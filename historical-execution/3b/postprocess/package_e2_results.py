"""Build result tables and two ESM-2 3B follow-up figures from verified analysis.
Presentation choices do not change models, blocks, controls, metrics or thresholds.
"""
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_experiments import ROOT, digest, write_json
GEN=['Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']
LABEL=['Y401N/Y403G','E404K','E423K','E425K','F428A/F432A']

def table(path,rows,fields):
    assert not path.exists(),path
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    # Verify serialized numbers and nulls survive export; independent raw arithmetic is verify_e2.py.
    reread=list(csv.DictReader(path.open(),delimiter='\t'))
    assert len(reread)==len(rows)
    for a,b in zip(rows,reread):
        for k,v in a.items():
            if v is None:assert b[k]==''
            elif isinstance(v,(float,int)):assert float(b[k])==v
            else:assert str(v)==b[k]

def build(analysis):
    r=json.loads(analysis.read_text());out=ROOT/'results';figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
    primary=[];controls=[];profiles=[];baselines=[]
    for model,layers in r['models'].items():
        for layer,b in layers.items():
            for g,m in b['genotype_native'].items():
                primary.append({'result_id':f'{model}/L{layer}/genotype_native/{g}','model':model,'layer':int(layer),'genotype':g,
                  **m,'eta':b['eta'],'denominator_floor':4*b['eta'],'separation_floor':8*b['eta'],
                  'scope':'primary' if model=='esm2_3b' and layer=='24' else 'sensitivity' if model=='esm2_3b' else 'secondary_partial_context',
                  'classification':b.get('interpretation',{}).get('classification','not_classified')})
            def walk(x,path):
                if isinstance(x,dict):
                    for key,v in x.items():
                        if key in ['profile','context_change_norm','CID_change_norm']:
                            positions=x.get('native_positions',list(range(337,431)))
                            for pos,val in zip(positions,v):profiles.append({'result_id':path+'/'+key,'model':model,'layer':int(layer),'native_position':pos,'value':val})
                        elif key not in ['seed','set','native_positions']:walk(v,path+'/'+key)
                elif isinstance(x,list):
                    for i,v in enumerate(x):walk(v,path+f'/{i}')
                elif isinstance(x,(float,int)) and not isinstance(x,bool):
                    controls.append({'result_id':path,'model':model,'layer':int(layer),'value':x,'value_status':'finite'})
                elif x is None:controls.append({'result_id':path,'model':model,'layer':int(layer),'value':None,'value_status':'undefined_denominator'})
            for key,x in b.items():
                if key not in ['genotype_native','eta']:walk(x,f'{model}/L{layer}/{key}')
    for ident,b in r['baselines'].items():
        for region,x in b.items():
            if region=='native_position_identity':
                baselines.append({'input_id':ident,'region':'whole','metric':region,'value':x});continue
            if x is None:continue
            for key,val in x.items():
                if key=='frequencies':
                    for aa,v in val.items():baselines.append({'input_id':ident,'region':region,'metric':'frequency_'+aa,'value':v})
                else:baselines.append({'input_id':ident,'region':region,'metric':key,'value':val})
    table(out/'primary-metrics.tsv',primary,list(primary[0]))
    table(out/'control-metrics.tsv',controls,list(controls[0]))
    table(out/'profiles.tsv',profiles,list(profiles[0]))
    table(out/'sequence-baselines.tsv',baselines,list(baselines[0]))
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':150})
    saved=[];datasets={}
    def save(fig,n,data):
        for suffix in ['pdf','png']:assert not (figdir/f'fig{n}.{suffix}').exists()
        fig.savefig(figdir/f'fig{n}.pdf',bbox_inches='tight',metadata={'Title':f'ETS E2 Figure {n}','Creator':'package_e2_results.py'})
        fig.savefig(figdir/f'fig{n}.png',bbox_inches='tight');plt.close(fig)
        write_json(figdir/f'fig{n}-data.json',{'analysis_path':str(analysis.relative_to(ROOT)),'analysis_sha256':digest(analysis),'data':data})
        saved.append(f'figures/fig{n}.pdf')
    b=r['models']['esm2_3b']['24'];fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for i,g in enumerate(GEN):
        vals=[v['M'] for v in b['genotype_shuffled'][g]]
        axs[0].scatter(np.full(20,i)+np.linspace(-.15,.15,20),vals,s=14,c='#9aa2ac',alpha=.8,label='20 fixed prefix shuffles' if i==0 else None)
        axs[0].scatter(i,b['genotype_native'][g]['M'],marker='D',s=45,c='#006f91',label='Native sequence' if i==0 else None,zorder=3)
    axs[0].set(ylabel='M: context-change RMS',title='A  Native response and prefix-order controls');axs[0].legend(fontsize=8)
    for i,g in enumerate(GEN):
        m=b['genotype_native'][g]
        if m['M_over_D'] is not None:
            axs[1].errorbar(i,m['M_over_D'],yerr=[[m['M_over_D']-m['L']],[m['U']-m['M_over_D']]],fmt='D',color='#006f91',capsize=4)
        else:axs[1].text(i,.02,'undefined',rotation=90,ha='center',transform=axs[1].get_xaxis_transform())
    axs[1].set(ylabel='M / D (CID-construct diagnostic)',title='B  Numerical sensitivity bounds')
    for ax in axs:ax.set_xticks(range(5),LABEL,rotation=30,ha='right');ax.grid(axis='y',alpha=.2)
    fig.suptitle('ESM-2 3B block 24 | '+b['interpretation']['classification'].replace('_',' '),fontsize=12)
    fig.supxlabel('One protein and five genotype comparisons. Bounds reflect eta, not biological confidence intervals.',fontsize=8)
    save(fig,1,{k:b[k] for k in ['eta','genotype_native','genotype_shuffled','interpretation']})
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for i,l in enumerate(['19','24','36']):
        bb=r['models']['esm2_3b'][l];vals=[x['rms'] for x in bb['factorial_shuffled']]
        axs[0,0].scatter(i+np.linspace(-.16,.16,20),vals,c='#9aa2ac',s=18)
        axs[0,0].scatter(i,bb['factorial_native']['rms'],marker='D',s=55,c='#006f91')
    axs[0,0].set(xticks=range(3),xticklabels=['19','24 (primary)','36'],ylabel='RMS(J)',title='A  WT factorial: native and 20 paired shuffles')
    pos=np.arange(337,431);control=np.array([x['profile'] for x in b['factorial_shuffled']])
    axs[0,1].fill_between(pos,control.min(0),control.max(0),color='#dbe0e6',label='Shuffle min-max')
    axs[0,1].plot(pos,b['factorial_native']['profile'],color='#006f91',label='Native')
    axs[0,1].set(xlabel='Native ETV4 residue',ylabel='Per-residue ||J||',title='B  Block 24 factorial profile');axs[0,1].legend(fontsize=8)
    for i,l in enumerate(['19','24','36']):
        bb=r['models']['esm2_3b'][l];vals=[x['rms'] for x in bb['N_only_shuffled_response']]
        axs[1,0].scatter(i+np.linspace(-.16,.16,20),vals,c='#9aa2ac',s=18)
        axs[1,0].scatter(i,bb['native_N_response_with_C']['rms'],marker='D',s=55,c='#006f91')
    axs[1,0].set(xticks=range(3),xticklabels=['19','24 (primary)','36'],ylabel='N-response RMS',title='C  N-order controls with the native C flank')
    axs[1,1].bar(range(4),[6.1,12,12,60],color=['#7d8c9e','#a1b5bd','#a1b5bd','#006f91'])
    axs[1,1].set(xticks=range(4),xticklabels=['0\n337-430','N\n165-430','C\n337-436','NC\n165-436'],ylabel='Reported mean KD (10^-11 M)',title=f"D  Source WT means | log ratio = {r['biological_factorial_log_ratio']:.4f}")
    for ax in axs.flat:ax.grid(axis='y',alpha=.2)
    fig.suptitle('WT flank interaction and order sensitivity',fontsize=13)
    fig.supxlabel('Grey: fixed computational controls; blue diamond: native. Source means: Currie et al., Table S3.\nNo biological uncertainty or calibration between KD and representation norms is inferred.',fontsize=8)
    save(fig,2,{'esm':{l:{k:r['models']['esm2_3b'][l][k] for k in ['factorial_native','factorial_shuffled','native_N_response_with_C','N_only_shuffled_response']} for l in ['19','24','36']},'KD_means':[6.1,12,12,60],'biological_log_ratio':r['biological_factorial_log_ratio']})
    record={'analysis_sha256':digest(analysis),'primary_rows':len(primary),'control_rows':len(controls),'profile_rows':len(profiles),'baseline_rows':len(baselines),'round_trip_tables_passed':True,'figures':saved}
    write_json(out/'presentation-verification.json',record);print(json.dumps(record,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--analysis',required=True);a=p.parse_args();build(Path(a.analysis).resolve())
