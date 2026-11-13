"""Build review figures from frozen summaries only. No model loading or inference."""
from pathlib import Path
import csv,hashlib,json,sys,os,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle
HERE=Path(__file__).resolve().parent
PILOT=HERE/'inputs'
for folder in ['data','figures']:
 (HERE/folder).mkdir(exist_ok=True)
G=['Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']
LABEL=[g.replace('_','/') for g in G]
SPECS=[('650M','esm2_650m','22','exports/results-v1/results/E2-analysis-v1.json','dc5d35b481c9de6ac8d6f88a3c2d76928c0043af0a327740ac5881b0914dc502'),('3B','esm2_3b','24','followups/esm2-3b-20261111/job-3647838/exports/results-v1/e2/results/E2-3B-analysis.json','d15441709fb0ebd85b302ffc182f5ae1620f2387b23fc87bc36f5f3242b69f2c')]
provenance=[]; models={}; source={}; primary=[]; shuffles=[]; wt=[]; margins=[]; allmetrics=[]
for name,key,layer,rel,expected in SPECS:
 p=PILOT/rel;h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==expected
 provenance.append({'path':rel,'sha256':h});d=json.loads(p.read_text());source[name]=d;models[name]=d['models'][key]
 for l,v in models[name].items():
  for g in G:
   row={'model':name,'block':int(l),'genotype':g,**v['genotype_native'][g]};allmetrics.append(row)
   assert abs(row['M']/row['D']-row['M_over_D'])<1e-12
   assert abs(max(0,row['M']-8*v['eta'])/(row['D']+4*v['eta'])-row['L'])<1e-12
   if l==layer:
    primary.append(row)
    for i,q in enumerate(v['genotype_shuffled'][g],1):shuffles.append({'model':name,'block':int(l),'genotype':g,'draw':i,'M':q['M']})
   count=sum(x['M']>=row['M'] for x in v['genotype_shuffled'][g]);assert count==v['interpretation']['shuffle_counts_at_or_above_native'][g]
  native=v['genotype_native'][G[0]]['M'];maximum=max(x['M'] for x in v['genotype_shuffled'][G[0]])
  margin=native-maximum;threshold=8*v['eta'];passed=margin>threshold
  assert passed==v['interpretation']['conditions']['native_order_specificity']
  margins.append({'model':name,'block':int(l),'primary':l==layer,'native_M':native,'max_shuffle_M':maximum,'margin':margin,'required_8eta':threshold,'passes':passed})
  if l==layer:
   wt.append({'model':name,'block':int(l),'kind':'native','draw':0,'rms_J':v['factorial_native']['rms']})
   wt.extend({'model':name,'block':int(l),'kind':'shuffle','draw':i,'rms_J':q['rms']} for i,q in enumerate(v['factorial_shuffled'],1))
prof=source['650M']['models']['profam1'];pr=[];ps=[];pc=[];po=[];pd=[]
for l,v in prof.items():
 for name,q in v['prompt_response'].items():
  kind='empty' if name=='empty' else name.rsplit('_',2)[0]
  pr.append({'block':int(l),'prompt':name,'kind':kind,'rms_q':q['rms']})
 for name,rows in v['N_shuffles'].items():
  ps.extend({'block':int(l),'prompt':name,'draw':i,'rms_q':q['rms']} for i,q in enumerate(rows,1))
 for q in v['prompt_contrasts']:pc.append({'block':int(l),**q})
 for name,value in v['order_contrasts'].items():po.append({'block':int(l),'prompt_pair':name,'rms_order_vector_difference':value})
 for s in range(1,6):
  for o in ['ab','ba']:
   h=v['prompt_response'][f'hom_{s}_{o}']['rms'];n=v['prompt_response'][f'nonhom_{s}_{o}']['rms']
   pd.append({'block':int(l),'set':s,'order':o,'hom_norm':h,'nonhom_norm':n,'difference_of_norms':h-n})
 assert v['factorial_native']['rms']==0
assert len(primary)==10 and len(shuffles)==200 and len(wt)==42 and len(pr)==93 and len(ps)==420
files={'esm-primary.tsv':primary,'esm-all-fixed-blocks.tsv':allmetrics,'esm-primary-shuffles.tsv':shuffles,'esm-wt-factorial.tsv':wt,'esm-native-order-margins.tsv':margins,'profam-prompt-responses.tsv':pr,'profam-NID-shuffles.tsv':ps,'profam-vector-contrasts.tsv':pc,'profam-order-contrasts.tsv':po,'profam-matched-norm-differences.tsv':pd}
for name,rows in files.items():
 with (HERE/'data'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
# Supporting identity and biological tables: no long sequences in presentation output.
for rel in ['exports/results-v1/data/derived/construct-registry.tsv','exports/results-v1/data/derived/measurement-registry.tsv','exports/results-v1/data/derived/prompt-panel.json','design/protocol-v1.md','sources/source-audit.md']:
 p=PILOT/rel;provenance.append({'path':rel,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
with (PILOT/'exports/results-v1/data/derived/measurement-registry.tsv').open() as f:
 measures=list(csv.DictReader(f,delimiter='\t'))
kd=[r for r in measures if r['endpoint']=='equilibrium_KD']
assert len(kd)==7
with (HERE/'data/published-WT-affinities.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['construct','KD_nM','SEM_nM','reported_n'],delimiter='\t');w.writeheader()
 for r in kd:w.writerow({'construct':r['construct_id'],'KD_nM':float(r['value'])*.01,'SEM_nM':float(r['uncertainty'])*.01,'reported_n':r['n']})
fontdir=Path(os.environ.get('LIBERTINUS_FONT_DIR','/usr/local/texlive/2025/texmf-dist/fonts/opentype/public/libertinus-fonts'))
if not fontdir.exists():
 found=subprocess.run(['kpsewhich','LibertinusSans-Regular.otf'],capture_output=True,text=True,check=True).stdout.strip()
 if not found: raise RuntimeError('Set LIBERTINUS_FONT_DIR to the directory containing LibertinusSans OTF fonts')
 fontdir=Path(found).parent
for p in fontdir.glob('LibertinusSans-*.otf'):fm.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Libertinus Sans','font.size':11,'axes.titlesize':12,'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10,'pdf.fonttype':3,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.titlelocation':'left','figure.facecolor':'white'})
BLUE='#243B53';RED='#A64B32';TEAL='#287D7A';GREY='#9BA7B2';COL={'650M':BLUE,'3B':RED};PCOL={'hom':TEAL,'nonhom':RED,'homshuffle':BLUE}
def clean(ax):
 ax.grid(axis='x',color='#E6E9EC',lw=.6);ax.set_axisbelow(True)
def save(fig,stem):
 for ax in fig.axes:
  if ax.axison:
   if stem=='fig4-profam-conditioning' and ax in [fig.axes[0],fig.axes[1],fig.axes[3]]:ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
   elif stem=='fig3-order-controls' and ax==fig.axes[2]:ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{1000*x:.3f}'))
   else:ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
 for ext in ['pdf','svg','png']:fig.savefig(HERE/'figures'/f'{stem}.{ext}',dpi=200,bbox_inches='tight')
 plt.close(fig)
def box(ax,x,y,w,h,text,color='#EDF0F2',size=11):
 ax.add_patch(Rectangle((x,y),w,h,facecolor=color,edgecolor='white',lw=1.5));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=size)
# Figure 1: regions and the two distinct comparisons.
fig=plt.figure(figsize=(9,8));gs=fig.add_gridspec(3,1,height_ratios=[1.05,1.35,1.65],hspace=.38)
a=fig.add_subplot(gs[0]);a.axis('off');a.set_xlim(0,1);a.set_ylim(0,1);a.set_title('A  Named regions of full-length ETV4',fontweight='bold')
regions=[('N-terminal\nextension','1–164','#E3E7EB',.20),('NID','165–336','#D3E6F3',.20),('Junction','337–338','#E3E7EB',.12),('ETS domain','339–420','#EBC8B9',.20),('CID segment','421–436','#CBE5DE',.15),('C-terminal\ntail','437–484','#E3E7EB',.13)]
x=0
for name,span,color,w in regions:box(a,x,.35,w,.50,name+'\n'+span,color,10.5);x+=w
a.text(0,.10,'Schematic widths are not sequence lengths. CID includes H4; CID is not wholly disordered.',fontsize=10)
b=fig.add_subplot(gs[1]);b.axis('off');b.set_title('B  Published mutant contrasts guide the ESM comparison',fontweight='bold')
table=b.table(cellText=[['Y401N/Y403G','ETS','Weaker DNA binding; relieves NID-associated inhibition'],['E404K','ETS','Stronger DNA binding; inhibition reported preserved'],['E423K; E425K','CID','Binding and inhibition reported preserved'],['F428A/F432A','CID','Binding and inhibition reported preserved']],colLabels=['Substitution(s)','Region','Published observation'],colWidths=[.22,.10,.68],bbox=[0,.03,1,.89],cellLoc='left',colLoc='left')
table.auto_set_font_size(False);table.set_fontsize(10)
for (r,c),cell in table.get_celld().items():
 cell.set_edgecolor('white');cell.set_facecolor('#E7EDF2' if r==0 else ('#F4F6F8' if r%2 else 'white'))
 if r==0:cell.set_text_props(weight='bold')
c=fig.add_subplot(gs[2]);c.axis('off');c.set_xlim(0,1);c.set_ylim(0,1);c.set_title('C  Two ways to vary the information available to the model',fontweight='bold')
box(c,0,.61,.49,.31,'ESM-2: mutation × sequence context\nETS + CID (337–436) versus full length (1–484)\nWT and the same five mutant genotypes','#F0E6E0',10.5)
box(c,.51,.61,.49,.31,'ProFam: NID response × protein prompt\nNID + ETS + CID (165–436) versus ETS + CID\nWT targets; the same prompt precedes both','#E2EFEB',10.5)
c.text(.01,.45,'ESM readout: 88 unchanged positions',fontweight='bold',fontsize=10.5)
c.text(.01,.31,'79 ETS + 7 CID + 2 junction positions\nMutation sites excluded with one common mask.',fontsize=10)
c.text(.52,.45,'ProFam WT readout: all 94 shared positions',fontweight='bold',fontsize=10.5)
c.text(.52,.31,'82 ETS + 10 CID + 2 junction positions\nPrompt: related / unrelated / shuffled / empty.',fontsize=10)
c.text(.01,.07,'Inputs include complete named constructs. Both readouts lie within 337–430.\nProFam is left-to-right: the readout cannot see the later CID residues 431–436.',fontsize=10)
save(fig,'fig1-design')
# Figure 2: all genotypes, both models, no inferential cross-model ranking.
fig,axes=plt.subplots(1,3,figsize=(9,4.6),sharey=True,gridspec_kw={'width_ratios':[1,1,1.25]})
for ax,key,title in zip(axes,['M','D','M_over_D'],['A  Mutation-dependent\ncontext effect','B  Mutation effect\nwithin ETS + CID','C  Original ratio\ndiagnostic']):
 for model,offset,marker in [('650M',-.12,'o'),('3B',.12,'s')]:
  rows=[r for r in primary if r['model']==model];ys=np.arange(5)+offset
  if key=='M_over_D':
   for y,r in zip(ys,rows):ax.plot([r['L'],r['U']],[y,y],color=COL[model],lw=1.7)
  ax.scatter([r[key] for r in rows],ys,c=COL[model],marker=marker,s=35,label=f'{model}: block '+('22' if model=='650M' else '24'),zorder=3)
 ax.set_title(title,fontsize=10.5,pad=13);ax.set_xlabel({'M':'M (RMS displacement)','D':'D (RMS displacement)','M_over_D':'M/D; numerical bounds'}[key]);ax.set_xlim(left=0);clean(ax)
axes[0].set_yticks(range(5),LABEL);axes[0].invert_yaxis();axes[1].set_xlim(0,.080);axes[0].set_xlim(0,.080);axes[0].set_xticks([0,.020,.040,.060,.080]);axes[1].set_xticks([0,.020,.040,.060,.080]);axes[2].set_xlim(0,.94)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.56,-.015));fig.subplots_adjust(left=.18,right=.99,bottom=.18,wspace=.28)
save(fig,'fig2-mutation-responses')
# Figure 3: native-order distributions and the two different WT/mutant controls.
fig,axes=plt.subplots(2,2,figsize=(9,7.7),gridspec_kw={'height_ratios':[1.2,1]});fig.subplots_adjust(left=.17,right=.98,top=.94,bottom=.10,wspace=.31,hspace=.55)
for ax,model,title in zip(axes[0],['650M','3B'],['A  650M, block 22','B  3B, block 24']):
 v=models[model]['22' if model=='650M' else '24']
 for y,g in enumerate(G):
  vals=[r['M'] for r in v['genotype_shuffled'][g]];ax.scatter(vals,y+np.linspace(-.16,.16,20),color=GREY,s=13,alpha=.85)
  native=v['genotype_native'][g]['M'];ax.scatter(native,y,color=COL[model],marker='D',s=45,zorder=4)
  count=sum(x>=native for x in vals);ax.text(.0382,y,f'{count}/20',va='center',fontsize=9)
 ax.set_xlim(0,.044);ax.set_xticks([0,.01,.02,.03]);ax.set_ylim(4.5,-.65);ax.set_yticks(range(5),LABEL if model=='650M' else ['']*5);ax.set_xlabel('M: native diamond; shuffled dots');ax.set_title(title,fontweight='bold');clean(ax)
ax=axes[1,0]
for y,r in enumerate(margins):
 ax.plot([r['margin'],r['required_8eta']],[y,y],color='#C6CDD3',lw=1.5)
 ax.scatter(r['margin'],y,color=COL[r['model']],marker='D',s=34)
 ax.scatter(r['required_8eta'],y,edgecolor=COL[r['model']],facecolor='white',marker='s',s=34,zorder=3)
ax.set_yticks(range(6),[r['model']+' / '+str(r['block'])+('*' if r['primary'] else '') for r in margins]);ax.invert_yaxis();ax.axvline(0,color='#777',lw=.6);ax.set_title('C  Fixed Y401N/Y403G order test',fontweight='bold');ax.set_xlabel('Margin (◆) and required 8η (□), in units of 10⁻³',fontsize=9);clean(ax)
ax=axes[1,1]
for y,model in enumerate(['650M','3B']):
 rows=[r for r in wt if r['model']==model];vals=[r['rms_J'] for r in rows if r['kind']=='shuffle'];native=next(r['rms_J'] for r in rows if r['kind']=='native')
 ax.scatter(vals,y+np.linspace(-.13,.13,20),color=GREY,s=16);ax.scatter(native,y,marker='D',color=COL[model],s=45,zorder=3);ax.text(.022,y,f'{sum(x>=native for x in vals)}/20',va='center',fontsize=9)
ax.set_xlim(0,.026);ax.set_ylim(1.65,-.65);ax.set_yticks([0,1],['650M / 22','3B / 24']);ax.set_title('D  WT factorial control',fontweight='bold');ax.set_xlabel('RMS(J): native diamond; shuffled dots');clean(ax)
save(fig,'fig3-order-controls')
# Figure 4: main-text ProFam candidate, with full primary and depth results.
fig,axes=plt.subplots(2,2,figsize=(9,7.8));fig.subplots_adjust(left=.10,right=.98,top=.93,bottom=.11,wspace=.32,hspace=.57)
v=prof['11'];ax=axes[0,0];settings=[(s,o) for s in range(1,6) for o in ['ab','ba']]
for cat,marker,label in [('hom','o','Homologous'),('nonhom','^','Nonhomologous'),('homshuffle','s','Shuffled homologues')]:
 vals=[v['prompt_response'][f'{cat}_{s}_{o}']['rms'] for s,o in settings];ax.scatter(np.arange(10),vals,color=PCOL[cat],marker=marker,s=26,label=label)
ax.axhline(v['prompt_response']['empty']['rms'],color='#777',ls='--',lw=1,label='Empty prompt');ax.set_xticks(range(10),[str(s)+o for s,o in settings],rotation=45);ax.set_ylim(0,.55);ax.set_title('A  Primary block 11: all prompt settings',fontweight='bold',fontsize=11);ax.set_ylabel('NID response RMS(q)');ax.set_xlabel('Matched set and sequence order');ax.grid(axis='y',color='#E6E9EC',lw=.5)
ax=axes[0,1]
for j,l in enumerate([8,11,16]):
 vals=[r['difference_of_norms'] for r in pd if r['block']==l];ax.scatter(j+np.linspace(-.16,.16,10),vals,c=TEAL,s=26)
ax.axhline(0,color='#555',lw=.8);ax.set_xticks(range(3),['8','11 (primary)','16']);ax.set_ylabel('RMS(q homologous) − RMS(q nonhomologous)',fontsize=9);ax.set_xlabel('Fixed block');ax.set_title('B  Matched ordering changes with depth',fontweight='bold',fontsize=11);ax.set_ylim(-.13,.04);ax.grid(axis='y',color='#E6E9EC',lw=.5)
ax=axes[1,0];names=['empty','hom_1_ab','hom_1_ba','nonhom_1_ab','nonhom_1_ba','homshuffle_1_ab','homshuffle_1_ba'];short=['Empty','Hom 1ab','Hom 1ba','Nonhom 1ab','Nonhom 1ba','Shuffled 1ab','Shuffled 1ba']
for y,name in enumerate(names):
 vals=[r['rms'] for r in v['N_shuffles'][name]];native=v['prompt_response'][name]['rms'];ax.scatter(vals,y+np.linspace(-.17,.17,20),c=GREY,s=10);ax.scatter(native,y,marker='D',color=TEAL,s=30,zorder=3);ax.text(.545,y,f'{sum(x>=native for x in vals)}/20',va='center',fontsize=8)
ax.set_yticks(range(7),short,fontsize=9);ax.set_ylim(6.6,-.6);ax.set_xlim(0,.615);ax.set_xlabel('RMS(q): native diamond; NID-shuffled dots');ax.set_title('C  Set-1 NID-order controls, block 11',fontweight='bold',fontsize=11);clean(ax)
ax=axes[1,1]
for i,(cat,label) in enumerate([('hom','Homologous'),('nonhom','Nonhomologous'),('homshuffle','Shuffled\nhomologues')]):
 vals=[v['order_contrasts'][f'{cat}_{s}'] for s in range(1,6)];ax.scatter(i+np.linspace(-.12,.12,5),vals,c=PCOL[cat],s=28)
ax.set_xticks(range(3),['Homologous','Nonhomologous','Shuffled\nhomologues'],fontsize=9);ax.set_ylim(0,.18);ax.set_ylabel('RMS(q order ab − q order ba)',fontsize=10);ax.set_title('D  Prompt-order effects, block 11',fontweight='bold',fontsize=11);ax.grid(axis='y',color='#E6E9EC',lw=.5)
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(.5,-.005),fontsize=9)
save(fig,'fig4-profam-conditioning')
# Round-trip check the values actually exported for the plots.
for name,rows in files.items():
 with (HERE/'data'/name).open() as f:back=list(csv.DictReader(f,delimiter='\t'))
 assert len(back)==len(rows)
 for a,b in zip(rows,back):
  for k,value in a.items():
   if isinstance(value,float): assert float(b[k])==value
proof={'sources':provenance,'scope':'Frozen-summary extraction and presentation only; no inference or threshold changes','source_rows':{k:len(v) for k,v in files.items()},'checks':['Source SHA256 identities','Ratios and lower numerical bound arithmetic','Native-order classification and all reported genotype counts','Table float round trips','Seven published WT KD rows converted from 1e-11 M to nM'],'figure_note':'New presentation under Checkpoint B; original exports and figures untouched','independent_review':False}
(HERE/'data-provenance.json').write_text(json.dumps(proof,indent=2)+'\n')
print('Built four figures in PDF/SVG/PNG and eleven source tables; presentation assertions passed.')
