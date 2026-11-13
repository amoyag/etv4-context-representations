"""Sequence-only completion of protocol baselines for exact compared inputs.
Adds region-wise identity, explicitly handling missing/nonoverlapping flanks.
No model, metric selection or inference.
"""
import csv,json
from pathlib import Path
from run_experiments import ROOT,read,write_json,digest
AA='ACDEFGHIKLMNPQRSTVWY'
GEN=['WT','Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']

def build():
    inputs=read('data/derived/model-inputs.json');prompts=read('data/derived/prompt-panel.json')
    pairs=set()
    def add(a,b,kind):pairs.add((a,b,kind))
    for g in GEN:
        add(f'ETV4_1_484_{g}',f'ETV4_337_436_{g}','FL_vs_CID')
        for bounds in ['1_484','337_436']:add(f'ETV4_{bounds}_{g}',f'ETV4_{bounds}_WT','genotype_vs_WT')
        for s in range(1,21):
            add(f'FLshuffle_{s}_{g}',f'ETV4_337_436_{g}','shuffled_FL_vs_CID')
            add(f'FLshuffle_{s}_{g}',f'ETV4_1_484_{g}','shuffled_vs_native_FL')
            add(f'FLshuffle_{s}_{g}',f'FLshuffle_{s}_WT','shuffled_genotype_vs_WT')
    for a,b in [('165_436','337_436'),('165_430','337_430'),('337_436','337_430'),('165_436','165_430'),('337_484','337_436'),('165_484','337_436'),('1_484','337_436')]:
        add('ETV4_'+a+'_WT','ETV4_'+b+'_WT','WT_boundary_or_factorial')
    for s in range(1,21):
        for a,b in [(f'Nshuffle_{s}_430','ETV4_337_430_WT'),(f'Cshuffle_{s}','ETV4_337_430_WT'),(f'NCshuffle_{s}',f'Nshuffle_{s}_430'),(f'NCshuffle_{s}',f'Cshuffle_{s}'),(f'Nshuffle_{s}_436','ETV4_337_436_WT'),(f'Nshuffle_{s}_436','ETV4_165_436_WT')]:add(a,b,'WT_shuffle_control')
    rows=[]
    regions={'whole':(1,484),'ETS_339_420':(339,420),'N_flank':(1,338),'C_flank':(421,484)}
    def subseq(x,lo,hi):
        start=x['native_start'];return {i:x['sequence'][i-start] for i in range(max(start,lo),min(start+x['length']-1,hi)+1)}
    for a,b,kind in sorted(pairs):
      for region,(lo,hi) in regions.items():
        x,y=subseq(inputs[a],lo,hi),subseq(inputs[b],lo,hi);shared=sorted(x.keys()&y.keys())
        row={'input_a':a,'input_b':b,'comparison':kind,'region':region,'length_a':len(x),'length_b':len(y),'aligned_positions':len(shared),
          'native_coordinate_identity':sum(x[i]==y[i] for i in shared)/len(shared) if shared else None,
          'identity_status':'defined_on_shared_positions' if shared else 'undefined_no_shared_positions'}
        for name,seq in [('a',''.join(x.values())),('b',''.join(y.values()))]:
            row['charge_'+name]=(seq.count('K')+seq.count('R')-seq.count('D')-seq.count('E'))/len(seq) if seq else None
            for aa in AA:row[f'frequency_{name}_{aa}']=seq.count(aa)/len(seq) if seq else None
        rows.append(row)
    def save(path,rows):
        assert not path.exists()
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        back=list(csv.DictReader(path.open(),delimiter='\t'));assert len(back)==len(rows)
        for original,parsed in zip(rows,back):
            for k,v in original.items():
                if v is None:assert parsed[k]==''
                elif isinstance(v,(int,float)):assert float(parsed[k])==v
                else:assert parsed[k]==v
    save(ROOT/'results/paired-sequence-baselines.tsv',rows)
    prows=[]
    for p,x in prompts.items():
        for i,seq in enumerate(x['sequences']):
            prows.append({'prompt_id':p,'member':i+1,'length':len(seq),'charge_proxy':(seq.count('K')+seq.count('R')-seq.count('D')-seq.count('E'))/len(seq),**{'frequency_'+a:seq.count(a)/len(seq) for a in AA}})
    save(ROOT/'results/prompt-sequence-baselines.tsv',prows)
    out={'pair_count':len(pairs),'region_rows':len(rows),'prompt_members':len(prows),'inputs_sha256':digest(ROOT/'data/derived/model-inputs.json'),'prompts_sha256':digest(ROOT/'data/derived/prompt-panel.json'),'round_trip_passed':True,
      'scope':'Protocol sequence baselines only; aligned-region identity is undefined without shared native positions. No new response metric or fitted model.'}
    write_json(ROOT/'results/paired-baseline-verification.json',out);print(json.dumps(out,indent=2))
if __name__=='__main__':build()
