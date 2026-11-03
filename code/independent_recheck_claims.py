"""Read-only audit of deposited sequence labels and roles; never trains models."""
from pathlib import Path
import csv, json, re
B = Path(__file__).resolve().parents[1]
O = Path(__file__).resolve().parent

def fasta(path):
    records = []; header = None; seq = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if header is not None: records.append((header, ''.join(seq).upper()))
            header = line; seq = []
        else: seq.append(re.sub(r'\s+', '', line))
    if header is not None: records.append((header, ''.join(seq).upper()))
    return records

D = {}
files = {
'H1': ('H1_iDPPIV-SCM/source/DPP-IV', ['train-positive.txt','train-negative.txt','test-positive.txt','test-negative.txt']),
'H2': ('H2_iDPPIV-SI/source/iDPPIV-SI/data', ['train_positive.txt','train_negative.txt','test_positive.txt','test_negative.txt']),
'H3': ('H3_iDPPIV-PC/source/iDPPIV/data', ['train_positive.txt','train_negative.txt','test_positive.txt','test_negative.txt']),
'H5': ('H5_IPPF-FE/source/IPPF-FE-main/Dipeptidyl_peptidase‑IV_inhibitory_peptides', ['train-positive.fasta','train-negative.fasta','test-positive.fasta','test-negative.fasta']),
'A3': ('A3_StackDPPIV/source', ['TR_pos_dppiv.txt','TR_neg_dppiv.txt','TS_pos_dppiv.txt','TS_neg_dppiv.txt']),
'A5': ('A5_DL设计规则/source/Fine_tune_data/DPP-IV_Dataset', ['train-positive.txt','train-negative.txt','test-positive.txt','test-negative.txt']),
'A13': ('A13_FedAntiDiaPep/source/FedAntiDiaPep-main', ['Train-Positive.fasta','Train-Negative.fasta','Test-Positive.fasta','Test-Negative.fasta'])}
for model,(folder,names) in files.items():
    for (split,label),name in zip([('train',1),('train',0),('test',1),('test',0)],names):
        D[model,split,label] = {s for h,s in fasta(B/folder/name)}
for split in ['train','test']:
    p = B/f'A4_StructuralDPPIV/source/Structural-DPP-IV-main/data/DPP-IV/{split}/{split}.tsv'
    rows = list(csv.DictReader(p.open(),delimiter='\t'))
    for label in [0,1]: D['A4',split,label] = {r['sequence'].upper() for r in rows if int(r['label'])==label}
for split,name in [('train','Layer1_training.txt'),('test','L1_Ind.txt')]:
    rows = fasta(B/'A2_ADP-Fuse/source'/name)
    for label,word in [(0,'negative'),(1,'positive')]: D['A2',split,label] = {s for h,s in rows if word in h.lower()}
rows = list(csv.DictReader((B/'A9_BertADP/source/BertADP-main/data/allData.csv').open()))
for split in ['train','test']:
    for label in [0,1]: D['A9',split,label] = {r['Sequence'] for r in rows if r['Source']==split and int(r['Label'])==label}
for split,name in [('train','Train.csv'),('test','Test.csv'),('external','ExternalValidation.csv')]:
    rows = list(csv.DictReader((B/'B1_Bilingua-DPP-IV/source/data'/name).open()))
    for label in [0,1]: D['B1',split,label] = {r['seq'] for r in rows if int(r['label'])==label}

def union(model,label): return set().union(*(s for (m,sp,l),s in D.items() if m==model and l==label))
conflicts = union('H1',0)&union('H1',1)
models = ['H1','H2','H3','H5','A3','A4','A5']
results = {
'negative_sets_and_deposited_splits_identical': {m:all(D[m,sp,0]==D['H1',sp,0] for sp in ['train','test']) for m in models},
'conflicts_by_model': {m:sorted(union(m,0)&union(m,1)) for m in models},
'cross_split_any_label': {m:sorted((D[m,'train',0]|D[m,'train',1])&(D[m,'test',0]|D[m,'test',1])) for m in models},
'A9_role_matrix_from_A2_negatives': {sp2:{sp9:len(D['A2',sp2,0]&D['A9',sp9,0]) for sp9 in ['train','test']} for sp2 in ['train','test']},
'A9_pool_external_intersection': {str(l):len(D['A9','train',l]&D['A9','test',l]) for l in [0,1]},
'A9_A13_published_sets_match': all(D['A9',sp,l]==D['A13',sp,l] for sp in ['train','test'] for l in [0,1]),
'B1_original_negative_sequences_removed': sorted(union('H1',0)-union('B1',0)),
'B1_removed_exactly_all_seven_conflicts': union('H1',0)-union('B1',0)==conflicts,
'B1_seven_conflicts_absent_both_labels': not (conflicts & (union('B1',0)|union('B1',1))),
'B1_external_overlap_H1_any_label': {str(l):len(D['B1','external',l]&(union('H1',0)|union('H1',1))) for l in [0,1]},
'limits': ['Exact full sequence identity only, not homology or assay truth.', 'A9 train field is the development pool; final training/validation assignment is not available.', 'Deposited files are not historical execution logs.', 'No performance effect estimated.']}
(O/'strict_direction_results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
with (O/'f665_conflicting_labels.csv').open('w',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['sequence','model','roles'])
    for seq in sorted(conflicts):
        for m in models:w.writerow([seq,m,';'.join(f'{sp}:{l}' for sp in ['train','test'] for l in [0,1] if seq in D[m,sp,l])])
print(json.dumps(results,ensure_ascii=False,indent=2))
