"""Audited full-sequence parser. Run from any working directory; never alters raw data."""
from pathlib import Path
from collections import Counter, defaultdict
import csv, hashlib, json, re
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '02_解析'
VERSION = 'evidence-20261102-v1'
AA = set('ACDEFGHIKLMNPQRSTVWY')
RECORDS = []
SOURCES = {}


def write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None: fields = list(rows[0])
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def track(path):
    rel = path.relative_to(ROOT).as_posix()
    if rel not in SOURCES:
        SOURCES[rel] = {'source_file': rel, 'bytes': path.stat().st_size,
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    return rel


def add(model, dataset, split, label, seq, path, record, line, depad=False):
    raw = str(seq); norm = re.sub(r'\s+', '', raw).upper()
    if depad: norm = norm.rstrip('X')
    if len(norm) < 2 or set(norm) - AA:
        raise ValueError(f'Invalid sequence: {path}:{line} {raw!r}')
    RECORDS.append(dict(model=model, dataset=dataset, split=split, label=int(label),
        sequence=norm, source_file=track(path), source_record=record,
        source_line=line, sequence_raw=raw,
        normalization='uppercase;whitespace_removed' + (';trailing_X_padding_removed' if depad else '')))


def fasta(path, model, dataset, split, fixed_label=None, pdf_text=False, depad=False):
    header = None; chunks = []; first_line = None
    def flush():
        if header is None: return
        if fixed_label is not None: lab = fixed_label
        elif pdf_text:
            match = re.fullmatch(r'>[^|]+\|([01])\|\w+', header)
            if not match: raise ValueError(f'Unrecognized A1 header {header}')
            lab = int(match[1])
        elif re.search(r'negative', header, re.I): lab = 0
        elif re.search(r'positive', header, re.I): lab = 1
        else: raise ValueError(f'Ambiguous FASTA label {header}')
        add(model,dataset,split,lab,''.join(chunks),path,header,first_line,depad)
    for n, line in enumerate(path.read_text().splitlines(),1):
        s=line.strip()
        if not s: continue
        if s.startswith('>'):
            flush(); header=s; chunks=[]; first_line=None
        elif pdf_text and not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',s):
            # PDF browser print headers, footers and page markers are not sequence lines.
            continue
        else:
            if header is None: raise ValueError(f'Sequence before header at {path}:{n}')
            if first_line is None: first_line=n
            chunks.append(s)
    flush()


def tabular(path, model, dataset, split, seqcol, labelcol, delimiter=',', rolecol=None, rolemap=None):
    with path.open(encoding='utf-8-sig',newline='') as f:
        for n,r in enumerate(csv.DictReader(f,delimiter=delimiter),2):
            role=rolemap[r[rolecol]] if rolecol else split
            add(model,dataset,role,r[labelcol],r[seqcol],path,n-1,n)


def parse():
    pairs={
     'H1':('H1_iDPPIV-SCM/source/DPP-IV','train-positive.txt','train-negative.txt','test-positive.txt','test-negative.txt'),
     'H2':('H2_iDPPIV-SI/source/iDPPIV-SI/data','train_positive.txt','train_negative.txt','test_positive.txt','test_negative.txt'),
     'H3':('H3_iDPPIV-PC/source/iDPPIV/data','train_positive.txt','train_negative.txt','test_positive.txt','test_negative.txt'),
     'H5':('H5_IPPF-FE/source/IPPF-FE-main/Dipeptidyl_peptidase‑IV_inhibitory_peptides','train-positive.fasta','train-negative.fasta','test-positive.fasta','test-negative.fasta'),
     'A3':('A3_StackDPPIV/source','TR_pos_dppiv.txt','TR_neg_dppiv.txt','TS_pos_dppiv.txt','TS_neg_dppiv.txt'),
     'A5':('A5_DL设计规则/source/Fine_tune_data/DPP-IV_Dataset','train-positive.txt','train-negative.txt','test-positive.txt','test-negative.txt'),
     'A13':('A13_FedAntiDiaPep/source/FedAntiDiaPep-main','Train-Positive.fasta','Train-Negative.fasta','Test-Positive.fasta','Test-Negative.fasta')}
    for m,(folder,*files) in pairs.items():
        for (split,label),fn in zip([('train',1),('train',0),('test',1),('test',0)],files):
            fasta(ROOT/folder/fn,m,'released_binary_dataset',split,label)
    fasta(ROOT/'A1_AntiDMPpred/source/A1_positive_negative.fasta.txt','A1','released_benchmark','development_pool',pdf_text=True)
    for m,folder,files in [('A2','A2_ADP-Fuse/source',['Layer1_training.txt','L1_Ind.txt']),
                           ('A14','A14_DQN/source/DQNCode-main/Dataset',['Layer 1 Training.fasta','Layer 1 Independent.fasta'])]:
        for split,fn in zip(['train','test'],files):fasta(ROOT/folder/fn,m,'released_L1_binary_dataset',split)
    for split in ['train','test']:
        tabular(ROOT/f'A4_StructuralDPPIV/source/Structural-DPP-IV-main/data/DPP-IV/{split}/{split}.tsv','A4','released_binary_dataset',split,'sequence','label','\t')
    for fn,split,label in [('Training dataset-ADP.fasta','train',1),('Training dataset-non_ADP.fasta','train',0),('independent_test_ADP.fasta','test',1),('independent_testing_non-ADP.fasta','test',0)]:
        fasta(ROOT/'A7_Yue/source/Yue-main'/fn,'A7','repository_version_paper_mapping_unresolved',split,label,depad=True)
    tabular(ROOT/'A9_BertADP/source/BertADP-main/data/allData.csv','A9','released_development_and_external',None,'Sequence','Label',rolecol='Source',rolemap={'train':'development_pool','test':'external_test'})
    tabular(ROOT/'H6_SeqDPPIV/source/SeqDPPIV-main/SeqDPPIV model/DPPIV sequence and descriptors.csv','H6','actual_model_input_664pos_629neg','model_input_pool','sequence','label')
    for fn,split in [('Train.csv','train'),('Test.csv','test'),('ExternalValidation.csv','external_test')]:
        tabular(ROOT/'B1_Bilingua-DPP-IV/source/data'/fn,'B1','released_recurated_dataset',split,'seq','label')
    p=ROOT/'H4_BBATProt/source/BBATProt-master/codes/DPPIV/train.tsv'
    for n,r in enumerate(csv.reader(p.open(),delimiter='\t'),1):add('H4','available_training_only','train',r[0],r[1],p,n,n)
    p=ROOT/'A12_DeNovo/source/Anti-Diabetic-Peptide-Prediction--master/negSamples.xlsx'
    wb=openpyxl.load_workbook(p,read_only=True,data_only=True);sheet=wb.active;rows=list(sheet.values);col=rows[0].index('Sequence')
    for n,r in enumerate(rows[1:],2):
        if r[col] is not None:add('A12','repository_null_paper_mapping_unresolved','available_only',0,r[col],p,f'{sheet.title}!A{n}',n)
    wb.close()


def generate():
    parse()
    by_model=defaultdict(list);sets=defaultdict(set)
    for r in RECORDS:
        by_model[r['model']].append(r);sets[r['model'],r['split'],r['label']].add(r['sequence'])
    def total(m,label=0):return set().union(*(s for (mm,sp,l),s in sets.items() if mm==m and l==label))
    models=sorted(by_model,key=lambda x:('AHB'.index(x[0]),int(x[1:])))
    ordered=sorted(RECORDS,key=lambda r:(models.index(r['model']),r['split'],-r['label'],r['source_file'],r['source_line']))
    write_csv(OUT/'negative_sets_long.csv',ordered)
    summ=[];counts=[]
    for m in models:
        rr=by_model[m];nn=[r for r in rr if r['label']==0];pp=[r for r in rr if r['label']==1]
        summ.append(dict(model=m,rows=len(rr),pos=len(pp),neg=len(nn),unique_neg=len(total(m)),unique_pos=len(total(m,1)),duplicate_neg=len(nn)-len(total(m)),duplicate_pos=len(pp)-len(total(m,1)),splits=';'.join(sorted({r['split'] for r in rr})),dataset_scope=rr[0]['dataset']))
        for sp in sorted({r['split'] for r in rr}):
            for lab in [0,1]:
                z=[r for r in rr if r['split']==sp and r['label']==lab]
                if z:counts.append(dict(model=m,split=sp,label=lab,records=len(z),unique_sequences=len(sets[m,sp,lab]),dataset_scope=rr[0]['dataset']))
    write_csv(OUT/'negative_sets_summary.csv',summ);write_csv(OUT/'partition_counts.csv',counts)
    matrix=[dict(model=m,**{n:len(total(m)&total(n)) for n in models}) for m in models]
    write_csv(OUT/'negative_overlap_matrix.csv',matrix,['model']+models)
    # Compatibility alias: generated from the same complete17-model matrix, not hand maintained.
    (OUT/'negative_overlap_matrix_v2.csv').write_bytes((OUT/'negative_overlap_matrix.csv').read_bytes())
    relations=[]
    def rel(k,v,scope='available_full_sequences'):relations.append(dict(relation=k,value=int(v),scope=scope))
    F=total('H1');G=total('A9');B=total('B1');C=F&total('H1',1)
    for k,v in [('F665 total',len(F)),('F966 total',len(G)),('F665 ∩ F966',len(F&G)),('H6 total',len(total('H6'))),('H6 ∩ F665',len(total('H6')&F)),('H4 ∩ F665',len(total('H4')&F)),('Bilingua total',len(B)),('Bilingua benchmark negatives',len(sets['B1','train',0]|sets['B1','test',0])),('Bilingua external negatives',len(sets['B1','external_test',0])),('Bilingua ∩ F665',len(B&F)),('Bilingua unique',len(B-F)),('F665 absent from Bilingua',len(F-B)),('Bilingua 658 ∩ A2 train 658',len((B&F)&sets['A2','train',0])),('F665 conflicting-label sequences',len(C))]:rel(k,v)
    for m in ['A2','A14']:
        for sp in ['train','test']:rel(f'{m} {sp} ∩ F966',len(sets[m,sp,0]&G))
        rel(f'{m} total ∩ F966',len(total(m)&G))
    for m in ['A7','A1','A12','H6']:
        for fn,fs in [('F665',F),('F966',G)]:rel(f'{m} ∩ {fn}',len(total(m)&fs),'repository_only' if m in ['A7','A12'] else 'full_available_version')
    rel('A7 ∩ A2',len(total('A7')&total('A2')),'repository_only');rel('A7 ∩ A1',len(total('A7')&total('A1')),'repository_only')
    for sp in ['train','test']:rel(f'A2 {sp} == A14 {sp}',sets['A2',sp,0]==sets['A14',sp,0])
    rel('A2 total == A14 total',total('A2')==total('A14'));write_csv(OUT/'set_relations.csv',relations)
    split=[]
    for m in ['H1','H2','H3','H5','A3','A4','A5']:
        for sp in ['train','test']:
            split.extend([dict(relation=f'{m} {sp} n',value=len(sets[m,sp,0]),scope='released_negative_members'),dict(relation=f'{m} {sp} == H1 {sp}',value=int(sets[m,sp,0]==sets['H1',sp,0]),scope='released_negative_members')])
    split.append(dict(relation='H4 train == H1 train',value=int(sets['H4','train',0]==sets['H1','train',0]),scope='train_only'))
    for sp in ['train','test']:split.append(dict(relation=f'A2 {sp} == A14 {sp}',value=int(sets['A2',sp,0]==sets['A14',sp,0]),scope='released_negative_members'))
    write_csv(OUT/'split_relations.csv',split)
    roles=[]
    for a in ['train','test']:
        for z in ['development_pool','external_test']:
            for lab in [0,1]:roles.append(dict(upstream='A2',upstream_role=a,downstream='A9',downstream_role=z,label=lab,shared=len(sets['A2',a,lab]&sets['A9',z,lab])))
    write_csv(OUT/'role_relations.csv',roles)
    details=[];loc=[];model_conf=[]
    for m in models:
        cc=total(m)&total(m,1);model_conf.append(dict(model=m,conflicting_sequences=len(cc),scope=by_model[m][0]['dataset']))
        for seq in sorted(cc):
            occurrences=[r for r in by_model[m] if r['sequence']==seq]
            roleset=sorted({f"{r['split']}:{r['label']}" for r in occurrences})
            details.append(dict(sequence=seq,model=m,roles=';'.join(roleset)))
            for r in occurrences:loc.append({k:r[k] for k in ['model','sequence','split','label','source_file','source_record','source_line']})
    write_csv(OUT/'label_conflicts_by_model.csv',model_conf)
    write_csv(ROOT/'05_标签冲突/conflict_occurrences.csv',loc)
    full=['H1','H2','H3','H5','A3','A4','A5']
    write_csv(ROOT/'05_标签冲突/f665_conflict_versions.csv',[r for r in details if r['model'] in full])
    brief=[]
    for seq in sorted(C):
        h=[r for r in by_model['H1'] if r['sequence']==seq]
        roleloc={(r['split'],r['label']):f"{Path(r['source_file']).name}:{r['source_line']}" for r in h}
        cross=len({r['split'] for r in h})>1
        brief.append(dict(sequence=seq,conflict_type='cross_train_test' if cross else 'same_train_opposite_labels',H1_train_pos=roleloc.get(('train',1),''),H1_train_neg=roleloc.get(('train',0),''),H1_test_pos=roleloc.get(('test',1),''),H1_test_neg=roleloc.get(('test',0),''),fixed_split_versions=';'.join(full),H6_status='same_label_conflict;own_random_split',Bilingua_status='absent_both_labels' if seq not in B|total('B1',1) else 'present'))
    write_csv(ROOT/'05_标签冲突/label_conflicts.csv',brief)
    source_rows=sorted(SOURCES.values(),key=lambda r:r['source_file']);write_csv(OUT/'pipeline_source_manifest.csv',source_rows)
    (OUT/'pipeline_run.json').write_text(json.dumps(dict(version=VERSION,models=models,model_count=len(models),record_count=len(RECORDS),source_file_count=len(SOURCES),source_manifest_sha256=hashlib.sha256((OUT/'pipeline_source_manifest.csv').read_bytes()).hexdigest(),normalization='uppercase; remove whitespace; A7 trailing X padding only; preserve valid NA; standard20AA; minimum length2',roles='train/test are published files; A9 development_pool/external_test; B1 external_test separate; H6 model_input_pool; A12 available_only',scope='17 retrievable records; excludes source-linked A6 from direct verification'),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'version':VERSION,'models':len(models),'records':len(RECORDS),'counts':{m:len(total(m)) for m in models}},ensure_ascii=False))

if __name__=='__main__':generate()
