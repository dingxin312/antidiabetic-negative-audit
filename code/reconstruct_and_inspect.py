"""Read-only check of H6 partition reconstruction and A4 pure data loaders."""
from pathlib import Path
import argparse,csv,json,ast,os
from sklearn.model_selection import train_test_split
p=argparse.ArgumentParser();p.add_argument('project_root');args=p.parse_args();B=Path(args.project_root).resolve()
rows=list(csv.DictReader((B/'H6_SeqDPPIV/source/SeqDPPIV-main/SeqDPPIV model/DPPIV sequence and descriptors.csv').open()))
checks=[]
for seed in range(1,11):
 tr,te=train_test_split(rows,test_size=.2,stratify=[r['label'] for r in rows],random_state=seed)
 trset={r['sequence'] for r in tr};teset={r['sequence'] for r in te}
 checks.append({'seed':seed,'negative_train':sum(r['label']=='0' for r in tr),'negative_test':sum(r['label']=='0' for r in te),'cross_split_exact':sorted(trset&teset)})
source=B/'A4_StructuralDPPIV/source/Structural-DPP-IV-main/util/util_file.py';tree=ast.parse(source.read_text());names={'load_tsv_format_data','load_DPP_IV_data','method_driver'}
funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
for n in funcs:
 n.returns=None
 for a in n.args.args:a.annotation=None
ns={};exec(compile(ast.Module(body=funcs,type_ignores=[]),str(source),'exec'),ns)
cwd=os.getcwd();os.chdir(source.parent.parent/'main')
try:
 v=ns['method_driver']('DPP-IV','valid',None);t=ns['method_driver']('DPP-IV','test',None)
finally:os.chdir(cwd)
print(json.dumps({'H6':checks,'A4':{'validation_records':len(v[0]),'test_records':len(t[0]),'same_ordered_sequences_and_labels':v==t}},indent=2))
