"""Meaningful scientific acceptance checks; does not train or fetch models."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,json,hashlib,subprocess,sys,argparse
from docx import Document
B=Path(__file__).resolve().parents[1];R=B/'07_封版';R.mkdir(exist_ok=True)
def read(rel):return list(csv.DictReader((B/rel).open(encoding='utf-8-sig')))
def run(materials_only=False):
 checks=[]
 def check(name,ok,detail=''):checks.append(dict(check=name,passed=bool(ok),detail=detail))
 records=read('02_解析/negative_sets_long.csv');S=defaultdict(set)
 for r in records:S[r['model'],r['split'],int(r['label'])].add(r['sequence'])
 def total(m,l=0):return set().union(*(v for (mm,sp,lab),v in S.items() if mm==m and lab==l))
 F=total('H1');G=total('A9');C=F&total('H1',1);seven=['H1','H2','H3','H5','A3','A4','A5']
 expected={'A1':236,'A2':1399,'A3':665,'A4':665,'A5':665,'A7':5250,'A9':966,'A12':200,'A13':966,'A14':1399,'H1':665,'H2':665,'H3':665,'H4':532,'H5':665,'H6':629,'B1':666}
 check('17_model_negative_counts',all(len(total(m))==n for m,n in expected.items()) and {r['model'] for r in records}==set(expected))
 for m in seven:
  for sp,n in [('train',532),('test',133)]:check(f'{m}_{sp}_membership_H1',len(S[m,sp,0])==n and S[m,sp,0]==S['H1',sp,0])
 check('H4_train_only',S['H4','train',0]==S['H1','train',0] and not S['H4','test',0])
 check('H6_actual_629_length_filter',total('H6')=={s for s in F if len(s)<=33})
 check('H6_not_upstream_fixed_split',not S['H6','train',0] and not S['H6','test',0] and len(S['H6','model_input_pool',0])==629)
 a=read('H6_SeqDPPIV/source/SeqDPPIV-main/SeqDPPIV model/DPPIV sequence and descriptors.csv');z=read('H6_SeqDPPIV/source/SeqDPPIV-main/SeqDPPIV model/DPPIV descriptors.csv')
 check('H6_features_align_to_sequences',len(a)==len(z)==1293 and all(all(abs(float(x[k])-float(y[k]))<1e-9 for k in y) for x,y in zip(a,z)))
 for m,sp in [('A2','train'),('A2','test')]:check(f'A2_A14_{sp}',S[m,sp,0]==S['A14',sp,0])
 check('A9_658_plus241_pool',len(S['A2','train',0]&S['A9','development_pool',0])==658 and len(S['A2','test',0]&S['A9','development_pool',0])==241 and len(S['A9','development_pool',0])==899)
 check('A9_67_A2test_external',len(S['A9','external_test',0])==67 and S['A9','external_test',0]<=S['A2','test',0])
 check('A2_remainder433',len(total('A2')-G)==433 and G<=total('A2'))
 for lab in [0,1]:
  check(f'A9_pool_external_disjoint_label{lab}',not S['A9','development_pool',lab]&S['A9','external_test',lab])
  check(f'A9_A13_published_roles_label{lab}',S['A9','development_pool',lab]==S['A13','train',lab] and S['A9','external_test',lab]==S['A13','test',lab])
 for m,n in [('H1',51),('A9',101),('A2',120),('A1',205)]:check('A7_depadding_overlap_'+m,len(total('A7')&total(m))==n)
 check('A7_padding_recorded',sum(r['model']=='A7' and r['split']=='train' and r['label']=='0' and r['sequence_raw'].endswith('X') for r in records)==4189)
 check('NA_valid_dipeptide_preserved','NA' in total('H1',1) and all(r['sequence']!='NAN' for r in records))
 expectedC={'FAQTQS','GPFPIIV','LHLPLP','LPVPQ','LTDVEN','LTQTPVVVPPF','VRGPFP'}
 check('seven_conflict_identities',C==expectedC)
 for m in seven:
  check(m+'_same_conflicting_sequences',total(m)&total(m,1)==C)
  check(m+'_three_cross_split_sequences',(S[m,'train',0]|S[m,'train',1])&(S[m,'test',0]|S[m,'test',1])=={'FAQTQS','LTDVEN','VRGPFP'})
 check('H6_retains_seven_conflicts',total('H6')&total('H6',1)==C)
 for m in ['A2','A14','A9','A13']:check(m+'_no_exact_crosslabel_conflict',not total(m)&total(m,1))
 check('B1_negative_roles',tuple(len(S['B1',sp,0]) for sp in ['train','test','external_test'])==(529,132,5))
 check('B1_external49',len(S['B1','external_test',1])==44 and len(S['B1','external_test',0])==5)
 check('B1_removed_exactly_conflicts',F-total('B1')==C and not C&(total('B1')|total('B1',1)))
 check('B1_benchmark_additions3',len((S['B1','train',0]|S['B1','test',0])-F)==3)
 check('two_658_direct_intersection',(total('B1')&F)&S['A2','train',0]=={'DFKLFAVYIKYR','KAFDITYVRLKF'})
 check('A12_available_only',len(S['A12','available_only',0])==200 and not S['A12','train',0])
 matrix={r['model']:r for r in read('02_解析/negative_overlap_matrix.csv')}
 check('matrix_all_pairs_from_canonical_sets',set(matrix)==set(expected) and all(int(matrix[m][n])==len(total(m)&total(n)) for m in expected for n in expected))
 prior=read('07_封版/pre_release_20261102_v1/02_解析/negative_overlap_matrix_v2.csv'); prior_key=list(prior[0])[0]
 check('matches_prior_independently_corrected15model_matrix',all(int(matrix[r[prior_key]][k])==int(v) for r in prior for k,v in r.items() if k!=prior_key))
 check('v2_generated_alias_same_bytes',(B/'02_解析/negative_overlap_matrix.csv').read_bytes()==(B/'02_解析/negative_overlap_matrix_v2.csv').read_bytes())
 # Compare prior independent evidence, not only implementation-derived expectations.
 independent=json.loads((B/'06_独立复核/strict_direction_results.json').read_text())
 check('independent_conflicts_agree',all(sorted(total(m)&total(m,1))==independent['conflicts_by_model'][m] for m in seven))
 check('independent_role_matrix_agrees',all(len(S['A2',sp,0]&S['A9',{'train':'development_pool','test':'external_test'}[z],0])==v for sp,row in independent['A9_role_matrix_from_A2_negatives'].items() for z,v in row.items()))
 master=read('03_Table2/Table2_master.csv');M={r['model']:r for r in master};count=Counter(r['sequence_verified'] for r in master)
 check('master21_status14Y3P4N',len(master)==21 and count=={'Y':14,'P':3,'N':4})
 check('overview_matches_master',all(r['sequence_verified']==M[r['model']]['sequence_verified'] and r['n_neg_total']==M[r['model']]['n_neg_total'] for r in read('00_索引/数据总览.csv')))
 check('H6_Y_subset_OWN_and_overlap',M['H6']['sequence_verified']=='Y' and M['H6']['split_identity']=='OWN' and M['H6']['overlap_F665']=='629/665' and M['H6']['overlap_F966']=='2/966')
 check('partial_total_not_available_count',all(M[m]['n_neg_total']=='' for m in ['H4','A7','A12']))
 check('N_reported_not_retrieved',M['A10']['n_neg_total']=='' and M['A10']['n_neg_reported']=='460' and M['A6']['citation_linked']=='Y')
 check('B1_count_fields_separate',M['B1']['n_neg_train']=='529' and M['B1']['n_neg_test']=='132' and M['B1']['n_neg_external']=='5')
 check('A9_final_train_unknown',M['A9']['n_neg_train']=='' and M['A9']['n_neg_development_pool']=='899' and M['A9']['n_neg_external']=='67')
 check('formal_overlaps_only_complete',all(not r['overlap_F665'] and not r['overlap_F966'] for r in master if r['sequence_verified']!='Y'))
 paper=read('03_Table2/Table2_paper.csv');family=read('03_Table2/Table2_family_master.csv')
 check('family_and_paper9rows',len(paper)==len(family)==9)
 if not materials_only:
  doc=Document(B/'03_Table2/Table2_paper.docx');check('docx9data_rows',len(doc.tables)==1 and len(doc.tables[0].rows)==10)
  def clean(s):return ''.join(s.split())
  check('docx_csv_model_alignment',all(clean(doc.tables[0].rows[i+1].cells[1].text)==clean(r['models']) for i,r in enumerate(paper)))
  check('docx_csv_status_alignment',all(doc.tables[0].rows[i+1].cells[2].text==r['sequence_verified'] for i,r in enumerate(paper)))
 check('no_false_F966_seven_conflicts',all('7 label conflicts retained' not in r['notes'] and '7 label conflicts retained' not in r.get('source_description','') for r in family if r['family_id']=='F966'))
 if not materials_only:
  fig=json.loads((B/'03_Table2/figure_sources.json').read_text());check('figure_cells_equal_canonical',all(int(matrix[a][c])==n for a,rr in fig['displayed_cells'].items() for c,n in rr.items()))
  check('figures_have21_model_nodes',set(fig['figure2_model_order'])==set(M))
 check('conflict_summary7_detail49',len(read('05_标签冲突/label_conflicts.csv'))==7 and len(read('05_标签冲突/f665_conflict_versions.csv'))==49)
 baseline=json.loads((B/'07_封版/pre_release_20261102_v1/raw_source_baseline.json').read_text());bad=[]
 for r in baseline:
  p=B/r['path']
  if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:bad.append(r['path'])
 check('all597_raw_source_files_unchanged',not bad,','.join(bad))
 result=dict(version='evidence-20261102-v1',scope='materials_only' if materials_only else 'materials_and_existing_displays',passed=sum(c['passed'] for c in checks),failed=sum(not c['passed'] for c in checks),checks=checks,limits=['Published data/implementation scope only','No performance retraining/effect estimate','H4test/A7paper mapping/A12full remain unresolved','A4fresh online code check previously failed; local saved source scope','No first-discovery or field prevalence claim'])
 (R/'validation_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
 print(json.dumps({'passed':result['passed'],'failed':result['failed'],'failures':[c for c in checks if not c['passed']]},ensure_ascii=False))
 if result['failed']:raise SystemExit(1)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--materials-only',action='store_true');args=parser.parse_args();run(materials_only=args.materials_only)
