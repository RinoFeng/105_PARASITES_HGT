from pathlib import Path
import pandas as pd,shutil,json,hashlib
from PIL import Image
from zipfile import ZipFile
O=Path(__file__).resolve().parents[1];R=O.parent;S=R/'manuscript_author_comments_20261004'
for folder in ['data','figures','tables','provenance']:shutil.copytree(S/folder,O/folder,dirs_exist_ok=True)
(O/'pdf').mkdir(exist_ok=True)
for f in list(S.glob('*.docx'))+list(S.glob('*.md')):shutil.copy2(f,O/f.name)
# Audited current data; do not rerun biological filtering.
T=R/'targeted_manuscript_revision_20260921/tables';d=pd.read_csv(T/'plotting_catalogue_metadata_8989.tsv',sep='\t');d=d[d.cohort_role!='certified_supplementary_candidate'];assert len(d)==8985;d.to_csv(O/'data/Figure2B_existing_records.tsv',sep='\t',index=False)
a=pd.read_csv(T/'figure2_abundance.tsv',sep='\t');assert a.records.sum()==8989;a.to_csv(O/'data/Figure2A_existing_counts.tsv',sep='\t',index=False)
shutil.copy2(R/'manuscript_refinement_20260923/tables/figure2_components48.tsv',O/'data/Figure2CD_existing_components.tsv')
# Lower-level lineage labels may be inherited only for the unchanged reference accession.
c=pd.read_excel(S/'tables/Table_S1_Catalogue.xlsx',sheet_name='Catalogue_8989')
labels=[]
for _,r in c.iterrows():
 domain=str(r.display_best_distant_domain);domain='Viruses' if domain=='VIRAL' else domain
 same=str(r.display_best_distant_accession).strip()==str(r.Target_accesson).strip()
 lineage=str(r.Target_lineage);reviewed=r.uniform_search_reviewed is True or r.uniform_search_reviewed==True
 if domain=='Viruses':name='Viruses';evidence='current_domain_only'
 elif reviewed and not same:name='Other '+domain;evidence='updated_reference;prior_lower_taxonomy_not_reused'
 else:
  terms=[x.strip() for x in lineage.split(';')]
  if domain=='Bacteria':
   names=['Pseudomonadota','Bacillota','Bacteroidota','Actinomycetota','Cyanobacteria','Cyanobacteriota','Spirochaetota','Chlamydiota','Verrucomicrobiota','Planctomycetota','Mycoplasmatota','Proteobacteria','Firmicutes','Bacteroidetes','Actinobacteria','Fusobacteriota','Acidobacteriota','Deinococcota']
  elif domain=='Eukaryota':names=['Ascomycota','Basidiomycota','Streptophyta','Chlorophyta','Chordata','Arthropoda','Nematoda','Apicomplexa','Ciliophora','Mollusca','Bacillariophyta','Oomycota','Amoebozoa','Euglenozoa','Discosea']
  elif domain=='Archaea':names=['Euryarchaeota','Thermoproteota','Crenarchaeota','Nitrososphaerota','Thaumarchaeota','Thermoplasmatota']
  else:raise ValueError(domain)
  found=[x for x in terms if x in names];name=found[0] if found else 'Other '+domain;evidence='preserved_reference_lineage'
 labels.append(dict(hgt_id=r.hgt_id,domain=domain,lineage_group=name,basis=evidence,current_reference=r.display_best_distant_accession))
l=pd.DataFrame(labels);assert len(l)==8989;l.to_csv(O/'data/Figure3A_ring_record_mapping.tsv',sep='\t',index=False)
v=l.groupby(['domain','lineage_group']).size().reset_index(name='n');out=[]
for domain,x in v.groupby('domain',sort=False):
 x=x.sort_values('n',ascending=False)
 for j,r in enumerate(x.itertuples()):out.append(dict(domain=domain,lineage_group=r.lineage_group if r.n/8989>=.005 and j<6 else 'Other '+domain,n=r.n))
v=pd.DataFrame(out).groupby(['domain','lineage_group'],sort=False).n.sum().reset_index();v.to_csv(O/'data/Figure3A_ring_counts.tsv',sep='\t',index=False)
assert v.groupby('domain').n.sum().to_dict()=={'Archaea':46,'Bacteria':7758,'Eukaryota':1131,'Viruses':54}
# Preserve current original panels byte-for-byte unless explicitly restored.
rows=json.loads((O/'checks/figure_inventory.json').read_text());audit=[]
pairs=[('current_main',2,'original_main',2),('current_main',5,'original_main',5),('current_SI',2,'original_SI',2),('current_SI',7,'original_SI',7),('current_SI',8,'original_SI',8)]
for ck,ci,ok,oi in pairs:
 import numpy as np
 x=np.array(Image.open(O/'reference_images'/f'{ck}_{ci}.png'));y=np.array(Image.open(O/'reference_images'/f'{ok}_{oi}.png'));audit.append(dict(current=f'{ck}_{ci}',original=f'{ok}_{oi}',pixel_identical=bool(x.shape==y.shape and np.array_equal(x,y))))
(O/'checks/unchanged_original_image_comparison.json').write_text(json.dumps(audit,indent=2))
print('source counts',a.to_dict('records'));print('preserved source images',audit)
