from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parents[1];S=O.parent/'manuscript_author_comments_20261004'
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
   names=['Pseudomonadota','Bacillota','Bacteroidota','Actinomycetota','Cyanobacteria','Cyanobacteriota','Spirochaetota','Chlamydiota','Verrucomicrobiota','Planctomycetota','Mycoplasmatota','Proteobacteria','Firmicutes','Bacteroidetes','Actinobacteria','Fusobacteriota','Acidobacteriota','Deinococcota','Chloroflexi','Chloroflexota','Acidobacteria','Chlamydiae','Verrucomicrobia','Spirochaetes','Gemmatimonadetes','Gemmatimonadota','Armatimonadetes','Deinococcus-Thermus','Fusobacteria','Elusimicrobia','Fibrobacteres','Nitrospirae','Nitrospirota','Ignavibacteriae','Balneolaeota']
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
