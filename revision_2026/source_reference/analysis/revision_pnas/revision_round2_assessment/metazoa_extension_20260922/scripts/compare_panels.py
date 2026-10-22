from pathlib import Path
from collections import defaultdict
import hashlib,json,bisect
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[1];R=B.parent;F=R/'innovation_final8989_20260915';D=R/'family_contribution_extension_20260922';T=B/'tables';P=B/'figures';CL=['Amoebozoa','Apicomplexa','Euglenozoa','Metazoa'];CATS=list('EFGH');PARTS=CATS+['Other','Unknown'];checks=[]
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(n,d):pd.DataFrame(d).to_csv(T/n,sep='\t',index=False,compression='infer')
def check(n,x):
 checks.append(dict(check=n,passed=bool(x)))
 if not x:raise AssertionError(n)
def truth(x):return str(x).lower() in {'true','1','yes'}
def abin(x):return '0' if x==0 else '1' if x==1 else 'low' if x<=.5 else 'high'
def key(r):return tuple(abin(r[c]) for c in ['eggnog_fraction','informative_cog_fraction','ko_fraction','pfam_fraction'])
def match(dat):
 rows=[];elig=[];pairs=[]
 for sp,z in dat.groupby('species'):
  pools=defaultdict(list)
  for b in z[z.role=='background'].to_dict('records'):pools[key(b)].append(b)
  for k in pools:pools[k].sort(key=lambda r:(r['length'],r['OG']))
  lens={k:[r['length'] for r in v] for k,v in pools.items()}
  for c in z[z.role=='candidate'].to_dict('records'):
   k=key(c);ll=lens.get(k,[]);pp=pools.get(k,[]);el=pp[bisect.bisect_left(ll,.8*c['length']):bisect.bisect_right(ll,1.25*c['length'])];cs=sorted(el,key=lambda b:(abs(np.log(b['length']/c['length'])),hashlib.sha256((c['OG']+'|'+b['OG']).encode()).hexdigest()))[:5]
   elig.append(dict(species=sp,clade=c['clade'],OG=c['OG'],eligible_controls=len(el),matched=len(cs)>=3))
   if len(cs)<3:continue
   for b in cs:pairs.append(dict(species=sp,candidate_OG=c['OG'],background_OG=b['OG'],length_ratio=b['length']/c['length']))
   for cat in PARTS:
    bg=np.mean([b[cat] for b in cs]);rows.append(dict(species=sp,clade=c['clade'],OG=c['OG'],category=cat,candidate_share=c[cat],background_share=bg,difference=c[cat]-bg,candidate_loci=c['units'],matched_controls=len(cs)))
 return pd.DataFrame(rows),pd.DataFrame(elig),pd.DataFrame(pairs)
old=read(D/'tables/family_profiles47.tsv.gz');old['role']=old.role.replace({'ordinary_background':'background'});oldscope=set(old.species);check('47_original_species',len(oldscope)==47)
u=read(T/'gene_units_recovered.tsv.gz',keep_default_na=False);newblock=set(read(T/'new_historical_blocked_OGs.tsv').OG);allblock=set(read(T/'all_historical_blocked_OGs.tsv').OG)
# Coarse labels only; original functional annotations reused exactly via the sequence hash.
vv=[]
for z in u.itertuples():
 ll=set(z.COG_category)&set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')-set('RS');v=[1/len(ll) if c in ll else 0 for c in CATS] if ll else [0.]*4;v+=[sum(c not in CATS for c in ll)/len(ll) if ll else 0,0 if ll else 1];vv.append(v)
u[PARTS]=np.asarray(vv);u['eggnog_fraction']=u.eggnog_present.map(truth).astype(float);u['ko_fraction']=u.KO_present.map(truth).astype(float);u['pfam_fraction']=u.PFAM_present.map(truth).astype(float);u['informative_cog_fraction']=1-u.Unknown;u['length']=u.aa_length_measured.astype(float)
uf=u[u.role.isin(['candidate','ordinary_background'])].copy();uf['role']=uf.role.replace({'ordinary_background':'background'});gp=uf.groupby(['species','reference_assembly','same_original_reference','OG','role'],as_index=False);new=gp[PARTS+['eggnog_fraction','ko_fraction','pfam_fraction','informative_cog_fraction']].mean();aux=gp.agg(units=('unit_id','nunique'),length=('length','median'));new=new.merge(aux,on=['species','reference_assembly','same_original_reference','OG','role'],validate='one_to_one');new['clade']='Metazoa';write('recovered_family_profiles.tsv.gz',new)
check('new_profile_sums',np.allclose(new[PARTS].sum(axis=1),1));check('no_historical_background',not new.loc[new.role=='background','OG'].isin(allblock).any())
# Descriptive quantity at each recovered reference, without merging assemblies.
rs=read(T/'species_recovery_summary.tsv');aa=read(T/'candidate_locus_audit.tsv');recover=[]
for z in rs.itertuples():
 q=new[(new.species==z.species)&(new.reference_assembly==z.reference_assembly)];c=q[q.role=='candidate'];bg=q[q.role=='background'];eff,el,pp=match(q);n=eff.OG.nunique() if len(eff) else 0
 recover.append(dict(species=z.species,reference_assembly=z.reference_assembly,same_original_reference=z.same_original_reference,final_records=z.final_records,full_CDS_strict_records=z.strict_final_records,candidate_loci=c.units.sum(),candidate_OGs=len(c),background_loci=bg.units.sum(),background_OGs=len(bg),fine_matched_candidate_OGs=n,eligible_for_augmented_panel=truth(z.same_original_reference) and n>=10,reason='eligible' if truth(z.same_original_reference) and n>=10 else 'alternative_reference_bridge_only' if not truth(z.same_original_reference) else 'fewer_than10_matched_families'))
recover=pd.DataFrame(recover);write('resource_to_panel_admission.tsv',recover)
# Original profiles reproduce the prior fine matching before any new exclusion is imposed.
baseeff,baseel,basepairs=match(old);saved=read(F/'ecology/matched_family_function_effects.tsv.gz');saved=saved[saved.species.isin(oldscope)];j=baseeff[baseeff.category.isin(CATS)].merge(saved,on=['species','OG','category'],suffixes=('_rebuilt','_saved'),validate='one_to_one');check('original_matching_row_count',len(j)==len(saved));check('original_candidate_profiles',np.allclose(j.candidate_share_rebuilt,j.candidate_share_saved));check('original_control_profiles',np.allclose(j.background_share_rebuilt,j.background_share_saved))
quarantine=old[~((old.role=='background')&old.OG.isin(newblock))].copy();added_species=set(recover[recover.eligible_for_augmented_panel].species)-oldscope;same=new[new.species.isin(added_species)&new.same_original_reference.map(truth)].copy();ext=pd.concat([quarantine,same],ignore_index=True)
scenarios={'original47':old,'original47_updated_quarantine':quarantine,'augmented_same_reference':ext};effrows=[];adm=[];qrows=[];allpairs=[];fullcounts=[];cladeeffects=[]
for scen,dat in scenarios.items():
 eff,el,pp=(baseeff,baseel,basepairs) if scen=='original47' else match(dat);n=eff[eff.category=='E'].groupby('species').OG.nunique();eligible_species=set(n[n>=10].index);valid=eff[eff.species.isin(eligible_species)];se=valid.groupby(['species','clade','category'],as_index=False).agg(candidate_share=('candidate_share','mean'),background_share=('background_share','mean'),effect=('difference','mean'),matched_families=('OG','nunique'));se['scenario']=scen;effrows.append(se);pp['scenario']=scen;allpairs.append(pp)
 for sp,z in dat.groupby('species'):
  c=z[z.role=='candidate'];b=z[z.role=='background'];num=int(n.get(sp,0));adm.append(dict(scenario=scen,species=sp,clade=z.clade.iloc[0],matched_families=num,in_comparison=sp in eligible_species));cc=len(c);bb=len(b);cn=c.units.sum();bn=b.units.sum();row=dict(scenario=scen,species=sp,clade=z.clade.iloc[0],C=cc,B=bb,N=cn,M=bn,family_ratio=cc/bb if bb else np.nan,locus_ratio=cn/bn if bn else np.nan,candidate_loci_per_family=cn/cc if cc else np.nan,background_loci_per_family=bn/bb if bb else np.nan,matched_families=num,entered_function_comparison=sp in eligible_species);qrows.append(row)
  if cc and bb:check('two_level_identity_'+scen+sp,np.isclose(row['family_ratio']*(cn/cc)/(bn/bb),row['locus_ratio']))
  for zz in c.itertuples():fullcounts.append(dict(scenario=scen,species=sp,clade=z.clade.iloc[0],role='candidate',OG=zz.OG,loci=zz.units))
  for zz in b.itertuples():fullcounts.append(dict(scenario=scen,species=sp,clade=z.clade.iloc[0],role='background',OG=zz.OG,loci=zz.units))
 ce=se.groupby(['clade','category'],as_index=False).agg(species=('species','nunique'),candidate_share=('candidate_share','mean'),background_share=('background_share','mean'),effect=('effect','mean'));ce['scenario']=scen;cladeeffects.append(ce)
 for cat in PARTS:
  z=ce[(ce.category==cat)&ce.clade.isin(CL)];cladeeffects.append(pd.DataFrame([dict(scenario=scen,clade='FOUR_CLADE_EQUAL',category=cat,species=z.species.sum(),candidate_share=z.candidate_share.mean(),background_share=z.background_share.mean(),effect=z.effect.mean())]))
se=pd.concat(effrows,ignore_index=True);ce=pd.concat(cladeeffects,ignore_index=True);qr=pd.DataFrame(qrows);write('species_function_effects.tsv',se);write('clade_function_effects.tsv',ce);write('matched_pairs.tsv.gz',pd.concat(allpairs,ignore_index=True));write('panel_admission.tsv',adm);write('species_quantity.tsv',qr);write('species_family_counts.tsv.gz',fullcounts)
# Fixed n=3 species draws, identical method for all panels. Spread is not a CI or p-value.
rng=np.random.default_rng(20260922);draws=[];contrasts=[]
for scen,dat in qr.groupby('scenario'):
 rng=np.random.default_rng(20260922)  # paired seed: unchanged panels have identical draws
 dat=dat[dat.entered_function_comparison];draw_by_clade={}
 for cl in CL:
  z=dat[dat.clade==cl].sort_values('species').reset_index(drop=True);check('three_species_'+scen+cl,len(z)>=3);arr=[]
  for rep in range(1000):
   take=z.iloc[rng.choice(len(z),3,replace=False)];c,b,n,m=take[['C','B','N','M']].sum();fr=c/b;lr=n/m;vr=(n/c)/(m/b);check('draw_identity_'+scen+cl+str(rep),np.isclose(fr*vr,lr));arr.append((fr,vr,lr));draws.append(dict(scenario=scen,clade=cl,replicate=rep,species=';'.join(take.species),family_ratio=fr,within_family_loci_ratio=vr,locus_ratio=lr))
  draw_by_clade[cl]=np.asarray(arr)
 for cl in CL[:-1]:
  rat=draw_by_clade[cl]/draw_by_clade['Metazoa'];check('contrast_identity_'+scen+cl,np.allclose(rat[:,0]*rat[:,1],rat[:,2]))
  for k,label in enumerate(['relative_family_abundance','relative_loci_per_family','relative_locus_abundance']):
   vals=rat[:,k];contrasts.append(dict(scenario=scen,clade=cl,reference='sampled_nematode_species',component=label,geometric_mean=np.exp(np.log(vals).mean()),q025=np.quantile(vals,.025),q975=np.quantile(vals,.975),draws=1000,species_per_draw=3,interpretation='Central95% species-subset spread, not CI or phylogenetic P-value'))
write('equal_species_draws.tsv.gz',draws);ct=pd.DataFrame(contrasts);write('quantity_contrasts.tsv',ct)
# Full catalogue organisation audit; whole-catalogue numbers never mix with strict-locus counts.
coverage=read(F/'macro/species_coverage_105.tsv');coverage['organisational_comparison_group']=np.where(coverage.clade=='Metazoa','multicellular_metazoan_parasites','sampled_unicellular_parasites');coverage['literal_comparison_included']=~coverage.species.isin(['Aphanomyces astaci','Pythium insidiosum','Chromera velia']);coverage['scope_note']='species-level catalogue observation, not true event/entry/retention rate';coverage.loc[coverage.species.isin(['Aphanomyces astaci','Pythium insidiosum']),'scope_note']='oomycete organisation outside this literal two-group contrast';coverage.loc[coverage.species=='Chromera velia','scope_note']='free-living reference; excluded from literal parasite contrast';write('catalogue105_organisational_scope.tsv',coverage)
catrows=[]
for label,z in [('literal102',coverage[coverage.literal_comparison_included]),('archived105',coverage)]:
 for group,q in z.groupby('organisational_comparison_group'):
  catrows.append(dict(scope=label,group=group if label=='literal102' else 'Metazoa' if group.startswith('multicellular') else 'non_Metazoa_not_all_unicellular_parasites',species=len(q),catalogue_records=q.catalogue_records.sum(),mean_records=q.catalogue_records.mean(),median_records=q.catalogue_records.median(),minimum_records=q.catalogue_records.min(),maximum_records=q.catalogue_records.max()))
write('catalogue_group_descriptions.tsv',catrows)
# Figures: all panels label the limited taxonomic range.
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
def save(fig,n):
 for ex in ['pdf','svg','png']:fig.savefig(P/(n+'.'+ex),dpi=240)
 plt.close(fig)
meta_n=se[(se.scenario=='augmented_same_reference')&(se.clade=='Metazoa')].species.nunique();fig,axs=plt.subplots(1,3,figsize=(165/25.4,83/25.4),layout='constrained')
for ax,cl in zip(axs,CL[:-1]):
 for i,sc in enumerate(['original47','augmented_same_reference']):
  z=ct[(ct.scenario==sc)&(ct.clade==cl)].set_index('component').loc[['relative_family_abundance','relative_loci_per_family']];y=np.log2(z.geometric_mean.to_numpy());lo=np.log2(z.q025.to_numpy());hi=np.log2(z.q975.to_numpy());ax.errorbar(np.arange(2)+(i-.5)*.12,y,yerr=[y-lo,hi-y],fmt='o',ms=4,capsize=2,label='Original panel' if i==0 else 'Augmented panel')
 ax.axhline(0,c='0.6',ls=':',lw=.8);ax.set_xticks([0,1],['Family\nabundance','Loci per\nfamily']);ax.set_title(cl);ax.spines[['top','right']].set_visible(False)
axs[0].set_ylabel('Contrast with sampled nematodes (log2)');axs[1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.22));save(fig,'Quantity_comparison_before_after')
fig,axs=plt.subplots(1,2,figsize=(165/25.4,85/25.4),layout='constrained');colors=['#4477AA','#EE6677','#228833','#AA3377']
for j,sc in enumerate(['original47','augmented_same_reference']):
 z=ce[(ce.scenario==sc)&(ce.clade=='Metazoa')].set_index('category').loc[CATS];axs[0].plot(np.arange(4)+(j-.5)*.08,z.effect*100,'o-',label='Original 5 nematodes' if j==0 else f'Augmented {meta_n} nematodes',ms=4)
axs[0].set_xticks(np.arange(4),CATS);axs[0].axhline(0,c='0.6',ls=':',lw=.8);axs[0].set_ylabel('Candidate minus matched background (pp)');axs[0].set_title('A  Nematode functional contrasts',loc='left');axs[0].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.12))
for cl,col in zip(CL,colors):
 z=ce[(ce.scenario=='augmented_same_reference')&(ce.clade==cl)].set_index('category').loc[CATS];axs[1].plot(np.arange(4),z.effect*100,'o-',label=cl if cl!='Metazoa' else f'Nematodes (n={meta_n})',color=col,ms=4)
axs[1].set_xticks(np.arange(4),CATS);axs[1].axhline(0,c='0.6',ls=':',lw=.8);axs[1].set_ylabel('Candidate minus matched background (pp)');axs[1].set_title('B  Separate lineage comparisons',loc='left');axs[1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2)
for ax in axs:ax.spines[['top','right']].set_visible(False)
save(fig,'Functional_comparison_before_after');write('comparison_checks.tsv',checks)
status=dict(status='complete',original_panel_species=47,added_same_reference_species=sorted(added_species),augmented_panel_species=se[se.scenario=='augmented_same_reference'].species.nunique(),augmented_metazoan_species=meta_n,metazoan_scope='nematodes only',strict_new_candidate_records=int(recover.full_CDS_strict_records.sum()),new_catalogue_records=0,quarantine_added_OGs=len(newblock),original_matched_results_reproduced=True,formal_pvalues=False,multicellularity_independent_effect_identified=False,checks=len(checks),all_checks_passed=all(c['passed'] for c in checks))
(B/'COMPARISON_STATUS.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status));print(recover.to_string(index=False));print(ce[(ce.scenario=='augmented_same_reference')&(ce.clade=='Metazoa')].to_string(index=False))
