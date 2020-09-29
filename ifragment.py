#!/usr/bin/env python
from __future__ import print_function
try:
	from ifragment_functions import *
except ImportError as e:
	sys.exit("Couldn't load \"ifragment_functions\". Check the files were downloaded correctly.")

if (len(sys.argv)<2):
	entid='HMDB0000252'
	print("WARNING! Running without arguments, using default values...")
else:entid=sys.argv[1]
if entid.startswith("c"):db="envipath"
elif entid.upper().startswith("CHEBI"):db="reactome"
elif entid.startswith("HMDB"):db="smpdb"
else:db=""
if db=='': sys.exit("Not valid compound ID, check the files in '{}' folder".format("entities_list/"))

nodef,df,rf=False,os.path.join(os.getcwd(),"ifragment_datafiles"),os.path.join(os.path.dirname(sys.argv[0]),"results")
fpfile,dbfile,null_model_file="{0}_fingerprints.csv {0}.tsv {0}_NullModel.tsv".format(os.path.join(df,db)).split()
results=os.path.join(rf,"scores_"+entid+".tsv")

#Load necesary files
null_model=load_nullmodel(null_model_file)
database=load_db(dbfile)
fingerprints,allents,nents,fpsize=load_fp(fpfile)
if not len(fingerprints)==nents:sys.exit("Entity count mismatching")
#Check saving path
if not os.path.isdir(rf):
	try:os.mkdir(rf)
	except OSError:
		nodef=True
		results=os.path.dirname(sys.argv[0])
		print('WARNING! Couldn\'t create folder "{}".\nSaving results into: "{}"'.format(rf,results),file=sys.stderr)
if nodef:print('Results will be saved into: "{}"'.format(results))
try:
	with open(results,"w") as fh:pass
except OSError as e:
	sys.exit('ERROR! Results file can\'t be saved into: "{}" check that you have permission to write in that folder.'.format(os.path.basename(results)))

inbg=np.bincount(np.concatenate(fingerprints),minlength=fpsize)
#Fingerprint of the compound selected
fp=fingerprints[(allents==entid).nonzero()[0]][0]

with open(results,"w") as fh:
	for (pathwayid,name),entities in database.items():
		ind=np.in1d(allents,entities).astype(int).nonzero()[0]
		pv1,pvind=fp_pv(fingerprints[ind],nents,inbg,fpsize)
		common_fpind,common_pvind=np.intersect1d(fp,pvind,assume_unique=True,return_indices=True)[1:]
		sc=0. if len(common_pvind)==0 else -np.log10(pv1[common_pvind]).sum()/len(fp)
		zr_mean,zr_std,location,scale=null_model[pathwayid]
		zsc=(sc-zr_mean)/zr_std
		pv=genextreme.sf(zsc,0,location,scale)
		print(pathwayid,name,entid,sc,zsc,pv,file=fh)
print("Saved",results.replace(os.getcwd(),""))
