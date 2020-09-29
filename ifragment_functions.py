import sys
import os
#grepf <(grep CSBG entities_list/smpdb.tsv|sort -u) /tmp/chemprof/smpdb27ac4col_conversion !!
try:
	import numpy as np
except ImportError as e:
	sys.exit("Couldn't load mandatory module: 'numpy'.")
try:
	from scipy.stats import hypergeom,genextreme
except ImportError as e:
	sys.exit("Couldn't load mandatory module: 'scipy'.\nInstructions for installing can be found at: https://scipy.org/install.html")

def fp_pv(fpset,nents,inbg,fpsize,bon=50.0,**kwargs):
	frset,inset=np.unique(np.concatenate(fpset),return_counts=True)
	pv=hypergeom.sf(inset-1,nents,inbg[frset],len(fpset))
	nozpv=pv.nonzero()[0]
	if not len(pv)==len(nozpv):
		m=pv[nozpv].min()
		minval=np.nextafter(0,m) if m/bon == 0 else m/bon
		pv[pv==0]=minval
	return(pv,frset)

def load_fp(fpf):
	fingerprints,allents,fpsize=[],[],0
	try:  
		with open(fpf,"r") as fh:
			for line in fh:
				l=line.strip().split(",")
				fingerprints.append(np.array(l[1:],dtype=int))
				allents.append(l[0])
	except OSError as e:sys.exit('ERROR! Couldn\'t open "{}".\nCheck that the file exists and you can open it.'.format(fpf))
	fpsize=max(map(int,[fr for fp in fingerprints for fr in fp]))+1
	return(np.array(fingerprints,dtype=object),np.array(allents),len(allents),fpsize)

def load_db(dbf):
	dbd={}
	try:
		with open(dbf,"r") as fh:
			for line in fh:
				pathwayid,entid,name=line.strip().split("\t")
				if (pathwayid,name) in dbd: dbd[(pathwayid,name)].append(entid)
				else: dbd[(pathwayid,name)]=[entid]
	except OSError as e:sys.exit('ERROR! Couldn\'t open "{}".\nCheck that the file exists and you can open it.'.format(dbf))
	return(dbd)

def load_nullmodel(nmf):
	try:
		with open(nmf,"r") as fh:
			nm=[line.strip().split("\t") for line in fh]
	except OSError as e:sys.exit('ERROR! Couldn\'t open "{}".\nCheck that the file exists and you can open it.'.format(nmf))  
	return(dict([(row[0],map(float,row[1:])) for row in nm]))
