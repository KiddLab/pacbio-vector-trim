# pacbio-vector-trim

This repository contains scripts for preprocessing of Pac Bio reads for clone assembly.

The problem is that multiple clones pooled into a single SMRT cell will all have the same
vector backbone.  This can create problems in the assembly.  There are several descriptions
of how to deal with this -- here is what we came up with that seems to work.


## Required Components
Python: h5py, genutils
Other: BLASR, access to SMRTPORTAL or other installation of assembly pipeline


## Description
The process begins with mapping reads to the vector backbone sequence.  The included
sequence, pCC1FOS.cloning.rc.fa, is the backbone for the pCC1 fosmid vector altered based
on the position of the standard cloning site.

```
blasr  \
m150606_042117_42131_c100805572550000001823179010081530_s1_p0.bas.h5 \
pCC1FOS.cloning.rc.fa \
-m 4 -header  -out m150606_042117_42131.pCC1FOS.blasr.out
```


Next, a new version (with same name, but please put in a NEW DIRECTORY) is made of
the original PacBio hdf5 files, but reads that span across the vector backbone are removed.
These are identified by reads with matches to within 200 bp of the start/end of the vector and
is done by editing the 'Regions' listing in the .bax files to edit high quality region
definitions.




```
python trim-vector-pacbio-hdf5.py \
--blasr m150606_042117_42131.pCC1FOS.blasr.out \
--original m150606_042117_42131_c100805572550000001823179010081530_s1_p0 \
--newdir filtered/ 
```

The resulting files can then be assembled in the usual way.  For pools of fosmid clones
we currently use HGAP.v2, being sure to set genome size appropriately based on the number of
clones in the pool. 


For more information, see

http://files.pacb.com/software/instrument/2.0.0/bas.h5%20Reference%20Guide.pdf

and discussions at
http://seqanswers.com/forums/archive/index.php/t-30439.html


http://www.freelists.org/post/mira_talk/BAC-vector-sequece-masking-for-de-novo-assembly-using-PacBio-C2,7


There are probably better ways of doing this, but this seems to work in our hands.

