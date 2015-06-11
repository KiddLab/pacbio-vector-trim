import os
import h5py
import sys
import glob
import genutils

from optparse import  OptionParser

###############################################################################
USAGE = """
python trim-vector-pacbio-hdf5.py --blasr <output of blasr mapping to vector> --vectorlen <length of vector>
                                  --original <path to original hd5 pacbio files>
                                  --newdir <directory for new hd5files>

  --blasr is output of blasr mapping of reads to vector backbone sequence.  Map with -m 4 -header
  --bectorlen is length of the vector sequence [default 8139]

Output is in bed format.
"""
parser = OptionParser(USAGE)
parser.add_option('--blasr',dest='blasrFile', help = 'blaster output file')
parser.add_option('--vectorlen',dest='vectorLen',type='int',default=8139, help= 'length of vector sequenc')
parser.add_option('--original',dest='originalPacBio',help= 'prefix to original pacbio files')
parser.add_option('--newdir',dest='newDir',help= 'directory for new pacbio files')

(options, args) = parser.parse_args()

if options.blasrFile is None:
    parser.error('blasr file not given')
if options.vectorLen is None:
    parser.error('vectorLen not given')
if options.originalPacBio is None:
    parser.error('original pacbio prefix not given')
if options.newDir is None:
    parser.error('new output directory not given')

if options.newDir[-1] != '/':
    options.newDir += '/'

###############################################################################

originalBaxFileNames = glob.glob('%s*bax.h5' % options.originalPacBio)
print 'Found %i bax files' % len(originalBaxFileNames)

originalBasFileNames = glob.glob('%s*bas.h5' % options.originalPacBio)
print 'Found %i bas files' % len(originalBasFileNames)

if os.path.isdir(options.newDir) is False:
    print 'Directory %s does not exist' % options.newDir
    print 'Please create it!'
    sys.exit()


# copy over new BAS file
if len(originalBasFileNames) != 1:
    print 'mulitple bas files, do not know what to do'
    sys.exit()
fileName = originalBasFileNames[0]
fileName = fileName.split('/')[-1]

newBasFileName = options.newDir + fileName
print 'New bas file:',newBasFileName
cmd = 'cp %s %s' % (originalBasFileNames[0],newBasFileName)
print cmd
genutils.runCMD(cmd)

#copy over baxFiles
baxFileNames = []
for fileName in originalBaxFileNames:
	newName = fileName.split('/')[-1]
	newName = options.newDir + newName
	baxFileNames.append(newName)
	cmd = 'cp %s %s' % (fileName,newName)
	print cmd
	genutils.runCMD(cmd)
	
#read in vector hits
vectorHits = {}
vectorLen = options.vectorLen

print 'Read in',options.blasrFile
inFile = open(options.blasrFile,'r')
for line in inFile:
	line = line.rstrip()
	line = line.split()
	n = line[0]
	if n == 'qName':
		continue
	parts = n.split('/')
	holeNum = int(parts[1])
	qStart = int(line[5])
	qEnd = int(line[6])
	
	tStart = int(line[9])
	tEnd = int(line[10])

	tHits = [tStart,tEnd]
	tHits.sort()
	# drop if both within 200
	if abs(tHits[0]-0) <=200 and abs(tHits[1]-vectorLen) <= 200:
		if holeNum not in vectorHits:
			vectorHits[holeNum] = []
		vectorHits[holeNum].append(qStart)            
inFile.close()
print 'num with vector hit: %i' % len(vectorHits)


numVectorHits = 0
for hd5FileName in baxFileNames:
    print hd5FileName
    hf = h5py.File(hd5FileName,'r+')
    regions = hf['/PulseData/Regions/']
    numRow = regions.shape[0]
    numCol = regions.shape[1]
    print 'numRow = %i numCol = %i' % (numRow,numCol)
    for row_n in range(numRow):
        row = regions[row_n]
        holeNum = row[0]
        regionTypes = row[1]
        # get ones that are selected vector hits
        if holeNum in vectorHits and regionTypes == 2:
            firstHit = vectorHits[holeNum]
            newEndHighQual = min(firstHit)
            row[3] = newEndHighQual
            if row[2] > row[3]:  # this should not be possible, but check anyway
                row[2] = 0
                row[3] = 0
                row[4] = 0
            regions[row_n] = row
            numVectorHits += 1
    hf.close()
    print 'Num vector hits filtered: %i' % numVectorHits
print 'DONE'
print 'Num vector hits filtered: %i' % numVectorHits



