import mstpython

pdbfile = "testfiles/2ZTA.pdb"
outPdbFile = "/tmp/out.pdb"

S = mstpython.Structure(pdbfile);
print "read %s, has %d chains, %d residues, and %d atoms" % (pdbfile, S.chainSize(), S.residueSize(), S.atomSize())

for i in range(S.chainSize()):
  C = S.getChain(i)
  print "\tchain %d, ID = %s, seg ID = %s, has %d residues and %d atoms" % (i+1, C.getID(), C.getSegID(), C.residueSize(), C.atomSize())
  newSegID = "S%d" % i
  print "\tchanging segid to %s" % newSegID
  C.setSegID(newSegID)

print "writing to %s.." % outPdbFile
S.writePDB(outPdbFile, "")
