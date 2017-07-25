import mstpython

pdbfile = "testfiles/2ZTA.pdb"
outPdbFile = "/tmp/out.pdb"

S = mstpython.Structure(pdbfile);
print "read %s, has %d chains, %d residues, and %d atoms" % (pdbfile, S.chainSize(), S.residueSize(), S.atomSize())

for i in range(S.chainSize()):
  C = S.getChain(i)
  print "\tchain %d, ID = %s, seg ID = %s, has %d residues and %d atoms" % (i+1, C.getID(), C.getSegID(), C.residueSize(), C.atomSize())
  for j in range(C.residueSize()):
    R = C.getResidue(j)
    print "\t\tresidue %d, name = %s, has %d atoms" % (i+1, R.getName(), R.atomSize())
    for k in range(R.atomSize()):
      A = R.getAtom(k)
      print "\t\t\tatom %d, name = %s, coor = [%8.3f %8.3f %8.3f]" % (i+1, A.getName(), A.getX(),  A.getY(),  A.getZ())

  newSegID = "S%d" % i
  print "\tchanging segid to %s" % newSegID
  C.setSegID(newSegID)

print "writing to %s.." % outPdbFile
S.writePDB(outPdbFile, "")
