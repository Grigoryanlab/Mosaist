use GENERAL;
use PDB;
use CHARMM;
use DEFINITIONS;
use SEQUENCE;

if (scalar(@ARGV) != 3) {
 die "\n<exec> [rotLibFile.dat] [chi_definition] [outputRotamerLibrary.bin]\n\n";
}

my $rotLibFile = shift;
my $chiDefFile = shift;
my $outRotLibFile = shift;
my $chidef = readChiDefinition($chiDefFile);
my $tmpBase = "/tmp/$$";

# read Dunbrack rotamer library
my $R = readDunbrackFormattedRotamerLibrary($rotLibFile);

# create entries for ALA and GLY (naturally, not listed in library)
my %b = ("phi" => 0, "psi" => 0);
my %r = ("chi" => GENERAL::arrayRef(), "chiSig" => GENERAL::arrayRef());
foreach my $aa ("ALA", "GLY") { push(@{$R->{$aa}{bins}}, \%b); }
foreach my $aa ("ALA", "GLY") { push(@{$R->{$aa}{bin}{$b{phi}}{$b{psi}}}, \%r); }
foreach my $aa ("ALA", "GLY") { $R->{$aa}{phiStep} = 360.0; }

# go over all amino acids and rewrite library into new (binary) form with explicit coordinates
my $ofh = GENERAL::GetOutFH($outRotLibFile);
for (my $i = 0; $i < 19; $i++) {
  my $aa = SEQUENCE::index2t($i);
  printf("$aa\n");

  # -- start CHARMM
  my $chrm = CHARMM::new($CHARMM_DEF);#, "tmp.inp", undef, 1);
  $chrm->{par}->{param} = 19; # choose whether to use united-hydrogen or explicit-hydrogen model
  $chrm->loadParm();

  # -- create peptide
  # place CA at origin, C no the X-axis, and N in the XY plane
  # IMPORTANT: it is key that the seed atoms be N, CA, C (for proper placement of the side-chain in
  # absolute 3D space; see function readRotamerCoordinates()). This happens to be the default, but
  # other possibilities can also be allowed based on the IC table, so coding it in explicitly.
  $chrm->setupPeptide("A", GENERAL::arrayRef($aa), undef, undef, undef, "none", "none", ("seed" => "1 n 1 ca 1 c", "ori" => 0));

  # -- dump coordinates of all rotamers for this amino acid
  my $k = 0;
  my $scCA = "select (.NOT. (" . CHARMM::backbone() . ")) .OR. type CA end";
  foreach my $b (@{$R->{$aa}{bins}}) {
    my $phi = $b->{phi};
    my $psi = $b->{psi};
    foreach my $rot (@{$R->{$aa}{bin}{$phi}{$psi}}) {
      placeRotamer($chrm, "A", 1, $chidef->{$aa}, $rot);
      $chrm->writePDB(tmpPdbFile($tmpBase, $aa, $k), $scCA);
      $k++;
    }
  }
  $chrm->finish();

  # -- populate binary rotamer library with this amino acid
  # start with the header info for this amino acid
  my ($first, $firstAtoms) = readRotamerCoordinates(tmpPdbFile($tmpBase, $aa, 0));
  $ofh->print(pack("Z", $aa));
  $ofh->print(pack("i", scalar(@{$chidef->{$aa}})));
  $ofh->print(pack("i", scalar(@$firstAtoms)));
  $ofh->print(pack("i", scalar(@{$R->{$aa}{bins}})));
  foreach my $chi (@{$chidef->{$aa}}) {
    $ofh->print(pack("Z", @$chi));
  }
  foreach my $a (@$firstAtoms) {
    $ofh->print(pack("Z", $a->{atomname}));
  }

  # -- write individual rotamers
  $k = 0;
  foreach my $b (@{$R->{$aa}{bins}}) {
    my $phi = $b->{phi};
    my $psi = $b->{psi};
    $ofh->print(pack("f", ($phi, $psi)));
    $ofh->print(pack("f", ($R->{$aa}{phiStep}, $R->{$aa}{psiStep})));
    $ofh->print(pack("i", scalar(@{$R->{$aa}{bin}{$phi}{$psi}})));
    foreach my $rot (@{$R->{$aa}{bin}{$phi}{$psi}}) {
      foreach my $f ("chi", "chiSig") {
        for (my $xi = 0; $xi < scalar(@{$rot->{chi}}); $xi++) {
          $ofh->print(pack("f", $rot->{$f}->[$xi]));
        }
      }
      my ($rotPDB, $atoms) = readRotamerCoordinates(tmpPdbFile($tmpBase, $aa, $k));
      GENERAL::assert(scalar(@$firstAtoms) == scalar(@$atoms), "the 0-th and $k-th rotamers of amino acid $aa have different numbers of atoms!");
      for (my $ai = 0; $ai < scalar(@$atoms); $ai++) {
        GENERAL::assert($atoms->[$ai]->{atomname} eq $firstAtoms->[$ai]->{atomname}, "the 0-th and $k-th rotamers of $aa were dumped with atoms in different orders!");
        $ofh->print(pack("f", ($atoms->[$ai]->{xcoor}, $atoms->[$ai]->{ycoor}, $atoms->[$ai]->{zcoor})));
      }
      GENERAL::crm(tmpPdbFile($tmpBase, $aa, $k));
      $k++;
    }
  }
}
close($ofh);


sub placeRotamer {
  my $self = shift;
  my $chain = shift;
  my $iresnum = shift;
  my $chidef = shift;
  my $rot = shift;

  my $fhandle = $self->{files}->{inp};
  $fhandle->print("coor init select .byres. atom $chain $iresnum * .and. .not. (" . CHARMM::backbone() . ") end \n");
  $fhandle->print("ic param \n");
  $fhandle->print("ic edit \n\n");
  my $site = "$chain $iresnum";
  my $angles = $rot->{chi};
  # specify the dihedral angles to be rebuilt
  my $k = 0;
  foreach my $chi (@$chidef) {
    GENERAL::error("dihe $site $chi->[0] $site $chi->[1] $site $chi->[2] $site $chi->[3]. Chi angle number " . ($k+1) .
                   " was not found.") if (!defined $angles->[$k]);
    $fhandle->print("dihe $site $chi->[0] $site $chi->[1] $site $chi->[2] $site $chi->[3] $angles->[$k] \n");
    $k++;
  }
  $fhandle->print("end\n");
  $fhandle->print("ic build\n\n");
}

sub readDunbrackFormattedRotamerLibrary {
  my $rlibFile = shift;
  my $ifh = GENERAL::GetInFH($rlibFile);
  my %R;

  # Parse rotamer library. Example line:
  #ARG  -180 -180    10     1  2  2  1  0.249730    62.5   176.9   176.6    85.7       6.9    11.1    10.5     9.9
  while (1) {
    # first read the header for the amino acid
    my ($bwLine1, $bwLine2, $aaLine);
    last if (!GENERAL::skipTo($ifh, "step, deg", \$bwLine1));
    GENERAL::assert(GENERAL::skipTo($ifh, "step, deg", \$bwLine2), "found the first, but not the second phi/psi interval line!");
    GENERAL::assert(GENERAL::skipTo($ifh, "Residue type", \$aaLine), "could not find the amino-acid header line");
    GENERAL::assert(($aaLine =~ /Residue type\s+(\S+)\s*$/) ? 1 : 0, "could not parse out amino acid name from '$aaLine'");
    my $aa = $1;
    foreach my $bwLine ($bwLine1, $bwLine2) {
      GENERAL::assert(($bwLine =~ /(phi|psi) step, deg\s+(\S+)\s*$/) ? 1 : 0, "could not parse out bin width from '$bwLine'");
      $R{$aa}{$1."Step"} = $2;
    }
    GENERAL::assert(defined($R{$aa}{phiStep}) && defined($R{$aa}{psiStep}), "either phi or psi step size was not found for amino acid $aa");
    $R{$aa}{bins} = GENERAL::arrayRef();
    $R{$aa}{bin} = {};
    $R{$aa}{N} = 0;

    # then read the rotamers
    my $rotLines = "";
    GENERAL::skipTo($ifh, "Backbone-dependent rotamer", undef, \$rotLines);
    my @rotLines = split("\n", $rotLines);
    foreach my $line (@rotLines) {
      $line =~ s/\#.*$//g;
      $line = GENERAL::Trim($line);
      next if ($line eq "");
      my @line = split(" ", $line);
      GENERAL::assert(scalar(@line) == 17, "unexpected number of entries in rotamer library line '$line'");
      GENERAL::assert($aa eq $line[0], "was expecting rotamers for $aa when found line '$line', containing a rotamer for $line[0]");
      my $phi = $line[1];
      my $psi = $line[2];
      my $pr = $line[8];
      my (@chi, @chiSig);
      for (my $i = 0; $i < 4; $i++) {
        next if ($line[$i+4] == 0); # if chi index is marked as zero, that means this residues does not have this chi angle
        push(@chi, $line[$i+9]);
        push(@chiSig, $line[$i+13]);
      }
      my %r = ("chi" => \@chi, "chiSig" => \@chiSig, "p" => $pr);
      if (!defined($R{$aa}{bin}{$phi}) || !defined($R{$aa}{bin}{$phi}{$psi})) { my %b = ("phi" => $phi, "psi" => $psi); push(@{$R{$aa}{bins}}, \%b); }
      push(@{$R{$aa}{bin}{$phi}{$psi}}, \%r);
      $R{$aa}{N}++;
    }
  }
  close($ifh);

  return \%R;
}

sub readChiDefinition {
  my $chiDefFile = shift;
  my $ifh = GENERAL::GetInFH($chiDefFile);
  my %C;
  my $aa = undef;

  foreach my $line (<$ifh>) {
    if ($line =~ /^\s*RESI\s+(\S+)\s*$/) {
      $aa = $1;
      $C{$aa} = GENERAL::arrayRef();
      next;
    }
    $line = GENERAL::Trim($line);
    next if ($line eq "");
    my @chi = split(" ", $line);
    push(@{$C{$aa}}, \@chi);
  }  
  close($ifh);

  return \%C;
}

sub tmpPdbFile {
  my $base = shift;
  my $aa = shift;
  my $k = shift;

  return lc($base . "-$aa-$k.pdb");
}

sub readRotamerCoordinates {
  my $pdbf = shift;
  my $pdb = PDB::new($pdbf);
  my $res = $pdb->{chain}->[0]->{res}->[0];
  my @atoms = $pdb->conAtoms();

  # align the residue such that CA is at the origin, N-CA is along the X axis, and N-CA-C is in the XY plane
  # IMPORTANT: because the seed atoms for generating the amino-acid coordinates in CHARMM were N, CA, C, that
  # means that N was placed at the origin, N-CA along the X-axis, and N-CA-C in the XY plane. So that means
  # all we have to do is move CA to the origin and everything will be as we want.
  my $CA = PDB::getAtomInRes($res, "CA");
  my $dx = $CA->{xcoor}; my $dy = $CA->{ycoor}; my $dz = $CA->{zcoor};
  foreach my $a (@atoms) {
    $a->{xcoor} -= $dx;
    $a->{ycoor} -= $dy;
    $a->{zcoor} -= $dz;
  }
  
  return ($pdb, \@atoms);
}
