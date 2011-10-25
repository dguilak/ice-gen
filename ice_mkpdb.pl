#!/usr/bin/perl

# process file from crystal generator - 1 value per line
open(IN, "fort.9") ||
  die("cannot open file fort.9 -- $!\n");

# make one-column list
$fields = 0;
my @list;
while (<IN>) {
  @input = split();
  foreach(@input) {
    $list[$fields] = $_ * 0.5291778;
    $fields++;
  }
}
close(IN);

# make triplets of data
my @coords;
my $tempstr;
$line = 0;
$fields = 0;
for $c (0..$#list) {
  if ($c <= 5) { next; }	# skip 6 lines
  $tempstr = join(' ', $tempstr, $list[$c]);
  $fields++;
  if ($fields == 3) {
    $coords[$line] = $tempstr;
    $line++;
    $tempstr = "";
    $fields = 0;
  }
}

undef(@list);
undef($tempstr);

# write PDB file
$natom = 0;
$nres = 1;
my @pdbline;
my @tempstr;

print("REMARK PBC $coords[0]\n");
shift(@coords);

$line = 0;
for $c (0..$#coords) {
  $_ = $coords[$c];
  @tempstr = split();
  $natom++;

  $pdbline[0] = "ATOM";
  $pdbline[1] = $natom;
  $pdbline[3] = "SPC";
  $pdbline[4] = $nres;
  $pdbline[5] = $tempstr[0];
  $pdbline[6] = $tempstr[1];
  $pdbline[7] = $tempstr[2];
  
  if ((($c + 1) % 3) == 1) {
    $pdbline[2] = "OW";
  } elsif ((($c + 1) % 3) == 2) {
    $pdbline[2] = "HW1";
  } else {
    $pdbline[2] = "HW2";
  }
  printf("%-6s%5d%4s  %3s  %4d    %8.3f%8.3f%8.3f\n", @pdbline);
  if ((($c + 1) % 3) == 0) {
    print("TER\n");
    $nres++;
  }
}

print("END\n");
