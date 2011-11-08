import iceIh

#roo, cutt (cutoff for potential [A])
print "Calling cons."
iceIh.cons(2.75,10)
### Arg list:
#nx, ny, nz: number of layers in x, y, z dir.0
print "Calling build:"
iceIh.build(8,4,10)
print "Calling anal"
iceIh.anal()
### Arg list:
#dseed (dipole seed?)
#nimp (number of improvement rounds)
#nmcma (number of steps in round)
#maxt (max tries)
#isto (after how many futile attempts to stop?)
#cri1 (e-)
#cri2 (dip. identity criterions)
print "Calling hydro"
iceIh.hydro()
#trans writes graphic file for oxygen lattice
#if 1 = oxygens only
#if 2 include hydrogens.
print "Calling trans ..."
iceIh.trans(2)
