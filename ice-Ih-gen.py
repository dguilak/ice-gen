import iceIh

#roo, cutt (cutoff for potential [A])
iceIh.cons(2.75,10)
### Arg list:
#nx, ny, nz: number of layers in x, y, z dir.0
iceIh.build(8,4,10)
iceIh.anal()
### Arg list:
#dseed (dipole seed?)
#nimp (number of improvement rounds)
#nmcma (number of steps in round)
#maxt (max tries)
#isto (after how many futile attempts to stop?)
#cri1 (e-)
#cri2 (dip. identity criterions)
iceIh.hydro(float(1), 1000, 10000, 400, 400, float(0.0001), float(0.0001))
#trans writes graphic file for oxygen lattice
#if 1 = oxygens only
#if 2 include hydrogens.
#write(99,777) *not* working here.
#writ needs to get called as well!
iceIh.trans(2)
