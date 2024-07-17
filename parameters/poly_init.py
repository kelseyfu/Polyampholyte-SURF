#!/usr/bin/python
 
# Script:  polyanion_init.py
# Purpose: Make input file
# Example: python polyanion_init.py
# Author:  Pierre Walker

# lammps types    
# 1: Cation bead on polymer
# 2: Anion bead on polymer
# 3: Neutral bead on polymer
# 4: Cation salt
# 5: Anion salt

import sys
from numpy import *
from random import *

#-------------------------------------------------------------------------
iseed = 9328                 # random number seed
npoly = NCHAIN                   # number of polymers
bond = 0.97                  # bond length. depends on bond potential, but close to 1 is good enough
minsep = 1.0                 # allowed separation in overlap check

cisize = 1.0                 # Ion size
z_c =  1                     # counterion valence                 
Lx   = LENGTHX               # box size in x
Ly   = LENGTHY               # box size in y
Lz   = LENGTHZ               # box size in z

minsep2 = minsep*minsep

sequence = SEQUENCE
nmonomers = len(sequence)
print("N monomers: "+str(nmonomers))

net_charge = sum(sequence)

ncounterions = int(abs(net_charge*npoly))
# END INPUT Parameters ------------------------------------------------------

INPUT_LAMMPS = open('input.data', 'w')

ntypes = 5 # number of atom types

# Initialize ncounterions to 0 for now, will calculate later

ntot = nmonomers*npoly + ncounterions # total number of particles
vol  = Lx*Ly*Lz # volume of the box

dens = ntot/vol # density

dim = int(ntot+1)

# simulation cell parameters
hx = Lx
hy = Ly
hz = Lz

hx2 = hx/2.
hy2 = hy/2.
hz2 = hz/2.

nbonds = npoly*nmonomers-npoly # number of bonds

print("\n")
print("Total number of particles: "+str(ntot)+"\n")
print("Density = "+str(dens)+"\n")
print("Number of chains = "+ str(npoly)+"\n")


print("Number of atoms types = "+str(ntypes)+"\n")
print("seed = "+str( iseed)+"\n")

# init position variables
xc=zeros(dim,float32)
yc=zeros(dim,float32)
zc=zeros(dim,float32)
cx=zeros(dim)
cy=zeros(dim)
cz=zeros(dim)
typeb=[0]*dim
molnum=[0]*dim
q=[0.0]*dim

# Build polymers
k=0
                            
for ix in range(npoly):
    lengthcurrentpoly = 0
    for iy in range(nmonomers):
        q[k] = sequence[iy]
        if sequence[iy] == 1.:
            typeb[k] = 1 # cation on polymer
        elif sequence[iy] == -1.:
            typeb[k] = 2 # anion on polymer
        else:
            typeb[k] = 3 # neutral on polymer

        molnum[k] = ix + 1 # Assign molecule number
        if iy == 0: # First monomer in polymer
            xc[k] = (random()-0.5)*hx # Random x position
            yc[k] = (random()-0.5)*hy # Random y position
            zc[k] = (random()-0.5)*hz # Random z position
        else:
            dx = random()-0.5 # Random x displacement
            dy = random()-0.5 # Random y displacement
            dz = random()-0.5 # Random z displacement
            r = sqrt(dx*dx+dy*dy+dz*dz) # Distance from previous monomer
            scale = bond/r # Scale factor to bond length
            dx = scale*dx # Scale x displacement
            dy = scale*dy # Scale y displacement
            dz = scale*dz # Scale z displacement 
            xc[k] = xc[k-1] + dx # New x position
            yc[k] = yc[k-1] + dy # New y position
            zc[k] = zc[k-1] + dz # New z position 
        k = k + 1   
        lengthcurrentpoly = lengthcurrentpoly + 1 # Increment polymer length
print(typeb[199:201])

# Adjust for periodic boundary conditions when bead is outside the box
for k in range(npoly*nmonomers):
    if (xc[k] > hx):
        cx[k] = int(xc[k]/hx)
        xc[k] = xc[k] - cx[k]*hx - hx2
    elif (xc[k] < 0.0):
        cx[k] = -int((-xc[k]+hx)/hx)
        xc[k] = xc[k] - cx[k]*hx - hx2
    else:
        cx[k] = 0
        xc[k] = xc[k] - hx2
    if (yc[k] > hy):
        cy[k] = int(yc[k]/hy)
        yc[k] = yc[k] - cy[k]*hy - hy2
    elif (yc[k] < 0.0):
        cy[k] = -int((-yc[k]+hy)/hy)
        yc[k] = yc[k] - cy[k]*hy - hy2
    else:
        cy[k] = 0
        yc[k] = yc[k] - hy2
    if (zc[k] > hz):
        cz[k] = int(zc[k]/hz)
        zc[k] = zc[k] - cz[k]*hz - hz2
    elif (zc[k] < 0.0):
        cz[k] = -int((-zc[k]+hz)/hz)
        zc[k] = zc[k] - cz[k]*hz - hz2
    else:
        cz[k] = 0
        zc[k] = zc[k] - hz2


print("Polymers built."+"\n")

# Add counterions to neutralise system
for ii in range(1,abs(ncounterions)+1):
    k = ii + ntot - ncounterions - 1
    xc[k] = (random()-0.5)*hx 
    yc[k] = (random()-0.5)*hy
    zc[k] = (random()-0.5)*hz
    typeb[k] = 4
    q[k] = sign(ncounterions)*-1

    

print("Counterions complete."+"\n")

# OUTPUT headers ---------------------------------------------------------------


# input.psf header line
##INPUT_PSF.write("PSF\n\n%8i !NATOM\n" % (ntot))
##INPUT_PDB.write("CRYST1  %7.3f  %7.3f  %7.3f %6.2f %6.2f %6.2f P 1           1\n" % (hx,hy,hz,90.0,90.0,90.0)) #Mark says this will allow use of pbctools in VMD

# input.lammps header 
INPUT_LAMMPS.write("# Polyanions PJW 8/2022\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" %     ntot)
INPUT_LAMMPS.write("%10i    bonds\n" %     nbonds)
INPUT_LAMMPS.write("%10i    angles\n" %     0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
##INPUT_LAMMPS.write("%10i    impropers\n" % npendant)
INPUT_LAMMPS.write("%10i    impropers\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % 1)
##INPUT_LAMMPS.write("%10i    improper types\n" % 1)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2,hx2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2,hy2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2,hz2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("1 1\n")
INPUT_LAMMPS.write("2 1\n")
INPUT_LAMMPS.write("3 1\n")
INPUT_LAMMPS.write("4 1\n")
INPUT_LAMMPS.write("5 1\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Pair Coeffs # lj/cut/coul/long\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("1 1 1\n")
INPUT_LAMMPS.write("2 1 1\n")
INPUT_LAMMPS.write("3 1 1\n")
INPUT_LAMMPS.write("4 1 1\n")
INPUT_LAMMPS.write("5 1 1\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bond Coeffs # fene\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("1 30 1.5 1 1\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

# END OUTPUT headers -----------------------------------------------------------

# Atoms output
mass = 1.0

# Polymers
i = 0
imol = 0

# positions 
for i in range(dim-1):
    itype = typeb[i]

    # could use a dictionary here between type and segname
    if itype != 4 and itype != 5:
        imol = molnum[i] #this implies the polymers must be placed in before the counterions
    else:
        #imol = npoly+1 #the molecule number for all counterions is the same; it's more like a group number
        imol = i-ntot+ncounterions+npoly+1 #LMH now each ion has its own molecule number
    INPUT_LAMMPS.write("%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i+1, imol, itype, q[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))
# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
#jbond1 = zeros(nbonds+1)
#jbond2 = zeros(nbonds+1)
pbond1 = zeros(nbonds+1)
pbond2 = zeros(nbonds+1)
ibond=0
pbond=0
i0 = 0
for i in range(ntot-abs(ncounterions)):
        #if not at the end of the polymer
        if molnum[i+1] == molnum[i]:
            ibond = ibond+1 #the bond number
            j=i+1
            INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond,i+1,j+1))

    
# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")

for ii in range(1,ntypes+1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)

#Close files
INPUT_LAMMPS.close()
##INPUT_PDB.close()
##INPUT_PSF.close()
##
print("LAMMPS output complete."+"\n")
