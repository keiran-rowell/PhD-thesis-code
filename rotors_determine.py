#Written by Keiran Rowell - PhD Candidate UNSW - May-June 2018
#Inspiration of problem set up from https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/

import argparse
import subprocess
import numpy
import scipy.constants
from scipy.linalg import eigh
import math

#Isotopic masses purloined from https://www.lfd.uci.edu/~gohlke/code/elements.py.html, which got them from NIST
#Masses above beyond period 3 have NOT been rigorously double-checked for entry errors
atomic_mass = dict(H=1.0078250321, He=4.0026032497, Li=7.016004, Be=9.0121821, B=11.0093055, C=12.0,
                   N=14.0030740052, O=15.9949146221, F=18.9984032, Ne=19.9924401759, Na=22.98976967, Mg=23.9850419,
                   Al=26.98153844, Si=27.9769265327, P=30.97376151, S=31.97207069, Cl=34.96885271, Ar=39.962383123,
                   K=38.9637069, Ca=39.9625912, Sc=44.9559102, Ti=47.9479471, V=50.9439637, Cr=51.9405119,
                   Mn=54.9380496, Fe=55.9349421, Co=58.9332002, Ni=57.9353479, Cu=62.9296011, Zn=63.9291466,
                   Ga=68.925581, Ge=73.9211782, As=74.9215964, Se=79.9165218, Br=78.9183376, Kr=83.911507,
                   Rb=84.9117893, Sr=87.9056143, Y=88.9058479, Zr=89.9047037, Nb=92.9063775, Mo=97.9054078,
                   Tc=97.907216, Ru=101.9043495, Rh=102.905504, Pd=105.903483, Ag=106.905093,
                   Cd=113.9033581, In=114.903878, Sn=119.9021966, Sb=120.903818, Te=129.9062228,
                   I=126.904468, Xe=131.9041545, Cs=132.905447, Ba=137.905241, La=138.906348,
                   Ce=139.905434, Pr=140.907648, Nd=141.907719, Pm=144.912744, Sm=151.919728,
                   Eu=152.921226, Gd=157.924101, Tb=158.925343, Dy=163.929171, Ho=164.930319,
                   Er=165.93029, Tm=168.934211, Yb=173.9388581, Lu=174.9407679, Hf=179.9465488,
                   Ta=180.947996, W=183.9509326, Re=186.9557508, Os=191.961479, Ir=192.962924,
                   Pt=194.964774, Au=196.966552, Hg=201.970626, Tl=204.974412, Pb=207.976636,
                   Bi=208.980383, Po=208.982416, At=209.987131, Rn=222.0175705, Fr=223.0197307,
                   Ra=226.0254026, Ac=227.027747, Th=232.0380504, Pa=231.0358789, U=238.0507826,
                   Np=237.0481673, Pu=244.064198, Am=243.0613727, Cm=247.070347, Bk=247.070299,
                   Cf=251.07958, Es=252.08297, Fm=257.095099, Md=258.098425, No=259.10102,
                   Lr=262.10969, Rf=261.10875, Db=262.11415, Sg=266.12193, Bh=264.12473,
                   Hs=269.13411, Mt=268.13882)

class Atom(object):
    def __init__(self, index, symbol, x_abs, y_abs, z_abs):
        self.index = index 
        self.symbol = symbol
        self.mass = atomic_mass[symbol]
        self.x_abs = x_abs
        self.y_abs = y_abs
        self.z_abs = z_abs
        self.x_rel = None
        self.y_rel = None
        self.z_rel = None

class COM(object):
    def __init__(self, x_abs, y_abs, z_abs):
        self.x_abs = x_abs
        self.y_abs = y_abs
        self.z_abs = z_abs

class Moms_inertia(object):
    def __init__(self, Ix, Iy, Iz, Ia, Ib, Ic):
        self.Ix = Ix
        self.Iy = Iy
        self.Iz = Iz
        self.Ia = Ia 
        self.Ib = Ib
        self.Ic = Ic
        self.I2d = None
        self.Ik = None

class Rot_consts(object):
    def __init__(self, Bx, By, Bz, Ba, Bb, Bc):
        self.Bx = Bx
        self.By = By
        self.Bz = Bz
        self.Ba = Ba
        self.Bb = Bb
        self.Bc = Bc
        self.B2d = None
        self.Bk = None
    
def parse_args():
    parser = argparse.ArgumentParser(description='Calculate moments of inertia, rotational constants, assign top type from .xyz file.')
    parser.add_argument('-f',"--filename", help=".xyz geometry files, with the 1st column the atomic symbol", type=str, required=True)
    parser.add_argument('-u',"--units", help="The units you want the output in of: [cm-1, GHz]", type=str, \
    choices=['cm-1', 'GHz'], default='GHz') 
    parser.add_argument('-I',"--moment_inertia", help="The moment(s) of inertia to be printed", type=str, \
    nargs = '*', choices=['Ix','Iy','Iz','Ia','Ib','Ic','I2d','Ik'])
    parser.add_argument('-B',"--rotational_const", help="The rotational constant(s) to be printed", type=str, \
    nargs = '*', choices=['Bx','By','Bz','Ba','Bb','Bc','B2d','Bk']) 
    args = parser.parse_args()
    return args

#Recognising equivalent moments of inertia after the trials of floating point calculations
inertia_tol = 1e-3
#Conversion factors
wavenum_to_GHz = 29.9792458    

def main():
    args = parse_args()
    fname = args.filename       
    with open(fname) as f: 
        file_contents = f.read().splitlines()[2:] #skip #atoms, comment line
    atom_list = assign_atoms(file_contents)
    xyz_mass = create_xyz_mass_matrix(atom_list)
    com = find_com(xyz_mass)
    atom_list = determine_relative_dists(atom_list,com)
    inertia_tensor = form_inertia_tensor(atom_list,com)
    moms_inertia = get_moments_of_inertia(inertia_tensor)
    top_type = assign_top_type(moms_inertia)
    moms_inertia = get_J_K_rotors(moms_inertia, top_type)
    rot_consts = get_rot_consts(moms_inertia, top_type)
    print_output(moms_inertia, rot_consts, args)

def assign_atoms(file_contents):
    atom_list = []

    for index, line in enumerate(file_contents):
        symbol = line.split()[0]
        x_abs = float(line.split()[1])
        y_abs = float(line.split()[2])
        z_abs = float(line.split()[3])
        atom = Atom(index, symbol, x_abs, y_abs, z_abs)
        atom_list.append(atom)

    return atom_list

def create_xyz_mass_matrix(atom_list):
    xyz_mass = numpy.zeros(shape=(0,4))

    for atom in atom_list:
        xyzm_vals = atom.x_abs, atom.y_abs, atom.z_abs, atom.mass
        xyz_mass = numpy.append(xyz_mass, [xyzm_vals] ,axis=0)

    return xyz_mass

def find_com(xyz_mass):
    CM = numpy.average(xyz_mass[:,:3], axis=0, weights=xyz_mass[:,3])
    com = COM(CM[0],CM[1],CM[2])

    return com

def determine_relative_dists(atom_list,com):
    for atom in atom_list: 
        atom.x_rel = atom.x_abs - com.x_abs
        atom.y_rel = atom.y_abs - com.y_abs
        atom.z_rel = atom.z_abs - com.z_abs

    return atom_list

def form_inertia_tensor(atom_list,com):
    inertia_tensor = numpy.zeros(shape=(3,3)) 
    
    Ixx, Ixy, Ixz, Iyz, Iyy, Iyz, Izx, Izy, Izz = [0] * 9

    for atom in atom_list:
        Ixx = Ixx + atom.mass * (atom.y_rel ** 2 + atom.z_rel ** 2)
        Ixy = Ixy - atom.mass * (atom.x_rel * atom.y_rel)  
        Ixz = Ixz - atom.mass * (atom.x_rel * atom.z_rel) 
        Iyy = Iyy + atom.mass * (atom.x_rel ** 2 + atom.z_rel ** 2) 
        Iyz = Iyz - atom.mass * (atom.y_rel * atom.z_rel) 
        Izz = Izz + atom.mass * (atom.x_rel ** 2 + atom.y_rel ** 2) 

    inertia_tensor[0,0], inertia_tensor[0,1], inertia_tensor[0,2] = Ixx, Ixy, Ixz
    inertia_tensor[1,0], inertia_tensor[1,1], inertia_tensor[1,2] = Ixy, Iyy, Iyz
    inertia_tensor[2,0], inertia_tensor[2,1], inertia_tensor[2,2] = Ixz, Iyz, Izz
        
    return inertia_tensor

def get_moments_of_inertia(inertia_tensor):
    w, v =  eigh(inertia_tensor)
    Ix, Iy, Iz = w[2], w[1], w[0]
    Ia, Ib, Ic = sorted(w)[0], sorted(w)[1], sorted(w)[2]

    moms_inertia = Moms_inertia(Ix, Iy, Iz, Ia, Ib, Ic)

    return moms_inertia

def get_rot_consts(moms_inertia, top_type):
    if moms_inertia.Ix <= inertia_tol: #I should have this iterating over Ix, Iy, Iz not 3 cases
        Bx = None
    else:
        Bx = rotconst_to_mominertia(moms_inertia.Ix)

    if moms_inertia.Iy <= inertia_tol:
        By = None
    else:
        By = rotconst_to_mominertia(moms_inertia.Iy)

    if moms_inertia.Iz <= inertia_tol:
        Bz = None
    else:
        Bz = rotconst_to_mominertia(moms_inertia.Iz)
    
    if top_type == "linear":
        Ba = None
    else:
        Ba = rotconst_to_mominertia(moms_inertia.Ia)
    Bb = rotconst_to_mominertia(moms_inertia.Ib)
    Bc = rotconst_to_mominertia(moms_inertia.Ic)
    B2d = rotconst_to_mominertia(moms_inertia.I2d)
    Bk = rotconst_to_mominertia(moms_inertia.Ik)
   
    rot_consts = Rot_consts(Bx, By, Bz, Ba, Bb, Bc,) 
    
    rot_consts.B2d = B2d
    rot_consts.Bk = Bk

    return rot_consts
  
def rotconst_to_mominertia(I):
    if I != 'None':
       B = ( (scipy.constants.h * 1E7) / (8*(scipy.constants.pi ** 2) * (scipy.constants.c * 1E2) * (I * 1.6605654723603E-40)) ) #cgs... 
    return B 

def assign_top_type(moms_inertia):
    if is_totally_symmetric(moms_inertia):
        top_type = "totally symmetric"
    elif is_linear(moms_inertia):
        top_type = "linear"
    elif is_oblate(moms_inertia):
        top_type = "oblate"
    elif is_prolate(moms_inertia):
        top_type = "prolate"
    elif is_asymmetric(moms_inertia):
        top_type = "asymmetric"

    return top_type 

#Have to check tolerance for comparisons due to floating point errors
#There's redundancy but I can write (and alter) each top criteria seperately
def is_totally_symmetric(moms_inertia):
    if abs(moms_inertia.Ia - moms_inertia.Ib) <= inertia_tol and  abs(moms_inertia.Ib - moms_inertia.Ic) <= inertia_tol \
    and abs(moms_inertia.Ia - moms_inertia.Ic) <= inertia_tol:
            return True
    else:
        return False

def is_oblate(moms_inertia):
    if abs(moms_inertia.Ia - moms_inertia.Ib) <= inertia_tol:
        if abs(moms_inertia.Ib - moms_inertia.Ic) > inertia_tol and moms_inertia.Ib < moms_inertia.Ic:
            return True
    else:
        return False
        
def is_prolate(moms_inertia):
    if abs(moms_inertia.Ia - moms_inertia.Ib) > inertia_tol and moms_inertia.Ia < moms_inertia.Ib:
        if abs(moms_inertia.Ib - moms_inertia.Ic) <= inertia_tol:
            return True
    else:
        return False
       
def is_asymmetric(moms_inertia):
    if abs(moms_inertia.Ia - moms_inertia.Ib) > inertia_tol and abs(moms_inertia.Ib - moms_inertia.Ic) > inertia_tol \
    and abs(moms_inertia.Ia - moms_inertia.Ic) > inertia_tol:
        return True
    else:
        return False
def is_linear(moms_inertia):
    if moms_inertia.Ia <= inertia_tol: #linear molecule has no moment of inertia along a-axis
        return True
    else:
        return False


def get_J_K_rotors(moms_inertia, top_type):
    if top_type == "prolate" or top_type == "oblate" or top_type == "asymmetric":       
    #Approximation formulas from the MultiWell Feb 13 2017 manual pg 141
        moms_inertia.I2d = math.sqrt(moms_inertia.Ib * moms_inertia.Ic)
        moms_inertia.Ik = ((moms_inertia.Ia ** -1) - (moms_inertia.I2d ** -1)) ** -1
    elif top_type == "totally symmetric": #somewhat redundant to compute b.c. of defintion       
        moms_inertia.I2d = math.sqrt(moms_inertia.Ib * moms_inertia.Ic)
        moms_inertia.Ik = moms_inertia.I2d #spherical so J = K etc
    
    return moms_inertia

def print_output(moms_inertia, rot_consts, args): #TODO: should check atrribute I'm accessing not 'None' (e.g. linear)
    if args.moment_inertia != None:
        for inertia_label in args.moment_inertia:
            if args.units == 'cm-1': #this got too messy when I put the condition in .format()
                print("{}: {}").format(inertia_label, getattr(moms_inertia, inertia_label)) #use getattr for attribute from string
            elif args.units == 'GHz': 
                print("{}: {}").format(inertia_label, getattr(moms_inertia, inertia_label)*wavenum_to_GHz) 
    if args.rotational_const != None:
        for rotational_label in args.rotational_const:
            if args.units == 'cm-1':
                print("{}: {}").format(rotational_label, getattr(rot_consts, rotational_label)) 
            elif args.units == 'GHz': 
                print("{}: {}").format(rotational_label, getattr(rot_consts, rotational_label)*wavenum_to_GHz) 


if __name__ == '__main__':
        main()
