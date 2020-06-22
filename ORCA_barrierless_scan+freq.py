import argparse
import subprocess

#Creates a series of xyz files which should got into their own folders to be submitted to a frequency calculation

def main():
    args = parse_args()
    write_scan_input_file(args)   #I'd like to have both scan and freq use same funciton with keyword arguements
    write_scan_pbs_script(args)
    write_freqs_pbs_script(args)
    write_processing_pbs_script(args)
    if args.qsub == 'y':
        subprocess.call(['qsub',args.filename+'.pbs']) #need to check this works and how to see when it terminates
    #now I need to parse the freq output into .vib format for MultiWell
    #Complete assignment still can't be automated? I think needs human decisions, still need the .vib files at hand

def parse_args(): #I need to clean up some of these unused options
    parser = argparse.ArgumentParser(description='Parse input options to automatically generate ORCA relaxed scan file.')

    parser.add_argument('-f',"--filename", help="The name of your scan input file", type=str, required=True)
    parser.add_argument('-com',"--comments", help="Comments in input file", default='off', choices=['on','off'], type=str)
    parser.add_argument('-nm',"--name", help="The Orca Job name", default='default name', type=str)
    parser.add_argument('-n',"--ncpus", help="The number of CPUS to run ORCA job with", default=4, type=int)
    parser.add_argument('-rst',"--restricted", help="Restricted (R) or unrestricted (U)", default='U', choices=['R','U'], type=str)
    parser.add_argument('-mtd',"--method", help="The ab initio method you wish to use", default='B2GP-PLYP', type=str)
    parser.add_argument('-bs',"--basis_set", help="The job's basis set (check ORCA basis-set library)", default='def2-TZVP', type=str)
    parser.add_argument('-d',"--dispersion", help="Dispersion correction", default='D3BJ', choices=['D2', 'D3ZERO', 'D3','D3BJ','none'], type=str) 
    parser.add_argument('-aprx',"--approx", help="The approximation speed-ups (e.g. RIJK, RIJCOSX etc)", default='RIJK', \
                       choices=['RIJK', 'RIJCOSX', 'RIJNOX', 'RI-J','none'], type=str)
    parser.add_argument('-ri',"--ri_aux", help="If an auxillary basis set for RI is needed", default='off', choices=['on','off'], type=str)
    parser.add_argument('-scf',"--scf_tol", help="The tolerance to SCF convergence (loose, tight) etc", default='tight', 
                        choices=['sloppy','loose','normal','strong','tight','verytight','extreme'], type=str)
    parser.add_argument('-g',"--grid", help="The grid size (Grid4, Grid6 etc.) used for calculation", default=5, choices=[1,2,3,4,5,6,7], type=int)
    parser.add_argument('-j',"--job_type", help="The type of ORCA job to run (usually opt)", default='opt', type=str)
    parser.add_argument('-a',"--atom_a", help="The 1st atom of the scanned bond coordinate", type=int, required=True)
    parser.add_argument('-b',"--atom_b", help="The 2nd atom of the scanned bond coordinate", type=int, required=True)
    parser.add_argument('-r_s',"--r_start", help="Interatomic distance and the start of scan", type=float, required=True)  
    parser.add_argument('-r_e',"--r_end", help="Interatomic distance to terminate scan at", type=float, required=True)
    parser.add_argument('-i',"--increments", help="How many steps in between start and end bond distance", type=int, required=True )
    parser.add_argument('-c',"--charge", help="Total molecular charge", default=0, type=int)
    parser.add_argument('-s',"--multiplicity", help="Total molecular multiplicity", type=int, required=True)
    parser.add_argument('-xyz',"--xyzfile", help="xyzfile for the molecule's coordinates at start of the scan", type=str, required=True)
    parser.add_argument('-clst',"--cluster", help="The computer cluster name to run on (changes PBS queue submission script)", default='katana', \
                        choices=['katana','raijin','caynow'], type=str)
    parser.add_argument('-mem',"--memory", help="The amount of memory for ORCA (in GB)", default=10, type=int)
    parser.add_argument('-t',"--time", help="The amount of walltime for ORCA (in hrs)", default=12,  type=int)
    parser.add_argument('-q',"--qsub", help="Automatically qsub generated files to the PBS queue (y/n)", default='n', choices=['y','n'])
    
    args = parser.parse_args()
    return args

def write_scan_input_file(args):
    with open(args.filename+'.inp','w') as f:  
        f.write('# '+args.name+" \n")
        f.write('! '),
        if args.ncpus != 1:
                f.write('PAL' + str(args.ncpus) + ' ') #often good to double int in freq jobs (e.g. PAL8), however easy to edit
        if args.restricted == 'U':
                f.write('UKS '),
        elif args.restricted == 'R':
                f.write('RKS '),
        if args.approx != 'none':
                f.write('RI-'+args.method+' ')
        else:
                f.write(args.method+' ')
        f.write(args.basis_set+' '),
        if args.dispersion != 'none':
                f.write(args.dispersion+' ')
        if args.comments != 'off':
            f.write('# Number of CPUs and level of theory'),
        f.write("\n")
        f.write('! '),
        if args.approx != 'none':
                f.write(args.approx+' '),
        if args.approx == 'RIJK':
                f.write('def2/JK '),
        if args.method in ['MP2','B2GP-PLYP', 'B2-PLYP', 'B2K-PLYP', 'B2T-PLYP'] and args.approx != 'none':
            args.ri_aux = 'on'
        if args.ri_aux == 'on':
                f.write(args.basis_set+'/C '),
        if args.comments != 'off':
                f.write('# Density fitting/resoultion-of-identity speed-ups'),
        f.write("\n")
        f.write('! '),
        f.write(args.scf_tol+'SCF '),
        f.write('Grid'+str(args.grid)+' '),
        if args.comments != 'off':
                f.write('# Convergence tolerence settings'),
        f.write("\n")
        f.write('! '),
        f.write(args.job_type+' '),
        f.write("\n")
        f.write('%geom scan')
        f.write("\n")
        f.write('B '+str(args.atom_a)+' '+str(args.atom_b)+' = '+str(args.r_start)+', '+str(args.r_end)+', '+str(args.increments))
        f.write("\n")
        f.write('end\n')
        f.write('end\n')
        f.write("\n")
        f.write('*xyzfile '+str(args.charge)+' '+str(args.multiplicity)+' '+args.xyzfile)
        f.write("\n")

def write_scan_pbs_script(args):
    with open(args.filename+'.pbs','w') as f:
        if args.cluster == "katana":
            f.write("#!/bin/bash\n")
            f.write("\n")
            f.write("#PBS -l nodes=1:ppn="+str(args.ncpus)+"\n")
            f.write("#PBS -l vmem="+str(args.memory)+"gb"+"\n")
            f.write("#PBS -l walltime="+str(args.time)+":00:00\n")
            f.write("#PBS -j oe\n") 
            f.write("\n")
            f.write("module add openmpi\n")
            f.write("module add orca/4.0.1.2\n")
            f.write("\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("\n")
            f.write("/share/apps/orca/4.0.1.2/bin/orca "+args.filename+".inp > "+args.filename+".out\n" )
            f.write("\n")
            point_radius = args.r_start
            increment_step = (args.r_end - args.r_start)/(args.increments-1)
            xyz_point_num = 1 
            for optismised_points in range(args.increments): #Ideally would be a bash loop
                xyz_point_val=str("{0:03}".format(xyz_point_num))
                dir='r_'+str("{0:.2f}".format(point_radius))
                subprocess.call(['mkdir',dir])  #Could check if dir exists to avoid complaints
                f.write("mv "+args.filename+"."+xyz_point_val+".xyz "+dir+"\n")
                f.write("mv "+args.filename+"."+xyz_point_val+".gbw "+dir+"\n")
                write_freq_input_file(args,dir,xyz_point_val)
                point_radius += increment_step
                xyz_point_num += 1
                f.write("\n")
            f.write("qsub "+args.filename+"_freqs.pbs\n")
            f.write("qsub "+args.filename+"_process.pbs\n")

def write_freqs_pbs_script(args):
    with open(args.filename+'_freqs.pbs','w') as f:
        if args.cluster == "katana":
            f.write("#!/bin/bash\n")
            f.write("\n")
            f.write("#PBS -l nodes=1:ppn="+str(args.ncpus)+"\n")
            f.write("#PBS -l vmem="+str(args.memory)+"gb"+"\n")
            f.write("#PBS -l walltime="+str(args.time)+":00:00\n")
            f.write("#PBS -j oe\n") 
            f.write("\n")
            f.write("#PBS -t 0-"+str(args.increments-1)+"\n") #Arrays start at 0
            f.write("\n")
            f.write("RADII=(") 
            point_radius = args.r_start #I should really pass XYZVALS, POINTS as an array
            increment_step = (args.r_end - args.r_start)/(args.increments-1)
            for point in range(args.increments):
                pnt_rad_formatted=str("{:.2f}".format(point_radius)) # trailing zeros for 2 d.p.
                f.write("'"+pnt_rad_formatted+"' ")
                point_radius += increment_step
            f.write(")\n") 
            f.write("\n")
            xyz_point_num = 1 #ORCA xyzfiles numbered from 1
            f.write("XYZ_VAL=(") 
            for point in range(args.increments):
                xyz_point_val=str("{0:03}".format(xyz_point_num)) #0 filled formatting
                f.write("'"+xyz_point_val+"' ")
                xyz_point_num += 1
            f.write(")\n") 
            f.write("\n")
            f.write("module add openmpi\n")
            f.write("module add orca/4.0.1.2\n")
            f.write("\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("\n")
            f.write("cd r_${RADII[$PBS_ARRAYID]}\n")
            f.write("/share/apps/orca/4.0.1.2/bin/orca "+args.filename+"_freq_${XYZ_VAL[$PBS_ARRAYID]}.inp"+ \
                    " > "+args.filename+"_freq_${XYZ_VAL[$PBS_ARRAYID]}.out \n" )
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("\n")

#All the frequency calculations should be done now

def write_processing_pbs_script(args):
    with open(args.filename+'_process.pbs','w') as f:
        if args.cluster == "katana":
            f.write("#!/bin/bash\n")
            f.write("\n")
            f.write("#PBS -l nodes=1:ppn=1\n") #Processing shouldn't take many resources
            f.write("#PBS -l vmem=1gb\n")
            f.write("#PBS -l walltime=00:10:00\n")
            f.write("#PBS -j oe\n") 
            f.write("\n")
            f.write("#PBS -t 0-"+str(args.increments-1)+"\n") #Arrays start at 0
            f.write("\n")
            f.write("RADII=(") 
            point_radius = args.r_start #I should really be passing these no recomputing
            increment_step = (args.r_end - args.r_start)/(args.increments-1)
            for point in range(args.increments):
                pnt_rad_formatted=str("{:.2f}".format(point_radius)) # trailing zerso 2 d.p.
                f.write("'"+pnt_rad_formatted+"' ")
                point_radius += increment_step
            f.write(")\n") 
            f.write("\n")
            xyz_point_num = 1 #ORCA xyzfiles numbered from 1
            f.write("XYZ_VAL=(") 
            for point in range(args.increments):
                xyz_point_val=str("{0:03}".format(xyz_point_num)) #0 filled formatting
                f.write("'"+xyz_point_val+"' ")
                xyz_point_num += 1
            f.write(")\n") 
            f.write("\n")
            f.write("module load python/2.7.9")
            f.write("\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("\n")
            f.write("freqs_to_multiwell.py -f r_${RADII[$PBS_ARRAYID]}/"+args.filename+"_freq_${XYZ_VAL[$PBS_ARRAYID]}.out"+ \
                    " > r_${RADII[$PBS_ARRAYID]}/"+args.filename+".${XYZ_VAL[$PBS_ARRAYID]}.freqs\n")
            f.write("rotors_determine.py -u GHz -B B2d Bk -f r_${RADII[$PBS_ARRAYID]}/"+args.filename+".${XYZ_VAL[$PBS_ARRAYID]}.xyz"+ \
                    " > r_${RADII[$PBS_ARRAYID]}/"+args.filename+".${XYZ_VAL[$PBS_ARRAYID]}.rotors\n")
            f.write("\n")

def write_freq_input_file(args,dir,xyz_point_val): 
    with open(dir+'/'+args.filename+'_freq_'+xyz_point_val+'.inp','w') as f:
        args.job_type = 'freq'
        f.write('# '+args.name+' frequency calculation\n')
        f.write('! '),
        f.write('PAL' + str(args.ncpus) + ' ')
        if args.restricted == 'U':
                f.write('UKS '),
        elif args.restricted == 'R':
                f.write('RKS '),
        if args.approx != 'none':
                f.write('RI-'+args.method+' ')
        else:
                f.write(args.method+' '),
        f.write(args.basis_set+' '),
        if args.dispersion != 'none':
                f.write(args.dispersion+' ')
        if args.comments != 'off':
            f.write('# Number of CPUs and level of theory'),
        f.write("\n")
        f.write('! '),
        if args.approx != 'none':
                f.write(args.approx+' '),
        if args.approx == 'RIJK':
                f.write('def2/JK '),
        if args.method in ['RI-MP2','RI-B2GP-PLYP', 'RI-B2-PLYP', 'RI-B2K-PLYP', 'RI-B2T-PLYP']:
            args.ri_aux = 'on'
        if args.ri_aux == 'on':
                f.write(args.basis_set+'/C '),
        if args.comments != 'off':
                f.write('# Density fitting/resoultion-of-identity speed-ups'),
        f.write('MOREAD ')
        f.write("\n")
        f.write('! '),
        f.write(args.scf_tol+'SCF '),
        f.write('Grid'+str(args.grid)+' '),
        if args.comments != 'off':
                f.write('# Convergence tolerence settings'),
        f.write("\n")
        f.write('! '),
        if args.method in ['B2GP-PLYP', 'B2-PLYP', 'B2K-PLYP', 'B2T-PLYP']: #analytical Hessians not implemented for double-hybrids
                f.write('num')
        f.write(args.job_type+' '),
        f.write("\n")
        f.write('%moinp "'+args.filename+'.'+xyz_point_val+'.gbw'+'"')
        f.write("\n")
        f.write("\n")
        f.write('*xyzfile '+str(args.charge)+' '+str(args.multiplicity)+' '+args.filename+'.'+xyz_point_val+'.xyz')
        f.write("\n")

if __name__ == '__main__':
    main()
