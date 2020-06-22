import argparse

eh_to_kjmol = 2625.499638

def main():
    args = parse_args()
    delh_list = get_delhs(args)    
  
    #Write out the relative energy in kJ/mol in each point's folder for each grabbing by the processing file
    for idx, entry in enumerate(delh_list):
        point_radius = float(entry[0])
        dir = ("r_" + "{0:.2f}").format(point_radius)
        point_index = str(("{0:03}").format(idx+1))
        with open(dir+"/"+args.filename.replace('out','')+point_index+".ediff", 'w') as delh_f:
            delh_f.write(str(entry[1]))
        
def parse_args():
    parser = argparse.ArgumentParser(description="Take the relaxed scan job summary, and convert into 'delh' energy change in kJ/mol which can be read in at each point folder to construct a VTST MultiWell data file line.")
    parser.add_argument('-f',"--filename", help="Output file containing the relaxed scan information. Currently only ORCA .out files supported.", type=str, required=True)
    parser.add_argument('-re',"--reactant_energy", help="The reactant minimum energy all VTST points will be reported relative to, in Hartree", type=float, required=True) 
    args = parser.parse_args()
    return args

def get_delhs(args):
    begin_delh_section="The Calculated Surface using the 'Actual Energy'"
    end_delh_section="The Calculated Surface using the SCF energy"
    is_delh_section = False
    raw_surface_energies = []
    
    f = open(args.filename)
    for line in f:
        if end_delh_section in line:
            is_delh_section = False
        if is_delh_section:
            if line.strip(): #Get rid of whitespace lines
                raw_surface_energies.append(line.split()) 
        if begin_delh_section in line:
            is_delh_section = True
    
    delh_list =  convert_to_rel_kjmol(raw_surface_energies, args)

    return delh_list 

def convert_to_rel_kjmol(raw_surface_energies, args):
    for row in range(len(raw_surface_energies)):
        raw_surface_energies[row][1] = (float(raw_surface_energies[row][1]) - args.reactant_energy) * eh_to_kjmol  
    return raw_surface_energies

if __name__ == '__main__':
    main()
