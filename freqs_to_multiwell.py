#!/usr/bin/env python
import argparse
import re

def main(): 
    args = parse_args()
    fname = args.filename
    cleaned_vibs = get_vibs(fname)

    print("{0}\tHAR\tghz").format(len(cleaned_vibs))
    for mode, vib_freq in enumerate(cleaned_vibs):
        print("{0}\tvib\t{1}\t\t0\t1").format(mode+1,vib_freq)

def parse_args():
    parser = argparse.ArgumentParser(description='Take the frequencies stored in the quantum chemistry output file and put them in MultiWell format')
    parser.add_argument('-f',"--filename", help="File to read frequency information from. Currently only ORCA .out files supported", type=str, required=True)
    args = parser.parse_args()
    return args

def get_vibs(fname):
    begin_vibs_section="VIBRATIONAL FREQUENCIES"
    end_vibs_section="NORMAL MODES"
    raw_vibs_section=[]
    is_vibs_section = False

    f = open(fname)
    for line in f:
        if end_vibs_section in line:
            is_vibs_section = False
        if is_vibs_section:
            raw_vibs_section.append(line.strip())
        if begin_vibs_section in line:
            is_vibs_section = True
    
    cleaned_vibs = clean_vibs_string(raw_vibs_section)

    return cleaned_vibs    

def clean_vibs_string(raw_vibs_section):
    cleaned_vibs = []

    for line in  raw_vibs_section:
        value_match = re.search('\d*\.\d*',line) #take only decimals
        if value_match != None:
            vib_value = value_match.group(0) #get value of match
            if vib_value != '0.00': #remove projected out translations/rotations
                cleaned_vibs.append(vib_value)

    return cleaned_vibs 

if __name__ == '__main__':
        main()
