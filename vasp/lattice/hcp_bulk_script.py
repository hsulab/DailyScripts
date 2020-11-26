# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 20:52:51 2017

@author: Alexander
"""
#!/usr/bin/env python3
###
def get_lattic_constant(file_name = r'./POSCAR'):
    with open(file_name, 'r') as file_open:
        lines = file_open.readlines()
        n = 0
        i = 0
        a = 0
        b = 0
        c = 0
        for line in lines:
            n = n + 1
            if not line:
                break
            else:
                if n == 3: 
                    i = line.split()[0]
                    a = line.split()[1]
                if n == 4: 
                    b = line.split()[1]
                if n == 5: 
                    c = line.split()[2]
    return [i,a,b,c]
# v=$(echo "scale=7;114.634/ 14.61^3 "| bc)
# 2.41 2.45 2.49 2.53 2.57 2.61 2.65 2.69 2.73 2.77 2.81
# a=$(echo "scale=15;$i/2.60673646539116* -1.505 "| bc)
# b=$(echo "scale=15;$i/2.60673646539116*3.01 "| bc) 
# c=$(echo "scale=15;$i/2.60673646539116*14.61 "| bc)
def string_switch(file_name,old_words,new_words,switch_option=1):
    with open(file_name, "r") as file_open:
        #readlines in the form of list
        lines = file_open.readlines()
 
    with open(file_name, "w") as file_write:
        # define a number to indicate where it is in the file
        line_number = 0
        # default option, only replace the first 
        if switch_option == 1:
            for line in lines:
                if old_words in line:
                    line = line.replace(old_words,new_words)
                    file_write.write(line)
                    line_number += 1
                    break
                file_write.write(line)
                line_number += 1
            # output the residue file
            for i in range(line_number,len(lines)):
                file_write.write(lines[i])
        # global match and replace
        elif switch_option == 'g':
            for line in lines:
                if old_words in line:
                    line = line.replace(old_words, new_words)
                file_write.write(line)
###
def create_bulk_script(vasp_bulk_script, POSCAR):
    lattice_constant = get_lattic_constant(POSCAR)
    #string_switch(r'vasp_bulk.script.', r'volume_location', str(), 1)
    ### POSCAR
    with open(r'./POSCAR', 'r') as file_open:
        string_switch( vasp_bulk_script , r'POSCAR_location', ''.join(file_open.readlines()), 1)
    ###
    string_switch(vasp_bulk_script , str(lattice_constant[0]), r'$i', 1)
    string_switch(vasp_bulk_script , str(lattice_constant[1]), r'$a', 1)
    string_switch(vasp_bulk_script , str(lattice_constant[2]), r'$b', 1)
    string_switch(vasp_bulk_script , str(lattice_constant[3]), r'$c', 1)
    ###
    i_range = []
    for i in range(11):
        epi = float(lattice_constant[0])/100.
        i_range.append(float(lattice_constant[0])+(-5+i)*epi)
    i_range = [ str(round(i,15)) for i in i_range]
    string_switch(vasp_bulk_script , r'i_range_location', ' '.join(i_range), 1)
    ### i, a, b, c
    string_switch(vasp_bulk_script , r'i_location', str(lattice_constant[0]), r'g')
    string_switch(vasp_bulk_script , r'a_location', str(lattice_constant[1]), 1)
    string_switch(vasp_bulk_script , r'b_location', str(lattice_constant[2]), 1)
    string_switch(vasp_bulk_script , r'c_location', str(lattice_constant[3]), r'g')
    ###
def main():
    create_bulk_script(r'./vasp_bulk.script', r'./POSCAR')
    print(get_lattic_constant())

if __name__ == "__main__":
    main()
