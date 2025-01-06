import glob
import sys

atom_line_length = [10, 10, 10, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
def format_atom_line(arr):
    ret = ''
    if arr[3] == 'R' or arr[3] == '*':     
        for i in range(len(atom_line_length)):
            if i == 3:
                # ret += ' ' + arr[i]
                ret += ' ' + 'e*'
            elif i == 4:
                l = len(arr[i])
                ret += ' ' * (atom_line_length[i] - 1 - l) + arr[i]
            else:
                l = len(arr[i])
                ret += ' ' * (atom_line_length[i] - l) + arr[i]
    else:
        for i in range(len(atom_line_length)):
            l = len(arr[i])
            ret += ' ' * (atom_line_length[i] - l) + arr[i]
    return ret

bond_line_length = [3, 3, 3, 3, 3, 3, 3]
def format_bond_line(arr):
    ret = ''
    for i in range(len(bond_line_length)):
        l = len(arr[i])
        ret += ' ' * (bond_line_length[i] - l) + arr[i]
    return ret

def main(argv):
    if len(argv) < 2:
        print("Please input SDF file name and output file name")
        print("e.g.) python polymer_converter_sdf_to_sdf.py input.sdf output")
        sys.exit()

    input_name = argv[1]
    output_name = argv[2]

    write_file = open(output_name + '.sdf', 'w')
    
    read_file = open(input_name, 'r')

    while 1:
        molecule_name = read_file.readline().rstrip()
        if molecule_name == '':
            break
        software_name = read_file.readline().rstrip()
        comment = read_file.readline().rstrip()
        count_line = read_file.readline().rstrip()
        atom_num = int(count_line[ : 3])
        bond_num = int(count_line[3 : 6])

        atoms = []
        atom_line_length = [10, 10, 10, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
        for i in range(atom_num):
            atom = read_file.readline().rstrip()
            tmp = []
            pos = 0
            if 'R' in atom or '*' in atom:
                for j in range(len(atom_line_length)):
                    if j == 3:
                        tmp.append(atom[pos : pos + atom_line_length[j] + 1].strip())
                        pos += atom_line_length[j] + 1
                    elif j == 4:
                        tmp.append(atom[pos : pos + atom_line_length[j] - 1].strip())
                        pos += atom_line_length[j] - 1
                    else:
                        tmp.append(atom[pos : pos + atom_line_length[j]].strip())
                        pos += atom_line_length[j]
            else:
                for j in range(len(atom_line_length)):
                    tmp.append(atom[pos : pos + atom_line_length[j]].strip())
                    pos += atom_line_length[j]
            atoms.append(tmp)

        bonds = []
        bond_line_length = [3, 3, 3, 3, 3, 3, 3]
        for i in range(bond_num):
            bond = read_file.readline().rstrip()
            tmp = []
            pos = 0
            for j in range(len(bond_line_length)):
                tmp.append(bond[pos : pos + bond_line_length[j]].strip())
                pos += bond_line_length[j]
            bonds.append(tmp)
            
        others = []
        charge_info_atom = []
        charge_info_chg = []

        while 1:
            tmp = read_file.readline().rstrip()

            ## need to deal with "M  CHG""
            if "M  CHG" in tmp:
                k = int(tmp[6:9])
                for i in range(k):
                    charge_info_atom.append(int(tmp[9 + 8 * i:13 + 8 * i]))
                    charge_info_chg.append(int(tmp[13 + 8 * i:17 + 8 * i]))
            else:
                others.append(tmp)

            if tmp == '$$$$':
                break

        e_star_num_1 = -1
        e_star_num_2 = -1
        for i, atom in enumerate(atoms, 1):
            if e_star_num_1 == -1:
                if atom[3] == 'R' or atom[3] == '*':
                    e_star_num_1 = i
            else:        
                if atom[3] == 'R' or atom[3] == '*':
                    e_star_num_2 = i

        real_n = atom_num if e_star_num_2 == -1 else atom_num - 1

        if e_star_num_2 != -1:
            molecule_name += "_p"
        else:
            molecule_name += "_m"

        write_file.write(molecule_name + '\n')
        write_file.write(software_name + '\n')
        write_file.write(comment + '\n')
        count_line = ' ' * (3 - len(str(real_n))) + str(real_n) + count_line[3 : ]
        write_file.write(count_line + '\n')

        for i, atom in enumerate(atoms, 1):
            if i == e_star_num_2:
                continue
            else:
                write_file.write(format_atom_line(atom) + '\n')

        for bond in bonds:
            if int(bond[0]) == e_star_num_2:
                bond[0] = str(e_star_num_1)
            elif int(bond[1]) == e_star_num_2:
                bond[1] = str(e_star_num_1)

            if int(bond[0]) > e_star_num_2 and e_star_num_2 != -1:
                bond[0] = str(int(bond[0]) - 1)
            if int(bond[1]) > e_star_num_2 and e_star_num_2 != -1:
                bond[1] = str(int(bond[1]) - 1)
            
            write_file.write(format_bond_line(bond) + '\n')

        if len(charge_info_atom) != 0:
            ### need to deal with charge information
            charge_line = "M  CHG"
            charge_num = int(len(charge_info_atom))
            charge_line += " " * (3 - len(str(charge_num))) + str(charge_num)
            for (atom, chg) in zip(charge_info_atom, charge_info_chg):
                if atom > e_star_num_2 and e_star_num_2 != -1:
                    atom -= 1
                charge_line += " " * (4 - len(str(atom))) + str(atom) + " " * (4 - len(str(chg))) + str(chg)
            # print(molecule_name, charge_line)
            write_file.write(charge_line + '\n')

        for other in others:
            write_file.write(other + '\n')

    read_file.close()

    write_file.close()

main(sys.argv)