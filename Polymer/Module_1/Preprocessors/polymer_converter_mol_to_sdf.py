import glob
import sys

def convert_mol_file(read_file, write_file):
    
    atom_line_length = [10, 10, 10, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    def format_atom_line(arr):
        ret = ''
        if arr[3] == 'e*':        
            for i in range(len(atom_line_length)):
                if i == 3:
                    ret += ' ' + arr[i]
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

    file_data = read_file.readlines()

    molecule_name = file_data[0].rstrip()
    software_name = file_data[1].rstrip()
    comment = file_data[2].rstrip()
    count_line = file_data[3].rstrip()
    atom_num = int(count_line[ : 3])
    bond_num = int(count_line[3 : 6])

    atoms = []
    atom_line_length = [10, 10, 10, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    for i in range(atom_num):
        atom = file_data[4 + i].rstrip()
        tmp = []
        pos = 0
        if 'e*' in atom:
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
        bond = file_data[4 + atom_num + i].rstrip()
        tmp = []
        pos = 0
        for j in range(len(bond_line_length)):
            tmp.append(bond[pos : pos + bond_line_length[j]].strip())
            pos += bond_line_length[j]
        bonds.append(tmp)

    others = file_data[4 + atom_num + bond_num : ]

    write_file.write(molecule_name + '\n')
    write_file.write(software_name + '\n')
    write_file.write(comment + '\n')
    count_line = ' ' * (3 - len(str(atom_num - 1))) + str(atom_num - 1) + count_line[3 : ]
    write_file.write(count_line + '\n')

    e_star_num_1 = -1
    e_star_num_2 = -1
    for i, atom in enumerate(atoms, 1):
        if e_star_num_1 == -1:
            if atom[3] == 'e*':
                e_star_num_1 = i
        else:        
            if atom[3] == 'e*':
                e_star_num_2 = i

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
        write_file.write(format_bond_line(bond) + '\n')

    for other in others:
        write_file.write(other)

def main(argv):
    if len(argv) < 2:
        print("Please input folder name and output file name")        
        print("e.g.) python polymer_converter_mol_to_sdf.py folder/name output")
        sys.exit()

    input_name = argv[1]
    output_name = argv[2]

    file_addr = input_name+'/*.mol'
    file_list = glob.glob(file_addr)

    write_file=open(output_name + '.sdf', 'w')
    
    for file_name in file_list:
        read_file=open(file_name, 'r')

        convert_mol_file(read_file, write_file)

        write_file.write('$$$$\n')            
        read_file.close()

    write_file.close()

main(sys.argv)