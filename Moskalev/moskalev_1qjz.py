import re
import csv
import numpy


name_of_file = '1qjz.pdb'
ln10 = 2.30258509299
pH = 7
pKa = {
    'ASP': 3.71,
    'GLU': 4.15,
    'TYR': 10.10,
    'ARG': 12.10,
    'HIS': 6.04,
    'LYS': 10.67,
    'CYS': 8.14,
}
m_atom = {
    'N': 14.0067,
    'C': 12.011,
    'O': 15.999,
    'H': 1.00784,
    'S': 32.065,
}
keys_of_m_atom = ['N', 'C', 'O', 'H', 'S']
data_chem_elem = []
data_amino_acid = []

x_i = 0
y_i = 0
z_i = 0
x_i_q = 0
y_i_q = 0
z_i_q = 0
q_i_o = 0
m_i = 0
a = 0
b = 0


def main_A(name_file: str):
    slice_file(name_file)
    for line in csv.DictReader(open(name_file[:4] + '_sliced.csv')):
        global x_i, y_i, z_i, m_i, x_i_q, y_i_q, z_i_q, q_i_o, a, b
        if line['chem_elem'] not in data_chem_elem:
            data_chem_elem.append(line['chem_elem'])
        if line['amino_acid'] not in data_amino_acid:
            data_amino_acid.append(line['amino_acid'])
        if line['name_protein'] == 'A':
            try:
                x_i += float(line['X']) * m_atom[line['chem_elem']]
                y_i += float(line['Y']) * m_atom[line['chem_elem']]
                z_i += float(line['Z']) * m_atom[line['chem_elem']]
                q_i_o += q_i(line['amino_acid'])
                if line['chem_elem'] in keys_of_m_atom:
                    m_i += m_atom[line['chem_elem']]
                else:
                    m_i += 0
                x_i_q += float(line['X']) * q_i(line['amino_acid'])
                y_i_q += float(line['Y']) * q_i(line['amino_acid'])
                z_i_q += float(line['Z']) * q_i(line['amino_acid'])
                '''print('Mass gravity:{} {} {} Charge center: X = {} Y = {} Z = \
                     {} qi: {} name of acid: {}'.format(line['X'], line['Y'], \
                     line['Z'], x_i_q, y_i_q, z_i_q, q_i(line['amino_acid']), \
                     line['amino_acid']))'''
            except KeyError:
                continue

    a = (x_i/m_i, y_i/m_i, z_i/m_i)
    b = (x_i_q/q_i_o, y_i_q/q_i_o, z_i_q/q_i_o)
    #return {'mass_gravity': a, 'charge_center': b, 'all_acid': data_amino_acid}
    return {'mass_gravity': a, 'charge_center': b}


def main_B(name_file: str):
    slice_file(name_file)
    for line in csv.DictReader(open(name_file[:4] + '_sliced.csv')):
        global x_i, y_i, z_i, m_i, x_i_q, y_i_q, z_i_q, q_i_o, a, b
        if line['chem_elem'] not in data_chem_elem:
            data_chem_elem.append(line['chem_elem'])
        if line['amino_acid'] not in data_amino_acid:
            data_amino_acid.append(line['amino_acid'])
        if line['name_protein'] == 'B':
            try:
                x_i += float(line['X']) * m_atom[line['chem_elem']]
                y_i += float(line['Y']) * m_atom[line['chem_elem']]
                z_i += float(line['Z']) * m_atom[line['chem_elem']]
                q_i_o += q_i(line['amino_acid'])
                if line['chem_elem'] in keys_of_m_atom:
                    m_i += m_atom[line['chem_elem']]
                else:
                    m_i += 0
                x_i_q += float(line['X']) * q_i(line['amino_acid'])
                y_i_q += float(line['Y']) * q_i(line['amino_acid'])
                z_i_q += float(line['Z']) * q_i(line['amino_acid'])
                '''print('Mass gravity:{} {} {} Charge center: X = {} Y = {} Z = \
                     {} qi: {} name of acid: {}'.format(line['X'], line['Y'], \
                     line['Z'], x_i_q, y_i_q, z_i_q, q_i(line['amino_acid']), \
                     line['amino_acid']))'''
            except KeyError:
                continue

    a = (x_i/m_i, y_i/m_i, z_i/m_i)
    b = (x_i_q/q_i_o, y_i_q/q_i_o, z_i_q/q_i_o)
    #return {'mass_gravity': a, 'charge_center': b, 'all_acid': data_amino_acid}
    return {'mass_gravity': a, 'charge_center': b}


def main_C(name_file: str):
    slice_file(name_file)
    for line in csv.DictReader(open(name_file[:4] + '_sliced.csv')):
        global x_i, y_i, z_i, m_i, x_i_q, y_i_q, z_i_q, q_i_o, a, b
        if line['chem_elem'] not in data_chem_elem:
            data_chem_elem.append(line['chem_elem'])
        if line['amino_acid'] not in data_amino_acid:
            data_amino_acid.append(line['amino_acid'])
        if line['name_protein'] == 'C':
            try:
                x_i += float(line['X']) * m_atom[line['chem_elem']]
                y_i += float(line['Y']) * m_atom[line['chem_elem']]
                z_i += float(line['Z']) * m_atom[line['chem_elem']]
                q_i_o += q_i(line['amino_acid'])
                if line['chem_elem'] in keys_of_m_atom:
                    m_i += m_atom[line['chem_elem']]
                else:
                    m_i += 0
                x_i_q += float(line['X']) * q_i(line['amino_acid'])
                y_i_q += float(line['Y']) * q_i(line['amino_acid'])
                z_i_q += float(line['Z']) * q_i(line['amino_acid'])
                '''print('Mass gravity:{} {} {} Charge center: X = {} Y = {} Z = \
                     {} qi: {} name of acid: {}'.format(line['X'], line['Y'], \
                     line['Z'], x_i_q, y_i_q, z_i_q, q_i(line['amino_acid']), \
                     line['amino_acid']))'''
            except KeyError:
                continue

    a = (x_i/m_i, y_i/m_i, z_i/m_i)
    b = (x_i_q/q_i_o, y_i_q/q_i_o, z_i_q/q_i_o)
    #return {'mass_gravity': a, 'charge_center': b, 'all_acid': data_amino_acid}
    return {'mass_gravity': a, 'charge_center': b}


def slice_file(name_file: str):
    f = open(name_file)
    f_new = open(name_file[:4] + '_sliced.csv', 'w')
    f_new.write('ATOM/HETATM,number,atom_type,amino_acid,name_protein,' +
                'residue_number,X,Y,Z,occupancy,B-factor,chem_elem\n')
    for i in f:
        if re.match(r'ATOM ', i):
            f_new.write(repl_spac_to_comma(i))
    f.close()
    f_new.close()


def repl_spac_to_comma(stroke: str):
    return stroke.replace('      ', ' ').replace('     ', ' ')\
                    .replace('    ', ' ').replace('   ', ' ')\
                    .replace('  ', ' ').replace(' ', ',').replace('-', ',-') \
                    .replace(',,', ',')


def q_i(amino_acid: str):
    if amino_acid in ['ARG', 'HIS', "LYS"]:
        return 1/(1 + numpy.exp(ln10*(pH - pKa[amino_acid])))
    elif amino_acid in ['ASP', 'GLU', 'TYR']:
        return (-1)/(1 + numpy.exp(ln10*(-(pH - pKa[amino_acid]))))
    else:
        return (-1)/(1 + numpy.exp(ln10*(pH)))


print('Name of virus ' + name_of_file[:-4])
print('For protein A: ', end='')
print(main_A(name_of_file))
print('For protein B:', end='')
print(main_B(name_of_file))
print('For protein C:', end='')
print(main_C(name_of_file))

if '__init__' == '__main__':
    pass
