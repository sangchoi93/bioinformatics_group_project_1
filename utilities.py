import pandas as pd
import numpy as np
import json
from scrape_pdb import pdb_parser
from matplotlib import pyplot as plt

# utilities to analyze the pdb data parsed out by pdb_parser class in scrape_pdb.py
class pdb_utilities:
    def __init__(self, df_atom, df_helix, df_sheet):
        self.df_atom = df_atom
        self.df_helix = df_helix
        self.df_sheet = df_sheet
        self.df_aa = pd.DataFrame()
        self.type_helix = {
            '1': 'right-handed alpha',
            '2': 'right-handed omega',
            '3': 'right-handed psi',
            '4': 'right-handed gama',
            '5': 'right-handed 3-10',
            '6': 'left-handed alpha',
            '7': 'left-handed omega',
            '8': 'left-handed gamma',
            '9': '2-7 ribbon/hex',
            '10': 'polyproline'
        }

    def find_coordinates_atom(self, protein_name: str, atom_name: str):
        l_atom_name = atom_name.split('.')
        found_atom = self.df_atom[(self.df_atom['res_seq'] == l_atom_name[0]) &
                            (self.df_atom['atom_name'] == l_atom_name[1])]
        if protein_name:
            found_atom = found_atom[found_atom['protein_name'] == protein_name]
        
        # print('{}\'s coordinate is ({}, {}, {})'.format(atom_name, 
        #                                                 float(found_atom['x']), 
        #                                                 float(found_atom['y']), 
        #                                                 float(found_atom['z'])))
        return {'x': float(found_atom['x']), 
                'y': float(found_atom['y']), 
                'z': float(found_atom['z'])}

    @staticmethod
    def listify_coordinates(d: dict):
        return np.array([d['x'], d['y'], d['z']])

    def calculate_angle(self, protein_name: str, seq_aa:str, typ: str):
        try:
            if typ == 'phi':
                coord_c = self.find_coordinates_atom(protein_name, str(int(seq_aa)-1) + '.C')
                coord_n = self.find_coordinates_atom(protein_name, seq_aa + '.N')
                coord_ca = self.find_coordinates_atom(protein_name, seq_aa + '.CA')
                coord_c_2 = self.find_coordinates_atom(protein_name, seq_aa+'.C')
                return self.calculate_dihedral([self.listify_coordinates(coord_c), 
                                        self.listify_coordinates(coord_n), 
                                        self.listify_coordinates(coord_ca), 
                                        self.listify_coordinates(coord_c_2)])
            elif typ == 'psi':
                coord_n = self.find_coordinates_atom(protein_name, seq_aa + '.N')
                coord_ca = self.find_coordinates_atom(protein_name, seq_aa + '.CA')
                coord_c = self.find_coordinates_atom(protein_name, seq_aa + '.C')
                coord_n_2 = self.find_coordinates_atom(protein_name, str(int(seq_aa)+1) + '.N')
                return self.calculate_dihedral([self.listify_coordinates(coord_n), 
                                                self.listify_coordinates(coord_ca), 
                                                self.listify_coordinates(coord_c), 
                                                self.listify_coordinates(coord_n_2)])
        
        # specify float exception in find_coordinate function if atom is not found
        except FloatingPointError:
            print('Floating point error occurred!! {} {} {}'.format(protein_name, seq_aa, typ))
            return None
        except TypeError:
            print('atom not found!', protein_name, seq_aa, typ)
            return None

    @staticmethod
    def calculate_dihedral(list_coords):
        v1 = list_coords[0] - list_coords[1]
        v2 = list_coords[2] - list_coords[1]
        v3 = list_coords[3] - list_coords[2]

        v1xv2 = np.cross(v1, v2)
        v2xv3 = np.cross(v3, v2)

        v1xv2_x_v2xv3 = np.cross(v1xv2, v2xv3)

        y = np.dot(v1xv2_x_v2xv3, v2)/np.linalg.norm(v2)
        x = np.dot(v1xv2, v2xv3)

        return round(np.degrees(np.arctan2(y, x)), 3)        

    def plot_ramachandran(self, df:pd.DataFrame):
        plt.grid()
        plt.scatter(df['psi'], df['phi'], alpha=0.5)

    def build_ramachandran_aa(self, res_name):
        df_aa = self.df_atom[self.df_atom['res_name'] == res_name][['protein_name', 'res_seq']].drop_duplicates(keep='first')
        df_aa['psi'] = df_aa.apply(lambda x: self.calculate_angle(x['protein_name'], x['res_seq'], 'psi'), axis=1)
        df_aa['phi'] = df_aa.apply(lambda x: self.calculate_angle(x['protein_name'], x['res_seq'], 'phi'), axis=1)
        self.plot_ramachandran(df_aa)
        return df_aa

    def build_ramachandran_helices(self):
        df_helices_ramanchandran = pd.DataFrame()
        l_dict_ramanchandra = list()

        for helix in self.df_helix[['protein_name', 'helix_class', 'init_seq_num', 'end_seq_num', 'init_chain_id', 'end_chain_id']].values:
            if helix[4] != helix[5]:
                print('end_chain_id and init_chain_id not matching!!!!!', helix[4], helix[5])
            else:
                for i in range(int(helix[2]), int(helix[3])+1):
                    dict_helix = {
                        'helix_type': self.type_helix[helix[1]],
                        'helix_type_code': helix[1],
                        'protein_name': helix[0],
                        'aa_seq_num': str(i),
                        'psi': self.calculate_angle(helix[0], str(i), 'psi'),
                        'phi': self.calculate_angle(helix[0], str(i), 'phi')
                    }
                    l_dict_ramanchandra += [dict_helix]
            
        df_helices_ramanchandran = pd.read_json(json.dumps(l_dict_ramanchandra))
        return df_helices_ramanchandran
   
if __name__ == '__main__':
    pdb_parser = pdb_parser(1)
    pdb_utilities = pdb_utilities(pdb_parser.df_atom, 
                                  pdb_parser.df_helix, 
                                  pdb_parser.df_sheet)
    # print(pdb_parser.df_atom[pdb_parser.df_atom['protein_name'] == '12AS'])
    # print(pdb_utilities.calculate_angle('12ASA', '327', 'psi'))
    pdb_utilities.build_ramachandran_aa('VAL')
