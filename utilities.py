import pandas as pd
import numpy as np
from scrape_pdb import pdb_parser

class pdb_utilities:
    def __init__(self, df_atom, df_helix, df_sheet):
        self.df_atom = df_atom
        self.df_helix = df_helix
        self.df_sheet = df_sheet
        self.df_aa = pd.DataFrame()

    def find_coordinates_atom(self, protein_name: str, atom_name: str):
        atom_name.split('.')[0]
        found_atom = self.df_atom[(self.df_atom['res_seq'] == atom_name.split('.')[0]) &
                            (self.df_atom['atom_name'] == atom_name.split('.')[1])]
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
            print('Type Error occurred!! {} {} {}'.format(protein_name, seq_aa, typ))


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

    def build_ramachandran_aa(self, res_name):
        df_aa = self.df_atom[self.df_atom['res_name'] == res_name][['protein_name', 'res_seq']].drop_duplicates(keep='first')
        self.df_aa = df_aa
        df_aa['psi'] = df_aa.apply(lambda x: self.calculate_angle(x['protein_name'], x['res_seq'], 'psi'), axis=1)
        df_aa['phi'] = df_aa.apply(lambda x: self.calculate_angle(x['protein_name'], x['res_seq'], 'phi'), axis=1)
        
        return df_aa

if __name__ == '__main__':
    pdb_parser = pdb_parser(10)
    pdb_utilities = pdb_utilities(pdb_parser.df_atom, 
                                  pdb_parser.df_helix, 
                                  pdb_parser.df_sheet)
    # print(pdb_parser.df_atom[pdb_parser.df_atom['protein_name'] == '12AS'])
    # print(pdb_utilities.calculate_angle('12ASA', '327', 'psi'))
    pdb_utilities.build_ramachandran_aa('VAL')
