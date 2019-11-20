import requests
import json
import numpy as np
import pandas as pd
import os

format_spacing = {
    'ATOM' : { 
        'spacing': [6, 11, 16, 17, 20, 22, 26, 27, 38, 46, 54, 60, 66, 78, 80],
        'label': ['record_name', 'serial_number', 
                  'atom_name', 'alt_loc', 'res_name', 
                  'chain_id', 'res_seq', 'iCode', 
                  'x', 'y', 'z', 'occupancy', 'tempFactor', 
                  'element', 'charge']
    },
    'HELIX' : {
        'spacing': [6, 10, 14, 18, 20, 25, 26, 30, 32, 37, 38, 40, 70, 76],
        'label': ['record_name', 'serial_number', 'helix_id', 'init_res_name',
                  'init_chain_id', 'init_seq_num', 'init_iCode', 'end_res_name',
                  'end_chain_id', 'end_seq_num', 'end_iCode', 'helix_class',
                  'comment', 'length']
    },
    'SHEET' : {
        'spacing': [6, 10, 14, 16, 20, 22, 26, 27, 31, 33, 37, 38, 40, 45, 48, 50, 54, 55, 60, 63, 65, 69, 70],
        'label': ['record_name', 'strand', 'sheet_id', 'num_strands', 'init_res_name',
                  'init_chain_id', 'init_seq_num', 'init_iCode', 'end_res_name',
                  'end_chain_id', 'end_seq_num', 'end_iCode', 'sense', 'cur_atom', 'cur_res_name',
                  'cur_chain_id', 'cur_res_seq', 'cur_iCode', 'prev_atom', 'prev_res_name', 
                  'prev_chain_id', 'prev_res_seq', 'rev_iCode']
    }
}

atom_backchain = ['N', 'CA', 'C', 'O']

class pdb_parser():
    def __init__(self, index_to_break: int=-1):
        self.df_atom = pd.DataFrame()
        self.df_sheet = pd.DataFrame()
        self.df_helix = pd.DataFrame()
        self.filename_list_pdb = 'cullpdb_pc30_res3.0_R1.0_d191017_chains18877.gz'
        self.df_pdb_list = self.parse_list_pdb(self.filename_list_pdb)
        self.pdb_dir = './pdb_data'

        if index_to_break == -1:
            print('parsing all pdb in {}'.format(self.filename_list_pdb))
        else:
            print('parsing first {} proteins in {}'.format(index_to_break, self.filename_list_pdb))

        for index, pdb in self.df_pdb_list.iterrows():
            if index == index_to_break:
                break
            pdb_id = pdb['IDs']
            print('parsing {}...'.format(pdb_id))
            self.process_pdb(pdb_id, pdb['length'])

    def parse_list_pdb(self, file_name: str):
        with open(file_name, 'r') as file_list_pdb:
            df_list_pdb = pd.read_csv(file_list_pdb, delimiter= ' ', skipinitialspace=True)
            return df_list_pdb

    def parse_pdb_data(self, str_pdb: str, keyword: str, pdb_name: str):

        if keyword in format_spacing.keys() and str_pdb.split(' ')[0] == keyword:
            l_label = format_spacing[keyword]['label']
            l_spacing = format_spacing[keyword]['spacing']
            if len(l_label) != len(l_spacing):
                raise Exception('length of label and spacing for {} not matching'.format(keyword))
            else:
                tmp_dict = dict()
                for i in range(len(l_spacing)):
                    if i == 0:
                        tmp_dict[l_label[i]] = str_pdb[0:l_spacing[i]].strip()
                    else:
                        tmp_dict[l_label[i]] = str_pdb[l_spacing[i-1]:l_spacing[i]].strip()
                
                tmp_dict['protein_name'] = pdb_name
                return tmp_dict

    def process_pdb(self, pdb_name, len_atom):
        import requests
        protein_name = pdb_name[:-1]
        protein_chain = pdb_name[-1]
        # print('Beginning {} pdb file download with requests'.format(pdb_name))

        if not(os.path.exists(self.pdb_dir + '/{}.pdb'.format(protein_name))):
            print('pdb file for {} not found. Downloading from protein data bank...'.format(protein_name))
            with open(self.pdb_dir + '/{}.pdb'.format(protein_name), 'wb') as f:
                url = 'https://files.rcsb.org/view/{}.pdb'.format(protein_name)
                r = requests.get(url)
                f.write(r.content)

        with open(self.pdb_dir + '/{}.pdb'.format(protein_name), 'r') as f:
            for line in f.readlines():
                dict_parsed_atom = self.parse_pdb_data(line, 'ATOM', pdb_name)
                dict_parsed_helix = self.parse_pdb_data(line, 'HELIX', pdb_name)
                dict_parsed_sheet = self.parse_pdb_data(line, 'SHEET', pdb_name)

                if dict_parsed_atom:
                    if dict_parsed_atom['chain_id'] == protein_chain:
                        self.df_atom = self.df_atom.append(dict_parsed_atom, ignore_index=True)

                elif dict_parsed_helix:
                    if dict_parsed_helix['init_chain_id'] == protein_chain or \
                            dict_parsed_helix['end_chain_id'] == protein_chain:
                        self.df_helix = self.df_helix.append(dict_parsed_helix, ignore_index=True)

                elif dict_parsed_sheet:
                    if dict_parsed_sheet['cur_chain_id'] == protein_chain:
                        self.df_sheet = self.df_sheet.append(dict_parsed_sheet, ignore_index=True)

    def print_stats(self):
        # print(self.df_atom[:5])
        print('# atoms:{} # helices:{} # sheets: {}'.format(str(len(self.df_atom)),
                                                            str(len(self.df_helix)), 
                                                            str(len(self.df_sheet))))

                
if __name__ == '__main__':
    parser = pdb_parser(3)
    # parser.process_pdb('3UTS')
    parser.print_stats()