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

# parser class that should read list of pdbs.
# It first collects pdb files and then parse them into df_atom, df_helix, and df_sheet
# df_atom contains all atom information including each of their coordinates
# df_helix contains helix information
# df_sheet contains sheet information
class pdb_parser():
    def __init__(self, index_to_break: int=-1, flag=True):
        self.df_atom = pd.DataFrame()
        self.df_sheet = pd.DataFrame()
        self.df_helix = pd.DataFrame()
        self.filename_list_pdb = 'cullpdb_pc30_res3.0_R1.0_d191017_chains18877.gz'
        self.df_pdb_list = self.parse_list_pdb(self.filename_list_pdb)
        self.pdb_dir = './pdb_data'
        self.l_atom = list()
        self.l_sheet = list()
        self.l_helix = list()

        # using multithreading to download pdb files...
        self.download_all_pdb(index_to_break)
        if index_to_break == -1:
            print('parsing all pdb in {}'.format(self.filename_list_pdb))
        else:
            print('parsing first {} proteins in {}'.format(index_to_break, self.filename_list_pdb))

        for index, pdb in enumerate(list(self.df_pdb_list['IDs'])):
            if index == index_to_break:
                break
            # self.process_pdb(pdb)
            # flushing out list for memory...
            if index%500 == 0 and index != 0:
                print('completed parsing {} pdbs'.format(index))
                self.df_atom = pd.concat([self.df_atom, pd.read_json(json.dumps(self.l_atom))], ignore_index=True)
                self.df_helix = pd.concat([self.df_helix, pd.read_json(json.dumps(self.l_helix))], ignore_index=True)
                self.df_sheet = pd.concat([self.df_sheet, pd.read_json(json.dumps(self.l_sheet))], ignore_index=True)
                self.l_atom = list()
                self.l_helix = list()
                self.l_sheet = list()

            self.process_pdb_new(pdb)
        
        self.df_atom = pd.concat([self.df_atom, pd.read_json(json.dumps(self.l_atom))], ignore_index=True).astype('str')
        self.df_helix = pd.concat([self.df_helix, pd.read_json(json.dumps(self.l_helix))], ignore_index=True).astype('str')
        self.df_sheet = pd.concat([self.df_sheet, pd.read_json(json.dumps(self.l_sheet))], ignore_index=True).astype('str')
        
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

    def process_pdb_new(self, pdb_name):
        protein_name = pdb_name[:-1]
        protein_chain = pdb_name[-1]
        # print('parsing {}...'.format(pdb_name))

        with open(self.pdb_dir + '/{}.pdb'.format(protein_name), 'r') as f:
            for line in f.readlines():
                dict_parsed_atom = self.parse_pdb_data(line, 'ATOM', pdb_name)
                dict_parsed_helix = self.parse_pdb_data(line, 'HELIX', pdb_name)
                dict_parsed_sheet = self.parse_pdb_data(line, 'SHEET', pdb_name)

                if dict_parsed_atom:
                    if dict_parsed_atom['chain_id'] == protein_chain:
                        self.l_atom += [dict_parsed_atom]

                elif dict_parsed_helix:
                    if dict_parsed_helix['init_chain_id'] == protein_chain or \
                            dict_parsed_helix['end_chain_id'] == protein_chain:
                        self.l_helix += [dict_parsed_helix]

                elif dict_parsed_sheet:
                    if dict_parsed_sheet['cur_chain_id'] == protein_chain:
                        self.l_sheet += [dict_parsed_sheet]


    def process_pdb(self, pdb_name):
        protein_name = pdb_name[:-1]
        protein_chain = pdb_name[-1]
        # print('parsing {}...'.format(pdb_name))

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

    def download_all_pdb(self, index_to_break):
        from multiprocessing import Pool
        p = Pool(2)
        p.map(self.download_pdb, list(self.df_pdb_list['IDs'])[:index_to_break]) 
        
    @staticmethod
    def download_pdb(protein_name):
        import requests
        pdb_dir = './pdb_data'
        protein_name = protein_name[:-1]
        if not(os.path.exists(pdb_dir + '/{}.pdb'.format(protein_name))):
#             print('pdb file for {} not found. Downloading from protein data bank...'.format(protein_name))
            
            with open(pdb_dir + '/{}.pdb'.format(protein_name), 'wb') as f:
                url = 'https://files.rcsb.org/view/{}.pdb'.format(protein_name)
                r = requests.get(url)
                f.write(r.content)

    @staticmethod
    def parse_list_pdb(file_name: str):
        with open(file_name, 'r') as file_list_pdb:
            df_list_pdb = pd.read_csv(file_list_pdb, delimiter= ' ', skipinitialspace=True)
            return df_list_pdb

if __name__ == '__main__':
    import time
    time_now = time.time()
    parser = pdb_parser(2000)
    parser.print_stats()
    print('----{}s----'.format(time.time() - time_now))