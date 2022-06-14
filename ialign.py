import subprocess
import time
from pathlib import Path
import linecache
import re
try:
    import pandas as pd
except:
    import os
    os.system('cmd /c "python -m pip install pandas"')

from numpy import log as ln
#
# Bash command arguments
#
# The first one is the path to the perl script of iAlign
# -p -> contains the pdb files "/NIL/lib"
# -w -> output files will be saved "/ialign_py/outputs_new"
# Output format 2 = detailed a = "2"
# list.lst contains the names of the pdb files of the library
# target.lst contains the names of the files we run against the library (here only 1r0r.pdb)
class ialign:
    def __init__(self):
        self.Bin=r'/home/jbrender/ialign/bin/ialign.pl'
        self.OutDir=r'/home/jbrender/ialign/iAlignOut'
        #self.Library=r'/home/jbrender/Interfaces'
        self.Library=r'/home/jbrender/interfaces-v2.0'
        #self.LibraryList=r'/home/jbrender/Interfaces/interfaces.lst'
        self.Library=r'/home/jbrender/interfaces-v2.0/interfaces.lst'
        self.TargetListDir=r'/home/jbrender/ialign'
        #self.PDBLib=r'/home/jbrender/Repaired'
        self.PDBLib=r'/home/jbrender/optimus/Repaired'
        self.OutTime=r'/home/jbrender/ialign/iAlignTime.txt'
        self.OutPDBDir=r'/home/jbrender/ialign/ParsedSKEMPI'
        
    # Function to change amino acid name i.e ALA -> A
    def changeAA(self,amino_acid):
        switcher = {
            "ALA": 'A',
            "ARG": 'R',
            "ASP": 'D',
            "ASN": 'N',
            "CYS": 'C',
            "GLN": 'Q',
            "GLY": 'G',
            "GLU": 'E',
            "HIS": 'H',
            "ILE": 'I',
            "LEU": 'L',
            "LYS": 'K',
            "MET": 'M',
            "PHE": 'F',
            "PRO": 'P',
            "SER": 'S',
            "THR": 'T',
            "TRP": 'W',
            "TYR": 'Y',
            "VAL": 'V',
        }
        return switcher.get(amino_acid, "0")
        
    def run_ialign(self,PDBCode):
        start = time.time()
        f = open(Path(self.TargetListDir+'/target.lst'), "w")
        PDBFileLoc=str(Path(self.PDBLib+'/Optimized_'+PDBCode+'_Repair.pdb'))
        f.write(PDBFileLoc+"\n")
        #print(str(Path(self.TargetListDir+'/target.lst'))+"\n")
        bashCommand = 'perl ' +self.Bin + ' -a 2 --mini 1 --minp 1 -w ' +self.OutDir + ' -p ' + self.Library + ' -l '  + self.LibraryList +' '+PDBFileLoc
        #bashCommand = "perl /home/andronikos/Documents/ialign/bin/ialign.pl -a 2 -mini 9 -w /home/andronikos/Documents/Optimus_Bind/ialign_py/test1_out -p /home/andronikos/Documents/Optimus_Bind/ialign_py/test1 1a4eABAB_int.pdb 1A4Y.pdb"
        outfile=self.OutDir+'/'+PDBCode+'.log'
        log = open(Path(outfile), 'w')
        process = subprocess.Popen(bashCommand.split(), stdout=log)
        output, error = process.communicate()
        end = time.time()
        f = open(Path(self.OutTime), "a")
        f.write(PDBCode+"\t"+str(end-start)+"\n")
        #print (bashCommand)
        print (PDBCode+"\t"+str(end-start)+"\n")
        return outfile
        
    def ialign_self(self,PDBCode):
        PDBFileLoc=str(Path(self.PDBLib+'/Optimized_'+PDBCode+'_Repair.pdb'))
        bashCommand = 'perl ' +self.Bin + ' -a 2 --mini 1 --minp 1 -do_self -w ' +self.OutPDBDir + ' -p ' +self.OutPDBDir + ' -p2 ' + self.PDBLib + ' -p1 '  + self.PDBLib+' '+PDBFileLoc
        outfile=self.OutPDBDir+'/'+PDBCode+'.log'
        log = open(Path(outfile), 'r')
        process = subprocess.Popen(bashCommand.split(), stdout=log)
        output, error = process.communicate()
        return outfile
            
    def get_interface(self,PDBCode):
        pattern=re.compile('\s+\d+\s+\w\s+\d+\s+\w{3}\s+\w\s+\d+\s+\w{3}\s+')
        # Loop through the file to find the lines using regex
        infile=r'C:\Users\jbren\OneDrive\Desktop\testing.log'
        log=open(infile,'r')
        interface={}
        lines=log.readlines()
        for line in lines:
            match=pattern.match(line)
            if match:
                line_list=line.split()
                # print(str(RowCt)+"\t"+str(+ChainCt)+"\n")
                if line_list[1] in interface:
                    interface[line_list[1]].append(line_list[2])
                else:
                    interface[line_list[1]]=[line_list[2]]
        print(interface)
        return interface
    
    def init_al(self,int_res):
        EmptyRow={}
        dict_temp = EmptyRow.fromkeys(int_res, '-')    # Use the positions as keys for the dictionary
        EmptyRow = {'pID': '-', 'is-score': 0,'rmsd': 0,'seq-id': 0}
        EmptyRow.update(dict_temp)
        return EmptyRow
                
    def read_al_file(self,PDBCode,ch_set):
        infile=r'C:\Users\jbren\OneDrive\Desktop\1A22.log'       
        interface=self.get_interface(PDBCode)
        #infile = open(Path(self.OutDir+'/'+PDBCode+'.log'), "r")       
        ial_df={}
        EmptyRow={}
        dict_list={}           
        n_align = 0         # Number of alignments (occured from running iAlign)
        for chain in interface:
            int_res=list(map(int,interface[chain]))
            EmptyRow[chain]=self.init_al(int_res)
            ial_df[chain]= pd.DataFrame(EmptyRow[chain],index=[0]) 
            dict_list[chain] = list()
        # The following lists will be used to store the lines that contain the information we want       
        is_lines = []       # Each item of the list contains the number of the line where the IS score is
        align_lines = []    # Each item of the list contains the number of the line where the residues' columns begin
        struct1_lines = []  # Each item of the list contains the number of the line where the name of the Protein interface is
        total_num_res = []        # List of the number of the residues
        rmsd_line_list = []
        # Regex patterns to search the output file of iAlign
        pattern1 = "\A>>>"
        pattern2 = "\AIS-score"
        pattern3 = "\A Index"
        pattern4 = "\ARMSD"
        # Loop through the file to find the lines using regex
        for i, line in enumerate(open(infile)):
            for match in re.finditer(pattern1, line):
                struct1_lines.append(i+1)
            for match in re.finditer(pattern2, line):
                is_lines.append(i+1)
            for match in re.finditer(pattern3, line):
                align_lines.append(i+2)
            for match in re.finditer(pattern4, line):
                rmsd_line_list.append(i+1)
                
        # The length of the two lists must be the same
        if (len(is_lines) == len(align_lines)):
            n_align = len(align_lines)
         
        # # loop over all proteins in the library (number=n_align)
        # # One dictionary storing the line for each side of the interface

        for i in range(n_align):      
            # Get the name of the target protein and use it as ID
            struct1_line = linecache.getline(infile, struct1_lines[i])
            struct_spl = struct1_line.split()
            stat_row={}
            pID = struct_spl[0]
            pID = pID[3:]
            # Get IS score
            is_score_line = linecache.getline(infile, is_lines[i])
            score_spl = is_score_line.split()
            isc = score_spl[2]
            if isc[-1] == ',':
                isc = isc[:-1]  # If the IS-Score ends with a comma remove it
            rmsd_line = linecache.getline(infile, rmsd_line_list[i])
            rmsd_spl = rmsd_line.split()
            rmsd = rmsd_spl[2]
            if rmsd[-1] == ',':
                rmsd = float(rmsd[:-1])  # If the IS-Score ends with a comma remove it
            seq_id = float(rmsd_spl[6])
            
            stat_row['is-score'] = float(isc)
            stat_row['pID'] = pID
            stat_row['rmsd'] = float(rmsd)
            stat_row['seq-id'] = seq_id
        
            # Get number of aligned residues for that protein (total_num_res)
            num_res_line = linecache.getline(infile, is_lines[i] + 1)
            num_spl = num_res_line.split()
            total_num_res.append(int(num_spl[5]))  # See how many lines you need to read
        
            # Get all aligned residues in list res1 and the positions in res_pos
            # This could probably be set off in a sub-function with the chain id as input
            res_pos = {}
            prot_row={}
            for chain in interface:
                res_pos[chain]={}
                prot_row[chain]={}

            for n in range(align_lines[i], align_lines[i] + total_num_res[i]):
                al_line = linecache.getline(infile, n).split()
                if al_line[4] in res_pos:
                    res_pos[al_line[4]][int(al_line[5])]=self.changeAA(al_line[3])
            for chain in interface:
                prot_row[chain].update(EmptyRow[chain])
                # Replace '-' with the aligned residue (if it exists) and update the dataframe
                prot_row[chain].update(res_pos[chain])
                prot_row[chain].update(stat_row)
                dict_list[chain].append(prot_row[chain])
        for chain in interface:    
            ial_df[chain] = ial_df[chain].append(dict_list[chain], ignore_index=True)
            ial_df[chain] = ial_df[chain].loc[ial_df[chain]['pID'] != '-']
        return (ial_df)
    
    def filter_al(self,ial_df,cutoff):
        filt_df = ial_df[ial_df['is-score'] > cutoff] 
        return filt_df
    
    def score_al(self,ial_df,pos,WTAA,mutAA):
        pseudo=0
        try:
            WT = ial_df[pos].value_counts()[WTAA]
        except:
            WT=0
        try:
            mut = ial_df[pos].value_counts()[mutAA]
        except:
            mut=0
        score=ln((WT+pseudo)/(mut+pseudo))
        return score

        # #print(df_1r0r)       
        
# #test
# iclass=ialign()
# iclass.name='1R0R'        
# iclass.run_ialign('1R0R')
