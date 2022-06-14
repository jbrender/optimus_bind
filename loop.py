# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:57:13 2022

@author: jbren
"""
PDB_Loc="C:\\Users\\jbren\\OneDrive\\Documents\\Optimus\\SKEMPI2_PDBs\\PDBs"


from SKEMPI import skempi_final as db
#from runfoldx_old import foldx as foldx
from protein import ProteinMethods
from ialign import ialign as ialign
from Mutation import MutationMethods
import regex as re
import numpy as np

def get_pdbs(db):
# Returns an array of the protein codes
    pdbs= db['#Pdb'].str.slice(start=0, stop=4, step=1)
    # Returns only the unique values from the array
    uniq_pdbs=pdbs.unique()   
    return uniq_pdbs

# =============================================================================
# Retrieve size of protein in residues
# Used to route larger proteins to node with more memory
def get_size(pdb_code):
    p=ProteinMethods()    
    p.PDB_folder=PDB_Loc
    x=p.getPDBfile(pdb_code)
    xx=p.getStructure(x)
    max_res=0
    for chain in xx:
        for residue in chain.get_residues():
            tag=residue.get_full_id()[3]
            if (tag[0])==" ":
                max_res=max_res+1
    return max_res
 
def rm_large_prot(db):
    pdb_list=get_pdbs(db)
    pdb_size=[get_size(pdb) for pdb in pdb_list]
    pdb_size=dict(zip(pdb_list,pdb_size))
    dict_large={key: value for key, value in pdb_size.items() if value>300}
    dict_small={key: value for key, value in pdb_size.items() if value<300}
    pdbs= db['#Pdb'].str.slice(start=0, stop=4, step=1)
    db_large=db[pdbs.isin(list(dict_large.keys()))]
    db_small=db[pdbs.isin(list(dict_small.keys()))]
    return db_large,db_small           
# =============================================================================
# Divide database into equal size sets regardless of size
# Div=number of divisions
# Offset=set slice  
def slice_db_row(db,div,offset):  
    row_nums=[db.index.get_loc(i) for i in db.index]
    row_set = [i for i in row_nums if i%div==offset]
    db_slice=db.iloc[row_set]
    return db_slice
# =============================================================================
    
def slice_db_prot(pdb_list,div,offset):
    prot_slice=pdb_list[offset::div]
    return prot_slice
    

def SplitByProt(db):
    slice_pdbs=[]
    pdbs=get_pdbs(db)
    for i in range(4):
        slice_pdbs.append(slice_db_prot(pdbs,4,i))
    return slice_pdbs
# =========================================================================
# Run the largest proteins on the 4GB server and split the remainder into three parts
def SplitByProtSize(db):
    split_db_size=rm_large_prot(db)
    slice_pdbs=[]
    pdb_large=get_pdbs(split_db_size[0])
    pdb_small=get_pdbs(split_db_size[1])
    for i in range(3):
        slice_pdbs.append(slice_db_prot(pdb_small,3,i))
    return slice_pdbs

def RunIAlignLoop(pdb):
    ialignDict[pdb] = ialign()
    try:
        ialignDict[pdb].run_ialign(pdb)
    except:
        file = open("errors_ialign.txt","a+")  
        file.write(pdb+"\n")      
        print("Error"+pdb+"\n")
        
def ReadIAlignLoop(pdb,chain):
    ialignDict[pdb] = ialign()
    try:
        ialignDict[pdb].read_al_file(pdb,chain)
    except:
        file = open("errors_read_ialign.txt","a+")  
        file.write(pdb+"\t"+chain+"\n")      
        print("Error"+pdb+"\n")
        
        
        
# =============================================================================
# Main
# =============================================================================
node=1
RunIAlignFlag=False
ReadIAlignFlag=True
# =============================================================================
# Loop over proteins
# =============================================================================
PDBSlice=SplitByProt(db)
ialignDict= {}
for pdb in PDBSlice[node]:
    if RunIAlignFlag==True:
        RunIAlignLoop(pdb)
    ialignDict[pdb] = ialign()
    try:
        interface=ialignDict[pdb].interface_residues(pdb)
        print(interface)
    except:
        print(pdb)
# =============================================================================
# Loop over chains        
# =============================================================================
db_subset=db[db['PDB'].isin(PDBSlice[node])]
uniq_chains= db_subset.groupby(['PDB', 'Prot1Chain','Prot2Chain']).size()   
for items in uniq_chains.iteritems():
    ChainSetA=items[0][1]
    ChainSetB=items[0][2]
    pdb=items[0][0]
    if ReadIAlignFlag==True:
        for chA in ChainSetA:
            ReadIAlignLoop(pdb,chA)
        for chB in ChainSetB:
            ReadIAlignLoop(pdb,chB)

        
# =============================================================================
# Loop over mutations
# =============================================================================

m=MutationMethods() 
for i in db.index:   
    # 
    db.loc[i,'Blosum']=m.getSMScore([db['Mutation(s)_cleaned'].loc[i]],'Blosum')
    db.loc[i,'ITSM']=m.getSMScore([db['Mutation(s)_cleaned'].loc[i]],'ITSM')
    db.loc[i,'Proline']=m.CheckProline([db['Mutation(s)_cleaned'].loc[i]])
    db.loc[i,'Glycine']=m.CheckGlycine([db['Mutation(s)_cleaned'].loc[i]])
    # print(i)


    


#        #dummy placeholder command
#  #       print (str(db.index[i])+"_"+pdb.iloc[i]+"_"+db['Mutation(s)_cleaned'].iloc[i])

        

# slices=SplitByProtSize(db)
# for pdb in slices[node]:
#     icode='i'+pdb
#     icode=ialign()
   # icode.run_ialign(pdb)       
    # i1A22.run_ialign('1A22')


# for i in range(3):
#     slice=slice_db_row(split_db_size[1],3,i)
#     slice_pdbs.append(get_proteins(slice))