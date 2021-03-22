# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:33:03 2021

@author: jbren
"""
# =============================================================================

# =============================================================================
from Bio.PDB import Selection
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast import NCBIXML
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
import re
        
class ProteinMethods:
    def __init__(self,PDB_folder_loc):
        self.PDB_folder=PDB_folder_loc
       
    def getPDBfile(self,pdbcode):
        PDBloc=self.PDB_folder+"\\"+pdbcode+".pdb"
        return PDBloc
    
    # Accepts a PDB file and returns a Biopython structure object    
    def getStructure(self,PDBfile):
        parser = PDBParser(PERMISSIVE=1)
        structure=parser.get_structure('x',PDBfile)
        return structure
        
    # Accepts a Biopython structure object and returns the sequence of the chain specified      
    def getSeq(self, PDBfile,chain):
        #for loop isn't really necessary here - should be only one record
        for record in SeqIO.parse(PDBfile, "pdb-atom"):
            if record.annotations['chain']==chain:
                seq=record._seq._data
        self.seq=seq
        return seq
    
    def writeSeq(self, seq,pdbcode,chain,fastaloc):
        fasta_file=fastaloc+"\\"+pdbcode+"_"+chain+".fasta"
        f = open(fasta_file, "w")
        f.write(">"+pdbcode+"\n")
        f.write(seq)
        f.close()
        return fasta_file
    
    def runPSIBlast(self,fasta_file,db_loc):
        cmdstr= 'psiblast -query '+fasta_file+' -db '+db_loc+' -evalue 0.0001 -outfmt 5 -out '+fasta_file[:-6]+'.xmld'
        #add check for completed process and error handling
        subprocess.run(cmdstr,shell=True)
        
    def readBlast(self,blast_file):
        MSA={}
        blastlist=[]
        with open(blast_file) as bf:
            blast_records = NCBIXML.parse(bf)
            for blast_record in blast_records:
                self.blast_record=blast_record
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        QueryGaps=[m.start() for m in re.finditer('-', hsp.query)]                        
                        AlignedStr=[letter for i, letter in enumerate(hsp.sbjct) if i not in QueryGaps]
                        AlignedStr=''.join(AlignedStr)
                        alignStart='-'*(hsp.query_start-1)
                        if (len(AlignedStr)<blast_record.query_length):
                            alignEnd='-'*(blast_record.query_length-hsp.query_end)
                        else:
                            alignEnd=''
                        MSA[alignment.accession]=alignStart+AlignedStr+alignEnd
                        blastentry = SeqRecord(Seq(alignStart+AlignedStr+alignEnd), id=alignment.accession)
                        blastlist.append(blastentry)
        self.MSA=MultipleSeqAlignment(blastlist)
        self.MSAInfo=AlignInfo.SummaryInfo(self.MSA)
        self.PSSM=self.MSAInfo.pos_specific_score_matrix(chars_to_ignore=['-'])
      
     
           
                
# ============================================================================+=
# test        
# =============================================================================
p=ProteinMethods("C:\\Users\\jbren\\OneDrive\\Documents\\Optimus\\SKEMPI2_PDBs\\PDBs")

x=p.readBlast('C://Users//jbren//OneDrive//Documents//Optimus//seq//1A22_A.xmld')
PDBfile=p.getPDBfile('1A22')
structure=p.getStructure(PDBfile)
seq=p.getSeq(PDBfile,'B')
# p.writeSeq(seq,'1A22','C:\\Users\\jbren\\OneDrive\\Documents\\Optimus')