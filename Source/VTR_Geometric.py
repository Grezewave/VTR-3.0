#V1.2
#Contacts match by geometric propertie, VMD(Vector Medium Distance)
#Some entry changes
#Write a plot .pml file
#Match skill improved and refined
#Most detailed analisys avaiable(by residue and color scale)
#Most Graphic detail and execution stats
#More run options(Chain filter)

import sys
import Classify
import Contacts
import OSfunct
import time
import VTR_Functions as vtr
import PymolGen as py
import os

class Match:
    #parameters---------------------------------------------------------------------
    def __init__(self, _protein1, _protein2, *args, **kwargs):
        self.rtt_protein = 0
        self.stc_protein = 0
        self.rtt_contacts = []
        self.stc_contacts = []
        self.matches = []
        self.rtt_dismatches = []
        self.stc_dismatches = []
        self.RMSD = False
        self.VTR = False
        self.mAVD = False
        self.number = 0
        self.time = 0
        self.type = ''
        self.chains = ''
        self.cutoff = 0
        start = time.time()
        _cutoff = kwargs.get('cutoff', 2)
        self.cutoff = _cutoff
        __type = kwargs.get('mtype', 'x')
        _chain11 = kwargs.get('chain11', '/')
        _chain12 = kwargs.get('chain12', '/')
        _chain21 = kwargs.get('chain21', '/')
        _chain22 = kwargs.get('chain22', '/')
        self.chains = _chain11 + '-' + _chain12 + ' ' + _chain21 + '-' + _chain22
        protein1 = _protein1
        protein2 = _protein2
        _type = ''

        if "d" == __type:
            _type += '-d'
        self.type = _type

        outname = "../Logs/"+protein1[protein1.rfind("/")+1:-4]+"x"+protein2[protein2.rfind("/")+1:-4]+_type+"Log.txt"
        out = open(outname, 'w')

    #build aligned protein and contacts list --------------------------------------------------------------------------------------------
        path = OSfunct.TMAlign(protein1,protein2)
        rtt_name = "../Data/" + path + protein1[protein1.rfind("/"):-4] + "_rotate.pdb"
        rtt_protein = Classify.classify(rtt_name)
        stc_protein = Classify.classify(protein2)
        self.rtt_protein = rtt_protein
        self.stc_protein = stc_protein
        rtt_contacts = Contacts.contacts(rtt_protein,protein1[protein1.rfind("/"):-4] + "_rotate",_chain11,_chain12)
        stc_contacts = Contacts.contacts(stc_protein,protein2[protein2.rfind("/"):-4],_chain21,_chain22)
        self.rtt_contacts = rtt_contacts
        self.stc_contacts = stc_contacts

    #match contacts ---------------------------------------------------------------------------------------------------------------------
        matches,rtt_dismatches,stc_dismatches = vtr.match_contacts(rtt_contacts,stc_contacts,int(_cutoff))
        self.matches = matches
        self.rtt_dismatches = rtt_dismatches
        self.stc_dismatches = stc_dismatches
        end = time.time()

    #write output -----------------------------------------------------------------------------------------------------------------------
        vtr.write_dismatch(protein1,protein2,rtt_dismatches,stc_dismatches,_type)
        result = str(len(matches)) + " matches found" + "\n"
        self.number = len(matches)
        out.write(result)
        out.write("Match execution time: " + str(round(end-start,0))+" seconds\n")
        self.time = end-start
        if (0 == len(matches)):
            out.write("\n")
            out.write("\n")
            out.write(protein1[protein1.rfind("/")+1:]+"\n")
            if "/" != _chain11:
                out.write(_chain11+"\n")
            else:
                out.write("All\n")
            if "/" != _chain12:
                out.write(_chain12+"\n")
            else:
                out.write("All\n")
            out.write(protein2[protein2.rfind("/")+1:]+"\n")
            if "/" != _chain21:
                out.write(_chain21+"\n")
            else:
                out.write("All\n")
            if "/" != _chain22:
                out.write(_chain22+"\n")
            else:
                out.write("All\n")
            out.write(str(_cutoff)+"\n")

            if "d" in _type:
                out.write("Detailed")
            else:
                out.write("Simple")
            out.close()
        out.write("Match RMSD = "+str(round(vtr.RMSD(matches, rtt_protein, stc_protein),2)) + "\n")
        self.RMSD = vtr.RMSD(matches, rtt_protein, stc_protein)
        VTR, mean_AVD = vtr.VTR(matches, rtt_contacts, stc_contacts ,len(rtt_dismatches),len(stc_dismatches),int(_cutoff))
        out.write("VTR = "+str(round(VTR,2)) + "\n")
        self.VTR = VTR
        out.write("Average AVD = "+str(round(mean_AVD,2)) + " A\n")
        self.mAVD = mean_AVD
        vtr.writer(protein1,protein2,rtt_protein,stc_protein,rtt_contacts,stc_contacts,matches,_type)

    #make plots -------------------------------------------------------------------------------------------------------------------------
        folder = py.detailed_ploter(rtt_name, protein2, matches,rtt_dismatches,stc_dismatches, int(_cutoff), _type)
        py.multi_ploter(rtt_name, protein2, matches, int(_cutoff), folder)

    #make graphs ------------------------------------------------------------------------------------------------------------------------
        folder = OSfunct.create_dir("../Graphs",rtt_name,protein2,_type)

        if 'd' in _type:
            vtr.freq_VMD(matches,int(_cutoff),"d",folder,protein1[protein1.rfind("/")+1:-4],protein2[protein2.rfind("/")+1:-4])
        else:
            vtr.freq_VMD(matches,int(_cutoff),"x",folder,protein1[protein1.rfind("/")+1:-4],protein2[protein2.rfind("/")+1:-4])

    #write recent parameters ------------------------------------------------------------------------------------------------------------
        out.write(protein1[protein1.rfind("/")+1:]+"\n")
        if "/" != _chain11:
            out.write(_chain11+"\n")
        else:
            out.write("All\n")
        if "/" != _chain12:
            out.write(_chain12+"\n")
        else:
            out.write("All\n")
        out.write(protein2[protein2.rfind("/")+1:]+"\n")
        if "/" != _chain21:
            out.write(_chain21+"\n")
        else:
            out.write("All\n")
        if "/" != _chain22:
            out.write(_chain22+"\n")
        else:
            out.write("All\n")
        out.write(str(_cutoff)+"\n")

        if "d" in _type:
            out.write("Detailed")
        else:
            out.write("Simple")
        out.close()

        output = open(outname, 'r')
        lines = output.readlines()
        for data in range(0,5):
            print(lines[data])
    def to_csv(self,_list, outname):
        out = open(outname, 'w')
        stack = ''''''
        for i in _list:
            if type(i) is Contacts.Contact:
                stack+=i.chain1.id + ',' + i.residue1.id + str(i.residue1.parameter) + ',' + i.atom1.type.strip() + ','
                stack+=i.chain2.id + ',' + i.residue2.id + str(i.residue2.parameter) + ',' + i.atom2.type.strip() + ',' 
                stack+=str(round(i.distance,2)) + ','
                for x in i.contacts:
                    stack += x + ' '
                stack += '''
'''
            elif type(i) is vtr.match:
                stack+=i.rtt_contact.chain1.id + ',' + i.rtt_contact.residue1.id + str(i.rtt_contact.residue1.parameter) + ',' + i.rtt_contact.atom1.type.strip() + ','
                stack+=i.rtt_contact.chain2.id + ',' + i.rtt_contact.residue2.id + str(i.rtt_contact.residue2.parameter) + ',' + i.rtt_contact.atom2.type.strip() + ',' 
                stack+=str(round(i.rtt_contact.distance,2)) + ','
                for x in i.rtt_contact.contacts:
                    stack += x + ' '
                stack += ','
                stack+=i.stc_contact.chain1.id + ',' + i.stc_contact.residue1.id + str(i.stc_contact.residue1.parameter) + ',' + i.stc_contact.atom1.type.strip() + ','
                stack+=i.stc_contact.chain2.id + ',' + i.stc_contact.residue2.id + str(i.stc_contact.residue2.parameter) + ',' + i.stc_contact.atom2.type.strip() + ',' 
                stack+=str(round(i.stc_contact.distance,2)) + ','
                for x in i.stc_contact.contacts:
                    stack += x + ' '
                stack += ','
                stack += str(round(i.VMD(),2)) + '''
'''
            else:
                print('Incorrect input(Must be a list of Contacts or vtr.matches)')
                return False
        out.write(stack)