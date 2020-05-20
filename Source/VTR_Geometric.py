#V3.0
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

def main(_protein1, _protein2, *args, **kwargs):
#parameters -------------------------------------------------------------------------------------------------------------------------
    start = time.time()
    _cutoff = kwargs.get('cutoff', 2)
    __type = kwargs.get('mtype', 'x')
    _chain11 = kwargs.get('chain11', '/')
    _chain12 = kwargs.get('chain12', '/')
    _chain21 = kwargs.get('chain21', '/')
    _chain22 = kwargs.get('chain22', '/')
    protein1 = _protein1
    protein2 = _protein2
    _type = ''

    if "d" == __type:
        _type += '-d'

    outname = "../Logs/"+protein1[protein1.rfind("/")+1:-4]+"x"+protein2[protein2.rfind("/")+1:-4]+_type+"Log.txt"
    out = open(outname, 'w')

#build aligned protein and contacts list --------------------------------------------------------------------------------------------
    path = OSfunct.TMAlign(protein1,protein2)
    rtt_name = "../Data/" + path + protein1[protein1.rfind("/"):-4] + "_rotate.pdb"
    rtt_protein = Classify.classify(rtt_name)
    stc_protein = Classify.classify(protein2)
    rtt_contacts = Contacts.contacts(rtt_protein,protein1[protein1.rfind("/"):-4] + "_rotate",_chain11,_chain12)
    stc_contacts = Contacts.contacts(stc_protein,protein2[protein2.rfind("/"):-4],_chain21,_chain22)

#match contacts ---------------------------------------------------------------------------------------------------------------------
    matches,rtt_dismatches,stc_dismatches = vtr.match_contacts(rtt_contacts,stc_contacts,int(_cutoff))
    end = time.time()

#write output -----------------------------------------------------------------------------------------------------------------------
    vtr.write_dismatch(protein1,protein2,rtt_dismatches,stc_dismatches,_type)
    result = str(len(matches)) + " matches found" + "\n"
    out.write(result)
    out.write("Match execution time: " + str(round(end-start,0))+" seconds\n")
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
        return (0,rtt_dismatches,stc_contacts)
    out.write("RMSD = "+str(round(vtr.RMSD(matches, rtt_protein, stc_protein),2)) + "\n")
    VTR, mean_AVD = vtr.VTR(matches, rtt_contacts, stc_contacts ,len(rtt_dismatches),len(stc_dismatches),int(_cutoff))
    out.write("VTR = "+str(round(VTR,2)) + "\n")
    out.write("Average AVD = "+str(round(mean_AVD,2)) + " A\n")
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
    return (matches, rtt_dismatches, stc_dismatches)