import OSfunct
import os

def colorscale(VMD, cutoff, out):
    #return a rgb code color between blue and red, based on cutoff and VMD, bigger the VMD, bluest the color
    Redest = [255,0,0]
    Bluest = [0,0,255]
    R = (((-255)/cutoff)*VMD)+255
    G = 0
    B = (((255)/cutoff)*VMD)
    if out == 'l':
        color = [int(R),int(G),int(B)]
    elif out == 't':
        color = (float(B)/255,float(G),float(R)/255)
    return(color)

def colorchange(index):
    #return the color of pymol type
    colors = ["aquamarine","blue","bluewhite","br0","br1","br2","br3","br4","br5","br6","br7","br8","br9","brightorange","brown","carbon","chartreuse","chocolate","cyan","darksalmon","dash","deepblue","deepolive","deeppurple","deepsalmon","deepteal","density","dirtyviolet","firebrick","forest","gray","green","greencyan","grey","hotpink","hydrogen","lightblue","lightmagenta","lightorange","lightpink","lightteal","lime","limegreen","limon","magenta","marine","nitrogen","olive","orange","oxygen","palecyan","palegreen","paleyellow","pink","purple","purpleblue","raspberry","red","ruby","salmon","sand","skyblue","slate","smudge","splitpea","sulfur","teal","tv_blue","tv_green","tv_orange","tv_red","tv_yellow","violet","violetpurple","warmpink","wheat","white","yellow","yelloworange"]
    return (colors[(index)%(len(colors)-1)])

def detailed_ploter(rtt_path, stc_path, matches, rtt_dismatches, stc_dismatches, cutoff, _type):
    #create a pml file for the full structure superposition with the matchs in highlights
    folder = OSfunct.create_dir("../Plots",rtt_path,stc_path,_type)
    #define ne file name
    pmlname = "../Plots/"+folder+"/c_scale" + rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + "_x_" + stc_path[stc_path.rfind("/")+1:-4] + ".pml"
    pml = open(pmlname,'w')
    pml.write("load ../../Data/" + rtt_path[rtt_path.rfind("/")+1:-11]+'x'+stc_path[stc_path.rfind("/")+1:-4] +'_align/'+rtt_path[rtt_path.rfind("/")+1:] + "\n")
    pml.write("load ../../Data/" + stc_path[stc_path.rfind("/")+1:] + "\n")
    x = 0
    #put each match in the visualization, with a color gradded equivalent to VMD
    for i in matches:
        color = colorscale(i.VMD(),cutoff,'l')
        entry = "set_color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", [" + str(color[0]) + "," + str(color[1]) + "," + str(color[2]) + "]\n"
        pml.write(entry)
        selection1 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom1.id)
        entry = "select " + selection1 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom1.id) + "\n"
        pml.write(entry)
        selection2 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom2.id)
        entry = "select " + selection2 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        selection1 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom1.id)
        entry = "select " + selection1 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom1.id) + "\n"
        pml.write(entry)
        selection2 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom2.id)
        entry = "select " + selection2 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        x += 1
    pml.write("set_color 125_125_125,[125,125,125]\n")
    entry = ' and('
    #Put the dismatches in a separated sellection,on gray color
    for i in rtt_dismatches:
        entry+=" id "+str(i.atom1.id)
        entry+=" or id "+str(i.atom2.id)+" or "
    pml.write('select '+rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")]+'_dismatches,model '+rtt_path[rtt_path.rfind("/")+1:-4]+entry[:-4]+')\n')
    pml.write("color 125_125_125, "+rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")]+"_dismatches\n")

    entry = ' and('
    for i in stc_dismatches:
        entry+=" id "+str(i.atom1.id)
        entry+=" or id "+str(i.atom2.id)+" or "
    pml.write('select '+stc_path[stc_path.rfind("/")+1:-4]+'_dismatches,model '+stc_path[stc_path.rfind("/")+1:-4]+entry[:-4]+')\n')
    pml.write("color 125_125_125, "+stc_path[stc_path.rfind("/")+1:-4]+"_dismatches\n")
    
    hideH = "sele resn HOH\nhide (sele)" 
    pml.write(hideH)
    pml.close()
    #return the folder name
    return(folder)

def multi_ploter(rtt_path, stc_path, matches, cutoff, folder):
    #create a visualization of each contact on Pymol
    for i in matches:
        #load the proteins, and hide everything
        pmlname = "../Plots/" + folder + "/" + i.rtt_contact.residue1.id + str(i.rtt_contact.residue1.parameter)+ i.rtt_contact.atom1.type+str(i.rtt_contact.atom1.id)+ "--" + i.rtt_contact.residue2.id + str(i.rtt_contact.residue2.parameter)+i.rtt_contact.atom2.type + str(i.rtt_contact.atom2.id)+ "_x_" + i.stc_contact.residue1.id + str(i.stc_contact.residue1.parameter)+i.stc_contact.atom1.type + str(i.stc_contact.atom1.id) + "--" + i.stc_contact.residue2.id + str(i.stc_contact.residue2.parameter)+i.stc_contact.atom2.type + str(i.stc_contact.atom2.id) + ".pml"
        pml = open(pmlname,'w')
        pml.write("load ../../Data/" + rtt_path[rtt_path.rfind("/")+1:-11]+'x'+stc_path[stc_path.rfind("/")+1:-4] +'_align/'+rtt_path[rtt_path.rfind("/")+1:] + "\n")
        pml.write("load ../../Data/" + stc_path[stc_path.rfind("/")+1:] + "\n")
        entry = "hide all\n"
        pml.write(entry)
        #show residues as sticks
        selection1 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + i.rtt_contact.residue1.id + str(i.rtt_contact.residue1.parameter)
        entry = "select " + selection1 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and chain " + i.rtt_contact.chain1.id + " and resi " + str(i.rtt_contact.residue1.parameter) + "\n"
        pml.write(entry)
        selection2 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + i.rtt_contact.residue2.id + str(i.rtt_contact.residue2.parameter)
        entry = "select " + selection2 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and chain " + i.rtt_contact.chain2.id + " and resi " + str(i.rtt_contact.residue2.parameter) + "\n"
        pml.write(entry)
        entry = "show sticks, " + selection1 + " " + selection2 + "\n"
        pml.write(entry)
        selection1 = stc_path[stc_path.rfind("/")+1:-4] + i.stc_contact.residue1.id + str(i.stc_contact.residue1.parameter)
        entry = "select " + selection1 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and chain " + i.stc_contact.chain1.id + " and resi " + str(i.stc_contact.residue1.parameter) + "\n"
        pml.write(entry)
        selection2 = stc_path[stc_path.rfind("/")+1:-4] + i.stc_contact.residue2.id + str(i.stc_contact.residue2.parameter)
        entry = "select " + selection2 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and chain " + i.stc_contact.chain2.id + " and resi " + str(i.stc_contact.residue2.parameter) + "\n"
        pml.write(entry)
        entry = "show sticks, " + selection1 + " " + selection2 + "\n"
        pml.write(entry)
        entry = "sele resn HOH\nhide (sele)\n"
        pml.write(entry)

        #set each match, criating a dashed line between atoms
        color = colorscale(i.VMD(),cutoff,'l')
        entry = "set_color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", [" + str(color[0]) + "," + str(color[1]) + "," + str(color[2]) + "]\n"
        pml.write(entry)
        if i.rtt_contact.centroid1:
            selection1 = 'rtt_centroid1'
            entry = 'pseudoatom ' + selection1 + ', pos=[' + str(i.rtt_contact.centroid1.x) + ', ' + str(i.rtt_contact.centroid1.y) + ', ' + str(i.rtt_contact.centroid1.z) +']\n'
        else:
            selection1 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom1.id)
            entry = "select " + selection1 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom1.id) + "\n"
        pml.write(entry)
        if i.rtt_contact.centroid2:
            selection2 = 'rtt_centroid2'
            entry = 'pseudoatom ' + selection2 + ', pos=[' + str(i.rtt_contact.centroid2.x) + ', ' + str(i.rtt_contact.centroid2.y) + ', ' + str(i.rtt_contact.centroid2.z) +']\n'
        else:
            selection2 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom2.id)
            entry = "select " + selection2 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        if i.stc_contact.centroid1:
            selection1 = 'stc_centroid1'
            entry = 'pseudoatom ' + selection1 + ', pos=[' + str(i.stc_contact.centroid1.x) + ', ' + str(i.stc_contact.centroid1.y) + ', ' + str(i.stc_contact.centroid1.z) +']\n'
        else:
            selection1 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom1.id)
            entry = "select " + selection1 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom1.id) + "\n"
        pml.write(entry)
        if i.stc_contact.centroid2:
            selection2 = 'stc_centroid2'
            entry = 'pseudoatom ' + selection2 + ', pos=[' + str(i.stc_contact.centroid2.x) + ', ' + str(i.stc_contact.centroid2.y) + ', ' + str(i.stc_contact.centroid2.z) +']\n'
        else:
            selection2 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom2.id)
            entry = "select " + selection2 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        pml.write('hide lines\n')
