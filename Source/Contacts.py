#v1.6
#This version analyzes contacts between diferent chains too, if there are.
#The progress can be verified by messages along the code(remove #)
#This version can be runned on cmd, just remove the '#' from the last line.
#Verify messages updated
#Bugs fixed
#Writing structure improved
#Better computational and biological eficience
#Better file read
#The input is a protein and a out name, not a file name
#Better class structure

#Contacts v1.6
import sys
import Classify

#This program take a protein object and find contacts
class Contact:
    def __init__(self,number,model,chain1,chain2,residue1,residue2,atom1,atom2,contacts, distance):
        self.number = number
        self.model = model
        self.chain1 = chain1
        self.chain2 = chain2
        self.residue1 = residue1
        self.residue2 = residue2
        self.atom1 = atom1
        self.atom2 = atom2
        self.contacts = contacts
        self.distance = distance
        self.centroid1 = False
        self.centroid2 = False

class State:
    def __init__(self, atom, residue, state):
        self.atom = atom
        self.residue = residue
        self.state = state

class Centroid:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

#Atom interactions type definitions    
hydrophobic = {"ALA" : ["CB  "],
               "ARG" : ["CB  ","CG  ","CD  "],
               "ASN" : ["CB  "],
               "ASP" : ["CB  "],
               "CYS" : ["CB  "],
               "GLN" : ["CB  ","CG  "],
               "GLU" : ["CB  ","CG  "],
               "HIS" : ["CB  ","CG  ","CD2 ","CE1 "],
               "ILE" : ["CB  ","CG1 ","CG2 ","CD1 "],
               "LEU" : ["CB  ","CG  ","CD1 ","CD2 "],
               "LYS" : ["CB  ","CG  ","CD  "],
               "MET" : ["CB  ","CG  ","CE  "],
               "PHE" : ["CB  ","CG  ","CD1 ","CD2 ","CE1 ","CE2 ","CZ  "],
               "PRO" : ["CB  ","CG  ","CD  "],
               "THR" : ["CG2 "],
               "TRP" : ["CB  ","CG  ","CD1 ","CD2 ","CE2 ","CE3 ","CH2 ","CZ  ","CZ2 ","CZ3 "],
               "TYR" : ["CB  ","CG  ","CD1 ","CD2 ","CE1 ","CE2 ","CZ  "],
               "VAL" : ["CB  ","CG1 ","CG2 "]}

positives = {"ARG" : ["NH1 ","NH2 "],
             "HIS" : ["ND1 ","NE2 "],
             "LYS" : ["NZ  "]}

negatives = {"ASP" : ["OD1 ","OD2 "],
             "GLU" : ["OE1 ","OE2 "]}

acceptors = {"ALA" : ["O   "],
             "ARG" : ["O   "],
             "ASN" : ["O   ","OD1 "],
             "ASP" : ["O   ","OD1 ","OD2 "],
             "CYS" : ["O   "],
             "GLN" : ["O   ","OE1 "],
             "GLU" : ["O   ","OE1 ","OE2 "],
             "GLY" : ["O   "],
             "HIS" : ["O   "],
             "ILE" : ["O   "],
             "LEU" : ["O   "],
             "LYS" : ["O   "],
             "MET" : ["O   "],
             "PHE" : ["O   "],
             "PRO" : ["O   "],
             "SER" : ["O   "],
             "THR" : ["O   "],
             "TRP" : ["O   "],
             "TYR" : ["O   "],
             "VAL" : ["O   "]}

donners = {"ALA" : ["N   "],
           "ARG" : ["N   ","NE  ","NH1 ","NH2 "],
           "ASN" : ["N   ","ND2 ","OD1 "],
           "ASP" : ["N   "],
           "CYS" : ["N   "],
           "GLN" : ["N   ","NE2 "],
           "GLU" : ["N   "],
           "GLY" : ["N   "],
           "HIS" : ["N   ","ND1 ","NE2 "],
           "ILE" : ["N   "],
           "LEU" : ["N   "],
           "LYS" : ["N   ","NZ  "],
           "MET" : ["N   "],
           "PHE" : ["N   "],
           "PRO" : ["N   "],
           "SER" : ["N   ","OG  "],
           "THR" : ["N   ","OG1 "],
           "TRP" : ["N   ","NE1 "],
           "TYR" : ["N   ","OH  "],
           "VAL" : ["N   "]}

aromathics = {"HIS" : ["CG  ","ND1 ","CD2 ","CE1 ","NE2 "],
              "PHE" : ["CG  ","CD1 ","CD2 ","CE1 ","CE2 ","CZ  "],
              "TRP" : ["CG  ","CD1 ","CD2 ","NE1 ","CE2 ","CE3 ","CZ2 ","CZ3 ","CH2 "],
              "TYR" : ["CD1 ","CD2 ","CE1 ","CE2 ","CG  ","CZ  "]}

sulfur =  {"CYS" : ["S   "]}

classes = [hydrophobic, positives, negatives, acceptors, donners, aromathics, sulfur]

def adistance(A,B):
    #Return the euclidean distance between two atoms
    D = (((A.x - B.x)**2)+((A.y - B.y)**2)+((A.z - B.z)**2))**(1/2)
    return D

def centroid_distance(A,B,resA,resB):
    #search for the aromatic ring and calculate distance between rings centroids to avoid redundance and miss matches in future
    ringMembersA = []
    ringMembersB = []
    if resA.id == 'PHE':
        for atom in resA.atoms:
            if atom.type == 'CZ  ':
                ringMembersA.append(atom)
                a = atom
            elif atom.type == 'CG  ':
                ringMembersA.append(atom)
                b = atom
            elif atom.type == 'CE1 ':
                ringMembersA.append(atom)
            elif atom.type == 'CE2 ':
                ringMembersA.append(atom)
            elif atom.type == 'CD1 ':
                ringMembersA.append(atom)
            elif atom.type == 'CD2 ':
                ringMembersA.append(atom)
        centroid1 = Centroid((a.x + b.x)/2,(a.y + b.y)/2,(a.z + b.z)/2)
    
    elif resA.id == 'TYR':
        for atom in resA.atoms:
            if atom.type == 'CZ  ':
                ringMembersA.append(atom)
                a = atom
            elif atom.type == 'CG  ':
                ringMembersA.append(atom)
                b = atom
            elif atom.type == 'CE1 ':
                ringMembersA.append(atom)
            elif atom.type == 'CE2 ':
                ringMembersA.append(atom)
            elif atom.type == 'CD1 ':
                ringMembersA.append(atom)
            elif atom.type == 'CD2 ':
                ringMembersA.append(atom)
        centroid1 = Centroid((a.x + b.x)/2,(a.y + b.y)/2,(a.z + b.z)/2)

    elif resA.id == 'HIS':
        for atom in resA.atoms:
            if atom.type == 'ND1 ':
                ringMembersA.append(atom)
                a = atom
            elif atom.type == 'CG  ':
                ringMembersA.append(atom)
                b = atom
            elif atom.type == 'CD2 ':
                ringMembersA.append(atom)
                c = atom
            elif atom.type == 'CE1 ':
                ringMembersA.append(atom)
                d = atom
            elif atom.type == 'NE2 ':
                ringMembersA.append(atom)
                e = atom
        centroid1 = Centroid((a.x + b.x + c.x + d.x + e.x)/5,(a.y + b.y + c.y + d.y + e.y)/5,(a.z + b.z + c.z + d.z + e.z)/5)

    elif resA.id == 'TRP':
        for atom in resA.atoms:
            if atom.type == 'CE2 ':
                ringMembersA.append(atom)
                a = atom
            elif atom.type == 'CD2 ':
                ringMembersA.append(atom)
                b = atom
            elif atom.type == 'CD1 ':
                ringMembersA.append(atom)
            elif atom.type == 'NE1 ':
                ringMembersA.append(atom)
            elif atom.type == 'CG  ':
                ringMembersA.append(atom)
            elif atom.type == 'CZ2 ':
                ringMembersA.append(atom)
            elif atom.type == 'CE3 ':
                ringMembersA.append(atom)
            elif atom.type == 'CH2 ':
                ringMembersA.append(atom)
            elif atom.type == 'CZ3 ':
                ringMembersA.append(atom)
        centroid1 = Centroid((a.x + b.x)/2,(a.y + b.y)/2,(a.z + b.z)/2)

#-------------------------------------------------------------------------------------------------------------------------------
    if resB.id == 'PHE':
        for atom in resB.atoms:
            if atom.type == 'CZ  ':
                ringMembersB.append(atom)
                a = atom
            elif atom.type == 'CG  ':
                ringMembersB.append(atom)
                b = atom
            elif atom.type == 'CE1 ':
                ringMembersB.append(atom)
            elif atom.type == 'CE2 ':
                ringMembersB.append(atom)
            elif atom.type == 'CD1 ':
                ringMembersB.append(atom)
            elif atom.type == 'CD2 ':
                ringMembersB.append(atom)
        centroid2 = Centroid((a.x + b.x)/2,(a.y + b.y)/2,(a.z + b.z)/2)
    
    elif resB.id == 'TYR':
        for atom in resB.atoms:
            if atom.type == 'CZ  ':
                ringMembersB.append(atom)
                a = atom
            elif atom.type == 'CG  ':
                ringMembersB.append(atom)
                b = atom
            elif atom.type == 'CE1 ':
                ringMembersB.append(atom)
            elif atom.type == 'CE2 ':
                ringMembersB.append(atom)
            elif atom.type == 'CD1 ':
                ringMembersB.append(atom)
            elif atom.type == 'CD2 ':
                ringMembersB.append(atom)
        centroid2 = Centroid((a.x + b.x)/2,(a.y + b.y)/2,(a.z + b.z)/2)

    elif resB.id == 'HIS':
        for atom in resB.atoms:
            if atom.type == 'ND1 ':
                ringMembersB.append(atom)
                a = atom
            elif atom.type == 'CG  ':
                ringMembersB.append(atom)
                b = atom
            elif atom.type == 'CD2 ':
                ringMembersB.append(atom)
                c = atom
            elif atom.type == 'CE1 ':
                ringMembersB.append(atom)
                d = atom
            elif atom.type == 'NE2 ':
                ringMembersB.append(atom)
                e = atom
        centroid2 = Centroid((a.x + b.x + c.x + d.x + e.x)/5,(a.y + b.y + c.y + d.y + e.y)/5,(a.z + b.z + c.z + d.z + e.z)/5)

    elif resB.id == 'TRP':
        for atom in resB.atoms:
            if atom.type == 'CE2 ':
                ringMembersB.append(atom)
                a = atom
            elif atom.type == 'CD2 ':
                ringMembersB.append(atom)
                b = atom
            elif atom.type == 'CD1 ':
                ringMembersB.append(atom)
            elif atom.type == 'NE1 ':
                ringMembersB.append(atom)
            elif atom.type == 'CG  ':
                ringMembersB.append(atom)
            elif atom.type == 'CZ2 ':
                ringMembersB.append(atom)
            elif atom.type == 'CE3 ':
                ringMembersB.append(atom)
            elif atom.type == 'CH2 ':
                ringMembersB.append(atom)
            elif atom.type == 'CZ3 ':
                ringMembersB.append(atom)
        centroid2 = Centroid((a.x + b.x)/2,(a.y + b.y)/2,(a.z + b.z)/2)

    D = adistance(centroid1, centroid2)    
    return D, ringMembersA, ringMembersB, centroid1, centroid2

def atomclass(residue, atom):
    stats = []
    #Return a list of interactions types of a atom
    for i in range(0,len(classes)):
        if residue in classes[i]:
            if atom in classes[i][residue]:
                stats.append(i)
    state = State(atom, residue, stats)
    return state

def search_stacking(atom1,stats1,atom2,stats2, res1, res2):
    ringMembersA = []
    ringMembersB = []
    contact = []
    distance = 0
    centroid1 = 0
    centroid2 = 0
    if (5 in stats1) and (5 in stats2):
        distance, ringMembersA, ringMembersB, centroid1, centroid2 = centroid_distance(atom1,atom2, res1, res2)
        if (distance >= 3) and (distance <= 6):
            contact.append("Aromatic Stacking")
    return contact, ringMembersA, ringMembersB, distance, centroid1, centroid2

def defcontact(atom1,stats1,atom2,stats2):
    #Define contacts between two atoms based on interactions and distance
    distance = adistance(atom1,atom2)
    contact = []
    if distance >=2:
        if (distance <= 3.8):
            if (0 in stats1) and (0 in stats2):
                contact.append("Hydrophobic")
                if (5 in stats1) and (5 in stats2):
                    contact.append("Aromatic Stacking")
        if (distance <= 6):
            if (1 in stats1) and (2 in stats2):
                contact.append("Attractive")
            elif (2 in stats1) and (1 in stats2):
                contact.append("Attractive")
            if (1 in stats1) and (1 in stats2):
                contact.append("Repulsive")
            elif (2 in stats1) and (2 in stats2):
                contact.append("Repulsive")
        if (distance <= 3.2):
            if (3 in stats1) and (4 in stats2):
                contact.append("Hydrogen Bonds")
            elif (4 in stats1) and (3 in stats2):
                contact.append("Hydrogen Bonds")
    if (distance >= 1.5) and (distance <= 2.8):
        if (6 in stats1) and (6 in stats2):
            contact.append("Disulphide Bonds")
    return contact, distance
        


def contacts(protein,outname,chain1,chain2):
    #define contacts
    contacts = []
    outfile = "../Results/Contacts/" + outname + "_Contacts.txt"
    out = open(outfile,'w')
    x = 0
    allow = 0
    #filter the chains and try to make intrachain contacts
    for a in protein.chains:
        if chain1 == '/':
            if chain2 == '/':
                allow = 1
            elif chain2 == a.id:
                allow = 1
                
        elif chain1 == a.id:
            if chain2 == '/':
                allow = 1
            elif chain2 == a.id:
                allow = 1

        #if the chains match, try to make contact in the chain
        if allow:
            for b in range(0, len(a.residues)):
                ringMembersA = []
                for c in range(b+2, len(a.residues)):
                    ringMembersB = []
                    for d in a.residues[b].atoms:
                        for e in a.residues[c].atoms:
                            #define the atom interaction types
                            stat1 = atomclass(a.residues[b].id,d.type)
                            stat2 = atomclass(a.residues[c].id,e.type)
                            #verify interactions of N and O atoms
                            
                            if ("O   " == stat1.atom) or ("N   " == stat1.atom):
                                if ("O   " == stat2.atom) or ("N   " == stat2.atom):
                                    if abs(a.residues[b].parameter - a.residues[c].parameter) <= 4:
                                        validation = 0
                                    else: 
                                        validation = 1
                                else:
                                    validation = 1
                            else:
                                validation = 1
                            if validation:
                                #try to define the contacts
                                if (d not in ringMembersA) and (e not in ringMembersB):
                                    stack, ringMembersA, ringMembersB, cdistance, centroid1, centroid2 = search_stacking(d,stat1.state,e,stat2.state, a.residues[b], a.residues[c])
                                else:
                                    stack = []
                                atomcontact, distance = defcontact(d,stat1.state,e,stat2.state)
                                if (atomcontact != []):
                                    #create a contact object and put in the list
                                    item = Contact(x,outname,a,a,a.residues[b],a.residues[c],d,e,atomcontact,distance)
                                    contacts.append(item)
                                    x+=1
                                if (stack != []):
                                    item = Contact(x,outname,a,a,a.residues[b],a.residues[c],d,e,stack,cdistance)
                                    item.centroid1 = centroid1
                                    item.centroid2 = centroid2
                                    contacts.append(item)
                                    x+=1                             


        allow = 0

    #filter the chains and try to make interchain contacts
    for f in range(0,len(protein.chains)):
        for g in range(f+1,len(protein.chains)):

            if chain1 == '/':
                if chain2 == '/':
                    allow = 1
                elif (chain2 == protein.chains[f].id) or (chain2 == protein.chains[g].id):
                    allow = 1
                
            elif chain1 == protein.chains[f].id:
                if chain2 == '/':
                    allow = 1
                elif chain2 == protein.chains[g].id:
                    allow = 1

            elif chain1 == protein.chains[g].id:
                if chain2 == '/':
                    allow = 1
                elif chain2 == protein.chains[f].id:
                    allow = 1

            #if the chains match, try to make contact in the chain
            if allow:
                for h in protein.chains[f].residues:
                    ringMembersA = []
                    for i in protein.chains[g].residues:
                        ringMembersB = []
                        for j in h.atoms:
                            for k in i.atoms:   
                                #define the atom interaction types
                                stat1 = atomclass(h.id,j.type)
                                stat2 = atomclass(i.id,k.type)
                                #try to define the contacts
                                if (j not in ringMembersA) and (k not in ringMembersB):
                                    stack, ringMembersA, ringMembersB, cdistance, centroid1, centroid2 = search_stacking(j,stat1.state,k,stat2.state, h, i)
                                else:
                                    stack = []
                                atomcontact, distance = defcontact(j,stat1.state,k,stat2.state)
                                if (atomcontact != []):
                                    #create a contact object and put in the list
                                    item = Contact(x,outname,protein.chains[f],protein.chains[g],h,i,j,k,atomcontact, distance)
                                    contacts.append(item)
                                    x+=1
                                if (stack != []):
                                    ringStack = []
                                    ringStack.append(stack[0] + '(')
                                    for atom in ringMembersA:
                                        ringStack[0] += atom.type + str(atom.id) + ','
                                    ringStack[0] += ' - '
                                    for atom in ringMembersB:
                                        ringStack[0] += atom.type + str(atom.id) + ','
                                    ringStack[0] += ')'
                                    item = Contact(x,outname,protein.chains[f],protein.chains[g],h,i,j,k,ringStack, cdistance)
                                    item.centroid1 = centroid1
                                    item.centroid2 = centroid2
                                    contacts.append(item)
                                    x+=1
            allow = 0

    #Write output
    out.write(protein.header)
    out.write("\n")
    out.write(protein.title.rstrip(" "))
    out.write(":")
    out.write("\n\n")
    for i in contacts:
        out.write(str(i.number))
        out.write(" - ")
        out.write(i.chain1.id)
        out.write(" ")
        out.write(i.residue1.id)
        out.write(" ")
        out.write(str(i.residue1.parameter))
        out.write(" ")
        out.write(i.atom1.type)
        out.write(" VS ")
        out.write(i.chain2.id)
        out.write(" ")
        out.write(i.residue2.id)
        out.write(" ")
        out.write(str(i.residue2.parameter))
        out.write(" ")
        out.write(i.atom2.type)
        out.write("\n")
        out.write("Distance: ")
        out.write("%f"%i.distance)
        out.write(" - ")
        for e in i.contacts:
            out.write(e)
            out.write(", ")
        out.write("\n")
    print("File ",outfile[outfile.rfind("/")+1:]," generated")
    out.close()
    return contacts
