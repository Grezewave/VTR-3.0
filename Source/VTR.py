# VTR 3.4 for Linux
import os
import VTR_Geometric

while (1):
     os.system("clear")
     os.system("echo ____________________________________________________________\n")
     os.system("echo VTR Geometric v 3.4\n")
     os.system("echo ____________________________________________________________\n")
     os.system("echo Proteins in directory: ")
     os.system("find ../Data -maxdepth 1 -not -type d")
     os.system("echo  ___________________________________________________________\n")
     rtt_protein = input("Choose the first:")
     while (not(os.path.exists("../Data/" + rtt_protein + ".pdb"))):
          print("This file does not exist...")
          rtt_protein = input("Choose the first:")
     stc_protein = input("Choose the second:")
     while (not(os.path.exists("../Data/" + stc_protein + ".pdb"))):
          print("This file does not exist...")
          stc_protein = input("Choose the second:")
     cutoff = input("Choose the cutoff:")
     chainrtt1 = input("Choose the first chain for " + rtt_protein +"(A, B ..., press enter for run all x all chains):")
     if chainrtt1 == "":
          chainrtt1 = "/"
     chainrtt2 = input("Choose the second chain for " + rtt_protein +"(A, B ..., press enter for run all x all chains):")
     if chainrtt2 == "":
          chainrtt2 = "/"
     chainstc1 = input("Choose the first chain for " + stc_protein +"(A, B ..., press enter for run all x all chains):")
     if chainstc1 == "":
          chainstc1 = "/"
     chainstc2 = input("Choose the second chain for " + stc_protein +"(A, B ..., press enter for run all x all chains):")
     if chainstc2 == "":
          chainstc2 = "/"

     question = 0
     selection = '''Please Select the Execution Type:
     (1) Include main chain Atoms(N and O) in contact calculation
     (2) Do not include'''
     
     while (question > 2 or question < 1):
          print(selection)
          question = int(input('   '))
     if (question == 2):
          VTR_Geometric.Match("../Data/" + rtt_protein + ".pdb", "../Data/" + stc_protein + ".pdb",cutoff = cutoff,mtype = 'x',chain11 = chainrtt1,chain12 = chainrtt2,chain21 = chainstc1,chain22 = chainstc2)
     elif (question == 1):
          VTR_Geometric.Match("../Data/" + rtt_protein + ".pdb", "../Data/" + stc_protein + ".pdb",cutoff = cutoff,mtype = 'd',chain11 = chainrtt1,chain12 = chainrtt2,chain21 = chainstc1,chain22 = chainstc2)
     input('Press Enter to Return')
     

