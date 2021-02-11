import os,math,csv
MAOA_Results = 'C:\\People\\Kaichen\\WATGEN\\Runs_4NM_fadFIXED\\Runs_4\\L1-357-001'
Lig = MAOA_Results + '\\Storage\\Lig'
Prot = MAOA_Results + '\\Storage\\Prot'
Prot_Lig = MAOA_Results + '\\Storage\\Prot_Lig'
Compare2 = MAOA_Results + '\\WaterAnalysisCompare2'


ThreeToOneLetter = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
                    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "LIG": "L"}
polar_C = { "ARG": "CZ", "ASN": "CG", "ASP": "CG", "GLU": "CD",  "GLN": "CD"}
Lig_polar = ['C13', 'C20','C25','C26','C28','C29','C30','C31','C32','C46','C47','C48']

#ligand_term, protein_term = 'Lig Final Waters.pdb','Prot Final Waters.pdb'

def aliphatic_C(carbon):
    if carbon[2][0] == 'C' and len(carbon[2]) > 1:
        return True
    else:
        return False 

def water_O(oxygen):
    if oxygen[2][0] == 'O' and oxygen[3] == 'WAT':
        return True
    else:
        return False

def HP_Interaction(ATM1,ATM2):
    coo1 = ATM1[-3:]
    coo2 = ATM2[-3:]
    distance = math.sqrt((coo1[0]-coo2[0])**2 + (coo1[1]-coo2[1])**2 + (coo1[2]-coo2[2])**2)
    if distance >= 3.5:
        return False
   
    else:
        return str(round(distance, 3))


def type_conv(split):
    new_line = split[0:6]
    new_line.append(float(split[6]))
    new_line.append(float(split[7]))
    new_line.append(float(split[8]))
    return new_line

#We need to do this three times:
#First, get ligand hydrophobic interaction(just ligand by itself)
#Second, get protein only hydrophobic interaction(with FAD inside pocket)
#Finally, get prot+lig hydrophobic interaction
#def read_LIG_MAOA(ligand_term,complex_term,folder):
    #We need to consider HP interactions between water-protein and water-ligand
def read_LIG_MAOA(ligand_term,folder1, complex_term, folder2, option):
    ##Option = 1, regular, option =2, read HP clashes between protein and ligand
    #If option =2, we should compare prot and proglig files in storage 
    pdb_list = []                
    alip_C, wat_O, lig_C = [], [], []
        
    os.chdir(folder1)
    for file in os.listdir(folder1):
        if file.endswith(ligand_term):
            pdb_list.append(file)
            
            
    
    if option == 1:
        os.chdir(folder2)
        for file in os.listdir(folder2):
            if file.endswith(complex_term) and file.startswith('new'):
                pdb_list.append(file)
        print("Two files to be read:", pdb_list)                
        watg1 = open(folder1+'\\'+pdb_list[0],'r+')
        
        lines0 = watg1.readlines()
        for line in lines0: 
            split = line.split()
            if split[3] != 'WAT':
                if  aliphatic_C(split) is True:
                    if split[3] in polar_C:
                        if split[2] != polar_C[split[3]]:                    
                            alip_C.append(type_conv(split))
                        #========================
                        
                    else:
                        if split[2] in Lig_polar:
                            pass
                        else:
                            alip_C.append(type_conv(split))
                    
        watg1.close()           
        watg4 = open(folder2+'\\'+pdb_list[1],'r+')
        lines1 = watg4.readlines()
        for line in lines1:
            split = line.split()
            if water_O(split) is True:
                wat_O.append(type_conv(split))            
           
        watg4.close()    
        return[alip_C, wat_O]
        
    elif option  == 2:
        FAD = []
        os.chdir(folder2)
        for file in os.listdir(folder2):
            if file.endswith(complex_term):
                pdb_list.append(file)
        print("Two files to be read:", pdb_list)                
        watg1 = open(folder1+'\\'+pdb_list[0],'r+')
        
        lines0 = watg1.readlines()
        for line in lines0: 
            split = line.split()
            if split[3] != 'WAT' :
                if  aliphatic_C(split) is True:
                    if split[3] != 'LIG': 
                        
                        if split[3] in polar_C:
                            if split[2] != polar_C[split[3]]:                
                                alip_C.append(type_conv(split))
                        else:
                            
                            alip_C.append(type_conv(split))
                        #========================
                        
                    else:
                        FAD.append(split[2])
#                        if split[2] in Lig_polar:
#                            pass
#                        else:
#                            alip_C.append(type_conv(split))
                    
        watg1.close() 
        print(FAD)          
        watg4 = open(folder2+'\\'+pdb_list[1],'r+')
        lines1 = watg4.readlines()        
        for k in range(len(lines1) - 500,len(lines1)):
            split = lines1[k].split()
            if split[3] == 'LIG':                
                if  aliphatic_C(split) is True:
                    
                    if split[2] not in FAD:
                        if split[2] not in Lig_polar:
                            lig_C.append(type_conv(split))
                    
#            if water_O(split) is True:
#                wat_O.append(type_conv(split))            
           
        watg4.close()

        return [alip_C, lig_C]
    else:
        print("Not an acceptable option. End reading PDB.")
        return 0





def selection_sort_alt(input_list):

    for idx in range(len(input_list)):

        min_idx = idx
        for j in range( idx +1, len(input_list)):
            if int(input_list[min_idx][0]) > int(input_list[j][0]):
                min_idx = j
# Swap the minimum value with the compared value

        input_list[idx], input_list[min_idx] = input_list[min_idx], input_list[idx]


def clash_rem(OC):#Take [C,O] as input
    clash_list =[]
    for carbon in OC[0]:
        for oxygen in OC[1]:
            if  HP_Interaction(carbon,oxygen) != False:
                water_num = oxygen[5]
                residue = carbon[5]
                carbon_name =  carbon[3] + '_' + carbon[2]
                oxygen_info =  oxygen[2]
                
                clash_list.append([water_num, oxygen_info,residue,carbon_name, HP_Interaction(carbon,oxygen)])
    selection_sort_alt(clash_list)

    return clash_list

def output_opts(clash,name,opt):
    terms = ['oxygen_name', 'carbon_name', 'distance']
    if opt == 'text' or 'txt':
        with open(name + '.txt', 'w') as txt_out:
            for elements in terms:
                txt_out.write(elements + '  ')
            txt_out.write('\n')
            for lines in clash:
                for info in lines:
                    txt_out.write(info + '\t')
                txt_out.write('\n')
        txt_out.close()
        return 0
    if opt == 'csv':
        with open(name + '.csv', 'w') as csv_out:
            csvwriter = csv.writer(csv_out)
            csvwriter.writerow(terms)
            for lines in clash:
                csvwriter.writerow(lines)
        csv_out.close()
        return 0
            
    else:
        print('Not a supported output')  
        return 0    
