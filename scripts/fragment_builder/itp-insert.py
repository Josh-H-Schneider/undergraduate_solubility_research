import os
import fileinput

nmol = 1
frag_prefix = "test"
fragitp = '"' + frag_prefix + '.itp"'
count = 0

project_combo_directory = './combo_folder/'

combo = os.path.join(project_combo_directory, f'combo.top')

   
for line in fileinput.FileInput(combo,inplace=1): 
    if count < 1:
        if "#include" in line:
            line=line.replace(line,line+"\n; Include fragment topology\n#include " + fragitp + "\n")
            count += 1
    print (line, end='')
    
    
with open(combo, 'a') as filedata:
    filedata.write("\n" + frag_prefix + "           "  + str(nmol))
    filedata.close()