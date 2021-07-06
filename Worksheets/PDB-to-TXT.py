import pandas as pd

# for pName in ['6acd', '7krq', '7lwt', '7lws', '7lww', '7lyn', '6zgh', '6zge',
#         '6zgi', '6xkl', '7lyl', '7kdk', '6zgg', '7m8k', '7mjg']:
#     with open(f'Proteins-PDB/{pName}.pdb', 'r') as f:
#         pdb_lines = f.readlines()

#     txt_lines = []
#     for line in pdb_lines:
#         words = line.split()
#         if not words[0] == 'ATOM': continue
#         if words[2] == 'CA':
#             if words[4][0] == 'A':
#                 txt_lines.append(line)

#     with open(f'Proteins-PDB/Txt-Files/{pName}.txt', 'w') as f:
#         for line in txt_lines:
#             f.writelines(line)

proteins = ['6acd', '7krq', '7lwt', '7lws', '7lww', '7lyn', '6zgh', '6zge',
        '6zgi', '6xkl', '7lyl', '7kdk', '6zgg', '7m8k', '7mjg']
df = pd.DataFrame()
data = dict()
for pName in proteins:
    amino_indices = []
    with open(f'Proteins-PDB/Txt-Files/{pName}.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            amino_indices.append([int(line[22:26])])
    data[pName] = pd.Series(amino_indices)

df = pd.concat(data, axis = 1)

with open("AminoAcid-Indices.csv", mode='a', newline='') as f:
    df.to_csv(f, header=f.tell()==0, index=False)
