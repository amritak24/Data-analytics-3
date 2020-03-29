#Output of the following code is the gene count for intersection with respective datasets 
#and data for genes for which the count for smokers increase or decrease 
#in smokers wrt non smokers in males and females respectively and the final set of genes
#whose value increase or decrease due to smoking for males and females respectively.
import numpy as np
import pandas as pd
import math
from  scipy import stats
import matplotlib.pyplot as plt
import statistics

df=pd.read_csv("./../data/Raw Data_GeneSpring.txt", sep="\t")

genecodes=df.iloc[:,0:1].values.tolist()
data=df.iloc[:,1:49]
#print(len(genecodes))

A=np.zeros(shape=(48,4))
#print(A)
for i in range(12):
  A[i][0]=1
  A[i+12][1]=1
  A[i+24][2]=1
  A[i+36][3]=1

A_null=np.zeros(shape=(48,4))
#print(A)
for i in range(12):
  A_null[i][0]=1
  A_null[i][2]=1
  A_null[i+12][0]=1
  A_null[i+12][3]=1
  A_null[i+24][1]=1
  A_null[i+24][2]=1
  A_null[i+36][1]=1
  A_null[i+36][3]=1

data_list=data.values.tolist()
#print(len(data_list[0]))

A= np.matrix(A)
A_null= np.matrix(A_null)
F_stat=[]
for row in data_list:
  rowm=np.matrix(row)
  x=np.matmul(np.matmul(A, np.linalg.pinv(np.matmul(A.T, A))),A.T)
  y=np.matmul(np.matmul(A_null, np.linalg.pinv(np.matmul(A_null.T, A_null))),A_null.T)
  F=np.matmul(np.matmul(rowm,(x-y)),rowm.T)/np.matmul(np.matmul(rowm,(np.identity(48)-x)),rowm.T)
  F_mod= float(F) * (48-np.linalg.matrix_rank(A))/(np.linalg.matrix_rank(A)-np.linalg.matrix_rank(A_null))
  F_stat.append(F_mod)
#print(len(F_stat))

dfn=(np.linalg.matrix_rank(A)-np.linalg.matrix_rank(A_null))
dfd=(48-np.linalg.matrix_rank(A))
p_val= 1 - stats.f.cdf(F_stat,dfn,dfd)
#print(len(p_val))

p_genes=[]
p_val_req=[]
i_list=[]
count=0
#print(1)
for i in range(len(p_val)):
  #print(2)
  if p_val[i] < 0.05 :
    #print(1)
    count=count+1
    i_list.append(i)
    p_genes.append(genecodes[i][0])
    p_val_req.append(p_val[i])
#print(count)

p_val_new=[]
for i in range(len(p_val)):
  if  not math.isnan(p_val[i]) :
    p_val_new.append(p_val[i])

#print(len(p_val) , len(p_val_new))

plt.hist(p_val_new,bins=20)
plt.xlabel('p_value')
plt.ylabel('Frequency')
plt.title('Histogram of p_value')
plt.savefig("p_value.jpg")

genesymbol=[]
genesymbols=df.iloc[:,49:50].values.tolist()
for i in genesymbols:
  genesymbol.append(i[0])
#print(genesymbol)

gs=[]
for i in i_list:
    gs.append(genesymbol[i])
#print((gs))

def find_mean(inter):
  male_up=[]
  male_down=[]
  female_up=[]
  female_down=[]
  for i in inter:
    a=np.mean(data_list[i][0:12])
    b=np.mean(data_list[i][12:24])
    c=np.mean(data_list[i][24:36])
    d=np.mean(data_list[i][36:])
    if(a<b):
      male_up.append(genesymbol[i])
    if(a>b):
      male_down.append(genesymbol[i])
    if(c<d):
      female_up.append(genesymbol[i])
    if(c>d):
      female_down.append(genesymbol[i])
  return male_up,male_down,female_up,female_down

new_i=[]
nkcell=pd.read_csv("./../data/NKCellCytotoxicity.txt")
nkcell=nkcell['KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY'].loc[1:].values.tolist()
#print(len(set(nkcell)))
inter_nkcell=set(gs).intersection(set(nkcell))
print("\ngene count NKCellCytotoxicity data : ", len(inter_nkcell))
print("intersecting genes for NKCellCytotoxicity data :", inter_nkcell)
 
with open("./../data/FreeRadicalResponse.txt", "r") as f:
  with open("newfile.txt", "w") as output: 
        next(f)
        next(f)
        for line in f:
                output.write(line)

new=pd.read_csv("newfile.txt",header=None)
#print(new[0].values)
inter_frad=set(gs).intersection(set(new[0].values))

print("\ngene count FreeRadicalResponse data : ", len(inter_frad)) 
print("intersecting genes for FreeRadicalResponse data :", inter_frad)


with open("./../data/DNARepair1.txt", "r") as f:
  with open("newfile.txt", "w") as output: 
        next(f)
        next(f)
        for line in f:
                output.write(line)

new=pd.read_csv("newfile.txt",header=None)
#print(new[0].values)
inter_dna_rep=set(gs).intersection(set(new[0].values))

print("\ngene count DNARepair1 data : ", len(inter_dna_rep)) 
print("intersecting genes for DNARepair1 data :", inter_dna_rep)


with open("./../data/XenobioticMetabolism1.txt", "r") as f:
  with open("newfile.txt", "w") as output: 
        next(f)
        next(f)
        for line in f:
                output.write(line)

new=pd.read_csv("newfile.txt",header=None)
#print(new[0].values)
inter_xeno=set(gs).intersection(set(new[0].values))

print("\ngene count XenobioticMetabolism1 data : ", len(inter_xeno)) 
print("intersecting genes for XenobioticMetabolism1 data :", inter_xeno)


intersect= list(inter_nkcell) + list(inter_frad) + list(inter_dna_rep) + list(inter_xeno)
print("\nTotal number of intersecting genes: ", len(list(set(intersect))))
for i in intersect:
    li= df.index[df['GeneSymbol']==i]
    #print(li)
    new_i.extend(list(li))
#print(inter_nkcell)
#print(len(inter_nkcell))
male_up,male_down,female_up,female_down= find_mean(new_i)
print("\nmale smoker value higher than male non smoker for these genes : \n",set(male_up))
print("\nmale smoker value lower than male non smoker for these genes : \n",set(male_down))
print("\nfemale smoker value higher than female non smoker for these genes : \n",set(female_up))
print("\nfemale smoker value lower than female non smoker for these genes : \n",set(female_down))

n0=1
q=0.05
l=0
p_sort=np.sort(p_val)
for i in range(len(p_val)):
  if (p_sort[i]*n0/i) <= q :
    l=i
#print(l)
