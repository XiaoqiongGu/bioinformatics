
import sys
import numpy as np

f1 = open(sys.argv[1])
lines = f1.readlines()
f1.close()

f2 = open(sys.argv[2],'w')
sig_char = lines[0].split('\t')
out_char1 = sig_char[0]+'\t'
for i in range(2,len(sig_char)-1):
    out_char1 += sig_char[i]+'\t'
out_char1 += sig_char[len(sig_char)-1]

f2.write(out_char1)
newlines = [lines[0]]



contig_num = len(lines)-1
sample_num = len(sig_char)-2


fin_mat = np.zeros((contig_num,sample_num+1))

for i in range(contig_num):
    sig_char = lines[i+1].split('\t')
    for j in range(sample_num+1):
        fin_mat[i,j] = int(sig_char[j+1])
        
length_vect = fin_mat[:,0]
length_mat = np.tile(length_vect,[sample_num,1]).T
sample_mat = fin_mat[:,1:sample_num+1]
sum_vect = np.sum(sample_mat,0)
sum_mat = np.tile(sum_vect,[contig_num,1])

out_mat = sample_mat*1e9/sum_mat/length_mat
out_mat[np.isnan(out_mat)] = 0

for i in range(1,contig_num+1):
    sig_char = lines[i].split('\t')
    out_char = sig_char[0]+'\t'
    for j in range(sample_num-1):
        out_char += str(out_mat[i-1,j])+'\t'
    out_char += str(out_mat[i-1,sample_num-1])+'\n'
    f2.write(out_char)
    
f2.close()
