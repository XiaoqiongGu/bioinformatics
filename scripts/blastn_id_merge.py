### this script is going to extract fasta ids containing both forward and reverse ids.
### python blastn_id_merge.py blast.output extracted.id.txt

import pandas as pd
import sys


blast = pd.read_csv(sys.argv[1],sep='\t',header=None)
blast.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','qcovs']
r1tem = blast[blast['qseqid'] == 'R1tem']['sseqid'].to_list()
f3tem = blast[blast['qseqid'] == 'F3tem']['sseqid'].to_list()

r1tem_f3tem = list(set(r1tem)&set(f3tem))

with open(sys.argv[2],'w') as f:
    for i in r1tem_f3tem:
        f.write(i + '\n')