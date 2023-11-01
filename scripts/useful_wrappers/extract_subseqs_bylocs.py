import sys


def produce_extract_loc_bash(input_file, output_file):
    candi_num = {}
    dic2 = {}
    n = 0
    with open(input_file, 'r') as f:
        for i in f.readlines():
            locs = i.split('\t')

            labs, loc1, loc2 = locs[1], int(locs[8]), int(locs[9])

            #if labs not in dup_set:
            #    continue

            if labs not in candi_num:
                candi_num[labs] = [loc1, loc2]
            else:
                if labs in dic2:
                    continue
                candi_num[labs].extend([loc1, loc2])
                max_num = max(candi_num[labs])
                min_num = min(candi_num[labs])
                dic2[labs] = [min_num, max_num]

    with open(output_file, 'w') as f:
        fasta_file = output_file.split('blast_')[1].replace(
            '.sh', '.linearized.fasta ')
        for key, value in dic2.items():
            f.write('samtools faidx ' + fasta_file + key + ':' +
                    str(value[0]) + '-' + str(value[1]) + ' >> ' + output_file.replace('sh', 'fasta') + '\n')


produce_extract_loc_bash(sys.argv[1], sys.argv[2])
