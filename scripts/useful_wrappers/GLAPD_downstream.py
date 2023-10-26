import sys
import csv

with open(sys.argv[1], 'r') as fr, open(sys.argv[2], 'w', newline='') as fw:
    writer = csv.writer(fw)
    head = ['Index', 'ID', 'Seq', 'Pos', 'Len']
    writer.writerow(head)
    index = 0
    for line in fr.readlines():
        line_list = line.strip('\t').strip('\n').strip().split(':')
        if '' in line_list:
            index += 1
        if len(line_list) == 2:
            if line_list[-1].strip() == "NULL":
                writer.writerow([index, line_list[0], 'NULL'])
        if len(line_list) == 5:
            id = line_list[0]
            pos = int(line_list[2].split(',')[0])
            length = int(line_list[3].split(' ')[0])
            seq = line_list[-1]
            writer.writerow([index, id, seq, pos, length])
        
        
        
        
