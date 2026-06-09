import sys, os
# fix gapps in snp-sites output from gubbins
f = open(sys.argv[1],'r')
for line in f.readlines():
    line = line.strip()
    if line.startswith('>'):
        print(line)
    else:
        subbed = line.replace("-","N")
        print(subbed)

f.close()
