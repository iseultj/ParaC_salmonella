#!/usr/bin/env python3

###USAGE: python3 completedeletioncode.py <fasta_filename> <percentage_of_deletion> <output_filename> <0 for variant sites only; anyother number above 0 for including invariant sites> '####
## some edits by iseult so partial deletion works

'''
Script filters sites that contain N in any seq in multifasta
4th argument set to 0, will filter in addition for variable sites only (i.e. >1 non-N character!)
'''

import argparse,sys

''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Filters MultiFasta for CompleteDeletion (i.e. removes sites w/ "N"), and for invariable sites (if specified)
                                ''',
                                epilog="Questions or comments to Adytia ;) ")
parser.add_argument("-f", dest="fasta_file", help="Input multifasta", type=argparse.FileType('rt'))
parser.add_argument("-o", "--out_fasta", dest="out", help="Output multifasta", type=argparse.FileType('w'),default=sys.stdout)
parser.add_argument("-l", "--seq_length", dest="seq_length", help="Number of bases to be printed per line in output fasta [100]", type=int,default=100)
parser.add_argument("-p", "--percentage", dest="percentage", help="Percentage for partial deletion, default [100] is complete deletion ", type=int,default=100)
parser.add_argument('-n', dest="allSites", help="Set flag to print non-varaible sites too [False]", action='store_true', default=False) # default= not necessary as implied by action
args = parser.parse_args()


out = args.out


## read the file and splits it and make a hash table####
with args.fasta_file as f:
    g=f.read().strip().split('>')

Alignment={}
for i in g[1:]:
    Alignment[i.split('\n')[0]]=[''.join(i.split('\n')[1:])]

### Calculate the length of Alignment ###

length=[]
for i in Alignment:
    length.append(len(Alignment[str(i)][0][0:]))
###select for seqeunces within the percentage of deletion ####

p=list(set(length))
if len(p)>1:
	sys.exit(' All the sequences in the alignment are of not same length')



nucbases=['A','T','G','C']
percentage=int(args.percentage)

def var_invar(s):
    '''Identifies if a site is Variable or Invariable and makes a count of the bases'''
    F=[]
    F.append(s)
    if (len(list(set(s))))>1:
        F.append('var')
        F.append(s.count('A')/len(s))
        F.append(s.count('T')/len(s))
        F.append(s.count('G')/len(s))
        F.append(s.count('C')/len(s))
        F.append(1-((s.count('A')/len(s))+(s.count('T')/len(s))+(s.count('G')/len(s))+(s.count('C')/len(s))))
    else:
        F.append('invar')
        F.append(s.count('A')/len(s))
        F.append(s.count('T')/len(s))
        F.append(s.count('G')/len(s))
        F.append(s.count('C')/len(s))
        F.append(1-((s.count('A')/len(s))+(s.count('T')/len(s))+(s.count('G')/len(s))+(s.count('C')/len(s))))
    return(F)



## build alignment matrix (lol) and filter for N (i.e. complete deletion)
n=[]
n1=[]
m1=[]
A=[0,0]
T=[0,0]
G=[0,0]
C=[0,0]
O=[0,0]
for j in range(0,p[0]):
	s=[]
	for i in Alignment:
		s.append(Alignment[i][0][j])
	if (float(((s.count('A')+s.count('T')+s.count('G')+s.count('C'))/len(s))*100))>=float(percentage):
		w1=var_invar(s)
		if w1[1]=='var':
			n.append(w1[0])
			A[0]=A[0]+w1[2]
			T[0]=T[0]+w1[3]
			G[0]=G[0]+w1[4]
			C[0]=C[0]+w1[5]
			O[0]=O[0]+w1[6]
			m1.append(w1[0])
		elif w1[1]=='invar':
			n1.append(w1[0])
			A[1]=A[1]+w1[2]
			T[1]=T[1]+w1[3]
			G[1]=G[1]+w1[4]
			C[1]=C[1]+w1[5]
			O[1]=O[1]+w1[6]
			m1.append(w1[0])
	else:
		pass


#    if 'N' in s:
#        pass
#    else:
#        m.append(s)

## filter for variable sites
#m1=[]
#if args.allSites==False:
#    for o in m:
#        if (len(list(set(o))))>1:
#            m1.append(o)
#else:
#    m1=m

                       #del m # remove m to save resources (in case alignment is large!)

## re-build fasta
					   #Alignment1={}
					   #a=0
					   #for k in Alignment:
					   #n2=[]
					   #    for l in m1:
					   #        n2.append(l[a])
					   #    Alignment1[k]=[''.join(n2)]
					   #    a=a+1
                       
					   #def rebuild(m1):
					   #    Alignment1={}
					   #    a=0
#    for k in Alignment:
#        n2=[]
#        for l in m1:
#            n2.append(l[a])
#		Alignment1[k]=[''.join(n2)]
#		a=a+1
#    return(Alignment1)

def rebuild(m1):
	Alignment1={}
	a=0
	for k in Alignment:
		n2=[]
		for l in m1:
			n2.append(l[a])
		Alignment1[k]=[''.join(n2)]
		a=a+1
	return(Alignment1)


if args.allSites==False:
	Alignment1=rebuild(n)
	print('printing alignment with only variable sites')
else:
	Alignment1=rebuild(m1)
	print('printing alignment with both variable and invariable sites')


##Reporting the Sequences and proportion of Nucleotides
print("Total Number of Sites in Alignment before deletion :",p[0])
print('\n')
print("Total Number of Sites in Alignment after "+str(percentage)+" deletion",len(m1))
print('\n')
print("Total Number of Variable sites : ",len(n))
print("Composition of A :",str(float((A[0]/len(n)))))
print("Composition of T :",str(float(T[0]/len(n))))
print("Composition of G :",str(float((G[0]/len(n)))))
print("Composition of C :",str(float((C[0]/len(n)))))
print("Composition of Others :",str(float((O[0]/len(n)))))
print('\n\n')
print("Total Number of Invariable sites : ",len(n1))
print("Number of A sites:",A[1])
print("Number of T sites:",T[1])
print("Number of G sites:",G[1])
print("Number of C sites:",C[1])
print("Number of Other sites:",O[1])

### Write to the file ###
seq_length = args.seq_length
o = args.out
for i in Alignment1:
    # o.write('>'+str(i)+'\n'+Alignment1[str(i)][0]+'\n')
    o.write('>'+str(i)+'\n')
    sequence = Alignment1[str(i)][0]
    while len(sequence) > 0:
        o.write(sequence[:seq_length]+'\n')
        sequence = sequence[seq_length:]

o.close()
print("Successfully completed writing the alignment to ",args.out.name)
        

### END #####






