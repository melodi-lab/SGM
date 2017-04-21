##read indent file and calcualte n at q=0.01 and scores
import csv
from math import *
#file_name="/n/trombone/s1/wrbai/codes/john_idents/drip-bullseyeWorm-ch12345-ident.txt"
def read_indent(FILE_NAME):
	input_file=csv.DictReader(open(FILE_NAME),delimiter='\t')
	print(file_name)
	target_decoy=[]
	for row in input_file:
		target_decoy.append((row["Kind"],int(row["Sid"]),row["Peptide"],float(row["Score"])))
	sorted_target_decoy=sorted(target_decoy,key=lambda x:-x[3])
#	print sorted_target_decoy[0:3]
	n_t=0;
	n_d=0;
	score=0;
	n=-1;
	for row in sorted_target_decoy:
		if row[0]=='t':
			n_t=n_t+1
		elif row[0]=='d':
			n_d=n_d+1
		if float(n_d)/float((n_t+n_d))>0.01:
			score=row[3]
			n=n_t+n_d
			break
	#print(score)
	#print(n)
	#print(n_t)
	#print(n_d)
	n_t=0;
	n_d=0;
	score=0;
	n=-1;
	for row in sorted_target_decoy:
		if row[0]=='t':
			n_t=n_t+1
		elif row[0]=='d':
			n_d=n_d+1
		if float(n_d)/float((n_t+n_d))>0.04:
			score=row[3]
			n=n_t+n_d
			break
	#print(score)
	#print(n)
	#print(n_t)
	#print(n_d)
	
	#print(sid)
	q_range=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
        q_range=range(100)
        q_range=[float(i)/1000 for i in q_range]
	n_range=[]
	score_range=[]
	for q in q_range:
		n_t=0;
		n_d=0;
		score=0;
		n=-1;
		for row in sorted_target_decoy:
			if row[0]=='t':
				n_t=n_t+1
			elif row[0]=='d':
				n_d=n_d+1
			if float(n_d)/float((n_t+n_d))>q:
				score=row[3]
				n=n_t+n_d
				break
		n_range.append(n_t)
		score_range.append(score)
#	print "q:\t",
#	for q in q_range:
#		print"%f\t"%(q),
	print "\nn:\t",
	for n in n_range:
		print"%f\t"%(n),
#	print "\nscore:\t",
#	for s in score_range:
#		print"%f\t"%(s),
		
	print "\narea=%d"%(sum(n_range)/11)
	#print(sid)

print("//////////////////////////////////////////////////////////////////////////////////////////")	
print("//////////////////////////////////////////////////////////////////////////////////////////")	
print("//////////////////////////////////////////////////////////////////////////////////////////")	
print("//////////////////////////////////////////////////////////////////////////////////////////")	

print("//////////////////////////////////////////////////////////////////////////////////////////")	
import sys
file_name= sys.argv[1]
read_indent(file_name)
