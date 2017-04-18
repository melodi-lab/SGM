##read indent file and calcualte n at q=0.01 and scores
import random
import csv
from math import *
file_name="/n/trombone/s1/wrbai/codes/bullseyeworm-ch2/xcorrIdent-bullseyeworm-0-ch2.txt"
file_name="/n/trombone/s1/wrbai/codes/output/bullseyeworm_ch2_by.txt.sumScoreIdent"
file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_bullseyeworm-ch2.txt.sumScoreIdent"
file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_bullseyeworm-ch3.txt.sumScoreIdent"
file_name="/s1/wrbai/codes/malaria-ch3/xcorrIdent-malaria-0-ch3.txt"
file_name="/n/trombone/s1/wrbai/codes/linfeng-ch3/xcorrIdent-linfeng-0-ch3.txt"
file_name="/n/trombone/s1/wrbai/codes/malaria-ch3/xcorrIdent-malaria-0-ch3.txt"
file_name="/n/trombone/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch3/malaria-ch3-ident.txt"
#file_name="/n/trombone/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch5/malaria-ch5-ident.txt"
file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch3.txt.sumScoreIdent"
#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch4.txt.sumScoreIdent"
#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch5.txt.sumScoreIdent"
#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_bullseyeworm-ch2.txt.sumScoreIdent"
#file_name="/n/trombone/s1/wrbai/codes/bullseyeworm-ch1/xcorrIdent-bullseyeworm-0-ch1.txt"
#file_name="/n/trombone/s1/wrbai/codes/john_idents/msgfPlusIdents/msgfPlus-bullseyeWorm-idents/bullseyeWorm-ch2-ident.txt"
#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2.txt.sumScoreIdent"
#file_name="/n/trombone/s1/wrbai/codes/malaria-ch3/xcorrIdent-malaria-0-ch3.txt"
#file_name="/n/trombone/s1/wrbai/codes/john_idents/drip-bullseyeWorm-ch12345-ident.txt"
def delted_some(target):
	if len(target)<=1:
		return target
	t=sorted(target,key=lambda x:x[1])
	t2=[]
	i=0
	for row in t:
		if not row[1]==i:
			t2.append(row)
			i=row[1]
			continue
		if row[3]>t2[-1][3]:
			a=t2.pop()
			print a 
			t2.append(row)
			i=row[1]
	return t2		

def read_indent(FILE_NAME,FILE_NAME2):
	input_file=csv.DictReader(open(FILE_NAME),delimiter='\t')
	print(file_name)
	target_decoy=[]
	for row in input_file:
		target_decoy.append((row["Kind"],int(row["Sid"]),row["Peptide"],float(row["Score"])))
#	target_decoy=sorted(target_decoy,key=lambda x:x[0])
#	target_decoy=sorted(target_decoy,key=lambda x:x[1])
	target=[x for x in target_decoy if x[0]=='t']
	decoy=[x for x in target_decoy if x[0]=='d']
	set1=set([x[1] for x in target])
	set2=set([x[1] for x in decoy])
	print len(set1)
	print len(set2)
	
	if not len(set1)==len(set2):
		target=delted_some(target)
		decoy=delted_some(decoy)
		set3=set1.union(set2);
		set_no_t=set3-set1;
		set_no_d=set3-set2;
		for row in set_no_t:
			target.append(('t',row,'AAAAAAA',-1000000))
		for row in set_no_d:	
			decoy.append(('d',row,'AAAAAAA',-1000000))
	target=sorted(target,key=lambda x:x[1])
	decoy=sorted(decoy,key=lambda x:x[1])
	print len(target)
	print len(decoy)
	target_decoy2=[]
	for row,row2 in zip(target,decoy):
		ind=0
		if row[3]-row2[3]>1e-7:
			ind=1
		elif row2[3]-row[3]>1e-7:
			ind=0
		else:
			ind=random.randint(0, 1)
			
		if ind:
			target_decoy2.append((row[0],row[1],row[2],row[3]))
			target_decoy2.append((row2[0],row2[1],row2[2],-1000000))
		else:
			target_decoy2.append((row[0],row[1],row[2],-1000000))
			target_decoy2.append((row2[0],row2[1],row2[2],row2[3]))
			
	f=open(FILE_NAME2,'w')		
	f.write('Kind\tSid\tPeptide\tScore\n')
	for row in target_decoy2:
		f.write('%s\t%s\t%s\t%.10f\n'%(row[0],row[1],row[2],row[3]))
	
#
file_name="malaria-1.txt"
file_name2="malaria-2.txt"
read_indent(file_name,file_name2)
exit()
name="plasm"
name2="plasm2"
file_name="output/bcb2016/%s_bipartite.txt"%name
file_name2="output/bcb2016/%s_bipartite.txt"%name2
read_indent(file_name,file_name2)
exit()
file_name="output/bcb2016/%s_msgfPlus.txt"%name
file_name2="output/bcb2016/%s_msgfPlus.txt"%name2
read_indent(file_name,file_name2)

file_name="output/bcb2016/%s_tide.txt"%name
file_name2="output/bcb2016/%s_tide.txt"%name2
read_indent(file_name,file_name2)

file_name="output/bcb2016/%s_pvalue.txt"%name
file_name2="output/bcb2016/%s_pvalue.txt"%name2
read_indent(file_name,file_name2)



#file_name="run_crux/kim-xcorr-concat4/kim_tide.txt"
#file_name2="output/ismb2016/kim_xcorr.txt"
#read_indent(file_name,file_name2)


#file_name="run_crux/bill-result/wenruo-result/Linfeng_tide.txt"
#file_name2="run_crux/bill-result/wenruo-result/Linfeng2_tide.txt"
#read_indent(file_name,file_name2)


#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng_080510_HapMap9_%d_ch%s_test.txt.sumScoreIdent"%(6,5)

#file_name="output/ismb2016/human2_msgfPlus.txt"
#file_name2="output/ismb2016/human3_msgfPlus.txt"
#read_indent(file_name,file_name2)
#
#
#file_name="output/ismb2016/human_tide.txt"
#file_name2="output/ismb2016/human1_tide.txt"
#read_indent(file_name,file_name2)
#
#file_name="output/ismb2016/human2_pvalue.txt"
#file_name2="output/ismb2016/human3_pvalue.txt"
#read_indent(file_name,file_name2)
#file_name="run_crux/malaria-TMT10_output_p-ch3/malaria_tide_pvalue.txt"
#file_name2="run_crux/malaria-TMT10_output_p-ch3/malaria_tide_pvalue2.txt"
#read_indent(file_name,file_name2)



#
#file_name="output/ismb2016/human_bipartite.txt"
#file_name2="output/ismb2016/human1_bipartite.txt"
#read_indent(file_name,file_name2)
#
#
##read_indent(file_name)
#file_name="output/ismb2016/malaria_bipartite.txt"
#file_name2="output/ismb2016/malaria2_bipartite.txt"
#read_indent(file_name,file_name2)
#
#file_name="output/ismb2016/malaria_msgfPlus.txt"
#file_name2="output/ismb2016/malaria2_msgfPlus.txt"
#read_indent(file_name,file_name2)
#
#
#file_name="output/ismb2016/malaria_tide.txt"
#file_name2="output/ismb2016/malaria2_tide.txt"
#read_indent(file_name,file_name2)
#
#file_name="output/ismb2016/malaria_pvalue.txt"
#file_name2="output/ismb2016/malaria2_pvalue.txt"
#read_indent(file_name,file_name2)
#

