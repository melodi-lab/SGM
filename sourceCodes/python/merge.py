import csv
from math import *
import os.path
from os import listdir
from os.path import isfile, join

#ch2file="test_scores_output_fast_yeast_01_ch2.txt.sumScoreIdent"
#ch3file="test_scores_output_fast_yeast_01_ch3.txt.sumScoreIdent"
#ch23file="test_scores_output_fast_yeast_01_ch23.txt.sumScoreIdent"
#ch2file="../output/test_scores_output_fast_bullseyeworm-ch4.txt.sumScoreIdent"
#ch3file="../output/test_scores_output_fast_bullseyeworm-ch3.txt.sumScoreIdent"
#ch3file="test_scores_output_fast_bullseyeworm-ch2.txt.sumScoreIdent2"
#ch23file="test_scores_output_fast_bullseyeworm-ch2345.txt.sumScoreIdent2"
#ch2file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
#ch3file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch4.txt.sumScoreIdent"
#ch23file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
#ch2file="/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch2/malaria-ch2-ident.txt"
#ch2file="malaria-ch2345-ident.txt"
#ch3file="/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch3/malaria-ch3-ident.txt"
#ch23file="malaria-ch2345-ident.txt"
def merge_indent(CHA,CHB,CHAB,fix_r=0):
	print("-------------------------");
	print("CHA   :"+CHA)
	print("CHB   :"+CHB)
	print("CHAB  :"+CHAB)
	input_file=csv.DictReader(open(CHA),delimiter='\t')
	spec2=[]
	for row in input_file:
		spec2.append((row["Kind"],int(row["Sid"]),row["Peptide"],float(row["Score"])))
	spec2sorted=sorted(spec2,key=lambda x:x[3],reverse=True)
	n_decoy=0;
	n_all=0;
	for row in spec2sorted:
		if row[0]=='t':
			n_all=n_all+1;
		else:
			n_all=n_all+1;
			n_decoy=n_decoy+1;
		if float(n_decoy)/float(n_all)>0.01:
			break
	n_all_A=n_all
	score2=(spec2sorted[n_all][3])
	
	input_file=csv.DictReader(open(CHB),delimiter='\t')
	spec3=[]
	for row in input_file:
		spec3.append((row["Kind"],int(row["Sid"]),row["Peptide"],float(row["Score"])))
	spec3sorted=sorted(spec3,key=lambda x:x[3],reverse=True)
	n_decoy=0;
	n_all=0;
	for row in spec3sorted:
		if row[0]=='t':
			n_all=n_all+1;
		else:
			n_all=n_all+1;
			n_decoy=n_decoy+1;
		if float(n_decoy)/float(n_all)>0.01:
			break
	n_all_B=n_all	
	score3=(spec3sorted[n_all][3])
	#print(score2/score3)
	#print(len(spec2sorted))
	#print(len(spec3sorted))
	
	sid_all=[]
	sid_2t=[]
	sid_2d=[]

	sid_3t=[]
	sid_3d=[]
	for row in spec2:
		sid_all.append(row[1])
		if row[0]=='t':
			sid_2t.append(row[1])
		else:
			sid_2d.append(row[1])
	for row in spec3:
		sid_all.append(row[1])
		if row[0]=='t':
			sid_3t.append(row[1])
		else:
			sid_3d.append(row[1])
	
#	sid_not_2=list(set(sid_all)-set(sid_2))
#	sid_not_3=list(set(sid_all)-set(sid_3))
	
	for row in list(set(sid_all)-set(sid_2t)):
		spec2.append(('t',row,'AAAA',-1000000000))
	for row in list(set(sid_all)-set(sid_2d)):
		spec2.append(('d',row,'AAAA',-1000000000))
	for row in list(set(sid_all)-set(sid_3t)):
		spec3.append(('t',row,'AAAA',-1000000000))
	for row in list(set(sid_all)-set(sid_3d)):
		spec3.append(('d',row,'AAAA',-1000000000))



	
	spec2=sorted(spec2,key=lambda x:x[0],reverse=True)
	spec2=sorted(spec2,key=lambda x:x[1])
	spec3=sorted(spec3,key=lambda x:x[0],reverse=True)
	spec3=sorted(spec3,key=lambda x:x[1])
#	print(spec2[0])
#	print(spec2[1])
#	print(spec2[2])
#	print(spec2[3])
#	print(spec2[4])
#	print(spec2[5])
#	print(spec3[0])
#	print(spec3[1])
#	print(spec3[2])
#	print(spec3[3])
#	print(spec3[4])
#	print(spec3[5])
	print(len(spec2))
	print(len(spec3))
	r1=1
	r2=1
	if score2*score3>1e-10:
		r1=score2/1000
		r2=score3/1000
	else:
		r1=1
		r2=1
	if fix_r:
		r1=1
		r2=1
	if len(spec3)<100:
		r1=1
		r2=1
	if len(spec2)<100:
		r1=1
		r2=1
	if r1<0.00000001:
		r1=1
		r2=1
		
	if r2<0.00000001:
		r1=1
		r2=1
	if r1/r2>100:
		r1=1
		r2=1
	if r2/r1>100:
		r1=1
		r2=1
	spec23=[]
	print r1,r2
	for row2,row3 in zip(spec2,spec3):
	#	if row2[0]!=row3[0]:
	#		print "asdfasdf"
		if row2[3]/r1>=row3[3]/r2:
			s2=row2[3]
			if s2>-1000000:
				s2=s2/r1
			row2r=(row2[0],row2[1],row2[2],s2)
			spec23.append(row2r);
	

		#	spec23.append(row2);
		if row2[3]/r1<row3[3]/r2:
			s3=row3[3]
			if s3>-1000000:
				s3=s3/r2

			row3r=(row3[0],row3[1],row3[2],s3)
			spec23.append(row3r);
	
	spec23sorted=sorted(spec23,key=lambda x:x[3],reverse=True)
	n_decoy=0;
	n_all=0;
	for row in spec23sorted:
		if row[0]=='t':
			n_all=n_all+1;
		else:
			n_all=n_all+1;
			n_decoy=n_decoy+1;
		if float(n_decoy)/float(n_all)>0.01:
			break
#	print ('n_A:%d n_b:%d n_ab:%d'%(n_all_A,n_all_B,n_all))
	
	
	
	f=open(CHAB,'w')
	f.write('Kind\tSid\tPeptide\tScore\n')
	for row in spec23:
		f.write('%s\t%s\t%s\t%.15f\n'%(row[0],row[1],row[2],row[3]))
def merge_mutilple_file(allfilenames,outputfile,fix_r=0):
	filenames=[]
	for row in allfilenames:
		if os.path.isfile(row) and os.path.getsize(row) > 0:
			filenames.append(row)
		else:
			print ("%s does not exist"%row)

	if len(filenames)<=1:
		print "less than two files-----------------------"
		return 
	
	merge_indent(filenames[0],filenames[1],outputfile,fix_r)		
	if len(filenames)==2:
		return
	for i in filenames[2:]:
		merge_indent(i,outputfile,outputfile,fix_r)		


#ch2file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2.txt.sumScoreIdent"
#ch3file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch3.txt.sumScoreIdent"
#ch4file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch4.txt.sumScoreIdent"
#ch5file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch5.txt.sumScoreIdent"
#ch6file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch6.txt.sumScoreIdent"
#ch7file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch7.txt.sumScoreIdent"
#ch2345file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
#merge_indent(ch2file,ch3file,ch2345file)		
#merge_indent(ch4file,ch2345file,ch2345file)		
#merge_indent(ch5file,ch2345file,ch2345file)		

#
#
#ch2file="/s1/wrbai/codes/msgfplus/malaria_TMT13/malaria-ch2/malaria-ch2-ident.txt"
#ch3file="/s1/wrbai/codes/msgfplus/malaria_TMT13/malaria-ch3/malaria-ch3-ident.txt"
#ch4file="/s1/wrbai/codes/msgfplus/malaria_TMT13/malaria-ch4/malaria-ch4-ident.txt"
#ch5file="/s1/wrbai/codes/msgfplus/malaria_TMT13/malaria-ch5/malaria-ch5-ident.txt"
#ch2345file="/s1/wrbai/codes/msgfplus/malaria_TMT13/malaria-ch2345-ident.txt"
#ch3file="/s1/wrbai/codes/msgfplus/malaria_jeffH/malaria-ch3/malaria-ch3-ident.txt"
#ch4file="/s1/wrbai/codes/msgfplus/malaria_jeffH/malaria-ch4/malaria-ch4-ident.txt"
#ch5file="/s1/wrbai/codes/msgfplus/malaria_jeffH/malaria-ch5/malaria-ch5-ident.txt"
#ch2345file="/s1/wrbai/codes/msgfplus/malaria_jeffH/malaria-ch2345-ident.txt"
#
#ch2file="/s1/wrbai/codes/msgfplus/linfeng_010211-2/linfeng-ch2/linfeng-ch2-ident.txt"
#ch3file="/s1/wrbai/codes/msgfplus/linfeng_010211-2/linfeng-ch2/linfeng-ch2-ident.txt"
#ch4file="/s1/wrbai/codes/msgfplus/linfeng_010211-2/linfeng-ch2/linfeng-ch2-ident.txt"
#ch5file="/s1/wrbai/codes/msgfplus/linfeng_010211-2/linfeng-ch2/linfeng-ch2-ident.txt"
#ch2345file="/s1/wrbai/codes/msgfplus/linfeng_010211-2/linfeng-ch2345-ident.txt"
#
#merge_indent(ch2file,ch3file,ch2345file)		
#merge_indent(ch4file,ch2345file,ch2345file)		
#merge_indent(ch5file,ch2345file,ch2345file)		
#merge_indent(ch6file,ch2345file,ch2345file)		
#merge_indent(ch7file,ch2345file,ch2345file)		

#ch2file="/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch2.txt.sumScoreIdent"
#ch3file="/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch3.txt.sumScoreIdent"
#ch4file="/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch4.txt.sumScoreIdent"
#ch5file="/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch5.txt.sumScoreIdent"
#
#
#ch2345file="/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
#merge_indent(ch2file,ch3file,ch2345file)		
#merge_indent(ch4file,ch2345file,ch2345file)		
#merge_indent(ch5file,ch2345file,ch2345file)		
#
#dir="/n/trombone"
#ch2file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch2_test.txt.sumScoreIdent"
#ch3file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch3_test.txt.sumScoreIdent"
#ch4file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch4_test.txt.sumScoreIdent"
#ch5file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch5_test.txt.sumScoreIdent"
#ch2345file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch2345_test.txt.sumScoreIdent"
##
##
#merge_indent(ch2file,ch3file,ch2345file)		
#merge_indent(ch4file,ch2345file,ch2345file)		
#for i in range(143,200):
#	id=i+1
#
#	dir="/n/trombone"
#	ch2file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng-%d_ch2_test.txt.sumScoreIdent"%(id)
#	ch3file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng-%d_ch3_test.txt.sumScoreIdent"%(id)
#	ch4file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng-%d_ch4_test.txt.sumScoreIdent"%(id)
#	ch2345file="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng-%d-ch2345_test.txt.sumScoreIdent"%(id)
#	merge_indent(ch2file,ch3file,ch2345file)		
#	merge_indent(ch4file,ch2345file,ch2345file)		
#
for i in range(20):

	id=i+1
	filenames=[]
	dir="/n/trombone"
	for j in range(4):
		c=j+2
		#filenames.append("/n/trombone/s1/wrbai/codes/msgfplus/linfeng_9_%d-1/linfeng-ch%d/linfeng.txt"%(id,c))
#		filenames.append("/n/trombone/s1/wrbai/codes/msgfplus/kim-%d/kim-ch%d/kim.txt"%(id,c))
	#	filenames.append("/n/trombone/s1/wrbai/codes/msgfplus/linfeng-%d/linfeng-ch%d/linfeng.txt"%(id,c))
	#	if id==17 and c==5:
	#		continue;
	#	if id==9 and c==5:
	#		continue
		#filenames.append("/n/trombone/s1/wrbai/codes/output/kim-random20-%d-0-ch%d.txt"%(id,c))
#		filenames.append("/n/trombone/s1/wrbai/codes/output/postbcb-mean/plasm-%d-ch%d-lambda1-7-lambda2-15.txt"%(id,c))
		filenames.append("/n/trombone/s1/wrbai/codes/output/postbcb/2016_6_29/Linfeng1/plasm-%d-ch%d-lambda1-7-lambda2-15.txt"%(id,c))
	#	filenames.append("/n/trombone/s1/wrbai/codes/msgfplus/plasm-%d/plasm-ch%d/plasm.txt"%(id,c))

#	outputfile="/n/trombone/s1/wrbai/codes/msgfplus/linfeng_9_%d-1/linfeng.txt"%(id)
	#outputfile="/n/trombone/s1/wrbai/codes/output/postbcb-mean/plasm-%d-lambda1-7-lambda2-15.txt"%(id)
	outputfile="/n/trombone/s1/wrbai/codes/output/postbcb/2016_6_29/Linfeng1/Linfeng-%d-lambda1-7-lambda2-15.txt"%(id)
#	outputfile="/n/trombone/s1/wrbai/codes/msgfplus/plasm-%d/plasm.txt"%(id)
#	outputfile="/n/trombone/s1/wrbai/codes/msgfplus/kim-%d/kim.txt"%(id)
#	outputfile="/n/trombone/s1/wrbai/codes/msgfplus/linfeng-3/linfeng.txt"
	merge_mutilple_file(filenames,outputfile,0)
##
##
#merge_indent(ch4file,ch2345file,ch2345file)		
#
#m
#merge_indent(ch5file,ch2345file,ch2345file)		
#
##
#
#
##ch2file="/s1/wrbai/codes/output/test_scores_output_fast_linfeng_112710_ch2_test.txt.sumScoreIdent"
##ch3file="/s1/wrbai/codes/output/test_scores_output_fast_linfeng_112710_ch3_test.txt.sumScoreIdent"
##ch4file="/s1/wrbai/codes/output/test_scores_output_fast_linfeng_112710_ch4_test.txt.sumScoreIdent"
##ch5file="/s1/wrbai/codes/output/test_scores_output_fast_linfeng_112710_ch5_test.txt.sumScoreIdent"
##ch2345file="/s1/wrbai/codes/output/test_scores_output_fast_linfeng_112710_ch2345_test.txt.sumScoreIdent"
##
##
##merge_indent(ch2file,ch3file,ch2345file)		
##merge_indent(ch4file,ch2345file,ch2345file)		
##merge_indent(ch5file,ch2345file,ch2345file)		
##
#

#ch2file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-TMT11_ch3_test.txt.sumScoreIdent"
#ch3file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-TMT11_ch3_test.txt.sumScoreIdent"
#ch4file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-TMT11_ch4_test.txt.sumScoreIdent"
#ch5file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-TMT11_ch5_test.txt.sumScoreIdent"
#ch2345file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-TMT11_ch2345_test.txt.sumScoreIdent"
#
#
#merge_indent(ch2file,ch3file,ch2345file)		
#merge_indent(ch4file,ch2345file,ch2345file)		
#merge_indent(ch5file,ch2345file,ch2345file)		
#






