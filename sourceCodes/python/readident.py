##read indent file and calcualte n at q=0.01 and scores
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

file_name="temp/malaria-2.txt"
read_indent(file_name)
exit()

file_name = "output/postbcb/plasm-1-ch3-lambda1-7-lambda2-15.txt"

read_indent(file_name)
exit()
file_name = "output/postbcb/2016_12_2/1/malarial/1.txt"

read_indent(file_name)
file_name = "output/postbcb/2016_11_8/malarial/1.txt"
read_indent(file_name)

file_name = "output/postbcb/2016_11_8/malarial/0.txt"

read_indent(file_name)

file_name = "output/postbcb/2016_11_8/malarial/1.txt"

read_indent(file_name)
file_name="msgfplus/plasm-10/plasm-ch3/plasm.txt"
read_indent(file_name)
exit()
file_name="matlab_codes/8_29_2016/plasm-10-ch3.txt"
print read_indent(file_name)
file_name="output/postbcb/2016_8_27/malaria/1.txt"
print read_indent(file_name)
file_name="output/postbcb/2016_8_25/malaria3/1.txt"
print read_indent(file_name)


file_name="backup05272014/310/test_scores_output_fast_Linfeng-10_ch3_test.txt.sumScoreIdent"
read_indent(file_name)

exit()
file_name="output/postbcb/2016_8_18/mann/4.txt"
read_indent(file_name)
flie_name="output/postbcb/2016_8_18/malaria/1.txt"

file_name="output/postbcb/2016_8_2/kim/plasm-10-ch3-lambda1-1.0-lambda2-1.0-power-2.5.txt"
read_indent(file_name)

file_name="output/postbcb/2016_8_2/kim/plasm-10-ch3-lambda1-1.0-lambda2-0.0-power-2.5.txt"
read_indent(file_name)
file_name="output/postbcb/2016_8_2/kim/plasm-10-ch3-lambda1-0.0-lambda2-1.0-power-2.5.txt"

read_indent(file_name)
file_name="output/postbcb/2016_7_27/malaria/plasm-10-ch3-lambda1-0.0-lambda2-1.0-power-2.0.txt"
read_indent(file_name)
file_name="output/postbcb/2016_7_27/malaria/plasm-10-ch3-lambda1-1.0-lambda2-0-power-2.0.txt"
read_indent(file_name)

exit()
print("-----------------------------------------------")	
for charge in [2,3,4,5,6,7]:
	file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(charge)
	read_indent(file_name)
file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(2345)
read_indent(file_name)

print("-----------------------------------------------")	
for charge in [2,3,4,5]:
	file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch%d_test.txt.sumScoreIdent"%(charge)
	read_indent(file_name)
file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_010211_ch%d_test.txt.sumScoreIdent"%(2345)
read_indent(file_name)

print("-----------------------------------------------")	
for charge in [2,3,4]:
	file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng_080510_HapMap9_1_ch%s_test.txt.sumScoreIdent"%(charge)
	read_indent(file_name)

#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_Linfeng_080510_HapMap9_1_ch%s_test.txt.sumScoreIdent"%(2345)
#read_indent(file_name)


file_name="/n/trombone/s1/wrbai/codes/run_crux/Linfeng_9_1-tide-Pvalue/Linfeng_tide_pvalue.txt"

read_indent(file_name)

file_name="/n/trombone/s1/wrbai/codes/run_crux/Linfeng_32_1-tide-Pvalue/Linfeng_32_1-ch2/Linfeng_tide_pvalue.txt"

read_indent(file_name)

file_name="/n/trombone/s1/wrbai/codes/run_crux/Linfeng_32_1-tide-Pvalue/Linfeng_32_1-ch3/Linfeng_tide_pvalue.txt"

read_indent(file_name)


file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_yeast_01_ch2.txt.sumScoreIdent"

read_indent(file_name)



#print("-----------------------------------------------")	
#for charge in [2,3,4,5]:
#	file_name="/n/trombone/s1/wrbai/codes/backup05272014/286/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(charge)
#	read_indent(file_name)


#print("-----------------------------------------------")	
#for charge in [2,3]:
#	file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_bullseyeworm-ch%d.txt.sumScoreIdent"%(charge)
#	read_indent(file_name)
#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(2345)
#read_indent(file_name)

#print("-----------------------------------------------")	
#for charge in [2,3,4,5,6,7]:
#	file_name="/n/trombone/s1/wrbai/codes/msgfplus/malaria/malaria-ch%d/malaria-ch%d-ident-20150821.txt"%(charge,charge)
#	read_indent(file_name)
#file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(2345)
#read_indent(file_name)
#print("-----------------------------------------------")	
#for charge in [2,3,4,5,6,7]:
#for charge in [2]:
#	file_name="/n/trombone/s1/wrbai/codes/msgfplus/malaria_decoy3/malaria-ch%d/malaria-ch%d-ident-20150821.txt"%(charge,charge)
#	read_indent(file_name)




#print("-------------------------------------")
#for charge in [2,3,4,5]:
#	file_name="/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch%d/malaria-ch%d-ident.txt"%(charge,charge)
#	read_indent(file_name)
#file_name="/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch2345-correct-ident.txt"
#
#read_indent(file_name)	
#file_name="/s1/wrbai/codes/merge_charge_ident/malaria-ch2345-ident.txt"
#read_indent(file_name)	


#for charge in [2,3,4,5,6,7]:
#	file_name="/n/trombone/s1/wrbai/codes/msgfplus/malaria_jeffH/malaria-ch%d/malaria-ch%d-ident.txt"%(charge,charge)
#	read_indent(file_name)
#print("-----------------------------------------------")	
#print("-----------------------------------------------")	
#for charge in [2,3,4,5,6,7]:
#	file_name="/n/trombone/s1/wrbai/codes/msgfplus/malaria_crux9_17/malaria-ch%d/malaria-ch%d-ident-20150821.txt"%(charge,charge)
#	read_indent(file_name)
#print("-----------------------------------------------")	

#for charge in [2,3,4,5,6,7]:
#	file_name="/n/trombone/s1/wrbai/codes/msgfplus/malaria_crux/malaria-ch%d/malaria-ch%d-ident-20150821.txt"%(charge,charge)
#	read_indent(file_name)




#file_name="/s1/wrbai/codes/msgfplus/malaria_jeffH/malaria-ch2345-ident.txt"
#read_indent(file_name)


#file_name="/s1/wrbai/codes/msgfplus/jeffH_malaria_ident.txt"
#read_indent(file_name)


#file_name="/s1/wrbai/codes/msgfplus/malaria_decoy1/malaria-ch2/malaria-ch2-ident-20150821.txt"
#read_indent(file_name)

#for charge in [2,3,4,5]:
#	file_name="/n/trombone/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(charge)
#	read_indent(file_name)
#file_name="/n/trombone/s1/wrbai/codes/output/calibration/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(2345)
#read_indent(file_name)
#
#file_name="/n/trombone/s1/wrbai/codes/backup05272014/293/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
#
#read_indent(file_name)
#
#
#print("-----------------------------------------------")	
#for charge in [2,3,4,5]:
#	file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch%d.txt.sumScoreIdent"%(charge)
#
#	read_indent(file_name)
#print("-----------------------------------------------")	
#
##for charge in [3]:
##	file_name="/n/trombone/s1/wrbai/codes/msgfplus/linfeng_010211/linfeng-ch3/linfeng-ch3-ident.txt"%(charge,charge)
##	read_indent(file_name)
#
#file_name="run_crux/Linfeng_112710_tide_pvalue.txt"
#
#read_indent(file_name)
#
file_name="result.txt"
read_indent(file_name)
exit()
file_name="run_crux/Linfeng010221-output-1-7-tide/Linfeng_tide.txt"
#
read_indent(file_name)
file_name="run_crux/Linfeng010221-output-1-7-tide-Pvalue/Linfeng_tide_pvalue.txt"
#
read_indent(file_name)
#
#
#
#print("-----------------------------------------------")	
#for charge in [2,3,4,5]:
#	file_name="/n/trombone/s1/wrbai/codes/output/test_scores_output_fast_linfeng_112710_ch%d_test.txt.sumScoreIdent"%(charge)
#
#	read_indent(file_name)
#print("-----------------------------------------------")	
#
#for charge in [2,3,4,5]:
#	file_name="output/test_scores_output_fast_malaria-TMT12_ch%d_test.txt.sumScoreIdent"%charge
#	read_indent(file_name)
#
#file_name="run_crux/malaria_TMT11_tide_pvalue.txt"
#
#read_indent(file_name)
#
#
#print("-----------------------------------------------")	
#for charge in [2,3,4,5,6,7]:
#	file_name="/n/trombone/s1/wrbai/codes/msgfplus/malaria_TMT11/malaria-ch%d/malaria-ch%d-ident.txt"%(charge,charge)
#	read_indent(file_name)
#
#
#file_name="msgfplus/malaria_TMT10/malaria-ch2345-ident.txt"
#
#read_indent(file_name)
#
#file_name="msgfplus/malaria_TMT11/malaria-ch2345-ident.txt"
#
#read_indent(file_name)
#
#file_name="run_crux/malaria-TMT10_output_exact-p/malaria_TMT10_tide_pvalue_T.txt"
#
#read_indent(file_name)
#
#file_name="run_crux/malaria-TMT12_output_tide/malaria_TMT12_tide.txt"
#read_indent(file_name)
#
#file_name="run_crux/malaria-TMT12_output_exact-pvalue//malaria_TMT12_tide_pvalue_T.txt"
#read_indent(file_name)

file_name="output/test_scores_output_fast_yeast_01_ch3.txt.sumScoreIdent"
read_indent(file_name)
file_name="backup05272014/305/test_scores_output_fast_yeast_01_ch2.txt.sumScoreIdent"
file_name="backup05272014/305/test_scores_output_fast_yeast_01_ch3.txt.sumScoreIdent"
file_name="run_crux/bowyer_output-tide/bowyer_tide.txt"
read_indent(file_name)
file_name="run_crux/bowyer_output_exact-p/bowyer_tide_pvalue.txt"
read_indent(file_name)
file_name="msgfplus/linfeng_010211-2/linfeng-ch2345-ident.txt"

read_indent(file_name)
