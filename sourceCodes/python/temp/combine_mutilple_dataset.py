import csv
import os.path
from math import *

def read_indent(FILE_NAME):
	input_file=csv.DictReader(open(FILE_NAME),delimiter='\t')
	print(FILE_NAME)
	target_decoy=[]
	sid=[]
	for row in input_file:
		target_decoy.append((row["Kind"],int(row["Sid"]),row["Peptide"],float(row["Score"])))	
		sid.append(int(row["Sid"]))
	spec2sorted=sorted(target_decoy,key=lambda x:x[3],reverse=True)
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
	


	return target_decoy, max(sid),score2
def combine_mutilple_dataset(file_name,output_name):
	print(output_name)
	s_max=0;
	t_all=[]
	for name in file_name:
		t,s,sc=read_indent(name)
		if s>40000:
			print s
		s=40000;
		sc=1
		for row in t:
			t_all.append((row[0],row[1]+s_max,row[2],(row[3]/(sc))))
		s_max=s_max+s	
	t_t=[x[1] for x in t_all if x[0]=='t']
	
	t_d=[x[1] for x in t_all if x[0]=='d']
	
	print len(t_t)-len(set(t_t))
	print len(t_d)-len(set(t_d))
	sid_no_t=set(t_d)-set(t_t)
	sid_no_d=set(t_t)-set(t_d)
	print(sid_no_t)
	print(sid_no_d)
	for i in sid_no_t:
		t_all.append(('t',i,"AAAAAA",-100000))
	
	for i in sid_no_d:
		t_all.append(('d',i,"AAAAAA",-100000))
	
	t_all=sorted(t_all,key=lambda x:x[0],reverse=True)
	t_all=sorted(t_all,key=lambda x:x[1])
	
	
	f=open(output_name,'w')
	f.write('Kind\tSid\tPeptide\tScore\n')
	for row in t_all:
		f.write('%s\t%s\t%s\t%.10f\n'%(row[0],row[1],row[2],row[3]))
	

org_id=2
method_id=1
orgs=["yeast","worm","malaria","human","human4","kim3","linfeng-lambda1-10-lambda2-10"]
org=orgs[org_id]
methods=["bipartite","msgfPlus","pvalue","tide"]
method=methods[method_id]
file_name=[]
#for id in range(3):
#	temp="merge_charge_ident/msgfPlusIdents/msgfPlus-%s-0%d-idents/msgfPlus-%s-0%d-ch23-ident.txt"%(org,id+1,org,id+1)
#	t="tide"
#	temp="/s1/wrbai/codes/john_idents/%sIdents/%s-%s-0%d-idents/%s-%s-0%d-ch123-ident.txt"%(t,t,org,id+1,t,org,id+1)
#	temp="/s1/wrbai/codes/john_idents/%sIdents/%s-%s-0%d-idents/%s-%s-0%d-ch123-ident.txt"%(t,t,org,id+1,t,org,id+1)
#	#file_name.append("/s1/wrbai/codes/merge_charge_ident/1/test_scores_output_fast_%s_0%d_ch23.txt.sumScoreIdent"%(org,id+1))
#	file_name.append(temp)
#TMT=[]
#TMT.append(["/s1/wrbai/codes/backup05272014/292/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent","/s1/wrbai/codes/msgfplus/jeffH_malaria_ident.txt" ,"/s1/wrbai/codes/run_crux/malaria-TMT10_output_exact-p/malaria_TMT10_tide_pvalue_T.txt","/s1/wrbai/codes/run_crux/malaria-TMT10_output-tide/malaria_TMT10_tide_pvalue_T.txt"])
#TMT.append(["/s1/wrbai/codes/backup05272014/300/test_scores_output_fast_malaria-TMT11_ch2345_test.txt.sumScoreIdent","/s1/wrbai/codes/msgfplus/malaria_TMT11/malaria-ch2345-ident.txt","/s1/wrbai/codes/run_crux/malaria_TMT11_tide_pvalue_T.txt","/s1/wrbai/codes/run_crux/malaria_TMT11_tide_pvalue.txt"])
#TMT.append(["/s1/wrbai/codes/backup05272014/301/test_scores_output_fast_malaria-TMT12_ch2345_test.txt.sumScoreIdent","/s1/wrbai/codes/msgfplus/malaria_TMT12/malaria-ch2345-ident.txt","/s1/wrbai/codes/run_crux/malaria-TMT12_output_exact-pvalue/malaria_TMT12_tide_pvalue_T.txt","/s1/wrbai/codes/run_crux/malaria-TMT12_output_tide/malaria_TMT12_tide.txt"])
#TMT.append(["/s1/wrbai/codes/backup05272014/301/test_scores_output_fast_malaria-TMT13_ch2345_test.txt.sumScoreIdent","/s1/wrbai/codes/msgfplus/malaria_TMT13/malaria-ch2345-ident.txt","/s1/wrbai/codes/run_crux/malaria-TMT13_output_exact-pvalue/malaria_TMT13_tide_pvalue_T.txt","/s1/wrbai/codes/run_crux/malaria-TMT13_output_tide/malaria_TMT13_tide_pvalue_T.txt"])
#for id in range(4):
#	file_name.append(TMT[id][method_id])
#for i in range(6):
#	id=i+1
	#file_name.append("/n/trombone/s1/wrbai/codes/backup05272014/310/test_scores_output_fast_Linfeng-%d-ch2345_test.txt.sumScoreIdent"%(id))
#	file_name.append("/n/trombone/s1/wrbai/codes/backup05272014/309/test_scores_output_fast_Linfeng_080510_HapMap9_%d_ch2345_test.txt.sumScoreIdent"%(id))
#for i in range(100):
#	id=i+1
	#file_name.append("/n/trombone/s1/wrbai/codes/backup05272014/310/test_scores_output_fast_Linfeng-%d-ch2345_test.txt.sumScoreIdent"%(id))
#	file_name.append("/n/trombone/s1/wrbai/codes/msgfplus/linfeng-%d/linfeng.txt"%(id))
	#file_name.append("/n/trombone/s1/wrbai/codes/msgfplus/kim-%d/kim.txt"%(id))
	#file_name.append("/n/trombone/s1/wrbai/codes/output/kim-random20-%d.txt"%id)
	#file_name.append("/n/trombone/s1/wrbai/codes/msgfplus/linfeng_9_%d-1/linfeng.txt"%(id))
for j in range(1):
	charge = j+2;
	file_name=[]
	for i in range(20):
		id = i+1

		#file_name.append("/s1/wrbai/codes/run_crux/plasm/plasm-%d-pvalue/plasm_tide_pvalue.txt"%id)
#		file_name.append("/s1/wrbai/codes/backup05272014/329/plasm/plasm-%d-ch%d.txt"%(id,charge))
#		aa1=("/n/trombone/s1/wrbai/codes/output/postbcb-mean/plasm-%d-lambda1-7-lambda2-15.txt")%id
	#	aa1=("/n/trombone/s1/wrbai/codes/output/postbcb/2016_6_29/Linfeng1/Linfeng-%d-lambda1-10-lambda2-10.txt")%id
		aa1=("../../MSGFPlus/plasm-%d/plasm.txt"%(id))
		#aa1=("/s1/wrbai/codes/run_crux/Linfeng-%d-tide/Linfeng_tide.txt"%(id))
#		aa1=("/s1/wrbai/codes/backup05272014/330/plasm/plasm-uncal-%d-ch%d.txt"%(id,charge))
#		aa1=("/s1/wrbai/codes/backup05272014/329/plasm/plasm-%d-ch%d.txt"%(id,charge))
		if os.path.isfile(aa1):
			num_lines = sum(1 for line in open(aa1))
			if num_lines>20:
				file_name.append(aa1)
	#file_name.append("/s1/wrbai/codes/msgfplus/plasm-%d/plasm.txt"%id)
#	file_name.append("/s1/wrbai/codes/run_crux/plasm/plasm-%d-tide/plasm_tide.txt"%id)
#for i in range(100):
#	id=i+1
#	file_name.append("/s1/wrbai/codes/run_crux/Linfeng-%d-tide-pvalue/Linfeng_tide_pvalue.txt"%(id))
#	file_name.append("/s1/wrbai/codes/run_crux/Linfeng-%d-tide/Linfeng_tide.txt"%(id))

	output_name="%s-%s.txt"%(org,method_id)	
	combine_mutilple_dataset(file_name,output_name)
