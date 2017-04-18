import csv
#target_file="Linfeng_112710_output/tide-search.target.txt"
import os.path
for j in range(20):
	id =j+1
	print id
	for i in range(4):
		charge=i+2
	#	target_file="linfeng_9_1/"
	#	target_file="linfeng_9_%d-1/linfeng-ch%d/linfeng-ch%d-targets.tsv"%(id,charge,charge)
#		target_file="linfeng-%d/linfeng-ch%d/linfeng-ch%d-targets.tsv"%(id,charge,charge)
		#target_file="kim-%d/kim-ch%d/kim-ch%d-targets.tsv"%(id,charge,charge)
		target_file="plasm-%d/plasm-ch%d/plasm-ch%d-targets.tsv"%(id,charge,charge)
		print target_file
		#target_file="Linfeng_010211_output_2/tide-search.target.txt"
		#decoy_file="Linfeng_112710_output/tide-search.decoy.txt"
#		decoy_file="linfeng_9_%d-1/linfeng-ch%d/linfeng-ch%d-decoys.tsv"%(id,charge,charge)
#		decoy_file="linfeng-%d/linfeng-ch%d/linfeng-ch%d-decoys.tsv"%(id,charge,charge)
		decoy_file="plasm-%d/plasm-ch%d/plasm-ch%d-decoys.tsv"%(id,charge,charge)
		#decoy_file="kim-%d/kim-ch%d/kim-ch%d-decoys.tsv"%(id,charge,charge)
	
#		output=open("linfeng-%d/linfenlinfeng.txt"%(id,charge),'w')		
		output=open("plasm-%d/plasm-ch%d/plasm.txt"%(id,charge),'w')		
#		output=open("kim-%d/kim-ch%d/kim.txt"%(id,charge),'w')		

		if not os.path.isfile(target_file):
			continue
		if not os.path.isfile(decoy_file):
			continue
		#decoy_file="Linfeng_010211_output_2/tide-search.decoy.txt"
		target=[]
		decoy=[]
		name="EValue"
	#	name="MSGFScore"
		set_t=[]
		with open(target_file) as csvfile:
			spamreader = csv.DictReader(csvfile, delimiter='\t')
			for row in spamreader:
				s=row["Peptide"]
				a=s.split("+229.163")
				b=''.join(a)
				c=b.split("+57.021")
				d=''.join(c)
				a1=row["ScanNum"]
				a2=row["Title"]
				a3=a2.split(".")
				a4=a3[-1]
				sid=int(a4)
			#	a1=row["SpecID"]
			#	a2=a1.split("index=")	
			#	a3=''.join(a2)	
				if sid not in set(set_t):
					set_t.append(sid)
					target.append((sid,d,float(row[name])))
		set_d=[]		
		with open(decoy_file) as csvfile:
			spamreader = csv.DictReader(csvfile, delimiter='\t')
			for row in spamreader:
				s=row["Peptide"]
				a=s.split("+229.163")
				b=''.join(a)
				c=b.split("+57.021")
				d=''.join(c)
				a1=row["ScanNum"]
				a2=row["Title"]
				a3=a2.split(".")
				a4=a3[-1]
				sid=int(a4)
			#	a1=row["SpecID"]
			#	a2=a1.split("index=")	
			#	a3=''.join(a2)	
				if sid not in set(set_d):
					set_d.append(sid)
					decoy.append((sid,d,float(row[name])))
		
		
		#		decoy.append((row["scan"],d,float(row[name])))
		#output=open("Linfeng_112710_tide_pvalue.txt",'w')		
		#output=open("Linfeng_010211_tide_pvalue_2.txt",'w')		
		idt=[x[0] for x in target]
		idd=[x[0] for x in decoy]
		id_t=list(set(idd)-set(idt))
		id_d=list(set(idt)-set(idd))
		print len(idt),len(idd)
		for i in id_t:
			target.append((i,"AAAAAA",10000000))
		for i in id_d:
	
			decoy.append((i,"AAAAAA",10000000))
		target=sorted(target,key=lambda x:x[0])	
		decoy=sorted(decoy,key=lambda x:x[0])	
		set1=set([])
		target2=[]
		
		for row in target:
			if row[0] not in set1:
				target2.append(row)
				set1.add(row[0])
		set2=set([])	
		decoy2=[]
		for row in decoy:
			if row[0] not in set2:
				decoy2.append(row)
				set2.add(row[0])
		target=target2		
		decoy=decoy2
		print "asdf",len(target),len(decoy)
		tt=[x[2]-y[2]<0 for x,y in zip(target,decoy)]
		print sum(tt),len(tt),float(sum(tt))/float(len(tt))
	
		output.write("Kind\tSid\tPeptide\tScore\n")
		for row1,row2 in zip(target,decoy):
			a=-row1[2]
			b=-row2[2]
			if a>b:
				b=-10000000
			else:
				a=-10000000

			output.write("t\t%s\t%s\t%.15f\n"%(row1[0],row1[1],a))
			output.write("d\t%s\t%s\t%.15f\n"%(row2[0],row2[1],b))
