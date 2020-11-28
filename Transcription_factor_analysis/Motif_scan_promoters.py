# -*- coding: utf-8 -*-
#Python script to extract motifs from promoter sequences
import re
#Input check path
infile=open("promoters.txt","r")
#Output
out=open("TF_targets.txt","w")
pattern1=re.compile("TFbindingmotif",flags=re.IGNORECASE)

first = 0;
for line in infile:
   line = line.strip("\n")
   if line.startswith('>'):
       if first == 1:
          print '%s:%s' %(name,count1)
          out.write ('%s:\t%s\n' %(name,count1))
         # out.write ('%s\n' %(count1))
       count1 = 0
       #count2 = 0
       name=line
   else:
      first = 1
      s = re.findall(pattern1,line)
      count1 = count1 + len(s);
      #s = re.findall(pattern2,line)
      #count2 = count2 + len(s);
print '%s:%s' %(name,len(s))
out.write ('%s:\t%s\n' %(name,len(s)))
