import string
import sys
import random
import os
import math


if __name__ == "__main__":
     random.seed()
     infileName = ""
     outfileName = ""
     rsize = 0
     rnumber = 0
     shuffleperseq = 1
     if len(sys.argv) < 2:
         print "\nUsage: <euler output file> <converted procrastAlign file>"
         sys.exit(1)
     else:
         infileName = sys.argv[1]
         outfileName = sys.argv[2]
         
     eulerout = open(infileName,'r')
     procrastout = open(outfileName,'w')
     lines = eulerout.readlines()
     lmas = {}
     lma_count = 0
     key = ""
     alignments = {}
     inaln = 0
     aln = ""
     alns = []
     for line in lines:
         if line.find("#") != -1:
             lma_count +=1
             if key != "":
                 alignments[key] = alns[::]
             alns = []
             aln = ""
             #print lma_count
             if lma_count < 10:
                 key = "#procrastAlignment 0%s"%lma_count
             else:
                 key = "#procrastAlignment %s"%lma_count
             lmas[key] = []
             alignments[key] = []
             
         elif line.find(">") != -1:
             if len(aln) > 0:
                 alns.append(aln)
             aln = ""
             sp = line.find(",")
             ep = line.find("\n")
             spep = line[sp+1:ep]
             guion = spep.find("-")
             spos = int(spep[:guion])
             epos = int(spep[guion+1:])
             lmas[key].append([spos, epos])
             #key = "rubbish"
         else:
             inaln = 1
             ep = line.find("\n")
             aln += line
             aln = aln.replace('\n','')
             aln = aln.replace(' ','')
    
     
     #pick up leftovers
     if len(alns) > 0 and key != "":
         lma_count +=1
         alignments[key] = alns[::]
     
     sorted_keys = lmas.keys()
     sorted_keys.sort()
     xmfafile = open(outfileName+".xmfa",'w')
     xmfafile.write("#FormatVersion Mauve1\n#Sequence1File    modified.fna\n#Sequence1Entry    1\n#Sequence1Format    FastA\n")
     for key in sorted_keys:
         startpos = lmas[key]
         ccnt = 0
         width = 80
         tuples = {}
         for spos in startpos:
             #print spos[0]
             pos = 0
             xmfafile.write("> %d:%d-%d + modified_%d.fna\n"%(ccnt+1,int(spos[0]),int(spos[1]),ccnt+1))
             while pos+width < len(alignments[key][ccnt-1]):
                 xmfafile.write(alignments[key][ccnt-1][pos:pos+width])
                 xmfafile.write("\n")
                 pos+= width
             xmfafile.write(alignments[key][ccnt-1][pos:])
             xmfafile.write("\n")
             ccnt +=1
         xmfafile.write("=\n")
     xmfafile.close()
     for key in sorted_keys:
         if len(lmas[key]) == 0:
             continue
         procrastout.write(key+"\n")
         spos_list = []
         len_list = []
         values = lmas[key]
         for value in values:
             spos_list.append(value[0])
             len_list.append(int(value[1])-int(value[0]))
         procrastout.write(`max(len_list)`)
         procrastout.write("\t")
         for value in spos_list:
             procrastout.write(`value`)
             procrastout.write("\t")
         procrastout.write("\nlens: ")
         for value in len_list:
             procrastout.write(`value`)
             procrastout.write("\t")
         procrastout.write("\n")
         
     print "finished converting %d lmas from euler to procrast format"%len(lmas.keys())
     procrastout.close()