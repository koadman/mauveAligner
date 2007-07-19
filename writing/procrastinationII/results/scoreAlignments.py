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
         print "\nUsage: scoreAlignments <correct alignment(tuples)> <calculated alignment(xmfa)> "
         sys.exit(1)
     else:
         correct = sys.argv[1]
         calculated = sys.argv[2]
         
     correct_tuples = {}
     correct_file = open(correct, 'r')
     lines = correct_file.readlines()
     correct_num = lines[0]
     correct_num = correct_num.replace("\n","")
     correct_num = int(correct_num)
     tcount = 0
     #parse correct alignment tuples
     for line in lines[1:]:
         line = line.replace("\n","")
         data = line.split('\t')
         ckey = "%d:%d"%(int(data[0]),int(data[1]))
         correct_tuples[ckey] = 1
         tcount +=1
     if tcount == correct_num:
         print "%d correct tuples read in OK! "%tcount 
     #read through calculated xmfa to get list of tuples
     calculated_file = open(calculated,'r')
     lines = calculated_file.readlines()
     lines = lines[4:]
     calculated_file.close()
     calculated_file = open(calculated,'r')
     alldata = calculated_file.read()
     num_lmas = alldata.count("=")
     lmacount = 0
     endlma = 0
     calc_tuples = []
     calc_dict = {}
     startpos = []
     alignment = []
     curI = 0
     numlines = len(lines)
     cI = 0
     #parse xmfa file
     while curI < len(lines):
         pI = curI-1
         if (curI*100) / numlines != (pI*100)/ numlines:
             print `(curI*100)/numlines`+"%...",
         if lines[curI].find('>') != -1 and endlma != 1:
             sI = lines[curI].find(':')+1
             eI = lines[curI].find('+')-1
             data = lines[curI][sI:eI]
             data = data.split('-')
             startpos.append(int(data[0]))
             length = int(data[1])-int(data[0])
             curI +=1
         elif lines[curI].find('>') == -1 and lines[curI].find('=') == -1:
             #this would be the alignment data
             data = ""
             while lines[curI].find('>') == -1:
                 if lines[curI].find('=') != -1:
                     break 
                 data += lines[curI].replace('\n','')
                 curI+=1
             
             alignment.append(data)
             #print data
         elif lines[curI].find('=') != -1:
             #end of lma
             endlma = 1
             ccnt = 1
             cI= 0
             #create tuple data from alignment
             for spos in startpos:
                 combinations = (len(startpos)*(len(startpos)-1))/2 
                 #print spos
                 pos = 0
                 for tpos in startpos[ccnt:]:
                     for site in xrange(len(alignment[ccnt-1])):
                         #print alignment[ccnt-1][site],alignment[ccnt][site]
                         if alignment[ccnt-1][site] != "-" and alignment[ccnt][site] != "-":
                             ttup = [spos+site,tpos+site]
                             ckey = "%d:%d"%(spos+site,tpos+site)
                             calc_dict[ckey] = 1
                             ckey = "%d:%d"%(tpos+site,spos+site)
                             calc_dict[ckey] = 1
                 ccnt +=1
             curI +=1
             endlma = 0
             startpos = []
             alignment = []
     
     truepos = 0
     falseneg = 0
     correct_keys = correct_tuples.keys()
     total_correct = len(correct_keys)
     #total_calc = len(calc_tuples)
     calc_keys = calc_dict.keys()
     total_calc = len(calc_keys)
     
     for key in correct_keys:
         #print key
         if calc_dict.has_key(key):
             truepos +=1
         else:
             falseneg +=1
     
     print "true positives: ", truepos
     print "possible: ", total_correct
     print "sensitivity: ", float(truepos)/float(total_correct)
     print "specificity: ", float(truepos)/float(total_calc)
     
      