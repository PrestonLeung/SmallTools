#!/usr/local/bin/python
#Ver. 0.041
#Pre.L

import argparse, re
import os, sys, csv
import scipy.stats as sci_stat
import numpy as np

#========================"Global Variables"=======================#

OutputDir = ''

#========================"Self defined Errors and Exceptions"=======================#

class Error(Exception):
    """Base class for exceptions."""
    pass
    
#Error classes for handling input errors
class InputError (Error):
    """Exception for errors in the input.
        Attributes:
        expr -- the input that caused the error
        msg  -- message giving details of the error
    """
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg    
    def __str__(self):
        return repr(str(self.msg)+str(self.expr))

#==========================="okFile"===========================#
#Check if the given file is a file.

def okFile(file):
    if(os.path.isfile(file)):
        return True
    else:
        return False
        
#==========================="Change Percentage to Float numbers"===========================#
# p is a string which needs to be converted to a float.
# sometimes they occur over a range and is represented by
# something like 0.07% -> 0.09%. We take the lower bound.
# This is chosen to reduce the chance of distribution of SNPs
# summing to > 1.

def stringToFloat(p):

    hasArrow = re.search('(\d+\.\d+\%) \-\> (\d+\.\d+\%)',p)
    if(hasArrow):
        result = float(hasArrow.group(1).strip('%'))/100
        #b = float(hasArrow.group(2).strip('%'))/100
        #result = np.average([a,b])
    else:
        result = float(p.strip('%'))/100

    return result

#==========================="Parse Geneious CSV"===========================#
# Geneious CSV example:
# Column01: NT Change           Column09: Variant Frequency
# Column02: Start               Column10: Variant P-Value (Approx.)
# Column03: End                 Column11: Amino Acid Change
# Column04: Length              Column12: CDS Codon Number
# Column05: Avg Quality         Column13: CDS Position
# Column06: Change From->To     Column14: Codon Change
# Column07: Coverage            Column15: Protein Effect
# Column08: Polymorphism Type

# s_column => NT change Start Column
# pt_column => Polymorphism Type column number
# pe_column => Protein Effect column number 
# Global (Non-Syn & Syn)

# This function parses the CSV file and then manipulates it into a dictionary
# Handles mutations by codons.

def parseCSV(filename,pt_column,pe_column, s_column, cc_column, aa_column,fq_column):
    headerCount = 1
    snpHash = {}

    with open(filename, 'rb') as csvfile:
        csvHandle = csv.reader(csvfile)
        for row in csvHandle:
            if(headerCount > 0):
                headerCount -= 1
                continue
            
            if(emptyArrowCodon(row[cc_column])):
                continue
            
            if(re.match('Substitution', row[pe_column]) or re.match('None', row[pe_column])):
                if(re.match('Substitution', row[pt_column])):
                    splitHash = splitCodon(row,pt_column,pe_column, s_column, cc_column, aa_column,fq_column)
                    #print "New SplitHash"
                    #for key in sorted(splitHash.keys()):
                    #    print "Key = {}.".format(key), splitHash[key]

                    for key in sorted(splitHash.keys()):
                        #print "Key = {}.".format(key), splitHash[key]
                        if(key in snpHash):
                            #print "Check Point 1.1"
                            tempArray = []
                            foundFlag = False
                            for item in snpHash[key]:

                                if(item[cc_column] == splitHash[key][cc_column]):
                                    #print "Check Point 1.1.1"
                                    item[fq_column] = str(float(re.sub('\%','',item[fq_column])) +
                                                          float(re.sub('\%','',splitHash[key][fq_column]))) + "%"
                                    foundFlag = True
                                    break
                                else:
                                    #print "Check Point 1.1.2"
                                    if(splitHash[key] not in tempArray):
                                        tempArray.append(splitHash[key])

                            if(foundFlag and len(tempArray) > 0):
                                tempArray.remove(splitHash[key])

                            if(len(tempArray)>0):
                                #print "tempArray: ", tempArray
                                for item in tempArray:
                                    snpHash[key].append(item)

                        else:
                            #print "Check Point 1.2"
                            snpHash[key] = [splitHash[key]]
                else:
                    #print "Check Point 1.3 - {}".format(row)
                    stringNum = int(re.sub(',','',row[s_column]))
                    codonNum = posToCodon(stringNum)

                    if(codonNum in snpHash):
                        snpHash[codonNum].append(row)
                    else:
                        snpHash[codonNum] = [row]

	return snpHash


#==========================="splitting very annoying codons from geneious"===========================#

def splitCodon(row,pt_column,pe_column, s_column, cc_column, aa_column,fq_column):

    splitHash = {}

    if(re.match('Substitution', row[pt_column])):
        wtCodonList = row[cc_column].split(' -> ')[0].split(',')
        mtCodonList = row[cc_column].split(' -> ')[1].split(',')

        assert(len(wtCodonList) == len(mtCodonList))
        aaCount = len(wtCodonList)

        if(aaCount > 1): #need to split the codons
            aaCountCP = 0
            #find the first codon
            codonNum = posToCodon(int(re.sub(',','',row[s_column])))

            while(aaCountCP < aaCount):
                codonNum2 = codonNum + aaCountCP
                tempRow = list(row)
                if(re.match('None', row[pe_column])):
                    tempRow[aa_column] = ""
                elif(re.match('Substitution', row[pe_column])):
                    wtCodonList2 = list(row[aa_column].split(' -> ')[0])
                    mtCodonList2 = list(row[aa_column].split(' -> ')[1])

                    if(wtCodonList2[aaCountCP] == mtCodonList2[aaCountCP]):
                        tempRow[aa_column] = ""
                        tempRow[pe_column] = "None"
                    else:
                        tempRow[aa_column] = "{} -> {}".format(wtCodonList2[aaCountCP],mtCodonList2[aaCountCP])
                        tempRow[pe_column] = "Substitution"

                tempRow[s_column] = "{}".format(codonNum2)
                tempRow[cc_column] = "{} -> {}".format(wtCodonList[aaCountCP],mtCodonList[aaCountCP])
                tempRow[pt_column] = "SNP (manual)"
                hasArrow = re.search('(\d+\.\d+\%) \-\> (\d+\.\d+\%)',tempRow[fq_column])
                
                if(hasArrow):
                    a = float(hasArrow.group(1).strip('%'))
                    b = float(hasArrow.group(2).strip('%'))
                    tempRow[fq_column] = str(min([a,b]))+"%"

                splitHash[codonNum2] = tempRow

                aaCountCP +=1

        else: #no need splitting
            stringNum = int(re.sub(',','',row[s_column]))
            codonNum = posToCodon(stringNum)
            splitHash[codonNum] = row

    return splitHash



#==========================="Convert nucleotide position to codon number"===========================#
# This is done because when we created the csv files from geneious not
# all of them had CDS codon number. So have to manually convert.

def posToCodon(ntNum):
    r = ntNum%3

    if(r == 0):
         codonNum = int(ntNum/3)
    else:
         codonNum = int(ntNum/3) + 1

    return codonNum


#==========================="Parse Geneious CSV nonSyn"===========================#
# s_column => NT change Start Column
# pt_column => Polymorphism Type column number
# pe_column => Protein Effect column number 
# Non-Syn mutations only

def parseCSV_nonSyn(filename,pt_column,pe_column, s_column,cc_column,aa_column,fq_column):
    headerCount = 1    
    snpHash = {}
    
    with open(filename, 'rb') as csvfile:
        csvHandle = csv.reader(csvfile)
        for row in csvHandle:
            if(headerCount > 0):
                headerCount -= 1
                continue

            if(emptyArrowCodon(row[cc_column])):
                continue

            if(re.match('Substitution', row[pe_column])):
                if(re.match('Substitution', row[pt_column])):

                    splitHash = splitCodon(row,pt_column,pe_column, s_column, cc_column, aa_column,fq_column)

                    for key in sorted(splitHash.keys()):
                        if(splitHash[key][pe_column] == "None"):
                            continue
                        if(key in snpHash):
                            tempArray = []
                            foundFlag = False
                            for item in snpHash[key]:
                                if(item[cc_column] == splitHash[key][cc_column]):
                                    item[fq_column] = str(float(re.sub('\%','',item[fq_column])) +
                                                          float(re.sub('\%','',splitHash[key][fq_column]))) + "%"
                                    foundFlag = True
                                    break
                                else:
                                    if(splitHash[key] not in tempArray):
                                        tempArray.append(splitHash[key])

                            if(foundFlag and len(tempArray) > 0):
                                tempArray.remove(splitHash[key])

                            if(len(tempArray)>0):
                                for item in tempArray:
                                    snpHash[key].append(item)

                        else:
                            snpHash[key] = [splitHash[key]]
                else:
                    stringNum = int(re.sub(',','',row[s_column]))
                    codonNum = posToCodon(stringNum)

                    if(codonNum in snpHash):
                        snpHash[codonNum].append(row)
                    else:
                        snpHash[codonNum] = [row]
	
    return snpHash

#==========================="Parse Geneious CSV Syn"===========================#
# s_column => NT change Start Column number
# pt_column => Polymorphism Type column number
# pe_column => Protein Effect column number 
# Syn mutations only

def parseCSV_syn(filename,pt_column,pe_column, s_column,cc_column,aa_column,fq_column):
    headerCount = 1    
    snpHash = {}
    
    with open(filename, 'rb') as csvfile:
        csvHandle = csv.reader(csvfile)        
        for row in csvHandle:
            if(headerCount > 0):
                headerCount -= 1
                continue

            if(emptyArrowCodon(row[cc_column])):
                continue

            if(re.match('Substitution', row[pe_column]) or re.match('None', row[pe_column])):
                if(re.match('Substitution', row[pt_column])):
                    splitHash = splitCodon(row,pt_column,pe_column, s_column, cc_column, aa_column,fq_column)
                    for key in sorted(splitHash.keys()):
                        if(splitHash[key][pe_column]=="Substitution"): #Only take in the entries with "None"
                            continue
                        if(key in snpHash):
                            tempArray = []
                            foundFlag = False
                            for item in snpHash[key]:
                                if(item[cc_column] == splitHash[key][cc_column]):
                                    item[fq_column] = str(float(re.sub('\%','',item[fq_column])) +
                                                          float(re.sub('\%','',splitHash[key][fq_column]))) + "%"
                                    foundFlag = True
                                    break
                                else:
                                    if(splitHash[key] not in tempArray):
                                        tempArray.append(splitHash[key])
                            
                            if(foundFlag and len(tempArray) > 0):
                                tempArray.remove(splitHash[key])
                            
                            if(len(tempArray)> 0):
                                for item in tempArray:
                                    snpHash[key].append(item)

                        else:
                            snpHash[key] = [splitHash[key]]
                elif(re.match('None', row[pe_column])):
                    stringNum = int(re.sub(',','',row[s_column]))
                    codonNum = posToCodon(stringNum)

                    if(codonNum in snpHash):
                        snpHash[codonNum].append(row)
                    else:
                        snpHash[codonNum] = [row]

	return snpHash

#==========================="Parse protein segment file"===========================#
# Breaks down the genome region into separate bits
# according to what was defined in the file
# ** needs to some how loop over the start stops to retreieve segments for the 
# ** confirmation (segements that are not continuous).

def readSegement(filename):

    chocolateHash = {}
    gLength = 0
    currMax = 0

    with open(filename, 'rb') as csvfile:
        csvHandle = csv.reader(csvfile)
        totalPos = 0
        for row in csvHandle: # this whole needs to be rewritten
            if(len(row)>=3):                
                protein = row[0]
                segment = row[1:len(row)]
                start = 0
                end = 0
                segLength = 0
                for i in range(1,len(row)):
                    if(i%2 == 1):
                        start = int(row[i])
                    elif(i%2 == 0):
                        end = int(row[i])
                        segLength += end - start + 1
                chocolateHash[tuple(segment)] = protein
                gLength += segLength
    return chocolateHash, gLength
	
#==========================="Calculate Shannon Entropy"===========================#
# Iterate through hash table, calculate SE for each position.
# Also output avg SE over given genome length.
# Note -> avgSE is determined by the regionHash. Even if SNP file
#         goes for full region, avgSE is the SE averaged over given
#         regions (ie. E1 E2 only if thats what's given).

def getSE(hashTable, f_col, gLength, regionHash,pt_col, cc_col):

    resultHash = {} #SE per position in hastTable
    regionSE = {}
    regionLength = {}
    avgSE_region = {} #avg SE per protein region
    std_dev = {} #standard deviation of SE

    #print "Check Point 2.1"
    #This loop is for getting SE of all positions (generally speaking)
    for key in sorted(hashTable.keys()):
        snpList = hashTable[key]
        #print "Check Point 2.1.1 - Key = {}, List = {}".format(key,snpList)
        freqList = appendFreq(snpList, f_col)
        resultHash[key] = sci_stat.entropy(freqList, base=2)

    #This loop splits SE into proteins regions
    #print "Check Point 2.2"
    for region in sorted(regionHash.keys()):
        
        # (start,end) = region #modify this
        protein = regionHash[region]
        posHash = {}
        for i in range(0,len(region)):            
            if(i%2 == 0):
                start = int(region[i])
            elif(i%2 == 1):
                end = int(region[i])
                if(protein in regionLength):
                    regionLength[protein] += (end - start + 1)
                else:
                    regionLength[protein] = (end - start + 1)        
                # this requires some mod to make it work
                for position in sorted(hashTable.keys()):
                    if(position > end):
                        break
                    if(position >= start and position <= end):
                        snpList = hashTable[position]
                        freqList = appendFreq(snpList, f_col)
                        if(position not in posHash):
                            posHash[position] = sci_stat.entropy(freqList, base=2)
                        else:
                            print "posHash should be not be pre-filled with position {}.".format(str(position))
        regionSE[protein] = posHash

    #Calculate avg SE over the entire region
    #print "Check Point 2.3"
    sumOfSE = 0
    for prot in regionSE:
        if(len(regionSE[prot].values())<1):
            sumOfSE += 0.0
        else:
            sumOfSE += sum(regionSE[prot].values())
    avgSE = sumOfSE / gLength

    #Calculating avg SE for each region
    #help chaturaka add here standard deviation of SE
    #print "Check Point 2.4"
    for prot in regionSE:
        #print "Check Point 2.4.1 - Protein:{}".format(prot)
        #print "Protein:{}, Protein Length: {}".format(prot,regionLength[prot])
        if(len(regionSE[prot].values())<1):            
            avgSE_region[prot] = 0
        else:            
            avgSE_region[prot] = sum(regionSE[prot].values())/regionLength[prot]
        #print "Check Point 2.4.2 - Protein:{}".format(prot)
        if(len(regionSE[prot].values())<1):
            std_dev[prot] = "N/A"
        else:
            #print "Check Point 2.4.3 - Protein:{}. Values  {}".format(prot,regionSE[prot].values())

            std_dev[prot] = np.std(regionSE[prot].values())
    #print "Check Point 2.Exit"
    return resultHash, avgSE, avgSE_region, std_dev

#==========================="appendFreq"===========================#

def appendFreq(theList, column):

    freqList = []

    for row in theList:
        freqList.append(stringToFloat(row[column]))

    if(1-sum(freqList) > 0):
        freqList.append(1-sum(freqList))

    return freqList

#==========================="print SE results"===========================#

def outputResults(hashTable, se_type):
    
    name = re.sub("/", "",OutputDir)
    
    with open(OutputDir+"ShannonEntropy_"+se_type+".txt", 'wb') as outputHandle:
        for key in sorted(hashTable.keys()):
            outputHandle.write("{}\t{}\t{}\n".format(name,key, hashTable[key]))
            
#==========================="main function"===========================#
# Checks if the codon change has proper codon change and not
# empty. Towards the end of the sequence, Geneious will still try
# to call SNPs, but it does this regardless if a codon is completed.
# e.g. of the 3 nucleotides, the 2nd or 3rd one is missing. What happens
# is that the SNP is called, but there's no information on the codon and
# Geneious ends up writing just an arrow, with  no before and after.
# If this happens, we skip that entry. This function checks if this case
# occurrs.

def emptyArrowCodon(codon):

    answer = False
    hasArrow = re.search(' -> ',codon)
    if(hasArrow):
        wt = codon.split(' -> ')[0]
        mt = codon.split(' -> ')[1]
        if(len(wt) < 1):
            answer = True
    return answer

#==========================="main function"===========================#

if __name__ == '__main__':  
    
    program_function = """
    *****   
        snpExtractorSEG.py is a SNP extractor tool that filters desired SNPS
        from Geneious (based on ver. R9.1) SNP files (usually in csv format) 
        and calculates Shannon Entropy (SE) by codons.

        Positions without recorded SNPs will be treated as 0.

        Note on Optional Arguments:

        -   By default this tool calculates SE with any nucleotide changes.
        -   Adding optional arguments -in and/or -is will produce extra files
            that filter SE calculated for non-synonymous and/or synonymous
            mutations separately.

        Note - ignores insertion/deletions/frameshifts/truncation
               in current version.

    	Ver. 0.05

	    Written by PresDawgz (~w~)v
    
    *****
    """
    parser = argparse.ArgumentParser()
    
    #Required Arguments
    parser.add_argument('-f', '--file', help='File to extract snp information from.',  required = True)
    parser.add_argument('-s', '--s_col',  help='The column position storing the positions of mutations.',  required = True,  type = int)
    parser.add_argument('-pt', '--pt_col',  help='The column position storing polymorphism type.',  required = True,  type = int)
    parser.add_argument('-pe', '--pe_col',  help='The column position storing protein effect.',  required = True,  type = int)
    parser.add_argument('-fq', '--freq_col',  help='The column position storing the frequencies.',  required = True,  type = int)
    parser.add_argument('-cc','--cdn_chg', help='The column position storing codon change.', required = True, type = int)
    parser.add_argument('-aa','--aa_chg', help='The column position storing amino acid change.', required = True, type = int)
    #parser.add_argument('-l', '--genome_len', help='Length of genome for calculating avg SE over all positions', required = True, type = int)
    parser.add_argument('-p', '--protein_segment', help='Protein region file for calculating avg SE over all positions')

    #Optional Arguments
    parser.add_argument('-is',  '--include_syn',  help='Include SE output for synonymous mutations only.',  action = 'store_true')
    parser.add_argument('-in',  '--include_nonsyn',  help='Include SE output for non-synonymous mutations only.',  action = 'store_true')        
    parser.add_argument('-d', '--directory', help='User specified directory to save results, otherwise saves at current location')
    parser.add_argument('-seo', '--SE_output', help="Retrieve the files containing SE calculations per position", action = 'store_true')
    
        
    if len(sys.argv) < 3:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()

    if(not okFile(args.file)):
        raise InputError(args.file,  'Cannot find file: ')
    
    if(not okFile(args.protein_segment)):
        raise InputError(args.protein_segment,  'Cannot find file: ')

    if(args.directory):
        if(os.path.exists(args.directory)):
            print "Using {} for output destination.".format(args.directory)
            pass
        else:
            print "Path: {} not found. Creating output destination.".format(args.directory)
            os.makedirs(args.directory)
        OutputDir = str(args.directory)+"/"

    # Read protein segment file
    pSegHash, genome_len = readSegement(args.protein_segment)

    # SE for Non-synonymous SNPs only.
    if(args.include_nonsyn):
        nonSynHash = parseCSV_nonSyn(args.file, (args.pt_col-1), (args.pe_col-1), (args.s_col-1),(args.cdn_chg-1),(args.aa_chg-1),(args.freq_col-1))
        #print "***NonSyn***"
        (result, avg_nonSyn_se, avg_nonSyn_se_Region, std_dv) = getSE(nonSynHash, (args.freq_col-1),genome_len, pSegHash,(args.pt_col-1),(args.cdn_chg-1))
        if(args.SE_output):
            outputResults(result, "Non_Synonymous")
        outputResults(avg_nonSyn_se_Region, "Non_Syn_Region")
        outputResults(std_dv, "Non_Syn_RegionSD")

    # SE for Synonymous SNPs only.
    if(args.include_syn):
        synHash = parseCSV_syn(args.file, (args.pt_col-1), (args.pe_col-1), (args.s_col-1),(args.cdn_chg-1),(args.aa_chg-1),(args.freq_col-1))
        #print "***Syn***"
        (result,avg_syn_se, avg_syn_se_Region,std_dv) = getSE(synHash, (args.freq_col-1), genome_len, pSegHash,(args.pt_col-1),(args.cdn_chg-1))
        if(args.SE_output):
            outputResults(result, "Synonymous")
        outputResults(avg_syn_se_Region, "Syn_Region")
        outputResults(std_dv, "Syn_RegionSD")

    # SE for All SNPs
    #print "Check Point 1."
    snpHash = parseCSV(args.file, (args.pt_col-1), (args.pe_col-1), (args.s_col-1),(args.cdn_chg-1),(args.aa_chg-1),(args.freq_col-1))

    #for key in sorted(snpHash.keys()):
    #    print "Key = {}.".format(key), snpHash[key]

    #print "Check Point 2."
    #print "***Global***"
    (result, avg_se, avg_se_Region,std_dv)= getSE(snpHash, (args.freq_col-1), genome_len, pSegHash,(args.pt_col-1),(args.cdn_chg-1))


    if(args.SE_output):
        outputResults(result, "Global")
    outputResults(avg_se_Region, "Global_Region")
    outputResults(std_dv, "Global_RegionSD")

    with open(OutputDir+"AverageSE.txt", 'wb') as outputHandle:
        outputHandle.write("Genome length (AA) = {}\n".format(genome_len))
        if(args.include_nonsyn):
            outputHandle.write("Average SE Non-synonymous SNPs only = {}\n".format(avg_nonSyn_se))
        if(args.include_syn):
            outputHandle.write("Average SE Synonymous SNPs only = {}\n".format(avg_syn_se))
        outputHandle.write("Average SE for all SNPs = {}\n".format(avg_se))
    
    print "Done."



    
            
        
        
        
