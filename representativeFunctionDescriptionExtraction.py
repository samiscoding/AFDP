'''
@auther: Samaneh

FINAL VERSION/ UPDATED ON 7.November.2016

'''
from Bio import SeqIO
import pandas as pd
import numpy as np
import xlsxwriter
import re
import multiprocessing
from multiprocessing.dummy import Pool
#from suffix_tree import SuffixTree
#from suffix_tree import GeneralisedSuffixTree 
import suffixtreeLibrary as st
import json
#from django.utils.encoding import smart_str

##################################################################################
##################################################################################
##################################################################################

class clusters():

	def __init__(self):

		self.suffixDic = dict()
		self.longest_common_string = ""
		self.prefix_longest_common_string = ""
		self. suffix_longest_common_string = ""


	def findMax(self,descList):

		maxOne = ""
		for s in descList:
			if len(s)>len(maxOne):
				maxOne = s
		return maxOne		

	def prefixCheck(self, ncp1):

		pPrefix = [] 
		pSuffix = []
		prefixOk = False
		prefixFlag = False

		prefixStree = st.STree(input=ncp1)
		self.prefix_longest_common_string, prefixWii = prefixStree.lcs()
		self.prefix_longest_common_string = self.prefix_longest_common_string.encode("ascii","ignore").lstrip(" ").rstrip(" ").lower()
		words = self.prefix_longest_common_string.split()
		for j in range(len(ncp1)):
			prefixOkCount = 0
			for w in words:
				if w in ncp1[j].split():
					prefixOkCount = prefixOkCount + 1
			if prefixOkCount == len(words):
				prefixOk = True  		

		if prefixOk == True:
			for i in range(len(ncp1)):
				prefix = suffix = ""
				plcsStartingIndx = ncp1[i].find(self.prefix_longest_common_string) 
				plcsEndingIndx = plcsStartingIndx + len(self.prefix_longest_common_string)
				if plcsStartingIndx > 0:
					pPrefix.append(ncp1[i][0:plcsStartingIndx])
				if plcsEndingIndx < len(ncp1[j]):
					pSuffix.append(ncp1[i][plcsEndingIndx:])
			if pPrefix:
				prefix = self.prefixCount(pPrefix)
			if pSuffix:
				suffix = self.suffixCount(pSuffix)
			self.longest_common_string = prefix + self.prefix_longest_common_string + suffix + " " + self.longest_common_string
			prefixFlag = True	
		
		return prefixFlag	


	def suffixCheck(self, ncp2):
		
		sPrefix = []
		sSuffix = []
		suffixOk =False
		suffixFlag = False

		suffixStree = st.STree(input=ncp2)
		self.suffix_longest_common_string, suffixWii = suffixStree.lcs()
		self.suffix_longest_common_string = self.suffix_longest_common_string.encode("ascii","ignore").lstrip(" ").rstrip(" ").lower()
		#print ncp2, self.suffix_longest_common_string

		words = self.suffix_longest_common_string.split()
		for j in range(len(ncp2)):
			suffixOkCount = 0
			for w in words:
				if w in ncp2[j].split():
					suffixOkCount = suffixOkCount + 1
			if suffixOkCount == len(words):
				suffixOk = True  

		if suffixOk == True:
			sLength = len(self.suffix_longest_common_string) 
			if sLength == 1:
				self.longest_common_string = self.longest_common_string.strip(" ") + "/" + self.longest_common_string.split()[-1] + self.suffix_longest_common_string
				suffixFlag = True
			else:
				for i in range(len(ncp2)):
					prefix = suffix = ""
					slcsStartingIndx = ncp2[i].find(self.suffix_longest_common_string) 
					slcsEndingIndx = slcsStartingIndx + len(self.suffix_longest_common_string)
					if slcsStartingIndx > 0:
						sPrefix.append(ncp2[i][0:slcsStartingIndx])
					if slcsEndingIndx < len(ncp2[i]):
						sSuffix.append(ncp2[i][slcsEndingIndx:])
				if sPrefix:
					prefix = self.prefixCount(sPrefix)
				if sSuffix:
					suffix = self.suffixCount(sSuffix)

				#print self.suffix_longest_common_string
				self.longest_common_string = self.longest_common_string + " " + prefix + " " + self.suffix_longest_common_string + " " + suffix
				suffixFlag = True

		return suffixFlag


	def prefixCount(self, ncp1):

		prefixKeyList = []
		prefixDict = {i:ncp1.count(i) for i in ncp1} #count repeated words of prefixes
		
		prefixMax = max(prefixDict.values())

		prefix = ""
		for key, value in prefixDict.iteritems():	
			if value == prefixMax:
				prefixKeyList.append(key)
		if prefixKeyList[0] != " ":						
			prefix = prefixKeyList[0]
		for k in range(1,len(prefixKeyList)):
			if prefixKeyList[k] != " ":
				prefix = prefix + "/" + prefixKeyList[k]	
		return prefix

	def suffixCount(self, ncp2):

		suffixKeyList = []
		suffixDict = {j:ncp2.count(j) for j in ncp2} #count repreated words of suffixes

		suffixMax = max(suffixDict.values())

		suffix = ""
		for key, value in suffixDict.iteritems():
			if value == suffixMax:
				suffixKeyList.append(key)	
		if suffixKeyList[0] != " ":
			suffix = suffixKeyList[0]
		for k in range(1,len(suffixKeyList)):
			if suffixKeyList[k] != " " and not re.match("\d",suffixKeyList[k]):	##@sam## numbers should not be considered as suffixes
				suffix = suffix + "/" + suffixKeyList[k]	

		return suffix

	def extractDescriptions(self,wii):

		suffList = []
		charList = []
		wiiStr = wii.encode("utf8","ignore")		
		##@sa## extracting all descriptions which are return from lcs functions of suffixtreeLibrary
		UPPAs = list(list(range(0xE000,0xF8FF+1)) + list(range(0xF0000,0xFFFFD+1)) + list(range(0x100000, 0x10FFFD+1)))
		for i in range(len(UPPAs)):
			y = unichr(UPPAs[i])
			charList.append(y)

		##@sam## a list of all descriptions
		suffList.append(wiiStr)

		##@sam## generating a list of all descriptions
		for c in charList:
			convertedC = c.encode("utf8","ignore")
			for j in range(len(suffList)):
				if convertedC in suffList[j]:
					temp = suffList[j]
					del suffList[j]	
					for token in temp.split(convertedC): 	
						token = token.replace("putative","").replace("(fragment)","").replace("(fragments)","").replace(" truncated","").replace("truncated","").replace("homolog","").replace("probable","").replace("(predicted)","")
						token = token.lstrip().rstrip()
						suffList.append(token)
		del suffList[len(suffList)-1]
		return suffList


	def pairwiseCheck(self, inputList):	##@sam## those clusters with sequences having so far descriptions from each other (according to the length of descriptions) will be removed at the begining  
		
		removeFlag = False
		base = max(inputList)
		threshold = (0.8*len(base))
		for i in range(len(inputList)):
			if len(inputList[i]) < threshold:
				removeFlag = True

	def makeSuffixTree(self, names, descriptions, resFile):

		suffixDict = dict()
		resultDic = dict() #dictionary of sequence ids with corresponding selected (representative) descriptions
		#longestLength = (l(descriptions))
		if not self.pairwiseCheck(descriptions):
			stree = st.STree(input=descriptions)
			self.longest_common_string, wii = stree.lcs()
			suffList = self.extractDescriptions(wii)

		##@sam## return all descriptions of a cluster and its selected longest common part of descriptions whose length is at least as long as 5% of the longest description of the cluster
			ncpDict = dict()	#dictionary whose keys are all non-conservative part (i.e prefix or suffix)(=NCP) and the values are the number of appearance in descriptions 
			compDict = dict()
			keysList = []

			self.longest_common_string = self.longest_common_string.encode("ascii","ignore").lstrip(" ").rstrip(" ") ##@sam## convert unicode to string
			minimalLength = len(self.findMax(descriptions))*0.2 ##@sam## define the minimal lenth of accepted common description
			if len(self.longest_common_string) >= minimalLength and self.longest_common_string != "protein":  ##@sam## clusters with totally different descriptions should be removed entirely
				ncp1 = []
				ncp2 = [] # list of prefixes and suffixes of longest_common_string of all descriptions
				wholeDes = ""
				allTheSame = False
				for su in suffList:
					cpStartingPoint = su.find(self.longest_common_string)	#finds the starting index of longest_common_string
					lcsLen = len(self.longest_common_string)	
					cpFinishingPoint = cpStartingPoint + lcsLen #finds the ending index of longest_commom_string
					suFinishingPoint = len(su)
					if cpStartingPoint > 0:
						ncp1.append(su[0:cpStartingPoint])
						
					if cpFinishingPoint < suFinishingPoint:
						
						ncp2.append(su[cpFinishingPoint:suFinishingPoint])

				### check prefixes and suffixes for possible longest_common_string extention 
				if ncp1:
					pCheck = self.prefixCheck(ncp1)	#if prefixCheck returns True i.e. self.longest_common_string has been extended
					if pCheck == False:
						prefix = self.prefixCount(ncp1)	#otherwise it tries to find words from prefixes to add to self.longest_common_string
						self.longest_common_string = prefix + " " + self.longest_common_string
				if ncp2:
					sCheck = self.suffixCheck(ncp2)	#if suffixCheck returns True i.e. self.longest_common_string has been extended	
					if sCheck == False:
						suffix = self.suffixCount(ncp2)	#otherwise it tries to find words from suffixes to add to self.longest_common_string	
						self.longest_common_string = self.longest_common_string + " " + suffix 
				
				self.longest_common_string = self.longest_common_string.replace("putative","").replace("(fragment)","").replace("(fragments)","").replace(" truncated","").replace("truncated","").replace("homolog","").replace("probable","").replace("(predicted)","")
			
				resFile.write("longest common:" + "\n" + self.longest_common_string + "\n" + "all descriptions:" + str(suffList) + "\n" + "*******" + "\n")


				for n in names:
					resultDic.update({n:self.longest_common_string})
				##@sam## write the dictionary into a json file:
				#outFile = open("/home/samaneh/AHRD/clustering/seqIDsAndDescriptions","w")
				#json.dump(resultDic, outFile)	
			return resultDic

########################################
	# @staticmethod ##@sam## create static method:
	def extractDesc(self,clstrsFile,resFile):
		
		tempFlag = False ### temporary flag until the file contains \n as the ending line
		seqNames = []
		names = []
		descs = []
		lines = clstrsFile.readlines()

		for i in range(len(lines)):
			if lines[i][0] == "#":
				tempFlag = False
				seqNames = []
				names = []
				descs = []
			elif lines[i][0] == "\n":
				if tempFlag == True:
					resultDictionary = self.makeSuffixTree(names, descs, resFile)	

			else:
				tempFlag = True
				seqNames = lines[i].split("|")
				for i in range(0,2):
					seqNames[i] = seqNames[i].strip(" ").strip("\n")
				seqNames[1] = seqNames[1].replace("putative","").replace("(fragment)","").replace("(fragments)","").replace(" truncated","").replace("truncated","").replace("homolog","").replace("probable","").replace("(predicted)","")
				seqNames[1] = seqNames[1].lstrip().rstrip()
				names.append(seqNames[0])
				descs.append(seqNames[1])	

		####from here you can start generateing new fata file using resultdictionary:


def handler():

	d = clusters()
	clustersFile = open("/home/samaneh/AHRD/clustering/clusteredDescriptions_final.txt","r")
	resultFile = open("/home/samaneh/AHRD/clustering/testClusterDescriptions_onSprot.txt","w")
	d.extractDesc(clustersFile,resultFile)

if __name__ == "__main__":
	handler()


