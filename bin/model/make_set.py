#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
import random
os.environ

opt = sys.argv
def make_set():
	if '-h' in opt:
		print 'Arguments:\n'
		print '-t: Cell type\n'
		print '-pp: List of intron pair with (close to) promoter\n'
		print '--circ: known circRNA BED intron pairs\n'
		print '-R: Ratio of negative VS positive\n'
		print '--sl: list of intron pair length to do stratified random sampling (comma seperated 3 number, defult: 10000 20000 30000 )\n'
		pass
	else:
		CELL = opt[opt.index('-t')+1]
		Pair_anno = CELL + '_anno_comb'
		PP_LIST = opt[opt.index('-pp')+1]
		Circ = opt[opt.index('--circ')+1]
		Train = CELL + '_train'
		Pred = CELL + '_pred'
		ratio = opt[opt.index('-R')+1]
		SL = opt[opt.index('--sl')+1].split(',')
		PART1 = CELL + "_IP_part1"
		PART2 = CELL + "_IP_part2"
		PART3 = CELL + "_IP_part3"
		PART4 = CELL + "_IP_part4"

		with open(Train, 'w') as a:
			with open(Pred, 'w') as b:
				with open(Pair_anno, 'r') as c:
					with open(Circ, 'r') as d:
						with open(PP_LIST, 'r') as f:
							with open(PART1, 'r') as P1:
								with open(PART2, 'r') as P2:
									with open(PART3, 'r') as P3:
										with open(PART4, 'r') as P4:


# Make intron pair list
											print 'List: ', SL, '\n'
											print ">>> Intron pair length to do stratified random sampling: 0 ~ %s, %s ~ %s, %s ~ %s, >= %s\n" % (SL[0], SL[0], SL[1], SL[1], SL[2], SL[2])
											print '    Start to make list of intron pairs, known circRNAs intron pairs and intron pairs with (close to) promoter ... ...\n'
				
											IP1 = []
											IP2 = []
											IP3 = []
											IP4 = []
											
											for pair in P1.readlines():
												p1 = pair.strip()
												IP1.append(p1)
												if len(IP1) % 10000 == 0:
													print '    %d intron pairs ...\n' % len(IP1)
											print '    %d intron pairs of PART 1.\n' % len(IP1)
											
											for pair in P2.readlines():
												p2 = pair.strip()
												IP2.append(p2)
												if len(IP2) % 10000 == 0:
													print '    %d intron pairs ...\n' % len(IP2)
											print '    %d intron pairs of PART 2.\n' % len(IP2)

											for pair in P3.readlines():
												p3 = pair.strip()
												IP3.append(p3)
												if len(IP3) % 10000 == 0:
													print '    %d intron pairs ...\n' % len(IP3)
											print '    %d intron pairs of PART 3.\n' % len(IP3)

											for pair in P4.readlines():
												p4 = pair.strip()
												IP4.append(p4)
												if len(IP4) % 10000 == 0:
													print '    %d intron pairs ...\n' % len(IP4)
											print '    %d intron pairs of PART 4.\n' % len(IP4)
											
											CIP = []
											for rd in d:
												cip = rd.strip()
												CIP.append(cip)
												if len(CIP) % 10000 == 0:
													print '    %d circRNA intron pairs ...\n' % len(CIP)
											print '    %d intron pairs mapped with known circRNA.\n' % len(CIP)

											PP = []
											for rf in f.readlines():
												PP.append(rf.strip())
											print '    Total of %d intron pairs with (close to) promoter.\n' % len(PP)
											

# Assign circRNA intron pairs to each part

											print '>>> Assign circRNA intron pairs to each part ... ... \n'
											print '[1] Part1:\n'
											CIP1_OUT = list(set(IP1) - set(CIP))
											CIP1_IN = list(set(IP1) - set(CIP1_OUT))
											print '    Number of PART 1 intron pairs in and out of circRNA intron pairs: (%d, %d)\n' % (len(CIP1_IN), len(CIP1_OUT))
											print '[2] Part2:\n'
											CIP2_OUT = list(set(IP2) - set(CIP))
											CIP2_IN = list(set(IP2) - set(CIP2_OUT))
											print '    Number of PART 2 intron pairs in and out of circRNA intron pairs: (%d, %d)\n' % (len(CIP2_IN), len(CIP2_OUT))
											print '[3] Part3:\n'
											CIP3_OUT = list(set(IP3) - set(CIP))
											CIP3_IN = list(set(IP3) - set(CIP3_OUT))
											print '    Number of PART 3 intron pairs in and out of circRNA intron pairs: (%d, %d)\n' % (len(CIP3_IN), len(CIP3_OUT))
											print '[4] Part4:\n'
											CIP4_OUT = list(set(IP4) - set(CIP))
											CIP4_IN = list(set(IP4) - set(CIP4_OUT))
											print '    Number of PART 4 intron pairs in and out of circRNA intron pairs: (%d, %d)\n' % (len(CIP4_IN), len(CIP4_OUT))
											print '    All circRNA intron pairs asigned.\n'


# Select intron pairs with (close to) or without (far from) promoter

											print '>>> Start to select intron pairs with (close to) or without (far from) promoter by part ... ...\n'

											print '[1] Part1:\n'
											NPCIP1_IN = list(set(CIP1_IN) - set(PP))
											PCIP1_IN = list(set(CIP1_IN) - set(NPCIP1_IN))
											NPCIP1_OUT = list(set(CIP1_OUT) - set(PP))
											PCIP1_OUT = list(set(CIP1_OUT) - set(NPCIP1_OUT))
											print '    Number of PART 1 circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP1_IN), len(NPCIP1_IN))
											print '    Number of PART 1 non-circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP1_OUT), len(NPCIP1_OUT))

											print '[2] Part2:\n'
											NPCIP2_IN = list(set(CIP2_IN) - set(PP))
											PCIP2_IN = list(set(CIP2_IN) - set(NPCIP2_IN))
											NPCIP2_OUT = list(set(CIP2_OUT) - set(PP))
											PCIP2_OUT = list(set(CIP2_OUT) - set(NPCIP2_OUT))
											print '    Number of PART 2 circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP2_IN), len(NPCIP2_IN))
											print '    Number of PART 2 non-circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP2_OUT), len(NPCIP2_OUT))

											print '[3] Part3:\n'
											NPCIP3_IN = list(set(CIP3_IN) - set(PP))
											PCIP3_IN = list(set(CIP3_IN) - set(NPCIP3_IN))
											NPCIP3_OUT = list(set(CIP3_OUT) - set(PP))
											PCIP3_OUT = list(set(CIP3_OUT) - set(NPCIP3_OUT))
											print '    Number of PART 3 circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP3_IN), len(NPCIP3_IN))
											print '    Number of PART 3 non-circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP3_OUT), len(NPCIP3_OUT))

											print '[4] Part4:\n'
											NPCIP4_IN = list(set(CIP4_IN) - set(PP))
											PCIP4_IN = list(set(CIP4_IN) - set(NPCIP4_IN))
											NPCIP4_OUT = list(set(CIP4_OUT) - set(PP))
											PCIP4_OUT = list(set(CIP4_OUT) - set(NPCIP4_OUT))
											print '    Number of PART 4 circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP4_IN), len(NPCIP4_IN))
											print '    Number of PART 4 non-circRNA intron pairs with and without promoter: (%d, %d)\n' % (len(PCIP4_OUT), len(NPCIP4_OUT))

											
# Random select intron pairs											
				
											print '>>> Do stratified random sampling of negative intron pairs with ratio of %s according to intron length and distance to promoter ... ...\n' % ratio

											random.shuffle(NPCIP1_OUT)
											random.shuffle(PCIP1_OUT)
											NEG_PCIP1 = PCIP1_OUT[0:int(len(PCIP1_IN) * float(ratio))]
											PRED_PCIP1 = PCIP1_OUT[int(len(PCIP1_IN) * float(ratio)):]
											NEG_NPCIP1 = NPCIP1_OUT[0:int(len(NPCIP1_IN) * float(ratio))]
											PRED_NPCIP1 = NPCIP1_OUT[int(len(NPCIP1_IN) * float(ratio)):]
											print '    Number of negative and predicted intron pairs with promoter in part 1: (%d, %d)\n' % (len(NEG_PCIP1), len(PRED_PCIP1))
											print '    Number of negative and predicted intron pairs without promoter in part 1: (%d, %d)\n' % (len(NEG_NPCIP1), len(PRED_NPCIP1))

											random.shuffle(NPCIP2_OUT)
											random.shuffle(PCIP2_OUT)
											NEG_PCIP2 = PCIP2_OUT[0:int(len(PCIP2_IN) * float(ratio))]
											PRED_PCIP2 = PCIP2_OUT[int(len(PCIP2_IN) * float(ratio)):]
											NEG_NPCIP2 = NPCIP2_OUT[0:int(len(NPCIP2_IN) * float(ratio))]
											PRED_NPCIP2 = NPCIP2_OUT[int(len(NPCIP2_IN) * float(ratio)):]
											print '    Number of negative and predicted intron pairs with promoter in part 2: (%d, %d)\n' % (len(NEG_PCIP2), len(PRED_PCIP2))
											print '    Number of negative and predicted intron pairs without promoter in part 2: (%d, %d)\n' % (len(NEG_NPCIP2), len(PRED_NPCIP2))

											random.shuffle(NPCIP3_OUT)
											random.shuffle(PCIP3_OUT)
											NEG_PCIP3 = PCIP3_OUT[0:int(len(PCIP3_IN) * float(ratio))]
											PRED_PCIP3 = PCIP3_OUT[int(len(PCIP3_IN) * float(ratio)):]
											NEG_NPCIP3 = NPCIP3_OUT[0:int(len(NPCIP3_IN) * float(ratio))]
											PRED_NPCIP3 = NPCIP3_OUT[int(len(NPCIP3_IN) * float(ratio)):]
											print '    Number of negative and predicted intron pairs with promoter in part 3: (%d, %d)\n' % (len(NEG_PCIP3), len(PRED_PCIP3))
											print '    Number of negative and predicted intron pairs without promoter in part 3: (%d, %d)\n' % (len(NEG_NPCIP3), len(PRED_NPCIP3))

											random.shuffle(NPCIP4_OUT)
											random.shuffle(PCIP4_OUT)
											NEG_PCIP4 = PCIP4_OUT[0:int(len(PCIP4_IN) * float(ratio))]
											PRED_PCIP4 = PCIP4_OUT[int(len(PCIP4_IN) * float(ratio)):]
											NEG_NPCIP4 = NPCIP4_OUT[0:int(len(NPCIP4_IN) * float(ratio))]
											PRED_NPCIP4 = NPCIP4_OUT[int(len(NPCIP4_IN) * float(ratio)):]
											print '    Number of negative and predicted intron pairs with promoter in part 4: (%d, %d)\n' % (len(NEG_PCIP4), len(PRED_PCIP4))
											print '    Number of negative and predicted intron pairs without promoter in part 4: (%d, %d)\n' % (len(NEG_NPCIP4), len(PRED_NPCIP4))
											
											NEG = NEG_PCIP1 + NEG_NPCIP1 + NEG_PCIP2 + NEG_NPCIP2 + NEG_PCIP3 + NEG_NPCIP3 + NEG_PCIP4 + NEG_NPCIP4
											PRED = PRED_PCIP1 + PRED_NPCIP1 + PRED_PCIP2 + PRED_NPCIP2 + PRED_PCIP3 + PRED_NPCIP3 + PRED_PCIP4 + PRED_NPCIP4
											print '    %d negative intron pairs.\n' % len(NEG)
											print '    %d intron pairs for prediction.\n' % len(PRED)
				
											print '>>> Start to build data set for training and prediction ... ...\n'
#											TRAIN = CIP + NEG
											header = c.readline()
											a.write(header.strip() + '\tType\n')
											b.write(header.strip() + '\tType\n')
											nip = 0
											npo = 0
											nne = 0
											npr = 0
											for rc in c.readlines():
												nip += 1
												if nip % 10000 == 0:
													print '    The %d intron pairs.\n' % nip
												lc = rc.strip().split('\t')
												ip = lc[3]
												if ip in CIP:
													npo += 1
													a.write('\t'.join(lc + ['T']) + '\n')
													if npo % 100 == 0:
														print '    The %d positive intron pairs.\n' % npo
												if ip in NEG:
													nne += 1
													a.write('\t'.join(lc + ['F']) + '\n')
													if nne % 100 == 0:
														print '    The %d negative intron pairs.\n' % nne
												if ip in PRED:
													npr += 1
													b.write('\t'.join(lc + ['P']) + '\n')
													if npr % 1000 == 0:
														print '    The %d intron pairs for prediction.\n' % npr
											
											print '>>> Total of %d intron pairs.\n' % nip
											print '    Number of intron pairs for traning: %d (%d positive, %d negative),\nand intron pairs for prediction is: %d.\n' % (npo + nne, npo, nne, npr)


if __name__ == '__main__':
	make_set()
