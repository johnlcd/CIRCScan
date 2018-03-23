#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
#import pp
import math
os.environ

opt = sys.argv
def anno_IP():
	if '-h' in opt:
		print '>>> Arguments:\n'
		print '-i: Input file (intron pair bedpe file)\n'
		print '-r: Reference file for intron and intron pair\n'
		pass
	else:
		IN = opt[opt.index('-i')+1]
		ANNO = 'pair_anno_alu'
		ANNO_M1 = 'pair_anno_alu_M1'
#		IL = 'IRAlu_intron_pair.list'
		REF = opt[opt.index('-r')+1]
#		FAIS = 'fw_alu_interval_score'
#		RAIS = 'rev_alu_interval_score'
		FAIO = 'fw_alu_intron_overlap'
		RAIO = 'rev_alu_intron_overlap'

		with open(ANNO, 'w') as a:
			with open(FAIO, 'r') as b:
				with open(RAIO, 'r') as c:
					with open(IN, 'r') as d:
						with open(REF, 'r') as e:
							with open(ANNO_M1, 'w') as f:

								IP = []
								INTRON = []
								DIPI1 = {}
								DIPI2 = {}
								DFI = {}
								DRI = {}
#								IP_POS = {}
								k = 0
								for re in e.readlines():
									k += 1
									le = re.strip().split('\t')
									i1 = le[1]
									i2 = le[2]
									ip = le[0]
									DIPI1[ip] = i1
									DIPI2[ip] = i2
									IP.append(ip)
									if k % 50000 == 0:
										print 'The %d intron pair:\t' % k, ip, ' [', i1, ', ', i2, ']\n'
									INTRON.append(i1)
									INTRON.append(i2)
								INTRON = list(set(INTRON))
								for i in INTRON:
									DFI[i] = 0
									DRI[i] = 0
	
	
								if len(DIPI1) == len(DIPI2):
									print 'Number of intron and score pairs: %d.\n' % len(DIPI1)
								else:
									print 'Number of forword and reverse intron (pairs) not consistent.\n'
	
								for rb in b.readlines():
									lb = rb.strip().split('\t')
									intron = lb[3]
#									print 'Intron is: %s\n' % intron
#									if intron not in INTRON:
#										print '%s not in intron list.\n' % intron
									fis = DFI.get(intron)
									DFI[intron] = fis + 1
	
								for rc in c.readlines():
									lc = rc.strip().split('\t')
									intron = lc[3]
									ris = DRI.get(intron)
									DRI[intron] = ris + 1
	
								m = 0
								n = 0
								a.write('Intron_pair' + '\t' + 'Alu' + '\n')
								f.write('Intron_pair' + '\t' + 'Alu' + '\n')
								for rd in d.readlines():
									m += 1
									ld = rd.strip().split('\t')
									ip = ld[6]
									if m % 10000 == 0:
										print 'The %d intron pair: %s\n' % (m, ip)
									i1 = DIPI1[ip]
									i2 = DIPI2[ip]
#									print 'Intron 1 position: (%s, %s)' % (ld[1], ld[2])
#									print 'Intron 2 position: (%s, %s)' % (ld[4], ld[5])
									l1 = int(ld[2]) - int(ld[1])
									l2 = int(ld[5]) - int(ld[4])
									if m % 10000 == 0:
										print 'Intron pair: %s\n' % ip
										print 'Introns: (%s,%s)\n' % (i1,i2)
										print 'Length of introns in intron pair %s: (%d, %d)' % (ip, l1, l2)
									fi1 = DFI.get(i1)
									fi2 = DFI.get(i2)
									ri1 = DRI.get(i1)
									ri2 = DRI.get(i2)
									iralu_across = fi1*ri2 + fi2*ri1
									iralu_within = fi1*ri1 + fi2*ri2
									iralu_score = '%.3f' % (float(iralu_across-iralu_within)*1000/math.sqrt(l1*l2))
									if m % 10000 == 0:
										print 'Number if IRAlus across and within: (%d, %d), length of intron 1 and intron 2: (%d, %d).\n' % (iralu_across, iralu_within, l1, l2)
#										print 'IRAlus score is: %s.\n' % iralu_score
									if m % 10000 == 0:
										print 'IRAlu score of intron pair (all) %s: %s.\n' % (ip, iralu_score)
									a.write(ip + '\t' + iralu_score + '\n')
									if iralu_across >= 1:
										n += 1
										f.write(ip + '\t' + iralu_score + '\n')
									if m % 10000 == 0:
										print 'IRAlu score of intron pair (with IRAlus across more than 1) %s: %s.\n' % (ip, iralu_score)

								print 'Number of total intron pairs: %d.\n' % m
								print 'Number of intron pairs with IRAlus across more than 1: %d.\n' % n


if __name__ == '__main__':
	anno_IP()
