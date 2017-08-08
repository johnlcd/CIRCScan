#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
import random
os.environ

opt = sys.argv
def make_set():
	if '-h' in opt:
		print 'Arguments:\n','-t: Cell type\n','-R: Ratio of negative VS positive\n'
		pass
	else:
		CELL = opt[opt.index('-t')+1]
		PS = CELL + '_circ_true'
		UK = CELL + '_unknown'
		Train = CELL + '_train'
		Pred = CELL + '_pred'
		ratio = float(opt[opt.index('-R')+1])

		with open(Train,'w') as a:
			with open(Pred,'w') as b:
				with open(PS,'r') as c:
					with open(UK,'r') as d:
						header = c.readline()
						a.write(header)
						b.write(header)
						ps = c.readlines()
						tpn = len(ps)
						a.writelines(ps)
						
						d.readline()
						uk = d.readlines()
						random.shuffle(uk)
						ns = uk[0:int(tpn*ratio)]
						m = 0
						n = 0
						for i in ns:
							a.write(i.replace('UK','F'))
							m+=1
							if m%10000 == 0:
								print '%d pisitive circRNAs.\n' % m
						pred = uk[int(tpn*ratio):]
						for j in pred:
							b.writelines(j.replace('UK','P'))
							n+=1
							if n%100000 == 0:
								print '%d intervals for prediction.\n' % n
						print 'Total number of positive set is: %d,\nand pred set is: %d.\n' % (m, n)


if __name__ == '__main__':
	make_set()
