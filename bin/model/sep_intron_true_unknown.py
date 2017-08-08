#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
#import pp
os.environ

opt = sys.argv
def sep_pairs():
	if '-h' in opt:
		print 'Arguments:\n','-t: Cell type\n','--circ: known circRNA BED file\n'
		pass
	else:
		CELL = opt[opt.index('-t')+1]
		Pair_anno = CELL + '_anno_comb'
		circ_bed = opt[opt.index('--circ')+1]

		with open(CELL + '_circ_true','w') as a:
			with open(CELL + '_unknown','w') as b:
				with open(Pair_anno,'r') as c:
					with open(circ_bed,'r') as d:
						CIRC_pos = []
						for rd in d:
							ld = rd.strip().split('\t')
							pos = ld[0:3]
							CIRC_pos.append(pos)
	
						head = c.readline()
						a.write(head.strip() + '\tType\n')
						b.write(head.strip() + '\tType\n')
						n = 1
						for rc in c.readlines():
							lc = rc.strip().split('\t')
							pos = lc[0:3]
							if pos in CIRC_pos:
								a.write(rc.strip() + '\tT\n')
							else:
								b.write(rc.strip() + '\tUK\n')
							if n%100000 == 0:
								print 'The %s intron pair.\n' % str(n)
							n+=1

if __name__ == '__main__':
	sep_pairs()
