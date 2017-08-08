#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
os.environ

opt = sys.argv

def comb_anno_pair():
	if '-h' in opt:
		print 'Arguments:\n','-t: Cell type\n','-f: Feature (exclude hi-c)\n','-r: reference intron pairs\n'
		pass
	else:
		CELL = opt[opt.index('-t')+1]
		FEA = opt[opt.index('-f')+1]
		RIP = opt[opt.index('-r')+1]
		print 'Reference file is:\n',RIP,'\n'

		with open('pair_anno_' + CELL + '_' + FEA,'w') as a:
			with open('intron_anno_' + CELL + '_' + FEA,'r') as b:
				with open(RIP,'r') as c:	
					Pair = []
					P1 = {}
					P2 = {}
					for rc in c:
						lc = rc.strip().split('\t')
						pair = lc[0]
						p1 = lc[1]
						p2 = lc[2]
						Pair.append(pair)
						P1[pair] = p1
						P2[pair] = p2
					print 'Total of %d intron pairs.\n' % len(Pair)
	
					L = {}
					r1 = b.readline()
					l1 = r1.strip().split('\t')
					Feature = l1[4:len(l1)]
					print 'Featues are:'
					print Feature,'\n'
					wl1 = ['Intron_pair'] + Feature
					print 'Header of intron annotaion file:'
					print wl1,'\n'
					a.write('\t'.join(wl1) + '\n')
					for rb in b.readlines():
						lb = rb.strip().split('\t')
						intron = lb[3]
						L[intron] = lb
					print 'Intron number:\t',len(L),'\n'
					
					n = 0
					for pair in Pair:
						wl = [pair] + ['0']*len(Feature)
						i1 = P1.get(pair)
						i2 = P2.get(pair)
						li1 = L.get(i1)
						li2 = L.get(i2)
						for f in Feature:
							value1 = li1[l1.index(f)]
							value2 = li2[l1.index(f)]
							value = int(value1) + int(value2)
							wl[wl1.index(f)] = str(value)
							n+=1
						if n%100000 == 0:
							print 'Intron pair number is:\n',n
							print 'Intron pairs are:\n',i1,',',i2
							print 'Pair name is:',pair,'\n'
						a.write('\t'.join(wl) + '\n')


if __name__ == '__main__':
	comb_anno_pair()
