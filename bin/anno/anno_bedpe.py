#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
os.environ

opt = sys.argv
def anno_pairs():
	if '-h' in opt:
		print 'Arguments:\n','-t: Cell type\n','-f: Feature\n','-r: reference intron pairs\n','-l: feature list\n'
		pass
	else:
		CELL = opt[opt.index('-t')+1]
		FEA = opt[opt.index('-f')+1]
		RIP = opt[opt.index('-r')+1]
		FL = opt[opt.index('-l')+1]
		print 'Reference file is:\n',RIP,'\n'
		
		with open('raw_pair_anno_' + CELL + '_' + FEA,'w') as a:
			with open('overlap_' + CELL + '_' + FEA,'r') as b:
				with open(RIP,'r') as c:
					with open(FL,'r') as d:
						Fea = []
						for rd in d:
							Fea.append(rd.strip())	
						print '>>> Features are:\n',Fea,'\n'
						
						Intron_pair = []
						for rc in c:
							Intron_pair.append(rc.strip().split('\t')[6])
						print '>>> The number of intron pairs:\n',len(Intron_pair),'\n'
						
						D = {}
						for i in Intron_pair:
							D[i] = []

						header = ['Intron_pair'] + Fea
						a.write('\t'.join(header) + '\n')
						for rb in b:
							lb = rb.strip().split('\t')
							fea = lb[16]
							intron_pair = lb[6]
							F = D.get(intron_pair)	
							if fea in Fea and fea not in F:	
								F.append(fea)
								D[intron_pair] = F	

						n = 0
						for i in Intron_pair:
							wl = [i] + ['0']*len(Fea)
							F = D.get(i)
							for fea in F:
								wl[header.index(fea)] = '1'
							a.write('\t'.join(wl) + '\n')
							n+=1
							if n%10000 == 0:
								print '>>> This is the %s intron pair in the list,' % str(n)
								print '    and the annotation type of intron pair %s is:\n' % i
								print wl[0:len(Fea)],'\n'
						print '>>> Total of %d intron pairs.\n' % n


if __name__ == '__main__':
	anno_pairs()
