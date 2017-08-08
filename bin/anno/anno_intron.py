#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
os.environ

opt = sys.argv
def anno_intron():
	if '-h' in opt:
		print 'Arguments:\n','-t: Cell type\n','-f: Feature\n','-r: reference intron\n','-l: feature list\n'
		pass
	else:
		CELL = opt[opt.index('-t')+1]
		FEA = opt[opt.index('-f')+1]
		RI = opt[opt.index('-r')+1]
		FL = opt[opt.index('-l')+1]

		with open('intron_anno_' + CELL + '_' + FEA,'w') as a:
			with open('overlap_' + CELL + '_' + FEA,'r') as b:
				with open(RI,'r') as c:
					with open(FL,'r') as d:
						Fea = []	
						for rd in d:
							Fea.append(rd.strip())	
						print 'Features are:\n',Fea,'\n'
	
						Intron = []
						for rc in c:
							Intron.append(rc.strip().split('\t')[3])
						print 'The number of intron:\n',len(Intron),'\n'
					
						D = {}
						for i in Intron:
							D[i] = []
						
						header = ['Chr','Start','End','INTRON'] + Fea 
						a.write('\t'.join(header) + '\n')
						hp = {}
						for rb in b:
							lb = rb.strip().split('\t')
							fea = lb[9]	
							pos = lb[0:3]
							intron = lb[3]
							hp[intron] = pos
							F = D.get(intron)	
							if fea in Fea and fea not in F:	
								F.append(fea)	
								D[intron] = F	
						
						n = 0
						for i in Intron:
							pos = hp.get(i)
							wl = pos + [i] + ['0']*len(Fea)
							F = D.get(i)
							for fea in F:
								wl[header.index(fea)] = '1'
							a.write('\t'.join(wl) + '\n')
							n+=1
							if n%10000 == 0:
								print 'This is the %s intron in the list,' % str(n)
								print 'and the annotation type of intron %s is:\n' % i
								print wl[3:len(Fea)+3],'\n'
						print 'Total of %d introns.\n' % n


if __name__ == '__main__':
	anno_intron()
