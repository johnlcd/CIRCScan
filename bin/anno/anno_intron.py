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
						print '>>> Make feature list ... ...\n'
						Fea = []	
						for rd in d:
							Fea.append(rd.strip())	
						print '    Features are:\n',Fea,'\n'
	
						print '>>> Make intron list ... ...\n'
						D_IL = {}
						D_IP = {}
						Intron = []
						for rc in c:
							lc = rc.strip().split('\t')
							i_chr = lc[0]
							i_sts = lc[1]
							i_end = lc[2]
							pos = lc[0:3]
							i_len = int(i_end) - int(i_sts)
							intron = rc.strip().split('\t')[3]
							Intron.append(intron)
							if len(Intron) % 10000 == 0:
								print '    The %d intron: %s\n' % (len(Intron), intron)
							D_IL[intron] = i_len
							D_IP[intron] = pos
						print '>>> The number of intron:\n    ', len(Intron), '\n'
					
						print '>>> Make intron feature pairs ... ...\n'
						D_IFL = {}
						for i in Intron:
							for fea in Fea:
								i_fea = i + '_' + fea
								D_IFL[i_fea] = 0
						print '    %d intron feature pairs.\n' % len(D_IFL)
						
						print '>>> Calculate the feature length in each intron ... ...\n'
						header = ['Chr','Start','End','INTRON'] + Fea 
						a.write('\t'.join(header) + '\n')
						nn = 0
						for rb in b:
							nn += 1
							lb = rb.strip().split('\t')
							intron = lb[3]
							if ( len(lb) == 11 ) or ( len(lb) == 13 ):
								fea = lb[9]
								o_len = int(lb[-1])
								i_fea = intron + '_' + fea
								len_sum = D_IFL.get(i_fea)	
								D_IFL[i_fea] = len_sum + o_len

						print '>>> Calculate the density of feature in intron region ... ...\n'
						n = 0
						for intron in Intron:
							n += 1
							if n % 10000 == 0:
								print '>>> This is the %s intron in the list.\n' % str(n)
							pos = D_IP.get(intron)
							wl = pos + [intron]
							for fea in Fea:
								i_fea = intron + '_' + fea
								o_len_sum = D_IFL.get(i_fea)
								i_len = D_IL.get(intron)
								if n % 10000 == 0:
									print '    Intron length of intron %s: %d, overlapped length with %s: %d. \n' % (intron, i_len, fea , o_len_sum)
								dense = "%.3f" % (float(o_len_sum)/float(i_len))
								wl = wl + [str(dense)]
							a.write('\t'.join(wl) + '\n')
							n+=1
							if n % 50000 == 0:
								print '    Annotation of %d intron %s is:\n' % (n, intron)
								print wl[3:len(Fea)+3],'\n'
						print '    Total of %d introns.\n' % n


if __name__ == '__main__':
	anno_intron()
