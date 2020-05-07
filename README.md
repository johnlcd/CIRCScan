
# CIRCScan  
> Tools for predicting circRNAs expression with ***`EPIGENETIC`*** features by ***`Machine Learning`***  
<br>

## Introduction  

> Circular RNAs ( circRNAs ) are an abundant class of noncoding RNAs with the widespread, cell/tissue specific pattern.  
> This tool, `CIRCScan`, is used for predicting circRNAs expression in a cell/tissue specific manner by machine learning based on epigenetic features.  
<br>
  
## License
> This software is distributed under the terms of GPL 2.0  
<br>

## Source
> [https://github.com/johnlcd/CIRCScan](https://github.com/johnlcd/CIRCScan)  
<br>

## Contact
### Author
> **Jia-Bin Chen**, **Shan-Shan Dong**, **Shi Yao**, **Yan Guo**, **Tie-Lin Yang**  
> Key Laboratory of Biomedical Information Engineering of Ministry of Education, School of Life Science and Technology, Xi'an Jiaotong University, Xi'an, Shaanxi Province, 710049, P. R. China  
> [yangtielin@mail.xjtu.edu.cn](yangtielin@mail.xjtu.edu.cn)  
<br>

## Maintainer
> **Jia-Bin Chen**  
> You can contact [johnlcd@stu.xjtu.edu.cn](johnlcd@stu.xjtu.edu.cn) when you have any questions, suggestions, comments, etc.  
> Please describe in details, and attach your command line and log messages if possible.  
<br>

## Requiremnets
- **bedtools** \( v2.25.0 \)
- **Python** \( recommended: python2.7 \)
- **R** \( >= 3.2.4 \)
	- R packages: caret (6.0-73), ggplot2 (2.2.0), doParallel (3.2.4), ROCR (1.0-7), etc. ( Dependent packages for different models ) 
		
> Check the log file ".out" to validate which package is required if get an error info
<br>

## Building `CIRCscan`

***CMD:***  

		git clone https://github.com/johnlcd/CIRCScan.git

<br>

## Directory catalog

- **bin**  	
	- **anno**  
		- alu_anno_IP.py
		- anno_bedpe.py  
		- anno_intron.py
		- bt_intersect_alu_intron.sh
		- bt_overlap
		- comb_pair_anno.py
	- anno_pair
	- circscan
	- merge_feature
	- **model**
		- Circ_pred.R
		- Eval_test_cor.R
		- Eval_test_perf.R
		- feature_selection.R
		- make_set.py
		- Model_train_FS_exp.R
		- Model_train_FS_exp_rep.R
		- Model_train_FS_exp_rmo.R
		- Model_train_FS_exp_rmo_rep.R
		- Model_train_pred_exp.R
		- Model_train.R
		- Model_train_rep.R
	- prepare_train_set
- **data**
	- **Alu**
		- alu_hg19.bed
		- pair_anno_alu
	- **anno**
		- GM12878_anno_comb.gz
		- H1-hESC_anno_comb.gz
		- HeLa-S3_anno_comb.gz
		- HepG2_anno_comb.gz
		- K562_anno_comb.gz
		- NHEK_anno_comb.gz
	- **DNaseI**
		- GM12878_dnase.bed
		- H1-hESC_dnase.bed
		- HeLa-S3_dnase.bed
		- HepG2_dnase.bed
		- HUVEC_dnase.bed
		- K562_dnase.bed
		- NHEK_dnase.bed
	- **histone**
		- A549_his.bed
		- GM12878_his.bed
		- H1-hESC_his.bed
		- HeLa-S3_his.bed
		- HepG2_his.bed
		- HUVEC_his.bed
		- K562_his.bed
		- NHEK_his.bed
	- **known_circ**
		- A549_circ_overlap.bed
		- GM12878_circ_overlap.bed
		- H1-hESC_circ_overlap.bed
		- HeLa-S3_circ_overlap.bed
		- HepG2_circ_overlap.bed
		- HMEC_circ_overlap.bed
		- HOB_circ_overlap.bed
		- HSMM_circ_overlap.bed
		- HUVEC_circ_overlap.bed
		- K562_circ_overlap.bed
		- NHA_circ_overlap.bed
		- NHDF_circ_overlap.bed
		- NHEK_circ_overlap.bed
		- NHLF_circ_overlap.bed
	- **pred_true_bycell**
		- GM12878_pred_true.bed.gz
		- H1-hESC_pred_true.bed.gz
		- HeLa-S3_pred_true.bed.gz
		- HepG2_pred_true.bed.gz
		- K562_pred_true.bed.gz
		- NHEK_pred_true.bed.gz
	- **raw_data**
		- DNaseI.txt.gz
		- Histone_part1.txt.gz
		- Histone_part2.txt.gz
		- select_cell.list
		- select_his.list
- **info**
	- models_ALL.txt
	- models_classification.list
	- models_test.list
	- models_test.txt
- README.md
- **sample**
	- **anno**
	- **feature**
	- **model**
<br>


## Runing preparation
- ### Set environment variables

***CMD:***  

		export PKG_DIR=/path/to/CIRCScan
		export PATH=$PKG_DIR/bin:PKG_DIR/bin/anno:$PKG_DIR/bin/model:$PATH

- ### Unzip data files

***CMD:***  

		cd $PKG_DIR/data
		tar -zxvf intron_pairs_data.tgz
		cd $PKG_DIR/data/raw_data
		gunzip *.gz
		cat Histone_part1.txt Histone_part2.txt > Histone.txt
<br>


## Work flow
<br>

### 1. Data preparation and feature generation
> Epigenetic data including DNaseI HS, Histone modification, downloaded from `ENCODE` [http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/)  
- ### Process the raw broadPeak/narrowPeak file of ENCODE

***CMD:***

For each histone mark: (e.g.: H3K36me3)
		
		zcat wgEncodeBroadHistoneK562H3k36me3StdPk.broadPeak.gz | awk -v OFS='\t' '{print "K562","None","H3K36me3",$1,$2,$3}' > K562_H3K36me3.txt

then merged all marks of selected cell lines:

		#rm -f Histon.txt
		cat <(echo "cell treatment antibody chr start end" | sed 's/ /\t/g') *.txt > Histon.txt (Put these '.txt' file in one directory)


For DHS:

		zcat wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak.gz | awk -v OFS='\t' '{print "K562","None",$1,$2,$3}' > K562_DNaseI.txt

then merged all cell line:

		#rm -f DNaseI.txt
		cat <(echo "cell treatment chr start end" | sed 's/ /\t/g') *.txt > DNaseI.txt
		

- ### Extract features data from `.txt` file, transform into `BED` fromate

***CMD:***  

		grep K562 Histon.txt | grep -f $PKG_DIR/data/raw_data/select_his.list | awk -v OFS='\t' '{print $4,$5,$6,$3}' > K562_his.bed

		grep K562 DNaseI.txt | awk -v OFS='\t' '{print $4,$5,$6,"DNaseI_HS"}' > K562_dnase.bed

> `"Histon.txt"`  

	cell	treatment	antibody	chr	start	end  
	GM12878	None	CTCF	chr22	16846634	16869580  
	GM12878	None	CTCF	chr22	16850639	16850924  
	GM12878	None	CTCF	chr22	16851700	16851834  
	GM12878	None	CTCF	chr22	16852344	16852458  
	GM12878	None	CTCF	chr22	16853076	16853192  
	GM12878	None	CTCF	chr22	16853755	16853871  
	GM12878	None	CTCF	chr22	16854517	16854638  
	GM12878	None	CTCF	chr22	16857119	16857231  
	GM12878	None	CTCF	chr22	16857764	16857871  
	...  

> `"DNaseI.txt"`  

	lab	cell	treatment	chr	start	end
	Duke	8988T	None	chr1	564665	564815
	Duke	8988T	None	chr1	565025	565175
	Duke	8988T	None	chr1	565865	566015
	Duke	8988T	None	chr1	714005	714155
	Duke	8988T	None	chr1	762785	762935
	Duke	8988T	None	chr1	766705	766855
	Duke	8988T	None	chr1	767945	768095
	Duke	8988T	None	chr1	794145	794295
	Duke	8988T	None	chr1	795945	796095
	...  

> `"K562_his.bed"`  

	chr1	10140	10374	H3K9me3
	chr1	118494	118714	H3K9ac
	chr1	118556	118713	H3K4me3
	chr1	137502	140080	H3K9ac
	chr1	138030	140084	H3K4me2
	chr1	138411	138738	H3K4me3
	chr1	138424	138651	H3K27ac
	chr1	138426	138651	H3K79me2
	chr1	138934	139174	H3K4me3
	chr1	138938	139177	CTCF
	...  

> `"K562_dnase.bed"`  

	chr1	115600	115750	DNaseI_HS
	chr1	136280	136430	DNaseI_HS
	chr1	138960	139110	DNaseI_HS
	chr1	235040	235190	DNaseI_HS
	chr1	235600	235750	DNaseI_HS
	chr1	237640	237790	DNaseI_HS
	chr1	521460	521610	DNaseI_HS
	chr1	564480	564630	DNaseI_HS
	chr1	565280	565430	DNaseI_HS
	chr1	565860	566010	DNaseI_HS
	...


- ### Feature generation and annotation:

***1. Histone modifications, DNaseI HS ... ( Feature types of `"bed"` format )***  

Make feature list, and overlap intron with feature, annotate intron by features, then combine intron annotation to pair ( `"anno_pair"` )

***CMD:***  

		anno_pair -t <cell_type> -f <feature (his, dnase ...)> [ --is (Ignor strands) ] --bed <feature.bed>

		e.g.:
		anno_pair -t K562 -f his --is --bed K562_his.bed

		anno_pair -t K562 -f dnase --is --bed K562_dnase.bed

Generate 4 files:  

> `"K562_his.list"`  

	CTCF  
	EZH2_(39875)  
	H2A.Z  
	H3K27ac  
	H3K27me3  
	H3K36me3  
	H3K4me1  
	H3K4me2  
	H3K4me3  
	H3K79me2  
	H3K9ac  
	H3K9me3  
	H4K20me1  

> `"overlap_K562_his"` 

	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	712769	712874	H3K79me2	105
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713056	713748	H3K4me3	607
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713188	713524	H3K27ac	336
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713195	713556	H3K4me2	361
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713199	713548	H3K79me2	349
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713205	713560	H3K9ac	355
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713575	713751	H3K27ac	88
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713578	713747	H3K4me2	85
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713578	713752	H3K79me2	85
	chr1	709660	713663	LOC100288069-1-1	1.00	-	chr1	713579	713759	H3K9ac	84
	...  

> `"intron_anno_K562_his"`  

	Chr	Start	End	INTRON	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2 H3K9ac	H3K9me3	H4K20me1
	chr1	12227	12612	DDX11L1-1-1	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	12721	13220	DDX11L1-1-2	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	15038	15795	WASH7P-1-9	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	15947	16606	WASH7P-1-8	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	18366	24737	WASH7P-1-2	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	24891	29320	WASH7P-1-1	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	700627	701708	LOC100288069-1-6	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	701767	703927	LOC100288069-1-5	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr1	703993	704876	LOC100288069-1-4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	...  

> `"pair_anno_K562_his"`  

	Intron_pair	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3 H4K20me1
	A1BG-1-4_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-5_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-5_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-6_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-6_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-6_5	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-7_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-7_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	A1BG-1-7_5	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	...  

> `"K562_dnase.list"`  

	DNaseI_HS  

> `"overlap_K562_dnase"`  

	chr1	705092	708355	LOC100288069-1-3	1.00	-	chr1	706265	706415	DNaseI_HS	150
	chr1	764484	776579	LINC01128-3-2	1.00	+	chr1	767140	767290	DNaseI_HS	150
	chr1	764484	783033	LINC01128-2-2	1.00	+	chr1	767140	767290	DNaseI_HS	150
	chr1	764484	783033	LINC01128-2-2	1.00	+	chr1	778220	778370	DNaseI_HS	150
	chr1	764484	787306	LINC01128-1-2	1.00	+	chr1	767140	767290	DNaseI_HS	150
	chr1	764484	787306	LINC01128-1-2	1.00	+	chr1	778220	778370	DNaseI_HS	150
	chr1	764484	787306	LINC01128-1-2	1.00	+	chr1	785040	785190	DNaseI_HS	150
	chr1	764484	787306	LINC01128-4-2	1.00	+	chr1	767140	767290	DNaseI_HS	150
	chr1	764484	787306	LINC01128-4-2	1.00	+	chr1	778220	778370	DNaseI_HS	150
	chr1	764484	787306	LINC01128-4-2	1.00	+	chr1	785040	785190	DNaseI_HS	150
	...  

> `"intron_anno_K562_dnase"`  

	Chr	Start	End	INTRON	DNaseI_HS
	chr1	12227	12612	DDX11L1-1-1	0.000
	chr1	12721	13220	DDX11L1-1-2	0.000
	chr1	15038	15795	WASH7P-1-9	0.000
	chr1	15947	16606	WASH7P-1-8	0.000
	chr1	18366	24737	WASH7P-1-2	0.000
	chr1	24891	29320	WASH7P-1-1	0.000
	chr1	700627	701708	LOC100288069-1-6	0.000
	chr1	701767	703927	LOC100288069-1-5	0.000
	chr1	703993	704876	LOC100288069-1-4	0.000
	...  

> `"pair_anno_K562_dnase"`  

	Intron_pair	DNaseI_HS
	A1BG-1-4_3	0.000
	A1BG-1-5_3	0.000
	A1BG-1-5_4	0.000
	A1BG-1-6_3	0.000
	A1BG-1-6_4	0.000
	A1BG-1-6_5	0.000
	A1BG-1-7_3	0.000
	A1BG-1-7_4	0.000
	A1BG-1-7_5	0.000
	...  


- ### Merge all features  

***CMD:***  

		merge_feature -t <cell_type>

		e.g.:
		merge_feature -t K562

Generate `"K562_anno_comb"`, e.g.:  

	Chr	Start	End	Intron_pair	Alu	DNaseI_HS	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3	H4K20me1
	chr19	58863648	58863921	A1BG-1-4_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58862756	58863921	A1BG-1-5_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58862756	58863053	A1BG-1-5_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58861735	58863921	A1BG-1-6_3	-1.985	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58861735	58863053	A1BG-1-6_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58861735	58862017	A1BG-1-6_5	-1.408	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58858718	58863921	A1BG-1-7_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58858718	58863053	A1BG-1-7_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	chr19	58858718	58862017	A1BG-1-7_5	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
	...  


>***Newly Updated:***  
>> Annotated intron (pair) using raw signal (reads count) from [ENCODE](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone)  `"bigWig"` file.  


- ### Data set preparation for model training, testing, validation and prediction

***CMD:***  

		prepare_train_set -t <cell_type> --circ <known_circ.bed> -R <ratio of negative VS positive> [ --sl < list of intron pair length ( sum of 2 introns ) to do stratified random sampling (space seperated 3 number, defult: 20000 30000 40000 ) > ]
		
		OR:
		
		prepare_train_set -t <cell_type> --circ <known_circ.bed (with expression (SRPBM) of 6 column)> ] --exp ( prepare data sets for expression prediction )

		e.g.:
		prepare_train_set -t K562 --circ K562_circ_overlap.bed -R 1 --sl 30000 50000 70000
		prepare_train_set -t K562 --circ K562_circ_overlap.bed --exp

Generate multiple files:   "K562_train", "K562_pred", "K562_circ_intron_pair", "K562_IP_part1", "K562_IP_part2", "K562_IP_part3", "K562_IP_part4"

OR:

"K562_exptrain", "K562_exp_pred"

`"K562_train"`, `"K562_pred"` used for modeling circRNAs expression status  
`"K562_exp_train"`, `"K562_exp_pred"` used for modeling circRNAs expression levels

> `"K562_train"` 

	Chr	Start	End	Intron_pair	Alu	DNaseI_HS	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3	H4K20me1	Type
	chr19	58861735	58863053	A1BG-1-6_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr10	52610424	52619745	A1CF-2-3_1	-0.464	0.000	0.000	0.335	0.062	0.000	0.744	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr10	52595833	52619745	A1CF-2-6_1	-0.850	0.000	0.000	0.482	0.000	0.000	0.716	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr10	52619601	52619745	A1CF-6-4_3	-0.148	0.000	0.000	0.000	0.000	0.000	0.300	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr12	9258831	9265132	A2M-1-10_2	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr12	9256834	9266139	A2M-1-11_1	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr12	9246060	9262930	A2M-1-18_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr12	9242497	9246175	A2M-1-21_17	0.000	0.000	0.000	0.000	0.000	0.000	1.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	chr12	9231839	9247680	A2M-1-25_16	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	F
	...  

> `"K562_pred"`  

	Chr	Start	End	Intron_pair	Alu	DNaseI_HS	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3	H4K20me1	Type
	chr19	58863648	58863921	A1BG-1-4_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58862756	58863921	A1BG-1-5_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58862756	58863053	A1BG-1-5_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58861735	58863921	A1BG-1-6_3	-1.985	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58861735	58862017	A1BG-1-6_5	-1.408	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58858718	58863921	A1BG-1-7_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58858718	58863053	A1BG-1-7_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58858718	58862017	A1BG-1-7_5	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	P
	chr19	58858718	58859006	A1BG-1-7_6	-2.130	0.000	0.000	0.000	0.519	0.000	0.000	0.768	0.000	0.729	0.533	0.000	0.279	0.000	0.000	P
	... 

> `"K562_exp_train"`

	Chr	Start	End	Intron_pair	Alu	DNaseI_HS	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3	H4K20me1	SRPBM
	chr12	125558421	125576069	AACS-1-1_5	0.105	0.015	0.032	1.000	0.041	0.000	0.227	1.000	0.111	0.053	0.000	0.0000.000	0.000	1.000	0.00936076724592655
	chr5	178199429	178203277	AACSP1-1-8_4	-6.708	0.000	0.000	0.000	0.000	0.000	0.000	0.454	0.000	0.000	0.000	0.0000.000	0.000	0.000	0.0561646034755593
	chr5	178199429	178203277	AACSP1-2-8_4	-6.711	0.000	0.000	0.000	0.000	0.000	0.000	0.454	0.000	0.000	0.000	0.0000.000	0.000	0.000	0.0561646034755593
	chr9	99413671	99413994	AAED1-1-4_2	0.000	0.000	0.000	1.000	0.000	0.000	0.000	1.000	0.330	0.000	0.000	0.0000.000	0.221	0.000	0.0280823017377796
	chr15	67528316	67529158	AAGAB-1-4_1	-1.888	0.000	0.000	1.000	0.000	0.000	0.000	1.000	0.000	0.000	0.000	1.0000.000	0.996	0.749	0.0374430689837062
	chr15	67524151	67529158	AAGAB-1-5_1	0.606	0.000	0.000	1.000	0.000	0.000	0.000	1.000	0.154	0.000	0.000	1.0000.242	0.938	0.492	0.636532172723005
	chr15	67500899	67501882	AAGAB-1-7_5	-1.614	0.000	0.000	1.000	0.000	0.000	0.000	1.000	0.168	0.000	0.000	0.9690.562	0.941	0.657	0.0187215344918531
	chr15	67528316	67529158	AAGAB-2-4_1	-1.851	0.000	0.000	1.000	0.000	0.000	0.000	1.000	0.000	0.000	0.000	1.0000.000	0.996	0.738	0.0374430689837062
	chr15	67524151	67529158	AAGAB-2-5_1	0.594	0.000	0.000	1.000	0.000	0.000	0.000	1.000	0.153	0.000	0.000	1.0000.262	0.938	0.485	0.636532172723005
	... 

> `"K562_exp_pred"`

	Chr	Start	End	Intron_pair	Alu	DNaseI_HS	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3	H4K20me1	SRPBM
	chr19	58863648	58863921	A1BG-1-4_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58862756	58863921	A1BG-1-5_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58862756	58863053	A1BG-1-5_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58861735	58863921	A1BG-1-6_3	-1.985	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58861735	58863053	A1BG-1-6_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58861735	58862017	A1BG-1-6_5	-1.408	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58858718	58863921	A1BG-1-7_3	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58858718	58863053	A1BG-1-7_4	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	chr19	58858718	58862017	A1BG-1-7_5	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.0000.000	0.000	0.000	EP
	... 


>***NOTE:***   
>> `known_circ_overlap.bed`: [Known circRNAs](https://github.com/johnlcd/CIRCScan/tree/master/data/known_circ) from RNA-seq data by overlapping the detection of multiple tools
<br>


### 2. Model traing and prediction

- ### Complete process for predicting circRNAs expression status 

***a). Model training ( Used for obtaining rank of importance )***  

***CMD:***  

		circscan --train -t <cell_type> -m <model> -s/-r <seed>/<run (1-5)> -n <cores>
		# "-n": used for models training by parellel
		# "-s": used for random sample training set ("seed" mode, multiple traning and reproducibility )
		# "-r": used for 5 independent CV run ("run" mode)

		e.g.:
		circscan --train -t K562 -m rf -s 111 -n 8
		or:
		circscan --train -t K562 -m rf -r 1 -n 8

Generate models and R data file `"K562_rf_train.RData"`, log file `"K562_rf_train.out"` with model evaluation  


***b). Feature selection***  

***CMD:***  

		circscan --fs -t <cell_type> -m <model> -n <cores> -l <all/feature_number_list> < --auc / --f1 (referenece index) > [ -pt (type of prediction) raw/prob (probabilities, default) ]	
		# "-n": used for models training by parellel
		# "-l": list of feature number for feature selection. If value is "all", then run feature selection with feature number from 1 to all, if is a list of feature number ( comma separsted ), for example: 1,2,3,4,5,10,15, then run feature selection with feature number you provide
		# "--auc / --f1": referenece index to evaluate model performance
		# "--pt": type of prediction, defult is 'prob' (probabilities), 'raw' is used for models without probabilities

		e.g.:
		circscan --fs -t K562 -m rf -n 8 -l all --auc

Generate R data file `"K562_rf_FS.RData"` of feature selection and log file `"K562_rf_FS.out"`, and `"K562_rf_perf_test.txt"` with results of feature selection	( Feature number with highest *F1* score )  

>***NOTE:***  
>> Feature selection is required to generate and select the best model for circRNAs prediction.  


***c). Model performance evaluation (testing data)***  

***CMD:***  

		circscan --eval -t <cell_type> -m <model> -n <cores>
		# "-n": used for models training by parellel

		e.g.:
		circscan --eval -t K562 -m rf -n 8

Generate R data file `"K562_rf_test_perf_eval.RData"` of model performance evaluation and log file `"eval_K562_rf_perf.out"`, and `"K562_rf_eval_test_perf"` with results of model performance in testing data  


***d). CircRNAs prediction and annotation***  

***CMD:***  

		circscan --pred -t <cell_type> -m <model> -n <cores>
		# "-n": used for models training by parellel

		e.g.:
		circscan --pred -t K562 -m rf -n 8

Generate predicted anaotated circRNAs file `"K562_rf_pred_true.bed"`  

<br>

***Mode 2). Model training, feature selection, and validation (known circRNAs) (Optional)***  

***CMD:***  

		circscan --exp -t <cell_type> -m <model> -n <cores> -sf < all/select_fea_list (comma separated)> -l <known circRNA_intron_FIP_list>	
		# "-n": used for models training by parellel
		# "-sf": list of selected feature list (comma separated, default: all)
		# "-l": circBase circRNAs FIP list file

		e.g.:
		circscan --exp -t K562 -m rf -n 8 -sf all/Alu,H3K36me3,... -l GM12878_circ_FIP.list


<br><br>

- ### Complete process for predicting circRNAs expression levels 

***a). Model training and feature selection( Used for obtaining rank of importance and select the best features)***  

***CMD:***  

		circscan --exp-fs -t <cell_type> -m <model> -s/-r <seed>/<run (1-5)> -n <cores> [ --RM (remove outlier) ]
		# "-n": used for models training by parellel
		# "-s": used for random sample training set ("seed" mode, multiple traning and reproducibility )
		# "-r": used for 5 independent CV run ("run" mode)
		# "--RM": whether remove outlier data points

		e.g.:
		circscan --exp-fs -t K562 -m rf -s 111 -n 8
		or:
		circscan --exp-fs -t K562 -m rf -r 1 -n 8

Generate models and R data file `"K562_rf_FS_exp.RData"`, log file `"K562_rf_FS_exp.out"`, result of model performance `"GM12878_rf_perf_test_reg.txt"` with model evaluation of feature selection  


***b). Model performance evaluation (testing data)***  

***CMD:***  

		circscan --exp-eval -t <cell_type> -m <model> -n <cores>
		# "-n": used for models training by parellel

		e.g.:
		circscan --exp-eval -t K562 -m rf -n 8

Generate models and R data file `"K562_rf_pred_exp_all.RData"`, log file `"K562_rf_eval_cor.out"`, result of model performance `"K562_rf_eval_test_exp_perf"` with model evaluation in testing data, results of predicted expression levels `"K562_rf_obs_pred_exp_test.all"`, `"K562_rf_obs_pred_exp_train.all"`  


***c). Predict expression***  

***CMD:***  

		circscan --exp-pred -t <cell_type> -m <model> -n <cores> -sf < all/select_fea_list (comma separated)> [ --RM (remove outlier) ]
		# "-n": used for models training by parellel
		# "-sf": select fea list (according to results of feature selection)
		# "--RM": whether remove outlier data points

		e.g.:
		circscan --exp-pred -t K562 -m rf -n 8 -sf all

Generate predicted circRNAs expression file `"K562_rf_train_pred_exp"`  

<br><br>


>## *All data are avaliable in `"$PKG_DIR/sample"` for performing and replicating the works in this Document.*

