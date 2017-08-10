
# CIRCScan  
Tools for predicting circRNAs with `EPIGENETIC` features by `Machine Learning`  
<br>

## Introduction
Circular RNAs (circRNAs) are an abundant class of noncoding RNAs with the widespread, cell/tissue specific pattern. This tool \(pipeline\), 
		 **CIRCScan**, is used to predict circRNAs in a cell/tissue specific manner by machine learning based on epigenetic features.  
<br>
  
## License
This software is distributed under the terms of GPL 2.0  
<br>

## Source
[https://github.com/johnlcd/CIRCScan](https://github.com/johnlcd/CIRCScan)  
<br>

## Contact
### Author
**Jia-Bin Chen**, **Shan-Shan Dong**, **Shi Yao**, **Yan Guo**, **Tie-Lin Yang**  
Key Laboratory of Biomedical Information Engineering of Ministry of Education, School of Life Science and Technology, Xi'an Jiaotong University, Xi'an, Shaanxi Province, 710049, P. R. China  
[yangtielin@mail.xjtu.edu.cn](yangtielin@mail.xjtu.edu.cn)
<br>

## Maintainer
**Jia-Bin Chen**  
You can contact [johnlcd@stu.xjtu.edu.cn](johnlcd@stu.xjtu.edu.cn) when you have any questions, suggestions, comments, etc.  
Please describe in details, and attach your command line and log messages if possible.
<br>

## Requiremnets
- **bedtools** \(v2.25.0\)
- **Python** \(recommended: python2.7\)
- **R** \(>= 3.2.4\)
		R packages: caret, ggplot, doParallel, ROCR, etc. (Dependent packages for different models) 
		
		Check the output file ".out" to validate which package is required if got an error
<br>

## Directory catalog
<br>

- **bin**  	
	- **anno**  
		- anno_bedpe.py  
		- anno_intron.py
		- bt_overlap
		- comb_pair_anno.py
	- anno_pair
	- circpred
	- fast_model
	- merge_feature
	- **model**
		- Circ_pred.R
		- fast_model.R
		- feature_selection.R
		- make_set.py
		- Model_train.R
		- sep_intron_true_unknown.py
	- prepare_train_set
- **data**
	- intron_intron-pairs.tgz
	- **pred_circ_bycell**
		- A549_pred_circ.bed.gz
		- GM12878_pred_circ.bed.gz
		- H1-hESC_pred_circ.bed.gz
		- HeLa-S3_pred_circ.bed.gz
		- HepG2_pred_circ.bed.gz
		- HMEC_pred_circ.bed.gz
		- HSMM_pred_circ.bed.gz
		- HUVEC_pred_circ.bed.gz
		- K562_pred_circ.bed.gz
		- NHEK_pred_circ.bed.gz
		- NHLF_pred_circ.bed.gz
	- **raw_data**
		- 4DGenome_HomoSapiens_hg19.txt.gz
		- ChromHMM.txt.gz
		- DNaseI.txt.gz
		- GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_new.txt.gz
		- GSE63525_HeLa_HiCCUPS_looplist_new.txt.gz
		- GSE63525_HMEC_HiCCUPS_looplist_new.txt.gz
		- GSE63525_K562_HiCCUPS_looplist_new.txt.gz
		- GSE63525_NHEK_HiCCUPS_looplist_new.txt.gz
		- Histone_part1.txt.gz
		- Histone_part2.txt.gz
		- RBP.txt.gz
		- select_cell.list
- **info**
	- models_ALL.list
	- models_classification.list
	- models_test.list
	- models_test.txt
- README.md
- **sample**
	- **anno**
		- intron_anno_K562_dnase.gz
		- ...
	- **fast_process**
		- K562_pred.gz
		- K562_train.gz
	- **feature**
		- **all_feature**
		- feature_select.list
		- K562_dnase.bed.gz
		- ...
	- K562_anno_comb.gz
	- K562_circ.bed.gz
	- K562_circ_true.gz
	- K562_unknown.gz
	- **model**
		- K562_pred.gz
		- K562_train.gz
<br>


## Runing preparation
### Set environment variables  
		export PKG_DIR=/path/to/tool_package
		export PATH=$PKG_DIR/bin:PKG_DIR/bin/anno:$PKG_DIR/bin/model:$PATH
<br>
### Unzip data files  
		cd $PKG_DIR/data
		tar -zxvf intron_intron-pairs.tgz
		cd $PKG_DIR/data/raw_data
		gunzip *.gz
		cat Histone_part1.txt Histone_part2.txt > Histone.txt
<br>

## Work flow


# PART 1: Data preparation and feature generation

	Prepare intron and intron pair data file.

	CMD:	cd $PKG_DIR/data; tar -zxvf intron_intron-pairs.tgz


## 1. Feature generation and annotation:

## RBPs, Histoen modifications, ChromHMM, DNaseI HS && other features of region ( Feature types of "bed" format )

	Make feature list, and overlap intron with feature, annotate intron by features, then combine intron annotation to pair ( "anno_pair" )

	CMD:	$PKG_DIR/bin/anno_pair -t <cell_type> -f <feature (rbp, his, hmm, dnase ...)> [ --is (Ignor strands) ] --bed <feature.bed>

	 generate 4 files:

	 a. "cell_feature.list"
#   CBP_(sc-369)
#   CBX2
#   CBX3_(SC-101004)
#   CBX8
#   CHD1_(A301-218A)
#   CHD4_Mi2
#   CHD7_(A301-223A-1)
#   CTCF
#   EZH2_(39875)
#   H2A.Z

	 b. "ovelapped_cell_feature"
#	chrY	16636816	16733888	NLGN4Y-1-1	1.00	+	chrY	16636680	16636830	dnase	14
#	chrY	16734471	16834996	NLGN4Y-1-2	1.00	+	chrY	16740180	16740330	dnase	150
#	chrY	16835149	16936067	NLGN4Y-1-3	1.00	+	chrY	16875640	16875790	dnase	150
#	chrY	16936253	16941609	NLGN4Y-1-4	1.00	+	.	-1	-1		.	0
#	chrY	16942399	16952292	NLGN4Y-1-5	1.00	+	.	-1	-1		.	0
#	chr7	28452536	28527792	CREB5-1-1	1.00	+	chr7	28512045	28512195	dnase	150
#	chr7	28527864	28534523	CREB5-1-2	1.00	+	.	-1	-1		.	0
#	chr7	28534617	28547233	CREB5-1-3	1.00	+	.	-1	-1		.	0
#	chr7	28547355	28609982	CREB5-1-4	1.00	+	chr7	28568865	28569015	dnase	150
#	chr7	28610155	28758369	CREB5-1-5	1.00	+	chr7	28666365	28666515	dnase	150

	c. "intron_anno_cell_feature"
#	Chr	Start	End	INTRON	DNaseI_HS
#	chrY	16636816	16733888	NLGN4Y-1-1	1
#	chrY	16734471	16834996	NLGN4Y-1-2	1
#	chrY	16835149	16936067	NLGN4Y-1-3	1
#	chrY	16936253	16941609	NLGN4Y-1-4	0
#	chrY	16942399	16952292	NLGN4Y-1-5	0
#	chr7	28452536	28527792	CREB5-1-1	1
#	chr7	28527864	28534523	CREB5-1-2	0
#	chr7	28534617	28547233	CREB5-1-3	0
#	chr7	28547355	28609982	CREB5-1-4	1

	d. "pair_anno_cell_feature"
#	Intron_pair	DNase_HS
#	A1BG-1-2_1	0
#	A1BG-1-3_1	1
#	A1BG-1-3_2	1
#	A1BG-1-4_1	0
#	A1BG-1-4_2	0
#	A1BG-1-4_3	1
#	A1BG-1-5_1	1
#	A1BG-1-5_2	1
#	A1BG-1-5_3	2

## Hi-C/pairs data ( Feature types of "bedpe" format )
	
	Overlap intron pairs with Hi-C pairs, annotate intron pairs ( "anno_pair" )

	CMD:	$PKG_DIR/bin/anno_pair -t <cell_type> -f <paired feature (hic ...)> [ --is (Ignor strands) ] --bedpe <feature.bedpe>

	generate 2 files:

	a. "overlapped_cell_hic"
#	chr1	1238661	1239465	chr1	1239523	1243148	ACAP3-1-2_1	1.00	-	-	chr1	1239001	1241001	chr1	1242947	1244947	Hi-C
#	chr1	1577362	1647784	chr1	1647917	1650766	CDK11B-5-5_4	1.00	-	-	chr1	1585000	1590000	chr1	1645000	1650000	Hi-C
#	chr1	1577869	1643702	chr1	1643839	1647784	CDK11B-2-6_5	1.00	-	-	chr1	1585000	1590000	chr1	1645000	1650000	Hi-C
#	chr1	1577869	1643702	chr1	1647917	1650766	CDK11B-2-6_4	1.00	-	-	chr1	1585000	1590000	chr1	1645000	1650000	Hi-C
#	chr1	1577869	1643702	chr1	1647917	1650796	CDK11B-3-6_4	1.00	-	-	chr1	1585000	1590000	chr1	1645000	1650000	Hi-C
#	chr1	1577869	1647784	chr1	1647917	1650766	CDK11B-1-5_4	1.00	-	-	chr1	1585000	1590000	chr1	1645000	1650000	Hi-C
#	chr1	1688749	1690539	chr1	1696885	1711343	NADK-2-5_1	1.00	-	-	chr1	1705401	1707401	chr1	1689081	1691081	Hi-C
#	chr1	1688749	1690539	chr1	1696885	1711343	NADK-2-5_1	1.00	-	-	chr1	1707001	1709001	chr1	1689081	1691081	Hi-C
#	chr1	1688749	1693390	chr1	1696885	1709727	NADK-4-3_1	1.00	-	-	chr1	1705401	1707401	chr1	1689081	1691081	Hi-C
#	chr1	1688749	1693390	chr1	1696885	1709727	NADK-4-3_1	1.00	-	-	chr1	1707001	1709001	chr1	1689081	1691081	Hi-C

	b. "pair_anno_cell_hic"
#	Intron_pair	Hi-C
#	DDX11L1-1-1_2	0
#	WASH7P-1-10_9	0
#	WASH7P-1-10_8	0
#	WASH7P-1-10_7	0
#	WASH7P-1-10_6	0
#	WASH7P-1-10_5	0
#	WASH7P-1-10_4	0
#	WASH7P-1-10_3	0
#	WASH7P-1-10_2	0


## 2. Merge all features

	CMD:	merge_feature -t <cell_type>

	generate "cell_anno_comb", e.g:
#	Chr	Start	End	Intron_pair	DNaseI-HS	Hi-C	CTCF	EZH2_(39875)	H2A.Z	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me2	H3K4me3	H3K79me2	H3K9ac	H3K9me3	H4K20me1
#	chr19	58864657	58864693	A1BG-1-2_1	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58864293	58864693	A1BG-1-3_1	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58864293	58864563	A1BG-1-3_2	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58863648	58864693	A1BG-1-4_1	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58863648	58864563	A1BG-1-4_2	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58863648	58863921	A1BG-1-4_3	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58862756	58864693	A1BG-1-5_1	0	0	0	0	0	0	0	0	0	0	0	0    0
#	chr19	58862756	58864563	A1BG-1-5_2	0	0	0	0	0	0	0	0	0	0	0	0    0


## 3. Data set preparation ( Preparation for training set and intron pairs for prediction )

	CMD:	prepare_train_set -t <cell_type> --circ <known_circ.bed> -R <ratio of negative VS positive> \
	[--re-sample (re-sample and overwrite data set if exsist)]	
	# "known_circ.bed":reported circRNA BED file (positive set)

	generate 4 files:	"cell_circ_true", "cell_unknown", "cell_train", "cell_pred"



# PART 2: Model traing and prediction


## 2-1. Fast process for model training and prediction ( Train and predict with H3K36me3 and H3K79me2 )

# NOTE: 
#	Data matrix for model training ( "cell_train" ) and prediction ( "cell_pred" ) contain only H3K36me3 and H3K79me2 but not other features, e.g
#	Chr	Start	End	Intron_pair	H3K36me3	H3K79me2	Type
#	chr12	53702065	53714476	AAAS-1-12_1	2	1	T
#	chr12	53701835	53709210	AAAS-1-13_3	2	0	T
#	chr12	125587224	125603311	AACS-1-5_10	2	1	T
#	chr12	125587224	125591814	AACS-1-5_8	2	1	T
#	chr12	125587224	125599103	AACS-1-5_9	2	1	T
#	chr12	125599022	125603311	AACS-1-8_10	2	0	T
#	chr15	67524151	67529158	AAGAB-1-5_1	2	2	T
#	chr15	67496381	67529158	AAGAB-1-8_1	2	2	T
#	chr15	67524151	67529158	AAGAB-2-5_1	2	2	T
	
	CMD:	fast_model -t <cell_type> -m <model> -n <cores>
	# "-n": used for models training by parellel
	
	generate models and R data file ( "cell_model_fast_model.RData" ), output file ( "cell_model_fast_model.out" ) and predicted circRNA file "cell_model_pred_true.bed"


## 2-2. Complete process for model training prediction ( With feature selection )

### 1. Model training ( Used for obtaining rank of importance )

	CMD:	circpred --train -t <cell_type> -m <model> -s <seed> -n <cores>
	# "-n": used for models training by parellel
	# "seed": used for random sample training set (multiple traning and reproducibility)

	generate models and R data file ( "cell_model_fast_model.RData" ), output file ( "cell_model_fast_model.out" ) with model evaluation


### 2. Feature selection
	
	CMD:	circpred --fs -t <cell_type> -m <model> -n <cores> -l <all/feature_number_list>	
	# "-n": used for models training by parellel
	# "-l": list of feature number for feature selection. If value is "all", then run feature selection with feature number from 1 to all, if is a list of feature number (comma separsted), for example: 1,2,3,4,5,10,15, then run feature selection with feature number you provided

	generate R data file "cell_model_FS.RData" of feature selection and output file ( "cell_model_FS.out" ) with results of feature selection	# Select feature number with best performance (F1 score)

## NOTE:
# Feature selection is required to generate and select the best model for circRNA prediction.


### 3. circRNA prediction and annotation with trained model

	CMD:	circpred --pred -t <cell_type> -m <model> -n <cores>
	# "-n": used for models training by parellel

	generate predicted anaotated circRNA file "cell_model_pred_true.bed"

-----------------------------------------------------------------


=================================================================

### END

=================================================================

