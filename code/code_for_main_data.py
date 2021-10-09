import numpy as np

###### ### 1), 2), 3) --> our feature, 4) --> control feature
### 1) conservation of RC
### 2) conservation of RE
### 3) expression conservation
### 4) conservation of RC (single gene)

input_path = "../data/"

### label - PS score
ps_dict = dict()
f = open(input_path+"PSscores.txt")
head_list = f.readline()
for line in f.readlines():
	line = line.strip().split("\t") 
	gene = line[0]
	ps = line[4]
	pg = line[5]
	ps_dict[gene] = (ps, pg)
f.close()

### features - 1) conservation of RC
rc_dict = dict()
f = open(input_path+"link_TF_conservation_score_50_GO_bp.txt")
head_list = f.readline()
for line in f.readlines():
	line = line.strip().split("\t")
	gene = line[0]
	ji = float(line[1])
	oc = float(line[2])
	rc_dict[gene] = (ji, oc)
f.close()

### features - 2) conservation of RE
re_pro_dict, re_en_dict = dict(), dict()
f = open(input_path+"RE_score_promoter.txt")
head_list = f.readline()
for line in f.readlines():
	line = line.strip().split("\t")
	gene = line[0]
	if line[1] != "nan":
		conserved = float(line[1])
		specific = float(line[2])
		re_pro_dict[gene] = (conserved, specific)
f.close()

f = open(input_path+"RE_score_enhancer.txt")
head_list = f.readline()
for line in f.readlines():
	line = line.strip().split("\t")
	gene = line[0]
	if line[1] != "nan":
		conserved = float(line[1])
		specific = float(line[2])
		re_en_dict[gene] = (conserved, specific)
f.close()

### features - 3) expression conservation --> pearson  
expr_dict = dict()
f = open(input_path+"module_expression_scores_GO_bp.txt")
head_list = f.readline()
for line in f.readlines():
	line = line.strip().split("\t")
	gene = line[0]
	if line[3] != "nan":
		expr_conv = float(line[3])
		expr_dict[gene] = expr_conv
f.close()

### feature - 4) conservation of RC single
rc_sing_dict = dict()
f = open(input_path+"single_gene_score_conservation_of_RC.txt")
head_list = f.readline()
for line in f.readlines():
    line = line.strip().split("\t")
    gene = line[0]
    ji = float(line[1])
    oc = float(line[2])
    rc_sing_dict[gene] = (ji, oc)
f.close()


###### graph presentation
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

labels = ["LPG", "HPG"]
### feature 1)
rc_lpg_score, rc_hpg_score = [], []
for gene in ps_dict:
    if gene in rc_dict:
        if ps_dict[gene][1] == "LPG":
            rc_lpg_score.append(rc_dict[gene][0])
        if ps_dict[gene][1] == "HPG":
            rc_hpg_score.append(rc_dict[gene][0])
rc_data = [rc_lpg_score, rc_hpg_score]

significance = mannwhitneyu(rc_data[0], rc_data[1])
pvalue = significance[1]

plt.title(pvalue)
plt.ylabel("Conservation of RC")
plt.boxplot(rc_data, labels = labels)
plt.show()


### feature 2) -> enhancer
re_lpg_score, re_hpg_score = [], []
for gene in ps_dict:
    if gene in re_en_dict:
        if ps_dict[gene][1] == "LPG":
            re_lpg_score.append(re_en_dict[gene][1])
        if ps_dict[gene][1] == "HPG":
            re_hpg_score.append(re_en_dict[gene][1])
re_data = [re_lpg_score, re_hpg_score]

significance = mannwhitneyu(re_data[0], re_data[1])
pvalue = significance[1]

plt.title(pvalue)
plt.ylabel("Species-specific RE\nin regulatory network")
plt.boxplot(re_data, labels = labels)
plt.show()


### feature 3)
expr_lpg_score, expr_hpg_score = [], []
for gene in ps_dict:
    if gene in expr_dict:
        if ps_dict[gene][1] == "LPG":
            expr_lpg_score.append(expr_dict[gene])
        if ps_dict[gene][1] == "HPG":
            expr_hpg_score.append(expr_dict[gene])
expr_data = [expr_lpg_score, expr_hpg_score]

significance = mannwhitneyu(expr_data[0], expr_data[1])
pvalue = significance[1]

plt.title(pvalue)
plt.ylabel("Expression conservation")
plt.boxplot(expr_data, labels = labels)
plt.show()


### feature 4)
sing_lpg_score, sing_hpg_score = [], []
for gene in ps_dict:
    if gene in rc_sing_dict:
        if ps_dict[gene][1] == "LPG":
            sing_lpg_score.append(rc_sing_dict[gene][0])
        if ps_dict[gene][1] == "HPG":
            sing_hpg_score.append(rc_sing_dict[gene][0])
sing_data = [sing_lpg_score, sing_hpg_score]

significance = mannwhitneyu(sing_data[0], sing_data[1])
pvalue = significance[1]

plt.title(pvalue)
plt.ylabel("Conservation of RC")
plt.boxplot(sing_data, labels = labels)
plt.show()





