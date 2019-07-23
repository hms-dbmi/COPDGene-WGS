
# Load Hail & packages
import hail as hl
hl.init()

from hail.plot import show
hl.plot.output_notebook()
from hail.expr import aggregators
from hail.expr.expressions import *
from hail.expr.expressions import Expression
from hail.typecheck import *
from hail import Table
from bokeh.plotting import figure, output_file, show, save
from pprint import pprint
from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import Span, ColumnDataSource
output_notebook()
from bokeh.models import *
import numpy as np
from numpy import array, empty
import pandas as pd
from pprint import pprint
from math import log, isnan, log10
from itertools import cycle
import os, time, sys

#Load COPD freeze 4 matrix and annotation table
mt = hl.read_matrix_table("path_to_vcf")
t  = hl.import_table("spath_to_table", impute = True).key_by('Samples')
mt = mt.annotate_cols(**t[mt.s])

# Speed up by dividing MT file into partitions that match the cloud cluster
# Number of repartitions depends on CPU per node and number of nodes.
CPU = 16
nodes = 20
mt = mt.repartition(4 * CPU * nodes)

# Remove monomorphic data (homozygous to reference), keeping alternate variants
# Keep genotype of variant when it appears more than 1
mt = mt.filter_rows(hl.agg.fraction(mt.GT.is_hom_ref()) < 1)




######## 1. SANITY CHECK
mt.describe()

# Show all possible alternate alleles
mt.row.alleles.show(5)

# Calculate all possible GT calls over entire dataset
possibleGT = mt.aggregate_entries(hl.agg.collect_as_set(mt.GT))
print(len(possibleGT))
possibleGT

#Show unique possible calls occuring in entire dataset
unique_allelecalls = mt.aggregate_rows(hl.struct(
    ref = hl.agg.collect_as_set(mt.alleles[0]),
    alt = hl.agg.collect_as_set(mt.alleles[1])))
pprint(unique_allelecalls)






######## 2. QUALITY CONTROL VARIANTS
######## 2.1 Optional: Pruning in Linkage disequilibrium

# Function works only on biallelic data
biallelic_mt = mt.filter_rows(hl.len(mt.alleles) == 2)
# Prune with window size of 44 kb for Caucasian and 22 kb for African
pruned_t = hl.ld_prune(mt.GT, r2 = 0.8, bp_window_size = 44000,  memory_per_core = 128)
mt = mt.filter_rows(hl.is_defined(pruned_t[mt.row_key]))

######### 2.2 Hardy-weinberg equilibrium
# HWE separately for each ethnic group, on only controls.
mt_NHW = mt.filter_cols(mt.Race == "Caucasian")
mt_NHW = mt_NHW.annotate_rows(hwe_ctrl = hl.agg.filter(mt_NHW.Affection == 'Control', hl.agg.hardy_weinberg_test(mt_NHW.GT)))
mt_NHW = mt_NHW.filter_rows(mt_NHW.hwe_ctrl.p_value > 10**-5)
mt_AA = mt.filter_cols(mt.Race == "African American")
mt_AA = mt_AA.annotate_rows(hwe_ctrl = hl.agg.filter(mt_AA.Affection == 'Control', hl.agg.hardy_weinberg_test(mt_AA.GT)))
mt_AA = mt_AA.filter_rows(mt_AA.hwe_ctrl.p_value > 10**-5)
# Merge of NHW and AA by columns (samples)
mt = mt_AA.union_cols(mt_NHW)

######## 2.3 Check for missingness
# Calculate variant statistics
mt = hl.variant_qc(mt)
# Keep variants with call rate (called/(non-called +called)) > 0.98.
mt = mt.filter_rows(mt.variant_qc.call_rate > 0.98)

######## 2.5 Remove samples with identical ID classifiers
# Rename the second duplicate ID <name> as <name>_1
t = hl.rename_duplicates(mt).cols()
t.key_by("s")
mt.key_cols_by("s")
mt = mt.filter_cols(hl.is_defined(t[mt.s]))

####### 2.6 Remove wrongly assigned variant types
# calculate rate of allele depth
ab = mt.AD[1] / hl.sum(mt.AD)
# compare allele depth rate of genotype (GT) to labels homozygous to reference, heterozygous and homozygous variant
filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                        (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                        (mt.GT.is_hom_var() & (ab >= 0.9)))
mt = mt.filter_entries(filter_condition_ab)

######## 2.7 Remove non-autosomes
# Note: This step needs to be executed after the sex check for samples, because impute_sex() function needs the sex chromosome. You can find this step after the sample filtering steps.

######## 2.8 Optional: Keep SNPs only (remove indels, CNVs etc.)
mt_snp = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))  #ref (0) to alternates(1)
#Check SNPs: unique possible reference (ref) and alternate allele calls (alt) from entire dataset (all samples)
unique_allelecalls = mt_snp.aggregate_rows(hl.struct(
    ref = hl.agg.collect_as_set(mt_snp.alleles[0]),
    alt = hl.agg.collect_as_set(mt_snp.alleles[1])))
pprint(unique_allelecalls)

#Check SNPs: shows all lenghts of vectors with possible allels (including ref, alternate)
a = mt_snp.aggregate_rows(hl.agg.collect_as_set(hl.len(mt_snp.alleles)))
pprint(a)
mt_AF = mt.filter_rows(mt.variant_qc.AF[1] >= 0.01)






######## 3. QUALITY CONTROL SAMPLES
######## 3.1 Filter samples for outliers more than (6 * SD) from mean (Part 1)
# Calculate sample statistics
mt = hl.sample_qc(mt)
# Calculate statistics on sample statistics
stats_singleton   = mt.aggregate_cols(hl.agg.stats(mt.sample_qc.n_singleton))
stats_ti_tv       = mt.aggregate_cols(hl.agg.stats(mt.sample_qc.r_ti_tv))
stats_het_hom_var = mt.aggregate_cols(hl.agg.stats(mt.sample_qc.r_het_hom_var))
stats_het         = mt.aggregate_cols(hl.agg.stats(mt.sample_qc.n_het))

######## 3.2 Sex check on chromosome X (inbreeding coefficient)
# Determine sex from GT calls in sex chromosomes
t = hl.impute_sex(mt.GT)
# Only keep those where genetic sex matches self-reported Sex
mt = mt.filter_cols(t[mt.s].is_female == mt.is_female)

######## 3.3 Check for genetic relationship / "duplicates"
# Calculate identity-by-descent matrix
mt_relatedness = hl.identity_by_descent(mt)
# keep pairs of samples with PI_HAT in [0.2, 1] using MAF computed from the dataset itself in row field panel_maf.
t_ibd = relatedness.filter(relatedness.ibd.PI_HAT > 0.2)
t_ibd.key_by('i')
mt.key_cols_by("s")
#Collect the IDs of the related samples in t_ibd
ibd_idx = t_ibd.aggregate(hl.agg.collect_as_set(t_ibd.i))
mt_ibd = mt.filter_cols(hl.is_defined(ibd_idx))

######### 3.3 Filter samples for outliers more than (6 * SD) from mean (Part 2)
# Number of singletons
mt = mt.filter_cols(mt.sample_qc.n_singleton < (stats_singleton.mean + (6 * stats_singleton.stdev)))
mt = mt.filter_cols(mt.sample_qc.n_singleton > (stats_singleton.mean - (6 * stats_singleton.stdev)))
#Ti/Tv ratio
mt = mt.filter_cols(mt.sample_qc.r_ti_tv < (stats_ti_tv.mean + (6 * stats_ti_tv.stdev)))
mt = mt.filter_cols(mt.sample_qc.r_ti_tv > (stats_ti_tv.mean - (6 * stats_ti_tv.stdev)))
#Het/hom ratio
mt = mt.filter_cols(mt.sample_qc.r_het_hom_var < (stats_het_hom_var.mean + (6 * stats_het_hom_var.stdev)))
mt = mt.filter_cols(mt.sample_qc.r_het_hom_var > (stats_het_hom_var.mean - (6 * stats_het_hom_var.stdev)))
#Number of heterozygous calls
mt = mt.filter_cols(mt.sample_qc.n_het < (stats_het.mean + (6 * stats_het.stdev)))
mt = mt.filter_cols(mt.sample_qc.n_het > (stats_het.mean - (6 * stats_het.stdev)))

######## 3.4 Remove non-autosomes(X, Y and MT DNA)
mt = mt.filter_rows(mt.locus.in_autosome())




######## 4. BASELINE CHARACTERISTICS QC-FILTERED DATA
# Summary on number of SNPs, indels and variants per chromosomes
hl.summarize_variants(mt)

#Partition data into cases (mt_case) and controls (mt_ctrl)
mt_case = mt.filter_cols(mt.Affection == 'Case')
mt_ctrl = mt.filter_cols(mt.Affection == 'Control')

#Calculate subject statistics
print('Age of cases =', mt_case.aggregate_cols(hl.agg.stats(mt_case.Age)))
print('Age of controls =', mt_ctrl.aggregate_cols(hl.agg.stats(mt_ctrl.Age)))

print('#Individuals of Cases:', mt_case.aggregate_cols(hl.agg.counter(mt_case.Race)))
print('#Individuals of Controls:', mt_ctrl.aggregate_cols(hl.agg.counter(mt_ctrl.Race)))

print('Gender  Cases:', mt_case.aggregate_cols(hl.agg.counter(mt_case.is_female)))
print('Gender Controls:', mt_ctrl.aggregate_cols(hl.agg.counter(mt_ctrl.is_female)))

print(mt_case.aggregate_cols(hl.agg.stats(mt_case.FEV1perc)))
print(mt_ctrl.aggregate_cols(hl.agg.stats(mt_ctrl.FEV1perc)))

print(mt_case.aggregate_cols(hl.agg.stats(mt_case.FEV1FVC)))
print(mt_ctrl.aggregate_cols(hl.agg.stats(mt_ctrl.FEV1FVC)))

print('Pack-years of cases =', mt_case.aggregate_cols(hl.agg.stats(hl.float(mt_case.PackYear))))
print('Pack-years of controls =', mt_ctrl.aggregate_cols(hl.agg.stats(hl.float(mt_ctrl.PackYear))))

print('#Current Smoker Cases:', mt_case.aggregate_cols(hl.agg.counter(mt_case.CurrentSmoker)))
print('#Current Smoker Controls:', mt_ctrl.aggregate_cols(hl.agg.counter(mt_ctrl.CurrentSmoker)))






######## 5 POPULATION STRATIFICATION CORRECTION

######## 5.1 Common variants to preserve power
mt = hl.variant_qc(mt)
mt_common = mt.filter_rows(mt.variant_qc.AF[1] > 0.05)

#Pruning. PCA is sensitive to linkage disequilibrium
#r2 = 0 means that variants are unrelated, and when 1 they are in full linkage disequilibrium. r = 0.5 means removing those with correlation between 0.5 and 1.
pruned_t = hl.ld_prune(mt_common.GT, r2 = 0.5, bp_window_size = 500000)
mt_pruned = mt_common.filter_rows(hl.is_defined(pruned_t[mt_common.row_key]))

# Calculate eigenvalues, PC scores (k=10) and loadings
eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt_common.GT, k = 10, compute_loadings = True)
mt = mt.annotate_cols(scores = scores[mt.s].scores)

# Plot first 2 PCs
pca = hl.plot.scatter(mt_common.scores[0],
                      mt_common.scores[1],
                      label = mt_common.Race,
                      title = 'PCA Caucasian', xlabel = 'PC1', ylabel = 'PC2')
show(pca)

######## 5.2 Optional: Project new samples on existing PCs
def pc_project(
        # reference: https://github.com/macarthur-lab/gnomad_hail/blob/master/utils/generic.py#L131
        mt: hl.MatrixTable,
        loadings_ht: hl.Table,
        loading_location: str = "loadings",
        af_location: str = "pca_af"
) -> hl.Table:
    n_variants = loadings_ht.count()
    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location]
    )
    mt = mt.filter_rows(hl.is_defined(mt.pca_loadings) & hl.is_defined(mt.pca_af) &
                        (mt.pca_af > 0) & (mt.pca_af < 1))
    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(n_variants * 2 * mt.pca_af * (1 - mt.pca_af))
    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))
    return mt.cols().select('scores')
eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt_common.GT, k = 10, compute_loadings = True)
mt_common = mt_common.annotate_cols(scores = scores[mt_common.s].scores)
mt_common = mt_common.annotate_rows(pca_af = hl.agg.mean(mt_common.GT.n_alt_alleles()) / 2)
loadings = loadings.annotate(pca_af = mt_common[loadings.key, :].pca_af)
related_scores = pc_project(related_mt, pca_loadings)

######## 5.3 Optional: Plot pca manually
x = pca_scores.scores[0]
y = pca_scores.scores[1]
label = mt.cols()[pca_scores.s].Race
collect_all = nullable(bool)

if isinstance(x, Expression) and isinstance(y, Expression):
        agg_f = x._aggregation_method()
        if isinstance(label, Expression):
            if collect_all:
                res = hail.tuple([x, y, label]).collect()
                label = [point[2] for point in res]
            else:
                res = agg_f(aggregators.downsample(x, y, label=label, n_divisions=n_divisions))
                label = [point[2][0] for point in res]
            x = [point[0] for point in res]
            y = [point[1] for point in res]
        else:
            if collect_all:
                res = hail.tuple([x, y]).collect()
            else:
                res = agg_f(aggregators.downsample(x, y, n_divisions=n_divisions))
            x = [point[0] for point in res]
            y = [point[1] for point in res]

arg = list(set(label))
color=[]
for i in label :
    if i == arg[1]:
        color.append('red')
    elif i == arg[2] :
        color.append('black')
    else :
        color.append('blue')

p = figure(title='PCA')
fields = dict(x=x, y= y, label=label, color=color)
source = ColumnDataSource(fields)
p.circle( x='x', y='y', legend='label', color = 'color', source=source)
show(p)






######## 6 GWAS
# Keep at leats allele frequency of 1% (common + rare)
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)

######## 6.1 Firth logistic regression ~ Affection status
# Covariates:
### 1.0 is input variable number of alternate alleles with input variable the genotype dosage derived from the PL field,
### number of pack years smoking,
### population structure by PCA scores,
gwas = hl.logistic_regression_rows(
            test = "firth",   #controls false positives
            y = hl.float(mt.AffectionBool),
            x = mt.GT.n_alt_alleles(),
            covariates = [1, hl.float(mt.PackYear),
                          mt.scores[0], mt.scores[1],
                         mt.scores[2], mt.scores[3],
                         mt.scores[4], mt.scores[5],
                         mt.scores[6], mt.scores[7],
                         mt.scores[8], mt.scores[9]])

######## 6.2 Firth linear regression ~ FEV1% predicted (lung function variable)
gwas = hl.linear_regression_rows(
            y = hl.float(mt.FEV1),
            x = mt.GT.n_alt_alleles(),
            covariates = [1, hl.float(mt.PackYear),
                          mt.scores[0], mt.scores[1],
                         mt.scores[2], mt.scores[3],
                         mt.scores[4], mt.scores[5],
                         mt.scores[6], mt.scores[7],
                         mt.scores[8], mt.scores[9]])

######## 6.3 Q-Q plot
qqplot = hl.plot.qq(gwas.p_value)
show(qqplot)

######## 6.4 Manhattan-like plots
#GWAS significanse level = 5.0 10e-8, suggestive: 5.0 10e-8 < P < 5.0 * 10e-6.

# Calculate Bonferroni based cut off lines
Bonferroni_line = -np.log10(0.05 / mt.count_rows())
Suggestive_line = -np.log10(1 / mt.count_rows())

#Plot manhattan with conventional GWAS significance line
manh = hl.plot.manhattan(gwas.p_value,
                         title = "Manhattan-like Plot",
                         size = 4,
                         significance_line = 5e-08)
# Add lines to plot
line1 = Span(location = Bonferroni_line, dimension = "width", line_color = "red", line_width = 1)
line2 = Span(location = Suggestive_line, dimension = "width", line_color = "orange", line_width = 1)
manh.renderers.extend([line1, line2])

show(manh)






######## 6.5 Optional: FDR correction on P-values

#Create pandas df of locus and P values from GWAS
pval = gwas.select(gwas.locus, gwas.p_value)
pval_pd = pval.to_pandas(flatten = True)
pval = pval_pd.iloc[:,2]

#FDR correction function
def FDR(pvalues):
    pvalues = array(pvalues)
    n = len(pvalues)
    new_pvalues = empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort() # is already sorted
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

adjPval = FDR(pval)
df = pd.concat([pval_pd, pd.DataFrame(adjPval)], axis = 1)
df.rename(index = str, columns = {0:"adjPval"}, inplace = True)

#Export to tsv format
df.to_csv("/home/hail/adjPval.tsv", sep = "\t", header=True)

#Convert tsv to Hail table format
df2 = hl.import_table("file:///home/hail/adjPval.tsv",
                      types = {"locus.contig": hl.tstr,
                               "locus.position": hl.tint32,
                               "p_value": hl.tfloat64,
                               "adjPval": hl.tfloat64})
df2 = df2.annotate(locus=hl.locus(df2["locus.contig"], df2["locus.position"]))
df2 = df2.select(df2.locus, df2.adjPval)

#Add adjusted P-values to MT file
mt = mt.annotate_rows(**df2[mt.locus])







######## 7. TOP HITS

######## 7.1 Look into top hits
# Order based on p-value
gwas_ordered = gwas.order_by(gwas.p_value)
gwas_ordered.p_value.show(25)
# Check allele frequencies (pathogenic variants are expected to be low-frequency)
mt_hits.variant_qc.AF[1].show()
# Check which ones are SNPs
mt_snp = mt_top.filter_rows(hl.is_snp(mt_top.alleles[0], mt_top.alleles[1]))

######## 7.2 Filter on suggestive significance level
signPval = gwas.filter(gwas.p_value < 5e-06)
signPval = signPval.key_by("locus")
mt_hits = mt.filter_rows(hl.is_defined(signPval[mt.locus]))
# Show locus, alleles, and allele frequency
mt_hits.variant_qc.AF[1].show(mt_hits.count())

######## 7.3 Check linkage disequilibrium by pruning
#The pruning algorithm finds correlation between genetic variants and keeps the first in chronological sequence order of the linkage disequilibrium block.
pruned_t = hl.ld_prune(mt_hits.GT, r2 = 0.1, bp_window_size = 1000000000)
mt_pruned = mt_hits.filter_rows(hl.is_defined(pruned_t[mt_hits.row_key]))
mt_pruned.variant_qc.AF[1].show()

######## 7.4 Samples per top hits
# Merge the data from the matrices globals, columns and rows into a single matrix
t = mt_hits.entries()
# Remove keys, because they are kept even after select()
t = t.key_by()
# Select the columns to keep
t_GT = t.select(t.GT, t.dbGaP_Subject_ID, t.locus)
# Check
t_GT.describe()

######## 7.4 Export as csv file
# First convert to pandas format
t_pd = t_GT.to_pandas(flatten = True)
# Write a csv in the local directory
t_pd.to_csv('<name-csv>', encoding = 'utf-8', index = False)
