.. _quick_tutorials:

===============
Quick tutorials
===============

We provide four quick tutorials; each shows one function of KGGSEE. In each tutorial, we provide the command line and a brief explanation of flags and output files. Please refer to :ref:`Detailed Document <detailed_document>` and :ref:`Options <options>` for details. The first tutorial (:ref:`GATES and ECS (gene-based association tests) <t1>`) should be done first, then you can run any of the following tutorials.

Make sure the KGGSEE Java archive ``kggsee.jar``, the running resource data folder ``resources/``, and the tutorial data folder ``tutorials/`` are under the same directory. Then, we suppose you are under the directory ``tutorials/``.


.. _t1:

GATES and ECS (gene-based association tests)
============================================

GATES and ECS are two statistical methods combining the p-values of a group of SNPs into one p-value. This analysis inputs p-values of SNPs and outputs p-values of genes. The tutorial command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --sum-file ./scz_gwas_eur_chr1.tsv.gz \
      --vcf-ref ./1kg_hg19_eur_chr1.vcf.gz \
      --keep-ref ./VCFRefhg19/ \
      --gene-assoc \
      --out t1


Options and input files
-----------------------
.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, columns of SNP coordinates and p-values (CHR, BP, and P by default) are needed.
    * - ``--vcf-ref``
      - Specifies a VCF file of genotypes sampled from a reference population. These genotypes are used to estimate LD correlation coefficients among SNPs.
    * - ``--keep-ref``
      - Keep the parsed VCF file in KGGSEE object format in the specified directory. KGGSEE will read these files in the following tutorials, which will be faster than parsing VCF files.
    * - ``--gene-assoc``
      - Triggers gene-based association tests.
    * - ``--out``
      - Specifies the prefix of output files.


Output files
------------
The numeric results of gene-based association tests are saved in ``t1.gene.pvalue.txt``. There are seven columns in the file:

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - Gene
      - Gene symbol
    * - #Var
      - Number of variants within a gene
    * - ECSP
      - p-value of ECS
    * - GATESP
      - p-value of GATES
    * - Chrom
      - Chromosome of the gene
    * - Pos
      - Coordinate of the variant with the lowest p-value within the gene
    * - GWAS_Var_P
      - p-value of the variant with the lowest p-value within the gene

The columns of ``t1.gene.var.pvalue.txt.gz`` are the same as ``t1.gene.pvalue.txt``. The difference is that, for each gene, in ``t1.gene.pvalue.txt``, only the variant with the lowest p-value is output, while in ``t1.gene.var.pvalue.txt.gz``, all variants are output.

The Q-Q plots for p-values of input GWAS file (inside or outside of each gene) and gene-based association tests by GATES or ECS are saved in ``t1.qq.pdf``.


.. _t2:

DESE (driver-tissue inference)
==============================
    
DESE performs phenotype-tissue association tests and conditional gene-based association tests at the same time. This analysis inputs p-values of a GWAS and expression profile of multiple tissues and outputs p-values of phenotype-tissue associations and conditional p-values of genes. The tutorial command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --sum-file ./scz_gwas_eur_chr1.tsv.gz \
      --saved-ref ./VCFRefhg19/ \
      --expression-file ../resources/GTEx_v8_TMM_all.gene.meanSE.txt.gz \
      --gene-assoc-condi \
      --out t2


Options and input files
-----------------------
.. list-table::
    :widths: 1 3
    :header-rows: 1
    :class: tight-table

    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, columns of SNP coordinates and p-values are needed.
    * - ``--saved-ref``
      - Specifies the directory of the genotypes of the reference population in KGGSEE object format, which is saved by the ``--keep-ref`` flag in the first tutorial.
    * - ``--expression-file``
      - Specifies a gene expression file that contains means and standard errors of gene expressions for tissues/cell types. Here ``GTEx_v8_TMM_all.gene.meanSE.txt`` is for gene-level DESE. Try ``GTEx_v8_TMM_all.transcript.meanSE.txt`` for transcript-level DESE.
    * - ``--gene-assoc-condi``
      - Triggers the DESE analysis.
    * - ``--out``
      - Specifies the prefix of the output files.


Output files
------------
The three files of ``t2.gene.pvalue.txt``, ``t2.gene.var.pvalue.txt.gz``, and ``t2.qq.pdf`` are the same as their counterparts with the same suffixes of the first tutorial. In addition, the results of conditional gene-based association tests are in ``t2.gene.assoc.condi.txt`` which contains nine columns:

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - Gene
      - Gene symbol
    * - Chrom
      - Chromosome of the gene
    * - StartPos
      - Start coordinate of the gene
    * - EndPos
      - End coordinate of the gene
    * - #Var
      - Number of variants within the gene
    * - Group
      - LD group number. Conditional ECS tests were performed for genes within the same LD group.
    * - ECSP
      - p-value of ECS
    * - CondiECSP
      - p-value of the conditional gene-based association tests by conditional ECS
    * - GeneScore
      - Gene's selective-expression score. A gene with a high score will be given higher priority to enter the conditioning procedure.
       

Results of driver-tissue prioritizations are in ``t2.celltype.txt``. This is a Wilcoxon rank-sum test which tests whether the selective expression median of the phenotype-associated genes is significantly higher than that of the other genes in the interrogated tissue. The file contains three columns:

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - TissueName
      - Name of the tissue being tested
    * - Unadjusted(p)
      - Unadjusted p-values for the tissue-phenotype associations
    * - Adjusted(p)
      - Adjusted p-values calculated by adjusting both selection bias and multiple testing
    * - Median(IQR)SigVsAll
      - Median (interquartile range) expression of the conditionally significant genes and all the background genes


.. _t3:

EMIC (gene-expression causal-effect inference)
==============================================

EMIC inferences gene expressions' causal effect on a complex phenotype with dependent expression quantitative loci by a robust median-based Mendelian randomization. SNPs with effects on both the phenotype and a gene are considered as instrumental variables (IVs) of the gene, which can be used to infer the gene's expression effect on the phenotype. This analysis uses effect sizes of SNPs on the phenotype and genes' expressions and outputs effect sizes and p-values of the expression effects on the phenotype. The tutorial command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --sum-file ./scz_gwas_eur_chr1.tsv.gz \
      --saved-ref ./VCFRefhg19/ \
      --eqtl-file ./GTEx_v8_gene_BrainBA9.eqtl.txt.gz \
      --emic-plot-p 0.01 \
      --beta-col OR \
      --beta-type 2 \
      --emic \
      --out t3


Options and input files
-----------------------
.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, in addition to the columns of SNP coordinates and p-values, two columns of SNP alleles (named A1 and A2 by default), a column of the effect allele (A1) frequency (named FRQ_U by default), and two columns of SNP effect sizes and their standard errors (named SE by default) are also needed.
    * - ``--saved-ref``
      - Specifies the directory of genotypes of reference population in KGGSEE object format, which is saved by the ``--keep-ref`` flag in the first tutorial.
    * - ``--eqtl-file``
      - Specifies a fasta-styled file of SNPs' effects on gene expressions. Here ``GTEx_v8_gene_BrainBA9.eqtl.txt.gz`` is for gene-level EMIC. You can try ``GTEx_v8_transcript_BrainBA9.eqtl.txt.gz`` for a transcript-level EMIC.
    * - ``--emic-plot-p``
      - Specifies the p-value threshold for plotting a scatter plot.
    * - ``--beta-col``
      - Specifies the column name of effect sizes in the GWAS file.
    * - ``--beta-type``
      - Specifies the type of the effect sizes; here ``2`` means that it is the odds ratio for a qualitative phenotype.
    * - ``--emic``
      - Triggers the EMIC analysis.
    * - ``--out``
      - Specifies the prefix of the output files.


Output files
------------
The numeric results of EMIC are saved in ``t3.emic.gene.txt``. There are nine columns in the file:

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - Gene
      - The gene symbol
    * - #Var
      - Number of IVs within the gene
    * - minP_EMIC
      - p-value of EMIC. When a transcript-level EMIC is performed, this is the minimum p-value among all transcripts of the gene.
    * - Details_EMIC
      - Each detailed result has four components in brackets: the number of IVs, the causal effect estimate and its standard error, and the p-value. When a transcript-level EMIC is performed, results for each transcript are listed.
    * - Chrom
      - Chromosome of a gene
    * - Pos
      - The coordinate of the IV with the lowest GWAS p-value
    * - GWAS_Var_P
      - GWAS p-value of an IV
    * - GWAS_Var_Beta
      - The phenotype association effect size of an IV
    * - GWAS_Var_SE
      - Standard error of an effect size


The columns of ``t3.emic.gene.var.tsv.gz`` are the same as ``t3.emic.gene.txt``. The difference is that, for each gene, in ``t3.emic.gene.txt``, only the eQTL with the lowest GWAS p-value is output, while in ``turorial_3.emic.gene.var.tsv.gz``, all eQTLs are output. In this tutorial, the file ``t3.emic.gene.PleiotropyFinemapping.txt`` is empty, we ignore it here.

File ``t3.qq.pdf`` saves the Q-Q plot for the GWAS p-values of IVs. File ``t3.emic.qq.pdf`` saves the Q-Q plot for the EMIC p-values. 

File ``t3.scatterplots.emic.pdf`` saves the scatter plots of the genetic association with gene expression. Each gene with an EMIC p-value lower than 2.5E-3 (default threshold) is saved on a separate page of the PDF. A filled rectangle on the plots denotes an IV. The red rectangle denotes the most significant GWAS variant among all the IVs of a gene. The slope of the line represents the estimated causal effect. The color of an IV denotes the degree of the LD between the IVs and the most significant GWAS variant. The error bar in a rectangle denotes the standard error of the coefficient estimate. File ``t3.scatterplots.emic.txt`` saves the numeric results of the scatter plots in ``t3.scatterplots.emic.pdf``.


.. _t4:

EHE (gene-based heritability estimation)
========================================
    
Heritability is a measure of how well differences in people's genes account for differences in their phenotypes. This tutorial estimates the heritability of each gene using GWAS summary statistics. The tutorial command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --sum-file ./scz_gwas_eur_chr1.tsv.gz \
      --saved-ref ./VCFRefhg19/ \
      --case-col Nca \
      --control-col Nco \
      --gene-herit \
      --out t4


Options and input files
-----------------------
.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, in addition to the columns of SNP coordinates and p-values, two columns of case and control sample sizes are also needed.
    * - ``--saved-ref``
      - Specifies the directory of the genotypes of the reference population in KGGSEE object format, which is saved by the ``--keep-ref`` flag in the first tutorial.
    * - ``--case-col``
      - Specifies the column name of the case sample size.
    * - ``--control-col``
      - Specifies the column name of the control sample size.
    * - ``--gene-herit``
      - Triggers gene-based association tests and estimations of gene heritability.
    * - ``--out``
      - Specifies the prefix of the output files.


Output files
------------
The output files are generally the same as the first tutorial, except that, in ``t4.gene.pvalue.txt``, ``t4.gene.var.pvalue.txt.gz``, there are two more columns named ``Herit`` and ``HeritSE``, which are the estimate and its standard error of a gene's heritability.

