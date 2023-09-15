.. _quick_tutorials:

===============
Quick tutorials
===============

We provide four quick tutorials; each shows one function of KGGSEE. In each tutorial, we provide the command line and a brief explanation of flags and output files. Please refer to :ref:`Detailed Document <detailed_document>` and :ref:`Options <options>` for details. The first tutorial (:ref:`Gene-based association tests <t1>`) should be done first, then you can run any of the following tutorials.

Make sure the KGGSEE Java archive ``kggsee.jar``, the running resource data folder ``resources/``, and the tutorial data folder ``tutorials/`` are under the same directory. Then, we suppose you are under the directory ``tutorials/``.


.. _t1:

Gene-based association tests
============================

GATES and ECS are two statistical methods combining the p-values of a group of SNPs into one p-value. This analysis inputs p-values of SNPs and outputs p-values of genes. The command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
    --sum-file scz_gwas_eur_chr1.tsv.gz \
    --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
    --keep-ref \
    --gene-assoc \
    --out t1


**Explanation of the flags and input files:**

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
      - Keep the parsed VCF file (KGGSEE object format) in a folder named ``VCFRefhg19`` under the output folder. KGGSEE will read these files in the following tutorials, which will be faster than parsing VCF files.
    * - ``--gene-assoc``
      - Triggers gene-based association tests.
    * - ``--out``
      - Specifies the prefix of output files.


**Explanation of the output files:**

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
      - Number of variants within the gene
    * - ECSP
      - p-value of ECS
    * - GATESP
      - p-value of GATES
    * - Chrom
      - Chromosome of the gene
    * - Pos
      - Coordinate of the variant with the lowest p-value within the gene
    * - GWAS_Var_P
      - p-value of the variant

The columns of ``t1.gene.var.pvalue.txt.gz`` are the same as ``t1.gene.pvalue.txt``. The difference is that, for each gene, in ``t1.gene.pvalue.txt``, only the variant with the lowest p-value is output, while in ``t1.gene.var.pvalue.txt.gz``, all variants are output.

The Q-Q plots for p-values of inputted GWAS file (inside or outside of gene) and gene-based association tests by GATES or ECS are saved in ``t1.qq.png``.


.. _t2:

DESE
====
    
DESE performs phenotype-tissue association tests and conditional gene-based association tests at the same time. This analysis inputs p-values of a GWAS and expression profile of multiple tissues; outputs p-values of phenotype-tissue associations and conditional p-values of genes. The command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
    --sum-file scz_gwas_eur_chr1.tsv.gz \
    --saved-ref VCFRefhg19 \
    --expression-file GTEx_v8_TMM.gene.meanSE.txt.gz \
    --gene-condi \
    --out t2


**Explanation of the flags and input files:**

.. list-table::
    :widths: 1 3
    :header-rows: 1
    :class: tight-table

    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, columns of SNP coordinates and p-values are needed.
    * - ``--saved-ref``
      - Specifies the folder of genotypes of reference population in KGGSEE object format, which is saved by the ``--keep-ref`` flag in the first tutorial.
    * - ``--expression-file``
      - Specifies a gene expression file that contains means and standard errors of gene expressions in multiple tissues/cell types. Here ``GTEx_v8_TMM.gene.meanSE.txt`` is for gene-level DESE. Try ``GTEx_v8_TMM.transcript.meanSE.txt`` for transcript-level DESE.
    * - ``--gene-condi``
      - Triggers the DESE analysis.
    * - ``--out``
      - Specifies the prefix of output files.


**Explanation of the output files:**

The three files of ``t2.gene.pvalue.txt``, ``t2.gene.var.pvalue.txt.gz``, and ``t2.qq.png`` are the same as their counterparts with the same suffixes of the first tutorial. In addition, the results of conditional gene-based association tests are in ``t2.finemapping.gene.ecs.txt`` which contains nine columns:

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
      - LD group number. Conditional ECS tests were performed for genes within a same LD group.
    * - ECSP
      - p-value of ECS
    * - CondiECSP
      - p-value of conditional gene-based association tests by conditional ECS
    * - GeneScore
      - The gene's selective expression score in all tissues. A gene with a high socre will be given higher priority to enter the conditioning procedure.
       

Results of driver-tissue prioritizations are in ``t2.celltype.txt``. This is basically a Wilcoxon rank-sum test which tests whether the selective expression median of the phenotype-associated genes is significantly higher than that of other genes in an interrogated tissue. The file contains three columns:

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - TissueName
      - Name of the tissue being tested
    * - p
      - These p-values are for tissue prioritization but NOT for hypothesis test.
    * - BHFDRq
      - The Benjamini-Hochberg adjusted p-values


.. _t3:

EMIC
====

EMIC inferences gene expressions' causal effect on a complex phenotype with dependent expression quantitative loci by a robust median-based Mendelian randomization. SNPs with effects on both the phenotype and a gene are considered instrumental variables (IVs) of the gene, which can be used to infer the gene's expression effect on the phenotype. This analysis inputs effect sizes of SNPs on the phenotype and genes' expressions; outputs effect sizes and p-values of genes' expression effects on the phenotype. The command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
    --sum-file scz_gwas_eur_chr1.tsv.gz \
    --saved-ref VCFRefhg19 \
    --eqtl-file GTEx_v8_gene_BrainBA9.eqtl.txt.gz \
    --beta-col OR \
    --beta-type 2 \
    --emic \
    --out t3


**Explanation of the flags and input files:**

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, in addition to the columns of SNP coordinates and p-values, two columns of SNP alleles (named A1 and A2 by default), a column of A1 allele frequency (named FRQ_U by default), and two columns of SNP effect sizes (no default header) and their standard errors (named SE by default) are also needed.
    * - ``--saved-ref``
      - Specifies the folder of genotypes of reference population in KGGSEE object format, which is saved by the ``--keep-ref`` flag in the first tutorial.
    * - ``--eqtl-file``
      - Specifies a fasta-styled file of SNPs' effects on gene expressions. Here ``GTEx_v8_gene_BrainBA9.eqtl.txt.gz`` for gene-level EMIC. Try ``GTEx_v8_transcript_BrainBA9.eqtl.txt.gz`` for transcript-level EMIC.
    * - ``--beta-col``
      - Specifies the column name of effect sizes in the GWAS file.
    * - ``--beta-type``
      - Specifies the type of the effect size; here ``2`` means that it is the odds ratio for a qualitative phenotype.
    * - ``--emic``
      - Triggers the EMIC analysis.
    * - ``--out``
      - Specifies the prefix of output files.


**Explanation of the output files:**

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
      - Chromosome of the gene
    * - Pos
      - The coordinate of the IV with the lowest GWAS p-value
    * - GWAS_Var_P
      - GWAS p-value of the IV
    * - GWAS_Var_Beta
      - The phenotype association effect size of the IV
    * - GWAS_Var_SE
      - Standard error of the effect size


The columns of ``t3.emic.gene.var.tsv.gz`` are the same as ``t3.emic.gene.txt``. The difference is that, for each gene, in ``t3.emic.gene.txt``, only the eQTL with the lowest GWAS p-value is output, while in ``turorial_3.emic.gene.var.tsv.gz``, all eQTLs are output. In this tutorial, the file ``t3.emic.gene.PleiotropyFinemapping.txt`` is empty, we ignore it here.

File ``t3.qq.png`` saves the Q-Q plot for GWAS p-values of IVs. File ``t3.emic.qq.png`` saves the Q-Q plot for EMIC p-values. 

File ``t3.scatterplots.emic.pdf`` saves the scatter plots of genetic association with gene expression. Each gene with an EMIC p-value lower than 2.5E-3 (default threshold) is saved on a separate page of the PDF. A filled rectangle on the plots denotes an IV. The red rectangle denotes the most significant GWAS variant among all the IVs of a gene. The slope of the line represents the estimated causal effect. The color of an IV denotes the degree of the LD between the IV and the most significant GWAS variant. The error bars in the rectangles denote the standard errors of the coefficient estimates.


.. _t4:

Gene-based heritability estimation
==================================
    
Heritability is a measure of how well differences in people's genes account for differences in their phenotypes. This tutorial estimates the heritability of each gene with GWAS summary statistics. The command is:

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
    --sum-file scz_gwas_eur_chr1.tsv.gz \
    --saved-ref VCFRefhg19 \
    --case-col Nca \
    --control-col Nco \
    --gene-herit \
    --out t4


**Explanation of the flags and input files:**

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Flag
      - Description
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics. In this analysis, in addition to the columns of SNP coordinates and p-values, two columns of case and control sample sizes are also needed.
    * - ``--saved-ref``
      - Specifies the folder of genotypes of reference population in KGGSEE object format, which is saved by the ``--keep-ref`` flag in the first tutorial.
    * - ``--case-col``
      - Specifies the column name of the case sample size.
    * - ``--control-col``
      - Specifies the column name of the control sample size.
    * - ``--gene-herit``
      - Triggers gene-based association tests and estimation of gene heritability.
    * - ``--out``
      - Specifies the prefix of output files.


**Explanation of the output files:**

The output files are generally the same as the first tutorial, except that, in ``t4.gene.pvalue.txt``, ``t4.gene.var.pvalue.txt.gz``, there are two more columns named SNPHerit and SNPHeritSE, which are the estimate and its standard error of the gene heritability.

