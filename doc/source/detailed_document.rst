.. _detailed_document:

=================
Detailed Document
=================

We first describe the general aspects of all analyses and then describe details for each analysis. KGGSEE performs analysis according to the following procedure:

1. Reads genotypes of an ancestrally matched reference population, e.g., a super population panel of the 1000 Genomes Project. The genotype VCF file is specified by ``--vcf-ref``, and if ``--keep-ref`` is used, KGGSEE saves the parsed VCF files in KGGSEE object format in the specified directory. In a following run, specifying the directory by ``--saved-ref``, KGGSEE reads genotypes from the object format files, which is faster than re-parsing the VCF files. KGGSEE calculates the minor allele frequency (MAF) of each SNP and filters out SNPs with an MAF lower than the threshold specified by ``--filter-maf-le`` (default: ``0.05``). KGGSEE also calculates the p-value of rejecting Hardy-Weinberg equilibrium for each SNP and filters out SNPs with a p-value lower than the threshold specified by ``--hwe-all`` (default: ``1E-5``). Only SNPs with genotypes of the reference population and who have passed the two filters will be considered in the following procedures.

2. Reads GWAS summary statistics from a whitespace delimited file specified by ``--sum-file``. Depending on the analysis performed, this file needs to have different columns, which will be described in the corresponding sections. For gene/transcript-expression causal-effect estimation, an eQTL summary statistic file has to be specified by ``--eqtl-file``. For the other gene-based analyses, a gene is a set of SNPs within a genomic region. By default, KGGSEE reads gene annotations (``resources/{hg19,hg38}/*{GEncode,refGene}.txt.gz``) to find out which SNPs should be considered as one SNP set. Instead, an eQTL summary statistic file can also be specified to let KGGSEE consider eQTLs of a gene/transcript as a SNP set. In addition to this, a user can also specify a genomic range file by ``--regions-bed`` to let KGGSEE consider SNPs in each range as a gene. We provide gene and transcript-based eQTL summary statistics for 49 GTEx v8 tissues available for downloading on our OneDrive (`hg19 <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EnhWhqLUNcpOrh6O3enFvCUBRvQ13v2970tcpOnNmmlKyg?e=1jkl06>`_ or `hg38 <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EtWxtqj5HTRHsEw4IiZ9xAMBu9S8Defi67pmL3_rNUjb9w?e=ufFapJ>`_).

3. Based on the flag specified, KGGSEE reads more needed files and performs the corresponding analysis.

    * :ref:`GATES and ECS (gene-based association tests) <detail_ECS>` trigered by ``--gene-assoc``;
    * :ref:`DESE (driver-tissue inference) <detail_DESE>` trigered by ``--gene-assoc-condi``;
    * :ref:`EMIC (gene-expression causal-effect inference) <detail_EMIC>` trigered by ``--emic``;
    * :ref:`EHE (gene-based heritability estimation) <detail_h2>` trigered by ``--gene-herit``.


.. _eqtl_file:

The KGGSEE format of eQTL summary statistics is fasta-styled. E.g.,

.. code::

    #symbol   id	chr	pos	ref	alt	altfreq	beta	se	p	neff	r2
    >WASH7P	ENSG00000227232	1
    52238	T	G	0.94	-1.77	0.28	5.1E-9	65	0.38
    74681	G	T	0.95	-1.45	0.33	1.1E-5	63	0.23
    92638	A	T	0.24	0.54	0.20	7.9E-3	53	0.12
    >MIR130	ENSG00000284557	1
    52238	T	G	0.94	-1.77	0.28	5.1E-9	65	0.38
    74681	G	T	0.95	-1.45	0.33	1.1E-5	63	0.23

The first row starting with ``#`` is the header line. Then, eQTLs of each gene/transcript are chunked. For each gene/transcript, the first row has three columns: (1) the gene symbol prefixed by ``>``, (2) Ensembl gene/transcript ID, and (3) chromosome. The second and following rows have nine columns: (1) the eQTL coordinate, (2) the reference allele, (3) the alternative allele, (4) the frequency of the alternative allele, (5) the effect size, (6) the standard error of the effect size, (7) the p-value of nonzero effect size, (8) the effective sample size, and (9) coefficient of determination.

A BED file specified by ``--regions-bed`` defines customized gene coordinates instead of the annotation from RefSeqGene or GENCODE. The first three columns of the BED file define gene coordinates and are mandatory; the fourth column defines gene names and is optional. When the fourth column is absent, a gene name of the format like chr1:100-200 will be allocated.

.. code::

    1       59091     80008     OR4F5
    1       850260    889955    SAMD11
    1       869584    904689    NOC2L
    1       885967    911095    KLHL17
    2       28814     56870     FAM110C
    2       207730    276398    SH3YL1
    2       254140    288283    ACP1


.. _detail_ECS:

GATES and ECS (gene-based association tests)
============================================

KGGSEE performs the gene-based association analysis by GATES (a rapid and powerful **G**\ ene-based **A**\ ssociation **T**\ est using **E**\ xtended **S**\ imes procedure) and ECS (an **E**\ ffective **C**\ hi-square **S** \tatistics). The ``--gene-assoc`` flag triggers both.

GATES (`Li et al. 2011 <https://doi.org/10.1016/j.ajhg.2011.01.019>`_) is basically an extension of the Simes procedure to dependent tests, as the individual GWAS tests are dependent due to LD. GATES calculates an effective number of independent p-values which is then used by a Simes procedure.

ECS (`Li et al. 2019 <https://doi.org/10.1093/bioinformatics/bty682>`_) first converts the p-values of a gene to chi-square statistics(one degree of freedom). Then, merges all chi-square statistics of a gene after correcting the redundancy of the statistics due to LD. The merged statistic is called an ECS which is used to calculate the p-value of the gene. 


Synopsis
--------

.. code:: shell

    java -Xms16g -Xmx16g -jar kggsee.jar
      --gene-assoc
      --out <prefix>
      --vcf-ref <file>
      --sum-file <file>
      --chrom-col <header>  # default: CHR
      --pos-col <header>  # default: BP
      --p-col <header>  # default: P 
      --neargene <basepair>  # default: 5000
      --eqtl-file <file>
      --filter-eqtl-p <pval>  # default: 0.01


The flag ``--gene-assoc`` triggers the gene-based association tests. ``--sum-file`` specifies a white space-delimited GWAS summary statistic file which must have three columns of the chromosome of SNP, coordinate of SNP, and p-value of SNP; headers of the three columns can be specified by ``--chrom-col``, ``--pos-col`` and ``--p-col`` separately. SNPs belonging to a gene can be defined either by SNPs close to the gene or by eQTLs of the gene. If ``--neargene`` is specified, KGGSEE reads gene annotations and considers SNPs inside a gene and its adjacent regions at a fixed number of basepairs on both sides to be a test unit. If ``--eqtl-file`` is specified, KGGSEE reads the eQTL summary statistic file and considers eQTLs of a gene or a transcript to be a test unit, and ``--neargene`` is overridden. When ``--eqtl-file`` is specified, ``--filter-eqtl-p`` can be used to specify a threshold of eQTL p-values. Only eQTLs with a p-value lower than the threshold will be considered. :ref:`A description of the eQTL file format <eqtl_file>` is near the beginning of the page.


Examples
--------


Test based on physical distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, SNPs inside a gene and its 10 kb adjacent regions will be grouped for association tests.

.. code:: shell

      java -Xmx4g -jar ../kggsee.jar \
        --gene-assoc \
        --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
        --sum-file scz_gwas_eur_chr1.tsv.gz \
        --neargene 10000 \
        --out t1.1


Test guided by eQTLs of a gene
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, eQTLs of a gene will be grouped for association tests.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --gene-assoc \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --eqtl-file GTEx_v8_gene_BrainBA9.eqtl.txt.gz \
      --out t1.2


Test guided by eQTLs of a transcript
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, eQTLs of a transcript will be grouped for association tests.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --gene-assoc \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --eqtl-file GTEx_v8_transcript_BrainBA9.eqtl.txt.gz \
      --out t1.3



Outputs
-------
The file with a suffix of ``.gene.pvalue.txt`` saves the results of gene-based association tests. The columns of the file are as follows:


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
      - The coordinate of the variant with the lowest p-value within the gene
    * - GWAS_Var_P
      - p-value of the variant


The columns of the file with the suffix of ``.gene.var.pvalue.txt.gz`` are the same as ``*.gene.pvalue.txt``. The difference is that, for each gene, in ``*.gene.pvalue.txt``, only the variant with the lowest p-value is output, while in ``*.gene.var.pvalue.txt.gz``, all variants are output. The file with the suffix of ``.qq.png`` is the Q-Q plots for p-values of GWAS summary statistics and gene-based association tests by GATES and ECS.


.. admonition:: Citation of GATES

    Miaoxin Li, Hong-Sheng Gui, Johnny Sheung Him Kwan and Pak Chung Sham. GATES: a rapid and powerful gene-based association test using extended Simes procedure. The American Journal of Human Genetics (2011). 88(3):283-293. https://doi.org/10.1016/j.ajhg.2011.01.019


.. admonition:: Citation of ECS

    Miaoxin Li, Lin Jiang, Timothy Shin Heng Mak, Johnny Sheung Him Kwan, Chao Xue, Peikai Chen, Henry Chi-Ming Leung, Liqian Cui, Tao Li and Pak Chung Sham. A powerful conditional gene-based association approach implicated functionally important genes for schizophrenia. Bioinformatics (2019). 35(4):628-635. https://doi.org/10.1093/bioinformatics/bty682


.. _detail_DESE:

DESE (driver-tissue inference)
==============================

DESE (**D**\ river-tissue **E**\ stimation by **S**\ elective **E**\ xpression; `Jiang et al. 2019 <https://doi.org/10.1186/s13059-019-1801-5>`_) estimates driver tissues by tissue-selective expression of phenotype-associated genes in GWAS. The assumption is that the tissue-selective expression of causal or susceptibility genes indicates the tissues where complex phenotypes happen primarily, which are called driver tissues. Therefore, a driver tissue is very likely to be enriched with selective expression of susceptibility genes of a phenotype. 

DESE initially performed the association analysis by mapping SNPs to genes according to their physical distance. We further demonstrated that grouping eQTLs of a gene or a transcript to perform the association analysis could be more powerful. We named the **e**\ QTL-guided **DESE** eDESE. KGGSEE implements DESE and eDESE with an improved effective chi-squared statistic to control type I error rates and remove redundant associations (`Li et al. 2022 <https://doi.org/10.7554/eLife.70779>`_).

.. note::
    We have developed an online service called PCGA (https://pmglab.top/pcga; `Xue et al. 2022 <https://doi.org/10.1093/nar/gkac425>`_), which implements DESE and integrates vast amounts of scRNA-seq datasets. PCGA collects and processes expression profile data from 54 tissue types and 6,598 cell types, enabling more convenient hierarchical estimation of the associated tissues and cell types of complex diseases. Additionally, PCGA has analyzed 1,871 public GWASs related to 1,588 unique phenotypes, which can be browsed and searched on the website.


Synopsis
--------

.. code:: shell

    java -Xms16g -Xmx16g -jar kggsee.jar
      --gene-assoc-condi
      --out <prefix>
      --vcf-ref <file>
      --sum-file <file>
      --chrom-col <header>  # default: CHR
      --pos-col <header>  # default: BP
      --p-col <header>  # default: P 
      --neargene <both-sides-bp|upstream-bp,downstream-bp>  # default: 5000
      --eqtl-file <file>
      --filter-eqtl-p <pval>  # default: 0.01
      --multiple-testing <bonf|benfdr|fixed>  # default: bonf
      --p-value-cutoff <pval>  # default: 0.05
      --top-gene <number>
      --expression-file <file>
      --geneset-db <cura|cgp|cano|cmop|onto|onco|immu>
      --geneset-file <file>
      --dese-permu-num <number>


The flag ``--gene-assoc-condi`` triggers DESE. First, KGGSEE performs gene-based association tests, which is the same as the analyses triggered by ``--gene-assoc``. ``--sum-file`` specifies a white space delimited GWAS summary statistic file which must have three columns of the chromosome of SNP, coordinate of SNP, and p-value of SNP; headers of the three columns can be specified by ``--chrom-col``, ``--pos-col`` and ``--p-col`` separately. SNPs belonging to a gene can be defined either by SNPs close to the gene or by eQTLs of the gene. If ``--neargene`` is specified by one number, KGGSEE reads gene annotations and considers SNPs inside a gene and its adjacent regions at a fixed number of basepairs on both sides to be a test unit. ``--neargene`` can also have two values to set an asymmetric boundary extension, e.g., 5 kb upstream and 15 kb downstream of a gene can be set by ``--neargene 5000,15000``. If ``--eqtl-file`` is specified, eDESE is evoked; KGGSEE reads the eQTL summary statistic file and considers eQTLs of a gene or a transcript to be a test unit, and ``--neargene`` is overridden. When ``--eqtl-file`` is specified, ``--filter-eqtl-p`` can be used to specify a threshold of eQTL p-values. Only eQTLs with a p-value lower than the threshold will be considered. :ref:`A description of the eQTL file format <eqtl_file>` is near the beginning of the page.

Second, after the gene-based association tests, significant genes by ECS are retained for fine-mapping. ``--multiple-testing`` specifies the method for multiple testing correction: ``bonf`` denotes Bonferroni correction; ``benfdr`` denotes Benjamini–Hochberg FDR; ``fixed`` denotes no correction. ``--p-value-cutoff`` specifies the threshold of the adjusted p-value. ``--top-gene`` specifies the maximum number of genes retained for fine-mapping. So, only genes (no more than the specified maximum number) with adjusted p-values lower than the specified threshold are retained for fine-mapping. Then, KGGSEE reads the expression file specified by ``--expression-file`` and performs iterative estimation of driver tissues. When ``--dese-permu-num`` is omitted, only unadjusted p-values are output. The unadjusted p-values are inflated due to selection bias in the iterations, which is only valid for tissue prioritization. For phenotype-tissue association tests, add ``--dese-permu-num 100`` for an adjustment by 100 permutations for selection bias and multiple testing.

Finally, if ``--geneset-db`` is specified, KGGSEE tests if the conditional significant genes are enriched in gene sets of `MSigDB <http://www.gsea-msigdb.org/gsea/msigdb/index.jsp>`_. The abbreviations of gene sets are as follows:

    | ``cura``: C2. curated gene sets;
    | ``cgp`` : C2. chemical and genetic perturbations;
    | ``cano``: C2. canonical pathways;
    | ``cmop``: C4. computational gene sets;
    | ``onto``: C5. ontology gene sets;
    | ``onco``: C6. oncogenic signature gene sets;
    | ``immu``: C7. immunologic signature gene sets.

Customized gene sets for enrichment tests can be specified by ``--geneset-file``. Please refer to ``resources/*.symbols.gmt.gz`` under the KGGSEE directory for file formats.


Expression files should be tab-delimitated. The first column is gene/transcript IDs. The IDs should be Ensembl gene IDs, Ensembl transcript IDs, or HGNC symbols. The version of Ensembl IDs will be trimmed by KGGSEE. For transcript-level expression profile,  a transcript label should be an Ensembl transcript ID and an ID of another type joint by ``:``.  Headers of the same tissue must have the same prefix. Headers of mean values must end with ``.mean``. Headers of SEs must end with ``.SE``. All SE values must be positive. The following columns are means and SEs of expression levels of genes or transcripts in multiple tissues. A gene-level expression file looks like this:

.. _expression_file:

.. code::

    Name               Tissue1.mean   Tissue1.SE     Tissue2.mean   Tissue2.SE     ...
    ENSG00000223972    0.0038016      0.00036668     0.0045709      0.00046303     ...
    ENSG00000227232    1.9911         0.030021       1.8841         0.040247       ...
    ENSG00000278267    0.00049215     0.00010645     0.00036466     9.2944E-05     ...
    ENSG00000243485    0.0047772      0.00038018     0.0067897      0.00074318     ...
    ENSG00000237613    0.0030462      0.00027513     0.0030465      0.00031694     ...
    ENSG00000268020    0.011766       0.00061769     0.013409       0.0011429      ...
    ENSG00000240361    0.017913       0.00093294     0.021833       0.001556       ...


A transcript-level expression file looks like this:

.. code:: 

    Name                               Tissue1.mean   Tissue1.SE     Tissue2.mean   Tissue2.SE     ...
    ENST00000373020:ENSG00000000003    35.06          0.52271        35.725         0.66812        ...
    ENST00000494424:ENSG00000000003    0.0034329      0.001209       0.0016207      0.0006441      ...
    ENST00000496771:ENSG00000000003    1.0462         0.019697       1.1043         0.02552        ...
    ENST00000612152:ENSG00000000003    2.5764         0.041124       2.4045         0.043626       ...
    ENST00000614008:ENSG00000000003    0.42826        0.01346        0.41354        0.01551        ...
    ENST00000373031:ENSG00000000005    15.215         0.58333        9.5993         0.49941        ...
    ENST00000485971:ENSG00000000005    1.0715         0.04074        1.1209         0.052269       ...


Examples
--------

DESE based on physical distance (eDESE:dist)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, SNPs inside a gene and its 10 kb adjacent regions will be considered as belonging to a gene. Significant genes by ECS with Bonferroni-adjusted p<0.05 will be retained for fine-mapping. Adjustment for selection bias and multiple testing will be carried out by 100 permutations. 

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --db-gene refgene,gencode \
      --only-hgnc-gene \
      --gene-assoc-condi \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --neargene 10000 \
      --multiple-testing bonf \
      --p-value-cutoff 0.05 \
      --expression-file ../resources/GTEx_v8_TMM_all.gene.meanSE.txt.gz \
      --dese-permu-num 100 \
      --out geneAssoc


DESE guided by eQTLs (eDESE:gene and eDESE:isoform)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform conditional gene-based association analysis using another two different strategies to map variants to genes, i.e., gene-level and isoform-level eQTLs (also are variants). The two strategies correspond to two models, i.e., eDESE:gene and eDESE:isoform, respectively.

eDESE:gene

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --db-gene refgene,gencode \
      --only-hgnc-gene \
      --gene-assoc-condi \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --eqtl-file GTEx_v8_gene_BrainBA9.eqtl.txt.gz \
      --filter-eqtl-p 0.01 \
      --multiple-testing bonf \
      --p-value-cutoff 0.05 \
      --expression-file ../resources/GTEx_v8_TMM_all.gene.meanSE.txt.gz \
      --out geneAssoceQTL

eDESE:isoform

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --db-gene refgene,gencode \
      --only-hgnc-gene \
      --gene-assoc-condi \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --eqtl-file GTEx_v8_transcript_BrainBA9.eqtl.txt.gz \
      --filter-eqtl-p 0.01 \
      --multiple-testing bonf \
      --p-value-cutoff 0.05 \
      --expression-file GTEx_v8_TMM.transcript.meanSE.txt.gz \
      --out geneAssocIsoformeQTL

DESE for drug repositioning
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, ``--expression-file`` specifies a customized file of the drug-induced gene-expression fold-change profile of which the format is the same as :ref:`the gene expression file <expression_file>`. DESE estimates the selective drug perturbation effect on the phenotype-associated genes' expression to aid the drug repositioning for complex diseases.


.. code:: shell

    java -Xmx10g -jar ../kggsee.jar \
      --db-gene refgene \
      --only-hgnc-gene \
      --gene-assoc-condi \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --neargene 5000 \
      --multiple-testing bonf \
      --p-value-cutoff 0.05 \
      --expression-file drug-induced_expression_change_profile \
      --dese-permu-num 100 \
      --out Selective_Perturbed_Drugs


Outputs
-------
The three files with suffixes of ``.gene.pvalue.txt``, ``.gene.var.pvalue.txt.gz``, and ``.qq.png`` are the same as their counterparts output by :ref:`Gene-based association tests <detail_ECS>`.

In addition, results of conditional gene-based association tests are saved in a file with a suffix of ``.finemapping.gene.ecs.txt``. The columns of the file are as follows:

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
      - Start position of the gene
    * - EndPos
      - End position of the gene
    * - #Var
      - Number of variants within the gene
    * - Group
      - LD group number. Conditional ECS tests were performed for genes within the same LD group.
    * - ECSP
      - p-value of ECS
    * - CondiECSP
      - p-value of conditional gene-based association tests by conditional ECS
    * - GeneScore
      - The gene's selective expression score in all tissues. A gene with a high score will be given higher priority to enter the conditioning procedure.


Results of phenotype-tissue associations are saved in a file with a suffix of ``.celltype.txt``. The columns of the file are as follows:

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


If ``--geneset-db`` or ``--geneset-file`` is specified, results of enrichment tests are saved in a file with a suffix of ``.geneset.txt``. The columns of the file are as follows:

.. list-table::
    :widths: 1 2
    :header-rows: 1
    :class: tight-table


    * - Header
      - Description
    * - GeneSet_ID
      - Gene-set ID in the first column of the gene-set file
    * - Enrichment_PValue_Hypergeometric
      - p-values of the hypergeometric tests.
    * - IsSignificant_Hypergeometric
      - If the conditional significant genes are significantly enriched in the gene set.
    * - Total_GeneSet_Gene#
      - The total number of genes in the gene set.
    * - GeneSet_URL
      - Gene-set URL in the second column of the gene-set file
    * - Gene_PValue
      - p-values of conditional significant genes within the gene set.


.. admonition:: Citation of DESE

    Lin Jiang, Chao Xue, Sheng Dai, Shangzhen Chen, Peikai Chen, Pak Chung Sham, Haijun Wang and Miaoxin Li. DESE: estimating driver tissues by selective expression of genes associated with complex diseases or traits. Genome Biology (2019). 20(1):1-19. https://doi.org/10.1186/s13059-019-1801-5


.. admonition:: Citation of eDESE

    Xiangyi Li, Lin Jiang, Chao Xue, Mulin Jun Li and Miaoxin Li. A conditional gene-based association framework integrating isoform-level eQTL data reveals new susceptibility genes for schizophrenia. Elife (2022). 10:e70779. https://doi.org/10.7554/elife.70779


.. admonition:: Citation of the enrichment analysis

    Hongsheng Gui, Johnny S. Kwan, Pak C. Sham, Stacey S. Cherny and Miaoxin Li. Sharing of Genes and Pathways Across Complex Phenotypes: A Multilevel Genome-Wide Analysis. Genetics (2017). 206(3):1601–1609. https://doi.org/10.1534/genetics.116.198150


.. _detail_EMIC:

EMIC (gene-expression causal-effect inference)
==============================================

EMIC (**E**\ ffective-median-based **M**\ endelian randomization framework for **I**\ nferring the **C**\ ausal genes of complex phenotypes) inferences gene expressions' causal effect on a complex phenotype with dependent expression quantitative loci by a robust median-based Mendelian randomization. The effective-median method solved the high false-positive issue in the existing MR methods due to either correlation among instrumental variables or noises in approximated linkage disequilibrium (LD). EMIC can further perform a pleiotropy fine-mapping analysis to remove possible false-positive estimates (`Jiang et al. 2022 <https://doi.org/10.1016/j.ajhg.2022.04.004>`_).


Synopsis
--------

.. code:: shell

    java -Xms16g -Xmx16g -jar kggsee.jar
      --emic
      --out <prefix>
      --vcf-ref <file>
      --sum-file <file>
      --chrom-col <header>  # default: CHR
      --pos-col <header>  # default: BP
      --a1-col <header>  # default: A1
      --a2-col <header>  # default: A2
      --freq-a1-col <header>  # default: FRQ_U
      --beta-col <header>
      --beta-type <0|1|2>
      --se-col <header>  # default: SE
      --eqtl-file <file>
      --filter-eqtl-p <pval>  # default: 1E-4
      --ld-pruning-mr  <r2>  # default: 0.5
      --emic-pfm-p <pval>  # default: 2.5E-6
      --emic-plot-p <pval>  # default: 2.5E-3


When performing EMIC (triggered by ``--emic``), a GWAS summary statistic file (specified by ``--sum-file``) and an eQTL summary statistic file (specified by ``eqtl-file``) are needed. The GWAS summary statistic file must have columns of SNP coordinates (specified by ``--chrom-col`` and ``--pos-col``), the two alleles (specified by ``--a1-col`` and ``--a2-col``), frequencies of the allele specified by ``--a1-col`` (specified by ``--freq-a1-col``), the effect sizes and its standard errors (specified by ``--beta-col`` and ``--se-col``). The type of effect sizes is specified by ``--beta-type`` (``0`` for the linear regression coefficients of a quantitative phenotype; ``1`` for the logarithm of odds ratio or logistic regression coefficient of a qualitative phenotype; ``2`` for an odds ratio of a qualitative phenotype). ``--filter-eqtl-p`` specifies the p-value threshold of eQTLs; only eQTLs with a p-value lower than the threshold will be considered; we note here that the default value is ``1E-4`` for EMIC, which is different from the other analyses. ``--ld-pruning-mr`` specifies the threshold of LD coefficient when pruning variants; for each gene or transcript, eQTLs with LD coefficients higher than the threshold will be pruned. ``--emic-pfm-p`` specifies the p-value threshold to further perform an EMIC pleiotropy fine-mapping (EMIC-PFM) analysis; if the EMIC p-value of a gene is lower than the threshold, an EMIC-PFM will be performed to control the false-positive caused by pleiotropy. ``--emic-plot-p`` specifies the p-value threshold for plotting a scatter plot; genes with an EMIC p-value lower than the threshold will be plotted. :ref:`A description of the eQTL file format <eqtl_file>` is near the beginning of the page.


Examples
--------

Gene-expression causal-effect inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an example of gene-level EMIC. Only eQTLs with a p-value lower than 1E-6 will be considered IVs. Genes with a p-value of EMIC lower than 0.05 will also undergo EMIC-PFM. Genes with a p-value of EMIC lower than 0.01 will be plotted.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --eqtl-file GTEx_v8_gene_BrainBA9.eqtl.txt.gz \
      --beta-col OR \
      --beta-type 2 \
      --emic \
      --filter-eqtl-p 1e-6 \
      --emic-pfm-p 0.05 \
      --emic-plot-p 0.01 \
      --out t3.1


Transcript-expression causal-effect inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an example of transcript-level EMIC. Only eQTLs with a p-value lower than 1E-6 will be considered IVs. Transcripts with a p-value of EMIC lower than 0.05 will also undergo EMIC-PFM. Transcripts with a p-value of EMIC lower than 0.01 will be plotted.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --eqtl-file GTEx_v8_transcript_BrainBA9.eqtl.txt.gz \
      --beta-col OR \
      --beta-type 2 \
      --emic \
      --filter-eqtl-p 1e-6 \
      --emic-pfm-p 0.05 \
      --emic-plot-p 0.01 \
      --out t3.2


Outputs
-------
The numeric results of EMIC are saved in a file with a suffix of ``.emic.gene.txt``. There are nine columns in the file:

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
      - Detailed results of EMIC-PFM separated by semicolons. Each result has four components in brackets: the number of IVs, the causal effect estimate and its standard error, and the p-value. When a transcript-level EMIC is performed, results for each transcript are listed.
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

The numeric results of EMIC-PFM are saved in a file with a suffix of ``.emic.gene.PleiotropyFinemapping.txt``. Only genes with a p-value lower than the threshold specified by ``--emic-pfm-p`` are saved. The file has thirteen columns, of which nine are the same as columns of ``*.emic.gene.txt``. The other four columns are:


.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - Group
      - IDs of a group of genes that share eQTLs.
    * - minP_EMIC_PFM
      - p-value of EMIC-PFM. When a transcript-level EMIC-PFM is performed, this is the minimum p-value among all transcripts of the gene.
    * - DetailsEMIC_PFM
      - Detailed results of EMIC-PFM separated by semicolons. Each result has four components in brackets: the number of IVs, the causal effect estimate and its standard error, and the p-value. When a transcript-level EMIC-PFM is performed, results for each transcript are listed.
    * - CochransQ
      - The p-value of an extended Cochran's Q test. The significance (p<1E-3) means that the causal effect is more likely to be false-positive. At this point, KGGSEE excludes its eQTLs which are also the eQTLs of other significant genes, and redoes EMIC. In this case, results in the columns of minP_EMIC_PFM and DetailsEMIC_PFM will be different from those in the columns of minP_EMIC and Details_EMIC.


The columns of the file with a suffix of ``.emic.gene.var.tsv.gz`` are the same as ``*.emic.gene.txt``. The difference is that, for each gene, in ``*.emic.gene.txt``, only the eQTL with the lowest GWAS p-value is output, while in ``*.emic.gene.var.tsv.gz``, all eQTLs are output. The file with a suffix of ``.qq.png`` saves the Q-Q plot for GWAS p-values of IVs. The file with a suffix of ``.emic.qq.png`` saves the Q-Q plot for EMIC p-values. The file with the suffix of ``.scatterplots.emic.pdf`` saves the scatter plots of genetic association with gene expression. Each gene with an EMIC p-value lower than the threshold specified by ``--emic-plot-p`` is saved on a separate page of the PDF. A filled rectangle on the plots denotes an IV. The red rectangle denotes the most significant GWAS variant among all the IVs of a gene. The slope of the line represents the estimated causal effect. The color of an IV denotes the degree of the LD between the IV and the most significant GWAS variant. The error bars in the rectangles denote the standard errors of the coefficient estimates. The file with the suffix of ``..scatterplots.emic.txt`` saves the numeric results of the scatter plots. For each gene, the first line is the gene ID after ``>`` which is followed the estimated causal effect size (``Eff.=``) and the p-value (``p=``). Then, a table of the information of each IV is shown:


.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - BP
      - The coordinate of the IV
    * - ExpBeta
      - The effect size of the IV to the expression
    * - ExpBetaSE
      - The SE of ExpBeta
    * - TraitBeta
      - The effect size of the IV to the trait. It is the linear regression coefficient for a quantitative trait; it is the logistic regression coefficient or the natural logarithm of the odds ratio for a dichotomous trait, which is depending on the value of ``--beta-type``.
    * - TraitBetaSE
      - The SE of TraitBeta


After the table, the LD coefficient matrix of the IVs is shown.



.. admonition:: Citation of EMIC

    Lin Jiang, Lin Miao, Guorong Yi, Xiangyi Li, Chao Xue, Mulin Jun Li, Hailiang Huang and Miaoxin Li. Powerful and robust inference of causal genes of complex phenotypes with dependent expression quantitative loci by a novel median-based Mendelian randomization. The American Journal of Human Genetics (2022). 109(5):838-856. https://doi.org/10.1016/j.ajhg.2022.04.004


.. _detail_h2:

EHE (gene-based heritability estimation)
========================================

This analysis estimates the heritability of each gene and performs gene-based association tests at the same time (`Miao et al. 2023 <https://doi.org/10.1016/j.ajhg.2023.08.006>`_).


Synopsis
--------

.. code:: shell

    java -Xms16g -Xmx16g -jar kggsee.jar
      --gene-herit
      --out <prefix>
      --vcf-ref <file>
      --sum-file <file>
      --chrom-col <header>  # default: CHR
      --pos-col <header>  # default: BP
      --p-col <header>  # default: P
      --nmiss-col <header>
      --case-col <header>
      --control-col <header>
      --prevalence <value>  # default: 0.01
      --neargene <basepair>  # default: 5000
      --eqtl-file <file>
      --filter-eqtl-p <pval>  # default: 0.01
      --gene-assoc-condi


``--gene-herit`` triggers gene-based association tests and estimation of gene heritability. ``--sum-file`` specifies a white space delimited GWAS summary statistic file which must have three columns of the chromosome of SNP, coordinate of SNP, and p-value of SNP; headers of the three columns can be specified by ``--chrom-col``, ``--pos-col`` and ``--p-col`` separately. In addition, for quantitative phenotype, a column of sample sizes is needed, and its header is specified by ``--nmiss-col``; for qualitative phenotype, two columns of case sample sizes and control sample sizes are needed, and their header is specified by ``--case-col`` and ``--control-col`` separately. SNPs belonging to a gene can be defined either by SNPs close to the gene or by eQTLs of the gene. If ``--neargene`` is specified, KGGSEE reads gene annotations and considers SNPs inside a gene and its adjacent regions at a fixed number of basepairs on both sides to be a test unit. If ``--eqtl-file`` is specified, KGGSEE reads the eQTL summary statistic file and considers eQTLs of a gene or a transcript to be a test unit, and ``--neargene`` is overridden. When ``--eqtl-file`` is specified, ``--filter-eqtl-p`` can be used to specify a threshold of eQTL p-values. Only eQTLs with a p-value lower than the threshold will be considered. When ``--gene-assoc-condi`` is specified, KGGSEE also calculates the conditional heritability of genes, and the flags of ``--multiple-testing``, ``--p-value-cutoff``, ``--top-gene`` and ``--expression-file`` have the same meaning as in :ref:`DESE <detail_DESE>`. :ref:`A description of the eQTL file format <eqtl_file>` is near the beginning of the page.


Examples
--------

Gene heritability based on physical distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, SNPs inside a gene and its 10 kb adjacent regions will be grouped to estimate heritability. The prevalence of affected individuals is set to 0.01.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --gene-herit \
      --prevalence 0.01 \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --case-col Nca \
      --control-col Nco \
      --neargene 10000 \
      --out t4.1

.. note::
    When ``--case-col`` and ``--control-col`` are specified, KGGSEE will regard the input as summary statistics from case/control samples and automatically adjust for the disease prevalence. On the other hand, if the ``--nmiss-col`` is specified, KGGSEE will regard the input as summary statistics for a quantitative trait.


Gene heritability guided by eQTLs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, eQTLs of a gene will be grouped to estimate heritability.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --gene-herit \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --case-col Nca \
      --control-col Nco \
      --eqtl-file GTEx_v8_gene_BrainBA9.eqtl.txt.gz \
      --out t4.2


Transcript heritability guided by eQTLs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, eQTLs of a transcript will be grouped to estimate heritability.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --gene-herit \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --case-col Nca \
      --control-col Nco \
      --eqtl-file GTEx_v8_transcript_BrainBA9.eqtl.txt.gz \
      --out t4.3
    

Gene's conditional heritability based on physical distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, SNPs inside a gene and its 10 kb adjacent regions will be grouped to estimate heritability. Significant genes by ECS with Bonferroni-adjusted p<0.05 will be retained for fine-mapping and then calculating conditional heritability.

.. code:: shell

    java -Xmx4g -jar ../kggsee.jar \
      --gene-herit \
      --prevalence 0.01 \
      --vcf-ref 1kg_hg19_eur_chr1.vcf.gz \
      --sum-file scz_gwas_eur_chr1.tsv.gz \
      --case-col Nca \
      --control-col Nco \
      --neargene 10000 \
      --multiple-testing bonf \
      --p-value-cutoff 0.05 \
      --expression-file ../resources/GTEx_v8_TMM_all.gene.meanSE.txt.gz \
      --out t4.4


Outputs
-------
The file with a suffix of ``.gene.pvalue.txt`` saves the results of gene-based heritability estimates and association tests. The columns of the file are as follows:


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
    * - Herit
      - Heritability estimate
    * - HeritSE
      - Standard error of the heritability estimate
    * - Chrom
      - Chromosome of the gene
    * - Pos
      - The coordinate of the variant with the lowest p-value within the gene
    * - GWAS_Var_P
      - p-value of the variant


The columns of the file with the suffix of ``.gene.var.pvalue.txt.gz`` are the same as ``*.gene.pvalue.txt``. The difference is that, for each gene, in ``*.gene.pvalue.txt``, only the variant with the lowest p-value is output, while in ``*.gene.var.pvalue.txt.gz``, all variants are output. The file with the suffix of ``.qq.png`` is the Q-Q plots for p-values of GWAS summary statistics and gene-based association tests by GATES and ECS.


When ``--gene-assoc-condi`` is specified, a file with a suffix of ``.finemapping.gene.ecs.txt`` is also output. This file has the following four more columns in addition to its counterpart output by :ref:`DESE <detail_DESE>`.

.. list-table::
    :widths: 1 4
    :header-rows: 1
    :class: tight-table

    * - Header
      - Description
    * - Herit
      - Unconditional heritability estimate
    * - HeritSE
      - Standard error of the unconditional heritability estimate
    * - CondiHerit
      - Conditional heritability estimate
    * - CondiHeritSE
      - Standard error of the conditional heritability estimate


.. admonition:: Citation of EHE

    Lin Miao, Lin Jiang, Bin Tang, Pak Chung Sham and Miaoxin Li. Dissecting the high-resolution genetic architecture of complex phenotypes by accurately estimating gene-based conditional heritability. The American Journal of Human Genetics (2023). 110(9):1534–1548. https://doi.org/10.1016/j.ajhg.2023.08.006

