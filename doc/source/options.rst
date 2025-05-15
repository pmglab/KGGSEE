.. _options:

=======
Options
=======

The options for :ref:`Reference population genotypes <option_vcf>`, :ref:`GWAS summary statistics <option_gwas>`, and :ref:`Global misc <option_misc>` act on all the analyses. For clarity, we have categorized the other parameters by :ref:`GATES and EHE <option_assoc>`, :ref:`DESE <option_dese>`, :ref:`EMIC <option_emic>` and :ref:`EHE <option_h2>`, although this has resulted in duplications.

In the "Default" columns of the following tables, "null" denotes that the flag works with an argument but there is no default value; "n/a" denotes that the flag works without any argument.


.. _option_vcf:

Reference population genotypes
==============================

These options work on the VCF file of reference population genotypes. Only SNPs that pass the filters will be used for subsequent analyses.


.. list-table:: 
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--vcf-ref``
      - Specifies a VCF file of genotypes sampled from the target population. These genotypes are used to estimate LD correlation coefficients between any pair of SNPs. For VCF files of separated chromosomes, use wildcards and quotes like ``"chr*.vcf.gz"``.
      - null
    * - ``--keep-ref``
      - Specifies a directory to keep the parsed VCF files in KGGSEE object format.
      - null
    * - ``--saved-ref``
      - Specifies the directory of the genotypes of the reference population in KGGSEE object format, which was saved by ``--keep-ref``. Reading KGGSEE object format files is faster than parsing VCF files.
      - null
    * - ``--filter-maf-le``
      - Filter SNPs with a minor allele frequency lower than the setting.
      - ``0.05``
    * - ``--hwe-all``
      - Filter SNPs with a p-value of rejecting Hardy-Weinberg equilibrium lower than the setting.
      - ``1E-5``
    * - ``--chrom``
      - Specify chromosome labels. e.g., ``--chrom NC_052532.1,NW_024095932.1,NW_024095933.1,NW_024095934.1``.
      - 1,...,22,X,Y,M
    * - ``--ld-block-max-r2``
      - Set the max tolerable LD coefficient between SNPs from two LD blocks. KGGSEE divides SNPs within a genomic region into LD blocks to improve computational efficiency. Any pairwise LD coefficient between SNPs from two LD blocks are always less than the number specified by ``--ld-block-max-r2``. A smaller number leads to larger LD blocks and in turn more RAM and CPU time usage; a larger number may increase the sampling error of LD coefficient estimates. We recommand that only increase the number when the partitioned LD blocks are too big to be analyzed.
      - ``0.15``
    * - ``ld-block-win``
      - Set the window size in base-pair for searching LD blocks on chromosomes
      - ``5000000``


.. _option_gwas:

GWAS summary statistics
=======================


.. list-table:: 
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--sum-file``
      - Specifies a whitespace delimitated file of GWAS summary statistics.
      - null
    * - ``--chrom-col``
      - Specifies the column header of chromosomes. 
      - ``CHR``
    * - ``--pos-col``
      - Specifies the column header of coordinates.
      - ``BP``
    * - ``--p-col``
      - Specifies the column header of p-values.
      - ``P``
    * - ``--a1-col``
      - Specifies the column header of the effect allele.
      - ``A1``
    * - ``--a2-col``
      - Specifies the column header of the other allele.
      - ``A2``
    * - ``--freq-a1-col``
      - Specifies the column header of the frequencies of the allele specified by ``--a1-col``.
      - ``FRQ_U``
    * - ``--beta-col``
      - Specifies the column header of effect sizes.
      - null
    * - ``--beta-type``
      - Specifies the type of effect sizes:  ``0`` for the linear regression coefficient of a quantitative trait; ``1`` for the logarithm of odds ratio or logistic regression coefficient of a dichotomous trait; ``2`` for an odds ratio of a dichotomous trait.
      - null
    * - ``--se-col``
      - Specifies the column header of standard errors of effect sizes. Note: even if the effect size is provided as an odds ratio, this is still the standard error of the logarithm (base e) of the odds ratio.
      - ``SE``
    * - ``--nmiss-col``
      - Specifies the column header of sample sizes for a quantitative trait.
      - null
    * - ``--case-col``
      - Specifies the column header of case sample sizes for a qualitative trait.
      - null
    * - ``--control-col``
      - Specifies the column header of control sample sizes for a dichotomous trait.
      - null
    * - ``--adjust-gc``
      - Ask KGGSEE to adjust the p-values and chi-square statistics using the genomic control factors from the input GWAS data before all follow-up analyses.
      - 1


.. _option_misc:

Global misc
===========


.. list-table::
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--nt``
      - Specifies the number of threads.
      - ``4``
    * - ``--lib-update``
      - Download ``kggsee.jar`` from http://pmglab.top/kggsee and replace the current running one.
      - n/a
    * - ``--buildver``
      - Specifies the reference genome version of the coordinates. The supported versions are ``hg19`` and ``hg38``.
      - ``hg19``
    * - ``--db-gene``
      - Specifies the database of gene annotations. ``refgene`` for RefSeq Genes; ``gencode`` for GENCODE; ``refgene,gencode`` for both.
      - ``refgene``
    * - ``--excel``
      - Output results in Excel format.
      - n/a
    * - ``--only-hgnc-gene``
      - Only genes with an HGNC-approved gene symbol are considered in the analysis.
      - n/a
    * - ``--out``
      - Specifies the output prefix of results.
      - ``./kggsee1``
    * - ``--regions-bed``
      - Specify a `BED file <https://en.wikipedia.org/wiki/BED_(file_format)>`_ to define customized gene coordinates instead of the annotation from RefSeqGene or GENCODE. The first three columns of the BED file define gene coordinates and are mandatory; the fourth column defines gene names and is optional. When the fourth column is absent, a gene name of the format like ``chr1:100-200`` will be allocated.
      - null
    * - ``--regions-out``
      - Specifies genomic regions to be excluded in the analysis, e.g., ``chr1,chr2:2323-34434,chr2:43455-345555``. 
      - null
    * - ``--resource``
      - Specifies the path KGGSEE running resource data.
      - ``path/to/kggsee.jar/resources/``


.. _option_assoc:

GATES and ECS
=============


.. list-table::
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--gene-assoc``
      - Triggers gene-based association tests.
      - n/a
    * - ``--neargene``
      - One number sets the basepair to extend at both sides of a gene, when considering SNPs belonging to the gene, e.g., ``--neargene 5000``. This flag can also have two values to set an asymmetric boundary extension, e.g., 5 kb upstream and 15 kb downstream of a gene can be set by ``--neargene 5000,15000``.
      - ``5000``
    * - ``--eqtl-file``
      - Specifies a fasta-styled file of eQTL summary statistics. If this flag is used, ``--neargene`` is overridden, and eQTLs of a gene or transcript will be grouped and tested.
      - null
    * - ``--filter-eqtl-p``
      - Specifies the threshold of eQTL p-values. Only eQTLs with a p-value lower than the threshold will be used. The default is ``0.01`` when performing gene-based association tests and heritability estimation.
      - ``0.01``


.. _option_dese:

DESE
====


.. list-table::
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--gene-assoc-condi``
      - Trigers the DESE, eDESE or SelDP.
      - n/a
    * - ``--expression-file``
      - Specifies a gene expression file that contains means and standard errors of gene expressions in multiple tissues.
      - null
    * - ``--multiple-testing``
      - Specifies the method for multiple testing correction. ``bonf`` denotes performing Bonferroni correction; ``benfdr`` denotes controlling false discovery rate by the Benjamini–Hochberg method; ``fixed`` denotes no correction.
      - ``bonf``
    * - ``--p-value-cutoff``
      - Specifies the threshold of the adjusted p-value for fine-mapping. Only genes with an adjusted p-value lower than the threshold will be retained for fine-mapping.
      - 0.05
    * - ``--top-gene``
      - Specifies the maximum number of genes with the smallest p-values that will be retained for fine-mapping.
      - null
    * - ``--geneset-db``
      - Specifies `MSigDB <http://www.gsea-msigdb.org/gsea/msigdb/index.jsp>`_ gene sets for enrichment analysis:
        
        ``cura``: C2. curated gene sets;
        
        ``cgp``: C2. chemical and genetic perturbations;
        
        ``cano``: C2. canonical pathways;
        
        ``cmop``: C4. computational gene sets;
        
        ``onto``: C5. ontology gene sets;
        
        ``onco``: C6. oncogenic signature gene sets;
        
        ``immu``: C7. immunologic signature gene sets.
      - null
    * - ``--geneset-file``
      - Specifies a user-defined file of gene sets for enrichment analysis.
      - null
    * - ``--neargene``
      - One number sets the basepair to extend at both sides of a gene when considering SNPs belonging to the gene, e.g., ``--neargene 5000``. This flag can also have two values to set an asymmetric boundary extension, e.g., 5 kb upstream and 15 kb downstream of a gene can be set by ``--neargene 5000,15000``.
      - ``5000``
    * - ``--eqtl-file``
      - Specifies a fasta-styled file of eQTL summary statistics. If this flag is used, ``--neargene`` is overridden, and eQTLs of a gene or transcript will be grouped and tested.
      - null
    * - ``--filter-eqtl-p``
      - Specifies the threshold of eQTL p-values. Only eQTLs with a p-value lower than the threshold will be used. The default is ``0.01`` when performing DESE.
      - ``0.01``
    * - ``--dese-permu-num``
      - The number of permutations for an adjustment of selection bias and multiple testing
      - null


.. _option_emic:

EMIC
====


.. list-table::
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--emic``
      - Triggers the EMIC.
      - n/a
    * - ``--eqtl-file``
      - Specifies a fasta-styled file of eQTL summary statistics.
      - null
    * - ``--filter-eqtl-p``
      - Specifies the threshold of eQTL p-values. Only eQTLs with a p-value lower than the threshold will be used. The default is ``1E-4`` when performing EMIC.
      - ``1E-4``
    * - ``--ld-pruning-mr``
      - Specifies the threshold of LD coefficients when pruning variants. For each gene or transcript, eQTLs with LD coefficients higher than the threshold will be pruned.
      - 0.5
    * - ``--emic-pfm-p``
      - Specifies the p-value threshold to further perform an EMIC pleiotropy fine-mapping (EMIC-PFM) analysis. If the EMIC p-value of a gene is lower than the threshold, an EMIC-PFM will be performed to control the false-positive caused by pleiotropy. 
      - ``2.5E-6``
    * - ``--emic-plot-p``
      - Specifies the p-value threshold for plotting a scatter plot. Genes with an EMIC p-value lower than the threshold will be plotted.
      - ``2.5E-3``      


.. _option_h2:

EHE
===


.. list-table::
    :widths: 3 8 2
    :header-rows: 1
    :class: tight-table


    * - Flag
      - Description
      - Default
    * - ``--gene-herit``
      - Triggers gene-based association tests and estimation of gene heritability. The flags of ``--neargene``, ``--eqtl-file`` and ``--filter-eqtl-p`` have the same meaning as in :ref:`GATES and ECS <option_assoc>`.
      - n/a
    * - ``--case-col``, ``--control-col``, ``--nmiss-col``
      - When ``--case-col`` and ``--control-col`` are specified, KGGSEE will regard the input as summary statistics from case/control samples and automatically adjust for the disease prevalence. On the other hand, if the ``--nmiss-col`` is specified, KGGSEE will regard the input as summary statistics for a quantitative trait.
      - null
    * - ``--gene-assoc-condi``
      - When ``--gene-assoc-condi`` is specified in addition to ``--gene-herit``, KGGSEE also calculates the conditional heritability, and the flags of ``--multiple-testing``, ``--p-value-cutoff``, ``--top-gene`` and ``--expression-file`` have the same meaning as in :ref:`DESE <option_dese>`.
      - n/a
    * - ``--prevalence``
      - Specifies the proportion of cases in the population when estimating the heritability of a dichotomous trait.
      - 0.01
