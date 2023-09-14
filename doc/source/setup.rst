.. _setup:

=====
Setup
=====


System requirements
===================

.. list-table::
    :widths: 3 7
    :header-rows: 0
    :class: tight-table

    * - Operating system
      - KGGSEE runs in a Java Virtual Machine. It does not matter which operating system it runs in.
    * - Java Runtime Environment
      - A Java SE Runtime Environment of version 1.8 or higher is needed.
    * - CPU
      - A CPU with four cores or more is recommended.
    * - Memory
      - 16 GB RAM or higher is recommended.
    * - Free space
      - KGGSEE and related datasets may take up to 10 GB.


Setup the Java Runtime Environment (JRE)
========================================

KGGSEE needs JRE 1.8 or higher to run. Both `Java(TM) SE JRE <https://java.com/en/download/manual.jsp>`_ and `OpenJDK JRE <https://openjdk.java.net/install>`_ are competent for KGGSEE.

After installing a JRE, check by entering ``java -version`` in a Terminal of Linux/MacOS, or a CMD/PowerShell of MS Windows. If it displays the JRE version like ``Java(TM) SE Runtime Environment (build 1.8.0_xxx)`` or ``OpenJDK Runtime Environment (build 1.8.0_xxx)``, it means the JRE has already been set up. Otherwise, check if JRE has been installed and if ``java`` is in ``$PATH``.


KGGSEE and its running resources
================================

.. list-table::
    :widths: 2 5 2 2
    :header-rows: 1
    :class: tight-table

    * - File
      - Description
      - Size
      - Updated
    * - `kggsee.jar <https://pmglab.top/kggsee/download/lib/v1/kggsee.jar>`_
      - The KGGSEE program
      - 46 MB
      - Sep 2023
    * - `resources/ <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EpXRqLXIToZItErUHiDNDO0BM29gbEn1-Grs14D_EqORJQ?e=0ZjvlN>`_
      - All running resource files in our OneDrive
      - 7.0 GB
      - Sep 2023
    * - `resources.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EYhQXE95WZFMqERo_xNOhZUB8lGeyTwPuiWM26AX8CHP8Q?e=PwbMoa>`_
      - Running resource files except for reference genotypes and eQTL summary statistics 
      - 362 MB
      - Sep 2023
    * - `tutorials.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EWqZHY25tT5Nq1GMwtl06ocBHoTAXGyBTH74zAp68dv5VA?e=tPtZ7B>`_
      - A tutorial dataset to run through :ref:`the four analyses <four_analyses>`
      - 155 MB
      - Apr 2022


KGGSEE is written in Java and distributed as a Java Archive ``kggsee.jar``. To perform an analysis, corresponding resource files are needed. For example, reference genotypes and gene annotations are needed for gene-based association test and heritability estimation; in addition, eQTL summary statistics is needed for gene-expression causal effect estimation.

A quick and easy way to set up an environment for :ref:`Quick tutorials <quick_tutorials>` is to (1) download ``kggsee.jar``, ``resources.zip`` and ``tutorials.zip``; (2) unzip ``resources.zip`` and ``tutorials.zip``; (3) put ``kggsee.jar``, ``resources/`` and ``tutorials/`` under one directory.

A straightforward way to set up an environment for customized analyses is to (1) download ``kggsee.jar`` and ``resources.zip``; (2) unzip ``resources.zip``, and put ``kggsee.jar`` and ``resources/`` under one directory; (3) download the reference genotypes of the population that matches your GWAS; (4) for gene/transcript expression causal effect estimation (EMIC), also download the eQTL summary statistics of phenotype-related tissues. To prepare customized resource files, refer to :ref:`Detailed Document <detailed_document>` for descriptions of the file formats.

The files in `resources/ <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EpXRqLXIToZItErUHiDNDO0Bk-jeiAtIlA-abGjOCdbqEw?e=3Jhy5g>`_ are described in the following, where gene annotations, genotypes of the 1000 genomes project and eQTL summary statistics are provided by both hg19 and hg38 coordinates.

.. list-table::
    :widths: 1 1
    :header-rows: 0
    :class: tight-table

    * - ``resources/{hg19,hg38}/kggseqv1.1_{hg19,hg38}_GEncode.txt.gz``
      - `GENCODE <https://www.gencodegenes.org>`_ annotations
    * - ``resources/{hg19,hg38}/kggseqv1.1_{hg19,hg38}_refGene.txt.gz``
      - `RefGene <https://www.ncbi.nlm.nih.gov/refseq/rsg>`_ annotations
    * - ``resources/HgncGene.txt.gz``
      - `HGNC <https://www.genenames.org>`_ gene ID
    * - ``resources/ENSTGene.gz``
      - `Ensembl <https://www.ensembl.org/index.html>`_ gene ID and transcript ID
    * - ``resources/*.symbols.gmt.gz``
      - `MSigDB <http://www.gsea-msigdb.org/gsea/msigdb/index.jsp>`_ gene sets
    * - ``resources/GTEx_v8_TMM.gene.meanSE.txt.gz``
      - The gene-level expression profile of the `GTEx v8 <https://www.gtexportal.org/home/>`_ tissues
    * - ``resources/GTEx_v8_TMM.transcript.meanSE.txt.gz``
      - The transcript-level expression profile of the `GTEx v8 <https://www.gtexportal.org/home/>`_ tissues
    * - ``resources/{hg19,hg38}/gty/1KG.{AFR,AMR,EAS,EUR,SAS}.{hg19,hg38}.vcf.gz`` (not included in ``resources.zip``)
      - VCF files of each super population panel of `1000 Genomes Project <https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/>`_ including biallelic variants with MAF>0.01
    * - ``resources/{hg19,hg38}/eqtl/*.{gene,transcript}.{hg19,hg38}.cov.eqtl.txt.gz`` (not included in ``resources.zip``)
      - cis-eQTL summary statistics calculated from the gene or transcript-level expression profile of the GTEx v8 dataset


On top of already downloading ``resources.zip``, ``tutorials/scz_gwas_eur_chr1.tsv.gz`` and ``tutorials/1kg_hg19_eur_chr1.vcf.gz`` are sufficient to go through gene-based association test (ECS and GATES) and heritability estimation (EHE). In addition, ``tutorials/GTEx_v8_gene_BrainBA9.eqtl.txt.gz`` is needed for driver-tissue estimation (DESE); ``tutorials/GTEx_v8_{gene,transcript}_BrainBA9.eqtl.txt.gz`` is needed for gene/transcript expression causal effect estimation (EMIC).

.. list-table::
    :widths: 1 1
    :header-rows: 0
    :class: tight-table
    
    * - ``tutorials/scz_gwas_eur_chr1.tsv.gz``
      - Summary statistics of chr1 SNPs for a GWAS of schizophrenia in the European population
    * - ``tutorials/1kg_hg19_eur_chr1.vcf.gz``
      - Genotypes of the chr1 SNPs of the 1000 Genomes Project European panel
    * - ``tutorials/GTEx_v8_gene_BrainBA9.eqtl.txt.gz``
      - eQTL summary statistics calculated from the gene-level brain BA9 expression profile of the GTEx v8 dataset
    * - ``tutorials/GTEx_v8_transcript_BrainBA9.eqtl.txt.gz``
      - eQTL summary statistics calculated from the transcript-level brain BA9 expression profile of the GTEx v8 dataset

