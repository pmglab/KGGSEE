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

KGGSEE needs JRE 1.8 or higher to run. Both `Java(TM) SE JRE <https://java.com/en/download/manual.jsp>`_ and `OpenJDK JRE <https://openjdk.java.net/install>`_ are competent for KGGSEE. After installing a JRE, check by entering ``java -version`` in a Terminal of Linux/MacOS, or CMD/PowerShell of MS Windows. If it displays the JRE version like ``Java(TM) SE Runtime Environment (build 1.8.0_xxx)`` or ``OpenJDK Runtime Environment (build 1.8.0_xxx)``, it means the JRE has already been set up. Otherwise, check if JRE has been installed and if Java is in the system PATH.


Setup KGGSEE
============

KGGSEE is written in Java and distributed as a Java Archive ``kggsee.jar``. In addition, resource datasets, such as gene annotations and eQTL summary statistics, are needed to perfrom the corresponding analyses. A tutorial dataset is also prepared to go through all four KGGSEE analyses.

From `the download page <http://pmglab.top/kggsee/#/download>`_, download ``kggsee.jar`` and ``kggsee_resource_tutorial.zip`` and unzip ``kggsee_resource_tutorial.zip``. Put ``kggsee.jar``, ``resources/`` and ``tutorials/`` in the same directory, and then, it's ready.


Files in ``resources/`` are:

.. list-table::
    :widths: 1 1
    :header-rows: 0
    :class: tight-table

    * - ``resources/hg19/kggseqv1.1_hg19_GEncode.txt.gz``
      - hg19 `GENCODE <https://www.gencodegenes.org>`_ annotation
    * - ``resources/hg19/kggseqv1.1_hg19_refGene.txt.gz``
      - hg19 `RefGene <https://www.ncbi.nlm.nih.gov/refseq/rsg>`_ annotation
    * - ``resources/hg38/kggseqv1.1_hg38_GEncode.txt.gz``
      - hg38 `GENCODE <https://www.gencodegenes.org>`_ annotation
    * - ``resources/hg38/kggseqv1.1_hg38_refGene.txt.gz``
      - hg38 `RefGene <https://www.ncbi.nlm.nih.gov/refseq/rsg>`_ annotation
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


Files in ``tutorials/`` are:

.. list-table::
    :widths: 1 1
    :header-rows: 0
    :class: tight-table
    
    * - ``tutorials/scz_gwas_eur_chr1.tsv.gz``
      - Summary statistics of chr1 SNPs for a GWAS of schizophrenia on the European population
    * - ``tutorials/1kg_hg19_eur_chr1.vcf.gz``
      - Genotypes of chr1 SNPs sampled from the 1000 Genome Project European population
    * - ``tutorials/GTEx_v8_gene_BrainBA9.eqtl.txt.gz``
      - Summary statistics of cis-eQTLs calculated from the gene-level brain BA9 expression profile of the GTEx v8 dataset
    * - ``tutorials/GTEx_v8_transcript_BrainBA9.eqtl.txt.gz``
      - Summary statistics of cis-eQTLs calculated from the transcript-level brain BA9 expression profile of the GTEx v8 dataset


For running customized analyses, the following data is needed, refer to :ref:`Detailed Document <detailed_document>` for descriptions of file formats.

* A file of GWAS summary statistics of the phenotype to be studied.
* VCF files of genotypes sampled from the population of the GWAS to be studied.
* A file of eQTL summary statistics calculated from target tissues may be used.

VCF files of the 1000 Genomes Project and eQTL summary statistics of GTEx v8 tissues are available in the ``resources/`` folder in `the download page <http://pmglab.top/kggsee/#/download>`_: 

.. list-table::
    :widths: 1 1
    :header-rows: 0
    :class: tight-table
    
    * - ``resources/{hg19,hg38}/*.vcf.gz``
      - Genotypes of super population panels of the 1000 Genomes Project including biallelic variants with MAF>0.01.
    * - ``resources/{hg19,hg38}/eqtl_gene/*.gene.{hg19,hg38}.cov.eqtl.txt.gz``
      - Summary statistics of cis-eQTLs calculated from the gene-level expression profile of the GTEx v8 dataset
    * - ``resources/{hg19,hg38}/eqtl_transcript/*.transcript.{hg19,hg38}.cov.eqtl.txt.gz``
      - Summary statistics of cis-eQTLs calculated from the transcript-level expression profile of the GTEx v8 dataset
