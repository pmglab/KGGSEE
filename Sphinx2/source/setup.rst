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

KGGSEE needs JRE 1.8 or higher to run. Both `Java(TM) SE JRE <https://java.com/en/download/manual.jsp>`_ and `OpenJDK JRE <https://openjdk.java.net/install>`_ are competent for KGGSEE. Please follow the instructions on the websites to complete the installation and also add Java to the system PATH.

Check the JRE by entering ``java -version`` in a Terminal of Linux or MacOS, or CMD or PowerShell of Windows. If it displays the JRE version like ``Java(TM) SE Runtime Environment (build 1.8.0_xxx)`` or ``OpenJDK Runtime Environment (build 1.8.0_xxx)``, it means the JRE has already been set up. Otherwise, check if JRE has been installed and if Java is in the system PATH.


Setup KGGSEE
============

1. Download the `KGGSEE Java archive <http://pmglab.top/kggsee/download/lib/v1/kggsee.jar>`_, and save it as ``kggsee.jar``)

2. Download the `running resource dataset <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/ESO0a08VRuFPq7YOxEdsMU8BaEX9c-wbhyBP537Qdg-EKw?e=tyQhhz>`_, unzip and save the folder as ``resources/`` (in the same folder with ``kggsee.jar``), the dataset includes:


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
 

3. For running Quick tutorials, download the `tutorial dataset <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EcQbomHELXVNrKdxqWvzOg4BVYoqsYpLU-wWy3_0ewgHiA?e=CRjYmq>`_, unzip and save the folder as ``tutorials/`` (in the same folder with ``kggsee.jar``), the dataset includes:


.. list-table::
   :widths: 1 1
   :header-rows: 0
   :class: tight-table
   
   * - ``tutorials/scz_gwas_eur_chr1.tsv.gz``
     - Summary statistics of chr1 SNPs for a GWAS of schizophrenia on the European population.
   
   * - ``tutorials/1kg_hg19_eur_chr1.vcf.gz``
     - Genotypes of chr1 SNPs sampled from the 1000 Genome Project European population.
   
   * - ``tutorials/GTEx_v8_gene_meanSE.tsv.gz``
     - The gene-leveled expression profile of GTEx v8 tissues.
   
   * - ``tutorials/GTEx_v8_transcript_meanSE.tsv.gz``
     - The transcript-leveled expression profile of GTEx v8 tissues.
   
   * - ``tutorials/GTEx_v8_gene_BrainBA9.eqtl.txt.gz``
     - Summary statistics of eQTLs calculated from gene-leveled expression profile of GTEx v8 brain BA9.

   * - ``tutorials/GTEx_v8_transcript_BrainBA9.eqtl.txt.gz``
     - Summary statistics of eQTLs calculated from transcript-leveled expression profile of GTEx v8 brain BA9.


4. For running customized analyses, the following data is needed, refer to :ref:`Detailed Document <detailed_document>` for descriptions of file formats.

   * GWAS summary statistics of the phenotype to be studied

   * VCF files of genotypes sampled from the population of the GWAS to be studied. Genotypes of the 1000 Genomes Project Phase3 v5 can be downloaded from the `NCBI FTP site <ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502>`_ or the `1000 Genomes Project FTP site <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502>`_.

   * If running DESE, a file of multiple-tissue gene expression profile is needed. We provide tissue-leveled and cell type-leveled expression profiles available for download on `OneDrive <https://mailsysueducn-my.sharepoint.com/personal/limiaoxin_mail_sysu_edu_cn/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Flimiaoxin%5Fmail%5Fsysu%5Fedu%5Fcn%2FDocuments%2Ftools%2Fkggsee%2Fresources&ga=1>`_.

   * For all analyses, a file of eQTL summary statistics calculated from target tissues may be used. We provide gene-based and transcript-based eQTL summary statistics for GTEx v8 tissues available for download on `OneDrive <https://mailsysueducn-my.sharepoint.com/personal/limiaoxin_mail_sysu_edu_cn/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Flimiaoxin%5Fmail%5Fsysu%5Fedu%5Fcn%2FDocuments%2Ftools%2Fkggsee%2Fresources&ga=1>`_.

