.. _setup:

=====
Setup
=====


System requirements
===================

.. list-table::
    :widths: 1 2
    :header-rows: 1
    :class: tight-table

    * - Hardware/Software
      - Requirement
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


Setup a Java Runtime Environment (JRE)
======================================

KGGSEE needs JRE 1.8 or higher. Both `Java(TM) SE JRE <https://java.com/en/download/manual.jsp>`_ and `OpenJDK JRE <https://openjdk.java.net/install>`_ are competent.

After installing a JRE, check by entering ``java -version`` in a Terminal of Linux/MacOS, or a CMD/PowerShell of MS Windows. If it displays the JRE version like ``Java(TM) SE Runtime Environment (build x)`` or ``OpenJDK Runtime Environment (build x)``, it means the JRE has already been set up. Otherwise, check if JRE has been installed and if ``java`` is in ``$PATH``.


KGGSEE and its running resources
================================

KGGSEE is written in Java and distributed as a Java Archive ``kggsee.jar``. To perform an analysis, corresponding running resources are also needed.  For example, reference genotypes and gene annotations are needed for gene-based association tests (GATES and ECS) and heritability estimations (EHE); in addition, eQTL summary statistics are needed for gene-expression causal-effect estimations (EMIC). Thus, ``kggsee.jar`` is always needed and which resource files are needed depends on the analysis. We provide the following download links.

.. list-table::
    :widths: 1 2 1 
    :header-rows: 1
    :class: tight-table

    * - File
      - Description
      - Size
    * - `kggsee.jar <https://pmglab.top/kggsee/download/lib/v1/kggsee.jar>`_
      - The KGGSEE program
      - 46 MB
    * - `resources/ <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EpXRqLXIToZItErUHiDNDO0Bk-jeiAtIlA-abGjOCdbqEw?e=3Jhy5g>`_
      - A OneDrive folder containing all running resource files provided by us
      - 
    * - `resources.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EYhQXE95WZFMqERo_xNOhZUB4mE73oh8Gs6ObS9aJe3XmA>`_
      - Running resource files except for reference genotypes and eQTL summary statistics 
      - 362 MB
    * - `tutorials.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EWqZHY25tT5Nq1GMwtl06ocBHoTAXGyBTH74zAp68dv5VA?e=tPtZ7B>`_
      - A tutorial dataset to run through :ref:`the four types of analyses <four_analyses>`
      - 155 MB


Set up an environment for the Quick tutorials 
=============================================

A quick and easy way to set up an environment for the :ref:`Quick tutorials <quick_tutorials>` is

#. Download `kggsee.jar <https://pmglab.top/kggsee/download/lib/v1/kggsee.jar>`_, `resources.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EYhQXE95WZFMqERo_xNOhZUB4mE73oh8Gs6ObS9aJe3XmA>`_ and `tutorials.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EWqZHY25tT5Nq1GMwtl06ocBHoTAXGyBTH74zAp68dv5VA?e=tPtZ7B>`_
#. Unzip ``resources.zip`` and ``tutorials.zip``
#. Put ``kggsee.jar``, ``resources/`` and ``tutorials/`` under one directory.

where `resources.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EYhQXE95WZFMqERo_xNOhZUB4mE73oh8Gs6ObS9aJe3XmA>`_ contains

.. list-table::
    :widths: 1 1
    :header-rows: 1
    :class: tight-table

    * - File
      - Description
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
    * - ``resources/GTEx_v8_TMM_all.gene.meanSE.txt.gz``
      - The gene-level expression profile of the `GTEx v8 <https://www.gtexportal.org/home/>`_ tissues
    * - ``resources/GTEx_v8_TMM_all.transcript.meanSE.txt.gz``
      - The transcript-level expression profile of the `GTEx v8 <https://www.gtexportal.org/home/>`_ tissues


and `tutorials.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EWqZHY25tT5Nq1GMwtl06ocBHoTAXGyBTH74zAp68dv5VA?e=tPtZ7B>`_ contains

.. list-table::
    :widths: 1 1
    :header-rows: 1
    :class: tight-table

    * - File
      - Description
    * - ``tutorials/scz_gwas_eur_chr1.tsv.gz``
      - Chromosome 1 summary statistics of a schizophrenia GWAS with a European sample.
    * - ``tutorials/1kg_hg19_eur_chr1.vcf.gz``
      - Chromosome 1 genotypes of the European panel of the 1000 Genomes Project
    * - ``tutorials/GTEx_v8_gene_BrainBA9.eqtl.txt.gz``
      - eQTL summary statistics calculated from the brain BA9 gene-level expression profile of GTEx v8
    * - ``tutorials/GTEx_v8_transcript_BrainBA9.eqtl.txt.gz``
      - eQTL summary statistics calculated from the brain BA9 transcript-level expression profile of GTEx v8


Set up an environment for customized analyses 
=============================================

In addition to the files packaged in `resources.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EYhQXE95WZFMqERo_xNOhZUB4mE73oh8Gs6ObS9aJe3XmA>`_, reference genotypes of five 1000 Genomes Project super populations and eQTL summary statistics of 49 GTEx v8 tissues are also available for downloading under `resources/ <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EpXRqLXIToZItErUHiDNDO0Bk-jeiAtIlA-abGjOCdbqEw?e=3Jhy5g>`_:

.. list-table::
    :widths: 1 2
    :header-rows: 1
    :class: tight-table

    * - File
      - Description
    * - `resources/hg19/gty/*.vcf.gz <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/Etg8dblAlUtGhtyN9RO49e0BvkXzgZj6Byy7PtNOUdLMMA?e=TftaGO>`_
      - VCF files of each super-population panel of the `1000 Genomes Project <https://www.internationalgenome.org/>`_ using hg19 coordinates. Each VCF file includes biallelic variants with MAF>0.01 of the super population. The VCF files include autosomes and chrX.
    * - `resources/hg38/gty/*.vcf.gz <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/Ep3EPaJSEqtAk_Eh7I7X4OwB9MDNe-LEwGUTFGC1V__O-A?e=sJyI59>`_
      - VCF files of each super-population panel of the `1000 Genomes Project <https://www.internationalgenome.org/>`_ using hg38 coordinates. Each VCF file includes biallelic variants with MAF>0.01 of the super population.  The VCF files include only autosomes.
    * - `resources/hg19/eqtl/*.eqtl.txt.gz <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EnhWhqLUNcpOrh6O3enFvCUBRvQ13v2970tcpOnNmmlKyg?e=JhXZh1>`_
      - cis-eQTL summary statistics using hg19 coordinates calculated from the gene or transcript-level expression profile of the `GTEx v8 <https://www.gtexportal.org/home/>`_ dataset
    * - `resources/hg38/eqtl/*.eqtl.txt.gz <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EtWxtqj5HTRHsEw4IiZ9xAMBu9S8Defi67pmL3_rNUjb9w?e=oCg45g>`_
      - cis-eQTL summary statistics using hg38 coordinates calculated from the gene or transcript-level expression profile of the `GTEx v8 <https://www.gtexportal.org/home/>`_ dataset


Then, a straightforward way to set up an environment for customized analyses is

#. Download `kggsee.jar <https://pmglab.top/kggsee/download/lib/v1/kggsee.jar>`_ and `resources.zip <https://mailsysueducn-my.sharepoint.com/:u:/g/personal/limiaoxin_mail_sysu_edu_cn/EYhQXE95WZFMqERo_xNOhZUB4mE73oh8Gs6ObS9aJe3XmA>`_
#. Unzip ``resources.zip``, and put ``kggsee.jar`` and ``resources/`` under one directory
#. Download the reference genotypes (`1kg_hg19 <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/Etg8dblAlUtGhtyN9RO49e0BvkXzgZj6Byy7PtNOUdLMMA?e=ks1hm1>`_ or `1kg_hg38 <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/Ep3EPaJSEqtAk_Eh7I7X4OwB9MDNe-LEwGUTFGC1V__O-A?e=d3KbyH>`_) of the population that matches your GWAS.
#. For running EMIC or eDESE, also download the eQTL summary statistics (`eqtl_hg19 <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EnhWhqLUNcpOrh6O3enFvCUBRvQ13v2970tcpOnNmmlKyg?e=1jkl06>`_ or `eqtl_hg38 <https://mailsysueducn-my.sharepoint.com/:f:/g/personal/limiaoxin_mail_sysu_edu_cn/EtWxtqj5HTRHsEw4IiZ9xAMBu9S8Defi67pmL3_rNUjb9w?e=ufFapJ>`_) of phenotype-associated tissues.
#. To prepare customized resource files, refer to :ref:`Detailed Document <detailed_document>` for descriptions of the file formats.

