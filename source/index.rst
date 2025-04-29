sqlranges Documentation
=======================

sqlranges is a Python library for performing **SQL-backed operations** on genomic intervals using **SQLite** or **DuckDB**. It allows seamless conversion of results to **PyRanges** for downstream analysis.

With sqlranges, you can easily interact with genomic data stored in SQL databases and perform operations such as exon counting, merging exon intervals, and finding overlapping genes.

This documentation will guide you through the installation, usage, and API of sqlranges.

Table of Contents
=================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   api
   modules

Installation
============

To install **sqlranges**, simply run the following command:

.. code-block:: bash

   pip install sqlranges

Usage
=====

Here's an example of how to use the `sqlranges` class:

.. code-block:: python

   from sqlranges import sqlranges

   # Initialize the sqlranges object
   sqlranges = sqlranges(input_file="Homo_sapiens.GRCh38.112.chr.gtf", 
                         table_name="human", 
                         db_name="db-human.duckdb", 
                         backend="duckdb")

   # Convert to PyRanges
   pyranges_obj = sqlranges.to_pyranges()

   # Count exons
   exon_counts = sqlranges.count_exons()

API Documentation
================

The following sections describe the available methods and their usage.

.. toctree::
   :maxdepth: 2
   :caption: API Reference:

   sqlranges

Modules
=======

The sqlranges package includes several modules that handle database operations, genomic interval manipulations, and PyRanges conversions.

.. toctree::
   :maxdepth: 1
   :caption: Modules:

   modules
