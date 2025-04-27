SQLRanges Documentation
=======================

SQLRanges is a Python library for performing **SQL-backed operations** on genomic intervals using **SQLite** or **DuckDB**. It allows seamless conversion of results to **PyRanges** for downstream analysis.

With SQLRanges, you can easily interact with genomic data stored in SQL databases and perform operations such as exon counting, merging exon intervals, and finding overlapping genes.

This documentation will guide you through the installation, usage, and API of SQLRanges.

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

To install **SQLRanges**, simply run the following command:

.. code-block:: bash

   pip install SQLRanges

Usage
=====

Here's an example of how to use the `SQLRanges` class:

.. code-block:: python

   from SQLRanges import SQLRanges

   # Initialize the SQLRanges object
   sqlranges = SQLRanges(input_file="Homo_sapiens.GRCh38.112.chr.gtf", 
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

   SQLRanges

Modules
=======

The SQLRanges package includes several modules that handle database operations, genomic interval manipulations, and PyRanges conversions.

.. toctree::
   :maxdepth: 1
   :caption: Modules:

   modules
