-- Create reference database for lipid network

-- Show databases and tables
.databases
.tables

-- open an existing database
.open lipidref.db

-- Drop database
drop database lipidref;

-- Add master list table
create table master(
  Source TEXT,
  Species TEXT,
  Organ TEXT,
  Matrix TEXT,
  Tx TEXT,
  standard_name TEXT,
  Name TEXT,
  Major Class TEXT,
  Class TEXT,
  Abbreviation TEXT,
  Saturation TEXT,
  Highest_Saturation_Degree TEXT,
  Cluster TEXT,
  Subclass TEXT,
  Chain_Length TEXT,
  Chain_Type TEXT,
  Subclass_2 TEXT,
  synthesis_pathway TEXT,
  label TEXT,
  label_saturation TEXT
);

-- Add lipid network nodes
create table lipid_nodes(
  Node TEXT,
  label TEXT,
  x_man TEXT,
  y_man TEXT,
  ID TEXT,
  ID_string TEXT,
  abbreviation TEXT,
  major_class TEXT,
  synthesis_pathway TEXT,
  PMIDs TEXT,
  URLs TEXT
);

-- Add lipid network edges
create table lipid_edges(
  'from' TEXT,
  'to' TEXT,
  synthesis_pathway TEXT,
  leading_name TEXT,
  enzyme_id TEXT,
  ensembl TEXT,
  hgnc_symbol TEXT,
  PMCIDs_DOI TEXT,
  URLs TEXT
);

-- import .txt file containing annotations into respective tables
.mode csv
.separator ","
.import --csv --skip 1 0_masterlist.csv master

-- nodes
.import --csv --skip 1 0_lipidnet_base_nodes.csv lipid_nodes

-- edges
.import --csv --skip 1 0_lipidnet_base_edges.csv lipid_edges

-- add rows to an existing table

-- Remove table
drop table tablename;

-- Export the results of a query (run line-by-line)
-- .show shows changes to the database parameters
.show
.headers on
.mode csv
.output test.csv
select * from lipid_nodes;
.output stdout

-- quit SQLite
.quit