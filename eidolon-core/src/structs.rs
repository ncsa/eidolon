//! Common structs are objects that we will build off to create the models and simulation. It
//! includes a struct for holding DNA information read in from a fasta file, objects for holding
//! data on sequencing errors and variants, and a transition matrix object that forms the basis for
//! most of the models in NEAT. Storing these in a common library will help us create other utilities
//! that will generate the data that generate-reads can then use to run simulations.
pub mod bed_record;
pub mod distributions;
pub mod haplotype_map;
pub mod mutated_map;
pub mod nucleotides;
pub mod read_record;
pub mod sequence_block;
pub mod sv_model;
pub mod transition_matrix;
pub mod variants;
