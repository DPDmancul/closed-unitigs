mod graph;

use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter};

fn main() {

  let input_file = &env::args().collect::<Vec<String>>()[1];
  let output_fasta = input_file.clone() + ".clo.fa";
  let output_counts = input_file.clone() + ".clo.counts";

  // Read BCALM FASTA file and generate graph
  let graph = graph::Graph::from(BufReader::new(File::open(input_file).unwrap()));
  // Close unitigs and write output files
  graph.close(&mut BufWriter::new(File::create(&output_fasta).unwrap()), &mut BufWriter::new(File::create(&output_counts).unwrap()));

}

