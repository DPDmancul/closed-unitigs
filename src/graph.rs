//! Represents a de Bruijn graph

#[path="./unitig.rs"]
mod unitig;

use snafu::Snafu;
use std::io::{BufRead, Write};
use regex::Regex;
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use unitig::*;

#[derive(Debug, Snafu)]
/// Describes and error on graph generation
enum GraphError {
  #[snafu(display("Unknown '{}' nucleotide", nucleo))]
  WrongNucleotide{nucleo: char},
  #[snafu(display("Unknown nucleotide into sequence \"{}\"", seq))]
  WrongNucleotideInto{seq: String}
}

impl From<UnitigError> for GraphError {
  /// Convert UnitigError in GraphError
  fn from(e: UnitigError) -> GraphError {
    match e {
      UnitigError::WrongNucleotide{nucleo} => GraphError::WrongNucleotide{nucleo}
    }
  }
}

/// Represents a graph node
#[derive(Debug)]
struct Edge{
  to: usize,
  start: bool,
  end: bool
}

/// Represents a graph node
#[derive(Debug)]
struct Node{
  /// Sequence of nucleotides
  kmer: Unitig,
  /// Sequence reverse complement
  complement: Unitig,
  /// kmer count
  count: u32,
  /// paths to other nodes
  out: Vec<Edge>,
  /// paths from other nodes
  into: Vec<Edge>
}

impl Node{
  /// Creates a new node containing the given sequence and vector of counts.
  fn new(kmer: Unitig, count: u32) -> Node {
    Node{
      complement: kmer.rev_compl(),
      kmer,
      count,
      out: vec![],
      into: vec![]
    }
  }
}

/// Represents a de Bruijn graph
pub struct Graph {
  /// List of nodes of the graph
  nodes: Vec<Node>,
  /// size of the k-mers
  k: usize
}

impl Graph {
  /// Builds an empty graph
  fn new(k: usize) -> Graph {
    Graph{
      nodes: Vec::new(),
      k
    }
  }

  /// Appends a new node to the graph
  fn append(&mut self, seq: String, count: u32) -> Result<(), GraphError>{
    self.nodes.push(Node::new(seq.try_into()?, count));
    Ok(())
  }

  /// Finds support of u, with memorization
  fn supp(u: &Unitig, k: usize, supp: &mut HashMap<Unitig, u32>) -> u32 {
    if let Some(&s) = supp.get(u) {
      // Use memorization
      return s;
    }
    // Compite support taking the minimum of k-mer counts
    let s = *(0..u.len()-k+1).into_iter()
      .map(|i| Unitig::try_from(&u[i..i+k]).unwrap()) // Safe because coming from an unitig
      .map(|u| supp.get(&u).unwrap_or(&0)) // k-mers counts must be already memorized; if the k-mer is not present its support is zero
      .min().unwrap_or(&0);
    supp.insert(u.clone(), s); // Memorize
    s
  }

  /// Finds closure of m
  fn closure<'a>(&'a self, m: &Unitig, first: (&'a Node, bool), last: (&'a Node, bool), k: usize, supp: &mut HashMap<Unitig, u32>, (is_closed, n_closed): (&mut HashMap<Unitig, bool>, &mut u32)) -> Unitig {
    let (mut m, mut first, mut last) = (m.clone(), first, last); // Make those mutable

    // Explore the graph trying to extend this unitig until support decreases
    'clo: loop {
      // dbg!(&m);
      let my_supp = Self::supp(&m, k, supp);

      // Try to extend to the right
      for Edge{to, start, end} in &last.0.out {
        if *start != last.1 {continue} // direction do not match
        let node = &self.nodes[*to]; // target node
        if m.contains(&node.kmer) || m.contains(&node.complement) {continue} // avoid loops
        let kmer = if *end {&node.kmer} else {&node.complement};
        // dbg!("out", node, to, start, end);
        let c = node.count;
        if c >= my_supp {
          if c == my_supp {
            // The closed unitig we are building is valid also for this k-mer
            is_closed.insert(kmer.clone(), true);
            *n_closed += 1;
          }
          m = &m + kmer; // Join
          last = (node, *end); // Extend
          continue 'clo
        } // elsewhere the support decreases and so we cannot extend
      }
      // Try to extend to the left
      for Edge{to, start, end} in &first.0.into {
        if *start != first.1 {continue} // direction do not match
        let node = &self.nodes[*to]; // target node
        if m.contains(&node.kmer) || m.contains(&node.complement) {continue} // avoid loops
        let kmer = if *end {&node.kmer} else {&node.complement};
        // dbg!("into", node, to, start, end);
        let c = node.count;
        if c >= my_supp {
          if c == my_supp {
            // The closed unitig we are building is valid also for this k-mer
            is_closed.insert(kmer.clone(), true);
            *n_closed += 1;
          }
          m = kmer + &m; // Join
          first = (node, *end); // Extend
          continue 'clo
        } // elsewhere the support decreases and so we cannot extend
      }
      break
    };
    is_closed.insert(m.clone(), true);
    *n_closed += 1;
    m //clo
  }

  /// Shrinks a closed unitig removing head and tail with higher support
  fn shrink(u: Unitig, k: usize, supp: &HashMap<Unitig, u32>) -> (Unitig, u32) {
    let (mut a, mut b) = (0, u.len()); // extremities
    let my_supp = supp.get(&u).unwrap();
    // Try shrink on left
    while a+k < b && supp.get(&u[a..a+k].try_into().unwrap()).unwrap() > my_supp { a += 1 }
    // Try shrink on right
    while b >= k && supp.get(&u[b-k..b].try_into().unwrap()).unwrap() > my_supp { b -= 1 }
    // Return shrunk closed unitig
    (u[a..b].try_into().unwrap(), *my_supp)
  }

  /// Finds closed unitigs
  pub fn close<T: Write, U: Write>(&self, fasta: &mut T, counts: &mut U) {
    let k = self.k;
    let mut closed = HashMap::<Unitig, u32>::new(); // using a map instead of a vector avoids duplicates

    {
      let mut supp = HashMap::<Unitig, u32>::new();
      let mut is_closed = HashMap::<Unitig, bool>::new();
      // Compute k-mers supports (equal to their counts)
      for node in &self.nodes {
        supp.insert(node.kmer.clone(), node.count);
        is_closed.insert(node.kmer.clone(), false);
      }

      let mut n_closed = 0;
      // Close and shrink all k-mers
      for node in &self.nodes {
        if is_closed[&node.kmer] {continue}
        print!("Closing {:?} ({:.2}%)\r", node.kmer, (1. + n_closed as f64)/self.nodes.len() as f64*100.);
        let close = self.closure(&node.kmer, (&node, true), (&node, true), k, &mut supp, (&mut is_closed, &mut n_closed));
        let (u, c) = Self::shrink(close, k, &supp);
        closed.insert(u, c);
      }
    }

    let mut closed: Vec<_> = closed.iter().collect();
    closed.sort_by_key(|(_, &c)| c); // Sort by count to reduce count differences
    for (u, c) in closed {
      writeln!(fasta, ">\n{}", u).unwrap();
      writeln!(counts, "{}", c).unwrap();
    }
  }

}

impl<T: BufRead> std::convert::From<T> for Graph {
  /// Build a de Bruijn graph from FASTA file
  fn from(buf: T) -> Graph {
    let mut graph = Graph::new(0); // temporary k = 0
    let mut nodes = Vec::<(usize, usize)>::new(); // left, right
    let mut edges = Vec::<((usize, bool), (usize, bool))>::new(); // from, to which are (node, dir)

    let count_re = Regex::new(r"ab:Z:(\d+(?: \d+)*)").unwrap();
    let link_re = Regex::new(r"L:([+-]):(\d+):([+-])").unwrap();

    let mut opt = String::new();

    for (index, line) in buf.lines().enumerate() {
      let line = line.unwrap();

      print!("Reading fasta file (line {})\r", index+1);

      // If line is even get options
      if index%2 == 0{
        if !opt.is_empty() || !line.starts_with('>') {
          panic!("Syntax error at line {}: \"{}\"", index+1, line);
        }
        opt = line;
        continue;
      }

      // Get counts
      let count: Vec<_> = match count_re.captures(&opt) {
          Some(r) => r.get(1).unwrap().as_str(),
          None => ""
        }.split(' ').map(|s| s.parse().unwrap()).collect();

      // Get k
      if graph.k == 0 {
        graph.k = line.len() - count.len() + 1; // line.len = count.len + k - 1
        println!("\x1B[2K\rk = {}", graph.k);
      }

      // Append this node
      let first = graph.nodes.len();
      for (i, &c) in count.iter().enumerate() {
        assert!(i + graph.k <= line.len());
        if let Err(e) = graph.append(String::from(&line[i..i+graph.k]), c) {
          panic!("{}; on line {}", e, index+1);
        }
      }
      let last = graph.nodes.len()-1;
      nodes.push((first, last));

      // Edges between k-mers of the unitig
      for i in first..last {
        graph.nodes[i].out.push(Edge{to: i+1, start: true, end: true});
        graph.nodes[i].into.push(Edge{to: i+1, start: false, end: false}); // Reverse complement
      }
      // Reverse direction
      for i in first+1..=last {
        graph.nodes[i].out.push(Edge{to: i-1, start: false, end: false}); // Reverse complement
        graph.nodes[i].into.push(Edge{to: i-1, start: true, end: true});
      }

      // Get edges between unitigs
      for group in link_re.captures_iter(&opt) {
        edges.push((
          (nodes.len()-1, group[1].starts_with('+')),
          (group[2].parse().unwrap(), group[3].starts_with('+'))
        ));
      }

      opt = String::new();
    }

    // Store edges
    for ((from, start), (to, end)) in edges {
      let (from, to) = (nodes[from], nodes[to]);

      // check direction
      let from = if start {from.1} else {from.0};
      let to = if end {to.0} else {to.1};

      if from == to {continue} // avoid self loops

      graph.nodes[from].out.push(Edge{to, start, end});
      graph.nodes[to].into.push(Edge{to: from, start: end, end: start}); // Reverse direction
    }

    // dbg!(graph.nodes);

    graph
  }
}
