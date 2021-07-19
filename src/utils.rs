
/// Gets the complement of a nucleotide.
/// Returns None if the input char is not a valid nucleotide.
pub fn complement(nucleo: char) -> Option<char>{
  match nucleo{
    'A' => Some('T'),
    'C' => Some('G'),
    'G' => Some('C'),
    'T' => Some('A'),
    _ => None
  }
}

/// Gets the reverse complement of a unitig.
/// Returns None if the input contains not a valid nucleotide.
pub fn rev_compl(seq: &str) -> Option<String>{
  seq.chars().map(complement).rev().collect()
}

/// Gets the normalized version of a unitig.
/// Returns None if the input contains not a valid nucleotide.
pub fn norm(seq: &str) -> Option<String>{
  Some(String::from(seq)).min(rev_compl(seq))
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test() {
    assert_eq!("ATGC", norm("GCAT").unwrap());
  }
}