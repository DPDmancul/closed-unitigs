
#[path="./utils.rs"]
mod utils;

use std::{
  ops::{Deref, Add},
  fmt::{self, Display},
  hash::{Hash, Hasher},
  cmp::Ordering,
  convert::{From, TryFrom, TryInto}
};
use snafu::Snafu;

#[derive(Debug, Snafu)]
/// Describes and error on unitig generation
pub enum UnitigError {
  #[snafu(display("Unknown '{}' nucleotide", nucleo))]
  WrongNucleotide{nucleo: char},
}

#[derive(Clone, Default, Debug)]
#[repr(transparent)]
/// Represents an unitig
pub struct Unitig(String);

impl Unitig {
  /// Returns the reverse complement of this unitig
  pub fn rev_compl(&self) -> Unitig {
    Unitig(utils::rev_compl(&self.0).unwrap())
  }

  /// Returns the normalized unitig (the lexicographically lower among itself and its reverse complement)
  pub fn norm(&self) -> Unitig {
    Unitig(utils::norm(&self.0).unwrap())
  }

  /// Check if this unitig contains as substring the given unitig
  pub fn contains(&self, x: &Unitig) -> bool {
    self.0.contains(&x.0)
  }
}

impl Add for &Unitig {
  type Output = Unitig;

  /// Concatenates two unitigs sharing a common tail-head
  fn add(self, other: Self) -> Self::Output {
      let common = self.0.len().min(other.0.len())-1;
      assert!(self.0[self.0.len()-common..] == other.0[..common], "The two Unitigs {:?} and {:?} are not joinable", self, other);
      Unitig(String::from(&self.0) + &other.0[common..])
  }
}

impl Display for Unitig {
  /// Displays an unitig
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "{}", self.0)
  }
}

// Casting

impl Deref for Unitig {
  type Target = str;

  /// Deref an unitig to String
  fn deref(&self) -> &Self::Target {
      &self.0
  }
}

impl From<Unitig> for String {
  /// Generates a String from a unitig
  fn from(u: Unitig) -> String {
    u.0
  }
}

impl TryFrom<&str> for Unitig {
  type Error = UnitigError;

  /// Generates an Unitig from a string slice
  fn try_from(u: &str) -> Result<Unitig, Self::Error> {
    String::from(u).try_into()
  }
}

impl TryFrom<String> for Unitig {
  type Error = UnitigError;

  /// Generates an Unitig from a String
  fn try_from(u: String) -> Result<Unitig, Self::Error> {
    let u = u.to_uppercase();
    for c in u.chars(){
      match c {
        'A' | 'C' | 'G' | 'T' => (),
        nucleo => return Err(UnitigError::WrongNucleotide{nucleo})
      }
    }
    Ok(Unitig(u))
  }
}

// Comparison

impl Ord for Unitig {
  /// Lexicographically compare two unitigs by normal form
  fn cmp(&self, other: &Self) -> Ordering {
    self.norm().0.cmp(&other.norm().0)
  }
}

impl PartialOrd for Unitig {
  /// Lexicographically compare two unitigs by normal form
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(&other))
  }
}

impl PartialEq for Unitig {
  /// Lexicographically compare two unitigs by normal form
  fn eq(&self, other: &Self) -> bool {
    self.norm().0 == other.norm().0
  }
}

impl Eq for Unitig {}

impl Hash for Unitig {
  /// Hash the normal form of this unitig
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.norm().0.hash(state);
  }
}

