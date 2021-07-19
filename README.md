# Closed unitigs

Generates the closed unitigs starting from a [BCALM](https://github.com/GATB/bcalm) generated FASTA file representing a de Bruijn graph of sequencing data.

[What are closed unitigs?](CLOSED_UNITIGS.md)


## Usage
```sh
closed-unitigs INPUT
```

### Workflow example

1. Generate the de Bruijn graph with BCALM:
```sh
bcalm -in list.fa -kmer-size 16 -max-memory 1000 -all-abundance-counts 
```

2. Generate the closed unitigs:
```sh
closed-unitigs list.unitigs.fa
```

## Download builds
  * [Linux (64 bit)](https://gitlab.com/DPDmancul/closed-unitigs/-/jobs/artifacts/main/raw/target/x86_64-unknown-linux-gnu/release/closed-unitigs?job=linux-gnu-64)
  * [Linux (32 bit)](https://gitlab.com/DPDmancul/closed-unitigs/-/jobs/artifacts/main/raw/target/i686-unknown-linux-gnu/release/closed-unitigs?job=linux-gnu-32)
  * [macOS (64 bit)](https://gitlab.com/DPDmancul/closed-unitigs/-/jobs/artifacts/main/raw/target/x86_64-apple-darwin/release/closed-unitigs?job=macos-64)
  * [Windows (64 bit)](https://gitlab.com/DPDmancul/closed-unitigs/-/jobs/artifacts/main/raw/target/x86_64-pc-windows-gnu/release/closed-unitigs.exe?job=windows-mingw-64)
  <!-- * [Windows (32 bit)](https://gitlab.com/DPDmancul/closed-unitigs/-/jobs/artifacts/main/raw/target/i686-pc-windows-gnu/release/closed-unitigs.exe?job=windows-mingw-32) -->
  <!-- * [Linux (armv7)](https://gitlab.com/DPDmancul/closed-unitigs/-/jobs/artifacts/main/raw/target/armv7-unknown-linux-gnueabihf/release/closed-unitigs?job=linux-arm) -->
