//! This example demonstrates using bio-seq with [noodles](https://crates.io/crates/noodles) to count amino acid k-mers.
//!
//! It encodes amino acids sequences in an efficient bit-packed representation that can be directly used to index an array of counts.
//!
//! Run with:
//! ```sh
//! cargo run --example aminokmers -- path/to/proteins.faa
//! ```

use std::fs::File;
use std::io::{self, BufReader};
use std::marker::PhantomData;
use std::path::PathBuf;

use bio_seq::prelude::*;
use noodles::fasta::Reader;

use clap::Parser;

#[derive(Parser)]
#[command(
    name = "aminokmers",
    about = "Example of using noodles with bio-seq to count amino acid kmers"
)]
struct Cli {
    #[arg(help = "Path to amino acid fasta file")]
    fasta: PathBuf,
}

/// Construct a histogram on the heap that's indexed directly with Kmers.
struct Histogram<C: Codec, const K: usize> {
    counts: Vec<usize>,
    _p: PhantomData<C>,
}

impl<C: Codec, const K: usize> Histogram<C, K> {
    fn new() -> Self {
        Self {
            counts: vec![0; 1 << (K * C::BITS as usize)],
            _p: PhantomData,
        }
    }

    fn add(&mut self, kmer: Kmer<C, K, usize>) {
        self.counts[usize::from(&kmer)] += 1;
    }

    /*
        fn get(&self, kmer: Kmer<C, K, usize>) -> usize {
            self.counts[usize::from(&kmer)]
        }
    */

    fn items(&self) -> impl Iterator<Item = (Kmer<C, K, usize>, usize)> + '_ {
        self.counts
            .iter()
            .enumerate()
            .filter(|(_, &count)| count > 0)
            .map(|(i, &count)| (Kmer::<C, K>::from(i), count))
    }
}

fn main() -> io::Result<()> {
    let args = Cli::parse();
    let fasta = File::open(&args.fasta)?;

    let mut reader = Reader::new(BufReader::new(fasta));
    let mut histogram = Histogram::<Amino, 3>::new();

    for result in reader.records() {
        let record = result?;
        let seq: Seq<Amino> = record
            .sequence()
            .as_ref()
            .try_into()
            .expect("could not parse fasta");

        for kmer in seq.kmers() {
            histogram.add(kmer);
        }
    }

    for (kmer, count) in histogram.items() {
        println!("{kmer}: {count}");
    }

    Ok(())
}
