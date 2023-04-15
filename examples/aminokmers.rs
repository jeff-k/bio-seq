extern crate bio_seq;
extern crate bio_streams;

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;

use bio_streams::fasta::Fasta;

use bio_seq::prelude::*;

#[derive(Parser)]
struct Cli {
    faa: PathBuf,
}

const K: usize = 4;

fn main() {
    let args = Cli::parse();

    let faa: Fasta<BufReader<File>, Seq<Amino>> =
        Fasta::new(BufReader::new(File::open(&args.faa).unwrap()));

    let mut count: usize = 0;
    let mut total: usize = 0;
    let mut histo = Box::new([0; 1 << (K * Amino::WIDTH as usize)]);

    for contig in faa {
        for kmer in contig.unwrap().seq.kmers::<K>() {
            histo[usize::from(kmer)] += 1;
            total += 1;
        }
        count += 1;
    }

    let n = 100;
    let mut empty_counts = 0;
    for i in 0..histo.len() {
        if histo[i] > n {
            println!("{}\t{:?}", Kmer::<Amino, K>::from(i), histo[i]);
        }
        if histo[i] == 0 {
            empty_counts += 1;
        }
    }
    println!(
        "Read {} contigs. Showing counts of {}-mers that occur more than {} times ({} total). {}/{} of possible {}-mers have a count of 0.",
        count, K, n, total, empty_counts, 1 << (K * Amino::WIDTH as usize), K,
    );
}
