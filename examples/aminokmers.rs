extern crate bio_seq;
extern crate bio_streams;

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;

use bio_streams::fasta::Fasta;

use bio_seq::codec::amino::Amino;
use bio_seq::{Kmer, Seq};

#[derive(Parser)]
struct Cli {
    faa: PathBuf,
}

fn main() {
    let args = Cli::parse();

    let faa: Fasta<BufReader<File>, Seq<Amino>> =
        Fasta::new(BufReader::new(File::open(&args.faa).unwrap()));

    let mut count: usize = 0;
    const K: usize = 3;
    let mut histo = [0; 1 << (K * 6)];

    for contig in faa {
        for kmer in contig.seq.kmers::<K>() {
            histo[usize::from(kmer)] += 1;
        }
        count += 1;
        println!("read contig {}", count);
    }

    for i in 0..histo.len() {
        if histo[i] > 0 {
            println!("{}\t{:?}", Kmer::<Amino, K>::from(i), histo[i]);
        }
    }
}
