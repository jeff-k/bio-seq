//! Codec processing benchmarks (work in progress)
//!
//! Run with:
//! ```sh
//! cargo run --example codec_bench -- path/to/input.fasta
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::PathBuf;

use bio_seq::codec;
use bio_seq::prelude::*;
use noodles::fasta::Reader;

use clap::{Parser, ValueEnum};

#[derive(Debug, Copy, Clone, ValueEnum)]
enum CodecType {
    Dna,
    Iupac,
    DegenerateWS,
    Amino,
    Text,
}

#[derive(Parser)]
#[command(name = "codec-bench", about = "Compare bio-seq codecs")]
struct Cli {
    #[arg(help = "Path to FASTA file")]
    fasta: PathBuf,

    #[arg(short, long, value_enum, help = "Encoding to use")]
    codec: CodecType,

    #[arg(short, long, default_value_t = 16, help = "Length of k-mers")]
    k: usize,
}

trait KmerCounter {
    fn add(&mut self, seq: &[u8]);
    fn get_stats(&self) -> (usize, f64);
}

/// Construct a histogram on the heap that's indexed directly with Kmers.
struct KmerCounts<C: Codec> {
    counts: HashMap<Seq<C>, usize>,
    k: usize,
    total: usize,
}

impl<C: Codec> KmerCounts<C> {
    fn new(k: usize) -> Self {
        Self {
            counts: HashMap::new(),
            k,
            total: 0,
        }
    }
}

impl<C: Codec> KmerCounter for KmerCounts<C> {
    fn add(&mut self, seq: &[u8]) {
        let seq: Seq<C> = seq.try_into().unwrap();
        for kmer in seq.windows(self.k) {
            *self.counts.entry(kmer.to_owned()).or_insert(0) += 1;
            self.total += 1;
        }
    }

    fn get_stats(&self) -> (usize, f64) {
        let total = self.total as f64;

        let entropy: f64 = self
            .counts
            .values()
            .map(|&count| {
                let p = count as f64 / total;
                -p * p.log2()
            })
            .sum();

        (self.total, entropy)
    }
}

fn main() -> io::Result<()> {
    let args = Cli::parse();
    let fasta = File::open(&args.fasta)?;
    let k = args.k;

    let mut reader = Reader::new(BufReader::new(fasta));

    let mut counts: Box<dyn KmerCounter> = match args.codec {
        CodecType::Dna => Box::new(KmerCounts::<Dna>::new(k)),
        CodecType::Iupac => Box::new(KmerCounts::<Iupac>::new(k)),
        CodecType::DegenerateWS => Box::new(KmerCounts::<codec::degenerate::WS>::new(k)),
        CodecType::Amino => Box::new(KmerCounts::<Amino>::new(k)),
        CodecType::Text => Box::new(KmerCounts::<codec::text::Dna>::new(k)),
    };

    for result in reader.records() {
        let record = result?;
        let seq: &[u8] = record.sequence().as_ref();

        counts.add(seq);
    }

    let (total, entropy) = counts.get_stats();
    // show stats
    println!("{k}, {total}, {entropy}");

    Ok(())
}
