[![Docs.rs](https://docs.rs/bio-seq/badge.svg)](https://docs.rs/bio-seq)
[![CI status](https://github.com/jeff-k/bio-seq/actions/workflows/rust.yml/badge.svg)](https://github.com/jeff-k/bio-seq/actions/workflows/rust.yml)

<div class="title-block" style="text-align: center;" align="center">

# bio-seq

### Bit-packed and well-typed biological sequences
</div>

This crate provides types and traits for sequences of genomic data. Common encodings are provided and can be extended with the `Codec` trait.

Short sequences of fixed length (kmers) are given special attention.

### Quick start

Add [bio-seq](https://crates.io/crates/bio-seq) to `Cargo.toml`:

```toml
[dependencies]
bio-seq = "0.13"
```

Iterating over the [kmer](https://docs.rs/bio-seq/latest/bio_seq/kmer)s of a [sequence](https://docs.rs/bio-seq/latest/bio_seq/seq):

```rust
use bio_seq::prelude::*;

// The `dna!` macro packs a static sequence with 2-bits per symbol at compile time
let seq = dna!("ATACGATCGATCGATCGATCCGT");

// iterate over the 8-mers of the reverse complement
for kmer in seq.revcomp().kmers::<8>() {
    println!("{kmer}");
}

// ACGGATCG
// CGGATCGA
// GGATCGAT
// GATCGATC
// ATCGATCG
// ...
```

Sequences are analogous to rust's string types and follow similar dereferencing conventions:

```rust
// Static sequences behave like static string literals:
let s: &'static str = "hello!";
let seq: &'static SeqSlice<Dna> = dna!("CGCTAGCTACGATCGCAT");

// Sequences can also be copied as `Kmer`s:
let kmer: Kmer<Dna, 18> = dna!("CGCTAGCTACGATCGCAT").into();
// or with the kmer! macro:
let kmer = kmer!("CGCTAGCTACGATCGCAT");

// `Seq`s are allocated on the heap like `String`s are:
let s: String = "hello!".into();
let seq: Seq<Dna> = dna!("CGCTAGCTACGATCGCAT").into();

// Alternatively, a `Seq` can be fallibly encoded at runtime:
let seq: Seq<Dna> = "CGCTAGCTACGATCGCAT".try_into().unwrap();

// `&SeqSlice`s are analogous to `&str`, `String` slices:
let slice: &str = &s[1..3];
let seqslice: &SeqSlice<Dna> = &seq[2..4];
```

Sequences can be read from popular third-party crates like [noodles](https://crates.io/crates/noodles):

```rust
let mut reader = noodles::fasta::Reader::new(BufReader::new(fasta));

for result in reader.records() {
    let record = result?;
    let seq: Seq<Dna> = record
        .sequence()
        .as_ref()
        .try_into()
        .unwrap()

    // ...
}
```

## Application examples

* [Saving packed sequences to binary files](https://github.com/jeff-k/bio-seq/blob/main/bio-seq/examples/seq2bin.rs)
* [Using noodles and counting kmers](https://github.com/jeff-k/bio-seq/blob/main/bio-seq/examples/aminokmers.rs)
* [Codec benchmarks](https://github.com/jeff-k/bio-seq/blob/main/bio-seq/examples/codec-bench.rs) comparing memory use and entropy of lossy/degenerate encodings

## Philosophy

Many bioinformatics crates implement their own kmer packing logic. This effort began as a way to define types and traits that allow kmer code to be shared between projects. It quickly became apparent that a kmer type doesn't make sense without being tightly coupled to a general type for sequences. The scope of this crate will be limited to operating on fixed and arbitrary length sequences with an emphasis on safety.

Some people like to engineer clever bit twiddling hacks to reverse complement a sequence and some people want to rapidly prototype succinct datastructures. Most people don't want to worry about endianess. The strength of rust is that we can safely abstract the science from the engineering to work towards both objectives.

## Contributing

Contributions and suggestions are very much welcome. Check out the [Roadmap](https://github.com/jeff-k/bio-seq/issues/10) to get started.

## Contents

* [Codec](#codecs): Coding/Decoding schemes for the symbols of a biological sequence
* [Seq](#sequences): A sequence of encoded symbols
* [Kmer](#kmers): A fixed size sequence of length `K`
* [Derivable codecs](#defining-new-codecs): This crate offers utilities for defining your own bit-level encodings

## [Sequences](https://docs.rs/bio-seq/latest/bio_seq/seq)

Strings of encoded symbols are packed into [`Seq`](https://docs.rs/bio-seq/latest/bio_seq/seq/struct.Seq.html). Slicing, chunking, and windowing return [`SeqSlice`](https://docs.rs/bio-seq/latest/bio_seq/seq/struct.SeqSlice.html). `Seq<A: Codec>` and `&SeqSlice<A: Codec>` are analogous to `String` and `&str`. As with the standard string types, these are stored on the heap and implement `Clone`.

## [Kmers](https://docs.rs/bio-seq/latest/bio_seq/kmer)

kmers are short sequences of length `k` that generally fit into a register (e.g. `usize`, or SIMD vector) and implement `Copy`. `k` is a compile-time constant.

All data is stored little-endian. This effects the order that sequences map to the integers:

```rust
for i in 0..=15 {
    println!("{}: {}", i, Kmer::<Dna, 5>::from(i));
}
```

```
0: AAAAA
1: CAAAA
2: GAAAA
3: TAAAA
4: ACAAA
5: CCAAA
6: GCAAA
7: TCAAA
8: AGAAA
9: CGAAA
10: GGAAA
11: TGAAA
12: ATAAA
13: CTAAA
14: GTAAA
15: TTAAA
```

### Succinct encodings

A lookup table can be indexed in constant time by treating kmers directly as `usize`:

```rust
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

    fn add(&mut self, kmer: Kmer<C, K>) {
        self.counts[usize::from(&kmer)] += 1;
    }

    fn get(&self, kmer: Kmer<C, K>) -> usize {
        self.counts[usize::from(&kmer)]
    }
}
```


## Sketching

### Minimising

The [2-bit representation](https://docs.rs/bio-seq/latest/bio_seq/codec/dna) of nucleotides is ordered `A < C < G < T`. Sequences and kmers are stored little-endian and are ordered "colexicographically". This means that `AAAA` < `CAAA` < `GAAA` < `...` < `AAAC` < `...` < `TTTT`:

```rust
let seq = dna!("GCTCGATCGTAAAAAATCGTATT");
let minimiser = seq.kmers::<8>().min().unwrap();

assert_eq!(minimiser, Kmer::from(dna!("GTAAAAAA")));
```

### Hashing

`Hash` is implemented for sequence and kmer types so equal values of these types will hash identically:

```rust
let seq_arr: &'static SeqSlice<Dna> = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCT");
let seq: Seq<Dna> = seq_arr.into();
let seq_slice: &SeqSlice<Dna> = &seq;
let kmer: Kmer<Dna, 32> = seq_arr.into();

assert_eq!(hash(seq_arr), hash(&seq));
assert_eq!(hash(&seq), hash(&seq_slice));
assert_eq!(hash(&seq_slice), hash(&kmer));
```

### Hashing minimisers

In practice we want to hash sequences that we minimise:

```rust
fn hash<T: Hash>(seq: T) -> u64 {
    let mut hasher = DefaultHasher::new();
    seq.hash(&mut hasher);
    hasher.finish()
}

let (minimiser, min_hash) = seq
    .kmers::<16>()
    .map(|kmer| (kmer, hash(&kmer)))
    .min_by_key(|&(_, hash)| hash)
    .unwrap();
```

### Canonical kmers:

To consider both the forward and reverse complement of kmers when minimising:

```rust
let (canonical_minimiser, canonical_hash) = seq
    .kmers::<16>()
    .map(|kmer| {
        let canonical_hash = hash(min(kmer, kmer.revcomp()));
        (kmer, canonical_hash)
    })
    .min_by_key(|&(_, hash)| hash)
    .unwrap();
```

Although it's more efficient to minimise `seq` and `seq.revcomp()` separately.

## Codecs

The `Codec` trait describes the coding/decoding process for the symbols of a biological sequence. This trait can be derived procedurally. There are four built-in codecs:

* `codec::Dna`, Using the lexicographically ordered 2-bit representation

* `codec::Iupac`, IUPAC nucleotide ambiguity codes are represented with 4 bits. This automatically gives us membership semantics for bitwise operations. Logical `or` is the union:

    ```rust
    assert_eq!(iupac!("AS-GYTNA") | iupac!("ANTGCAT-"), iupac!("ANTGYWNA"));
    ```

    Logical `and` is the intersection of two iupac sequences:

    ```rust
    assert_eq!(iupac!("ACGTSWKM") & iupac!("WKMSTNNA"), iupac!("A----WKA"));
    ```

* `codec::Text`, utf-8 strings that are read directly from common plain-text file formats can be treated as sequences. Additional logic can be defined to ensure that `'a' == 'A'` and for handling `'N'`.

* `codec::Amino`, Amino acid sequences are represented with 6 bits. The representation of amino acids is designed to be easy to coerce from sequences of 2-bit encoded DNA.

## Defining new codecs

Custom codecs can be defined by implementing the `Codec` trait.

In simple cases the `Codec` trait can be derived from the variant names and discriminants of enum types:

```rust
use bio_seq_derive::Codec;
use bio_seq::codec::Codec;

#[derive(Clone, Copy, Debug, PartialEq, Codec)]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}
```

Note that you need to explicitly provide a "discriminant" (e.g. `0b00`) in the enum.

A `#[width(n)]` attribute specifies how many bits the encoding requires per symbol. The maximum supported is 8. If this attribute isn't specified then the optimal width will be chosen.

`#[alt(...,)]` and `#[display('x')]` attributes can be used to define alternative representations or display the item with a special character. Here is the definition for the stop codon in `codec::Amino`:

```rust
pub enum Amino {
    #[display('*')] // print the stop codon as a '*'
    #[alt(0b001011, 0b100011)] // TGA, TAG
    X = 0b000011, // TAA (stop)
```

## Sequence conversions

### Translation table traits 

[Translation tables](https://docs.rs/bio-seq/latest/bio_seq/translation) provide methods for translating codons into amino acids.

Enable the translation feature in `Cargo.toml`:

```
[dependencies]
bio-seq = { version="0.13", features=["translation"] }
```

```rust
pub trait TranslationTable<A: Codec, B: Codec> {
    fn to_amino(&self, codon: &SeqSlice<A>) -> B;
    fn to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError>;
}

/// A partial translation table where not all triples of characters map to amino acids
pub trait PartialTranslationTable<A: Codec, B: Codec> {
    fn try_to_amino(&self, codon: &SeqSlice<A>) -> Result<B, TranslationError>;
    fn try_to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError>;
}
```

The standard genetic code is provided as a `translation::STANDARD` constant:

```rust
use crate::prelude::*;
use crate::translation::STANDARD;
use crate::translation::TranslationTable;

let seq = dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");

let aminos: Seq<Amino> = seq
    .windows(3)
    .map(|codon| STANDARD.to_amino(&codon))
    .collect::<Seq<Amino>>();

assert_eq!(
    aminos,
    Seq<Amino>::try_from("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK").unwrap()
);
```

### Custom translation tables

Instantiate a translation table from a type that implements `Into<HashMap<Seq<A>, B>>`:

```rust
let codon_mapping: [(Seq<Dna>, Amino); 6] = [
    (dna!("AAA"), Amino::A),
    (dna!("ATG"), Amino::A),
    (dna!("CCC"), Amino::C),
    (dna!("GGG"), Amino::E),
    (dna!("TTT"), Amino::D),
    (dna!("TTA"), Amino::F),
];

let table = CodonTable::from_map(codon_mapping);

let seq: Seq<Dna> = dna!("AAACCCGGGTTTTTATTAATG");
let mut amino_seq: Seq<Amino> = Seq::new();

for codon in seq.chunks(3) {
    amino_seq.push(table.try_to_amino(codon).unwrap());
}
```

Implementing the `TranslationTable` trait directly:

```rust
struct Mitochondria;

impl TranslationTable<Dna, Amino> for Mitochondria {
    fn to_amino(&self, codon: &SeqSlice<Dna>) -> Amino {
        if *codon == dna!("AGA") {
            Amino::X
        } else if *codon == dna!("AGG") {
            Amino::X
        } else if *codon == dna!("ATA") {
            Amino::M
        } else if *codon == dna!("TGA") {
            Amino::W
        } else {
            Amino::unsafe_from_bits(Into::<u8>::into(codon))
        }
    }

    fn to_codon(&self, _amino: Amino) -> Result<Seq<Dna>, TranslationError> {
        unimplemented!()
    }
}
```


