<div class="title-block" style="text-align: center;" align="center">

# bio-seq

### Bit-packed and well-typed biological sequences
</div>

```rust
use bio_seq::prelude::*;

let seq = dna!("ATACGATCGATCGATCGATCCGT");

// iterate over the 8-mers of the reverse complement
for kmer in seq.revcomp().kmers::<8>() {
    println!("{}", kmer);
}
```

The IUPAC nucleotide ambiguity codes naturally encode a set of bases for each position:

```rust
use bio_seq::prelude::*;

let seq = iupac!("AGCTNNCAGTCGACGTATGTA");
let pattern = Seq::<Iupac>::from_str("AYG").unwrap();

for slice in seq.windows(pattern.len()) {
    if pattern.contains(slice) {
        println!("{} matches pattern", slice);
    }
}
```

The primary design goal of this crate is to make translating between biological sequence types safe and convenient:

```rust
// debruijn sequence of all 3-mers:
let seq: Seq<Dna> =
    dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
let aminos: Seq<Amino> = Seq::from_iter(seq.kmers().map(|kmer| kmer.into()));
assert_eq!(
    aminos,
    amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
);
```

## Contents

* [Codec](#codecs): Coding/Decoding schemes for the characters of a biological sequence
* [Seq](#sequences): A sequence of encoded characters
* [Kmer](#kmers): A fixed size sequence of length `K`
* [Derivable codecs](#derivable-codecs): This crate offers utilities for defining your own bit-level encodings
* [Safe conversion](#sequence-conversion) between sequences

## Codecs

The `Codec` trait describes the coding/decoding process for the characters of a biological sequence. This trait can be derived procedurally. There are three built-in codecs:

### `codec::Dna`
Using the lexicographically ordered 2-bit representation

### `codec::Iupac`
IUPAC  nucleotide ambiguity codes are represented with 4 bits. This supports membership resolution with bitwise operations. Logical `or` is the union:

```rust
assert_eq!(iupac!("AS-GYTNA") | iupac!("ANTGCAT-"), iupac!("ANTGYWNA"));
```

Logical `and` is the intersection of two iupac sequences:

```rust
assert_eq!(iupac!("ACGTSWKM") & iupac!("WKMSTNNA"), iupac!("A----WKA"));
```

### `codec::Amino`
Amino acid sequences are represented with 6 bits. The representation of amino acids is designed to be easy to coerce from sequences of 2-bit encoded DNA.

## Sequences

Strings of encoded characters are packed into `Seq`s. Slicing, chunking, and windowing return `SeqSlice`s. `Seq<A: Codec>`/`&SeqSlice<A: Codec>` are analogous to `String`/`&str`.

All data is stored little-endian. This effects the order that sequences map to the integers ("colexicographic" order):

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

## Kmers

kmers are sequences with a fixed size that can fit into a register. these are implemented with const generics.

### Dense encodings

For dense encodings, a lookup table can be populated and indexed in constant time with the `usize` representation:

```rust
fn kmer_histogram<C: Codec, const K: usize>(seq: &SeqSlice<C>) -> Vec<usize> {
    // For dna::Dna our histogram will need 4^4
    // bins to count every possible 4-mer.
    let mut histo = vec![0; 1 << (C::WIDTH * K as u8)];

    for kmer in seq.kmers::<K>() {
        histo[usize::from(kmer)] += 1;
    }

    histo
}
```

This example builds a histogram of kmer occurences

## Sketching

### Hashing

The `Hash` trait is implemented for Kmers

### Canonical Kmers

Depending on the application, it may be permissible to superimpose the forward and reverse complements of a kmer:

```rust
k = kmer!("ACGTGACGT");
let canonical = k ^ k.revcomp(); // TODO: implement ReverseComplement for Kmer
```

### Kmer minimisers

The 2-bit representation of nucleotides is ordered `A < C < G < T`. Sequences and kmers are stored in little-endian and are ordered "colexicographically". This means that `AAAA < CAAA < GAAA < ... < AAAC < ... < TTTT`


```rust
fn minimise(seq: Seq<Dna>) -> Option<Kmer::<Dna, 8>> {
    seq.kmers().min()
}
```

### Example: Hashing minimiser of canonical Kmers

```rust
for ckmer in seq.window(8).map(|kmer| hash(kmer ^ kmer.revcomp())) {
    // TODO: example
    ...
}
```

## Derivable codecs

Sequence coding/decoding is derived from the variant names and discriminants of enum types:

```rust
use bio_seq_derive::Codec;
use bio_seq::codec::Codec;

#[derive(Clone, Copy, Debug, PartialEq, Codec)]
#[width = 2]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl From<Dna> for u8 {
    fn from(dna: Dna) -> Self {
        dna as u8
    }
}
```

The `width` attribute specifies how many bits the encoding requires per symbol. The maximum supported is 8. If this attribute isn't specified then the optimal width will be chosen.

## Sequence conversions

`Iupac` from `Dna`; `Seq<Iupac>` from `Seq<Dna>`

`Amino` from `Kmer<3>`; `Seq<Amino>` from `Seq<Dna>`
  * Sequence length not a multiple of 3 is an error

`Seq<Iupac>` from `Amino`; `Seq<Iupac>` from `Seq<Amino>` (TODO)

`Vec<Seq<Dna>>` from `Seq<Iupac>`: A sequence of IUPAC codes can generate a list of DNA sequences of the same length. (TODO)
