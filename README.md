<div class="title-block" style="text-align: center;" align="center">

# `bio-seq`

### Bit-packed and well-typed biological sequences
</div>

```rust
use bio_seq::{dna, Seq, FromStr};
use bio_seq::codec::{dna::Dna, ReverseComplement};

let seq = dna!("ATACGATCGATCGATCGATCCGT");

// iterate over the 8-mers of the reverse complement
for kmer in seq.revcomp().kmers::<8>() {
    println!("{}", kmer);
}
```

The IUPAC nucleotide ambiguity codes naturally encode a set of bases for each position:

```rust
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
// debruijn sequence of order 3
let seq: Seq<Dna> =
    dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
let aminos: Seq<Amino> = Seq::from_vec(seq.kmers().map(|kmer| kmer.into()).collect());
assert_eq!(
    aminos,
    amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
);
```

## Contents

* [Codec](#codecs): Encoding scheme for the 'characters' of a biological sequence
* [Seq](#sequences): A sequence of encoded characters
* [Kmer](#kmers): A fixed size sequence of length `k`
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

Strings of encoded biological characters are packed into `Seq`s. Slicing, chunking, and windowing return `SeqSlice`s. `Seq<A: Codec>`/`&SeqSlice<A: Codec>` are analogous to `String`/`&str`.

## Kmers

kmers are sequences with a fixed size that can fit into a register. these are implemented with const generics.

### Dense encodings

For dense encodings, a lookup table can be populated and indexed in constant time with the `usize` representation:

```rust
let mut histogram = vec![0; 1 << C::WIDTH * K];
```

### Hashing

The `Hash` trait is implemented for Kmers

### Canonical Kmers

Depending on the application, it may be permissible to superimpose the forward and reverse complements of a kmer:

```rust
k = kmer!("ACGTGACGT");
let canonical = k ^ k.revcomp(); // TODO: implement ReverseComplement for Kmer
```

### Kmer minimisers

The 2-bit representation of DNA sequences is lexicographically ordered:

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
use bio_seq::codec::{Codec, ParseBioErr};

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

The `width` attribute specifies how many bits the encoding requires per symbol. The maximum supported is 8.

Kmers are stored as `usize`s with the least significant bit first.

## Sequence conversions

`Iupac` from `Dna`; `Seq<Iupac>` from `Seq<Dna>`

`Amino` from `Kmer<3>`; `Seq<Amino>` from `Seq<Dna>`
  * Sequence length not a multiple of 3 is an error

`Seq<Iupac>` from `Amino`; `Seq<Iupac>` from `Seq<Amino>` (TODO)

`Vec<Seq<Dna>>` from `Seq<Iupac>`: A sequence of IUPAC codes can generate a list of DNA sequences of the same length. (TODO)
