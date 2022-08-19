<div class="title-block" style="text-align: center;" align="center">

# `bio-seq`

### Bit-packed and well-typed biological sequences
</div>

```rust
use bio_seq::{Seq, FromStr};
use bio_seq::codec::{dna::Dna, ReverseComplement};
 
fn main() {
    // declare a sequence of 2-bit encoded DNA
    let seq = Seq::<Dna>::from_str("TACGATCGATCGATCGATC").unwrap();

    // iterate over the 8-mers of the reverse complement of our sequence
    for kmer in seq.revcomp().kmers::<8>() {
        println!("{}", kmer);
    }
}
```

## Contents

* [Codec](#codecs): Encoding scheme for the 'characters' of a biological sequence
* [Seq](#sequences): A sequence of encoded characters
* [Kmer](#kmers): A fixed size sequence of length `k`
* [Derivable codecs](#derivable-codecs): This crate offers utilities for defining your own bit-level encodings
* [Safe conversion](#sequence-conversion) between sequences

## Codecs

The `Codec` trait desribes the coding/decoding process for the characters of a biological sequence. There are four built-in codecs:

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

### TODO `codec::ascii::Dna`
for the 8-bit ascii representation of IUPAC ambiguity codes. This is intended to be compatible with existing bioinformatics packages such as `rust-bio`.

## Sequences

Strings of encoded biological characters are packed into `Seq`s. This are allocated on the heap and may be mutable. Slicing, chunking, and windowing return instances of `SeqSlice`.

## Kmers

kmers are sequences with a fixed size that can fit into a register. these are implemented with const generics.

`k * Codec::width` must fit in a `usize` (i.e. 64). for larger kmers use `bigk::kmer`: TODO

### Dense encodings

For dense encodings, a lookup table can be populated and indexed in constant time with the `usize` representation:

TODO: finish example
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

`Amino` from `Kmer<3>`; `Seq<Amino>` from `Seq<Dna>` (TODO)
  * Sequence length not a multiple of 3 is an error

`Seq<Iupac>` from `Amino`; `Seq<Iupac>` from `Seq<Amino>` (TODO)

`Vec<Seq<Dna>>` from `Seq<Iupac>`: A sequence of IUPAC codes can generate a list of DNA sequences of the same length. (TODO)

TODO: deal with alternate (e.g. mamalian mitochondrial) translation codes
