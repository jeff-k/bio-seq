<div class="title-block" style="text-align: center;" align="center">

# `bio-seq`

### Bit packed and well-typed biological sequences
</div>

```rust
use bio_seq::*;

let seq = dna!("ACTGCTAGCA");

for kmer in seq.kmers::<8>() {
	println!("{}", kmer);
}
```

* `bio_seq::Dna`: DNA use the lexicographically ordered 2-bit representation

* `bio_seq::Iupac`: IUPAC  nucleotide ambiguity codes are represented with 4 bits

	```
	  A C G T
	  -------
	S 0 1 1 0
	- 0 0 0 0
	C 0 1 0 0
	N 1 1 1 1
	B 0 1 1 1
	  ... etc.
	```
	This supports membership resolution with bitwise operations:

	```rust
    assert_eq!(
        format!("{}", iupac!("AS-GYTNA") | iupac!("ANTGCAT-")),
        "ANTGYWNA"
    );
    assert_eq!(
        format!("{}", iupac!("ACGTSWKM") & iupac!("WKMSTNNA")),
        "A----WKA"
    );
	```
The Iupac struct implements `From<Dna>`

* `bio_seq::Amino`: Amino acid sequences are represented with 6 bits.

   The representation of amino acids is designed to be easy to coerce from sequences of 2-bit encoded DNA.
   TODO: deal with alternate (e.g. mamalian mitochondrial) translation codes

## Kmers

Kmers are sequences with a fixed size. These are implemented with const generics.

`K * Codec::WIDTH` must fit in a `usize` (i.e. 64). For larger Kmers use `bigk::Kmer`: (TODO)

### Minimisers for free

The 2-bit representation of DNA sequences is lexicographically ordered:

```rust
// find the lexicographically minimum 8-mer
fn minimise(seq: Seq<Dna>) -> Option<Kmer::<8>> {
    seq.kmers::<8>().min()
}
```

## Derived codecs

Alphabet coding/decoding is derived from the variant names and discriminants of enum types:

```rust
#[derive(Clone, Copy, Debug, PartialEq, Codec)]
#[width = 2]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}
```

The `width` attribute specifies how many bits the encoding requires per symbol.

## Little endian

Kmers are represented stored as `usize`s with the least significant bit first.

```rust
dna!("C") == 0b01 // not 0b0100_0000
dna!("CT") == 0b11_01
```

## Conversion with `From` and `Into`

`Iupac` from `Dna`; `Seq<Iupac>` from `Seq<Dna>`

`Amino` from `Kmer<3>`; `Seq<Amino>` from `Seq<Dna>` (TODO)
  * Sequence length not a multiple of 3 is an error

`Seq<Iupac>` from `Amino`; `Seq<Iupac>` from `Seq<Amino>` (TODO)

`Vec<Seq<Dna>>` from `Seq<Iupac>`: A sequence of IUPAC codes can generate a list of DNA sequences of the same length. (TODO)

### Deref coercion

TODO: find out if `Kmer<Dna, K>` -> `Kmer<Amino, K/3>` is possible

## Drop-in compatibility with `rust-bio`

meant to replace Text/TextSlice

## TODO

* benchmarking
* macros for defining alphabet codecs more concisely
* wider SIMD-sized Kmers
