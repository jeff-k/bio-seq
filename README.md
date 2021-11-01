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

* `bio_seq::alphabet::Dna`: DNA use the lexicographically ordered 2-bit representation

* `bio_seq::alphabet::Iupac`: IUPAC  nucleotide ambiguity codes are represented with 4 bits

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

* `bio_seq::alphabet::Amino`: Amino acid sequences are represented with 6 bits.

TODO: deal with alternate (e.g. mamalian mitochondrial) codes

## Kmers

Kmers are sequences of DNA with a fixed size. These are implemented with const generics.

## Minimisers for free

The 2-bit representation of DNA sequences is lexicographically ordered:

```rust
// find the lexicographically minimum 8-mer
fn minimise(seq: Seq<Dna>) -> Option<Kmer::<8>> {
    seq.kmers::<8>().min()
}
```

## Conversion with `From` and `Into`

`Iupac` from `Dna`; `Seq<Iupac>` from `Seq<Dna>`

`Amino` from `Kmer<3>`; `Seq<Amino>` from `Seq<Dna>`
  * Sequence length not a multiple of 3 is an error

`Seq<Iupac>` from `Amino`; `Seq<Iupac>` from `Seq<Amino>`

`Vec<Seq<Dna>>` from `Seq<Iupac>`: A sequence of IUPAC codes can generate a list of DNA sequences of the same length.


## Drop-in compatibility with `rust-bio`

meant to replace Text/TextSlice

## TODO

* benchmarking
* clever bit twiddling hacks for stuff like converting from `u8` to 2-bit representation
* macros for defining alphabet bit strings more concisely
