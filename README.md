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

TODO: deal with alternate (e.g. mamalian mitochondrial) codes

## Kmers

Kmers are sequences of DNA with a fixed size. These are implemented with const generics.

K is meant to fit in a `usize`. For larger Kmers or amino acid sequences use SeqSlices:

```rust


```

## Minimisers for free

The 2-bit representation of DNA sequences is lexicographically ordered:

```rust
// find the lexicographically minimum 8-mer
fn minimise(seq: Seq<Dna>) -> Option<Kmer::<8>> {
    seq.kmers::<8>().min()
}
```

## little end first

Kmers can be casted to `usize` indexes

```rust
dna!"C" == 0b01 // not 0b0100_0000
dna!"CT" == 0b0111
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
* macros for defining alphabet codecs more concisely
* wider SIMD-sized Kmers
