<div class="title-block" style="text-align: center;" align="center">

# `bio-seq`

### Bit packed and well-typed biological sequences
</div>

```rust
use bioseq::*;

let seq = dna!("ACTGCTAGCA");

for kmer in seq.kmers::<8>() {
	println!("{}", kmer);
}
```

* `bioseq::alphabet::Dna`: DNA use the lexicographically ordered 2-bit representation

* `bioseq::alphabet::Iupac`: IUPAC  nucleotide ambiguity codes are represented with 4 bits

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
	let a = iupac!("ASGYTNA");
	let b = iupac!("ANTGCAT");

	assert_eq!(a.union(b), iupac!(_));
	assert_eq!(a.intersection(b), iupac!(_));
	```

* *TODO* `bioseq::alphabet::amino`: Amino acid sequences

## Kmers

Kmers are sequences of DNA with a fixed size. These are implemented with const generics.

## Minimisers for free

```rust
// find the lexicographically minimum 8-mer
fn minimise(seq: Seq) -> Kmer::<8> {
    seq.kmers::<8>().min()
}
```

## Drop-in compatibility with `rust-bio`

meant to replace Text/TextSlice

## TODO

* benchmarking
* clever bit twiddling hacks for stuff like converting from `u8` to 2-bit representation
