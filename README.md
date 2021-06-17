<div class="title-block" style="text-align: center;" align="center">

# bio-seq <!-- omit in toc -->

## Biological sequences backed by bit vectors <!-- omit in toc -->
</div>

```rust
let seq = dna!("ACTGCTAGCA");
let y: Seq<Iupac> = Iupac::from_vec([A, N, A, T, N, N, G]);

assert_eq!(seq.to_usize(), 0b00011011);
```

IUPAC nucleotide ambiguity codes are represented with 4 bits:

```
  A C G T
  ------
S 0 1 1 0
- 0 0 0 0
C 0 1 0 0
N 1 1 1 1
B 0 1 1 1
  ... etc.
```

```rust
let seq = dna!("ACGCTACGA");

let seq = Seq<Dna>::from("ACTGCAT");

assert_eq!(seq.rc(), dna!"TCGTAGCGT");
```

```rust
let mut s: Seq<Dna> = Seq::new();
let y: Seq<Iupac> = iupac![A, N, A, T, N, N, G];
```

## Kmers

Kmers are slices of Seqs.

## Minimisers for free

```rust
// find the lexicographically minimum 8-mer
fn minimise(seq: Seq) -> Kmer<8> {
    seq.iter().min()
}
```

## Drop-in compatibility with `rust-bio`

replaces TextSlice

## TODO

* benchmarking
* clever bit twiddling hacks for stuff like converting from `u8` to 2-bit representation
