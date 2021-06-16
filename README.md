<div class="title-block" style="text-align: center;" align="center">

# bio-seq <!-- omit in toc -->

## Biological sequences backed by bit vectors <!-- omit in toc -->
</div>


```
A: 00, C: 01, G: 10, T: 11
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
let seq = dna!"ACGCTACGA";

assert_eq!(seq.rc(), dna!"TCGTAGCGT");
```

```rust
let mut s: Seq<Dna> = Seq::new();
let y: Seq<Iupac> = iupac![A, N, A, T, N, N, G];
```
