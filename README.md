# bio-seq
Types for biological sequences

```rust
let seq = dna!"ACGCTACGA";

assert_eq!(seq.rc(), dna!"TCGTAGCGT");
```
