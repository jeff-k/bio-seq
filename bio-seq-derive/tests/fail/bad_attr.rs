use bio_seq::prelude::*;

#[derive(Codec, PartialEq, Debug, Hash, Eq, Copy, Clone)]
enum BadAttr {
    A = 0b000,
    #[display('x')]
    #[alt(15)]
    #[asdf]
    C = 0b001,
    G = 0b010,
    T = 0b100,
    X = 0b011,
}

fn main() {
    let _b = BadAttr::C;
}
