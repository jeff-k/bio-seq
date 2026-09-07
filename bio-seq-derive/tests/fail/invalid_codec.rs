use bio_seq::prelude::*;

#[derive(Codec)]
struct NotEnum {
    a: (),
    c: (),
    g: (),
    t: (),
}

#[derive(Codec)]
enum NoDiscriminants {
    A,
    C,
    G,
    T,
}

#[derive(Codec, PartialEq, Debug, Hash)]
#[bits(9)]
enum ZeroWidth {}

#[derive(Codec, PartialEq, Debug, Hash)]
#[bits("nine")]
enum WeirdWidth {}

#[derive(Codec, PartialEq, Debug, Hash)]
#[bits(2)]
enum BadWidth {
    A = 0b000,
    C = 0b001,
    G = 0b010,
    T = 0b100,
    X = 0b011,
}

#[derive(Codec, PartialEq, Debug, Hash, Eq, Copy, Clone)]
#[bits(3)]
enum BadDiscriminant {
    A = 0b000,
    C = 3,
    G = 0b010,
    T = 0x04,
    X = 1.0,
}

// Test whether codec is #[repr(u8)]
/*
#[derive(Codec, Eq, Debug, Hash, PartialEq, Copy, Clone)]
enum NoRepr {
    A = 0b11,
    C = 0b01,
    G = 0b00,
    T = 0b10,
}
*/

#[derive(Codec, PartialEq, Debug, Hash, Eq, Copy, Clone)]
#[repr(u8)]
enum ActuallyGood {
    A = 0b001,
    C = b'\0',
    G = 0b0_10,
    T = 3,
    #[alt(16)]
    X = 4,
}

fn main() {
    let _c = NotEnum {
        a: (),
        c: (),
        g: (),
        t: (),
    };

    let _x = ActuallyGood::C;

    assert_eq!(ActuallyGood::BITS, 5);
}
