pub enum Dna {
    A,
    C,
    G,
    T,
}

pub enum Iupac {
    A,
    C,
    G,
    T,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    N,
}

pub enum Amino {
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
}

pub trait Complementary {
    fn complement(base: Self) -> Self;
}

impl Iupac {
    fn from_char(b: char) -> Self {
        match b {
            'N' => Iupac::N,
            _ => panic!(),
        }
    }
}

impl Dna {
    fn from_u8(b: u8) {
        //      let c = b & 0x3; // mask the last 2 bits
        //  let d[0] = (b ^ 0x2) ^ (b 0x1);
    }
}
