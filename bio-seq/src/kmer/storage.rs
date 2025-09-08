const fn make_2bit_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut i: usize = 0;
    while i < 256 {
        let b0: u8 = (i as u8 & 0b11_00_00_00) >> 6;
        let b1: u8 = (i as u8 & 0b00_11_00_00) >> 2;
        let b2: u8 = (i as u8 & 0b00_00_11_00) << 2;
        let b3: u8 = (i as u8 & 0b00_00_00_11) << 6;

        table[i] = b3 | b2 | b1 | b0;
        i += 1;
    }
    table
}

pub(crate) const REV_2BIT: [u8; 256] = make_2bit_table();

pub(crate) mod sealed {
    use crate::Bs;

    pub trait KmerStorage: Copy + Clone + PartialEq + std::fmt::Debug {
        const BITS: usize;
        type BaN: AsRef<Bs> + AsMut<Bs>;

        fn to_bitarray(self) -> Self::BaN;
        fn from_bitslice(bs: &Bs) -> Self;

        //        fn rotate_left(self, n: u32) -> Self;
        //        fn rotate_right(self, n: u32) -> Self;

        fn shiftr(&mut self, n: u32);

        fn shiftl(&mut self, n: u32);

        fn mask(&mut self, bits: usize);

        fn complement(&mut self, mask: usize);
        fn rev_blocks_2(&mut self);
    }
}

pub trait KmerStorage: sealed::KmerStorage {}

impl KmerStorage for usize {}

impl KmerStorage for u64 {}

impl KmerStorage for u128 {}
