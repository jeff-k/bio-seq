impl SeqStorage for Vec<u8> {
    type Slice<'a> = &'a [u8];
    type Array<'a> = todo!();

    fn new() -> Self {
        todo!()
    }

    fn len(&self) -> usize {
        todo!()
    }

    fn with_capacity(_len: usize) -> Self {
        todo!()
    }

    fn is_empty(&self) -> bool {
        todo!()
    }

    fn as_slice(&self) -> &Self::Slice<'_> {
        todo!()
    }

    fn push(&mut self, _bits: u8) {
        todo!()
    }

   fn to_usize(&self) -> usize {
       todo!()
    }

    fn clear(&mut self) {
        todo!()
    }
}
