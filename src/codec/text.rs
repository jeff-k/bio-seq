use crate::seq::Seq;
use core::marker::PhantomData;
use core::ops::Deref;

use crate::codec::{Codec, Complement};
use crate::ParseBioError;

pub struct Text<A: Codec = Dna>(Vec<u8>, PhantomData<A>);

impl<A: Codec> Deref for Text<A> {
    type Target = Seq<A>;

    fn deref(&self) -> &Self::Target {
        //        &Seq::<A>::from(self.0.clone())
        //unsafe { std::mem::transmute(&self.0) }
        unimplemented!()
    }
}

type Dna = u8;

impl Codec for Dna {
    type Error = ParseBioError;

    const WIDTH: u8 = 8;

    fn unsafe_from_bits(b: u8) -> Self {
        b
    }

    fn try_from_bits(b: u8) -> Result<Self, Self::Error> {
        Ok(b)
    }

    fn from_char(_c: char) -> Result<Self, Self::Error> {
        unimplemented!()
    }

    fn to_char(self) -> char {
        self.into()
    }
}

impl Complement for Dna {
    fn comp(self) -> Self {
        unimplemented!()
    }
}

/*
pub struct Amino;

impl Codec for Amino {

}
*/
