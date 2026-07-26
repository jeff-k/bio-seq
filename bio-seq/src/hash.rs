// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Platform-stable hashing for bit-packed sequences.

use crate::Bs;
use bitvec::field::BitField;
use core::hash::Hasher;

/// Feed a bit slice into a hasher as a stable, platform-independent byte stream.
///
/// Bits are packed least-significant-first into bytes, matching the in-memory
/// `Lsb0` layout. Sidestepping `bitvec`'s `Hash` impl (which routes bits
/// through `bool::hash`, one hasher write per bit) and the `usize` storage
/// width keeps the byte stream identical on 32- and 64-bit targets.
#[inline]
pub(crate) fn hash_bits<H: Hasher>(bs: &Bs, state: &mut H) {
    // Extract whole bytes via `BitField::load_le::<u8>()`, which lets bitvec
    // do the bit-shuffle on the underlying storage word. Buffer them so the
    // hasher sees one bulk write per 64-byte block instead of one per byte.
    // `load_le` zero-pads the high bits of any partial trailing chunk.
    let mut buf = [0u8; 64];
    let mut len = 0;
    for chunk in bs.chunks(8) {
        buf[len] = chunk.load_le::<u8>();
        len += 1;
        if len == buf.len() {
            state.write(&buf);
            len = 0;
        }
    }
    if len > 0 {
        state.write(&buf[..len]);
    }
}
