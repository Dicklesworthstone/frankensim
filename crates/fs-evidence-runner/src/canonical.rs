//! Small bounded canonical-frame builder used by pure base-schema projections.

use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use fs_blake3::{ContentHash, hash_domain};

/// Canonical field sink shared by the count-only preflight and byte encoder.
///
/// Callers that need a complete-size preflight should express their field
/// sequence once and pass it to [`CanonicalFrameV1::preflighted`]. The closure
/// must be a pure description of the frame: it is invoked first against a
/// nonallocating [`CanonicalFrameSizeV1`] and then against the byte encoder.
pub(crate) trait CanonicalFrameSinkV1 {
    fn push_u8(&mut self, field: &'static str, value: u8) -> Result<(), ConstructionErrorV2>;

    fn push_u16(&mut self, field: &'static str, value: u16) -> Result<(), ConstructionErrorV2>;

    fn push_u32(&mut self, field: &'static str, value: u32) -> Result<(), ConstructionErrorV2>;

    fn push_bytes(&mut self, field: &'static str, value: &[u8]) -> Result<(), ConstructionErrorV2>;

    fn push_fixed_bytes_32(
        &mut self,
        field: &'static str,
        value: &[u8],
    ) -> Result<(), ConstructionErrorV2>;

    fn push_str(&mut self, field: &'static str, value: &str) -> Result<(), ConstructionErrorV2>;

    fn push_presence(
        &mut self,
        field: &'static str,
        present: bool,
    ) -> Result<(), ConstructionErrorV2>;
}

/// Count-only canonical frame sink.
///
/// This type retains only the checked encoded length and the inclusive frame
/// bound. It does not retain field values or allocate a caller-sized buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CanonicalFrameSizeV1 {
    encoded_len: usize,
    max_bytes: usize,
}

impl CanonicalFrameSizeV1 {
    pub(crate) fn new(magic: &'static [u8], max_bytes: usize) -> Result<Self, ConstructionErrorV2> {
        if magic.len() > max_bytes {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "canonical.magic",
                "magic fits the frame bound",
                magic.len(),
            ));
        }
        Ok(Self {
            encoded_len: magic.len(),
            max_bytes,
        })
    }

    /// Exact encoded length accumulated so far, including the magic.
    #[must_use]
    pub(crate) const fn exact_len(self) -> usize {
        self.encoded_len
    }

    fn extend(
        &mut self,
        field: &'static str,
        additional: usize,
    ) -> Result<(), ConstructionErrorV2> {
        self.encoded_len =
            checked_canonical_frame_length_v1(field, self.encoded_len, additional, self.max_bytes)?;
        Ok(())
    }
}

#[derive(Debug)]
pub(crate) struct CanonicalFrameV1 {
    bytes: Vec<u8>,
    max_bytes: usize,
}

impl CanonicalFrameV1 {
    pub(crate) fn new(magic: &'static [u8], max_bytes: usize) -> Result<Self, ConstructionErrorV2> {
        if magic.len() > max_bytes {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "canonical.magic",
                "magic fits the frame bound",
                magic.len(),
            ));
        }
        Ok(Self {
            bytes: magic.to_vec(),
            max_bytes,
        })
    }

    /// Preflight and encode one canonical frame from a single field sequence.
    ///
    /// The first closure invocation uses a count-only sink. Only after that
    /// pass succeeds is an exact-capacity byte buffer allocated. The second
    /// invocation is bounded by the preflighted length, and the final equality
    /// check refuses a shorter second pass.
    pub(crate) fn preflighted<F>(
        magic: &'static [u8],
        max_bytes: usize,
        encode_fields: F,
    ) -> Result<Self, ConstructionErrorV2>
    where
        F: Fn(&mut dyn CanonicalFrameSinkV1) -> Result<(), ConstructionErrorV2>,
    {
        let exact_len = Self::preflight_length(magic, max_bytes, &encode_fields)?;

        let mut bytes = Vec::with_capacity(exact_len);
        bytes.extend_from_slice(magic);
        let mut frame = Self {
            bytes,
            // A divergent larger second pass must refuse before Vec can grow
            // beyond the complete size established by the first pass.
            max_bytes: exact_len,
        };
        encode_fields(&mut frame)?;
        if frame.bytes.len() != exact_len {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "canonical.preflight",
                "encoded length exactly matches the size preflight",
                frame.bytes.len(),
            ));
        }
        Ok(frame)
    }

    /// Count one complete canonical frame without allocating its byte buffer.
    ///
    /// Composite constructors use this before any caller-cardinality-sized
    /// validation index. The later [`Self::preflighted`] call repeats the
    /// count from the same field encoder immediately before materialization,
    /// so semantic validation cannot make a stale length authoritative.
    pub(crate) fn preflight_length<F>(
        magic: &'static [u8],
        max_bytes: usize,
        encode_fields: &F,
    ) -> Result<usize, ConstructionErrorV2>
    where
        F: Fn(&mut dyn CanonicalFrameSinkV1) -> Result<(), ConstructionErrorV2> + ?Sized,
    {
        let mut size = CanonicalFrameSizeV1::new(magic, max_bytes)?;
        encode_fields(&mut size)?;
        Ok(size.exact_len())
    }

    pub(crate) fn push_u8(
        &mut self,
        field: &'static str,
        value: u8,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &[value])
    }

    pub(crate) fn push_i8(
        &mut self,
        field: &'static str,
        value: i8,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_u16(
        &mut self,
        field: &'static str,
        value: u16,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_u32(
        &mut self,
        field: &'static str,
        value: u32,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_u64(
        &mut self,
        field: &'static str,
        value: u64,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_i64(
        &mut self,
        field: &'static str,
        value: i64,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_u128(
        &mut self,
        field: &'static str,
        value: u128,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_i16(
        &mut self,
        field: &'static str,
        value: i16,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_i32(
        &mut self,
        field: &'static str,
        value: i32,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_i128(
        &mut self,
        field: &'static str,
        value: i128,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, &value.to_be_bytes())
    }

    pub(crate) fn push_bytes(
        &mut self,
        field: &'static str,
        value: &[u8],
    ) -> Result<(), ConstructionErrorV2> {
        let length = checked_u32_length(field, value.len())?;
        let additional = core::mem::size_of::<u32>()
            .checked_add(value.len())
            .ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    field,
                    "checked length prefix plus payload length",
                    value.len(),
                )
            })?;
        checked_canonical_frame_length_v1(field, self.bytes.len(), additional, self.max_bytes)?;

        // Capacity was checked for the complete field before either component
        // is appended, so a refusal cannot leave a partial length prefix.
        self.bytes.extend_from_slice(&length.to_be_bytes());
        self.bytes.extend_from_slice(value);
        Ok(())
    }

    pub(crate) fn push_fixed_bytes_32(
        &mut self,
        field: &'static str,
        value: &[u8],
    ) -> Result<(), ConstructionErrorV2> {
        if value.len() != 32 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                field,
                "exactly 32 unprefixed bytes",
                value.len(),
            ));
        }
        self.extend(field, value)
    }

    pub(crate) fn push_str(
        &mut self,
        field: &'static str,
        value: &str,
    ) -> Result<(), ConstructionErrorV2> {
        self.push_bytes(field, value.as_bytes())
    }

    pub(crate) fn push_presence(
        &mut self,
        field: &'static str,
        present: bool,
    ) -> Result<(), ConstructionErrorV2> {
        self.push_u8(field, u8::from(present))
    }

    pub(crate) fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Transfer the exactly preflighted canonical bytes without copying them.
    ///
    /// Component constructors call [`Self::root`] before this transfer when
    /// they also retain a content identity. Consuming the frame here preserves
    /// the one-allocation materialization contract of [`Self::preflighted`].
    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }

    pub(crate) fn root(&self, domain: &'static str) -> ContentHash {
        hash_domain(domain, &self.bytes)
    }

    fn extend(&mut self, field: &'static str, value: &[u8]) -> Result<(), ConstructionErrorV2> {
        checked_canonical_frame_length_v1(field, self.bytes.len(), value.len(), self.max_bytes)?;
        self.bytes.extend_from_slice(value);
        Ok(())
    }
}

impl CanonicalFrameSinkV1 for CanonicalFrameV1 {
    fn push_u8(&mut self, field: &'static str, value: u8) -> Result<(), ConstructionErrorV2> {
        Self::push_u8(self, field, value)
    }

    fn push_u16(&mut self, field: &'static str, value: u16) -> Result<(), ConstructionErrorV2> {
        Self::push_u16(self, field, value)
    }

    fn push_u32(&mut self, field: &'static str, value: u32) -> Result<(), ConstructionErrorV2> {
        Self::push_u32(self, field, value)
    }

    fn push_bytes(&mut self, field: &'static str, value: &[u8]) -> Result<(), ConstructionErrorV2> {
        Self::push_bytes(self, field, value)
    }

    fn push_fixed_bytes_32(
        &mut self,
        field: &'static str,
        value: &[u8],
    ) -> Result<(), ConstructionErrorV2> {
        Self::push_fixed_bytes_32(self, field, value)
    }

    fn push_str(&mut self, field: &'static str, value: &str) -> Result<(), ConstructionErrorV2> {
        Self::push_str(self, field, value)
    }

    fn push_presence(
        &mut self,
        field: &'static str,
        present: bool,
    ) -> Result<(), ConstructionErrorV2> {
        Self::push_presence(self, field, present)
    }
}

impl CanonicalFrameSinkV1 for CanonicalFrameSizeV1 {
    fn push_u8(&mut self, field: &'static str, _value: u8) -> Result<(), ConstructionErrorV2> {
        self.extend(field, core::mem::size_of::<u8>())
    }

    fn push_u16(&mut self, field: &'static str, _value: u16) -> Result<(), ConstructionErrorV2> {
        self.extend(field, core::mem::size_of::<u16>())
    }

    fn push_u32(&mut self, field: &'static str, _value: u32) -> Result<(), ConstructionErrorV2> {
        self.extend(field, core::mem::size_of::<u32>())
    }

    fn push_bytes(&mut self, field: &'static str, value: &[u8]) -> Result<(), ConstructionErrorV2> {
        let _length = checked_u32_length(field, value.len())?;
        let additional = core::mem::size_of::<u32>()
            .checked_add(value.len())
            .ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    field,
                    "checked length prefix plus payload length",
                    value.len(),
                )
            })?;
        self.extend(field, additional)
    }

    fn push_fixed_bytes_32(
        &mut self,
        field: &'static str,
        value: &[u8],
    ) -> Result<(), ConstructionErrorV2> {
        if value.len() != 32 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                field,
                "exactly 32 unprefixed bytes",
                value.len(),
            ));
        }
        self.extend(field, 32)
    }

    fn push_str(&mut self, field: &'static str, value: &str) -> Result<(), ConstructionErrorV2> {
        self.push_bytes(field, value.as_bytes())
    }

    fn push_presence(
        &mut self,
        field: &'static str,
        _present: bool,
    ) -> Result<(), ConstructionErrorV2> {
        self.extend(field, core::mem::size_of::<u8>())
    }
}

fn checked_u32_length(field: &'static str, length: usize) -> Result<u32, ConstructionErrorV2> {
    u32::try_from(length).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            field,
            "length representable as u32",
            length,
        )
    })
}

/// Checked encoded-length seam used by count-only and byte-emitting frames.
///
/// Keeping this callable within the crate lets bounded declaration tests
/// exercise arithmetic overflow without constructing an impossible allocation.
pub(crate) fn checked_canonical_frame_length_v1(
    field: &'static str,
    current: usize,
    additional: usize,
    maximum: usize,
) -> Result<usize, ConstructionErrorV2> {
    let next = current.checked_add(additional).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            field,
            "checked frame length",
            additional,
        )
    })?;
    if next > maximum {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            field,
            "value fits the canonical frame bound",
            next,
        ));
    }
    Ok(next)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::{
        CanonicalFrameSinkV1, CanonicalFrameSizeV1, CanonicalFrameV1,
        checked_canonical_frame_length_v1, checked_u32_length,
    };
    use crate::construction::ConstructionErrorKindV2;

    #[test]
    fn integer_and_presence_fields_have_independent_known_big_endian_bytes() {
        let mut frame = CanonicalFrameV1::new(b"\xAA\x55", 256).expect("fixture frame");
        frame.push_u8("u8", 0x81).expect("u8");
        frame.push_i8("i8", -2).expect("i8");
        frame.push_u16("u16", 0x0102).expect("u16");
        frame.push_i16("i16", -258).expect("i16");
        frame.push_u32("u32", 0x0304_0506).expect("u32");
        frame.push_i32("i32", -50_594_827).expect("i32");
        frame.push_u64("u64", 0x0708_090A_0B0C_0D0E).expect("u64");
        frame.push_i64("i64", -0x0102_0304_0506_0708).expect("i64");
        frame
            .push_u128("u128", 0x0F10_1112_1314_1516_1718_191A_1B1C_1D1E)
            .expect("u128");
        frame
            .push_i128("i128", -0x0102_0304_0506_0708_090A_0B0C_0D0E_0F10)
            .expect("i128");
        frame.push_presence("absent", false).expect("absence");
        frame.push_presence("present", true).expect("presence");

        let mut expected = vec![0xAA, 0x55, 0x81, 0xFE, 0x01, 0x02, 0xFE, 0xFE];
        expected.extend_from_slice(&[0x03, 0x04, 0x05, 0x06]);
        expected.extend_from_slice(&(-50_594_827_i32).to_be_bytes());
        expected.extend_from_slice(&[0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E]);
        expected.extend_from_slice(&(-0x0102_0304_0506_0708_i64).to_be_bytes());
        expected.extend_from_slice(&[
            0x0F, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x1B, 0x1C,
            0x1D, 0x1E,
        ]);
        expected
            .extend_from_slice(&(-0x0102_0304_0506_0708_090A_0B0C_0D0E_0F10_i128).to_be_bytes());
        expected.extend_from_slice(&[0, 1]);
        assert_eq!(frame.as_bytes(), expected);
    }

    #[test]
    fn byte_and_string_fields_use_exact_u32_length_prefixes() {
        let mut frame = CanonicalFrameV1::new(b"M", 32).expect("fixture frame");
        frame
            .push_bytes("bytes", &[0x00, 0x80, 0xFF])
            .expect("bytes");
        frame.push_str("text", "é").expect("UTF-8 text");
        assert_eq!(
            frame.as_bytes(),
            &[b'M', 0, 0, 0, 3, 0, 0x80, 0xFF, 0, 0, 0, 2, 0xC3, 0xA9]
        );
    }

    #[test]
    fn exact_frame_bound_accepts_and_one_over_refuses_atomically() {
        let mut frame = CanonicalFrameV1::new(b"M", 8).expect("fixture frame");
        frame.push_bytes("exact", &[1, 2, 3]).expect("exact fit");
        assert_eq!(frame.as_bytes(), &[b'M', 0, 0, 0, 3, 1, 2, 3]);

        let before = frame.as_bytes().to_vec();
        let error = frame.push_u8("one_over", 9).expect_err("must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(error.field(), "one_over");
        assert_eq!(error.observed(), "9");
        assert_eq!(frame.as_bytes(), before);

        let mut partial_prefix = CanonicalFrameV1::new(b"M", 6).expect("fixture frame");
        let error = partial_prefix
            .push_bytes("payload", &[1, 2])
            .expect_err("prefix plus payload exceeds the bound");
        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(partial_prefix.as_bytes(), b"M");
    }

    #[test]
    fn count_and_frame_length_overflow_helpers_refuse_precisely() {
        if usize::BITS > u32::BITS {
            let too_many = usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize");
            let error = checked_u32_length("count", too_many).expect_err("u32 overflow");
            assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
            assert_eq!(error.field(), "count");
            assert_eq!(error.expected(), "length representable as u32");
        }

        let error = checked_canonical_frame_length_v1("frame", usize::MAX, 1, usize::MAX)
            .expect_err("usize overflow");
        assert_eq!(error.kind(), ConstructionErrorKindV2::ArithmeticOverflow);
        assert_eq!(error.field(), "frame");
        assert_eq!(error.expected(), "checked frame length");
    }

    #[test]
    fn count_only_preflight_matches_the_exact_encoded_length() {
        let mut size = CanonicalFrameSizeV1::new(b"FRAME", 128).expect("size preflight");
        size.push_u16("code", 7).expect("u16");
        size.push_u32("count", 3).expect("u32");
        size.push_str("name", "alpha").expect("string");
        size.push_presence("optional", true).expect("presence");
        size.push_fixed_bytes_32("root", &[0xA5; 32])
            .expect("fixed root");

        let frame = CanonicalFrameV1::preflighted(b"FRAME", 128, |sink| {
            sink.push_u16("code", 7)?;
            sink.push_u32("count", 3)?;
            sink.push_str("name", "alpha")?;
            sink.push_presence("optional", true)?;
            sink.push_fixed_bytes_32("root", &[0xA5; 32])
        })
        .expect("preflighted frame");

        assert_eq!(size.exact_len(), frame.as_bytes().len());
        assert_eq!(size.exact_len(), 5 + 2 + 4 + 4 + 5 + 1 + 32);
    }

    #[test]
    fn preflight_refusal_happens_before_the_encoding_pass() {
        let invocations = Cell::new(0_u8);
        let error = CanonicalFrameV1::preflighted(b"M", 7, |sink| {
            invocations.set(invocations.get() + 1);
            sink.push_bytes("payload", &[1, 2, 3])
        })
        .expect_err("complete size exceeds the inclusive bound");

        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(error.field(), "payload");
        assert_eq!(invocations.get(), 1);
    }

    #[test]
    fn preflighted_frame_refuses_a_divergent_second_pass() {
        let invocations = Cell::new(0_u8);
        let error = CanonicalFrameV1::preflighted(b"M", 32, |sink| {
            let invocation = invocations.get();
            invocations.set(invocation + 1);
            if invocation == 0 {
                sink.push_u16("value", 7)
            } else {
                sink.push_u8("value", 7)
            }
        })
        .expect_err("second pass must have the preflighted exact length");

        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "canonical.preflight");
        assert_eq!(
            error.expected(),
            "encoded length exactly matches the size preflight"
        );
        assert_eq!(error.observed(), "2");
        assert_eq!(invocations.get(), 2);
    }

    #[test]
    fn preflighted_frame_refuses_growth_beyond_the_exact_capacity() {
        let invocations = Cell::new(0_u8);
        let error = CanonicalFrameV1::preflighted(b"M", 32, |sink| {
            let invocation = invocations.get();
            invocations.set(invocation + 1);
            if invocation == 0 {
                sink.push_u8("value", 7)
            } else {
                sink.push_u16("value", 7)
            }
        })
        .expect_err("second pass cannot grow beyond its preflighted capacity");

        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(error.field(), "value");
        assert_eq!(error.expected(), "value fits the canonical frame bound");
        assert_eq!(error.observed(), "3");
        assert_eq!(invocations.get(), 2);
    }

    #[test]
    fn magic_is_part_of_the_bound() {
        let error = CanonicalFrameV1::new(b"AB", 1).expect_err("oversized magic");
        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(error.field(), "canonical.magic");
        assert_eq!(error.observed(), "2");
    }
}
