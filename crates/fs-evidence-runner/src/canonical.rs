//! Small bounded canonical-frame builder used by pure base-schema projections.

use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use fs_blake3::{ContentHash, hash_domain};

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
        let length = u32::try_from(value.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                field,
                "length representable as u32",
                value.len(),
            )
        })?;
        self.push_u32(field, length)?;
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

    pub(crate) fn root(&self, domain: &'static str) -> ContentHash {
        hash_domain(domain, &self.bytes)
    }

    fn extend(&mut self, field: &'static str, value: &[u8]) -> Result<(), ConstructionErrorV2> {
        let next = self.bytes.len().checked_add(value.len()).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::ArithmeticOverflow,
                field,
                "checked frame length",
                value.len(),
            )
        })?;
        if next > self.max_bytes {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                field,
                "value fits the canonical frame bound",
                next,
            ));
        }
        self.bytes.extend_from_slice(value);
        Ok(())
    }
}
