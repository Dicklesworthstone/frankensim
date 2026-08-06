//! In-house OpenEXR writer + reader (plan §10.5), spec-conformant subset:
//! single-part scanline files, NONE compression, HALF/FLOAT channels,
//! multi-channel AOVs (beauty/albedo/normal/depth) with the spec's
//! alphabetical channel ordering, and opaque custom header attributes. The
//! reader covers our writer's subset (round-trips + ledger artifacts) with
//! structured rejections beyond it.
//! The bounded structural inspector validates that exact subset, including
//! offsets and scanline framing, without materializing image planes.
//!
//! Determinism: byte-exact encodes (pure integer/bit code).

use crate::ImgError;
use std::collections::BTreeSet;
use std::mem::size_of;

/// Canonical custom EXR attribute name for the rendered source artifact hash.
///
/// `fs-img` treats the value as opaque bytes.  The L6 composition layer owns
/// hash syntax and ledger lookup; this L5 crate only preserves the bytes.
pub const SOURCE_ARTIFACT_HASH_ATTRIBUTE: &str = "frankensim.sourceArtifactHash";

/// Channel sample type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PixelType {
    /// 16-bit half float.
    Half,
    /// 32-bit float.
    Float,
}

impl PixelType {
    fn code(self) -> u32 {
        match self {
            PixelType::Half => 1,
            PixelType::Float => 2,
        }
    }

    fn bytes(self) -> usize {
        match self {
            PixelType::Half => 2,
            PixelType::Float => 4,
        }
    }
}

/// One named planar channel (row-major f32 samples; converted on write).
#[derive(Debug, Clone, PartialEq)]
pub struct Channel {
    /// Channel name (e.g. "R", "albedo.G", "depth.Z"); 1..=31 UTF-8
    /// bytes with no NULs for the short-name-only version-2 subset.
    pub name: String,
    /// Storage type.
    pub ty: PixelType,
    /// Row-major samples (width × height).
    pub data: Vec<f32>,
}

/// One opaque custom OpenEXR header attribute.
///
/// Attribute names and type names are UTF-8 strings admitted by the writer;
/// values are preserved byte-for-byte and may contain NUL or non-UTF-8 data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExrAttribute {
    /// Attribute name, unique within the custom attribute set.
    pub name: String,
    /// OpenEXR attribute type name, such as `string`.
    pub ty: String,
    /// Raw attribute payload bytes.
    pub value: Vec<u8>,
}

/// Caller-owned ceilings for the allocations performed by the EXR writer.
///
/// Scratch covers logical reference slots used to sort channels and
/// attributes. Output covers the exact encoded artifact length. Neither value
/// claims allocator bookkeeping or a process-RSS ceiling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExrWriteLimits {
    /// Maximum logical bytes in sorting-reference storage.
    pub max_scratch_bytes: u64,
    /// Maximum exact encoded output bytes.
    pub max_output_bytes: u64,
}

impl ExrWriteLimits {
    /// Unlimited logical ceilings used by the compatibility writer APIs.
    pub const UNBOUNDED: Self = Self {
        max_scratch_bytes: u64::MAX,
        max_output_bytes: u64::MAX,
    };
}

/// Exact logical storage required by one deterministic EXR encode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExrWriteRequirements {
    /// Reference-vector payload used to establish canonical ordering.
    pub scratch_bytes: u64,
    /// Complete encoded artifact length, including headers and offset table.
    pub output_bytes: u64,
}

/// f32 → f16 bits with round-to-nearest-even (subnormals + specials).
#[must_use]
pub fn f32_to_f16_bits(x: f32) -> u16 {
    let bits = x.to_bits();
    let sign = ((bits >> 16) & 0x8000) as u16;
    let exp = (bits >> 23) & 0xFF;
    let man = bits & 0x007F_FFFF;
    if exp == 0xFF {
        // Inf / NaN (keep a NaN payload bit).
        return sign | 0x7C00 | u16::from(man != 0) << 9;
    }
    let unbiased = exp.cast_signed() - 127;
    if unbiased > 15 {
        return sign | 0x7C00; // overflow → ±inf
    }
    if unbiased >= -14 {
        // Normal half: 10-bit mantissa with RNE.
        let mut half = sign | (((unbiased + 15) as u16) << 10) | ((man >> 13) as u16);
        let rem = man & 0x1FFF;
        if rem > 0x1000 || (rem == 0x1000 && (half & 1) == 1) {
            half += 1; // carries correctly into the exponent
        }
        return half;
    }
    if unbiased < -25 {
        return sign; // underflow → ±0
    }
    // Subnormal half.
    let full = man | 0x0080_0000; // implicit bit
    let shift = (-14 - unbiased + 13) as u32;
    let mut half = sign | ((full >> shift) as u16);
    let rem = full & ((1u32 << shift) - 1);
    let halfway = 1u32 << (shift - 1);
    if rem > halfway || (rem == halfway && (half & 1) == 1) {
        half += 1;
    }
    half
}

/// f16 bits → f32.
#[must_use]
pub fn f16_bits_to_f32(h: u16) -> f32 {
    let sign = u32::from(h >> 15) << 31;
    let exp = u32::from(h >> 10) & 0x1F;
    let man = u32::from(h) & 0x3FF;
    let bits = match (exp, man) {
        (0, 0) => sign,
        (0, m) => {
            // Subnormal half = m·2⁻²⁴: normalize around the highest bit.
            let h = m.ilog2();
            let e = 103 + h; // 127 + (h − 24)
            let mant = (m << (23 - h)) & 0x007F_FFFF;
            sign | (e << 23) | mant
        }
        (0x1F, 0) => sign | 0x7F80_0000,
        (0x1F, m) => sign | 0x7F80_0000 | (m << 13),
        (e, m) => sign | ((e + 127 - 15) << 23) | (m << 13),
    };
    f32::from_bits(bits)
}

const MAGIC: [u8; 4] = [0x76, 0x2F, 0x31, 0x01];

fn push_attr(out: &mut Vec<u8>, name: &str, ty: &str, value: &[u8]) {
    out.extend_from_slice(name.as_bytes());
    out.push(0);
    out.extend_from_slice(ty.as_bytes());
    out.push(0);
    out.extend_from_slice(
        &i32::try_from(value.len())
            .expect("EXR attribute length was admitted before encoding")
            .to_le_bytes(),
    );
    out.extend_from_slice(value);
}

/// Encode channels as a single-part scanline EXR (NONE compression), without
/// custom attributes. Channels are stored in the spec's alphabetical order
/// regardless of the argument order.
///
/// # Errors
/// [`ImgError`] on shape/name defects.
pub fn write_exr(width: u32, height: u32, channels: &[Channel]) -> Result<Vec<u8>, ImgError> {
    write_exr_with_attributes(width, height, channels, &[])
}

/// Encode channels and opaque custom header attributes as a single-part
/// scanline EXR (NONE compression).
///
/// Channels and custom attributes are each stored in alphabetical name order.
/// Passing an empty attribute slice is byte-identical to [`write_exr`].
/// Built-in EXR attribute names are reserved so custom values cannot shadow
/// the structural header.
///
/// # Errors
/// [`ImgError`] on shape, name, duplicate, reserved-name, or size defects.
pub fn write_exr_with_attributes(
    width: u32,
    height: u32,
    channels: &[Channel],
    attributes: &[ExrAttribute],
) -> Result<Vec<u8>, ImgError> {
    write_exr_with_attributes_budgeted(
        width,
        height,
        channels,
        attributes,
        ExrWriteLimits::UNBOUNDED,
    )
}

/// Validate an EXR request and return its exact logical scratch and encoded
/// output requirements without allocating image-sized storage.
///
/// # Errors
/// [`ImgError`] on the same shape, name, duplicate, reserved-name, or checked
/// size defects as [`write_exr_with_attributes`].
pub fn exr_write_requirements(
    width: u32,
    height: u32,
    channels: &[Channel],
    attributes: &[ExrAttribute],
) -> Result<ExrWriteRequirements, ImgError> {
    let pixel_count = checked_pixel_count(width, height, channels.len())?;
    for channel in channels {
        if channel.data.len() != pixel_count {
            return Err(ImgError::Shape {
                expected: pixel_count,
                got: channel.data.len(),
                context: "channel sample count",
            });
        }
    }
    exr_layout_requirements_impl(
        width,
        height,
        channels.len(),
        |index| (channels[index].name.as_str(), channels[index].ty),
        attributes,
    )
}

/// Validate a channel layout and return exact EXR requirements without owning
/// or allocating any pixel planes. This is the admission seam for producers
/// that must reject an output ceiling before constructing image-sized data.
///
/// # Errors
/// [`ImgError`] on invalid dimensions, names, duplicates, attributes, or
/// checked size overflow.
pub fn exr_write_requirements_for_layout(
    width: u32,
    height: u32,
    channels: &[(&str, PixelType)],
    attributes: &[ExrAttribute],
) -> Result<ExrWriteRequirements, ImgError> {
    checked_pixel_count(width, height, channels.len())?;
    exr_layout_requirements_impl(
        width,
        height,
        channels.len(),
        |index| channels[index],
        attributes,
    )
}

fn checked_pixel_count(width: u32, height: u32, channel_count: usize) -> Result<usize, ImgError> {
    if width == 0 || height == 0 || channel_count == 0 {
        return Err(ImgError::Shape {
            expected: 1,
            got: 0,
            context: "write_exr needs a nonempty image and channel set",
        });
    }
    if i32::try_from(width).is_err() || i32::try_from(height).is_err() {
        return Err(ImgError::Malformed {
            what: format!("dimensions {width}x{height} exceed the EXR i32 data window"),
        });
    }
    let pixel_count_u64 =
        u64::from(width)
            .checked_mul(u64::from(height))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR pixel count",
            })?;
    let pixel_count = usize::try_from(pixel_count_u64).map_err(|_| ImgError::SizeOverflow {
        context: "EXR pixel count on this platform",
    })?;
    Ok(pixel_count)
}

#[allow(clippy::too_many_lines)] // one checked-size derivation mirrors the wire layout
fn exr_layout_requirements_impl<'a>(
    width: u32,
    height: u32,
    channel_count: usize,
    mut channel_at: impl FnMut(usize) -> (&'a str, PixelType),
    attributes: &[ExrAttribute],
) -> Result<ExrWriteRequirements, ImgError> {
    let mut channel_value_bytes = 1_u64;
    let mut line_bytes = 0_u64;
    for index in 0..channel_count {
        let (name, ty) = channel_at(index);
        validate_attribute_name("channel", name)?;
        if (0..index).any(|prior| channel_at(prior).0 == name) {
            return Err(ImgError::Malformed {
                what: format!("duplicate channel {name:?}"),
            });
        }
        channel_value_bytes = channel_value_bytes
            .checked_add(
                u64::try_from(name.len()).map_err(|_| ImgError::SizeOverflow {
                    context: "EXR channel name length",
                })?,
            )
            .and_then(|bytes| bytes.checked_add(17))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR channel-list bytes",
            })?;
        line_bytes = line_bytes
            .checked_add(u64::from(width).checked_mul(ty.bytes() as u64).ok_or(
                ImgError::SizeOverflow {
                    context: "EXR scanline bytes",
                },
            )?)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR scanline bytes",
            })?;
    }
    if i32::try_from(channel_value_bytes).is_err() || i32::try_from(line_bytes).is_err() {
        return Err(ImgError::SizeOverflow {
            context: "EXR signed attribute or scanline length",
        });
    }

    let mut custom_attribute_bytes = 0_u64;
    for (index, attribute) in attributes.iter().enumerate() {
        validate_attribute_name("attribute", &attribute.name)?;
        validate_attribute_name("attribute type", &attribute.ty)?;
        if is_builtin_attribute(&attribute.name) {
            return Err(ImgError::Malformed {
                what: format!("custom attribute {:?} shadows a built-in", attribute.name),
            });
        }
        if i32::try_from(attribute.value.len()).is_err() {
            return Err(ImgError::Shape {
                expected: usize::try_from(i32::MAX).unwrap_or(usize::MAX),
                got: attribute.value.len(),
                context: "EXR attribute payload bytes (maximum)",
            });
        }
        if attributes[..index]
            .iter()
            .any(|prior| prior.name == attribute.name)
        {
            return Err(ImgError::Malformed {
                what: format!("duplicate custom attribute {:?}", attribute.name),
            });
        }
        custom_attribute_bytes = custom_attribute_bytes
            .checked_add(encoded_attribute_len(
                &attribute.name,
                &attribute.ty,
                u64::try_from(attribute.value.len()).map_err(|_| ImgError::SizeOverflow {
                    context: "EXR custom attribute payload",
                })?,
            )?)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR custom attributes",
            })?;
    }

    let mut output_bytes = 8_u64;
    output_bytes = checked_add_size(
        output_bytes,
        encoded_attribute_len("channels", "chlist", channel_value_bytes)?,
        "EXR header",
    )?;
    for (name, ty, value_bytes) in [
        ("compression", "compression", 1_u64),
        ("dataWindow", "box2i", 16),
        ("displayWindow", "box2i", 16),
        ("lineOrder", "lineOrder", 1),
        ("pixelAspectRatio", "float", 4),
        ("screenWindowCenter", "v2f", 8),
        ("screenWindowWidth", "float", 4),
    ] {
        output_bytes = checked_add_size(
            output_bytes,
            encoded_attribute_len(name, ty, value_bytes)?,
            "EXR header",
        )?;
    }
    output_bytes = checked_add_size(output_bytes, custom_attribute_bytes, "EXR header")?;
    output_bytes = checked_add_size(output_bytes, 1, "EXR header terminator")?;
    output_bytes = checked_add_size(
        output_bytes,
        u64::from(height)
            .checked_mul(8)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR offset table",
            })?,
        "EXR offset table",
    )?;
    output_bytes = checked_add_size(
        output_bytes,
        u64::from(height)
            .checked_mul(line_bytes.checked_add(8).ok_or(ImgError::SizeOverflow {
                context: "EXR scanline block",
            })?)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR scanline blocks",
            })?,
        "EXR scanline blocks",
    )?;
    usize::try_from(output_bytes).map_err(|_| ImgError::SizeOverflow {
        context: "EXR output on this platform",
    })?;

    let scratch_bytes = u64::try_from(channel_count)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<&Channel>() as u64))
        .and_then(|bytes| {
            u64::try_from(attributes.len())
                .ok()
                .and_then(|count| count.checked_mul(size_of::<&ExrAttribute>() as u64))
                .and_then(|attribute_bytes| bytes.checked_add(attribute_bytes))
        })
        .ok_or(ImgError::SizeOverflow {
            context: "EXR ordering scratch",
        })?;
    Ok(ExrWriteRequirements {
        scratch_bytes,
        output_bytes,
    })
}

/// Encode one deterministic EXR after admitting all writer-owned logical
/// scratch and exact output bytes against caller ceilings.
///
/// Input channel planes and custom-attribute payloads remain caller-owned and
/// are therefore outside these two limits.
///
/// # Errors
/// [`ImgError`] on input defects, ceiling refusal, allocation refusal, or
/// checked size overflow.
#[allow(clippy::too_many_lines)] // one linear encoder keeps admission and wire order reviewable
pub fn write_exr_with_attributes_budgeted(
    width: u32,
    height: u32,
    channels: &[Channel],
    attributes: &[ExrAttribute],
    limits: ExrWriteLimits,
) -> Result<Vec<u8>, ImgError> {
    let requirements = exr_write_requirements(width, height, channels, attributes)?;
    if requirements.scratch_bytes > limits.max_scratch_bytes {
        return Err(ImgError::ResourceLimit {
            resource: "EXR ordering scratch",
            requested: requirements.scratch_bytes,
            limit: limits.max_scratch_bytes,
        });
    }
    if requirements.output_bytes > limits.max_output_bytes {
        return Err(ImgError::ResourceLimit {
            resource: "EXR encoded output",
            requested: requirements.output_bytes,
            limit: limits.max_output_bytes,
        });
    }

    let wi = i32::try_from(width).map_err(|_| ImgError::SizeOverflow {
        context: "EXR width",
    })?;
    let hi = i32::try_from(height).map_err(|_| ImgError::SizeOverflow {
        context: "EXR height",
    })?;
    let mut sorted = Vec::new();
    sorted
        .try_reserve_exact(channels.len())
        .map_err(|_| ImgError::AllocationRefused {
            resource: "EXR channel ordering scratch",
            requested: u64::try_from(channels.len())
                .unwrap_or(u64::MAX)
                .saturating_mul(size_of::<&Channel>() as u64),
        })?;
    sorted.extend(channels);
    sorted.sort_unstable_by(|left, right| left.name.cmp(&right.name));

    let mut sorted_attributes = Vec::new();
    sorted_attributes
        .try_reserve_exact(attributes.len())
        .map_err(|_| ImgError::AllocationRefused {
            resource: "EXR attribute ordering scratch",
            requested: u64::try_from(attributes.len())
                .unwrap_or(u64::MAX)
                .saturating_mul(size_of::<&ExrAttribute>() as u64),
        })?;
    sorted_attributes.extend(attributes);
    sorted_attributes.sort_unstable_by(|left, right| left.name.cmp(&right.name));

    let output_capacity =
        usize::try_from(requirements.output_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR output capacity",
        })?;
    let mut out = Vec::new();
    out.try_reserve_exact(output_capacity)
        .map_err(|_| ImgError::AllocationRefused {
            resource: "EXR encoded output",
            requested: requirements.output_bytes,
        })?;
    out.extend_from_slice(&MAGIC);
    out.extend_from_slice(&2u32.to_le_bytes()); // version 2, no flags

    let channel_value_bytes = sorted.iter().try_fold(1_u64, |bytes, channel| {
        bytes
            .checked_add(channel.name.len() as u64)
            .and_then(|value| value.checked_add(17))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR channel-list bytes",
            })
    })?;
    push_attr_prefix(
        &mut out,
        "channels",
        "chlist",
        i32::try_from(channel_value_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR channel-list signed length",
        })?,
    );
    for c in &sorted {
        out.extend_from_slice(c.name.as_bytes());
        out.push(0);
        out.extend_from_slice(&c.ty.code().to_le_bytes());
        out.extend_from_slice(&[0, 0, 0, 0]); // pLinear + reserved
        out.extend_from_slice(&1u32.to_le_bytes()); // xSampling
        out.extend_from_slice(&1u32.to_le_bytes()); // ySampling
    }
    out.push(0);
    push_attr(&mut out, "compression", "compression", &[0]); // NONE
    let mut window = [0_u8; 16];
    window[8..12].copy_from_slice(&(wi - 1).to_le_bytes());
    window[12..16].copy_from_slice(&(hi - 1).to_le_bytes());
    push_attr(&mut out, "dataWindow", "box2i", &window);
    push_attr(&mut out, "displayWindow", "box2i", &window);
    push_attr(&mut out, "lineOrder", "lineOrder", &[0]); // increasing y
    push_attr(&mut out, "pixelAspectRatio", "float", &1.0f32.to_le_bytes());
    push_attr(&mut out, "screenWindowCenter", "v2f", &[0u8; 8]);
    push_attr(
        &mut out,
        "screenWindowWidth",
        "float",
        &1.0f32.to_le_bytes(),
    );
    for attribute in &sorted_attributes {
        push_attr(&mut out, &attribute.name, &attribute.ty, &attribute.value);
    }
    out.push(0); // end of header

    // Scanline offset table placeholder.
    let table_pos = out.len();
    out.resize(out.len() + 8 * height as usize, 0);

    let line_bytes = sorted.iter().try_fold(0_usize, |bytes, channel| {
        (width as usize)
            .checked_mul(channel.ty.bytes())
            .and_then(|row| bytes.checked_add(row))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR scanline bytes",
            })
    })?;
    let line_bytes_i32 = i32::try_from(line_bytes).map_err(|_| ImgError::SizeOverflow {
        context: "EXR scanline signed length",
    })?;
    for y in 0..height as usize {
        let offset = out.len() as u64;
        out[table_pos + 8 * y..table_pos + 8 * (y + 1)].copy_from_slice(&offset.to_le_bytes());
        out.extend_from_slice(
            &i32::try_from(y)
                .expect("y < height <= i32::MAX")
                .to_le_bytes(),
        );
        out.extend_from_slice(&line_bytes_i32.to_le_bytes());
        for c in &sorted {
            let row = &c.data[y * width as usize..(y + 1) * width as usize];
            match c.ty {
                PixelType::Half => {
                    for &v in row {
                        out.extend_from_slice(&f32_to_f16_bits(v).to_le_bytes());
                    }
                }
                PixelType::Float => {
                    for &v in row {
                        out.extend_from_slice(&v.to_le_bytes());
                    }
                }
            }
        }
    }
    if out.len() != output_capacity {
        return Err(ImgError::SizeOverflow {
            context: "EXR encoded-length invariant",
        });
    }
    Ok(out)
}

fn checked_add_size(left: u64, right: u64, context: &'static str) -> Result<u64, ImgError> {
    left.checked_add(right)
        .ok_or(ImgError::SizeOverflow { context })
}

fn encoded_attribute_len(name: &str, ty: &str, value_bytes: u64) -> Result<u64, ImgError> {
    if value_bytes > i32::MAX as u64 {
        return Err(ImgError::SizeOverflow {
            context: "EXR attribute signed length",
        });
    }
    u64::try_from(name.len())
        .ok()
        .and_then(|bytes| bytes.checked_add(1))
        .and_then(|bytes| bytes.checked_add(u64::try_from(ty.len()).ok()?))
        .and_then(|bytes| bytes.checked_add(1 + 4))
        .and_then(|bytes| bytes.checked_add(value_bytes))
        .ok_or(ImgError::SizeOverflow {
            context: "EXR attribute encoding",
        })
}

fn push_attr_prefix(out: &mut Vec<u8>, name: &str, ty: &str, value_len: i32) {
    out.extend_from_slice(name.as_bytes());
    out.push(0);
    out.extend_from_slice(ty.as_bytes());
    out.push(0);
    out.extend_from_slice(&value_len.to_le_bytes());
}

fn validate_attribute_name(context: &str, name: &str) -> Result<(), ImgError> {
    if name.is_empty() || name.contains('\0') || name.len() > 31 {
        return Err(ImgError::Malformed {
            what: format!("{context} name {name:?} must contain 1..=31 UTF-8 bytes and no NUL"),
        });
    }
    Ok(())
}

fn is_builtin_attribute(name: &str) -> bool {
    matches!(
        name,
        "channels"
            | "compression"
            | "dataWindow"
            | "displayWindow"
            | "lineOrder"
            | "pixelAspectRatio"
            | "screenWindowCenter"
            | "screenWindowWidth"
    )
}

/// A decoded EXR (our writer's subset): alphabetical channels, f32 data
/// (HALF widened losslessly).
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedExr {
    /// Pixel width.
    pub width: u32,
    /// Pixel height.
    pub height: u32,
    /// Channels in file (alphabetical) order.
    pub channels: Vec<Channel>,
    /// Opaque custom header attributes in file order.
    pub attributes: Vec<ExrAttribute>,
}

/// Caller-owned byte ceilings for structural EXR inspection.
///
/// Inspection never materializes channel planes. `max_decoded_bytes` bounds
/// the logical uncompressed sample payload described by the file, while
/// `max_metadata_bytes` bounds the two small borrowed-descriptor vectors
/// returned to the caller.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExrInspectLimits {
    /// Maximum encoded artifact bytes accepted at entry.
    pub max_input_bytes: u64,
    /// Maximum bytes through the header terminator.
    pub max_header_bytes: u64,
    /// Maximum logical uncompressed channel-sample bytes.
    pub max_decoded_bytes: u64,
    /// Maximum logical bytes in returned borrowed metadata descriptors.
    pub max_metadata_bytes: u64,
}

impl ExrInspectLimits {
    /// Unlimited ceilings for callers that already enforce an outer budget.
    pub const UNBOUNDED: Self = Self {
        max_input_bytes: u64::MAX,
        max_header_bytes: u64::MAX,
        max_decoded_bytes: u64::MAX,
        max_metadata_bytes: u64::MAX,
    };
}

/// One channel descriptor borrowed directly from an inspected EXR header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExrInspectedChannel<'a> {
    /// Canonical short channel name.
    pub name: &'a str,
    /// On-wire scalar type.
    pub ty: PixelType,
}

/// One custom attribute borrowed directly from an inspected EXR header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExrInspectedAttribute<'a> {
    /// Canonical short attribute name.
    pub name: &'a str,
    /// Canonical short attribute type name.
    pub ty: &'a str,
    /// Opaque attribute payload.
    pub value: &'a [u8],
}

/// Structural facts for an EXR in the exact subset emitted by this crate.
///
/// Channel samples are not decoded or retained. The borrowed channel and
/// custom-attribute views remain valid for the lifetime of `bytes` supplied to
/// [`inspect_exr`] or [`inspect_exr_with_poll`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExrInspection<'a> {
    /// Pixel width.
    pub width: u32,
    /// Pixel height.
    pub height: u32,
    /// Alphabetically ordered channel schema.
    pub channels: Vec<ExrInspectedChannel<'a>>,
    /// Alphabetically ordered custom attributes.
    pub attributes: Vec<ExrInspectedAttribute<'a>>,
    /// Complete encoded artifact length.
    pub input_bytes: u64,
    /// Bytes through the header terminator.
    pub header_bytes: u64,
    /// Bytes in one uncompressed scanline payload.
    pub scanline_bytes: u64,
    /// Total logical uncompressed channel-sample bytes.
    pub decoded_bytes: u64,
    /// Logical bytes in the returned descriptor vectors.
    pub metadata_bytes: u64,
}

#[derive(Debug, Clone, Copy)]
struct ExrAttributeRef<'a> {
    name: &'a str,
    ty: &'a str,
    value: &'a [u8],
}

#[derive(Debug, Clone, Copy)]
struct StrictExrHeader<'a> {
    channels_value: &'a [u8],
    width: u32,
    height: u32,
    header_bytes: u64,
    scanline_bytes: u64,
    decoded_bytes: u64,
    channel_count: usize,
    custom_attribute_count: usize,
    metadata_bytes: u64,
}

fn inspect_continue(
    poll: &mut impl FnMut() -> bool,
    operation: &'static str,
) -> Result<(), ImgError> {
    if poll() {
        Ok(())
    } else {
        Err(ImgError::Cancelled { operation })
    }
}

fn admit_inspect_bytes(resource: &'static str, requested: u64, limit: u64) -> Result<(), ImgError> {
    if requested > limit {
        Err(ImgError::ResourceLimit {
            resource,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

fn read_short_cstr<'a>(
    bytes: &'a [u8],
    pos: &mut usize,
    context: &str,
) -> Result<&'a str, ImgError> {
    let start = *pos;
    let remaining = bytes.get(start..).ok_or_else(|| ImgError::Malformed {
        what: format!("truncated {context} at byte {start}"),
    })?;
    let searchable = &remaining[..remaining.len().min(32)];
    let Some(length) = searchable.iter().position(|&byte| byte == 0) else {
        let what = if remaining.len() <= 31 {
            format!("unterminated {context} at byte {start}")
        } else {
            format!("{context} at byte {start} exceeds the version-2 31-byte limit")
        };
        return Err(ImgError::Malformed { what });
    };
    if length == 0 {
        return Err(ImgError::Malformed {
            what: format!("empty {context} at byte {start}"),
        });
    }
    let value = std::str::from_utf8(&remaining[..length]).map_err(|_| ImgError::Malformed {
        what: format!("non-UTF-8 {context} at byte {start}"),
    })?;
    *pos = start
        .checked_add(length + 1)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR header string offset",
        })?;
    Ok(value)
}

fn read_attribute_ref<'a>(
    bytes: &'a [u8],
    pos: &mut usize,
    max_header_bytes: u64,
) -> Result<ExrAttributeRef<'a>, ImgError> {
    let name = read_short_cstr(bytes, pos, "EXR attribute name")?;
    let ty = read_short_cstr(bytes, pos, "EXR attribute type")?;
    let signed_size = i32::from_le_bytes(take(bytes, *pos, 4)?.try_into().expect("4 bytes"));
    *pos = pos.checked_add(4).ok_or(ImgError::SizeOverflow {
        context: "EXR attribute payload offset",
    })?;
    let size = usize::try_from(signed_size).map_err(|_| ImgError::Malformed {
        what: format!("negative payload length for EXR attribute {name:?}"),
    })?;
    let end = pos.checked_add(size).ok_or(ImgError::SizeOverflow {
        context: "EXR attribute payload offset",
    })?;
    admit_inspect_bytes(
        "EXR header",
        u64::try_from(end).map_err(|_| ImgError::SizeOverflow {
            context: "EXR header bytes",
        })?,
        max_header_bytes,
    )?;
    let value = take(bytes, *pos, size)?;
    *pos = end;
    Ok(ExrAttributeRef { name, ty, value })
}

fn require_builtin_attribute(
    attribute: ExrAttributeRef<'_>,
    expected_name: &str,
    expected_type: &str,
) -> Result<(), ImgError> {
    if attribute.name != expected_name || attribute.ty != expected_type {
        return Err(ImgError::Malformed {
            what: format!(
                "canonical EXR header expected {expected_name:?}/{expected_type:?}, found {:?}/{:?}",
                attribute.name, attribute.ty
            ),
        });
    }
    Ok(())
}

fn inspect_chlist(value: &[u8], poll: &mut impl FnMut() -> bool) -> Result<(usize, u64), ImgError> {
    let mut pos = 0usize;
    let mut count = 0usize;
    let mut bytes_per_pixel = 0u64;
    let mut previous_name: Option<&str> = None;
    loop {
        inspect_continue(poll, "EXR structural inspection")?;
        let Some(&next) = value.get(pos) else {
            return Err(ImgError::Malformed {
                what: "EXR channel list lacks its terminal NUL".to_string(),
            });
        };
        if next == 0 {
            pos += 1;
            if pos != value.len() {
                return Err(ImgError::Malformed {
                    what: "trailing bytes after EXR channel-list terminator".to_string(),
                });
            }
            break;
        }
        let name = read_short_cstr(value, &mut pos, "EXR channel name")?;
        if previous_name.is_some_and(|previous| previous >= name) {
            return Err(ImgError::Malformed {
                what: format!("EXR channels are not in strict alphabetical order at {name:?}"),
            });
        }
        let descriptor = take(value, pos, 16)?;
        let code = u32::from_le_bytes(descriptor[..4].try_into().expect("4 bytes"));
        let ty = match code {
            1 => PixelType::Half,
            2 => PixelType::Float,
            other => {
                return Err(ImgError::Unsupported {
                    what: format!("pixel type {other}"),
                });
            }
        };
        if descriptor[4..8] != [0, 0, 0, 0]
            || u32::from_le_bytes(descriptor[8..12].try_into().expect("4 bytes")) != 1
            || u32::from_le_bytes(descriptor[12..16].try_into().expect("4 bytes")) != 1
        {
            return Err(ImgError::Unsupported {
                what: format!(
                    "channel {name:?} uses noncanonical pLinear/reserved or subsampling fields"
                ),
            });
        }
        pos += 16;
        count = count.checked_add(1).ok_or(ImgError::SizeOverflow {
            context: "EXR channel count",
        })?;
        bytes_per_pixel =
            bytes_per_pixel
                .checked_add(ty.bytes() as u64)
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR bytes per pixel",
                })?;
        previous_name = Some(name);
    }
    if count == 0 {
        return Err(ImgError::Malformed {
            what: "EXR channel list is empty".to_string(),
        });
    }
    Ok((count, bytes_per_pixel))
}

fn canonical_window(value: &[u8]) -> Result<(u32, u32), ImgError> {
    if value.len() != 16 {
        return Err(ImgError::Malformed {
            what: "canonical EXR dataWindow must contain 16 bytes".to_string(),
        });
    }
    let min_x = i32::from_le_bytes(value[0..4].try_into().expect("4 bytes"));
    let min_y = i32::from_le_bytes(value[4..8].try_into().expect("4 bytes"));
    let max_x = i32::from_le_bytes(value[8..12].try_into().expect("4 bytes"));
    let max_y = i32::from_le_bytes(value[12..16].try_into().expect("4 bytes"));
    if min_x != 0 || min_y != 0 || max_x < 0 || max_y < 0 || max_x == i32::MAX || max_y == i32::MAX
    {
        return Err(ImgError::Unsupported {
            what: "EXR data window must start at (0, 0) within writer-admitted i32 dimensions"
                .to_string(),
        });
    }
    let width = u32::try_from(max_x)
        .ok()
        .and_then(|extent| extent.checked_add(1))
        .ok_or(ImgError::SizeOverflow {
            context: "EXR data-window width",
        })?;
    let height = u32::try_from(max_y)
        .ok()
        .and_then(|extent| extent.checked_add(1))
        .ok_or(ImgError::SizeOverflow {
            context: "EXR data-window height",
        })?;
    Ok((width, height))
}

#[allow(clippy::too_many_lines)]
fn inspect_strict_header<'a>(
    bytes: &'a [u8],
    limits: ExrInspectLimits,
    poll: &mut impl FnMut() -> bool,
) -> Result<StrictExrHeader<'a>, ImgError> {
    let mut pos = 8usize;
    let channels = read_attribute_ref(bytes, &mut pos, limits.max_header_bytes)?;
    require_builtin_attribute(channels, "channels", "chlist")?;
    let (channel_count, bytes_per_pixel) = inspect_chlist(channels.value, poll)?;

    inspect_continue(poll, "EXR structural inspection")?;
    let compression = read_attribute_ref(bytes, &mut pos, limits.max_header_bytes)?;
    require_builtin_attribute(compression, "compression", "compression")?;
    if compression.value != [0] {
        return Err(ImgError::Unsupported {
            what: "EXR compression must be canonical NONE".to_string(),
        });
    }

    let data_window = read_attribute_ref(bytes, &mut pos, limits.max_header_bytes)?;
    require_builtin_attribute(data_window, "dataWindow", "box2i")?;
    let (width, height) = canonical_window(data_window.value)?;

    let display_window = read_attribute_ref(bytes, &mut pos, limits.max_header_bytes)?;
    require_builtin_attribute(display_window, "displayWindow", "box2i")?;
    if display_window.value != data_window.value {
        return Err(ImgError::Unsupported {
            what: "EXR displayWindow must equal dataWindow".to_string(),
        });
    }

    for (name, ty, expected_value) in [
        ("lineOrder", "lineOrder", &[0][..]),
        ("pixelAspectRatio", "float", &1.0f32.to_le_bytes()[..]),
        ("screenWindowCenter", "v2f", &[0; 8][..]),
        ("screenWindowWidth", "float", &1.0f32.to_le_bytes()[..]),
    ] {
        inspect_continue(poll, "EXR structural inspection")?;
        let attribute = read_attribute_ref(bytes, &mut pos, limits.max_header_bytes)?;
        require_builtin_attribute(attribute, name, ty)?;
        if attribute.value != expected_value {
            return Err(ImgError::Unsupported {
                what: format!("EXR {name} is outside the canonical writer subset"),
            });
        }
    }

    let mut custom_attribute_count = 0usize;
    let mut previous_custom_name: Option<&str> = None;
    loop {
        inspect_continue(poll, "EXR structural inspection")?;
        let Some(&next) = bytes.get(pos) else {
            return Err(ImgError::Malformed {
                what: "EXR header lacks its terminal NUL".to_string(),
            });
        };
        if next == 0 {
            pos += 1;
            break;
        }
        let attribute = read_attribute_ref(bytes, &mut pos, limits.max_header_bytes)?;
        if is_builtin_attribute(attribute.name) {
            return Err(ImgError::Malformed {
                what: format!(
                    "duplicate or misplaced built-in EXR attribute {:?}",
                    attribute.name
                ),
            });
        }
        if previous_custom_name.is_some_and(|previous| previous >= attribute.name) {
            return Err(ImgError::Malformed {
                what: format!(
                    "custom EXR attributes are not in strict alphabetical order at {:?}",
                    attribute.name
                ),
            });
        }
        custom_attribute_count =
            custom_attribute_count
                .checked_add(1)
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR custom attribute count",
                })?;
        previous_custom_name = Some(attribute.name);
    }

    let header_bytes = u64::try_from(pos).map_err(|_| ImgError::SizeOverflow {
        context: "EXR header bytes",
    })?;
    admit_inspect_bytes("EXR header", header_bytes, limits.max_header_bytes)?;
    let scanline_bytes =
        u64::from(width)
            .checked_mul(bytes_per_pixel)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR scanline bytes",
            })?;
    if scanline_bytes > i32::MAX as u64 {
        return Err(ImgError::Malformed {
            what: "EXR scanline payload exceeds its signed 32-bit field".to_string(),
        });
    }
    let decoded_bytes =
        scanline_bytes
            .checked_mul(u64::from(height))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR decoded sample bytes",
            })?;
    admit_inspect_bytes(
        "EXR decoded samples",
        decoded_bytes,
        limits.max_decoded_bytes,
    )?;
    let metadata_bytes = u64::try_from(channel_count)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<ExrInspectedChannel<'static>>() as u64))
        .and_then(|channel_bytes| {
            u64::try_from(custom_attribute_count)
                .ok()
                .and_then(|count| {
                    count.checked_mul(size_of::<ExrInspectedAttribute<'static>>() as u64)
                })
                .and_then(|attribute_bytes| channel_bytes.checked_add(attribute_bytes))
        })
        .ok_or(ImgError::SizeOverflow {
            context: "EXR inspection metadata",
        })?;
    admit_inspect_bytes(
        "EXR inspection metadata",
        metadata_bytes,
        limits.max_metadata_bytes,
    )?;

    let offset_table_bytes = u64::from(height)
        .checked_mul(8)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR offset-table bytes",
        })?;
    let blocks_bytes = scanline_bytes
        .checked_add(8)
        .and_then(|bytes| bytes.checked_mul(u64::from(height)))
        .ok_or(ImgError::SizeOverflow {
            context: "EXR scanline-block bytes",
        })?;
    let expected_input = header_bytes
        .checked_add(offset_table_bytes)
        .and_then(|bytes| bytes.checked_add(blocks_bytes))
        .ok_or(ImgError::SizeOverflow {
            context: "EXR encoded bytes",
        })?;
    let actual_input = u64::try_from(bytes.len()).map_err(|_| ImgError::SizeOverflow {
        context: "EXR input bytes",
    })?;
    if actual_input != expected_input {
        return Err(ImgError::Malformed {
            what: format!(
                "EXR encoded length is {actual_input} bytes; canonical layout requires {expected_input}"
            ),
        });
    }

    Ok(StrictExrHeader {
        channels_value: channels.value,
        width,
        height,
        header_bytes,
        scanline_bytes,
        decoded_bytes,
        channel_count,
        custom_attribute_count,
        metadata_bytes,
    })
}

fn collect_inspected_channels<'a>(
    value: &'a [u8],
    count: usize,
    poll: &mut impl FnMut() -> bool,
) -> Result<Vec<ExrInspectedChannel<'a>>, ImgError> {
    let requested = u64::try_from(count)
        .unwrap_or(u64::MAX)
        .saturating_mul(size_of::<ExrInspectedChannel<'static>>() as u64);
    let mut channels = Vec::new();
    channels
        .try_reserve_exact(count)
        .map_err(|_| ImgError::AllocationRefused {
            resource: "EXR inspection channel descriptors",
            requested,
        })?;
    let mut pos = 0usize;
    while value.get(pos) != Some(&0) {
        inspect_continue(poll, "EXR structural inspection")?;
        let name = read_short_cstr(value, &mut pos, "EXR channel name")?;
        let descriptor = take(value, pos, 16)?;
        let ty = match u32::from_le_bytes(descriptor[..4].try_into().expect("4 bytes")) {
            1 => PixelType::Half,
            2 => PixelType::Float,
            other => {
                return Err(ImgError::Unsupported {
                    what: format!("pixel type {other}"),
                });
            }
        };
        pos += 16;
        channels.push(ExrInspectedChannel { name, ty });
    }
    Ok(channels)
}

fn collect_inspected_attributes<'a>(
    bytes: &'a [u8],
    count: usize,
    max_header_bytes: u64,
    poll: &mut impl FnMut() -> bool,
) -> Result<Vec<ExrInspectedAttribute<'a>>, ImgError> {
    let requested = u64::try_from(count)
        .unwrap_or(u64::MAX)
        .saturating_mul(size_of::<ExrInspectedAttribute<'static>>() as u64);
    let mut attributes = Vec::new();
    attributes
        .try_reserve_exact(count)
        .map_err(|_| ImgError::AllocationRefused {
            resource: "EXR inspection attribute descriptors",
            requested,
        })?;
    let mut pos = 8usize;
    for _ in 0..8 {
        inspect_continue(poll, "EXR structural inspection")?;
        let _ = read_attribute_ref(bytes, &mut pos, max_header_bytes)?;
    }
    while bytes.get(pos) != Some(&0) {
        inspect_continue(poll, "EXR structural inspection")?;
        let attribute = read_attribute_ref(bytes, &mut pos, max_header_bytes)?;
        attributes.push(ExrInspectedAttribute {
            name: attribute.name,
            ty: attribute.ty,
            value: attribute.value,
        });
    }
    Ok(attributes)
}

/// Strictly inspect an EXR emitted by this crate without decoding image
/// planes. This convenience entry point never observes cancellation.
///
/// # Errors
/// Returns [`ImgError`] for budget refusal, allocation refusal, malformed or
/// noncanonical structure, or an unsupported EXR feature.
pub fn inspect_exr(bytes: &[u8], limits: ExrInspectLimits) -> Result<ExrInspection<'_>, ImgError> {
    inspect_exr_with_poll(bytes, limits, || true)
}

/// Strictly inspect an EXR emitted by this crate without decoding image
/// planes, polling before bounded header/channel/scanline units. `poll`
/// returns `true` to continue and `false` to cancel.
///
/// The exact canonical writer subset is enforced: fixed built-in header
/// values and order, sorted channel/custom-attribute schemas, exact offset
/// table entries, increasing scanline indices, exact declared payload sizes,
/// and exact EOF.
///
/// # Errors
/// Returns [`ImgError`] for cancellation, budget/allocation refusal,
/// malformed or noncanonical structure, or an unsupported EXR feature.
#[allow(clippy::too_many_lines)] // one read-only transaction keeps all exact-offset checks visible
pub fn inspect_exr_with_poll(
    bytes: &[u8],
    limits: ExrInspectLimits,
    mut poll: impl FnMut() -> bool,
) -> Result<ExrInspection<'_>, ImgError> {
    inspect_continue(&mut poll, "EXR structural inspection")?;
    let input_bytes = u64::try_from(bytes.len()).map_err(|_| ImgError::SizeOverflow {
        context: "EXR input bytes",
    })?;
    admit_inspect_bytes("EXR input", input_bytes, limits.max_input_bytes)?;
    if take(bytes, 0, 4)? != MAGIC {
        return Err(ImgError::Malformed {
            what: "missing EXR magic".to_string(),
        });
    }
    let version = u32::from_le_bytes(take(bytes, 4, 4)?.try_into().expect("4 bytes"));
    if version != 2 {
        return Err(ImgError::Unsupported {
            what: format!("EXR version/flags {version:#x} (single-part v2 only)"),
        });
    }

    let header = inspect_strict_header(bytes, limits, &mut poll)?;
    inspect_continue(&mut poll, "EXR structural inspection")?;
    let channels =
        collect_inspected_channels(header.channels_value, header.channel_count, &mut poll)?;
    let attributes = collect_inspected_attributes(
        bytes,
        header.custom_attribute_count,
        limits.max_header_bytes,
        &mut poll,
    )?;

    let table_start = usize::try_from(header.header_bytes).map_err(|_| ImgError::SizeOverflow {
        context: "EXR offset-table start",
    })?;
    let table_bytes =
        usize::try_from(u64::from(header.height) * 8).map_err(|_| ImgError::SizeOverflow {
            context: "EXR offset-table bytes",
        })?;
    let mut block_pos = table_start
        .checked_add(table_bytes)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR first scanline offset",
        })?;
    let payload_bytes =
        usize::try_from(header.scanline_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR scanline bytes on this platform",
        })?;
    for y in 0..header.height {
        inspect_continue(&mut poll, "EXR structural inspection")?;
        let table_offset = usize::try_from(y)
            .map_err(|_| ImgError::SizeOverflow {
                context: "EXR offset-table index on this platform",
            })?
            .checked_mul(8)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR offset-table entry",
            })?;
        let table_pos = table_start
            .checked_add(table_offset)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR offset-table entry",
            })?;
        let observed_offset =
            u64::from_le_bytes(take(bytes, table_pos, 8)?.try_into().expect("8 bytes"));
        let expected_offset = u64::try_from(block_pos).map_err(|_| ImgError::SizeOverflow {
            context: "EXR scanline offset",
        })?;
        if observed_offset != expected_offset {
            return Err(ImgError::Malformed {
                what: format!(
                    "EXR offset table entry {y} is {observed_offset}; expected {expected_offset}"
                ),
            });
        }
        let observed_y =
            i32::from_le_bytes(take(bytes, block_pos, 4)?.try_into().expect("4 bytes"));
        let expected_y = i32::try_from(y).map_err(|_| ImgError::SizeOverflow {
            context: "EXR scanline index",
        })?;
        if observed_y != expected_y {
            return Err(ImgError::Malformed {
                what: format!("EXR scanline order broke at y={y}: found {observed_y}"),
            });
        }
        let declared_pos = block_pos.checked_add(4).ok_or(ImgError::SizeOverflow {
            context: "EXR scanline-size offset",
        })?;
        let declared_size =
            i32::from_le_bytes(take(bytes, declared_pos, 4)?.try_into().expect("4 bytes"));
        if declared_size < 0 || declared_size as u64 != header.scanline_bytes {
            return Err(ImgError::Malformed {
                what: format!(
                    "EXR scanline {y} declares {declared_size} bytes; expected {}",
                    header.scanline_bytes
                ),
            });
        }
        let payload_start = block_pos.checked_add(8).ok_or(ImgError::SizeOverflow {
            context: "EXR scanline payload offset",
        })?;
        let _ = take(bytes, payload_start, payload_bytes)?;
        block_pos = payload_start
            .checked_add(payload_bytes)
            .ok_or(ImgError::SizeOverflow {
                context: "EXR next scanline offset",
            })?;
    }
    if block_pos != bytes.len() {
        return Err(ImgError::Malformed {
            what: format!(
                "EXR parser stopped at byte {block_pos}, before exact EOF {}",
                bytes.len()
            ),
        });
    }

    Ok(ExrInspection {
        width: header.width,
        height: header.height,
        channels,
        attributes,
        input_bytes,
        header_bytes: header.header_bytes,
        scanline_bytes: header.scanline_bytes,
        decoded_bytes: header.decoded_bytes,
        metadata_bytes: header.metadata_bytes,
    })
}

const EXR_FLOAT_VERIFY_POLL_SAMPLES: usize = 4096;

/// Verify that every sample in one named `FLOAT` channel has exactly the
/// expected IEEE-754 bit pattern.
///
/// This is the non-cancellable convenience form of
/// [`verify_exr_float_channel_constant_with_poll`]. It first applies the same
/// strict writer-subset and caller-budget checks as [`inspect_exr`], then scans
/// only the requested channel in-place. No image plane is decoded or
/// allocated; the only successful-path allocations are the bounded borrowed
/// metadata descriptor vectors returned internally by the structural
/// inspector.
///
/// # Errors
/// Returns [`ImgError`] when the artifact is malformed, outside the strict
/// subset, over budget, lacks the requested channel, declares that channel as
/// a non-`FLOAT` type, or contains a sample whose bits differ from `expected`.
pub fn verify_exr_float_channel_constant(
    bytes: &[u8],
    channel_name: &str,
    expected: f32,
    limits: ExrInspectLimits,
) -> Result<(), ImgError> {
    verify_exr_float_channel_constant_with_poll(bytes, channel_name, expected, limits, || true)
}

/// Verify exact constant `FLOAT` channel samples with bounded cancellation
/// polling and without materializing image planes.
///
/// Structural validation polls at the units documented by
/// [`inspect_exr_with_poll`]. The sample scan polls at every scanline boundary
/// and after at most 4,096 samples within a scanline. A mismatch reports the
/// deterministic first `(x, y)` coordinate in row-major order and compares
/// [`f32::to_bits`] values, so signed zero and NaN payload differences are not
/// hidden by floating-point equality.
///
/// # Errors
/// Returns the same structural, budget, allocation, and cancellation errors as
/// [`inspect_exr_with_poll`]. A missing required channel or mismatching sample
/// is [`ImgError::Malformed`]; a present non-`FLOAT` channel is
/// [`ImgError::Unsupported`].
#[allow(clippy::too_many_lines)] // checked byte-offset arithmetic remains explicit and auditable
pub fn verify_exr_float_channel_constant_with_poll(
    bytes: &[u8],
    channel_name: &str,
    expected: f32,
    limits: ExrInspectLimits,
    mut poll: impl FnMut() -> bool,
) -> Result<(), ImgError> {
    let inspection = inspect_exr_with_poll(bytes, limits, &mut poll)?;
    let width = usize::try_from(inspection.width).map_err(|_| ImgError::SizeOverflow {
        context: "EXR FLOAT verification width on this platform",
    })?;
    let height = inspection.height;
    let header_bytes =
        usize::try_from(inspection.header_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR FLOAT verification header offset",
        })?;
    let scanline_bytes =
        usize::try_from(inspection.scanline_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR FLOAT verification scanline bytes on this platform",
        })?;

    let mut channel_row_offset = 0usize;
    let mut requested_type = None;
    for channel in &inspection.channels {
        if channel.name == channel_name {
            requested_type = Some(channel.ty);
            break;
        }
        let channel_row_bytes =
            width
                .checked_mul(channel.ty.bytes())
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR FLOAT verification channel-row bytes",
                })?;
        channel_row_offset =
            channel_row_offset
                .checked_add(channel_row_bytes)
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR FLOAT verification channel offset",
                })?;
    }
    let Some(requested_type) = requested_type else {
        return Err(ImgError::Malformed {
            what: format!("required EXR FLOAT channel {channel_name:?} is missing"),
        });
    };
    if requested_type != PixelType::Float {
        return Err(ImgError::Unsupported {
            what: format!(
                "required EXR channel {channel_name:?} is {requested_type:?}; FLOAT is required"
            ),
        });
    }

    let table_bytes = usize::try_from(u64::from(height).checked_mul(8).ok_or(
        ImgError::SizeOverflow {
            context: "EXR FLOAT verification offset table",
        },
    )?)
    .map_err(|_| ImgError::SizeOverflow {
        context: "EXR FLOAT verification offset table on this platform",
    })?;
    let scanline_block_bytes = scanline_bytes
        .checked_add(8)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR FLOAT verification scanline block",
        })?;
    let first_block = header_bytes
        .checked_add(table_bytes)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR FLOAT verification first scanline",
        })?;
    let expected_bits = expected.to_bits();

    for y in 0..height {
        let y_usize = usize::try_from(y).map_err(|_| ImgError::SizeOverflow {
            context: "EXR FLOAT verification row index on this platform",
        })?;
        let block_offset =
            y_usize
                .checked_mul(scanline_block_bytes)
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR FLOAT verification scanline offset",
                })?;
        let row_start = first_block
            .checked_add(block_offset)
            .and_then(|offset| offset.checked_add(8))
            .and_then(|offset| offset.checked_add(channel_row_offset))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR FLOAT verification channel row",
            })?;

        for x in 0..width {
            if x % EXR_FLOAT_VERIFY_POLL_SAMPLES == 0 {
                inspect_continue(&mut poll, "EXR FLOAT channel constant verification")?;
            }
            let sample_offset = x
                .checked_mul(size_of::<f32>())
                .and_then(|offset| row_start.checked_add(offset))
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR FLOAT verification sample offset",
                })?;
            let observed_bits =
                u32::from_le_bytes(take(bytes, sample_offset, 4)?.try_into().expect("4 bytes"));
            if observed_bits != expected_bits {
                return Err(ImgError::Malformed {
                    what: format!(
                        "EXR FLOAT channel {channel_name:?} sample at ({x}, {y}) has bits \
                         {observed_bits:#010x}; expected {expected_bits:#010x}"
                    ),
                });
            }
        }
    }
    Ok(())
}

const MAX_CONSECUTIVE_F32_INTEGER: f32 = 16_777_216.0;
const RAW_VALIDITY_PRIMARY: u32 = 1 << 0;
const RAW_VALIDITY_OBJECT_ID: u32 = 1 << 4;
const RAW_VALIDITY_MATERIAL_ID: u32 = 1 << 5;
const RAW_VALIDITY_CONTRIBUTION_SPLIT: u32 = 1 << 6;

/// Caller-owned semantic limits for a raw cinematic AOV payload.
///
/// Palette indices are inclusive one-based maxima; zero remains the canonical
/// background/unavailable value. Both maxima must be no larger than 2^24, the
/// largest consecutive integer representable by an EXR `FLOAT` sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExrRawFrameSemanticLimits {
    /// Complete validity-bit domain admitted by the bound render semantics.
    allowed_validity_bits: u32,
    /// Largest admitted nonzero `id.object` palette index.
    maximum_object_palette_index: u32,
    /// Largest admitted nonzero `id.material` palette index.
    maximum_material_palette_index: u32,
}

impl ExrRawFrameSemanticLimits {
    /// Admit an exact validity domain and inclusive one-based palette maxima.
    ///
    /// # Errors
    /// Returns [`ImgError::Unsupported`] when either palette maximum is above
    /// the largest consecutive integer representable by an EXR `FLOAT`.
    pub fn try_new(
        allowed_validity_bits: u32,
        maximum_object_palette_index: u32,
        maximum_material_palette_index: u32,
    ) -> Result<Self, ImgError> {
        for (kind, maximum) in [
            ("object", maximum_object_palette_index),
            ("material", maximum_material_palette_index),
        ] {
            if maximum > MAX_CONSECUTIVE_F32_INTEGER as u32 {
                return Err(ImgError::Unsupported {
                    what: format!(
                        "raw {kind} palette maximum {maximum} exceeds the consecutive-EXR-FLOAT \
                         index limit {}",
                        MAX_CONSECUTIVE_F32_INTEGER as u32
                    ),
                });
            }
        }
        Ok(Self {
            allowed_validity_bits,
            maximum_object_palette_index,
            maximum_material_palette_index,
        })
    }

    /// Preserve the historical payload validator's full binary32 ID range.
    #[must_use]
    pub const fn full_f32_index_range(allowed_validity_bits: u32) -> Self {
        Self {
            allowed_validity_bits,
            maximum_object_palette_index: MAX_CONSECUTIVE_F32_INTEGER as u32,
            maximum_material_palette_index: MAX_CONSECUTIVE_F32_INTEGER as u32,
        }
    }

    /// Complete admitted validity-bit domain.
    #[must_use]
    pub const fn allowed_validity_bits(self) -> u32 {
        self.allowed_validity_bits
    }

    /// Inclusive largest admitted nonzero `id.object` index.
    #[must_use]
    pub const fn maximum_object_palette_index(self) -> u32 {
        self.maximum_object_palette_index
    }

    /// Inclusive largest admitted nonzero `id.material` index.
    #[must_use]
    pub const fn maximum_material_palette_index(self) -> u32 {
        self.maximum_material_palette_index
    }
}

fn is_raw_semantic_channel(name: &str) -> bool {
    matches!(
        name,
        "id.object"
            | "id.material"
            | "samples"
            | "diagnostic.validity"
            | "primary.coverage"
            | "variance.Y"
    )
}

fn exact_nonnegative_f32_integer(value: f32) -> Option<u32> {
    if (0.0..=MAX_CONSECUTIVE_F32_INTEGER).contains(&value) && value.fract() == 0.0 {
        Some(value as u32)
    } else {
        None
    }
}

fn raw_float_semantic_violation(
    channel_name: &str,
    value: f32,
    semantic_limits: ExrRawFrameSemanticLimits,
) -> Option<String> {
    if !value.is_finite() {
        return Some("must be finite".to_string());
    }
    match channel_name {
        "id.object" | "id.material" => match exact_nonnegative_f32_integer(value) {
            None => Some(format!(
                "must be a nonnegative exact integer in the consecutive-f32 range \
                 [0, {}]",
                MAX_CONSECUTIVE_F32_INTEGER as u32
            )),
            Some(index) => {
                let maximum = if channel_name == "id.object" {
                    semantic_limits.maximum_object_palette_index
                } else {
                    semantic_limits.maximum_material_palette_index
                };
                (index > maximum).then(|| {
                    format!(
                        "palette index {index} exceeds the admitted one-based maximum {maximum}"
                    )
                })
            }
        },
        "samples" => match exact_nonnegative_f32_integer(value) {
            Some(integer) if integer > 0 => None,
            _ => Some(format!(
                "must be a positive exact integer in the consecutive-f32 range \
                 [1, {}]",
                MAX_CONSECUTIVE_F32_INTEGER as u32
            )),
        },
        "diagnostic.validity" => match exact_nonnegative_f32_integer(value) {
            Some(bits) if bits & !semantic_limits.allowed_validity_bits == 0 => None,
            Some(bits) => Some(format!(
                "contains unknown bits {:#010x}; caller allows {:#010x}",
                bits & !semantic_limits.allowed_validity_bits,
                semantic_limits.allowed_validity_bits,
            )),
            None => Some(format!(
                "must be a nonnegative exact integer bit mask in the consecutive-f32 range \
                 [0, {}]",
                MAX_CONSECUTIVE_F32_INTEGER as u32
            )),
        },
        "primary.coverage" if !(0.0..=1.0).contains(&value) => {
            Some("must lie in the closed interval [0, 1]".to_string())
        }
        "variance.Y" if value < 0.0 => Some("must be nonnegative".to_string()),
        _ => None,
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct RawSemanticChannelOffsets {
    validity: Option<usize>,
    object_id: Option<usize>,
    material_id: Option<usize>,
    primary_coverage: Option<usize>,
    samples: Option<usize>,
}

impl RawSemanticChannelOffsets {
    fn record(&mut self, name: &str, row_offset: usize) {
        let slot = match name {
            "diagnostic.validity" => &mut self.validity,
            "id.object" => &mut self.object_id,
            "id.material" => &mut self.material_id,
            "primary.coverage" => &mut self.primary_coverage,
            "samples" => &mut self.samples,
            _ => return,
        };
        *slot = Some(row_offset);
    }

    const fn correlated(self) -> Option<[usize; 5]> {
        match (
            self.validity,
            self.object_id,
            self.material_id,
            self.primary_coverage,
            self.samples,
        ) {
            (Some(validity), Some(object), Some(material), Some(coverage), Some(samples)) => {
                Some([validity, object, material, coverage, samples])
            }
            _ => None,
        }
    }
}

/// Validate every `FLOAT` sample in a strict raw-frame EXR without decoding
/// image planes.
///
/// This convenience form never observes cancellation. It applies the same
/// strict structure and caller-owned byte ceilings as [`inspect_exr`], then
/// scans every `FLOAT` channel in canonical wire order. Every value must be
/// finite. When present, the canonical raw-AOV semantic channels receive
/// additional checks:
///
/// - `id.object` and `id.material` are nonnegative exact integers no larger
///   than 2^24, the consecutive-integer limit of binary32;
/// - `samples` is a positive exact integer under the same limit;
/// - `diagnostic.validity` is such an integer mask and contains no bit outside
///   `allowed_validity_bits`;
/// - `primary.coverage` lies in `[0, 1]`; and
/// - `variance.Y` is nonnegative.
///
/// When `diagnostic.validity` is present, the final-diagnostic per-pixel
/// relations documented by [`validate_exr_raw_frame_payload_against`] are also
/// enforced. This compatibility entry point admits the full exactly
/// representable binary32 palette-index range; callers with an exact palette
/// should use the typed variant instead.
///
/// Other `FLOAT` channels have only the universal finite-value requirement.
/// Non-`FLOAT` channels are structurally validated but otherwise skipped;
/// canonical semantic channels are rejected if declared with another type.
///
/// # Errors
/// Returns [`ImgError`] for structure, budget, allocation, or semantic
/// failures. Semantic errors identify the deterministic first channel and
/// `(x, y)` sample in scanline/channel/x wire order.
pub fn validate_exr_raw_frame_payload(
    bytes: &[u8],
    allowed_validity_bits: u32,
    limits: ExrInspectLimits,
) -> Result<(), ImgError> {
    validate_exr_raw_frame_payload_against(
        bytes,
        ExrRawFrameSemanticLimits::full_f32_index_range(allowed_validity_bits),
        limits,
    )
}

/// Cancellable form of [`validate_exr_raw_frame_payload`].
///
/// The payload scan polls at every `FLOAT` channel row and after at most 4,096
/// samples within that row. Together with [`inspect_exr_with_poll`], this
/// bounds cancellation latency without allocating any image planes.
///
/// # Errors
/// Returns the errors documented by [`validate_exr_raw_frame_payload`], plus
/// [`ImgError::Cancelled`] when `poll` returns `false`.
#[allow(clippy::too_many_lines)] // explicit checked wire offsets keep the validator auditable
pub fn validate_exr_raw_frame_payload_with_poll(
    bytes: &[u8],
    allowed_validity_bits: u32,
    limits: ExrInspectLimits,
    poll: impl FnMut() -> bool,
) -> Result<(), ImgError> {
    validate_exr_raw_frame_payload_against_with_poll(
        bytes,
        ExrRawFrameSemanticLimits::full_f32_index_range(allowed_validity_bits),
        limits,
        poll,
    )
}

/// Validate a strict raw-frame EXR against exact palette and validity
/// semantics without allocating image planes.
///
/// In addition to [`validate_exr_raw_frame_payload`], this entry point rejects
/// object or material IDs above the caller's exact one-based palette maxima.
/// When `diagnostic.validity` is present, the final-diagnostic correlation
/// channels must also be present and every pixel must satisfy all four frozen
/// relations:
///
/// - `OBJECT_ID` is set exactly when `id.object != 0`;
/// - `MATERIAL_ID` is set exactly when `id.material != 0`;
/// - `PRIMARY` is set exactly when `primary.coverage > 0`; and
/// - `CONTRIBUTION_SPLIT` is set exactly when `samples > 0`.
///
/// A nonzero object or material ID also requires positive primary coverage;
/// the converse is intentionally not required because a primary can lack a
/// stable categorical identity.
///
/// Profiles without `diagnostic.validity` retain their profile-local checks;
/// for example, daily-core coverage is range-checked without inventing a
/// final-diagnostic validity plane.
///
/// # Errors
/// Returns [`ImgError`] for invalid semantic limits, structure, budget,
/// allocation, or payload-semantic failures.
pub fn validate_exr_raw_frame_payload_against(
    bytes: &[u8],
    semantic_limits: ExrRawFrameSemanticLimits,
    limits: ExrInspectLimits,
) -> Result<(), ImgError> {
    validate_exr_raw_frame_payload_against_with_poll(bytes, semantic_limits, limits, || true)
}

/// Cancellable form of [`validate_exr_raw_frame_payload_against`].
///
/// Structural and channel-local scans retain their existing poll bounds. The
/// correlation pass additionally polls before every row and after at most
/// 4,096 pixels within a row.
///
/// # Errors
/// Returns the errors documented by
/// [`validate_exr_raw_frame_payload_against`], plus [`ImgError::Cancelled`]
/// when `poll` returns `false`.
#[allow(clippy::too_many_lines)] // explicit checked wire offsets keep the validator auditable
pub fn validate_exr_raw_frame_payload_against_with_poll(
    bytes: &[u8],
    semantic_limits: ExrRawFrameSemanticLimits,
    limits: ExrInspectLimits,
    mut poll: impl FnMut() -> bool,
) -> Result<(), ImgError> {
    let inspection = inspect_exr_with_poll(bytes, limits, &mut poll)?;
    let width = usize::try_from(inspection.width).map_err(|_| ImgError::SizeOverflow {
        context: "EXR raw-frame validation width on this platform",
    })?;
    let mut semantic_offsets = RawSemanticChannelOffsets::default();
    let mut channel_row_offset = 0usize;
    for channel in &inspection.channels {
        inspect_continue(&mut poll, "EXR raw-frame payload validation")?;
        if is_raw_semantic_channel(channel.name) && channel.ty != PixelType::Float {
            return Err(ImgError::Unsupported {
                what: format!(
                    "raw-frame semantic channel {:?} is {:?}; FLOAT is required",
                    channel.name, channel.ty
                ),
            });
        }
        semantic_offsets.record(channel.name, channel_row_offset);
        channel_row_offset =
            channel_row_offset
                .checked_add(width.checked_mul(channel.ty.bytes()).ok_or(
                    ImgError::SizeOverflow {
                        context: "EXR raw-frame validation channel-row bytes",
                    },
                )?)
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR raw-frame validation channel offset",
                })?;
    }
    if semantic_offsets.validity.is_some() && semantic_offsets.correlated().is_none() {
        return Err(ImgError::Malformed {
            what: "EXR diagnostic.validity requires id.object, id.material, \
                   primary.coverage, and samples correlation channels"
                .to_string(),
        });
    }
    let height = inspection.height;
    let header_bytes =
        usize::try_from(inspection.header_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR raw-frame validation header offset",
        })?;
    let scanline_bytes =
        usize::try_from(inspection.scanline_bytes).map_err(|_| ImgError::SizeOverflow {
            context: "EXR raw-frame validation scanline bytes on this platform",
        })?;
    let table_bytes = usize::try_from(u64::from(height).checked_mul(8).ok_or(
        ImgError::SizeOverflow {
            context: "EXR raw-frame validation offset table",
        },
    )?)
    .map_err(|_| ImgError::SizeOverflow {
        context: "EXR raw-frame validation offset table on this platform",
    })?;
    let scanline_block_bytes = scanline_bytes
        .checked_add(8)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR raw-frame validation scanline block",
        })?;
    let first_block = header_bytes
        .checked_add(table_bytes)
        .ok_or(ImgError::SizeOverflow {
            context: "EXR raw-frame validation first scanline",
        })?;

    for y in 0..height {
        let y_usize = usize::try_from(y).map_err(|_| ImgError::SizeOverflow {
            context: "EXR raw-frame validation row index on this platform",
        })?;
        let block_offset =
            y_usize
                .checked_mul(scanline_block_bytes)
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR raw-frame validation scanline offset",
                })?;
        let payload_start = first_block
            .checked_add(block_offset)
            .and_then(|offset| offset.checked_add(8))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR raw-frame validation payload offset",
            })?;
        let mut channel_row_offset = 0usize;
        for channel in &inspection.channels {
            let channel_row_bytes =
                width
                    .checked_mul(channel.ty.bytes())
                    .ok_or(ImgError::SizeOverflow {
                        context: "EXR raw-frame validation channel-row bytes",
                    })?;
            if channel.ty == PixelType::Float {
                let row_start = payload_start.checked_add(channel_row_offset).ok_or(
                    ImgError::SizeOverflow {
                        context: "EXR raw-frame validation channel row",
                    },
                )?;
                for x in 0..width {
                    if x % EXR_FLOAT_VERIFY_POLL_SAMPLES == 0 {
                        inspect_continue(&mut poll, "EXR raw-frame payload validation")?;
                    }
                    let sample_offset = x
                        .checked_mul(size_of::<f32>())
                        .and_then(|offset| row_start.checked_add(offset))
                        .ok_or(ImgError::SizeOverflow {
                            context: "EXR raw-frame validation sample offset",
                        })?;
                    let observed_bits = u32::from_le_bytes(
                        take(bytes, sample_offset, 4)?.try_into().expect("4 bytes"),
                    );
                    let observed = f32::from_bits(observed_bits);
                    if let Some(reason) =
                        raw_float_semantic_violation(channel.name, observed, semantic_limits)
                    {
                        return Err(ImgError::Malformed {
                            what: format!(
                                "EXR raw FLOAT channel {:?} sample at ({x}, {y}) has bits \
                                 {observed_bits:#010x}: {reason}",
                                channel.name
                            ),
                        });
                    }
                }
            }
            channel_row_offset = channel_row_offset.checked_add(channel_row_bytes).ok_or(
                ImgError::SizeOverflow {
                    context: "EXR raw-frame validation channel offset",
                },
            )?;
        }
    }
    if let Some([validity, object, material, coverage, samples]) = semantic_offsets.correlated() {
        verify_raw_semantic_relations(
            bytes,
            width,
            height,
            first_block,
            scanline_block_bytes,
            [validity, object, material, coverage, samples],
            &mut poll,
        )?;
    }
    Ok(())
}

fn verify_raw_semantic_relations(
    bytes: &[u8],
    width: usize,
    height: u32,
    first_block: usize,
    scanline_block_bytes: usize,
    channel_offsets: [usize; 5],
    poll: &mut impl FnMut() -> bool,
) -> Result<(), ImgError> {
    let [
        validity_offset,
        object_offset,
        material_offset,
        coverage_offset,
        samples_offset,
    ] = channel_offsets;
    for y in 0..height {
        let y_usize = usize::try_from(y).map_err(|_| ImgError::SizeOverflow {
            context: "EXR raw semantic relation row index on this platform",
        })?;
        let payload_start = y_usize
            .checked_mul(scanline_block_bytes)
            .and_then(|offset| first_block.checked_add(offset))
            .and_then(|offset| offset.checked_add(8))
            .ok_or(ImgError::SizeOverflow {
                context: "EXR raw semantic relation row",
            })?;
        for x in 0..width {
            if x % EXR_FLOAT_VERIFY_POLL_SAMPLES == 0 {
                inspect_continue(poll, "EXR raw-frame semantic relation validation")?;
            }
            let sample_offset = x
                .checked_mul(size_of::<f32>())
                .ok_or(ImgError::SizeOverflow {
                    context: "EXR raw semantic relation sample offset",
                })?;
            let read = |channel_offset: usize| -> Result<f32, ImgError> {
                let offset = payload_start
                    .checked_add(channel_offset)
                    .and_then(|offset| offset.checked_add(sample_offset))
                    .ok_or(ImgError::SizeOverflow {
                        context: "EXR raw semantic relation channel sample",
                    })?;
                Ok(f32::from_bits(u32::from_le_bytes(
                    take(bytes, offset, 4)?.try_into().expect("4 bytes"),
                )))
            };
            let validity =
                exact_nonnegative_f32_integer(read(validity_offset)?).ok_or_else(|| {
                    ImgError::Malformed {
                        what: format!(
                            "EXR diagnostic.validity sample at ({x}, {y}) changed during validation"
                        ),
                    }
                })?;
            let object = exact_nonnegative_f32_integer(read(object_offset)?).ok_or_else(|| {
                ImgError::Malformed {
                    what: format!("EXR id.object sample at ({x}, {y}) changed during validation"),
                }
            })?;
            let material =
                exact_nonnegative_f32_integer(read(material_offset)?).ok_or_else(|| {
                    ImgError::Malformed {
                        what: format!(
                            "EXR id.material sample at ({x}, {y}) changed during validation"
                        ),
                    }
                })?;
            let coverage = read(coverage_offset)?;
            let samples =
                exact_nonnegative_f32_integer(read(samples_offset)?).ok_or_else(|| {
                    ImgError::Malformed {
                        what: format!("EXR samples value at ({x}, {y}) changed during validation"),
                    }
                })?;

            let relation = [
                (
                    "OBJECT_ID iff id.object != 0",
                    validity & RAW_VALIDITY_OBJECT_ID != 0,
                    object != 0,
                ),
                (
                    "MATERIAL_ID iff id.material != 0",
                    validity & RAW_VALIDITY_MATERIAL_ID != 0,
                    material != 0,
                ),
                (
                    "PRIMARY iff primary.coverage > 0",
                    validity & RAW_VALIDITY_PRIMARY != 0,
                    coverage > 0.0,
                ),
                (
                    "CONTRIBUTION_SPLIT iff samples > 0",
                    validity & RAW_VALIDITY_CONTRIBUTION_SPLIT != 0,
                    samples > 0,
                ),
            ]
            .into_iter()
            .find(|(_, bit_is_set, predicate_is_true)| bit_is_set != predicate_is_true);
            if let Some((relation, bit_is_set, predicate_is_true)) = relation {
                return Err(ImgError::Malformed {
                    what: format!(
                        "EXR raw semantic relation {relation} failed at ({x}, {y}): \
                         validity bit is {bit_is_set}, payload predicate is {predicate_is_true}"
                    ),
                });
            }
            if (object != 0 || material != 0) && coverage <= 0.0 {
                return Err(ImgError::Malformed {
                    what: format!(
                        "EXR raw semantic relation nonzero object/material ID requires \
                         primary.coverage > 0 failed at ({x}, {y})"
                    ),
                });
            }
        }
    }
    inspect_continue(poll, "EXR raw-frame semantic relation validation")
}

fn take(bytes: &[u8], pos: usize, n: usize) -> Result<&[u8], ImgError> {
    let end = pos.checked_add(n).ok_or_else(|| ImgError::Malformed {
        what: format!("byte range overflow at byte {pos}"),
    })?;
    bytes.get(pos..end).ok_or_else(|| ImgError::Malformed {
        what: format!("truncated at byte {pos}"),
    })
}

fn read_cstr(bytes: &[u8], pos: &mut usize) -> Result<String, ImgError> {
    let start = *pos;
    while *pos < bytes.len() && bytes[*pos] != 0 {
        *pos += 1;
    }
    if *pos >= bytes.len() {
        return Err(ImgError::Malformed {
            what: "unterminated string".to_string(),
        });
    }
    let s = std::str::from_utf8(&bytes[start..*pos])
        .map_err(|_| ImgError::Malformed {
            what: format!("non-UTF-8 header string at byte {start}"),
        })?
        .to_owned();
    *pos += 1;
    Ok(s)
}

fn parse_chlist(value: &[u8]) -> Result<Vec<(String, PixelType)>, ImgError> {
    let mut specs = Vec::new();
    let mut cp = 0usize;
    while cp < value.len() && value[cp] != 0 {
        let start = cp;
        while cp < value.len() && value[cp] != 0 {
            cp += 1;
        }
        let cname = String::from_utf8_lossy(&value[start..cp]).into_owned();
        cp += 1;
        let code = u32::from_le_bytes(
            value
                .get(cp..cp + 4)
                .ok_or_else(|| ImgError::Malformed {
                    what: "chlist truncated".to_string(),
                })?
                .try_into()
                .expect("4 bytes"),
        );
        let ty = match code {
            1 => PixelType::Half,
            2 => PixelType::Float,
            other => {
                return Err(ImgError::Unsupported {
                    what: format!("pixel type {other}"),
                });
            }
        };
        cp += 16; // type + pLinear/reserved + samplings
        specs.push((cname, ty));
    }
    Ok(specs)
}

/// Parsed header: (channel specs, width, height, custom attributes).
type HeaderInfo = (Vec<(String, PixelType)>, u32, u32, Vec<ExrAttribute>);

/// Parse the header attributes; returns [`HeaderInfo`] and leaves `pos`
/// just past the header terminator.
fn parse_header(bytes: &[u8], pos: &mut usize) -> Result<HeaderInfo, ImgError> {
    let mut specs: Vec<(String, PixelType)> = Vec::new();
    let mut window = (0u32, 0u32);
    let mut compression_seen = false;
    let mut attributes = Vec::new();
    let mut seen_names = BTreeSet::new();
    loop {
        if bytes.get(*pos) == Some(&0) {
            *pos += 1;
            break; // end of header
        }
        let name = read_cstr(bytes, pos)?;
        if !seen_names.insert(name.clone()) {
            return Err(ImgError::Malformed {
                what: format!("duplicate header attribute {name:?}"),
            });
        }
        let ty = read_cstr(bytes, pos)?;
        let size = u32::from_le_bytes(take(bytes, *pos, 4)?.try_into().expect("4 bytes")) as usize;
        *pos += 4;
        let value = take(bytes, *pos, size)?.to_vec();
        *pos += size;
        match name.as_str() {
            "channels" => specs = parse_chlist(&value)?,
            "compression" => {
                compression_seen = true;
                match value.as_slice() {
                    [0] => {}
                    [] => {
                        return Err(ImgError::Malformed {
                            what: "empty compression attribute".to_string(),
                        });
                    }
                    [code, ..] => {
                        return Err(ImgError::Unsupported {
                            what: format!("compression {code} (NONE only)"),
                        });
                    }
                }
            }
            "dataWindow" => {
                if value.len() != 16 {
                    return Err(ImgError::Malformed {
                        what: "box2i size".to_string(),
                    });
                }
                let x2 = i32::from_le_bytes(value[8..12].try_into().expect("4"));
                let y2 = i32::from_le_bytes(value[12..16].try_into().expect("4"));
                if x2 < 0 || y2 < 0 || x2 == i32::MAX || y2 == i32::MAX {
                    return Err(ImgError::Malformed {
                        what: "negative dataWindow extent".to_string(),
                    });
                }
                window = ((x2 + 1).cast_unsigned(), (y2 + 1).cast_unsigned());
            }
            "displayWindow" | "lineOrder" | "pixelAspectRatio" | "screenWindowCenter"
            | "screenWindowWidth" => {}
            _ => attributes.push(ExrAttribute { name, ty, value }),
        }
    }
    if !compression_seen || specs.is_empty() || window.0 == 0 || window.1 == 0 {
        return Err(ImgError::Malformed {
            what: "missing required header attributes".to_string(),
        });
    }
    Ok((specs, window.0, window.1, attributes))
}

/// Decode an EXR produced by [`write_exr`] or
/// [`write_exr_with_attributes`]. Structured rejection outside that subset.
///
/// # Errors
/// [`ImgError::Malformed`] / [`ImgError::Unsupported`].
pub fn read_exr(bytes: &[u8]) -> Result<DecodedExr, ImgError> {
    if take(bytes, 0, 4)? != MAGIC {
        return Err(ImgError::Malformed {
            what: "missing EXR magic".to_string(),
        });
    }
    let version = u32::from_le_bytes(take(bytes, 4, 4)?.try_into().expect("4 bytes"));
    if version != 2 {
        return Err(ImgError::Unsupported {
            what: format!("EXR version/flags {version:#x} (single-part v2 only)"),
        });
    }
    let mut pos = 8usize;
    let (specs, width, height, attributes) = parse_header(bytes, &mut pos)?;
    // Skip the offset table; read blocks sequentially (our writer's order).
    pos += 8 * height as usize;
    let n = width as usize * height as usize;
    let mut channels: Vec<Channel> = specs
        .iter()
        .map(|(name, ty)| Channel {
            name: name.clone(),
            ty: *ty,
            data: vec![0.0; n],
        })
        .collect();
    for y in 0..height as usize {
        let block_y = usize::try_from(i32::from_le_bytes(
            take(bytes, pos, 4)?.try_into().expect("4 bytes"),
        ))
        .map_err(|_| ImgError::Malformed {
            what: "negative scanline y".to_string(),
        })?;
        pos += 8; // y + declared size
        if block_y != y {
            return Err(ImgError::Malformed {
                what: format!("scanline order broke at y={y}"),
            });
        }
        for c in &mut channels {
            for x in 0..width as usize {
                let v = match c.ty {
                    PixelType::Half => {
                        let b = take(bytes, pos, 2)?;
                        pos += 2;
                        f16_bits_to_f32(u16::from_le_bytes([b[0], b[1]]))
                    }
                    PixelType::Float => {
                        let b = take(bytes, pos, 4)?;
                        pos += 4;
                        f32::from_le_bytes(b.try_into().expect("4 bytes"))
                    }
                };
                c.data[y * width as usize + x] = v;
            }
        }
    }
    Ok(DecodedExr {
        width,
        height,
        channels,
        attributes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f16_conversion_known_answers_and_round_trip() {
        assert_eq!(f32_to_f16_bits(0.0), 0x0000);
        assert_eq!(f32_to_f16_bits(-0.0), 0x8000);
        assert_eq!(f32_to_f16_bits(1.0), 0x3C00);
        assert_eq!(f32_to_f16_bits(-2.0), 0xC000);
        assert_eq!(f32_to_f16_bits(65504.0), 0x7BFF); // max half
        assert_eq!(f32_to_f16_bits(1e6), 0x7C00); // overflow → inf
        assert_eq!(f32_to_f16_bits(f32::INFINITY), 0x7C00);
        assert!(f16_bits_to_f32(f32_to_f16_bits(f32::NAN)).is_nan());
        // Smallest subnormal half.
        assert_eq!(f32_to_f16_bits(5.960_464_5e-8), 0x0001);
        // Every finite half survives f16 → f32 → f16 exactly.
        for h in 0..=0x7BFFu16 {
            let back = f32_to_f16_bits(f16_bits_to_f32(h));
            assert_eq!(back, h, "half round-trip broke at {h:#06x}");
        }
    }

    #[test]
    fn exr_aov_round_trip_is_lossless() {
        let (w, h) = (6u32, 4u32);
        let n = (w * h) as usize;
        let ch = |name: &str, ty: PixelType, k: f32| Channel {
            name: name.to_string(),
            ty,
            data: (0..n).map(|i| (i as f32) * k - 3.0).collect(),
        };
        let channels = vec![
            ch("R", PixelType::Float, 0.25),
            ch("G", PixelType::Float, 0.5),
            ch("B", PixelType::Float, 0.75),
            ch("albedo.R", PixelType::Half, 0.03125),
            ch("depth.Z", PixelType::Float, 1.5),
        ];
        let bytes = write_exr(w, h, &channels).unwrap();
        assert_eq!(
            bytes,
            write_exr(w, h, &channels).unwrap(),
            "byte-exact determinism"
        );
        let decoded = read_exr(&bytes).unwrap();
        assert_eq!((decoded.width, decoded.height), (w, h));
        // Alphabetical order per spec.
        let names: Vec<&str> = decoded.channels.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, vec!["B", "G", "R", "albedo.R", "depth.Z"]);
        // FLOAT channels round-trip exactly; HALF values chosen on the
        // half grid (multiples of 2⁻⁵) round-trip exactly too.
        for c in &decoded.channels {
            let orig = channels.iter().find(|o| o.name == c.name).expect("name");
            assert_eq!(c.data, orig.data, "channel {} drifted", c.name);
        }
        assert!(decoded.attributes.is_empty());
    }

    #[test]
    fn exr_custom_attributes_round_trip_and_empty_path_stays_exact() {
        let channel = Channel {
            name: "R".to_string(),
            ty: PixelType::Float,
            data: vec![0.25, 0.5],
        };
        let empty_legacy = write_exr(2, 1, std::slice::from_ref(&channel)).unwrap();
        let empty_explicit =
            write_exr_with_attributes(2, 1, std::slice::from_ref(&channel), &[]).unwrap();
        assert_eq!(empty_legacy, empty_explicit);

        let attributes = vec![
            ExrAttribute {
                name: SOURCE_ARTIFACT_HASH_ATTRIBUTE.to_string(),
                ty: "string".to_string(),
                value: b"blake3:0123456789abcdef".to_vec(),
            },
            ExrAttribute {
                name: "frankensim.binaryReceipt".to_string(),
                ty: "opaque".to_string(),
                value: vec![0, 0xFF, 0, 7],
            },
            ExrAttribute {
                name: "frankensim.emptyReceipt".to_string(),
                ty: "opaque".to_string(),
                value: Vec::new(),
            },
        ];
        let encoded =
            write_exr_with_attributes(2, 1, std::slice::from_ref(&channel), &attributes).unwrap();
        assert_eq!(
            encoded,
            write_exr_with_attributes(2, 1, std::slice::from_ref(&channel), &attributes).unwrap()
        );
        let decoded = read_exr(&encoded).unwrap();
        assert_eq!(
            decoded
                .attributes
                .iter()
                .map(|attribute| attribute.name.as_str())
                .collect::<Vec<_>>(),
            vec![
                "frankensim.binaryReceipt",
                "frankensim.emptyReceipt",
                SOURCE_ARTIFACT_HASH_ATTRIBUTE,
            ]
        );
        assert_eq!(decoded.attributes[0].value.as_slice(), &[0, 0xFF, 0, 7]);
        assert!(decoded.attributes[1].value.is_empty());
        assert_eq!(
            decoded.attributes[2].value.as_slice(),
            b"blake3:0123456789abcdef"
        );
        assert_eq!(
            encoded,
            write_exr_with_attributes(
                decoded.width,
                decoded.height,
                &decoded.channels,
                &decoded.attributes,
            )
            .unwrap(),
            "our metadata-bearing subset re-encodes byte exactly"
        );
    }

    fn inspected_fixture() -> Vec<u8> {
        let channels = vec![
            Channel {
                name: "depth.Z".to_string(),
                ty: PixelType::Half,
                data: vec![1.0; 12],
            },
            Channel {
                name: "R".to_string(),
                ty: PixelType::Float,
                data: vec![0.25; 12],
            },
        ];
        let attributes = vec![
            ExrAttribute {
                name: "frankensim.zeta".to_string(),
                ty: "opaque".to_string(),
                value: vec![0, 0xFF, 7],
            },
            ExrAttribute {
                name: "frankensim.alpha".to_string(),
                ty: "string".to_string(),
                value: b"bounded".to_vec(),
            },
        ];
        write_exr_with_attributes(4, 3, &channels, &attributes).unwrap()
    }

    const RAW_TEST_ALLOWED_VALIDITY_BITS: u32 = RAW_VALIDITY_PRIMARY
        | RAW_VALIDITY_OBJECT_ID
        | RAW_VALIDITY_MATERIAL_ID
        | RAW_VALIDITY_CONTRIBUTION_SPLIT;

    fn raw_payload_fixture() -> Vec<u8> {
        let channels = vec![
            Channel {
                name: "B".to_string(),
                ty: PixelType::Float,
                data: vec![0.25; 6],
            },
            Channel {
                name: "diagnostic.validity".to_string(),
                ty: PixelType::Float,
                data: vec![
                    RAW_VALIDITY_CONTRIBUTION_SPLIT as f32,
                    RAW_TEST_ALLOWED_VALIDITY_BITS as f32,
                    RAW_TEST_ALLOWED_VALIDITY_BITS as f32,
                    RAW_TEST_ALLOWED_VALIDITY_BITS as f32,
                    RAW_TEST_ALLOWED_VALIDITY_BITS as f32,
                    RAW_TEST_ALLOWED_VALIDITY_BITS as f32,
                ],
            },
            Channel {
                name: "id.material".to_string(),
                ty: PixelType::Float,
                data: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            },
            Channel {
                name: "id.object".to_string(),
                ty: PixelType::Float,
                data: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            },
            Channel {
                name: "primary.coverage".to_string(),
                ty: PixelType::Float,
                data: vec![-0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
            },
            Channel {
                name: "samples".to_string(),
                ty: PixelType::Float,
                data: vec![1.0, 2.0, 4.0, 8.0, 16.0, MAX_CONSECUTIVE_F32_INTEGER],
            },
            Channel {
                name: "variance.Y".to_string(),
                ty: PixelType::Float,
                data: vec![0.0, 0.01, 0.1, 1.0, 10.0, -0.0],
            },
        ];
        write_exr(3, 2, &channels).unwrap()
    }

    fn overwrite_float_sample(
        bytes: &mut [u8],
        channel_name: &str,
        x: usize,
        y: usize,
        value: f32,
    ) {
        let inspection = inspect_exr(bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        let width = inspection.width as usize;
        let mut channel_row_offset = 0usize;
        for channel in &inspection.channels {
            if channel.name == channel_name {
                assert_eq!(channel.ty, PixelType::Float);
                break;
            }
            channel_row_offset += width * channel.ty.bytes();
        }
        let sample_offset = inspection.header_bytes as usize
            + inspection.height as usize * 8
            + y * (inspection.scanline_bytes as usize + 8)
            + 8
            + channel_row_offset
            + x * size_of::<f32>();
        drop(inspection);
        bytes[sample_offset..sample_offset + 4].copy_from_slice(&value.to_bits().to_le_bytes());
    }

    #[test]
    fn strict_exr_inspection_is_plane_free_and_exactly_budgeted() {
        let bytes = inspected_fixture();
        let inspected = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        assert_eq!((inspected.width, inspected.height), (4, 3));
        assert_eq!(
            inspected
                .channels
                .iter()
                .map(|channel| (channel.name, channel.ty))
                .collect::<Vec<_>>(),
            vec![("R", PixelType::Float), ("depth.Z", PixelType::Half)]
        );
        assert_eq!(
            inspected
                .attributes
                .iter()
                .map(|attribute| attribute.name)
                .collect::<Vec<_>>(),
            vec!["frankensim.alpha", "frankensim.zeta"]
        );
        assert_eq!(inspected.scanline_bytes, 4 * (4 + 2));
        assert_eq!(inspected.decoded_bytes, 3 * inspected.scanline_bytes);
        assert_eq!(inspected.input_bytes, bytes.len() as u64);

        let exact = ExrInspectLimits {
            max_input_bytes: inspected.input_bytes,
            max_header_bytes: inspected.header_bytes,
            max_decoded_bytes: inspected.decoded_bytes,
            max_metadata_bytes: inspected.metadata_bytes,
        };
        let input_bytes = inspected.input_bytes;
        let header_bytes = inspected.header_bytes;
        let decoded_bytes = inspected.decoded_bytes;
        let metadata_bytes = inspected.metadata_bytes;
        drop(inspected);
        assert!(
            inspect_exr(&bytes, exact).is_ok(),
            "equality must be admitted"
        );

        for (limits, resource) in [
            (
                ExrInspectLimits {
                    max_input_bytes: input_bytes - 1,
                    ..exact
                },
                "EXR input",
            ),
            (
                ExrInspectLimits {
                    max_header_bytes: header_bytes - 1,
                    ..exact
                },
                "EXR header",
            ),
            (
                ExrInspectLimits {
                    max_decoded_bytes: decoded_bytes - 1,
                    ..exact
                },
                "EXR decoded samples",
            ),
            (
                ExrInspectLimits {
                    max_metadata_bytes: metadata_bytes - 1,
                    ..exact
                },
                "EXR inspection metadata",
            ),
        ] {
            assert!(matches!(
                inspect_exr(&bytes, limits),
                Err(ImgError::ResourceLimit {
                    resource: observed,
                    ..
                }) if observed == resource
            ));
        }
    }

    #[test]
    fn strict_exr_inspection_rejects_offsets_sizes_and_trailing_bytes() {
        let bytes = inspected_fixture();
        let inspected = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        let table_start = inspected.header_bytes as usize;
        let blocks_start = table_start + inspected.height as usize * 8;
        drop(inspected);

        let mut bad_offset = bytes.clone();
        let offset =
            u64::from_le_bytes(bad_offset[table_start..table_start + 8].try_into().unwrap());
        bad_offset[table_start..table_start + 8].copy_from_slice(&(offset + 1).to_le_bytes());
        assert!(matches!(
            inspect_exr(&bad_offset, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("offset table")
        ));

        let mut bad_size = bytes.clone();
        let declared = i32::from_le_bytes(
            bad_size[blocks_start + 4..blocks_start + 8]
                .try_into()
                .unwrap(),
        );
        bad_size[blocks_start + 4..blocks_start + 8].copy_from_slice(&(declared - 1).to_le_bytes());
        assert!(matches!(
            inspect_exr(&bad_size, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("declares")
        ));

        let mut bad_y = bytes.clone();
        bad_y[blocks_start..blocks_start + 4].copy_from_slice(&1i32.to_le_bytes());
        assert!(matches!(
            inspect_exr(&bad_y, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("scanline order")
        ));

        let mut bad_sampling = bytes.clone();
        let channels_prefix = b"channels\0chlist\0";
        let channels_pos = bad_sampling
            .windows(channels_prefix.len())
            .position(|window| window == channels_prefix)
            .unwrap();
        let channel_value_start = channels_pos + channels_prefix.len() + 4;
        let channel_name_end = bad_sampling[channel_value_start..]
            .iter()
            .position(|&byte| byte == 0)
            .unwrap();
        let x_sampling = channel_value_start + channel_name_end + 1 + 8;
        bad_sampling[x_sampling..x_sampling + 4].copy_from_slice(&2u32.to_le_bytes());
        assert!(matches!(
            inspect_exr(&bad_sampling, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Unsupported { what }) if what.contains("subsampling")
        ));

        let mut bad_line_order = bytes.clone();
        let line_order_prefix = b"lineOrder\0lineOrder\0";
        let line_order_pos = bad_line_order
            .windows(line_order_prefix.len())
            .position(|window| window == line_order_prefix)
            .unwrap();
        bad_line_order[line_order_pos + line_order_prefix.len() + 4] = 1;
        assert!(matches!(
            inspect_exr(&bad_line_order, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Unsupported { what }) if what.contains("lineOrder")
        ));

        let mut trailing = bytes;
        trailing.push(0);
        assert!(matches!(
            inspect_exr(&trailing, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("encoded length")
        ));
    }

    #[test]
    fn strict_exr_inspection_cancels_at_bounded_units() {
        let bytes = inspected_fixture();
        let mut total_polls = 0usize;
        inspect_exr_with_poll(&bytes, ExrInspectLimits::UNBOUNDED, || {
            total_polls += 1;
            true
        })
        .unwrap();
        assert!(total_polls > 10, "fixture must cross multiple poll units");

        let cancel_at = total_polls / 2;
        let mut observed = 0usize;
        let result = inspect_exr_with_poll(&bytes, ExrInspectLimits::UNBOUNDED, || {
            observed += 1;
            observed < cancel_at
        });
        assert!(matches!(
            result,
            Err(ImgError::Cancelled {
                operation: "EXR structural inspection"
            })
        ));
        assert_eq!(observed, cancel_at);
    }

    #[test]
    fn strict_exr_inspection_rejects_every_truncated_prefix() {
        let bytes = inspected_fixture();
        for end in 0..bytes.len() {
            assert!(
                inspect_exr(&bytes[..end], ExrInspectLimits::UNBOUNDED).is_err(),
                "truncated prefix of {end} bytes was admitted"
            );
        }
    }

    #[test]
    fn exact_float_channel_verifier_accepts_bits_and_reports_first_mismatch() {
        let bytes = inspected_fixture();
        verify_exr_float_channel_constant(&bytes, "R", 0.25, ExrInspectLimits::UNBOUNDED).unwrap();

        let inspected = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        let blocks_start = inspected.header_bytes as usize + inspected.height as usize * 8;
        let scanline_block_bytes = inspected.scanline_bytes as usize + 8;
        drop(inspected);
        let mut mismatching = bytes;
        let sample_offset = blocks_start + scanline_block_bytes + 8 + 2 * 4;
        mismatching[sample_offset..sample_offset + 4]
            .copy_from_slice(&0.5f32.to_bits().to_le_bytes());

        assert!(matches!(
            verify_exr_float_channel_constant(
                &mismatching,
                "R",
                0.25,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { what })
                if what.contains("sample at (2, 1)")
                    && what.contains("0x3f000000")
                    && what.contains("0x3e800000")
        ));
    }

    #[test]
    fn exact_float_channel_verifier_rejects_missing_and_non_float_channels() {
        let bytes = inspected_fixture();
        assert!(matches!(
            verify_exr_float_channel_constant(
                &bytes,
                "missing",
                0.25,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { what }) if what.contains("is missing")
        ));
        assert!(matches!(
            verify_exr_float_channel_constant(
                &bytes,
                "depth.Z",
                1.0,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Unsupported { what })
                if what.contains("depth.Z") && what.contains("FLOAT is required")
        ));
    }

    #[test]
    fn exact_float_channel_verifier_preserves_structure_and_budget_refusals() {
        let bytes = inspected_fixture();
        let inspected = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        let decoded_bytes = inspected.decoded_bytes;
        let table_start = inspected.header_bytes as usize;
        drop(inspected);

        assert!(matches!(
            verify_exr_float_channel_constant(
                &bytes,
                "R",
                0.25,
                ExrInspectLimits {
                    max_decoded_bytes: decoded_bytes - 1,
                    ..ExrInspectLimits::UNBOUNDED
                },
            ),
            Err(ImgError::ResourceLimit {
                resource: "EXR decoded samples",
                requested,
                limit,
            }) if requested == decoded_bytes && limit == decoded_bytes - 1
        ));

        let mut invalid_offset = bytes.clone();
        invalid_offset[table_start..table_start + 8].copy_from_slice(&0_u64.to_le_bytes());
        assert!(matches!(
            verify_exr_float_channel_constant(
                &invalid_offset,
                "R",
                0.25,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { what }) if what.contains("offset table")
        ));
        assert!(matches!(
            verify_exr_float_channel_constant(
                &bytes[..bytes.len() - 1],
                "R",
                0.25,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { .. })
        ));
    }

    #[test]
    fn exact_float_channel_verifier_cancels_during_bounded_sample_scan() {
        let bytes = inspected_fixture();
        let mut structural_polls = 0usize;
        inspect_exr_with_poll(&bytes, ExrInspectLimits::UNBOUNDED, || {
            structural_polls += 1;
            true
        })
        .unwrap();

        let mut polls = 0usize;
        let result = verify_exr_float_channel_constant_with_poll(
            &bytes,
            "R",
            0.25,
            ExrInspectLimits::UNBOUNDED,
            || {
                polls += 1;
                polls <= structural_polls
            },
        );
        assert!(matches!(
            result,
            Err(ImgError::Cancelled {
                operation: "EXR FLOAT channel constant verification"
            })
        ));
        assert_eq!(polls, structural_polls + 1);
    }

    #[test]
    fn raw_frame_payload_validator_accepts_all_canonical_semantic_boundaries() {
        let mut bytes = raw_payload_fixture();
        overwrite_float_sample(&mut bytes, "id.object", 2, 1, MAX_CONSECUTIVE_F32_INTEGER);
        overwrite_float_sample(&mut bytes, "id.material", 2, 1, MAX_CONSECUTIVE_F32_INTEGER);
        validate_exr_raw_frame_payload(
            &bytes,
            RAW_TEST_ALLOWED_VALIDITY_BITS,
            ExrInspectLimits::UNBOUNDED,
        )
        .unwrap();
    }

    #[test]
    fn raw_frame_payload_validator_enforces_exact_palette_maxima() {
        let bytes = raw_payload_fixture();
        let exact = ExrRawFrameSemanticLimits::try_new(RAW_TEST_ALLOWED_VALIDITY_BITS, 5, 5)
            .expect("fixture palette limits are exactly representable");
        assert_eq!(
            exact.allowed_validity_bits(),
            RAW_TEST_ALLOWED_VALIDITY_BITS
        );
        assert_eq!(exact.maximum_object_palette_index(), 5);
        assert_eq!(exact.maximum_material_palette_index(), 5);
        assert!(matches!(
            ExrRawFrameSemanticLimits::try_new(
                RAW_TEST_ALLOWED_VALIDITY_BITS,
                MAX_CONSECUTIVE_F32_INTEGER as u32 + 1,
                5,
            ),
            Err(ImgError::Unsupported { what })
                if what.contains("object palette maximum")
                    && what.contains("consecutive-EXR-FLOAT")
        ));
        validate_exr_raw_frame_payload_against(&bytes, exact, ExrInspectLimits::UNBOUNDED).unwrap();

        for (channel, field) in [("id.object", "object"), ("id.material", "material")] {
            let mut hostile = bytes.clone();
            overwrite_float_sample(&mut hostile, channel, 1, 0, 6.0);
            assert!(matches!(
                validate_exr_raw_frame_payload_against(
                    &hostile,
                    exact,
                    ExrInspectLimits::UNBOUNDED,
                ),
                Err(ImgError::Malformed { what })
                    if what.contains(channel)
                        && what.contains("sample at (1, 0)")
                        && what.contains("palette index 6")
                        && what.contains("maximum 5")
                        && what.contains(field)
            ));
        }
    }

    #[test]
    fn raw_frame_payload_validator_enforces_per_pixel_validity_relations() {
        let fixture = raw_payload_fixture();
        let exact = ExrRawFrameSemanticLimits::try_new(RAW_TEST_ALLOWED_VALIDITY_BITS, 5, 5)
            .expect("fixture palette limits are exactly representable");
        for (channel, value, relation) in [
            ("id.object", 0.0, "OBJECT_ID iff id.object != 0"),
            ("id.material", 0.0, "MATERIAL_ID iff id.material != 0"),
            ("primary.coverage", 0.0, "PRIMARY iff primary.coverage > 0"),
            (
                "diagnostic.validity",
                (RAW_TEST_ALLOWED_VALIDITY_BITS & !RAW_VALIDITY_CONTRIBUTION_SPLIT) as f32,
                "CONTRIBUTION_SPLIT iff samples > 0",
            ),
        ] {
            let mut hostile = fixture.clone();
            overwrite_float_sample(&mut hostile, channel, 1, 0, value);
            assert!(matches!(
                validate_exr_raw_frame_payload_against(
                    &hostile,
                    exact,
                    ExrInspectLimits::UNBOUNDED,
                ),
                Err(ImgError::Malformed { what })
                    if what.contains("at (1, 0)") && what.contains(relation)
            ));
        }

        let mut orphaned_id = fixture;
        overwrite_float_sample(&mut orphaned_id, "primary.coverage", 1, 0, 0.0);
        overwrite_float_sample(
            &mut orphaned_id,
            "diagnostic.validity",
            1,
            0,
            (RAW_TEST_ALLOWED_VALIDITY_BITS & !RAW_VALIDITY_PRIMARY) as f32,
        );
        assert!(matches!(
            validate_exr_raw_frame_payload_against(
                &orphaned_id,
                exact,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { what })
                if what.contains("nonzero object/material ID requires")
                    && what.contains("at (1, 0)")
        ));
    }

    #[test]
    fn raw_frame_payload_validator_rejects_reserved_validity_bit() {
        let mut hostile = raw_payload_fixture();
        overwrite_float_sample(
            &mut hostile,
            "diagnostic.validity",
            1,
            0,
            (RAW_TEST_ALLOWED_VALIDITY_BITS | (1 << 2)) as f32,
        );
        let exact = ExrRawFrameSemanticLimits::try_new(RAW_TEST_ALLOWED_VALIDITY_BITS, 5, 5)
            .expect("fixture palette limits are exactly representable");
        assert!(matches!(
            validate_exr_raw_frame_payload_against(
                &hostile,
                exact,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { what })
                if what.contains("diagnostic.validity")
                    && what.contains("sample at (1, 0)")
                    && what.contains("unknown bits 0x00000004")
        ));
    }

    #[test]
    fn raw_frame_payload_validator_reports_first_hostile_sample_exactly() {
        let fixture = raw_payload_fixture();
        for (channel, value, reason) in [
            ("B", f32::INFINITY, "must be finite"),
            ("id.object", -1.0, "nonnegative exact integer"),
            ("id.material", 0.5, "nonnegative exact integer"),
            ("id.object", 16_777_218.0, "consecutive-f32 range"),
            ("samples", 0.0, "positive exact integer"),
            ("samples", 1.5, "positive exact integer"),
            ("diagnostic.validity", 4.0, "unknown bits 0x00000004"),
            ("diagnostic.validity", 0.5, "exact integer bit mask"),
            ("primary.coverage", 1.25, "closed interval [0, 1]"),
            ("variance.Y", -0.25, "must be nonnegative"),
        ] {
            let mut hostile = fixture.clone();
            overwrite_float_sample(&mut hostile, channel, 1, 1, value);
            assert!(matches!(
                validate_exr_raw_frame_payload(
                    &hostile,
                    RAW_TEST_ALLOWED_VALIDITY_BITS,
                    ExrInspectLimits::UNBOUNDED,
                ),
                Err(ImgError::Malformed { what })
                    if what.contains(channel)
                        && what.contains("sample at (1, 1)")
                        && what.contains(reason)
            ));
        }
    }

    #[test]
    fn raw_frame_payload_validator_rejects_non_float_semantic_channel() {
        let bytes = write_exr(
            1,
            1,
            &[
                Channel {
                    name: "R".to_string(),
                    ty: PixelType::Float,
                    data: vec![0.25],
                },
                Channel {
                    name: "samples".to_string(),
                    ty: PixelType::Half,
                    data: vec![1.0],
                },
            ],
        )
        .unwrap();
        assert!(matches!(
            validate_exr_raw_frame_payload(&bytes, 0, ExrInspectLimits::UNBOUNDED),
            Err(ImgError::Unsupported { what })
                if what.contains("samples") && what.contains("FLOAT is required")
        ));
    }

    #[test]
    fn raw_frame_payload_validator_preserves_budget_and_structure_refusals() {
        let bytes = raw_payload_fixture();
        let inspection = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        let decoded_bytes = inspection.decoded_bytes;
        drop(inspection);
        assert!(matches!(
            validate_exr_raw_frame_payload(
                &bytes,
                RAW_TEST_ALLOWED_VALIDITY_BITS,
                ExrInspectLimits {
                    max_decoded_bytes: decoded_bytes - 1,
                    ..ExrInspectLimits::UNBOUNDED
                },
            ),
            Err(ImgError::ResourceLimit {
                resource: "EXR decoded samples",
                requested,
                limit,
            }) if requested == decoded_bytes && limit == decoded_bytes - 1
        ));
        assert!(matches!(
            validate_exr_raw_frame_payload(
                &bytes[..bytes.len() - 1],
                RAW_TEST_ALLOWED_VALIDITY_BITS,
                ExrInspectLimits::UNBOUNDED,
            ),
            Err(ImgError::Malformed { .. })
        ));
    }

    #[test]
    fn raw_frame_payload_validator_cancels_in_the_payload_scan() {
        let bytes = raw_payload_fixture();
        let mut structural_polls = 0usize;
        inspect_exr_with_poll(&bytes, ExrInspectLimits::UNBOUNDED, || {
            structural_polls += 1;
            true
        })
        .unwrap();
        let semantic_channels = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED)
            .unwrap()
            .channels
            .len();

        let mut polls = 0usize;
        let result = validate_exr_raw_frame_payload_with_poll(
            &bytes,
            RAW_TEST_ALLOWED_VALIDITY_BITS,
            ExrInspectLimits::UNBOUNDED,
            || {
                polls += 1;
                polls <= structural_polls + semantic_channels
            },
        );
        assert!(matches!(
            result,
            Err(ImgError::Cancelled {
                operation: "EXR raw-frame payload validation"
            })
        ));
        assert_eq!(polls, structural_polls + semantic_channels + 1);

        let inspection = inspect_exr(&bytes, ExrInspectLimits::UNBOUNDED).unwrap();
        let payload_polls = inspection.channels.len() * inspection.height as usize;
        let correlation_start = structural_polls + inspection.channels.len() + payload_polls;
        let exact = ExrRawFrameSemanticLimits::try_new(RAW_TEST_ALLOWED_VALIDITY_BITS, 5, 5)
            .expect("fixture palette limits are exactly representable");
        let mut polls = 0usize;
        let result = validate_exr_raw_frame_payload_against_with_poll(
            &bytes,
            exact,
            ExrInspectLimits::UNBOUNDED,
            || {
                polls += 1;
                polls <= correlation_start
            },
        );
        assert!(matches!(
            result,
            Err(ImgError::Cancelled {
                operation: "EXR raw-frame semantic relation validation"
            })
        ));
        assert_eq!(polls, correlation_start + 1);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one table proves the complete frozen 4K AOV layout
    fn exr_requirements_are_exact_and_budgets_refuse_before_encoding() {
        let channels = vec![
            Channel {
                name: "R".to_string(),
                ty: PixelType::Float,
                data: vec![0.25, 0.5],
            },
            Channel {
                name: "depth.Z".to_string(),
                ty: PixelType::Half,
                data: vec![1.0, 2.0],
            },
        ];
        let attributes = vec![ExrAttribute {
            name: "frankensim.test".to_string(),
            ty: "string".to_string(),
            value: b"bounded".to_vec(),
        }];
        let requirements = exr_write_requirements(2, 1, &channels, &attributes).unwrap();
        let expected_output = 259_u64
            + channels
                .iter()
                .map(|channel| channel.name.len() as u64 + 17)
                .sum::<u64>()
            + attributes
                .iter()
                .map(|attribute| {
                    attribute.name.len() as u64
                        + attribute.ty.len() as u64
                        + 6
                        + attribute.value.len() as u64
                })
                .sum::<u64>()
            + 16
            + 2 * (4 + 2);
        assert_eq!(requirements.output_bytes, expected_output);
        assert_eq!(
            requirements.scratch_bytes,
            2 * size_of::<&Channel>() as u64 + size_of::<&ExrAttribute>() as u64
        );

        let legacy = write_exr_with_attributes(2, 1, &channels, &attributes).unwrap();
        assert_eq!(legacy.len() as u64, requirements.output_bytes);
        assert_eq!(
            write_exr_with_attributes_budgeted(
                2,
                1,
                &channels,
                &attributes,
                ExrWriteLimits {
                    max_scratch_bytes: requirements.scratch_bytes,
                    max_output_bytes: requirements.output_bytes,
                },
            )
            .unwrap(),
            legacy,
            "adequate ceilings must preserve the legacy deterministic bytes"
        );
        assert!(matches!(
            write_exr_with_attributes_budgeted(
                2,
                1,
                &channels,
                &attributes,
                ExrWriteLimits {
                    max_scratch_bytes: requirements.scratch_bytes - 1,
                    max_output_bytes: requirements.output_bytes,
                },
            ),
            Err(ImgError::ResourceLimit {
                resource: "EXR ordering scratch",
                ..
            })
        ));
        assert!(matches!(
            write_exr_with_attributes_budgeted(
                2,
                1,
                &channels,
                &attributes,
                ExrWriteLimits {
                    max_scratch_bytes: requirements.scratch_bytes,
                    max_output_bytes: requirements.output_bytes - 1,
                },
            ),
            Err(ImgError::ResourceLimit {
                resource: "EXR encoded output",
                ..
            })
        ));

        let final_layout = [
            ("R", PixelType::Float),
            ("G", PixelType::Float),
            ("B", PixelType::Float),
            ("albedo.R", PixelType::Float),
            ("albedo.G", PixelType::Float),
            ("albedo.B", PixelType::Float),
            ("normal.X", PixelType::Float),
            ("normal.Y", PixelType::Float),
            ("normal.Z", PixelType::Float),
            ("depth.Z", PixelType::Float),
            ("primary.coverage", PixelType::Float),
            ("variance.Y", PixelType::Float),
            ("motion.prev.X", PixelType::Float),
            ("motion.prev.Y", PixelType::Float),
            ("normal_geom.X", PixelType::Float),
            ("normal_geom.Y", PixelType::Float),
            ("normal_geom.Z", PixelType::Float),
            ("direct.R", PixelType::Float),
            ("direct.G", PixelType::Float),
            ("direct.B", PixelType::Float),
            ("indirect.R", PixelType::Float),
            ("indirect.G", PixelType::Float),
            ("indirect.B", PixelType::Float),
            ("emission.R", PixelType::Float),
            ("emission.G", PixelType::Float),
            ("emission.B", PixelType::Float),
            ("id.object", PixelType::Float),
            ("id.material", PixelType::Float),
            ("samples", PixelType::Float),
            ("diagnostic.validity", PixelType::Float),
        ];
        let four_k = exr_write_requirements_for_layout(3_840, 2_160, &final_layout, &[]).unwrap();
        assert_eq!(four_k.output_bytes, 995_363_608);
        assert_eq!(four_k.scratch_bytes, 30 * size_of::<&Channel>() as u64);
    }

    #[test]
    fn exr_v2_short_names_and_signed_lengths_enforce_format_boundaries() {
        let name_31 = "x".repeat(31);
        let name_32 = "x".repeat(32);
        assert!(
            exr_write_requirements_for_layout(1, 1, &[(name_31.as_str(), PixelType::Float)], &[],)
                .is_ok(),
            "version 2 permits a 31-byte short channel name"
        );
        assert!(matches!(
            exr_write_requirements_for_layout(1, 1, &[(name_32.as_str(), PixelType::Float)], &[],),
            Err(ImgError::Malformed { .. })
        ));

        // Exercise the scanline limit through the allocation-free layout API:
        // 536,870,911 FLOAT samples occupy i32::MAX - 3 bytes, while one more
        // sample would set the sign bit in the EXR scanline-size field.
        let channel = [("R", PixelType::Float)];
        assert!(exr_write_requirements_for_layout(536_870_911, 1, &channel, &[]).is_ok());
        assert!(matches!(
            exr_write_requirements_for_layout(536_870_912, 1, &channel, &[]),
            Err(ImgError::SizeOverflow {
                context: "EXR signed attribute or scanline length"
            })
        ));

        assert!(encoded_attribute_len("a", "b", i32::MAX as u64).is_ok());
        assert!(matches!(
            encoded_attribute_len("a", "b", i32::MAX as u64 + 1),
            Err(ImgError::SizeOverflow {
                context: "EXR attribute signed length"
            })
        ));
    }

    #[test]
    fn malformed_and_unsupported_reject() {
        assert!(read_exr(b"nope").is_err());
        let ch = Channel {
            name: "R".to_string(),
            ty: PixelType::Float,
            data: vec![0.0; 4],
        };
        let mut bytes = write_exr(2, 2, std::slice::from_ref(&ch)).unwrap();
        // Flip compression byte to ZIP: structured Unsupported.
        let pos = bytes
            .windows(12)
            .position(|w| w.starts_with(b"compression\0"))
            .expect("attr present");
        // name + NUL + type("compression") + NUL + size(4) → value byte.
        let value_at = pos + 12 + 12 + 4;
        bytes[value_at] = 3;
        assert!(matches!(
            read_exr(&bytes),
            Err(ImgError::Unsupported { .. })
        ));
        // Duplicate channel names refuse at write time.
        assert!(write_exr(2, 2, &[ch.clone(), ch]).is_err());

        let reserved = ExrAttribute {
            name: "channels".to_string(),
            ty: "string".to_string(),
            value: Vec::new(),
        };
        assert!(
            write_exr_with_attributes(
                2,
                2,
                std::slice::from_ref(&Channel {
                    name: "R".to_string(),
                    ty: PixelType::Float,
                    data: vec![0.0; 4],
                }),
                std::slice::from_ref(&reserved)
            )
            .is_err()
        );
    }
}
