//! In-house OpenEXR writer + reader (plan §10.5), spec-conformant subset:
//! single-part scanline files, NONE compression, HALF/FLOAT channels,
//! multi-channel AOVs (beauty/albedo/normal/depth) with the spec's
//! alphabetical channel ordering, and opaque custom header attributes. The
//! reader covers our writer's subset (round-trips + ledger artifacts) with
//! structured rejections beyond it.
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

    #[test]
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
