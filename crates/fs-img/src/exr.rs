//! In-house OpenEXR writer + reader (plan §10.5), spec-conformant subset:
//! single-part scanline files, NONE compression, HALF/FLOAT channels,
//! multi-channel AOVs (beauty/albedo/normal/depth) with the spec's
//! alphabetical channel ordering. The reader covers our writer's subset
//! (round-trips + ledger artifacts) with structured rejections beyond it.
//!
//! Determinism: byte-exact encodes (pure integer/bit code).

use crate::ImgError;
use std::collections::BTreeMap;

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
    /// Channel name (e.g. "R", "albedo.G", "depth.Z"); no NULs.
    pub name: String,
    /// Storage type.
    pub ty: PixelType,
    /// Row-major samples (width × height).
    pub data: Vec<f32>,
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
    out.extend_from_slice(&(value.len() as u32).to_le_bytes());
    out.extend_from_slice(value);
}

/// Encode channels as a single-part scanline EXR (NONE compression).
/// Channels are stored in the spec's alphabetical order regardless of the
/// argument order.
///
/// # Errors
/// [`ImgError`] on shape/name defects.
pub fn write_exr(width: u32, height: u32, channels: &[Channel]) -> Result<Vec<u8>, ImgError> {
    if width == 0 || height == 0 || channels.is_empty() {
        return Err(ImgError::Shape {
            expected: 1,
            got: 0,
            context: "write_exr needs a nonempty image and channel set",
        });
    }
    let (Ok(wi), Ok(hi)) = (i32::try_from(width), i32::try_from(height)) else {
        return Err(ImgError::Malformed {
            what: format!("dimensions {width}x{height} exceed the EXR i32 data window"),
        });
    };
    let n = width as usize * height as usize;
    let mut sorted: BTreeMap<&str, &Channel> = BTreeMap::new();
    for c in channels {
        if c.name.is_empty() || c.name.contains('\0') {
            return Err(ImgError::Malformed {
                what: format!("channel name {:?} (empty/NUL)", c.name),
            });
        }
        if c.data.len() != n {
            return Err(ImgError::Shape {
                expected: n,
                got: c.data.len(),
                context: "channel sample count",
            });
        }
        if sorted.insert(&c.name, c).is_some() {
            return Err(ImgError::Malformed {
                what: format!("duplicate channel {:?}", c.name),
            });
        }
    }

    let mut out = Vec::new();
    out.extend_from_slice(&MAGIC);
    out.extend_from_slice(&2u32.to_le_bytes()); // version 2, no flags

    // chlist attribute value.
    let mut chlist = Vec::new();
    for (name, c) in &sorted {
        chlist.extend_from_slice(name.as_bytes());
        chlist.push(0);
        chlist.extend_from_slice(&c.ty.code().to_le_bytes());
        chlist.extend_from_slice(&[0, 0, 0, 0]); // pLinear + reserved
        chlist.extend_from_slice(&1u32.to_le_bytes()); // xSampling
        chlist.extend_from_slice(&1u32.to_le_bytes()); // ySampling
    }
    chlist.push(0);
    push_attr(&mut out, "channels", "chlist", &chlist);
    push_attr(&mut out, "compression", "compression", &[0]); // NONE
    let mut window = Vec::with_capacity(16);
    window.extend_from_slice(&0i32.to_le_bytes());
    window.extend_from_slice(&0i32.to_le_bytes());
    window.extend_from_slice(&(wi - 1).to_le_bytes());
    window.extend_from_slice(&(hi - 1).to_le_bytes());
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
    out.push(0); // end of header

    // Scanline offset table placeholder.
    let table_pos = out.len();
    out.resize(out.len() + 8 * height as usize, 0);

    let line_bytes: usize = sorted.values().map(|c| width as usize * c.ty.bytes()).sum();
    for y in 0..height as usize {
        let offset = out.len() as u64;
        out[table_pos + 8 * y..table_pos + 8 * (y + 1)].copy_from_slice(&offset.to_le_bytes());
        out.extend_from_slice(&i32::try_from(y).expect("y < height <= i32::MAX").to_le_bytes());
        out.extend_from_slice(&(line_bytes as u32).to_le_bytes());
        for c in sorted.values() {
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
    Ok(out)
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
}

fn take(bytes: &[u8], pos: usize, n: usize) -> Result<&[u8], ImgError> {
    bytes.get(pos..pos + n).ok_or_else(|| ImgError::Malformed {
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
    let s = String::from_utf8_lossy(&bytes[start..*pos]).into_owned();
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

/// Parsed header: (channel specs, width, height).
type HeaderInfo = (Vec<(String, PixelType)>, u32, u32);

/// Parse the header attributes; returns [`HeaderInfo`] and leaves `pos`
/// just past the header terminator.
fn parse_header(bytes: &[u8], pos: &mut usize) -> Result<HeaderInfo, ImgError> {
    let mut specs: Vec<(String, PixelType)> = Vec::new();
    let mut window = (0u32, 0u32);
    let mut compression_seen = false;
    loop {
        if bytes.get(*pos) == Some(&0) {
            *pos += 1;
            break; // end of header
        }
        let name = read_cstr(bytes, pos)?;
        let _ty = read_cstr(bytes, pos)?;
        let size = u32::from_le_bytes(take(bytes, *pos, 4)?.try_into().expect("4 bytes")) as usize;
        *pos += 4;
        let value = take(bytes, *pos, size)?.to_vec();
        *pos += size;
        match name.as_str() {
            "channels" => specs = parse_chlist(&value)?,
            "compression" => {
                compression_seen = true;
                if value != [0] {
                    return Err(ImgError::Unsupported {
                        what: format!("compression {} (NONE only)", value[0]),
                    });
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
            _ => {}
        }
    }
    if !compression_seen || specs.is_empty() || window.0 == 0 || window.1 == 0 {
        return Err(ImgError::Malformed {
            what: "missing required header attributes".to_string(),
        });
    }
    Ok((specs, window.0, window.1))
}

/// Decode an EXR produced by [`write_exr`]. Structured rejection outside
/// that subset.
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
    let (specs, width, height) = parse_header(bytes, &mut pos)?;
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
    }
}
