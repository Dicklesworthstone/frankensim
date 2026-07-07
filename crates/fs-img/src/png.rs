//! In-house PNG writer + reader (plan §10.5), spec-conformant subset:
//! 8/16-bit grayscale/RGB/RGBA, sRGB chunk, None filters, zlib streams
//! built from STORED deflate blocks (universally decodable; compression
//! ratio is an explicit no-claim — renders ship in EXR, PNG is the
//! preview/report format).
//!
//! Determinism: byte-exact encodes — same pixels, same bytes, every run,
//! every ISA (pure integer code; golden-hashed in conformance).
//!
//! The reader covers exactly OUR writer's subset (round-trips + ledger
//! artifact loading) and rejects everything else with structured errors —
//! it is not a general PNG decoder (documented no-claim).

use crate::ImgError;

/// PNG color layouts this writer speaks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PngColor {
    /// 1 channel.
    Gray,
    /// 3 channels.
    Rgb,
    /// 4 channels.
    Rgba,
}

impl PngColor {
    fn channels(self) -> usize {
        match self {
            PngColor::Gray => 1,
            PngColor::Rgb => 3,
            PngColor::Rgba => 4,
        }
    }

    fn type_byte(self) -> u8 {
        match self {
            PngColor::Gray => 0,
            PngColor::Rgb => 2,
            PngColor::Rgba => 6,
        }
    }

    fn from_type_byte(b: u8) -> Option<PngColor> {
        match b {
            0 => Some(PngColor::Gray),
            2 => Some(PngColor::Rgb),
            6 => Some(PngColor::Rgba),
            _ => None,
        }
    }
}

/// CRC-32 (IEEE 802.3, reflected 0xEDB88320) — the PNG chunk checksum.
#[must_use]
pub fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &b in bytes {
        crc ^= u32::from(b);
        for _ in 0..8 {
            let mask = 0u32.wrapping_sub(crc & 1);
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    !crc
}

/// Adler-32 — the zlib stream checksum.
#[must_use]
pub fn adler32(bytes: &[u8]) -> u32 {
    const MOD: u32 = 65_521;
    let (mut a, mut b) = (1u32, 0u32);
    for chunk in bytes.chunks(5000) {
        for &x in chunk {
            a += u32::from(x);
            b += a;
        }
        a %= MOD;
        b %= MOD;
    }
    (b << 16) | a
}

/// Wrap raw bytes in a zlib stream of STORED deflate blocks.
fn zlib_stored(raw: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(raw.len() + raw.len() / 65_535 * 5 + 16);
    out.extend_from_slice(&[0x78, 0x01]); // CMF/FLG (32K window, no dict)
    let mut chunks = raw.chunks(65_535).peekable();
    if raw.is_empty() {
        out.extend_from_slice(&[0x01, 0x00, 0x00, 0xFF, 0xFF]); // final empty block
    }
    while let Some(chunk) = chunks.next() {
        let bfinal = u8::from(chunks.peek().is_none());
        let len = chunk.len() as u16;
        out.push(bfinal); // BTYPE=00 stored, byte-aligned
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(&(!len).to_le_bytes());
        out.extend_from_slice(chunk);
    }
    out.extend_from_slice(&adler32(raw).to_be_bytes());
    out
}

/// Un-wrap a zlib stream of STORED blocks (our writer's subset).
fn unzlib_stored(z: &[u8]) -> Result<Vec<u8>, ImgError> {
    if z.len() < 6 {
        return Err(ImgError::Malformed {
            what: "zlib stream too short".to_string(),
        });
    }
    if z[0] & 0x0F != 8 {
        return Err(ImgError::Malformed {
            what: "not a deflate zlib stream".to_string(),
        });
    }
    let mut pos = 2usize;
    let mut out = Vec::new();
    loop {
        let Some(&header) = z.get(pos) else {
            return Err(ImgError::Malformed {
                what: "truncated deflate block".to_string(),
            });
        };
        if header & 0x06 != 0 {
            return Err(ImgError::Unsupported {
                what: "compressed deflate blocks (this reader covers our stored-block \
                       writer subset)"
                    .to_string(),
            });
        }
        let bfinal = header & 1 == 1;
        let len_bytes = z.get(pos + 1..pos + 5).ok_or_else(|| ImgError::Malformed {
            what: "truncated block header".to_string(),
        })?;
        let len = u16::from_le_bytes([len_bytes[0], len_bytes[1]]) as usize;
        let nlen = u16::from_le_bytes([len_bytes[2], len_bytes[3]]);
        if nlen != !(len as u16) {
            return Err(ImgError::Malformed {
                what: "stored block NLEN mismatch".to_string(),
            });
        }
        let data = z
            .get(pos + 5..pos + 5 + len)
            .ok_or_else(|| ImgError::Malformed {
                what: "stored block data truncated".to_string(),
            })?;
        out.extend_from_slice(data);
        pos += 5 + len;
        if bfinal {
            break;
        }
    }
    let adler = z.get(pos..pos + 4).ok_or_else(|| ImgError::Malformed {
        what: "missing adler32 trailer".to_string(),
    })?;
    if u32::from_be_bytes([adler[0], adler[1], adler[2], adler[3]]) != adler32(&out) {
        return Err(ImgError::Malformed {
            what: "adler32 mismatch (corrupt data)".to_string(),
        });
    }
    Ok(out)
}

fn push_chunk(out: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    out.extend_from_slice(&(data.len() as u32).to_be_bytes());
    out.extend_from_slice(kind);
    out.extend_from_slice(data);
    let mut crc_input = Vec::with_capacity(4 + data.len());
    crc_input.extend_from_slice(kind);
    crc_input.extend_from_slice(data);
    out.extend_from_slice(&crc32(&crc_input).to_be_bytes());
}

const SIGNATURE: [u8; 8] = [0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A];

/// Encode 8-bit samples (row-major, interleaved channels) as a PNG with
/// an sRGB chunk. Byte-exact deterministic.
///
/// # Errors
/// [`ImgError::Shape`] when the buffer disagrees with width × height ×
/// channels.
pub fn write_png8(
    width: u32,
    height: u32,
    color: PngColor,
    samples: &[u8],
) -> Result<Vec<u8>, ImgError> {
    let expected = width as usize * height as usize * color.channels();
    if samples.len() != expected || width == 0 || height == 0 {
        return Err(ImgError::Shape {
            expected,
            got: samples.len(),
            context: "write_png8 samples",
        });
    }
    let row = width as usize * color.channels();
    let mut raw = Vec::with_capacity((row + 1) * height as usize);
    for y in 0..height as usize {
        raw.push(0); // filter: None
        raw.extend_from_slice(&samples[y * row..(y + 1) * row]);
    }
    Ok(assemble(width, height, 8, color, &raw))
}

/// Encode 16-bit samples (row-major, interleaved; big-endian per spec).
///
/// # Errors
/// [`ImgError::Shape`] on buffer/shape disagreement.
pub fn write_png16(
    width: u32,
    height: u32,
    color: PngColor,
    samples: &[u16],
) -> Result<Vec<u8>, ImgError> {
    let expected = width as usize * height as usize * color.channels();
    if samples.len() != expected || width == 0 || height == 0 {
        return Err(ImgError::Shape {
            expected,
            got: samples.len(),
            context: "write_png16 samples",
        });
    }
    let row = width as usize * color.channels();
    let mut raw = Vec::with_capacity((row * 2 + 1) * height as usize);
    for y in 0..height as usize {
        raw.push(0);
        for &s in &samples[y * row..(y + 1) * row] {
            raw.extend_from_slice(&s.to_be_bytes());
        }
    }
    Ok(assemble(width, height, 16, color, &raw))
}

fn assemble(width: u32, height: u32, depth: u8, color: PngColor, raw: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(raw.len() + 128);
    out.extend_from_slice(&SIGNATURE);
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.push(depth);
    ihdr.push(color.type_byte());
    ihdr.extend_from_slice(&[0, 0, 0]); // deflate, adaptive filters, no interlace
    push_chunk(&mut out, b"IHDR", &ihdr);
    push_chunk(&mut out, b"sRGB", &[0]); // perceptual intent
    push_chunk(&mut out, b"IDAT", &zlib_stored(raw));
    push_chunk(&mut out, b"IEND", &[]);
    out
}

/// A decoded PNG (our writer's subset).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedPng {
    /// Pixel width.
    pub width: u32,
    /// Pixel height.
    pub height: u32,
    /// Bit depth (8 or 16).
    pub depth: u8,
    /// Color layout.
    pub color: PngColor,
    /// Interleaved samples as bytes (16-bit stays big-endian pairs; use
    /// [`DecodedPng::samples16`] for typed access).
    pub bytes: Vec<u8>,
}

impl DecodedPng {
    /// 16-bit samples (only valid when `depth == 16`).
    #[must_use]
    pub fn samples16(&self) -> Vec<u16> {
        self.bytes
            .chunks_exact(2)
            .map(|p| u16::from_be_bytes([p[0], p[1]]))
            .collect()
    }
}

/// Decode a PNG produced by [`write_png8`]/[`write_png16`]. Structured
/// rejection on anything outside that subset (fuzz-tested totality).
///
/// # Errors
/// [`ImgError::Malformed`] / [`ImgError::Unsupported`].
pub fn read_png(bytes: &[u8]) -> Result<DecodedPng, ImgError> {
    if bytes.len() < 8 || bytes[..8] != SIGNATURE {
        return Err(ImgError::Malformed {
            what: "missing PNG signature".to_string(),
        });
    }
    let mut pos = 8usize;
    let mut header: Option<(u32, u32, u8, PngColor)> = None;
    let mut idat = Vec::new();
    loop {
        let len_bytes = bytes.get(pos..pos + 8).ok_or_else(|| ImgError::Malformed {
            what: "truncated chunk header".to_string(),
        })?;
        let len =
            u32::from_be_bytes([len_bytes[0], len_bytes[1], len_bytes[2], len_bytes[3]]) as usize;
        let kind = &len_bytes[4..8];
        let data = bytes
            .get(pos + 8..pos + 8 + len)
            .ok_or_else(|| ImgError::Malformed {
                what: "truncated chunk data".to_string(),
            })?;
        let crc = bytes
            .get(pos + 8 + len..pos + 12 + len)
            .ok_or_else(|| ImgError::Malformed {
                what: "truncated chunk crc".to_string(),
            })?;
        let mut crc_input = Vec::with_capacity(4 + len);
        crc_input.extend_from_slice(kind);
        crc_input.extend_from_slice(data);
        if u32::from_be_bytes([crc[0], crc[1], crc[2], crc[3]]) != crc32(&crc_input) {
            return Err(ImgError::Malformed {
                what: "chunk crc mismatch".to_string(),
            });
        }
        match kind {
            b"IHDR" => {
                if data.len() != 13 {
                    return Err(ImgError::Malformed {
                        what: "IHDR length".to_string(),
                    });
                }
                let w = u32::from_be_bytes([data[0], data[1], data[2], data[3]]);
                let h = u32::from_be_bytes([data[4], data[5], data[6], data[7]]);
                let depth = data[8];
                let color =
                    PngColor::from_type_byte(data[9]).ok_or_else(|| ImgError::Unsupported {
                        what: format!("color type {}", data[9]),
                    })?;
                if depth != 8 && depth != 16 {
                    return Err(ImgError::Unsupported {
                        what: format!("bit depth {depth}"),
                    });
                }
                if data[12] != 0 {
                    return Err(ImgError::Unsupported {
                        what: "interlacing".to_string(),
                    });
                }
                header = Some((w, h, depth, color));
            }
            b"IDAT" => idat.extend_from_slice(data),
            b"IEND" => break,
            _ => {} // ancillary chunks (sRGB…) pass through
        }
        pos += 12 + len;
    }
    let Some((width, height, depth, color)) = header else {
        return Err(ImgError::Malformed {
            what: "no IHDR before IEND".to_string(),
        });
    };
    let raw = unzlib_stored(&idat)?;
    let bpp = color.channels() * (depth as usize / 8);
    let row = width as usize * bpp;
    let expected = (row + 1) * height as usize;
    if raw.len() != expected {
        return Err(ImgError::Shape {
            expected,
            got: raw.len(),
            context: "decoded scanlines",
        });
    }
    let mut out = Vec::with_capacity(row * height as usize);
    for y in 0..height as usize {
        let line = &raw[y * (row + 1)..(y + 1) * (row + 1)];
        if line[0] != 0 {
            return Err(ImgError::Unsupported {
                what: format!("filter type {} (our writer emits None)", line[0]),
            });
        }
        out.extend_from_slice(&line[1..]);
    }
    Ok(DecodedPng {
        width,
        height,
        depth,
        color,
        bytes: out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc_and_adler_known_answers() {
        // Standard test vectors.
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(adler32(b"Wikipedia"), 0x11E6_0398);
        assert_eq!(adler32(b""), 1);
    }

    #[test]
    fn png8_round_trips_bit_exactly() {
        let (w, h) = (5u32, 3u32);
        let px: Vec<u8> = (0..w * h * 3).map(|i| (i * 7 % 251) as u8).collect();
        let bytes = write_png8(w, h, PngColor::Rgb, &px).unwrap();
        let again = write_png8(w, h, PngColor::Rgb, &px).unwrap();
        assert_eq!(bytes, again, "byte-exact determinism");
        let decoded = read_png(&bytes).unwrap();
        assert_eq!((decoded.width, decoded.height, decoded.depth), (w, h, 8));
        assert_eq!(decoded.bytes, px, "pixel round-trip");
    }

    #[test]
    fn png16_round_trips() {
        let (w, h) = (4u32, 2u32);
        let px: Vec<u16> = (0..w * h * 4).map(|i| (i * 6151 % 65_521) as u16).collect();
        let bytes = write_png16(w, h, PngColor::Rgba, &px).unwrap();
        let decoded = read_png(&bytes).unwrap();
        assert_eq!(decoded.samples16(), px);
    }

    #[test]
    fn shape_and_malformed_rejections_teach() {
        assert!(matches!(
            write_png8(4, 4, PngColor::Rgb, &[0u8; 5]),
            Err(ImgError::Shape {
                expected: 48,
                got: 5,
                ..
            })
        ));
        assert!(read_png(b"not a png").is_err());
        // Corrupt one IDAT byte: crc must catch it.
        let px = vec![7u8; 12];
        let mut bytes = write_png8(2, 2, PngColor::Rgb, &px).unwrap();
        let idx = bytes.len() - 30;
        bytes[idx] ^= 0xFF;
        assert!(
            read_png(&bytes).is_err(),
            "corruption must not decode silently"
        );
    }
}
