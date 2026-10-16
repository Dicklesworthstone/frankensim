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
//! The bounded structural inspector verifies the canonical chunk and stored-
//! deflate grammar without retaining decoded pixels.

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
    !crc32_update(0xFFFF_FFFF, bytes)
}

fn crc32_update(mut crc: u32, bytes: &[u8]) -> u32 {
    for &b in bytes {
        crc ^= u32::from(b);
        for _ in 0..8 {
            let mask = 0u32.wrapping_sub(crc & 1);
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    crc
}

/// Adler-32 — the zlib stream checksum.
#[must_use]
pub fn adler32(bytes: &[u8]) -> u32 {
    let mut state = Adler32State::new();
    state.update(bytes);
    state.finish()
}

#[derive(Debug, Clone, Copy)]
struct Adler32State {
    a: u32,
    b: u32,
}

impl Adler32State {
    const fn new() -> Self {
        Self { a: 1, b: 0 }
    }

    fn update(&mut self, bytes: &[u8]) {
        const MOD: u32 = 65_521;
        for chunk in bytes.chunks(5000) {
            for &x in chunk {
                self.a += u32::from(x);
                self.b += self.a;
            }
            self.a %= MOD;
            self.b %= MOD;
        }
    }

    const fn finish(self) -> u32 {
        (self.b << 16) | self.a
    }
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
    let header = u16::from_be_bytes([z[0], z[1]]);
    if !header.is_multiple_of(31) {
        return Err(ImgError::Malformed {
            what: "zlib header check bits mismatch".to_string(),
        });
    }
    if z[1] & 0x20 != 0 {
        return Err(ImgError::Unsupported {
            what: "zlib preset dictionaries".to_string(),
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
    if pos + 4 != z.len() {
        return Err(ImgError::Malformed {
            what: "trailing bytes after zlib adler32".to_string(),
        });
    }
    Ok(out)
}

fn push_chunk(out: &mut Vec<u8>, kind: [u8; 4], data: &[u8]) {
    assert!(
        u32::try_from(data.len()).is_ok(),
        "PNG chunk exceeds u32 length field"
    );
    out.extend_from_slice(&(data.len() as u32).to_be_bytes());
    out.extend_from_slice(&kind);
    out.extend_from_slice(data);
    let mut crc_input = Vec::with_capacity(4 + data.len());
    crc_input.extend_from_slice(&kind);
    crc_input.extend_from_slice(data);
    out.extend_from_slice(&crc32(&crc_input).to_be_bytes());
}

const SIGNATURE: [u8; 8] = [0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A];

fn checked_sample_count(
    width: u32,
    height: u32,
    channels: usize,
    got: usize,
    context: &'static str,
) -> Result<usize, ImgError> {
    if width == 0 || height == 0 {
        return Err(ImgError::Shape {
            expected: 0,
            got,
            context,
        });
    }
    (width as usize)
        .checked_mul(height as usize)
        .and_then(|px| px.checked_mul(channels))
        .ok_or(ImgError::Shape {
            expected: usize::MAX,
            got,
            context,
        })
}

fn checked_scanline_capacity(
    row_bytes: usize,
    height: u32,
    context: &'static str,
) -> Result<usize, ImgError> {
    row_bytes
        .checked_add(1)
        .and_then(|row| row.checked_mul(height as usize))
        .ok_or(ImgError::Shape {
            expected: usize::MAX,
            got: 0,
            context,
        })
}

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
    let channels = color.channels();
    let expected =
        checked_sample_count(width, height, channels, samples.len(), "write_png8 samples")?;
    if samples.len() != expected {
        return Err(ImgError::Shape {
            expected,
            got: samples.len(),
            context: "write_png8 samples",
        });
    }
    let row = width as usize * channels;
    let mut raw = Vec::with_capacity(checked_scanline_capacity(
        row,
        height,
        "write_png8 scanlines",
    )?);
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
    let channels = color.channels();
    let expected = checked_sample_count(
        width,
        height,
        channels,
        samples.len(),
        "write_png16 samples",
    )?;
    if samples.len() != expected {
        return Err(ImgError::Shape {
            expected,
            got: samples.len(),
            context: "write_png16 samples",
        });
    }
    let row = width as usize * channels;
    let row_bytes = row.checked_mul(2).ok_or(ImgError::Shape {
        expected: usize::MAX,
        got: samples.len(),
        context: "write_png16 scanlines",
    })?;
    let mut raw = Vec::with_capacity(checked_scanline_capacity(
        row_bytes,
        height,
        "write_png16 scanlines",
    )?);
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
    push_chunk(&mut out, *b"IHDR", &ihdr);
    push_chunk(&mut out, *b"sRGB", &[0]); // perceptual intent
    push_chunk(&mut out, *b"IDAT", &zlib_stored(raw));
    push_chunk(&mut out, *b"IEND", &[]);
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

/// Caller-owned byte ceilings for structural PNG inspection.
///
/// Inspection retains no decoded pixels. `max_decoded_bytes` bounds the
/// logical interleaved sample bytes described by IHDR, excluding the one-byte
/// None-filter prefix on each encoded scanline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PngInspectLimits {
    /// Maximum encoded artifact bytes accepted at entry.
    pub max_input_bytes: u64,
    /// Maximum logical decoded interleaved sample bytes.
    pub max_decoded_bytes: u64,
}

impl PngInspectLimits {
    /// Unlimited ceilings for callers that already enforce an outer budget.
    pub const UNBOUNDED: Self = Self {
        max_input_bytes: u64::MAX,
        max_decoded_bytes: u64::MAX,
    };
}

/// Structural facts for a PNG in the exact subset emitted by this crate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PngInspection {
    /// Pixel width.
    pub width: u32,
    /// Pixel height.
    pub height: u32,
    /// Sample bit depth, 8 or 16.
    pub depth: u8,
    /// Packed color layout.
    pub color: PngColor,
    /// Complete encoded artifact length.
    pub input_bytes: u64,
    /// Interleaved sample bytes in one decoded row.
    pub scanline_bytes: u64,
    /// Total logical decoded interleaved sample bytes.
    pub decoded_bytes: u64,
    /// Total uncompressed zlib payload bytes, including row filter prefixes.
    pub filtered_bytes: u64,
}

impl DecodedPng {
    /// 16-bit samples (only valid when `depth == 16`).
    #[must_use]
    pub fn samples16(&self) -> Vec<u16> {
        assert_eq!(self.depth, 16, "samples16 requires a 16-bit PNG");
        let (samples, remainder) = self.bytes.as_chunks::<2>();
        assert!(remainder.is_empty(), "16-bit PNG payload must be even");
        samples.iter().map(|&p| u16::from_be_bytes(p)).collect()
    }
}

fn checked_range<'a>(
    bytes: &'a [u8],
    start: usize,
    len: usize,
    what: &'static str,
) -> Result<&'a [u8], ImgError> {
    let end = start.checked_add(len).ok_or_else(|| ImgError::Malformed {
        what: what.to_string(),
    })?;
    bytes.get(start..end).ok_or_else(|| ImgError::Malformed {
        what: what.to_string(),
    })
}

fn read_chunk(bytes: &[u8], pos: usize) -> Result<(&[u8], &[u8], usize), ImgError> {
    let len_bytes = checked_range(bytes, pos, 8, "truncated chunk header")?;
    let len = u32::from_be_bytes([len_bytes[0], len_bytes[1], len_bytes[2], len_bytes[3]]) as usize;
    let kind = &len_bytes[4..8];
    let data_start = pos.checked_add(8).ok_or_else(|| ImgError::Malformed {
        what: "chunk offset overflow".to_string(),
    })?;
    let data = checked_range(bytes, data_start, len, "truncated chunk data")?;
    let crc_start = data_start
        .checked_add(len)
        .ok_or_else(|| ImgError::Malformed {
            what: "chunk length overflow".to_string(),
        })?;
    let crc = checked_range(bytes, crc_start, 4, "truncated chunk crc")?;
    let mut crc_input = Vec::with_capacity(4 + len);
    crc_input.extend_from_slice(kind);
    crc_input.extend_from_slice(data);
    if u32::from_be_bytes([crc[0], crc[1], crc[2], crc[3]]) != crc32(&crc_input) {
        return Err(ImgError::Malformed {
            what: "chunk crc mismatch".to_string(),
        });
    }
    let next_pos = crc_start
        .checked_add(4)
        .ok_or_else(|| ImgError::Malformed {
            what: "chunk offset overflow".to_string(),
        })?;
    Ok((kind, data, next_pos))
}

fn parse_ihdr(data: &[u8]) -> Result<(u32, u32, u8, PngColor), ImgError> {
    if data.len() != 13 {
        return Err(ImgError::Malformed {
            what: "IHDR length".to_string(),
        });
    }
    let w = u32::from_be_bytes([data[0], data[1], data[2], data[3]]);
    let h = u32::from_be_bytes([data[4], data[5], data[6], data[7]]);
    if w == 0 || h == 0 {
        return Err(ImgError::Malformed {
            what: "PNG dimensions must be nonzero".to_string(),
        });
    }
    let depth = data[8];
    let color = PngColor::from_type_byte(data[9]).ok_or_else(|| ImgError::Unsupported {
        what: format!("color type {}", data[9]),
    })?;
    if depth != 8 && depth != 16 {
        return Err(ImgError::Unsupported {
            what: format!("bit depth {depth}"),
        });
    }
    if data[10] != 0 {
        return Err(ImgError::Unsupported {
            what: format!("compression method {}", data[10]),
        });
    }
    if data[11] != 0 {
        return Err(ImgError::Unsupported {
            what: format!("filter method {}", data[11]),
        });
    }
    if data[12] != 0 {
        return Err(ImgError::Unsupported {
            what: "interlacing".to_string(),
        });
    }
    Ok((w, h, depth, color))
}

#[derive(Debug, Clone, Copy)]
struct PngChunkRef<'a> {
    kind: &'a [u8],
    data: &'a [u8],
    expected_crc: u32,
    next_pos: usize,
}

fn png_inspect_continue(
    poll: &mut impl FnMut() -> bool,
    operation: &'static str,
) -> Result<(), ImgError> {
    if poll() {
        Ok(())
    } else {
        Err(ImgError::Cancelled { operation })
    }
}

fn admit_png_inspect_bytes(
    resource: &'static str,
    requested: u64,
    limit: u64,
) -> Result<(), ImgError> {
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

fn read_chunk_ref(bytes: &[u8], pos: usize) -> Result<PngChunkRef<'_>, ImgError> {
    let header = checked_range(bytes, pos, 8, "truncated chunk header")?;
    let data_len = u32::from_be_bytes(header[..4].try_into().expect("4 bytes")) as usize;
    let kind = &header[4..8];
    let data_start = pos.checked_add(8).ok_or_else(|| ImgError::Malformed {
        what: "PNG chunk offset overflow".to_string(),
    })?;
    let data = checked_range(bytes, data_start, data_len, "truncated chunk data")?;
    let crc_start = data_start
        .checked_add(data_len)
        .ok_or_else(|| ImgError::Malformed {
            what: "PNG chunk length overflow".to_string(),
        })?;
    let crc_bytes = checked_range(bytes, crc_start, 4, "truncated chunk crc")?;
    let next_pos = crc_start
        .checked_add(4)
        .ok_or_else(|| ImgError::Malformed {
            what: "PNG chunk offset overflow".to_string(),
        })?;
    Ok(PngChunkRef {
        kind,
        data,
        expected_crc: u32::from_be_bytes(crc_bytes.try_into().expect("4 bytes")),
        next_pos,
    })
}

fn require_png_chunk(chunk: PngChunkRef<'_>, expected: [u8; 4]) -> Result<(), ImgError> {
    if chunk.kind == expected {
        Ok(())
    } else {
        Err(ImgError::Malformed {
            what: format!(
                "canonical PNG expected chunk {:?}, found {:?}",
                String::from_utf8_lossy(&expected),
                String::from_utf8_lossy(chunk.kind)
            ),
        })
    }
}

fn verify_chunk_crc_with_poll(
    chunk: PngChunkRef<'_>,
    poll: &mut impl FnMut() -> bool,
) -> Result<(), ImgError> {
    const POLL_BYTES: usize = 64 * 1024;
    png_inspect_continue(poll, "PNG structural inspection")?;
    let mut crc = crc32_update(0xFFFF_FFFF, chunk.kind);
    for bytes in chunk.data.chunks(POLL_BYTES) {
        png_inspect_continue(poll, "PNG structural inspection")?;
        crc = crc32_update(crc, bytes);
    }
    let actual = !crc;
    if actual != chunk.expected_crc {
        return Err(ImgError::Malformed {
            what: format!(
                "PNG {:?} chunk CRC mismatch",
                String::from_utf8_lossy(chunk.kind)
            ),
        });
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn inspect_canonical_zlib(
    zlib: &[u8],
    expected_filtered_bytes: u64,
    row_stride: u64,
    poll: &mut impl FnMut() -> bool,
) -> Result<(), ImgError> {
    const POLL_BYTES: usize = 64 * 1024;
    if take_png(zlib, 0, 2)? != [0x78, 0x01] {
        return Err(ImgError::Unsupported {
            what: "PNG IDAT must use the canonical 0x78 0x01 zlib header".to_string(),
        });
    }
    let mut pos = 2usize;
    let mut raw_offset = 0u64;
    let mut adler = Adler32State::new();
    loop {
        png_inspect_continue(poll, "PNG structural inspection")?;
        let header = *zlib.get(pos).ok_or_else(|| ImgError::Malformed {
            what: "truncated deflate block".to_string(),
        })?;
        if header & 0x06 != 0 {
            return Err(ImgError::Unsupported {
                what: "compressed deflate blocks are outside the canonical writer subset"
                    .to_string(),
            });
        }
        if header & 0xF8 != 0 {
            return Err(ImgError::Malformed {
                what: format!("nonzero padding bits in stored-deflate header {header:#04x}"),
            });
        }
        let bfinal = header & 1 == 1;
        let block_header = take_png(zlib, pos + 1, 4)?;
        let len = u16::from_le_bytes(block_header[..2].try_into().expect("2 bytes"));
        let nlen = u16::from_le_bytes(block_header[2..].try_into().expect("2 bytes"));
        if nlen != !len {
            return Err(ImgError::Malformed {
                what: "stored-deflate NLEN mismatch".to_string(),
            });
        }
        let remaining =
            expected_filtered_bytes
                .checked_sub(raw_offset)
                .ok_or(ImgError::SizeOverflow {
                    context: "PNG filtered-byte accounting",
                })?;
        let canonical_len = remaining.min(u64::from(u16::MAX));
        if u64::from(len) != canonical_len {
            return Err(ImgError::Malformed {
                what: format!(
                    "stored-deflate block at raw byte {raw_offset} has length {len}; canonical length is {canonical_len}"
                ),
            });
        }
        if bfinal != (canonical_len == remaining) {
            return Err(ImgError::Malformed {
                what: format!(
                    "stored-deflate final-block marker is noncanonical at raw byte {raw_offset}"
                ),
            });
        }
        let data_start = pos.checked_add(5).ok_or(ImgError::SizeOverflow {
            context: "PNG deflate block offset",
        })?;
        let data = take_png(zlib, data_start, usize::from(len))?;
        let block_end = raw_offset
            .checked_add(u64::from(len))
            .ok_or(ImgError::SizeOverflow {
                context: "PNG filtered-byte accounting",
            })?;

        let remainder = raw_offset % row_stride;
        let mut filter_offset = if remainder == 0 {
            raw_offset
        } else {
            raw_offset
                .checked_add(row_stride - remainder)
                .ok_or(ImgError::SizeOverflow {
                    context: "PNG filter-byte offset",
                })?
        };
        while filter_offset < block_end {
            let in_block = usize::try_from(filter_offset - raw_offset).map_err(|_| {
                ImgError::SizeOverflow {
                    context: "PNG filter-byte offset on this platform",
                }
            })?;
            if data[in_block] != 0 {
                return Err(ImgError::Unsupported {
                    what: format!(
                        "filter type {} on scanline {} (canonical writer emits None)",
                        data[in_block],
                        filter_offset / row_stride
                    ),
                });
            }
            filter_offset =
                filter_offset
                    .checked_add(row_stride)
                    .ok_or(ImgError::SizeOverflow {
                        context: "PNG filter-byte offset",
                    })?;
        }

        for bytes in data.chunks(POLL_BYTES) {
            png_inspect_continue(poll, "PNG structural inspection")?;
            adler.update(bytes);
        }
        raw_offset = block_end;
        pos = data_start
            .checked_add(usize::from(len))
            .ok_or(ImgError::SizeOverflow {
                context: "PNG deflate block offset",
            })?;
        if bfinal {
            break;
        }
    }
    if raw_offset != expected_filtered_bytes {
        return Err(ImgError::Malformed {
            what: format!(
                "PNG filtered payload is {raw_offset} bytes; expected {expected_filtered_bytes}"
            ),
        });
    }
    let trailer = take_png(zlib, pos, 4)?;
    let expected_adler = u32::from_be_bytes(trailer.try_into().expect("4 bytes"));
    if adler.finish() != expected_adler {
        return Err(ImgError::Malformed {
            what: "PNG zlib Adler-32 mismatch".to_string(),
        });
    }
    if pos + 4 != zlib.len() {
        return Err(ImgError::Malformed {
            what: "trailing bytes after PNG zlib Adler-32".to_string(),
        });
    }
    Ok(())
}

fn take_png(bytes: &[u8], pos: usize, len: usize) -> Result<&[u8], ImgError> {
    checked_range(bytes, pos, len, "truncated PNG stored-deflate stream")
}

/// Strictly inspect a PNG emitted by this crate without retaining decoded
/// pixels. This convenience entry point never observes cancellation.
///
/// # Errors
/// Returns [`ImgError`] for budget refusal, malformed or noncanonical
/// structure, or an unsupported PNG feature.
pub fn inspect_png(bytes: &[u8], limits: PngInspectLimits) -> Result<PngInspection, ImgError> {
    inspect_png_with_poll(bytes, limits, || true)
}

/// Strictly inspect a PNG emitted by this crate without retaining decoded
/// pixels, polling before bounded chunk and stored-deflate units. `poll`
/// returns `true` to continue and `false` to cancel.
///
/// The exact canonical writer subset is enforced: `IHDR`, `sRGB`, one
/// `IDAT`, and `IEND` in that order; exact CRCs; canonical stored-deflate
/// framing; None filters; exact Adler-32; and exact EOF.
///
/// # Errors
/// Returns [`ImgError`] for cancellation, budget refusal, malformed or
/// noncanonical structure, or an unsupported PNG feature.
pub fn inspect_png_with_poll(
    bytes: &[u8],
    limits: PngInspectLimits,
    mut poll: impl FnMut() -> bool,
) -> Result<PngInspection, ImgError> {
    png_inspect_continue(&mut poll, "PNG structural inspection")?;
    let input_bytes = u64::try_from(bytes.len()).map_err(|_| ImgError::SizeOverflow {
        context: "PNG input bytes",
    })?;
    admit_png_inspect_bytes("PNG input", input_bytes, limits.max_input_bytes)?;
    if bytes.get(..8) != Some(&SIGNATURE) {
        return Err(ImgError::Malformed {
            what: "missing PNG signature".to_string(),
        });
    }

    let ihdr = read_chunk_ref(bytes, 8)?;
    require_png_chunk(ihdr, *b"IHDR")?;
    verify_chunk_crc_with_poll(ihdr, &mut poll)?;
    let (width, height, depth, color) = parse_ihdr(ihdr.data)?;
    let bytes_per_sample = u64::from(depth / 8);
    let scanline_bytes = u64::from(width)
        .checked_mul(color.channels() as u64)
        .and_then(|bytes| bytes.checked_mul(bytes_per_sample))
        .ok_or(ImgError::SizeOverflow {
            context: "PNG decoded scanline bytes",
        })?;
    let decoded_bytes =
        scanline_bytes
            .checked_mul(u64::from(height))
            .ok_or(ImgError::SizeOverflow {
                context: "PNG decoded sample bytes",
            })?;
    admit_png_inspect_bytes(
        "PNG decoded samples",
        decoded_bytes,
        limits.max_decoded_bytes,
    )?;
    let row_stride = scanline_bytes
        .checked_add(1)
        .ok_or(ImgError::SizeOverflow {
            context: "PNG filtered scanline bytes",
        })?;
    let filtered_bytes =
        row_stride
            .checked_mul(u64::from(height))
            .ok_or(ImgError::SizeOverflow {
                context: "PNG filtered bytes",
            })?;

    let srgb = read_chunk_ref(bytes, ihdr.next_pos)?;
    require_png_chunk(srgb, *b"sRGB")?;
    if srgb.data != [0] {
        return Err(ImgError::Unsupported {
            what: "canonical PNG sRGB rendering intent must be perceptual (0)".to_string(),
        });
    }
    verify_chunk_crc_with_poll(srgb, &mut poll)?;

    let idat = read_chunk_ref(bytes, srgb.next_pos)?;
    require_png_chunk(idat, *b"IDAT")?;
    let iend = read_chunk_ref(bytes, idat.next_pos)?;
    require_png_chunk(iend, *b"IEND")?;
    if !iend.data.is_empty() {
        return Err(ImgError::Malformed {
            what: "PNG IEND length must be zero".to_string(),
        });
    }
    if iend.next_pos != bytes.len() {
        return Err(ImgError::Malformed {
            what: "trailing bytes after canonical PNG IEND".to_string(),
        });
    }

    verify_chunk_crc_with_poll(idat, &mut poll)?;
    inspect_canonical_zlib(idat.data, filtered_bytes, row_stride, &mut poll)?;
    verify_chunk_crc_with_poll(iend, &mut poll)?;
    Ok(PngInspection {
        width,
        height,
        depth,
        color,
        input_bytes,
        scanline_bytes,
        decoded_bytes,
        filtered_bytes,
    })
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
    let mut saw_idat = false;
    loop {
        let (kind, data, next_pos) = read_chunk(bytes, pos)?;
        match kind {
            b"IHDR" => {
                if pos != 8 || header.is_some() {
                    return Err(ImgError::Malformed {
                        what: "IHDR must be the first and only header chunk".to_string(),
                    });
                }
                header = Some(parse_ihdr(data)?);
            }
            b"IDAT" => {
                if header.is_none() {
                    return Err(ImgError::Malformed {
                        what: "IDAT before IHDR".to_string(),
                    });
                }
                saw_idat = true;
                idat.extend_from_slice(data);
            }
            b"IEND" => {
                if data.is_empty() {
                    pos = next_pos;
                    break;
                }
                return Err(ImgError::Malformed {
                    what: "IEND length must be zero".to_string(),
                });
            }
            _ => {
                if header.is_none() {
                    return Err(ImgError::Malformed {
                        what: "ancillary chunk before IHDR".to_string(),
                    });
                }
            }
        }
        pos = next_pos;
    }
    if pos != bytes.len() {
        return Err(ImgError::Malformed {
            what: "trailing bytes after IEND".to_string(),
        });
    }
    if !saw_idat {
        return Err(ImgError::Malformed {
            what: "missing IDAT".to_string(),
        });
    }
    let Some((width, height, depth, color)) = header else {
        return Err(ImgError::Malformed {
            what: "no IHDR before IEND".to_string(),
        });
    };
    let raw = unzlib_stored(&idat)?;
    let bpp = color.channels() * (depth as usize / 8);
    let row = (width as usize).checked_mul(bpp).ok_or(ImgError::Shape {
        expected: usize::MAX,
        got: raw.len(),
        context: "decoded scanline width",
    })?;
    let expected = checked_scanline_capacity(row, height, "decoded scanlines")?;
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

    fn inspected_png_fixture() -> Vec<u8> {
        let (width, height) = (17u32, 5u32);
        let samples = (0..width * height * 4)
            .map(|index| (index.wrapping_mul(257) % 65_521) as u16)
            .collect::<Vec<_>>();
        write_png16(width, height, PngColor::Rgba, &samples).unwrap()
    }

    fn find_chunk_pos(bytes: &[u8], expected: [u8; 4]) -> usize {
        let mut pos = 8usize;
        loop {
            let chunk = read_chunk_ref(bytes, pos).unwrap();
            if chunk.kind == expected {
                return pos;
            }
            assert_ne!(chunk.kind, b"IEND", "chunk not found");
            pos = chunk.next_pos;
        }
    }

    fn rewrite_chunk_crc(bytes: &mut [u8], pos: usize) {
        let len = u32::from_be_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        let crc_start = pos + 8 + len;
        let crc = crc32(&bytes[pos + 4..crc_start]).to_be_bytes();
        bytes[crc_start..crc_start + 4].copy_from_slice(&crc);
    }

    #[test]
    fn strict_png_inspection_retains_no_pixels_and_is_exactly_budgeted() {
        let bytes = inspected_png_fixture();
        let inspected = inspect_png(&bytes, PngInspectLimits::UNBOUNDED).unwrap();
        assert_eq!(
            (
                inspected.width,
                inspected.height,
                inspected.depth,
                inspected.color
            ),
            (17, 5, 16, PngColor::Rgba)
        );
        assert_eq!(inspected.scanline_bytes, 17 * 4 * 2);
        assert_eq!(inspected.decoded_bytes, inspected.scanline_bytes * 5);
        assert_eq!(inspected.filtered_bytes, (inspected.scanline_bytes + 1) * 5);
        assert_eq!(inspected.input_bytes, bytes.len() as u64);

        let exact = PngInspectLimits {
            max_input_bytes: inspected.input_bytes,
            max_decoded_bytes: inspected.decoded_bytes,
        };
        assert!(
            inspect_png(&bytes, exact).is_ok(),
            "equality must be admitted"
        );
        assert!(matches!(
            inspect_png(
                &bytes,
                PngInspectLimits {
                    max_input_bytes: inspected.input_bytes - 1,
                    ..exact
                }
            ),
            Err(ImgError::ResourceLimit {
                resource: "PNG input",
                ..
            })
        ));
        assert!(matches!(
            inspect_png(
                &bytes,
                PngInspectLimits {
                    max_decoded_bytes: inspected.decoded_bytes - 1,
                    ..exact
                }
            ),
            Err(ImgError::ResourceLimit {
                resource: "PNG decoded samples",
                ..
            })
        ));
    }

    #[test]
    fn strict_png_inspection_rejects_noncanonical_blocks_filters_and_eof() {
        let bytes = inspected_png_fixture();
        let idat_pos = find_chunk_pos(&bytes, *b"IDAT");
        let srgb_pos = find_chunk_pos(&bytes, *b"sRGB");
        let idat_data_start = idat_pos + 8;

        let mut wrong_chunk_order = bytes.clone();
        wrong_chunk_order[srgb_pos + 4..srgb_pos + 8].copy_from_slice(b"gAMA");
        rewrite_chunk_crc(&mut wrong_chunk_order, srgb_pos);
        assert!(matches!(
            inspect_png(&wrong_chunk_order, PngInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("expected chunk")
        ));

        let mut bad_block = bytes.clone();
        let len_pos = idat_data_start + 3;
        let len = u16::from_le_bytes(bad_block[len_pos..len_pos + 2].try_into().unwrap());
        let shorter = len - 1;
        bad_block[len_pos..len_pos + 2].copy_from_slice(&shorter.to_le_bytes());
        bad_block[len_pos + 2..len_pos + 4].copy_from_slice(&(!shorter).to_le_bytes());
        rewrite_chunk_crc(&mut bad_block, idat_pos);
        assert!(matches!(
            inspect_png(&bad_block, PngInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("canonical length")
        ));

        let mut bad_filter = bytes.clone();
        let first_raw_byte = idat_data_start + 2 + 5;
        bad_filter[first_raw_byte] = 1;
        rewrite_chunk_crc(&mut bad_filter, idat_pos);
        assert!(matches!(
            inspect_png(&bad_filter, PngInspectLimits::UNBOUNDED),
            Err(ImgError::Unsupported { what }) if what.contains("filter type 1")
        ));

        let mut bad_adler = bytes.clone();
        bad_adler[first_raw_byte + 1] ^= 1;
        rewrite_chunk_crc(&mut bad_adler, idat_pos);
        assert!(matches!(
            inspect_png(&bad_adler, PngInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("Adler-32")
        ));

        let mut trailing = bytes;
        trailing.push(0);
        assert!(matches!(
            inspect_png(&trailing, PngInspectLimits::UNBOUNDED),
            Err(ImgError::Malformed { what }) if what.contains("trailing bytes")
        ));
    }

    #[test]
    fn strict_png_inspection_cancels_at_bounded_units() {
        let samples = vec![7u8; 30_000 * 3];
        let bytes = write_png8(30_000, 1, PngColor::Rgb, &samples).unwrap();
        let mut total_polls = 0usize;
        inspect_png_with_poll(&bytes, PngInspectLimits::UNBOUNDED, || {
            total_polls += 1;
            true
        })
        .unwrap();
        assert!(total_polls > 10, "fixture must cross multiple poll units");

        let cancel_at = total_polls / 2;
        let mut observed = 0usize;
        let result = inspect_png_with_poll(&bytes, PngInspectLimits::UNBOUNDED, || {
            observed += 1;
            observed < cancel_at
        });
        assert!(matches!(
            result,
            Err(ImgError::Cancelled {
                operation: "PNG structural inspection"
            })
        ));
        assert_eq!(observed, cancel_at);
    }

    #[test]
    fn strict_png_inspection_rejects_every_truncated_prefix() {
        let bytes = inspected_png_fixture();
        for end in 0..bytes.len() {
            assert!(
                inspect_png(&bytes[..end], PngInspectLimits::UNBOUNDED).is_err(),
                "truncated prefix of {end} bytes was admitted"
            );
        }
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

    fn rewrite_ihdr_crc(bytes: &mut [u8]) {
        let crc = crc32(&bytes[12..29]).to_be_bytes();
        bytes[29..33].copy_from_slice(&crc);
    }

    #[test]
    fn png_subset_validation_rejects_bad_headers_and_trailing_bytes() {
        let px = vec![1u8; 12];
        let good = write_png8(2, 2, PngColor::Rgb, &px).unwrap();

        let mut bad_compression = good.clone();
        bad_compression[26] = 1;
        rewrite_ihdr_crc(&mut bad_compression);
        assert!(matches!(
            read_png(&bad_compression),
            Err(ImgError::Unsupported { .. })
        ));

        let mut bad_filter = good.clone();
        bad_filter[27] = 1;
        rewrite_ihdr_crc(&mut bad_filter);
        assert!(matches!(
            read_png(&bad_filter),
            Err(ImgError::Unsupported { .. })
        ));

        let mut trailing = good;
        trailing.push(0);
        assert!(matches!(
            read_png(&trailing),
            Err(ImgError::Malformed { .. })
        ));
    }

    #[test]
    fn png_writer_rejects_zero_dimensions() {
        assert!(matches!(
            write_png8(0, 2, PngColor::Rgb, &[]),
            Err(ImgError::Shape { .. })
        ));
        assert!(matches!(
            write_png16(2, 0, PngColor::Gray, &[]),
            Err(ImgError::Shape { .. })
        ));
    }
}
