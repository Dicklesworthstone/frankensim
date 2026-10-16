//! Export backends beyond the mesh formats: 3MF (additive manufacturing;
//! a minimal STORED-entry ZIP container with the 3D-model XML), glTF as
//! GLB (binary container: JSON chunk + BIN chunk, previews), and legacy
//! VTK (unstructured grid + optional point field, scientific interop).
//! All exports are byte-deterministic.

use fs_rep_mesh::Soup;
use std::fmt::Write as _;

// ------------------------------------------------------------------ zip

/// CRC-32 (IEEE, reflected 0xEDB88320) — the ZIP entry checksum.
fn crc32(bytes: &[u8]) -> u32 {
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

struct ZipEntry {
    name: &'static str,
    offset: u32,
    crc: u32,
    len: u32,
}

/// Append one STORED zip entry; returns bookkeeping for the directory.
fn push_entry(out: &mut Vec<u8>, name: &'static str, data: &[u8]) -> ZipEntry {
    let offset = u32::try_from(out.len()).expect("bounded archive");
    let crc = crc32(data);
    let len = u32::try_from(data.len()).expect("bounded entry");
    out.extend_from_slice(&[0x50, 0x4B, 0x03, 0x04]); // local header
    out.extend_from_slice(&20u16.to_le_bytes()); // version needed
    out.extend_from_slice(&0u16.to_le_bytes()); // flags
    out.extend_from_slice(&0u16.to_le_bytes()); // method: STORED
    out.extend_from_slice(&[0, 0, 0, 0]); // mod time/date (fixed: determinism)
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&len.to_le_bytes()); // compressed
    out.extend_from_slice(&len.to_le_bytes()); // uncompressed
    out.extend_from_slice(&u16::try_from(name.len()).expect("short name").to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes()); // extra len
    out.extend_from_slice(name.as_bytes());
    out.extend_from_slice(data);
    ZipEntry {
        name,
        offset,
        crc,
        len,
    }
}

fn finish_zip(out: &mut Vec<u8>, entries: &[ZipEntry]) {
    let dir_start = u32::try_from(out.len()).expect("bounded archive");
    for e in entries {
        out.extend_from_slice(&[0x50, 0x4B, 0x01, 0x02]); // central header
        out.extend_from_slice(&20u16.to_le_bytes()); // version made by
        out.extend_from_slice(&20u16.to_le_bytes()); // version needed
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes()); // STORED
        out.extend_from_slice(&[0, 0, 0, 0]);
        out.extend_from_slice(&e.crc.to_le_bytes());
        out.extend_from_slice(&e.len.to_le_bytes());
        out.extend_from_slice(&e.len.to_le_bytes());
        out.extend_from_slice(&u16::try_from(e.name.len()).expect("short").to_le_bytes());
        // extra-field-len(2) + comment-len(2) + disk-start(2) + internal-attrs(2)
        // = 8 bytes. Writing 12 made each central-directory record 50 bytes
        // instead of the fixed 46, misaligning every following record's
        // PK\x01\x02 signature and the EOCD's offsets — an invalid ZIP/OPC
        // package no conformant reader can open.
        out.extend_from_slice(&[0u8; 8]);
        out.extend_from_slice(&[0u8; 4]); // external attrs
        out.extend_from_slice(&e.offset.to_le_bytes());
        out.extend_from_slice(e.name.as_bytes());
    }
    let dir_len = u32::try_from(out.len()).expect("bounded") - dir_start;
    out.extend_from_slice(&[0x50, 0x4B, 0x05, 0x06]); // EOCD
    out.extend_from_slice(&[0, 0, 0, 0]); // disk numbers
    let n = u16::try_from(entries.len()).expect("few entries");
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&dir_len.to_le_bytes());
    out.extend_from_slice(&dir_start.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes()); // comment len
}

// ------------------------------------------------------------------ 3mf

/// Export as a minimal 3MF package (OPC ZIP: content types, rels, and
/// the 3D-model XML in millimeter units).
#[must_use]
pub fn export_3mf(soup: &Soup) -> Vec<u8> {
    let content_types = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
        <Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\
        <Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\
        <Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/>\
        </Types>";
    let rels = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
        <Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\
        <Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" \
        Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/>\
        </Relationships>";
    let mut model = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
         <model unit=\"millimeter\" xml:lang=\"en-US\" \
         xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\
         <resources><object id=\"1\" type=\"model\"><mesh><vertices>",
    );
    for p in &soup.positions {
        let _ = write!(model, "<vertex x=\"{}\" y=\"{}\" z=\"{}\"/>", p.x, p.y, p.z);
    }
    model.push_str("</vertices><triangles>");
    for t in &soup.triangles {
        let _ = write!(
            model,
            "<triangle v1=\"{}\" v2=\"{}\" v3=\"{}\"/>",
            t[0], t[1], t[2]
        );
    }
    model.push_str(
        "</triangles></mesh></object></resources>\
         <build><item objectid=\"1\"/></build></model>",
    );
    let mut out = Vec::new();
    let entries = vec![
        push_entry(&mut out, "[Content_Types].xml", content_types.as_bytes()),
        push_entry(&mut out, "_rels/.rels", rels.as_bytes()),
        push_entry(&mut out, "3D/3dmodel.model", model.as_bytes()),
    ];
    finish_zip(&mut out, &entries);
    out
}

// ------------------------------------------------------------------ glb

/// Export as GLB (glTF 2.0 binary container: one mesh, f32 positions +
/// u32 indices; preview-grade by design).
#[must_use]
pub fn export_glb(soup: &Soup) -> Vec<u8> {
    // BIN chunk: positions then indices, 4-byte aligned by construction.
    let mut bin: Vec<u8> = Vec::new();
    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for p in &soup.positions {
        #[allow(clippy::cast_possible_truncation)]
        let v = [p.x as f32, p.y as f32, p.z as f32];
        for (k, x) in v.iter().enumerate() {
            min[k] = min[k].min(*x);
            max[k] = max[k].max(*x);
            bin.extend_from_slice(&x.to_le_bytes());
        }
    }
    let pos_len = bin.len();
    for t in &soup.triangles {
        for &i in t {
            bin.extend_from_slice(&i.to_le_bytes());
        }
    }
    while !bin.len().is_multiple_of(4) {
        bin.push(0);
    }
    let idx_len = soup.triangles.len() * 12;
    let json = format!(
        "{{\"asset\":{{\"version\":\"2.0\",\"generator\":\"fs-io\"}},\
         \"buffers\":[{{\"byteLength\":{}}}],\
         \"bufferViews\":[\
         {{\"buffer\":0,\"byteOffset\":0,\"byteLength\":{pos_len},\"target\":34962}},\
         {{\"buffer\":0,\"byteOffset\":{pos_len},\"byteLength\":{idx_len},\"target\":34963}}],\
         \"accessors\":[\
         {{\"bufferView\":0,\"componentType\":5126,\"count\":{},\"type\":\"VEC3\",\
         \"min\":[{},{},{}],\"max\":[{},{},{}]}},\
         {{\"bufferView\":1,\"componentType\":5125,\"count\":{},\"type\":\"SCALAR\"}}],\
         \"meshes\":[{{\"primitives\":[{{\"attributes\":{{\"POSITION\":0}},\"indices\":1}}]}}],\
         \"nodes\":[{{\"mesh\":0}}],\"scenes\":[{{\"nodes\":[0]}}],\"scene\":0}}",
        bin.len(),
        soup.positions.len(),
        min[0],
        min[1],
        min[2],
        max[0],
        max[1],
        max[2],
        soup.triangles.len() * 3,
    );
    let mut json_bytes = json.into_bytes();
    while !json_bytes.len().is_multiple_of(4) {
        json_bytes.push(b' ');
    }
    let total = 12 + 8 + json_bytes.len() + 8 + bin.len();
    let mut out = Vec::with_capacity(total);
    out.extend_from_slice(b"glTF");
    out.extend_from_slice(&2u32.to_le_bytes());
    out.extend_from_slice(&u32::try_from(total).expect("bounded").to_le_bytes());
    out.extend_from_slice(
        &u32::try_from(json_bytes.len())
            .expect("bounded")
            .to_le_bytes(),
    );
    out.extend_from_slice(b"JSON");
    out.extend_from_slice(&json_bytes);
    out.extend_from_slice(&u32::try_from(bin.len()).expect("bounded").to_le_bytes());
    out.extend_from_slice(b"BIN\0");
    out.extend_from_slice(&bin);
    out
}

// ------------------------------------------------------------------ vtk

/// Export as legacy-ASCII VTK unstructured grid, optionally with one
/// scalar point field (scientific-viz interop).
#[must_use]
pub fn export_vtk(soup: &Soup, field: Option<(&str, &[f64])>) -> String {
    let mut out = String::from("# vtk DataFile Version 3.0\nfs-io export\nASCII\n");
    out.push_str("DATASET UNSTRUCTURED_GRID\n");
    let _ = writeln!(out, "POINTS {} double", soup.positions.len());
    for p in &soup.positions {
        let _ = writeln!(out, "{} {} {}", p.x, p.y, p.z);
    }
    let n = soup.triangles.len();
    let _ = writeln!(out, "CELLS {n} {}", n * 4);
    for t in &soup.triangles {
        let _ = writeln!(out, "3 {} {} {}", t[0], t[1], t[2]);
    }
    let _ = writeln!(out, "CELL_TYPES {n}");
    for _ in 0..n {
        out.push_str("5\n"); // VTK_TRIANGLE
    }
    if let Some((name, values)) = field
        && values.len() == soup.positions.len()
    {
        let _ = writeln!(out, "POINT_DATA {}", values.len());
        let _ = writeln!(out, "SCALARS {name} double 1\nLOOKUP_TABLE default");
        for v in values {
            let _ = writeln!(out, "{v}");
        }
    }
    out
}
