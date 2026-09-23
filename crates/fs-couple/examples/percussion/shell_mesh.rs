//! Cold explicit shell-mesh input. All mechanics and radiation use fs-plate's
//! original geometry, section, reduction and extrusion owners.
use super::Error;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::survey::{IsotropicMaterial, MeshShell};
use std::collections::{BTreeMap, BTreeSet};

pub const MAX_BYTES: usize = 8_388_608;
pub const HEADER: &str = "frankensim-shell-mesh-v1";

pub fn parse(text: &str) -> Result<(MeshShell, [f64;2], [f64;2]), Error> {
    if text.len() > MAX_BYTES { return Err("shell mesh exceeds the 8 MiB input limit".into()); }
    let budget = super::super::mesh_budget();
    let mut nodes = BTreeMap::new(); let mut faces = BTreeMap::new();
    let mut materials = BTreeMap::new();
    let mut band = None; let mut strike = None; let mut header = false;
    for (line, raw) in text.lines().enumerate() {
        let data = raw.split('#').next().unwrap_or("").trim();
        if data.is_empty() { continue; }
        if !header {
            if data != HEADER { return Err("shell mesh needs frankensim-shell-mesh-v1 header".into()); }
            header = true; continue;
        }
        let f: Vec<_> = data.split(',').map(str::trim).collect();
        let invalid = || format!("shell mesh line {}: unknown, duplicated or malformed record",line+1);
        let number = |index: usize| -> Result<f64,Error> {
            let x: f64 = f.get(index).ok_or_else(invalid)?.parse().map_err(|_| invalid())?;
            if !x.is_finite() { return Err(invalid().into()); }
            Ok(x)
        };
        let id = |index: usize| -> Result<usize,Error> {
            Ok(f.get(index).ok_or_else(invalid)?.parse::<usize>().map_err(|_| invalid())?)
        };
        match f[0] {
            "band_hz" if f.len() == 3 && band.is_none() => {
                let v = [number(1)?,number(2)?];
                if v[0] <= 0. || v[1] <= v[0] { return Err("shell mesh needs 0 < lower < upper Hz".into()); }
                band = Some(v);
            }
            "strike" if f.len() == 3 && strike.is_none() => { strike = Some([number(1)?,number(2)?]); }
            "material" if f.len() == 5 && materials.len() < budget.max_triangles => {
                let key = id(1)?;
                let m = IsotropicMaterial { young_pa: number(2)?,poisson:number(3)?,density_kg_m3:number(4)? };
                if materials.insert(key,m).is_some() { return Err(invalid().into()); }
            }
            "node" if f.len() == 6 && nodes.len() < budget.max_nodes => {
                let key = id(1)?; let p = [number(2)?,number(3)?,number(4)?]; let h = number(5)?;
                if nodes.insert(key,(p,h)).is_some() { return Err(invalid().into()); }
            }
            "triangle" if f.len() == 6 && faces.len() < budget.max_triangles => {
                let key = id(1)?; let tri = [id(2)?,id(3)?,id(4)?]; let material = id(5)?;
                if faces.insert(key,(tri,material)).is_some() { return Err(invalid().into()); }
            }
            _ => return Err(invalid().into()),
        }
    }
    if !header || nodes.is_empty() || faces.is_empty() || materials.is_empty() {
        return Err("shell mesh requires vertices, triangles and explicit physical materials".into());
    }
    let band = band.ok_or("shell mesh requires band_hz")?;
    let strike = strike.ok_or("shell mesh requires an explicit default strike XY station")?;
    // IDs may be sparse and records reordered. Dense numbering is deterministic
    // by source ID; every position, triangle corner and material assignment survives.
    let mut map = BTreeMap::new(); let mut positions = Vec::with_capacity(nodes.len());
    let mut thickness = Vec::with_capacity(nodes.len());
    for (source,(p,h)) in nodes { map.insert(source,positions.len()); positions.push(p); thickness.push(h); }
    let mut triangles = Vec::with_capacity(faces.len()); let mut assigned = Vec::with_capacity(faces.len());
    let mut used_materials = BTreeSet::new();
    for (_, (corners, material)) in faces {
        let mut tri = [0;3];
        for k in 0..3 { tri[k] = *map.get(&corners[k]).ok_or("shell triangle references a missing vertex")?; }
        triangles.push(tri);
        assigned.push(*materials.get(&material).ok_or("shell triangle references a missing physical material")?);
        used_materials.insert(material);
    }
    if used_materials.len() != materials.len() { return Err("shell mesh has unused material records".into()); }
    let mesh = ShellMesh::new(positions,triangles)?;
    // The percussion host's stick, stand and mute coordinates are vertical XY
    // stations. A locally folded/vertical surface cannot use that chart.
    for e in 0..mesh.tris.len() {
        if mesh.facet(e)?.frame[2][2] <= 1e-6 {
            return Err("percussion shell mesh needs upward, nonvertical facets in its declared XY frame".into());
        }
    }
    let shell = MeshShell::new(mesh,thickness,&assigned,budget)?;
    super::super::playing::shell_location(&shell.mesh.nodes,&shell.mesh.tris,strike)?;
    Ok((shell,band,strike))
}

#[cfg(test)]
#[path = "shell_mesh/tests.rs"]
mod tests;
