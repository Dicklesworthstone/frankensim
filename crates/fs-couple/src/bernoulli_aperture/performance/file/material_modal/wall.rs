//! Source material walls on either original duct; no second wall evolution.
use super::{PlateValveInputError, Reader, ValveSide, bad, checkpoint};
use super::plate::{MaterialCache, relative_path};
use crate::bernoulli_aperture::wall::material::{
    AxialConstraint, MaterialWallAssignment, MaterialWallInput,
};
use crate::bernoulli_aperture::wall::material::source::WallConstitutiveSource;
use fs_exec::CancelGate;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use std::sync::Arc;

/// One solid source at an original tract or bore node, before numerical lowering.
#[derive(Clone, Debug)]
pub struct ModalValveWallRecord {
    side: ValveSide,
    node: usize,
    path: String,
    pack: Arc<NormalizedMaterialCardPack>,
}
impl ModalValveWallRecord {
    /// The original acoustic side; no shared index between unrelated ducts.
    pub const fn side(&self) -> ValveSide { self.side }
    /// Original physical graph node, unchanged by fractional/loss subdivision.
    pub const fn node(&self) -> usize { self.node }
    /// Caller-relative normalized material artifact.
    pub fn path(&self) -> &str { &self.path }
    /// Exact original source, shared with other walls/plates using the same path.
    pub fn pack(&self) -> &NormalizedMaterialCardPack { &self.pack }
}

type Walls = ([Option<MaterialWallInput>; 2], Vec<ModalValveWallRecord>);

pub(super) fn read(rows: &mut Reader<'_>, maximum_frequency_hz: f64,
    cache: &mut MaterialCache,
    load: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    gate: &CancelGate,
) -> Result<Walls, PlateValveInputError> {
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        != Some("walls") {
        return Ok(([None, None], Vec::new()));
    }
    let mut row = rows.row("walls")?;
    if row.fields.clone().next() == Some("none") {
        row.word()?; row.finish()?;
        return Ok(([None, None], Vec::new()));
    }
    let count = row.count(128)?;
    let band_hz = [row.scalar()?, row.scalar()?];
    row.finish()?;
    if count == 0 || band_hz[0] <= 0.0 || band_hz[0] >= band_hz[1]
        || band_hz[1] < maximum_frequency_hz {
        return Err(bad(rows.line, "walls require 1..=128 assignments and a positive ordered band covering propagation and observation"));
    }
    let mut assignments: [Vec<MaterialWallAssignment>; 2] = [Vec::new(), Vec::new()];
    let mut records = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("wall")?;
        let (side, index) = match row.word()? {
            "upstream" => (ValveSide::Upstream, 0),
            "downstream" => (ValveSide::Downstream, 1),
            _ => return Err(bad(row.line, "wall must select upstream or downstream")),
        };
        let node = row.parse()?;
        if assignments[index].len() >= 64 || assignments[index].iter().any(|a| a.node == node) {
            return Err(bad(row.line, "each side admits at most 64 distinct original wall nodes"));
        }
        let orthotropic = match row.word()? {
            "isotropic" => false,
            "orthotropic" => true,
            _ => return Err(bad(row.line, "wall requires an explicit isotropic or principal-aligned orthotropic membrane")),
        };
        let path = row.word()?.to_string();
        relative_path(&path, row.line)?;
        let temperature = row.scalar()?;
        let omega = row.scalar()?;
        let thickness_m = row.scalar()?;
        let axial = match row.word()? {
            "stress-free" => AxialConstraint::StressFree,
            "strain-fixed" => AxialConstraint::StrainFixed,
            _ => return Err(bad(row.line, "wall axial restraint must be stress-free or strain-fixed")),
        };
        let constitutive = WallConstitutiveSource::read(orthotropic, &mut row.fields)
            .map_err(|what| bad(row.line, what))?;
        let line = row.line;
        row.finish()?;
        if temperature <= 0.0 || omega < 0.0 || thickness_m <= 0.0 {
            return Err(bad(line, "wall needs positive solid temperature/thickness and nonnegative source frequency"));
        }
        let pack = cache.load(&path, line, load, gate)?;
        let point = QueryPoint::new().with("T", temperature)
            .and_then(|p| p.with("omega", omega))
            .map_err(|e| PlateValveInputError::MaterialSource { line, source: e.to_string() })?;
        let required = constitutive.requirements().map_err(PlateValveInputError::Physics)?;
        let material = resolve_material_state_point(pack.card(), &point, &required,
            MaterialPropertySelection::SingleClaimOnly)
            .map_err(|source| PlateValveInputError::Material { line, source })?;
        let (model, loss) = constitutive.selection(&material, thickness_m, axial, band_hz)
            .map_err(PlateValveInputError::Physics)?;
        assignments[index].push(MaterialWallAssignment { node, material, model, thickness_m, axial, loss });
        records.push(ModalValveWallRecord { side, node, path, pack });
    }
    Ok((assignments.map(|assignments| (!assignments.is_empty()).then(|| MaterialWallInput {
        max_walls: assignments.len(), assignments, band_hz,
    })), records))
}
