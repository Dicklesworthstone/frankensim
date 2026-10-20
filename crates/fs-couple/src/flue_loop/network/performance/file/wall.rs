//! Addressed material-card ingress for the existing local cylindrical wall law.
//! Source geometry and regional gases stay with the shared duct description.
use super::{FlueInputError, Row, Rows, bad, poll};
use fs_exec::CancelGate;
use crate::bernoulli_aperture::wall::material::{
    AxialConstraint, MaterialWallAssignment, MaterialWallInput, WallMaterialLoss,
    YOUNG_VISCOSITY_PROPERTY,
};
use crate::material_assembly::MaterialSource;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, resolve_material_state_point,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use std::sync::Arc;
use super::material::Sources;

/// File-ingress cap for one complete normalized material-card source.
pub const MAX_FLUE_WALL_PACK_BYTES: usize = super::material::MAX_PACK_BYTES;
const MAX_WALLS: usize = 64;

/// Exact normalized source for one original wall node.
/// Shared paths retain one immutable artifact; state queries remain independent.
#[derive(Clone, Debug)]
pub struct FlueWallRecord {
    node: usize,
    path: String,
    pack: Arc<NormalizedMaterialCardPack>,
}
impl FlueWallRecord {
    /// Original graph node, unchanged by loss or fractional lowering.
    #[must_use] pub const fn node(&self) -> usize { self.node }
    /// Caller-relative artifact path as supplied by the performance data.
    #[must_use] pub fn path(&self) -> &str { &self.path }
    /// Complete normalized source, named material state and source provenance.
    #[must_use] pub fn pack(&self) -> &NormalizedMaterialCardPack { &self.pack }
}

struct ParsedWall {
    node: usize,
    line: usize,
    path: String,
    temperature_k: f64,
    angular_frequency_rad_s: f64,
    thickness_m: f64,
    axial: AxialConstraint,
    loss: WallMaterialLoss,
}
pub(super) struct ParsedWalls {
    line: usize,
    band_hz: [f64; 2],
    assignments: Vec<ParsedWall>,
}
impl ParsedWalls {
    pub fn read(rows: &mut Rows<'_>) -> Result<Option<Self>, FlueInputError> {
        if !rows.next_is("walls") { return Ok(None); }
        let declaration = rows.record("walls")?;
        let line = rows.line;
        if declaration == "none" { return Ok(None); }
        let mut row = Row { fields: declaration.split_ascii_whitespace(), line };
        let count: usize = row.number()?;
        let band_hz = [row.scalar()?, row.scalar()?];
        row.end()?;
        if !(1..=MAX_WALLS).contains(&count) || band_hz[0] <= 0.0 || band_hz[0] >= band_hz[1] {
            return Err(bad(line, "walls needs 1..=64 assignments and a positive ordered frequency band"));
        }
        let mut assignments = Vec::<ParsedWall>::with_capacity(count);
        for _ in 0..count {
            let mut row = rows.row("wall")?;
            let node = row.number()?;
            if assignments.iter().any(|a| a.node == node) {
                return Err(bad(row.line, "wall must name a distinct original graph node"));
            }
            if row.fields.next() != Some("isotropic") {
                return Err(bad(row.line, "wall requires NODE isotropic PACK T OMEGA H AXIAL LOSS"));
            }
            let path = row.fields.next().ok_or_else(|| bad(row.line, "wall needs a material-card pack path"))?;
            // Relative on every supported host. Only the caller opens the file.
            super::material::relative_path(path, row.line)?;
            let temperature_k = row.scalar()?;
            let angular_frequency_rad_s = row.scalar()?;
            let thickness_m = row.scalar()?;
            let axial = match row.fields.next() {
                Some("stress-free") => AxialConstraint::StressFree,
                Some("strain-fixed") => AxialConstraint::StrainFixed,
                _ => return Err(bad(row.line, "wall axial constraint must be stress-free or strain-fixed")),
            };
            let loss = match row.fields.next() {
                Some("elastic") => WallMaterialLoss::Elastic,
                Some("kelvin-voigt") => WallMaterialLoss::KelvinVoigt,
                _ => return Err(bad(row.line, "wall loss must be elastic or sourced kelvin-voigt")),
            };
            let source_line = row.line;
            row.end()?;
            if temperature_k <= 0.0 || angular_frequency_rad_s < 0.0 || thickness_m <= 0.0 {
                return Err(bad(source_line, "wall needs a positive material temperature/thickness and nonnegative angular frequency"));
            }
            assignments.push(ParsedWall { node, line: source_line, path: path.into(),
                temperature_k, angular_frequency_rad_s, thickness_m, axial, loss });
        }
        Ok(Some(Self { line, band_hz, assignments }))
    }

    pub fn resolve(self, maximum_frequency_hz: f64, sources: &mut Sources,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate)
        -> Result<(MaterialWallInput, Vec<FlueWallRecord>), FlueInputError>
    {
        if self.band_hz[1] < maximum_frequency_hz {
            return Err(bad(self.line, "wall material band must cover the fractional propagation band"));
        }
        let mut assignments = Vec::with_capacity(self.assignments.len());
        let mut records = Vec::with_capacity(self.assignments.len());
        for wall in self.assignments {
            poll(gate)?;
            let pack = sources.load(&wall.path, wall.line, loader, gate)?;
            let point = QueryPoint::new().with("T", wall.temperature_k)
                .and_then(|point| point.with("omega", wall.angular_frequency_rad_s))
                .map_err(|source| FlueInputError::MaterialSource { line: wall.line, source: source.to_string() })?;
            let source = MaterialSource { card: pack.card(), point: &point,
                selection: MaterialPropertySelection::SingleClaimOnly };
            let material = (|| {
                let mut required = vec![
                    Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
                    Requirement::try_new("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?,
                    Requirement::try_new("poisson_ratio", Dims::NONE,
                        Admit::OpenInterval { lower: -1.0, upper: 0.5 })?,
                ];
                if matches!(wall.loss, WallMaterialLoss::KelvinVoigt) {
                    required.push(Requirement::try_new(YOUNG_VISCOSITY_PROPERTY,
                        DynViscosity::DIMS, Admit::NonNegative)?);
                }
                resolve_material_state_point(source.card, source.point, &required, source.selection)
            })().map_err(|source| FlueInputError::Material { line: wall.line, source })?;
            assignments.push(MaterialWallAssignment { node: wall.node, material,
                thickness_m: wall.thickness_m, axial: wall.axial, loss: wall.loss });
            records.push(FlueWallRecord { node: wall.node, path: wall.path, pack });
        }
        let max_walls = assignments.len();
        Ok((MaterialWallInput { assignments, band_hz: self.band_hz, max_walls }, records))
    }
}
