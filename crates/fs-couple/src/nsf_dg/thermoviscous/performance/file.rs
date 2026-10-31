//! Native resolved NSF source -> complete thermoviscous field -> original DG duct.
//!
//! This reader deliberately embeds the existing resolved-source schema instead
//! of cloning it. The embedded source is fully admitted first, then its cold NSF
//! and DG owners are transferred through the transactional cold handoff. The DG
//! load becomes a duct attached to the complete three-dimensional field; the NSF
//! zone couples to a separate declared field port. No source state is rebuilt.
//!
//! First-rung scope is intentionally strict: the embedded direct source must not
//! already own upstream DG, finite-reservoir, reciprocal-mouth or prepared
//! observer histories, because those owners cannot be silently discarded during
//! a topology replacement. The field itself is a two-port volume: one NSF port
//! and the original primary DG port. Additional field ducts and source-owned
//! molecular/wall memories remain subsequent explicit extensions.

use super::NsfThermoviscousDgPerformance;
use crate::cavity_modal::interior::thermoviscous::{ThermoviscousInterior, source_medium};
use crate::cavity_modal::interior::thermoviscous::dg::{ThermoviscousDg, ThermoviscousDgLimits};
use crate::cavity_modal::interior::thermoviscous::dg::load::ThermoviscousDgLoad;
use crate::nsf_dg::performance::file::{ResolvedDgFilePerformance, ResolvedDgInputError};
use crate::nsf_dg::performance::NsfDgColdPerformance;
use crate::nsf_dg::{NsfDgColdParts, NsfDgError};
use crate::pcm_wav::observation::PressureRenderer;
use fs_blake3::{ContentHash, hash_domain};
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_feec::acoustics::thermoviscous::{ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits};
use fs_feec::acoustics::thermoviscous::time::ThermoviscousTimeLimits;
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::NsfGrid3;
use fs_rep_mesh::TetComplex;
use std::collections::BTreeSet;
use std::sync::Arc;
use fs_exec::CancelGate;

pub const RESOLVED_THERMOVISCOUS_DG_PERFORMANCE_SCHEMA: &str =
    "frankensim-resolved-thermoviscous-dg-performance-v1";
pub const MAX_RESOLVED_THERMOVISCOUS_BYTES: usize = 2 * 1024 * 1024;
const MAX_VERTICES: usize = 16_384;
const MAX_TETS: usize = 16_384;
const MAX_FACES: usize = 32_768;

#[derive(Debug)]
pub enum ResolvedThermoviscousInputError {
    Input(&'static str),
    Source(ResolvedDgInputError),
    Build(String),
}
impl core::fmt::Display for ResolvedThermoviscousInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved thermoviscous input: {self:?}")
    }
}
impl std::error::Error for ResolvedThermoviscousInputError {}
impl From<ResolvedDgInputError> for ResolvedThermoviscousInputError {
    fn from(value: ResolvedDgInputError) -> Self { Self::Source(value) }
}
fn build(error: impl core::fmt::Display) -> ResolvedThermoviscousInputError {
    ResolvedThermoviscousInputError::Build(error.to_string())
}
fn bad(what: &'static str) -> ResolvedThermoviscousInputError {
    ResolvedThermoviscousInputError::Input(what)
}
fn canonical(mut face: [u32; 3]) -> [u32; 3] { face.sort_unstable(); face }

#[derive(Clone, Copy, Debug)]
pub struct ResolvedThermoviscousDgInfo {
    pub input_hash: ContentHash,
    pub vertices: usize,
    pub tetrahedra: usize,
    pub wall_faces: usize,
    pub isothermal_faces: usize,
    pub source_port: usize,
    pub primary_duct_port: usize,
    pub band_hz: [f64; 2],
    pub bulk_to_shear_ratio: f64,
    pub sample_rate_hz: u32,
    pub samples: u64,
    pub full_scale_pa: f64,
}

pub struct ResolvedThermoviscousDgInput {
    info: ResolvedThermoviscousDgInfo,
    points: Vec<[f64; 3]>,
    mesh: TetComplex,
    walls: Vec<[u32; 3]>,
    isothermal: Vec<[u32; 3]>,
    ports: Vec<Vec<[u32; 3]>>,
    source: ResolvedDgFilePerformance,
}

struct Reader<'a> { lines: std::str::Lines<'a>, line: usize }
impl<'a> Reader<'a> {
    fn row(&mut self, key: &'static str) -> Result<Vec<&'a str>, ResolvedThermoviscousInputError> {
        self.line += 1;
        let text = self.lines.next().ok_or_else(|| bad("missing field record"))?;
        let mut fields = text.split_ascii_whitespace();
        if fields.next() != Some(key) { return Err(bad(key)); }
        Ok(fields.collect())
    }
}
fn usize_field(s: &str) -> Result<usize, ResolvedThermoviscousInputError> {
    s.parse().map_err(|_| bad("invalid integer field"))
}
fn scalar(s: &str) -> Result<f64, ResolvedThermoviscousInputError> {
    let value: f64 = s.parse().map_err(|_| bad("invalid scalar field"))?;
    if value.is_finite() { Ok(value) } else { Err(bad("field scalars must be finite")) }
}
fn face(fields: &[&str], vertices: usize) -> Result<[u32; 3], ResolvedThermoviscousInputError> {
    if fields.len() != 3 { return Err(bad("field face requires three vertex indices")); }
    let values = [usize_field(fields[0])?, usize_field(fields[1])?, usize_field(fields[2])?];
    if values.iter().any(|&v| v >= vertices) || values[0] == values[1]
        || values[0] == values[2] || values[1] == values[2]
    { return Err(bad("field face requires three distinct in-range vertices")); }
    Ok(canonical(values.map(|v| v as u32)))
}

impl ResolvedThermoviscousDgInput {
    pub fn from_bytes(bytes: &[u8], maximum_block: usize, gate: Arc<CancelGate>)
        -> Result<Self, ResolvedThermoviscousInputError>
    {
        if bytes.len() > MAX_RESOLVED_THERMOVISCOUS_BYTES || gate.is_requested() {
            return Err(bad("resolved thermoviscous source exceeds its byte budget or was cancelled"));
        }
        let text = core::str::from_utf8(bytes).map_err(|_| bad("source must be UTF-8"))?;
        let (schema, rest) = text.split_once('\n').ok_or_else(|| bad("missing source schema"))?;
        if schema != RESOLVED_THERMOVISCOUS_DG_PERFORMANCE_SCHEMA {
            return Err(bad("unsupported resolved thermoviscous source schema"));
        }
        let marker = "\nresolved_source\n";
        let (field_text, source_text) = rest.split_once(marker)
            .ok_or_else(|| bad("field source requires one embedded resolved_source"))?;
        if source_text.contains(marker) { return Err(bad("resolved_source may occur only once")); }
        let source = ResolvedDgFilePerformance::from_bytes(source_text.as_bytes(), maximum_block,
            Arc::clone(&gate))?;
        let source_info = source.info();

        let mut reader = Reader { lines: field_text.lines(), line: 1 };
        let band = reader.row("field_band")?;
        if band.len() != 2 { return Err(bad("field_band requires FMIN FMAX")); }
        let band_hz = [scalar(band[0])?, scalar(band[1])?];
        if band_hz[0] <= 0.0 || band_hz[1] <= band_hz[0] {
            return Err(bad("field band must be positive and increasing"));
        }
        let bulk = reader.row("field_bulk_viscosity")?;
        if bulk.len() != 2 || bulk[1].is_empty() {
            return Err(bad("field_bulk_viscosity requires RATIO SOURCE_TOKEN"));
        }
        let bulk_to_shear_ratio = scalar(bulk[0])?;
        if bulk_to_shear_ratio < 0.0 { return Err(bad("field bulk/shear ratio must be nonnegative")); }

        let vertices = reader.row("field_vertices")?;
        if vertices.len() != 1 { return Err(bad("field_vertices requires a count")); }
        let vertex_count = usize_field(vertices[0])?;
        if !(4..=MAX_VERTICES).contains(&vertex_count) { return Err(bad("field vertex count exceeds its budget")); }
        let mut points = Vec::with_capacity(vertex_count);
        for _ in 0..vertex_count {
            if gate.is_requested() { return Err(bad("field assembly cancelled")); }
            let row = reader.row("field_vertex")?;
            if row.len() != 3 { return Err(bad("field_vertex requires X Y Z in metres")); }
            points.push([scalar(row[0])?, scalar(row[1])?, scalar(row[2])?]);
        }
        let tets = reader.row("field_tets")?;
        if tets.len() != 1 { return Err(bad("field_tets requires a count")); }
        let tet_count = usize_field(tets[0])?;
        if !(1..=MAX_TETS).contains(&tet_count) { return Err(bad("field tetrahedron count exceeds its budget")); }
        let mut tetrahedra = Vec::with_capacity(tet_count);
        let mut unique_tets = BTreeSet::new();
        for _ in 0..tet_count {
            let row = reader.row("field_tet")?;
            if row.len() != 4 { return Err(bad("field_tet requires four vertex indices")); }
            let mut tet = [0u32; 4];
            for (slot, value) in tet.iter_mut().zip(row) {
                let index = usize_field(value)?;
                if index >= vertex_count { return Err(bad("field tetrahedron vertex lies outside the mesh")); }
                *slot = index as u32;
            }
            let mut identity = tet; identity.sort_unstable();
            if identity.windows(2).any(|x| x[0] == x[1]) || !unique_tets.insert(identity) {
                return Err(bad("field tetrahedra require distinct vertices and unique cells"));
            }
            tetrahedra.push(tet);
        }
        let mesh = TetComplex::from_tets(vertex_count, tetrahedra);

        let walls_row = reader.row("field_walls")?;
        if walls_row.len() != 1 { return Err(bad("field_walls requires a count")); }
        let wall_count = usize_field(walls_row[0])?;
        if wall_count > MAX_FACES { return Err(bad("field wall count exceeds its budget")); }
        let mut walls = Vec::with_capacity(wall_count);
        let mut isothermal = Vec::new();
        let mut owned = BTreeSet::new();
        for _ in 0..wall_count {
            let row = reader.row("field_wall")?;
            if row.len() != 4 { return Err(bad("field_wall requires adiabatic|isothermal I J K")); }
            let boundary = face(&row[1..], vertex_count)?;
            if !owned.insert(boundary) { return Err(bad("field boundary face assigned more than once")); }
            match row[0] {
                "adiabatic" => {},
                "isothermal" => isothermal.push(boundary),
                _ => return Err(bad("field wall requires adiabatic or isothermal")),
            }
            walls.push(boundary);
        }

        let ports_row = reader.row("field_ports")?;
        if ports_row.len() != 1 { return Err(bad("field_ports requires a count")); }
        let port_count = usize_field(ports_row[0])?;
        if !(2..=32).contains(&port_count) { return Err(bad("field needs 2..32 pressure ports")); }
        let mut ports = Vec::with_capacity(port_count);
        for expected in 0..port_count {
            let row = reader.row("field_port")?;
            if row.len() != 2 || usize_field(row[0])? != expected {
                return Err(bad("field ports must be declared in consecutive zero-based order"));
            }
            let count = usize_field(row[1])?;
            if count == 0 || count > MAX_FACES { return Err(bad("field port requires a bounded nonempty face list")); }
            let mut selected = Vec::with_capacity(count);
            for _ in 0..count {
                let row = reader.row("field_face")?;
                let boundary = face(&row, vertex_count)?;
                if !owned.insert(boundary) { return Err(bad("field wall/port faces must be disjoint")); }
                selected.push(boundary);
            }
            ports.push(selected);
        }
        let connect = reader.row("field_connect")?;
        if connect.len() != 2 { return Err(bad("field_connect requires SOURCE_PORT PRIMARY_DG_PORT")); }
        let source_port = usize_field(connect[0])?;
        let primary_duct_port = usize_field(connect[1])?;
        if source_port >= port_count || primary_duct_port >= port_count || source_port == primary_duct_port {
            return Err(bad("field source and primary duct must be distinct declared ports"));
        }
        if reader.lines.next().is_some() { return Err(bad("unexpected field records before resolved_source")); }

        Ok(Self {
            info: ResolvedThermoviscousDgInfo {
                input_hash: hash_domain("org.frankensim.fs-couple.resolved-thermoviscous-dg-performance.v1", bytes),
                vertices: vertex_count, tetrahedra: tet_count, wall_faces: wall_count,
                isothermal_faces: isothermal.len(), source_port, primary_duct_port,
                band_hz, bulk_to_shear_ratio,
                sample_rate_hz: source_info.sample_rate_hz,
                samples: source_info.samples,
                full_scale_pa: source_info.full_scale_pa,
            }, points, mesh, walls, isothermal, ports, source,
        })
    }

    pub const fn info(&self) -> ResolvedThermoviscousDgInfo { self.info }

    /// Build the complete borrowed field only for the duration of the callback.
    /// This avoids self-referential storage and keeps the original source/mesh
    /// lifetimes honest. The callback receives the ordinary pressure renderer.
    pub fn with_renderer<T>(self,
        callback: impl FnOnce(&mut dyn PressureRenderer, ResolvedThermoviscousDgInfo)
            -> Result<T, ResolvedThermoviscousInputError>)
        -> Result<T, ResolvedThermoviscousInputError>
    {
        let info = self.info;
        let cold = self.source.into_cold_performance().map_err(build)?;
        let NsfDgColdPerformance {
            system: NsfDgColdParts { zone, load: primary, limits, high }, boundary,
            observation, sample_rate_hz, output_samples, maximum_block, gate,
        } = cold;
        validate_source_geometry(&zone, high, &self.points, &self.ports[info.source_port])?;
        if gate.is_requested() { return Err(bad("field assembly cancelled")); }
        let gas = *primary.duct().inlet_gas();
        let medium = source_medium(&gas, info.bulk_to_shear_ratio * gas.dynamic_viscosity)
            .map_err(build)?;
        let boundary_geometry = ThermoviscousBoundary {
            no_slip_faces: self.walls.clone(), isothermal_faces: self.isothermal.clone(),
            ports: self.ports.iter().map(|faces| BoundaryMode::uniform(faces)).collect(),
        };
        let field = ThermoviscousFe::new(&self.mesh, &self.points, medium, &boundary_geometry,
            ThermoviscousLimits::default(), &mut || gate.is_requested()).map_err(build)?;
        let source = ThermoviscousInterior::from_gas(field, gas, info.band_hz).map_err(build)?;
        let prepared = source.prepare_time(primary.propagator().dt(), ThermoviscousTimeLimits::default(),
            &mut || gate.is_requested()).map_err(build)?;
        let acoustic = ThermoviscousDg::new(&source, prepared, primary, info.primary_duct_port,
            ThermoviscousDgLimits { maximum_steps: boundary.len() as u64,
                ..ThermoviscousDgLimits::default() }, &gate).map_err(build)?;
        let load = ThermoviscousDgLoad::new(acoustic, info.source_port, &gate).map_err(build)?;
        let downstream_patch = zone.characteristic_patch_on(high).is_some();
        let source_patch = zone.characteristic_patch_on(!high).is_some();
        let system = match (downstream_patch, source_patch) {
            (false, false) => crate::nsf_dg::thermoviscous::NsfThermoviscousDg::new(zone, load, high, limits, &gate),
            (true, false) => crate::nsf_dg::thermoviscous::NsfThermoviscousDg::new_on_patch(zone, load, high, limits, &gate),
            (false, true) => crate::nsf_dg::thermoviscous::NsfThermoviscousDg::new_with_source_patch(zone, load, high, limits, &gate),
            (true, true) => crate::nsf_dg::thermoviscous::NsfThermoviscousDg::new_on_patches(zone, load, high, limits, &gate),
        }.map_err(build)?;
        let mut renderer = NsfThermoviscousDgPerformance::new_resampled_with_drives(
            system, info.primary_duct_port, observation, sample_rate_hz, boundary,
            output_samples, maximum_block, Arc::clone(&gate),
        ).map_err(build)?;
        callback(&mut renderer, info)
    }
}

fn validate_source_geometry(zone: &NsfGrid3<D3q27>, high: bool, points: &[[f64; 3]],
    port: &[[u32; 3]]) -> Result<(), ResolvedThermoviscousInputError>
{
    let open = zone.characteristic_axis().ok_or_else(|| bad("embedded resolved source lost its characteristic axis"))?;
    let traces = if zone.characteristic_patch_on(high).is_some() {
        zone.open_patch_traces_on(high)
    } else { zone.open_boundary_traces(high) }
        .ok_or_else(|| bad("embedded source has no coupled interface traces"))?;
    let first = traces.first().ok_or_else(|| bad("embedded source interface is empty"))?;
    let axis = open.axis;
    let dx = zone.law().dx_m;
    let edge = first.cell[axis];
    if traces.iter().any(|trace| trace.cell[axis] != edge) {
        return Err(bad("resolved interface traces do not lie on one physical plane"));
    }
    let coordinate = if high { (edge + 1) as f64 * dx } else { edge as f64 * dx };
    let tolerance = 256.0 * f64::EPSILON * coordinate.abs().max(dx);
    if points.iter().any(|point| if high { point[axis] < coordinate - tolerance }
        else { point[axis] > coordinate + tolerance })
    {
        return Err(bad("thermoviscous field volume overlaps the resolved NSF volume"));
    }
    let mut vertices = BTreeSet::new();
    for face in port { vertices.extend(face.iter().copied()); }
    if vertices.is_empty() || vertices.iter().any(|&vertex|
        (points[vertex as usize][axis] - coordinate).abs() > tolerance)
    {
        return Err(bad("field source port must lie exactly on the selected NSF interface plane"));
    }
    let tangential: Vec<_> = (0..3).filter(|&a| a != axis).collect();
    for &a in &tangential {
        let low = traces.iter().map(|trace| trace.cell[a]).min().unwrap() as f64 * dx;
        let high_bound = (traces.iter().map(|trace| trace.cell[a]).max().unwrap() + 1) as f64 * dx;
        if vertices.iter().any(|&vertex| {
            let value = points[vertex as usize][a];
            value < low - tolerance || value > high_bound + tolerance
        }) {
            return Err(bad("field source port extends beyond the selected NSF face or patch"));
        }
    }
    Ok(())
}
