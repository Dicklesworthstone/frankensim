//! Source-backed radiation on exterior solid faces. The original material,
//! gray-radiation, half-cell exchange and heat-storage owners do the physics.

use std::collections::BTreeSet;
use fs_conduction::radiation::{AmbientRadiationConfig, AmbientRadiationPatch,
    SurfaceEmissivity, SURFACE_EMISSIVITY_PROPERTY};
use fs_exec::CancelGate;
use fs_lbm::conjugate::boundary_exchange::{BoundaryExchange, PreparedBoundaryExchange};
use fs_lbm::conjugate::{ChtError, EnergySolution, ThermalFace, ThermalSetup, VoxelDomain};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_material::thermal::source::ThermalMaterialPoint;
use super::{J, array, budget, fields, integer, positive, quote, required, scalar, string};
use super::source::{Problem, Source, invalid, poll};

pub(super) struct Radiation {
    pub config: AmbientRadiationConfig,
    faces: Vec<Face>,
}
struct Face {
    name: String,
    cell: usize,
    face: usize,
    material: usize,
    surroundings_k: f64,
    convection_h: f64,
    convection_k: f64,
    initial_k: f64,
}

pub(super) struct FaceReport {
    cell: usize,
    face: usize,
    pub surface_temperature_k: f64,
    emissivity: f64,
    pub radiated_w: f64,
    pub convective_w: f64,
    applied_w: f64,
    pub mismatch_w: f64,
    source: ResolvedMaterialStatePoint,
}
impl FaceReport {
    pub fn outward_w(&self) -> f64 { self.radiated_w + self.convective_w }
}

pub(super) fn load(
    root: &J, domain: &VoxelDomain, materials: &[Source], indices: &[usize], thermal: &ThermalSetup,
) -> Result<Option<Radiation>, String> {
    let Some(input) = root.get("radiation") else { return Ok(None); };
    fields(input, &["surfaces", "max_iterations", "temperature_tolerance_k", "power_tolerance_w",
        "relative_tolerance", "relaxation"])?;
    let config = AmbientRadiationConfig {
        max_iterations: budget(input, "max_iterations", 1000)?,
        temperature_tolerance_k: positive(input, "temperature_tolerance_k")?,
        balance_tolerance_w: scalar(input, "power_tolerance_w")?,
        balance_relative_tolerance: scalar(input, "relative_tolerance")?,
        relaxation: positive(input, "relaxation")?,
    };
    config.validate().map_err(|e| e.to_string())?;
    if config.balance_relative_tolerance >= 1.0 { return Err("radiation relative tolerance must be less than one".into()); }
    let mut faces = Vec::new();
    let mut used = BTreeSet::new();
    let mut names = BTreeSet::new();
    for surface in array(input, "surfaces", 64)? {
        fields(surface, &["name", "material_index", "faces", "surroundings_k",
            "initial_surface_temperature_k", "convection"])?;
        let name = string(surface, "name")?.to_owned();
        if !names.insert(name.clone()) { return Err("radiating surface names must be unique".into()); }
        let material = integer(required(surface, "material_index")?)?;
        let source = materials.get(material).ok_or("radiating surface names an unknown material")?;
        // A single immutable claim prevents selecting a different curve on a
        // temperature-domain exit. No surface emissivity default is allowed.
        if source.pack.card().claims_for(SURFACE_EMISSIVITY_PROPERTY).len() != 1 {
            return Err("radiating material requires exactly one total-emissivity claim".into());
        }
        let surroundings_k = positive(surface, "surroundings_k")?;
        let initial_k = positive(surface, "initial_surface_temperature_k")?;
        let (convection_h, convection_k) = match surface.get("convection") {
            None => (0.0, surroundings_k),
            Some(c) => { fields(c, &["h_w_m2_k", "ambient_k"])?;
                (positive(c, "h_w_m2_k")?, positive(c, "ambient_k")?) },
        };
        let assigned = array(surface, "faces", super::MAX_CELLS)?;
        if assigned.is_empty() { return Err("radiating surface must bind at least one face".into()); }
        for pair in assigned {
            let pair = pair.as_array().filter(|p| p.len() == 2).ok_or("radiating face needs [cell, local_face]")?;
            let (cell, face) = (integer(&pair[0])?, integer(&pair[1])?);
            if cell >= domain.cell_count() || face >= 6 || domain.is_fluid(cell) || indices[cell] != material {
                return Err("radiating face must belong to the declared solid material".into());
            }
            // neighbour() is intentionally internal to fs-lbm; exterior
            // geometry is unambiguous from the public Cartesian coordinates.
            let at = domain.coords(cell); let axis = face / 2;
            if (face % 2 == 0 && at[axis] != 0)
                || (face % 2 == 1 && at[axis] + 1 != domain.dims()[axis])
                || thermal.faces[face] != ThermalFace::Adiabatic || !used.insert((cell, face))
            { return Err("radiation needs a unique exterior face delegated by an adiabatic boundary".into()); }
            if faces.len() >= super::MAX_CELLS { return Err("radiation face budget exceeded".into()); }
            faces.push(Face { name: name.clone(), cell, face, material, surroundings_k,
                convection_h, convection_k, initial_k });
        }
    }
    if faces.is_empty() { return Err("radiation requires at least one surface".into()); }
    Ok(Some(Radiation { config, faces }))
}

impl Face {
    fn at(&self, materials: &[Source], temperature: f64)
        -> Result<(AmbientRadiationPatch, ResolvedMaterialStatePoint), ChtError>
    {
        let source = &materials[self.material];
        let point = source.point_at(temperature)?;
        let (emissivity, resolved) = SurfaceEmissivity::from_condition(&self.name,
            ThermalMaterialPoint { material: source.pack.card(), context: &point,
                temperature_axis: &source.temperature_axis, liquid_fraction_axis: None },
            MaterialPropertySelection::SingleClaimOnly,
        ).map_err(|e| invalid("radiation.material", e))?;
        let patch = AmbientRadiationPatch::new(&self.name, emissivity, self.surroundings_k)
            .map_err(|e| invalid("radiation.material", e))?;
        Ok((patch, resolved))
    }
}
impl Radiation {
    pub fn initial_temperatures(&self) -> Vec<f64> { self.faces.iter().map(|f| f.initial_k).collect() }

    pub fn prepare(&self, problem: &Problem, setup: &ThermalSetup,
        conductivity: &[[f64; 3]], temperatures: &[f64], gate: &CancelGate,
    ) -> Result<PreparedBoundaryExchange, ChtError> {
        if temperatures.len() != self.faces.len() { return Err(invalid("radiation.state", "surface state shape mismatch")); }
        let mut boundaries = Vec::with_capacity(self.faces.len());
        for (i, (face, &t)) in self.faces.iter().zip(temperatures).enumerate() {
            if i.is_multiple_of(256) { poll(gate)?; }
            let (patch, _) = face.at(&problem.materials, t)?;
            let radiative = patch.secant_coefficient_w_m2_k(t).map_err(|e| invalid("radiation.secant", e))?;
            // Different physical mechanisms, ONE half-cell conductive path.
            let h = radiative + face.convection_h;
            let ambient = (face.convection_h / h).mul_add(
                face.convection_k - face.surroundings_k, face.surroundings_k);
            boundaries.push(BoundaryExchange { cell: face.cell, face: face.face,
                coefficient_w_m2_k: h, reservoir_temperature_k: ambient });
        }
        PreparedBoundaryExchange::new(&problem.domain, setup, conductivity, &boundaries, gate)
    }

    pub fn audit(&self, problem: &Problem, prepared: &PreparedBoundaryExchange,
        solution: &EnergySolution, driving: &[f64], gate: &CancelGate,
    ) -> Result<(Vec<FaceReport>, f64), ChtError> {
        let samples = prepared.report(solution, gate)?;
        if samples.len() != self.faces.len() || driving.len() != self.faces.len() {
            return Err(invalid("radiation.state", "prepared face count mismatch"));
        }
        let (mut total_defect, mut total_scale) = (0.0, 0.0);
        let mut defect = 0.0_f64;
        let mut reports = Vec::with_capacity(samples.len());
        for (i, ((face, sample), &old)) in self.faces.iter().zip(samples).zip(driving).enumerate() {
            if i.is_multiple_of(256) { poll(gate)?; }
            let t = sample.surface_temperature_k;
            // Re-query at the ACTUAL surface: the frozen coefficient's
            // validity alone is not permission to reuse its old emissivity.
            let (patch, source) = face.at(&problem.materials, t)?;
            let radiated_w = sample.area_m2 * patch.heat_flux_w_m2(t)
                .map_err(|e| invalid("radiation.flux", e))?;
            let convective_w = sample.area_m2 * face.convection_h * (t - face.convection_k);
            let actual = radiated_w + convective_w;
            let mismatch_w = (actual - sample.outward_w).abs();
            let scale = actual.abs().max(sample.outward_w.abs());
            let limit = self.config.balance_tolerance_w + self.config.balance_relative_tolerance * scale;
            if [radiated_w, convective_w, actual, mismatch_w, limit].iter().any(|v| !v.is_finite()) {
                return Err(invalid("radiation.flux", "nonfinite surface heat or tolerance"));
            }
            defect = defect.max(scaled(mismatch_w, limit))
                .max((t - old).abs() / self.config.temperature_tolerance_k);
            total_defect += mismatch_w; total_scale += scale;
            reports.push(FaceReport { cell: face.cell, face: face.face, surface_temperature_k: t,
                emissivity: patch.emissivity().value(), radiated_w, convective_w,
                applied_w: sample.outward_w, mismatch_w, source });
        }
        let limit = self.config.balance_tolerance_w + self.config.balance_relative_tolerance * total_scale;
        if !total_defect.is_finite() || !limit.is_finite() { return Err(invalid("radiation.flux", "aggregate surface heat overflow")); }
        // Absolute local defects cannot cancel in a globally balanced report.
        defect = defect.max(scaled(total_defect, limit));
        poll(gate)?;
        Ok((reports, defect))
    }
}

fn scaled(defect: f64, limit: f64) -> f64 {
    if defect == 0.0 { 0.0 } else { defect / limit }
}

pub(super) fn render(reports: &[FaceReport]) -> String {
    let rows = reports.iter().map(|r| format!(concat!(
        "{{\"cell\":{},\"face\":{},\"surface_temperature_k\":{},\"emissivity\":{},",
        "\"radiated_w\":{},\"convective_w\":{},\"applied_w\":{},\"heat_defect_w\":{},",
        "\"material_state_identity\":{}}}"), r.cell, r.face, r.surface_temperature_k,
        r.emissivity, r.radiated_w, r.convective_w, r.applied_w, r.mismatch_w,
        quote(&r.source.identity().to_hex()))).collect::<Vec<_>>().join(",");
    format!("[{rows}]")
}
