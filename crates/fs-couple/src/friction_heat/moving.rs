//! Moving, distributed friction heat on the original finite thermal network.
//!
//! Two planar material strips translate along one common tangent. Their original
//! cells retain all heat and conduct it through the existing network. Every
//! intersected cell pair has its own massless common-temperature interface.
//! The overlap area is integrated exactly along a linear path, including every
//! cell-edge crossing; no endpoint sampling, nearest-cell deposition, global
//! isothermal patch, body heat fraction or additional thermal store is used.
//!
//! The footprint has fixed width and length. Flux/temperature/material data are
//! endpoint values over the step, so backward-Euler time refinement is required
//! independently of mesh refinement. Uniform prism faces are explicit geometry,
//! not inferred from contact coefficients. Curved surfaces, changing width,
//! large rotations and nonuniform transverse temperature require a spatially
//! resolved surface owner; this strip model does not silently approximate them.
use super::{FrictionHeatEntropy, FrictionHeatTolerances, entropy::InterfaceEntropy, finite, invalid};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, StepCoefficients,
};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, EnergyAudit, MaterialId,
};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};

/// Material-coordinate partial-slip history and its original moving heat network.
pub mod partial_slip;

/// One fixed plane and a unit tangent; normal points from body A to body B.
#[derive(Clone, Copy, Debug)]
pub struct StripFrame {
    pub origin_m: [f64; 3],
    pub tangent: [f64; 3],
    pub normal: [f64; 3],
    /// Common transverse extent of every face and of the contact footprint.
    pub width_m: f64,
}
/// Actual surface face interval in its body's material coordinates.
#[derive(Clone, Copy, Debug)]
pub struct StripSurfaceCell {
    pub cell: CellId,
    pub span_m: [f64; 2],
}
#[derive(Clone, Copy, Debug)]
pub struct MovingFrictionBudget {
    pub maximum_surface_cells: usize,
    pub maximum_overlap_tests: usize,
    pub maximum_patches: usize,
}
/// Accepted or proposed tangential placement of the unchanged material grids.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingContactPosition {
    pub body_offsets_m: [f64; 2],
    pub footprint_offset_m: f64,
}
/// Linear interpolation between these placements over this thermal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingContactPath {
    pub start: MovingContactPosition,
    pub end: MovingContactPosition,
}
#[derive(Clone, Debug)]
struct BoundCell {
    face: StripSurfaceCell,
    geometry: CellGeometry,
    material: MaterialId,
    mass: f64,
    body: u64,
    distance: f64,
}
/// Original material face geometry; there is no mutable thermal history here.
#[derive(Clone, Debug)]
pub struct MovingFrictionRoute {
    frame: StripFrame,
    surfaces: [Vec<BoundCell>; 2],
    footprint_m: [f64; 2],
    budget: MovingFrictionBudget,
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 256.0 * f64::EPSILON * a.abs().max(b.abs())
}
// One pass over original links; calling exposed_area_m2 for every surface cell
// would turn a resolved moving band into a quadratic mesh traversal.
fn exposed_faces(
    thermal: &DepositionSimulation<'_>,
    ids: impl Iterator<Item = CellId>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, DepositionError> {
    let mut areas = BTreeMap::new();
    for id in ids {
        checkpoint()?;
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(invalid("inactive moving surface cell"))?;
        areas.insert(
            id,
            cell.geometry.surface_area_m2 - cell.boundary.reservoir.map_or(0.0, |r| r.area_m2),
        );
    }
    for link in thermal.links().values() {
        checkpoint()?;
        let covered = thermal
            .contact_patch(link.id)
            .ok_or(DepositionError::Arithmetic)?
            .covered_area_m2;
        for (id, area) in [(link.a, covered[0]), (link.b, covered[1])] {
            if let Some(value) = areas.get_mut(&id) {
                *value = finite(*value - area)?;
            }
        }
    }
    if areas.values().any(|a| *a < 0.0) {
        return Err(invalid("original contacts exceed surface area"));
    }
    Ok(areas)
}
impl MovingFrictionRoute {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        frame: StripFrame,
        mut surfaces: [Vec<StripSurfaceCell>; 2],
        footprint_m: [f64; 2],
        budget: MovingFrictionBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if frame
            .origin_m
            .iter()
            .chain(&frame.tangent)
            .chain(&frame.normal)
            .any(|v| !v.is_finite())
            || !frame.width_m.is_finite()
            || frame.width_m <= 0.0
            || !close(dot(frame.tangent, frame.tangent), 1.0)
            || !close(dot(frame.normal, frame.normal), 1.0)
            || dot(frame.tangent, frame.normal).abs() > 256.0 * f64::EPSILON
            || footprint_m.iter().any(|v| !v.is_finite())
            || footprint_m[1] <= footprint_m[0]
            || budget.maximum_surface_cells == 0
            || budget.maximum_overlap_tests == 0
            || budget.maximum_patches == 0
            || surfaces.iter().any(Vec::is_empty)
            || surfaces[0]
                .len()
                .checked_add(surfaces[1].len())
                .is_none_or(|n| n > budget.maximum_surface_cells)
        {
            return Err(invalid(
                "moving friction requires finite orthonormal strip geometry and bounded nonempty faces",
            ));
        }
        let exposed = exposed_faces(
            thermal,
            surfaces.iter().flatten().map(|f| f.cell),
            &mut checkpoint,
        )?;
        let mut ids = BTreeSet::new();
        let mut bound: [Vec<BoundCell>; 2] = std::array::from_fn(|_| Vec::new());
        for side in 0..2 {
            surfaces[side].sort_by(|a, b| a.span_m[0].total_cmp(&b.span_m[0]));
            for (i, face) in surfaces[side].iter().enumerate() {
                checkpoint()?;
                if face.span_m.iter().any(|v| !v.is_finite())
                    || face.span_m[1] <= face.span_m[0]
                    || !ids.insert(face.cell)
                    || (i != 0 && surfaces[side][i - 1].span_m[1] != face.span_m[0])
                {
                    return Err(invalid(
                        "moving surface faces must be unique, contiguous and nonoverlapping",
                    ));
                }
                let cell = thermal
                    .cells()
                    .get(&face.cell)
                    .ok_or(invalid("inactive moving surface cell"))?;
                if i != 0 && cell.body_id != bound[side][0].body {
                    return Err(invalid(
                        "one moving surface must belong to one original body",
                    ));
                }
                let expected: [f64; 3] = std::array::from_fn(|j| {
                    frame.origin_m[j]
                        + frame.tangent[j] * f64::midpoint(face.span_m[0], face.span_m[1])
                });
                let scale = frame
                    .origin_m
                    .into_iter()
                    .chain(expected)
                    .chain(cell.geometry.start_m)
                    .chain(cell.geometry.end_m)
                    .map(f64::abs)
                    .fold(1.0_f64, f64::max);
                let near = |point: [f64; 3]| {
                    difference(point, expected)
                        .iter()
                        .all(|v| v.abs() <= 256.0 * f64::EPSILON * scale)
                };
                let other = match (near(cell.geometry.start_m), near(cell.geometry.end_m)) {
                    (true, false) => cell.geometry.end_m,
                    (false, true) => cell.geometry.start_m,
                    _ => {
                        return Err(invalid(
                            "surface cell centreline must meet its declared face centre once",
                        ));
                    }
                };
                let depth = dot(difference(other, expected), frame.normal);
                let sign = if side == 0 { -1.0 } else { 1.0 };
                let length = finite(sign * depth)?;
                let residual = difference(
                    other,
                    std::array::from_fn(|j| expected[j] + depth * frame.normal[j]),
                );
                let area = finite((face.span_m[1] - face.span_m[0]) * frame.width_m)?;
                if length <= 0.0
                    || residual
                        .iter()
                        .any(|v| v.abs() > 256.0 * f64::EPSILON * scale)
                    || !close(cell.geometry.volume_m3, finite(area * length)?)
                    || exposed[&face.cell] + 256.0 * f64::EPSILON * area < area
                {
                    return Err(invalid(
                        "moving heat requires exposed, normal, uniform prismatic material faces",
                    ));
                }
                bound[side].push(BoundCell {
                    face: *face,
                    geometry: cell.geometry,
                    material: cell.material_id,
                    mass: cell.mass_kg,
                    body: cell.body_id,
                    distance: 0.5 * length,
                });
            }
        }
        if bound[0][0].body == bound[1][0].body {
            return Err(invalid(
                "moving friction surfaces must name distinct original bodies",
            ));
        }
        let result = Self {
            frame,
            surfaces: bound,
            footprint_m,
            budget,
        };
        result.validate(thermal, &mut checkpoint)?;
        if finite(result.footprint_area_m2())? <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        Ok(result)
    }
    pub fn footprint_area_m2(&self) -> f64 {
        (self.footprint_m[1] - self.footprint_m[0]) * self.frame.width_m
    }
    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        let mut faces = BTreeSet::new();
        for face in self.surfaces.iter().flatten() {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&face.face.cell)
                .ok_or(invalid("moving surface cell disappeared"))?;
            if cell.geometry != face.geometry
                || cell.material_id != face.material
                || cell.mass_kg != face.mass
                || cell.body_id != face.body
            {
                return Err(invalid(
                    "moving surface changed original geometry, material, mass or body",
                ));
            }
            faces.insert(face.face.cell);
        }
        for link in thermal.links().values() {
            checkpoint()?;
            if (faces.contains(&link.a) || faces.contains(&link.b))
                && thermal.cells()[&link.a].body_id != thermal.cells()[&link.b].body_id
            {
                return Err(invalid(
                    "moving contact must not duplicate an authored cross-body thermal link",
                ));
            }
        }
        Ok(())
    }
}

/// One cell pair's EXACT mean overlap area over the complete path.
#[derive(Clone, Debug)]
pub struct SweptContactPatch {
    pub cells: [CellId; 2],
    pub mean_area_m2: f64,
    pub distances_m: [f64; 2],
}
/// Linear endpoint values at normalized step time t in [0,1].
fn at(line: [f64; 2], t: f64) -> f64 {
    line[0] + t * (line[1] - line[0])
}
fn overlap(lines: &[[f64; 2]; 6], t: f64) -> f64 {
    let left = at(lines[0], t).max(at(lines[2], t)).max(at(lines[4], t));
    let right = at(lines[1], t).min(at(lines[3], t)).min(at(lines[5], t));
    (right - left).max(0.0)
}
// Max/min of affine edges is piecewise affine. Split at all edge crossings,
// then integrate each positive overlap exactly by the trapezoid rule.
fn swept_overlap(lines: [[f64; 2]; 6]) -> Result<f64, DepositionError> {
    let mut times = vec![0.0, 1.0];
    for i in 0..6 {
        for j in i + 1..6 {
            let d0 = finite(lines[i][0] - lines[j][0])?;
            let d1 = finite(lines[i][1] - lines[j][1])?;
            if d0 != d1 {
                let t = finite(d0 / finite(d0 - d1)?)?;
                if t > 0.0 && t < 1.0 {
                    times.push(t);
                }
            }
        }
    }
    times.sort_by(f64::total_cmp);
    times.dedup();
    let mut sum = 0.0;
    for w in times.windows(2) {
        sum = finite(sum + 0.5 * (w[1] - w[0]) * (overlap(&lines, w[0]) + overlap(&lines, w[1])))?;
    }
    Ok(sum)
}
fn edges(span: [f64; 2], start: f64, end: f64) -> Result<[[f64; 2]; 2], DepositionError> {
    Ok([
        [finite(span[0] + start)?, finite(span[0] + end)?],
        [finite(span[1] + start)?, finite(span[1] + end)?],
    ])
}
fn bounds(edges: [[f64; 2]; 2]) -> [f64; 2] {
    [edges[0][0].min(edges[0][1]), edges[1][0].max(edges[1][1])]
}

/// Complete endpoint source: local original friction work and original material data.
pub struct MovingFrictionInput<P> {
    /// One physical work value per swept patch; never numerical damping.
    pub dissipated_j: Vec<f64>,
    /// Each touched cell's material conductivity at that cell's endpoint state.
    pub conductivity_w_m_k: BTreeMap<CellId, f64>,
    /// Other sources/conduction only; friction is inserted by this owner.
    pub transport: StepCoefficients<P>,
}
#[derive(Clone, Debug)]
pub struct MovingPatchHeat {
    pub cells: [CellId; 2],
    pub mean_area_m2: f64,
    pub interface_temperature_k: f64,
    pub conductances_w_k: [f64; 2],
    pub body_heat_j: [f64; 2],
    pub dissipated_j: f64,
    pub interface_residual_j: f64,
}
#[derive(Clone, Debug)]
pub struct MovingFrictionHeatReport {
    pub patches: Vec<MovingPatchHeat>,
    pub body_heat_j: [f64; 2],
    pub dissipated_j: f64,
    pub external_heat_j: f64,
    pub thermal: EnergyAudit,
    pub total_energy_residual_j: f64,
    /// Candidate placement; retain only after accepting the complete transaction.
    pub position: MovingContactPosition,
    /// Optional entropy of this same trial, accepted only with all energy gates.
    /// Moving coverage removes the original ambient entropy exchange as well as
    /// its heat. Shared cells contribute storage once, not once per patch.
    pub entropy: Option<FrictionHeatEntropy>,
}
pub struct MovingFrictionHeatStep<'s, 'm> {
    thermal: BackwardEulerStep<'s, 'm>,
    patches: Vec<SweptContactPatch>,
    area: f64,
    end: MovingContactPosition,
    covered: BTreeMap<CellId, f64>,
}
impl<'s, 'm> MovingFrictionHeatStep<'s, 'm> {
    pub fn new(
        thermal: &'s DepositionSimulation<'m>,
        target_s: f64,
        route: &MovingFrictionRoute,
        path: MovingContactPath,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        route.validate(thermal, &mut checkpoint)?;
        let exposed = exposed_faces(
            thermal,
            route.surfaces.iter().flatten().map(|f| f.face.cell),
            &mut checkpoint,
        )?;
        for point in [path.start, path.end] {
            if point
                .body_offsets_m
                .iter()
                .chain(std::iter::once(&point.footprint_offset_m))
                .any(|v| !v.is_finite())
            {
                return Err(invalid("finite moving contact placements required"));
            }
            let f = edges(
                route.footprint_m,
                point.footprint_offset_m,
                point.footprint_offset_m,
            )?;
            for side in 0..2 {
                let cells = &route.surfaces[side];
                let lo = finite(cells[0].face.span_m[0] + point.body_offsets_m[side])?;
                let hi = finite(cells.last().unwrap().face.span_m[1] + point.body_offsets_m[side])?;
                if f[0][0] < lo || f[1][0] > hi {
                    return Err(invalid(
                        "the complete moving footprint must remain on both material surfaces",
                    ));
                }
            }
        }
        let foot = edges(
            route.footprint_m,
            path.start.footprint_offset_m,
            path.end.footprint_offset_m,
        )?;
        let fb = bounds(foot);
        let a = route.surfaces[0]
            .iter()
            .map(|c| {
                edges(
                    c.face.span_m,
                    path.start.body_offsets_m[0],
                    path.end.body_offsets_m[0],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let b = route.surfaces[1]
            .iter()
            .map(|c| {
                edges(
                    c.face.span_m,
                    path.start.body_offsets_m[1],
                    path.end.body_offsets_m[1],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut patches = Vec::new();
        let mut tested = 0usize;
        let mut first = 0;
        for (i, ae) in a.into_iter().enumerate() {
            checkpoint()?;
            let ab = bounds(ae);
            if ab[1] <= fb[0] || ab[0] >= fb[1] {
                continue;
            }
            while first < b.len() && bounds(b[first])[1] <= ab[0] {
                first += 1;
            }
            for (j, &be) in b.iter().enumerate().skip(first) {
                let bb = bounds(be);
                if bb[0] >= ab[1] {
                    break;
                }
                if bb[1] <= fb[0] || bb[0] >= fb[1] {
                    continue;
                }
                checkpoint()?;
                tested += 1;
                if tested > route.budget.maximum_overlap_tests {
                    return Err(DepositionError::Budget("moving footprint overlap tests"));
                }
                let area = finite(
                    route.frame.width_m
                        * swept_overlap([ae[0], ae[1], be[0], be[1], foot[0], foot[1]])?,
                )?;
                if area == 0.0 {
                    continue;
                }
                if patches.len() >= route.budget.maximum_patches {
                    return Err(DepositionError::Budget("moving footprint patches"));
                }
                patches.push(SweptContactPatch {
                    cells: [
                        route.surfaces[0][i].face.cell,
                        route.surfaces[1][j].face.cell,
                    ],
                    mean_area_m2: area,
                    distances_m: [route.surfaces[0][i].distance, route.surfaces[1][j].distance],
                });
            }
        }
        let total = patches
            .iter()
            .try_fold(0.0, |s, p| finite(s + p.mean_area_m2))?;
        let area = route.footprint_area_m2();
        if patches.is_empty()
            || (total - area).abs() > 1024.0 * f64::EPSILON * area * (patches.len() as f64).sqrt()
        {
            return Err(invalid(
                "swept contact area failed full-footprint conservation",
            ));
        }
        // Geometry is fixed for this residual. Prepare coverage once instead
        // of rebuilding it at every nonlinear source evaluation.
        let mut covered = BTreeMap::<CellId, f64>::new();
        for patch in &patches {
            checkpoint()?;
            for id in patch.cells {
                let value = covered.entry(id).or_default();
                *value = finite(*value + patch.mean_area_m2)?;
            }
        }
        for (&id, &covered_area) in &covered {
            checkpoint()?;
            if covered_area > exposed[&id] * (1.0 + 1024.0 * f64::EPSILON) {
                return Err(invalid("moving contact overlaps an already covered surface"));
            }
        }
        Ok(Self {
            thermal: BackwardEulerStep::new(thermal, target_s, checkpoint)?,
            patches,
            area,
            end: path.end,
            covered,
        })
    }
    pub fn patches(&self) -> &[SweptContactPatch] {
        &self.patches
    }
    pub fn cell_ids(&self) -> &[CellId] {
        self.thermal.cell_ids()
    }
    pub fn duration_s(&self) -> f64 {
        self.thermal.duration_s()
    }
    pub fn footprint_area_m2(&self) -> f64 {
        self.area
    }

    /// Explicit uniform work-density distribution in space and time. This is
    /// geometric deposition, never a heat split between the two bodies.
    pub fn uniform_work_j(&self, work_j: f64) -> Result<Vec<f64>, DepositionError> {
        if !work_j.is_finite() || work_j < 0.0 {
            return Err(invalid("finite nonnegative friction work required"));
        }
        let mut result = Vec::with_capacity(self.patches.len());
        let mut remaining = work_j;
        for (i, p) in self.patches.iter().enumerate() {
            let value = if i + 1 == self.patches.len() {
                remaining
            } else {
                finite(work_j * (p.mean_area_m2 / self.area))?
            };
            if value < 0.0 {
                return Err(DepositionError::Arithmetic);
            }
            result.push(value);
            remaining = finite(remaining - value)?;
        }
        Ok(result)
    }

    pub fn evaluate<P, E: From<DepositionError>>(
        &self,
        energy_j: &[f64],
        temperatures_k: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            &[SweptContactPatch],
            &[f64],
        ) -> Result<MovingFrictionInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<MovingFrictionHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_j, temperatures_k, resolve, false, checkpoint)
    }

    /// Evaluate the original moving-contact candidate with caloric entropy.
    /// Every swept patch uses its own solved interface temperature. The original
    /// source callback runs once, and coverage cancels ambient heat and entropy
    /// together. No extra solve, interface heat store or body heat split is used.
    /// Every active thermal law must provide its exact entropy integral. Other
    /// physical sources remain unclassified; numerical defects are not heating.
    /// Use on a final candidate to keep entropy work outside nonlinear iterations.
    pub fn evaluate_with_entropy<P, E: From<DepositionError>>(
        &self,
        energy_j: &[f64],
        temperatures_k: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            &[SweptContactPatch],
            &[f64],
        ) -> Result<MovingFrictionInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<MovingFrictionHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_j, temperatures_k, resolve, true, checkpoint)
    }

    fn evaluate_entropy_option<P, E: From<DepositionError>>(
        &self,
        energy_j: &[f64],
        temperatures_k: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            &[SweptContactPatch],
            &[f64],
        ) -> Result<MovingFrictionInput<P>, E>,
        with_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<MovingFrictionHeatTrial<'m, P>, E> {
        if temperatures_k.len() != self.patches.len()
            || temperatures_k.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(invalid(
                "one finite positive interface temperature per swept patch required",
            )
            .into());
        }
        let checkpoint = RefCell::new(checkpoint);
        let check = || (*checkpoint.borrow_mut())();
        let dt = self.duration_s();
        let mut entropy = None;
        let coefficients=|state: &DepositionSimulation<'m>|->Result<StepCoefficients<(P,MovingFrictionHeatReport)>,E> {
            let MovingFrictionInput { dissipated_j,conductivity_w_m_k,transport }=
                resolve(state,&self.patches,temperatures_k)?;
            if dissipated_j.len()!=self.patches.len()
                || dissipated_j.iter().any(|h| !h.is_finite() || *h<0.0)
            { return Err(invalid("one original nonnegative physical work value per swept patch required").into()); }
            for (&id,&k) in &conductivity_w_m_k {
                if !state.cells().contains_key(&id) || !k.is_finite() || k<0.0 {
                    return Err(invalid("moving contact conductivity must resolve active original material cells").into());
                }
            }
            let StepCoefficients { mut sources_w,conductances_w_k,payload }=transport;
            let mut interface_entropy = if with_entropy {
                Some(InterfaceEntropy::new(state, &sources_w, dt, &check)?)
            } else { None };
            // The footprint removes its mean covered area from the ORIGINAL
            // ambient boundary, without persisting a fictitious fixed contact.
            for (&id,&area) in &self.covered {
                check()?;
                let cell=&state.cells()[&id];
                let body_temperature_k = state.temperature_k(id)?;
                let correction=finite(area*cell.boundary.htc_w_m2_k
                    *(body_temperature_k-cell.boundary.ambient_k))?;
                if let Some(entropy) = &mut interface_entropy {
                    entropy.correct_reservoir_transfer(finite(dt * correction)?,
                        body_temperature_k, cell.boundary.ambient_k)?;
                }
                let source=sources_w.entry(id).or_default();
                *source=finite(*source+correction)?;
            }
            let external_heat_j=finite(dt*state.external_step_power_w(&sources_w,&conductances_w_k,&check)?)?;
            let mut report=MovingFrictionHeatReport { patches:Vec::with_capacity(self.patches.len()),
                body_heat_j:[0.0;2],dissipated_j:0.0,external_heat_j,
                thermal:EnergyAudit { incoming_mass_kg:0.0,incoming_energy_j:0.0,
                    external_heat_j:0.0,stored_energy_change_j:0.0,residual_j:0.0,substeps:0,events:0 },
                total_energy_residual_j:0.0,position:self.end,entropy:None };
            for ((patch,&temperature),&work) in self.patches.iter().zip(temperatures_k).zip(&dissipated_j) {
                check()?;
                let mut g=[0.0;2];
                let mut heat=[0.0;2];
                for side in 0..2 {
                    let id=patch.cells[side];
                    let k=*conductivity_w_m_k.get(&id).ok_or(invalid("missing moving surface conductivity"))?;
                    g[side]=finite(patch.mean_area_m2*k/patch.distances_m[side])?;
                    if k>0.0 && g[side]==0.0 { return Err(DepositionError::Arithmetic.into()); }
                    let body_temperature_k = state.temperature_k(id)?;
                    let power=finite(g[side]*(temperature-body_temperature_k))?;
                    heat[side]=finite(dt*power)?;
                    if let Some(entropy) = &mut interface_entropy {
                        entropy.add_transfer(heat[side], body_temperature_k, temperature)?;
                    }
                    let source=sources_w.entry(id).or_default();
                    *source=finite(*source+power)?;
                    report.body_heat_j[side]=finite(report.body_heat_j[side]+heat[side])?;
                }
                if finite(g[0]+g[1])?<=0.0 {
                    return Err(invalid("every moving interface needs a conducting body").into());
                }
                report.dissipated_j=finite(report.dissipated_j+work)?;
                report.patches.push(MovingPatchHeat { cells:patch.cells,mean_area_m2:patch.mean_area_m2,
                    interface_temperature_k:temperature,conductances_w_k:g,body_heat_j:heat,
                    dissipated_j:work,interface_residual_j:finite(heat[0]+heat[1]-work)? });
            }
            entropy = interface_entropy;
            Ok(StepCoefficients { sources_w,conductances_w_k,payload:(payload,report) })
        };
        let trial = if with_entropy {
            self.thermal.evaluate_with_entropy(energy_j, coefficients, &check)?
        } else {
            self.thermal.evaluate(energy_j, coefficients, &check)?
        };
        let mut report = trial.payload().1.clone();
        if let Some(interface) = entropy {
            report.entropy = Some(FrictionHeatEntropy::from_interfaces(
                trial.entropy_audit().ok_or(DepositionError::Arithmetic)?,
                interface,
                report.patches.iter().map(|patch| (patch.dissipated_j,
                    patch.interface_temperature_k, patch.interface_residual_j)),
                &check,
            )?);
            check()?;
        }
        report.thermal = *trial.audit();
        report.total_energy_residual_j = finite(
            report.thermal.stored_energy_change_j - report.external_heat_j - report.dissipated_j,
        )?;
        let mut residual = trial.residual_j().to_vec();
        residual.extend(report.patches.iter().map(|p| p.interface_residual_j));
        Ok(MovingFrictionHeatTrial {
            trial,
            residual,
            report,
        })
    }
}
pub struct MovingFrictionHeatTrial<'m, P> {
    trial: BackwardEulerTrial<'m, (P, MovingFrictionHeatReport)>,
    residual: Vec<f64>,
    report: MovingFrictionHeatReport,
}
impl<'m, P> MovingFrictionHeatTrial<'m, P> {
    pub fn residual_j(&self) -> &[f64] {
        &self.residual
    }
    pub fn payload(&self) -> &P {
        &self.trial.payload().0
    }
    pub fn report(&self) -> &MovingFrictionHeatReport {
        &self.report
    }
    pub fn accept(
        self,
        tolerance: FrictionHeatTolerances,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, P, MovingFrictionHeatReport), DepositionError> {
        checkpoint()?;
        if [
            tolerance.cell_j,
            tolerance.thermal_j,
            tolerance.interface_j,
            tolerance.total_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "finite nonnegative moving-heat tolerances required",
            ));
        }
        for patch in &self.report.patches {
            if patch.interface_residual_j.abs() > tolerance.interface_j {
                return Err(DepositionError::EnergyBalance(patch.interface_residual_j));
            }
        }
        if self.report.total_energy_residual_j.abs() > tolerance.total_j {
            return Err(DepositionError::EnergyBalance(
                self.report.total_energy_residual_j,
            ));
        }
        let (thermal, (payload, _), _) =
            self.trial
                .accept(tolerance.cell_j, tolerance.thermal_j, &mut checkpoint)?;
        checkpoint()?;
        Ok((thermal, payload, self.report))
    }
}
