//! A conforming tetrahedral liquid snapshot enters the original sparse VOF.
//!
//! The geometry owner constructs actual tet/cell intersections. Their moments
//! integrate a P1 velocity field; a mass-weighted rigid correction preserves
//! EACH face-connected piece's linear and angular momentum at receiving cell
//! centres. Disconnected fragments cannot exchange impulse or torque through
//! remapping, and aliasing them into one receiving cell requires refinement.
//! Both signs
//! of the resulting kinetic/potential energy change are numerical projection
//! error, bounded by the caller, never heat. The returned report describes
//! ingress only; subsequent evolution belongs to the original VOF solver.
//!
//! The base entry point is an already-liquid, stress-free, uniform-material
//! isothermal snapshot with capillarity explicitly omitted. The `enthalpy`
//! extension carries a P1 caloric field into the original thermal VOF owner.
//! Neither path melts a solid, discards plastic history, or performs a pressure
//! solve.
//! The supplied tetrahedra must form a nonoverlapping mesh. Overfilled target
//! cells refuse, but that check is not a global mesh-intersection certificate.

use std::collections::BTreeMap;

/// Conservative P1 enthalpy ingress into the original moving-liquid owner.
pub mod enthalpy;

mod components;
pub use components::TetLiquidComponentTransfer;

use fs_exec::Cx;
use fs_la::factor::{cholesky, lu};
use fs_mesh::overlap::{OverlapError, tetrahedron_box_moments};

use super::{
    MaterialFreeSurface3, MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3,
    MaterialFreeSurfaceScaling, MaterialLbmError, admit_domain,
};

/// Geometry and the existing P1 velocity field of one already-liquid body.
#[derive(Clone, Copy, Debug)]
pub struct TetrahedralLiquidState<'a> {
    /// Physical vertex positions [m], in the same frame as `grid_origin_m`.
    pub positions_m: &'a [[f64; 3]],
    /// Conforming, nonoverlapping tetrahedra, indexing `positions_m`.
    pub tetrahedra: &'a [[u32; 4]],
    /// Velocity [m/s] at each vertex; interpolated affinely within each tet.
    pub velocities_m_s: &'a [[f64; 3]],
}

/// Work and measured projection-error bounds. No tolerance manufactures heat.
#[derive(Clone, Copy, Debug)]
pub struct TetLiquidTransferLimits {
    /// Maximum number of input tetrahedra.
    pub max_tetrahedra: usize,
    /// Maximum tet/Cartesian-cell intersection evaluations.
    pub max_overlap_tests: usize,
    /// Maximum occupied target cells, also bounded by the native domain cap.
    pub max_wet_cells: usize,
    /// Relative volume partition/overfill roundoff tolerance, in [0, 1e-6].
    pub relative_volume_tolerance: f64,
    /// Sum of absolute per-component mass discrepancies [kg].
    pub mass_tolerance_kg: f64,
    /// Sum of per-component Euclidean linear-momentum discrepancies [kg m/s].
    pub momentum_tolerance_kg_m_s: f64,
    /// Sum of per-component angular-momentum discrepancy norms [kg m²/s],
    /// all about the same grid origin.
    pub angular_momentum_tolerance_kg_m2_s: f64,
    /// Maximum displacement of EACH piece's mass centroid due to cell-centre representation [m].
    pub maximum_centroid_shift_m: f64,
    /// Bound on the SUM over pieces of |change in kinetic energy| +
    /// |change in gravity potential| [J]. Opposite errors cannot cancel.
    pub maximum_projection_energy_change_j: f64,
}

/// Independently integrated source and actual initialized native-cell moments.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LiquidTransferMoments {
    /// Total physical liquid mass [kg].
    pub mass_kg: f64,
    /// First mass moment about the declared grid origin [kg m].
    pub first_mass_moment_kg_m: [f64; 3],
    /// Integrated linear momentum [kg m/s].
    pub momentum_kg_m_s: [f64; 3],
    /// Angular momentum about the declared grid origin [kg m²/s].
    pub angular_momentum_kg_m2_s: [f64; 3],
    /// Resolved kinetic energy [J], without subcell energy invented as heat.
    pub kinetic_energy_j: f64,
}

/// Frozen ingress accounting, not a conservation claim for subsequent steps.
#[derive(Clone, Debug, PartialEq)]
pub struct TetLiquidTransferReport {
    /// Origin used by every reported first/angular moment [m].
    pub grid_origin_m: [f64; 3],
    /// Analytic P1 integrals over the source tetrahedra.
    pub source: LiquidTransferMoments,
    /// Moments read back from the initialized native VOF, at cell centres.
    pub target: LiquidTransferMoments,
    /// Target minus source centroid [m].
    pub centroid_shift_m: [f64; 3],
    /// Signed target minus source kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Signed target minus source potential energy for uniform gravity [J].
    pub gravitational_energy_change_j: f64,
    /// Actual geometric intersection evaluations.
    pub overlap_tests: usize,
    /// Actual occupied target cells.
    pub wet_cells: usize,
    /// Original face-connected source pieces and their actual native cells.
    /// Every piece independently preserves mass and both momenta. These are
    /// ingress lineage and audits, not persistent labels after native merging.
    pub components: Vec<TetLiquidComponentTransfer>,
}

/// A geometric, physical, projection, budget or cancellation refusal.
#[derive(Debug)]
pub enum TetLiquidTransferError {
    /// Invalid source or a projection outside its admitted physical bounds.
    Invalid(&'static str),
    /// The existing geometry owner refused an intersection.
    Geometry(OverlapError),
    /// The existing material/native flow owner refused initialization.
    Flow(MaterialLbmError),
    /// The original execution context observed cancellation.
    Cancelled,
}
impl core::fmt::Display for TetLiquidTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral liquid transfer: {self:?}")
    }
}
impl std::error::Error for TetLiquidTransferError {}
impl From<OverlapError> for TetLiquidTransferError {
    fn from(error: OverlapError) -> Self {
        Self::Geometry(error)
    }
}
impl From<MaterialLbmError> for TetLiquidTransferError {
    fn from(error: MaterialLbmError) -> Self {
        Self::Flow(error)
    }
}

#[derive(Default)]
struct CellIntegral {
    component: usize,
    volume: f64,
    momentum: [f64; 3],
    enthalpy_j: f64,
}

#[derive(Default)]
struct TransferredEnthalpy {
    by_cell: BTreeMap<[usize; 3], f64>,
    source_j: f64,
    source_by_component_j: Vec<f64>,
}

fn invalid(what: &'static str) -> TetLiquidTransferError {
    TetLiquidTransferError::Invalid(what)
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), TetLiquidTransferError> {
    cx.checkpoint()
        .map_err(|_| TetLiquidTransferError::Cancelled)
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    a.map(|v| v * s)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn centre(c: [usize; 3], dx: f64) -> [f64; 3] {
    c.map(|v| (v as f64 + 0.5) * dx)
}
fn finite_moments(m: &LiquidTransferMoments) -> bool {
    m.mass_kg.is_finite()
        && m.mass_kg > 0.0
        && m.kinetic_energy_j.is_finite()
        && m.kinetic_energy_j >= 0.0
        && [
            m.first_mass_moment_kg_m,
            m.momentum_kg_m_s,
            m.angular_momentum_kg_m2_s,
        ]
        .iter()
        .flatten()
        .all(|v| v.is_finite())
}

impl<'a> MaterialFreeSurfaceScaling<'a> {
    /// Transfer actual geometry and P1 motion into this material's native flow.
    ///
    /// Source density is the SAME immutable liquid density used for native mass,
    /// inertia and viscosity scaling. No material name or extra density is accepted.
    /// Source tets must already be liquid and stress-free; this API does not
    /// authorize deleting a solid's internal energy or history. Surface tension
    /// must be explicitly neglected because a conservative surface-energy
    /// projection is not represented here. Initial numerical pressure is zero
    /// relative to the existing gas reference. Pressure initialization is not
    /// an incompressible projection or an initially equilibrated free surface.
    ///
    /// Cancellation and every refusal leave the input snapshot untouched and
    /// return no candidate. Native material speed/Mach gates also apply after
    /// the angular correction. Ill-resolved spin, excessive centroid/energy
    /// error, geometry outside the domain and target overfill require refinement.
    /// Face connectivity uses the original shared vertex indices. Each
    /// disconnected piece is corrected and audited independently; pieces that
    /// occupy the same receiving cell refuse. Source-tet and target-cell lineage
    /// describes ingress only; subsequent coalescence belongs to native VOF.
    pub fn free_surface3_from_tetrahedra(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        grid_origin_m: [f64; 3],
        source: TetrahedralLiquidState<'_>,
        limits: TetLiquidTransferLimits,
        cx: &Cx<'_>,
    ) -> Result<(MaterialFreeSurface3<'a>, TetLiquidTransferReport), TetLiquidTransferError> {
        let (flow, report, _) = self.transfer_tetrahedra(
            domain, grid_origin_m, source, None, limits, cx,
        )?;
        Ok((flow, report))
    }

    // One geometry/momentum owner for isothermal and enthalpy-bearing ingress.
    #[allow(clippy::too_many_arguments)]
    fn transfer_tetrahedra(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        grid_origin_m: [f64; 3],
        source: TetrahedralLiquidState<'_>,
        source_h: Option<&[f64]>,
        limits: TetLiquidTransferLimits,
        cx: &Cx<'_>,
    ) -> Result<(MaterialFreeSurface3<'a>, TetLiquidTransferReport, TransferredEnthalpy), TetLiquidTransferError> {
        checkpoint(cx)?;
        self.require_frozen_source()?;
        admit_domain(&domain, &|_, _, _| false)?;
        if self.material().surface_tension_n_m().is_some() || self.planar_wetting().is_some() {
            return Err(invalid(
                "tet transfer requires explicitly neglected capillarity",
            ));
        }
        admit_limits(limits)?;
        if source.tetrahedra.is_empty()
            || source.tetrahedra.len() > limits.max_tetrahedra
            || source.positions_m.len() < 4
            || source.positions_m.len() != source.velocities_m_s.len()
            || source_h.is_some_and(|h| h.len() != source.positions_m.len())
            || grid_origin_m.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("invalid tetrahedral source or input count limit"));
        }
        let dx = self.cell_size_m();
        let cell_volume = dx * dx * dx;
        let rho = self.material().density_kg_m3();
        let upper: [f64; 3] =
            core::array::from_fn(|i| grid_origin_m[i] + domain.cells[i] as f64 * dx);
        if !cell_volume.is_finite()
            || cell_volume <= 0.0
            || (0..3).any(|i| !upper[i].is_finite() || upper[i] <= grid_origin_m[i])
        {
            return Err(invalid("unrepresentable physical grid extent"));
        }
        let (component_for_tet, mut components) = components::partition(source, cx)?;
        let mut cells: BTreeMap<[usize; 3], CellIntegral> = BTreeMap::new();
        let mut original = LiquidTransferMoments::default();
        let mut thermal = TransferredEnthalpy::default();
        thermal.source_by_component_j.resize(components.len(), 0.0);
        let mut overlap_tests = 0usize;
        for (tet_index, indices) in source.tetrahedra.iter().enumerate() {
            checkpoint(cx)?;
            let component = component_for_tet[tet_index];
            if indices
                .iter()
                .any(|&i| i as usize >= source.positions_m.len())
            {
                return Err(invalid("tetrahedron index lies outside the source field"));
            }
            let p = indices.map(|i| source.positions_m[i as usize]);
            let v = indices.map(|i| source.velocities_m_s[i as usize]);
            if p.iter().any(|x| {
                (0..3).any(|a| !x[a].is_finite() || x[a] < grid_origin_m[a] || x[a] > upper[a])
            }) || v.iter().flatten().any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "source field is nonfinite or outside the receiving grid",
                ));
            }
            let lower_tet: [f64; 3] =
                core::array::from_fn(|a| p.iter().map(|x| x[a]).fold(f64::INFINITY, f64::min));
            let upper_tet: [f64; 3] =
                core::array::from_fn(|a| p.iter().map(|x| x[a]).fold(f64::NEG_INFINITY, f64::max));
            overlap_tests = overlap_tests
                .checked_add(1)
                .filter(|n| *n <= limits.max_overlap_tests)
                .ok_or_else(|| invalid("tet/cell intersection budget exhausted"))?;
            let whole = tetrahedron_box_moments(p, lower_tet, upper_tet, p[0], cx)?;
            if !whole.volume.is_finite() || whole.volume <= 0.0 {
                return Err(invalid(
                    "source tetrahedron has no representable positive volume",
                ));
            }
            integrate_source(&mut original, p, v, whole.volume * rho, grid_origin_m);
            integrate_source(&mut components[component].source, p, v, whole.volume * rho, grid_origin_m);
            let gradient = field_gradient(p, v)?;
            let h = source_h.map(|values| indices.map(|i| values[i as usize]));
            let h_gradient = h.map(|values| field_gradient(p, values.map(|value| [value])))
                .transpose()?;
            if let Some(values) = h {
                let energy = rho * whole.volume
                    * (values.into_iter().map(|value| 0.25 * value).sum::<f64>());
                thermal.source_j += energy;
                thermal.source_by_component_j[component] += energy;
            }
            let first: [usize; 3] = core::array::from_fn(|a| {
                (((lower_tet[a] - grid_origin_m[a]) / dx).floor() as usize).min(domain.cells[a] - 1)
            });
            let end: [usize; 3] = core::array::from_fn(|a| {
                (((upper_tet[a] - grid_origin_m[a]) / dx).ceil() as usize).min(domain.cells[a])
            });
            let mut covered = 0.0;
            for z in first[2]..end[2] {
                for y in first[1]..end[1] {
                    for x in first[0]..end[0] {
                        checkpoint(cx)?;
                        overlap_tests = overlap_tests
                            .checked_add(1)
                            .filter(|n| *n <= limits.max_overlap_tests)
                            .ok_or_else(|| invalid("tet/cell intersection budget exhausted"))?;
                        let c = [x, y, z];
                        let low = core::array::from_fn(|a| grid_origin_m[a] + c[a] as f64 * dx);
                        let high =
                            core::array::from_fn(|a| grid_origin_m[a] + (c[a] + 1) as f64 * dx);
                        let overlap = tetrahedron_box_moments(p, low, high, p[0], cx)?;
                        if overlap.volume == 0.0 {
                            continue;
                        }
                        covered += overlap.volume;
                        if !cells.contains_key(&c)
                            && (cells.len() >= limits.max_wet_cells
                                || cells.len() >= domain.max_cells)
                        {
                            return Err(invalid("occupied target-cell budget exhausted"));
                        }
                        let target = cells.entry(c).or_insert_with(|| CellIntegral {
                            component, ..CellIntegral::default()
                        });
                        if target.component != component {
                            return Err(invalid(
                                "disconnected source pieces share a receiving cell; refine before transfer",
                            ));
                        }
                        target.volume += overlap.volume;
                        if let (Some(values), Some(gradient)) = (h, h_gradient) {
                            target.enthalpy_j += rho
                                * (values[0] * overlap.volume + dot(gradient[0], overlap.first));
                        }
                        for a in 0..3 {
                            target.momentum[a] +=
                                rho * (v[0][a] * overlap.volume + dot(gradient[a], overlap.first));
                        }
                    }
                }
            }
            if !covered.is_finite()
                || (covered - whole.volume).abs() > limits.relative_volume_tolerance * whole.volume
            {
                return Err(invalid(
                    "geometric cell intersections do not partition source volume",
                ));
            }
        }
        if !finite_moments(&original) || cells.is_empty() {
            return Err(invalid(
                "unrepresentable source mass, momentum or kinetic energy",
            ));
        }
        let mut initial = Vec::with_capacity(cells.len());
        let mut component_for_cell = Vec::with_capacity(cells.len());
        for (coordinates, integral) in cells {
            checkpoint(cx)?;
            let fill = integral.volume / cell_volume;
            if !fill.is_finite() || fill <= 0.0 || fill > 1.0 + limits.relative_volume_tolerance {
                return Err(invalid("source mesh overfills a receiving Cartesian cell"));
            }
            // A roundoff overfill is clipped only inside the explicit bound;
            // the actual resulting mass discrepancy is audited below as well.
            let fill_fraction = fill.min(1.0);
            let mass = rho * cell_volume * fill_fraction;
            if source_h.is_some() {
                let h = integral.enthalpy_j / mass;
                if !h.is_finite() || !thermal.source_j.is_finite() {
                    return Err(invalid("transferred specific or extensive enthalpy is not finite"));
                }
                // Use the actual admitted native mass, including bounded
                // geometric roundoff clipping; no thermal energy is discarded.
                thermal.by_cell.insert(coordinates, h);
            }
            initial.push(MaterialFreeSurfaceInitialCell3 {
                coordinates,
                fill_fraction,
                velocity_m_s: scale(integral.momentum, 1.0 / mass),
            });
            components[integral.component].cell_indices.push(initial.len() - 1);
            component_for_cell.push(integral.component);
        }
        components::correct(&mut initial, &components, rho * cell_volume, dx, limits, cx)?;
        let predicted = integrate_cells(&initial, rho * cell_volume, dx, cx)?;
        check_projection(&original, &predicted, domain.gravity_m_s2, limits)?;
        let flow = self.free_surface3_from_cells(domain, &initial, cx)?;
        let mut actual = LiquidTransferMoments::default();
        let mut actual_components = vec![LiquidTransferMoments::default(); components.len()];
        for (index, cell) in initial.iter().enumerate() {
            checkpoint(cx)?;
            let [x, y, z] = cell.coordinates;
            let sample = flow
                .sample_at_si(x, y, z)?
                .ok_or_else(|| invalid("native initializer lost an occupied cell"))?;
            accumulate_cell(
                &mut actual,
                sample.density_kg_m3 * sample.fill_fraction * cell_volume,
                centre(cell.coordinates, dx),
                sample.velocity_m_s,
            );
            accumulate_cell(
                &mut actual_components[component_for_cell[index]],
                sample.density_kg_m3 * sample.fill_fraction * cell_volume,
                centre(cell.coordinates, dx), sample.velocity_m_s,
            );
        }
        if (flow.ledger_mass_kg()? - actual.mass_kg).abs() > limits.mass_tolerance_kg {
            return Err(invalid(
                "native mass ledger and initialized physical cells disagree",
            ));
        }
        let (centroid_shift_m, kinetic_energy_change_j, gravitational_energy_change_j) =
            check_projection(&original, &actual, domain.gravity_m_s2, limits)?;
        let components = components::reports(
            components, actual_components, &initial, domain.gravity_m_s2, limits, cx,
        )?;
        checkpoint(cx)?;
        Ok((
            flow,
            TetLiquidTransferReport {
                grid_origin_m,
                source: original,
                target: actual,
                centroid_shift_m,
                kinetic_energy_change_j,
                gravitational_energy_change_j,
                overlap_tests,
                wet_cells: initial.len(),
                components,
            },
            thermal,
        ))
    }
}

fn admit_limits(limits: TetLiquidTransferLimits) -> Result<(), TetLiquidTransferError> {
    if limits.max_tetrahedra == 0
        || limits.max_overlap_tests == 0
        || limits.max_wet_cells == 0
        || !limits.relative_volume_tolerance.is_finite()
        || !(0.0..=1e-6).contains(&limits.relative_volume_tolerance)
        || [
            limits.mass_tolerance_kg,
            limits.momentum_tolerance_kg_m_s,
            limits.angular_momentum_tolerance_kg_m2_s,
            limits.maximum_centroid_shift_m,
            limits.maximum_projection_energy_change_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(invalid(
            "invalid geometric, momentum, energy or count bounds",
        ));
    }
    Ok(())
}

fn field_gradient<const N: usize>(
    p: [[f64; 3]; 4],
    v: [[f64; N]; 4],
) -> Result<[[f64; 3]; N], TetLiquidTransferError> {
    let rows: [[f64; 3]; 3] = core::array::from_fn(|i| sub(p[i + 1], p[0]));
    let matrix: [f64; 9] = core::array::from_fn(|i| rows[i / 3][i % 3]);
    let factor =
        lu(&matrix, 3).map_err(|_| invalid("source tetrahedron has singular P1 geometry"))?;
    let mut gradient = [[0.0; 3]; N];
    for a in 0..N {
        gradient[a] = core::array::from_fn(|i| v[i + 1][a] - v[0][a]);
        factor.solve(&mut gradient[a]);
    }
    if gradient.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable P1 field gradient"));
    }
    Ok(gradient)
}

fn integrate_source(
    total: &mut LiquidTransferMoments,
    p: [[f64; 3]; 4],
    v: [[f64; 3]; 4],
    mass: f64,
    origin: [f64; 3],
) {
    // Exact P1 simplex integrals: integral(N_i N_j) = V*(1+delta_ij)/20.
    let r = p.map(|x| sub(x, origin));
    let sum_r = r.into_iter().fold([0.0; 3], add);
    let sum_v = v.into_iter().fold([0.0; 3], add);
    let diagonal_cross = (0..4).map(|i| cross(r[i], v[i])).fold([0.0; 3], add);
    total.mass_kg += mass;
    total.first_mass_moment_kg_m = add(total.first_mass_moment_kg_m, scale(sum_r, mass / 4.0));
    total.momentum_kg_m_s = add(total.momentum_kg_m_s, scale(sum_v, mass / 4.0));
    total.angular_momentum_kg_m2_s = add(
        total.angular_momentum_kg_m2_s,
        scale(add(cross(sum_r, sum_v), diagonal_cross), mass / 20.0),
    );
    total.kinetic_energy_j +=
        mass / 40.0 * (dot(sum_v, sum_v) + v.into_iter().map(|x| dot(x, x)).sum::<f64>());
}

fn accumulate_cell(
    total: &mut LiquidTransferMoments,
    mass: f64,
    position: [f64; 3],
    velocity: [f64; 3],
) {
    let momentum = scale(velocity, mass);
    total.mass_kg += mass;
    total.first_mass_moment_kg_m = add(total.first_mass_moment_kg_m, scale(position, mass));
    total.momentum_kg_m_s = add(total.momentum_kg_m_s, momentum);
    total.angular_momentum_kg_m2_s = add(total.angular_momentum_kg_m2_s, cross(position, momentum));
    total.kinetic_energy_j += 0.5 * mass * dot(velocity, velocity);
}
fn integrate_cells(
    cells: &[MaterialFreeSurfaceInitialCell3],
    mass_unit: f64,
    dx: f64,
    cx: &Cx<'_>,
) -> Result<LiquidTransferMoments, TetLiquidTransferError> {
    let mut total = LiquidTransferMoments::default();
    for cell in cells {
        checkpoint(cx)?;
        accumulate_cell(
            &mut total,
            mass_unit * cell.fill_fraction,
            centre(cell.coordinates, dx),
            cell.velocity_m_s,
        );
    }
    Ok(total)
}

fn correct_momenta(
    cells: &mut [MaterialFreeSurfaceInitialCell3],
    mass_unit: f64,
    dx: f64,
    source: &LiquidTransferMoments,
    limits: TetLiquidTransferLimits,
    cx: &Cx<'_>,
) -> Result<(), TetLiquidTransferError> {
    let before = integrate_cells(cells, mass_unit, dx, cx)?;
    if !finite_moments(&before)
        || (before.mass_kg - source.mass_kg).abs() > limits.mass_tolerance_kg
    {
        return Err(invalid(
            "geometric transfer does not preserve physical mass",
        ));
    }
    let centroid = scale(before.first_mass_moment_kg_m, 1.0 / before.mass_kg);
    let delta_v = scale(
        sub(source.momentum_kg_m_s, before.momentum_kg_m_s),
        1.0 / before.mass_kg,
    );
    let mut inertia = [0.0; 9];
    for cell in cells.iter_mut() {
        checkpoint(cx)?;
        cell.velocity_m_s = add(cell.velocity_m_s, delta_v);
        let r = sub(centre(cell.coordinates, dx), centroid);
        let r2 = dot(r, r);
        let mass = cell.fill_fraction * mass_unit;
        for a in 0..3 {
            for b in 0..3 {
                inertia[a * 3 + b] += mass * (if a == b { r2 } else { 0.0 } - r[a] * r[b]);
            }
        }
    }
    let translated = integrate_cells(cells, mass_unit, dx, cx)?;
    let mut omega = sub(
        source.angular_momentum_kg_m2_s,
        translated.angular_momentum_kg_m2_s,
    );
    if norm(omega) > limits.angular_momentum_tolerance_kg_m2_s {
        // Normalize before the native dense solve to retain density
        // similitude and avoid putting the physical mass scale in its pivots.
        for value in &mut inertia {
            *value /= before.mass_kg;
        }
        omega = scale(omega, 1.0 / before.mass_kg);
        let factor = cholesky(&inertia, 3)
            .map_err(|_| invalid("target grid cannot resolve source angular momentum; refine"))?;
        factor.solve(&mut omega);
        if omega.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable angular momentum correction"));
        }
        for cell in cells {
            checkpoint(cx)?;
            cell.velocity_m_s = add(
                cell.velocity_m_s,
                cross(omega, sub(centre(cell.coordinates, dx), centroid)),
            );
        }
    }
    Ok(())
}

fn check_projection(
    source: &LiquidTransferMoments,
    target: &LiquidTransferMoments,
    gravity: [f64; 3],
    limits: TetLiquidTransferLimits,
) -> Result<([f64; 3], f64, f64), TetLiquidTransferError> {
    if !finite_moments(target)
        || (target.mass_kg - source.mass_kg).abs() > limits.mass_tolerance_kg
        || norm(sub(target.momentum_kg_m_s, source.momentum_kg_m_s))
            > limits.momentum_tolerance_kg_m_s
        || norm(sub(
            target.angular_momentum_kg_m2_s,
            source.angular_momentum_kg_m2_s,
        )) > limits.angular_momentum_tolerance_kg_m2_s
    {
        return Err(invalid(
            "target mass, linear or angular momentum exceeds the transfer bound",
        ));
    }
    let shift = sub(
        scale(target.first_mass_moment_kg_m, 1.0 / target.mass_kg),
        scale(source.first_mass_moment_kg_m, 1.0 / source.mass_kg),
    );
    let kinetic = target.kinetic_energy_j - source.kinetic_energy_j;
    let potential = -dot(
        gravity,
        sub(target.first_mass_moment_kg_m, source.first_mass_moment_kg_m),
    );
    if !norm(shift).is_finite()
        || norm(shift) > limits.maximum_centroid_shift_m
        || !kinetic.is_finite()
        || !potential.is_finite()
        || kinetic.abs() + potential.abs() > limits.maximum_projection_energy_change_j
    {
        return Err(invalid(
            "cell projection exceeds the centroid or mechanical-energy error bound; refine",
        ));
    }
    Ok((shift, kinetic, potential))
}

#[cfg(test)]
mod tests;
