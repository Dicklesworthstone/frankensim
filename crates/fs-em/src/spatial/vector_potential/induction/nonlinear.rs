//! Reversible saturation in the original coupled 3D Whitney A–phi evolution.
//!
//! Backward Euler uses H(B_next) and
//! J=sigma*(-delta_A/dt-grad(phi)+v cross B_next). The original material's
//! exact differential reluctivity enters bounded Newton/FGMRES; the shared
//! induction owner supplies relative tree gauge, floating electric components,
//! source loading and exact P1 current/Ohmic/Lorentz transfers.
//!
//! For convex psi, D=H_next.delta_B-delta_psi is nonnegative numerical
//! dissipation. It is evaluated from original positive source slopes and
//! angular convexity, never by cancellation of large magnetic energies.
//! delta_psi+Ohmic_heat+Lorentz_work+D=impressed_work is independently gated.
//! D is NOT magnetic hysteresis or thermal heat. This first-order fixed-domain
//! law supports reversible isotropic saturation and fixed recoil. Conductor
//! coefficients may be refreshed from the same source at an owned temperature;
//! the original magnetic laws and conducting occupancy remain fixed. No ALE,
//! hysteresis, magnetic-temperature entropy, terminals or coupled motion is
//! inferred. A thermal consumer owns local Joule deposition and temperature.

use super::super::super::{
    SpatialMagneticError, dot, finite, invalid, nonlinear::SpatialNonlinearMaterial, poll,
};
use super::super::{
    apply,
    nonlinear::{
        NonlinearVectorError, SpatialNonlinearVectorField, SpatialNonlinearVectorMagnetostatics,
    },
};
use super::{
    Increment, Jacobi, SpatialInductionCell, SpatialInductionControls, SpatialInductionField,
    SpatialInductionLimits,
    coupling::{Coupling, Geometry, admit_step},
    diagonal,
};
use fs_feec::{
    assembly::stiffness,
    weighted::{CellWeight, WeightedError, weighted_mass_matrix},
};
use fs_material::{
    conductor::ResolvedConductorState,
    field_storage::{admit_source_temperature_refresh, linear_magnetic_permeability},
};
use fs_qty::Dims;
use fs_solver::{FgmresState, LinearOp, ResidualClaim, SolveReport, SolverRunProgress};
use fs_sparse::{Csr, ops};
use std::{
    cell::{Cell, RefCell},
    sync::Arc,
};

#[derive(Clone, Copy, Debug)]
pub struct SpatialNonlinearInductionControls {
    pub induction: SpatialInductionControls,
    pub maximum_newton_iterations: usize,
    pub maximum_backtracks: usize,
}
#[derive(Clone, Debug)]
pub struct SpatialNonlinearInductionReport {
    pub newton_iterations: usize,
    pub linear_iterations: usize,
    pub backtracks: usize,
    pub residual_norm_a: f64,
    pub residual_target_a: f64,
}
#[derive(Debug)]
pub enum SpatialNonlinearInductionError {
    Spatial(SpatialMagneticError),
    Magnetic(NonlinearVectorError),
    Nonconvergence(SpatialNonlinearInductionReport),
}
impl From<SpatialMagneticError> for SpatialNonlinearInductionError {
    fn from(error: SpatialMagneticError) -> Self {
        Self::Spatial(error)
    }
}
impl From<NonlinearVectorError> for SpatialNonlinearInductionError {
    fn from(error: NonlinearVectorError) -> Self {
        Self::Magnetic(error)
    }
}
impl core::fmt::Display for SpatialNonlinearInductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear spatial induction: {self:?}")
    }
}
impl std::error::Error for SpatialNonlinearInductionError {}
type Result<T> = core::result::Result<T, SpatialNonlinearInductionError>;

#[derive(Clone, Debug)]
pub struct SpatialNonlinearInductionMagneticCell {
    pub cell: usize,
    pub volume_m3: f64,
    pub flux_density_t: [[f64; 3]; 2],
    pub free_energy_j: [f64; 2],
    /// Stable signed change of original magnetic free energy.
    pub energy_change_j: f64,
    /// Backward-Euler convexity loss; never deposited as constitutive heat.
    pub numerical_dissipation_j: f64,
}
pub struct SpatialNonlinearInductionStep {
    pub cells: Vec<SpatialInductionCell>,
    pub magnetic_cells: Vec<SpatialNonlinearInductionMagneticCell>,
    pub endpoint_flux_density_t: Vec<[f64; 3]>,
    pub endpoint_field_strength_a_m: Vec<[f64; 3]>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    pub heat_j: f64,
    pub mechanical_work_j: f64,
    pub numerical_dissipation_j: f64,
    pub impressed_work_j: f64,
    pub energy_residual_j: f64,
    pub maximum_ampere_residual_a: f64,
    pub maximum_charge_residual_a: f64,
    pub maximum_cell_magnetic_peclet: f64,
    pub report: SpatialNonlinearInductionReport,
}
#[derive(Clone)]
pub struct SpatialNonlinearInduction<'m> {
    problem: &'m SpatialNonlinearVectorMagnetostatics,
    prepared: Arc<Prepared>,
    potential_wb: Vec<f64>,
    impressed_current_a_m2: Vec<[f64; 3]>,
    time_s: f64,
}
struct Prepared {
    coupling: Coupling,
    permeability_bound: Vec<f64>,
}
pub struct SpatialNonlinearInductionCandidate<'m> {
    state: SpatialNonlinearInduction<'m>,
    step: SpatialNonlinearInductionStep,
}
impl<'m> SpatialNonlinearInductionCandidate<'m> {
    pub fn state(&self) -> &SpatialNonlinearInduction<'m> {
        &self.state
    }
    pub fn step(&self) -> &SpatialNonlinearInductionStep {
        &self.step
    }
    pub fn into_parts(self) -> (SpatialNonlinearInduction<'m>, SpatialNonlinearInductionStep) {
        (self.state, self.step)
    }
}
fn geometry(problem: &SpatialNonlinearVectorMagnetostatics) -> Geometry<'_> {
    Geometry {
        geometry: &problem.geometry,
        cells: &problem.cells,
        curl: &problem.curl,
    }
}
struct Magnetic {
    field: SpatialInductionField,
    energy_j: Vec<f64>,
    effort_a: Vec<f64>,
    tangent: Vec<CellWeight>,
}
fn magnetic(
    problem: &SpatialNonlinearVectorMagnetostatics,
    potential: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Magnetic> {
    let g = geometry(problem);
    let b = g.flux(potential, cancelled)?;
    let mut h = Vec::with_capacity(b.len());
    let mut energy_j = Vec::with_capacity(b.len());
    let mut effort_a = vec![0.0; problem.mesh().edges.len()];
    let mut tangent = Vec::with_capacity(b.len());
    let mut stored_energy_j = 0.0;
    for (cell, &b) in b.iter().enumerate() {
        poll(cancelled)?;
        let response = problem.response(cell, b)?;
        let volume = g.geometry.geometry.vol_signed[cell].abs();
        let energy = finite(volume * response.energy)?;
        stored_energy_j = finite(stored_energy_j + energy)?;
        energy_j.push(energy);
        h.push(response.h);
        tangent.push(CellWeight::Tensor(response.tangent));
        for local in 0..6 {
            let e = problem.cells[cell].edges[local];
            effort_a[e] =
                finite(effort_a[e] + volume * dot(problem.cells[cell].curls[local], response.h))?;
        }
    }
    Ok(Magnetic {
        field: SpatialInductionField {
            flux_density_t: b,
            field_strength_a_m: h,
            stored_energy_j,
        },
        energy_j,
        effort_a,
        tangent,
    })
}
struct Evaluation {
    magnetic: Magnetic,
    potential: Vec<f64>,
    residual: Vec<f64>,
    norm: f64,
}
#[allow(clippy::too_many_arguments)]
fn evaluate(
    problem: &SpatialNonlinearVectorMagnetostatics,
    c: &Coupling,
    old: &[f64],
    x: &mut [f64],
    dt: f64,
    motion: &Csr,
    drive: &[f64],
    load: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Evaluation> {
    let mut potential = old.to_vec();
    for (e, index) in c.edge_unknown.iter().enumerate() {
        poll(cancelled)?;
        if let Some(i) = *index {
            potential[e] = finite(old[e] + x[i])?;
            // All nonlinear and final physical gates see the retained cochain.
            x[i] = finite(potential[e] - old[e])?;
        }
    }
    let magnetic = magnetic(problem, &potential, cancelled)?;
    let mass = apply(&c.mass, x, cancelled)?;
    let movement = apply(motion, x, cancelled)?;
    let mut residual: Vec<_> = mass
        .iter()
        .zip(&movement)
        .zip(drive)
        .map(|((m, v), d)| finite(m / dt - v - d))
        .collect::<core::result::Result<_, _>>()?;
    for (e, index) in c.edge_unknown.iter().enumerate() {
        if let Some(i) = *index {
            residual[i] = finite(residual[i] + magnetic.effort_a[e] - load[e])?;
        }
    }
    let norm = finite(residual.iter().fold(0.0_f64, |n, v| n.hypot(*v)))?;
    Ok(Evaluation {
        magnetic,
        potential,
        residual,
        norm,
    })
}
impl<'m> SpatialNonlinearInduction<'m> {
    /// Continue an admitted saturated/current-carrying field. Its actual fixed
    /// tangential trace and original constitutive sources remain authoritative.
    pub fn new(
        initial: &SpatialNonlinearVectorField<'m>,
        sources: Vec<Option<ResolvedConductorState>>,
        limits: SpatialInductionLimits,
        flux_tolerance_wb: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let problem = initial.problem();
        let entry_bound = problem
            .mesh()
            .tets
            .len()
            .checked_mul(36)
            .ok_or_else(|| invalid("nonlinear induction assembly bound overflows"))?;
        let (coupling, potential_wb) = Coupling::prepare(
            geometry(problem),
            initial.potential_wb(),
            initial.prescribed_edges().to_vec(),
            sources,
            entry_bound,
            limits,
            flux_tolerance_wb,
            |cell, source| {
                problem.materials()[cell]
                    .admit_fixed_material_refresh(source.material())
                    .map_err(SpatialMagneticError::Material)
            },
            cancelled,
        )?;
        let mut permeability_bound = Vec::with_capacity(problem.materials().len());
        for material in problem.materials() {
            poll(cancelled)?;
            let mu = match material {
                SpatialNonlinearMaterial::Linear(m) => {
                    linear_magnetic_permeability(m).map_err(SpatialMagneticError::Material)?
                }
                SpatialNonlinearMaterial::Recoil(m) => m.material().permeability_h_m(),
                SpatialNonlinearMaterial::Reversible(law) => {
                    // Both radial and tangential reluctivity exceed the least
                    // positive original segment slope on the admitted range.
                    let mut minimum = f64::INFINITY;
                    for segment in law.source_path().knots.windows(2) {
                        poll(cancelled)?;
                        let slope = (segment[1].1.evidence.value.value
                            - segment[0].1.evidence.value.value)
                            / (segment[1].0 - segment[0].0);
                        minimum = minimum.min(slope);
                    }
                    minimum.recip()
                }
            };
            if !super::positive(mu) {
                return Err(invalid("finite source permeability bound required").into());
            }
            permeability_bound.push(mu);
        }
        magnetic(problem, &potential_wb, cancelled)?;
        poll(cancelled)?;
        Ok(Self {
            problem,
            prepared: Arc::new(Prepared {
                coupling,
                permeability_bound,
            }),
            potential_wb,
            impressed_current_a_m2: initial.impressed_current_density_a_m2().to_vec(),
            time_s: 0.0,
        })
    }
    pub fn problem(&self) -> &'m SpatialNonlinearVectorMagnetostatics {
        self.problem
    }
    pub fn potential_wb(&self) -> &[f64] {
        &self.potential_wb
    }
    pub fn conductors(&self) -> &[Option<ResolvedConductorState>] {
        &self.prepared.coupling.sources
    }
    pub fn impressed_current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.impressed_current_a_m2
    }
    pub fn tree_edges(&self) -> &[usize] {
        &self.prepared.coupling.tree_edges
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        Ok(magnetic(self.problem, &self.potential_wb, cancelled)?.field)
    }

    /// First-order implicit A–phi step at held nodal velocity and endpoint
    /// impressed current. Source/event reversals require separate intervals.
    /// Every Newton trial uses the same accepted previous A and material laws.
    /// No source, field, clock or conductor state changes before all gates pass.
    pub fn step_candidate(
        &self,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        settings: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearInductionCandidate<'m>> {
        self.step_prepared_candidate(
            Arc::clone(&self.prepared),
            endpoint_current_a_m2,
            velocity_m_s,
            duration_s,
            settings,
            cancelled,
        )
    }

    /// Continue the SAME accepted field with source-resolved endpoint
    /// conductivity, retaining every original conductor claim and quantity.
    /// Only its owned temperature may change. The complete original magnetic
    /// law must remain unchanged, including its full saturation/recoil range;
    /// this is not a zero-work replacement of magnetic storage.
    ///
    /// The existing Whitney/charge operator is reassembled on its original
    /// geometry, relative tree gauge and conducting support. No field is
    /// re-equilibrated or restarted. A finite thermal owner must use the
    /// returned local Joule heat; numerical dissipation remains separate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_conductors_candidate(
        &self,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<ResolvedConductorState>>,
        settings: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearInductionCandidate<'m>> {
        let coupling = self.prepared.coupling.refresh_sources(
            geometry(self.problem),
            conductors,
            |cell, source| {
                let original = self.conductors()[cell]
                    .as_ref()
                    .expect("unchanged conducting support");
                admit_source_temperature_refresh(original.material(), source.material())
                    .map_err(SpatialMagneticError::Material)?;
                self.problem.materials()[cell]
                    .admit_fixed_material_refresh(source.material())
                    .map_err(SpatialMagneticError::Material)?;
                if original.source_basis() != source.source_basis()
                    || original.source_key() != source.source_key()
                    || original
                        .with_material(source.material().clone())
                        .map_err(|_| {
                            invalid("thermal induction left its original electrical source")
                        })?
                        != *source
                {
                    return Err(invalid(
                        "thermal induction changed its original electrical basis or source",
                    ));
                }
                Ok(())
            },
            cancelled,
        )?;
        self.step_prepared_candidate(
            Arc::new(Prepared {
                coupling,
                permeability_bound: self.prepared.permeability_bound.clone(),
            }),
            endpoint_current_a_m2,
            velocity_m_s,
            duration_s,
            settings,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared_candidate(
        &self,
        prepared: Arc<Prepared>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        settings: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearInductionCandidate<'m>> {
        poll(cancelled)?;
        let controls = settings.induction;
        let c = &prepared.coupling;
        let g = geometry(self.problem);
        let (time_s, restart) = admit_step(
            g,
            c,
            self.time_s,
            endpoint_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
        )?;
        if settings
            .maximum_backtracks
            .checked_add(1)
            .and_then(|trials| settings.maximum_newton_iterations.checked_mul(trials))
            .is_none()
            || settings
                .maximum_newton_iterations
                .checked_mul(controls.maximum_cycles)
                .and_then(|n| n.checked_mul(restart))
                .is_none()
        {
            return Err(invalid("nonlinear induction work bound overflows").into());
        }
        let load = c.current_load(
            g,
            endpoint_current_a_m2,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let before = magnetic(self.problem, &self.potential_wb, cancelled)?;
        let (motion, drive, maximum_cell_magnetic_peclet) = c.motion(
            g,
            velocity_m_s,
            &before.field.flux_density_t,
            &prepared.permeability_bound,
            controls.maximum_cell_magnetic_peclet,
            cancelled,
        )?;
        let transpose = ops::transpose(&motion);
        let mut x = vec![0.0; c.mass.nrows()];
        let mut value = evaluate(
            self.problem,
            c,
            &self.potential_wb,
            &mut x,
            duration_s,
            &motion,
            &drive,
            &load,
            cancelled,
        )?;
        let initial_norm = value.norm;
        let target = controls
            .absolute_residual_a
            .max(finite(controls.relative_tolerance * initial_norm)?);
        let mut report = SpatialNonlinearInductionReport {
            newton_iterations: 0,
            linear_iterations: 0,
            backtracks: 0,
            residual_norm_a: value.norm,
            residual_target_a: target,
        };
        while value.norm > target {
            poll(cancelled)?;
            if report.newton_iterations >= settings.maximum_newton_iterations {
                return Err(SpatialNonlinearInductionError::Nonconvergence(report));
            }
            // Original FEEC tensor Hodge and curl-curl owner, with the exact
            // original nonlinear source tangent; no fitted secant stiffness.
            let star = weighted_mass_matrix(
                self.problem.mesh(),
                &self.problem.geometry.geometry,
                2,
                &value.magnetic.tangent,
                Dims([-1, -1, 2, 0, 2, 0]),
                self.problem.limits,
                cancelled,
            )
            .map_err(|e| {
                if e == WeightedError::Cancelled {
                    SpatialMagneticError::Cancelled
                } else {
                    SpatialMagneticError::Assembly(e)
                }
            })?;
            let magnetic = c.restrict(&stiffness(&self.problem.curl, star.matrix()), cancelled)?;
            let inverse = (0..x.len())
                .map(|i| {
                    let d = finite(diagonal(&magnetic, i) + diagonal(&c.mass, i) / duration_s)?;
                    if !super::positive(d) || !super::positive(d.recip()) {
                        Err(invalid(
                            "positive representable nonlinear induction Jacobi scale required",
                        ))
                    } else {
                        Ok(d.recip())
                    }
                })
                .collect::<super::Result<Vec<_>>>()?;
            let rhs: Vec<_> = value.residual.iter().map(|r| -r).collect();
            // Bounded inexact Newton: the linear allowance shrinks with the
            // physical nonlinear residual, and is always checked independently.
            let forcing = 0.1_f64.min((value.norm / initial_norm).sqrt());
            let linear_target = (0.1 * target).max(forcing * value.norm);
            let (direction, linear) = linear_solve(
                &magnetic,
                &c.mass,
                &motion,
                &transpose,
                duration_s,
                &rhs,
                inverse,
                linear_target,
                restart,
                controls.maximum_cycles,
                cancelled,
            )?;
            report.linear_iterations = report
                .linear_iterations
                .checked_add(linear.iters)
                .ok_or_else(|| invalid("nonlinear induction linear work count overflows"))?;
            let mut alpha = 1.0_f64;
            let mut accepted = None;
            let mut source_refusal = None;
            for backtrack in 0..=settings.maximum_backtracks {
                poll(cancelled)?;
                let mut trial = x
                    .iter()
                    .zip(&direction)
                    .map(|(a, d)| finite(alpha.mul_add(*d, *a)))
                    .collect::<super::Result<Vec<_>>>()?;
                match evaluate(
                    self.problem,
                    c,
                    &self.potential_wb,
                    &mut trial,
                    duration_s,
                    &motion,
                    &drive,
                    &load,
                    cancelled,
                ) {
                    Ok(next)
                        if next.norm <= target
                            || next.norm <= (1.0 - 1e-4 * alpha) * value.norm =>
                    {
                        accepted = Some((trial, next, backtrack));
                        break;
                    }
                    Ok(_) => {}
                    Err(
                        error @ SpatialNonlinearInductionError::Magnetic(
                            NonlinearVectorError::FluxRange { .. },
                        ),
                    )
                    | Err(
                        error @ SpatialNonlinearInductionError::Magnetic(
                            NonlinearVectorError::Spatial(SpatialMagneticError::Recoil { .. }),
                        ),
                    ) => {
                        source_refusal = Some(error);
                    }
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            let Some((next_x, next, backtracks)) = accepted else {
                return Err(source_refusal
                    .unwrap_or(SpatialNonlinearInductionError::Nonconvergence(report)));
            };
            x = next_x;
            value = next;
            report.newton_iterations += 1;
            report.backtracks = report
                .backtracks
                .checked_add(backtracks)
                .ok_or_else(|| invalid("nonlinear induction backtrack count overflows"))?;
            report.residual_norm_a = value.norm;
        }
        let delta = value
            .potential
            .iter()
            .zip(&self.potential_wb)
            .map(|(a, b)| finite(a - b))
            .collect::<super::Result<Vec<_>>>()?;
        let current = c.recover(
            g,
            &delta,
            &x,
            &value.magnetic.field.flux_density_t,
            velocity_m_s,
            duration_s,
            cancelled,
        )?;
        let mut maximum_ampere_residual_a = 0.0_f64;
        for e in 0..delta.len() {
            poll(cancelled)?;
            if !c.fixed_edges[e] {
                let residual =
                    finite(value.magnetic.effort_a[e] - load[e] - current.edge_current_a[e])?;
                maximum_ampere_residual_a = maximum_ampere_residual_a.max(residual.abs());
            }
        }
        if maximum_ampere_residual_a > controls.current_tolerance_a
            || current.maximum_charge_residual_a > controls.current_tolerance_a
        {
            return Err(invalid(
                "full nonlinear Ampere or floating-component charge equation failed",
            )
            .into());
        }
        let face_flux = apply(&self.problem.curl, &value.potential, cancelled)?;
        for row in self.problem.mesh().d2().rows {
            poll(cancelled)?;
            let flux = finite(row.iter().map(|&(f, s)| f64::from(s) * face_flux[f]).sum())?;
            if flux.abs() > controls.flux_tolerance_wb {
                return Err(
                    invalid("nonlinear induction violates discrete div-curl allowance").into(),
                );
            }
        }
        let mut magnetic_cells = Vec::with_capacity(g.cells.len());
        let mut stored_change_j = 0.0;
        let mut numerical_dissipation_j = 0.0;
        for cell in 0..g.cells.len() {
            poll(cancelled)?;
            let old = before.field.flux_density_t[cell];
            let new = value.magnetic.field.flux_density_t[cell];
            let db = std::array::from_fn(|j| new[j] - old[j]);
            let volume = g.geometry.geometry.vol_signed[cell].abs();
            let loss = match &self.problem.materials()[cell] {
                SpatialNonlinearMaterial::Reversible(law) => law
                    .bregman_density_3d_j_m3(old, new)
                    .map_err(SpatialMagneticError::Material)?,
                SpatialNonlinearMaterial::Linear(m) => {
                    0.5 * dot(db, db)
                        / linear_magnetic_permeability(m).map_err(SpatialMagneticError::Material)?
                }
                SpatialNonlinearMaterial::Recoil(m) => {
                    0.5 * dot(db, db) / m.material().permeability_h_m()
                }
            };
            let loss = finite(volume * loss)?;
            if loss < 0.0 {
                return Err(invalid("negative nonlinear magnetic numerical dissipation").into());
            }
            let change =
                finite(volume * dot(value.magnetic.field.field_strength_a_m[cell], db) - loss)?;
            stored_change_j = finite(stored_change_j + change)?;
            numerical_dissipation_j = finite(numerical_dissipation_j + loss)?;
            magnetic_cells.push(SpatialNonlinearInductionMagneticCell {
                cell,
                volume_m3: volume,
                flux_density_t: [old, new],
                free_energy_j: [before.energy_j[cell], value.magnetic.energy_j[cell]],
                energy_change_j: change,
                numerical_dissipation_j: loss,
            });
        }
        let impressed_work_j = finite(load.iter().zip(&delta).map(|(j, a)| j * a).sum())?;
        let energy_residual_j = finite(
            stored_change_j + current.heat_j + current.mechanical_work_j + numerical_dissipation_j
                - impressed_work_j,
        )?;
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "nonlinear magnetic, Ohmic, Lorentz and numerical work do not close",
            )
            .into());
        }
        poll(cancelled)?;
        Ok(SpatialNonlinearInductionCandidate {
            state: Self {
                problem: self.problem,
                prepared,
                potential_wb: value.potential,
                impressed_current_a_m2: endpoint_current_a_m2.to_vec(),
                time_s,
            },
            step: SpatialNonlinearInductionStep {
                cells: current.cells,
                magnetic_cells,
                endpoint_flux_density_t: value.magnetic.field.flux_density_t,
                endpoint_field_strength_a_m: value.magnetic.field.field_strength_a_m,
                stored_before_j: before.field.stored_energy_j,
                stored_after_j: value.magnetic.field.stored_energy_j,
                stored_change_j,
                heat_j: current.heat_j,
                mechanical_work_j: current.mechanical_work_j,
                numerical_dissipation_j,
                impressed_work_j,
                energy_residual_j,
                maximum_ampere_residual_a,
                maximum_charge_residual_a: current.maximum_charge_residual_a,
                maximum_cell_magnetic_peclet,
                report,
            },
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn linear_solve(
    magnetic: &Csr,
    mass: &Csr,
    motion: &Csr,
    transpose: &Csr,
    dt: f64,
    rhs: &[f64],
    inverse: Vec<f64>,
    target: f64,
    restart: usize,
    maximum_cycles: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, SolveReport)> {
    let cancel = RefCell::new(cancelled);
    let stopped = Cell::new(false);
    let check = || {
        if stopped.get() {
            true
        } else {
            let yes = (*cancel.borrow_mut())();
            stopped.set(yes);
            yes
        }
    };
    let operator = Increment {
        magnetic,
        mass,
        motion,
        transpose,
        dt,
        magnetic_scale: 1.0,
        motion_scale: 1.0,
        check: &check,
    };
    let norm = finite(rhs.iter().fold(0.0_f64, |n, v| n.hypot(*v)))?;
    let tolerance = finite((0.25 * target / norm).min(0.5))?;
    if !super::positive(tolerance) {
        return Err(invalid("unrepresentable Newton linear residual target").into());
    }
    let mut solver = FgmresState::new(rhs, restart);
    let result = solver.run_controlled(
        &operator,
        &Jacobi(inverse),
        rhs,
        tolerance,
        maximum_cycles,
        &mut || check(),
    );
    if result.progress != SolverRunProgress::Complete || check() {
        return Err(SpatialMagneticError::Cancelled.into());
    }
    let mut ax = vec![0.0; rhs.len()];
    operator.apply(&solver.x, &mut ax);
    if check() {
        return Err(SpatialMagneticError::Cancelled.into());
    }
    let residual = finite(rhs.iter().zip(ax).fold(0.0_f64, |n, (b, a)| n.hypot(b - a)))?;
    let mut history = result.report.history;
    history.push(residual / norm);
    let report = SolveReport::from_claim(
        result.report.iters,
        ResidualClaim::TrueEuclidean(residual / norm),
        target / norm,
        history,
    );
    if !report.converged_euclidean() {
        return Err(SpatialMagneticError::Nonconvergence(report).into());
    }
    Ok((solver.x, report))
}
