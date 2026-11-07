//! Stationary three-dimensional induction with retained vector magnetic memory.
//!
//! The original Whitney A–phi coupling owns conducting components, the relative
//! tree gauge, impressed-current loading and exact local Joule integration.
//! Its shared bounded Newton/FGMRES kernel uses the native vector return map's
//! frozen-history tangent. Every nonlinear trial starts from the same accepted
//! cell polarizations; no trial history is published or substituted by a scalar
//! loop. Source-selected intrinsic viscosities use the same physical duration
//! in every trial. Backward Euler closes
//! delta W + Q_pinning + Q_viscous + Q_Joule + D = W_impressed + W_terminal.
//! Only physical heat terms can enter a thermal receiver. D remains endpoint
//! numerical dissipation, including during source-free relaxation.
//!
//! Geometry, occupancy, magnetic coefficients and their physical frame are
//! fixed. Conductivity can refresh from its original admitted material source.
//! This owner supplies no motion/transport law, anisotropic hysteresis fit,
//! empirical excess-loss correction or experimental material qualification.
//! Intrinsic viscosity fits must exclude the eddy currents resolved here.

/// Physical voltage and total-current electrodes with retained vector memory.
pub mod terminals;
/// Native graph terminals, stranded copper and retained vector magnetic memory.
pub mod winding;

use super::super::super::{
    SpatialGeometry, SpatialMagneticError, dot, finite, invalid, poll, prepare_geometry,
};
use super::super::{CellVectors, VectorTopology, prepare_cells};
use super::{
    SpatialInductionCell, SpatialInductionField, SpatialInductionLimits,
    coupling::{Coupling, Geometry, admit_step},
    kernel::{Kernel, Magnetic, Physical},
    nonlinear::{SpatialNonlinearInductionControls, SpatialNonlinearInductionReport},
    terminals::{SpatialInductionTerminalTransfer, TerminalDrive},
};
use fs_feec::{
    assembly::incidence_to_csr,
    weighted::{CellWeight, WeightedAssemblyLimits},
};
use fs_material::{
    conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial},
    field_storage::{
        admit_source_temperature_refresh,
        hysteresis::{
            source::ResolvedVectorHysteresis,
            vector::{
                VectorHysteresisControls, VectorHysteresisError, VectorHysteresisState,
                VectorHysteresisStep,
            },
        },
    },
};
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
pub struct SpatialHystereticInductionControls {
    pub induction: SpatialNonlinearInductionControls,
    pub material: VectorHysteresisControls,
}

#[derive(Debug)]
pub enum SpatialHystereticInductionError {
    Spatial(SpatialMagneticError),
    Material(VectorHysteresisError),
    Nonconvergence(SpatialNonlinearInductionReport),
}
impl From<SpatialMagneticError> for SpatialHystereticInductionError {
    fn from(error: SpatialMagneticError) -> Self {
        Self::Spatial(error)
    }
}
impl From<VectorHysteresisError> for SpatialHystereticInductionError {
    fn from(error: VectorHysteresisError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for SpatialHystereticInductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hysteretic spatial induction: {self:?}")
    }
}
impl std::error::Error for SpatialHystereticInductionError {}
type Result<T> = core::result::Result<T, SpatialHystereticInductionError>;

/// Original tetrahedral geometry and one explicit vector-law source per cell.
/// Branchless background cells represent air or nonmagnetic conductors without
/// inventing a hysteresis fit from an ordinary permeability property.
pub struct SpatialHystereticProblem {
    geometry: SpatialGeometry,
    cells: Vec<CellVectors>,
    curl: Csr,
    materials: Vec<ResolvedVectorHysteresis>,
    limits: WeightedAssemblyLimits,
}
impl SpatialHystereticProblem {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedVectorHysteresis>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if materials.len() != tetrahedra.len() {
            return Err(
                invalid("one original vector hysteresis source per tetrahedron required").into(),
            );
        }
        let geometry = prepare_geometry(positions_m, tetrahedra, limits, cancelled)?;
        let cells = prepare_cells(&geometry, cancelled)?;
        for material in &materials {
            poll(cancelled)?;
            finite(material.law().background_permeability_h_m().recip())?;
        }
        let curl = incidence_to_csr(&geometry.mesh.d1());
        poll(cancelled)?;
        Ok(Self {
            geometry,
            cells,
            curl,
            materials,
            limits,
        })
    }
    pub fn materials(&self) -> &[ResolvedVectorHysteresis] {
        &self.materials
    }
    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    pub fn cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.geometry.geometry.vol_signed.get(cell).map(|v| v.abs())
    }
    pub fn exterior_faces(&self) -> Vec<usize> {
        self.geometry
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, face)| (face.cells.len() == 1).then_some(i))
            .collect()
    }
    pub fn face_edges(&self, face: usize) -> Option<[usize; 3]> {
        VectorTopology {
            geometry: &self.geometry,
            cells: &self.cells,
        }
        .face_edges(face)
    }
    fn geometry(&self) -> Geometry<'_> {
        Geometry {
            geometry: &self.geometry,
            cells: &self.cells,
            curl: &self.curl,
        }
    }
    fn admit_conductor<C: SpatialConductorMaterial>(
        &self,
        cell: usize,
        source: &C,
    ) -> super::Result<()> {
        let original = &self.materials[cell];
        let material = source.material();
        if original.material().card_identity() != material.card_identity()
            || original.material().material() != material.material()
        {
            return Err(invalid(
                "magnetic fit and conductor must retain the same original material card",
            ));
        }
        let checked = original.at_query_point(&material.query())?;
        admit_source_temperature_refresh(original.material(), checked.material())?;
        for ((axis, before), (_, after)) in original
            .material()
            .query_point()
            .iter()
            .zip(checked.material().query_point())
        {
            if before.to_bits() != after.to_bits() {
                original
                    .admit_temperature_interval(axis, [before.min(*after), before.max(*after)])?;
            }
        }
        Ok(())
    }
}

struct Stationary {
    motion: Csr,
    drive: Vec<f64>,
    velocity: Vec<[f64; 3]>,
}
struct Prepared<C: SpatialConductorMaterial> {
    coupling: Coupling<C>,
    stationary: Arc<Stationary>,
}

/// One accepted cochain, actual vector memories, original conductors and clock.
#[derive(Clone)]
pub struct SpatialHystereticInduction<'m, C: SpatialConductorMaterial = ResolvedConductorState> {
    problem: &'m SpatialHystereticProblem,
    prepared: Arc<Prepared<C>>,
    potential_wb: Vec<f64>,
    cell_states: Vec<VectorHysteresisState>,
    impressed_current_a_m2: Vec<[f64; 3]>,
    time_s: f64,
}

#[derive(Clone, Debug)]
pub struct SpatialHystereticInductionMagneticCell {
    pub cell: usize,
    pub volume_m3: f64,
    /// Actual curl(A) before and after the accepted step.
    pub flux_density_t: [[f64; 3]; 2],
    pub stored_energy_j: [f64; 2],
    /// Stable increment from the original vector constitutive owner.
    pub energy_change_j: f64,
    /// Complete physical constitutive heat: dry pinning plus intrinsic viscosity.
    pub irreversible_heat_j: f64,
    pub pinning_heat_j: f64,
    pub viscous_heat_j: f64,
    pub numerical_dissipation_j: f64,
    /// Endpoint H paired with the actual curl increment minus its accounts.
    pub energy_residual_j: f64,
}
pub struct SpatialHystereticInductionStep {
    pub cells: Vec<SpatialInductionCell>,
    pub magnetic_cells: Vec<SpatialHystereticInductionMagneticCell>,
    pub endpoint_flux_density_t: Vec<[f64; 3]>,
    pub endpoint_field_strength_a_m: Vec<[f64; 3]>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    pub joule_heat_j: f64,
    /// Complete physical constitutive heat: pinning_heat_j + viscous_heat_j.
    pub hysteresis_heat_j: f64,
    pub pinning_heat_j: f64,
    pub viscous_heat_j: f64,
    /// Physical Joule plus dry/viscous switching heat, payable exactly once.
    pub heat_j: f64,
    /// Endpoint numerical loss, excluded from heat_j.
    pub numerical_dissipation_j: f64,
    pub impressed_work_j: f64,
    /// Endpoint electrode work, distinct from the imposed closed-current load.
    pub terminal_work_j: f64,
    pub terminals: Vec<SpatialInductionTerminalTransfer>,
    pub energy_residual_j: f64,
    pub maximum_ampere_residual_a: f64,
    pub maximum_charge_residual_a: f64,
    pub maximum_constitutive_residual_a_m: f64,
    pub maximum_flux_inverse_residual_t: f64,
    pub report: SpatialNonlinearInductionReport,
}
pub struct SpatialHystereticInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialHystereticInduction<'m, C>,
    step: SpatialHystereticInductionStep,
}
impl<'m, C: SpatialConductorMaterial> SpatialHystereticInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialHystereticInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialHystereticInductionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialHystereticInduction<'m, C>,
        SpatialHystereticInductionStep,
    ) {
        (self.state, self.step)
    }
}

impl<'m> SpatialHystereticInduction<'m> {
    /// Start at zero A and the actual demagnetized vector state. Prepare a
    /// magnetized or remanent initial condition through driven candidate steps.
    pub fn new(
        problem: &'m SpatialHystereticProblem,
        sources: Vec<Option<ResolvedConductorState>>,
        prescribed_edge_indices: &[usize],
        limits: SpatialInductionLimits,
        flux_tolerance_wb: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_conductors(
            problem,
            sources,
            prescribed_edge_indices,
            limits,
            flux_tolerance_wb,
            cancelled,
        )
    }
}
impl<'m, C: SpatialConductorMaterial> SpatialHystereticInduction<'m, C> {
    /// Scalar or source-bearing tensor conductivity on the same native A–phi
    /// owner. Prescribed edge traces remain zero throughout this history.
    pub fn new_with_conductors(
        problem: &'m SpatialHystereticProblem,
        sources: Vec<Option<C>>,
        prescribed_edge_indices: &[usize],
        limits: SpatialInductionLimits,
        flux_tolerance_wb: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let mut fixed = vec![false; problem.mesh().edges.len()];
        for &edge in prescribed_edge_indices {
            poll(cancelled)?;
            if edge >= fixed.len() || fixed[edge] {
                return Err(invalid("unique original prescribed edge indices required").into());
            }
            fixed[edge] = true;
        }
        let entry_bound = problem
            .cells
            .len()
            .checked_mul(36)
            .ok_or_else(|| invalid("magnetic assembly bound overflows"))?;
        let (coupling, potential_wb) = Coupling::prepare(
            problem.geometry(),
            &vec![0.0; fixed.len()],
            fixed,
            sources,
            entry_bound,
            limits,
            flux_tolerance_wb,
            |cell, source| problem.admit_conductor(cell, source),
            cancelled,
        )?;
        let n = coupling.mass.nrows();
        let stationary = Arc::new(Stationary {
            motion: Csr::from_parts(n, n, vec![0; n + 1], vec![], vec![]),
            drive: vec![0.0; n],
            velocity: vec![[0.0; 3]; problem.mesh().vertex_count],
        });
        let mut cell_states = Vec::with_capacity(problem.cells.len());
        for material in &problem.materials {
            poll(cancelled)?;
            cell_states.push(material.law().initial_state());
        }
        poll(cancelled)?;
        Ok(Self {
            problem,
            prepared: Arc::new(Prepared {
                coupling,
                stationary,
            }),
            potential_wb,
            cell_states,
            impressed_current_a_m2: vec![[0.0; 3]; problem.cells.len()],
            time_s: 0.0,
        })
    }
    pub fn problem(&self) -> &'m SpatialHystereticProblem {
        self.problem
    }
    pub fn potential_wb(&self) -> &[f64] {
        &self.potential_wb
    }
    pub fn cell_states(&self) -> &[VectorHysteresisState] {
        &self.cell_states
    }
    pub fn conductors(&self) -> &[Option<C>] {
        &self.prepared.coupling.sources
    }
    pub fn tree_edges(&self) -> &[usize] {
        &self.prepared.coupling.tree_edges
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn impressed_current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.impressed_current_a_m2
    }

    /// Recover the accepted curl field and its retained constitutive H without
    /// advancing, re-equilibrating or resetting any cell memory.
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        let flux_density_t = self
            .problem
            .geometry()
            .flux(&self.potential_wb, cancelled)?;
        let mut field_strength_a_m = Vec::with_capacity(self.cell_states.len());
        let mut stored_energy_j = 0.0;
        for (cell, state) in self.cell_states.iter().enumerate() {
            poll(cancelled)?;
            field_strength_a_m.push(state.field_a_m());
            stored_energy_j = finite(
                stored_energy_j
                    + self.problem.geometry.geometry.vol_signed[cell].abs()
                        * state.stored_energy_j_m3(),
            )?;
        }
        Ok(SpatialInductionField {
            flux_density_t,
            field_strength_a_m,
            stored_energy_j,
        })
    }

    pub fn step_candidate(
        &self,
        endpoint_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticInductionCandidate<'m, C>> {
        self.step_prepared_candidate(
            self.prepared.clone(),
            endpoint_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Refresh only conductivity from the same original electrical sources and
    /// fixed magnetic fit. The original occupancy and accepted memories remain.
    pub fn step_with_conductors_candidate(
        &self,
        endpoint_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticInductionCandidate<'m, C>> {
        let prepared = self.prepare_conductors(conductors, cancelled)?;
        self.step_prepared_candidate(
            prepared,
            endpoint_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    // Terminal Newton trials share this conductivity update while each starts
    // from the SAME accepted vector memories and edge cochains.
    fn prepare_conductors(
        &self,
        conductors: Vec<Option<C>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Arc<Prepared<C>>> {
        poll(cancelled)?;
        if conductors == self.prepared.coupling.sources {
            return Ok(self.prepared.clone());
        }
        let coupling = self.prepared.coupling.refresh_sources(
            self.problem.geometry(),
            conductors,
            |cell, source| {
                let original = self.conductors()[cell]
                    .as_ref()
                    .expect("unchanged conducting support");
                admit_source_temperature_refresh(original.material(), source.material())?;
                self.problem.admit_conductor(cell, source)?;
                if !original.has_same_electrical_source(source)
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
        Ok(Arc::new(Prepared {
            coupling,
            stationary: self.prepared.stationary.clone(),
        }))
    }

    fn magnetic(
        &self,
        potential: &[f64],
        duration_s: f64,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Magnetic<Vec<VectorHysteresisStep>>> {
        let b = self.problem.geometry().flux(potential, cancelled)?;
        let mut h = Vec::with_capacity(b.len());
        let mut energy_j = Vec::with_capacity(b.len());
        let mut effort_a = vec![0.0; self.problem.mesh().edges.len()];
        let mut tangent = Vec::with_capacity(b.len());
        let mut states = Vec::with_capacity(b.len());
        let mut stored_energy_j = 0.0;
        for (cell, &b) in b.iter().enumerate() {
            poll(cancelled)?;
            let step = self.problem.materials[cell]
                .law()
                .candidate_flux_density_timed(
                    &self.cell_states[cell],
                    b,
                    duration_s,
                    controls,
                    cancelled,
                )?;
            let h_cell = step.state.field_a_m();
            let volume = self.problem.geometry.geometry.vol_signed[cell].abs();
            let energy = finite(volume * step.state.stored_energy_j_m3())?;
            stored_energy_j = finite(stored_energy_j + energy)?;
            energy_j.push(energy);
            h.push(h_cell);
            tangent.push(CellWeight::Tensor(step.incremental_reluctivity_m_h()?));
            let vectors = &self.problem.cells[cell];
            for local in 0..6 {
                let edge = vectors.edges[local];
                effort_a[edge] =
                    finite(effort_a[edge] + volume * dot(vectors.curls[local], h_cell))?;
            }
            states.push(step);
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
            state: states,
        })
    }

    fn step_prepared_candidate(
        &self,
        prepared: Arc<Prepared<C>>,
        endpoint_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        settings: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticInductionCandidate<'m, C>> {
        self.step_driven_candidate(
            prepared,
            endpoint_current_a_m2,
            duration_s,
            settings,
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_driven_candidate(
        &self,
        prepared: Arc<Prepared<C>>,
        endpoint_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        settings: SpatialHystereticInductionControls,
        terminal_drive: Option<&TerminalDrive<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticInductionCandidate<'m, C>> {
        poll(cancelled)?;
        let g = self.problem.geometry();
        let c = &prepared.coupling;
        let stationary = &prepared.stationary;
        let controls = settings.induction.induction;
        let (time_s, restart) = admit_step(
            g,
            c,
            self.time_s,
            endpoint_current_a_m2,
            &stationary.velocity,
            duration_s,
            controls,
        )?;
        let before = self.field(cancelled)?;
        let load = c.current_load(
            g,
            endpoint_current_a_m2,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let drive = if let Some(terminal) = terminal_drive {
            let mut drive = stationary.drive.clone();
            terminal.add_rhs(c, &mut drive, cancelled)?;
            std::borrow::Cow::Owned(drive)
        } else {
            std::borrow::Cow::Borrowed(stationary.drive.as_slice())
        };
        let kernel = Kernel {
            geometry: g,
            coupling: c,
            accepted_potential: &self.potential_wb,
            duration_s,
            motion: &stationary.motion,
            drive: &drive,
            load: &load,
            assembly: self.problem.limits,
            controls: settings.induction,
            restart,
        };
        // The original A–phi Newton kernel globalizes on its complete residual
        // norm. Its effort and tangent come from the native incremental law,
        // whose potential includes HALF the physical viscous heat; the energy
        // retained above is the original storage, without any dissipative term.
        let value = kernel.solve(
            |potential, mut stop| {
                self.magnetic(potential, duration_s, settings.material, &mut stop)
            },
            |error| {
                matches!(
                    error,
                    SpatialHystereticInductionError::Material(
                        VectorHysteresisError::FluxNonconvergence { .. }
                            | VectorHysteresisError::FluxRange { .. }
                    ) | SpatialHystereticInductionError::Material(VectorHysteresisError::Invalid(
                        "vector magnetic field exceeds its original magnitude domain"
                    ))
                )
            },
            SpatialHystereticInductionError::Nonconvergence,
            cancelled,
        )?;
        let Physical {
            delta,
            current,
            maximum_ampere_residual_a,
            maximum_charge_residual_a,
            terminals,
            terminal_work_j,
        } = kernel.recover(&value, &stationary.velocity, terminal_drive, cancelled)?;
        let mut magnetic_cells = Vec::with_capacity(g.cells.len());
        let mut stored_change_j = 0.0;
        let mut hysteresis_heat_j = 0.0;
        let mut pinning_heat_j = 0.0;
        let mut viscous_heat_j = 0.0;
        let mut numerical_dissipation_j = 0.0;
        let mut constitutive_work_defect_j = 0.0;
        let mut maximum_constitutive_residual_a_m = 0.0_f64;
        let mut maximum_flux_inverse_residual_t = 0.0_f64;
        for (cell, step) in value.magnetic.state.iter().enumerate() {
            poll(cancelled)?;
            let volume = g.geometry.geometry.vol_signed[cell].abs();
            let old = before.flux_density_t[cell];
            let new = value.magnetic.field.flux_density_t[cell];
            let db = std::array::from_fn(|j| new[j] - old[j]);
            let change = finite(volume * step.energy_change_j_m3)?;
            let heat = finite(volume * step.irreversible_heat_j_m3)?;
            let pinning = finite(volume * step.pinning_heat_j_m3)?;
            let viscous = finite(volume * step.viscous_heat_j_m3)?;
            let loss = finite(volume * step.numerical_dissipation_j_m3)?;
            let defect = finite(volume * dot(step.state.field_a_m(), db) - change - heat - loss)?;
            // Preserve actual local memory and independently admit its inverse
            // error against curl(A); never overwrite the returned material B.
            constitutive_work_defect_j = finite(constitutive_work_defect_j + defect.abs())?;
            stored_change_j = finite(stored_change_j + change)?;
            hysteresis_heat_j = finite(hysteresis_heat_j + heat)?;
            pinning_heat_j = finite(pinning_heat_j + pinning)?;
            viscous_heat_j = finite(viscous_heat_j + viscous)?;
            numerical_dissipation_j = finite(numerical_dissipation_j + loss)?;
            maximum_constitutive_residual_a_m =
                maximum_constitutive_residual_a_m.max(step.maximum_constitutive_residual_a_m);
            maximum_flux_inverse_residual_t =
                maximum_flux_inverse_residual_t.max(step.inverse_residual_t);
            magnetic_cells.push(SpatialHystereticInductionMagneticCell {
                cell,
                volume_m3: volume,
                flux_density_t: [old, new],
                stored_energy_j: [
                    finite(volume * self.cell_states[cell].stored_energy_j_m3())?,
                    value.magnetic.energy_j[cell],
                ],
                energy_change_j: change,
                irreversible_heat_j: heat,
                pinning_heat_j: pinning,
                viscous_heat_j: viscous,
                numerical_dissipation_j: loss,
                energy_residual_j: defect,
            });
        }
        let heat_j = finite(current.heat_j + hysteresis_heat_j)?;
        let impressed_work_j = finite(load.iter().zip(&delta).map(|(j, a)| j * a).sum())?;
        let energy_residual_j = finite(
            stored_change_j + heat_j + numerical_dissipation_j - impressed_work_j - terminal_work_j,
        )?;
        if constitutive_work_defect_j > controls.energy_tolerance_j
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "hysteretic magnetic, Joule, constitutive, terminal and numerical work do not close",
            )
            .into());
        }
        poll(cancelled)?;
        Ok(SpatialHystereticInductionCandidate {
            state: Self {
                problem: self.problem,
                prepared,
                potential_wb: value.potential,
                cell_states: value
                    .magnetic
                    .state
                    .into_iter()
                    .map(|step| step.state)
                    .collect(),
                impressed_current_a_m2: endpoint_current_a_m2.to_vec(),
                time_s,
            },
            step: SpatialHystereticInductionStep {
                cells: current.cells,
                magnetic_cells,
                endpoint_flux_density_t: value.magnetic.field.flux_density_t,
                endpoint_field_strength_a_m: value.magnetic.field.field_strength_a_m,
                stored_before_j: before.stored_energy_j,
                stored_after_j: value.magnetic.field.stored_energy_j,
                stored_change_j,
                joule_heat_j: current.heat_j,
                hysteresis_heat_j,
                pinning_heat_j,
                viscous_heat_j,
                heat_j,
                numerical_dissipation_j,
                impressed_work_j,
                terminal_work_j,
                terminals,
                energy_residual_j,
                maximum_ampere_residual_a,
                maximum_charge_residual_a,
                maximum_constitutive_residual_a_m,
                maximum_flux_inverse_residual_t,
                report: value.report,
            },
        })
    }
}
