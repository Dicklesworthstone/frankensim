//! Spatial vector magnetic history on the original planar P1 field solver.
//!
//! Every triangle owns its original source-selected vector play law and accepted
//! polarization. Newton trials freeze that history and minimize the incremental
//! potential W_new + Q, not storage alone. The existing nonlinear field owner
//! supplies the FEEC tangent, checked CG and bounded globalization. See Egger &
//! Engertsberger (2025), section II, <https://arxiv.org/abs/2507.14521> for the
//! coupled variational formulation; this implementation reuses the existing
//! quadratic branch potential and saturation-ball model, without regularization.
//!
//! This opt-in fixed-geometry, fixed-condition model resolves spatial remanence,
//! minor loops and rotating-field memory. Branchless sourced cells represent a
//! linear background. Every physical reversal must be resolved in the loading
//! history. Numerical endpoint dissipation is reported separately from physical
//! pinning heat. It supplies no eddy current, motion, anisotropy, thermal feedback
//! or measured-material qualification, and no path-independent coenergy claim.

use super::{
    MagneticError, NonlinearMagneticControls, NonlinearMagneticError, NonlinearMagneticReport,
    PlanarMagnetostatics, VariationalEvaluation, finite, gradients, poll, solve_variational,
    weighted_error,
};
use fs_feec::{
    assembly::incidence_to_csr,
    weighted::WeightedAssemblyLimits,
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::field_storage::hysteresis::{
    source::ResolvedVectorHysteresis,
    vector::{
        VectorHysteresisControls, VectorHysteresisError, VectorHysteresisState, VectorHysteresisStep,
    },
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_sparse::Csr;

type Result<T> = core::result::Result<T, NonlinearMagneticError>;

/// Independent field, constitutive and complete work-balance acceptance budgets.
#[derive(Clone, Copy, Debug)]
pub struct PlanarHystereticControls {
    pub nonlinear: NonlinearMagneticControls,
    pub material: VectorHysteresisControls,
    /// Absolute tolerance on the assembled endpoint work identity [J].
    pub energy_tolerance_j: f64,
}

/// Prepared geometry and source laws. Shape gradients and cell measures are
/// cached once; no material card is re-resolved during a field iteration.
pub struct PlanarHystereticMagnetostatics {
    mesh: TriComplex2,
    materials: Vec<ResolvedVectorHysteresis>,
    shapes: Vec<[[f64; 2]; 3]>,
    volumes_m3: Vec<f64>,
    incidence: Csr,
    limits: WeightedAssemblyLimits,
}

/// One accepted field plus ALL of its internal magnetic history. Cloning is a
/// reversible checkpoint; failed or discarded candidates never alter this state.
#[derive(Clone)]
pub struct PlanarHystereticState<'a> {
    problem: &'a PlanarHystereticMagnetostatics,
    potential_wb_m: Vec<f64>,
    magnetic_flux_density_t: Vec<[f64; 2]>,
    cell_states: Vec<VectorHysteresisState>,
    stored_energy_j: f64,
}
impl<'a> PlanarHystereticState<'a> {
    pub fn problem(&self) -> &'a PlanarHystereticMagnetostatics {
        self.problem
    }
    pub fn potential_wb_m(&self) -> &[f64] {
        &self.potential_wb_m
    }
    /// The actual discrete curl of the accepted nodal field, not a substituted
    /// constitutive flux. The inverse residual is independently admitted.
    pub fn magnetic_flux_density_t(&self) -> &[[f64; 2]] {
        &self.magnetic_flux_density_t
    }
    pub fn cell_states(&self) -> &[VectorHysteresisState] {
        &self.cell_states
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
}

/// Local physical and numerical transfers, in original mesh order.
#[derive(Clone, Debug)]
pub struct PlanarHystereticCellStep {
    pub triangle: usize,
    pub volume_m3: f64,
    pub stored_energy_change_j: f64,
    /// Original material pinning heat, available once to a thermal consumer [J].
    pub irreversible_heat_j: f64,
    /// Endpoint-integration loss, never deposited as thermal energy [J].
    pub numerical_dissipation_j: f64,
    /// H_new dot the actual FE curl increment, integrated over this cell [J].
    pub field_work_j: f64,
    pub energy_residual_j: f64,
    /// Difference between the retained material B and actual curl(A) [T].
    pub flux_residual_t: f64,
}

/// Work supplied by impressed currents AND prescribed-boundary reactions.
#[derive(Clone, Debug)]
pub struct PlanarHystereticStep {
    pub cells: Vec<PlanarHystereticCellStep>,
    pub stored_energy_change_j: f64,
    pub irreversible_heat_j: f64,
    pub numerical_dissipation_j: f64,
    pub impressed_current_work_j: f64,
    pub boundary_work_j: f64,
    pub supplied_work_j: f64,
    pub field_work_j: f64,
    /// Supplied endpoint work minus storage change, pinning heat and numerical
    /// dissipation. A small algebraic norm alone cannot accept this identity.
    pub energy_residual_j: f64,
    pub nonlinear: NonlinearMagneticReport,
}

pub struct PlanarHystereticCandidate<'a> {
    state: PlanarHystereticState<'a>,
    step: PlanarHystereticStep,
}
impl<'a> PlanarHystereticCandidate<'a> {
    pub fn state(&self) -> &PlanarHystereticState<'a> {
        &self.state
    }
    pub fn step(&self) -> &PlanarHystereticStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarHystereticState<'a>, PlanarHystereticStep) {
        (self.state, self.step)
    }
}

fn history_error(triangle: usize, error: VectorHysteresisError) -> NonlinearMagneticError {
    match error {
        VectorHysteresisError::Cancelled => MagneticError::Cancelled.into(),
        error => NonlinearMagneticError::History { triangle, error },
    }
}
fn invalid(message: &'static str) -> NonlinearMagneticError {
    MagneticError::Invalid(message).into()
}

struct Evaluation {
    steps: Vec<VectorHysteresisStep>,
    flux: Vec<[f64; 2]>,
    residual: Vec<f64>,
    tangent_weights: Vec<CellWeight2>,
    objective: f64,
}
impl VariationalEvaluation for Evaluation {
    fn residual(&self) -> &[f64] {
        &self.residual
    }
    fn tangent_weights(&self) -> &[CellWeight2] {
        &self.tangent_weights
    }
    fn objective(&self) -> f64 {
        self.objective
    }
}

impl PlanarHystereticMagnetostatics {
    /// Explicit per-triangle sources, including branchless background cells.
    /// Isotropy is selected before construction by ResolvedVectorHysteresis;
    /// neither a material name nor an ordinary permeability implies a loop fit.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedVectorHysteresis>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Planar { .. }) {
            return Err(invalid("planar magnetic history requires the planar A_z metric"));
        }
        if materials.len() != mesh.faces().len()
            || mesh.vertices().len() > limits.max_dofs
            || mesh.edges().len() > limits.max_dofs
            || mesh.faces().len() > limits.max_cells
            || mesh.faces().len().checked_mul(9).is_none_or(|n| n > limits.max_triplets)
        {
            return Err(invalid(
                "one hysteresis source per cell and sufficient assembly limits required",
            ));
        }
        let weights: Vec<_> = materials
            .iter()
            .map(|material| {
                CellWeight2::Scalar(material.law().background_permeability_h_m().recip())
            })
            .collect();
        // Validate the original mesh/Hodge even when later fields have no free
        // coordinates. Actual Newton weights use the full history tangent.
        weighted_mass_matrix_2d(
            &mesh, 1, &weights, Dims([-1, -1, 2, 0, 2, 0]), limits, cancelled,
        )
        .map_err(weighted_error)?;
        let mut shapes = Vec::with_capacity(materials.len());
        let mut volumes_m3 = Vec::with_capacity(materials.len());
        for (triangle, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            shapes.push(gradients(face.map(|v| mesh.vertices()[v as usize]))?);
            volumes_m3.push(
                mesh.face_measure(triangle).ok_or_else(|| invalid("missing cell measure"))?,
            );
        }
        let incidence = incidence_to_csr(&mesh.d0());
        Ok(Self { mesh, materials, shapes, volumes_m3, incidence, limits })
    }

    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    pub fn materials(&self) -> &[ResolvedVectorHysteresis] {
        &self.materials
    }

    /// The virgin zero-field state at an arbitrary finite gauge. Previously
    /// magnetized specimens must be prepared through their actual load history.
    pub fn demagnetized(
        &self,
        gauge_wb_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarHystereticState<'_>> {
        finite(gauge_wb_m)?;
        let mut states = Vec::with_capacity(self.materials.len());
        for material in &self.materials {
            poll(cancelled)?;
            states.push(material.law().initial_state());
        }
        poll(cancelled)?;
        Ok(PlanarHystereticState {
            problem: self,
            potential_wb_m: vec![gauge_wb_m; self.mesh.vertices().len()],
            magnetic_flux_density_t: vec![[0.0; 2]; self.materials.len()],
            cell_states: states,
            stored_energy_j: 0.0,
        })
    }

    fn evaluate(
        &self,
        previous: &PlanarHystereticState<'_>,
        potential: &[f64],
        rhs: &[f64],
        references: &[usize],
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        let mut evaluation = Evaluation {
            steps: Vec::with_capacity(self.materials.len()),
            flux: Vec::with_capacity(self.materials.len()),
            residual: rhs.iter().map(|x| -x).collect(),
            tangent_weights: Vec::with_capacity(self.materials.len()),
            objective: 0.0,
        };
        for (triangle, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let shape = self.shapes[triangle];
            let reference = potential[face[0] as usize];
            let gradient = [0, 1].map(|axis| {
                (potential[face[1] as usize] - reference).mul_add(
                    shape[1][axis],
                    (potential[face[2] as usize] - reference) * shape[2][axis],
                )
            });
            let b = [finite(gradient[1])?, finite(-gradient[0])?];
            let step = self.materials[triangle]
                .law()
                .candidate_flux_density(
                    &previous.cell_states[triangle],
                    [b[0], b[1], 0.0],
                    controls,
                    cancelled,
                )
                .map_err(|e| history_error(triangle, e))?;
            let h = step.state.field_a_m();
            let nu = step.incremental_reluctivity_m_h().map_err(|e| history_error(triangle, e))?;
            let volume = self.volumes_m3[triangle];
            // Relative to accepted storage; this constant changes neither H
            // nor the minimizer and retains small increments around a bias.
            let density = step
                .incremental_potential_change_j_m3()
                .map_err(|e| history_error(triangle, e))?;
            evaluation.objective = finite(volume.mul_add(density, evaluation.objective))?;
            for (local, &vertex) in face.iter().enumerate() {
                let pairing = h[0].mul_add(shape[local][1], -h[1] * shape[local][0]);
                evaluation.residual[vertex as usize] =
                    finite(volume.mul_add(pairing, evaluation.residual[vertex as usize]))?;
            }
            // R^T (dH/dB) R for B=(d_y A,-d_x A). In-plane reduction is
            // valid because isotropic history starts and stays in-plane.
            evaluation.tangent_weights.push(CellWeight2::Tensor([
                [nu[1][1], -nu[1][0]],
                [-nu[0][1], nu[0][0]],
            ]));
            evaluation.flux.push(b);
            evaluation.steps.push(step);
        }
        for (v, &load) in rhs.iter().enumerate() {
            poll(cancelled)?;
            let reference = references[v];
            // Remove the accepted source-work constant as well as storage.
            // Both terms then retain small increments about a biased state.
            let delta = finite(
                (potential[v] - potential[reference])
                    - (previous.potential_wb_m[v] - previous.potential_wb_m[reference]),
            )?;
            evaluation.objective = finite((-load).mul_add(delta, evaluation.objective))?;
        }
        Ok(evaluation)
    }

    /// Solve one physical load increment without changing accepted history.
    /// Source J and prescribed A are endpoint values. Their endpoint work
    /// includes a separately reported numerical remainder; resolve/refine
    /// reversal and rotation histories for continuous-path convergence.
    /// The numerical guess retains the interior field while lifting each
    /// component's gauge, then replaces prescribed values. Use
    /// candidate_from_initial when that guess leaves a source's field domain
    /// despite an admissible equilibrium.
    pub fn candidate<'a>(
        &'a self,
        previous: &PlanarHystereticState<'a>,
        current_density_a_m2: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        controls: PlanarHystereticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarHystereticCandidate<'a>> {
        if !std::ptr::eq(self, previous.problem) {
            return Err(invalid("original magnetic history is required"));
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_boundary(
            &self.mesh, boundary_potential_wb_m, cancelled,
        )?;
        let mut initial = previous.potential_wb_m.clone();
        let references =
            PlanarMagnetostatics::mesh_gauge_references(&self.mesh, &prescribed, cancelled)?;
        for (v, value) in prescribed.iter().enumerate() {
            poll(cancelled)?;
            let reference = references[v];
            initial[v] = if let Some(value) = value {
                *value
            } else {
                finite(
                    (previous.potential_wb_m[v] - previous.potential_wb_m[reference])
                        + prescribed[reference].expect("admitted component gauge"),
                )?
            };
        }
        self.candidate_from_initial(
            previous, current_density_a_m2, boundary_potential_wb_m,
            &initial, controls, cancelled,
        )
    }

    /// Supply a complete admissible NUMERICAL initial field without altering
    /// physical history. Every constitutive trial still starts from previous.
    /// This permits large boundary changes whose equilibrium is in-domain even
    /// when retaining the old interior potential would not be.
    pub fn candidate_from_initial<'a>(
        &'a self,
        previous: &PlanarHystereticState<'a>,
        current_density_a_m2: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        initial_potential_wb_m: &[f64],
        controls: PlanarHystereticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarHystereticCandidate<'a>> {
        poll(cancelled)?;
        if !std::ptr::eq(self, previous.problem)
            || current_density_a_m2.len() != self.materials.len()
            || current_density_a_m2.iter().any(|x| !x.is_finite())
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
        {
            return Err(invalid(
                "original magnetic history, complete finite current and positive work budget required",
            ));
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_boundary(
            &self.mesh, boundary_potential_wb_m, cancelled,
        )?;
        let references =
            PlanarMagnetostatics::mesh_gauge_references(&self.mesh, &prescribed, cancelled)?;
        let mut rhs = vec![0.0; self.mesh.vertices().len()];
        for (triangle, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let load = finite(current_density_a_m2[triangle] * (self.volumes_m3[triangle] / 3.0))?;
            for &vertex in face {
                rhs[vertex as usize] = finite(rhs[vertex as usize] + load)?;
            }
        }
        let (potential, evaluation, nonlinear) = solve_variational(
            &self.mesh, &self.incidence, self.limits, &rhs,
            boundary_potential_wb_m, initial_potential_wb_m, None, controls.nonlinear,
            &mut |potential, _gauge, cancelled| {
                self.evaluate(
                    previous, potential, &rhs, &references, controls.material, cancelled,
                )
            },
            cancelled,
        )?;
        let mut step = PlanarHystereticStep {
            cells: Vec::with_capacity(self.materials.len()),
            stored_energy_change_j: 0.0,
            irreversible_heat_j: 0.0,
            numerical_dissipation_j: 0.0,
            impressed_current_work_j: 0.0,
            boundary_work_j: 0.0,
            supplied_work_j: 0.0,
            field_work_j: 0.0,
            energy_residual_j: 0.0,
            nonlinear,
        };
        let mut cell_states = Vec::with_capacity(self.materials.len());
        let mut energy = 0.0;
        for (triangle, native) in evaluation.steps.into_iter().enumerate() {
            poll(cancelled)?;
            let volume = self.volumes_m3[triangle];
            let h = native.state.field_a_m();
            let before = previous.magnetic_flux_density_t[triangle];
            let after = evaluation.flux[triangle];
            let work = finite(
                volume * h[0].mul_add(after[0] - before[0], h[1] * (after[1] - before[1])),
            )?;
            let change = finite(volume * native.energy_change_j_m3)?;
            let heat = finite(volume * native.irreversible_heat_j_m3)?;
            let numerical = finite(volume * native.numerical_dissipation_j_m3)?;
            let residual = finite(work - change - heat - numerical)?;
            step.stored_energy_change_j = finite(step.stored_energy_change_j + change)?;
            step.irreversible_heat_j = finite(step.irreversible_heat_j + heat)?;
            step.numerical_dissipation_j = finite(step.numerical_dissipation_j + numerical)?;
            step.field_work_j = finite(step.field_work_j + work)?;
            energy = finite(volume.mul_add(native.state.stored_energy_j_m3(), energy))?;
            step.cells.push(PlanarHystereticCellStep {
                triangle,
                volume_m3: volume,
                stored_energy_change_j: change,
                irreversible_heat_j: heat,
                numerical_dissipation_j: numerical,
                field_work_j: work,
                energy_residual_j: residual,
                flux_residual_t: native.inverse_residual_t,
            });
            cell_states.push(native.state);
        }
        for v in 0..potential.len() {
            poll(cancelled)?;
            let reference = references[v];
            // Uniform gauge changes do no physical work. Both current and
            // boundary reaction pairings use the SAME component reference.
            let delta = finite(
                (potential[v] - potential[reference])
                    - (previous.potential_wb_m[v] - previous.potential_wb_m[reference]),
            )?;
            step.impressed_current_work_j =
                finite(rhs[v].mul_add(delta, step.impressed_current_work_j))?;
            if prescribed[v].is_some() {
                step.boundary_work_j =
                    finite(evaluation.residual[v].mul_add(delta, step.boundary_work_j))?;
            }
        }
        step.supplied_work_j = finite(step.impressed_current_work_j + step.boundary_work_j)?;
        step.energy_residual_j = finite(
            step.supplied_work_j - step.stored_energy_change_j
                - step.irreversible_heat_j - step.numerical_dissipation_j,
        )?;
        let local_defect = step.cells.iter().try_fold(0.0, |sum, c| {
            finite(sum + c.energy_residual_j.abs())
        })?;
        if step.energy_residual_j.abs() > controls.energy_tolerance_j
            || local_defect > controls.energy_tolerance_j
            || (step.supplied_work_j - step.field_work_j).abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "spatial hysteresis work, storage, heat and numerical remainder do not close",
            ));
        }
        poll(cancelled)?;
        Ok(PlanarHystereticCandidate {
            state: PlanarHystereticState {
                problem: self,
                potential_wb_m: potential,
                magnetic_flux_density_t: evaluation.flux,
                cell_states,
                stored_energy_j: energy,
            },
            step,
        })
    }
}
