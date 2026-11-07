//! Discrete adjoints of the authentic two-stage primary-enthalpy SDIRK2 step.
//!
//! Write a=gamma*dt and r=(1-gamma)*dt/a. The accepted first-stage
//! all-node rate, including prescribed-node reaction, is (E(h1)-E(h0))/a.
//! Eliminating that rate gives the second free-row residual
//! E(h2)-(1-r)*E(h0)-r*E(h1)-a*Q2(h2), with E=m*h+Cfilm*T(h).
//! Its block transpose needs the actual storage slopes at h0, h1 and h2;
//! the interstage carry uses E'(h1), while history uses E'(h0). Neither a
//! fictitious old enthalpy nor a backward-Euler endpoint receipt is needed.
//!
//! Geometry, reference masses, original constitutive data and the time grid
//! stay fixed. The separate stage source and prescribed-temperature controls
//! are independent samples; callers contract them with their forcing model.

use super::super::sdirk::{
    EnthalpySdirk2Input, EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA,
};
use super::*;
use crate::interface::temperature::TemperatureContactEvaluation;

/// Complete two-stage implicit pullback. Array order is first, second stage.
#[derive(Debug)]
pub struct EnthalpySdirk2Gradient {
    pub previous_specific_enthalpy: Vec<f64>,
    /// P1 nodal volumetric-source-density derivatives at each tableau time.
    pub source_density: [Vec<f64>; 2],
    /// Bath-temperature derivatives, zero at each stage's unconstrained nodes.
    pub prescribed_temperature: [Vec<f64>; 2],
    /// Multipliers of the eliminated block residual documented by this module.
    pub stage_adjoint: [Vec<f64>; 2],
    pub relative_residual: [f64; 2],
    /// Krylov columns consumed at each stage, each bounded by max_iterations.
    pub iterations: [usize; 2],
}

/// Elementwise conductivity controls of both genuine SDIRK2 stages.
#[derive(Debug, Clone, PartialEq)]
pub struct EnthalpySdirk2ConductivityGradient {
    /// Independent full-law scale derivatives at the first and second stage.
    /// Each vector uses mesh element order and already includes its complete
    /// eliminated-block adjoint and diagonal tableau time factor.
    pub stage_scales: [Vec<f64>; 2],
    /// Derivative when the same element scale controls both stage laws.
    /// Exactly the elementwise sum of `stage_scales`.
    pub shared_scales: Vec<f64>,
}

#[derive(Debug)]
struct StageTangent {
    tangent: Csr,
    slopes: Vec<f64>,
    storage_diagonal: Vec<f64>,
    prescribed: Vec<bool>,
    inverse_diagonal: Vec<f64>,
    // Original stage K_eff*T only, excluding storage and all external loads.
    conductive_loads: Vec<[f64; 4]>,
    radiation: Option<PointwiseBoundaryAssembly>,
}

/// Both accepted physical stages and their exact spatial h tangents.
/// State/contact/radiation feedback is retained even on latent plateaus.
/// This owner never exposes a stage as a backward-Euler accepted interval.
#[derive(Debug)]
pub struct EnthalpySdirk2Linearization<'m> {
    primal: EnthalpySdirk2Step,
    mesh: &'m ConductionMesh,
    masses: Vec<f64>,
    previous_storage_diagonal: Vec<f64>,
    diagonal_dt: f64,
    carry_ratio: f64,
    stages: [StageTangent; 2],
}

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Solve the actual SDIRK2 interval, then bind both original stage tangents.
    ///
    /// Supports the primal's fixed interfaces, smooth source-dependent contacts
    /// and pointwise/material ambient or enclosure radiation. Strict original
    /// chart/source derivatives must exist at BOTH stages; a forward generalized
    /// derivative at a corner does not authorize a classical adjoint.
    /// `max_radiation_entries` caps the SUM of retained radiation work allowances
    /// and may be zero without radiation. No trajectory/checkpoint allocation
    /// is implicit; callers can compose accepted steps on a fixed time grid.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_sdirk2(
        &self,
        cx: &Cx<'_>,
        inputs: [EnthalpySdirk2Input<'_>; 2],
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        max_radiation_entries: usize,
    ) -> Result<EnthalpySdirk2Linearization<'m>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        let mut entries = [0; 2];
        for (index, input) in inputs.iter().enumerate() {
            entries[index] = radiation_entries(input.problem, self.sdirk_radiation(cx, *input)?)?;
        }
        let total =
            entries[0]
                .checked_add(entries[1])
                .ok_or(EnthalpyAdjointError::InvalidInput(
                    "SDIRK2 radiation workspace overflow",
                ))?;
        if total > max_radiation_entries {
            return Err(EnthalpyAdjointError::InvalidInput(
                "SDIRK2 radiation adjoints exceed their combined entry allowance",
            ));
        }
        total
            .checked_mul(std::mem::size_of::<f64>())
            .filter(|&bytes| isize::try_from(bytes).is_ok())
            .ok_or(EnthalpyAdjointError::InvalidInput(
                "SDIRK2 radiation address space overflow",
            ))?;
        let previous_storage_diagonal = storage_history_derivatives(cx, self, old_h)?;
        let primal = self.advance_sdirk2(cx, inputs, old_h, dt_s, config)?;
        let diagonal_dt = finite(SDIRK2_GAMMA * dt_s)?;
        let carry_ratio = finite(((1.0 - SDIRK2_GAMMA) * dt_s) / diagonal_dt)?;
        let first = self.sdirk_stage_tangent(
            cx,
            inputs[0],
            old_h,
            diagonal_dt,
            &primal.first_stage,
            entries[0],
        )?;
        let second = self.sdirk_stage_tangent(
            cx,
            inputs[1],
            old_h,
            diagonal_dt,
            &primal.second_stage,
            entries[1],
        )?;
        poll(cx, self.masses.len())?;
        Ok(EnthalpySdirk2Linearization {
            primal,
            mesh: self.mesh,
            masses: self.masses.clone(),
            previous_storage_diagonal,
            diagonal_dt,
            carry_ratio,
            stages: [first, second],
        })
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn sdirk_stage_tangent(
        &self,
        cx: &Cx<'_>,
        input: EnthalpySdirk2Input<'_>,
        old_h: &[f64],
        dt: f64,
        primal: &EnthalpySdirk2Stage,
        radiation_entries: usize,
    ) -> Result<StageTangent, EnthalpyAdjointError> {
        let h = &primal.specific_enthalpy_j_kg;
        check_chart(cx, self, h)?;
        if !crate::adjoint::robin::nonlinear::material_is_smooth(
            cx,
            input.problem,
            &primal.temperature,
        )? {
            return Err(EnthalpyAdjointError::MaterialKink);
        }
        let radiation = self.sdirk_radiation(cx, input)?;
        let context = StepContext {
            storage: self,
            cx,
            problem: input.problem,
            interfaces: primal
                .physical_contacts
                .as_ref()
                .map(TemperatureContactEvaluation::frozen_interfaces)
                .or(input.interfaces),
            old: old_h,
            dt,
            pointwise_radiation: radiation,
        };
        // Stage construction depends on the actual h, never an extrapolated
        // history. We retain only its physical tangent, not its BE residual.
        let mut stage = context.stage(h)?;
        match (radiation, &mut stage.radiation) {
            (
                Some(PointwiseBoundary::MaterialAmbient(owner)),
                Some(PointwiseBoundaryAssembly::Ambient(assembly)),
            ) => {
                if primal.pointwise_radiation.as_ref() != Some(&assembly.report) {
                    return Err(EnthalpyAdjointError::InvalidInput(
                        "SDIRK2 material radiation stage receipts differ",
                    ));
                }
                owner.bind_adjoint(cx, h, assembly, radiation_entries)?;
            }
            (
                Some(PointwiseBoundary::MaterialEnclosure(owner)),
                Some(PointwiseBoundaryAssembly::MaterialEnclosure(assembly)),
            ) => {
                if primal.pointwise_enclosure_radiation.as_ref() != Some(&assembly.report) {
                    return Err(EnthalpyAdjointError::InvalidInput(
                        "SDIRK2 material enclosure stage receipts differ",
                    ));
                }
                owner.bind_adjoint(cx, h, assembly, radiation_entries)?;
            }
            _ => {}
        }
        if let Some(radiation) = &stage.radiation {
            radiation.admit_adjoint(cx)?;
        }
        if let Some((owner, config)) = input.contacts {
            owner.check_mesh(input.problem.mesh)?;
            let contact = owner.linearize_enthalpy_at_fields(
                cx,
                input.problem.boundary,
                h,
                &primal.temperature,
                &stage.slope,
                config.max_face_pairs,
                |vertex| Some(self.curve_for_vertex(vertex)),
            )?;
            let original =
                primal
                    .physical_contacts
                    .as_ref()
                    .ok_or(EnthalpyAdjointError::InvalidInput(
                        "SDIRK2 physical contact stage missing",
                    ))?;
            let bound = contact.physical_contacts();
            if bound.states() != original.states()
                || bound.quadrature_states() != original.quadrature_states()
                || bound.material_node_states() != original.material_node_states()
                || bound.material_trace_states() != original.material_trace_states()
            {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "SDIRK2 original contact stage receipts differ",
                ));
            }
            let mut coo = Coo::new(h.len(), h.len());
            for row in 0..h.len() {
                if row % ASSEMBLY_TILE == 0 {
                    poll(cx, row)?;
                }
                let (columns, entries) = stage.tangent.row(row);
                for (&column, &entry) in columns.iter().zip(entries) {
                    coo.push(row, column, entry);
                }
                let (columns, entries) = contact.feedback().row(row);
                for (&column, &entry) in columns.iter().zip(entries) {
                    coo.push(row, column, finite(entry * stage.slope[column])?);
                }
                // Material fraction feedback is already in h coordinates.
                if let Some(phase) = contact.phase_feedback() {
                    let (columns, entries) = phase.row(row);
                    for (&column, &entry) in columns.iter().zip(entries) {
                        coo.push(row, column, entry);
                    }
                }
            }
            stage.tangent = coo.assemble();
        }
        let mut prescribed = vec![false; h.len()];
        for &(vertex, _) in input.problem.boundary.dirichlet() {
            prescribed[vertex] = true;
        }
        let mut inverse_diagonal = Vec::with_capacity(h.len());
        for row in 0..h.len() {
            if row % ASSEMBLY_TILE == 0 {
                poll(cx, row)?;
            }
            let diagonal = if prescribed[row] {
                self.masses[row]
            } else {
                let radiation_diagonal = stage
                    .radiation
                    .as_ref()
                    .map_or(0.0, |radiation| radiation.preconditioner_diagonal(row));
                finite(
                    stage.storage_diagonal[row]
                        + dt * finite(
                            stage.tangent.get(row, row) + radiation_diagonal * stage.slope[row],
                        )?,
                )?
            };
            if diagonal == 0.0 {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "nonzero SDIRK2 stage Jacobian diagonal required",
                ));
            }
            inverse_diagonal.push(finite(1.0 / diagonal)?);
        }
        let conductive_loads = {
            let conductivity = context.phase_scales(h)?;
            super::conductivity::element_loads(cx, input.problem, &primal.temperature, conductivity.as_ref())?
        };
        poll(cx, h.len())?;
        Ok(StageTangent {
            tangent: stage.tangent,
            slopes: stage.slope,
            storage_diagonal: stage.storage_diagonal,
            prescribed,
            inverse_diagonal,
            conductive_loads,
            radiation: stage.radiation,
        })
    }
}

impl EnthalpySdirk2Linearization<'_> {
    /// Genuine interval receipt with separate stage powers and tableau heat.
    #[must_use]
    pub const fn primal(&self) -> &EnthalpySdirk2Step {
        &self.primal
    }

    fn operator(&self, stage: usize) -> EnthalpyTangent<'_, '_> {
        let stage = &self.stages[stage];
        EnthalpyTangent {
            mesh: self.mesh,
            tangent: &stage.tangent,
            slopes: &stage.slopes,
            masses: &self.masses,
            storage_diagonal: &stage.storage_diagonal,
            previous_storage_diagonal: &self.previous_storage_diagonal,
            prescribed: &stage.prescribed,
            inverse_diagonal: &stage.inverse_diagonal,
            dt: self.diagonal_dt,
            feedback: &[],
            pointwise_radiation: stage.radiation.as_ref(),
        }
    }

    /// Convert two independent stage temperature cotangents to h coordinates.
    /// Add direct stage-enthalpy terms before calling `pullback_stages`.
    pub fn temperature_pullback(
        &self,
        cx: &Cx<'_>,
        seeds: [&[f64]; 2],
    ) -> Result<[Vec<f64>; 2], EnthalpyAdjointError> {
        let mut result = [Vec::new(), Vec::new()];
        for stage in 0..2 {
            vector(cx, seeds[stage], self.masses.len())?;
            result[stage].reserve(self.masses.len());
            for (vertex, (&seed, &slope)) in seeds[stage]
                .iter()
                .zip(&self.stages[stage].slopes)
                .enumerate()
            {
                if vertex % ASSEMBLY_TILE == 0 {
                    poll(cx, vertex)?;
                }
                result[stage].push(finite(seed * slope)?);
            }
        }
        poll(cx, self.masses.len())?;
        Ok(result)
    }

    /// Pull the complete two-stage adjoint back to element conductivity scales.
    ///
    /// At each stage the control is `K_eff,e(T,f) -> s_e*K_eff,e(T,f)` at
    /// `s_e=1`, equivalently a log scale at zero. The retained load includes
    /// the original directional temperature law and any admitted phase
    /// multiplier or source-resolved tensor. Neither storage, source, boundary
    /// nor contact coefficients are scaled.
    ///
    /// Pass the unchanged `stage_adjoint` from [Self::pullback] or
    /// [Self::pullback_stages]. Each contribution is
    /// `-gamma*dt*lambda_stage^T*(K_eff*T)_stage`, with that stage's
    /// prescribed rows removed. The eliminated block transpose ALREADY carries
    /// first-stage influence: adding a separate off-diagonal tableau term or
    /// another k'(T)/phase chain rule would count it twice.
    ///
    /// Stage contributions describe independent controls at the two tableau
    /// times. `shared_scales` sums them for a material control held throughout
    /// the interval; a trajectory additionally propagates the existing history
    /// pullback and sums accepted interval contributions. Explicit objective
    /// dependence stays caller-owned.
    ///
    /// No state, material query, sparse assembly or transpose solve is repeated.
    /// The linearization retains 64 bytes per element for these two loads;
    /// this call allocates 24 bytes per element for the two contributions and
    /// shared sum. `max_elements` bounds output/traversal, not total RSS.
    ///
    /// # Errors
    /// Wrong or nonfinite full-stage multipliers, exhausted element/address
    /// allowance, allocation failure, overflow, and cancellation return no
    /// partial gradient or changed accepted state.
    pub fn conductivity_scale_pullback(
        &self,
        cx: &Cx<'_>,
        full_stage_joule_multipliers: [&[f64]; 2],
        max_elements: usize,
    ) -> Result<EnthalpySdirk2ConductivityGradient, EnthalpyAdjointError> {
        poll(cx, 0)?;
        let elements = self.mesh.element_count();
        if elements > max_elements
            || self.stages.iter().any(|stage| stage.conductive_loads.len() != elements)
            || elements
                .checked_mul(std::mem::size_of::<[f64; 3]>())
                .is_none_or(|bytes| isize::try_from(bytes).is_err())
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "SDIRK2 conductivity pullback exceeds its element allowance",
            ));
        }
        for multiplier in full_stage_joule_multipliers {
            vector(cx, multiplier, self.masses.len())?;
        }
        let mut stage_scales = [Vec::new(), Vec::new()];
        let mut shared_scales = Vec::new();
        for output in stage_scales.iter_mut().chain(std::iter::once(&mut shared_scales)) {
            output.try_reserve_exact(elements).map_err(|_| {
                EnthalpyAdjointError::InvalidInput("SDIRK2 conductivity gradient allocation refused")
            })?;
        }
        for (element, &tet) in self.mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let mut contributions = [0.0; 2];
            for stage in 0..2 {
                contributions[stage] = stage_conductivity_contraction(
                    self.diagonal_dt,
                    self.stages[stage].conductive_loads[element],
                    tet.map(|v| full_stage_joule_multipliers[stage][v as usize]),
                    tet.map(|v| self.stages[stage].prescribed[v as usize]),
                )
                .ok_or(EnthalpyAdjointError::NonFiniteArithmetic)?;
                stage_scales[stage].push(contributions[stage]);
            }
            shared_scales.push(finite(contributions[0] + contributions[1])?);
        }
        poll(cx, elements)?;
        Ok(EnthalpySdirk2ConductivityGradient { stage_scales, shared_scales })
    }

    /// Apply the complete block Jacobian of the eliminated two-stage balance.
    /// Directions and outputs are ordered [first stage, second stage].
    pub fn apply_jacobian(
        &self,
        cx: &Cx<'_>,
        directions: [&[f64]; 2],
    ) -> Result<[Vec<f64>; 2], EnthalpyAdjointError> {
        self.action(cx, directions, false)
    }

    /// Exact transpose of `apply_jacobian`, including the interstage mass carry.
    pub fn apply_jacobian_transpose(
        &self,
        cx: &Cx<'_>,
        directions: [&[f64]; 2],
    ) -> Result<[Vec<f64>; 2], EnthalpyAdjointError> {
        self.action(cx, directions, true)
    }

    fn action(
        &self,
        cx: &Cx<'_>,
        directions: [&[f64]; 2],
        transpose: bool,
    ) -> Result<[Vec<f64>; 2], EnthalpyAdjointError> {
        let mut output = [
            self.operator(0).action(cx, directions[0], transpose)?,
            self.operator(1).action(cx, directions[1], transpose)?,
        ];
        for (vertex, &capacity) in self.stages[0].storage_diagonal.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if !self.stages[1].prescribed[vertex] {
                let (target, source) = if transpose { (0, 1) } else { (1, 0) };
                output[target][vertex] = finite(
                    output[target][vertex]
                        - finite(self.carry_ratio * finite(capacity * directions[source][vertex])?)?,
                )?;
            }
        }
        poll(cx, self.masses.len())?;
        Ok(output)
    }

    /// Endpoint-only objective/carry pullback; the second stage is the endpoint.
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        config: LinearConfig,
    ) -> Result<EnthalpySdirk2Gradient, EnthalpyAdjointError> {
        self.pullback_stages(cx, [&vec![0.0; self.masses.len()], seed], config)
    }

    /// Pull back arbitrary h objectives at both original physical stages.
    /// Each transpose uses the shared checked FGMRES implementation and its
    /// own `max_iterations` column cap; no gradient is returned if either fails.
    /// The history chain includes stage-one prescribed nodes that become free
    /// in stage two. Direct objective dependence on input controls is caller-owned.
    pub fn pullback_stages(
        &self,
        cx: &Cx<'_>,
        seeds: [&[f64]; 2],
        config: LinearConfig,
    ) -> Result<EnthalpySdirk2Gradient, EnthalpyAdjointError> {
        vector(cx, seeds[0], self.masses.len())?;
        vector(cx, seeds[1], self.masses.len())?;
        config
            .max_iterations
            .checked_mul(2)
            .ok_or(EnthalpyAdjointError::InvalidInput(
                "SDIRK2 transpose work budget overflow",
            ))?;
        let second = self.operator(1).pullback(cx, seeds[1], config)?;
        let mut first_seed = Vec::with_capacity(self.masses.len());
        for (vertex, (&seed, &lambda)) in seeds[0]
            .iter()
            .zip(&second.adjoint)
            .enumerate()
        {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let carry = if self.stages[1].prescribed[vertex] {
                0.0
            } else {
                finite(self.stages[0].storage_diagonal[vertex] * lambda)?
            };
            first_seed.push(finite(self.carry_ratio.mul_add(carry, seed))?);
        }
        let first = self.operator(0).pullback(cx, &first_seed, config)?;
        let mut previous_specific_enthalpy = Vec::with_capacity(self.masses.len());
        for (vertex, (&first, &second)) in first
            .previous_specific_enthalpy
            .iter()
            .zip(&second.previous_specific_enthalpy)
            .enumerate()
        {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            previous_specific_enthalpy
                .push(finite((1.0 - self.carry_ratio).mul_add(second, first))?);
        }
        poll(cx, self.masses.len())?;
        Ok(EnthalpySdirk2Gradient {
            previous_specific_enthalpy,
            source_density: [first.source_density, second.source_density],
            prescribed_temperature: [first.prescribed_temperature, second.prescribed_temperature],
            stage_adjoint: [first.adjoint, second.adjoint],
            relative_residual: [first.relative_residual, second.relative_residual],
            iterations: [first.iterations, second.iterations],
        })
    }
}

// Local contraction of the admitted joule-residual multiplier. Keeping the
// prescribed-row mask at its own physical stage is essential when baths change.
fn stage_conductivity_contraction(
    diagonal_dt: f64,
    loads: [f64; 4],
    multipliers: [f64; 4],
    prescribed: [bool; 4],
) -> Option<f64> {
    let mut contraction = 0.0_f64;
    for row in 0..4 {
        if !prescribed[row] {
            contraction = multipliers[row].mul_add(loads[row], contraction);
            if !contraction.is_finite() {
                return None;
            }
        }
    }
    let value = -diagonal_dt * contraction;
    value.is_finite().then_some(value)
}

fn radiation_entries(
    problem: ConductionProblem<'_>,
    radiation: Option<PointwiseBoundary<'_>>,
) -> Result<usize, EnthalpyAdjointError> {
    let workspace = |patches, bands| {
        pointwise::workspace_entries(
            problem.mesh.vertex_count(),
            problem.mesh.boundary().len(),
            patches,
        )?
        .checked_mul(bands)
        .ok_or(EnthalpyAdjointError::InvalidInput(
            "SDIRK2 spectral radiation workspace overflow",
        ))
    };
    match radiation {
        None => Ok(0),
        Some(PointwiseBoundary::Ambient(owner)) => workspace(owner.patch_count(), 1),
        Some(PointwiseBoundary::Enclosure(owner)) => {
            workspace(owner.patch_count(), owner.band_count())
        }
        Some(PointwiseBoundary::MaterialAmbient(owner)) => Ok(owner.adjoint_workspace_entries()?),
        Some(PointwiseBoundary::MaterialEnclosure(owner)) => Ok(owner.adjoint_workspace_entries()?),
    }
}

