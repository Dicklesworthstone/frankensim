//! Finite load jumps and explicit material reconfiguration at accepted times.
//! A load jump preserves velocity. A material rebind preserves REDUCED
//! canonical momentum and reports the work required by the new Hamiltonian.
use super::*;
use crate::linear3::TetMaterialField;

/// Independent physical admission limits for a zero-duration material change.
#[derive(Debug, Clone, Copy)]
pub struct MaterialRebindConfig3 {
    /// Bound on the absolute signed kinetic plus shifted-elastic work [J].
    /// Zero admits only an exactly zero-work replacement.
    pub maximum_absolute_work_j: f64,
    /// Absolute free canonical-momentum residual allowance [kg m/s].
    /// The run's relative linear-solve tolerance also applies to momentum norm.
    pub absolute_momentum_tolerance_kg_m_s: f64,
    /// Bound on both total and elastic Mandel strains, in (0,1).
    pub maximum_strain: f64,
}

/// Material replacement work at fixed displacement and accepted physical time.
/// Every two-entry array is ordered [before, after].
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialRebind3 {
    pub time_s: f64,
    pub step: usize,
    pub reference_mass_kg: [f64; 2],
    pub kinetic_energy_j: [f64; 2],
    pub elastic_energy_j: [f64; 2],
    pub kinetic_work_j: f64,
    pub elastic_work_j: f64,
    /// Signed work supplied TO the solid by the external reconfiguration.
    /// Negative work must be received by that owner; neither sign is heat.
    pub supplied_work_j: f64,
    /// Norm of M_after_ff*v_after - M_before_ff*v_before.
    /// Stationary constraints may exchange FULL nodal momentum. This is only
    /// conservation of the free generalized coordinates' canonical momentum.
    pub reduced_momentum_residual_kg_m_s: f64,
    /// Original tetrahedron order; each pair retains the exact source binding.
    pub material_state_identities: Vec<[[u8; 32]; 2]>,
    pub thermal_state_identities: Vec<[Option<ContentHash>; 2]>,
}

impl<'g, 'm> TetTransientRun3<'g, 'm, '_, '_> {
    /// Replace source-resolved density, stiffness and stress-free strain at the
    /// SAME accepted time, retaining displacement and p_free = M_ff*v_free.
    /// The native initializer rebuilds the actual operators, solves velocity
    /// from that momentum, then re-equilibrates acceleration under the explicit
    /// new load and thermal field. The clock and consumed step count are kept.
    ///
    /// This is externally supplied material reconfiguration, not mass transport,
    /// continuous parameter evolution or a caloric heat transfer. Required work
    /// is Delta(1/2 v'Mv) + Delta(sum V/2 e'C e), e = strain - strain_free.
    /// With stationary constraints, preserving reduced momentum does NOT imply
    /// conservation of full nodal or total physical momentum; the constraints
    /// can exchange momentum. No support impulse is inferred or reported here.
    ///
    /// The reference mesh, connectivity and fixed DOFs must match exactly.
    /// Moving supports refuse. The returned run borrows the actual new problem;
    /// later motion and stress recovery consume its new sources and operators.
    /// Any source, solve, momentum, strain, work or cancellation refusal leaves
    /// the original run unchanged. This event consumes no advancement slot.
    #[allow(clippy::too_many_arguments)]
    pub fn rebind_material_candidate<'p, 's>(
        &self, problem_after: &'p TetLinearElasticProblem<'s>,
        thermal_after: Option<TetThermalStrainField<'_>>,
        load_after: &TetLoadCase3<'g, 'm>, mesh_frame: ContentHash,
        controls: MaterialRebindConfig3, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, MaterialRebind3), TransientError3> {
        check(cx)?;
        if !controls.maximum_absolute_work_j.is_finite() || controls.maximum_absolute_work_j < 0.0
            || !controls.absolute_momentum_tolerance_kg_m_s.is_finite()
            || controls.absolute_momentum_tolerance_kg_m_s <= 0.0
            || !controls.maximum_strain.is_finite() || controls.maximum_strain <= 0.0
            || controls.maximum_strain >= 1.0 {
            return Err(TransientError3::Invalid("finite material work, positive momentum and small-strain allowances required"));
        }
        if self.prescribed.is_some() || problem_after.nodes_m != self.problem.nodes_m
            || problem_after.tetrahedra != self.problem.tetrahedra
            || problem_after.fixed_dofs != self.problem.fixed_dofs {
            return Err(TransientError3::Invalid("material reconfiguration requires the same mesh and stationary zero supports"));
        }
        let momentum = self.apply(1.0, 0.0, &self.state.velocity, cx)?;
        let motion = self.nodal_motion(cx)?;
        let before = self.recover_response(mesh_frame, controls.maximum_strain, cx)?;
        let mut candidate = load_after.begin_transient_impl(
            problem_after, thermal_after, self.state.time_s, &motion[0], &motion[1],
            None, Some(&momentum), self.config, cx)?;
        if candidate.free_dofs() != self.free_dofs()
            || candidate.state.displacement != self.state.displacement {
            return Err(TransientError3::Invalid("material reconfiguration changed the physical displacement coordinates"));
        }
        candidate.state.step = self.state.step;
        let next_momentum = candidate.apply(1.0, 0.0, &candidate.state.velocity, cx)?;
        let mut residual = 0.0_f64;
        for (&old, &next) in momentum.iter().zip(&next_momentum) {
            check(cx)?; residual = number(residual.hypot(number(next-old)?))?;
        }
        let allowance = number(controls.absolute_momentum_tolerance_kg_m_s
            + self.config.solve.relative_residual_tolerance*number(length(&momentum))?)?;
        if residual > allowance {
            return Err(TransientError3::Momentum {
                residual_kg_m_s: residual, allowance_kg_m_s: allowance });
        }
        let after = candidate.recover_response(mesh_frame, controls.maximum_strain, cx)?;
        let kinetic = [self.state.kinetic_energy_j, candidate.state.kinetic_energy_j];
        let elastic = [before.stress.elastic_energy_j(), after.stress.elastic_energy_j()];
        let kinetic_work = number(kinetic[1]-kinetic[0])?;
        let elastic_work = number(elastic[1]-elastic[0])?;
        let supplied_work = number(kinetic_work+elastic_work)?;
        if supplied_work.abs() > controls.maximum_absolute_work_j {
            return Err(TransientError3::Invalid("material reconfiguration exceeds its external work allowance"));
        }
        let mut materials = Vec::with_capacity(self.problem.tetrahedra.len());
        let mut thermal = Vec::with_capacity(self.problem.tetrahedra.len());
        for element in 0..self.problem.tetrahedra.len() {
            check(cx)?;
            let material_id = |field: TetMaterialField<'_>| match field {
                TetMaterialField::Uniform(material) => material.material_state_identity,
                TetMaterialField::PerElement(values) => values[element].material_state_identity,
            };
            let thermal_id = |field: Option<TetThermalStrainField<'_>>| match field {
                None => None,
                Some(TetThermalStrainField::Uniform(state)) => Some(state.identity()),
                Some(TetThermalStrainField::PerElement(values)) => Some(values[element].identity()),
            };
            materials.push([material_id(self.problem.materials), material_id(problem_after.materials)]);
            thermal.push([thermal_id(self.thermal_strain_field()), thermal_id(candidate.thermal_strain_field())]);
        }
        let report = MaterialRebind3 { time_s: self.state.time_s, step: self.state.step,
            reference_mass_kg: [self.operators.total_mass_kg, candidate.operators.total_mass_kg],
            kinetic_energy_j: kinetic, elastic_energy_j: elastic,
            kinetic_work_j: kinetic_work, elastic_work_j: elastic_work, supplied_work_j: supplied_work,
            reduced_momentum_residual_kg_m_s: residual,
            material_state_identities: materials, thermal_state_identities: thermal };
        check(cx)?;
        Ok((candidate, report))
    }
}

/// Diagnostics for changing a finite force at the current accepted time.
/// No time, step, displacement, velocity or energy is changed by this operation.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadJump3 {
    pub time_s: f64,
    pub iterations: usize,
    pub force_residual_n: f64,
    /// Euclidean norm over free xyz acceleration components; not a time error.
    pub acceleration_change_l2_m_s2: f64,
}

/// One original Newmark interval to the event's LEFT load, followed by the
/// zero-duration switch to its RIGHT load. Both operations publish atomically.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientEventStep3 {
    pub interval: TransientStep3,
    pub jump: LoadJump3,
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Apply a finite load discontinuity at the CURRENT time. Solve
    /// M*a_plus = f_plus + f_thermal - C*v - K*u using the original mass solve.
    /// A finite jump changes acceleration, not momentum: q/v and both energy
    /// fields are preserved bit-for-bit. No impulse or artificial work is added.
    /// The full new load, including fixed-node forces, becomes the accepted
    /// response input and the next Newmark interval's starting load.
    ///
    /// No step slot is consumed because time does not advance. Work remains
    /// bounded by the existing linear-iteration allowance. Source/solve failure
    /// or cancellation publishes neither the new force nor acceleration.
    pub fn apply_load_jump(&mut self, load_after: &TetLoadCase3<'_, '_>, cx: &Cx<'_>)
        -> Result<LoadJump3, TransientError3> {
        check(cx)?;
        if self.prescribed.is_some() {
            return Err(TransientError3::Invalid("finite load jumps with prescribed support motion are not admitted"));
        }
        let force = self.reduce_force(load_after, cx)?;
        let ku = self.apply(0.0, 1.0, &self.state.displacement, cx)?;
        let cv = self.apply(self.config.mass_damping_per_s, self.config.stiffness_damping_s,
            &self.state.velocity, cx)?;
        let mut rhs = Vec::with_capacity(force.len());
        for i in 0..force.len() {
            if i % 256 == 0 { check(cx)?; }
            rhs.push(number(force[i] + self.thermal[i] - ku[i] - cv[i])?);
        }
        let (a, iterations, force_residual_n) = solve(&self.operator(1.0, 0.0, cx), &rhs, self.config)?;
        let mut change = 0.0_f64;
        for (i, (&next, &old)) in a.iter().zip(&self.state.acceleration).enumerate() {
            if i % 256 == 0 { check(cx)?; }
            change = number(change.hypot(number(next-old)?))?;
        }
        let mut full_force = Vec::with_capacity(load_after.nodal_forces_n().len());
        for &(d, f) in load_after.nodal_forces_n() { check(cx)?; full_force.push((d, f)); }
        let report = LoadJump3 { time_s: self.state.time_s, iterations, force_residual_n,
            acceleration_change_l2_m_s2: change };
        check(cx)?;
        self.state.acceleration = a;
        self.force = force;
        self.full_force = full_force;
        Ok(report)
    }

    /// Advance to an explicitly resolved event, then switch finite loads there.
    /// `load_before` and `load_after` are the physical LEFT and RIGHT limits at
    /// the SAME endpoint. Use the left load for the preceding interval; using
    /// the right load there smears a step backwards through the entire interval.
    /// Unchanged Newmark quadrature applies between the retained starting load
    /// and `load_before`. Any intervening discontinuities must be split by the
    /// caller. This does not discover knots, approximate left limits by epsilon,
    /// wrap frame clocks, or claim to detect unknown sub-step pulses.
    ///
    /// If either the interval or right-limit mass solve fails, neither is
    /// published. Forking shares the assembled operators and frozen thermal data.
    /// A successful event consumes exactly one ordinary advancement step slot.
    pub fn advance_to_load_jump(&mut self, load_before: &TetLoadCase3<'_, '_>,
        load_after: &TetLoadCase3<'_, '_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientEventStep3, TransientError3> {
        check(cx)?;
        let mut candidate = self.clone();
        let interval = candidate.advance(load_before, h_s, cx)?;
        let jump = candidate.apply_load_jump(load_after, cx)?;
        check(cx)?;
        *self = candidate;
        Ok(TransientEventStep3 { interval, jump })
    }
}
