//! First-order force/contact/drift composition on the ORIGINAL numerical owners.
//!
//! Hyper3 evaluates original elastic, held referential and uniform closed-pressure forces.
//! The existing impulse owner applies h(F_external - F_internal) using the full
//! supported consistent mass. A caller then resolves simultaneous contact at
//! this unchanged configuration. Finally the original generalized-alpha owner
//! advances the mass-only kinetic subflow, with zero subflow acceleration and
//! zero resistance: q_next = q + h v_contact, v_next = v_contact.
//!
//! The retained native history entry is the ACTUAL implicit solve of that
//! kinetic subflow; it is not a full-solid midpoint-equilibrium receipt. The
//! endpoint acceleration is subsequently recomputed by the original mass solve
//! and original elastic/current-pressure forces. It describes the smooth force component only;
//! contact impulses remain measures and are not disguised as smooth forces.
//!
//! The caller must admit the whole coupled geometric path, including omitted
//! contact features and element orientation. These staged primitives provide no
//! collision or Moreau-Jean existence claim and no contact adjoint. Explicit
//! elastic splitting requires timestep refinement; its signed and noncancelling
//! work defects are reported and gated, never allocated to physical heat.

use std::sync::Arc;

use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};

use super::impulse::{
    FiniteSolidImpulse, FiniteSolidImpulseControls, FiniteSolidMomentum, ImpulseMass, moments,
};
use super::{
    Cx, DynamicOperator, FiniteSolidBody, FiniteSolidControls, FiniteSolidError,
    FiniteSolidMaterial, ImplicitStepTelemetry, MassDiagonal, OperatorGeneralizedAlpha,
    PreparedHyperTetProblem, SolveReport, checkpoint, dot, inflate, second_order_acceleration,
};

/// Native force kick and immutable binding to its original configuration/load.
/// The receipt may cross any number of instantaneous contact impulses, but not
/// another configuration, material binding, physical time or accepted step.
#[derive(Clone, Debug)]
pub struct FiniteSolidForceKick<'a> {
    pub duration_s: f64,
    /// Original consistent-mass impulse receipt, including stationary supports.
    pub impulse: FiniteSolidImpulse,
    /// Held external force times h, at the actual initial current positions.
    pub source_impulse: FiniteSolidMomentum,
    /// Minus original elastic force times h. Objectivity requires zero resultant
    /// and moment; their measured residuals pass the requested impulse gates.
    pub internal_impulse: FiniteSolidMomentum,
    /// Total held-load kick work, including the separately retained pressure part.
    pub source_impulse_work_j: f64,
    /// Actual pressure impulse paired with the mean velocity across this kick.
    pub pressure_impulse_work_j: f64,
    pub internal_impulse_work_j: f64,
    /// Native combined-impulse work minus separately accumulated source and
    /// internal work. This measured floating-point decomposition defect is
    /// gated independently and retained in the complete splitting account.
    pub work_decomposition_residual_j: f64,
    /// Native combined nodal-impulse resultant minus separate source/internal
    /// resultants. Both momentum decompositions have independent absolute gates.
    pub linear_decomposition_residual_kg_m_s: [f64; 3],
    pub angular_decomposition_residual_kg_m2_s: [f64; 3],
    pub strain_before_j: f64,
    prepared: Arc<PreparedHyperTetProblem<'a, FiniteSolidMaterial>>,
    initial_q: Vec<f64>,
    initial_t: f64,
    initial_steps: usize,
    initial_history_len: usize,
    initial_last_step: Option<ImplicitStepTelemetry>,
    external_force_n: Vec<f64>,
    pressure: Option<(Arc<ClosedTetPressureSurface>, f64)>,
    controls: FiniteSolidImpulseControls,
    // Public diagnostic fields may be edited by callers; only this immutable
    // binding supplies the numerical continuation and its physical accounting.
    binding: ForceKickBinding,
}

#[derive(Clone, Debug)]
struct ForceKickBinding {
    duration_s: f64,
    origin_m: [f64; 3],
    source_work_j: f64,
    pressure_work_j: f64,
    internal_work_j: f64,
    strain_before_j: f64,
    kick_work_residual_j: f64,
    decomposition_residual_j: f64,
}

/// Original kinetic-subflow motion and complete force-splitting accounting.
#[derive(Clone, Debug)]
pub struct FiniteSolidContactDrift {
    /// The ACTUAL native solve of the mass-only kinetic subflow, retained in
    /// SecondOrderState.history. This is not a full elastic midpoint step.
    pub time: ImplicitStepTelemetry,
    /// Original mass solve for smooth endpoint acceleration under the same held
    /// loads and endpoint elastic/pressure forces; excludes the contact impulse measure.
    pub endpoint_acceleration: SolveReport,
    pub before: FiniteSolidMomentum,
    pub after: FiniteSolidMomentum,
    /// Actual drift increments, evaluated without subtracting large momenta.
    pub linear_residual_kg_m_s: [f64; 3],
    pub angular_residual_kg_m2_s: [f64; 3],
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Average(v) dot M delta-v during kinetic drift, normally exactly zero.
    pub kinetic_change_j: f64,
    pub strain_before_j: f64,
    pub strain_after_j: f64,
    /// Held referential load work plus original uniform-pressure -p Delta V work.
    pub source_displacement_work_j: f64,
    /// Original closed-pressure work, already included in total source work.
    pub pressure_work_j: f64,
    /// All external kick work minus all external displacement work.
    pub source_splitting_work_j: f64,
    /// Pressure kick work minus -p Delta V, gated separately from dead-load work.
    pub pressure_splitting_work_j: f64,
    /// Strain-energy change plus minus-elastic-force kick work.
    pub elastic_splitting_work_j: f64,
    /// Sum of the two signed split defects, actual kick work residual and drift
    /// kinetic change. Contact storage/loss belongs to the separate impact map.
    pub integration_defect_j: f64,
    /// Sum of their absolute magnitudes, separating dead-load and pressure work.
    /// Opposing source/pressure/elastic errors cannot
    /// cancel this per-solid admission gate.
    pub absolute_integration_defect_j: f64,
}

impl<'a> FiniteSolidBody<'a> {
    /// Stage a force kick without advancing q, physical time or native history.
    ///
    /// Held loads and original internal forces act at the initial current
    /// configuration. This primitive must be followed by the contact map and
    /// `contact_drift_candidate`; it is not a complete mechanical timestep.
    /// Pressure-bound bodies retain their explicit-pressure refusal on this entry.
    #[allow(clippy::too_many_arguments)]
    pub fn force_kick_candidate(
        &self,
        duration_s: f64,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        origin_m: [f64; 3],
        controls: FiniteSolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidForceKick<'a>), FiniteSolidError> {
        self.force_kick_with_pressure(
            duration_s,
            nodal_forces_n,
            body_force_n_m3,
            None,
            origin_m,
            controls,
            cx,
        )
    }

    /// Stage the original first-order force kick with uniform pressure held on
    /// the bound closed material surface. Current pressure forces act at the
    /// same old geometry as the contact impulses. The matching drift retains
    /// original -p Delta V work and bounds pressure splitting error separately.
    /// This is not the midpoint energy/momentum map or an exact pressure flow.
    #[allow(clippy::too_many_arguments)]
    pub fn force_kick_candidate_with_pressure(
        &self,
        duration_s: f64,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_pa: f64,
        origin_m: [f64; 3],
        controls: FiniteSolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidForceKick<'a>), FiniteSolidError> {
        self.force_kick_with_pressure(
            duration_s,
            nodal_forces_n,
            body_force_n_m3,
            Some(pressure_pa),
            origin_m,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn force_kick_with_pressure(
        &self,
        duration_s: f64,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        uniform_pressure_pa: Option<f64>,
        origin_m: [f64; 3],
        controls: FiniteSolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidForceKick<'a>), FiniteSolidError> {
        checkpoint(cx)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(invalid("force kick requires a finite positive interval"));
        }
        let pressure = self.admit_split_pressure(uniform_pressure_pa, cx)?;
        let displacement = inflate(&self.state.q, &self.free, &self.fixed);
        let strain_before_j = self.prepared.stored_energy_j(&displacement, cx)?;
        let external = self
            .prepared
            .external_force_n(nodal_forces_n, body_force_n_m3, cx)?;
        let internal = self.prepared.internal_force(&displacement, cx)?;
        let pressure_impulse = pressure
            .as_ref()
            .map(|(surface, value)| {
                let mut impulse = surface.force_n(&displacement, *value, cx)?;
                for (dof, value) in impulse.iter_mut().enumerate() {
                    if dof % 4096 == 0 {
                        checkpoint(cx)?;
                    }
                    *value *= duration_s;
                }
                finite(&impulse)?;
                Ok::<_, FiniteSolidError>(impulse)
            })
            .transpose()?;
        let mut source = vec![0.0; external.len()];
        let mut elastic = vec![0.0; external.len()];
        let mut nodal_impulses = Vec::with_capacity(self.problem.nodes_m.len());
        for node in 0..self.problem.nodes_m.len() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            let mut impulse = [0.0; 3];
            for (axis, value) in impulse.iter_mut().enumerate() {
                let dof = 3 * node + axis;
                source[dof] = duration_s * external[dof];
                elastic[dof] = -duration_s * internal[dof];
                let applied_source = match &pressure_impulse {
                    Some(pressure) => source[dof] + pressure[dof],
                    None => source[dof],
                };
                *value = applied_source + elastic[dof];
            }
            finite(&impulse)?;
            nodal_impulses.push((node, impulse));
        }
        finite(&source)?;
        finite(&elastic)?;
        let (candidate, impulse) =
            self.nodal_impulse_candidate(&nodal_impulses, origin_m, controls, cx)?;
        let positions = self.current_positions(cx)?;
        let mut source_impulse = moments(&positions, &source, origin_m, cx)?;
        if let Some(pressure) = &pressure_impulse {
            let momentum = moments(&positions, pressure, origin_m, cx)?;
            if magnitude(momentum.linear_kg_m_s) > controls.linear_tolerance_kg_m_s
                || magnitude(momentum.angular_kg_m2_s) > controls.angular_tolerance_kg_m2_s
            {
                return Err(FiniteSolidError::ImpulseBalance {
                    linear_residual_kg_m_s: momentum.linear_kg_m_s,
                    angular_residual_kg_m2_s: momentum.angular_kg_m2_s,
                    work_residual_j: 0.0,
                });
            }
            for axis in 0..3 {
                source_impulse.linear_kg_m_s[axis] += momentum.linear_kg_m_s[axis];
                source_impulse.angular_kg_m2_s[axis] += momentum.angular_kg_m2_s[axis];
            }
            finite(&source_impulse.linear_kg_m_s)?;
            finite(&source_impulse.angular_kg_m2_s)?;
        }
        let internal_impulse = moments(&positions, &elastic, origin_m, cx)?;
        if magnitude(internal_impulse.linear_kg_m_s) > controls.linear_tolerance_kg_m_s
            || magnitude(internal_impulse.angular_kg_m2_s) > controls.angular_tolerance_kg_m2_s
        {
            return Err(FiniteSolidError::ImpulseBalance {
                linear_residual_kg_m_s: internal_impulse.linear_kg_m_s,
                angular_residual_kg_m2_s: internal_impulse.angular_kg_m2_s,
                work_residual_j: 0.0,
            });
        }
        let mut average_v = vec![0.0; external.len()];
        for ((&dof, &before), &after) in self.free.iter().zip(&self.state.v).zip(&candidate.state.v)
        {
            average_v[dof] = 0.5 * before + 0.5 * after;
        }
        let dead_load_impulse_work_j = dot(&source, &average_v);
        let pressure_impulse_work_j = pressure_impulse
            .as_ref()
            .map_or(0.0, |pressure| dot(pressure, &average_v));
        let source_impulse_work_j = if pressure_impulse.is_some() {
            dead_load_impulse_work_j + pressure_impulse_work_j
        } else {
            dead_load_impulse_work_j
        };
        let internal_impulse_work_j = dot(&elastic, &average_v);
        let work_decomposition_residual_j =
            impulse.impulse_work_j - (source_impulse_work_j + internal_impulse_work_j);
        let linear_decomposition_residual_kg_m_s = std::array::from_fn(|axis| {
            impulse.applied.linear_kg_m_s[axis]
                - (source_impulse.linear_kg_m_s[axis] + internal_impulse.linear_kg_m_s[axis])
        });
        let angular_decomposition_residual_kg_m2_s = std::array::from_fn(|axis| {
            impulse.applied.angular_kg_m2_s[axis]
                - (source_impulse.angular_kg_m2_s[axis] + internal_impulse.angular_kg_m2_s[axis])
        });
        finite(&[
            source_impulse_work_j,
            pressure_impulse_work_j,
            dead_load_impulse_work_j,
            internal_impulse_work_j,
            work_decomposition_residual_j,
        ])?;
        finite(&linear_decomposition_residual_kg_m_s)?;
        finite(&angular_decomposition_residual_kg_m2_s)?;
        if work_decomposition_residual_j.abs() > controls.work_tolerance_j
            || magnitude(linear_decomposition_residual_kg_m_s) > controls.linear_tolerance_kg_m_s
            || magnitude(angular_decomposition_residual_kg_m2_s)
                > controls.angular_tolerance_kg_m2_s
        {
            return Err(FiniteSolidError::ImpulseBalance {
                linear_residual_kg_m_s: linear_decomposition_residual_kg_m_s,
                angular_residual_kg_m2_s: angular_decomposition_residual_kg_m2_s,
                work_residual_j: work_decomposition_residual_j,
            });
        }
        let binding = ForceKickBinding {
            duration_s,
            origin_m,
            source_work_j: dead_load_impulse_work_j,
            pressure_work_j: pressure_impulse_work_j,
            internal_work_j: internal_impulse_work_j,
            strain_before_j,
            kick_work_residual_j: impulse.work_residual_j,
            decomposition_residual_j: work_decomposition_residual_j,
        };
        let report = FiniteSolidForceKick {
            duration_s,
            impulse,
            source_impulse,
            internal_impulse,
            source_impulse_work_j,
            pressure_impulse_work_j,
            internal_impulse_work_j,
            work_decomposition_residual_j,
            linear_decomposition_residual_kg_m_s,
            angular_decomposition_residual_kg_m2_s,
            strain_before_j,
            prepared: Arc::clone(&self.prepared),
            initial_q: self.state.q.clone(),
            initial_t: self.state.t,
            initial_steps: self.state.steps,
            initial_history_len: self.state.history.len(),
            initial_last_step: self.state.history.last().cloned(),
            external_force_n: external,
            pressure,
            controls,
            binding,
        };
        checkpoint(cx)?;
        Ok((candidate, report))
    }

    /// Admit the original pressure binding before a trajectory starts contact
    /// work. Zero pressure takes the unchanged no-pressure numerical path.
    pub(in crate::finite_solid) fn admit_split_pressure(
        &self,
        pressure_pa: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<Option<(Arc<ClosedTetPressureSurface>, f64)>, FiniteSolidError> {
        match pressure_pa {
            Some(value) => {
                let surface = self
                    .pressure_surface
                    .as_ref()
                    .ok_or_else(|| invalid("body was not initialized with a pressure surface"))?;
                if surface.validate_pressure(TetPressure::Uniform(value), cx)? {
                    Ok(None)
                } else {
                    Ok(Some((Arc::clone(surface), value)))
                }
            }
            None if self.pressure_surface.is_some() => Err(invalid(
                "contact splitting requires explicit pressure support",
            )),
            None => Ok(None),
        }
    }

    /// Stage the original kinetic subflow after the matching kick/contact map.
    ///
    /// The matching receipt fixes the material, initial configuration, duration,
    /// held loads, mass solve controls and accounting origin. Arbitrary velocity
    /// changes between kick and drift are the caller's separately audited contact
    /// impulses. No configuration projection or implicit full-solid solve occurs.
    /// A successful candidate still needs the caller's whole-path admission.
    pub fn contact_drift_candidate(
        &self,
        kick: &FiniteSolidForceKick<'a>,
        controls: FiniteSolidControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidContactDrift), FiniteSolidError> {
        checkpoint(cx)?;
        if !Arc::ptr_eq(&self.prepared, &kick.prepared)
            || self.state.q != kick.initial_q
            || self.state.t != kick.initial_t
            || self.state.steps != kick.initial_steps
            || self.state.history.len() != kick.initial_history_len
            || self.state.history.last() != kick.initial_last_step.as_ref()
        {
            return Err(invalid(
                "contact drift requires its original unchanged kick configuration",
            ));
        }
        if !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
            || controls.solve.max_newton_iterations == 0
            || self.state.steps >= controls.maximum_steps
        {
            return Err(invalid("invalid kinetic drift controls or step budget"));
        }
        let binding = &kick.binding;
        let origin = binding.origin_m;
        let before_u = inflate(&self.state.q, &self.free, &self.fixed);
        let before_positions = self.current_positions(cx)?;
        let before = self.momentum_about(origin, cx)?;
        let kinetic_before_j = self.kinetic_energy_j()?;
        let mut candidate = self.clone();
        // Acceleration belongs to this declared zero-force subflow only. The
        // actual smooth endpoint acceleration is restored by its native solve
        // below; the accepted input and prior telemetry are never changed.
        candidate.state.a.fill(0.0);
        let zero = vec![0.0; self.free.len()];
        let integrator =
            OperatorGeneralizedAlpha::new(self.free.len(), binding.duration_s, 1.0, controls.solve);
        let time = integrator.step_controlled(
            &mut candidate.state,
            &ImpulseMass {
                mass: &self.free_mass,
            },
            &zero,
            &mut || cx.checkpoint().is_err(),
        )?;
        let after_u = inflate(&candidate.state.q, &self.free, &self.fixed);
        let strain_after_j = self.prepared.stored_energy_j(&after_u, cx)?;
        let forcing: Vec<f64> = self
            .free
            .iter()
            .map(|&dof| kick.external_force_n[dof])
            .collect();
        let operator = DynamicOperator::new(
            &self.prepared,
            &self.free_mass,
            &self.free,
            &self.fixed,
            kick.pressure
                .as_ref()
                .map(|(surface, value)| (surface.as_ref(), TetPressure::Uniform(*value))),
            cx,
        );
        let acceleration = second_order_acceleration(
            &operator,
            &candidate.state.q,
            &candidate.state.v,
            &forcing,
            &MassDiagonal {
                mass: &self.free_mass,
            },
            kick.controls.mass_solve.solve,
            kick.controls.mass_solve.maximum_workspace_components,
            &mut || cx.checkpoint().is_err(),
        );
        if let Some(error) = operator.error.take() {
            return Err(error.into());
        }
        let acceleration = acceleration?;
        candidate.state.a = acceleration.value;

        let after_positions = candidate.current_positions(cx)?;
        let mut before_v = vec![0.0; self.mass.nrows()];
        let mut delta_v = vec![0.0; self.mass.nrows()];
        let mut average_v = vec![0.0; self.mass.nrows()];
        for ((&dof, &old), &new) in self.free.iter().zip(&self.state.v).zip(&candidate.state.v) {
            before_v[dof] = old;
            delta_v[dof] = new - old;
            average_v[dof] = 0.5 * old + 0.5 * new;
        }
        let mut before_p = vec![0.0; self.mass.nrows()];
        let mut delta_p = vec![0.0; self.mass.nrows()];
        self.mass.spmv(&before_v, &mut before_p);
        self.mass.spmv(&delta_v, &mut delta_p);
        let average_p: Vec<f64> = before_p
            .iter()
            .zip(&delta_p)
            .map(|(&p, &dp)| p + 0.5 * dp)
            .collect();
        let average_positions: Vec<[f64; 3]> = before_positions
            .iter()
            .zip(&after_positions)
            .map(|(old, new)| std::array::from_fn(|axis| 0.5 * old[axis] + 0.5 * new[axis]))
            .collect();
        let delta_positions: Vec<[f64; 3]> = before_positions
            .iter()
            .zip(&after_positions)
            .map(|(old, new)| std::array::from_fn(|axis| new[axis] - old[axis]))
            .collect();
        let momentum_change = moments(&average_positions, &delta_p, origin, cx)?;
        let position_change = moments(&delta_positions, &average_p, [0.0; 3], cx)?;
        let angular_residual_kg_m2_s: [f64; 3] = std::array::from_fn(|axis| {
            momentum_change.angular_kg_m2_s[axis] + position_change.angular_kg_m2_s[axis]
        });
        let kinetic_change_j = dot(&average_v, &delta_p);
        if magnitude(momentum_change.linear_kg_m_s) > kick.controls.linear_tolerance_kg_m_s
            || magnitude(angular_residual_kg_m2_s) > kick.controls.angular_tolerance_kg_m2_s
            || kinetic_change_j.abs() > kick.controls.work_tolerance_j
        {
            return Err(FiniteSolidError::ImpulseBalance {
                linear_residual_kg_m_s: momentum_change.linear_kg_m_s,
                angular_residual_kg_m2_s,
                work_residual_j: kinetic_change_j,
            });
        }
        let mut source_displacement_work_j = 0.0;
        for (dof, ((&old, &new), &force)) in before_u
            .iter()
            .zip(&after_u)
            .zip(&kick.external_force_n)
            .enumerate()
        {
            if dof % 4096 == 0 {
                checkpoint(cx)?;
            }
            source_displacement_work_j = force.mul_add(new - old, source_displacement_work_j);
        }
        let dead_load_splitting_work_j = binding.source_work_j - source_displacement_work_j;
        let pressure_work_j = match &kick.pressure {
            Some((surface, value)) => surface.work_j(&before_u, &after_u, *value, cx)?,
            None => 0.0,
        };
        let pressure_splitting_work_j = binding.pressure_work_j - pressure_work_j;
        let source_splitting_work_j = if kick.pressure.is_some() {
            source_displacement_work_j += pressure_work_j;
            dead_load_splitting_work_j + pressure_splitting_work_j
        } else {
            dead_load_splitting_work_j
        };
        let elastic_splitting_work_j =
            (strain_after_j - binding.strain_before_j) + binding.internal_work_j;
        let integration_defect_j = source_splitting_work_j
            + elastic_splitting_work_j
            + binding.kick_work_residual_j
            + binding.decomposition_residual_j
            + kinetic_change_j;
        let absolute_integration_defect_j = dead_load_splitting_work_j.abs()
            + pressure_splitting_work_j.abs()
            + elastic_splitting_work_j.abs()
            + binding.kick_work_residual_j.abs()
            + binding.decomposition_residual_j.abs()
            + kinetic_change_j.abs();
        finite(&[
            kinetic_change_j,
            source_displacement_work_j,
            pressure_work_j,
            source_splitting_work_j,
            pressure_splitting_work_j,
            dead_load_splitting_work_j,
            elastic_splitting_work_j,
            integration_defect_j,
            absolute_integration_defect_j,
        ])?;
        finite(&angular_residual_kg_m2_s)?;
        if absolute_integration_defect_j > controls.energy_tolerance_j {
            return Err(FiniteSolidError::Energy {
                defect_j: absolute_integration_defect_j,
                allowance_j: controls.energy_tolerance_j,
            });
        }
        let report = FiniteSolidContactDrift {
            time,
            endpoint_acceleration: acceleration.report,
            before,
            after: candidate.momentum_about(origin, cx)?,
            linear_residual_kg_m_s: momentum_change.linear_kg_m_s,
            angular_residual_kg_m2_s,
            kinetic_before_j,
            kinetic_after_j: candidate.kinetic_energy_j()?,
            kinetic_change_j,
            strain_before_j: binding.strain_before_j,
            strain_after_j,
            source_displacement_work_j,
            pressure_work_j,
            source_splitting_work_j,
            pressure_splitting_work_j,
            elastic_splitting_work_j,
            integration_defect_j,
            absolute_integration_defect_j,
        };
        checkpoint(cx)?;
        Ok((candidate, report))
    }
}

fn invalid(message: &'static str) -> FiniteSolidError {
    FiniteSolidError::Invalid(message)
}

fn finite(values: &[f64]) -> Result<(), FiniteSolidError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(invalid("solid contact splitting is unrepresentable"))
    }
}

fn magnitude(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}
