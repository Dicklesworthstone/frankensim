//! Original Newmark motion and independently reconstructed P1 traction work.
use super::*;
use crate::linear3::{TetLinearElasticProblem, TetThermalStrainField};
use crate::loads3::transient::{
    TetTransientRun3, TransientConfig3, TransientError3, TransientStep3,
};

/// Subdivisions of the native mechanical-work account, never an extra debit.
/// Includes pressure work on moving constrained nodes. `native.support_work_j`
/// separately retains the work of support reactions ON the solid.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveP1Step3 {
    pub native: TransientStep3,
    pub transfer: AdhesiveP1Transfer3,
    /// Average endpoint original traction integrated against actual P1 travel.
    /// Explicit full-stick projection includes both tangential components.
    pub spatial_traction_work_j: f64,
    /// Average represented nodal adhesive correction dotted with that travel.
    pub traction_work_on_solid_j: f64,
    pub other_external_work_j: f64,
    pub work_residual_j: f64,
    pub work_tolerance_j: f64,
    pub quadrature_work_estimate_j: f64,
}

/// Absolute world-component support endpoint on the original fixed DOF set.
/// Components must be finite, sorted and unique; the native owner checks all
/// motion, solve and displacement bounds. Both variants use the same Newmark
/// support kinematics, consistent mass, damping and elastic coupling.
#[derive(Clone, Copy, Debug)]
pub enum AdhesiveP1SupportEndpoint3<'a> {
    /// Absolute support displacement [m]. An omitted fixed DOF is zero at this
    /// endpoint. Support velocity and acceleration follow the Newmark map;
    /// sufficiently smooth motion and adequate time resolution are required.
    Displacement(&'a [(usize, f64)]),
    /// Endpoint support acceleration [m/s²]. An omitted fixed DOF has zero
    /// acceleration; its accepted displacement and velocity continue. No
    /// impulse or velocity-jump convention is inferred.
    Acceleration(&'a [(usize, f64)]),
}

/// Fixed coplanar pressure patch on the existing complete sparse Newmark owner.
/// Matrices remain shared by the native candidate clone. Projection storage is
/// reused; the native trajectory keeps its existing bounded vector allocation.
/// Supports may be stationary or explicitly prescribed. Pressure stays on the
/// reference material trace with its original fixed frame. Explicit full-stick
/// selection carries both tangential components through the same native force
/// balance and work account. Support movement does not create follower loads,
/// contact migration or a new adhesion law.
/// Read-only native access prevents an untracked load/support change; consume
/// this wrapper to continue elsewhere.
pub struct AdhesiveP1Transient3<'g, 'm, 'p, 's> {
    patch: AdhesiveP1Patch3<'g, 'm>,
    native: TetTransientRun3<'g, 'm, 'p, 's>,
    cells: Vec<SectorAdhesiveTraction>,
    vector_cells: Option<Vec<SectorAdhesiveVectorTraction>>,
    represented: Vec<Dd>,
    background: Vec<f64>,
    delta: Vec<f64>,
    native_config: TransientConfig3,
    current: AdhesiveP1Transfer3,
    last: Option<AdhesiveP1Step3>,
}

impl<'g, 'm> AdhesiveP1Patch3<'g, 'm> {
    /// Start the actual native tetrahedral trajectory with this original source
    /// in its INITIAL force/acceleration balance. Material, consistent mass,
    /// thermal strain, damping and supports enter the existing native owner.
    /// This constructor requires stationary zero supports. Use
    /// `begin_transient_with_supports` for explicitly prescribed motion.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient<'p, 's>(
        self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        config: TransientConfig3,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<AdhesiveP1Transient3<'g, 'm, 'p, 's>, AdhesiveP1Error3> {
        self.begin_impl(
            profile,
            other,
            problem,
            thermal,
            time_s,
            displacement_m,
            velocity_m_s,
            None,
            config,
            gate,
            cx,
        )
    }

    /// Start with prescribed supports in the ORIGINAL full mass/stiffness and
    /// force balance. Initial u/v include free and fixed components; supplied
    /// acceleration is nonzero only on fixed DOFs. Free acceleration is solved
    /// by the native owner, including M_fc, C_fc and K_fc coupling.
    ///
    /// Every subsequent step must supply an explicit support endpoint using
    /// `advance_with_supports` or its admission variant. The constrained DOF set
    /// and reference pressure trace remain frozen.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient_with_supports<'p, 's>(
        self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        prescribed_acceleration_m_s2: &[[f64; 3]],
        config: TransientConfig3,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<AdhesiveP1Transient3<'g, 'm, 'p, 's>, AdhesiveP1Error3> {
        self.begin_impl(
            profile,
            other,
            problem,
            thermal,
            time_s,
            displacement_m,
            velocity_m_s,
            Some(prescribed_acceleration_m_s2),
            config,
            gate,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn begin_impl<'p, 's>(
        mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        prescribed_acceleration_m_s2: Option<&[[f64; 3]]>,
        config: TransientConfig3,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<AdhesiveP1Transient3<'g, 'm, 'p, 's>, AdhesiveP1Error3> {
        let (load, current) = self.assemble(profile, other, gate, cx)?;
        let native = match prescribed_acceleration_m_s2 {
            Some(acceleration) => load.begin_transient_with_supports(
                problem,
                thermal,
                time_s,
                displacement_m,
                velocity_m_s,
                acceleration,
                config,
                cx,
            )?,
            None => load.begin_transient(
                problem,
                thermal,
                time_s,
                displacement_m,
                velocity_m_s,
                config,
                cx,
            )?,
        };
        let mut cells = vec![SectorAdhesiveTraction::default(); self.cells.len()];
        let mut represented = vec![Dd::ZERO; self.represented.len()];
        let mut background = vec![0.0; self.background.len()];
        let delta = vec![0.0; self.background.len()];
        std::mem::swap(&mut self.cells, &mut cells);
        std::mem::swap(&mut self.represented, &mut represented);
        std::mem::swap(&mut self.background, &mut background);
        let vector_cells = self.vector.as_mut().map(|scratch| {
            let mut accepted = vec![SectorAdhesiveVectorTraction::default(); scratch.cells.len()];
            std::mem::swap(&mut scratch.cells, &mut accepted);
            accepted
        });
        poll(gate, cx)?;
        Ok(AdhesiveP1Transient3 {
            patch: self,
            native,
            cells,
            vector_cells,
            represented,
            background,
            delta,
            native_config: config,
            current,
            last: None,
        })
    }
}

impl<'g, 'm, 'p, 's> AdhesiveP1Transient3<'g, 'm, 'p, 's> {
    pub const fn native(&self) -> &TetTransientRun3<'g, 'm, 'p, 's> {
        &self.native
    }
    pub fn into_native(self) -> TetTransientRun3<'g, 'm, 'p, 's> {
        self.native
    }
    pub const fn current_transfer(&self) -> &AdhesiveP1Transfer3 {
        &self.current
    }
    pub const fn last_step(&self) -> Option<&AdhesiveP1Step3> {
        self.last.as_ref()
    }

    /// Advance one pressure endpoint through the original Newmark scheme.
    /// Work uses average endpoint tractions on actual material displacement.
    /// Resolve source jumps and short pulses explicitly; changing a pressure
    /// endpoint is not an exact continuous adhesion trajectory or snap event.
    pub fn advance(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        step_s: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<&AdhesiveP1Step3, AdhesiveP1Error3> {
        let (candidate, report) = self.prepare(profile, other, None, step_s, gate, cx)?;
        self.publish(candidate, report);
        Ok(self.last.as_ref().unwrap())
    }

    /// Stage another receiver against actual complete native motion. Its result
    /// must be infallibly installable after success. Refusal or late cancellation
    /// preserves every accepted native coordinate, force endpoint and work report.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_admission<R, E>(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        step_s: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
        admit: impl FnOnce(&AdhesiveP1Step3, &TetTransientRun3<'g, 'm, 'p, 's>) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<AdhesiveP1Error3>,
    {
        let (candidate, report) = self
            .prepare(profile, other, None, step_s, gate, cx)
            .map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(gate, cx).map_err(E::from)?;
        self.publish(candidate, report);
        Ok(receiving)
    }

    /// Advance pressure and support motion together through the existing native
    /// Newmark owner. Full nodal travel includes constrained pressure nodes.
    /// Source-pressure work and support-reaction work remain separate physical
    /// inputs. The fixed reference patch does not follow a moving obstacle.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_supports(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        supports: AdhesiveP1SupportEndpoint3<'_>,
        step_s: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<&AdhesiveP1Step3, AdhesiveP1Error3> {
        let (candidate, report) = self.prepare(profile, other, Some(supports), step_s, gate, cx)?;
        self.publish(candidate, report);
        Ok(self.last.as_ref().unwrap())
    }

    /// Admit another receiver against the same complete free/support candidate.
    /// Refusal or late cancellation preserves accepted support motion, native
    /// forces, pressure endpoints and all work accounts together.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_supports_and_admission<R, E>(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        supports: AdhesiveP1SupportEndpoint3<'_>,
        step_s: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
        admit: impl FnOnce(&AdhesiveP1Step3, &TetTransientRun3<'g, 'm, 'p, 's>) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<AdhesiveP1Error3>,
    {
        let (candidate, report) = self
            .prepare(profile, other, Some(supports), step_s, gate, cx)
            .map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(gate, cx).map_err(E::from)?;
        self.publish(candidate, report);
        Ok(receiving)
    }

    fn publish(&mut self, candidate: TetTransientRun3<'g, 'm, 'p, 's>, report: AdhesiveP1Step3) {
        self.native = candidate;
        std::mem::swap(&mut self.cells, &mut self.patch.cells);
        if let (Some(accepted), Some(scratch)) = (&mut self.vector_cells, &mut self.patch.vector) {
            std::mem::swap(accepted, &mut scratch.cells);
        }
        std::mem::swap(&mut self.represented, &mut self.patch.represented);
        std::mem::swap(&mut self.background, &mut self.patch.background);
        self.current = report.transfer;
        self.last = Some(report);
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        supports: Option<AdhesiveP1SupportEndpoint3<'_>>,
        step_s: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, AdhesiveP1Step3), AdhesiveP1Error3> {
        let (load, transfer) = self.patch.assemble(profile, other, gate, cx)?;
        let mut candidate = self.native.clone();
        let native = match supports {
            Some(AdhesiveP1SupportEndpoint3::Displacement(values)) => {
                candidate.advance_with_supports(&load, values, step_s, cx)?
            }
            Some(AdhesiveP1SupportEndpoint3::Acceleration(values)) => {
                candidate.advance_with_support_accelerations(&load, values, step_s, cx)?
            }
            None => candidate.advance(&load, step_s, cx)?,
        };
        if supports.is_some() {
            // Use the original support owner's Newmark travel h*v_mid on ALL
            // components. Absolute u1-u0 can round away an acceleration-driven
            // increment about a nonzero support datum. Native work retains it.
            // Extract only velocity; the other temporary full fields drop here.
            let [_, before, _] = self.native.nodal_motion(cx)?;
            for (node, row) in before.iter().enumerate() {
                if node % 256 == 0 {
                    poll(gate, cx)?;
                }
                for c in 0..3 {
                    self.delta[3 * node + c] = 0.5 * row[c];
                }
            }
            drop(before);
            let [_, after, _] = candidate.nodal_motion(cx)?;
            for (node, row) in after.iter().enumerate() {
                if node % 256 == 0 {
                    poll(gate, cx)?;
                }
                for c in 0..3 {
                    self.delta[3 * node + c] =
                        finite(step_s * finite(self.delta[3 * node + c] + 0.5 * row[c])?)?;
                }
            }
        } else {
            self.delta.fill(0.0);
            for ((&d, &after), &before) in candidate
                .free_dofs()
                .iter()
                .zip(&candidate.state().displacement)
                .zip(&self.native.state().displacement)
            {
                poll(gate, cx)?;
                self.delta[d] = finite(after - before)?;
            }
        }
        let mut actual = Dd::ZERO;
        let mut other_work = Dd::ZERO;
        for (d, &change) in self.delta.iter().enumerate() {
            if d % 256 == 0 {
                poll(gate, cx)?;
            }
            actual =
                actual + (self.represented[d] + self.patch.represented[d]) * product(0.5, change);
            other_work = other_work
                + (Dd::from_f64(self.background[d]) + Dd::from_f64(self.patch.background[d]))
                    * product(0.5, change);
        }
        let mut spatial = Dd::ZERO;
        let mut scale = Dd::ZERO;
        let mut estimate = Dd::ZERO;
        for (cell_index, ((trace, before), after)) in self
            .patch
            .traces
            .iter()
            .zip(&self.cells)
            .zip(&self.patch.cells)
            .enumerate()
        {
            poll(gate, cx)?;
            let mut travel = [Dd::ZERO; 3];
            for (&node, coefficient) in trace.nodes.iter().zip(trace.coefficient) {
                let normal_travel = (0..3).fold(Dd::ZERO, |sum, axis| {
                    sum + product(self.patch.normal[axis], self.delta[3 * node + axis])
                });
                for i in 0..3 {
                    travel[i] = travel[i] + Dd::from_f64(coefficient[i]) * normal_travel;
                }
            }
            let before = [
                before.normal_force_n,
                before.first_moment_n_m[0],
                before.first_moment_n_m[1],
                before.moment_error_estimate_n_m,
            ];
            let after = [
                after.normal_force_n,
                after.first_moment_n_m[0],
                after.first_moment_n_m[1],
                after.moment_error_estimate_n_m,
            ];
            for i in 0..3 {
                let work = (Dd::from_f64(before[i]) + Dd::from_f64(after[i]))
                    * Dd::from_f64(0.5)
                    * travel[i];
                spatial = spatial + work;
                scale = scale + Dd::from_f64(finite(work.to_f64())?.abs());
            }
            let gradient = finite(travel[1].to_f64())?.hypot(finite(travel[2].to_f64())?);
            estimate = estimate
                + (Dd::from_f64(before[3]) + Dd::from_f64(after[3])) * product(0.5, gradient);
            if let (Some(accepted), Some(scratch)) = (&self.vector_cells, &self.patch.vector) {
                let before = &accepted[cell_index];
                let after = &scratch.cells[cell_index];
                for (component, direction) in [self.patch.frame.x_axis, self.patch.frame.y_axis]
                    .into_iter().enumerate()
                {
                    let mut travel = [Dd::ZERO; 3];
                    for (&node, coefficient) in trace.nodes.iter().zip(trace.coefficient) {
                        let tangential_travel = (0..3).fold(Dd::ZERO, |sum, axis| {
                            sum + product(direction[axis], self.delta[3 * node + axis])
                        });
                        for i in 0..3 {
                            travel[i] = travel[i] + Dd::from_f64(coefficient[i]) * tangential_travel;
                        }
                    }
                    let before = [before.tangential_force_n[component],
                        before.tangential_first_moment_n_m[component][0],
                        before.tangential_first_moment_n_m[component][1]];
                    let after = [after.tangential_force_n[component],
                        after.tangential_first_moment_n_m[component][0],
                        after.tangential_first_moment_n_m[component][1]];
                    for i in 0..3 {
                        let work = (Dd::from_f64(before[i]) + Dd::from_f64(after[i]))
                            * Dd::from_f64(0.5) * travel[i];
                        spatial = spatial + work;
                        scale = scale + Dd::from_f64(finite(work.to_f64())?.abs());
                    }
                }
            }
        }
        let config = self.patch.controls;
        let work_residual_j = finite((actual - spatial).to_f64())?;
        let work_tolerance_j = finite(
            config.work_absolute_tolerance_j
                + config.work_relative_tolerance
                    * finite(scale.to_f64())?.max(finite(actual.to_f64())?.abs()),
        )?;
        if work_residual_j.abs() > work_tolerance_j {
            return Err(AdhesiveP1Error3::Work {
                residual_j: work_residual_j,
                tolerance_j: work_tolerance_j,
            });
        }
        let quadrature_work_estimate_j = finite(estimate.to_f64())?;
        if quadrature_work_estimate_j > config.maximum_quadrature_work_estimate_j {
            return Err(AdhesiveP1Error3::QuadratureWork {
                estimate_j: quadrature_work_estimate_j,
                maximum_j: config.maximum_quadrature_work_estimate_j,
            });
        }
        let defect =
            finite((Dd::from_f64(native.mechanical_work_j) - actual - other_work).to_f64())?;
        let allowance = finite(
            self.native_config.absolute_energy_tolerance_j
                + self.native_config.relative_energy_tolerance
                    * native
                        .mechanical_work_j
                        .abs()
                        .max(finite(actual.to_f64())?.abs())
                        .max(finite(other_work.to_f64())?.abs()),
        )?;
        if defect.abs() > allowance {
            return Err(TransientError3::Energy {
                defect_j: defect,
                allowance_j: allowance,
            }
            .into());
        }
        let report = AdhesiveP1Step3 {
            native,
            transfer,
            spatial_traction_work_j: finite(spatial.to_f64())?,
            traction_work_on_solid_j: finite(actual.to_f64())?,
            other_external_work_j: finite(other_work.to_f64())?,
            work_residual_j,
            work_tolerance_j,
            quadrature_work_estimate_j,
        };
        poll(gate, cx)?;
        Ok((candidate, report))
    }
}
