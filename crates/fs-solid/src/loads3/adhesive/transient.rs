//! Original Newmark motion and independently reconstructed P1 pressure work.
use super::*;
use crate::linear3::{TetLinearElasticProblem, TetThermalStrainField};
use crate::loads3::transient::{
    TetTransientRun3, TransientConfig3, TransientError3, TransientStep3,
};

/// Subdivisions of the native mechanical-work account, never an extra debit.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveP1Step3 {
    pub native: TransientStep3,
    pub transfer: AdhesiveP1Transfer3,
    /// Average endpoint original pressure integrated against actual P1 travel.
    pub spatial_traction_work_j: f64,
    /// Average represented nodal adhesive correction dotted with that travel.
    pub traction_work_on_solid_j: f64,
    pub other_external_work_j: f64,
    pub work_residual_j: f64,
    pub work_tolerance_j: f64,
    pub quadrature_work_estimate_j: f64,
}

/// Fixed coplanar pressure patch on the existing complete sparse Newmark owner.
/// Matrices remain shared by the native candidate clone. Projection storage is
/// reused; the native trajectory keeps its existing bounded vector allocation.
/// Supports are stationary and zero. Read-only native access prevents an
/// untracked load/support change; consume this wrapper to continue elsewhere.
pub struct AdhesiveP1Transient3<'g, 'm, 'p, 's> {
    patch: AdhesiveP1Patch3<'g, 'm>,
    native: TetTransientRun3<'g, 'm, 'p, 's>,
    cells: Vec<SectorAdhesiveTraction>,
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
    /// Nonzero prescribed support motion has its own native path and is not
    /// silently omitted by this fixed-reference pressure adapter.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient<'p, 's>(
        mut self,
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
        let (load, current) = self.assemble(profile, other, gate, cx)?;
        let native = load.begin_transient(
            problem,
            thermal,
            time_s,
            displacement_m,
            velocity_m_s,
            config,
            cx,
        )?;
        let mut cells = vec![SectorAdhesiveTraction::default(); self.cells.len()];
        let mut represented = vec![Dd::ZERO; self.represented.len()];
        let mut background = vec![0.0; self.background.len()];
        let delta = vec![0.0; self.background.len()];
        std::mem::swap(&mut self.cells, &mut cells);
        std::mem::swap(&mut self.represented, &mut represented);
        std::mem::swap(&mut self.background, &mut background);
        poll(gate, cx)?;
        Ok(AdhesiveP1Transient3 {
            patch: self,
            native,
            cells,
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
        let (candidate, report) = self.prepare(profile, other, step_s, gate, cx)?;
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
            .prepare(profile, other, step_s, gate, cx)
            .map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(gate, cx).map_err(E::from)?;
        self.publish(candidate, report);
        Ok(receiving)
    }

    fn publish(&mut self, candidate: TetTransientRun3<'g, 'm, 'p, 's>, report: AdhesiveP1Step3) {
        self.native = candidate;
        std::mem::swap(&mut self.cells, &mut self.patch.cells);
        std::mem::swap(&mut self.represented, &mut self.patch.represented);
        std::mem::swap(&mut self.background, &mut self.patch.background);
        self.current = report.transfer;
        self.last = Some(report);
    }

    fn prepare(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        step_s: f64,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, AdhesiveP1Step3), AdhesiveP1Error3> {
        let (load, transfer) = self.patch.assemble(profile, other, gate, cx)?;
        let mut candidate = self.native.clone();
        let native = candidate.advance(&load, step_s, cx)?;
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
        // Fixed supports have exactly zero actual travel; their forces remain
        // in the original native load and later reaction recovery, not dropped.
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
        for ((trace, before), after) in self
            .patch
            .traces
            .iter()
            .zip(&self.cells)
            .zip(&self.patch.cells)
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
