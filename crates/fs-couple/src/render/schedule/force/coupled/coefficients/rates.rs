//! Storage-preserving rate changes on the original accepted material history.
//!
//! A Maxwell arm stores (sqrt(k) b.q-z)^2/2, independently of tau. Changing
//! only tau changes mobility, not energy, projection, stress or internal strain.
//! The surrounding source owner still owns temperature and material validity.
use super::*;

/// Replacement relaxation time of an existing material coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialRelaxationTime {
    pub branch: usize,
    pub relaxation_time_s: f64,
}

/// Boundary observation at unchanged q, velocity and internal strain z.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialRelaxationTimeChange {
    pub branch: usize,
    pub previous_time_s: f64,
    pub relaxation_time_s: f64,
    pub stored_energy_j: f64,
    /// Instantaneous powers, NOT integrated heat over a time interval.
    pub power_before_w: f64,
    pub power_after_w: f64,
}

/// Only genuinely changed coordinates are returned. An unchanged batch is an
/// allocation-free no-op. Parameter work is exactly zero for this rate-only law.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MaterialRelaxationTimeReport {
    pub changes: Vec<MaterialRelaxationTimeChange>,
}

/// Native preparation retained until the enclosing physical transaction accepts.
pub(crate) struct PreparedMaterialRelaxationTimes {
    replacement: Option<(Vec<fs_phs::RelaxationBranch>, preparation::PreparedOperator)>,
    report: MaterialRelaxationTimeReport,
}

impl CoupledModalSystem {
    /// Change selected Maxwell relaxation times at one accepted boundary.
    ///
    /// Indices must be unique and strictly increasing; times and reciprocal
    /// rates must be positive finite. Stiffness, projections, every original
    /// mechanical/material/fluid/acoustic state, clock, source port and previous
    /// interval report stay unchanged. No history is re-relaxed or zeroed.
    ///
    /// Reuse the original ordinary/prepared/harmonic-only construction choice
    /// and its complete budgets. Rebuild the actual propagator/loss integrals
    /// BEFORE publication; any source, solver or cancellation failure leaves
    /// the original owner intact. An unchanged batch runs no preparation.
    ///
    /// Source-bound tensor points retain their own relaxation law and therefore
    /// refuse changes to their protected memories. This API supplies no thermal
    /// law, temperature interpolation, modulus change, caloric heat, or source
    /// substitution. Callers must hold the selected times over each interval,
    /// split prescribed events at boundaries and refine smooth changes in time.
    pub fn set_material_relaxation_times(
        &mut self,
        updates: &[MaterialRelaxationTime],
        gate: &CancelGate,
    ) -> Result<MaterialRelaxationTimeReport, ModalCouplingError> {
        let prepared = self.prepare_rates(updates, CoefficientEndpoint::Accepted, gate)?;
        poll(Some(gate))?;
        Ok(self.install_material_relaxation_times(prepared))
    }

    pub(crate) fn prepare_staged_material_relaxation_times(
        &self, updates: &[MaterialRelaxationTime], gate: &CancelGate,
    ) -> Result<PreparedMaterialRelaxationTimes, ModalCouplingError> {
        self.prepare_rates(updates, CoefficientEndpoint::Staged, gate)
    }

    fn prepare_rates(
        &self, updates: &[MaterialRelaxationTime], endpoint: CoefficientEndpoint, gate: &CancelGate,
    ) -> Result<PreparedMaterialRelaxationTimes, ModalCouplingError> {
        poll(Some(gate))?;
        if updates.len() > self.relaxation_branches.len() {
            return Err(invalid("relaxation time updates exceed original branch count"));
        }
        let mut previous = None;
        let mut changed = 0usize;
        for update in updates {
            poll(Some(gate))?;
            let original = self.relaxation_branches.get(update.branch)
                .ok_or_else(|| invalid("relaxation time update names an absent branch"))?;
            let rate = 1.0 / update.relaxation_time_s;
            if previous.is_some_and(|index| update.branch <= index)
                || !update.relaxation_time_s.is_finite() || update.relaxation_time_s <= 0.0
                || !rate.is_finite() || rate <= 0.0
            {
                return Err(invalid("ordered unique branches and representable positive relaxation times required"));
            }
            previous = Some(update.branch);
            if original.relaxation_time_s != update.relaxation_time_s {
                if self.plane_stress_material_points.iter().any(|point| point.memories.contains(&update.branch)) {
                    return Err(invalid("tensor material histories require their original source-law update owner"));
                }
                changed += 1;
            }
        }
        if changed == 0 {
            poll(Some(gate))?;
            return Ok(PreparedMaterialRelaxationTimes { replacement: None, report: MaterialRelaxationTimeReport::default() });
        }
        let (models, memory) = match endpoint {
            CoefficientEndpoint::Accepted => (&self.models, &self.material_memory),
            CoefficientEndpoint::Staged => (&self.candidates, &self.candidate_material_memory),
        };
        let mut branches = self.relaxation_branches.clone();
        let mut report = MaterialRelaxationTimeReport { changes: Vec::with_capacity(changed) };
        for update in updates {
            poll(Some(gate))?;
            let branch = &mut branches[update.branch];
            if branch.relaxation_time_s == update.relaxation_time_s { continue; }
            let mismatch = finite(super::super::relaxation::relaxed_memory(std::slice::from_ref(branch), models)?[0]
                - memory[update.branch])?;
            let stored_energy_j = finite(0.5 * mismatch * mismatch)?;
            let power_before_w = finite(mismatch * (mismatch / branch.relaxation_time_s))?;
            let power_after_w = finite(mismatch * (mismatch / update.relaxation_time_s))?;
            report.changes.push(MaterialRelaxationTimeChange { branch: update.branch,
                previous_time_s: branch.relaxation_time_s,
                relaxation_time_s: update.relaxation_time_s,
                stored_energy_j, power_before_w, power_after_w });
            branch.relaxation_time_s = update.relaxation_time_s;
        }
        // No stiffness or energy chart changes: the existing complete Nyquist
        // bound and all stored momenta/internal strains remain the same. The
        // original preparation still admits pole rates and solver resources.
        let operator = preparation::prepare(&self.models, &self.columns, &self.connections,
            self.mode_count(), self.dt, self.preparation(), &branches,
            self.radiation.as_ref(), gate)?;
        poll(Some(gate))?;
        Ok(PreparedMaterialRelaxationTimes { replacement: Some((branches, operator)), report })
    }

    // Publication is infallible and belongs immediately after joint admission.
    pub(crate) fn install_material_relaxation_times(
        &mut self, prepared: PreparedMaterialRelaxationTimes,
    ) -> MaterialRelaxationTimeReport {
        if let Some((branches, operator)) = prepared.replacement {
            self.exact = operator.exact;
            self.harmonic_damping_columns = operator.harmonic_damping_columns;
            self.relaxation_branches = branches;
        }
        prepared.report
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_phs::RelaxationBranch;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use fs_math::c64::C64;

    fn config() -> ModalCouplingConfig {
        ModalCouplingConfig { max_modes: 4, max_connections: 4, max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9 }
    }
    fn owner(tau: f64, q: f64, v: f64, memory: f64, coupled: bool,
        preparation: ModalPreparation) -> CoupledModalSystem {
        let gate = CancelGate::new_clock_free();
        let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(1000,
            vec![ModalAcousticMode { angular_frequency_rad_s: 20.0, damping_ratio: 0.01,
                pressure_per_modal_velocity: C64::ZERO }],
            ModalAcousticTimeBudget::audible_reference()).unwrap();
        model.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v }]).unwrap();
        CoupledModalSystem::new_with_preparation(vec![model], vec![], &[], config(), preparation, &gate)
            .unwrap().with_material_relaxation(vec![RelaxationBranch {
                projection: vec![if coupled { 1.0 } else { 0.0 }, 0.0],
                stiffness: 100.0, relaxation_time_s: tau }], Some(vec![memory]), &gate).unwrap()
    }
    fn near(a: f64, b: f64) {
        assert!((a-b).abs() < 2e-10*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}");
    }
    fn prepared() -> PreparedModalBudget {
        PreparedModalBudget { max_modes: 4, max_dense_bytes: 8*1024*1024,
            max_multiply_terms: 1_000_000_000 }
    }

    #[test]
    fn changing_rate_preserves_stored_history_and_matches_piecewise_analytic_decay() {
        let gate = CancelGate::new_clock_free();
        for preparation in [ModalPreparation::Ordinary, ModalPreparation::Prepared(prepared())] {
            let mut native = owner(0.04, 0.0, 0.0, 0.2, false, preparation);
            let initial = native.total_energy_j().unwrap();
            let mut heat = 0.0;
            for _ in 0..5 { heat += native.step_under_gate(&[0.0], &gate).unwrap().material_dissipation_j; }
            let before_memory = native.material_memory().to_vec();
            let before_states = native.components()[0].states().to_vec();
            let before_energy = native.total_energy_j().unwrap();
            let before_clock = native.samples_rendered();
            let report = native.set_material_relaxation_times(&[MaterialRelaxationTime {
                branch: 0, relaxation_time_s: 0.01 }], &gate).unwrap();
            assert_eq!(native.material_memory(), before_memory);
            assert_eq!(native.components()[0].states(), before_states);
            assert_eq!(native.total_energy_j().unwrap(), before_energy);
            assert_eq!(native.samples_rendered(), before_clock);
            near(report.changes[0].power_after_w, 4.0*report.changes[0].power_before_w);
            for tick in 1..=10 {
                let frame = native.step_under_gate(&[0.0], &gate).unwrap().clone();
                heat += frame.material_dissipation_j;
                near(native.material_memory()[0], 0.2*(-0.005_f64/0.04-f64::from(tick)*0.001/0.01).exp());
                near(native.total_energy_j().unwrap()+heat, initial);
            }
        }
    }

    #[test]
    fn warm_coupled_history_uses_new_propagator_and_same_native_harmonic_operator() {
        let gate = CancelGate::new_clock_free();
        for preparation in [ModalPreparation::Ordinary, ModalPreparation::Prepared(prepared()),
            ModalPreparation::HarmonicOnly(prepared())] {
            let mut native = owner(0.04, 0.01, -0.03, 0.007, true, preparation);
            if !preparation.is_harmonic_only() {
                for _ in 0..3 { native.step_under_gate(&[0.02], &gate).unwrap(); }
            }
            let current = native.components()[0].states()[0];
            let memory = native.material_memory()[0];
            let count = native.samples_rendered();
            let energy = native.total_energy_j().unwrap();
            native.set_material_relaxation_times(&[MaterialRelaxationTime { branch: 0,
                relaxation_time_s: 0.009 }], &gate).unwrap();
            let mut reference = owner(0.009, current.displacement_m_sqrt_kg,
                current.velocity_m_sqrt_kg_per_s, memory, true, preparation);
            assert_eq!(native.samples_rendered(), count);
            assert_eq!(native.total_energy_j().unwrap(), energy);
            for omega in [1.0, 17.0, 50.0] {
                let a = native.harmonic_response(omega, &[C64::new(0.1, -0.2)], 100_000, &gate).unwrap();
                let b = reference.harmonic_response(omega, &[C64::new(0.1, -0.2)], 100_000, &gate).unwrap();
                near(a.displacement_m_sqrt_kg[0].re, b.displacement_m_sqrt_kg[0].re);
                near(a.displacement_m_sqrt_kg[0].im, b.displacement_m_sqrt_kg[0].im);
                near(a.material_dissipation_w, b.material_dissipation_w);
            }
            if !preparation.is_harmonic_only() {
                for _ in 0..8 {
                    let a = native.step_under_gate(&[0.03], &gate).unwrap().clone();
                    let b = reference.step_under_gate(&[0.03], &gate).unwrap().clone();
                    near(a.total_energy_j(), b.total_energy_j());
                    near(native.material_memory()[0], reference.material_memory()[0]);
                    near(native.components()[0].states()[0].displacement_m_sqrt_kg,
                        reference.components()[0].states()[0].displacement_m_sqrt_kg);
                }
            }
        }
    }

    #[test]
    fn no_op_refusal_and_cancel_preserve_clock_state_rates_and_retry() {
        let gate = CancelGate::new_clock_free();
        let mut native = owner(0.04, 0.01, -0.03, 0.007, true, ModalPreparation::Ordinary);
        native.step_under_gate(&[0.02], &gate).unwrap();
        let memory = native.material_memory().to_vec();
        let state = native.components()[0].states().to_vec();
        let clock = native.samples_rendered();
        assert!(native.set_material_relaxation_times(&[MaterialRelaxationTime {
            branch: 0, relaxation_time_s: 0.04 }], &gate).unwrap().changes.is_empty());
        for tau in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1), f64::MIN_POSITIVE] {
            assert!(native.set_material_relaxation_times(&[MaterialRelaxationTime {
                branch: 0, relaxation_time_s: tau }], &gate).is_err());
        }
        assert!(native.set_material_relaxation_times(&[MaterialRelaxationTime {branch: 1,
            relaxation_time_s: 0.01}], &gate).is_err());
        assert!(native.set_material_relaxation_times(&[MaterialRelaxationTime {branch: 0,
            relaxation_time_s: 0.01}; 2], &gate).is_err());
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        assert!(native.set_material_relaxation_times(&[MaterialRelaxationTime {branch: 0,
            relaxation_time_s: 0.01}], &cancelled).is_err());
        assert_eq!(native.material_memory(), memory);
        assert_eq!(native.components()[0].states(), state);
        assert_eq!(native.samples_rendered(), clock);
        assert_eq!(native.material_relaxation_branches()[0].relaxation_time_s, 0.04);
        native.set_material_relaxation_times(&[MaterialRelaxationTime {branch: 0,
            relaxation_time_s: 0.01}], &gate).unwrap();
        native.step_under_gate(&[0.02], &gate).unwrap();
        assert_eq!(native.samples_rendered(), clock+1);
    }
}
