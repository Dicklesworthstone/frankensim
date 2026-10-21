//! Endpoint acceleration readout of the COMPLETE staged mechanical solve.
//! No differentiation of sampled velocities or narrow-band phase substitution.

use super::*;

impl FrictionModalSystem {
    /// Observe physical endpoint acceleration rather than a narrow-band transfer.
    ///
    /// `weights[i]` maps the original mass-normalized modal acceleration to
    /// pressure [Pa s² / (m sqrt(kg))], in component-then-mode order. A compact
    /// baffled radiator supplies rho*area/(2*pi*r*sqrt(modal_mass)). Signed areas
    /// and exact zeros remain meaningful. The producer owns units, geometry,
    /// source validity and the radiation approximation; this API invents none.
    ///
    /// Acceleration uses the SAME held forces that advanced the candidate,
    /// including friction, material-memory and all bilateral reactions, minus its elastic and
    /// viscous restoring loads. Thus a static force-balanced deflection emits
    /// no pressure. No radiation reaction, delay, directivity or pressure-field
    /// energy is inferred from these observation coefficients.
    ///
    /// All original component pressure transfers must be zero, preventing
    /// double radiation and preserving their original pressure budgets. Bind
    /// before the first sample; live observer/medium changes are not supported.
    /// Nonfinite or over-budget pressure refuses before publishing ANY state,
    /// frame or clock. The original network pressure ceiling remains in force.
    pub fn with_acceleration_pressure(mut self, weights: Vec<f64>) -> Result<Self, ModalFrictionError> {
        if self.samples_rendered() != 0 || self.acceleration_pressure.is_some()
            || weights.len() != self.mode_count() || weights.iter().any(|x| !x.is_finite())
            || self.components().iter().flat_map(|m| m.modes()).any(|m|
                m.pressure_per_modal_velocity != fs_math::c64::C64::ZERO) {
            return Err(invalid("acceleration pressure requires an unstepped network, one finite coefficient per mode and no existing pressure transfers").into());
        }
        self.acceleration_pressure = Some(weights);
        Ok(self)
    }

    /// Exact admitted observation coefficients, without source-authority changes.
    #[must_use]
    pub fn acceleration_pressure_weights(&self) -> Option<&[f64]> {
        self.acceleration_pressure.as_deref()
    }

    pub(super) fn staged_observer_pressure(&self, gate: Option<&CancelGate>) -> Result<f64, ModalFrictionError> {
        let Some(weights) = &self.acceleration_pressure else {
            return Ok(self.network.staged_frame.observer_pressure_pa);
        };
        let mut pressure = 0.0;
        let mut index = 0;
        for component in &self.network.candidates {
            poll(gate)?;
            for (mode, state) in component.modes().iter().zip(component.states()) {
                let weight = weights[index];
                if weight != 0.0 {
                    let acceleration = finite(self.network.forces[index]
                        - 2.0 * mode.damping_ratio * mode.angular_frequency_rad_s
                            * state.velocity_m_sqrt_kg_per_s
                        - mode.angular_frequency_rad_s * mode.angular_frequency_rad_s
                            * state.displacement_m_sqrt_kg)?;
                    pressure = finite(pressure + weight * acceleration)?;
                }
                index += 1;
            }
        }
        limit("acceleration observer pressure", pressure.abs(), self.network.config.maximum_abs_pressure_pa)?;
        Ok(pressure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget};
    use crate::render::schedule::force::coupled::{ModalConnection, ModalCouplingConfig};
    use fs_math::c64::C64;

    fn system(displacement: f64, pressure_limit: f64) -> FrictionModalSystem {
        let gate = CancelGate::new_clock_free();
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
            angular_frequency_rad_s: 100.0, damping_ratio: 0.1,
            pressure_per_modal_velocity: C64::ZERO,
        }], ModalAcousticTimeBudget::audible_reference()).unwrap();
        model.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: displacement, velocity_m_sqrt_kg_per_s: 0.0,
        }]).unwrap();
        let attachment = ModalAttachment { component: 0, shapes: vec![1.0] };
        let network = CoupledModalSystem::new(vec![model], vec![ModalConnection {
            left: attachment.clone(), right: ModalAttachment { component: 0, shapes: vec![0.0] },
            stiffness_n_m: 100.0, damping_n_s_m: 0.1, rest_extension_m: 0.0,
        }], ModalCouplingConfig {
            max_modes: 4, max_connections: 2, max_setup_terms: 128,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
            maximum_abs_pressure_pa: pressure_limit, maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1.0e-10, energy_absolute_tolerance_j: 1.0e-12,
            energy_relative_tolerance: 1.0e-8,
        }, &gate).unwrap();
        FrictionModalSystem::new(network, ModalFrictionPort {
            left: attachment, right: None,
            law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 },
        }, ModalFrictionConfig {
            max_iterations: 96, maximum_force_n: 100.0, maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1.0e-10, force_relative_tolerance: 1.0e-9,
            sticking_tolerance_m: 1.0e-12,
        }, &gate).unwrap()
    }
    fn released() -> FrictionDrive { FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 } }

    #[test]
    fn static_deflection_does_not_radiate_and_force_feedthrough_is_present() {
        let mut equilibrium = system(0.001, 1000.0).with_acceleration_pressure(vec![2.0]).unwrap();
        // Spring-to-ground adds 0.1 N to the oscillator's 10 N restoring load.
        for _ in 0..32 {
            let p = equilibrium.step(&[10.1], released()).unwrap().observer_pressure_pa;
            assert!(p.abs() < 1.0e-10, "static body radiates {p}");
        }
        let mut driven = system(0.0, 1000.0).with_acceleration_pressure(vec![2.0]).unwrap();
        let frame = *driven.step(&[1.0], released()).unwrap();
        let state = driven.components()[0].states()[0];
        let reaction = driven.network.frame.connection_forces_n[0];
        let expected = 2.0 * (1.0 + reaction - 20.0 * state.velocity_m_sqrt_kg_per_s
            - 10_000.0 * state.displacement_m_sqrt_kg);
        assert!((frame.observer_pressure_pa - expected).abs() < 1.0e-14);
        assert!(frame.observer_pressure_pa > 1.9, "force feedthrough disappeared");
    }

    #[test]
    fn friction_and_connection_reactions_are_in_the_same_observed_acceleration() {
        let mut model = system(0.0, 1000.0).with_acceleration_pressure(vec![-0.25]).unwrap();
        for _ in 0..64 {
            let frame = *model.step(&[0.0], FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 }).unwrap();
            let s = model.components()[0].states()[0];
            let acceleration = frame.traction_n + model.network.frame.connection_forces_n[0]
                - 20.0 * s.velocity_m_sqrt_kg_per_s - 10_000.0 * s.displacement_m_sqrt_kg;
            assert!((frame.observer_pressure_pa + 0.25 * acceleration).abs() < 1.0e-14);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    }

    #[test]
    fn pressure_refusal_retains_previous_state_frame_and_clock() {
        let mut model = system(0.0, 0.01).with_acceleration_pressure(vec![2.0]).unwrap();
        model.step(&[0.0], released()).unwrap();
        let before = model.components()[0].states().to_vec();
        let frame = *model.last_frame().unwrap();
        assert!(model.step(&[1.0], released()).is_err());
        assert_eq!(model.components()[0].states(), before);
        assert_eq!(model.last_frame(), Some(&frame));
        assert_eq!(model.samples_rendered(), 1);
        assert_eq!(model.step(&[0.0], released()).unwrap().observer_pressure_pa, 0.0);
        assert_eq!(model.samples_rendered(), 2);
    }

    #[test]
    fn signed_observation_changes_only_pressure_not_mechanics_or_work() {
        let mut a = system(0.0, 1000.0).with_acceleration_pressure(vec![1.0]).unwrap();
        let mut b = system(0.0, 1000.0).with_acceleration_pressure(vec![-2.0]).unwrap();
        for _ in 0..32 {
            let fa = *a.step(&[0.5], released()).unwrap();
            let fb = *b.step(&[0.5], released()).unwrap();
            assert_eq!(fb.observer_pressure_pa, -2.0 * fa.observer_pressure_pa);
            assert_eq!(a.components()[0].states(), b.components()[0].states());
            assert_eq!(fa.network_energy_j, fb.network_energy_j);
            assert_eq!(fa.external_work_j, fb.external_work_j);
        }
    }

    #[test]
    fn malformed_duplicate_and_live_observers_refuse() {
        for weights in [vec![], vec![1.0, 2.0], vec![f64::NAN], vec![f64::INFINITY]] {
            assert!(system(0.0, 1000.0).with_acceleration_pressure(weights).is_err());
        }
        let bound = system(0.0, 1000.0).with_acceleration_pressure(vec![1.0]).unwrap();
        assert!(bound.with_acceleration_pressure(vec![1.0]).is_err());
        let mut live = system(0.0, 1000.0);
        live.step(&[0.0], released()).unwrap();
        assert!(live.with_acceleration_pressure(vec![1.0]).is_err());
    }
}
