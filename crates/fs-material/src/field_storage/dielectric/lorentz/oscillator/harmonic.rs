//! Continuous electric-field tones through the shared exact oscillator owner.

use super::{LorentzError, LorentzPole, LorentzState, LorentzStep, Result, finite};
use fs_modal::time::{ModalAcousticState, PreparedHarmonicForce};

/// Immutable interval preparation, retaining physical polarization storage and
/// independently integrated supplied work and frictional heat.
#[derive(Clone, Debug)]
pub struct PreparedLorentzHarmonic {
    pole: LorentzPole,
    transition: PreparedHarmonicForce,
}

impl LorentzPole {
    /// Prepare `E(t)=Ec*cos(Omega*t)+Es*sin(Omega*t)` on one positive interval.
    /// Frequencies are rad/s, damping remains the original physical rate, and
    /// the field coefficients are local to the initial state's clock. Lossless
    /// exact resonance retains secular growth; no steady-state inverse is used.
    /// Zero driving frequency is the static cosine field. Observation and field
    /// bandwidth remain the caller's responsibility; no Nyquist cutoff is added.
    pub fn prepare_harmonic(
        self,
        drive_angular_frequency_rad_s: f64,
        duration_s: f64,
    ) -> Result<PreparedLorentzHarmonic> {
        Ok(PreparedLorentzHarmonic {
            pole: self,
            transition: PreparedHarmonicForce::try_new_with_damping_rate(
                self.omega,
                self.damping,
                drive_angular_frequency_rad_s,
                duration_s,
            )?,
        })
    }
}

impl PreparedLorentzHarmonic {
    /// Advance the original P/J vectors under the continuous field, without
    /// allocation or transcendental calls. For the parent's `exp(+i Omega t)`
    /// susceptibility, these same coefficients correspond to `Ec - i*Es`.
    /// No accepted state or phase is retained; every failed trial is immutable.
    pub fn advance<const N: usize>(
        &self,
        before: LorentzState<N>,
        cosine_field_v_m: [f64; N],
        sine_field_v_m: [f64; N],
    ) -> Result<LorentzStep<N>> {
        let stored_before = self.pole.stored_energy_j_m3(before)?;
        let mut state = LorentzState::default();
        let mut work = 0.0;
        let mut heat = 0.0;
        let mut tolerance = 0.0;
        for i in 0..N {
            let motion = self.transition.advance_with_energy(
                ModalAcousticState {
                    displacement_m_sqrt_kg: finite(
                        before.polarization_c_m2[i] / self.pole.coupling,
                    )?,
                    velocity_m_sqrt_kg_per_s: finite(before.current_a_m2[i] / self.pole.coupling)?,
                },
                finite(self.pole.coupling * cosine_field_v_m[i])?,
                finite(self.pole.coupling * sine_field_v_m[i])?,
            )?;
            state.polarization_c_m2[i] =
                finite(self.pole.coupling * motion.state.displacement_m_sqrt_kg)?;
            state.current_a_m2[i] =
                finite(self.pole.coupling * motion.state.velocity_m_sqrt_kg_per_s)?;
            work = finite(work + motion.input_work_j)?;
            heat = finite(heat + motion.viscous_dissipation_j)?;
            tolerance = finite(tolerance + motion.dissipation_roundoff_tolerance_j)?;
        }
        let stored_after = self.pole.stored_energy_j_m3(state)?;
        Ok(LorentzStep {
            state,
            stored_before_j_m3: stored_before,
            stored_after_j_m3: stored_after,
            supplied_work_j_m3: work,
            heat_j_m3: heat,
            heat_roundoff_tolerance_j_m3: tolerance,
            energy_residual_j_m3: finite((stored_after - stored_before) + (heat - work))?,
        })
    }

    /// Rotate the original field coefficients to the next interval's clock.
    /// Apply only after accepting the associated motion. Repeated rotations
    /// accumulate floating-point phase error; no amplitude renormalization,
    /// material change, or phase reset is performed.
    pub fn advance_field_phase<const N: usize>(
        &self,
        cosine_field_v_m: [f64; N],
        sine_field_v_m: [f64; N],
    ) -> Result<([f64; N], [f64; N])> {
        if !(1..=3).contains(&N) {
            return Err(LorentzError::Invalid(
                "one to three electric field components required",
            ));
        }
        let mut cosine = [0.0; N];
        let mut sine = [0.0; N];
        for i in 0..N {
            (cosine[i], sine[i]) = self
                .transition
                .motion()
                .advance_force_phase(cosine_field_v_m[i], sine_field_v_m[i])?;
        }
        Ok((cosine, sine))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::c64::C64;

    fn near(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 8e-11 * scale.max(1e-300),
            "actual {actual:.17e}, expected {expected:.17e}, scale {scale:.17e}"
        );
    }

    #[test]
    fn g1_undamped_resonant_field_has_analytic_secular_growth_and_work() {
        let (delta, omega, duration) = (2e-6, 10.0, 2.37);
        let pole = LorentzPole::try_new(delta, omega, 0.0).unwrap();
        let cosine = [400.0, -70.0, 0.0];
        let sine = [-90.0, 30.0, 180.0];
        let frame = pole
            .prepare_harmonic(omega, duration)
            .unwrap()
            .advance(LorentzState::default(), cosine, sine)
            .unwrap();
        let angle = omega * duration;
        let (sn, cs) = angle.sin_cos();
        for i in 0..3 {
            let p = 0.5 * delta * (cosine[i] * angle * sn + sine[i] * (sn - angle * cs));
            let j = 0.5 * delta * omega * (cosine[i] * (sn + angle * cs) + sine[i] * angle * sn);
            near(frame.state.polarization_c_m2[i], p, 0.02);
            near(frame.state.current_a_m2[i], j, 0.2);
        }
        assert_eq!(frame.heat_j_m3, 0.0);
        near(
            frame.supplied_work_j_m3,
            frame.stored_after_j_m3,
            frame.stored_after_j_m3,
        );
        near(frame.energy_residual_j_m3, 0.0, frame.stored_after_j_m3);
        assert!(matches!(
            pole.susceptibility_f_m(omega),
            Err(LorentzError::SingularResonance)
        ));
    }

    #[test]
    fn g1_periodic_phasor_state_and_cycle_heat_match_the_original_lorentz_law() {
        let (delta, omega): (f64, f64) = (4e-6, 10.0);
        let cosine = [400.0, -80.0, 30.0];
        let sine = [-60.0, 90.0, 100.0];
        for gamma in [0.4, 20.0, 35.0] {
            let pole = LorentzPole::try_new(delta, omega, gamma).unwrap();
            for drive in [5.0, 10.0, 17.0] {
                let duration = core::f64::consts::TAU / drive;
                let susceptibility = C64::from_re(delta * omega * omega)
                    / C64::new(omega * omega - drive * drive, gamma * drive);
                let mut initial = LorentzState::default();
                let mut expected_heat = 0.0;
                for i in 0..3 {
                    let p = susceptibility * C64::new(cosine[i], -sine[i]);
                    initial.polarization_c_m2[i] = p.re;
                    initial.current_a_m2[i] = -drive * p.im;
                    expected_heat += 0.5 * duration * gamma * drive * drive * p.norm_sq()
                        / (delta * omega * omega);
                }
                let frame = pole
                    .prepare_harmonic(drive, duration)
                    .unwrap()
                    .advance(initial, cosine, sine)
                    .unwrap();
                for i in 0..3 {
                    near(
                        frame.state.polarization_c_m2[i],
                        initial.polarization_c_m2[i],
                        0.1,
                    );
                    near(frame.state.current_a_m2[i], initial.current_a_m2[i], 1.0);
                }
                assert!(frame.heat_j_m3 > 0.0);
                near(frame.heat_j_m3, expected_heat, expected_heat);
                near(frame.supplied_work_j_m3, expected_heat, expected_heat);
                near(
                    frame.energy_residual_j_m3,
                    0.0,
                    frame.stored_before_j_m3 + expected_heat,
                );
            }
        }
    }

    #[test]
    fn g3_phase_correct_subdivision_preserves_vector_state_work_and_heat() {
        for gamma in [0.0, 0.2, 20.0, 35.0] {
            let pole = LorentzPole::try_new(4e-6, 10.0, gamma).unwrap();
            let initial = LorentzState {
                polarization_c_m2: [1e-3, -2e-4, 7e-4],
                current_a_m2: [-2e-3, 3e-3, 1e-3],
            };
            let (cosine, sine) = ([400.0, -20.0, 80.0], [-20.0, 800.0, -500.0]);
            let whole = pole
                .prepare_harmonic(7.3, 0.24)
                .unwrap()
                .advance(initial, cosine, sine)
                .unwrap();
            let part = pole.prepare_harmonic(7.3, 0.08).unwrap();
            let (mut c, mut s) = (cosine, sine);
            let mut state = initial;
            let (mut work, mut heat) = (0.0, 0.0);
            for _ in 0..3 {
                let frame = part.advance(state, c, s).unwrap();
                state = frame.state;
                work += frame.supplied_work_j_m3;
                heat += frame.heat_j_m3;
                (c, s) = part.advance_field_phase(c, s).unwrap();
            }
            for i in 0..3 {
                near(
                    whole.state.polarization_c_m2[i],
                    state.polarization_c_m2[i],
                    0.01,
                );
                near(whole.state.current_a_m2[i], state.current_a_m2[i], 0.1);
            }
            near(whole.supplied_work_j_m3, work, 10.0);
            near(whole.heat_j_m3, heat, 10.0);
            near(whole.energy_residual_j_m3, 0.0, 10.0);
            assert_eq!(
                whole,
                pole.prepare_harmonic(7.3, 0.24)
                    .unwrap()
                    .advance(initial, cosine, sine)
                    .unwrap()
            );
        }
    }

    #[test]
    fn g0_static_limit_and_transactional_input_refusals() {
        let pole = LorentzPole::try_new(2e-6, 10.0, 2.0).unwrap();
        let initial = LorentzState::<2>::default();
        let static_frame = pole
            .prepare_harmonic(0.0, 0.1)
            .unwrap()
            .advance(initial, [300.0, -100.0], [50.0, 70.0])
            .unwrap();
        let held = pole
            .prepare_step(0.1)
            .unwrap()
            .advance(initial, [300.0, -100.0], [300.0, -100.0])
            .unwrap();
        for i in 0..2 {
            near(
                static_frame.state.polarization_c_m2[i],
                held.state.polarization_c_m2[i],
                0.001,
            );
            near(
                static_frame.state.current_a_m2[i],
                held.state.current_a_m2[i],
                0.01,
            );
        }
        near(static_frame.heat_j_m3, held.heat_j_m3, 1.0);
        near(
            static_frame.supplied_work_j_m3,
            held.supplied_work_j_m3,
            1.0,
        );
        let plan = pole.prepare_harmonic(7.0, 0.1).unwrap();
        assert!(plan.advance(initial, [1.0, f64::NAN], [0.0; 2]).is_err());
        assert!(
            plan.advance_field_phase([0.0; 2], [0.0, f64::INFINITY])
                .is_err()
        );
        assert!(plan.advance_field_phase::<0>([], []).is_err());
        assert!(plan.advance(LorentzState::<0>::default(), [], []).is_err());
        assert!(pole.prepare_harmonic(-1.0, 0.1).is_err());
        assert!(pole.prepare_harmonic(7.0, 0.0).is_err());
        assert_eq!(initial, LorentzState::default());
    }
}
