use super::*;
use crate::phs_orbit::{PhsFloquetControls, PhsNumericalStability};
use fs_math::{c64::C64, det};
use fs_orbit::ShootBudget;
use fs_phs::{RelaxationBranch, Storage};

struct EuclideanStorage;
impl Storage for EuclideanStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * x.iter().map(|v| v * v).sum::<f64>()
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(x);
    }
}

fn oscillator(resistance: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![resistance, 0.0, 0.0, resistance],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(EuclideanStorage),
    )
    .unwrap()
}

// Authored external source: the complete closed ODE is a radial oscillator,
// r' = gain*(1-r^2)*r, theta'=-1. All original component resistance remains.
struct RadialDrive {
    gain: f64,
    resistance: f64,
    memory_stiffness: f64,
}
impl AutonomousPortFeedback for RadialDrive {
    fn port_count(&self) -> usize {
        2
    }
    fn input(&self, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        let radial = self.gain * (1.0 - x[0] * x[0] - x[1] * x[1]);
        out[0] = self.resistance * flows[0] + radial * x[0];
        out[1] = self.resistance * flows[1] + radial * x[1];
        if self.memory_stiffness != 0.0 {
            // Explicit actuator compensation, accounted as supplied work.
            // The memory state still evolves and dissipates its original heat.
            out[1] += self.memory_stiffness * x[0] - det::sqrt(self.memory_stiffness) * x[2];
        }
        true
    }
}

fn hb(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 60,
        tolerance: 1e-11,
        deflate_equilibrium: true,
    }
}
fn controls(n: usize) -> PhsAutonomousControls {
    PhsAutonomousControls {
        floquet: PhsFloquetControls {
            shooting: ShootBudget {
                steps_per_period: 256,
                max_newton: 30,
                tolerance: 1e-10,
            },
            state_tolerances: vec![1e-5; n],
            spectral_radius_tolerance: 2e-5,
            unit_circle_guard: 1e-2,
        },
        relative_period_tolerance: 1e-5,
        phase_multiplier_tolerance: 1e-3,
        audit_samples: 256,
    }
}
fn rotate(candidate: &HbOrbit, phase: f64) -> HbOrbit {
    let mut shifted = candidate.clone();
    for row in &mut shifted.coeffs {
        for (n, coefficient) in row.iter_mut().enumerate().skip(1) {
            let angle = n as f64 * phase;
            *coefficient = *coefficient * C64::new(det::cos(angle), det::sin(angle));
        }
    }
    shifted
}

#[test]
fn physical_feedback_finds_stable_and_unstable_cycles_with_correct_phase_sections() {
    let system = oscillator(0.2);
    for gain in [0.1, -0.1] {
        let feedback = RadialDrive {
            gain,
            resistance: 0.2,
            memory_stiffness: 0.0,
        };
        let candidate =
            solve_autonomous_feedback_harmonic_balance(&system, &feedback, 1.05, 0.8, &hb(3))
                .unwrap();
        assert!((candidate.omega - 1.0).abs() < 1e-8);
        assert!((candidate.first_harmonic_amplitude(0) - 1.0).abs() < 1e-8);
        for phase in [0.0, core::f64::consts::FRAC_PI_2] {
            let shifted = rotate(&candidate, phase);
            let result =
                analyze_autonomous_feedback_stability(&system, &feedback, &shifted, &controls(2))
                    .unwrap();
            assert_eq!(
                result.classification,
                if gain > 0.0 {
                    PhsNumericalStability::Stable
                } else {
                    PhsNumericalStability::Unstable
                }
            );
            assert_eq!(result.section_coordinate, if phase == 0.0 { 1 } else { 0 });
            assert_eq!(
                result.section_value.to_bits(),
                shifted.sample(result.section_coordinate, 0.0).to_bits()
            );
            assert!(
                (result.dominant_transverse_growth_rate_per_s.unwrap() + 2.0 * gain).abs() < 2e-5
            );
            for level in &result.levels {
                assert!((level.orbit.period - core::f64::consts::TAU).abs() < 1e-5);
                assert_eq!(level.orbit.multipliers.len(), 2);
                assert!(level.phase_multiplier_index.is_some());
                let expected = det::exp(-2.0 * gain * level.orbit.period);
                assert!((level.transverse_spectral_radius.unwrap() - expected).abs() < 2e-5);
                for i in 0..2 {
                    assert!((level.orbit.x0[i] - shifted.sample(i, 0.0)).abs() < 1e-5);
                }
            }
            assert!((result.candidate_audit.mean_supplied_power_w - 0.2).abs() < 1e-9);
            assert!((result.candidate_audit.mean_resistance_loss_w - 0.2).abs() < 1e-9);
            assert!(result.candidate_audit.relative_power_residual < 1e-9);
        }
    }
}

#[test]
fn adaptive_feedback_keeps_maxwell_memory_multiplier_and_original_source_heat() {
    let system = oscillator(0.2)
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 0.3,
            relaxation_time_s: 2.0,
        }])
        .unwrap();
    let feedback = RadialDrive {
        gain: 0.1,
        resistance: 0.2,
        memory_stiffness: 0.3,
    };
    let refinement = HbRefinementBudget {
        max_harmonics: 3,
        residual_tolerances: vec![1e-8; 3],
        state_tolerances: vec![1e-7; 3],
        frequency_tolerance: 1e-8,
    };
    let adaptive = solve_autonomous_feedback_harmonic_balance_adaptive(
        &system,
        &feedback,
        1.05,
        0.8,
        &hb(1),
        &refinement,
    )
    .unwrap();
    assert_eq!(adaptive.orbit.coeffs.len(), 3);
    assert!(adaptive.orbit.coeffs[0].len() > 2);
    let result =
        analyze_autonomous_feedback_stability(&system, &feedback, &adaptive.orbit, &controls(3))
            .unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Stable);
    for level in &result.levels {
        assert_eq!(level.orbit.multipliers.len(), 3);
        let phase = level.phase_multiplier_index.unwrap();
        let mut moduli: Vec<f64> = level
            .orbit
            .multipliers
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != phase)
            .map(|(_, &(re, im))| C64::new(re, im).abs())
            .collect();
        moduli.sort_by(|a, b| b.total_cmp(a));
        assert!((moduli[0] - det::exp(-0.2 * level.orbit.period)).abs() < 2e-5);
        assert!((moduli[1] - det::exp(-0.5 * level.orbit.period)).abs() < 2e-5);
        assert!((level.orbit.x0[2] - 0.2 * det::sqrt(0.3)).abs() < 1e-5);
    }
    // q=cos(t), p=-sin(t), z/sqrt(kb)=.2*cos(t)+.4*sin(t). Original R loss is
    // .2*mean([(1.24*cos-.12*sin)^2+sin^2]) + .3/2*mean((q-z)^2).
    let expected_loss = 0.1 * (1.24_f64.powi(2) + 0.12_f64.powi(2) + 1.0) + 0.06;
    assert!((result.candidate_audit.mean_stored_energy_j - 0.56).abs() < 1e-8);
    assert!((result.candidate_audit.mean_resistance_loss_w - expected_loss).abs() < 1e-8);
    assert!((result.candidate_audit.mean_supplied_power_w - expected_loss).abs() < 1e-8);
    assert!(result.candidate_audit.relative_power_residual < 1e-8);
}

struct VanDerPolDrive;
impl AutonomousPortFeedback for VanDerPolDrive {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        out[0] = (0.1 + 1.0 - x[0] * x[0]) * flows[0];
        true
    }
}
#[test]
fn self_excited_force_on_an_admitted_oscillator_matches_van_der_pol_period() {
    let system = fs_phs::mass_spring_damper(1.0, 1.0, 0.1).unwrap();
    let coarse =
        solve_autonomous_feedback_harmonic_balance(&system, &VanDerPolDrive, 0.94, 2.0, &hb(7))
            .unwrap();
    let candidate =
        solve_autonomous_feedback_harmonic_balance(&system, &VanDerPolDrive, 0.94, 2.0, &hb(17))
            .unwrap();
    let mut c = controls(2);
    c.floquet.shooting.steps_per_period = 512;
    c.floquet.state_tolerances = vec![5e-4; 2];
    c.relative_period_tolerance = 2e-4;
    // Accurate frequency alone does not resolve the physical waveform. The
    // N=7 seed has a small period error but fails the strict section-state check.
    assert!((core::f64::consts::TAU / coarse.omega - 6.663_286_859_323_130).abs() < 1e-6);
    assert!(matches!(
        analyze_autonomous_feedback_stability(&system, &VanDerPolDrive, &coarse, &c),
        Err(PhsOrbitError::SeedDisagreement { .. })
    ));
    let result =
        analyze_autonomous_feedback_stability(&system, &VanDerPolDrive, &candidate, &c).unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Stable);
    // Same independent published mu=1 period pinned in fs-orbit conformance.
    let period = 6.663_286_859_323_130;
    assert!((result.levels[2].orbit.period - period).abs() < 2e-7);
    assert!(result.levels[2].transverse_spectral_radius.unwrap() < 0.001);
    assert!(result.candidate_audit.mean_supplied_power_w > 0.0);
    assert!(result.candidate_audit.relative_power_residual < 1e-8);
    assert_eq!(system.structure().1[3], 0.1);
}

struct NoDrive;
impl AutonomousPortFeedback for NoDrive {
    fn port_count(&self) -> usize {
        0
    }
    fn input(&self, _: &[f64], _: &[f64], out: &mut [f64]) -> bool {
        out.is_empty()
    }
}
#[test]
fn conservative_one_to_one_family_keeps_all_neutral_directions_unresolved() {
    let mut j = vec![0.0; 16];
    for k in [0, 2] {
        j[k * 4 + k + 1] = 1.0;
        j[(k + 1) * 4 + k] = -1.0;
    }
    let system =
        PortHamiltonian::new(4, 0, j, vec![0.0; 16], vec![], Box::new(EuclideanStorage)).unwrap();
    let candidate = HbOrbit {
        omega: 1.0,
        coeffs: vec![
            vec![C64::ZERO, C64::from_re(0.5)],
            vec![C64::ZERO, C64::new(0.0, 0.5)],
            vec![C64::ZERO, C64::new(0.2, 0.1)],
            vec![C64::ZERO, C64::new(-0.1, 0.2)],
        ],
        residual_trace: vec![0.0],
        residual: 0.0,
        aft_samples: 16,
        sampling_defect: 0.0,
    };
    let mut c = controls(4);
    c.floquet.shooting.steps_per_period = 1024;
    let result = analyze_autonomous_feedback_stability(&system, &NoDrive, &candidate, &c).unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Unresolved);
    assert!(result.dominant_transverse_growth_rate_per_s.is_none());
    for level in &result.levels {
        assert_eq!(level.orbit.multipliers.len(), 4);
        assert!(level.phase_multiplier_index.is_none());
        assert!(level.transverse_spectral_radius.is_none());
    }
    assert_eq!(result.candidate_audit.mean_supplied_power_w, 0.0);
    assert_eq!(result.candidate_audit.mean_resistance_loss_w, 0.0);
}

#[test]
fn feedback_domain_incomplete_output_and_unresolved_section_are_real_refusals() {
    struct Drive(u8);
    impl AutonomousPortFeedback for Drive {
        fn port_count(&self) -> usize {
            2
        }
        fn input(&self, x: &[f64], _: &[f64], out: &mut [f64]) -> bool {
            match self.0 {
                0 => false,
                1 => {
                    out[0] = 0.0;
                    true
                }
                _ => {
                    out[0] = 0.2 * x[0] - x[1];
                    out[1] = x[0] + 0.2 * x[1];
                    true
                }
            }
        }
    }
    let system = oscillator(0.2);
    let candidate = HbOrbit {
        omega: 1.0,
        coeffs: vec![
            vec![C64::ZERO, C64::from_re(0.5)],
            vec![C64::ZERO, C64::new(0.0, 0.5)],
        ],
        residual_trace: vec![0.0],
        residual: 0.0,
        aft_samples: 16,
        sampling_defect: 0.0,
    };
    for mode in [0, 1] {
        assert_eq!(
            solve_autonomous_feedback_harmonic_balance(&system, &Drive(mode), 1.0, 1.0, &hb(1))
                .unwrap_err(),
            PhsOrbitError::SourceRefused
        );
        assert_eq!(
            analyze_autonomous_feedback_stability(&system, &Drive(mode), &candidate, &controls(2))
                .unwrap_err(),
            PhsOrbitError::SourceRefused
        );
        assert_eq!(
            audit_autonomous_feedback_orbit(&system, &Drive(mode), &candidate, 64).unwrap_err(),
            PhsOrbitError::SourceRefused
        );
    }
    assert!(matches!(
        analyze_autonomous_feedback_stability(&system, &Drive(2), &candidate, &controls(2)),
        Err(PhsOrbitError::Invalid { .. })
    ));
}

#[test]
fn autonomous_admission_precedes_feedback_sampling() {
    use std::cell::Cell;
    struct Counted {
        ports: usize,
        calls: Cell<usize>,
    }
    impl AutonomousPortFeedback for Counted {
        fn port_count(&self) -> usize {
            self.ports
        }
        fn input(&self, _: &[f64], _: &[f64], out: &mut [f64]) -> bool {
            self.calls.set(self.calls.get() + 1);
            out.fill(0.0);
            true
        }
    }
    let system = oscillator(0.2);
    let source = Counted {
        ports: 2,
        calls: Cell::new(0),
    };
    let candidate = HbOrbit {
        omega: 1.0,
        coeffs: vec![
            vec![C64::ZERO, C64::from_re(0.5)],
            vec![C64::ZERO, C64::new(0.0, 0.5)],
        ],
        residual_trace: vec![0.0],
        residual: 0.0,
        aft_samples: 16,
        sampling_defect: 0.0,
    };
    for which in 0..4 {
        let mut c = controls(2);
        match which {
            0 => c.floquet.shooting.steps_per_period = usize::MAX,
            1 => c.relative_period_tolerance = 0.5,
            2 => c.phase_multiplier_tolerance = c.floquet.unit_circle_guard,
            _ => c.floquet.state_tolerances = vec![1e-5],
        }
        assert!(analyze_autonomous_feedback_stability(&system, &source, &candidate, &c).is_err());
    }
    assert!(
        solve_autonomous_feedback_harmonic_balance(&system, &source, 1.0, 1.0, &hb(0)).is_err()
    );
    assert!(audit_autonomous_feedback_orbit(&system, &source, &candidate, 0).is_err());
    assert_eq!(source.calls.get(), 0);
    let wrong = Counted {
        ports: 1,
        calls: Cell::new(0),
    };
    assert!(solve_autonomous_feedback_harmonic_balance(&system, &wrong, 1.0, 1.0, &hb(1)).is_err());
    assert_eq!(wrong.calls.get(), 0);
}
