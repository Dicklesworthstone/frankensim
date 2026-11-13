use super::*;
use fs_phs::Storage;

fn hb(harmonics: usize) -> HbBudget {
    HbBudget { harmonics, max_newton: 60, tolerance: 1e-11, deflate_equilibrium: false }
}
fn seed(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega, coeffs, residual_trace: vec![], residual: f64::NAN,
        aft_samples: 0, sampling_defect: f64::NAN,
    }
}
fn rows(harmonics: usize) -> Vec<Vec<C64>> {
    vec![vec![C64::ZERO; harmonics + 1]; 2]
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:?} != {b:?}");
}

#[test]
fn period_embedding_preserves_every_physical_port_phase_and_dc() {
    let mut a = vec![C64::ZERO; MAX_INPUT_HARMONICS + 1];
    let mut b = a.clone();
    a[0] = C64::from_re(0.3);
    a[1] = C64::new(0.2, -0.7);
    a[7] = C64::new(-0.4, 0.1);
    b[0] = C64::from_re(-0.5);
    b[3] = C64::new(0.3, 0.6);
    let original = PeriodicPortInput::new(2.5, vec![a.clone(), b.clone()]).unwrap();
    let unchanged = original.over_periods(1).unwrap();
    assert_eq!(unchanged.coefficients, original.coefficients);
    for multiplier in [2, 3, 5] {
        let remapped = original.over_periods(multiplier).unwrap();
        assert_eq!(remapped.highest_nonzero_harmonic(), 7 * multiplier);
        assert_eq!(remapped.harmonic_count(), 7 * multiplier);
        for (old, new) in original.coefficients.iter().zip(&remapped.coefficients) {
            assert_eq!(new[0], old[0]);
            for n in 1..=7 {
                assert_eq!(new[n * multiplier], old[n]);
            }
            assert!(new.iter().enumerate().all(|(n, z)| n % multiplier == 0 || *z == C64::ZERO));
        }
        for time in [0.0, 0.013, 0.71, 2.3, 17.0] {
            let (mut before, mut after) = ([0.0; 2], [0.0; 2]);
            assert!(original.sample_into(time, &mut before));
            assert!(remapped.sample_into(time, &mut after));
            for i in 0..2 { assert!((before[i] - after[i]).abs() < 2e-12); }
        }
    }
    assert_eq!(original.coefficients, vec![a, b]);
}

#[test]
fn multiplier_caps_and_unrepresentable_periods_never_drop_authored_harmonics() {
    let mut row = vec![C64::ZERO; 130];
    row[129] = C64::from_re(f64::from_bits(1));
    let high = PeriodicPortInput::new(2.0, vec![row]).unwrap();
    assert!(high.over_periods(2).is_err());
    assert!(high.over_periods(0).is_err());
    assert!(high.over_periods(usize::MAX).is_err());
    let slow = PeriodicPortInput::new(f64::MIN_POSITIVE, vec![vec![C64::ONE]]).unwrap();
    assert!(slow.over_periods(256).is_err());
    let constant = PeriodicPortInput::new(2.0, vec![vec![C64::from_re(0.4)]]).unwrap();
    assert_eq!(constant.over_periods(256).unwrap().coefficients, constant.coefficients);
}

#[test]
fn complete_seed_solves_linear_phases_and_all_original_memory_coordinates() {
    let system = fs_phs::mass_spring_damper(2.0, 17.0, 0.3).unwrap()
        .with_relaxation_branches(vec![fs_phs::RelaxationBranch {
            projection: vec![1.0, 0.0], stiffness: 4.0, relaxation_time_s: 0.7,
        }]).unwrap();
    let force = C64::new(0.5, -0.2);
    let drive = PeriodicPortInput::new(2.0, vec![vec![C64::from_re(0.8), force]]).unwrap();
    let mut coefficients = vec![vec![C64::ZERO; 4]; 3];
    coefficients[0][1] = C64::new(0.08, 0.02);
    coefficients[1][1] = C64::new(-0.03, 0.1);
    coefficients[2][1] = C64::new(0.02, -0.01);
    let initial = seed(2.0, coefficients.clone());
    let orbit = drive.solve_harmonic_balance_seeded(&system, &initial, &hb(3)).unwrap();
    let iw = C64::new(0.0, 2.0);
    let modulus = (iw.scale(0.7) / (C64::ONE + iw.scale(0.7))).scale(4.0);
    let displacement = force / (C64::new(17.0 - 2.0 * 4.0, 0.3 * 2.0) + modulus);
    close(orbit.coeffs[0][0], C64::from_re(0.8 / 17.0), 1e-9);
    close(orbit.coeffs[0][1], displacement, 1e-9);
    close(orbit.coeffs[1][1], iw.scale(2.0) * displacement, 1e-9);
    assert!(orbit.coeffs[2][1].abs() > 1e-3);
    assert_eq!(initial.coeffs, coefficients);
    assert_eq!(orbit.omega.to_bits(), drive.omega_rad_s().to_bits());
}

struct Duffing;
impl Storage for Duffing {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * 0.625 * x[0] * x[0] + 0.125 * x[0].powi(4) + 0.5 * x[1] * x[1]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 0.625 * x[0] + 0.5 * x[0].powi(3);
        out[1] = x[1];
    }
}
fn duffing(damping: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2, 1, vec![0.0, 1.0, -1.0, 0.0], vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 1.0], Box::new(Duffing),
    ).unwrap()
}
fn third_harmonic_drive() -> PeriodicPortInput {
    PeriodicPortInput::new(3.0, vec![vec![C64::ZERO, C64::from_re(0.0625)]]).unwrap()
}

#[test]
fn exact_period_three_duffing_branch_keeps_drive_frequency_and_response_phase() {
    // q=cos(t), p=-sin(t) solves q''+.625q+.5q^3=.125cos(3t).
    // The original drive is 3 rad/s; its response is genuinely 1 rad/s.
    let drive = third_harmonic_drive();
    let chart = drive.over_periods(3).unwrap();
    let mut coefficients = rows(3);
    coefficients[0][1] = C64::from_re(0.5);
    coefficients[1][1] = C64::new(0.0, 0.5);
    let initial = seed(1.0, coefficients);
    let system = duffing(0.0);
    let orbit = chart.solve_harmonic_balance_seeded(&system, &initial, &hb(3)).unwrap();
    assert_eq!(drive.omega_rad_s(), 3.0);
    assert_eq!(orbit.omega, 1.0);
    close(orbit.coeffs[0][1], C64::from_re(0.5), 1e-10);
    close(orbit.coeffs[1][1], C64::new(0.0, 0.5), 1e-10);
    let refinement = HbRefinementBudget {
        max_harmonics: 7, residual_tolerances: vec![1e-8; 2],
        state_tolerances: vec![1e-8; 2], frequency_tolerance: 1e-10,
    };
    let refined = chart.solve_harmonic_balance_adaptive_seeded(
        &system, &initial, &hb(3), &refinement,
    ).unwrap();
    assert!(refined.refinement_trace.len() >= 2);
    assert_eq!(refined.orbit.coeffs[0].len(), 8);
    close(refined.orbit.coeffs[0][1], C64::from_re(0.5), 1e-9);
    assert!(refined.orbit.coeffs[0].iter().skip(2).all(|c| c.abs() < 1e-9));
}

#[test]
fn malformed_full_seeds_refuse_before_any_original_storage_evaluation() {
    use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};
    struct Counted(Arc<AtomicUsize>);
    impl Storage for Counted {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            self.0.fetch_add(1, Ordering::Relaxed);
            0.5 * x.iter().map(|v| v * v).sum::<f64>()
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            self.0.fetch_add(1, Ordering::Relaxed);
            out.copy_from_slice(x);
        }
    }
    let calls = Arc::new(AtomicUsize::new(0));
    let system = PortHamiltonian::new(
        2, 1, vec![0.0, 1.0, -1.0, 0.0], vec![0.0, 0.0, 0.0, 0.1],
        vec![0.0, 1.0], Box::new(Counted(calls.clone())),
    ).unwrap();
    calls.store(0, Ordering::Relaxed);
    let chart = third_harmonic_drive().over_periods(3).unwrap();
    for which in 0..7 {
        let mut initial = seed(1.0, rows(3));
        let mut options = hb(3);
        match which {
            0 => initial.omega = 3.0,
            1 => initial.coeffs[1].pop().map(|_| ()).unwrap(),
            2 => initial.coeffs[0][0].im = 1.0,
            3 => initial.coeffs[1][1].re = f64::NAN,
            4 => options.harmonics = usize::MAX,
            5 => options.harmonics = 4,
            _ => { initial.coeffs = rows(1); options.harmonics = 1; }
        }
        assert!(chart.solve_harmonic_balance_seeded(&system, &initial, &options).is_err());
    }
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn shorter_cycles_and_numerical_tails_cannot_masquerade_as_period_six() {
    let drive = third_harmonic_drive();
    let system = duffing(0.005);
    let mut waveform = seed(0.5, rows(6)); // six cycles of the 3 rad/s drive
    waveform.coeffs[0][3] = C64::from_re(0.5); // really a period-two response
    waveform.coeffs[1][3] = C64::new(0.0, 0.75);
    assert!(drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 2]).is_err());
    waveform.coeffs = rows(6);
    waveform.coeffs[0][2] = C64::from_re(0.5); // really a period-three response
    waveform.coeffs[1][2] = C64::new(0.0, 0.5);
    assert!(drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 2]).is_err());
    waveform.coeffs = rows(6);
    waveform.coeffs[0][6] = C64::from_re(0.5); // ordinary drive-periodic response
    waveform.coeffs[0][1] = C64::new(1e-14, -1e-14); // numerical tail, not a new branch
    assert!(drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 2]).is_err());
    waveform.coeffs[0][1] = C64::new(1e-4, -1e-4);
    waveform.coeffs[0][0] = C64::from_re(1e100); // no tolerance dilution by preload
    drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 2]).unwrap();
}

#[test]
fn original_memory_can_resolve_period_even_at_a_structural_observation_node() {
    let drive = third_harmonic_drive();
    let system = duffing(0.005).with_relaxation_branches(vec![fs_phs::RelaxationBranch {
        projection: vec![1.0, 0.0], stiffness: 0.3, relaxation_time_s: 2.0,
    }]).unwrap();
    // This test isolates waveform admission, not the equations of an orbit.
    let mut waveform = seed(0.5, vec![vec![C64::ZERO; 7]; 3]);
    waveform.coeffs[0][6] = C64::from_re(0.5);
    assert!(drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 3]).is_err());
    waveform.coeffs[2][1] = C64::new(0.01, -0.03);
    drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 3]).unwrap();
    // Changing phase cannot change the RMS repeat distance.
    let before = repeat_rms(&waveform.coeffs[2], 6, 3).unwrap();
    for (n, coefficient) in waveform.coeffs[2].iter_mut().enumerate().skip(1) {
        let phase = n as f64 * 0.73;
        *coefficient = *coefficient * C64::new(det::cos(phase), det::sin(phase));
    }
    let after = repeat_rms(&waveform.coeffs[2], 6, 3).unwrap();
    assert!((before - after).abs() < 1e-14);
    drive.require_subharmonic_separation(&system, 6, &waveform, &[1e-8; 3]).unwrap();
}

#[test]
fn parseval_repeat_distance_matches_direct_waveform_quadrature() {
    let mut coefficients = vec![C64::ZERO; 10];
    coefficients[0] = C64::from_re(0.3);
    for n in 1..10 { coefficients[n] = C64::new(0.2 / n as f64, -0.1 / (n * n) as f64); }
    let waveform = seed(1.0, vec![coefficients.clone()]);
    for periods in [2, 3, 4, 6, 12] {
        for repeat in 1..periods {
            if periods % repeat != 0 { continue; }
            let shift = core::f64::consts::TAU * repeat as f64 / periods as f64;
            let mean_square = (0..4096).map(|j| {
                let phase = core::f64::consts::TAU * j as f64 / 4096.0;
                (waveform.sample(0, phase + shift) - waveform.sample(0, phase)).powi(2)
            }).sum::<f64>() / 4096.0;
            assert!((repeat_rms(&coefficients, periods, repeat).unwrap() - det::sqrt(mean_square)).abs() < 1e-12);
        }
    }
    let mut commensurate = vec![C64::ZERO; 13];
    commensurate[6] = C64::new(1e100, -1e100);
    commensurate[12] = C64::new(-1e100, 1e100);
    assert_eq!(repeat_rms(&commensurate, 6, 1).unwrap(), 0.0);
}

#[test]
fn ambiguous_drive_clocks_and_wrong_state_tolerances_refuse() {
    let system = duffing(0.005);
    let mut waveform = seed(1.5, rows(4));
    waveform.coeffs[0][1] = C64::from_re(0.5);
    for coefficients in [
        vec![C64::ONE],
        vec![C64::ZERO, C64::ZERO, C64::ONE],
    ] {
        let drive = PeriodicPortInput::new(3.0, vec![coefficients]).unwrap();
        assert!(drive.require_subharmonic_separation(&system, 2, &waveform, &[1e-8; 2]).is_err());
    }
    let drive = third_harmonic_drive();
    assert!(drive.require_subharmonic_separation(&system, 2, &waveform, &[1e-8]).is_err());
    assert!(drive.require_subharmonic_separation(&system, 2, &waveform, &[0.0; 2]).is_err());
    assert!(drive.require_subharmonic_separation(&system, 2, &waveform, &[f64::NAN; 2]).is_err());
    waveform.omega = 3.0;
    assert!(drive.require_subharmonic_separation(&system, 2, &waveform, &[1e-8; 2]).is_err());
}

#[test]
fn damped_period_three_power_and_floquet_map_use_the_whole_response_cycle() {
    use crate::phs_orbit::PhsNumericalStability;
    use fs_orbit::ShootBudget;
    let drive = third_harmonic_drive();
    let chart = drive.over_periods(3).unwrap();
    let system = duffing(0.005);
    let mut coefficients = rows(13);
    coefficients[0][1] = C64::new(0.47, -0.16);
    coefficients[1][1] = C64::new(0.16, 0.47);
    coefficients[0][3] = C64::new(-0.0035, -0.0066);
    coefficients[1][3] = C64::new(0.0198, -0.0105);
    let initial = seed(1.0, coefficients);
    let orbit = chart.solve_harmonic_balance_seeded(&system, &initial, &hb(13)).unwrap();
    let tolerances = [1e-6; 2];
    drive.require_subharmonic_separation(&system, 3, &orbit, &tolerances).unwrap();
    // Independent high-accuracy ODE shooting of the explicit Duffing equation:
    // integrate T=2*pi from (q0,p0), solve endpoint-initial=0, then integrate
    // its analytical variational matrix [[0,1],[-.625-1.5*q*q,-.005]].
    assert!((orbit.sample(0, 0.0) - 0.938_957_088_939_154).abs() < 2e-7);
    assert!((orbit.sample(1, 0.0) - 0.352_519_600_812_236).abs() < 2e-7);
    assert!((orbit.sample(0, core::f64::consts::TAU / 3.0) - orbit.sample(0, 0.0)).abs() > 0.1);
    let a = drive.audit_subharmonic_orbit(&system, 3, &orbit, &tolerances, 512).unwrap();
    let b = drive.audit_subharmonic_orbit(&system, 3, &orbit, &tolerances, 1024).unwrap();
    assert!(a.relative_power_residual < 1e-7);
    assert!((a.mean_resistance_loss_w - 0.002_487_010_550_273).abs() < 2e-9);
    assert!((a.mean_supplied_power_w - b.mean_supplied_power_w).abs() < 1e-11);
    let controls = PhsFloquetControls {
        shooting: ShootBudget { steps_per_period: 512, max_newton: 30, tolerance: 1e-10 },
        state_tolerances: vec![1e-5; 2], spectral_radius_tolerance: 2e-5,
        unit_circle_guard: 1e-3,
    };
    let stability = drive.analyze_subharmonic_orbit_stability(
        &system, 3, &orbit, &tolerances, &controls,
    ).unwrap();
    assert_eq!(stability.classification, PhsNumericalStability::Stable);
    for level in &stability.levels {
        assert!((level.orbit.period - core::f64::consts::TAU).abs() < 1e-13);
        assert_eq!(level.orbit.multipliers.len(), 2); // forced; delete no phase direction
        // Complex-conjugate multipliers have modulus sqrt(exp(-damping*T)).
        assert!((level.spectral_radius - det::exp(-0.005 * level.orbit.period / 2.0)).abs() < 2e-5);
    }
    assert!((stability.dominant_growth_rate_per_s.unwrap() + 0.0025).abs() < 1e-5);
}
