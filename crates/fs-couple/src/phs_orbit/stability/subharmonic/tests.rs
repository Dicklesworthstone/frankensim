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
