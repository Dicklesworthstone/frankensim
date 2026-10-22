//! Branch continuation in the coupled set-valued friction owner (IRP pexfp.152).
//!
//! A light attachment on one mode is driven at constant speed into stick-slip.
//! Before each accepted step the test recomputes the held-force load line from
//! the network's own prediction: `v_free = (u dt - B^T dq_free)/dt` and
//! `m = c/dt`. An independent dense scan of `g(v) = v - v_free + m N mu(v)`,
//! plus the static cone, gives the intersections. The frame must disclose
//! that count, its committed traction must be the selected intersection, and
//! three-root steps must be taken stuck while sticking and on the high branch
//! while sliding.

use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_math::c64::C64;

fn system(law: FrictionLaw) -> FrictionModalSystem {
    let model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * 200.0, damping_ratio: 0.01,
        pressure_per_modal_velocity: C64::ZERO,
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    let config = ModalCouplingConfig { max_modes: 8, max_connections: 8, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 };
    let gate = CancelGate::new_clock_free();
    let network = CoupledModalSystem::new(vec![model], vec![], config, &gate).unwrap();
    // 1/sqrt(1e-4 kg): a light attachment, so the per-step mobility is large
    // enough for a three-root window at a 1 N load.
    FrictionModalSystem::new(network, ModalFrictionPort {
        left: ModalAttachment { component: 0, shapes: vec![100.0] }, right: None, law },
        ModalFrictionConfig { max_iterations: 96, maximum_force_n: 50.0, maximum_slip_speed_m_s: 100.0,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10, sticking_tolerance_m: 1e-12 },
        &gate).unwrap()
}

/// Sliding roots of g on (0, s] by a dense scan, each refined by bisection.
fn scan_roots(s: f64, m: f64, normal: f64, law: &FrictionLaw) -> Vec<f64> {
    let g = |v: f64| {
        let mu = if v == 0.0 { law.zero_speed_kinetic_coefficient().unwrap() } else { law.kinetic_coefficient(v).unwrap() };
        v - s + m * normal * mu
    };
    let cells = 50_000_u32;
    let mut roots = Vec::new();
    let mut prev = (0.0, g(0.0));
    for i in 1..=cells {
        let v = s * f64::from(i) / f64::from(cells);
        let cur = g(v);
        if (prev.1 < 0.0) != (cur < 0.0) {
            let (mut lo, mut hi, below) = (prev.0, v, prev.1 < 0.0);
            loop {
                let mid = lo + 0.5 * (hi - lo);
                if mid <= lo || mid >= hi { break; }
                if (g(mid) < 0.0) == below { lo = mid } else { hi = mid }
            }
            roots.push(hi);
        }
        prev = (v, cur);
    }
    roots
}

fn check_law(law: FrictionLaw) -> (usize, usize) {
    let mut system = system(law.clone());
    let drive = FrictionDrive { speed_m_s: 0.2, normal_force_n: 1.0 };
    let (dt, c) = (system.network.dt, system.compliance);
    let m = c / dt;
    let static_mu = law.static_coefficient().unwrap();
    let (mut stuck_three, mut high_three) = (0, 0);
    for _ in 0..6_000 {
        system.network.prepare_prediction(&[0.0], None).unwrap();
        let free = drive.speed_m_s * dt - dot(&system.column, &system.network.free_delta).unwrap();
        let v_free = free / dt;
        let frame = *system.step(&[0.0], drive).unwrap();
        let roots = frame.roots.expect("a loaded Stribeck/hyperbolic contact discloses roots");
        let cone = v_free.abs() / m <= static_mu * drive.normal_force_n;
        let sliding = scan_roots(v_free.abs(), m, drive.normal_force_n, &law);
        assert_eq!(usize::from(roots.found), usize::from(cone) + sliding.len(), "intersections at v_free {v_free}");
        assert!(roots.selected < roots.found);
        if frame.regime == ModalFrictionRegime::Sticking {
            assert!(cone && roots.selected == 0);
        } else {
            let root = sliding[usize::from(roots.selected) - usize::from(cone)];
            let expected = v_free.signum() * law.kinetic_coefficient(root).unwrap() * drive.normal_force_n;
            assert!((frame.traction_n - expected).abs() <= 1e-12 * expected.abs(), "traction {} vs {expected}", frame.traction_n);
        }
        if roots.found == 3 {
            if roots.selected == 0 { stuck_three += 1 } else if roots.selected == 2 { high_three += 1 } else {
                panic!("the unstable middle root was selected");
            }
        }
    }
    (stuck_three, high_three)
}

#[test]
fn g1_coupled_owner_follows_branch_continuation_and_discloses_roots() {
    for law in [
        FrictionLaw::Stribeck { static_mu: 1.2, kinetic_mu: 0.35, characteristic_speed: 0.05, viscous_per_speed: 0.0 },
        FrictionLaw::Hyperbolic { static_mu: 1.2, kinetic_mu: 0.35, characteristic_speed: 0.02 },
    ] {
        let (stuck, high) = check_law(law.clone());
        println!("{{\"test\":\"coupled_branch_continuation\",\"law\":\"{law:?}\",\"three_root_steps_stuck\":{stuck},\"three_root_steps_high\":{high}}}");
        assert!(stuck > 0 && high > 0, "both hysteresis branches must be exercised: {stuck} stuck, {high} high");
    }
}

#[test]
fn g0_other_laws_and_release_disclose_nothing() {
    let mut coulomb = system(FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 });
    let frame = *coulomb.step(&[0.0], FrictionDrive { speed_m_s: 0.2, normal_force_n: 1.0 }).unwrap();
    assert_eq!(frame.roots, None);
    let mut released = system(FrictionLaw::Stribeck { static_mu: 1.2, kinetic_mu: 0.35, characteristic_speed: 0.05, viscous_per_speed: 0.0 });
    let frame = *released.step(&[0.0], FrictionDrive { speed_m_s: 0.2, normal_force_n: 0.0 }).unwrap();
    assert_eq!((frame.regime, frame.roots), (ModalFrictionRegime::Released, None));
}
