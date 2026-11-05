//! Spatial basis refinement preserves physical motion, force and Airy storage.
use fs_couple::thin_plate::{VkBody, VkToneConfig};
use fs_nlmodal::SineMode;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_scenario::ThinPlate;
use std::f64::consts::PI;

fn plate(n_modes: usize) -> ThinPlate {
    ThinPlate {
        length_m: 0.6, width_m: 0.4, thickness_m: 0.001,
        density_kg_m3: 7800.0, e1_pa: 2.0e11, e2_pa: 2.0e11,
        nu12: 0.3, g12_pa: 2.0e11 / 2.6, material_angle_rad: 0.0,
        damping_ratio: 0.001, thermoelastic: None, kelvin_voigt_bending: None,
        n_modes, geometric_nonlinearity: true, pretension_n_m: 0.0, clamped: false,
    }
}

fn config(dt: f64) -> VkToneConfig {
    VkToneConfig {
        time_step_s: dt,
        nyquist_guard_fraction: 0.8,
        dense_budget: DenseExponentialBudget {
            max_dimension: 256,
            max_dense_bytes: 64 * 1024 * 1024,
            max_multiply_terms: 10_000_000_000,
        },
    }
}

#[test]
fn explicit_basis_retains_dark_motion_and_signed_physical_force_in_caller_order() {
    let p = plate(4);
    let displacement = [
        SineMode { m: 2, n: 1 }, // zero integrated area, nonzero point drive
        SineMode { m: 1, n: 3 }, // negative point drive at the physical bridge
        SineMode { m: 1, n: 1 },
        SineMode { m: 3, n: 1 },
    ];
    let stress = [SineMode { m: 1, n: 1 }, SineMode { m: 2, n: 1 }];
    let body = VkBody::from_sine_basis(p, &displacement, &stress, 32).unwrap();
    assert_eq!(body.state().len(), 8);
    let norm = (4.0 / (p.density_kg_m3 * p.thickness_m * p.length_m * p.width_m)).sqrt();
    assert_eq!(body.areas[0], 0.0);
    for (index, mode) in displacement.iter().enumerate().skip(1) {
        let expected = norm * 4.0 * p.length_m * p.width_m
            / (mode.m as f64 * mode.n as f64 * PI * PI);
        assert!((body.areas[index] - expected).abs() < 2e-15);
    }
    let dt = 1e-7;
    let force = 2.0;
    let mut tones = body.into_prepared_tones(&[0.0], config(dt), &[]).unwrap();
    let frame = tones.drive_and_radiate(&[force], &[0.0], 1.2, 1.0).unwrap();
    let traces = [1.0, -std::f64::consts::FRAC_1_SQRT_2,
        std::f64::consts::FRAC_1_SQRT_2, std::f64::consts::FRAC_1_SQRT_2];
    for (mode, trace) in traces.iter().enumerate() {
        // At rest, p_dot = physical modal force; the restoring correction is O(h²).
        let velocity = tones.state()[2 * mode + 1];
        assert!((velocity / dt - force * norm * trace).abs() < 3e-6);
    }
    assert!(tones.state()[0] > 0.0 && tones.state()[1] > 0.0);
    assert!(frame.work.supplied > 0.0);
    assert!(frame.work.balance_residual().abs() < 1e-13);
}

// Independent area quadrature of Psi * [Phi,Phi], with mass-normalized
// displacement and unit-L2 stress shapes. No production coupling helper is used.
fn independently_integrated_cubic_stiffness(p: ThinPlate, stress: &[SineMode]) -> f64 {
    let intervals = 128;
    let phi = (4.0 / (p.density_kg_m3 * p.thickness_m * p.length_m * p.width_m)).sqrt();
    let psi = (4.0 / (p.length_m * p.width_m)).sqrt();
    let ax = PI / p.length_m;
    let ay = PI / p.width_m;
    let weight = |i| if i == 0 || i == intervals { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
    stress.iter().map(|mode| {
        let mut integral = 0.0;
        for i in 0..=intervals {
            let x = p.length_m * i as f64 / intervals as f64;
            let (sx, cx) = (ax * x).sin_cos();
            for j in 0..=intervals {
                let y = p.width_m * j as f64 / intervals as f64;
                let (sy, cy) = (ay * y).sin_cos();
                let wxx = -phi * ax * ax * sx * sy;
                let wyy = -phi * ay * ay * sx * sy;
                let wxy = phi * ax * ay * cx * cy;
                let airy = psi * (mode.m as f64 * ax * x).sin()
                    * (mode.n as f64 * ay * y).sin();
                integral += weight(i) * weight(j) * airy * 2.0 * (wxx * wyy - wxy * wxy);
            }
        }
        integral *= p.length_m * p.width_m / (9.0 * (intervals * intervals) as f64);
        let xi2 = (mode.m as f64 * ax).powi(2) + (mode.n as f64 * ay).powi(2);
        p.e1_pa * p.thickness_m / (2.0 * xi2 * xi2) * integral * integral
    }).sum()
}

#[test]
fn independent_airy_refinement_reaches_original_quartic_storage_and_trajectory() {
    let p = plate(1);
    let displacement = [SineMode { m: 1, n: 1 }];
    let stress = [SineMode { m: 1, n: 1 }, SineMode { m: 3, n: 1 }, SineMode { m: 1, n: 3 }];
    let dt = 5e-5;
    let mut trajectories = Vec::new();
    for count in [1, 3] {
        let cubic = independently_integrated_cubic_stiffness(p, &stress[..count]);
        let mut tones = VkBody::from_sine_basis(p, &displacement, &stress[..count], count)
            .unwrap().into_prepared_tones(&[0.0], config(dt), &[]).unwrap();
        for _ in 0..200 {
            tones.drive_and_radiate(&[80.0], &[0.0], 1.2, 1.0).unwrap();
        }
        let (q, momentum) = (tones.state()[0], tones.state()[1]);
        let k2 = (PI / p.length_m).powi(2) + (PI / p.width_m).powi(2);
        let bending = p.e1_pa * p.thickness_m.powi(3) / (12.0 * (1.0 - p.nu12.powi(2)));
        let omega2 = (bending * k2 * k2 + p.pretension_n_m * k2)
            / (p.density_kg_m3 * p.thickness_m);
        let harmonic = 0.5 * (momentum * momentum + omega2 * q * q);
        let quartic = 0.25 * cubic * q.powi(4);
        assert!(quartic > 1e-4 * harmonic, "test must exercise nonlinear stretching");
        assert!((tones.stored_energy_j() - harmonic - quartic).abs() < 3e-6 * quartic);
        trajectories.push((q, cubic));
    }
    assert!(trajectories[1].1 > trajectories[0].1 * 1.01);
    assert!((trajectories[1].0 - trajectories[0].0).abs() > 1e-7);
}

#[test]
fn complete_spatial_basis_admission_refuses_budget_aliases_and_unsupported_physics() {
    let displacement = [SineMode { m: 1, n: 1 }, SineMode { m: 2, n: 1 }];
    let stress = [SineMode { m: 1, n: 1 }, SineMode { m: 3, n: 1 }];
    assert!(VkBody::from_sine_basis(plate(2), &displacement, &stress, 8).is_ok());
    assert!(VkBody::from_sine_basis(plate(2), &displacement, &stress, 7).is_err());
    assert!(VkBody::from_sine_basis(plate(1), &displacement, &stress, 8).is_err());
    assert!(VkBody::from_sine_basis(plate(0), &[], &stress, 8).is_err());
    assert!(VkBody::from_sine_basis(plate(2), &displacement, &[], 8).is_err());
    assert!(VkBody::from_sine_basis(plate(2), &[displacement[0]; 2], &stress, 8).is_err());
    assert!(VkBody::from_sine_basis(plate(2), &displacement, &[stress[0]; 2], 8).is_err());
    for field in [0, 1, 2] {
        let mut invalid = plate(2);
        match field { 0 => invalid.clamped = true, 1 => invalid.e2_pa *= 0.5, _ => invalid.thickness_m = 0.0 }
        assert!(VkBody::from_sine_basis(invalid, &displacement, &stress, 8).is_err());
    }
}
