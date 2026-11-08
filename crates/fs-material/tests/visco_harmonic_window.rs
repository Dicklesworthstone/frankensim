//! Finite-window original material storage, work and heat from real-field ODEs.

use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::{
    Mandel6, PronyBranch, TensorProny, TensorPronyHarmonicWindowPoint, ViscoError,
};
use fs_math::c64::C64;

fn outer(diagonal: f64, a: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| a[i] * a[j] + if i == j { diagonal } else { 0.0 })
    })
}

fn law() -> TensorProny {
    TensorProny::new_mandel(
        outer(8.0, [1.0, -0.5, 0.2, 0.4, -0.3, 0.6]),
        vec![
            PronyBranch::new(outer(0.0, [1.0, 2.0, -0.5, 0.8, -0.4, 0.3]), 0.2).unwrap(),
            PronyBranch::new(outer(0.4, [-0.2, 0.3, 1.0, -0.6, 0.7, 0.8]), 1.4).unwrap(),
        ],
        0.05,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(outer(0.01, [0.1, -0.2, 0.05, 0.08, -0.03, 0.04]))
    .unwrap()
    .with_band((0.1, 100.0))
    .unwrap()
}

fn strain() -> [C64; 6] {
    core::array::from_fn(|i| C64::new((i as f64 - 2.0) * 0.001, (3.0 - i as f64) * 0.0007))
}

fn apply(matrix: &Mandel6, vector: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * vector[j]).sum())
}

fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn near(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}, relative={relative:e}"
    );
}

// Solve each Maxwell ODE in real cosine/sine coordinates, independently of
// the production spectrum and temporal Gram kernels. q=A*cos(wt)+B*sin(wt)
// satisfies [1,x;-x,1]*[A;B]=[x*epsilon_im;-x*epsilon_re], x=w*tau.
// Return instantaneous [storage, Kelvin power, Maxwell power, stress power].
fn real_observation(law: &TensorProny, omega: f64, e: &[C64; 6], t: f64) -> [f64; 4] {
    let (s, c) = (omega * t).sin_cos();
    let eps = e.map(|z| z.re * c + z.im * s);
    let rate = e.map(|z| omega * (z.im * c - z.re * s));
    let mut stress = apply(law.c_eq(), &eps);
    let mut storage = 0.5 * dot(&eps, &stress);
    let mut kelvin = 0.0;
    if let Some(viscosity) = law.kv_pa_s() {
        let viscous = apply(viscosity, &rate);
        kelvin = dot(&rate, &viscous);
        for i in 0..6 {
            stress[i] += viscous[i];
        }
    }
    let mut maxwell = 0.0;
    for branch in law.branches() {
        let x = omega * branch.tau_s;
        let denominator = 1.0 + x * x;
        let a = e.map(|z| (x * z.im + x * x * z.re) / denominator);
        let b = e.map(|z| (x * x * z.im - x * z.re) / denominator);
        let q = core::array::from_fn(|i| a[i] * c + b[i] * s);
        let q_rate: [f64; 6] = core::array::from_fn(|i| omega * (b[i] * c - a[i] * s));
        for i in 0..6 {
            near(
                q_rate[i] + q[i] / branch.tau_s,
                rate[i],
                omega * e[i].abs().max(1e-12),
                3e-14,
            );
        }
        let spring_stress = apply(&branch.c_pa, &q);
        let squared = dot(&q, &spring_stress);
        storage += 0.5 * squared;
        maxwell += squared / branch.tau_s;
        for i in 0..6 {
            stress[i] += spring_stress[i];
        }
    }
    [storage, kelvin, maxwell, dot(&stress, &rate)]
}

fn simpson(law: &TensorProny, omega: f64, strain: &[C64; 6], window: HarmonicWindow) -> [f64; 3] {
    let panels = 4096;
    let mut integrated = [0.0; 3];
    for i in 0..=panels {
        let t = window.start_s() + window.duration_s() * i as f64 / panels as f64;
        let sample = real_observation(law, omega, strain, t);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        for k in 0..3 {
            integrated[k] += weight * sample[k + 1];
        }
    }
    integrated.map(|value| value * window.duration_s() / (3.0 * panels as f64))
}

fn quantities(p: &TensorPronyHarmonicWindowPoint) -> [f64; 8] {
    [
        p.start_stored_energy_j_m3,
        p.end_stored_energy_j_m3,
        p.stored_energy_change_j_m3,
        p.kelvin_heat_j_m3,
        p.maxwell_heat_j_m3,
        p.material_heat_j_m3,
        p.stress_work_j_m3,
        p.energy_balance_residual_j_m3,
    ]
}

#[test]
fn anisotropic_subcycle_storage_heat_and_work_match_real_maxwell_odes_and_quadrature() {
    let law = law();
    let original = law.clone();
    let omega = 2.7;
    let strain = strain();
    let material = law.prepare_harmonic_material(omega).unwrap();
    for (start, duration) in [(0.17, 0.07), (0.19, 0.37), (-0.2, 0.83)] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let prepared = material.prepare_window(window).unwrap();
        let point = prepared.evaluate(&strain).unwrap();
        let start = real_observation(&law, omega, &strain, window.start_s())[0];
        let end = real_observation(&law, omega, &strain, window.start_s() + duration)[0];
        let [kelvin, maxwell, work] = simpson(&law, omega, &strain, window);
        near(point.start_stored_energy_j_m3, start, start, 3e-13);
        near(point.end_stored_energy_j_m3, end, end, 3e-13);
        near(
            point.stored_energy_change_j_m3,
            end - start,
            start.max(end),
            3e-13,
        );
        near(point.kelvin_heat_j_m3, kelvin, kelvin, 3e-12);
        near(point.maxwell_heat_j_m3, maxwell, maxwell, 3e-12);
        near(
            point.material_heat_j_m3,
            kelvin + maxwell,
            kelvin + maxwell,
            3e-12,
        );
        near(
            point.stress_work_j_m3,
            work,
            start.max(end).max(work.abs()),
            3e-12,
        );
        near(
            point.energy_balance_residual_j_m3,
            0.0,
            start.max(end),
            3e-13,
        );
        assert!(point.kelvin_heat_j_m3 > 0.0 && point.maxwell_heat_j_m3 > 0.0);
        assert_eq!(point, material.observe_window(&strain, window).unwrap());
    }
    // This fixture needs mutual tensor entries: dropping them changes the
    // equilibrium energy by far more than the quadrature/roundoff tolerance.
    let eps = strain.map(|z| z.re);
    let complete = dot(&eps, &apply(law.c_eq(), &eps));
    let diagonal: f64 = (0..6).map(|i| eps[i] * law.c_eq()[i][i] * eps[i]).sum();
    assert!((complete - diagonal).abs() > 1e-3 * complete);
    assert_eq!(law, original);
}

#[test]
fn full_period_recovers_original_cycle_heat_and_storage_returns_to_its_start() {
    let law = law();
    let omega = 3.4;
    let e = strain();
    let material = law.prepare_harmonic_material(omega).unwrap();
    let mean = material.evaluate(&e).unwrap();
    let duration = core::f64::consts::TAU / omega;
    let point = material
        .observe_window(&e, HarmonicWindow::new(0.23, duration).unwrap())
        .unwrap();
    near(
        point.kelvin_heat_j_m3 / duration,
        mean.mean_kelvin_heat_w_m3,
        mean.mean_kelvin_heat_w_m3,
        4e-13,
    );
    near(
        point.maxwell_heat_j_m3 / duration,
        mean.mean_maxwell_heat_w_m3,
        mean.mean_maxwell_heat_w_m3,
        4e-13,
    );
    near(
        point.material_heat_j_m3 / duration,
        mean.mean_material_heat_w_m3,
        mean.mean_material_heat_w_m3,
        4e-13,
    );
    near(
        point.stress_work_j_m3 / duration,
        mean.mean_stress_power_w_m3,
        mean.mean_material_heat_w_m3,
        4e-13,
    );
    near(
        point.end_stored_energy_j_m3,
        point.start_stored_energy_j_m3,
        point.start_stored_energy_j_m3,
        4e-13,
    );
    near(
        point.stored_energy_change_j_m3,
        0.0,
        point.start_stored_energy_j_m3,
        4e-13,
    );
}

#[test]
fn partition_and_clock_phase_changes_preserve_the_same_physical_exposure() {
    let law = law();
    let omega = 5.2;
    let e = strain();
    let material = law.prepare_harmonic_material(omega).unwrap();
    let start = 0.1875;
    let duration = 0.4375;
    let whole = material
        .observe_window(&e, HarmonicWindow::new(start, duration).unwrap())
        .unwrap();
    let mut time = start;
    let mut sums = [0.0; 5];
    let mut previous = whole.start_stored_energy_j_m3;
    for width in [0.125, 0.25, 0.0625] {
        let part = material
            .observe_window(&e, HarmonicWindow::new(time, width).unwrap())
            .unwrap();
        near(part.start_stored_energy_j_m3, previous, previous, 3e-13);
        previous = part.end_stored_energy_j_m3;
        for (total, value) in sums.iter_mut().zip([
            part.stored_energy_change_j_m3,
            part.kelvin_heat_j_m3,
            part.maxwell_heat_j_m3,
            part.material_heat_j_m3,
            part.stress_work_j_m3,
        ]) {
            *total += value;
        }
        time += width;
    }
    let scale = whole
        .start_stored_energy_j_m3
        .max(whole.end_stored_energy_j_m3);
    near(previous, whole.end_stored_energy_j_m3, scale, 3e-13);
    for (actual, expected) in sums.into_iter().zip([
        whole.stored_energy_change_j_m3,
        whole.kelvin_heat_j_m3,
        whole.maxwell_heat_j_m3,
        whole.material_heat_j_m3,
        whole.stress_work_j_m3,
    ]) {
        near(actual, expected, scale.max(expected.abs()), 3e-13);
    }

    let shift = 0.734;
    let (s, c) = (omega * shift).sin_cos();
    let translated = e.map(|z| z * C64::new(c, s));
    let point = material
        .observe_window(
            &translated,
            HarmonicWindow::new(start + shift, duration).unwrap(),
        )
        .unwrap();
    for (actual, expected) in quantities(&point).into_iter().zip(quantities(&whole)) {
        near(actual, expected, scale.max(expected.abs()), 4e-13);
    }
}

#[test]
fn short_window_at_a_rate_zero_retains_cubic_heat_and_sub_ulp_storage_change() {
    let viscosity = 0.7;
    let law = TensorProny::new_mandel(outer(1.0, [0.0; 6]), vec![], 0.01)
        .unwrap()
        .with_kelvin_viscosity_mandel(outer(viscosity, [0.0; 6]))
        .unwrap();
    let (omega, amplitude, duration) = (1.9_f64, 0.003_f64, 1e-9_f64);
    let mut e = [C64::ZERO; 6];
    e[0] = C64::from_re(amplitude);
    let point = law
        .prepare_harmonic_material(omega)
        .unwrap()
        .observe_window(&e, HarmonicWindow::new(0.0, duration).unwrap())
        .unwrap();
    let heat = viscosity * amplitude.powi(2) * omega.powi(4) * duration.powi(3) / 3.0;
    let change = -0.5 * amplitude.powi(2) * omega.powi(2) * duration.powi(2);
    near(point.kelvin_heat_j_m3, heat, heat, 4e-13);
    assert!(point.kelvin_heat_j_m3 > 0.0);
    assert_eq!(point.maxwell_heat_j_m3, 0.0);
    assert_eq!(point.start_stored_energy_j_m3, point.end_stored_energy_j_m3);
    near(point.stored_energy_change_j_m3, change, change.abs(), 4e-13);
    near(point.stress_work_j_m3, change + heat, change.abs(), 4e-13);
    near(point.energy_balance_residual_j_m3, 0.0, change.abs(), 4e-13);
}

#[test]
fn maxwell_heat_survives_an_unrepresentable_original_spring_strain() {
    let law = TensorProny::new_mandel(
        outer(1.0, [0.0; 6]),
        vec![PronyBranch::new(outer(1e100, [0.0; 6]), 1e-300).unwrap()],
        0.01,
    )
    .unwrap();
    let material = law.prepare_harmonic_material(1.0).unwrap();
    let mut e = [C64::ZERO; 6];
    e[0] = C64::from_re(1e-40);
    assert_eq!(
        material
            .branch_response(0, &e)
            .unwrap()
            .spring_strain_mandel,
        [C64::ZERO; 6]
    );
    let start = 0.3_f64;
    let duration = 0.4_f64;
    let point = material
        .observe_window(&e, HarmonicWindow::new(start, duration).unwrap())
        .unwrap();
    let integral_sin_squared =
        duration * 0.5 - ((2.0 * (start + duration)).sin() - (2.0 * start).sin()) * 0.25;
    let heat = 1e-280 * integral_sin_squared;
    near(point.maxwell_heat_j_m3, heat, heat, 5e-13);
    assert!(point.maxwell_heat_j_m3 > 0.0);
    assert_eq!(point.kelvin_heat_j_m3, 0.0);
    assert_eq!(point.material_heat_j_m3, point.maxwell_heat_j_m3);
}

#[test]
fn long_interval_work_keeps_original_loss_factors_when_the_dynamic_loss_tensor_underflows() {
    // Both eta*w and C_branch*w*tau are 1e-400 and round out of the cached
    // dynamic loss tensor. Their actual positive work over this admitted long
    // interval remains finite: coefficient*w²*amplitude²*T/2 = 5e-297.
    for maxwell in [false, true] {
        let branches = if maxwell {
            vec![PronyBranch::new(outer(1e-200, [0.0; 6]), 1.0).unwrap()]
        } else {
            vec![]
        };
        let mut law = TensorProny::new_mandel(outer(1e-300, [0.0; 6]), branches, 0.05).unwrap();
        if !maxwell {
            law = law
                .with_kelvin_viscosity_mandel(outer(1e-200, [0.0; 6]))
                .unwrap();
        }
        let mut e = [C64::ZERO; 6];
        e[0] = C64::from_re(0.01);
        let omega = 1e-200;
        assert_eq!(law.complex_modulus(omega).unwrap().1[0][0], 0.0);
        let point = law
            .prepare_harmonic_material(omega)
            .unwrap()
            .observe_window(&e, HarmonicWindow::new(-5e307, 1e308).unwrap())
            .unwrap();
        let heat = 5e-297;
        near(point.material_heat_j_m3, heat, heat, 5e-13);
        near(point.stress_work_j_m3, heat, heat, 5e-13);
        assert!(point.stress_work_j_m3 > 0.0);
        near(point.stored_energy_change_j_m3, 0.0, heat, 5e-13);
        near(point.energy_balance_residual_j_m3, 0.0, heat, 5e-13);
        if maxwell {
            assert_eq!(point.kelvin_heat_j_m3, 0.0);
            near(point.maxwell_heat_j_m3, heat, heat, 5e-13);
        } else {
            assert_eq!(point.maxwell_heat_j_m3, 0.0);
            near(point.kelvin_heat_j_m3, heat, heat, 5e-13);
        }
    }
}

#[test]
fn zero_channels_are_exact_and_original_limits_refuse_without_mutating_the_law() {
    let law = law();
    let original = law.clone();
    assert!(matches!(
        law.prepare_harmonic_material(0.05),
        Err(ViscoError::OutOfBand { .. })
    ));
    let material = law.prepare_harmonic_material(2.7).unwrap();
    let window = HarmonicWindow::new(0.0, 1e-9).unwrap();
    let zero = material.observe_window(&[C64::ZERO; 6], window).unwrap();
    assert_eq!(quantities(&zero), [0.0; 8]);
    let triangular =
        HarmonicWindow::with_profile(0.0, 0.1, HarmonicWindowProfile::Triangular).unwrap();
    assert!(material.prepare_window(triangular).is_err());
    // The actual strain is near zero inside this tiny window, but its original
    // whole-cycle amplitude exceeds the source law's small-strain limit.
    let mut excessive = [C64::ZERO; 6];
    excessive[0] = C64::new(0.0, 0.1);
    assert!(material.observe_window(&excessive, window).is_err());
    assert!(
        material
            .observe_window(&[C64::from_re(f64::NAN); 6], window)
            .is_err()
    );
    assert_eq!(law, original);

    let elastic = TensorProny::new_mandel(outer(1.0, [0.0; 6]), vec![], 0.01).unwrap();
    let point = elastic
        .prepare_harmonic_material(2.7)
        .unwrap()
        .observe_window(&strain(), HarmonicWindow::new(0.13, 0.17).unwrap())
        .unwrap();
    assert_eq!(point.kelvin_heat_j_m3, 0.0);
    assert_eq!(point.maxwell_heat_j_m3, 0.0);
    assert_eq!(point.material_heat_j_m3, 0.0);
    near(
        point.stress_work_j_m3,
        point.stored_energy_change_j_m3,
        point
            .start_stored_energy_j_m3
            .max(point.end_stored_energy_j_m3),
        3e-13,
    );
}
