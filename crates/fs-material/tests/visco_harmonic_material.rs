//! Original harmonic stress/storage/heat against full tensor and history oracles.
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco::{Mandel6, PronyBranch, TensorProny, ViscoError};
use fs_math::{c64::C64, det};

fn outer(d: f64, a: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| a[i] * a[j] + if i == j { d } else { 0.0 }))
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
fn apply(a: &Mandel6, x: &[C64; 6]) -> [C64; 6] {
    core::array::from_fn(|i| (0..6).fold(C64::ZERO, |s, j| s + x[j].scale(a[i][j])))
}
fn norm(a: &Mandel6, x: &[C64; 6]) -> f64 {
    let ax = apply(a, x);
    (0..6).map(|i| (x[i].conj() * ax[i]).re).sum()
}
fn memory(w: f64, tau: f64, e: &[C64; 6]) -> [C64; 6] {
    let s = C64::new(0.0, -w * tau);
    let transfer = s / (C64::ONE + s);
    e.map(|v| transfer * v)
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1e-12), "{a:e} != {b:e}");
}

#[test]
fn full_tensor_memory_stress_and_independent_energies_match() {
    let law = law();
    let w = 2.7;
    let e = strain();
    let prepared = law.prepare_harmonic_material(w).unwrap();
    let point = prepared.evaluate(&e).unwrap();
    let mut stress = apply(law.c_eq(), &e);
    let rate = e.map(|v| C64::new(0.0, -w) * v);
    let kv = apply(law.kv_pa_s().unwrap(), &rate);
    for i in 0..6 {
        stress[i] = stress[i] + kv[i];
    }
    let mut energy = 0.25 * norm(law.c_eq(), &e);
    let mut heat = 0.5 * norm(law.kv_pa_s().unwrap(), &rate);
    for (i, b) in law.branches().iter().enumerate() {
        let q = memory(w, b.tau_s, &e);
        let branch = prepared.branch_response(i, &e).unwrap();
        let branch_stress = apply(&b.c_pa, &q);
        for j in 0..6 {
            assert!((branch.spring_strain_mandel[j] - q[j]).abs() < 2e-17);
            let residual = (C64::ONE + C64::new(0.0, -w * b.tau_s))
                * branch.spring_strain_mandel[j]
                - C64::new(0.0, -w * b.tau_s) * e[j];
            assert!(residual.abs() < 2e-17);
            assert!((branch.stress_mandel_pa[j] - branch_stress[j]).abs() < 2e-16);
            stress[j] = stress[j] + branch_stress[j];
        }
        near(
            branch.mean_stored_energy_j_m3,
            0.25 * norm(&b.c_pa, &q),
            2e-13,
        );
        near(
            branch.mean_heat_w_m3,
            0.5 * norm(&b.c_pa, &q) / b.tau_s,
            2e-13,
        );
        energy += 0.25 * norm(&b.c_pa, &q);
        heat += 0.5 * norm(&b.c_pa, &q) / b.tau_s;
    }
    for i in 0..6 {
        assert!((point.stress_mandel_pa[i] - stress[i]).abs() < 3e-16);
    }
    near(point.mean_stored_energy_j_m3, energy, 2e-13);
    near(point.mean_material_heat_w_m3, heat, 2e-13);
    near(point.mean_stress_power_w_m3, heat, 2e-13);
    assert!(point.power_balance_residual_w_m3.abs() < 3e-13 * heat);
    let phase = C64::new(0.6, 0.8);
    let rotated = prepared.evaluate(&e.map(|v| v * phase)).unwrap();
    near(rotated.mean_stored_energy_j_m3, energy, 2e-13);
    near(rotated.mean_material_heat_w_m3, heat, 2e-13);
}

#[test]
fn original_piecewise_linear_history_converges_to_the_same_cycle_energy_and_heat() {
    let law = law();
    let w = 2.7;
    let e = strain();
    let period = core::f64::consts::TAU / w;
    let exact = law
        .prepare_harmonic_material(w)
        .unwrap()
        .evaluate(&e)
        .unwrap();
    let run = |n: usize| {
        let mut state = law.state();
        state.strain = e.map(|v| v.re);
        state.q = law
            .branches()
            .iter()
            .map(|b| memory(w, b.tau_s, &e).map(|v| v.re))
            .collect();
        let real = state.strain.map(C64::from_re);
        state.stored = 0.5 * norm(law.c_eq(), &real);
        for (b, q) in law.branches().iter().zip(&state.q) {
            state.stored += 0.5 * norm(&b.c_pa, &q.map(C64::from_re));
        }
        state.stress = exact.stress_mandel_pa.map(|v| v.re);
        let mut energy = 0.5 * state.stored;
        for step in 1..=n {
            let (sin, cos) = det::sin_cos(core::f64::consts::TAU * step as f64 / n as f64);
            let eps = e.map(|v| v.re * cos + v.im * sin);
            law.step_with_mean(&mut state, eps, period / n as f64)
                .unwrap();
            energy += if step == n {
                0.5 * state.stored
            } else {
                state.stored
            };
        }
        (
            (energy / n as f64 - exact.mean_stored_energy_j_m3).abs(),
            (state.dissipated / period - exact.mean_material_heat_w_m3).abs(),
        )
    };
    let coarse = run(128);
    let fine = run(512);
    assert!(fine.0 < coarse.0 / 8.0);
    assert!(fine.1 < coarse.1 / 8.0);
    assert!(fine.0 < 1e-4 * exact.mean_stored_energy_j_m3);
    assert!(fine.1 < 1e-4 * exact.mean_material_heat_w_m3);
}

#[test]
fn constrained_material_observation_keeps_all_transverse_stress_zero() {
    let law = law();
    let w = 4.0;
    let section = law
        .prepare_plane_stress_harmonic(MaterialPlane::Axes13, w)
        .unwrap();
    let strain = section
        .recover_strain(&[
            C64::new(0.001, 0.002),
            C64::new(-0.002, 0.001),
            C64::new(0.003, -0.001),
        ])
        .unwrap();
    let point = law
        .prepare_harmonic_material(w)
        .unwrap()
        .evaluate(&strain)
        .unwrap();
    for i in [1, 3, 4] {
        assert!(point.stress_mandel_pa[i].abs() < 5e-16);
    }
    assert!(point.mean_kelvin_heat_w_m3 > 0.0 && point.mean_maxwell_heat_w_m3 > 0.0);
    near(
        point.mean_stress_power_w_m3,
        point.mean_material_heat_w_m3,
        2e-13,
    );
}

#[test]
fn original_frequency_strain_and_branch_domains_refuse_without_changing_source() {
    let law = law();
    let original = law.clone();
    assert!(law.prepare_harmonic_material(0.0).is_err());
    assert!(matches!(
        law.prepare_harmonic_material(0.05),
        Err(ViscoError::OutOfBand { .. })
    ));
    let prepared = law.prepare_harmonic_material(3.0).unwrap();
    assert!(prepared.evaluate(&[C64::new(0.04, 0.04); 6]).is_err());
    assert!(prepared.evaluate(&[C64::from_re(f64::NAN); 6]).is_err());
    assert!(prepared.branch_response(2, &strain()).is_err());
    assert_eq!(
        prepared
            .evaluate(&[C64::ZERO; 6])
            .unwrap()
            .mean_material_heat_w_m3,
        0.0
    );
    assert_eq!(law, original);
}

#[test]
fn original_stress_and_heat_survive_spring_strain_underflow() {
    let stiffness =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1e100 } else { 0.0 }));
    let law = TensorProny::new_mandel(
        outer(1.0, [0.0; 6]),
        vec![PronyBranch::new(stiffness, 1e-300).unwrap()],
        0.01,
    )
    .unwrap();
    let prepared = law.prepare_harmonic_material(1.0).unwrap();
    let mut strain = [C64::ZERO; 6];
    strain[0] = C64::from_re(1e-40);
    let branch = prepared.branch_response(0, &strain).unwrap();
    assert_eq!(branch.spring_strain_mandel, [C64::ZERO; 6]);
    assert!((branch.stress_mandel_pa[0].im / -1e-240 - 1.0).abs() < 2e-14);
    assert!((branch.mean_heat_w_m3 / 5e-281 - 1.0).abs() < 2e-14);
    let point = prepared.evaluate(&strain).unwrap();
    assert!((point.mean_material_heat_w_m3 / 5e-281 - 1.0).abs() < 2e-14);
    assert!((point.mean_stress_power_w_m3 / 5e-281 - 1.0).abs() < 2e-14);
}
