//! Independent synthetic joint creep/yield oracles, not measured calibration.

use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{
    FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State,
};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.1,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 3.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}
fn backbone(hardening: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: hardening,
        },
        limits(),
    )
    .unwrap()
}
fn kinetics(eta: f64, n: f64) -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: eta,
        stress_scale_pa: 2e6,
        exponent: n,
    }
}
fn law(eta: f64, n: f64, dt: f64) -> FiniteNortonJ2 {
    FiniteNortonJ2::new(backbone(H), kinetics(eta, n), dt, 3.0).unwrap()
}
fn coaxial(mean: f64, axial: f64) -> [f64; 9] {
    [
        (mean + axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
    ]
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}, difference={:.8e}",
        (actual - expected).abs()
    );
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn check_tangent(law: &FiniteNortonJ2, f: &[f64; 9], old: &FiniteJ2State) {
    let r = law.response(f, old).unwrap();
    for column in 0..9 {
        let h = 2e-6;
        let mut plus = *f;
        plus[column] += h;
        let mut minus = *f;
        minus[column] -= h;
        let (pp, wp) = law.stress_and_potential(&plus, old).unwrap();
        let (pm, wm) = law.stress_and_potential(&minus, old).unwrap();
        near((wp - wm) / (2.0 * h), r.first_piola_pa[column], 0.25, 3e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * h),
                r.tangent_pa[row][column],
                2.0,
                3e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                2e-12,
            );
        }
    }
}

#[test]
fn independent_linear_joint_return_separates_mechanisms_storage_heat_and_potential() {
    let (eta, dt, mean) = (10e6, 0.1, 0.03);
    let law = law(eta, 1.0, dt);
    let old = law.initial_state();
    for (axial, yielding) in [(0.004, false), (0.12, true)] {
        let f = coaxial(mean, axial);
        let r = law.response(&f, &old).unwrap();
        let h = 3.0 * MU;
        let v = eta / dt;
        let trial = h * axial;
        let q = if yielding {
            let beta = H / (h + H);
            (Y + beta * (trial - Y)) / (1.0 + h * beta / v)
        } else {
            trial / (1.0 + h / v)
        };
        let dc = q / v;
        let dp = if yielding {
            (trial - Y - h * dc) / (h + H)
        } else {
            0.0
        };
        let elastic = axial - dp - dc;
        let stored =
            0.5 * K * (3.0 * mean).powi(2) + 1.5 * MU * elastic * elastic + 0.5 * H * dp * dp;
        near(r.plastic_increment, dp, 1e-14, 2e-13);
        near(r.creep_increment, dc, 1e-14, 2e-13);
        near(r.state.alpha, dp, 1e-14, 2e-13);
        near(r.state.creep_alpha, dc, 1e-14, 2e-13);
        near(r.equivalent_kirchhoff_pa, q, 2e-6, 2e-13);
        near(r.jacobian, (3.0 * mean).exp(), 1e-14, 1e-14);
        near(
            r.state.inverse_plastic_cauchy_green[0],
            (-2.0 * (dp + dc)).exp(),
            1e-14,
            2e-13,
        );
        for i in [1, 2] {
            near(
                r.state.inverse_plastic_cauchy_green[i],
                (dp + dc).exp(),
                1e-14,
                2e-13,
            );
        }
        for (i, tau) in [
            (0, 3.0 * K * mean + 2.0 * MU * elastic),
            (4, 3.0 * K * mean - MU * elastic),
            (8, 3.0 * K * mean - MU * elastic),
        ] {
            near(r.kirchhoff_pa[i], tau, 2e-6, 2e-13);
            near(r.first_piola_pa[i], tau / f[i], 2e-6, 2e-13);
        }
        near(r.stored_energy_j_m3, stored, 2e-7, 2e-13);
        near(r.plastic_dissipation_j_m3, Y * dp, 2e-7, 2e-13);
        near(r.creep_dissipation_j_m3, q * dc, 2e-7, 2e-13);
        near(r.inelastic_dissipation_j_m3(), Y * dp + q * dc, 2e-7, 2e-13);
        near(
            r.incremental_potential_j_m3,
            stored + Y * dp + 0.5 * q * dc,
            2e-7,
            2e-13,
        );
        assert_eq!(r.plastic_increment > 0.0, yielding);
        assert!(r.creep_increment > 0.0);
        assert_eq!(
            law.stress_and_potential(&f, &old).unwrap(),
            (r.first_piola_pa, r.incremental_potential_j_m3)
        );
        law.validate_committed(&f, &r.state, 1e-4).unwrap();
        check_tangent(&law, &f, &old);
    }
    assert_eq!(old, law.initial_state());
}

#[test]
fn sustained_subyield_load_uses_only_creep_bound_and_never_expands_yield_or_stores_hardening() {
    let (q, eta, dt) = (0.8e6, 8e6, 0.2);
    let native = FiniteJ2::new(
        *backbone(H).parameters(),
        FiniteJ2Limits {
            maximum_accumulated_plastic_strain: 0.02,
            ..limits()
        },
    )
    .unwrap();
    let law = FiniteNortonJ2::new(native, kinetics(eta, 1.0), dt, 1.0).unwrap();
    let mut state = law.initial_state();
    let dc = dt * q / eta;
    let mut final_f = ID;
    for step in 1..=20 {
        let accumulated = f64::from(step) * dc;
        let f = coaxial(0.01, accumulated + q / (3.0 * MU));
        let r = law.response(&f, &state).unwrap();
        near(r.state.creep_alpha, accumulated, 2e-13, 2e-12);
        near(r.creep_increment, dc, 2e-13, 2e-12);
        near(r.equivalent_kirchhoff_pa, q, 3e-5, 2e-12);
        assert_eq!(r.plastic_increment, 0.0);
        assert_eq!(r.plastic_dissipation_j_m3, 0.0);
        assert_eq!(r.state.alpha, 0.0);
        near(r.creep_dissipation_j_m3, q * dc, 2e-6, 2e-11);
        near(
            r.stored_energy_j_m3,
            0.5 * K * 0.03_f64.powi(2) + q * q / (6.0 * MU),
            2e-6,
            2e-11,
        );
        state = r.state;
        final_f = f;
    }
    assert!(state.creep_alpha > 0.39 && final_f[0] > 1.5);
    law.validate_committed(&final_f, &state, 1e-4).unwrap();
    let mut closed_cycle = state.clone();
    closed_cycle.creep_alpha += 0.1;
    assert_eq!(
        law.stored_energy(&final_f, &state).unwrap(),
        law.stored_energy(&final_f, &closed_cycle).unwrap()
    );
    // Existing families cannot reinterpret the new, separate mechanism history.
    assert!(backbone(H).validate_history(&state).is_err());
    assert!(
        FiniteOverstressJ2::new(backbone(H), kinetics(eta, 1.0), dt)
            .unwrap()
            .validate_history(&state)
            .is_err()
    );
    assert!(
        FiniteNortonCreep::new(MU, K, limits(), kinetics(eta, 1.0), dt)
            .unwrap()
            .validate_history(&state)
            .is_err()
    );
}

#[test]
fn initial_yielding_then_shape_hold_relaxes_below_yield_without_more_plastic_hardening() {
    let eta = 400e6;
    let f = coaxial(0.0, 0.12);
    let initial_law = law(eta, 1.0, 0.02);
    let initial = initial_law
        .response(&f, &initial_law.initial_state())
        .unwrap();
    assert!(initial.plastic_increment > 0.0 && initial.creep_increment > 0.0);
    // After the initial yield event, a sustained load at that same radius is
    // carried by creep alone. Reducing the load retains all plastic hardening.
    let loaded = law(eta, 1.0, 0.2);
    let mut loaded_state = initial.state.clone();
    for factor in [1.0, 1.0, 0.4] {
        let q = factor * initial.equivalent_kirchhoff_pa;
        let dc = loaded.duration_s() * q / eta;
        let axial = loaded_state.alpha + loaded_state.creep_alpha + dc + q / (3.0 * MU);
        let next = loaded
            .response(&coaxial(0.0, axial), &loaded_state)
            .unwrap();
        assert_eq!(next.plastic_increment, 0.0);
        assert_eq!(next.state.alpha, initial.state.alpha);
        near(next.creep_increment, dc, 1e-14, 2e-12);
        near(next.equivalent_kirchhoff_pa, q, 3e-5, 2e-12);
        loaded_state = next.state;
    }
    let exact = initial.equivalent_kirchhoff_pa * (-3.0 * MU * 2.0 / eta).exp();
    let mut previous_error = f64::INFINITY;
    for count in [10, 20, 40] {
        let dt = 2.0 / f64::from(count);
        let law = law(eta, 1.0, dt);
        let mut state = initial.state.clone();
        let mut q = initial.equivalent_kirchhoff_pa;
        let mut energy = law.stored_energy(&f, &state).unwrap();
        for _ in 0..count {
            let next = law.response(&f, &state).unwrap();
            near(
                next.equivalent_kirchhoff_pa,
                q / (1.0 + 3.0 * MU * dt / eta),
                3e-5,
                2e-12,
            );
            assert_eq!(next.plastic_increment, 0.0);
            assert_eq!(next.state.alpha, initial.state.alpha);
            assert!(next.state.creep_alpha > state.creep_alpha);
            assert!(next.equivalent_kirchhoff_pa < Y + H * next.state.alpha);
            let defect = energy - next.stored_energy_j_m3 - next.inelastic_dissipation_j_m3();
            near(defect, 1.5 * MU * next.creep_increment.powi(2), 2e-6, 2e-10);
            assert!(defect > 0.0 && next.creep_dissipation_j_m3 > 0.0);
            q = next.equivalent_kirchhoff_pa;
            energy = next.stored_energy_j_m3;
            state = next.state;
        }
        let error = (q - exact).abs();
        if previous_error.is_finite() {
            assert!(error < 0.55 * previous_error);
        }
        previous_error = error;
    }
}

#[test]
fn perfectly_plastic_limit_jointly_flows_and_resolves_the_creep_only_transition() {
    let (eta, dt) = (10e6, 0.1);
    let law = FiniteNortonJ2::new(backbone(0.0), kinetics(eta, 1.0), dt, 3.0).unwrap();
    let old = law.initial_state();
    let f = coaxial(0.0, 0.12);
    let r = law.response(&f, &old).unwrap();
    let dc = Y * dt / eta;
    near(r.equivalent_kirchhoff_pa, Y, 2e-6, 2e-13);
    near(r.creep_increment, dc, 1e-14, 2e-13);
    near(
        r.plastic_increment,
        0.12 - Y / (3.0 * MU) - dc,
        1e-14,
        2e-13,
    );
    near(r.stored_energy_j_m3, Y * Y / (6.0 * MU), 2e-7, 2e-13);
    check_tangent(&law, &f, &old);
    let switch = coaxial(0.0, Y / (3.0 * MU) + dc);
    let boundary = law.response(&switch, &old).unwrap();
    assert_eq!(boundary.plastic_increment, 0.0);
    near(boundary.creep_increment, dc, 1e-14, 2e-13);
    near(boundary.equivalent_kirchhoff_pa, Y, 2e-6, 2e-13);
    let held = law.response(&f, &r.state).unwrap();
    assert_eq!(held.plastic_increment, 0.0);
    assert_eq!(held.state.alpha, r.state.alpha);
    assert!(held.equivalent_kirchhoff_pa < Y && held.creep_increment > 0.0);
}

#[test]
fn nonlinear_joint_return_has_objective_history_reversal_and_full_potential_tangents() {
    let (eta, n, dt) = (20e6, 2.7, 0.3);
    let law = law(eta, n, dt);
    let initial = law.initial_state();
    let old = law.response(&coaxial(0.01, 0.07), &initial).unwrap().state;
    let f = [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96];
    let r = law.response(&f, &old).unwrap();
    assert!(r.plastic_increment > 0.0 && r.creep_increment > 0.0);
    let q = 2e6 * (eta * r.creep_increment / (dt * 2e6)).powf(1.0 / n);
    near(r.equivalent_kirchhoff_pa, q, 2e-6, 2e-12);
    near(q, Y + H * r.state.alpha, 2e-6, 2e-12);
    near(
        r.plastic_dissipation_j_m3,
        Y * r.plastic_increment,
        2e-7,
        2e-12,
    );
    near(r.creep_dissipation_j_m3, q * r.creep_increment, 2e-7, 2e-12);
    near(
        r.incremental_potential_j_m3,
        r.stored_energy_j_m3 + Y * r.plastic_increment + (n / (n + 1.0)) * q * r.creep_increment,
        2e-7,
        2e-12,
    );
    check_tangent(&law, &f, &old);
    let subyield = coaxial(0.01, old.alpha + old.creep_alpha + 0.001);
    let inactive = law.response(&subyield, &old).unwrap();
    assert_eq!(inactive.plastic_increment, 0.0);
    assert!(inactive.creep_increment > 0.0);
    check_tangent(&law, &subyield, &old);
    let (s, c) = 1.15_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = law.response(&multiply(&rotation, &f), &old).unwrap();
    near(rotated.state.alpha, r.state.alpha, 2e-13, 2e-12);
    near(rotated.state.creep_alpha, r.state.creep_alpha, 2e-13, 2e-12);
    for i in 0..6 {
        near(
            rotated.state.inverse_plastic_cauchy_green[i],
            r.state.inverse_plastic_cauchy_green[i],
            2e-13,
            2e-12,
        );
    }
    let expected = multiply(&rotation, &r.first_piola_pa);
    for i in 0..9 {
        near(rotated.first_piola_pa[i], expected[i], 3e-5, 3e-12);
    }
    let reverse = law.response(&coaxial(0.0, -0.10), &old).unwrap();
    assert!(reverse.state.alpha > old.alpha && reverse.state.creep_alpha > old.creep_alpha);
    assert!(reverse.state.inverse_plastic_cauchy_green[0] > old.inverse_plastic_cauchy_green[0]);
    assert!(reverse.kirchhoff_pa[0] < 0.0);
    assert!(reverse.plastic_dissipation_j_m3 > 0.0 && reverse.creep_dissipation_j_m3 > 0.0);
    law.validate_committed(&f, &r.state, 1e-4).unwrap();
}

#[test]
fn zero_stress_has_the_exact_linear_relaxation_tangent_and_positive_resolvable_creep() {
    let (eta, dt) = (60e6, 0.4);
    let law = law(eta, 1.0, dt);
    let r = law.response(&ID, &law.initial_state()).unwrap();
    let ratio = (eta / dt) / (3.0 * MU + eta / dt);
    assert_eq!(r.state, law.initial_state());
    assert_eq!(r.first_piola_pa, [0.0; 9]);
    assert_eq!(r.inelastic_dissipation_j_m3(), 0.0);
    for ij in 0..9 {
        let (i, j) = (ij / 3, ij % 3);
        for kl in 0..9 {
            let (k, l) = (kl / 3, kl % 3);
            let delta = |a: usize, b: usize| if a == b { 1.0 } else { 0.0 };
            let expected = K * delta(i, j) * delta(k, l)
                + MU * ratio
                    * (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k)
                        - (2.0 / 3.0) * delta(i, j) * delta(k, l));
            near(r.tangent_pa[ij][kl], expected, 2e-6, 2e-14);
        }
    }
    let tiny = coaxial(0.0, 1e-15);
    let logs = [tiny[0].ln(), tiny[4].ln(), tiny[8].ln()];
    let mean = logs.iter().sum::<f64>() / 3.0;
    let dev = logs.map(|x| x - mean);
    let q_trial = (6.0 * dev.iter().map(|x| x * x).sum::<f64>()).sqrt() * MU;
    let small = law.response(&tiny, &law.initial_state()).unwrap();
    assert!(small.creep_increment > 0.0 && small.plastic_increment == 0.0);
    near(
        small.creep_increment,
        q_trial / (3.0 * MU + eta / dt),
        1e-29,
        1e-11,
    );
    let nonlinear = FiniteNortonJ2::new(backbone(H), kinetics(eta, 3.0), dt, 3.0).unwrap();
    let n = nonlinear.response(&ID, &nonlinear.initial_state()).unwrap();
    near(n.tangent_pa[1][1], MU, 1e-7, 1e-14);
}

#[test]
fn separate_history_bounds_and_unresolved_returns_refuse_without_mutation() {
    for maximum in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(FiniteNortonJ2::new(backbone(H), kinetics(60e6, 1.0), 0.4, maximum).is_err());
    }
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(FiniteNortonJ2::new(backbone(H), kinetics(60e6, 1.0), dt, 3.0).is_err());
    }
    let ordinary = law(60e6, 1.0, 0.4);
    let f = coaxial(0.0, 0.12);
    let old = ordinary
        .response(&f, &ordinary.initial_state())
        .unwrap()
        .state;
    let saved = old.clone();
    let unresolved = law(1e30, 1.0, 1.0);
    assert!(matches!(
        unresolved.response(&f, &old),
        Err(FiniteJ2Error::Numerical(_))
    ));
    let bounded = FiniteNortonJ2::new(
        backbone(H),
        kinetics(60e6, 1.0),
        0.4,
        old.creep_alpha + 1e-5,
    )
    .unwrap();
    assert!(matches!(
        bounded.response(&f, &old),
        Err(FiniteJ2Error::Domain { .. })
    ));
    for invalid_value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = old.clone();
        invalid.creep_alpha = invalid_value;
        assert!(ordinary.validate_history(&invalid).is_err());
    }
    let mut unreachable = old.clone();
    unreachable.alpha = 0.0;
    unreachable.creep_alpha = 1e-6;
    assert!(ordinary.validate_history(&unreachable).is_err());
    assert!(ordinary.response(&[0.0; 9], &old).is_err());
    assert!(
        ordinary
            .validate_committed(&coaxial(0.0, 0.4), &old, 1e-4)
            .is_err()
    );
    ordinary.validate_committed(&f, &old, 1e-4).unwrap();
    // Inspection remains available, but exposing the backbone cannot silently
    // drop the selected creep kinetics through the static or Perzyna adapter.
    let exposed = ordinary.backbone();
    exposed.validate_history(&old).unwrap();
    assert!(exposed.response(&f, &old).is_err());
    assert!(exposed.stress_and_potential(&f, &old).is_err());
    let wrong_law = FiniteOverstressJ2::new(exposed.clone(), kinetics(60e6, 1.0), 0.4).unwrap();
    assert!(wrong_law.response(&f, &old).is_err());
    assert!(wrong_law.stress_and_potential(&f, &old).is_err());
    assert!(
        wrong_law
            .validate_committed(&coaxial(0.0, 0.4), &old, 1e-4)
            .is_err()
    );
    wrong_law.validate_committed(&f, &old, 1e-4).unwrap();
    assert_eq!(old, saved);
}
