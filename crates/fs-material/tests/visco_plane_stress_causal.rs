//! Original anisotropic plane-stress histories under coherent harmonic loading.
//! The oracle integrates the six-dimensional source equation in physical time;
//! it does not use the production Gram modes, harmonic moments or transfer map.
use fs_fft::harmonic_window::phase_angular_interval_from_clock_dd;
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco_tensor::stress_control::{
    PreparedTensorMixedControl as Mixed, TensorPlaneStressHarmonicBudget as HarmonicBudget,
    TensorStressControlBudget,
};
use fs_material::visco_tensor::{
    TensorProny, TensorPronyBranch, TensorPronyError, TensorPronyHistory,
};
use fs_material::{Tangent6, Voigt};
use fs_math::{c64::C64, dd::Dd};

const TAUS: [f64; 2] = [0.23, 1.1];
const WEIGHTS: Voigt = [
    1.0,
    1.0,
    1.0,
    std::f64::consts::SQRT_2,
    std::f64::consts::SQRT_2,
    std::f64::consts::SQRT_2,
];
const FREQUENCIES: [f64; 2] = [2.3, 3.9];
const STRAINS: [[C64; 3]; 2] = [
    [
        C64::new(0.0012, -0.0007),
        C64::new(-0.0005, 0.0009),
        C64::new(0.0006, 0.0011),
    ],
    [
        C64::new(-0.0008, 0.0004),
        C64::new(0.0007, -0.0003),
        C64::new(-0.0009, 0.0002),
    ],
];
const BIAS: [f64; 3] = [0.001, -0.0005, 0.0007];

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        maximum_payload_bytes: 4_000_000,
        maximum_scalar_work: 5_000_000,
    }
}
fn axes(plane: MaterialPlane) -> ([usize; 3], [usize; 3]) {
    match plane {
        MaterialPlane::Axes12 => ([0, 1, 3], [2, 4, 5]),
        MaterialPlane::Axes21 => ([1, 0, 3], [2, 4, 5]),
        MaterialPlane::Axes13 => ([0, 2, 5], [1, 3, 4]),
        MaterialPlane::Axes23 => ([1, 2, 4], [0, 3, 5]),
    }
}
fn mask(plane: MaterialPlane) -> [bool; 6] {
    let (_, outside) = axes(plane);
    core::array::from_fn(|i| outside.contains(&i))
}
fn prescribed(plane: MaterialPlane, engineering: [f64; 3]) -> Voigt {
    let (inside, _) = axes(plane);
    let mut values = [0.0; 6];
    for j in 0..3 {
        values[inside[j]] = engineering[j] / if j == 2 { 2.0 } else { 1.0 };
    }
    values
}
fn diagonal(values: Voigt) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn apply(a: &Tangent6, x: Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| a[i][j] * x[j]).sum())
}
fn dot(a: Voigt, b: Voigt) -> f64 {
    (0..6).map(|i| a[i] * b[i]).sum()
}
fn mandel(v: Voigt) -> Voigt {
    core::array::from_fn(|i| v[i] * WEIGHTS[i])
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

// Noncommuting anisotropic tensors: both branch factors mix normal and shear
// strains, and one is rank one. No simultaneous scalar material reduction exists.
fn coupled() -> (TensorProny, OriginalOde) {
    let mut eq = diagonal([1200.0, 1700.0, 1000.0, 900.0, 700.0, 1300.0]);
    eq[0][2] = -130.0;
    eq[2][0] = -130.0;
    eq[1][4] = 90.0;
    eq[4][1] = 90.0;
    let mut factors = [[[0.0; 6]; 6]; 2];
    factors[0][0] = [17.0, -3.0, 14.0, 7.0, 1.0, -2.0];
    factors[1][0] = [1.0, 12.0, -2.0, 6.0, 9.0, -4.0];
    factors[1][1] = [-3.0, 2.0, 8.0, -4.0, 0.0, 11.0];
    let branches = factors.map(|b| {
        core::array::from_fn(|i| core::array::from_fn(|j| (0..6).map(|k| b[k][i] * b[k][j]).sum()))
    });
    let c0 = core::array::from_fn(|i| {
        core::array::from_fn(|j| eq[i][j] + branches[0][i][j] + branches[1][i][j])
    });
    let arms = factors
        .into_iter()
        .zip(TAUS)
        .map(|(b, tau)| TensorPronyBranch::from_gram_rows(b, tau).unwrap())
        .collect();
    let law = TensorProny::new_mandel(eq, arms, 0.1)
        .unwrap()
        .with_band((0.0, 100.0))
        .unwrap();
    (law, OriginalOde { eq, branches, c0 })
}

fn waveform(
    frequencies: &[f64],
    strains: &[[C64; 3]],
    bias: [f64; 3],
    t: f64,
) -> ([f64; 3], [f64; 3]) {
    let mut value = bias;
    let mut rate = [0.0; 3];
    for (&omega, strain) in frequencies.iter().zip(strains) {
        let (sine, cosine) = (omega * t).sin_cos();
        for i in 0..3 {
            value[i] += strain[i].re * cosine + strain[i].im * sine;
            rate[i] += omega * (-strain[i].re * sine + strain[i].im * cosine);
        }
    }
    (value, rate)
}
fn at_clock(frequencies: &[f64], strains: &[[C64; 3]], bias: [f64; 3], start: Dd) -> [f64; 3] {
    let mut value = bias.map(Dd::from_f64);
    for (&omega, strain) in frequencies.iter().zip(strains) {
        let phase = phase_angular_interval_from_clock_dd(omega, start, 0.125).unwrap();
        for i in 0..3 {
            value[i] = value[i] + Dd::from_f64(strain[i].re) * phase.start[0]
                - Dd::from_f64(strain[i].im) * phase.start[1];
        }
    }
    value.map(Dd::to_f64)
}

// Seed actual original branch memory through a finite strain ramp, then use the
// original instantaneous mixed controller to establish zero transverse stress.
// Raw branch strains are tracked independently through those two known events.
fn seed(
    law: &TensorProny,
    plane: MaterialPlane,
    initial: [f64; 3],
) -> (Mixed, TensorPronyHistory, [f64; 20]) {
    let mut history = law.history();
    let preload = [0.002, -0.001, 0.0015, 0.0007, -0.0004, 0.0005];
    let duration = 0.17;
    history.advance(preload, duration, || false).unwrap();
    let previous = mandel(history.strain());
    let mut y = [0.0; 20];
    for (j, tau) in TAUS.into_iter().enumerate() {
        let r = duration / tau;
        let gain = -(-r).exp_m1() / r;
        for i in 0..6 {
            y[6 + 6 * j + i] = gain * previous[i];
        }
    }
    let controller = Mixed::new(
        &history,
        mask(plane),
        TensorStressControlBudget::default(),
        || false,
    )
    .unwrap();
    controller
        .advance(&mut history, prescribed(plane, initial), 0.0, || false)
        .unwrap();
    let next = mandel(history.strain());
    y[..6].copy_from_slice(&next);
    for j in 0..2 {
        for i in 0..6 {
            y[6 + 6 * j + i] += next[i] - previous[i];
        }
    }
    (controller, history, y)
}

// Independent partial-pivot elimination of C0_FF, written as six physical
// constraints. Production uses its original QR/SVD spectral realization.
fn solve(mut a: Tangent6, mut rhs: Voigt) -> Voigt {
    for j in 0..6 {
        let pivot = (j..6)
            .max_by(|&i, &k| a[i][j].abs().total_cmp(&a[k][j].abs()))
            .unwrap();
        a.swap(j, pivot);
        rhs.swap(j, pivot);
        let row = a[j];
        for i in j + 1..6 {
            let factor = a[i][j] / row[j];
            for k in j..6 {
                a[i][k] -= factor * row[k];
            }
            rhs[i] -= factor * rhs[j];
        }
    }
    let mut x = [0.0; 6];
    for i in (0..6).rev() {
        x[i] = (rhs[i] - (i + 1..6).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
struct OriginalOde {
    eq: Tangent6,
    branches: [Tangent6; 2],
    c0: Tangent6,
}
impl OriginalOde {
    fn stress(&self, y: &[f64; 20]) -> Voigt {
        let mut stress = apply(&self.eq, core::array::from_fn(|i| y[i]));
        for j in 0..2 {
            let arm = apply(
                &self.branches[j],
                core::array::from_fn(|i| y[6 + 6 * j + i]),
            );
            for i in 0..6 {
                stress[i] += arm[i];
            }
        }
        stress
    }
    fn branch_energy(&self, y: &[f64; 20]) -> f64 {
        (0..2)
            .map(|j| {
                let q = core::array::from_fn(|i| y[6 + 6 * j + i]);
                0.5 * dot(q, apply(&self.branches[j], q))
            })
            .sum()
    }
    fn stored(&self, y: &[f64; 20]) -> f64 {
        let e = core::array::from_fn(|i| y[i]);
        0.5 * dot(e, apply(&self.eq, e)) + self.branch_energy(y)
    }
    fn heat_rate(&self, y: &[f64; 20]) -> f64 {
        (0..2)
            .map(|j| {
                let q = core::array::from_fn(|i| y[6 + 6 * j + i]);
                dot(q, apply(&self.branches[j], q)) / TAUS[j]
            })
            .sum()
    }
    fn derivative(&self, y: &[f64; 20], plane: MaterialPlane, imposed: Voigt) -> [f64; 20] {
        let mask = mask(plane);
        let mut rhs = imposed;
        for (j, tau) in TAUS.into_iter().enumerate() {
            let stress = apply(
                &self.branches[j],
                core::array::from_fn(|i| y[6 + 6 * j + i]),
            );
            for i in 0..6 {
                if mask[i] {
                    rhs[i] += stress[i] / tau;
                }
            }
        }
        let mut a = self.c0;
        for i in 0..6 {
            if !mask[i] {
                a[i] = [0.0; 6];
                a[i][i] = 1.0;
            }
        }
        let rate = solve(a, rhs);
        let mut result = [0.0; 20];
        result[..6].copy_from_slice(&rate);
        for (j, tau) in TAUS.into_iter().enumerate() {
            for i in 0..6 {
                result[6 + 6 * j + i] = rate[i] - y[6 + 6 * j + i] / tau;
            }
        }
        result[18] = dot(self.stress(y), rate);
        result[19] = self.heat_rate(y);
        result
    }
    fn integrate(
        &self,
        mut y: [f64; 20],
        plane: MaterialPlane,
        start: f64,
        duration: f64,
        frequencies: &[f64],
        strains: &[[C64; 3]],
        bias: [f64; 3],
    ) -> [f64; 20] {
        y[18] = 0.0;
        y[19] = 0.0;
        let steps = 4096;
        let h = duration / steps as f64;
        let derivative = |state: &[f64; 20], t: f64| {
            let (value, rate) = waveform(frequencies, strains, bias, t);
            let imposed = mandel(prescribed(plane, rate));
            let mut state = *state;
            let e = mandel(prescribed(plane, value));
            for i in axes(plane).0 {
                state[i] = e[i];
            }
            self.derivative(&state, plane, imposed)
        };
        for step in 0..steps {
            let t = start + step as f64 * h;
            let k1 = derivative(&y, t);
            let k2 = derivative(
                &core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]),
                t + 0.5 * h,
            );
            let k3 = derivative(
                &core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]),
                t + 0.5 * h,
            );
            let k4 = derivative(&core::array::from_fn(|i| y[i] + h * k3[i]), t + h);
            for i in 0..20 {
                y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
            }
        }
        let value = mandel(prescribed(
            plane,
            waveform(frequencies, strains, bias, start + duration).0,
        ));
        for i in axes(plane).0 {
            y[i] = value[i];
        }
        y
    }
    fn compare_history(&self, history: &TensorPronyHistory, y: &[f64; 20]) {
        let stress = self.stress(y);
        for i in 0..6 {
            close(history.strain()[i], y[i] / WEIGHTS[i], 5e-13);
            close(history.stress()[i], stress[i] / WEIGHTS[i], 2e-10);
        }
        close(history.stored_energy(), self.stored(y), 2e-12);
        close(history.branch_stored_energy(), self.branch_energy(y), 2e-12);
    }
}

#[test]
fn retained_plane_stress_harmonics_match_original_full_tensor_ode_and_handoff() {
    let (law, ode) = coupled();
    for plane in [
        MaterialPlane::Axes12,
        MaterialPlane::Axes21,
        MaterialPlane::Axes13,
        MaterialPlane::Axes23,
    ] {
        let (start, duration) = (0.37, 0.29);
        let initial = waveform(&FREQUENCIES, &STRAINS, BIAS, start).0;
        let (control, history, y) = seed(&law, plane, initial);
        let original = history.clone();
        let interval = control
            .prepare_plane_stress_harmonic_window(
                plane,
                &FREQUENCIES,
                Dd::from_f64(start),
                duration,
                budget(),
                || false,
            )
            .unwrap();
        let result = interval
            .evaluate(&STRAINS, &BIAS, &history, || false)
            .unwrap();
        let reference = ode.integrate(y, plane, start, duration, &FREQUENCIES, &STRAINS, BIAS);
        ode.compare_history(&result.history, &reference);
        close(result.material.material_heat_j_m3, reference[19], 2e-12);
        close(result.material.stress_work_j_m3, reference[18], 2e-12);
        close(
            result.material.stored_energy_change_j_m3,
            ode.stored(&reference) - ode.stored(&y),
            2e-12,
        );
        assert!(result.material.material_heat_j_m3 > 0.0);
        close(result.material.energy_balance_residual_j_m3, 0.0, 2e-12);
        assert_eq!(history, original);
        for i in axes(plane).1 {
            assert_eq!(result.history.stress()[i], 0.0);
        }
        for i in 0..6 {
            assert!(result.strain_bounds[0][i] <= history.strain()[i]);
            assert!(result.strain_bounds[1][i] >= result.history.strain()[i]);
        }
        // Resume with the unmodified original mixed controller. Matching future
        // creep and heat tests the returned private memory, not just its endpoint.
        let held = waveform(&FREQUENCIES, &STRAINS, BIAS, start + duration).0;
        let mut resumed = result.history.clone();
        let report = control
            .advance(&mut resumed, prescribed(plane, held), 0.13, || false)
            .unwrap();
        let expected = ode.integrate(reference, plane, 0.0, 0.13, &[], &[], held);
        ode.compare_history(&resumed, &expected);
        close(report.dissipated, expected[19], 2e-12);
        close(report.work, 0.0, 2e-12);
    }
}

#[test]
fn held_memory_keeps_heat_below_rounded_endpoint_energy_resolution() {
    let (law, ode) = coupled();
    let plane = MaterialPlane::Axes12;
    let bias = [0.001, -0.0004, 0.0];
    let (control, history, y) = seed(&law, plane, bias);
    let h = 1e-20;
    let interval = control
        .prepare_plane_stress_harmonic_window(plane, &[], Dd::from_f64(2.0), h, budget(), || false)
        .unwrap();
    let result = interval.evaluate(&[], &bias, &history, || false).unwrap();
    let expected = ode.heat_rate(&y) * h;
    assert!(expected > 0.0);
    close(
        result.material.material_heat_j_m3,
        expected,
        5e-12 * expected,
    );
    close(
        result.material.stored_energy_change_j_m3,
        -expected,
        5e-12 * expected,
    );
    assert_eq!(result.material.stress_work_j_m3, 0.0);
    assert_eq!(
        result.material.start_stored_energy_j_m3,
        result.material.end_stored_energy_j_m3
    );
    assert!(result.strain_increment.iter().any(|v| *v != 0.0));
}

#[test]
fn dd_clock_partition_preserves_complete_constrained_history() {
    let (law, _) = coupled();
    let plane = MaterialPlane::Axes13;
    let start = Dd::from_pair(1e16, 0.3);
    let initial = at_clock(&FREQUENCIES, &STRAINS, BIAS, start);
    let (control, history, _) = seed(&law, plane, initial);
    let whole = control
        .prepare_plane_stress_harmonic_window(plane, &FREQUENCIES, start, 0.375, budget(), || false)
        .unwrap()
        .evaluate(&STRAINS, &BIAS, &history, || false)
        .unwrap();
    let first = control
        .prepare_plane_stress_harmonic_window(plane, &FREQUENCIES, start, 0.125, budget(), || false)
        .unwrap()
        .evaluate(&STRAINS, &BIAS, &history, || false)
        .unwrap();
    let second = control
        .prepare_plane_stress_harmonic_window(
            plane,
            &FREQUENCIES,
            start + Dd::from_f64(0.125),
            0.25,
            budget(),
            || false,
        )
        .unwrap()
        .evaluate(&STRAINS, &BIAS, &first.history, || false)
        .unwrap();
    for i in 0..6 {
        close(whole.history.strain()[i], second.history.strain()[i], 2e-14);
        close(whole.history.stress()[i], second.history.stress()[i], 2e-11);
    }
    close(
        whole.history.branch_stored_energy(),
        second.history.branch_stored_energy(),
        2e-13,
    );
    close(
        whole.material.material_heat_j_m3,
        first.material.material_heat_j_m3 + second.material.material_heat_j_m3,
        2e-13,
    );
    close(
        whole.material.stress_work_j_m3,
        first.material.stress_work_j_m3 + second.material.stress_work_j_m3,
        2e-13,
    );
    close(
        whole.material.stored_energy_change_j_m3,
        first.material.stored_energy_change_j_m3 + second.material.stored_energy_change_j_m3,
        2e-13,
    );
    close(whole.history.work(), second.history.work(), 2e-13);
    close(
        whole.history.dissipated(),
        second.history.dissipated(),
        2e-13,
    );
}

#[test]
fn equal_carriers_cancel_before_transverse_recovery_without_erasing_memory() {
    let (law, _) = coupled();
    let plane = MaterialPlane::Axes23;
    let (control, history, _) = seed(&law, plane, BIAS);
    let huge = [C64::new(1e100, -2e100); 3];
    let opposite = huge.map(|v| C64::new(-v.re, -v.im));
    let cancelled = control
        .prepare_plane_stress_harmonic_window(
            plane,
            &[3.0, 3.0],
            Dd::from_f64(0.2),
            0.4,
            budget(),
            || false,
        )
        .unwrap()
        .evaluate(&[huge, opposite], &BIAS, &history, || false)
        .unwrap();
    let held = control
        .prepare_plane_stress_harmonic_window(plane, &[], Dd::from_f64(0.2), 0.4, budget(), || {
            false
        })
        .unwrap()
        .evaluate(&[], &BIAS, &history, || false)
        .unwrap();
    assert_eq!(cancelled.history, held.history);
    assert_eq!(cancelled.strain_increment, held.strain_increment);
    assert_eq!(cancelled.strain_bounds, held.strain_bounds);
    let physical = |p: &fs_material::visco::TensorPronyHarmonicWindowPoint| {
        [
            p.start_stored_energy_j_m3,
            p.end_stored_energy_j_m3,
            p.stored_energy_change_j_m3,
            p.kelvin_heat_j_m3,
            p.maxwell_heat_j_m3,
            p.material_heat_j_m3,
            p.stress_work_j_m3,
        ]
    };
    assert_eq!(physical(&cancelled.material), physical(&held.material));
    // A zero-amplitude carrier and an empty bank may sum DD diagnostic words
    // differently while every reported physical account remains identical.
    for material in [&cancelled.material, &held.material] {
        close(
            material.energy_balance_residual_j_m3,
            0.0,
            64.0 * f64::EPSILON * material.material_heat_j_m3,
        );
    }
    assert!(cancelled.material.material_heat_j_m3 > 0.0);
    assert_eq!(cancelled.material.stress_work_j_m3, 0.0);
}

#[test]
fn startup_admits_actual_transverse_strain_envelope() {
    let eq = diagonal([1.0, 1.0, 0.01, 1.0, 1.0, 1.0]);
    let mut factor = [[0.0; 6]; 6];
    factor[0] = [50.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let law = TensorProny::new_mandel(
        eq,
        vec![TensorPronyBranch::from_gram_rows(factor, 0.2).unwrap()],
        0.05,
    )
    .unwrap();
    let history = law.history();
    let control = Mixed::new(
        &history,
        mask(MaterialPlane::Axes12),
        TensorStressControlBudget::default(),
        || false,
    )
    .unwrap();
    let interval = control
        .prepare_plane_stress_harmonic_window(
            MaterialPlane::Axes12,
            &[3.0],
            Dd::ZERO,
            0.6,
            budget(),
            || false,
        )
        .unwrap();
    let strain = [[C64::new(0.0, 0.002), C64::ZERO, C64::ZERO]];
    let error = interval
        .evaluate(&strain, &[0.0; 3], &history, || false)
        .unwrap_err();
    assert_eq!(
        error,
        TensorPronyError::Invalid(
            "transient plane-stress strain envelope exceeds the original material domain"
        )
    );
    assert_eq!(history.strain(), [0.0; 6]);
    assert_eq!(history.dissipated(), 0.0);
}

#[test]
fn constrained_harmonic_refusals_preserve_the_original_history() {
    let (law, _) = coupled();
    let plane = MaterialPlane::Axes12;
    let start = Dd::from_f64(0.37);
    let (control, history, _) = seed(&law, plane, at_clock(&FREQUENCIES, &STRAINS, BIAS, start));
    let unchanged = history.clone();
    let interval = control
        .prepare_plane_stress_harmonic_window(plane, &FREQUENCIES, start, 0.2, budget(), || false)
        .unwrap();
    let foreign = law.history();
    assert!(
        interval
            .evaluate(&STRAINS, &BIAS, &foreign, || false)
            .is_err()
    );
    assert!(
        interval
            .evaluate(&STRAINS, &[0.02, 0.0, 0.0], &history, || false)
            .is_err()
    );
    assert!(
        control
            .prepare_plane_stress_harmonic_window(
                MaterialPlane::Axes13,
                &FREQUENCIES,
                start,
                0.2,
                budget(),
                || false
            )
            .is_err()
    );
    for limits in [
        HarmonicBudget {
            maximum_payload_bytes: 0,
            ..budget()
        },
        HarmonicBudget {
            maximum_scalar_work: 0,
            ..budget()
        },
    ] {
        assert!(
            control
                .prepare_plane_stress_harmonic_window(
                    plane,
                    &FREQUENCIES,
                    start,
                    0.2,
                    limits,
                    || false
                )
                .is_err()
        );
    }
    assert!(
        control
            .prepare_plane_stress_harmonic_window(plane, &[101.0], start, 0.2, budget(), || false)
            .is_err()
    );
    assert!(matches!(
        control.prepare_plane_stress_harmonic_window(
            plane,
            &FREQUENCIES,
            start,
            0.2,
            budget(),
            || true
        ),
        Err(TensorPronyError::Cancelled)
    ));
    let mut polls = 0;
    assert_eq!(
        interval
            .evaluate(&STRAINS, &BIAS, &history, || {
                polls += 1;
                polls > 3
            })
            .unwrap_err(),
        TensorPronyError::Cancelled
    );
    assert_eq!(history, unchanged);
    // Fixed source DC/frequency admission also applies to the constrained path.
    let restricted = law.clone().with_band((1.0, 10.0)).unwrap();
    let (restricted_control, restricted_history, _) = seed(&restricted, plane, BIAS);
    let dc = restricted_control
        .prepare_plane_stress_harmonic_window(plane, &[], start, 0.2, budget(), || false)
        .unwrap();
    assert!(
        dc.evaluate(&[], &BIAS, &restricted_history, || false)
            .is_err()
    );
    let kelvin = law
        .with_kelvin_viscosity_mandel(diagonal([0.1; 6]))
        .unwrap()
        .history();
    assert!(
        Mixed::new(
            &kelvin,
            mask(plane),
            TensorStressControlBudget::default(),
            || false
        )
        .is_err()
    );
}
