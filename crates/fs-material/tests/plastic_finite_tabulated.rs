//! G0/G1/G3: independent scalar consistency, work, finite tensor derivatives and
//! loading histories for explicitly tabulated Kirchhoff hardening. Synthetic
//! data test the numerical model, not the calibration of a particular metal.

use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{
    FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2Response,
};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const KNOTS: [(f64, f64); 5] = [(0.0, Y), (0.04, 5e6), (0.12, 6e6), (0.3, 6e6), (0.8, 20e6)];

fn table() -> TabulatedHardening {
    TabulatedHardening::new(
        KNOTS
            .iter()
            .map(|&(plastic_strain, yield_stress_pa)| HardeningPoint {
                plastic_strain,
                yield_stress_pa,
            })
            .collect(),
    )
    .unwrap()
}

fn base(hardening: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: hardening,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 0.8,
        },
    )
    .unwrap()
}
fn law() -> FiniteJ2 {
    base(0.0).with_tabulated_hardening(table()).unwrap()
}
fn kinetics(exponent: f64) -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: 16e6,
        stress_scale_pa: 2e6,
        exponent,
    }
}
fn coaxial(axial: f64, mean: f64) -> [f64; 9] {
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
fn radius(alpha: f64) -> f64 {
    for pair in KNOTS.windows(2) {
        let [(a, ya), (b, yb)] = [pair[0], pair[1]];
        if alpha <= b {
            return ya + (yb - ya) * ((alpha - a) / (b - a));
        }
    }
    panic!("oracle outside measured support: {alpha}")
}
fn storage(alpha: f64) -> f64 {
    let mut total = 0.0;
    for pair in KNOTS.windows(2) {
        let [(a, ya), (b, yb)] = [pair[0], pair[1]];
        if alpha <= a {
            break;
        }
        let dx = alpha.min(b) - a;
        total += (ya - Y) * dx + 0.5 * ((yb - ya) / (b - a)) * dx * dx;
    }
    total
}
fn bisect(mut a: f64, mut b: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(a) <= 0.0 && f(b) >= 0.0);
    for _ in 0..160 {
        let x = 0.5 * a + 0.5 * b;
        if f(x) < 0.0 {
            a = x;
        } else {
            b = x;
        }
    }
    0.5 * a + 0.5 * b
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "actual={a:.16e}, expected={b:.16e}, difference={:.8e}",
        (a - b).abs()
    );
}
fn derivatives(r: &FiniteJ2Response, f: &[f64; 9], update: impl Fn(&[f64; 9]) -> ([f64; 9], f64)) {
    for column in 0..9 {
        let step = 1e-6;
        let mut plus = *f;
        let mut minus = *f;
        plus[column] += step;
        minus[column] -= step;
        let (pp, wp) = update(&plus);
        let (pm, wm) = update(&minus);
        near(
            (wp - wm) / (2.0 * step),
            r.first_piola_pa[column],
            0.3,
            5e-6,
        );
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * step),
                r.tangent_pa[row][column],
                3.0,
                5e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                3e-12,
            );
        }
    }
}

#[test]
fn tabulated_finite_compression_crosses_knots_and_preserves_work_partition() {
    let material = law();
    for axial in [-0.025_f64, -0.09, -0.23, -0.5] {
        let mean = -0.02;
        let trial = 3.0 * MU * axial.abs();
        let gamma = bisect(0.0, axial.abs(), |g| 3.0 * MU * g + radius(g) - trial);
        let r = material
            .response(&coaxial(axial, mean), &material.initial_state())
            .unwrap();
        near(r.state.alpha, gamma, 2e-14, 3e-13);
        near(r.equivalent_kirchhoff_pa, radius(gamma), 3e-6, 3e-13);
        near(
            r.state.inverse_plastic_cauchy_green[0],
            (2.0 * gamma).exp(),
            5e-14,
            3e-13,
        );
        let stored = 0.5 * K * (3.0 * mean).powi(2)
            + 1.5 * MU * (axial.abs() - gamma).powi(2)
            + storage(gamma);
        near(r.stored_energy_j_m3, stored, 3e-7, 3e-13);
        near(r.plastic_dissipation_j_m3, Y * gamma, 3e-7, 3e-13);
        near(
            r.incremental_potential_j_m3,
            stored + Y * gamma,
            3e-7,
            3e-13,
        );
        near(r.jacobian, (3.0 * mean).exp(), 3e-14, 3e-13);
        material
            .validate_committed(&coaxial(axial, mean), &r.state, 1e-4)
            .unwrap();
    }
}

#[test]
fn tabulated_finite_unloading_and_reverse_flow_keep_the_accumulated_history() {
    let material = law();
    let mut state = material.initial_state();
    let mut plastic = 0.0;
    let mut alpha = 0.0;
    for axial in [0.1, 0.15, 0.14, 0.0, -0.1, -0.16, 0.0, 0.1] {
        let trial: f64 = 3.0 * MU * (axial - plastic);
        let increment = if trial.abs() <= radius(alpha) {
            0.0
        } else {
            bisect(0.0, 0.8 - alpha, |g| {
                3.0 * MU * g + radius(alpha + g) - trial.abs()
            })
        };
        alpha += increment;
        plastic += trial.signum() * increment;
        let f = coaxial(axial, 0.0);
        let r = material.response(&f, &state).unwrap();
        near(r.state.alpha, alpha, 3e-14, 5e-13);
        near(r.kirchhoff_pa[0], 2.0 * MU * (axial - plastic), 4e-6, 5e-13);
        near(
            r.state.inverse_plastic_cauchy_green[0],
            (-2.0 * plastic).exp(),
            5e-14,
            5e-13,
        );
        let held = material.response(&f, &r.state).unwrap();
        assert_eq!(held.state, r.state);
        assert_eq!(held.plastic_dissipation_j_m3, 0.0);
        state = r.state;
    }
}

#[test]
fn tabulated_full_tensor_tangent_and_potential_survive_noncoaxial_history_and_rotation() {
    let material = law();
    let old = material
        .response(&coaxial(0.07, 0.0), &material.initial_state())
        .unwrap()
        .state;
    let f = [1.3, 0.15, 0.03, 0.02, 0.9, 0.08, 0.01, -0.03, 0.93];
    let r = material.response(&f, &old).unwrap();
    derivatives(&r, &f, |f| material.stress_and_potential(f, &old).unwrap());
    let (s, c) = 0.67_f64.sin_cos();
    let rot = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let multiply = |a: &[f64; 9], b: &[f64; 9]| -> [f64; 9] {
        core::array::from_fn(|ij| {
            (0..3)
                .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
                .sum()
        })
    };
    let turned = material.response(&multiply(&rot, &f), &old).unwrap();
    for (a, b) in turned
        .first_piola_pa
        .into_iter()
        .zip(multiply(&rot, &r.first_piola_pa))
    {
        near(a, b, 5e-6, 5e-12);
    }
    near(turned.state.alpha, r.state.alpha, 3e-14, 3e-12);
    near(turned.stored_energy_j_m3, r.stored_energy_j_m3, 3e-7, 3e-12);
    for i in 0..6 {
        near(
            turned.state.inverse_plastic_cauchy_green[i],
            r.state.inverse_plastic_cauchy_green[i],
            5e-14,
            3e-12,
        );
    }
}

#[test]
fn tabulated_perzyna_matches_independent_rate_equation_and_full_tangent() {
    for n in [1.0, 3.2] {
        let material = FiniteOverstressJ2::new(law(), kinetics(n), 0.3).unwrap();
        let old = material.initial_state();
        let f = coaxial(0.23, 0.015);
        let viscous = |g: f64| 2e6 * (16e6 * g / (0.3 * 2e6)).powf(1.0 / n);
        let gamma = bisect(0.0, 0.23, |g| {
            3.0 * MU * g + radius(g) + viscous(g) - 3.0 * MU * 0.23
        });
        let r = material.response(&f, &old).unwrap();
        near(r.state.alpha, gamma, 3e-14, 5e-12);
        near(
            r.equivalent_kirchhoff_pa,
            radius(gamma) + viscous(gamma),
            4e-6,
            5e-12,
        );
        near(
            r.plastic_dissipation_j_m3,
            (Y + viscous(gamma)) * gamma,
            3e-7,
            5e-12,
        );
        near(
            r.flow_potential_j_m3,
            (Y + n / (n + 1.0) * viscous(gamma)) * gamma,
            3e-7,
            5e-12,
        );
        derivatives(&r, &f, |f| material.stress_and_potential(f, &old).unwrap());
    }
}

#[test]
fn tabulated_joint_creep_retains_distinct_plastic_and_creep_accumulations() {
    for n in [1.0, 2.5] {
        let material = FiniteNortonJ2::new(law(), kinetics(n), 0.1, 2.0).unwrap();
        let old = material.initial_state();
        let f = coaxial(0.3, 0.01);
        let creep_at = |q: f64| 0.1 * (2e6 / 16e6) * (q / 2e6).powf(n);
        let plastic = bisect(0.0, 0.3, |p| {
            radius(p) + 3.0 * MU * (p + creep_at(radius(p))) - 3.0 * MU * 0.3
        });
        let creep = creep_at(radius(plastic));
        let r = material.response(&f, &old).unwrap();
        near(r.plastic_increment, plastic, 3e-14, 5e-12);
        near(r.creep_increment, creep, 3e-14, 5e-12);
        near(r.equivalent_kirchhoff_pa, radius(plastic), 4e-6, 5e-12);
        near(r.plastic_dissipation_j_m3, Y * plastic, 3e-7, 5e-12);
        near(
            r.creep_dissipation_j_m3,
            radius(plastic) * creep,
            3e-7,
            5e-12,
        );
        derivatives(&r, &f, |f| material.stress_and_potential(f, &old).unwrap());
    }
}

#[test]
fn tabulated_affine_curves_match_native_linear_law_and_refuse_extrapolation() {
    for h in [0.0, 8e6] {
        let curve = TabulatedHardening::new(vec![
            HardeningPoint {
                plastic_strain: 0.0,
                yield_stress_pa: Y,
            },
            HardeningPoint {
                plastic_strain: 0.8,
                yield_stress_pa: Y + 0.8 * h,
            },
        ])
        .unwrap();
        let original = base(h);
        let table = base(0.0).with_tabulated_hardening(curve).unwrap();
        let f = coaxial(-0.15, 0.025);
        let a = original.response(&f, &original.initial_state()).unwrap();
        let b = table.response(&f, &table.initial_state()).unwrap();
        near(a.state.alpha, b.state.alpha, 3e-15, 3e-13);
        near(a.stored_energy_j_m3, b.stored_energy_j_m3, 3e-7, 3e-13);
        for i in 0..9 {
            near(a.first_piola_pa[i], b.first_piola_pa[i], 3e-6, 3e-13);
            for j in 0..9 {
                near(a.tangent_pa[i][j], b.tangent_pa[i][j], 3e-5, 3e-13);
            }
        }
    }
    let material = law();
    let old = material.initial_state();
    assert!(material.response(&coaxial(1.1, 0.0), &old).is_err());
    assert_eq!(old, material.initial_state());
    assert!(base(1.0).with_tabulated_hardening(table()).is_err());
    assert!(material.clone().with_voce_hardening(1.0, 1.0).is_err());
    let rate = FiniteOverstressJ2::new(material, kinetics(1.0), 100.0).unwrap();
    assert!(rate.response(&coaxial(1.1, 0.0), &old).is_err());
}

#[test]
fn tabulated_primary_creep_brackets_without_extrapolating_and_keeps_its_native_potential() {
    use fs_material::creep::StrainHardeningCreep;
    use fs_material::plastic::finite::FiniteJ2State;

    let (eta, dt, reference, mean) = (40e6, 0.2, 0.015, 0.015);
    let (old_plastic, old_creep, increment) = (0.02_f64, 0.035_f64, 0.025_f64);
    let old_total = old_plastic + old_creep;
    let old = FiniteJ2State {
        inverse_plastic_cauchy_green: [
            (-2.0 * old_total).exp(),
            old_total.exp(),
            old_total.exp(),
            0.0,
            0.0,
            0.0,
        ],
        alpha: old_plastic,
        creep_alpha: old_creep,
    };
    // Only a narrow measured interval remains. The initial primary-creep
    // bracket would imply alpha > 0.046; the actual coupled root is 0.045.
    // Refusing an out-of-domain bracket sample would reject this valid step.
    let native = base(0.0);
    let bounded = FiniteJ2::new(
        *native.parameters(),
        FiniteJ2Limits {
            maximum_accumulated_plastic_strain: 0.046,
            ..*native.limits()
        },
    )
    .unwrap()
    .with_tabulated_hardening(
        TabulatedHardening::new(
            [0.0, 0.04, 0.046]
                .into_iter()
                .map(|plastic_strain| HardeningPoint {
                    plastic_strain,
                    yield_stress_pa: radius(plastic_strain),
                })
                .collect(),
        )
        .unwrap(),
    )
    .unwrap();
    let material = FiniteNortonJ2::new(
        bounded,
        FiniteOverstressParameters {
            viscosity_pa_s: eta,
            stress_scale_pa: 2e6,
            exponent: 1.0,
        },
        dt,
        2.0,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: reference,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap();
    material.validate_history(&old).unwrap();

    // Independently prescribed plastic endpoint. The primary multiplier
    // ref/(ref+old_creep+dc) gives a quadratic scalar equation for creep.
    let q = radius(old_plastic + increment);
    let shifted = reference + old_creep;
    let creep_at = |stress: f64| {
        let drive = reference * dt * stress / eta;
        2.0 * drive / (shifted + (shifted * shifted + 4.0 * drive).sqrt())
    };
    let creep = creep_at(q);
    let f = coaxial(old_total + increment + creep + q / (3.0 * MU), mean);
    let response = material.response(&f, &old).unwrap();
    near(response.state.alpha, old_plastic + increment, 2e-13, 2e-11);
    near(response.state.creep_alpha, old_creep + creep, 2e-13, 2e-11);
    near(response.equivalent_kirchhoff_pa, q, 3e-5, 2e-11);
    near(
        response.plastic_dissipation_j_m3,
        Y * increment,
        2e-5,
        2e-11,
    );
    near(response.creep_dissipation_j_m3, q * creep, 2e-5, 2e-11);
    let creep_potential = eta / dt
        * (0.5 * (1.0 + old_creep / reference) * creep * creep + creep.powi(3) / (3.0 * reference));
    let stored =
        0.5 * K * (3.0 * mean).powi(2) + q * q / (6.0 * MU) + storage(old_plastic + increment);
    near(
        response.flow_potential_j_m3,
        Y * increment + creep_potential,
        2e-5,
        2e-10,
    );
    near(
        response.incremental_potential_j_m3,
        stored + Y * increment + creep_potential,
        2e-5,
        2e-10,
    );
    assert!(creep_potential < 0.5 * q * creep);
    derivatives(&response, &f, |f| {
        material.stress_and_potential(f, &old).unwrap()
    });
    material
        .validate_committed(&f, &response.state, 1e-4)
        .unwrap();

    // A root that actually requires extrapolation still refuses without
    // changing the borrowed plastic metric or either accumulated mechanism.
    let outside_plastic = 0.048;
    let outside_q = radius(outside_plastic);
    let outside_f = coaxial(
        old_total + (outside_plastic - old_plastic) + creep_at(outside_q) + outside_q / (3.0 * MU),
        mean,
    );
    let before = old.clone();
    assert!(material.response(&outside_f, &old).is_err());
    assert_eq!(old, before);
}
