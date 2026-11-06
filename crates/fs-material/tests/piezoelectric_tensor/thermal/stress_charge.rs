//! G0/G1/G3 synthetic checks of actual inverse-tensor thermal derivatives.
//! Closed scalar potentials and an independent dense solve are the oracles.
use super::*;

fn stress_source(
    c: C,
    e: E,
    p: P,
    cp: C,
    ep: E,
    pp: P,
    mut alter: impl FnMut(&mut PropertyClaim),
) -> MaterialCard {
    source(c, e, p, |claim| {
        let PropertyValue::Scalar { value, dims } = claim.value else {
            unreachable!()
        };
        let name = claim.key.name();
        let slope = if name == "density" {
            0.0
        } else {
            let suffix = name.rsplit_once('_').unwrap().1.as_bytes();
            let (i, j) = ((suffix[0] - b'1') as usize, (suffix[1] - b'1') as usize);
            if name.starts_with("piezoelectric_stiffness_e_mandel") {
                cp[i][j]
            } else if name.starts_with("piezoelectric_stress_charge_mandel") {
                ep[i][j]
            } else {
                pp[i][j]
            }
        };
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, value - 20.0 * slope), (340.0, value + 40.0 * slope)],
            dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claim.validity =
            ValidityDomain::unconstrained().with_quantity("T", temperature(), 280.0, 340.0);
        alter(claim);
    })
}
fn stress_tensor(card: &MaterialCard) -> TensorPiezoelectric {
    TensorPiezoelectric::resolve(
        card,
        &thermal_point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap()
}
fn bind(card: &MaterialCard, n: [f64; 3]) -> ThermalLongitudinalPiezoelectric {
    ThermalLongitudinalPiezoelectric::from_stress_charge_tensor(
        card,
        &stress_tensor(card),
        n,
        "T",
        [280.0, 340.0],
    )
    .unwrap()
}
fn scalar_source(cp33: f64, ep33: f64, pp33: f64) -> MaterialCard {
    let mut c = [[0.0; 6]; 6];
    let mut p = [[0.0; 3]; 3];
    for (i, row) in c.iter_mut().enumerate() {
        row[i] = 1e11;
    }
    for (i, row) in p.iter_mut().enumerate() {
        row[i] = 1.1e-8;
    }
    let mut e = [[0.0; 6]; 3];
    e[2][2] = 14.0;
    let mut cp = [[0.0; 6]; 6];
    let mut ep = [[0.0; 6]; 3];
    let mut pp = [[0.0; 3]; 3];
    cp[2][2] = cp33;
    ep[2][2] = ep33;
    pp[2][2] = pp33;
    stress_source(c, e, p, cp, ep, pp, |_| {})
}

#[test]
fn g1_affine_stress_charge_uses_exact_inverse_derivatives_and_scalar_potential() {
    let (cp, ep, pp) = (-2e8, 0.04, 2e-11);
    let card = scalar_source(cp, ep, pp);
    let law = bind(&card, [0.0, 0.0, 1.0]);
    let (strain, charge, volume) = (2e-5, 3e-4, 1e-5);
    let lower = law
        .minimum_heat_capacity_j_k(strain, charge, volume)
        .unwrap();
    for t in [280.0, 300.0, 317.0, 340.0] {
        let c = 1e11 + cp * (t - 300.0);
        let e = 14.0 + ep * (t - 300.0);
        let p = 1.1e-8 + pp * (t - 300.0);
        let at = law.at_temperature(t).unwrap();
        law.at_material(at.material()).unwrap();
        assert_eq!(at.material().material().properties().len(), 64);
        assert_eq!(
            at.material().material().card_identity(),
            card.content_hash()
        );
        near(
            at.material().compliance_at_fixed_field_pa_inv(),
            1.0 / c,
            1e-26,
            2e-13,
        );
        near(
            at.material().strain_charge_coefficient_m_v(),
            e / c,
            1e-25,
            2e-13,
        );
        near(
            at.material().permittivity_at_fixed_stress_f_m(),
            p + e * e / c,
            1e-22,
            2e-13,
        );
        let first = [
            -cp / (c * c),
            ep / c - e * cp / (c * c),
            pp + 2.0 * e * ep / c - e * e * cp / (c * c),
        ];
        let second = [
            2.0 * cp * cp / (c * c * c),
            2.0 * e * cp * cp / (c * c * c) - 2.0 * ep * cp / (c * c),
            2.0 * (ep - e * cp / c).powi(2) / c,
        ];
        for (actual, expected) in at
            .coefficient_first_derivatives_per_k()
            .into_iter()
            .zip(first)
        {
            near(actual, expected, 1e-28, 3e-12);
        }
        for (actual, expected) in at
            .coefficient_second_derivatives_per_k2()
            .unwrap()
            .into_iter()
            .zip(second)
        {
            near(actual, expected, 1e-30, 3e-12);
        }
        // Independent fixed-strain/charge Helmholtz potential:
        // psi=V*(C strain^2+(D-e strain)^2/epsilonS)/2.
        let z = charge - e * strain;
        let zp = -ep * strain;
        let free = 0.5 * volume * (c * strain * strain + z * z / p);
        let entropy =
            -volume * (0.5 * cp * strain * strain + z * zp / p - 0.5 * z * z * pp / (p * p));
        let capacity = -t * volume / p * (zp - z * pp / p).powi(2);
        let state = at.state(strain, charge, volume).unwrap();
        near(state.free_energy_j, free, 1e-18, 3e-12);
        near(state.entropy_j_k, entropy, 1e-21, 3e-12);
        near(state.internal_energy_j, free + t * entropy, 1e-18, 3e-12);
        near(state.heat_capacity_j_k, capacity, 1e-21, 3e-12);
        assert!(state.heat_capacity_j_k >= lower);
    }
    assert_ne!(
        law.at_temperature(280.0)
            .unwrap()
            .coefficient_first_derivatives_per_k(),
        law.at_temperature(340.0)
            .unwrap()
            .coefficient_first_derivatives_per_k()
    );
}

#[test]
fn g1_affine_elastic_stiffness_has_zero_field_capacity_despite_inverse_curvature() {
    let card = scalar_source(2e8, 0.0, 0.0);
    let law = bind(&card, [0.0, 0.0, 1.0]);
    for t in [280.0, 300.0, 340.0] {
        let at = law.at_temperature(t).unwrap();
        assert!(at.coefficient_second_derivatives_per_k2().unwrap()[0] > 0.0);
        for (strain, charge) in [(2e-5, 3e-4), (-2e-5, 3e-4), (0.0, 3e-4)] {
            let state = at.state(strain, charge, 1e-5).unwrap();
            near(
                state.entropy_j_k,
                -0.5 * 1e-5 * 2e8 * strain * strain,
                1e-20,
                3e-12,
            );
            // Omitting A'' gives a false negative capacity here. Exact
            // cancellation of the two finite terms must also remain admissible.
            near(state.heat_capacity_j_k, 0.0, 1e-18, 0.0);
        }
    }
}

#[test]
fn g1_transverse_thermal_relaxation_can_have_positive_field_heat_capacity() {
    let mut c = [[0.0; 6]; 6];
    let mut cp = [[0.0; 6]; 6];
    let mut p = [[0.0; 3]; 3];
    for (i, row) in c.iter_mut().enumerate() {
        row[i] = 1e11;
    }
    for (i, row) in p.iter_mut().enumerate() {
        row[i] = 1e-8;
    }
    let a = 0.003;
    cp[0][0] = 1e11 * a;
    cp[1][1] = -1e11 * a;
    let card = stress_source(
        c,
        [[0.0; 6]; 3],
        p,
        cp,
        [[0.0; 6]; 3],
        [[0.0; 3]; 3],
        |_| {},
    );
    let n = [
        core::f64::consts::FRAC_1_SQRT_2,
        core::f64::consts::FRAC_1_SQRT_2,
        0.0,
    ];
    let law = bind(&card, n);
    let (strain, volume) = (1e-5, 1e-5);
    let at = law.at_temperature(300.0).unwrap();
    // s(T)=(1/2+1/[2(1-a^2(T-300)^2)])/C0;
    // s'=0 and s''=a^2/C0 at the reference temperature.
    near(at.coefficient_first_derivatives_per_k()[0], 0.0, 1e-28, 0.0);
    near(
        at.coefficient_second_derivatives_per_k2().unwrap()[0],
        a * a / 1e11,
        1e-30,
        3e-12,
    );
    let state = at.state(strain, 0.0, volume).unwrap();
    let expected = 0.5 * 300.0 * volume * 1e11 * a * a * strain * strain;
    assert!(state.heat_capacity_j_k > 0.0);
    near(state.heat_capacity_j_k, expected, 1e-20, 3e-12);
    let h = 0.1;
    let free = |t| {
        law.at_temperature(t)
            .unwrap()
            .state(strain, 0.0, volume)
            .unwrap()
            .free_energy_j
    };
    near(
        state.heat_capacity_j_k,
        -300.0 * (free(300.0 + h) - 2.0 * free(300.0) + free(300.0 - h)) / (h * h),
        1e-16,
        3e-6,
    );
    assert!(state.heat_capacity_j_k >= law.minimum_heat_capacity_j_k(strain, 0.0, volume).unwrap());
}

// Independent partial-pivot elimination, never the material's Cholesky owner.
fn solve(mut a: C, mut b: [f64; 6]) -> [f64; 6] {
    for j in 0..6 {
        let pivot = (j..6)
            .max_by(|&i, &k| a[i][j].abs().total_cmp(&a[k][j].abs()))
            .unwrap();
        a.swap(j, pivot);
        b.swap(j, pivot);
        for i in j + 1..6 {
            let ratio = a[i][j] / a[j][j];
            for k in j..6 {
                let pivot = a[j][k];
                a[i][k] -= ratio * pivot;
            }
            let pivot = b[j];
            b[i] -= ratio * pivot;
        }
    }
    let mut x = [0.0; 6];
    for i in (0..6).rev() {
        x[i] = (b[i] - (i + 1..6).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
fn dense_source() -> (C, E, P, C, E, P) {
    let (c, e, p) = coefficients();
    let mut cp: C =
        core::array::from_fn(|i| core::array::from_fn(|j| 0.0002 * (i + j + 1) as f64 * c[i][j]));
    cp[0][4] = 2e7;
    cp[4][0] = cp[0][4];
    cp[2][5] = -4e7;
    cp[5][2] = cp[2][5];
    let ep = core::array::from_fn(|i| {
        core::array::from_fn(|j| 0.001 * (2 + i) as f64 * e[i][j] + 0.003 * (j as f64 - i as f64))
    });
    let pp =
        core::array::from_fn(|i| core::array::from_fn(|j| 0.0007 * (i + j + 1) as f64 * p[i][j]));
    (c, e, p, cp, ep, pp)
}
fn independent_coefficients(t: f64, n: [f64; 3]) -> [f64; 3] {
    let (c, e, p, cp, ep, pp) = dense_source();
    let c = core::array::from_fn(|i| core::array::from_fn(|j| c[i][j] + (t - 300.0) * cp[i][j]));
    let e: E = core::array::from_fn(|i| core::array::from_fn(|j| e[i][j] + (t - 300.0) * ep[i][j]));
    let p: P = core::array::from_fn(|i| core::array::from_fn(|j| p[i][j] + (t - 300.0) * pp[i][j]));
    let r = [
        n[0] * n[0],
        n[1] * n[1],
        n[2] * n[2],
        2.0_f64.sqrt() * n[0] * n[1],
        2.0_f64.sqrt() * n[1] * n[2],
        2.0_f64.sqrt() * n[2] * n[0],
    ];
    let h = core::array::from_fn(|j| (0..3).map(|i| n[i] * e[i][j]).sum());
    let v = solve(c, r);
    let w = solve(c, h);
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    [
        dot(&r, &v),
        dot(&h, &v),
        dot(&n, &apply(&p, &n)) + dot(&h, &w),
    ]
}
fn independent_free(t: f64, n: [f64; 3], strain: f64, charge: f64, volume: f64) -> f64 {
    let [s, d, p] = independent_coefficients(t, n);
    0.5 * volume * (p * strain * strain - 2.0 * d * strain * charge + s * charge * charge)
        / (s * p - d * d)
}

#[test]
fn g1_g3_full_stress_charge_thermal_path_matches_independent_tensor_inverse_and_corotation() {
    let (c, e, p, cp, ep, pp) = dense_source();
    let card = stress_source(c, e, p, cp, ep, pp, |_| {});
    let tensor = stress_tensor(&card);
    let n = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    let rotated = tensor.rotated(ContentHash([9; 32]), q).unwrap();
    let base = bind(&card, n);
    let moved = ThermalLongitudinalPiezoelectric::from_stress_charge_tensor(
        &card,
        &rotated,
        apply(&q, &n),
        "T",
        [280.0, 340.0],
    )
    .unwrap();
    let reversed = bind(&card, n.map(|v| -v));
    let (strain, charge, volume) = (2e-5, 3e-4, 1e-5);
    let lower = base
        .minimum_heat_capacity_j_k(strain, charge, volume)
        .unwrap();
    for t in [280.0, 295.0, 300.0, 320.0, 340.0] {
        let at = base.at_temperature(t).unwrap();
        let coefficients = independent_coefficients(t, n);
        for (actual, expected) in [
            at.material().compliance_at_fixed_field_pa_inv(),
            at.material().strain_charge_coefficient_m_v(),
            at.material().permittivity_at_fixed_stress_f_m(),
        ]
        .into_iter()
        .zip(coefficients)
        {
            near(actual, expected, 1e-25, 3e-12);
        }
        base.at_material(at.material()).unwrap();
        let b = moved.at_temperature(t).unwrap();
        let r = reversed.at_temperature(t).unwrap();
        moved.at_material(b.material()).unwrap();
        assert!(base.at_material(b.material()).is_err());
        assert!(base.at_material(r.material()).is_err());
        for (s, d) in [
            (strain, charge),
            (-strain, charge),
            (strain, -charge),
            (-strain, -charge),
        ] {
            let state = at.state(s, d, volume).unwrap();
            near(
                state.free_energy_j,
                independent_free(t, n, s, d, volume),
                1e-18,
                3e-12,
            );
            assert!(state.heat_capacity_j_k >= lower);
            for other in [
                b.state(s, d, volume).unwrap(),
                r.state(s, -d, volume).unwrap(),
            ] {
                near(state.free_energy_j, other.free_energy_j, 1e-18, 3e-12);
                near(state.entropy_j_k, other.entropy_j_k, 1e-21, 3e-12);
                near(
                    state.heat_capacity_j_k,
                    other.heat_capacity_j_k,
                    1e-21,
                    3e-12,
                );
            }
        }
    }
    let t = 315.0;
    let h = 0.05;
    let at = base.at_temperature(t).unwrap();
    let state = at.state(strain, charge, volume).unwrap();
    let free = |t| independent_free(t, n, strain, charge, volume);
    near(
        state.entropy_j_k,
        -(free(t + h) - free(t - h)) / (2.0 * h),
        1e-20,
        3e-7,
    );
    near(
        state.heat_capacity_j_k,
        -t * (free(t + h) - 2.0 * free(t) + free(t - h)) / (h * h),
        1e-16,
        3e-6,
    );
    near(
        state.internal_energy_j,
        free(t) + t * state.entropy_j_k,
        1e-18,
        3e-12,
    );
    let minus = independent_coefficients(t - h, n);
    let zero = independent_coefficients(t, n);
    let plus = independent_coefficients(t + h, n);
    for (i, actual) in at
        .coefficient_first_derivatives_per_k()
        .into_iter()
        .enumerate()
    {
        near(actual, (plus[i] - minus[i]) / (2.0 * h), 1e-27, 3e-7);
    }
    for (i, actual) in at
        .coefficient_second_derivatives_per_k2()
        .unwrap()
        .into_iter()
        .enumerate()
    {
        near(
            actual,
            (plus[i] - 2.0 * zero[i] + minus[i]) / (h * h),
            1e-25,
            3e-5,
        );
    }
}

#[test]
fn g0_stress_charge_thermal_refuses_nonaffine_sources_and_unexcited_instability() {
    let (c, e, p, cp, ep, pp) = dense_source();
    let card = stress_source(c, e, p, cp, ep, pp, |_| {});
    let law = bind(&card, [0.0, 0.0, 1.0]);
    assert!(law.at_temperature(279.0).is_err());
    assert!(law.at_temperature(f64::NAN).is_err());
    assert!(law.minimum_heat_capacity_j_k(-1.0, 1.0, 1.0).is_err());
    let nonlinear = stress_source(c, e, p, cp, ep, pp, |claim| {
        if claim.key.name() == "piezoelectric_stiffness_e_mandel_44" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                knots.insert(1, (300.0, 5.1e10));
            }
        }
    });
    assert!(
        ThermalLongitudinalPiezoelectric::from_stress_charge_tensor(
            &nonlinear,
            &stress_tensor(&nonlinear),
            [0.0, 0.0, 1.0],
            "T",
            [280.0, 340.0],
        )
        .is_err()
    );
    let unstable = stress_source(c, e, p, cp, ep, pp, |claim| {
        if [
            "piezoelectric_stiffness_e_mandel_45",
            "piezoelectric_stiffness_e_mandel_54",
        ]
        .contains(&claim.key.name())
        {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                *knots = vec![(280.0, -1.2e11), (340.0, 2.4e11)];
            }
        }
    });
    assert!(
        ThermalLongitudinalPiezoelectric::from_stress_charge_tensor(
            &unstable,
            &stress_tensor(&unstable),
            [0.0, 0.0, 1.0],
            "T",
            [280.0, 340.0],
        )
        .is_err()
    );
    let density = stress_source(c, e, p, cp, ep, pp, |claim| {
        if claim.key.name() == "density" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                knots[1].1 += 100.0;
            }
        }
    });
    assert!(
        ThermalLongitudinalPiezoelectric::from_stress_charge_tensor(
            &density,
            &stress_tensor(&density),
            [0.0, 0.0, 1.0],
            "T",
            [280.0, 340.0],
        )
        .is_err()
    );
    let wrong = thermal_source(|_| {});
    assert!(
        ThermalLongitudinalPiezoelectric::from_stress_charge_tensor(
            &wrong,
            &tensor(&wrong),
            [0.0, 0.0, 1.0],
            "T",
            [280.0, 340.0],
        )
        .is_err()
    );
}
