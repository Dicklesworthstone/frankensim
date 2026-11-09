//! Original continuous CTE tables through the same 3-D potential/exchange.
//! Manufactured coefficients and independent physical-tensor references.
use super::*;

const TABLES: [[(f64, f64); 4]; 3] = [
    [
        (280.0, 8e-5),
        (310.0, 1e-4),
        (330.0, 1.8e-4),
        (360.0, 1.4e-4),
    ],
    [
        (280.0, -3e-5),
        (310.0, -8e-5),
        (330.0, -2e-5),
        (360.0, 4e-5),
    ],
    [(280.0, 2e-4), (310.0, 1e-4), (330.0, 5e-5), (360.0, 3e-4)],
];

fn tabulated_card(varying_elasticity: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, knots: Vec<(f64, f64)>| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored piecewise CTE tensor fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let affine = |base: f64, slope: f64| {
        let slope = if varying_elasticity { slope } else { 0.0 };
        vec![(280.0, base - 20.0 * slope), (360.0, base + 60.0 * slope)]
    };
    for i in 0..3 {
        insert(
            ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i],
            Pressure::DIMS,
            affine(E[i], ET[i]),
        );
        insert(
            ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i],
            Dims::NONE,
            affine(NU[i], NUT[i]),
        );
        insert(
            ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i],
            Pressure::DIMS,
            affine(G[i], GT[i]),
        );
        insert(
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i],
            Dims([0, 0, 0, -1, 0, 0]),
            TABLES[i].to_vec(),
        );
    }
    insert(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, affine(E[0], ET[0]));
    insert(POISSON_RATIO_PROPERTY, Dims::NONE, affine(NU[0], NUT[0]));
    insert(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        Dims([0, 0, 0, -1, 0, 0]),
        TABLES[0].to_vec(),
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic table expansion".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn prepared(card: &MaterialCard, isotropic: bool, range: [f64; 2]) -> AffineThermoelastic {
    let (elastic, expansion) = sources(card, 300.0, isotropic, isotropic);
    let stiffness = if isotropic {
        ThermoelasticSymmetry::Isotropic
    } else {
        ThermoelasticSymmetry::OrientedOrthotropic {
            principal_to_target: QS,
        }
    };
    let expansion_frame = if isotropic {
        ThermoelasticSymmetry::Isotropic
    } else {
        ThermoelasticSymmetry::OrientedOrthotropic {
            principal_to_target: QA,
        }
    };
    AffineThermoelastic::new(
        card,
        &elastic,
        &expansion,
        stiffness,
        expansion_frame,
        range,
        0.1,
    )
    .unwrap()
}

fn integral(axis: usize, t: f64) -> f64 {
    let (lo, hi, sign) = if t >= 300.0 {
        (300.0, t, 1.0)
    } else {
        (t, 300.0, -1.0)
    };
    let mut sum = 0.0;
    for pair in TABLES[axis].windows(2) {
        let [(a, x), (b, y)] = [pair[0], pair[1]];
        let start = lo.max(a);
        let end = hi.min(b);
        if end > start {
            let slope = (y - x) / (b - a);
            let width = end - start;
            sum += (x + slope * (start - a)) * width + 0.5 * slope * width * width;
        }
    }
    sign * sum
}

// Reuse the independent compliance-elimination oracle of the parent fixture,
// shifting only its analytic free strain to this table's exact integral.
fn expected(t: f64, total: [f64; 6], isotropic: bool, volume: f64) -> (f64, [f64; 6]) {
    if isotropic {
        let e = E[0] + ET[0] * (t - 300.0);
        let nu = NU[0] + NUT[0] * (t - 300.0);
        let bulk = e / (3.0 * (1.0 - 2.0 * nu));
        let shear = e / (2.0 * (1.0 + nu));
        let mut elastic = total;
        for x in &mut elastic[..3] {
            *x -= integral(0, t);
        }
        let trace = elastic[..3].iter().sum::<f64>();
        let mut deviator = elastic;
        for x in &mut deviator[..3] {
            *x -= trace / 3.0;
        }
        let stress = core::array::from_fn(|i| {
            2.0 * shear * deviator[i] + if i < 3 { bulk * trace } else { 0.0 }
        });
        let psi = volume
            * (0.5 * bulk * trace * trace + shear * deviator.iter().map(|x| x * x).sum::<f64>());
        return (psi, stress);
    }
    let dt = t - 300.0;
    let mut correction = [[0.0; 3]; 3];
    for i in 0..3 {
        correction[i][i] = ALPHA[i] * dt + 0.5 * ALPHAT[i] * dt * dt - integral(i, t);
    }
    let correction = mandel(rotate(correction, QA));
    oracle(
        t,
        core::array::from_fn(|i| total[i] + correction[i]),
        QS,
        QA,
        volume,
    )
}

#[test]
fn tabulated_isotropic_and_rotated_tensor_states_follow_original_potential() {
    let card = tabulated_card(true);
    let total = [0.025, -0.012, 0.018, 0.004, -0.006, 0.003];
    for isotropic in [true, false] {
        let law = prepared(&card, isotropic, [280.0, 360.0]);
        let floor = law.minimum_heat_capacity_j_k(0.7).unwrap();
        for t in [287.0, 304.0, 323.0, 347.0] {
            let at = law.at_temperature(t).unwrap();
            let state = at.state(total, 0.7).unwrap();
            let (psi, stress) = expected(t, total, isotropic, 0.7);
            close(state.free_energy_j, psi, 2e-13);
            for (a, b) in state.stress_mandel_pa.into_iter().zip(stress) {
                close(a, b, 3e-12);
            }
            let h = 0.01;
            let samples =
                [-2.0, -1.0, 0.0, 1.0, 2.0].map(|x| expected(t + x * h, total, isotropic, 0.7).0);
            let entropy =
                -(samples[0] - 8.0 * samples[1] + 8.0 * samples[3] - samples[4]) / (12.0 * h);
            let capacity = -t
                * (-samples[0] + 16.0 * samples[1] - 30.0 * samples[2] + 16.0 * samples[3]
                    - samples[4])
                / (12.0 * h * h);
            close(state.entropy_j_k, entropy, 3e-10);
            close(state.heat_capacity_j_k, capacity, 3e-6);
            assert!(state.heat_capacity_j_k >= floor);
            assert!(at.temperature_is_smooth());
            at.temperature_derivatives().unwrap();
            let free = at.free_strain_mandel();
            close(at.state(free, 0.7).unwrap().free_energy_j, 0.0, 1e-25);
        }
    }
}

#[test]
fn cte_knots_keep_continuous_energy_and_explicit_one_sided_capacity() {
    let card = tabulated_card(true);
    let total = [0.025, -0.012, 0.018, 0.004, -0.006, 0.003];
    for isotropic in [true, false] {
        let law = prepared(&card, isotropic, [280.0, 360.0]);
        for t in [310.0, 330.0] {
            let at = law.at_temperature(t).unwrap();
            assert!(!at.temperature_is_smooth());
            assert!(at.temperature_derivatives().is_err());
            let state = at.state(total, 0.7).unwrap();
            let h = 1e-4;
            let left = law
                .at_temperature(t - h)
                .unwrap()
                .state(total, 0.7)
                .unwrap();
            let right = law
                .at_temperature(t + h)
                .unwrap()
                .state(total, 0.7)
                .unwrap();
            close(
                state.entropy_j_k,
                0.5 * (left.entropy_j_k + right.entropy_j_k),
                1e-7,
            );
            close(
                state.internal_energy_j,
                0.5 * (left.internal_energy_j + right.internal_energy_j),
                4e-5,
            );
            close(
                state.heat_capacity_j_k,
                t * (right.entropy_j_k - state.entropy_j_k) / h,
                2e-6,
            );
            assert!(
                at.adiabatic_acoustics(
                    super::super::super::acoustic::AdiabaticThermoelasticInput {
                        reference_density_kg_m3: 2.0,
                        caloric_heat_capacity_j_m3_k: 1000.0,
                        total_strain_mandel: total,
                        mechanical_material_identity: ContentHash([1; 32]),
                        caloric_state_identity: ContentHash([2; 32]),
                    }
                )
                .is_err()
            );
            law.prepare_exchange(t, t + 1.0)
                .unwrap()
                .evaluate(total, total, 0.7)
                .unwrap();
        }
        let bounded = prepared(&card, isotropic, [280.0, 330.0]);
        let at = bounded.at_temperature(330.0).unwrap();
        let state = at.state(total, 0.7).unwrap();
        let left = bounded
            .at_temperature(330.0 - 1e-4)
            .unwrap()
            .state(total, 0.7)
            .unwrap();
        close(
            state.heat_capacity_j_k,
            330.0 * (state.entropy_j_k - left.entropy_j_k) / 1e-4,
            2e-6,
        );
    }
}

#[test]
fn finite_exchange_integrates_every_source_segment_and_reverses() {
    let card = tabulated_card(true);
    let e0 = [0.025, -0.012, 0.018, 0.004, -0.006, 0.003];
    let e1 = [0.021, -0.009, 0.014, 0.006, -0.003, 0.002];
    for isotropic in [true, false] {
        let law = prepared(&card, isotropic, [280.0, 360.0]);
        let prepared = law.prepare_exchange(287.0, 347.0).unwrap();
        let exchange = prepared.evaluate(e0, e1, 0.7).unwrap();
        let reverse = law.exchange(347.0, e1, 287.0, e0, 0.7).unwrap();
        close(
            exchange.free_energy_change_j,
            exchange.after.free_energy_j - exchange.before.free_energy_j,
            2e-12,
        );
        close(
            exchange.entropy_change_j_k,
            exchange.after.entropy_j_k - exchange.before.entropy_j_k,
            2e-14,
        );
        close(
            exchange.internal_energy_change_j,
            exchange.after.internal_energy_j - exchange.before.internal_energy_j,
            2e-12,
        );
        let s0 = expected(287.0, e0, isotropic, 0.7).1;
        let s1 = expected(347.0, e1, isotropic, 0.7).1;
        let work = 0.35
            * (0..6)
                .map(|i| (s0[i] + s1[i]) * (e1[i] - e0[i]))
                .sum::<f64>();
        close(exchange.mechanical_work_j, work, 2e-13);
        close(
            exchange.reversible_heat_j,
            work - (exchange.after.internal_energy_j - exchange.before.internal_energy_j),
            2e-12,
        );
        for (a, b) in [
            exchange.free_energy_change_j,
            exchange.entropy_change_j_k,
            exchange.internal_energy_change_j,
            exchange.reversible_heat_j,
        ]
        .into_iter()
        .zip([
            reverse.free_energy_change_j,
            reverse.entropy_change_j_k,
            reverse.internal_energy_change_j,
            reverse.reversible_heat_j,
        ]) {
            close(a, -b, 3e-12);
        }
    }
}

#[test]
fn one_ulp_temperature_increments_retain_prestress_heat_on_both_knot_sides() {
    let card = tabulated_card(false);
    let law = prepared(&card, true, [280.0, 360.0]);
    let strain = [0.025, -0.012, 0.018, 0.004, -0.006, 0.003];
    for (t0, t1) in [
        (319.0_f64, 319.0_f64.next_up()),
        (310.0, 310.0_f64.next_up()),
        (310.0, 310.0_f64.next_down()),
    ] {
        let exchange = law.exchange(t0, strain, t1, strain, 0.7).unwrap();
        let side = if t1 < t0 { t1 } else { t0 };
        let capacity = law
            .at_temperature(side)
            .unwrap()
            .state(strain, 0.7)
            .unwrap()
            .heat_capacity_j_k;
        let expected = -capacity * (t1 - t0);
        assert!(expected.abs() > 1e-20);
        close(exchange.reversible_heat_j, expected, expected.abs() * 2e-10);
        assert_eq!(exchange.mechanical_work_j, 0.0);
        assert_eq!(
            exchange.retained_reversible_heat_j(),
            exchange.retained_mechanical_work_j() - exchange.retained_internal_energy_change_j()
        );
    }
}
