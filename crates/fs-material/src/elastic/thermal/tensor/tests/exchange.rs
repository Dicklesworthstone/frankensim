//! Finite exchange on the same synthetic original tensor/source fixtures.
use super::*;

fn bulk_shear_law(bulk: f64, shear: f64, alpha: f64) -> AffineThermoelastic {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (BULK_MODULUS_PROPERTY, Pressure::DIMS, bulk),
        (SHEAR_MODULUS_PROPERTY, Pressure::DIMS, shear),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, Dims([0, 0, 0, -1, 0, 0]), alpha),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", quantity(), 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite thermal exchange".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "independent bulk and shear".into(),
            phase: "solid".into(),
            process: "finite exchange test".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
    let requirements = [BULK_MODULUS_PROPERTY, SHEAR_MODULUS_PROPERTY].map(|name| {
        ScalarPropertyRequirement::try_new(
            name, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive,
        ).unwrap()
    });
    let expansion = ScalarPropertyRequirement::try_new(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        Dims([0, 0, 0, -1, 0, 0]),
        ScalarAdmissibility::Finite,
    ).unwrap();
    let elastic = resolve_material_state_point(
        &card, &point(300.0), &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let expansion = resolve_material_state_point(
        &card, &point(300.0), &[expansion], MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    AffineThermoelastic::new(
        &card, &elastic, &expansion,
        ThermoelasticSymmetry::IsotropicBulkShear,
        ThermoelasticSymmetry::Isotropic, [280.0, 360.0], 0.05,
    ).unwrap()
}

fn entropy_oracle(t: f64, strain: [f64; 6], volume: f64) -> f64 {
    let h = 0.02;
    let values = [-2.0, -1.0, 1.0, 2.0]
        .map(|offset| oracle(t + h * offset, strain, QS, QA, volume).0);
    -(values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h)
}

#[test]
fn finite_exchange_retains_original_rotated_compliance_potential_and_work() {
    let material = card(false);
    let source = law(&material, QS, QA);
    let temperatures = [289.0, 342.0];
    let strains = [
        [0.013, -0.004, 0.002, 0.007, -0.003, 0.005],
        [0.009, 0.002, -0.003, -0.004, 0.006, 0.003],
    ];
    let volume = 0.7;
    let exchange = source.exchange(
        temperatures[0], strains[0], temperatures[1], strains[1], volume,
    ).unwrap();
    let [(psi0, stress0), (psi1, stress1)] = [0, 1].map(|i| {
        oracle(temperatures[i], strains[i], QS, QA, volume)
    });
    let [s0, s1] = [0, 1].map(|i| entropy_oracle(temperatures[i], strains[i], volume));
    let work = 0.5 * volume * (0..6).map(|i| {
        (stress0[i] + stress1[i]) * (strains[1][i] - strains[0][i])
    }).sum::<f64>();
    let internal = psi1 - psi0 + temperatures[1] * s1 - temperatures[0] * s0;
    close(exchange.mechanical_work_j, work, 2e-12);
    close(exchange.free_energy_change_j, psi1 - psi0, 2e-12);
    close(exchange.entropy_change_j_k, s1 - s0, 4e-10);
    close(exchange.internal_energy_change_j, internal, 3e-7);
    close(exchange.reversible_heat_j, work - internal, 3e-7);
    close(exchange.free_energy_balance_residual_j, 0.0, 2e-12);
    assert_eq!(
        exchange.before,
        source.at_temperature(temperatures[0]).unwrap().state(strains[0], volume).unwrap(),
    );
    assert_eq!(
        exchange.after,
        source.at_temperature(temperatures[1]).unwrap().state(strains[1], volume).unwrap(),
    );

    // Reversing the endpoints reverses every finite transfer. An observer
    // rotation moves both independent source frames and all actual strains.
    let reverse = source.exchange(
        temperatures[1], strains[1], temperatures[0], strains[0], volume,
    ).unwrap();
    let q = transpose(QS);
    let rotated = law(&material, product(q, QS), product(q, QA));
    let transformed = strains.map(|strain| mandel(rotate(tensor(strain), q)));
    let observed = rotated.exchange(
        temperatures[0], transformed[0], temperatures[1], transformed[1], volume,
    ).unwrap();
    for values in [
        [exchange.mechanical_work_j, reverse.mechanical_work_j, observed.mechanical_work_j],
        [exchange.free_energy_change_j, reverse.free_energy_change_j, observed.free_energy_change_j],
        [exchange.entropy_change_j_k, reverse.entropy_change_j_k, observed.entropy_change_j_k],
        [exchange.internal_energy_change_j, reverse.internal_energy_change_j, observed.internal_energy_change_j],
        [exchange.reversible_heat_j, reverse.reversible_heat_j, observed.reversible_heat_j],
    ] {
        close(values[0], -values[1], 2e-12);
        close(values[0], values[2], 2e-11);
    }
}

#[test]
fn finite_exchange_preserves_one_ulp_motion_under_nonzero_thermal_prestress() {
    let (bulk, shear, alpha, temperature, volume) = (1800.0, 600.0, 1e-4, 350.0, 0.7);
    let source = bulk_shear_law(bulk, shear, alpha);
    let before = [0.0001_f64, -0.002, 0.001, 0.003, 0.0, 0.0];
    let mut after = before;
    after[0] = after[0].next_up();
    let increment = after[0] - before[0];
    let exchange = source.exchange(temperature, before, temperature, after, volume).unwrap();
    assert_eq!(exchange.before.elastic_strain_mandel, exchange.after.elastic_strain_mandel);
    assert_eq!(exchange.before.free_energy_j, exchange.after.free_energy_j);
    assert!(exchange.free_energy_change_j != 0.0);

    // Constant K/G imply Delta S = 3 V K alpha Delta tr(epsilon).
    // Poisson's ratio is nonzero (0.35), and free strain is about 0.005.
    let entropy = 3.0 * volume * bulk * alpha * increment;
    let free = alpha * (temperature - 300.0);
    let trace = before[..3].iter().sum::<f64>() - 3.0 * free;
    let sigma_x = bulk * trace + 2.0 * shear * (before[0] - free - trace / 3.0);
    let work = volume * sigma_x * increment;
    close(exchange.mechanical_work_j, work, work.abs() * 2e-12);
    close(exchange.free_energy_change_j, work, work.abs() * 2e-12);
    close(exchange.entropy_change_j_k, entropy, entropy.abs() * 2e-12);
    close(exchange.internal_energy_change_j, work + temperature * entropy,
        (work + temperature * entropy).abs() * 2e-12);
    close(exchange.reversible_heat_j, -temperature * entropy,
        (temperature * entropy).abs() * 2e-12);
}

#[test]
fn finite_exchange_keeps_source_weak_bulk_mode_in_thermal_heat() {
    let (bulk, shear, alpha, volume) = (1e-5, 1e9, 2e-5, 0.7);
    let source = bulk_shear_law(bulk, shear, alpha);
    let strain = [0.002, 0.002, 0.002, 0.0, 0.0, 0.0];
    let (t0, t1) = (310.0, 340.0);
    let exchange = source.exchange(t0, strain, t1, strain, volume).unwrap();
    let e0 = strain[0] - alpha * (t0 - 300.0);
    let e1 = strain[0] - alpha * (t1 - 300.0);
    let free_change = 4.5 * volume * bulk * (e1 - e0) * (e1 + e0);
    let entropy_change = 9.0 * volume * bulk * alpha * (e1 - e0);
    let internal_change = -4.5 * volume * bulk * alpha * alpha * (t1 - t0) * (t1 + t0);
    assert_eq!(exchange.mechanical_work_j, 0.0);
    close(exchange.free_energy_change_j, free_change, free_change.abs() * 2e-12);
    close(exchange.entropy_change_j_k, entropy_change, entropy_change.abs() * 2e-12);
    close(exchange.internal_energy_change_j, internal_change, internal_change.abs() * 2e-12);
    close(exchange.reversible_heat_j, -internal_change, internal_change.abs() * 2e-12);
    assert!(exchange.reversible_heat_j > 0.0);
}

#[test]
fn finite_exchange_refines_to_original_mixed_thermal_derivative() {
    let material = card(false);
    let source = law(&material, QS, QA);
    let temperature = 320.0;
    let strain = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    let rate = [0.0003, -0.0002, 0.0001, -0.0004, 0.0002, 0.0005];
    let temperature_rate = 4.0;
    let volume = 0.7;
    let state = source.at_temperature(temperature).unwrap().state(strain, volume).unwrap();
    let expected = temperature * volume * (0..6).map(|i| {
        state.stress_temperature_derivative_mandel_pa_k[i] * rate[i]
    }).sum::<f64>() - state.heat_capacity_j_k * temperature_rate;
    let mut errors = Vec::new();
    for duration in [0.02, 0.01, 0.005] {
        let after = core::array::from_fn(|i| strain[i] + duration * rate[i]);
        let exchange = source.exchange(
            temperature, strain, temperature + duration * temperature_rate, after, volume,
        ).unwrap();
        errors.push((exchange.reversible_heat_j / duration - expected).abs());
    }
    assert!(errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1]);
}

#[test]
fn finite_exchange_admits_both_original_source_and_strain_endpoints() {
    let material = card(false);
    let source = law(&material, QS, QA);
    let strain = [0.0; 6];
    for (t0, e0, t1, e1, volume) in [
        (279.0, strain, 320.0, strain, 1.0),
        (320.0, strain, 361.0, strain, 1.0),
        (320.0, [0.2; 6], 330.0, strain, 1.0),
        (320.0, strain, 330.0, [0.2; 6], 1.0),
        (320.0, strain, 330.0, strain, 0.0),
        (320.0, strain, f64::NAN, strain, 1.0),
    ] {
        assert!(source.exchange(t0, e0, t1, e1, volume).is_err());
    }
    let stationary = source.exchange(320.0, strain, 320.0, strain, 1.0).unwrap();
    assert_eq!(stationary.mechanical_work_j, 0.0);
    assert_eq!(stationary.free_energy_change_j, 0.0);
    assert_eq!(stationary.entropy_change_j_k, 0.0);
    assert_eq!(stationary.internal_energy_change_j, 0.0);
    assert_eq!(stationary.reversible_heat_j, 0.0);
}

#[test]
fn finite_exchange_retains_one_ulp_temperature_change_under_thermal_prestress() {
    let (bulk, shear, alpha, volume) = (1800.0, 600.0, 1e-4, 0.7);
    let source = bulk_shear_law(bulk, shear, alpha);
    let t0 = 350.0_f64;
    let t1 = t0.next_up();
    let dt = t1 - t0;
    let strain = [0.0001, -0.002, 0.001, 0.003, 0.0, 0.0];
    let free0 = source.at_temperature(t0).unwrap().free_strain_mandel()[0];
    let free1 = source.at_temperature(t1).unwrap().free_strain_mandel()[0];
    assert_ne!(free1 - free0, alpha * dt);
    let exchange = source.exchange(t0, strain, t1, strain, volume).unwrap();
    let entropy = -9.0 * volume * bulk * alpha * alpha * dt;
    let internal = -4.5 * volume * bulk * alpha * alpha * (t1 + t0) * dt;
    assert_eq!(exchange.mechanical_work_j, 0.0);
    assert_eq!(exchange.stiffness_work_j, 0.0);
    close(exchange.entropy_change_j_k, entropy, entropy.abs() * 2e-12);
    close(exchange.internal_energy_change_j, internal, internal.abs() * 2e-12);
    close(exchange.reversible_heat_j, -internal, internal.abs() * 2e-12);
    close(exchange.expansion_work_j, exchange.free_energy_change_j,
        exchange.free_energy_change_j.abs() * 2e-12);
}

#[test]
fn affine_young_modulus_has_no_spurious_internal_heat_when_endpoint_stiffness_rounds_equal() {
    let original = card(false);
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        let replacement = match claim.key.name() {
            POISSON_RATIO_PROPERTY => Some((NU[0], Dims::NONE)),
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY => {
                Some((0.0, Dims([0, 0, 0, -1, 0, 0])))
            }
            _ => None,
        };
        if let Some((value, dims)) = replacement {
            claim.value = PropertyValue::Scalar { value, dims };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    let material = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let (elastic, expansion) = sources(&material, 300.0, true, true);
    let source = AffineThermoelastic::new(
        &material, &elastic, &expansion,
        ThermoelasticSymmetry::Isotropic, ThermoelasticSymmetry::Isotropic,
        [280.0, 360.0], 0.1,
    ).unwrap();
    let t0 = 350.0_f64;
    let t1 = t0.next_up();
    let strain = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    let before = source.at_temperature(t0).unwrap();
    let after = source.at_temperature(t1).unwrap();
    assert_eq!(before.stiffness_mandel_pa(), after.stiffness_mandel_pa());
    let exchange = source.exchange(t0, strain, t1, strain, 0.7).unwrap();
    let expected = exchange.before.free_energy_j * ET[0]
        / (E[0] + ET[0] * (t0 - 300.0)) * (t1 - t0);
    assert!(expected != 0.0);
    close(exchange.free_energy_change_j, expected, expected.abs() * 2e-12);
    close(exchange.stiffness_work_j, expected, expected.abs() * 2e-12);
    assert_eq!(exchange.mechanical_work_j, 0.0);
    assert_eq!(exchange.expansion_work_j, 0.0);
    assert_eq!(exchange.entropy_change_j_k, 0.0);
    // Psi is affine in T and S=-Psi_T, so U=Psi+T*S is constant.
    // The allowance covers arithmetic in the native rounded entropy snapshot.
    close(exchange.internal_energy_change_j, 0.0, expected.abs() * 2e-12);
    close(exchange.reversible_heat_j, 0.0, expected.abs() * 2e-12);
}

#[test]
fn tiny_temperature_exchange_retains_original_compliance_derivative_increment() {
    let material = card(false);
    let strain = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    let t0 = 321.0_f64;
    let t1 = t0.next_up();
    for elastic_isotropic in [false, true] {
        for expansion_isotropic in [false, true] {
            let (elastic, expansion) = sources(&material, 300.0,
                elastic_isotropic, expansion_isotropic);
            let source = AffineThermoelastic::new(
                &material, &elastic, &expansion,
                if elastic_isotropic {
                    ThermoelasticSymmetry::Isotropic
                } else {
                    ThermoelasticSymmetry::OrientedOrthotropic { principal_to_target: QS }
                },
                if expansion_isotropic {
                    ThermoelasticSymmetry::Isotropic
                } else {
                    ThermoelasticSymmetry::OrientedOrthotropic { principal_to_target: QA }
                },
                [280.0, 360.0], 0.1,
            ).unwrap();
            let exchange = source.exchange(t0, strain, t1, strain, 0.7).unwrap();
            let expected = -exchange.before.heat_capacity_j_k * (t1 - t0);
            assert!(expected.abs() > 1e-30);
            close(exchange.reversible_heat_j, expected, expected.abs() * 5e-10);
            let retained = exchange.retained_mechanical_work_j()
                - exchange.retained_internal_energy_change_j();
            assert_eq!(retained, exchange.retained_reversible_heat_j());
        }
    }
}
