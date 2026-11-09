//! Synthetic source-backed regressions for separated bulk/shear thermal modes.
use super::*;
use crate::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ))
}

fn source_law(
    symmetry: ThermoelasticSymmetry,
    coefficients: [[f64; 2]; 2],
    expansion: [[f64; 2]; 3],
    expansion_symmetry: ThermoelasticSymmetry,
) -> AffineThermoelastic {
    let elastic_names = match symmetry {
        ThermoelasticSymmetry::Isotropic => [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (POISSON_RATIO_PROPERTY, Dims::NONE),
        ],
        ThermoelasticSymmetry::IsotropicBulkShear => [
            (BULK_MODULUS_PROPERTY, Pressure::DIMS),
            (SHEAR_MODULUS_PROPERTY, Pressure::DIMS),
        ],
        ThermoelasticSymmetry::IsotropicYoungShear => [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (SHEAR_MODULUS_PROPERTY, Pressure::DIMS),
        ],
        _ => panic!("fixture requires explicitly isotropic stiffness"),
    };
    let expansion_names = if expansion_symmetry == ThermoelasticSymmetry::Isotropic {
        vec![LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY]
    } else {
        ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES.to_vec()
    };
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, coefficient: [f64; 2]| {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, coefficient[0] - 20.0 * coefficient[1]),
                    (360.0, coefficient[0] + 60.0 * coefficient[1]),
                ],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", quantity(), 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic separated thermal modes".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    };
    for (i, (name, dims)) in elastic_names.iter().enumerate() {
        insert(name, *dims, coefficients[i]);
    }
    for (i, name) in expansion_names.iter().enumerate() {
        insert(name, Dims([0, 0, 0, -1, 0, 0]), expansion[i]);
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic isotropic thermal solid".into(),
        phase: "solid".into(),
        process: "weak-mode numerical regression".into(),
        revision: 0,
    }, claims, vec![]).unwrap();
    let elastic_requirements: Vec<_> = elastic_names.into_iter().map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
    }).collect();
    let expansion_requirements: Vec<_> = expansion_names.into_iter().map(|name| {
        ScalarPropertyRequirement::try_new(
            name, Dims([0, 0, 0, -1, 0, 0]), ScalarAdmissibility::Finite,
        ).unwrap()
    }).collect();
    let point = QueryPoint::new().with_quantity("T", quantity(), 300.0).unwrap();
    let elastic = resolve_material_state_point(
        &card, &point, &elastic_requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let thermal = resolve_material_state_point(
        &card, &point, &expansion_requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    AffineThermoelastic::new(
        &card, &elastic, &thermal, symmetry, expansion_symmetry, [280.0, 360.0], 0.1,
    ).unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-30);
    assert!((actual - expected).abs() <= tolerance * scale,
        "{actual:.16e} != {expected:.16e}");
}

// Differentiate the real source-resolved potential at FIXED total strain.
// Elastic capacity alone can be negative; no caloric model is invented here.
fn check_temperature_derivatives(
    law: &AffineThermoelastic, total: [f64; 6], t: f64, volume: f64,
) {
    let h = 1.0 / 128.0;
    let state = law.at_temperature(t).unwrap().state(total, volume).unwrap();
    let low = law.at_temperature(t - h).unwrap().state(total, volume).unwrap();
    let high = law.at_temperature(t + h).unwrap().state(total, volume).unwrap();
    relative(state.entropy_j_k, -(high.free_energy_j - low.free_energy_j) / (2.0 * h), 2e-6);
    relative(state.heat_capacity_j_k, t * (high.entropy_j_k - low.entropy_j_k) / (2.0 * h), 2e-6);
    relative(state.internal_energy_j, state.free_energy_j + t * state.entropy_j_k, 2e-14);
    for i in 0..6 {
        relative(state.stress_temperature_derivative_mandel_pa_k[i],
            (high.stress_mandel_pa[i] - low.stress_mandel_pa[i]) / (2.0 * h), 2e-6);
    }
}

#[test]
fn weak_bulk_source_potential_retains_hydrostatic_work_and_thermal_derivatives() {
    let bulk = 1e-9;
    let bulk_first = 1e-12;
    let alpha = 1.0 / 16384.0;
    let alpha_first = 1.0 / 1048576.0;
    let law = source_law(
        ThermoelasticSymmetry::IsotropicBulkShear,
        [[bulk, bulk_first], [400.0, 0.25]],
        [[alpha, alpha_first]; 3], ThermoelasticSymmetry::Isotropic,
    );
    let total = [0.013, 0.013, 0.013, 0.0, 0.0, 0.0];
    let t = 321.0;
    let volume = 0.7;
    let at = law.at_temperature(t).unwrap();
    let state = at.state(total, volume).unwrap();
    let dt = t - 300.0;
    let free = alpha * dt + 0.5 * alpha_first * dt * dt;
    let trace = 3.0 * (total[0] - free);
    let k = bulk + bulk_first * dt;
    let a = alpha + alpha_first * dt;
    for i in 0..3 {
        relative(state.stress_mandel_pa[i], k * trace, 2e-12);
        relative(state.stress_temperature_derivative_mandel_pa_k[i],
            bulk_first * trace - 3.0 * k * a, 2e-12);
    }
    assert_eq!(&state.stress_mandel_pa[3..], &[0.0; 3]);
    relative(state.free_energy_j, 0.5 * volume * k * trace * trace, 2e-12);
    relative(state.entropy_j_k,
        volume * (3.0 * k * a * trace - 0.5 * bulk_first * trace * trace), 2e-12);
    relative(state.heat_capacity_j_k, -t * volume * (
        -6.0 * bulk_first * a * trace + 9.0 * k * a * a - 3.0 * k * alpha_first * trace
    ), 2e-12);
    check_temperature_derivatives(&law, total, t, volume);
    assert_eq!(at.state(at.free_strain_mandel(), volume).unwrap().free_energy_j, 0.0);
}

#[test]
fn weak_shear_with_independently_rotated_expansion_retains_strength_and_capacity() {
    let alpha = 1.0 / 4096.0;
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let law = source_law(
        ThermoelasticSymmetry::IsotropicBulkShear,
        [[1e20, 0.0], [1.0, 1.0 / 1024.0]],
        [[alpha, 0.0], [-alpha, 0.0], [0.0, 0.0]],
        ThermoelasticSymmetry::OrientedOrthotropic {
            principal_to_target: [[q, -q, 0.0], [q, q, 0.0], [0.0, 0.0, 1.0]],
        },
    );
    let total = [1.0 / 64.0, -1.0 / 64.0, 0.0, 1.0 / 128.0, 0.0, 0.0];
    let volume = 0.7;
    let at = law.at_temperature(300.0).unwrap();
    // The stored dense tensor is deliberately retained unchanged, even when
    // its rounded normal entries cannot encode the original weak shear mode.
    assert_eq!(at.stiffness_mandel_pa()[0][0], at.stiffness_mandel_pa()[0][1]);
    let state = at.state(total, volume).unwrap();
    let alpha_target = at.temperature_derivatives().unwrap().free_strain_mandel_per_k;
    assert!(alpha_target[3].abs() > alpha);
    relative(state.free_energy_j, volume * total.iter().map(|x| x * x).sum::<f64>(), 2e-14);
    for i in 0..6 {
        relative(state.stress_mandel_pa[i], 2.0 * total[i], 2e-14);
        relative(state.stress_temperature_derivative_mandel_pa_k[i],
            2.0 * total[i] / 1024.0 - 2.0 * alpha_target[i], 2e-14);
    }
    let work = total.into_iter().zip(alpha_target).map(|(e, a)| e * a).sum::<f64>();
    let alpha_squared = alpha_target.iter().map(|a| a * a).sum::<f64>();
    relative(state.heat_capacity_j_k,
        -300.0 * volume * (-4.0 * work / 1024.0 + 2.0 * alpha_squared), 2e-14);
    check_temperature_derivatives(&law, total, 300.0, volume);
}

#[test]
fn weak_bulk_young_shear_and_original_young_poisson_share_stable_source_potential() {
    let alpha = 1.0 / 16384.0;
    let alpha_first = 1.0 / 1048576.0;
    let total = [0.013, 0.013, 0.013, 0.0, 0.0, 0.0];
    let t = 321.0;
    let volume = 0.7;
    for (symmetry, coefficients) in [
        (ThermoelasticSymmetry::IsotropicYoungShear, [[1e-9, 1e-12], [400.0, 0.25]]),
        (ThermoelasticSymmetry::Isotropic, [[3.0, 0.01], [-1.0 + 2.0_f64.powi(-40), 0.0]]),
    ] {
        let law = source_law(symmetry, coefficients, [[alpha, alpha_first]; 3],
            ThermoelasticSymmetry::Isotropic);
        let at = law.at_temperature(t).unwrap();
        let state = at.state(total, volume).unwrap();
        let dt = t - 300.0;
        let e = coefficients[0][0] + coefficients[0][1] * dt;
        let second = coefficients[1][0] + coefficients[1][1] * dt;
        let bulk = if symmetry == ThermoelasticSymmetry::IsotropicYoungShear {
            1.0 / (9.0 / e - 3.0 / second)
        } else { e / (3.0 * (1.0 - 2.0 * second)) };
        let free = alpha * dt + 0.5 * alpha_first * dt * dt;
        let trace = 3.0 * (total[0] - free);
        for stress in &state.stress_mandel_pa[..3] {
            relative(*stress, bulk * trace, 2e-12);
        }
        relative(state.free_energy_j, 0.5 * volume * bulk * trace * trace, 2e-12);
        check_temperature_derivatives(&law, total, t, volume);
    }
}
