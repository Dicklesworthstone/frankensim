use super::*;
use crate::state_point::{
    MaterialPropertySelection, ScalarPropertyRequirement, YIELD_STRESS_PROPERTY,
    resolve_isotropic_solid_state_point, resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Density, Dims, Pressure};

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn card(name: &str, e: [(f64, f64); 3], nu: [f64; 2], yield_stress: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, knots) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, e.to_vec()),
        (
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            vec![(200.0, nu[0]), (600.0, nu[1])],
        ),
        (
            YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            vec![(200.0, yield_stress[0]), (600.0, yield_stress[1])],
        ),
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            vec![(200.0, 1.0), (600.0, 1.0)],
        ),
    ] {
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
                    point(300.0).axis_quantities()["T"],
                    200.0,
                    600.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored thermal contact reference".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn affine(name: &str, e: f64, slope: f64, nu: f64) -> MaterialCard {
    card(
        name,
        [200.0, 400.0, 600.0].map(|t| (t, slope.mul_add(t - 300.0, e))),
        [nu; 2],
        [1e12; 2],
    )
}
fn state(card: &MaterialCard, t: f64) -> IsotropicSolidStatePoint {
    resolve_isotropic_solid_state_point(card, &point(t), MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}
fn pair(cards: &[MaterialCard; 2]) -> ThermalHertzPair {
    ThermalHertzPair::new(
        [&cards[0], &cards[1]],
        [&state(&cards[0], 300.0), &state(&cards[1], 320.0)],
        [[200.0, 600.0]; 2],
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn original_compliance_and_all_thermal_derivatives_match_independent_differences() {
    let cards = [affine("a", 2e9, -2e6, 0.25), affine("b", 3e9, -1e6, 0.3)];
    let model = pair(&cards);
    let t = [360.0, 420.0];
    let at = model.at_temperatures(t).unwrap();
    let independent = |t: [f64; 2]| {
        ((1.0 - 0.25 * 0.25) / (2e9 - 2e6 * (t[0] - 300.0))
            + (1.0 - 0.3 * 0.3) / (3e9 - 1e6 * (t[1] - 300.0)))
            .recip()
    };
    near(at.reduced_modulus.value_pa, independent(t), 3e-7);
    for i in 0..2 {
        let h = 0.05;
        let mut plus = t;
        plus[i] += h;
        let mut minus = t;
        minus[i] -= h;
        near(
            at.gradient_pa_k[i],
            (independent(plus) - independent(minus)) / (2.0 * h),
            0.002,
        );
        near(
            at.hessian_pa_k2[i][i],
            (independent(plus) - 2.0 * independent(t) + independent(minus)) / (h * h),
            0.001,
        );
        assert_eq!(
            at.surfaces[i].resolved().query_point(),
            state(&cards[i], t[i]).resolved().query_point()
        );
    }
    let h = 0.1;
    let mixed = (independent([t[0] + h, t[1] + h])
        - independent([t[0] + h, t[1] - h])
        - independent([t[0] - h, t[1] + h])
        + independent([t[0] - h, t[1] - h]))
        / (4.0 * h * h);
    near(at.hessian_pa_k2[0][1], mixed, 0.001);
    assert_eq!(at.hessian_pa_k2[0][1], at.hessian_pa_k2[1][0]);
    assert!(at.hessian_pa_k2[0][0] < 0.0 && at.hessian_pa_k2[1][1] < 0.0);
    // A hypothetical native Hertz geometry factor G is fixed in this check.
    // U=psi+sum(Ts), not psi, and dU/dTi=-G sum_j Tj*d2E/dTj dTi.
    let g = 4e-10;
    let energy = |t: [f64; 2]| {
        let s = model.at_temperatures(t).unwrap();
        g * (s.reduced_modulus.value_pa - t[0] * s.gradient_pa_k[0] - t[1] * s.gradient_pa_k[1])
    };
    for i in 0..2 {
        let mut plus = t;
        plus[i] += h;
        let mut minus = t;
        minus[i] -= h;
        near(
            (energy(plus) - energy(minus)) / (2.0 * h),
            -g * (t[0] * at.hessian_pa_k2[0][i] + t[1] * at.hessian_pa_k2[1][i]),
            1e-11,
        );
    }
}

#[test]
fn surface_exchange_and_constant_limit_preserve_original_values() {
    let cards = [affine("a", 2e9, -2e6, 0.25), affine("b", 3e9, 1e6, 0.3)];
    let model = pair(&cards);
    let swap = ThermalHertzPair::new(
        [&cards[1], &cards[0]],
        [&state(&cards[1], 320.0), &state(&cards[0], 300.0)],
        [[200.0, 600.0]; 2],
    )
    .unwrap();
    let a = model.at_temperatures([350.0, 410.0]).unwrap();
    let b = swap.at_temperatures([410.0, 350.0]).unwrap();
    near(a.reduced_modulus.value_pa, b.reduced_modulus.value_pa, 0.0);
    for i in 0..2 {
        near(a.gradient_pa_k[i], b.gradient_pa_k[1 - i], 1e-9);
        for j in 0..2 {
            near(a.hessian_pa_k2[i][j], b.hessian_pa_k2[1 - i][1 - j], 1e-11);
        }
    }
    assert_ne!(model.identity(), swap.identity());
    let fixed = [affine("c", 2e9, 0.0, 0.25), affine("d", 3e9, 0.0, 0.3)];
    let response = pair(&fixed).at_temperatures([200.0, 600.0]).unwrap();
    assert_eq!(response.gradient_pa_k, [0.0; 2]);
    assert_eq!(response.hessian_pa_k2, [[0.0; 2]; 2]);
    assert_eq!(
        response.reduced_modulus.value_pa,
        state(&fixed[0], 300.0)
            .reduced_modulus_against(&state(&fixed[1], 300.0))
            .unwrap()
            .value_pa
    );
}

#[test]
fn original_source_domains_claims_and_differentiability_are_enforced() {
    let a = affine("a", 2e9, -2e6, 0.25);
    let b = affine("b", 3e9, -1e6, 0.3);
    let model = pair(&[a.clone(), b.clone()]);
    assert!(model.at_temperatures([199.9, 300.0]).is_err());
    assert!(model.at_temperatures([300.0, f64::NAN]).is_err());
    let refs = [state(&a, 300.0), state(&b, 320.0)];
    assert!(ThermalHertzPair::new([&b, &b], [&refs[0], &refs[1]], [[200.0, 600.0]; 2]).is_err());
    let kink = card(
        "kink",
        [(200.0, 2.2e9), (400.0, 1.8e9), (600.0, 1.5e9)],
        [0.25; 2],
        [1e12; 2],
    );
    let nu = card(
        "nu",
        [(200.0, 2.2e9), (400.0, 1.8e9), (600.0, 1.4e9)],
        [0.25, 0.26],
        [1e12; 2],
    );
    let yield_change = card(
        "yield",
        [(200.0, 2.2e9), (400.0, 1.8e9), (600.0, 1.4e9)],
        [0.25; 2],
        [1e12, 9e11],
    );
    for unsupported in [&kink, &yield_change] {
        assert!(
            ThermalHertzPair::new(
                [unsupported, &b],
                [&state(unsupported, 300.0), &refs[1]],
                [[200.0, 600.0]; 2]
            )
            .is_err()
        );
    }
    assert!(
        ThermalHertzPair::new(
            [&nu, &b],
            [&state(&nu, 300.0), &refs[1]],
            [[200.0, 600.0]; 2]
        )
        .is_ok(),
        "an affine Poisson source has an explicit thermal derivative"
    );
    assert!(
        ThermalHertzPair::new(
            [&a, &b],
            [&refs[0], &refs[1]],
            [[200.0, 601.0], [200.0, 600.0]]
        )
        .is_err()
    );
    // A rejected query cannot mutate source selection or the next result.
    assert_eq!(
        model.at_temperatures([300.0, 320.0]).unwrap(),
        model.at_temperatures([300.0, 320.0]).unwrap()
    );
}

fn selected_curve(
    name: &str,
    dims: Dims,
    value: f64,
    slope: f64,
    support: [f64; 2],
) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: [200.0, 400.0, 600.0]
                .map(|t| (t, slope.mul_add(t - 300.0, value)))
                .to_vec(),
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            point(300.0).axis_quantities()["T"],
            support[0],
            support[1],
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured independent thermal modulus sources".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn selected(
    name: &str,
    bulk_shear: bool,
    second_slope: f64,
    auxiliary_support: [f64; 2],
) -> (MaterialCard, IsotropicSolidStatePoint) {
    let requirement = |key, dims| {
        ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    };
    let density = requirement("rho-source", Density::DIMS);
    let yield_limit = requirement("yield-source", Pressure::DIMS);
    let first_name = if bulk_shear { "K-source" } else { "E-source" };
    let first_value = if bulk_shear { 1.6e9 } else { 2e9 };
    let first_slope = if bulk_shear { -4e5 } else { -1e6 };
    let mut claims = ClaimSet::new();
    for claim in [
        selected_curve("rho-source", Density::DIMS, 1000.0, 0.0, [200.0, 600.0]),
        selected_curve(
            first_name,
            Pressure::DIMS,
            first_value,
            first_slope,
            [200.0, 600.0],
        ),
        selected_curve(
            "G-source",
            Pressure::DIMS,
            8e8,
            second_slope,
            [200.0, 600.0],
        ),
        selected_curve("yield-source", Pressure::DIMS, 1e12, 0.0, [200.0, 600.0]),
        selected_curve("retained-optical", Dims::NONE, 1.5, 1e-4, auxiliary_support),
        // Unselected canonical answers deliberately disagree with this law.
        selected_curve(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            9e9,
            0.0,
            [200.0, 600.0],
        ),
        selected_curve(POISSON_RATIO_PROPERTY, Dims::NONE, 0.1, 0.0, [200.0, 600.0]),
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic sources".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let first = requirement(first_name, Pressure::DIMS);
    let shear = requirement("G-source", Pressure::DIMS);
    let parameters = if bulk_shear {
        IsotropicElasticParameters::BulkShear {
            bulk_modulus: first,
            shear_modulus: shear,
        }
    } else {
        IsotropicElasticParameters::YoungShear {
            young_modulus: first,
            shear_modulus: shear,
        }
    };
    let mut requirements = vec![
        density.clone(),
        yield_limit.clone(),
        requirement("retained-optical", Dims::NONE),
    ];
    requirements.extend(parameters.requirements().into_iter().cloned());
    let resolved = resolve_material_state_point(
        &material,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let state = IsotropicSolidStatePoint::try_from_resolved_with_parameters(
        resolved,
        &density,
        &parameters,
        &yield_limit,
    )
    .unwrap();
    (material, state)
}

fn check_derivatives(model: &ThermalHertzPair, t: [f64; 2], independent: impl Fn([f64; 2]) -> f64) {
    let at = model.at_temperatures(t).unwrap();
    near(at.reduced_modulus.value_pa, independent(t), 2e-6);
    let h = 0.1;
    for i in 0..2 {
        let mut plus = t;
        plus[i] += h;
        let mut minus = t;
        minus[i] -= h;
        near(
            at.gradient_pa_k[i],
            (independent(plus) - independent(minus)) / (2.0 * h),
            0.02,
        );
        near(
            at.hessian_pa_k2[i][i],
            (independent(plus) - 2.0 * independent(t) + independent(minus)) / (h * h),
            0.02,
        );
    }
    let mixed = (independent([t[0] + h, t[1] + h])
        - independent([t[0] + h, t[1] - h])
        - independent([t[0] - h, t[1] + h])
        + independent([t[0] - h, t[1] - h]))
        / (4.0 * h * h);
    near(at.hessian_pa_k2[0][1], mixed, 0.02);
    assert_eq!(at.hessian_pa_k2[0][1], at.hessian_pa_k2[1][0]);
    // Differentiate the actual reciprocal internal energy at fixed indentation.
    let geometry = 4e-10;
    let internal = |temperatures: [f64; 2]| {
        let value = model.at_temperatures(temperatures).unwrap();
        geometry
            * (value.reduced_modulus.value_pa
                - temperatures[0] * value.gradient_pa_k[0]
                - temperatures[1] * value.gradient_pa_k[1])
    };
    for i in 0..2 {
        let mut plus = t;
        plus[i] += h;
        let mut minus = t;
        minus[i] -= h;
        near(
            (internal(plus) - internal(minus)) / (2.0 * h),
            -geometry * (t[0] * at.hessian_pa_k2[0][i] + t[1] * at.hessian_pa_k2[1][i]),
            2e-10,
        );
    }
}

#[test]
fn independent_pairs_drive_exact_two_temperature_contact_without_canonical_reselection() {
    let (a, sa) = selected("young-shear", false, -2e5, [200.0, 600.0]);
    let (b, sb) = selected("bulk-shear", true, 3e5, [200.0, 600.0]);
    let model = ThermalHertzPair::new([&a, &b], [&sa, &sb], [[200.0, 600.0]; 2]).unwrap();
    let independent = |t: [f64; 2]| {
        let young = 2e9 - 1e6 * (t[0] - 300.0);
        let ga = 8e8 - 2e5 * (t[0] - 300.0);
        let bulk = 1.6e9 - 4e5 * (t[1] - 300.0);
        let gb = 8e8 + 3e5 * (t[1] - 300.0);
        let nua = young / (2.0 * ga) - 1.0;
        let eb = 9.0 * bulk * gb / (3.0 * bulk + gb);
        let nub = (3.0 * bulk - 2.0 * gb) / (2.0 * (3.0 * bulk + gb));
        ((1.0 - nua * nua) / young + (1.0 - nub * nub) / eb).recip()
    };
    check_derivatives(&model, [360.0, 420.0], independent);
    let response = model.at_temperatures([360.0, 420.0]).unwrap();
    for (original, refreshed) in [&sa, &sb].into_iter().zip(&response.surfaces) {
        assert_eq!(
            original.resolved().properties().len(),
            refreshed.resolved().properties().len()
        );
        assert!(
            refreshed
                .resolved()
                .property(YOUNG_MODULUS_PROPERTY)
                .is_none()
        );
        assert!(
            refreshed
                .resolved()
                .property(POISSON_RATIO_PROPERTY)
                .is_none()
        );
        assert!(refreshed.resolved().property("retained-optical").is_some());
        assert_eq!(original.yield_requirement(), refreshed.yield_requirement());
        assert_eq!(
            original.elastic_state().derivation().unwrap().parameters(),
            refreshed.elastic_state().derivation().unwrap().parameters()
        );
        for property in original.resolved().properties() {
            assert_eq!(
                property.answer().receipt.selected,
                refreshed
                    .resolved()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    assert_ne!(response.surfaces[0].poisson_ratio(), sa.poisson_ratio());
    let swapped = ThermalHertzPair::new([&b, &a], [&sb, &sa], [[200.0, 600.0]; 2]).unwrap();
    let reverse = swapped.at_temperatures([420.0, 360.0]).unwrap();
    near(
        response.reduced_modulus.value_pa,
        reverse.reduced_modulus.value_pa,
        0.0,
    );
    near(response.gradient_pa_k[0], reverse.gradient_pa_k[1], 1e-9);
    assert_ne!(model.identity(), swapped.identity());
}

#[test]
fn affine_poisson_and_young_sources_keep_their_full_reciprocal_derivatives() {
    let a = card(
        "varying-nu",
        [(200.0, 2.2e9), (400.0, 1.8e9), (600.0, 1.4e9)],
        [0.15, 0.35],
        [1e12; 2],
    );
    let b = affine("fixed-nu", 3e9, -1e6, 0.3);
    let model = pair(&[a, b]);
    check_derivatives(&model, [360.0, 420.0], |t| {
        let ea = 2e9 - 2e6 * (t[0] - 300.0);
        let nua = 0.15 + 0.0005 * (t[0] - 200.0);
        ((1.0 - nua * nua) / ea + (1.0 - 0.3 * 0.3) / (3e9 - 1e6 * (t[1] - 300.0))).recip()
    });
}

#[test]
fn independent_source_intervals_and_auxiliary_support_are_admitted_together() {
    let (b, sb) = selected("bulk-shear", true, 3e5, [200.0, 600.0]);
    let (narrow, narrow_state) = selected("narrow-auxiliary", false, -2e5, [250.0, 450.0]);
    assert!(
        ThermalHertzPair::new([&narrow, &b], [&narrow_state, &sb], [[200.0, 600.0]; 2]).is_err()
    );
    let (unstable, unstable_state) = selected("unstable-young-shear", false, -2e6, [200.0, 600.0]);
    assert!(
        ThermalHertzPair::new([&unstable, &b], [&unstable_state, &sb], [[200.0, 600.0]; 2])
            .is_err()
    );
    let (a, sa) = selected("young-shear", false, -2e5, [200.0, 600.0]);
    assert!(ThermalHertzPair::new([&b, &b], [&sa, &sb], [[200.0, 600.0]; 2]).is_err());
    let model = ThermalHertzPair::new([&a, &b], [&sa, &sb], [[200.0, 600.0]; 2]).unwrap();
    let before = model.at_temperatures([350.0, 450.0]).unwrap();
    assert!(model.at_temperatures([601.0, 450.0]).is_err());
    assert_eq!(model.at_temperatures([350.0, 450.0]).unwrap(), before);
}

#[test]
fn affine_poisson_zero_crossing_preserves_near_reference_queries() {
    let a = card(
        "nu-crossing-zero",
        [(200.0, 2e9), (400.0, 2e9), (600.0, 2e9)],
        [0.1, -0.3],
        [1e12; 2],
    );
    let b = affine("fixed-neighbor", 3e9, 0.0, 0.3);
    let model = pair(&[a, b]);
    for temperature in [300.0 - 1e-10, 300.0, 300.0 + 1e-10] {
        let response = model.at_temperatures([temperature, 320.0]).unwrap();
        let nu = -0.001 * (temperature - 300.0);
        near(response.surfaces[0].poisson_ratio(), nu, 1e-14);
        let expected = ((1.0 - nu * nu) / 2e9 + (1.0 - 0.3 * 0.3) / 3e9).recip();
        near(response.reduced_modulus.value_pa, expected, 2e-6);
        assert!(response.gradient_pa_k.iter().all(|value| value.is_finite()));
        assert!(
            response
                .hessian_pa_k2
                .iter()
                .flatten()
                .all(|value| value.is_finite())
        );
    }
}
