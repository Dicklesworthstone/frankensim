use super::*;
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
    for unsupported in [&kink, &nu, &yield_change] {
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
