//! Analytic synthetic expansion paths; these are not measured material data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("p", pressure)
        .unwrap()
}
fn claim(name: &str, scale: f64, upper: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300., scale), (340., 3. * scale), (upper, 2. * scale)],
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with("T", 300., upper)
            .with("p", 1e5, 2e5),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic directional expansion; not measurements".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(count: usize, upper: f64, dims: Dims) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, scale) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .zip([1e-5, 2e-5, -1e-5])
        .take(count)
    {
        claims
            .insert_claim(claim(name, scale, upper, dims))
            .unwrap();
    }
    // An isotropic claim must never complete a missing directional law.
    claims
        .insert_claim(claim(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            1e-5,
            400.,
            INVERSE_TEMPERATURE_DIMS,
        ))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic solid".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn path(card: &MaterialCard, from: f64, to: f64) -> IntegratedOrthotropicThermalExpansion {
    integrate_orthotropic_thermal_expansion(
        card,
        &point(from, 1e5),
        &point(to, 1e5),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-12 * a.abs().max(b.abs()).max(1.),
        "{a} != {b}"
    );
}

#[test]
fn exact_piecewise_integral_preserves_signed_axes_and_isotropic_limit() {
    let card = card(3, 400., INVERSE_TEMPERATURE_DIMS);
    let result = path(&card, 300., 400.);
    // Independent trapezoids: 40*(1+3)/2 + 60*(3+2)/2 = 230.
    for (actual, scale) in result
        .free_principal_strains()
        .iter()
        .zip([1e-5, 2e-5, -1e-5])
    {
        near(*actual, 230. * scale);
        assert!((*actual - 200. * scale).abs() > 1e-5);
    }
    let isotropic = integrate_isotropic_thermal_expansion(
        &card,
        &point(300., 1e5),
        &point(400., 1e5),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        result.free_principal_strains()[0].to_bits(),
        isotropic.free_linear_strain().to_bits()
    );
    for (name, selected) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .zip(result.selected_claims())
    {
        assert_eq!(
            *selected,
            result
                .reference()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
        assert_eq!(
            *selected,
            result
                .current()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    assert_eq!(result, path(&card, 300., 400.));
}

#[test]
fn changing_reference_is_additive_and_cooling_reverses_the_path() {
    let card = card(3, 400., INVERSE_TEMPERATURE_DIMS);
    let all = path(&card, 300., 400.);
    let early = path(&card, 300., 325.);
    let late = path(&card, 325., 400.);
    let cooling = path(&card, 400., 300.);
    for i in 0..3 {
        near(
            all.free_principal_strains()[i],
            early.free_principal_strains()[i] + late.free_principal_strains()[i],
        );
        near(
            cooling.free_principal_strains()[i],
            -all.free_principal_strains()[i],
        );
    }
    assert_eq!(path(&card, 350., 350.).free_principal_strains(), &[0.; 3]);
    assert_ne!(all.identity(), cooling.identity());
    assert_ne!(all.identity(), late.identity());
}

#[test]
fn missing_wrong_units_unsupported_endpoint_and_changed_context_refuse() {
    for card in [
        card(2, 400., INVERSE_TEMPERATURE_DIMS),
        card(3, 400., Dims::NONE),
        card(3, 350., INVERSE_TEMPERATURE_DIMS),
    ] {
        assert!(
            integrate_orthotropic_thermal_expansion(
                &card,
                &point(300., 1e5),
                &point(400., 1e5),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let card = card(3, 400., INVERSE_TEMPERATURE_DIMS);
    // Both points independently fit every source domain; the context change
    // still cannot be interpreted as a one-coordinate temperature integral.
    assert!(
        integrate_orthotropic_thermal_expansion(
            &card,
            &point(300., 1e5),
            &point(400., 2e5),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    assert!(
        integrate_orthotropic_thermal_expansion(
            &card,
            &point(300., 1e5),
            &point(401., 1e5),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn competing_claims_require_a_single_explicit_source_for_each_axis() {
    let source = card(3, 400., INVERSE_TEMPERATURE_DIMS);
    let mut claims = source.claims().clone();
    let name = ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[0];
    let alternate = claims
        .insert_claim(claim(name, 7e-5, 400., INVERSE_TEMPERATURE_DIMS))
        .unwrap();
    let ambiguous = MaterialCard::assemble(source.id().clone(), claims, vec![]).unwrap();
    assert!(
        integrate_orthotropic_thermal_expansion(
            &ambiguous,
            &point(300., 1e5),
            &point(400., 1e5),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let original = path(&source, 300., 400.);
    let pinned = ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .zip(original.selected_claims().iter().copied())
        .map(|(n, id)| (n.to_owned(), if n == name { alternate } else { id }))
        .collect();
    let selected = integrate_orthotropic_thermal_expansion(
        &ambiguous,
        &point(300., 1e5),
        &point(400., 1e5),
        MaterialPropertySelection::PinnedByProperty(pinned),
    )
    .unwrap();
    near(selected.free_principal_strains()[0], 230. * 7e-5);
    near(
        selected.free_principal_strains()[1],
        original.free_principal_strains()[1],
    );
    assert_eq!(selected.selected_claims()[0], alternate);
}

fn surface_card(gap: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (axis, name) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .enumerate()
    {
        let (samples, triangles) = if gap && axis == 1 {
            (
                vec![
                    ([300., 1e5], 1e-6),
                    ([350., 1e5], 2e-6),
                    ([300., 2e5], 2e-6),
                    ([350., 2e5], 3e-6),
                    ([450., 1e5], 3e-6),
                    ([500., 1e5], 4e-6),
                    ([450., 2e5], 5e-6),
                    ([500., 2e5], 6e-6),
                ],
                vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
            )
        } else {
            (
                vec![
                    ([300., 1e5], 1e-6),
                    ([500., 1e5], 3e-6),
                    ([300., 2e5], 2e-6),
                    ([500., 2e5], 6e-6),
                ],
                vec![[0, 1, 2], [1, 3, 2]],
            )
        };
        let mut source = claim(name, 1e-6, 500., INVERSE_TEMPERATURE_DIMS);
        source.value = PropertyValue::Surface2D {
            axes: ["T".into(), "p".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples,
            triangles,
            dims: INVERSE_TEMPERATURE_DIMS,
        };
        source.validity = ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.,
                500.,
            )
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 1e5, 2e5);
        claims.insert_claim(source).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic surface expansion".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn surface_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 1.5e5)
        .unwrap()
}

#[test]
fn each_axis_uses_the_source_surface_integral_and_a_gap_refuses_all_axes() {
    let complete = surface_card(false);
    let result = integrate_orthotropic_thermal_expansion(
        &complete,
        &surface_point(300.),
        &surface_point(500.),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    // At this fixed pressure the source edge is crossed at 400 K. The
    // endpoint and edge values are 1.5, 2.5 and 4.5 microstrain/K.
    for strain in result.free_principal_strains() {
        near(*strain, 550e-6);
    }
    let gap = surface_card(true);
    for t in [300., 500.] {
        integrate_orthotropic_thermal_expansion(
            &gap,
            &surface_point(t),
            &surface_point(t),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    }
    for (from, to) in [(300., 500.), (500., 300.)] {
        assert!(matches!(
            integrate_orthotropic_thermal_expansion(
                &gap,
                &surface_point(from),
                &surface_point(to),
                MaterialPropertySelection::SingleClaimOnly
            ),
            Err(MaterialStatePointError::Support { .. })
        ));
    }
}
