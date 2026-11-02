//! Manufactured nonlinear diffusivity laws; not measured porous-material data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};

const D: Dims = Dims([2, 0, -1, 0, 0, 0]);
const NAME: &str = "moisture_diffusivity";

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(m: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "moisture_content",
            quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
            m,
        )
        .unwrap()
}
fn card(value: PropertyValue, interpolation: InterpolationPolicy) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(NAME, value.dims()),
            value,
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    290.0,
                    310.0,
                )
                .with_quantity(
                    "moisture_content",
                    quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
                    0.0,
                    1.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: vec![],
            provenance: Provenance {
                source: "authored nonlinear diffusion verification law".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured sorbing solid".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(knots: Vec<(f64, f64)>) -> MaterialCard {
    card(
        PropertyValue::Curve {
            abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE,
            knots,
            dims: D,
        },
        InterpolationPolicy::LinearInside,
    )
}
fn binding(c: &MaterialCard, p: &QueryPoint) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        c,
        p,
        &[ScalarPropertyRequirement::try_new(NAME, D, ScalarAdmissibility::NonNegative).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn chart(c: &MaterialCard, bounds: [f64; 2]) -> KirchhoffMoistureDiffusivity {
    KirchhoffMoistureDiffusivity::from_binding(c, &binding(c, &point(bounds[0])), NAME, bounds)
        .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn degenerate_linear_diffusivity_has_nonzero_wetting_flux_and_exact_potential() {
    // D(M)=2e-9 M, K(M)=1e-9 M². The dry endpoint is not a barrier.
    let c = curve(vec![(0.0, 0.0), (1.0, 2e-9)]);
    let k = chart(&c, [0.0, 1.0]);
    assert_eq!(k.diffusivity(0.0).unwrap(), 0.0);
    near(k.mean_diffusivity(0.0, 0.3).unwrap(), 3e-10, 1e-25);
    near(k.integral(0.0, 0.3).unwrap(), 9e-11, 1e-25);
    near(k.integral(0.3, 0.1).unwrap(), -8e-11, 1e-25);
    assert_eq!(k.integral(0.1, 0.1).unwrap(), 0.0);
    assert_eq!(k.mean_diffusivity(0.1, 0.3), k.mean_diffusivity(0.3, 0.1));
    assert_eq!(k, chart(&c, [0.0, 1.0]));
    assert_eq!(k.moisture_bounds(), [0.0, 1.0]);
    assert_eq!(k.knot_count(), 2);
}

#[test]
fn every_crossed_source_segment_contributes_to_the_flux() {
    let c = curve(vec![(0.0, 1e-9), (0.25, 5e-9), (0.5, 1e-9), (1.0, 3e-9)]);
    let k = chart(&c, [0.0, 1.0]);
    // Two triangular ramps followed by a linear rise: 3e-9*.5 + 2e-9*.5.
    near(k.integral(0.0, 1.0).unwrap(), 2.5e-9, 1e-24);
    near(k.mean_diffusivity(0.125, 0.375).unwrap(), 4e-9, 1e-24);
    near(
        k.integral(0.0, 0.375).unwrap() + k.integral(0.375, 1.0).unwrap(),
        k.integral(0.0, 1.0).unwrap(),
        1e-24,
    );
    let narrow = f64::from_bits(0.9f64.to_bits() + 1);
    near(
        k.mean_diffusivity(0.9, narrow).unwrap(),
        k.diffusivity(0.9).unwrap(),
        1e-24,
    );
    assert!(k.integral(0.9, narrow).unwrap() > 0.0);
}

#[test]
fn constant_impermeability_and_zero_width_paths_have_exact_limits() {
    for value in [0.0, 1e-9] {
        let c = card(
            PropertyValue::Scalar { value, dims: D },
            InterpolationPolicy::ConstantWithinValidity,
        );
        let full = chart(&c, [0.0, 1.0]);
        let point = chart(&c, [0.25, 0.25]);
        assert_eq!(full.mean_diffusivity(0.1, 0.9).unwrap(), value);
        assert_eq!(point.mean_diffusivity(0.25, 0.25).unwrap(), value);
        assert_eq!(point.integral(0.25, 0.25).unwrap(), 0.0);
        assert_eq!(point.knot_count(), 1);
        assert!(point.diffusivity(0.2500001).is_err());
    }
}

#[test]
fn negative_intervening_law_wrong_card_and_extrapolation_are_refused() {
    let c = curve(vec![(0.0, 1e-9), (0.5, -1e-9), (1.0, 1e-9)]);
    let bound = binding(&c, &point(0.0));
    assert!(matches!(
        KirchhoffMoistureDiffusivity::from_binding(&c, &bound, NAME, [0.0, 1.0]),
        Err(MaterialStatePointError::OutsideConsumerDomain { .. })
    ));
    let other = curve(vec![(0.0, 1e-9), (1.0, 1e-9)]);
    assert!(KirchhoffMoistureDiffusivity::from_binding(&other, &bound, NAME, [0.0, 1.0]).is_err());
    let k = chart(&other, [0.2, 0.8]);
    for m in [0.19, 0.81, f64::NAN, f64::INFINITY] {
        assert!(k.diffusivity(m).is_err());
        assert!(k.integral(0.3, m).is_err());
    }
    let untyped = QueryPoint::new()
        .with("temperature", 300.0)
        .unwrap()
        .with("moisture_content", 0.2)
        .unwrap();
    // Even when legacy source admission accepts dimensional coordinates, this
    // dry-basis physical owner does not infer their meaning.
    if let Ok(bound) = resolve_material_state_point(
        &other,
        &untyped,
        &[ScalarPropertyRequirement::try_new(NAME, D, ScalarAdmissibility::NonNegative).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    ) {
        assert!(
            KirchhoffMoistureDiffusivity::from_binding(&other, &bound, NAME, [0.2, 0.8]).is_err()
        );
    }
}

fn surface(gap: bool) -> MaterialCard {
    let (samples, triangles) = if gap {
        (
            vec![
                ([0.0, 290.0], 1e-9),
                ([0.25, 290.0], 2e-9),
                ([0.0, 310.0], 2e-9),
                ([0.25, 310.0], 3e-9),
                ([0.75, 290.0], 3e-9),
                ([1.0, 290.0], 4e-9),
                ([0.75, 310.0], 5e-9),
                ([1.0, 310.0], 6e-9),
            ],
            vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
        )
    } else {
        (
            vec![
                ([0.0, 290.0], 1e-9),
                ([1.0, 290.0], 3e-9),
                ([0.0, 310.0], 2e-9),
                ([1.0, 310.0], 6e-9),
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        )
    };
    card(
        PropertyValue::Surface2D {
            axes: ["moisture_content".into(), "temperature".into()],
            axis_dims: [Dims::NONE, Temperature::DIMS],
            samples,
            triangles,
            dims: D,
        },
        InterpolationPolicy::LinearInside,
    )
}

#[test]
fn held_temperature_surface_integrates_crossed_cells_and_cannot_bridge_a_hole() {
    let complete = surface(false);
    let k = chart(&complete, [0.0, 1.0]);
    // At T=300, the actual triangulated slice is (0,1.5),(.5,2.5),(1,4.5) e-9.
    assert_eq!(k.knot_count(), 3);
    near(k.integral(0.0, 1.0).unwrap(), 2.75e-9, 1e-24);
    let gap = surface(true);
    let bound = binding(&gap, &point(0.0));
    binding(&gap, &point(1.0)); // Both endpoints are valid independently.
    assert!(matches!(
        KirchhoffMoistureDiffusivity::from_binding(&gap, &bound, NAME, [0.0, 1.0]),
        Err(MaterialStatePointError::Support { .. })
    ));
}
