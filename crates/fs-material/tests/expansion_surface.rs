//! Face geometry follows the same source, reference datum and principal axes.
#[path = "expansion_surface/oriented.rs"]
mod oriented;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
    surface::PreparedExpansionSurface,
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(), t).unwrap()
}
fn source(reference: f64, values: &[Vec<(f64, f64)>]) -> PreparedFreeExpansion {
    source_in_span(reference, values, [300.0, 400.0])
}
fn source_in_span(
    reference: f64,
    values: &[Vec<(f64, f64)>],
    span: [f64; 2],
) -> PreparedFreeExpansion {
    let mut claims = ClaimSet::new();
    let mut keys = Vec::new();
    for (axis, knots) in values.iter().enumerate() {
        let key = PropertyKey::new(format!("engineering_length_change_{axis}"), Dims::NONE);
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: knots.clone(),
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    300.0,
                    400.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic principal face fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        keys.push(key);
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let strains = if values.len() == 1 {
        FreeExpansionStrainSources::Isotropic {
            strain: &keys[0],
            selection,
        }
    } else {
        FreeExpansionStrainSources::FixedPrincipal {
            frame: ContentHash([19; 32]),
            strains: [&keys[0], &keys[1], &keys[2]],
            selections: [selection; 3],
        }
    };
    PreparedFreeExpansion::try_new(
        &card,
        FreeExpansionGeometryQuery {
            reference: &point(reference),
            strains,
            temperature_axis: "T",
            lower: &point(span[0]),
            upper: &point(span[1]),
        },
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}

fn relative(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-13 * b.abs(), "{a} != {b}");
}

#[test]
fn g1_native_expansion_tangent_survives_equal_rounded_endpoint_lengths() {
    let lower = 350.0_f64;
    let middle = lower.next_up();
    let upper = middle.next_up();
    let change = 2.0_f64.powi(-20);
    let path = source_in_span(
        lower,
        &[vec![(300.0, 0.125), (400.0, 0.125 + change)]],
        [lower, upper],
    );
    let expected = (change / 100.0) / (1.125 + 0.5 * change);
    // Geometry legitimately rounds equal; its differential must not be inferred
    // as zero from that rounding. This also feeds the original face-area owner.
    assert_eq!(
        path.at_temperature(upper).unwrap().principal_stretches(),
        [1.0; 3]
    );
    for t in [lower, middle, upper] {
        let at = path.at_temperature_with_tangent(t).unwrap();
        for actual in at.principal_stretch_derivatives_per_k() {
            relative(actual, expected);
        }
        relative(at.relative_volume_derivative_per_k(), 3.0 * expected);
    }
    let smooth = path.at_temperature_with_smooth_tangent(middle).unwrap();
    relative(smooth.relative_volume_derivative_per_k(), 3.0 * expected);
    let surface = PreparedExpansionSurface::try_new(path, [2.0, 3.0, 5.0]).unwrap();
    relative(
        surface
            .at_temperature(middle)
            .unwrap()
            .total_area_derivative_m2_k(),
        20.0 * expected,
    );
}

#[test]
fn g3_inserted_reference_preserves_signed_principal_tangents_and_real_corners() {
    let change = 2.0_f64.powi(-20);
    let reference = 300.0_f64.next_up();
    let path = source(
        reference,
        &[
            vec![
                (300.0, 0.125),
                (350.0, 0.125 + change),
                (400.0, 0.125 + 3.0 * change),
            ],
            vec![(300.0, 0.25), (400.0, 0.25 - change)],
            vec![(300.0, 0.0), (400.0, 0.0)],
        ],
    );
    // The artificial reference knot is not a discontinuity of any source law.
    let at = path.at_temperature_with_smooth_tangent(reference).unwrap();
    let expected = [change / 50.0 / 1.125, -change / 100.0 / 1.25, 0.0];
    for (actual, expected) in at
        .principal_stretch_derivatives_per_k()
        .into_iter()
        .zip(expected)
    {
        relative(actual, expected);
    }
    // The actual source slope changes at 350 K; the ordinary right derivative
    // remains available while the two-sided thermoelastic derivative refuses.
    assert!(path.at_temperature_with_smooth_tangent(350.0).is_err());
    relative(
        path.at_temperature_with_tangent(350.0)
            .unwrap()
            .principal_stretch_derivatives_per_k()[0],
        2.0 * change / 50.0 / 1.125,
    );
}

#[test]
fn g1_native_surface_expansion_tangent_needs_no_representable_interior_sample() {
    let change = 2.0_f64.powi(-20);
    let key = PropertyKey::new("engineering_length_change", Dims::NONE);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "composition".into()],
                axis_dims: [Temperature::DIMS, Dims::NONE],
                samples: vec![
                    ([300.0, 0.0], 0.125),
                    ([400.0, 0.0], 0.125 + change),
                    ([400.0, 1.0], 0.1875 + change),
                    ([300.0, 1.0], 0.1875),
                ],
                triangles: vec![[0, 1, 2], [0, 2, 3]],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", quantity(), 300.0, 400.0)
                .with_quantity(
                    "composition",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                    1.0,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured affine thermal strain at fixed composition".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = |t| {
        point(t)
            .with_quantity("composition", QuantitySpec::dimensional(Dims::NONE), 0.25)
            .unwrap()
    };
    let lower = 350.0_f64;
    let upper = lower.next_up();
    let path = PreparedFreeExpansion::try_new(
        &card,
        FreeExpansionGeometryQuery {
            reference: &query(lower),
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &key,
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            temperature_axis: "T",
            lower: &query(lower),
            upper: &query(upper),
        },
    )
    .unwrap();
    let expected = (change / 100.0) / (1.140625 + 0.5 * change);
    for t in [lower, upper] {
        let at = path.at_temperature_with_tangent(t).unwrap();
        for actual in at.principal_stretch_derivatives_per_k() {
            relative(actual, expected);
        }
    }
}

#[test]
fn g0_unrepresentable_native_tangent_does_not_remove_valid_geometry() {
    let lower = 300.0_f64;
    let upper = lower.next_up();
    let path = source_in_span(
        lower,
        &[
            vec![(lower, 0.0), (upper, 1.0e307)],
            vec![(300.0, 0.0), (400.0, 0.0)],
            vec![(300.0, 0.0), (400.0, 0.0)],
        ],
        [lower, upper],
    );
    assert_eq!(
        path.at_temperature(upper).unwrap().relative_volume(),
        1.0e307
    );
    assert!(path.at_temperature_with_tangent(lower).is_err());
    assert!(path.at_temperature_with_tangent(upper).is_err());
}

#[test]
fn g1_source_normalization_anisotropic_face_areas_and_tangents_follow_the_same_axes() {
    let path = source(
        300.0,
        &[
            vec![(300.0, 0.2), (350.0, 0.5), (400.0, 1.1)],
            vec![(300.0, 0.1), (350.0, 0.2), (400.0, 0.3)],
            vec![(300.0, -0.1), (350.0, 0.0), (400.0, 0.3)],
        ],
    );
    let reference = [2.0, 3.0, 5.0];
    let surface = PreparedExpansionSurface::try_new(path.clone(), reference).unwrap();
    assert_eq!(surface.source(), &path);
    assert_eq!(surface.reference_areas_m2(), reference);
    assert_eq!(surface.reference_total_area_m2(), 10.0);
    assert_eq!(
        surface.at_temperature(300.0).unwrap().principal_areas_m2(),
        reference
    );
    for (t, e, de) in [
        (340.0, [0.44, 0.18, -0.02], [0.006, 0.002, 0.002]),
        (350.0, [0.5, 0.2, 0.0], [0.012, 0.002, 0.006]),
        (400.0, [1.1, 0.3, 0.3], [0.012, 0.002, 0.006]),
    ] {
        let s: [f64; 3] = core::array::from_fn(|i| (1.0 + e[i]) / [1.2, 1.1, 0.9][i]);
        let ds: [f64; 3] = core::array::from_fn(|i| de[i] / [1.2, 1.1, 0.9][i]);
        let sample = surface.at_temperature(t).unwrap();
        let mut total = 0.0;
        let mut tangent = 0.0;
        for i in 0..3 {
            let [j, k] = [(i + 1) % 3, (i + 2) % 3];
            let area = reference[i] * s[j] * s[k];
            let derivative = reference[i] * (ds[j] * s[k] + s[j] * ds[k]);
            near(sample.principal_areas_m2()[i], area);
            near(sample.principal_area_derivatives_m2_k()[i], derivative);
            total += area;
            tangent += derivative;
        }
        near(sample.total_area_m2(), total);
        near(sample.total_area_derivative_m2_k(), tangent);
        assert_eq!(sample.geometry().path_identity(), path.identity());
        assert_eq!(
            surface
                .from_geometry(path.at_temperature_with_tangent(t).unwrap())
                .unwrap(),
            sample
        );
    }
    let dt = 0.001;
    let finite_difference = (surface.at_temperature(340.0 + dt).unwrap().total_area_m2()
        - surface.at_temperature(340.0 - dt).unwrap().total_area_m2())
        / (2.0 * dt);
    assert!(
        (finite_difference
            - surface
                .at_temperature(340.0)
                .unwrap()
                .total_area_derivative_m2_k())
        .abs()
            < 1e-9
    );
}

#[test]
fn g3_patch_selection_zero_area_and_source_derivative_cancellation_are_physical() {
    let path = source(
        300.0,
        &[
            vec![(300.0, 0.0), (400.0, 1.0)],
            vec![(300.0, 0.0), (400.0, 0.0)],
            vec![(300.0, 0.0), (400.0, 0.0)],
        ],
    );
    let face = PreparedExpansionSurface::try_new(path.clone(), [2.0, 3.0, 0.0]).unwrap();
    let at = face.at_temperature(400.0).unwrap();
    assert_eq!(at.principal_areas_m2(), [2.0, 6.0, 0.0]);
    near(at.principal_area_derivatives_m2_k()[0], 0.0);
    near(at.principal_area_derivatives_m2_k()[1], 0.03);
    let empty = PreparedExpansionSurface::try_new(path, [0.0; 3]).unwrap();
    assert_eq!(empty.at_temperature(350.0).unwrap().total_area_m2(), 0.0);
    assert!(empty.at_temperature(401.0).is_err());
    let cancelling = source(
        300.0,
        &[
            vec![(300.0, 0.0), (400.0, 0.0)],
            vec![(300.0, -0.5), (400.0, 0.5)],
            vec![(300.0, 0.5), (400.0, -0.5)],
        ],
    );
    let selected = PreparedExpansionSurface::try_new(cancelling, [1.0, 0.0, 0.0]).unwrap();
    near(
        selected
            .at_temperature(350.0)
            .unwrap()
            .total_area_derivative_m2_k(),
        0.0,
    );
}

#[test]
fn g4_foreign_reference_invalid_areas_and_source_domain_never_supply_a_surface() {
    let values = [vec![(300.0, 0.2), (400.0, 0.4)]];
    let path = source(300.0, &values);
    for area in [[-1.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0], [f64::MAX; 3]] {
        assert!(PreparedExpansionSurface::try_new(path.clone(), area).is_err());
    }
    let face = PreparedExpansionSurface::try_new(path.clone(), [1.0, 2.0, 3.0]).unwrap();
    let foreign = source(325.0, &values);
    assert!(
        face.from_geometry(foreign.at_temperature_with_tangent(350.0).unwrap())
            .is_err()
    );
    assert!(face.at_temperature(299.0).is_err());
    assert!(face.at_temperature(f64::NAN).is_err());
    assert_eq!(face.source(), &path);
}
