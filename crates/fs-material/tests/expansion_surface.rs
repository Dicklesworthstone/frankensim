//! Face geometry follows the same source, reference datum and principal axes.
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
            lower: &point(300.0),
            upper: &point(400.0),
        },
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
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
