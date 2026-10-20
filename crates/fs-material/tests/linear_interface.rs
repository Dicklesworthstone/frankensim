//! One source-derived bilateral law shared by physical joint consumers.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ResolvedLinearInterface, ScalarAdmissibility, ScalarPropertyRequirement,
    interface_stiffness_per_area_dims, interface_viscosity_per_area_dims,
    resolve_interface_state_point,
};
use fs_qty::semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0)
            .with("omega", 0.0, 200_000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic scalar interface regression, not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn properties(k: f64, c: f64) -> Vec<PropertyClaim> {
    vec![
        claim("stiffness", interface_stiffness_per_area_dims(), k),
        claim("viscosity", interface_viscosity_per_area_dims(), c),
    ]
}

fn resolve(
    properties: Vec<PropertyClaim>,
    point: &QueryPoint,
    history: &str,
) -> ResolvedInterfaceStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for property in properties {
        requirements.push(
            ScalarPropertyRequirement::try_with_key(&property.key, ScalarAdmissibility::Finite)
                .unwrap(),
        );
        claims.insert_claim(property).unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic surface preparation".into(),
            revision: 1,
        },
        texture_frame: format!("original {name} surface coordinates"),
    };
    let card = InterfaceSystemCard::assemble(
        surface("first"),
        surface("second"),
        SystemContext {
            medium: "synthetic bilateral layer".into(),
            third_body: None,
            environment: "explicit fixed-state fixture".into(),
            history: history.into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_interface_state_point(
        &card,
        point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 100.0)
        .unwrap()
}

fn law(
    source: ResolvedInterfaceStatePoint,
    area: f64,
) -> Result<ResolvedLinearInterface, MaterialStatePointError> {
    ResolvedLinearInterface::new(source, "stiffness", "viscosity", area)
}

#[test]
fn original_ordered_source_and_area_control_one_passive_physical_joint() {
    assert_eq!(
        interface_stiffness_per_area_dims(),
        Dims([-2, 1, -2, 0, 0, 0])
    );
    assert_eq!(
        interface_viscosity_per_area_dims(),
        Dims([-2, 1, -1, 0, 0, 0])
    );
    let source = resolve(properties(1e8, 2e4), &point(300.0), "unloaded source A");
    let result = law(source.clone(), 1e-4).unwrap();
    assert_eq!(result.resolved(), &source);
    assert_eq!(
        (result.stiffness_property(), result.viscosity_property()),
        ("stiffness", "viscosity")
    );
    assert_eq!(result.area_m2(), 1e-4);
    assert_eq!((result.stiffness_n_m(), result.damping_n_s_m()), (1e4, 2.0));
    assert_eq!(result.omega_band_rad_s(), (0.0, 200_000.0));
    let larger = law(source.clone(), 2e-4).unwrap();
    assert_eq!(larger.stiffness_n_m(), 2.0 * result.stiffness_n_m());
    assert_eq!(larger.damping_n_s_m(), 2.0 * result.damping_n_s_m());
    assert_eq!(larger.resolved().identity(), result.resolved().identity());
    let pair = source.surface_materials();
    result
        .require_surface_materials(&pair[0], &pair[1])
        .unwrap();
    assert!(
        result
            .require_surface_materials(&pair[1], &pair[0])
            .is_err()
    );
    let mut revised = pair[0].clone();
    revised.revision += 1;
    assert!(
        result
            .require_surface_materials(&revised, &pair[1])
            .is_err()
    );
    let changed = law(
        resolve(
            properties(1e8, 2e4),
            &point(300.0),
            "different preparation history",
        ),
        1e-4,
    )
    .unwrap();
    assert_eq!(result.stiffness_n_m(), changed.stiffness_n_m());
    assert_ne!(result.resolved().identity(), changed.resolved().identity());
    for (k, c) in [(1e8, 0.0), (0.0, 2e4)] {
        let pure = law(
            resolve(properties(k, c), &point(300.0), "one passive mechanism"),
            1e-4,
        )
        .unwrap();
        assert_eq!(pure.stiffness_n_m(), k * 1e-4);
        assert_eq!(pure.damping_n_s_m(), c * 1e-4);
    }
}

#[test]
fn temperature_requery_changes_coefficients_but_frequency_curves_are_not_frozen() {
    let mut input = properties(1e8, 2e4);
    input[0].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 1e8), (320.0, 3e8)],
        dims: interface_stiffness_per_area_dims(),
    };
    input[0].interpolation = InterpolationPolicy::LinearInside;
    let cool = resolve(input.clone(), &point(280.0), "fixed temperature source");
    let warm = resolve(input.clone(), &point(300.0), "fixed temperature source");
    let a = law(cool.clone(), 1e-4).unwrap();
    let b = law(warm.clone(), 1e-4).unwrap();
    assert_eq!(a.stiffness_n_m(), 1e4);
    assert_eq!(b.stiffness_n_m(), 2e4);
    assert_eq!(a.damping_n_s_m(), b.damping_n_s_m());
    assert_eq!(a.resolved(), &cool);
    assert_eq!(b.resolved(), &warm);
    assert_ne!(a.resolved().identity(), b.resolved().identity());
    input[0].value = PropertyValue::Curve {
        abscissa: "omega".into(),
        abscissa_dims: fs_qty::Frequency::DIMS,
        knots: vec![(0.0, 1e8), (200_000.0, 1e8)],
        dims: interface_stiffness_per_area_dims(),
    };
    let frequency = resolve(input, &point(300.0), "frequency-dependent model");
    assert!(
        law(frequency, 1e-4).is_err(),
        "even equal samples do not change the declared model"
    );
}

#[test]
fn typed_cyclic_frequency_support_is_intersected_and_checked_in_physical_rad_s() {
    let cyclic = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Cyclic),
        ValueForm::Static,
    ));
    let mut input = properties(1e8, 2e4);
    for (p, (lo, hi)) in input.iter_mut().zip([(20.0, 200.0), (30.0, 100.0)]) {
        p.validity = ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0)
            .with_quantity("frequency", cyclic, lo, hi);
    }
    let query = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity("frequency", cyclic, 50.0)
        .unwrap();
    let source = resolve(input, &query, "cyclic-frequency applicability");
    let result = law(source.clone(), 1e-4).unwrap();
    let factor = 2.0 * core::f64::consts::PI;
    let band = (30.0 * factor, 100.0 * factor);
    assert_eq!(result.omega_band_rad_s(), band);
    assert_eq!(result.resolved(), &source);
    for valid in [band.0, 50.0 * factor, band.1] {
        result.require_angular_frequency_rad_s(valid).unwrap();
    }
    for invalid in [
        50.0,
        band.0 - 1.0,
        band.1 + 1.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(result.require_angular_frequency_rad_s(invalid).is_err());
    }
}

#[test]
fn missing_wrong_inactive_unbounded_and_unrepresentable_laws_refuse() {
    let source = resolve(properties(1e8, 2e4), &point(300.0), "admitted source");
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(law(source.clone(), area).is_err());
    }
    let missing = resolve(
        vec![claim("stiffness", interface_stiffness_per_area_dims(), 1e8)],
        &point(300.0),
        "missing viscosity",
    );
    assert!(law(missing, 1e-4).is_err());
    let mut wrong = properties(1e8, 2e4);
    wrong[0] = claim("stiffness", Pressure::DIMS, 1e8);
    assert!(matches!(
        law(
            resolve(
                wrong,
                &point(300.0),
                "bulk modulus is not interface stiffness"
            ),
            1e-4
        ),
        Err(MaterialStatePointError::DimensionMismatch { .. })
    ));
    for (k, c, area) in [
        (0.0, 0.0, 1e-4),
        (-1.0, 2e4, 1e-4),
        (1e8, -1.0, 1e-4),
        (1e300, 1.0, 1e300),
        (1e-300, 1.0, 1e-200),
        (1.0, 1e-300, 1e-200),
    ] {
        assert!(
            law(
                resolve(
                    properties(k, c),
                    &point(300.0),
                    "invalid physical resultant"
                ),
                area
            )
            .is_err()
        );
    }
    let mut unbounded = properties(1e8, 2e4);
    unbounded[1].validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    assert!(
        law(
            resolve(unbounded, &point(300.0), "no frequency support"),
            1e-4
        )
        .is_err()
    );
}
