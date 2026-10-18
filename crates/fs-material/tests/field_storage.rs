//! G0/G1 uniform-field reductions from synthetic data, not measurements.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Temperature};

fn card(key: &str, dims: Dims, value: PropertyValue, validity: ValidityDomain) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic uniform-field storage fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic linear material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn constant(key: &str, dims: Dims, value: f64) -> MaterialCard {
    card(
        key,
        dims,
        PropertyValue::Scalar { value, dims },
        ValidityDomain::unconstrained(),
    )
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1e-300), "{a} != {b}");
}

#[test]
fn geometry_and_material_set_capacitance_and_the_same_field_energy() {
    let source = constant(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        4e-11,
    );
    let capacitor = ResolvedCapacitor::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
        0.001,
    )
    .unwrap();
    near(capacitor.capacitance_f(), 8e-10);
    let voltage = 12.0;
    let q = capacitor.capacitance_f() * voltage;
    let field = voltage / capacitor.gap_m();
    // Independent field integral: integral(epsilon E^2/2 dV).
    near(
        capacitor.stored_energy_j(q).unwrap(),
        0.5 * 4e-11 * field * field * capacitor.area_m2() * capacitor.gap_m(),
    );
    near(
        capacitor.stored_energy_j(-q).unwrap(),
        capacitor.stored_energy_j(q).unwrap(),
    );
    assert_eq!(capacitor.material().card_identity(), source.content_hash());
    let larger =
        ResolvedCapacitor::try_from_material(capacitor.material().clone(), 0.04, 0.001).unwrap();
    near(larger.capacitance_f(), 2.0 * capacitor.capacitance_f());
    assert_eq!(larger.material(), capacitor.material());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ResolvedCapacitor::try_from_material(capacitor.material().clone(), bad, 1.0).is_err()
        );
        assert!(
            ResolvedCapacitor::try_from_material(capacitor.material().clone(), 1.0, bad).is_err()
        );
        if !bad.is_finite() {
            assert!(capacitor.stored_energy_j(bad).is_err());
        }
    }
    assert_eq!(capacitor.stored_energy_j(0.0).unwrap(), 0.0);
}

#[test]
fn magnetic_material_and_turns_set_flux_linkage_energy() {
    let source = constant(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        2e-4,
    );
    let inductor = ResolvedInductor::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.001,
        0.1,
        50,
    )
    .unwrap();
    near(inductor.inductance_h(), 0.005);
    let current = 3.0;
    let linkage = inductor.inductance_h() * current;
    let h = f64::from(inductor.turns()) * current / inductor.path_length_m();
    near(
        inductor.stored_energy_j(linkage).unwrap(),
        0.5 * 2e-4 * h * h * inductor.area_m2() * inductor.path_length_m(),
    );
    let doubled =
        ResolvedInductor::try_from_material(inductor.material().clone(), 0.001, 0.1, 100).unwrap();
    near(doubled.inductance_h(), 4.0 * inductor.inductance_h());
    assert_eq!(doubled.material(), inductor.material());
    assert!(
        ResolvedInductor::try_from_material(inductor.material().clone(), 0.001, 0.1, 0).is_err()
    );
}

#[test]
fn temperature_curve_refresh_keeps_source_and_refuses_extrapolation() {
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let source = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 4e-11), (400.0, 8e-11)],
            dims,
        },
        ValidityDomain::unconstrained().with("T", 300.0, 400.0),
    );
    let resolve = |t| {
        ResolvedCapacitor::resolve(
            &source,
            &QueryPoint::new().with("T", t).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.02,
            0.001,
        )
    };
    let cold = resolve(300.0).unwrap();
    let warm = resolve(350.0).unwrap();
    near(warm.capacitance_f() / cold.capacitance_f(), 1.5);
    assert_eq!(
        cold.material().properties()[0].answer().receipt.selected,
        warm.material().properties()[0].answer().receipt.selected
    );
    assert_ne!(cold.material(), warm.material());
    assert!(resolve(401.0).is_err());
}

#[test]
fn incompatible_quantities_dispersive_curves_and_field_support_refuse() {
    let relative = constant(ELECTRIC_PERMITTIVITY_PROPERTY, Dims([0; 6]), 4.0);
    assert!(
        ResolvedCapacitor::resolve(
            &relative,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0
        )
        .is_err()
    );
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let dispersive = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Curve {
            abscissa: "omega".into(),
            abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
            knots: vec![(1.0, 4e-11), (100.0, 8e-11)],
            dims,
        },
        ValidityDomain::unconstrained().with("omega", 1.0, 100.0),
    );
    assert!(
        ResolvedCapacitor::resolve(
            &dispersive,
            &QueryPoint::new().with("omega", 10.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0
        )
        .is_err()
    );
    let restricted = card(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        PropertyValue::Scalar {
            value: 2e-4,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained().with("B", 0.0, 1.0),
    );
    assert!(
        ResolvedInductor::resolve(
            &restricted,
            &QueryPoint::new().with("B", 0.5).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
            1
        )
        .is_err()
    );
    let frequency_limited = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Scalar { value: 4e-11, dims },
        ValidityDomain::unconstrained().with("frequency", 100.0, 1000.0),
    );
    assert!(
        ResolvedCapacitor::resolve(
            &frequency_limited,
            &QueryPoint::new().with("frequency", 500.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0
        )
        .is_err()
    );
}
