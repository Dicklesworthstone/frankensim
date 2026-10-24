//! G0/G1/G3 synthetic field thermodynamics; no measured material validation.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
    thermal::{ThermalFieldElement, ThermalFieldStorage},
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

#[test]
fn spatial_magnetic_law_uses_actual_b_volume_and_shared_thermal_derivatives() {
    use fs_material::field_storage::thermal::magnetic::ThermalMagneticMaterial;
    for slope in [-1e-6,0.0,1e-6] {
        let source=curve(MAGNETIC_PERMEABILITY_PROPERTY,MAGNETIC_PERMEABILITY_DIMS,
            vec![(200.0,5e-4-100.0*slope),(500.0,5e-4+200.0*slope)]);
        let original=ResolvedInductor::resolve(&source,&point(300.0),
            MaterialPropertySelection::SingleClaimOnly,0.001,0.1,50).unwrap();
        let law=ThermalMagneticMaterial::new(&source,original.material(),"T",[200.0,500.0]).unwrap();
        let t=350.0; let volume=0.002; let b=0.3;
        let at=law.at_temperature(t).unwrap();
        let actual=at.state(b,volume).unwrap();
        let mu=5e-4+(t-300.0)*slope;
        let psi=volume*b*b/(2.0*mu);
        close(actual.free_energy_j,psi,1e-12);
        if slope!=0.0 {
            close(actual.entropy_j_k,psi*slope/mu,1e-12);
            close(actual.heat_capacity_j_k,-t*volume*b*b*slope*slope/mu.powi(3),1e-12);
            let lo=law.at_temperature(t-0.01).unwrap().state(b,volume).unwrap();
            let hi=law.at_temperature(t+0.01).unwrap().state(b,volume).unwrap();
            close(-(hi.free_energy_j-lo.free_energy_j)/0.02,actual.entropy_j_k,1e-8);
        } else {assert_eq!(actual.entropy_j_k,0.0); assert_eq!(actual.heat_capacity_j_k,0.0);}
        close(actual.internal_energy_j,psi+t*actual.entropy_j_k,1e-12);
        close(at.state(b,3.0*volume).unwrap().free_energy_j,3.0*psi,1e-12);
        assert_eq!(at.material().card_identity(),source.content_hash());
        let bound=law.minimum_heat_capacity_j_k(b,volume).unwrap();
        for k in 0..=30 {
            let sampled=law.at_temperature(200.0+10.0*k as f64).unwrap().state(b,volume).unwrap();
            assert!(sampled.heat_capacity_j_k>=bound);
        }
        assert!(law.at_temperature(501.0).is_err());
        assert!(at.state(-b,volume).is_err());
        assert!(at.state(b,0.0).is_err());
    }
}

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
}
fn card(name: &str, dims: Dims, value: PropertyValue, validity: ValidityDomain) -> MaterialCard {
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored affine field thermodynamics fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric or magnetic medium".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(name: &str, dims: Dims, knots: Vec<(f64, f64)>) -> MaterialCard {
    card(
        name,
        dims,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims,
        },
        ValidityDomain::unconstrained().with_quantity(
            "T",
            quantity(QuantityKind::AbsoluteTemperature),
            200.0,
            500.0,
        ),
    )
}
fn capacitor(source: &MaterialCard, p: &QueryPoint) -> ResolvedCapacitor {
    ResolvedCapacitor::resolve(
        source,
        p,
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
        0.001,
    )
    .unwrap()
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-100),
        "{a} != {b}"
    );
}

#[test]
fn g1_affine_dielectric_uses_helmholtz_entropy_and_internal_energy() {
    // epsilon(T) = 1e-11 + 1e-13 T [F/m], C=20 epsilon.
    let source = curve(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        vec![(200.0, 3e-11), (300.0, 4e-11), (500.0, 6e-11)],
    );
    let original = capacitor(&source, &point(300.0));
    let law = ThermalFieldStorage::capacitor(&source, &original, "T", [200.0, 500.0]).unwrap();
    let t = 350.0;
    let q = 2e-9;
    let at = law.at_temperature(t).unwrap();
    let state = at.state(q).unwrap();
    let c = 20.0 * (1e-11 + 1e-13 * t);
    let dc = 20.0 * 1e-13;
    let expected_f = q * q / (2.0 * c);
    let expected_s = q * q * dc / (2.0 * c * c);
    close(state.free_energy_j, expected_f, 1e-13);
    close(state.entropy_j_k, expected_s, 1e-13);
    close(state.internal_energy_j, expected_f + t * expected_s, 1e-13);
    close(
        state.heat_capacity_j_k,
        -t * q * q * dc * dc / (c * c * c),
        1e-13,
    );
    assert_eq!(at.state(-q).unwrap(), state);
    assert_eq!(at.state(0.0).unwrap().internal_energy_j, 0.0);
    assert!(state.heat_capacity_j_k < 0.0);
    match at.element() {
        ThermalFieldElement::Capacitor(native) => {
            close(native.capacitance_f(), c, 1e-13);
            assert_eq!(native.stored_energy_j(q).unwrap(), state.free_energy_j);
            assert_eq!(native.material().card_identity(), source.content_hash());
            assert_eq!(
                native.material().properties()[0].answer().receipt.selected,
                original.material().properties()[0]
                    .answer()
                    .receipt
                    .selected
            );
        }
        _ => panic!("capacitor changed kind"),
    }
    // Independent centered derivatives of source-resolved free and internal energy.
    let dt = 0.01;
    let lo = law.at_temperature(t - dt).unwrap().state(q).unwrap();
    let hi = law.at_temperature(t + dt).unwrap().state(q).unwrap();
    close(
        -(hi.free_energy_j - lo.free_energy_j) / (2.0 * dt),
        state.entropy_j_k,
        1e-8,
    );
    close(
        (hi.internal_energy_j - lo.internal_energy_j) / (2.0 * dt),
        state.heat_capacity_j_k,
        1e-8,
    );
    assert_eq!(law.material(), original.material());
    assert_eq!(law.temperature_range_k(), [200.0, 500.0]);
    assert_eq!(law.original_temperature_k(), 300.0);
    assert_eq!(law.card(), &source);
}

#[test]
fn g1_magnetic_storage_keeps_flux_linkage_and_reversible_thermal_sign() {
    // mu(T)=8e-4 - 1e-6 T, L = 25 mu. Negative slope gives negative field entropy.
    let source = curve(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        vec![(200.0, 6e-4), (350.0, 4.5e-4), (500.0, 3e-4)],
    );
    let original = ResolvedInductor::resolve(
        &source,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        0.001,
        0.1,
        50,
    )
    .unwrap();
    let law = ThermalFieldStorage::inductor(&source, &original, "T", [200.0, 500.0]).unwrap();
    let linkage = 0.01;
    let t = 360.0;
    let state = law.at_temperature(t).unwrap().state(linkage).unwrap();
    let l = 25.0 * (8e-4 - 1e-6 * t);
    let dl = -25e-6;
    close(state.free_energy_j, linkage * linkage / (2.0 * l), 1e-13);
    close(
        state.entropy_j_k,
        linkage * linkage * dl / (2.0 * l * l),
        1e-13,
    );
    close(
        state.internal_energy_j,
        state.free_energy_j + t * state.entropy_j_k,
        1e-13,
    );
    assert!(state.entropy_j_k < 0.0);
    assert!(state.heat_capacity_j_k < 0.0);
    // This is flux LINKAGE; the N^2 geometry is present exactly once.
    match law.at_temperature(t).unwrap().element() {
        ThermalFieldElement::Inductor(native) => close(native.inductance_h(), l, 1e-13),
        _ => panic!("inductor changed kind"),
    }
}

#[test]
fn g3_constant_source_has_no_field_entropy_or_hidden_caloric_exchange() {
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let source = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Scalar { value: 4e-11, dims },
        ValidityDomain::unconstrained(),
    );
    let original = capacitor(&source, &point(300.0));
    let law = ThermalFieldStorage::capacitor(&source, &original, "T", [200.0, 500.0]).unwrap();
    let a = law.at_temperature(200.0).unwrap().state(3e-9).unwrap();
    let b = law.at_temperature(500.0).unwrap().state(3e-9).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.entropy_j_k, 0.0);
    assert_eq!(a.heat_capacity_j_k, 0.0);
    assert_eq!(a.free_energy_j, a.internal_energy_j);
    assert_eq!(law.minimum_heat_capacity_j_k(1.0).unwrap(), 0.0);
}

#[test]
fn g3_capacity_bound_covers_every_supported_temperature_and_coordinate() {
    for values in [
        vec![(200.0, 3e-11), (500.0, 6e-11)],
        vec![(200.0, 6e-11), (500.0, 3e-11)],
    ] {
        let source = curve(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            values,
        );
        let original = capacitor(&source, &point(300.0));
        let law = ThermalFieldStorage::capacitor(&source, &original, "T", [200.0, 500.0]).unwrap();
        let qmax = 1e-6;
        let bound = law.minimum_heat_capacity_j_k(qmax).unwrap();
        assert!(bound < 0.0);
        for index in 0..=30 {
            let at = law.at_temperature(200.0 + 10.0 * f64::from(index)).unwrap();
            for fraction in [-1.0, -0.5, 0.0, 0.25, 1.0] {
                assert!(at.state(fraction * qmax).unwrap().heat_capacity_j_k >= bound);
            }
        }
        assert_eq!(law.minimum_heat_capacity_j_k(0.0).unwrap(), 0.0);
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(law.minimum_heat_capacity_j_k(bad).is_err());
        }
    }
}

#[test]
fn g3_non_affine_and_incomplete_support_refuse_even_with_equal_endpoints() {
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let kinked = curve(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        vec![(200.0, 4e-11), (350.0, 8e-11), (500.0, 4e-11)],
    );
    let original = capacitor(&kinked, &point(300.0));
    assert!(ThermalFieldStorage::capacitor(&kinked, &original, "T", [200.0, 500.0]).is_err());
    let affine = curve(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        vec![(250.0, 4e-11), (450.0, 8e-11)],
    );
    let original = capacitor(&affine, &point(300.0));
    assert!(ThermalFieldStorage::capacitor(&affine, &original, "T", [200.0, 500.0]).is_err());
    for range in [
        [300.0, 300.0],
        [500.0, 200.0],
        [0.0, 500.0],
        [f64::NAN, 500.0],
        [200.0, f64::INFINITY],
        [350.0, 400.0],
    ] {
        assert!(ThermalFieldStorage::capacitor(&affine, &original, "T", range).is_err());
    }
    let admitted = ThermalFieldStorage::capacitor(&affine, &original, "T", [250.0, 450.0]).unwrap();
    for bad in [249.0, 451.0, f64::NAN, f64::INFINITY] {
        assert!(admitted.at_temperature(bad).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(admitted.at_temperature(300.0).unwrap().state(bad).is_err());
    }
    // Same material label does not authorize a different source card.
    assert!(ThermalFieldStorage::capacitor(&kinked, &original, "T", [250.0, 450.0]).is_err());
    let untyped = capacitor(&affine, &QueryPoint::new().with("T", 300.0).unwrap());
    assert!(ThermalFieldStorage::capacitor(&affine, &untyped, "T", [250.0, 450.0]).is_err());
}

#[test]
fn g1_affine_surface_keeps_fixed_pressure_and_rejects_a_second_thermal_owner() {
    let pressure = quantity(QuantityKind::Pressure);
    let temperature = quantity(QuantityKind::AbsoluteTemperature);
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let source = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Surface2D {
            axes: ["pressure".into(), "T".into()],
            axis_dims: [pressure.dims(), temperature.dims()],
            samples: vec![
                ([1e5, 200.0], 3e-11),
                ([3e5, 200.0], 5e-11),
                ([1e5, 500.0], 6e-11),
                ([3e5, 500.0], 8e-11),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims,
        },
        ValidityDomain::unconstrained()
            .with_quantity("pressure", pressure, 1e5, 3e5)
            .with_quantity("T", temperature, 200.0, 500.0),
    );
    let context = point(300.0)
        .with_quantity("pressure", pressure, 2e5)
        .unwrap();
    let original = capacitor(&source, &context);
    let law = ThermalFieldStorage::capacitor(&source, &original, "T", [200.0, 500.0]).unwrap();
    let at = law.at_temperature(400.0).unwrap();
    assert_eq!(
        at.element().material().query().axes().get("pressure"),
        Some(&2e5)
    );
    close(at.element().coefficient(), 20.0 * 6e-11, 1e-13);
    close(law.coefficient_slope_per_k(), 20.0 * 1e-13, 1e-13);

    let second = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        dims,
        PropertyValue::Curve {
            abscissa: "other_temperature".into(),
            abscissa_dims: temperature.dims(),
            knots: vec![(200.0, 3e-11), (500.0, 6e-11)],
            dims,
        },
        ValidityDomain::unconstrained().with_quantity(
            "other_temperature",
            temperature,
            200.0,
            500.0,
        ),
    );
    let context = point(300.0)
        .with_quantity("other_temperature", temperature, 300.0)
        .unwrap();
    let original = capacitor(&second, &context);
    assert!(ThermalFieldStorage::capacitor(&second, &original, "T", [200.0, 500.0]).is_err());
}

