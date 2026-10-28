//! G0/G1/G4 uniform-field reductions from synthetic data, not measurements.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{
        CompositionBasis, FrequencyConvention, MoistureBasis, QuantityKind, SemanticType,
        StrainBasis, StrainComponent, ValueForm,
    },
};

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

mod nonlinear_magnetics {
    use super::*;
    use fs_material::field_storage::{
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    };

    fn point() -> QueryPoint {
        QueryPoint::new()
            .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
            .unwrap()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap()
    }

    fn curve(knots: Vec<(f64, f64)>) -> MaterialCard {
        card(
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots,
                dims: FIELD_STRENGTH_DIMS,
            },
            ValidityDomain::unconstrained(),
        )
    }

    fn resolve(source: &MaterialCard, maximum: f64) -> ReversibleMagneticMaterial {
        ReversibleMagneticMaterial::resolve(
            source,
            &point(),
            "B",
            maximum,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    fn check_source_and_inverse(source: &MaterialCard, maximum: f64) {
        let law = resolve(source, maximum);
        for i in 0..=100 {
            let b = maximum * f64::from(i) / 100.0;
            let response = law.evaluate(b).unwrap();
            let source_point = point().with("B", b).unwrap();
            let canonical = law
                .material()
                .at_query_point(source, &source_point)
                .unwrap();
            near(
                response.field_strength_a_m,
                canonical
                    .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
                    .unwrap()
                    .value_si(),
            );
            near(
                law.flux_density_at_field_strength(response.field_strength_a_m)
                    .unwrap(),
                b,
            );
            near(
                response.energy_density_j_m3 + response.coenergy_density_j_m3,
                b * response.field_strength_a_m,
            );
        }
        // Exact original knots keep their original canonical value, and the
        // inverse uses these same knots rather than another source fit.
        for (b, answer) in &law.source_path().knots {
            let h = answer.evidence.value.value;
            assert_eq!(law.evaluate(*b).unwrap().field_strength_a_m, h);
            assert_eq!(law.flux_density_at_field_strength(h).unwrap(), *b);
        }
    }

    #[test]
    fn compiled_saturation_curve_retains_source_inverse_and_legendre_derivatives() {
        let source = curve(vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]);
        check_source_and_inverse(&source, 4.0);
        let law = resolve(&source, 4.0);
        let response = law.evaluate(1.25).unwrap();
        near(response.energy_density_j_m3, 0.875);
        near(response.coenergy_density_j_m3, 1.625);
        near(response.differential_reluctivity_m_h, 4.0);
        let db = 1e-5;
        let dh = 4.0 * db;
        let derivative = (law.evaluate(1.25 + db).unwrap().energy_density_j_m3
            - law.evaluate(1.25 - db).unwrap().energy_density_j_m3)
            / (2.0 * db);
        assert!((derivative - 2.0).abs() < 1e-8);
        let dual = |h| {
            law.evaluate(law.flux_density_at_field_strength(h).unwrap())
                .unwrap()
                .coenergy_density_j_m3
        };
        assert!(((dual(2.0 + dh) - dual(2.0 - dh)) / (2.0 * dh) - 1.25).abs() < 1e-8);
        near(
            law.bregman_density_j_m3([0.5, 0.0], [2.0, 0.0]).unwrap(),
            4.125,
        );
        near(
            law.bregman_density_j_m3([2.0, 0.0], [0.5, 0.0]).unwrap(),
            2.625,
        );
        near(
            law.bregman_density_j_m3([0.0, 2.0], [2.0, 0.0]).unwrap(),
            10.0,
        );
        assert_eq!(
            law.bregman_density_j_m3([1.25, 0.0], [1.25, 0.0]).unwrap(),
            0.0
        );
    }

    #[test]
    fn compiled_fixed_temperature_surface_slice_matches_original_source() {
        let source = card(
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            PropertyValue::Surface2D {
                axes: ["B".into(), "T".into()],
                axis_dims: [FLUX_DENSITY_DIMS, Temperature::DIMS],
                samples: vec![
                    ([0.0, 200.0], 0.0),
                    ([1.0, 200.0], 2.0),
                    ([2.0, 200.0], 8.0),
                    ([0.0, 400.0], 0.0),
                    ([1.0, 400.0], 4.0),
                    ([2.0, 400.0], 12.0),
                ],
                triangles: vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
                dims: FIELD_STRENGTH_DIMS,
            },
            ValidityDomain::unconstrained(),
        );
        check_source_and_inverse(&source, 2.0);
        let law = resolve(&source, 2.0);
        assert!(law.source_path().knots.len() > 3);
        assert_eq!(law.material().card_identity(), source.content_hash());
    }

    #[test]
    fn tiny_coenergy_survives_subtraction_cancellation_in_high_susceptibility_span() {
        let source = curve(vec![(0.0, 0.0), (1e-30, 1.0), (1.0, 1.0 + f64::EPSILON)]);
        let law = resolve(&source, 1.0);
        let response = law.evaluate(0.5).unwrap();
        // The original B*H-psi subtraction rounds to zero here. The integral
        // of B dH remains positive and analytically known for this source law.
        assert_eq!(
            0.5 * response.field_strength_a_m - response.energy_density_j_m3,
            0.0
        );
        near(
            response.coenergy_density_j_m3,
            0.125 * f64::EPSILON + 0.5e-30,
        );
        assert!(response.differential_reluctivity_m_h > 0.0);
    }

    #[test]
    fn inverse_and_compiled_evaluation_keep_original_field_support() {
        let source = curve(vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]);
        let law = resolve(&source, 2.0);
        assert_eq!(law.maximum_field_strength_a_m(), 5.0);
        assert_eq!(law.flux_density_at_field_strength(0.0).unwrap(), 0.0);
        assert_eq!(law.flux_density_at_field_strength(5.0).unwrap(), 2.0);
        for bad in [-1.0, f64::NAN, f64::INFINITY, 5.0_f64.next_up()] {
            assert!(law.flux_density_at_field_strength(bad).is_err());
        }
        for bad in [-1.0, f64::NAN, f64::INFINITY, 2.0_f64.next_up()] {
            assert!(law.evaluate(bad).is_err());
        }
    }
}

mod permanent_recoil_refresh {
    use super::*;
    use fs_material::{
        conductor::{ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY},
        field_storage::permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
        state_point::{
            ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
        },
    };

    const RECOIL_PROPERTIES: [(&str, Dims, f64); 5] = [
        (
            RECOIL_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.5e-6,
        ),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, 0.8),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            -1e6,
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            1e6,
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            1e5,
        ),
    ];

    fn point(temperature_k: f64) -> QueryPoint {
        QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                temperature_k,
            )
            .unwrap()
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                1e5,
            )
            .unwrap()
    }

    fn source(temperature_dependent_recoil: Option<&str>) -> MaterialCard {
        let mut set = ClaimSet::new();
        for (name, dims, value) in RECOIL_PROPERTIES.into_iter().chain([(
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            2e6,
        )]) {
            let is_recoil_curve = temperature_dependent_recoil == Some(name);
            let is_curve = is_recoil_curve || name == ELECTRICAL_CONDUCTIVITY_PROPERTY;
            let payload = if is_curve {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: if is_recoil_curve {
                        vec![(300.0, value), (350.0, 1.1 * value), (400.0, value)]
                    } else {
                        vec![(300.0, value), (400.0, 0.5 * value)]
                    },
                    dims,
                }
            } else {
                PropertyValue::Scalar { value, dims }
            };
            set.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: payload,
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                interpolation: if is_curve {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured conducting recoil material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic conducting magnet".into(),
                phase: "solid".into(),
                process: "fixed-recoil refresh fixture".into(),
                revision: 0,
            },
            set,
            vec![],
        )
        .unwrap()
    }

    fn recoil(source: &MaterialCard, include_conductor: bool) -> ResolvedRecoilMaterial {
        let mut properties = RECOIL_PROPERTIES.to_vec();
        if include_conductor {
            properties.push((
                ELECTRICAL_CONDUCTIVITY_PROPERTY,
                ELECTRICAL_CONDUCTIVITY_DIMS,
                2e6,
            ));
        }
        let requirements: Vec<_> = properties
            .into_iter()
            .map(|(name, dims, _)| {
                ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
            })
            .collect();
        ResolvedRecoilMaterial::try_from_material(
            resolve_material_state_point(
                source,
                &point(300.0),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn fixed_recoil_admits_native_conductivity_temperature_refresh() {
        let source = source(None);
        let law = recoil(&source, true);
        let warm = law
            .material()
            .at_query_point(&source, &point(400.0))
            .unwrap();
        assert_eq!(warm.properties().len(), 6);
        near(
            warm.property(ELECTRICAL_CONDUCTIVITY_PROPERTY)
                .unwrap()
                .value_si(),
            1e6,
        );
        law.admit_fixed_material_refresh(&warm).unwrap();
        for (name, _, _) in RECOIL_PROPERTIES {
            assert_eq!(
                law.material().property(name).unwrap().value_si().to_bits(),
                warm.property(name).unwrap().value_si().to_bits()
            );
        }
        assert_eq!(
            law.material().query_point(),
            &[("T".into(), 300.0), ("pressure".into(), 1e5)]
        );
    }

    #[test]
    fn every_temperature_dependent_recoil_property_refuses_even_equal_endpoints() {
        for (name, _, _) in RECOIL_PROPERTIES {
            let source = source(Some(name));
            let law = recoil(&source, true);
            law.admit_fixed_material_refresh(law.material()).unwrap();
            let warm = law
                .material()
                .at_query_point(&source, &point(400.0))
                .unwrap();
            assert_eq!(
                law.material().property(name).unwrap().value_si().to_bits(),
                warm.property(name).unwrap().value_si().to_bits()
            );
            assert!(law.admit_fixed_material_refresh(&warm).is_err(), "{name}");
            let midway = law
                .material()
                .at_query_point(&source, &point(350.0))
                .unwrap();
            assert!(law.admit_fixed_material_refresh(&midway).is_err(), "{name}");
        }
    }

    #[test]
    fn recoil_refresh_preserves_complete_bundle_card_and_nonthermal_context() {
        let original = source(None);
        let law = recoil(&original, true);
        let reduced = recoil(&original, false);
        assert!(
            law.admit_fixed_material_refresh(reduced.material())
                .is_err()
        );
        let other_card = source(Some(REMANENT_FLUX_DENSITY_PROPERTY));
        let substituted = recoil(&other_card, true);
        assert!(
            law.admit_fixed_material_refresh(substituted.material())
                .is_err()
        );
        let changed_pressure = law
            .material()
            .at_query_point(&original, &point(300.0).with("pressure", 2e5).unwrap())
            .unwrap();
        assert!(law.admit_fixed_material_refresh(&changed_pressure).is_err());
    }
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

fn static_quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

#[test]
fn g1_typed_fixed_external_curves_set_magnetic_storage_without_losing_sources() {
    for kind in [
        QuantityKind::Pressure,
        QuantityKind::Composition(CompositionBasis::MassFraction),
        QuantityKind::Moisture(MoistureBasis::RelativeHumidity),
        QuantityKind::MassConcentration,
        QuantityKind::AmountConcentration,
    ] {
        let quantity = static_quantity(kind);
        let source = card(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            PropertyValue::Curve {
                abscissa: "environment".into(),
                abscissa_dims: quantity.dims(),
                knots: vec![(0.25, 2e-4), (0.75, 4e-4)],
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
            ValidityDomain::unconstrained().with_quantity("environment", quantity, 0.25, 0.75),
        );
        let point = QueryPoint::new()
            .with_quantity("environment", quantity, 0.5)
            .unwrap();
        let first = ResolvedInductor::resolve(
            &source,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.001,
            0.1,
            50,
        )
        .unwrap();
        near(first.inductance_h(), 0.0075);
        near(
            linear_magnetic_permeability(first.material()).unwrap(),
            3e-4,
        );
        assert_eq!(first.material().card_identity(), source.content_hash());
        assert_eq!(
            first.material().axis_quantities().get("environment"),
            Some(&quantity)
        );
        assert_eq!(first.material().properties()[0].source_axes().len(), 1);
        assert!(first.material().properties()[0].is_static_at_fixed_external_conditions());

        // Another fixed physical context is a new snapshot, not a live zero-work update.
        let next_point = point.clone().with("environment", 0.625).unwrap();
        let refreshed = first
            .material()
            .at_query_point(&source, &next_point)
            .unwrap();
        let second = ResolvedInductor::try_from_material(refreshed, 0.001, 0.1, 50).unwrap();
        near(second.inductance_h(), 0.00875);
        near(first.inductance_h(), 0.0075);
        assert_eq!(
            first.material().properties()[0].answer().receipt.selected,
            second.material().properties()[0].answer().receipt.selected
        );
        assert!(
            ResolvedInductor::resolve(
                &source,
                &point.clone().with("environment", 0.8).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                0.001,
                0.1,
                50,
            )
            .is_err()
        );
    }
}

#[test]
fn g1_fixed_pressure_composition_and_temperature_surfaces_set_electric_storage() {
    let pressure = static_quantity(QuantityKind::Pressure);
    for (kind, low, high) in [
        (
            QuantityKind::Composition(CompositionBasis::MoleFraction),
            0.0,
            1.0,
        ),
        (QuantityKind::AbsoluteTemperature, 300.0, 400.0),
    ] {
        let context = static_quantity(kind);
        let source = card(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            PropertyValue::Surface2D {
                axes: ["pressure".into(), "context".into()],
                axis_dims: [pressure.dims(), context.dims()],
                samples: vec![
                    ([1e5, low], 4e-11),
                    ([3e5, low], 8e-11),
                    ([1e5, high], 6e-11),
                    ([3e5, high], 10e-11),
                ],
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            },
            ValidityDomain::unconstrained()
                .with_quantity("pressure", pressure, 1e5, 3e5)
                .with_quantity("context", context, low, high),
        );
        let point = QueryPoint::new()
            .with_quantity("pressure", pressure, 2e5)
            .unwrap()
            .with_quantity("context", context, f64::midpoint(low, high))
            .unwrap();
        let first = ResolvedCapacitor::resolve(
            &source,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.01,
            0.001,
        )
        .unwrap();
        // Independent affine field coefficient: epsilon=4e-11+4e-11*p_hat+2e-11*c_hat.
        near(first.capacitance_f(), 7e-10);
        near(first.stored_energy_j(2e-9).unwrap(), 0.5 * 4e-18 / 7e-10);
        assert_eq!(first.material().properties()[0].source_axes().len(), 2);
        let changed = first
            .material()
            .at_query_point(&source, &point.clone().with("pressure", 2.5e5).unwrap())
            .unwrap();
        let second = ResolvedCapacitor::try_from_material(changed, 0.01, 0.001).unwrap();
        near(second.capacitance_f(), 8e-10);
        assert_eq!(
            first.material().card_identity(),
            second.material().card_identity()
        );
        assert_eq!(
            first.material().properties()[0].answer().receipt.selected,
            second.material().properties()[0].answer().receipt.selected
        );
    }
}

#[test]
fn g4_dynamic_or_untyped_dependencies_and_applicability_bounds_are_not_fixed_context() {
    let cases = [
        (
            "frequency",
            static_quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
        ),
        ("stress", static_quantity(QuantityKind::Stress)),
        (
            "strain",
            static_quantity(QuantityKind::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Normal,
            }),
        ),
        ("B", QuantitySpec::dimensional(Dims([0, 1, -2, 0, -1, 0]))),
        ("E", QuantitySpec::dimensional(Dims([1, 1, -3, 0, -1, 0]))),
        (
            "strain_rate",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
        ),
        (
            "pressure",
            QuantitySpec::dimensional(Dims([-1, 1, -2, 0, 0, 0])),
        ),
        ("composition", QuantitySpec::dimensional(Dims([0; 6]))),
    ];
    for (axis, quantity) in cases {
        for varying in [false, true] {
            let dims = ELECTRIC_PERMITTIVITY_DIMS;
            let value = if varying {
                PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: quantity.dims(),
                    knots: vec![(0.25, 4e-11), (0.75, 8e-11)],
                    dims,
                }
            } else {
                PropertyValue::Scalar { value: 4e-11, dims }
            };
            let source = card(
                ELECTRIC_PERMITTIVITY_PROPERTY,
                dims,
                value,
                ValidityDomain::unconstrained().with_quantity(axis, quantity, 0.25, 0.75),
            );
            let result = ResolvedCapacitor::resolve(
                &source,
                &QueryPoint::new()
                    .with_quantity(axis, quantity, 0.5)
                    .unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                1.0,
                1.0,
            );
            assert!(
                matches!(result, Err(FieldStorageError::Invalid(_))),
                "unsupported {axis} must be refused by field-storage admission; varying={varying}"
            );
        }
    }
}
