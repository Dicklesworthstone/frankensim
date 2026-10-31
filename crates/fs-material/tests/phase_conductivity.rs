//! Source paths and analytic latent/sensible derivatives on the original chart.
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::conduction::phase::{PhaseConductivityError, PreparedPhaseConductivity};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionLaw;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

fn claim(value: PropertyValue, fraction_bounds: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("k", KD),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                800.0,
            )
            .with_quantity(
                "f",
                QuantitySpec::dimensional(Dims::NONE),
                fraction_bounds[0],
                fraction_bounds[1],
            ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: fs_matdb::UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured constitutive fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn phase_values(knots: &[(f64, f64)]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "f".into(),
        abscissa_dims: Dims::NONE,
        knots: knots.to_vec(),
        dims: KD,
    }
}
fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic phase path".into(),
            phase: "declared equilibrium chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn curve(identity: ContentHash, states: &[(f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        identity,
        states
            .iter()
            .map(|&(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 1000.0,
            })
            .collect(),
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    curve(
        card.content_hash(),
        &[
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (50.0, 400.0, 0.5),
            (100.0, 400.0, 1.0),
            (500.0, 800.0, 1.0),
        ],
    )
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
}
fn prepare<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    range: [f64; 2],
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<PreparedPhaseConductivity<'a>, PhaseConductivityError> {
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context(),
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        range,
        checkpoint,
    )
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn every_source_and_caloric_knot_survives_with_correct_one_sided_latent_tangents() {
    let card = card(vec![claim(
        phase_values(&[(0.0, 2.0), (0.25, 4.0), (0.75, 1.0), (1.0, 3.0)]),
        [0.0, 1.0],
    )]);
    let curve = caloric(&card);
    let path = prepare(&card, &curve, [-100.0, 500.0], || Ok(())).unwrap();
    assert_eq!(
        path.knots().iter().map(|k| k[0]).collect::<Vec<_>>(),
        vec![-100.0, 0.0, 25.0, 50.0, 75.0, 100.0, 500.0]
    );
    for (h, k, dk, dt) in [
        (-100.0, 2.0, 0.0, 1.0),
        (-20.0, 2.0, 0.0, 1.0),
        (0.0, 2.0, 0.08, 0.0),
        (12.5, 3.0, 0.08, 0.0),
        (25.0, 4.0, -0.06, 0.0),
        (40.0, 3.1, -0.06, 0.0),
        (75.0, 1.0, 0.08, 0.0),
        (100.0, 3.0, 0.0, 1.0),
        (500.0, 3.0, 0.0, 1.0),
    ] {
        let sample = path.state_at_specific_enthalpy(h).unwrap();
        close(sample.conductivity_w_m_k(), k);
        close(sample.conductivity_enthalpy_derivative_kg_w_m_k_j(), dk);
        assert_eq!(sample.temperature_enthalpy_derivative_kg_k_j(), dt);
        let original = curve.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(
            sample.temperature_k().to_bits(),
            original.temperature_k().to_bits()
        );
        assert_eq!(
            sample.liquid_mass_fraction().to_bits(),
            original.liquid_mass_fraction().to_bits()
        );
        assert_eq!(sample.phase(), original.phase());
        assert_eq!(
            ConductivityLaw::conductivity_w_m_k(&path, h).unwrap(),
            sample.conductivity_w_m_k()
        );
        let native = ImplicitConductionLaw::state_at_specific_enthalpy(&path, h).unwrap();
        assert_eq!(
            native.conductivity_slope_kg_w_m_k_j,
            sample.conductivity_enthalpy_derivative_kg_w_m_k_j()
        );
        assert_eq!(native.temperature_slope_kg_k_j, dt);
        assert!(sample.source_line_index() < path.source_lines().len());
    }
    let clipped = prepare(&card, &curve, [10.0, 100.0], || Ok(())).unwrap();
    assert_eq!(clipped.enthalpy_bounds_j_kg(), [10.0, 100.0]);
    close(
        clipped
            .state_at_specific_enthalpy(100.0)
            .unwrap()
            .conductivity_enthalpy_derivative_kg_w_m_k_j(),
        0.08,
    );
    assert_eq!(
        clipped
            .state_at_specific_enthalpy(100.0)
            .unwrap()
            .temperature_enthalpy_derivative_kg_k_j(),
        0.0
    );
    assert!(
        clipped
            .state_at_specific_enthalpy(100.0_f64.next_up())
            .is_err()
    );
    assert_ne!(path.identity(), clipped.identity());
    assert!(std::ptr::eq(path.caloric_curve(), &curve));
    assert!(std::ptr::addr_eq(
        ImplicitConductionLaw::caloric_law(&path) as *const dyn ThermalLaw,
        &curve as *const EquilibriumEnthalpyPhaseCurve
    ));
    assert!(std::ptr::addr_eq(
        ImplicitConductionLaw::conductivity_law(&path) as *const dyn ConductivityLaw,
        &path as *const PreparedPhaseConductivity<'_>
    ));
    assert!(
        path.source_lines()
            .iter()
            .flat_map(|line| &line.knots)
            .all(|(_, answer)| answer.receipt.selected == path.selected_claim())
    );
    let observable = ImplicitConductionLaw::time_observables(&path, 60.0).unwrap();
    assert_eq!(observable.temperature_k, 400.0);
    close(observable.phase_mass_fraction.unwrap(), 0.6);
}

#[test]
fn one_two_coordinate_source_drives_sensible_and_latent_transport_without_a_mixing_rule() {
    // k=2+0.01(T-300)+3f, supplied as two explicit affine source cells.
    // The latent path T=400 crosses their diagonal at f=.2 (h=20).
    let card = card(vec![claim(
        PropertyValue::Surface2D {
            axes: ["T".into(), "f".into()],
            axis_dims: [Temperature::DIMS, Dims::NONE],
            samples: vec![
                ([300.0, 0.0], 2.0),
                ([800.0, 0.0], 7.0),
                ([800.0, 1.0], 10.0),
                ([300.0, 1.0], 5.0),
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3]],
            dims: KD,
        },
        [0.0, 1.0],
    )]);
    let curve = caloric(&card);
    let path = prepare(&card, &curve, [-100.0, 500.0], || Ok(())).unwrap();
    assert!(path.knots().iter().any(|k| (k[0] - 20.0).abs() < 1e-12));
    for h in [-80.0, -1.0, 10.0, 20.0, 50.0, 90.0, 120.0, 300.0, 500.0] {
        let sample = path.state_at_specific_enthalpy(h).unwrap();
        let expected =
            2.0 + 0.01 * (sample.temperature_k() - 300.0) + 3.0 * sample.liquid_mass_fraction();
        close(sample.conductivity_w_m_k(), expected);
        close(
            sample.conductivity_enthalpy_derivative_kg_w_m_k_j(),
            if h < 0.0 || h >= 100.0 { 0.01 } else { 0.03 },
        );
    }
    // Original numeric and identity-bearing caloric queries agree exactly.
    for h in [-100.0, 0.0, 37.0, 40.0, 100.0, 299.0, 500.0] {
        let values = curve.values_at_specific_enthalpy(h).unwrap();
        let state = curve.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(
            values.temperature_k.to_bits(),
            state.temperature_k().to_bits()
        );
        assert_eq!(
            values.liquid_mass_fraction.to_bits(),
            state.liquid_mass_fraction().to_bits()
        );
        assert_eq!(
            values.bulk_density_kg_m3.to_bits(),
            state.bulk_density_kg_m3().to_bits()
        );
    }
}

#[test]
fn source_changes_holes_negative_interior_and_unowned_phase_paths_refuse() {
    let split = card(vec![
        claim(phase_values(&[(0.0, 2.0), (1.0, 3.0)]), [0.0, 0.5]),
        claim(phase_values(&[(0.0, 4.0), (1.0, 5.0)]), [0.5, 1.0]),
    ]);
    let chart = caloric(&split);
    assert!(matches!(
        prepare(&split, &chart, [-100.0, 500.0], || Ok(())),
        Err(PhaseConductivityError::Source(_))
    ));
    let negative = card(vec![claim(
        phase_values(&[(0.0, 2.0), (0.5, -1.0), (1.0, 3.0)]),
        [0.0, 1.0],
    )]);
    assert!(matches!(
        prepare(&negative, &caloric(&negative), [-100.0, 500.0], || Ok(())),
        Err(PhaseConductivityError::Invalid(_))
    ));
    let ordinary = card(vec![claim(
        phase_values(&[(0.0, 2.0), (1.0, 3.0)]),
        [0.0, 1.0],
    )]);
    let mushy = curve(
        ordinary.content_hash(),
        &[(0.0, 300.0, 0.0), (100.0, 800.0, 1.0)],
    );
    assert!(matches!(
        prepare(&ordinary, &mushy, [0.0, 100.0], || Ok(())),
        Err(PhaseConductivityError::Invalid(_))
    ));
    let mixed_sensible = curve(
        ordinary.content_hash(),
        &[
            (0.0, 300.0, 0.0),
            (25.0, 400.0, 0.5),
            (50.0, 500.0, 0.5),
            (100.0, 800.0, 1.0),
        ],
    );
    assert!(matches!(
        prepare(&ordinary, &mixed_sensible, [25.0, 50.0], || Ok(())),
        Err(PhaseConductivityError::Invalid(_))
    ));
    let ordinary_chart = caloric(&ordinary);
    let missing_axis = PreparedPhaseConductivity::new(
        &ordinary_chart,
        ThermalMaterialPoint {
            material: &ordinary,
            context: &context(),
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [0.0, 100.0],
        || Ok(()),
    );
    assert!(matches!(
        missing_axis,
        Err(PhaseConductivityError::Invalid(_))
    ));
    let foreign = curve(
        ContentHash([0xaa; 32]),
        &[
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 400.0, 1.0),
            (500.0, 800.0, 1.0),
        ],
    );
    assert!(matches!(
        prepare(&ordinary, &foreign, [-100.0, 500.0], || Ok(())),
        Err(PhaseConductivityError::Material(_))
    ));
}

#[test]
fn explicit_insulation_domains_and_late_cancellation_preserve_the_original_law() {
    let zero = card(vec![claim(
        phase_values(&[(0.0, 0.0), (1.0, 0.0)]),
        [0.0, 1.0],
    )]);
    let curve = caloric(&zero);
    let mut calls = 0;
    let path = prepare(&zero, &curve, [-100.0, 500.0], || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    for h in [-100.0, 0.0, 50.0, 100.0, 500.0] {
        let sample = path.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(sample.conductivity_w_m_k(), 0.0);
        assert_eq!(sample.conductivity_enthalpy_derivative_kg_w_m_k_j(), 0.0);
    }
    let mut retry_calls = 0;
    assert!(matches!(
        prepare(&zero, &curve, [-100.0, 500.0], || {
            retry_calls += 1;
            if retry_calls == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PhaseConductivityError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(
        prepare(&zero, &curve, [-100.0, 500.0], || Ok(()))
            .unwrap()
            .identity(),
        path.identity()
    );
    for h in [-101.0, 501.0, f64::NAN] {
        assert!(matches!(
            path.state_at_specific_enthalpy(h),
            Err(PhaseConductivityError::OutsideEnthalpyDomain { .. })
        ));
        assert_eq!(
            ConductivityLaw::conductivity_w_m_k(&path, h),
            Err(DepositionError::MaterialDomain)
        );
    }
    assert_eq!(
        path.state_at_specific_enthalpy(50.0).unwrap().phase(),
        SolidLiquidPhase::SolidLiquid
    );
}
