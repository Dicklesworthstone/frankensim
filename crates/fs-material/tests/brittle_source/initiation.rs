//! G0/G1/G3/G4: original sources reach native tensile onset and fracture history.
use super::*;
use fs_material::Voigt;
use fs_material::brittle::initiation::{TensileInitiationError, TensileInitiationSettings};
use fs_material::brittle::source::initiation::ResolvedTensileInitiationError;
use fs_material::tensor::contract;

const E: f64 = 1e9;
const NU: f64 = 0.2;
const FT: f64 = 1e6;

fn resolution() -> TensileInitiationSettings {
    TensileInitiationSettings {
        maximum_strain: 0.05,
        relative_stress_tolerance: 1e-10,
        fraction_tolerance: 1e-12,
        maximum_bisections: 64,
    }
}

fn uniaxial(axis: usize, stress: f64) -> Voigt {
    let mut strain = [0.0; 6];
    strain[..3].fill(-NU * stress / E);
    strain[axis] = stress / E;
    strain
}

#[test]
fn discovered_plane_retains_complete_source_and_seeds_native_fracture_without_new_work() {
    let mut claims = properties();
    claims.push(scalar(DENSITY_PROPERTY, Density::DIMS, 2300.0));
    let source = card(claims);
    let mut requirements = brittle_crack_band_requirements().unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let resolved = ResolvedBrittleCrackBand::try_from_material(bundle, settings(0.01)).unwrap();
    let detector = resolved.tensile_initiator(resolution()).unwrap();
    assert!(core::ptr::eq(detector.source(), &resolved));
    assert!(core::ptr::eq(
        detector.native().elastic(),
        resolved.law().elastic()
    ));

    // The original source band is x-normal; the actual tensile load is along y.
    let start = uniaxial(1, 0.2 * FT);
    let end = uniaxial(1, 4.0 * FT);
    let event = detector.between(&start, &end, None).unwrap().unwrap();
    assert!(core::ptr::eq(event.source(), &resolved));
    let [low, high] = event.native().fraction_bracket();
    let expected = 0.8 / 3.8;
    assert!(low <= expected && expected <= high);
    assert!(event.native().normal()[1].abs() > 1.0 - 1e-14);

    let receiving = event.resolve_band(0.008).unwrap();
    assert_eq!(receiving.material(), resolved.material());
    assert_eq!(receiving.material().properties().len(), 8);
    assert_eq!(receiving.law().parameters(), resolved.law().parameters());
    assert_eq!(receiving.law().settings().maximum_strain, 0.05);
    assert_eq!(receiving.law().settings().energy_tolerance_j_m3, 1e-6);
    assert_eq!(receiving.law().band_width_m(), 0.008);
    near(
        receiving
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        2300.0,
    );
    let seed = event.seed(&receiving).unwrap();
    let stress = resolved.law().elastic().stress(&start, &());
    let existing_storage = 0.5 * contract(&stress, &start);
    near(seed.elastic_storage_j_m3, existing_storage);
    near(seed.state.stored_energy_j_m3(), existing_storage);
    assert_eq!(seed.state.strain(), &start);
    assert_eq!(seed.state.damage(), 0.0);
    assert_eq!(seed.state.dissipated_j_m3(), 0.0);
    let step = receiving.law().trial(&end, &seed.state).unwrap();
    assert!(step.state.damage() > 0.0);
    near(
        step.energy.work,
        step.state.stored_energy_j_m3() - existing_storage + step.state.dissipated_j_m3(),
    );
    // The original source remains available; geometry selection did not mutate it.
    assert_eq!(resolved.law().normal(), [1.0, 0.0, 0.0]);
    assert_eq!(resolved.law().band_width_m(), 0.01);
}

#[test]
fn equal_numeric_coefficients_do_not_authorize_another_source_or_context() {
    let source = card(properties());
    let resolved = resolve(&source, 0.01);
    let event = resolved
        .tensile_initiator(resolution())
        .unwrap()
        .between(&uniaxial(0, 0.2 * FT), &uniaxial(0, 2.0 * FT), None)
        .unwrap()
        .unwrap();
    let mut other_properties = properties();
    other_properties[6].provenance.source = "another specimen with equal numbers".into();
    let other = resolve(&card(other_properties), 0.01);
    assert_eq!(other.law().parameters(), resolved.law().parameters());
    // The unsourced numerical law can match; the sourced bridge refuses it.
    assert!(event.native().seed(other.law()).is_ok());
    assert!(matches!(
        event.seed(&other),
        Err(ResolvedTensileInitiationError::SourceMismatch)
    ));
    let at_temperature = |temperature| {
        resolve_brittle_crack_band(
            &source,
            &QueryPoint::new()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    temperature,
                )
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .unwrap()
    };
    let cold = at_temperature(300.0);
    let warm = at_temperature(400.0);
    let thermal_event = cold
        .tensile_initiator(resolution())
        .unwrap()
        .between(&uniaxial(0, 0.2 * FT), &uniaxial(0, 2.0 * FT), None)
        .unwrap()
        .unwrap();
    assert_eq!(cold.law().parameters(), warm.law().parameters());
    assert_eq!(
        cold.material().card_identity(),
        warm.material().card_identity()
    );
    assert!(thermal_event.native().seed(warm.law()).is_ok());
    assert!(matches!(
        thermal_event.seed(&warm),
        Err(ResolvedTensileInitiationError::SourceMismatch)
    ));
}

#[test]
fn original_temperature_dependent_strength_changes_the_actual_onset() {
    let mut claims = properties();
    claims[2].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, FT), (600.0, 0.5 * FT)],
        dims: Pressure::DIMS,
    };
    claims[2].validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        600.0,
    );
    claims[2].interpolation = InterpolationPolicy::LinearInside;
    let source = card(claims);
    for (temperature, expected_fraction) in [(300.0, 0.5), (600.0, 0.25)] {
        let resolved = resolve_brittle_crack_band(
            &source,
            &QueryPoint::new()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    temperature,
                )
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .unwrap();
        let event = resolved
            .tensile_initiator(resolution())
            .unwrap()
            .between(&[0.0; 6], &uniaxial(0, 2.0 * FT), None)
            .unwrap()
            .unwrap();
        let [before, after] = event.native().fraction_bracket();
        assert!(before <= expected_fraction && expected_fraction <= after);
        assert_eq!(
            event.source().material().state_coordinate("T"),
            Some(temperature)
        );
        let receiving = event.resolve_band(0.012).unwrap();
        assert_eq!(receiving.material(), resolved.material());
        let seed = event.seed(&receiving).unwrap();
        assert_eq!(seed.state.damage(), 0.0);
        let next = receiving
            .law()
            .trial(&uniaxial(0, 2.0 * FT), &seed.state)
            .unwrap();
        assert!(next.state.damage() > 0.0);
    }
}

#[test]
fn receiving_controls_geometry_and_cancelled_detection_refuse_without_losing_source() {
    let source = card(properties());
    let resolved = resolve(&source, 0.01);
    let detector = resolved.tensile_initiator(resolution()).unwrap();
    let start = uniaxial(0, 0.2 * FT);
    let end = uniaxial(0, 2.0 * FT);
    assert!(matches!(
        detector.between_controlled(&start, &end, None, &mut || true),
        Err(ResolvedTensileInitiationError::Native(
            TensileInitiationError::Cancelled
        ))
    ));
    assert!(
        detector
            .between(&[0.0; 6], &uniaxial(0, -FT), None)
            .unwrap()
            .is_none()
    );
    let event = detector.between(&start, &end, None).unwrap().unwrap();
    for width in [0.0, -0.01, f64::INFINITY, 5.0] {
        assert!(matches!(
            event.resolve_band(width),
            Err(ResolvedTensileInitiationError::Source(_))
        ));
    }
    let mut changed = settings(0.01);
    changed.energy_tolerance_j_m3 *= 2.0;
    let altered =
        ResolvedBrittleCrackBand::try_from_material(resolved.material().clone(), changed).unwrap();
    assert!(matches!(
        event.seed(&altered),
        Err(ResolvedTensileInitiationError::ControlsMismatch)
    ));
    let receiving = event.resolve_band(0.01).unwrap();
    let first = event.seed(&receiving).unwrap();
    let repeated = event.seed(&receiving).unwrap();
    assert_eq!(first, repeated);
    assert_eq!(receiving.material(), resolved.material());
}
