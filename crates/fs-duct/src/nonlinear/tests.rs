use super::*;
use crate::RiccatiSteps;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_math::det;
use fs_qty::{MolarMass, Pressure, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};
use std::sync::Arc;

fn gas(cp_linear: f64, specific_gas_constant: f64) -> FrozenGas {
    // Manufactured source cards exercise the production NASA-9 owner and its
    // bounded inversion. They are analytic fixtures, not measured gas data.
    let mut a = [0.0; 9];
    a[2] = 3.5;
    a[3] = cp_linear;
    a[7] = 1e12;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 0.001, 5_000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured nonlinear-duct caloric test".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("nonlinear-duct-fixture").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("test-gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / specific_gas_constant),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    FrozenGas::new(
        source,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(5_000.0)),
        || false,
    )
    .unwrap()
}
fn rest(_: f64) -> Result<InitialGasState, GasError> {
    Ok(InitialGasState {
        density_kg_m3: 1.2,
        velocity_m_s: 0.0,
        pressure_pa: 101_325.0,
    })
}
fn closed() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}
fn open() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
    }
}
fn cylinder(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn run(duct: &mut NonlinearDuct, end: f64) {
    let mut steps = 0;
    while duct.time().value() < end {
        let report = duct
            .step(
                Time::new(end - duct.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        assert!(report.maximum_cfl <= 0.40000000000001);
        assert!(report.conservation_residual.mass.abs() < 1e-12);
        assert!(report.conservation_residual.energy.abs() < 2e-9);
        steps += 1;
        assert!(steps < 20_000);
    }
}

#[test]
fn g0_segment_faces_and_exact_conical_volume_are_retained() {
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.2,
            },
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.03,
                length: 0.7,
            },
        ],
    };
    let duct = NonlinearDuct::new(&geometry, gas(0.0, 287.0), &[3, 11], closed(), rest).unwrap();
    assert_eq!(duct.geometry(), &geometry);
    assert_eq!(duct.solver().mesh().cell_count(), 14);
    assert_eq!(duct.solver().mesh().positions()[3], 0.2);
    let volume: f64 = duct.solver().mesh().volumes().iter().sum();
    let exact = core::f64::consts::PI
        * (0.01_f64.powi(2) * 0.2
            + 0.7 / 3.0 * (0.01_f64.powi(2) + 0.01 * 0.03 + 0.03_f64.powi(2)));
    relative(volume, exact, 2.0e-14);
    let sample = duct.sample(7).unwrap();
    let local_radius = 0.01 + 0.02 * (sample.axial_position_m - 0.2) / 0.7;
    relative(
        sample.area_m2,
        core::f64::consts::PI * local_radius * local_radius,
        2.0e-14,
    );
    assert_eq!(sample.volume_velocity_m3_s, 0.0);
}

#[test]
fn g1_profile_volume_matches_independent_integral_and_spline_cone_limit() {
    let profile = RadiusProfile::exponential(0.01, 3.0, 0.8).unwrap();
    let geometry = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(profile),
            steps: RiccatiSteps(2),
        }],
    };
    let duct = NonlinearDuct::new(&geometry, gas(0.0, 287.0), &[7], closed(), rest).unwrap();
    let actual: f64 = duct.solver().mesh().volumes().iter().sum();
    let exact = core::f64::consts::PI * 0.01_f64.powi(2) * (det::exp(2.4) - 1.0) / 3.0;
    relative(actual, exact, 4.0e-14);
    // Unequally spaced original stations deliberately fall inside solver cells.
    let profile = RadiusProfile::spline(
        vec![(2.0, 0.01), (2.17, 0.0134), (2.39, 0.0178), (3.0, 0.03)],
        Some([0.02, 0.02]),
    )
    .unwrap();
    let geometry = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(profile),
            steps: RiccatiSteps(2),
        }],
    };
    let duct = NonlinearDuct::new(&geometry, gas(0.0, 287.0), &[5], closed(), rest).unwrap();
    let actual: f64 = duct.solver().mesh().volumes().iter().sum();
    let exact = core::f64::consts::PI / 3.0 * (0.01_f64.powi(2) + 0.01 * 0.03 + 0.03_f64.powi(2));
    relative(actual, exact, 5.0e-14);
}

#[test]
fn g0_tapered_duct_remains_at_uniform_pressure_rest() {
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.2,
            },
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.035,
                length: 0.8,
            },
        ],
    };
    let mut duct =
        NonlinearDuct::new(&geometry, gas(0.001, 287.0), &[12, 32], closed(), rest).unwrap();
    let before = duct.solver().states().to_vec();
    run(&mut duct, 0.001);
    for (&a, b) in before.iter().zip(duct.solver().states()) {
        relative(b.density, a.density, 1.0e-12);
        relative(b.energy, a.energy, 1.0e-12);
        assert!(b.momentum.abs() < 1.0e-10);
    }
    assert_eq!(duct.time().value(), 0.001);
}

#[test]
fn g2_pressure_jump_propagates_through_the_actual_duct_adapter() {
    let mut duct = NonlinearDuct::new(&cylinder(0.02, 1.0), gas(0.0, 1.0), &[240], open(), |x| {
        Ok(if x < 0.5 {
            InitialGasState {
                density_kg_m3: 1.0,
                velocity_m_s: 0.0,
                pressure_pa: 1.0,
            }
        } else {
            InitialGasState {
                density_kg_m3: 0.125,
                velocity_m_s: 0.0,
                pressure_pa: 0.1,
            }
        })
    })
    .unwrap();
    run(&mut duct, 0.1);
    let samples = duct.samples().unwrap();
    let shocked = samples
        .iter()
        .find(|sample| sample.axial_position_m > 0.62)
        .unwrap();
    assert!(
        shocked.pressure_pa > 0.2 && shocked.velocity_m_s > 0.5,
        "{shocked:?}"
    );
    let untouched = samples
        .iter()
        .find(|sample| sample.axial_position_m > 0.85)
        .unwrap();
    relative(untouched.pressure_pa, 0.1, 1.0e-6);
    assert!(
        samples
            .iter()
            .all(|sample| sample.pressure_pa > 0.0 && sample.temperature_k > 0.0)
    );
}

#[test]
fn g2_finite_amplitude_wave_reflects_and_returns_with_reversed_flow() {
    let p0 = 101_325.0;
    let rho0 = p0 / (287.0 * 300.0);
    let c0 = det::sqrt(1.4 * 287.0 * 300.0);
    let mut duct = NonlinearDuct::new(
        &cylinder(0.02, 1.0),
        gas(0.0, 287.0),
        &[240],
        closed(),
        |x| {
            let perturbation = 1_000.0 * det::exp(-0.5 * ((x - 0.7) / 0.04).powi(2));
            let ratio = 1.0 + perturbation / p0;
            Ok(InitialGasState {
                density_kg_m3: rho0 * det::exp(det::ln(ratio) / 1.4),
                velocity_m_s: 2.0 * c0 / 0.4 * (det::exp(det::ln(ratio) * (0.4 / 2.8)) - 1.0),
                pressure_pa: p0 + perturbation,
            })
        },
    )
    .unwrap();
    assert!(
        duct.samples()
            .unwrap()
            .iter()
            .map(|s| s.velocity_m_s)
            .sum::<f64>()
            > 0.0
    );
    run(&mut duct, 0.0013);
    let samples = duct.samples().unwrap();
    let peak = samples
        .iter()
        .max_by(|a, b| a.pressure_pa.total_cmp(&b.pressure_pa))
        .unwrap();
    assert!(
        peak.axial_position_m > 0.78 && peak.axial_position_m < 0.93,
        "{peak:?}"
    );
    assert!(peak.pressure_pa > p0 + 400.0, "{peak:?}");
    assert!(
        peak.volume_velocity_m3_s < 0.0,
        "reflected wave did not reverse: {peak:?}"
    );
}

#[test]
fn g2_temperature_dependent_caloric_source_survives_shock_propagation() {
    let source = gas(0.002, 287.0);
    let identity = *source.source().components()[0].regions()[0].card_identity();
    let mut duct = NonlinearDuct::new(&cylinder(0.02, 1.0), source, &[120], open(), |x| {
        Ok(InitialGasState {
            density_kg_m3: 1.0,
            velocity_m_s: 0.0,
            pressure_pa: 287.0 * if x < 0.5 { 1_800.0 } else { 300.0 },
        })
    })
    .unwrap();
    run(&mut duct, 0.0002);
    assert_eq!(
        *duct.solver().gas().source().components()[0].regions()[0].card_identity(),
        identity
    );
    let samples = duct.samples().unwrap();
    assert!(samples.iter().any(|sample| sample.velocity_m_s > 100.0));
    for sample in &samples {
        let expected_gamma =
            (3.5 + 0.002 * sample.temperature_k) / (2.5 + 0.002 * sample.temperature_k);
        let actual_gamma = sample.sound_speed_m_s.powi(2) / (287.0 * sample.temperature_k);
        relative(actual_gamma, expected_gamma, 3.0e-12);
    }
    let cold_gamma = (3.5 + 0.002 * 300.0) / (2.5 + 0.002 * 300.0);
    assert!(samples.iter().any(|sample| {
        let gamma = sample.sound_speed_m_s.powi(2) / (287.0 * sample.temperature_k);
        (gamma - cold_gamma).abs() > 0.05
    }));
}

#[test]
fn g0_unsupported_junctions_and_invalid_meshes_refuse() {
    assert!(
        NonlinearDuct::new(&cylinder(0.01, 1.0), gas(0.0, 287.0), &[], closed(), rest).is_err()
    );
    assert!(
        NonlinearDuct::new(&cylinder(0.01, 1.0), gas(0.0, 287.0), &[0], closed(), rest).is_err()
    );
    assert!(
        NonlinearDuct::new(
            &cylinder(0.01, 1.0),
            gas(0.0, 287.0),
            &[MAX_NONLINEAR_DUCT_CELLS + 1],
            closed(),
            rest
        )
        .is_err()
    );
    let jump = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 1.0,
            },
            Segment::Cylinder {
                radius: 0.02,
                length: 1.0,
            },
        ],
    };
    assert!(matches!(
        NonlinearDuct::new(&jump, gas(0.0, 287.0), &[1, 1], closed(), rest),
        Err(NonlinearDuctError::UnsupportedGeometry(_))
    ));
    let hole = Duct {
        segments: vec![Segment::ToneHole {
            hole_radius: 0.001,
            chimney_height: 0.001,
            bore_radius: 0.01,
            state: crate::HoleState::Open,
        }],
    };
    assert!(matches!(
        NonlinearDuct::new(&hole, gas(0.0, 287.0), &[1], closed(), rest),
        Err(NonlinearDuctError::UnsupportedGeometry(_))
    ));
    let duct =
        NonlinearDuct::new(&cylinder(0.01, 1.0), gas(0.0, 287.0), &[3], closed(), rest).unwrap();
    assert!(duct.sample(3).is_err());
}

#[test]
fn g4_cancellation_publishes_no_partial_duct_or_state() {
    let geometry = cylinder(0.01, 1.0);
    assert!(matches!(
        NonlinearDuct::new_controlled(&geometry, gas(0.0, 287.0), &[8], closed(), rest, || true),
        Err(NonlinearDuctError::Gas(GasError::Cancelled))
    ));
    let mut duct = NonlinearDuct::new(&geometry, gas(0.0, 287.0), &[8], closed(), rest).unwrap();
    let before = duct.solver().states().to_vec();
    let mut polls = 0;
    let result = duct.step(Time::new(0.001), StepControl::default(), || {
        polls += 1;
        polls > 2
    });
    assert!(matches!(
        result,
        Err(NonlinearDuctError::Gas(GasError::Cancelled))
    ));
    assert_eq!(duct.solver().states(), before.as_slice());
    assert_eq!(duct.time().value(), 0.0);
}
