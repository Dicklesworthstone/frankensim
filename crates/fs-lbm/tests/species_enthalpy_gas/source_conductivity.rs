//! Original transport-card integrals in the existing nonlinear gas heat owner.
use super::*;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_lbm::conjugate::{
    TensorNonlinearConfig, ThermalFace, ThermalSetup, ThermochemicalConductionControls,
    ThermochemicalConductionScheme, ThermochemicalConductivityModel,
    ThermochemicalSpeciesConductionStep,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{
    SpeciesTransportMaterial, conductivity::GasConductivityQuadrature,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, QuantitySpec, Temperature};

fn transport_model(datum: f64, conductivity: &[(f64, f64)]) -> ThermochemicalGasModel {
    let original = source::constant_cp_model(datum, 1.0);
    let [lo, hi] = [conductivity[0].0, conductivity.last().unwrap().0];
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let artifact = ContentHash([0x6B; 32]);
    let provenance = Provenance {
        source: "authored synthetic conductivity coupling fixture; not measured".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic pure gas".into(),
            method: "authored algebra".into(),
            artifact,
            caveats: "not empirical".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, knots) in [
        (
            "dynamic_viscosity",
            DynViscosity::DIMS,
            vec![(lo, 1.8e-5), (hi, 1.8e-5)],
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            conductivity.to_vec(),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    temperature,
                    lo,
                    hi,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "gas".into(),
            process: "authored temperature law".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ThermochemicalGasModel::new(
        original.thermodynamics().clone(),
        vec![SpeciesTransportMaterial {
            species: SpeciesId::new("manufactured").unwrap(),
            material,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        None,
    )
    .unwrap()
}

fn heat_controls(scheme: ThermochemicalConductionScheme) -> ThermochemicalConductionControls {
    let mut thermal = TensorNonlinearConfig::default();
    thermal.linear.tolerance = 1e-10;
    ThermochemicalConductionControls {
        scheme,
        thermal,
        maximum_caloric_queries: 100_000,
        maximum_enthalpy_defect_j: 1e-3,
        conductivity: ThermochemicalConductivityModel::SourceIntegral {
            quadrature: GasConductivityQuadrature {
                absolute_tolerance_w_m_k: 1e-12,
                relative_tolerance: 1e-11,
                maximum_subdivisions: 128,
                maximum_depth: 8,
            },
            maximum_evaluations: 100_000,
        },
    }
}
fn pure_h(model: &ThermochemicalGasModel, t: f64) -> f64 {
    model
        .evaluate(&source::point(t, 100_000.0))
        .unwrap()
        .specific_enthalpy_j_kg()
}

#[allow(clippy::too_many_arguments)]
fn advance<'a>(
    accepted: &MaterialAdvectionStep,
    domain: &VoxelDomain,
    inventory: &[[f64; 1]],
    h: &[f64],
    dt: f64,
    gas: &GasSpeciesEnthalpySource<'a, 1>,
    controls: ThermochemicalConductionControls,
    gate: &CancelGate,
) -> Result<ThermochemicalSpeciesConductionStep<'a, 1>, ChtError> {
    accepted.conduct_thermochemical_species(
        domain,
        inventory,
        h,
        &ThermalSetup::new([ThermalFace::Adiabatic; 6]),
        0.0,
        dt,
        gas,
        controls,
        gate,
    )
}

#[test]
fn original_kirchhoff_integral_crosses_source_knots_and_preserves_energy_references() {
    let domain = VoxelDomain::graded([vec![0.5, 1.5], vec![1.0], vec![1.0]]).unwrap();
    let masses = [1.0, 2.0];
    let accepted = carrier(&domain, &masses);
    let before = accepted.clone();
    let inventory = masses.map(|m| [m]);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let gate = CancelGate::new_clock_free();
    let d = [0.0; 2];
    let initial_t = [480.0, 320.0];
    let cp = 3.0 * R / MOLAR_MASS[0];
    let capacity = masses.map(|m| m * cp);
    let dt = 10_000.0;
    for extra_slope in [0.0, 2e-4] {
        let knots = [
            (200.0, 0.02),
            (400.0, 0.04),
            (700.0, 0.07 + 300.0 * extra_slope),
        ];
        let primitive = |t: f64| 0.5e-4 * t * t + 0.5 * extra_slope * (t - 400.0).max(0.0).powi(2);
        // Independent scalar BE equation. Both materials share this source,
        // so q=dt*A/(dL+dR)*(Phi(TL)-Phi(TR)), with dL+dR=1 m.
        let equilibrium =
            (capacity[0] * initial_t[0] + capacity[1] * initial_t[1]) / (capacity[0] + capacity[1]);
        let (mut low, mut high) = (0.0, capacity[0] * (initial_t[0] - equilibrium));
        for _ in 0..100 {
            let q = 0.5 * (low + high);
            let t = [
                initial_t[0] - q / capacity[0],
                initial_t[1] + q / capacity[1],
            ];
            if q > dt * (primitive(t[0]) - primitive(t[1])) {
                high = q;
            } else {
                low = q;
            }
        }
        let q = 0.5 * (low + high);
        let expected = [
            initial_t[0] - q / capacity[0],
            initial_t[1] + q / capacity[1],
        ];
        assert!(expected[0] > 400.0 && expected[1] < 400.0);
        let mut reference: Option<[f64; 2]> = None;
        for datum in [-4000.0, 0.0, 10000.0] {
            let model = transport_model(datum, &knots);
            let gas = GasSpeciesEnthalpySource::new(&model, &names, &context, &d, energy(), &gate)
                .unwrap();
            let h = initial_t.map(|t| pure_h(&model, t));
            let controls = heat_controls(ThermochemicalConductionScheme::BackwardEuler);
            let result = advance(
                &accepted, &domain, &inventory, &h, dt, &gas, controls, &gate,
            )
            .unwrap();
            for cell in 0..2 {
                near(result.state.temperatures_k()[cell], expected[cell], 2e-6);
                near(
                    masses[cell] * (result.specific_enthalpy_j_kg()[cell] - h[cell]),
                    if cell == 0 { -q } else { q },
                    1e-3,
                );
            }
            let actual = [
                result.state.temperatures_k()[0],
                result.state.temperatures_k()[1],
            ];
            if let Some(old) = reference {
                for cell in 0..2 {
                    near(actual[cell], old[cell], 2e-6);
                }
            }
            reference = Some(actual);
            near(result.thermal.record().stored_energy_change_j, 0.0, 1e-3);
            assert!(result.report.absolute_enthalpy_defect_j <= 1e-3);
            assert!(result.report.conductivity_evaluations > 0);
            assert!(result.report.maximum_estimated_conductivity_error_w_m_k <= 4e-12);
            let held = advance(
                &accepted,
                &domain,
                &inventory,
                &h,
                dt,
                &gas,
                ThermochemicalConductionControls {
                    conductivity: ThermochemicalConductivityModel::HeldAtOpening,
                    ..controls
                },
                &gate,
            )
            .unwrap();
            let held_heat = masses[0] * (h[0] - held.specific_enthalpy_j_kg()[0]);
            assert!((held_heat - q).abs() > 10.0);
        }
    }
    assert_eq!(accepted, before);
    assert_eq!(inventory, [[1.0], [2.0]]);
}

#[test]
fn constant_original_conductivity_matches_native_be_and_sdirk_held_laws() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let masses = [1.0, 2.0];
    let accepted = carrier(&domain, &masses);
    let inventory = masses.map(|m| [m]);
    let model = transport_model(10000.0, &[(200.0, 0.04), (700.0, 0.04)]);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&model, &names, &context, &d, energy(), &gate).unwrap();
    let h = [480.0, 320.0].map(|t| pure_h(&model, t));
    let capacity = masses.map(|m| m * 3.0 * R / MOLAR_MASS[0]);
    let dt = 10_000.0;
    let z = dt * 0.04 * (1.0 / capacity[0] + 1.0 / capacity[1]);
    let mean = (capacity[0] * 480.0 + capacity[1] * 320.0) / (capacity[0] + capacity[1]);
    for scheme in [
        ThermochemicalConductionScheme::BackwardEuler,
        ThermochemicalConductionScheme::Sdirk2,
    ] {
        let controls = heat_controls(scheme);
        let source = advance(
            &accepted, &domain, &inventory, &h, dt, &gas, controls, &gate,
        )
        .unwrap();
        let held = advance(
            &accepted,
            &domain,
            &inventory,
            &h,
            dt,
            &gas,
            ThermochemicalConductionControls {
                conductivity: ThermochemicalConductivityModel::HeldAtOpening,
                ..controls
            },
            &gate,
        )
        .unwrap();
        let factor = match scheme {
            ThermochemicalConductionScheme::BackwardEuler => 1.0 / (1.0 + z),
            ThermochemicalConductionScheme::Sdirk2 => {
                let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
                (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
            }
        };
        let expected = [
            mean + 2.0 / 3.0 * 160.0 * factor,
            mean - 1.0 / 3.0 * 160.0 * factor,
        ];
        for cell in 0..2 {
            near(source.state.temperatures_k()[cell], expected[cell], 2e-6);
            near(
                source.state.temperatures_k()[cell],
                held.state.temperatures_k()[cell],
                2e-6,
            );
        }
        near(source.thermal.record().closure_j, 0.0, 1e-3);
        assert!(source.report.conductivity_evaluations > 0);
    }
}

#[test]
fn conductivity_work_source_bracket_and_cancellation_refusals_are_atomic() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let accepted = carrier(&domain, &[1.0]);
    let before = accepted.clone();
    let inventory = [[1.0]];
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0];
    let gate = CancelGate::new_clock_free();
    let model = transport_model(0.0, &[(200.0, 0.02), (700.0, 0.07)]);
    let gas = GasSpeciesEnthalpySource::new(&model, &names, &context, &d, energy(), &gate).unwrap();
    let h = [pure_h(&model, 400.0)];
    let mut limited = heat_controls(ThermochemicalConductionScheme::BackwardEuler);
    if let ThermochemicalConductivityModel::SourceIntegral {
        maximum_evaluations,
        ..
    } = &mut limited.conductivity
    {
        *maximum_evaluations = 1;
    }
    assert!(
        advance(
            &accepted, &domain, &inventory, &h, 1.0, &gas, limited, &gate,
        )
        .is_err()
    );
    // A valid opening transport sample does not authorize conductivity outside
    // the original card's narrower support across the requested caloric branch.
    let narrow_model = transport_model(0.0, &[(300.0, 0.03), (450.0, 0.045)]);
    let narrow =
        GasSpeciesEnthalpySource::new(&narrow_model, &names, &context, &d, energy(), &gate)
            .unwrap();
    let controls = heat_controls(ThermochemicalConductionScheme::BackwardEuler);
    assert!(
        advance(
            &accepted, &domain, &inventory, &h, 1.0, &narrow, controls, &gate,
        )
        .is_err()
    );
    let held = advance(
        &accepted,
        &domain,
        &inventory,
        &h,
        1.0,
        &narrow,
        ThermochemicalConductionControls {
            conductivity: ThermochemicalConductivityModel::HeldAtOpening,
            ..controls
        },
        &gate,
    )
    .unwrap();
    assert_eq!(held.specific_enthalpy_j_kg(), h);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        advance(
            &accepted, &domain, &inventory, &h, 1.0, &gas, controls, &cancelled,
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(accepted, before);
    assert_eq!(inventory, [[1.0]]);
    assert_eq!(h, [pure_h(&model, 400.0)]);
}
