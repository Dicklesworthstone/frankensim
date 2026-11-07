//! Source-state consistency, physical strip integration and original coating heat/mass.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::fluid_transport::ThermalFluidTransportState;
use fs_material::state_point::MaterialPropertySelection;
use fs_thermal::sublimation::{
    MOLAR_GAS_CONSTANT, SolverConfig, SublimationError, SublimationLaw, SublimationProperties,
    VaporReservoir,
    coating::{SupportedSurfaceState, advance_supported_surface_with_heat},
};

fn claim(axis: &str, knots: Vec<(f64, f64)>) -> PropertyClaim {
    let temperature_range = if axis == "temperature" {
        [knots[0].0, knots[knots.len() - 1].0]
    } else {
        [250.0, 350.0]
    };
    PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: if axis == "pressure" {
                Pressure::DIMS
            } else {
                Temperature::DIMS
            },
            knots,
            dims: BINARY_DIFFUSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                temperature_range[0],
                temperature_range[1],
            )
            .with_quantity("pressure", typed(QuantityKind::Pressure), 50000.0, 150000.0)
            .with("carrier_fraction", 0.99, 0.99),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured binary gas film fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "declared A/B pair".into(),
            phase: "gas".into(),
            process: "synthetic analytic fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn constant_card(d: f64) -> MaterialCard {
    card(vec![claim("temperature", vec![(250.0, d), (350.0, d)])])
}
fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            280.0,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 60000.0)
        .unwrap()
        .with("carrier_fraction", 0.99)
        .unwrap()
}
fn pair<'a>(material: &'a MaterialCard, context: &'a QueryPoint) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial {
        species: [SpeciesId::new("B").unwrap(), SpeciesId::new("A").unwrap()],
        material,
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn fluid(temperature: f64) -> ThermalFluidTransportState {
    ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(1.0, 1000.0, 0.1, 4.0e-5).unwrap(),
        temperature,
        None,
    )
    .unwrap()
}
fn conditions() -> LaminarPlateFilmConditions {
    LaminarPlateFilmConditions {
        start: Length::new(0.0),
        end: Length::new(0.1),
        speed_m_s: 1.0,
        temperature: Temperature::new(300.0),
        total_pressure: Pressure::new(100000.0),
    }
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= 1e-20 + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn film_queries_actual_total_pressure_and_retains_pair_source_and_fluid() {
    let material = card(vec![claim(
        "pressure",
        vec![(50000.0, 1e-5), (150000.0, 3e-5)],
    )]);
    let context = point();
    let binding = pair(&material, &context);
    let gas = fluid(300.0);
    let vapor = SpeciesId::new("A").unwrap();
    let film =
        evaluate_laminar_plate_gas_film(&gas, &binding, &vapor, conditions(), || false).unwrap();
    assert!(core::ptr::eq(film.fluid(), &gas));
    assert!(core::ptr::eq(film.diffusion_source().material, &material));
    assert_eq!(film.vapor_species(), &vapor);
    assert_eq!(film.carrier_species().as_str(), "B");
    assert_eq!(film.diffusion().species[0], vapor);
    near(film.diffusion().diffusivity_m2_s, 2e-5, 1e-15);
    let resolved = &film.diffusion().material;
    assert_eq!(resolved.state_coordinate("temperature"), Some(300.0));
    assert_eq!(resolved.state_coordinate("pressure"), Some(100000.0));
    assert_eq!(resolved.state_coordinate("carrier_fraction"), Some(0.99));
    let receipt = &resolved
        .property(BINARY_DIFFUSIVITY_PROPERTY)
        .unwrap()
        .answer()
        .receipt;
    assert!(!receipt.source_hashes.is_empty());
    assert_eq!(
        receipt.selected,
        material.claims_for(BINARY_DIFFUSIVITY_PROPERTY)[0].0
    );
    near(film.reynolds(), 2500.0, 1e-15);
    near(film.schmidt(), 2.0, 1e-15);
    let expected_sh = 0.664 * 50.0 * 2.0_f64.powf(1.0 / 3.0);
    near(film.downstream_mean_sherwood(), expected_sh, 2e-13);
    near(
        film.mass_transfer_coefficient_m_s(),
        expected_sh * 2e-5 / 0.1,
        2e-13,
    );
    assert_eq!(binding.context, &context);
}

#[test]
fn molecular_diffusivity_is_independent_of_heat_diffusivity_and_scales_physically() {
    let a = constant_card(2e-5);
    let b = constant_card(4e-5);
    let context = point();
    let gas = fluid(300.0);
    let vapor = SpeciesId::new("A").unwrap();
    let first =
        evaluate_laminar_plate_gas_film(&gas, &pair(&a, &context), &vapor, conditions(), || false)
            .unwrap();
    let faster =
        evaluate_laminar_plate_gas_film(&gas, &pair(&b, &context), &vapor, conditions(), || false)
            .unwrap();
    near(
        faster.mass_transfer_coefficient_m_s() / first.mass_transfer_coefficient_m_s(),
        2.0_f64.powf(2.0 / 3.0),
        2e-13,
    );
    let other_heat = ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(1.0, 4000.0, 0.01, 4e-5).unwrap(),
        300.0,
        None,
    )
    .unwrap();
    let same = evaluate_laminar_plate_gas_film(
        &other_heat,
        &pair(&a, &context),
        &vapor,
        conditions(),
        || false,
    )
    .unwrap();
    assert_eq!(
        same.mass_transfer_coefficient_m_s(),
        first.mass_transfer_coefficient_m_s()
    );
}

#[test]
fn strip_refinement_preserves_original_leading_edge_transfer_and_thin_strip_limit() {
    let material = constant_card(2e-5);
    let context = point();
    let gas = fluid(300.0);
    let vapor = SpeciesId::new("A").unwrap();
    let binding = pair(&material, &context);
    let whole =
        evaluate_laminar_plate_gas_film(&gas, &binding, &vapor, conditions(), || false).unwrap();
    let mut integrated = 0.0;
    for [a, b] in [[0.0, 0.001], [0.001, 0.037], [0.037, 0.1]] {
        let section = whole
            .at_film_state(
                &gas,
                LaminarPlateFilmConditions {
                    start: Length::new(a),
                    end: Length::new(b),
                    ..conditions()
                },
                || false,
            )
            .unwrap();
        integrated += section.mass_transfer_coefficient_m_s() * (b - a);
    }
    near(
        integrated,
        whole.mass_transfer_coefficient_m_s() * 0.1,
        3e-13,
    );
    let thin = whole
        .at_film_state(
            &gas,
            LaminarPlateFilmConditions {
                start: Length::new(0.1_f64.next_down()),
                ..conditions()
            },
            || false,
        )
        .unwrap();
    near(
        thin.mass_transfer_coefficient_m_s(),
        0.5 * whole.mass_transfer_coefficient_m_s(),
        2e-15,
    );
}

#[test]
fn refreshed_film_keeps_the_original_selected_law_and_refuses_a_source_switch() {
    let material = card(vec![
        claim("temperature", vec![(250.0, 1e-5), (320.0, 2e-5)]),
        claim("temperature", vec![(320.0, 3e-5), (350.0, 4e-5)]),
    ]);
    let context = point();
    let binding = pair(&material, &context);
    let gas = fluid(300.0);
    let vapor = SpeciesId::new("A").unwrap();
    let film =
        evaluate_laminar_plate_gas_film(&gas, &binding, &vapor, conditions(), || false).unwrap();
    let warm = fluid(310.0);
    let next = film
        .at_film_state(
            &warm,
            LaminarPlateFilmConditions {
                temperature: Temperature::new(310.0),
                ..conditions()
            },
            || false,
        )
        .unwrap();
    near(
        next.diffusion().diffusivity_m2_s,
        1e-5 + 60.0 / 70.0 * 1e-5,
        2e-15,
    );
    let selected = |film: &SourceLaminarPlateGasFilm<'_, '_, ThermalFluidTransportState>| {
        film.diffusion()
            .material
            .property(BINARY_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected
    };
    assert_eq!(selected(&next), selected(&film));
    let hot = fluid(330.0);
    let hot_conditions = LaminarPlateFilmConditions {
        temperature: Temperature::new(330.0),
        ..conditions()
    };
    assert!(matches!(
        film.at_film_state(&hot, hot_conditions, || false),
        Err(GasFilmError::Source(_))
    ));
    let newly_selected =
        evaluate_laminar_plate_gas_film(&hot, &binding, &vapor, hot_conditions, || false).unwrap();
    assert_ne!(selected(&newly_selected), selected(&film));
}

#[test]
fn invalid_pair_state_regime_alias_and_cancellation_refuse_without_extrapolation() {
    let material = constant_card(2e-5);
    let context = point();
    let gas = fluid(300.0);
    let binding = pair(&material, &context);
    let vapor = SpeciesId::new("A").unwrap();
    assert!(matches!(
        evaluate_laminar_plate_gas_film(
            &gas,
            &binding,
            &SpeciesId::new("foreign").unwrap(),
            conditions(),
            || false
        ),
        Err(GasFilmError::Invalid(_))
    ));
    let mut duplicate = binding.clone();
    duplicate.species[1] = duplicate.species[0].clone();
    assert!(matches!(
        evaluate_laminar_plate_gas_film(&gas, &duplicate, &vapor, conditions(), || false),
        Err(GasFilmError::Invalid(_))
    ));
    for changed in [
        LaminarPlateFilmConditions {
            speed_m_s: 0.0,
            ..conditions()
        },
        LaminarPlateFilmConditions {
            start: Length::new(0.1),
            ..conditions()
        },
        LaminarPlateFilmConditions {
            temperature: Temperature::new(301.0),
            ..conditions()
        },
        LaminarPlateFilmConditions {
            total_pressure: Pressure::new(f64::NAN),
            ..conditions()
        },
    ] {
        assert!(matches!(
            evaluate_laminar_plate_gas_film(&gas, &binding, &vapor, changed, || false),
            Err(GasFilmError::Invalid(_))
        ));
    }
    assert!(matches!(
        evaluate_laminar_plate_gas_film(
            &gas,
            &binding,
            &vapor,
            LaminarPlateFilmConditions {
                speed_m_s: 201.0,
                ..conditions()
            },
            || false
        ),
        Err(GasFilmError::OutsideCorrelation { .. })
    ));
    let too_diffusive = constant_card(1e-4);
    assert!(matches!(
        evaluate_laminar_plate_gas_film(
            &gas,
            &pair(&too_diffusive, &context),
            &vapor,
            conditions(),
            || false
        ),
        Err(GasFilmError::OutsideCorrelation { .. })
    ));
    let alias = context
        .clone()
        .with_quantity(
            "film_temperature",
            typed(QuantityKind::AbsoluteTemperature),
            280.0,
        )
        .unwrap();
    assert!(matches!(
        evaluate_laminar_plate_gas_film(
            &gas,
            &pair(&material, &alias),
            &vapor,
            conditions(),
            || false
        ),
        Err(GasFilmError::Invalid(_))
    ));
    for (axis, value) in [("T", 280.0), ("p", 60000.0)] {
        let legacy = context.clone().with(axis, value).unwrap();
        assert!(matches!(
            evaluate_laminar_plate_gas_film(
                &gas,
                &pair(&material, &legacy),
                &vapor,
                conditions(),
                || false
            ),
            Err(GasFilmError::Invalid(_))
        ));
    }
    let mut calls = 0;
    assert!(matches!(
        evaluate_laminar_plate_gas_film(&gas, &binding, &vapor, conditions(), || {
            calls += 1;
            calls == 3
        }),
        Err(GasFilmError::Cancelled)
    ));
    assert_eq!(binding.context, &context);
}

struct Deposit;
impl SublimationLaw for Deposit {
    fn temperature_range_k(&self) -> [f64; 2] {
        [250.0, 350.0]
    }
    fn properties(&self, _: f64) -> Result<SublimationProperties, SublimationError> {
        Ok(SublimationProperties {
            saturation_pressure_pa: 10.0,
            saturation_pressure_derivative_pa_k: 0.0,
            latent_heat_j_kg: 2.8e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}

#[test]
fn source_film_drives_bare_deposition_with_original_latent_heat_and_mass_conservation() {
    let material = constant_card(2e-5);
    let context = point();
    let gas = fluid(300.0);
    let film = evaluate_laminar_plate_gas_film(
        &gas,
        &pair(&material, &context),
        &SpeciesId::new("A").unwrap(),
        conditions(),
        || false,
    )
    .unwrap();
    let state = SupportedSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 500.0 * 300.0,
        deposit_mass_kg: 0.0,
        area_m2: 0.01,
    };
    let kinetic = (0.018 / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * 300.0)).sqrt();
    let film_pressure = film.mass_transfer_coefficient_m_s() / (MOLAR_GAS_CONSTANT / 0.018 * 300.0);
    // Independent two-resistance solution at 300 K, with supplied cooling
    // exactly removing deposition heat. No numerical deposit seed is present.
    let expected_mass = -10.0 / (1.0 / kinetic + 1.0 / film_pressure) * state.area_m2 * 10.0;
    let supplied_heat = expected_mass * 2.8e6;
    let step = advance_supported_surface_with_heat(
        &Deposit,
        film.surface_kinetics(1.0).unwrap(),
        state,
        VaporReservoir::Prescribed {
            partial_pressure_pa: 20.0,
            temperature_k: 300.0,
        },
        10.0,
        supplied_heat,
        SolverConfig::default(),
        |h| Ok(h / 500.0),
        |t| Ok(2000.0 * t),
        || Ok(()),
    )
    .unwrap();
    near(step.phase.transferred_mass_kg, expected_mass, 2e-12);
    near(step.state.deposit_mass_kg, -expected_mass, 2e-12);
    near(step.phase.surface_temperature_k, 300.0, 2e-13);
    assert_eq!(step.state.substrate_mass_kg, state.substrate_mass_kg);
    near(
        step.substrate_energy_change_j + step.deposit_energy_change_j + step.phase.vapor_energy_j,
        supplied_heat,
        2e-10,
    );
    assert!(step.phase.recession_m < 0.0);
    assert_eq!(
        film.surface_kinetics(0.0)
            .unwrap()
            .accommodation_coefficient,
        0.0
    );
    assert!(film.surface_kinetics(1.01).is_err());
}
