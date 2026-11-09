//! Accepted-temperature optical feedback and transactional source support.

#[path = "thermo_optical/analytic.rs"]
mod analytic;

#[path = "thermo_optical/schott.rs"]
mod schott;
#[path = "thermo_optical/spectral.rs"]
mod spectral;

use super::*;
use fs_couple::layered_wave_heating::{ThermoOpticalLayer, ThermoOpticalWaveHeating};
use fs_em::plane_wave::{OpticalConductorResponse, OpticalMaterialResponse};
use fs_material::state_point::{
    VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, VISIBLE_COMPLEX_IOR_K_PROPERTIES,
    VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
    VISIBLE_DIELECTRIC_IOR_PROPERTIES, VisibleConductorStatePoint, VisibleOpticalStatePoint,
    resolve_visible_optical_state_point,
};

const LIGHT_SPEED: f64 = 299_792_458.0;
const OPTICAL_OMEGA: f64 = core::f64::consts::TAU * LIGHT_SPEED / 550e-9;
const WIDTH: f64 = 70e-9;

fn optical_card(name: &str, upper: f64) -> MaterialCard {
    let base = card(name, 1.0, 0.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    for (names, endpoints) in [
        (VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, [0.3, 1.5]),
        (VISIBLE_COMPLEX_IOR_K_PROPERTIES, [3.2, 1.6]),
    ] {
        for name in names {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, Dims::NONE),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(280.0, endpoints[0]), (360.0, endpoints[1])],
                        dims: Dims::NONE,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "T",
                            quantity(QuantityKind::AbsoluteTemperature),
                            280.0,
                            upper,
                        )
                        .with_quantity(
                            "pressure",
                            quantity(QuantityKind::Pressure),
                            90_000.0,
                            200_000.0,
                        ),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "authored thermo-optical source/heat oracle".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn optical_source(card: &MaterialCard, t: f64, pressure: f64) -> VisibleConductorStatePoint {
    let mut requirements: Vec<_> = VISIBLE_COMPLEX_IOR_ETA_PROPERTIES
        .into_iter()
        .chain(VISIBLE_COMPLEX_IOR_K_PROPERTIES)
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                Dims::NONE,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        })
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_new(MU_NAME, MU, ScalarAdmissibility::StrictlyPositive)
            .unwrap(),
    );
    VisibleConductorStatePoint::try_from_material(
        resolve_material_state_point(
            card,
            &at(t, pressure),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
    .unwrap()
}

fn glass_card(upper: f64) -> MaterialCard {
    let base = card("authored thermo-optic glass", 1.0, 0.0);
    let mut claims = base.claims().clone();
    for (names, quantity, values) in [
        (
            VISIBLE_DIELECTRIC_IOR_PROPERTIES,
            QuantitySpec::dimensional(Dims::NONE),
            [1.4, 1.56],
        ),
        (
            VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
            VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
            [5e5, 2.1e6],
        ),
    ] {
        for name in names {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, quantity),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(280.0, values[0]), (360.0, values[1])],
                        dims: quantity.dims(),
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "T",
                            super::quantity(QuantityKind::AbsoluteTemperature),
                            280.0,
                            upper,
                        )
                        .with_quantity(
                            "pressure",
                            super::quantity(QuantityKind::Pressure),
                            90_000.0,
                            200_000.0,
                        ),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "authored spectral glass heat oracle".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn glass_source(card: &MaterialCard) -> VisibleOpticalStatePoint {
    resolve_visible_optical_state_point(
        card,
        &at(300.0, 100_000.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_additional_requirements(
        card,
        &[
            ScalarPropertyRequirement::try_new(MU_NAME, MU, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn exterior() -> ResolvedMaterialStatePoint {
    state(
        &card(
            "optical exterior",
            1.0 / (M * LIGHT_SPEED * LIGHT_SPEED * E),
            0.0,
        ),
        false,
        100_000.0,
    )
}

fn initial_at<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperatures: &[f64],
) -> DepositionSimulation<'a> {
    let n = temperatures.len();
    let events = bound
        .birth_events(0.0, temperatures, conduction(), || Ok(()))
        .unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: n,
            max_links: n - 1,
            max_events: 2 * n - 1,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn native_nonuniform_temperatures_refresh_each_original_optical_cell_and_feed_next_step() {
    let card = optical_card("temperature-dependent film", 360.0);
    let source = optical_source(&card, 300.0, 100_000.0);
    let outside = exterior();
    let native = body(&[&card], &[WIDTH], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial_temperatures = [290.0, 330.0];
    let initial = initial_at(&bound, &initial_temperatures);
    let initial_temperatures = bound
        .cell_ids()
        .into_iter()
        .map(|id| initial.temperature_k(id).unwrap())
        .collect::<Vec<_>>();
    let before = initial.cells().clone();
    let owner = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        0.8,
        Polarization::Tm,
        &[ThermoOpticalLayer {
            source: &source,
            card: &card,
        }],
        || false,
    )
    .unwrap();
    let input = [C64::new(2e6, 0.0), C64::new(3e5, 2e5)];
    let drive = LayeredWaveDrive::Coherent(input);
    let dt = 1e-11;
    let refreshed: Vec<_> = initial_temperatures
        .iter()
        .map(|&t| source.at_query_point(&card, &at(t, 100_000.0)).unwrap())
        .collect();
    let responses: Vec<_> = refreshed
        .iter()
        .map(|s| {
            DielectricResponse::OpticalConductor(OpticalConductorResponse::try_new(s).unwrap())
        })
        .collect();
    let reference = LayeredWaveHeating::try_new_cellwise(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        &responses,
        || false,
    )
    .unwrap();
    let mut powers = [0.0; 2];
    reference
        .prepare(0.8, Polarization::Tm, || false)
        .unwrap()
        .cell_power_into(drive, &mut powers, || false)
        .unwrap();
    let frozen = LayeredWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        &[DielectricResponse::OpticalConductor(
            OpticalConductorResponse::try_new(&source).unwrap(),
        )],
        || false,
    )
    .unwrap();
    let mut frozen_powers = [0.0; 2];
    frozen
        .prepare(0.8, Polarization::Tm, || false)
        .unwrap()
        .cell_power_into(drive, &mut frozen_powers, || false)
        .unwrap();
    assert!(
        powers
            .iter()
            .zip(frozen_powers)
            .any(|(&a, b)| (a - b).abs() > 0.01 * a)
    );

    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let controls = LayeredThermalSolveControls {
        cell_tolerance_j: 1e-12,
        total_tolerance_j: 1e-12,
        ..LayeredThermalSolveControls::default()
    };
    let mut total_polls = 0;
    let result = owner
        .advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();

    // Independent constant-capacity two-cell backward-Euler heat equation.
    // The nonuniform initial state conducts as well as absorbing local light.
    let c = RHO * AREA * (WIDTH / 2.0) * CP;
    let g_dt = AREA * K / (WIDTH / 2.0) * dt;
    let mean = 0.5 * (initial_temperatures[0] + initial_temperatures[1])
        + 0.5 * dt * (powers[0] + powers[1]) / c;
    let contrast = (0.5 * c * (initial_temperatures[0] - initial_temperatures[1])
        + 0.5 * dt * (powers[0] - powers[1]))
        / (c + 2.0 * g_dt);
    let expected = [mean + contrast, mean - contrast];
    assert_eq!(result.applied_optical_states, refreshed);
    for (i, id) in bound.cell_ids().into_iter().enumerate() {
        let temperature = result.step.thermal.temperature_k(id).unwrap();
        near(temperature, expected[i], 3e-9);
        near(result.step.cell_energy_j[i].1, powers[i] * dt, 1e-13);
        assert_eq!(result.step.cell_energy_j[i].0, id);
        assert_eq!(
            result.endpoint_optical_states[i],
            source
                .at_query_point(&card, &at(temperature, 100_000.0))
                .unwrap()
        );
    }
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.optical.absorbed_j,
        3e-12,
    );
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);

    let next = owner
        .advance_to(
            &result.step.thermal,
            2.0 * dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
    assert_ne!(next.applied_optical_states, result.applied_optical_states);

    // Interrupt the final checkpoint after both source refreshes and the heat
    // solve: the caller must still receive no partially committed endpoint.
    let mut polls = 0;
    assert!(matches!(
        owner.advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || {
                polls += 1;
                polls == total_polls
            },
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(initial.cells(), &before);
}

#[test]
fn optical_source_substitution_and_endpoint_domain_exit_refuse_without_thermal_commit() {
    let card = optical_card("narrow optical support", 301.0);
    let source = optical_source(&card, 300.0, 100_000.0);
    let other_card = optical_card("different provenance same values", 301.0);
    let wrong_pressure = optical_source(&card, 300.0, 150_000.0);
    let outside = exterior();
    let native = body(&[&card], &[WIDTH], &[1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    for binding in [
        ThermoOpticalLayer {
            source: &source,
            card: &other_card,
        },
        ThermoOpticalLayer {
            source: &wrong_pressure,
            card: &card,
        },
    ] {
        assert!(
            ThermoOpticalWaveHeating::try_new(
                &bound,
                &outside,
                DielectricResponse::Lossless(&outside),
                OPTICAL_OMEGA,
                1.0,
                Polarization::Te,
                &[binding],
                || false,
            )
            .is_err()
        );
    }
    let owner = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        1.0,
        Polarization::Te,
        &[ThermoOpticalLayer {
            source: &source,
            card: &card,
        }],
        || false,
    )
    .unwrap();
    let initial = initial(&bound, &[0.0]);
    let before = initial.cells().clone();
    let drive = LayeredWaveDrive::Coherent([C64::new(1.5e6, 0.0), C64::ZERO]);
    let dt = 1e-10;
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    // Establish that the original thermal owner admits the candidate while
    // only the original optical claim's tighter temperature support refuses.
    let frozen = LayeredWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        &[DielectricResponse::OpticalConductor(
            OpticalConductorResponse::try_new(&source).unwrap(),
        )],
        || false,
    )
    .unwrap();
    let candidate = frozen
        .prepare(1.0, Polarization::Te, || false)
        .unwrap()
        .advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let candidate_temperature = candidate
        .thermal
        .temperature_k(bound.cell_ids()[0])
        .unwrap();
    assert!(candidate_temperature > 301.0 && candidate_temperature < 360.0);
    assert!(matches!(
        owner.advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        ),
        Err(LayeredWaveHeatingError::Material(_))
    ));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(source.resolved().state_coordinate("T"), Some(300.0));
}

#[test]
fn mixed_glass_and_metal_refresh_original_spatial_optics_and_native_heat_together() {
    let glass_card = glass_card(360.0);
    let metal_card = optical_card("metal behind glass", 360.0);
    let glass = glass_source(&glass_card);
    let metal = VisibleOpticalStatePoint::Conductor(optical_source(&metal_card, 300.0, 100_000.0));
    let outside = exterior();
    let native = body(&[&glass_card, &metal_card], &[420e-9, 35e-9], &[2, 1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[290.0, 330.0, 300.0]);
    let before = initial.cells().clone();
    let sources = [
        ThermoOpticalLayer {
            source: &glass,
            card: &glass_card,
        },
        ThermoOpticalLayer {
            source: &metal,
            card: &metal_card,
        },
    ];
    let owner: ThermoOpticalWaveHeating<'_, '_, '_, VisibleOpticalStatePoint> =
        ThermoOpticalWaveHeating::try_new(
            &bound,
            &outside,
            DielectricResponse::Lossless(&outside),
            OPTICAL_OMEGA,
            0.9,
            Polarization::Tm,
            &sources,
            || false,
        )
        .unwrap();
    let original_cells = [
        (&glass, &glass_card),
        (&glass, &glass_card),
        (&metal, &metal_card),
    ];
    let refreshed: Vec<_> = bound
        .cell_ids()
        .into_iter()
        .zip(original_cells)
        .map(|(id, (source, card))| {
            source
                .at_query_point(card, &at(initial.temperature_k(id).unwrap(), 100_000.0))
                .unwrap()
        })
        .collect();
    let responses: Vec<_> = refreshed
        .iter()
        .map(|source| {
            DielectricResponse::OpticalMaterial(OpticalMaterialResponse::try_new(source).unwrap())
        })
        .collect();
    let reference = LayeredWaveHeating::try_new_cellwise(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        &responses,
        || false,
    )
    .unwrap();
    let drive = LayeredWaveDrive::Coherent([C64::new(1e6, 0.0), C64::new(2e5, -1e5)]);
    let mut powers = [0.0; 3];
    reference
        .prepare(0.9, Polarization::Tm, || false)
        .unwrap()
        .cell_power_into(drive, &mut powers, || false)
        .unwrap();
    assert!(powers.iter().all(|p| *p > 0.0));
    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    let controls = LayeredThermalSolveControls {
        cell_tolerance_j: 1e-11,
        total_tolerance_j: 1e-11,
        ..LayeredThermalSolveControls::default()
    };
    let dt = 1e-11;
    let result = owner
        .advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(result.applied_optical_states, refreshed);
    assert!(matches!(
        result.applied_optical_states[0],
        VisibleOpticalStatePoint::TabulatedDielectric(_)
    ));
    assert!(matches!(
        result.applied_optical_states[2],
        VisibleOpticalStatePoint::Conductor(_)
    ));
    for (i, (id, energy)) in result.step.cell_energy_j.iter().enumerate() {
        near(*energy, powers[i] * dt, 1e-13);
        let (source, card) = original_cells[i];
        assert_eq!(
            result.endpoint_optical_states[i],
            source
                .at_query_point(
                    card,
                    &at(result.step.thermal.temperature_k(*id).unwrap(), 100_000.0)
                )
                .unwrap()
        );
    }
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.optical.absorbed_j,
        3e-11,
    );
    let next = owner
        .advance_to(
            &result.step.thermal,
            2.0 * dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn glass_endpoint_outside_original_spectral_source_support_refuses_atomically() {
    let card = glass_card(301.0);
    let source = glass_source(&card);
    let outside = exterior();
    let native = body(&[&card], &[420e-9], &[1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[0.0]);
    let before = initial.cells().clone();
    let owner = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        1.0,
        Polarization::Te,
        &[ThermoOpticalLayer {
            source: &source,
            card: &card,
        }],
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let result = owner.advance_to(
        &initial,
        1e-9,
        LayeredWaveDrive::Coherent([C64::new(1.5e6, 0.0), C64::ZERO]),
        conduction(),
        LayeredThermalSolveControls::default(),
        &mut workspace,
        || false,
    );
    assert!(matches!(result, Err(LayeredWaveHeatingError::Material(_))));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(source.resolved().state_coordinate("T"), Some(300.0));
}
