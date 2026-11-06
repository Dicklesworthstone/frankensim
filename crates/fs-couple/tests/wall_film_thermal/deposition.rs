//! Native spatial storage, phase heat and one finite gas share one endpoint.
use super::fixture::*;
use std::{cell::Cell, collections::BTreeMap};

use fs_blake3::ContentHash;
use fs_couple::finite_gas_thermal::wall_film::{
    deposition::{
        DepositionWallFilmBinding, DepositionWallFilmControls, advance_deposition_with_wall_film,
    },
    implicit::WallFilmImplicitControls,
};
use fs_material::{
    evaporation::gas::{
        GasLiquidBinding, implicit::ImplicitGasError, wall_film::GasWallFilmInventory,
    },
    phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve},
    thermochemical_gas::energy::parcel::ClosedGasParcel,
};
use fs_thermal::{
    enthalpy::{
        Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
        DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
        backward_euler::{
            BoundaryExchangeStepCoefficients, BoundaryStepCoefficients, StepCoefficients,
            TransportStepCoefficients,
        },
    },
    evaporation::EvaporationLaw,
};

pub(super) const MASS: f64 = 0.002;
pub(super) const CONTACT_AREA: f64 = 0.001;
pub(super) const LINK: u64 = 9;

pub(super) fn cell(id: u64, material: u64, h: f64, area: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id + 10,
        body_id: 1,
        material_id: MaterialId(material),
        mass_kg: MASS,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 2e-6,
            surface_area_m2: area,
        },
        boundary: Boundary {
            ambient_k: 280.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
pub(super) fn contact(g: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(1),
            b: CellId(2),
            area_m2: CONTACT_AREA,
            conductance_w_k: g,
        }),
    }
}
pub(super) fn network<'a>(
    cells: Vec<DepositedCell>,
    laws: &[(MaterialId, &'a dyn ThermalLaw)],
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    events.extend(later);
    DepositionSimulation::new(
        events,
        laws,
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
pub(super) fn coefficients<T>(payload: T) -> BoundaryExchangeStepCoefficients<T> {
    TransportStepCoefficients::from(BoundaryStepCoefficients::from(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload,
    }))
    .into()
}
pub(super) fn empty_heat(
    _: &DepositionSimulation<'_>,
    _: &ClosedGasParcel,
) -> Result<BoundaryExchangeStepCoefficients<()>, ImplicitGasError> {
    Ok(coefficients(()))
}
pub(super) fn patch(
    binding: &GasLiquidBinding,
    film: GasWallFilmInventory,
    area: f64,
    htc: f64,
) -> DepositionWallFilmBinding<'_> {
    let mut boundary = boundary(binding);
    boundary.area_m2 = area;
    DepositionWallFilmBinding {
        cell_id: CellId(1),
        film,
        boundary,
        heat_transfer_coefficient_w_m2_k: htc,
    }
}
pub(super) fn policy() -> DepositionWallFilmControls {
    let mut result = DepositionWallFilmControls::new(controls().phase);
    result.implicit.thermal.newton.absolute_tolerance = 1e-11;
    result.implicit.thermal.newton.relative_tolerance = 0.0;
    result
}
pub(super) fn balances(
    old: &DepositionSimulation<'_>,
    gas: &ClosedGasParcel,
    film: GasWallFilmInventory,
    binding: &GasLiquidBinding,
    next: (
        &DepositionSimulation<'_>,
        &ClosedGasParcel,
        GasWallFilmInventory,
    ),
    external: f64,
) {
    let (thermal, next_gas, next_film) = next;
    let wall_delta: f64 = old
        .cells()
        .iter()
        .map(|(id, cell)| {
            assert_eq!(thermal.cells()[id].mass_kg, cell.mass_kg);
            cell.mass_kg * (thermal.cells()[id].specific_h_j_kg - cell.specific_h_j_kg)
        })
        .sum();
    let gas_delta = next_gas.mass_kg() * next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    near(
        wall_delta + gas_delta + film_energy(next_film, binding) - film_energy(film, binding),
        external,
        1e-7,
    );
    near(
        next_gas.species_mass_kg(&water()).unwrap() + next_film.mass_kg(),
        gas.species_mass_kg(&water()).unwrap() + film.mass_kg(),
        1e-13,
    );
    near(
        next_gas.species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-17,
    );
    near(
        next_gas.mass_kg() + next_film.mass_kg(),
        gas.mass_kg() + film.mass_kg(),
        1e-13,
    );
    assert_eq!(next_gas.volume_m3(), GAS_VOLUME);
    assert_eq!(next_gas.point().axes().get("wall-film-lot"), Some(&13.0));
}

#[test]
fn spatial_neighbor_and_native_wall_freezing_supply_full_film_evaporation() {
    const WALL_LATENT: f64 = 100_000.0;
    let plateau = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x87; 32]),
        [
            (280.0, -20.0 * CP_WALL, 0.0),
            (300.0, 0.0, 0.0),
            (300.0, WALL_LATENT, 1.0),
            (340.0, WALL_LATENT + 40.0 * CP_WALL, 1.0),
        ]
        .map(
            |(temperature_k, specific_enthalpy_j_kg, liquid_mass_fraction)| EnthalpyPhaseKnot {
                temperature_k,
                specific_enthalpy_j_kg,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1000.0,
            },
        )
        .to_vec(),
    )
    .unwrap();
    let sensible = curve();
    let thermal = network(
        vec![
            cell(1, 1, 0.5 * WALL_LATENT, AREA + CONTACT_AREA),
            cell(2, 2, 10.0 * CP_WALL, AREA + CONTACT_AREA),
        ],
        &[(MaterialId(1), &plateau), (MaterialId(2), &sensible)],
        vec![contact(0.1)],
    );
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let binding = binding(&source);
    let film = GasWallFilmInventory::Wet {
        mass_kg: 1e-7,
        specific_enthalpy_j_kg: 0.0,
    };
    let next = advance_deposition_with_wall_film(
        &thermal,
        patch(
            &binding,
            film,
            thermal.exposed_area_m2(CellId(1)).unwrap(),
            0.0,
        ),
        &gas,
        &source,
        0.1,
        0.0,
        policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let capacity = MASS * CP_WALL;
    let neighbor_t = (capacity * 310.0 + 0.01 * 300.0) / (capacity + 0.01);
    let conducted = capacity * (310.0 - neighbor_t);
    let released = film.mass_kg() * LATENT;
    near(
        next.thermal().temperature_k(CellId(2)).unwrap(),
        neighbor_t,
        1e-7,
    );
    assert_eq!(next.thermal().temperature_k(CellId(1)).unwrap(), 300.0);
    let native = plateau
        .state_at_specific_enthalpy(next.thermal().cells()[&CellId(1)].specific_h_j_kg)
        .unwrap();
    near(
        native.specific_enthalpy_j_kg(),
        0.5 * WALL_LATENT + (conducted - released) / MASS,
        1e-6,
    );
    near(
        native.liquid_mass_fraction(),
        0.5 + (conducted - released) / (MASS * WALL_LATENT),
        1e-11,
    );
    assert!(conducted > 0.0 && conducted < released);
    assert_eq!(next.film(), GasWallFilmInventory::Dry);
    assert_eq!(next.report().phase.transferred_mass_kg, film.mass_kg());
    near(next.report().wall_heat_into_phase_system_j, released, 1e-9);
    let cg = 2.5 * 100_000.0 * GAS_VOLUME / 290.0;
    let tg =
        290.0 + film.mass_kg() * (R * 290.0 + 3.5 * R * 10.0) / (cg + film.mass_kg() * 2.5 * R);
    near(next.gas().state().state().temperature, tg, 2e-7);
    assert_eq!(next.report().sensible_heat_into_gas_j, 0.0);
    near(next.report().external_thermal_heat_j, 0.0, 1e-12);
    balances(
        &thermal,
        &gas,
        film,
        &binding,
        (next.thermal(), next.gas(), next.film()),
        0.0,
    );
}

#[test]
fn dry_condensation_warms_a_spatial_wall_and_matches_the_original_uniform_limit() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let one = network(
        vec![cell(1, 1, 0.0, AREA)],
        &[(MaterialId(1), &curve)],
        vec![],
    );
    let interface = patch(&binding, GasWallFilmInventory::Dry, AREA, 0.0);
    let uniform = advance_deposition_with_wall_film(
        &one,
        interface,
        &gas,
        &source,
        0.1,
        0.0,
        policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let body = body(0.0);
    let lumped = scope(|cx, _| {
        environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body)
            .advance_implicit(
                cx,
                input(&curve, 300.0, 0.1),
                WallFilmImplicitControls::default(),
            )
            .unwrap()
    });
    near(
        uniform.thermal().temperature_k(CellId(1)).unwrap(),
        lumped.state.temperature_k(),
        1e-7,
    );
    near(
        uniform.gas().state().state().temperature,
        lumped.report.next_gas.state().state().temperature,
        1e-7,
    );
    near(
        uniform.film().mass_kg(),
        lumped.report.next_film.mass_kg(),
        1e-14,
    );

    let thermal = network(
        vec![
            cell(1, 1, 0.0, AREA + CONTACT_AREA),
            cell(2, 1, -10.0 * CP_WALL, AREA + CONTACT_AREA),
        ],
        &[(MaterialId(1), &curve)],
        vec![contact(0.1)],
    );
    let interface = patch(
        &binding,
        GasWallFilmInventory::Dry,
        thermal.exposed_area_m2(CellId(1)).unwrap(),
        0.0,
    );
    let next = advance_deposition_with_wall_film(
        &thermal,
        interface,
        &gas,
        &source,
        0.1,
        0.0,
        policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let tw = next.thermal().temperature_k(CellId(1)).unwrap();
    let tn = next.thermal().temperature_k(CellId(2)).unwrap();
    let tg = next.gas().state().state().temperature;
    let condensed = next.film().mass_kg();
    assert!(next.report().phase.film_formed && condensed > 0.0);
    assert!(tn > 290.0 && tn < tw);
    assert!(tw < uniform.thermal().temperature_k(CellId(1)).unwrap());
    let conducted = 0.01 * (tw - tn);
    near(MASS * CP_WALL * (tn - 290.0), conducted, 1e-8);
    near(
        (MASS * CP_WALL + condensed * CP_LIQUID) * (tw - 300.0),
        condensed * (LATENT + 3.5 * R * (tg - 300.0)) - conducted,
        1e-8,
    );
    // Vapor leaves a fixed-volume gas at its actual endpoint enthalpy.
    let cg = 2.5 * 100_000.0 * GAS_VOLUME / 320.0;
    near(tg, cg * 320.0 / (cg + condensed * R), 2e-7);
    balances(
        &thermal,
        &gas,
        GasWallFilmInventory::Dry,
        &binding,
        (next.thermal(), next.gas(), next.film()),
        0.0,
    );
}

#[test]
fn partial_patch_keeps_remaining_ambient_area_and_endpoint_source_contact_feedback() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let mut wall = cell(1, 1, 0.0, AREA + CONTACT_AREA);
    wall.boundary.htc_w_m2_k = 20.0;
    wall.boundary.source_w = 0.7;
    let thermal = network(
        vec![wall, cell(2, 1, 10.0 * CP_WALL, AREA + CONTACT_AREA)],
        &[(MaterialId(1), &curve)],
        vec![contact(0.5)],
    );
    let dt = 0.1;
    let (tw, tn, tg, condensed) = (301.0, 309.5, 319.0, 1e-7);
    let area = 0.5 * AREA;
    let carrier_r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
    let vapor_mass = gas.species_mass_kg(&water()).unwrap() - condensed;
    let pressure =
        (gas.species_mass_kg(&carrier()).unwrap() * carrier_r + vapor_mass * R) * tg / GAS_VOLUME;
    let partial = vapor_mass * R * tg / GAS_VOLUME;
    let saturation = PHASE.properties(tw).unwrap().saturation_pressure_pa;
    let unit_flux = pressure / (R * tg) * ((pressure - partial) / (pressure - saturation)).ln();
    let mut interface = patch(&binding, GasWallFilmInventory::Dry, area, 20.0);
    interface.boundary.mass_transfer_coefficient_m_s = -condensed / (dt * area * unit_flux);
    let sensible_gas = dt * area * 20.0 * (tw - tg);
    let cg = 2.5 * 100_000.0 * GAS_VOLUME / 320.0;
    let gas_heat = cg * (tg - 320.0) + condensed * R * tg - sensible_gas;
    let contact_heat = MASS * CP_WALL * (310.0 - tn);
    let conductance = contact_heat / (dt * (tn - tw));
    let ambient_loss = dt * (AREA - area) * 20.0 * (tw - 280.0);
    let phase_heat = -condensed * (LATENT + 3.5 * R * (tg - 300.0) - CP_LIQUID * (tw - 300.0));
    let source_power =
        (MASS * CP_WALL * (tw - 300.0) - contact_heat + ambient_loss + sensible_gas + phase_heat)
            / dt
            - 0.7;
    let next = advance_deposition_with_wall_film(
        &thermal,
        interface,
        &gas,
        &source,
        dt,
        gas_heat,
        policy(),
        |thermal, gas| {
            let t = thermal.temperature_k(CellId(1))?;
            let neighbor = thermal.temperature_k(CellId(2))?;
            let gas_t = gas.state().state().temperature;
            let g = conductance * (1.0 + 0.001 * (t - tw) + 0.001 * (neighbor - tn));
            let power = source_power + 0.02 * (t - tw) + 0.01 * (gas_t - tg);
            let mut result = coefficients([t, neighbor, gas_t, g, power]);
            result
                .transport
                .heat
                .bulk
                .sources_w
                .insert(CellId(1), power);
            result.transport.heat.bulk.conductances_w_k.insert(LINK, g);
            Ok(result)
        },
        || Ok(()),
    )
    .unwrap();
    near(next.thermal().temperature_k(CellId(1)).unwrap(), tw, 1e-7);
    near(next.thermal().temperature_k(CellId(2)).unwrap(), tn, 1e-7);
    near(next.gas().state().state().temperature, tg, 1e-7);
    near(next.film().mass_kg(), condensed, 1e-14);
    near(next.report().sensible_heat_into_gas_j, sensible_gas, 1e-8);
    near(
        next.report().wall_heat_into_phase_system_j,
        phase_heat,
        1e-8,
    );
    let external_thermal = dt * (source_power + 0.7) - ambient_loss;
    near(
        next.report().external_thermal_heat_j,
        external_thermal,
        1e-8,
    );
    for (actual, expected) in next
        .payload()
        .iter()
        .zip([tw, tn, tg, conductance, source_power])
    {
        near(*actual, expected, 1e-7);
    }
    // Step-local constitutive results do not overwrite authored boundaries.
    assert_eq!(next.thermal().cells()[&CellId(1)].boundary.htc_w_m2_k, 20.0);
    assert_eq!(next.thermal().cells()[&CellId(1)].boundary.source_w, 0.7);
    assert_eq!(next.thermal().links()[&LINK].conductance_w_k, 0.5);
    balances(
        &thermal,
        &gas,
        GasWallFilmInventory::Dry,
        &binding,
        (next.thermal(), next.gas(), next.film()),
        external_thermal + gas_heat,
    );
}

#[test]
fn pending_events_occupied_area_late_sources_work_and_cancellation_preserve_inputs() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let initial = cell(1, 1, 0.0, AREA);
    let boundary_event = TimedEvent {
        time_s: 0.1,
        event: Event::Boundary(
            CellId(1),
            Boundary {
                source_w: 4.0,
                ..initial.boundary
            },
        ),
    };
    let thermal = network(
        vec![initial],
        &[(MaterialId(1), &curve)],
        vec![boundary_event],
    );
    let original_cells = thermal.cells().clone();
    let interface = patch(&binding, GasWallFilmInventory::Dry, AREA, 0.0);
    assert!(
        advance_deposition_with_wall_film(
            &thermal,
            interface,
            &gas,
            &source,
            0.2,
            0.0,
            policy(),
            empty_heat,
            || Ok(())
        )
        .is_err()
    );
    let queries = Cell::new(0usize);
    let polls = Cell::new(0usize);
    let accepted = advance_deposition_with_wall_film(
        &thermal,
        interface,
        &gas,
        &source,
        0.1,
        0.0,
        policy(),
        |_, _| {
            queries.set(queries.get() + 1);
            Ok(coefficients(()))
        },
        || {
            polls.set(polls.get() + 1);
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(accepted.thermal().time_s(), 0.1);
    assert_eq!(
        accepted.thermal().cells()[&CellId(1)].boundary.source_w,
        0.0
    );
    assert!(
        advance_deposition_with_wall_film(
            accepted.thermal(),
            patch(&binding, accepted.film(), AREA, 0.0),
            accepted.gas(),
            &source,
            0.2,
            0.0,
            policy(),
            empty_heat,
            || Ok(())
        )
        .is_err()
    );
    let processed = accepted
        .thermal()
        .advance_to(
            0.1,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-6,
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(processed.1.events, 1);
    assert_eq!(processed.0.cells()[&CellId(1)].boundary.source_w, 4.0);

    assert!(
        advance_deposition_with_wall_film(
            &thermal,
            patch(&binding, GasWallFilmInventory::Dry, 2.0 * AREA, 0.0),
            &gas,
            &source,
            0.1,
            0.0,
            policy(),
            empty_heat,
            || Ok(())
        )
        .is_err()
    );
    assert!(
        advance_deposition_with_wall_film(
            &thermal,
            interface,
            &gas,
            &source,
            0.1,
            0.0,
            policy(),
            |_, _| {
                let mut result = coefficients(());
                result
                    .transport
                    .heat
                    .boundaries
                    .ambient_w_k
                    .insert(CellId(1), 0.1);
                Ok(result)
            },
            || Ok(())
        )
        .is_err()
    );
    let mut limited = policy();
    limited.implicit.thermal.maximum_evaluations = accepted.report().residual_evaluations - 1;
    assert!(
        advance_deposition_with_wall_film(
            &thermal,
            interface,
            &gas,
            &source,
            0.1,
            0.0,
            limited,
            empty_heat,
            || Ok(())
        )
        .is_err()
    );

    let last_query = queries.get();
    queries.set(0);
    assert!(matches!(
        advance_deposition_with_wall_film(
            &thermal,
            interface,
            &gas,
            &source,
            0.1,
            0.0,
            policy(),
            |_, _| {
                queries.set(queries.get() + 1);
                if queries.get() == last_query {
                    Err(ImplicitGasError::Invalid(
                        "last spatial heat source refused",
                    ))
                } else {
                    Ok(coefficients(()))
                }
            },
            || Ok(())
        ),
        Err(ImplicitGasError::Invalid(
            "last spatial heat source refused"
        ))
    ));
    assert_eq!(queries.get(), last_query);
    let last_poll = polls.get();
    polls.set(0);
    assert!(
        advance_deposition_with_wall_film(
            &thermal,
            interface,
            &gas,
            &source,
            0.1,
            0.0,
            policy(),
            empty_heat,
            || {
                polls.set(polls.get() + 1);
                if polls.get() == last_poll {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        )
        .is_err()
    );
    assert_eq!(polls.get(), last_poll);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(gas, original_gas);
    assert_eq!(interface.film, GasWallFilmInventory::Dry);
}
