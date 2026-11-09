//! G0/G1/G3/G4: conservative fine gas-network trajectories and source feedback.
use super::*;
use crate::graph::thermal::medium::mass_transfer::evaporation::droplet::chamber::transport::{
    GasDropletTransportCandidate,
    adaptive::{
        AdaptiveGasDropletTransportControls, AdaptiveGasDropletTransportError,
        advance_droplets_with_gas_transport_adaptive,
        advance_droplets_with_gas_transport_and_heat_adaptive,
    },
};
use fs_material::evaporation::gas::{
    GasEvaporationError, droplet::chamber::GasDropletSurfaceStep,
};

fn check_closed(
    initial_droplets: &[SourceGasDropletSurface<'_, '_>],
    initial_gases: &BTreeMap<u64, ClosedGasParcel>,
    final_droplets: &[GasDropletSurfaceStep],
    final_gases: &BTreeMap<u64, ClosedGasParcel>,
    mass_tolerance: f64,
    energy_tolerance: f64,
) {
    assert_eq!(initial_droplets.len(), final_droplets.len());
    assert!(initial_gases.keys().eq(final_gases.keys()));
    let mut liquid_mass = 0.0;
    let mut liquid_energy = 0.0;
    for next in final_droplets {
        let old = initial_droplets.iter().find(|row| row.id == next.id).unwrap();
        liquid_mass += next.droplet.mass_kg - old.droplet.mass_kg;
        liquid_energy += next.droplet.mass_kg
            * (next.droplet.specific_enthalpy_j_kg + old.binding.liquid_enthalpy_offset_j_kg)
            - old.droplet.mass_kg
                * (old.droplet.specific_enthalpy_j_kg + old.binding.liquid_enthalpy_offset_j_kg);
    }
    let vapor = SpeciesId::new("B").unwrap();
    let carrier = SpeciesId::new("A").unwrap();
    let mut gas_mass = 0.0;
    let mut gas_energy = 0.0;
    let mut vapor_mass = 0.0;
    let mut carrier_mass = 0.0;
    for (id, old) in initial_gases {
        let next = &final_gases[id];
        assert_eq!(next.volume_m3(), old.volume_m3());
        assert_eq!(next.point().axes().get("local-lot"), Some(&9.0));
        gas_mass += next.mass_kg() - old.mass_kg();
        gas_energy += next.mass_kg() * next.state().specific_internal_energy_j_kg()
            - old.mass_kg() * old.state().specific_internal_energy_j_kg();
        vapor_mass += next.species_mass_kg(&vapor).unwrap()
            - old.species_mass_kg(&vapor).unwrap();
        carrier_mass += next.species_mass_kg(&carrier).unwrap()
            - old.species_mass_kg(&carrier).unwrap();
    }
    near(liquid_mass + gas_mass, 0.0, mass_tolerance, 0.0);
    near(liquid_mass + vapor_mass, 0.0, mass_tolerance, 0.0);
    near(carrier_mass, 0.0, mass_tolerance, 0.0);
    near(liquid_energy + gas_energy, 0.0, energy_tolerance, 0.0);
}

fn maximum_gas_temperature_error(
    actual: &BTreeMap<u64, ClosedGasParcel>,
    reference: &BTreeMap<u64, ClosedGasParcel>,
) -> f64 {
    actual.iter().map(|(id, gas)| {
        (gas.state().state().temperature - reference[id].state().state().temperature).abs()
    }).fold(0.0, f64::max)
}

// Manual refinement calls only the original full-network interval owner.
// Every step receives the complete previous gas map and the actual droplets.
#[allow(clippy::too_many_arguments)]
fn fixed_reference<'catalog>(
    initial_droplets: &[SourceGasDropletSurface<'_, '_>],
    initial_gases: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    count: usize,
    solve: GasDropletTransportSolve,
    caloric: &LiquidCaloric,
    with_heat: bool,
) -> GasDropletTransportCandidate<'catalog> {
    let mut rows = initial_droplets.to_vec();
    rows.sort_unstable_by_key(|row| row.id);
    let bounds = endpoints(initial_gases);
    let mut gases = initial_gases.clone();
    let mut last = None;
    for _ in 0..count {
        let next = if with_heat {
            let warm: Vec<_> = rows.iter().map(|&surface| SourceGasDropletHeatSurface {
                surface,
                heat_transfer_coefficient_w_m2_k: 80.0,
                heat_controls: heat_controls(caloric),
            }).collect();
            advance_droplets_with_gas_transport_and_heat(
                &warm, REGION, &gases, catalog, paths, &bounds,
                duration_s / count as f64, solve, joint_heat(), || false,
            )
        } else {
            advance_droplets_with_gas_transport(
                &rows, REGION, &gases, catalog, paths, &bounds,
                duration_s / count as f64, solve, || false,
            )
        }.unwrap();
        for (row, step) in rows.iter_mut().zip(next.droplets()) {
            assert_eq!(row.id, step.id);
            row.droplet = step.droplet;
        }
        gases = next.parcels().clone();
        last = Some(next);
    }
    last.unwrap()
}

fn time_controls(span: f64, loose: bool) -> AdaptiveGasDropletTransportControls {
    let mut controls = AdaptiveGasDropletTransportControls::new(interval_controls());
    let time = &mut controls.droplets;
    time.initial_dt_s = span;
    time.minimum_dt_s = 1e-8;
    time.maximum_dt_s = span;
    time.maximum_comparisons = 4096;
    time.maximum_reductions = 32;
    time.enthalpy_absolute_tolerance_j_kg = if loose { 1e4 } else { 1.0 };
    time.temperature_absolute_tolerance_k = if loose { 10.0 } else { 1e-3 };
    time.mass_absolute_tolerance_kg = if loose { 1e-6 } else { 1e-10 };
    time.radius_absolute_tolerance_m = if loose { 1e-3 } else { 1e-8 };
    time.relative_tolerance = 0.0;
    controls
}

#[test]
fn phase_and_heat_accept_only_the_exact_original_fine_network_pair_and_last_source_receipts() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let parcels = BTreeMap::from([
        (REGION, parcel(&catalog, 400.0, 100000.0, 0.003, 1e-5)),
        (DONOR, parcel(&catalog, 440.0, 130000.0, 0.04, 1e-4)),
        (IDLE, parcel(&catalog, 410.0, 90000.0, 0.01, 1e-5)),
    ]);
    let saved = parcels.clone();
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let pinned = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &parcels[&REGION],
        || false,
    ).unwrap();
    let forced = SourceForcedSphericalDropletFilm::try_new(
        &pinned, ForcedSphereMassTransfer::Frossling1938, 1.0, 0.05,
    ).unwrap();
    let temperature = |h| inverse(&caloric, h);
    let rows = [
        SourceGasDropletSurface {
            id: 29,
            law: &LiquidPhase,
            binding: &binding,
            droplet: sphere(&caloric, 410.0, 8e-6),
            temperature_from_enthalpy: &temperature,
            film: SourceDropletFilm::Forced(forced),
        },
        SourceGasDropletSurface {
            id: 11,
            law: &LiquidPhase,
            binding: &binding,
            droplet: sphere(&caloric, 390.0, 1e-6),
            temperature_from_enthalpy: &temperature,
            film: SourceDropletFilm::Stagnant(&pinned),
        },
    ];
    let initial = rows.map(|row| row.droplet);
    let bounds = endpoints(&parcels);
    let links = [nozzle(1e-8)];
    let paths = GasTransportPaths { nozzles: &links, ..Default::default() };
    let solve = interval_controls();
    let controls = time_controls(DURATION, true);
    for with_heat in [false, true] {
        let half = fixed_reference(
            &rows, &parcels, &catalog, paths, DURATION / 2.0, 1, solve,
            &caloric, with_heat,
        );
        let updated: Vec<_> = rows.iter().map(|&row| SourceGasDropletSurface {
            droplet: half.droplets().iter().find(|step| step.id == row.id).unwrap().droplet,
            ..row
        }).collect();
        let fine = fixed_reference(
            &updated, half.parcels(), &catalog, paths, DURATION / 2.0, 1,
            solve, &caloric, with_heat,
        );
        let warm = rows.map(|surface| SourceGasDropletHeatSurface {
            surface,
            heat_transfer_coefficient_w_m2_k: 80.0,
            heat_controls: heat_controls(&caloric),
        });
        let next = if with_heat {
            advance_droplets_with_gas_transport_and_heat_adaptive(
                &warm, REGION, &parcels, &catalog, paths, &bounds,
                DURATION, controls, joint_heat(), || false,
            )
        } else {
            advance_droplets_with_gas_transport_adaptive(
                &rows, REGION, &parcels, &catalog, paths, &bounds,
                DURATION, controls, || false,
            )
        }.unwrap();
        assert_eq!(next.endpoint(), &fine);
        assert_eq!(next.endpoint().films(), fine.films());
        assert_eq!(next.endpoint().report().duration_s, DURATION / 2.0);
        assert_eq!(next.report().duration_s, DURATION);
        assert_eq!(next.report().attempted_steps, 3);
        assert_eq!(next.report().comparisons, 1);
        assert_eq!(next.report().accepted_comparisons, 1);
        assert_eq!(next.report().accepted_intervals, 2);
        assert_eq!(next.report().exchange.accepted_intervals, 2);
        near(next.report().exchange.duration_s, DURATION, 1e-17, 0.0);
        let audit = &next.report().exchange.exchange;
        near(
            audit.transferred_mass_kg,
            half.report().coupled.exchange.transferred_mass_kg
                + fine.report().coupled.exchange.transferred_mass_kg,
            1e-24, 1e-13,
        );
        near(
            audit.sensible_heat_into_liquid_j,
            half.report().coupled.heat.map_or(0.0, |h| h.heat_into_liquid_j)
                + fine.report().coupled.heat.map_or(0.0, |h| h.heat_into_liquid_j),
            1e-16, 1e-13,
        );
        assert_eq!(next.report().exchange.parcels.len(), parcels.len());
        assert!(audit.sum_absolute_mass_residual_kg <= solve.mass_tolerance_kg);
        assert!(audit.sum_absolute_energy_residual_j <= solve.energy_tolerance_j);
        assert!(std::ptr::eq(next.endpoint().catalog(), &catalog));
        assert_eq!(next.endpoint().parcels()[&IDLE], parcels[&IDLE]);
        check_closed(
            &rows, &parcels, next.endpoint().droplets(), next.endpoint().parcels(),
            solve.mass_tolerance_kg, solve.energy_tolerance_j,
        );

        let permuted = if with_heat {
            advance_droplets_with_gas_transport_and_heat_adaptive(
                &[warm[1], warm[0]], REGION, &parcels, &catalog, paths, &bounds,
                DURATION, controls, joint_heat(), || false,
            )
        } else {
            advance_droplets_with_gas_transport_adaptive(
                &[rows[1], rows[0]], REGION, &parcels, &catalog, paths, &bounds,
                DURATION, controls, || false,
            )
        }.unwrap();
        assert_eq!(next.endpoint(), permuted.endpoint());
        assert_eq!(next.report(), permuted.report());
        assert_eq!(parcels, saved);
        assert_eq!(rows.map(|row| row.droplet), initial);
    }
}

#[test]
fn refined_ventilation_tracks_all_finite_gases_against_original_small_intervals() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let parcels = BTreeMap::from([
        (REGION, parcel(&catalog, 400.0, 100000.0, 0.003, 1e-4)),
        (DONOR, parcel(&catalog, 440.0, 130000.0, 0.04, 1e-4)),
        (IDLE, parcel(&catalog, 360.0, 90000.0, 0.001, 1e-4)),
    ]);
    let saved = parcels.clone();
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let pinned = SourceSphericalDropletFilm::new(
        &source(&card, &context).source, &binding.species,
        &parcels[&REGION], || false,
    ).unwrap();
    let temperature = |h| inverse(&caloric, h);
    let rows = [SourceGasDropletSurface {
        id: 11,
        law: &LiquidPhase,
        binding: &binding,
        droplet: sphere(&caloric, 390.0, 1e-6),
        temperature_from_enthalpy: &temperature,
        film: SourceDropletFilm::Stagnant(&pinned),
    }];
    let warm = rows.map(|surface| SourceGasDropletHeatSurface {
        surface,
        heat_transfer_coefficient_w_m2_k: 80.0,
        heat_controls: heat_controls(&caloric),
    });
    let links = [
        nozzle(5e-8),
        SnapshotNozzleLink {
            id: 32,
            parcels: [DONOR, IDLE],
            ..nozzle(5e-8)
        },
    ];
    let paths = GasTransportPaths { nozzles: &links, ..Default::default() };
    let bounds = endpoints(&parcels);
    let duration = 0.4;
    let run = |controls| advance_droplets_with_gas_transport_and_heat_adaptive(
        &warm, REGION, &parcels, &catalog, paths, &bounds,
        duration, controls, joint_heat(), || false,
    ).unwrap();
    let loose = run(time_controls(duration, true));
    let tight = run(time_controls(duration, false));
    let reference = fixed_reference(
        &rows, &parcels, &catalog, paths, duration, 256,
        interval_controls(), &caloric, true,
    );
    let loose_error = maximum_gas_temperature_error(
        loose.endpoint().parcels(), reference.parcels(),
    );
    let tight_error = maximum_gas_temperature_error(
        tight.endpoint().parcels(), reference.parcels(),
    );
    assert!(
        tight_error < 0.35 * loose_error,
        "tight={tight_error:e}, loose={loose_error:e}",
    );
    assert!(tight.report().accuracy_rejections > 0);
    assert!(tight.report().accepted_intervals > loose.report().accepted_intervals);
    assert!(tight.report().maximum_accepted_normalized_error <= 1.0);
    assert_ne!(
        tight.endpoint().parcels()[&DONOR].state().state().temperature,
        parcels[&DONOR].state().state().temperature,
    );
    assert_ne!(
        tight.endpoint().parcels()[&IDLE].state().state().temperature,
        parcels[&IDLE].state().state().temperature,
    );
    assert_ne!(
        tight.endpoint().parcels()[&IDLE].species_mass_kg(&binding.species).unwrap(),
        parcels[&IDLE].species_mass_kg(&binding.species).unwrap(),
    );
    let last = &tight.endpoint().films()[0];
    let SourceDropletFilmSample::Stagnant(final_sample) = &last.films[1] else {
        panic!("the explicitly stagnant droplet changed its transport model");
    };
    let refreshed = pinned.resolve(
        tight.endpoint().droplets()[0].geometry[1],
        &tight.endpoint().parcels()[&REGION],
        || false,
    ).unwrap();
    assert_eq!(final_sample, &refreshed);
    assert!(tight.report().exchange.exchange.transferred_mass_kg > 0.0);
    assert!(tight.report().exchange.exchange.sensible_heat_into_liquid_j > 0.0);
    check_closed(
        &rows, &parcels, tight.endpoint().droplets(), tight.endpoint().parcels(),
        interval_controls().mass_tolerance_kg, interval_controls().energy_tolerance_j,
    );
    assert_eq!(parcels, saved);
}

#[test]
fn an_empty_droplet_roster_still_refines_remote_only_gas_transport() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let parcels = BTreeMap::from([
        (REGION, parcel(&catalog, 400.0, 100000.0, 0.003, 1e-5)),
        (DONOR, parcel(&catalog, 440.0, 130000.0, 0.04, 1e-4)),
        (IDLE, parcel(&catalog, 360.0, 90000.0, 0.001, 1e-4)),
    ]);
    let saved = parcels.clone();
    let links = [SnapshotNozzleLink {
        id: 32,
        parcels: [DONOR, IDLE],
        ..nozzle(1e-7)
    }];
    let paths = GasTransportPaths { nozzles: &links, ..Default::default() };
    let bounds = endpoints(&parcels);
    let duration = 0.4;
    let controls = time_controls(duration, false);
    let next = advance_droplets_with_gas_transport_adaptive(
        &[], REGION, &parcels, &catalog, paths, &bounds,
        duration, controls, || false,
    ).unwrap();
    let reference = fixed_reference(
        &[], &parcels, &catalog, paths, duration, 256,
        interval_controls(), &caloric, false,
    );
    let loose = advance_droplets_with_gas_transport_adaptive(
        &[], REGION, &parcels, &catalog, paths, &bounds,
        duration, time_controls(duration, true), || false,
    ).unwrap();
    // The selected atmosphere is completely static. Only the other two gas
    // nodes can cause these accuracy rejections and refine the actual flow.
    assert_eq!(next.endpoint().parcels()[&REGION], parcels[&REGION]);
    assert!(next.endpoint().droplets().is_empty());
    assert!(next.endpoint().films().is_empty());
    assert!(next.report().accuracy_rejections > 0);
    assert!(next.report().accepted_comparisons > 1);
    assert!(next.report().accepted_intervals > loose.report().accepted_intervals);
    assert_ne!(next.endpoint().parcels()[&DONOR], parcels[&DONOR]);
    assert_ne!(next.endpoint().parcels()[&IDLE], parcels[&IDLE]);
    assert!(
        maximum_gas_temperature_error(next.endpoint().parcels(), reference.parcels())
            < 0.35 * maximum_gas_temperature_error(loose.endpoint().parcels(), reference.parcels()),
    );
    assert_eq!(next.report().exchange.exchange.transferred_mass_kg, 0.0);
    assert_eq!(next.report().exchange.exchange.sensible_heat_into_liquid_j, 0.0);
    assert_eq!(next.report().exchange.parcels.len(), parcels.len());
    near(next.report().exchange.duration_s, duration, 1e-14, 0.0);
    check_closed(
        &[], &parcels, next.endpoint().droplets(), next.endpoint().parcels(),
        interval_controls().mass_tolerance_kg, interval_controls().energy_tolerance_j,
    );
    assert!(matches!(
        advance_droplets_with_gas_transport_adaptive(
            &[], REGION, &parcels, &catalog, paths, &bounds, 0.0, controls, || false,
        ),
        Err(AdaptiveGasDropletTransportError::Invalid(_)),
    ));
    assert_eq!(parcels, saved);
}

#[test]
fn bounded_phase_retry_late_cancellation_and_terminal_source_refusal_preserve_original_inputs() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let parcels = BTreeMap::from([
        (REGION, parcel(&catalog, 400.0, 100000.0, 0.003, 1e-5)),
        (DONOR, parcel(&catalog, 440.0, 130000.0, 0.04, 1e-4)),
        (IDLE, parcel(&catalog, 410.0, 90000.0, 0.01, 1e-5)),
    ]);
    let saved = parcels.clone();
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let pinned = SourceSphericalDropletFilm::new(
        &source(&card, &context).source, &binding.species,
        &parcels[&REGION], || false,
    ).unwrap();
    let temperature = |h| inverse(&caloric, h);
    let rows = [SourceGasDropletSurface {
        id: 11,
        law: &LiquidPhase,
        binding: &binding,
        droplet: sphere(&caloric, 390.0, 1e-6),
        temperature_from_enthalpy: &temperature,
        film: SourceDropletFilm::Stagnant(&pinned),
    }];
    let initial = rows[0].droplet;
    let links = [nozzle(1e-8)];
    let paths = GasTransportPaths { nozzles: &links, ..Default::default() };
    let bounds = endpoints(&parcels);
    let coarse = advance_droplets_with_gas_transport(
        &rows, REGION, &parcels, &catalog, paths, &bounds,
        DURATION, interval_controls(), || false,
    ).unwrap();
    let requested_fraction = coarse.droplets()[0].phase.transferred_mass_kg
        / rows[0].droplet.mass_kg;
    assert!(requested_fraction > 0.0);
    let mut physical = time_controls(DURATION, true);
    // The original owner supplies this fixture's actual withdrawal scale.
    // A tighter authored limit must cause time refinement, not clipped mass.
    physical.droplets.solve.surface.maximum_inventory_fraction = requested_fraction * 0.5;
    let refined = advance_droplets_with_gas_transport_adaptive(
        &rows, REGION, &parcels, &catalog, paths, &bounds,
        DURATION, physical, || false,
    ).unwrap();
    assert!(refined.report().physical_rejections > 0);
    assert!(refined.report().accepted_intervals > 2);
    assert_eq!(refined.report().duration_s, DURATION);
    check_closed(
        &rows, &parcels, refined.endpoint().droplets(), refined.endpoint().parcels(),
        interval_controls().mass_tolerance_kg, interval_controls().energy_tolerance_j,
    );

    let fixed = time_controls(DURATION, true);
    let mut exhausted = fixed;
    exhausted.droplets.maximum_comparisons = 1;
    assert!(matches!(
        advance_droplets_with_gas_transport_adaptive(
            &rows, REGION, &parcels, &catalog, paths, &bounds,
            2.0 * DURATION, exhausted, || false,
        ),
        Err(AdaptiveGasDropletTransportError::ComparisonBudget),
    ));
    let mut calls = 0;
    let accepted = advance_droplets_with_gas_transport_adaptive(
        &rows, REGION, &parcels, &catalog, paths, &bounds,
        2.0 * DURATION, fixed, || { calls += 1; false },
    ).unwrap();
    let final_gate = calls;
    let mut calls = 0;
    let cancelled = advance_droplets_with_gas_transport_adaptive(
        &rows, REGION, &parcels, &catalog, paths, &bounds,
        2.0 * DURATION, fixed, || {
            calls += 1;
            calls == final_gate
        },
    );
    assert!(matches!(
        cancelled,
        Err(AdaptiveGasDropletTransportError::Material(GasEvaporationError::Cancelled))
            | Err(AdaptiveGasDropletTransportError::Step {
                error: DropletGasTransportError::Cancelled, ..
            }),
    ));
    let retry = advance_droplets_with_gas_transport_adaptive(
        &rows, REGION, &parcels, &catalog, paths, &bounds,
        2.0 * DURATION, fixed, || false,
    ).unwrap();
    assert_eq!(retry.endpoint(), accepted.endpoint());
    assert_eq!(retry.report(), accepted.report());

    let restricted_card = material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "temperature",
        vec![(200.0, 1.6e-5), (400.0001, 1.6e-5)], [200.0, 400.0001],
    )]);
    let restricted = SourceSphericalDropletFilm::new(
        &source(&restricted_card, &context).source, &binding.species,
        &parcels[&REGION], || false,
    ).unwrap();
    let selected = restricted.selected_diffusivity().clone();
    let unsupported = SourceGasDropletSurface {
        film: SourceDropletFilm::Stagnant(&restricted),
        ..rows[0]
    };
    // Original gas transport leaves this pinned diffusivity's T domain. The
    // exact original source error is terminal at the first full trial span.
    assert!(matches!(
        advance_droplets_with_gas_transport_adaptive(
            &[unsupported], REGION, &parcels, &catalog, paths, &bounds,
            DURATION, fixed, || false,
        ),
        Err(AdaptiveGasDropletTransportError::Step {
            time_s, duration_s,
            error: DropletGasTransportError::Droplets(
                SphericalDropletFilmError::Film(GasFilmError::Source(_))
            ),
        }) if time_s == 0.0 && duration_s == DURATION,
    ));
    assert_eq!(restricted.selected_diffusivity(), &selected);
    assert_eq!(rows[0].droplet, initial);
    assert_eq!(parcels, saved);
}
