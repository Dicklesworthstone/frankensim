//! Original loading clocks, finite thermal transport, and contact disappearance.

use super::{MASS, binding, close, driver, mixed_source};
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode,
    combined::{JkrCombinedBoundary, loading::JkrCombinedLoadingState},
};
use fs_couple::friction_heat::{
    FrictionHeatRoute,
    adhesive::snap::combined_loading::{
        JkrCombinedHeat,
        timed::{
            JkrCombinedBoundaryPolicy as Policy, JkrTimedCombinedHeatControls,
            footprint::JkrCombinedHeatFootprint,
        },
    },
};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLink, TimedEvent,
};

const LINK: u64 = 31;
const STEP: f64 = 1.0 / 16.0;
const FOOTPRINT_STEP: f64 = 0.25;

fn controls() -> JkrTimedCombinedHeatControls {
    JkrTimedCombinedHeatControls {
        heat: super::controls().heat,
        transport: StepConfig {
            max_dt_s: STEP,
            energy_tolerance_j: 2e-22,
        },
        maximum_transport_substeps: 64,
        maximum_footprint_intervals: 16,
        maximum_footprint_dt_s: FOOTPRINT_STEP,
        total_energy_tolerance_j: 1e-20,
    }
}

// Manufactured constant finite heat capacity, bound to the original two cards.
fn caloric(source: &SourceJkrMixedMode, capacity_j_k: f64) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let elastic = source.normal_source().normal_state().elastic();
    [
        elastic.surface_a_state().resolved().card_identity(),
        elastic.surface_b_state().resolved().card_identity(),
    ]
    .map(|card| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card,
            fs_material::phase::SolidLiquidPhase::Solid,
            [0.0, 1000.0]
                .map(|rise| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: rise * capacity_j_k / MASS,
                    temperature_k: 290.0 + rise,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                })
                .to_vec(),
        )
        .unwrap()
    })
}

fn geometry(face_area_m2: f64, side: usize) -> CellGeometry {
    let volume = MASS / 1000.0;
    let length = volume / face_area_m2;
    CellGeometry {
        start_m: [if side == 0 { -length } else { 0.0 }, 0.0, 0.0],
        end_m: [if side == 0 { 0.0 } else { length }, 0.0, 0.0],
        volume_m3: volume,
        surface_area_m2: 2.0 * face_area_m2 + 4.0 * face_area_m2.sqrt() * length,
    }
}

fn ambient(capacity: f64, full_area: f64, cooling_rate_per_s: f64) -> Boundary {
    Boundary {
        ambient_k: 290.0,
        htc_w_m2_k: cooling_rate_per_s * capacity / full_area,
        reservoir: None,
        source_w: 0.0,
    }
}

fn timeline<'m>(
    caloric: &'m [EquilibriumEnthalpyPhaseCurve; 2],
    initial: JkrCombinedLoadingState,
    capacity: f64,
    face_area: f64,
    initial_rise: [f64; 2],
    cooling_rate: f64,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'m> {
    let mut schedule: Vec<_> = (0..2)
        .map(|side| {
            let geometry = geometry(face_area, side);
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(side as u64),
                    segment_id: side as u64,
                    body_id: side as u64,
                    material_id: MaterialId(side as u64),
                    mass_kg: MASS,
                    specific_h_j_kg: initial_rise[side] * capacity / MASS,
                    geometry,
                    boundary: ambient(capacity, geometry.surface_area_m2, cooling_rate),
                }),
            }
        })
        .collect();
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(0),
            b: CellId(1),
            area_m2: initial.contact().unwrap().contact_area_m2,
            conductance_w_k: capacity,
        }),
    });
    schedule.extend(later);
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &caloric[0]), (MaterialId(1), &caloric[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 16,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(0.0, controls().transport, || Ok(()))
    .unwrap()
    .0
}

fn footprint(
    heat: &JkrCombinedHeat<'_>,
    thermal: &DepositionSimulation<'_>,
) -> JkrCombinedHeatFootprint {
    let cell = &thermal.cells()[&CellId(0)];
    let length = cell.geometry.end_m[0] - cell.geometry.start_m[0];
    let link = thermal.links()[&LINK];
    // Two original equal half-lengths in series: G=A*k/L.
    let k = link.conductance_w_k * length / link.area_m2;
    let route = FrictionHeatRoute::from_link(thermal, LINK, [k; 2]).unwrap();
    JkrCombinedHeatFootprint::new(heat, thermal, route).unwrap()
}

#[test]
fn boundary_release_heats_at_its_clock_then_exposes_both_faces_to_native_cooling() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let limit = loading.limit().unwrap();
    let original_fold = loading
        .advance_to_load_factor(initial, 2.0 * limit.load_factor, &gate)
        .unwrap();
    let release = original_fold.state().stored_energy_j();
    let caloric = caloric(&source, release);
    let face = 2.0 * initial.contact().unwrap().contact_area_m2;
    let surface = geometry(face, 0).surface_area_m2;
    let change = (0..2)
        .map(|side| TimedEvent {
            time_s: 0.5,
            event: Event::Boundary(CellId(side), ambient(release, surface, 2.0)),
        })
        .collect();
    let mut thermal = timeline(&caloric, initial, release, face, [0.0; 2], 1.0, change);
    let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
    let mut footprint = footprint(&heat, &thermal);
    let ramp = loading
        .advance_linear_ramp(initial, 0.0, 1.0, 2.0 * limit.load_factor, 1e-15, &gate)
        .unwrap();
    assert_eq!(ramp.accepted_final_time_s(), 0.5);
    let report = heat
        .apply_timed_ramp(
            &mut thermal,
            heat.stamp(),
            &ramp,
            &mut footprint,
            Policy::ReleaseIfSeparatingFold,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.boundary_time_s, Some(0.5));
    assert_eq!(report.release_time_s, Some(0.5));
    assert_eq!(report.mechanics[0].released_energy_j, 0.0);
    close(report.mechanics[1].released_energy_j, release);
    assert_eq!(report.final_time_s, 1.0);
    assert_eq!(thermal.time_s(), 1.0);
    assert_eq!(report.transport.substeps, 16);
    assert_eq!(report.transport.events, 2);
    assert_eq!(heat.state().contact(), None);
    assert_eq!(heat.state().load_factor(), 2.0 * limit.load_factor);
    assert_eq!(footprint.area_m2(), None);
    assert!(!thermal.links().contains_key(&LINK));
    close(heat.ledger().dissipated_work_j(), release);
    let cooling = (1.0 - 2.0 * STEP).powi(8);
    for (id, share) in [(CellId(0), 0.25), (CellId(1), 0.5)] {
        close(thermal.exposed_area_m2(id).unwrap(), surface);
        close(
            MASS * thermal.cells()[&id].specific_h_j_kg,
            share * release * cooling,
        );
    }
    close(
        report.transport.external_heat_j,
        0.75 * release * (cooling - 1.0),
    );
    close(
        report.thermal.stored_energy_change_j,
        0.75 * release * cooling,
    );
    assert!(report.combined_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
    assert!(report.absolute_stage_energy_residual_j <= controls().total_energy_tolerance_j);
}

#[test]
fn held_and_reversing_contacts_conduct_with_the_original_area_and_guard_the_peak_area() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.6, &gate).unwrap();
    let initial_area = initial.contact().unwrap().contact_area_m2;
    let zero_area = loading
        .initial_state(0.0, &gate)
        .unwrap()
        .contact()
        .unwrap()
        .contact_area_m2;
    let capacity = initial.stored_energy_j().abs();
    let caloric = caloric(&source, capacity);
    for endpoint in [0.6, -0.6] {
        let mut thermal = timeline(
            &caloric,
            initial,
            capacity,
            2.0 * zero_area,
            [1.0, 0.0],
            0.0,
            vec![],
        );
        let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
        let mut footprint = footprint(&heat, &thermal);
        let ramp = loading
            .advance_linear_ramp(initial, 0.0, 1.0, endpoint, 1e-15, &gate)
            .unwrap();
        let report = heat
            .apply_timed_ramp(
                &mut thermal,
                heat.stamp(),
                &ramp,
                &mut footprint,
                Policy::Stop,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        // Equal finite capacities have conserved sum and decaying difference.
        // Integrate the declared left-footprint split independently from the
        // original source areas; there is no contact heat production here.
        let mut difference = capacity;
        for i in 0..4 {
            let load = 0.6 + (endpoint - 0.6) * i as f64 / 4.0;
            let area = loading
                .initial_state(load, &gate)
                .unwrap()
                .contact()
                .unwrap()
                .contact_area_m2;
            difference *= (1.0 - 2.0 * (area / initial_area) * STEP).powi(4);
        }
        close(
            MASS * thermal.cells()[&CellId(0)].specific_h_j_kg,
            0.5 * (capacity + difference),
        );
        close(
            MASS * thermal.cells()[&CellId(1)].specific_h_j_kg,
            0.5 * (capacity - difference),
        );
        close(thermal.links()[&LINK].area_m2, initial_area);
        close(
            thermal.links()[&LINK].conductance_w_k / thermal.links()[&LINK].area_m2,
            capacity / initial_area,
        );
        assert_eq!(footprint.area_m2(), Some(thermal.links()[&LINK].area_m2));
        assert_eq!(report.release_time_s, None);
        assert_eq!(report.footprint_intervals, 4);
        assert_eq!(report.zero_crossing_time_s, (endpoint < 0.0).then_some(0.5));
        assert!(
            report.maximum_admitted_area_m2
                >= if endpoint < 0.0 {
                    zero_area
                } else {
                    initial_area
                }
        );
        assert_eq!(heat.ledger().dissipated_work_j(), 0.0);
        assert_eq!(report.transport.external_heat_j, 0.0);
    }

    // Endpoints fit this face exactly, but their intervening zero-load patch
    // does not. Even one broad footprint interval must refuse that geometry.
    assert!(zero_area > initial_area);
    let thermal = timeline(
        &caloric,
        initial,
        capacity,
        initial_area,
        [0.0; 2],
        0.0,
        vec![],
    );
    let heat = binding(loading, initial, &caloric, &thermal, &gate);
    let footprint = footprint(&heat, &thermal);
    let ramp = loading
        .advance_linear_ramp(initial, 0.0, 1.0, -0.6, 1e-15, &gate)
        .unwrap();
    assert!(
        heat.propose_timed_ramp(
            &thermal,
            heat.stamp(),
            &ramp,
            &footprint,
            Policy::Stop,
            JkrTimedCombinedHeatControls {
                maximum_footprint_dt_s: 1.0,
                ..controls()
            },
            &gate,
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heat.state(), initial);
}

#[test]
fn explicit_stop_and_nonseparating_boundaries_return_the_actual_untraversed_tail() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    for (d, b, t, policy, kind) in [
        (
            -1.0,
            4.0,
            8.0,
            Policy::Stop,
            JkrCombinedBoundary::RadialFold,
        ),
        (
            1.0,
            8.0,
            2.0,
            Policy::ReleaseIfSeparatingFold,
            JkrCombinedBoundary::OpeningEdge,
        ),
        (
            1.0,
            2.0,
            8.0,
            Policy::ReleaseIfSeparatingFold,
            JkrCombinedBoundary::RadialFold,
        ),
    ] {
        let (loading, _) = driver(&source, d, b, t, &gate);
        let initial = loading.initial_state(0.0, &gate).unwrap();
        let limit = loading.limit().unwrap();
        assert_eq!(limit.kind, kind);
        let capacity = initial.stored_energy_j().abs();
        let caloric = caloric(&source, capacity);
        let thermal = timeline(
            &caloric,
            initial,
            capacity,
            2.0 * initial.contact().unwrap().contact_area_m2,
            [0.0; 2],
            0.0,
            vec![],
        );
        let heat = binding(loading, initial, &caloric, &thermal, &gate);
        let footprint = footprint(&heat, &thermal);
        let ramp = loading
            .advance_linear_ramp(initial, 0.0, 1.0, 2.0 * limit.load_factor, 1e-15, &gate)
            .unwrap();
        let candidate = heat
            .propose_timed_ramp(
                &thermal,
                heat.stamp(),
                &ramp,
                &footprint,
                policy,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        let report = candidate.report();
        assert_eq!(report.final_time_s, 0.5);
        assert_eq!(report.untraversed_time_s, 0.5);
        assert_eq!(report.boundary_time_s, Some(0.5));
        assert_eq!(report.release_time_s, None);
        assert_eq!(candidate.thermal().time_s(), 0.5);
        assert_eq!(
            candidate.history().state().contact().unwrap().boundary,
            Some(kind)
        );
        assert!(candidate.thermal().links().contains_key(&LINK));
        assert!(candidate.footprint().area_m2().is_some());
        assert_eq!(candidate.history().ledger().dissipated_work_j(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(heat.state(), initial);
    }
}

#[test]
fn source_clock_and_shared_transport_budgets_refuse_without_partial_publication() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let capacity = initial.stored_energy_j().abs();
    let caloric = caloric(&source, capacity);
    let mut thermal = timeline(
        &caloric,
        initial,
        capacity,
        2.0 * initial.contact().unwrap().contact_area_m2,
        [0.0; 2],
        0.0,
        vec![],
    );
    let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
    let mut footprint = footprint(&heat, &thermal);
    let stamp = heat.stamp();
    let cells = thermal.cells().clone();
    let links = thermal.links().clone();
    let area = footprint.area_m2();
    let normal = loading
        .advance_linear_ramp(initial, 0.0, 1.0, 0.5, 1e-15, &gate)
        .unwrap();
    let wrong_clock = loading
        .advance_linear_ramp(initial, 0.25, 1.0, 0.5, 1e-15, &gate)
        .unwrap();
    let (foreign, _) = driver(&source, -1.0, 4.0, 9.0, &gate);
    let wrong_source = foreign
        .advance_linear_ramp(
            foreign.initial_state(0.0, &gate).unwrap(),
            0.0,
            1.0,
            0.5,
            1e-15,
            &gate,
        )
        .unwrap();
    for ramp in [&wrong_clock, &wrong_source] {
        assert!(
            heat.apply_timed_ramp(
                &mut thermal,
                stamp,
                ramp,
                &mut footprint,
                Policy::Stop,
                controls(),
                &gate,
                || Ok(())
            )
            .is_err()
        );
    }
    for limited in [
        JkrTimedCombinedHeatControls {
            maximum_transport_substeps: 15,
            ..controls()
        },
        JkrTimedCombinedHeatControls {
            maximum_footprint_intervals: 1,
            ..controls()
        },
    ] {
        assert!(
            heat.apply_timed_ramp(
                &mut thermal,
                stamp,
                &normal,
                &mut footprint,
                Policy::Stop,
                limited,
                &gate,
                || Ok(())
            )
            .is_err()
        );
    }
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.links(), &links);
    assert_eq!(footprint.area_m2(), area);
    let accepted = heat
        .apply_timed_ramp(
            &mut thermal,
            stamp,
            &normal,
            &mut footprint,
            Policy::Stop,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(accepted.transport.substeps, 16);
    assert_eq!(thermal.time_s(), 1.0);
}

#[test]
fn scheduled_link_replacement_recontact_and_final_cancellation_keep_all_three_owners() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let limit = loading.limit().unwrap();
    let capacity = loading
        .advance_to_load_factor(initial, 2.0 * limit.load_factor, &gate)
        .unwrap()
        .state()
        .stored_energy_j();
    let caloric = caloric(&source, capacity);
    let area = initial.contact().unwrap().contact_area_m2;
    let replacement = ThermalLink {
        id: LINK,
        a: CellId(0),
        b: CellId(1),
        area_m2: area,
        conductance_w_k: 2.0 * capacity,
    };
    let ramp = loading
        .advance_linear_ramp(initial, 0.0, 1.0, 2.0 * limit.load_factor, 1e-15, &gate)
        .unwrap();
    for later in [
        vec![
            TimedEvent {
                time_s: 0.25,
                event: Event::Separation(LINK),
            },
            TimedEvent {
                time_s: 0.25,
                event: Event::Contact(replacement),
            },
        ],
        vec![TimedEvent {
            time_s: 0.75,
            event: Event::Contact(replacement),
        }],
    ] {
        let mut thermal = timeline(
            &caloric,
            initial,
            capacity,
            2.0 * area,
            [0.0; 2],
            0.0,
            later,
        );
        let cells = thermal.cells().clone();
        let links = thermal.links().clone();
        let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
        let stamp = heat.stamp();
        let mut footprint = footprint(&heat, &thermal);
        assert!(
            heat.apply_timed_ramp(
                &mut thermal,
                stamp,
                &ramp,
                &mut footprint,
                Policy::ReleaseIfSeparatingFold,
                controls(),
                &gate,
                || Ok(())
            )
            .is_err()
        );
        assert_eq!(heat.stamp(), stamp);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.links(), &links);
        assert_eq!(footprint.area_m2(), Some(area));
    }

    let unrelated = ThermalLink {
        id: LINK + 1,
        a: CellId(0),
        b: CellId(1),
        area_m2: area / 16.0,
        conductance_w_k: 0.25 * capacity,
    };
    let original = timeline(
        &caloric,
        initial,
        capacity,
        2.0 * area,
        [0.0; 2],
        0.0,
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(unrelated),
        }],
    );
    let mut reference_thermal = original.clone();
    let mut reference = binding(loading, initial, &caloric, &original, &gate);
    let mut reference_footprint = footprint(&reference, &original);
    let mut polls = 0;
    let expected = reference
        .apply_timed_ramp(
            &mut reference_thermal,
            reference.stamp(),
            &ramp,
            &mut reference_footprint,
            Policy::ReleaseIfSeparatingFold,
            controls(),
            &gate,
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut thermal = original.clone();
    let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
    let mut footprint = footprint(&heat, &thermal);
    let stamp = heat.stamp();
    let mut calls = 0;
    assert!(
        heat.apply_timed_ramp(
            &mut thermal,
            stamp,
            &ramp,
            &mut footprint,
            Policy::ReleaseIfSeparatingFold,
            controls(),
            &gate,
            || {
                calls += 1;
                if calls == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        )
        .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), original.cells());
    assert_eq!(thermal.links(), original.links());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(footprint.area_m2(), Some(area));
    let retry = heat
        .apply_timed_ramp(
            &mut thermal,
            stamp,
            &ramp,
            &mut footprint,
            Policy::ReleaseIfSeparatingFold,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(thermal.cells(), reference_thermal.cells());
    assert_eq!(thermal.links(), reference_thermal.links());
    assert_eq!(thermal.links().len(), 1);
    assert_eq!(thermal.links()[&unrelated.id], unrelated);
    assert_eq!(footprint.area_m2(), None);
}
