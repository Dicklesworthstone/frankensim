//! Original normal snap work, material calorics and native heat impulses.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrEnvelope, JkrJumpKind, JkrState, SourceJkrSphereContact,
};
use fs_couple::friction_heat::adhesive::snap::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_thermal::enthalpy::*;
use fs_tribo::{DissipationStep, HeatPartition};

const MASS: f64 = 1e-6;

fn source(adhesion: f64) -> SourceJkrSphereContact {
    let a = common::solid("snap-surface-a", 1000.0, 1e7);
    let b = common::solid("snap-surface-b", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, adhesion).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}

fn curves(source: &SourceJkrSphereContact, span: f64) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let e = source.normal_state().elastic();
    let cards = [
        e.surface_a_state().resolved().card_identity(),
        e.surface_b_state().resolved().card_identity(),
    ];
    // An explicitly manufactured latent plateau; no fitted material claim.
    let h = -source.evaluate_attached(0.0).unwrap().reversible_energy_j / MASS;
    cards.map(|card| {
        EquilibriumEnthalpyPhaseCurve::try_new(
            card,
            [
                (0.0, 290.0, 0.0),
                (0.1, 300.0, 0.0),
                (0.9, 300.0, 1.0),
                (10.0, 400.0, 1.0),
            ]
            .map(
                |(v, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: v * h * span,
                    temperature_k,
                    liquid_mass_fraction,
                    bulk_density_kg_m3: 1000.0,
                },
            )
            .to_vec(),
        )
        .unwrap()
    })
}

fn thermal(
    curves: &[EquilibriumEnthalpyPhaseCurve; 2],
    mass_scale: f64,
) -> DepositionSimulation<'_> {
    let schedule = (0..2)
        .map(|side| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(side),
                segment_id: side,
                body_id: side,
                material_id: MaterialId(side),
                mass_kg: MASS * mass_scale,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [side as f64, 0.0, 0.0],
                    end_m: [side as f64 + 1e-3, 0.0, 0.0],
                    volume_m3: 1e-9,
                    surface_area_m2: 6e-6,
                },
                boundary: Boundary {
                    ambient_k: 290.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &curves[0]), (MaterialId(1), &curves[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 10,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-20,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn controls() -> JkrSnapHeatControls {
    JkrSnapHeatControls {
        maximum_segments: 16,
        cell_energy_tolerance_j: 1e-22,
        thermal_energy_tolerance_j: 2e-22,
        transfer_energy_tolerance_j: 4e-22,
    }
}

fn end_state(source: &SourceJkrSphereContact, mut state: JkrState, endpoints: &[f64]) -> JkrState {
    for &end in endpoints {
        state = source
            .advance_monotone(state, end)
            .unwrap()
            .response
            .state();
    }
    state
}

fn binding<'s>(
    source: &'s SourceJkrSphereContact,
    curves: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
    state: JkrState,
) -> JkrSnapHeat<'s> {
    JkrSnapHeat::new(
        source,
        state,
        thermal,
        [CellId(0), CellId(1)],
        [&curves[0], &curves[1]],
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("resolved elastic waves".into()),
    )
    .unwrap()
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn attachment_heats_original_latent_storage_once_with_explicit_outside_work() {
    let source = source(0.05);
    let caloric = curves(&source, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let before = thermal.cells().clone();
    let initial = source.detached_state(-1e-6).unwrap();
    let mut heat = binding(&source, &caloric, &thermal, initial);
    let expected = heat.stamp();
    let final_state = end_state(&source, initial, &[0.0]);
    let report = heat
        .apply(
            &mut thermal,
            expected,
            &[0.0],
            final_state,
            controls(),
            || Ok(()),
        )
        .unwrap();
    let release = -source.evaluate_attached(0.0).unwrap().reversible_energy_j;
    assert_eq!(report.jumps.len(), 1);
    assert_eq!(report.jumps[0].kind, JkrJumpKind::Attach);
    assert_eq!(report.jumps[0].indentation_m, 0.0);
    assert_eq!(heat.state().branch(), JkrBranch::Attached);
    close(report.work.dissipated_work_j(), release);
    close(report.work.surface_a_heat_j(), 0.25 * release);
    close(report.work.surface_b_heat_j(), 0.5 * release);
    close(report.work.other_work_j(), 0.25 * release);
    close(report.thermal.stored_energy_change_j, 0.75 * release);
    assert!(report.transfer_residual_j.abs() <= controls().transfer_energy_tolerance_j);
    assert_eq!(thermal.time_s(), 0.0);
    for id in [CellId(0), CellId(1)] {
        close(thermal.temperature_k(id).unwrap(), 300.0);
        assert_eq!(thermal.cells()[&id].mass_kg, before[&id].mass_kg);
        assert_eq!(thermal.cells()[&id].geometry, before[&id].geometry);
        assert!(thermal.cells()[&id].specific_h_j_kg > 0.0);
    }
    let accepted = thermal.cells().clone();
    assert_eq!(
        heat.apply(
            &mut thermal,
            expected,
            &[0.0],
            final_state,
            controls(),
            || Ok(())
        ),
        Err(JkrSnapHeatError::HistoryMismatch)
    );
    assert_eq!(thermal.cells(), &accepted);
    let no_jump = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[0.0],
            final_state,
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(no_jump.work.dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &accepted);
}

#[test]
fn batched_adhesive_cycle_balances_original_work_and_detaches_at_exact_fold() {
    let source = source(0.05);
    let caloric = curves(&source, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let initial = source.detached_state(-1e-6).unwrap();
    let fold = source.minimum_attached_indentation_m();
    let endpoints = [0.0, 1e-6, 0.5 * fold, fold, -1e-6];
    let mut state = initial;
    let mut prescribed = 0.0;
    for &end in &endpoints {
        let transition = source.advance_monotone(state, end).unwrap();
        prescribed += transition.prescribed_work_j;
        state = transition.response.state();
    }
    let mut heat = binding(&source, &caloric, &thermal, initial);
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &endpoints,
            state,
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.jumps.len(), 2);
    assert_eq!(report.jumps[1].kind, JkrJumpKind::Detach);
    assert_eq!(report.jumps[1].indentation_m, fold);
    assert_eq!(report.final_state.branch(), JkrBranch::Detached);
    close(report.work.dissipated_work_j(), prescribed);
    close(
        report.thermal.stored_energy_change_j + report.work.other_work_j(),
        prescribed,
    );
    close(heat.ledger().dissipated_work_j(), prescribed);
    assert_eq!(report.thermal.substeps, 0);
    assert_eq!(report.thermal.events, 0);
}

#[test]
fn reversible_tensile_and_compressive_work_never_heats_the_cells() {
    let source = source(0.05);
    let caloric = curves(&source, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let before = thermal.cells().clone();
    let initial = source.evaluate_attached(0.0).unwrap().state();
    let endpoints = [0.5 * source.minimum_attached_indentation_m(), 1e-6, 0.0];
    let mut heat = binding(&source, &caloric, &thermal, initial);
    let final_state = end_state(&source, initial, &endpoints);
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &endpoints,
            final_state,
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert!(report.jumps.is_empty());
    assert_eq!(report.work.dissipated_work_j(), 0.0);
    assert_eq!(report.thermal.stored_energy_change_j, 0.0);
    assert_eq!(thermal.cells(), &before);
}

#[test]
fn original_source_path_material_objects_and_finite_receivers_remain_bound() {
    let source = source(0.05);
    let other_source = self::source(0.06);
    let caloric = curves(&source, 1.0);
    let thermal = thermal(&caloric, 1.0);
    let initial = source.detached_state(-1e-6).unwrap();
    let heat = binding(&source, &caloric, &thermal, initial);
    let foreign = other_source.evaluate_attached(0.0).unwrap().state();
    assert!(matches!(
        heat.propose(&thermal, heat.stamp(), &[0.0], foreign, controls(), || Ok(
            ()
        )),
        Err(JkrSnapHeatError::HistoryMismatch)
    ));
    assert!(matches!(
        JkrSnapHeat::new(
            &source,
            foreign,
            &thermal,
            [CellId(0), CellId(1)],
            [&caloric[0], &caloric[1]],
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            None
        ),
        Err(JkrSnapHeatError::HistoryMismatch)
    ));
    let identical_copy = caloric.clone();
    assert!(
        JkrSnapHeat::new(
            &source,
            initial,
            &thermal,
            [CellId(0), CellId(1)],
            [&identical_copy[0], &identical_copy[1]],
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            None
        )
        .is_err()
    );
    assert!(
        JkrSnapHeat::new(
            &source,
            initial,
            &thermal,
            [CellId(0), CellId(1)],
            [&caloric[1], &caloric[0]],
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            None
        )
        .is_err()
    );
    assert!(
        JkrSnapHeat::new(
            &source,
            initial,
            &thermal,
            [CellId(0), CellId(1)],
            [&caloric[0], &caloric[1]],
            HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
            None
        )
        .is_err()
    );
    let changed_mass = self::thermal(&caloric, 2.0);
    let final_state = end_state(&source, initial, &[0.0]);
    assert!(matches!(
        heat.propose(
            &changed_mass,
            heat.stamp(),
            &[0.0],
            final_state,
            controls(),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Invalid(_))
    ));
}

#[test]
fn cancellation_at_every_phase_and_caloric_refusal_leave_both_owners_retryable() {
    let source = source(0.05);
    let caloric = curves(&source, 1.0);
    let original = thermal(&caloric, 1.0);
    let initial = source.detached_state(-1e-6).unwrap();
    let endpoints = [0.0, source.minimum_attached_indentation_m()];
    let final_state = end_state(&source, initial, &endpoints);
    let reference = binding(&source, &caloric, &original, initial);
    let mut calls = 0;
    let proposal = reference
        .propose(
            &original,
            reference.stamp(),
            &endpoints,
            final_state,
            controls(),
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    for stop in 1..=calls {
        let mut thermal = original.clone();
        let mut heat = binding(&source, &caloric, &thermal, initial);
        let stamp = heat.stamp();
        let mut polls = 0;
        let result = heat.apply(
            &mut thermal,
            stamp,
            &endpoints,
            final_state,
            controls(),
            || {
                polls += 1;
                if polls == stop {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(
            result,
            Err(JkrSnapHeatError::Thermal(DepositionError::Cancelled))
        );
        assert_eq!(heat.stamp(), stamp);
        assert_eq!(heat.ledger().dissipated_work_j(), 0.0);
        assert_eq!(thermal.cells(), original.cells());
        let retry = heat
            .apply(
                &mut thermal,
                stamp,
                &endpoints,
                final_state,
                controls(),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(&retry, proposal.report());
        assert_eq!(thermal.cells(), proposal.thermal().cells());
    }
    // The first cell accepts its heat before the second cell refuses. Neither
    // private caloric change nor accepted branch history may escape that failure.
    let mut narrow = curves(&source, 1e-4);
    narrow[0] = caloric[0].clone();
    let mut thermal = thermal(&narrow, 1.0);
    let before = thermal.cells().clone();
    let mut heat = binding(&source, &narrow, &thermal, initial);
    let stamp = heat.stamp();
    assert_eq!(
        heat.apply(
            &mut thermal,
            stamp,
            &endpoints,
            final_state,
            controls(),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Thermal(DepositionError::MaterialDomain))
    );
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), &before);
}

#[test]
fn explicit_transfer_gate_refuses_an_overallocated_declared_partition() {
    let source = source(0.05);
    let caloric = curves(&source, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let before = thermal.cells().clone();
    let initial = source.detached_state(-1e-6).unwrap();
    // This is within the partition owner's dimensionless roundoff allowance,
    // but the caller's dimensional transfer budget still independently applies.
    let partition = HeatPartition::new(0.25, 0.5, 0.25 + 32.0 * f64::EPSILON).unwrap();
    let mut heat = JkrSnapHeat::new(
        &source,
        initial,
        &thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        partition,
        Some("resolved elastic waves".into()),
    )
    .unwrap();
    let stamp = heat.stamp();
    let final_state = end_state(&source, initial, &[0.0]);
    let strict = JkrSnapHeatControls {
        transfer_energy_tolerance_j: 0.0,
        ..controls()
    };
    assert!(matches!(
        heat.apply(&mut thermal, stamp, &[0.0], final_state, strict, || Ok(())),
        Err(JkrSnapHeatError::TransferBalance(_))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.ledger().dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &before);
}

#[test]
fn already_integrated_work_needs_no_artificial_power_or_duration() {
    let partition = HeatPartition::new(0.25, 0.5, 0.25).unwrap();
    let integrated = DissipationStep::from_work(7.0, partition).unwrap();
    assert_eq!(
        integrated,
        DissipationStep::from_power(14.0, 0.5, partition).unwrap()
    );
    assert_eq!(
        DissipationStep::from_work(0.0, partition)
            .unwrap()
            .total_work_j(),
        0.0
    );
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(DissipationStep::from_work(invalid, partition).is_err());
    }
}
