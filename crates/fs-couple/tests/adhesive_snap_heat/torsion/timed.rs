//! G1/G3/G4: original twist loss and snap release during finite heat transport.

use super::{
    CRITICAL_TWIST, MASS, binding, close, controls, events, ledgers, loaded, source, thermal,
};
use fs_contact::normal_patch::adhesive::{JkrBranch, SourceJkrSphereContact};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    torsion::{JkrTorsionHeatError, timed::JkrTimedTorsionHeatControls},
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_thermal::enthalpy::{
    Boundary, CellId, DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId,
    Reservoir, StepConfig, ThermalLink, TimedEvent,
};

const STEP: f64 = 1.0 / 16.0;

fn timed_controls() -> JkrTimedTorsionHeatControls {
    JkrTimedTorsionHeatControls {
        heat: controls(),
        transport: StepConfig {
            max_dt_s: STEP,
            energy_tolerance_j: 2e-22,
        },
        maximum_transport_substeps: 64,
        total_energy_tolerance_j: 8e-22,
    }
}

// Original material identities with a manufactured constant finite capacity.
fn linear_curves(
    normal: &SourceJkrSphereContact,
    capacity_j_k: f64,
) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let elastic = normal.normal_state().elastic();
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

fn boundary(conductance_w_k: f64) -> Boundary {
    Boundary {
        ambient_k: 290.0,
        htc_w_m2_k: 0.0,
        reservoir: Some(Reservoir {
            area_m2: 1e-6,
            conductance_w_k,
            temperature_k: 290.0,
        }),
        source_w: 0.0,
    }
}

// Reuse the original finite cells, masses, material laws and geometry, then
// declare the independent conduction/cooling experiment on the native timeline.
fn timeline(
    caloric: &[EquilibriumEnthalpyPhaseCurve; 2],
    reservoir_w_k: f64,
    link_w_k: f64,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'_> {
    let initial = thermal(caloric, 1.0);
    let mut schedule: Vec<_> = initial
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.boundary = boundary(reservoir_w_k);
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 0,
            a: CellId(0),
            b: CellId(1),
            area_m2: 1e-6,
            conductance_w_k: link_w_k,
        }),
    });
    schedule.extend(later);
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &caloric[0]), (MaterialId(1), &caloric[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, timed_controls().transport, || Ok(()))
    .unwrap()
    .0
}

#[test]
fn plastic_twist_power_follows_independent_conduction_and_cooling_modes() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let delta = initial.normal_response().indentation_m;
    let prior_loss = initial.accumulated_torsion_loss_j();
    let angle = 2.0 * CRITICAL_TWIST;
    let original = events.advance_torsion(&initial, angle).unwrap();
    let loss = original.report().irreversible_torsion_loss_j;
    assert!(loss > 0.0);
    let caloric = linear_curves(&normal, loss);
    // C = loss J/K, reservoir G = C/s and link G = C/s. Sum and
    // difference modes decay at 1/s and 3/s under the original explicit steps.
    let mut thermal = timeline(&caloric, loss, loss, vec![]);
    let mut heat = binding(events, initial, &caloric, &thermal);
    let report = heat
        .apply_interval(
            &mut thermal,
            heat.stamp(),
            delta,
            angle,
            1.0,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(thermal.time_s(), 1.0);
    assert_eq!(report.transport.substeps, 16);
    assert_eq!(heat.state(), original.state());
    assert_eq!(report.normal.work.dissipated_work_j(), 0.0);
    close(heat.state().accumulated_torsion_loss_j(), prior_loss + loss);
    close(heat.torsion_loss_ledger().dissipated_work_j(), loss);
    close(report.injected_surface_heat_j[0], 0.5 * loss);
    close(report.injected_surface_heat_j[1], 0.25 * loss);
    let sum = 0.75 * loss * (1.0 - (1.0 - STEP).powi(16));
    let difference = 0.25 * loss * (1.0 - (1.0 - 3.0 * STEP).powi(16)) / 3.0;
    for (side, expected) in [(0, (sum + difference) / 2.0), (1, (sum - difference) / 2.0)] {
        close(
            MASS * thermal.cells()[&CellId(side)].specific_h_j_kg,
            expected,
        );
    }
    close(report.external_transport_heat_j, sum - 0.75 * loss);
    assert!(report.external_transport_heat_j < 0.0);
    assert!(sum < 0.7 * 0.75 * loss);
    assert!(report.absolute_stage_energy_residual_j <= timed_controls().total_energy_tolerance_j);
    let next = heat
        .apply_interval(
            &mut thermal,
            heat.stamp(),
            delta,
            0.0,
            0.5,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.initial_time_s, 1.0);
    assert_eq!(thermal.time_s(), 1.5);
    assert_eq!(next.torsion_loss.dissipated_work_j(), 0.0);
    close(heat.torsion_loss_ledger().dissipated_work_j(), loss);
    close(
        MASS * thermal.cells()[&CellId(0)].specific_h_j_kg
            + MASS * thermal.cells()[&CellId(1)].specific_h_j_kg,
        sum * (1.0 - STEP).powi(8),
    );
}

#[test]
fn fold_spring_and_normal_release_precede_the_scheduled_cooling_change() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let prior_loss = initial.accumulated_torsion_loss_j();
    let fold = normal.minimum_attached_indentation_m();
    let original = events.advance_normal(&initial, fold).unwrap();
    let normal_release = original.report().normal_released_energy_j;
    // Circular-punch torsional stiffness scales with a^3, and the original
    // displacement fold has a_fold^3 = a_zero_normal_force^3 / 36.
    let spring_release = initial.torsion_storage_j() / 36.0;
    close(original.report().torsion_released_energy_j, spring_release);
    let release = normal_release + spring_release;
    let caloric = linear_curves(&normal, release);
    let changes = (0..2)
        .map(|side| TimedEvent {
            time_s: 0.5,
            event: Event::Boundary(CellId(side), boundary(2.0 * release)),
        })
        .collect();
    let mut thermal = timeline(&caloric, 0.0, 0.0, changes);
    let mut heat = binding(events, initial, &caloric, &thermal);
    let report = heat
        .apply_interval(
            &mut thermal,
            heat.stamp(),
            fold,
            0.4,
            1.0,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.normal.segments[0], *original.report());
    close(
        report.normal.normal_release.dissipated_work_j(),
        normal_release,
    );
    close(
        report.normal.torsion_release.dissipated_work_j(),
        spring_release,
    );
    close(report.normal.thermal.stored_energy_change_j, 0.75 * release);
    assert_eq!(report.normal.thermal.substeps, 0);
    assert_eq!(report.transport.events, 2);
    assert_eq!(report.torsion_loss.dissipated_work_j(), 0.0);
    assert_eq!(report.injected_surface_heat_j, [0.0; 2]);
    assert_eq!(heat.state().normal_response().branch, JkrBranch::Detached);
    assert_eq!(heat.state().torsion_storage_j(), 0.0);
    assert_eq!(heat.state().accumulated_torsion_loss_j(), prior_loss);
    for (side, fraction) in [(0, 0.25), (1, 0.5)] {
        close(
            MASS * thermal.cells()[&CellId(side)].specific_h_j_kg,
            fraction * release * (1.0 - 2.0 * STEP).powi(8),
        );
    }
    close(
        report.thermal.external_heat_j,
        report.thermal.stored_energy_change_j,
    );
    assert!(report.external_transport_heat_j < 0.0);
    assert_eq!(thermal.time_s(), 1.0);
}

#[test]
fn late_refusals_preserve_all_owners_and_retry_accepts_the_original_proposal() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let fold = normal.minimum_attached_indentation_m();
    let release = events
        .advance_normal(&initial, fold)
        .unwrap()
        .report()
        .released_energy_j;
    let caloric = linear_curves(&normal, release);
    let mut thermal = timeline(&caloric, release, 0.0, vec![]);
    let mut heat = binding(events.clone(), initial.clone(), &caloric, &thermal);
    let stamp = heat.stamp();
    let cells = thermal.cells().clone();
    let links = thermal.links().clone();
    let original_ledgers = ledgers(&heat);
    let mut bounded = timed_controls();
    bounded.maximum_transport_substeps = 1;
    assert!(matches!(
        heat.apply_interval(&mut thermal, stamp.clone(), fold, 0.0, 1.0, bounded, || Ok(
            ()
        ),),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Budget(_)
        )))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.links(), &links);
    assert_eq!(thermal.time_s(), 0.0);
    let mut proposal_calls = 0;
    let proposal = heat
        .propose_interval(
            &thermal,
            stamp.clone(),
            fold,
            0.0,
            1.0,
            timed_controls(),
            || {
                proposal_calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        heat.apply_interval(
            &mut thermal,
            stamp.clone(),
            fold,
            0.0,
            1.0,
            timed_controls(),
            || {
                calls += 1;
                if calls == proposal_calls + 1 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(calls, proposal_calls + 1);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.links(), &links);
    assert_eq!(thermal.time_s(), 0.0);
    let mut advanced = thermal
        .advance_to(STEP, timed_controls().transport, || Ok(()))
        .unwrap()
        .0;
    assert!(matches!(
        heat.apply_interval(
            &mut advanced,
            stamp.clone(),
            fold,
            0.0,
            1.0,
            timed_controls(),
            || Ok(()),
        ),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::HistoryMismatch))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    assert_eq!(advanced.cells(), &cells);
    assert_eq!(advanced.time_s(), STEP);
    let accepted = heat
        .apply_interval(
            &mut thermal,
            stamp.clone(),
            fold,
            0.0,
            1.0,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(&accepted, proposal.report());
    assert_eq!(heat.stamp(), proposal.history().stamp());
    assert_eq!(ledgers(&heat), ledgers(proposal.history()));
    assert_eq!(thermal.cells(), proposal.thermal().cells());
    assert_eq!(thermal.links(), proposal.thermal().links());
    assert_eq!(thermal.time_s(), proposal.thermal().time_s());
    assert!(matches!(
        heat.apply_interval(
            &mut thermal,
            stamp,
            fold,
            0.0,
            1.0,
            timed_controls(),
            || Ok(()),
        ),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::HistoryMismatch))
    ));
    assert_eq!(heat.stamp(), proposal.history().stamp());
    assert_eq!(thermal.cells(), proposal.thermal().cells());

    // A separate narrow material domain refuses after staging plastic motion
    // and before any proposed caloric transport or source history can publish.
    let delta = initial.normal_response().indentation_m;
    let angle = 2.0 * CRITICAL_TWIST;
    let loss = events
        .advance_torsion(&initial, angle)
        .unwrap()
        .report()
        .irreversible_torsion_loss_j;
    let narrow = linear_curves(&normal, loss * 1e-6);
    let mut refused_thermal = timeline(&narrow, 0.0, 0.0, vec![]);
    let refused_cells = refused_thermal.cells().clone();
    let mut refused = binding(events, initial, &narrow, &refused_thermal);
    let refused_stamp = refused.stamp();
    let refused_ledgers = ledgers(&refused);
    assert!(matches!(
        refused.apply_interval(
            &mut refused_thermal,
            refused_stamp.clone(),
            delta,
            angle,
            1.0,
            timed_controls(),
            || Ok(()),
        ),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(refused.stamp(), refused_stamp);
    assert_eq!(ledgers(&refused), refused_ledgers);
    assert_eq!(refused_thermal.cells(), &refused_cells);
    assert_eq!(refused_thermal.time_s(), 0.0);
}
