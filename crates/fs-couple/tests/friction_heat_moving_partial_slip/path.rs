//! G1/G3/G4: accelerated residence, cyclic shear and complete path rollback.
use super::*;
use fs_couple::friction_heat::moving::partial_slip::path::*;
use fs_couple::friction_heat::moving::trajectory::QuadraticContactPath;

fn limits(dt: f64) -> MovingPartialSlipPathControls {
    MovingPartialSlipPathControls {
        step: controls(),
        maximum_step_s: dt,
        maximum_steps: 128,
        maximum_trials: 20000,
        total_energy_tolerance_j: 2e-8,
    }
}
fn excursion() -> QuadraticContactPath {
    QuadraticContactPath {
        start: position(0.0, 0.0, 0.0),
        midpoint: position(0.04, 0.0, 0.04),
        end: position(0.0, 0.0, 0.0),
    }
}
fn run(
    contact: &MovingPartialSlip,
    thermal: &DepositionSimulation<'static>,
    target: f64,
    path: QuadraticContactPath,
    dt: f64,
    mu: f64,
) -> MovingPartialSlipPathCandidate<'static, ()> {
    contact
        .advance_quadratic_candidate(
            thermal,
            target,
            path,
            &work(),
            limits(dt),
            300.0,
            |_, t, r, _, _| Ok::<_, Infallible>(source(t, r, mu)),
            || false,
        )
        .unwrap()
}

#[test]
fn closed_sliding_excursion_retains_hysteresis_and_physical_heat() {
    let (thermal, contact) = fixture(&[-1.0, 2.0], 1.0, [0.0; 2]);
    let candidate = run(&contact, &thermal, 1.0, excursion(), 1.0, 1.0);
    let r = candidate.report();
    assert_eq!(
        r.steps.len(),
        2,
        "relative motion reversal is a required event"
    );
    assert_eq!(r.steps[0].phase.end_s, 0.5);
    near(r.steps[0].phase.relative_velocity_m_s[0], 0.16, 1e-15);
    near(r.steps[0].phase.relative_velocity_m_s[1], 0.0, 1e-15);
    near(r.steps[1].phase.relative_velocity_m_s[1], -0.16, 1e-15);
    assert_eq!(candidate.contact().position(), contact.position());
    near(r.relative_travel_m, 0.08, 1e-15);
    // Independent two-spring cycle: yield strains C/k; forward plastic travel
    // D-C/k, reverse plastic travel D-2C/k for D=.04 on both branches.
    let expected_heat: f64 = [(0.3, 60.0), (0.7, 40.0)]
        .iter()
        .map(|&(c, k)| c * (0.08 - 3.0 * c / k))
        .sum();
    near(expected_heat, 0.03875, 1e-15);
    near(r.physical_heat_j, expected_heat, 3e-14);
    near(r.storage_after_j, 0.006875, 2e-14);
    near(r.numerical_dissipation_j, 0.034375, 3e-14);
    near(r.mechanical_work_j, 0.08, 3e-14);
    near(r.thermal_energy_change_j, expected_heat, 2e-9);
    near(r.elastic_outflow_j, 0.0, 1e-14);
    near(r.total_energy_residual_j, 0.0, 2e-9);
    for parcel in candidate.contact().parcels() {
        let z = parcel.state().region_elastic_displacement_m();
        near(z[0][0], -0.005, 1e-14);
        near(z[1][0], -0.0175, 1e-14);
    }
    // Subintervals retain distinct work ownership, and no caller state moved.
    let ids: Vec<_> = r
        .steps
        .iter()
        .map(|s| {
            s.report.regions[0]
                .response
                .generalized_work
                .ownership
                .interval_id()
        })
        .collect();
    assert_eq!(ids, ["interval/segment/0", "interval/segment/1"]);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.parcels()[0].state(), &PartialSlipState::zero());
}

#[test]
fn accelerated_sliding_deposits_heat_by_quadratic_residence_not_endpoint_chord() {
    let (thermal, contact) = fixture(&[-1.0, 0.0, 1.0, 2.0], 1.0, [0.005, 0.0175]);
    let path = QuadraticContactPath {
        start: contact.position(),
        midpoint: position(0.025, 0.0, 0.025),
        end: position(0.1, 0.0, 0.1),
    };
    let candidate = run(&contact, &thermal, 1.0, path, 1.0, 1.0);
    assert_eq!(candidate.report().steps.len(), 1);
    let r = &candidate.report().steps[0].report;
    near(r.heat.dissipated_j, 0.1, 2e-14);
    // Footprint [d,1+d], d=.1t². Its residence in B's [1,2] cell
    // integrates to .1/3 m², compared with .1/2 for the endpoint chord.
    let right: Vec<_> = r
        .heat
        .patches
        .iter()
        .filter(|p| p.cells[1] == CellId(5))
        .collect();
    near(right.iter().map(|p| p.mean_area_m2).sum(), 0.1 / 3.0, 2e-15);
    near(
        right.iter().map(|p| p.dissipated_j).sum(),
        0.01 / 3.0,
        2e-14,
    );
    near(
        candidate.report().physical_heat_j,
        candidate.report().thermal_energy_change_j,
        2e-9,
    );
    // The native step holds each regional heat increment at its mean power.
    // Refine that time split against the continuous power/residence integral:
    // integral_0^1 (.2t)*(.1t²) dt = .005 J in the right-hand B cell.
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let refined = run(&contact, &thermal, 1.0, path, 1.0 / n as f64, 1.0);
        let work: f64 = refined
            .report()
            .steps
            .iter()
            .flat_map(|s| &s.report.heat.patches)
            .filter(|p| p.cells[1] == CellId(5))
            .map(|p| p.dissipated_j)
            .sum();
        errors.push((work - 0.005).abs());
        near(refined.report().physical_heat_j, 0.1, 3e-13);
    }
    assert!(errors.windows(2).all(|w| w[1] < 0.3 * w[0]), "{errors:?}");
}

#[test]
fn returning_footprint_does_not_resurrect_departed_shear_memory() {
    let (thermal, contact) = fixture(&[-1.0, 2.0], 1000.0, [0.01, 0.02]);
    let path = QuadraticContactPath {
        start: contact.position(),
        midpoint: position(0.0, 0.0, 0.25),
        end: contact.position(),
    };
    let candidate = run(&contact, &thermal, 1.0, path, 1.0, 1000.0);
    let r = candidate.report();
    assert_eq!(r.steps.len(), 2);
    near(r.relative_travel_m, 0.0, 0.0);
    near(r.storage_before_j, 0.011, 1e-14);
    near(r.elastic_outflow_j, 0.00275, 1e-14);
    near(r.storage_after_j, 0.00825, 1e-14);
    near(r.physical_heat_j, 0.0, 1e-14);
    assert_eq!(candidate.contact().parcels()[0].span_m, [0.0, 0.25]);
    assert_eq!(
        candidate.contact().parcels()[0]
            .state()
            .region_elastic_displacement_m(),
        [[0.0; 2]; 2]
    );
    near(r.total_energy_residual_j, 0.0, 1e-13);
}

#[test]
fn quadratic_shared_motion_preserves_memory_and_actual_mean_shear_moment() {
    let state = PartialSlipState::from_region_checkpoint(
        [[0.0, 0.01]; 2],
        [0.0; 2],
        [false; 2],
        Some(frame()),
        4,
    )
    .unwrap();
    let (thermal, contact) = fixture_state(&[-3.0, 5.0], 1000.0, state);
    let candidate = run(
        &contact,
        &thermal,
        1.0,
        QuadraticContactPath {
            start: contact.position(),
            midpoint: position(0.075, 0.075, 0.075),
            end: position(0.3, 0.3, 0.3),
        },
        1.0,
        1000.0,
    );
    let r = &candidate.report().steps[0].report;
    near(r.force_on_a_n[2], 1.0, 2e-14);
    // Integral(.5+.3t²)dt = .6, with force +z and moment -y.
    near(r.moment_on_a_about_origin_nm[1], -0.6, 2e-14);
    near(r.storage_after_j, 0.005, 2e-14);
    near(r.heat.dissipated_j, 0.0, 1e-14);
    near(r.elastic_outflow_j, 0.0, 1e-14);
}

#[test]
fn cyclic_refinement_removes_numerical_loss_without_heating_it() {
    let (thermal, contact) = fixture(&[-1.0, 2.0], 1.0, [0.0; 2]);
    let mut defects = Vec::new();
    for n in [8, 16, 32] {
        let candidate = run(&contact, &thermal, 1.0, excursion(), 1.0 / n as f64, 1.0);
        let r = candidate.report();
        near(r.relative_travel_m, 0.08, 2e-15);
        near(r.physical_heat_j, 0.03875, 2e-13);
        near(r.thermal_energy_change_j, r.physical_heat_j, 2e-8);
        near(r.storage_after_j, 0.006875, 2e-13);
        near(
            r.mechanical_work_j - r.physical_heat_j - r.storage_after_j,
            r.numerical_dissipation_j,
            2e-13,
        );
        defects.push(r.numerical_dissipation_j);
    }
    assert!(defects.windows(2).all(|w| w[1] < 0.7 * w[0]), "{defects:?}");
}

#[test]
fn late_refusal_shared_budgets_cancellation_and_hidden_excursion_are_transactional() {
    let (thermal, contact) = fixture(&[-1.0, 2.0], 1.0, [0.0; 2]);
    let before = contact.parcels().to_vec();
    let mut limit = limits(1.0);
    limit.maximum_steps = 1;
    let no_source = |_: &MovingPartialSlipPathPhase,
                     _: &DepositionSimulation<'static>,
                     _: &[MovingSlipRegion],
                     _: &[SweptContactPatch],
                     _: &[f64]|
     -> Result<MovingPartialSlipSource<()>, &'static str> {
        panic!("preflight source call")
    };
    assert!(matches!(
        contact.advance_quadratic_candidate(
            &thermal,
            1.0,
            excursion(),
            &work(),
            limit,
            300.0,
            no_source,
            || false
        ),
        Err(MovingPartialSlipError::StepBudget)
    ));
    let late = contact.advance_quadratic_candidate(
        &thermal,
        1.0,
        excursion(),
        &work(),
        limits(1.0),
        300.0,
        |phase, t, r, _, _| {
            if phase.index == 1 {
                Err("late material refusal")
            } else {
                Ok(source(t, r, 1.0))
            }
        },
        || false,
    );
    assert!(matches!(
        late,
        Err(MovingPartialSlipError::Source("late material refusal"))
    ));
    let seen_second = Cell::new(false);
    let cancelled = contact.advance_quadratic_candidate(
        &thermal,
        1.0,
        excursion(),
        &work(),
        limits(1.0),
        300.0,
        |phase, t, r, _, _| {
            if phase.index == 1 {
                seen_second.set(true);
            }
            Ok::<_, Infallible>(source(t, r, 1.0))
        },
        || seen_second.get(),
    );
    assert!(matches!(
        cancelled,
        Err(MovingPartialSlipError::Thermal(DepositionError::Cancelled))
    ));
    let complete = run(&contact, &thermal, 1.0, excursion(), 1.0, 1.0);
    limit = limits(1.0);
    limit.maximum_trials = complete.report().steps[0].report.trials + 1;
    assert!(matches!(
        contact.advance_quadratic_candidate(
            &thermal,
            1.0,
            excursion(),
            &work(),
            limit,
            300.0,
            |_, t, r, _, _| Ok::<_, Infallible>(source(t, r, 1.0)),
            || false
        ),
        Err(MovingPartialSlipError::TrialBudget)
    ));
    assert_eq!(contact.parcels(), before);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    let retry = run(&contact, &thermal, 1.0, excursion(), 1.0, 1.0);
    assert_eq!(retry.contact().parcels(), complete.contact().parcels());
    assert_eq!(retry.thermal().cells(), complete.thermal().cells());
    let (short_thermal, short) = fixture(&[0.0, 1.24], 1.0, [0.0; 2]);
    // All three supplied placements fit; the actual interior vertex at t=.625
    // reaches .25 and leaves the material surface. Endpoint checks cannot see it.
    assert!(matches!(
        short.advance_quadratic_candidate(
            &short_thermal,
            1.0,
            QuadraticContactPath {
                start: short.position(),
                midpoint: position(0.0, 0.0, 0.24),
                end: position(0.0, 0.0, 0.16)
            },
            &work(),
            limits(1.0),
            300.0,
            no_source,
            || false
        ),
        Err(MovingPartialSlipError::Invalid(_))
    ));
}
