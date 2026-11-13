//! Actual native DG transients, positive material losses and finite heat cells.
//! The manufactured pulses exercise the admitted realization's algebra, not
//! a new material-validity claim for broadband content outside its source band.
use super::*;
use fs_couple::porous_heating::transient::TransientPorousHeating;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgPorousMaterial, PortId, PortInput};
use fs_duct::{Duct, Segment};
use fs_material::porous::diffusive::PorousDiffusive;
use fs_phs::linear_exact::DenseExponentialBudget;

const PORT: PortId = PortId(81);
const ENDS: [f64; 2] = [LENGTH / 2.0, LENGTH];
fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 128 << 20,
        max_multiply_terms: 80_000_000_000,
    }
}
fn duct(inlet: DgBoundary, outlet: DgBoundary) -> DgDuct<DgPorousMaterial<GasState>> {
    let fluid = gas();
    let model = PorousDiffusive::new(
        &material(),
        &fluid,
        8,
        [
            100.0 * core::f64::consts::TAU,
            4000.0 * core::f64::consts::TAU,
        ],
        5e-5,
    )
    .unwrap();
    DgDuct::from_porous_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: (AREA / core::f64::consts::PI).sqrt(),
                length: LENGTH,
            }],
        },
        &model,
        &fluid,
        DgMesh {
            elements: 2,
            order: 1,
        },
        inlet,
        outlet,
    )
    .unwrap()
}

#[test]
fn g1_native_pulse_heats_its_original_elements_and_preserves_mechanics() {
    let duct = duct(DgBoundary::Closed, DgBoundary::Closed);
    let propagator = duct.propagator(2e-6, 1e-11).unwrap();
    let acoustic =
        duct.project_initial_state(|x| if x < LENGTH / 2.0 { 200.0 } else { 20.0 }, |_| 0.0);
    let thermal = initial(events(&ENDS, &[0.0; 2], &[0.0; 2], 0.0), &SENSIBLE);
    let original_acoustic = acoustic.clone();
    let original_cells = thermal.cells().clone();
    let heating = TransientPorousHeating::try_new(
        &duct,
        &propagator,
        &thermal,
        &[id(0), id(1)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance(
            &acoustic,
            0.0,
            &thermal,
            &[],
            config(propagator.dt() / 7.0),
            || Ok(()),
        )
        .unwrap();
    let mut native = acoustic.clone();
    let native_step = propagator.step(&mut native, &[]).unwrap();
    assert_eq!(result.acoustic, native);
    assert_eq!(result.acoustic_step, native_step);
    assert_eq!(heating.duct().operator_identity(), duct.operator_identity());
    assert!(result.generated.viscous_j > 0.0 && result.generated.thermal_j > 0.0);
    assert!(result.cells[0].generated.total_j() > 20.0 * result.cells[1].generated.total_j());
    let total = result.generated.total_j();
    for cell in &result.cells {
        close(cell.deposited_j, cell.generated.total_j(), 2e-13 * total);
        close(
            result.thermal.cells()[&cell.cell].specific_h_j_kg * original_cells[&cell.cell].mass_kg,
            cell.generated.total_j(),
            3e-13 * total,
        );
        assert_eq!(result.thermal.cells()[&cell.cell].boundary.source_w, 0.0);
    }
    close(result.deposited_j, total, 2e-13 * total);
    close(result.heat_report.external_heat_j, total, 2e-13 * total);
    close(result.heat_report.residual_j, 0.0, 2e-13 * total);
    close(
        acoustic.stored_energy_j() - result.acoustic.stored_energy_j(),
        total,
        3e-8 * total,
    );
    close(result.other_acoustic_loss_and_residual_j, 0.0, 3e-8 * total);
    assert_eq!(result.thermal.time_s(), propagator.dt());
    assert_eq!(acoustic, original_acoustic);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g2_held_port_heat_retains_original_conduction_sources_and_events() {
    let duct = duct(DgBoundary::Port(PORT), DgBoundary::Closed);
    let propagator = duct.propagator(2e-5, 1e-11).unwrap();
    let acoustic = duct.zero_state();
    let mut schedule = events(&ENDS, &[10.0, 2.0], &[2.0, 0.0], 4.0);
    schedule.push(TimedEvent {
        time_s: 0.4 * propagator.dt(),
        event: Event::Boundary(
            id(0),
            Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 3.0,
            },
        ),
    });
    let thermal = initial(schedule, &SENSIBLE);
    let heating = TransientPorousHeating::try_new(
        &duct,
        &propagator,
        &thermal,
        &[id(0), id(1)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    let ports = [PortInput {
        id: PORT,
        effort: 1e-4,
    }];
    let result = heating
        .advance(
            &acoustic,
            0.0,
            &thermal,
            &ports,
            config(propagator.dt() / 9.0),
            || Ok(()),
        )
        .unwrap();
    let native_work = result.acoustic_step.ledger.port_work_step_j;
    close(
        native_work,
        result.acoustic.stored_energy_j() + result.generated.total_j(),
        3e-8 * native_work,
    );
    close(
        result.heat_report.external_heat_j,
        (0.4 * 2.0 + 0.6 * 3.0) * propagator.dt() + result.deposited_j,
        1e-13,
    );
    close(
        result.heat_report.stored_energy_change_j,
        result.heat_report.external_heat_j,
        1e-13,
    );
    assert_eq!(result.heat_report.events, 1);
    assert_eq!(result.thermal.cells()[&id(0)].boundary.source_w, 3.0);
    assert_eq!(thermal.cells()[&id(0)].boundary.source_w, 2.0);
    assert!(result.thermal.cells()[&id(1)].specific_h_j_kg > 2.0);
    assert!(
        result
            .cells
            .iter()
            .all(|cell| cell.generated.total_j() > 0.0)
    );
    let second = heating
        .advance(
            &result.acoustic,
            result.thermal.time_s(),
            &result.thermal,
            &ports,
            config(propagator.dt()),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(second.acoustic.step_count(), 2);
    close(
        second.heat_report.external_heat_j,
        3.0 * propagator.dt() + second.deposited_j,
        1e-13,
    );
}

#[test]
fn g3_matched_boundary_loss_never_becomes_porous_heat() {
    let duct = duct(DgBoundary::Matched, DgBoundary::Matched);
    let propagator = duct.propagator(4e-6, 1e-11).unwrap();
    let acoustic = duct.project_initial_state(|_| 200.0, |_| 0.0);
    let thermal = initial(events(&ENDS, &[150.0; 2], &[0.0; 2], 0.0), &LATENT);
    let heating = TransientPorousHeating::try_new(
        &duct,
        &propagator,
        &thermal,
        &[id(0), id(1)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance(
            &acoustic,
            0.0,
            &thermal,
            &[],
            config(propagator.dt() / 3.0),
            || Ok(()),
        )
        .unwrap();
    assert!(result.other_acoustic_loss_and_residual_j > 0.0);
    let deficit = acoustic.stored_energy_j() - result.acoustic.stored_energy_j();
    assert!(deficit > result.deposited_j);
    close(
        result.heat_report.external_heat_j,
        result.generated.total_j(),
        2e-13 * deficit,
    );
    close(
        deficit,
        result.generated.total_j() + result.other_acoustic_loss_and_residual_j,
        1e-14 * deficit,
    );
    for cell in &result.cells {
        close(
            (result.thermal.cells()[&cell.cell].specific_h_j_kg - 150.0)
                * thermal.cells()[&cell.cell].mass_kg,
            cell.generated.total_j(),
            2e-15,
        );
        assert_eq!(result.thermal.temperature_k(cell.cell).unwrap(), 300.5);
    }
}

#[test]
fn g4_shifted_clocks_keep_exact_heat_and_underresolved_intervals_refuse() {
    let duct = duct(DgBoundary::Closed, DgBoundary::Closed);
    let propagator = duct.propagator(2e-6, 1e-11).unwrap();
    let acoustic = duct.project_initial_state(|_| 100.0, |_| 0.0);
    let thermal = initial(events(&ENDS, &[0.0; 2], &[0.0; 2], 0.0), &SENSIBLE);
    let heating = TransientPorousHeating::try_new(
        &duct,
        &propagator,
        &thermal,
        &[id(0), id(1)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    let shifted = thermal.advance_to(1.0, config(1.0), || Ok(())).unwrap().0;
    let baseline = heating
        .advance(
            &acoustic,
            0.0,
            &thermal,
            &[],
            config(propagator.dt() / 7.0),
            || Ok(()),
        )
        .unwrap();
    let result = heating
        .advance(
            &acoustic,
            shifted.time_s(),
            &shifted,
            &[],
            config(propagator.dt() / 7.0),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(result.generated, baseline.generated);
    assert_eq!(result.acoustic, baseline.acoustic);
    close(
        result.deposited_j,
        result.generated.total_j(),
        3e-13 * result.generated.total_j(),
    );
    close(
        result.heat_report.external_heat_j,
        result.deposited_j,
        3e-13 * result.deposited_j,
    );
    assert_eq!(result.thermal.time_s(), 1.0 + propagator.dt());

    let epoch = 17_179_869_184.0; // 2^34 s: a positive rounded dt is >90% too long.
    let large = thermal
        .advance_to(epoch, config(epoch), || Ok(()))
        .unwrap()
        .0;
    let represented = (epoch + propagator.dt()) - epoch;
    assert!(represented > 1.9 * propagator.dt());
    let original = acoustic.clone();
    assert!(matches!(
        heating.advance(
            &acoustic,
            epoch,
            &large,
            &[],
            config(propagator.dt()),
            || Ok(()),
        ),
        Err(PorousHeatingError::Invalid(_))
    ));
    assert_eq!(large.time_s(), epoch);
    assert_eq!(large.cells(), thermal.cells());
    assert_eq!(acoustic, original);
}

#[test]
fn g4_physical_binding_clock_budget_and_late_refusals_are_atomic() {
    let duct = duct(DgBoundary::Closed, DgBoundary::Closed);
    let propagator = duct.propagator(2e-6, 1e-11).unwrap();
    let acoustic = duct.project_initial_state(|_| 100.0, |_| 0.0);
    let thermal = initial(events(&ENDS, &[0.0; 2], &[0.0; 2], 0.0), &SENSIBLE);
    for destinations in [[id(0), id(0)], [id(1), id(0)]] {
        assert!(
            TransientPorousHeating::try_new(
                &duct,
                &propagator,
                &thermal,
                &destinations,
                budget(),
                || Ok(()),
            )
            .is_err()
        );
    }
    assert!(
        TransientPorousHeating::try_new(
            &duct,
            &propagator,
            &thermal,
            &[id(0), id(1)],
            DenseExponentialBudget {
                max_dense_bytes: 0,
                ..budget()
            },
            || Ok(()),
        )
        .is_err()
    );
    let heating = TransientPorousHeating::try_new(
        &duct,
        &propagator,
        &thermal,
        &[id(0), id(1)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        heating
            .advance(
                &acoustic,
                1.0,
                &thermal,
                &[],
                config(propagator.dt()),
                || Ok(())
            )
            .is_err()
    );
    let mut changed = events(&ENDS, &[0.0; 2], &[0.0; 2], 0.0);
    if let Event::Birth(cell) = &mut changed[0].event {
        cell.mass_kg *= 1.01;
    }
    assert!(
        heating
            .advance(
                &acoustic,
                0.0,
                &initial(changed, &SENSIBLE),
                &[],
                config(propagator.dt()),
                || Ok(()),
            )
            .is_err()
    );
    let original_acoustic = acoustic.clone();
    let original_cells = thermal.cells().clone();
    let mut calls = 0;
    heating
        .advance(
            &acoustic,
            0.0,
            &thermal,
            &[],
            config(propagator.dt()),
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut cancelled_calls = 0;
    let refused = heating.advance(
        &acoustic,
        0.0,
        &thermal,
        &[],
        config(propagator.dt()),
        || {
            cancelled_calls += 1;
            if cancelled_calls == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(
        cancelled_calls, calls,
        "cancel after both candidates have been calculated"
    );
    assert_eq!(acoustic, original_acoustic);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);

    let full = initial(events(&ENDS, &[300.0; 2], &[0.0; 2], 0.0), &LATENT);
    let heating = TransientPorousHeating::try_new(
        &duct,
        &propagator,
        &full,
        &[id(0), id(1)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        heating.advance(&acoustic, 0.0, &full, &[], config(propagator.dt()), || Ok(
            ()
        ),),
        Err(PorousHeatingError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(acoustic, original_acoustic);
    assert_eq!(full.cells()[&id(0)].specific_h_j_kg, 300.0);
}
