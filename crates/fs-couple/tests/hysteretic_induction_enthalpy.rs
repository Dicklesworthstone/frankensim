//! Independent P1/play balances through actual source-backed field/heat owners.
//! Manufactured fit coefficients establish accounting, not catalog validity.
#![cfg(feature = "magnetic-hysteresis")]

#[path = "hysteretic_induction/common.rs"]
mod common;
use common::*;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    magnetic::{MagneticRegion, MagneticThermalError},
    material::MaterialNetwork,
};
use fs_em::planar::eddy::hysteretic::PlanarHystereticEddy;
use fs_matdb::{MaterialCard, MaterialStateId};
use fs_material::conductor::ElectricalPropertyBasis;

struct Reference {
    t: f64,
    h: f64,
    z: f64,
    b: f64,
    joule: f64,
    pin: f64,
    numerical: f64,
}
// Analytic P1 moments give projected conductor mass sigma*volume/18. The
// nested scalar solves do not call the production field, material or Newton code.
fn reference(t0: f64, old_h: f64, old_z: f64, old_b: f64) -> Reference {
    let at = |t| {
        let sigma = 1.0 / resistivity(t);
        let h = bisect(-30.0, 30.0, |h| {
            h + sigma / (72.0 * DT) * (MU * h + play(old_z, h) - old_b)
        });
        let z = play(old_z, h);
        let b = MU * h + z;
        (
            h,
            z,
            b,
            sigma * VOLUME / 18.0 * (0.5 * (b - old_b)).powi(2) / DT,
            VOLUME * PIN * (z - old_z).abs(),
        )
    };
    let t = bisect(
        t0,
        t0 + VOLUME * storage(old_h, old_z) / (CP * VOLUME) + 1.0,
        |t| {
            let (_, _, _, joule, pin) = at(t);
            CP * VOLUME * (t - t0) - joule - pin
        },
    );
    let (h, z, b, joule, pin) = at(t);
    let numerical = VOLUME * (h * (b - old_b) - storage(h, z) + storage(old_h, old_z)) - pin;
    Reference {
        t,
        h,
        z,
        b,
        joule,
        pin,
        numerical,
    }
}

#[test]
fn finite_temperature_refresh_preserves_play_history_and_deposits_only_physical_heat() {
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(4);
    let net = network(&cells, &routes);
    let owner = problem(&card, square());
    let mut state = retained(&card, &owner);
    let mut heat = thermal(&curve, VOLUME);
    let (mut old_h, mut old_z, mut old_b) = (12.0, 1.0, 1.12);
    for step in 1..=16 {
        let expected = reference(heat.temperature_k(CELL).unwrap(), old_h, old_z, old_b);
        let before = (
            heat.cells().clone(),
            state.potential_wb_m().to_vec(),
            state.magnetic_state().cell_states().to_vec(),
        );
        let candidate = net
            .step_hysteretic_eddy_candidate(
                &heat,
                step as f64 * DT,
                &state,
                &[0.0; 4],
                controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        near(candidate.thermal().temperature_k(CELL).unwrap(), expected.t);
        near(candidate.field().potential_wb_m()[4], 0.5 * expected.b);
        near(report.electrical.joule_heat_j, expected.joule);
        near(report.electrical.hysteresis_heat_j, expected.pin);
        near(
            report.electrical.numerical_dissipation_j,
            expected.numerical,
        );
        near(report.deposited_joule_j, expected.joule + expected.pin);
        near(
            report.thermal.stored_energy_change_j,
            expected.joule + expected.pin,
        );
        near(report.total_energy_residual_j, 0.0);
        near(report.non_electrical_heat_j, 0.0);
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.potential_wb_m(), before.1);
        assert_eq!(state.magnetic_state().cell_states(), before.2);
        for (cell, direction) in DIRECTIONS.iter().enumerate() {
            let native = &candidate.field().magnetic_state().cell_states()[cell];
            for axis in 0..2 {
                near(native.field_a_m()[axis], expected.h * direction[axis]);
                near(
                    native.polarizations_t()[0][axis],
                    expected.z * direction[axis],
                );
            }
        }
        for source in candidate.field().conductors() {
            near(
                source.conductor.resistivity_ohm_m(),
                resistivity(expected.t),
            );
            assert_eq!(
                source.conductor.source_basis(),
                ElectricalPropertyBasis::Resistivity
            );
        }
        if step == 1 {
            assert!(expected.joule > 0.01 && expected.pin > 0.01 && expected.numerical > 0.01);
            assert!(
                (candidate.thermal().energy_j()
                    - heat.energy_j()
                    - expected.joule
                    - expected.pin
                    - expected.numerical)
                    .abs()
                    > 0.01
            );
            assert!((1.0 / resistivity(expected.t) - 18.0).abs() > 0.1);
        }
        let (next_heat, next_state, _) = candidate.into_parts();
        heat = next_heat;
        state = next_state;
        (old_h, old_z, old_b) = (expected.h, expected.z, expected.b);
    }
    let memory = &state.magnetic_state().cell_states()[0];
    assert!(memory.field_a_m()[0].abs() < 1e-6);
    near(memory.polarizations_t()[0][0], CHI * PIN);
    assert!(state.potential_wb_m()[4] > 0.19);
}

#[test]
fn insulating_hysteresis_receives_pinning_heat_and_rejects_an_unowned_receiver() {
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let owner = problem(&card, paired_squares());
    let zero = owner.demagnetized(0.0, &mut || false).unwrap();
    let state = PlanarHystereticEddy::new(
        &zero,
        vec![],
        &(0..8).collect::<Vec<_>>(),
        &boundary(6),
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap();
    let cells = cells(&card, &curve, &context);
    let mut routes = routes(8);
    let heat = thermal(&curve, 2.0 * VOLUME);
    let drive = [72.0, 72.0, 72.0, 72.0, -72.0, -72.0, -72.0, -72.0];
    let next = network(&cells, &routes)
        .step_hysteretic_eddy_candidate(
            &heat,
            DT,
            &state,
            &drive,
            controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    let expected_pin = 2.0 * VOLUME * PIN;
    near(
        next.thermal().temperature_k(CELL).unwrap(),
        300.0 + PIN / CP,
    );
    near(next.report().electrical.joule_heat_j, 0.0);
    near(next.report().electrical.hysteresis_heat_j, expected_pin);
    near(next.report().deposited_joule_j, expected_pin);
    near(next.report().thermal.stored_energy_change_j, expected_pin);
    near(next.report().total_energy_residual_j, 0.0);
    assert!(next.report().electrical.numerical_dissipation_j > 0.1);
    assert!(next.report().electrical.cells.is_empty());
    assert_eq!(next.report().electrical.magnetic.cells.len(), 8);
    routes[0] = MagneticRegion::Fixed { triangle: 0 };
    // Match only mapped volume so refusal specifically reaches the lossy
    // fixed magnetic cell, rather than an unrelated volume mismatch.
    let incomplete_heat = thermal(&curve, 1.75 * VOLUME);
    assert!(matches!(
        network(&cells, &routes).step_hysteretic_eddy_candidate(
            &incomplete_heat,
            DT,
            &state,
            &drive,
            controls(),
            field_controls(),
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::Source(MagneticThermalError::Invalid(
            "every dissipating magnetic triangle needs an actual thermal heat destination"
        )))
    ));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
    assert!(state.potential_wb_m().iter().all(|a| *a == 0.0));
}

#[test]
fn unsupported_temperature_fit_card_and_late_cancellation_publish_neither_history() {
    for upper in [None, Some(300.0), Some(500.0)] {
        let card = card(upper);
        let curve = caloric(&card);
        let context = point(300.0);
        let owner = problem(&card, square());
        let state = retained(&card, &owner);
        let heat = thermal(&curve, VOLUME);
        let cells = cells(&card, &curve, &context);
        let routes = routes(4);
        let net = network(&cells, &routes);
        let before = (
            heat.cells().clone(),
            state.potential_wb_m().to_vec(),
            state.magnetic_state().cell_states().to_vec(),
        );
        let mut polls = 0;
        let run = net.step_hysteretic_eddy_candidate(
            &heat,
            DT,
            &state,
            &[0.0; 4],
            controls(),
            field_controls(),
            &[0.0],
            || {
                polls += 1;
                false
            },
        );
        if upper != Some(500.0) {
            assert!(run.is_err());
        } else {
            let accepted = run.unwrap();
            assert!(polls > 20);
            let mut seen = 0;
            let cancelled = net.step_hysteretic_eddy_candidate(
                &heat,
                DT,
                &state,
                &[0.0; 4],
                controls(),
                field_controls(),
                &[0.0],
                || {
                    seen += 1;
                    seen >= polls
                },
            );
            assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
            let retry = net
                .step_hysteretic_eddy_candidate(
                    &heat,
                    DT,
                    &state,
                    &[0.0; 4],
                    controls(),
                    field_controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            assert_eq!(retry.thermal().cells(), accepted.thermal().cells());
            assert_eq!(
                retry.field().potential_wb_m(),
                accepted.field().potential_wb_m()
            );
            assert_eq!(
                retry.field().magnetic_state().cell_states(),
                accepted.field().magnetic_state().cell_states()
            );
            let wrong = MaterialCard::assemble(
                MaterialStateId {
                    revision: 1,
                    ..card.id().clone()
                },
                card.claims().clone(),
                vec![],
            )
            .unwrap();
            let wrong_cells = self::cells(&wrong, &curve, &context);
            assert!(
                network(&wrong_cells, &routes)
                    .step_hysteretic_eddy_candidate(
                        &heat,
                        DT,
                        &state,
                        &[0.0; 4],
                        controls(),
                        field_controls(),
                        &[0.0],
                        || false
                    )
                    .is_err()
            );
        }
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.potential_wb_m(), before.1);
        assert_eq!(state.magnetic_state().cell_states(), before.2);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(state.time_s(), 0.0);
    }
}

#[test]
fn complete_heat_provider_keeps_the_native_field_and_both_physical_losses() {
    use fs_couple::circuit_enthalpy::{
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        eddy::thermal::hysteretic_heat_assembly_candidate,
    };
    use fs_exec::{CancelGate, Cx};
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let owner = problem(&card, square());
    let state = retained(&card, &owner);
    let heat = thermal(&curve, VOLUME);
    let cells = cells(&card, &curve, &context);
    let routes = routes(4);
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 59,
                kernel_id: 73,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = assembly.prepare(&cx, &heat).unwrap();
        let next = hysteretic_heat_assembly_candidate(
            &prepared,
            &cx,
            &heat,
            DT,
            &state,
            &[0.0; 4],
            &routes,
            1e-12,
            controls(),
            field_controls(),
            &[0.0],
        )
        .unwrap();
        let expected = reference(300.0, 12.0, 1.0, 1.12);
        near(next.thermal().temperature_k(CELL).unwrap(), expected.t);
        near(next.field().potential_wb_m()[4], 0.5 * expected.b);
        near(
            next.report().deposited_joule_j,
            expected.joule + expected.pin,
        );
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.network.material.cells.len(), 1);
        gate.request();
        assert!(matches!(
            hysteretic_heat_assembly_candidate(
                &prepared,
                &cx,
                &heat,
                DT,
                &state,
                &[0.0; 4],
                &routes,
                1e-12,
                controls(),
                field_controls(),
                &[0.0]
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(heat.time_s(), 0.0);
    });
}
