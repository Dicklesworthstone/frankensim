//! Original heat transport while the finite atmosphere has no active interface.
use super::*;
use fs_thermal::enthalpy::evaporation::implicit::EvaporatingCellBinding;

fn heat_payload(
    thermal: &DepositionSimulation<'_>,
    liquids: &BTreeMap<CellId, LiquidSurfaceState>,
) -> Result<BoundaryExchangeStepCoefficients<(f64, f64, f64)>, ThermalEvaporationError> {
    assert!(liquids.is_empty());
    Ok(coefficients((
        thermal.time_s(),
        thermal.temperature_k(CellId(1))?,
        thermal.temperature_k(CellId(2))?,
    )))
}

#[test]
fn empty_phase_roster_matches_original_nonlinear_heat_network_and_freezes_complete_gas_state() {
    let source = catalog(-1e7);
    let initial_gas = gas(&source, 340.0, 5000.0, 0.001);
    let dt = 20.0;
    let required = [
        LIQUID_MASS * delta_h(300.0, 305.0),
        LIQUID_MASS * delta_h(320.0, 315.0),
    ];
    let contact_w = 2.0 * (315.0 - 305.0);
    let old = thermal(
        vec![
            cell(1, 300.0, required[0] / dt - contact_w, 0.0012),
            cell(2, 320.0, required[1] / dt + contact_w, 0.0012),
        ],
        Some((0.0002, 2.0)),
    );
    let snapshot = old.clone();
    let gas_snapshot = initial_gas.clone();
    let mut config = controls();
    // No active gas/liquid interface exists, so phase temperature and gas-root
    // work limits must not constrain the independently admitted native heat.
    config.phase.maximum_temperature_increment_k = 1e-3;
    config.phase.gas.gas.maximum_iterations = 1;
    let no_bindings: [EvaporatingCellBinding<'_>; 0] = [];
    let (native, native_payload, native_report) = old
        .solve_implicit_evaporating_network_step(
            &no_bindings,
            dt,
            config.thermal,
            heat_payload,
            || Ok(()),
        )
        .unwrap();
    let candidate = solve_deposition_with_gas(
        &old,
        &[],
        &initial_gas,
        &source,
        dt,
        config,
        |thermal, liquids, gas| {
            assert_eq!(gas, &initial_gas);
            heat_payload(thermal, liquids).map_err(ImplicitGasError::from)
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(candidate.thermal().cells(), native.cells());
    assert_eq!(candidate.thermal().links(), native.links());
    assert_eq!(candidate.payload(), &native_payload);
    assert_eq!(candidate.report().thermal.thermal, native_report.physical);
    near(
        candidate.thermal().temperature_k(CellId(1)).unwrap(),
        305.0,
        1e-9,
    );
    near(
        candidate.thermal().temperature_k(CellId(2)).unwrap(),
        315.0,
        1e-9,
    );
    near(
        candidate.report().external_heat_j,
        required[0] + required[1],
        1e-8,
    );
    near(
        candidate.thermal().energy_j() - old.energy_j(),
        required[0] + required[1],
        1e-8,
    );
    assert_eq!(candidate.thermal().mass_kg(), old.mass_kg());
    assert_eq!(candidate.thermal().time_s(), dt);
    assert_eq!(candidate.gas(), &initial_gas);
    assert_eq!(
        candidate.gas().species_inventory_kg().unwrap(),
        initial_gas.species_inventory_kg().unwrap()
    );
    assert!(std::ptr::eq(candidate.catalog(), &source));
    assert!(candidate.report().thermal.thermal.surfaces.is_empty());
    assert!(
        candidate
            .report()
            .thermal
            .vapor_enthalpy_correction_j
            .is_empty()
    );
    assert!(candidate.report().sensible_heat_into_liquid_j.is_empty());
    assert_eq!(candidate.report().gas_iterations, 0);
    assert_eq!(candidate.report().gas_mass_residual_kg, 0.0);
    assert_eq!(candidate.report().gas_energy_residual_j, 0.0);
    let gas_step = candidate.gas_exchange().step();
    assert!(gas_step.species.is_empty());
    assert_eq!(gas_step.mass_into_gas_kg, 0.0);
    assert_eq!(gas_step.energy_into_gas_j, 0.0);
    assert_eq!(gas_step.stored_energy_change_j, 0.0);
    assert_eq!(gas_step.iterations, 0);
    unchanged(&old, &snapshot, &initial_gas, &gas_snapshot);
}

#[test]
fn idle_phase_keeps_catalog_identity_source_freshness_work_and_cancellation_gates() {
    let source = catalog(0.0);
    let changed_source = catalog(1000.0);
    let initial_gas = gas(&source, 300.0, 1000.0, 0.001);
    let old = thermal(vec![cell(1, 300.0, 5.0, 0.001)], None);
    let snapshot = old.clone();
    let gas_snapshot = initial_gas.clone();
    let calls = Cell::new(0usize);
    let result = solve_deposition_with_gas(
        &old,
        &[],
        &initial_gas,
        &changed_source,
        1.0,
        controls(),
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Phase(GasEvaporationError::Gas(
            GasEnergyError::Invalid(
                "boundary catalog substitutes an active caloric or transport source"
            )
        )))
    ));
    assert_eq!(calls.get(), 0);

    let mut polls = 0usize;
    let accepted = solve_deposition_with_gas(
        &old,
        &[],
        &initial_gas,
        &source,
        1.0,
        controls(),
        |_, liquids, gas| {
            assert!(liquids.is_empty());
            assert_eq!(gas, &initial_gas);
            calls.set(calls.get() + 1);
            Ok(coefficients(calls.get()))
        },
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(*accepted.payload(), calls.get());
    assert_eq!(accepted.report().gas_iterations, 0);
    let final_source_call = calls.get();
    calls.set(0);
    let result = solve_deposition_with_gas(
        &old,
        &[],
        &initial_gas,
        &source,
        1.0,
        controls(),
        |_, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() == final_source_call {
                Err(DepositionError::InvalidInput("final original idle heat source refused").into())
            } else {
                Ok(coefficients(()))
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("final original idle heat source refused")
        )))
    ));
    assert_eq!(calls.get(), final_source_call);

    let mut limited = controls();
    limited.thermal.maximum_evaluations = 1;
    assert!(matches!(
        solve_deposition_with_gas(
            &old,
            &[],
            &initial_gas,
            &source,
            1.0,
            limited,
            empty_heat,
            || Ok(())
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::Budget("implicit finite-gas residual evaluations")
        )))
    ));
    let mut count = 0usize;
    assert!(matches!(
        solve_deposition_with_gas(
            &old,
            &[],
            &initial_gas,
            &source,
            1.0,
            controls(),
            empty_heat,
            || {
                count += 1;
                if count == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(count, polls);
    unchanged(&old, &snapshot, &initial_gas, &gas_snapshot);
}

#[test]
fn idle_phase_retains_native_empty_geometry_and_pending_event_semantics() {
    let source = catalog(0.0);
    let initial_gas = gas(&source, 300.0, 1000.0, 0.001);
    let budget = DepositionBudget {
        max_cells: 1,
        max_links: 0,
        max_events: 2,
        max_substeps: 8,
    };
    let empty =
        DepositionSimulation::new(vec![], &[(MaterialId(1), &ORIGINAL_LIQUID)], budget).unwrap();
    // fs-thermal backward_euler.rs:173-182 refuses the empty active-cell roster with
    // `DepositionError::InvalidInput`; ef5ab372af reworded its text. No refusal-specific
    // variant exists: `InvalidInput` also carries, on this same path, the pending-event
    // refusal of endpoint.rs:83-86, so the text still names which refusal this is.
    assert!(matches!(
        solve_deposition_with_gas(
            &empty,
            &[],
            &initial_gas,
            &source,
            1.0,
            controls(),
            empty_heat,
            || Ok(())
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(
                "implicit step needs active cells, positive duration and a matching later target"
            )
        )))
    ));
    let born = cell(1, 300.0, 5.0, 0.001);
    let mut changed = born.boundary;
    changed.source_w = 9.0;
    let original = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(born),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), changed),
            },
        ],
        &[(MaterialId(1), &ORIGINAL_LIQUID)],
        budget,
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
    .0;
    let snapshot = original.clone();
    let gas_snapshot = initial_gas.clone();
    assert!(matches!(
        solve_deposition_with_gas(
            &original,
            &[],
            &initial_gas,
            &source,
            2.0,
            controls(),
            empty_heat,
            || Ok(())
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(
                "process pending events or split implicit step at the next event"
            )
        )))
    ));
    let candidate = solve_deposition_with_gas(
        &original,
        &[],
        &initial_gas,
        &source,
        1.0,
        controls(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(
        candidate.thermal().pending_events(),
        original.pending_events()
    );
    assert_eq!(
        candidate.thermal().cells()[&CellId(1)].boundary.source_w,
        5.0
    );
    assert_eq!(candidate.gas(), &initial_gas);
    near(candidate.report().external_heat_j, 5.0, 1e-10);
    assert!(matches!(
        solve_deposition_with_gas(
            candidate.thermal(),
            &[],
            candidate.gas(),
            &source,
            2.0,
            controls(),
            empty_heat,
            || Ok(())
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(
                "process pending events or split implicit step at the next event"
            )
        )))
    ));
    unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
}
