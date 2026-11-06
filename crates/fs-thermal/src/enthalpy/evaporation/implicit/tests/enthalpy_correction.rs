//! Complete vapor enthalpy on the original conservative thermal residual.
use super::*;
use crate::enthalpy::evaporation::implicit::network::VaporEnthalpyStepCoefficients;
use std::cell::Cell;

fn binding<'l>(
    id: u64,
    law: &'l Interface,
    film: GasFilm,
    datum: f64,
) -> EvaporatingCellBinding<'l, Interface> {
    EvaporatingCellBinding {
        id: CellId(id),
        law,
        film,
        pure_liquid_enthalpy_range_j_kg: range(datum),
    }
}

fn unchanged(actual: &DepositionSimulation<'_>, original: &DepositionSimulation<'_>) {
    assert_eq!(actual.time_s, original.time_s);
    assert_eq!(actual.cursor, original.cursor);
    assert_eq!(actual.cells, original.cells);
    assert_eq!(actual.links, original.links);
    assert_eq!(actual.covered_areas, original.covered_areas);
    assert_eq!(actual.schedule, original.schedule);
}

#[test]
fn complete_vapor_streams_preserve_latent_heat_contacts_and_large_caloric_datums() {
    let dt = 500.0;
    let temperatures = [305.0, 306.0];
    let films = BTreeMap::from([(CellId(1), film(0.0)), (CellId(2), film(8000.0))]);
    let areas = [0.007, 0.008];
    let increments = [0.1 * CP * 5.0, -0.2 * CP * 4.0];
    // Evaporation has an admitted negative caloric-source mismatch. Incoming
    // vapor also carries its sensible enthalpy from a hotter upstream gas.
    let excess = [-1250.0, 2000.0 * (320.0 - temperatures[1]) + 450.0];
    let dm = [
        transfer(&VARYING, temperatures[0], films[&CellId(1)], areas[0], dt),
        transfer(&VARYING, temperatures[1], films[&CellId(2)], areas[1], dt),
    ];
    assert!(dm[0] > 0.0 && dm[1] < 0.0);
    let latent = [LATENT - 2000.0 * 5.0, LATENT - 2000.0 * 6.0];
    let required = [
        increments[0] + dm[0] * (latent[0] + excess[0]),
        increments[1] + dm[1] * (latent[1] + excess[1]),
    ];
    let mut reference = None;
    for datum in [0.0, 1e12] {
        let caloric = Caloric {
            datum,
            fixed_mass: None,
        };
        let mut a = cell(1, 0.1, 300.0, datum, required[0] / dt - (-1.05 + 3.0 + 0.4));
        a.boundary.ambient_k = 290.0;
        a.boundary.htc_w_m2_k = 10.0;
        a.boundary.reservoir = Some(Reservoir {
            area_m2: 0.001,
            temperature_k: 320.0,
            conductance_w_k: 0.2,
        });
        let mut b = cell(2, 0.2, 310.0, datum, required[1] / dt - (-0.96 - 0.4));
        b.boundary.htc_w_m2_k = 20.0;
        let old = simulation(&caloric, vec![a, b], Some(0.4), vec![]);
        let step = EvaporatingNetworkBackwardEulerStep::new(
            &old,
            &[
                binding(1, &VARYING, films[&CellId(1)], datum),
                binding(2, &VARYING, films[&CellId(2)], datum),
            ],
            dt,
            || Ok(()),
        )
        .unwrap();
        let calls = Cell::new(0);
        let trial = step
            .evaluate_with_vapor_enthalpy(
                &increments,
                |_| Ok::<_, ThermalEvaporationError>(films.clone()),
                |thermal, liquids, resolved| -> Result<_, ThermalEvaporationError> {
                    calls.set(calls.get() + 1);
                    assert_eq!(resolved, &films);
                    for i in 0..2 {
                        let id = CellId(i as u64 + 1);
                        close(thermal.temperature_k(id)?, temperatures[i], 1e-12);
                        close(liquids[&id].area_m2, areas[i], 1e-17);
                        assert_eq!(thermal.cells()[&id].mass_kg, old.cells()[&id].mass_kg);
                        close(
                            liquids[&id].mass_kg,
                            old.cells()[&id].mass_kg - dm[i],
                            1e-15,
                        );
                    }
                    Ok(VaporEnthalpyStepCoefficients {
                        heat: coefficients("complete vapor stream"),
                        vapor_enthalpy_excess_j_kg: BTreeMap::from([
                            (CellId(1), excess[0]),
                            (CellId(2), excess[1]),
                        ]),
                    })
                },
                || Ok(()),
            )
            .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(trial.payload(), &"complete vapor stream");
        assert!(trial.maximum_cell_residual_j() < 1e-7);
        let (next, payload, audit) = trial.accept(1e-7, 1e-6, 1e-12, || Ok(())).unwrap();
        assert_eq!(payload, "complete vapor stream");
        for i in 0..2 {
            let id = CellId(i as u64 + 1);
            let receipt = &audit.thermal.surfaces[&id];
            close(receipt.latent_energy_j, dm[i] * latent[i], 1e-7);
            close(
                audit.vapor_enthalpy_correction_j[&id],
                dm[i] * excess[i],
                1e-9,
            );
            let vapor = receipt.transferred_mass_kg
                * (datum + CP * (temperatures[i] - 300.0) + latent[i] + excess[i]);
            close(
                receipt.vapor_energy_j,
                vapor,
                1e-7 + 8.0 * f64::EPSILON * vapor.abs(),
            );
            close(receipt.supplied_heat_j, required[i], 1e-7);
            close(receipt.energy_residual_j, 0.0, 1e-7);
            close(next.temperature_k(id).unwrap(), temperatures[i], 1e-12);
            assert_eq!(next.cells()[&id].geometry, old.cells()[&id].geometry);
            assert_eq!(next.cells()[&id].boundary, old.cells()[&id].boundary);
        }
        close(audit.thermal.external_heat_j, required.iter().sum(), 1e-7);
        close(audit.thermal.energy_residual_j, 0.0, 1e-7);
        close(audit.thermal.mass_balance_residual_kg, 0.0, 1e-14);
        let invariant = (
            audit.thermal.external_heat_j,
            audit.thermal.energy_residual_j,
            audit.vapor_enthalpy_correction_j,
        );
        if let Some(reference) = &reference {
            assert_eq!(&invariant, reference);
        } else {
            reference = Some(invariant);
        }
    }
}

#[test]
fn hot_incoming_vapor_releases_its_sensible_energy_in_the_adiabatic_endpoint() {
    let law = Interface {
        pressure_slope: 0.0,
        latent_slope: -2000.0,
    };
    let mut humid = film(6000.0);
    humid.temperature_k = 330.0;
    let dt = 500.0;
    let mass = 0.1;
    let raw_dm = transfer(&law, 300.0, humid, 0.01, dt);
    let dm = mass - (mass - raw_dm);
    assert!(dm < 0.0);
    let gas_cp = 2000.0;
    let mismatch = 250.0;
    // Kirchhoff-consistent L'(T)=Cp_v-Cp_l and incoming h_v(Tgas) make the
    // closed adiabatic endpoint an analytic mixing-enthalpy equation.
    let rise = -dm * (LATENT + gas_cp * (330.0 - 300.0) + mismatch) / ((mass - dm) * CP);
    let old = simulation(&CALORIC, vec![cell(1, mass, 300.0, 0.0, 0.0)], None, vec![]);
    let step =
        EvaporatingNetworkBackwardEulerStep::new(&old, &[binding(1, &law, humid, 0.0)], dt, || {
            Ok(())
        })
        .unwrap();
    let increments = [mass * CP * rise];
    let trial = step
        .evaluate_with_vapor_enthalpy(
            &increments,
            |_| Ok::<_, ThermalEvaporationError>(BTreeMap::from([(CellId(1), humid)])),
            |thermal, _, _| -> Result<_, ThermalEvaporationError> {
                Ok(VaporEnthalpyStepCoefficients {
                    heat: coefficients(()),
                    vapor_enthalpy_excess_j_kg: BTreeMap::from([(
                        CellId(1),
                        gas_cp * (330.0 - thermal.temperature_k(CellId(1))?) + mismatch,
                    )]),
                })
            },
            || Ok(()),
        )
        .unwrap();
    let (next, (), audit) = trial.accept(1e-7, 1e-7, 1e-12, || Ok(())).unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 300.0 + rise, 1e-11);
    let expected_extra = dm * (gas_cp * (330.0 - 300.0 - rise) + mismatch);
    close(
        audit.vapor_enthalpy_correction_j[&CellId(1)],
        expected_extra,
        1e-9,
    );
    assert!(expected_extra < -1.0);
    let receipt = &audit.thermal.surfaces[&CellId(1)];
    close(receipt.latent_energy_j, dm * (LATENT - 2000.0 * rise), 1e-8);
    close(
        audit.thermal.stored_energy_change_j + receipt.vapor_energy_j,
        0.0,
        1e-7,
    );
    assert_eq!(audit.thermal.external_heat_j, 0.0);

    let missing_sensible = step
        .evaluate(
            &increments,
            |_, _| Ok::<_, ThermalEvaporationError>(coefficients(())),
            || Ok(()),
        )
        .unwrap();
    close(
        missing_sensible.audit().energy_residual_j,
        -expected_extra,
        1e-7,
    );
    assert!(matches!(
        missing_sensible.accept(1e-6, 1e-6, 1e-12, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::EnergyBalance(_)
        ))
    ));
}

#[test]
fn exact_correction_maps_and_individual_energy_gates_cannot_be_bypassed() {
    let dt = 1000.0;
    let films = BTreeMap::from([(CellId(1), film(0.0)), (CellId(2), film(7840.0))]);
    let cells = (1..=2)
        .map(|id| {
            let raw_dm = transfer(&CONSTANT, 300.0, films[&CellId(id)], 0.01, dt);
            let dm = 0.1 - (0.1 - raw_dm);
            cell(id, 0.1, 300.0, 0.0, dm * LATENT / dt)
        })
        .collect();
    let old = simulation(&CALORIC, cells, None, vec![]);
    let before = old.clone();
    let step = EvaporatingNetworkBackwardEulerStep::new(
        &old,
        &[
            binding(1, &CONSTANT, films[&CellId(1)], 0.0),
            binding(2, &CONSTANT, films[&CellId(2)], 0.0),
        ],
        dt,
        || Ok(()),
    )
    .unwrap();
    let evaluate = |excess| {
        step.evaluate_with_vapor_enthalpy(
            &[0.0, 0.0],
            |_| Ok::<_, ThermalEvaporationError>(films.clone()),
            |_, _, _| -> Result<_, ThermalEvaporationError> {
                Ok(VaporEnthalpyStepCoefficients {
                    heat: coefficients(()),
                    vapor_enthalpy_excess_j_kg: excess,
                })
            },
            || Ok(()),
        )
    };
    for invalid in [
        BTreeMap::new(),
        BTreeMap::from([(CellId(1), 0.0)]),
        BTreeMap::from([(CellId(1), 0.0), (CellId(3), 0.0)]),
        BTreeMap::from([(CellId(1), 0.0), (CellId(2), 0.0), (CellId(3), 0.0)]),
        BTreeMap::from([(CellId(1), f64::NAN), (CellId(2), 0.0)]),
        BTreeMap::from([(CellId(1), 0.0), (CellId(2), f64::INFINITY)]),
        BTreeMap::from([(CellId(1), f64::NEG_INFINITY), (CellId(2), 0.0)]),
    ] {
        assert!(matches!(
            evaluate(invalid),
            Err(ThermalEvaporationError::Thermal(
                DepositionError::InvalidInput(_)
            ))
        ));
    }
    let corrected_zero = evaluate(BTreeMap::from([(CellId(1), 0.0), (CellId(2), 0.0)])).unwrap();
    let ordinary = step
        .evaluate_with_films(
            &[0.0, 0.0],
            |_| Ok::<_, ThermalEvaporationError>(films.clone()),
            |_, _, _| Ok::<_, ThermalEvaporationError>(coefficients(())),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(corrected_zero.residual_j(), ordinary.residual_j());
    assert_eq!(&corrected_zero.audit().thermal, ordinary.audit());
    let (corrected, (), audit) = corrected_zero.accept(1e-7, 1e-7, 1e-12, || Ok(())).unwrap();
    let (uncorrected, (), _) = ordinary.accept(1e-7, 1e-7, 1e-12, || Ok(())).unwrap();
    assert_eq!(corrected.cells(), uncorrected.cells());
    assert!(
        audit
            .vapor_enthalpy_correction_j
            .values()
            .all(|q| *q == 0.0)
    );

    let opposing = evaluate(BTreeMap::from([(CellId(1), 1e5), (CellId(2), 1e5)])).unwrap();
    close(opposing.audit().thermal.energy_residual_j, 0.0, 1e-8);
    assert!(opposing.maximum_cell_residual_j() > 100.0);
    assert!(matches!(
        opposing.accept(1e-6, 1e-6, 1e-12, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::EnergyBalance(_)
        ))
    ));
    unchanged(&old, &before);
}

#[test]
fn source_errors_original_materials_pending_events_and_final_cancel_remain_atomic() {
    let dt = 10.0;
    let raw_dm = transfer(&CONSTANT, 300.0, film(0.0), 0.01, dt);
    let dm = 0.1 - (0.1 - raw_dm);
    let liquid = cell(1, 0.1, 300.0, 0.0, dm * (LATENT + 1000.0) / dt);
    let mut future_boundary = liquid.boundary;
    future_boundary.source_w = 99.0;
    let old = simulation(
        &CALORIC,
        vec![liquid.clone()],
        None,
        vec![TimedEvent {
            time_s: dt,
            event: Event::Boundary(CellId(1), future_boundary),
        }],
    );
    let before = old.clone();
    let selected = [binding(1, &CONSTANT, film(0.0), 0.0)];
    let films = BTreeMap::from([(CellId(1), film(0.0))]);
    let step = EvaporatingNetworkBackwardEulerStep::new(&old, &selected, dt, || Ok(())).unwrap();
    let source_error = step.evaluate_with_vapor_enthalpy(
        &[0.0],
        |_| Ok::<_, ThermalEvaporationError>(films.clone()),
        |_, _, _| -> Result<VaporEnthalpyStepCoefficients<()>, ThermalEvaporationError> {
            Err(DepositionError::InvalidInput("original vapor enthalpy source").into())
        },
        || Ok(()),
    );
    assert!(matches!(
        source_error,
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("original vapor enthalpy source")
        ))
    ));
    let evaluate = || {
        step.evaluate_with_vapor_enthalpy(
            &[0.0],
            |_| Ok::<_, ThermalEvaporationError>(films.clone()),
            |_, _, _| -> Result<_, ThermalEvaporationError> {
                Ok(VaporEnthalpyStepCoefficients {
                    heat: coefficients(()),
                    vapor_enthalpy_excess_j_kg: BTreeMap::from([(CellId(1), 1000.0)]),
                })
            },
            || Ok(()),
        )
        .unwrap()
    };
    let mut checks = 0;
    let (accepted, (), _) = evaluate()
        .accept(1e-7, 1e-7, 1e-12, || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(accepted.time_s(), dt);
    assert_eq!(accepted.cursor, old.cursor);
    assert_eq!(accepted.cells()[&CellId(1)].boundary, liquid.boundary);
    assert!(accepted.cells()[&CellId(1)].mass_kg < liquid.mass_kg);
    let mut cancelled_checks = 0;
    let cancelled = evaluate().accept(1e-7, 1e-7, 1e-12, || {
        cancelled_checks += 1;
        if cancelled_checks == checks {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(ThermalEvaporationError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(cancelled_checks, checks);
    unchanged(&old, &before);

    let fixed = Caloric {
        datum: 0.0,
        fixed_mass: Some(0.1),
    };
    let fixed_old = simulation(&fixed, vec![liquid], None, vec![]);
    let fixed_step =
        EvaporatingNetworkBackwardEulerStep::new(&fixed_old, &selected, dt, || Ok(())).unwrap();
    let mut correction_calls = 0;
    let refused = fixed_step.evaluate_with_vapor_enthalpy(
        &[0.0],
        |_| Ok::<_, ThermalEvaporationError>(films.clone()),
        |_, _, _| -> Result<_, ThermalEvaporationError> {
            correction_calls += 1;
            Ok(VaporEnthalpyStepCoefficients {
                heat: coefficients(()),
                vapor_enthalpy_excess_j_kg: BTreeMap::from([(CellId(1), 1000.0)]),
            })
        },
        || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("original fixed inventory changed")
        ))
    ));
    assert_eq!(correction_calls, 0);
    assert_eq!(fixed_old.cells()[&CellId(1)].mass_kg, 0.1);
}
