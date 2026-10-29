use super::*;
use crate::enthalpy::{
    CellGeometry, DepositionBudget, Reservoir, StepConfig, ThermalLink, TimedEvent,
};

struct Caloric {
    cp: f64,
    offset: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.offset) / self.cp;
        if (200.0..=1200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
impl CondensedCaloric for Caloric {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if (200.0..=1200.0).contains(&t) {
            Ok(self.offset + self.cp * (t - 300.0))
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}
static VIRGIN: Caloric = Caloric {
    cp: 2000.0,
    offset: 0.0,
};
static CHAR: Caloric = Caloric {
    cp: 1200.0,
    offset: 400_000.0,
};

fn setup<'a>(
    law: &'a dyn CondensedCaloric,
    gas: &'a dyn VolatileEnthalpy,
    temperatures: &[f64],
    masses: &[f64],
    power: f64,
    conductance: f64,
    char_law: Option<&'a dyn CondensedCaloric>,
    yield_char: f64,
    rate: f64,
) -> PyrolysisSimulation<'a> {
    let mut events = Vec::new();
    for (i, (&t, &mass)) in temperatures.iter().zip(masses).enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: mass,
                specific_h_j_kg: law.specific_enthalpy_j_kg(t).unwrap(),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: mass / 1000.0,
                    surface_area_m2: 0.02,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: if i == 0 { power } else { 0.0 },
                },
            }),
        });
    }
    for i in 1..temperatures.len() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                area_m2: 0.005,
                conductance_w_k: conductance,
            }),
        });
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 128,
            max_links: 127,
            max_events: 255,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: law,
            char_caloric: char_law,
            volatile: gas,
            char_mass_yield: yield_char,
            kinetics: Arrhenius {
                prefactor_s_inv: rate,
                activation_energy_j_mol: 0.0,
            },
            temperature_range_k: [250.0, 1000.0],
            source: "independent manufactured backward-Euler heat/reaction oracle",
        }],
        || Ok(()),
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} differs from {b:.17e} by {:.3e} > {tolerance:.3e}",
        (a - b).abs()
    );
}

fn advance<'a>(
    state: &PyrolysisSimulation<'a>,
    end: f64,
    maximum_dt: f64,
) -> (PyrolysisSimulation<'a>, ImplicitReport) {
    let mut workspace = state.implicit_workspace().unwrap();
    state
        .advance_implicit_to(
            end,
            PyrolysisControls {
                maximum_dt_s: maximum_dt,
                ..Default::default()
            },
            ImplicitSolve::default(),
            &mut workspace,
            || Ok(()),
        )
        .unwrap()
}

#[test]
fn tiny_cell_stiff_conduction_matches_independent_two_capacity_solution() {
    let gas = |t| Ok(2000.0 * (t - 300.0));
    let mut initial = setup(
        &VIRGIN,
        &gas,
        &[300.0, 900.0],
        &[1e-8, 0.01],
        0.0,
        1000.0,
        None,
        0.0,
        0.0,
    );
    initial.thermal.budget.max_substeps = 4;
    assert!(matches!(
        initial.advance_to(
            0.1,
            PyrolysisControls {
                maximum_dt_s: 0.1,
                ..Default::default()
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    let (next, report) = advance(&initial, 0.1, 0.1);
    let [c0, c1] = [1e-8 * 2000.0, 0.01 * 2000.0];
    let difference = 600.0 / (1.0 + 0.1 * 1000.0 * (1.0 / c0 + 1.0 / c1));
    let t0 = (c0 * 300.0 + c1 * 900.0 - c1 * difference) / (c0 + c1);
    close(next.temperature_k(CellId(1)).unwrap(), t0, 1e-8);
    close(
        next.temperature_k(CellId(2)).unwrap(),
        t0 + difference,
        1e-8,
    );
    assert_eq!(report.pyrolysis.thermal.substeps, 1);
    assert!(report.linear_solves > 0 && report.linear_solves < 8);
    close(next.energy_j(), initial.energy_j(), 1e-8);
    assert_eq!(report.pyrolysis.volatile_mass_kg, 0.0);
    assert_eq!(report.pyrolysis.thermal.external_heat_j, 0.0);
}

#[test]
fn stiff_neighbor_transfer_heating_and_reaction_close_one_endpoint_system() {
    let latent = 200_000.0;
    let gas = |t| Ok(2000.0 * (t - 300.0) + latent);
    let masses = [1e-5, 0.01];
    let initial = setup(
        &VIRGIN,
        &gas,
        &[900.0, 300.0],
        &masses,
        5000.0,
        100.0,
        None,
        0.0,
        2.0,
    );
    let (next, report) = advance(&initial, 0.1, 0.1);
    // Constant rate and equal virgin/volatile Cp make the full coupled problem
    // linear: original initial capacity carries retained + emitted sensible heat.
    let c = masses.map(|m| 2000.0 * m);
    let consumed = masses.map(|m| m - m / 1.2);
    let source = [500.0 - latent * consumed[0], -latent * consumed[1]];
    let gdt = 10.0;
    let b = [c[0] * 900.0 + source[0], c[1] * 300.0 + source[1]];
    let determinant = c[0] * c[1] + gdt * (c[0] + c[1]);
    let expected = [
        (b[0] * (c[1] + gdt) + gdt * b[1]) / determinant,
        (b[1] * (c[0] + gdt) + gdt * b[0]) / determinant,
    ];
    for (i, expected) in expected.into_iter().enumerate() {
        close(
            next.temperature_k(CellId(i as u64 + 1)).unwrap(),
            expected,
            2e-8,
        );
        close(
            next.inventories()[&CellId(i as u64 + 1)].virgin_mass_kg,
            masses[i] / 1.2,
            1e-17,
        );
    }
    close(report.pyrolysis.thermal.external_heat_j, 500.0, 1e-11);
    close(
        next.energy_j() - initial.energy_j() + report.pyrolysis.volatile_enthalpy_j,
        500.0,
        1e-8,
    );
    close(
        next.mass_kg() + report.pyrolysis.volatile_mass_kg,
        initial.mass_kg(),
        1e-17,
    );
    assert!(report.pyrolysis.maximum_cell_energy_residual_j <= 1e-8);
}

#[test]
fn distinct_char_calorics_and_datum_shift_preserve_physical_solution() {
    let solve = |offset: f64| {
        let virgin = Caloric { cp: 2000.0, offset };
        let char_law = Caloric {
            cp: 1200.0,
            offset: offset + 400_000.0,
        };
        let gas = |t| Ok(offset + 2000.0 * (t - 300.0) + 800_000.0);
        let initial = setup(
            &virgin,
            &gas,
            &[600.0, 700.0],
            &[0.01, 0.01],
            100.0,
            20.0,
            Some(&char_law),
            0.2,
            2.0,
        );
        let (next, report) = advance(&initial, 0.1, 0.1);
        let endpoint = next.inventories().values().copied().collect::<Vec<_>>();
        for cell in &endpoint {
            assert!(cell.char_mass_kg > 0.0 && cell.virgin_mass_kg > 0.0);
        }
        close(
            next.mass_kg() + report.pyrolysis.volatile_mass_kg,
            initial.mass_kg(),
            1e-17,
        );
        close(
            next.energy_j() - initial.energy_j() + report.pyrolysis.volatile_enthalpy_j,
            10.0,
            2e-8,
        );
        (endpoint, report.pyrolysis)
    };
    let a = solve(0.0);
    let b = solve(1e8);
    for (x, y) in a.0.iter().zip(&b.0) {
        close(x.temperature_k, y.temperature_k, 2e-8);
        close(x.virgin_mass_kg, y.virgin_mass_kg, 1e-17);
        close(x.char_mass_kg, y.char_mass_kg, 1e-17);
    }
    close(
        b.1.volatile_enthalpy_j - a.1.volatile_enthalpy_j,
        1e8 * a.1.volatile_mass_kg,
        1e-8,
    );
}

#[test]
fn robin_and_reservoir_heat_use_actual_covered_area_and_endpoint_temperature() {
    let gas = |t| Ok(2000.0 * (t - 300.0));
    let mut initial = setup(&VIRGIN, &gas, &[300.0], &[0.01], 0.0, 0.0, None, 0.0, 0.0);
    initial
        .set_boundary(
            CellId(1),
            Boundary {
                ambient_k: 900.0,
                htc_w_m2_k: 100.0,
                reservoir: Some(Reservoir {
                    area_m2: 0.005,
                    conductance_w_k: 3.0,
                    temperature_k: 500.0,
                }),
                source_w: 100.0,
            },
            || Ok(()),
        )
        .unwrap();
    let (next, report) = advance(&initial, 0.5, 0.5);
    let ambient_g = (0.02 - 0.005) * 100.0;
    let expected = (20.0 * 300.0 + 0.5 * (ambient_g * 900.0 + 3.0 * 500.0 + 100.0))
        / (20.0 + 0.5 * (ambient_g + 3.0));
    close(next.temperature_k(CellId(1)).unwrap(), expected, 1e-9);
    close(
        report.pyrolysis.thermal.external_heat_j,
        20.0 * (expected - 300.0),
        1e-8,
    );
}

#[test]
fn temporal_refinement_converges_to_independent_conduction_and_reaction_trajectory() {
    let gas = |t| Ok(2000.0 * (t - 300.0));
    let initial = setup(
        &VIRGIN,
        &gas,
        &[400.0, 800.0],
        &[0.01, 0.01],
        0.0,
        2.0,
        None,
        0.0,
        0.0,
    );
    let exact_t0 = 600.0 - 200.0 * (-0.2_f64).exp();
    let errors = [0.1, 0.05, 0.025].map(|dt| {
        let (next, _) = advance(&initial, 1.0, dt);
        (next.temperature_k(CellId(1)).unwrap() - exact_t0).abs()
    });
    assert!(errors[0] / errors[1] > 1.97 && errors[1] / errors[2] > 1.98);
    let reactive = setup(
        &VIRGIN,
        &gas,
        &[600.0],
        &[0.01],
        0.0,
        0.0,
        Some(&CHAR),
        0.0,
        2.0,
    );
    let mass_errors = [0.1, 0.05, 0.025].map(|dt| {
        let (next, report) = advance(&reactive, 1.0, dt);
        close(next.temperature_k(CellId(1)).unwrap(), 600.0, 1e-9);
        close(
            next.mass_kg() + report.pyrolysis.volatile_mass_kg,
            0.01,
            1e-16,
        );
        next.mass_kg() - 0.01 * (-2.0_f64).exp()
    });
    assert!(mass_errors[0] / mass_errors[1] > 1.9 && mass_errors[1] / mass_errors[2] > 1.95);
}

#[test]
fn source_domain_work_and_late_cancellation_refuse_atomically_and_workspace_replays() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let initial = setup(
        &VIRGIN,
        &gas,
        &[600.0, 900.0],
        &[0.01, 0.01],
        300.0,
        5.0,
        None,
        0.0,
        0.1,
    );
    let original_cells = initial.cells().clone();
    let original_inventory = initial.inventories().clone();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.1,
        ..Default::default()
    };
    let mut workspace = initial.implicit_workspace().unwrap();
    assert!(matches!(
        initial.advance_implicit_to(
            1.0,
            controls,
            ImplicitSolve {
                maximum_iterations: 1,
                ..Default::default()
            },
            &mut workspace,
            || Ok(())
        ),
        Err(ImplicitError::Iterations)
    ));
    let mut polls = 0;
    assert!(matches!(
        initial.advance_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            &mut workspace,
            || {
                polls += 1;
                if polls == 250 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(ImplicitError::Cancelled)
    ));
    assert_eq!(polls, 250);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original_cells);
    assert_eq!(initial.inventories(), &original_inventory);
    let replay = initial
        .advance_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let direct = advance(&initial, 1.0, 0.1);
    assert_eq!(replay.0.cells(), direct.0.cells());
    assert_eq!(replay.0.inventories(), direct.0.inventories());
    assert_eq!(replay.1, direct.1);
    let mut overheated = initial.clone();
    overheated
        .set_boundary(
            CellId(1),
            Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 1e12,
            },
            || Ok(()),
        )
        .unwrap();
    assert!(
        overheated
            .advance_implicit_to(
                1.0,
                controls,
                ImplicitSolve::default(),
                &mut workspace,
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(overheated.time_s(), 0.0);
}

#[test]
fn activated_reaction_and_conduction_refine_toward_original_explicit_transport() {
    let law = Caloric {
        cp: 1500.0,
        offset: 0.0,
    };
    let gas = |t| Ok(1500.0 * (t - 300.0) + 1_000_000.0);
    let seed = setup(
        &law,
        &gas,
        &[600.0, 500.0],
        &[0.002, 0.002],
        400.0,
        1.0,
        None,
        0.0,
        1.0,
    );
    let material = PyrolysisMaterial {
        kinetics: Arrhenius {
            prefactor_s_inv: 1e8,
            activation_energy_j_mol: 100_000.0,
        },
        ..seed.materials[&MaterialId(1)]
    };
    let mut initial =
        PyrolysisSimulation::new(seed.thermal.clone(), &[material], || Ok(())).unwrap();
    initial.thermal.budget.max_substeps = 10_000;
    let exact = initial
        .advance_to(
            0.5,
            PyrolysisControls {
                maximum_dt_s: 0.0001,
                cell_energy_tolerance_j: 1e-10,
                total_energy_tolerance_j: 1e-5,
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let errors = [0.04, 0.02, 0.01].map(|dt| {
        let (next, report) = advance(&initial, 0.5, dt);
        assert!(report.pyrolysis.volatile_mass_kg > 0.0);
        close(
            next.mass_kg() + report.pyrolysis.volatile_mass_kg,
            initial.mass_kg(),
            1e-16,
        );
        close(
            next.energy_j() - initial.energy_j() + report.pyrolysis.volatile_enthalpy_j,
            200.0,
            1e-6,
        );
        [CellId(1), CellId(2)]
            .into_iter()
            .map(|id| (next.temperature_k(id).unwrap() - exact.temperature_k(id).unwrap()).abs())
            .fold(0.0, f64::max)
    });
    assert!(
        errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8,
        "{errors:?}"
    );
}
