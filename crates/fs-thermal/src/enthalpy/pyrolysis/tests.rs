use super::*;
use crate::enthalpy::{CellGeometry, DepositionBudget, StepConfig, ThermalLink, TimedEvent};

struct Caloric {
    cp: f64,
    offset: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.offset) / self.cp;
        if (100.0..=2000.0).contains(&t) {
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
        if (100.0..=2000.0).contains(&t) {
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

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn thermal<'a>(
    law: &'a dyn CondensedCaloric,
    temperatures: &[f64],
    power: f64,
    conductance: f64,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, &t) in temperatures.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 0.01,
                specific_h_j_kg: law.specific_enthalpy_j_kg(t).unwrap(),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 0.02,
                },
                boundary: boundary(if i == 0 { power } else { 0.0 }),
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
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 100_000,
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
    .0
}
fn model<'a>(
    virgin: &'a dyn CondensedCaloric,
    char_law: Option<&'a dyn CondensedCaloric>,
    gas: &'a dyn VolatileEnthalpy,
    yield_char: f64,
    rate: f64,
) -> PyrolysisMaterial<'a> {
    PyrolysisMaterial {
        material_id: MaterialId(1),
        virgin,
        char_caloric: char_law,
        volatile: gas,
        char_mass_yield: yield_char,
        kinetics: Arrhenius {
            prefactor_s_inv: rate,
            activation_energy_j_mol: 0.0,
        },
        temperature_range_k: [250.0, 1000.0],
        source: "analytic manufactured calorimetry",
    }
}
fn close(a: f64, b: f64, eps: f64) {
    assert!((a - b).abs() <= eps, "{a:.17e} != {b:.17e}");
}

#[test]
fn isothermal_backward_euler_has_actual_positive_inventory_and_first_order_refinement() {
    let law = Arrhenius {
        prefactor_s_inv: 2.0,
        activation_energy_j_mol: 0.0,
    };
    let step = law.isothermal_step(3.0, 600.0, 0.25).unwrap();
    assert_eq!(step.virgin_mass_kg, 2.0);
    assert_eq!(step.consumed_mass_kg, 1.0);
    let solve = |n: usize| {
        let mut mass = 1.0;
        for _ in 0..n {
            mass = law
                .isothermal_step(mass, 600.0, 1.0 / n as f64)
                .unwrap()
                .virgin_mass_kg;
        }
        mass
    };
    let exact = (-2.0_f64).exp();
    let errors = [100, 200, 400].map(|n| solve(n) - exact);
    assert!(errors.iter().all(|e| *e > 0.0));
    assert!(errors[0] / errors[1] > 1.95 && errors[1] / errors[2] > 1.97);
    let huge = law.isothermal_step(1.0, 600.0, 1e100).unwrap();
    assert!(huge.virgin_mass_kg > 0.0);
}

#[test]
fn small_rate_products_preserve_representable_mass_without_fictitious_release() {
    let law = Arrhenius {
        prefactor_s_inv: 9e-17,
        activation_energy_j_mol: 0.0,
    };
    assert_eq!(1.0 + 9e-17, 1.0);
    let step = law.isothermal_step(1.5, 600.0, 1.0).unwrap();
    assert_eq!(step.virgin_mass_kg, 1.5_f64.next_down());
    assert_eq!(step.consumed_mass_kg, 1.5 - step.virgin_mass_kg);
    assert!(step.consumed_mass_kg > 0.0);
    close(
        step.inventory_roundoff_kg,
        step.consumed_mass_kg - 1.5 * 9e-17,
        1e-32,
    );
    let tiny = law.isothermal_step(1.5, 600.0, 0.001).unwrap();
    assert_eq!(tiny.virgin_mass_kg, 1.5);
    assert_eq!(tiny.consumed_mass_kg, 0.0);
    assert!(tiny.inventory_roundoff_kg < 0.0);
}

#[test]
fn unresolved_reaction_heat_accumulates_in_the_budget_without_becoming_a_source() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let initial = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[600.0], 0.0, 0.0),
        &[model(&VIRGIN, None, &gas, 0.0, 1e-16)],
        || Ok(()),
    )
    .unwrap();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.01,
        cell_energy_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        ..Default::default()
    };
    let (next, audit) = initial.advance_to(1.0, controls, || Ok(())).unwrap();
    assert_eq!(next.mass_kg(), initial.mass_kg());
    assert_eq!(audit.volatile_mass_kg, 0.0);
    assert_eq!(audit.volatile_enthalpy_j, 0.0);
    assert!(audit.kinetic_mass_roundoff_kg > 0.0);
    close(audit.kinetic_energy_roundoff_j, 8e-13, 1e-26);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert!(matches!(
        initial.advance_to(
            1.0,
            PyrolysisControls {
                total_energy_tolerance_j: 1e-13,
                ..controls
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(
            "pyrolysis accumulated reaction-energy roundoff"
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn noncharring_adiabatic_energy_uses_volatile_enthalpy_once() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let initial = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[600.0], 0.0, 0.0),
        &[model(&VIRGIN, None, &gas, 0.0, 2.0)],
        || Ok(()),
    )
    .unwrap();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.1,
        ..Default::default()
    };
    let (next, audit) = initial.advance_to(0.1, controls, || Ok(())).unwrap();
    let consumed = 0.01 * (1.0 - 1.0 / 1.2);
    let expected_t = 600.0 - consumed * 800_000.0 / (0.01 * 2000.0);
    close(next.temperature_k(CellId(1)).unwrap(), expected_t, 1e-9);
    close(audit.volatile_mass_kg, consumed, 1e-17);
    close(
        next.mass_kg() + audit.volatile_mass_kg,
        initial.mass_kg(),
        1e-17,
    );
    close(
        next.energy_j() + audit.volatile_enthalpy_j,
        initial.energy_j(),
        1e-8,
    );
    close(audit.thermal.external_heat_j, 0.0, 0.0);
    assert!(audit.thermal.residual_j.abs() < 1e-8);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 600.0);
}

#[test]
fn coupled_temperature_and_inventory_refine_to_the_adiabatic_exact_solution() {
    // Equal vapor/solid Cp and constant k imply dm/dt=-k*m and
    // m*Cp*dT/dt=-k*m*L, so T decreases LINEARLY while m decays exponentially.
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let initial = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[600.0], 0.0, 0.0),
        &[model(&VIRGIN, None, &gas, 0.0, 0.2)],
        || Ok(()),
    )
    .unwrap();
    let solve = |n| {
        initial
            .advance_to(
                1.0,
                PyrolysisControls {
                    maximum_dt_s: 1.0 / n as f64,
                    ..Default::default()
                },
                || Ok(()),
            )
            .unwrap()
            .0
    };
    let states = [10, 20, 40].map(solve);
    let exact_t = 600.0 - 0.2 * 800_000.0 / 2000.0;
    let exact_m = 0.01 * (-0.2_f64).exp();
    let te = states
        .each_ref()
        .map(|s| s.temperature_k(CellId(1)).unwrap() - exact_t);
    let me = states.each_ref().map(|s| s.mass_kg() - exact_m);
    assert!(te[0] / te[1] > 1.97 && te[1] / te[2] > 1.98);
    assert!(me[0] / me[1] > 1.97 && me[1] / me[2] > 1.98);
}

#[test]
fn distinct_char_calorics_and_external_heat_close_independent_algebra() {
    let gas = |t| Ok(2400.0 * (t - 300.0) + 1_000_000.0);
    let initial = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[600.0], 100.0, 0.0),
        &[model(&VIRGIN, Some(&CHAR), &gas, 0.25, 1.0)],
        || Ok(()),
    )
    .unwrap();
    let dt = 0.2;
    let (next, audit) = initial
        .advance_to(
            dt,
            PyrolysisControls {
                maximum_dt_s: dt,
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
    let remaining = 0.01 / (1.0 + dt);
    let consumed = 0.01 - remaining;
    let char_mass = 0.25 * consumed;
    let gas_mass = 0.75 * consumed;
    let sensible =
        0.01 * 2000.0 * 300.0 + 100.0 * dt - char_mass * 400_000.0 - gas_mass * 1_000_000.0;
    let expected = 300.0 + sensible / (remaining * 2000.0 + char_mass * 1200.0 + gas_mass * 2400.0);
    close(next.temperature_k(CellId(1)).unwrap(), expected, 1e-8);
    close(
        next.inventories()[&CellId(1)].char_mass_kg,
        char_mass,
        1e-17,
    );
    close(next.mass_kg() + gas_mass, initial.mass_kg(), 1e-17);
    close(
        next.energy_j() + audit.volatile_enthalpy_j - initial.energy_j(),
        20.0,
        1e-8,
    );
    close(audit.thermal.external_heat_j, 20.0, 1e-12);
}

#[test]
fn disabling_reaction_reuses_original_thermal_transport_exactly_to_roundoff() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let thermal = thermal(&VIRGIN, &[350.0, 500.0, 300.0], 100.0, 4.0);
    let expected = thermal
        .advance_to(
            1.0,
            StepConfig {
                max_dt_s: 0.02,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
    let state =
        PyrolysisSimulation::new(thermal, &[model(&VIRGIN, None, &gas, 0.0, 0.0)], || Ok(()))
            .unwrap();
    let (next, audit) = state
        .advance_to(1.0, PyrolysisControls::default(), || Ok(()))
        .unwrap();
    for &id in next.cells().keys() {
        close(
            next.temperature_k(id).unwrap(),
            expected.0.temperature_k(id).unwrap(),
            1e-8,
        );
    }
    close(
        audit.thermal.external_heat_j,
        expected.1.external_heat_j,
        1e-9,
    );
    assert_eq!(audit.volatile_mass_kg, 0.0);
    assert_eq!(audit.thermal.substeps, expected.1.substeps);
}

#[test]
fn hot_face_drives_resolved_mass_release_with_internal_heat_conservation() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let mut m = model(&VIRGIN, Some(&CHAR), &gas, 0.15, 1e7);
    m.kinetics.activation_energy_j_mol = 90_000.0;
    let state =
        PyrolysisSimulation::new(thermal(&VIRGIN, &[500.0; 4], 200.0, 1.0), &[m], || Ok(()))
            .unwrap();
    let (next, audit) = state
        .advance_to(
            2.0,
            PyrolysisControls {
                total_energy_tolerance_j: 1e-5,
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
    assert!(next.temperature_k(CellId(1)).unwrap() > next.temperature_k(CellId(4)).unwrap());
    assert!(audit.releases[&CellId(1)].mass_kg > audit.releases[&CellId(4)].mass_kg);
    close(
        next.mass_kg() + audit.volatile_mass_kg,
        state.mass_kg(),
        1e-14,
    );
    close(
        next.energy_j() + audit.volatile_enthalpy_j - state.energy_j(),
        400.0,
        1e-5,
    );
    close(audit.thermal.external_heat_j, 400.0, 1e-9);
}

#[test]
fn common_formation_datum_shift_preserves_temperature_and_reaction() {
    let a = Caloric {
        cp: 2000.0,
        offset: 0.0,
    };
    let b = Caloric {
        cp: 2000.0,
        offset: 1_000_000.0,
    };
    let ga = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let gb = |t| Ok(2000.0 * (t - 300.0) + 1_800_000.0);
    let x = PyrolysisSimulation::new(
        thermal(&a, &[600.0], 80.0, 0.0),
        &[model(&a, None, &ga, 0.0, 1.0)],
        || Ok(()),
    )
    .unwrap();
    let y = PyrolysisSimulation::new(
        thermal(&b, &[600.0], 80.0, 0.0),
        &[model(&b, None, &gb, 0.0, 1.0)],
        || Ok(()),
    )
    .unwrap();
    let (xn, xa) = x
        .advance_to(0.2, PyrolysisControls::default(), || Ok(()))
        .unwrap();
    let (yn, ya) = y
        .advance_to(0.2, PyrolysisControls::default(), || Ok(()))
        .unwrap();
    close(
        xn.temperature_k(CellId(1)).unwrap(),
        yn.temperature_k(CellId(1)).unwrap(),
        1e-9,
    );
    close(xa.volatile_mass_kg, ya.volatile_mass_kg, 1e-16);
    close(
        ya.volatile_enthalpy_j - xa.volatile_enthalpy_j,
        1_000_000.0 * xa.volatile_mass_kg,
        1e-8,
    );
}

#[test]
fn cancellation_work_limit_domain_and_boundary_refusals_are_atomic() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let state = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[600.0, 500.0], 100.0, 1.0),
        &[model(&VIRGIN, None, &gas, 0.0, 1.0)],
        || Ok(()),
    )
    .unwrap();
    let cells = state.cells().clone();
    let mut polls = 0;
    let error = state.advance_to(1.0, PyrolysisControls::default(), || {
        polls += 1;
        if polls == 50 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(error, Err(DepositionError::Cancelled)));
    assert_eq!(state.cells(), &cells);
    assert_eq!(state.time_s(), 0.0);
    let mut bounded = state.clone();
    bounded.thermal.budget.max_substeps = 1;
    assert!(matches!(
        bounded.advance_to(1.0, PyrolysisControls::default(), || Ok(())),
        Err(DepositionError::Budget(_))
    ));
    assert_eq!(bounded.cells(), &cells);
    let mut invalid = state.clone();
    assert!(
        invalid
            .set_boundary(CellId(1), boundary(f64::NAN), || Ok(()))
            .is_err()
    );
    assert_eq!(invalid.cells(), &cells);
    invalid
        .set_boundary(CellId(1), boundary(-1e8), || Ok(()))
        .unwrap();
    assert!(matches!(
        invalid.advance_to(0.02, PyrolysisControls::default(), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(invalid.time_s(), 0.0);
}

#[test]
fn law_identity_residue_calorics_and_exothermic_model_are_not_guessed() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 800_000.0);
    let counterfeit = Caloric {
        cp: 2000.0,
        offset: 0.0,
    };
    assert!(
        PyrolysisSimulation::new(
            thermal(&VIRGIN, &[500.0], 0.0, 0.0),
            &[model(&counterfeit, None, &gas, 0.0, 1.0)],
            || Ok(())
        )
        .is_err()
    );
    assert!(
        PyrolysisSimulation::new(
            thermal(&VIRGIN, &[500.0], 0.0, 0.0),
            &[model(&VIRGIN, None, &gas, 0.2, 1.0)],
            || Ok(())
        )
        .is_err()
    );
    let exothermic = |t| Ok(2000.0 * (t - 300.0) - 100.0);
    assert!(
        PyrolysisSimulation::new(
            thermal(&VIRGIN, &[500.0], 0.0, 0.0),
            &[model(&VIRGIN, None, &exothermic, 0.0, 1.0)],
            || Ok(())
        )
        .is_err()
    );
}
