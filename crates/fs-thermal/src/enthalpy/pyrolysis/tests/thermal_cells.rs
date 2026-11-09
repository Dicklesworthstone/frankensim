use super::*;
use crate::enthalpy::pyrolysis::implicit::ImplicitSolve;

struct HeatSink {
    basis: ThermalStorageBasis,
    latent: bool,
}
impl ThermalLaw for HeatSink {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, e: f64) -> Result<f64, DepositionError> {
        if !(-50_000.0..=200_000.0).contains(&e) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if !self.latent {
            500.0 + e / 1000.0
        } else if e < 0.0 {
            600.0 + e / 1000.0
        } else if e <= 100_000.0 {
            600.0
        } else {
            600.0 + (e - 100_000.0) / 1000.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 1000.0
    }
}

fn thermal_pair<'a>(sink: &'a dyn ThermalLaw, sink_energy: f64) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (index, material, mass, energy) in [
        (
            1,
            MaterialId(1),
            0.01,
            VIRGIN.specific_enthalpy_j_kg(650.0).unwrap(),
        ),
        (2, MaterialId(2), 0.02, sink_energy),
    ] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(index),
                segment_id: index,
                body_id: index,
                material_id: material,
                mass_kg: mass,
                specific_h_j_kg: energy,
                geometry: CellGeometry {
                    start_m: [index as f64, 0.0, 0.0],
                    end_m: [index as f64 + 1.0, 0.0, 0.0],
                    volume_m3: mass / 1000.0,
                    surface_area_m2: 0.02,
                },
                boundary: boundary(0.0),
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: CellId(1),
            b: CellId(2),
            area_m2: 0.005,
            conductance_w_k: 10.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &VIRGIN), (MaterialId(2), sink)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
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
    .0
}

fn advance_pair<'a>(
    state: &PyrolysisSimulation<'a>,
    implicit: bool,
) -> (PyrolysisSimulation<'a>, PyrolysisAudit) {
    let controls = PyrolysisControls {
        maximum_dt_s: 0.1,
        ..Default::default()
    };
    if implicit {
        let mut workspace = state.implicit_workspace().unwrap();
        let (next, report) = state
            .advance_implicit_to(
                0.1,
                controls,
                ImplicitSolve::default(),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        (next, report.pyrolysis)
    } else {
        state.advance_to(0.1, controls, || Ok(())).unwrap()
    }
}

fn check_account(
    start: &PyrolysisSimulation<'_>,
    next: &PyrolysisSimulation<'_>,
    audit: &PyrolysisAudit,
) {
    close(
        next.energy_j() + audit.volatile_enthalpy_j,
        start.energy_j(),
        2e-7,
    );
    close(
        next.mass_kg() + audit.volatile_mass_kg,
        start.mass_kg(),
        1e-15,
    );
    assert!(audit.absolute_energy_residual_j < 2e-7);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.releases[&CellId(2)], VolatileRelease::default());
    assert_eq!(next.inventories().len(), 1);
    assert!(!next.inventories().contains_key(&CellId(2)));
    assert_eq!(
        next.cells()[&CellId(2)].mass_kg,
        start.cells()[&CellId(2)].mass_kg
    );
    assert_eq!(
        next.cells()[&CellId(2)].geometry,
        start.cells()[&CellId(2)].geometry
    );
}

#[test]
fn g1_reactive_body_and_finite_heat_sink_follow_both_discrete_energy_equations() {
    let gas = |t| Ok(2000.0 * (t - 300.0) + 200_000.0);
    // Closed-mass heat exchange does not require converting u or a declared
    // storage coordinate into h. Only the emitting reactive body transports h.
    for basis in [
        ThermalStorageBasis::Enthalpy,
        ThermalStorageBasis::InternalEnergy,
        ThermalStorageBasis::DeclaredSensibleStorage,
    ] {
        let sink = HeatSink {
            basis,
            latent: false,
        };
        let state = PyrolysisSimulation::new_with_thermal_materials(
            thermal_pair(&sink, 0.0),
            &[model(&VIRGIN, None, &gas, 0.0, 1.0)],
            &[MaterialId(2)],
            || Ok(()),
        )
        .unwrap();
        let reaction = (0.01 * 0.1 / 1.1) * 200_000.0;
        for implicit in [false, true] {
            let (next, audit) = advance_pair(&state, implicit);
            let (hot, cold) = if implicit {
                // Independently eliminate a 2x2 BE heat-capacity system.
                let hot = ((20.0 * 650.0 - reaction) * 21.0 + 20.0 * 500.0) / 440.0;
                (hot, (20.0 * 500.0 + hot) / 21.0)
            } else {
                (650.0 + (-150.0 - reaction) / 20.0, 500.0 + 150.0 / 20.0)
            };
            close(next.temperature_k(CellId(1)).unwrap(), hot, 2e-8);
            close(next.temperature_k(CellId(2)).unwrap(), cold, 2e-8);
            check_account(&state, &next, &audit);
        }
    }
}

#[test]
fn g1_latent_sink_retains_exact_plateau_and_absorbed_energy_during_decomposition() {
    let sink = HeatSink {
        basis: ThermalStorageBasis::Enthalpy,
        latent: true,
    };
    let gas = |t| Ok(2000.0 * (t - 300.0) + 200_000.0);
    let state = PyrolysisSimulation::new_with_thermal_materials(
        thermal_pair(&sink, 50_000.0),
        &[model(&VIRGIN, None, &gas, 0.0, 1.0)],
        &[MaterialId(2)],
        || Ok(()),
    )
    .unwrap();
    let reaction = (0.01 * 0.1 / 1.1) * 200_000.0;
    for implicit in [false, true] {
        let (next, audit) = advance_pair(&state, implicit);
        let hot = if implicit {
            (20.0 * 650.0 + 600.0 - reaction) / 21.0
        } else {
            650.0 + (-50.0 - reaction) / 20.0
        };
        let absorbed = if implicit { hot - 600.0 } else { 50.0 };
        close(next.temperature_k(CellId(1)).unwrap(), hot, 2e-8);
        assert_eq!(next.temperature_k(CellId(2)).unwrap(), 600.0);
        close(
            next.cells()[&CellId(2)].specific_h_j_kg,
            50_000.0 + absorbed / 0.02,
            1e-6,
        );
        assert!(next.cells()[&CellId(2)].specific_h_j_kg > 50_000.0);
        check_account(&state, &next, &audit);
    }
}

fn latent_buffer_reference() -> [f64; 3] {
    // Independent continuous balances: m Cp Tdot = Q - k(T)m L,
    // mdot = -k(T)m, E_buffer_dot = -Q. The original outgoing h carries
    // sensible energy; only its difference L remains after cancelling mdot*h.
    let derivative = |y: [f64; 3]| {
        let [temperature, mass, _] = y;
        let consumption = 2.0 * (-4000.0 / (GAS_CONSTANT * temperature)).exp() * mass;
        let power = 10.0 * (600.0 - temperature);
        [
            (power - 200_000.0 * consumption) / (2000.0 * mass),
            -consumption,
            -power,
        ]
    };
    let mut exact = [650.0, 0.01, 1000.0];
    let dt = 0.5 / 4096.0;
    for _ in 0..4096 {
        let k1 = derivative(exact);
        let k2 = derivative(std::array::from_fn(|i| exact[i] + 0.5 * dt * k1[i]));
        let k3 = derivative(std::array::from_fn(|i| exact[i] + 0.5 * dt * k2[i]));
        let k4 = derivative(std::array::from_fn(|i| exact[i] + dt * k3[i]));
        exact = std::array::from_fn(|i| {
            exact[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    assert!(
        exact[0] < 600.0,
        "the initially cooling buffer must later feed reaction heat"
    );
    exact
}

#[test]
fn g1_activated_reaction_and_latent_buffer_converge_to_independent_energy_odes() {
    let sink = HeatSink {
        basis: ThermalStorageBasis::Enthalpy,
        latent: true,
    };
    let gas = |t| Ok(2000.0 * (t - 300.0) + 200_000.0);
    let mut material = model(&VIRGIN, None, &gas, 0.0, 2.0);
    material.kinetics.activation_energy_j_mol = 4000.0;
    let state = PyrolysisSimulation::new_with_thermal_materials(
        thermal_pair(&sink, 50_000.0),
        &[material],
        &[MaterialId(2)],
        || Ok(()),
    )
    .unwrap();
    let exact = latent_buffer_reference();
    for implicit in [false, true] {
        let mut errors = [0.0; 2];
        for (rung, steps) in [50, 100].into_iter().enumerate() {
            let controls = PyrolysisControls {
                maximum_dt_s: 0.5 / steps as f64,
                ..Default::default()
            };
            let (next, audit) = if implicit {
                let mut workspace = state.implicit_workspace().unwrap();
                let (next, report) = state
                    .advance_implicit_to(
                        0.5,
                        controls,
                        ImplicitSolve::default(),
                        &mut workspace,
                        || Ok(()),
                    )
                    .unwrap();
                (next, report.pyrolysis)
            } else {
                state.advance_to(0.5, controls, || Ok(())).unwrap()
            };
            check_account(&state, &next, &audit);
            assert_eq!(next.temperature_k(CellId(2)).unwrap(), 600.0);
            let observed = [
                next.temperature_k(CellId(1)).unwrap(),
                next.cells()[&CellId(1)].mass_kg,
                next.cells()[&CellId(2)].mass_kg * next.cells()[&CellId(2)].specific_h_j_kg,
            ];
            // Normalize each independently measured physical coordinate.
            errors[rung] = (observed[0] - exact[0]).abs()
                + 10000.0 * (observed[1] - exact[1]).abs()
                + 0.01 * (observed[2] - exact[2]).abs();
        }
        assert!(
            errors[1] < 0.55 * errors[0],
            "implicit={implicit}: {errors:?}"
        );
        assert!(errors[1] < 0.4, "implicit={implicit}: {errors:?}");
    }
}

#[test]
fn adaptive_activated_reaction_resolves_latent_storage_and_reciprocal_heat() {
    use crate::enthalpy::pyrolysis::implicit::adaptive::AdaptiveControls;

    let sink = HeatSink {
        basis: ThermalStorageBasis::Enthalpy,
        latent: true,
    };
    let gas = |t| Ok(2000.0 * (t - 300.0) + 200_000.0);
    let mut material = model(&VIRGIN, None, &gas, 0.0, 2.0);
    material.kinetics.activation_energy_j_mol = 4000.0;
    let state = PyrolysisSimulation::new_with_thermal_materials(
        thermal_pair(&sink, 50_000.0),
        &[material],
        &[MaterialId(2)],
        || Ok(()),
    )
    .unwrap();
    let exact = latent_buffer_reference();
    let mut workspace = state.implicit_workspace().unwrap();
    let mut errors = [0.0; 2];
    let mut work = [0; 2];
    for (rung, scale) in [1.0, 0.01].into_iter().enumerate() {
        let (next, report) = state
            .advance_adaptive_implicit_to(
                0.5,
                PyrolysisControls {
                    maximum_dt_s: 0.5,
                    cell_energy_tolerance_j: 1e-9,
                    ..Default::default()
                },
                ImplicitSolve::default(),
                AdaptiveControls {
                    initial_dt_s: 0.5,
                    storage_absolute_tolerance_j_kg: 100.0 * scale,
                    temperature_absolute_tolerance_k: 0.05 * scale,
                    mass_absolute_tolerance_kg: 1e-6 * scale,
                    volatile_energy_absolute_tolerance_j: scale,
                    relative_tolerance: 0.0,
                    ..Default::default()
                },
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        let audit = &report.physical.pyrolysis;
        check_account(&state, &next, audit);
        assert_eq!(next.temperature_k(CellId(2)).unwrap(), 600.0);
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        let observed = [
            next.temperature_k(CellId(1)).unwrap(),
            next.cells()[&CellId(1)].mass_kg,
            next.cells()[&CellId(2)].mass_kg * next.cells()[&CellId(2)].specific_h_j_kg,
        ];
        errors[rung] = (observed[0] - exact[0]).abs()
            + 10000.0 * (observed[1] - exact[1]).abs()
            + 0.01 * (observed[2] - exact[2]).abs();
        work[rung] = audit.thermal.substeps;
        assert!(next.temperature_k(CellId(1)).unwrap() < 600.0);
    }
    assert!(
        errors[1] < 0.2 * errors[0],
        "errors={errors:?}, steps={work:?}"
    );
    assert!(errors[1] < 0.15, "errors={errors:?}, steps={work:?}");
    assert!(work[1] > work[0]);
}

#[test]
fn heat_only_materials_are_explicit_and_reacting_calorics_require_enthalpy() {
    let sink = HeatSink {
        basis: ThermalStorageBasis::Enthalpy,
        latent: true,
    };
    let gas = |t| Ok(2000.0 * (t - 300.0) + 200_000.0);
    let materials = [model(&VIRGIN, None, &gas, 0.0, 1.0)];
    assert!(
        PyrolysisSimulation::new(thermal_pair(&sink, 50_000.0), &materials, || Ok(())).is_err()
    );
    for declared in [
        &[MaterialId(1)][..],
        &[MaterialId(9)][..],
        &[MaterialId(2), MaterialId(2)][..],
    ] {
        assert!(
            PyrolysisSimulation::new_with_thermal_materials(
                thermal_pair(&sink, 50_000.0),
                &materials,
                declared,
                || Ok(())
            )
            .is_err()
        );
    }
    struct OtherBasis(ThermalStorageBasis);
    impl ThermalLaw for OtherBasis {
        fn storage_basis(&self) -> ThermalStorageBasis {
            self.0
        }
        fn temperature_k(&self, e: f64) -> Result<f64, DepositionError> {
            VIRGIN.temperature_k(e)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            VIRGIN.maximum_temperature_slope()
        }
    }
    impl CondensedCaloric for OtherBasis {
        fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
            VIRGIN.specific_enthalpy_j_kg(t)
        }
    }
    for basis in [
        ThermalStorageBasis::InternalEnergy,
        ThermalStorageBasis::DeclaredSensibleStorage,
    ] {
        let other = OtherBasis(basis);
        assert!(
            PyrolysisSimulation::new(
                thermal(&other, &[650.0], 0.0, 0.0),
                &[model(&other, None, &gas, 0.0, 1.0)],
                || Ok(())
            )
            .is_err()
        );
        assert!(
            PyrolysisSimulation::new(
                thermal(&VIRGIN, &[650.0], 0.0, 0.0),
                &[model(&VIRGIN, Some(&other), &gas, 0.2, 1.0)],
                || Ok(())
            )
            .is_err()
        );
    }
}

#[test]
fn heat_only_caloric_refusal_and_cancellation_leave_reactive_history_unchanged() {
    let sink = HeatSink {
        basis: ThermalStorageBasis::Enthalpy,
        latent: true,
    };
    let gas = |t| Ok(2000.0 * (t - 300.0) + 200_000.0);
    let mut state = PyrolysisSimulation::new_with_thermal_materials(
        thermal_pair(&sink, 199_999.0),
        &[model(&VIRGIN, None, &gas, 0.0, 1.0)],
        &[MaterialId(2)],
        || Ok(()),
    )
    .unwrap();
    let mut hot_boundary = boundary(1e6);
    hot_boundary.ambient_k = 700.0;
    state
        .set_boundary(CellId(2), hot_boundary, || Ok(()))
        .unwrap();
    let cells = state.cells().clone();
    let inventories = state.inventories().clone();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.1,
        ..Default::default()
    };
    assert!(state.advance_to(0.1, controls, || Ok(())).is_err());
    let mut workspace = state.implicit_workspace().unwrap();
    assert!(
        state
            .advance_implicit_to(
                0.1,
                controls,
                ImplicitSolve::default(),
                &mut workspace,
                || Ok(())
            )
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        state
            .advance_to(0.1, controls, || {
                polls += 1;
                if polls == 4 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .err(),
        Some(DepositionError::Cancelled)
    );
    assert_eq!(state.cells(), &cells);
    assert_eq!(state.inventories(), &inventories);
    assert_eq!(state.time_s(), 0.0);
}
