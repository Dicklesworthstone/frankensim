//! Independent manufactured physical references for live frame thermodynamics.
use super::*;

fn live_config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_connections: 18,
        max_setup_terms: 32_768,
        ..config()
    }
}
fn live_specs() -> [Spec; 2] {
    let mut specs = specimens();
    specs[0].modulus_slope = -2.0;
    specs[1].modulus_slope = -3.0;
    for spec in &mut specs {
        spec.cp *= 0.01;
    }
    specs
}

#[test]
fn affine_modulus_changes_axial_and_bending_motion_with_exact_caloric_internal_energy() {
    for bending in [false, true] {
        let mut endpoints = Vec::new();
        for cached in [false, true] {
            for temperature_dependent in [false, true] {
                let specs = live_specs().map(|mut s| {
                    if !temperature_dependent {
                        s.modulus_slope = 0.0;
                    }
                    s
                });
                fixture(
                    specs,
                    bending,
                    false,
                    0.0,
                    |component, mut thermal, network| {
                        let gate = CancelGate::new_clock_free();
                        let mut owner = component
                            .into_thermoelastic_system(
                                &thermal,
                                network,
                                &DESTINATIONS,
                                0.05,
                                live_config(),
                                cached.then(prepared),
                                &gate,
                            )
                            .unwrap();
                        let mass = RHO * areas().iter().sum::<f64>() * LENGTH
                            / if bending { 105.0 } else { 3.0 };
                        let geometry = if bending {
                            moments().map(|i| 4.0 * i / LENGTH)
                        } else {
                            areas().map(|a| a / LENGTH)
                        };
                        let viscosity: [f64; 2] = core::array::from_fn(|i| {
                            geometry[i]
                                * if bending {
                                    specs[i].bending
                                } else {
                                    specs[i].axial
                                }
                        });
                        let damping: f64 = viscosity.iter().sum();
                        let (mut q, mut v, mut temperatures) = (0.01, 0.03, [300.0; 2]);
                        let dof = if bending { 5 } else { 3 };
                        let shape = owner.basis().modes()[0].displacement[dof];
                        let initial_total = 0.5 * mass * v * v
                            + (0..2)
                                .map(|i| {
                                    RHO * areas()[i]
                                        * LENGTH
                                        * specs[i].cp
                                        * (temperatures[i] - 200.0)
                                        + 0.5
                                            * (YOUNG - 300.0 * specs[i].modulus_slope)
                                            * geometry[i]
                                            * q
                                            * q
                                })
                                .sum::<f64>();
                        for _ in 0..32 {
                            let e: [f64; 2] = core::array::from_fn(|i| {
                                YOUNG + specs[i].modulus_slope * (temperatures[i] - 300.0)
                            });
                            let k = (0..2).map(|i| e[i] * geometry[i]).sum::<f64>();
                            let (next_q, next_v, next_energy) =
                                oscillator_from(mass, k, damping, 1.0 / f64::from(RATE), q, v);
                            let loss = 0.5 * (mass * v * v + k * q * q) - next_energy;
                            let report = owner
                                .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                                .unwrap();
                            for i in 0..2 {
                                let heat = loss * viscosity[i] / damping;
                                let reversible = 0.5
                                    * temperatures[i]
                                    * specs[i].modulus_slope
                                    * geometry[i]
                                    * (next_q * next_q - q * q);
                                temperatures[i] +=
                                    (heat + reversible) / (RHO * areas()[i] * LENGTH * specs[i].cp);
                                near(report.layers[i].material_heat_j, heat, 2e-11);
                                near(report.layers[i].reversible_heat_j, reversible, 2e-11);
                                near(
                                    thermal.temperature_k(CELLS[i]).unwrap(),
                                    temperatures[i],
                                    2e-7,
                                );
                                let coefficient =
                                    YOUNG + specs[i].modulus_slope * (temperatures[i] - 300.0);
                                let state = report.layers[i].states[1];
                                near(
                                    state.free_energy_j,
                                    0.5 * coefficient * geometry[i] * next_q * next_q,
                                    2e-11,
                                );
                                near(
                                    state.entropy_j_k,
                                    -0.5 * specs[i].modulus_slope * geometry[i] * next_q * next_q,
                                    2e-13,
                                );
                                near(
                                    state.internal_energy_j,
                                    0.5 * (YOUNG - 300.0 * specs[i].modulus_slope)
                                        * geometry[i]
                                        * next_q
                                        * next_q,
                                    2e-11,
                                );
                                let ka = owner
                                    .mechanical()
                                    .connection_coefficients(3 + 3 * i)
                                    .unwrap()
                                    .stiffness_n_m;
                                near(ka, coefficient * areas()[i] / LENGTH, 2e-8);
                            }
                            q = next_q;
                            v = next_v;
                            let actual = owner.mechanical().components()[0].states()[0];
                            near(shape * actual.displacement_m_sqrt_kg, q, 2e-9);
                            near(shape * actual.velocity_m_sqrt_kg_per_s, v, 2e-8);
                            let thermal_u = (0..2)
                                .map(|i| {
                                    let cell = &thermal.cells()[&CELLS[i]];
                                    cell.mass_kg * cell.specific_h_j_kg
                                })
                                .sum::<f64>();
                            near(
                                owner.mechanical_internal_energy_j().unwrap() + thermal_u,
                                initial_total,
                                2e-9,
                            );
                            near(report.total_energy_residual_j, 0.0, 2e-10);
                            near(report.entropy.accounting_residual_j_k, 0.0, 2e-12);
                        }
                        endpoints.push((
                            shape
                                * owner.mechanical().components()[0].states()[0]
                                    .displacement_m_sqrt_kg,
                            thermal.temperature_k(CELLS[0]).unwrap(),
                        ));
                    },
                );
            }
        }
        // DATA alone changes both real motion and the same caloric history.
        assert!((endpoints[0].0 - endpoints[1].0).abs() > if bending { 1e-8 } else { 1e-6 });
        near(endpoints[0].0, endpoints[2].0, 2e-9);
        near(endpoints[1].0, endpoints[3].0, 2e-9);
    }
}

fn chain_specs() -> [Spec; 4] {
    let base = specimens();
    core::array::from_fn(|i| Spec {
        modulus_slope: [-2.0, 1.0, -1.0, -2.5][i],
        alpha: [1e-4, -2e-4, 3e-4, -1e-4][i],
        cp: [0.2, 0.25, 0.15, 0.3][i],
        axial: [2.0, 8.0, 5.0, 3.0][i],
        curved_eta: i == 2,
        ..base[i % 2]
    })
}
const POWERS: [f64; 4] = [0.1, 0.2, 0.3, 0.4];
fn chain<R>(
    run: impl FnOnce(
        FrameVibrationComponent,
        DepositionSimulation<'_>,
        &MaterialNetwork<'_>,
        &[FrameLayerCell],
    ) -> R,
) -> R {
    let specs = chain_specs();
    let sources: Vec<_> = specs
        .iter()
        .enumerate()
        .map(|(i, &s)| material(i, s))
        .collect();
    let context = point(300.0);
    let requirements = [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("axial_eta", DynViscosity::DIMS),
        ("bending_eta", DynViscosity::DIMS),
        (
            fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
        ),
    ]
    .map(|(n, d)| ScalarPropertyRequirement::try_new(n, d, ScalarAdmissibility::Finite).unwrap());
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    for i in 0..3 {
        frame.add_node(i as f64 * LENGTH, 0.0);
        frame.fix(i, [i == 0, true, true]);
    }
    for element in 0..2 {
        let layers = (0..2)
            .map(|i| ConcentricMaterialLayer {
                outer_radius_m: RADIUS[i],
                material: resolve_material_state_point(
                    &sources[2 * element + i].0,
                    &context,
                    &requirements,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap(),
            })
            .collect();
        let section = FrameViscousSection::try_new(
            ConcentricElasticSection::try_new(0.0, layers).unwrap(),
            "axial_eta",
            "bending_eta",
        )
        .unwrap();
        frame
            .add_material_viscous_element(element, element + 1, 5, section)
            .unwrap();
    }
    let basis = FrameVibrationBasis::from_frame(
        &frame,
        2,
        FrameVibrationBudget {
            max_frame_dofs: 9,
            max_dense_entries: 81,
        },
        || false,
    )
    .unwrap();
    let mut component = basis
        .coupled_component(
            RATE,
            &[0.0; 2],
            ModalAcousticTimeBudget::audible_reference(),
            72,
            || false,
        )
        .unwrap();
    let mu = RHO * areas().iter().sum::<f64>() * LENGTH / 6.0;
    let project = |u: [f64; 2], phi: &[f64]| {
        phi[3] * mu * (4.0 * u[0] + u[1]) + phi[6] * mu * (u[0] + 2.0 * u[1])
    };
    component
        .restore_states(
            &basis
                .modes()
                .iter()
                .map(|m| ModalAcousticState {
                    displacement_m_sqrt_kg: project([0.008, 0.003], &m.displacement),
                    velocity_m_sqrt_kg_per_s: project([0.01, -0.015], &m.displacement),
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let events = (0..4)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(i as u64),
                mass_kg: RHO * areas()[i % 2] * LENGTH,
                specific_h_j_kg: 100.0 * specs[i].cp,
                geometry: CellGeometry {
                    start_m: [(i / 2) as f64 * LENGTH, 0.0, 0.0],
                    end_m: [(i / 2 + 1) as f64 * LENGTH, 0.0, 0.0],
                    volume_m3: areas()[i % 2] * LENGTH,
                    surface_area_m2: 2.0 * core::f64::consts::PI * RADIUS[i % 2] * LENGTH,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: POWERS[i],
                },
            }),
        })
        .collect();
    let caloric: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, (_, curve))| (MaterialId(i as u64), curve.curve() as &dyn ThermalLaw))
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &caloric,
        DepositionBudget {
            max_cells: 4,
            max_links: 0,
            max_events: 4,
            max_substeps: 1024,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let cells: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, (card, curve))| MaterialCell {
            cell: CellId(i as u64),
            material_id: MaterialId(i as u64),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect();
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let destinations: Vec<_> = (0..4)
        .map(|i| FrameLayerCell {
            element: i / 2,
            layer: i % 2,
            cell: CellId(i as u64),
        })
        .collect();
    run(component, thermal, &network, &destinations)
}
fn potential(spec: Spec, layer: usize, t: f64, extension: f64) -> [f64; 3] {
    let e = extension / LENGTH - spec.alpha * (t - 300.0);
    let young = YOUNG + spec.modulus_slope * (t - 300.0);
    let volume = areas()[layer] * LENGTH;
    let psi = 0.5 * volume * young * e * e;
    let entropy = volume * (young * spec.alpha * e - 0.5 * spec.modulus_slope * e * e);
    [psi, entropy, psi + t * entropy]
}
fn chain_reference(qv: [f64; 4], t: [f64; 4]) -> [f64; 6] {
    let specs = chain_specs();
    let mut k = [0.0; 2];
    let mut c = [0.0; 2];
    let mut load = [0.0; 2];
    for i in 0..4 {
        let young = YOUNG + specs[i].modulus_slope * (t[i] - 300.0);
        k[i / 2] += young * areas()[i % 2] / LENGTH;
        let eta = specs[i].axial
            * if specs[i].curved_eta {
                0.005 * t[i] - 0.5
            } else {
                1.0
            };
        c[i / 2] += eta * areas()[i % 2] / LENGTH;
        load[i / 2] += young * areas()[i % 2] * specs[i].alpha * (t[i] - 300.0);
    }
    let mu = RHO * areas().iter().sum::<f64>() * LENGTH / 6.0;
    let derivative = |s: [f64; 6]| -> [f64; 6] {
        let f0 = -k[0] * s[0] - k[1] * (s[0] - s[1]) - c[0] * s[2] - c[1] * (s[2] - s[3]) + load[0]
            - load[1];
        let f1 = -k[1] * (s[1] - s[0]) - c[1] * (s[3] - s[2]) + load[1];
        [
            s[2],
            s[3],
            (2.0 * f0 - f1) / (7.0 * mu),
            (-f0 + 4.0 * f1) / (7.0 * mu),
            c[0] * s[2] * s[2],
            c[1] * (s[3] - s[2]).powi(2),
        ]
    };
    let mut y = [qv[0], qv[1], qv[2], qv[3], 0.0, 0.0];
    let h = 1.0 / f64::from(RATE) / 64.0;
    for _ in 0..64 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + h * c[i]));
        for i in 0..6 {
            y[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn unequal_layer_expansion_and_modulus_refresh_full_two_dof_stiffness_without_reset() {
    let mut endpoints = Vec::new();
    for cached in [false, true] {
        chain(|component, mut thermal, network, destinations| {
            let gate = CancelGate::new_clock_free();
            let mut owner = component
                .into_thermoelastic_system(
                    &thermal,
                    network,
                    destinations,
                    0.05,
                    live_config(),
                    cached.then(prepared),
                    &gate,
                )
                .unwrap();
            let specs = chain_specs();
            let mut temperatures = [300.0; 4];
            let mut qv = [0.008, 0.003, 0.01, -0.015];
            let mut controls = solve();
            controls.maximum_cells = 4;
            let initial = owner.mechanical_internal_energy_j().unwrap()
                + thermal
                    .cells()
                    .values()
                    .map(|c| c.mass_kg * c.specific_h_j_kg)
                    .sum::<f64>();
            for step in 0..24 {
                let next = chain_reference(qv, temperatures);
                let old_extensions = [qv[0], qv[1] - qv[0]];
                let new_extensions = [next[0], next[1] - next[0]];
                let eta: [f64; 4] = core::array::from_fn(|i| {
                    specs[i].axial
                        * if specs[i].curved_eta {
                            0.005 * temperatures[i] - 0.5
                        } else {
                            1.0
                        }
                });
                let report = owner
                    .step(&mut thermal, &[0.0; 2], controls, &[0.0; 4], &gate)
                    .unwrap();
                for i in 0..4 {
                    let e = i / 2;
                    let layer = i % 2;
                    let initial = potential(specs[i], layer, temperatures[i], old_extensions[e]);
                    let frozen = potential(specs[i], layer, temperatures[i], new_extensions[e]);
                    let work = frozen[0] - initial[0];
                    let loss = next[4 + e] * eta[i] * areas()[layer]
                        / (eta[2 * e] * areas()[0] + eta[2 * e + 1] * areas()[1]);
                    let capacity = RHO * areas()[layer] * LENGTH * specs[i].cp;
                    let residual = |t: f64| {
                        capacity * (t - temperatures[i])
                            + potential(specs[i], layer, t, new_extensions[e])[2]
                            - initial[2]
                            - work
                            - loss
                            - POWERS[i] / f64::from(RATE)
                    };
                    let mut temperature = temperatures[i];
                    for _ in 0..8 {
                        let slope =
                            (residual(temperature + 0.001) - residual(temperature - 0.001)) / 0.002;
                        temperature -= residual(temperature) / slope;
                    }
                    temperatures[i] = temperature;
                    let expected = potential(specs[i], layer, temperature, new_extensions[e]);
                    let actual = report.layers[i].states[1];
                    near(
                        thermal.temperature_k(CellId(i as u64)).unwrap(),
                        temperature,
                        3e-7,
                    );
                    near(report.layers[i].stress_work_j, work, 3e-10);
                    near(report.layers[i].material_heat_j, loss, 3e-10);
                    near(actual.free_energy_j, expected[0], 3e-10);
                    near(actual.entropy_j_k, expected[1], 3e-12);
                    near(actual.internal_energy_j, expected[2], 3e-10);
                    let port = owner
                        .mechanical()
                        .connection_coefficients(6 + 3 * i)
                        .unwrap();
                    near(
                        port.rest_extension_m,
                        specs[i].alpha * (temperature - 300.0) * LENGTH,
                        2e-10,
                    );
                }
                qv = [next[0], next[1], next[2], next[3]];
                let actual: [f64; 4] = core::array::from_fn(|i| {
                    owner
                        .basis()
                        .modes()
                        .iter()
                        .zip(owner.mechanical().components()[0].states())
                        .map(|(mode, state)| {
                            mode.displacement[if i % 2 == 0 { 3 } else { 6 }]
                                * if i < 2 {
                                    state.displacement_m_sqrt_kg
                                } else {
                                    state.velocity_m_sqrt_kg_per_s
                                }
                        })
                        .sum()
                });
                for i in 0..4 {
                    near(actual[i], qv[i], 3e-8);
                }
                assert_eq!(owner.mechanical().samples_rendered(), step + 1);
                let total = owner.mechanical_internal_energy_j().unwrap()
                    + thermal
                        .cells()
                        .values()
                        .map(|c| c.mass_kg * c.specific_h_j_kg)
                        .sum::<f64>();
                near(
                    total,
                    initial + (step + 1) as f64 / f64::from(RATE) * POWERS.iter().sum::<f64>(),
                    3e-9,
                );
                near(report.entropy.accounting_residual_j_k, 0.0, 3e-12);
                near(report.total_energy_residual_j, 0.0, 2e-10);
            }
            assert!((temperatures[0] - temperatures[1]).abs() > 0.01);
            near(
                owner
                    .mechanical()
                    .connection_coefficients(3)
                    .unwrap()
                    .damping_n_s_m,
                (specs[2].axial * (0.005 * temperatures[2] - 0.5) * areas()[0]
                    + specs[3].axial * areas()[1])
                    / LENGTH,
                2e-8,
            );
            // Unequal members create a nonzero off-diagonal K in the INITIAL
            // mass basis. A scalar frequency retune would drop this term.
            let phi = owner.basis().modes();
            let k: [f64; 2] = core::array::from_fn(|e| {
                (0..2)
                    .map(|i| {
                        (YOUNG + specs[2 * e + i].modulus_slope * (temperatures[2 * e + i] - 300.0))
                            * areas()[i]
                            / LENGTH
                    })
                    .sum()
            });
            let cross = k[0] * phi[0].displacement[3] * phi[1].displacement[3]
                + k[1]
                    * (phi[0].displacement[6] - phi[0].displacement[3])
                    * (phi[1].displacement[6] - phi[1].displacement[3]);
            assert!(cross.abs() > 0.01);
            endpoints.push((
                owner.mechanical().components()[0].states().to_vec(),
                thermal
                    .cells()
                    .values()
                    .map(|c| c.specific_h_j_kg)
                    .collect::<Vec<_>>(),
            ));
        });
    }
    for (a, b) in endpoints[0].0.iter().zip(&endpoints[1].0) {
        near(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg, 3e-10);
        near(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s, 3e-9);
    }
    for (a, b) in endpoints[0].1.iter().zip(&endpoints[1].1) {
        near(*a, *b, 1e-8);
    }
}

#[test]
fn thermal_strain_and_cancelled_refusals_preserve_live_coefficients_and_clean_retry() {
    fixture(
        live_specs(),
        false,
        false,
        0.0,
        |component, mut thermal, network| {
            let gate = CancelGate::new_clock_free();
            let mut other = component
                .basis()
                .coupled_component(
                    RATE,
                    &[0.0],
                    ModalAcousticTimeBudget::audible_reference(),
                    18,
                    || false,
                )
                .unwrap();
            other.restore_states(component.states()).unwrap();
            let mut baseline_thermal = thermal.clone();
            let mut baseline = other
                .into_thermoelastic_system(
                    &baseline_thermal,
                    network,
                    &DESTINATIONS,
                    0.011,
                    live_config(),
                    None,
                    &gate,
                )
                .unwrap();
            let mut owner = component
                .into_thermoelastic_system(
                    &thermal,
                    network,
                    &DESTINATIONS,
                    0.011,
                    live_config(),
                    None,
                    &gate,
                )
                .unwrap();
            for _ in 0..4 {
                owner
                    .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .unwrap();
                baseline
                    .step(&mut baseline_thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .unwrap();
            }
            let states = owner.mechanical().components()[0].states().to_vec();
            let frame = owner.mechanical().last_frame().unwrap().clone();
            let coefficients: Vec<_> = (0..9)
                .map(|i| owner.mechanical().connection_coefficients(i).unwrap())
                .collect();
            let energies: Vec<_> = thermal
                .cells()
                .values()
                .map(|c| c.specific_h_j_kg)
                .collect();
            let clock = thermal.time_s();
            let mut refused = solve();
            refused.maximum_material_trials = 1;
            assert!(
                owner
                    .step(&mut thermal, &[0.0], refused, &[0.0; 2], &gate)
                    .is_err()
            );
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(matches!(
                owner.step(&mut thermal, &[0.0], solve(), &[0.0; 2], &cancelled),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            let force = 100.0 * owner.basis().modes()[0].displacement[3];
            assert!(matches!(
                owner.step(&mut thermal, &[force], solve(), &[0.0; 2], &gate),
                Err(CircuitEnthalpyError::Source(
                    fs_couple::frame_vibration::thermal::FrameThermalSourceError::Material(_)
                ))
            ));
            assert_eq!(owner.mechanical().components()[0].states(), states);
            assert_eq!(owner.mechanical().last_frame(), Some(&frame));
            assert_eq!(thermal.time_s(), clock);
            assert_eq!(
                thermal
                    .cells()
                    .values()
                    .map(|c| c.specific_h_j_kg)
                    .collect::<Vec<_>>(),
                energies
            );
            assert_eq!(
                (0..9)
                    .map(|i| owner.mechanical().connection_coefficients(i).unwrap())
                    .collect::<Vec<_>>(),
                coefficients
            );
            for _ in 0..8 {
                owner
                    .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .unwrap();
                baseline
                    .step(&mut baseline_thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .unwrap();
            }
            assert_eq!(
                owner.mechanical().components()[0].states(),
                baseline.mechanical().components()[0].states()
            );
            assert_eq!(thermal.cells(), baseline_thermal.cells());
        },
    );
    let mut invalid_capacity = live_specs();
    invalid_capacity[0].alpha = 0.1;
    fixture(
        invalid_capacity,
        false,
        false,
        0.0,
        |component, thermal, network| {
            assert!(
                component
                    .into_thermoelastic_system(
                        &thermal,
                        network,
                        &DESTINATIONS,
                        0.05,
                        live_config(),
                        None,
                        &CancelGate::new_clock_free()
                    )
                    .is_err()
            );
            assert_eq!(thermal.time_s(), 0.0);
        },
    );
    let source_exit = specimens().map(|mut spec| {
        spec.modulus_slope = 2.0;
        spec.mechanical_upper = 300.000001;
        spec.cp *= 0.01;
        spec
    });
    fixture(
        source_exit,
        false,
        false,
        0.0,
        |component, mut thermal, network| {
            let gate = CancelGate::new_clock_free();
            let mut owner = component
                .into_thermoelastic_system(
                    &thermal,
                    network,
                    &DESTINATIONS,
                    0.05,
                    live_config(),
                    None,
                    &gate,
                )
                .unwrap();
            let before = owner.mechanical().components()[0].states().to_vec();
            let coefficients: Vec<_> = (0..9)
                .map(|i| owner.mechanical().connection_coefficients(i))
                .collect();
            assert!(
                owner
                    .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .is_err()
            );
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(owner.mechanical().samples_rendered(), 0);
            assert_eq!(owner.mechanical().components()[0].states(), before);
            assert_eq!(
                (0..9)
                    .map(|i| owner.mechanical().connection_coefficients(i))
                    .collect::<Vec<_>>(),
                coefficients
            );
        },
    );
}
