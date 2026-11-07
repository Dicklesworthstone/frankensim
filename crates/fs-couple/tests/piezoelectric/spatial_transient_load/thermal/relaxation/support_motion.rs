//! Independent affine-strain, Maxwell/Debye and finite-cell heat oracles.
//! All material constants are authored fixtures, not empirical calibration.
use super::*;

const MAXWELL_FRACTION: f64 = 0.5;
const MAXWELL_TIME_S: f64 = 0.7e-6;
const RESISTANCE_ALPHA_K_INV: f64 = 0.02;
const PIEZO_CHARGE_C: f64 = 12.0 * AREA;

struct ElectricalReference {
    midpoint_voltage: f64,
    voltage: f64,
    polarization: f64,
    dielectric_heat: f64,
    resistor_heat: f64,
}

// For affine axial strain, Q = C*V + P - e33*A*strain.
// Independently eliminate Qdot = -V/R and tau*Pdot + P = Cp*V.
fn electrical_reference(
    voltage: f64,
    polarization: f64,
    midpoint_strain_rate: f64,
    dt: f64,
    resistance: Option<f64>,
    memory: bool,
) -> ElectricalReference {
    let h = 0.5 * dt;
    let conductance = resistance.map_or(0.0, |r| 1.0 / r);
    let cp = if memory {
        POLE_INCREMENT_F_M * LENGTH
    } else {
        0.0
    };
    let vm = (CAPACITANCE * voltage / h
        + PIEZO_CHARGE_C * midpoint_strain_rate
        + polarization / (POLE_TIME_S + h))
        / (CAPACITANCE / h + conductance + cp / (POLE_TIME_S + h));
    let pm = (POLE_TIME_S * polarization + h * cp * vm) / (POLE_TIME_S + h);
    let next_p = 2.0 * pm - polarization;
    ElectricalReference {
        midpoint_voltage: vm,
        voltage: 2.0 * vm - voltage,
        polarization: next_p,
        dielectric_heat: if memory {
            POLE_TIME_S * (next_p - polarization).powi(2) / (cp * dt)
        } else {
            0.0
        },
        resistor_heat: dt * conductance * vm * vm,
    }
}

fn support_heating(memory: bool, with_resistor: bool) {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (ceramic_card, ceramic_curve, base, loss) =
            ceramic_with_mechanical_memory(false, Some((MAXWELL_FRACTION, MAXWELL_TIME_S)));
        let materials = [base];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|dof| (dof, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let mut body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            field_controls(),
            cx,
        )
        .unwrap();
        if memory {
            body = body
                .with_relaxation(vec![loss], MEMORY_BUDGET, field_controls(), cx)
                .unwrap();
        }
        let (resistor_card, resistor_curve) =
            thermal_material(2.0, 4.0, RESISTANCE_ALPHA_K_INV, 400.0);
        let cards = [ceramic_card, resistor_card];
        let curves = [ceramic_curve, resistor_curve];
        let context = point(300.0);
        let cells = caloric_cells(&cards, &curves, &context);
        let resistors = if with_resistor {
            vec![CellResistor {
                destination: ResistorCell {
                    branch: LOAD.resistor_branch,
                    cell: CELLS[1],
                },
                selection: MaterialPropertySelection::SingleClaimOnly,
                length_m: LENGTH,
                area_m2: AREA,
            }]
        } else {
            vec![]
        };
        let bindings = if with_resistor { vec![LOAD] } else { vec![] };
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let mut thermal = isolated_cells(&curves);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let routes = if memory {
            ceramic_routes(cube.tets.len())
        } else {
            vec![]
        };
        let initial_energy = thermal.energy_j();
        let dt = 0.4 * POLE_TIME_S;
        let h = 0.5 * dt;
        let delta_c = MAXWELL_FRACTION * STIFFNESS[2];
        let cp = POLE_INCREMENT_F_M * LENGTH;
        let mass = CERAMIC_DENSITY_KG_M3 * VOLUME;
        let mut strain = 0.0_f64;
        let mut rate = 0.0_f64;
        let mut voltage = 0.0_f64;
        let mut polarization = 0.0_f64;
        let mut maxwell_stress = 0.0_f64;
        let mut ceramic_heat = 0.0_f64;
        let mut supplied = 0.0_f64;
        let mut returned_work = false;
        for (step, next_rate) in [80.0, 40.0, -60.0, 0.0].into_iter().enumerate() {
            let target = thermal.time_s() + dt;
            let midpoint_rate = 0.5 * (rate + next_rate);
            let next_strain = strain + dt * midpoint_rate;
            let supports = (0..3 * n)
                .map(|dof| {
                    (
                        dof,
                        if dof % 3 == 2 {
                            next_rate * cube.nodes[dof / 3][2]
                        } else {
                            0.0
                        },
                    )
                })
                .collect::<Vec<_>>();
            let before_t = CELLS.map(|cell| thermal.temperature_k(cell).unwrap());
            let before_u = body.state().displacement_m().to_vec();
            let before_v = body.state().velocity_m_s().to_vec();
            let before_time = body.state().time_s();
            let before_memory = body.state().material_memory_energy_j();
            // The finite conductor's endpoint temperature determines R(T).
            // Solve its scalar caloric equation independently by bisection.
            let mut resistor_temperature = before_t[1];
            if with_resistor {
                let (mut low, mut high) = (before_t[1], 400.0);
                for _ in 0..80 {
                    let t = 0.5 * (low + high);
                    let r = RESISTANCE * (1.0 + RESISTANCE_ALPHA_K_INV * (t - 300.0));
                    let reference = electrical_reference(
                        voltage,
                        polarization,
                        midpoint_rate,
                        dt,
                        Some(r),
                        memory,
                    );
                    if reference.resistor_heat > 2.0 * VOLUME * (t - before_t[1]) {
                        low = t;
                    } else {
                        high = t;
                    }
                }
                resistor_temperature = 0.5 * (low + high);
            }
            let resistance = with_resistor.then_some(
                RESISTANCE * (1.0 + RESISTANCE_ALPHA_K_INV * (resistor_temperature - 300.0)),
            );
            let reference =
                electrical_reference(voltage, polarization, midpoint_rate, dt, resistance, memory);
            // Maxwell stress obeys tau*sdot + s = delta_C*tau*strain_dot.
            let sm = if memory {
                (MAXWELL_TIME_S * maxwell_stress + h * delta_c * MAXWELL_TIME_S * midpoint_rate)
                    / (MAXWELL_TIME_S + h)
            } else {
                0.0
            };
            let next_stress = 2.0 * sm - maxwell_stress;
            let mechanical_heat = if memory {
                VOLUME * dt * sm * sm / (delta_c * MAXWELL_TIME_S)
            } else {
                0.0
            };
            let retained_memory = if memory {
                0.5 * reference.polarization.powi(2) / cp
                    + 0.5 * VOLUME * next_stress.powi(2) / delta_c
            } else {
                0.0
            };
            let stored = 0.5 * VOLUME * STIFFNESS[2] * next_strain.powi(2)
                + 0.5 * CAPACITANCE * reference.voltage.powi(2)
                + retained_memory;
            let kinetic = mass * LENGTH.powi(2) * next_rate.powi(2) / 6.0;
            let midpoint_stress = 0.5 * STIFFNESS[2] * (strain + next_strain)
                + sm
                + 12.0 * reference.midpoint_voltage / LENGTH;
            let expected_work = VOLUME * midpoint_stress * (next_strain - strain) + kinetic
                - mass * LENGTH.powi(2) * rate.powi(2) / 6.0;
            let advance = |trace: &[(usize, f64)], solve| {
                if memory {
                    body.step_with_material_thermal_support_motion_candidate(
                        &thermal,
                        target,
                        trace,
                        &prepared,
                        &bindings,
                        &routes,
                        "T",
                        &[Drive::Voltage(0.0), Drive::Current(0.0)],
                        &zero,
                        &charge,
                        field_controls(),
                        solve,
                        &[0.0; 2],
                        cx,
                    )
                } else {
                    body.step_with_thermal_support_motion_candidate(
                        &thermal,
                        target,
                        trace,
                        &prepared,
                        &bindings,
                        &[Drive::Voltage(0.0), Drive::Current(0.0)],
                        &zero,
                        &charge,
                        field_controls(),
                        solve,
                        &[0.0; 2],
                        cx,
                    )
                }
            };
            if step == 0 && memory && with_resistor {
                // Refusal after field evaluation cannot publish a partial heat
                // or support-motion history; the same accepted state can retry.
                assert!(matches!(
                    advance(
                        &supports,
                        CircuitEnthalpySolve {
                            maximum_circuit_trials: 1,
                            ..thermal_controls()
                        }
                    ),
                    Err(CircuitEnthalpyError::CircuitBudget)
                ));
                let mut missing = supports.clone();
                missing.pop();
                assert!(matches!(
                    advance(&missing, thermal_controls()),
                    Err(CircuitEnthalpyError::Source(
                        LoadThermalError::Piezoelectric(_)
                    ))
                ));
            }
            if step == 1 && memory && with_resistor {
                // A stationary call must not silently erase accepted support
                // velocity and discard the mechanical stopping work.
                assert!(matches!(
                    body.step_with_material_thermal_loads_candidate(
                        &thermal,
                        target,
                        &prepared,
                        &bindings,
                        &routes,
                        "T",
                        &[Drive::Voltage(0.0), Drive::Current(0.0)],
                        &zero,
                        &charge,
                        field_controls(),
                        thermal_controls(),
                        &[0.0; 2],
                        cx,
                    ),
                    Err(CircuitEnthalpyError::Source(
                        LoadThermalError::Piezoelectric(_)
                    ))
                ));
            }
            let candidate = advance(&supports, thermal_controls()).unwrap();
            let report = candidate.report();
            let native = &report.electrical.native;
            near(native.support_work_j, expected_work, 2e-16, 4e-8);
            assert_eq!(native.applied_force_work_j, 0.0);
            near(native.mechanical_work_j, expected_work, 2e-16, 4e-8);
            near(
                native.mechanical_material_heat_j,
                mechanical_heat,
                2e-16,
                4e-8,
            );
            near(
                native.dielectric_material_heat_j,
                reference.dielectric_heat,
                2e-16,
                4e-8,
            );
            near(
                report.electrical.loads.joule_heat_j,
                reference.resistor_heat,
                2e-16,
                4e-8,
            );
            near(
                report.deposited_joule_j,
                mechanical_heat + reference.dielectric_heat + reference.resistor_heat,
                2e-16,
                4e-8,
            );
            near(
                report.electrical.loads.electrical_source_work_j,
                0.0,
                2e-16,
                0.0,
            );
            near(report.non_electrical_heat_j, 0.0, 2e-16, 0.0);
            near(report.total_energy_residual_j, 0.0, 2e-16, 0.0);
            assert_eq!(report.electrical.ceramic_destinations, routes);
            let reaction_work: f64 = report
                .electrical
                .midpoint_support_force_n
                .iter()
                .zip(
                    candidate
                        .body()
                        .state()
                        .displacement_m()
                        .iter()
                        .zip(&before_u),
                )
                .map(|(force, (next, old))| {
                    (0..3)
                        .map(|axis| force[axis] * (next[axis] - old[axis]))
                        .sum::<f64>()
                })
                .sum();
            near(reaction_work, expected_work, 2e-16, 4e-8);
            near(
                candidate.body().state().kinetic_energy_j(),
                kinetic,
                2e-16,
                4e-8,
            );
            near(
                candidate.body().state().stored_energy_j(),
                stored,
                2e-16,
                4e-8,
            );
            near(
                candidate.body().state().material_memory_energy_j(),
                retained_memory,
                2e-16,
                4e-8,
            );
            near(
                candidate.body().state().electrodes()[1].voltage_v,
                reference.voltage,
                2e-8,
                4e-8,
            );
            near(
                candidate.body().state().electrodes()[1].charge_c,
                CAPACITANCE * reference.voltage + reference.polarization
                    - PIEZO_CHARGE_C * next_strain,
                2e-18,
                4e-8,
            );
            for (node, coordinates) in cube.nodes.iter().enumerate() {
                near(
                    candidate.body().state().displacement_m()[node][2],
                    next_strain * coordinates[2],
                    1e-18,
                    4e-8,
                );
                assert_eq!(
                    candidate.body().state().velocity_m_s()[node][2],
                    supports[3 * node + 2].1
                );
            }
            ceramic_heat += mechanical_heat + reference.dielectric_heat;
            supplied += expected_work;
            returned_work |= expected_work < 0.0;
            near(
                candidate.thermal().temperature_k(CELLS[0]).unwrap(),
                300.0 + ceramic_heat / (mass * CERAMIC_CP_J_KG_K),
                2e-7,
                0.0,
            );
            near(
                candidate.thermal().temperature_k(CELLS[1]).unwrap(),
                resistor_temperature,
                2e-7,
                0.0,
            );
            near(
                candidate.thermal().energy_j() + candidate.body().state().total_energy_j(),
                initial_energy + supplied,
                8e-16,
                0.0,
            );
            assert_eq!(body.state().time_s(), before_time);
            assert_eq!(body.state().material_memory_energy_j(), before_memory);
            assert_eq!(body.state().displacement_m(), before_u.as_slice());
            assert_eq!(body.state().velocity_m_s(), before_v.as_slice());
            assert_eq!(
                CELLS.map(|cell| thermal.temperature_k(cell).unwrap()),
                before_t
            );
            assert_eq!(candidate.thermal().time_s(), target);
            assert_eq!(candidate.body().state().time_s(), target);
            (thermal, body, _) = candidate.into_parts();
            strain = next_strain;
            rate = next_rate;
            voltage = reference.voltage;
            polarization = reference.polarization;
            maxwell_stress = next_stress;
        }
        assert!(returned_work);
        if memory {
            assert!(thermal.temperature_k(CELLS[0]).unwrap() > 300.0);
        } else {
            assert_eq!(thermal.temperature_k(CELLS[0]).unwrap(), 300.0);
        }
        if with_resistor {
            assert!(thermal.temperature_k(CELLS[1]).unwrap() > 300.0);
        } else {
            assert_eq!(thermal.temperature_k(CELLS[1]).unwrap(), 300.0);
        }
        let accepted_time = thermal.time_s();
        let accepted_energy = thermal.energy_j() + body.state().total_energy_j();
        gate.request();
        assert!(matches!(
            body.step_with_material_thermal_support_motion_candidate(
                &thermal,
                accepted_time + dt,
                &fixed,
                &prepared,
                &bindings,
                &routes,
                "T",
                &[Drive::Voltage(0.0), Drive::Current(0.0)],
                &zero,
                &charge,
                field_controls(),
                thermal_controls(),
                &[0.0; 2],
                cx,
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.time_s(), accepted_time);
        assert_eq!(body.state().time_s(), accepted_time);
        assert_eq!(
            thermal.energy_j() + body.state().total_energy_j(),
            accepted_energy
        );
    });
}

#[test]
fn g1_support_work_heats_a_temperature_dependent_physical_shunt() {
    support_heating(false, true);
}

#[test]
fn g1_support_work_heats_maxwell_and_debye_memory_without_external_resistors() {
    support_heating(true, false);
}

#[test]
fn g1_g4_support_work_couples_separate_ceramic_and_resistor_heat_transactionally() {
    support_heating(true, true);
}

#[test]
fn g3_finite_ceramic_and_resistor_heat_ignore_unchanged_rigid_kinetic_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (ceramic_card, ceramic_curve, base, memory) =
            ceramic_with_mechanical_memory(false, Some((MAXWELL_FRACTION, MAXWELL_TIME_S)));
        let materials = [base];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let fixed = (0..3 * n).map(|dof| (dof, 0.0)).collect::<Vec<_>>();
        let original = Body::new(&model, &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero, &zero, &charge, field_controls(), cx).unwrap()
            .with_relaxation(vec![memory], MEMORY_BUDGET, field_controls(), cx).unwrap();
        let (resistor_card, resistor_curve) =
            thermal_material(2.0, 4.0, RESISTANCE_ALPHA_K_INV, 400.0);
        let cards = [ceramic_card, resistor_card];
        let curves = [ceramic_curve, resistor_curve];
        let context = point(300.0);
        let cells = caloric_cells(&cards, &curves, &context);
        let resistors = [CellResistor {
            destination: ResistorCell { branch: LOAD.resistor_branch, cell: CELLS[1] },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: LENGTH, area_m2: AREA,
        }];
        let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
        let routes = ceramic_routes(cube.tets.len());
        let dt = 2.0_f64.powi(-20);
        let rate = 0.25;
        let expected_kinetic = CERAMIC_DENSITY_KG_M3 * LENGTH.powi(5) * rate * rate / 6.0;
        let mut reference = Vec::new();
        for boost in [0.0, 2.0_f64.powi(16)] {
            let mut thermal = isolated_cells(&curves);
            let prepared = network.prepare(&thermal, || false).unwrap();
            let mut body = original.clone();
            for (index, next_rate) in [0.0, rate, 0.0].into_iter().enumerate() {
                let supports = (0..3 * n).map(|dof| (dof, match dof % 3 {
                    0 => boost,
                    2 => next_rate * cube.nodes[dof / 3][2],
                    _ => 0.0,
                })).collect::<Vec<_>>();
                let candidate = body.step_with_material_thermal_support_motion_candidate(
                    &thermal, thermal.time_s() + dt, &supports, &prepared, &[LOAD],
                    &routes, "T", &[Drive::Voltage(0.0), Drive::Current(0.0)],
                    &zero, &charge, field_controls(), thermal_controls(), &[0.0; 2], cx,
                ).unwrap();
                let report = candidate.report();
                if index > 0 {
                    let native = &report.electrical.native;
                    near(native.kinetic_energy_change_j,
                        if index == 1 { expected_kinetic } else { -expected_kinetic },
                        1e-24, 1e-10);
                    assert!(report.deposited_joule_j > 0.0);
                    assert!(native.material_heat_j > 0.0);
                    assert!(report.electrical.loads.joule_heat_j > 0.0);
                    assert!(report.total_energy_residual_j.abs() < 1e-16);
                    assert_eq!(report.electrical.ceramic_destinations, routes);
                    let observed = [native.support_work_j, native.material_heat_j,
                        report.electrical.loads.joule_heat_j, report.deposited_joule_j,
                        candidate.body().state().electrodes()[1].voltage_v,
                        candidate.thermal().temperature_k(CELLS[0]).unwrap(),
                        candidate.thermal().temperature_k(CELLS[1]).unwrap()];
                    if boost == 0.0 {
                        reference.push(observed);
                    } else {
                        for (j, (actual, expected)) in observed.into_iter()
                            .zip(reference[index - 1]).enumerate() {
                            near(actual, expected, if j >= 5 { 1e-10 } else { 1e-23 }, 1e-9);
                        }
                    }
                }
                (thermal, body, _) = candidate.into_parts();
            }
            assert!(thermal.temperature_k(CELLS[0]).unwrap() > 300.0);
            assert!(thermal.temperature_k(CELLS[1]).unwrap() > 300.0);
            assert!(body.state().velocity_m_s().iter().all(|v| v[0] == boost));
        }
    });
}
