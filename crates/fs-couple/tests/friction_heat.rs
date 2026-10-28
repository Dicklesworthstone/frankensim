//! Independent finite-volume references for general friction heat coupling.
//! Synthetic caloric and conductivity data are analytical fixtures, not a
//! material calibration. Every mechanics evaluation uses the real dry law.
use fs_couple::friction_heat::*;
use fs_la::prepared_lu::LuWorkspace;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use std::collections::BTreeMap;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static SENSIBLE: Sensible = Sensible;
const CONTACT: u64 = 0;
fn tolerances() -> FrictionHeatTolerances {
    FrictionHeatTolerances {
        cell_j: 1e-10,
        thermal_j: 1e-10,
        interface_j: 1e-10,
        total_j: 1e-10,
    }
}
fn network(
    n: usize,
    lengths: [f64; 2],
    k: [f64; 2],
    capacity_per_volume: [f64; 2],
    cooling: [f64; 2],
    initial_rise: [f64; 2],
) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for side in 0..2 {
        let dx = lengths[side] / n as f64;
        for i in 0..n {
            let id = CellId((side * n + i + 1) as u64);
            let surface = 2.0 + 4.0 * dx;
            let covered = if n == 1 { 1.0 } else { 2.0 };
            let (start, end) = if side == 0 {
                (-((i + 1) as f64) * dx, -(i as f64) * dx)
            } else {
                (i as f64 * dx, (i + 1) as f64 * dx)
            };
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: side as u64,
                    material_id: MaterialId(1),
                    mass_kg: capacity_per_volume[side] * dx,
                    specific_h_j_kg: initial_rise[side],
                    geometry: CellGeometry {
                        start_m: [start, 0.0, 0.0],
                        end_m: [end, 0.0, 0.0],
                        volume_m3: dx,
                        surface_area_m2: surface,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: if i == 0 {
                            cooling[side] / (surface - covered)
                        } else {
                            0.0
                        },
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
    }
    for side in 0..2 {
        for i in 0..n - 1 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: (side * n + i + 1) as u64,
                    a: CellId((side * n + i + 1) as u64),
                    b: CellId((side * n + i + 2) as u64),
                    area_m2: 1.0,
                    conductance_w_k: k[side] / (lengths[side] / n as f64),
                }),
            });
        }
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: CONTACT,
            a: CellId(1),
            b: CellId((n + 1) as u64),
            area_m2: 1.0,
            conductance_w_k: series_conductance_w_k(
                1.0,
                [lengths[0] / (2 * n) as f64, lengths[1] / (2 * n) as f64],
                k,
            )
            .unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &SENSIBLE as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2 * n,
            max_links: 2 * n,
            max_events: 4 * n + 4,
            max_substeps: 10000,
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
    .0
}
fn input(dt: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 10.0,
        slip_increment_m: [dt, 0.0],
        dt_s: dt,
        contact_temperature_k: None,
        nominal_pressure_pa: Some(10.0),
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
    }
}
fn source<P>(heat: f64, k: [f64; 2], payload: P) -> FrictionHeatInput<P> {
    FrictionHeatInput {
        dissipated_j: heat,
        conductivity_w_m_k: k,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload,
        },
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.16e} != {b:.16e}");
}

// Independent backward-Euler matrix in TEMPERATURE coordinates. The production
// owner works in joules; this uses actual fs-la LU, not a duplicate heat step.
fn factor(thermal: &DepositionSimulation<'_>, dt: f64, half_g: [f64; 2], beta: f64) -> LuWorkspace {
    let count = thermal.cells().len();
    let n = count + 1;
    let mut matrix = vec![0.0; n * n];
    for (i, cell) in thermal.cells().values().enumerate() {
        matrix[i * n + i] = cell.mass_kg / dt
            + cell.boundary.htc_w_m2_k * thermal.exposed_area_m2(cell.id).unwrap();
    }
    for link in thermal.links().values().filter(|link| link.id != CONTACT) {
        let a = link.a.0 as usize - 1;
        let b = link.b.0 as usize - 1;
        let g = link.conductance_w_k;
        matrix[a * n + a] += g;
        matrix[b * n + b] += g;
        matrix[a * n + b] -= g;
        matrix[b * n + a] -= g;
    }
    let ports = [0, count / 2];
    for side in 0..2 {
        let i = ports[side];
        let g = half_g[side];
        matrix[i * n + i] += g;
        matrix[i * n + count] -= g;
        matrix[count * n + i] -= g;
        matrix[count * n + count] += g;
    }
    matrix[count * n + count] += beta;
    let mut lu = LuWorkspace::new(n).unwrap();
    lu.factor(&matrix).unwrap();
    lu
}
fn endpoint(
    thermal: &DepositionSimulation<'_>,
    lu: &mut LuWorkspace,
    dt: f64,
    p0: f64,
) -> (Vec<f64>, f64) {
    let mut rhs: Vec<_> = thermal
        .cells()
        .values()
        .map(|cell| cell.mass_kg / dt * cell.specific_h_j_kg)
        .collect();
    rhs.push(p0);
    let mut temperatures = vec![0.0; rhs.len()];
    lu.solve_factored_into(&rhs, &mut temperatures).unwrap();
    let energy = thermal
        .cells()
        .values()
        .zip(&temperatures)
        .map(|(cell, &t)| cell.mass_kg * (t - cell.specific_h_j_kg))
        .collect();
    (energy, 300.0 + temperatures.last().unwrap())
}

#[test]
fn g0_friction_work_enters_existing_network_once_with_unequal_half_cells() {
    let old = network(1, [1.0, 3.0], [2.0, 9.0], [2.0, 2.0], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [2.0, 9.0]).unwrap();
    assert_eq!(route.distances_m(), [0.5, 1.5]);
    let dt = 0.2;
    let mut lu = factor(&old, dt, [4.0, 6.0], 0.0);
    let (energy, t) = endpoint(&old, &mut lu, dt, 10.0);
    let law = FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 1.0,
    };
    let accepted = law.initial_state();
    let step = FrictionHeatStep::new(&old, dt, &route, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &energy,
            t,
            |_, temperature| {
                let friction =
                    tangential_heat_candidate(&law, &accepted, input(dt), temperature, 1e-12)
                        .unwrap();
                Ok::<_, DepositionError>(source(friction.dissipated_j, [2.0, 9.0], friction))
            },
            || Ok(()),
        )
        .unwrap();
    assert!(trial.residual_j().iter().all(|r| r.abs() < 1e-12));
    let (next, friction, report) = trial.accept(tolerances(), || Ok(())).unwrap();
    // Effective endpoint admittances include each ORIGINAL finite heat capacity.
    let a = 4.0 * 2.0 / (2.0 + dt * 4.0);
    let b = 6.0 * 6.0 / (6.0 + dt * 6.0);
    close(report.body_heat_j[0], 2.0 * a / (a + b), 1e-12);
    close(
        report.body_heat_j.iter().sum(),
        friction.dissipated_j,
        1e-12,
    );
    close(next.energy_j() - old.energy_j(), 2.0, 1e-12);
    assert_eq!(next.contact_patch(CONTACT), old.contact_patch(CONTACT));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(accepted.slip_m(), [0.0; 2]);
    assert_eq!(friction.state.contact_temperature_k(), Some(t));
    assert_eq!(next.cells()[&CellId(1)].boundary.source_w, 0.0);
}

#[test]
fn g0_zero_slip_conduction_is_not_duplicated_or_labelled_friction_heat() {
    let old = network(1, [1.0, 3.0], [2.0, 9.0], [2.0, 2.0], [0.0; 2], [10.0, 0.0]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [2.0, 9.0]).unwrap();
    let dt = 0.2;
    let mut lu = factor(&old, dt, [4.0, 6.0], 0.0);
    let (energy, t) = endpoint(&old, &mut lu, dt, 0.0);
    let step = FrictionHeatStep::new(&old, dt, &route, || Ok(())).unwrap();
    let (next, (), report) = step
        .evaluate(
            &energy,
            t,
            |_, _| Ok::<_, DepositionError>(source(0.0, [2.0, 9.0], ())),
            || Ok(()),
        )
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    let exchange = 2.4 * 10.0 * dt / (1.0 + dt * 2.4 * (0.5 + 1.0 / 6.0));
    close(report.body_heat_j[0], -exchange, 1e-12);
    close(report.body_heat_j[1], exchange, 1e-12);
    close(next.energy_j(), old.energy_j(), 1e-12);
    assert_eq!(report.dissipated_j, 0.0);
}

#[test]
fn g1_decreasing_mu_feedback_reaches_the_analytic_steady_temperature() {
    let mut thermal = network(1, [1.0; 2], [1.0, 2.0], [1.0; 2], [3.0, 5.0], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0, 2.0]).unwrap();
    let dt = 0.25;
    let mut lu = factor(&thermal, dt, [2.0, 4.0], 0.1);
    let cold = FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.8,
    };
    let mut history = cold.initial_state();
    let mut last_temperature = 300.0;
    let mut last_force = 0.0;
    for _ in 0..80 {
        let (energy, t) = endpoint(&thermal, &mut lu, dt, 8.0);
        let step =
            FrictionHeatStep::new(&thermal, thermal.time_s() + dt, &route, || Ok(())).unwrap();
        let trial = step
            .evaluate(
                &energy,
                t,
                |_, temperature| {
                    // Declared analytical material law; genuine Coulomb return map at
                    // this trial's interface temperature, with zero reversible storage.
                    let mu = 0.8 - 0.01 * (temperature - 300.0);
                    let law = FrictionLaw::Coulomb {
                        static_mu: mu,
                        kinetic_mu: mu,
                    };
                    let friction =
                        tangential_heat_candidate(&law, &history, input(dt), temperature, 1e-12)
                            .unwrap();
                    Ok::<_, DepositionError>(source(friction.dissipated_j, [1.0, 2.0], friction))
                },
                || Ok(()),
            )
            .unwrap();
        let (next, friction, report) = trial.accept(tolerances(), || Ok(())).unwrap();
        assert!(report.total_energy_residual_j.abs() < 1e-12);
        last_temperature = report.interface_temperature_k;
        last_force = friction.traction_n[0];
        history = friction.state;
        thermal = next;
    }
    let effective = 2.0 * 3.0 / (2.0 + 3.0) + 4.0 * 5.0 / (4.0 + 5.0);
    let expected = 300.0 + 8.0 / (effective + 0.1);
    close(last_temperature, expected, 1e-9);
    assert!(last_force < 8.0);
    close(last_force, 8.0 - 0.1 * (expected - 300.0), 1e-9);
    assert_eq!(history.contact_temperature_k(), Some(last_temperature));
}

#[test]
fn g0_refusal_publishes_neither_thermal_nor_tangential_state() {
    let old = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let law = FrictionLaw::Coulomb {
        static_mu: 0.5,
        kinetic_mu: 0.4,
    };
    let history = law.initial_state();
    let original = history.clone();
    let original_cells = old.cells().clone();
    let step = FrictionHeatStep::new(&old, 1.0, &route, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &[0.0; 2],
            301.0,
            |_, temperature| {
                let friction =
                    tangential_heat_candidate(&law, &history, input(1.0), temperature, 1e-12)
                        .unwrap();
                Ok::<_, DepositionError>(source(friction.dissipated_j, [1.0; 2], friction))
            },
            || Ok(()),
        )
        .unwrap();
    assert!(trial.accept(tolerances(), || Ok(())).is_err());
    let refused = step.evaluate::<(), DepositionError>(
        &[0.0; 2],
        300.0,
        |_, _| Err(DepositionError::MaterialDomain),
        || Ok(()),
    );
    assert!(matches!(refused, Err(DepositionError::MaterialDomain)));
    let trial = step
        .evaluate(
            &[0.0; 2],
            300.0,
            |_, _| Ok::<_, DepositionError>(source(0.0, [1.0; 2], ())),
            || Ok(()),
        )
        .unwrap();
    let mut checks = 0;
    assert!(matches!(
        trial.accept(tolerances(), || {
            checks += 1;
            if checks == 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.cells(), &original_cells);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(history, original);
}

#[test]
fn g0_material_geometry_and_exclusive_heat_ownership_are_admitted() {
    let old = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    assert!(FrictionHeatRoute::from_link(&old, CONTACT, [2.0; 2]).is_err());
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let step = FrictionHeatStep::new(&old, 1.0, &route, || Ok(())).unwrap();
    for heat in [-1.0, f64::NAN] {
        assert!(
            step.evaluate(
                &[0.0; 2],
                300.0,
                |_, _| Ok::<_, DepositionError>(source(heat, [1.0; 2], ())),
                || Ok(())
            )
            .is_err()
        );
    }
    assert!(
        step.evaluate(
            &[0.0; 2],
            300.0,
            |_, _| {
                let mut s = source(0.0, [1.0; 2], ());
                s.transport.conductances_w_k.insert(CONTACT, 1.0);
                Ok::<_, DepositionError>(s)
            },
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn g1_stationary_partition_and_interface_temperature_converge_on_three_grids() {
    // Constant surface generation between two initially isothermal half spaces:
    // share_A=e_A/(e_A+e_B), dT=2 P sqrt(t)/(sqrt(pi)*(e_A+e_B)).
    // Here A=1, rho*cp=1, k=[1,4], so e=[1,2]. The finite outer
    // boundaries are far beyond the diffusion depth over this interval.
    let mut errors = Vec::new();
    let dt = 0.25 / 512.0;
    let law = FrictionLaw::Coulomb {
        static_mu: 0.3,
        kinetic_mu: 0.3,
    };
    for n in [16, 32, 64] {
        let mut thermal = network(n, [4.0; 2], [1.0, 4.0], [1.0; 2], [0.0; 2], [0.0; 2]);
        let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0, 4.0]).unwrap();
        let half_g = [2.0 * n as f64 / 4.0, 8.0 * n as f64 / 4.0];
        let mut lu = factor(&thermal, dt, half_g, 0.0);
        let mut history = law.initial_state();
        let mut received = [0.0; 2];
        let mut last_t = 300.0;
        for _ in 0..512 {
            let (energy, t) = endpoint(&thermal, &mut lu, dt, 3.0);
            let step =
                FrictionHeatStep::new(&thermal, thermal.time_s() + dt, &route, || Ok(())).unwrap();
            let trial = step
                .evaluate(
                    &energy,
                    t,
                    |_, temperature| {
                        let friction = tangential_heat_candidate(
                            &law,
                            &history,
                            input(dt),
                            temperature,
                            1e-12,
                        )
                        .unwrap();
                        Ok::<_, DepositionError>(source(
                            friction.dissipated_j,
                            [1.0, 4.0],
                            friction,
                        ))
                    },
                    || Ok(()),
                )
                .unwrap();
            let (next, friction, report) = trial.accept(tolerances(), || Ok(())).unwrap();
            for side in 0..2 {
                received[side] += report.body_heat_j[side];
            }
            last_t = report.interface_temperature_k;
            history = friction.state;
            thermal = next;
        }
        close(received.iter().sum(), 0.75, 1e-10);
        close(thermal.energy_j(), 0.75, 1e-10);
        let partition_error = (received[0] / 0.75 - 1.0 / 3.0).abs();
        let temperature_error = ((last_t - 300.0) - 1.0 / core::f64::consts::PI.sqrt()).abs();
        println!(
            "n={n} partition_error={partition_error:e} temperature_error={temperature_error:e}"
        );
        errors.push((partition_error, temperature_error));
    }
    for pair in errors.windows(2) {
        assert!(pair[1].0 < pair[0].0);
        assert!(pair[1].1 < pair[0].1);
    }
    assert!(errors[2].0 / (1.0 / 3.0) < 0.01);
    assert!(errors[2].1 / (1.0 / core::f64::consts::PI.sqrt()) < 0.02);
}
