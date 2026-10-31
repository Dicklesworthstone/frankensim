//! Manufactured fixed-domain conductor motion through source-normalized coils.
//! The exact center-hat mode separates diffusion, transport drag and back-EMF.
use super::*;
use fs_em::{
    planar::eddy::{moving::MovingEddyControls, nonlinear::moving::NonlinearMovingEddyControls},
    winding::nonlinear::NonlinearWindingControls,
};

fn moving_controls() -> NonlinearWindingControls<NonlinearMovingEddyControls> {
    let base = controls();
    NonlinearWindingControls {
        field: NonlinearMovingEddyControls {
            moving: MovingEddyControls {
                relative_tolerance: 1e-12,
                absolute_residual_a_m: 1e-13,
                restart: 8,
                max_cycles: 50,
                max_krylov_scalars: 10_000,
                max_cell_magnetic_peclet: 2.0,
                energy_tolerance_j: 1e-9,
                current_balance_tolerance_a: 1e-9,
            },
            absolute_residual_tolerance_a_m: 1e-11,
            relative_residual_tolerance: 1e-11,
            max_newton_iterations: 30,
            max_backtracks: 30,
        },
        maximum_field_trials: base.maximum_field_trials,
        maximum_newton_iterations: base.maximum_newton_iterations,
        difference_step_a: base.difference_step_a,
        terminal_voltage_tolerance_v: base.terminal_voltage_tolerance_v,
        power_tolerance_w: base.power_tolerance_w,
        energy_tolerance_j: base.energy_tolerance_j,
        newton: base.newton,
    }
}
fn velocities(problem: &PlanarNonlinearMagnetostatics, speed: f64) -> Vec<[f64; 2]> {
    // Constant through-flow within each homogeneous conductor. Outgoing and
    // return material move oppositely; their fixed boundary data supply the
    // invariant Eulerian occupancy. No finite conductor interface is translated.
    problem
        .mesh()
        .vertices()
        .iter()
        .map(|p| {
            [
                if p[0] <= 1.0 {
                    speed
                } else if p[0] >= 2.0 {
                    -speed
                } else {
                    0.0
                },
                0.0,
            ]
        })
        .collect()
}
fn uniform_sources(card: &MaterialCard, temperature: f64) -> Vec<EddyConductor> {
    let copper =
        common::material::conductor(card, temperature, ElectricalPropertyBasis::Resistivity);
    sources(card, None)
        .into_iter()
        .map(|mut source| {
            source.conductor = copper.clone();
            source
        })
        .collect()
}

#[test]
fn zero_velocity_retains_the_same_nonlinear_solid_diffusion_and_heat() {
    for saturating in [false, true] {
        let card = common::source(saturating, SIGMA);
        let problem = common::material::problem(&card, mesh());
        let original = driven(&problem, &card, 0.0, 0.1);
        let stationary = original
            .step_to_candidate(
                0.2,
                sources(&card, None),
                vec![],
                &[2.0],
                controls(),
                &mut || false,
            )
            .unwrap();
        let moving = original
            .step_with_eulerian_velocity_candidate(
                0.2,
                &velocities(&problem, 0.0),
                &[2.0],
                moving_controls(),
                &mut || false,
            )
            .unwrap();
        near(
            moving.step().terminal_voltage_v[0],
            stationary.step().terminal_voltage_v[0],
        );
        near(
            moving.step().inductive_voltage_v[0],
            stationary.step().inductive_voltage_v[0],
        );
        near(moving.step().motional_voltage_v[0], 0.0);
        near(moving.step().heat_j, stationary.step().heat_j);
        near(
            moving.step().numerical_dissipation_j,
            stationary.step().numerical_dissipation_j,
        );
        near(moving.step().mechanical_work_j, 0.0);
        near(
            moving.step().field.uniform_drive_work_j,
            stationary.step().field.uniform_drive_work_j,
        );
        near(moving.step().field.motional_drive_work_j, 0.0);
        for (&a, &b) in moving
            .circuit()
            .field()
            .potential_wb_m()
            .iter()
            .zip(stationary.circuit().field().potential_wb_m())
        {
            near(a, b);
        }
        for (a, b) in moving
            .step()
            .field
            .cells
            .iter()
            .zip(&stationary.step().field.cells)
        {
            near(a.joule.heat_j, b.heat_j);
            for k in 0..3 {
                near(a.joule.current_density_a_m2[k], b.current_density_a_m2[k]);
                near(a.laboratory_electric_field_v_m[k], b.electric_field_v_m[k]);
            }
        }
    }
}

#[test]
fn moving_solid_back_emf_survives_zero_current_and_closes_motor_generator_work() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let bias = 0.2;
    let end = 0.04;
    let dt = 0.1;
    let dc_resistance = 72.0 * DEPTH / SIGMA;
    for speed in [-0.08_f64, 0.08] {
        let original = driven(&problem, &card, bias, end);
        for current in [0.0, 0.4 * speed.signum(), -0.01 * speed.signum()] {
            let a = current / (2.0 + SIGMA / (36.0 * dt));
            let emf = 12.0 * DEPTH * speed * bias;
            let inductive = (4.0 * DEPTH * a + end * current) / dt;
            let voltage = dc_resistance * current + inductive + emf;
            let dc_heat = dc_resistance * current * current * dt;
            let diffusion_heat = SIGMA * DEPTH * a * a / (9.0 * dt);
            // Integral (dN/dx)^2=2 on each unit square. The central Galerkin
            // convection term is zero because integral N*dN/dx=0, but local
            // transport current, drag and heat remain nonzero and reciprocal.
            let drag_heat = 4.0 * SIGMA * DEPTH * dt * speed * speed * a * a;
            let mechanical = current * emf * dt - drag_heat;
            let candidate = original
                .step_with_eulerian_velocity_candidate(
                    dt,
                    &velocities(&problem, speed),
                    &[current],
                    moving_controls(),
                    &mut || false,
                )
                .unwrap();
            let step = candidate.step();
            near(candidate.circuit().currents_a()[0], current);
            near(
                candidate.circuit().field().potential_wb_m()[8],
                0.5 * bias + a,
            );
            near(
                candidate.circuit().field().potential_wb_m()[10],
                2.5 * bias - a,
            );
            near(step.inductive_voltage_v[0], inductive);
            near(step.motional_voltage_v[0], emf);
            near(step.terminal_voltage_v[0], voltage);
            near(step.field.uniform_drive_work_j, dc_heat);
            near(step.field.motional_drive_work_j, current * emf * dt);
            near(step.heat_j, dc_heat + diffusion_heat + drag_heat);
            near(step.mechanical_work_j, mechanical);
            near(step.supplied_work_j, dt * current * voltage);
            near(step.stored_before_j, 1.5 * DEPTH * bias * bias);
            near(
                step.stored_after_j,
                1.5 * DEPTH * bias * bias + 4.0 * DEPTH * a * a + 0.5 * end * current * current,
            );
            near(
                step.numerical_dissipation_j,
                4.0 * DEPTH * a * a + 0.5 * end * current * current,
            );
            assert!(step.energy_residual_j.abs() < moving_controls().energy_tolerance_j);
            for &(region, mean) in &step.field.region_motional_electric_field_v_m {
                near(
                    mean,
                    if region == 10 {
                        -speed * bias
                    } else {
                        speed * bias
                    },
                );
            }
            for cell in &step.field.cells {
                let source = &cell.joule.source;
                let polarity = if source.region == 10 { 1.0 } else { -1.0 };
                let nx = [0.0, -2.0, 0.0, 2.0][source.triangle % 4];
                let ny = [2.0, 0.0, -2.0, 0.0][source.triangle % 4];
                let gradient = [bias + polarity * a * nx, polarity * a * ny];
                let offset = polarity * speed * bias + speed * a * nx;
                for (k, &vertex) in problem.mesh().faces()[source.triangle].iter().enumerate() {
                    let n = if vertex >= 8 { 1.0 } else { 0.0 };
                    let electric = polarity * (6.0 * current / SIGMA - a * (n - 1.0 / 3.0) / dt)
                        - speed * a * nx;
                    near(cell.joule.electric_field_v_m[k], electric);
                    near(cell.joule.current_density_a_m2[k], SIGMA * electric);
                    near(cell.laboratory_electric_field_v_m[k], electric + offset);
                    near(cell.motional_electric_field_v_m[k], -offset);
                    for axis in 0..2 {
                        near(
                            cell.lorentz_force_density_n_m3[k][axis],
                            SIGMA * electric * gradient[axis],
                        );
                    }
                }
                near(
                    cell.mechanical_work_j,
                    dt * cell.joule.volume_m3
                        * offset
                        * cell.joule.current_density_a_m2.iter().sum::<f64>()
                        / 3.0,
                );
            }
            if current == 0.0 {
                assert!(step.motional_voltage_v[0].abs() > 0.05);
                near(step.heat_j, 0.0);
                near(step.mechanical_work_j, 0.0);
            } else if current * speed < 0.0 {
                assert!(step.supplied_work_j < 0.0 && step.mechanical_work_j < 0.0);
            } else {
                assert!(step.mechanical_work_j > 0.0);
            }
        }
    }
}

#[test]
fn moving_solid_source_refresh_retains_diffusion_and_recomputes_actual_heat() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let bias = 0.2;
    let speed = 0.08;
    let original = driven(&problem, &card, bias, 0.0);
    let first = original
        .step_with_eulerian_velocity_candidate(
            0.1,
            &velocities(&problem, speed),
            &[0.4],
            moving_controls(),
            &mut || false,
        )
        .unwrap();
    let accepted = first.into_parts().0;
    let before_a = 0.4 / (2.0 + SIGMA / 3.6);
    let sigma = 12.0_f64; // Original rho(T) at 400 K, with fixed H(B).
    let dt = 0.15;
    let current = 0.2;
    let q = sigma / (36.0 * dt);
    let a = (current + q * before_a) / (2.0 + q);
    let da = a - before_a;
    let dc = 72.0 * DEPTH / sigma;
    let emf = 12.0 * DEPTH * speed * bias;
    let drag_heat = 4.0 * sigma * DEPTH * dt * speed * speed * a * a;
    let candidate = accepted
        .step_with_eulerian_velocity_and_sources_candidate(
            0.25,
            &velocities(&problem, speed),
            uniform_sources(&card, 400.0),
            vec![],
            &[current],
            moving_controls(),
            &mut || false,
        )
        .unwrap();
    let step = candidate.step();
    near(
        candidate.circuit().field().potential_wb_m()[8],
        0.5 * bias + a,
    );
    near(candidate.circuit().resistance_ohm()[0], dc);
    near(step.inductive_voltage_v[0], 4.0 * DEPTH * da / dt);
    near(step.motional_voltage_v[0], emf);
    near(
        step.terminal_voltage_v[0],
        dc * current + 4.0 * DEPTH * da / dt + emf,
    );
    near(
        step.heat_j,
        dc * current * current * dt + sigma * DEPTH * da * da / (9.0 * dt) + drag_heat,
    );
    near(step.mechanical_work_j, current * emf * dt - drag_heat);
    assert!((a - current / (2.0 + q)).abs() > 0.01);
    for cell in &step.field.cells {
        near(
            cell.joule.source.conductor.resistivity_ohm_m(),
            sigma.recip(),
        );
    }
    for (&actual, &expected) in candidate.circuit().profiles_per_m2()[0]
        .iter()
        .zip(&profiles()[0])
    {
        near(actual, expected);
    }
    near(accepted.time_s(), 0.1);
    near(accepted.resistance_ohm()[0], 72.0 * DEPTH / SIGMA);
    near(accepted.field().potential_wb_m()[8], 0.5 * bias + before_a);
}

#[test]
fn motion_drives_a_loaded_native_circuit_without_external_electrical_supply() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let bias = 0.2;
    let speed = 0.08;
    let dt = 0.1;
    let load = 0.7;
    let end = 0.04;
    let original = NonlinearWindingCircuit::new_with_solid_regions(
        initial(
            &problem,
            sources(&card, None),
            &[4, 5, 6, 7],
            &boundary(&problem, bias),
        ),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 777.0 }),
                (2, 0, Branch::Resistor { ohms: load }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles(),
            end_inductance_h: vec![end],
        },
        vec![],
        bindings(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap();
    let q = SIGMA / (36.0 * dt);
    let emf = 12.0 * DEPTH * speed * bias;
    let current = -emf / (72.0 * DEPTH / SIGMA + load + end / dt + 4.0 * DEPTH / (dt * (2.0 + q)));
    let candidate = original
        .step_with_eulerian_velocity_candidate(
            dt,
            &velocities(&problem, speed),
            &[0.0],
            moving_controls(),
            &mut || false,
        )
        .unwrap();
    let step = candidate.step();
    near(candidate.circuit().currents_a()[0], current);
    near(step.terminal_voltage_v[0], -load * current);
    near(
        step.graph_resistors[0].heat_j,
        load * current * current * dt,
    );
    near(step.supplied_work_j, 0.0);
    assert!(step.mechanical_work_j < 0.0 && step.graph_resistors[0].heat_j > 0.0);
    near(
        step.stored_after_j - step.stored_before_j
            + step.heat_j
            + step.numerical_dissipation_j
            + step.mechanical_work_j,
        0.0,
    );
    let mut bad = moving_controls();
    bad.maximum_field_trials = 1;
    assert!(matches!(
        original.step_with_eulerian_velocity_candidate(
            dt,
            &velocities(&problem, speed),
            &[0.0],
            bad,
            &mut || false,
        ),
        Err(NonlinearWindingError::TrialBudget)
    ));
    assert!(
        original
            .step_with_eulerian_velocity_candidate(
                dt,
                &velocities(&problem, 100.0),
                &[0.0],
                moving_controls(),
                &mut || false,
            )
            .is_err()
    );
    let mut foreign = sources(&card, None);
    foreign[0].region = 99;
    assert!(
        original
            .step_with_eulerian_velocity_and_sources_candidate(
                dt,
                &velocities(&problem, speed),
                foreign,
                vec![],
                &[0.0],
                moving_controls(),
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        original.step_with_eulerian_velocity_candidate(
            dt,
            &velocities(&problem, speed),
            &[0.0],
            moving_controls(),
            &mut || true,
        ),
        Err(NonlinearWindingError::Cancelled)
    ));
    near(original.currents_a()[0], 0.0);
    near(original.time_s(), 0.0);
}
