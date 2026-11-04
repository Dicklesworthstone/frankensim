//! G1/G3/G4 dynamic checks on the original spatial tensor/mesh operators.
use super::*;
use fs_couple::piezoelectric::spatial::transient::TetPiezoelectricTransient as Body;

fn zero_source(n: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
    (vec![[0.0; 3]; n], vec![0.0; n])
}

#[test]
fn g1_open_electrode_oscillator_converges_second_order_to_exact_sinusoid() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed = (0..3 * n)
            .filter(|&i| i != dof)
            .map(|i| (i, 0.0))
            .collect::<Vec<_>>();
        // One unconstrained mechanical DOF and one floating electrode give an
        // exact semidiscrete oscillator. Eliminate its scalar charge equation
        // analytically, then compare time integration against sin/cos, not
        // against another midpoint recurrence or its conserved energy.
        let diagonal = |matrix: &fs_sparse::Csr| {
            let (cols, values) = matrix.row(dof);
            values[cols.iter().position(|&i| i == dof).unwrap()]
        };
        let mass = diagonal(&model.solid().mass);
        let stiffness = diagonal(&model.solid().stiffness);
        let (cols, values) = model.coupling().row(dof);
        let coupling = cols
            .iter()
            .zip(values)
            .filter(|(col, _)| cube.top.contains(col))
            .map(|(_, value)| value)
            .sum::<f64>();
        let capacitance = cube
            .top
            .iter()
            .map(|&row| {
                let (cols, values) = model.dielectric().stiffness().row(row);
                cols.iter()
                    .zip(values)
                    .filter(|(col, _)| cube.top.contains(col))
                    .map(|(_, value)| value)
                    .sum::<f64>()
            })
            .sum::<f64>();
        assert!(coupling.abs() > 0.0);
        let omega = ((stiffness + coupling * coupling / capacitance) / mass).sqrt();
        let duration = 0.7 / omega;
        let amplitude = 1e-9;
        let mut initial_u = zero.clone();
        initial_u[node][2] = amplitude;
        let conditions = [Condition::Voltage(0.0), Condition::Charge(0.0)];
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let mut body = Body::new(
                &model,
                &fixed,
                &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &initial_u,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            for _ in 0..steps {
                body = body
                    .step_candidate(
                        duration / steps as f64,
                        &conditions,
                        &zero,
                        &charge,
                        controls(),
                        cx,
                    )
                    .unwrap()
                    .into_parts()
                    .0;
            }
            let displacement_error =
                body.state().displacement_m()[node][2] / amplitude - 0.7_f64.cos();
            let velocity_error =
                body.state().velocity_m_s()[node][2] / (amplitude * omega) + 0.7_f64.sin();
            errors.push(displacement_error.hypot(velocity_error));
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{errors:?}"
            );
        }
        assert!(errors[2] < 4e-5, "{errors:?}");
    });
}

#[test]
fn g1_free_body_has_consistent_inertia_and_exact_uniform_acceleration() {
    with_cx(|cx, _| {
        let cube = cube(2);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let v0 = [0.012, -0.007, 0.003];
        let velocity = vec![v0; n];
        let a = [2e4, -1e4, 5e3];
        // Independent P1 uniform-body-force integration: integral N_i dV=V/4.
        let mut force = zero.clone();
        let cell_mass = 7500.0 * LENGTH.powi(3) / (6.0 * 2.0_f64.powi(3));
        for tet in &cube.tets {
            for &node in tet {
                for axis in 0..3 {
                    force[node][axis] += 0.25 * cell_mass * a[axis];
                }
            }
        }
        let body = Body::new(
            &model,
            &[],
            &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
            &zero,
            &velocity,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let dt = 1e-7;
        let candidate = body
            .step_candidate(
                dt,
                &[Condition::Voltage(0.0); 2],
                &force,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let state = candidate.body().state();
        for node in 0..n {
            for axis in 0..3 {
                near(
                    state.displacement_m()[node][axis],
                    dt * v0[axis] + 0.5 * dt * dt * a[axis],
                    1e-19,
                    2e-9,
                );
                near(
                    state.velocity_m_s()[node][axis],
                    v0[axis] + dt * a[axis],
                    1e-12,
                    2e-9,
                );
            }
        }
        let expected_ke = 0.5
            * 7500.0
            * LENGTH.powi(3)
            * (0..3)
                .map(|axis| (v0[axis] + dt * a[axis]).powi(2))
                .sum::<f64>();
        near(state.kinetic_energy_j(), expected_ke, 1e-22, 2e-10);
        assert!(state.stored_energy_j() < 1e-25);
        assert!(candidate.step().energy_residual_j.abs() < 1e-20);
        assert!(candidate.step().electrical_work_j.abs() < 1e-25);
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(body.state().velocity_m_s(), &velocity);
    });
}

#[test]
fn g3_short_and_open_circuit_free_ringing_conserve_energy_and_reverse_time() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let initial_u = cube
            .nodes
            .iter()
            .map(|p| [0.0, 0.0, 1e-6 * p[2]])
            .collect::<Vec<_>>();
        let mut energies = Vec::new();
        for top in [Condition::Voltage(0.0), Condition::Charge(0.0)] {
            let electrodes = electrodes(&cube, 0.0, top);
            let conditions = [Condition::Voltage(0.0), top];
            let mut body = Body::new(
                &model,
                &cube.supports,
                &electrodes,
                &initial_u,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let initial_energy = body.state().total_energy_j();
            energies.push(initial_energy);
            let mut peak_kinetic: f64 = 0.0;
            for _ in 0..40 {
                let (next, step) = body
                    .step_candidate(1e-8, &conditions, &zero, &charge, controls(), cx)
                    .unwrap()
                    .into_parts();
                assert_eq!(step.mechanical_work_j, 0.0);
                near(step.electrical_work_j, 0.0, 1e-22, 0.0);
                near(next.state().total_energy_j(), initial_energy, 1e-22, 2e-8);
                peak_kinetic = peak_kinetic.max(next.state().kinetic_energy_j());
                body = next;
            }
            assert!(peak_kinetic > 0.3 * initial_energy);
            assert_ne!(body.state().displacement_m(), &initial_u);
            if matches!(top, Condition::Charge(_)) {
                assert!(body.state().electrodes()[1].voltage_v.abs() > 1e-4);
                near(body.state().electrodes()[1].charge_c, 0.0, 1e-18, 0.0);
            }
            let reversed_v = body
                .state()
                .velocity_m_s()
                .iter()
                .map(|v| v.map(|x| -x))
                .collect::<Vec<_>>();
            let mut reversed = Body::new(
                &model,
                &cube.supports,
                &electrodes,
                body.state().displacement_m(),
                &reversed_v,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            for _ in 0..40 {
                reversed = reversed
                    .step_candidate(1e-8, &conditions, &zero, &charge, controls(), cx)
                    .unwrap()
                    .into_parts()
                    .0;
            }
            for node in 0..n {
                for axis in 0..3 {
                    near(
                        reversed.state().displacement_m()[node][axis],
                        initial_u[node][axis],
                        2e-17,
                        2e-7,
                    );
                    near(reversed.state().velocity_m_s()[node][axis], 0.0, 2e-9, 0.0);
                }
            }
        }
        // Uniform axial strain: isolated electrodes add e33^2/epsilon33 to C33.
        near(
            energies[1] / energies[0],
            1.0 + 12.0_f64.powi(2) / (EPSILON[2] * STIFFNESS[2]),
            1e-10,
            1e-10,
        );
    });
}

#[test]
fn g1_voltage_and_charge_drives_pay_for_strain_electric_and_kinetic_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        for (initial, target) in [
            (Condition::Voltage(0.0), Condition::Voltage(-10.0)),
            (Condition::Charge(0.0), Condition::Charge(-1e-10)),
        ] {
            let body = Body::new(
                &model,
                &cube.supports,
                &electrodes(&cube, 0.0, initial),
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let candidate = body
                .step_candidate(
                    4e-8,
                    &[Condition::Voltage(0.0), target],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let state = candidate.body().state();
            assert!(state.kinetic_energy_j() > 0.0);
            assert!(state.stored_energy_j() > 0.0);
            assert!(
                cube.top
                    .iter()
                    .map(|&i| state.displacement_m()[i][2])
                    .sum::<f64>()
                    > 0.0
            );
            near(
                state.total_energy_j(),
                candidate.step().electrical_work_j,
                1e-22,
                2e-10,
            );
            assert_eq!(candidate.step().mechanical_work_j, 0.0);
            match target {
                Condition::Voltage(v) => assert_eq!(state.electrodes()[1].voltage_v, v),
                Condition::Charge(q) => near(state.electrodes()[1].charge_c, q, 1e-19, 1e-9),
            }
        }
    });
}

#[test]
fn g1_injected_volume_charge_has_its_own_electrical_supply() {
    with_cx(|cx, _| {
        let cube = cube(2);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let center = cube
            .nodes
            .iter()
            .position(|p| p.iter().all(|&x| x == 0.5 * LENGTH))
            .unwrap();
        let mut next_charge = charge.clone();
        next_charge[center] = 1e-11;
        let candidate = body
            .step_candidate(
                2e-8,
                &[Condition::Voltage(0.0); 2],
                &zero,
                &next_charge,
                controls(),
                cx,
            )
            .unwrap();
        let state = candidate.body().state();
        assert_eq!(state.displacement_m(), &zero);
        assert_eq!(state.velocity_m_s(), &zero);
        assert!(state.voltage_v()[center] > 0.0);
        near(
            state.total_energy_j(),
            0.5 * state.voltage_v()[center] * next_charge[center],
            1e-23,
            1e-10,
        );
        near(
            candidate.step().electrical_work_j,
            state.total_energy_j(),
            1e-23,
            1e-10,
        );
    });
}

#[test]
fn g3_changing_a_large_common_voltage_preserves_the_dynamic_fields() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let body = Body::new(
            &model,
            &cube.supports,
            &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let plain = body
            .step_candidate(
                4e-8,
                &[Condition::Voltage(0.0), Condition::Voltage(-10.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let shifted = body
            .step_candidate(
                4e-8,
                &[Condition::Voltage(1e12), Condition::Voltage(1e12 - 10.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            plain.body().state().displacement_m(),
            shifted.body().state().displacement_m()
        );
        assert_eq!(
            plain.body().state().velocity_m_s(),
            shifted.body().state().velocity_m_s()
        );
        assert_eq!(
            plain.step().electrical_work_j,
            shifted.step().electrical_work_j
        );
        for (a, b) in plain
            .body()
            .state()
            .cells()
            .iter()
            .zip(shifted.body().state().cells())
        {
            assert_eq!(a.response, b.response);
        }
    });
}

#[test]
fn g4_invalid_or_cancelled_candidates_preserve_the_accepted_endpoint() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let electrodes = electrodes(&cube, 0.0, Condition::Voltage(0.0));
        let body = Body::new(
            &model,
            &cube.supports,
            &electrodes,
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let conditions = [Condition::Voltage(0.0), Condition::Voltage(-100.0)];
        let mut narrow = controls();
        narrow.maximum_electric_field_v_m = 10.0;
        assert!(
            body.step_candidate(1e-8, &conditions, &zero, &charge, narrow, cx)
                .is_err()
        );
        assert!(
            body.step_candidate(0.0, &conditions, &zero, &charge, controls(), cx)
                .is_err()
        );
        assert!(
            body.step_candidate(
                1e-8,
                &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                &zero,
                &charge,
                controls(),
                cx
            )
            .is_err()
        );
        let mut invalid_velocity = zero.clone();
        invalid_velocity[0][0] = 1.0;
        assert!(
            Body::new(
                &model,
                &cube.supports,
                &electrodes,
                &zero,
                &invalid_velocity,
                &charge,
                controls(),
                cx
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            body.step_candidate(1e-8, &conditions, &zero, &charge, controls(), cx),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(body.state().displacement_m(), &zero);
        assert_eq!(body.state().velocity_m_s(), &zero);
        assert_eq!(body.state().total_energy_j(), 0.0);
    });
}
