//! Original nonlinear moving laws and independently reduced signed work.
use super::*;
use fs_em::conduction::motion::ConductionMotion;

#[test]
fn moving_field_source_separates_short_circuit_generation_and_motor_work() {
    let law = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.25,
        10.0,
    );
    let model = specimen(law.clone(), law, 0.2);
    let before = model.conductors().to_vec();
    for sign in [-1.0_f64, 1.0] {
        let speed = sign * 2.0;
        let velocity = [[0.0, speed]; 4];
        for voltage in [0.0, 1.0, 3.0] {
            let candidate = model
                .solve_with_motion_candidate(
                    &ends(sign * voltage),
                    &[],
                    ConductionMotion {
                        velocity_m_s: &velocity,
                        normal_flux_density_t: &[1.0; 4],
                    },
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let report = candidate.report();
            // Independent uniform law J=(1+|E'|/4)E', E'=u_y B_z-V/L.
            let ep = sign * (2.0 - voltage);
            let j = (1.0 + 0.25 * ep.abs()) * ep;
            for cell in 0..4 {
                near(report.field.electric_field_v_m[cell][0], -sign * voltage);
                near(report.motion.conductor_field_v_m[cell][0], ep);
                near(report.field.current_density_a_m2[cell][0], j);
                near(report.field.current_density_a_m2[cell][1], 0.0);
                near(report.field.cell_power_w[cell], 0.05 * ep * j);
                near(
                    report.motion.cell_mechanical_power_w[cell],
                    -0.05 * speed * j,
                );
                near(report.motion.lorentz_force_density_n_m3[cell][1], -j);
                let source = &report.field.evaluated[cell];
                near(
                    source.conductor.material().query().axes()["operating"],
                    ep.abs(),
                );
                near(source.power_defect_w, 0.0);
                assert_eq!(
                    source.conductor.source_property().answer().receipt.selected,
                    model.conductors()[cell]
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                );
            }
            near(report.field.power_w, 0.2 * ep * j);
            near(report.motion.mechanical_power_w, -0.2 * speed * j);
            near(report.motion.terminal_power_w, -0.2 * sign * voltage * j);
            near(
                report.motion.electrical_power_w,
                report.motion.terminal_power_w,
            );
            near(report.motion.contact_power_w, 0.0);
            near(report.field.power_balance_residual_w, 0.0);
            near(report.motion.work_pairing_residual_w, 0.0);
            near(report.field.boundary_current_a, 0.0);
            assert!(report.field.field_solves > 1);
            assert!(report.field.power_w > 0.0);
            if voltage == 0.0 {
                near(report.motion.terminal_power_w, 0.0);
                assert!(report.motion.mechanical_power_w < 0.0);
            } else if voltage < 2.0 {
                assert!(report.motion.terminal_power_w < 0.0);
                assert!(report.motion.mechanical_power_w < 0.0);
            } else {
                assert!(report.motion.terminal_power_w > 0.0);
                assert!(report.motion.mechanical_power_w > 0.0);
            }
        }
    }
    assert_eq!(model.conductors(), before);
}

#[test]
fn moving_current_source_and_trace_work_survive_an_unresolvable_common_voltage() {
    let law = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.0,
        0.25,
        20.0,
    );
    let model = specimen(law.clone(), law, 0.2)
        .with_relaxation(0.5)
        .unwrap();
    let velocity = [[0.0, 2.0]; 4];
    let motion = ConductionMotion {
        velocity_m_s: &velocity,
        normal_flux_density_t: &[1.0; 4],
    };
    // Under voltage drive, rho(J)J=E'=1, independently inverted as a quadratic.
    let accepted = model
        .solve_with_motion_candidate(
            &ends(1.0),
            &[],
            motion,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    let expected = 2.0 / (1.0 + 2.0_f64.sqrt());
    for cell in 0..4 {
        near(
            accepted.report().field.current_density_a_m2[cell][0],
            expected,
        );
        near(
            accepted.report().field.evaluated[cell]
                .conductor
                .material()
                .query()
                .axes()["operating"],
            expected,
        );
    }
    for j in [-0.8_f64, 0.8] {
        let ep = (1.0 + 0.25 * j.abs()) * j;
        for bias in [0.0, 2.0_f64.powi(54)] {
            let candidate = model
                .solve_with_motion_candidate(
                    &[(0, bias), (3, bias)],
                    &[ConductionCurrentTrace {
                        vertices: [2, 5],
                        inward_current_density_a_m2: -j,
                    }],
                    motion,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let report = candidate.report();
            assert!(report.field.reference_voltage_v.iter().all(|&v| v == bias));
            near(report.field.relative_voltage_v[2], 2.0 - ep);
            for cell in 0..4 {
                near(report.field.current_density_a_m2[cell][0], j);
                near(report.field.electric_field_v_m[cell][0], ep - 2.0);
                near(report.motion.conductor_field_v_m[cell][0], ep);
                near(report.field.cell_power_w[cell], 0.05 * ep * j);
            }
            near(report.field.power_w, 0.2 * ep * j);
            near(report.motion.mechanical_power_w, -0.4 * j);
            near(report.motion.terminal_power_w, -0.2 * j * (2.0 - ep));
            near(report.field.power_balance_residual_w, 0.0);
            near(report.motion.work_pairing_residual_w, 0.0);
        }
    }
}

#[test]
fn open_circuit_screens_motion_and_zero_velocity_reproduces_stationary_sources() {
    let law = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.25,
        10.0,
    );
    let model = specimen(law.clone(), law, 0.2);
    let open = model
        .solve_with_motion_candidate(
            &[(0, 0.0), (3, 0.0)],
            &[],
            ConductionMotion {
                velocity_m_s: &[[0.0, 2.0]; 4],
                normal_flux_density_t: &[1.0; 4],
            },
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    near(open.report().field.voltage_v[2], 2.0);
    for cell in 0..4 {
        near(open.report().field.electric_field_v_m[cell][0], -2.0);
        near(open.report().motion.conductor_field_v_m[cell][0], 0.0);
        near(open.report().field.current_density_a_m2[cell][0], 0.0);
        near(open.report().field.cell_power_w[cell], 0.0);
        near(
            open.report().field.evaluated[cell]
                .conductor
                .material()
                .query()
                .axes()["operating"],
            0.0,
        );
    }
    near(open.report().motion.mechanical_power_w, 0.0);
    near(open.report().motion.terminal_power_w, 0.0);
    let boundary = ends(1.5);
    let stationary = model
        .solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            false
        })
        .unwrap();
    let zero = model
        .solve_with_motion_candidate(
            &boundary,
            &[],
            ConductionMotion {
                velocity_m_s: &[[0.0; 2]; 4],
                normal_flux_density_t: &[3.0; 4],
            },
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(zero.report().field.voltage_v, stationary.report().voltage_v);
    assert_eq!(
        zero.report().field.cell_power_w,
        stationary.report().cell_power_w
    );
    assert_eq!(zero.report().field.evaluated, stationary.report().evaluated);
    assert_eq!(
        zero.report().field.field_solves,
        stationary.report().field_solves
    );
    assert_eq!(
        zero.report().field.source_evaluations,
        stationary.report().source_evaluations
    );
    near(zero.report().motion.mechanical_power_w, 0.0);
    near(
        zero.report().motion.terminal_power_w,
        stationary.report().power_w,
    );
}

#[test]
fn moving_trials_retain_domain_iteration_and_cancellation_gates() {
    let law = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.25,
        10.0,
    );
    let model = specimen(law.clone(), law, 0.2);
    let before = model.conductors().to_vec();
    let motion = ConductionMotion {
        velocity_m_s: &[[0.0, 2.0]; 4],
        normal_flux_density_t: &[1.0; 4],
    };
    let mut calls = 0;
    let baseline = model
        .solve_with_motion_candidate(
            &ends(1.0),
            &[],
            motion,
            limits(),
            controls(),
            policy(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut interrupted = 0;
    assert!(matches!(
        model.solve_with_motion_candidate(
            &ends(1.0),
            &[],
            motion,
            limits(),
            controls(),
            policy(),
            &mut || {
                interrupted += 1;
                interrupted == calls / 2
            },
        ),
        Err(ConductionResponseError::Cancelled),
    ));
    assert!(matches!(
        model.solve_with_motion_candidate(
            &ends(1.0),
            &[],
            motion,
            limits(),
            controls(),
            ConductionResponseControls {
                maximum_iterations: 1,
                ..policy()
            },
            &mut || false,
        ),
        Err(ConductionResponseError::NotConverged { .. }),
    ));
    let bounded = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.25,
        0.5,
    );
    let bounded = specimen(bounded.clone(), bounded, 0.2);
    let bounded_before = bounded.conductors().to_vec();
    // Laboratory E is zero, but the actual material field E'=2 is unsupported.
    assert!(
        bounded
            .solve_with_motion_candidate(
                &ends(0.0),
                &[],
                motion,
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(bounded.conductors(), bounded_before);
    assert_eq!(model.conductors(), before);
    let retry = model
        .solve_with_motion_candidate(
            &ends(1.0),
            &[],
            motion,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        retry.report().field.evaluated,
        baseline.report().field.evaluated
    );
    assert_eq!(
        retry.report().field.voltage_v,
        baseline.report().field.voltage_v
    );
    assert_eq!(
        retry.report().motion.cell_mechanical_power_w,
        baseline.report().motion.cell_mechanical_power_w
    );
}

fn joint_specimen(depth: f64) -> PlanarConductionResponse {
    let law = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.2,
        10.0,
    );
    let contact = super::super::contacts::interface(law.original(), law.original(), 1e5, 1.0);
    let mesh = TriComplex2::from_indexed_triangles(
        "nonlinear moving bulk and stationary joint",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [2.0, 1.0],
        ],
        vec![
            [0, 1, 2],
            [0, 2, 3],
            [1, 4, 5],
            [1, 5, 2],
            [6, 7, 8],
            [6, 8, 9],
        ],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    PlanarConductionResponse::new_with_phase_contacts(
        mesh,
        vec![law; 6],
        vec![fs_em::conduction::ConductionContact {
            first_vertices: [4, 5],
            second_vertices: [6, 9],
            source: contact,
        }],
        vec![[None, None]],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn moving_nonlinear_bulk_keeps_stationary_joint_heat_and_refuses_sliding_sides() {
    let model = joint_specimen(0.2);
    let velocity = [
        [0.0, 2.0],
        [0.0, 2.0],
        [0.0; 2],
        [0.0; 2],
        [0.0; 2],
        [0.0; 2],
    ];
    for voltage in [0.0_f64, 1.0, 3.0] {
        let candidate = model
            .solve_with_motion_candidate(
                &[(0, 0.0), (3, 0.0), (7, voltage), (8, voltage)],
                &[],
                ConductionMotion {
                    velocity_m_s: &velocity,
                    normal_flux_density_t: &[1.0; 6],
                },
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .unwrap();
        let report = candidate.report();
        // Three unit lengths share J=(1+0.2|E'|)E' and the joint R''=0.5.
        // Therefore |2-V|=3.5|E'|+0.1|E'|², independently of the trial solver.
        let drive = 2.0 - voltage;
        let ep = drive.signum() * 2.0 * drive.abs()
            / (3.5 + (3.5_f64.powi(2) + 0.4 * drive.abs()).sqrt());
        let j = (1.0 + 0.2 * ep.abs()) * ep;
        for cell in 0..6 {
            near(report.field.current_density_a_m2[cell][0], j);
            near(report.motion.conductor_field_v_m[cell][0], ep);
            near(
                report.field.electric_field_v_m[cell][0],
                ep - if cell < 2 { 2.0 } else { 0.0 },
            );
            near(report.field.cell_power_w[cell], 0.1 * ep * j);
            near(
                report.field.evaluated[cell]
                    .conductor
                    .material()
                    .query()
                    .axes()["operating"],
                ep.abs(),
            );
        }
        near(report.field.contact_response[0].current_a, 0.2 * j);
        near(report.field.contact_response[0].power_w, 0.1 * j * j);
        near(report.motion.contact_power_w, 0.1 * j * j);
        near(report.field.power_w, 0.6 * ep * j + 0.1 * j * j);
        near(report.motion.mechanical_power_w, -0.4 * j);
        near(report.motion.electrical_power_w, 0.6 * ep * j - 0.4 * j);
        near(report.motion.terminal_power_w, -0.2 * voltage * j);
        near(report.field.power_balance_residual_w, 0.0);
    }
    for cell in [2, 5] {
        let mut moving_side = velocity;
        moving_side[cell] = [0.0, 1.0];
        assert!(matches!(
            model.solve_with_motion_candidate(
                &[(0, 0.0), (3, 0.0), (7, 1.0), (8, 1.0)],
                &[],
                ConductionMotion {
                    velocity_m_s: &moving_side,
                    normal_flux_density_t: &[0.0; 6],
                },
                limits(),
                controls(),
                policy(),
                &mut || false,
            ),
            Err(ConductionResponseError::Field(ConductionError::Invalid(
                "moving electrical contacts require an owned surface work and heat law"
            )))
        ));
    }
}
