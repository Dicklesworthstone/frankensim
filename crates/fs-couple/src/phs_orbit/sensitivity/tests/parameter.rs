//! G1/G3: complete physical directions and the prescribed-clock chain.
use super::*;
use crate::phs_orbit::PhsPowerObjective;
use fs_phs::{PhsStructureParameterDirection, QuadraticStorage};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

#[test]
fn simultaneous_structure_storage_ports_drive_and_clock_match_original_mechanics() {
    let (mass, stiffness, loss, coupling, port, omega) = (2.0, 10.0, 1.3, 0.8, 1.2, 1.1);
    let (dm, dk, dr, dj, dg, dw) = (0.3, -0.4, 0.2, 0.07, -0.09, 0.31);
    let source = PortHamiltonian::new(
        2,
        1,
        vec![0.0, coupling, -coupling, 0.0],
        vec![0.0, 0.0, 0.0, loss],
        vec![0.0, port],
        Box::new(QuadraticStorage::new(vec![stiffness, 0.0, 0.0, 1.0 / mass], 2).unwrap()),
    )
    .unwrap();
    let forces = vec![
        C64::from_re(0.4),
        C64::new(0.5, -0.2),
        C64::new(-0.12, 0.07),
    ];
    let df = vec![
        C64::from_re(-0.1),
        C64::new(0.05, 0.08),
        C64::new(0.03, -0.02),
    ];
    let input = PeriodicPortInput::new(omega, vec![forces.clone()]).unwrap();
    let derivative = PeriodicPortInput::new(omega, vec![df.clone()]).unwrap();
    let orbit = solve_forced_harmonic_balance(&source, &input, 0.1, &hb(2)).unwrap();
    let mut prepared =
        PreparedPhsForcedSensitivity::new(&source, &input, &orbit, controls(), &scales(2), || {
            false
        })
        .unwrap();
    // The original drive is retained without changing the constructor lifetime.
    drop(input);
    let j = [0.0, dj, -dj, 0.0];
    let r = [0.0, 0.0, 0.0, dr];
    let g = [0.0, dg];
    let storage = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[dk * x[0], -dm * x[1] / mass / mass]);
        true
    };
    let physical = || PhsForcedParameterDirection {
        structure: PhsStructureParameterDirection {
            interconnection: Some(&j),
            dissipation: Some(&r),
            ports: Some(&g),
        },
        storage: Some(&storage),
        input: Some(&derivative),
    };
    let chart = ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 0.7,
        angular_frequency_derivative: dw,
    };
    let response = prepared
        .solve_parameter(chart, physical(), || false)
        .unwrap();
    close(response.angular_frequency_derivative, dw, 0.0);
    let (mut explicit, mut total) = (0.0, 0.0);
    let mut powers = [0.0; 4]; // supplied, resistance, and their total derivatives
    for harmonic in 0..=2 {
        let w = harmonic as f64 * omega;
        let wp = harmonic as f64 * dw;
        let z = C64::new(coupling * coupling * stiffness - mass * w * w, loss * w);
        let zp = C64::new(
            2.0 * coupling * dj * stiffness + coupling * coupling * dk
                - dm * w * w
                - 2.0 * mass * w * wp,
            dr * w + loss * wp,
        );
        let numerator = forces[harmonic].scale(coupling * port);
        let np =
            forces[harmonic].scale(dj * port + coupling * dg) + df[harmonic].scale(coupling * port);
        let q = numerator / z;
        let qp = (np - q * zp) / z;
        let p = q * C64::new(0.0, mass * w / coupling);
        let pp = qp * C64::new(0.0, mass * w / coupling)
            + q * C64::new(
                0.0,
                (dm * w + mass * wp) / coupling - mass * w * dj / coupling / coupling,
            );
        // Independent complex mechanics for the actual conjugate port flow.
        // DC has multiplicity one, positive-frequency half amplitudes two.
        let velocity = p.scale(1.0 / mass);
        let velocity_p = pp.scale(1.0 / mass) - p.scale(dm / mass / mass);
        let dot = |a: C64, b: C64| a.re * b.re + a.im * b.im;
        let multiplicity = if harmonic == 0 { 1.0 } else { 2.0 };
        powers[0] += multiplicity * port * dot(velocity, forces[harmonic]);
        powers[1] += multiplicity * loss * dot(velocity, velocity);
        powers[2] += multiplicity
            * (dg * dot(velocity, forces[harmonic])
                + port * dot(velocity_p, forces[harmonic])
                + port * dot(velocity, df[harmonic]));
        powers[3] +=
            multiplicity * (dr * dot(velocity, velocity) + 2.0 * loss * dot(velocity, velocity_p));
        for (actual, expected) in [
            (response.coeffs[0][harmonic], qp),
            (response.coeffs[1][harmonic], pp),
        ] {
            close(actual.re, expected.re, 3e-12);
            close(actual.im, expected.im, 3e-12);
        }
        let weight = if harmonic == 0 { 0.5 } else { 1.0 };
        let direct = weight
            * (dk * (q.re * q.re + q.im * q.im) - dm / mass / mass * (p.re * p.re + p.im * p.im));
        explicit += direct;
        total += direct
            + 2.0
                * weight
                * (stiffness * (q.re * qp.re + q.im * qp.im)
                    + (p.re * pp.re + p.im * pp.im) / mass);
    }
    let objective = prepared
        .mean_stored_energy(objective_controls(), || false)
        .unwrap();
    let mut adjoint = prepared
        .prepare_adjoint(&objective.gradient, 1.0, || false)
        .unwrap();
    close(
        adjoint
            .gradient_parameter(chart, explicit, physical(), || false)
            .unwrap()
            .derivative,
        total,
        3e-12,
    );
    drop(adjoint);
    // Per-cycle stored-energy integral: explicit period change must appear once.
    let period = core::f64::consts::TAU / omega;
    let gradient: Vec<_> = objective.gradient.iter().map(|x| period * x).collect();
    let direct = period * explicit - period * dw / omega * objective.value;
    let mut integral = prepared
        .prepare_adjoint(&gradient, period, || false)
        .unwrap();
    close(
        integral
            .gradient_parameter(chart, direct, physical(), || false)
            .unwrap()
            .derivative,
        period * total - period * dw / omega * objective.value,
        5e-12,
    );
    drop(integral);
    // Every built-in objective now assembles its own direct J/R/G/storage/input
    // partial and prescribed period chain. No manual explicit q_p is supplied.
    for (kind, energy, index) in [
        (PhsPowerObjective::MeanSuppliedPower, false, 0),
        (PhsPowerObjective::MeanResistancePower, false, 1),
        (PhsPowerObjective::SuppliedEnergyPerCycle, true, 0),
        (PhsPowerObjective::ResistanceEnergyPerCycle, true, 1),
    ] {
        let mut power = prepared
            .prepare_power_adjoint(kind, objective_controls(), || false)
            .unwrap();
        let value = if energy {
            period * powers[index]
        } else {
            powers[index]
        };
        let expected = if energy {
            period * (powers[index + 2] - dw / omega * powers[index])
        } else {
            powers[index + 2]
        };
        close(power.value(), value, 4e-12);
        close(
            power
                .gradient_parameter(chart, physical(), || false)
                .unwrap()
                .derivative,
            expected,
            8e-12,
        );
    }
}

#[test]
fn nonlinear_memory_and_driving_clock_keep_generated_harmonics() {
    let (db, dj, dr, dg, dw, du) = (0.2, 0.04, -0.03, 0.06, 0.13, -0.08);
    let system = |p: f64| {
        PortHamiltonian::new(
            3,
            1,
            vec![
                0.0,
                1.0 + dj * p,
                0.0,
                -1.0 - dj * p,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![0.0, 0.0, 0.0, 0.0, 0.4 + dr * p, 0.0, 0.0, 0.0, 2.0],
            vec![0.0, 1.0 + dg * p, 0.0],
            Box::new(NonlinearMemory { beta: 1.4 + db * p }),
        )
        .unwrap()
    };
    let input = |p: f64| {
        PeriodicPortInput::new(
            0.85 + dw * p,
            vec![vec![
                C64::from_re(0.05 * p),
                C64::new(0.7, -0.2).scale(1.0 + du * p),
            ]],
        )
        .unwrap()
    };
    let source = system(0.0);
    let base_input = input(0.0);
    let orbit = memory_orbit(&source, &base_input);
    let mut prepared = PreparedPhsForcedSensitivity::new(
        &source,
        &base_input,
        &orbit,
        controls(),
        &scales(3),
        || false,
    )
    .unwrap();
    let j = [0.0, dj, 0.0, -dj, 0.0, 0.0, 0.0, 0.0, 0.0];
    let r = [0.0, 0.0, 0.0, 0.0, dr, 0.0, 0.0, 0.0, 0.0];
    let g = [0.0, dg, 0.0];
    let forcing = PeriodicPortInput::new(
        0.85,
        vec![vec![C64::from_re(0.05), C64::new(0.7, -0.2).scale(du)]],
    )
    .unwrap();
    let storage = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[db * x[0].powi(3), 0.0, 0.0]);
        true
    };
    let physical = || PhsForcedParameterDirection {
        structure: PhsStructureParameterDirection {
            interconnection: Some(&j),
            dissipation: Some(&r),
            ports: Some(&g),
        },
        storage: Some(&storage),
        input: Some(&forcing),
    };
    let chart = ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 1.0,
        angular_frequency_derivative: dw,
    };
    let response = prepared
        .solve_parameter(chart, physical(), || false)
        .unwrap();
    let step = 1e-4;
    let lo = memory_orbit(&system(-step), &input(-step));
    let hi = memory_orbit(&system(step), &input(step));
    for state in 0..3 {
        for harmonic in 0..=9 {
            let reference =
                (hi.coeffs[state][harmonic] - lo.coeffs[state][harmonic]).scale(0.5 / step);
            close(response.coeffs[state][harmonic].re, reference.re, 3e-9);
            close(response.coeffs[state][harmonic].im, reference.im, 3e-9);
        }
    }
    assert!(response.coeffs[2][3].abs() > 1e-5);
    let mean = prepared
        .mean_stored_energy(objective_controls(), || false)
        .unwrap();
    let count = 1009;
    let explicit = (0..count)
        .map(|i| {
            0.25 * db
                * orbit
                    .sample(0, core::f64::consts::TAU * i as f64 / count as f64)
                    .powi(4)
        })
        .sum::<f64>()
        / count as f64;
    let energy = |source: &PortHamiltonian, orbit: &HbOrbit| {
        (0..count)
            .map(|i| {
                let theta = core::f64::consts::TAU * i as f64 / count as f64;
                source.hamiltonian(&[
                    orbit.sample(0, theta),
                    orbit.sample(1, theta),
                    orbit.sample(2, theta),
                ])
            })
            .sum::<f64>()
            / count as f64
    };
    let mut adjoint = prepared
        .prepare_adjoint(&mean.gradient, 1.0, || false)
        .unwrap();
    let derivative = adjoint
        .gradient_parameter(chart, explicit, physical(), || false)
        .unwrap();
    close(
        derivative.derivative,
        (energy(&system(step), &hi) - energy(&system(-step), &lo)) / (2.0 * step),
        5e-9,
    );
    drop(adjoint);
    // Independent dense phase quadrature of the original physical powers,
    // including the actual relaxing-state resistance, in separately solved
    // neighboring systems. Neither oracle calls the new power/adjoint APIs.
    let power_mean = |p: f64, cycle: &HbOrbit, supplied: bool| {
        (0..count)
            .map(|i| {
                let theta = core::f64::consts::TAU * i as f64 / count as f64;
                let q = cycle.sample(0, theta);
                let momentum = cycle.sample(1, theta);
                let memory = cycle.sample(2, theta);
                if supplied {
                    let force = 0.05 * p
                        + (1.0 + du * p)
                            * (1.4 * fs_math::det::cos(theta) + 0.4 * fs_math::det::sin(theta));
                    (1.0 + dg * p) * momentum * force
                } else {
                    (0.4 + dr * p) * momentum * momentum + 0.5 * (memory - q).powi(2)
                }
            })
            .sum::<f64>()
            / count as f64
    };
    for (kind, supplied, energy) in [
        (PhsPowerObjective::MeanSuppliedPower, true, false),
        (PhsPowerObjective::MeanResistancePower, false, false),
        (PhsPowerObjective::SuppliedEnergyPerCycle, true, true),
        (PhsPowerObjective::ResistanceEnergyPerCycle, false, true),
    ] {
        let reference = |p: f64, cycle: &HbOrbit| {
            power_mean(p, cycle, supplied)
                * if energy {
                    core::f64::consts::TAU / cycle.omega
                } else {
                    1.0
                }
        };
        let mut power = prepared
            .prepare_power_adjoint(kind, objective_controls(), || false)
            .unwrap();
        close(power.value(), reference(0.0, &orbit), 2e-11);
        close(
            power
                .gradient_parameter(chart, physical(), || false)
                .unwrap()
                .derivative,
            (reference(step, &hi) - reference(-step, &lo)) / (2.0 * step),
            2e-8,
        );
    }
}

#[test]
fn complete_directions_refuse_bad_sources_and_cancel_without_losing_reuse() {
    let source = fs_phs::mass_spring_damper(2.0, 10.0, 1.3).unwrap();
    let input = PeriodicPortInput::new(1.1, vec![vec![C64::ZERO, C64::new(0.5, -0.2)]]).unwrap();
    let orbit = solve_forced_harmonic_balance(&source, &input, 0.1, &hb(2)).unwrap();
    let mut prepared =
        PreparedPhsForcedSensitivity::new(&source, &input, &orbit, controls(), &scales(2), || {
            false
        })
        .unwrap();
    let physical = || PhsForcedParameterDirection {
        structure: PhsStructureParameterDirection::default(),
        storage: None,
        input: None,
    };
    let chart = ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 1.0,
        angular_frequency_derivative: 0.3,
    };
    let expected = prepared
        .solve_parameter(chart, physical(), || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.solve_parameter(chart, physical(), || {
            polls += 1;
            polls > 10
        }),
        Err(PhsForcedSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let bad_storage = |_: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(
        prepared
            .solve_parameter(
                chart,
                PhsForcedParameterDirection {
                    storage: Some(&bad_storage),
                    ..physical()
                },
                || false
            )
            .is_err()
    );
    let bad_j = [0.0, 1.0, 1.0, 0.0];
    assert!(matches!(
        prepared.solve_parameter(
            chart,
            PhsForcedParameterDirection {
                structure: PhsStructureParameterDirection {
                    interconnection: Some(&bad_j),
                    ..Default::default()
                },
                ..physical()
            },
            || false
        ),
        Err(PhsForcedSensitivityError::Physical(_))
    ));
    let bad_clock = input.at_frequency(1.2).unwrap();
    assert!(
        prepared
            .solve_parameter(
                chart,
                PhsForcedParameterDirection {
                    input: Some(&bad_clock),
                    ..physical()
                },
                || false
            )
            .is_err()
    );
    let wide = PeriodicPortInput::new(
        1.1,
        vec![vec![
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::from_re(f64::from_bits(1)),
        ]],
    )
    .unwrap();
    assert!(
        prepared
            .solve_parameter(
                chart,
                PhsForcedParameterDirection {
                    input: Some(&wide),
                    ..physical()
                },
                || false
            )
            .is_err()
    );
    assert_eq!(
        expected.coeffs,
        prepared
            .solve_parameter(chart, physical(), || false)
            .unwrap()
            .coeffs
    );
    let mut power = prepared
        .prepare_power_adjoint(
            PhsPowerObjective::SuppliedEnergyPerCycle,
            objective_controls(),
            || false,
        )
        .unwrap();
    let complete = power
        .gradient_parameter(chart, physical(), || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        power.gradient_parameter(chart, physical(), || {
            polls += 1;
            polls > 10
        }),
        Err(PhsForcedSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    assert!(
        power
            .gradient_parameter(
                chart,
                PhsForcedParameterDirection {
                    storage: Some(&bad_storage),
                    ..physical()
                },
                || false,
            )
            .is_err()
    );
    assert!(
        power
            .gradient_parameter(
                chart,
                PhsForcedParameterDirection {
                    input: Some(&bad_clock),
                    ..physical()
                },
                || false,
            )
            .is_err()
    );
    assert_eq!(
        complete,
        power
            .gradient_parameter(chart, physical(), || false)
            .unwrap()
    );
    drop(power);
    let mut tiny = controls();
    tiny.maximum_payload_bytes = 1;
    assert!(matches!(
        PreparedPhsForcedSensitivity::new(&source, &input, &orbit, tiny, &scales(2), || false),
        Err(PhsForcedSensitivityError::Sensitivity(
            ForcedSensitivityError::Budget { .. }
        ))
    ));
}
