//! Complete physical work derivatives against independent neighboring cycles.

use super::*;
use crate::phs_orbit::PhsPowerObjective;
use std::cell::Cell;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

// Independent exact circular solutions, including the original relaxing pair.
// Supply is evaluated from its own radial/rotation/memory law, not replaced by
// loss or inferred from energy balance, even though their settled means agree.
fn physical_powers(cycle: &Cycle) -> [f64; 2] {
    let a = cycle.radius_squared();
    let w = cycle.frequency();
    let gamma = cycle.eta * cycle.mu;
    let denominator = gamma * gamma + w * w;
    let b = a * (cycle.h * cycle.alpha).powi(2) / denominator;
    let y2 = (cycle.g * cycle.k).powi(2) * a;
    let correlation = a * cycle.h * cycle.alpha * gamma / denominator;
    [
        (cycle.gain - cycle.saturation * y2) * y2 + cycle.mu * cycle.h * cycle.alpha * correlation,
        cycle.r * cycle.k * cycle.k * a + cycle.eta * cycle.mu * cycle.mu * b,
    ]
}

#[test]
fn automatic_complete_work_gradients_include_material_feedback_memory_and_period() {
    let cycle = Cycle::source();
    let slope = Cycle::slope();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian_calls = Cell::new(0);
    let hessian_calls = Cell::new(0);
    let jacobian = |x: &[f64], y: &[f64], out: &mut [f64]| {
        jacobian_calls.set(jacobian_calls.get() + 1);
        cycle.jacobian(x, y, out)
    };
    let hessian = |_: &[f64], v: &[f64], out: &mut [f64]| {
        hessian_calls.set(hessian_calls.get() + 1);
        for i in 0..4 {
            out[i] = v[i] * if i < 2 { cycle.k } else { cycle.mu };
        }
        true
    };
    let storage = |x: &[f64], out: &mut [f64]| {
        for i in 0..4 {
            out[i] = x[i] * if i < 2 { slope.k } else { slope.mu };
        }
        true
    };
    let feedback =
        |x: &[f64], y: &[f64], yp: &[f64], out: &mut [f64]| cycle.parameter(&slope, x, y, yp, out);
    let mut jp = vec![0.0; 16];
    jp[1] = -slope.j;
    jp[4] = slope.j;
    let rp = diagonal([slope.r, slope.r, slope.eta, slope.eta]);
    let gp = diagonal([slope.g, slope.g, slope.h, slope.h]);
    let direction = || PhsAutonomousParameterDirection {
        structure: PhsStructureParameterDirection {
            interconnection: Some(&jp),
            dissipation: Some(&rp),
            ports: Some(&gp),
        },
        storage: Some(&storage),
        feedback: &feedback,
    };
    let epsilon = 2e-5;
    let lower = cycle.shifted(&slope, -epsilon);
    let upper = cycle.shifted(&slope, epsilon);
    let lo = physical_powers(&lower);
    let hi = physical_powers(&upper);
    let base = physical_powers(&cycle);
    let period = core::f64::consts::TAU / orbit.omega;
    let mut derivatives = [0.0; 4];
    for scaled_units in [false, true] {
        let units = if scaled_units {
            HbScales {
                states: vec![2.5, 0.4, 7.0, 0.3],
                equations: vec![0.3, 1.8, 0.6, 3.1],
                angular_frequency: 0.7,
            }
        } else {
            scales()
        };
        let mut prepared = PreparedPhsAutonomousSensitivity::new(
            &system,
            &cycle,
            &jacobian,
            &orbit,
            &reference(&orbit),
            controls(),
            &units,
            || false,
        )
        .unwrap();
        for (index, objective) in [
            PhsPowerObjective::MeanSuppliedPower,
            PhsPowerObjective::MeanResistancePower,
            PhsPowerObjective::SuppliedEnergyPerCycle,
            PhsPowerObjective::ResistanceEnergyPerCycle,
        ]
        .into_iter()
        .enumerate()
        {
            let mut objective_units = objective_controls();
            objective_units.objective_scale = if scaled_units { 7.0 } else { 0.13 };
            let mut adjoint = if scaled_units {
                prepared.prepare_power_adjoint_analytic(
                    objective,
                    objective_units,
                    &jacobian,
                    &hessian,
                    || false,
                )
            } else {
                prepared.prepare_power_adjoint(objective, objective_units, &jacobian, || false)
            }
            .unwrap();
            let factor = if index >= 2 { period } else { 1.0 };
            assert_eq!(adjoint.objective(), objective);
            close(adjoint.value(), factor * base[index % 2], 4e-13);
            close(
                adjoint.angular_frequency_partial(),
                if index >= 2 {
                    -adjoint.value() / orbit.omega
                } else {
                    0.0
                },
                2e-13,
            );
            assert!(adjoint.value_sampling_defect() < 1e-10);
            assert!(adjoint.gradient_sampling_defect() < 1e-10);
            let calls = (jacobian_calls.get(), hessian_calls.get());
            let derivative = adjoint
                .gradient_parameter(0.2, direction(), || false)
                .unwrap();
            let repeat = adjoint
                .gradient_parameter(3.0, direction(), || false)
                .unwrap();
            assert_eq!(calls, (jacobian_calls.get(), hessian_calls.get()));
            close(derivative.derivative, repeat.derivative, 3e-12);
            let expected = if index >= 2 {
                core::f64::consts::TAU
                    * (hi[index % 2] / upper.frequency() - lo[index % 2] / lower.frequency())
                    / (2.0 * epsilon)
            } else {
                (hi[index % 2] - lo[index % 2]) / (2.0 * epsilon)
            };
            close(derivative.derivative, expected, 4e-8);
            if scaled_units {
                close(derivative.derivative, derivatives[index], 5e-12);
            }
            derivatives[index] = derivative.derivative;
        }
    }
    // This direction changes the actual period substantially; treating cycle
    // energy as constant-period mean power would fail this independent oracle.
    assert!((derivatives[2] - period * derivatives[0]).abs() > 0.01);
    assert!(hessian_calls.get() > 0);
}

#[test]
fn automatic_objective_refusals_and_cancellation_preserve_cached_dual() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], y: &[f64], out: &mut [f64]| cycle.jacobian(x, y, out);
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let mut adjoint = prepared
        .prepare_power_adjoint(
            PhsPowerObjective::SuppliedEnergyPerCycle,
            objective_controls(),
            &jacobian,
            || false,
        )
        .unwrap();
    let feedback = |_: &[f64], y: &[f64], yp: &[f64], out: &mut [f64]| {
        assert!(yp.iter().all(|v| *v == 0.0));
        out.copy_from_slice(&[y[0], y[1], 0.0, 0.0]);
        true
    };
    let direction = || PhsAutonomousParameterDirection {
        structure: Default::default(),
        storage: None,
        feedback: &feedback,
    };
    let baseline = adjoint
        .gradient_parameter(1.0, direction(), || false)
        .unwrap();
    let incomplete_storage = |_: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        adjoint.gradient_parameter(
            1.0,
            PhsAutonomousParameterDirection {
                storage: Some(&incomplete_storage),
                ..direction()
            },
            || false
        ),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::ObjectiveRefused { .. }
        ))
    ));
    let incomplete_feedback = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        adjoint.gradient_parameter(
            1.0,
            PhsAutonomousParameterDirection {
                feedback: &incomplete_feedback,
                ..direction()
            },
            || false
        ),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::ObjectiveRefused { .. }
        ))
    ));
    assert!(matches!(
        adjoint.gradient_parameter(
            1.0,
            PhsAutonomousParameterDirection {
                structure: PhsStructureParameterDirection {
                    ports: Some(&[0.0; 2]),
                    ..Default::default()
                },
                ..direction()
            },
            || false
        ),
        Err(PhsAutonomousSensitivityError::Physical(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        adjoint.gradient_parameter(1.0, direction(), || {
            polls += 1;
            polls > 100
        }),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let repeat = adjoint
        .gradient_parameter(1.0, direction(), || false)
        .unwrap();
    assert_eq!(baseline.derivative, repeat.derivative);
}

#[test]
fn source_parameter_power_retains_weak_signed_work_and_nonsymmetric_resistance() {
    let big = 2.0_f64.powi(54);
    for g in [vec![big, 1.0, -big, 0.0], vec![big, -big, 1.0, 0.0]] {
        let system = PortHamiltonian::from_raw_parts(
            2,
            2,
            vec![0.0; 4],
            vec![0.0; 4],
            g.clone(),
            Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
        );
        let partial = system
            .prepare_parameter_direction(PhsStructureParameterDirection {
                ports: Some(&g),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            partial
                .supplied_power_from_effort(&[1.0; 2], &[1.0; 2], &[1.0; 2], &[1.0; 2])
                .unwrap(),
            3.0
        );
        assert!(
            partial
                .supplied_power_from_effort(&[1.0], &[1.0; 2], &[1.0; 2], &[1.0; 2])
                .is_err()
        );
    }
    let raw = |r| {
        PortHamiltonian::from_raw_parts(
            2,
            0,
            vec![0.0; 4],
            r,
            vec![],
            Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
        )
    };
    let system = raw(vec![big, 1.0, -big, 0.0]);
    let rp = [0.5, 0.0, 0.0, 0.5];
    let partial = system
        .prepare_parameter_direction(PhsStructureParameterDirection {
            dissipation: Some(&rp),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        partial
            .resistance_power_from_effort(&[1.0; 2], &[1.0; 2])
            .unwrap(),
        3.0
    );
    let system = raw(vec![1.0, 2.0, 0.0, 3.0]);
    let rp = [0.3, 0.2, 0.2, -0.7];
    let partial = system
        .prepare_parameter_direction(PhsStructureParameterDirection {
            dissipation: Some(&rp),
            ..Default::default()
        })
        .unwrap();
    close(
        partial
            .resistance_power_from_effort(&[2.0, -1.0], &[-0.5, 4.0])
            .unwrap(),
        -9.3,
        2e-15,
    );
    assert!(
        partial
            .resistance_power_from_effort(&[1e-300, 0.0], &[1e-300, 0.0])
            .is_err()
    );
}
