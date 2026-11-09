//! Original source support, finite driver power, and reciprocal Miller storage.

use super::*;

#[path = "gate/current.rs"]
mod current;
use fs_phs::circuit::converter::commutation::capacitance::gate::{
    GateCommutationParameters, GateCommutationStorage,
};
use fs_phs::circuit::converter::commutation::gate::{GateConductanceCurve, GateDriver};

const CAP: [f64; 4] = [0.02, 0.03, 0.025, 0.015];
const GS: [f64; 4] = [0.003, 0.004, 0.005, 0.006];
const GD: [f64; 4] = [0.001, 0.002, 0.0015, 0.0025];
const DRIVER: GateDriver = GateDriver {
    on_voltage_v: 3.0,
    off_voltage_v: 0.0,
    on_resistance_ohm: 2.0,
    off_resistance_ohm: 3.0,
    maximum_current_a: 2.0,
};
fn law() -> GateConductanceCurve {
    GateConductanceCurve::try_new(vec![
        (-4.0, 0.025),
        (0.0, 0.025),
        (1.0, 0.025),
        (2.0, 5.0),
        (4.0, 5.0),
    ])
    .unwrap()
}
fn gate_parameters() -> GateCommutationParameters {
    GateCommutationParameters {
        converter: parameters([0.04, 0.05], CAP).converter,
        intrinsic_drain_charge_curves: CAP.map(curve),
        gate_source_capacitance_f: GS,
        gate_drain_capacitance_f: GD,
        channel_conductance: core::array::from_fn(|_| law()),
        gate_drivers: [DRIVER; 4],
        maximum_cells_per_chord: 64,
    }
}
fn authored_energy(p: &GateCommutationParameters, v: [f64; 9]) -> f64 {
    let vd = [v[1] - v[3], v[3], v[2] - v[4], v[4]];
    let mut result = 0.5 * p.converter.inductance_h * v[0] * v[0];
    for k in 0..2 {
        result += 0.5 * p.converter.bus_capacitance_f[k] * v[1 + k] * v[1 + k];
    }
    for k in 0..4 {
        result += energy(CAP[k], vd[k]);
        result += 0.5 * GS[k] * v[5 + k] * v[5 + k];
        result += 0.5 * GD[k] * (v[5 + k] - vd[k]).powi(2);
    }
    result
}

#[test]
fn reciprocal_miller_charge_and_energy_match_independent_physical_electrodes() {
    let p = gate_parameters();
    let va = [0.3, 0.8, 0.7, 0.3, 0.2, 0.6, 0.4, 0.8, 0.5];
    let vb = [0.1, 0.9, 0.6, 0.35, 0.25, 0.7, 0.3, 0.9, 0.4];
    let a = p
        .state(
            va[0],
            [va[1], va[2]],
            [va[3], va[4]],
            va[5..].try_into().unwrap(),
        )
        .unwrap();
    let b = p
        .state(
            vb[0],
            [vb[1], vb[2]],
            [vb[3], vb[4]],
            vb[5..].try_into().unwrap(),
        )
        .unwrap();
    let vd = [va[1] - va[3], va[3], va[2] - va[4], va[4]];
    let drain =
        core::array::from_fn::<_, 4, _>(|k| charge(CAP[k], vd[k]) + GD[k] * (vd[k] - va[5 + k]));
    for k in 0..2 {
        close(
            a[1 + k],
            p.converter.bus_capacitance_f[k] * va[1 + k] + drain[2 * k],
            2e-16,
        );
        close(a[3 + k], drain[2 * k + 1] - drain[2 * k], 2e-16);
    }
    for k in 0..4 {
        close(a[5 + k], (GS[k] + GD[k]) * va[5 + k] - GD[k] * vd[k], 2e-17);
    }
    let storage = GateCommutationStorage::new(&p).unwrap();
    let mut gradient = [0.0; 9];
    storage.gradient(&a, &mut gradient);
    for k in 0..9 {
        close(gradient[k], va[k], 2e-12);
    }
    let ha = authored_energy(&p, va);
    let hb = authored_energy(&p, vb);
    close(storage.hamiltonian(&a), ha, 2e-13);
    close(storage.hamiltonian(&b), hb, 2e-13);
    let midpoint: [f64; 9] = core::array::from_fn(|k| f64::midpoint(a[k], b[k]));
    let mut mean = [0.0; 9];
    ConstitutiveDiscreteGradient::evaluate(&storage, &a, &b, &midpoint, &mut mean);
    for k in 0..9 {
        close(mean[k], f64::midpoint(va[k], vb[k]), 2e-12);
    }
    close(
        (0..9).map(|k| mean[k] * (b[k] - a[k])).sum(),
        hb - ha,
        3e-13,
    );
    let mut bus_direction = [0.0; 9];
    bus_direction[1] = 1.0;
    let mut gate_direction = [0.0; 9];
    gate_direction[5] = 1.0;
    let mut bus_tangent = [0.0; 9];
    let mut gate_tangent = [0.0; 9];
    assert!(ContinuousStorageTangent::apply(
        &storage,
        &a,
        &bus_direction,
        &mut bus_tangent
    ));
    assert!(ContinuousStorageTangent::apply(
        &storage,
        &a,
        &gate_direction,
        &mut gate_tangent
    ));
    let alpha = GD[0] / (GS[0] + GD[0]);
    let ch = 3.0 * CAP[0] + GS[0] * alpha;
    let cl = 3.0 * CAP[1] + GS[1] * GD[1] / (GS[1] + GD[1]);
    let det = p.converter.bus_capacitance_f[0] * (ch + cl) + ch * cl;
    let reciprocal = alpha * cl / det;
    assert!(reciprocal > 0.0);
    close(bus_tangent[5], reciprocal, 2e-12);
    close(gate_tangent[1], reciprocal, 2e-12);
    let mut raised = [0.6, 0.4, 0.8, 0.5];
    raised[0] += alpha * 0.1;
    let drain_step = p.state(0.3, [0.9, 0.7], [0.3, 0.2], raised).unwrap();
    close(drain_step[5], a[5], 2e-17);
    storage.gradient(&drain_step, &mut gradient);
    close(gradient[5] - va[5], alpha * 0.1, 2e-12);
}

#[test]
fn finite_driver_power_and_original_channel_tangent_preserve_signed_work_and_limits() {
    let law = law();
    let point = law.evaluate(1.5).unwrap();
    close(point.conductance_s, 2.5125, 1e-14);
    close(point.slope_s_per_v, 4.975, 1e-14);
    let epsilon = 1e-6;
    close(
        (law.evaluate(1.5 + epsilon).unwrap().conductance_s
            - law.evaluate(1.5 - epsilon).unwrap().conductance_s)
            / (2.0 * epsilon),
        point.slope_s_per_v,
        2e-9,
    );
    assert!(law.evaluate(4.01).is_err());
    assert!(law.evaluate(f64::NAN).is_err());
    assert!(GateConductanceCurve::try_new(vec![(0.0, 0.0), (1.0, 1.0)]).is_err());
    assert!(GateConductanceCurve::try_new(vec![(0.0, 1.0), (0.0, 2.0)]).is_err());
    for on in [true, false] {
        let point = DRIVER.evaluate(on, 0.7).unwrap();
        let r = if on {
            DRIVER.on_resistance_ohm
        } else {
            DRIVER.off_resistance_ohm
        };
        let u = if on {
            DRIVER.on_voltage_v
        } else {
            DRIVER.off_voltage_v
        };
        let current = (u - 0.7) / r;
        close(point.current_a, current, 1e-15);
        close(point.source_power_w, u * current, 2e-15);
        close(point.terminal_power_w, 0.7 * current, 2e-15);
        close(point.resistor_power_w, r * current * current, 2e-15);
        close(
            point.source_power_w - point.terminal_power_w,
            point.resistor_power_w,
            2e-15,
        );
        assert_eq!(point.current_a.is_sign_positive(), on);
    }
    let returning = GateDriver {
        off_voltage_v: 1.0,
        ..DRIVER
    }
    .evaluate(false, 2.0)
    .unwrap();
    assert!(returning.source_power_w < 0.0 && returning.resistor_power_w > 0.0);
    assert!(
        GateDriver {
            maximum_current_a: 0.01,
            ..DRIVER
        }
        .evaluate(true, 0.0)
        .is_err()
    );
}

#[test]
fn gate_source_support_and_capacitance_or_chord_budgets_refuse_without_regularization() {
    let p = gate_parameters();
    assert!(p.state(0.0, [0.0; 2], [0.0; 2], [4.01; 4]).is_err());
    let mut invalid = p.clone();
    invalid.gate_source_capacitance_f[0] = 0.0;
    assert!(invalid.validate().is_err());
    invalid = p.clone();
    invalid.gate_drain_capacitance_f[0] = -0.001;
    assert!(invalid.validate().is_err());
    let mut zero = p.clone();
    zero.gate_drain_capacitance_f = [0.0; 4];
    assert!(zero.validate().is_ok());
    let origin = zero.state(0.0, [0.0; 2], [0.0; 2], [0.0; 4]).unwrap();
    assert_eq!(origin, [0.0; 9]);
    let zero_storage = GateCommutationStorage::new(&zero).unwrap();
    assert_eq!(zero_storage.hamiltonian(&origin), 0.0);
    let mut origin_effort = [1.0; 9];
    zero_storage.gradient(&origin, &mut origin_effort);
    assert_eq!(origin_effort, [0.0; 9]);
    let a = p.state(0.0, [4.0, -2.0], [-0.6, 3.0], [0.0; 4]).unwrap();
    let b = p.state(0.0, [-3.0, 5.0], [2.0, -2.5], [0.0; 4]).unwrap();
    let midpoint: [f64; 9] = core::array::from_fn(|k| f64::midpoint(a[k], b[k]));
    let storage = GateCommutationStorage::new(&p).unwrap();
    let mut output = [0.0; 9];
    ConstitutiveDiscreteGradient::evaluate(&storage, &a, &b, &midpoint, &mut output);
    assert!(output.iter().all(|x| x.is_finite()));
    let mut limited = p.clone();
    limited.maximum_cells_per_chord = 1;
    let limited = GateCommutationStorage::new(&limited).unwrap();
    ConstitutiveDiscreteGradient::evaluate(&limited, &a, &b, &midpoint, &mut output);
    assert!(!output.iter().all(|x| x.is_finite()));
}
