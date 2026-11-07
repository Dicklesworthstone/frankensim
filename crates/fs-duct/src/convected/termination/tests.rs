//! G1/G3: independent boundary equations, stationary parity and physical work.

use super::*;
use crate::Segment;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use std::cell::Cell;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn cylinder() -> Segment {
    Segment::Cylinder {
        radius: 0.012,
        length: 0.83,
    }
}

fn omega() -> f64 {
    core::f64::consts::TAU * 723.0
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale.max(1e-30),
        "actual {actual:e}, expected {expected:e}, scale {scale:e}"
    );
}

fn complex_close(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale.max(1e-30),
        "actual {actual:?}, expected {expected:?}, scale {scale:e}"
    );
}

fn loads(reference: f64) -> [FlowImpedance; 2] {
    [C64::new(0.7, 0.2), C64::new(1.6, -0.4)]
        .map(|z| FlowImpedance::try_new(z.scale(reference)).unwrap())
}

#[test]
fn zero_flow_drives_recover_original_pressure_volume_flow_transfer() {
    let gas = gas();
    let geometry = cylinder();
    let duct = UniformFlowCylinder::try_new(&gas, &geometry, omega(), 0.0).unwrap();
    let reference = duct
        .characteristic_flow_impedance()
        .unwrap()
        .value_j_s_kg2()
        .re;
    let loads = loads(reference);
    let prepared = duct.prepare_terminations(loads, 1e8).unwrap();
    let drives = [C64::new(1.2, -0.3), C64::new(-0.2, 0.4)];
    let actual = prepared.solve(drives).unwrap();
    // Independent stationary p/Q calculation through the original ABCD owner.
    // H=p/rho, m_out=+/-rho Q, hence Z_pQ=rho² Z_Hm and p_drive=rho H_drive.
    let wave = crate::fluid::segment_wave(&gas, duct.radius_m(), omega()).unwrap();
    let matrix = crate::segment_matrix(&geometry, &wave).unwrap();
    let z = loads.map(|load| load.value_j_s_kg2().scale(gas.density * gas.density));
    let source = drives.map(|h| h.scale(gas.density));
    let a = matrix[0] + z[0] * matrix[2];
    let b = matrix[1] + z[0] * matrix[3];
    let flow_right = (source[0] - a * source[1]) / (a * z[1] + b);
    let pressure_right = source[1] + z[1] * flow_right;
    let pressure_left = matrix[0] * pressure_right + matrix[1] * flow_right;
    let flow_left = matrix[2] * pressure_right + matrix[3] * flow_right;
    let left = actual.at_position(0.0).unwrap();
    let right = actual.at_position(duct.length_m()).unwrap();
    for (actual, expected) in [
        (left.pressure_pa, pressure_left),
        (right.pressure_pa, pressure_right),
        (left.volume_flow_m3_s, flow_left),
        (right.volume_flow_m3_s, flow_right),
    ] {
        complex_close(actual, expected, expected.abs());
    }
    let absorbed = 0.5 * (z[0].re * flow_left.norm_sq() + z[1].re * flow_right.norm_sq());
    let work = 0.5 * (source[0] * flow_left.conj() - source[1] * flow_right.conj()).re;
    close(
        actual.power().absorbed_disturbance_power_w,
        absorbed,
        absorbed,
    );
    close(actual.power().source_work_into_duct_w, work, absorbed);
    close(actual.power().balance_residual_w, 0.0, absorbed);
}

#[test]
fn original_pressure_velocity_fields_satisfy_both_mean_flow_loads_and_power() {
    let gas = gas();
    for mach in [-0.65, -0.2, 0.3, 0.65] {
        let velocity = mach * gas.sound_speed;
        let duct = UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), velocity).unwrap();
        let reference = duct
            .characteristic_flow_impedance()
            .unwrap()
            .value_j_s_kg2()
            .re;
        let loads = loads(reference);
        let prepared = duct.prepare_terminations(loads, 1e8).unwrap();
        let drives = [C64::new(0.9, -0.2), C64::new(-0.4, 0.6)];
        let result = prepared.solve(drives).unwrap();
        let mut absorbed = 0.0;
        let mut source_work = 0.0;
        let mut outward_flux = 0.0;
        for (i, x) in [0.0, duct.length_m()].into_iter().enumerate() {
            let field = result.at_position(x).unwrap();
            // Independent physical definitions, without using the reported H/m.
            let enthalpy = field.pressure_pa.scale(1.0 / gas.density)
                + field.particle_velocity_m_s.scale(velocity);
            let mass = (field.particle_velocity_m_s.scale(gas.density)
                + field
                    .pressure_pa
                    .scale(velocity / (gas.sound_speed * gas.sound_speed)))
            .scale(duct.area_m2() * if i == 0 { -1.0 } else { 1.0 });
            let load = loads[i].value_j_s_kg2();
            complex_close(enthalpy, load * mass + drives[i], drives[i].abs());
            let loss = 0.5 * load.re * mass.norm_sq();
            let work = -0.5 * (drives[i] * mass.conj()).re;
            let flux = 0.5 * (enthalpy * mass.conj()).re;
            let end = result.terminals()[i];
            complex_close(end.enthalpy_j_kg, enthalpy, enthalpy.abs());
            complex_close(end.outward_mass_flow_kg_s, mass, mass.abs());
            complex_close(end.boundary_residual_j_kg, C64::ZERO, drives[i].abs());
            close(end.absorbed_disturbance_power_w, loss, loss);
            close(end.source_work_into_duct_w, work, loss + work.abs());
            close(end.outward_disturbance_power_w, flux, loss + work.abs());
            assert!(end.absorbed_disturbance_power_w >= 0.0);
            absorbed += loss;
            source_work += work;
            outward_flux += flux;
        }
        assert!(absorbed > 0.0);
        close(source_work, absorbed, absorbed);
        close(outward_flux, 0.0, absorbed);
        close(result.power().balance_residual_w, 0.0, absorbed);
        close(result.power().interior_flux_residual_w, 0.0, absorbed);
    }
}

#[test]
fn passive_load_can_reflect_pressure_gain_without_disturbance_power_gain() {
    let gas = gas();
    let mach = 0.5;
    let duct =
        UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), mach * gas.sound_speed).unwrap();
    let matched = duct.characteristic_flow_impedance().unwrap();
    let reactive = FlowImpedance::try_new(C64::new(0.0, matched.value_j_s_kg2().re)).unwrap();
    let prepared = duct.prepare_terminations([matched, reactive], 1e8).unwrap();
    let result = prepared
        .solve([C64::from_re(2.0 * (1.0 + mach) / gas.density), C64::ZERO])
        .unwrap();
    let forward = result.field().boundary_power().outgoing_pressure_pa[1];
    let reflected = result.at_position(duct.length_m()).unwrap().pressure_pa - forward;
    // Z=iR gives enthalpy reflection +i. Pressure reflection additionally
    // carries (1+M)/(1-M)=3, although reflected physical power is unchanged.
    complex_close(prepared.enthalpy_reflection()[1], C64::new(0.0, 1.0), 1.0);
    complex_close(reflected / forward, C64::new(0.0, 3.0), 3.0);
    let factor = duct.area_m2() / (2.0 * gas.density * gas.sound_speed);
    let incident_power = factor * (1.0 + mach).powi(2) * forward.norm_sq();
    let reflected_power = factor * (1.0 - mach).powi(2) * reflected.norm_sq();
    close(incident_power, reflected_power, incident_power);
    assert_eq!(result.terminals()[1].absorbed_disturbance_power_w, 0.0);
    close(
        result.terminals()[1].outward_disturbance_power_w,
        0.0,
        incident_power,
    );
    close(
        result.power().balance_residual_w,
        0.0,
        result.power().absorbed_disturbance_power_w,
    );
}

#[test]
fn coherent_drives_retain_interference_in_field_and_terminal_work() {
    let gas = gas();
    let duct =
        UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), 0.3 * gas.sound_speed).unwrap();
    let reference = duct
        .characteristic_flow_impedance()
        .unwrap()
        .value_j_s_kg2()
        .re;
    let prepared = duct.prepare_terminations(loads(reference), 1e8).unwrap();
    let h = [C64::new(0.7, 0.3), C64::new(-0.4, 0.6)];
    let left = prepared.solve([h[0], C64::ZERO]).unwrap();
    let right = prepared.solve([C64::ZERO, h[1]]).unwrap();
    let together = prepared.solve(h).unwrap();
    for x in [0.0, 0.19, 0.47, duct.length_m()] {
        let a = left.at_position(x).unwrap();
        let b = right.at_position(x).unwrap();
        let sum = together.at_position(x).unwrap();
        complex_close(
            sum.pressure_pa,
            a.pressure_pa + b.pressure_pa,
            a.pressure_pa.abs() + b.pressure_pa.abs(),
        );
        complex_close(
            sum.particle_velocity_m_s,
            a.particle_velocity_m_s + b.particle_velocity_m_s,
            a.particle_velocity_m_s.abs() + b.particle_velocity_m_s.abs(),
        );
    }
    let independent_power =
        left.power().absorbed_disturbance_power_w + right.power().absorbed_disturbance_power_w;
    assert!(
        (together.power().absorbed_disturbance_power_w - independent_power).abs()
            > 0.02 * independent_power,
        "coherent load work must retain a live cross term"
    );
    close(together.power().balance_residual_w, 0.0, independent_power);
}

#[test]
fn flow_reversal_and_end_exchange_preserve_fields_and_physical_work() {
    let gas = gas();
    let forward =
        UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), 0.42 * gas.sound_speed).unwrap();
    let backward =
        UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), -0.42 * gas.sound_speed).unwrap();
    let reference = forward
        .characteristic_flow_impedance()
        .unwrap()
        .value_j_s_kg2()
        .re;
    let z = loads(reference);
    let h = [C64::new(0.8, -0.1), C64::new(-0.3, 0.7)];
    let original = forward
        .prepare_terminations(z, 1e8)
        .unwrap()
        .solve(h)
        .unwrap();
    let reversed = backward
        .prepare_terminations([z[1], z[0]], 1e8)
        .unwrap()
        .solve([h[1], h[0]])
        .unwrap();
    for x in [0.0, 0.17, 0.52, forward.length_m()] {
        let a = original.at_position(x).unwrap();
        let b = reversed.at_position(forward.length_m() - x).unwrap();
        complex_close(a.pressure_pa, b.pressure_pa, a.pressure_pa.abs());
        complex_close(
            a.particle_velocity_m_s,
            -b.particle_velocity_m_s,
            a.particle_velocity_m_s.abs(),
        );
        complex_close(
            a.mass_flow_perturbation_kg_s,
            -b.mass_flow_perturbation_kg_s,
            a.mass_flow_perturbation_kg_s.abs(),
        );
        complex_close(
            a.stagnation_enthalpy_perturbation_j_kg,
            b.stagnation_enthalpy_perturbation_j_kg,
            a.stagnation_enthalpy_perturbation_j_kg.abs(),
        );
        close(
            a.disturbance_power_w,
            -b.disturbance_power_w,
            original.power().absorbed_disturbance_power_w,
        );
    }
    for i in 0..2 {
        close(
            original.terminals()[i].absorbed_disturbance_power_w,
            reversed.terminals()[1 - i].absorbed_disturbance_power_w,
            original.power().absorbed_disturbance_power_w,
        );
        close(
            original.terminals()[i].source_work_into_duct_w,
            reversed.terminals()[1 - i].source_work_into_duct_w,
            original.power().absorbed_disturbance_power_w,
        );
    }
}

#[test]
fn round_trip_resonance_shift_and_declared_condition_limit_are_enforced() {
    let gas = gas();
    let mach = 0.4;
    let length = 0.83;
    // H'=0 at both ends gives r_left*r_right=1. The exact round-trip
    // phase is 2 omega L/[c(1-M²)], so omega_n=n pi c(1-M²)/L.
    let resonance = core::f64::consts::PI * gas.sound_speed * (1.0 - mach * mach) / length;
    let zero = FlowImpedance::try_new(C64::ZERO).unwrap();
    let resonant =
        UniformFlowCylinder::try_new(&gas, &cylinder(), resonance, mach * gas.sound_speed).unwrap();
    assert!(matches!(
        resonant.prepare_terminations([zero; 2], 1e8),
        Err(DuctError::Singular)
    ));
    let nearby = UniformFlowCylinder::try_new(
        &gas,
        &cylinder(),
        resonance * 1.0001,
        mach * gas.sound_speed,
    )
    .unwrap();
    let admitted = nearby.prepare_terminations([zero; 2], 1e5).unwrap();
    assert!(admitted.condition_estimate() > 1e3);
    assert!(matches!(
        nearby.prepare_terminations([zero; 2], 1e3),
        Err(DuctError::Singular)
    ));
    for limit in [0.0, 0.9, f64::NAN, f64::INFINITY] {
        assert!(nearby.prepare_terminations([zero; 2], limit).is_err());
    }
    for impedance in [
        C64::new(-0.1, 0.0),
        C64::new(0.1, f64::NAN),
        C64::new(f64::INFINITY, 0.0),
    ] {
        assert!(FlowImpedance::try_new(impedance).is_err());
    }
}

#[test]
fn prepared_loads_retain_original_owner_and_reject_invalid_drives() {
    #[derive(Debug)]
    struct Counted {
        state: FluidAcousticState,
        calls: Cell<usize>,
    }
    impl FluidAcoustics for Counted {
        fn acoustic_state(&self) -> Result<FluidAcousticState, fs_material::MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.state)
        }
    }
    let source = Counted {
        state: gas().acoustic_state().unwrap(),
        calls: Cell::new(0),
    };
    let duct = UniformFlowCylinder::try_new(&source, &cylinder(), omega(), 60.0).unwrap();
    let matched = duct.characteristic_flow_impedance().unwrap();
    let prepared = duct.prepare_terminations([matched; 2], 1e8).unwrap();
    assert!(core::ptr::eq(prepared.duct().source(), &source));
    assert_eq!(prepared.impedances(), &[matched; 2]);
    assert_eq!(prepared.enthalpy_reflection(), [C64::ZERO; 2]);
    for i in 0..20 {
        let result = prepared
            .solve([C64::new(f64::from(i) / 40.0, -0.2), C64::new(0.3, 0.1)])
            .unwrap();
        let first = result.at_position(0.21).unwrap();
        result.at_position(0.63).unwrap();
        assert_eq!(result.at_position(0.21).unwrap(), first);
    }
    assert_eq!(source.calls.get(), 1);
    let silence = prepared.solve([C64::ZERO; 2]).unwrap();
    assert_eq!(silence.at_position(0.4).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(silence.power().absorbed_disturbance_power_w, 0.0);
    assert_eq!(silence.power().source_work_into_duct_w, 0.0);
    for source in [
        C64::new(f64::NAN, 0.0),
        C64::new(0.0, f64::INFINITY),
        C64::from_re(f64::MAX),
    ] {
        assert!(prepared.solve([source, C64::ZERO]).is_err());
    }
    // A finite load with an otherwise overflowing Z+R intermediate is admitted
    // through the scaled reflection calculation, without substituting infinity.
    let huge = FlowImpedance::try_new(C64::new(f64::MAX, f64::MAX)).unwrap();
    let large = duct.prepare_terminations([huge, matched], 1e8).unwrap();
    assert!(large.enthalpy_reflection().into_iter().all(finite));
    assert_eq!(source.calls.get(), 1);
}
