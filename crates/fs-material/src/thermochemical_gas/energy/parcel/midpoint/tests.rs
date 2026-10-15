//! Public constitutive steps with synthetic source laws, not empirical gas claims.
use super::*;
use crate::fluid_transport::FluidTransport;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn parcel(x: f64, slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], slope),
        point(400.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.001).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn combined_heat_and_volume_follow_the_closed_form_constant_capacity_equation() {
    let old = parcel(0.75, 0.0);
    let medium = old.state().state();
    let cv = old.state().thermodynamics().properties().mass_specific_cv().value();
    for volume in [0.0006, 0.0009, 0.001, 0.0014, 0.002] {
        for heat in [-10.0, 0.0, 30.0] {
            let trial = old.midpoint_heat_volume(volume, heat, controls(), || false).unwrap();
            let vm = 0.5*(old.volume_m3()+volume);
            let ratio = (volume-old.volume_m3())/vm;
            let tm = (cv*400.0+heat/(2.0*old.mass_kg()))/(cv+0.5*ratio*medium.specific_gas_constant);
            near(trial.midpoint().state().temperature, tm, 2e-11);
            near(trial.endpoint().state().state().temperature, 2.0*tm-400.0, 2e-11);
            let work = old.mass_kg()*medium.specific_gas_constant*tm*ratio;
            near(trial.step().work_by_gas_j, work, 2e-10);
            near(trial.step().stored_energy_change_j+work, heat, 2e-9);
            assert_eq!(trial.endpoint().mass_kg(), old.mass_kg());
            assert_eq!(trial.endpoint().point().axes().get("material-lot"), Some(&42.0));
            assert!(trial.step().energy_residual_j.abs() <= old.mass_kg()*controls().energy_tolerance_j_kg);
        }
    }
}

#[test]
fn nonlinear_heat_and_composition_update_the_same_mechanical_thermal_acoustic_states() {
    let old = parcel(0.75, 0.002);
    let heated = old.midpoint_heat_volume(0.0007, 40.0, controls(), || false).unwrap();
    let final_state = heated.endpoint().state();
    let tm = heated.midpoint().state().temperature;
    let t1 = final_state.state().temperature;
    // Midpoint is in INTERNAL ENERGY, not in Kelvin for nonlinear heat capacity.
    assert!((tm-0.5*(400.0+t1)).abs() > 0.05);
    near(2.0*(heated.midpoint().specific_internal_energy_j_kg()-old.state().specific_internal_energy_j_kg()),
        final_state.specific_internal_energy_j_kg()-old.state().specific_internal_energy_j_kg(), 1e-10);
    let gas = final_state.state();
    assert!(gas.sound_speed > old.state().state().sound_speed);
    assert!(gas.dynamic_viscosity > old.state().state().dynamic_viscosity);
    assert!(gas.thermal_conductivity > old.state().state().thermal_conductivity);
    near(final_state.transport_state().unwrap().density_kg_m3(), old.mass_kg()/0.0007, 1e-12);
    near(gas.sound_speed.powi(2), gas.gamma*gas.pressure/gas.density, 1e-12);
    assert_ne!(final_state.components(), old.state().components());
    let other = parcel(0.25, 0.002).midpoint_heat_volume(0.0007, 40.0, controls(), || false).unwrap();
    assert_ne!(other.midpoint().state().pressure, heated.midpoint().state().pressure);
    assert_ne!(other.endpoint().state().state().sound_speed, gas.sound_speed);
}

#[test]
fn pressure_derivatives_match_fresh_nonlinear_steps_in_both_control_directions() {
    for (x, slope) in [(0.75, 0.0), (0.75, 0.002), (0.25, 0.001)] {
        let old = parcel(x, slope);
        for volume in [0.0007, 0.001, 0.0015] {
            let trial = old.midpoint_heat_volume(volume, 20.0, controls(), || false).unwrap();
            let eval = |v, q| old.midpoint_heat_volume(v, q, controls(), || false).unwrap().midpoint().state().pressure;
            let dv = volume*1e-4;
            near(trial.step().pressure_final_volume_derivative_pa_m3,
                (eval(volume+dv, 20.0)-eval(volume-dv, 20.0))/(2.0*dv), 1e-7);
            let dq = 0.001;
            near(trial.step().pressure_heat_derivative_pa_j,
                (eval(volume, 20.0+dq)-eval(volume, 20.0-dq))/(2.0*dq), 1e-7);
        }
    }
}

#[test]
fn fixed_volume_preserves_existing_heating_behavior_and_no_op_is_bitwise() {
    let old = parcel(0.75, 0.001);
    let zero = old.midpoint_heat_volume(old.volume_m3(), 0.0, controls(), || false).unwrap();
    assert_eq!(zero.endpoint(), &old);
    assert_eq!(zero.midpoint(), old.state());
    assert_eq!(zero.step().iterations, 0);
    for heat in [-10.0, 20.0] {
        let mut reference = old.clone();
        reference.add_heat(heat, controls(), || false).unwrap();
        let trial = old.midpoint_heat_volume(old.volume_m3(), heat, controls(), || false).unwrap();
        near(trial.endpoint().state().state().temperature, reference.state().state().temperature, 1e-10);
        near(trial.endpoint().state().state().pressure, reference.state().state().pressure, 1e-10);
        assert_eq!(trial.step().work_by_gas_j, 0.0);
    }
}

#[test]
fn adiabatic_refinement_converges_without_claiming_exact_finite_step_isentropy() {
    let old = parcel(0.75, 0.0);
    let exact = 400.0*2.0_f64.powf(old.state().state().gamma-1.0);
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut gas = old.clone();
        let mut work = 0.0;
        for i in 1..=count {
            let v = 0.001-0.0005*f64::from(i)/f64::from(count);
            let trial = gas.midpoint_heat_volume(v, 0.0, controls(), || false).unwrap();
            work += trial.step().work_by_gas_j;
            gas = trial.into_parts().0;
        }
        errors.push((gas.state().state().temperature-exact).abs());
        near(gas.mass_kg()*(gas.state().specific_internal_energy_j_kg()-old.state().specific_internal_energy_j_kg())+work,
            0.0, 1e-8);
    }
    assert!(errors[0]/errors[1] > 3.8);
    assert!(errors[1]/errors[2] > 3.8);
}

#[test]
fn failed_midpoint_endpoint_and_every_cancellation_boundary_leave_input_unchanged() {
    let old = parcel(0.75, 0.0);
    let saved = old.clone();
    for (v, q) in [(0.0, 0.0), (f64::NAN, 0.0), (0.001, f64::INFINITY),
        (0.001, 1e12), (0.001, -1e12), (1e-12, 0.0)] {
        assert!(old.midpoint_heat_volume(v, q, controls(), || false).is_err());
        assert_eq!(old, saved);
    }
    // The midpoint (~629 K) is supported but the endpoint (~857 K) is not.
    let narrow = GasEnergySolve { maximum_temperature_k: 700.0, ..controls() };
    assert!(old.midpoint_heat_volume(0.0001, 0.0, narrow, || false).is_err());
    let mut total_polls = 0;
    old.midpoint_heat_volume(0.0008, 20.0, controls(), || { total_polls += 1; false }).unwrap();
    for stop in 1..=total_polls {
        let mut polls = 0;
        assert!(old.midpoint_heat_volume(0.0008, 20.0, controls(), || { polls += 1; polls == stop }).is_err());
        assert_eq!(old, saved);
    }
    let tiny = GasEnergySolve { maximum_iterations: 1, ..controls() };
    assert!(parcel(0.75, 0.002).midpoint_heat_volume(0.0008, 20.0, tiny, || false).is_err());
}
