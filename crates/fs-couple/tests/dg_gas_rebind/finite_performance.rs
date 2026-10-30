//! Authored heat reaches the original finite gas, source and retained DG histories.
#[path = "finite_performance/transport.rs"]
mod transport;

use super::*;
use super::finite_gas::{channel_load, controls, energy_solve, gas};
use fs_couple::bernoulli_aperture::performance::{
    ApertureObservation, AperturePerformance, AperturePerformanceConfig, CoupledAperture,
    finite_gas::{FiniteGasHeatProgram, GasHeatControl},
};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::RenderError;
use fs_exec::CancelGate;
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};
use std::sync::Arc;

fn config(block: usize) -> AperturePerformanceConfig {
    AperturePerformanceConfig { sample_rate_hz: 96_000, samples: 48, max_block: block,
        max_compile_work: 64, max_controls: 8 }
}
fn source(cold: &AxialGasField) -> ApertureDg {
    let load = channel_load(cold);
    ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
}
fn performer(cold: &AxialGasField, config: AperturePerformanceConfig) -> AperturePerformance {
    let phrase = GestureSchedule::try_new(1, vec![GestureTrack {
        id: "mouth".into(), target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(30.0), events: vec![],
    }]).unwrap();
    AperturePerformance::new(CoupledAperture::Dg(source(cold)),
        ApertureObservation::Inlet, phrase, config).unwrap()
}
fn program() -> FiniteGasHeatProgram {
    FiniteGasHeatProgram {
        controls: [(0, 0.0), (7, 6.0), (19, -3.0), (31, 0.0)].into_iter()
            .map(|(sample, power_w)| GasHeatControl { sample, power_w }).collect(),
        caloric: energy_solve(), coupling: controls(),
    }
}
fn retained(player: &AperturePerformance) -> &ApertureDg {
    let CoupledAperture::Dg(source) = player.system() else { unreachable!() };
    source
}
fn same_state(actual: &ApertureDg, expected: &ApertureDg) {
    assert_eq!(actual.aperture().state(), expected.aperture().state());
    assert_eq!(actual.aperture().spec(), expected.aperture().spec());
    assert_eq!(actual.aperture().accepted_steps(), expected.aperture().accepted_steps());
    assert_eq!(actual.duct_state(), expected.duct_state());
    assert_eq!(actual.load().duct().gas_field(), expected.load().duct().gas_field());
    assert_eq!(actual.load().terminal().unwrap().pole_states(),
        expected.load().terminal().unwrap().pole_states());
    assert_eq!(actual.load().last_terminal_channel_frames(),
        expected.load().last_terminal_channel_frames());
    assert_eq!(actual.load().terminal_channels()[0].channel().unwrap().velocities_m_s(),
        expected.load().terminal_channels()[0].channel().unwrap().velocities_m_s());
    assert_eq!(actual.stored_energy_j(), expected.stored_energy_j());
}

#[test]
fn scheduled_heat_matches_original_finite_gas_transaction_at_every_source_sample() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let mut direct = source(&cold);
    let mut inventory = gas(&cold);
    let initial = inventory.clone();
    let heat = program();
    let gate = Arc::new(CancelGate::new());
    let mut player = performer(&cold, config(1))
        .with_finite_gas_heat(gas(&cold), heat.clone(), gate.clone()).unwrap();
    assert_eq!(player.finite_gas_heat_program(), Some(&heat));
    assert!(Arc::ptr_eq(player.finite_gas_cancel_gate().unwrap(), &gate));
    let mut absolute_work = 0.0;
    for sample in 0..48 {
        let power = heat.controls.iter().rev().find(|event| event.sample <= sample).unwrap().power_w;
        let expected = direct.step_with_finite_gas_heat_controlled(&mut inventory,
            TubeDrive { upstream_pressure_pa: 30.0, body_flow_m3_s: 0.0 }, 0.0,
            power * DT, energy_solve(), controls(), &mut || false).unwrap();
        let mut pressure = [f64::NAN];
        player.block(&mut pressure).unwrap();
        assert_eq!(pressure[0].to_bits(), expected.coupled.aperture.bore_pressure_pa.to_bits());
        let actual = player.last_finite_gas_step().unwrap();
        assert_eq!(actual.coupled, expected.coupled);
        assert_eq!(actual.gas.refresh, expected.gas.refresh);
        assert_eq!(actual.gas.background, expected.gas.background);
        assert_eq!(actual.gas.supplied_heat_j, power * DT);
        assert_eq!(actual.gas.transported_enthalpy_j, 0.0);
        assert!(actual.gas.work_residual_j.abs() <= controls().work_tolerance_j);
        assert!(actual.gas.gas_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(actual.combined_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(actual.coupled.balance_residual_j().abs() <= controls().energy_tolerance_j);
        absolute_work += actual.gas.refresh.parameter_work_j.abs();
        assert_eq!(player.finite_gas(), Some(&inventory));
        assert_eq!(inventory.mass_kg(), initial.mass_kg());
        assert_eq!(inventory.volume_m3(), initial.volume_m3());
        assert_eq!(direct.load().duct().gas_field().unwrap().uniform_state(), Some(inventory.state()));
        assert_eq!(direct.load().terminal_channels()[0].gas(), inventory.state().state());
        assert_eq!(direct.aperture().spec().density_kg_m3, inventory.state().state().density);
        assert_eq!(direct.aperture().spec().impedance_pa_s_m3, direct.load().impedance_pa_s_m3());
        same_state(retained(&player), &direct);
        let accepted = heat.controls.iter().take_while(|event| event.sample <= sample).count();
        assert_eq!(player.applied_heat_controls(), &heat.controls[..accepted]);
    }
    assert!(absolute_work > 0.0);
    assert!(inventory.state().state().temperature > 300.0);
    assert_eq!(player.samples_rendered(), 48);
}

#[test]
fn heat_playback_is_block_invariant_and_zero_power_preserves_the_frozen_solution() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let mut reference: Option<([f64; 48], AperturePerformance)> = None;
    for block in [1, 7, 48] {
        let mut player = performer(&cold, config(block)).with_finite_gas_heat(
            gas(&cold), program(), Arc::new(CancelGate::new())).unwrap();
        let mut pressure = [0.0; 48];
        for chunk in pressure.chunks_mut(block) { player.block(chunk).unwrap(); }
        if let Some((values, previous)) = &reference {
            assert_eq!(&pressure, values);
            assert_eq!(player.finite_gas(), previous.finite_gas());
            assert_eq!(player.applied_heat_controls(), previous.applied_heat_controls());
            same_state(retained(&player), retained(previous));
        }
        reference = Some((pressure, player));
    }
    let mut zero = program();
    zero.controls.truncate(1);
    let mut finite = performer(&cold, config(48)).with_finite_gas_heat(
        gas(&cold), zero, Arc::new(CancelGate::new())).unwrap();
    let mut frozen = performer(&cold, config(48));
    let (mut heated, mut ordinary) = ([0.0; 48], [0.0; 48]);
    finite.block(&mut heated).unwrap();
    frozen.block(&mut ordinary).unwrap();
    for (&actual, &expected) in heated.iter().zip(&ordinary) { close(actual, expected, expected.abs()); }
    for (&actual, &expected) in retained(&finite).duct_state().pressure().iter()
        .zip(retained(&frozen).duct_state().pressure()) { close(actual, expected, expected.abs()); }
    close(retained(&finite).stored_energy_j(), retained(&frozen).stored_energy_j(),
        retained(&frozen).stored_energy_j());
    assert!((finite.finite_gas().unwrap().state().state().temperature - 300.0).abs() < 2e-9);
    assert!(finite.last_finite_gas_step().unwrap().gas.refresh.parameter_work_j.abs() < 1e-15);
    assert_ne!(reference.unwrap().0, ordinary, "heat must change the actual continuing pressure");
}

#[test]
fn refused_heat_and_cancellation_preserve_the_accepted_sample_and_control_cursors() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let mut reference = performer(&cold, config(7)).with_finite_gas_heat(
        gas(&cold), program(), Arc::new(CancelGate::new())).unwrap();
    reference.block(&mut [0.0; 7]).unwrap();
    for cancel in [false, true] {
        let mut heat = program();
        if !cancel { heat.coupling.maximum_temperature_increment_k = 0.001; }
        let gate = Arc::new(CancelGate::new());
        let mut attempted = performer(&cold, config(7))
            .with_finite_gas_heat(gas(&cold), heat, gate.clone()).unwrap();
        attempted.block(&mut [0.0; 7]).unwrap();
        if cancel { gate.request(); }
        let mut output = [123.0];
        let error = attempted.block(&mut output).unwrap_err();
        if cancel { assert!(matches!(error, RenderError::Voice(
            fs_couple::acoustic_realize::AcousticRealizeError::Cancelled))); }
        assert_eq!(output, [123.0]);
        assert_eq!(attempted.samples_rendered(), 7);
        assert_eq!(attempted.finite_gas(), reference.finite_gas());
        assert_eq!(attempted.applied_heat_controls(), reference.applied_heat_controls());
        assert_eq!(attempted.applied_controls(), reference.applied_controls());
        assert_eq!(attempted.last_finite_gas_step().unwrap().coupled,
            reference.last_finite_gas_step().unwrap().coupled);
        same_state(retained(&attempted), retained(&reference));
        assert!(matches!(attempted.block(&mut output), Err(RenderError::Poisoned)));
    }
}

#[test]
fn heat_admission_enforces_original_clock_shared_storage_work_and_solver_domains() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let attach = |config, heat| performer(&cold, config).with_finite_gas_heat(
        gas(&cold), heat, Arc::new(CancelGate::new()));
    // One pressure assignment and four heat assignments share five slots;
    // pressure compilation costs two, independently of four heat validations.
    assert!(attach(AperturePerformanceConfig { max_controls: 5, max_compile_work: 4,
        ..config(48) }, program()).is_ok());
    for limited in [AperturePerformanceConfig { max_controls: 4, ..config(48) },
        AperturePerformanceConfig { max_compile_work: 3, ..config(48) }] {
        assert!(attach(limited, program()).is_err());
    }
    for bad in [vec![], vec![(1, 0.0)], vec![(0, 0.0), (0, 1.0)],
        vec![(0, 0.0), (48, 1.0)], vec![(0, f64::NAN)],
        vec![(0, f64::INFINITY)], vec![(0, f64::from_bits(1))]] {
        let mut heat = program();
        heat.controls = bad.into_iter().map(|(sample, power_w)| GasHeatControl { sample, power_w }).collect();
        assert!(attach(config(48), heat).is_err());
    }
    let mut heat = program();
    heat.caloric.minimum_temperature_k = 310.0;
    assert!(attach(config(48), heat).is_err());
    let mut heat = program();
    heat.coupling.maximum_iterations = 0;
    assert!(attach(config(48), heat).is_err());
    let gate = Arc::new(CancelGate::new());
    gate.request();
    assert!(performer(&cold, config(48)).with_finite_gas_heat(gas(&cold), program(), gate).is_err());
    let already = attach(config(48), program()).unwrap();
    assert!(already.with_finite_gas_heat(gas(&cold), program(), Arc::new(CancelGate::new())).is_err());
    let mut advanced = performer(&cold, config(48));
    advanced.block(&mut [0.0]).unwrap();
    assert!(advanced.with_finite_gas_heat(gas(&cold), program(), Arc::new(CancelGate::new())).is_err());
}

#[path = "finite_performance/wall.rs"]
mod wall;
