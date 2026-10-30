//! One shared finite inventory pays the work of both original loaded tracts.
use super::finite_gas::{channel_load, compression_energy, controls, energy_solve};
use super::*;
use fs_couple::dg_port::finite_gas::{FiniteGasDgBackground, FiniteGasDgError};
use fs_material::thermochemical_gas::energy::{GasEnergyError, parcel::ClosedGasParcel};

const UP_LENGTH: f64 = 0.035;
const GAP_VOLUME: f64 = 0.002 * 4e-5 * 0.001;
const HEAT: f64 = 0.003;
const FORCES: [f64; 2] = [2e-5, -3e-5];

fn volume() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS * (UP_LENGTH + LENGTH) + 2.0 * GAP_VOLUME
}

fn warm() -> (LoadedApertureDg, ClosedGasParcel) {
    let cold = [
        field(300.0, 1.0, 101325.0, UP_LENGTH),
        field(300.0, 1.0, 101325.0, LENGTH),
    ];
    let up = channel_load(&cold[0]);
    let down = channel_load(&cold[1]);
    let (_, impedance) = series_inlet_state(&up, &down).unwrap();
    let mut model =
        LoadedApertureDg::new_loaded_with_upstream(paired_source(&down, impedance), up, down)
            .unwrap();
    let base = super::finite_gas::gas(&cold[1]);
    let parcel =
        ClosedGasParcel::new(base.model().clone(), base.point().clone(), volume()).unwrap();
    for sample in 0..24 {
        model.step(drive(sample)).unwrap();
    }
    (model, parcel)
}

fn loads(model: &LoadedApertureDg) -> [&DgPortLoad; 2] {
    [model.upstream_load().unwrap(), model.load()]
}

// Actual quadratic bore fields and uniform-volume FV channel cells, independent
// of the production mass matrices and complete-load storage accumulator.
fn represented_energy(model: &LoadedApertureDg) -> f64 {
    pair(model).stored_energy_j()
        + loads(model)
            .iter()
            .map(|load| {
                let channel = load.terminal_channels()[0].channel().unwrap();
                wave_energy(
                    load.state(),
                    load.duct().inlet_gas(),
                    load.duct().total_length(),
                ) + 0.5 * channel.gas().density * GAP_VOLUME / channel.velocities_m_s().len() as f64
                    * channel.velocities_m_s().iter().map(|v| v * v).sum::<f64>()
            })
            .sum::<f64>()
}

fn same_load(actual: &DgPortLoad, expected: &DgPortLoad) {
    assert_eq!(actual.state(), expected.state());
    assert_eq!(actual.duct().gas_field(), expected.duct().gas_field());
    assert_eq!(
        actual.duct().operator_identity(),
        expected.duct().operator_identity()
    );
    assert_eq!(
        actual.port_response().unwrap(),
        expected.port_response().unwrap()
    );
    assert_eq!(
        actual.terminal().unwrap().pole_states(),
        expected.terminal().unwrap().pole_states()
    );
    assert_eq!(
        actual.last_terminal_frames(),
        expected.last_terminal_frames()
    );
    assert_eq!(
        actual.last_terminal_channel_frames(),
        expected.last_terminal_channel_frames()
    );
    assert_eq!(
        actual.cumulative_parameter_work_j(),
        expected.cumulative_parameter_work_j()
    );
    let (a, b) = (
        &actual.terminal_channels()[0],
        &expected.terminal_channels()[0],
    );
    assert_eq!(a.gas(), b.gas());
    assert_eq!(
        a.channel().unwrap().velocities_m_s(),
        b.channel().unwrap().velocities_m_s()
    );
    assert_eq!(
        a.channel().unwrap().last_step(),
        b.channel().unwrap().last_step()
    );
    assert_eq!(
        a.cumulative_parameter_work_j(),
        b.cumulative_parameter_work_j()
    );
    assert_eq!(
        a.cumulative_pressure_work_j(),
        b.cumulative_pressure_work_j()
    );
    assert_eq!(
        a.cumulative_viscous_dissipation_j(),
        b.cumulative_viscous_dissipation_j()
    );
    assert_eq!(
        a.cumulative_numerical_dissipation_j(),
        b.cumulative_numerical_dissipation_j()
    );
}

#[test]
fn g1_shared_heating_pays_both_bore_work_accounts_and_preserves_both_blades() {
    let (mut live, mut gas) = warm();
    let (mut frozen, _) = warm();
    let old = gas.clone();
    let initial_energy = represented_energy(&live);
    let reference = frozen.step_with_forces(drive(24), FORCES).unwrap();
    let compression = loads(&frozen).map(|load| {
        compression_energy(load.state(), old.state().state()) * load.duct().total_length() / LENGTH
    });
    assert!(compression.iter().all(|energy| *energy > 0.0));
    let result = live
        .step_with_shared_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            FORCES,
            [0.0; 2],
            HEAT,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    let r = old.state().state().specific_gas_constant;
    let rho = old.state().state().density;
    let modulus = |t: f64| rho * r * t * (3.5 + 0.001 * t) / (2.5 + 0.001 * t);
    let caloric =
        |t: f64| old.mass_kg() * r * (2.5 * (t - 300.0) + 0.0005 * (t * t - 300.0 * 300.0));
    let balance = |t: f64| {
        caloric(t) + compression.iter().sum::<f64>() * (modulus(t) / modulus(300.0) - 1.0) - HEAT
    };
    let (mut lo, mut hi) = (300.0, 380.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if balance(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let temperature = gas.state().state().temperature;
    assert!((temperature - 0.5 * (lo + hi)).abs() < 2e-7);
    assert_eq!(gas.volume_m3(), volume());
    assert_eq!(gas.mass_kg(), old.mass_kg());
    let mut work = 0.0;
    let mut exterior_work = 0.0;
    for (i, (next, previous)) in loads(&live).into_iter().zip(loads(&frozen)).enumerate() {
        let expected_work = compression[i] * (modulus(temperature) / modulus(300.0) - 1.0);
        assert!((result.refresh[i].parameter_work_j - expected_work).abs() < 2e-14);
        work += result.refresh[i].parameter_work_j;
        assert_eq!(
            next.duct().gas_field().unwrap().uniform_state(),
            Some(gas.state())
        );
        assert_eq!(next.terminal_channels()[0].gas(), gas.state().state());
        assert_eq!(next.state().step_count(), 25);
        assert_eq!(
            next.terminal().unwrap().pole_states(),
            previous.terminal().unwrap().pole_states()
        );
        assert_eq!(next.last_terminal_frames(), previous.last_terminal_frames());
        assert_eq!(
            next.last_terminal_channel_frames(),
            previous.last_terminal_channel_frames()
        );
        exterior_work += next
            .last_terminal_frames()
            .iter()
            .map(|frame| frame.load_work_j)
            .sum::<f64>();
    }
    assert_eq!(result.coupled.aperture, reference.aperture);
    assert_eq!(pair(&live).blade_states(), pair(&frozen).blade_states());
    assert_ne!(pair(&live).blade_states()[0], pair(&live).blade_states()[1]);
    assert_eq!(pair(&live).gas(), Some(gas.state().state()));
    assert_eq!(
        pair(&live).spec().impedance_pa_s_m3,
        live.upstream_load().unwrap().impedance_pa_s_m3() + live.load().impedance_pa_s_m3()
    );
    assert_eq!(result.coupled.parameter_work_j, work);
    assert!((result.work_by_gas_j - work).abs() <= controls().work_tolerance_j);
    close(
        represented_energy(&live),
        live.stored_energy_j(),
        live.stored_energy_j(),
    );
    let residual = caloric(temperature) + represented_energy(&live) - initial_energy
        + result.coupled.inferred_total_dissipation_j
        + exterior_work
        - result.coupled.upstream_work_j
        - result.coupled.body_work_j
        - result.coupled.aperture.mechanical_work_j
        - HEAT;
    assert!(
        residual.abs() <= controls().energy_tolerance_j,
        "{residual:e}"
    );
    assert!(result.combined_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert_eq!(
        result.coupled.upstream.as_ref().unwrap().inlet_work_j,
        reference.upstream.as_ref().unwrap().inlet_work_j
    );
    let next = live.step(drive(25)).unwrap();
    assert_eq!(next.parameter_work_j, 0.0);
    assert_eq!(
        next.duct.ledger.cumulative_parameter_work_j,
        result.refresh[1].parameter_work_j
    );
    assert_eq!(
        next.upstream
            .unwrap()
            .duct
            .ledger
            .cumulative_parameter_work_j,
        result.refresh[0].parameter_work_j
    );
}

#[test]
fn g3_shared_inventory_admits_actual_summed_volume_and_one_complete_source() {
    let (mut model, gas) = warm();
    let old_loads = loads(&model).map(Clone::clone);
    let motion = *pair(&model).blade_states();
    for (model_source, inventory_volume) in [
        (gas.model().clone(), volume() - 2.0 * GAP_VOLUME),
        (
            gas.model().clone(),
            core::f64::consts::PI * RADIUS * RADIUS * LENGTH + GAP_VOLUME,
        ),
        (
            gas.model().with_composition(composition(0.5)).unwrap(),
            volume(),
        ),
    ] {
        let mut wrong =
            ClosedGasParcel::new(model_source, gas.point().clone(), inventory_volume).unwrap();
        let before = wrong.clone();
        let mut called = false;
        let result = model.step_with_shared_finite_gas_controlled::<(), GasEnergyError, _>(
            &mut wrong,
            drive(24),
            FORCES,
            [0.0; 2],
            controls(),
            |_, _, _, _| {
                called = true;
                Err(GasEnergyError::Invalid(
                    "must not evaluate an unowned inventory",
                ))
            },
            &mut || false,
        );
        assert!(matches!(result, Err(FiniteGasDgError::Invalid(_))));
        assert!(!called);
        assert_eq!(wrong, before);
        assert_eq!(pair(&model).blade_states(), &motion);
        for (load, before) in loads(&model).into_iter().zip(&old_loads) {
            same_load(load, before);
        }
    }
}

#[test]
fn g4_shared_source_refusal_and_final_cancellation_leave_both_tracts_retryable() {
    let (mut model, mut gas) = warm();
    let (mut clean, mut expected_gas) = warm();
    let original_gas = gas.clone();
    let old_loads = loads(&model).map(Clone::clone);
    let motion = *pair(&model).blade_states();
    let spec = *pair(&model).spec();
    let viscosity = pair(&model).viscous_slit();
    let unchanged = |model: &LoadedApertureDg, gas: &ClosedGasParcel| {
        assert_eq!(gas, &original_gas);
        assert_eq!(pair(model).blade_states(), &motion);
        assert_eq!(pair(model).spec(), &spec);
        assert_eq!(pair(model).gas(), Some(original_gas.state().state()));
        assert_eq!(pair(model).viscous_slit(), viscosity);
        assert_eq!(model.aperture().accepted_steps(), 24);
        for (load, old) in loads(model).into_iter().zip(&old_loads) {
            same_load(load, old);
        }
    };
    let failed = model.step_with_shared_finite_gas_controlled(
        &mut gas,
        drive(24),
        FORCES,
        [0.0; 2],
        controls(),
        |_, _, _, _| {
            Err::<FiniteGasDgBackground<()>, _>(GasEnergyError::Invalid(
                "source refused after both acoustic previews",
            ))
        },
        &mut || false,
    );
    assert!(matches!(failed, Err(FiniteGasDgError::Background(_))));
    unchanged(&model, &gas);
    let mut polls = 0;
    let expected = clean
        .step_with_shared_finite_gas_heat_controlled(
            &mut expected_gas,
            drive(24),
            FORCES,
            [0.0; 2],
            HEAT,
            energy_solve(),
            controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        model.step_with_shared_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            FORCES,
            [0.0; 2],
            HEAT,
            energy_solve(),
            controls(),
            &mut || {
                calls += 1;
                calls == polls
            },
        ),
        Err(FiniteGasDgError::Cancelled)
    ));
    assert_eq!(calls, polls);
    unchanged(&model, &gas);
    let retry = model
        .step_with_shared_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            FORCES,
            [0.0; 2],
            HEAT,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.coupled, expected.coupled);
    assert_eq!(retry.refresh, expected.refresh);
    assert_eq!(gas, expected_gas);
    for sample in 25..32 {
        assert_eq!(
            model.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
}
