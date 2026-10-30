//! Finite thermodynamic energy reaches original moving-aperture trajectories.
use super::*;
use fs_couple::dg_port::finite_gas::{FiniteGasDgError, FiniteGasDgLoss, FiniteGasDgSolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn controls() -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        work_bracket_j: [-1e-5, 1e-5],
        work_tolerance_j: 1e-15,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 80,
        maximum_temperature_increment_k: 80.0,
        loss: FiniteGasDgLoss::Lossless,
    }
}

fn energy_solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 900.0,
        energy_tolerance_j_kg: 1e-8,
        maximum_iterations: 128,
    }
}

fn kind() -> DgTerminalChannelKind {
    DgTerminalChannelKind::Slit {
        width_m: 0.002,
        height_m: 4e-5,
        length_m: 0.001,
        transverse_cells: 8,
        maximum_mach: 0.2,
    }
}

fn volume() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS * LENGTH + 0.002 * 4e-5 * 0.001
}

fn gas(cold: &AxialGasField) -> ClosedGasParcel {
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            101_325.0,
        )
        .unwrap();
    ClosedGasParcel::new(
        cold.model().with_composition(composition(1.0)).unwrap(),
        point,
        volume(),
    )
    .unwrap()
}

fn channel_load(cold: &AxialGasField) -> DgPortLoad {
    load(cold, false)
        .with_terminal_channels(
            vec![DgTerminalChannelInput {
                port: PortId(1),
                kind: kind(),
            }],
            8,
        )
        .unwrap()
}

fn scalar(cold: &AxialGasField) -> ApertureDg {
    let load = channel_load(cold);
    ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
}

// Independent Gauss integration of the actual quadratic pressure polynomial.
fn compression_energy(state: &DgState, gas: &GasState) -> f64 {
    let node = (3.0_f64 / 5.0).sqrt();
    let p = state.pressure();
    assert_eq!(p.len(), 3);
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    [(-node, 5.0 / 9.0), (0.0, 8.0 / 9.0), (node, 5.0 / 9.0)]
        .iter()
        .map(|&(x, weight)| {
            let value =
                0.5 * x * (x - 1.0) * p[0] + (1.0 - x * x) * p[1] + 0.5 * x * (x + 1.0) * p[2];
            LENGTH * weight * area * value * value / (4.0 * bulk(gas))
        })
        .sum()
}

#[test]
fn g1_finite_heat_pays_wave_work_and_reaches_the_same_moving_source() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let (mut live, mut frozen) = (scalar(&cold), scalar(&cold));
    for sample in 0..24 {
        assert_eq!(
            live.step(drive(sample)).unwrap(),
            frozen.step(drive(sample)).unwrap()
        );
    }
    let mut gas = gas(&cold);
    let initial_gas = gas.clone();
    let initial_energy = live.stored_energy_j();
    let frozen_frame = frozen.step_with_force(drive(24), 2e-5).unwrap();
    let pressure_energy = compression_energy(frozen.duct_state(), initial_gas.state().state());
    assert!(pressure_energy > 0.0);
    let heat = 0.003;
    let result = live
        .step_with_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            2e-5,
            heat,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();

    // The synthetic NASA source has Cp/R = 3.5 + 0.001*T. At fixed mass
    // and volume, kinetic channel/wave energy is unchanged by heating;
    // compression energy changes with bulk modulus at retained compression.
    let r = initial_gas.state().state().specific_gas_constant;
    let density = initial_gas.state().state().density;
    let modulus = |t: f64| density * r * t * (3.5 + 0.001 * t) / (2.5 + 0.001 * t);
    let caloric =
        |t: f64| initial_gas.mass_kg() * r * (2.5 * (t - 300.0) + 0.0005 * (t * t - 300.0 * 300.0));
    let balance =
        |t: f64| caloric(t) + pressure_energy * (modulus(t) / modulus(300.0) - 1.0) - heat;
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
    let expected_work = pressure_energy * (modulus(temperature) / modulus(300.0) - 1.0);
    assert!((result.gas.refresh.parameter_work_j - expected_work).abs() < 2e-14);
    assert!(result.gas.refresh.parameter_work_j > 0.0);
    let mut heat_only = initial_gas.clone();
    heat_only.add_heat(heat, energy_solve(), || false).unwrap();
    assert!(temperature < heat_only.state().state().temperature);

    assert_eq!(result.coupled.aperture, frozen_frame.aperture);
    assert_eq!(live.aperture().state(), frozen.aperture().state());
    assert_eq!(live.aperture().accepted_steps(), 25);
    assert_eq!(live.duct_state().step_count(), 25);
    assert_eq!(
        live.load().terminal().unwrap().pole_states(),
        frozen.load().terminal().unwrap().pole_states()
    );
    assert_eq!(
        live.load().last_terminal_channel_frames(),
        frozen.load().last_terminal_channel_frames()
    );
    assert_eq!(
        live.aperture().spec().density_kg_m3,
        gas.state().state().density
    );
    assert_eq!(
        live.aperture().spec().impedance_pa_s_m3,
        live.load().impedance_pa_s_m3()
    );
    assert_eq!(
        live.load().duct().gas_field().unwrap().uniform_state(),
        Some(gas.state())
    );
    assert_eq!(
        live.load().terminal_channels()[0].gas(),
        gas.state().state()
    );
    assert_eq!(result.coupled.duct, result.gas.acoustic);
    assert_eq!(
        result.coupled.parameter_work_j,
        result.gas.refresh.parameter_work_j
    );
    assert_eq!(result.coupled.stored_energy_j, live.stored_energy_j());
    assert_eq!(
        result.coupled.storage_change_j,
        live.stored_energy_j() - initial_energy
    );
    assert_eq!(
        result.coupled.duct.port_outputs[0].delta_x,
        frozen_frame.duct.port_outputs[0].delta_x
    );
    assert_eq!(
        result.coupled.duct.port_outputs[0].end_conjugate,
        live.load()
            .duct()
            .port_pressures(live.duct_state())
            .unwrap()[0]
            .1
    );
    assert!(result.coupled.balance_residual_j().abs() <= controls().energy_tolerance_j);
    assert!(result.combined_energy_residual_j.abs() <= controls().energy_tolerance_j);

    let mut difference = 0.0;
    for sample in 25..40 {
        let next = live.step(drive(sample)).unwrap();
        let reference = frozen.step(drive(sample)).unwrap();
        assert_eq!(next.parameter_work_j, 0.0);
        assert_eq!(
            next.duct.ledger.cumulative_parameter_work_j,
            result.gas.refresh.parameter_work_j
        );
        difference +=
            (next.aperture.bore_pressure_pa - reference.aperture.bore_pressure_pa).powi(2);
    }
    assert!(difference > 0.0);
}

#[test]
fn g4_finite_background_failure_and_final_cancellation_preserve_mechanics_gas_and_load() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let (mut attempted, mut clean) = (scalar(&cold), scalar(&cold));
    for sample in 0..24 {
        assert_eq!(
            attempted.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
    let mut gas = gas(&cold);
    let original_gas = gas.clone();
    let state = attempted.duct_state().clone();
    let motion = attempted.aperture().state();
    let spec = *attempted.aperture().spec();
    let terminal = attempted.load().terminal().unwrap().pole_states().to_vec();
    let channels = attempted.load().last_terminal_channel_frames().to_vec();
    let velocity = attempted.load().terminal_channels()[0]
        .channel()
        .unwrap()
        .velocities_m_s()
        .to_vec();
    let unchanged = |model: &ApertureDg, gas: &ClosedGasParcel| {
        assert_eq!(gas, &original_gas);
        assert_eq!(model.duct_state(), &state);
        assert_eq!(model.aperture().state(), motion);
        assert_eq!(model.aperture().spec(), &spec);
        assert_eq!(model.aperture().accepted_steps(), 24);
        assert_eq!(model.load().terminal().unwrap().pole_states(), terminal);
        assert_eq!(model.load().last_terminal_channel_frames(), channels);
        assert_eq!(
            model.load().terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s(),
            velocity
        );
    };
    let failed = attempted.step_with_finite_gas_controlled::<(), GasEnergyError, _>(
        &mut gas,
        drive(24),
        0.0,
        controls(),
        |_, _, _, _| {
            Err(GasEnergyError::Invalid(
                "original finite background refused after acoustic preparation",
            ))
        },
        &mut || false,
    );
    assert!(matches!(failed, Err(FiniteGasDgError::Background(_))));
    unchanged(&attempted, &gas);
    assert!(
        attempted
            .step_with_finite_gas_heat_controlled(
                &mut gas,
                drive(24),
                0.0,
                0.003,
                energy_solve(),
                FiniteGasDgSolve {
                    work_bracket_j: [0.0, 0.0],
                    ..controls()
                },
                &mut || false
            )
            .is_err()
    );
    unchanged(&attempted, &gas);

    let mut clean_gas = original_gas.clone();
    let mut polls = 0;
    let expected = clean
        .step_with_finite_gas_heat_controlled(
            &mut clean_gas,
            drive(24),
            0.0,
            0.003,
            energy_solve(),
            controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        attempted.step_with_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            0.0,
            0.003,
            energy_solve(),
            controls(),
            &mut || {
                count += 1;
                count == polls
            }
        ),
        Err(FiniteGasDgError::Cancelled)
    ));
    assert_eq!(count, polls);
    unchanged(&attempted, &gas);
    let retry = attempted
        .step_with_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            0.0,
            0.003,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.coupled, expected.coupled);
    assert_eq!(retry.gas.refresh, expected.gas.refresh);
    assert_eq!(gas, clean_gas);
    for sample in 25..40 {
        assert_eq!(
            attempted.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
}

#[test]
fn g3_finite_heat_preserves_independent_blades_and_rejects_an_unowned_second_inventory() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let build = || {
        let load = channel_load(&cold);
        LoadedApertureDg::new_loaded(paired_source(&load, load.impedance_pa_s_m3()), load).unwrap()
    };
    let (mut live, mut frozen) = (build(), build());
    for sample in 0..24 {
        assert_eq!(
            live.step(drive(sample)).unwrap(),
            frozen.step(drive(sample)).unwrap()
        );
    }
    let mut gas = gas(&cold);
    let forces = [2e-5, -3e-5];
    let original = frozen.step_with_forces(drive(24), forces).unwrap();
    let result = live
        .step_with_finite_gas_heat_controlled(
            &mut gas,
            drive(24),
            forces,
            [0.0; 2],
            0.003,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.coupled.aperture, original.aperture);
    assert_eq!(pair(&live).blade_states(), pair(&frozen).blade_states());
    assert_ne!(pair(&live).blade_states()[0], pair(&live).blade_states()[1]);
    assert_eq!(pair(&live).gas(), Some(gas.state().state()));
    assert_eq!(
        pair(&live).spec().impedance_pa_s_m3,
        live.load().impedance_pa_s_m3()
    );
    assert!(result.combined_energy_residual_j.abs() <= controls().energy_tolerance_j);

    let up = load(&cold, false);
    let down = load(&cold, false);
    let (_, impedance) = series_inlet_state(&up, &down).unwrap();
    let mut two =
        LoadedApertureDg::new_loaded_with_upstream(paired_source(&down, impedance), up, down)
            .unwrap();
    let old_gas = gas.clone();
    assert!(matches!(
        two.step_with_finite_gas_heat_controlled(
            &mut gas,
            drive(0),
            [0.0; 2],
            [0.0; 2],
            0.003,
            energy_solve(),
            controls(),
            &mut || false
        ),
        Err(FiniteGasDgError::Invalid(_))
    ));
    assert_eq!(gas, old_gas);
    assert_eq!(two.aperture().accepted_steps(), 0);
    assert_eq!(two.duct_state().step_count(), 0);
    assert_eq!(two.upstream_load().unwrap().state().step_count(), 0);
}
