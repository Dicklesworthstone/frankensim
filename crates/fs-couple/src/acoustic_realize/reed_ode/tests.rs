//! G0/G1/G3: actual reed-driven Foster ducts retain complete vibrating sources
//! and close moving apertures against the current physical interval.
use super::*;
use crate::thin_plate::CompactBody;
use fs_nlmodal::SosModalStorage;

const RATE: u32 = 48_000;
const STEPS: usize = 64;

fn gas() -> GasState {
    GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        293.15,
        101_325.0,
    )
    .unwrap()
}

fn sections() -> [AcousticSection; 1] {
    [AcousticSection {
        length: 0.3,
        radius: 0.01,
        outlet_radius: 0.01,
        cells: 2,
    }]
}

fn reed(mass: f64, pressure: f64) -> BeatingReed {
    BeatingReed {
        rest_opening_m: 1e-3,
        width_m: 0.01,
        closing_pressure_pa: 1e6,
        blowing_pressure_pa: pressure,
        attack_s: 0.0,
        mass_kg: mass,
        stiffness_n_m: 100.0,
        damping_ratio: 0.01,
    }
}

fn compact(drive: f64, area: f64) -> CompactBody {
    let mut body = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 0.003,
        mass_kg: 0.02,
        frequency_hz: 400.0 / core::f64::consts::TAU,
        damping_ratio: 0.02,
    })
    .unwrap();
    body.area_m2 = area;
    body.drive_participation = drive;
    body
}

fn bank(drive: [f64; 2], area: [f64; 2], mixed: bool) -> PlateBank {
    let body = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![600.0, 900.0],
            channels: vec![],
        },
        // Positive definite, but not diagonally dominant. Its full material
        // block accompanies the actual duct's singular Foster resistance.
        &[60.0, -70.0, -70.0, 90.0],
        &drive,
        area.to_vec(),
    )
    .unwrap()
    .with_relaxation_branches(
        vec![fs_phs::RelaxationBranch {
            projection: vec![0.7, 0.0, -0.4, 0.0],
            stiffness: 1600.0,
            relaxation_time_s: 0.0008,
        }],
        1,
    )
    .unwrap();
    let mut bank = PlateBank::default();
    bank.vk = Some(body);
    if mixed {
        bank.linear.push(compact(drive[0], 0.003));
    }
    bank
}

fn render(bank: &mut PlateBank, reed: BeatingReed, n: usize) -> Vec<f64> {
    realize_reed_ode(
        &sections(),
        ChainOutlet::PressureRelease,
        &[],
        reed,
        bank,
        &gas(),
        1.0,
        RATE,
        n,
        None,
    )
    .unwrap()
}

fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    let scale = a.iter().chain(b).map(|x| x * x).sum::<f64>().sqrt();
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
        / scale.max(f64::MIN_POSITIVE)
}

#[test]
fn both_reed_selections_keep_mixed_modes_and_independent_surface_coupling() {
    for mass in [0.0, 1e-4] {
        let reed = reed(mass, 100.0);
        let response = render(&mut bank([0.0, 0.0], [0.015, -0.02], true), reed, STEPS);
        let other_drive = render(&mut bank([1.4, -0.6], [0.015, -0.02], true), reed, STEPS);
        let other_surface = render(&mut bank([0.0, 0.0], [0.015, 0.02], true), reed, STEPS);
        let vk_only = render(&mut bank([0.0, 0.0], [0.015, -0.02], false), reed, STEPS);
        let bare = render(&mut PlateBank::default(), reed, STEPS);
        assert!(response.iter().all(|pressure| pressure.is_finite()));
        assert!(relative_difference(&response, &other_drive) < 1e-9);
        assert!(relative_difference(&response, &other_surface) > 1e-5);
        assert!(relative_difference(&response, &vk_only) > 1e-4);
        assert!(relative_difference(&vk_only, &bare) > 1e-3);
    }
}

#[test]
fn no_hidden_primer_or_render_length_switch_changes_the_accepted_history() {
    for mass in [0.0, 1e-4] {
        let quiet = render(
            &mut bank([0.7, -0.4], [0.015, -0.02], true),
            reed(mass, 0.0),
            STEPS,
        );
        assert!(quiet.iter().all(|pressure| *pressure == 0.0));
        let mut warm = bank([0.7, -0.4], [0.015, -0.02], true);
        for _ in 0..24 {
            warm.advance_controlled(0.7, 1.0 / f64::from(RATE), gas().density, 1.0, || false)
                .unwrap();
        }
        let initial = warm
            .physical_image(&gas(), RATE, &fs_exec::CancelGate::default())
            .unwrap()
            .unwrap()
            .state()
            .to_vec();
        let response = render(&mut warm, reed(mass, 0.0), STEPS);
        assert!(response.iter().any(|pressure| pressure.abs() > 1e-5));
        let preserved = warm
            .physical_image(&gas(), RATE, &fs_exec::CancelGate::default())
            .unwrap()
            .unwrap();
        assert_eq!(preserved.state(), initial);
        let short = render(&mut PlateBank::default(), reed(mass, 100.0), 2);
        let long = render(&mut PlateBank::default(), reed(mass, 100.0), 8);
        assert_eq!(short, long[..2]);
    }
}

#[test]
fn quasistatic_first_window_closes_against_independent_physical_port_integration() {
    let gas = gas();
    let reed = reed(0.0, 100.0);
    let observed = render(&mut PlateBank::default(), reed, 1)[0];
    assert!(observed > 0.0 && observed < reed.blowing_pressure_pa);
    let drop = reed.blowing_pressure_pa - observed;
    let opening = quasistatic_aperture_opening(reed.rest_opening_m, reed.closing_pressure_pa, drop);
    let valve_flow = bernoulli_volume_flow(reed.width_m, opening, drop, gas.density);
    let duct = ode_bore_chain(
        &sections(),
        ChainOutlet::PressureRelease,
        1,
        &[],
        &gas,
        None,
    )
    .unwrap();
    let duration = 1.0 / f64::from(RATE);
    let (_, pressure, _) =
        held_flow_reference(&duct, &vec![0.0; duct.state_dim()], 0, valve_flow, duration);
    assert!((pressure / observed - 1.0).abs() < 2e-9);
}

// Integrate the original complete body at a fixed physical input. Pressure and
// original R heat are independent accumulated observables, not energy remainders.
fn held_flow_reference(
    system: &fs_phs::PortHamiltonian,
    initial: &[f64],
    source: usize,
    flow: f64,
    duration: f64,
) -> (Vec<f64>, f64, f64) {
    let n = system.state_dim();
    let m = system.port_dim();
    let (j, r, g) = system.structure();
    let rhs = |state: &[f64]| {
        let effort = system.effort(&state[..n]);
        let pressure: f64 = (0..n).map(|row| g[row * m + source] * effort[row]).sum();
        let mut result = vec![0.0; n + 2];
        for row in 0..n {
            result[row] = (0..n)
                .map(|column| (j[row * n + column] - r[row * n + column]) * effort[column])
                .sum::<f64>()
                + g[row * m + source] * flow;
            result[n + 1] += effort[row]
                * (0..n)
                    .map(|column| r[row * n + column] * effort[column])
                    .sum::<f64>();
        }
        result[n] = pressure;
        result
    };
    let mut state = vec![0.0; n + 2];
    state[..n].copy_from_slice(initial);
    let h = duration / 2048.0;
    for _ in 0..2048 {
        let a = rhs(&state);
        let b = rhs(&state
            .iter()
            .zip(&a)
            .map(|(x, d)| x + 0.5 * h * d)
            .collect::<Vec<_>>());
        let c = rhs(&state
            .iter()
            .zip(&b)
            .map(|(x, d)| x + 0.5 * h * d)
            .collect::<Vec<_>>());
        let d = rhs(&state
            .iter()
            .zip(&c)
            .map(|(x, d)| x + h * d)
            .collect::<Vec<_>>());
        for k in 0..state.len() {
            state[k] += h * (a[k] + 2.0 * b[k] + 2.0 * c[k] + d[k]) / 6.0;
        }
    }
    let pressure = state[n] / duration;
    let heat = state[n + 1];
    state.truncate(n);
    (state, pressure, heat)
}

#[test]
fn quasistatic_warm_complete_body_closes_mouth_jet_and_original_loss_work() {
    let gas = gas();
    let duration = 1.0 / f64::from(RATE);
    let mut warm = bank([0.7, -0.4], [0.015, -0.02], true);
    for _ in 0..24 {
        warm.advance_controlled(0.7, duration, gas.density, 1.0, || false)
            .unwrap();
    }
    assert!(warm.vk.as_ref().unwrap().state().last().unwrap().abs() > 0.0);
    let line = || {
        ode_bore_chain(
            &sections(),
            ChainOutlet::PressureRelease,
            2,
            &[],
            &gas,
            None,
        )
        .unwrap()
    };
    let (reference, initial, observation) =
        compose_plate_chain(&mut warm, line(), 1, &gas, RATE).unwrap();
    assert!(observation.is_some());
    let source = 1;
    let before = reference.hamiltonian(&initial);
    assert!(before > 0.0);
    for mouth in [0.0, 100.0] {
        let reed = reed(0.0, mouth);
        let (system, state, _) = compose_plate_chain(&mut warm, line(), 1, &gas, RATE).unwrap();
        assert_eq!(
            state, initial,
            "the physical bank's complete warm history changed"
        );
        // Observe the fluid port directly. Omitting the separate microphone
        // observation leaves every physical surface and feedback path in place.
        let pressure =
            realize_quasistatic(system, state, None, source, reed, &gas, 1.0, duration, 1).unwrap()
                [0];
        let drop = mouth - pressure;
        let opening =
            quasistatic_aperture_opening(reed.rest_opening_m, reed.closing_pressure_pa, drop);
        let flow = bernoulli_volume_flow(reed.width_m, opening, drop, gas.density);
        let (next, mean_pressure, heat) =
            held_flow_reference(&reference, &initial, source, flow, duration);
        assert!((pressure - mean_pressure).abs() < 2e-9 * pressure.abs().max(mean_pressure.abs()));
        let after = reference.hamiltonian(&next);
        let jet_heat = duration * drop * flow;
        let mouth_work = duration * mouth * flow;
        let bore_work = duration * mean_pressure * flow;
        let scale = before.abs()
            + after.abs()
            + mouth_work.abs()
            + bore_work.abs()
            + heat.abs()
            + jet_heat.abs();
        assert!(heat >= 0.0 && jet_heat >= 0.0);
        assert!((after - before + heat - bore_work).abs() <= 2e-9 * scale);
        assert!((after - before + heat + jet_heat - mouth_work).abs() <= 2e-9 * scale);
        if mouth == 0.0 {
            assert!(
                pressure.abs() > 1e-5,
                "complete warm body must drive the valve"
            );
            assert!(
                after < before,
                "unpowered warm body must pay for both loss channels"
            );
        }
    }
}

#[test]
fn nonlinear_body_uses_the_exponential_owner_without_weakening_midpoint_admission() {
    let mut bank = PlateBank::default();
    bank.vk = Some(
        VkBody::from_modal_storage_with_damping(
            SosModalStorage {
                omegas: vec![600.0],
                channels: vec![fs_nlmodal::StressChannel {
                    coefficient: 2e5,
                    coupling: vec![1.0],
                }],
            },
            &[60.0],
            &[0.0],
            vec![0.015],
        )
        .unwrap(),
    );
    let gas = gas();
    let line = ode_bore_chain(
        &sections(),
        ChainOutlet::PressureRelease,
        2,
        &[],
        &gas,
        None,
    )
    .unwrap();
    let (original, initial, _) = compose_plate_chain(&mut bank, line, 1, &gas, RATE).unwrap();
    let split = prepare_render_linear(&original).unwrap();
    assert!(!split.is_exact_quadratic());
    assert!(matches!(
        fs_phs::linear_exact::step_render(
            &original,
            &split,
            &initial,
            &[0.0, 0.0],
            1.0 / f64::from(RATE)
        ),
        Err(fs_phs::PhsError::TimeResolution(_))
    ));
    for rate in [48_000, 96_000, 192_000] {
        for mass in [0.0, 1e-4] {
            let pressure = realize_reed_ode(
                &sections(),
                ChainOutlet::PressureRelease,
                &[],
                reed(mass, 100.0),
                &mut bank,
                &gas,
                1.0,
                rate,
                8,
                None,
            )
            .unwrap();
            assert!(pressure.iter().all(|p| p.is_finite()));
            assert!(pressure.iter().any(|p| p.abs() > 1e-4));
        }
    }
    assert!(
        bank.vk
            .as_ref()
            .unwrap()
            .state()
            .iter()
            .all(|value| *value == 0.0)
    );
}

// Independently integrate the continuous coupled equations. Only the physical
// bore matrices/storage are shared with the renderer: no prepared transition,
// mean-port map, reed junction, contact solve or production stepper enters this
// smooth, noncontact reference. The final coordinate integrates bore pressure.
fn continuous_massive_reference(
    reed: BeatingReed,
    rate: u32,
    windows: usize,
    subdivisions: usize,
) -> Vec<f64> {
    let gas = gas();
    let bore = ode_bore_chain(
        &sections(),
        ChainOutlet::PressureRelease,
        1,
        &[],
        &gas,
        None,
    )
    .unwrap();
    continuous_complete_reference(
        &bore,
        &vec![0.0; bore.state_dim()],
        None,
        0,
        reed,
        rate,
        windows,
        subdivisions,
    )
}

#[allow(clippy::too_many_arguments)]
fn continuous_complete_reference(
    bore: &fs_phs::PortHamiltonian,
    initial: &[f64],
    observation: Option<&PlateBankObservation>,
    source: usize,
    reed: BeatingReed,
    rate: u32,
    windows: usize,
    subdivisions: usize,
) -> Vec<f64> {
    let gas = gas();
    let n = bore.state_dim();
    let m = bore.port_dim();
    let (j, r, g) = bore.structure();
    let face = reed.stiffness_n_m * reed.rest_opening_m / reed.closing_pressure_pa;
    let damping = 2.0 * reed.damping_ratio * (reed.stiffness_n_m * reed.mass_kg).sqrt();
    let rhs = |state: &[f64]| {
        let effort = bore.effort(&state[..n]);
        let pressure: f64 = (0..n).map(|row| g[row * m + source] * effort[row]).sum();
        let drop = reed.blowing_pressure_pa - pressure;
        let opening = if reed.mass_kg > 0.0 {
            let opening = reed.rest_opening_m + state[n];
            assert!(
                opening > 0.25 * reed.rest_opening_m,
                "smooth reference entered contact"
            );
            opening
        } else {
            reed.rest_opening_m * (1.0 - drop / reed.closing_pressure_pa).clamp(0.0, 1.0)
        };
        let jet = if reed.mass_kg == 0.0 && drop.abs() < 1e-12 {
            0.0
        } else {
            reed.width_m * opening * drop.signum() * (2.0 * drop.abs() / gas.density).sqrt()
        };
        let flow = jet
            - if reed.mass_kg > 0.0 {
                face * state[n + 1]
            } else {
                0.0
            };
        let mut result = vec![0.0; n + 3];
        for row in 0..n {
            result[row] = (0..n)
                .map(|column| (j[row * n + column] - r[row * n + column]) * effort[column])
                .sum::<f64>()
                + g[row * m + source] * flow;
        }
        if reed.mass_kg > 0.0 {
            result[n] = state[n + 1];
            result[n + 1] = (-reed.stiffness_n_m * state[n] - damping * state[n + 1] - face * drop)
                / reed.mass_kg;
        }
        result[n + 2] = pressure;
        result
    };
    let duration = 1.0 / f64::from(rate);
    let h = duration / subdivisions as f64;
    let mut state = vec![0.0; n + 3];
    state[..n].copy_from_slice(initial);
    let mut pressure = Vec::with_capacity(windows);
    for _ in 0..windows {
        let before_volume = observation.map_or(0.0, |o| o.volume_velocity(&state[..n]).unwrap());
        state[n + 2] = 0.0;
        for _ in 0..subdivisions {
            let a = rhs(&state);
            let b = rhs(&state
                .iter()
                .zip(&a)
                .map(|(x, d)| x + 0.5 * h * d)
                .collect::<Vec<_>>());
            let c = rhs(&state
                .iter()
                .zip(&b)
                .map(|(x, d)| x + 0.5 * h * d)
                .collect::<Vec<_>>());
            let d = rhs(&state
                .iter()
                .zip(&c)
                .map(|(x, d)| x + h * d)
                .collect::<Vec<_>>());
            for k in 0..state.len() {
                state[k] += h * (a[k] + 2.0 * b[k] + 2.0 * c[k] + d[k]) / 6.0;
            }
        }
        assert!(state.iter().all(|x| x.is_finite()));
        let after_volume = observation.map_or(0.0, |o| o.volume_velocity(&state[..n]).unwrap());
        pressure.push(
            state[n + 2] / duration
                + gas.density * (after_volume - before_volume)
                    / (core::f64::consts::TAU * duration),
        );
    }
    pressure
}

fn resonant_reed(mass: f64) -> BeatingReed {
    BeatingReed {
        rest_opening_m: 4e-4,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa: 300.0,
        attack_s: 0.0,
        mass_kg: mass,
        stiffness_n_m: 500.0,
        damping_ratio: 0.04,
    }
}

fn massive_trace(reed: BeatingReed, rate: u32, windows: usize) -> Vec<f64> {
    realize_reed_ode(
        &sections(),
        ChainOutlet::PressureRelease,
        &[],
        reed,
        &mut PlateBank::default(),
        &gas(),
        1.0,
        rate,
        windows,
        None,
    )
    .unwrap()
}

#[test]
fn massive_complete_bore_refines_to_the_simultaneous_continuous_equations() {
    let reed = resonant_reed(1e-5);
    let rate = 24_000;
    let windows = 24;
    let reference = continuous_massive_reference(reed, rate, windows, 512);
    let mut errors = Vec::new();
    for refinement in [1, 2, 4] {
        let pressure = massive_trace(reed, rate * refinement as u32, windows * refinement);
        // Compare the same physical windows, not a finer endpoint against a
        // coarser interval mean. This also tests the renderer's output meaning.
        let averaged: Vec<f64> = pressure
            .chunks_exact(refinement)
            .map(|chunk| chunk.iter().sum::<f64>() / refinement as f64)
            .collect();
        let error = relative_difference(&averaged, &reference);
        assert!(error.is_finite());
        errors.push(error);
    }
    assert!(
        errors[2] < 2e-3,
        "complete coupled pressure errors: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[0] > 3.0 * pair[1],
            "same-window coupling must refine: {errors:?}"
        );
    }
}

#[test]
fn equal_static_compliance_with_distinct_reed_mass_changes_actual_bore_motion() {
    let mut traces = Vec::new();
    for mass in [1e-5, 4e-5] {
        let reed = resonant_reed(mass);
        let reference = continuous_massive_reference(reed, RATE, STEPS, 256);
        let actual = massive_trace(reed, RATE, STEPS);
        let error = relative_difference(&actual, &reference);
        assert!(error < 4e-3, "mass={mass}, coupled bore error={error}");
        traces.push(actual);
    }
    let difference = relative_difference(&traces[0], &traces[1]);
    assert!(
        difference > 1e-3,
        "mass affected only an unobserved reed: {difference}"
    );
}

fn nonlinear_bank(coefficient: f64) -> PlateBank {
    let body = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![600.0, 900.0],
            channels: vec![fs_nlmodal::StressChannel {
                coefficient,
                coupling: vec![1.0, 0.3, 0.3, 0.5],
            }],
        },
        &[60.0, -70.0, -70.0, 90.0],
        &[0.7, -0.4],
        vec![0.015, -0.02],
    )
    .unwrap()
    .with_relaxation_branches(
        vec![fs_phs::RelaxationBranch {
            projection: vec![0.7, 0.0, -0.4, 0.0],
            stiffness: 1600.0,
            relaxation_time_s: 0.0008,
        }],
        1,
    )
    .unwrap();
    let mut bank = PlateBank::default();
    bank.linear.push(compact(0.7, 0.003));
    bank.vk = Some(body);
    bank
}

fn nonlinear_image(
    coefficient: f64,
    rate: u32,
) -> (
    fs_phs::PortHamiltonian,
    Vec<f64>,
    Option<PlateBankObservation>,
) {
    let gas = gas();
    let line = ode_bore_chain(
        &sections(),
        ChainOutlet::PressureRelease,
        2,
        &[],
        &gas,
        None,
    )
    .unwrap();
    let (system, mut state, observation) =
        compose_plate_chain(&mut nonlinear_bank(coefficient), line, 1, &gas, rate).unwrap();
    assert_eq!(observation.as_ref().unwrap().state_dim(), 7);
    // Original order: compact q/p, two Vk q/p pairs, then Maxwell memory.
    // Both coefficient selections start with exactly the same physical state.
    state[..7].copy_from_slice(&[2e-6, 1e-5, 5e-6, 0.01, -3e-6, -0.006, 8e-5]);
    (system, state, observation)
}

#[test]
fn nonlinear_complete_reed_body_retains_finite_amplitude_stress_and_refines() {
    let coefficient = 1e17;
    let windows = 8;
    let (system, initial, observation) = nonlinear_image(coefficient, RATE);
    let split = prepare_render_linear(&system).unwrap();
    assert!(!split.is_exact_quadratic());
    let n = system.state_dim();
    let q = split.quadratic_matrix();
    let quadratic: f64 = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| 0.5 * initial[i] * q[i * n + j] * initial[j])
                .sum::<f64>()
        })
        .sum();
    assert!(
        system.hamiltonian(&initial) > 1.05 * quadratic,
        "fixture must carry material nonlinear storage, not an infinitesimal channel"
    );
    for mass in [0.0, 1e-5] {
        let reed = resonant_reed(mass);
        let reference = continuous_complete_reference(
            &system,
            &initial,
            observation.as_ref(),
            1,
            reed,
            RATE,
            windows,
            1024,
        );
        let mut errors = Vec::new();
        let mut finest = Vec::new();
        for refinement in [1, 2, 4] {
            let rate = RATE * refinement as u32;
            let (system, state, observer) = nonlinear_image(coefficient, rate);
            let pressure = nonlinear::realize_nonlinear(
                system,
                state,
                observer,
                1,
                reed,
                &gas(),
                1.0,
                1.0 / f64::from(rate),
                windows * refinement,
            )
            .unwrap();
            let average: Vec<f64> = pressure
                .chunks_exact(refinement)
                .map(|chunk| chunk.iter().sum::<f64>() / refinement as f64)
                .collect();
            errors.push(relative_difference(&average, &reference));
            finest = average;
        }
        assert!(
            errors[2] < 2e-3,
            "mass={mass}: finite-amplitude errors {errors:?}"
        );
        for pair in errors.windows(2) {
            assert!(
                pair[0] > 2.5 * pair[1],
                "mass={mass}: original nonlinear flow did not refine {errors:?}"
            );
        }
        let (linear, same_state, observer) = nonlinear_image(0.0, 4 * RATE);
        assert_eq!(same_state, initial);
        let linear_pressure = if mass == 0.0 {
            realize_quasistatic(
                linear,
                same_state,
                observer,
                1,
                reed,
                &gas(),
                1.0,
                1.0 / f64::from(4 * RATE),
                4 * windows,
            )
            .unwrap()
        } else {
            realize_massive(
                linear,
                same_state,
                observer,
                1,
                reed,
                &gas(),
                1.0,
                1.0 / f64::from(4 * RATE),
                4 * windows,
            )
            .unwrap()
        };
        let linear_average: Vec<f64> = linear_pressure
            .chunks_exact(4)
            .map(|chunk| chunk.iter().sum::<f64>() / 4.0)
            .collect();
        assert!(
            relative_difference(&finest, &linear_average) > 1e-3,
            "quartic material stress must change the actual coupled pressure"
        );
    }
}

#[test]
fn refused_nonlinear_observation_preserves_every_original_bank_history() {
    let gas = gas();
    let mut bank = nonlinear_bank(1e17);
    for _ in 0..24 {
        bank.advance_controlled(0.7, 1.0 / f64::from(RATE), gas.density, 1.0, || false)
            .unwrap();
    }
    let initial = bank
        .physical_image(&gas, RATE, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap()
        .state()
        .to_vec();
    assert!(initial.last().unwrap().abs() > 0.0);
    for mass in [0.0, 1e-5] {
        let result = realize_reed_ode(
            &sections(),
            ChainOutlet::PressureRelease,
            &[],
            resonant_reed(mass),
            &mut bank,
            &gas,
            0.0,
            RATE,
            8,
            None,
        );
        assert!(matches!(
            result,
            Err(AcousticRealizeError::InvalidDescription { .. })
        ));
        let retained = bank
            .physical_image(&gas, RATE, &fs_exec::CancelGate::default())
            .unwrap()
            .unwrap();
        assert_eq!(retained.state(), initial);
    }
}
