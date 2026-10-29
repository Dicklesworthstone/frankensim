//! G0/G1/G3/G4: a physical image retains every original law, state and port.
use super::*;
use fs_material::gas::GasSpec;
use fs_nlmodal::{SosModalStorage, StressChannel};

const AREAS: [f64; 4] = [0.03, -0.015, -0.02, 0.045];
const DRIVES: [f64; 4] = [-0.8, 0.35, 0.6, -0.9];

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn flow_storage(mass: f64, stiffness: f64, resistance: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, resistance],
        vec![1.0, 0.0],
        Box::new(
            fs_phs::QuadraticStorage::new(vec![stiffness, 0.0, 0.0, mass.recip()], 2).unwrap(),
        ),
    )
    .unwrap()
}

fn compact(mass: f64, omega: f64, zeta: f64, area: f64, drive: f64) -> CompactBody {
    let mut body = CompactBody::from_radiator(fs_scenario::RadiatingPlate {
        area_m2: area.abs(),
        mass_kg: mass,
        frequency_hz: omega / core::f64::consts::TAU,
        damping_ratio: zeta,
    })
    .unwrap();
    body.area_m2 = area;
    body.drive_participation = drive;
    body
}

fn mixed() -> PlateBank {
    let mut first = compact(2.0, 18.0, 0.04, AREAS[0], DRIVES[0]);
    let mut second = compact(3.0, 30.0, 0.07, AREAS[1], DRIVES[1]);
    first.y = 0.01;
    first.v = -0.025;
    second.y = -0.004;
    second.v = 0.015;
    let mut vk = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![25.0, 37.0],
            channels: Vec::<StressChannel>::new(),
        },
        &[1.4, -0.35, -0.35, 0.8],
        &DRIVES[2..],
        AREAS[2..].to_vec(),
    )
    .unwrap()
    .with_relaxation_branches(
        vec![fs_phs::RelaxationBranch {
            projection: vec![0.7, 0.0, -0.4, 0.0],
            stiffness: 16.0,
            relaxation_time_s: 0.08,
        }],
        1,
    )
    .unwrap();
    vk.x.copy_from_slice(&[0.008, -0.03, -0.006, 0.04, 0.017]);
    let sys = flow_storage(500.0, 20_000.0, 25.0);
    let cavity = PlateCavity {
        linear: prepare_render_linear(&sys).unwrap(),
        sys,
        x: vec![0.0003, -0.02],
    };
    PlateBank {
        linear: vec![first, second],
        vk: Some(vk),
        cavity: Some(cavity),
        ..PlateBank::default()
    }
}

fn physical_state(bank: &PlateBank) -> Vec<f64> {
    let mut x = Vec::new();
    for b in &bank.linear {
        x.extend_from_slice(&[b.y, b.v]);
    }
    x.extend_from_slice(bank.vk.as_ref().unwrap().state());
    x.extend_from_slice(&bank.cavity.as_ref().unwrap().x);
    x
}

fn close(a: f64, b: f64, relative: f64, floor: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(floor),
        "{a} != {b}"
    );
}

#[test]
fn warmed_image_keeps_complete_state_independent_ports_and_observation() {
    let gas = air();
    let gate = CancelGate::default();
    let dt = 1.0 / 48_000.0;
    let mut bank = mixed();
    for _ in 0..7 {
        bank.advance_controlled(0.6, dt, gas.density, 1.7, || false)
            .unwrap();
    }
    let accepted = bank.runtime.as_ref().unwrap().body.state().to_vec();
    let physical = physical_state(&bank);
    let energy = bank.stored_energy_j();
    let image = bank.physical_image(&gas, 48_000, &gate).unwrap().unwrap();
    assert_eq!(image.state_dim(), 11); // Four mechanical pairs, one memory and one cavity pair.
    assert_eq!(image.state(), accepted);
    assert_eq!(image.system().port_dim(), 2);
    assert_eq!(image.system().hamiltonian(image.state()), energy);
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), accepted);
    let outputs = image.system().output(image.state());
    let velocity = [physical[1], physical[3], physical[5], physical[7]];
    let force_flow: f64 = DRIVES.iter().zip(velocity).map(|(a, v)| a * v).sum();
    let volume: f64 = AREAS.iter().zip(velocity).map(|(a, v)| a * v).sum();
    close(outputs[0], force_flow, 2e-14, 1e-12);
    close(outputs[1], volume, 2e-14, 1e-12);

    let (sys, initial, observer) = image.into_parts();
    assert_eq!(observer.state_dim(), initial.len());
    close(
        observer.volume_velocity(&initial).unwrap(),
        volume,
        2e-14,
        1e-12,
    );
    assert!(
        observer
            .volume_velocity(&initial[..initial.len() - 1])
            .is_err()
    );
    assert!(
        observer
            .pressure_between(&initial, &initial, 0.0, gas.density, 1.7)
            .is_err()
    );
    assert!(
        observer
            .pressure_between(&initial, &initial, dt, gas.density, 0.0)
            .is_err()
    );
    let split = prepare_render_linear(&sys).unwrap();
    let step = fs_phs::linear_exact::step_render(&sys, &split, &initial, &[0.4, 0.0], dt).unwrap();
    let live = bank
        .advance_controlled(0.4, dt, gas.density, 1.7, || false)
        .unwrap();
    for (a, b) in step
        .x
        .iter()
        .zip(bank.runtime.as_ref().unwrap().body.state())
    {
        close(*a, *b, 2e-12, 1e-10);
    }
    close(step.dissipated, live.work.dissipated, 2e-11, 1e-12);
    close(step.supplied, live.work.supplied, 2e-11, 1e-12);
    close(
        observer
            .pressure_between(&initial, &step.x, dt, gas.density, 1.7)
            .unwrap(),
        live.pressure_pa,
        2e-10,
        1e-10,
    );

    let current = bank.runtime.as_ref().unwrap().body.state().to_vec();
    let cancelled = CancelGate::default();
    cancelled.request();
    assert!(matches!(
        bank.physical_image(&gas, 48_000, &cancelled),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), current);
    bank.linear[0].drive_participation *= 2.0;
    assert!(bank.physical_image(&gas, 48_000, &gate).is_err());
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), current);
}

fn rhs(x: &[f64], force: f64) -> Vec<f64> {
    let mut dx = vec![0.0; 15];
    let cavity_pressure = 20_000.0 * x[9];
    let exterior_pressure = 40_000.0 * x[11];
    let net_pressure = cavity_pressure - exterior_pressure;
    let volume = AREAS[0] * x[1] + AREAS[1] * x[3] + AREAS[2] * x[5] + AREAS[3] * x[7];
    let strain = 4.0 * (0.7 * x[4] - 0.4 * x[6]) - x[8];
    dx[0] = x[1];
    dx[1] = (DRIVES[0] * force + AREAS[0] * net_pressure) / 2.0 - 324.0 * x[0] - 1.44 * x[1];
    dx[2] = x[3];
    dx[3] = (DRIVES[1] * force + AREAS[1] * net_pressure) / 3.0 - 900.0 * x[2] - 4.2 * x[3];
    dx[4] = x[5];
    dx[5] = DRIVES[2] * force + AREAS[2] * net_pressure - 625.0 * x[4] - 1.4 * x[5] + 0.35 * x[7]
        - 2.8 * strain;
    dx[6] = x[7];
    dx[7] = DRIVES[3] * force + AREAS[3] * net_pressure - 1369.0 * x[6] + 0.35 * x[5] - 0.8 * x[7]
        + 1.6 * strain;
    dx[8] = strain / 0.08;
    dx[9] = x[10] / 500.0 - volume;
    dx[10] = -cavity_pressure - 25.0 * x[10] / 500.0;
    dx[11] = x[12] / 200.0 + volume;
    dx[12] = -exterior_pressure - 10.0 * x[12] / 200.0;
    dx[13] = force * (DRIVES[0] * x[1] + DRIVES[1] * x[3] + DRIVES[2] * x[5] + DRIVES[3] * x[7]);
    dx[14] = 2.88 * x[1] * x[1] + 12.6 * x[3] * x[3] + 1.4 * x[5] * x[5] - 0.7 * x[5] * x[7]
        + 0.8 * x[7] * x[7]
        + strain * strain / 0.08
        + 25.0 * (x[10] / 500.0).powi(2)
        + 10.0 * (x[12] / 200.0).powi(2);
    dx
}

fn rk4(x: &mut [f64], dt: f64, force: f64) {
    let k1 = rhs(x, force);
    let stage = |k: &[f64], a: f64| {
        x.iter()
            .zip(k)
            .map(|(x, k)| x + a * dt * k)
            .collect::<Vec<_>>()
    };
    let k2 = rhs(&stage(&k1, 0.5), force);
    let k3 = rhs(&stage(&k2, 0.5), force);
    let k4 = rhs(&stage(&k3, 1.0), force);
    for i in 0..x.len() {
        x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

#[test]
fn a_complete_image_couples_to_an_external_flow_load_through_the_actual_surface() {
    let gas = air();
    let gate = CancelGate::default();
    let mut bank = mixed();
    let source_state = physical_state(&bank);
    let mut reference = source_state.clone();
    reference.extend_from_slice(&[0.0001, 0.02, 0.0, 0.0]);
    let (body, mut state, observer) = bank
        .physical_image(&gas, 48_000, &gate)
        .unwrap()
        .unwrap()
        .into_parts();
    // This external pressure opposes OUTWARD motion; the original internal
    // cavity uses the opposite orientation and still retains its own history.
    let sys = fs_phs::transformer(body, flow_storage(200.0, 40_000.0, 10.0), 1, 0, -1.0).unwrap();
    assert_eq!(sys.port_dim(), 1);
    assert_eq!(sys.state_dim(), 13);
    state.extend_from_slice(&[0.0001, 0.02]);
    let start = sys.hamiltonian(&state);
    let split = prepare_render_linear(&sys).unwrap();
    let (force, dt) = (0.7, 0.0004);
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..16 {
        let step = fs_phs::linear_exact::step_render(&sys, &split, &state, &[force], dt).unwrap();
        work += step.supplied;
        heat += step.dissipated;
        state = step.x;
        for _ in 0..64 {
            rk4(&mut reference, dt / 64.0, force);
        }
    }
    let mut physical = state.clone();
    for (i, mass) in [2.0_f64, 3.0].into_iter().enumerate() {
        physical[2 * i] /= mass.sqrt();
        physical[2 * i + 1] /= mass.sqrt();
    }
    let weights = [
        648.0,
        2.0,
        2700.0,
        3.0,
        625.0,
        1.0,
        1369.0,
        1.0,
        1.0,
        20_000.0,
        1.0 / 500.0,
        40_000.0,
        1.0 / 200.0,
    ];
    let error = physical
        .iter()
        .zip(&reference)
        .zip(weights)
        .map(|((a, b), w)| w * (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(error < 3e-11, "complete physical state error {error}");
    assert!((work - reference[13]).abs() < 3e-13);
    assert!((heat - reference[14]).abs() < 3e-13);
    assert!((sys.hamiltonian(&state) - start + heat - work).abs() < 3e-12);
    let volume: f64 = AREAS
        .iter()
        .enumerate()
        .map(|(i, a)| a * physical[2 * i + 1])
        .sum();
    close(
        observer
            .volume_velocity(&state[..observer.state_dim()])
            .unwrap(),
        volume,
        2e-14,
        1e-12,
    );
    assert_eq!(physical_state(&bank), source_state);
}

#[test]
fn shared_face_radiation_is_prepared_once_and_its_warm_memory_survives_export() {
    let gas = air();
    let gate = CancelGate::default();
    let area = core::f64::consts::PI * 0.01_f64.powi(2);
    let face = std::sync::Arc::new(());
    let mut first = compact(2e-5, 600.0, 0.04, area, 1.0);
    let mut second = compact(3e-5, 900.0, 0.04, -0.6 * area, 0.0);
    for body in [&mut first, &mut second] {
        body.piston_face = Some(face.clone());
        body.piston_area_m2 = area;
    }
    let mut bank = PlateBank {
        linear: vec![first, second],
        ..PlateBank::default()
    };
    bank.declare_radiation_loads(&gas, 48_000).unwrap();
    let cold = bank.physical_image(&gas, 48_000, &gate).unwrap().unwrap();
    assert!(
        bank.runtime.is_none(),
        "export does not prepare an unused time integrator"
    );
    assert_eq!(bank.radiation.len(), 1);
    assert_eq!(
        cold.state_dim(),
        4 + 2 * bank.radiation[0].source.model().poles.len()
    );
    assert!(cold.state().iter().all(|v| *v == 0.0));
    let dt = 1.0 / 48_000.0;
    for _ in 0..8 {
        bank.advance_controlled(2e-4, dt, gas.density, 1.0, || false)
            .unwrap();
    }
    let accepted = bank.runtime.as_ref().unwrap().body.state().to_vec();
    assert!(accepted[4..].iter().any(|v| *v != 0.0));
    let image = bank.physical_image(&gas, 48_000, &gate).unwrap().unwrap();
    assert_eq!(image.state(), accepted);
    assert_eq!(
        image.system().hamiltonian(image.state()),
        bank.stored_energy_j()
    );
    let (sys, initial, observer) = image.into_parts();
    let split = prepare_render_linear(&sys).unwrap();
    let step = fs_phs::linear_exact::step_render(&sys, &split, &initial, &[0.0, 0.0], dt).unwrap();
    let live = bank
        .advance_controlled(0.0, dt, gas.density, 1.0, || false)
        .unwrap();
    for (a, b) in step
        .x
        .iter()
        .zip(bank.runtime.as_ref().unwrap().body.state())
    {
        close(*a, *b, 1e-9, 1e-16);
    }
    close(step.dissipated, live.work.dissipated, 1e-9, 1e-20);
    close(
        observer.volume_velocity(&step.x).unwrap(),
        live.volume_velocity_m3_s,
        1e-10,
        1e-16,
    );
    let other = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202_650.0).unwrap();
    let current = bank.runtime.as_ref().unwrap().body.state().to_vec();
    assert!(bank.physical_image(&other, 48_000, &gate).is_err());
    assert!(bank.physical_image(&gas, 0, &gate).is_err());
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), current);
}
