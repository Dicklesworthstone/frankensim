//! G0/G1/G3: original material history and reciprocal surface coupling survive
//! the actual ODE assembly adapter; references use physical equations.
use super::*;
use crate::thin_plate::CompactBody;
use fs_nlmodal::SosModalStorage;
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

const AREA: [f64; 3] = [0.03, -0.02, 0.045];
const DRIVE: [f64; 3] = [-0.8, 0.6, -0.9];

fn gas() -> GasState {
    GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        293.15,
        101_325.0,
    )
    .unwrap()
}

fn duct() -> PortHamiltonian {
    // Flow-driven compliance and inertance; both openings share the physical
    // pressure effort. Port zero is the source, port one faces the structure.
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 25.0],
        vec![1.0, 1.0, 0.0, 0.0],
        Box::new(QuadraticStorage::new(vec![20_000.0, 0.0, 0.0, 1.0 / 500.0], 2).unwrap()),
    )
    .unwrap()
}

fn mixed_bank() -> PlateBank {
    let mut compact = CompactBody::from_radiator(RadiatingPlate {
        area_m2: AREA[0],
        mass_kg: 2.0,
        frequency_hz: 18.0 / core::f64::consts::TAU,
        damping_ratio: 0.04,
    })
    .unwrap();
    compact.drive_participation = DRIVE[0];
    let vk = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![25.0, 37.0],
            channels: vec![],
        },
        &[1.4, -0.35, -0.35, 0.8],
        &DRIVE[1..],
        AREA[1..].to_vec(),
    )
    .unwrap()
    .with_relaxation_branches(
        vec![RelaxationBranch {
            projection: vec![0.7, 0.0, -0.4, 0.0],
            stiffness: 16.0,
            relaxation_time_s: 0.08,
        }],
        1,
    )
    .unwrap();
    let mut bank = PlateBank::default();
    bank.linear.push(compact);
    bank.vk = Some(vk);
    bank
}

fn rk4(x: &mut [f64], dt: f64, rhs: impl Fn(&[f64]) -> Vec<f64>) {
    let k1 = rhs(x);
    let stage = |k: &[f64], scale: f64| {
        x.iter()
            .zip(k)
            .map(|(a, b)| a + scale * dt * b)
            .collect::<Vec<_>>()
    };
    let k2 = rhs(&stage(&k1, 0.5));
    let k3 = rhs(&stage(&k2, 0.5));
    let k4 = rhs(&stage(&k3, 1.0));
    for i in 0..x.len() {
        x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

fn physical_rhs(x: &[f64], force: f64, flow_source: f64) -> Vec<f64> {
    // Original compact [q,p]/sqrt(m), two mass-normalized Vk pairs,
    // Maxwell sqrt-energy coordinate, then actual fluid [compression,impulse].
    let mass_root = 2.0_f64.sqrt();
    let pressure = 20_000.0 * x[7];
    let flow = x[8] / 500.0;
    let strain = 4.0 * (0.7 * x[2] - 0.4 * x[4]) - x[6];
    vec![
        x[1],
        DRIVE[0] / mass_root * force + AREA[0] / mass_root * pressure - 324.0 * x[0] - 1.44 * x[1],
        x[3],
        DRIVE[1] * force + AREA[1] * pressure - 625.0 * x[2] - 1.4 * x[3] + 0.35 * x[5]
            - 2.8 * strain,
        x[5],
        DRIVE[2] * force + AREA[2] * pressure - 1369.0 * x[4] + 0.35 * x[3] - 0.8 * x[5]
            + 1.6 * strain,
        strain / 0.08,
        flow_source + flow - AREA[0] / mass_root * x[1] - AREA[1] * x[3] - AREA[2] * x[5],
        -pressure - 25.0 * flow,
        force * (DRIVE[0] / mass_root * x[1] + DRIVE[1] * x[3] + DRIVE[2] * x[5])
            + flow_source * pressure,
        1.44 * x[1] * x[1] + 1.4 * x[3] * x[3] - 0.7 * x[3] * x[5]
            + 0.8 * x[5] * x[5]
            + strain * strain / 0.08
            + 25.0 * flow * flow,
    ]
}

#[test]
fn complete_mixed_plate_and_duct_retain_warm_memory_and_independent_ports() {
    let gas = gas();
    let mut bank = mixed_bank();
    for _ in 0..12 {
        bank.advance_controlled(0.7, 0.0004, gas.density, 1.0, || false)
            .unwrap();
    }
    let source = bank
        .physical_image(&gas, 48_000, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap();
    let accepted = source.state().to_vec();
    assert_eq!(accepted.len(), 7);
    assert!(
        accepted[6].abs() > 1e-8,
        "pre-existing Maxwell history is live"
    );
    let energy = bank.stored_energy_j();
    let (system, mut x, observation) =
        compose_plate_chain(&mut bank, duct(), 1, &gas, 48_000).unwrap();
    let observation = observation.unwrap();
    assert_eq!(observation.state_dim(), 7);
    assert_eq!(x.len(), 9);
    assert_eq!(&x[..7], accepted);
    assert_eq!(system.port_dim(), 2);
    assert!((system.hamiltonian(&x) - energy).abs() < 1e-14);
    let (force, flow_source, dt) = (0.7, 1e-4, 0.0004);
    let linear = prepare_render_linear(&system).unwrap();
    let exact = prepare_chain_exact(&system, &linear, dt).unwrap();
    assert!(exact.is_some());
    assert!(
        step_chain_flow(
            &system,
            &linear,
            exact.as_ref(),
            &x,
            &[force, flow_source],
            dt * 2.0
        )
        .is_err()
    );
    let mut reference = x.clone();
    reference.extend_from_slice(&[0.0, 0.0]);
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..16 {
        let old_flow = observation.volume_velocity(&x[..7]).unwrap();
        let record = step_chain_flow(
            &system,
            &linear,
            exact.as_ref(),
            &x,
            &[force, flow_source],
            dt,
        )
        .unwrap();
        let pressure =
            plate_chain_radiation(Some(&observation), &record.x, &x, &gas, 1.7, dt).unwrap();
        let new_flow =
            AREA[0] / 2.0_f64.sqrt() * record.x[1] + AREA[1] * record.x[3] + AREA[2] * record.x[5];
        let expected = gas.density * (new_flow - old_flow) / (core::f64::consts::TAU * 1.7 * dt);
        assert!((pressure - expected).abs() < 1e-13);
        work += record.supplied;
        heat += record.dissipated;
        x = record.x;
        for _ in 0..64 {
            rk4(&mut reference, dt / 64.0, |state| {
                physical_rhs(state, force, flow_source)
            });
        }
    }
    let weights: [f64; 9] = [
        324.0,
        1.0,
        625.0,
        1.0,
        1369.0,
        1.0,
        1.0,
        20_000.0,
        1.0 / 500.0,
    ];
    let error = x
        .iter()
        .zip(&reference)
        .zip(weights)
        .map(|((a, b), w)| w * (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(error < 2e-11, "{error}");
    assert!((work - reference[9]).abs() < 2e-13);
    assert!((heat - reference[10]).abs() < 2e-13);
    assert!((system.hamiltonian(&x) - energy + heat - work).abs() < 2e-12);
    assert!(plate_chain_radiation(Some(&observation), &x[..6], &x, &gas, 1.0, dt).is_err());
}

#[test]
fn attached_radiation_remains_in_the_ode_prefix_and_same_source_history() {
    let gas = gas();
    let mut bank = PlateBank::default();
    bank.linear.push(
        CompactBody::from_radiator(RadiatingPlate {
            area_m2: core::f64::consts::PI * 0.01_f64.powi(2),
            mass_kg: 2e-5,
            frequency_hz: 600.0 / core::f64::consts::TAU,
            damping_ratio: 0.01,
        })
        .unwrap(),
    );
    bank.declare_radiation_loads(&gas, 48_000).unwrap();
    for _ in 0..12 {
        bank.advance_controlled(2e-4, 1.0 / 48_000.0, gas.density, 1.0, || false)
            .unwrap();
    }
    let image = bank
        .physical_image(&gas, 48_000, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap();
    let accepted = image.state().to_vec();
    assert!(accepted.len() > 2 && accepted[2..].iter().any(|x| *x != 0.0));
    let energy = bank.stored_energy_j();
    let (system, x, observation) = compose_plate_chain(&mut bank, duct(), 1, &gas, 48_000).unwrap();
    let observation = observation.unwrap();
    assert_eq!(observation.state_dim(), accepted.len());
    assert_eq!(&x[..accepted.len()], accepted);
    assert_eq!(x.len(), accepted.len() + 2);
    assert!((system.hamiltonian(&x) - energy).abs() < 1e-13 * energy);
    let linear = prepare_render_linear(&system).unwrap();
    let exact = prepare_chain_exact(&system, &linear, 1.0 / 48_000.0).unwrap();
    let record = step_chain_flow(
        &system,
        &linear,
        exact.as_ref(),
        &x,
        &[0.0, 0.0],
        1.0 / 48_000.0,
    )
    .unwrap();
    let pressure =
        plate_chain_radiation(Some(&observation), &record.x, &x, &gas, 1.0, 1.0 / 48_000.0)
            .unwrap();
    assert!(pressure.is_finite());
    assert!(record.dissipated > 0.0);
    assert!((record.delta_h + record.dissipated).abs() < 1e-11 * energy);
}

#[test]
fn empty_bank_preserves_the_original_duct_and_has_no_structural_observer() {
    let gas = gas();
    let original = duct();
    let mut bank = PlateBank::default();
    let (system, x, observation) = compose_plate_chain(&mut bank, duct(), 1, &gas, 48_000).unwrap();
    assert_eq!(system.structure(), original.structure());
    assert_eq!(x, [0.0, 0.0]);
    assert!(observation.is_none());
    assert_eq!(
        plate_chain_radiation(None, &x, &x, &gas, 1.0, 0.001).unwrap(),
        0.0
    );
}
