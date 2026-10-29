use super::*;
use fs_material::gas::GasSpec;
use fs_nlmodal::{SosModalStorage, StressChannel};

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn compact(mass: f64, omega: f64, area: f64, drive: f64) -> CompactBody {
    let mut body = CompactBody::from_radiator(fs_scenario::RadiatingPlate {
        area_m2: area.abs(),
        mass_kg: mass,
        frequency_hz: omega / core::f64::consts::TAU,
        damping_ratio: 0.04,
    })
    .unwrap();
    body.area_m2 = area;
    body.drive_participation = drive;
    body
}

const AREAS: [f64; 3] = [0.03, -0.02, 0.045];
const DRIVE: [f64; 3] = [-0.8, 0.6, -0.9];
const CAVITY_MASS: f64 = 500.0;
const CAVITY_STIFFNESS: f64 = 20_000.0;
const CAVITY_RESISTANCE: f64 = 25.0;

fn mixed(quartic: f64) -> PlateBank {
    let mut linear = compact(2.0, 18.0, AREAS[0], DRIVE[0]);
    linear.y = 0.01;
    linear.v = -0.025;
    let mut vk = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![25.0, 37.0],
            channels: vec![StressChannel {
                coefficient: quartic,
                coupling: vec![1.0, 0.2, 0.2, 0.5],
            }],
        },
        &[1.4, -0.35, -0.35, 0.8],
        &DRIVE[1..],
        AREAS[1..].to_vec(),
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
    // Flow-driven compliance/inertance, the same canonical cavity law, with
    // independently chosen coefficients so its force and flow signs are visible.
    let sys = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, CAVITY_RESISTANCE],
        vec![1.0, 0.0],
        Box::new(
            fs_phs::QuadraticStorage::new(vec![CAVITY_STIFFNESS, 0.0, 0.0, CAVITY_MASS.recip()], 2)
                .unwrap(),
        ),
    )
    .unwrap();
    let cavity = PlateCavity {
        linear: prepare_render_linear(&sys).unwrap(),
        sys,
        x: vec![0.0003, -0.02],
    };
    PlateBank {
        linear: vec![linear],
        vk: Some(vk),
        cavity: Some(cavity),
        ..PlateBank::default()
    }
}

fn physical_state(bank: &PlateBank) -> Vec<f64> {
    let b = &bank.linear[0];
    let mut x = vec![b.y, b.v];
    x.extend_from_slice(bank.vk.as_ref().unwrap().state());
    x.extend_from_slice(&bank.cavity.as_ref().unwrap().x);
    x
}

fn mixed_rhs(x: &[f64], force: f64, quartic: f64) -> Vec<f64> {
    let mut dx = vec![0.0; 11];
    let pressure = CAVITY_STIFFNESS * x[7];
    let flow = x[8] / CAVITY_MASS;
    let strain = 4.0 * (0.7 * x[2] - 0.4 * x[4]) - x[6];
    let form = x[2] * x[2] + 0.4 * x[2] * x[4] + 0.5 * x[4] * x[4];
    dx[0] = x[1];
    dx[1] = (DRIVE[0] * force + AREAS[0] * pressure) / 2.0 - 18.0_f64.powi(2) * x[0] - 1.44 * x[1];
    dx[2] = x[3];
    dx[3] = DRIVE[1] * force + AREAS[1] * pressure - 625.0 * x[2] - 1.4 * x[3] + 0.35 * x[5]
        - 2.8 * strain
        - quartic * form * (x[2] + 0.2 * x[4]);
    dx[4] = x[5];
    dx[5] = DRIVE[2] * force + AREAS[2] * pressure - 1369.0 * x[4] + 0.35 * x[3] - 0.8 * x[5]
        + 1.6 * strain
        - quartic * form * (0.2 * x[2] + 0.5 * x[4]);
    dx[6] = strain / 0.08;
    dx[7] = flow - AREAS[0] * x[1] - AREAS[1] * x[3] - AREAS[2] * x[5];
    dx[8] = -pressure - CAVITY_RESISTANCE * flow;
    dx[9] = force * (DRIVE[0] * x[1] + DRIVE[1] * x[3] + DRIVE[2] * x[5]);
    dx[10] = 2.88 * x[1] * x[1] + 1.4 * x[3] * x[3] - 0.7 * x[3] * x[5]
        + 0.8 * x[5] * x[5]
        + strain * strain / 0.08
        + CAVITY_RESISTANCE * flow * flow;
    dx
}

fn rk4(x: &mut [f64], dt: f64, rhs: impl Fn(&[f64]) -> Vec<f64>) {
    let k1 = rhs(x);
    let stage = |k: &[f64], s: f64| {
        x.iter()
            .zip(k)
            .map(|(a, b)| a + s * dt * b)
            .collect::<Vec<_>>()
    };
    let k2 = rhs(&stage(&k1, 0.5));
    let k3 = rhs(&stage(&k2, 0.5));
    let k4 = rhs(&stage(&k3, 1.0));
    for i in 0..x.len() {
        x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

fn mixed_error(actual: &[f64], expected: &[f64]) -> f64 {
    let weights = [
        2.0 * 324.0,
        2.0,
        625.0,
        1.0,
        1369.0,
        1.0,
        1.0,
        CAVITY_STIFFNESS,
        CAVITY_MASS.recip(),
    ];
    actual
        .iter()
        .zip(expected)
        .zip(weights)
        .map(|((a, b), w)| w * (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn complete_mixed_bank_matches_physical_equations_and_clock_refinement() {
    let (force, dt, steps) = (0.7, 0.0004, 16);
    let mut bank = mixed(0.0);
    let mut half = mixed(0.0);
    let initial_energy = bank.stored_energy_j();
    let mut reference = physical_state(&bank);
    reference.extend_from_slice(&[0.0, 0.0]);
    let (mut work, mut heat) = (0.0, 0.0);
    for step in 0..steps {
        let old_flow = bank.volume_velocity();
        let record = bank
            .advance_controlled(force, dt, 1.2, 1.7, || false)
            .unwrap();
        work += record.work.supplied;
        heat += record.work.dissipated;
        let expected_pressure =
            1.2 * (bank.volume_velocity() - old_flow) / (core::f64::consts::TAU * 1.7 * dt);
        assert!((record.pressure_pa - expected_pressure).abs() < 1e-13);
        assert!((record.volume_velocity_m3_s - bank.volume_velocity()).abs() < 1e-15);
        assert!(record.work.dissipated >= 0.0);
        // Start with one full clock, then refine a moving body with nonzero
        // Maxwell and cavity histories. No constructor resets are available.
        if step == 0 {
            half.advance_controlled(force, dt, 1.2, 1.7, || false)
                .unwrap();
        } else {
            for _ in 0..2 {
                half.advance_controlled(force, dt / 2.0, 1.2, 1.7, || false)
                    .unwrap();
            }
        }
        for _ in 0..64 {
            rk4(&mut reference, dt / 64.0, |x| mixed_rhs(x, force, 0.0));
        }
    }
    assert!(mixed_error(&physical_state(&bank), &reference) < 2e-11);
    assert!(mixed_error(&physical_state(&bank), &physical_state(&half)) < 2e-11);
    assert!((work - reference[9]).abs() < 2e-13);
    assert!((heat - reference[10]).abs() < 2e-13);
    assert!((bank.stored_energy_j() - initial_energy + heat - work).abs() < 2e-12);
}

#[test]
fn nonlinear_mixed_bank_retains_cross_damping_memory_and_second_order_motion() {
    let (quartic, force, duration) = (300_000.0, 0.7, 0.0064);
    let mut reference = physical_state(&mixed(quartic));
    reference.extend_from_slice(&[0.0, 0.0]);
    for _ in 0..2048 {
        rk4(&mut reference, duration / 2048.0, |x| {
            mixed_rhs(x, force, quartic)
        });
    }
    let mut errors = Vec::new();
    let mut heat_errors = Vec::new();
    let mut work_errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut bank = mixed(quartic);
        let start = bank.stored_energy_j();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..steps {
            let record = bank
                .advance_controlled(force, duration / f64::from(steps), 1.2, 1.0, || false)
                .unwrap();
            heat += record.work.dissipated;
            work += record.work.supplied;
        }
        assert!((bank.stored_energy_j() - start + heat - work).abs() < 3e-12);
        assert!(heat > 0.0);
        errors.push(mixed_error(&physical_state(&bank), &reference));
        work_errors.push((work - reference[9]).abs());
        heat_errors.push((heat - reference[10]).abs());
    }
    for values in [&errors, &work_errors, &heat_errors] {
        assert!(
            values[0] / values[1] > 3.8 && values[1] / values[2] > 3.8,
            "{values:?}"
        );
    }
    assert!(errors[2] < 2e-7, "{errors:?}");
    assert!(work_errors[2] < 1e-9);
    assert!(heat_errors[2] < 1.1e-5 * reference[10]);
}

#[test]
fn rejected_observation_late_cancellation_and_topology_change_publish_no_history() {
    let mut bank = mixed(0.0);
    bank.advance_controlled(0.7, 0.0004, 1.2, 1.0, || false)
        .unwrap();
    let accepted = physical_state(&bank);
    let complete = bank.runtime.as_ref().unwrap().body.state().to_vec();
    let energy = bank.stored_energy_j();
    assert!(
        bank.advance_controlled(0.7, 0.0004, 1.2, 0.0, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        bank.advance_controlled(0.7, 0.0004, 1.2, 1.0, || {
            polls += 1;
            polls == 3
        }),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(polls, 3); // The prepared exact owner has calculated the candidate.
    assert_eq!(physical_state(&bank), accepted);
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), complete);
    assert_eq!(bank.stored_energy_j(), energy);
    bank.linear[0].drive_participation *= 2.0;
    assert!(
        bank.advance_controlled(0.7, 0.0004, 1.2, 1.0, || false)
            .is_err()
    );
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), complete);
    bank.linear.clear();
    bank.vk = None;
    assert!(
        bank.advance_controlled(0.7, 0.0004, 1.2, 1.0, || false)
            .is_err()
    );
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), complete);
}

#[test]
fn a_shared_piston_face_reacts_on_every_mode_through_one_complete_history() {
    let gas = air();
    let area = core::f64::consts::PI * 0.01_f64.powi(2);
    let face = std::sync::Arc::new(());
    let mut first = compact(2e-5, 600.0, area, 1.0);
    let mut second = compact(3e-5, 900.0, -0.6 * area, 0.0);
    for body in [&mut first, &mut second] {
        body.piston_face = Some(face.clone());
        body.piston_area_m2 = area;
    }
    let mut bank = PlateBank {
        linear: vec![first, second],
        ..PlateBank::default()
    };
    bank.attach_radiation_loads(&gas, 48_000).unwrap();
    assert_eq!(bank.radiation.len(), 1);
    assert_eq!(bank.radiation[0].members, [0, 1]);
    let poles = bank.radiation[0].source.model().poles.clone();
    let n = 4 + 2 * poles.len();
    let mut reference = vec![0.0; n + 2];
    let (force, dt) = (2e-4, 1.0 / 48_000.0);
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..12 {
        let record = bank
            .advance_controlled(force, dt, gas.density, 1.0, || false)
            .unwrap();
        work += record.work.supplied;
        heat += record.work.dissipated;
        for _ in 0..100 {
            rk4(&mut reference, dt / 100.0, |x| {
                let mut dx = vec![0.0; n + 2];
                let flow = area * x[1] - 0.6 * area * x[3];
                let pressure: f64 = poles
                    .iter()
                    .enumerate()
                    .map(|(i, p)| p.coupling[0] * x[5 + 2 * i])
                    .sum();
                for (i, mass, omega, a, f) in [
                    (0, 2e-5, 600.0, area, force),
                    (1, 3e-5, 900.0, -0.6 * area, 0.0),
                ] {
                    let damping = 2.0 * 0.04 * omega * mass;
                    dx[2 * i] = x[2 * i + 1];
                    dx[2 * i + 1] = (f - a * pressure - damping * x[2 * i + 1]) / mass
                        - omega * omega * x[2 * i];
                    dx[n + 1] += damping * x[2 * i + 1] * x[2 * i + 1];
                }
                for (i, p) in poles.iter().enumerate() {
                    let q = 4 + 2 * i;
                    let damping = 2.0 * p.zeta * p.omega;
                    dx[q] = x[q + 1];
                    dx[q + 1] =
                        p.coupling[0] * flow - p.omega * p.omega * x[q] - damping * x[q + 1];
                    dx[n + 1] += damping * x[q + 1] * x[q + 1];
                }
                dx[n] = force * x[1];
                dx
            });
        }
    }
    let actual = bank.runtime.as_ref().unwrap().body.state();
    assert_eq!(actual.len(), n);
    for (i, body) in bank.linear.iter().enumerate() {
        assert!((body.y - reference[2 * i]).abs() < 1e-13);
        assert!((body.v - reference[2 * i + 1]).abs() < 2e-10);
    }
    assert!(
        bank.linear[1].v.abs() > 1e-6,
        "an undriven mode must feel the shared pressure"
    );
    let mut memory_error = 0.0;
    for (i, p) in poles.iter().enumerate() {
        let q = 4 + 2 * i;
        memory_error += (p.omega * (actual[q] - reference[q])).powi(2)
            + (actual[q + 1] - reference[q + 1]).powi(2);
    }
    assert!(memory_error.sqrt() < 2e-8 * work.sqrt());
    assert!((heat - reference[n + 1]).abs() < 2e-8 * work);
    assert!((work - reference[n]).abs() < 2e-8 * work);
    assert!((bank.stored_energy_j() + heat - work).abs() < 2e-11 * work);
    let accepted = actual.to_vec();
    bank.linear[0].omega *= 1.1;
    assert!(
        bank.advance_controlled(force, dt, gas.density, 1.0, || false)
            .is_err()
    );
    assert_eq!(bank.runtime.as_ref().unwrap().body.state(), accepted);
}

#[test]
fn an_uncoupled_large_bank_keeps_the_existing_exact_scalar_path() {
    let mut bank = PlateBank::default();
    for i in 0..128 {
        bank.linear.push(compact(
            0.02,
            300.0 + f64::from(i),
            0.001,
            if i % 2 == 0 { 0.8 } else { -0.3 },
        ));
    }
    let mut reference = bank.linear.clone();
    let expected: f64 = reference
        .iter_mut()
        .map(|b| {
            b.drive_and_radiate(0.7 * b.drive_participation, 1.0 / 48_000.0, 1.2, 1.0)
                .unwrap()
        })
        .sum();
    assert_eq!(
        bank.drive_and_radiate(0.7, 1.0 / 48_000.0, 1.2, 1.0)
            .unwrap(),
        expected
    );
    assert!(bank.runtime.is_none());
    for (a, b) in bank.linear.iter().zip(reference) {
        assert_eq!((a.y, a.v), (b.y, b.v));
    }
}
