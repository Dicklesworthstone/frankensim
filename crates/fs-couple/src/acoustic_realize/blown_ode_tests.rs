//! G0/G1/G3: actual prescribed-flow rendering excites the signed physical
//! surface map, retains complete Vk history and reuses the original ODE owner.
use super::*;
use crate::thin_plate::CompactBody;
use fs_nlmodal::SosModalStorage;

const RATE: u32 = 48_000;
const STEPS: usize = 96;

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

fn render(bank: &mut PlateBank, peak: f64) -> Vec<f64> {
    realize_blown_ode(
        &sections(),
        ChainOutlet::PressureRelease,
        &[],
        VolumeVelocityPulse {
            peak_m3_s: peak,
            duration_s: 0.001,
        },
        bank,
        &gas(),
        1.0,
        RATE,
        STEPS,
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

fn compact(drive: f64, area: f64) -> PlateBank {
    let mut body = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 0.003,
        mass_kg: 0.02,
        frequency_hz: 400.0 / core::f64::consts::TAU,
        damping_ratio: 0.02,
    })
    .unwrap();
    body.area_m2 = area;
    body.drive_participation = drive;
    let mut bank = PlateBank::default();
    bank.linear.push(body);
    bank
}

fn vk(drive: [f64; 2], area: [f64; 2]) -> PlateBank {
    let body = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![600.0, 900.0],
            channels: vec![],
        },
        // Positive definite (determinant 500), but not diagonally dominant:
        // the material block requires its general physical loss factor while
        // the duct's singular Foster blocks retain their passive graph roots.
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
    bank
}

#[test]
fn prescribed_flow_uses_surface_area_even_when_the_point_drive_is_zero() {
    let zero_drive = render(&mut compact(0.0, 0.003), 1e-5);
    let drive = render(&mut compact(2.5, 0.003), 1e-5);
    let reverse_drive = render(&mut compact(-1.7, 0.003), 1e-5);
    let silent_surface = render(&mut compact(2.5, 0.0), 1e-5);
    assert!(zero_drive.iter().all(|p| p.is_finite()));
    assert!(relative_difference(&zero_drive, &drive) < 1e-10);
    assert!(relative_difference(&zero_drive, &reverse_drive) < 1e-10);
    assert!(
        relative_difference(&zero_drive, &silent_surface) > 1e-3,
        "pressure-driven body motion must remain present with zero point-force participation"
    );
}

#[test]
fn vk_only_and_mixed_banks_retain_signed_surface_coupling_and_material_loss() {
    let mut only_vk = vk([0.0, 0.0], [0.015, -0.02]);
    let response = render(&mut only_vk, 1e-5);
    let other_drive = render(&mut vk([1.4, -0.6], [0.015, -0.02]), 1e-5);
    let reversed_surface = render(&mut vk([0.0, 0.0], [0.015, 0.02]), 1e-5);
    let bare = render(&mut PlateBank::default(), 1e-5);
    assert!(relative_difference(&response, &other_drive) < 1e-10);
    assert!(
        relative_difference(&response, &bare) > 1e-3,
        "Vk-only must be an actual moving boundary, not the empty-bank route"
    );
    assert!(
        relative_difference(&response, &reversed_surface) > 1e-5,
        "relative signed areas interact with the original nonproportional material law"
    );
    let mut mixed = vk([0.0, 0.0], [0.015, -0.02]);
    mixed.linear = compact(0.0, 0.003).linear;
    let mixed_response = render(&mut mixed, 1e-5);
    assert!(
        relative_difference(&mixed_response, &response) > 1e-3,
        "an additional compact mode must survive alongside Vk and its memory"
    );
}

#[test]
fn a_warm_material_source_remains_audible_without_a_new_blow() {
    let gas = gas();
    let mut warm = vk([0.7, -0.4], [0.015, -0.02]);
    for _ in 0..24 {
        warm.advance_controlled(0.7, 1.0 / f64::from(RATE), gas.density, 1.0, || false)
            .unwrap();
    }
    let image = warm
        .physical_image(&gas, RATE, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap();
    assert_eq!(image.state_dim(), 5);
    assert!(
        image.state()[4].abs() > 1e-8,
        "the source has live Maxwell history"
    );
    let response = render(&mut warm, 0.0);
    let quiescent = render(&mut vk([0.7, -0.4], [0.015, -0.02]), 0.0);
    assert!(quiescent.iter().all(|p| *p == 0.0));
    assert!(response.iter().all(|p| p.is_finite()));
    assert!(
        response.iter().any(|p| p.abs() > 1e-5),
        "accepted vibration must not be zeroed by ODE routing"
    );
}

fn nonlinear_warm_bank(coefficient: f64, force: f64) -> PlateBank {
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
            relaxation_time_s: 2e-5,
        }],
        1,
    )
    .unwrap();
    let mut bank = compact(0.7, 0.003);
    bank.vk = Some(body);
    // The original bank owner admits this short preparation clock. Every
    // renderer below receives exactly this accepted history, independent of
    // its own later sample rate. No test writes a physical state directly.
    for _ in 0..256 {
        bank.advance_controlled(force, 5e-7, gas().density, 1.0, || false)
            .unwrap();
    }
    bank
}

#[test]
fn nonlinear_prescribed_flow_retains_warm_memory_and_amplitude_dependent_motion() {
    let mut strong = nonlinear_warm_bank(1e17, 200.0);
    let initial = strong
        .physical_image(&gas(), RATE, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap()
        .state()
        .to_vec();
    assert_eq!(initial.len(), 7);
    assert!(
        initial[6].abs() > 1e-8,
        "rapid material memory must already be populated"
    );
    let strain =
        initial[2] * initial[2] + 0.6 * initial[2] * initial[4] + 0.5 * initial[4] * initial[4];
    let nonlinear_energy = 0.25 * 1e17 * strain * strain;
    let linear_strain_energy =
        0.5 * (600.0_f64.powi(2) * initial[2].powi(2) + 900.0_f64.powi(2) * initial[4].powi(2));
    assert!(nonlinear_energy > 0.01 * linear_strain_energy);
    let ringing = render(&mut strong, 0.0);
    let driven = render(&mut strong, 1e-5);
    let weak = render(&mut nonlinear_warm_bank(1e17, 100.0), 0.0);
    let doubled_weak: Vec<_> = weak.iter().map(|p| 2.0 * p).collect();
    let quadratic = render(&mut nonlinear_warm_bank(0.0, 200.0), 0.0);
    assert!(ringing.iter().all(|p| p.is_finite()));
    assert!(ringing.iter().any(|p| p.abs() > 1e-5));
    assert!(
        relative_difference(&ringing, &driven) > 1e-4,
        "new physical flow must reach the retained nonlinear body"
    );
    assert!(
        relative_difference(&ringing, &quadratic) > 1e-3,
        "the original quartic force must remain active"
    );
    assert!(
        relative_difference(&ringing, &doubled_weak) > 1e-3,
        "nonlinear motion must depend on physical amplitude"
    );
    let preserved = strong
        .physical_image(&gas(), RATE, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap();
    assert_eq!(preserved.state(), initial);
    assert!(
        render(&mut nonlinear_warm_bank(1e17, 0.0), 0.0)
            .iter()
            .all(|p| *p == 0.0)
    );
}

#[test]
fn zero_flow_nonlinear_body_refines_against_independent_complete_equations() {
    let gas = gas();
    let mut bank = nonlinear_warm_bank(1e17, 200.0);
    let line = ode_bore_chain(
        &sections(),
        ChainOutlet::PressureRelease,
        2,
        &[],
        &gas,
        None,
    )
    .unwrap();
    let (system, initial, _) = compose_plate_chain(&mut bank, line, 1, &gas, RATE).unwrap();
    let n = system.state_dim();
    let ports = system.port_dim();
    let (j, r, g) = system.structure();
    let rhs = |state: &[f64]| {
        let effort = system.effort(&state[..n]);
        let mut result = vec![0.0; n + 1];
        for row in 0..n {
            result[row] = (0..n)
                .map(|column| (j[row * n + column] - r[row * n + column]) * effort[column])
                .sum();
        }
        result[n] = (0..n).map(|row| g[row * ports + 1] * effort[row]).sum();
        result
    };
    // Independent original pressure observation: the compact mode is mass
    // normalized, while the supplied Vk surface areas are already modal.
    let surface_flow = |x: &[f64]| 0.003 / 0.02_f64.sqrt() * x[1] + 0.015 * x[3] - 0.02 * x[5];
    let duration = 1.0 / f64::from(RATE);
    let windows = 24;
    let subdivisions = 256;
    let h = duration / f64::from(subdivisions);
    let mut x = vec![0.0; n + 1];
    x[..n].copy_from_slice(&initial);
    let mut reference = Vec::with_capacity(windows);
    for _ in 0..windows {
        let before_flow = surface_flow(&x);
        x[n] = 0.0;
        for _ in 0..subdivisions {
            let a = rhs(&x);
            let b = rhs(&x
                .iter()
                .zip(&a)
                .map(|(v, d)| v + 0.5 * h * d)
                .collect::<Vec<_>>());
            let c = rhs(&x
                .iter()
                .zip(&b)
                .map(|(v, d)| v + 0.5 * h * d)
                .collect::<Vec<_>>());
            let d = rhs(&x.iter().zip(&c).map(|(v, d)| v + h * d).collect::<Vec<_>>());
            for i in 0..x.len() {
                x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
        reference.push(
            x[n] / duration
                + gas.density * (surface_flow(&x) - before_flow)
                    / (core::f64::consts::TAU * duration),
        );
    }
    let mut errors = Vec::new();
    for refinement in [1, 2, 4] {
        let response = realize_blown_ode(
            &sections(),
            ChainOutlet::PressureRelease,
            &[],
            VolumeVelocityPulse {
                peak_m3_s: 0.0,
                duration_s: 0.001,
            },
            &mut bank,
            &gas,
            1.0,
            RATE * refinement as u32,
            windows * refinement,
            None,
        )
        .unwrap();
        let averaged: Vec<f64> = response
            .chunks_exact(refinement)
            .map(|window| window.iter().sum::<f64>() / refinement as f64)
            .collect();
        errors.push(relative_difference(&averaged, &reference));
    }
    assert!(
        errors[2] < 2e-3,
        "complete nonlinear pressure errors: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[0] > 2.5 * pair[1],
            "nonlinear physical clock refinement: {errors:?}"
        );
    }
}
