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
