use super::*;
use fs_duct::dg::DgMesh;
use fs_duct::{Duct, Segment};

const RHO: f64 = 1000.0;
const SOUND: f64 = 1500.0;
const REFERENCE: f64 = 100_000.0;
const VAPOR: f64 = 2000.0;
const RADIUS: f64 = 0.01;
fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}
fn material() -> CavitationState {
    CavitationState::try_new(293.15, RHO, 0.0, 0.0, VAPOR, SOUND).unwrap()
}
fn limits() -> DuctCavityLimits {
    DuctCavityLimits {
        reference_pressure_pa: REFERENCE,
        maximum_cavity_volume_m3: 3e-6,
        maximum_liquid_strain: 0.001,
        maximum_interface_mach: 0.01,
        maximum_vapor_density_ratio: 0.01,
    }
}
fn numerics() -> DuctCavityNumerics {
    DuctCavityNumerics {
        relative_tolerance: 1e-10,
        absolute_energy_tolerance_j: 1e-11,
        absolute_pressure_tolerance_pa: 1e-3,
        absolute_volume_tolerance_m3: 1e-14,
        maximum_step_s: 1e-5,
        event_time_tolerance_s: 1e-13,
        event_scan_substeps: 8,
        maximum_attempts: 30_000,
        maximum_events: 32,
    }
}
fn native(outlet: DgBoundary) -> DgDuct<DgFluidMaterial<CavitationState>> {
    DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: 1.0,
            }],
        },
        &material(),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(8)),
        outlet,
    )
    .unwrap()
}
fn spec(port: u32, volume: f64) -> DuctPocket {
    DuctPocket {
        port: PortId(port),
        initial_volume_m3: volume,
        maximum_volume_m3: 2e-6,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e} ({tolerance:e})");
}
fn gas(port: u32, volume: f64, pressure: f64, minimum: f64) -> TrappedGasPocket {
    TrappedGasPocket {
        port: PortId(port),
        gas: BubbleGas::from_reference_volume_m3(volume, pressure, 1.4).unwrap(),
        minimum_volume_m3: minimum,
    }
}

#[test]
fn gas_cushion_oscillates_at_the_liquid_column_frequency_and_conserves_work() {
    let volume = 1e-6;
    let inventory = gas(8, volume, REFERENCE - VAPOR, 0.2 * volume);
    let stiffness = inventory
        .gas
        .volume_response(volume)
        .unwrap()
        .stiffness_pa_per_m3;
    // Distributed lossless column: z*tan(z)=K_g*A*L/(rho*c^2), z=omega*L/c.
    // This analytic transmission-line oracle is independent of the DG operator.
    let target = stiffness * area() / (RHO * SOUND * SOUND);
    let (mut lo, mut hi) = (0.0f64, core::f64::consts::FRAC_PI_2);
    for _ in 0..70 {
        let mid = 0.5 * (lo + hi);
        if mid * mid.tan() < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let z = 0.5 * (lo + hi);
    let omega = z * SOUND;
    let amplitude = volume * 1e-4;
    let q0 = omega * amplitude;
    let duct = native(DgBoundary::IdealOpen);
    let initial = duct.project_initial_state(|_| 0.0, |x| q0 * (z * (1.0 - x)).cos() / z.cos());
    let energy = initial.stored_energy_j()
        + (REFERENCE - VAPOR) * volume
        + inventory
            .gas
            .volume_response(volume)
            .unwrap()
            .mechanical_potential_j;
    let mut model = CavitatingDuctPockets::new_with_trapped_gas(
        duct,
        initial,
        0.0,
        &[spec(8, volume)],
        &[inventory],
        limits(),
    )
    .unwrap();
    for quarter in 1..=4 {
        let phase = quarter as f64 * core::f64::consts::FRAC_PI_2;
        let report = model
            .advance(phase / omega, &numerics(), &mut || false)
            .unwrap();
        let pocket = report.pockets[0];
        near(
            pocket.volume_m3,
            volume + amplitude * phase.sin(),
            amplitude * 0.003,
        );
        near(pocket.end_inflow_m3_s, q0 * phase.cos(), q0 * 0.003);
        let response = inventory.gas.volume_response(pocket.volume_m3).unwrap();
        near(
            pocket.end_absolute_pressure_pa,
            VAPOR + response.pressure_pa,
            numerics().absolute_pressure_tolerance_pa,
        );
        near(report.total_stored_energy_j, energy, 2e-11);
        assert_eq!(pocket.collapse_events, 0);
        assert_eq!(pocket.mode, DuctCavityMode::VaporCavity);
        assert!(pocket.volume_m3 > inventory.minimum_volume_m3);
    }
    assert_eq!(model.trapped_gas(PortId(8)), Some(inventory));
}

#[test]
fn zero_original_inventory_recovers_the_exact_pure_vapor_path() {
    let duct = native(DgBoundary::IdealOpen);
    let initial =
        duct.project_initial_state(|x| -(REFERENCE - VAPOR) * (1.0 - x), |_| area() * 0.1);
    let volume = 1e-8;
    let empty = gas(8, volume, 0.0, 1e-12);
    let mut vapor = CavitatingDuctPockets::new(
        duct.clone(),
        initial.clone(),
        0.0,
        &[spec(8, volume)],
        limits(),
    )
    .unwrap();
    let mut empty_gas = CavitatingDuctPockets::new_with_trapped_gas(
        duct,
        initial,
        0.0,
        &[spec(8, volume)],
        &[empty],
        limits(),
    )
    .unwrap();
    let a = vapor.advance(5e-4, &numerics(), &mut || false).unwrap();
    let b = empty_gas.advance(5e-4, &numerics(), &mut || false).unwrap();
    assert_eq!(vapor.state(), empty_gas.state());
    assert_eq!(a.pockets, b.pockets);
    assert_eq!(a.total_stored_energy_j, b.total_stored_energy_j);
    assert_eq!(empty_gas.trapped_gas(PortId(8)), Some(empty));
}

#[test]
fn gas_stays_open_while_a_neighboring_vapor_pocket_collapses() {
    let duct = native(DgBoundary::Port(PortId(17)));
    let initial = duct.project_initial_state(|x| -(REFERENCE - VAPOR) * x, |_| area() * 0.01);
    let inventory = gas(8, 1e-6, REFERENCE - VAPOR, 2e-7);
    let mut model = CavitatingDuctPockets::new_with_trapped_gas(
        duct,
        initial,
        0.0,
        &[spec(17, 1e-10), spec(8, 1e-6)],
        &[inventory],
        limits(),
    )
    .unwrap();
    let report = model.advance(5e-5, &numerics(), &mut || false).unwrap();
    assert_eq!(report.pockets[0].port, PortId(8));
    assert_eq!(report.pockets[0].collapse_events, 0);
    assert_eq!(report.pockets[0].mode, DuctCavityMode::VaporCavity);
    assert_eq!(report.pockets[1].collapse_events, 1);
    assert_eq!(report.pockets[1].mode, DuctCavityMode::LiquidFilled);
    assert!(report.pockets[0].volume_m3 > 1e-6);
    near(
        report.pockets[0].end_absolute_pressure_pa,
        VAPOR
            + inventory
                .gas
                .volume_response(report.pockets[0].volume_m3)
                .unwrap()
                .pressure_pa,
        1e-3,
    );
    assert!(report.energy_balance_residual_j.abs() < 1e-11);
}

#[test]
fn gas_minimum_volume_source_pressure_and_cancellation_refuse_without_bounce() {
    let volume = 1e-6;
    let inventory = gas(8, volume, REFERENCE - VAPOR, 0.2 * volume);
    let duct = native(DgBoundary::IdealOpen);
    let initial = duct.project_initial_state(|_| 0.0, |_| -area() * 0.01);
    let mut model = CavitatingDuctPockets::new_with_trapped_gas(
        duct,
        initial,
        0.0,
        &[spec(8, volume)],
        &[inventory],
        limits(),
    )
    .unwrap();
    let before = model.state().clone();
    let mut polls = 0;
    assert_eq!(
        model
            .advance(1e-4, &numerics(), &mut || {
                polls += 1;
                polls == 7
            })
            .unwrap_err(),
        DuctCavityError::Cancelled
    );
    assert_eq!(model.state(), &before);

    // The refusal boundary is a separate trajectory; do not expect a late
    // cancellation request to supersede an earlier physical domain refusal.
    let duct = native(DgBoundary::IdealOpen);
    let initial = duct.project_initial_state(|_| 0.0, |_| -area() * 0.01);
    let inventory = gas(8, volume, REFERENCE - VAPOR, volume * (1.0 - 1e-6));
    let mut model = CavitatingDuctPockets::new_with_trapped_gas(
        duct,
        initial,
        0.0,
        &[spec(8, volume)],
        &[inventory],
        limits(),
    )
    .unwrap();
    let before = model.state().clone();
    assert!(matches!(
        model.advance(1e-4, &numerics(), &mut || false),
        Err(DuctCavityError::OutsideDomain(_))
    ));
    assert_eq!(model.state(), &before);

    let duct = native(DgBoundary::IdealOpen);
    let initial = duct.zero_state();
    let wrong = gas(8, volume, 0.5 * (REFERENCE - VAPOR), 0.2 * volume);
    assert!(
        CavitatingDuctPockets::new_with_trapped_gas(
            duct,
            initial,
            0.0,
            &[spec(8, volume)],
            &[wrong],
            limits()
        )
        .is_err()
    );
}
