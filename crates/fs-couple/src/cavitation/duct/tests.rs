//! Manufactured inviscid water coefficients test the coupled equations, not
//! source-data accuracy. Every duct below is the production native DG owner.
use super::*;
use fs_duct::dg::DgMesh;
use fs_duct::{Duct, Segment};
use fs_material::cavitation::SaturatedVaporState;

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 1.0;
const RHO: f64 = 1000.0;
const VAPOR: f64 = 2000.0;
const REFERENCE: f64 = 100_000.0;
const PORT: PortId = PortId(17);

fn material() -> CavitationState {
    CavitationState::try_new(293.15, RHO, 0.0, 0.0, VAPOR, 1500.0).unwrap()
}
fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}
fn limits() -> DuctCavityLimits {
    DuctCavityLimits {
        reference_pressure_pa: REFERENCE,
        maximum_cavity_volume_m3: 0.01 * area() * LENGTH,
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
        maximum_step_s: 2e-5,
        event_time_tolerance_s: 1e-13,
        event_scan_substeps: 8,
        maximum_attempts: 20_000,
        maximum_events: 32,
    }
}
fn duct(
    material: CavitationState,
    outlet_cavity: bool,
) -> DgDuct<DgFluidMaterial<CavitationState>> {
    DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &material,
        DgMesh {
            elements: 2,
            order: 2,
        },
        if outlet_cavity {
            DgBoundary::IdealOpen
        } else {
            DgBoundary::Port(PORT)
        },
        if outlet_cavity {
            DgBoundary::Port(PORT)
        } else {
            DgBoundary::IdealOpen
        },
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tol {tolerance:e})"
    );
}

#[test]
fn liquid_filled_limit_reproduces_native_wave_propagation() {
    let native = duct(material(), false);
    let initial = native.project_initial_state(
        |x| 1000.0 * (0.5 * core::f64::consts::PI * x / LENGTH).cos(),
        |_| 0.0,
    );
    let dt = 1e-4;
    let propagator = native.propagator(dt, 1e-11).unwrap();
    let mut exact = initial.clone();
    propagator.step(&mut exact, &[]).unwrap();
    propagator.sync_nodal(&mut exact).unwrap();
    let mut model = CavitatingDuct::new(native, initial, 0.0, 0.0, limits()).unwrap();
    let mut accuracy = numerics();
    accuracy.absolute_energy_tolerance_j = 1e-14;
    let report = model.advance(dt, &accuracy, &mut || false).unwrap();
    assert_eq!(report.mode, DuctCavityMode::LiquidFilled);
    assert_eq!(report.rupture_events + report.collapse_events, 0);
    assert_eq!(report.cavity_volume_m3, 0.0);
    assert_eq!(report.duct.ledger.port_work_step_j, 0.0);
    for (actual, expected) in model.state().duct().pressure().iter().zip(exact.pressure()) {
        near(*actual, *expected, 2e-6);
    }
    for (actual, expected) in model
        .state()
        .duct()
        .volume_flow()
        .iter()
        .zip(exact.volume_flow())
    {
        near(*actual, *expected, 2e-12);
    }
    assert!(report.energy_balance_residual_j.abs() < 1e-11);
}

#[test]
fn matched_boundary_loss_is_integrated_from_original_generator() {
    let native = DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &material(),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PORT),
        DgBoundary::Matched,
    )
    .unwrap();
    let initial = native.project_initial_state(|_| 500.0, |_| 0.0);
    let energy = initial.stored_energy_j();
    let mut model = CavitatingDuct::new(native, initial, 0.0, 0.0, limits()).unwrap();
    let report = model.advance(5e-4, &numerics(), &mut || false).unwrap();
    assert!(report.duct.ledger.cumulative_dissipated_j > 0.01 * energy);
    near(
        report.total_stored_energy_j + report.duct.ledger.cumulative_dissipated_j,
        energy,
        1e-11,
    );
    assert_eq!(report.duct.ledger.port_work_step_j, 0.0);
}

#[test]
fn exact_accelerating_column_conserves_volume_and_work_at_either_endpoint() {
    for outlet in [false, true] {
        let native = duct(material(), outlet);
        let speed = 0.1;
        let q0 = area() * speed;
        let acceleration = area() * (REFERENCE - VAPOR) / (RHO * LENGTH);
        let initial_volume = 1e-8;
        let collapse_time =
            (q0 + (q0 * q0 + 2.0 * acceleration * initial_volume).sqrt()) / acceleration;
        let initial = native.project_initial_state(
            |x| -(REFERENCE - VAPOR) * if outlet { x / LENGTH } else { 1.0 - x / LENGTH },
            |_| if outlet { -q0 } else { q0 },
        );
        let initial_energy = initial.stored_energy_j() + (REFERENCE - VAPOR) * initial_volume;
        let mut model =
            CavitatingDuct::new(native, initial, 0.0, initial_volume, limits()).unwrap();
        let t = 0.4 * collapse_time;
        let expected_q = q0 - acceleration * t;
        let expected_volume = initial_volume + q0 * t - 0.5 * acceleration * t * t;
        let report = model.advance(t, &numerics(), &mut || false).unwrap();
        assert_eq!(report.mode, DuctCavityMode::VaporCavity);
        near(report.end_inflow_m3_s, expected_q, 2e-12);
        near(report.cavity_volume_m3, expected_volume, 1e-14);
        near(
            report.displaced_liquid_volume_m3,
            expected_volume - initial_volume,
            1e-14,
        );
        near(report.end_absolute_pressure_pa, VAPOR, 1e-6);
        near(report.total_stored_energy_j, initial_energy, 1e-11);
        near(
            report.duct.ledger.port_work_step_j,
            (VAPOR - REFERENCE) * (expected_volume - initial_volume),
            1e-11,
        );
        let original_port = &report.duct.port_outputs[0];
        near(original_port.conjugate_rate, VAPOR - REFERENCE, 1e-6);
        near(
            original_port.end_conjugate + REFERENCE,
            report.end_absolute_pressure_pa,
            1e-8,
        );
        for q in model.state().duct().volume_flow() {
            near(*q, if outlet { -expected_q } else { expected_q }, 2e-12);
        }
        let collapsed = model
            .advance(1.01 * collapse_time, &numerics(), &mut || false)
            .unwrap();
        assert_eq!(collapsed.collapse_events, 1);
        assert_eq!(collapsed.mode, DuctCavityMode::LiquidFilled);
        near(
            collapsed.last_collapse_time_s.unwrap(),
            collapse_time,
            2e-11,
        );
        assert!(collapsed.cavity_volume_m3.abs() < numerics().absolute_volume_tolerance_m3);
        assert!(collapsed.end_absolute_pressure_pa >= VAPOR);
        near(collapsed.total_stored_energy_j, initial_energy, 2e-11);
        near(
            model.duct().port_pressures(model.state().duct()).unwrap()[0].1 + REFERENCE,
            collapsed.end_absolute_pressure_pa,
            1e-8,
        );
    }
}

#[test]
fn incident_liquid_motion_ruptures_then_feeds_back_vapor_pressure() {
    let native = duct(material(), false);
    let initial = native.project_initial_state(
        |x| -0.5 * (REFERENCE - VAPOR) * (1.0 - x / LENGTH),
        |_| area() * 0.1,
    );
    let energy = initial.stored_energy_j();
    let mut model = CavitatingDuct::new(native, initial, 0.0, 0.0, limits()).unwrap();
    let report = model.advance(2e-4, &numerics(), &mut || false).unwrap();
    assert_eq!(report.rupture_events, 1);
    assert_eq!(report.mode, DuctCavityMode::VaporCavity);
    assert!(report.last_rupture_time_s.unwrap() > 0.0);
    assert!(report.cavity_volume_m3 > 0.0);
    near(
        report.end_absolute_pressure_pa,
        VAPOR,
        numerics().absolute_pressure_tolerance_pa,
    );
    near(report.total_stored_energy_j, energy, 2e-11);
    assert!(report.duct.ledger.port_work_step_j < 0.0);
    assert!(report.duct.ledger.balance_checked);
}

#[test]
fn cancellation_attempt_limit_and_domain_failure_are_atomic() {
    let native = duct(material(), false);
    let initial = native.project_initial_state(
        |x| -(REFERENCE - VAPOR) * (1.0 - x / LENGTH),
        |_| area() * 0.1,
    );
    let mut model = CavitatingDuct::new(native, initial, 0.0, 1e-8, limits()).unwrap();
    let before = model.state().clone();
    let mut polls = 0;
    assert_eq!(
        model
            .advance(1e-3, &numerics(), &mut || {
                polls += 1;
                polls == 9
            })
            .unwrap_err(),
        DuctCavityError::Cancelled
    );
    assert_eq!(model.state(), &before);
    let mut budget = numerics();
    budget.maximum_attempts = 1;
    assert_eq!(
        model.advance(1e-3, &budget, &mut || false).unwrap_err(),
        DuctCavityError::AttemptLimit
    );
    assert_eq!(model.state(), &before);
    model.limits.maximum_cavity_volume_m3 = 1.001e-8;
    assert!(matches!(
        model.advance(1e-3, &numerics(), &mut || false),
        Err(DuctCavityError::OutsideDomain(_))
    ));
    assert_eq!(model.state(), &before);
}

#[test]
fn incompatible_source_initial_pressure_and_dense_vapor_refuse() {
    let foreign = duct(
        CavitationState::try_new(293.15, 900.0, 0.0, 0.0, VAPOR, 1500.0).unwrap(),
        false,
    );
    let native = duct(material(), false);
    assert!(matches!(
        CavitatingDuct::new(native, foreign.zero_state(), 0.0, 0.0, limits()),
        Err(DuctCavityError::Duct(_))
    ));
    let native = duct(material(), false);
    let initial = native.zero_state();
    assert!(matches!(
        CavitatingDuct::new(native, initial, 0.0, 1e-8, limits()),
        Err(DuctCavityError::InvalidInput(_))
    ));
    let native = duct(material(), false);
    let interior_void = native.project_initial_state(
        |x| -2.0 * (REFERENCE - VAPOR) * (core::f64::consts::PI * x / LENGTH).sin().powi(2),
        |_| 0.0,
    );
    assert!(matches!(
        CavitatingDuct::new(native, interior_void, 0.0, 0.0, limits()),
        Err(DuctCavityError::OutsideDomain(_))
    ));
    for (density, admitted) in [(1.0, true), (100.0, false)] {
        let source = material()
            .with_saturated_vapor(SaturatedVaporState::try_new(density, 2e6).unwrap())
            .unwrap();
        let native = duct(source, false);
        let initial = native.zero_state();
        assert_eq!(
            CavitatingDuct::new(native, initial, 0.0, 0.0, limits()).is_ok(),
            admitted
        );
    }
}
