use super::*;
use fs_duct::dg::{DgMesh, DgSidePort};
use fs_duct::{Duct, Segment};

const RADIUS: f64 = 0.01;
const RHO: f64 = 1000.0;
const REFERENCE: f64 = 100_000.0;
const VAPOR: f64 = 2000.0;
fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}
fn material() -> CavitationState {
    CavitationState::try_new(293.15, RHO, 0.0, 0.0, VAPOR, 1500.0).unwrap()
}
fn limits() -> DuctCavityLimits {
    DuctCavityLimits {
        reference_pressure_pa: REFERENCE,
        maximum_cavity_volume_m3: 2e-6,
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
fn native(
    elements: usize,
    inlet: DgBoundary,
    outlet: DgBoundary,
    sides: &[DgSidePort],
) -> DgDuct<DgFluidMaterial<CavitationState>> {
    DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: 1.0,
            }],
        },
        &material(),
        DgMesh { elements, order: 2 },
        inlet,
        outlet,
    )
    .unwrap()
    .with_side_ports(sides)
    .unwrap()
}
fn pocket(id: u32, volume: f64) -> DuctPocket {
    DuctPocket {
        port: PortId(id),
        initial_volume_m3: volume,
        maximum_volume_m3: 1e-6,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn two_end_cavities_exchange_liquid_through_a_coupled_gram_system() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(9)),
        DgBoundary::Port(PortId(2)),
        &[],
    );
    let q = area() * 0.1;
    let initial = duct.project_initial_state(|_| VAPOR - REFERENCE, |_| q);
    let v0 = 1e-7;
    let initial_energy = initial.stored_energy_j() + (REFERENCE - VAPOR) * 2.0 * v0;
    // Input declaration order does not replace the original port ordering.
    let mut model = CavitatingDuctPockets::new(
        duct,
        initial,
        0.0,
        &[pocket(2, v0), pocket(9, v0)],
        limits(),
    )
    .unwrap();
    assert!(
        model.gram[1].abs() > 0.01 * (model.gram[0] * model.gram[3]).sqrt(),
        "fixture must exercise a nondiagonal pressure Gram matrix"
    );
    let dt = 1e-3;
    let report = model.advance(dt, &numerics(), &mut || false).unwrap();
    for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
        near(report.pockets[i].end_inflow_m3_s, sign * q, 2e-12);
        near(report.pockets[i].volume_m3, v0 + sign * q * dt, 1e-14);
        near(report.pockets[i].end_absolute_pressure_pa, VAPOR, 1e-6);
        near(
            report.duct.port_outputs[i].end_conjugate + REFERENCE,
            report.pockets[i].end_absolute_pressure_pa,
            1e-8,
        );
        assert_eq!(report.pockets[i].mode, DuctCavityMode::VaporCavity);
    }
    near(
        report
            .pockets
            .iter()
            .map(|p| p.displaced_liquid_volume_m3)
            .sum(),
        0.0,
        1e-14,
    );
    near(report.duct.ledger.port_work_step_j, 0.0, 1e-11);
    near(report.total_stored_energy_j, initial_energy, 1e-11);
    assert_eq!(model.state().pockets().len(), 2);
    assert_eq!(model.state().duct().step_count(), 1);
}

#[test]
fn interior_pocket_recovers_exact_two_column_acceleration_and_collapse() {
    let duct = native(
        2,
        DgBoundary::IdealOpen,
        DgBoundary::IdealOpen,
        &[DgSidePort {
            id: PortId(4),
            downstream_element: 1,
        }],
    );
    let q = area() * 0.1;
    let initial = duct.project_initial_state(
        |x| -(REFERENCE - VAPOR) * (1.0 - 2.0 * (x - 0.5).abs()),
        |x| if x <= 0.5 { -q } else { q },
    );
    let volume = 1e-8;
    let q0 = 2.0 * q;
    let acceleration = 4.0 * area() * (REFERENCE - VAPOR) / RHO;
    let collapse = (q0 + (q0 * q0 + 2.0 * acceleration * volume).sqrt()) / acceleration;
    let energy = initial.stored_energy_j() + (REFERENCE - VAPOR) * volume;
    let mut model =
        CavitatingDuctPockets::new(duct, initial, 0.0, &[pocket(4, volume)], limits()).unwrap();
    let t = 0.4 * collapse;
    let report = model.advance(t, &numerics(), &mut || false).unwrap();
    near(
        report.pockets[0].volume_m3,
        volume + q0 * t - 0.5 * acceleration * t * t,
        1e-14,
    );
    near(
        report.pockets[0].end_inflow_m3_s,
        q0 - acceleration * t,
        2e-12,
    );
    near(report.total_stored_energy_j, energy, 1e-11);
    let final_step = model
        .advance(1.01 * collapse, &numerics(), &mut || false)
        .unwrap();
    assert_eq!(final_step.pockets[0].collapse_events, 1);
    near(
        final_step.pockets[0].last_collapse_time_s.unwrap(),
        collapse,
        2e-11,
    );
    assert_eq!(final_step.pockets[0].mode, DuctCavityMode::LiquidFilled);
    near(final_step.total_stored_energy_j, energy, 2e-11);
}

#[test]
fn simultaneous_internal_collapses_commit_together_or_not_at_all() {
    let duct = native(
        3,
        DgBoundary::IdealOpen,
        DgBoundary::IdealOpen,
        &[
            DgSidePort {
                id: PortId(101),
                downstream_element: 1,
            },
            DgSidePort {
                id: PortId(7),
                downstream_element: 2,
            },
        ],
    );
    let left = 1.0 / 3.0;
    let right = 2.0 / 3.0;
    let q = area() * 0.1;
    let initial = duct.project_initial_state(
        |x| {
            -(REFERENCE - VAPOR)
                * if x <= left {
                    x / left
                } else if x <= right {
                    1.0
                } else {
                    (1.0 - x) / (1.0 - right)
                }
        },
        |x| {
            if x <= left {
                -q
            } else if x <= right {
                0.0
            } else {
                q
            }
        },
    );
    let volume = 1e-8;
    let acceleration = 3.0 * area() * (REFERENCE - VAPOR) / RHO;
    let collapse = (q + (q * q + 2.0 * acceleration * volume).sqrt()) / acceleration;
    let energy = initial.stored_energy_j() + 2.0 * (REFERENCE - VAPOR) * volume;
    let mut model = CavitatingDuctPockets::new(
        duct,
        initial,
        0.0,
        &[pocket(7, volume), pocket(101, volume)],
        limits(),
    )
    .unwrap();
    let before = model.state().clone();
    let mut limited = numerics();
    limited.maximum_events = 1;
    assert_eq!(
        model
            .advance(1.005 * collapse, &limited, &mut || false)
            .unwrap_err(),
        DuctCavityError::EventLimit
    );
    assert_eq!(model.state(), &before);
    let report = model
        .advance(
            collapse + 4.0 * numerics().event_time_tolerance_s,
            &numerics(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(report.simultaneous_groups, 1);
    for pocket in &report.pockets {
        assert_eq!(pocket.collapse_events, 1);
        assert_eq!(pocket.mode, DuctCavityMode::LiquidFilled);
        near(pocket.last_collapse_time_s.unwrap(), collapse, 2e-11);
        assert!(pocket.volume_m3.abs() <= numerics().absolute_volume_tolerance_m3);
    }
    assert_eq!(
        report.pockets[0].last_collapse_time_s,
        report.pockets[1].last_collapse_time_s
    );
    near(report.total_stored_energy_j, energy, 2e-11);
    // This zero-pressure-margin middle column develops an undershoot under
    // the native dispersive collision. A pocket at another station would be
    // needed; the coupling refuses it instead of clipping bulk pressure.
    let at_collision = model.state().clone();
    assert!(matches!(
        model.advance(1.005 * collapse, &numerics(), &mut || false),
        Err(DuctCavityError::OutsideDomain(_))
    ));
    assert_eq!(model.state(), &at_collision);
}

#[test]
fn empty_contacts_select_both_consistent_growth_constraints() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
        &[],
    );
    let initial =
        duct.project_initial_state(|_| VAPOR - REFERENCE, |x| area() * 0.01 * (1.0 - 2.0 * x));
    let energy = initial.stored_energy_j();
    let mut model = CavitatingDuctPockets::new(
        duct,
        initial,
        0.0,
        &[pocket(0, 0.0), pocket(1, 0.0)],
        limits(),
    )
    .unwrap();
    assert!(
        model
            .state()
            .pockets()
            .iter()
            .all(|p| p.mode == DuctCavityMode::VaporCavity)
    );
    let report = model.advance(1e-5, &numerics(), &mut || false).unwrap();
    for pocket in &report.pockets {
        assert!(pocket.volume_m3 > 0.0);
        assert!(pocket.end_inflow_m3_s > 0.0);
        near(pocket.end_absolute_pressure_pa, VAPOR, 1e-6);
    }
    near(report.total_stored_energy_j, energy, 1e-11);
}

#[test]
fn internal_geometry_checks_both_liquid_interfaces_even_when_net_flow_is_zero() {
    let duct = DgDuct::from_material(
        &Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: RADIUS,
                    length: 0.5,
                },
                Segment::Cylinder {
                    radius: 2.0 * RADIUS,
                    length: 0.5,
                },
            ],
        },
        &material(),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::IdealOpen,
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: PortId(4),
        downstream_element: 1,
    }])
    .unwrap();
    let geometry = duct.side_port_geometry(PortId(4)).unwrap();
    near(geometry.station_m, 0.5, 0.0);
    near(geometry.upstream_area_m2, area(), 1e-18);
    near(geometry.downstream_area_m2, 4.0 * area(), 1e-18);
    let q = area() * material().liquid_sound_speed_m_s() * 0.02;
    let initial = duct.project_initial_state(|_| VAPOR - REFERENCE, |_| q);
    assert!(matches!(
        CavitatingDuctPockets::new(duct, initial, 0.0, &[pocket(4, 1e-8)], limits()),
        Err(DuctCavityError::OutsideDomain(_))
    ));
}
