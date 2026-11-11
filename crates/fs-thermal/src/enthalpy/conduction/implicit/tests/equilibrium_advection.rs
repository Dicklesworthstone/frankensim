//! G1/G3/G4: prescribed throughflow steady balance and original continuous response.
use super::super::equilibrium::{ThermalOperatingControls, ThermalOperatingError};
use super::super::harmonic::ThermalImpedanceBudget;
use super::*;
use crate::enthalpy::advection::reconstruction::{AxialCellStencil, LimitedEnthalpyReconstruction};
use crate::enthalpy::advection::schedule::{AdvectionSchedule, AdvectionScheduleKnot};
use crate::enthalpy::advection::{Advected, MassBoundary, MassFlow};
use fs_math::c64::C64;

fn stationary_controls() -> ThermalOperatingControls {
    ThermalOperatingControls {
        maximum_newton_iterations: 12,
        maximum_backtracks: 12,
        maximum_residual_evaluations: 64,
        cell_tolerance_w: 1e-10,
        total_tolerance_w: 1e-10,
    }
}

fn impedance_budget() -> ThermalImpedanceBudget {
    ThermalImpedanceBudget {
        maximum_cells: 8,
        maximum_dense_bytes: 100_000,
        maximum_dense_work: 100_000,
        equilibrium_tolerance_w: 1e-9,
        relative_residual_tolerance: 1e-12,
    }
}

fn throughflow(cells: usize, signed_rate: f64, inlet_h: f64) -> Advected<()> {
    let (inlet, outlet) = if signed_rate < 0.0 {
        (cells - 1, 0)
    } else {
        (0, cells - 1)
    };
    Advected::new(
        (),
        (0..cells - 1)
            .map(|i| MassFlow {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                kg_s: signed_rate,
            })
            .collect(),
        cells + 1,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: cells as u64,
                cell: CellId(inlet as u64),
                kg_s: signed_rate.abs(),
                material_id: MaterialId(0),
                specific_h_j_kg: inlet_h,
            },
            MassBoundary::Outlet {
                id: cells as u64 + 1,
                cell: CellId(outlet as u64),
                kg_s: signed_rate.abs(),
            },
        ],
        cells + 1,
    )
    .unwrap()
}

fn complex_near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 1e-12,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn throughflow_chain_and_graph_close_original_rows_in_both_directions() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    for reversed in [false, true] {
        let (mut thermal, conduction) = setup(&law, &[450.0, 225.0], None, vec![]);
        for (i, power) in (if reversed { [30.0, 10.0] } else { [10.0, 30.0] })
            .into_iter()
            .enumerate()
        {
            thermal
                .cells
                .get_mut(&CellId(i as u64))
                .unwrap()
                .boundary
                .source_w = power;
        }
        let before = thermal.cells().clone();
        let flow = throughflow(2, if reversed { -2.0 } else { 2.0 }, 300.0);
        for graph in [false, true] {
            let system = ImplicitConductionChain::prepare(
                &thermal,
                &conduction,
                &[(MaterialId(0), &law)],
                if graph {
                    LinearTopology::Graph(ImplicitConductionLinearControls::default())
                } else {
                    LinearTopology::Chain
                },
                || Ok(()),
            )
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
            let point = system
                .solve_operating_point(
                    &thermal,
                    stationary_controls(),
                    &mut ImplicitConductionWorkspace::new(2).unwrap(),
                    || Ok(()),
                )
                .unwrap();
            // P_downstream=3*(h_upstream-h_downstream)+30=0 and
            // P_total=2*(300-h_downstream)+40=0, with G=1 and Cp=m=1.
            for (actual, expected) in point.specific_enthalpies_j_kg().iter().zip(if reversed {
                [320.0, 310.0]
            } else {
                [310.0, 320.0]
            }) {
                near(*actual, expected, 1e-10);
            }
            assert_eq!(point.report().nonlinear_iterations, 1);
            assert!(point.net_powers_w().iter().all(|p| p.abs() <= 1e-10));
        }
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}

#[test]
fn throughflow_impedance_preserves_directed_mutual_response_at_dc_and_frequency() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (mut thermal, conduction) = setup(&law, &[310.0, 320.0], None, vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 10.0;
    thermal.cells.get_mut(&CellId(1)).unwrap().boundary.source_w = 30.0;
    let flow = throughflow(2, 2.0, 300.0);
    let system =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
    let model = system
        .prepare_thermal_impedance(&thermal, impedance_budget(), || Ok(()))
        .unwrap();
    for omega in [0.0, 0.2, 4.0] {
        // The original outward tangent is [[3,-1],[-3,3]], not symmetric.
        let a = C64::new(3.0, -omega);
        let determinant = a * a - C64::new(3.0, 0.0);
        let factor = model.factor(omega, || Ok(())).unwrap();
        let left = factor.response(&[C64::ONE, C64::ZERO], || Ok(())).unwrap();
        let right = factor.response(&[C64::ZERO, C64::ONE], || Ok(())).unwrap();
        complex_near(left.energy_j[0], a / determinant);
        complex_near(left.energy_j[1], C64::new(3.0, 0.0) / determinant);
        complex_near(right.energy_j[0], C64::ONE / determinant);
        complex_near(right.energy_j[1], a / determinant);
        complex_near(left.temperature_k[1], right.temperature_k[0].scale(3.0));
        assert!(left.relative_residual <= 1e-12 && right.relative_residual <= 1e-12);
    }
}

#[test]
fn throughflow_resolves_latent_storage_without_inventing_temperature_capacity() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    let (mut thermal, conduction) = setup(&law, &[420.0, 450.0], None, vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 20.0;
    thermal.cells.get_mut(&CellId(1)).unwrap().boundary.source_w = 40.0;
    let flow = throughflow(2, 2.0, 350.0);
    let system =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
    let point = system
        .solve_operating_point(
            &thermal,
            stationary_controls(),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    near(point.specific_enthalpies_j_kg()[0], 360.0, 1e-11);
    near(point.specific_enthalpies_j_kg()[1], 380.0, 1e-11);
    assert_eq!(point.temperatures_k(), &[300.0, 300.0]);
    let model = point
        .prepare_thermal_impedance(impedance_budget(), || Ok(()))
        .unwrap();
    for omega in [0.0, 4.0] {
        let response = model
            .factor(omega, || Ok(()))
            .unwrap()
            .response(&[C64::ONE, C64::ZERO], || Ok(()))
            .unwrap();
        let a = C64::new(2.0, -omega);
        complex_near(response.energy_j[0], C64::ONE / a);
        complex_near(response.energy_j[1], C64::new(2.0, 0.0) / (a * a));
        assert_eq!(response.temperature_k, vec![C64::ZERO; 2]);
    }
}

#[test]
fn stationary_schedule_uses_frozen_rates_inlets_mass_and_clock() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (mut thermal, conduction) = setup(
        &law,
        &[250.0],
        None,
        vec![TimedEvent {
            time_s: 1.0,
            event: Event::Boundary(CellId(0), boundary(Some(900.0))),
        }],
    );
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 50.0;
    let before = thermal.cells().clone();
    let flow = throughflow(1, 2.0, 300.0);
    let schedule = AdvectionSchedule::new(
        vec![
            AdvectionScheduleKnot {
                time_s: 0.0,
                flow_multiplier: 0.5,
                inlet_enthalpies_j_kg: BTreeMap::from([(1, 400.0)]),
            },
            AdvectionScheduleKnot {
                time_s: 1.0,
                flow_multiplier: 2.0,
                inlet_enthalpies_j_kg: BTreeMap::from([(1, 600.0)]),
            },
        ],
        2,
        1,
        || Ok(()),
    )
    .unwrap();
    let system =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(()))
            .unwrap();
    let point = system
        .solve_operating_point(
            &thermal,
            stationary_controls(),
            &mut ImplicitConductionWorkspace::new(1).unwrap(),
            || Ok(()),
        )
        .unwrap();
    near(point.temperatures_k()[0], 450.0, 1e-11);
    assert_eq!(point.operating_time_s(), 0.0);
    let model = point
        .prepare_thermal_impedance(impedance_budget(), || Ok(()))
        .unwrap();
    let response = model
        .factor(3.0, || Ok(()))
        .unwrap()
        .response(&[C64::ONE], || Ok(()))
        .unwrap();
    let energy = C64::ONE / C64::new(0.5, -3.0);
    complex_near(response.energy_j[0], energy);
    complex_near(response.temperature_k[0], energy.scale(0.5));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
    assert_eq!(schedule.next_transition_time_s(0.0), Some(1.0));
}

#[test]
fn original_latent_corners_allow_steady_storage_but_refuse_two_sided_harmonics() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    for h in [300.0, 400.0, 500.0] {
        let (thermal, conduction) = setup(&law, &[h], None, vec![]);
        let flow = throughflow(1, 1.0, h);
        let system = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap()
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
        let point = system
            .solve_operating_point(
                &thermal,
                stationary_controls(),
                &mut ImplicitConductionWorkspace::new(1).unwrap(),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(point.specific_enthalpies_j_kg(), &[h]);
        let response = point.prepare_thermal_impedance(impedance_budget(), || Ok(()));
        if h == 400.0 {
            assert!(response.is_ok());
        } else {
            assert!(response.is_err());
        }
    }
}

#[test]
fn stopped_and_closed_flows_do_not_invent_an_energy_reference() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[300.0, 300.0], None, vec![]);
    let circulation = Advected::new(
        (),
        vec![
            MassFlow {
                id: 0,
                a: CellId(0),
                b: CellId(1),
                kg_s: 2.0,
            },
            MassFlow {
                id: 1,
                a: CellId(1),
                b: CellId(0),
                kg_s: 2.0,
            },
        ],
        2,
    )
    .unwrap();
    // A zero-rate inlet or outlet is refused, which keeps a stopped port from
    // anchoring energy. A stopped duct is therefore a zero interior flow alone.
    let stopped = || {
        Advected::new(
            (),
            vec![MassFlow {
                id: 0,
                a: CellId(0),
                b: CellId(1),
                kg_s: 0.0,
            }],
            2,
        )
        .unwrap()
    };
    for port in [
        MassBoundary::Inlet {
            id: 1,
            cell: CellId(0),
            kg_s: 0.0,
            material_id: MaterialId(0),
            specific_h_j_kg: 300.0,
        },
        MassBoundary::Outlet {
            id: 1,
            cell: CellId(1),
            kg_s: 0.0,
        },
    ] {
        assert!(matches!(
            stopped().with_boundaries(vec![port], 2),
            Err(DepositionError::InvalidInput(
                "invalid or duplicate mass-flow boundary"
            ))
        ));
    }
    for flow in [stopped(), circulation] {
        let system = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap()
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
        assert!(matches!(
            system.solve_operating_point(
                &thermal,
                stationary_controls(),
                &mut ImplicitConductionWorkspace::new(2).unwrap(),
                || Ok(())
            ),
            Err(ThermalOperatingError::Unanchored { .. })
        ));
        let model = system
            .prepare_thermal_impedance(&thermal, impedance_budget(), || Ok(()))
            .unwrap();
        assert!(model.factor(0.0, || Ok(())).is_err());
    }
}

#[test]
fn throughflow_budget_cancellation_and_unsupported_reconstruction_publish_no_field() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[200.0, 250.0, 300.0], None, vec![]);
    let before = thermal.cells().clone();
    let flow = throughflow(3, 1.0, 350.0);
    let system = ImplicitConductionChain::prepare(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        LinearTopology::Graph(ImplicitConductionLinearControls::default()),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(&thermal, &flow, || Ok(()))
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let mut calls = 0;
    let expected = system
        .solve_operating_point(&thermal, stationary_controls(), &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap()
        .specific_enthalpies_j_kg()
        .to_vec();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(
            system
                .solve_operating_point(&thermal, stationary_controls(), &mut workspace, || {
                    count += 1;
                    if count == stop {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
    }
    let mut bounded = stationary_controls();
    bounded.maximum_residual_evaluations = 1;
    assert!(
        system
            .solve_operating_point(&thermal, bounded, &mut workspace, || Ok(()))
            .is_err()
    );
    let retry = system
        .solve_operating_point(&thermal, stationary_controls(), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.specific_enthalpies_j_kg(), &expected);
    let reconstruction = LimitedEnthalpyReconstruction::try_new(
        &thermal,
        &[AxialCellStencil {
            previous: CellId(0),
            cell: CellId(1),
            next: CellId(2),
        }],
        3,
        || Ok(()),
    )
    .unwrap();
    let limited = ImplicitConductionChain::prepare(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        LinearTopology::Graph(ImplicitConductionLinearControls::default()),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection_reconstruction(&thermal, &flow, &reconstruction, || Ok(()))
    .unwrap();
    assert!(
        limited
            .solve_operating_point(&thermal, stationary_controls(), &mut workspace, || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(flow.audit().substeps, 0);
}
