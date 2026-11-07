//! Prescribed pump ramps reach all original graph thermal owners.

use super::*;
use crate::graph::thermal::transport::transient::implicit::SparseAdvectionControls;
use crate::graph::thermal::transport::transient::nonlinear::ImplicitTransportControls;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::advection::schedule::AdvectionScheduleKnot;
use fs_thermal::enthalpy::backward_euler::{BoundaryStepCoefficients, StepCoefficients};

fn history(datum: f64) -> AdvectionSchedule {
    AdvectionSchedule::new(
        [(0.0, 0.0), (0.2, 0.0), (0.4, 1.0), (0.6, 1.0)]
            .into_iter()
            .map(|(time_s, flow_multiplier)| AdvectionScheduleKnot {
                time_s,
                flow_multiplier,
                inlet_enthalpies_j_kg: [(20, datum + 100.0)].into(),
            })
            .collect(),
        4,
        1,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_scheduled_explicit_graph_resolves_pump_start_and_two_tank_residence_time() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let fixed = binding(cx, &open_flow(cx, false), &old, streams(0.0));
        let scheduled = fixed.with_schedule(cx, &old, history(0.0)).unwrap();
        assert!(scheduled.schedule().is_some());
        assert!(scheduled.wrap_observer(()).unwrap().schedule().is_some());
        let mut errors = Vec::new();
        for dt in [0.02, 0.01] {
            let (next, energy, transport) = scheduled
                .advance_to(
                    cx,
                    &old,
                    0.6,
                    StepConfig {
                        max_dt_s: dt,
                        energy_tolerance_j: 1e-9,
                    },
                )
                .unwrap();
            // Operational residence time is integral r(t) dt = 0.3. The
            // independent exact solution for equal serial tanks and h_in=100
            // is h0=100(1-exp(-tau)), h1=100(1-(1+tau)exp(-tau)).
            let decay = (-0.3_f64).exp();
            let expected = [100.0 * (1.0 - decay), 100.0 * (1.0 - 1.3 * decay)];
            let actual = [
                next.cells()[&CellId(0)].specific_h_j_kg,
                next.cells()[&CellId(1)].specific_h_j_kg,
            ];
            errors.push((actual[0] - expected[0]).abs() + (actual[1] - expected[1]).abs());
            near(energy.stored_energy_change_j, actual[0] + actual[1], 1e-9);
            near(energy.external_heat_j, transport.external_enthalpy_j, 1e-9);
            near(energy.residual_j, 0.0, 1e-9);
            near(transport.inlet_mass_kg, transport.outlet_mass_kg, 1e-14);
            assert!(transport.inlet_mass_kg > 0.25 && transport.inlet_mass_kg < 0.31);
            assert!(transport.maximum_courant <= dt);
        }
        assert!(errors[1] < 0.6 * errors[0], "errors={errors:?}");
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 0.0);
    });
}

#[test]
fn g1_dense_sparse_and_nonlinear_steps_share_the_scheduled_endpoint() {
    with_cx(|cx| {
        let law = Law {
            datum: 1e6,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let scheduled = binding(cx, &open_flow(cx, false), &old, streams(law.datum))
            .with_schedule(cx, &old, history(law.datum))
            .unwrap();
        let (stopped, stopped_energy) = scheduled
            .advance_advection_implicit(cx, &old, 0.2, 1e-8, 1e-8)
            .unwrap();
        assert_eq!(stopped.cells(), old.cells());
        assert_eq!(stopped_energy.external_heat_j, 0.0);
        let (dense, dense_energy) = scheduled
            .advance_advection_implicit(cx, &stopped, 0.4, 1e-8, 1e-8)
            .unwrap();
        let (sparse, sparse_report) = scheduled
            .advance_advection_sparse_implicit(
                cx,
                &stopped,
                0.4,
                SparseAdvectionControls {
                    maximum_cells: 2,
                    maximum_entries: 32,
                    restart: 2,
                    maximum_cycles: 8,
                    maximum_corrections: 3,
                    relative_linear_tolerance: 1e-12,
                    cell_tolerance_j: 1e-8,
                    total_tolerance_j: 1e-8,
                },
            )
            .unwrap();
        let (nonlinear, nonlinear_report) = scheduled
            .advance_implicit_with_heat(
                cx,
                &stopped,
                0.4,
                ImplicitTransportControls {
                    solver: NewtonKrylovConfig {
                        absolute_tolerance: 1e-10,
                        relative_tolerance: 0.0,
                        ..NewtonKrylovConfig::default()
                    },
                    maximum_cells: 2,
                    maximum_newton_iterations: 16,
                    maximum_endpoint_evaluations: 256,
                    energy_scale_j: 100.0,
                    finite_difference_step: 1e-5,
                    cell_tolerance_j: 1e-7,
                    total_tolerance_j: 1e-7,
                },
                &[0.0; 2],
                |_, _| {
                    Ok::<BoundaryStepCoefficients<()>, DepositionError>(
                        StepCoefficients {
                            sources_w: BTreeMap::new(),
                            conductances_w_k: BTreeMap::new(),
                            payload: (),
                        }
                        .into(),
                    )
                },
            )
            .unwrap();
        // BE at the common endpoint has r=1, h_in=100 and dt=0.2:
        // 1.2*h0=20; 1.2*h1=0.2*h0. The old inlet value 120 is unused.
        for state in [&dense, &sparse, &nonlinear] {
            near(
                state.cells()[&CellId(0)].specific_h_j_kg - law.datum,
                100.0 / 6.0,
                1e-7,
            );
            near(
                state.cells()[&CellId(1)].specific_h_j_kg - law.datum,
                100.0 / 36.0,
                1e-7,
            );
            assert_eq!(state.time_s(), 0.4);
        }
        for energy in [
            dense_energy,
            sparse_report.thermal,
            nonlinear_report.thermal,
        ] {
            near(
                energy.stored_energy_change_j,
                100.0 / 6.0 + 100.0 / 36.0,
                1e-7,
            );
            near(energy.external_heat_j, energy.stored_energy_change_j, 1e-7);
        }
        assert_eq!(stopped.time_s(), 0.2);
    });
}

#[test]
fn g0_schedule_knots_and_inlet_identity_refuse_without_advancing_graph_state() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let fixed = binding(cx, &open_flow(cx, false), &old, streams(0.0));
        let scheduled = fixed.clone().with_schedule(cx, &old, history(0.0)).unwrap();
        assert!(
            matches!(scheduled.advance_advection_implicit(cx, &old, 0.4, 1e-8, 1e-8),
            Err(TransientTransportError::Invalid(message)) if message.contains("next advection schedule knot"))
        );
        let invalid = AdvectionSchedule::new(
            vec![AdvectionScheduleKnot {
                time_s: 0.0,
                flow_multiplier: 1.0,
                inlet_enthalpies_j_kg: [(21, 100.0)].into(),
            }],
            1,
            1,
            || Ok(()),
        )
        .unwrap();
        assert!(fixed.with_schedule(cx, &old, invalid).is_err());
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 0.0);
        assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    });
}
