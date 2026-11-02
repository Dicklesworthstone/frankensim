//! Larger native transport graphs against independent backward-Euler solutions.
use super::*;
use crate::graph::tests::{edge, with_cx};
use crate::graph::{GraphBranch, LossGraph};
use fs_qty::Density;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};

struct Caloric {
    datum: f64,
    curved: bool,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let excess = h - self.datum;
        if !excess.is_finite() || !(-1e-6..=80_000.0).contains(&excess) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if self.curved {
            300.0 + 10.0 * ((1.0 + excess / 1000.0).sqrt() - 1.0)
        } else {
            300.0 + excess / 1000.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        if self.curved { 0.01 } else { 0.001 }
    }
}
fn inventory<'a>(law: &'a Caloric, values: &[f64]) -> DepositionSimulation<'a> {
    let schedule = values
        .iter()
        .enumerate()
        .map(|(i, &excess)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: law.datum + excess,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: values.len(),
            max_links: 2 * values.len(),
            max_events: values.len(),
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn controls(cells: usize) -> SparseAdvectionControls {
    SparseAdvectionControls {
        maximum_cells: cells,
        maximum_entries: 8 * cells,
        restart: 32,
        maximum_cycles: 128,
        maximum_corrections: 3,
        relative_linear_tolerance: 1e-13,
        cell_tolerance_j: 2e-6,
        total_tolerance_j: 2e-5,
    }
}
fn binding(
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    closed: bool,
    q: f64,
    inlet: f64,
) -> GraphEnthalpyTransport {
    let n = thermal.cells().len();
    let count = if closed { n } else { n - 1 };
    let branches: Vec<GraphBranch> = (0..count)
        .map(|i| edge(&format!("edge-{i}"), i, (i + 1) % n, 1.0))
        .collect();
    let graph = LossGraph::new(n, branches).unwrap();
    let cells: Vec<_> = (0..n).map(|i| CellId(i as u64)).collect();
    let boundaries = if closed {
        vec![]
    } else {
        let (first, last) = if q > 0.0 { (0, n - 1) } else { (n - 1, 0) };
        vec![
            MassBoundary::Inlet {
                id: (n + 1) as u64,
                cell: cells[first],
                kg_s: q.abs(),
                material_id: MaterialId(1),
                specific_h_j_kg: inlet,
            },
            MassBoundary::Outlet {
                id: (n + 2) as u64,
                cell: cells[last],
                kg_s: q.abs(),
            },
        ]
    };
    GraphEnthalpyTransport::from_volume_flows(
        cx,
        &graph,
        &vec![q; count],
        Density::new(1.0),
        &cells,
        thermal,
        boundaries,
        n + 2,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g1_sparse_chain_above_dense_limit_preserves_native_enthalpy_and_nonlinear_caloric_law() {
    with_cx(|cx| {
        let n = 96;
        for curved in [false, true] {
            let law = Caloric { datum: 0.0, curved };
            let old = inventory(&law, &vec![0.0; n]);
            let transport = binding(cx, &old, false, 9.0, 60_000.0);
            assert!(
                transport
                    .advance_advection_implicit(cx, &old, 1.0, 1e-6, 1e-5)
                    .is_err()
            );
            let (next, report) = transport
                .advance_advection_sparse_implicit(cx, &old, 1.0, controls(n))
                .unwrap();
            // Original mass=1, dt=1 and mass rate=9: h_i = 60000*(9/10)^(i+1).
            for i in 0..n {
                let expected = 60_000.0 * 0.9_f64.powi((i + 1) as i32);
                near(
                    next.cells()[&CellId(i as u64)].specific_h_j_kg,
                    expected,
                    2e-6,
                );
                near(
                    next.temperature_k(CellId(i as u64)).unwrap(),
                    law.temperature_k(expected).unwrap(),
                    2e-8,
                );
            }
            let expected_energy = 540_000.0 * (1.0 - 0.9_f64.powi(n as i32));
            near(report.thermal.stored_energy_change_j, expected_energy, 2e-5);
            near(report.thermal.external_heat_j, expected_energy, 2e-5);
            assert!(report.matrix_entries <= 4 * n);
            assert!(!report.linear.is_empty());
            assert!(report.linear.iter().all(SolveReport::converged_euclidean));
            assert_eq!(old.time_s(), 0.0);
            assert_eq!(next.time_s(), 1.0);
        }
    });
}

#[test]
fn g1_sparse_closed_transport_preserves_energy_direction_and_formation_datum() {
    with_cx(|cx| {
        let n = 96;
        for datum in [0.0, 1e8] {
            for q in [-1.0, 1.0] {
                let law = Caloric {
                    datum,
                    curved: false,
                };
                let mut initial = vec![0.0; n];
                initial[0] = 10_000.0;
                let old = inventory(&law, &initial);
                let transport = binding(cx, &old, true, q, datum);
                let (next, report) = transport
                    .advance_advection_sparse_implicit(cx, &old, 9.0, controls(n))
                    .unwrap();
                // Periodic bidiagonal inverse: h_i = h0/(1+a) * r^distance/(1-r^N).
                let amplitude = 1000.0 / (1.0 - 0.9_f64.powi(n as i32));
                for i in 0..n {
                    let distance = if q > 0.0 || i == 0 { i } else { n - i };
                    let expected = amplitude * 0.9_f64.powi(distance as i32);
                    near(
                        next.cells()[&CellId(i as u64)].specific_h_j_kg - datum,
                        expected,
                        2e-6,
                    );
                }
                let total: f64 = next
                    .cells()
                    .values()
                    .map(|c| c.specific_h_j_kg - datum)
                    .sum();
                near(total, 10_000.0, 2e-5);
                near(report.thermal.external_heat_j, 0.0, 1e-12);
                near(report.thermal.stored_energy_change_j, 0.0, 2e-5);
                assert!(report.maximum_cell_residual_j <= controls(n).cell_tolerance_j);
            }
        }
    });
}

#[test]
fn g0_sparse_limits_refuse_without_changing_inventory() {
    with_cx(|cx| {
        let n = 96;
        let law = Caloric {
            datum: 0.0,
            curved: false,
        };
        let old = inventory(&law, &vec![0.0; n]);
        let transport = binding(cx, &old, false, 9.0, 60_000.0);
        let mut limited = controls(n);
        limited.maximum_cells = n - 1;
        assert!(matches!(
            transport.advance_advection_sparse_implicit(cx, &old, 1.0, limited),
            Err(SparseAdvectionError::Invalid(_))
        ));
        let mut limited = controls(n);
        limited.maximum_entries = 1;
        assert!(matches!(
            transport.advance_advection_sparse_implicit(cx, &old, 1.0, limited),
            Err(SparseAdvectionError::Invalid(_))
        ));
        let mut limited = controls(n);
        limited.restart = 1;
        limited.maximum_cycles = 1;
        assert!(matches!(
            transport.advance_advection_sparse_implicit(cx, &old, 1.0, limited),
            Err(SparseAdvectionError::Linear(_))
        ));
        assert_eq!(old.time_s(), 0.0);
        assert!(old.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    });
}
