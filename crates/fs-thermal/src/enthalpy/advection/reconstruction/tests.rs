//! G1/G3/G4: original endpoint transport and implicit thermal acceptance.

use super::*;
use crate::enthalpy::advection::{MassBoundary, MassFlow};
use crate::enthalpy::backward_euler::BackwardEulerStep;
use crate::enthalpy::{
    Boundary, DepositedCell, DepositionBudget, Event, StepConfig, ThermalLaw, TimedEvent,
};

struct Law {
    datum: f64,
    hole: bool,
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let value = h - self.datum;
        if !value.is_finite()
            || !(-100.0..=1000.0).contains(&value)
            || (self.hole && (140.0..160.0).contains(&value))
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + value)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}

fn cell(law: &Law, id: usize, h: f64, a: [f64; 3], b: [f64; 3], mass: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id as u64),
        segment_id: id as u64,
        body_id: 0,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: law.datum + h,
        geometry: CellGeometry {
            start_m: a,
            end_m: b,
            volume_m3: mass,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn prepare<'a>(law: &'a Law, cells: Vec<DepositedCell>) -> DepositionSimulation<'a> {
    let n = cells.len();
    DepositionSimulation::new(
        cells
            .into_iter()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect(),
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: n,
            max_links: 2 * n + 4,
            max_events: n,
            max_substeps: 64,
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

fn ring<'a>(law: &'a Law, h: &[f64]) -> DepositionSimulation<'a> {
    let n = h.len();
    // The authored polygon has perimeter one, with exact shared vertex values.
    // Axial distances follow its original centreline, not a chord between cells.
    let radius = 1.0 / (2.0 * n as f64 * (std::f64::consts::PI / n as f64).sin());
    let vertices: Vec<_> = (0..n)
        .map(|i| {
            let angle = std::f64::consts::TAU * i as f64 / n as f64;
            [radius * angle.cos(), radius * angle.sin(), 0.0]
        })
        .collect();
    prepare(
        law,
        h.iter()
            .enumerate()
            .map(|(i, &h)| {
                cell(
                    law,
                    i,
                    h,
                    vertices[i],
                    vertices[(i + 1) % n],
                    1.0 / n as f64,
                )
            })
            .collect(),
    )
}

fn ring_flow(n: usize, rate: f64) -> Advected<()> {
    Advected::new(
        (),
        (0..n)
            .map(|i| MassFlow {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(((i + 1) % n) as u64),
                kg_s: rate,
            })
            .collect(),
        n,
    )
    .unwrap()
}

fn ring_plan(state: &DepositionSimulation<'_>) -> LimitedEnthalpyReconstruction {
    let n = state.cells().len();
    LimitedEnthalpyReconstruction::try_new(
        state,
        &(0..n)
            .map(|i| AxialCellStencil {
                previous: CellId(((i + n - 1) % n) as u64),
                cell: CellId(i as u64),
                next: CellId(((i + 1) % n) as u64),
            })
            .collect::<Vec<_>>(),
        n,
        || Ok(()),
    )
    .unwrap()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn nonuniform_affine_traces_use_actual_distances_and_reverse_the_donor_face() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = Law { datum, hole: false };
        let state = prepare(
            &law,
            vec![
                cell(&law, 0, 1.0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),
                cell(&law, 1, 5.0, [1.0, 0.0, 0.0], [4.0, 0.0, 0.0], 3.0),
                cell(&law, 2, 10.0, [4.0, 0.0, 0.0], [6.0, 0.0, 0.0], 2.0),
            ],
        );
        let plan = LimitedEnthalpyReconstruction::try_new(
            &state,
            &[AxialCellStencil {
                previous: CellId(0),
                cell: CellId(1),
                next: CellId(2),
            }],
            1,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(plan.len(), 1);
        assert!(!plan.is_empty());
        for rate in [2.0, -2.0] {
            let (inlet, outlet, inlet_h) = if rate > 0.0 {
                (0, 2, 0.0)
            } else {
                (2, 0, 12.0)
            };
            let flow = Advected::new(
                (),
                vec![
                    MassFlow {
                        id: 0,
                        a: CellId(0),
                        b: CellId(1),
                        kg_s: rate,
                    },
                    MassFlow {
                        id: 1,
                        a: CellId(1),
                        b: CellId(2),
                        kg_s: rate,
                    },
                ],
                2,
            )
            .unwrap()
            .with_boundaries(
                vec![
                    MassBoundary::Inlet {
                        id: 2,
                        cell: CellId(inlet),
                        kg_s: rate.abs(),
                        material_id: MaterialId(1),
                        specific_h_j_kg: datum + inlet_h,
                    },
                    MassBoundary::Outlet {
                        id: 3,
                        cell: CellId(outlet),
                        kg_s: rate.abs(),
                    },
                ],
                4,
            )
            .unwrap();
            let donor = flow.endpoint_transport(&state, || Ok(())).unwrap();
            let limited = flow
                .endpoint_transport_with_reconstruction(&state, &plan, || Ok(()))
                .unwrap();
            assert_eq!(donor.heat.bulk.sources_w, limited.heat.bulk.sources_w);
            let face = limited
                .transfers
                .iter()
                .find(|face| face.a == CellId(1))
                .unwrap();
            assert_eq!(face.b, CellId(if rate > 0.0 { 2 } else { 0 }));
            // h(x)=2x: the actual face values are 8 and 2, independently of
            // unequal cell lengths. Native fluxes subtract component datum h0=1.
            near(
                face.power_w / rate.abs() + 1.0,
                if rate > 0.0 { 8.0 } else { 2.0 },
                1e-12,
            );
            let empty = LimitedEnthalpyReconstruction::try_new(&state, &[], 0, || Ok(())).unwrap();
            let original = flow
                .endpoint_transport_with_reconstruction(&state, &empty, || Ok(()))
                .unwrap();
            assert_eq!(donor.transfers, original.transfers);
            assert_eq!(donor.heat.bulk.sources_w, original.heat.bulk.sources_w);
        }
    }
}

#[test]
fn smooth_periodic_face_transport_refines_and_is_less_diffusive_than_donor_cells() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let dx = 1.0 / n as f64;
        let k = std::f64::consts::TAU;
        // Independent exact cell averages of h(s)=100+10 sin(2*pi*s).
        let averages: Vec<_> = (0..n)
            .map(|i| {
                let left = i as f64 * dx;
                100.0 + 10.0 * ((k * left).cos() - (k * (left + dx)).cos()) / (k * dx)
            })
            .collect();
        let state = ring(&law, &averages);
        let plan = ring_plan(&state);
        let flow = ring_flow(n, 1.0);
        let low = flow.endpoint_transport(&state, || Ok(())).unwrap();
        let high = flow
            .endpoint_transport_with_reconstruction(&state, &plan, || Ok(()))
            .unwrap();
        let mut low_error = 0.0;
        let mut high_error = 0.0;
        let mut rates = vec![0.0; n];
        for (i, face) in high.transfers.iter().enumerate() {
            let exact = 100.0 + 10.0 * (k * (i + 1) as f64 * dx).sin();
            let trace = face.power_w + averages[0];
            high_error += dx * (trace - exact).abs();
            low_error += dx * (low.transfers[i].power_w + averages[0] - exact).abs();
            let donor = averages[i];
            let adjacent = averages[(i + 1) % n];
            assert!(trace >= donor.min(adjacent) - 1e-12 && trace <= donor.max(adjacent) + 1e-12);
            rates[face.a.0 as usize] -= face.power_w;
            rates[face.b.0 as usize] += face.power_w;
        }
        near(rates.iter().sum(), 0.0, 1e-11);
        assert!(high.heat.bulk.sources_w.is_empty());
        assert!(
            high_error < 0.2 * low_error,
            "n={n}: high={high_error}, donor={low_error}"
        );
        errors.push(high_error);
    }
    // This checks spatial FACE reconstruction, independently of time accuracy.
    // It does not claim this rate through extrema of every coupled trajectory.
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
}

#[test]
fn original_backward_euler_accepts_bounded_limited_transport_above_explicit_cfl() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = Law { datum, hole: false };
        let old = ring(&law, &[0.0, 9.0, 0.0]);
        let plan = ring_plan(&old);
        let old_cells = old.cells().clone();
        for rate in [1.0, -1.0] {
            let flow = ring_flow(3, rate);
            let dt = 1.0 / 3.0;
            let step = BackwardEulerStep::new(&old, dt, || Ok(())).unwrap();
            // Independent BE balance at lambda=rate*dt/m=1. The maximum and
            // minimum cells have zero slope; the intermediate cell transports
            // its neighbor mean. Solving the three balances gives [1,5,3]
            // (or its reversed ordering), instead of donor-cell [9,36,18]/7.
            let expected = if rate > 0.0 {
                [1.0, 5.0, 3.0]
            } else {
                [3.0, 5.0, 1.0]
            };
            let delta: Vec<_> = expected
                .into_iter()
                .zip([0.0, 9.0, 0.0])
                .map(|(new, old)| (new - old) / 3.0)
                .collect();
            let trial = step
                .evaluate_with_transport(
                    &delta,
                    |candidate| {
                        flow.endpoint_transport_with_reconstruction(candidate, &plan, || Ok(()))
                    },
                    || Ok(()),
                )
                .unwrap();
            assert!(trial.maximum_cell_residual_j() < 1e-10);
            let (next, (), audit) = trial.accept(1e-9, 1e-9, || Ok(())).unwrap();
            near(next.time_s(), dt, 0.0);
            near(audit.external_heat_j, 0.0, 0.0);
            near(audit.stored_energy_change_j, 0.0, 1e-10);
            for (i, expected) in expected.into_iter().enumerate() {
                near(
                    next.temperature_k(CellId(i as u64)).unwrap(),
                    300.0 + expected,
                    1e-10,
                );
                let h = next.cells()[&CellId(i as u64)].specific_h_j_kg - datum;
                assert!((0.0..=9.0).contains(&h));
            }
            assert_eq!(old.cells(), &old_cells);
            assert_eq!(flow.audit(), Default::default());
        }
    }
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let sharp = ring(&law, &[0.0, 0.0, 0.0, 100.0, 100.0, 100.0]);
    let plan = ring_plan(&sharp);
    for rate in [1.0, -1.0] {
        let flow = ring_flow(6, rate);
        let low = flow.endpoint_transport(&sharp, || Ok(())).unwrap();
        let high = flow
            .endpoint_transport_with_reconstruction(&sharp, &plan, || Ok(()))
            .unwrap();
        // A discontinuity is limited to donor traces; no new extremum is made.
        assert_eq!(low.transfers, high.transfers);
    }
}

#[test]
fn original_geometry_mass_topology_and_stencil_budgets_cannot_be_substituted() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let state = ring(&law, &[0.0, 100.0, 200.0]);
    let plan = ring_plan(&state);
    let flow = ring_flow(3, 0.0);
    assert!(
        flow.endpoint_transport_with_reconstruction(&state, &plan, || Ok(()))
            .unwrap()
            .transfers
            .is_empty()
    );
    let one = AxialCellStencil {
        previous: CellId(0),
        cell: CellId(1),
        next: CellId(2),
    };
    assert!(LimitedEnthalpyReconstruction::try_new(&state, &[one], 0, || Ok(())).is_err());
    assert!(LimitedEnthalpyReconstruction::try_new(&state, &[one, one], 2, || Ok(())).is_err());
    let bad = AxialCellStencil {
        previous: CellId(1),
        cell: CellId(1),
        next: CellId(2),
    };
    assert!(LimitedEnthalpyReconstruction::try_new(&state, &[bad], 1, || Ok(())).is_err());
    for change in 0..3 {
        let mut changed = state.clone();
        let cell = changed.cells.get_mut(&CellId(1)).unwrap();
        match change {
            0 => cell.geometry.volume_m3 *= 2.0,
            1 => cell.mass_kg *= 2.0,
            _ => cell.material_id = MaterialId(2),
        }
        assert!(
            flow.endpoint_transport_with_reconstruction(&changed, &plan, || Ok(()))
                .is_err()
        );
    }
    let missing = Advected::new(
        (),
        vec![MassFlow {
            id: 0,
            a: CellId(0),
            b: CellId(1),
            kg_s: 0.0,
        }],
        1,
    )
    .unwrap();
    assert!(
        missing
            .endpoint_transport_with_reconstruction(&state, &plan, || Ok(()))
            .is_err()
    );
    let stream = flow
        .clone()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 4,
                    cell: CellId(1),
                    kg_s: 1.0,
                    material_id: MaterialId(1),
                    specific_h_j_kg: 0.0,
                },
                MassBoundary::Outlet {
                    id: 5,
                    cell: CellId(1),
                    kg_s: 1.0,
                },
            ],
            5,
        )
        .unwrap();
    assert!(
        stream
            .endpoint_transport_with_reconstruction(&state, &plan, || Ok(()))
            .is_err()
    );
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(flow.audit(), Default::default());
}

#[test]
fn reconstructed_face_source_support_and_cancellation_are_transactional() {
    let law = Law {
        datum: 0.0,
        hole: true,
    };
    let state = ring(&law, &[0.0, 100.0, 200.0]);
    let plan = ring_plan(&state);
    let flow = ring_flow(3, 1.0);
    let original = state.cells().clone();
    flow.endpoint_transport(&state, || Ok(())).unwrap();
    // Supported donor states do not authorize the reconstructed h=150 face.
    assert!(matches!(
        flow.endpoint_transport_with_reconstruction(&state, &plan, || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    let mut calls = 0;
    assert!(matches!(
        flow.endpoint_transport_with_reconstruction(&state, &plan, || {
            calls += 1;
            if calls == 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(state.cells(), &original);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(flow.audit(), Default::default());
}
