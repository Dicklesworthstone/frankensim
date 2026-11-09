//! Real axial pipe inventories without artificial junction storage.

use super::super::transient::{
    TransientRegionalCarreauControls, TransientRegionalCarreauTransport,
};
use super::*;
use crate::graph::rheology::transient::{
    TransientCarreauControls, TransientCarreauGraph, TransientCarreauStep,
};

const PATH: [CellId; 4] = [CellId(0), CellId(20), CellId(10), CellId(2)];

fn cells(law: &Caloric, temperatures: [f64; 4], split: bool) -> Vec<DepositedCell> {
    let pipe_mass = RHO * core::f64::consts::PI * RADIUS * RADIUS * LENGTH;
    let geometry = [
        ([-0.1, 0.0, 0.0], [0.0; 3], 0.04),
        (
            [0.0; 3],
            [if split { LENGTH } else { 2.0 * LENGTH }, 0.0, 0.0],
            if split { pipe_mass } else { 2.0 * pipe_mass },
        ),
        ([2.0 * LENGTH, 0.0, 0.0], [LENGTH, 0.0, 0.0], pipe_mass),
        (
            [2.0 * LENGTH, 0.0, 0.0],
            [2.0 * LENGTH + 0.1, 0.0, 0.0],
            0.03,
        ),
    ];
    (0..4)
        .filter(|&i| split || i != 2)
        .map(|i| {
            let (start_m, end_m, mass_kg) = geometry[i];
            DepositedCell {
                id: PATH[i],
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg,
                specific_h_j_kg: law.datum + law.sensible(temperatures[i]),
                geometry: CellGeometry {
                    start_m,
                    end_m,
                    volume_m3: mass_kg / RHO,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }
        })
        .collect()
}

fn initial<'m>(law: &'m Caloric, cells: Vec<DepositedCell>) -> DepositionSimulation<'m> {
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
            max_cells: 8,
            max_links: 8,
            max_events: 8,
            max_substeps: 100,
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

fn split_owner(cx: &Cx<'_>, state: &DepositionSimulation<'_>) -> RegionalCarreauTransport {
    RegionalCarreauTransport::new_with_series_splices(
        cx,
        graph(false),
        Density::new(RHO),
        &[Some(CellId(0)), None, Some(CellId(2))],
        &[CellId(20), CellId(10)],
        state,
        8,
    )
    .unwrap()
}

// Independent serial backward-Euler control volumes. There are FOUR physical
// inventories, not the five used by the finite-junction fixture in our parent.
fn reference(
    law: &Caloric,
    old: &DepositionSimulation<'_>,
    dt: f64,
    dp: f64,
    inlet: f64,
) -> (f64, [f64; 4]) {
    let order = if dp > 0.0 { [0, 1, 2, 3] } else { [3, 2, 1, 0] };
    let march = |q: f64| {
        let mut next = [0.0; 4];
        let mut donor = law.sensible(inlet);
        for i in order {
            let cell = &old.cells()[&PATH[i]];
            let old_h = cell.specific_h_j_kg - law.datum;
            let t = if i == 1 || i == 2 {
                bisect(280.0, 360.0, |t| {
                    let h = law.sensible(t);
                    cell.mass_kg * (h - old_h)
                        - dt * (RHO * q * (donor - h) + resistance(t) * q * q)
                })
            } else {
                law.temperature(
                    (cell.mass_kg * old_h + dt * RHO * q * donor) / (cell.mass_kg + dt * RHO * q),
                )
            };
            next[i] = t;
            donor = law.sensible(t);
        }
        next
    };
    let q = bisect(0.0, 0.002, |q| {
        let t = march(q);
        q * (resistance(t[1]) + resistance(t[2])) - dp.abs()
    });
    (q.copysign(dp), march(q))
}

#[test]
fn g1_series_splice_transports_hot_front_and_reversal_with_only_real_inventory() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 1e8,
            slope: 4.0,
        };
        let mut old = initial(&law, cells(&law, [300.0, 310.0, 315.0, 298.0], true));
        let owner = split_owner(cx, &old);
        assert_eq!(owner.node_cells(), &[CellId(0), CellId(2)]);
        assert_eq!(
            owner.node_cell_bindings(),
            &[Some(CellId(0)), None, Some(CellId(2))]
        );
        assert_eq!(owner.node_cell(1), None);
        assert_eq!(owner.node_cell(2), Some(CellId(2)));
        let mass: f64 = old.cells().values().map(|cell| cell.mass_kg).sum();
        near(
            mass,
            0.07 + RHO * core::f64::consts::PI * RADIUS.powi(2) * 2.0 * LENGTH,
            1e-15,
            0.0,
        );
        for (target, dp, inlet) in [(0.5, 20000.0, 340.0), (0.75, -16000.0, 290.0)] {
            let dt = target - old.time_s();
            let (q, expected) = reference(&law, &old, dt, dp, inlet);
            let (next, report) = owner
                .advance(
                    cx,
                    &old,
                    target,
                    &pressures(dp),
                    &[],
                    hydraulics(),
                    controls(),
                    &[0.0; 4],
                    |_, _, flow| Ok::<_, ()>(exchange(&law, inlet, flow)),
                )
                .unwrap();
            assert_eq!(next.cells().len(), 4);
            assert!(!next.cells().contains_key(&CellId(1)));
            for (&id, &t) in PATH.iter().zip(&expected) {
                near(next.temperature_k(id).unwrap(), t, 4e-6, 0.0);
            }
            for (index, signed_q) in [q, -q].into_iter().enumerate() {
                let branch = &report.heat.hydraulics.branches[index];
                near(
                    branch.profile.flow().volumetric_flow_m3_s,
                    signed_q,
                    1e-12,
                    3e-8,
                );
                near(
                    branch.profile.pipe().temperature_k(),
                    expected[1 + index],
                    4e-6,
                    0.0,
                );
                near(
                    report.heat.branch_viscous_heat_w[index],
                    resistance(expected[1 + index]) * q * q,
                    2e-8,
                    4e-8,
                );
            }
            let outlet = if dp > 0.0 { 3 } else { 0 };
            let streams = RHO * q.abs() * (law.sensible(inlet) - law.sensible(expected[outlet]));
            near(report.heat.external_stream_enthalpy_w, streams, 2e-5, 4e-8);
            near(
                report.thermal.stored_energy_change_j,
                dt * (streams + dp * q),
                2e-5,
                4e-8,
            );
            assert_eq!(old.time_s(), target - dt);
            old = next;
        }
    });
}

fn radial_controls() -> TransientRegionalCarreauControls {
    let mut momentum = TransientCarreauControls {
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-20),
        relative_flow_tolerance: 2e-15,
        ..TransientCarreauControls::default()
    };
    momentum.radial.absolute_velocity_tolerance_m_s = 1e-13;
    momentum.radial.relative_velocity_tolerance = 1e-13;
    TransientRegionalCarreauControls {
        hydraulics: momentum,
        thermal: controls(),
        maximum_total_branch_solves: 10000,
        mechanical_energy_tolerance_j: 1e-10,
    }
}

fn thermostat(
    regional: &RegionalCarreauTransport,
    flow: &TransientCarreauStep,
) -> CarreauTransportExchange<()> {
    let mut boundaries = Vec::new();
    for (node, rate) in flow.node_outflows.iter().enumerate() {
        let Some(cell) = regional.node_cell(node) else {
            continue;
        };
        let kg_s = RHO * rate.value();
        if kg_s > 0.0 {
            boundaries.push(MassBoundary::Inlet {
                id: 100 + node as u64,
                cell,
                kg_s,
                material_id: MaterialId(1),
                specific_h_j_kg: 0.0,
            });
        } else if kg_s < 0.0 {
            boundaries.push(MassBoundary::Outlet {
                id: 100 + node as u64,
                cell,
                kg_s: -kg_s,
            });
        }
    }
    // Authored active cooling cancels actual local viscous heat. Uniform 300 K
    // then isolates partition-invariant original momentum and work, without
    // disguising that energy as an advective or hydraulic source.
    let sources_w = regional
        .pipe_cells()
        .iter()
        .copied()
        .zip(
            flow.branches
                .iter()
                .map(|branch| -branch.viscous_heat_j / flow.duration_s),
        )
        .collect();
    CarreauTransportExchange {
        mass_boundaries: boundaries,
        heat: StepCoefficients {
            sources_w,
            conductances_w_k: BTreeMap::new(),
            payload: (),
        }
        .into(),
    }
}

#[test]
fn g1_splitting_a_pipe_preserves_native_radial_history_work_and_mass_in_both_directions() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 0.0,
            slope: 0.0,
        };
        let mut split_thermal = initial(&law, cells(&law, [300.0; 4], true));
        let mut whole_thermal = initial(&law, cells(&law, [300.0; 4], false));
        let split = split_owner(cx, &split_thermal);
        let pipe = split.graph().branches()[0].pipe;
        let whole_graph = CarreauPipeGraph::new(
            2,
            vec![CarreauPipeBranch {
                name: "whole".into(),
                from: 0,
                to: 1,
                pipe: pipe.with_length(2.0 * LENGTH).unwrap(),
            }],
        )
        .unwrap();
        let whole = RegionalCarreauTransport::new(
            cx,
            whole_graph,
            Density::new(RHO),
            &[CellId(0), CellId(2)],
            &[CellId(20)],
            &whole_thermal,
            8,
        )
        .unwrap();
        assert_eq!(
            whole.node_cell_bindings(),
            &[Some(CellId(0)), Some(CellId(2))]
        );
        let split_momentum =
            TransientCarreauGraph::new(split.graph().clone(), RHO, &[2, 2], vec![None; 3]).unwrap();
        let whole_momentum =
            TransientCarreauGraph::new(whole.graph().clone(), RHO, &[2], vec![None; 2]).unwrap();
        let mut split_owner =
            TransientRegionalCarreauTransport::new(cx, split, split_momentum).unwrap();
        let mut whole_owner =
            TransientRegionalCarreauTransport::new(cx, whole, whole_momentum).unwrap();
        let mut split_state = split_owner.momentum().rest_state(&[None; 3], cx).unwrap();
        let mut whole_state = whole_owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let mut split_work = split_owner.momentum().workspace().unwrap();
        let mut whole_work = whole_owner.momentum().workspace().unwrap();
        for (target, dp) in [(0.02, 1000.0), (0.06, -5000.0)] {
            let split = split_owner
                .advance(
                    cx,
                    &split_thermal,
                    &split_state,
                    &mut split_work,
                    target,
                    &pressures(dp),
                    &[],
                    radial_controls(),
                    &[0.0; 4],
                    |_, _, flow| Ok::<_, ()>(thermostat(split_owner.regional(), flow)),
                )
                .unwrap();
            let whole_ports = [
                FixedPressure {
                    node: 0,
                    pressure: Pressure::new(dp),
                },
                FixedPressure {
                    node: 1,
                    pressure: Pressure::new(0.0),
                },
            ];
            let whole = whole_owner
                .advance(
                    cx,
                    &whole_thermal,
                    &whole_state,
                    &mut whole_work,
                    target,
                    &whole_ports,
                    &[],
                    radial_controls(),
                    &[0.0; 3],
                    |_, _, flow| Ok::<_, ()>(thermostat(whole_owner.regional(), flow)),
                )
                .unwrap();
            let reference = &whole.report.hydraulics;
            let actual = &split.report.hydraulics;
            assert_eq!(
                actual.branches[0]
                    .midpoint_volume_flow_m3_s
                    .is_sign_positive(),
                dp > 0.0
            );
            for (branch, sign) in [(0, 1.0), (1, -1.0)] {
                for (&u, &v) in split.momentum.branches()[branch]
                    .velocities_m_s()
                    .iter()
                    .zip(whole.momentum.branches()[0].velocities_m_s())
                {
                    near(u, sign * v, 2e-12, 2e-10);
                }
                near(
                    actual.branches[branch].midpoint_volume_flow_m3_s,
                    sign * reference.branches[0].midpoint_volume_flow_m3_s,
                    2e-15,
                    2e-10,
                );
            }
            for (a, b) in [
                (actual.viscous_heat_j, reference.viscous_heat_j),
                (
                    actual.kinetic_energy_change_j,
                    reference.kinetic_energy_change_j,
                ),
                (actual.reservoir_work_j, reference.reservoir_work_j),
            ] {
                near(a, b, 2e-12, 2e-10);
            }
            for cell in split.thermal.cells().keys() {
                near(
                    split.thermal.temperature_k(*cell).unwrap(),
                    300.0,
                    1e-8,
                    0.0,
                );
            }
            near(
                split.report.other_heat_j,
                -actual.viscous_heat_j,
                1e-12,
                1e-12,
            );
            assert!(
                split.report.combined_energy_residual_j.abs()
                    <= split.report.combined_energy_tolerance_j
            );
            split_owner = split.next;
            split_state = split.momentum;
            split_thermal = split.thermal;
            whole_owner = whole.next;
            whole_state = whole.momentum;
            whole_thermal = whole.thermal;
        }
        near(
            split_thermal
                .cells()
                .values()
                .map(|cell| cell.mass_kg)
                .sum(),
            whole_thermal
                .cells()
                .values()
                .map(|cell| cell.mass_kg)
                .sum(),
            1e-15,
            0.0,
        );
    });
}

#[test]
fn g0_splices_refuse_false_geometry_ports_and_unresolved_continuity_without_changing_inventory() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 0.0,
            slope: 0.0,
        };
        let old = initial(&law, cells(&law, [300.0; 4], true));
        let owner = split_owner(cx, &old);
        let before = old.cells().clone();
        let bind = |graph, nodes: &[Option<CellId>], state: &DepositionSimulation<'_>| {
            RegionalCarreauTransport::new_with_series_splices(
                cx,
                graph,
                Density::new(RHO),
                nodes,
                &[CellId(20), CellId(10)],
                state,
                8,
            )
        };
        assert!(bind(graph(false), &[None, None, Some(CellId(2))], &old).is_err());
        let mut bent = cells(&law, [300.0; 4], true);
        bent[2].geometry.start_m = [LENGTH, LENGTH, 0.0];
        assert!(
            bind(
                graph(false),
                &[Some(CellId(0)), None, Some(CellId(2))],
                &initial(&law, bent)
            )
            .is_err()
        );
        let mut unequal = cells(&law, [300.0; 4], true);
        let mut branches = graph(false).branches().to_vec();
        let pipe = branches[1].pipe;
        branches[1].pipe = LaminarCarreauPipe::new(
            *pipe.law(),
            300.0,
            1.1 * RADIUS,
            LENGTH,
            pipe.maximum_laminar_wall_rate_per_s(),
        )
        .unwrap();
        unequal[2].mass_kg *= 1.21;
        unequal[2].geometry.volume_m3 *= 1.21;
        assert!(
            bind(
                CarreauPipeGraph::new(3, branches).unwrap(),
                &[Some(CellId(0)), None, Some(CellId(2))],
                &initial(&law, unequal)
            )
            .is_err()
        );
        for (ports, injections) in [
            (
                vec![FixedPressure {
                    node: 1,
                    pressure: Pressure::new(0.0),
                }],
                vec![],
            ),
            (
                pressures(1000.0).to_vec(),
                vec![FlowInjection {
                    node: 1,
                    flow: VolumetricFlowRate::new(0.0),
                }],
            ),
        ] {
            assert!(matches!(
                owner.advance(
                    cx,
                    &old,
                    0.1,
                    &ports,
                    &injections,
                    hydraulics(),
                    controls(),
                    &[0.0; 4],
                    |_, _, _| -> Result<CarreauTransportExchange<()>, ()> {
                        panic!("a splice port must refuse before exchange")
                    },
                ),
                Err(ImplicitTransportError::Invalid(
                    "zero-volume Carreau splice cannot carry pressure or injection ports"
                ))
            ));
        }
        for rates in [[1e-4, -1.0001e-4], [1e-4, 1e-4], [0.0, 1e-30]] {
            assert!(matches!(
                owner
                    .series
                    .validate_flows::<()>(cx, Density::new(RHO), &rates),
                Err(ImplicitTransportError::Invalid(
                    "zero-volume Carreau splice requires native fixed-mass flow continuity"
                ))
            ));
        }
        for rates in [[1e-4, -1e-4], [-1e-4, 1e-4], [0.0, 0.0]] {
            owner
                .series
                .validate_flows::<()>(cx, Density::new(RHO), &rates)
                .unwrap();
        }
        assert_eq!(old.cells(), &before);
        assert_eq!(old.time_s(), 0.0);
    });
}
