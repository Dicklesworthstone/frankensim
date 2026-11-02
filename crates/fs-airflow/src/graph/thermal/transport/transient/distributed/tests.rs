//! Physical residence, retained pipe heat and geometry-consistent inertia.

use super::*;
use crate::graph::gravity::GraphGravitySection;
use crate::graph::inertial::InertialGraph;
use crate::graph::tests::{boundary, edge, with_cx};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::{GasSpec, GasState};
use fs_thermal::enthalpy::{
    Boundary, DepositedCell, DepositionBudget, DepositionError, Event, StepConfig, ThermalLaw,
    TimedEvent,
};

const NODES: [CellId; 2] = [CellId(1000), CellId(2000)];

struct Caloric {
    cp: f64,
    datum: f64,
    identity: u8,
}
impl ThermalLaw for Caloric {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.identity; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.datum) / self.cp;
        if t.is_finite() && (250.0..=500.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 101325.0).unwrap()
}
fn law(datum: f64) -> Caloric {
    Caloric {
        cp: gas().specific_heat_cp,
        datum,
        identity: 17,
    }
}
fn cells(
    law: &Caloric,
    n: usize,
    pipe_volume: f64,
    node_volume: f64,
) -> (Vec<DepositedCell>, Vec<CellId>) {
    let mut data = Vec::new();
    let mut path = Vec::new();
    let make = |id: CellId, start, end, volume| DepositedCell {
        id,
        segment_id: id.0,
        body_id: 0,
        material_id: MaterialId(1),
        mass_kg: gas().density * volume,
        specific_h_j_kg: law.datum,
        geometry: CellGeometry {
            start_m: [start, 0.0, 0.0],
            end_m: [end, 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    data.push(make(NODES[0], -node_volume, 0.0, node_volume));
    data.push(make(
        NODES[1],
        pipe_volume,
        pipe_volume + node_volume,
        node_volume,
    ));
    for i in 0..n {
        // Physical traversal deliberately differs from native sorted CellId order.
        let id = CellId(10 + 7 * (n - i) as u64);
        path.push(id);
        data.push(make(
            id,
            pipe_volume * i as f64 / n as f64,
            pipe_volume * (i + 1) as f64 / n as f64,
            pipe_volume / n as f64,
        ));
    }
    (data, path)
}
fn simulation<'a>(law: &'a Caloric, cells: &[DepositedCell]) -> DepositionSimulation<'a> {
    let events = cells
        .iter()
        .cloned()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: cells.len(),
            max_links: 8,
            max_events: cells.len() + 2,
            max_substeps: 2000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-5,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn graph(cx: &Cx<'_>, length: f64) -> LossGraph {
    LossGraph::new(2, vec![edge("pipe", 0, 1, 1.0)])
        .unwrap()
        .with_gravity(&gas(), &[[0.0; 3], [length, 0.0, 0.0]], [0.0; 3], cx)
        .unwrap()
}
fn streams(
    law: &Caloric,
    volume_rate: f64,
    inlet: CellId,
    outlet: CellId,
    temperature_rise: f64,
) -> Vec<MassBoundary> {
    let mass = gas().density * volume_rate.abs();
    vec![
        MassBoundary::Inlet {
            id: 100_000,
            cell: inlet,
            kg_s: mass,
            material_id: MaterialId(1),
            specific_h_j_kg: law.datum + law.cp * temperature_rise,
        },
        MassBoundary::Outlet {
            id: 100_001,
            cell: outlet,
            kg_s: mass,
        },
    ]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g0_integrated_gravity_must_follow_the_original_stored_bend_and_partition() {
    with_cx(|cx| {
        let law = law(0.0);
        let (mut data, path) = cells(&law, 2, 2.0, 1.0);
        data[2].geometry.end_m = [1.0, 0.0, 1.0];
        data[3].geometry.start_m = [1.0, 0.0, 1.0];
        let thermal = simulation(&law, &data);
        let graph = graph(cx, 2.0);
        let layout = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &thermal,
            16,
        )
        .unwrap();
        let sections: Vec<_> = path
            .iter()
            .map(|id| {
                let geometry = thermal.cells()[id].geometry;
                GraphGravitySection {
                    start_m: geometry.start_m,
                    end_m: geometry.end_m,
                    transport: gas().transport_state().unwrap(),
                }
            })
            .collect();
        let bind = |sections: Vec<GraphGravitySection>| {
            graph
                .clone()
                .with_distributed_boussinesq_gravity(
                    &gas(),
                    &[sections],
                    &[[0.0; 3], [2.0, 0.0, 0.0]],
                    [0.0, 0.0, -8.0],
                    0.2,
                    4,
                    cx,
                )
                .unwrap()
        };
        let correct = bind(sections.clone());
        layout
            .admit_gravity(cx, correct.gravity().unwrap())
            .unwrap();
        layout.admit_storage(cx, &correct, &thermal).unwrap();
        let mut refreshed = sections.clone();
        refreshed[0].transport = GasState::try_new(&GasSpec::dry_air_ussa1976(), 330.0, 101325.0)
            .unwrap()
            .transport_state()
            .unwrap();
        layout
            .admit_storage(cx, &bind(refreshed), &thermal)
            .unwrap();

        let mut other_bend = sections.clone();
        other_bend[0].end_m[2] = 2.0;
        other_bend[1].start_m[2] = 2.0;
        let changed = bind(other_bend);
        assert!(
            layout
                .admit_gravity(cx, changed.gravity().unwrap())
                .is_err()
        );
        assert!(layout.admit_storage(cx, &changed, &thermal).is_err());
        let merged = bind(vec![GraphGravitySection {
            end_m: sections[1].end_m,
            ..sections[0]
        }]);
        assert!(layout.admit_gravity(cx, merged.gravity().unwrap()).is_err());
        assert!(
            DistributedGraphInventory::new(
                cx,
                &changed,
                Density::new(gas().density),
                &NODES,
                &[path],
                &thermal,
                16
            )
            .is_err()
        );

        // A direct connection explicitly owns no interior storage to compare.
        let direct = simulation(&law, &data[..2]);
        let empty = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[vec![]],
            &direct,
            16,
        )
        .unwrap();
        empty.admit_gravity(cx, changed.gravity().unwrap()).unwrap();
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn g1_native_pipe_storage_delays_heat_and_reversal_moves_the_stored_profile() {
    with_cx(|cx| {
        let law = law(1e6);
        let volume = 1.0 / gas().density;
        let (data, path) = cells(&law, 2, 2.0 * volume, volume);
        let old = simulation(&law, &data);
        let graph = graph(cx, 2.0 * volume);
        let layout = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &old,
            16,
        )
        .unwrap();
        let q = volume; // exactly one reference kilogram per second
        let binding = layout
            .bind(
                cx,
                &graph,
                &[q],
                &old,
                streams(&law, q, NODES[0], NODES[1], 60.0),
                16,
            )
            .unwrap();
        assert_eq!(binding.node_cells(), &NODES);
        assert_eq!(layout.fluid_cell_count(), 4);
        near(
            layout.nominal_residence_times_s(cx, &[q]).unwrap()[0].unwrap(),
            2.0,
            1e-14,
        );
        let (forward, audit) = binding
            .advance_advection_implicit(cx, &old, 1.0, 1e-5, 1e-5)
            .unwrap();
        // Four equal native tanks in series. BE with dt*mdot/M=1 gives
        // successive half factors: the pipe volume is not an instant edge.
        let order = [NODES[0], path[0], path[1], NODES[1]];
        for (i, id) in order.into_iter().enumerate() {
            near(
                forward.temperature_k(id).unwrap() - 300.0,
                60.0 / 2.0_f64.powi(i as i32 + 1),
                2e-8,
            );
        }
        near(audit.stored_energy_change_j, audit.external_heat_j, 1e-5);
        assert_eq!(old.time_s(), 0.0);
        // Reverse the SAME physical pipe. New cold inlet is at the old outlet;
        // hot material already in its interior remains in the original cells.
        let reverse = layout
            .bind(
                cx,
                &graph,
                &[-q],
                &forward,
                streams(&law, q, NODES[1], NODES[0], 0.0),
                16,
            )
            .unwrap();
        let (back, audit) = reverse
            .advance_advection_implicit(cx, &forward, 2.0, 1e-5, 1e-5)
            .unwrap();
        let mut entering = 0.0;
        for id in order.into_iter().rev() {
            let old_rise = forward.temperature_k(id).unwrap() - 300.0;
            entering = 0.5 * (old_rise + entering);
            near(back.temperature_k(id).unwrap() - 300.0, entering, 3e-8);
        }
        assert!(audit.external_heat_j < 0.0);
        near(audit.stored_energy_change_j, audit.external_heat_j, 1e-5);
        // An explicitly empty path keeps the historical two-junction model.
        let direct = simulation(&law, &data[..2]);
        let empty = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[vec![]],
            &direct,
            16,
        )
        .unwrap();
        assert!(empty.inertial_branches()[0].is_none());
        assert_eq!(
            empty.nominal_residence_times_s(cx, &[q]).unwrap(),
            vec![Some(0.0)]
        );
        let direct_binding = empty
            .bind(
                cx,
                &graph,
                &[q],
                &direct,
                streams(&law, q, NODES[0], NODES[1], 60.0),
                16,
            )
            .unwrap();
        let (instant, _) = direct_binding
            .advance_advection_implicit(cx, &direct, 1.0, 1e-5, 1e-5)
            .unwrap();
        near(instant.temperature_k(NODES[1]).unwrap() - 300.0, 15.0, 2e-8);
        assert!(
            instant.temperature_k(NODES[1]).unwrap() > forward.temperature_k(NODES[1]).unwrap()
        );
    });
}

#[test]
fn g2_axial_refinement_sharpens_residence_without_changing_hydraulic_inertia_or_volume() {
    with_cx(|cx| {
        let law = law(0.0);
        let graph = graph(cx, 0.96);
        let mut outlet = Vec::new();
        let mut inertias = Vec::new();
        for n in [4, 16] {
            let (data, path) = cells(&law, n, 0.96, 0.02);
            let old = simulation(&law, &data);
            let layout = DistributedGraphInventory::new(
                cx,
                &graph,
                Density::new(gas().density),
                &NODES,
                &[path],
                &old,
                64,
            )
            .unwrap();
            let equivalent = layout.inertial_branches()[0].unwrap();
            let momentum = InertialGraph::try_new(
                cx,
                &graph,
                Density::new(gas().density),
                &[equivalent],
                &[boundary(0, 0.0)],
                &[0.0],
                0.0,
                8,
            )
            .unwrap();
            layout.admit_momentum(cx, &momentum).unwrap();
            inertias.push(momentum.inertances_pa_s2_m3()[0]);
            near(layout.branch_volumes_m3()[0], 0.96, 1e-14);
            near(
                layout.nominal_residence_times_s(cx, &[1.0]).unwrap()[0].unwrap(),
                0.96,
                1e-14,
            );
            let binding = layout
                .bind(
                    cx,
                    &graph,
                    &[1.0],
                    &old,
                    streams(&law, 1.0, NODES[0], NODES[1], 60.0),
                    64,
                )
                .unwrap();
            let config = StepConfig {
                max_dt_s: 1.0 / 512.0,
                energy_tolerance_j: 1e-5,
            };
            let (early, early_energy, _) = binding.advance_to(cx, &old, 0.5, config).unwrap();
            let (late, late_energy, _) = binding.advance_to(cx, &early, 1.5, config).unwrap();
            assert!(early_energy.residual_j.abs() < 1e-5);
            assert!(late_energy.residual_j.abs() < 1e-5);
            outlet.push((
                (early.temperature_k(NODES[1]).unwrap() - 300.0) / 60.0,
                (late.temperature_k(NODES[1]).unwrap() - 300.0) / 60.0,
            ));
        }
        near(inertias[0], inertias[1], 1e-13);
        // Fixed total volume=1 m3 and Q=1 m3/s: refining the same pipe
        // suppresses premature breakthrough and the delayed smeared tail.
        assert!(outlet[1].0 < 0.25 * outlet[0].0, "{outlet:?}");
        assert!(outlet[1].1 > outlet[0].1 + 0.05, "{outlet:?}");
        assert!(outlet[1].0 < 0.02 && outlet[1].1 > 0.95, "{outlet:?}");
        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.branches().len(), 1);
    });
}

#[test]
fn g1_variable_cross_section_inertia_comes_from_the_actual_stored_volume() {
    with_cx(|cx| {
        let law = law(0.0);
        let (mut data, path) = cells(&law, 2, 4.0, 0.02);
        data[2].geometry.end_m = [1.0, 0.0, 0.0];
        data[3].geometry.start_m = [1.0, 0.0, 0.0];
        data[3].geometry.end_m = [3.0, 0.0, 0.0];
        data[1].geometry.start_m = [3.0, 0.0, 0.0];
        data[1].geometry.end_m = [3.02, 0.0, 0.0];
        let old = simulation(&law, &data);
        let graph = graph(cx, 3.0);
        let layout = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path],
            &old,
            16,
        )
        .unwrap();
        let equivalent = layout.inertial_branches()[0].unwrap();
        // (length, area) are (1,2) then (2,1): sum(length/area)=2.5.
        near(equivalent.length_m, 3.0, 1e-14);
        near(equivalent.flow_area_m2, 1.2, 1e-14);
        let momentum = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(gas().density),
            &[equivalent],
            &[boundary(0, 0.0)],
            &[0.0],
            0.0,
            8,
        )
        .unwrap();
        near(
            momentum.inertances_pa_s2_m3()[0],
            gas().density * 2.5,
            1e-13,
        );
        layout.admit_momentum(cx, &momentum).unwrap();
        let wrong = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(gas().density),
            &[InertialBranch {
                length_m: 3.0,
                flow_area_m2: 4.0 / 3.0,
            }],
            &[boundary(0, 0.0)],
            &[0.0],
            0.0,
            8,
        )
        .unwrap();
        assert!(layout.admit_momentum(cx, &wrong).is_err());
    });
}

#[test]
fn g3_pipe_inventory_keeps_native_heating_at_zero_flow() {
    with_cx(|cx| {
        let law = law(1e6);
        let volume = 1.0 / gas().density;
        let (mut data, path) = cells(&law, 2, 2.0 * volume, volume);
        data[2].boundary.source_w = law.cp * 5.0;
        let old = simulation(&law, &data);
        let graph = graph(cx, 2.0 * volume);
        let layout = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &old,
            16,
        )
        .unwrap();
        assert_eq!(
            layout.nominal_residence_times_s(cx, &[0.0]).unwrap(),
            vec![None]
        );
        let binding = layout.bind(cx, &graph, &[0.0], &old, vec![], 16).unwrap();
        let (next, audit) = binding
            .advance_advection_implicit(cx, &old, 1.0, 1e-5, 1e-5)
            .unwrap();
        near(next.temperature_k(path[0]).unwrap(), 305.0, 2e-8);
        for id in [NODES[0], NODES[1], path[1]] {
            near(next.temperature_k(id).unwrap(), 300.0, 1e-12);
        }
        near(audit.external_heat_j, law.cp * 5.0, 1e-7);
        near(audit.stored_energy_change_j, audit.external_heat_j, 1e-5);
    });
}

#[test]
fn g0_changed_source_geometry_interior_ports_and_expanded_budgets_refuse() {
    with_cx(|cx| {
        let law = law(0.0);
        let (data, path) = cells(&law, 2, 1.0, 0.02);
        let old = simulation(&law, &data);
        let graph = graph(cx, 1.0);
        let layout = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &old,
            16,
        )
        .unwrap();
        let retained = layout.bind(cx, &graph, &[0.0], &old, vec![], 16).unwrap();
        let other = Caloric {
            cp: law.cp,
            datum: law.datum,
            identity: 18,
        };
        let changed = simulation(&other, &data);
        assert!(retained.endpoint_transport(cx, &changed).is_err());
        assert!(
            layout
                .bind(cx, &graph, &[0.0], &changed, vec![], 16)
                .is_err()
        );
        let mut broken = data.clone();
        broken[3].geometry.start_m[0] += 0.1;
        let broken = simulation(&law, &broken);
        assert!(retained.endpoint_transport(cx, &broken).is_err());
        assert!(
            DistributedGraphInventory::new(
                cx,
                &graph,
                Density::new(gas().density),
                &NODES,
                &[path.clone()],
                &broken,
                16
            )
            .is_err()
        );
        assert!(
            layout
                .bind(cx, &graph, &[0.0], &broken, vec![], 16)
                .is_err()
        );
        let interior = streams(&law, 1.0, path[0], NODES[1], 60.0);
        assert!(layout.bind(cx, &graph, &[1.0], &old, interior, 16).is_err());
        // Three interior transport segments plus two actual reservoir streams.
        assert_eq!(layout.interior_connection_count(), 3);
        assert!(
            layout
                .bind(
                    cx,
                    &graph,
                    &[1.0],
                    &old,
                    streams(&law, 1.0, NODES[0], NODES[1], 60.0),
                    4
                )
                .is_err()
        );
        assert_eq!(old.time_s(), 0.0);
        for cell in old.cells().values() {
            assert_eq!(cell.specific_h_j_kg, 0.0);
        }
    });
}
#[path = "tests/reconstruction.rs"]
mod reconstruction;
