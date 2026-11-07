//! Authored nonuniform pipe geometry reaches the native explicit heat timeline.

use super::*;

#[test]
fn reconstructed_graph_moves_affine_heat_with_native_sources_in_both_directions() {
    with_cx(|cx| {
        let law = law(1e9);
        let lengths = [0.1, 0.15, 0.2, 0.125, 0.175, 0.25];
        let n = lengths.len();
        let (mut data, path) = cells(&law, n, 1.0, 0.1);
        let mut start = 0.0;
        for (cell, length) in data[2..].iter_mut().zip(lengths) {
            let end = start + length;
            cell.geometry.start_m[0] = start;
            cell.geometry.end_m[0] = end;
            cell.geometry.volume_m3 = length;
            cell.mass_kg = gas().density * length;
            start = end;
        }
        for cell in &mut data {
            let centre = 0.5 * (cell.geometry.start_m[0] + cell.geometry.end_m[0]);
            cell.specific_h_j_kg = law.datum + law.cp * (20.0 + 20.0 * centre);
            // Native external heating supplies 2 K/s independently of transport.
            cell.boundary.source_w = cell.mass_kg * law.cp * 2.0;
        }
        let source_w: f64 = data.iter().map(|cell| cell.boundary.source_w).sum();
        let initial = simulation(&law, &data);
        let graph = graph(cx, 1.0);
        let inventory = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &initial,
            32,
        )
        .unwrap()
        .with_limited_reconstruction(cx, &graph, &initial, n - 2)
        .unwrap();
        let dt = 0.01;
        for q in [0.2, -0.2, 0.0] {
            let external = if q == 0.0 {
                vec![]
            } else {
                let (inlet, outlet) = if q > 0.0 {
                    (NODES[0], NODES[1])
                } else {
                    (NODES[1], NODES[0])
                };
                streams(
                    &law,
                    q,
                    inlet,
                    outlet,
                    initial.temperature_k(inlet).unwrap() - 300.0,
                )
            };
            let binding = inventory
                .bind(cx, &graph, &[q], &initial, external, 32)
                .unwrap();
            let observer = binding.wrap_observer(()).unwrap();
            assert_eq!(observer.reconstruction().unwrap().len(), n - 2);
            let (next, energy, transport) = binding
                .advance_to(
                    cx,
                    &initial,
                    dt,
                    StepConfig {
                        // The timeline clips to dt AFTER the combined CFL.
                        // Keep its proposal above this short exact interval.
                        max_dt_s: 1.0,
                        energy_tolerance_j: 1e-5,
                    },
                )
                .unwrap();
            assert_eq!(energy.substeps, 1);
            // One explicit step is exact for an interior affine enthalpy field:
            // dT/dt = source/(m Cp) - signed_speed * dT/dx. The different cell
            // lengths and deliberately unsorted IDs exercise physical geometry.
            for &id in &path[2..n - 2] {
                let expected = initial.temperature_k(id).unwrap() + dt * (2.0 - q * 20.0);
                near(next.temperature_k(id).unwrap(), expected, 1e-8);
            }
            near(energy.residual_j, 0.0, 1e-5);
            near(
                energy.external_heat_j,
                dt * source_w + transport.external_enthalpy_j,
                1e-6,
            );
            near(transport.inlet_mass_kg, dt * gas().density * q.abs(), 1e-14);
            near(transport.outlet_mass_kg, transport.inlet_mass_kg, 1e-14);
            assert_eq!(initial.time_s(), 0.0);
        }
    });
}
