use super::*;
use fs_thermal::enthalpy::DepositionObserver;

#[derive(Clone, Default)]
struct AcceptedDurations(Vec<f64>);
impl DepositionObserver for AcceptedDurations {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.0.push(dt);
        Ok(())
    }
}

#[test]
fn g1_affine_pipe_reversal_retains_spatial_inventory_reconstruction_and_native_cfl() {
    with_cx(|cx| {
        let law = law(0.0);
        let (mut data, path) = cells(&law, 3, 3.0, 1.0);
        data.iter_mut()
            .find(|cell| cell.id == path[1])
            .unwrap()
            .specific_h_j_kg = law.cp * 100.0;
        let old = simulation(&law, &data);
        let graph = graph(cx, 3.0);
        let layout = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &old,
            12,
        )
        .unwrap()
        .with_limited_reconstruction(cx, &graph, &old, 1)
        .unwrap();
        let mut boundaries = Vec::new();
        for (node, &cell) in NODES.iter().enumerate() {
            boundaries.push(MassBoundary::Inlet {
                id: 100_000 + 2 * node as u64,
                cell,
                kg_s: 0.5 * gas().density,
                material_id: MaterialId(1),
                specific_h_j_kg: 0.0,
            });
            boundaries.push(MassBoundary::Outlet {
                id: 100_001 + 2 * node as u64,
                cell,
                kg_s: 0.5 * gas().density,
            });
        }
        let transport = layout
            .bind_affine(cx, &graph, &[2.0], &[0.0], &[-2.0], &old, boundaries, 12)
            .unwrap();
        assert_eq!(transport.node_cells(), &NODES);
        assert_eq!(transport.reconstruction().unwrap().len(), 1);
        let observer = transport
            .wrap_observer(AcceptedDurations::default())
            .unwrap();
        let (next, observer, energy) = old
            .advance_observed(
                &observer,
                0.1,
                StepConfig {
                    max_dt_s: 0.05,
                    energy_tolerance_j: 1e-6,
                },
                || Ok(()),
            )
            .unwrap();
        let advection = observer.audit();
        let durations = &observer.inner().0;
        // Independent five-volume balance on the ORIGINAL accepted durations.
        // Native advection/reconstruction combine their timestep restrictions;
        // max_dt is a ceiling, not a promise of two equal .05-second steps.
        // Each directional volume rate is .5, all volumes are one, and the
        // sole limited interior slope stays zero at the central maximum.
        let mut rise = [0.0, 0.0, 100.0, 0.0, 0.0];
        let mut external_heat = 0.0;
        for &dt in durations {
            external_heat -= dt * 0.5 * gas().density * law.cp * (rise[0] + rise[4]);
            let previous = rise;
            for i in 0..rise.len() {
                let left = if i > 0 { previous[i - 1] } else { 0.0 };
                let right = if i + 1 < rise.len() {
                    previous[i + 1]
                } else {
                    0.0
                };
                rise[i] += dt * (0.5 * left + 0.5 * right - previous[i]);
            }
        }
        for (cell, expected) in [NODES[0], path[0], path[1], path[2], NODES[1]]
            .into_iter()
            .zip(rise)
        {
            near(next.temperature_k(cell).unwrap(), 300.0 + expected, 1e-10);
        }
        near(durations.iter().sum(), 0.1, 1e-15);
        assert_eq!(advection.substeps, durations.len());
        near(energy.external_heat_j, external_heat, 1e-7);
        near(energy.stored_energy_change_j, energy.external_heat_j, 1e-6);
        near(advection.inlet_mass_kg, 0.1 * gas().density, 1e-12);
        near(advection.outlet_mass_kg, advection.inlet_mass_kg, 1e-12);
        assert!(advection.maximum_courant > 0.0);
        assert_eq!(old.time_s(), 0.0);
    });
}
