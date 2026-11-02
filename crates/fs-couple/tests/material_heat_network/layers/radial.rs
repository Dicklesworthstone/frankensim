//! Physical radial widths reach the original complete circuit/heat owner.
use super::*;

#[test]
fn radial_source_integrals_reach_coupled_heat_with_physical_widths_in_both_orders() {
    for cylinder in [false, true] {
        for reversed in [false, true] {
            let radii = if reversed {
                [0.125_f64, 0.03125]
            } else {
                [0.03125_f64, 0.125]
            };
            let contact = 0.0625;
            let length = 0.25;
            let geometry = if cylinder {
                LayerGeometry::cylindrical(radii, contact, length).unwrap()
            } else {
                LayerGeometry::spherical(radii, contact).unwrap()
            };
            let widths = radii.map(|radius| (radius - contact).abs());
            assert_ne!(widths, geometry.equivalent_distances_m());
            // Independent integral(dr/A(r)) between the original traces.
            let low = radii[0].min(radii[1]);
            let high = radii[0].max(radii[1]);
            let resistance_metric = if cylinder {
                (high / low).ln() / (2.0 * core::f64::consts::PI * length)
            } else {
                (1.0 / low - 1.0 / high) / (4.0 * core::f64::consts::PI)
            };
            with_layer_geometry(
                1.0,
                false,
                Some(geometry),
                widths,
                |cx, _, assembly, layers, mut thermal, mut circuit| {
                    let prepared = assembly
                        .prepare(cx, &thermal)
                        .unwrap()
                        .with_integrated_conduction(cx, &thermal, layers, layer_controls())
                        .unwrap();
                    let initial = prepared.resolve(cx, &thermal).unwrap();
                    near(initial.bulk.conductances_w_k[&7], 1.25 / resistance_metric);
                    assert_eq!(
                        initial.bulk.source.integrated_conduction[0].geometry,
                        geometry
                    );
                    let total = thermal.energy_j() + 12.5;
                    for step in 1..=2 {
                        let candidate = prepared
                            .step_candidate(
                                cx,
                                &thermal,
                                step as f64 * 0.0005,
                                &circuit,
                                &[],
                                controls(),
                                &[0.0; 2],
                            )
                            .unwrap();
                        let source = &candidate.report().source.material;
                        let path = &source.integrated_conduction[0];
                        let ta = candidate.thermal().temperature_k(CellId(30)).unwrap();
                        let tb = candidate.thermal().temperature_k(CellId(10)).unwrap();
                        let mean_k = 1.0 + ((ta + tb) * 0.5 - 300.0) / 200.0;
                        let conductance = mean_k / resistance_metric;
                        assert_eq!(path.geometry, geometry);
                        near(path.response.conductance_w_k, conductance);
                        near(path.applied_power_a_to_b_w, conductance * (ta - tb));
                        near(source.material.links[0].conductance_w_k, conductance);
                        near(
                            candidate.thermal().energy_j()
                                + candidate.report().electrical.stored_after_j,
                            total,
                        );
                        let (next_thermal, next_circuit, _) = candidate.into_parts();
                        thermal = next_thermal;
                        circuit = next_circuit;
                    }
                },
            );
        }
    }
}

#[test]
fn radial_resistance_coordinates_and_wrong_physical_widths_refuse() {
    for geometry in [
        LayerGeometry::cylindrical([0.03125, 0.125], 0.0625, 0.25).unwrap(),
        LayerGeometry::spherical([0.03125, 0.125], 0.0625).unwrap(),
    ] {
        for authored_distances in [
            geometry.equivalent_distances_m(),
            [0.0625, 0.0625],
            [0.0625, 0.03125],
        ] {
            with_layer_geometry(
                1.0,
                false,
                Some(geometry),
                authored_distances,
                |cx, _, assembly, layers, thermal, _| {
                    assert!(matches!(
                        assembly
                            .prepare(cx, &thermal)
                            .unwrap()
                            .with_integrated_conduction(cx, &thermal, layers, layer_controls()),
                        Err(HeatNetworkError::Invalid(
                            "integrated conduction must retain original physical layer distances",
                        ))
                    ));
                },
            );
        }
    }
}
