//! Manufactured anisotropic links through the actual circuit/contact/radiation solve.
use super::*;

#[test]
fn g1_g3_material_orientation_changes_finite_heat_and_the_same_circuit_history() {
    let run = |direction: [f64; 3]| with_case(0.2, 0.2,
        |cx, _, original, mut thermal, mut circuit| {
            // No cell supplies a legacy scalar conductivity role. The same
            // original cards, caloric objects and electrical recipients remain.
            let cells = original.electrical.bulk.cells.iter().map(|cell| MaterialCell {
                cell: cell.cell, material_id: cell.material_id, curve: cell.curve,
                point: cell.point, conductivity: None,
            }).collect::<Vec<_>>();
            let bulk = MaterialNetwork { cells: &cells,
                resistors: original.electrical.bulk.resistors,
                links: original.electrical.bulk.links };
            let electrical = ElectricalContactNetwork { bulk: &bulk,
                contacts: original.electrical.contacts,
                electrical_bases: original.electrical.electrical_bases };
            let assembly = MaterialHeatNetwork { electrical: &electrical,
                thermal_contacts: original.thermal_contacts, radiation: original.radiation };
            let side = DirectionalConductionSide {
                principal_properties: ["k_principal_1", "k_principal_2", "k_principal_3"]
                    .map(|name| PropertyKey::new(name, THERMAL_CONDUCTIVITY_DIMS)),
                selection: MaterialPropertySelection::SingleClaimOnly,
                direction_in_material_axes: direction,
                material_frame: fs_blake3::ContentHash([19; 32]),
            };
            let directions = [DirectionalConductionLink { link: 7, sides: [side.clone(), side] }];
            assert!(assembly.prepare(cx, &thermal).is_err());
            let selected = assembly.prepare_with_directional_conduction(cx, &thermal, &directions).unwrap();
            let first = selected.resolve(cx, &thermal).unwrap();
            let ta = thermal.temperature_k(CellId(30)).unwrap();
            let tb = thermal.temperature_k(CellId(10)).unwrap();
            let projection = direction[0].powi(2) + 4.0 * direction[1].powi(2)
                + 9.0 * direction[2].powi(2);
            let ka = projection * (1.0 + (ta - 300.0) / 200.0);
            let kb = projection * (1.0 + (tb - 300.0) / 200.0);
            let expected_bulk = 0.1 / (0.5 / ka + 0.5 / kb);
            let resistance = 0.2 * (1.0 + (ta - 300.0) / 200.0);
            let expected_contact = 1.0 / (1.0 / expected_bulk + resistance / 0.1);
            near(first.bulk.source.material.links[0].conductance_w_k, expected_bulk);
            near(first.bulk.conductances_w_k[&7], expected_contact);
            assert_eq!(first.contacts.len(), 1);
            assert_eq!(first.bulk.source.thermal_contacts.len(), 1);
            assert!(first.bulk.source.radiation.is_some());
            for cell in &first.bulk.source.material.cells {
                let source = cell.conductivity.as_ref().unwrap();
                assert!(source.property("thermal_conductivity").is_none());
                assert!(source.property("k_principal_1").is_some());
                assert!(source.property("k_principal_2").is_some());
                assert!(source.property("k_principal_3").is_some());
                assert!(source.property("electrical_resistivity").is_some());
            }
            // Initialize at the unchanged charge/flux coordinates with this
            // complete principal/electrical source bundle, before any motion.
            circuit = MaterialCircuit::new_with_contacts(circuit.graph().clone(),
                first.bulk.resistors, first.contacts, circuit.fluxes_wb().to_vec(),
                circuit.charges_c().to_vec(), &[]).unwrap();
            for n in 1..=8 {
                let proposal = selected.step_candidate(cx, &thermal, f64::from(n) * 0.0005,
                    &circuit, &[], controls(), &[0.0; 2]).unwrap();
                near(proposal.thermal().energy_j() + proposal.report().electrical.stored_after_j, 150.0);
                let (t, c, _) = proposal.into_parts();
                thermal = t;
                circuit = c;
            }
            [thermal.temperature_k(CellId(30)).unwrap(), thermal.temperature_k(CellId(10)).unwrap(),
                circuit.charges_c()[0]]
        });
    let along = run([1.0, 0.0, 0.0]);
    let reverse = run([-1.0, 0.0, 0.0]);
    let across = run([0.0, 1.0, 0.0]);
    let rotated = run([core::f64::consts::FRAC_1_SQRT_2, core::f64::consts::FRAC_1_SQRT_2, 0.0]);
    assert_eq!(along, reverse);
    assert!(across[0] < rotated[0] && rotated[0] < along[0]);
    assert!(across[1] > rotated[1] && rotated[1] > along[1]);
    assert!((across[2] - along[2]).abs() > 1e-12);
}
