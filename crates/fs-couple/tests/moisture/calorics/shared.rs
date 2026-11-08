//! Shared finite vapor and distinct bound-water energy stores on original sources.
use super::*;
use fs_couple::moisture::atmosphere::{
    CartesianMoistureParticipant, DryBasisSorption, MoistureParticipant,
};

fn energy_controls() -> SorptionEnergyControls {
    SorptionEnergyControls {
        // Native water admission is 1e-13 kg; the original vapor energy datum
        // can convert that residual into a few microjoules.
        energy_balance_tolerance_j: 1e-5,
        ..super::energy_controls()
    }
}

fn participant<'a>(
    id: u64,
    field: &'a MoistureSlab,
    faces: &'a [Option<DryBasisSorption>; 2],
) -> MoistureParticipant<'a> {
    MoistureParticipant {
        id,
        slab: field,
        faces,
        transport: config(),
    }
}

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn spatial(material: MaterialCard, contents: Vec<f64>, cx: &Cx<'_>) -> MoistureCartesian {
    MoistureCartesian::new(
        material,
        point(0.1),
        ["moisture_diffusivity"; 3],
        MaterialPropertySelection::SingleClaimOnly,
        MoistureCartesianGeometry {
            lengths_m: [0.1; 3],
            cells: [2, 1, 1],
            frame: ContentHash([4; 32]),
            maximum_moisture: 0.3,
            maximum_cells: 2,
        },
        contents,
        cx,
    )
    .unwrap()
    .with_kirchhoff_diffusivity()
}

#[test]
fn unlike_specimens_release_heat_while_the_shared_vapor_gains_no_net_water() {
    with_cx(|cx| {
        let material_a = source([3e6; 2]);
        let material_b = source([2e6; 2]);
        let heat_a = law(&material_a, [0.0, 0.3]);
        let heat_b = law(&material_b, [0.0, 0.3]);
        let a = slab(material_a, 1.0, vec![0.05], cx);
        let b = slab(material_b, 1.0, vec![0.15], cx);
        let gas = air::atmosphere(0.4, 1.0);
        let faces = [Some(air::sorption(0.25)), Some(air::sorption(0.25))];
        let participants = [participant(42, &b, &faces), participant(7, &a, &faces)];
        let bindings = [
            SorptionEnergyBinding {
                id: 7,
                law: &heat_a,
            },
            SorptionEnergyBinding {
                id: 42,
                law: &heat_b,
            },
        ];
        let native = gas
            .advance_shared(&participants, 2.0, 2, air::controls(), cx)
            .unwrap();
        let result = gas
            .advance_shared_with_sorption_energy(
                &participants,
                &bindings,
                2.0,
                2,
                air::controls(),
                energy_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.native().atmosphere(), native.atmosphere());
        for (with_energy, original) in result.native().solids().iter().zip(native.solids()) {
            assert_eq!(with_energy.id, original.id);
            assert_eq!(
                with_energy.moisture.moisture(),
                original.moisture.moisture()
            );
        }
        // Independent one-cell backward Euler: equal opposing inventories keep
        // RH=0.4, but transfer water from the weakly bound to strongly bound solid.
        let expected_a = (0.05 + 0.08 * 0.1) / 1.08;
        let water = a.dry_mass_kg() * (expected_a - 0.05);
        near(
            result.native().solids()[0].moisture.moisture()[0],
            expected_a,
            2e-13,
        );
        near(result.native().step().water_to_solids_kg, 0.0, 2e-13);
        let energy = result.energy();
        assert_eq!(
            energy
                .specimens
                .iter()
                .map(|body| body.id)
                .collect::<Vec<_>>(),
            [7, 42]
        );
        near(energy.specimens[0].storage.water_change_kg, water, 2e-13);
        near(energy.specimens[1].storage.water_change_kg, -water, 2e-13);
        near(energy.thermostat_heat_j, -1e6 * water, 2e-7);
        assert!(energy.thermostat_heat_j < -1000.0);
        assert_eq!(energy.source_segments, 2);
        assert_eq!(
            energy.gas_internal_energy_change_j,
            native.step().gas_internal_energy_change_j
        );
        near(
            energy.gas_internal_energy_change_j + energy.bound_internal_energy_change_j
                - energy.thermostat_heat_j,
            0.0,
            energy_controls().energy_balance_tolerance_j,
        );
        let reordered = gas
            .advance_shared_with_sorption_energy(
                &[participant(7, &a, &faces), participant(42, &b, &faces)],
                &[bindings[1], bindings[0]],
                2.0,
                2,
                air::controls(),
                energy_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            reordered.energy().thermostat_heat_j.to_bits(),
            energy.thermostat_heat_j.to_bits()
        );
        for (left, right) in reordered.energy().specimens.iter().zip(&energy.specimens) {
            assert_eq!(left.id, right.id);
            assert_eq!(
                left.storage.cell_binding_energy_change_j,
                right.storage.cell_binding_energy_change_j
            );
        }
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.15]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn cartesian_sdirk_accounts_original_to_final_storage_once_for_every_specimen() {
    with_cx(|cx| {
        let material_a = source([3e6, 2.4e6]);
        let material_b = source([2e6; 2]);
        let heat_a = law(&material_a, [0.0, 0.3]);
        let heat_b = law(&material_b, [0.0, 0.3]);
        let a = spatial(material_a, vec![0.04, 0.06], cx);
        let b = spatial(material_b, vec![0.16, 0.12], cx);
        let faces = std::array::from_fn(|_| Some(air::sorption(0.25)));
        let participants = [
            CartesianMoistureParticipant {
                id: 5,
                moisture: &a,
                faces: &faces,
                transport: config(),
                nonlinear: nonlinear(),
            },
            CartesianMoistureParticipant {
                id: 9,
                moisture: &b,
                faces: &faces,
                transport: config(),
                nonlinear: nonlinear(),
            },
        ];
        let bindings = [
            SorptionEnergyBinding {
                id: 9,
                law: &heat_b,
            },
            SorptionEnergyBinding {
                id: 5,
                law: &heat_a,
            },
        ];
        let gas = air::atmosphere(0.8, 1.0);
        let native = gas
            .advance_shared_cartesian_sdirk2(&participants, 0.005, 4, air::controls(), cx)
            .unwrap();
        let result = gas
            .advance_shared_cartesian_sdirk2_with_sorption_energy(
                &participants,
                &bindings,
                0.005,
                4,
                air::controls(),
                energy_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.native().atmosphere(), native.atmosphere());
        let mut expected = 0.0;
        for (index, (before, after)) in [&a, &b]
            .into_iter()
            .zip(result.native().solids())
            .enumerate()
        {
            assert_eq!(
                after.moisture.moisture(),
                native.solids()[index].moisture.moisture()
            );
            let mass = before.dry_mass_kg() / before.moisture().len() as f64;
            for (cell, (&old, &new)) in before
                .moisture()
                .iter()
                .zip(after.moisture.moisture())
                .enumerate()
            {
                // Independent primitive of RT-q(M), integrated over the whole
                // accepted interval rather than summing both stage diagnostics.
                let heat = mass
                    * if index == 0 {
                        (rt(&gas) - 3e6) * (new - old) + 1e6 * (new * new - old * old)
                    } else {
                        (rt(&gas) - 2e6) * (new - old)
                    };
                near(
                    result.energy().specimens[index]
                        .storage
                        .cell_binding_energy_change_j[cell],
                    heat,
                    2e-9,
                );
                expected += heat;
            }
        }
        near(result.energy().thermostat_heat_j, expected, 4e-9);
        assert_eq!(result.energy().source_segments, 4);
        assert_eq!(result.energy().interval_s, [0.0, 0.005]);
        assert_eq!(
            result.energy().gas_internal_energy_change_j,
            native.step().gas_internal_energy_change_j
        );
        let single = gas
            .advance_cartesian_sdirk2_with_sorption_energy(
                &a,
                0.005,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                &heat_a,
                energy_controls(),
                cx,
            )
            .unwrap();
        let shared_single = gas
            .advance_shared_cartesian_sdirk2_with_sorption_energy(
                &participants[..1],
                &bindings[1..],
                0.005,
                2,
                air::controls(),
                energy_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            single.native().moisture().moisture(),
            shared_single.native().solids()[0].moisture.moisture()
        );
        assert_eq!(
            single.energy().thermostat_heat_j.to_bits(),
            shared_single.energy().thermostat_heat_j.to_bits()
        );

        let first_order = gas
            .advance_shared_cartesian_with_sorption_energy(
                &participants,
                &bindings,
                0.005,
                4,
                air::controls(),
                energy_controls(),
                cx,
            )
            .unwrap();
        let first_native = gas
            .advance_shared_cartesian(&participants, 0.005, 4, air::controls(), cx)
            .unwrap();
        assert_eq!(first_order.native().atmosphere(), first_native.atmosphere());
        for (a, b) in first_order
            .native()
            .solids()
            .iter()
            .zip(first_native.solids())
        {
            assert_eq!(a.moisture.moisture(), b.moisture.moisture());
        }
    });
}

#[test]
fn source_binding_and_aggregate_work_refusals_publish_no_partial_shared_endpoint() {
    with_cx(|cx| {
        let material = source([2.5e6; 2]);
        let heat = law(&material, [0.0, 0.3]);
        let narrow = law(&material, [0.0, 0.10000001]);
        let a = slab(material.clone(), 1.0, vec![0.1], cx);
        let b = slab(material, 1.0, vec![0.1], cx);
        let gas = air::atmosphere(0.8, 1.0);
        let saved_gas = gas.clone();
        let faces = [Some(air::sorption(0.25)), Some(air::sorption(0.25))];
        let participants = [participant(1, &a, &faces), participant(2, &b, &faces)];
        let binding = [
            SorptionEnergyBinding { id: 1, law: &heat },
            SorptionEnergyBinding { id: 2, law: &heat },
        ];
        let run = |binding: &[SorptionEnergyBinding<'_>], controls| {
            gas.advance_shared_with_sorption_energy(
                &participants,
                binding,
                2.0,
                2,
                air::controls(),
                controls,
                cx,
            )
        };
        assert!(matches!(
            run(
                &binding,
                SorptionEnergyControls {
                    maximum_source_segments: 1,
                    ..energy_controls()
                }
            ),
            Err(SorptionEnergyError::SourceBudget {
                required: 2,
                maximum: 1
            })
        ));
        assert!(matches!(
            run(
                &[
                    binding[0],
                    SorptionEnergyBinding {
                        id: 2,
                        law: &narrow
                    },
                ],
                energy_controls()
            ),
            Err(SorptionEnergyError::Moisture(MoistureError::Material(_)))
        ));
        assert!(matches!(
            run(&[binding[0], binding[0]], energy_controls()),
            Err(SorptionEnergyError::Invalid(_))
        ));
        assert!(matches!(
            run(
                &[binding[0], SorptionEnergyBinding { id: 3, law: &heat },],
                energy_controls()
            ),
            Err(SorptionEnergyError::Invalid(_))
        ));
        let accepted = run(&binding, energy_controls()).unwrap();
        assert_eq!(accepted.energy().source_segments, 2);
        assert_eq!(gas, saved_gas);
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.1]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}
