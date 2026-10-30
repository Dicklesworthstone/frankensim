//! Declared chemistry reaches spatial heat and the native solid response.
//! Synthetic isomers demonstrate coupling; this is not a combustion fit.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_couple::cavity_phs::piston::reactive::material_rate::ACTIVITY_RATE_SCALE_PROPERTY;
use fs_couple::cavity_phs::piston::reactive::{GasReaction, GasReactionError, GasReactionSolve};
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::solid_enthalpy::gas::reactive::*;
use fs_couple::solid_enthalpy::gas::*;
use fs_couple::solid_enthalpy::*;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::PropertyKey;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergySolve, IsochoricGasParcel, parcel::ClosedGasParcel},
};
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};
use fs_solid::loads3::VolumeLoad3;
use fs_thermal::enthalpy::*;

#[path = "solid_reactive_gas/fixture.rs"]
mod fixture;
use fixture::*;

#[path = "solid_reactive_gas/coupled.rs"]
mod coupled;

#[path = "solid_reactive_gas/network.rs"]
mod network;

#[test]
fn kinetic_caloric_and_transport_data_drive_spatial_heat_deformation_and_stress() {
    with_cx(|cx, _| {
        let mut observed = Vec::new();
        for (rate, release, product_k_factor) in [
            (1.0, 20000.0, 1.0),
            (2.0, 20000.0, 1.0),
            (1.0, 40000.0, 1.0),
            (1.0, 20000.0, 8.0),
        ] {
            let law = reaction(rate, release, product_k_factor);
            let gas = gas(&law);
            let saved_gas = gas.clone();
            let (card, caloric) = wall(1.0, false);
            let thermal = initial(caloric.curve());
            let context = point(300.0);
            let model = solid(&thermal, &card, caloric.curve(), &context, cx);
            let next = model
                .step_with_reacting_isochoric_gas(
                    &thermal,
                    &gas,
                    &law,
                    1.0,
                    0.0,
                    &contacts(),
                    None,
                    chemical(),
                    solve(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();

            // Initially pure A has unit activity and no reverse contribution.
            // Equal molar masses/viscosities make the original WMS mixture law
            // reduce to mole-weighting; cp/R=3.5 gives cv/R=2.5 for both gases.
            let extent = rate * VOLUME;
            let amount =
                100000.0 * VOLUME / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 300.0);
            let gas_capacity = 2.5 * 100000.0 * VOLUME / 300.0;
            let chemical_temperature = 300.0 + release * extent / gas_capacity;
            let gas_k = 0.03 * (1.0 + extent / amount * (product_k_factor - 1.0));
            let expected = reference(chemical_temperature, gas_capacity, gas_k);
            near(next.report().reaction.extent_mol, extent, 1e-13);
            near(
                next.report().reaction.species[0].final_mol,
                amount - extent,
                1e-12,
            );
            near(next.report().reaction.species[1].initial_mol, 0.0, 1e-14);
            near(next.report().reaction.species[1].final_mol, extent, 1e-12);
            near(
                next.gas().state().state().thermal_conductivity,
                gas_k,
                1e-12,
            );
            near(next.gas().state().state().temperature, expected[2], 2e-7);
            near(
                next.gas().state().state().pressure,
                100000.0 * expected[2] / 300.0,
                1e-4,
            );
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                expected[0],
                2e-7,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                expected[1],
                2e-7,
            );
            assert!(expected[1] > expected[0] && expected[0] > 300.0);
            near(next.report().duration_s, 1.0, 0.0);
            assert!(next.report().reaction.entropy_change_j_k >= -chemical().entropy_tolerance_j_k);
            assert!(next.report().reaction.stored_energy_change_j.abs() < 2e-7);
            near(next.report().thermal.external_gas_heat_j, 0.0, 0.0);
            near(next.report().thermal.external_thermal_heat_j, 0.0, 0.0);
            assert!(
                next.report().combined_energy_residual_j.abs() < solve().total_energy_tolerance_j
            );
            assert!(next.report().thermal.solid_to_gas_heat_j < 0.0);
            near(
                next.thermal().energy_j() - thermal.energy_j()
                    + gas.mass_kg()
                        * (next.gas().state().specific_internal_energy_j_kg()
                            - gas.state().specific_internal_energy_j_kg()),
                0.0,
                2e-7,
            );

            // Solve the ORIGINAL 3-D tetrahedral problem. For nu=0 the common
            // face takes the mean in-plane strain, each apex its own strain.
            let strain = [expected[0], expected[1]].map(|t| 1e-5 * (t - 300.0));
            let loads = model
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let free = next.state().solve(&loads, &FREE, config(), cx).unwrap();
            let displacement = free.solution.displacement.displacement_m();
            near(displacement[1][0], 0.5 * (strain[0] + strain[1]), 1e-10);
            near(displacement[2][1], 0.5 * (strain[0] + strain[1]), 1e-10);
            near(displacement[3][2], strain[0], 1e-10);
            near(displacement[4][2], -strain[1], 1e-10);
            let fixed: Vec<_> = (0..15).collect();
            let clamped = next.state().solve(&loads, &fixed, config(), cx).unwrap();
            for (element, expansion) in clamped.solution.stress.elements().iter().zip(strain) {
                for stress in &element.stress_mandel_pa()[..3] {
                    near(*stress, -1000.0 * expansion, 1e-8);
                }
            }
            assert_eq!(next.gas().mass_kg(), gas.mass_kg());
            assert_eq!(next.gas().volume_m3(), gas.volume_m3());
            assert_eq!(next.gas().point().axes().get("catalyst-lot"), Some(&7.0));
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(gas, saved_gas);
            observed.push((
                expected[1] - expected[0],
                displacement[3][2],
                clamped.solution.stress.elements()[0].stress_mandel_pa()[0],
            ));
        }
        // Changing only a declared rate, formation datum or product transport
        // card reaches both the spatial gradient and actual solid observables.
        for changed in &observed[1..] {
            assert!((changed.0 - observed[0].0).abs() > 1e-4);
            assert!((changed.1 - observed[0].1).abs() > 1e-8);
            assert!((changed.2 - observed[0].2).abs() > 1e-5);
        }
    });
}

#[test]
fn continuation_uses_cooled_species_inventory_and_the_original_nonlinear_material_network() {
    with_cx(|cx, _| {
        let law = reaction(1.5, 20000.0, 4.0);
        let original_gas = gas(&law);
        let mut gas = original_gas.clone();
        let (card, caloric) = wall(1.0, true);
        let mut thermal = initial(caloric.curve());
        let initial_solid_energy = thermal.energy_j();
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = [30, 10].map(|id| MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        });
        let links = [ConductionLink {
            link: 7,
            distance_a_m: 0.03,
            distance_b_m: 0.07,
        }];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &links,
        };
        let transport = network.prepare(&thermal, || false).unwrap();
        let loads = model
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let mut previous_product = 0.0;
        let mut previous_displacement = 0.0;
        let mut previous_uncooled_temperature = 300.0;
        for interval in 1..=2 {
            let initial_temperature = gas.state().state().temperature;
            if interval == 2 {
                assert!(initial_temperature < previous_uncooled_temperature);
            }
            let expected = law
                .parcel_candidate(
                    &ClosedGasParcel::from_isochoric(gas.clone()),
                    1.0,
                    chemical(),
                    || false,
                )
                .unwrap();
            let next = model
                .step_with_reacting_isochoric_gas(
                    &thermal,
                    &gas,
                    &law,
                    interval as f64,
                    0.0,
                    &contacts(),
                    Some(&transport),
                    chemical(),
                    solve(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            // One physical chemical stage irrespective of thermal Newton work.
            assert_eq!(&next.report().reaction, expected.step());
            let rate = next.report().reaction.rate_material.as_ref().unwrap();
            near(
                rate.property(ACTIVITY_RATE_SCALE_PROPERTY)
                    .unwrap()
                    .value_si(),
                1.5 * initial_temperature / 300.0,
                1e-12,
            );
            near(
                next.report().reaction.species[1].initial_mol,
                previous_product,
                1e-12,
            );
            assert!(next.report().reaction.species[1].final_mol > previous_product);
            let temperatures =
                [CellId(30), CellId(10)].map(|id| next.thermal().temperature_k(id).unwrap());
            let sample = next.report().thermal.transport.as_ref().unwrap();
            let link = &sample.links[0];
            near(link.conductivity_w_m_k[0], temperatures[0] / 300.0, 1e-10);
            near(link.conductivity_w_m_k[1], temperatures[1] / 300.0, 1e-10);
            near(
                link.conductance_w_k,
                0.5 / (0.03 / (temperatures[0] / 300.0) + 0.07 / (temperatures[1] / 300.0)),
                1e-10,
            );
            assert!(link.conductance_w_k > 5.0); // authored timeline value was 0.5.
            for patch in &next.report().thermal.contacts {
                let internal = sample
                    .cells
                    .iter()
                    .find(|cell| cell.cell == patch.cell)
                    .unwrap()
                    .conductivity
                    .as_ref()
                    .unwrap();
                assert_eq!(internal.query(), patch.source.query());
                assert_eq!(internal.card_identity(), patch.source.card_identity());
            }
            near(
                next.thermal().energy_j() - initial_solid_energy
                    + original_gas.mass_kg()
                        * (next.gas().state().specific_internal_energy_j_kg()
                            - original_gas.state().specific_internal_energy_j_kg()),
                0.0,
                4e-7,
            );
            let equilibrium = next.state().solve(&loads, &FREE, config(), cx).unwrap();
            let displacement = equilibrium.solution.displacement.displacement_m()[3][2];
            assert!(displacement > previous_displacement);
            assert_eq!(next.gas().mass_kg(), original_gas.mass_kg());
            previous_product = next.report().reaction.species[1].final_mol;
            previous_displacement = displacement;
            previous_uncooled_temperature = expected.gas().state().state().temperature;
            let (accepted_thermal, accepted_gas, _, _) = next.into_parts();
            thermal = accepted_thermal;
            gas = accepted_gas;
        }
        assert_eq!(thermal.time_s(), 2.0);
    });
}

#[test]
fn chemistry_and_late_thermal_refusals_or_cancellation_leave_both_histories_unchanged() {
    with_cx(|cx, gate| {
        let law = reaction(1.0, 20000.0, 4.0);
        let gas = gas(&law);
        let saved_gas = gas.clone();
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let saved_cells = thermal.cells().clone();
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let patches = contacts();
        let attempt = |reaction: &GasReaction, controls: SolidGasSolve| {
            model.step_with_reacting_isochoric_gas(
                &thermal,
                &gas,
                reaction,
                1.0,
                0.0,
                &patches,
                None,
                chemical(),
                controls,
                &[0.0; 3],
                cx,
            )
        };
        let expected = law
            .parcel_candidate(
                &ClosedGasParcel::from_isochoric(gas.clone()),
                1.0,
                chemical(),
                || false,
            )
            .unwrap();
        assert!(expected.step().extent_mol > 0.0);
        // The chemistry succeeds, but the original heat solve exhausts its
        // trial budget. Its temporary reacted gas must not become accepted.
        let limited = SolidGasSolve {
            maximum_trials: 1,
            ..solve()
        };
        assert!(matches!(
            attempt(&law, limited),
            Err(ReactiveSolidGasError::Thermal(_))
        ));
        let narrow = SolidGasSolve {
            maximum_gas_temperature_increment_k: 1.0,
            ..solve()
        };
        assert!(matches!(
            attempt(&law, narrow),
            Err(ReactiveSolidGasError::Invalid(_))
        ));
        let rapid = reaction(1000.0, 20000.0, 4.0);
        assert!(matches!(
            attempt(&rapid, solve()),
            Err(ReactiveSolidGasError::Reaction(
                GasReactionError::ConsumptionLimit { .. }
            ))
        ));
        assert!(
            model
                .step_with_reacting_isochoric_gas(
                    &thermal,
                    &gas,
                    &law,
                    1.0,
                    f64::NAN,
                    &patches,
                    None,
                    chemical(),
                    solve(),
                    &[0.0; 3],
                    cx
                )
                .is_err()
        );
        assert!(
            model
                .step_with_reacting_isochoric_gas(
                    &thermal,
                    &gas,
                    &law,
                    1.0,
                    0.0,
                    &patches,
                    None,
                    chemical(),
                    solve(),
                    &[0.0; 2],
                    cx
                )
                .is_err()
        );
        // Replaying after discarded proposals advances chemistry exactly once.
        let accepted = attempt(&law, solve()).unwrap();
        let replay = attempt(&law, solve()).unwrap();
        assert_eq!(&accepted.report().reaction, expected.step());
        assert_eq!(accepted.gas(), replay.gas());
        assert_eq!(accepted.thermal().cells(), replay.thermal().cells());
        gate.request();
        assert!(matches!(
            attempt(&law, solve()),
            Err(ReactiveSolidGasError::Cancelled)
        ));
        assert_eq!(gas, saved_gas);
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
