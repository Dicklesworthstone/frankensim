//! G1/G3/G4: reacting sound and native species/enthalpy flow share one sample.

use super::super::exchange::quiet_body;
use super::super::flow::{flow_solve, nozzle};
use super::*;
use fs_couple::acoustic_realize::finite_gas::flow::FiniteGasStringFlowDirection;
use fs_couple::cavity_phs::piston::thermal::{
    HeatExchangeSolve, layers::MaterialConductionLayers,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;

fn reservoir(model: ThermochemicalGasModel, pressure: f64) -> ClosedGasParcel {
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    ClosedGasParcel::new(
        model,
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                300.0,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
            .unwrap(),
        0.001,
    )
    .unwrap()
}

fn reactant_reservoir(reaction: &GasReaction, pressure: f64) -> ClosedGasParcel {
    reservoir(
        reaction
            .catalog()
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap(),
            )
            .unwrap(),
        pressure,
    )
}

// Independent extensive inventory for the equal-molar-mass A/B fixture.
fn amounts(gas: &ClosedGasParcel) -> [f64; 2] {
    let total = gas.mass_kg() / 0.028;
    ["A", "B"].map(|name| {
        gas.model()
            .thermodynamics()
            .components()
            .iter()
            .zip(gas.model().thermodynamics().mole_composition().fractions())
            .find(|(component, _)| component.species() == &SpeciesId::new(name).unwrap())
            .map_or(0.0, |(_, fraction)| total * fraction)
    })
}

#[test]
fn g1_reacting_inlet_outlet_and_closed_nozzle_preserve_elements_energy_and_next_inventory() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        for (pressure, area) in [(200000.0, 1e-4), (50000.0, 1e-4), (200000.0, 0.0)] {
            let (mut live, mut gas) = prepared(cx, &reaction);
            let mut reservoir = reactant_reservoir(&reaction, pressure);
            let before = live.clone();
            let old_gas = gas.clone();
            let old_reservoir = reservoir.clone();
            let local_before = amounts(&ClosedGasParcel::from_isochoric(gas.clone()));
            let reservoir_before = amounts(&reservoir);
            let heat = FiniteGasStringHeat {
                string_heat_j: 5e-8,
                gas_heat_j: 0.01,
            };
            let geometry = GasNozzle {
                area_m2: area,
                ..nozzle()
            };
            // Reference: the existing chemical owner once, then the existing
            // native transfer once over exactly the same sample duration.
            let mut thermal = before.clone();
            let mut thermal_gas = old_gas.clone();
            let heated = thermal
                .step_with_reacting_isochoric_gas(
                    cx,
                    8,
                    &[0.0; 2],
                    &mut thermal_gas,
                    heat,
                    &reaction,
                    controls(),
                )
                .unwrap();
            let local = ClosedGasParcel::from_isochoric(thermal_gas.clone());
            let into = pressure >= local.state().state().pressure;
            let (donor, receiver) = if into {
                (&reservoir, &local)
            } else {
                (&local, &reservoir)
            };
            let (proposal, expected_flow) = donor
                .through_nozzle_to(receiver, geometry, DT, flow_solve(), || false)
                .unwrap();
            let (donor, receiver, expected_transfer) = proposal.into_parts();
            let (expected_gas, expected_reservoir) = if into {
                (receiver, donor)
            } else {
                (donor, receiver)
            };
            let result = live
                .step_with_reacting_isochoric_gas_nozzle(
                    cx,
                    8,
                    &[0.0; 2],
                    &mut gas,
                    &mut reservoir,
                    heat,
                    &reaction,
                    controls(),
                    geometry,
                    flow_solve(),
                )
                .unwrap();
            assert_eq!(result.thermal.chemistry, heated.chemistry);
            assert_eq!(result.flow, expected_flow);
            assert_eq!(result.transfer, expected_transfer);
            assert_eq!(gas, expected_gas.into_isochoric());
            assert_eq!(reservoir, expected_reservoir);
            assert_eq!(live.thermal(), thermal.thermal());
            assert_eq!(live.mechanical().states(), thermal.mechanical().states());
            assert_eq!(live.mechanical().ambient(), gas.state().state());
            assert_eq!(live.mechanical().epoch(), 9);
            assert_eq!(live.mechanical().accepted_samples(), 9);
            assert_eq!(
                result.direction,
                if into {
                    FiniteGasStringFlowDirection::IntoSurroundings
                } else {
                    FiniteGasStringFlowDirection::OutOfSurroundings
                }
            );
            if area == 0.0 {
                assert_eq!(result.transfer.transferred_mass_kg, 0.0);
                assert_eq!(result.thermal.coupled.string, heated.coupled.string);
                assert_eq!(gas, thermal_gas);
                assert_eq!(reservoir, old_reservoir);
            } else {
                assert!(result.transfer.transferred_mass_kg > 0.0);
            }
            let local_after = amounts(&ClosedGasParcel::from_isochoric(gas.clone()));
            let reservoir_after = amounts(&reservoir);
            let extent = result.thermal.chemistry.reaction.extent_mol;
            for (i, change) in [-extent, extent].into_iter().enumerate() {
                near(
                    local_after[i] + reservoir_after[i],
                    local_before[i] + reservoir_before[i] + change,
                    2e-12,
                );
            }
            near(
                local_after.iter().sum::<f64>() + reservoir_after.iter().sum::<f64>(),
                local_before.iter().sum::<f64>() + reservoir_before.iter().sum::<f64>(),
                2e-12,
            );
            let storage = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
                - old_gas.mass_kg() * old_gas.state().specific_internal_energy_j_kg()
                + reservoir.mass_kg() * reservoir.state().specific_internal_energy_j_kg()
                - old_reservoir.mass_kg()
                    * old_reservoir.state().specific_internal_energy_j_kg()
                + live.mechanical().specimen().mass_kg()
                    * (live.thermal().specific_enthalpy_j_kg()
                        - before.thermal().specific_enthalpy_j_kg())
                + live.mechanical_internal_energy_j().unwrap()
                - before.mechanical_internal_energy_j().unwrap();
            let vibration = &result.thermal.coupled.string.vibration;
            let residual =
                heat.string_heat_j + heat.gas_heat_j + vibration.acoustic.input_work_j
                    - storage
                    - vibration.dissipation.authored_loss_j
                    - vibration.dissipation.roundoff_residual_j;
            assert!(residual.abs() <= result.energy_tolerance_j);
            assert!(result.energy_balance_residual_j.abs() <= result.energy_tolerance_j);

            // The next chemistry solve must start from the transferred
            // inventory; reusing the pre-transfer report loses the inlet.
            let next = live
                .step_with_reacting_isochoric_gas_nozzle(
                    cx,
                    9,
                    &[0.0; 2],
                    &mut gas,
                    &mut reservoir,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    controls(),
                    geometry,
                    flow_solve(),
                )
                .unwrap();
            for (species, initial) in next
                .thermal
                .chemistry
                .reaction
                .species
                .iter()
                .zip(local_after)
            {
                near(species.initial_mol, initial, 1e-12);
            }
            assert_eq!(live.mechanical().epoch(), 10);
            assert_eq!(live.mechanical().accepted_samples(), 10);
        }
    });
}

#[test]
fn g3_reacting_conduction_keeps_its_endpoint_receipts_through_flow() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        for area in [0.0, 1e-4] {
            let mut gas = gas(&reaction);
            let mut reservoir = reactant_reservoir(&reaction, 200000.0);
            let (mut live, card, point, key) = quiet_body(cx, gas.state().state(), 300.0);
            let radius = live.mechanical().specimen().radius_m();
            let layers = MaterialConductionLayers {
                material: &card,
                context: &point,
                temperature_axis: "T",
                liquid_fraction_axis: None,
                conductivity: &key,
                selection: MaterialPropertySelection::SingleClaimOnly,
                area_m2: 2.0
                    * core::f64::consts::PI
                    * radius
                    * (live.mechanical().specimen().string().length_m + radius),
                gas_distance_m: 0.001,
                body_distance_m: radius / 2.0,
            };
            let exchange = HeatExchangeSolve {
                heat_bracket_j: [-8e-5, 8e-5],
                source_tolerance_j: 1e-13,
                energy_tolerance_j: 2e-10,
                maximum_iterations: 64,
            };
            let mut reference = live.clone();
            let mut reference_gas = gas.clone();
            let expected = reference
                .step_with_reacting_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    &mut reference_gas,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    controls(),
                    exchange,
                    &layers,
                )
                .unwrap();
            let result = live
                .step_with_reacting_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    &mut gas,
                    &mut reservoir,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    controls(),
                    exchange,
                    &layers,
                    GasNozzle {
                        area_m2: area,
                        ..nozzle()
                    },
                    flow_solve(),
                )
                .unwrap();
            assert_eq!(result.thermal.chemistry, expected.chemistry);
            assert_eq!(
                result.thermal.coupled.source.card_identity(),
                expected.coupled.source.card_identity()
            );
            assert_eq!(
                result.thermal.coupled.source.query_point(),
                expected.coupled.source.query_point()
            );
            assert_eq!(
                result.thermal.coupled.solid_to_gas_heat_j,
                expected.coupled.solid_to_gas_heat_j
            );
            assert_eq!(
                result.thermal.coupled.conductance_w_k,
                expected.coupled.conductance_w_k
            );
            assert_eq!(live.thermal(), reference.thermal());
            assert!(live.thermal().temperature_k() > 300.0);
            assert_eq!(live.mechanical().ambient(), gas.state().state());
            assert_eq!(live.mechanical().accepted_samples(), 1);
            assert!(result.energy_balance_residual_j.abs() <= result.energy_tolerance_j);
            if area == 0.0 {
                assert_eq!(gas, reference_gas);
                assert_eq!(
                    result.thermal.coupled.coupled.string,
                    expected.coupled.coupled.string
                );
            } else {
                assert_ne!(gas, reference_gas);
                assert!(result.transfer.transferred_mass_kg > 0.0);
            }
        }
    });
}

#[test]
fn g4_inlet_outside_reaction_catalog_refuses_after_transfer_without_publishing_or_losing_retry() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let (mut live, mut gas) = prepared(cx, &reaction);
        let before = live.clone();
        let old_gas = gas.clone();
        let (species, transport) = source::component("C", 0.028, 0.0);
        let foreign = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![species],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            )
            .unwrap(),
            vec![transport],
            None,
        )
        .unwrap();
        let mut reservoir = reservoir(foreign, 200000.0);
        let old_reservoir = reservoir.clone();
        // The original nozzle admits this real mixture; only the reacting
        // continuation cannot account for C. This exercises the late boundary.
        let mut reference = live.clone();
        let mut reference_gas = gas.clone();
        reference
            .step_with_reacting_isochoric_gas(
                cx,
                8,
                &[0.0; 2],
                &mut reference_gas,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
            )
            .unwrap();
        let (transfer, _) = reservoir
            .through_nozzle_to(
                &ClosedGasParcel::from_isochoric(reference_gas),
                nozzle(),
                DT,
                flow_solve(),
                || false,
            )
            .unwrap();
        assert!(transfer.into_parts().2.transferred_mass_kg > 0.0);
        let error = live
            .step_with_reacting_isochoric_gas_nozzle(
                cx,
                8,
                &[0.0; 2],
                &mut gas,
                &mut reservoir,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
                nozzle(),
                flow_solve(),
            )
            .unwrap_err();
        assert!(
            error.to_string().contains("outside the reaction source catalog"),
            "{error}"
        );
        assert_eq!(gas, old_gas);
        assert_eq!(reservoir, old_reservoir);
        assert_eq!(live.thermal(), before.thermal());
        assert_eq!(live.mechanical().states(), before.mechanical().states());
        assert_eq!(live.mechanical().ambient(), before.mechanical().ambient());
        assert_eq!(live.mechanical().epoch(), 8);
        assert_eq!(live.mechanical().accepted_samples(), 8);

        let mut reservoir = reactant_reservoir(&reaction, 200000.0);
        let mut reference_reservoir = reservoir.clone();
        let mut reference = before;
        let mut reference_gas = old_gas;
        let retry = live
            .step_with_reacting_isochoric_gas_nozzle(
                cx,
                8,
                &[0.0; 2],
                &mut gas,
                &mut reservoir,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
                nozzle(),
                flow_solve(),
            )
            .unwrap();
        let expected = reference
            .step_with_reacting_isochoric_gas_nozzle(
                cx,
                8,
                &[0.0; 2],
                &mut reference_gas,
                &mut reference_reservoir,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
                nozzle(),
                flow_solve(),
            )
            .unwrap();
        assert_eq!(retry.thermal.chemistry, expected.thermal.chemistry);
        assert_eq!(retry.thermal.coupled.string, expected.thermal.coupled.string);
        assert_eq!(retry.transfer, expected.transfer);
        assert_eq!(gas, reference_gas);
        assert_eq!(reservoir, reference_reservoir);
    });
}
