//! Manufactured chemistry drives phase transfer on the original native clock.

#[path = "reactive/conduction.rs"]
mod conduction;
#[path = "reactive/momentum.rs"]
mod momentum;
#[path = "reactive/pressure.rs"]
mod pressure;

use super::*;
use fs_couple::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use fs_couple::liquid_gas::evaporation::multiple::{
    DistributedLiquidGasSolve, LiquidGasHeatSurface, advance_liquid_surfaces_with_gas,
    advance_liquid_surfaces_with_gas_and_heat,
};
use fs_couple::liquid_gas::evaporation::reactive::{
    ReactiveLiquidPhaseCandidate, ReactiveLiquidPhaseError, ReactiveLiquidPhaseSolve,
};
use fs_material::evaporation::gas::chamber::{
    GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve,
};
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::*;

#[allow(dead_code)]
#[path = "../../material_exchange/gas.rs"]
mod gas_data;

// Equal-mass manufactured isomers. The volatile uses the existing authored
// phase fixture; this is neither an empirical water mechanism nor combustion.
pub(super) fn reaction(rate: f64, release_j_mol: f64) -> GasReaction {
    let reference = -13.4e6 * M;
    let (a, ta) = gas_data::component("carrier", M, reference + release_j_mol);
    let (b, tb) = gas_data::component("water", M, reference);
    let species = vec![a.species().clone(), b.species().clone()];
    let id = ReactionId::new("manufactured-volatile-production").unwrap();
    let matrix =
        StoichiometricMatrix::new(species.clone(), vec![id.clone()], vec![vec![-1], vec![1]])
            .unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let certificate = verify_conservation(&elements, &matrix, &charge).unwrap();
    let equilibrium =
        IdealGasReactionEquilibriumV1::new(matrix, certificate, id, vec![a.clone(), b.clone()])
            .unwrap();
    let kinetics = DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(rate)).unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    );
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![ta, tb],
        Some(&water()),
    )
    .unwrap();
    GasReaction::new(catalog, kinetics).unwrap()
}

pub(super) fn reactive_controls() -> ReactiveLiquidPhaseSolve {
    ReactiveLiquidPhaseSolve {
        reaction: GasReactionSolve {
            gas: controls().phase.gas.gas,
            maximum_consumed_fraction: 0.2,
            species_tolerance_mol: 1e-15,
            entropy_tolerance_j_k: 1e-12,
        },
        maximum_gas_temperature_increment_k: 10.0,
        mass_tolerance_kg: 1e-15,
        energy_tolerance_j: 1e-9,
    }
}

fn complete_inventory<R>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &ClosedGasParcel,
    binding: &GasLiquidBinding,
    next: &ReactiveLiquidPhaseCandidate<'_, '_, '_, R>,
    native_supplied_j: f64,
) {
    let dm =
        next.liquid().flow().ledger_mass_kg().unwrap() - liquid.flow().ledger_mass_kg().unwrap();
    let dh = next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
    let du = next.gas().mass_kg() * next.gas().state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    near(
        dm + next.gas().mass_kg() - gas.mass_kg(),
        0.0,
        reactive_controls().mass_tolerance_kg,
    );
    near(
        dh + binding.liquid_enthalpy_offset_j_kg * dm + du,
        native_supplied_j,
        reactive_controls().energy_tolerance_j,
    );
    let chemical_volatile = next
        .report()
        .reaction
        .species
        .iter()
        .find(|row| row.species == water())
        .unwrap();
    near(
        dm + next.gas().species_mass_kg(&water()).unwrap()
            - gas.species_mass_kg(&water()).unwrap()
            - M * chemical_volatile.stoichiometric_change_mol,
        0.0,
        reactive_controls().mass_tolerance_kg,
    );
    assert!(next.report().absolute_mass_defect_kg <= reactive_controls().mass_tolerance_kg);
    assert!(next.report().absolute_energy_defect_j <= reactive_controls().energy_tolerance_j);
    assert_eq!(next.liquid().accepted_steps(), liquid.accepted_steps() + 1);
    assert_eq!(next.gas().volume_m3(), gas.volume_m3());
}

#[test]
fn g1_chemical_vapor_production_reverses_native_phase_flux_without_extra_heat_or_clock() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let reaction = reaction(2000.0, 75.0);
        let source = reaction.catalog();
        let gas = gas(source, 300.0, 1000.0);
        let binding = binding(source, 0.0);
        let patch = surface(&Phase, &binding);
        let initial_bits = liquid.flow().lattice().grid().state_bits();
        let original_gas = gas.clone();
        let frozen =
            advance_liquid_with_gas(liquid, &gas, source, patch, controls(), gate).unwrap();
        assert!(frozen.report().phase.transferred_mass_kg > 0.0);
        let chemical = reaction
            .parcel_candidate(&gas, DT, reactive_controls().reaction, || false)
            .unwrap();
        assert!(chemical.gas().state().state().temperature > gas.state().state().temperature);
        assert!(
            chemical
                .gas()
                .species_partial_pressure_pa(&water())
                .unwrap()
                > 3536.8
        );
        let direct =
            advance_liquid_with_gas(liquid, chemical.gas(), source, patch, controls(), gate)
                .unwrap();
        let next = reaction
            .evaporating_liquid_candidate(
                liquid,
                &gas,
                patch,
                controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        let report = next.report();
        assert!(report.reaction.extent_mol > 0.0);
        assert!(report.phase.phase.transferred_mass_kg < 0.0);
        assert!(report.phase.mass_exchange.heat_j > 0.0);
        assert_eq!(&report.reaction, chemical.step());
        assert_eq!(next.gas(), direct.gas());
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            direct.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(report.duration_s, DT);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert_eq!(
            report.phase.exchange.gas_temperature_k[0],
            chemical.gas().state().state().temperature
        );
        assert!(core::ptr::eq(next.catalog(), source));
        let flow = report.phase.flow;
        complete_inventory(
            liquid,
            &gas,
            &binding,
            &next,
            flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
        );
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), initial_bits);
        assert_eq!(liquid.accepted_steps(), 0);
    });
}

#[test]
fn g3_reactive_sensible_heat_and_continuation_keep_original_sources_and_remaining_species() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let reaction = reaction(2000.0, 75.0);
        let source = reaction.catalog();
        let gas = gas(source, 305.0, 1000.0);
        let binding = binding(source, 0.0);
        let patch = surface(&Phase, &binding);
        let chemical = reaction
            .parcel_candidate(&gas, DT, reactive_controls().reaction, || false)
            .unwrap();
        let direct = advance_liquid_with_gas_and_heat(
            liquid,
            chemical.gas(),
            source,
            patch,
            500.0,
            controls(),
            heat_controls(),
            gate,
        )
        .unwrap();
        let next = reaction
            .evaporating_liquid_heat_candidate(
                liquid,
                &gas,
                patch,
                500.0,
                controls(),
                heat_controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        assert_eq!(next.gas(), direct.gas());
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert!(next.report().phase.heat.unwrap().heat_into_liquid_j > 0.0);
        let flow = next.report().phase.flow;
        complete_inventory(
            liquid,
            &gas,
            &binding,
            &next,
            flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
        );
        let second = reaction
            .evaporating_liquid_heat_candidate(
                next.liquid(),
                next.gas(),
                patch,
                500.0,
                controls(),
                heat_controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        assert_eq!(second.liquid().accepted_steps(), 2);
        assert_eq!(second.liquid().physical_time_s(), 2.0 * DT);
        assert!(second.report().reaction.extent_mol < next.report().reaction.extent_mol);
        let before = next.gas().species_mass_kg(&water()).unwrap() / M;
        let volatile = second
            .report()
            .reaction
            .species
            .iter()
            .find(|row| row.species == water())
            .unwrap();
        near(volatile.initial_mol, before, 1e-20);
        assert_eq!(second.gas().point().axes().get("specimen-lot"), Some(&27.0));
        let flow = second.report().phase.flow;
        complete_inventory(
            next.liquid(),
            next.gas(),
            &binding,
            &second,
            flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
        );
    });
}

#[test]
fn g3_zero_reaction_extent_and_zero_sensible_heat_recover_original_phase_owner() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let reaction = reaction(2000.0, 0.0);
        let source = reaction.catalog();
        let gas = gas(source, 300.0, 50000.0);
        let binding = binding(source, 0.0);
        let patch = surface(&Phase, &binding);
        let original =
            advance_liquid_with_gas(liquid, &gas, source, patch, controls(), gate).unwrap();
        let next = reaction
            .evaporating_liquid_candidate(
                liquid,
                &gas,
                patch,
                controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        let zero_heat = reaction
            .evaporating_liquid_heat_candidate(
                liquid,
                &gas,
                patch,
                0.0,
                controls(),
                heat_controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        assert_eq!(next.report().reaction.extent_mol, 0.0);
        assert_eq!(next.gas(), original.gas());
        assert_eq!(next.gas(), zero_heat.gas());
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            original.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            zero_heat.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            zero_heat.report().phase.heat.unwrap().heat_into_liquid_j,
            0.0
        );
    });
}

#[test]
fn g3_reacting_multiple_surfaces_share_one_gas_for_heat_and_phase_in_canonical_order() {
    let coordinates = [[3, 2, 5], [3, 4, 2]];
    with_liquid_interfaces(
        0.0,
        &[(coordinates[0], 0.4, 295.0), (coordinates[1], 0.4, 305.0)],
        |liquid, gate| {
            let reaction = reaction(500.0, 75.0);
            let source = reaction.catalog();
            let gas = gas(source, 300.0, 2500.0);
            let binding = binding(source, 0.0);
            let surfaces = coordinates.map(|liquid| LiquidGasEvaporationSurface {
                liquid,
                gas: [4, liquid[1], liquid[2]],
                ..surface(&Phase, &binding)
            });
            let single = controls();
            let joint = DistributedLiquidGasSolve {
                chamber: GasLiquidChamberSolve {
                    surface: single.phase,
                    maximum_surfaces: 2,
                },
                native: single.native,
                mass_tolerance_kg: single.mass_tolerance_kg,
                energy_tolerance_j: single.energy_tolerance_j,
            };
            let chemical = reaction
                .parcel_candidate(&gas, DT, reactive_controls().reaction, || false)
                .unwrap();
            let direct = advance_liquid_surfaces_with_gas(
                liquid,
                chemical.gas(),
                source,
                &surfaces,
                joint,
                gate,
            )
            .unwrap();
            let next = reaction
                .evaporating_liquid_surfaces_candidate(
                    liquid,
                    &gas,
                    &surfaces,
                    joint,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert!(next.report().phase.surfaces[0].phase.transferred_mass_kg < 0.0);
            assert!(next.report().phase.surfaces[1].phase.transferred_mass_kg > 0.0);
            let flow = next.report().phase.flow;
            complete_inventory(
                liquid,
                &gas,
                &binding,
                &next,
                flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
            );
            let heated = surfaces.map(|surface| LiquidGasHeatSurface {
                surface,
                heat_transfer_coefficient_w_m2_k: 1000.0,
                heat_controls: heat_controls(),
            });
            let joint_heat = GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            };
            let direct = advance_liquid_surfaces_with_gas_and_heat(
                liquid,
                chemical.gas(),
                source,
                &heated,
                joint,
                joint_heat,
                gate,
            )
            .unwrap();
            let next = reaction
                .evaporating_liquid_surfaces_heat_candidate(
                    liquid,
                    &gas,
                    &heated,
                    joint,
                    joint_heat,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            let permuted = reaction
                .evaporating_liquid_surfaces_heat_candidate(
                    liquid,
                    &gas,
                    &[heated[1], heated[0]],
                    joint,
                    joint_heat,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(next.gas(), permuted.gas());
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                permuted.liquid().flow().lattice().grid().state_bits()
            );
            let flow = next.report().phase.flow;
            complete_inventory(
                liquid,
                &gas,
                &binding,
                &next,
                flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
            );
        },
    );
}

#[test]
fn g4_reaction_source_temperature_and_late_phase_refusals_preserve_both_original_owners() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let reaction = reaction(2000.0, 75.0);
        let source = reaction.catalog();
        let gas = gas(source, 300.0, 1000.0);
        let binding = binding(source, 0.0);
        let original_gas = gas.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let mut narrow = reactive_controls();
        narrow.maximum_gas_temperature_increment_k = 1e-12;
        assert!(
            reaction
                .evaporating_liquid_candidate(
                    liquid,
                    &gas,
                    surface(&Phase, &binding),
                    controls(),
                    narrow,
                    gate,
                )
                .is_err()
        );
        let other = self::reaction(2000.0, 76.0);
        assert!(
            other
                .evaporating_liquid_candidate(
                    liquid,
                    &gas,
                    surface(&Phase, &binding),
                    controls(),
                    reactive_controls(),
                    gate,
                )
                .is_err()
        );

        struct Refuse;
        impl EvaporationLaw for Refuse {
            fn temperature_range_k(&self) -> [f64; 2] {
                [280.0, 340.0]
            }
            fn properties(
                &self,
                _: f64,
            ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
            {
                Err(EvaporationError::Invalid(
                    "original phase source refuses after chemical staging",
                ))
            }
        }
        assert!(
            reaction
                .evaporating_liquid_candidate(
                    liquid,
                    &gas,
                    surface(&Refuse, &binding),
                    controls(),
                    reactive_controls(),
                    gate,
                )
                .is_err()
        );

        struct Cancel<'a>(&'a CancelGate);
        impl EvaporationLaw for Cancel<'_> {
            fn temperature_range_k(&self) -> [f64; 2] {
                [280.0, 340.0]
            }
            fn properties(
                &self,
                t: f64,
            ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
            {
                self.0.request();
                Phase.properties(t)
            }
        }
        let error = reaction.evaporating_liquid_candidate(
            liquid,
            &gas,
            surface(&Cancel(gate), &binding),
            controls(),
            reactive_controls(),
            gate,
        );
        assert!(matches!(error, Err(ReactiveLiquidPhaseError::Cancelled)));
        assert!(gate.is_requested());
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
