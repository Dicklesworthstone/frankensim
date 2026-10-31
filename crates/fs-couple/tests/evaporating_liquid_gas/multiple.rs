use std::cell::Cell;

use fs_couple::liquid_gas::evaporation::multiple::{
    DistributedLiquidGasCandidate, DistributedLiquidGasSolve, advance_liquid_surfaces_with_gas,
};
use fs_material::evaporation::gas::chamber::GasLiquidChamberSolve;

use super::*;

#[path = "multiple/heat.rs"]
mod heat;

const COLD: [usize; 3] = [3, 2, 5];
const WARM: [usize; 3] = [3, 4, 2];
const INITIAL: [([usize; 3], f64, f64); 2] = [(COLD, 0.4, 295.0), (WARM, 0.4, 305.0)];

fn joint_controls() -> DistributedLiquidGasSolve {
    let single = controls();
    DistributedLiquidGasSolve {
        chamber: GasLiquidChamberSolve {
            surface: single.phase,
            maximum_surfaces: 2,
        },
        native: single.native,
        mass_tolerance_kg: single.mass_tolerance_kg,
        energy_tolerance_j: single.energy_tolerance_j,
    }
}

fn patch<'a>(
    coordinates: [usize; 3],
    law: &'a dyn EvaporationLaw,
    binding: &'a GasLiquidBinding,
) -> LiquidGasEvaporationSurface<'a> {
    LiquidGasEvaporationSurface {
        liquid: coordinates,
        gas: [coordinates[0] + 1, coordinates[1], coordinates[2]],
        ..surface(law, binding)
    }
}

fn complete_inventory(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &ClosedGasParcel,
    binding: &GasLiquidBinding,
    next: &DistributedLiquidGasCandidate<'_, '_, '_>,
) {
    let liquid_mass_change =
        next.liquid().flow().ledger_mass_kg().unwrap() - liquid.flow().ledger_mass_kg().unwrap();
    let liquid_h_change =
        next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
    let gas_u_change = next.gas().mass_kg() * next.gas().state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let report = next.report();
    near(
        liquid_mass_change + next.gas().mass_kg() - gas.mass_kg(),
        0.0,
        joint_controls().mass_tolerance_kg,
    );
    near(
        liquid_mass_change + next.gas().species_mass_kg(&water()).unwrap()
            - gas.species_mass_kg(&water()).unwrap(),
        0.0,
        joint_controls().mass_tolerance_kg,
    );
    near(
        next.gas().species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-22,
    );
    near(
        liquid_h_change + binding.liquid_enthalpy_offset_j_kg * liquid_mass_change + gas_u_change,
        report.flow.external_heat_j + report.flow.viscous_heat_j + report.flow.body_heat_j,
        joint_controls().energy_tolerance_j,
    );
    near(
        liquid_h_change,
        report.native_liquid_enthalpy_change_j,
        1e-12,
    );
    assert!(report.absolute_mass_defect_kg <= joint_controls().mass_tolerance_kg);
    assert!(report.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
}

#[test]
fn g1_joint_native_cold_condensation_and_warm_evaporation_share_one_gas_and_clock() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let binding = binding(&source, 0.0);
        let surfaces = [patch(COLD, &Phase, &binding), patch(WARM, &Phase, &binding)];
        let next = advance_liquid_surfaces_with_gas(
            liquid,
            &gas,
            &source,
            &surfaces,
            joint_controls(),
            gate,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.surfaces.len(), 2);
        assert_eq!(report.surfaces[0].coordinates, COLD);
        assert_eq!(report.surfaces[1].coordinates, WARM);
        assert!(report.surfaces[0].phase.transferred_mass_kg < 0.0);
        assert!(report.surfaces[1].phase.transferred_mass_kg > 0.0);
        assert!(report.surfaces[0].native.heat_j > 0.0);
        assert!(report.surfaces[1].native.heat_j < 0.0);
        assert!(report.exchange.evaporated_mass_kg > 0.0);
        assert!(report.exchange.condensed_mass_kg > 0.0);
        for (index, row) in report.surfaces.iter().enumerate() {
            assert_eq!(row.phase.id, index as u64);
            assert_eq!(row.coordinates, row.native.coordinates);
            assert_eq!(row.native.mass_change_kg, -row.phase.transferred_mass_kg);
            assert_eq!(
                row.phase.vapor_partial_pressure_pa[0],
                gas.species_partial_pressure_pa(&water()).unwrap(),
            );
            near(row.phase.surface_temperature_k[0], INITIAL[index].2, 1e-11);
        }
        assert_eq!(report.mass_exchange.physical_time_s, 0.0);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        complete_inventory(liquid, &gas, &binding, &next);

        let reversed = advance_liquid_surfaces_with_gas(
            liquid,
            &gas,
            &source,
            &[surfaces[1], surfaces[0]],
            joint_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(next.gas(), reversed.gas());
        assert_eq!(report.exchange, reversed.report().exchange);
        assert_eq!(report.mass_exchange, reversed.report().mass_exchange);
        assert_eq!(report.flow, reversed.report().flow);
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            reversed.liquid().flow().lattice().grid().state_bits(),
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            reversed.liquid().total_enthalpy_j().unwrap(),
        );
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn g1_joint_native_continuation_retains_actual_shared_gas_and_advected_liquid() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let binding = binding(&source, 0.0);
        let surfaces = [patch(COLD, &Phase, &binding), patch(WARM, &Phase, &binding)];
        let first = advance_liquid_surfaces_with_gas(
            liquid,
            &gas,
            &source,
            &surfaces,
            joint_controls(),
            gate,
        )
        .unwrap();
        let second = advance_liquid_surfaces_with_gas(
            first.liquid(),
            first.gas(),
            first.catalog(),
            &surfaces,
            joint_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(second.liquid().accepted_steps(), 2);
        assert_eq!(second.liquid().physical_time_s(), 2.0 * DT);
        assert_eq!(second.report().mass_exchange.physical_time_s, DT);
        for row in &second.report().surfaces {
            assert_eq!(
                row.phase.vapor_partial_pressure_pa[0],
                first.gas().species_partial_pressure_pa(&water()).unwrap(),
            );
            assert_eq!(
                row.native.specific_enthalpy_j_kg[0],
                first
                    .liquid()
                    .specific_enthalpy_at(row.coordinates)
                    .unwrap(),
            );
        }
        complete_inventory(first.liquid(), first.gas(), &binding, &second);
        assert_eq!(first.liquid().physical_time_s(), DT);
        assert_eq!(gas.point().axes().get("specimen-lot"), Some(&27.0));
        assert_eq!(second.gas().point().axes().get("specimen-lot"), Some(&27.0));
        assert_eq!(second.gas().volume_m3(), gas.volume_m3());
    });
}

#[test]
fn g3_joint_native_singleton_preserves_existing_phase_only_endpoints() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let single =
            advance_liquid_with_gas(liquid, &gas, &source, face, controls(), gate).unwrap();
        let joint = advance_liquid_surfaces_with_gas(
            liquid,
            &gas,
            &source,
            &[face],
            joint_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(single.gas(), joint.gas());
        assert_eq!(single.report().mass_exchange, joint.report().mass_exchange);
        assert_eq!(single.report().flow, joint.report().flow);
        assert_eq!(
            single.liquid().flow().lattice().grid().state_bits(),
            joint.liquid().flow().lattice().grid().state_bits(),
        );
        assert_eq!(
            single.liquid().total_enthalpy_j().unwrap(),
            joint.liquid().total_enthalpy_j().unwrap(),
        );
        complete_inventory(liquid, &gas, &binding, &joint);
    });
}

#[test]
fn g4_joint_native_duplicate_work_budget_late_source_and_cancel_publish_neither_owner() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let binding = binding(&source, 0.0);
        let faces = [patch(COLD, &Phase, &binding), patch(WARM, &Phase, &binding)];
        let original_gas = gas.clone();
        let original_bits = liquid.flow().lattice().grid().state_bits();
        let original_h = liquid.total_enthalpy_j().unwrap();
        assert!(
            advance_liquid_surfaces_with_gas(
                liquid,
                &gas,
                &source,
                &[faces[0], faces[0]],
                joint_controls(),
                gate,
            )
            .is_err()
        );
        let mut bounded = joint_controls();
        bounded.chamber.maximum_surfaces = 1;
        assert!(
            advance_liquid_surfaces_with_gas(liquid, &gas, &source, &faces, bounded, gate,)
                .is_err()
        );
        let mut other_binding = binding.clone();
        other_binding.liquid_enthalpy_offset_j_kg += 1.0;
        assert!(
            advance_liquid_surfaces_with_gas(
                liquid,
                &gas,
                &source,
                &[faces[0], patch(WARM, &Phase, &other_binding)],
                joint_controls(),
                gate,
            )
            .is_err()
        );

        struct CountedPhase<'a>(&'a Cell<usize>);
        impl EvaporationLaw for CountedPhase<'_> {
            fn temperature_range_k(&self) -> [f64; 2] {
                Phase.temperature_range_k()
            }
            fn properties(
                &self,
                t: f64,
            ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
            {
                self.0.set(self.0.get() + 1);
                Phase.properties(t)
            }
        }
        struct MissingPhase;
        impl EvaporationLaw for MissingPhase {
            fn temperature_range_k(&self) -> [f64; 2] {
                Phase.temperature_range_k()
            }
            fn properties(
                &self,
                _: f64,
            ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
            {
                Err(EvaporationError::Invalid(
                    "later surface source unavailable",
                ))
            }
        }
        let calls = Cell::new(0);
        let counted = CountedPhase(&calls);
        assert!(
            advance_liquid_surfaces_with_gas(
                liquid,
                &gas,
                &source,
                &[
                    patch(COLD, &counted, &binding),
                    patch(WARM, &MissingPhase, &binding)
                ],
                joint_controls(),
                gate,
            )
            .is_err()
        );
        assert!(
            calls.get() > 0,
            "first surface must have been prepared before later refusal"
        );
        gate.request();
        assert!(matches!(
            advance_liquid_surfaces_with_gas(liquid, &gas, &source, &faces, joint_controls(), gate,),
            Err(EvaporatingLiquidGasError::Cancelled),
        ));
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
