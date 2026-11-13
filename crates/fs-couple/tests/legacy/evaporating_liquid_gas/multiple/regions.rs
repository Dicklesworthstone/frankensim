//! Distinct finite gas regions retain their own humidity and caloric accounts.
use super::*;

#[path = "regions/pressure.rs"]
mod pressure;
use fs_couple::liquid_gas::evaporation::multiple::{
    LiquidGasHeatSurface,
    regions::{
        RegionalLiquidGasHeatSurface, RegionalLiquidGasRegionReport, RegionalLiquidGasSurface,
        advance_liquid_regions_with_gas, advance_liquid_regions_with_gas_and_heat,
    },
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;
use std::collections::BTreeMap;

const DRY: u64 = 10;
const HUMID: u64 = 20;
const IDLE: u64 = 90;
const CELLS: [([usize; 3], f64, f64); 2] = [(COLD, 0.4, 300.0), (WARM, 0.4, 300.0)];

fn parcels(source: &ThermochemicalGasModel) -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([
        (DRY, gas(source, 300.0, 1000.0)),
        (HUMID, gas(source, 310.0, 6000.0)),
        (IDLE, gas(source, 290.0, 2500.0)),
    ])
}
fn surfaces(binding: &GasLiquidBinding) -> [RegionalLiquidGasSurface<'_>; 2] {
    [
        RegionalLiquidGasSurface {
            region_id: DRY,
            surface: patch(COLD, &Phase, binding),
        },
        RegionalLiquidGasSurface {
            region_id: HUMID,
            surface: patch(WARM, &Phase, binding),
        },
    ]
}
fn region_accounts(
    original: &ClosedGasParcel,
    next: &ClosedGasParcel,
    report: &RegionalLiquidGasRegionReport,
) {
    let transferred: f64 = report
        .surfaces
        .iter()
        .map(|row| row.phase.transferred_mass_kg)
        .sum();
    near(
        next.species_mass_kg(&water()).unwrap() - original.species_mass_kg(&water()).unwrap(),
        transferred,
        controls().phase.gas.mass_tolerance_kg,
    );
    near(
        next.species_mass_kg(&carrier()).unwrap(),
        original.species_mass_kg(&carrier()).unwrap(),
        1e-22,
    );
    let stored = next.mass_kg() * next.state().specific_internal_energy_j_kg()
        - original.mass_kg() * original.state().specific_internal_energy_j_kg();
    near(
        stored,
        report.exchange.gas_internal_energy_change_j,
        controls().phase.gas.energy_tolerance_j,
    );
    near(
        report.exchange.liquid_enthalpy_change_j + stored,
        0.0,
        controls().phase.gas.energy_tolerance_j,
    );
    near(
        report.exchange.energy_residual_j,
        0.0,
        controls().phase.gas.energy_tolerance_j,
    );
    assert_eq!(next.volume_m3(), original.volume_m3());
    assert_eq!(next.point().axes().get("specimen-lot"), Some(&27.0));
    assert_eq!(report.exchange.duration_s, DT);
    for row in &report.surfaces {
        assert_eq!(row.native.mass_change_kg, -row.phase.transferred_mass_kg);
    }
}

#[test]
fn separate_humidity_regions_evaporate_and_condense_with_one_native_step_and_canonical_order() {
    with_liquid_interfaces(0.0, &CELLS, |liquid, gate| {
        let source = source(-13.4e6);
        let original = parcels(&source);
        let saved = original.clone();
        let binding = binding(&source, 0.0);
        let rows = surfaces(&binding);
        let next = advance_liquid_regions_with_gas(
            liquid,
            &original,
            &source,
            &rows,
            joint_controls(),
            gate,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(
            report.regions.keys().copied().collect::<Vec<_>>(),
            vec![DRY, HUMID]
        );
        assert_eq!(report.mass_exchange.samples.len(), 2);
        assert_eq!(report.mass_exchange.physical_time_s, 0.0);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert_eq!(next.parcels()[&IDLE], original[&IDLE]);
        assert!(std::ptr::eq(next.catalog(), &source));
        for (id, coordinates, evaporation) in [(DRY, COLD, true), (HUMID, WARM, false)] {
            let regional = &report.regions[&id];
            assert_eq!(regional.surfaces.len(), 1);
            let row = &regional.surfaces[0];
            assert_eq!(row.coordinates, coordinates);
            assert_eq!(row.phase.transferred_mass_kg > 0.0, evaporation);
            assert_eq!(
                row.phase.vapor_partial_pressure_pa[0],
                original[&id].species_partial_pressure_pa(&water()).unwrap()
            );
            near(row.phase.surface_temperature_k[0], 300.0, 1e-11);
            region_accounts(&original[&id], &next.parcels()[&id], regional);
        }
        let liquid_mass = next.liquid().flow().ledger_mass_kg().unwrap()
            - liquid.flow().ledger_mass_kg().unwrap();
        let gas_mass: f64 = original
            .iter()
            .map(|(id, old)| next.parcels()[id].mass_kg() - old.mass_kg())
            .sum();
        near(
            liquid_mass + gas_mass,
            0.0,
            joint_controls().mass_tolerance_kg,
        );
        let gas_energy: f64 = original
            .iter()
            .map(|(id, old)| {
                next.parcels()[id].mass_kg()
                    * next.parcels()[id].state().specific_internal_energy_j_kg()
                    - old.mass_kg() * old.state().specific_internal_energy_j_kg()
            })
            .sum();
        let liquid_energy = next.liquid().total_enthalpy_j().unwrap()
            - liquid.total_enthalpy_j().unwrap()
            + binding.liquid_enthalpy_offset_j_kg * liquid_mass;
        near(
            liquid_energy + gas_energy,
            report.flow.external_heat_j + report.flow.viscous_heat_j + report.flow.body_heat_j,
            joint_controls().energy_tolerance_j,
        );
        assert!(report.absolute_mass_defect_kg <= joint_controls().mass_tolerance_kg);
        assert!(report.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);

        let reversed = advance_liquid_regions_with_gas(
            liquid,
            &original,
            &source,
            &[rows[1], rows[0]],
            joint_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(next.parcels(), reversed.parcels());
        assert_eq!(report.mass_exchange, reversed.report().mass_exchange);
        assert_eq!(report.flow, reversed.report().flow);
        for id in [DRY, HUMID] {
            assert_eq!(
                report.regions[&id].phase,
                reversed.report().regions[&id].phase
            );
            assert_eq!(
                report.regions[&id].exchange,
                reversed.report().regions[&id].exchange
            );
        }
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            reversed.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            reversed.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(original, saved);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn each_region_supplies_its_own_reciprocal_heat_before_its_phase_exchange() {
    with_liquid_interfaces(0.0, &CELLS, |liquid, gate| {
        let source = source(-13.4e6);
        let original = BTreeMap::from([
            (DRY, gas(&source, 320.0, 1000.0)),
            (HUMID, gas(&source, 290.0, 6000.0)),
        ]);
        let binding = binding(&source, 0.0);
        let rows = surfaces(&binding).map(|row| RegionalLiquidGasHeatSurface {
            region_id: row.region_id,
            surface: LiquidGasHeatSurface {
                surface: row.surface,
                heat_transfer_coefficient_w_m2_k: 50.0,
                heat_controls: heat_controls(),
            },
        });
        let next = advance_liquid_regions_with_gas_and_heat(
            liquid,
            &original,
            &source,
            &rows,
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        for (id, warming) in [(DRY, true), (HUMID, false)] {
            let report = &next.report().regions[&id];
            let heat = report.heat.as_ref().unwrap();
            let row = &report.surfaces[0];
            let surface_heat = row.heat.as_ref().unwrap();
            assert_eq!(heat.heat_into_liquid_j > 0.0, warming);
            assert_eq!(row.phase.surface_temperature_k[0] > 300.0, warming);
            assert_eq!(row.phase.transferred_mass_kg > 0.0, warming);
            assert_eq!(
                surface_heat.temperature_k[1],
                row.phase.surface_temperature_k[0]
            );
            assert_eq!(heat.gas_temperature_k[1], report.phase.gas_temperature_k[0]);
            near(
                heat.heat_into_liquid_j,
                50.0 * 1e-6
                    * DT
                    * (report.phase.gas_temperature_k[0] - row.phase.surface_temperature_k[0]),
                1e-12,
            );
            near(
                report.exchange.gas_internal_energy_change_j,
                report.exchange.energy_into_gas_j - heat.heat_into_liquid_j,
                controls().phase.gas.energy_tolerance_j,
            );
            region_accounts(&original[&id], &next.parcels()[&id], report);
        }
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert_eq!(next.report().mass_exchange.samples.len(), 2);
        assert!(next.report().absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
    });
}

#[test]
fn missing_repeated_and_incompatible_region_bindings_refuse_before_native_publication() {
    let corner = [3, 6, 6];
    with_liquid_interfaces(
        0.0,
        &[(corner, 0.4, 300.0), (COLD, 0.4, 300.0)],
        |liquid, gate| {
            let source = source(-13.4e6);
            let original = parcels(&source);
            let saved = original.clone();
            let binding = binding(&source, 0.0);
            let baseline = [
                RegionalLiquidGasSurface {
                    region_id: DRY,
                    surface: patch(corner, &Phase, &binding),
                },
                RegionalLiquidGasSurface {
                    region_id: HUMID,
                    surface: patch(COLD, &Phase, &binding),
                },
            ];
            let bits = liquid.flow().lattice().grid().state_bits();
            let old_h = liquid.total_enthalpy_j().unwrap();
            let mut incompatible = binding.clone();
            incompatible.liquid_enthalpy_offset_j_kg += 1.0;
            for change in 0..4 {
                let mut rows = baseline;
                match change {
                    0 => rows[1].region_id = 999,
                    1 => rows[1] = rows[0],
                    2 => {
                        // Both faces are actual gas neighbors, but one native
                        // liquid inventory cannot be assigned two regional owners.
                        rows[1].surface = rows[0].surface;
                        rows[1].surface.gas = [3, 7, 6];
                    }
                    _ => rows[1].surface.binding = &incompatible,
                }
                assert!(
                    advance_liquid_regions_with_gas(
                        liquid,
                        &original,
                        &source,
                        &rows,
                        joint_controls(),
                        gate,
                    )
                    .is_err()
                );
            }
            assert_eq!(original, saved);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(liquid.physical_time_s(), 0.0);
            let retry = advance_liquid_regions_with_gas(
                liquid,
                &original,
                &source,
                &baseline,
                joint_controls(),
                gate,
            )
            .unwrap();
            assert_eq!(retry.liquid().accepted_steps(), 1);
        },
    );
}

#[test]
fn a_late_second_region_source_refusal_discards_both_gas_updates_and_retries_once() {
    struct EndpointRefusal<'a>(&'a Cell<usize>);
    impl EvaporationLaw for EndpointRefusal<'_> {
        fn temperature_range_k(&self) -> [f64; 2] {
            Phase.temperature_range_k()
        }
        fn properties(
            &self,
            t: f64,
        ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError> {
            if t > 300.0 + 1e-9 {
                self.0.set(self.0.get() + 1);
                return Err(EvaporationError::Invalid(
                    "regional condensation endpoint unavailable",
                ));
            }
            Phase.properties(t)
        }
    }
    with_liquid_interfaces(0.0, &CELLS, |liquid, gate| {
        let source = source(-13.4e6);
        let original = parcels(&source);
        let saved = original.clone();
        let binding = binding(&source, 0.0);
        let baseline = surfaces(&binding);
        let refusals = Cell::new(0);
        let failing = EndpointRefusal(&refusals);
        let mut rows = baseline;
        rows[1].surface.law = &failing;
        let bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        assert!(
            advance_liquid_regions_with_gas(
                liquid,
                &original,
                &source,
                &rows,
                joint_controls(),
                gate,
            )
            .is_err()
        );
        assert!(
            refusals.get() > 0,
            "refusal must reach the actual condensation endpoint"
        );
        assert_eq!(original, saved);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
        let retry = advance_liquid_regions_with_gas(
            liquid,
            &original,
            &source,
            &baseline,
            joint_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(retry.liquid().accepted_steps(), 1);
        assert_eq!(retry.liquid().physical_time_s(), DT);
        assert!(retry.report().regions[&DRY].exchange.transferred_mass_kg > 0.0);
        assert!(retry.report().regions[&HUMID].exchange.transferred_mass_kg < 0.0);
    });
}
