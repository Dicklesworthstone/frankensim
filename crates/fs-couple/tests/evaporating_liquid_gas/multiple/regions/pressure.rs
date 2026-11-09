//! Regional post-phase EOS pressures reach the original native gas-link owner.
use super::*;
use fs_couple::liquid_gas::evaporation::multiple::regions::pressure::{
    RegionalGasPressure, advance_liquid_regions_with_gas_and_heat_and_pressure,
    advance_liquid_regions_with_gas_and_pressure,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

const REFERENCE: f64 = 100000.0;
const MAXIMUM: f64 = 0.05;

fn gas_regions() -> BTreeMap<[usize; 3], u64> {
    let mut result = BTreeMap::new();
    // Cover every original D3Q19 gas neighbor, including the diagonal links
    // and gas cells without an authored evaporation face.
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if x > 3 || y == 0 || y == 7 || z == 0 || z == 7 {
                    result.insert([x, y, z], if z >= 4 { DRY } else { HUMID });
                }
            }
        }
    }
    result
}

fn atmospheres(source: &ThermochemicalGasModel) -> BTreeMap<u64, ClosedGasParcel> {
    [(DRY, 320.0, 1000.0), (HUMID, 290.0, 6000.0)]
        .map(|(id, t, vapor)| {
            let initial = gas(source, t, vapor);
            let parcel = ClosedGasParcel::new(
                initial.model().clone(),
                initial.point().clone(),
                3e-7,
            )
            .unwrap();
            (id, parcel)
        })
        .into_iter()
        .collect()
}

fn pressure_from_inventory(gas: &ClosedGasParcel) -> f64 {
    let amount = gas.species_mass_kg(&carrier()).unwrap() / 0.028
        + gas.species_mass_kg(&water()).unwrap() / M;
    amount * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * gas.state().state().temperature
        / gas.volume_m3()
}

fn normal_velocity(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>, cell: [usize; 3]) -> f64 {
    liquid
        .flow()
        .sample_at_si(cell[0], cell[1], cell[2])
        .unwrap()
        .unwrap()
        .velocity_m_s[0]
}

#[test]
fn distinct_post_phase_and_heat_pressures_reach_native_links_without_extra_energy_or_steps() {
    for with_heat in [false, true] {
        with_liquid_interfaces_pressure(0.0, &CELLS, Some(REFERENCE), |liquid, gate| {
            let source = source(-13.4e6);
            let original = atmospheres(&source);
            let preserved = original.clone();
            let old_bits = liquid.flow().lattice().grid().state_bits();
            let old_h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(&source, 0.0);
            let rows = surfaces(&binding);
            let heated = rows.map(|row| RegionalLiquidGasHeatSurface {
                region_id: row.region_id,
                surface: LiquidGasHeatSurface {
                    surface: row.surface,
                    heat_transfer_coefficient_w_m2_k: 50.0,
                    heat_controls: heat_controls(),
                },
            });
            let membership = gas_regions();
            let pressure = RegionalGasPressure {
                gas_regions: &membership,
                maximum_density_excursion: MAXIMUM,
            };
            let thermal = GasLiquidChamberHeatSolve { maximum_iterations: 128 };
            let next = if with_heat {
                advance_liquid_regions_with_gas_and_heat_and_pressure(
                    liquid, &original, &source, &heated, joint_controls(), thermal, pressure, gate,
                )
            } else {
                advance_liquid_regions_with_gas_and_pressure(
                    liquid, &original, &source, &rows, joint_controls(), pressure, gate,
                )
            }
            .unwrap();
            let unloaded = if with_heat {
                advance_liquid_regions_with_gas_and_heat(
                    liquid, &original, &source, &heated, joint_controls(), thermal, gate,
                )
            } else {
                advance_liquid_regions_with_gas(
                    liquid, &original, &source, &rows, joint_controls(), gate,
                )
            }
            .unwrap();
            let report = next.report();
            let receipt = report.pressure.as_ref().unwrap();
            assert!(unloaded.report().pressure.is_none());
            assert_eq!(next.parcels(), unloaded.parcels());
            assert_eq!(report.mass_exchange, unloaded.report().mass_exchange);
            assert_eq!(receipt.regions.len(), 2);
            for (id, cell) in [(DRY, COLD), (HUMID, WARM)] {
                let load = receipt.regions[&id];
                assert_eq!(load.initial_pressure_pa, original[&id].state().state().pressure);
                assert_eq!(load.applied_pressure_pa, next.parcels()[&id].state().state().pressure);
                assert_eq!(load.volume_m3, original[&id].volume_m3());
                assert_eq!(load.maximum_density_excursion, MAXIMUM);
                assert!(receipt.sampled_links[&id] > 5, "unpatched interfaces also load");
                assert_eq!(report.regions[&id].phase, unloaded.report().regions[&id].phase);
                assert_eq!(report.regions[&id].exchange, unloaded.report().regions[&id].exchange);
                assert_eq!(report.regions[&id].heat, unloaded.report().regions[&id].heat);
                assert_eq!(report.regions[&id].heat.is_some(), with_heat);
                region_accounts(&original[&id], &next.parcels()[&id], &report.regions[&id]);

                // Independent species inventory -> EOS -> flat-front moment.
                // All five missing links at this cell belong to this region.
                // Their weights sum to 1/6; rho=1000 and dx/dt=1. The fixture's
                // linear tangential shear cancels in the normal moment.
                let expected_pressure = pressure_from_inventory(&next.parcels()[&id]);
                near(load.applied_pressure_pa, expected_pressure, 2e-10);
                let d = (expected_pressure - REFERENCE) / 1000.0;
                near(normal_velocity(next.liquid(), cell), -d / (1.0 + d), 2e-11);
                near(normal_velocity(unloaded.liquid(), cell), 0.0, 2e-13);
            }
            assert_ne!(
                receipt.regions[&DRY].applied_pressure_pa,
                receipt.regions[&HUMID].applied_pressure_pa
            );
            assert_ne!(
                next.liquid().flow().lattice().grid().state_bits(),
                unloaded.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_eq!(next.liquid().physical_time_s(), DT);
            assert_eq!(report.mass_exchange.physical_time_s, 0.0);
            assert_eq!(report.flow.external_heat_j, 0.0, "pressure is not a second heat source");
            let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                - liquid.flow().ledger_mass_kg().unwrap();
            let mut gas_mass = 0.0;
            let mut gas_energy = 0.0;
            for (id, old) in &original {
                let new = &next.parcels()[id];
                gas_mass += new.mass_kg() - old.mass_kg();
                gas_energy += new.mass_kg() * new.state().specific_internal_energy_j_kg()
                    - old.mass_kg() * old.state().specific_internal_energy_j_kg();
                assert_eq!(new.volume_m3(), old.volume_m3());
            }
            near(dm + gas_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(
                next.liquid().total_enthalpy_j().unwrap() - old_h
                    + binding.liquid_enthalpy_offset_j_kg * dm + gas_energy,
                report.flow.viscous_heat_j + report.flow.body_heat_j,
                joint_controls().energy_tolerance_j,
            );
            assert!(report.absolute_mass_defect_kg <= joint_controls().mass_tolerance_kg);
            assert!(report.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
            assert_eq!(original, preserved);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
            assert_eq!(liquid.accepted_steps(), 0);
        });
    }
}

#[test]
fn incomplete_conflicting_and_excessive_pressure_loads_refuse_atomically_and_retry() {
    with_liquid_interfaces_pressure(0.0, &CELLS, Some(REFERENCE), |liquid, gate| {
        let source = source(-13.4e6);
        let original = atmospheres(&source);
        let preserved = original.clone();
        let old_bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        let binding = binding(&source, 0.0);
        let rows = surfaces(&binding);
        let complete = gas_regions();
        let run = |map: &BTreeMap<[usize; 3], u64>, maximum| {
            advance_liquid_regions_with_gas_and_pressure(
                liquid,
                &original,
                &source,
                &rows,
                joint_controls(),
                RegionalGasPressure {
                    gas_regions: map,
                    maximum_density_excursion: maximum,
                },
                gate,
            )
        };
        let mut missing_link = complete.clone();
        // This diagonal is not an authored phase face, so admission proceeds
        // until the original native pressure callback actually samples it.
        missing_link.remove(&[4, 3, 5]);
        assert!(matches!(run(&missing_link, MAXIMUM), Err(EvaporatingLiquidGasError::Invalid(_))));
        let mut contradictory = complete.clone();
        contradictory.insert([4, 2, 5], HUMID);
        assert!(matches!(run(&contradictory, MAXIMUM), Err(EvaporatingLiquidGasError::Invalid(_))));
        let mut unknown = complete.clone();
        unknown.insert([4, 2, 5], 999);
        assert!(matches!(run(&unknown, MAXIMUM), Err(EvaporatingLiquidGasError::Invalid(_))));
        let mut absent = original.clone();
        absent.remove(&HUMID);
        assert!(advance_liquid_regions_with_gas_and_pressure(
            liquid,
            &absent,
            &source,
            &rows,
            joint_controls(),
            RegionalGasPressure {
                gas_regions: &complete,
                maximum_density_excursion: MAXIMUM,
            },
            gate,
        ).is_err());
        for maximum in [0.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(run(&complete, maximum), Err(EvaporatingLiquidGasError::Invalid(_))));
        }
        // A finite positive limit passes admission, then rejects the actual
        // post-phase EOS load after source/native mass and h proposals exist.
        assert!(matches!(run(&complete, 1e-12), Err(EvaporatingLiquidGasError::Liquid(_))));
        assert_eq!(original, preserved);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
        let retry = run(&complete, MAXIMUM).unwrap();
        assert_eq!(retry.liquid().accepted_steps(), 1);
        assert_eq!(retry.liquid().physical_time_s(), DT);
        gate.request();
        assert!(matches!(run(&complete, MAXIMUM), Err(EvaporatingLiquidGasError::Cancelled)));
        assert_eq!(original, preserved);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
    });
    with_liquid_interfaces(0.0, &CELLS, |liquid, gate| {
        let source = source(-13.4e6);
        let original = atmospheres(&source);
        let binding = binding(&source, 0.0);
        let membership = gas_regions();
        assert!(matches!(
            advance_liquid_regions_with_gas_and_pressure(
                liquid,
                &original,
                &source,
                &surfaces(&binding),
                joint_controls(),
                RegionalGasPressure {
                    gas_regions: &membership,
                    maximum_density_excursion: MAXIMUM,
                },
                gate,
            ),
            Err(EvaporatingLiquidGasError::Invalid(_))
        ));
        assert_eq!(liquid.accepted_steps(), 0);
    });
}
