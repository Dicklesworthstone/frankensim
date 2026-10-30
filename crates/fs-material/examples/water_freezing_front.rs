//! A computed ice front in a one-dimensional water layer, using original calorics.
//!
//! Run: `cargo run -p fs-material --example water_freezing_front`.
//! A 10 mm reference layer begins just above its actual atmospheric melting
//! point. The left surface is held at 253.15 K; the right and transverse faces
//! are insulated. Cell-to-cell heat flow, temperature gradients and phase
//! inventories are computed by the existing thermal owner. No front speed,
//! apparent heat capacity or separate latent heat source is prescribed.
//!
//! Liquid conductivity is IAPWS SR6-08(2011); ordinary freshwater-ice conductivity
//! is the named empirical Yen (1981) correlation. Mixed cells explicitly use a
//! series-layer closure with original phase VOLUME fractions. Geometry remains
//! fixed in its reference configuration while source expansion and p*Delta V
//! are reported separately. This is a conduction approximation, not a moving
//! fluid/free-surface or mechanically constrained freezing calculation.

use fs_material::phase::water::freezing::spatial::{
    FreezingSlab, FreezingSlabGeometry, FreezingSlabReservoir,
};
use fs_material::phase::water::freezing::thermal::{
    FreezingThermalControls, FreezingWaterThermalLaw,
};
use fs_material::phase::water::freezing::transport::{
    FreezingWaterConductivity, PhaseConductivityMixture,
};
use fs_material::phase::water::freezing::{FreezingWater, FreezingWaterControls};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{DepositionBudget, StepConfig};
use fs_thermochem::water::ice::IceMeltingSolve;
use fs_thermochem::water::ice::bounds::IceIhHeatCapacityBoundSolve;
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chart = FreezingWater::new(
        Pressure::new(101_325.0),
        FreezingWaterControls {
            minimum_temperature_k: 243.15,
            maximum_temperature_k: 303.15,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            maximum_gibbs_mismatch_j_kg: 30.0,
            melting: IceMeltingSolve {
                pressure_tolerance_pa: 1e-5,
                temperature_tolerance_k: 1e-10,
                maximum_iterations: 96,
            },
        },
        || false,
    )?;
    let thermal_law = FreezingWaterThermalLaw::new(
        chart,
        FreezingThermalControls {
            ice_heat_capacity: IceIhHeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
            liquid_heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )?;
    let conductivity =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Series, || false)?;
    let initial = chart.state_at_temperature(Temperature::new(273.16), || false)?;
    let mut slab = FreezingSlab::new(
        &thermal_law,
        &conductivity,
        FreezingSlabGeometry {
            thickness_m: 0.01,
            area_m2: 0.001,
            cells: 24,
        },
        initial,
        [
            Some(FreezingSlabReservoir {
                temperature_k: 253.15,
                surface_resistance_m2_k_w: 0.0,
            }),
            None,
        ],
        &[],
        DepositionBudget {
            max_cells: 24,
            max_links: 23,
            max_events: 47,
            max_substeps: 100_000,
        },
        || Ok(()),
    )?;
    let config = StepConfig {
        max_dt_s: 0.2,
        energy_tolerance_j: 2e-6,
    };
    let mut heat = 0.0;
    let mut enthalpy = 0.0;
    let mut internal_energy = 0.0;
    let mut pressure_work = 0.0;
    let mut entropy = [0.0; 4];
    let mut max_residual = 0.0_f64;
    let mut steps = 0;
    println!("# 10 mm reference layer, 24 cells, area=0.001 m2, pressure=101325 Pa");
    println!(
        "# original IAPWS ice/IF97 calorics; SR6-08 liquid and Yen1981 ice conductivity; series mixed-cell closure"
    );
    println!("# initial_temperature=273.16 K; left_surface=253.15 K; other faces insulated");
    println!(
        "# melting_temperature={:.12} K; fusion_enthalpy={:.9} J/kg; reference_mass={:.12} kg",
        chart.melting().temperature().value(),
        chart.melting().latent_enthalpy_j_kg(),
        slab.mass_kg()
    );
    println!(
        "time_s,frozen_reference_equiv_mm,ice_kg,liquid_kg,free_ice_m3,free_liquid_m3,min_temperature_K,max_temperature_K,heat_J,delta_H_J,delta_U_J,pressure_work_J,body_entropy_J_K,reservoir_entropy_input_J_K,conductive_entropy_J_K,temporal_entropy_defect_J_K"
    );
    for time in [0.0, 20.0, 60.0, 120.0] {
        if time > slab.time_s() {
            let (next, audit) = slab.advance_to(time, config, || Ok(()))?;
            heat += audit.thermal.external_heat_j;
            enthalpy += audit.source_enthalpy_change_j;
            internal_energy += audit.internal_energy_change_j;
            pressure_work += audit.pressure_work_j;
            steps += audit.thermal.substeps;
            for value in [
                audit.thermal.residual_j,
                audit.source_enthalpy_residual_j,
                audit.internal_energy_residual_j,
            ] {
                max_residual = max_residual.max(value.abs());
            }
            for (total, value) in entropy.iter_mut().zip([
                audit.entropy.thermal_storage_change_j_k,
                audit.entropy.reservoir_entropy_input_j_k,
                audit.entropy.conduction_generation_j_k,
                audit.entropy.temporal_defect_j_k,
            ]) {
                *total += value;
            }
            if audit.entropy.has_unclassified_sources {
                return Err(
                    "surface-conduction example acquired an unclassified heat source".into(),
                );
            }
            slab = next;
        }
        let inventory = slab.inventory(|| Ok(()))?;
        println!(
            "{time:.0},{:.9},{:.12},{:.12},{:.12},{:.12},{:.9},{:.9},{heat:.9},{enthalpy:.9},{internal_energy:.9},{pressure_work:.9},{:.9},{:.9},{:.9},{:.9}",
            1000.0 * inventory.equivalent_frozen_reference_length_m,
            inventory.ice_mass_kg,
            inventory.liquid_mass_kg,
            inventory.free_ice_volume_m3,
            inventory.free_liquid_volume_m3,
            inventory.minimum_temperature_k,
            inventory.maximum_temperature_k,
            entropy[0],
            entropy[1],
            entropy[2],
            entropy[3]
        );
    }
    println!(
        "# final spatial field: location is REFERENCE cell centre; phase volumes are source values"
    );
    println!(
        "reference_x_mm,temperature_K,liquid_mass_fraction,liquid_volume_fraction,conductivity_W_m_K,specific_enthalpy_J_kg"
    );
    let mut mixed_cells = 0;
    for cell in slab.profile(|| Ok(()))? {
        let state = cell.material;
        let k = conductivity.state_at_material_state(state, || false)?;
        let x = state.liquid_mass_fraction();
        mixed_cells += usize::from(x > 0.0 && x < 1.0);
        println!(
            "{:.9},{:.9},{x:.9},{:.9},{:.9},{:.9}",
            1000.0 * cell.reference_centre_m,
            state.temperature_k(),
            state.liquid_volume_fraction(),
            k.conductivity_w_m_k(),
            cell.stored_specific_enthalpy_j_kg
        );
    }
    let inventory = slab.inventory(|| Ok(()))?;
    if inventory.ice_mass_kg <= 0.0
        || inventory.liquid_mass_kg <= 0.0
        || mixed_cells == 0
        || (enthalpy - heat).abs() > 2e-6
        || (internal_energy + pressure_work - heat).abs() > 2e-6
        || entropy[0] - entropy[1] <= 0.0
    {
        return Err("spatial freezing did not produce a balanced resolved phase transition".into());
    }
    println!(
        "# substeps={steps}; maximum_energy_residual={max_residual:.9e} J; body_plus_reservoir_entropy={:.9} J/K; temporal_entropy_defect={:.9e} J/K",
        entropy[0] - entropy[1],
        entropy[3]
    );
    Ok(())
}
