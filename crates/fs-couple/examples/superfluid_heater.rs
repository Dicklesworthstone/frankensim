//! Finite heater -> thermal counterflow -> first/second sound at sourced 1.8 K.
//!
//! cargo run -p fs-couple --features superfluid --example superfluid_heater
//!
//! Helium properties come from the existing Donnelly–Barenghi reference. Heater
//! capacitance, interface conductance, geometry, power and frozen vortex tangle
//! are authored demonstration inputs, not calibrated apparatus measurements.

use fs_couple::superfluid::counterflow::{
    CounterflowHeatControls, CounterflowThermalContact, advance_counterflow_heat,
};
use fs_duct::superfluid::{ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel};
use fs_material::superfluid::helium::HeliumSvpReference;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const T0: f64 = 1.8;
const ID: CellId = CellId(1);
const MASS_KG: f64 = 1e-6;
const CP_J_KG_K: f64 = 50.0;
const SOURCE_W: f64 = 1e-7;

struct Heater;
impl ThermalLaw for Heater {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        let t = T0 + u / CP_J_KG_K;
        if t.is_finite() && (1.79..=1.81).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP_J_KG_K
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference = HeliumSvpReference::at_1_8_k()?;
    let mut channel = TwoFluidChannel::try_new(
        *reference.state(),
        0.05,
        1e-6,
        32,
        LinearizationLimits {
            max_pressure_perturbation_pa: 5.0,
            max_temperature_perturbation_k: 0.001,
            max_component_speed_m_s: 0.001,
        },
    )?;
    let thermal_config = StepConfig {
        max_dt_s: 1e-4,
        energy_tolerance_j: 1e-18,
    };
    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: ID,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: MASS_KG,
                specific_h_j_kg: CP_J_KG_K * 0.001,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1e-3, 0.0, 0.0],
                    volume_m3: 1e-9,
                    surface_area_m2: 6e-6,
                },
                boundary: Boundary {
                    ambient_k: T0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: SOURCE_W,
                },
            }),
        }],
        &[(MaterialId(1), &Heater)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 4,
        },
    )?;
    let mut thermal = initial.advance_to(0.0, thermal_config, || Ok(()))?.0;
    let initial_energy = thermal.energy_j();
    let contact = CounterflowThermalContact {
        end: ChannelEnd::Left,
        cell: ID,
        conductance_w_k: 0.002,
    };
    let controls = CounterflowHeatControls {
        thermal: thermal_config,
        max_iterations: 64,
        max_abs_heat_power_w: 1e-4,
        max_temperature_departure_k: 0.005,
        contact_power_tolerance_w: 1e-14,
        cell_energy_tolerance_j: 1e-20,
        transfer_energy_tolerance_j: 1e-19,
    };
    let friction = MutualFriction::IsotropicVortexTangle {
        circulation_quantum_m2_s: 9.97e-8,
        line_density_per_m2: 1e8,
        dissipative_b: 1.0,
    };
    let dt = 0.5 * channel.maximum_step_s();
    let mut transferred_j = 0.0;
    let mut maximum_transfer_defect: f64 = 0.0;
    let modes = reference.state().sound_modes();
    println!(
        "sourced helium: T0={T0} K, fast={:.6} m/s, slow={:.6} m/s",
        modes.fast_m_s, modes.slow_m_s
    );
    println!(
        "time_s,heater_K,inlet_temperature_perturbation_K,mid_temperature_perturbation_K,heat_into_helium_W"
    );
    for step in 0..2000 {
        let report = advance_counterflow_heat(
            &mut channel,
            &mut thermal,
            contact,
            dt,
            friction,
            controls,
            || false,
        )?;
        transferred_j -= report.heat_to_thermal_j;
        maximum_transfer_defect =
            maximum_transfer_defect.max(report.reference_heat_transfer_residual_j.abs());
        if step % 400 == 399 {
            println!(
                "{:.9e},{:.9},{:.9e},{:.9e},{:.9e}",
                channel.time_s(),
                thermal.temperature_k(ID)?,
                channel.cell(0).unwrap().temperature_perturbation_k,
                channel.cell(16).unwrap().temperature_perturbation_k,
                -report.heat_outward_w
            );
        }
    }
    let volume = channel.area_m2() * channel.length_m() / channel.cells() as f64;
    let reference_heat: f64 = (0..channel.cells())
        .map(|i| {
            let c = channel.cell(i).unwrap();
            T0 * volume
                * (reference.state().density_kg_m3() * c.specific_entropy_perturbation_j_kg_k
                    + reference.state().entropy_j_kg_k() * c.density_perturbation_kg_m3)
        })
        .sum();
    let thermal_balance =
        thermal.energy_j() - initial_energy + transferred_j - SOURCE_W * thermal.time_s();
    assert!((reference_heat - transferred_j).abs() < 1e-16);
    assert!(thermal_balance.abs() < 1e-16);
    println!("reference_heat_into_helium_J={transferred_j:.9e}");
    println!("thermal_balance_residual_J={thermal_balance:.3e}");
    println!("max_reference_heat_transfer_residual_J={maximum_transfer_defect:.3e}");
    println!(
        "exported_vortex_friction_heat_J={:.9e}",
        channel.cumulative_mutual_friction_heat_j()
    );
    println!("wave_availability_J={:.9e}", channel.wave_energy_j()?);
    Ok(())
}
