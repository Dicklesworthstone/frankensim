//! Finite heater -> viscous helium capillary -> finite thermal sink.
//!
//! cargo run -p fs-couple --features superfluid --example superfluid_heat_link
//!
//! Helium is the sourced 1.8 K SVP reference. Geometry, interface conductances
//! and finite-cell calorics are authored examples, not an apparatus calibration.
//! A quasisteady initial profile avoids imposing a sudden transverse startup
//! on the reduced low-Womersley normal-wall closure.

use fs_couple::superfluid::counterflow::{CounterflowHeatControls, CounterflowThermalContact};
use fs_couple::superfluid::counterflow_pair::advance_counterflow_heat_pair;
use fs_duct::superfluid::{ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel};
use fs_material::superfluid::NormalViscousWall;
use fs_material::superfluid::helium::HeliumSvpReference;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const T0: f64 = 1.8;
const MASS: f64 = 1e-9;
const CP: f64 = 1.0;
const G: f64 = 1e-8;
const THETA: f64 = 0.001;
const CELLS: usize = 48;
const LENGTH: f64 = 1.0;
const RADIUS: f64 = 1e-6;
const IDS: [CellId; 2] = [CellId(1), CellId(2)];

struct FiniteCaloricCell;
impl ThermalLaw for FiniteCaloricCell {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        let t = T0 + u / CP;
        if t.is_finite() && (1.79..=1.81).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference = HeliumSvpReference::at_1_8_k()?;
    let material = *reference.state();
    let wall = NormalViscousWall::circular_capillary(RADIUS)?;
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let mut channel = TwoFluidChannel::try_new_with_normal_wall(
        material,
        LENGTH,
        area,
        CELLS,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1e4,
            max_temperature_perturbation_k: 0.002,
            max_component_speed_m_s: 0.001,
        },
        wall,
    )?;
    // Independent steady two-fluid equations: mu'=p'/rho-s*T'=0,
    // gamma*vn=-rho*s*grad(T'), q=rho*s*T0*vn at zero mass current.
    // The explicit ports sit at the end cells, with N-1 resolved dual volumes.
    let gamma = material.normal_viscosity_pa_s() / wall.permeability_m2();
    let entropy_density = material.density_kg_m3() * material.entropy_j_kg_k();
    let conductivity = entropy_density * entropy_density * T0 / gamma;
    let resolved_length = LENGTH * (CELLS - 1) as f64 / CELLS as f64;
    let initial_power = 2.0 * THETA / (2.0 / G + resolved_length / (conductivity * area));
    let left_theta = THETA - initial_power / G;
    let temperature: Vec<_> = (0..CELLS)
        .map(|i| left_theta * (1.0 - 2.0 * i as f64 / (CELLS - 1) as f64))
        .collect();
    let pressure: Vec<_> = temperature
        .iter()
        .map(|theta| entropy_density * theta)
        .collect();
    let normal = initial_power / (area * entropy_density * T0);
    let superfluid = -material.normal_fraction() / material.superfluid_fraction() * normal;
    channel.initialize(
        &pressure,
        &temperature,
        &vec![normal; CELLS + 1],
        &vec![superfluid; CELLS + 1],
        || false,
    )?;

    let thermal_config = StepConfig {
        max_dt_s: 0.001,
        energy_tolerance_j: 1e-23,
    };
    let events = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: CP * if i == 0 { THETA } else { -THETA },
                geometry: CellGeometry {
                    start_m: [i as f64 * LENGTH, 0.0, 0.0],
                    end_m: [i as f64 * LENGTH + 1e-6, 0.0, 0.0],
                    volume_m3: 1e-12,
                    surface_area_m2: 1e-8,
                },
                boundary: Boundary {
                    ambient_k: T0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let original = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &FiniteCaloricCell)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 4,
        },
    )?;
    let mut thermal = original.advance_to(0.0, thermal_config, || Ok(()))?.0;
    let contacts = [ChannelEnd::Left, ChannelEnd::Right].map(|end| CounterflowThermalContact {
        end,
        cell: IDS[usize::from(end == ChannelEnd::Right)],
        conductance_w_k: G,
    });
    let controls = CounterflowHeatControls {
        thermal: thermal_config,
        max_iterations: 96,
        max_abs_heat_power_w: 1e-9,
        max_temperature_departure_k: 0.005,
        contact_power_tolerance_w: 1e-22,
        cell_energy_tolerance_j: 1e-25,
        transfer_energy_tolerance_j: 1e-23,
    };
    let dt = 0.5 * channel.maximum_step_s();
    let mut transfer = [0.0; 2];
    let mut maximum_transfer_residual: f64 = 0.0;
    println!("sourced_helium_T0_K={T0}, authored_radius_m={RADIUS}, authored_length_m={LENGTH}");
    println!(
        "normal_relaxation_rate_per_s={:.9e}",
        material.normal_wall_decay_rate_s_inv(wall)?
    );
    println!("initial_steady_heat_power_W={initial_power:.9e}");
    println!("time_s,heater_K,sink_K,heat_from_heater_W,heat_to_sink_W");
    for step in 0..1000 {
        let report = advance_counterflow_heat_pair(
            &mut channel,
            &mut thermal,
            contacts,
            dt,
            MutualFriction::None,
            controls,
            || false,
        )?;
        for i in 0..2 {
            transfer[i] += report.heat_to_thermal_j[i];
        }
        maximum_transfer_residual =
            maximum_transfer_residual.max(report.reference_heat_transfer_residual_j.abs());
        if step % 200 == 199 {
            println!(
                "{:.9e},{:.9},{:.9},{:.9e},{:.9e}",
                channel.time_s(),
                thermal.temperature_k(IDS[0])?,
                thermal.temperature_k(IDS[1])?,
                -report.heat_outward_w[0],
                report.heat_outward_w[1]
            );
        }
    }
    for i in 0..2 {
        let initial_h = CP * if i == 0 { THETA } else { -THETA };
        let residual = MASS * (thermal.cells()[&IDS[i]].specific_h_j_kg - initial_h) - transfer[i];
        assert!(residual.abs() < 1e-23);
        println!("thermal_cell_{i}_balance_residual_J={residual:.3e}");
    }
    assert!(transfer[0] < 0.0 && transfer[1] > 0.0);
    assert!(channel.cumulative_normal_viscous_heat_j() > 0.0);
    println!("max_reference_heat_transfer_residual_J={maximum_transfer_residual:.3e}");
    println!(
        "exported_normal_viscous_heat_J={:.9e}",
        channel.cumulative_normal_viscous_heat_j()
    );
    println!("wave_availability_J={:.9e}", channel.wave_energy_j()?);
    Ok(())
}
