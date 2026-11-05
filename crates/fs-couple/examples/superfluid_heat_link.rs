//! Finite heater -> unsteady viscous helium capillary -> finite thermal sink.
//!
//! cargo run -p fs-couple --features superfluid --example superfluid_heat_link
//!
//! Helium is the sourced 1.8 K SVP reference. Geometry, Kapitza coefficients
//! and finite-cell calorics are authored examples, not an apparatus calibration.
//! Transverse Stokes modes retain normal-profile memory and kinetic energy.
//! Explicit local caloric receivers collect distributed dissipation. Their
//! storage is authored; it does not relinearize the frozen helium reference.

use fs_couple::superfluid::counterflow::{CounterflowHeatControls, CounterflowThermalContact};
use fs_couple::superfluid::counterflow_pair::dissipation::advance_counterflow_heat_pair_with_dissipation;
use fs_duct::superfluid::{ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel};
use fs_material::superfluid::NormalViscousSpectrum;
use fs_material::superfluid::helium::HeliumSvpReference;
use fs_material::superfluid::kapitza::{KapitzaInterface, KapitzaParameters};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const T0: f64 = 1.8;
const MASS: f64 = 1e-9;
const CP: f64 = 1.0;
const INTERFACE_COEFFICIENT_W_M2_K4: f64 = 135.0;
const THETA: f64 = 0.001;
const CELLS: usize = 48;
const LENGTH: f64 = 1.0;
const RADIUS: f64 = 1e-6;
const IDS: [CellId; 2] = [CellId(1), CellId(2)];
const RESOLVED_PROFILE_MODES: usize = 8;
const RECEIVER_MASS: f64 = 1e-12;

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
    let spectrum =
        NormalViscousSpectrum::circular_capillary(material, RADIUS, RESOLVED_PROFILE_MODES)?;
    let wall = spectrum.wall();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    // Authored surface coefficient and actual endpoint area. This example
    // exercises nonlinear boundary transfer; it is not a surface calibration.
    let interface = KapitzaInterface::try_new(KapitzaParameters {
        area_m2: area,
        coefficient_w_m2_k_n: INTERFACE_COEFFICIENT_W_M2_K4,
        exponent: 4.0,
        min_temperature_k: 1.79,
        max_temperature_k: 1.81,
    })?;
    let reference_conductance = interface.tangent_from_w_k(T0)?;
    let mut channel = TwoFluidChannel::try_new_with_normal_profile(
        material,
        LENGTH,
        area,
        CELLS,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1e4,
            max_temperature_perturbation_k: 0.002,
            max_component_speed_m_s: 0.001,
        },
        spectrum,
    )?;
    // Independent steady two-fluid equations: mu'=p'/rho-s*T'=0,
    // gamma*vn=-rho*s*grad(T'), q=rho*s*T0*vn at zero mass current.
    // The explicit ports sit at the end cells, with N-1 resolved dual volumes.
    // Only this initialization uses the reference interface tangent; every
    // transient step solves both original nonlinear surface laws together.
    let gamma = material.normal_viscosity_pa_s() / wall.permeability_m2();
    let entropy_density = material.density_kg_m3() * material.entropy_j_kg_k();
    let conductivity = entropy_density * entropy_density * T0 / gamma;
    let resolved_length = LENGTH * (CELLS - 1) as f64 / CELLS as f64;
    let initial_power =
        2.0 * THETA / (2.0 / reference_conductance + resolved_length / (conductivity * area));
    let left_theta = THETA - initial_power / reference_conductance;
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
    // Seed the developed modal shape without changing either section mean.
    // Its positive excess kinetic energy is part of the channel's storage.
    channel.initialize_normal_profile_steady(|| false)?;

    let thermal_config = StepConfig {
        max_dt_s: 0.001,
        energy_tolerance_j: 1e-23,
    };
    let destinations: Vec<_> = (0..CELLS).map(|i| CellId(100 + i as u64)).collect();
    let events = (0..CELLS + 2)
        .map(|i| {
            let endpoint = i < 2;
            let position = if endpoint {
                i as f64 * LENGTH
            } else {
                (i - 2) as f64 * LENGTH / CELLS as f64
            };
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: if endpoint {
                        IDS[i]
                    } else {
                        destinations[i - 2]
                    },
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(1),
                    mass_kg: if endpoint { MASS } else { RECEIVER_MASS },
                    specific_h_j_kg: CP
                        * match i {
                            0 => THETA,
                            1 => -THETA,
                            _ => 0.0,
                        },
                    geometry: CellGeometry {
                        start_m: [position, 0.0, 0.0],
                        end_m: [position + 1e-6, 0.0, 0.0],
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
            }
        })
        .collect();
    let original = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &FiniteCaloricCell)],
        DepositionBudget {
            max_cells: CELLS + 2,
            max_links: 0,
            max_events: CELLS + 2,
            max_substeps: 4,
        },
    )?;
    let mut thermal = original.advance_to(0.0, thermal_config, || Ok(()))?.0;
    let contacts = [ChannelEnd::Left, ChannelEnd::Right].map(|end| {
        CounterflowThermalContact::kapitza(
            end,
            IDS[usize::from(end == ChannelEnd::Right)],
            interface,
        )
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
    let mut received_dissipation = vec![0.0; CELLS];
    let mut maximum_dissipation_residual: f64 = 0.0;
    println!("sourced_helium_T0_K={T0}, authored_radius_m={RADIUS}, authored_length_m={LENGTH}");
    println!(
        "normal_relaxation_rate_per_s={:.9e}",
        material.normal_wall_decay_rate_s_inv(wall)?
    );
    println!("authored_kapitza_coefficient_W_m2_K4={INTERFACE_COEFFICIENT_W_M2_K4}");
    println!("interface_area_m2={area:.9e}, reference_tangent_W_K={reference_conductance:.9e}");
    println!("linearized_initial_heat_power_W={initial_power:.9e}");
    println!("resolved_transverse_modes={RESOLVED_PROFILE_MODES}, matched_tail_modes=1");
    println!("time_s,heater_K,sink_K,heat_from_heater_W,heat_to_sink_W");
    for step in 0..1000 {
        let step_report = advance_counterflow_heat_pair_with_dissipation(
            &mut channel,
            &mut thermal,
            contacts,
            &destinations,
            dt,
            MutualFriction::None,
            controls,
            || false,
        )?;
        let report = step_report.contacts;
        maximum_dissipation_residual =
            maximum_dissipation_residual.max(step_report.dissipation_transfer_residual_j.abs());
        for (i, received) in received_dissipation.iter_mut().enumerate() {
            *received +=
                channel.last_normal_viscous_heat_j()[i] + channel.last_mutual_friction_heat_j()[i];
        }
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
    let mut stored_dissipation = 0.0;
    for (i, &id) in destinations.iter().enumerate() {
        let stored = RECEIVER_MASS * thermal.cells()[&id].specific_h_j_kg;
        assert!((stored - received_dissipation[i]).abs() < 1e-23);
        stored_dissipation += stored;
    }
    assert!((stored_dissipation - channel.cumulative_normal_viscous_heat_j()).abs() < 1e-23);
    println!("max_reference_heat_transfer_residual_J={maximum_transfer_residual:.3e}");
    println!(
        "normal_viscous_heat_J={:.9e}",
        channel.cumulative_normal_viscous_heat_j()
    );
    println!("native_distributed_heat_storage_J={stored_dissipation:.9e}");
    println!("max_dissipation_transfer_residual_J={maximum_dissipation_residual:.3e}");
    println!(
        "profile_excess_kinetic_energy_J={:.9e}",
        channel.normal_profile_excess_energy_j()?
    );
    println!("wave_availability_J={:.9e}", channel.wave_energy_j()?);
    Ok(())
}
