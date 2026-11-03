//! Isolated helium-II waves, two native pistons, and original thermal cells.
//!
//! Run: cargo run -p fs-couple --features superfluid --example superfluid_piston_pair
//! Optional first argument: channel cells (default 64). The 1.8 K helium
//! reference is sourced; geometry, excitation, piston dimensions and friction
//! are authored demonstration inputs, not a calibrated cryogenic apparatus.
//! Reference pressure is balanced outside both pistons. A local isochoric
//! caloric tangent receives dissipated energy, without updating the frozen
//! wave reference or claiming a nonlinear helium equation of state.

use fs_couple::superfluid::FrictionHeatControls;
use fs_couple::superfluid::piston::{BalancedPistonBoundary, PistonControls};
use fs_couple::superfluid::piston_pair::advance_balanced_piston_pair_with_heat;
use fs_duct::superfluid::{
    BodyForce, ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::helium::HeliumSvpReference;
use fs_math::{dd::Dd, det};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

struct IsochoricReference {
    temperature_k: f64,
    cv_j_kg_k: f64,
}
impl ThermalLaw for IsochoricReference {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, specific_energy_j_kg: f64) -> Result<f64, DepositionError> {
        if !specific_energy_j_kg.is_finite()
            || !(0.0..=self.cv_j_kg_k * 1e-3).contains(&specific_energy_j_kg)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(self.temperature_k + specific_energy_j_kg / self.cv_j_kg_k)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cv_j_kg_k
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cells = std::env::args().nth(1).map_or(Ok(64usize), |s| s.parse())?;
    if !(8..=4096).contains(&cells) {
        return Err("cell count must be between 8 and 4096".into());
    }
    let reference = HeliumSvpReference::at_1_8_k()?;
    let material = *reference.state();
    let length_m = 0.1;
    let area_m2 = 1e-4;
    let dx_m = length_m / cells as f64;
    let radius_m = det::sqrt(area_m2 / core::f64::consts::PI);
    let mut channel = TwoFluidChannel::try_new(
        material,
        length_m,
        area_m2,
        cells,
        LinearizationLimits {
            max_pressure_perturbation_pa: 100.0,
            max_temperature_perturbation_k: 1e-3,
            max_component_speed_m_s: 0.1,
        },
    )?;
    let pressure: Vec<_> = (0..cells)
        .map(|i| det::cos(core::f64::consts::PI * (i as f64 + 0.5) / cells as f64))
        .collect();
    let temperature: Vec<_> = (0..cells)
        .map(|i| 2e-4 * det::cos(core::f64::consts::TAU * (i as f64 + 0.5) / cells as f64))
        .collect();
    let still = vec![0.0; cells + 1];
    channel.initialize(&pressure, &temperature, &still, &still, || false)?;

    // Two unequal centered disks; both positions and momenta evolve natively.
    let masses_kg = [0.002, 0.003];
    let thickness_m = [0.001, 0.0015];
    let properties = [0, 1].map(|i| {
        let axial = 0.5 * masses_kg[i] * radius_m * radius_m;
        let transverse =
            masses_kg[i] * (3.0 * radius_m * radius_m + thickness_m[i] * thickness_m[i]) / 12.0;
        MassProperties::new(
            masses_kg[i],
            Vec3::ZERO,
            Vec3::new(axial, transverse, transverse),
        )
    });
    let [left_properties, right_properties] = properties;
    let properties = [left_properties?, right_properties?];
    let positions = [Vec3::ZERO, Vec3::new(length_m, 0.0, 0.0)];
    let mut bodies = [
        RigidBodyState::new(
            Pose::new(positions[0], UnitQuaternion::IDENTITY)?,
            Vec3::ZERO,
            Vec3::ZERO,
        )?,
        RigidBodyState::new(
            Pose::new(positions[1], UnitQuaternion::IDENTITY)?,
            Vec3::ZERO,
            Vec3::ZERO,
        )?,
    ];
    let boundaries = [0, 1].map(|i| BalancedPistonBoundary {
        end: if i == 0 {
            ChannelEnd::Left
        } else {
            ChannelEnd::Right
        },
        positive_axis_world: Vec3::new(1.0, 0.0, 0.0),
        equilibrium_position_world: positions[i],
        max_displacement_m: 0.001,
    });
    let duration_s = 0.005;
    let steps = (duration_s / (0.8 * channel.maximum_step_s())).ceil() as usize;
    let dt_s = duration_s / steps as f64;
    let controls = PistonControls {
        momentum_tolerance_n_s: 1e-14,
        energy_tolerance_j: 1e-16,
        velocity_tolerance_m_s: 1e-11,
    };
    let heat_controls = FrictionHeatControls {
        thermal: StepConfig {
            max_dt_s: dt_s,
            energy_tolerance_j: 1e-18,
        },
        cell_heat_tolerance_j: 1e-20,
        heat_energy_tolerance_j: 1e-18,
        combined_energy_tolerance_j: 1e-16,
    };
    let caloric = IsochoricReference {
        temperature_k: material.temperature_k(),
        cv_j_kg_k: material.specific_heat_cv_j_kg_k(),
    };
    let destinations: Vec<_> = (0..cells).map(|i| CellId(i as u64 + 1)).collect();
    let cell_mass_kg = material.density_kg_m3() * area_m2 * dx_m;
    let schedule = destinations
        .iter()
        .enumerate()
        .map(|(i, &id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: cell_mass_kg,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [i as f64 * dx_m, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * dx_m, 0.0, 0.0],
                    volume_m3: area_m2 * dx_m,
                    surface_area_m2: core::f64::consts::TAU * radius_m * dx_m,
                },
                boundary: Boundary {
                    ambient_k: material.temperature_k(),
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let original_thermal = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), &caloric as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: cells,
            max_links: 0,
            max_events: cells,
            max_substeps: 32,
        },
    )?;
    let mut thermal = original_thermal
        .advance_to(0.0, heat_controls.thermal, || Ok(()))?
        .0;
    let initial_energy_j = channel.wave_energy_j()?;
    let initial_momentum = channel.resolved_momentum_kg_m_s()?;
    let half_mass = channel.boundary_half_mass_kg();
    let mut transferred_heat_j = Dd::ZERO;
    for _ in 0..steps {
        let report = advance_balanced_piston_pair_with_heat(
            &mut channel,
            &mut bodies,
            properties,
            boundaries,
            &mut thermal,
            &destinations,
            dt_s,
            BodyForce::None,
            [0.0; 2],
            MutualFriction::Linear {
                gamma_kg_m3_s: 1000.0,
            },
            controls,
            heat_controls,
            || false,
        )?;
        transferred_heat_j =
            transferred_heat_j + Dd::from_f64(report.heat.heat_intake.external_heat_j);
    }

    // Reconstruct cumulative physical balances from the final original states,
    // independently of the step residuals and without a fictitious support load.
    let mut energy = Dd::from_f64(channel.wave_energy_j()?);
    let mut momentum = Dd::from_f64(channel.resolved_momentum_kg_m_s()?);
    for i in 0..2 {
        let velocity = bodies[i].center_of_mass_velocity_world(properties[i])?.x;
        energy = energy
            + Dd::from_f64(bodies[i].kinetic_energy(properties[i])?)
            + Dd::from_f64(0.5 * half_mass * velocity * velocity);
        momentum = momentum
            + Dd::from_f64(bodies[i].linear_momentum_world().x)
            + Dd::from_f64(half_mass * velocity);
        println!("piston_{i}_mass_kg={:.12e}", masses_kg[i]);
        println!(
            "piston_{i}_displacement_x_m={:.12e}",
            bodies[i].pose().position_world().x - positions[i].x
        );
        println!("piston_{i}_velocity_x_m_s={velocity:.12e}");
    }
    let thermal_energy = thermal.cells().values().fold(Dd::ZERO, |sum, cell| {
        sum + Dd::from_f64(cell.mass_kg * cell.specific_h_j_kg)
    });
    let maximum_dissipative_temperature_rise_k = thermal
        .cells()
        .values()
        .map(|cell| cell.specific_h_j_kg / caloric.cv_j_kg_k)
        .fold(0.0_f64, f64::max);
    let energy_defect = (energy + thermal_energy - Dd::from_f64(initial_energy_j)).to_f64();
    let momentum_defect = (momentum - Dd::from_f64(initial_momentum)).to_f64();
    let transfer_defect = (thermal_energy - transferred_heat_j).to_f64();
    let momentum_scale = det::sqrt(
        2.0 * initial_energy_j
            * (material.density_kg_m3() * area_m2 * length_m + masses_kg.iter().sum::<f64>()),
    );
    if energy_defect.abs() > 1e-8 * initial_energy_j
        || momentum_defect.abs() > 1e-10 * momentum_scale
        || transfer_defect.abs() > 1e-9 * transferred_heat_j.to_f64()
        || thermal.time_s() != channel.time_s()
    {
        return Err("cumulative momentum, energy, heat transfer or clock balance failed".into());
    }
    println!(
        "helium_reference_temperature_k={:.9}",
        material.temperature_k()
    );
    println!(
        "cells={cells} steps={steps} simulated_time_s={:.12e}",
        channel.time_s()
    );
    println!(
        "fast_sound_m_s={:.9} slow_sound_m_s={:.9}",
        material.sound_modes().fast_m_s,
        material.sound_modes().slow_m_s
    );
    println!("boundary_fluid_half_mass_kg={half_mass:.12e}");
    println!("deposited_friction_heat_j={:.12e}", thermal_energy.to_f64());
    println!(
        "maximum_dissipative_temperature_rise_k={maximum_dissipative_temperature_rise_k:.12e}"
    );
    println!(
        "relative_combined_energy_defect={:.12e}",
        energy_defect / initial_energy_j
    );
    println!("combined_momentum_defect_n_s={momentum_defect:.12e}");
    println!("heat_transfer_defect_j={transfer_defect:.12e}");
    Ok(())
}
