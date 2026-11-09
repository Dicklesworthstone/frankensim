//! Passive reservoir withdrawal/return with native liquid-column separation.
//!
//! cargo run -p fs-couple --features cavitation --example cavitating_reservoir
//!
//! IF97 water at 293.15 K retains its original density, sound speed, vapor
//! pressure and vapor-density admission. The one-metre, radius-0.01-m column
//! starts from an authored quadratic pressure profile with one endpoint at
//! saturation. Two constant-pressure reservoirs are selected at t=1 us through
//! the SAME declared resistance. This ideal switching event is split exactly;
//! Q=(p_res-p_port)/R responds to every actual liquid stage and cavity event.
//!
//! Geometry, initial pressure profile and resistance are demonstration inputs,
//! not measured valve data or a pump curve. The original small-strain DG model
//! retains equilibrium, isothermal, dilute pure-vapor cavities. Resistance loss
//! is reported as an external heat port and is not added to liquid temperature.
//! The model does not resolve valve inertia, nucleation delay or latent heat.

use fs_couple::cavitation::duct::pockets::driven::DuctReservoirResistance;
use fs_couple::cavitation::duct::pockets::{CavitatingDuctPockets, DuctPocket};
use fs_couple::cavitation::duct::{DuctCavityLimits, DuctCavityNumerics};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId};
use fs_duct::{Duct, Segment};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().nth(1).is_some() {
        return Err("usage: cavitating_reservoir (no arguments)".into());
    }
    let water = SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(293.15))?,
        0.072,
    )?;
    let vapor = water.state().saturation_pressure_pa();
    let reference = 100_000.0;
    let radius = 0.01;
    let area = std::f64::consts::PI * radius * radius;
    let duct = DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 1.0,
            }],
        },
        &water,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(1)),
        DgBoundary::Port(PortId(2)),
    )?;
    let initial =
        duct.project_initial_state(|x| vapor - reference + 100_000.0 * x * (2.0 - x), |_| 0.0);
    let mut coupled = CavitatingDuctPockets::new_with_driven_port(
        duct,
        initial,
        0.0,
        &[DuctPocket {
            port: PortId(1),
            initial_volume_m3: 0.0,
            maximum_volume_m3: 1e-6,
        }],
        &[],
        &[],
        PortId(2),
        DuctCavityLimits {
            reference_pressure_pa: reference,
            maximum_cavity_volume_m3: 2e-6,
            maximum_liquid_strain: 0.001,
            maximum_interface_mach: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
    )?;
    let controls = DuctCavityNumerics {
        relative_tolerance: 1e-10,
        absolute_energy_tolerance_j: 1e-11,
        absolute_pressure_tolerance_pa: 1e-3,
        absolute_volume_tolerance_m3: 1e-14,
        maximum_step_s: 1e-7,
        event_time_tolerance_s: 1e-13,
        event_scan_substeps: 8,
        maximum_attempts: 30_000,
        maximum_events: 32,
    };
    let resistance = 72_000.0 / (area * 0.1);
    let initial_drive_pressure = vapor + 100_000.0;
    let withdrawal =
        DuctReservoirResistance::try_new(initial_drive_pressure - 72_000.0, resistance)?;
    let return_flow =
        DuctReservoirResistance::try_new(initial_drive_pressure + 72_000.0, resistance)?;
    let mut volume = 0.0;
    let mut acoustic_work = 0.0;
    let mut reservoir_work = 0.0;
    let mut resistance_loss = 0.0;
    let mut ruptures = 0;
    let mut collapses = 0;
    println!(
        "time_s,reservoir_pressure_pa,port_pressure_pa,inflow_m3_s,cavity_volume_m3,mode,transported_volume_m3,acoustic_work_j,reservoir_absolute_work_j,resistance_loss_j,ruptures,collapses,liquid_balance_residual_j,reservoir_balance_residual_j"
    );
    for index in 1..=36 {
        let time = f64::from(index) * (1e-6 / 12.0);
        let reservoir = if index <= 12 { withdrawal } else { return_flow };
        let step = coupled.advance_with_reservoir(time, &controls, reservoir, &mut || false)?;
        let drive = step
            .drive
            .as_ref()
            .ok_or("declared hydraulic port is missing")?;
        let balance = drive
            .reservoir
            .as_ref()
            .ok_or("reservoir balance is missing")?;
        volume += drive.transported_volume_m3;
        acoustic_work += drive.acoustic_work_j;
        reservoir_work += balance.supplied_absolute_work_j;
        resistance_loss += balance.dissipated_energy_j;
        ruptures += step.pockets[0].rupture_events;
        collapses += step.pockets[0].collapse_events;
        println!(
            "{time:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:?},{volume:.12e},{acoustic_work:.12e},{reservoir_work:.12e},{resistance_loss:.12e},{ruptures},{collapses},{:.12e},{:.12e}",
            reservoir.absolute_pressure_pa(),
            drive.end_absolute_pressure_pa,
            drive.inflow_m3_s,
            step.pockets[0].volume_m3,
            step.pockets[0].mode,
            step.energy_balance_residual_j,
            balance.energy_balance_residual_j
        );
    }
    Ok(())
}
