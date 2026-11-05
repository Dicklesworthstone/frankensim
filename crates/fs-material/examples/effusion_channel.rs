//! Synthetic finite-reservoir molecular-channel exchange with a held wall.
//!
//! Run: cargo run -p fs-material --example effusion_channel -- 450 200 0.01
//! Compare wall temperatures 250 and 900 K; halve dt for time refinement.
//! Optional fourth argument adds a SEPARATE thin-aperture area [m2].
//! Geometry, gas coefficients and the minimum mean-free-path bound are authored
//! demonstration inputs, not measurements. No collision cross sections are
//! inferred. This is a fixed-step, first-order model with a prescribed thermostat;
//! cumulative wall heat is an account, not a whole-window error certificate.

#[path = "material_gas_spring/source.rs"]
mod source;

use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        channel::{CylindricalMolecularGeometry, DiffuseMolecularChannel},
    },
    network::{
        SpeciesNetworkSolve,
        transport::{
            GasTransportPaths, GasTransportSolve,
            effusion::{
                EffusionNetworkLink,
                channel::{DiffuseChannelNetworkLink, step_candidate},
            },
        },
    },
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_math::dd::Dd;
use std::{collections::BTreeMap, io::Write};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "usage: effusion_channel [wall_temperature_k=450] [steps=200] [dt_s=0.01] [separate_aperture_area_m2=0]";
const DIAMETER_M: f64 = 20e-6;
const LENGTH_M: f64 = 100e-6;
const AREA_M2: f64 = 1e-6;
const LAMBDA_MIN_M: f64 = 1e-3;
const QUALIFICATION: &str = "SYNTHETIC: stationary Maxwellian gases; collisionless circular bores; lambda_min=0.001 m across all retained mixtures and authored T/p support; no empirical validation";

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(5).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() > 4 {
        return Err(USAGE.into());
    }
    let wall_temperature_k: f64 = args.first().map_or(Ok(450.0), |v| v.parse())?;
    let steps: u32 = args.get(1).map_or(Ok(200), |v| v.parse())?;
    let dt: f64 = args.get(2).map_or(Ok(0.01), |v| v.parse())?;
    let parallel_area: f64 = args.get(3).map_or(Ok(0.0), |v| v.parse())?;
    if !(200.0..=1800.0).contains(&wall_temperature_k)
        || !(1..=10_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || dt * f64::from(steps) > 3600.0
        || !parallel_area.is_finite()
        || parallel_area < 0.0
    {
        return Err("wall T must be 200..1800 K, steps 1..10000, dt positive with duration <=3600 s, and separate area finite/nonnegative".into());
    }

    let model = source::load_with_pressure_support(
        include_str!("material_gas_spring/mixture-a.tsv"),
        "built-in synthetic mixture-a molecular-channel demonstration",
        (1e-3, 2.0),
    )?;
    let transmissions = model
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                1.0,
                "synthetic open entrance mouth; channel geometry owns transmission",
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let mouth = EffusiveInterface::try_new(
        model.clone(),
        transmissions,
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1e-3, 2.0),
        },
        EffusiveRegime::try_new(DIAMETER_M, 0.0, LAMBDA_MIN_M, QUALIFICATION)?,
    )?;
    let channel = DiffuseMolecularChannel::try_new(
        mouth,
        CylindricalMolecularGeometry::try_new(DIAMETER_M, LENGTH_M, LAMBDA_MIN_M, QUALIFICATION)?,
    )?;
    let mut gases = BTreeMap::from([
        (
            0,
            ClosedGasParcel::new(model.clone(), source::point(600.0, 1.0)?, 0.001)?,
        ),
        (
            1,
            ClosedGasParcel::new(model, source::point(300.0, 0.5)?, 0.001)?,
        ),
    ]);
    let energy = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let endpoints = BTreeMap::from([(0, energy), (1, energy)]);
    let controls = GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-16,
            energy_tolerance_j: 1e-12,
        },
        entropy_tolerance_j_k: 1e-13,
    };
    let separate_aperture = [EffusionNetworkLink {
        id: 100,
        parcels: [0, 1],
        interface: channel.mouth(),
        area_m2: parallel_area,
    }];
    let paths = GasTransportPaths {
        effusion: &separate_aperture,
        ..GasTransportPaths::default()
    };
    let links = [DiffuseChannelNetworkLink {
        id: 1,
        parcels: [0, 1],
        channel: &channel,
        area_m2: AREA_M2,
        wall_temperature_k,
    }];

    eprintln!("SYNTHETIC two-reservoir example; each rigid volume is 0.001 m3.");
    eprintln!("D={DIAMETER_M} m, L={LENGTH_M} m, aggregate mouth area={AREA_M2} m2.");
    eprintln!(
        "Estimated transmission={:.9}; held wall={wall_temperature_k} K.",
        channel.geometry().transmission_probability()
    );
    eprintln!("Qualification: {QUALIFICATION}");
    eprintln!("Positive wall heat enters the thermostat. Refusals publish no step.");
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    writeln!(
        output,
        "time_s,left_pressure_pa,left_temperature_k,left_mass_kg,right_pressure_pa,right_temperature_k,right_mass_kg,step_heat_into_wall_j,total_heat_into_wall_j,step_energy_residual_j,step_total_entropy_change_j_k"
    )?;
    let mut total_wall_heat = Dd::ZERO;
    for index in 0..=steps {
        let (mut wall_heat, mut residual, mut entropy) = (0.0, 0.0, 0.0);
        if index > 0 {
            let candidate = step_candidate(
                &gases, paths, &links, dt, &endpoints, controls, || false,
            ).map_err(|error| format!(
                "step {index} refused: {error}; refine dt for a time-step refusal; source/regime limits require admissible physical inputs"
            ))?;
            let (next, step) = candidate.into_parts();
            wall_heat = step.heat_into_walls_j;
            residual = step.energy_residual_j;
            entropy = step.total_entropy_change_j_k;
            total_wall_heat = total_wall_heat + Dd::from_f64(wall_heat);
            gases = next;
        }
        let left = gases[&0].state().state();
        let right = gases[&1].state().state();
        writeln!(
            output,
            "{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            f64::from(index) * dt,
            left.pressure,
            left.temperature,
            gases[&0].mass_kg(),
            right.pressure,
            right.temperature,
            gases[&1].mass_kg(),
            wall_heat,
            total_wall_heat.to_f64(),
            residual,
            entropy,
        )?;
    }
    output.flush()?;
    Ok(())
}
