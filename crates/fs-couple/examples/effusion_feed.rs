//! Finite vessel filling and discharge through distinct qualified molecular ports.
//!
//! Run from the repository root:
//! `cargo run -p fs-couple --example effusion_feed --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv
//! 0.01 1e-7 1e-7 1e-9 0 1e-7 800 300000 100 0.01`
//!
//! Set vacuum area to zero for filling. With equal areas, source T=1600 K and
//! p=400000 Pa give initial zero net mass but nonzero heating and recoil.
//! The source uses the SAME original gas catalog and remains INSTANTANEOUSLY
//! COMOVING with the vessel. It is prescribed, not a depleted finite donor or a
//! stationary inertial ambient. Both distinct aperture groups point +x; exterior
//! wall pressure/support forces and gravity are absent. Both opposing molecular
//! beams push the vessel toward -x. No second thermalization heat is added.
//!
//! GAS.tsv is the existing SYNTHETIC seven-column source format: species, mole
//! fraction, molar mass [kg/mol], Cp/R constant, Cp/R slope [1/K], viscosity
//! [Pa s], conductivity [W/(m K)] at 200 K. No measured accuracy is claimed.
//! Initial volume is 0.001 m3. The authored mean-free-path bound must cover every
//! evolving mixture over 200..1800 K and 1000..1e7 Pa; unit angle-independent,
//! nondeflecting transmission is assumed. Fixed steps remain first order and
//! require time refinement. A refused interval emits no new endpoint row.
//! CSV gross streams remain separate; authoritative net energy is evaluated
//! before subtracting gross formation-energy reports. Totals are compensated.

#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
use fs_couple::rocket_blowdown::{
    GasRocketState,
    effusion::{
        EffusionOutlet, EffusionRocketEnvironment, EffusionRocketSolve,
        exchange::{EffusionExchangeStep, EffusionExchangeStream, EffusionReservoirPort},
    },
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        momentum::EffusiveMomentumInterface,
    },
    network::{SpeciesNetworkSolve, transport::GasTransportSolve},
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_mbd::Vec3;
use std::io::{Read, Write};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "usage: effusion_feed GAS.tsv dry_mass_kg vacuum_area_m2 source_area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m [source_temperature_k=800] [source_pressure_pa=300000] [steps=100] [dt_s=0.01] [initial_temperature_k=400] [initial_pressure_pa=100000]";
const METRICS: [&str; 9] = [
    "mass_kg",
    "relative_momentum_x_kg_m_s",
    "relative_momentum_y_kg_m_s",
    "relative_momentum_z_kg_m_s",
    "relative_energy_j",
    "momentum_world_x_kg_m_s",
    "momentum_world_y_kg_m_s",
    "momentum_world_z_kg_m_s",
    "energy_world_j",
];

fn stream(row: &EffusionExchangeStream) -> [f64; 9] {
    let r = row.relative_momentum_kg_m_s;
    let p = row.momentum_world_kg_m_s;
    [
        row.mass_kg,
        r.x,
        r.y,
        r.z,
        row.relative_energy_j,
        p.x,
        p.y,
        p.z,
        row.energy_world_j,
    ]
}

fn accounts(step: &EffusionExchangeStep) -> [f64; 27] {
    let mut values = [0.0; 27];
    values[..9].copy_from_slice(&stream(&step.incoming));
    values[9..18].copy_from_slice(&stream(&step.outgoing));
    let r = step
        .incoming
        .relative_momentum_kg_m_s
        .sub(step.outgoing.relative_momentum_kg_m_s);
    let p = step.net_momentum_into_vessel_world_kg_m_s;
    values[18..].copy_from_slice(&[
        step.net_mass_into_gas_kg,
        r.x,
        r.y,
        r.z,
        step.net_relative_energy_into_gas_j,
        p.x,
        p.y,
        p.z,
        step.net_energy_into_vessel_world_j,
    ]);
    values
}

fn accumulate(sum: &mut [f64; 27], correction: &mut [f64; 27], terms: [f64; 27]) -> Result<()> {
    for k in 0..27 {
        let increment = terms[k] - correction[k];
        let next = sum[k] + increment;
        correction[k] = (next - sum[k]) - increment;
        sum[k] = next;
        if !next.is_finite() || !correction[k].is_finite() {
            return Err("cumulative particle account is not representable".into());
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(14).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !(7..=13).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let dry_mass: f64 = args[1].parse()?;
    let vacuum_area: f64 = args[2].parse()?;
    let source_area: f64 = args[3].parse()?;
    let diameter: f64 = args[4].parse()?;
    let thickness: f64 = args[5].parse()?;
    let mean_free_path: f64 = args[6].parse()?;
    let source_t: f64 = args.get(7).map_or(Ok(800.0), |v| v.parse())?;
    let source_p: f64 = args.get(8).map_or(Ok(300_000.0), |v| v.parse())?;
    let steps: u32 = args.get(9).map_or(Ok(100), |v| v.parse())?;
    let dt: f64 = args.get(10).map_or(Ok(0.01), |v| v.parse())?;
    let initial_t: f64 = args.get(11).map_or(Ok(400.0), |v| v.parse())?;
    let initial_p: f64 = args.get(12).map_or(Ok(100_000.0), |v| v.parse())?;
    if !dry_mass.is_finite()
        || dry_mass <= 0.0
        || [vacuum_area, source_area]
            .iter()
            .any(|a| !a.is_finite() || *a < 0.0)
        || !(1..=10_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || dt * f64::from(steps) > 3600.0
        || [source_t, initial_t]
            .iter()
            .any(|t| !(200.0..=1800.0).contains(t))
        || [source_p, initial_p]
            .iter()
            .any(|p| !(1000.0..=1e7).contains(p))
    {
        return Err("finite positive dry mass/dt, nonnegative areas, 1..=10000 steps, at most 3600 s, T in 200..1800 K and p in 1000..1e7 Pa required".into());
    }
    let mut text = String::new();
    std::fs::File::open(&args[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let model = source::load(&text, &args[0])?; // Original 64 KiB / 32-species bounds.
    let transmissions = model
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                1.0,
                "authored reciprocal thin aperture independent of molecular energy",
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let interface = EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            model.clone(),
            transmissions,
            EffusiveSupport {
                temperature_k: (200.0, 1800.0),
                pressure_pa: (1000.0, 1e7),
            },
            EffusiveRegime::try_new(
                diameter,
                thickness,
                mean_free_path,
                format!(
                    "authored lambda_min={mean_free_path:.17e} m for all evolving source mixtures at 200..1800 K and 1000..1e7 Pa; no empirical validation"
                ),
            )?,
        )?,
        "authored uncollimated angle-independent nondeflecting transmission; no wall scattering or measured validation",
    )?;
    // This source's volume only constructs its thermochemical state. The API
    // explicitly prescribes its distribution instead of spending its inventory.
    let supply = ClosedGasParcel::new(model.clone(), source::point(source_t, source_p)?, 0.001)?;
    let gas = ClosedGasParcel::new(model, source::point(initial_t, initial_p)?, 0.001)?;
    let mut vessel = GasRocketState::new(gas, dry_mass, Vec3::ZERO, Vec3::ZERO, 0)?;
    let vacuum = [EffusionOutlet {
        id: 0,
        interface: &interface,
        area_m2: vacuum_area,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }];
    let ports = [EffusionReservoirPort {
        id: 1,
        source: &supply,
        interface: &interface,
        area_m2: source_area,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }];
    let environment = EffusionRocketEnvironment {
        inertial_frame: 0,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
    };
    let endpoint = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let controls = EffusionRocketSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.1,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-8,
            },
            entropy_tolerance_j_k: 1e-10,
        },
        momentum_tolerance_n_s: 1e-12,
    };
    eprintln!(
        "SYNTHETIC source {}; finite 0.001 m3 gas at {initial_t} K/{initial_p} Pa; prescribed COMOVING source at {source_t} K/{source_p} Pa.",
        args[0]
    );
    eprintln!(
        "Distinct +x vacuum/source areas {vacuum_area}/{source_area} m2; no external wall force or gravity; fixed first-order dt={dt} s."
    );
    eprintln!(
        "Qualification: {}",
        interface.interface().regime().provenance()
    );
    eprintln!("Angular qualification: {}", interface.provenance());
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    write!(
        output,
        "time_s,temperature_k,pressure_pa,gas_mass_kg,position_x_m,position_y_m,position_z_m,velocity_x_m_s,velocity_y_m_s,velocity_z_m_s"
    )?;
    for window in ["step", "total"] {
        for direction in ["incoming", "outgoing", "net_into_vessel"] {
            for metric in METRICS {
                write!(output, ",{window}_{direction}_{metric}")?;
            }
        }
    }
    write!(
        output,
        ",step_mass_residual_kg,step_momentum_residual_x_kg_m_s,step_momentum_residual_y_kg_m_s,step_momentum_residual_z_kg_m_s,step_energy_residual_j"
    )?;
    let species = interface
        .interface()
        .gas_catalog()
        .thermodynamics()
        .components();
    for row in species {
        write!(
            output,
            ",\"mole_fraction_{}\"",
            row.species().as_str().replace('"', "\"\"")
        )?;
    }
    writeln!(output)?;
    let (mut totals, mut correction) = ([0.0; 27], [0.0; 27]);
    for i in 0..=steps {
        let (mut exchanged, mut residual) = ([0.0; 27], [0.0; 5]);
        if i > 0 {
            let candidate = vessel
                .advance_effusion_exchange(
                    environment,
                    &vacuum,
                    &ports,
                    dt,
                    endpoint,
                    controls,
                    || false,
                )
                .map_err(|e| {
                    format!(
                        "joint interval {i} refused at t={} s: {e}",
                        f64::from(i - 1) * dt
                    )
                })?;
            let (next, step) = candidate.into_parts();
            exchanged = accounts(&step);
            accumulate(&mut totals, &mut correction, exchanged)?;
            let p = step.momentum_residual_kg_m_s;
            residual = [
                step.gas.exchange.total_mass_residual_kg,
                p.x,
                p.y,
                p.z,
                step.energy_residual_j,
            ];
            vessel = next;
        }
        let gas = vessel.reservoir();
        let state = gas.state().state();
        let x = vessel.vehicle().position_m;
        let v = vessel.vehicle().velocity_m_s;
        write!(output, "{:.12e}", f64::from(i) * dt)?;
        for value in [
            state.temperature,
            state.pressure,
            gas.mass_kg(),
            x.x,
            x.y,
            x.z,
            v.x,
            v.y,
            v.z,
        ]
        .into_iter()
        .chain(exchanged)
        .chain(totals)
        .chain(residual)
        {
            write!(output, ",{value:.12e}")?;
        }
        let active = gas.model().thermodynamics();
        for row in species {
            let fraction = active
                .components()
                .binary_search_by(|s| s.species().cmp(row.species()))
                .map(|k| active.mole_composition().fractions()[k])
                .unwrap_or(0.0);
            write!(output, ",{fraction:.12e}")?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
