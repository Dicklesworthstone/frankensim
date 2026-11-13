//! Finite vessel filling and discharge through distinct qualified molecular ports.
//!
//! Run from the repository root:
//! `cargo run -p fs-couple --example effusion_feed --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv
//! 0.01 1e-7 1e-7 1e-9 0 1e-7 800 300000 100 0.01`
//! Append `--inertial-source 0 0 0 2000 --vessel-velocity 150 30 0 --adaptive`
//! for a stationary source, an initially moving vessel and bounded refinement.
//!
//! Set vacuum area to zero for filling. With equal areas, source T=1600 K and
//! p=400000 Pa give initial zero net mass but nonzero heating and recoil.
//! The source uses the SAME original gas catalog and is COMOVING by default.
//! `--inertial-source` instead fixes its bulk velocity in frame 0 and explicitly
//! qualifies relative speeds up to MAX_RELATIVE_SPEED [m/s]. It is prescribed,
//! not a depleted finite donor. Both distinct aperture groups point +x; exterior
//! wall pressure/support forces and gravity are absent. Normal molecular moments
//! push toward -x; ambient slip also transfers tangential momentum and energy.
//! No second thermalization heat is added.
//!
//! GAS.tsv is the existing SYNTHETIC seven-column source format: species, mole
//! fraction, molar mass [kg/mol], Cp/R constant, Cp/R slope [1/K], viscosity
//! [Pa s], conductivity [W/(m K)] at 200 K. No measured accuracy is claimed.
//! Initial volume is 0.001 m3. The authored mean-free-path bound must cover every
//! evolving mixture over 200..1800 K and 1000..1e7 Pa; unit angle-independent,
//! nondeflecting transmission is assumed. Drift qualification additionally covers
//! grazing transmission and counterstream collisions at every admitted speed.
//! All finite steps remain first order. `--adaptive` refines each CSV interval
//! internally using authored local-error estimates, not a global certificate.
//! Speeds are admitted at sampled states and endpoints, without a continuous
//! trajectory enclosure. A refused interval emits no new endpoint row.
//! CSV gross streams remain separate; authoritative net energy is evaluated
//! before subtracting gross formation-energy reports. Totals are compensated.

#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
use fs_couple::rocket_blowdown::{
    GasRocketState,
    effusion::{
        EffusionOutlet, EffusionRocketEnvironment, EffusionRocketSolve,
        exchange::{EffusionAmbientPort, EffusionExchangeStream, EffusionReservoirPort},
        trajectory::EffusionAdaptiveSolve,
    },
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        drift::EffusiveDriftInterface, momentum::EffusiveMomentumInterface,
    },
    network::{SpeciesNetworkSolve, transport::GasTransportSolve},
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_mbd::Vec3;
use std::io::{Read, Write};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "usage: effusion_feed GAS.tsv dry_mass_kg vacuum_area_m2 source_area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m [source_temperature_k=800] [source_pressure_pa=300000] [steps=100] [dt_s=0.01] [initial_temperature_k=400] [initial_pressure_pa=100000] [--inertial-source VX VY VZ MAX_RELATIVE_SPEED] [--vessel-velocity VX VY VZ] [--adaptive]";
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

fn accounts(
    incoming: &EffusionExchangeStream,
    outgoing: &EffusionExchangeStream,
    net_mass: f64,
    net_relative_energy: f64,
    net_momentum: Vec3,
    net_energy: f64,
) -> [f64; 27] {
    let mut values = [0.0; 27];
    values[..9].copy_from_slice(&stream(incoming));
    values[9..18].copy_from_slice(&stream(outgoing));
    let r = incoming
        .relative_momentum_kg_m_s
        .sub(outgoing.relative_momentum_kg_m_s);
    let p = net_momentum;
    values[18..].copy_from_slice(&[
        net_mass,
        r.x,
        r.y,
        r.z,
        net_relative_energy,
        p.x,
        p.y,
        p.z,
        net_energy,
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

fn velocity(values: &[String], flag: &str) -> Result<Vec3> {
    let vector = Vec3::new(values[0].parse()?, values[1].parse()?, values[2].parse()?);
    if !vector.is_finite() {
        return Err(format!("{flag} requires three finite SI velocity components").into());
    }
    Ok(vector)
}

fn main() -> Result<()> {
    let raw: Vec<_> = std::env::args().skip(1).take(24).collect();
    if raw == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if raw.len() > 23 {
        return Err(USAGE.into());
    }
    let mut args = Vec::with_capacity(13);
    let (mut inertial_source, mut initial_velocity) = (None, None);
    let (mut adaptive_mode, mut index) = (false, 0);
    while index < raw.len() {
        match raw[index].as_str() {
            "--inertial-source" => {
                if inertial_source.is_some() {
                    return Err("duplicate --inertial-source".into());
                }
                let values = raw
                    .get(index + 1..index + 5)
                    .ok_or("--inertial-source requires VX VY VZ MAX_RELATIVE_SPEED")?;
                let speed: f64 = values[3].parse()?;
                if !speed.is_finite() || speed < 0.0 {
                    return Err("MAX_RELATIVE_SPEED must be finite and nonnegative in m/s".into());
                }
                inertial_source = Some((velocity(values, "--inertial-source")?, speed));
                index += 5;
            }
            "--vessel-velocity" => {
                if initial_velocity.is_some() {
                    return Err("duplicate --vessel-velocity".into());
                }
                let values = raw
                    .get(index + 1..index + 4)
                    .ok_or("--vessel-velocity requires VX VY VZ")?;
                initial_velocity = Some(velocity(values, "--vessel-velocity")?);
                index += 4;
            }
            "--adaptive" => {
                if adaptive_mode {
                    return Err("duplicate --adaptive".into());
                }
                adaptive_mode = true;
                index += 1;
            }
            value if value.starts_with("--") => {
                return Err(format!("unknown option {value}; {USAGE}").into());
            }
            value => {
                if args.len() == 13 {
                    return Err(USAGE.into());
                }
                args.push(value);
                index += 1;
            }
        }
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
    std::fs::File::open(args[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let model = source::load(&text, args[0])?; // Original 64 KiB / 32-species bounds.
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
    let drifting_interface = inertial_source.map(|(_, speed)| {
        EffusiveDriftInterface::try_new(interface.clone(), speed, format!(
            "authored collisionless Maxwellian incidence, angle/energy-independent nondeflecting transmission including grazing paths and counterstreams, for every retained mixture/T/p and relative speed <= {speed:.17e} m/s; no empirical validation"
        ))
    }).transpose()?;
    // This source's volume only constructs its thermochemical state. The API
    // explicitly prescribes its distribution instead of spending its inventory.
    let supply = ClosedGasParcel::new(model.clone(), source::point(source_t, source_p)?, 0.001)?;
    let gas = ClosedGasParcel::new(model, source::point(initial_t, initial_p)?, 0.001)?;
    let mut vessel = GasRocketState::new(
        gas,
        dry_mass,
        Vec3::ZERO,
        initial_velocity.unwrap_or(Vec3::ZERO),
        0,
    )?;
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
    let comoving_ports: &[EffusionReservoirPort<'_>] = if inertial_source.is_some() {
        &[]
    } else {
        &ports
    };
    let ambient =
        drifting_interface
            .as_ref()
            .zip(inertial_source)
            .map(|(law, (source_velocity, _))| EffusionAmbientPort {
                id: 1,
                source: &supply,
                velocity_world_m_s: source_velocity,
                inertial_frame: 0,
                interface: law,
                area_m2: source_area,
                outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
            });
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
    let adaptive = EffusionAdaptiveSolve {
        initial_step_s: dt,
        minimum_step_s: dt * 2.0_f64.powi(-20),
        maximum_step_s: dt,
        maximum_attempts: 512,
        maximum_accepted_intervals: 128,
        relative_tolerance: 1e-5,
        temperature_tolerance_k: 1e-4,
        pressure_tolerance_pa: 0.025,
        species_mass_tolerance_kg: 1e-13,
        position_tolerance_m: 1e-8,
        velocity_tolerance_m_s: 1e-7,
    };
    eprintln!(
        "SYNTHETIC source {}; finite 0.001 m3 gas at {initial_t} K/{initial_p} Pa; prescribed source at {source_t} K/{source_p} Pa.",
        args[0]
    );
    eprintln!(
        "Distinct +x vacuum/source areas {vacuum_area}/{source_area} m2; no external wall force or gravity; CSV interval={dt} s."
    );
    eprintln!(
        "Qualification: {}",
        interface.interface().regime().provenance()
    );
    eprintln!("Angular qualification: {}", interface.provenance());
    if let Some(source) = ambient {
        eprintln!(
            "FIXED INERTIAL source velocity {:?} m/s in frame 0; initial vessel velocity {:?} m/s.",
            source.velocity_world_m_s,
            vessel.vehicle().velocity_m_s
        );
        eprintln!("Drift qualification: {}", source.interface.provenance());
    } else {
        eprintln!(
            "Source is INSTANTANEOUSLY COMOVING; initial vessel velocity {:?} m/s.",
            vessel.vehicle().velocity_m_s
        );
    }
    if adaptive_mode {
        eprintln!(
            "Each CSV interval is refined internally; accepted method remains first order. Local estimates, not global accuracy bounds: {adaptive:?}"
        );
    } else {
        eprintln!(
            "Fixed first-order steps freeze sampled relative moments; refine dt to assess accuracy."
        );
    }
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
    writeln!(
        output,
        ",accepted_substeps,attempted_intervals,rejected_intervals"
    )?;
    let (mut totals, mut correction) = ([0.0; 27], [0.0; 27]);
    for i in 0..=steps {
        let (mut exchanged, mut residual) = ([0.0; 27], [0.0; 5]);
        let mut counts = [0_usize; 3];
        if i > 0 {
            let context = |error: &dyn std::fmt::Display| {
                format!(
                    "joint interval {i} refused at t={} s: {error}",
                    f64::from(i - 1) * dt
                )
            };
            let next = if adaptive_mode {
                let candidate = vessel
                    .advance_effusion_exchange_with_ambient_adaptive(
                        environment,
                        &vacuum,
                        comoving_ports,
                        ambient.as_slice(),
                        dt,
                        endpoint,
                        controls,
                        adaptive,
                        || false,
                    )
                    .map_err(|e| context(&e))?;
                let (next, step) = candidate.into_parts();
                exchanged = accounts(
                    &step.total_incoming,
                    &step.total_outgoing,
                    step.total_net_mass_into_gas_kg,
                    step.total_net_relative_energy_into_gas_j,
                    step.total_net_momentum_into_vessel_world_kg_m_s,
                    step.total_net_energy_into_vessel_world_j,
                );
                let p = step.total_momentum_residual_kg_m_s;
                residual = [
                    step.total_mass_residual_kg,
                    p.x,
                    p.y,
                    p.z,
                    step.total_energy_residual_j,
                ];
                counts = [
                    step.accepted_steps.len(),
                    step.attempted_intervals,
                    step.rejected_intervals,
                ];
                next
            } else {
                let candidate = vessel
                    .advance_effusion_exchange_with_ambient(
                        environment,
                        &vacuum,
                        comoving_ports,
                        ambient.as_slice(),
                        dt,
                        endpoint,
                        controls,
                        || false,
                    )
                    .map_err(|e| context(&e))?;
                let (next, step) = candidate.into_parts();
                exchanged = accounts(
                    &step.incoming,
                    &step.outgoing,
                    step.net_mass_into_gas_kg,
                    step.net_relative_energy_into_gas_j,
                    step.net_momentum_into_vessel_world_kg_m_s,
                    step.net_energy_into_vessel_world_j,
                );
                let p = step.momentum_residual_kg_m_s;
                residual = [
                    step.gas.exchange.total_mass_residual_kg,
                    p.x,
                    p.y,
                    p.z,
                    step.energy_residual_j,
                ];
                counts = [1, 1, 0];
                next
            };
            accumulate(&mut totals, &mut correction, exchanged)?;
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
        writeln!(output, ",{},{},{}", counts[0], counts[1], counts[2])?;
    }
    output.flush()?;
    Ok(())
}
