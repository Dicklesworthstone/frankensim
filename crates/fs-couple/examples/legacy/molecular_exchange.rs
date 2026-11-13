//! Two finite homogeneous gases exchange particles through one shared aperture.
//!
//! Run from the repository root:
//! `cargo run -p fs-couple --example molecular_exchange --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv --adaptive`
//! Use `--help` for explicit SI inputs. Defaults are authored SYNTHETIC examples:
//! left=(400 K, 100000 Pa, 0.001 m³, [0,100,0] m/s),
//! right=(800 K, 200000 Pa, 0.001 m³, [0,-100,0] m/s), aperture at rest, normal +x,
//! area=1e-7 m², pore diameter=1e-9 m, wall thickness=0 m, qualified minimum
//! mean free path=1e-7 m, maximum relative speed=2000 m/s, 100 intervals of 0.01 s.
//!
//! GAS.tsv uses the existing seven-column synthetic NASA/transport source format.
//! Both finite inventories use that SAME original catalog. Unit, angle/energy-
//! independent nondeflecting transmission is explicitly authored over every
//! evolving mixture at 200..1800 K, 1000..1e7 Pa and the selected speed bound.
//! Mean-free-path qualification must cover grazing paths and counterstreams;
//! these constants are not measured material or aperture validation.
//!
//! The zero-storage aperture transfers the same particles between fixed-volume
//! homogeneous gas owners. Their bulk velocity is gas momentum, not container
//! motion. Other-wall pressure/support forces, moving-volume work and particle
//! flight between separated vessels are not inferred. Equal counterflow can
//! transfer momentum with zero net mass. All finite steps remain first order;
//! optional refinement controls local error, not a global accuracy certificate
//! or a continuous source-domain enclosure. Refused intervals emit no endpoint.
//! CSV interval and cumulative accounts retain both gross directions separately.

#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
use fs_couple::{
    liquid_gas::evaporation::momentum::{MovingGasParcel, boundary::MovingGasBoundarySolve},
    molecular_exchange::{
        MolecularExchangeControls, MolecularExchangeSurface, step_candidate,
        trajectory::{MolecularExchangeAdaptiveControls, advance_candidate},
    },
};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::species::{
            boundary::GasSpeciesBoundarySolve,
            effusion::{
                EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
                drift::EffusiveDriftInterface, momentum::EffusiveMomentumInterface,
            },
        },
    },
};
use fs_math::dd::Dd;
use std::{
    collections::BTreeSet,
    io::{Read, Write},
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "usage: molecular_exchange GAS.tsv [--left T_K P_PA V_M3 VX VY VZ] [--right T_K P_PA V_M3 VX VY VZ] [--aperture-velocity VX VY VZ] [--normal NX NY NZ] [--area M2] [--geometry DIAMETER_M WALL_M MIN_MEAN_FREE_PATH_M] [--maximum-relative-speed M_PER_S] [--steps N] [--dt SECONDS] [--adaptive]\nVelocities are m/s in common inertial axes; normal is unit left-to-right. Defaults: left 400 100000 0.001 0 100 0; right 800 200000 0.001 0 -100 0; aperture velocity 0 0 0; normal 1 0 0; area 1e-7; geometry 1e-9 0 1e-7; maximum-relative-speed 2000; steps 100; dt 0.01. Authored synthetic qualifications, no empirical validation.";
const METRICS: [&str; 13] = [
    "gross_left_to_right_mass_kg",
    "gross_right_to_left_mass_kg",
    "net_left_to_right_mass_kg",
    "net_left_to_right_momentum_x_kg_m_s",
    "net_left_to_right_momentum_y_kg_m_s",
    "net_left_to_right_momentum_z_kg_m_s",
    "mass_residual_kg",
    "momentum_residual_x_kg_m_s",
    "momentum_residual_y_kg_m_s",
    "momentum_residual_z_kg_m_s",
    "stream_energy_residual_j",
    "aperture_energy_residual_j",
    "laboratory_energy_residual_j",
];

struct Input {
    file: String,
    gases: [[f64; 6]; 2],
    surface: MolecularExchangeSurface,
    geometry: [f64; 3],
    speed: f64,
    steps: u32,
    dt: f64,
    adaptive: bool,
}

fn numbers<const N: usize>(args: &[String], index: usize) -> Result<[f64; N]> {
    let fields = args
        .get(index + 1..index + 1 + N)
        .ok_or_else(|| format!("{} requires {N} numeric values", args[index]))?;
    let mut values = [0.0_f64; N];
    for (value, text) in values.iter_mut().zip(fields) {
        *value = text.parse()?;
        if !value.is_finite() {
            return Err(format!("{} requires finite SI values", args[index]).into());
        }
    }
    Ok(values)
}

fn input() -> Result<Option<Input>> {
    let args: Vec<_> = std::env::args().skip(1).take(37).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(None);
    }
    if args.is_empty()
        || args.len() > 36
        || args.iter().any(|s| s.len() > 1024)
        || args[0].starts_with("--")
    {
        return Err(USAGE.into());
    }
    let mut value = Input {
        file: args[0].clone(),
        gases: [
            [400.0, 100000.0, 0.001, 0.0, 100.0, 0.0],
            [800.0, 200000.0, 0.001, 0.0, -100.0, 0.0],
        ],
        surface: MolecularExchangeSurface {
            area_m2: 1e-7,
            normal_left_to_right: [1.0, 0.0, 0.0],
            velocity_world_m_s: [0.0; 3],
        },
        geometry: [1e-9, 0.0, 1e-7],
        speed: 2000.0,
        steps: 100,
        dt: 0.01,
        adaptive: false,
    };
    let (mut index, mut seen) = (1, BTreeSet::new());
    while index < args.len() {
        let flag = args[index].as_str();
        if !seen.insert(flag) {
            return Err(format!("duplicate option {flag}").into());
        }
        index += match flag {
            "--left" => {
                value.gases[0] = numbers(&args, index)?;
                7
            }
            "--right" => {
                value.gases[1] = numbers(&args, index)?;
                7
            }
            "--aperture-velocity" => {
                value.surface.velocity_world_m_s = numbers(&args, index)?;
                4
            }
            "--normal" => {
                value.surface.normal_left_to_right = numbers(&args, index)?;
                4
            }
            "--area" => {
                value.surface.area_m2 = numbers::<1>(&args, index)?[0];
                2
            }
            "--geometry" => {
                value.geometry = numbers(&args, index)?;
                4
            }
            "--maximum-relative-speed" => {
                value.speed = numbers::<1>(&args, index)?[0];
                2
            }
            "--dt" => {
                value.dt = numbers::<1>(&args, index)?[0];
                2
            }
            "--steps" => {
                value.steps = args
                    .get(index + 1)
                    .ok_or("--steps requires an integer")?
                    .parse()?;
                2
            }
            "--adaptive" => {
                value.adaptive = true;
                1
            }
            _ => return Err(format!("unknown option {flag}; {USAGE}").into()),
        };
    }
    let normal_squared: f64 = value
        .surface
        .normal_left_to_right
        .iter()
        .map(|v| v * v)
        .sum();
    if value.gases.iter().any(|g| {
        !(200.0..=1800.0).contains(&g[0]) || !(1000.0..=1e7).contains(&g[1]) || g[2] <= 0.0
    }) || value.surface.area_m2 < 0.0
        || value.speed < 0.0
        || value.dt <= 0.0
        || !(1..=10_000).contains(&value.steps)
        || value.dt * f64::from(value.steps) > 3600.0
        || !normal_squared.is_finite()
        || (normal_squared - 1.0).abs() > 64.0 * f64::EPSILON
        || (value.adaptive && value.dt * 2.0_f64.powi(-20) == 0.0)
    {
        return Err("require T in 200..1800 K, p in 1000..1e7 Pa, positive volumes/dt, nonnegative area/speed, unit normal, 1..=10000 intervals and at most 3600 s".into());
    }
    Ok(Some(value))
}

fn accounts(
    gross: [f64; 2],
    mass: f64,
    impulse: [f64; 3],
    mass_residual: f64,
    momentum_residual: [f64; 3],
    energy_residual: [f64; 3],
) -> [f64; 13] {
    [
        gross[0],
        gross[1],
        mass,
        impulse[0],
        impulse[1],
        impulse[2],
        mass_residual,
        momentum_residual[0],
        momentum_residual[1],
        momentum_residual[2],
        energy_residual[0],
        energy_residual[1],
        energy_residual[2],
    ]
}

fn main() -> Result<()> {
    let Some(input) = input()? else {
        return Ok(());
    };
    let mut text = String::new();
    std::fs::File::open(&input.file)?
        .take(65_537)
        .read_to_string(&mut text)?;
    let model = source::load(&text, &input.file)?;
    let transmission = model
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                1.0,
                "authored energy-independent unit transmission; synthetic thin-aperture example",
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let [diameter, thickness, mean_free_path] = input.geometry;
    let interface = EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(
                model.clone(),
                transmission,
                EffusiveSupport {
                    temperature_k: (200.0, 1800.0),
                    pressure_pa: (1000.0, 1e7),
                },
                EffusiveRegime::try_new(
                    diameter,
                    thickness,
                    mean_free_path,
                    "authored synthetic mean-free-path bound over all retained mixtures at 200..1800 K and 1000..1e7 Pa; no empirical calibration",
                )?,
            )?,
            "authored angle-independent nondeflecting molecular transmission",
        )?,
        input.speed,
        format!(
            "authored collisionless Maxwellian incidence including grazing and counterstream paths at relative speeds <= {} m/s; no empirical validation",
            input.speed
        ),
    )?;
    let make_gas = |g: [f64; 6]| -> Result<MovingGasParcel> {
        Ok(MovingGasParcel::new(
            ClosedGasParcel::new(model.clone(), source::point(g[0], g[1])?, g[2])?,
            [g[3], g[4], g[5]],
        )?)
    };
    let mut gases = [make_gas(input.gases[0])?, make_gas(input.gases[1])?];
    let boundary = MovingGasBoundarySolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 1800.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        momentum_absolute_n_s: 1e-12,
        momentum_relative_tolerance: 1e-12,
    };
    let controls = MolecularExchangeControls {
        gas: [boundary; 2],
        maximum_duration_s: input.dt,
        maximum_species_fraction: 0.05,
        mass_absolute_kg: 1e-13,
        momentum_absolute_n_s: 1e-12,
        energy_absolute_j: 5e-8,
        relative_tolerance: 1e-10,
    };
    let adaptive = MolecularExchangeAdaptiveControls {
        initial_step_s: input.dt,
        minimum_step_s: input.dt * 2.0_f64.powi(-20),
        maximum_step_s: input.dt,
        maximum_attempts: 512,
        maximum_accepted_intervals: 128,
        relative_tolerance: 1e-5,
        temperature_tolerance_k: 1e-4,
        pressure_tolerance_pa: 0.025,
        species_mass_tolerance_kg: 1e-13,
        velocity_tolerance_m_s: 1e-7,
    };
    eprintln!(
        "SYNTHETIC source {}; left/right [T_K,p_Pa,V_m3,vx,vy,vz]: {:?}; shared aperture {:?}.",
        input.file, input.gases, input.surface
    );
    eprintln!(
        "Authored diameter/wall/minimum-mean-free-path {:?} m; maximum relative speed {} m/s; fixed gas volumes, no other-wall pressure or particle-flight model.",
        input.geometry, input.speed
    );
    eprintln!(
        "First-order {} intervals of {} s; optional local-error refinement={}, no global accuracy or continuous-domain certificate.",
        input.steps, input.dt, input.adaptive
    );
    eprintln!(
        "Closed-vessel equilibration and entropy/passivity require the remaining wall/face momentum operators; reducing dt does not supply them. This example advances only the shared aperture exchange."
    );
    if input.adaptive {
        eprintln!("Adaptive local controls: {adaptive:?}");
    }

    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    write!(output, "time_s")?;
    for side in ["left", "right"] {
        for field in [
            "temperature_k",
            "pressure_pa",
            "volume_m3",
            "mass_kg",
            "velocity_x_m_s",
            "velocity_y_m_s",
            "velocity_z_m_s",
        ] {
            write!(output, ",{side}_{field}")?;
        }
    }
    for prefix in ["interval", "cumulative"] {
        for metric in METRICS {
            write!(output, ",{prefix}_{metric}")?;
        }
    }
    for side in ["left", "right"] {
        for species in model.thermodynamics().components() {
            write!(
                output,
                ",\"{side}_mole_fraction_{}\"",
                species.species().as_str().replace('"', "\"\"")
            )?;
        }
    }
    writeln!(
        output,
        ",accepted_native_substeps,attempted_intervals,rejected_intervals"
    )?;
    let mut totals = [Dd::ZERO; 13];
    for index in 0..=input.steps {
        let (mut interval, mut counts) = ([0.0; 13], [0_usize; 3]);
        if index > 0 {
            let fail = |e: &dyn std::fmt::Display| {
                format!(
                    "interval {index} refused at t={} s: {e}",
                    f64::from(index - 1) * input.dt
                )
            };
            let (next, values, accepted) = if input.adaptive {
                let (next, report) = advance_candidate(
                    [&gases[0], &gases[1]],
                    &interface,
                    input.surface,
                    input.dt,
                    controls,
                    adaptive,
                    || false,
                )
                .map_err(|e| fail(&e))?
                .into_parts();
                let t = &report.totals;
                (
                    next,
                    accounts(
                        t.gross_mass_kg,
                        t.net_mass_left_to_right_kg,
                        t.net_momentum_left_to_right_n_s,
                        t.mass_residual_kg,
                        t.momentum_residual_n_s,
                        [
                            t.stream_energy_residual_j,
                            t.energy_residual_j,
                            t.laboratory_energy_residual_j,
                        ],
                    ),
                    [
                        report.accepted_steps.len(),
                        report.attempted_intervals,
                        report.rejected_intervals,
                    ],
                )
            } else {
                let (next, report) = step_candidate(
                    [&gases[0], &gases[1]],
                    &interface,
                    input.surface,
                    input.dt,
                    controls,
                    || false,
                )
                .map_err(|e| fail(&e))?
                .into_parts();
                let mut gross = [Dd::ZERO; 2];
                for species in &report.exchange.species {
                    for side in 0..2 {
                        gross[side] = gross[side] + Dd::from_f64(species.mass_kg[side]);
                    }
                }
                (
                    next,
                    accounts(
                        gross.map(|m| m.to_f64()),
                        report.net_mass_left_to_right_kg,
                        report.net_momentum_left_to_right_n_s,
                        report.mass_residual_kg,
                        report.momentum_residual_n_s,
                        [
                            report.stream_energy_residual_j,
                            report.energy_residual_j,
                            report.laboratory_energy_residual_j,
                        ],
                    ),
                    [1, 1, 0],
                )
            };
            for (total, value) in totals.iter_mut().zip(values) {
                *total = *total + Dd::from_f64(value);
                if !total.to_f64().is_finite() {
                    return Err("cumulative particle account is unrepresentable".into());
                }
            }
            gases = next;
            interval = values;
            counts = accepted;
        }
        write!(output, "{:.12e}", f64::from(index) * input.dt)?;
        for gas in &gases {
            let state = gas.parcel().state().state();
            let v = gas.velocity_m_s();
            for value in [
                state.temperature,
                state.pressure,
                gas.parcel().volume_m3(),
                gas.parcel().mass_kg(),
                v[0],
                v[1],
                v[2],
            ] {
                write!(output, ",{value:.12e}")?;
            }
        }
        for value in interval
            .into_iter()
            .chain(totals.iter().map(|x| x.to_f64()))
        {
            write!(output, ",{value:.12e}")?;
        }
        for gas in &gases {
            let mixture = gas.parcel().model().thermodynamics();
            for species in model.thermodynamics().components() {
                let fraction = mixture
                    .components()
                    .iter()
                    .position(|s| s.species() == species.species())
                    .map_or(0.0, |i| mixture.mole_composition().fractions()[i]);
                write!(output, ",{fraction:.12e}")?;
            }
        }
        writeln!(output, ",{},{},{}", counts[0], counts[1], counts[2])?;
    }
    output.flush()?;
    Ok(())
}
