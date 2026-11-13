//! Finite source gas -> molecular vacuum escape -> vessel recoil and exhaust.
//!
//! Run from the repository root:
//! `cargo run -p fs-couple --example effusion_recoil --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv
//! 0.01 1e-7 1e-9 0 1e-7 100 0.01`
//! Append `adaptive 0,-9.81,0 0.002,0,0` for adaptive motion under gravity
//! and a 0.002 N force along +x.
//!
//! GAS.tsv uses the existing SYNTHETIC seven-column example format: species,
//! mole fraction, molar mass [kg/mol], Cp/R constant, Cp/R slope [1/K],
//! viscosity [Pa s], conductivity [W/(m K)] at 200 K. No measured accuracy is
//! claimed. Initial gas is 0.001 m3 at 400 K and 100000 Pa; the vessel starts at
//! rest, with a +x outlet. Optional held world gravity and external force are
//! comma-separated x,y,z vectors; their defaults are zero. Recoil points -x.
//! The authored minimum mean free path must cover every evolving mixture over
//! 200..1800 K and 1000..1e7 Pa. Unit, angle-independent, nondeflecting
//! transmission is assumed. Optional `adaptive` mode uses bounded step doubling
//! for cooling, species depletion and motion; the accepted method remains first
//! order. Its authored local tolerances are estimates, not global certificates.
//! Refused intervals emit no endpoint row. All mechanical accounts use the
//! vessel's declared inertial frame, including external impulse and work.

#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
use fs_couple::rocket_blowdown::{
    GasRocketState,
    effusion::{
        EffusionOutlet, EffusionRocketEnvironment, EffusionRocketSolve, EffusionRocketStep,
        trajectory::EffusionAdaptiveSolve,
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
const USAGE: &str = "usage: effusion_recoil GAS.tsv dry_mass_kg exposed_area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m [steps=100] [dt_s=0.01] [mode=fixed|adaptive] [gravity_world_m_s2=0,0,0] [external_force_world_n=0,0,0]";

fn vector_argument(value: Option<&String>) -> Result<Vec3> {
    let Some(value) = value else {
        return Ok(Vec3::ZERO);
    };
    let mut entries = value.split(',');
    let vector = Vec3::new(
        entries.next().ok_or("expected x,y,z vector")?.parse()?,
        entries.next().ok_or("expected x,y,z vector")?.parse()?,
        entries.next().ok_or("expected x,y,z vector")?.parse()?,
    );
    if entries.next().is_some() || !vector.is_finite() {
        return Err(
            "gravity and force each require three finite comma-separated components".into(),
        );
    }
    Ok(vector)
}

fn accumulate<const N: usize>(
    sum: &mut [f64; N],
    correction: &mut [f64; N],
    terms: [f64; N],
) -> Result<()> {
    for k in 0..N {
        let increment = terms[k] - correction[k];
        let next = sum[k] + increment;
        correction[k] = (next - sum[k]) - increment;
        sum[k] = next;
        if !next.is_finite() || !correction[k].is_finite() {
            return Err("cumulative physical account is not representable".into());
        }
    }
    Ok(())
}

fn accounts(step: &EffusionRocketStep) -> [f64; 9] {
    let p = step.exhaust_momentum_world_kg_m_s;
    let impulse = step.exhaust.external_impulse_n_s;
    [
        step.expelled_mass_kg,
        p.x,
        p.y,
        p.z,
        step.exhaust_energy_world_j,
        impulse.x,
        impulse.y,
        impulse.z,
        step.exhaust.external_work_j,
    ]
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(12).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !(6..=11).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let dry_mass: f64 = args[1].parse()?;
    let area: f64 = args[2].parse()?;
    let diameter: f64 = args[3].parse()?;
    let thickness: f64 = args[4].parse()?;
    let mean_free_path: f64 = args[5].parse()?;
    let steps: u32 = args.get(6).map_or(Ok(100), |v| v.parse())?;
    let dt: f64 = args.get(7).map_or(Ok(0.01), |v| v.parse())?;
    let adaptive_mode = match args.get(8).map_or("fixed", String::as_str) {
        "fixed" => false,
        "adaptive" => true,
        _ => return Err("mode must be fixed or adaptive".into()),
    };
    let environment = EffusionRocketEnvironment {
        inertial_frame: 0,
        gravity_world_m_s2: vector_argument(args.get(9))?,
        external_force_world_n: vector_argument(args.get(10))?,
    };
    if !dry_mass.is_finite()
        || dry_mass <= 0.0
        || !area.is_finite()
        || area < 0.0
        || !(1..=10_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || dt * f64::from(steps) > 3600.0
    {
        return Err("positive finite dry mass, nonnegative finite area, 1..=10000 steps, positive finite dt and at most 3600 s required".into());
    }
    let mut text = String::new();
    std::fs::File::open(&args[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let model = source::load(&text, &args[0])?; // Original loader enforces 64 KiB / 32 species.
    let transmissions = model
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                1.0,
                "authored open thin apertures; reciprocal and independent of molecular energy",
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let qualification = format!(
        "authored CLI bound lambda_min={mean_free_path:.17e} m over all source mixtures at 200..1800 K, 1000..1e7 Pa; no empirical validation"
    );
    let interface = EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            model.clone(),
            transmissions,
            EffusiveSupport {
                temperature_k: (200.0, 1800.0),
                pressure_pa: (1000.0, 1e7),
            },
            EffusiveRegime::try_new(diameter, thickness, mean_free_path, qualification)?,
        )?,
        "authored uncollimated angle-independent nondeflecting transmission; no wall scattering or measured validation",
    )?;
    let gas = ClosedGasParcel::new(model, source::point(400.0, 100_000.0)?, 0.001)?;
    let mut vessel = GasRocketState::new(gas, dry_mass, Vec3::ZERO, Vec3::ZERO, 0)?;
    let outlets = [EffusionOutlet {
        id: 0,
        interface: &interface,
        area_m2: area,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }];
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
        "SYNTHETIC source: {}; 0.001 m3, 400 K, 100000 Pa; +x vacuum emission, dry mass {dry_mass} kg.",
        args[0]
    );
    eprintln!(
        "Qualification: {}",
        interface.interface().regime().provenance()
    );
    eprintln!("Angular qualification: {}", interface.provenance());
    eprintln!("Held environment: {environment:?}");
    if adaptive_mode {
        eprintln!("First-order adaptive local error controls: {adaptive:?}");
    }
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    write!(
        output,
        concat!(
            "time_s,temperature_k,pressure_pa,gas_mass_kg,",
            "position_x_m,position_y_m,position_z_m,velocity_x_m_s,velocity_y_m_s,velocity_z_m_s,",
            "step_escaped_mass_kg,total_escaped_mass_kg,",
            "step_escaped_momentum_x_kg_m_s,step_escaped_momentum_y_kg_m_s,step_escaped_momentum_z_kg_m_s,",
            "total_escaped_momentum_x_kg_m_s,total_escaped_momentum_y_kg_m_s,total_escaped_momentum_z_kg_m_s,",
            "step_escaped_energy_j,total_escaped_energy_j,",
            "step_external_impulse_x_n_s,step_external_impulse_y_n_s,step_external_impulse_z_n_s,",
            "total_external_impulse_x_n_s,total_external_impulse_y_n_s,total_external_impulse_z_n_s,",
            "step_external_work_j,total_external_work_j,",
            "step_mass_residual_kg,step_momentum_residual_x_kg_m_s,step_momentum_residual_y_kg_m_s,step_momentum_residual_z_kg_m_s,step_energy_residual_j,",
            "step_gravity_impulse_reserve_x_n_s,step_gravity_impulse_reserve_y_n_s,step_gravity_impulse_reserve_z_n_s,step_gravity_work_reserve_j,",
            "accepted_substeps,attempted_intervals,rejected_intervals"
        )
    )?;
    for species in vessel.reservoir().model().thermodynamics().components() {
        write!(
            output,
            ",\"mole_fraction_{}\"",
            species.species().as_str().replace('"', "\"\"")
        )?;
    }
    writeln!(output)?;
    let (mut totals, mut compensation) = ([0.0; 9], [0.0; 9]);
    for i in 0..=steps {
        let (mut escaped, mut residuals) = ([0.0; 9], [0.0; 9]);
        let (mut escaped_correction, mut residual_correction) = ([0.0; 9], [0.0; 9]);
        let mut counts = [0_usize; 3];
        if i > 0 {
            let context = |error: &dyn std::fmt::Display| {
                format!(
                    "joint gas/vessel interval {i} refused at t={} s: {error}",
                    f64::from(i - 1) * dt
                )
            };
            let (next, accepted_steps) = if adaptive_mode {
                let candidate = vessel
                    .advance_effusion_adaptive(
                        environment,
                        &outlets,
                        dt,
                        endpoint,
                        controls,
                        adaptive,
                        || false,
                    )
                    .map_err(|error| context(&error))?;
                let (next, report) = candidate.into_parts();
                counts = [
                    report.accepted_steps.len(),
                    report.attempted_intervals,
                    report.rejected_intervals,
                ];
                (next, report.accepted_steps)
            } else {
                let candidate = vessel
                    .advance_effusion_with_environment(
                        environment,
                        &outlets,
                        dt,
                        endpoint,
                        controls,
                        || false,
                    )
                    .map_err(|error| context(&error))?;
                let (next, step) = candidate.into_parts();
                counts = [1, 1, 0];
                (next, vec![step])
            };
            for step in &accepted_steps {
                accumulate(&mut escaped, &mut escaped_correction, accounts(step))?;
                let p = step.momentum_residual_kg_m_s;
                let reserve = step.gravity_impulse_uncertainty_n_s;
                accumulate(
                    &mut residuals,
                    &mut residual_correction,
                    [
                        step.gas.exchange.total_mass_residual_kg,
                        p.x,
                        p.y,
                        p.z,
                        step.energy_residual_j,
                        reserve.x,
                        reserve.y,
                        reserve.z,
                        step.gravity_work_uncertainty_j,
                    ],
                )?;
            }
            accumulate(&mut totals, &mut compensation, escaped)?;
            vessel = next;
        }
        let gas = vessel.reservoir();
        let state = gas.state().state();
        let row = [
            f64::from(i) * dt,
            state.temperature,
            state.pressure,
            gas.mass_kg(),
            vessel.vehicle().position_m.x,
            vessel.vehicle().position_m.y,
            vessel.vehicle().position_m.z,
            vessel.vehicle().velocity_m_s.x,
            vessel.vehicle().velocity_m_s.y,
            vessel.vehicle().velocity_m_s.z,
            escaped[0],
            totals[0],
            escaped[1],
            escaped[2],
            escaped[3],
            totals[1],
            totals[2],
            totals[3],
            escaped[4],
            totals[4],
            escaped[5],
            escaped[6],
            escaped[7],
            totals[5],
            totals[6],
            totals[7],
            escaped[8],
            totals[8],
            residuals[0],
            residuals[1],
            residuals[2],
            residuals[3],
            residuals[4],
            residuals[5],
            residuals[6],
            residuals[7],
            residuals[8],
        ];
        write!(output, "{:.12e}", row[0])?;
        for value in &row[1..] {
            write!(output, ",{value:.12e}")?;
        }
        for count in counts {
            write!(output, ",{count}")?;
        }
        for fraction in gas.model().thermodynamics().mole_composition().fractions() {
            write!(output, ",{fraction:.12e}")?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
