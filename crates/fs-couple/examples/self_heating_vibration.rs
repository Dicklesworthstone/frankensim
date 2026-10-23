//! A source-card specimen rings down, heats, and changes its own relaxation.
//!
//! Uses the existing creep example's bounded material-card/table ingress. The
//! physical solver is the general inertial material owner, not a fitted audio
//! oscillator. Fixed reference geometry: L=0.1 m, A=1e-5 m2, square section.
//! All spring/relaxation/caloric/transport coefficients come from the input.
//!
//! Example (manufactured data, not a calibrated material):
//! cargo run -p fs-couple --example self_heating_vibration -- --table \
//! crates/fs-couple/examples/self_heating_vibration/material.txt 0.5 0 0.04 800 0
//!
//! Arguments: --table FILE VELOCITY PAYLOAD_KG DURATION STEPS FILM_COEFFICIENT
//! or: --pack FILE BRANCH_COUNT VELOCITY PAYLOAD_KG DURATION STEPS FILM_COEFFICIENT.
//! Output is the accepted displacement, velocity, temperature and energy history.
//! Add --finite-gas VOLUME_M3 TMIN_K TMAX_K MAX_EXCHANGE_J with the shared
//! --gas-record/--gas-component/--gas-water/--gas-phase source options to evolve
//! a finite gas reservoir. See self_heating_vibration/README.md for a complete run.
#[path = "self_heating_vibration/finite_gas.rs"]
mod finite_gas;
#[path = "self_heating_creep/source.rs"]
mod source;

use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_couple::thermovisco_dynamics::{Controls, Geometry, Oscillator};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::thermovisco::thermal::FeedbackConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const LENGTH_M: f64 = 0.1;
const AREA_M2: f64 = 1e-5;
const VOLUME_M3: f64 = LENGTH_M * AREA_M2;
fn surface_area() -> f64 {
    4.0 * AREA_M2.sqrt() * LENGTH_M + 2.0 * AREA_M2
}

fn run(args: &[String]) -> Result<()> {
    let (args, finite) = finite_gas::Options::extract(args)?;
    let mode = args.first().map(String::as_str).unwrap_or("");
    let offset = match (mode, args.len()) {
        ("--table", 7) => 0, ("--pack", 8) => 1,
        _ => return Err("usage: self_heating_vibration --table FILE VELOCITY PAYLOAD_KG DURATION STEPS FILM_COEFFICIENT; or --pack FILE BRANCH_COUNT VELOCITY PAYLOAD_KG DURATION STEPS FILM_COEFFICIENT".into()),
    };
    let velocity: f64 = args[2 + offset].parse()?;
    let payload: f64 = args[3 + offset].parse()?;
    let duration: f64 = args[4 + offset].parse()?;
    let count: usize = args[5 + offset].parse()?;
    let film: f64 = args[6 + offset].parse()?;
    if !velocity.is_finite()
        || !payload.is_finite()
        || payload < 0.0
        || !duration.is_finite()
        || duration <= 0.0
        || duration > 1000.0
        || !(1..=100_000).contains(&count)
        || !film.is_finite()
        || film < 0.0
    {
        return Err("finite velocity, nonnegative payload/film, duration in (0,1000] s and 1..100000 steps required".into());
    }
    if let Some(options) = finite {
        return finite_gas::run(
            mode,
            &args[1],
            if offset == 1 {
                Some(args[2].parse()?)
            } else {
                None
            },
            finite_gas::Motion {
                velocity,
                payload,
                duration,
                count,
                film,
            },
            options,
        );
    }
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let (description, records) = pool.scope(|arena| -> Result<_> {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 0, kernel_id: 0x544d44, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let input = if mode == "--pack" { source::Source::pack(&args[1], args[2].parse()?, &gate)? }
            else { source::Source::table(&args[1], &gate)? };
        let mut state = Oscillator::try_new(&input.mechanical, &input.caloric, Geometry {
            length_m: LENGTH_M, cross_section_m2: AREA_M2, thermal_surface_m2: surface_area(), end_payload_kg: payload,
        }, 300.0, 0.0, velocity)?;
        let description = format!("material_card={}, caloric_curve={}, branches={}, source_states={}; L={LENGTH_M} m, A={AREA_M2} m2; specimen mass={} kg; uniform-mode inertia=specimen_mass/3+payload={} kg; initial velocity={velocity} m/s; film={film} W/(m2 K); ambient=300 K; radiation disabled. Fixed-geometry, small-strain, single-phase axial reduction; temperature changes relaxation. Source provenance does not establish measured material or vibration accuracy.",
            input.card.content_hash().to_hex(), input.caloric.curve().identity().to_hex(),
            input.mechanical.law().branch_count(), input.mechanical.samples().len(),
            state.point().mass_kg(), state.effective_mass_kg());
        let environment = LumpedThermalEnvironment {
            temperature_k: 300.0, radiation_temperature_k: 300.0, convection_w_per_m2_k: film,
            transport: input.transport.clone(), enthalpy_tolerance_j_kg: 1e-12, maximum_thermal_residual_j: 1e-10,
        };
        let controls = Controls { feedback: FeedbackConfig { max_iterations: 100,
            temperature_tolerance_k: 1e-8, energy_tolerance_j: 1e-9, relaxation: 1.0 },
            force_tolerance_n: 1e-6, energy_tolerance_j: 2e-9 };
        let initial_energy = state.kinetic_energy_j();
        let initial_enthalpy = state.point().caloric_state().specific_enthalpy_j_kg();
        let (mut numerical, mut heat, mut boundary_heat) = (0.0, 0.0, 0.0);
        let mut records = Vec::<[f64; 9]>::with_capacity(count);
        for index in 1..=count {
            let target = duration*(index as f64/count as f64);
            let step = state.step(0.0, target-state.point().time_s(), controls,
                |input| environment.advance(&cx, input), &gate)?;
            numerical += step.numerical_loss_j;
            heat += step.material_heat_j;
            boundary_heat += step.transport.external_heat_j;
            state = step.next;
            let stored = state.kinetic_energy_j()+state.point().volume_m3()*state.point().mechanical().stored_energy_j_m3()
                + state.point().mass_kg()*(state.point().caloric_state().specific_enthalpy_j_kg()-initial_enthalpy);
            let defect = stored+numerical-initial_energy-boundary_heat;
            if !defect.is_finite() || defect.abs() > index as f64*controls.energy_tolerance_j {
                return Err("cumulative mechanical/thermal energy balance refused".into());
            }
            records.push([state.point().time_s(), state.extension_m(), state.velocity_m_s(),
                state.point().caloric_state().temperature_k(), heat, numerical, boundary_heat, stored, defect]);
        }
        Ok((description, records))
    })?;
    eprintln!("{description}");
    println!(
        "time_s,extension_m,velocity_m_s,temperature_K,physical_material_heat_J,numerical_loss_J,boundary_heat_J,stored_energy_J,cumulative_energy_defect_J"
    );
    for r in records {
        println!(
            "{},{},{},{},{},{},{},{},{}",
            r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]
        );
    }
    Ok(())
}
fn main() -> Result<()> {
    run(&std::env::args().skip(1).collect::<Vec<_>>())
}
