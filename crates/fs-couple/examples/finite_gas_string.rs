//! Source packs -> electrical/thermal feedback -> continuous string motion and sound.
//!
//! The existing string, conductor, caloric, finite-gas and conduction owners
//! advance one shared accepted trajectory. CSV contains physical pressure and
//! independent cumulative storage accounting; WAV has an explicit fixed scale.
//! See finite_gas_string/README.md for source compilation and the bounded input.

#[path = "finite_gas_string/input.rs"]
mod input;
#[allow(dead_code)]
#[path = "../src/bin/music_render/json.rs"]
mod json;
#[allow(dead_code)]
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
use json::json_string;

use fs_couple::acoustic_realize::{
    LinearMaterialStringRuntime, StringPointForce, ThermalMaterialStringRuntime,
    finite_gas::{FiniteGasStringHeat, FiniteGasStringSolve},
};
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::string_specimen::{
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, StringGeometryConstraint, StringPrestress,
    with_uniform_circular_material_and_constraints,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::gas_species::input::{GasRecord, read_isochoric_parcel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_qty::{Density, DynViscosity, Pressure};
use fs_scenario::PrestressedString;
use std::{
    io::{Read, Write},
    path::Path,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn bounded_read(path: &Path, limit: usize) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take((limit + 1) as u64)
        .read_to_end(&mut bytes)?;
    if bytes.len() > limit {
        return Err(format!("{} exceeds {limit} bytes", path.display()).into());
    }
    Ok(bytes)
}

fn material(pack: &NormalizedPack) -> Result<MaterialCard> {
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("source-pack:{}", pack.pack_id()),
            phase: "declared single solid phase".into(),
            process: "explicit finite-gas string source association".into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )?)
}

fn load_gas(options: &thermochemical::Options, volume_m3: f64) -> Result<IsochoricGasParcel> {
    let loaded = options.load()?;
    let bindings = loaded.bindings();
    let registry = loaded.registry(&bindings).ok_or(
        "finite gas requires explicit component packs, water identity and frozen phase policy",
    )?;
    let path = options
        .record
        .as_deref()
        .ok_or("--gas-record is required")?;
    let bytes = bounded_read(path, 1 << 20)?;
    let mut fields = std::str::from_utf8(&bytes)?.split_ascii_whitespace();
    let gas = read_isochoric_parcel(&mut fields, &registry, volume_m3)?;
    if fields.next().is_some() {
        return Err("gas record contains trailing fields".into());
    }
    Ok(gas)
}

fn context(temperature_k: f64, omega_rad_s: f64) -> Result<QueryPoint> {
    Ok(QueryPoint::new()
        .with("T", temperature_k)?
        .with("omega", omega_rad_s)?)
}

fn build(
    cx: &Cx<'_>,
    spec: &input::Spec,
    card: &MaterialCard,
    gas: &IsochoricGasParcel,
) -> Result<ThermalMaterialStringRuntime> {
    spec.validate()?;
    if !(spec.gas_solve.minimum_temperature_k..=spec.gas_solve.maximum_temperature_k)
        .contains(&gas.state().state().temperature)
    {
        return Err("initial gas temperature is outside its inversion bracket".into());
    }
    let point = context(spec.solid_temperature_k, spec.query_omega_rad_s)?;
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            "young_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
    ];
    if spec.drive.is_some() {
        requirements.push(ScalarPropertyRequirement::try_new(
            ELECTRICAL_RESISTIVITY_PROPERTY,
            ELECTRICAL_RESISTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    let resolved = resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let specimen = with_uniform_circular_material_and_constraints(
        PrestressedString {
            length_m: spec.length_m,
            tension_n: spec.tension_n,
            lin_density_kg_m: 1.0,
            axial_stiffness_n: 0.0,
            width_m: 1.0,
            n_modes: spec.modes,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
        },
        StringGeometryConstraint::FixedRadius(spec.radius_m),
        &resolved,
        StringPrestress::FixedTension(spec.tension_n),
    )?
    .with_kelvin_voigt_bending_loss()?;
    let mut mechanical = LinearMaterialStringRuntime::try_new_with_gas(
        cx,
        specimen,
        None,
        gas.state().state(),
        spec.listener_m,
        spec.sample_rate_hz,
    )?;
    mechanical.initialize_point_load_equilibrium(
        cx,
        &[StringPointForce {
            station_frac: spec.station,
            force_n: spec.force_n,
        }],
    )?;
    let lower = context(spec.solid_minimum_k, spec.query_omega_rad_s)?;
    let upper = context(spec.solid_maximum_k, spec.query_omega_rad_s)?;
    let heat_capacity = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let caloric = ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat_capacity,
            density: &density,
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        spec.enthalpy_interpolation_j_kg,
    )?;
    let initial = caloric.state_at_temperature(spec.solid_temperature_k)?;
    Ok(ThermalMaterialStringRuntime::try_new(
        cx,
        mechanical,
        card.clone(),
        caloric.curve().clone(),
        initial.specific_enthalpy_j_kg(),
        spec.maximum_solid_increment_k,
    )?)
}

fn mechanical_energy(live: &ThermalMaterialStringRuntime) -> f64 {
    live.mechanical()
        .modes()
        .iter()
        .zip(live.mechanical().states())
        .map(|(mode, state)| {
            0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2))
        })
        .sum()
}

struct History {
    description: String,
    pressure_pa: Vec<f64>,
    records: Vec<[f64; 22]>,
}

fn simulate(
    cx: &Cx<'_>,
    spec: &input::Spec,
    card: &MaterialCard,
    mut gas: IsochoricGasParcel,
) -> Result<History> {
    let mut live = build(cx, spec, card, &gas)?;
    let initial_mechanical_j = mechanical_energy(&live);
    let initial_enthalpy_j_kg = live.thermal().specific_enthalpy_j_kg();
    let initial_gas_energy_j_kg = gas.state().specific_internal_energy_j_kg();
    let solid_mass_kg = live.mechanical().specimen().mass_kg();
    let gas_mass_kg = gas.mass_kg();
    let description = format!(
        "material_card={}, caloric_curve={}, solid_mass_kg={solid_mass_kg:e}, gas_mass_kg={gas_mass_kg:e}, gas_volume_m3={:e}, initial_mechanical_energy_j={initial_mechanical_j:e}; fully solid, fixed geometry/tension/E/density, temperature-dependent loss/resistivity/conduction; first-order frozen-law samples; finite gas is a thermal inventory, not a chamber; diagnostic compact pressure supplies no radiation loading; source identity is not measured validation{}",
        card.content_hash().to_hex(),
        live.thermal().phase_curve_identity().to_hex(),
        gas.volume_m3(),
        thermochemical::provenance(&[GasRecord::Thermochemical(gas.state().clone())]),
    );
    let point = context(spec.solid_temperature_k, spec.query_omega_rad_s)?;
    let quantity = card
        .claims()
        .registered_quantity("thermal_conductivity")
        .ok_or("solid claims must contain thermal_conductivity")?;
    let conductivity = PropertyKey::with_quantity("thermal_conductivity", quantity);
    let layers = MaterialConductionLayers {
        material: card,
        context: &point,
        temperature_axis: "T",
        liquid_fraction_axis: None,
        conductivity: &conductivity,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: spec.transfer_area_m2,
        gas_distance_m: spec.gas_distance_m,
        body_distance_m: spec.solid_distance_m,
    };
    let controls = FiniteGasStringSolve {
        gas: spec.gas_solve,
        maximum_gas_temperature_increment_k: spec.maximum_gas_increment_k,
    };
    let dt = 1.0 / f64::from(spec.sample_rate_hz);
    let heat = FiniteGasStringHeat {
        string_heat_j: spec.solid_power_w * dt,
        gas_heat_j: spec.gas_power_w * dt,
    };
    let mut records = Vec::with_capacity(spec.samples);
    let mut pressure_pa = Vec::with_capacity(spec.samples);
    let (mut material_heat, mut drag_heat, mut exchange_heat, mut joule_heat) =
        (0.0, 0.0, 0.0, 0.0);
    let (mut solid_heat, mut gas_heat, mut work, mut other_loss, mut allowance) =
        (0.0, 0.0, 0.0, 0.0, 0.0);
    for sample in 0..spec.samples {
        cx.checkpoint().map_err(|_| "finite gas string cancelled")?;
        let forces = live.mechanical().project_point_forces(
            cx,
            &[StringPointForce {
                station_frac: spec.station,
                force_n: if sample < spec.release_sample {
                    spec.force_n
                } else {
                    0.0
                },
            }],
        )?;
        let epoch = live.mechanical().epoch();
        let (frame, supplied_joule, current, final_current, resistance, final_resistance) =
            if let Some(drive) = spec.drive {
                let report = live.step_with_drive_and_isochoric_gas_conduction(
                    cx,
                    epoch,
                    &forces,
                    drive,
                    &mut gas,
                    heat,
                    controls,
                    spec.exchange,
                    &layers,
                )?;
                let values = (
                    report.joule_heat_j,
                    report.current_a,
                    report.final_current_a,
                    report.conductor_before.resistance_ohm(),
                    report.conductor_after.resistance_ohm(),
                );
                (
                    report.coupled,
                    values.0,
                    values.1,
                    values.2,
                    values.3,
                    values.4,
                )
            } else {
                (
                    live.step_with_isochoric_gas_conduction(
                        cx,
                        epoch,
                        &forces,
                        &mut gas,
                        heat,
                        controls,
                        spec.exchange,
                        &layers,
                    )?,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                )
            };
        let coupled = &frame.coupled;
        let vibration = &coupled.string.vibration;
        material_heat += vibration.dissipation.material_heat_j;
        drag_heat += coupled.drag_heat_j;
        exchange_heat += frame.solid_to_gas_heat_j;
        joule_heat += supplied_joule;
        solid_heat += heat.string_heat_j;
        gas_heat += heat.gas_heat_j;
        work += vibration.acoustic.input_work_j;
        other_loss +=
            vibration.dissipation.authored_loss_j + vibration.dissipation.roundoff_residual_j;
        allowance += coupled.energy_tolerance_j;
        // Measure all three actual storage changes from their accepted states,
        // independently of the step reports' claimed energy residuals.
        let stored_change = mechanical_energy(&live) - initial_mechanical_j
            + solid_mass_kg * (live.thermal().specific_enthalpy_j_kg() - initial_enthalpy_j_kg)
            + gas_mass_kg * (gas.state().specific_internal_energy_j_kg() - initial_gas_energy_j_kg);
        let defect = stored_change + other_loss - work - joule_heat - solid_heat - gas_heat;
        let tolerance = allowance
            + 64.0
                * f64::EPSILON
                * (sample + 1) as f64
                * (stored_change.abs()
                    + other_loss.abs()
                    + work.abs()
                    + joule_heat.abs()
                    + solid_heat.abs()
                    + gas_heat.abs()
                    + initial_mechanical_j);
        if !defect.is_finite()
            || !tolerance.is_finite()
            || defect.abs() > tolerance
            || vibration.epoch != live.mechanical().epoch()
            || live.mechanical().accepted_samples() != (sample + 1) as u64
            || live.mechanical().ambient() != gas.state().state()
        {
            return Err("accepted live string/gas trajectory violates its shared epoch or cumulative energy balance".into());
        }
        let motion = live.mechanical().point_motion(cx, spec.station)?;
        let pressure = vibration.acoustic.observer_pressure_pa;
        let record = [
            (sample + 1) as f64 * dt,
            pressure,
            motion.displacement_m,
            motion.velocity_m_s,
            live.thermal().temperature_k(),
            gas.state().state().temperature,
            gas.state().state().pressure,
            material_heat,
            drag_heat,
            exchange_heat,
            joule_heat,
            solid_heat,
            gas_heat,
            work,
            stored_change,
            defect,
            tolerance,
            current,
            final_current,
            resistance,
            final_resistance,
            live.mechanical().modes()[0].damping_ratio,
        ];
        if record.iter().any(|value| !value.is_finite()) {
            return Err("nonfinite trajectory observation".into());
        }
        records.push(record);
        pressure_pa.push(pressure);
    }
    cx.checkpoint().map_err(|_| "finite gas string cancelled")?;
    Ok(History {
        description,
        pressure_pa,
        records,
    })
}

fn run(args: &[String]) -> Result<()> {
    let (args, gas_options) = thermochemical::Options::extract(args)?;
    if args.len() != 3 {
        return Err("usage: finite_gas_string SCENARIO SOLID_CLAIMS.pack OUTPUT.wav --gas-record FILE --gas-component SPECIES.pack NASA9.pack TRANSPORT.pack [repeat] --gas-water none|ID --gas-phase frozen".into());
    }
    let bytes = bounded_read(Path::new(&args[0]), input::MAX_INPUT_BYTES)?;
    let spec = input::parse(std::str::from_utf8(&bytes)?)?;
    let pack = NormalizedPack::from_bytes(&bounded_read(Path::new(&args[1]), 16 << 20)?)?;
    let card = material(&pack)?;
    let gas = load_gas(&gas_options, spec.gas_volume_m3)?;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let history = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x464753,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        simulate(&cx, &spec, &card, gas)
    })?;
    let (wav, clips) = encode_pcm16_wav(
        &history.pressure_pa,
        spec.sample_rate_hz,
        spec.full_scale_pa,
    )?;
    // Physics admission completes before either output is exposed. Creating a
    // new file also protects an earlier accepted render from accidental reuse.
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[2])?;
    file.write_all(&wav)?;
    file.flush()?;
    eprintln!(
        "{}; sample_rate_hz={}; full_scale_pa={:e}; clipped_samples={clips}",
        history.description, spec.sample_rate_hz, spec.full_scale_pa
    );
    let mut output = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "time_s,pressure_pa,displacement_m,velocity_m_s,solid_temperature_k,gas_temperature_k,gas_pressure_pa,material_heat_j,drag_heat_j,solid_to_gas_heat_j,joule_heat_j,prescribed_solid_heat_j,prescribed_gas_heat_j,mechanical_work_j,storage_change_j,energy_defect_j,energy_tolerance_j,current_a,endpoint_current_a,resistance_ohm,endpoint_resistance_ohm,first_mode_damping_ratio"
    )?;
    for record in &history.records {
        for (i, value) in record.iter().enumerate() {
            if i != 0 {
                write!(output, ",")?;
            }
            write!(output, "{value:.17e}")?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    run(&std::env::args().skip(1).collect::<Vec<_>>())
}

#[cfg(test)]
#[path = "finite_gas_string/tests.rs"]
mod tests;
