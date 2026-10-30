//! Source packs -> electrical/thermal feedback -> continuous string motion and sound.
//!
//! The existing string, conductor, caloric, finite-gas and conduction owners
//! advance one shared accepted trajectory. CSV contains physical pressure and
//! independent cumulative storage accounting; WAV has an explicit fixed scale.
//! See finite_gas_string/README.md for source compilation and the bounded input.

#[path = "finite_gas_string/input.rs"]
mod input;
#[path = "finite_gas_string/layered.rs"]
mod layered;
#[path = "finite_gas_string/flow.rs"]
mod flow;
#[path = "finite_gas_string/output.rs"]
mod output;
#[path = "../src/bin/music_render/wind_input/attachments.rs"]
mod attachments;
#[path = "material_strings/native_input.rs"]
mod native_sources;
#[allow(dead_code)]
#[path = "../src/bin/music_render/json.rs"]
mod json;
#[allow(dead_code)]
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
use json::json_string;

use fs_couple::acoustic_realize::{
    LinearMaterialStringRuntime, StringPointForce, ThermalMaterialStringRuntime,
    finite_gas::{
        FiniteGasStringHeat, FiniteGasStringSolve, flow::FiniteGasStringFlowDirection,
    },
};
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::string_specimen::{
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, StringGeometryConstraint, StringPrestress,
    with_uniform_circular_material_and_constraints, with_uniform_circular_thermal_extension,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::gas_species::input::{GasRecord, read_isochoric_parcel};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point, integrate_isotropic_thermal_expansion,
};
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_qty::{Density, DynViscosity, Pressure};
use fs_scenario::PrestressedString;
use std::{
    io::{Read, Write},
    path::Path,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
// The existing native string adapter and this runner share one source reader.
use bounded_read as read_bounded;

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

fn load_gas(
    options: &thermochemical::Options, spec: &input::Spec, directory: &Path,
    mut read: impl FnMut(&Path, usize) -> Result<Vec<u8>>,
) -> Result<(IsochoricGasParcel, Option<ClosedGasParcel>)> {
    spec.validate()?;
    let loaded = options.load_with(|path| read(path, 16 << 20).map_err(|e| e.to_string()))?;
    let bindings = loaded.bindings();
    let registry = loaded.registry(&bindings).ok_or(
        "finite gas requires explicit component packs, water identity and frozen phase policy",
    )?;
    let mut parcel = |path: &Path, volume_m3| -> Result<IsochoricGasParcel> {
        let bytes = read(path, 1 << 20)?;
        if bytes.len() > 1 << 20 { return Err("gas record exceeds 1 MiB".into()); }
        let mut fields = std::str::from_utf8(&bytes)?.split_ascii_whitespace();
        let gas = read_isochoric_parcel(&mut fields, &registry, volume_m3)?;
        if fields.next().is_some() { return Err("gas record contains trailing fields".into()); }
        Ok(gas)
    };
    let gas = parcel(options.record.as_deref().ok_or("--gas-record is required")?,
        spec.gas_volume_m3)?;
    let reservoir = spec.nozzle.as_ref().map(|nozzle|
        parcel(&directory.join(&nozzle.reservoir_record), nozzle.reservoir_volume_m3)
            .map(ClosedGasParcel::from_isochoric)).transpose()?;
    Ok((gas, reservoir))
}

fn context(temperature_k: f64, omega_rad_s: f64) -> Result<QueryPoint> {
    Ok(QueryPoint::new()
        .with("T", temperature_k)?
        .with("omega", omega_rad_s)?)
}

fn string_template(spec: &input::Spec) -> PrestressedString {
    PrestressedString {
        length_m: spec.length_m,
        // The material/prestress binding replaces all template properties.
        tension_n: 0.0,
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
        tension_loss_factor: None,
    }
}

fn caloric_curve(spec: &input::Spec, card: &MaterialCard) -> Result<ResolvedHeatCapacityCurve> {
    let lower = context(spec.solid_minimum_k, spec.query_omega_rad_s)?;
    let upper = context(spec.solid_maximum_k, spec.query_omega_rad_s)?;
    let heat_capacity = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    Ok(ResolvedHeatCapacityCurve::try_new(
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
    )?)
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
    let template = string_template(spec);
    let geometry = StringGeometryConstraint::FixedRadius(spec.radius_m);
    let specimen = match (spec.prestress, spec.thermoelastic) {
        (input::Prestress::FixedTension(tension), None) =>
            with_uniform_circular_material_and_constraints(template, geometry, &resolved,
                StringPrestress::FixedTension(tension))?,
        (input::Prestress::FixedSupports, Some(thermal)) => {
            let reference = context(thermal.reference_temperature_k, spec.query_omega_rad_s)?;
            let expansion = integrate_isotropic_thermal_expansion(card, &reference, &point,
                MaterialPropertySelection::SingleClaimOnly)?;
            with_uniform_circular_thermal_extension(template, geometry, &resolved, &expansion,
                thermal.reference_length_m, thermal.strain_limit)?
        }
        _ => return Err("inconsistent string prestress and thermoelastic input".into()),
    }.with_kelvin_voigt_bending_loss()?;
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
    let caloric = caloric_curve(spec, card)?;
    let initial = caloric.state_at_temperature(spec.solid_temperature_k)?;
    Ok(if let Some(thermal) = spec.thermoelastic {
        ThermalMaterialStringRuntime::try_new_thermoelastic(cx, mechanical, card.clone(),
            caloric.curve().clone(), initial.specific_enthalpy_j_kg(),
            spec.maximum_solid_increment_k, thermal.solve)?
    } else {
        ThermalMaterialStringRuntime::try_new(cx, mechanical, card.clone(),
            caloric.curve().clone(), initial.specific_enthalpy_j_kg(),
            spec.maximum_solid_increment_k)?
    })
}

struct History {
    description: String,
    pressure_pa: Vec<f64>,
    records: Vec<[f64; 27]>,
    layer_temperatures: Vec<Vec<f64>>,
    gas_flow: Option<flow::History>,
    electrical_feedback: Option<Vec<[f64; 3]>>,
}

fn simulate(
    cx: &Cx<'_>,
    spec: &input::Spec,
    card: &MaterialCard,
    mut gas: IsochoricGasParcel,
    mut reservoir: Option<ClosedGasParcel>,
) -> Result<History> {
    spec.validate()?;
    let valid_reservoir = match (&spec.nozzle, &reservoir) {
        (None, None) => true,
        (Some(nozzle), Some(parcel)) => {
            parcel.volume_m3().to_bits() == nozzle.reservoir_volume_m3.to_bits()
                && (spec.gas_solve.minimum_temperature_k..=spec.gas_solve.maximum_temperature_k)
                    .contains(&parcel.state().state().temperature)
        }
        _ => false,
    };
    if !valid_reservoir {
        return Err("nozzle data requires its original finite reservoir and admitted initial temperature".into());
    }
    let mut live = build(cx, spec, card, &gas)?;
    let initial_mechanical_j = live.mechanical_internal_energy_j()?;
    let initial_enthalpy_j_kg = live.thermal().specific_enthalpy_j_kg();
    let initial_gas_energy_j_kg = gas.state().specific_internal_energy_j_kg();
    let solid_mass_kg = live.mechanical().specimen().mass_kg();
    let gas_mass_kg = gas.mass_kg();
    let prestress = if spec.thermoelastic.is_some() {
        "fixed supports, sourced thermal expansion and affine stiffness with reversible heat"
    } else { "fixed tension and constant elastic modulus" };
    let electrical_scope = if spec.electrical_feedback.is_some() {
        "first-order frozen mechanics with solved endpoint Joule heat"
    } else {
        "first-order frozen-law samples"
    };
    let mut description = format!(
        "material_card={}, caloric_curve={}, solid_mass_kg={solid_mass_kg:e}, gas_mass_kg={gas_mass_kg:e}, gas_volume_m3={:e}, initial_mechanical_internal_energy_j={initial_mechanical_j:e}; fully solid, fixed geometry/density, {prestress}, temperature-dependent loss/resistivity/conduction; {electrical_scope}; finite gas is a thermal inventory, not a chamber; diagnostic compact pressure supplies no radiation loading; source identity is not measured validation{}",
        card.content_hash().to_hex(),
        live.thermal().phase_curve_identity().to_hex(),
        gas.volume_m3(),
        thermochemical::provenance(&[GasRecord::Thermochemical(gas.state().clone())]),
    );
    let mut gas_flow = reservoir
        .as_ref()
        .map(|parcel| flow::History::new(&gas, parcel, spec.samples));
    if let (Some(nozzle), Some(parcel)) = (&spec.nozzle, &reservoir) {
        description.push_str(&format!(
            "; nozzle={:?}, reservoir_volume_m3={:e}; first-order heat-then-flow, finite well-mixed volumes, no jet momentum or chamber load; reservoir{}",
            nozzle.geometry,
            parcel.volume_m3(),
            thermochemical::provenance(&[GasRecord::Thermochemical(parcel.state().clone())]),
        ));
    }
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
    let mut electrical_feedback = spec.electrical_feedback.map(|_| Vec::with_capacity(spec.samples));
    let (mut electrical_source_work, mut electrical_source_residual) = (0.0, 0.0);
    let (mut material_heat, mut drag_heat, mut exchange_heat, mut joule_heat) =
        (0.0, 0.0, 0.0, 0.0);
    let (mut solid_heat, mut gas_heat, mut work, mut other_loss, mut allowance) =
        (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut reversible_heat, mut axial_change, mut entropy_energy_change) = (0.0, 0.0, 0.0);
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
        let (frame, electrical, flux, source) =
            if let (Some(nozzle), Some(parcel)) = (&spec.nozzle, reservoir.as_mut()) {
                let (report, electrical, source) = if let Some(drive) = spec.drive {
                    let (report, source) = if let Some(feedback) = spec.electrical_feedback {
                        let report = live.step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
                            cx, epoch, &forces, drive, &mut gas, parcel, heat, controls,
                            spec.exchange, &layers, nozzle.geometry, nozzle.controls(spec.gas_solve),
                            feedback,
                        )?;
                        let source = (report.source_work_j, report.source_residual_j,
                            report.trials, report.energy_tolerance_j);
                        (report.electrical, Some(source))
                    } else {
                        (live.step_with_drive_and_isochoric_gas_nozzle_conduction(
                            cx, epoch, &forces, drive, &mut gas, parcel, heat, controls,
                            spec.exchange, &layers, nozzle.geometry, nozzle.controls(spec.gas_solve),
                        )?, None)
                    };
                    let electrical = [
                        report.joule_heat_j,
                        report.current_a,
                        report.final_current_a,
                        report.conductor_before.resistance_ohm(),
                        report.conductor_after.resistance_ohm(),
                    ];
                    (report.coupled, electrical, source)
                } else {
                    (
                        live.step_with_isochoric_gas_nozzle_conduction(
                            cx, epoch, &forces, &mut gas, parcel, heat, controls,
                            spec.exchange, &layers, nozzle.geometry, nozzle.controls(spec.gas_solve),
                        )?,
                        [0.0; 5],
                        None,
                    )
                };
                let sign = match report.direction {
                    FiniteGasStringFlowDirection::IntoSurroundings => 1.0,
                    FiniteGasStringFlowDirection::OutOfSurroundings => -1.0,
                };
                let flux = (
                    sign * report.flow.mass_flow_kg_s,
                    sign * report.transfer.transferred_mass_kg,
                    sign * report.transfer.transferred_enthalpy_j,
                    report.energy_tolerance_j,
                );
                (report.thermal, electrical, Some(flux), source)
            } else if let Some(drive) = spec.drive {
                let (report, source) = if let Some(feedback) = spec.electrical_feedback {
                    let report = live.step_with_implicit_drive_and_isochoric_gas_conduction(
                        cx, epoch, &forces, drive, &mut gas, heat, controls,
                        spec.exchange, &layers, feedback,
                    )?;
                    let source = (report.source_work_j, report.source_residual_j,
                        report.trials, report.energy_tolerance_j);
                    (report.electrical, Some(source))
                } else {
                    (live.step_with_drive_and_isochoric_gas_conduction(
                        cx, epoch, &forces, drive, &mut gas, heat, controls,
                        spec.exchange, &layers,
                    )?, None)
                };
                let electrical = [
                    report.joule_heat_j,
                    report.current_a,
                    report.final_current_a,
                    report.conductor_before.resistance_ohm(),
                    report.conductor_after.resistance_ohm(),
                ];
                (report.coupled, electrical, None, source)
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
                    [0.0; 5],
                    None,
                    None,
                )
            };
        let [supplied_joule, current, final_current, resistance, final_resistance] = electrical;
        let coupled = &frame.coupled;
        let vibration = &coupled.string.vibration;
        material_heat += vibration.dissipation.material_heat_j;
        drag_heat += coupled.drag_heat_j;
        exchange_heat += frame.solid_to_gas_heat_j;
        joule_heat += supplied_joule;
        electrical_source_work += source.as_ref().map_or(supplied_joule, |source| source.0);
        solid_heat += heat.string_heat_j;
        gas_heat += heat.gas_heat_j;
        work += vibration.acoustic.input_work_j;
        other_loss +=
            vibration.dissipation.authored_loss_j + vibration.dissipation.roundoff_residual_j;
        allowance += source.as_ref().map_or_else(
            || flux.as_ref().map_or(coupled.energy_tolerance_j, |flux| flux.3),
            |source| source.3,
        );
        reversible_heat += coupled.string.reversible_heat_j;
        axial_change += coupled.string.axial_free_energy_change_j;
        entropy_energy_change += coupled.string.mechanical_entropy_energy_change_j;
        // Measure accepted storage independently of the step reports. With
        // transport, both extensive inventories replace the fixed-mass term.
        let solid_and_motion = live.mechanical_internal_energy_j()? - initial_mechanical_j
            + solid_mass_kg * (live.thermal().specific_enthalpy_j_kg() - initial_enthalpy_j_kg);
        let (stored_change, gas_reference_scale) =
            if let (Some(history), Some(parcel)) = (&gas_flow, &reservoir) {
                let changes = history.energy_changes(&gas, parcel);
                (
                    solid_and_motion + changes[0] + changes[1],
                    history.reference_energy_scale(&gas, parcel),
                )
            } else {
                (
                    solid_and_motion
                        + gas_mass_kg * (gas.state().specific_internal_energy_j_kg() - initial_gas_energy_j_kg),
                    0.0,
                )
            };
        let defect = stored_change + other_loss - work - electrical_source_work - solid_heat - gas_heat;
        let tolerance = allowance
            + 64.0
                * f64::EPSILON
                * (sample + 1) as f64
                * (stored_change.abs()
                    + other_loss.abs()
                    + work.abs()
                    + electrical_source_work.abs()
                    + solid_heat.abs()
                    + gas_heat.abs()
                    + initial_mechanical_j.abs()
                    + gas_reference_scale);
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
            live.mechanical().specimen().string().tension_n,
            live.mechanical().modes()[0].angular_frequency_rad_s / (2.0 * std::f64::consts::PI),
            reversible_heat,
            axial_change,
            entropy_energy_change,
        ];
        if record.iter().any(|value| !value.is_finite()) {
            return Err("nonfinite trajectory observation".into());
        }
        if let (Some(history), Some(parcel), Some(flux)) = (&mut gas_flow, &reservoir, flux) {
            history.record(&gas, parcel, flux.0, flux.1, flux.2)?;
        }
        if let (Some(feedback), Some(source)) = (&mut electrical_feedback, source) {
            electrical_source_residual += source.1;
            let row = [electrical_source_work, electrical_source_residual, source.2 as f64];
            if row.iter().any(|value| !value.is_finite()) {
                return Err("nonfinite electrical feedback observation".into());
            }
            feedback.push(row);
        }
        records.push(record);
        pressure_pa.push(pressure);
    }
    cx.checkpoint().map_err(|_| "finite gas string cancelled")?;
    Ok(History {
        description,
        pressure_pa,
        records,
        layer_temperatures: Vec::new(),
        gas_flow,
        electrical_feedback,
    })
}

fn run(args: &[String]) -> Result<()> {
    let (args, gas_options) = thermochemical::Options::extract(args)?;
    let (args, output_options) = output::Options::extract(&args)?;
    let (args, source_options) = native_sources::SourceOptions::extract(&args)?;
    if args.len() != 3 {
        return Err("usage: finite_gas_string SCENARIO SOLID_CLAIMS.pack OUTPUT.wav --gas-record FILE --gas-component SPECIES.pack NASA9.pack TRANSPORT.pack [repeat] --gas-water none|ID --gas-phase frozen [--self-contained --attached-source RELATIVE_PATH ...] [--history-output NEW.csv] [--require-sample-rate HZ]".into());
    }
    if source_options.decimate {
        return Err("finite_gas_string does not support --decimate; author the physical sample rate explicitly".into());
    }
    let bytes = bounded_read(Path::new(&args[0]), input::MAX_INPUT_BYTES)?;
    let spec = input::parse(std::str::from_utf8(&bytes)?)?;
    output_options.admit(spec.sample_rate_hz, Path::new(&args[2]))?;
    let directory = Path::new(&args[0]).parent().unwrap_or_else(|| Path::new("."));
    let sources = source_options.load(directory)?;
    let pack = NormalizedPack::from_bytes(&sources.read_option(Path::new(&args[1]), 16 << 20)?)?;
    let mut cards = vec![material(&pack)?];
    if let Some(bonded) = &spec.bonded {
        for sleeve in &bonded.sleeves {
            let bytes = sources.read_relative(&sleeve.pack, 16 << 20)?;
            cards.push(material(&NormalizedPack::from_bytes(&bytes)?)?);
        }
    }
    // Attached names remain exact logical paths; ordinary CLI record names
    // retain their historical scenario-relative reservoir resolution.
    let gas_directory = if sources.is_restricted() { Path::new("") } else { directory };
    let (gas, reservoir) = load_gas(&gas_options, &spec, gas_directory,
        |path, limit| sources.read_option(path, limit))?;
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
        if spec.bonded.is_some() {
            layered::simulate(&cx, &spec, &cards, gas, reservoir)
        } else {
            simulate(&cx, &spec, &cards[0], gas, reservoir)
        }
    })?;
    let (wav, clips) = encode_pcm16_wav(
        &history.pressure_pa,
        spec.sample_rate_hz,
        spec.full_scale_pa,
    )?;
    let csv = output_options.history.as_ref()
        .map(|_| output::csv_bytes(&history, output::MAX_HISTORY_BYTES)).transpose()?;
    // Physics admission completes before either output is exposed. Creating a
    // new file also protects an earlier accepted render from accidental reuse.
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[2])?;
    file.write_all(&wav)?;
    file.flush()?;
    if let (Some(path), Some(csv)) = (&output_options.history, csv) {
        let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(path)?;
        file.write_all(&csv)?;
        file.flush()?;
    } else {
        output::write_csv(&history, std::io::BufWriter::new(std::io::stdout().lock()))?;
    }
    eprintln!(
        "{}; sample_rate_hz={}; full_scale_pa={:e}; clipped_samples={clips}",
        history.description, spec.sample_rate_hz, spec.full_scale_pa
    );
    // The same calibrated receipt grammar as the existing material performer.
    eprintln!("samples={}; rate_Hz={}; clips={clips}; full_scale_Pa={:e}",
        history.pressure_pa.len(), spec.sample_rate_hz, spec.full_scale_pa);
    Ok(())
}

fn main() -> Result<()> {
    run(&std::env::args().skip(1).collect::<Vec<_>>())
}

#[cfg(test)]
#[path = "finite_gas_string/tests.rs"]
mod tests;
