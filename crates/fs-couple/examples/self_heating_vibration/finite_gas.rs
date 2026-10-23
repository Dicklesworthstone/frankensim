//! Finite surrounding gas on the existing source-bound vibration example.
//! All reservoir proposals remain inside the original material feedback loop.

#[path = "../../src/bin/music_render/json.rs"]
mod json;
use json::json_string;
#[allow(dead_code)]
#[path = "../../src/bin/music_render/thermochemical.rs"]
mod thermochemical;

use super::{AREA_M2, LENGTH_M, Result, source, surface_area};
use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::finite_gas_thermal::FiniteGasThermalEnvironment;
use fs_couple::thermovisco_dynamics::{Controls, Geometry, Oscillator};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::gas_species::input::{GasRecord, read_isochoric_parcel};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_material::thermovisco::thermal::FeedbackConfig;
use std::io::Read;

#[derive(Clone, Copy)]
pub(super) struct Motion {
    pub velocity: f64,
    pub payload: f64,
    pub duration: f64,
    pub count: usize,
    pub film: f64,
}

#[derive(Clone, Copy, Debug)]
struct Reservoir {
    volume_m3: f64,
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
    maximum_exchange_j: f64,
    power_w: f64,
    solid_temperature_k: f64,
}

#[derive(Debug)]
pub(super) struct Options {
    sources: thermochemical::Options,
    reservoir: Reservoir,
}

impl Options {
    /// Preserve the old positional command and its output when no finite
    /// reservoir is requested. Missing numerical/source policy refuses before I/O.
    pub fn extract(args: &[String]) -> Result<(Vec<String>, Option<Self>)> {
        if !args.iter().any(|arg| {
            arg == "--finite-gas" || arg == "--solid-temperature" || arg.starts_with("--gas-")
        }) {
            return Ok((args.to_vec(), None));
        }
        let (remaining, sources) = thermochemical::Options::extract(args)?;
        let mut words = remaining.iter();
        let mut positional = Vec::new();
        let mut finite = None;
        let mut power = None;
        let mut temperature = None;
        while let Some(word) = words.next() {
            match word.as_str() {
                "--finite-gas" => {
                    if finite.is_some() {
                        return Err("--finite-gas may only be supplied once".into());
                    }
                    let mut number = || -> Result<f64> {
                        Ok(words
                            .next()
                            .ok_or("--finite-gas needs VOLUME_M3 TMIN_K TMAX_K MAX_EXCHANGE_J")?
                            .parse()?)
                    };
                    finite = Some([number()?, number()?, number()?, number()?]);
                }
                "--gas-power" => {
                    if power.is_some() {
                        return Err("--gas-power may only be supplied once".into());
                    }
                    power = Some(
                        words
                            .next()
                            .ok_or("--gas-power needs signed watts")?
                            .parse::<f64>()?,
                    );
                }
                "--solid-temperature" => {
                    if temperature.is_some() {
                        return Err("--solid-temperature may only be supplied once".into());
                    }
                    temperature = Some(
                        words
                            .next()
                            .ok_or("--solid-temperature needs kelvin")?
                            .parse::<f64>()?,
                    );
                }
                _ => positional.push(word.clone()),
            }
        }
        let [volume, minimum, maximum, heat] =
            finite.ok_or("gas source/thermal options require --finite-gas")?;
        let power_w = power.unwrap_or(0.0);
        let solid_temperature_k = temperature.unwrap_or(300.0);
        if ![volume, minimum, maximum, heat, solid_temperature_k]
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
            || maximum <= minimum
            || !power_w.is_finite()
        {
            return Err("finite gas requires positive volume, temperatures and heat bracket; ordered temperature bounds and finite gas power".into());
        }
        if sources.record.is_none() {
            return Err(
                "--finite-gas requires --gas-record and explicit thermochemical component sources"
                    .into(),
            );
        }
        Ok((
            positional,
            Some(Self {
                sources,
                reservoir: Reservoir {
                    volume_m3: volume,
                    minimum_temperature_k: minimum,
                    maximum_temperature_k: maximum,
                    maximum_exchange_j: heat,
                    power_w,
                    solid_temperature_k,
                },
            }),
        ))
    }

    fn load(&self) -> Result<IsochoricGasParcel> {
        let loaded = self.sources.load()?;
        let bindings = loaded.bindings();
        let registry = loaded.registry(&bindings).ok_or(
            "finite gas requires explicit --gas-component, --gas-water and --gas-phase frozen",
        )?;
        let path = self.sources.record.as_ref().ok_or("missing gas record")?;
        let limit = fs_material::gas_species::input::packs::MAX_RECORD_BYTES;
        let mut bytes = Vec::new();
        std::fs::File::open(path)?
            .take((limit + 1) as u64)
            .read_to_end(&mut bytes)?;
        if bytes.len() > limit {
            return Err("gas record exceeds 8 KiB".into());
        }
        let text = std::str::from_utf8(&bytes)?;
        let mut fields = text.split_ascii_whitespace();
        let gas = read_isochoric_parcel(&mut fields, &registry, self.reservoir.volume_m3)?;
        if fields.next().is_some() {
            return Err("unexpected fields after finite gas declaration".into());
        }
        admit_initial(&gas, self.reservoir)?;
        Ok(gas)
    }
}

fn admit_initial(gas: &IsochoricGasParcel, reservoir: Reservoir) -> Result<()> {
    let temperature = gas.state().state().temperature;
    if temperature < reservoir.minimum_temperature_k
        || temperature > reservoir.maximum_temperature_k
    {
        return Err(
            "initial gas temperature is outside the declared energy-inversion bracket".into(),
        );
    }
    Ok(())
}

struct History {
    description: String,
    records: Vec<[f64; 14]>,
}

/// This driver composes the existing mechanical and thermal owners. It does
/// not introduce a different oscillator, caloric law or convergence procedure.
fn simulate(
    cx: &Cx<'_>,
    gate: &CancelGate,
    input: &source::Source,
    mut gas: IsochoricGasParcel,
    motion: Motion,
    reservoir: Reservoir,
) -> Result<History> {
    admit_initial(&gas, reservoir)?;
    let mut state = Oscillator::try_new(
        &input.mechanical,
        &input.caloric,
        Geometry {
            length_m: LENGTH_M,
            cross_section_m2: AREA_M2,
            thermal_surface_m2: surface_area(),
            end_payload_kg: motion.payload,
        },
        reservoir.solid_temperature_k,
        0.0,
        motion.velocity,
    )?;
    let source_details =
        thermochemical::provenance(&[GasRecord::Thermochemical(gas.state().clone())]);
    let description = format!(
        "material_card={}, caloric_curve={}, branches={}; L={LENGTH_M} m, A={AREA_M2} m2; specimen_mass={} kg; modal_mass={} kg; gas_mass={} kg; gas_volume={} m3; initial_solid_temperature={} K; initial_gas_temperature={} K; film={} W/(m2 K); external_gas_power={} W; source_details={{{}}}. Fixed geometry, single solid phase, uniform axial reduction and well-mixed fixed-volume frozen gas; radiation disabled. Film coefficient is declared, not inferred from source transport or geometry. Source data do not establish measured vibration accuracy.",
        input.card.content_hash().to_hex(),
        input.caloric.curve().identity().to_hex(),
        input.mechanical.law().branch_count(),
        state.point().mass_kg(),
        state.effective_mass_kg(),
        gas.mass_kg(),
        gas.volume_m3(),
        reservoir.solid_temperature_k,
        gas.state().state().temperature,
        motion.film,
        reservoir.power_w,
        source_details.trim_start_matches(','),
    );
    let body = LumpedThermalEnvironment {
        temperature_k: gas.state().state().temperature,
        radiation_temperature_k: 300.0,
        convection_w_per_m2_k: motion.film,
        transport: input.transport.clone(),
        enthalpy_tolerance_j_kg: 1e-12,
        maximum_thermal_residual_j: 1e-10,
    };
    let controls = Controls {
        feedback: FeedbackConfig {
            max_iterations: 100,
            temperature_tolerance_k: 1e-8,
            energy_tolerance_j: 1e-9,
            relaxation: 1.0,
        },
        force_tolerance_n: 1e-6,
        energy_tolerance_j: 2e-9,
    };
    let gas_solve = GasEnergySolve {
        minimum_temperature_k: reservoir.minimum_temperature_k,
        maximum_temperature_k: reservoir.maximum_temperature_k,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let exchange = HeatExchangeSolve {
        heat_bracket_j: [-reservoir.maximum_exchange_j, reservoir.maximum_exchange_j],
        source_tolerance_j: 1e-12,
        energy_tolerance_j: 2e-10,
        maximum_iterations: 100,
    };
    let initial_energy = state.kinetic_energy_j()
        + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3();
    let initial_h = state.point().caloric_state().specific_enthalpy_j_kg();
    let initial_u = gas.state().specific_internal_energy_j_kg();
    let gas_mass = gas.mass_kg();
    let (mut numerical, mut material_heat, mut transferred, mut gas_source, mut radiation) =
        (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut records = Vec::with_capacity(motion.count);
    for index in 1..=motion.count {
        cx.checkpoint()
            .map_err(|_| "cancelled finite-gas vibration")?;
        let target = motion.duration * (index as f64 / motion.count as f64);
        let dt = target - state.point().time_s();
        let gas_heat_j = reservoir.power_w * dt;
        if !gas_heat_j.is_finite() || (reservoir.power_w != 0.0 && gas_heat_j == 0.0) {
            return Err("prescribed gas heat is unrepresentable".into());
        }
        let environment = FiniteGasThermalEnvironment {
            gas: &gas,
            body: &body,
            gas_heat_j,
            gas_solve,
            exchange,
            maximum_gas_temperature_increment_k: reservoir.maximum_temperature_k
                - reservoir.minimum_temperature_k,
        };
        let step = state.step(
            0.0,
            dt,
            controls,
            |thermal| environment.advance(cx, thermal),
            gate,
        )?;
        let next = &step.next;
        let transport = &step.transport.report;
        numerical += step.numerical_loss_j;
        material_heat += step.material_heat_j;
        transferred += transport.solid_to_gas_heat_j;
        gas_source += transport.external_gas_heat_j;
        radiation += transport.radiation_heat_j;
        let final_u = transport.next_gas.state().specific_internal_energy_j_kg();
        let stored = next.kinetic_energy_j()
            + next.point().volume_m3() * next.point().mechanical().stored_energy_j_m3()
            + next.point().mass_kg()
                * (next.point().caloric_state().specific_enthalpy_j_kg() - initial_h)
            + gas_mass * (final_u - initial_u);
        // Independently measured accepted storage: Q and material loss are
        // internal transfers, and backward-Euler loss is never gas/solid heat.
        let defect = stored + numerical - initial_energy - gas_source - radiation;
        let scale = stored.abs()
            + numerical.abs()
            + initial_energy.abs()
            + gas_source.abs()
            + radiation.abs()
            + gas_mass * (initial_u.abs() + final_u.abs());
        let tolerance = index as f64
            * (controls.energy_tolerance_j + gas_mass * gas_solve.energy_tolerance_j_kg)
            + 64.0 * f64::EPSILON * scale;
        if !defect.is_finite() || !tolerance.is_finite() || defect.abs() > tolerance {
            return Err("cumulative specimen/inertia/gas energy balance refused".into());
        }
        records.push([
            next.point().time_s(),
            next.extension_m(),
            next.velocity_m_s(),
            next.point().caloric_state().temperature_k(),
            transport.next_gas.state().state().temperature,
            transport.next_gas.state().state().pressure,
            material_heat,
            numerical,
            transferred,
            gas_source,
            radiation,
            stored,
            defect,
            tolerance,
        ]);
        cx.checkpoint()
            .map_err(|_| "cancelled finite-gas publication")?;
        state = step.next;
        gas = step.transport.report.next_gas;
    }
    Ok(History {
        description,
        records,
    })
}

pub(super) fn run(
    mode: &str,
    path: &str,
    branches: Option<usize>,
    motion: Motion,
    options: Options,
) -> Result<()> {
    let gas = options.load()?;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let history = pool.scope(|arena| -> Result<_> {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x544d44,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let input = if mode == "--pack" {
            source::Source::pack(path, branches.ok_or("missing branch count")?, &gate)?
        } else {
            source::Source::table(path, &gate)?
        };
        simulate(&cx, &gate, &input, gas, motion, options.reservoir)
    })?;
    // Publish only a complete accepted run. Admission, source, convergence and
    // cancellation failures cannot leave a misleading partial CSV trajectory.
    eprintln!("{}", history.description);
    println!(
        "time_s,extension_m,velocity_m_s,solid_temperature_K,gas_temperature_K,gas_pressure_Pa,physical_material_heat_J,numerical_loss_J,solid_to_gas_heat_J,external_gas_heat_J,radiation_heat_J,stored_energy_J,cumulative_energy_defect_J,cumulative_energy_tolerance_J"
    );
    for row in history.records {
        println!(
            "{}",
            row.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
        );
    }
    Ok(())
}

#[cfg(test)]
#[path = "finite_gas_tests.rs"]
mod tests;
