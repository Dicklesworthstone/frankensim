//! Authored finite wall enthalpy and reciprocal gas/wall heat.
//!
//! The original DG wall declaration owns the material and complete condition
//! point. This adapter selects its caloric properties, finite inventory and
//! external powers; existing material and coupled owners compute every state.
//!
//! After optional finite-gas records and before observation, the bounded grammar
//! is (SI units throughout; samples use the mechanical clock):
//! ```text
//! finite_wall solid
//! wall_caloric CP_PROPERTY DENSITY_PROPERTY T TMIN TMAX HREF MAX_INTERPOLATION_ERROR single-claim
//! wall_inventory MASS_KG VOLUME_M3 AREA_M2 INITIAL_SPECIFIC_ENTHALPY_J_KG
//! wall_feedback WMIN WMAX SOURCE_TOL ENERGY_TOL MAX_ITERS MAX_ABS_PARAMETER_WORK
//! wall_dissipation none
//! wall_heat COUNT
//! wall_power SAMPLE POWER_W
//! ```
//! Repeat wall_power COUNT times. Native mechanical loss heating instead uses
//! `wall_dissipation native MAX_DIM MAX_BYTES MAX_MULTIPLY_TERMS`.
//! With finite_gas, append both reciprocal interface records:
//! ```text
//! gas_wall_conduction K_PROPERTY AREA_M2 GAS_DISTANCE_M BODY_DISTANCE_M
//! gas_wall_feedback QMIN QMAX SOURCE_TOL ENERGY_TOL MAX_ITERS MAX_MATERIAL_TRIALS
//! ```
//! The initial enthalpy reconstructs the original wall temperature through its
//! source-derived all-solid chart. HREF is the specific enthalpy at TMIN. Signed
//! external powers are held between events; paired events consume two shared
//! control/work slots. This model fixes the existing geometry and inventories.

use super::finite_gas::FiniteGasInput;
use super::{AperturePerformance, MAX_CONTROLS, PlateValveInputError, Reader, bad, checkpoint};
use crate::bernoulli_aperture::dg::thermal_wall::FiniteWallThermalSolve;
use crate::bernoulli_aperture::dg::thermal_wall::gas::FiniteGasWallThermalSolve;
use crate::bernoulli_aperture::performance::finite_gas::GasHeatControl;
use crate::bernoulli_aperture::performance::gas_wall::{
    FiniteGasWallHeatProgram, GasWallConduction, GasWallHeatControl,
};
use crate::bernoulli_aperture::performance::thermal_wall::{
    FiniteWallHeatProgram, WallHeatControl,
};
use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, feedback};
use crate::dg_port::input::wall::regions::DgMaterialWallLayout;
use fs_exec::CancelGate;
use fs_matdb::{MaterialCard, PropertyKey, SelectionPolicy};
use fs_material::gas_species::input::GasRecord;
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_phs::linear_exact::{DenseExponentialBudget, MAX_PREPARED_EXPONENTIAL_DIM};
use fs_qty::Density;
use std::sync::Arc;

struct JointInput {
    conduction: GasWallConduction,
    exchange: HeatExchangeSolve,
    maximum_material_trials: usize,
}

pub(super) struct FiniteWallInput {
    line: usize,
    material: Arc<MaterialCard>,
    curve: EquilibriumEnthalpyPhaseCurve,
    initial: EquilibriumPhaseState,
    temperature_axis: String,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    program: FiniteWallHeatProgram,
    joint: Option<JointInput>,
}

fn exchange(row: &mut super::Row<'_>) -> Result<HeatExchangeSolve, PlateValveInputError> {
    let controls = HeatExchangeSolve {
        heat_bracket_j: [row.scalar()?, row.scalar()?],
        source_tolerance_j: row.scalar()?,
        energy_tolerance_j: row.scalar()?,
        maximum_iterations: row.count(128)?,
    };
    feedback::admit::<core::convert::Infallible>(controls).map_err(|_| {
        bad(
            row.line,
            "invalid signed heat/work bracket or feedback budget",
        )
    })?;
    Ok(controls)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl FiniteWallInput {
    pub(super) fn read(
        rows: &mut Reader<'_>,
        layout: Option<&DgMaterialWallLayout>,
        gate: &CancelGate,
    ) -> Result<Self, PlateValveInputError> {
        checkpoint(gate)?;
        let mut row = rows.row("finite_wall")?;
        let line = row.line;
        if row.word()? != "solid" {
            return Err(bad(
                line,
                "finite_wall requires the explicit fixed-density solid model",
            ));
        }
        row.finish()?;
        let layout = layout
            .ok_or_else(|| bad(line, "finite wall heat requires dg_wall material sources"))?;
        let mut records = layout
            .uniform_record()
            .into_iter()
            .chain(layout.regions().iter().map(|region| region.source()));
        let first = records
            .next()
            .ok_or_else(|| bad(line, "a rigid DG wall has no finite material caloric owner"))?;
        let card = first.pack().card();
        let context = first.wall().material().query().clone();
        // One uniform caloric body cannot silently average distinct source
        // conditions or materials. Thickness and mechanical restraint can vary.
        for record in records {
            checkpoint(gate)?;
            if record.pack().card().content_hash() != card.content_hash()
                || record.wall().material().query() != context
            {
                return Err(bad(
                    line,
                    "one finite wall body requires one original card and condition point across its main-bore regions",
                ));
            }
        }

        let mut row = rows.row("wall_caloric")?;
        let caloric_line = row.line;
        let heat_capacity = PropertyKey::new(row.word()?, SPECIFIC_HEAT_CAPACITY_DIMS);
        let density = PropertyKey::new(row.word()?, Density::DIMS);
        let temperature_axis = row.word()?.to_owned();
        let minimum_temperature_k = row.scalar()?;
        let maximum_temperature_k = row.scalar()?;
        let reference_enthalpy_j_kg = row.scalar()?;
        let interpolation_error_j_kg = row.scalar()?;
        if row.word()? != "single-claim" {
            return Err(bad(
                caloric_line,
                "wall caloric sources require explicit single-claim selection",
            ));
        }
        row.finish()?;
        // The shared authored DG wall grammar fixes T as absolute temperature;
        // omega and every other condition must remain unchanged.
        if temperature_axis != "T" || !context.axes().contains_key("T") {
            return Err(bad(
                caloric_line,
                "wall caloric temperature axis must be the original T axis",
            ));
        }
        let lower = context
            .clone()
            .with(&temperature_axis, minimum_temperature_k)
            .map_err(|_| bad(caloric_line, "invalid lower caloric temperature"))?;
        let upper = context
            .clone()
            .with(&temperature_axis, maximum_temperature_k)
            .map_err(|_| bad(caloric_line, "invalid upper caloric temperature"))?;
        let source = ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &heat_capacity,
                density: &density,
                temperature_axis: &temperature_axis,
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            reference_enthalpy_j_kg,
            interpolation_error_j_kg,
        )
        .map_err(|source| PlateValveInputError::Phase {
            line: caloric_line,
            source,
        })?;
        checkpoint(gate)?;

        let mut row = rows.row("wall_inventory")?;
        let inventory_line = row.line;
        let mass_kg = row.scalar()?;
        let volume_m3 = row.scalar()?;
        let surface_area_m2 = row.scalar()?;
        let initial_enthalpy_j_kg = row.scalar()?;
        row.finish()?;
        if ![mass_kg, volume_m3, surface_area_m2]
            .into_iter()
            .all(positive)
        {
            return Err(bad(
                inventory_line,
                "finite wall mass, volume and area must be positive",
            ));
        }
        let initial = source
            .curve()
            .state_at_specific_enthalpy(initial_enthalpy_j_kg)
            .map_err(|source| PlateValveInputError::Phase {
                line: inventory_line,
                source,
            })?;
        if context.axes()[&temperature_axis].to_bits() != initial.temperature_k().to_bits() {
            return Err(bad(
                inventory_line,
                "initial wall enthalpy must reconstruct the original mechanical wall temperature",
            ));
        }

        let mut row = rows.row("wall_feedback")?;
        let wall_exchange = exchange(&mut row)?;
        let maximum_absolute_parameter_work_j = row.scalar()?;
        if maximum_absolute_parameter_work_j < 0.0 {
            return Err(bad(
                row.line,
                "wall parameter-work allowance must be nonnegative",
            ));
        }
        row.finish()?;
        let mut row = rows.row("wall_dissipation")?;
        let dissipation_heat = match row.word()? {
            "none" => None,
            "native" => {
                let budget = DenseExponentialBudget {
                    max_dimension: row.count(MAX_PREPARED_EXPONENTIAL_DIM)?,
                    max_dense_bytes: row.parse()?,
                    max_multiply_terms: row.parse()?,
                };
                if budget.max_dimension == 0
                    || budget.max_dense_bytes == 0
                    || budget.max_multiply_terms == 0
                {
                    return Err(bad(
                        row.line,
                        "native wall heat requires positive dimension, byte and work budgets",
                    ));
                }
                Some(budget)
            }
            _ => {
                return Err(bad(
                    row.line,
                    "wall_dissipation requires none or native with explicit budgets",
                ));
            }
        };
        row.finish()?;
        let count = rows.count("wall_heat", MAX_CONTROLS)?;
        if count == 0 {
            return Err(bad(
                rows.line,
                "wall_heat requires an explicit initial power",
            ));
        }
        let mut controls = Vec::new();
        controls
            .try_reserve_exact(count)
            .map_err(|_| bad(rows.line, "cannot reserve wall heat controls"))?;
        let mut previous = None;
        for _ in 0..count {
            checkpoint(gate)?;
            let mut row = rows.row("wall_power")?;
            let sample = row.parse()?;
            let power_w = row.scalar()?;
            if previous.map_or(sample != 0, |last| sample <= last) {
                return Err(bad(
                    row.line,
                    "wall powers must start at zero and increase strictly",
                ));
            }
            row.finish()?;
            controls.push(WallHeatControl { sample, power_w });
            previous = Some(sample);
        }
        let joint = if rows
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next())
            == Some("gas_wall_conduction")
        {
            let mut row = rows.row("gas_wall_conduction")?;
            let conductivity = PropertyKey::new(row.word()?, THERMAL_CONDUCTIVITY_DIMS);
            let area_m2 = row.scalar()?;
            let gas_distance_m = row.scalar()?;
            let body_distance_m = row.scalar()?;
            if area_m2 < 0.0 || !positive(gas_distance_m) || !positive(body_distance_m) {
                return Err(bad(
                    row.line,
                    "gas/wall conduction requires nonnegative area and positive distances",
                ));
            }
            row.finish()?;
            let mut row = rows.row("gas_wall_feedback")?;
            let exchange = exchange(&mut row)?;
            let maximum_material_trials = row.count(1_000_000)?;
            if maximum_material_trials == 0 {
                return Err(bad(
                    row.line,
                    "gas/wall feedback requires a positive total material-trial budget",
                ));
            }
            row.finish()?;
            Some(JointInput {
                conduction: GasWallConduction {
                    context,
                    conductivity,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    liquid_fraction_axis: None,
                    area_m2,
                    gas_distance_m,
                    body_distance_m,
                },
                exchange,
                maximum_material_trials,
            })
        } else {
            None
        };
        checkpoint(gate)?;
        Ok(Self {
            line,
            material: Arc::new(card.clone()),
            curve: source.curve().clone(),
            initial,
            temperature_axis,
            mass_kg,
            volume_m3,
            surface_area_m2,
            program: FiniteWallHeatProgram {
                controls,
                coupling: FiniteWallThermalSolve {
                    exchange: wall_exchange,
                    maximum_absolute_parameter_work_j,
                    dissipation_heat,
                },
            },
            joint,
        })
    }

    pub(super) fn requires_gas(&self) -> bool {
        self.joint.is_some()
    }

    pub(super) fn bind(
        self,
        renderer: AperturePerformance,
        gas: Option<FiniteGasInput>,
        gate: &CancelGate,
        records: &mut Vec<GasRecord>,
    ) -> Result<AperturePerformance, PlateValveInputError> {
        checkpoint(gate)?;
        let body = UniformEnthalpyStepInput {
            curve: &self.curve,
            initial: self.initial,
            mass_kg: self.mass_kg,
            volume_m3: self.volume_m3,
            surface_area_m2: self.surface_area_m2,
            internal_heat_j: 0.0,
            duration_s: renderer.system().aperture().spec().time_step_s,
        };
        let playback_gate = Arc::new(CancelGate::new_clock_free());
        let renderer = match (gas, self.joint) {
            (None, None) => renderer
                .with_finite_wall_heat(
                    body,
                    self.material,
                    self.temperature_axis,
                    self.program,
                    playback_gate,
                )
                .map_err(PlateValveInputError::Gesture)?,
            (Some(gas), Some(joint)) => {
                let (inventory, gas_program) = gas.into_inventory(&renderer, gate)?;
                let initial = GasRecord::Thermochemical(inventory.state().clone());
                let controls = merge_powers(
                    &renderer,
                    &gas_program.controls,
                    &self.program.controls,
                    self.line,
                    gate,
                )?;
                let program = FiniteGasWallHeatProgram {
                    controls,
                    coupling: FiniteGasWallThermalSolve {
                        gas: gas_program.coupling,
                        caloric: gas_program.caloric,
                        wall: self.program.coupling,
                        exchange: joint.exchange,
                        maximum_material_trials: joint.maximum_material_trials,
                    },
                };
                let renderer = renderer
                    .with_finite_gas_wall_heat(
                        inventory,
                        body,
                        self.material,
                        self.temperature_axis,
                        joint.conduction,
                        program,
                        playback_gate,
                    )
                    .map_err(PlateValveInputError::Gesture)?;
                records.push(initial);
                renderer
            }
            _ => {
                return Err(bad(
                    self.line,
                    "finite gas and wall heat require their explicit reciprocal conduction; conduction also requires finite gas",
                ));
            }
        };
        checkpoint(gate)?;
        Ok(renderer)
    }
}

/// Validate each original list before merging; sorting could silently repair an
/// invalid schedule. Every union sample retains both actual held powers.
fn merge_powers(
    renderer: &AperturePerformance,
    gas: &[GasHeatControl],
    wall: &[WallHeatControl],
    line: usize,
    gate: &CancelGate,
) -> Result<Vec<GasWallHeatControl>, PlateValveInputError> {
    let samples = renderer.config().samples;
    let dt = renderer.system().aperture().spec().time_step_s;
    validate_powers(
        gas.iter().map(|event| (event.sample, event.power_w)),
        samples,
        dt,
        line,
        gate,
    )?;
    validate_powers(
        wall.iter().map(|event| (event.sample, event.power_w)),
        samples,
        dt,
        line,
        gate,
    )?;
    let (assignments, work) = renderer
        .remaining_control_budgets()
        .map_err(PlateValveInputError::Gesture)?;
    let capacity = (assignments / 2)
        .min(usize::try_from(work / 2).unwrap_or(usize::MAX))
        .min(MAX_CONTROLS / 2);
    if gas.len().max(wall.len()) > capacity {
        return Err(bad(
            line,
            "paired gas/wall powers exceed the shared control or work budget",
        ));
    }
    let mut output = Vec::new();
    output
        .try_reserve_exact(gas.len().saturating_add(wall.len()).min(capacity))
        .map_err(|_| bad(line, "cannot reserve paired thermal powers"))?;
    let (mut g, mut w) = (0, 0);
    let mut held = [0.0; 2];
    while g < gas.len() || w < wall.len() {
        checkpoint(gate)?;
        if output.len() == capacity {
            return Err(bad(
                line,
                "union of thermal event samples exceeds paired control capacity",
            ));
        }
        let sample = match (gas.get(g), wall.get(w)) {
            (Some(a), Some(b)) => a.sample.min(b.sample),
            (Some(a), None) => a.sample,
            (None, Some(b)) => b.sample,
            (None, None) => unreachable!(),
        };
        if gas.get(g).is_some_and(|event| event.sample == sample) {
            held[0] = gas[g].power_w;
            g += 1;
        }
        if wall.get(w).is_some_and(|event| event.sample == sample) {
            held[1] = wall[w].power_w;
            w += 1;
        }
        output.push(GasWallHeatControl {
            sample,
            gas_power_w: held[0],
            wall_power_w: held[1],
        });
    }
    Ok(output)
}

fn validate_powers(
    events: impl Iterator<Item = (u64, f64)>,
    samples: u64,
    dt: f64,
    line: usize,
    gate: &CancelGate,
) -> Result<(), PlateValveInputError> {
    let mut previous = None;
    for (sample, power) in events {
        checkpoint(gate)?;
        let heat = power * dt;
        if sample >= samples
            || previous.map_or(sample != 0, |last| sample <= last)
            || !power.is_finite()
            || !heat.is_finite()
            || (power != 0.0 && heat == 0.0)
        {
            return Err(bad(
                line,
                "thermal powers need increasing in-window samples and representable signed energy",
            ));
        }
        previous = Some(sample);
    }
    if previous.is_none() {
        return Err(bad(
            line,
            "both thermal inventories need initial power assignments",
        ));
    }
    Ok(())
}
