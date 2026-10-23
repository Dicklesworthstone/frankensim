//! Physical requirements for the standalone CSV scan-fitting consumer.
//!
//! Parsing introduces no numerical model or penalty objective. Absolute times
//! are inserted before measurement quadrature is built, then lowered to the
//! existing scan constraints on that exact grid. The program stays borrowed so
//! source/knot identities cannot change between parsing and binding.
use crate::{heaters::ScanHeatingProgram, specimen::Result};
use fs_ascent::scan_design::constraints::{ScanActuatorBounds, ScanControlBound,
    ScanEnergyCap, ScanStateBound, ScanPeakBound, ScanPowerSlewBound, ScanPathSpeedBound,
    ScanResponse, ScanResponseBudget, ScanPeakResponse};
use fs_conduction::ConductionMesh;
use fs_conduction::bc::history::march::moving::motion::velocity::MovingHeatPathBudget;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::ScanObservable;
use fs_exec::Cx;
use std::collections::BTreeSet;

/// Independent input, dense-optimizer and response-work allowances.
#[derive(Debug, Clone, Copy)]
pub struct RequirementBudget {
    pub max_bytes: usize,
    pub max_declarations: usize,
    pub max_kkt_dimension: usize,
    pub max_peak_samples: usize,
    pub max_gradient_entries: usize,
    pub max_transpose_solves: usize,
}

#[derive(Clone)]
enum Declaration {
    Power(ScanControlBound),
    Energy(ScanEnergyCap),
    State { time: f64, vertex: usize, quantity: ScanObservable,
        lower: Option<f64>, upper: Option<f64>, scale: f64 },
    Peak { start: f64, end: f64, quantity: ScanObservable,
        maximum: f64, scale: f64, smoothing: f64, vertices: Vec<usize> },
    Slew(ScanPowerSlewBound),
    Speed(ScanPathSpeedBound),
}

/// Declarations bound to an immutable heater program and original mesh.
pub struct ScanRequirements<'a> {
    program: &'a ScanHeatingProgram,
    mesh: &'a ConductionMesh,
    declarations: Vec<Declaration>,
    times_s: Vec<f64>,
    budget: RequirementBudget,
}

/// Native declarations and their admitted exact-grid work/storage shapes.
/// Values and gradients are evaluated by ConstrainedScanStudy, not this loader.
pub struct BoundScanRequirements {
    pub powers: Vec<ScanControlBound>,
    pub energy: Vec<ScanEnergyCap>,
    pub states: Vec<ScanStateBound>,
    pub peaks: Vec<ScanPeakBound>,
    pub actuators: ScanActuatorBounds,
    pub response_budget: ScanResponseBudget,
    pub max_peak_samples: usize,
    pub kkt_dimension: usize,
}

impl<'a> ScanRequirements<'a> {
    /// No extra physical requirements; the native nonnegative-power rows remain.
    pub fn none(program: &'a ScanHeatingProgram, mesh: &'a ConductionMesh,
        budget: RequirementBudget) -> Self
    {
        Self { program, mesh, declarations: Vec::new(), times_s: Vec::new(), budget }
    }

    /// Read FS_SCAN_REQUIREMENTS_V1. All limits/scales have explicit SI units;
    /// `*` denotes an absent lower/upper side, never an infinite number.
    /// Fixed power values are not controls: use an energy cap to include all
    /// heaters, or change the heater declaration intentionally before fitting.
    pub fn from_text(cx: &Cx<'_>, program: &'a ScanHeatingProgram, mesh: &'a ConductionMesh,
        text: &str, budget: RequirementBudget) -> Result<Self>
    {
        poll(cx)?;
        if text.len() > budget.max_bytes { return Err("requirements exceed their byte allowance".into()); }
        let mut lines = text.lines().enumerate().filter_map(|(i, line)| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i + 1, line))
        });
        if lines.next().map(|(_, line)| line) != Some("FS_SCAN_REQUIREMENTS_V1") {
            return Err("requirements must begin with FS_SCAN_REQUIREMENTS_V1".into());
        }
        let mut out = Self::none(program, mesh, budget);
        let mut selected = BTreeSet::new();
        let mut slew_sources = BTreeSet::new();
        let mut speed_sources = BTreeSet::new();
        for (line_number, line) in lines {
            poll(cx)?;
            let result = (|| -> Result<Declaration> {
                if out.declarations.len() >= budget.max_declarations {
                    return Err("physical requirement declaration allowance exceeded".into());
                }
                let fields: Vec<_> = line.split_whitespace().take(10).collect();
                let source = |text: &str| -> Result<usize> {
                    let id: usize = text.parse()?;
                    if id >= program.sources().len() { return Err("requirement names an absent heater".into()); }
                    Ok(id)
                };
                let time = |text: &str, zero: bool| -> Result<f64> {
                    let value = number(text)?;
                    if value < 0.0 || (!zero && value == 0.0) || value > program.end_time_s() {
                        return Err("requirement time is outside the physical observation horizon".into());
                    }
                    Ok(value)
                };
                Ok(match fields.as_slice() {
                    ["power", id, knot, lo, hi, scale] => {
                        let id = source(id)?;
                        let knot: usize = knot.parse()?;
                        let coordinate = program.controls().iter().position(|c| c.source == id && c.value == knot)
                            .ok_or("power bound must name a selected heater/knot, not a fixed value")?;
                        if !selected.insert(coordinate) { return Err("duplicate power bound; combine its two sides".into()); }
                        let (lower, upper, scale) = interval(lo, hi, scale)?;
                        if lower.into_iter().chain(upper).any(|v| v < 0.0) {
                            return Err("heater power limits must be nonnegative watts".into());
                        }
                        Declaration::Power(ScanControlBound { coordinate, lower, upper, scale })
                    }
                    ["energy", start, end, maximum, scale] => {
                        let start = time(start, true)?; let end = time(end, false)?;
                        if end <= start { return Err("energy window needs start < end".into()); }
                        let maximum_j = nonnegative(maximum)?;
                        Declaration::Energy(ScanEnergyCap { start_time_s: start, end_time_s: end,
                            maximum_j, scale_j: positive(scale)? })
                    }
                    ["state", at, vertex, observable, lo, hi, scale] => {
                        let time = time(at, false)?; let vertex: usize = vertex.parse()?;
                        if vertex >= mesh.vertex_count() { return Err("state bound names an absent mesh vertex".into()); }
                        let quantity = quantity(observable)?;
                        let (lower, upper, scale) = interval(lo, hi, scale)?;
                        if quantity == ScanObservable::LiquidMassFraction
                            && lower.into_iter().chain(upper).any(|v| !(0.0..=1.0).contains(&v)) {
                            return Err("state fraction limits must be in [0,1], not percentages".into());
                        }
                        Declaration::State { time, vertex, quantity, lower, upper, scale }
                    }
                    ["peak", start, end, observable, maximum, scale, smoothing, selection] => {
                        let start = time(start, false)?; let end = time(end, false)?;
                        if end < start { return Err("peak window needs start <= end at accepted endpoints".into()); }
                        let vertices = vertices(cx, selection, mesh.vertex_count(), budget.max_peak_samples)?;
                        Declaration::Peak { start, end, quantity: quantity(observable)?,
                            maximum: number(maximum)?, scale: positive(scale)?,
                            smoothing: positive(smoothing)?, vertices }
                    }
                    ["slew", id, rise, fall, scale] => {
                        let source = source(id)?;
                        if !slew_sources.insert(source) { return Err("duplicate heater slew requirement".into()); }
                        Declaration::Slew(ScanPowerSlewBound { source,
                            maximum_rise_w_s: nonnegative(rise)?, maximum_fall_w_s: nonnegative(fall)?,
                            scale_w_s: positive(scale)? })
                    }
                    ["speed", id, maximum, scale] => {
                        let source = source(id)?;
                        if !speed_sources.insert(source) { return Err("duplicate heater speed requirement".into()); }
                        Declaration::Speed(ScanPathSpeedBound { source,
                            maximum_m_s: nonnegative(maximum)?, scale_m_s: positive(scale)? })
                    }
                    _ => return Err("expected power, energy, state, peak, slew or speed declaration with its documented fields".into()),
                })
            })();
            let declaration = result.map_err(|e| format!("requirement line {line_number}: {e}"))?;
            match &declaration {
                Declaration::State { time, .. } => out.times_s.push(*time),
                Declaration::Peak { start, end, .. } => out.times_s.extend([*start, *end]),
                _ => (),
            }
            out.declarations.push(declaration);
        }
        if out.declarations.is_empty() { return Err("an explicit requirements file must contain a requirement".into()); }
        out.times_s.sort_by(f64::total_cmp);
        out.times_s.dedup();
        poll(cx)?;
        Ok(out)
    }

    /// Merge these exact times with measurement/source events BEFORE refinement
    /// and construction of temporal readouts. Energy is integrated analytically
    /// from source schedules, so its window does not require a physical stop.
    pub fn event_times_s(&self) -> &[f64] { &self.times_s }

    /// Lower onto the final simulation grid, without constructing a physical run.
    /// Every refined endpoint in a peak window participates. Budgets cover those
    /// samples and the additional original adjoint sweeps, not just file rows.
    pub fn bind(&self, cx: &Cx<'_>, endpoints: &[f64]) -> Result<BoundScanRequirements> {
        poll(cx)?;
        let mut previous = 0.0;
        for &time in endpoints {
            poll(cx)?;
            if !time.is_finite() || time <= previous || !(time - previous).is_finite() {
                return Err("constraint binding needs a finite increasing physical grid".into());
            }
            previous = time;
        }
        if endpoints.is_empty() || previous.to_bits() != self.program.end_time_s().to_bits() {
            return Err("constraint grid must end at the original heater-program horizon".into());
        }
        let endpoint = |time: f64| -> Result<usize> {
            endpoints.binary_search_by(|t| t.total_cmp(&time)).map_err(|_| "physical grid omits an exact constraint time".into())
        };
        let mut out = BoundScanRequirements {
            powers: Vec::new(), energy: Vec::new(), states: Vec::new(), peaks: Vec::new(),
            actuators: ScanActuatorBounds::default(),
            response_budget: ScanResponseBudget { max_responses: 0, max_gradient_entries: 0,
                max_transpose_solves: endpoints.len() }, max_peak_samples: 0, kkt_dimension: 0,
        };
        let mut rows = self.program.controls().len();
        for declaration in &self.declarations {
            poll(cx)?;
            let added = match declaration {
                Declaration::Power(bound) => {
                    out.powers.push(*bound);
                    usize::from(bound.lower.is_some()) + usize::from(bound.upper.is_some())
                }
                Declaration::Energy(cap) => { out.energy.push(*cap); 1 }
                Declaration::State { time, vertex, quantity, lower, upper, scale } => {
                    let at = endpoint(*time)?;
                    out.states.push(ScanStateBound { response: ScanResponse {
                        endpoint: at, vertex: *vertex, observable: *quantity }, lower: *lower, upper: *upper, scale: *scale });
                    out.response_budget.max_transpose_solves = add(out.response_budget.max_transpose_solves, at + 1)?;
                    usize::from(lower.is_some()) + usize::from(upper.is_some())
                }
                Declaration::Peak { start, end, quantity, maximum, scale, smoothing, vertices } => {
                    let first = endpoint(*start)?; let last = endpoint(*end)?;
                    let count = (last - first + 1).checked_mul(vertices.len()).ok_or("peak sample count overflow")?;
                    out.max_peak_samples = add(out.max_peak_samples, count)?;
                    if out.max_peak_samples > self.budget.max_peak_samples { return Err("refined peak samples exceed their allowance".into()); }
                    out.peaks.push(ScanPeakBound { response: ScanPeakResponse {
                        first_endpoint: first, last_endpoint: last, vertices: vertices.clone(),
                        observable: *quantity, smoothing_scale: *smoothing }, maximum: *maximum, scale: *scale });
                    out.response_budget.max_transpose_solves = add(out.response_budget.max_transpose_solves, last + 1)?;
                    1
                }
                Declaration::Slew(bound) => {
                    let cap = self.budget.max_kkt_dimension.saturating_sub(rows) / 2;
                    let count = self.program.sources()[bound.source].absorbed_power_slew_segments(cx, cap)?.len();
                    out.actuators.max_power_segments = add(out.actuators.max_power_segments, count)?;
                    out.actuators.power_slew.push(*bound);
                    count.checked_mul(2).ok_or("slew row count overflow")?
                }
                Declaration::Speed(bound) => {
                    let max_knots = self.program.authored_control_count();
                    let count = self.program.sources()[bound.source].path_velocity_segments(cx, MovingHeatPathBudget {
                        max_knots, max_segments: self.budget.max_kkt_dimension.saturating_sub(rows) })?.len();
                    out.actuators.path_budget.max_knots_per_source = max_knots;
                    out.actuators.path_budget.max_segments = add(out.actuators.path_budget.max_segments, count)?;
                    out.actuators.path_speed.push(*bound);
                    count
                }
            };
            rows = add(rows, added)?;
            out.kkt_dimension = add(self.program.controls().len(), rows)?;
            if out.kkt_dimension > self.budget.max_kkt_dimension { return Err("physical requirements exceed the dense SQP allowance".into()); }
        }
        out.kkt_dimension = add(self.program.controls().len(), rows)?;
        out.response_budget.max_responses = add(out.states.len(), out.peaks.len())?;
        out.response_budget.max_gradient_entries = add(self.mesh.vertex_count(), self.program.authored_control_count())?
            .checked_mul(out.response_budget.max_responses).ok_or("response gradient shape overflow")?;
        if out.kkt_dimension > self.budget.max_kkt_dimension
            || out.response_budget.max_gradient_entries > self.budget.max_gradient_entries
            || out.response_budget.max_transpose_solves > self.budget.max_transpose_solves {
            return Err("constraint matrix, response gradients or transpose work exceed their allowances".into());
        }
        poll(cx)?;
        Ok(out)
    }
}

fn number(text: &str) -> Result<f64> {
    let v: f64 = text.parse()?;
    if !v.is_finite() { return Err("requirement numbers must be finite physical values".into()); }
    Ok(v)
}
fn positive(text: &str) -> Result<f64> {
    let v = number(text)?;
    if v <= 0.0 { return Err("requirement scales must be positive".into()); }
    Ok(v)
}
fn nonnegative(text: &str) -> Result<f64> {
    let v = number(text)?;
    if v < 0.0 { return Err("requirement must be nonnegative".into()); }
    Ok(v)
}
fn interval(lo: &str, hi: &str, scale: &str) -> Result<(Option<f64>, Option<f64>, f64)> {
    let lower = if lo == "*" { None } else { Some(number(lo)?) };
    let upper = if hi == "*" { None } else { Some(number(hi)?) };
    if (lower.is_none() && upper.is_none()) || lower.zip(upper).is_some_and(|(a,b)| a >= b) {
        return Err("bound needs at least one side and strictly ordered finite limits".into());
    }
    Ok((lower, upper, positive(scale)?))
}
fn quantity(text: &str) -> Result<ScanObservable> {
    match text {
        "temperature_k" => Ok(ScanObservable::Temperature),
        "specific_enthalpy_j_kg" => Ok(ScanObservable::SpecificEnthalpy),
        "liquid_mass_fraction" => Ok(ScanObservable::LiquidMassFraction),
        _ => Err("quantity must name temperature_k, specific_enthalpy_j_kg or liquid_mass_fraction".into()),
    }
}
fn vertices(cx: &Cx<'_>, text: &str, count: usize, maximum: usize) -> Result<Vec<usize>> {
    if text == "all" {
        if count == 0 || count > maximum { return Err("peak vertex selection exceeds its allowance".into()); }
        return Ok((0..count).collect());
    }
    let mut out = Vec::new();
    for item in text.split(';') {
        poll(cx)?;
        let v: usize = item.parse()?;
        if out.len() >= maximum || v >= count || out.last().is_some_and(|&old| old >= v) {
            return Err("peak vertices must be increasing distinct mesh IDs within their allowance".into());
        }
        out.push(v);
    }
    Ok(out)
}
fn add(a: usize, b: usize) -> Result<usize> { a.checked_add(b).ok_or_else(|| "constraint work/shape count overflow".into()) }
fn poll(cx: &Cx<'_>) -> Result<()> { cx.checkpoint().map_err(|_| "constraint loading cancelled".into()) }

#[cfg(test)]
mod tests;
