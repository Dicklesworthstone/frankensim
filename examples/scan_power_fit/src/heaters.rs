//! Declared multi-heater schedules for the standalone measurement-fitting consumer.
//!
//! Powers marked `fit` become existing nonnegative-linear scan controls. Everything
//! else stays fixed, including unselected heaters and exact zero shutdowns. This
//! file format is not a second solver or a general FrankenScript interface.
use crate::specimen::{Result, ScanSpecimen};
use fs_ascent::scan_design::ScanPowerControl;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::{ConductionMesh, COMPONENT_POWER_DIMS};
use fs_exec::Cx;
use fs_qty::Dims;
use fs_scenario::signal::{Interp, TimeSignal};

/// Retained input-shape allowances, not total memory or accuracy bounds.
#[derive(Debug, Clone, Copy)]
pub struct HeatingProgramBudget {
    pub max_bytes: usize,
    pub max_sources: usize,
    pub max_knots: usize,
    pub max_fitted_powers: usize,
    pub max_authored_controls: usize,
}

/// One complete exported source schedule, including all fixed/off values.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatingSchedule {
    pub source: usize,
    pub interpolation: Interp,
    pub times_s: Vec<f64>,
    pub powers_w: Vec<f64>,
    pub fitted: Vec<bool>,
}

/// Immutable source templates, control identities and their exact event union.
/// Every source is clamped to zero after its last fixed knot. The observation
/// horizon can be later than shutdown, so post-heating cooling is retained.
pub struct ScanHeatingProgram {
    sources: Vec<MovingHeatSource>,
    controls: Vec<ScanPowerControl>,
    schedules: Vec<HeatingSchedule>,
    events_s: Vec<f64>,
    end_time_s: f64,
    authored_controls: usize,
}

struct DeclaredSource {
    interpolation: Interp,
    half_width_m: [f64; 3],
    times: Vec<f64>,
    centers: Vec<[f64; 3]>,
    powers: Vec<f64>,
    references: Vec<Option<f64>>,
}

impl ScanHeatingProgram {
    /// Read FS_SCAN_HEATERS_V1. Heater IDs are contiguous declaration-order
    /// indices; knot records may interleave sources but must increase in time
    /// within each source. All positions/widths are metres and all powers watts.
    /// Each source starts at zero seconds and ends at a FIXED zero power no later
    /// than end_time_s. A linear source may start on; zero startup is explicit.
    /// There must be at least one fitted value across the complete program.
    pub fn from_text(cx: &Cx<'_>, mesh: &ConductionMesh, end_time_s: f64,
        text: &str, budget: HeatingProgramBudget) -> Result<Self>
    {
        poll(cx)?;
        if text.len() > budget.max_bytes || !end_time_s.is_finite() || end_time_s <= 0.0 {
            return Err("heater input needs a bounded file and finite positive observation horizon".into());
        }
        let mut lines = text.lines().enumerate().filter_map(|(i, line)| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i + 1, line))
        });
        if lines.next().map(|(_, line)| line) != Some("FS_SCAN_HEATERS_V1") {
            return Err("heater input must begin with FS_SCAN_HEATERS_V1".into());
        }
        let mut declared = Vec::<DeclaredSource>::new();
        let mut total_knots = 0usize;
        let mut selected = 0usize;
        for (line_number, line) in lines {
            poll(cx)?;
            let result = (|| -> Result<()> {
                let fields: Vec<_> = line.split_whitespace().take(10).collect();
                match fields.first().copied() {
                    Some("heater") if fields.len() == 6 => {
                        let id: usize = fields[1].parse()?;
                        if id != declared.len() || declared.len() >= budget.max_sources {
                            return Err("heater IDs must be contiguous from zero within the source allowance".into());
                        }
                        let interpolation = match fields[2] {
                            "hold" => Interp::Hold,
                            "linear" => Interp::Linear,
                            _ => return Err("heater interpolation must be hold or linear".into()),
                        };
                        let half_width_m = [number(fields[3])?, number(fields[4])?, number(fields[5])?];
                        if half_width_m.iter().any(|&x| x <= 0.0) {
                            return Err("heater half widths must be positive metres".into());
                        }
                        declared.push(DeclaredSource { interpolation, half_width_m,
                            times: Vec::new(), centers: Vec::new(), powers: Vec::new(), references: Vec::new() });
                    }
                    Some("knot") if matches!(fields.len(), 8 | 9) => {
                        let id: usize = fields[1].parse()?;
                        let source = declared.get_mut(id).ok_or("knot names an undeclared heater")?;
                        let time = number(fields[2])?;
                        let center = [number(fields[3])?, number(fields[4])?, number(fields[5])?];
                        let power = number(fields[6])?;
                        let reference = match (fields[7], fields.len()) {
                            ("fixed", 8) => None,
                            ("fit", 9) => {
                                let scale = number(fields[8])?;
                                if scale <= 0.0 || !(power / scale).is_finite()
                                    || (power > 0.0 && power / scale == 0.0) {
                                    return Err("fitted powers need positive watt scales and representable coordinates".into());
                                }
                                if selected >= budget.max_fitted_powers {
                                    return Err("fitted powers exceed their explicit dimension allowance".into());
                                }
                                Some(scale)
                            }
                            _ => return Err("knot ends with fixed or fit REFERENCE_W".into()),
                        };
                        if power < 0.0 || time < 0.0 || time > end_time_s
                            || source.times.last().is_some_and(|&old| time <= old || !(time - old).is_finite())
                            || (source.times.is_empty() && time != 0.0)
                        { return Err("knots require nonnegative watts and increasing seconds from zero through the observation horizon".into()); }
                        if total_knots >= budget.max_knots { return Err("heater knots exceed their combined allowance".into()); }
                        total_knots += 1;
                        selected += usize::from(reference.is_some());
                        source.times.push(time);
                        source.centers.push(center);
                        source.powers.push(power);
                        source.references.push(reference);
                    }
                    _ => return Err("expected heater ID hold|linear HX HY HZ or knot ID TIME X Y Z W fixed|fit [REFERENCE_W]".into()),
                }
                Ok(())
            })();
            result.map_err(|error| format!("heater line {line_number}: {error}"))?;
        }
        if declared.is_empty() || selected == 0 { return Err("heater program needs sources and at least one fitted power".into()); }
        // Each table contributes all position/power values, not just the fitted
        // ones. Admit this before constructing the original source objects.
        let authored = total_knots.checked_mul(4).and_then(|n|
            declared.len().checked_mul(3).and_then(|widths| n.checked_add(widths)))
            .filter(|&n| n <= budget.max_authored_controls)
            .ok_or("heater declarations exceed their combined authored-control allowance")?;
        let mut sources = Vec::with_capacity(declared.len());
        let mut controls = Vec::with_capacity(selected);
        let mut schedules = Vec::with_capacity(declared.len());
        let mut events_s = Vec::with_capacity(total_knots + 1);
        for (id, source) in declared.into_iter().enumerate() {
            poll(cx)?;
            if source.times.len() < 2 || source.powers.last() != Some(&0.0)
                || source.references.last().is_none_or(|r| r.is_some())
            { return Err(format!("heater {id} needs at least two knots and a fixed zero-power final knot").into()); }
            for (value, reference) in source.references.iter().enumerate() {
                if let Some(reference_w) = reference {
                    controls.push(ScanPowerControl { source: id, value, reference_w: *reference_w });
                }
            }
            let center = std::array::from_fn(|axis| TimeSignal::Table {
                times: source.times.clone(), values: source.centers.iter().map(|p| p[axis]).collect(),
                dims: Dims([1, 0, 0, 0, 0, 0]), interp: source.interpolation,
            });
            sources.push(MovingHeatSource::new(center, source.half_width_m,
                TimeSignal::Table { times: source.times.clone(), values: source.powers.clone(),
                    dims: COMPONENT_POWER_DIMS, interp: source.interpolation }, 0.0,
                MovingHeatBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count(),
                    min_active_vertices: 1, power_relative_tolerance: 1e-11 })?);
            events_s.extend(source.times.iter().copied().filter(|&t| t > 0.0));
            schedules.push(HeatingSchedule { source: id, interpolation: source.interpolation,
                times_s: source.times, powers_w: source.powers,
                fitted: source.references.iter().map(Option::is_some).collect() });
        }
        events_s.push(end_time_s);
        events_s.sort_by(f64::total_cmp);
        events_s.dedup();
        poll(cx)?;
        Ok(Self { sources, controls, schedules, events_s, end_time_s, authored_controls: authored })
    }

    /// Preserve the exact historical single-source/two-power declaration.
    pub fn from_specimen(cx: &Cx<'_>, specimen: &ScanSpecimen, initial_w: [f64; 2]) -> Result<Self> {
        poll(cx)?;
        let source = specimen.source(initial_w)?;
        let authored_controls = source.schedule_control_count()?;
        Ok(Self { sources: vec![source], controls: vec![
                ScanPowerControl { source: 0, value: 0, reference_w: 1.0 },
                ScanPowerControl { source: 0, value: 1, reference_w: 1.0 }],
            schedules: vec![HeatingSchedule { source: 0, interpolation: Interp::Hold,
                times_s: vec![0.0, specimen.switch_time_s, specimen.end_time_s],
                powers_w: vec![initial_w[0], initial_w[1], 0.0], fitted: vec![true, true, false] }],
            events_s: vec![specimen.switch_time_s, specimen.end_time_s],
            end_time_s: specimen.end_time_s, authored_controls })
    }

    pub fn sources(&self) -> &[MovingHeatSource] { &self.sources }
    pub fn controls(&self) -> &[ScanPowerControl] { &self.controls }
    pub fn event_times_s(&self) -> &[f64] { &self.events_s }
    pub fn end_time_s(&self) -> f64 { self.end_time_s }
    pub fn authored_control_count(&self) -> usize { self.authored_controls }

    /// Check every unselected value, original path, clock, interpolation and
    /// footprint before presenting a fitted schedule. No clipping or rounding
    /// to zero occurs. This check alone asserts neither fit nor feasibility.
    pub fn export(&self, cx: &Cx<'_>, sources: &[MovingHeatSource]) -> Result<Vec<HeatingSchedule>> {
        poll(cx)?;
        if sources.len() != self.sources.len() { return Err("export changed the heater count".into()); }
        let mut out = Vec::with_capacity(sources.len());
        for (i, ((source, template), schedule)) in sources.iter().zip(&self.sources).zip(&self.schedules).enumerate() {
            poll(cx)?;
            let powers = source.absorbed_power_values(schedule.times_s.len())?;
            if powers.len() != schedule.powers_w.len() { return Err("export changed the heater knot count".into()); }
            for ((&power, &fixed), &fitted) in powers.iter().zip(&schedule.powers_w).zip(&schedule.fitted) {
                if !fitted && power.to_bits() != fixed.to_bits() { return Err("export changed a fixed heater power or shutdown".into()); }
            }
            if template.with_absorbed_power_values(&powers, 0.0, powers.len())? != *source {
                return Err("export changed a heater path, clock, footprint or sampling policy".into());
            }
            out.push(HeatingSchedule { source: i, interpolation: schedule.interpolation,
                times_s: schedule.times_s.clone(), powers_w: powers, fitted: schedule.fitted.clone() });
        }
        poll(cx)?;
        Ok(out)
    }
}

fn number(text: &str) -> Result<f64> {
    let value: f64 = text.parse()?;
    if !value.is_finite() { return Err("heater numbers must be finite SI values".into()); }
    Ok(value)
}
fn poll(cx: &Cx<'_>) -> Result<()> { cx.checkpoint().map_err(|_| "heater program cancelled".into()) }

#[cfg(test)]
mod tests;
