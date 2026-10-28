//! Contract the accepted trajectory's absolute-watt controls with the ORIGINAL
//! scheduled workload. dJ/da=sum_i (dJ/dP_i) P_original,i. The physical adjoint
//! already owns dt and all storage/boundary feedback; neither dt nor duty is
//! applied here. Original watts also make the derivative well-defined at a=0.
use fs_project::{ProjectSpec, TransientRegionPower};
use super::{J, Result, array, invalid, named, number, policy};

pub(super) fn coefficient(
    report: &J, original: &ProjectSpec, sampled: &ProjectSpec,
    region: &str, amplitude: f64, prefix_steps: usize,
) -> Result<f64> {
    let original_time = policy(original).ok_or_else(|| invalid("missing original workload time policy"))?;
    let sampled_time = policy(sampled).ok_or_else(|| invalid("missing sampled workload time policy"))?;
    let original_grid = original_time.coarse_step_ends_s()
        .map_err(|e| invalid(format!("{}: {}", e.code, e.what)))?;
    let sampled_grid = sampled_time.coarse_step_ends_s()
        .map_err(|e| invalid(format!("{}: {}", e.code, e.what)))?;
    if original_grid.len() != sampled_grid.len()
        || original_grid.iter().zip(&sampled_grid).any(|(a, b)| a.to_bits() != b.to_bits())
    {
        return Err(invalid("schedule amplitude cannot change the original event-aligned time grid"));
    }
    let find = |schedules: &[TransientRegionPower]| -> Result<usize> {
        let mut found = schedules.iter().enumerate().filter(|(_, s)| s.region == region);
        let (index, _) = found.next().ok_or_else(|| invalid("missing amplitude-controlled schedule"))?;
        if found.next().is_some() { return Err(invalid("ambiguous amplitude-controlled schedule")); }
        Ok(index)
    };
    let original = &original_time.power_schedules[find(&original_time.power_schedules)?];
    let sampled = &sampled_time.power_schedules[find(&sampled_time.power_schedules)?];
    contract(report, original, sampled, &original_grid, amplitude, prefix_steps)
}

fn contract(
    report: &J, original: &TransientRegionPower, sampled: &TransientRegionPower,
    coarse: &[f64], amplitude: f64, prefix_steps: usize,
) -> Result<f64> {
    let count = coarse.len().checked_mul(2)
        .filter(|&count| count > 0 && count <= 10_000)
        .ok_or_else(|| invalid("schedule adjoint grid exceeds the native step allowance"))?;
    if !amplitude.is_finite() || amplitude < 0.0 || prefix_steps > count
        || original.region != sampled.region || original.source != sampled.source
        || original.steps.len() != sampled.steps.len() || original.steps.is_empty()
    {
        return Err(invalid("schedule amplitude lost its original region, source, shape or selected prefix"));
    }
    let mut start = 0.0;
    for (old, actual) in original.steps.iter().zip(&sampled.steps) {
        let watts = old.watts.value;
        let scaled = watts * amplitude;
        if !watts.is_finite() || watts < 0.0 || !scaled.is_finite()
            || (watts > 0.0 && amplitude > 0.0 && scaled == 0.0)
            || actual.watts.value.to_bits() != scaled.to_bits()
            || old.watts.dims != fs_project::spec::dims::POWER
            || actual.watts.dims != old.watts.dims
            || old.until.dims != fs_project::spec::dims::TIME
            || actual.until.dims != old.until.dims
            || !old.until.value.is_finite() || old.until.value <= start
            || actual.until.value.to_bits() != old.until.value.to_bits()
        {
            return Err(invalid("sampled schedule is not the original delivered-watt history times its amplitude"));
        }
        start = old.until.value;
    }
    if start.to_bits() != coarse[coarse.len() - 1].to_bits() {
        return Err(invalid("schedule amplitude does not cover the original complete horizon"));
    }
    let intervals = array(report, "intervals", 10_000)?;
    if intervals.len() != count {
        return Err(invalid("schedule sensitivity requires every original fine-grid interval"));
    }
    let mut schedule_index = 0;
    let mut row_index = 0;
    let mut coarse_start = 0.0;
    let mut start = 0.0;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for &coarse_end in coarse {
        // Exact arithmetic of the native scheduled-workload TimeGrids owner.
        // This validates its retained grid; it does not choose a new grid.
        let middle = coarse_start + 0.5 * (coarse_end - coarse_start);
        if !(coarse_end.is_finite() && middle > coarse_start && middle < coarse_end) {
            return Err(invalid("schedule sensitivity requires representable positive half steps"));
        }
        for end in [middle, coarse_end] {
            let row = &intervals[row_index];
            if row.f64_field("step") != Some((row_index + 1) as f64)
                || number(row, "start_s")?.to_bits() != start.to_bits()
                || number(row, "end_s")?.to_bits() != end.to_bits()
            {
                return Err(invalid("schedule sensitivity interval differs from the original accepted fine grid"));
            }
            // Schedule events take effect on the FOLLOWING interval. The
            // ending interval retains its left-limit delivered power.
            while schedule_index + 1 < original.steps.len()
                && start >= original.steps[schedule_index].until.value
            {
                schedule_index += 1;
            }
            let segment = &original.steps[schedule_index];
            if end > segment.until.value {
                return Err(invalid("schedule sensitivity interval crosses an authored workload switch"));
            }
            let derivative = number(named(row, "regional_power", "region", &original.region)?,
                "dtemperature_ddelivered_watts_k_per_w")?;
            if row_index >= prefix_steps && derivative != 0.0 {
                return Err(invalid("future scheduled heating cannot influence the selected earlier peak"));
            }
            let term = derivative * segment.watts.value;
            let increment = term - correction;
            let next = sum + increment;
            let next_correction = (next - sum) - increment;
            if !term.is_finite() || !increment.is_finite() || !next.is_finite()
                || !next_correction.is_finite()
                || (derivative != 0.0 && segment.watts.value > 0.0 && term == 0.0)
            {
                return Err(invalid("schedule sensitivity contraction is not representable"));
            }
            sum = next;
            correction = next_correction;
            start = end;
            row_index += 1;
        }
        coarse_start = coarse_end;
    }
    Ok(sum)
}

#[cfg(test)]
mod tests;
