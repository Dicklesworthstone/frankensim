//! Bounded modal output FIRs using the original frequency and time-image owners.
use super::*;
use fs_duct::{DuctError, modal::outlet::{self, MouthResponse}};
use fs_math::c64::C64;

pub(super) fn load_at<'a>(loads: MmBankLoad<'a>, combo: usize) -> MmLoad<'a> {
    match loads { MmBankLoad::Shared(load) => load, MmBankLoad::PerCombo(loads) => loads[combo] }
}

pub(super) fn admit(ducts: &[Duct], labels: &[&str], gas: &GasState,
    loads: MmBankLoad<'_>, line: &MmLineConfig, output: MouthFlowConfig)
    -> Result<(), RenderError>
{
    if ducts.is_empty() || ducts.len() != labels.len() || line.sample_rate_hz == 0
        || !(1..=fs_duct::modal::MAX_MODES).contains(&line.n_modes)
        || !(512..=65_536).contains(&output.max_taps) || !output.max_taps.is_power_of_two()
        || !output.tail_tolerance.is_finite() || !(0.0..1.0).contains(&output.tail_tolerance)
        || output.tail_tolerance == 0.0
        || !output.tap_change_tolerance.is_finite()
        || !(0.0..1.0).contains(&output.tap_change_tolerance) || output.tap_change_tolerance == 0.0
    {
        return Err(invalid("modal mouth realization requires valid paths, modes, clock and finite tolerances"));
    }
    let total = ducts.len().checked_mul(line.n_modes)
        .and_then(|n| n.checked_mul(output.max_taps));
    if total.is_none_or(|n| n > output.max_total_taps) {
        return Err(invalid("modal mouth realization exceeds its total ceiling-tap budget"));
    }
    if let MmBankLoad::PerCombo(list) = loads {
        if list.len() != ducts.len() {
            return Err(invalid("modal mouth output requires one load per original valve path"));
        }
    }
    let dt = 1.0 / f64::from(line.sample_rate_hz);
    for (combo, duct) in ducts.iter().enumerate() {
        let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
        let transit = length / gas.sound_speed / dt;
        if !length.is_finite() || length <= 0.0 || !transit.is_finite()
            || transit < 1.0 || 4.0 * transit > output.max_taps as f64
        {
            return Err(invalid("modal mouth transit does not fit the declared output ladder"));
        }
        match load_at(loads, combo) {
            MmLoad::Analytic(_) => {},
            MmLoad::Tabulated(table) => {
                let radius = duct.segments.last().ok_or_else(|| invalid("empty mouth path"))?
                    .outlet_radius();
                table.validate_mouth(length, radius, gas).map_err(physical)?;
                for bin in [1, output.max_taps / 2] {
                    table.z_at(core::f64::consts::TAU * bin as f64
                        / (output.max_taps as f64 * dt)).map_err(physical)?;
                }
            }
            MmLoad::FlangedModal(exterior) => {
                // The input bank's own preflight (mm_line::validate_mouth_load): one static
                // pressure across the interface, and the exterior ka ceiling must cover the
                // whole Nyquist band, since every later DFT rung evaluates this same matrix.
                if gas.pressure != exterior.pressure {
                    return Err(invalid("flanged modal mouth needs the bore and exterior at one static pressure"));
                }
                let radius = duct.segments.last().ok_or_else(|| invalid("empty mouth path"))?
                    .outlet_radius();
                fs_duct::flanged::flanged_modal_impedance(exterior, radius,
                    core::f64::consts::PI * f64::from(line.sample_rate_hz), 1).map_err(physical)?;
            }
            _ => return Err(invalid("modal mouth output needs an analytic, flanged modal or strict tabulated load")),
        }
    }
    Ok(())
}

fn response(duct: &Duct, gas: &GasState, load: MmLoad<'_>, line: &MmLineConfig,
    omega: f64, dc: bool) -> Result<MouthResponse, RenderError>
{
    let evaluate = |loss| match load {
        MmLoad::Analytic(termination) => outlet::mouth_response(
            duct, gas, omega, loss, termination, line.n_modes, line.extra_slices),
        MmLoad::Tabulated(table) if dc => {
            let mouth = table.mouth().ok_or(DuctError::BadParameter {
                what: "strict modal output needs an explicit static mouth load",
            })?;
            outlet::mouth_response_load(duct, gas, omega, loss,
                C64::from_re(mouth.dc_impedance), line.n_modes, line.extra_slices)
        }
        MmLoad::Tabulated(table) => outlet::mouth_response_tabulated(
            duct, gas, omega, loss, table, line.n_modes, line.extra_slices),
        // Same coupled matrix as the input bank (mm_line::resolve_load), including its quasi-static
        // probe: no scalar fallback at DC, so `dc` needs no separate arm.
        MmLoad::FlangedModal(exterior) => outlet::mouth_response_flanged(
            duct, gas, exterior, omega, loss, line.n_modes, line.extra_slices),
        _ => Err(DuctError::BadParameter { what: "unsupported modal output load" }),
    };
    // This is the original input bank's loss admission, not an analytic
    // radiation substitution. Strict table coverage never falls back.
    match evaluate(line.loss) {
        Err(DuctError::TooNarrow { .. }) if matches!(line.loss, LossModel::WideTube) =>
            evaluate(LossModel::Bessel).map_err(physical),
        result => result.map_err(physical),
    }
}

// Scale before squaring: weak modes must not vanish merely because h*h
// underflows. Every nonzero mode meets its OWN relative convergence gates.
fn metrics(ir: &[f64], previous: &[f64]) -> (f64, f64) {
    let scale = ir.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return (0.0, if previous.is_empty() || previous.iter().any(|x| *x != 0.0) {
            f64::INFINITY
        } else { 0.0 });
    }
    let total: f64 = ir.iter().map(|h| (h / scale).powi(2)).sum();
    let tail: f64 = ir[3 * (ir.len() / 4)..].iter().map(|h| (h / scale).powi(2)).sum();
    let change = if previous.is_empty() { f64::INFINITY } else {
        let difference: f64 = ir.iter().zip(previous)
            .map(|(a, b)| ((a - b) / scale).powi(2)).sum();
        (difference / total).sqrt()
    };
    (tail / total, change)
}

/// `observed` is how many LEADING modes get a FIR and a convergence gate (mode zero is total flow). The whole
/// `line.n_modes` image is still solved at every bin, so each observed flow is the coupled one; only the unobserved
/// modes' time images are skipped. A sequence observes every mode. A caller whose receiver cannot see the higher
/// modes (an on-axis point in the far field: every mode n >= 1 integrates to zero over the aperture) observes one.
#[allow(clippy::too_many_arguments)]
pub(super) fn path(duct: &Duct, label: &str, gas: &GasState, load: MmLoad<'_>,
    line: &MmLineConfig, reference: f64, output: MouthFlowConfig, observed: usize)
    -> Result<(Vec<DelayedFilter>, MouthFlowReceipt), RenderError>
{
    if !(1..=line.n_modes).contains(&observed) {
        return Err(invalid("modal mouth output observes between one and all retained modes"));
    }
    let dt = 1.0 / f64::from(line.sample_rate_hz);
    let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
    let transit = length / gas.sound_speed / dt;
    let mut n_fft = ((4.0 * transit).ceil() as usize).max(256).next_power_of_two()
        .min(output.max_taps / 2);
    // Same quasi-static probe as the input bank. A strict table's explicit
    // static impedance is independent of its positive-frequency support.
    let dc = response(duct, gas, load, line, 1e-3, true)?;
    let radius = dc.mouth_radius_m;
    let static_flow = dc.per_incident_plane_pressure(reference).map_err(physical)?;
    let mut spectra: Vec<Vec<C64>> = Vec::new();
    let mut previous = vec![Vec::<f64>::new(); observed];
    let mut rungs = Vec::new();
    loop {
        let half = n_fft / 2;
        let mut next = vec![vec![C64::ZERO; half + 1]; observed];
        if !spectra.is_empty() {
            for (new, old) in next.iter_mut().zip(&spectra) {
                for (k, &value) in old.iter().enumerate() { new[2*k] = value; }
            }
        }
        for k in 1..=half {
            if !spectra.is_empty() && k % 2 == 0 { continue; }
            let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
            let flow = response(duct, gas, load, line, omega, false)?
                .per_incident_plane_pressure(reference).map_err(physical)?;
            for (mode, &value) in next.iter_mut().zip(&flow) { mode[k] = value; }
        }
        for (mode, dc) in next.iter_mut().zip(&static_flow) { mode[0] = C64::from_re(dc.re); }
        let taps: Vec<Vec<f64>> = next.iter().map(|spectrum|
            crate::driving_point::realize_reflectance_ir_from_positive_spectrum(
                spectrum, n_fft, transit, true)).collect();
        if taps.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("modal mouth time image left the finite set"));
        }
        let (tail_energy, tap_change): (Vec<_>, Vec<_>) = taps.iter().zip(&previous)
            .map(|(ir, old)| metrics(ir, old)).unzip();
        let converged = n_fft as f64 > transit + 2.0
            && tail_energy.iter().all(|&e| e < output.tail_tolerance)
            && tap_change.iter().all(|&e| e < output.tap_change_tolerance);
        rungs.push(MouthFlowRung { n_fft, tail_energy, tap_change });
        if converged {
            let filters = taps.into_iter().map(|ir|
                DelayedFilter::from_impulse_response(dt, ir)
                    .map_err(|e| physical(format!("modal mouth FIR refused: {e:?}"))))
                .collect::<Result<_, _>>()?;
            return Ok((filters, MouthFlowReceipt { label: label.to_string(),
                mouth_radius_m: radius, n_fft, rungs }));
        }
        if n_fft == output.max_taps {
            return Err(physical(format!(
                "modal mouth output for {label:?} did not converge at {} taps: {:?}",
                output.max_taps, rungs.last())));
        }
        previous = taps.iter().map(|ir| ir[..half].to_vec()).collect();
        spectra = next;
        n_fft *= 2;
    }
}

#[cfg(test)]
mod tests {
    use super::metrics;
    #[test]
    fn convergence_is_scale_invariant_and_does_not_silence_weak_modes() {
        let ir = [0.0, 1.0, 0.5, 0.0, 0.0, 0.0, 1e-6, 0.0];
        let old = [0.0, 1.0, 0.49, 0.0];
        let reference = metrics(&ir, &old);
        for scale in [1e-200, 1.0, 1e200] {
            let scaled = ir.map(|v| v * scale);
            let previous = old.map(|v| v * scale);
            let actual = metrics(&scaled, &previous);
            assert!((actual.0-reference.0).abs() < 1e-25);
            assert!((actual.1-reference.1).abs() < 1e-15);
        }
        assert_eq!(metrics(&[0.0; 8], &[0.0; 4]), (0.0, 0.0));
        assert!(metrics(&[0.0; 8], &[]).1.is_infinite());
        assert!(metrics(&[0.0; 8], &[1e-200; 4]).1.is_infinite());
    }
}
