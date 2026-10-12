//! Multi-partial octave tuning in the existing pinned stiff-string model.
//!
//! A supplied weight vector selects the lower:upper partial pairs 2:1, 4:2,
//! 6:3 and 8:4. Mixed policies minimize weighted squared beat mismatches divided
//! by the upper partial number. This is a declared local octave objective, not
//! a psychoacoustic fit, measured Railsback curve or global temperament search.
//! Only physical tension changes. The source's mass, EI and geometry stay fixed.
use super::{Course, PI, TAU, at_partial, det, tune_octaves};

/// Nonnegative relative weights for 2:1, 4:2, 6:3 and 8:4, respectively.
/// At least one must be positive. Scaling every weight leaves the policy alone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OctaveWeights([f64; 4]);
impl OctaveWeights {
    pub const TWO_ONE: Self = Self([1.0, 0.0, 0.0, 0.0]);

    pub fn new(weights: [f64; 4]) -> Result<Self, String> {
        if weights.iter().any(|w| !w.is_finite() || *w < 0.0) {
            return Err("octave weights must be finite and nonnegative".into());
        }
        let scale = weights.iter().copied().fold(0.0, f64::max);
        if scale == 0.0 { return Err("at least one octave weight must be positive".into()); }
        let normalized = weights.map(|w| w / scale);
        if weights.iter().zip(normalized).any(|(&w, v)| w > 0.0 && v == 0.0) {
            return Err("relative octave weights are not representable".into());
        }
        Ok(Self(normalized))
    }

    pub fn parse(text: &str) -> Result<Self, String> {
        let mut fields = text.split(',');
        let mut weights = [0.0; 4];
        for w in &mut weights {
            *w = fields.next().ok_or("octave weights need four comma-separated values")?
                .trim().parse().map_err(|_| "invalid octave weight")?;
        }
        if fields.next().is_some() { return Err("octave weights need exactly four values".into()); }
        Self::new(weights)
    }

    pub fn values(self) -> [f64; 4] { self.0 }
}

#[derive(Clone, Copy, Debug)]
pub struct Profile {
    /// Policy below the fixed C4..B4 temperament octave.
    pub bass: OctaveWeights,
    /// Policy above the fixed C4..B4 temperament octave.
    pub treble: OctaveWeights,
    /// Optional ceiling on the already-tuned neighbour's target partial [Hz].
    /// Excluded pairs are reported, not replaced by an invented lower partial.
    pub max_partial_hz: Option<f64>,
}
impl Default for Profile {
    fn default() -> Self {
        Self { bass: OctaveWeights::TWO_ONE, treble: OctaveWeights::TWO_ONE,
            max_partial_hz: None }
    }
}

#[derive(Clone, Debug)]
pub struct OctaveReport {
    pub lower_midi: u8,
    /// Lower partial minus upper partial, in physical Hz, for every active pair.
    /// None means a zero weight or a target above the explicitly supplied band.
    pub signed_beats_hz: [Option<f64>; 4],
    /// sqrt(sum(w_n * (beat_n/n)^2) / sum(w_n)), in upper-fundamental Hz.
    /// This is the declared octave objective, not loudness-weighted roughness.
    pub weighted_rms_hz: f64,
}

/// Use the existing anchor/traversal/rollback owner, replacing only its local
/// octave equation. The default profile preserves the original 2:1 path exactly.
/// The complete scale is unchanged on ANY admission, solve or report failure.
pub fn tune(courses: &mut [Course], concert_a_hz: f64, profile: Profile)
    -> Result<Vec<OctaveReport>, String> {
    if profile.max_partial_hz.is_some_and(|hz| !hz.is_finite() || hz <= 0.0) {
        return Err("tuning band must be finite and positive".into());
    }
    let mut report = Vec::with_capacity(courses.len().min(88));
    tune_octaves(courses, concert_a_hz, |course, neighbour, upper| {
        let weights = if upper { profile.treble } else { profile.bass };
        let (tuned, row) = match_octave(course, neighbour, upper, weights, profile.max_partial_hz)?;
        report.push(row);
        Ok(tuned)
    })?;
    report.sort_by_key(|r| r.lower_midi);
    Ok(report)
}

#[derive(Clone, Copy, Default)]
struct Pair { partial: u8, target_hz: f64, weight: f64 }

fn match_octave(course: Course, neighbour: Course, upper: bool,
    weights: OctaveWeights, band: Option<f64>) -> Result<(Course, OctaveReport), String> {
    let mut pairs = [Pair::default(); 4];
    let mut count = 0;
    for (i, &weight) in weights.0.iter().enumerate() {
        if weight == 0.0 { continue; }
        let n = i + 1;
        let target = neighbour.partial_hz(if upper { 2*n } else { n }, neighbour.tension_n);
        if !target.is_finite() || target <= 0.0 {
            return Err(format!("key {} has an unrepresentable octave target", course.midi));
        }
        if band.is_some_and(|maximum| target > maximum) { continue; }
        pairs[i] = Pair { partial: (if upper { n } else { 2*n }) as u8,
            target_hz: target, weight };
        count += 1;
    }
    if count == 0 {
        return Err(format!("key {} has no weighted octave partial in the supplied tuning band", course.midi));
    }
    // The band may exclude the largest original weight. Renormalize admitted
    // weights before arithmetic so a retained tiny-weight group cannot underflow
    // its objective or falsely satisfy the stationarity test with zero work.
    let active_scale = pairs.iter().map(|p| p.weight).fold(0.0_f64, f64::max);
    for pair in &mut pairs { pair.weight /= active_scale; }
    // Exact single-pair solve, including the bit-preserving original 2:1 image.
    let tuned = if count == 1 {
        let pair = pairs.iter().find(|p| p.weight > 0.0).expect("one admitted pair");
        at_partial(course, pair.partial, pair.target_hz)?
    } else {
        minimize(course, &pairs)?
    };
    let mut beats = [None; 4];
    let mut squared = 0.0;
    let mut total_weight = 0.0;
    for (i, pair) in pairs.iter().enumerate().filter(|(_, p)| p.weight > 0.0) {
        let actual = tuned.partial_hz(usize::from(pair.partial), tuned.tension_n);
        let beat = if upper { pair.target_hz - actual } else { actual - pair.target_hz };
        beats[i] = Some(beat);
        squared += pair.weight * (beat / (i+1) as f64).powi(2);
        total_weight += pair.weight;
    }
    let weighted_rms_hz = det::sqrt(squared / total_weight);
    if !weighted_rms_hz.is_finite() || beats.iter().flatten().any(|v| !v.is_finite()) {
        return Err(format!("key {} octave beat report is not representable", course.midi));
    }
    Ok((tuned, OctaveReport { lower_midi: if upper { neighbour.midi } else { course.midi },
        signed_beats_hz: beats, weighted_rms_hz }))
}

/// For variable string partial n, f_n/n = sqrt(A*(T+b_n)), A=1/(4*mu*L^2),
/// b_n=EI*(n*pi/L)^2. With s_n=(target_n/n)^2/A the objective derivative is
/// proportional to g(T)=sum w_n*(1-sqrt(s_n/(T+b_n))). Each derivative term is
/// strictly increasing: g'(T)=sum w_n*sqrt(s_n)/(2*(T+b_n)^(3/2)) > 0.
/// Hence its positive root is the unique constrained interior minimum, not a
/// local optimizer guess. Dividing lower-string mismatches by 2*n rather than
/// n changes the complete objective by 1/4, never its minimizer.
fn minimize(mut course: Course, pairs: &[Pair; 4]) -> Result<Course, String> {
    let mut bend = [0.0_f64; 4];
    let mut target = [0.0_f64; 4];
    let mut scale = 0.0_f64;
    for (i, p) in pairs.iter().enumerate().filter(|(_, p)| p.weight > 0.0) {
        let k = f64::from(p.partial) * PI / course.length_m;
        bend[i] = course.flexural_rigidity_nm2 * k * k;
        target[i] = (TAU * p.target_hz / k).powi(2) * course.linear_density_kg_m;
        if !bend[i].is_finite() || !target[i].is_finite() || target[i] <= 0.0 {
            return Err(format!("key {} multi-partial tension coefficients are not representable", course.midi));
        }
        scale = scale.max(bend[i]).max(target[i]);
    }
    for (i, p) in pairs.iter().enumerate().filter(|(_, p)| p.weight > 0.0) {
        let b = bend[i] / scale;
        let s = target[i] / scale;
        if (bend[i] > 0.0 && b == 0.0) || s == 0.0 {
            return Err(format!("key {} scaled octave coefficients underflow", course.midi));
        }
        bend[i] = b; target[i] = s;
        // p is used by the filter: keep all explicitly positive weights.
        debug_assert!(p.weight > 0.0);
    }
    let gradient = |t: f64| -> f64 {
        pairs.iter().enumerate().filter(|(_, p)| p.weight > 0.0).map(|(i, p)| {
            // Separate roots avoid overflow of the positive ratio near T=0.
            p.weight * (1.0 - det::sqrt(target[i]) / det::sqrt(t + bend[i]))
        }).sum()
    };
    let mut lo = 0.0;
    let mut hi = pairs.iter().enumerate().filter(|(_, p)| p.weight > 0.0)
        .map(|(i, _)| target[i] - bend[i]).fold(0.0_f64, f64::max);
    // A bending-only optimum is not an admitted tensioned string. Do not clamp.
    let at_zero = gradient(0.0);
    if at_zero.is_nan() || at_zero >= 0.0 || hi <= 0.0 {
        return Err(format!("key {} multi-partial optimum requires nonpositive tension", course.midi));
    }
    let weight: f64 = pairs.iter().map(|p| p.weight).sum();
    let tolerance = 128.0 * f64::EPSILON * weight;
    let mut root = hi;
    let mut converged = false;
    for _ in 0..96 {
        let g = gradient(root);
        if !g.is_finite() { return Err("nonfinite multi-partial tuning gradient".into()); }
        if g.abs() <= tolerance { converged = true; break; }
        if g > 0.0 { hi = root; } else { lo = root; }
        let next = lo + 0.5 * (hi - lo);
        if next == lo || next == hi { break; }
        root = next;
    }
    if !converged { return Err(format!("key {} multi-partial tension solve did not converge", course.midi)); }
    course.tension_n = root * scale;
    if !course.tension_n.is_finite() || course.tension_n <= 0.0 {
        return Err(format!("key {} multi-partial tension is not positive and representable", course.midi));
    }
    course.validate().map_err(|e| format!("multi-partial tuning key {}: {e}", course.midi))?;
    // Verify through the played scale's frequency owner, not just scaled algebra.
    let residual: f64 = pairs.iter().filter(|p| p.weight > 0.0).map(|p| {
        let actual = course.partial_hz(usize::from(p.partial), course.tension_n);
        p.weight * (1.0 - p.target_hz / actual)
    }).sum();
    if !residual.is_finite() || residual.abs() > 1e-11 * weight {
        return Err(format!("key {} failed its physical multi-partial stationarity residual", course.midi));
    }
    Ok(course)
}

#[cfg(test)]
#[path = "partial_tuning_tests.rs"]
mod tests;
