//! The plane-mode equivalent of the complete flanged mouth.
//!
//! A plane-wave interior (the nonlinear Euler bore) carries one mode, so it can
//! take a mouth load only as one impedance. The multimodal flanged mouth of
//! [`mm_input_impedance_flanged`] is richer: a plane wave arriving at the mouth
//! also excites the higher modes there, and the bore's own mode conversion
//! (radius steps, flare) feeds them back. [`equivalent_plane_mouth_impedance`]
//! folds all of that into the one load the plane mode sees at the mouth
//!
//! `Z_eq = p_0 / U_0`,  `p = Z_rad U_mouth`,  `U_mouth = flow_map[:, 0]`
//!
//! for a plane inlet flow (higher inlet flows zero, the same `Z00` source
//! closure as the line bank), with `Z_rad` the complete flanged matrix at the
//! ACTUAL last staircase station radius and `flow_map` the original mouth map
//! of [`outlet::MouthResponse`]. Convention: `e^{-i omega t}`, like every
//! modal and scalar impedance here, so a flanged piston has negative reactance.
//!
//! [`plane_mouth_kernel`] realizes `Z_eq` as the causal real weights a
//! time-domain end uses (`p = Z_inf (U_t - sum_j w_j U_{t-j})`, the form of the
//! nonlinear bore's piston radiation). `Z_eq` is not a certified-passive
//! impedance: the tests report its resistance, and the kernel is built from
//! that resistance alone, so it is causal by construction.
use super::*;

/// The largest `ka` the equivalent-load kernel is evaluated to: the band of its recorded measurements (24), not the flanged load's admitted
/// cap (raised to 28, MB 13958), so the kernels, their claims and the bits of every bore that uses this builder stay as measured.
const EVALUATED_MAX_KA: f64 = 24.0;

/// Plane-mode equivalent of the complete flanged mouth at one frequency.
///
/// `state` is the bore gas, `exterior` the room; they must share a static
/// pressure. `n_modes` and `extra_slices` are the multimodal solve's own.
///
/// # Errors
/// As [`mm_input_impedance_flanged`], or a singular mouth flow (a plane inlet
/// flow that produces no plane mouth flow has no equivalent load).
pub fn equivalent_plane_mouth_impedance(
    duct: &Duct,
    state: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    n_modes: usize,
    extra_slices: usize,
) -> Result<C64, DuctError> {
    if n_modes == 0 {
        return Err(DuctError::BadParameter {
            what: "at least one mode (the plane wave)",
        });
    }
    if n_modes > MAX_MODES {
        return Err(DuctError::BadParameter {
            what: "mode count exceeds MAX_MODES",
        });
    }
    let geometry = PreparedModalDuct::new(duct, n_modes, extra_slices)?;
    equivalent_prepared(&geometry, state, exterior, omega, loss)
}

fn equivalent_prepared(
    geometry: &PreparedModalDuct<'_>,
    state: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
) -> Result<C64, DuctError> {
    if !(omega.is_finite() && omega > 0.0) {
        return Err(DuctError::BadParameter {
            what: "omega must be finite and positive",
        });
    }
    let response = mm_evaluate_prepared(
        geometry,
        state,
        omega,
        loss,
        MmLoad::Flanged(exterior),
        true,
    )?;
    let n = response.input.n_modes;
    let radiation =
        crate::flanged::flanged_modal_impedance(exterior, response.mouth_radius_m, omega, n)?;
    // Column 0 of the row-major map: mouth modal flows per unit plane inlet flow.
    let flow: Vec<C64> = (0..n).map(|row| response.flow_map[row * n]).collect();
    if !(flow[0].abs() > 0.0 && flow.iter().all(|u| u.re.is_finite() && u.im.is_finite())) {
        return Err(DuctError::Singular);
    }
    let pressure = (0..n).fold(C64::ZERO, |sum, j| sum + radiation[j] * flow[j]);
    let z = pressure * flow[0].recip();
    if !(z.re.is_finite() && z.im.is_finite()) {
        return Err(DuctError::Singular);
    }
    Ok(z)
}

/// How [`plane_mouth_kernel`] samples and realizes the equivalent load.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MouthKernelSpec {
    /// Time step of the realized weights [s].
    pub step_s: f64,
    /// Requested frequency spacing of the evaluation grid [Hz]; the kernel is
    /// `1 / (bin_hz step_s)` taps long before its tail is dropped.
    pub bin_hz: f64,
    /// Top of the band over which the realized load is claimed to match `Z_eq` [Hz]. The load is evaluated far beyond it, up to the
    /// flanged radiation limit, so that the resistance omitted above the evaluation does not leak into the reactance inside the claim;
    /// the step must be fine enough for its Nyquist frequency to allow that (the evaluation must reach 0.9 of the limit), and the claim must lie below a third of it.
    pub band_hz: f64,
    /// Largest fraction of the kernel's energy dropped from its tail, in `(0, 1)`.
    pub tail_tolerance: f64,
}

/// Causal real weights of a plane-mode mouth load, in the form of the
/// nonlinear bore's piston radiation: `p = z_inf (U_t - sum_j weights[j] U_{t-j})`.
#[derive(Debug, Clone, PartialEq)]
pub struct MouthKernel {
    /// Time step of the weights [s].
    pub step_s: f64,
    /// `rho c / (pi a^2)` of the BORE gas (`state`) at the actual mouth radius [Pa s/m^3]: the characteristic impedance the time-domain bore end is matched to
    /// (`NonlinearBore::with_mouth_kernel` checks it against the gas the bore holds). With one gas for bore and exterior it is the exterior's, as it always was, bit
    /// for bit; with two it differs from the exterior's `rho c / (pi a^2)` by their `rho c` ratio.
    pub z_inf: f64,
    /// The actual last staircase station radius the load was built for [m].
    pub mouth_radius_m: f64,
    /// `w_0 ..`; they sum to one (zero load at zero frequency) before the tail is dropped.
    pub weights: Vec<f64>,
    /// Fraction of the kernel's energy in the dropped tail (at most the declared tolerance).
    pub dropped_tail_energy: f64,
    /// Sum of the dropped weights: the error of the kept weights' sum, hence of the load at zero frequency.
    pub dropped_sum: f64,
    /// Frequency spacing of the evaluation grid actually used [Hz].
    pub bin_hz: f64,
    /// Top of the evaluated band [Hz]: the load is `Z_inf` above it.
    pub evaluated_hz: f64,
}

impl MouthKernel {
    /// The realized load at `omega`, in the `e^{-i omega t}` convention of
    /// [`equivalent_plane_mouth_impedance`]: `z_inf (1 - conj(sum_j w_j e^{-i omega dt j}))`.
    #[must_use]
    pub fn impedance_at(&self, omega: f64) -> C64 {
        let mut transfer = C64::ZERO;
        for (j, &w) in self.weights.iter().enumerate() {
            let phase = -omega * self.step_s * j as f64;
            transfer = transfer + C64::new(w * det::cos(phase), w * det::sin(phase));
        }
        (C64::ONE - transfer.conj()).scale(self.z_inf)
    }
}

/// Realize the plane-mode equivalent mouth load as causal weights.
///
/// The target resistance `Re Z_eq` is sampled at `omega_k = 2 pi k / (N step_s)`
/// for `k = 1 ..=` the evaluated top (the flanged radiation limit, or just below
/// the step's Nyquist frequency), with `Z_eq(0) = 0` (radiation vanishes at
/// zero frequency), `N = 1 / (bin_hz step_s)` rounded to an even count. The
/// weights are the causal doubling of the inverse cosine transform of
/// `1 - Re Z_eq / Z_inf` (the even part of the impulse response determines the
/// causal one), so causality holds exactly; the reactance follows by
/// Kramers-Kronig and the tests compare it with `Z_eq`. The resistance omitted
/// above the evaluated top would otherwise leak into the reactance as an error
/// growing with frequency; evaluating to the radiation limit, where the resistance
/// is within a few percent of `Z_inf` and oscillating, keeps that error far below
/// the claim (offline on a bell bore: 5e-5 of `Z_inf` at a claim a third of the top).
///
/// `state` is the bore gas and `exterior` the room (one static pressure). `Re Z_eq` is the two-gas flanged law of [`equivalent_plane_mouth_impedance`]; its reference
/// impedance `Z_inf` is the BORE gas's, because the time-domain bore end is matched to the gas it holds. Above the evaluated top the realized load is therefore the
/// bore gas's `Z_inf` where the flanged radiation resistance tends to the exterior's: a mismatch of `(rho c)_exterior / (rho c)_bore - 1` in a band beyond the claim
/// (1.3 % for the 27 C saturated breath state against dry 20 C air). With equal gases nothing differs from the one-gas builder, bit for bit.
///
/// # Errors
/// A bad spec, any frequency's refusal from the multimodal solve, or a grid
/// whose band reaches the Nyquist frequency.
pub fn plane_mouth_kernel(
    duct: &Duct,
    state: &GasState,
    exterior: &GasState,
    loss: LossModel,
    n_modes: usize,
    extra_slices: usize,
    spec: MouthKernelSpec,
) -> Result<MouthKernel, DuctError> {
    let bad = |what: &'static str| DuctError::BadParameter { what };
    if !(spec.step_s.is_finite() && spec.step_s > 0.0) {
        return Err(bad("mouth kernel needs a positive finite step"));
    }
    if !(spec.bin_hz.is_finite()
        && spec.bin_hz > 0.0
        && spec.band_hz.is_finite()
        && spec.band_hz > 2.0 * spec.bin_hz)
    {
        return Err(bad(
            "mouth kernel needs a positive bin spacing and a band of at least two bins",
        ));
    }
    if !(spec.tail_tolerance > 0.0 && spec.tail_tolerance < 1.0) {
        return Err(bad("mouth kernel tail tolerance must lie in (0, 1)"));
    }
    if n_modes == 0 || n_modes > MAX_MODES {
        return Err(bad("mouth kernel mode count outside 1..=MAX_MODES"));
    }
    let mut n = (1.0 / (spec.bin_hz * spec.step_s)).round() as usize;
    n += n % 2;
    let bin_hz = 1.0 / (n as f64 * spec.step_s);
    let geometry = PreparedModalDuct::new(duct, n_modes, extra_slices)?;
    let mouth_radius_m = geometry.mouth_radius_m();
    // The bore end is matched to the gas the bore holds: the bore gas's characteristic impedance, which is the exterior's when the two gases are one.
    let z_inf = state.density * state.sound_speed
        / (core::f64::consts::PI * mouth_radius_m * mouth_radius_m);
    if !(z_inf.is_finite() && z_inf > 0.0) {
        return Err(DuctError::Singular);
    }
    // Evaluate to the flanged radiation limit or just below Nyquist, whichever is lower.
    let radiation_limit_hz = EVALUATED_MAX_KA * exterior.sound_speed
        / (core::f64::consts::TAU * mouth_radius_m);
    let evaluated_hz = (0.999 * radiation_limit_hz).min(0.98 * 0.5 / spec.step_s);
    let top = (evaluated_hz / bin_hz).floor() as usize;
    if n < 8 || top < 8 || top + 2 >= n / 2 {
        return Err(bad(
            "mouth kernel step is too coarse for the evaluated band",
        ));
    }
    // The resistance omitted above the evaluation leaks into the reactance, and its size falls only as the radiation
    // limit is approached: a step whose Nyquist frequency stops the evaluation short of 0.9 of it cannot meet the claim.
    if evaluated_hz < 0.9 * radiation_limit_hz {
        return Err(bad(
            "mouth kernel step is too coarse to evaluate the load up to the flanged radiation limit",
        ));
    }
    if spec.band_hz * 3.0 > evaluated_hz {
        return Err(bad(
            "mouth kernel band must lie below a third of the evaluated band (the flanged radiation limit and the step's Nyquist frequency)",
        ));
    }
    // Re H_k = 1 - Re Z_eq / Z_inf, tapered over the last tenth of the band.
    let taper_from = (0.9 * top as f64).floor();
    let mut spectrum = vec![0.0; top + 1];
    spectrum[0] = 1.0;
    for (k, slot) in spectrum.iter_mut().enumerate().skip(1) {
        let omega = core::f64::consts::TAU * k as f64 * bin_hz;
        let z = equivalent_prepared(&geometry, state, exterior, omega, loss)?;
        let taper = if (k as f64) <= taper_from {
            1.0
        } else {
            0.5 * (1.0
                + det::cos(
                    core::f64::consts::PI * (k as f64 - taper_from) / (top as f64 - taper_from),
                ))
        };
        *slot = (1.0 - z.re / z_inf) * taper;
    }
    // Inverse cosine transform, then double the positive times (causal).
    let half = n / 2;
    let mut weights = Vec::with_capacity(half);
    for j in 0..half {
        let theta = core::f64::consts::TAU * j as f64 / n as f64;
        let c1 = det::cos(theta);
        // cos(k theta) for k = 1, 2, ... by the three-term recurrence
        let (mut previous, mut current) = (1.0, c1);
        let mut sum = spectrum[0];
        for &value in spectrum.iter().skip(1) {
            sum += 2.0 * value * current;
            let next = 2.0 * c1 * current - previous;
            previous = current;
            current = next;
        }
        let s = sum / n as f64;
        weights.push(if j == 0 { s } else { 2.0 * s });
    }
    // Drop the tail whose energy is within the declared fraction.
    let total: f64 = weights.iter().map(|w| w * w).sum();
    let mut kept = weights.len();
    let (mut dropped, mut dropped_sum) = (0.0, 0.0);
    while kept > 1 && dropped + weights[kept - 1] * weights[kept - 1] <= spec.tail_tolerance * total
    {
        dropped += weights[kept - 1] * weights[kept - 1];
        dropped_sum += weights[kept - 1];
        kept -= 1;
    }
    weights.truncate(kept);
    Ok(MouthKernel {
        step_s: spec.step_s,
        z_inf,
        mouth_radius_m,
        weights,
        dropped_tail_energy: if total > 0.0 { dropped / total } else { 0.0 },
        dropped_sum,
        bin_hz,
        evaluated_hz: top as f64 * bin_hz,
    })
}

#[cfg(test)]
mod tests;
