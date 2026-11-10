//! Volume flow at the terminating plane of a duct, from the wave incident on
//! its inlet.
//!
//! A characteristic line returns the wave a duct reflects to its inlet. The
//! same linear duct also carries flow to its mouth. With the inlet split into
//! an incident wave `p+` and a returned wave `p-` referenced to the port
//! impedance `zc` (`p = p+ + p-`, `U = (p+ - p-) / zc`), the reflectance
//! `R = (Z - zc) / (Z + zc)` gives `U_in = 2 p+ / (Z + zc)`, so
//!
//! `T(w) = U_mouth / p+ = (U_mouth / U_in) * 2 / (Z + zc)`
//!
//! from the SAME transfer-matrix response the line was sampled from (same
//! viscothermal loss, wall and radiation load). No loss, radiation or
//! geometry law is added here; this only reads a second output of the
//! existing chain. The radiated pressure of that flow belongs to the baffled
//! receiver, not to this module.

use super::{
    DrivingPointError, LineLadder, duct_axial_length, realize_reflectance_ir_from_positive_spectrum,
    validate_ladder,
};
use fs_duct::{Duct, LossModel, MouthLoad, Termination, WallModel, input_impedance_wall};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

/// Angular frequencies [rad/s] probed for the DC bin, lowest first. The bin must
/// be the true limit `T(0+)`: a value taken from a higher frequency steps
/// against bin 1 and puts a constant `J/N` into every tap, a tail that falls
/// only as `1/N` and that no ladder can converge. (The first version took
/// 20 Hz, where the tube's inertance has already moved the real part by
/// ~0.4%, and missed the production tolerance at 65536 taps.) The reflectance
/// line probes 1e-3 rad/s the same way; a refused probe is answered from the
/// next decade up instead of a made-up value.
const DC_PROBE_RAD_S: [f64; 5] = [1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0];

fn dc_limit(
    response: &impl Fn(f64) -> Result<C64, DrivingPointError>,
    probes: &[f64],
) -> Result<f64, DrivingPointError> {
    let mut refusal = None;
    for &omega in probes {
        match response(omega) {
            Ok(value) => return Ok(value.re),
            Err(error) => refusal = Some(error),
        }
    }
    Err(refusal.unwrap_or(DrivingPointError::Invalid {
        what: "the DC probe list is empty",
    }))
}

/// One rung of the mouth-flow ladder.
#[derive(Debug, Clone, PartialEq)]
pub struct MouthFlowRung {
    /// DFT length of this rung.
    pub n_fft: usize,
    /// Relative energy in the last quarter of the taps.
    pub tail_energy: f64,
    /// `||h_n[..m] - h_{n/2}[..m]|| / ||h_n||` over the first `m = n / 4`
    /// taps, the part of the previous rung free of its own wrap-around
    /// (`+inf` on the first rung).
    pub tap_change: f64,
}

/// Receipt of a converged mouth-flow line.
#[derive(Debug, Clone, PartialEq)]
pub struct MouthFlowRealization {
    /// Converged tap count.
    pub n_fft: usize,
    /// Every rung traversed, in order.
    pub rungs: Vec<MouthFlowRung>,
}

/// `U_mouth / p+` at one angular frequency, `e^{-iwt}` convention.
///
/// # Errors
/// TMM refusal, including a lossless pole where `Z` is unrepresentable.
pub fn mouth_flow_per_incident_pressure(
    physics: &Duct,
    gas: &GasState,
    omega: f64,
    termination: Termination,
    wall: Option<&WallPin>,
    zc: f64,
) -> Result<C64, DrivingPointError> {
    let response = input_impedance_wall(physics, gas, omega, LossModel::Bessel, termination, wall)
        .map_err(DrivingPointError::Duct)?;
    let denominator = C64::new(response.impedance.re + zc, response.impedance.im);
    Ok((response.u_mouth_over_u_in * C64::new(2.0, 0.0)) / denominator)
}

/// [`mouth_flow_per_incident_pressure`] against an explicit mouth load `Z_L` [Pa s/m^3] (a referenced table's row): the
/// same chain the tabulated reflectance line samples ([`super::baked`]), so both lines read one load.
fn mouth_flow_against_load(
    physics: &Duct,
    gas: &GasState,
    omega: f64,
    load: C64,
    wall: Option<&WallPin>,
    zc: f64,
) -> Result<C64, DrivingPointError> {
    let response = super::baked::response(physics, gas, omega, load, WallModel::from_pin(wall))
        .map_err(DrivingPointError::Duct)?;
    let denominator = C64::new(response.impedance.re + zc, response.impedance.im);
    Ok((response.u_mouth_over_u_in * C64::new(2.0, 0.0)) / denominator)
}

/// Causal impulse response of [`mouth_flow_per_incident_pressure`], sized by a
/// convergence ladder: the DFT length doubles until the wrap-around tail and
/// the change of the shared taps between rungs are both below the ladder's
/// tolerance (`tail_tol` and its square root), or the ladder refuses.
///
/// # Errors
/// Empty geometry, non-positive rate or port impedance, TMM refusal, or
/// non-convergence at `ladder.max_taps`. A tabulated mouth load must carry its
/// mouth reference (plane, radius, gas) and cover the ladder's top rung; its DC
/// bin is the duct against the table's declared static load, never a table row.
pub fn realize_mouth_flow_ir<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<(Vec<f64>, MouthFlowRealization), DrivingPointError> {
    realize_with_dc_probes(
        physics, gas, gas, termination, sample_rate_hz, zc, wall, ladder, &DC_PROBE_RAD_S,
    )
}

// The DC probe list is a parameter only so a test can rebuild the stepped-DC
// variant beside the true-limit one and measure what the step costs. `exterior` is the gas a tabulated load must have
// been computed in; every bore quantity is `gas`'s.
#[allow(clippy::too_many_arguments)]
fn realize_with_dc_probes<'a>(
    physics: &Duct,
    gas: &GasState,
    exterior: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
    dc_probes: &[f64],
) -> Result<(Vec<f64>, MouthFlowRealization), DrivingPointError> {
    let length = duct_axial_length(physics);
    if !(length > 0.0 && length.is_finite()) {
        return Err(DrivingPointError::Invalid {
            what: "mouth-flow line needs positive axial length",
        });
    }
    if sample_rate_hz == 0 || !(zc.is_finite() && zc > 0.0) {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate and finite positive port impedance required",
        });
    }
    let max_taps = validate_ladder(ladder)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let transit = length / gas.sound_speed / dt;
    let demand = (4.0 * transit).ceil();
    if !demand.is_finite() {
        return Err(DrivingPointError::Invalid {
            what: "transit delay is not finite",
        });
    }
    let load = termination.into();
    // A declared termination takes exactly its old path; a table plays only through its
    // own mouth and must cover the top rung before any bin is sampled.
    let static_load = match load {
        MouthLoad::Analytic(_) => None,
        MouthLoad::Tabulated(table) => {
            let radius = physics
                .segments
                .last()
                .ok_or(DrivingPointError::Duct(fs_duct::DuctError::EmptyDuct))?
                .outlet_radius();
            let mouth = table.validate_mouth(length, radius, exterior).map_err(DrivingPointError::Duct)?;
            super::baked::validate_band(table, sample_rate_hz, max_taps)?;
            Some(C64::new(mouth.dc_impedance, 0.0))
        }
        MouthLoad::Matched => return Err(super::baked::matched_refusal()),
    };
    let response = |omega: f64| match load {
        MouthLoad::Tabulated(table) => {
            let z = table.z_at(omega).map_err(DrivingPointError::Duct)?;
            mouth_flow_against_load(physics, gas, omega, z, wall, zc)
        }
        MouthLoad::Analytic(termination) => {
            mouth_flow_per_incident_pressure(physics, gas, omega, termination, wall, zc)
        }
        MouthLoad::Matched => Err(super::baked::matched_refusal()),
    };
    let dc = match static_load {
        Some(z) => dc_limit(&|omega| mouth_flow_against_load(physics, gas, omega, z, wall, zc), dc_probes)?,
        None => dc_limit(&response, dc_probes)?,
    };
    transfer_ladder(&response, dc, transit, dt, max_taps, ladder)
}

/// [`dc_limit`] on the default probe list, for another transfer than a single duct's (a network port's flow).
pub(super) fn transfer_dc_limit(
    response: &impl Fn(f64) -> Result<C64, DrivingPointError>,
) -> Result<f64, DrivingPointError> {
    dc_limit(response, &DC_PROBE_RAD_S)
}

/// The mouth-flow ladder on a sampled transfer: `response` at the positive bins, `dc` (the true limit `T(0+)`) at bin
/// 0, `transit` the one-way delay in samples, `max_taps` the validated ladder cap. The DFT length starts at the
/// 4-transit power-of-two rung and doubles until the wrap-around tail and the change of the shared taps are both under
/// the ladder's tolerance, or refuses at the cap. A single duct and a network port use this same arithmetic.
pub(super) fn transfer_ladder(
    response: &impl Fn(f64) -> Result<C64, DrivingPointError>,
    dc: f64,
    transit: f64,
    dt: f64,
    max_taps: usize,
    ladder: &LineLadder,
) -> Result<(Vec<f64>, MouthFlowRealization), DrivingPointError> {
    let demand = (4.0 * transit).ceil();
    if !demand.is_finite() {
        return Err(DrivingPointError::Invalid {
            what: "transit delay is not finite",
        });
    }
    // Start at the 4-transit power-of-two rung, leaving at least one doubling.
    let mut n_fft = (demand.min(max_taps as f64) as usize)
        .max(256)
        .next_power_of_two()
        .min(max_taps);
    if max_taps >= 512 {
        n_fft = n_fft.min(max_taps / 2).max(256);
    }
    let tap_tol = ladder.tail_tol.sqrt();
    let mut t_pos: Vec<C64> = Vec::new();
    let mut previous: Vec<f64> = Vec::new();
    let mut rungs: Vec<MouthFlowRung> = Vec::new();
    loop {
        let half = n_fft / 2;
        let mut next = vec![C64::ZERO; half + 1];
        if t_pos.is_empty() {
            for (k, slot) in next.iter_mut().enumerate().skip(1) {
                *slot = response(core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt))?;
            }
        } else {
            // Even bins of this rung are the bins of the previous one.
            for (k_prev, &value) in t_pos.iter().enumerate() {
                next[2 * k_prev] = value;
            }
            for k in (1..half).step_by(2) {
                next[k] = response(core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt))?;
            }
        }
        next[0] = C64::new(dc, 0.0);
        t_pos = next;
        let ir = realize_reflectance_ir_from_positive_spectrum(&t_pos, n_fft, transit, true);
        let total: f64 = ir.iter().map(|h| h * h).sum();
        let tail: f64 = ir[3 * (n_fft / 4)..].iter().map(|h| h * h).sum();
        let tail_energy = if total > 0.0 { tail / total } else { 0.0 };
        let tap_change = if previous.is_empty() {
            f64::INFINITY
        } else {
            let diff: f64 = ir[..previous.len()]
                .iter()
                .zip(&previous)
                .map(|(a, b)| (a - b) * (a - b))
                .sum();
            if total > 0.0 { (diff / total).sqrt() } else { 0.0 }
        };
        rungs.push(MouthFlowRung { n_fft, tail_energy, tap_change });
        if (n_fft as f64) > transit + 2.0 && tail_energy < ladder.tail_tol && tap_change < tap_tol {
            return Ok((ir, MouthFlowRealization { n_fft, rungs }));
        }
        if n_fft >= max_taps {
            return Err(DrivingPointError::TransferLadderUnconverged {
                max_taps,
                tail_energy,
                tail_tol: ladder.tail_tol,
                tap_change,
                tap_tol,
            });
        }
        previous = ir[..half].to_vec();
        n_fft *= 2;
    }
}

/// A causal FIR line carrying incident inlet pressure to mouth volume flow:
/// `push(p+)` returns `U_mouth` for the same sample.
///
/// # Errors
/// As [`realize_mouth_flow_ir`], plus filter construction.
pub fn mouth_flow_line<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<(DelayedFilter, MouthFlowRealization), DrivingPointError> {
    let (ir, receipt) =
        realize_mouth_flow_ir(physics, gas, termination, sample_rate_hz, zc, wall, ladder)?;
    let line = DelayedFilter::from_impulse_response(1.0 / f64::from(sample_rate_hz), ir)
        .map_err(DrivingPointError::Discrete)?;
    Ok((line, receipt))
}

/// [`mouth_flow_line`] for a bore filled with `gas` that radiates into a different `exterior`: only a tabulated mouth
/// load is taken, checked against the `exterior` it was computed in. Every bore quantity (the transit, `zc`, the
/// propagation and its loss) stays `gas`'s, and the load is absolute, so nothing is converted at the mouth. With
/// `exterior` equal to `gas` this is [`mouth_flow_line`] on that table exactly.
///
/// # Errors
/// As [`mouth_flow_line`], or a table whose mouth or exterior gas does not match.
#[allow(clippy::too_many_arguments)]
pub fn mouth_flow_line_exterior(
    physics: &Duct,
    gas: &GasState,
    exterior: &GasState,
    table: &fs_duct::TabulatedLoad,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<(DelayedFilter, MouthFlowRealization), DrivingPointError> {
    let (ir, receipt) = realize_with_dc_probes(
        physics, gas, exterior, MouthLoad::Tabulated(table), sample_rate_hz, zc, wall, ladder, &DC_PROBE_RAD_S,
    )?;
    let line = DelayedFilter::from_impulse_response(1.0 / f64::from(sample_rate_hz), ir)
        .map_err(DrivingPointError::Discrete)?;
    Ok((line, receipt))
}

/// The two lines of a series pressure source `p_s` at a junction inside the duct (a constriction's turbulence, MB 15368),
/// with the inlet port held at `zc` exactly as the island's line and [`mouth_flow_line`] hold it, so the island's own
/// waves and the source's add by superposition on the linear duct.
#[derive(Debug, Clone)]
pub struct JunctionSourceLines {
    /// `push(p_s)` returns the mouth volume flow [m^3/s] of the same sample, to add to [`mouth_flow_line`]'s.
    pub to_mouth_flow: DelayedFilter,
    /// `push(p_s)` returns the wave [Pa] the source sends to the inlet port, referenced to `zc`, to add to the island's
    /// returned wave in the island's own convention (its line's output is read on the next step).
    pub to_inlet_wave: DelayedFilter,
    /// Receipts of the two ladders, mouth first.
    pub receipts: [MouthFlowRealization; 2],
}

/// `(U_mouth / p_s, p-_inlet / p_s)` at one angular frequency (`e^{-iwt}`) for a series source between segment `k` and
/// segment `k + 1`, raising the pressure on the mouth side. Looking upstream from the junction, the segments `k, ..., 0`
/// end in the port's `zc` (`Z_U`, inlet-flow transfer `T_U`); downstream, `k + 1, ...` end in the termination (`Z_D`,
/// `T_D`). The loop flow is `U = p_s / (Z_U + Z_D)` toward the mouth, so `U_mouth = T_D U`; the flow `-T_U U` reaches the
/// port, whose matched load turns it into the arriving wave `p- = -zc T_U U`.
///
/// # Errors
/// A TMM refusal on either side.
pub fn junction_source_transfer(
    upstream_reversed: &Duct,
    downstream: &Duct,
    gas: &GasState,
    omega: f64,
    termination: Termination,
    wall: Option<&WallPin>,
    zc: f64,
) -> Result<(C64, C64), DrivingPointError> {
    let up = fs_duct::input_impedance_load(upstream_reversed, gas, omega, LossModel::Bessel, C64::new(zc, 0.0), wall)
        .map_err(DrivingPointError::Duct)?;
    let down = input_impedance_wall(downstream, gas, omega, LossModel::Bessel, termination, wall)
        .map_err(DrivingPointError::Duct)?;
    let loop_impedance = up.impedance + down.impedance;
    let to_mouth = down.u_mouth_over_u_in / loop_impedance;
    let to_inlet = (up.u_mouth_over_u_in * C64::new(-zc, 0.0)) / loop_impedance;
    Ok((to_mouth, to_inlet))
}

/// [`junction_source_transfer`] with the mouth closed by the absolute load `z_load` [Pa s/m^3] (a tabulated load's value
/// at `omega`) in place of a declared termination: the same loop, with the downstream side ending in that load through
/// the seam the mouth-flow line plays a table through ([`mouth_flow_line`], `baked::response`).
///
/// # Errors
/// A TMM refusal on either side.
fn junction_source_transfer_against_load(
    upstream_reversed: &Duct,
    downstream: &Duct,
    gas: &GasState,
    omega: f64,
    z_load: C64,
    wall: Option<&WallPin>,
    zc: f64,
) -> Result<(C64, C64), DrivingPointError> {
    let up = fs_duct::input_impedance_load(upstream_reversed, gas, omega, LossModel::Bessel, C64::new(zc, 0.0), wall)
        .map_err(DrivingPointError::Duct)?;
    let down = super::baked::response(downstream, gas, omega, z_load, WallModel::from_pin(wall))
        .map_err(DrivingPointError::Duct)?;
    let loop_impedance = up.impedance + down.impedance;
    let to_mouth = down.u_mouth_over_u_in / loop_impedance;
    let to_inlet = (up.u_mouth_over_u_in * C64::new(-zc, 0.0)) / loop_impedance;
    Ok((to_mouth, to_inlet))
}

/// The upstream (reversed, segments `k, ..., 0`) and downstream (`k + 1, ...`) ducts of a junction after segment `k`.
///
/// # Errors
/// A junction at or after the last segment (nothing downstream), or a segment that is not a cylinder (reversal is only
/// exact for cylinders here), by name.
pub fn split_at_junction(physics: &Duct, k: usize) -> Result<(Duct, Duct), DrivingPointError> {
    if k + 1 >= physics.segments.len() {
        return Err(DrivingPointError::Invalid {
            what: "a junction source needs at least one segment downstream of its section",
        });
    }
    if !physics.segments.iter().all(|s| matches!(s, fs_duct::Segment::Cylinder { .. })) {
        return Err(DrivingPointError::Invalid {
            what: "a junction source needs a duct of cylinders",
        });
    }
    let upstream_reversed = Duct { segments: physics.segments[..=k].iter().rev().cloned().collect() };
    let downstream = Duct { segments: physics.segments[k + 1..].to_vec() };
    Ok((upstream_reversed, downstream))
}

/// [`JunctionSourceLines`] for a series source after segment `k` of `physics` (a duct of cylinders), on the same
/// transfer ladder, sampled bins and DC probes as [`mouth_flow_line`]; both lines start at the whole duct's 4-transit
/// rung, since the loop rings with the whole duct's modes. The mouth load is the island's own, taken as
/// [`mouth_flow_line`] takes it: a declared termination exactly as before, or a tabulated load checked against the whole
/// duct's mouth in `gas` and covering the top rung, with its DC bin against the table's declared static load.
///
/// # Errors
/// As [`split_at_junction`] and [`mouth_flow_line`].
#[allow(clippy::too_many_arguments)]
pub fn junction_source_lines<'a>(
    physics: &Duct,
    k: usize,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<JunctionSourceLines, DrivingPointError> {
    junction_lines_in(physics, k, gas, gas, termination.into(), sample_rate_hz, zc, wall, ladder)
}

/// [`junction_source_lines`] for a bore filled with `gas` that radiates into a different `exterior`, as
/// [`mouth_flow_line_exterior`]: only a tabulated load is taken, checked against the `exterior` it was computed in; every
/// bore quantity stays `gas`'s.
///
/// # Errors
/// As [`junction_source_lines`], or a table whose mouth or exterior gas does not match.
#[allow(clippy::too_many_arguments)]
pub fn junction_source_lines_exterior(
    physics: &Duct,
    k: usize,
    gas: &GasState,
    exterior: &GasState,
    table: &fs_duct::TabulatedLoad,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<JunctionSourceLines, DrivingPointError> {
    junction_lines_in(physics, k, gas, exterior, MouthLoad::Tabulated(table), sample_rate_hz, zc, wall, ladder)
}

#[allow(clippy::too_many_arguments)]
fn junction_lines_in(
    physics: &Duct,
    k: usize,
    gas: &GasState,
    exterior: &GasState,
    load: MouthLoad<'_>,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<JunctionSourceLines, DrivingPointError> {
    let (up, down) = split_at_junction(physics, k)?;
    let length = duct_axial_length(physics);
    if !(length > 0.0 && length.is_finite()) {
        return Err(DrivingPointError::Invalid {
            what: "a junction source needs positive axial length",
        });
    }
    if sample_rate_hz == 0 || !(zc.is_finite() && zc > 0.0) {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate and finite positive port impedance required",
        });
    }
    let max_taps = validate_ladder(ladder)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let transit = length / gas.sound_speed / dt;
    let demand = (4.0 * transit).ceil();
    if !demand.is_finite() {
        return Err(DrivingPointError::Invalid {
            what: "transit delay is not finite",
        });
    }
    // A declared termination takes exactly its old path; a table plays only through its own mouth (the whole duct's,
    // as the island's line and the mouth-flow line check it) and must cover the top rung before any bin is sampled.
    let static_load = match load {
        MouthLoad::Analytic(_) => None,
        MouthLoad::Tabulated(table) => {
            let radius = physics
                .segments
                .last()
                .ok_or(DrivingPointError::Duct(fs_duct::DuctError::EmptyDuct))?
                .outlet_radius();
            let mouth = table.validate_mouth(length, radius, exterior).map_err(DrivingPointError::Duct)?;
            super::baked::validate_band(table, sample_rate_hz, max_taps)?;
            Some(C64::new(mouth.dc_impedance, 0.0))
        }
        MouthLoad::Matched => return Err(super::baked::matched_refusal()),
    };
    let transfer = |omega: f64| match load {
        MouthLoad::Analytic(termination) => junction_source_transfer(&up, &down, gas, omega, termination, wall, zc),
        MouthLoad::Tabulated(table) => {
            let z = table.z_at(omega).map_err(DrivingPointError::Duct)?;
            junction_source_transfer_against_load(&up, &down, gas, omega, z, wall, zc)
        }
        MouthLoad::Matched => Err(super::baked::matched_refusal()),
    };
    let to_mouth = |omega: f64| transfer(omega).map(|t| t.0);
    let to_inlet = |omega: f64| transfer(omega).map(|t| t.1);
    let dc = |part: fn((C64, C64)) -> C64| match static_load {
        Some(z) => dc_limit(
            &|omega| junction_source_transfer_against_load(&up, &down, gas, omega, z, wall, zc).map(part),
            &DC_PROBE_RAD_S,
        ),
        None => dc_limit(&|omega| transfer(omega).map(part), &DC_PROBE_RAD_S),
    };
    // Both lines start at the whole duct's 4-transit rung: `transit` is the whole duct's.
    let (mouth_ir, mouth_receipt) = transfer_ladder(&to_mouth, dc(|t| t.0)?, transit, dt, max_taps, ladder)?;
    let (inlet_ir, inlet_receipt) = transfer_ladder(&to_inlet, dc(|t| t.1)?, transit, dt, max_taps, ladder)?;
    Ok(JunctionSourceLines {
        to_mouth_flow: DelayedFilter::from_impulse_response(dt, mouth_ir).map_err(DrivingPointError::Discrete)?,
        to_inlet_wave: DelayedFilter::from_impulse_response(dt, inlet_ir).map_err(DrivingPointError::Discrete)?,
        receipts: [mouth_receipt, inlet_receipt],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_duct::Segment;
    use fs_material::gas::GasSpec;

    const RATE: u32 = 48_000;

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    fn tube() -> (Duct, f64) {
        let radius = 0.01;
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length: 0.17 }],
        };
        let gas = air();
        let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        (duct, zc)
    }

    #[test]
    fn realized_fir_reproduces_the_sampled_mouth_flow_response_between_bins() {
        // The IR is built from samples on a DFT grid; at frequencies that are
        // NOT grid bins its transform must still equal the direct TMM response,
        // because the converged tail leaves nothing to alias. 1 and 7 Hz sit
        // below the first bins and would expose a DC bin that is not the true
        // T(0+): the first version took T(20 Hz), which missed the production
        // tail tolerance at 65536 taps (2.87e-10 against 1e-10).
        let (duct, zc) = tube();
        let gas = air();
        let (ir, receipt) = realize_mouth_flow_ir(
            &duct, &gas, Termination::FlangedOpen, RATE, zc, None, &LineLadder::default(),
        )
        .expect("ladder converges");
        assert!(receipt.rungs.len() >= 2, "a ladder needs at least two rungs");
        let dt = 1.0 / f64::from(RATE);
        let frequencies = [1.0, 7.0, 310.0, 1013.0, 2207.0, 3119.0];
        let mut worst = 0.0_f64;
        let mut scale = 0.0_f64;
        for f in frequencies {
            let omega = core::f64::consts::TAU * f;
            let direct = mouth_flow_per_incident_pressure(
                &duct, &gas, omega, Termination::FlangedOpen, None, zc,
            )
            .expect("direct response");
            let (mut re, mut im) = (0.0_f64, 0.0_f64);
            for (n, h) in ir.iter().enumerate() {
                let phase = omega * n as f64 * dt;
                re += h * phase.cos();
                im += h * phase.sin();
            }
            let error = ((re - direct.re).powi(2) + (im - direct.im).powi(2)).sqrt();
            worst = worst.max(error);
            scale = scale.max(direct.abs());
        }
        println!(
            "IRPVOICE_FIR_ERROR taps={} worst={worst:e} scale={scale:e} ratio={:e}",
            receipt.n_fft,
            worst / scale
        );
        assert!(
            worst < 1.0e-3 * scale,
            "FIR departs from the TMM mouth-flow response by {worst:e} of scale {scale:e}"
        );
    }

    // Measurement, not a gate on a number. An unreachable tail tolerance makes
    // each ladder run to its cap and report the tail energy it reached there,
    // so the tail can be read against tap count: exponential ringing falls by
    // orders of magnitude per doubling, an algebraic (boundary-layer) tail by
    // a fixed factor. Output lines are tagged IRPVOICE_TAIL.
    fn never_converges(cap: usize) -> LineLadder {
        LineLadder { tail_tol: 1.0e-300, peak_tol_cents: 0.05, max_taps: cap }
    }

    fn tail_of<T>(result: Result<T, DrivingPointError>) -> f64 {
        match result {
            Err(
                DrivingPointError::LadderUnconverged { tail_energy, .. }
                | DrivingPointError::TransferLadderUnconverged { tail_energy, .. },
            ) => tail_energy,
            Err(other) => panic!("unexpected refusal: {other}"),
            Ok(_) => panic!("a ladder with an unreachable tolerance cannot converge"),
        }
    }

    /// A transfer ladder that reaches its ceiling prints what it gated on, the measured tap change and its tolerance.
    #[test]
    fn a_transfer_ladder_refusal_prints_the_tap_change() {
        // One pole decaying over 2400 samples: the 1024-tap ceiling cannot hold it.
        let sigma = 20.0;
        let response = |w: f64| Ok(C64::new(1.0, 0.0) / C64::new(sigma, -w));
        let ladder = LineLadder { tail_tol: 1.0e-10, peak_tol_cents: 0.05, max_taps: 1024 };
        let refusal = transfer_ladder(&response, 1.0 / sigma, 8.0, 1.0 / f64::from(RATE), 1024, &ladder)
            .map(|_| ())
            .expect_err("the ceiling is far shorter than the decay");
        let message = refusal.to_string();
        assert!(message.contains("tap_change=") && !message.contains("NaN"), "{message}");
        match refusal {
            DrivingPointError::TransferLadderUnconverged { max_taps, tail_tol, tap_change, tap_tol, .. } => {
                assert_eq!(max_taps, 1024);
                assert_eq!(tap_tol, tail_tol.sqrt());
                assert!(tap_change.is_finite() && tap_change > 0.0, "tap change {tap_change}");
            }
            other => panic!("unexpected refusal: {other}"),
        }
    }

    fn tails(
        name: &str, duct: &Duct, gas: &GasState, rate: u32, zc: f64, wall: Option<&WallPin>,
    ) -> Vec<(f64, f64)> {
        let mut series = Vec::new();
        for cap in [4096_usize, 16384, 65536] {
            let ladder = never_converges(cap);
            let reflectance = tail_of(crate::driving_point::realize_line(
                duct, gas, Termination::FlangedOpen, rate, zc, wall, &ladder,
            ));
            let mouth_flow = tail_of(realize_mouth_flow_ir(
                duct, gas, Termination::FlangedOpen, rate, zc, wall, &ladder,
            ));
            println!(
                "IRPVOICE_TAIL case={name} rate_hz={rate} taps={cap} \
                 reflectance_tail={reflectance:e} mouth_flow_tail={mouth_flow:e}"
            );
            series.push((reflectance, mouth_flow));
        }
        series
    }

    fn description_case(file: &str) -> (Duct, GasState, u32, f64, Option<WallPin>) {
        let path = format!(
            "{}/../../data/instruments/descriptions/voice/{file}",
            env!("CARGO_MANIFEST_DIR")
        );
        let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let description = crate::voice_performance::VoiceDescription::from_bytes(&bytes)
            .expect("description parses");
        let chart = crate::tract::TractChart::try_new(
            description.sections.clone(),
            &description.tract_source,
            &description.tract_licence,
        )
        .expect("chart admitted");
        let rate = description.output_rate_hz * description.mechanics_per_output;
        let gas = *description.ambient_gas();
        let zc = gas.density * gas.sound_speed / description.sections[0].area_m2;
        (chart.to_duct(), gas, rate, zc, description.wall)
    }

    #[test]
    fn a_dc_bin_taken_from_20_hz_leaves_a_tail_the_true_limit_does_not() {
        // Controlled A/B on the same tube at the same cap: the only difference
        // is where the DC bin is read. A jump J between bin 0 and bin 1 adds a
        // constant J/N to every tap, so its tail energy falls only as 1/N.
        let (duct, zc) = tube();
        let gas = air();
        let ladder = never_converges(65536);
        let tail = |probes: &[f64]| {
            tail_of(realize_with_dc_probes(
                &duct, &gas, &gas, Termination::FlangedOpen, RATE, zc, None, &ladder, probes,
            ))
        };
        let stepped = tail(&[core::f64::consts::TAU * 20.0]);
        let limit = tail(&DC_PROBE_RAD_S);
        println!("IRPVOICE_DC_AB taps=65536 stepped_dc_tail={stepped:e} true_limit_tail={limit:e} ratio={:e}", stepped / limit);
        assert!(stepped.is_finite() && limit.is_finite());
        assert!(
            limit <= stepped,
            "reading the DC bin at its true limit must not leave a longer tail: {limit:e} vs {stepped:e}"
        );
    }

    #[test]
    fn tail_energy_against_tap_count_for_both_lines() {
        let (tube_duct, tube_zc) = tube();
        let mut cases = vec![(
            "uniform-tube".to_string(),
            tails("uniform-tube", &tube_duct, &air(), RATE, tube_zc, None),
        )];
        for file in [
            "story2008-a-Ds3.performance",
            "story2008-u-Ds3.performance",
            "story2008-u-Ds3-r1.performance",
        ] {
            let (duct, gas, rate, zc, wall) = description_case(file);
            cases.push((file.to_string(), tails(file, &duct, &gas, rate, zc, wall.as_ref())));
        }
        for (name, series) in &cases {
            for (reflectance, mouth_flow) in series {
                assert!(
                    reflectance.is_finite() && mouth_flow.is_finite()
                        && *reflectance >= 0.0 && *mouth_flow >= 0.0,
                    "{name}: tail energies must be finite and non-negative"
                );
            }
            let (first, last) = (series[0], series[series.len() - 1]);
            assert!(
                last.0 < first.0 && last.1 < first.1,
                "{name}: tail energy must fall with tap count, got {first:?} then {last:?}"
            );
        }
    }

    #[test]
    fn flow_gain_matches_the_inlet_flow_in_the_low_frequency_limit() {
        // Below the first resonance the duct is compact: the flow leaving the
        // mouth equals the flow entering the inlet, so
        // |T| = |U_in / p+| = 2 / |Z + zc| with |U_mouth / U_in| -> 1.
        let (duct, zc) = tube();
        let gas = air();
        let omega = core::f64::consts::TAU * 20.0;
        let t = mouth_flow_per_incident_pressure(&duct, &gas, omega, Termination::FlangedOpen, None, zc)
            .expect("response");
        let z = input_impedance_wall(&duct, &gas, omega, LossModel::Bessel, Termination::FlangedOpen, None)
            .expect("impedance");
        let inlet = C64::new(2.0, 0.0) / C64::new(z.impedance.re + zc, z.impedance.im);
        assert!(
            (t.abs() - inlet.abs()).abs() < 5.0e-3 * inlet.abs(),
            "|T| {} versus inlet flow gain {}",
            t.abs(),
            inlet.abs()
        );
    }

    /// A series source after segment `k` of a uniform tube, inlet held at `zc`: the inlet wave cannot arrive before the
    /// upstream transit `x_k / c` and peaks there (the direct wave; every echo has further to go and a mouth reflection to
    /// lose), the mouth flow likewise after `(L - x_k) / c`. The share of each line's energy more than 4 samples before its
    /// transit is no more than the share the realization's own kernel puts there for a unit delay at the same transit and
    /// rung (MB 16979, 16990, 17012): the kernel's Nyquist rolloff and precursor fold give even a perfect delay early taps,
    /// so the bound is computed from the kernel, never from the line, and a wave arriving early would put a share of its
    /// whole energy there. A share is free of the line's units (the mouth flow's `1 / zc`). On a tube whose transits are
    /// whole numbers of samples, and on one whose inlet transit is 9.5. At DC the loop is the port's resistor, so the inlet
    /// wave is `-p_s` and the mouth flow `p_s / zc`, to the tube's viscous share. A junction with nothing downstream and a
    /// duct that is not all cylinders are refused by name.
    #[test]
    fn junction_source_lines_are_causal_from_the_junction_and_see_the_port_resistor_at_dc() {
        let radius = 0.01;
        let gas = air();
        let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let samples = |metres: f64| metres / gas.sound_speed * f64::from(RATE);
        let tube = |segment: f64| Duct { segments: vec![Segment::Cylinder { radius, length: segment }; 10] };
        let lines_of = |duct: &Duct| {
            junction_source_lines(duct, 3, &gas, Termination::FlangedOpen, RATE, zc, None, &LineLadder::default())
                .expect("junction lines")
        };
        // The share of its energy more than 4 samples before `transit` (taps `..floor(transit - 4)`).
        let early_share = |taps: &[f64], transit: f64| {
            let total: f64 = taps.iter().map(|h| h * h).sum();
            taps[..(transit - 4.0).floor() as usize].iter().map(|h| h * h).sum::<f64>() / total
        };
        // A unit delay of `transit` samples realized as every transfer line of this module is at rung `n_fft` (acoustic
        // e^{-iwt} bins, as `transfer_ladder` passes them).
        let kernel = |transit: f64, n_fft: usize| {
            let delay: Vec<C64> = (0..=n_fft / 2)
                .map(|k| {
                    let phase = core::f64::consts::TAU * k as f64 / n_fft as f64 * transit;
                    C64::new(phase.cos(), phase.sin())
                })
                .collect();
            realize_reflectance_ir_from_positive_spectrum(&delay, n_fft, transit, true)
        };
        let whole = 2.5 * gas.sound_speed / f64::from(RATE);
        for segment in [whole, 0.017] {
            let lines = lines_of(&tube(segment));
            let (_, mouth) = lines.to_mouth_flow.quiescent_fir().expect("mouth FIR");
            let (_, inlet) = lines.to_inlet_wave.quiescent_fir().expect("inlet FIR");
            for (name, taps, transit, rung) in [
                ("inlet", inlet, samples(4.0 * segment), lines.receipts[1].n_fft),
                ("mouth", mouth, samples(6.0 * segment), lines.receipts[0].n_fft),
            ] {
                let (share, bound) = (early_share(taps, transit), early_share(&kernel(transit, rung), transit));
                assert!(
                    share <= bound,
                    "{name}: {share:e} of its energy more than 4 samples before the transit {transit:.3}, above the \
                     kernel's own {bound:e} for a unit delay at rung {rung}"
                );
                let peak = taps.iter().enumerate().max_by(|a, b| a.1.abs().total_cmp(&b.1.abs())).expect("taps").0;
                assert!((peak as f64 - transit).abs() < 4.0, "{name}: peak at tap {peak}, transit {transit:.3} samples");
            }
        }
        let duct = tube(0.017);
        let lines = lines_of(&duct);
        let (_, mouth) = lines.to_mouth_flow.quiescent_fir().expect("mouth FIR");
        let (_, inlet) = lines.to_inlet_wave.quiescent_fir().expect("inlet FIR");
        let (inlet_dc, mouth_dc): (f64, f64) = (inlet.iter().sum(), mouth.iter().sum());
        assert!((inlet_dc + 1.0).abs() < 1.0e-2, "inlet wave at DC {inlet_dc}, the port's resistor gives -1");
        assert!((mouth_dc * zc - 1.0).abs() < 1.0e-2, "mouth flow at DC {mouth_dc}, the port's resistor gives 1 / zc");
        let refused = |duct: &Duct, k: usize| match junction_source_lines(duct, k, &gas, Termination::FlangedOpen, RATE, zc, None, &LineLadder::default()) {
            Err(DrivingPointError::Invalid { what }) => what,
            other => panic!("expected a refusal by name, got {:?}", other.map(|_| ())),
        };
        assert!(refused(&duct, 9).contains("downstream"));
        let mut coned = duct.clone();
        coned.segments[5] = Segment::Cone { inlet_radius: radius, outlet_radius: 1.2 * radius, length: 0.017 };
        assert!(refused(&coned, 3).contains("cylinders"));
    }
}
