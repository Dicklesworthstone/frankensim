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
use fs_duct::{Duct, LossModel, Termination, input_impedance_wall};
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

/// Causal impulse response of [`mouth_flow_per_incident_pressure`], sized by a
/// convergence ladder: the DFT length doubles until the wrap-around tail and
/// the change of the shared taps between rungs are both below the ladder's
/// tolerance (`tail_tol` and its square root), or the ladder refuses.
///
/// # Errors
/// Empty geometry, non-positive rate or port impedance, TMM refusal, or
/// non-convergence at `ladder.max_taps`.
pub fn realize_mouth_flow_ir(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<(Vec<f64>, MouthFlowRealization), DrivingPointError> {
    realize_with_dc_probes(
        physics, gas, termination, sample_rate_hz, zc, wall, ladder, &DC_PROBE_RAD_S,
    )
}

// The DC probe list is a parameter only so a test can rebuild the stepped-DC
// variant beside the true-limit one and measure what the step costs.
#[allow(clippy::too_many_arguments)]
fn realize_with_dc_probes(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
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
    // Start at the 4-transit power-of-two rung, leaving at least one doubling.
    let mut n_fft = (demand.min(max_taps as f64) as usize)
        .max(256)
        .next_power_of_two()
        .min(max_taps);
    if max_taps >= 512 {
        n_fft = n_fft.min(max_taps / 2).max(256);
    }
    let tap_tol = ladder.tail_tol.sqrt();
    let response = |omega: f64| {
        mouth_flow_per_incident_pressure(physics, gas, omega, termination, wall, zc)
    };
    let dc = dc_limit(&response, dc_probes)?;
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
            return Err(DrivingPointError::LadderUnconverged {
                max_taps,
                tail_energy,
                tail_tol: ladder.tail_tol,
                // Peak shift is a reflectance-line criterion; the transfer
                // ladder gates on the tap change instead.
                peak_shift_cents: f64::NAN,
                peak_tol_cents: ladder.peak_tol_cents,
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
pub fn mouth_flow_line(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
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
            Err(DrivingPointError::LadderUnconverged { tail_energy, .. }) => tail_energy,
            Err(other) => panic!("unexpected refusal: {other}"),
            Ok(_) => panic!("a ladder with an unreachable tolerance cannot converge"),
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
                &duct, &gas, Termination::FlangedOpen, RATE, zc, None, &ladder, probes,
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
}
