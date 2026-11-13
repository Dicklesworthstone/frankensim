//! Physical output harmonics of a full-state periodic pHS waveform.
//!
//! State coefficients are not output coefficients when storage is nonlinear:
//! `y(t) = G^T grad H(x(t))` can contain DC, even harmonics and harmonics above
//! the retained state order. This observer evaluates that original law, then
//! uses the shared FFT owner to project its output. It neither linearizes H
//! nor invents a harmonic envelope. The state and port order stay unchanged.

use super::port_flows_into;
use fs_fft::{C64 as FftComplex, Fft};
use fs_math::{c64::C64, det};
use fs_orbit::HbOrbit;
use fs_phs::PortHamiltonian;

/// Physical output accuracy and explicit cold-observation resource limits.
#[derive(Clone, Debug)]
pub struct PortSpectrumControls {
    /// Number of positive output harmonics, independent of the state order.
    pub harmonics: usize,
    /// First power-of-two grid, at least 32 and strictly above twice both
    /// the state and output orders. No Nyquist state is silently discarded.
    pub initial_samples: usize,
    /// Largest power-of-two grid, at most 65536 and at least four times the
    /// initial grid. At least two successive grid comparisons must pass.
    pub max_samples: usize,
    /// Maximum total evaluations of the original storage gradient. Its own
    /// internal cost is source-defined. FFT and port work additionally scale
    /// with the original state/port counts and admitted sample grids.
    pub max_source_evaluations: usize,
    /// Conservative peak numeric-array bytes, including FFT plans/scratch,
    /// sampled states/ports, current/previous spectra and the returned trace.
    /// Excludes the caller-owned system and orbit and allocator metadata.
    pub max_workspace_bytes: usize,
    /// Positive absolute RMS tolerance per physical port, in that output's
    /// units. Sampling change and observed omitted-band RMS SHARE this budget.
    pub rms_tolerances: Vec<f64>,
}

/// Finite-grid evidence, not an enclosure of the infinite harmonic tail.
#[derive(Clone, Debug, PartialEq)]
pub struct PortSpectrumCheck {
    /// Complete original-storage evaluations on this grid.
    pub samples: usize,
    /// RMS in measured output bins above the requested order, per port.
    /// Includes the Nyquist bin once and each interior conjugate pair twice.
    pub omitted_rms: Vec<f64>,
    /// RMS change of the retained DC and complex coefficients from the prior
    /// grid. None on the first grid. Every grid uses the SAME phase origin.
    pub sampling_change_rms: Option<Vec<f64>>,
}

/// Physical port outputs ready for existing negative-time harmonic consumers.
#[derive(Clone, Debug)]
pub struct PortSpectrum {
    /// Original orbit fundamental [rad/s].
    pub omega_rad_s: f64,
    /// Physical DC output per original port. DC is real and is never doubled.
    pub dc: Vec<f64>,
    /// `[port][harmonic - 1]` PEAK phasors under `exp(-i n omega t)`.
    /// These are `2 conj(c_n)` of the positive-time Fourier coefficients.
    /// Mechanical force, velocity, current, pressure, etc. retain their source
    /// units; a consumer must use the corresponding physical coupling map.
    pub peak_phasors: Vec<Vec<C64>>,
    /// All evaluated grids, including two accepted successive comparisons.
    pub refinement_trace: Vec<PortSpectrumCheck>,
    /// Total evaluations of the original nonlinear storage gradient.
    pub source_evaluations: usize,
}

/// Observation refuses atomically; it never changes the original system/orbit.
#[derive(Clone, Debug, PartialEq)]
pub enum PortSpectrumError {
    Invalid(&'static str),
    Budget(&'static str),
    SourceRefused,
    Cancelled,
    /// The admitted ladder did not resolve both sampling and output truncation.
    Unresolved {
        last: PortSpectrumCheck,
    },
}
impl core::fmt::Display for PortSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical port spectrum: {self:?}")
    }
}
impl std::error::Error for PortSpectrumError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PortSpectrumError> {
    if cancelled() {
        Err(PortSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Evaluate all original physical port outputs, with measured grid refinement.
///
/// FFT synthesis reconstructs every retained state, including material memory,
/// before calling the original gradient once per sample. Output projection is
/// O((states + ports) N log N), and no trigonometric loop over every state,
/// sample and harmonic is needed. Successive grids double N and shift their
/// sample nodes; projection removes that shift so all reported phases refer
/// to the original orbit clock. This helps expose folded nonlinear harmonics.
///
/// Two successive comparisons must satisfy, on EVERY port,
/// `sampling_change_rms + omitted_rms <= rms_tolerance`. This is a measured
/// finite-grid criterion: arbitrary unresolved content can still evade finite
/// sampling. It proves neither the orbit's equation residual, convergence in
/// state order, stability, a continuum error bound nor an infinite-tail bound.
/// It is equally applicable to forced and autonomous candidate waveforms.
/// Source gradients must obey their original deterministic/read-only contract.
///
/// Cancellation is checked before each state FFT, source evaluation and output
/// FFT; an individual FFT has at most the admitted 65536 nodes. Failed work
/// publishes no partial output spectrum. DC is retained separately because
/// positive-frequency acoustic radiation APIs must not silently radiate it.
pub fn observe_port_spectrum(
    system: &PortHamiltonian,
    orbit: &HbOrbit,
    controls: &PortSpectrumControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PortSpectrum, PortSpectrumError> {
    poll(&mut cancelled)?;
    admit(system, orbit, controls)?;
    let mut previous: Option<Vec<Vec<C64>>> = None;
    let mut trace = Vec::new();
    let mut samples = controls.initial_samples;
    let mut evaluations = 0_usize;
    let mut consecutive = 0;
    loop {
        poll(&mut cancelled)?;
        evaluations = evaluations
            .checked_add(samples)
            .filter(|&n| n <= controls.max_source_evaluations)
            .ok_or(PortSpectrumError::Budget(
                "original storage evaluation budget exhausted",
            ))?;
        admit_storage(system, controls, samples)?;
        // A third-of-a-cell shift changes nodes without changing the clock.
        let shift = core::f64::consts::TAU * (trace.len() % 3) as f64 / (3.0 * samples as f64);
        let (coefficients, omitted_rms) = project(
            system,
            orbit,
            controls.harmonics,
            samples,
            shift,
            &mut cancelled,
        )?;
        let change = previous.as_ref().map(|old| {
            coefficients
                .iter()
                .zip(old)
                .map(|(new, old)| coefficient_change(new, old))
                .collect::<Vec<_>>()
        });
        let resolved = change.as_ref().is_some_and(|change| {
            change
                .iter()
                .zip(&omitted_rms)
                .zip(&controls.rms_tolerances)
                .all(|((&delta, &tail), &tolerance)| delta + tail <= tolerance)
        });
        consecutive = if resolved { consecutive + 1 } else { 0 };
        let check = PortSpectrumCheck {
            samples,
            omitted_rms,
            sampling_change_rms: change,
        };
        trace.push(check.clone());
        if consecutive >= 2 {
            let dc = coefficients.iter().map(|row| row[0].re).collect();
            let mut peak_phasors = Vec::with_capacity(system.port_dim());
            for row in &coefficients {
                let peaks: Vec<_> = row[1..].iter().map(|z| z.conj().scale(2.0)).collect();
                if peaks.iter().any(|&z| !finite(z)) {
                    return Err(PortSpectrumError::Invalid("physical peak phasor overflow"));
                }
                peak_phasors.push(peaks);
            }
            poll(&mut cancelled)?;
            return Ok(PortSpectrum {
                omega_rad_s: orbit.omega,
                dc,
                peak_phasors,
                refinement_trace: trace,
                source_evaluations: evaluations,
            });
        }
        if samples == controls.max_samples {
            return Err(PortSpectrumError::Unresolved { last: check });
        }
        previous = Some(coefficients);
        samples *= 2;
    }
}

fn admit(
    system: &PortHamiltonian,
    orbit: &HbOrbit,
    c: &PortSpectrumControls,
) -> Result<(), PortSpectrumError> {
    let n = system.state_dim();
    let m = system.port_dim();
    let width = orbit.coeffs.first().map_or(0, Vec::len);
    let highest = width.saturating_sub(1).max(c.harmonics);
    let (_, _, g) = system.structure();
    if n == 0
        || m == 0
        || width == 0
        || orbit.coeffs.len() != n
        || !orbit.omega.is_finite()
        || orbit.omega <= 0.0
        || !((highest.max(1) as f64) * orbit.omega).is_finite()
        || orbit
            .coeffs
            .iter()
            .any(|row| row.len() != width || row[0].im != 0.0 || row.iter().any(|&z| !finite(z)))
        || n.checked_mul(m) != Some(g.len())
        || g.iter().any(|v| !v.is_finite())
        || c.rms_tolerances.len() != m
        || c.rms_tolerances.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || c.initial_samples < 32
        || !c.initial_samples.is_power_of_two()
        || !c.max_samples.is_power_of_two()
        || c.max_samples > 65_536
        || c.initial_samples
            .checked_mul(4)
            .is_none_or(|v| v > c.max_samples)
        || highest
            .checked_mul(2)
            .is_none_or(|v| v >= c.initial_samples)
    {
        return Err(PortSpectrumError::Invalid(
            "complete finite original-state rows, positive per-port tolerances and three resolved sample grids required",
        ));
    }
    Ok(())
}

fn admit_storage(
    system: &PortHamiltonian,
    c: &PortSpectrumControls,
    samples: usize,
) -> Result<(), PortSpectrumError> {
    // Sample matrices, FFT plans and complex scratch; three spectra and the
    // bounded (at most 12-grid) evidence trace plus gradient/port scratch.
    let n = system.state_dim() as u128;
    let m = system.port_dim() as u128;
    let bytes =
        8 * (n + m + 12) * samples as u128 + 48 * m * (c.harmonics as u128 + 1) + 8 * (n + 32 * m);
    if bytes > c.max_workspace_bytes as u128 {
        return Err(PortSpectrumError::Budget(
            "physical spectrum numeric workspace exceeds its byte budget",
        ));
    }
    Ok(())
}

fn coefficient_change(new: &[C64], old: &[C64]) -> f64 {
    let mut rms = (new[0].re - old[0].re).abs();
    for (a, b) in new[1..].iter().zip(&old[1..]) {
        let difference = *a - *b;
        // Count the conjugate pair without squaring physical amplitudes.
        rms = det::hypot(rms, difference.abs());
        rms = det::hypot(rms, difference.abs());
    }
    rms
}

fn project(
    system: &PortHamiltonian,
    orbit: &HbOrbit,
    harmonics: usize,
    samples: usize,
    shift: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<Vec<C64>>, Vec<f64>), PortSpectrumError> {
    let n = system.state_dim();
    let m = system.port_dim();
    let fft = Fft::new(samples);
    let mut data = vec![FftComplex::default(); samples];
    let mut scratch = data.clone();
    let mut states = vec![0.0; n * samples];
    let mut outputs = vec![0.0; m * samples];
    for (coordinate, row) in orbit.coeffs.iter().enumerate() {
        poll(cancelled)?;
        data.fill(FftComplex::default());
        data[0].re = row[0].re;
        for (k, &coefficient) in row.iter().enumerate().skip(1) {
            let angle = k as f64 * shift;
            let c = coefficient * C64::new(det::cos(angle), det::sin(angle));
            // Forward FFT of the conjugated positive-frequency coefficients
            // synthesizes x without an overflow-prone N/inverse-N round trip.
            data[k] = FftComplex::new(c.re, -c.im);
            data[samples - k] = FftComplex::new(c.re, c.im);
        }
        fft.forward(&mut data, &mut scratch);
        for (j, value) in data.iter().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(PortSpectrumError::Invalid(
                    "original state waveform overflow",
                ));
            }
            states[j * n + coordinate] = value.re;
        }
    }
    let mut effort = vec![0.0; n];
    let mut flows = vec![0.0; m];
    for j in 0..samples {
        poll(cancelled)?;
        effort.fill(f64::NAN);
        system
            .storage()
            .gradient(&states[j * n..(j + 1) * n], &mut effort);
        if effort.iter().any(|v| !v.is_finite()) || !port_flows_into(system, &effort, &mut flows) {
            return Err(PortSpectrumError::SourceRefused);
        }
        for p in 0..m {
            outputs[p * samples + j] = flows[p];
        }
    }
    let mut coefficients = Vec::with_capacity(m);
    let mut tails = Vec::with_capacity(m);
    for values in outputs.chunks_exact(samples) {
        poll(cancelled)?;
        let scale = values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let scale = if scale == 0.0 { 1.0 } else { scale };
        for (z, &value) in data.iter_mut().zip(values) {
            *z = FftComplex::new((value / scale) / samples as f64, 0.0);
        }
        fft.forward(&mut data, &mut scratch);
        let mut row = Vec::with_capacity(harmonics + 1);
        let mut tail = 0.0;
        for (k, &z) in data[..=samples / 2].iter().enumerate() {
            let angle = -(k as f64) * shift;
            let c = C64::new(z.re, z.im).scale(scale) * C64::new(det::cos(angle), det::sin(angle));
            if !finite(c) {
                return Err(PortSpectrumError::Invalid(
                    "physical output Fourier coefficient overflow",
                ));
            }
            if k <= harmonics {
                row.push(if k == 0 { C64::from_re(c.re) } else { c });
            } else {
                tail = det::hypot(tail, c.abs());
                if k != samples / 2 {
                    tail = det::hypot(tail, c.abs());
                }
            }
        }
        if !tail.is_finite() {
            return Err(PortSpectrumError::Invalid(
                "physical output RMS is not representable",
            ));
        }
        coefficients.push(row);
        tails.push(tail);
    }
    Ok((coefficients, tails))
}
