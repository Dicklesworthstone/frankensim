//! Delayed free-field mutual radiation between a flue's compact outlets.
//!
//! Each outlet's exterior pressure gains the free-field pressure of every
//! OTHER outlet's volume flow, p_m(t) = sum_n rho/(4 pi d_mn) dq_n/dt(t - d_mn/c)
//! (the off-diagonal of Z_mn = j w rho exp(-j k d)/(4 pi d)), realized
//! explicitly from accepted flows:
//!
//! p_m[k] = sum_n rho/(4 pi d) [s_n(k-N) - s_n(k-N-2)]/(2 dt),
//! s_n(j) = (1-f) q_n[j] + f q_n[j-1], N + f = d/(c dt) - 1.
//!
//! A central difference (zero at Nyquist, one sample of group delay) on a
//! linearly interpolated delay. With N >= 1 it reads only accepted flows, so
//! it enters each outlet's own terminal solve as a KNOWN series pressure and
//! no simultaneous multiport solve exists. The exact continuous term is not
//! passive together with the realized compact self loads (above about 120 kHz
//! at 288 kHz its resistance outgrows theirs); the central difference's
//! sin(w dt)/dt roll-off is the band limit, and passivity is checked, not
//! assumed: before any sample, Re Z = diag(realized self) + Re M must have a
//! Cholesky factor at every 5 Hz from 5 Hz to Nyquist, and M must stay within
//! 5 % of the exact term over the receiver's band.
//!
//! Not modeled: scattering by the instrument body (the outlets sit in the
//! card's flat plane at their declared centers, in free field), the image of
//! a baffle, the mouth (it is not an outlet), and any mean flow.
use crate::flue_loop::FlueError;
use fs_math::c64::C64;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationSampling};
use fs_vfit::waveguide::WaveguideError;
use fs_vfit::waveguide::network::WaveguideNetwork;

/// Admission grid step [Hz].
pub const OUTLET_MUTUAL_GRID_HZ: f64 = 5.0;
/// Largest admitted relative error of the realized mutual term in band.
pub const OUTLET_MUTUAL_BAND_ERROR: f64 = 0.05;
/// Largest number of coupled outlets.
pub const MAX_MUTUAL_OUTLETS: usize = 64;
/// Largest retained flow history per outlet [samples].
pub const MAX_MUTUAL_HISTORY: usize = 1 << 16;
// Largest admitted whole delay, N + 3 <= MAX_MUTUAL_HISTORY - 1.
const MAX_DELAY: f64 = 65_532.0;

fn invalid(what: &'static str) -> FlueError { FlueError::Invalid { what } }

/// One coupled pair, in the declared outlet order (`outlets[0] < outlets[1]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutletMutualPair {
    /// Positions in the declared outlet list.
    pub outlets: [usize; 2],
    /// Center distance in the card's plane [m].
    pub distance_m: f64,
    /// Whole samples N of the interpolated delay, N >= 1.
    pub delay_samples: usize,
    /// Fraction f in [0, 1) of the interpolated delay.
    pub delay_fraction: f64,
    /// rho / (4 pi d) / (2 dt) [Pa s/m^3 per sample].
    gain: f64,
}

/// Cold admission numbers, retained for the render's record.
#[derive(Debug, Clone, PartialEq)]
pub struct OutletMutualRecord {
    /// Exterior density of the self loads [kg/m^3].
    pub density_kg_m3: f64,
    /// Exterior sound speed of the self loads [m/s].
    pub sound_speed_m_s: f64,
    /// Mechanical step [s].
    pub time_step_s: f64,
    /// Receiver band the 5 % error gate covers [Hz].
    pub band_hz: f64,
    /// Number of admission frequencies (5 Hz to below Nyquist).
    pub grid_points: usize,
    /// Smallest whole delay over the pairs [samples].
    pub min_delay_samples: usize,
    /// Largest relative error of the realized term against the exact one in band.
    pub max_band_error: f64,
    /// Frequency of that error [Hz].
    pub max_band_error_hz: f64,
    /// Smallest Cholesky pivot over its own diagonal entry, over the grid.
    pub min_pivot_ratio: f64,
    /// Frequency of that pivot [Hz].
    pub min_pivot_hz: f64,
}

/// Explicit delayed mutual radiation: one flow ring per outlet.
#[derive(Debug, Clone)]
pub struct OutletMutual {
    nodes: Vec<usize>,
    pairs: Vec<OutletMutualPair>,
    len: usize,
    head: usize,
    rings: Vec<f64>,
    scratch: Vec<f64>,
    record: OutletMutualRecord,
}

/// Bilinear-realized self resistance Re Z(w) of an admitted relaxation load.
fn self_resistance(load: &RelaxationImpedanceSpec, warped: f64) -> f64 {
    let w2 = warped * warped;
    let mut r = load.base().resistance_pa_s_m3;
    for term in load.terms() {
        let p2 = term.rate_per_s * term.rate_per_s;
        r += term.resistance_pa_s_m3 * (w2 / (w2 + p2));
    }
    for term in load.compliance_terms() {
        let p2 = term.rate_per_s * term.rate_per_s;
        r += term.resistance_pa_s_m3 * (p2 / (w2 + p2));
    }
    r
}

impl OutletMutualPair {
    /// The realized kernel at angular frequency `omega`, under exp(+j w t):
    /// gain (1 - z^-2) z^-N ((1-f) + f z^-1), z^-1 = exp(-j w dt).
    #[must_use]
    pub fn realized(&self, omega: f64, dt: f64) -> C64 {
        let theta = omega * dt;
        let (s1, c1) = fs_math::det::sin_cos(theta);
        let (s2, c2) = fs_math::det::sin_cos(2.0 * theta);
        #[allow(clippy::cast_precision_loss)] // delays are bounded by MAX_MUTUAL_HISTORY
        let (sn, cn) = fs_math::det::sin_cos(self.delay_samples as f64 * theta);
        let difference = C64::new(1.0 - c2, s2);
        let delay = C64::new(cn, -sn);
        let interp = C64::new((1.0 - self.delay_fraction) + self.delay_fraction * c1, -self.delay_fraction * s1);
        (difference * delay * interp).scale(self.gain)
    }

    /// The exact free-field term j w rho exp(-j w d/c)/(4 pi d).
    #[must_use]
    pub fn exact(&self, omega: f64, density: f64, sound_speed: f64) -> C64 {
        let (s, c) = fs_math::det::sin_cos(omega * self.distance_m / sound_speed);
        let k = omega * density / (4.0 * core::f64::consts::PI * self.distance_m);
        C64::new(k * s, k * c)
    }
}

impl OutletMutual {
    /// Admit the coupling of these outlets (terminal `nodes`, plane `centers`,
    /// realized self `loads`, all in the declared outlet order) for exterior
    /// `density` and `sound_speed`, step `dt` and receiver band `band_hz`.
    ///
    /// # Errors
    /// Fewer than two or more than 64 outlets, mismatched lists, a repeated
    /// node, nonfinite or nonpositive inputs, coincident centers, a pair closer
    /// than two samples of travel (not explicit), a load not sampled
    /// bilinear-midpoint, a history beyond 65536 samples, a band at or above
    /// Nyquist, a realized term more than 5 % from the exact one in band, or
    /// a Re Z without a Cholesky factor at some grid frequency (not passive).
    pub fn admit(nodes: &[usize], centers: &[[f64; 2]], loads: &[RelaxationImpedanceSpec],
        density: f64, sound_speed: f64, dt: f64, band_hz: f64, gate: &fs_exec::CancelGate,
    ) -> Result<Self, FlueError> {
        let n = nodes.len();
        if n < 2 || n > MAX_MUTUAL_OUTLETS {
            return Err(invalid("outlet_mutual couples two to 64 outlets"));
        }
        if centers.len() != n || loads.len() != n {
            return Err(invalid("outlet_mutual needs one center and one self load per outlet"));
        }
        if (0..n).any(|i| nodes[i + 1..].contains(&nodes[i])) {
            return Err(invalid("outlet_mutual outlets must be distinct terminals"));
        }
        if ![density, sound_speed, dt, band_hz].iter().all(|v| v.is_finite() && *v > 0.0)
            || centers.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("outlet_mutual needs a finite positive gas, clock and band and finite centers"));
        }
        let nyquist = 0.5 / dt;
        if band_hz >= nyquist {
            return Err(invalid("outlet_mutual band must lie below Nyquist"));
        }
        if loads.iter().any(|l| l.sampling() != RelaxationSampling::BilinearMidpoint) {
            return Err(invalid("outlet_mutual admits bilinear-midpoint outlet loads only"));
        }
        let mut pairs = Vec::new();
        pairs.try_reserve_exact(n * (n - 1) / 2).map_err(|_| invalid("outlet_mutual allocation failed"))?;
        let mut longest = 0usize;
        for i in 0..n {
            for j in i + 1..n {
                let d = fs_math::det::hypot(centers[i][0] - centers[j][0], centers[i][1] - centers[j][1]);
                if !d.is_finite() || d <= 0.0 {
                    return Err(invalid("outlet_mutual outlets must have distinct centers"));
                }
                let rest = d / (sound_speed * dt) - 1.0;
                if !rest.is_finite() || rest < 1.0 {
                    return Err(invalid("outlet_mutual pair is closer than two samples of travel: the kernel would not be explicit"));
                }
                if rest >= MAX_DELAY {
                    return Err(invalid("outlet_mutual pair delay exceeds the retained history"));
                }
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)] // 1 <= rest < 65532
                let whole = rest.floor() as usize;
                #[allow(clippy::cast_precision_loss)]
                let fraction = rest - whole as f64;
                longest = longest.max(whole);
                pairs.push(OutletMutualPair { outlets: [i, j], distance_m: d, delay_samples: whole,
                    delay_fraction: fraction,
                    gain: density / (4.0 * core::f64::consts::PI * d) / (2.0 * dt) });
            }
        }
        let len = longest + 3;
        // Grid: 5 Hz steps strictly below Nyquist.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let points = ((nyquist / OUTLET_MUTUAL_GRID_HZ).ceil() as usize).saturating_sub(1);
        let mut matrix = vec![0.0; n * n];
        let mut factor = vec![0.0; n * n];
        let (mut worst_err, mut worst_err_hz) = (0.0f64, 0.0);
        let (mut worst_pivot, mut worst_pivot_hz) = (f64::INFINITY, 0.0);
        for k in 1..=points {
            if k % 4096 == 0 && gate.is_requested() {
                return Err(invalid("outlet_mutual admission cancelled"));
            }
            #[allow(clippy::cast_precision_loss)]
            let hz = OUTLET_MUTUAL_GRID_HZ * k as f64;
            let omega = core::f64::consts::TAU * hz;
            let warped = 2.0 * fs_math::det::tan(0.5 * omega * dt) / dt;
            for (i, load) in loads.iter().enumerate() {
                matrix[i * n + i] = self_resistance(load, warped);
            }
            for pair in &pairs {
                let m = pair.realized(omega, dt);
                let [i, j] = pair.outlets;
                matrix[i * n + j] = m.re;
                matrix[j * n + i] = m.re;
                if hz <= band_hz {
                    let exact = pair.exact(omega, density, sound_speed);
                    let e = (m - exact).abs() / exact.abs();
                    if !e.is_finite() || e > worst_err { (worst_err, worst_err_hz) = (e, hz); }
                }
            }
            // Cholesky of the symmetric Re Z; every pivot must be positive.
            for j in 0..n {
                let mut pivot = matrix[j * n + j];
                for q in 0..j { pivot -= factor[j * n + q] * factor[j * n + q]; }
                let ratio = pivot / matrix[j * n + j];
                if !(pivot > 0.0) || !ratio.is_finite() {
                    return Err(invalid("outlet_mutual Re Z has no Cholesky factor at some frequency: not passive with these self loads"));
                }
                if ratio < worst_pivot { (worst_pivot, worst_pivot_hz) = (ratio, hz); }
                let diagonal = fs_math::det::sqrt(pivot);
                factor[j * n + j] = diagonal;
                for i in j + 1..n {
                    let mut v = matrix[i * n + j];
                    for q in 0..j { v -= factor[i * n + q] * factor[j * n + q]; }
                    factor[i * n + j] = v / diagonal;
                }
            }
        }
        if !(worst_err <= OUTLET_MUTUAL_BAND_ERROR) {
            return Err(invalid("outlet_mutual realized term is more than 5 % from the exact free-field term in band"));
        }
        let mut rings = Vec::new();
        rings.try_reserve_exact(n * len).map_err(|_| invalid("outlet_mutual allocation failed"))?;
        rings.resize(n * len, 0.0);
        let min_delay = pairs.iter().map(|p| p.delay_samples).min().unwrap_or(0);
        Ok(Self { nodes: nodes.to_vec(), pairs, len, head: 0, rings, scratch: vec![0.0; n],
            record: OutletMutualRecord { density_kg_m3: density, sound_speed_m_s: sound_speed, time_step_s: dt,
                band_hz, grid_points: points, min_delay_samples: min_delay, max_band_error: worst_err,
                max_band_error_hz: worst_err_hz, min_pivot_ratio: worst_pivot, min_pivot_hz: worst_pivot_hz } })
    }

    /// Terminal nodes, in the declared outlet order.
    #[must_use]
    pub fn nodes(&self) -> &[usize] { &self.nodes }
    /// Coupled pairs with their realized delays.
    #[must_use]
    pub fn pairs(&self) -> &[OutletMutualPair] { &self.pairs }
    /// Cold admission numbers.
    #[must_use]
    pub const fn record(&self) -> &OutletMutualRecord { &self.record }

    /// Flow q[now - lag] of outlet `outlet`, lag in 1..=len; zero before the start.
    fn past(&self, outlet: usize, lag: usize) -> f64 {
        self.rings[outlet * self.len + (self.head + self.len - lag) % self.len]
    }

    /// Exterior pressure of every outlet for the NEXT step [Pa], from the
    /// accepted flows only. `out` holds one value per outlet.
    pub fn exterior_pressures(&self, out: &mut [f64]) {
        out.fill(0.0);
        for pair in &self.pairs {
            let (n, f) = (pair.delay_samples, pair.delay_fraction);
            let [a, b] = pair.outlets;
            let slope = |k: usize| (1.0 - f) * (self.past(k, n) - self.past(k, n + 2))
                + f * (self.past(k, n + 1) - self.past(k, n + 3));
            out[a] += pair.gain * slope(b);
            out[b] += pair.gain * slope(a);
        }
    }

    /// Retain one accepted step's outlet flows [m^3/s], in outlet order.
    pub fn push_flows(&mut self, flows: &[f64]) {
        for (k, &q) in flows.iter().enumerate().take(self.nodes.len()) {
            self.rings[k * self.len + self.head] = q;
        }
        self.head = (self.head + 1) % self.len;
    }

    /// Set every outlet's exterior pressure on the network for its next step.
    ///
    /// # Errors
    /// A nonfinite pressure or a terminal the network did not admit.
    pub(crate) fn apply(&mut self, physics: &mut WaveguideNetwork) -> Result<(), WaveguideError> {
        let mut out = core::mem::take(&mut self.scratch);
        self.exterior_pressures(&mut out);
        let result = self.nodes.iter().zip(&out)
            .try_for_each(|(&node, &p)| physics.set_terminal_exterior_pressure(node, p));
        self.scratch = out;
        result
    }

    /// Retain the flows the network just accepted at the outlets.
    pub(crate) fn record_step(&mut self, physics: &WaveguideNetwork) {
        for (k, &node) in self.nodes.iter().enumerate() {
            let q = physics.node_frame(node).map_or(0.0, |f| f.net_flow_into_node_m3_s);
            self.rings[k * self.len + self.head] = q;
        }
        self.head = (self.head + 1) % self.len;
    }
}

#[cfg(test)]
mod tests;
