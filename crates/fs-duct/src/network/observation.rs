//! Finite-window observations of already solved pressure/volume-flow ports.
//!
//! Calling preparation declares that all inputs describe the SAME physical
//! ports, orientation and source clock. Matching labels checks array layout;
//! NetworkCoherentResponse retains no original topology/material owner and
//! this module cannot certify that declaration. It neither solves nor changes
//! a duct, source law or BEM load. Every input model limitation still applies.

use super::{MAX_NETWORK_PORTS, NetworkCoherentResponse};
use crate::DuctError;
use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow, ProperHarmonicProduct, WindowError};
use fs_math::{c64::C64, dd::Dd};

/// Bound the uninterruptible in-place grouping sort.
pub const MAX_OBSERVATION_LINES: usize = 4096;
/// Empty factors also count against this bound.
pub const MAX_OBSERVATION_FACTORS: usize = 256;

/// Caller-authored physical correspondence, not a verified solver identity.
#[derive(Clone, Copy, Debug)]
pub struct NetworkObservationScope<'a> {
    /// Common physical ports, with every flow directed into the same network.
    /// Each input must have exactly this label order. Labels alone do not prove
    /// common topology, material, load, location or orientation.
    pub ports: &'a [String],
    /// Nonempty identifier for the authored phase origin, in seconds. All
    /// amplitudes are peak phasors Re(A exp(-i omega t)) on this clock.
    /// No clock alignment or phase origin is inferred from this identifier.
    pub source_clock: &'a str,
}

/// Admission bounds covering preparation and one observation.
#[derive(Clone, Copy, Debug)]
pub struct NetworkObservationLimits {
    pub maximum_factors: usize,
    /// Total supplied lines, including the optional deterministic mean.
    pub maximum_lines: usize,
    /// Port count times total supplied lines.
    pub maximum_port_line_pairs: usize,
    /// Port count times the sum of SQUARED supplied line counts per bank.
    /// This conservative bound is checked before allocation; equal-frequency
    /// coalescing reduces the actual observation work.
    pub maximum_port_frequency_pairs: usize,
    /// Logical heap payload: prepared fields, line/factor headers, owned labels
    /// and clock, grouping indices, query accumulators and one returned result.
    /// Borrowed inputs, stack temporaries, previous results and allocator
    /// overhead are excluded. Capacities are admitted before allocation.
    pub maximum_workspace_bytes: usize,
}

/// Averages of actual real quantities at one named physical port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkPortObservation {
    pub mean_pressure_pa: f64,
    pub mean_inflow_m3_s: f64,
    /// Mean square, without subtracting the deterministic window mean.
    pub mean_square_pressure_pa2: f64,
    pub mean_square_inflow_m6_s2: f64,
    /// <p(t) U(t)> [W], positive into the network. This can be negative.
    /// Finite-window boundary work includes storage exchange; it is not heat.
    pub mean_boundary_power_w: f64,
}

/// One finite observation, borrowing its port/clock declaration from the owner.
#[derive(Debug)]
pub struct WindowedNetworkObservation<'a> {
    window: HarmonicWindow,
    ports: &'a [String],
    source_clock: &'a str,
    values: Vec<NetworkPortObservation>,
    boundary_power_w: f64,
}

impl WindowedNetworkObservation<'_> {
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }
    pub fn ports(&self) -> &[String] {
        self.ports
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn port_values(&self) -> &[NetworkPortObservation] {
        &self.values
    }
    /// Sum of signed powers into the declared ports. No loss attribution.
    pub fn boundary_power_w(&self) -> f64 {
        self.boundary_power_w
    }
}

#[derive(Debug)]
struct Line {
    omega: f64,
    // Pressure real/imaginary, volume inflow real/imaginary. Retain the low
    // components after coalescing, including through every later product.
    fields: Vec<[Dd; 4]>,
}

#[derive(Debug)]
struct Bank {
    lines: Vec<Line>,
}

/// Optional deterministic mean plus authored proper complex fluctuation factors.
///
/// The physical field is Re(sum_j mean_j exp(-i omega_j t)
/// + sum_r xi_r sum_j factor[r,j] exp(-i omega[r,j] t)), with JOINT moments
/// E[xi_r]=0, E[xi_r conj(xi_s)]=delta_rs and E[xi_r xi_s]=0.
/// These second moments need no Gaussian distribution. Square-root statistical
/// weights are already in the supplied peak amplitudes; these are not PSD/Hz.
///
/// Frequencies within a factor interfere; different factors have no cross
/// terms, including with the deterministic mean. Properness removes fluctuation
/// sum-frequency terms even for a subperiod shutter. The deterministic mean
/// retains both sum and difference terms. Cross-frequency factor coherence
/// can make the ensemble nonstationary; no stationarity is inferred.
///
/// Preparation copies only omega and physical p/U fields from the solved
/// responses. Their cycle-averaged power summaries are not finite-window data
/// and are deliberately not used. Equal frequencies coalesce in source order
/// using double-double sums; merely close frequencies remain distinct.
#[derive(Debug)]
pub struct PreparedNetworkObservation {
    ports: Vec<String>,
    source_clock: String,
    mean: Bank,
    factors: Vec<Bank>,
    workspace_bytes: usize,
    port_frequency_pairs: usize,
}

impl PreparedNetworkObservation {
    pub fn new(
        scope: NetworkObservationScope<'_>,
        deterministic_lines: &[NetworkCoherentResponse],
        proper_factors: &[&[NetworkCoherentResponse]],
        limits: NetworkObservationLimits,
    ) -> Result<Self, DuctError> {
        Self::new_cancellable(scope, deterministic_lines, proper_factors, limits, || false)
    }

    /// As new, with bounded cancellation checks and no partially prepared bank.
    /// Calling this declares the physical correspondence described by scope;
    /// no original network identity is established by matching labels.
    pub fn new_cancellable(
        scope: NetworkObservationScope<'_>,
        deterministic_lines: &[NetworkCoherentResponse],
        proper_factors: &[&[NetworkCoherentResponse]],
        limits: NetworkObservationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(&mut cancelled, 0)?;
        let n = scope.ports.len();
        if n == 0 || n > MAX_NETWORK_PORTS || scope.source_clock.trim().is_empty() {
            return Err(bad("observation needs physical ports and an authored source clock"));
        }
        for (i, label) in scope.ports.iter().enumerate() {
            if label.trim().is_empty() || scope.ports[..i].contains(label) {
                return Err(bad("observation port labels must be nonempty and unique"));
            }
        }
        if proper_factors.len() > MAX_OBSERVATION_FACTORS
            || proper_factors.len() > limits.maximum_factors
        {
            return Err(bad("observation factor count exceeds its budget"));
        }
        let mut total_lines = 0usize;
        let mut pairs = 0usize;
        let mut largest_bank = 0usize;
        for bank in core::iter::once(deterministic_lines).chain(proper_factors.iter().copied()) {
            total_lines = total_lines.checked_add(bank.len()).ok_or_else(overflow)?;
            pairs = bank.len().checked_mul(bank.len())
                .and_then(|count| pairs.checked_add(count)).ok_or_else(overflow)?;
            largest_bank = largest_bank.max(bank.len());
        }
        let port_lines = n.checked_mul(total_lines).ok_or_else(overflow)?;
        let port_pairs = n.checked_mul(pairs).ok_or_else(overflow)?;
        if total_lines > MAX_OBSERVATION_LINES || total_lines > limits.maximum_lines
            || port_lines > limits.maximum_port_line_pairs
            || port_pairs > limits.maximum_port_frequency_pairs
        {
            return Err(bad("observation line/pair work exceeds its budget"));
        }
        let mut identity_bytes = n.checked_mul(core::mem::size_of::<String>())
            .and_then(|bytes| bytes.checked_add(scope.source_clock.len())).ok_or_else(overflow)?;
        for label in scope.ports {
            identity_bytes = identity_bytes.checked_add(label.len()).ok_or_else(overflow)?;
        }
        let per_line = n.checked_mul(core::mem::size_of::<[Dd; 4]>())
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<Line>()))
            .ok_or_else(overflow)?;
        let workspace_bytes = total_lines.checked_mul(per_line)
            .and_then(|bytes| bytes.checked_add(identity_bytes))
            .and_then(|bytes| bytes.checked_add(
                proper_factors.len().checked_mul(core::mem::size_of::<Bank>())?))
            .and_then(|bytes| bytes.checked_add(
                largest_bank.checked_mul(core::mem::size_of::<usize>())?))
            .and_then(|bytes| bytes.checked_add(n.checked_mul(
                core::mem::size_of::<[Dd; 5]>() + core::mem::size_of::<NetworkPortObservation>())?))
            .ok_or_else(overflow)?;
        if workspace_bytes > limits.maximum_workspace_bytes {
            return Err(bad("observation workspace exceeds its budget"));
        }

        let mut completed = 0usize;
        for bank in core::iter::once(deterministic_lines).chain(proper_factors.iter().copied()) {
            for line in bank {
                poll(&mut cancelled, completed)?;
                if line.ports.as_slice() != scope.ports || line.port_pressures.len() != n
                    || line.port_inflows.len() != n
                    || !line.omega.is_finite() || line.omega <= 0.0
                    || line.port_pressures.iter().chain(&line.port_inflows)
                        .any(|z| !z.re.is_finite() || !z.im.is_finite())
                {
                    return Err(bad("observation line has different layout or invalid p/U/frequency"));
                }
                completed += 1;
            }
        }
        let mut ports = storage(n)?;
        for label in scope.ports {
            ports.push(owned_text(label)?);
        }
        let source_clock = owned_text(scope.source_clock)?;
        let mean = prepare_bank(deterministic_lines, n, &mut cancelled)?;
        let mut factors = storage(proper_factors.len())?;
        for lines in proper_factors {
            factors.push(prepare_bank(lines, n, &mut cancelled)?);
        }
        let distinct_pairs = core::iter::once(&mean).chain(factors.iter())
            .try_fold(0usize, |sum, bank| {
                bank.lines.len().checked_mul(bank.lines.len())
                    .and_then(|count| sum.checked_add(count))
            }).ok_or_else(overflow)?;
        let port_frequency_pairs = n.checked_mul(distinct_pairs).ok_or_else(overflow)?;
        poll(&mut cancelled, completed)?;
        Ok(Self { ports, source_clock, mean, factors, workspace_bytes, port_frequency_pairs })
    }

    pub fn ports(&self) -> &[String] {
        &self.ports
    }
    pub fn source_clock(&self) -> &str {
        &self.source_clock
    }
    pub fn factor_count(&self) -> usize {
        self.factors.len()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Actual ordered point-pair count after coalescing within each bank.
    pub fn port_frequency_pairs(&self) -> usize {
        self.port_frequency_pairs
    }

    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedNetworkObservation<'_>, DuctError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Shared kernels are prepared once per frequency pair and reused at all
    /// ports. No network solve, time sampling or carrier-cycle average occurs.
    /// Cancelled/nonfinite/negative-self-moment results publish nothing and leave
    /// the prepared owner reusable. No relative error bound at arbitrary nulls
    /// or interval enclosure is claimed.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedNetworkObservation<'_>, DuctError> {
        poll(&mut cancelled, 0)?;
        let n = self.ports.len();
        let mut sums = storage(n)?;
        sums.resize(n, [Dd::ZERO; 5]);
        let mut completed = 0usize;
        for line in &self.mean.lines {
            poll(&mut cancelled, completed)?;
            let phase = window.mean_phase(line.omega).map_err(window_error)?;
            for (sum, field) in sums.iter_mut().zip(&line.fields) {
                sum[0] = sum[0] + real_mean([field[0], field[1]], phase);
                sum[1] = sum[1] + real_mean([field[2], field[3]], phase);
            }
            completed += 1;
        }
        for (bank_index, bank) in core::iter::once(&self.mean).chain(self.factors.iter()).enumerate() {
            poll(&mut cancelled, completed)?;
            for (i, a) in bank.lines.iter().enumerate() {
                for (j, b) in bank.lines[..=i].iter().enumerate() {
                    poll(&mut cancelled, completed)?;
                    let kernel = if bank_index == 0 {
                        Kernel::Deterministic(window.product_kernel(a.omega, b.omega)
                            .map_err(window_error)?)
                    } else {
                        Kernel::Proper(window.proper_product_kernel(a.omega, b.omega)
                            .map_err(window_error)?)
                    };
                    for ((sum, af), bf) in sums.iter_mut().zip(&a.fields).zip(&b.fields) {
                        let (ap, au) = ([af[0], af[1]], [af[2], af[3]]);
                        let (bp, bu) = ([bf[0], bf[1]], [bf[2], bf[3]]);
                        let pp = kernel.product(ap, bp)?;
                        let uu = kernel.product(au, bu)?;
                        sum[2] = sum[2] + pp;
                        sum[3] = sum[3] + uu;
                        sum[4] = sum[4] + kernel.product(ap, bu)?;
                        if i != j {
                            sum[2] = sum[2] + pp;
                            sum[3] = sum[3] + uu;
                            // Reversed frequency pair, with the original
                            // kernel orientation retained by swapping fields.
                            sum[4] = sum[4] + kernel.product(au, bp)?;
                        }
                    }
                    completed += 1;
                }
            }
        }
        let mut values = storage(n)?;
        let mut boundary_power = Dd::ZERO;
        for sum in sums {
            if sum.iter().any(|value| !finite(*value)) {
                return Err(bad("finite-window port observation is unrepresentable"));
            }
            let [p, u, pp, uu, power] = sum.map(Dd::to_f64);
            if pp < 0.0 || uu < 0.0 {
                return Err(bad("finite-window port mean square lost numerical positivity"));
            }
            boundary_power = boundary_power + sum[4];
            values.push(NetworkPortObservation {
                mean_pressure_pa: p,
                mean_inflow_m3_s: u,
                mean_square_pressure_pa2: pp,
                mean_square_inflow_m6_s2: uu,
                mean_boundary_power_w: power,
            });
        }
        if !finite(boundary_power) {
            return Err(bad("finite-window total boundary power is unrepresentable"));
        }
        poll(&mut cancelled, completed)?;
        Ok(WindowedNetworkObservation {
            window, ports: &self.ports, source_clock: &self.source_clock,
            values, boundary_power_w: boundary_power.to_f64(),
        })
    }
}

fn prepare_bank(
    input: &[NetworkCoherentResponse],
    ports: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Bank, DuctError> {
    let mut indices = storage(input.len())?;
    indices.extend(0..input.len());
    indices.sort_unstable_by(|a, b| input[*a].omega.total_cmp(&input[*b].omega)
        .then_with(|| a.cmp(b)));
    let mut lines: Vec<Line> = storage(input.len())?;
    for (completed, index) in indices.into_iter().enumerate() {
        poll(cancelled, completed)?;
        let source = &input[index];
        if lines.last().is_none_or(|line| line.omega.to_bits() != source.omega.to_bits()) {
            let mut fields = storage(ports)?;
            fields.resize(ports, [Dd::ZERO; 4]);
            lines.push(Line { omega: source.omega, fields });
        }
        let fields = &mut lines.last_mut().expect("group was just inserted").fields;
        for (port, field) in fields.iter_mut().enumerate() {
            let p = source.port_pressures[port];
            let u = source.port_inflows[port];
            for (sum, value) in field.iter_mut().zip([p.re, p.im, u.re, u.im]) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
    }
    if lines.iter().flat_map(|line| &line.fields).flatten().any(|value| !finite(*value)) {
        return Err(bad("coalesced observation field is unrepresentable"));
    }
    Ok(Bank { lines })
}

#[derive(Clone, Copy)]
enum Kernel {
    Deterministic(HarmonicProduct),
    Proper(ProperHarmonicProduct),
}

impl Kernel {
    fn leaf(self, a: C64, b: C64) -> Result<Dd, DuctError> {
        match self {
            Self::Deterministic(kernel) => kernel.mean_real_product_dd(a, b),
            Self::Proper(kernel) => kernel.mean_real_product_dd(a, b),
        }.map_err(window_error)
    }

    // Preserve grouped low components through shared C64-only kernel APIs.
    // Ordinary single-line amplitudes take just the high/high contraction.
    fn product(self, a: [Dd; 2], b: [Dd; 2]) -> Result<Dd, DuctError> {
        let ah = C64::new(a[0].hi, a[1].hi);
        let al = C64::new(a[0].lo, a[1].lo);
        let bh = C64::new(b[0].hi, b[1].hi);
        let bl = C64::new(b[0].lo, b[1].lo);
        let mut value = self.leaf(ah, bh)?;
        if al != C64::ZERO {
            value = value + self.leaf(al, bh)?;
        }
        if bl != C64::ZERO {
            value = value + self.leaf(ah, bl)?;
            if al != C64::ZERO {
                value = value + self.leaf(al, bl)?;
            }
        }
        if finite(value) { Ok(value) } else { Err(bad("port product is unrepresentable")) }
    }
}

fn real_mean(a: [Dd; 2], phase: C64) -> Dd {
    a[0] * Dd::from_f64(phase.re) - a[1] * Dd::from_f64(phase.im)
}

fn storage<T>(count: usize) -> Result<Vec<T>, DuctError> {
    let mut output = Vec::new();
    output.try_reserve_exact(count).map_err(|_| bad("observation allocation failed"))?;
    Ok(output)
}

fn owned_text(text: &str) -> Result<String, DuctError> {
    let mut output = String::new();
    output.try_reserve_exact(text.len()).map_err(|_| bad("observation label allocation failed"))?;
    output.push_str(text);
    Ok(output)
}

fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

fn overflow() -> DuctError {
    bad("observation work/storage arithmetic overflow")
}

fn window_error(_: WindowError) -> DuctError {
    bad("observation temporal moment is unrepresentable")
}

fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() { Err(DuctError::Cancelled { completed_steps }) } else { Ok(()) }
}
