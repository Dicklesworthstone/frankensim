//! Propagate original nonlinear port harmonics and probe sidebands through BEM.
//!
//! The source observer owns nonlinear generation and its truncation evidence;
//! BEM owns each frequency's actual propagation. This bridge only converts the
//! declared phasor convention and physical source coordinates. It does not add
//! radiation loading to the orbit, infer port units, or synthesize an envelope.

use super::{sideband::PhsPortSidebandResponse, spectrum::PortSpectrum};
use fs_aeroac::{
    AeroacError,
    observation::{AngularHarmonicFieldLine, HarmonicObservationBudget, PreparedHarmonicField},
};
use fs_bem::observation_transfers::{
    HarmonicField, PreparedHarmonicObservation, ReceiverSpec, TransferError,
};
use fs_math::{c64::C64, dd::Dd};

/// Bounds apply separately to spectrum preparation and complete propagation.
#[derive(Clone, Copy, Debug)]
pub struct PortAcousticLimits {
    /// Original physical output ports.
    pub maximum_ports: usize,
    /// Counts authored signed bands before folding, including a DC band.
    pub maximum_lines: usize,
    /// Complete physical receiver array, not a per-source allowance.
    pub maximum_receivers: usize,
    /// Conservative work score for folding, physical maps, BEM application and
    /// receiver-line grouping. No sparsity or zero-amplitude discount is used.
    pub maximum_products: usize,
    /// Peak additional numeric payload, including the returned value. Borrowed
    /// spectra, BEM banks/maps, receiver names and allocator metadata are excluded.
    /// An additional 256 MiB hard cap applies; no field solve is performed here.
    pub maximum_workspace_bytes: usize,
}

/// A failed preparation or propagation publishes no partial physical field.
#[derive(Debug)]
pub enum PortAcousticError {
    /// Incompatible original shapes, clocks, coordinates or nonfinite data.
    Invalid(&'static str),
    /// The complete operation exceeds an explicit resource allowance.
    Budget(&'static str),
    /// The caller cancelled before publication.
    Cancelled,
    /// The original propagation owner refused its physical field.
    Bem(TransferError),
    /// The existing acoustic observation owner refused the resulting lines.
    Field(AeroacError),
}
impl core::fmt::Display for PortAcousticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical port acoustic propagation: {self:?}")
    }
}
impl std::error::Error for PortAcousticError {}
type Result<T> = core::result::Result<T, PortAcousticError>;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PortAcousticError::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(message: &'static str) -> PortAcousticError {
    PortAcousticError::Invalid(message)
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn reserved<T>(count: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| PortAcousticError::Budget("physical spectrum allocation failed"))?;
    Ok(values)
}
fn admit(
    limits: PortAcousticLimits,
    ports: usize,
    lines: usize,
    receivers: usize,
    work: u128,
    bytes: u128,
) -> Result<()> {
    if ports > limits.maximum_ports
        || lines > limits.maximum_lines
        || receivers > limits.maximum_receivers
        || lines > 65_536
    {
        return Err(PortAcousticError::Budget(
            "port, line or receiver count exceeds its limit",
        ));
    }
    if work > limits.maximum_products as u128 {
        return Err(PortAcousticError::Budget(
            "physical propagation exceeds its work limit",
        ));
    }
    if bytes > limits.maximum_workspace_bytes as u128 || bytes > 256 * 1024 * 1024 {
        return Err(PortAcousticError::Budget(
            "physical propagation exceeds its workspace limit",
        ));
    }
    Ok(())
}
fn admit_spectrum(limits: PortAcousticLimits, ports: usize, lines: usize) -> Result<()> {
    // Public response metadata can contain arbitrary usize counts without a
    // backing allocation. Bound lines before forming even the u128 work score.
    admit(limits, ports, lines, 0, 0, 0)?;
    let (p, l) = (ports as u128, lines as u128);
    admit(
        limits,
        ports,
        lines,
        0,
        l * (l + 8 * p) + p,
        32 * l * p + 64 * l + 40 * p,
    )?;
    if ports == 0 {
        return Err(invalid("an original physical port is required"));
    }
    Ok(())
}

/// A real physical port waveform, retaining its original units and phase clock.
///
/// Positive lines use peak amplitudes in `Re(A exp(-i omega t))`. Only exactly
/// coincident physical frequencies combine; close lines retain their beats.
/// DC stays in original port units, separate from the positive-frequency BEM
/// field. Preparation neither solves nor validates the source orbit: its native
/// residual, truncation and small-probe assumptions remain authoritative.
pub struct PortAcousticSpectrum {
    dc: Vec<f64>,
    omega: Vec<f64>,
    // Frequency-major, then original physical port.
    amplitudes: Vec<C64>,
}

impl PortAcousticSpectrum {
    /// Consume the original nonlinear-storage observer's already negative-time
    /// peak phasors. No conjugation, doubling or amplitude threshold is applied.
    pub fn from_harmonics(
        source: &PortSpectrum,
        limits: PortAcousticLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let ports = source.dc.len();
        let lines = source.peak_phasors.first().map_or(0, Vec::len);
        admit_spectrum(limits, ports, lines)?;
        if !source.omega_rad_s.is_finite()
            || source.omega_rad_s <= 0.0
            || source.peak_phasors.len() != ports
            || source.dc.iter().any(|value| !value.is_finite())
            || source.peak_phasors.iter().any(|row| row.len() != lines)
        {
            return Err(invalid("complete finite physical port spectrum required"));
        }
        let mut omega = reserved(lines)?;
        let mut amplitudes = reserved(ports * lines)?;
        for line in 0..lines {
            poll(&mut cancelled)?;
            let frequency = (line + 1) as f64 * source.omega_rad_s;
            if !frequency.is_finite() || frequency <= 0.0 {
                return Err(invalid("physical harmonic frequency is not representable"));
            }
            omega.push(frequency);
            for (index, row) in source.peak_phasors.iter().enumerate() {
                if index % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                if !finite(row[line]) {
                    return Err(invalid("physical port amplitude must be finite"));
                }
                amplitudes.push(row[line]);
            }
        }
        let mut dc = reserved(ports)?;
        dc.extend_from_slice(&source.dc);
        poll(&mut cancelled)?;
        Ok(Self {
            dc,
            omega,
            amplitudes,
        })
    }

    /// Fold the REAL infinitesimal probe, not its complex envelope, into the
    /// existing negative-time acoustic convention. A positive signed physical
    /// frequency conjugates its peak coefficient; a negative one does not.
    /// At zero frequency only the real coefficient is physical DC. There is
    /// never a factor of two. Coincident opposite bands add coherently in DD
    /// BEFORE any source map or quadratic acoustic observation is evaluated.
    ///
    /// Signed frequencies use the native sideband owner's fused
    /// `k.mul_add(omega, carrier)`. The returned ordered frequencies can be used
    /// directly to prepare the required BEM operators, including subharmonics.
    pub fn from_sidebands(
        source: &PhsPortSidebandResponse,
        limits: PortAcousticLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let ports = source.port_count;
        let lines = source
            .output_harmonics
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| invalid("sideband count overflow"))?;
        admit_spectrum(limits, ports, lines)?;
        if !source.omega_rad_s.is_finite()
            || source.omega_rad_s <= 0.0
            || !source.carrier_angular_frequency.is_finite()
            || ports.checked_mul(lines) != Some(source.port_coefficients.len())
        {
            return Err(invalid("complete finite physical sideband output required"));
        }
        let mut bands = reserved(lines)?;
        let mut dc_sum = reserved(ports)?;
        dc_sum.resize(ports, Dd::ZERO);
        for (index, row) in source.port_coefficients.chunks_exact(ports).enumerate() {
            poll(&mut cancelled)?;
            if row.iter().any(|&value| !finite(value)) {
                return Err(invalid("physical sideband amplitude must be finite"));
            }
            let order = index as isize - source.output_harmonics as isize;
            let signed =
                (order as f64).mul_add(source.omega_rad_s, source.carrier_angular_frequency);
            if !signed.is_finite() {
                return Err(invalid("physical sideband frequency is not representable"));
            }
            if signed == 0.0 {
                for (sum, coefficient) in dc_sum.iter_mut().zip(row) {
                    *sum = *sum + Dd::from_f64(coefficient.re);
                }
            } else {
                bands.push((signed.abs(), index, signed > 0.0));
            }
        }
        bands.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        poll(&mut cancelled)?;
        let mut omega = reserved(bands.len())?;
        let mut amplitudes = reserved(ports * bands.len())?;
        let mut sums = reserved(ports)?;
        sums.resize(ports, [Dd::ZERO; 2]);
        let mut start = 0;
        while start < bands.len() {
            poll(&mut cancelled)?;
            let frequency = bands[start].0;
            sums.fill([Dd::ZERO; 2]);
            let mut end = start;
            while end < bands.len() && bands[end].0 == frequency {
                let (_, index, conjugate) = bands[end];
                let row = &source.port_coefficients[index * ports..(index + 1) * ports];
                for (port, (sum, coefficient)) in sums.iter_mut().zip(row).enumerate() {
                    if port % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    sum[0] = sum[0] + Dd::from_f64(coefficient.re);
                    sum[1] = sum[1]
                        + Dd::from_f64(if conjugate {
                            -coefficient.im
                        } else {
                            coefficient.im
                        });
                }
                end += 1;
            }
            omega.push(frequency);
            for sum in &sums {
                amplitudes.push(complex(*sum)?);
            }
            start = end;
        }
        let mut dc = reserved(ports)?;
        for value in dc_sum {
            let value = value.to_f64();
            if !value.is_finite() {
                return Err(invalid("physical DC overflow"));
            }
            dc.push(value);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            dc,
            omega,
            amplitudes,
        })
    }

    /// Number of original physical ports, including any silent coordinates.
    pub fn port_count(&self) -> usize {
        self.dc.len()
    }
    /// Actual positive physical frequencies, strictly increasing, in rad/s.
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.omega
    }
    /// Static original-port output, not a guessed steady acoustic pressure.
    pub fn dc_port_outputs(&self) -> &[f64] {
        &self.dc
    }
    /// Complete original port row for one positive frequency, in source order.
    pub fn peak_port_outputs(&self, line: usize) -> Option<&[C64]> {
        let start = line.checked_mul(self.port_count())?;
        self.amplitudes
            .get(start..start.checked_add(self.port_count())?)
    }

    /// Apply exact-frequency native BEM fields to the complete real spectrum.
    ///
    /// There must be one binding for every retained positive line, in the same
    /// order as `angular_frequencies_rad_s()`, including zero-amplitude lines.
    /// All native banks must have exactly the supplied receivers/geometry. No
    /// nearest-frequency match, transfer interpolation, plane-wave intensity
    /// shortcut or material/source-unit inference is performed. A DC-only
    /// source takes no bindings and returns silent acoustic receivers plus DC.
    ///
    /// The returned fields go directly to the existing finite-window observer
    /// or `pcm_wav::observation::harmonic::HarmonicPressureRenderer`. Their
    /// common clock is the source's original clock; BEM already includes the
    /// propagation phase. Sampling bandwidth is admitted by that audio owner.
    pub fn propagate<'a>(
        &self,
        receivers: &'a [ReceiverSpec],
        bindings: &[PortAcousticBinding<'a>],
        limits: PortAcousticLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PortAcousticObservation<'a>> {
        poll(&mut cancelled)?;
        let (ports, lines, count) = (self.port_count(), self.omega.len(), receivers.len());
        if bindings.len() != lines || count == 0 {
            return Err(invalid(
                "one exact-frequency binding per line and physical receivers required",
            ));
        }
        let (p, l, r) = (ports as u128, lines as u128, count as u128);
        let mut work = r * l * l + r + l * p;
        let mut max_sources = 0;
        for binding in bindings {
            let sources = binding.observation.sources().len();
            work += sources as u128 * (p + 4 * r);
            max_sources = max_sources.max(sources);
        }
        let bytes = 256 * r * l + 160 * r + 32 * max_sources as u128 + 144 * l + 32 * p;
        admit(limits, ports, lines, count, work, bytes)?;
        if receivers.iter().any(|receiver| {
            receiver.id.is_empty() || receiver.position_m.iter().any(|value| !value.is_finite())
        }) {
            return Err(invalid("finite identified physical receivers required"));
        }
        for (binding, &omega) in bindings.iter().zip(&self.omega) {
            poll(&mut cancelled)?;
            if binding.observation.angular_frequency_rad_s() != omega
                || binding.observation.receivers() != receivers
                || binding.observation.sources().len().checked_mul(ports)
                    != Some(binding.source_from_ports.len())
                || binding
                    .source_from_ports
                    .iter()
                    .any(|&value| !finite(value))
            {
                return Err(invalid(
                    "BEM frequency, receivers or physical source map differs from the spectrum",
                ));
            }
        }
        let mut received = reserved(count)?;
        for _ in 0..count {
            received.push(reserved(lines)?);
        }
        let mut fields = reserved(count)?;
        fields.resize(count, HarmonicField::default());
        let mut source = reserved(max_sources)?;
        let mut observations = reserved(lines)?;
        for (line, binding) in bindings.iter().enumerate() {
            poll(&mut cancelled)?;
            source.clear();
            let ports = self
                .peak_port_outputs(line)
                .expect("admitted original line");
            for row in binding.source_from_ports.chunks_exact(self.port_count()) {
                poll(&mut cancelled)?;
                let mut sum = [Dd::ZERO; 2];
                for (index, (&map, &port)) in row.iter().zip(ports).enumerate() {
                    if index % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let (a, b, c, d) = (
                        Dd::from_f64(map.re),
                        Dd::from_f64(map.im),
                        Dd::from_f64(port.re),
                        Dd::from_f64(port.im),
                    );
                    sum[0] = sum[0] + a * c - b * d;
                    sum[1] = sum[1] + a * d + b * c;
                }
                source.push(complex(sum)?);
            }
            binding
                .observation
                .evaluate_coherent(&source, &mut fields, &mut cancelled)
                .map_err(|error| match error {
                    TransferError::Cancelled => PortAcousticError::Cancelled,
                    other => PortAcousticError::Bem(other),
                })?;
            for (receiver, field) in received.iter_mut().zip(&fields) {
                receiver.push(AngularHarmonicFieldLine {
                    angular_frequency_rad_s: self.omega[line],
                    pressure_pa: field.pressure_pa,
                    particle_velocity_m_s: field.particle_velocity_m_s,
                });
            }
            observations.push(binding.observation);
        }
        let mut acoustic = reserved(count)?;
        let budget = HarmonicObservationBudget {
            max_lines: lines,
            max_line_pairs: lines * lines,
        };
        for receiver in received {
            acoustic.push(
                PreparedHarmonicField::new_angular_cancellable(&receiver, budget, &mut cancelled)
                    .map_err(|error| match error {
                    AeroacError::Cancelled { .. } => PortAcousticError::Cancelled,
                    other => PortAcousticError::Field(other),
                })?,
            );
        }
        let mut dc = reserved(self.port_count())?;
        dc.extend_from_slice(&self.dc);
        poll(&mut cancelled)?;
        Ok(PortAcousticObservation {
            receivers,
            fields: acoustic,
            dc,
            observations,
        })
    }
}

fn complex(sum: [Dd; 2]) -> Result<C64> {
    let value = C64::new(sum[0].to_f64(), sum[1].to_f64());
    if !finite(value) {
        return Err(invalid("coherent physical source overflow"));
    }
    Ok(value)
}

/// Explicit physical law at this bank's positive angular frequency.
#[derive(Clone, Copy)]
pub struct PortAcousticBinding<'a> {
    /// Already solved physical source-to-receiver propagation, borrowed intact.
    pub observation: &'a PreparedHarmonicObservation,
    /// Source-major, then original port. Units are the declared BEM source's
    /// SI unit per original port unit, in the common exp(-i omega t) convention.
    /// The caller owns this law, source orientation and radiation loading;
    /// complex entries may retain physical delay or a derivative conversion.
    pub source_from_ports: &'a [C64],
}

/// Propagated positive-frequency field with the static physical output retained.
/// Borrowed native owners preserve their source identities, receiver geometry
/// and producer uncertainty; this composition does not upgrade those claims.
pub struct PortAcousticObservation<'a> {
    receivers: &'a [ReceiverSpec],
    fields: Vec<PreparedHarmonicField<3>>,
    dc: Vec<f64>,
    observations: Vec<&'a PreparedHarmonicObservation>,
}
impl<'a> PortAcousticObservation<'a> {
    /// Common receiver identities and metre-valued coordinates, in field order.
    pub fn receivers(&self) -> &'a [ReceiverSpec] {
        self.receivers
    }
    /// Original calibrated acoustic field owner for each physical receiver.
    pub fn fields(&self) -> &[PreparedHarmonicField<3>] {
        &self.fields
    }
    /// Original-port DC, before any undefined zero-frequency acoustic map.
    pub fn dc_port_outputs(&self) -> &[f64] {
        &self.dc
    }
    /// Native positive-frequency propagation owners, in spectrum order.
    pub fn observations(&self) -> &[&'a PreparedHarmonicObservation] {
        &self.observations
    }
}
