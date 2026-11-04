//! Spatial finite-window moments of already propagated acoustic ensembles.
//!
//! Every receiver shares one authored source clock and ordered fluctuation
//! coordinates. A coordinate names the SAME random coefficient at all receivers;
//! equal factor counts or laws alone never establish physical correlation.
//! The binding is supplied by the caller, not inferred from distance or phase.
//! Original pressure/velocity fields carry all propagation and source limits.

use super::*;
use fs_material::fluid_acoustics::FluidAcousticState;

mod readout;
pub use readout::{
    HarmonicReadoutBudget, HarmonicReadoutTap, PreparedHarmonicReadout, WindowedHarmonicReadoutPair,
};

/// One explicitly shared zero-mean coefficient in the physical source model.
/// All receiver ensembles bind their factor at this index to this identity.
/// Distinct coordinates retain the original ensemble's uncorrelated-factor
/// assumptions, including vanishing real/proper mixed moments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HarmonicFactorCoordinate<'a> {
    /// Caller-authored physical coefficient identity, unique within the array.
    pub id: &'a str,
    /// Unit-variance real or unit-mean-square proper complex coefficient law.
    pub statistics: HarmonicFactorStatistics,
}

/// One physical receiver's original field and local microphone calibration.
/// Positions and velocity components use a single common Cartesian frame.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicArrayReceiver<'a, const D: usize> {
    /// Position [m] in the common frame; no propagation is inferred from it.
    pub position_m: [f64; D],
    /// Original local state supplying rho*c for directional microphone response.
    pub acoustic: FluidAcousticState,
    /// Already propagated fields with factors in the declared coordinate order.
    pub ensemble: &'a PreparedHarmonicEnsemble<D>,
}

/// Bounds both retained receiver metadata and all cross-receiver queries.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicArrayBudget {
    /// Maximum retained receiver count.
    pub max_receivers: usize,
    /// Maximum common fluctuation-coordinate count, including silent factors.
    pub max_factors: usize,
    /// Total already-grouped lines over every receiver, mean and factor.
    pub max_total_lines: usize,
    /// Sum, over mean/factor banks, of the squared line count across receivers.
    /// Bounds all spatial frequency pairs; each has at most sixteen p/u products.
    pub max_line_pairs: usize,
}

/// Real cross moments between two locations, in left/right receiver order.
/// Cross-location p/u products have intensity units but are NOT local flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticCrossMoments<const D: usize> {
    /// Pressure at the left receiver times pressure at the right [Pa²].
    pub pressure_pressure_pa2: f64,
    /// Left pressure times each right velocity component [W/m²].
    pub pressure_velocity_w_m2: [f64; D],
    /// Each left velocity component times right pressure [W/m²].
    pub velocity_pressure_w_m2: [f64; D],
    /// Left velocity i times right velocity j [m²/s²].
    pub velocity_velocity_m2_s2: [[f64; D]; D],
}

/// Two different spatial statistics on the same physical detector interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowedHarmonicReceiverPair<const D: usize> {
    /// Common physical observer interval [s].
    pub window: HarmonicWindow,
    /// Ordered left/right receiver indices.
    pub receivers: [usize; 2],
    /// <E[x_left(t) x_right(t)]>, including the deterministic mean waveform.
    pub second_moments: AcousticCrossMoments<D>,
    /// Cov(<x_left>, <x_right>) over the authored source ensemble. Only
    /// fluctuation factors contribute; no large mean-square subtraction occurs.
    /// This differs from the time average of instantaneous covariance.
    pub averaged_covariance: AcousticCrossMoments<D>,
}

/// Borrowed original receiver ensembles with an explicit common factor binding.
/// Preparation copies only receiver metadata; it does not rebuild propagation
/// or duplicate spectral coefficients. Pair queries allocate nothing and use
/// the original DD harmonic-window kernels for deterministic, real and proper
/// factors. Close frequencies retain their beats and real factors retain sum
/// frequencies. Covariance does not assume Gaussianity or stationarity.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicArray<'a, const D: usize> {
    source_clock: &'a str,
    coordinates: &'a [HarmonicFactorCoordinate<'a>],
    receivers: Vec<HarmonicArrayReceiver<'a, D>>,
}

impl<'a, const D: usize> PreparedHarmonicArray<'a, D> {
    /// Bind each receiver's factor k to `coordinates[k]` on `source_clock`.
    /// The caller establishes actual shared coefficient identity and clock;
    /// this constructor checks the declared shape/law, never estimates it.
    /// Dimension is explicitly 2 or 3. A source can be silent at a receiver
    /// through an empty factor, which must retain its original coordinate slot.
    pub fn new(
        source_clock: &'a str,
        coordinates: &'a [HarmonicFactorCoordinate<'a>],
        receivers: &[HarmonicArrayReceiver<'a, D>],
        budget: HarmonicArrayBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(source_clock, coordinates, receivers, budget, || false)
    }

    /// Admit the common binding with bounded work and cancellation checks.
    pub fn new_cancellable(
        source_clock: &'a str,
        coordinates: &'a [HarmonicFactorCoordinate<'a>],
        receivers: &[HarmonicArrayReceiver<'a, D>],
        budget: HarmonicArrayBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if !(D == 2 || D == 3)
            || source_clock.trim().is_empty()
            || receivers.is_empty()
            || receivers.len() > budget.max_receivers
            || coordinates.len() > budget.max_factors
        {
            return Err(invalid(
                "acoustic array needs a clock, bounded receivers and 2D/3D frame",
            ));
        }
        let mut ids = Vec::new();
        ids.try_reserve_exact(coordinates.len())
            .map_err(|_| invalid("acoustic factor identity allocation failed"))?;
        for coordinate in coordinates {
            poll(&mut cancelled)?;
            if coordinate.id.trim().is_empty() {
                return Err(invalid("acoustic array factor identity is empty"));
            }
            ids.push(coordinate.id);
        }
        ids.sort_unstable();
        if ids.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid("acoustic array factor identities must be unique"));
        }
        for receiver in receivers {
            poll(&mut cancelled)?;
            if receiver.position_m.iter().any(|x| !x.is_finite())
                || receiver.ensemble.factors.len() != coordinates.len()
                || receiver
                    .ensemble
                    .factors
                    .iter()
                    .zip(coordinates)
                    .any(|((law, _), coordinate)| *law != coordinate.statistics)
            {
                return Err(invalid(
                    "acoustic receiver geometry or common factor binding differs",
                ));
            }
        }
        let mut total = 0usize;
        let mut pairs = 0usize;
        for bank in 0..=coordinates.len() {
            let mut lines = 0usize;
            for receiver in receivers {
                poll(&mut cancelled)?;
                lines = lines
                    .checked_add(field(receiver.ensemble, bank).lines.len())
                    .ok_or_else(|| invalid("acoustic array line count overflows"))?;
            }
            total = total
                .checked_add(lines)
                .ok_or_else(|| invalid("acoustic array line count overflows"))?;
            pairs = lines
                .checked_mul(lines)
                .and_then(|n| pairs.checked_add(n))
                .ok_or_else(|| invalid("acoustic array pair work overflows"))?;
        }
        if total > budget.max_total_lines || pairs > budget.max_line_pairs {
            return Err(invalid("acoustic array exceeds line/pair work budget"));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(receivers.len())
            .map_err(|_| invalid("acoustic receiver allocation failed"))?;
        retained.extend_from_slice(receivers);
        poll(&mut cancelled)?;
        Ok(Self {
            source_clock,
            coordinates,
            receivers: retained,
        })
    }

    /// Authored common physical time origin identity.
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    /// Ordered shared physical coefficients, preserving silent coordinates.
    pub fn factor_coordinates(&self) -> &[HarmonicFactorCoordinate<'a>] {
        self.coordinates
    }
    /// Original receiver fields and metadata in admitted order.
    pub fn receivers(&self) -> &[HarmonicArrayReceiver<'a, D>] {
        &self.receivers
    }

    /// Original local moments, including signed intensity, without recomputing
    /// a propagator or interpreting cross-location products as local power.
    pub fn observe_receiver(
        &self,
        receiver: usize,
        window: HarmonicWindow,
        cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticEnsemble<D>, AeroacError> {
        self.receiver(receiver)?
            .ensemble
            .observe_window_cancellable(window, cancelled)
    }

    /// Cross-receiver second moments and covariance of averaged signals.
    pub fn observe_pair(
        &self,
        receivers: [usize; 2],
        window: HarmonicWindow,
    ) -> Result<WindowedHarmonicReceiverPair<D>, AeroacError> {
        self.observe_pair_cancellable(receivers, window, || false)
    }

    /// Allocation-free cross-location products, with a poll per frequency pair.
    /// Reversing the receiver order transposes p/u order and the velocity tensor.
    pub fn observe_pair_cancellable(
        &self,
        receivers: [usize; 2],
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedHarmonicReceiverPair<D>, AeroacError> {
        poll(&mut cancelled)?;
        let left = self.receiver(receivers[0])?.ensemble;
        let right = self.receiver(receivers[1])?.ensemble;
        let (second_moments, averaged_covariance) = pair(
            left,
            right,
            window,
            receivers[0] == receivers[1],
            &mut cancelled,
        )?;
        Ok(WindowedHarmonicReceiverPair {
            window,
            receivers,
            second_moments,
            averaged_covariance,
        })
    }

    fn receiver(&self, index: usize) -> Result<&HarmonicArrayReceiver<'a, D>, AeroacError> {
        self.receivers
            .get(index)
            .ok_or_else(|| invalid("acoustic array receiver index"))
    }
}

fn field<const D: usize>(
    ensemble: &PreparedHarmonicEnsemble<D>,
    bank: usize,
) -> &PreparedHarmonicField<D> {
    if bank == 0 {
        &ensemble.mean
    } else {
        &ensemble.factors[bank - 1].1
    }
}

fn kernel(
    window: HarmonicWindow,
    kind: FactorKind,
    a: f64,
    b: f64,
) -> Result<MomentKernel, AeroacError> {
    Ok(match kind {
        FactorKind::Deterministic => {
            MomentKernel::Deterministic(window.product_kernel_hz(a, b).map_err(window_error)?)
        }
        FactorKind::Proper => MomentKernel::Proper(
            window
                .proper_product_kernel_hz(a, b)
                .map_err(window_error)?,
        ),
    })
}

fn pair<const D: usize>(
    left: &PreparedHarmonicEnsemble<D>,
    right: &PreparedHarmonicEnsemble<D>,
    window: HarmonicWindow,
    self_pair: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(AcousticCrossMoments<D>, AcousticCrossMoments<D>), AeroacError> {
    let mut moments = CrossSums::new();
    let mut covariance = CrossSums::new();
    for bank in 0..=left.factors.len() {
        poll(cancelled)?;
        let kind = if bank == 0 || left.factors[bank - 1].0 == HarmonicFactorStatistics::Real {
            FactorKind::Deterministic
        } else {
            FactorKind::Proper
        };
        let a = field(left, bank);
        let b = field(right, bank);
        if a.lines.is_empty() || b.lines.is_empty() {
            continue;
        }
        for x in &a.lines {
            for y in &b.lines {
                poll(cancelled)?;
                moments.add(kernel(window, kind, x.frequency_hz, y.frequency_hz)?, x, y)?;
            }
        }
        if bank != 0 {
            // Covariance of time-averaged signals: integrate each factor FIRST.
            // This is not the average instantaneous covariance above.
            let a = averaged_amplitudes(a, window, cancelled)?;
            let b = averaged_amplitudes(b, window, cancelled)?;
            covariance.add(kernel(window, kind, 0.0, 0.0)?, &a, &b)?;
        }
    }
    poll(cancelled)?;
    Ok((moments.finish(self_pair)?, covariance.finish(self_pair)?))
}

fn averaged_amplitudes<const D: usize>(
    field: &PreparedHarmonicField<D>,
    window: HarmonicWindow,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedLine<D>, AeroacError> {
    let mut result = PreparedLine {
        frequency_hz: 0.0,
        pressure: ComplexSum::default(),
        velocity: [ComplexSum::default(); D],
    };
    for line in &field.lines {
        poll(cancelled)?;
        let integral = window
            .product_kernel_hz(line.frequency_hz, 0.0)
            .map_err(window_error)?;
        let add = |sum: &mut ComplexSum, value: ComplexSum| -> Result<(), AeroacError> {
            let [re, im] = value.parts();
            sum.re = sum.re
                + integral
                    .mean_real_product_parts_dd([re, im], [Dd::ONE, Dd::ZERO])
                    .map_err(window_error)?;
            sum.im = sum.im
                + integral
                    .mean_real_product_parts_dd([im, -re], [Dd::ONE, Dd::ZERO])
                    .map_err(window_error)?;
            Ok(())
        };
        add(&mut result.pressure, line.pressure)?;
        for (out, &value) in result.velocity.iter_mut().zip(&line.velocity) {
            add(out, value)?;
        }
    }
    Ok(result)
}

struct CrossSums<const D: usize> {
    pp: Dd,
    pu: [Dd; D],
    up: [Dd; D],
    uu: [[Dd; D]; D],
}
impl<const D: usize> CrossSums<D> {
    fn new() -> Self {
        Self {
            pp: Dd::ZERO,
            pu: [Dd::ZERO; D],
            up: [Dd::ZERO; D],
            uu: [[Dd::ZERO; D]; D],
        }
    }
    fn add(
        &mut self,
        kernel: MomentKernel,
        a: &PreparedLine<D>,
        b: &PreparedLine<D>,
    ) -> Result<(), AeroacError> {
        self.pp = self.pp + kernel.product(a.pressure, b.pressure)?;
        for i in 0..D {
            self.pu[i] = self.pu[i] + kernel.product(a.pressure, b.velocity[i])?;
            self.up[i] = self.up[i] + kernel.product(a.velocity[i], b.pressure)?;
            for j in 0..D {
                self.uu[i][j] = self.uu[i][j] + kernel.product(a.velocity[i], b.velocity[j])?;
            }
        }
        Ok(())
    }
    fn finish(self, self_pair: bool) -> Result<AcousticCrossMoments<D>, AeroacError> {
        let value = AcousticCrossMoments {
            pressure_pressure_pa2: self.pp.to_f64(),
            pressure_velocity_w_m2: self.pu.map(Dd::to_f64),
            velocity_pressure_w_m2: self.up.map(Dd::to_f64),
            velocity_velocity_m2_s2: self.uu.map(|row| row.map(Dd::to_f64)),
        };
        if !value.pressure_pressure_pa2.is_finite()
            || value
                .pressure_velocity_w_m2
                .iter()
                .chain(&value.velocity_pressure_w_m2)
                .chain(value.velocity_velocity_m2_s2.iter().flatten())
                .any(|x| !x.is_finite())
        {
            return Err(nonfinite("finite-window acoustic cross moments"));
        }
        if self_pair
            && (value.pressure_pressure_pa2 < 0.0
                || (0..D).any(|i| value.velocity_velocity_m2_s2[i][i] < 0.0))
        {
            return Err(invalid(
                "acoustic array self moment lost numerical positivity",
            ));
        }
        Ok(value)
    }
}
