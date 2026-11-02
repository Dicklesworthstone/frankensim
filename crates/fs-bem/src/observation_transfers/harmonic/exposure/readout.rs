//! Joint uncertainty of linearly averaged physical sensor readouts.
//!
//! Averaging a fluctuating field precedes the covariance here. This is
//! different from averaging instantaneous field products (energy/intensity).
//! A coherent proper factor is integrated across all frequencies before its
//! outer product, retaining cross-frequency and cross-receiver interference.

use super::*;

/// Fixed ceiling on requested receiver pairs, including repeats.
pub const MAX_READOUT_RECEIVER_PAIRS: usize = 65_536;

/// Ordered physical receivers; reversing a pair transposes its field covariance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HarmonicReceiverPair {
    pub first: usize,
    pub second: usize,
}

/// Per-query bounds. The prepared field's existing storage is not counted again.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReadoutLimits {
    pub maximum_receiver_pairs: usize,
    /// Eight DD real contractions per retained receiver/frequency and sixteen
    /// proper covariance contractions per requested pair/factor.
    pub maximum_products: usize,
    /// Logical heap payload of the temporary integrated factor, covariance
    /// accumulators, returned means and returned pair covariances.
    pub maximum_workspace_bytes: usize,
}

/// Ensemble mean of the normalized time-averaged pressure/velocity readout.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicAveragedField {
    pub pressure_pa: f64,
    pub particle_velocity_m_s: [f64; 3],
}

/// CENTERED covariance of two time-averaged real sensor readouts.
///
/// The two pressure/velocity blocks need not agree at different receivers.
/// Their units equal intensity units, but these are NOT energy fluxes:
/// Cov(average(p), average(u)) differs from average(p*u). Deterministic fields
/// contribute to the means only; proper fluctuation factors determine these
/// covariances without subtracting a large squared deterministic mean.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicReadoutCovariance {
    pub receivers: HarmonicReceiverPair,
    pub pressure_covariance_pa2: f64,
    pub first_pressure_second_velocity_covariance_w_m2: [f64; 3],
    pub first_velocity_second_pressure_covariance_w_m2: [f64; 3],
    pub velocity_covariance_m2_s2: [[f64; 3]; 3],
}

/// Complete readout result retaining the original physical array and clock.
pub struct WindowedHarmonicReadout<'a> {
    window: HarmonicWindow,
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    means: Vec<HarmonicAveragedField>,
    covariances: Vec<HarmonicReadoutCovariance>,
}

impl WindowedHarmonicReadout<'_> {
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        &self.reference.receivers
    }
    pub fn sources(&self) -> &[SourceId] {
        &self.reference.sources
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.reference.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        ProducerUncertainty::Unresolved {
            source_hash: self.reference.producer_hash,
        }
    }
    /// Receiver order is exactly the prepared physical array order.
    pub fn means(&self) -> &[HarmonicAveragedField] {
        &self.means
    }
    /// Requested order is retained, including reversed or repeated pairs.
    pub fn covariances(&self) -> &[HarmonicReadoutCovariance] {
        &self.covariances
    }
}

impl PreparedHarmonicExposure<'_> {
    /// Observe uncertainty of a linear pressure/velocity detector that reports
    /// the normalized window average. A long observation can cancel the
    /// readout variance even while the instantaneous acoustic energy is nonzero.
    ///
    /// The existing explicitly proper source-factor assumptions apply. No
    /// Gaussian law, detector electronics/noise, calibration or sampling
    /// bandwidth is invented. The receiver correlations support authored real
    /// array combinations: a pressure readout with weights w has mean w^T m
    /// and variance w^T C w when its required receiver pairs are requested.
    ///
    /// Work is O(receivers * retained lines + requested pairs * factors);
    /// there is no quadratic frequency-pair loop and no BEM solve. Exactly
    /// repeated frequencies were already combined by field preparation.
    pub fn observe_readout(
        &self,
        window: HarmonicWindow,
        pairs: &[HarmonicReceiverPair],
        limits: HarmonicReadoutLimits,
    ) -> Result<WindowedHarmonicReadout<'_>, TransferError> {
        self.observe_readout_cancellable(window, pairs, limits, || false)
    }

    /// Refusal or cancellation publishes no partial array and leaves the
    /// immutable prepared field reusable.
    pub fn observe_readout_cancellable(
        &self,
        window: HarmonicWindow,
        pairs: &[HarmonicReceiverPair],
        limits: HarmonicReadoutLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedHarmonicReadout<'_>, TransferError> {
        poll(&mut cancelled)?;
        let receivers = self.reference.receivers.len();
        admit(pairs.len(), MAX_READOUT_RECEIVER_PAIRS)?;
        admit(pairs.len(), limits.maximum_receiver_pairs)?;
        for pair in pairs {
            poll(&mut cancelled)?;
            if pair.first >= receivers || pair.second >= receivers {
                return Err(bad(
                    "BEM readout receiver index is outside the physical array",
                ));
            }
        }
        let lines = core::iter::once(&self.mean)
            .chain(self.factors.iter())
            .try_fold(0usize, |n, bank| n.checked_add(bank.lines.len()))
            .ok_or_else(overflow)?;
        let products = receivers
            .checked_mul(lines)
            .and_then(|n| n.checked_mul(8))
            .and_then(|n| {
                n.checked_add(
                    pairs
                        .len()
                        .checked_mul(self.factors.len())?
                        .checked_mul(16)?,
                )
            })
            .ok_or_else(overflow)?;
        admit(products, MAX_PRODUCTS)?;
        admit(products, limits.maximum_products)?;
        let workspace = receivers
            .checked_mul(
                core::mem::size_of::<[Dd; 8]>() + core::mem::size_of::<HarmonicAveragedField>(),
            )
            .and_then(|n| {
                n.checked_add(pairs.len().checked_mul(
                    core::mem::size_of::<[Dd; 16]>()
                        + core::mem::size_of::<HarmonicReadoutCovariance>(),
                )?)
            })
            .ok_or_else(overflow)?;
        admit(workspace, limits.maximum_workspace_bytes)?;

        let mut integrated = storage(receivers)?;
        integrated.resize(receivers, [Dd::ZERO; 8]);
        integrate_bank(&self.mean, window, &mut integrated, &mut cancelled)?;
        let mut means = storage(receivers)?;
        for field in &integrated {
            poll(&mut cancelled)?;
            means.push(HarmonicAveragedField {
                pressure_pa: readout_value(field[0])?,
                particle_velocity_m_s: [
                    readout_value(field[2])?,
                    readout_value(field[4])?,
                    readout_value(field[6])?,
                ],
            });
        }
        let mut sums = storage(pairs.len())?;
        sums.resize(pairs.len(), [Dd::ZERO; 16]);
        // Once integrated, each factor is a single proper complex amplitude at
        // each receiver. The zero-frequency kernel supplies the physical 1/2.
        let covariance = window
            .proper_product_kernel(0.0, 0.0)
            .map_err(window_error)?;
        for bank in &self.factors {
            poll(&mut cancelled)?;
            integrated.fill([Dd::ZERO; 8]);
            integrate_bank(bank, window, &mut integrated, &mut cancelled)?;
            for (pair, sum) in pairs.iter().zip(&mut sums) {
                poll(&mut cancelled)?;
                let a = &integrated[pair.first];
                let b = &integrated[pair.second];
                for row in 0..4 {
                    for col in 0..4 {
                        sum[4 * row + col] = sum[4 * row + col]
                            + covariance
                                .mean_real_product_parts_dd(component(a, row), component(b, col))
                                .map_err(window_error)?;
                    }
                }
            }
        }
        let mut covariances = storage(pairs.len())?;
        for (&pair, sum) in pairs.iter().zip(sums) {
            poll(&mut cancelled)?;
            let mut values = [0.0; 16];
            for (out, value) in values.iter_mut().zip(sum) {
                *out = readout_value(value)?;
            }
            if pair.first == pair.second && (0..4).any(|c| values[5 * c] < 0.0) {
                return Err(bad("BEM averaged readout variance became negative"));
            }
            covariances.push(HarmonicReadoutCovariance {
                receivers: pair,
                pressure_covariance_pa2: values[0],
                first_pressure_second_velocity_covariance_w_m2: [values[1], values[2], values[3]],
                first_velocity_second_pressure_covariance_w_m2: [values[4], values[8], values[12]],
                velocity_covariance_m2_s2: core::array::from_fn(|row| {
                    core::array::from_fn(|col| values[4 * (row + 1) + col + 1])
                }),
            });
        }
        poll(&mut cancelled)?;
        Ok(WindowedHarmonicReadout {
            window,
            reference: self.reference,
            source_clock: self.source_clock,
            means,
            covariances,
        })
    }
}

fn integrate_bank(
    bank: &Bank,
    window: HarmonicWindow,
    integrated: &mut [[Dd; 8]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    for line in &bank.lines {
        poll(cancelled)?;
        let kernel = window
            .product_kernel(line.omega, 0.0)
            .map_err(window_error)?;
        for (r, (sum, field)) in integrated.iter_mut().zip(&line.fields).enumerate() {
            if r % 64 == 0 {
                poll(cancelled)?;
            }
            for c in 0..4 {
                let a = component(field, c);
                sum[2 * c] = sum[2 * c]
                    + kernel
                        .mean_real_product_parts_dd(a, [Dd::ONE, Dd::ZERO])
                        .map_err(window_error)?;
                // Im(A Phi) = Re((-i A) Phi), with both low components kept.
                sum[2 * c + 1] = sum[2 * c + 1]
                    + kernel
                        .mean_real_product_parts_dd([a[1], -a[0]], [Dd::ONE, Dd::ZERO])
                        .map_err(window_error)?;
            }
        }
    }
    for field in integrated {
        if field.iter().any(|&value| !finite(value)) {
            return Err(bad("BEM integrated readout amplitude is unrepresentable"));
        }
    }
    Ok(())
}

fn readout_value(value: Dd) -> Result<f64, TransferError> {
    let rounded = value.to_f64();
    if finite(value) && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(bad("BEM averaged readout moment is unrepresentable"))
    }
}

#[cfg(test)]
mod tests;
