//! Spatial incident-power absorption after multiple scattering.
//!
//! A history samples the existing exact extinction free flight. At a collision
//! it deposits `P * kappa_a / (kappa_a + kappa_s)` in that original cell and
//! continues `P * kappa_s / (kappa_a + kappa_s)` through the native phase law.
//! This is implicit capture: scattering changes the next path, never creates
//! heat. Each wavelength has independent monochromatic histories, avoiding
//! spectral mixture weights that can exceed the source's physical power.
//!
//! This open, index-one grid has no reflecting surfaces. The source is a narrow
//! collimated beam with already quadrature-weighted watts per lane. LTE emission
//! is independent and is not traced, subtracted or counted as incident heating.
//! A collision cap or unrepresentable split retains surviving power as unresolved radiation. The mean
//! unresolved power estimates the omitted tail; it is not a certified bound on
//! that expectation. Sample standard errors are fixed-sample estimates, not
//! confidence intervals or optional-stopping guarantees. In particular, an
//! unobserved rare event can have a zero empirical error estimate.

use super::super::{cells::validate_ray, checkpoint, deposition::removed_power};
use super::prepared::{self, Monochromatic, PaletteOptics, PreparedOptics, UniformOptics};
use super::{ScatteringGrid, VolumeTransportError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;
use fs_rand::StreamKey;
use std::collections::HashMap;

/// Caller-owned sampling, work and retained spatial-output limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ScatteringDepositionConfig {
    /// Independent histories per positive, scattering spectral lane; at least 2.
    pub histories_per_lane: u32,
    /// Capture at most this many collisions; retain the remaining power as
    /// unresolved radiation. No Russian roulette or unreported cutoff is used.
    pub max_collisions_per_history: u32,
    /// Hard limit on all sampled free flights plus deterministic lane traversals.
    /// Exhaustion refuses the complete operation, without a partial result.
    pub max_total_segments: u64,
    /// Maximum distinct absorbing grid cells retained, across all lanes.
    /// Accumulator storage is O(this limit times the wavelength packet width).
    pub max_deposition_cells: usize,
    /// Study seed for the native counter-based generator.
    pub seed: u64,
    /// Stable source identity. Distinct sources should use distinct values.
    pub source_id: u64,
}

/// Fixed-sample power estimate and estimated standard error of its mean, in W.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerEstimate<const N: usize> {
    pub mean_w: [f64; N],
    pub standard_error_w: [f64; N],
}

impl<const N: usize> PowerEstimate<N> {
    fn zero() -> Self {
        Self {
            mean_w: [0.0; N],
            standard_error_w: [0.0; N],
        }
    }
}

/// Heat delivered to an original X-major density-cell index.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellAbsorbedPower<const N: usize> {
    pub cell_index: usize,
    /// Repeat visits within a history are reduced before computing variance.
    pub absorbed: PowerEstimate<N>,
}

/// Conservative sampled power balance for the full open grid.
///
/// Each lane satisfies incident = ballistic escape + diffuse escape + absorbed
/// + unresolved, to floating-point roundoff, for every retained sample average.
/// The absorbed total's error uses history totals, preserving correlations
/// between cells; summing cell variances would give an incorrect error.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringPowerDeposition<const N: usize> {
    pub medium_identity: ContentHash,
    pub wavelength_bounds_nm: (f64, f64),
    pub wavelengths_nm: [f64; N],
    pub incident_power_w: [f64; N],
    /// Power leaving without any scattering collision.
    pub ballistic_escape: PowerEstimate<N>,
    /// Power leaving after at least one scattering collision; this is not heat.
    pub diffuse_escape: PowerEstimate<N>,
    pub absorbed: PowerEstimate<N>,
    /// Power still in transport at the collision cap or when a positive split
    /// falls below representable watts. Neither is silently discarded as heat.
    pub unresolved: PowerEstimate<N>,
    /// Ascending original cell indices. Unobserved cells have no retained record.
    pub cells: Vec<CellAbsorbedPower<N>>,
    /// Zero means exact deterministic absorption-only transport in that lane.
    pub sampled_histories: [u32; N],
    /// Histories retained unresolved because a positive split underflowed.
    pub precision_limited_histories: [u32; N],
    pub traced_segments: u64,
    pub config: ScatteringDepositionConfig,
}

impl ScatteringGrid<'_> {
    /// Estimate spatial heating including reabsorption of scattered radiation.
    ///
    /// The ray and every continuation run until they leave this grid. Density,
    /// coefficients and phase law are the original immutable medium's owners.
    /// Exact absorption-only lanes reuse `deposit_beam_power` with no sampling
    /// error. The remaining lanes use independent seed/source/lane/history
    /// streams and deterministic reductions; hash-table iteration cannot affect
    /// the results. Every free flight retains the grid's DDA work budget and
    /// cancellation checkpoints. No report is published on any refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_multiple_scattering_power<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, unit_direction, f64::INFINITY)?;
        if self.absorption.absorption_palette().is_some() || self.scattering_palette().is_some() {
            let medium = PaletteOptics::try_new(self, *wavelengths_nm, &mut || checkpoint(cx))?;
            return deposit_prepared(
                &medium,
                &mut || checkpoint(cx),
                origin,
                unit_direction,
                incident_power_w,
                config,
            );
        }
        let (absorption, scattering, extinction) = self.coefficients(wavelengths_nm)?;
        let medium = UniformOptics {
            medium: self,
            wavelengths: *wavelengths_nm,
            absorption,
            scattering,
            extinction,
        };
        deposit_prepared(
            &medium,
            &mut || checkpoint(cx),
            origin,
            unit_direction,
            incident_power_w,
            config,
        )
    }
}

/// One implicit-capture estimator for uniform and genuinely cell-local optics.
pub(super) fn deposit_prepared<const N: usize>(
    medium: &impl PreparedOptics<N>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    origin: Point3,
    unit_direction: Vec3,
    incident_power_w: &[f64; N],
    config: ScatteringDepositionConfig,
) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
    let wavelengths_nm = medium.wavelengths();
    check()?;
    validate_ray(origin, unit_direction, f64::INFINITY)?;
    if config.histories_per_lane < 2 || config.max_collisions_per_history == 0 {
        return Err(VolumeTransportError::Budget);
    }
    for (lane, &power) in incident_power_w.iter().enumerate() {
        if !power.is_finite() || power < 0.0 {
            return Err(VolumeTransportError::BeamPower { lane });
        }
    }
    let mut result = ScatteringPowerDeposition {
        medium_identity: medium.identity(),
        wavelength_bounds_nm: medium.bounds(),
        wavelengths_nm: *wavelengths_nm,
        incident_power_w: *incident_power_w,
        ballistic_escape: PowerEstimate::zero(),
        diffuse_escape: PowerEstimate::zero(),
        absorbed: PowerEstimate::zero(),
        unresolved: PowerEstimate::zero(),
        cells: Vec::new(),
        sampled_histories: [0; N],
        precision_limited_histories: [0; N],
        traced_segments: 0,
        config,
    };
    let mut spatial = SpatialTallies::<N>::new(config.max_deposition_cells);
    for lane in 0..N {
        check()?;
        let source_w = incident_power_w[lane];
        if source_w == 0.0 {
            continue;
        }
        let wavelength = [wavelengths_nm[lane]];
        let monochromatic = Monochromatic::<_, N> {
            medium,
            lane,
            wavelength,
        };
        if !medium.lane_has_scattering(lane) {
            count_segment(&mut result.traced_segments, config.max_total_segments)?;
            let exact = super::super::deposition::deposit_prepared(
                check,
                medium.geometry(),
                (
                    origin,
                    unit_direction,
                    if medium.lane_has_extinction(lane) {
                        f64::INFINITY
                    } else {
                        0.0
                    },
                ),
                &wavelength,
                &[source_w],
                config.max_deposition_cells,
                (medium.identity(), medium.bounds()),
                |cell| monochromatic.coefficients(cell),
            )?;
            result.ballistic_escape.mean_w[lane] = exact.transmitted_power_w[0];
            result.absorbed.mean_w[lane] = exact.absorbed_power_w[0];
            for cell in exact.cells {
                spatial.cell(cell.cell_index)?.lanes[lane]
                    .moments
                    .observe(cell.absorbed_power_w[0])?;
            }
            continue;
        }
        result.sampled_histories[lane] = config.histories_per_lane;
        let seed = lane_seed(config.seed, config.source_id, lane);
        let mut totals = [Moments::default(); 4];
        for history in 0..config.histories_per_lane {
            check()?;
            let mut random = StreamKey {
                seed,
                kernel: 0x4142_5348, // "ABSH": absorbed source-power histories.
                tile: history,
            }
            .stream();
            let history_result = trace_history(
                &monochromatic,
                check,
                (origin, unit_direction, source_w),
                (config.max_collisions_per_history, config.max_total_segments),
                &mut result.traced_segments,
                &mut || random.next_f64(),
                &mut |cell, heat| spatial.cell(cell)?.lanes[lane].capture(history, heat),
            )?;
            result.precision_limited_histories[lane] += u32::from(history_result.precision_limited);
            let outcome = history_result.values;
            for (moment, value) in totals.iter_mut().zip(outcome) {
                moment.observe(value)?;
            }
        }
        for (estimate, moment) in [
            &mut result.ballistic_escape,
            &mut result.diffuse_escape,
            &mut result.absorbed,
            &mut result.unresolved,
        ]
        .into_iter()
        .zip(totals)
        {
            (estimate.mean_w[lane], estimate.standard_error_w[lane]) = moment.estimate()?;
        }
    }
    result.cells = spatial.finish(check, &result.sampled_histories)?;
    check()?;
    Ok(result)
}

/// One complete implicit-capture history. Input packet units may be watts for
/// incident beams or dimensionless for reciprocal emission probabilities.
/// Source selection happens before this call, so existing beam draws are exact.
pub(super) struct HistoryOutcome {
    /// Ballistic escape, diffuse escape, absorption, unresolved continuation.
    pub values: [f64; 4],
    pub precision_limited: bool,
}
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_history(
    medium: &impl PreparedOptics<1>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    source: (Point3, Vec3, f64),
    limits: (u32, u64),
    traced_segments: &mut u64,
    random: &mut impl FnMut() -> f64,
    capture: &mut impl FnMut(usize, f64) -> Result<(), VolumeTransportError>,
) -> Result<HistoryOutcome, VolumeTransportError> {
    let mut precision_limited = false;
    let (mut position, mut direction, source_w) = source;
    let mut power = source_w;
    let mut survival = PowerSum {
        sum: source_w,
        correction: 0.0,
    };
    let mut captured = PowerSum::default();
    let mut outcome = [0.0; 4]; // ballistic, diffuse, heat, unresolved
    for collision in 0..limits.0 {
        count_segment(traced_segments, limits.1)?;
        let flight = prepared::sample_extinction(
            medium,
            check,
            position,
            direction,
            f64::INFINITY,
            0.0,
            random(),
            false,
        )?;
        let Some(event) = flight.event else {
            outcome[usize::from(collision != 0)] = power;
            break;
        };
        let (a, s, t) = medium.coefficient(event.cell_index, 0)?;
        let heat = split_power(power, a, t)?;
        let continuation = split_power(power, s, t)?;
        let (Some(heat), Some(continuation)) = (heat, continuation) else {
            // Retain the entire representable packet before any
            // partial capture, including the below-resolution tail.
            outcome[3] = power;
            precision_limited = true;
            break;
        };
        let roundoff = 32.0 * f64::EPSILON * power + 8.0 * f64::from_bits(1);
        if ((power - heat) - continuation).abs() > roundoff {
            return Err(VolumeTransportError::Unrepresentable);
        }
        if heat > 0.0 {
            captured.add(heat)?;
            capture(event.cell_index, heat)?;
        }
        survival.continue_after_capture(heat, continuation)?;
        power = survival.value();
        if power == 0.0 {
            break;
        }
        if collision + 1 == limits.0 {
            outcome[3] = power;
            break;
        }
        // A monochromatic matched phase sample has physical/PDF
        // weight one, including interpolated spectral phase tables.
        direction = medium
            .phase(event.cell_index)
            .sample_packet(medium.wavelengths(), direction, random(), random())?
            .direction;
        position = event.position;
    }
    outcome[2] = captured.value();
    Ok(HistoryOutcome {
        values: outcome,
        precision_limited,
    })
}

fn count_segment(count: &mut u64, maximum: u64) -> Result<(), VolumeTransportError> {
    if *count == maximum {
        return Err(VolumeTransportError::Budget);
    }
    *count += 1;
    Ok(())
}

fn lane_seed(seed: u64, source: u64, lane: usize) -> u64 {
    let mut hash = DomainHasher::new("org.frankensim.render.absorbed-power-stream.v1");
    hash.update(&seed.to_le_bytes());
    hash.update(&source.to_le_bytes());
    hash.update(&(lane as u64).to_le_bytes());
    let identity = hash.finalize();
    u64::from_le_bytes(
        identity.as_bytes()[..8]
            .try_into()
            .expect("eight hash bytes"),
    )
}

fn split_power(
    power: f64,
    coefficient: f64,
    extinction: f64,
) -> Result<Option<f64>, VolumeTransportError> {
    // Reuse the same stable coefficient split as exact beam absorption,
    // including its logarithmic fallback for very small coefficient ratios.
    let value = removed_power(power, 1.0, 1.0, 1.0, coefficient, extinction);
    if !value.is_finite() {
        Err(VolumeTransportError::Unrepresentable)
    } else if coefficient > 0.0 && value == 0.0 {
        Ok(None)
    } else {
        Ok(Some(value))
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct PowerSum {
    sum: f64,
    correction: f64,
}
impl PowerSum {
    fn add(&mut self, value: f64) -> Result<(), VolumeTransportError> {
        let increment = value - self.correction;
        let next = self.sum + increment;
        if !next.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        self.correction = (next - self.sum) - increment;
        self.sum = next;
        Ok(())
    }
    fn value(self) -> f64 {
        self.sum
    }

    fn continue_after_capture(
        &mut self,
        heat: f64,
        continuation: f64,
    ) -> Result<(), VolumeTransportError> {
        if heat <= 0.5 * self.sum {
            // Near unit albedo, multiplication may round the continuation
            // back to its incoming value. Retain small successive captures
            // with compensated subtraction of exactly the heat we tallied.
            self.add(-heat)
        } else {
            // When most energy is absorbed, direct subtraction could erase
            // a representable tail. Keep the stable scattering-ratio result.
            self.sum = continuation;
            self.correction = 0.0;
            Ok(())
        }
    }
}

/// Welford moments normalized by the largest observed value, so neither large
/// source watts nor very small deposits overflow/underflow the squared moment.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Moments {
    count: u32,
    scale: f64,
    mean: f64,
    m2: f64,
}
impl Moments {
    pub(super) fn observe(&mut self, value: f64) -> Result<(), VolumeTransportError> {
        if !value.is_finite() || value < 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        if value > self.scale {
            let factor = self.scale / value;
            self.mean *= factor;
            self.m2 *= factor * factor;
            self.scale = value;
        }
        let normalized = if self.scale == 0.0 {
            0.0
        } else {
            value / self.scale
        };
        self.count += 1;
        let delta = normalized - self.mean;
        self.mean += delta / f64::from(self.count);
        self.m2 += delta * (normalized - self.mean);
        Ok(())
    }
    fn zeroes_until(&mut self, total: u32) {
        if total > self.count {
            let old = f64::from(self.count);
            let zeroes = f64::from(total - self.count);
            let new = f64::from(total);
            self.m2 += self.mean * self.mean * old * (zeroes / new);
            self.mean *= old / new;
            self.count = total;
        }
    }
    pub(super) fn estimate(self) -> Result<(f64, f64), VolumeTransportError> {
        let mean = self.mean * self.scale;
        let error = if self.count < 2 {
            0.0
        } else {
            self.scale
                * det::sqrt(self.m2.max(0.0) / f64::from(self.count) / f64::from(self.count - 1))
        };
        if mean.is_finite() && error.is_finite() {
            Ok((mean, error))
        } else {
            Err(VolumeTransportError::Unrepresentable)
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct CellStatistics {
    moments: Moments,
    history: Option<u32>,
    power: PowerSum,
}
impl CellStatistics {
    fn capture(&mut self, history: u32, power: f64) -> Result<(), VolumeTransportError> {
        if self.history != Some(history) {
            self.flush()?;
            self.history = Some(history);
        }
        self.power.add(power)
    }
    fn flush(&mut self) -> Result<(), VolumeTransportError> {
        if let Some(history) = self.history.take() {
            self.moments.zeroes_until(history);
            self.moments.observe(self.power.value())?;
            self.power = PowerSum::default();
        }
        Ok(())
    }
}

struct CellTally<const N: usize> {
    index: usize,
    lanes: [CellStatistics; N],
}
pub(super) struct SpatialTallies<const N: usize> {
    cells: Vec<CellTally<N>>,
    // Lookup only: randomized hash iteration never controls reduction order.
    indices: HashMap<usize, usize>,
    limit: usize,
}
impl<const N: usize> SpatialTallies<N> {
    pub(super) fn capture(
        &mut self,
        cell: usize,
        lane: usize,
        history: u32,
        value: f64,
    ) -> Result<(), VolumeTransportError> {
        self.cell(cell)?.lanes[lane].capture(history, value)
    }
    pub(super) fn new(limit: usize) -> Self {
        Self {
            cells: Vec::new(),
            indices: HashMap::new(),
            limit,
        }
    }
    fn cell(&mut self, index: usize) -> Result<&mut CellTally<N>, VolumeTransportError> {
        if let Some(&slot) = self.indices.get(&index) {
            return Ok(&mut self.cells[slot]);
        }
        if self.cells.len() == self.limit {
            return Err(VolumeTransportError::Budget);
        }
        if self.cells.len() == self.cells.capacity() {
            let amount = self
                .cells
                .capacity()
                .max(16)
                .min(self.limit - self.cells.len());
            self.cells
                .try_reserve_exact(amount)
                .map_err(|_| VolumeTransportError::Budget)?;
            self.indices
                .try_reserve(amount)
                .map_err(|_| VolumeTransportError::Budget)?;
        }
        let slot = self.cells.len();
        self.cells.push(CellTally {
            index,
            lanes: [CellStatistics::default(); N],
        });
        self.indices.insert(index, slot);
        Ok(&mut self.cells[slot])
    }
    pub(super) fn finish(
        mut self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        samples: &[u32; N],
    ) -> Result<Vec<CellAbsorbedPower<N>>, VolumeTransportError> {
        check()?;
        self.cells.sort_unstable_by_key(|cell| cell.index);
        let mut output = Vec::new();
        output
            .try_reserve_exact(self.cells.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for mut cell in self.cells {
            check()?;
            let mut absorbed = PowerEstimate::zero();
            for lane in 0..N {
                cell.lanes[lane].flush()?;
                cell.lanes[lane].moments.zeroes_until(samples[lane]);
                (absorbed.mean_w[lane], absorbed.standard_error_w[lane]) =
                    cell.lanes[lane].moments.estimate()?;
            }
            output.push(CellAbsorbedPower {
                cell_index: cell.index,
                absorbed,
            });
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests;
