//! Prepared total outgoing fields and positive angular radiation integrals.
//! Source columns are contracted before squaring, including quiet cancellations;
//! no source-power Gram matrix replaces the original complex fields.

use super::harmonic::{MAX_MODES, MAX_PRODUCTS, Sum};
use super::*;
use fs_exec::CancelGate;
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

/// Explicit cold-preparation storage and work. The quadrature ladder uses the
/// original Gauss–Legendre/uniform-longitude grids, growing latitude count by
/// n/2+1 and alternating half-cell longitude offsets. At least three grids are
/// required. Runtime convergence is checked for each actual drive/coherence,
/// never inferred from separate source self-powers.
#[derive(Clone, Copy, Debug)]
pub struct PreparedRadiationOptions {
    pub angular: helmholtz::DirectivityProjectionOptions,
    /// At least eight latitudes on the first grid; independent of SH degree.
    pub initial_latitudes: usize,
    pub refinements: usize,
    /// All source/direction kernel evaluations, including requested observation
    /// directions and every cached angular grid, admitted before evaluation.
    pub maximum_kernel_evaluations: usize,
    /// Retained complex source-column entries across all grids and directions;
    /// capped additionally by MAX_OBSERVATION_ENTRIES. This is not an RSS limit.
    pub maximum_source_entries: usize,
}

/// Full outgoing time/ensemble-averaged power at one physical frequency.
/// For peak exp(-iωt) amplitudes, dP/dΩ=E[|F|²]/(2ρc), p~F exp(ikr)/r.
/// These are line powers, not spectral densities. Wall absorption and source
/// self-loading are separate quantities and are not debited here.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiatedPowerMoments {
    pub total_power_w: f64,
    /// Integral direction_i dP [W], signed. Dividing by c gives the net
    /// outgoing acoustic momentum rate, not a finite-body force solution.
    pub first_directional_moment_w: [f64; 3],
    /// Integral direction_i direction_j dP [W], symmetric; its trace is P.
    pub second_directional_moment_w: [[f64; 3]; 3],
    pub evaluated_directions: usize,
    /// Last grid-to-grid change divided by the larger total power. Two
    /// successive changes must pass. These estimates do not bound continuum,
    /// panel-integration or boundary-solve error.
    pub relative_power_change: f64,
    pub relative_moment_change: f64,
}

#[derive(Clone, Debug)]
struct Grid {
    directions: Vec<[f64; 3]>,
    weights: Vec<f64>,
    columns: Vec<C64>,
}

/// Immutable source angular columns of one existing radiation owner. Arbitrary
/// requested directions are separate from the positive full-sphere quadrature.
/// Repeated contractions allocate nothing and perform no solve/kernel work.
/// Each result keeps the original source SI units and common phase reference.
#[derive(Clone, Debug)]
pub struct PreparedRadiationObservation {
    omega: f64,
    medium: Medium,
    sources: Vec<SourceId>,
    producer_hash: ContentHash,
    directions: Vec<[f64; 3]>,
    columns: Vec<C64>,
    grids: Vec<Grid>,
    angular_directions: usize,
    relative_tolerance: f64,
    kernel_evaluations: usize,
    power_scale: f64,
}

impl IterativeBoundaryResponse<'_, '_> {
    /// Reuse retained surface/compact-source solutions, including exact direct
    /// monopoles and force dipoles. No new FMM setup or source solve occurs.
    pub fn prepare_radiation(
        &self,
        directions: &[[f64; 3]],
        options: PreparedRadiationOptions,
        gate: &CancelGate,
    ) -> Result<PreparedRadiationObservation, IterativeObservationError> {
        let nodes = if self.surface().triangles().is_some() {
            16
        } else {
            1
        };
        let mut per_direction = self
            .surface()
            .areas()
            .len()
            .checked_mul(nodes)
            .and_then(|n| n.checked_mul(self.sources().len()))
            .ok_or(TransferError::Invalid("radiation kernel work overflow"))?;
        for index in 0..self.sources().len() {
            per_direction += usize::from(self.source_incident(index)?.is_some());
        }
        PreparedRadiationObservation::prepare(
            self.omega_rad_s(),
            self.medium(),
            self.sources(),
            self.operator_hash(),
            directions,
            options,
            per_direction,
            |directions| {
                self.far_field(directions, options.maximum_kernel_evaluations, gate)
                    .map(|sample| sample.amplitudes)
            },
            &|| gate.is_requested(),
        )
    }
}

impl BemObservation<'_> {
    /// Dense retained-solution counterpart, with the same exact source-trace,
    /// surface, medium and frequency admission as retained harmonic observation.
    /// Uses the original integrated far-field owner; no second solver is added.
    pub fn prepare_radiation_from_solutions(
        &self,
        omega: f64,
        solutions: &[helmholtz::RadiationSolution],
        directions: &[[f64; 3]],
        options: PreparedRadiationOptions,
        gate: &CancelGate,
    ) -> Result<PreparedRadiationObservation, TransferError> {
        self.admit_retained_solutions(omega, solutions, &mut || gate.is_requested())?;
        let nodes = if self.surface.triangles().is_some() {
            16
        } else {
            1
        };
        let per_direction = self
            .surface
            .areas()
            .len()
            .checked_mul(nodes)
            .and_then(|n| n.checked_mul(solutions.len()))
            .ok_or(TransferError::Invalid("radiation kernel work overflow"))?;
        PreparedRadiationObservation::prepare(
            omega,
            self.medium,
            &self.sources,
            self.operator_hash(),
            directions,
            options,
            per_direction,
            |directions| {
                let mut columns = reserve(solutions.len() * directions.len())?;
                for solution in solutions {
                    columns.extend(helmholtz::far_field_integrated_cancellable(
                        self.surface,
                        solution,
                        self.medium,
                        directions,
                        options.maximum_kernel_evaluations,
                        || gate.is_requested(),
                    )?);
                }
                Ok(columns)
            },
            &|| gate.is_requested(),
        )
    }
}

impl PreparedRadiationObservation {
    #[allow(clippy::too_many_arguments)]
    fn prepare<E: From<TransferError> + From<HelmholtzError>>(
        omega: f64,
        medium: Medium,
        sources: &[SourceId],
        producer_hash: ContentHash,
        directions: &[[f64; 3]],
        options: PreparedRadiationOptions,
        per_direction_work: usize,
        mut evaluate: impl FnMut(&[[f64; 3]]) -> Result<Vec<C64>, E>,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self, E> {
        let mut cancellation = || cancelled();
        check_cancel(&mut cancellation)?;
        harmonic::admit_frequency(omega, medium)?;
        helmholtz::admit_angular_options(options.angular)?;
        cap(sources.len(), MAX_OBSERVATION_INPUTS)?;
        cap(directions.len(), helmholtz::MAX_SH_QUADRATURE_DIRECTIONS)?;
        if sources.is_empty()
            || options.initial_latitudes < 8
            || !(3..=16).contains(&options.refinements)
        {
            return Err(TransferError::Invalid(
                "radiation bank needs sources and at least three angular grids",
            )
            .into());
        }
        for (index, &direction) in directions.iter().enumerate() {
            if index % 256 == 0 {
                check_cancel(&mut cancellation)?;
            }
            let length = det::sqrt(direction.iter().map(|v| v * v).sum::<f64>());
            if !length.is_finite() || (length - 1.0).abs() > 1e-10 {
                return Err(TransferError::Invalid(
                    "radiation observation directions must be finite unit vectors",
                )
                .into());
            }
        }
        let mut latitudes = reserve(options.refinements)?;
        let mut latitude = options.initial_latitudes;
        let mut angular_directions = 0usize;
        for _ in 0..options.refinements {
            let count = latitude
                .checked_mul(latitude.saturating_mul(2).saturating_add(1))
                .ok_or(TransferError::Invalid("angular grid shape overflow"))?;
            angular_directions = angular_directions
                .checked_add(count)
                .ok_or(TransferError::Invalid("angular ladder shape overflow"))?;
            cap(angular_directions, options.angular.maximum_directions)?;
            latitudes.push(latitude);
            latitude = latitude
                .checked_add(latitude / 2 + 1)
                .ok_or(TransferError::Invalid("angular latitude overflow"))?;
        }
        let all_directions = angular_directions + directions.len();
        let entries = all_directions
            .checked_mul(sources.len())
            .ok_or(TransferError::Invalid(
                "radiation source bank extent overflow",
            ))?;
        cap(
            entries,
            options.maximum_source_entries.min(MAX_OBSERVATION_ENTRIES),
        )?;
        let work = per_direction_work.checked_mul(all_directions);
        if work.is_none_or(|n| n > options.maximum_kernel_evaluations) {
            return Err(HelmholtzError::FarFieldBudget {
                required_evaluations: work.unwrap_or(usize::MAX),
                maximum_evaluations: options.maximum_kernel_evaluations,
            }
            .into());
        }
        let columns = if directions.is_empty() {
            Vec::new()
        } else {
            evaluate(directions)?
        };
        admit_columns(&columns, sources.len() * directions.len())?;
        let mut grids = reserve(options.refinements)?;
        for (index, latitude) in latitudes.into_iter().enumerate() {
            check_cancel(&mut cancellation)?;
            let grid = helmholtz::angular_grid(latitude, index % 2 == 1, cancelled)?;
            let columns = evaluate(&grid.directions)?;
            admit_columns(&columns, sources.len() * grid.directions.len())?;
            let mut weights = reserve(grid.directions.len())?;
            let longitude_weight = core::f64::consts::TAU / grid.phases.len() as f64;
            for weight in grid.weights {
                weights.extend(std::iter::repeat_n(
                    weight * longitude_weight,
                    grid.phases.len(),
                ));
            }
            grids.push(Grid {
                directions: grid.directions,
                weights,
                columns,
            });
        }
        check_cancel(&mut cancellation)?;
        Ok(Self {
            omega,
            medium,
            sources: sources.to_vec(),
            producer_hash,
            directions: directions.to_vec(),
            columns,
            grids,
            angular_directions,
            relative_tolerance: options.angular.relative_tolerance,
            kernel_evaluations: work.unwrap(),
            power_scale: core::f64::consts::FRAC_1_SQRT_2
                / det::sqrt(medium.density * medium.sound_speed),
        })
    }

    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn medium(&self) -> Medium {
        self.medium
    }
    pub fn sources(&self) -> &[SourceId] {
        &self.sources
    }
    pub fn directions(&self) -> &[[f64; 3]] {
        &self.directions
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.producer_hash
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.kernel_evaluations
    }
    pub fn cached_angular_directions(&self) -> usize {
        self.angular_directions
    }

    /// Requested-direction peak amplitudes only, without integrating power or
    /// requiring angular convergence. O(sources × requested directions), no
    /// allocation or kernel evaluation. Discard output after any error.
    pub fn observe_coherent(
        &self,
        drives: &[C64],
        budget: HarmonicWorkBudget,
        output: &mut [C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.admit(
            Input::Modes {
                factors: drives,
                modes: 1,
            },
            budget,
            output.len(),
            false,
        )?;
        for (direction, out) in output.iter_mut().enumerate() {
            check_cancel(&mut cancelled)?;
            *out = self.field(&self.columns, self.directions.len(), direction, drives)?;
        }
        check_cancel(&mut cancelled)
    }

    /// Requested-direction dP/dΩ only, with the same peak coherence factors as
    /// `evaluate_coherence_modes`. No full-sphere quadrature is evaluated.
    pub fn observe_coherence_modes(
        &self,
        factors: &[C64],
        modes: usize,
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.observe_power(
            Input::Modes { factors, modes },
            budget,
            output_w_per_sr,
            &mut cancelled,
        )
    }

    /// Requested-direction independent-source dP/dΩ in O(sources × directions).
    /// No diagonal factor matrix or angular power integration is constructed.
    pub fn observe_independent_sources(
        &self,
        peak_second_moments: &[f64],
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.observe_power(
            Input::Independent(peak_second_moments),
            budget,
            output_w_per_sr,
            &mut cancelled,
        )
    }

    /// Requested-direction radiation from measured source coherence in one
    /// interior Welch bin. Output is this bin's W/sr, not W/(sr Hz).
    /// Channels retain the original source kinds, order, SI units and clock.
    /// No kernel evaluation, boundary solve or full-sphere integration occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn observe_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[SourceId],
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        check_cancel(&mut cancelled)?;
        let view =
            harmonic::welch::source_bin(spectrum, bin, self.omega, &self.sources, recorded_order)?;
        self.observe_power(Input::Welch(view), budget, output_w_per_sr, &mut cancelled)
    }

    fn observe_power(
        &self,
        input: Input<'_>,
        budget: HarmonicWorkBudget,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.admit(input, budget, output.len(), false)?;
        for (direction, out) in output.iter_mut().enumerate() {
            check_cancel(cancelled)?;
            *out = self.power_density(
                &self.columns,
                self.directions.len(),
                direction,
                input,
                cancelled,
            )?;
        }
        check_cancel(cancelled)
    }

    /// Coherent peak far-field amplitudes in requested-direction order plus
    /// total time-averaged radiated power. No allocation after preparation.
    /// On any error, discard output; no converged integrated result is returned.
    pub fn evaluate_coherent(
        &self,
        drives: &[C64],
        budget: HarmonicWorkBudget,
        output: &mut [C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiatedPowerMoments, TransferError> {
        let input = Input::Modes {
            factors: drives,
            modes: 1,
        };
        self.admit(input, budget, output.len(), true)?;
        let power = self.integrate(input, &mut cancelled)?;
        for (direction, out) in output.iter_mut().enumerate() {
            check_cancel(&mut cancelled)?;
            *out = self.field(&self.columns, self.directions.len(), direction, drives)?;
        }
        check_cancel(&mut cancelled)?;
        Ok(power)
    }

    /// B is mode-major/source-minor; C=B B*. Each entry uses its source's SI
    /// coordinate. Modes are independent, sources within each mode interfere.
    /// Output is dP/dΩ [W/sr], never |sum of source powers| or a PSD in W/Hz.
    /// A deterministic peak drive is one mode; zero modes are silence.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[C64],
        modes: usize,
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiatedPowerMoments, TransferError> {
        self.evaluate_power(
            Input::Modes { factors, modes },
            budget,
            output_w_per_sr,
            &mut cancelled,
        )
    }

    /// Diagonal peak source second moments in O(sources × directions), without
    /// constructing a diagonal factor matrix. Entries use squared source SI units.
    pub fn evaluate_independent_sources(
        &self,
        peak_second_moments: &[f64],
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiatedPowerMoments, TransferError> {
        self.evaluate_power(
            Input::Independent(peak_second_moments),
            budget,
            output_w_per_sr,
            &mut cancelled,
        )
    }

    /// Full outgoing power and angular moments of one measured spectral bin.
    ///
    /// Retained source fields are combined coherently for each measured factor
    /// before positive angular power is formed. The ORIGINAL angular refinement
    /// and two-successive-grid convergence gate run on this actual covariance;
    /// converged source self-powers do not qualify a different interference
    /// pattern. The cache is reused without solves, quadrature kernels or query
    /// allocation. Work includes recorded-factor conversion and every cached
    /// direction needed by the original convergence ladder.
    ///
    /// Returns bin watts and W/sr, with bandwidth and peak/time normalization
    /// each applied once. DC/Nyquist, frequency and source-order mismatches
    /// refuse. The record does not validate stationarity, unobserved bands,
    /// independent segments, source physics or a continuum radiation error.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[SourceId],
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiatedPowerMoments, TransferError> {
        check_cancel(&mut cancelled)?;
        let view =
            harmonic::welch::source_bin(spectrum, bin, self.omega, &self.sources, recorded_order)?;
        self.evaluate_power(Input::Welch(view), budget, output_w_per_sr, &mut cancelled)
    }

    fn evaluate_power(
        &self,
        input: Input<'_>,
        budget: HarmonicWorkBudget,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RadiatedPowerMoments, TransferError> {
        self.admit(input, budget, output.len(), true)?;
        let power = self.integrate(input, cancelled)?;
        for (direction, out) in output.iter_mut().enumerate() {
            check_cancel(cancelled)?;
            *out = self.power_density(
                &self.columns,
                self.directions.len(),
                direction,
                input,
                cancelled,
            )?;
        }
        check_cancel(cancelled)?;
        Ok(power)
    }

    fn admit(
        &self,
        input: Input<'_>,
        budget: HarmonicWorkBudget,
        output: usize,
        integrate_power: bool,
    ) -> Result<(), TransferError> {
        if output != self.directions.len() {
            return Err(TransferError::Invalid("radiation output shape"));
        }
        cap(budget.maximum_products, MAX_PRODUCTS)?;
        let sources = self.sources.len();
        let products = match input {
            Input::Modes { factors, modes } => {
                cap(modes, MAX_MODES)?;
                if factors.len() != modes * sources || factors.iter().any(|&z| !finite(z)) {
                    return Err(TransferError::Invalid(
                        "radiation coherence-factor shape or value",
                    ));
                }
                modes * (sources + 1)
            }
            Input::Independent(powers) => {
                if powers.len() != sources || powers.iter().any(|&p| !p.is_finite() || p < 0.0) {
                    return Err(TransferError::Invalid(
                        "radiation source peak second moments",
                    ));
                }
                2 * sources
            }
            Input::Welch(view) => {
                cap(view.segment_count(), MAX_MODES)?;
                if view.channel_count() != sources {
                    return Err(TransferError::Invalid("radiation recorded source shape"));
                }
                // Include both scalar components of the checked peak conversion.
                view.segment_count() * (3 * sources + 1)
            }
        };
        let directions = self.directions.len()
            + if integrate_power {
                self.angular_directions
            } else {
                0
            };
        let work = directions
            .checked_mul(products + 10)
            .ok_or(TransferError::Invalid(
                "radiation contraction work overflow",
            ))?;
        cap(work, budget.maximum_products)
    }

    fn field(
        &self,
        columns: &[C64],
        directions: usize,
        direction: usize,
        drives: &[C64],
    ) -> Result<C64, TransferError> {
        self.field_with(columns, directions, direction, drives.len(), |source| {
            Ok(drives[source])
        })
    }

    fn field_with(
        &self,
        columns: &[C64],
        directions: usize,
        direction: usize,
        sources: usize,
        mut drive: impl FnMut(usize) -> Result<C64, TransferError>,
    ) -> Result<C64, TransferError> {
        let mut sum = Sum::default();
        for source in 0..sources {
            sum.add(columns[source * directions + direction] * drive(source)?);
        }
        let value = sum.finish();
        if !finite(value) {
            return Err(TransferError::Invalid("coherent outgoing field overflow"));
        }
        Ok(value)
    }

    fn power_density(
        &self,
        columns: &[C64],
        directions: usize,
        direction: usize,
        input: Input<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, TransferError> {
        let mut sum = Sum::default();
        match input {
            Input::Modes { factors, .. } => {
                for mode in factors.chunks_exact(self.sources.len()) {
                    check_cancel(cancelled)?;
                    let value = self
                        .field(columns, directions, direction, mode)?
                        .scale(self.power_scale);
                    sum.add(C64::from_re(value.norm_sq()));
                }
            }
            Input::Independent(powers) => {
                for (source, &power) in powers.iter().enumerate() {
                    if source % 256 == 0 {
                        check_cancel(cancelled)?;
                    }
                    let value = columns[source * directions + direction]
                        .scale(det::sqrt(power))
                        .scale(self.power_scale);
                    sum.add(C64::from_re(value.norm_sq()));
                }
            }
            Input::Welch(view) => {
                for factor in view.factors() {
                    check_cancel(cancelled)?;
                    let value = self
                        .field_with(columns, directions, direction, factor.len(), |source| {
                            if source % 256 == 0 {
                                check_cancel(cancelled)?;
                            }
                            harmonic::welch::peak(view, factor[source])
                        })?
                        .scale(self.power_scale);
                    sum.add(C64::from_re(value.norm_sq()));
                }
            }
        }
        let power = sum.finish().re;
        if !power.is_finite() || power < 0.0 {
            return Err(TransferError::Invalid("outgoing angular power overflow"));
        }
        Ok(power)
    }

    fn integrate(
        &self,
        input: Input<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RadiatedPowerMoments, TransferError> {
        check_cancel(cancelled)?;
        let mut previous: Option<[f64; 10]> = None;
        let mut settled = 0;
        let mut evaluated = 0;
        let mut changes = None;
        for grid in &self.grids {
            let mut sums = [Sum::default(); 10];
            for (index, (&d, &weight)) in grid.directions.iter().zip(&grid.weights).enumerate() {
                check_cancel(cancelled)?;
                let power = self.power_density(
                    &grid.columns,
                    grid.directions.len(),
                    index,
                    input,
                    cancelled,
                )? * weight;
                let factors = [
                    1.0,
                    d[0],
                    d[1],
                    d[2],
                    d[0] * d[0],
                    d[0] * d[1],
                    d[0] * d[2],
                    d[1] * d[1],
                    d[1] * d[2],
                    d[2] * d[2],
                ];
                for (sum, factor) in sums.iter_mut().zip(factors) {
                    sum.add(C64::from_re(power * factor));
                }
            }
            evaluated += grid.directions.len();
            let current = sums.map(|sum| sum.finish().re);
            if current.iter().any(|v| !v.is_finite()) || current[0] < 0.0 {
                return Err(TransferError::Invalid("integrated outgoing power overflow"));
            }
            if let Some(old) = previous {
                let scale = current[0].max(old[0]);
                let (power, moments) = if scale == 0.0 {
                    (0.0, 0.0)
                } else {
                    (
                        (current[0] - old[0]).abs() / scale,
                        (1..10)
                            .map(|i| (current[i] / scale - old[i] / scale).abs())
                            .fold(0.0_f64, f64::max),
                    )
                };
                changes = Some((power, moments));
                settled = if power <= self.relative_tolerance && moments <= self.relative_tolerance
                {
                    settled + 1
                } else {
                    0
                };
                if settled >= 2 {
                    check_cancel(cancelled)?;
                    return Ok(RadiatedPowerMoments {
                        total_power_w: current[0],
                        first_directional_moment_w: [current[1], current[2], current[3]],
                        second_directional_moment_w: [
                            [current[4], current[5], current[6]],
                            [current[5], current[7], current[8]],
                            [current[6], current[8], current[9]],
                        ],
                        evaluated_directions: evaluated,
                        relative_power_change: power,
                        relative_moment_change: moments,
                    });
                }
            }
            previous = Some(current);
        }
        Err(TransferError::AngularPowerNotConverged {
            evaluated_directions: evaluated,
            power_change: changes.map(|c| c.0),
            moment_change: changes.map(|c| c.1),
        })
    }
}

#[derive(Clone, Copy)]
enum Input<'a> {
    Modes { factors: &'a [C64], modes: usize },
    Independent(&'a [f64]),
    Welch(WelchPeakBin<'a>),
}

fn admit_columns(values: &[C64], expected: usize) -> Result<(), TransferError> {
    if values.len() != expected || values.iter().any(|&z| !finite(z)) {
        Err(TransferError::Invalid(
            "incomplete or nonfinite outgoing source columns",
        ))
    } else {
        Ok(())
    }
}
fn reserve<T>(count: usize) -> Result<Vec<T>, TransferError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| HelmholtzError::AllocationFailed {
            what: "prepared radiation bank",
        })?;
    Ok(values)
}
fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "radiation/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "radiation/welch_tests.rs"]
mod welch_tests;
