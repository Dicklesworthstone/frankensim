//! Directed transport of nonthermal, isotropic volumetric emission.
//!
//! Each row launches unit-energy packets uniformly within one original cell,
//! using the existing exact extinction/free-flight and implicit-capture owner.
//! Absorption, including return to the emitting cell, is retained separately
//! from escape and unresolved radiation. These dimensionless fractions apply
//! to line/bin-integrated watts or joules; there is no Planck factor, wavelength
//! quadrature, reciprocity symmetrization, or assumption that a source absorbs.
//!
//! The operator describes one passage through an immutable, open, index-one
//! optical snapshot. Scattering is elastic. The caller owns excitation,
//! wavelength redistribution, delayed reemission and conversion to heat. Photon
//! flight time is neglected. Uniform cell emission is the spatial assumption,
//! including when a receiver is subsequently used as a recycling source.
//!
//! Every sampled history conserves energy to roundoff. Fixed-sample standard
//! errors are diagnostics, not confidence bounds; an unobserved rare outcome
//! can have zero empirical error. Neither unresolved radiation nor sampling
//! error may be silently assigned to heat. Receiver links within a row have
//! correlated errors. A row's total absorption error is measured from
//! history totals, not obtained by adding receiver variances.

use super::super::{VolumeTransportError, checkpoint, mixture::PositiveSum};
use super::ScatteringGrid;
use super::deposition::{Moments, SpatialTallies, trace_history};
use super::exchange::emitted_ray;
use super::field::CellScatteringGrid;
use super::prepared::{Monochromatic, PaletteOptics, PreparedOptics, UniformOptics};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_rand::StreamKey;

/// Sampling and storage limits for the complete operator, not renewed per row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LuminescenceTransportConfig {
    /// Independent histories for every supplied emitting cell and lane; >= 2.
    pub histories_per_emitter_lane: u32,
    /// Surviving packets at this cap remain unresolved, never heat.
    pub max_collisions_per_history: u32,
    /// All sampled free flights over all emitting cells and lanes.
    pub max_total_segments: u64,
    /// All positive-density cells, even those with no absorption or emission.
    pub max_material_cells: usize,
    pub max_emitter_cells: usize,
    /// All distinct directed emitter/receiver links, summed over rows. A link
    /// holds all spectral lanes, and a self-absorption link counts once.
    pub max_absorption_links: usize,
    pub seed: u64,
    pub source_id: u64,
}

/// Dimensionless fraction of emitted radiant energy and estimated error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FractionEstimate<const M: usize> {
    pub mean: [f64; M],
    pub standard_error: [f64; M],
}

impl<const M: usize> FractionEstimate<M> {
    fn zero() -> Self {
        Self {
            mean: [0.0; M],
            standard_error: [0.0; M],
        }
    }

    fn validate(&self, allowance: f64) -> Result<(), VolumeTransportError> {
        if self
            .mean
            .iter()
            .any(|value| !nonnegative(*value) || *value - 1.0 > allowance)
            || self.standard_error.iter().any(|value| !nonnegative(*value))
        {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LuminescenceAbsorption<const M: usize> {
    /// Original X-major grid index, including the source for self-absorption.
    pub cell_index: usize,
    pub fraction: FractionEstimate<M>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LuminescenceTransferRow<const M: usize> {
    pub emitter_cell: usize,
    /// Ascending original receiver indices. Unobserved cells have no record.
    pub absorption: Vec<LuminescenceAbsorption<M>>,
    pub absorbed: FractionEstimate<M>,
    /// Ballistic plus scattered radiation leaving the open grid.
    pub escape: FractionEstimate<M>,
    /// Surviving radiation at a collision cap or unrepresentable positive split.
    pub unresolved: FractionEstimate<M>,
    pub precision_limited_histories: [u32; M],
}

#[derive(Clone, Debug, PartialEq)]
pub struct LuminescenceTransferOperator<const M: usize> {
    /// Full absorption/scattering/phase snapshot identity.
    pub medium_identity: ContentHash,
    /// Original AbsorbingGrid owner. A CellScatteringGrid may replace its
    /// spectra; in that case medium_identity binds those replacement optics.
    pub absorption_identity: ContentHash,
    pub wavelength_bounds_nm: (f64, f64),
    pub wavelengths_nm: [f64; M],
    pub grid_cell_count: usize,
    /// Ascending original positive-density cells, including non-emitting cells
    /// whose temperature can invalidate the fixed optical snapshot.
    pub material_cells: Vec<usize>,
    /// Ascending original emitter indices, independent of local absorptivity.
    pub rows: Vec<LuminescenceTransferRow<M>>,
    pub traced_segments: u64,
    pub config: LuminescenceTransportConfig,
}

impl<const M: usize> LuminescenceTransferOperator<M> {
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        &self.wavelengths_nm
    }

    /// Structural and conservative admission for a retained/public report.
    /// Consumers must also match identities and original material-cell bindings
    /// to their optical snapshot. This does not certify Monte Carlo accuracy.
    ///
    /// `fraction_roundoff` is the caller's absolute dimensionless allowance for
    /// row closure and the independently reduced receiver/total absorption.
    /// It is never enlarged or replaced by statistical uncertainty. Fractions
    /// are never clipped, renormalized or symmetrized during admission.
    pub fn validate_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        fraction_roundoff: f64,
    ) -> Result<(), VolumeTransportError> {
        check()?;
        validate_wavelengths(&self.wavelengths_nm, self.wavelength_bounds_nm)?;
        if self.medium_identity == ContentHash([0; 32])
            || self.absorption_identity == ContentHash([0; 32])
        {
            return Err(VolumeTransportError::SourceIdentity);
        }
        if !nonnegative(fraction_roundoff) || self.grid_cell_count == 0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        validate_config(self.config)?;
        if self.material_cells.len() > self.config.max_material_cells
            || self.rows.len() > self.config.max_emitter_cells
            || self.traced_segments > self.config.max_total_segments
        {
            return Err(VolumeTransportError::Budget);
        }
        let minimum = minimum_segments::<M>(self.rows.len(), self.config)?;
        if self.traced_segments < minimum
            || self.traced_segments
                > minimum.saturating_mul(u64::from(self.config.max_collisions_per_history))
        {
            return Err(VolumeTransportError::Budget);
        }
        let mut previous = None;
        for &cell in &self.material_cells {
            check()?;
            validate_cell(cell, previous, self.grid_cell_count)?;
            previous = Some(cell);
        }
        let mut previous = None;
        let mut links = 0_usize;
        for row in &self.rows {
            check()?;
            validate_cell(row.emitter_cell, previous, self.grid_cell_count)?;
            previous = Some(row.emitter_cell);
            links = links
                .checked_add(row.absorption.len())
                .ok_or(VolumeTransportError::Budget)?;
            if links > self.config.max_absorption_links {
                return Err(VolumeTransportError::Budget);
            }
            for fraction in [&row.absorbed, &row.escape, &row.unresolved] {
                fraction.validate(fraction_roundoff)?;
            }
            let mut sums: [PositiveSum; M] = core::array::from_fn(|_| PositiveSum::default());
            let mut previous_receiver = None;
            for receiver in &row.absorption {
                check()?;
                validate_cell(receiver.cell_index, previous_receiver, self.grid_cell_count)?;
                if self
                    .material_cells
                    .binary_search(&receiver.cell_index)
                    .is_err()
                {
                    return Err(VolumeTransportError::Grid);
                }
                previous_receiver = Some(receiver.cell_index);
                receiver.fraction.validate(fraction_roundoff)?;
                for (sum, value) in sums.iter_mut().zip(receiver.fraction.mean) {
                    sum.add(value)?;
                }
            }
            for (lane, sum) in sums.into_iter().enumerate() {
                if row.precision_limited_histories[lane] > self.config.histories_per_emitter_lane {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                if (sum.finish()? - row.absorbed.mean[lane]).abs() > fraction_roundoff {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                let mut closure = PositiveSum::default();
                for value in [
                    row.absorbed.mean[lane],
                    row.escape.mean[lane],
                    row.unresolved.mean[lane],
                ] {
                    closure.add(value)?;
                }
                if (closure.finish()? - 1.0).abs() > fraction_roundoff {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
        }
        check()
    }
}

impl ScatteringGrid<'_> {
    /// Prepare reusable directed probabilities for arbitrary nonthermal source
    /// cells. Source order is immaterial; duplicate/out-of-grid indices refuse.
    /// Zero-density or transparent cells can emit. Each supplied cell and lane
    /// uses independent counter-based histories, with deterministic reductions.
    pub fn estimate_luminescence_transfer<const M: usize>(
        &self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64; M],
        emitter_cells: &[usize],
        config: LuminescenceTransportConfig,
    ) -> Result<LuminescenceTransferOperator<M>, VolumeTransportError> {
        checkpoint(cx)?;
        if self.absorption.absorption_palette().is_some() || self.scattering_palette().is_some() {
            let prepared = PaletteOptics::try_new(self, *wavelengths_nm, &mut || checkpoint(cx))?;
            return estimate(&prepared, &mut || checkpoint(cx), emitter_cells, config);
        }
        let (absorption, scattering, extinction) = self.coefficients(wavelengths_nm)?;
        estimate(
            &UniformOptics {
                medium: self,
                wavelengths: *wavelengths_nm,
                absorption,
                scattering,
                extinction,
            },
            &mut || checkpoint(cx),
            emitter_cells,
            config,
        )
    }
}

impl<const M: usize> CellScatteringGrid<'_, '_, M> {
    pub fn estimate_luminescence_transfer(
        &self,
        cx: &Cx<'_>,
        emitter_cells: &[usize],
        config: LuminescenceTransportConfig,
    ) -> Result<LuminescenceTransferOperator<M>, VolumeTransportError> {
        self.estimate_luminescence_transfer_with_checkpoint(
            &mut || checkpoint(cx),
            emitter_cells,
            config,
        )
    }

    /// Hosting solvers may supply their own cancellation/budget checkpoint.
    /// A refusal never publishes a partially traced operator.
    pub fn estimate_luminescence_transfer_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        emitter_cells: &[usize],
        config: LuminescenceTransportConfig,
    ) -> Result<LuminescenceTransferOperator<M>, VolumeTransportError> {
        estimate(self, check, emitter_cells, config)
    }
}

fn estimate<const M: usize>(
    medium: &impl PreparedOptics<M>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    emitter_cells: &[usize],
    config: LuminescenceTransportConfig,
) -> Result<LuminescenceTransferOperator<M>, VolumeTransportError> {
    check()?;
    validate_wavelengths(medium.wavelengths(), medium.bounds())?;
    validate_config(config)?;
    if emitter_cells.len() > config.max_emitter_cells
        || minimum_segments::<M>(emitter_cells.len(), config)? > config.max_total_segments
    {
        return Err(VolumeTransportError::Budget);
    }
    let geometry = medium.geometry();
    let mut emitters = Vec::new();
    emitters
        .try_reserve_exact(emitter_cells.len())
        .map_err(|_| VolumeTransportError::Budget)?;
    emitters.extend_from_slice(emitter_cells);
    emitters.sort_unstable();
    let mut previous = None;
    for &cell in &emitters {
        check()?;
        validate_cell(cell, previous, geometry.cell_count())?;
        previous = Some(cell);
    }
    let mut output = LuminescenceTransferOperator {
        medium_identity: medium.identity(),
        absorption_identity: geometry.content_identity(),
        wavelength_bounds_nm: medium.bounds(),
        wavelengths_nm: *medium.wavelengths(),
        grid_cell_count: geometry.cell_count(),
        material_cells: Vec::new(),
        rows: Vec::new(),
        traced_segments: 0,
        config,
    };
    for cell in 0..geometry.cell_count() {
        check()?;
        if geometry.grid.data[cell] > 0.0 {
            reserve_one(&mut output.material_cells, config.max_material_cells)?;
            output.material_cells.push(cell);
        }
    }
    output
        .rows
        .try_reserve_exact(emitters.len())
        .map_err(|_| VolumeTransportError::Budget)?;
    let mut retained_links = 0;
    for emitter_cell in emitters {
        check()?;
        let mut row = LuminescenceTransferRow {
            emitter_cell,
            absorption: Vec::new(),
            absorbed: FractionEstimate::zero(),
            escape: FractionEstimate::zero(),
            unresolved: FractionEstimate::zero(),
            precision_limited_histories: [0; M],
        };
        let mut spatial = SpatialTallies::<M>::new(config.max_absorption_links - retained_links);
        for lane in 0..M {
            let mono = Monochromatic::<_, M> {
                medium,
                lane,
                wavelength: [output.wavelengths_nm[lane]],
            };
            let mut totals = [Moments::default(); 3];
            let seed = emitter_seed(config.seed, config.source_id, emitter_cell, lane);
            for history in 0..config.histories_per_emitter_lane {
                check()?;
                let mut stream = StreamKey {
                    seed,
                    kernel: 0x4C55_4D45, // LUME
                    tile: history,
                }
                .stream();
                let (position, direction) =
                    emitted_ray(geometry, emitter_cell, &mut || stream.next_f64())?;
                let outcome = trace_history(
                    &mono,
                    check,
                    (position, direction, 1.0),
                    (config.max_collisions_per_history, config.max_total_segments),
                    &mut output.traced_segments,
                    &mut || stream.next_f64(),
                    &mut |cell, capture| spatial.capture(cell, lane, history, capture),
                )?;
                row.precision_limited_histories[lane] += u32::from(outcome.precision_limited);
                for (moment, value) in totals.iter_mut().zip([
                    outcome.values[2],
                    outcome.values[0] + outcome.values[1],
                    outcome.values[3],
                ]) {
                    moment.observe(value)?;
                }
            }
            for (estimate, moment) in [&mut row.absorbed, &mut row.escape, &mut row.unresolved]
                .into_iter()
                .zip(totals)
            {
                (estimate.mean[lane], estimate.standard_error[lane]) = moment.estimate()?;
            }
        }
        let receivers = spatial.finish(check, &[config.histories_per_emitter_lane; M])?;
        row.absorption
            .try_reserve_exact(receivers.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for receiver in receivers {
            check()?;
            row.absorption.push(LuminescenceAbsorption {
                cell_index: receiver.cell_index,
                fraction: FractionEstimate {
                    mean: receiver.absorbed.mean_w,
                    standard_error: receiver.absorbed.standard_error_w,
                },
            });
        }
        retained_links += row.absorption.len();
        output.rows.push(row);
    }
    check()?;
    Ok(output)
}

fn validate_config(config: LuminescenceTransportConfig) -> Result<(), VolumeTransportError> {
    if config.histories_per_emitter_lane < 2 || config.max_collisions_per_history == 0 {
        Err(VolumeTransportError::Budget)
    } else {
        Ok(())
    }
}

fn validate_wavelengths<const M: usize>(
    wavelengths: &[f64; M],
    (lo, hi): (f64, f64),
) -> Result<(), VolumeTransportError> {
    if M == 0 {
        return Err(VolumeTransportError::Sample);
    }
    if !lo.is_finite()
        || !hi.is_finite()
        || lo <= 0.0
        || lo >= hi
        || wavelengths
            .iter()
            .any(|value| !value.is_finite() || !(lo..=hi).contains(value))
    {
        return Err(VolumeTransportError::Wavelength);
    }
    Ok(())
}

fn validate_cell(
    cell: usize,
    previous: Option<usize>,
    count: usize,
) -> Result<(), VolumeTransportError> {
    if cell >= count || previous.is_some_and(|previous| previous >= cell) {
        Err(VolumeTransportError::Grid)
    } else {
        Ok(())
    }
}

fn minimum_segments<const M: usize>(
    cells: usize,
    config: LuminescenceTransportConfig,
) -> Result<u64, VolumeTransportError> {
    u64::try_from(cells)
        .ok()
        .and_then(|cells| cells.checked_mul(u64::try_from(M).ok()?))
        .and_then(|lanes| lanes.checked_mul(u64::from(config.histories_per_emitter_lane)))
        .ok_or(VolumeTransportError::Budget)
}

fn reserve_one<T>(values: &mut Vec<T>, limit: usize) -> Result<(), VolumeTransportError> {
    if values.len() == limit {
        return Err(VolumeTransportError::Budget);
    }
    if values.len() == values.capacity() {
        values
            .try_reserve_exact(values.capacity().max(16).min(limit - values.len()))
            .map_err(|_| VolumeTransportError::Budget)?;
    }
    Ok(())
}

fn emitter_seed(seed: u64, source: u64, cell: usize, lane: usize) -> u64 {
    let mut hash = DomainHasher::new("org.frankensim.render.luminescence-transfer-stream.v1");
    for value in [seed, source, cell as u64, lane as u64] {
        hash.update(&value.to_le_bytes());
    }
    u64::from_le_bytes(
        hash.finalize().as_bytes()[..8]
            .try_into()
            .expect("eight hash bytes"),
    )
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

#[cfg(test)]
mod tests;
