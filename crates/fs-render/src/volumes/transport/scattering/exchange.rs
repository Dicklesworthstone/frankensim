//! Reciprocal volume radiation exchange on one prescribed, index-one grid.
//!
//! For each wavelength, an isotropic unit packet is emitted at a uniform point
//! in each absorbing cell. The shared implicit-capture history gives directed
//! absorption probabilities p_ij and escape probabilities. With
//! A_i = 4π rho_i kappa_a_i V_i, reciprocity gives A_i p_ij = A_j p_ji.
//! The retained finite estimate C_ij = (A_i p_ij + A_j p_ji)/2 is nonnegative
//! and symmetric. Multiplying C_ij by a Planck-radiance difference and a positive
//! wavelength weight yields exactly equal/opposite heat transfers. Isothermal
//! equilibrium therefore has zero exchange even with finite photon samples.
//!
//! C_i,infinity = A_i p_i,escape couples a cell to an isotropic external bath.
//! No boundary photons or independent emitted-minus-absorbed subtraction are
//! needed. Self absorption cancels, scattering is elastic, and no dielectric
//! boundary, polarization or wavelength redistribution is inferred.
//!
//! This volume-launched reciprocal operator uses the net-exchange formulation
//! discussed by Eymet et al., JQSRT 91 (2005), doi:10.1016/j.jqsrt.2004.05.049;
//! their boundary-based acceleration for thick media is not implemented here.
//! Symmetrization does NOT normalize noisy rows:
//! exact finite-sample pair conservation differs from emissive-capacity closure
//! in expectation. Collision tails are retained explicitly. Empirical standard
//! errors are fixed-sample diagnostics, not confidence intervals, certified
//! omission bounds or guarantees that rare paths were observed.

use super::super::{AbsorbingGrid, VolumeTransportError, checkpoint};
use super::ScatteringGrid;
use super::deposition::{Moments, SpatialTallies, trace_history};
use super::field::CellScatteringGrid;
use super::prepared::{Monochromatic, PaletteOptics, PreparedOptics, UniformOptics};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;
use fs_rand::StreamKey;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RadiativeExchangeConfig {
    /// Independent uniform-position/isotropic-direction histories, at least 2,
    /// for EVERY cell/lane whose emission capacity is positive.
    pub histories_per_emitter_lane: u32,
    pub max_collisions_per_history: u32,
    /// All source cells and lanes together, not renewed per emitter.
    pub max_total_segments: u64,
    /// All positive-density material cells, including transparent/scattering
    /// cells whose temperature can invalidate the prescribed optical snapshot.
    pub max_material_cells: usize,
    pub max_emitter_cells: usize,
    /// Maximum distinct unordered interacting cell pairs retained.
    pub max_pairs: usize,
    pub seed: u64,
    pub source_id: u64,
}

/// Multiply by W/(m² sr nm) and quadrature nm to obtain watts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExchangeCoefficientEstimate<const N: usize> {
    pub mean_m2_sr: [f64; N],
    pub standard_error_m2_sr: [f64; N],
}
impl<const N: usize> ExchangeCoefficientEstimate<N> {
    fn zero() -> Self {
        Self {
            mean_m2_sr: [0.0; N],
            standard_error_m2_sr: [0.0; N],
        }
    }
    fn set(
        &mut self,
        lane: usize,
        fraction: (f64, f64),
        capacity: f64,
    ) -> Result<(), VolumeTransportError> {
        self.mean_m2_sr[lane] = scale(fraction.0, capacity)?;
        self.standard_error_m2_sr[lane] = scale(fraction.1, capacity)?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RadiativeExchangeCell<const N: usize> {
    pub cell_index: usize,
    /// Exact supplied-model A_i = 4π rho_i kappa_a_i V_i; no Planck factor.
    pub emission_capacity_m2_sr: [f64; N],
    /// Complete-history total directed absorption, INCLUDING self return.
    pub absorbed: ExchangeCoefficientEstimate<N>,
    pub self_return: ExchangeCoefficientEstimate<N>,
    pub escape: ExchangeCoefficientEstimate<N>,
    /// Surviving coefficient at a collision cap or below-resolution split.
    /// Omitted net heat depends on Planck differences, not emitted power alone.
    pub unresolved: ExchangeCoefficientEstimate<N>,
    pub sampled_histories: [u32; N],
    pub precision_limited_histories: [u32; N],
}

#[derive(Clone, Debug, PartialEq)]
pub struct RadiativeExchangePair<const N: usize> {
    pub cell_a: usize,
    pub cell_b: usize,
    /// One coefficient for BOTH directions, so positivity and reciprocity hold
    /// in the retained sample. Edges sharing an emitter row are correlated.
    pub exchange: ExchangeCoefficientEstimate<N>,
    /// Directed a→b minus b→a before symmetrization. Its standard error is
    /// twice `exchange.standard_error_m2_sr`, from independent emitter streams.
    pub reciprocity_residual_m2_sr: [f64; N],
}

#[derive(Clone, Debug, PartialEq)]
pub struct RadiativeExchangeOperator<const N: usize> {
    /// Binds the original density, coefficients, phases and source conditions.
    pub medium_identity: ContentHash,
    pub wavelength_bounds_nm: (f64, f64),
    pub wavelengths_nm: [f64; N],
    pub grid_cell_count: usize,
    pub cell_volume_m3: f64,
    /// Ascending original positive-density cells. These all require the fixed
    /// optical model's temperature-validity binding, even when they emit no light.
    pub material_cells: Vec<usize>,
    /// Ascending ORIGINAL grid indices with positive capacity in any lane.
    pub cells: Vec<RadiativeExchangeCell<N>>,
    /// Ascending (cell_a, cell_b), with cell_a < cell_b.
    pub pairs: Vec<RadiativeExchangePair<N>>,
    pub traced_segments: u64,
    pub config: RadiativeExchangeConfig,
}

impl<const N: usize> RadiativeExchangeOperator<N> {
    /// Admission policy for consumers of a retained/public report. Directed
    /// rows are reconstructed BEFORE any thermal use; symmetric row sums are
    /// intentionally not forced to equal emitter capacities. The supplied
    /// coefficient roundoff allowance is absolute m² sr and is not enlarged.
    pub fn validate_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        coefficient_roundoff_m2_sr: f64,
    ) -> Result<(), VolumeTransportError> {
        check()?;
        let allowance = coefficient_roundoff_m2_sr;
        let (lo, hi) = self.wavelength_bounds_nm;
        if N == 0 {
            return Err(VolumeTransportError::Sample);
        }
        if !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || lo >= hi
            || self
                .wavelengths_nm
                .iter()
                .any(|v| !v.is_finite() || !(lo..=hi).contains(v))
        {
            return Err(VolumeTransportError::Wavelength);
        }
        if !nonnegative(allowance)
            || !self.cell_volume_m3.is_finite()
            || self.cell_volume_m3 <= 0.0
            || self.grid_cell_count == 0
        {
            return Err(VolumeTransportError::Unrepresentable);
        }
        if self.config.histories_per_emitter_lane < 2
            || self.config.max_collisions_per_history == 0
            || self.material_cells.len() > self.config.max_material_cells
            || self.cells.len() > self.config.max_emitter_cells
            || self.pairs.len() > self.config.max_pairs
            || self.traced_segments > self.config.max_total_segments
        {
            return Err(VolumeTransportError::Budget);
        }
        let mut previous_material = None;
        for &cell in &self.material_cells {
            check()?;
            if cell >= self.grid_cell_count
                || previous_material.is_some_and(|previous| previous >= cell)
            {
                return Err(VolumeTransportError::Grid);
            }
            previous_material = Some(cell);
        }
        let mut sums = Vec::new();
        sums.try_reserve_exact(self.cells.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        let mut previous = None;
        let mut minimum_segments = 0_u64;
        for cell in &self.cells {
            check()?;
            if cell.cell_index >= self.grid_cell_count
                || self.material_cells.binary_search(&cell.cell_index).is_err()
                || previous.is_some_and(|id| id >= cell.cell_index)
                || !cell.emission_capacity_m2_sr.iter().all(|a| nonnegative(*a))
                || cell.emission_capacity_m2_sr.iter().all(|a| *a == 0.0)
            {
                return Err(VolumeTransportError::Unrepresentable);
            }
            previous = Some(cell.cell_index);
            for estimate in [
                &cell.absorbed,
                &cell.self_return,
                &cell.escape,
                &cell.unresolved,
            ] {
                if !estimate
                    .mean_m2_sr
                    .iter()
                    .chain(estimate.standard_error_m2_sr.iter())
                    .all(|v| nonnegative(*v))
                {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
            let mut row = [CoefficientSum::default(); N];
            for lane in 0..N {
                let capacity = cell.emission_capacity_m2_sr[lane];
                let histories = if capacity > 0.0 {
                    self.config.histories_per_emitter_lane
                } else {
                    0
                };
                if cell.sampled_histories[lane] != histories
                    || cell.precision_limited_histories[lane] > histories
                {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                minimum_segments = minimum_segments
                    .checked_add(u64::from(histories))
                    .ok_or(VolumeTransportError::Budget)?;
                for estimate in [
                    &cell.absorbed,
                    &cell.self_return,
                    &cell.escape,
                    &cell.unresolved,
                ] {
                    if estimate.mean_m2_sr[lane] - capacity > allowance
                        || (capacity == 0.0
                            && (estimate.mean_m2_sr[lane] != 0.0
                                || estimate.standard_error_m2_sr[lane] != 0.0))
                    {
                        return Err(VolumeTransportError::Unrepresentable);
                    }
                }
                if cell.self_return.mean_m2_sr[lane] - cell.absorbed.mean_m2_sr[lane] > allowance {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                let mut balance = CoefficientSum::default();
                for value in [
                    cell.absorbed.mean_m2_sr[lane],
                    cell.escape.mean_m2_sr[lane],
                    cell.unresolved.mean_m2_sr[lane],
                    -capacity,
                ] {
                    balance.add(value)?;
                }
                if balance.value().abs() > allowance {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                row[lane].add(cell.self_return.mean_m2_sr[lane])?;
            }
            sums.push(row);
        }
        if self.traced_segments < minimum_segments
            || self.traced_segments
                > minimum_segments.saturating_mul(u64::from(self.config.max_collisions_per_history))
            || (self.cells.is_empty() && self.traced_segments != 0)
        {
            return Err(VolumeTransportError::Budget);
        }
        let mut previous = None;
        for pair in &self.pairs {
            check()?;
            let key = (pair.cell_a, pair.cell_b);
            if pair.cell_a >= pair.cell_b
                || previous.is_some_and(|k| k >= key)
                || !pair
                    .exchange
                    .mean_m2_sr
                    .iter()
                    .chain(pair.exchange.standard_error_m2_sr.iter())
                    .all(|v| nonnegative(*v))
                || !pair
                    .reciprocity_residual_m2_sr
                    .iter()
                    .all(|v| v.is_finite())
            {
                return Err(VolumeTransportError::Unrepresentable);
            }
            previous = Some(key);
            let a = self
                .cells
                .binary_search_by_key(&pair.cell_a, |c| c.cell_index)
                .map_err(|_| VolumeTransportError::Grid)?;
            let b = self
                .cells
                .binary_search_by_key(&pair.cell_b, |c| c.cell_index)
                .map_err(|_| VolumeTransportError::Grid)?;
            for lane in 0..N {
                let half_residual = 0.5 * pair.reciprocity_residual_m2_sr[lane];
                let forward = pair.exchange.mean_m2_sr[lane] + half_residual;
                let reverse = pair.exchange.mean_m2_sr[lane] - half_residual;
                if !forward.is_finite()
                    || !reverse.is_finite()
                    || forward < -allowance
                    || reverse < -allowance
                    || ((self.cells[a].emission_capacity_m2_sr[lane] == 0.0
                        || self.cells[b].emission_capacity_m2_sr[lane] == 0.0)
                        && (pair.exchange.mean_m2_sr[lane] != 0.0
                            || pair.exchange.standard_error_m2_sr[lane] != 0.0
                            || pair.reciprocity_residual_m2_sr[lane] != 0.0))
                {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                sums[a][lane].add(forward)?;
                sums[b][lane].add(reverse)?;
            }
        }
        for (cell, sum) in self.cells.iter().zip(sums) {
            check()?;
            for lane in 0..N {
                if (sum[lane].value() - cell.absorbed.mean_m2_sr[lane]).abs() > allowance {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
        }
        check()
    }
}

#[derive(Clone, Copy, Default)]
struct CoefficientSum {
    sum: f64,
    correction: f64,
}
impl CoefficientSum {
    fn add(&mut self, value: f64) -> Result<(), VolumeTransportError> {
        let next = self.sum + value;
        if !next.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        if !self.correction.is_finite() || !self.value().is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(())
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

impl ScatteringGrid<'_> {
    /// Estimate temperature-independent exchange coefficients for a fixed
    /// optical snapshot. Source and receiver coefficients are the SAME local
    /// absorption owner used by rendering's LTE source and by incident heating.
    pub fn estimate_thermal_exchange<const N: usize>(
        &self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64; N],
        config: RadiativeExchangeConfig,
    ) -> Result<RadiativeExchangeOperator<N>, VolumeTransportError> {
        checkpoint(cx)?;
        if self.absorption.absorption_palette().is_some() || self.scattering_palette().is_some() {
            let prepared = PaletteOptics::try_new(self, *wavelengths_nm, &mut || checkpoint(cx))?;
            return estimate(&prepared, &mut || checkpoint(cx), config);
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
            config,
        )
    }
}

impl<const N: usize> CellScatteringGrid<'_, '_, N> {
    pub fn estimate_thermal_exchange(
        &self,
        cx: &Cx<'_>,
        config: RadiativeExchangeConfig,
    ) -> Result<RadiativeExchangeOperator<N>, VolumeTransportError> {
        self.estimate_thermal_exchange_with_checkpoint(&mut || checkpoint(cx), config)
    }
    pub fn estimate_thermal_exchange_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        config: RadiativeExchangeConfig,
    ) -> Result<RadiativeExchangeOperator<N>, VolumeTransportError> {
        estimate(self, check, config)
    }
}

fn estimate<const N: usize>(
    medium: &impl PreparedOptics<N>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    config: RadiativeExchangeConfig,
) -> Result<RadiativeExchangeOperator<N>, VolumeTransportError> {
    check()?;
    if N == 0 {
        return Err(VolumeTransportError::Sample);
    }
    if config.histories_per_emitter_lane < 2 || config.max_collisions_per_history == 0 {
        return Err(VolumeTransportError::Budget);
    }
    let geometry = medium.geometry();
    let cell_volume_m3 = positive_product(&geometry.grid.cell)?;
    let mut output = RadiativeExchangeOperator {
        medium_identity: medium.identity(),
        wavelength_bounds_nm: medium.bounds(),
        wavelengths_nm: *medium.wavelengths(),
        grid_cell_count: geometry.cell_count(),
        cell_volume_m3,
        material_cells: Vec::new(),
        cells: Vec::new(),
        pairs: Vec::new(),
        traced_segments: 0,
        config,
    };
    let mut minimum_segments = 0_u64;
    for cell in 0..geometry.cell_count() {
        check()?;
        let density = geometry.grid.data[cell];
        if density == 0.0 {
            continue;
        }
        reserve(&mut output.material_cells, config.max_material_cells)?;
        output.material_cells.push(cell);
        let mut capacity = [0.0; N];
        for lane in 0..N {
            let absorption = medium.coefficient(cell, lane)?.0;
            if absorption > 0.0 {
                capacity[lane] = positive_product(&[
                    4.0 * core::f64::consts::PI,
                    density,
                    absorption,
                    cell_volume_m3,
                ])?;
                minimum_segments = minimum_segments
                    .checked_add(u64::from(config.histories_per_emitter_lane))
                    .ok_or(VolumeTransportError::Budget)?;
            }
        }
        if capacity.iter().all(|a| *a == 0.0) {
            continue;
        }
        reserve(&mut output.cells, config.max_emitter_cells)?;
        output.cells.push(RadiativeExchangeCell {
            cell_index: cell,
            emission_capacity_m2_sr: capacity,
            absorbed: ExchangeCoefficientEstimate::zero(),
            self_return: ExchangeCoefficientEstimate::zero(),
            escape: ExchangeCoefficientEstimate::zero(),
            unresolved: ExchangeCoefficientEstimate::zero(),
            sampled_histories: [0; N],
            precision_limited_histories: [0; N],
        });
    }
    if minimum_segments > config.max_total_segments {
        return Err(VolumeTransportError::Budget);
    }
    let mut pairs = PairTallies::<N>::new(config.max_pairs);
    for row in &mut output.cells {
        check()?;
        let mut spatial = SpatialTallies::<N>::new(config.max_emitter_cells);
        for lane in 0..N {
            let capacity = row.emission_capacity_m2_sr[lane];
            if capacity == 0.0 {
                continue;
            }
            row.sampled_histories[lane] = config.histories_per_emitter_lane;
            let mono = Monochromatic::<_, N> {
                medium,
                lane,
                wavelength: [output.wavelengths_nm[lane]],
            };
            let mut totals = [Moments::default(); 3]; // absorption, escape, unresolved
            let seed = emitter_seed(config.seed, config.source_id, row.cell_index, lane);
            for history in 0..config.histories_per_emitter_lane {
                check()?;
                let mut stream = StreamKey {
                    seed,
                    kernel: 0x4E45_5843,
                    tile: history,
                }
                .stream(); // NEXC
                let (position, direction) =
                    emitted_ray(geometry, row.cell_index, &mut || stream.next_f64())?;
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
                estimate.set(lane, moment.estimate()?, capacity)?;
            }
        }
        for receiver in spatial.finish(check, &row.sampled_histories)? {
            check()?;
            let mut directed = ExchangeCoefficientEstimate::zero();
            for lane in 0..N {
                directed.set(
                    lane,
                    (
                        receiver.absorbed.mean_w[lane],
                        receiver.absorbed.standard_error_w[lane],
                    ),
                    row.emission_capacity_m2_sr[lane],
                )?;
            }
            if receiver.cell_index == row.cell_index {
                row.self_return = directed;
            } else {
                pairs.insert(row.cell_index, receiver.cell_index, directed)?;
            }
        }
    }
    output.pairs = pairs.finish(check)?;
    check()?;
    Ok(output)
}

/// No rejection resampling or shifted origin: an unrepresentable source point
/// is a named refusal rather than being moved into another physical cell.
fn emitted_ray(
    geometry: &AbsorbingGrid<'_>,
    cell: usize,
    random: &mut impl FnMut() -> f64,
) -> Result<(Point3, Vec3), VolumeTransportError> {
    let grid = &geometry.grid;
    let lattice = [
        cell % grid.dims[0],
        (cell / grid.dims[0]) % grid.dims[1],
        cell / (grid.dims[0] * grid.dims[1]),
    ];
    let mut coordinates = [0.0; 3];
    for axis in 0..3 {
        let lo = (lattice[axis] as f64).mul_add(grid.cell[axis], grid.origin[axis]);
        let hi = ((lattice[axis] + 1) as f64).mul_add(grid.cell[axis], grid.origin[axis]);
        coordinates[axis] = random().mul_add(hi - lo, lo);
        if !(coordinates[axis] > lo && coordinates[axis] < hi) {
            return Err(VolumeTransportError::Unrepresentable);
        }
    }
    let z = 1.0 - 2.0 * random();
    let radius = det::sqrt(((1.0 - z) * (1.0 + z)).max(0.0));
    let (sin, cos) = det::sin_cos(core::f64::consts::TAU * random());
    Ok((
        Point3::new(coordinates[0], coordinates[1], coordinates[2]),
        Vec3::new(radius * cos, radius * sin, z),
    ))
}

fn emitter_seed(seed: u64, source: u64, cell: usize, lane: usize) -> u64 {
    let mut hash = DomainHasher::new("org.frankensim.render.reciprocal-exchange-stream.v1");
    for value in [seed, source, cell as u64, lane as u64] {
        hash.update(&value.to_le_bytes());
    }
    u64::from_le_bytes(
        hash.finalize().as_bytes()[..8]
            .try_into()
            .expect("eight hash bytes"),
    )
}

struct DirectedPair<const N: usize> {
    cells: (usize, usize),
    forward: ExchangeCoefficientEstimate<N>,
    reverse: ExchangeCoefficientEstimate<N>,
}
struct PairTallies<const N: usize> {
    pairs: Vec<DirectedPair<N>>,
    indices: HashMap<(usize, usize), usize>,
    limit: usize,
}
impl<const N: usize> PairTallies<N> {
    fn new(limit: usize) -> Self {
        Self {
            pairs: Vec::new(),
            indices: HashMap::new(),
            limit,
        }
    }
    fn insert(
        &mut self,
        emitter: usize,
        receiver: usize,
        value: ExchangeCoefficientEstimate<N>,
    ) -> Result<(), VolumeTransportError> {
        let cells = (emitter.min(receiver), emitter.max(receiver));
        let slot = if let Some(&slot) = self.indices.get(&cells) {
            slot
        } else {
            reserve(&mut self.pairs, self.limit)?;
            self.indices
                .try_reserve(1)
                .map_err(|_| VolumeTransportError::Budget)?;
            let slot = self.pairs.len();
            self.pairs.push(DirectedPair {
                cells,
                forward: ExchangeCoefficientEstimate::zero(),
                reverse: ExchangeCoefficientEstimate::zero(),
            });
            self.indices.insert(cells, slot);
            slot
        };
        if emitter == cells.0 {
            self.pairs[slot].forward = value;
        } else {
            self.pairs[slot].reverse = value;
        }
        Ok(())
    }
    fn finish(
        mut self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    ) -> Result<Vec<RadiativeExchangePair<N>>, VolumeTransportError> {
        check()?;
        self.pairs.sort_unstable_by_key(|p| p.cells);
        let mut result = Vec::new();
        result
            .try_reserve_exact(self.pairs.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for pair in self.pairs {
            check()?;
            let mut exchange = ExchangeCoefficientEstimate::zero();
            let mut residual = [0.0; N];
            for lane in 0..N {
                let a = pair.forward.mean_m2_sr[lane];
                let b = pair.reverse.mean_m2_sr[lane];
                let sum = a + b;
                exchange.mean_m2_sr[lane] = if sum.is_finite() {
                    0.5 * sum
                } else {
                    0.5 * a + 0.5 * b
                };
                if (a > 0.0 || b > 0.0) && exchange.mean_m2_sr[lane] == 0.0 {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                residual[lane] = a - b;
                let sa = pair.forward.standard_error_m2_sr[lane];
                let sb = pair.reverse.standard_error_m2_sr[lane];
                let maximum = sa.max(sb);
                exchange.standard_error_m2_sr[lane] = if maximum == 0.0 {
                    0.0
                } else {
                    maximum * (0.5 * det::sqrt((sa / maximum).powi(2) + (sb / maximum).powi(2)))
                };
            }
            for lane in 0..N {
                if (pair.forward.standard_error_m2_sr[lane] > 0.0
                    || pair.reverse.standard_error_m2_sr[lane] > 0.0)
                    && exchange.standard_error_m2_sr[lane] == 0.0
                {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
            result.push(RadiativeExchangePair {
                cell_a: pair.cells.0,
                cell_b: pair.cells.1,
                exchange,
                reciprocity_residual_m2_sr: residual,
            });
        }
        Ok(result)
    }
}

fn reserve<T>(values: &mut Vec<T>, limit: usize) -> Result<(), VolumeTransportError> {
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
fn positive_product(factors: &[f64]) -> Result<f64, VolumeTransportError> {
    let direct: f64 = factors.iter().product();
    let result = if direct.is_finite() && direct >= f64::MIN_POSITIVE {
        direct
    } else {
        det::exp(factors.iter().map(|v| det::ln(*v)).sum())
    };
    if result.is_finite() && result > 0.0 {
        Ok(result)
    } else {
        Err(VolumeTransportError::Unrepresentable)
    }
}
fn scale(value: f64, capacity: f64) -> Result<f64, VolumeTransportError> {
    if value == 0.0 || capacity == 0.0 {
        return Ok(0.0);
    }
    let result = value * capacity;
    if result.is_finite() && result > 0.0 {
        Ok(result)
    } else {
        Err(VolumeTransportError::Unrepresentable)
    }
}

#[cfg(test)]
mod tests;
