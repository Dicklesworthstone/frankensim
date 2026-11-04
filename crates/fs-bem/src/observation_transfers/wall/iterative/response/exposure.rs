//! Finite-window work from the original retained boundary fields.
//!
//! The wall's frequency-dependent admittance supplies its actual inward
//! velocity. Finite-window pressure work includes reversible wall storage;
//! it must not be mislabeled irreversible heat or cycle-averaged absorption.

use super::*;
use crate::observation_transfers::HarmonicExposureFactorStatistics;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::dd::Dd;

/// One original frequency response and its physical source-rate coordinates.
/// All lines use the same exp(-i omega t) emission clock. Propagation phase is
/// already in the retained pressure/velocity solutions.
#[derive(Clone, Copy)]
pub struct BoundaryWorkExposureLine<'a, 'o, 's> {
    pub response: &'a IterativeBoundaryResponse<'o, 's>,
    pub rates: &'a [C64],
}

/// One zero-mean coordinate shared by all its sources and frequencies.
/// Distinct factors have zero conjugate and unconjugated cross moments, as in
/// the original harmonic exposure owner. Their statistics are authored, never
/// inferred from the number of lines or the boundary material.
#[derive(Clone, Copy)]
pub struct BoundaryWorkExposureFactor<'a, 'o, 's> {
    pub lines: &'a [BoundaryWorkExposureLine<'a, 'o, 's>],
    pub statistics: HarmonicExposureFactorStatistics,
}

#[derive(Clone, Copy, Debug)]
pub struct BoundaryWorkExposureLimits {
    /// Includes all supplied lines, before frequency grouping; hard cap 4096.
    pub maximum_lines: usize,
    /// Includes empty factors; hard cap 256.
    pub maximum_factors: usize,
    /// Three coherent trace products per source, panel and supplied line.
    pub maximum_source_panel_products: usize,
    /// Panels times sum of squared supplied line counts over all banks.
    pub maximum_panel_frequency_pairs: usize,
    /// Retained DD fields/headers, sort indices and material query buffers.
    /// Borrowed inputs, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Expected normalized-window boundary work [W], not integrated energy [J].
/// The original physical finite-source map and discrete boundary orientation
/// apply. Signed values are preserved, including release of stored energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowedBoundaryWork {
    /// Integral of p * prescribed outward normal velocity over the surface.
    pub input_power_w: f64,
    /// Integral of p * actual outward normal velocity from the retained solve.
    pub outward_power_w: f64,
    /// Work into the local wall: p * (Y_wall p), including its reactive storage.
    /// Only a complete cycle/stationary average yields absorbed wall power.
    pub wall_power_w: f64,
    /// Input minus outward minus wall work; measured, never forced to zero.
    pub balance_w: f64,
}

// p, actual outward velocity, prescribed velocity, inward wall velocity.
// Keep their complex components in DD through frequency grouping and products.
type Field = [[Dd; 2]; 4];
struct Line {
    omega: f64,
    fields: Vec<Field>,
}
struct Bank {
    lines: Vec<Line>,
    statistics: HarmonicExposureFactorStatistics,
}

/// Reusable work observer for one original finite-source boundary owner.
/// Preparation contracts retained source fields once; every later shutter
/// query allocates nothing and performs no BEM solve or material query.
///
/// Input and outward work use distinct retained traces. Wall work uses the
/// original local law directly, preserving weak wall exchange beside large
/// input/outward powers. Their residual therefore still exposes a boundary
/// solve defect. Far-field radiation, volume storage and irreversible wall
/// heat are not inferred from this finite-surface observation.
pub struct PreparedBoundaryWorkExposure<'a> {
    source_clock: &'a str,
    sources: &'a [SourceId],
    areas: &'a [f64],
    operator_hash: ContentHash,
    mean: Bank,
    factors: Vec<Bank>,
    workspace_bytes: usize,
    panel_frequency_pairs: usize,
}

impl<'a> PreparedBoundaryWorkExposure<'a> {
    /// Bind frequencies from the SAME immutable iterative boundary owner.
    /// A matching hash alone is insufficient: source geometry, local material
    /// laws and numerical settings must share the original owner. Point-source
    /// banks have no finite self-work surface and refuse, even when dark.
    pub fn new(
        reference: &'a IterativeBoundaryResponse<'_, '_>,
        source_clock: &'a str,
        mean: &[BoundaryWorkExposureLine<'_, '_, '_>],
        factors: &[BoundaryWorkExposureFactor<'_, '_, '_>],
        limits: BoundaryWorkExposureLimits,
        gate: &CancelGate,
    ) -> Result<Self, IterativeObservationError> {
        poll(gate)?;
        let Excitation::Surface { owner, .. } = &reference.observer.observation.excitation else {
            return Err(invalid("boundary exposure requires finite surface sources"));
        };
        if source_clock.trim().is_empty() {
            return Err(invalid(
                "boundary exposure requires an explicit source clock",
            ));
        }
        let panels = reference.surface().areas().len();
        let sources = reference.sources().len();
        let (mut lines, mut pairs, mut largest) = (0usize, 0usize, 0usize);
        admit(factors.len(), limits.maximum_factors.min(256))?;
        for bank in core::iter::once(mean).chain(factors.iter().map(|f| f.lines)) {
            lines = lines.checked_add(bank.len()).ok_or_else(overflow)?;
            pairs = bank
                .len()
                .checked_mul(bank.len())
                .and_then(|n| pairs.checked_add(n))
                .ok_or_else(overflow)?;
            largest = largest.max(bank.len());
        }
        admit(lines, limits.maximum_lines.min(4096))?;
        let source_products = panels
            .checked_mul(sources)
            .and_then(|n| n.checked_mul(lines))
            .and_then(|n| n.checked_mul(3))
            .ok_or_else(overflow)?;
        let panel_frequency_pairs = panels.checked_mul(pairs).ok_or_else(overflow)?;
        admit(source_products, limits.maximum_source_panel_products)?;
        admit(panel_frequency_pairs, limits.maximum_panel_frequency_pairs)?;
        let per_line = panels
            .checked_mul(core::mem::size_of::<Field>())
            .and_then(|n| n.checked_add(core::mem::size_of::<Line>()))
            .ok_or_else(overflow)?;
        let workspace_bytes = lines.checked_mul(per_line)
            .and_then(|n| n.checked_add(factors.len().checked_mul(core::mem::size_of::<Bank>())?))
            .and_then(|n| n.checked_add(largest.checked_mul(core::mem::size_of::<usize>())?))
            // The local law owner retains both law and panel coefficients;
            // every law is used, so two panel-sized buffers bound both.
            .and_then(|n| n.checked_add(panels.checked_mul(2 * core::mem::size_of::<C64>())?))
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;
        // Check the complete supplied bank before any field/material allocation.
        for bank in core::iter::once(mean).chain(factors.iter().map(|f| f.lines)) {
            for line in bank {
                poll(gate)?;
                if !core::ptr::eq(line.response.observer, reference.observer)
                    || line.rates.len() != sources
                    || line.rates.iter().any(|&rate| !finite(rate))
                {
                    return Err(invalid(
                        "boundary exposure source/owner mismatch or nonfinite rate",
                    ));
                }
            }
        }
        let mean = prepare_bank(mean, HarmonicExposureFactorStatistics::Real, owner, gate)?;
        let mut prepared_factors = storage(factors.len())?;
        for factor in factors {
            prepared_factors.push(prepare_bank(factor.lines, factor.statistics, owner, gate)?);
        }
        let actual_pairs = core::iter::once(&mean)
            .chain(prepared_factors.iter())
            .map(|bank| bank.lines.len() * bank.lines.len())
            .sum::<usize>();
        poll(gate)?;
        Ok(Self {
            source_clock,
            sources: reference.sources(),
            areas: reference.surface().areas(),
            operator_hash: reference.operator_hash(),
            mean,
            factors: prepared_factors,
            workspace_bytes,
            panel_frequency_pairs: panels * actual_pairs,
        })
    }

    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn sources(&self) -> &[SourceId] {
        self.sources
    }
    pub fn operator_hash(&self) -> ContentHash {
        self.operator_hash
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn panel_frequency_pairs(&self) -> usize {
        self.panel_frequency_pairs
    }

    /// Exact shared harmonic-window contraction of deterministic/real/proper
    /// banks, including carrier-sum and cross-frequency beat terms where their
    /// authored statistics require them. No clamp or implicit cycle average.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
        gate: &CancelGate,
    ) -> Result<WindowedBoundaryWork, IterativeObservationError> {
        poll(gate)?;
        let mut sums = [Dd::ZERO; 3];
        for bank in core::iter::once(&self.mean).chain(self.factors.iter()) {
            for a in &bank.lines {
                for b in &bank.lines {
                    poll(gate)?;
                    let kernel = bank.statistics.kernel(window, a.omega, b.omega)?;
                    for (panel, &area) in self.areas.iter().enumerate() {
                        if panel % 64 == 0 {
                            poll(gate)?;
                        }
                        for (sum, component) in sums.iter_mut().zip([2, 1, 3]) {
                            *sum = *sum
                                + kernel.product(a.fields[panel][0], b.fields[panel][component])?
                                    * Dd::from_f64(area);
                        }
                    }
                }
            }
        }
        let balance = sums[0] - sums[1] - sums[2];
        if sums
            .into_iter()
            .chain([balance])
            .any(|value| !finite_dd(value))
        {
            return Err(invalid("boundary exposure work is unrepresentable"));
        }
        poll(gate)?;
        Ok(WindowedBoundaryWork {
            input_power_w: sums[0].to_f64(),
            outward_power_w: sums[1].to_f64(),
            wall_power_w: sums[2].to_f64(),
            balance_w: balance.to_f64(),
        })
    }
}

fn prepare_bank(
    input: &[BoundaryWorkExposureLine<'_, '_, '_>],
    statistics: HarmonicExposureFactorStatistics,
    owner: &BemObservation<'_>,
    gate: &CancelGate,
) -> Result<Bank, IterativeObservationError> {
    let mut indices = storage(input.len())?;
    indices.extend(0..input.len());
    indices.sort_unstable_by(|&a, &b| {
        input[a]
            .response
            .omega
            .total_cmp(&input[b].response.omega)
            .then_with(|| a.cmp(&b))
    });
    let mut lines: Vec<Line> = storage(input.len())?;
    for index in indices {
        poll(gate)?;
        let line = input[index];
        if lines
            .last()
            .is_none_or(|last| last.omega.to_bits() != line.response.omega.to_bits())
        {
            let mut fields = storage(owner.surface.areas().len())?;
            fields.resize(owner.surface.areas().len(), [[Dd::ZERO; 2]; 4]);
            lines.push(Line {
                omega: line.response.omega,
                fields,
            });
        }
        let fields = &mut lines.last_mut().expect("frequency group exists").fields;
        for (panel, field) in fields.iter_mut().enumerate() {
            poll(gate)?;
            for (source, rate) in line.rates.iter().enumerate() {
                if source % 32 == 0 {
                    poll(gate)?;
                }
                let solution = &line.response.solutions[source];
                for (target, value) in field.iter_mut().take(3).zip([
                    solution.pressure[panel],
                    solution.velocity[panel],
                    owner.fields[source][panel],
                ]) {
                    add_product(target, value, *rate);
                }
            }
        }
    }
    // Evaluate each bank's distinct frequencies once, after coherent grouping.
    // This original local law is not replaced by drive-minus-solved-velocity.
    if let Some(reference) = input.first() {
        let observation = &reference.response.observer.observation;
        let rho_c = reference.response.medium().density * reference.response.medium().sound_speed;
        for line in &mut lines {
            let beta = observation
                .boundary
                .normalized_admittance(line.omega, || gate.is_requested())?;
            for (panel, field) in line.fields.iter_mut().enumerate() {
                if panel % 64 == 0 {
                    poll(gate)?;
                }
                let br = Dd::from_f64(beta[panel].re) / Dd::from_f64(rho_c);
                let bi = Dd::from_f64(beta[panel].im) / Dd::from_f64(rho_c);
                field[3] = [
                    field[0][0] * br - field[0][1] * bi,
                    field[0][0] * bi + field[0][1] * br,
                ];
            }
        }
    }
    if lines
        .iter()
        .flat_map(|line| &line.fields)
        .flatten()
        .flatten()
        .any(|&x| !finite_dd(x))
    {
        return Err(invalid(
            "boundary exposure coherent field is unrepresentable",
        ));
    }
    Ok(Bank { lines, statistics })
}

fn add_product(target: &mut [Dd; 2], a: C64, b: C64) {
    let [ar, ai, br, bi] = [a.re, a.im, b.re, b.im].map(Dd::from_f64);
    target[0] = target[0] + ar * br - ai * bi;
    target[1] = target[1] + ar * bi + ai * br;
}
fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite()
}
fn invalid(message: &'static str) -> IterativeObservationError {
    TransferError::Invalid(message).into()
}
fn overflow() -> IterativeObservationError {
    invalid("boundary exposure work/storage overflow")
}
fn admit(requested: usize, cap: usize) -> Result<(), IterativeObservationError> {
    if requested > cap {
        Err(TransferError::Cap { requested, cap }.into())
    } else {
        Ok(())
    }
}
fn storage<T>(count: usize) -> Result<Vec<T>, IterativeObservationError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("boundary exposure allocation failed"))?;
    Ok(values)
}

#[cfg(test)]
mod tests;
