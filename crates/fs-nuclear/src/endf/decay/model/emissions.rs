//! Energy-resolved emitted particles from the original integrated decay counts.
//!
//! ENDF-102 (2023), §8.4, Eq. (8.1) gives the parent-wide spectral intensity
//! `FC * RP(E) + FD * sum(RI * delta(E-ER))`. Mode branching is already included
//! in those normalizations and is never applied again here. Prepared group
//! integrals are reusable across inventory histories and interval durations.
//! Histogram and linear continua are integrated analytically, including their
//! first energy moment; an energy-group midpoint is never substituted for it.
//!
//! Discrete beta records describe endpoints, not monoenergetic particles. They
//! require a separate beta-shape model and refuse here; a supplied beta-minus
//! continuum can be used when it has no positive endpoint contribution.
//! Electron-capture/positron continua also refuse: ENDF §8.4.2 rule 5(c) requires
//! a resolved positron-intensity normalization that is not supplied by the
//! retained per-line RIS fields alone. Fission-fragment and
//! unidentified radiation semantics are not inferred. Missing spectral records
//! remain explicit missing coverage, even when the aggregate heat data exist.
//!
//! Emitted particles and their energy are source observables, not extra local
//! heat. The original average-energy routing, thermal transaction, and material
//! remain unchanged. No transport, dose, covariance propagation, or independent
//! nuclear-data validation is claimed.
//!
//! Per-decay coefficients are retained as f64. A nonzero subnormal coefficient
//! cannot be amplified by a transformation count above one: rounding in that
//! coefficient could otherwise become a large relative error in a normal-range
//! output. Such applications refuse rather than implying preserved precision.

use std::{collections::BTreeMap, sync::Arc};

use super::EvaluatedDecayAssembly;
use crate::burnup::transport::BurnupTransportStep;
use crate::burnup::{BurnupModel, BurnupRate, BurnupStep};
use crate::decay::{DecayModel, DecayStep};
use crate::endf::decay::spectra::{DecayRadiation, DecaySpectrum, DecaySpectrumPoint};
use crate::endf::spectrum::PositiveIntegral;
use crate::endf::{EndfError, budget, convex};
use crate::irradiation::NuclideId;
use crate::{Accumulator, FissionError, nonnegative};

/// Explicit allocation and work bounds for preparation and source evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DecayEmissionLimits {
    pub max_groups: usize,
    /// Count and energy scalar values together, for prepared coefficients or
    /// the group arrays returned by a compartment-report evaluation.
    pub max_group_values: usize,
    /// Bounded source records, segment/group intersections, and coefficient
    /// applications. A compartment call shares one work allowance.
    pub max_work: usize,
}

impl Default for DecayEmissionLimits {
    fn default() -> Self {
        Self {
            max_groups: 4096,
            max_group_values: 524_288,
            max_work: 4_000_000,
        }
    }
}

/// Coverage for one original-order nuclide and the selected radiation kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DecayEmissionCoverage {
    /// The explicitly stable law has no transformations.
    Stable,
    /// The evaluation contains the selected radiation spectrum.
    Reported,
    /// Absence of a record does not prove absence of emitted particles.
    Unreported,
}

/// Actual transformations for which the selected spectrum was not supplied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnreportedDecayEmission {
    pub nuclide: NuclideId,
    pub decays: f64,
}

/// Integrated expected emitted particles and their energy, in original group
/// order. Internal edges belong to the upper bin; the last upper edge is closed.
/// Values include only reported spectra. Inspect `unreported_decays` before
/// interpreting the result as a complete source.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayEmissionStep {
    pub duration_s: f64,
    pub radiation: DecayRadiation,
    pub group_edges_j: Arc<[f64]>,
    pub emitted_particles: Vec<f64>,
    pub emitted_energy_j: Vec<f64>,
    pub unreported_decays: Vec<UnreportedDecayEmission>,
}

/// Immutable per-parent response. Clone shares the evaluation, group grid and
/// coefficient arrays; it never shares mutable inventory or emission history.
#[derive(Clone, Debug)]
pub struct PreparedDecayEmission {
    assembly: Arc<EvaluatedDecayAssembly>,
    radiation: DecayRadiation,
    edges: Arc<[f64]>,
    coverage: Arc<[DecayEmissionCoverage]>,
    particles: Arc<[f64]>,
    energy: Arc<[f64]>,
    limits: DecayEmissionLimits,
}

impl PreparedDecayEmission {
    /// Admit one selected radiation and a strictly increasing SI energy grid.
    /// Every positive retained line/continuum segment must fit in the grid;
    /// neither omitted tails nor extrapolated spectral values are invented.
    pub fn new(
        assembly: Arc<EvaluatedDecayAssembly>,
        radiation: DecayRadiation,
        group_edges_j: Vec<f64>,
        limits: DecayEmissionLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        let groups = group_edges_j
            .len()
            .checked_sub(1)
            .filter(|&n| n > 0)
            .ok_or_else(|| invalid("at least two emission group edges are required"))?;
        if groups > limits.max_groups {
            return Err(budget("decay emission groups"));
        }
        check_values(assembly.nuclides().len(), groups, limits)?;
        if !group_edges_j.iter().copied().all(nonnegative)
            || !group_edges_j.windows(2).all(|pair| pair[0] < pair[1])
        {
            return Err(invalid(
                "finite nonnegative ascending emission energies required",
            ));
        }
        let mut particles = Vec::with_capacity(assembly.nuclides().len() * groups);
        let mut energy = Vec::with_capacity(particles.capacity());
        let mut coverage = Vec::with_capacity(assembly.nuclides().len());
        let mut work = 0;
        for evaluation in assembly.evaluations() {
            visit(&mut work, limits, &mut checkpoint)?;
            let mut counts = vec![PositiveIntegral::new(); groups];
            let mut energies = counts.clone();
            let status = if evaluation.law().decay_constant().value() == 0.0 {
                DecayEmissionCoverage::Stable
            } else if let Some(spectrum) = evaluation
                .spectra()
                .iter()
                .find(|spectrum| spectrum.radiation == radiation)
            {
                integrate_spectrum(
                    spectrum,
                    &group_edges_j,
                    &mut counts,
                    &mut energies,
                    &mut work,
                    limits,
                    &mut checkpoint,
                )?;
                DecayEmissionCoverage::Reported
            } else {
                DecayEmissionCoverage::Unreported
            };
            coverage.push(status);
            for (count, joules) in counts.into_iter().zip(energies) {
                visit(&mut work, limits, &mut checkpoint)?;
                // A positive coefficient outside f64 range refuses explicitly;
                // it cannot become an invented zero before inventory scaling.
                particles.push(count.value()?);
                energy.push(joules.value()?);
            }
        }
        checkpoint()?;
        Ok(Self {
            assembly,
            radiation,
            edges: group_edges_j.into(),
            coverage: coverage.into(),
            particles: particles.into(),
            energy: energy.into(),
            limits,
        })
    }

    pub fn assembly(&self) -> &Arc<EvaluatedDecayAssembly> {
        &self.assembly
    }
    pub fn radiation(&self) -> DecayRadiation {
        self.radiation
    }
    pub fn group_edges_j(&self) -> &[f64] {
        &self.edges
    }
    pub fn coverage(&self) -> &[DecayEmissionCoverage] {
        &self.coverage
    }

    /// One original-order parent's expected particle counts per transformation.
    pub fn particles_per_decay(&self, parent: usize) -> Option<&[f64]> {
        let groups = self.edges.len() - 1;
        (parent < self.coverage.len())
            .then(|| &self.particles[parent * groups..(parent + 1) * groups])
    }

    /// Exact group energy integrals of the admitted interpolation, in J/decay.
    pub fn energy_per_decay_j(&self, parent: usize) -> Option<&[f64]> {
        let groups = self.edges.len() - 1;
        (parent < self.coverage.len()).then(|| &self.energy[parent * groups..(parent + 1) * groups])
    }

    /// Use actual integrated parent transformations, including daughters born
    /// and decayed within the interval. A report carries no model identity, so
    /// the caller must supply its producing inventory's original model Arc.
    pub fn from_decay_step(
        &self,
        model: &Arc<DecayModel>,
        step: &DecayStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayEmissionStep, EndfError> {
        checkpoint()?;
        if !Arc::ptr_eq(model, self.assembly.model()) {
            return Err(invalid(
                "emission report belongs to a different decay model",
            ));
        }
        let mut work = 0;
        self.apply(step.duration_s, &step.decays, &mut work, &mut checkpoint)
    }

    /// Bind the original evaluated decay channels inside a full burnup model.
    /// Capture/fission channels may be added and reactions may be reordered.
    /// Nuclide order/laws and every original decay branch, product, source label
    /// and energy route must still agree. Matching labels alone are insufficient.
    pub fn bind_burnup(
        &self,
        model: Arc<BurnupModel>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupDecayEmission, EndfError> {
        checkpoint()?;
        if model.nuclides().len() != self.assembly.nuclides().len()
            || model
                .nuclides()
                .iter()
                .zip(self.assembly.nuclides())
                .zip(self.assembly.model().nuclides())
                .any(|((actual, &id), original)| actual.id != id || actual.law != original.law)
        {
            return Err(invalid(
                "burnup identities or half-lives differ from the emission evaluation",
            ));
        }
        let mut work = 0;
        let mut expected = BTreeMap::new();
        for (branch, &(record, mode)) in self
            .assembly
            .model()
            .branches()
            .iter()
            .zip(&self.assembly.mode_sources)
        {
            visit(&mut work, self.limits, &mut checkpoint)?;
            let source = &self.assembly.evaluations()[record];
            let label = format!(
                "{} | {} | MAT{} MF8/MT457 mode{}",
                self.assembly.model().data_id(),
                source.data_id(),
                source.metadata().material_number,
                mode
            );
            if expected.insert((branch.parent, label), branch).is_some() {
                return Err(invalid("ambiguous evaluated decay channel identity"));
            }
        }
        let mut decay_reactions = Vec::with_capacity(expected.len());
        for (index, reaction) in model.reactions().iter().enumerate() {
            visit(&mut work, self.limits, &mut checkpoint)?;
            if let BurnupRate::Decay { branch_fraction } = reaction.rate {
                let branch = expected
                    .remove(&(reaction.parent, reaction.data_id.clone()))
                    .ok_or_else(|| {
                        invalid("burnup contains an unbound or duplicate decay channel")
                    })?;
                if branch_fraction != branch.probability
                    || reaction.energy_per_event != branch.energy_per_decay
                    || reaction.products.len() != 1
                    || reaction.products[0].nuclide != branch.daughter
                    || reaction.products[0].yield_per_event != 1.0
                {
                    return Err(invalid(
                        "burnup decay branch differs from its evaluated source",
                    ));
                }
                decay_reactions.push((index, reaction.parent));
            }
        }
        if !expected.is_empty() {
            return Err(invalid("burnup omits an evaluated decay branch"));
        }
        checkpoint()?;
        Ok(BurnupDecayEmission {
            source: self.clone(),
            model,
            decay_reactions: decay_reactions.into(),
        })
    }

    fn apply(
        &self,
        duration_s: f64,
        decays: &[f64],
        work: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayEmissionStep, EndfError> {
        if !nonnegative(duration_s)
            || decays.len() != self.coverage.len()
            || !decays.iter().copied().all(nonnegative)
            || (duration_s == 0.0 && decays.iter().any(|&count| count != 0.0))
        {
            return Err(invalid(
                "matching finite integrated decay counts and duration required",
            ));
        }
        let groups = self.edges.len() - 1;
        let mut particles = vec![PositiveIntegral::new(); groups];
        let mut energy = particles.clone();
        let mut unreported_decays = Vec::new();
        for (parent, &count) in decays.iter().enumerate() {
            visit(work, self.limits, checkpoint)?;
            if count == 0.0 {
                continue;
            }
            match self.coverage[parent] {
                DecayEmissionCoverage::Stable => {
                    return Err(invalid("a stable nuclide cannot have decay events"));
                }
                DecayEmissionCoverage::Unreported => {
                    unreported_decays.push(UnreportedDecayEmission {
                        nuclide: self.assembly.nuclides()[parent],
                        decays: count,
                    });
                }
                DecayEmissionCoverage::Reported => {
                    let offset = parent * groups;
                    for group in 0..groups {
                        visit(work, self.limits, checkpoint)?;
                        if count > 1.0
                            && (self.particles[offset + group].is_subnormal()
                                || self.energy[offset + group].is_subnormal())
                        {
                            return Err(FissionError::Arithmetic.into());
                        }
                        particles[group].add_product(
                            1.0,
                            count,
                            self.particles[offset + group],
                            1.0,
                        )?;
                        energy[group].add_product(1.0, count, self.energy[offset + group], 1.0)?;
                    }
                }
            }
        }
        let emitted_particles = particles
            .into_iter()
            .map(PositiveIntegral::value)
            .collect::<Result<_, _>>()?;
        let emitted_energy_j = energy
            .into_iter()
            .map(PositiveIntegral::value)
            .collect::<Result<_, _>>()?;
        checkpoint()?;
        Ok(DecayEmissionStep {
            duration_s,
            radiation: self.radiation,
            group_edges_j: Arc::clone(&self.edges),
            emitted_particles,
            emitted_energy_j,
            unreported_decays,
        })
    }
}

/// Decay emission response bound to an admitted full reaction model. Nuclear
/// capture, fission, feed and movement remain owned by the original solvers.
#[derive(Clone, Debug)]
pub struct BurnupDecayEmission {
    source: PreparedDecayEmission,
    model: Arc<BurnupModel>,
    decay_reactions: Arc<[(usize, usize)]>,
}

impl BurnupDecayEmission {
    pub fn source(&self) -> &PreparedDecayEmission {
        &self.source
    }
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }

    /// Select physical decay reactions only; capture/fission removals and
    /// transport counts never become emissions. Branch probabilities have
    /// already been applied by the original reaction owner.
    pub fn from_burnup_step(
        &self,
        model: &Arc<BurnupModel>,
        step: &BurnupStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayEmissionStep, EndfError> {
        checkpoint()?;
        self.check_model(model)?;
        self.apply(step, &mut 0, &mut checkpoint)
    }

    /// Preserve the location of actual events in original compartment order.
    /// Reported movement does not relocate already emitted radiation. All
    /// candidates are returned together or refused together.
    pub fn from_transport_step(
        &self,
        model: &Arc<BurnupModel>,
        step: &BurnupTransportStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<DecayEmissionStep>, EndfError> {
        checkpoint()?;
        self.check_model(model)?;
        check_values(
            step.burnup.len(),
            self.source.edges.len() - 1,
            self.source.limits,
        )?;
        if !nonnegative(step.duration_s)
            || step.burnup.is_empty()
            || step
                .burnup
                .iter()
                .any(|local| local.duration_s != step.duration_s)
        {
            return Err(invalid(
                "matching nonempty compartment decay-report intervals required",
            ));
        }
        let mut work = 0;
        let mut result = Vec::with_capacity(step.burnup.len());
        for local in &step.burnup {
            result.push(self.apply(local, &mut work, &mut checkpoint)?);
        }
        checkpoint()?;
        Ok(result)
    }

    fn check_model(&self, model: &Arc<BurnupModel>) -> Result<(), EndfError> {
        if Arc::ptr_eq(model, &self.model) {
            Ok(())
        } else {
            Err(invalid(
                "emission report belongs to a different burnup model",
            ))
        }
    }

    fn apply(
        &self,
        step: &BurnupStep,
        work: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayEmissionStep, EndfError> {
        if step.reaction_counts.len() != self.model.reactions().len()
            || !step.reaction_counts.iter().copied().all(nonnegative)
        {
            return Err(invalid(
                "matching finite original-order reaction counts required",
            ));
        }
        let mut counts = vec![Accumulator::new(0.0); self.model.nuclides().len()];
        for &(reaction, parent) in self.decay_reactions.iter() {
            visit(work, self.source.limits, checkpoint)?;
            counts[parent].add(step.reaction_counts[reaction])?;
        }
        let decays: Vec<_> = counts.into_iter().map(|value| value.value()).collect();
        self.source
            .apply(step.duration_s, &decays, work, checkpoint)
    }
}

fn integrate_spectrum(
    spectrum: &DecaySpectrum,
    edges: &[f64],
    counts: &mut [PositiveIntegral],
    energies: &mut [PositiveIntegral],
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    let fd = spectrum.discrete_normalization.value;
    if fd > 0.0 {
        for line in spectrum.discrete.iter() {
            visit(work, limits, checkpoint)?;
            let ri = line.intensity.value;
            let positrons = spectrum.radiation == DecayRadiation::ElectronCapturePositron
                && line.pair_or_positron.is_some_and(|v| v.value > 0.0);
            if ri == 0.0 && !positrons {
                continue;
            }
            if !matches!(
                spectrum.radiation,
                DecayRadiation::Gamma
                    | DecayRadiation::Alpha
                    | DecayRadiation::Neutron
                    | DecayRadiation::Proton
                    | DecayRadiation::Electron
                    | DecayRadiation::XRay
            ) {
                return Err(EndfError::Unsupported(
                    "this discrete radiation record is not an admitted monoenergetic particle",
                ));
            }
            let e = line.energy_j.value;
            check_domain(edges, e, e)?;
            let group = edges
                .partition_point(|&edge| edge <= e)
                .saturating_sub(1)
                .min(counts.len() - 1);
            counts[group].add_product(1.0, fd, ri, 1.0)?;
            energies[group].add_triple_product(1.0, fd, ri, e, 1.0)?;
        }
    }
    let fc = spectrum.continuum_normalization.value;
    if fc == 0.0 {
        return Ok(());
    }
    let continuum = spectrum
        .continuum
        .as_ref()
        .ok_or_else(|| invalid("positive continuum normalization requires a retained spectrum"))?;
    let mut region = 0;
    let mut group = 0;
    for (index, points) in continuum.points.windows(2).enumerate() {
        visit(work, limits, checkpoint)?;
        let [a, b] = [points[0], points[1]];
        if a.energy_j == b.energy_j {
            continue;
        }
        while continuum.interpolation[region].end_point < index + 2 {
            region += 1;
        }
        let law = continuum.interpolation[region].law;
        if spectrum.radiation == DecayRadiation::ElectronCapturePositron
            && (a.density_per_j > 0.0 || (law != 1 && b.density_per_j > 0.0))
        {
            return Err(EndfError::Unsupported(
                "electron-capture/positron continuum requires resolved positron-intensity normalization",
            ));
        }
        if a.density_per_j == 0.0 && (b.density_per_j == 0.0 || law == 1) {
            continue;
        }
        if matches!(
            spectrum.radiation,
            DecayRadiation::FissionFragment | DecayRadiation::Other(_)
        ) {
            return Err(EndfError::Unsupported(
                "unresolved fragment or unidentified radiation spectrum",
            ));
        }
        if !matches!(law, 1 | 2) {
            return Err(EndfError::Unsupported(
                "decay emission grouping supports histogram and linear spectra",
            ));
        }
        check_domain(edges, a.energy_j, b.energy_j)?;
        while group + 1 < counts.len() && edges[group + 1] <= a.energy_j {
            group += 1;
        }
        let mut left = a.energy_j;
        while left < b.energy_j {
            visit(work, limits, checkpoint)?;
            let right = b.energy_j.min(edges[group + 1]);
            let p = density_at(a, b, law, left)?;
            let q = density_at(a, b, law, right)?;
            let width = right - left;
            counts[group].add_product(width, fc, p, 2.0)?;
            counts[group].add_product(width, fc, q, 2.0)?;
            // Integral of the product of two linear Bernstein polynomials,
            // E and RP(E). Every term is positive; no slope cancellation.
            energies[group].add_triple_product(width, fc, p, left, 3.0)?;
            energies[group].add_triple_product(width, fc, q, right, 3.0)?;
            energies[group].add_triple_product(width, fc, p, right, 6.0)?;
            energies[group].add_triple_product(width, fc, q, left, 6.0)?;
            left = right;
            if left < b.energy_j {
                group += 1;
            }
        }
    }
    Ok(())
}

fn density_at(
    a: DecaySpectrumPoint,
    b: DecaySpectrumPoint,
    law: u8,
    energy: f64,
) -> Result<f64, EndfError> {
    if law == 1 || energy == a.energy_j {
        return Ok(a.density_per_j);
    }
    if energy == b.energy_j {
        return Ok(b.density_per_j);
    }
    let width = b.energy_j - a.energy_j;
    let lower = (b.energy_j - energy) / width;
    let upper = (energy - a.energy_j) / width;
    if lower == 0.0 || upper == 0.0 {
        return Err(FissionError::Arithmetic.into());
    }
    convex(a.density_per_j, b.density_per_j, lower, upper)
}

fn check_domain(edges: &[f64], first: f64, last: f64) -> Result<(), EndfError> {
    if first < edges[0] || last > edges[edges.len() - 1] {
        Err(invalid(
            "positive decay spectrum lies outside the requested energy groups",
        ))
    } else {
        Ok(())
    }
}

fn check_values(rows: usize, groups: usize, limits: DecayEmissionLimits) -> Result<(), EndfError> {
    let values = rows
        .checked_mul(groups)
        .and_then(|n| n.checked_mul(2))
        .ok_or_else(|| budget("decay emission group values"))?;
    if values > limits.max_group_values {
        Err(budget("decay emission group values"))
    } else {
        Ok(())
    }
}

fn visit(
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    if *work >= limits.max_work {
        return Err(budget("decay emission work"));
    }
    *work += 1;
    if *work % 64 == 0 {
        checkpoint()?;
    }
    Ok(())
}

fn invalid(message: &'static str) -> EndfError {
    FissionError::InvalidInput(message).into()
}

#[cfg(test)]
mod tests;
