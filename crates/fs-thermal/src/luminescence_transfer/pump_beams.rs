//! Simultaneous, photon-conservative pumping of shared finite populations.
//!
//! Beams are incoherent or phase averaged: their absorption hazards add. Each
//! original cell has one ground/singlet/triplet inventory, even when several
//! directed paths cross it or traverse the cells in opposite orders. Molecular
//! depths are cross section times segment length divided by ORIGINAL cell
//! volume, per ground emitter. Cell averaging and quasistatic photon flight are
//! the same physical closure as [`super::pump_column`].
//!
//! Acyclic path dependencies admit an exact ordered sweep of the local residence
//! closures. Cyclic dependencies use bounded physical source iteration on the
//! common mean ground field. Every trial kinetic history restarts from the same
//! substep populations; no beam advances an independent copy. Acceptance checks
//! common-field residence and each beam's photon balance. Temporal refinement
//! discards unconverged grids and publishes only a converged finer history.
//! Phase interference, reflection, scattering, stimulated emission, excited
//! state absorption and active recycling of emitted light are outside scope.

use super::pump_column::{
    KineticTotals, PumpColumnLimits, Sum, admit_energy, admit_photons,
    cap, cell, energies, photon_balance, poll, product3, values,
};
use super::saturation::{SaturableSingletTripletPlan, SaturableSingletTripletStep};
use super::{SingletTripletDecay, SingletTripletError, finite, nonnegative, positive_product};

/// One original cell visited once by this beam. Different beams may share it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturablePumpSegment<const B: usize> {
    pub cell_index: usize,
    /// Lane then destination `[singlet, triplet]`, dimensionless per emitter.
    pub optical_depth_per_ground: [[f64; 2]; B],
    pub passive_optical_depth: [f64; B],
}

/// Entry-to-exit order and original photon energies of one incoherent beam.
#[derive(Clone, Debug, PartialEq)]
pub struct SaturablePumpPath<const B: usize> {
    pub photon_energies_j: [f64; B],
    pub segments: Vec<SaturablePumpSegment<B>>,
}

/// Original scalar/kinetic/time limits plus bounded cyclic source iteration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpBeamsLimits {
    pub column: PumpColumnLimits,
    /// Per temporal substep, including rejected coarse grids. Aggregate kinetic
    /// work remains bounded by `column.maximum_cell_evaluations`.
    pub maximum_coupling_iterations: usize,
}

impl Default for PumpBeamsLimits {
    fn default() -> Self {
        Self {
            column: PumpColumnLimits::default(),
            maximum_coupling_iterations: 64,
        }
    }
}

/// External photons and energy of one beam; internal segment crossings are not
/// extra sources. Arrays retain this beam's original spectral-lane order.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpBeamStep<const B: usize> {
    pub incident_photons: [f64; B],
    pub transmitted_photons: [f64; B],
    pub absorbed_photons: [f64; B],
    pub passive_absorbed_photons: [f64; B],
    pub incident_energy_j: [f64; B],
    pub transmitted_energy_j: [f64; B],
    pub absorbed_energy_j: [f64; B],
    pub passive_absorbed_energy_j: [f64; B],
    pub photon_residual: [f64; B],
}

/// One candidate shared history, in original model and beam order.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpBeamsStep<const B: usize> {
    pub cells: Vec<SaturableSingletTripletStep>,
    pub passive_heat_j: Vec<f64>,
    pub beams: Vec<PumpBeamStep<B>>,
    pub residual_j: f64,
    pub substeps: usize,
    /// Includes all local scalar, coupling and discarded temporal work.
    pub cell_evaluations: usize,
    pub coupling_iterations: usize,
    pub max_energy_discrepancy_j: f64,
}

#[derive(Clone, Copy, Debug)]
struct Visit {
    beam: usize,
    segment: usize,
}

/// Immutable geometry, spectra and constitutive models. Populations and the
/// incident fields belong to the caller and are never modified or cached.
#[derive(Clone, Debug)]
pub struct SaturablePumpBeams<const B: usize> {
    models: Vec<SingletTripletDecay>,
    paths: Vec<SaturablePumpPath<B>>,
    visits: Vec<Vec<Visit>>,
    topological: Option<Vec<usize>>,
}

impl<const B: usize> SaturablePumpBeams<B> {
    pub fn try_new(
        models: &[SingletTripletDecay],
        paths: &[SaturablePumpPath<B>],
        maximum_cells: usize,
        maximum_beams: usize,
        maximum_segments: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        poll(&mut cancelled)?;
        if B == 0 || B > 64 || maximum_cells == 0
            || maximum_beams == 0 || maximum_segments == 0
        {
            return Err(SingletTripletError::Invalid("pump beam lanes or admission limits"));
        }
        cap("pump beam cells", models.len(), maximum_cells)?;
        cap("pump beams", paths.len(), maximum_beams)?;
        let mut count = 0usize;
        for path in paths {
            poll(&mut cancelled)?;
            count = count.checked_add(path.segments.len())
                .ok_or(SingletTripletError::Unrepresentable("pump beam segment count"))?;
            cap("pump beam segments", count, maximum_segments)?;
        }
        count.checked_mul(B)
            .ok_or(SingletTripletError::Unrepresentable("pump beam channel count"))?;
        let mut visits = vec![Vec::new(); models.len()];
        let mut successors: Vec<Vec<usize>> = vec![Vec::new(); models.len()];
        let mut indegree = vec![0usize; models.len()];
        let mut last_visit = vec![None; models.len()];
        for (beam, path) in paths.iter().enumerate() {
            poll(&mut cancelled)?;
            if !path.photon_energies_j.iter().all(|&e| e.is_finite() && e > 0.0) {
                return Err(SingletTripletError::Invalid("pump beam photon energies"));
            }
            let mut previous: Option<usize> = None;
            for (ordinal, segment) in path.segments.iter().enumerate() {
                poll(&mut cancelled)?;
                let index = segment.cell_index;
                if index >= models.len() || last_visit[index] == Some(beam) {
                    return Err(SingletTripletError::Invalid("pump beam cell index or repeated cell"));
                }
                last_visit[index] = Some(beam);
                if !segment.optical_depth_per_ground.iter().flatten().copied().all(nonnegative)
                    || !segment.passive_optical_depth.iter().copied().all(nonnegative)
                {
                    return Err(SingletTripletError::Invalid("pump beam optical depths"));
                }
                for lane in 0..B {
                    finite(
                        segment.optical_depth_per_ground[lane].iter().sum(),
                        "summed active depth",
                    )?;
                    for target in 0..2 {
                        if segment.optical_depth_per_ground[lane][target] > 0.0
                            && path.photon_energies_j[lane]
                                < models[index].excitation_energies_j()[target]
                        {
                            return Err(SingletTripletError::Invalid(
                                "pump photon cannot populate its destination",
                            ));
                        }
                    }
                }
                visits[index].push(Visit { beam, segment: ordinal });
                if let Some(parent) = previous {
                    successors[parent].push(index);
                    indegree[index] += 1;
                }
                previous = Some(index);
            }
        }
        let mut order: Vec<_> = indegree.iter().enumerate()
            .filter_map(|(i, &degree)| (degree == 0).then_some(i)).collect();
        let mut next = 0;
        while next < order.len() {
            poll(&mut cancelled)?;
            for &child in &successors[order[next]] {
                indegree[child] -= 1;
                if indegree[child] == 0 {
                    order.push(child);
                }
            }
            next += 1;
        }
        let topological = (order.len() == models.len()).then_some(order);
        poll(&mut cancelled)?;
        Ok(Self {
            models: models.to_vec(),
            paths: paths.to_vec(),
            visits,
            topological,
        })
    }

    pub fn models(&self) -> &[SingletTripletDecay] {
        &self.models
    }

    pub fn paths(&self) -> &[SaturablePumpPath<B>] {
        &self.paths
    }

    /// Constant simultaneous incident photon rates during this interval. Models,
    /// geometry and spectral coefficients are fixed until the caller refreshes
    /// them. Every error leaves all caller populations and accounts untouched.
    pub fn advance(
        &self,
        populations: &[[f64; 3]],
        incident_photon_rates_s: &[[f64; B]],
        duration_s: f64,
        limits: PumpBeamsLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PumpBeamsStep<B>, SingletTripletError> {
        poll(&mut cancelled)?;
        let column = limits.column;
        if populations.len() != self.models.len()
            || incident_photon_rates_s.len() != self.paths.len()
            || !populations.iter().flatten().copied().all(nonnegative)
            || !incident_photon_rates_s.iter().flatten().copied().all(nonnegative)
            || !nonnegative(duration_s) || limits.maximum_coupling_iterations == 0
            || column.maximum_root_iterations < 2 || column.maximum_substeps == 0
            || column.maximum_cell_evaluations == 0 || !nonnegative(column.relative_tolerance)
            || column.relative_tolerance >= 1.0 || !nonnegative(column.absolute_energy_tolerance_j)
            || (column.relative_tolerance == 0.0 && column.absolute_energy_tolerance_j == 0.0)
        {
            return Err(SingletTripletError::Invalid(
                "pump beam populations, illumination or work limits",
            ));
        }
        SaturableSingletTripletPlan::prepare(
            SingletTripletDecay::try_new([1.0; 2], [0.0; 2], [0.0; 2], [0.0; 2])?,
            [0.0; 2], [0.0; 2], 0.0, column.kinetic, &mut cancelled,
        )?;
        let nonlinear = duration_s > 0.0
            && self.paths.iter().zip(incident_photon_rates_s).any(|(path, rates)| {
                path.segments.iter().any(|segment| {
                    (0..B).any(|lane| rates[lane] > 0.0
                        && segment.optical_depth_per_ground[lane].iter().any(|&a| a > 0.0))
                })
            });
        let mut evaluations = 0usize;
        let mut coupling_iterations = 0usize;
        let mut substeps = 1usize;
        let mut previous: Option<Vec<f64>> = None;
        loop {
            poll(&mut cancelled)?;
            cap("pump beam temporal substeps", substeps, column.maximum_substeps)?;
            let candidate = self.uniform(
                populations, incident_photon_rates_s, duration_s, substeps, nonlinear,
                limits, &mut evaluations, &mut coupling_iterations, &mut cancelled,
            )?;
            if let Some(mut candidate) = candidate {
                candidate.cell_evaluations = evaluations;
                candidate.coupling_iterations = coupling_iterations;
                if !nonlinear {
                    return Ok(candidate);
                }
                let components = self.components(&candidate, duration_s)?;
                if let Some(old) = &previous {
                    let mut converged = true;
                    let mut discrepancy = 0.0_f64;
                    for (&before, &after) in old.iter().zip(&components) {
                        let difference = (after - before).abs();
                        discrepancy = discrepancy.max(difference);
                        converged &= difference <= column.absolute_energy_tolerance_j
                            + column.relative_tolerance * before.abs().max(after.abs());
                    }
                    if converged {
                        candidate.max_energy_discrepancy_j = discrepancy;
                        poll(&mut cancelled)?;
                        return Ok(candidate);
                    }
                }
                previous = Some(components);
            } else {
                // A failed common-field solve is not a physical history and
                // cannot be the coarse side of an accuracy comparison.
                previous = None;
            }
            substeps = substeps.checked_mul(2)
                .ok_or(SingletTripletError::Unrepresentable("pump beam refinement counter"))?;
        }
    }
}

#[path = "pump_beams/evolution.rs"]
mod evolution;
#[path = "pump_beams/accounts.rs"]
mod accounts;
#[cfg(test)]
#[path = "pump_beams/tests.rs"]
mod tests;
#[cfg(test)]
#[path = "pump_beams/oracles.rs"]
mod oracles;
