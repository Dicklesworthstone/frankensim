//! State-refreshed exponential midpoint on the original species/enthalpy fluxes.
//!
//! Each predictor supplies midpoint properties, while the corrector advances
//! FROM the original substep state. Complete intervals at n and 2n substeps
//! estimate smooth second-order time error; the accepted state is the unchanged
//! fine solution. Source discontinuities and phase corners require resolution
//! and are not covered by a certified error bound.

use super::super::{
    MAX_MATERIAL_SPECIES, SpeciesBoundaryDiffusionReport, SpeciesDiffusionBoundary, ratio,
};
use super::boundary::BoundaryEnthalpyAccount;
use super::{
    CancelGate, ChtError, EnergySum, MaterialAdvectionStep, SpeciesDiffusionControls,
    SpeciesEnthalpyDiffusionStep, VoxelDomain, nonnegative, poll, poll_every, refuse,
};
use std::cell::Cell;

mod boundary;
pub use boundary::{
    AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls,
    AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionStep,
};
use boundary::{BoundaryHistory, BoundaryOptions};

/// One original caloric/transport owner for the complete supplied state.
/// Preparations must own their state-dependent data, retain the same species
/// order and enthalpy convention, and poll the supplied cancellation gate.
/// Properties are functions of this state; no unrepresented time dependence
/// or mutable constitutive history may be advanced by a query.
pub trait SpeciesEnthalpyDiffusionModel<const N: usize> {
    type Prepared;

    fn prepare(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<Self::Prepared, ChtError>;
    fn temperature_k<'p>(&self, prepared: &'p Self::Prepared) -> &'p [f64];
    fn diffusivity_m2_s<'p>(&self, prepared: &'p Self::Prepared) -> &'p [f64];
    fn face_partial_enthalpy_j_kg(
        &self,
        prepared: &Self::Prepared,
        domain: &VoxelDomain,
        cells: [usize; 2],
        axis: usize,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError>;
    /// Original partial specific enthalpies at a declared physical reservoir face.
    /// The descriptor and its canonical index refer to the SAME supplied boundary
    /// list. A model must explicitly bind this thermal law; reservoir composition
    /// alone does not supply a face temperature or an extra heat source.
    fn boundary_partial_enthalpy_j_kg(
        &self,
        _prepared: &Self::Prepared,
        _domain: &VoxelDomain,
        _boundary_index: usize,
        _boundary: &SpeciesDiffusionBoundary<N>,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        poll(gate)?;
        Err(refuse(
            "adaptive reservoirs require an original boundary partial-enthalpy law",
        ))
    }
    /// Conservative owned payload of one preparation, including its fixed
    /// value. Borrowed original source tables are caller-owned and additional.
    fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize>;
    /// Actual mass-weighted absolute caloric-inversion residual [J]. Analytic
    /// closures may return zero; a numerical source must report its residual.
    fn caloric_defect_j(&self, prepared: &Self::Prepared) -> Result<f64, ChtError>;
}

/// Independent accuracy, physical-account and total resource allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSpeciesEnthalpyDiffusionControls<const N: usize> {
    /// Original inner face-sweep controls. Pair work and storage bound the
    /// WHOLE adaptive call, including predictors and rejected resolutions.
    /// Species mass allowances cover cumulative accepted-substep defects and
    /// the final initial-to-endpoint closed-domain audit.
    pub solver: SpeciesDiffusionControls<N>,
    /// Initial number of physical property-refresh substeps. The original
    /// solver's `substeps` independently resolves each frozen face composition.
    pub initial_substeps: usize,
    /// Number of fine solves; each doubles the previous outer resolution.
    pub maximum_refinements: usize,
    /// Total preparations PLUS face-property queries, including rejected work.
    pub maximum_source_evaluations: usize,
    pub absolute_temperature_tolerance_k: f64,
    /// Absolute J/kg only: formation energy never scales the thermal tolerance.
    /// Compare both represented endpoints and retained interval heat / mass.
    pub absolute_specific_enthalpy_tolerance_j_kg: f64,
    pub absolute_mass_fraction_tolerance: [f64; N],
    pub relative_mass_fraction_tolerance: f64,
    /// Cumulative accepted-corrector defects, the complete interval audit and
    /// original/final caloric closure must fit this absolute joule allowance.
    /// Every source-prepared stage must also satisfy it independently.
    pub maximum_enthalpy_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveSpeciesEnthalpyDiffusionStep<const N: usize, P = ()> {
    /// Accepted endpoint source image, including its original transport/EOS
    /// fields. Reuse it without another whole-domain caloric inversion.
    /// Editing returned inventories requires a fresh source preparation.
    pub state: P,
    pub species_mass_kg: Vec<[f64; N]>,
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// Actual final source inversion, in K; solid entries are zero.
    pub temperature_k: Vec<f64>,
    /// Integrated paired heat already included in the accepted enthalpy [J].
    /// This scalar readback rounds once; internal interval audits retain tails.
    pub diffusive_enthalpy_change_j: Vec<f64>,
    pub report: AdaptiveSpeciesEnthalpyDiffusionReport<N>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveSpeciesEnthalpyDiffusionReport<const N: usize> {
    pub duration_s: f64,
    pub accepted_substeps: usize,
    pub refinements: usize,
    pub total_pair_updates: usize,
    pub total_component_updates: usize,
    pub source_preparations: usize,
    pub face_source_evaluations: usize,
    pub estimated_temperature_error_k: f64,
    pub estimated_specific_enthalpy_error_j_kg: f64,
    pub estimated_mass_fraction_error: [f64; N],
    pub maximum_error_ratio: f64,
    /// Full-outer face relaxation on substeps with sampled changing D or
    /// relative partial enthalpies. Admission requires at most one half;
    /// identical original coefficients retain their exact stiff fast path.
    pub maximum_changing_source_relaxation: f64,
    /// Full-INNER relaxation where a reservoir competes with another reservoir
    /// or an interior map on its cell. Admission requires at most one half,
    /// even for unchanged sources: saturated maps can misassign external heat.
    /// Explicitly resolved inner sweeps and isolated exact maps retain their
    /// stiff fast paths; changing sources independently require the outer cap.
    pub maximum_competing_boundary_relaxation: f64,
    pub absolute_mass_defect_kg: [f64; N],
    pub stored_enthalpy_change_j: f64,
    pub transported_enthalpy_change_j: f64,
    /// Sum of the original physical correctors' local/global/frame defects.
    /// Predictors and rejected resolutions are not physical intervals.
    pub accepted_substep_enthalpy_defect_j: f64,
    /// Independent original-to-final represented storage versus paired heat.
    pub interval_enthalpy_defect_j: f64,
    pub initial_caloric_defect_j: f64,
    pub endpoint_caloric_defect_j: f64,
    /// Largest source closure residual on this accepted trajectory, including
    /// its property predictors. Rejected trajectories are excluded.
    pub maximum_stage_caloric_defect_j: f64,
    pub absolute_enthalpy_defect_j: f64,
}

struct Work {
    pair_updates: usize,
    component_updates: usize,
    source_preparations: Cell<usize>,
    face_queries: Cell<usize>,
    maximum_sources: usize,
}

impl Work {
    fn source(&self, preparation: bool) -> Result<(), ChtError> {
        let used = self
            .source_preparations
            .get()
            .checked_add(self.face_queries.get())
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| refuse("adaptive caloric source count overflows"))?;
        if used > self.maximum_sources {
            return Err(refuse("adaptive caloric source budget exhausted"));
        }
        let count = if preparation {
            &self.source_preparations
        } else {
            &self.face_queries
        };
        count.set(count.get() + 1);
        Ok(())
    }
}

struct FrozenBoundary<const N: usize> {
    mass: SpeciesBoundaryDiffusionReport<N>,
    account: BoundaryEnthalpyAccount<N>,
}

struct FrozenSolve<const N: usize> {
    step: SpeciesEnthalpyDiffusionStep<N>,
    increments: Vec<EnergySum>,
    boundary: Option<FrozenBoundary<N>>,
}

struct Refined<const N: usize, P> {
    step: AdaptiveSpeciesEnthalpyDiffusionStep<N, P>,
    boundary: Option<BoundaryHistory<N>>,
    boundary_mass_error: [f64; N],
    boundary_heat_error: f64,
}

struct Trial<const N: usize> {
    species: Vec<[f64; N]>,
    h: Vec<f64>,
    temperature: Vec<f64>,
    heat: Vec<EnergySum>,
    boundary: Option<BoundaryHistory<N>>,
    mass_defect: [f64; N],
    stored: f64,
    transported: f64,
    substep_defect: f64,
    interval_defect: f64,
    endpoint_caloric: f64,
    maximum_stage_caloric: f64,
    maximum_source_relaxation: f64,
    maximum_boundary_relaxation: f64,
    total_defect: f64,
}

// Compose positive face relaxations without forming 1-r and losing a tiny r.
// The exponent includes ALL inner sweeps: more frozen sweeps cannot resolve
// an outer interval whose caloric properties have not yet been refreshed.
fn outer_relaxation(mut value: f64, mut sweeps: u64) -> f64 {
    let mut result = 0.0;
    while sweeps != 0 {
        if sweeps & 1 != 0 {
            result = value.mul_add(1.0 - result, result);
        }
        sweeps >>= 1;
        if sweeps != 0 {
            value = value.mul_add(1.0 - value, value);
        }
    }
    result
}

fn changed_contrasts<const N: usize>(a: [f64; N], b: [f64; N]) -> Result<bool, ChtError> {
    if a.iter().chain(&b).any(|h| !h.is_finite()) {
        return Err(refuse(
            "stiff-source comparison requires finite partial enthalpy",
        ));
    }
    let mut changed = false;
    for k in 0..N {
        let left = fs_math::eft::two_sum(a[k], -a[0]);
        let right = fs_math::eft::two_sum(b[k], -b[0]);
        if !left.0.is_finite()
            || !left.1.is_finite()
            || !right.0.is_finite()
            || !right.1.is_finite()
        {
            return Err(refuse("stiff-source enthalpy contrast is unrepresentable"));
        }
        changed |= left != right;
    }
    Ok(changed)
}

fn positive(value: f64) -> Result<f64, ChtError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(refuse(
            "adaptive diffusion requires representable positive controls",
        ))
    }
}

fn prepared<const N: usize, M: SpeciesEnthalpyDiffusionModel<N>>(
    model: &M,
    state: &M::Prepared,
    domain: &VoxelDomain,
    maximum_defect_j: f64,
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    let temperature = model.temperature_k(state);
    let diffusivity = model.diffusivity_m2_s(state);
    if temperature.len() != domain.cell_count() || diffusivity.len() != domain.cell_count() {
        return Err(refuse(
            "prepared caloric and diffusion fields must match the domain",
        ));
    }
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        nonnegative(diffusivity[c])?;
        if domain.is_fluid(c) {
            positive(temperature[c])?;
        } else if temperature[c] != 0.0 || diffusivity[c] != 0.0 {
            return Err(refuse(
                "solid entries of prepared thermal fields must be zero",
            ));
        }
    }
    let defect = model.caloric_defect_j(state)?;
    nonnegative(defect)?;
    if defect > maximum_defect_j {
        return Err(refuse(
            "prepared caloric inversion exceeds its absolute joule budget",
        ));
    }
    poll(gate)?;
    Ok(defect)
}

// Three prepared images can coexist: original, current and midpoint/endpoint.
// Six complete state payloads cover retained coarse, live current, predictor,
// corrector, comparison/output and their temporary vectors conservatively.
fn solver_storage<const N: usize, M: SpeciesEnthalpyDiffusionModel<N>>(
    model: &M,
    cells: usize,
    total: usize,
) -> Result<usize, ChtError> {
    let bytes = (|| {
        let state = core::mem::size_of::<[f64; N]>().checked_add(16)?;
        let cell = state
            .checked_mul(6)?
            .checked_add(3 * core::mem::size_of::<EnergySum>() + 8)?;
        let fixed = core::mem::size_of::<AdaptiveSpeciesEnthalpyDiffusionStep<N>>()
            .checked_add(6 * core::mem::size_of::<Trial<N>>())?
            .checked_add(8 * core::mem::size_of::<[EnergySum; N]>())?
            .checked_add(4096)?;
        cells
            .checked_mul(cell)?
            .checked_add(model.maximum_prepared_bytes(cells)?.checked_mul(3)?)?
            .checked_add(fixed)
    })()
    .filter(|&n| isize::try_from(n).is_ok())
    .ok_or_else(|| refuse("adaptive caloric storage extent overflows"))?;
    total
        .checked_sub(bytes)
        .ok_or_else(|| refuse("adaptive caloric live-state storage budget exhausted"))
}

#[allow(clippy::too_many_arguments)]
fn frozen<const N: usize, M: SpeciesEnthalpyDiffusionModel<N>>(
    carrier: &MaterialAdvectionStep,
    domain: &VoxelDomain,
    model: &M,
    state: &M::Prepared,
    species: &[[f64; N]],
    h: &[f64],
    dt: f64,
    controls: AdaptiveSpeciesEnthalpyDiffusionControls<N>,
    storage: usize,
    work: &mut Work,
    gate: &CancelGate,
    comparison: Option<(&M::Prepared, &Cell<bool>)>,
    boundary: Option<BoundaryOptions<'_, N>>,
) -> Result<FrozenSolve<N>, ChtError> {
    let mut solver = controls.solver;
    solver.maximum_storage_bytes = storage;
    solver.maximum_pair_updates = solver
        .maximum_pair_updates
        .checked_sub(work.pair_updates)
        .ok_or_else(|| refuse("adaptive diffusion total pair budget exhausted"))?;
    let interior = |cells, axis| {
        poll(gate)?;
        work.source(false)?;
        let partial = model.face_partial_enthalpy_j_kg(state, domain, cells, axis, gate)?;
        if let Some((start, changed)) = comparison
            && !changed.get()
        {
            work.source(false)?;
            let original = model.face_partial_enthalpy_j_kg(start, domain, cells, axis, gate)?;
            changed.set(changed_contrasts(original, partial)?);
        }
        Ok(partial)
    };
    let result = if let Some(boundary) = boundary {
        let (open, mut account) = carrier.diffuse_species_with_enthalpy_and_boundaries_accounted(
            domain,
            species,
            h,
            model.diffusivity_m2_s(state),
            boundary.boundaries,
            dt,
            interior,
            |index| {
                poll(gate)?;
                work.source(false)?;
                let face = &boundary.boundaries[index];
                let partial =
                    model.boundary_partial_enthalpy_j_kg(state, domain, index, face, gate)?;
                if let Some((start, changed)) = comparison
                    && !changed.get()
                {
                    work.source(false)?;
                    let original =
                        model.boundary_partial_enthalpy_j_kg(start, domain, index, face, gate)?;
                    changed.set(changed_contrasts(original, partial)?);
                }
                Ok(partial)
            },
            solver,
            controls.maximum_enthalpy_defect_j,
            gate,
        )?;
        let increments = core::mem::take(&mut account.increments);
        FrozenSolve {
            step: open.step,
            increments,
            boundary: Some(FrozenBoundary {
                mass: open.boundary,
                account,
            }),
        }
    } else {
        let (step, increments) = carrier.diffuse_species_with_enthalpy_accounted(
            domain,
            species,
            h,
            model.diffusivity_m2_s(state),
            dt,
            interior,
            solver,
            controls.maximum_enthalpy_defect_j,
            gate,
        )?;
        FrozenSolve {
            step,
            increments,
            boundary: None,
        }
    };
    work.pair_updates = work
        .pair_updates
        .checked_add(result.step.species.report.pair_updates)
        .ok_or_else(|| refuse("adaptive diffusion total pair extent overflows"))?;
    work.component_updates = work
        .pair_updates
        .checked_mul(N)
        .ok_or_else(|| refuse("adaptive diffusion total component extent overflows"))?;
    Ok(result)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn integrate<const N: usize, M: SpeciesEnthalpyDiffusionModel<N>>(
    carrier: &MaterialAdvectionStep,
    domain: &VoxelDomain,
    model: &M,
    original: &M::Prepared,
    original_caloric: f64,
    initial_species: &[[f64; N]],
    initial_h: &[f64],
    duration: f64,
    substeps: usize,
    controls: AdaptiveSpeciesEnthalpyDiffusionControls<N>,
    storage: usize,
    work: &mut Work,
    gate: &CancelGate,
    boundary: Option<BoundaryOptions<'_, N>>,
) -> Result<(Trial<N>, M::Prepared), ChtError> {
    if substeps == 0 || substeps as u64 > (1u64 << 53) {
        return Err(refuse("adaptive thermal substep count is unrepresentable"));
    }
    let dt = positive(duration / substeps as f64)?;
    let half = positive(0.5 * dt)?;
    let mut species = initial_species.to_vec();
    let mut h = initial_h.to_vec();
    let mut heat = vec![EnergySum::default(); h.len()];
    let mut boundary_history = boundary.map(|b| BoundaryHistory::new(b.boundaries));
    let mut current = None;
    let mut substep_defect = EnergySum::default();
    let mut mass_defect = [EnergySum::default(); N];
    let mut endpoint_caloric = original_caloric;
    let mut maximum_stage_caloric = original_caloric;
    let mut maximum_source_relaxation: f64 = 0.0;
    let mut maximum_boundary_relaxation: f64 = 0.0;
    for _ in 0..substeps {
        poll(gate)?;
        let start = current.as_ref().unwrap_or(original);
        let FrozenSolve {
            step: predictor,
            increments: predictor_heat,
            boundary: predictor_boundary,
        } = frozen(
            carrier, domain, model, start, &species, &h, half, controls, storage, work, gate, None,
            boundary,
        )?;
        drop(predictor_heat);
        drop(predictor_boundary);
        work.source(true)?;
        let midpoint = model.prepare(
            domain,
            carrier,
            &predictor.species.species_mass_kg,
            &predictor.specific_enthalpy_j_kg,
            gate,
        )?;
        maximum_stage_caloric = maximum_stage_caloric.max(prepared(
            model,
            &midpoint,
            domain,
            controls.maximum_enthalpy_defect_j,
            gate,
        )?);
        let changed = Cell::new(false);
        for (c, (&old, &new)) in model
            .diffusivity_m2_s(start)
            .iter()
            .zip(model.diffusivity_m2_s(&midpoint))
            .enumerate()
        {
            poll_every(gate, c)?;
            changed.set(changed.get() || old != new);
        }
        // The predictor spans dt/2, with two half sweeps per inner substep.
        // Four times that count recovers the full outer face relaxation.
        let start_relaxation = outer_relaxation(
            predictor.species.report.maximum_pair_relaxation,
            4 * controls.solver.substeps as u64,
        );
        let comparison = (!changed.get() && start_relaxation > 0.5).then_some((start, &changed));
        // The predictor locates properties only. Its heat is not another source
        // and its state is NEVER the starting state of this full corrector.
        let FrozenSolve {
            step: corrected,
            increments: corrected_heat,
            boundary: corrected_boundary,
        } = frozen(
            carrier, domain, model, &midpoint, &species, &h, dt, controls, storage, work, gate,
            comparison, boundary,
        )?;
        match (&mut boundary_history, corrected_boundary) {
            (Some(history), Some(current)) => {
                history.accumulate(&current.mass, &current.account, gate)?
            }
            (None, None) => {}
            _ => return Err(refuse("adaptive reservoir account lost its original owner")),
        }
        if let Some(options) = boundary {
            // Fixed laws can still alias when a saturated map consumes the
            // inventory before another face can exchange it. Resolve the INNER
            // splitting interval; unlike source refresh, more native sweeps do
            // resolve this competition. The predictor spans only dt/2.
            let start_inner = outer_relaxation(predictor.species.report.maximum_pair_relaxation, 4);
            let midpoint_inner =
                outer_relaxation(corrected.species.report.maximum_pair_relaxation, 2);
            if start_inner > maximum_boundary_relaxation
                && options.competing_maps(domain, model.diffusivity_m2_s(start), gate)?
            {
                maximum_boundary_relaxation = start_inner;
            }
            if midpoint_inner > maximum_boundary_relaxation
                && options.competing_maps(domain, model.diffusivity_m2_s(&midpoint), gate)?
            {
                maximum_boundary_relaxation = midpoint_inner;
            }
        }
        if changed.get() {
            maximum_source_relaxation =
                maximum_source_relaxation
                    .max(start_relaxation)
                    .max(outer_relaxation(
                        corrected.species.report.maximum_pair_relaxation,
                        2 * controls.solver.substeps as u64,
                    ));
        }
        drop(midpoint);
        drop(predictor);
        work.source(true)?;
        let endpoint = model.prepare(
            domain,
            carrier,
            &corrected.species.species_mass_kg,
            &corrected.specific_enthalpy_j_kg,
            gate,
        )?;
        endpoint_caloric = prepared(
            model,
            &endpoint,
            domain,
            controls.maximum_enthalpy_defect_j,
            gate,
        )?;
        maximum_stage_caloric = maximum_stage_caloric.max(endpoint_caloric);
        substep_defect.add(corrected.report.absolute_enthalpy_defect_j);
        for (k, sum) in mass_defect.iter_mut().enumerate() {
            sum.add(corrected.species.report.absolute_mass_defect_kg[k]);
        }
        for (cell, value) in corrected_heat.iter().enumerate() {
            poll_every(gate, cell)?;
            heat[cell].add(value.sum);
            heat[cell].add(value.correction);
        }
        drop(corrected_heat);
        species = corrected.species.species_mass_kg;
        h = corrected.specific_enthalpy_j_kg;
        current = Some(endpoint);
    }
    let mut stored = EnergySum::default();
    let mut transported = EnergySum::default();
    let mut interval_defect = EnergySum::default();
    let mut mass_change = [EnergySum::default(); N];
    for c in 0..h.len() {
        poll_every(gate, c)?;
        let (delta, tail) = fs_math::eft::two_sum(h[c], -initial_h[c]);
        let mut local = EnergySum::default();
        for dh in [delta, tail] {
            let (value, error) = fs_math::eft::two_prod(carrier.mass_kg[c], dh);
            for term in [value, error] {
                stored.add(term);
                local.add(term);
            }
        }
        for term in [heat[c].sum, heat[c].correction] {
            transported.add(term);
            local.add(-term);
        }
        interval_defect.add(local.value().abs());
        for k in 0..N {
            let (value, error) = fs_math::eft::two_sum(species[c][k], -initial_species[c][k]);
            mass_change[k].add(value);
            mass_change[k].add(error);
        }
    }
    let external_mass = boundary_history
        .as_ref()
        .map(BoundaryHistory::external_mass)
        .unwrap_or([EnergySum::default(); N]);
    let mass_defect = core::array::from_fn(|k| {
        let mut residual = mass_change[k];
        residual.add(-external_mass[k].sum);
        residual.add(-external_mass[k].correction);
        mass_defect[k].value() + residual.value().abs()
    });
    for k in 0..N {
        nonnegative(mass_defect[k])?;
        if mass_defect[k] > controls.solver.maximum_mass_defect_kg[k] {
            return Err(refuse(
                "adaptive interval exceeds its constituent-mass budget",
            ));
        }
    }
    let external = boundary_history
        .as_ref()
        .map(BoundaryHistory::external_energy)
        .unwrap_or_default();
    let mut stored_residual = stored;
    let mut transported_residual = transported;
    for residual in [&mut stored_residual, &mut transported_residual] {
        residual.add(-external.sum);
        residual.add(-external.correction);
    }
    let interval_defect = interval_defect.value()
        + stored_residual.value().abs()
        + transported_residual.value().abs();
    let total_defect =
        substep_defect.value() + interval_defect + original_caloric + endpoint_caloric;
    nonnegative(total_defect)?;
    if total_defect > controls.maximum_enthalpy_defect_j {
        return Err(refuse(
            "adaptive interval exceeds its physical and caloric joule budget",
        ));
    }
    let endpoint =
        current.ok_or_else(|| refuse("adaptive interval has no endpoint preparation"))?;
    let temperature = model.temperature_k(&endpoint).to_vec();
    poll(gate)?;
    Ok((
        Trial {
            species,
            h,
            temperature,
            heat,
            boundary: boundary_history,
            mass_defect,
            stored: stored.value(),
            transported: transported.value(),
            substep_defect: substep_defect.value(),
            interval_defect,
            endpoint_caloric,
            maximum_stage_caloric,
            maximum_source_relaxation,
            maximum_boundary_relaxation,
            total_defect,
        },
        endpoint,
    ))
}

impl MaterialAdvectionStep {
    /// Refresh one original caloric and transport law through positive
    /// exponential midpoint, then refine the complete interval globally.
    ///
    /// Coarse and fine trajectories start from the SAME original inventory.
    /// Each local predictor runs the native frozen solve for dt/2; its source
    /// image drives a native full-dt corrector from the original local state.
    /// Endpoint properties become the next substep's start properties.
    /// The previous fine trajectory is reused as the next coarse comparison.
    ///
    /// Error estimates are |fine-coarse|/3 for smooth second-order evolution,
    /// tested independently in K, absolute J/kg and every species mass fraction.
    /// This is not a certified bound at discontinuities or latent/source corners.
    /// Saturated coarse/fine agreement can also hide a fast smooth thermal
    /// transient. Changed sampled D or face-relative partial enthalpies therefore
    /// require full-outer face relaxation <=1/2 before acceptance. Extra original
    /// face queries are budgeted; common enthalpy shifts do not change contrasts.
    /// Unchanged sampled coefficients retain the stiff frozen-property path.
    /// This guard resolves saturation, not an unsampled constitutive error bound.
    /// No extrapolation, state normalization or invented caloric law is applied.
    /// Original source, mass, storage, positivity and cancellation refusals
    /// propagate; a refused intermediate state is not silently extrapolated.
    /// All resource costs include discarded predictors and resolutions. Physical
    /// defect accounts contain only the accepted correctors and an independent
    /// original-to-final audit, with source inversion defects kept explicit.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn diffuse_species_with_enthalpy_adaptive<
        const N: usize,
        M: SpeciesEnthalpyDiffusionModel<N>,
    >(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        duration_s: f64,
        model: &M,
        controls: AdaptiveSpeciesEnthalpyDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<AdaptiveSpeciesEnthalpyDiffusionStep<N, M::Prepared>, ChtError> {
        self.diffuse_species_with_enthalpy_refined(
            domain,
            species_mass_kg,
            specific_enthalpy_j_kg,
            duration_s,
            model,
            controls,
            gate,
            None,
        )
        .map(|result| result.step)
    }

    /// Refresh original cell and face laws while advancing reservoir exchange.
    ///
    /// Reuses the closed midpoint/refinement algorithm and original positive
    /// reservoir maps. Accepted physical histories retain every face's signed
    /// and gross species/enthalpy transfers; independent transfer tolerances
    /// prevent opposing boundaries from hiding behind stationary endpoint fields.
    /// Predictors and rejected trials consume work, but never external heat.
    ///
    /// The model must bind the supplied boundary descriptors to its original
    /// physical face caloric law. Coefficients are refreshed at each midpoint.
    /// Source-changing stiff boundary rates obey the SAME full-outer relaxation
    /// guard as interior rates. Competing boundary/interior maps additionally
    /// resolve full-inner relaxation, even for fixed coefficients: identical
    /// saturated trials can otherwise allocate all exchange to the first face.
    /// Final mass/energy audits subtract retained actual
    /// external transfers before rounding. Returned heat is already in h.
    /// Accuracy remains estimated; the caller owns external reservoir storage,
    /// thermal interface closure and any pressure/volume evolution.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive<
        const N: usize,
        M: SpeciesEnthalpyDiffusionModel<N>,
    >(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        boundaries: &[SpeciesDiffusionBoundary<N>],
        duration_s: f64,
        model: &M,
        controls: AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionStep<N, M::Prepared>, ChtError> {
        let options = BoundaryOptions {
            boundaries,
            mass_tolerance: controls.absolute_boundary_transfer_tolerance_kg,
            heat_tolerance: controls.absolute_boundary_enthalpy_tolerance_j,
        };
        let result = self.diffuse_species_with_enthalpy_refined(
            domain,
            species_mass_kg,
            specific_enthalpy_j_kg,
            duration_s,
            model,
            controls.integrator,
            gate,
            Some(options),
        )?;
        let (boundary, thermal_boundary) = result
            .boundary
            .ok_or_else(|| refuse("accepted reservoir trajectory lost its external account"))?
            .finish(gate)?;
        poll(gate)?;
        Ok(AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionStep {
            step: result.step,
            boundary,
            thermal_boundary,
            estimated_boundary_transfer_error_kg: result.boundary_mass_error,
            estimated_boundary_enthalpy_error_j: result.boundary_heat_error,
        })
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn diffuse_species_with_enthalpy_refined<
        const N: usize,
        M: SpeciesEnthalpyDiffusionModel<N>,
    >(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        duration_s: f64,
        model: &M,
        controls: AdaptiveSpeciesEnthalpyDiffusionControls<N>,
        gate: &CancelGate,
        boundary: Option<BoundaryOptions<'_, N>>,
    ) -> Result<Refined<N, M::Prepared>, ChtError> {
        poll(gate)?;
        let transport = self.accepted_mass_transport(gate)?;
        let cells = domain.cell_count();
        if N == 0
            || N > MAX_MATERIAL_SPECIES
            || controls.initial_substeps == 0
            || controls.initial_substeps as u64 > (1u64 << 53)
            || controls.maximum_refinements == 0
            || controls.solver.substeps == 0
            || controls.solver.substeps as u64 > (1u64 << 53)
            || self.projection.field.dims() != domain.dims()
            || self.mass_kg.len() != cells
            || species_mass_kg.len() != cells
            || specific_enthalpy_j_kg.len() != cells
        {
            return Err(refuse(
                "adaptive thermal diffusion shapes or work controls are invalid",
            ));
        }
        positive(duration_s)?;
        positive(controls.absolute_temperature_tolerance_k)?;
        positive(controls.absolute_specific_enthalpy_tolerance_j_kg)?;
        nonnegative(controls.relative_mass_fraction_tolerance)?;
        nonnegative(controls.maximum_enthalpy_defect_j)?;
        for tolerance in controls.absolute_mass_fraction_tolerance {
            positive(tolerance)?;
        }
        for allowance in controls.solver.maximum_mass_defect_kg {
            nonnegative(allowance)?;
        }
        nonnegative(controls.solver.maximum_partition_defect_kg)?;
        let mut initial_partition = EnergySum::default();
        for c in 0..cells {
            poll_every(gate, c)?;
            if domain.is_fluid(c) != (transport.row_of[c] != usize::MAX)
                || !specific_enthalpy_j_kg[c].is_finite()
                || (!domain.is_fluid(c)
                    && (specific_enthalpy_j_kg[c] != 0.0 || species_mass_kg[c] != [0.0; N]))
            {
                return Err(refuse(
                    "adaptive thermal state differs from its accepted carrier domain",
                ));
            }
            if domain.is_fluid(c) {
                positive(self.mass_kg[c])?;
                let mut partition = EnergySum::default();
                partition.add(-self.mass_kg[c]);
                for &mass in &species_mass_kg[c] {
                    nonnegative(mass)?;
                    ratio(mass, self.mass_kg[c])?;
                    partition.add(mass);
                }
                initial_partition.add(partition.value().abs());
            }
        }
        nonnegative(initial_partition.value())?;
        if initial_partition.value() > controls.solver.maximum_partition_defect_kg {
            return Err(refuse(
                "adaptive initial species do not partition the accepted carrier",
            ));
        }
        let extra = if let Some(boundary) = boundary {
            boundary.validate()?;
            boundary
                .retained_bytes()
                .ok_or_else(|| refuse("adaptive reservoir history storage extent overflows"))?
        } else {
            0
        };
        let remaining = controls
            .solver
            .maximum_storage_bytes
            .checked_sub(extra)
            .ok_or_else(|| refuse("adaptive reservoir history exceeds storage budget"))?;
        let storage = solver_storage(model, cells, remaining)?;
        let mut work = Work {
            pair_updates: 0,
            component_updates: 0,
            source_preparations: Cell::new(0),
            face_queries: Cell::new(0),
            maximum_sources: controls.maximum_source_evaluations,
        };
        work.source(true)?;
        let original =
            model.prepare(domain, self, species_mass_kg, specific_enthalpy_j_kg, gate)?;
        let original_caloric = prepared(
            model,
            &original,
            domain,
            controls.maximum_enthalpy_defect_j,
            gate,
        )?;
        let mut substeps = controls.initial_substeps;
        let (mut coarse, coarse_prepared) = integrate(
            self,
            domain,
            model,
            &original,
            original_caloric,
            species_mass_kg,
            specific_enthalpy_j_kg,
            duration_s,
            substeps,
            controls,
            storage,
            &mut work,
            gate,
            boundary,
        )?;
        // Comparison needs the copied T/h/Y, not another live source image.
        // Drop it before a fine solve to preserve the three-preparation bound.
        drop(coarse_prepared);
        for refinement in 1..=controls.maximum_refinements {
            poll(gate)?;
            substeps = substeps
                .checked_mul(2)
                .ok_or_else(|| refuse("adaptive thermal refinement count overflows"))?;
            let (fine, endpoint) = integrate(
                self,
                domain,
                model,
                &original,
                original_caloric,
                species_mass_kg,
                specific_enthalpy_j_kg,
                duration_s,
                substeps,
                controls,
                storage,
                &mut work,
                gate,
                boundary,
            )?;
            let mut temperature_error: f64 = 0.0;
            let mut h_error: f64 = 0.0;
            let mut species_error = [0.0_f64; N];
            let mut maximum_error_ratio: f64 = 0.0;
            for c in 0..cells {
                poll_every(gate, c)?;
                if !domain.is_fluid(c) {
                    continue;
                }
                let t = (fine.temperature[c] - coarse.temperature[c]).abs() / 3.0;
                // The physical storage budget may intentionally admit rounding
                // below a large formation datum. It must not erase the separate
                // time-error estimate: compare retained transported heat too.
                let mut heat_delta = fine.heat[c];
                heat_delta.add(-coarse.heat[c].sum);
                heat_delta.add(-coarse.heat[c].correction);
                let heat_delta = heat_delta.value().abs();
                nonnegative(heat_delta)?;
                let h = (fine.h[c] - coarse.h[c])
                    .abs()
                    .max(heat_delta / self.mass_kg[c])
                    / 3.0;
                temperature_error = temperature_error.max(t);
                h_error = h_error.max(h);
                maximum_error_ratio = maximum_error_ratio
                    .max(t / controls.absolute_temperature_tolerance_k)
                    .max(h / controls.absolute_specific_enthalpy_tolerance_j_kg);
                for (k, error) in species_error.iter_mut().enumerate() {
                    let value =
                        (fine.species[c][k] - coarse.species[c][k]).abs() / self.mass_kg[c] / 3.0;
                    let scale = (species_mass_kg[c][k] / self.mass_kg[c])
                        .max(fine.species[c][k] / self.mass_kg[c]);
                    let tolerance = positive(
                        controls.absolute_mass_fraction_tolerance[k]
                            + controls.relative_mass_fraction_tolerance * scale,
                    )?;
                    *error = error.max(value);
                    maximum_error_ratio = maximum_error_ratio.max(value / tolerance);
                }
            }
            let (boundary_mass_error, boundary_heat_error, boundary_ratio) =
                match (&coarse.boundary, &fine.boundary, boundary) {
                    (Some(coarse), Some(fine), Some(options)) => {
                        coarse.compare(fine, options, gate)?
                    }
                    (None, None, None) => ([0.0; N], 0.0, 0.0),
                    _ => {
                        return Err(refuse(
                            "adaptive reservoir comparison lost its physical history",
                        ));
                    }
                };
            maximum_error_ratio = maximum_error_ratio.max(boundary_ratio);
            if maximum_error_ratio <= 1.0
                && fine.maximum_source_relaxation <= 0.5
                && fine.maximum_boundary_relaxation <= 0.5
            {
                let changes = fine.heat.iter().map(|sum| sum.value()).collect();
                let report = AdaptiveSpeciesEnthalpyDiffusionReport {
                    duration_s,
                    accepted_substeps: substeps,
                    refinements: refinement,
                    total_pair_updates: work.pair_updates,
                    total_component_updates: work.component_updates,
                    source_preparations: work.source_preparations.get(),
                    face_source_evaluations: work.face_queries.get(),
                    estimated_temperature_error_k: temperature_error,
                    estimated_specific_enthalpy_error_j_kg: h_error,
                    estimated_mass_fraction_error: species_error,
                    maximum_error_ratio,
                    maximum_changing_source_relaxation: fine.maximum_source_relaxation,
                    maximum_competing_boundary_relaxation: fine.maximum_boundary_relaxation,
                    absolute_mass_defect_kg: fine.mass_defect,
                    stored_enthalpy_change_j: fine.stored,
                    transported_enthalpy_change_j: fine.transported,
                    accepted_substep_enthalpy_defect_j: fine.substep_defect,
                    interval_enthalpy_defect_j: fine.interval_defect,
                    initial_caloric_defect_j: original_caloric,
                    endpoint_caloric_defect_j: fine.endpoint_caloric,
                    maximum_stage_caloric_defect_j: fine.maximum_stage_caloric,
                    absolute_enthalpy_defect_j: fine.total_defect,
                };
                poll(gate)?;
                return Ok(Refined {
                    step: AdaptiveSpeciesEnthalpyDiffusionStep {
                        state: endpoint,
                        species_mass_kg: fine.species,
                        specific_enthalpy_j_kg: fine.h,
                        temperature_k: fine.temperature,
                        diffusive_enthalpy_change_j: changes,
                        report,
                    },
                    boundary: fine.boundary,
                    boundary_mass_error,
                    boundary_heat_error,
                });
            }
            drop(endpoint);
            coarse = fine;
        }
        Err(refuse(
            "adaptive thermal diffusion accuracy, source relaxation or boundary competition is unresolved at its refinement limit",
        ))
    }
}
