//! Simultaneous diffusion and competitive finite-site storage in a rigid matrix.
//!
//! Each species retains its original Henry storage and layered D*S face law.
//! Fluxes are relative to the stationary polymer, not a zero-sum gas-mixture
//! reference frame. All adsorbing channels in a layer consume ONE common site
//! inventory. The existing competitive site owner eliminates its implicit
//! endpoint, including off-diagonal displacement derivatives. Native scalar
//! tridiagonal factors and bounded block Gauss-Seidel solve Newton corrections
//! in O(cells*species) work per sweep, without a dense global Jacobian.
//!
//! This numerical owner does not infer mixture-qualified material sources,
//! gas EOS/calorics, swelling, cross diffusion or heat of sorption. Temperature,
//! geometry, Henry coefficients, kinetics and site capacities remain fixed.

use super::finite::MolarPressureCapacity;
use super::trapping::{LangmuirMembraneTrap, TrappingMembraneStep, TrappingReservoirBalance};
use super::*;
use crate::adsorption::{
    AdsorptionError,
    competitive::{
        MAX_COMPETITIVE_SPECIES,
        sites::{CompetitiveSiteBuffer, backward_euler_sites_into},
    },
};
use fs_thermochem::SpeciesId;

mod solve;

/// One original mobile field, per-layer trap laws and explicit retained loading.
/// Non-None trap capacities must agree exactly across species within each layer;
/// they describe one common population and are never added together.
#[derive(Clone, Debug)]
pub struct CompetitiveMembraneChannel {
    pub species: SpeciesId,
    pub mobile: MembraneDiffusionState,
    pub traps: Vec<Option<LangmuirMembraneTrap>>,
    pub initial_trapped_moles: Vec<f64>,
}
#[derive(Clone, Debug, PartialEq)]
struct Channel {
    species: SpeciesId,
    mobile: MembraneDiffusionState,
    traps: Arc<[Option<LangmuirMembraneTrap>]>,
    trapped: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompetitiveMembraneSolve {
    pub balance: MembraneDiffusionSolve,
    pub maximum_iterations: usize,
    pub maximum_backtracks: usize,
    pub maximum_linear_sweeps: usize,
    /// Every correction must satisfy its original Jacobian residual, relative
    /// to max(current nonlinear residual, caller's amount allowance).
    pub linear_relative_tolerance: f64,
}
impl Default for CompetitiveMembraneSolve {
    fn default() -> Self {
        Self {
            balance: MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: 1e-14,
                relative_amount_tolerance: 1e-10,
            },
            maximum_iterations: 64,
            maximum_backtracks: 32,
            maximum_linear_sweeps: 256,
            linear_relative_tolerance: 0.05,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub enum CompetitiveMembraneError {
    Invalid(&'static str),
    Diffusion(MembraneDiffusionError),
    Sites(AdsorptionError),
    NotConverged {
        iterations: usize,
        residual_ratio: f64,
    },
    LinearWork {
        sweeps: usize,
        residual_ratio: f64,
    },
    Cancelled,
}
impl core::fmt::Display for CompetitiveMembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "competitive membrane diffusion: {self:?}")
    }
}
impl std::error::Error for CompetitiveMembraneError {}
impl From<MembraneDiffusionError> for CompetitiveMembraneError {
    fn from(e: MembraneDiffusionError) -> Self {
        if e == MembraneDiffusionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Diffusion(e)
        }
    }
}
impl From<TridiagonalError> for CompetitiveMembraneError {
    fn from(e: TridiagonalError) -> Self {
        MembraneDiffusionError::from(e).into()
    }
}
impl From<AdsorptionError> for CompetitiveMembraneError {
    fn from(e: AdsorptionError) -> Self {
        if e == AdsorptionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Sites(e)
        }
    }
}
type Result<T> = core::result::Result<T, CompetitiveMembraneError>;
fn invalid(message: &'static str) -> CompetitiveMembraneError {
    CompetitiveMembraneError::Invalid(message)
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveMembraneDiffusionState {
    channels: Vec<Channel>,
    capacities: Arc<[f64]>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveMembraneStep {
    /// Canonical SpeciesId order. Each row independently audits storage, all
    /// original faces, boundary streams and optional finite bath conservation.
    pub species: Vec<(SpeciesId, TrappingMembraneStep)>,
    pub vacancy_before_mol: Vec<f64>,
    pub vacancy_after_mol: Vec<f64>,
    pub nonlinear_iterations: usize,
    pub linear_sweeps: usize,
    pub residual_evaluations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveMembraneCandidate<S = CompetitiveMembraneDiffusionState> {
    state: S,
    step: CompetitiveMembraneStep,
}
impl<S> CompetitiveMembraneCandidate<S> {
    pub const fn state(&self) -> &S {
        &self.state
    }
    pub const fn step(&self) -> &CompetitiveMembraneStep {
        &self.step
    }
    pub fn into_parts(self) -> (S, CompetitiveMembraneStep) {
        (self.state, self.step)
    }
}

/// Prepared matrix and kinetics storage. No nonlinear evaluation, correction
/// sweep or line-search trial allocates; only returned physical candidates and
/// reports own new arrays. The explicit byte ceiling covers retained buffers.
#[derive(Debug)]
pub struct CompetitiveMembraneWorkspace {
    rows: usize,
    species: usize,
    factors: Vec<TridiagonalWorkspace>,
    current: Vec<f64>,
    trial: Vec<f64>,
    correction: Vec<f64>,
    diagonal: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    rhs: Vec<f64>,
    scale: Vec<f64>,
    linear_allowance: Vec<f64>,
    occupied: Vec<f64>,
    adsorbed: Vec<f64>,
    desorbed: Vec<f64>,
    site_residual: Vec<f64>,
    site_diagonal: Vec<f64>,
    coupling_row: Vec<f64>,
    coupling_column: Vec<f64>,
    cross: Vec<f64>,
    cross_compensation: Vec<f64>,
    linear_rhs: Vec<f64>,
    next_correction: Vec<f64>,
    vacancy: Vec<f64>,
}
fn zeros(n: usize) -> Result<Vec<f64>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| invalid("competitive workspace allocation refused"))?;
    v.resize(n, 0.0);
    Ok(v)
}
impl CompetitiveMembraneWorkspace {
    fn required_bytes(rows: usize, species: usize) -> Option<usize> {
        // 16 channel arrays, four temporary rows, one vacancy row, plus each
        // native factor's five f64 arrays, bool pivots and allocation headers.
        rows.checked_mul(species)?
            .checked_mul(22 * size_of::<f64>() + size_of::<bool>())?
            .checked_add(rows.checked_mul(5 * size_of::<f64>())?)?
            .checked_add(species.checked_mul(size_of::<TridiagonalWorkspace>())?)?
            .checked_add(size_of::<Self>())
    }
    fn new(rows: usize, species: usize, maximum_bytes: usize) -> Result<Self> {
        let bytes = Self::required_bytes(rows, species)
            .ok_or(invalid("competitive workspace extent overflow"))?;
        if bytes > maximum_bytes || bytes > isize::MAX as usize {
            return Err(invalid("competitive workspace byte budget exceeded"));
        }
        let n = rows * species;
        let mut factors = Vec::new();
        factors
            .try_reserve_exact(species)
            .map_err(|_| invalid("competitive factor allocation refused"))?;
        for _ in 0..species {
            factors.push(TridiagonalWorkspace::new(rows)?);
        }
        Ok(Self {
            rows,
            species,
            factors,
            current: zeros(n)?,
            trial: zeros(n)?,
            correction: zeros(n)?,
            diagonal: zeros(n)?,
            lower: zeros(n)?,
            upper: zeros(n)?,
            rhs: zeros(n)?,
            scale: zeros(n)?,
            linear_allowance: zeros(n)?,
            occupied: zeros(n)?,
            adsorbed: zeros(n)?,
            desorbed: zeros(n)?,
            site_residual: zeros(n)?,
            site_diagonal: zeros(n)?,
            coupling_row: zeros(n)?,
            coupling_column: zeros(n)?,
            cross: zeros(rows)?,
            cross_compensation: zeros(rows)?,
            linear_rhs: zeros(rows)?,
            next_correction: zeros(rows)?,
            vacancy: zeros(rows)?,
        })
    }
}
#[derive(Clone, Copy)]
enum Boundary<'a> {
    Prescribed(&'a [[f64; 2]]),
    Finite {
        pressures: &'a [[f64; 2]],
        capacities: [f64; 2],
    },
}
impl<'a> Boundary<'a> {
    fn pressures(self) -> &'a [[f64; 2]] {
        match self {
            Self::Prescribed(p) | Self::Finite { pressures: p, .. } => p,
        }
    }
    fn offset(self) -> usize {
        usize::from(matches!(self, Self::Finite { .. }))
    }
}

impl CompetitiveMembraneDiffusionState {
    pub fn new(
        mut inputs: Vec<CompetitiveMembraneChannel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if inputs.is_empty() || inputs.len() > MAX_COMPETITIVE_SPECIES {
            return Err(invalid("competitive membrane needs 1..32 explicit species"));
        }
        inputs.sort_by(|a, b| a.species.cmp(&b.species));
        let model = &inputs[0].mobile.model;
        let n = model.cells.len();
        let layers = model.layers.len();
        let time = inputs[0].mobile.elapsed_s;
        let mut capacity_density = vec![0.0; layers];
        for (s, input) in inputs.iter().enumerate() {
            poll(&mut cancelled)?;
            let other = &input.mobile.model;
            if (s > 0 && input.species == inputs[s - 1].species)
                || input.mobile.elapsed_s != time
                || other.area_m2 != model.area_m2
                || other.layer_ranges != model.layer_ranges
                || other.cells.len() != n
                || other.layers.len() != layers
                || other
                    .cells
                    .iter()
                    .zip(&model.cells)
                    .any(|(a, b)| a.width_m != b.width_m || a.layer != b.layer)
                || input.traps.len() != layers
                || input.initial_trapped_moles.len() != n
            {
                return Err(invalid(
                    "unique channels must share the same ordered geometry/time and explicit trap/loading extents",
                ));
            }
            for (i, trap) in input.traps.iter().enumerate() {
                if let Some(t) = trap {
                    let capacity = t.capacity_mol_m3();
                    if capacity_density[i] > 0.0 && capacity_density[i] != capacity {
                        return Err(invalid(
                            "species must share ONE identical layer site capacity",
                        ));
                    }
                    capacity_density[i] = capacity;
                }
            }
            for (cell, &initial) in input.initial_trapped_moles.iter().enumerate() {
                if !nonnegative(initial)
                    || (input.traps[model.cells[cell].layer].is_none() && initial != 0.0)
                {
                    return Err(invalid(
                        "trapped loading needs finite nonnegative amount and an explicit species site law",
                    ));
                }
            }
        }
        let capacities = model
            .cells
            .iter()
            .map(|cell| finite(capacity_density[cell.layer] * (cell.width_m * model.area_m2)))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut sum = Sum::default();
            for input in &inputs {
                sum.add(input.initial_trapped_moles[i])?;
            }
            if sum.value > capacities[i]
                || (capacity_density[model.cells[i].layer] > 0.0 && capacities[i] == 0.0)
            {
                return Err(invalid(
                    "combined species loading exceeds shared capacity or capacity underflowed",
                ));
            }
        }
        let channels = inputs
            .into_iter()
            .map(|i| Channel {
                species: i.species,
                mobile: i.mobile,
                traps: i.traps.into(),
                trapped: i.initial_trapped_moles,
            })
            .collect();
        let state = Self {
            channels,
            capacities: capacities.into(),
        };
        for s in 0..state.channels.len() {
            poll(&mut cancelled)?;
            state.species_amount_mol(s)?;
        }
        poll(&mut cancelled)?;
        Ok(state)
    }
    pub fn species(&self) -> impl ExactSizeIterator<Item = &SpeciesId> {
        self.channels.iter().map(|c| &c.species)
    }
    pub fn mobile(&self, species: usize) -> Option<&MembraneDiffusionState> {
        self.channels.get(species).map(|c| &c.mobile)
    }
    pub fn trapped_moles(&self, species: usize) -> Option<&[f64]> {
        self.channels.get(species).map(|c| c.trapped.as_slice())
    }
    pub fn site_capacities_mol(&self) -> &[f64] {
        &self.capacities
    }
    pub fn elapsed_s(&self) -> f64 {
        self.channels[0].mobile.elapsed_s
    }
    pub fn species_amount_mol(&self, species: usize) -> Result<f64> {
        let c = self
            .channels
            .get(species)
            .ok_or(invalid("species outside competitive membrane roster"))?;
        let mut sum = Sum::default();
        sum.add(c.mobile.stored_amount_mol()?)?;
        for &value in &c.trapped {
            sum.add(value)?;
        }
        Ok(sum.value)
    }
    pub fn workspace(&self, maximum_bytes: usize) -> Result<CompetitiveMembraneWorkspace> {
        CompetitiveMembraneWorkspace::new(self.capacities.len(), self.channels.len(), maximum_bytes)
    }
    /// Boundary arrays follow the canonical order returned by `species()`.
    pub fn step_candidate(
        &self,
        partial_pressures_pa: &[[f64; 2]],
        duration_s: f64,
        controls: CompetitiveMembraneSolve,
        work: &mut CompetitiveMembraneWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompetitiveMembraneCandidate> {
        self.solve(
            Boundary::Prescribed(partial_pressures_pa),
            duration_s,
            controls,
            work,
            cancelled,
        )
    }
}

/// Each bath supplies one fixed whole-reservoir molar pressure capacity. All
/// species partial pressures and both kinds of membrane storage evolve in the
/// same coupled solve. No gas mixture EOS or energy law is inferred here.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteCompetitiveMembraneDiffusionState {
    membrane: CompetitiveMembraneDiffusionState,
    capacities: [MolarPressureCapacity; 2],
    pressures: Vec<[f64; 2]>,
}
impl FiniteCompetitiveMembraneDiffusionState {
    pub fn new(
        membrane: CompetitiveMembraneDiffusionState,
        capacities: [MolarPressureCapacity; 2],
        partial_pressures_pa: Vec<[f64; 2]>,
    ) -> Result<Self> {
        if partial_pressures_pa.len() != membrane.channels.len()
            || partial_pressures_pa
                .iter()
                .flatten()
                .any(|&p| !nonnegative(p))
        {
            return Err(invalid(
                "finite competitive baths need canonical nonnegative partial pressures for every species",
            ));
        }
        let state = Self {
            membrane,
            capacities,
            pressures: partial_pressures_pa,
        };
        for s in 0..state.membrane.channels.len() {
            state.species_amount_mol(s)?;
        }
        Ok(state)
    }
    pub const fn membrane(&self) -> &CompetitiveMembraneDiffusionState {
        &self.membrane
    }
    pub fn reservoir_partial_pressures_pa(&self) -> &[[f64; 2]] {
        &self.pressures
    }
    pub const fn capacities(&self) -> [MolarPressureCapacity; 2] {
        self.capacities
    }
    pub fn species_amount_mol(&self, species: usize) -> Result<f64> {
        let mut sum = Sum::default();
        sum.add(self.membrane.species_amount_mol(species)?)?;
        for side in 0..2 {
            sum.add(finite(
                self.capacities[side].value_si() * self.pressures[species][side],
            )?)?;
        }
        Ok(sum.value)
    }
    pub fn workspace(&self, maximum_bytes: usize) -> Result<CompetitiveMembraneWorkspace> {
        CompetitiveMembraneWorkspace::new(
            self.membrane.capacities.len() + 2,
            self.membrane.channels.len(),
            maximum_bytes,
        )
    }
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: CompetitiveMembraneSolve,
        work: &mut CompetitiveMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompetitiveMembraneCandidate<Self>> {
        let (membrane, step) = self
            .membrane
            .solve(
                Boundary::Finite {
                    pressures: &self.pressures,
                    capacities: self.capacities.map(MolarPressureCapacity::value_si),
                },
                duration_s,
                controls,
                work,
                &mut cancelled,
            )?
            .into_parts();
        let pressures = step
            .species
            .iter()
            .map(|(_, r)| {
                r.finite_reservoirs
                    .as_ref()
                    .expect("finite bath report")
                    .pressures_after_pa
            })
            .collect();
        let state = Self::new(membrane, self.capacities, pressures)?;
        poll(&mut cancelled)?;
        Ok(CompetitiveMembraneCandidate { state, step })
    }
}

#[cfg(test)]
mod tests;
