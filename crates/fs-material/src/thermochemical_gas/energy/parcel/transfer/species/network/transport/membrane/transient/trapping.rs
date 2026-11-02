//! Original-source immobile Langmuir traps on the mobile Henry membrane owner.
//!
//! The adsorbent supplies its existing capacity, calorically consistent affinity
//! and desorption rate. Dry skeleton density is explicit per membrane volume;
//! this binding does not infer porosity, swelling, mobile diffusivity or heat.
//! The fixed-temperature thermostat is external. Actual old/new cell and layer
//! surface pressure ranges must remain within BOTH original source selections.
//! Distinct mobile-polymer and trap material cards declare a composite layer;
//! their complete typed physical contexts must match at each admitted pressure.
//! This pure-permeant bridge does not infer the meaning of mixture/composition
//! or concentration axes; sources carrying those coordinates require an owner
//! that explicitly maps and evolves them. Other fixed conditions stay matched.

use super::slab::finite::MolarPressureCapacity;
use super::slab::trapping::{
    FiniteTrappingMembraneDiffusionState, LangmuirMembraneTrap, TrappingMembraneCandidate,
    TrappingMembraneDiffusionState, TrappingMembraneError, TrappingMembraneSolve,
    TrappingMembraneStep, TrappingMembraneWorkspace,
};
use super::*;
use crate::adsorption::{AdsorptionError, LangmuirAdsorbent, LangmuirProperties};
use fs_qty::semantic::QuantityKind;

/// Accuracy-controlled windows with original-source admission on every probe.
pub mod adaptive;
/// Original finite gas parcels coupled to both membrane inventories.
pub mod reservoir;

#[derive(Clone, Debug)]
pub struct LangmuirMembraneTrapSource {
    pub adsorbent: Arc<LangmuirAdsorbent>,
    /// kg dry adsorbent skeleton per m3 of the complete membrane layer.
    pub dry_skeleton_density_kg_m3: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ResolvedTrappingMembraneError {
    Invalid(&'static str),
    Membrane(ResolvedMembraneDiffusionError),
    Trapping(TrappingMembraneError),
    Adsorption(AdsorptionError),
    Cancelled,
}
impl core::fmt::Display for ResolvedTrappingMembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound kinetic membrane traps: {self:?}")
    }
}
impl std::error::Error for ResolvedTrappingMembraneError {}
impl From<ResolvedMembraneDiffusionError> for ResolvedTrappingMembraneError {
    fn from(error: ResolvedMembraneDiffusionError) -> Self {
        if error == ResolvedMembraneDiffusionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Membrane(error)
        }
    }
}
impl From<MembraneDiffusionError> for ResolvedTrappingMembraneError {
    fn from(error: MembraneDiffusionError) -> Self {
        ResolvedMembraneDiffusionError::from(error).into()
    }
}
impl From<TrappingMembraneError> for ResolvedTrappingMembraneError {
    fn from(error: TrappingMembraneError) -> Self {
        if error == TrappingMembraneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Trapping(error)
        }
    }
}
impl From<AdsorptionError> for ResolvedTrappingMembraneError {
    fn from(error: AdsorptionError) -> Self {
        if error == AdsorptionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Adsorption(error)
        }
    }
}
type Result<T> = core::result::Result<T, ResolvedTrappingMembraneError>;
use ResolvedTrappingMembraneError as Error;

#[derive(Clone, Debug)]
pub struct ResolvedTrappingMembraneDiffusion<'a> {
    bindings: Arc<Bindings<'a>>,
    trap_bindings: Arc<[Option<LangmuirMembraneTrapSource>]>,
    state: TrappingMembraneDiffusionState,
    support: MembraneDiffusionSourceSupport,
    trap_support: Vec<Option<[LangmuirProperties; 2]>>,
}
#[derive(Clone, Debug)]
pub struct ResolvedTrappingMembraneCandidate<'a> {
    state: ResolvedTrappingMembraneDiffusion<'a>,
    step: TrappingMembraneStep,
}
impl<'a> ResolvedTrappingMembraneCandidate<'a> {
    pub const fn state(&self) -> &ResolvedTrappingMembraneDiffusion<'a> {
        &self.state
    }
    pub const fn step(&self) -> &TrappingMembraneStep {
        &self.step
    }
    pub fn into_parts(self) -> (ResolvedTrappingMembraneDiffusion<'a>, TrappingMembraneStep) {
        (self.state, self.step)
    }
}

fn admit_traps(
    species: &SpeciesId,
    temperature: f64,
    traps: &[Option<LangmuirMembraneTrapSource>],
    support: &MembraneDiffusionSourceSupport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(
    Vec<Option<LangmuirMembraneTrap>>,
    Vec<Option<[LangmuirProperties; 2]>>,
)> {
    if traps.len() != support.layers.len() {
        return Err(Error::Invalid(
            "one original trap source slot per membrane layer required",
        ));
    }
    let mut coefficients = Vec::with_capacity(traps.len());
    let mut receipts = Vec::with_capacity(traps.len());
    for (binding, layer) in traps.iter().zip(&support.layers) {
        poll(cancelled)?;
        let Some(binding) = binding else {
            coefficients.push(None);
            receipts.push(None);
            continue;
        };
        if binding.adsorbent.species() != species
            || !binding.dry_skeleton_density_kg_m3.is_finite()
            || binding.dry_skeleton_density_kg_m3 <= 0.0
        {
            return Err(Error::Invalid(
                "trap source must name the same permeant and positive explicit dry skeleton density",
            ));
        }
        // Reuse the original mobile owner's actual cell/ordered-face envelope;
        // no independent pressure interpolation or source selection is made.
        let low = layer[0].state_coordinate("pressure").ok_or(Error::Invalid(
            "membrane support lacks its actual pressure coordinate",
        ))?;
        let high = layer[1].state_coordinate("pressure").ok_or(Error::Invalid(
            "membrane support lacks its actual pressure coordinate",
        ))?;
        let a = binding.adsorbent.properties(temperature, low)?;
        poll(cancelled)?;
        let b = if high == low {
            a.clone()
        } else {
            binding.adsorbent.properties(temperature, high)?
        };
        for (mobile, trap) in layer.iter().zip([&a.source, &b.source]) {
            if mobile.query_point() != trap.query_point()
                || mobile.axis_quantities() != trap.axis_quantities()
            {
                return Err(Error::Invalid(
                    "mobile polymer and trap sources must share the same complete typed physical context",
                ));
            }
            if mobile.axis_quantities().values().any(|q| {
                q.semantic_type().is_some_and(|q| {
                    matches!(
                        q.kind(),
                        QuantityKind::Composition(_)
                            | QuantityKind::MassConcentration
                            | QuantityKind::AmountConcentration
                    )
                })
            }) {
                return Err(Error::Invalid(
                    "pure-permeant trap bridge cannot freeze unowned composition or concentration coordinates",
                ));
            }
        }
        if a.affinity_pa_inverse != b.affinity_pa_inverse
            || a.desorption_rate_s_inverse != b.desorption_rate_s_inverse
        {
            return Err(Error::Invalid(
                "fixed-temperature trap coefficients changed across the admitted pressure interval",
            ));
        }
        coefficients.push(Some(LangmuirMembraneTrap::try_new(
            binding.dry_skeleton_density_kg_m3 * binding.adsorbent.capacity_mol_kg_dry(),
            a.affinity_pa_inverse,
            a.desorption_rate_s_inverse,
        )?));
        receipts.push(Some([a, b]));
    }
    Ok((coefficients, receipts))
}

// Both boundary owners admit the exact same old/new mobile/face envelopes and
// original trap sources. A joint finite-reservoir step therefore cannot replace
// or extend a constitutive law when its reservoirs deplete or pressurize.
fn admit_step(
    bindings: &Bindings<'_>,
    traps: &[Option<LangmuirMembraneTrapSource>],
    before: &TrappingMembraneDiffusionState,
    before_faces: &[[Option<f64>; 2]],
    after: &TrappingMembraneDiffusionState,
    step: &TrappingMembraneStep,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(
    MembraneDiffusionSourceSupport,
    Vec<Option<[LangmuirProperties; 2]>>,
)> {
    let support = bindings.admit(
        [
            (before.mobile(), before_faces),
            (
                after.mobile(),
                step.transport.face_partial_pressures_pa.as_slice(),
            ),
        ],
        cancelled,
    )?;
    let (coefficients, trap_support) = admit_traps(
        &bindings.species,
        bindings.temperature_k,
        traps,
        &support,
        cancelled,
    )?;
    if coefficients.as_slice() != after.traps() {
        return Err(Error::Invalid(
            "candidate changed its original trap coefficients or density binding",
        ));
    }
    poll(cancelled)?;
    Ok((support, trap_support))
}

impl<'a> ResolvedTrappingMembraneDiffusion<'a> {
    pub fn new(
        mobile: ResolvedMembraneDiffusion<'a>,
        traps: &[Option<LangmuirMembraneTrapSource>],
        initial_trapped_moles: &[f64],
        initial_reservoirs_pa: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let flux = mobile
            .state
            .sample_fluxes(initial_reservoirs_pa, &mut cancelled)?;
        let sample = (&mobile.state, flux.face_partial_pressures_pa.as_slice());
        let support = mobile.bindings.admit([sample, sample], &mut cancelled)?;
        let (coefficients, trap_support) = admit_traps(
            &mobile.bindings.species,
            mobile.bindings.temperature_k,
            traps,
            &support,
            &mut cancelled,
        )?;
        let state = TrappingMembraneDiffusionState::new(
            mobile.state,
            &coefficients,
            initial_trapped_moles,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            bindings: mobile.bindings,
            trap_bindings: Arc::from(traps),
            state,
            support,
            trap_support,
        })
    }
    pub const fn state(&self) -> &TrappingMembraneDiffusionState {
        &self.state
    }
    pub const fn source_support(&self) -> &MembraneDiffusionSourceSupport {
        &self.support
    }
    pub fn trap_source_support(&self) -> &[Option<[LangmuirProperties; 2]>] {
        &self.trap_support
    }
    pub fn species(&self) -> &SpeciesId {
        &self.bindings.species
    }
    pub fn temperature_k(&self) -> f64 {
        self.bindings.temperature_k
    }
    pub fn workspace(&self) -> Result<TrappingMembraneWorkspace> {
        Ok(self.state.workspace()?)
    }

    /// Simultaneous mobile diffusion/site kinetics, then original source
    /// admission of both endpoint inventories and their ordered face traces.
    /// Every refusal leaves both stores, the clock and all receipts unchanged.
    pub fn step_candidate(
        &self,
        reservoirs_pa: [f64; 2],
        duration_s: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedTrappingMembraneCandidate<'a>> {
        poll(&mut cancelled)?;
        let before = self
            .state
            .mobile()
            .sample_fluxes(reservoirs_pa, &mut cancelled)?;
        let candidate: TrappingMembraneCandidate =
            self.state
                .step_candidate(reservoirs_pa, duration_s, controls, work, &mut cancelled)?;
        let (state, step) = candidate.into_parts();
        let (support, trap_support) = admit_step(
            &self.bindings,
            &self.trap_bindings,
            &self.state,
            &before.face_partial_pressures_pa,
            &state,
            &step,
            &mut cancelled,
        )?;
        Ok(ResolvedTrappingMembraneCandidate {
            state: Self {
                bindings: self.bindings.clone(),
                trap_bindings: self.trap_bindings.clone(),
                state,
                support,
                trap_support,
            },
            step,
        })
    }
}

/// Source-qualified mobile/trapped storage between two finite isothermal baths.
///
/// Each bath supplies its explicit whole-reservoir dn/dp [mol/Pa]. All four
/// inventories participate in the existing native implicit solve; source
/// support follows their actual old/new pressures. This owner does not infer
/// reservoir gas composition, an EOS, caloric storage or thermostat work.
#[derive(Clone, Debug)]
pub struct ResolvedFiniteTrappingMembraneDiffusion<'a> {
    bindings: Arc<Bindings<'a>>,
    trap_bindings: Arc<[Option<LangmuirMembraneTrapSource>]>,
    state: FiniteTrappingMembraneDiffusionState,
    support: MembraneDiffusionSourceSupport,
    trap_support: Vec<Option<[LangmuirProperties; 2]>>,
}

#[derive(Clone, Debug)]
pub struct ResolvedFiniteTrappingMembraneCandidate<'a> {
    state: ResolvedFiniteTrappingMembraneDiffusion<'a>,
    step: TrappingMembraneStep,
}
impl<'a> ResolvedFiniteTrappingMembraneCandidate<'a> {
    pub const fn state(&self) -> &ResolvedFiniteTrappingMembraneDiffusion<'a> {
        &self.state
    }
    pub const fn step(&self) -> &TrappingMembraneStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedFiniteTrappingMembraneDiffusion<'a>,
        TrappingMembraneStep,
    ) {
        (self.state, self.step)
    }
}

impl<'a> ResolvedFiniteTrappingMembraneDiffusion<'a> {
    /// Consume the original source-bound film and bind its actual initial
    /// finite reservoir pressures. Construction rechecks every original source
    /// because these baths need not match the prior prescribed boundaries.
    pub fn new(
        membrane: ResolvedTrappingMembraneDiffusion<'a>,
        capacities: [MolarPressureCapacity; 2],
        reservoir_pressures_pa: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let flux = membrane
            .state
            .mobile()
            .sample_fluxes(reservoir_pressures_pa, &mut cancelled)?;
        let sample = (
            membrane.state.mobile(),
            flux.face_partial_pressures_pa.as_slice(),
        );
        let support = membrane.bindings.admit([sample, sample], &mut cancelled)?;
        let (coefficients, trap_support) = admit_traps(
            &membrane.bindings.species,
            membrane.bindings.temperature_k,
            &membrane.trap_bindings,
            &support,
            &mut cancelled,
        )?;
        if coefficients.as_slice() != membrane.state.traps() {
            return Err(Error::Invalid(
                "finite reservoir binding changed original trap coefficients or density",
            ));
        }
        let state = FiniteTrappingMembraneDiffusionState::new(
            membrane.state,
            capacities,
            reservoir_pressures_pa,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            bindings: membrane.bindings,
            trap_bindings: membrane.trap_bindings,
            state,
            support,
            trap_support,
        })
    }
    pub const fn state(&self) -> &FiniteTrappingMembraneDiffusionState {
        &self.state
    }
    pub const fn membrane(&self) -> &TrappingMembraneDiffusionState {
        self.state.membrane()
    }
    pub const fn capacities(&self) -> [MolarPressureCapacity; 2] {
        self.state.capacities()
    }
    pub const fn reservoir_pressures_pa(&self) -> [f64; 2] {
        self.state.reservoir_pressures_pa()
    }
    pub const fn source_support(&self) -> &MembraneDiffusionSourceSupport {
        &self.support
    }
    pub fn trap_source_support(&self) -> &[Option<[LangmuirProperties; 2]>] {
        &self.trap_support
    }
    pub fn species(&self) -> &SpeciesId {
        &self.bindings.species
    }
    pub fn temperature_k(&self) -> f64 {
        self.bindings.temperature_k
    }
    pub fn workspace(&self) -> Result<TrappingMembraneWorkspace> {
        Ok(self.state.workspace()?)
    }

    /// Propose one simultaneous finite-bath/mobile/site update. Unsupported
    /// old or new cell/face conditions refuse the whole candidate atomically.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedFiniteTrappingMembraneCandidate<'a>> {
        poll(&mut cancelled)?;
        let before = self
            .membrane()
            .mobile()
            .sample_fluxes(self.reservoir_pressures_pa(), &mut cancelled)?;
        let (state, step) = self
            .state
            .step_candidate(duration_s, controls, work, &mut cancelled)?
            .into_parts();
        let (support, trap_support) = admit_step(
            &self.bindings,
            &self.trap_bindings,
            self.membrane(),
            &before.face_partial_pressures_pa,
            state.membrane(),
            &step,
            &mut cancelled,
        )?;
        Ok(ResolvedFiniteTrappingMembraneCandidate {
            state: Self {
                bindings: self.bindings.clone(),
                trap_bindings: self.trap_bindings.clone(),
                state,
                support,
                trap_support,
            },
            step,
        })
    }
}
