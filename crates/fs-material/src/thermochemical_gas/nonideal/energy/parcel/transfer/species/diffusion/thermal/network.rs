//! Nonisothermal PR faces on the original simultaneous and adaptive owners.
//! Species enthalpy, reciprocal heat and donor-to-face correction enter ONE
//! conservative native reconstruction; no endpoint is reheated afterward.
use super::super::network::{
    MolecularFace, RealGasDiffusionNetworkCandidate, RealGasDiffusionNetworkLink,
    RealGasDiffusionNetworkStep, step_counted_candidate,
};
use super::{
    GasDiffusionGeometry, GasEnergyError, IsochoricRealGasParcel, RealGasThermodiffusionFlux,
    RealGasThermodiffusionModel, RealGasThermodiffusionSolve, RealGasThermodiffusionWorkspace,
    SpeciesDiffusionFlux, SpeciesId, positive,
};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::transfer::species::{
        diffusion::network::{
            DiffusionNetworkSolve,
            adaptive::{
                Candidate, DiffusionNetworkAdvanceCandidate, DiffusionNetworkAdvanceControls,
                DiffusionNetworkAdvanceError, DiffusionNetworkAdvanceReport, Transfer,
                TrialFailure, advance_with_kernel, thermal::integral_error,
            },
        },
        network::SpeciesNetworkHeatTransfer,
    },
};
use std::collections::BTreeMap;

pub use super::super::network::{
    RealGasDiffusionNetworkError as RealGasThermodiffusionNetworkError,
    RealGasDiffusionNetworkSolve as RealGasThermodiffusionNetworkSolve,
};
pub type RealGasThermodiffusionLink<'a> =
    RealGasDiffusionNetworkLink<'a, RealGasThermodiffusionModel<'a>>;
pub type RealGasThermodiffusionNetworkStep =
    RealGasDiffusionNetworkStep<RealGasThermodiffusionFlux>;
pub type RealGasThermodiffusionNetworkCandidate =
    RealGasDiffusionNetworkCandidate<RealGasThermodiffusionFlux>;
pub type RealGasThermodiffusionNetworkAdvanceCandidate =
    DiffusionNetworkAdvanceCandidate<RealGasThermodiffusionNetworkStep, IsochoricRealGasParcel>;
pub type RealGasThermodiffusionNetworkAdvanceReport =
    DiffusionNetworkAdvanceReport<RealGasThermodiffusionNetworkStep>;

/// Endpoint accuracy plus independent per-face signed/gross species and heat
/// accuracy. Opposing links or flow reversals cannot hide transferred amounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasThermodiffusionNetworkAdvanceControls {
    pub temporal: DiffusionNetworkAdvanceControls,
    /// All actual original face calls, including failed/rejected partial probes.
    pub maximum_face_evaluations: usize,
    pub absolute_species_transfer_tolerance_kg: f64,
    /// Additional heat only; original donor species enthalpy remains separate.
    pub absolute_additional_heat_tolerance_j: f64,
}

impl MolecularFace for RealGasThermodiffusionModel<'_> {
    type Flux = RealGasThermodiffusionFlux;
    fn species(&self) -> &[SpeciesId] {
        RealGasThermodiffusionModel::species(self)
    }
    fn workspace(&self) -> Result<RealGasThermodiffusionWorkspace, GasEnergyError> {
        RealGasThermodiffusionModel::workspace(self)
    }
    fn sample(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        geometry: GasDiffusionGeometry,
        controls: RealGasThermodiffusionSolve,
        work: &mut RealGasThermodiffusionWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Flux, GasEnergyError> {
        self.face_flux(left, right, geometry, controls, work, cancelled)
    }
    fn admit(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        controls: RealGasThermodiffusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        self.admit_endpoints(left, right, controls, cancelled)
    }
    fn species_flux(flux: &Self::Flux) -> &[SpeciesDiffusionFlux] {
        &flux.species
    }
    fn heat_rate(flux: &Self::Flux) -> f64 {
        flux.additional_heat_rate_left_to_right_w
    }
    fn frame_residual(flux: &Self::Flux) -> f64 {
        flux.mass_frame_residual_kg_m2_s
    }
}

/// Explicit first-order interval using every ORIGINAL face snapshot. Thermal
/// faces already include ordinary molecular transport and Fourier conduction;
/// do not add duplicate faces for those mechanisms on the same exposure.
/// Actual source support and whole-network species/U/S are admitted together.
pub fn step_candidate(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasThermodiffusionLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: RealGasThermodiffusionNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasThermodiffusionNetworkCandidate, RealGasThermodiffusionNetworkError> {
    step_counted_candidate(
        parcels,
        links,
        duration_s,
        endpoints,
        controls,
        &mut BTreeMap::new(),
        &mut 0,
        usize::MAX,
        cancelled,
    )
}

struct ThermalCandidate {
    candidate: RealGasThermodiffusionNetworkCandidate,
    controls: RealGasThermodiffusionNetworkAdvanceControls,
}
impl Candidate for ThermalCandidate {
    type Parcel = IsochoricRealGasParcel;
    type Receipt = RealGasThermodiffusionNetworkStep;
    fn parcels(&self) -> &BTreeMap<u64, Self::Parcel> {
        self.candidate.parcels()
    }
    fn entropy_change(&self) -> f64 {
        self.candidate.step().entropy_change_j_k
    }
    fn energy_residual(&self) -> f64 {
        self.candidate.step().exchange.total_energy_residual_j
    }
    fn transfers(&self) -> impl Iterator<Item = Transfer<'_>> {
        self.candidate
            .step()
            .exchange
            .fluxes
            .iter()
            .map(|f| Transfer {
                donor: f.transfer.donor,
                receiver: f.transfer.receiver,
                species: &f.transfer.species,
                mass_kg: f.transfer.mass_kg,
                enthalpy_j: f.carried_enthalpy_j,
            })
    }
    fn heat(&self) -> &[SpeciesNetworkHeatTransfer] {
        &self.candidate.step().exchange.heat
    }
    fn transfer_error_ratio(
        &self,
        first: &Self,
        second: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<f64>, DiffusionNetworkAdvanceError> {
        let stages = [
            self.candidate.step(),
            first.candidate.step(),
            second.candidate.step(),
        ];
        let mut maximum = 0.0_f64;
        if stages
            .iter()
            .any(|s| s.links.len() != stages[0].links.len())
        {
            return Err(DiffusionNetworkAdvanceError::Invalid(
                "real-gas thermal link identity changed",
            ));
        }
        for i in 0..stages[0].links.len() {
            if cancelled() {
                return Err(DiffusionNetworkAdvanceError::Cancelled);
            }
            let links = stages.map(|s| &s.links[i]);
            if links.iter().any(|l| {
                l.id != links[0].id
                    || l.parcels != links[0].parcels
                    || l.flux.is_some() != links[0].flux.is_some()
            }) {
                return Err(DiffusionNetworkAdvanceError::Invalid(
                    "real-gas thermal face support changed",
                ));
            }
            if links[0].flux.is_none() {
                continue;
            }
            let flux = links.map(|l| l.flux.as_ref().unwrap());
            if flux
                .iter()
                .any(|f| f.species.len() != flux[0].species.len())
            {
                return Err(DiffusionNetworkAdvanceError::Invalid(
                    "real-gas thermal species axis changed",
                ));
            }
            for k in 0..flux[0].species.len() {
                if cancelled() {
                    return Err(DiffusionNetworkAdvanceError::Cancelled);
                }
                if flux
                    .iter()
                    .any(|f| f.species[k].species != flux[0].species[k].species)
                {
                    return Err(DiffusionNetworkAdvanceError::Invalid(
                        "real-gas thermal species order changed",
                    ));
                }
                let amounts = std::array::from_fn::<_, 3, _>(|j| {
                    flux[j].species[k].mass_rate_left_to_right_kg_s * stages[j].duration_s
                });
                let error = integral_error(amounts[0], [amounts[1], amounts[2]]).ok_or(
                    DiffusionNetworkAdvanceError::Invalid(
                        "unrepresentable real-gas thermal mass integral",
                    ),
                )?;
                maximum = maximum.max(error / self.controls.absolute_species_transfer_tolerance_kg);
            }
            let heat = std::array::from_fn::<_, 3, _>(|j| {
                flux[j].additional_heat_rate_left_to_right_w * stages[j].duration_s
            });
            let error = integral_error(heat[0], [heat[1], heat[2]]).ok_or(
                DiffusionNetworkAdvanceError::Invalid(
                    "unrepresentable real-gas thermal heat integral",
                ),
            )?;
            maximum = maximum.max(error / self.controls.absolute_additional_heat_tolerance_j);
        }
        Ok(Some(maximum))
    }
    fn into_parts(self) -> (BTreeMap<u64, Self::Parcel>, Self::Receipt) {
        self.candidate.into_parts()
    }
}

/// Complete a source-refreshed nonisothermal window on the existing temporal
/// controller. Accepted fine steps alone contribute mass, donor enthalpy and
/// additional heat. Matrix scratch and actual face-call accounting span every
/// rejected/failed probe. The original whole-window PR audit and fixed physical
/// budgets remain distinct from local endpoint and per-face transfer estimates.
/// No implicit stability, phase-stability or global/spatial error claim follows.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    original: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasThermodiffusionLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: RealGasThermodiffusionNetworkSolve,
    controls: RealGasThermodiffusionNetworkAdvanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasThermodiffusionNetworkAdvanceCandidate, DiffusionNetworkAdvanceError> {
    if cancelled() {
        return Err(DiffusionNetworkAdvanceError::Cancelled);
    }
    if ![
        physical.exchange.maximum_species_fraction,
        physical.exchange.mass_tolerance_kg,
        physical.exchange.energy_tolerance_j,
        physical.entropy_tolerance_j_k,
        controls.absolute_species_transfer_tolerance_kg,
        controls.absolute_additional_heat_tolerance_j,
    ]
    .into_iter()
    .all(positive)
        || physical.exchange.maximum_species_fraction >= 1.0
        || !(1..=1_000_000).contains(&controls.maximum_face_evaluations)
    {
        return Err(DiffusionNetworkAdvanceError::Invalid(
            "positive native PR thermal budgets and bounded work required",
        ));
    }
    let mut workspaces = BTreeMap::new();
    advance_with_kernel(
        original,
        links.len(),
        duration_s,
        endpoints,
        DiffusionNetworkSolve {
            exchange: physical.exchange,
            entropy_tolerance_j_k: physical.entropy_tolerance_j_k,
        },
        controls.temporal,
        None,
        Some(controls.maximum_face_evaluations),
        |state, dt, work, cancelled| {
            step_counted_candidate(
                state,
                links,
                dt,
                endpoints,
                physical,
                &mut workspaces,
                &mut work.face_evaluations,
                controls.maximum_face_evaluations,
                cancelled,
            )
            .map(|candidate| ThermalCandidate {
                candidate,
                controls,
            })
            .map_err(|error| match error {
                RealGasThermodiffusionNetworkError::FaceWorkLimit { .. } => TrialFailure::FaceWork,
                other => TrialFailure::RealGas(other),
            })
        },
        cancelled,
    )
}
