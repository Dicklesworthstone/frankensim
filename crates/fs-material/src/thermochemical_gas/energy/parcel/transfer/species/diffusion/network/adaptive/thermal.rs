//! Source-conjugate thermal candidates on the original temporal controller.

use super::super::super::thermal::{
    GasThermodiffusionWorkspace,
    network::{
        self, GasThermodiffusionLink, GasThermodiffusionNetworkAdvanceCandidate,
        GasThermodiffusionNetworkAdvanceControls, GasThermodiffusionNetworkCandidate,
        GasThermodiffusionNetworkError, GasThermodiffusionNetworkStep,
    },
};
use super::{
    Candidate, ClosedGasParcel, DiffusionNetworkError, DiffusionNetworkSolve, Error,
    GasEnergyError, GasEnergySolve, GasTransportError, GasTransportSolve, SpeciesNetworkError,
    SpeciesNetworkHeatTransfer, Transfer, TrialFailure, advance_with_kernel, poll,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;
use fs_math::dd::Dd;
use std::collections::BTreeMap;

struct ThermalCandidate {
    candidate: GasThermodiffusionNetworkCandidate,
    controls: GasThermodiffusionNetworkAdvanceControls,
}

impl Candidate for ThermalCandidate {
    type Receipt = GasThermodiffusionNetworkStep;
    fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        self.candidate.parcels()
    }
    fn entropy_change(&self) -> f64 {
        self.candidate.transport().step().entropy_change_j_k
    }
    fn energy_residual(&self) -> f64 {
        self.candidate
            .transport()
            .step()
            .exchange
            .total_energy_residual_j
    }
    fn transfers(&self) -> impl Iterator<Item = Transfer<'_>> {
        self.candidate
            .transport()
            .step()
            .exchange
            .fluxes
            .iter()
            .map(|row| Transfer {
                donor: row.transfer.donor,
                receiver: row.transfer.receiver,
                species: &row.transfer.species,
                mass_kg: row.transfer.mass_kg,
                enthalpy_j: row.carried_enthalpy_j,
            })
    }
    fn heat(&self) -> &[SpeciesNetworkHeatTransfer] {
        &self.candidate.transport().step().exchange.heat
    }
    fn transfer_error_ratio(
        &self,
        first: &Self,
        second: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<f64>, Error> {
        let stages = [&self.candidate, &first.candidate, &second.candidate];
        let durations = stages.map(|stage| stage.transport().step().duration_s);
        let mut maximum = 0.0_f64;
        if stages
            .iter()
            .any(|stage| stage.links().len() != stages[0].links().len())
        {
            return Err(Error::Invalid("thermal comparison changed link identities"));
        }
        for index in 0..stages[0].links().len() {
            poll(cancelled)?;
            let links = stages.map(|stage| &stage.links()[index]);
            if links.iter().any(|link| {
                link.id != links[0].id
                    || link.parcels != links[0].parcels
                    || link.flux.species.len() != links[0].flux.species.len()
            }) {
                return Err(Error::Invalid(
                    "thermal comparison changed link or species identities",
                ));
            }
            for k in 0..links[0].flux.species.len() {
                poll(cancelled)?;
                let rows = links.map(|link| &link.flux.species[k]);
                if rows.iter().any(|row| row.species != rows[0].species) {
                    return Err(Error::Invalid(
                        "thermal comparison changed constituent order",
                    ));
                }
                // Match the actual original transaction's f64 rate*time values.
                // Keep the fine sum's low part until after the coarse subtraction.
                let amounts = core::array::from_fn::<_, 3, _>(|i| {
                    rows[i].mass_rate_left_to_right_kg_s * durations[i]
                });
                let error = integral_error(amounts[0], [amounts[1], amounts[2]]).ok_or(
                    Error::Invalid("unrepresentable thermal mass-transfer comparison"),
                )?;
                maximum = maximum.max(error / self.controls.absolute_species_transfer_tolerance_kg);
            }
            let energies = core::array::from_fn::<_, 3, _>(|i| {
                links[i].flux.additional_heat_rate_left_to_right_w * durations[i]
            });
            let error = integral_error(energies[0], [energies[1], energies[2]]).ok_or(
                Error::Invalid("unrepresentable thermal heat-transfer comparison"),
            )?;
            maximum = maximum.max(error / self.controls.absolute_additional_heat_tolerance_j);
        }
        Ok(Some(maximum))
    }
    fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, Self::Receipt) {
        let (candidate, links) = self.candidate.into_parts();
        let (parcels, transport) = candidate.into_parts();
        (parcels, GasThermodiffusionNetworkStep { transport, links })
    }
}

// An actual source/domain refusal is never reclassified by parsing its text.
// Only the native finite-step physical failures have a retry interpretation.
pub(super) fn retryable(error: &GasThermodiffusionNetworkError) -> bool {
    matches!(
        error,
        GasThermodiffusionNetworkError::Transport(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
                | SpeciesNetworkError::Parcel {
                    error: GasEnergyError::OutsideBracket { .. },
                    ..
                }
        )) | GasThermodiffusionNetworkError::Transport(GasTransportError::Diffusion(
            DiffusionNetworkError::EntropyDecrease { .. }
        ))
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn advance_candidate(
    original: &BTreeMap<u64, ClosedGasParcel>,
    links: &[GasThermodiffusionLink<'_>],
    workspaces: &mut [GasThermodiffusionWorkspace],
    other_paths: GasTransportPaths<'_>,
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: GasTransportSolve,
    controls: GasThermodiffusionNetworkAdvanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasThermodiffusionNetworkAdvanceCandidate, Error> {
    poll(&mut cancelled)?;
    if ![
        physical.entropy_tolerance_j_k,
        physical.exchange.mass_tolerance_kg,
        physical.exchange.energy_tolerance_j,
        physical.exchange.maximum_species_fraction,
        controls.absolute_species_transfer_tolerance_kg,
        controls.absolute_additional_heat_tolerance_j,
    ]
    .into_iter()
    .all(|x| x.is_finite() && x > 0.0)
        || physical.exchange.maximum_species_fraction >= 1.0
    {
        return Err(Error::Invalid(
            "positive original thermal physical and transfer budgets required",
        ));
    }
    advance_with_kernel(
        original,
        &[],
        duration_s,
        endpoints,
        DiffusionNetworkSolve {
            exchange: physical.exchange,
            entropy_tolerance_j_k: physical.entropy_tolerance_j_k,
        },
        controls.temporal,
        None,
        |state, dt, _, cancelled| {
            network::step_candidate(
                state,
                links,
                workspaces,
                other_paths,
                dt,
                endpoints,
                physical,
                cancelled,
            )
            .map(|candidate| ThermalCandidate {
                candidate,
                controls,
            })
            .map_err(TrialFailure::Thermal)
        },
        cancelled,
    )
}

// Largest signed, gross-forward and gross-reverse coarse/fine discrepancy.
// The local method is first order, so the step-doubling divisor is one.
fn integral_error(coarse: f64, fine: [f64; 2]) -> Option<f64> {
    if ![coarse, fine[0], fine[1]].into_iter().all(f64::is_finite) {
        return None;
    }
    let mut maximum = 0.0_f64;
    for [c, a, b] in [
        [coarse, fine[0], fine[1]],
        [coarse.max(0.0), fine[0].max(0.0), fine[1].max(0.0)],
        [(-coarse).max(0.0), (-fine[0]).max(0.0), (-fine[1]).max(0.0)],
    ] {
        let difference = (Dd::from_f64(a) - Dd::from_f64(c)) + Dd::from_f64(b);
        let error = difference.to_f64().abs();
        if !difference.hi.is_finite() || !difference.lo.is_finite() || !error.is_finite() {
            return None;
        }
        maximum = maximum.max(error);
    }
    Some(maximum)
}

#[cfg(test)]
mod tests {
    use super::integral_error;

    #[test]
    fn g1_signed_and_both_gross_directions_resolve_reversal() {
        assert_eq!(integral_error(10.0, [4.0, 6.0]), Some(0.0));
        assert_eq!(integral_error(-10.0, [-4.0, -6.0]), Some(0.0));
        assert_eq!(integral_error(0.0, [3.0, -3.0]), Some(3.0));
        assert_eq!(integral_error(2.0, [5.0, -3.0]), Some(3.0));
        assert_eq!(integral_error(-2.0, [-5.0, 3.0]), Some(3.0));
    }

    #[test]
    fn g1_transport_error_retains_low_parts_and_trace_scales() {
        assert_eq!(1e16 + 1.0, 1e16);
        assert_eq!(integral_error(1e16, [1e16, 1.0]), Some(1.0));
        assert_eq!(integral_error(-1e16, [-1e16, -1.0]), Some(1.0));
        for scale in [1e-250, 1.0, 1e250] {
            let value = integral_error(2.0 * scale, [5.0 * scale, -3.0 * scale]).unwrap();
            assert!((value / scale - 3.0).abs() < 2e-15);
        }
    }

    #[test]
    fn g0_unrepresentable_integrals_cannot_appear_accurate() {
        assert_eq!(integral_error(f64::NAN, [0.0, 0.0]), None);
        assert_eq!(integral_error(0.0, [f64::INFINITY, 0.0]), None);
        assert_eq!(integral_error(-f64::MAX, [f64::MAX, 0.0]), None);
    }
}
