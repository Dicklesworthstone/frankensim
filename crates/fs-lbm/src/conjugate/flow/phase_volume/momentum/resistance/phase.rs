//! Original material sources inside the coupled resistance/heat/volume step.

use super::super::super::refuse;
use super::super::ThermalVolumeMomentumInput;
use super::{
    ResistanceResponse, ResistiveThermalVolumeMomentumControls, ResistiveThermalVolumeMomentumStep,
};
use crate::conjugate::flow::density_projection::BoundaryPressure;
use crate::conjugate::flow::{FlowFace, FlowField, StaggeredBoundaryPressure, poll_every};
use crate::conjugate::{ChtError, VoxelDomain, finite, finite_positive, poll};
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::mushy::{MushyPhaseLaw, MushyPhaseState, PhaseMobility};
use fs_material::phase::transport::ResolvedPhaseTransportState;
use std::borrow::Cow;

/// One accepted source-resolved resistance, phase, volume and momentum endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedResistiveThermalVolumeMomentumStep {
    /// Original M/P inertia, accepted material routes and physical drag heat.
    /// Applied beta is retained alongside its physical endpoint mismatch gates.
    pub evolution: ResistiveThermalVolumeMomentumStep,
    /// Sources resolved at the RETURNED endpoint phases, one per fluid cell.
    /// Their beta agrees with the applied field within the native accepted
    /// impulse/heat budgets; these receipts are not stale predictor sources.
    pub endpoint_sources: Vec<Option<ResolvedPhaseTransportState>>,
}

impl FlowField {
    /// Close phase-dependent resistance and physical drag heat through the
    /// original chart's density/volume response and conservative material flow.
    ///
    /// Each bounded trial starts from the original M, h and P. The native
    /// resistance owner supplies physical cell heat; the thermal-volume owner
    /// transports that heat and admits the actual returned occupied volumes.
    /// `resolve` then evaluates each fluid cell's returned phase through the
    /// original material source, and mu*K^-1 supplies the next trial's beta.
    /// The common owner gates both caloric heat and constitutive impulse/heat
    /// residuals on the same accepted velocities before publishing any state.
    ///
    /// Resolve the supplied state with
    /// `resolve_phase_transport_state_for_volume_change`: true phase VOLUME
    /// fraction controls permeability, and a kinematic-viscosity source uses
    /// its chart-matching current density. Complete phase/law/source identity
    /// and current source/transport density remain mandatory. Native inertia
    /// still uses authoritative original M/V. No density substitution occurs.
    ///
    /// Original, incoming, endpoint and accepted caloric-intermediate states
    /// must remain mobile under this law. Exact solids require a separate
    /// immobile momentum owner and refuse before their liquid source query.
    /// Source queries occur for original and trial endpoint states; accepted
    /// caloric-intermediate admission adds no liquid-property query.
    ///
    /// Mechanical pressure work remains in the native account and numerical
    /// damping never becomes heat. This bounded split retains the original
    /// fixed thermodynamic-pressure/composition convention. Shear coordinates
    /// are caller-declared; this entry does not infer shear or move a skeleton.
    ///
    /// # Errors
    /// Any source/phase inconsistency, unsupported solid state, unrepresentable
    /// beta, native refusal, exhausted coupling budget or cancellation returns
    /// no candidate and leaves all supplied physical inventories unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_resolved_phase_resistance_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedResistiveThermalVolumeMomentumStep, ChtError> {
        self.advance_resolved_phase_resistance_thermal_volume_with_reservoirs(
            domain, faces, None, law, input, controls, gate, resolve,
        )
    }

    /// The same coupled step with the original staggered reservoir pressure
    /// profile. Each trial uses one native pressure field for volume, drag,
    /// momentum and material routes; source thermodynamic pressure is retained.
    ///
    /// # Errors
    /// Invalid profile or any source, phase, balance or coupling refusal from
    /// [`Self::advance_resolved_phase_resistance_thermal_volume_momentum`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_resolved_phase_resistance_thermal_volume_momentum_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedResistiveThermalVolumeMomentumStep, ChtError> {
        self.advance_resolved_phase_resistance_thermal_volume_with_reservoirs(
            domain,
            faces,
            Some(BoundaryPressure::Profile(profile)),
            law,
            input,
            controls,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_resolved_phase_resistance_thermal_volume_with_reservoirs(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: Option<BoundaryPressure<'_>>,
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedResistiveThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        for phase in input.volume.boundary_inflow.into_iter().flatten() {
            admit_mobile_phase(law, &phase)?;
        }
        let (evolution, endpoint_sources) = self.advance_thermal_volume_resistance_with_response(
            domain,
            faces,
            boundary_pressure,
            law.curve(),
            input,
            controls,
            gate,
            |phases| {
                let n = domain.cell_count();
                if phases.len() != n {
                    return Err(refuse(
                        "phase resistance requires one state per original cell",
                    ));
                }
                let mut beta = vec![0.0; n];
                let mut sources = Vec::with_capacity(n);
                for cell in 0..n {
                    poll_every(gate, cell)?;
                    if !domain.is_fluid(cell) {
                        if phases[cell].is_some() {
                            return Err(refuse("fixed grid solids cannot contain phase inventory"));
                        }
                        sources.push(None);
                        continue;
                    }
                    let phase = phases[cell]
                        .ok_or_else(|| refuse("fluid resistance needs an original chart state"))?;
                    let mobile = admit_mobile_phase(law, &phase)?;
                    let source = resolve(cell, &phase)?;
                    poll(gate)?;
                    if *source.mobility() != PhaseMobility::Mobile(mobile)
                        || *source.phase_state() != phase
                        || source.source_density_kg_m3().to_bits()
                            != phase.bulk_density_kg_m3().to_bits()
                    {
                        return Err(refuse(
                            "resistance source must match its exact phase, law and current density",
                        ));
                    }
                    let transport = source.mobile_transport().ok_or_else(|| {
                        refuse("mobile phase resistance needs original fluid transport")
                    })?;
                    if transport.density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits() {
                        return Err(refuse(
                            "phase viscosity must retain the current source density",
                        ));
                    }
                    let viscosity = transport.dynamic_viscosity_pa_s();
                    finite_positive("flow.thermal_resistance_dynamic_viscosity", viscosity)?;
                    let inverse = mobile.inverse_permeability_per_m2;
                    beta[cell] = if inverse == 0.0 {
                        0.0
                    } else {
                        viscosity * inverse
                    };
                    finite("flow.thermal_phase_resistance_beta", beta[cell])?;
                    if inverse != 0.0 && beta[cell] <= 0.0 {
                        return Err(refuse("nonzero phase resistance is not representable"));
                    }
                    sources.push(Some(source));
                }
                Ok(ResistanceResponse {
                    resistance_pa_s_m2: Cow::Owned(beta),
                    sources,
                })
            },
        )?;
        // The accepted caloric split can differ from its transported endpoint.
        // Admit that state too, without spending another source query or
        // treating incoming mobile material as a hidden solid-release owner.
        for cell in 0..domain.cell_count() {
            poll_every(gate, cell)?;
            if !domain.is_fluid(cell) {
                continue;
            }
            let original = input.volume.phase_states[cell]
                .ok_or_else(|| refuse("mobile thermal resistance needs its original phase"))?;
            let external = input.volume.net_enthalpy_increment_j[cell];
            let heat = evolution.caloric_resistance_heat_j[cell];
            let increment = if heat == 0.0 {
                external
            } else {
                external + heat
            };
            let heated = law
                .curve()
                .advance_specific_energy(original, increment / input.volume.mass_kg[cell])
                .map_err(|error| refuse(error.to_string()))?;
            admit_mobile_phase(law, &heated)?;
        }
        poll(gate)?;
        Ok(ResolvedResistiveThermalVolumeMomentumStep {
            evolution,
            endpoint_sources,
        })
    }
}

fn admit_mobile_phase(
    law: &MushyPhaseLaw<'_>,
    phase: &EquilibriumPhaseState,
) -> Result<MushyPhaseState, ChtError> {
    if phase.phase_curve_identity() != law.curve().identity()
        || phase.material_card_identity() != law.curve().material_card_identity()
    {
        return Err(refuse(
            "thermal resistance needs the original phase chart and material",
        ));
    }
    let mobility = law
        .classify_at_specific_enthalpy(phase.specific_enthalpy_j_kg())
        .map_err(|error| refuse(error.to_string()))?;
    if mobility.phase() != *phase {
        return Err(refuse(
            "thermal resistance needs the exact h-resolved phase",
        ));
    }
    match mobility {
        PhaseMobility::Mobile(mobile) => Ok(mobile),
        PhaseMobility::ImmobileSolid { .. } => Err(refuse(
            "exact solid phase needs an explicit immobile momentum owner",
        )),
    }
}

#[cfg(test)]
mod tests;
