//! Preserve finite electrical supply when a spatial birth batch is refused.
//! The existing electrical/beam interval and thermal scheduler remain owners.

use core::convert::Infallible;
use crate::material_circuit::MaterialCircuit;
use crate::circuit_enthalpy::material::{MaterialNetworkSample, PreparedMaterialNetwork};
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkError, HeatNetworkSample,
    MaterialHeatNetwork, PreparedHeatNetwork};
use fs_exec::Cx;
use crate::circuit_enthalpy::contacts::material::PreparedElectricalContactNetwork;
use fs_solid::bridge as solid;
use super::{BirthAudit, GrowingBridge, poll};
use super::super::{BridgeEnthalpyError, invalid, electrical as feedback};

#[derive(Debug)]
pub enum ElectricalGrowthError {
    Interval(feedback::ElectricalError),
    Birth(BridgeEnthalpyError<Infallible>),
    /// The physical interval/births passed, but their declared successor did not.
    Successor(HeatNetworkError),
}
impl core::fmt::Display for ElectricalGrowthError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrical growing bridge: {self:?}")
    }
}
impl std::error::Error for ElectricalGrowthError {}

pub struct ElectricalGrowthReport<S = MaterialNetworkSample> {
    /// Original left-limit electrical, mechanical and thermal interval accounts.
    pub interval: feedback::ElectricalReport<S>,
    /// Actual zero-duration arrivals, including passive-body incoming enthalpy.
    pub births: BirthAudit,
    /// New material has not equilibrated under the enlarged load field yet.
    pub birth_equilibrium_pending: bool,
    /// Electrical + elastic + thermal storage and numerical loss, less electrical
    /// supply, load work, independent boundary heat and incoming enthalpy, J.
    pub total_energy_residual_j: f64,
}

impl<'t,'m> GrowingBridge<'t,'m> {
    /// Solve the existing electrical/beam feedback, then accept the actual birth
    /// batch before returning ANY endpoint. A failure in the last birth discards
    /// the solved interval, circuit proposal and earlier successful appends.
    ///
    /// Inputs describe the pre-event topology; split at events and source/load
    /// discontinuities. No circuit branch is inferred from new material. Charges,
    /// fluxes and electrical topology survive zero-duration births unchanged.
    /// The next interval needs a complete material-network binding for the new
    /// thermal cells; explicit preparation/source pins remain the caller's choice.
    /// This does not migrate circuit topology, create an ideal supply, or bypass
    /// the bulk electrical path's refusal of unrouted series contacts.
    #[allow(clippy::too_many_arguments)]
    pub fn step_electrical_candidate(
        &self, circuit: &MaterialCircuit, transport: &PreparedMaterialNetwork<'_,'_>,
        ports: &[f64], target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
        config: feedback::ElectricalConfig, guess_j: &[f64], mut cancelled: impl FnMut()->bool,
    ) -> Result<(Self,MaterialCircuit,ElectricalGrowthReport),ElectricalGrowthError> {
        let candidate=feedback::step_candidate(&self.thermal,&self.bridge,&self.bindings,
            circuit,transport,ports,target_s,gravity_m_s2,loads_n_m,config,guess_j,&mut cancelled)
            .map_err(ElectricalGrowthError::Interval)?;
        self.finish_electrical_candidate(candidate,config,&mut cancelled)
    }

    /// Accept contact-heated circuit, beam and thermal histories together with
    /// the scheduled birth batch. Both contact-side heat destinations and source
    /// selections remain owned by the prepared electrical contact network.
    ///
    /// This differs from `step_electrical_candidate` only in the interval owner:
    /// it includes the original contact-current and explicit A/B heat-partition
    /// calculation. Births still run after the interval, never inside Newton
    /// trials. A refused final birth also discards the finite electrical debit
    /// and the passive contact receiver's proposed temperature change.
    ///
    /// Inputs and contact geometry describe the pre-event topology. No newborn
    /// circuit branch, interface, heat fraction or post-birth electrical solve
    /// is inferred; prepare a complete successor binding before the next step.
    #[allow(clippy::too_many_arguments)]
    pub fn step_electrical_with_contacts_candidate(
        &self, circuit: &MaterialCircuit,
        transport: &PreparedElectricalContactNetwork<'_,'_>,
        ports: &[f64], target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
        config: feedback::ElectricalConfig, guess_j: &[f64], mut cancelled: impl FnMut()->bool,
    ) -> Result<(Self,MaterialCircuit,ElectricalGrowthReport),ElectricalGrowthError> {
        let candidate=feedback::step_with_contacts_candidate(
            &self.thermal,&self.bridge,&self.bindings,circuit,transport,ports,
            target_s,gravity_m_s2,loads_n_m,config,guess_j,&mut cancelled)
            .map_err(ElectricalGrowthError::Interval)?;
        self.finish_electrical_candidate(candidate,config,&mut cancelled)
    }

    /// Material thermal contacts and closed radiation share the existing
    /// electrical/beam interval and the SAME atomic birth-completion path.
    /// A late birth failure discards radiative changes in passive receivers as
    /// well as the finite circuit debit, beam update and all earlier appends.
    ///
    /// The network describes pre-event geometry. After accepted births, supply
    /// a complete successor network before continuing; no newborn radiator,
    /// conductive contact, circuit branch or view factor is inferred here.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network_candidate(
        &self, cx: &Cx<'_>, circuit: &MaterialCircuit,
        transport: &PreparedHeatNetwork<'_,'_>, ports: &[f64],
        target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
        config: feedback::ElectricalConfig, guess_j: &[f64],
    ) -> Result<(Self,MaterialCircuit,ElectricalGrowthReport<HeatNetworkSample>),ElectricalGrowthError> {
        let candidate=feedback::step_heat_network_candidate(cx,&self.thermal,&self.bridge,&self.bindings,
            circuit,transport,ports,target_s,gravity_m_s2,loads_n_m,config,guess_j)
            .map_err(ElectricalGrowthError::Interval)?;
        self.finish_electrical_candidate(candidate,config,&mut ||cx.checkpoint().is_err())
    }

    /// Accept the physical interval, birth batch AND source-preserving successor
    /// heat network together. A late source or exposed-surface refusal returns
    /// none of the circuit, beam, heat field or successor binding.
    ///
    /// `transport` describes the old topology; `successor` describes the complete
    /// post-event topology. The surface map has the exact semantics of
    /// [`PreparedHeatNetwork::extend_fixed_circuit`]. Surviving source claims are
    /// retained before new coefficients are evaluated. Returned bindings borrow
    /// only the authored successor, not the returned mutable physical state.
    ///
    /// This performs no extra physical step or source debit. Use the returned
    /// circuit and prepared network for the next interval. New surfaces, view
    /// factors and thermal links must be authored; the circuit stays unchanged.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_heat_network_with_successor_candidate<'n, 'a>(
        &self, cx: &Cx<'_>, circuit: &MaterialCircuit,
        transport: &PreparedHeatNetwork<'_,'_>, successor: &'n MaterialHeatNetwork<'a>,
        surface_successors: &[Option<usize>], ports: &[f64],
        target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
        config: feedback::ElectricalConfig, guess_j: &[f64],
    ) -> Result<(Self,MaterialCircuit,PreparedHeatNetwork<'n,'a>,
        ElectricalGrowthReport<HeatNetworkSample>),ElectricalGrowthError> {
        let (next,next_circuit,report)=self.step_heat_network_candidate(
            cx,circuit,transport,ports,target_s,gravity_m_s2,loads_n_m,config,guess_j)?;
        let prepared=transport.extend_fixed_circuit(cx,successor,next.thermal(),surface_successors)
            .map_err(ElectricalGrowthError::Successor)?;
        // The extension reuses original source/geometry admission and its final
        // checkpoint. Neither it nor this composition advances a second clock.
        Ok((next,next_circuit,prepared,report))
    }

    // One birth/energy/publication path for every admitted heat-network payload.
    // The candidate carries the original admitted circuit, not a reconstructed
    // charge/flux state or a second electrical advance during deposition.
    fn finish_electrical_candidate<S>(
        &self, candidate: feedback::ElectricalCandidate<'t,'m,S>,
        config: feedback::ElectricalConfig, cancelled: &mut impl FnMut()->bool,
    ) -> Result<(Self,MaterialCircuit,ElectricalGrowthReport<S>),ElectricalGrowthError> {
        let (thermal,bridge,next_circuit,interval)=candidate.into_parts();
        let staged=Self {thermal,bridge,bindings:self.bindings.clone(),future:self.future.clone(),
            body:self.body,pending:false};
        let c=config.bridge.feedback;
        let mechanics=solid::Config {section:c.mechanics,
            maximum_deflection_m:config.bridge.maximum_deflection_m};
        let (next,births)=staged.synchronize(mechanics,c.thermal_energy_tolerance_j,cancelled)
            .map_err(ElectricalGrowthError::Birth)?;
        // Reconstruct from the original accepted storage and the fully born
        // endpoint, not by summing two reports that might share an internal heat.
        let total=(interval.electrical.stored_after_j-interval.electrical.stored_before_j)
            +(next.bridge.stored_energy_j()-self.bridge.stored_energy_j())
            +(next.thermal.energy_j()-self.thermal.energy_j())
            +interval.mechanical.energy_j.numerical_dissipation
            -interval.electrical.native.supplied-interval.mechanical.load_work_j
            -interval.other_external_heat_j-births.thermal.incoming_energy_j
            -births.thermal.external_heat_j;
        if !total.is_finite() || total.abs()>c.total_energy_tolerance_j {
            return Err(ElectricalGrowthError::Birth(invalid(
                "electrical storage, spatial growth and incoming enthalpy do not close")));
        }
        poll(cancelled).map_err(|error|ElectricalGrowthError::Birth(error.into()))?;
        let pending=next.pending;
        Ok((next,next_circuit,ElectricalGrowthReport {interval,births,
            birth_equilibrium_pending:pending,total_energy_residual_j:total}))
    }
}
