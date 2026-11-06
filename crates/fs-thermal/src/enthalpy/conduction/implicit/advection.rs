//! A prepared directed tangent for the original endpoint upwind flux owner.
use super::*;
use crate::enthalpy::advection::Advected;
use crate::enthalpy::advection::entropy::EndpointAdvectionEntropy;
use crate::enthalpy::advection::reconstruction::LimitedEnthalpyReconstruction;
use crate::enthalpy::advection::schedule::AdvectionSchedule;
use std::borrow::Cow;
use crate::enthalpy::entropy::AdvectiveEntropyContribution;

#[derive(Clone, Copy)]
pub(super) struct ClassifiedAdvectionEntropy {
    pub(super) advection: AdvectiveEntropyContribution,
    pub(super) source_entropy_j_k: f64,
    pub(super) has_unclassified_sources: bool,
}

impl ClassifiedAdvectionEntropy {
    pub(super) fn sample(
        thermal: &DepositionSimulation<'_>,
        sources_w: &BTreeMap<CellId, f64>,
        transfers: &[HeatTransfer],
        endpoint: EndpointAdvectionEntropy,
        duration_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        let mut source = Sum::default();
        let mut has_unclassified_sources = false;
        let temperature = |id| -> Result<f64, DepositionError> {
            let t = *endpoint.temperatures.get(&id).ok_or(DepositionError::InvalidInput(
                "inactive entropy source cell",
            ))?;
            if positive(t) { Ok(t) } else { Err(DepositionError::MaterialDomain) }
        };
        for (i, (&id, cell)) in thermal.cells.iter().enumerate() {
            if i % 256 == 0 { checkpoint()?; }
            let authored = cell.boundary.source_w;
            let generated = sources_w.get(&id).copied().unwrap_or(0.0);
            has_unclassified_sources |= authored != 0.0 || generated != 0.0;
            source.add(finite(finite(authored + generated)? / temperature(id)?)?)?;
        }
        for (i, transfer) in transfers.iter().enumerate() {
            if i % 256 == 0 { checkpoint()?; }
            has_unclassified_sources |= transfer.power_w != 0.0;
            source.add(finite(-transfer.power_w / temperature(transfer.a)?)?)?;
            source.add(finite(transfer.power_w / temperature(transfer.b)?)?)?;
        }
        Ok(Self {
            advection: endpoint.rates.integrated(duration_s)?,
            source_entropy_j_k: finite(duration_s * source.total()?)?,
            has_unclassified_sources,
        })
    }
}

pub(super) struct ImplicitAdvection {
    pub(super) flow: Advected<()>,
    // d(heat into row)/d(cell energy in column), in 1/s. Fixed mass and
    // prescribed flow make this constant even for nonlinear/latent T(h).
    pub(super) tangent: Vec<(usize, usize, f64)>,
    // The same compiled connectivity in carried-enthalpy coordinates, divided
    // by original cell masses. Liquid mobility changes its column slopes at
    // each endpoint, so these entries must never enter the static Jacobian.
    liquid_tangent: Vec<(usize, CellId, f64)>,
    reconstruction: Option<LimitedEnthalpyReconstruction>,
    entropy: bool,
    schedule: Option<AdvectionSchedule>,
}

pub(super) struct ImplicitAdvectionStage<'a> {
    flow: Cow<'a, Advected<()>>,
    pub(super) multiplier: f64,
}

impl ImplicitAdvection {
    pub(super) fn stage(
        &self,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitAdvectionStage<'_>, DepositionError> {
        if let Some(schedule) = &self.schedule {
            let (flow, multiplier) = schedule.flow_at(&self.flow, time_s, checkpoint)?;
            Ok(ImplicitAdvectionStage { flow: Cow::Owned(flow), multiplier })
        } else {
            Ok(ImplicitAdvectionStage { flow: Cow::Borrowed(&self.flow), multiplier: 1.0 })
        }
    }

    pub(super) fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.schedule.as_ref().and_then(|schedule| schedule.next_transition_time_s(after_s))
    }

    pub(super) fn endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        stage: &ImplicitAdvectionStage<'_>,
        index: &BTreeMap<CellId, usize>,
        with_entropy: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (TransportStepCoefficients<()>, Vec<(usize, usize, f64)>, Option<EndpointAdvectionEntropy>),
        ImplicitConductionError,
    >
    {
        let flow = stage.flow.as_ref();
        // Liquid upwind columns remain dynamic with or without reconstruction.
        // Face corrections below augment this same carrier tangent; they must
        // neither discard it nor reintroduce the bulk-storage derivative.
        let mut tangent = Vec::new();
        tangent
            .try_reserve_exact(self.liquid_tangent.len())
            .map_err(|_| TridiagonalError::Capacity)?;
        if stage.multiplier != 0.0 {
            for &(row, cell, coefficient) in &self.liquid_tangent {
                checkpoint()?;
                let slope = flow.implicit_carried_enthalpy_slope(&thermal.cells[&cell])?;
                if slope != 0.0 {
                    let derivative = scaled_tangent(finite(coefficient * slope)?, stage.multiplier)?;
                    tangent.push((row, index[&cell], derivative));
                }
            }
        }
        let reconstruction = self.reconstruction.as_ref().filter(|_| stage.multiplier != 0.0);
        let Some(reconstruction) = reconstruction else {
            if self.entropy && with_entropy {
                let (transport, entropy) =
                    flow.endpoint_transport_entropy_state(thermal, checkpoint)?;
                return Ok((transport, tangent, Some(entropy)));
            }
            return Ok((
                flow.endpoint_transport(thermal, checkpoint)?,
                tangent,
                None,
            ));
        };
        let (transport, rates) =
            flow.reconstructed_endpoint(thermal, reconstruction, true, &mut checkpoint)?;
        tangent
            .try_reserve_exact(rates.len())
            .map_err(|_| TridiagonalError::Capacity)?;
        for (row, col, rate) in rates {
            checkpoint()?;
            let derivative = finite(rate / thermal.cells[&col].mass_kg)?;
            if derivative == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            tangent.push((index[&row], index[&col], derivative));
        }
        Ok((transport, tangent, None))
    }
}

impl ImplicitConductionChain<'_> {
    pub(super) fn next_advection_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.advection.as_ref().and_then(|flow| flow.next_transition_time_s(after_s))
    }

    /// Solve prescribed upwind enthalpy and nonlinear conduction together.
    ///
    /// This binds an immutable snapshot of the original `Advected` graph. Its
    /// endpoint query retains native continuity, material, caloric-datum and
    /// inlet/outlet admission. Every Newton residual transports the proposed
    /// endpoint h, including latent energy at constant temperature. The same
    /// flow must not also be supplied as an explicit source observer.
    ///
    /// The exact directed tangent is compiled once in joule coordinates. A
    /// direct chain needs a tridiagonal combined operator; other flow graphs
    /// require `ImplicitConductionGraph`. A reference-coordinate derivative is
    /// retained even for admitted continuity roundoff. No symmetric diffusion
    /// surrogate or numerical enthalpy differences are used.
    ///
    /// Ordinary fixed-mass bulk enthalpy and the original stationary-skeleton
    /// liquid carrier are admitted. For liquid flow, the carrier tangent is zero
    /// on its latent plateau and one on its fully liquid branch, reevaluated at
    /// every trial. Thus liquid circulation does not advect stationary-solid
    /// latent energy. Fully solid flow still refuses; no mobility is inferred.
    ///
    /// Stationary species may use distinct immutable mixture caloric laws.
    /// Their original inlet-fed species balances, exact cell compositions and
    /// common enthalpy reference remain authoritative at every endpoint. The
    /// compiled tangent transports stored mixture h, so unlike heat capacities
    /// are never averaged into a temperature flux. Species do not evolve here.
    /// Selected outlets retain their prescribed h and composition. Their exact
    /// tangent has zero outgoing-state slope while retaining the component
    /// reference derivative; they are not ordinary local-enthalpy outlets.
    /// Source law bindings are optional for energy and required for entropy.
    /// An attached explicit flow entropy history still refuses.
    /// The native thermal entropy account otherwise
    /// leaves mass-carried transfers unclassified unless
    /// [`Self::with_implicit_advection_entropy`] is selected. The supplied
    /// observer's history is unchanged, including on success.
    ///
    /// No explicit flow CFL is imposed on this implicit term. The selected
    /// backward-Euler or SDIRK2 time error and upwind spatial diffusion still
    /// require refinement; the existing adaptive h/T/phase controllers can
    /// refine these same endpoints. SDIRK2 re-evaluates transport at both stages.
    /// Other observed sources retain all their original timestep restrictions.
    pub fn with_implicit_advection(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        self.bind_advection(thermal, flow, None, checkpoint)
    }

    /// Bind a prescribed pump and inlet history to the same native flow.
    ///
    /// One immutable flow snapshot is sampled at each actual BE or SDIRK stage,
    /// then reused by every Newton/line-search trial. A common nonnegative flow
    /// multiplier preserves topology, cell masses and frozen stationary species.
    /// Inlet values retain their original caloric law and common enthalpy datum.
    /// The native time owners clip to schedule knots, including adaptive trials.
    ///
    /// Zero-flow stages carry no mass and require no liquid mobility. Later
    /// nonzero stages still require native per-cell continuity, original liquid
    /// support and source-domain admission. This does not infer a pump momentum
    /// law or silently evolve composition/density. The base flow is not also
    /// supplied as an explicit source observer.
    pub fn with_scheduled_implicit_advection(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        schedule: &AdvectionSchedule,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        self.bind_advection(thermal, flow, Some(schedule.clone()), checkpoint)
    }

    /// Add a history to an already attached upwind or reconstructed flow.
    /// Existing entropy/reconstruction selections retain their original rules.
    pub fn with_implicit_advection_schedule(
        mut self,
        thermal: &DepositionSimulation<'_>,
        schedule: &AdvectionSchedule,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        self.admit(thermal, &mut checkpoint)?;
        let advection = self.advection.as_mut().ok_or(ImplicitConductionError::Invalid(
            "advection schedule requires an attached implicit flow",
        ))?;
        if advection.schedule.is_some() {
            return Err(ImplicitConductionError::Invalid("implicit advection already has a schedule"));
        }
        schedule.admit(&advection.flow, thermal, &mut checkpoint)?;
        advection.schedule = Some(schedule.clone());
        let stage = advection.stage(thermal.time_s(), &mut checkpoint)?;
        advection.endpoint(thermal, &stage, &self.index, advection.entropy, &mut checkpoint)?;
        drop(stage);
        checkpoint()?;
        Ok(self)
    }

    fn bind_advection(
        mut self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        schedule: Option<AdvectionSchedule>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if self.advection.is_some() {
            return Err(ImplicitConductionError::Invalid(
                "implicit advection already bound",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let rates = flow.implicit_enthalpy_tangent(&mut checkpoint)?;
        // Admit the actual initial prescribed flow with its original owner.
        // A stopped liquid pump does not require mobile liquid at preparation.
        if let Some(schedule) = &schedule {
            schedule.admit(flow, thermal, &mut checkpoint)?;
            let (initial, _) = schedule.flow_at(flow, thermal.time_s(), &mut checkpoint)?;
            initial.endpoint_transport(thermal, &mut checkpoint)?;
        } else {
            flow.endpoint_transport(thermal, &mut checkpoint)?;
        }
        let mut tangent = Vec::with_capacity(rates.len());
        for (row, col, rate) in rates {
            checkpoint()?;
            let row = *self.index.get(&row).ok_or(ImplicitConductionError::Invalid(
                "inactive implicit advection row",
            ))?;
            let col = *self.index.get(&col).ok_or(ImplicitConductionError::Invalid(
                "inactive implicit advection column",
            ))?;
            if matches!(self.topology, LinearTopology::Chain) && row.abs_diff(col) > 1 {
                return Err(ImplicitConductionError::Invalid(
                    "nontridiagonal enthalpy transport requires the graph solver",
                ));
            }
            let derivative = finite(rate / thermal.cells[&self.ids[col]].mass_kg)?;
            if derivative == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            tangent.push((row, col, derivative));
        }
        let liquid_tangent = if flow.liquid_carrier().is_some() {
            std::mem::take(&mut tangent)
                .into_iter()
                .map(|(row, col, coefficient)| (row, self.ids[col], coefficient))
                .collect()
        } else {
            Vec::new()
        };
        self.advection = Some(ImplicitAdvection {
            flow: flow.clone(),
            tangent,
            liquid_tangent,
            reconstruction: None,
            entropy: false,
            schedule,
        });
        checkpoint()?;
        Ok(self)
    }

    /// Classify entropy of an already attached bulk, liquid or stationary-mixture
    /// upwind flow with its original source entropy closure.
    ///
    /// The existing caloric integral supplies actual storage, and the original
    /// flow owner supplies inlet/outlet entropy and upwind mixing. The native
    /// audit retains them in `entropy.advection`, excluding them from genuine
    /// unclassified heat/work sources. Backward Euler and SDIRK2 use exactly
    /// their native enthalpy-flux quadrature; adaptive comparison/rejected
    /// trajectories publish no entropy. Numerical time/algebraic defects remain
    /// separate from conductive production and spatial upwind mixing.
    ///
    /// Every flowing material must supply its original caloric entropy.
    /// Stationary mixtures additionally require source entropy anchors in one
    /// shared species-reduced convention. Selective outlets additionally bind
    /// their original outgoing law, composition and caloric reference. Their
    /// signed removal term includes the source h/s once, with no extra latent
    /// heat, phase equilibrium or separator work. Mixing/pressure contributions
    /// retain their sign; no pressure-work or full fluid dissipation law is
    /// inferred. A stationary-skeleton liquid carrier uses its original liquid
    /// h/s at every stage while storage uses mixture entropy; no additional
    /// latent-energy source is introduced. Original phase support, source
    /// identity and fixed-density admission remain authoritative.
    /// Nonempty face reconstruction requires its own entropy closure and
    /// refuses this opt-in; its existing energy-only API is unchanged.
    pub fn with_implicit_advection_entropy(
        mut self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        self.admit(thermal, &mut checkpoint)?;
        let advection = self.advection.as_mut().ok_or(ImplicitConductionError::Invalid(
            "implicit entropy requires an attached upwind flow",
        ))?;
        if advection.reconstruction.is_some() {
            return Err(ImplicitConductionError::Invalid(
                "reconstructed enthalpy transport requires its own entropy closure",
            ));
        }
        let stage = advection.stage(thermal.time_s(), &mut checkpoint)?;
        stage.flow.endpoint_transport_with_entropy(thermal, &mut checkpoint)?;
        drop(stage);
        advection.entropy = true;
        checkpoint()?;
        Ok(self)
    }

    pub(super) fn with_implicit_advection_reconstruction(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        reconstruction: &LimitedEnthalpyReconstruction,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !matches!(self.topology, LinearTopology::Graph(_)) {
            return Err(ImplicitConductionError::Invalid(
                "limited enthalpy reconstruction requires the graph solver",
            ));
        }
        let mut result = self.with_implicit_advection(thermal, flow, &mut checkpoint)?;
        if !reconstruction.is_empty() {
            // Existing source, geometry, topology and face-domain admission is
            // authoritative. No reconstruction is frozen at the initial state.
            flow.endpoint_transport_with_reconstruction(thermal, reconstruction, &mut checkpoint)?;
            result
                .advection
                .as_mut()
                .ok_or(DepositionError::Arithmetic)?
                .reconstruction = Some(reconstruction.clone());
        }
        checkpoint()?;
        Ok(result)
    }
}

pub(super) fn scaled_tangent(value: f64, multiplier: f64) -> Result<f64, DepositionError> {
    let scaled = finite(value * multiplier)?;
    if value != 0.0 && multiplier != 0.0 && scaled == 0.0 {
        return Err(DepositionError::Arithmetic);
    }
    Ok(scaled)
}
