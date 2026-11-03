//! Original isobaric, fixed-composition gas sources on finite graph storage.
//!
//! The exact registered GasEnthalpyLaw owns NASA h(T), T(h), caloric support and
//! formation-energy datum. Its original pressure, mole composition and source
//! selections supply each local transport state. This is a low-Mach,
//! fixed-reference-inventory approximation, not compressible pressure evolution.

use super::conduction::{GraphConductionError, GraphConductiveLink, ResolvedGraphConduction};
use super::convection::{GraphConvectionError, GraphConvectiveLink, GraphSourceHeat};
use super::convection::reservoir::{GraphConvectiveReservoir, GraphReservoirConvectionError, GraphReservoirHeat};
use super::*;
use fs_material::thermochemical_gas::{
    ResolvedThermochemicalGasState, energy::enthalpy::GasEnthalpyLaw,
};

pub struct GasGraphSources<'a> {
    source: &'a GasEnthalpyLaw,
    reference: ResolvedThermochemicalGasState,
    maximum_mach: f64,
    core: BoundGraphCore<'a>,
}

fn resolve_gas(
    source: &GasEnthalpyLaw,
    temperature: Temperature,
) -> Result<ResolvedThermochemicalGasState, SourceGraphError> {
    let point = source
        .context()
        .clone()
        .with("temperature", temperature.value())
        .map_err(|_| SourceGraphError::Invalid("invalid original gas temperature coordinate"))?;
    source
        .model()
        .evaluate(&point)
        .map_err(SourceGraphError::Material)
}

impl<'a> GasGraphSources<'a> {
    /// Reference rho comes from the unchanged original source context and must
    /// equal the distributed inventory's conserved reference density exactly.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        source: &'a GasEnthalpyLaw,
        thermal: &DepositionSimulation<'_>,
        binding: SourceGraphBinding<'a>,
        limits: SourceGraphLimits,
        maximum_mach: f64,
    ) -> Result<Self, SourceGraphError> {
        poll(cx)?;
        if !positive(maximum_mach) || maximum_mach > 0.3 {
            return Err(SourceGraphError::Invalid(
                "source gas graph requires a Mach limit in (0,0.3]",
            ));
        }
        let reference = source
            .model()
            .evaluate(source.context())
            .map_err(SourceGraphError::Material)?;
        source
            .enthalpy_at_temperature(reference.state().temperature)
            .map_err(SourceGraphError::GasCaloric)?;
        let transport = reference
            .transport_state()
            .map_err(SourceGraphError::Material)?;
        let core = BoundGraphCore::try_new(cx, source, transport, thermal, binding, limits)?;
        Ok(Self {
            source,
            reference,
            maximum_mach,
            core,
        })
    }

    pub const fn source(&self) -> &'a GasEnthalpyLaw {
        self.source
    }
    pub const fn reference(&self) -> &ResolvedThermochemicalGasState {
        &self.reference
    }
    pub fn inventory(&self) -> &'a DistributedGraphInventory {
        self.core.inventory
    }
    pub const fn maximum_mach(&self) -> f64 {
        self.maximum_mach
    }

    /// Axial conduction uses both actual bulk source states at the native
    /// thermal candidate, retaining original gas pressure and composition.
    /// Existing reciprocal links own heat independently of advective direction;
    /// no stagnant-flow conductivity floor or extra caloric law is introduced.
    pub fn evaluate_conduction_links(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        bindings: &[GraphConductiveLink],
        maximum_links: usize,
    ) -> Result<ResolvedGraphConduction<ResolvedThermochemicalGasState>, GraphConductionError> {
        conduction::evaluate_links(&self.core, cx, thermal, sources, bindings, maximum_links)
    }

    fn admit_mach(
        &self,
        cx: &Cx<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        flows: &[f64],
    ) -> Result<(), SourceGraphError> {
        if flows.len() != self.core.section_areas.len() {
            return Err(SourceGraphError::Invalid(
                "complete original gas flow vector required",
            ));
        }
        for (branch, sections) in sources.sections.iter().enumerate() {
            for (section, source) in sections.iter().enumerate() {
                poll(cx)?;
                let sound = source.momentum.medium.state().sound_speed;
                let speed = flows[branch].abs() / self.core.section_areas[branch][section];
                let mach = speed / sound;
                if !positive(sound) || !nonnegative(mach) || mach > self.maximum_mach {
                    return Err(SourceGraphError::GasMach {
                        branch,
                        section,
                        value: mach,
                        maximum: self.maximum_mach,
                    });
                }
            }
        }
        Ok(())
    }

    /// Resolve original source conditions at actual native caloric temperatures.
    /// The trial viscosity, Boussinesq departure, Reynolds and Mach limits all
    /// apply. Source failures are terminal proposals, never absent flow roots.
    pub fn resolve(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        trial_volume_flows_m3_s: &[f64],
    ) -> Result<ResolvedSourceGraph<ResolvedThermochemicalGasState>, SourceGraphError> {
        let resolved = self
            .core
            .resolve(cx, thermal, trial_volume_flows_m3_s, |t| {
                resolve_gas(self.source, t)
            })?;
        self.admit_mach(cx, &resolved.sources, trial_volume_flows_m3_s)?;
        Ok(resolved)
    }

    /// Re-admit original section laws at actual average, old and new flow,
    /// including zero-crossing support. The largest absolute old/new/average
    /// flow independently bounds Mach and omitted mechanical work. Neither
    /// lower Reynolds limits nor fast endpoints can hide behind an average.
    pub fn admit_actual(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        report: &InertialGraphReport,
    ) -> Result<GraphSourceAdmission, SourceGraphError> {
        self.admit_flow(cx, thermal, sources, GraphFlow::Inertial(report))
    }

    /// Admit a native quasi-steady endpoint with its actual thermal duration.
    /// The rate is held for this interval; no momentum history is inferred.
    /// Original section Reynolds, gas Mach, source identity and non-cancelling
    /// local work caps all apply before its heat coefficients can be accepted.
    pub fn admit_quasi_steady(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        hydraulics: &GraphSolution,
        duration_s: f64,
    ) -> Result<GraphSourceAdmission, SourceGraphError> {
        self.admit_flow(
            cx,
            thermal,
            sources,
            GraphFlow::QuasiSteady {
                hydraulics,
                duration_s,
            },
        )
    }

    fn admit_flow(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        flow: GraphFlow<'_>,
    ) -> Result<GraphSourceAdmission, SourceGraphError> {
        let admission = self.core.admit_flow(cx, thermal, sources, flow)?;
        self.admit_mach(cx, sources, &admission.speed_envelope_m3_s)?;
        Ok(admission)
    }

    /// Reciprocal wall heat uses original gas k, mu and Cp at the declared
    /// bulk/film temperature, retaining this law's pressure and composition.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_links(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        report: &InertialGraphReport,
        bindings: &[GraphConvectiveLink],
        maximum_links: usize,
    ) -> Result<GraphSourceHeat<ResolvedThermochemicalGasState>, GraphConvectionError> {
        self.evaluate_link_flow(
            cx,
            thermal,
            sources,
            GraphFlow::Inertial(report),
            bindings,
            maximum_links,
        )
    }

    /// Source-retaining reciprocal wall heat at a solved quasi-steady endpoint.
    /// The original gas pressure, composition, caloric object and source claims
    /// remain fixed while local temperatures and flow change.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_quasi_steady_links(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        hydraulics: &GraphSolution,
        duration_s: f64,
        bindings: &[GraphConvectiveLink],
        maximum_links: usize,
    ) -> Result<GraphSourceHeat<ResolvedThermochemicalGasState>, GraphConvectionError> {
        self.evaluate_link_flow(
            cx,
            thermal,
            sources,
            GraphFlow::QuasiSteady {
                hydraulics,
                duration_s,
            },
            bindings,
            maximum_links,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_link_flow(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        flow: GraphFlow<'_>,
        bindings: &[GraphConvectiveLink],
        maximum_links: usize,
    ) -> Result<GraphSourceHeat<ResolvedThermochemicalGasState>, GraphConvectionError> {
        let admission = self
            .admit_flow(cx, thermal, sources, flow)
            .map_err(GraphConvectionError::Source)?;
        let convection = convection::evaluate_links(
            &self.core,
            cx,
            thermal,
            flow,
            bindings,
            maximum_links,
            |t| resolve_gas(self.source, t),
        )?;
        Ok(GraphSourceHeat {
            admission,
            convection,
        })
    }
    /// Convective heat from already-declared prescribed wall reservoirs at the
    /// actual native inertial endpoint. The original wall temperatures, areas,
    /// fluid source, caloric identity and physical flow limits remain authoritative.
    /// Apply the returned boundary conductances once inside the thermal residual;
    /// the diagnostic heat powers are not additional cell sources.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_reservoirs(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        report: &InertialGraphReport,
        bindings: &[GraphConvectiveReservoir],
        maximum_reservoirs: usize,
    ) -> Result<GraphReservoirHeat<ResolvedThermochemicalGasState>, GraphReservoirConvectionError> {
        self.evaluate_reservoir_flow(
            cx, thermal, sources, GraphFlow::Inertial(report), bindings, maximum_reservoirs,
        )
    }

    /// Original prescribed-wall convection with the solved quasi-steady flow.
    /// The duration is the actual native thermal interval for source work limits.
    /// Reservoir temperature is held by its original boundary event, while
    /// film/bulk properties and conductance are recomputed at each thermal trial.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_quasi_steady_reservoirs(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        hydraulics: &GraphSolution,
        duration_s: f64,
        bindings: &[GraphConvectiveReservoir],
        maximum_reservoirs: usize,
    ) -> Result<GraphReservoirHeat<ResolvedThermochemicalGasState>, GraphReservoirConvectionError> {
        self.evaluate_reservoir_flow(
            cx,
            thermal,
            sources,
            GraphFlow::QuasiSteady { hydraulics, duration_s },
            bindings,
            maximum_reservoirs,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_reservoir_flow(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedThermochemicalGasState>,
        flow: GraphFlow<'_>,
        bindings: &[GraphConvectiveReservoir],
        maximum_reservoirs: usize,
    ) -> Result<GraphReservoirHeat<ResolvedThermochemicalGasState>, GraphReservoirConvectionError> {
        let admission = self.admit_flow(cx, thermal, sources, flow)
            .map_err(GraphReservoirConvectionError::Source)?;
        let convection = convection::reservoir::evaluate(
            &self.core,
            cx,
            thermal,
            flow,
            bindings,
            maximum_reservoirs,
            |t| resolve_gas(self.source, t),
        )?;
        Ok(GraphReservoirHeat { admission, convection })
    }

}
