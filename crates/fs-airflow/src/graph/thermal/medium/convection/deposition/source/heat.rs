//! Source-resolved fluid cooling on the complete material/circuit heat solve.
//!
//! Convection, bulk and contact Joule heat, conduction and enclosure radiation
//! all see one proposed enthalpy endpoint. The original thermal and circuit
//! drivers remain the only state and time-integration owners.

use super::super::ConvectionSourceError;
use super::links::{ResolvedLinkConvection, SourceConvectiveLink, evaluate_source_link_convection};
use super::{
    ResolvedCellConvectionSources, SourceCellConvectionError, SourceConvectiveCell,
    evaluate_source_cell_convection,
};
use fs_conduction::transient::deposition::DepositionSimulation;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyCandidate, CircuitEnthalpyError, CircuitEnthalpyInput, CircuitEnthalpySolve,
    CircuitThermalCoefficients,
    contacts::{
        self, ContactThermalCoefficients, ContactThermalSample,
        heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork},
    },
};
use fs_couple::material_circuit::MaterialCircuit;
use fs_exec::Cx;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Temperature;

/// Original source failures and the complete closed-system energy gate.
#[derive(Debug)]
pub enum ConvectionHeatError<E> {
    Heat(HeatNetworkError),
    Convection(SourceCellConvectionError<E>),
    /// Closed enclosure radiation is internal exchange. Its net numerical
    /// defect may not masquerade as external supply in the combined balance.
    ClosedEnergy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl<E> From<HeatNetworkError> for ConvectionHeatError<E> {
    fn from(error: HeatNetworkError) -> Self {
        match error {
            HeatNetworkError::Cancelled => Self::Cancelled,
            other => Self::Heat(other),
        }
    }
}
impl<E> From<SourceCellConvectionError<E>> for ConvectionHeatError<E> {
    fn from(error: SourceCellConvectionError<E>) -> Self {
        match error {
            SourceCellConvectionError::Convection(ConvectionSourceError::Cancelled) => {
                Self::Cancelled
            }
            other => Self::Convection(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ConvectionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material convection/heat network: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ConvectionHeatError<E> {}

/// Every heat-transfer source at the same accepted endpoint.
#[derive(Debug)]
pub struct ConvectionHeatSample<F, L = F> {
    pub heat: HeatNetworkSample,
    /// Exchange with prescribed reservoirs, already in the source-power map.
    pub convection: ResolvedCellConvectionSources<F>,
    /// Reciprocal exchange with actual finite fluid cells, already in the
    /// native link-conductance map. This is never external heat supply.
    pub finite_convection: Option<ResolvedLinkConvection<L>>,
}

/// A prepared heat assembly and physical fluid boundaries on its actual cells.
///
/// Keep the original fluid models in `resolve_fluid` across nonlinear trials
/// and accepted steps. For example, capture a prepared liquid-temperature law
/// or a source-pinned gas-mixture model. Temperature is supplied by this adapter;
/// pressure, composition, geometry and relative flow remain explicitly authored
/// reservoir conditions. No fluid type or material name selects a correlation.
///
/// The callback sees the same proposed thermal state as the complete heat
/// network. It must be repeatable and preserve material identity, selected
/// claims and source support. Its retained return value can be any original
/// `ThermalFluidTransport` owner; no duplicate fluid evaluator is introduced.
pub struct SourceConvectionHeatNetwork<'p, 'n, 'a, 'b, R> {
    pub heat: &'p PreparedHeatNetwork<'n, 'a>,
    pub cells: &'b [SourceConvectiveCell],
    pub maximum_cells: usize,
    pub resolve_fluid: R,
}

impl<R> SourceConvectionHeatNetwork<'_, '_, '_, '_, R> {
    /// Assemble all mechanisms on one proposed enthalpy state.
    ///
    /// This can also serve an existing nonlinear thermal/coupling residual.
    /// Original exposed-area coverage is authoritative. A convecting cell must
    /// have zero authored ambient h so this source cannot duplicate a Robin
    /// boundary. Contact/radiation/source failures return no partial assembly.
    ///
    /// Convection is external heat from prescribed reservoirs. Closed enclosure
    /// radiation remains internal exchange; consumers must exclude its signed
    /// roundoff net when closing the full energy account, as `step_candidate`
    /// does. Returned powers already include convection exactly once.
    pub fn resolve<F: ThermalFluidTransport, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<ConvectionHeatSample<F>>, ConvectionHeatError<E>>
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        self.resolve_combined::<F, F, E>(cx, thermal, None)
    }

    /// Resolve finite fluid links AND prescribed fluid reservoirs on the same
    /// complete heat-network trial. A wall can exchange with either or both on
    /// distinct native interface patches; original contact coverage remains
    /// authoritative for its exposed reservoir area.
    ///
    /// Prepare `heat` with `prepare_with_link_conductances` using the original
    /// links' initial resolved G. Its Fourier bindings must exclude these links.
    /// Keep the source-pinned fluid laws in `resolve_links`; every trial gets
    /// the actual wall/fluid enthalpy temperatures and the declared bulk/film
    /// property temperature. The two finite cells keep their original caloric
    /// laws, inventories and energy ownership. No bath or duplicate fluid state
    /// is introduced. Empty `cells` permits a completely closed fluid network.
    pub fn resolve_with_links<F: ThermalFluidTransport, L: ThermalFluidTransport, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        links: &[SourceConvectiveLink],
        maximum_links: usize,
        resolve_links: impl Fn(
            &SourceConvectiveLink, &DepositionSimulation<'_>, Temperature,
        ) -> Result<L, E>,
    ) -> Result<ContactThermalCoefficients<ConvectionHeatSample<F, L>>, ConvectionHeatError<E>>
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        cx.checkpoint().map_err(|_| ConvectionHeatError::Cancelled)?;
        let finite = evaluate_source_link_convection(
            thermal, links, maximum_links, resolve_links, || cx.checkpoint().is_err(),
        )?;
        self.resolve_combined(cx, thermal, Some(finite))
    }

    fn resolve_combined<F: ThermalFluidTransport, L, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        finite_convection: Option<ResolvedLinkConvection<L>>,
    ) -> Result<ContactThermalCoefficients<ConvectionHeatSample<F, L>>, ConvectionHeatError<E>>
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        cx.checkpoint()
            .map_err(|_| ConvectionHeatError::Cancelled)?;
        let mut coefficients = match &finite_convection {
            Some(finite) => self.heat.resolve_with_link_conductances(
                cx, thermal, finite.conductances_w_k(),
            ),
            None => self.heat.resolve(cx, thermal),
        }?;
        let convection = evaluate_source_cell_convection(
            thermal,
            self.cells,
            self.maximum_cells,
            |binding, temperature| (self.resolve_fluid)(binding, thermal, temperature),
            || cx.checkpoint().is_err(),
        )?;
        convection
            .convection()
            .add_to(&mut coefficients.bulk.sources_w, || {
                cx.checkpoint().is_err()
            })
            .map_err(SourceCellConvectionError::from)?;
        cx.checkpoint()
            .map_err(|_| ConvectionHeatError::Cancelled)?;
        let bulk = coefficients.bulk;
        Ok(ContactThermalCoefficients {
            contacts: coefficients.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: bulk.resistors,
                sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: ConvectionHeatSample {
                    heat: bulk.source,
                    convection,
                    finite_convection,
                },
            },
        })
    }

    /// Evolve the original circuit and finite thermal cells with material cooling.
    ///
    /// Every residual re-evaluates the body's source-bound electrical resistance,
    /// bulk/contact/layer conduction, emissivity and the fluid's bulk/film
    /// properties at one candidate endpoint. Thus a fluid data change can change
    /// cooling, body temperature, electrical current and stored circuit energy.
    /// Original contact heat partitions and source pins are retained.
    ///
    /// The original report counts convection in `non_electrical_heat_j`.
    /// True external heat is that field minus dt times the retained enclosure's
    /// `radiation_net_power_w()`; a second gate below checks the corrected total
    /// residual against the SAME total-energy allowance. No internal heat is
    /// deposited twice and no residual is relabeled as physical supply.
    ///
    /// This is the existing endpoint-coupled, generally first-order scheme.
    /// Refine time steps for physical accuracy; its convergence gates are not
    /// temporal-error estimates. Reservoirs have prescribed states and no finite
    /// inventory here. No mass transfer, resolved flow, pressure work, moving
    /// geometry or phase-front motion is inferred. No inputs advance on refusal.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_candidate<'m, F: ThermalFluidTransport, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
    ) -> Result<
        CircuitEnthalpyCandidate<'m, ContactThermalSample<ConvectionHeatSample<F>>>,
        CircuitEnthalpyError<ConvectionHeatError<E>>,
    >
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        self.step_resolved(cx, thermal, target_s, circuit, ports, controls, guess_j,
            |state| self.resolve(cx, state))
    }

    /// Evolve one original circuit and the same finite wall/fluid enthalpy
    /// network, including all other authored heat paths and reservoirs.
    ///
    /// Finite convection is internal reciprocal exchange. Its G(Twall,Tfluid)
    /// is re-resolved in every native residual alongside R(T), Joule/contact
    /// heat, conduction and radiation. Thus changing fluid transport or caloric
    /// data changes both temperatures and the actual circuit state.
    /// All accepted inputs remain unchanged on source, budget or energy refusal.
    ///
    /// Cell masses, fixed geometry and prescribed relative flow are retained.
    /// This does not infer fluid pressure work, advection, momentum or phase
    /// front motion. The native endpoint scheme remains generally first order.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_candidate_with_links<'m, F: ThermalFluidTransport, L: ThermalFluidTransport, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        links: &[SourceConvectiveLink],
        maximum_links: usize,
        resolve_links: impl Fn(
            &SourceConvectiveLink, &DepositionSimulation<'_>, Temperature,
        ) -> Result<L, E>,
    ) -> Result<
        CircuitEnthalpyCandidate<'m, ContactThermalSample<ConvectionHeatSample<F, L>>>,
        CircuitEnthalpyError<ConvectionHeatError<E>>,
    >
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        self.step_resolved(cx, thermal, target_s, circuit, ports, controls, guess_j,
            |state| self.resolve_with_links(cx, state, links, maximum_links, &resolve_links))
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn step_resolved<'m, F, L, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'_>)
            -> Result<ContactThermalCoefficients<ConvectionHeatSample<F, L>>, ConvectionHeatError<E>>,
    ) -> Result<
        CircuitEnthalpyCandidate<'m, ContactThermalSample<ConvectionHeatSample<F, L>>>,
        CircuitEnthalpyError<ConvectionHeatError<E>>,
    > {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let routes = self
            .heat
            .heat_destinations(thermal)
            .map_err(ConvectionHeatError::from)
            .map_err(circuit_error)?;
        let destinations = self.heat.resistor_destinations();
        let proposal = contacts::step_candidate(
            thermal,
            target_s,
            CircuitEnthalpyInput {
                circuit,
                ports,
                destinations: &destinations,
            },
            &routes,
            controls,
            guess_j,
            resolve,
            || cx.checkpoint().is_err(),
        );
        let proposal = match proposal {
            Err(CircuitEnthalpyError::Source(ConvectionHeatError::Cancelled)) => {
                return Err(CircuitEnthalpyError::Cancelled);
            }
            other => other?,
        };
        let report = proposal.report();
        let internal_heat_j =
            (target_s - thermal.time_s()) * report.source.material.heat.radiation_net_power_w();
        let residual_j = report.total_energy_residual_j + internal_heat_j;
        if !internal_heat_j.is_finite()
            || !residual_j.is_finite()
            || residual_j.abs() > controls.total_energy_tolerance_j
        {
            return Err(CircuitEnthalpyError::Source(
                ConvectionHeatError::ClosedEnergy {
                    residual_j,
                    tolerance_j: controls.total_energy_tolerance_j,
                },
            ));
        }
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(proposal)
    }
}

fn circuit_error<E>(error: ConvectionHeatError<E>) -> CircuitEnthalpyError<ConvectionHeatError<E>> {
    match error {
        ConvectionHeatError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(other),
    }
}
