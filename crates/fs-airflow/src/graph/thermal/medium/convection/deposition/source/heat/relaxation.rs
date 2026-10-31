//! Original source convection on the complete material-memory trajectory.
//! This only composes existing source evaluations with the existing joint step.
/// Pressure-derived flow, convection and hydraulic heat on the same endpoint.
pub mod pressure;

use super::*;
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::{
    FrameRelaxationThermalSystem, RelaxationHeatReport,
    electrical::RelaxationElectricalInput,
    transport::{RelaxationHeatInput, RelaxationTransportError},
};
use fs_couple::frame_vibration::thermal::FrameThermalSolve;
use fs_exec::CancelGate;
use fs_conduction::transient::deposition::advection::Advected;

/// The accepted reservoir-fluid sources accompany the original complete
/// mechanical/electrical/caloric report in `transport_source`. They are the
/// final accepted trial, not a fresh post-publication evaluation.
pub type ConvectiveRelaxationReport<F> = RelaxationHeatReport<ResolvedCellConvectionSources<F>>;
pub type ConvectiveRelaxationError<E> = RelaxationTransportError<ConvectionHeatError<E>>;

impl<R> SourceConvectionHeatNetwork<'_, '_, '_, '_, R> {
    /// Advance material relaxation, fluid memory, finite-layer enthalpy and
    /// optional circuit/field storage with source-resolved reservoir convection.
    ///
    /// EVERY caloric trial evaluates the original fluid at the declared bulk
    /// or film temperature. Its actual rho/cp/k/mu (and expansion for buoyancy)
    /// feed the existing correlation and original exposed-area mapping. The
    /// resulting heat changes local temperature, R(T)/C(T)/L(T) where bound,
    /// and the same layers' source tau(T), which acts next mechanical interval.
    /// No independent wall/coolant temperature, heat step or damping fit is added.
    ///
    /// `M` is the retained mechanical-fluid source type; `F` is the heat-transfer
    /// source type. They may differ physically. This entry does NOT silently
    /// retune frozen Stokes inertia/drag with a film temperature or heat-source
    /// swap. Bind each material/geometry explicitly according to its own model.
    ///
    /// The reservoirs, relative flow and pressure/composition are prescribed;
    /// this is not a finite-fluid energy, momentum or moving-boundary solver.
    /// Correlation support and heat-surface coverage remain original admission
    /// gates. An authored ambient h on a source-convecting cell refuses rather
    /// than contributing duplicate cooling. Source pins belong to the caller's
    /// retained fluid law, not a fabricated anonymous transport record.
    ///
    /// Convection is external heat, radiation is internal transfer, and Joule,
    /// Maxwell and reversible field heat retain separate accounts. The generic
    /// caloric entropy report leaves callback heat unclassified at the receiver;
    /// `transport_source` retains the physical reservoir and heat-transfer data.
    /// It does not claim that Q/T_cell is an irreversible reservoir entropy.
    /// Source failures retain E, and cancellation/budget/domain failures publish
    /// none of the mechanical, thermal, source, electrical or fluid histories.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_relaxation<F: ThermalFluidTransport, M, E>(
        &self,
        cx: &Cx<'_>,
        owner: &mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        electrical: Option<RelaxationElectricalInput<'_>>,
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<ConvectiveRelaxationReport<F>, ConvectiveRelaxationError<E>>
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        owner.step_heat_network_with(
            thermal, external_modal_forces, RelaxationHeatInput { network: self.heat, cx },
            electrical, controls, energy_guess_j, gate, |trial| {
                let ContactThermalCoefficients { bulk, contacts } = self.resolve(cx, trial)?;
                let ConvectionHeatSample { heat, convection, .. } = bulk.source;
                Ok(ContactThermalCoefficients { contacts, bulk: CircuitThermalCoefficients {
                    resistors: bulk.resistors, sources_w: bulk.sources_w,
                    conductances_w_k: bulk.conductances_w_k, source: (heat, convection),
                } })
            },
        )
    }
}

/// Original accepted sources: reservoirs exchange external heat; finite links
/// exchange internal heat between existing caloric inventories. Neither payload
/// introduces a second temperature or advances after the joint publication.
#[derive(Debug)]
pub struct FiniteRelaxationConvection<F, L = F> {
    pub reservoirs: ResolvedCellConvectionSources<F>,
    pub finite: ResolvedLinkConvection<L>,
}
pub type FiniteConvectiveRelaxationReport<F, L = F> = RelaxationHeatReport<FiniteRelaxationConvection<F, L>>;

impl<R> SourceConvectionHeatNetwork<'_, '_, '_, '_, R> {
    /// Continue the complete material/circuit/frame trajectory with FINITE
    /// coolant cells, optionally alongside separately prescribed reservoirs.
    /// Every trial reuses `resolve_with_links`: each wall-fluid interface sees
    /// both actual caloric temperatures and its declared bulk/film convention.
    ///
    /// Prepare the thermal owner and heat network with their existing
    /// `*_with_link_conductances` entry points, using the original initial
    /// source evaluation. Fourier bindings must exclude these links. Their
    /// native areas, coverage, equal/opposite heat, cell masses and caloric
    /// laws stay authoritative. Several walls can share one finite coolant.
    /// Empty reservoir `cells` gives a closed finite heat-exchange system.
    ///
    /// Changing coolant capacity or transport changes its evolving temperature
    /// and thus the solids' temperature, source relaxation and circuit response.
    /// Finite exchange enters only the native conductance map, never a duplicate
    /// pair of heater powers or an external bath. Native conductive entropy
    /// accounting already includes this reciprocal transfer once.
    ///
    /// Retain the original coolant laws and contexts in `resolve_links`; generic
    /// transport cannot prove their correspondence to a caloric material. Flow,
    /// pressure/composition and geometry are prescribed: no advection, mass or
    /// momentum transfer, phase-front motion or variable Stokes loading is
    /// inferred. The original rate/thermal split still requires time refinement.
    /// Typed source failures and either cancellation authority preserve ALL
    /// accepted histories, including finite coolant and field charge/flux.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_relaxation_with_links<F: ThermalFluidTransport, L: ThermalFluidTransport, M, E>(
        &self, cx:&Cx<'_>,
        owner:&mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal:&mut DepositionSimulation<'_>, external_modal_forces:&[f64],
        electrical:Option<RelaxationElectricalInput<'_>>, controls:FrameThermalSolve,
        energy_guess_j:&[f64], gate:&CancelGate, links:&[SourceConvectiveLink],
        maximum_links:usize,
        resolve_links:impl Fn(&SourceConvectiveLink,&DepositionSimulation<'_>,Temperature)->Result<L,E>,
    )->Result<FiniteConvectiveRelaxationReport<F,L>,ConvectiveRelaxationError<E>>
    where R:Fn(&SourceConvectiveCell,&DepositionSimulation<'_>,Temperature)->Result<F,E> {
        self.step_relaxation_links_flow(cx,owner,thermal,external_modal_forces,electrical,
            controls,energy_guess_j,gate,links,maximum_links,resolve_links,None)
    }

    /// The same finite coolant/convection/circuit/material trajectory with
    /// prescribed balanced mass-carried enthalpy. `flow` uses original finite
    /// coolant cell IDs and caloric material identities. Keep this graph and
    /// its inlet states fixed over the interval. No accepted observer is stepped.
    ///
    /// Both convective transfer and advective donor enthalpy see the SAME trial
    /// temperatures. Recirculation is internal; only inlet/outlet enthalpy is
    /// external supply. The native implicit residual applies every exchange
    /// once. Flow through bound mechanical material inventories refuses.
    ///
    /// Mass rates, relative surface speeds and pressure/composition remain
    /// prescribed and must represent the same physical flow where applicable.
    /// This does not solve pumps, pressure, fluid momentum or moving mass.
    /// It preserves the original upwind spatial and first-order time limits.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_advected_relaxation_with_links<F: ThermalFluidTransport, L: ThermalFluidTransport, M, E>(
        &self, cx:&Cx<'_>,
        owner:&mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal:&mut DepositionSimulation<'_>, external_modal_forces:&[f64],
        electrical:Option<RelaxationElectricalInput<'_>>, controls:FrameThermalSolve,
        energy_guess_j:&[f64], gate:&CancelGate, links:&[SourceConvectiveLink],
        maximum_links:usize,
        resolve_links:impl Fn(&SourceConvectiveLink,&DepositionSimulation<'_>,Temperature)->Result<L,E>,
        flow:&Advected<()>,
    )->Result<FiniteConvectiveRelaxationReport<F,L>,ConvectiveRelaxationError<E>>
    where R:Fn(&SourceConvectiveCell,&DepositionSimulation<'_>,Temperature)->Result<F,E> {
        self.step_relaxation_links_flow(cx,owner,thermal,external_modal_forces,electrical,
            controls,energy_guess_j,gate,links,maximum_links,resolve_links,Some(flow))
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn step_relaxation_links_flow<F: ThermalFluidTransport, L: ThermalFluidTransport, M, E>(
        &self, cx:&Cx<'_>,
        owner:&mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal:&mut DepositionSimulation<'_>, external_modal_forces:&[f64],
        electrical:Option<RelaxationElectricalInput<'_>>, controls:FrameThermalSolve,
        energy_guess_j:&[f64], gate:&CancelGate, links:&[SourceConvectiveLink],
        maximum_links:usize,
        resolve_links:impl Fn(&SourceConvectiveLink,&DepositionSimulation<'_>,Temperature)->Result<L,E>,
        flow:Option<&Advected<()>>,
    )->Result<FiniteConvectiveRelaxationReport<F,L>,ConvectiveRelaxationError<E>>
    where R:Fn(&SourceConvectiveCell,&DepositionSimulation<'_>,Temperature)->Result<F,E> {
        let resolve = |trial:&DepositionSimulation<'_>| {
                let ContactThermalCoefficients {bulk,contacts}=self.resolve_with_links(
                    cx,trial,links,maximum_links,&resolve_links)?;
                let ConvectionHeatSample {heat,convection,finite_convection}=bulk.source;
                let finite=finite_convection.ok_or(ConvectionHeatError::Convection(
                    SourceCellConvectionError::Convection(ConvectionSourceError::Invalid(
                        "complete finite convection evaluation is absent"))))?;
                Ok(ContactThermalCoefficients {contacts,bulk:CircuitThermalCoefficients {
                    resistors:bulk.resistors,sources_w:bulk.sources_w,conductances_w_k:bulk.conductances_w_k,
                    source:(heat,FiniteRelaxationConvection {reservoirs:convection,finite}),
                }})
            };
        match flow {
            Some(flow) => owner.step_heat_network_with_advection(thermal,external_modal_forces,
                RelaxationHeatInput {network:self.heat,cx},electrical,controls,energy_guess_j,gate,flow,resolve),
            None => owner.step_heat_network_with(thermal,external_modal_forces,
                RelaxationHeatInput {network:self.heat,cx},electrical,controls,energy_guess_j,gate,resolve),
        }
    }
}
