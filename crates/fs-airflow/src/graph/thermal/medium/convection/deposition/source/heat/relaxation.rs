//! Original source convection on the complete material-memory trajectory.
//! This only composes existing source evaluations with the existing joint step.
use super::*;
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::{
    FrameRelaxationThermalSystem, RelaxationHeatReport,
    electrical::RelaxationElectricalInput,
    transport::{RelaxationHeatInput, RelaxationTransportError},
};
use fs_couple::frame_vibration::thermal::FrameThermalSolve;
use fs_exec::CancelGate;

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
