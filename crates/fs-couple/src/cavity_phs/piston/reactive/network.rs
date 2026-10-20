//! One reaction update followed by the existing reciprocal circuit/wall/gas solve.
//!
//! Chemical formation energy stays in the original gas internal-energy store.
//! Every nonlinear thermal/electrical trial sees the same reacted composition;
//! it never consumes another portion of the accepted chemical inventory.

use super::{
    GasPiston, GasReaction, GasReactionError, GasReactionSolve, GasReactionStep, PistonInput,
};
use crate::cavity_phs::piston::thermal::ConductanceSample;
use crate::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use crate::cavity_phs::piston::thermal::boundary::network::{
    NetworkGasCandidate, NetworkGasDrive, NetworkGasError, NetworkGasReport, NetworkGasSolve,
};
use crate::circuit_enthalpy::material::{
    MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyInput, CircuitThermalCoefficients,
};
use crate::material_circuit::MaterialCircuit;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_thermal::enthalpy::DepositionSimulation;

/// Original owners' refusals, before publishing any chemical or physical history.
#[derive(Debug)]
pub enum ReactiveNetworkError<E, F> {
    Invalid(&'static str),
    Reaction(GasReactionError),
    Network(CircuitEnthalpyError<NetworkGasError<E, F>>),
    Cancelled,
}
impl<E, F> From<GasReactionError> for ReactiveNetworkError<E, F> {
    fn from(error: GasReactionError) -> Self {
        if error == GasReactionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Reaction(error)
        }
    }
}
impl<E, F> From<CircuitEnthalpyError<NetworkGasError<E, F>>> for ReactiveNetworkError<E, F> {
    fn from(error: CircuitEnthalpyError<NetworkGasError<E, F>>) -> Self {
        if matches!(error, CircuitEnthalpyError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Network(error)
        }
    }
}
impl<E: core::fmt::Debug, F: core::fmt::Debug> core::fmt::Display for ReactiveNetworkError<E, F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reactive gas/network coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static, F: std::error::Error + 'static> std::error::Error
    for ReactiveNetworkError<E, F>
{
}

/// Original material/chemical records and an independently measured full balance.
pub struct ReactiveNetworkReport<S, I> {
    /// Fixed-volume reaction at the accepted initial state, exactly once per dt.
    pub reaction: GasReactionStep,
    /// Original coupled electrical, thermal, gas-boundary and motion reports.
    /// Their stage account starts after the chemical reconstruction.
    pub network: NetworkGasReport<S, I>,
    /// Gas internal-energy change across chemistry AND the coupled solve [J].
    pub gas_energy_change_j: f64,
    /// dU_gas + dH_cells + dH_circuit + dK + exterior pressure work,
    /// minus external force work, electrical supply and external heat [J].
    /// Chemical formation energy, Joule heat and reciprocal wall heat are
    /// internal transfers and receive no second supply credit.
    pub total_energy_residual_j: f64,
}

/// The original finite thermal, circuit and gas/motion owners returned together.
pub struct ReactiveNetworkCandidate<'m, S, I, P = GasPiston> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    piston: P,
    report: ReactiveNetworkReport<S, I>,
}
impl<'m, S, I, P> ReactiveNetworkCandidate<'m, S, I, P> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    pub const fn report(&self) -> &ReactiveNetworkReport<S, I> {
        &self.report
    }
    /// Publish all endpoints; all reported internal heat is already applied.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        MaterialCircuit,
        P,
        ReactiveNetworkReport<S, I>,
    ) {
        (self.thermal, self.circuit, self.piston, self.report)
    }
}
impl<'m, S, I> ReactiveNetworkCandidate<'m, S, I> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston) -> P,
    ) -> ReactiveNetworkCandidate<'m, S, I, P> {
        ReactiveNetworkCandidate {
            thermal: self.thermal,
            circuit: self.circuit,
            piston: bind(self.piston),
            report: self.report,
        }
    }
}

impl GasPiston {
    fn network_reaction<E, F>(
        &self,
        reaction: &GasReaction,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
        drive: NetworkGasDrive,
        kinetics: GasReactionSolve,
        controls: NetworkGasSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, GasReactionStep), ReactiveNetworkError<E, F>> {
        super::poll(cancelled)?;
        if !super::positive(drive.input.duration_s)
            || !target_s.is_finite()
            || drive.input.duration_s != target_s - thermal.time_s()
            || !super::positive(controls.total_energy_tolerance_j)
        {
            return Err(ReactiveNetworkError::Invalid(
                "reaction and finite network need one positive physical interval and energy tolerance",
            ));
        }
        let (gas, report) = reaction
            .parcel_candidate(&self.gas, drive.input.duration_s, kinetics, cancelled)?
            .into_parts();
        let mut reacted = self.clone();
        reacted.gas = gas;
        Ok((reacted, report))
    }

    fn finish_reactive_network<'m, S, I, E, F>(
        &self,
        initial: &DepositionSimulation<'_>,
        reaction: GasReactionStep,
        candidate: NetworkGasCandidate<'m, S, I>,
        drive: NetworkGasDrive,
        controls: NetworkGasSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkCandidate<'m, S, I>, ReactiveNetworkError<E, F>> {
        let (thermal, circuit, mut piston, network) = candidate.into_parts();
        let gas_energy = self.gas.mass_kg()
            * (piston.gas.state().specific_internal_energy_j_kg()
                - self.gas.state().specific_internal_energy_j_kg());
        let cell_energy = thermal.energy_j() - initial.energy_j();
        let electrical = &network.network.electrical;
        let motion = &network.boundary.piston;
        let total = ((electrical.stored_after_j - electrical.stored_before_j) + cell_energy)
            + gas_energy
            + motion.kinetic_energy_change_j
            + motion.exterior_pressure_work_j
            - motion.external_force_work_j
            - electrical.native.supplied
            - network.coupled.external_network_heat_j
            - drive.input.heat_j;
        if !gas_energy.is_finite()
            || !cell_energy.is_finite()
            || !total.is_finite()
            || total.abs() > controls.total_energy_tolerance_j
        {
            return Err(ReactiveNetworkError::Invalid(
                "reaction, circuit, finite thermal cells, gas and motion do not close total energy",
            ));
        }
        // Replace the tentative stage's account with the full chemical interval.
        // The chamber heat includes exchanged wall heat; the enclosing balance
        // above includes the reciprocal thermal debit instead of a supply credit.
        let gas_input = PistonInput {
            heat_j: motion.gas.supplied_heat_j,
            ..drive.input
        };
        self.record_reactive_account(
            &mut piston,
            gas_input,
            motion,
            &reaction,
            gas_energy,
            controls.total_energy_tolerance_j,
        )?;
        super::poll(cancelled)?;
        Ok(ReactiveNetworkCandidate {
            thermal,
            circuit,
            piston,
            report: ReactiveNetworkReport {
                reaction,
                network,
                gas_energy_change_j: gas_energy,
                total_energy_residual_j: total,
            },
        })
    }

    /// React once and couple the products to the existing circuit/enthalpy solver.
    ///
    /// Actual product gas heat capacity, pressure and conductivity drive the
    /// original wall heat/momentum solve. Thermal endpoints update the original
    /// material resistance and conductivity on every network residual trial.
    /// All trials start from the same accepted electrical/thermal histories and
    /// the same reacted species: chemistry is never repeated by Newton probes.
    ///
    /// Both stages approximate ONE duration. This is first-order homogeneous
    /// reaction/network splitting, not a flame, shock, stiff mechanism or engine
    /// cycle model. The thermal geometry, circuit topology/C/L and rigid piston
    /// inertia retain the original owners' fixed-state assumptions. No external
    /// heat argument may include reaction energy or already deposited Joule heat.
    /// Original source support, consumption, entropy, work and cancellation gates
    /// remain in force; failure publishes none of the endpoints.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn reactive_network_candidate<'m, S, I, E, F>(
        &self,
        reaction: &GasReaction,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        electrical: CircuitEnthalpyInput<'_>,
        drive: NetworkGasDrive,
        kinetics: GasReactionSolve,
        controls: NetworkGasSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
        resolve_port: impl Fn(
            &DepositionSimulation<'m>,
            &S,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, F>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkCandidate<'m, S, I>, ReactiveNetworkError<E, F>> {
        let (reacted, chemistry) = self.network_reaction(
            reaction,
            thermal,
            target_s,
            drive,
            kinetics,
            controls,
            &mut cancelled,
        )?;
        let network = reacted.network_candidate(
            thermal,
            target_s,
            electrical,
            drive,
            controls,
            guess_j,
            resolve,
            resolve_port,
            &mut cancelled,
        )?;
        self.finish_reactive_network(thermal, chemistry, network, drive, controls, &mut cancelled)
    }

    /// Same-card electrical, caloric and gas-wall transport without a custom callback.
    ///
    /// The network pins the original conductivity/resistivity claims and reuses
    /// the complete endpoint bundle for wall conduction. The gas conductivity
    /// comes from the actual reacted mixture's caloric midpoint. Source identity,
    /// physical state and energy therefore travel together through all owners.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn reactive_material_network_candidate<'m>(
        &self,
        reaction: &GasReaction,
        binding: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        drive: NetworkGasDrive,
        kinetics: GasReactionSolve,
        controls: NetworkGasSolve,
        geometry: GasConductionGeometry,
        guess_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        ReactiveNetworkCandidate<'m, MaterialNetworkSample, ResolvedMaterialStatePoint>,
        ReactiveNetworkError<MaterialNetworkError, MaterialNetworkError>,
    > {
        let (reacted, chemistry) = self.network_reaction(
            reaction,
            thermal,
            target_s,
            drive,
            kinetics,
            controls,
            &mut cancelled,
        )?;
        let network = reacted.material_network_candidate_with_bases(
            binding,
            thermal,
            target_s,
            circuit,
            ports,
            electrical_bases,
            drive,
            controls,
            geometry,
            guess_j,
            &mut cancelled,
        )?;
        self.finish_reactive_network(thermal, chemistry, network, drive, controls, &mut cancelled)
    }
}
