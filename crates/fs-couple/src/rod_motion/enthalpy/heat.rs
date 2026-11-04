//! Complete source heat transfer on the original native rod transaction.
//!
//! This adapter replaces only the transport evaluation in the existing finite
//! enthalpy residual. Native motion, full layer Psi/S/U, reversible exchange,
//! intrinsic axial/bending/torsional heat and endpoint source refresh keep their one owner.
use super::*;
use crate::circuit_enthalpy::{CircuitThermalCoefficients, contacts::heat::{
    HeatAssembly, HeatNetworkError,
}};
use fs_conduction::transient::deposition::radiation::RadiationCellSources;
use fs_exec::Cx;
use std::cell::RefCell;

/// Retain the concrete complete-assembly error and every original rod error.
/// Cancellation remains `Rod(CircuitEnthalpyError::Cancelled)`.
#[derive(Debug)]
pub enum RodHeatError<E = HeatNetworkError> {
    Heat(E),
    Rod(RodThermoelasticError),
}
impl<E> From<RodThermoelasticError> for RodHeatError<E> {
    fn from(error: RodThermoelasticError) -> Self { Self::Rod(error) }
}
impl<E: core::fmt::Debug> core::fmt::Display for RodHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rod and complete material heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RodHeatError<E> {}

impl<I: RodMotionInertia> RodThermoelasticSystem<'_, '_, '_, I> {
    /// Advance actual rod dynamics with complete source contact, layer,
    /// radiation, finite-coolant and reservoir heat on the SAME finite cells.
    /// Each original thermal trial resolves the assembly once, retains its
    /// full sample, and refreshes E/G/integrated expansion at the native endpoint.
    /// Optional intrinsic axial, bending and torsional viscosity enter those
    /// same cells once. Pass [`RodMaterialViscosity`] for any combination.
    /// Existing [`RodAngularViscosity`] and `Option<RodTorsionalViscosity>`
    /// arguments remain supported.
    ///
    /// The prepared assembly must borrow the original material cell/link slices.
    /// Additional finite links require cold admission through
    /// [`Self::try_new_with_link_conductances`]; preparation conductances are
    /// never reused as evolving physics. Electrical storage, resistor and
    /// electrical-contact work refuse because this path has no circuit owner.
    ///
    /// Thermal contact, fluid-boundary and radiation geometry remains explicitly
    /// fixed at its reference configuration. Native rod displacement does not
    /// infer changed area, view factors or coolant flow. The original first-order
    /// thermoelastic split requires time refinement. A source, convergence,
    /// energy or cancellation failure publishes neither mechanical nor heat state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_assembly<'v, A: HeatAssembly + ?Sized>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: &[RodNodalWrench],
        assembly: &A,
        cx: &Cx<'_>,
        viscosity: impl Into<RodMaterialViscosity<'v>>,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> core::result::Result<RodThermoelasticReport<A::Sample>, RodHeatError<A::Error>> {
        self.step_heat_with_loads(thermal, target_s, StageLoads::World(loads), assembly, cx,
            viscosity.into(), mechanical_controls, controls, energy_guess_j, gate)
    }

    /// The same complete finite-heat transaction with original material-point
    /// forces and body/world-held free couples. Spatial load work and native
    /// integration error retain their existing accounts; neither becomes heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_heat_assembly<'v, A: HeatAssembly + ?Sized>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: &PreparedRodLoads,
        assembly: &A,
        cx: &Cx<'_>,
        viscosity: impl Into<RodMaterialViscosity<'v>>,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> core::result::Result<RodThermoelasticReport<A::Sample>, RodHeatError<A::Error>> {
        self.step_heat_with_loads(thermal, target_s, StageLoads::Spatial(loads), assembly, cx,
            viscosity.into(), mechanical_controls, controls, energy_guess_j, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_heat_with_loads<A: HeatAssembly + ?Sized>(
        &mut self, thermal: &mut DepositionSimulation<'_>, target_s: f64,
        loads: StageLoads<'_>, assembly: &A, cx: &Cx<'_>,
        viscosity: RodMaterialViscosity<'_>,
        mechanical_controls: RodMotionControls, controls: RodThermoelasticSolve,
        energy_guess_j: &[f64], gate: &CancelGate,
    ) -> core::result::Result<RodThermoelasticReport<A::Sample>, RodHeatError<A::Error>> {
        let trial = AssemblyTrial { assembly, cx, failure: RefCell::new(None) };
        trial.checkpoint(gate)?;
        let heat = assembly.prepared_heat_network();
        let material = heat.field_storage_material_network();
        if !core::ptr::eq(self.network_source.cells, material.cells)
            || !core::ptr::eq(self.network_source.links, material.links)
            || !material.resistors.is_empty()
            || !heat.resistor_destinations().is_empty() {
            return Err(invalid("rod heat assembly needs original cell/link slices and no unowned electrical graph").into());
        }
        let routes = heat.heat_destinations(thermal).map_err(|error| match error {
            HeatNetworkError::Cancelled => {
                gate.request();
                RodHeatError::Rod(CircuitEnthalpyError::Cancelled)
            }
            other => RodHeatError::Heat(A::Error::from(other)),
        })?;
        if !routes.is_empty() {
            return Err(invalid("rod heat assembly cannot discard electrical contact work").into());
        }
        let result = self.step_with_transport(thermal, target_s, loads,
            viscosity, &trial,
            mechanical_controls, controls, energy_guess_j, gate);
        match result {
            Err(CircuitEnthalpyError::Invalid(SOURCE_REFUSAL)) => match trial.failure.into_inner() {
                Some(error) => Err(RodHeatError::Heat(error)),
                None => Err(invalid("complete rod heat source refusal lost its original cause").into()),
            },
            other => other.map_err(RodHeatError::from),
        }
    }
}

// Only transport and its retained receipt vary. The parent owns the sole
// mechanical step, cell-energy solve, thermodynamic audit and publication.
pub(super) trait Transport {
    type Sample;
    fn checkpoint(&self, gate: &CancelGate) -> Result<()> { check(gate) }
    fn resolve(&self, prepared: &PreparedMaterialNetwork<'_, '_>,
        thermal: &DepositionSimulation<'_>, gate: &CancelGate,
    ) -> Result<CircuitThermalCoefficients<Self::Sample>>;
    fn radiation<'s>(&self, _source: &'s Self::Sample) -> Option<&'s RadiationCellSources> { None }
}
pub(super) struct Plain;
impl Transport for Plain {
    type Sample = MaterialNetworkSample;
    fn resolve(&self, prepared: &PreparedMaterialNetwork<'_, '_>,
        thermal: &DepositionSimulation<'_>, gate: &CancelGate,
    ) -> Result<CircuitThermalCoefficients<Self::Sample>> {
        prepared.resolve(thermal, || gate.is_requested()).map_err(network)
    }
}

const SOURCE_REFUSAL: &str = "original complete heat source refused rod trial";
struct AssemblyTrial<'p, 'scope, A: HeatAssembly + ?Sized> {
    assembly: &'p A,
    cx: &'p Cx<'scope>,
    failure: RefCell<Option<A::Error>>,
}
impl<A: HeatAssembly + ?Sized> Transport for AssemblyTrial<'_, '_, A> {
    type Sample = A::Sample;
    fn checkpoint(&self, gate: &CancelGate) -> Result<()> {
        if self.cx.checkpoint().is_err() { gate.request(); }
        check(gate)
    }
    fn resolve(&self, _prepared: &PreparedMaterialNetwork<'_, '_>,
        thermal: &DepositionSimulation<'_>, gate: &CancelGate,
    ) -> Result<CircuitThermalCoefficients<Self::Sample>> {
        self.checkpoint(gate)?;
        self.failure.borrow_mut().take();
        let resolved = match self.assembly.resolve_heat(self.cx, thermal) {
            Ok(resolved) => resolved,
            Err(error) => match assembly_error(error, gate) {
                RodHeatError::Rod(error) => return Err(error),
                RodHeatError::Heat(error) => {
                    *self.failure.borrow_mut() = Some(error);
                    return Err(invalid(SOURCE_REFUSAL));
                }
            },
        };
        self.checkpoint(gate)?;
        if !resolved.contacts.is_empty() || !resolved.bulk.resistors.is_empty() {
            return Err(invalid("rod thermal trial cannot discard electrical storage or Joule/contact work"));
        }
        Ok(resolved.bulk)
    }
    fn radiation<'s>(&self, source: &'s Self::Sample) -> Option<&'s RadiationCellSources> {
        source.as_ref().radiation.as_ref()
    }
}

fn assembly_error<E>(error: CircuitEnthalpyError<E>, gate: &CancelGate) -> RodHeatError<E> {
    let native = match error {
        CircuitEnthalpyError::Source(error) => return RodHeatError::Heat(error),
        CircuitEnthalpyError::Invalid(message) => CircuitEnthalpyError::Invalid(message),
        CircuitEnthalpyError::Circuit(error) => error.into(),
        CircuitEnthalpyError::Thermal(error) => error.into(),
        CircuitEnthalpyError::Newton(error) => CircuitEnthalpyError::Newton(error),
        CircuitEnthalpyError::NotConverged(report) => CircuitEnthalpyError::NotConverged(report),
        CircuitEnthalpyError::Interrupted(progress) => CircuitEnthalpyError::Interrupted(progress),
        CircuitEnthalpyError::CircuitBudget => CircuitEnthalpyError::CircuitBudget,
        CircuitEnthalpyError::Cancelled => CircuitEnthalpyError::Cancelled,
    };
    if matches!(&native, CircuitEnthalpyError::Cancelled) { gate.request(); }
    RodHeatError::Rod(native)
}

pub(super) fn cell_power_w(radiation: Option<&RadiationCellSources>, cell: CellId) -> f64 {
    radiation.and_then(|r| r.sources_w().get(&cell)).copied().unwrap_or(0.0)
}
pub(super) fn exchange(radiation: Option<&RadiationCellSources>,
    thermal: &DepositionSimulation<'_>, dt: f64, gate: &CancelGate,
) -> Result<(f64, f64)> {
    let mut entropy = 0.0;
    if let Some(radiation) = radiation {
        for (&cell, &power) in radiation.sources_w() {
            check(gate)?;
            entropy = finite(entropy + dt * power / thermal.temperature_k(cell)?)?;
        }
    }
    Ok((finite(dt * radiation.map_or(0.0, |r| r.net_heat_w()))?, entropy))
}
