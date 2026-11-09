//! Nonlinear spatial current and field storage through the shared heat owner.
mod motion;
pub use motion::step_spatial_moving_conduction_candidate;

use super::*;
use crate::circuit_enthalpy::contacts::heat::HeatAssembly;

/// Original spatial/field failures and the complete assembly's typed failures.
#[derive(Debug)]
pub enum ConductionHeatError<E> {
    Conduction(ConductionThermalError),
    Heat(E),
}
impl<E> From<ConductionThermalError> for ConductionHeatError<E> {
    fn from(error: ConductionThermalError) -> Self {
        Self::Conduction(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ConductionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial field and heat assembly: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ConductionHeatError<E> {}

fn assembly_heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> CircuitEnthalpyError<ConductionHeatError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionHeatError::Heat(other.into())),
    }
}

/// Close native tetrahedral current, electrical-contact work and finite caloric
/// storage with any complete [`HeatAssembly`] on the SAME proposed endpoint.
///
/// The original prepared network retains the material cards, finite inventory
/// and heat destinations. The complete assembly supplies its original sourced
/// thermal contacts, layers, radiation, reservoirs and reciprocal native links.
/// Passive finite cells remain in that inventory: a source-fluid convection
/// assembly can warm a coolant and thereby change actual bulk conductivity,
/// contact resistance and current, without another thermal timeline or a fitted
/// electrical resistor. Its entire typed source sample and errors are retained.
///
/// Native tetrahedral E.J and exact surface voltage-jump work enter the original
/// enthalpy residual once. Every electrical contact requires explicit two-sided
/// source coordinates and heat ownership. Reciprocal finite links and enclosure
/// radiation remain internal exchange; external non-electrical heat is
/// `report.non_electrical_heat_j - dt * report.source.as_ref().radiation_net_power_w()`.
///
/// The drive accepts existing voltage/current-density traces or equipotential
/// total-current electrodes. Under total-current drive, the original refreshed
/// bulk and contact sources determine voltage and local current on every trial;
/// their actual I²R work joins the same finite-cell residual and heat assembly.
///
/// This keeps the existing fixed geometry, matching electrical traces, confined
/// phase charts and field-independent thermal-source requirements. An assembly
/// must resolve the supplied trial state without publishing side effects. Source
/// failure, rejected trials and cancellation publish neither endpoint. Graph
/// electrical bindings refuse because this native field has no owner for their
/// work. This first-order coupling requires temporal and spatial refinement.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_spatial_conduction_candidate<
    't,
    'd,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &SpatialConduction<C>,
    drive: impl Into<SpatialConductionDrive<'d>>,
    contacts: &[contacts::ConductionContactThermal<'_>],
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialConductionEnthalpyCandidate<'t, C, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared
        .thermal_material_network()
        .map_err(assembly_heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(assembly_heat_error)?;
    let layout = SpatialConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let next = layout.step_resolved_candidate(
        thermal,
        target_s,
        problem,
        drive,
        contacts,
        assembly,
        field_controls,
        controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(map_conduction_error(invalid(
                    "3D current cannot silently omit graph electrical work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    let report = next.report();
    closed_balance(
        target_s - thermal.time_s(),
        report.total_energy_residual_j,
        report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(assembly_heat_error)?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}

/// Evolve the full multi-electrode tetrahedral conductor, native fixed C/L
/// history and complete heat assembly at one finite-enthalpy endpoint.
///
/// Original scalar or principal conductivity, electrical contacts and ordinary
/// graph resistors refresh from the SAME proposed material temperature/phase
/// as thermal interfaces, layers, baths, finite links and enclosure radiation.
/// Each native midpoint field/contact/branch loss enters its original receiving
/// volume once. Every trial starts from the accepted charges and fluxes.
///
/// Graph-resistor heat routes belong to the prepared assembly and must remain
/// disjoint from the field's receiving volumes. Spatial contacts retain their
/// explicit two-sided sources and heat partition. Graph series contacts refuse
/// because this circuit has no owner for their additional work. Source-dependent
/// C/L still requires the original reversible field-storage owner.
///
/// The complete typed heat sample and failures are retained. Reciprocal finite
/// links and enclosure radiation remain internal: external non-electrical heat
/// is `report.non_electrical_heat_j - dt * report.source.as_ref().radiation_net_power_w()`.
/// Original closed-enclosure balance, source domains, budgets, event cuts and
/// cancellation apply before either history is published. Geometry is fixed,
/// phase support remains confined to admitted charts, thermal boundary sources
/// must be field-independent, and temporal/spatial refinement is still required.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_spatial_conduction_circuit_candidate<
    't,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &SpatialConductionCircuit<C>,
    ports: &[f64],
    contacts: &[contacts::ConductionContactThermal<'_>],
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    response: SpatialConductanceControls,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialConductionCircuitEnthalpyCandidate<'t, C, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared
        .field_material_network()
        .map_err(assembly_heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(assembly_heat_error)?;
    let layout = SpatialConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let destinations = prepared.resistor_destinations();
    let next = layout.step_circuit_resolved_candidate(
        thermal,
        target_s,
        circuit,
        ports,
        &destinations,
        contacts,
        response,
        assembly,
        field_controls,
        controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(map_conduction_error(invalid(
                    "3D finite circuit cannot omit graph series electrical-contact work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    let report = next.report();
    closed_balance(
        target_s - thermal.time_s(),
        report.total_energy_residual_j,
        report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(assembly_heat_error)?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}

/// Evolve distributed nonlinear current, native electrical contacts, source
/// R(T)/C(T)/L(T) and finite caloric storage with any complete [`HeatAssembly`].
///
/// The assembly's original prepared network supplies source identities and
/// heat destinations. Its complete resolver supplies thermal contacts, layers,
/// radiation, prescribed reservoirs and additional reciprocal native links at
/// the SAME proposed enthalpy endpoint. This includes the original source-fluid
/// convection assembly; finite coolant warming can change actual spatial
/// current, contact loss and subsequent charge/flux history without another
/// thermal or electrical integrator.
///
/// Local E.J and surface Joule loss remain owned by the original spatial solve.
/// Source field rebinding and reversible entropy heat remain owned by the same
/// field driver. The accepted report retains the entire typed assembly sample;
/// physical field U = psi + T s plus all caloric inventories close against
/// external supply. Finite links and enclosure radiation remain internal.
/// No iteration error or enclosure roundoff is relabeled as physical heat.
///
/// The original field-independent thermal-source and boundary requirements
/// also apply to added heat laws: no unresolved E/J/B/H may affect their response
/// or source applicability. Fixed geometry, complete disjoint link ownership,
/// solid affine field support, positive capacity/coordinate envelopes, explicit
/// event cuts and first-order temporal-refinement requirements still apply.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_response_field_circuit_candidate<'t, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &ConductionResponseCircuit,
    ports: &[f64],
    contacts: &[contacts::ConductionContactThermal<'_>],
    fields: &[FieldThermalBinding],
    regions: &[ConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: ConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    ConductionResponseFieldThermalCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let next = super::super::circuit::response::thermal::step_response_circuit_with_heat(
        heat,
        cx,
        thermal,
        target_s,
        circuit,
        ports,
        None,
        &[],
        contacts,
        Some(fields),
        regions,
        volume_tolerance_m3,
        assembly,
        field_controls,
        response,
        controls,
        initial_energy_increments_j,
    )?;
    let internal_heat_j = finite(
        (target_s - thermal.time_s())
            * next
                .report()
                .source
                .network
                .as_ref()
                .radiation_net_power_w(),
    )?;
    let next = next
        .close_field_report(internal_heat_j, controls.total_energy_tolerance_j)
        .map_err(map_conduction_error)?;
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    Ok(next)
}
