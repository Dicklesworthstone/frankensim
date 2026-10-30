//! Nonlinear spatial current and field storage through the shared heat owner.
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
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared
        .field_material_network()
        .map_err(assembly_heat_error)?;
    let layout = ConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    if regions.len() != circuit.sources().len() || regions.len() > assembly.max_cells {
        return Err(map_conduction_error(invalid(
            "one bounded heat route per original circuit triangle required",
        )));
    }
    let mut field_axes = Vec::new();
    for region in regions {
        check(&cancel).map_err(map_conduction_error)?;
        let source = circuit
            .sources()
            .get(region.triangle)
            .ok_or_else(|| invalid("thermal route names an absent circuit conductor"))
            .map_err(map_conduction_error)?;
        for (axis, _) in source.axes() {
            field_axes.push((region.cell, axis.as_str()));
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &field_axes)
        .map_err(assembly_heat_error)?;
    let destinations = prepared.resistor_destinations();
    let next = layout.step_response_circuit_resolved_bound(
        thermal,
        target_s,
        circuit,
        ports,
        &destinations,
        bulk.resistors,
        contacts,
        Some(fields),
        assembly,
        field_controls,
        response,
        controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(map_conduction_error(invalid(
                    "spatial field-storage circuit cannot omit series electrical-contact work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
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
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}
