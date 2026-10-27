//! Voltage-driven source permeability, field entropy and finite caloric storage.
//!
//! The native circuit/field step owns Maxwell/Kirchhoff evolution and exposes
//! signed magnetic parameter work P. The shared magnetic caloric owner adds
//! Q_rev=-P-delta(T*S) to the original thermal residual. No empirical core loss,
//! surrogate inductance or second thermal integrator is introduced.

use super::*;
use crate::circuit_enthalpy::{
    finite,
    magnetic::{MagneticRoutes, MagneticThermalNetwork, MagneticThermalSample},
    material::MaterialNetworkSample,
};
use fs_em::planar::nonlinear::PlanarNonlinearMagnetostatics;
use fs_material::{
    conductor::ResolvedConductorState, field_storage::admit_magnetic_temperature_refresh,
};

/// Physical internal-energy account, independent of the native free-energy
/// balance. Joule heating and reversible magnetic heat are internal transfers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MagneticCaloricBalance {
    /// Authored boundary, reservoir and independent source heat [J].
    pub external_heat_j: f64,
    /// Delta(field U + native circuit/end storage + caloric energy) + native
    /// numerical dissipation - external electrical work - external heat [J].
    pub total_energy_residual_j: f64,
}

impl WindingThermalNetwork<'_> {
    /// Evolve voltage-driven winding current, sourced mu(T), conductor R(T)
    /// and the original finite caloric cells in one joint residual.
    ///
    /// Every thermally mapped magnetic triangle needs its original derivative
    /// law and a field envelope with positive combined caloric/field heat
    /// capacity. All magnetic regions use linear isotropic permeability with
    /// the shared law's positive affine mu(T) over its declared interval here;
    /// H(B,T), hysteresis and remanence require their own thermal derivatives.
    /// Fixed external regions must have temperature-independent permeability.
    /// Geometry, end inductances, topology, source selections and fixed material
    /// coordinates remain unchanged. The native retained-field solver rebuilds
    /// its operator at each proposed temperature, retaining the accepted nodal
    /// field, circuit charge/flux, original material identities and clock.
    ///
    /// `graph_material` optionally binds all ordinary graph resistors through
    /// the existing material transport owner; otherwise their authored values
    /// remain fixed. Its electrical bases and destinations have the same
    /// semantics as `step_nonlinear_material_candidate`.
    ///
    /// This is first-order backward Euler. Its numerical loss is reported but
    /// never deposited as heat. Each trial starts from the same accepted state;
    /// a source refusal, cancellation or failed local/global balance returns no
    /// partial physical endpoint. `caloric_balance()` reports physical closure;
    /// the ordinary report retains native free-energy/parameter-work closure.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_thermal_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &NonlinearWindingCircuit<'m>,
        external_inputs: &[f64],
        bindings: &[MagneticThermalLaw],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearWindingEnthalpyCandidate<'t, 'm>, Error> {
        if graph_material.is_none() && !circuit.graph_resistor_materials().is_empty() {
            return Err(invalid(
                "source-bound graph resistors require their complete thermal material bindings",
            ));
        }
        self.step_nonlinear_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            graph_material,
            Some((bindings, assembly)),
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }
}

// The same caloric, volume, conductivity and complete-source refresh as every
// other magnetic consumer. Only magnetic storage admission changes: its
// temperature work is now owned by the caller's shared caloric binding.
#[allow(clippy::too_many_arguments)]
pub(super) fn materials_at(
    layout: &MagneticThermalNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    problem: &PlanarNonlinearMagnetostatics,
    conductors: &BTreeMap<usize, &ResolvedConductorState>,
    sample: &MaterialNetworkSample,
    routes: &MagneticRoutes,
    cancelled: &RefCell<impl FnMut() -> bool>,
) -> Result<
    (
        Vec<(usize, ResolvedConductorState)>,
        Vec<MagneticThermalSample>,
    ),
    Error,
> {
    let originals: Vec<_> = problem
        .materials()
        .iter()
        .map(IsotropicMagneticMaterial::material)
        .collect();
    layout.materials_at_using(
        thermal,
        &originals,
        conductors,
        sample,
        routes,
        cancelled,
        |_, original, snapshot, _| {
            admit_magnetic_temperature_refresh(original, snapshot)
                .map_err(|e| CircuitEnthalpyError::Source(WindingThermalError::Magnetic(e)))
        },
    )
}

pub(super) fn balance(
    report: &CircuitEnthalpyReport<WindingThermalSample, NonlinearWindingStep>,
    tolerance_j: f64,
) -> Result<Option<MagneticCaloricBalance>, Error> {
    let Some(source) = &report.source.magnetic_caloric else {
        return Ok(None);
    };
    let step = &report.electrical;
    // The complete native circuit store already includes field free energy,
    // end inductance, capacitors and any original unbound inductors. Replace
    // only the represented field contribution by its physical internal energy.
    let correction = finite(
        finite(source.field_internal_energy_j[1] - step.field.stored_after_j)?
            - finite(source.field_internal_energy_j[0] - step.field.stored_before_j)?,
    )?;
    let external_heat_j = finite(report.non_electrical_heat_j - source.reversible_heat_j)?;
    let total_energy_residual_j = finite(
        finite(step.stored_after_j - step.stored_before_j)?
            + correction
            + step.numerical_dissipation_j
            + report.thermal.stored_energy_change_j
            - step.supplied_work_j
            - external_heat_j,
    )?;
    if total_energy_residual_j.abs() > tolerance_j {
        return Err(invalid(
            "native circuit, magnetic internal energy and caloric storage do not close external supply",
        ));
    }
    Ok(Some(MagneticCaloricBalance {
        external_heat_j,
        total_energy_residual_j,
    }))
}
