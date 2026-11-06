//! Physical solid-conductor electrodes, retained vector memory and finite heat.
//!
//! Voltage control solves the current distribution; total-current control solves
//! each equipotential electrode voltage in the same A–phi system. Every thermal
//! trial restarts the accepted electrode, magnetic and caloric histories. Native
//! electrode supply, impressed-current work, eddy heat and irreversible pinning
//! enter one balance. Numerical dissipation remains outside physical heating.

use super::super::{
    SpatialInductionHeatError, SpatialMagneticRegion,
    heat::{HeatError, close_heat_balance, heat_error},
    terminals::SpatialTerminalThermalCandidate,
};
use super::*;
use crate::circuit_enthalpy::contacts::heat::HeatAssembly;
use fs_em::spatial::vector_potential::induction::hysteretic::{
    SpatialHystereticProblem,
    terminals::{
        SpatialHystereticCurrentInduction, SpatialHystereticCurrentInductionCandidate,
        SpatialHystereticTerminalInduction, SpatialHystereticTerminalInductionCandidate,
    },
};
use fs_exec::Cx;

// Dispatch only the existing native electrode causality. Material queries,
// thermal residuals, heat deposition and accepted-memory handling stay shared.
trait Electrode<C: SpatialConductorMaterial>: HystereticInductionHistory<C> + Sized {
    type Candidate: HystereticInductionDriver<Model = Self, Step = SpatialHystereticInductionStep>;
    #[allow(clippy::too_many_arguments)]
    fn advance(
        &self,
        endpoint_values: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Candidate, SpatialHystereticInductionError>;
}

// The two native owners deliberately retain different stepping interfaces.
// Only this private adapter shares their identical read/receipt operations.
macro_rules! electrode {
    ($model:ident, $candidate:ident, $advance:ident) => {
        impl<C: SpatialConductorMaterial> HystereticInductionHistory<C> for $model<'_, C> {
            fn time_s(&self) -> f64 { $model::time_s(self) }
            fn problem(&self) -> &SpatialHystereticProblem { $model::problem(self) }
            fn conductors(&self) -> &[Option<C>] { $model::conductors(self) }
        }
        impl<'m, C: SpatialConductorMaterial> ThermalDriver for $candidate<'m, C> {
            type Model = $model<'m, C>;
            type Step = SpatialHystereticInductionStep;
            fn account(&self) -> DriverAccount<'_> { Self::receipt_account(self.step()) }
            fn into_parts(self) -> (Self::Model, Self::Step) { $candidate::into_parts(self) }
        }
        impl<C: SpatialConductorMaterial> HystereticInductionDriver for $candidate<'_, C> {
            fn induction_step(&self) -> &SpatialHystereticInductionStep { self.step() }
            fn receipt_account(step: &Self::Step) -> DriverAccount<'_> {
                <SpatialHystereticInductionCandidate<'_, C> as HystereticInductionDriver>::receipt_account(step)
            }
            fn field_trials(&self) -> usize { 1 }
        }
        impl<'m, C: SpatialConductorMaterial> Electrode<C> for $model<'m, C> {
            type Candidate = $candidate<'m, C>;
            fn advance(
                &self,
                endpoint_values: &[f64],
                endpoint_impressed_current_a_m2: &[[f64; 3]],
                duration_s: f64,
                conductors: Vec<Option<C>>,
                controls: SpatialHystereticInductionControls,
                cancelled: &mut impl FnMut() -> bool,
            ) -> Result<Self::Candidate, SpatialHystereticInductionError> {
                self.$advance(endpoint_values, endpoint_impressed_current_a_m2,
                    duration_s, conductors, controls, cancelled)
            }
        }
    };
}
electrode!(
    SpatialHystereticTerminalInduction,
    SpatialHystereticTerminalInductionCandidate,
    step_voltage_with_conductors_candidate
);
electrode!(
    SpatialHystereticCurrentInduction,
    SpatialHystereticCurrentInductionCandidate,
    step_current_with_conductors_candidate
);

impl SpatialMagneticThermalNetwork<'_> {
    /// Solve physical endpoint electrode voltages with the same finite caloric
    /// inventory, original conductivity sources and accepted magnetic memories.
    ///
    /// The original electrode solve determines all spatial currents and reaction
    /// work. Core eddy and pinning heat enter their actual whole-volume receivers
    /// once. Coefficients of each source-selected magnetic fit remain fixed over
    /// the complete supported thermal interval. Geometry and material are stationary.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_voltage_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialHystereticTerminalInduction<'m, C>,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialHystereticInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalThermalCandidate<
            't,
            SpatialHystereticTerminalInduction<'m, C>,
            SpatialHystereticInductionStep,
        >,
        Error,
    > {
        step_bound(
            self,
            thermal,
            target_s,
            field,
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            controls,
            field_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Prescribe total incoming current on each physical electrode and solve its
    /// voltage together with induction, vector memory and finite temperature.
    /// Current balances separately on each original conducting component.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_current_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialHystereticCurrentInduction<'m, C>,
        endpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialHystereticInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalThermalCandidate<
            't,
            SpatialHystereticCurrentInduction<'m, C>,
            SpatialHystereticInductionStep,
        >,
        Error,
    > {
        step_bound(
            self,
            thermal,
            target_s,
            field,
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            controls,
            field_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn step_bound<'t, C: SpatialConductorMaterial, M: Electrode<C>>(
    layout: &SpatialMagneticThermalNetwork<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &M,
    endpoint_values: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHystereticInductionControls,
    initial_energy_increments_j: &[f64],
    cancelled: impl FnMut() -> bool,
) -> Result<SpatialTerminalThermalCandidate<'t, M, SpatialHystereticInductionStep>, Error> {
    let cancel = RefCell::new(cancelled);
    check(&cancel)?;
    crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let network = MaterialNetwork {
        cells: layout.cells,
        resistors: &[],
        links: layout.links,
    };
    let prepared = network
        .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
        .map_err(material)?;
    step_resolved(
        layout,
        thermal,
        target_s,
        field,
        endpoint_values,
        endpoint_impressed_current_a_m2,
        controls,
        field_controls,
        initial_energy_increments_j,
        |state| {
            prepared
                .resolve(state, || (*cancel.borrow_mut())())
                .map_err(material)
        },
        |sample| sample,
        &cancel,
    )
}

#[allow(clippy::too_many_arguments)]
fn step_resolved<
    't,
    S,
    E: From<MagneticThermalError>,
    C: SpatialConductorMaterial,
    M: Electrode<C>,
>(
    layout: &SpatialMagneticThermalNetwork<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &M,
    endpoint_values: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHystereticInductionControls,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'_>,
    ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
    sample: impl Fn(&S) -> &MaterialNetworkSample,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<
    SpatialTerminalThermalCandidate<'t, M, SpatialHystereticInductionStep, S>,
    CircuitEnthalpyError<E>,
> {
    let (thermal, field, report) = layout.step_hysteretic_driver_candidate(
        thermal,
        target_s,
        field,
        controls,
        initial_energy_increments_j,
        resolve,
        sample,
        |_, _, conductors, _| {
            field
                .advance(
                    endpoint_values,
                    endpoint_impressed_current_a_m2,
                    target_s - thermal.time_s(),
                    conductors,
                    field_controls,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(induction)
                .map_err(native_error)
        },
        |_, _, _, _| Ok(()),
        cancel,
    )?;
    Ok(SpatialTerminalThermalCandidate {
        thermal,
        field,
        report,
    })
}

/// Continue physical voltage electrodes and retained magnetic memory with the
/// original complete contact/radiation/coolant heat assembly and typed sources.
#[allow(clippy::too_many_arguments)]
pub fn voltage_heat_assembly_candidate<
    't,
    'm,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &SpatialHystereticTerminalInduction<'m, C>,
    endpoint_voltage_v: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHystereticInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalThermalCandidate<
        't,
        SpatialHystereticTerminalInduction<'m, C>,
        SpatialHystereticInductionStep,
        A::Sample,
    >,
    HeatError<A::Error>,
> {
    heat_candidate(
        provider,
        cx,
        thermal,
        target_s,
        field,
        endpoint_voltage_v,
        endpoint_impressed_current_a_m2,
        regions,
        volume_tolerance_m3,
        controls,
        field_controls,
        initial_energy_increments_j,
    )
}

/// Continue physical total-current electrodes through the same complete finite
/// heat assembly; solved voltage supplies the actual original terminal work.
#[allow(clippy::too_many_arguments)]
pub fn current_heat_assembly_candidate<
    't,
    'm,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &SpatialHystereticCurrentInduction<'m, C>,
    endpoint_incoming_current_a: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHystereticInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalThermalCandidate<
        't,
        SpatialHystereticCurrentInduction<'m, C>,
        SpatialHystereticInductionStep,
        A::Sample,
    >,
    HeatError<A::Error>,
> {
    heat_candidate(
        provider,
        cx,
        thermal,
        target_s,
        field,
        endpoint_incoming_current_a,
        endpoint_impressed_current_a_m2,
        regions,
        volume_tolerance_m3,
        controls,
        field_controls,
        initial_energy_increments_j,
    )
}

#[allow(clippy::too_many_arguments)]
fn heat_candidate<'t, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial, M: Electrode<C>>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &M,
    endpoint_values: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHystereticInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalThermalCandidate<'t, M, SpatialHystereticInductionStep, A::Sample>,
    HeatError<A::Error>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error)?;
    let layout = SpatialMagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: regions,
        volume_tolerance_m3,
    };
    let mut next = step_resolved(
        &layout,
        thermal,
        target_s,
        field,
        endpoint_values,
        endpoint_impressed_current_a_m2,
        controls,
        field_controls,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|error| error.map_source(SpatialInductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(native_error(invalid(
                    "electrode induction cannot discard graph electrical work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    close_heat_balance(
        &mut next.report,
        target_s - thermal.time_s(),
        controls.total_energy_tolerance_j,
    )?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}
