//! Finite temperature, native 3D induction and finite rotor momentum together.
//!
//! The original finite-enthalpy residual owns temperatures and material
//! resolution. Each temperature trial restarts the same accepted field and
//! mechanical history, then closes motional EMF and rotor speed. Only local
//! conductor-frame Joule heat and the existing reversible magnetic/caloric
//! exchange enter finite caloric storage. Lorentz work is an internal transfer.

use super::{
    LinearInductionThermalDriver, SpatialInductionStep, SpatialInductionThermalReport,
    SpatialMagneticRegion, SpatialMagneticThermalLaw, SpatialMagneticThermalNetwork, check,
    invalid, material, native_error,
};
use crate::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve, CircuitThermalCoefficients, DriverAccount,
        ThermalDriver,
        contacts::heat::{HeatAssembly, HeatNetworkError, HeatNetworkSample},
        magnetic::MagneticThermalError,
        material::{MaterialNetwork, MaterialNetworkSample},
    },
    induction_rotor::{RotorBalance, RotorControls},
    spatial_induction::{SpatialInductionRotor, SpatialRotorCandidate, SpatialRotorError},
};
use fs_em::spatial::vector_potential::{
    VectorMagneticBoundary, induction::SpatialInductionControls,
};
use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_math::dd::Dd;
use fs_mbd::articulated::PrescribedBaseImpulseResponse;
use fs_thermal::enthalpy::DepositionSimulation;
use std::cell::RefCell;

#[derive(Debug)]
pub enum SpatialRotorThermalError<E = HeatNetworkError> {
    Native(MagneticThermalError),
    Rotor(SpatialRotorError),
    Heat(E),
}
impl<E> From<MagneticThermalError> for SpatialRotorThermalError<E> {
    fn from(error: MagneticThermalError) -> Self {
        Self::Native(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpatialRotorThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial induction rotor/thermal coupling: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpatialRotorThermalError<E> {}
pub type SpatialRotorEnthalpyError<E = HeatNetworkError> =
    CircuitEnthalpyError<SpatialRotorThermalError<E>>;

/// Actual accepted field and momentum accounts at the converged temperature.
pub struct SpatialRotorThermalStep {
    pub induction: SpatialInductionStep,
    pub mechanics: PrescribedBaseImpulseResponse,
    pub balance: RotorBalance,
}

/// The physical residual includes rotor kinetic and field internal energy,
/// finite caloric storage, external shaft/impressed work and external heat.
pub type SpatialRotorThermalReport<S = MaterialNetworkSample> =
    SpatialInductionThermalReport<SpatialRotorThermalStep, S>;

pub struct SpatialRotorThermalCandidate<'t, 'm, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    rotor: SpatialInductionRotor<'m>,
    report: SpatialRotorThermalReport<S>,
}
impl<'t, 'm, S> SpatialRotorThermalCandidate<'t, 'm, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn rotor(&self) -> &SpatialInductionRotor<'m> {
        &self.rotor
    }
    pub fn report(&self) -> &SpatialRotorThermalReport<S> {
        &self.report
    }
    /// Publish all three histories together. Reported Joule/reversible heat
    /// has already entered the returned caloric state.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialInductionRotor<'m>,
        SpatialRotorThermalReport<S>,
    ) {
        (self.thermal, self.rotor, self.report)
    }
}

struct RotorDriver<'m>(SpatialRotorCandidate<'m>);
impl<'m> ThermalDriver for RotorDriver<'m> {
    type Model = SpatialInductionRotor<'m>;
    type Step = SpatialRotorThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        let field = &self.0.induction;
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: field.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: (Dd::from_f64(field.stored_after_j)
                - Dd::from_f64(field.stored_before_j)
                - Dd::from_f64(field.parameter_work_j)
                + Dd::from_f64(self.0.mechanics.kinetic_energy_change))
            .to_f64(),
            supplied_j: (Dd::from_f64(field.impressed_work_j)
                + Dd::from_f64(field.magnetic_boundary_work_j)
                + Dd::from_f64(self.0.balance.shaft_work_j))
            .to_f64(),
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        let (rotor, induction, mechanics, balance) = self.0.into_parts();
        (
            rotor,
            SpatialRotorThermalStep {
                induction,
                mechanics,
                balance,
            },
        )
    }
}
impl LinearInductionThermalDriver for RotorDriver<'_> {
    fn induction_step(&self) -> &SpatialInductionStep {
        &self.0.induction
    }
    fn induction_receipt(step: &Self::Step) -> &SpatialInductionStep {
        &step.induction
    }
    fn field_trials(&self) -> usize {
        self.0.balance.iterations
    }
    fn additional_stored_change_j(step: &Self::Step) -> f64 {
        step.mechanics.kinetic_energy_change
    }
    fn additional_supplied_j(step: &Self::Step) -> f64 {
        step.balance.shaft_work_j
    }
    fn internal_mechanical_work_j(step: &Self::Step) -> f64 {
        step.induction.mechanical_work_j
    }
}

fn rotor_error<E>(error: SpatialRotorError) -> SpatialRotorEnthalpyError<E> {
    match error {
        SpatialRotorError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialRotorThermalError::Rotor(other)),
    }
}
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> SpatialRotorEnthalpyError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialRotorThermalError::Heat(other.into())),
    }
}

impl SpatialMagneticThermalNetwork<'_> {
    /// Solve temperature-dependent induction, Joule heating and finite rotor
    /// inertia in the original finite-enthalpy residual. Rotor conductor
    /// temperature must remain homogeneous; mapping its complete volume to
    /// one finite caloric cell is the usual realization. Additional thermal
    /// bodies, links and reservoirs retain their original material ownership.
    ///
    /// Magnetic entropy/parameter work use the original source-bound law.
    /// Rotating material remains scalar, homogeneous and fully solid, with
    /// uniform permeability throughout the field domain. All field solves
    /// inside speed closure consume the shared thermal trial budget. Endpoint
    /// temperature coefficients make this coupled method generally first order.
    #[allow(clippy::too_many_arguments)]
    pub fn step_rotor_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        rotor: &SpatialInductionRotor<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        shaft_torque_nm: f64,
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        rotor_controls: RotorControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialRotorThermalCandidate<'t, 'm>, SpatialRotorEnthalpyError> {
        self.step_rotor_boundary_candidate(
            thermal,
            target_s,
            rotor,
            endpoint_current_a_m2,
            None,
            shaft_torque_nm,
            bindings,
            controls,
            field_controls,
            rotor_controls,
            assembly,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Apply an affine magnetic boundary history to the same jointly solved
    /// field, finite rotor momentum and finite enthalpy. Each thermal/speed
    /// trial reuses the accepted initial trace and the supplied endpoint.
    /// Magnetic boundary reaction work enters external supply once; Lorentz
    /// work remains an internal field/rotor transfer, and only original Joule
    /// and reversible material heat enter caloric storage.
    #[allow(clippy::too_many_arguments)]
    pub fn step_rotor_with_boundary_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        rotor: &SpatialInductionRotor<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        endpoint_boundary: &VectorMagneticBoundary,
        shaft_torque_nm: f64,
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        rotor_controls: RotorControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialRotorThermalCandidate<'t, 'm>, SpatialRotorEnthalpyError> {
        self.step_rotor_boundary_candidate(
            thermal,
            target_s,
            rotor,
            endpoint_current_a_m2,
            Some(endpoint_boundary),
            shaft_torque_nm,
            bindings,
            controls,
            field_controls,
            rotor_controls,
            assembly,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_rotor_boundary_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        rotor: &SpatialInductionRotor<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        endpoint_boundary: Option<&VectorMagneticBoundary>,
        shaft_torque_nm: f64,
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        rotor_controls: RotorControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialRotorThermalCandidate<'t, 'm>, SpatialRotorEnthalpyError> {
        let cancel = RefCell::new(cancelled);
        check(&cancel).map_err(native_error)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)
            .map_err(native_error)?;
        self.step_rotor_resolved_candidate(
            thermal,
            target_s,
            rotor,
            endpoint_current_a_m2,
            endpoint_boundary,
            shaft_torque_nm,
            bindings,
            controls,
            field_controls,
            rotor_controls,
            assembly,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
                    .map_err(native_error)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_rotor_resolved_candidate<'t, 'm, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        rotor: &SpatialInductionRotor<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        endpoint_boundary: Option<&VectorMagneticBoundary>,
        shaft_torque_nm: f64,
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        rotor_controls: RotorControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, SpatialRotorEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialRotorThermalCandidate<'t, 'm, S>, SpatialRotorEnthalpyError<E>> {
        crate::induction_rotor::admit_controls(
            target_s - thermal.time_s(),
            shaft_torque_nm,
            rotor_controls,
        )
        .map_err(SpatialRotorError::from)
        .map_err(rotor_error)?;
        let (thermal, rotor, report) = self.step_linear_driver_candidate(
            thermal,
            target_s,
            rotor.field(),
            bindings,
            controls,
            initial_energy_increments_j,
            resolve,
            sample,
            |_, _, materials, conductors, remaining| {
                if remaining == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let mut bounded = rotor_controls;
                bounded.maximum_iterations = bounded.maximum_iterations.min(remaining);
                match endpoint_boundary {
                    Some(boundary) => rotor.step_with_materials_and_boundary_candidate(
                        target_s,
                        endpoint_current_a_m2,
                        boundary,
                        shaft_torque_nm,
                        materials,
                        conductors,
                        assembly,
                        field_controls,
                        bounded,
                        &mut || (*cancel.borrow_mut())(),
                    ),
                    None => rotor.step_with_materials_candidate(
                        target_s,
                        endpoint_current_a_m2,
                        shaft_torque_nm,
                        materials,
                        conductors,
                        assembly,
                        field_controls,
                        bounded,
                        &mut || (*cancel.borrow_mut())(),
                    ),
                }
                .map(RotorDriver)
                .map_err(|error| match error {
                    SpatialRotorError::Nonconvergence { iterations, .. }
                        if iterations == remaining
                            && remaining < rotor_controls.maximum_iterations =>
                    {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    other => rotor_error(other),
                })
            },
            |_, _, _, _| Ok(()),
            cancel,
        )?;
        check(cancel).map_err(native_error)?;
        Ok(SpatialRotorThermalCandidate {
            thermal,
            rotor,
            report,
        })
    }
}

/// Apply the same rotor/temperature residual to a complete heat assembly.
/// Contact, layer, radiation and finite coolant payloads and typed failures
/// remain those of the provider. Electrical resistor/contact work and field-
/// dependent thermal boundaries refuse without an owner in this formulation.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    rotor: &SpatialInductionRotor<'m>,
    endpoint_current_a_m2: &[[f64; 3]],
    shaft_torque_nm: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    bindings: &[SpatialMagneticThermalLaw],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialInductionControls,
    rotor_controls: RotorControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialRotorThermalCandidate<'t, 'm, A::Sample>, SpatialRotorEnthalpyError<A::Error>> {
    heat_assembly_boundary_candidate(
        provider,
        cx,
        thermal,
        target_s,
        rotor,
        endpoint_current_a_m2,
        None,
        shaft_torque_nm,
        regions,
        volume_tolerance_m3,
        bindings,
        controls,
        field_controls,
        rotor_controls,
        assembly,
        initial_energy_increments_j,
    )
}

/// Drive the original rotor with an external magnetic boundary history while
/// the same complete heat assembly resolves all finite thermal transport.
/// The original sources, geometry, complete trace support and jointly accepted
/// field/rotor/caloric histories remain owned by the existing coupled residual.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_with_boundary_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    rotor: &SpatialInductionRotor<'m>,
    endpoint_current_a_m2: &[[f64; 3]],
    endpoint_boundary: &VectorMagneticBoundary,
    shaft_torque_nm: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    bindings: &[SpatialMagneticThermalLaw],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialInductionControls,
    rotor_controls: RotorControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialRotorThermalCandidate<'t, 'm, A::Sample>, SpatialRotorEnthalpyError<A::Error>> {
    heat_assembly_boundary_candidate(
        provider,
        cx,
        thermal,
        target_s,
        rotor,
        endpoint_current_a_m2,
        Some(endpoint_boundary),
        shaft_torque_nm,
        regions,
        volume_tolerance_m3,
        bindings,
        controls,
        field_controls,
        rotor_controls,
        assembly,
        initial_energy_increments_j,
    )
}

#[allow(clippy::too_many_arguments)]
fn heat_assembly_boundary_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    rotor: &SpatialInductionRotor<'m>,
    endpoint_current_a_m2: &[[f64; 3]],
    endpoint_boundary: Option<&VectorMagneticBoundary>,
    shaft_torque_nm: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    bindings: &[SpatialMagneticThermalLaw],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialInductionControls,
    rotor_controls: RotorControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialRotorThermalCandidate<'t, 'm, A::Sample>, SpatialRotorEnthalpyError<A::Error>> {
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
    let mut next = layout.step_rotor_resolved_candidate(
        thermal,
        target_s,
        rotor,
        endpoint_current_a_m2,
        endpoint_boundary,
        shaft_torque_nm,
        bindings,
        controls,
        field_controls,
        rotor_controls,
        assembly,
        initial_energy_increments_j,
        |trial| {
            let coefficients = provider
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(SpatialRotorThermalError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(native_error(invalid(
                    "rotor heat cannot discard unowned electrical work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    close_heat_balance::<_, A::Sample, A::Error>(
        &mut next.report,
        target_s - thermal.time_s(),
        controls.total_energy_tolerance_j,
    )?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}

fn close_heat_balance<D, S: AsRef<HeatNetworkSample>, E: From<HeatNetworkError>>(
    report: &mut SpatialInductionThermalReport<D, S>,
    duration_s: f64,
    tolerance_j: f64,
) -> Result<(), SpatialRotorEnthalpyError<E>> {
    super::heat::close_heat_balance::<D, S, E>(report, duration_s, tolerance_j).map_err(|error| {
        error.map_source(|source| match source {
            super::heat::SpatialInductionHeatError::Native(error) => {
                SpatialRotorThermalError::Native(error)
            }
            super::heat::SpatialInductionHeatError::Heat(error) => {
                SpatialRotorThermalError::Heat(error)
            }
        })
    })
}
