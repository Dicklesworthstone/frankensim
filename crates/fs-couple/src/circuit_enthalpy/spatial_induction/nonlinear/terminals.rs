//! Physical electrode work, nonlinear magnetic storage and finite heat.
//!
//! Voltage and total-current electrodes retain their original charge topology
//! through every thermal trial. The existing constitutive/enthalpy driver owns
//! source refresh, optional magnetic entropy, local Joule deposition and atomic
//! acceptance. A stationary body exchanges electrical work with its electrodes;
//! numerical magnetic dissipation remains separate from deposited heat.

pub use super::super::terminals::SpatialTerminalThermalCandidate;
use super::super::{
    SpatialMagneticRegion,
    heat::{HeatError, SpatialInductionHeatError, close_heat_balance, heat_error},
};
use super::*;
use crate::circuit_enthalpy::contacts::heat::HeatAssembly;
use fs_em::spatial::vector_potential::induction::nonlinear::terminals::{
    SpatialNonlinearCurrentInduction, SpatialNonlinearCurrentInductionCandidate,
    SpatialNonlinearTerminalInduction, SpatialNonlinearTerminalInductionCandidate,
};
use fs_exec::Cx;

impl<C: SpatialConductorMaterial> NonlinearInductionHistory<C>
    for SpatialNonlinearTerminalInduction<'_, C>
{
    fn time_s(&self) -> f64 {
        SpatialNonlinearTerminalInduction::time_s(self)
    }
    fn problem(&self) -> &SpatialNonlinearVectorMagnetostatics {
        SpatialNonlinearTerminalInduction::problem(self)
    }
    fn conductors(&self) -> &[Option<C>] {
        SpatialNonlinearTerminalInduction::conductors(self)
    }
    fn field(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionField, SpatialNonlinearInductionError> {
        SpatialNonlinearTerminalInduction::field(self, cancelled)
    }
}

impl<C: SpatialConductorMaterial> NonlinearInductionHistory<C>
    for SpatialNonlinearCurrentInduction<'_, C>
{
    fn time_s(&self) -> f64 {
        SpatialNonlinearCurrentInduction::time_s(self)
    }
    fn problem(&self) -> &SpatialNonlinearVectorMagnetostatics {
        SpatialNonlinearCurrentInduction::problem(self)
    }
    fn conductors(&self) -> &[Option<C>] {
        SpatialNonlinearCurrentInduction::conductors(self)
    }
    fn field(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionField, SpatialNonlinearInductionError> {
        SpatialNonlinearCurrentInduction::field(self, cancelled)
    }
}

// Only electrode causality differs. Source preparation, stationary motion and
// every thermal/caloric account remain in the original shared driver.
trait NonlinearElectrode<C: SpatialConductorMaterial>: NonlinearInductionHistory<C> + Sized {
    type Candidate: InductionThermalDriver<Model = Self, Step = SpatialNonlinearInductionStep>;

    #[allow(clippy::too_many_arguments)]
    fn drive(
        &self,
        endpoint_values: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Option<Vec<SpatialNonlinearMaterial>>,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Candidate, SpatialNonlinearInductionError>;
}

impl<'m, C: SpatialConductorMaterial> NonlinearElectrode<C>
    for SpatialNonlinearTerminalInduction<'m, C>
{
    type Candidate = SpatialNonlinearTerminalInductionCandidate<'m, C>;

    fn drive(
        &self,
        values: &[f64],
        impressed: &[[f64; 3]],
        velocity: &[[f64; 3]],
        dt: f64,
        materials: Option<Vec<SpatialNonlinearMaterial>>,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Candidate, SpatialNonlinearInductionError> {
        match materials {
            Some(materials) => self.step_voltage_with_materials_candidate(
                values, impressed, velocity, dt, materials, conductors, controls, cancelled,
            ),
            None => self.step_voltage_with_conductors_candidate(
                values, impressed, velocity, dt, conductors, controls, cancelled,
            ),
        }
    }
}

impl<'m, C: SpatialConductorMaterial> NonlinearElectrode<C>
    for SpatialNonlinearCurrentInduction<'m, C>
{
    type Candidate = SpatialNonlinearCurrentInductionCandidate<'m, C>;

    fn drive(
        &self,
        values: &[f64],
        impressed: &[[f64; 3]],
        velocity: &[[f64; 3]],
        dt: f64,
        materials: Option<Vec<SpatialNonlinearMaterial>>,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Candidate, SpatialNonlinearInductionError> {
        match materials {
            Some(materials) => self.step_current_with_materials_candidate(
                values, impressed, velocity, dt, materials, conductors, controls, cancelled,
            ),
            None => self.step_current_with_conductors_candidate(
                values, impressed, velocity, dt, conductors, controls, cancelled,
            ),
        }
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver
    for SpatialNonlinearTerminalInductionCandidate<'m, C>
{
    type Model = SpatialNonlinearTerminalInduction<'m, C>;
    type Step = SpatialNonlinearInductionStep;

    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialNonlinearTerminalInductionCandidate::into_parts(self)
    }
}

impl<C: SpatialConductorMaterial> InductionThermalDriver
    for SpatialNonlinearTerminalInductionCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialNonlinearInductionStep {
        self.step()
    }
    fn induction_receipt(step: &Self::Step) -> &SpatialNonlinearInductionStep {
        step
    }
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_> {
        <SpatialNonlinearInductionCandidate<'_, C> as InductionThermalDriver>::receipt_account(step)
    }
    fn field_trials(&self) -> usize {
        1
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver
    for SpatialNonlinearCurrentInductionCandidate<'m, C>
{
    type Model = SpatialNonlinearCurrentInduction<'m, C>;
    type Step = SpatialNonlinearInductionStep;

    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialNonlinearCurrentInductionCandidate::into_parts(self)
    }
}

impl<C: SpatialConductorMaterial> InductionThermalDriver
    for SpatialNonlinearCurrentInductionCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialNonlinearInductionStep {
        self.step()
    }
    fn induction_receipt(step: &Self::Step) -> &SpatialNonlinearInductionStep {
        step
    }
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_> {
        <SpatialNonlinearInductionCandidate<'_, C> as InductionThermalDriver>::receipt_account(step)
    }
    fn field_trials(&self) -> usize {
        1
    }
}

impl SpatialMagneticThermalNetwork<'_> {
    /// Couple physical electrode voltages to stationary nonlinear induction and
    /// finite enthalpy. Local conductor losses enter their original receivers;
    /// terminal V*I work and impressed work fund storage, heat and numerical loss.
    ///
    /// Conductivity follows each actual endpoint temperature. Explicit magnetic
    /// laws enable the original source-separable entropy exchange; without them,
    /// magnetic coefficients must remain fixed over the admitted temperature
    /// interval. Every trial restarts from the accepted field and caloric state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_voltage_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearTerminalInduction<'m, C>,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalThermalCandidate<
            't,
            SpatialNonlinearTerminalInduction<'m, C>,
            SpatialNonlinearInductionStep,
        >,
        Error,
    > {
        let cancel = RefCell::new(cancelled);
        self.step_nonlinear_terminal_candidate(
            thermal,
            target_s,
            field,
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            magnetic_laws,
            controls,
            field_controls,
            initial_energy_increments_j,
            &cancel,
        )
    }

    /// Couple signed incoming total electrode currents to stationary nonlinear
    /// induction and finite heat. The native grouped charge rows determine the
    /// electrode voltages and retain the current-driven reference convention.
    ///
    /// Original conductor/material identity, whole-volume solid receivers and
    /// optional magnetic entropy have the same admission as voltage drive.
    /// Cancellation or a failed source, solve or balance publishes neither state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_current_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearCurrentInduction<'m, C>,
        endpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalThermalCandidate<
            't,
            SpatialNonlinearCurrentInduction<'m, C>,
            SpatialNonlinearInductionStep,
        >,
        Error,
    > {
        let cancel = RefCell::new(cancelled);
        self.step_nonlinear_terminal_candidate(
            thermal,
            target_s,
            field,
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            magnetic_laws,
            controls,
            field_controls,
            initial_energy_increments_j,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_nonlinear_terminal_candidate<
        't,
        C: SpatialConductorMaterial,
        O: NonlinearElectrode<C>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &O,
        endpoint_values: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialTerminalThermalCandidate<'t, O, SpatialNonlinearInductionStep>, Error> {
        check(cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_nonlinear_terminal_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_values,
            endpoint_impressed_current_a_m2,
            magnetic_laws,
            controls,
            field_controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_nonlinear_terminal_resolved_candidate<
        't,
        C: SpatialConductorMaterial,
        O: NonlinearElectrode<C>,
        S,
        E: From<MagneticThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &O,
        endpoint_values: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        snapshot: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        SpatialTerminalThermalCandidate<'t, O, SpatialNonlinearInductionStep, S>,
        CircuitEnthalpyError<E>,
    > {
        check(cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        // Thermal cells are stationary; magnetic work must not silently stand
        // in for absent heat advection from a prescribed moving body.
        let velocity = vec![[0.0; 3]; field.problem().positions_m().len()];
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = self.step_constitutive_driver_candidate(
            thermal,
            target_s,
            field,
            magnetic_laws,
            controls,
            initial_energy_increments_j,
            resolve,
            snapshot,
            |_, _, materials, conductors, _| {
                field
                    .drive(
                        endpoint_values,
                        endpoint_impressed_current_a_m2,
                        &velocity,
                        dt,
                        materials,
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
}

/// Continue voltage-driven nonlinear induction with the complete original heat
/// assembly. Thermal contacts, radiation and finite coolant share the same
/// endpoint temperatures, source payload and optional magnetic entropy exchange.
/// Electrical graph work requires its own owner and refuses in this field path.
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
    field: &SpatialNonlinearTerminalInduction<'m, C>,
    endpoint_voltage_v: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialNonlinearInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalThermalCandidate<
        't,
        SpatialNonlinearTerminalInduction<'m, C>,
        SpatialNonlinearInductionStep,
        A::Sample,
    >,
    HeatError<A::Error>,
> {
    terminal_heat_assembly_candidate(
        provider,
        cx,
        thermal,
        target_s,
        field,
        endpoint_voltage_v,
        endpoint_impressed_current_a_m2,
        regions,
        volume_tolerance_m3,
        magnetic_laws,
        controls,
        field_controls,
        initial_energy_increments_j,
    )
}

/// Continue total-current-driven nonlinear induction with the complete original
/// heat assembly. Native electrode voltages supply the signed terminal work;
/// only physical local heat is deposited. Closed radiation and coolant exchange
/// retain their original internal accounting and finite histories.
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
    field: &SpatialNonlinearCurrentInduction<'m, C>,
    endpoint_incoming_current_a: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialNonlinearInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalThermalCandidate<
        't,
        SpatialNonlinearCurrentInduction<'m, C>,
        SpatialNonlinearInductionStep,
        A::Sample,
    >,
    HeatError<A::Error>,
> {
    terminal_heat_assembly_candidate(
        provider,
        cx,
        thermal,
        target_s,
        field,
        endpoint_incoming_current_a,
        endpoint_impressed_current_a_m2,
        regions,
        volume_tolerance_m3,
        magnetic_laws,
        controls,
        field_controls,
        initial_energy_increments_j,
    )
}

#[allow(clippy::too_many_arguments)]
fn terminal_heat_assembly_candidate<
    't,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
    O: NonlinearElectrode<C>,
>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &O,
    endpoint_values: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialNonlinearInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalThermalCandidate<'t, O, SpatialNonlinearInductionStep, A::Sample>,
    HeatError<A::Error>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
    let mut magnetic_axes = Vec::new();
    for region in regions {
        check(&cancel).map_err(native_error)?;
        if let SpatialMagneticRegion::Cell { tetrahedron, cell } = *region {
            if let Some(SpatialNonlinearMaterial::Reversible(law)) =
                field.problem().materials().get(tetrahedron)
            {
                magnetic_axes.push((cell, law.source_path().axis.as_str()));
            }
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &magnetic_axes)
        .map_err(heat_error)?;
    let layout = SpatialMagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: regions,
        volume_tolerance_m3,
    };
    let mut next = layout.step_nonlinear_terminal_resolved_candidate(
        thermal,
        target_s,
        field,
        endpoint_values,
        endpoint_impressed_current_a_m2,
        magnetic_laws,
        controls,
        field_controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = provider
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(SpatialInductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(native_error(invalid(
                    "spatial electrode heat cannot discard graph electrical work",
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
