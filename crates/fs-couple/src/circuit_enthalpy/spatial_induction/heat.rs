//! Complete material heat transport on the retained three-dimensional fields.
//! The existing interface/layer/radiosity and coolant owners evaluate the same
//! endpoint temperatures as conductivity and magnetic storage; no second heat
//! path is stepped and no closed-enclosure transfer becomes external supply.

use super::nonlinear::SpatialNonlinearInductionThermalCandidate;
use super::*;
use crate::circuit_enthalpy::{
    contacts::heat::{
        HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance,
    },
    magnetic::{MagneticThermalError, caloric::MagneticConstitutiveThermalLaw},
};
use fs_em::spatial::{
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::induction::nonlinear::{
        SpatialNonlinearInduction, SpatialNonlinearInductionControls,
    },
};
use fs_exec::Cx;

/// Original field/material refusal or the complete heat provider's typed error.
/// Native cancellation, work exhaustion and thermal/numerical failures remain
/// in their original outer `CircuitEnthalpyError` variants.
#[derive(Debug)]
pub enum SpatialInductionHeatError<E = HeatNetworkError> {
    Native(MagneticThermalError),
    Heat(E),
}
impl<E> From<MagneticThermalError> for SpatialInductionHeatError<E> {
    fn from(error: MagneticThermalError) -> Self {
        Self::Native(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpatialInductionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial induction/material heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpatialInductionHeatError<E> {}
type HeatError<E> = CircuitEnthalpyError<SpatialInductionHeatError<E>>;

fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> HeatError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialInductionHeatError::Heat(other.into())),
    }
}

// Existing concrete entry points retain their exact public error type while
// delegating to the same complete assembly path.
fn concrete_error(error: SpatialInductionHeatError<HeatNetworkError>) -> MagneticThermalError {
    match error {
        SpatialInductionHeatError::Native(error) => error,
        SpatialInductionHeatError::Heat(error) => MagneticThermalError::Heat(error),
    }
}

fn close_heat_balance<D, S: AsRef<HeatNetworkSample>, E: From<HeatNetworkError>>(
    report: &mut SpatialInductionThermalReport<D, S>,
    dt: f64,
    tolerance: f64,
) -> Result<(), HeatError<E>> {
    let sample = report.coupled.source.network.as_ref();
    // Both the original native/free-energy account and the physical magnetic
    // internal-energy account must exclude the same internal radiation net.
    closed_balance(
        dt,
        report.coupled.total_energy_residual_j,
        sample,
        tolerance,
    )
    .map_err(heat_error)?;
    closed_balance(dt, report.total_energy_residual_j, sample, tolerance).map_err(heat_error)?;
    let internal_net = finite(dt * sample.radiation_net_power_w())?;
    report.external_heat_j = finite(report.external_heat_j - internal_net)?;
    report.total_energy_residual_j = finite(report.total_energy_residual_j + internal_net)?;
    Ok(())
}

/// Continue the original three-dimensional A–phi and finite-enthalpy histories
/// with any complete material heat assembly, including finite material coolant.
///
/// Each trial resolves the SAME proposed thermal state, then the native field
/// restarts from the same accepted magnetic history. The original conductor
/// owner supplies sigma(T) and local conductor-frame Joule heat; the original
/// magnetic caloric law supplies affine-mu(T) entropy, internal energy and
/// reversible heat. Finite coolant exchange stays internal. Prescribed material
/// reservoirs, thermal contacts, integrated layers and radiation retain their
/// full accepted source samples in `report().coupled.source.network`.
///
/// `regions` and `bindings` retain the original whole-volume, solid-chart
/// and magnetic-capacity admissions. Passive caloric cells need not belong to
/// the magnetic mesh. The assembly must own all native thermal links; no
/// initial cooling coefficient is frozen for later steps. Graph electrical
/// bindings and resolved resistor/contact work refuse because they have no
/// electrical owner in this field formulation. Thermal boundaries cannot
/// depend on unowned E/J/B/H.
///
/// Geometry, conducting occupancy, contact areas and enclosure view factors
/// stay fixed. Prescribed velocity exports signed Lorentz work; it does not
/// evolve mechanics. The original first-order thermal approximation and
/// linear-magnetic source domain remain. Both native and physical energy gates
/// exclude closed-radiation roundoff from external supply. A refusal publishes
/// neither thermal nor field history and retains the provider's typed error.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &SpatialInduction<'m, C>,
    endpoint_current_a_m2: &[[f64; 3]],
    velocity_m_s: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    bindings: &[SpatialMagneticThermalLaw],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialInductionControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialInductionThermalCandidate<'t, 'm, A::Sample, C>, HeatError<A::Error>> {
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
    let mut next = layout.step_resolved_candidate(
        thermal,
        target_s,
        field,
        endpoint_current_a_m2,
        velocity_m_s,
        bindings,
        controls,
        field_controls,
        assembly,
        initial_energy_increments_j,
        |trial| {
            let coefficients = provider
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(SpatialInductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(native_error(invalid(
                    "spatial induction cannot discard graph electrical work",
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

/// Continue native saturating/fixed-recoil induction with the SAME complete
/// material heat assembly and original finite thermal inventory.
///
/// Conductivity follows each actual caloric temperature while the complete
/// original magnetic law remains fixed over its admitted field range. Only
/// physical local Joule heat enters the finite stores; backward-Euler magnetic
/// loss remains a separate numerical export. Coolant, source selection,
/// field-independent boundaries, complete link ownership, typed refusals and
/// closed-radiation accounting follow `heat_assembly_candidate`.
///
/// No restarted field, fitted core loss, separate coolant temperature history
/// or duplicated heat is introduced. Existing fixed-volume fully-solid field
/// receivers, source domains and first-order thermal coupling remain. A
/// temperature-dependent nonlinear magnetic law still requires its own
/// field-entropy receiver and refuses here.
#[allow(clippy::too_many_arguments)]
pub fn nonlinear_heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &SpatialNonlinearInduction<'m, C>,
    endpoint_current_a_m2: &[[f64; 3]],
    velocity_m_s: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialNonlinearInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, A::Sample, C>, HeatError<A::Error>> {
    nonlinear_heat_assembly_with_laws_candidate(
        provider,
        cx,
        thermal,
        target_s,
        field,
        endpoint_current_a_m2,
        velocity_m_s,
        regions,
        volume_tolerance_m3,
        None,
        controls,
        field_controls,
        initial_energy_increments_j,
    )
}

/// Continue sourced nonlinear magnetic temperature laws with the SAME complete
/// heat assembly, including finite material coolant.
///
/// Explicit original linear mu(T) or separable H(B,T)=s(T)H_ref(B) bindings
/// supply the field-entropy receiver. Conductivity, magnetic source response,
/// boundary heat and coolant resolve at each proposed thermal endpoint. The
/// original field retains A and current history; its parameter work and the
/// original caloric owner's reversible heat close U=psi+T*S without double
/// deposition. Numerical magnetic loss remains a separate export.
///
/// The source, field-envelope, combined-capacity and whole-volume admissions
/// of `SpatialMagneticThermalNetwork::step_nonlinear_thermal_candidate` remain
/// authoritative. Fixed geometry, sourced scalar or tensor conductivity and an affine positive
/// nonlinear field scale are retained. General H(B,T), recoil entropy,
/// hysteresis and phase change require additional physical owners. Full
/// provider samples and typed failures follow `heat_assembly_candidate`.
#[allow(clippy::too_many_arguments)]
pub fn nonlinear_thermal_heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &SpatialNonlinearInduction<'m, C>,
    endpoint_current_a_m2: &[[f64; 3]],
    velocity_m_s: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    bindings: &[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>],
    controls: CircuitEnthalpySolve,
    field_controls: SpatialNonlinearInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, A::Sample, C>, HeatError<A::Error>> {
    nonlinear_heat_assembly_with_laws_candidate(
        provider,
        cx,
        thermal,
        target_s,
        field,
        endpoint_current_a_m2,
        velocity_m_s,
        regions,
        volume_tolerance_m3,
        Some(bindings),
        controls,
        field_controls,
        initial_energy_increments_j,
    )
}

// Fixed and temperature-dependent magnetic material paths share complete heat
// ownership, source admission, coolant resolution and physical energy checks.
#[allow(clippy::too_many_arguments)]
fn nonlinear_heat_assembly_with_laws_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    field: &SpatialNonlinearInduction<'m, C>,
    endpoint_current_a_m2: &[[f64; 3]],
    velocity_m_s: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    bindings: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialNonlinearInductionControls,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, A::Sample, C>, HeatError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
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
    let mut next = layout.step_constitutive_resolved_candidate(
        thermal,
        target_s,
        field,
        endpoint_current_a_m2,
        velocity_m_s,
        bindings,
        controls,
        field_controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = provider
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(SpatialInductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(native_error(invalid(
                    "spatial induction cannot discard graph electrical work",
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

impl PreparedHeatNetwork<'_, '_> {
    /// Continue three-dimensional induction with this prepared heat network
    /// through the same complete `heat_assembly_candidate` driver.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_induction_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialInductionThermalCandidate<'t, 'm, HeatNetworkSample, C>, Error> {
        heat_assembly_candidate(
            self,
            cx,
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            regions,
            volume_tolerance_m3,
            bindings,
            controls,
            field_controls,
            assembly,
            initial_energy_increments_j,
        )
        .map_err(|error| error.map_source(concrete_error))
    }

    /// Continue original sourced nonlinear magnetic temperature response
    /// through the same complete `nonlinear_thermal_heat_assembly_candidate`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_nonlinear_induction_thermal_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, HeatNetworkSample, C>, Error> {
        nonlinear_thermal_heat_assembly_candidate(
            self,
            cx,
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            regions,
            volume_tolerance_m3,
            bindings,
            controls,
            field_controls,
            initial_energy_increments_j,
        )
        .map_err(|error| error.map_source(concrete_error))
    }

    /// Continue native saturating/fixed-recoil induction through the same
    /// complete `nonlinear_heat_assembly_candidate` driver.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_nonlinear_induction_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, HeatNetworkSample, C>, Error> {
        nonlinear_heat_assembly_candidate(
            self,
            cx,
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            regions,
            volume_tolerance_m3,
            controls,
            field_controls,
            initial_energy_increments_j,
        )
        .map_err(|error| error.map_source(concrete_error))
    }
}
