//! Complete shared heat transport on the original planar induction histories.
//! All providers resolve at the same thermal endpoint as conductivity and
//! magnetic response. Native field, caloric storage and heat-transfer owners
//! retain their states, local work accounts and source refusals.

use super::*;
use crate::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::heat::{
        HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance,
    },
    eddy::{
        EddyEnthalpyCandidate,
        nonlinear::{NonlinearEddyEnthalpyCandidate, thermal::NonlinearEddyCaloricCandidate},
    },
    magnetic::{MagneticRegion, MagneticThermalError, caloric::MagneticConstitutiveThermalLaw},
};
use fs_em::planar::{
    eddy::nonlinear::{NonlinearEddyControls, PlanarNonlinearEddy},
    nonlinear::IsotropicMagneticMaterial,
};
use fs_exec::Cx;

/// Native magnetic/material failure or the complete heat provider's error.
/// Cancellation, budgets and numerical failures retain their original outer
/// `CircuitEnthalpyError` classification.
#[derive(Debug)]
pub enum PlanarInductionHeatError<E = HeatNetworkError> {
    Native(MagneticThermalError),
    Heat(E),
}
impl<E> From<MagneticThermalError> for PlanarInductionHeatError<E> {
    fn from(error: MagneticThermalError) -> Self { Self::Native(error) }
}
impl<E: core::fmt::Debug> core::fmt::Display for PlanarInductionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar induction/material heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for PlanarInductionHeatError<E> {}
type HeatError<E> = CircuitEnthalpyError<PlanarInductionHeatError<E>>;

fn assembly_heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> HeatError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(PlanarInductionHeatError::Heat(other.into())),
    }
}
fn concrete_error(error: PlanarInductionHeatError<HeatNetworkError>) -> MagneticThermalError {
    match error {
        PlanarInductionHeatError::Native(error) => error,
        PlanarInductionHeatError::Heat(error) => MagneticThermalError::Heat(error),
    }
}

fn close_assembly_heat_balance<D, S: AsRef<HeatNetworkSample>, E: From<HeatNetworkError>>(
    report: &mut EddyCaloricReport<S, D>, dt: f64, tolerance: f64,
) -> Result<(), HeatError<E>> {
    let sample = report.coupled.source.material.network.as_ref();
    closed_balance(dt, report.coupled.total_energy_residual_j, sample, tolerance)
        .map_err(assembly_heat_error)?;
    closed_balance(dt, report.total_energy_residual_j, sample, tolerance)
        .map_err(assembly_heat_error)?;
    let internal_net = finite(dt * sample.radiation_net_power_w())?;
    report.external_heat_j = finite(report.external_heat_j - internal_net)?;
    report.total_energy_residual_j = finite(report.total_energy_residual_j + internal_net)?;
    Ok(())
}

// A complete provider can add physical thermal links and reservoirs, but this
// field formulation has no owner for an independent graph's electrical work.
fn resolve_heat<A: HeatAssembly + ?Sized>(
    provider: &A, cx: &Cx<'_>, trial: &DepositionSimulation<'_>,
) -> Result<CircuitThermalCoefficients<A::Sample>, HeatError<A::Error>> {
    let coefficients = provider.resolve_heat(cx, trial)
        .map_err(|error| error.map_source(PlanarInductionHeatError::Heat))?;
    if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
        return Err(native_error(invalid("planar induction cannot discard graph electrical work")));
    }
    Ok(coefficients.bulk)
}

fn nonlinear_axes<'a, E: From<HeatNetworkError>>(
    cx: &Cx<'_>, field: &'a PlanarNonlinearEddy<'_>, regions: &[MagneticRegion],
) -> Result<Vec<(fs_thermal::enthalpy::CellId, &'a str)>, HeatError<E>> {
    let mut axes = Vec::new();
    for region in regions {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        if let MagneticRegion::Cell { triangle, cell } = *region {
            if let Some(IsotropicMagneticMaterial::Reversible(law)) =
                field.problem().materials().get(triangle) {
                axes.push((cell, law.source_path().axis.as_str()));
            }
        }
    }
    Ok(axes)
}

/// Continue original midpoint induction with fixed permeability and any complete
/// material heat assembly, including finite coolant. The provider re-resolves
/// the same trial temperatures as sigma(T); no initial conductance is frozen.
/// All original source identities, geometry, magnetic regions, caloric cells
/// and native local Joule heat retain their owners. Passive coolant cells may
/// lie outside the magnetic mesh.
///
/// Graph electrical work and field-dependent thermal boundaries refuse. The
/// original midpoint impressed-current work convention remains unchanged.
/// Full provider samples are in `report().source.network`. As in the standard
/// heat-network report, true external heat excludes
/// `dt * report().source.network.as_ref().radiation_net_power_w()` from
/// `non_electrical_heat_j`; the additional gate enforces that closed balance.
/// Source domains, fully solid field receivers and first-order thermal
/// coupling remain. A refusal advances neither accepted history.
#[allow(clippy::too_many_arguments)]
pub fn fixed_heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A, cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
    field: &PlanarEddy<'m>, impressed_jz_a_m2: &[f64], regions: &[MagneticRegion],
    volume_tolerance_m3: f64, controls: CircuitEnthalpySolve, field_controls: EddyControls,
    initial_energy_increments_j: &[f64],
) -> Result<EddyEnthalpyCandidate<'t, 'm, A::Sample>, HeatError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(assembly_heat_error)?;
    prepared.admit_field_independent_boundaries(cx, &[]).map_err(assembly_heat_error)?;
    let layout = MagneticThermalNetwork {
        cells: bulk.cells, links: bulk.links, magnetic_regions: regions, volume_tolerance_m3,
    };
    let next = layout.step_eddy_resolved_candidate(
        thermal, target_s, field, impressed_jz_a_m2, controls, field_controls,
        initial_energy_increments_j,
        |trial| resolve_heat(provider, cx, trial),
        |sample: &A::Sample| &sample.as_ref().material, &cancel,
    )?;
    closed_balance(target_s - thermal.time_s(), next.report.total_energy_residual_j,
        next.report.source.network.as_ref(), controls.total_energy_tolerance_j)
        .map_err(assembly_heat_error)?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}

/// Continue the original fixed nonlinear H(B) law with the same complete heat
/// provider. Its backward-Euler field step and explicit numerical loss remain
/// distinct from physical Joule heat. Temperature-dependent nonlinear laws
/// still require the explicit entropy receiver below; this entry cannot hide
/// their reversible heat. Reporting, finite coolant and ownership follow
/// `fixed_heat_assembly_candidate`.
#[allow(clippy::too_many_arguments)]
pub fn fixed_nonlinear_heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A, cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
    field: &PlanarNonlinearEddy<'m>, impressed_jz_a_m2: &[f64], regions: &[MagneticRegion],
    volume_tolerance_m3: f64, controls: CircuitEnthalpySolve,
    field_controls: NonlinearEddyControls, initial_energy_increments_j: &[f64],
) -> Result<NonlinearEddyEnthalpyCandidate<'t, 'm, A::Sample>, HeatError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(assembly_heat_error)?;
    let axes = nonlinear_axes(cx, field, regions)?;
    prepared.admit_field_independent_boundaries(cx, &axes).map_err(assembly_heat_error)?;
    let layout = MagneticThermalNetwork {
        cells: bulk.cells, links: bulk.links, magnetic_regions: regions, volume_tolerance_m3,
    };
    let next = layout.step_nonlinear_eddy_resolved_candidate(
        thermal, target_s, field, impressed_jz_a_m2, controls, field_controls,
        initial_energy_increments_j,
        |trial| resolve_heat(provider, cx, trial),
        |sample: &A::Sample| &sample.as_ref().material, &cancel,
    )?;
    closed_balance(target_s - thermal.time_s(), next.report.total_energy_residual_j,
        next.report.source.network.as_ref(), controls.total_energy_tolerance_j)
        .map_err(assembly_heat_error)?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}

/// Continue sourced mu(T), original endpoint-effort induction and the same
/// finite enthalpy inventory with any complete material heat assembly.
///
/// Conductivity, permeability, coolant transport, interfaces, integrated
/// layers and radiation resolve at each proposed endpoint. Native Joule heat,
/// magnetic parameter work and signed reversible heat retain distinct accounts;
/// no transport or field loss is deposited twice. Full provider samples are
/// in `report().coupled.source.material.network`.
///
/// The original magnetic derivative laws, field envelopes, whole volumes and
/// positive combined-capacity admissions remain. Passive coolant cells need
/// no magnetic region. All thermal links require their physical owner; no
/// standalone coolant trajectory or lagged coefficient is introduced. Graph
/// electrical work and field-dependent thermal boundaries refuse. Geometry
/// and source domains stay fixed; the original first-order temporal
/// approximation requires refinement. Both energy gates exclude closed
/// radiation from exterior supply. Typed provider failures publish no state.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A, cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
    field: &PlanarEddy<'m>, impressed_jz_a_m2: &[f64], regions: &[MagneticRegion],
    volume_tolerance_m3: f64, bindings: &[MagneticThermalLaw],
    controls: CircuitEnthalpySolve, field_controls: EddyControls,
    assembly: WeightedAssemblyLimits, initial_energy_increments_j: &[f64],
) -> Result<EddyCaloricCandidate<'t, 'm, A::Sample>, HeatError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(assembly_heat_error)?;
    prepared.admit_field_independent_boundaries(cx, &[]).map_err(assembly_heat_error)?;
    let layout = MagneticThermalNetwork {
        cells: bulk.cells, links: bulk.links, magnetic_regions: regions, volume_tolerance_m3,
    };
    let mut next = layout.step_eddy_thermal_resolved_candidate(
        thermal, target_s, field, impressed_jz_a_m2, bindings, controls,
        field_controls, assembly, initial_energy_increments_j,
        |trial| resolve_heat(provider, cx, trial),
        |sample: &A::Sample| &sample.as_ref().material, &cancel,
    )?;
    close_assembly_heat_balance(&mut next.report, target_s - thermal.time_s(),
        controls.total_energy_tolerance_j)?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}

/// Continue mixed linear mu(T) and explicitly separable H(B,T)=s(T)H_ref(B)
/// laws with the same complete heat assembly. Each original thermal magnetic
/// law owns its field entropy and finite receiving cells. The native nonlinear
/// solver retains its potential/current history, parameter work and separate
/// numerical loss. Positive affine temperature scaling is an authored model,
/// not an inferred general H(B,T) surface. General recoil entropy, hysteresis,
/// changing geometry and phase change remain outside these admitted laws.
/// Coolant, reports, error propagation and energy ownership follow
/// `heat_assembly_candidate`.
#[allow(clippy::too_many_arguments)]
pub fn nonlinear_heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A, cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
    field: &PlanarNonlinearEddy<'m>, impressed_jz_a_m2: &[f64], regions: &[MagneticRegion],
    volume_tolerance_m3: f64,
    bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>],
    controls: CircuitEnthalpySolve, field_controls: NonlinearEddyControls,
    assembly: WeightedAssemblyLimits, initial_energy_increments_j: &[f64],
) -> Result<NonlinearEddyCaloricCandidate<'t, 'm, A::Sample>, HeatError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(assembly_heat_error)?;
    let axes = nonlinear_axes(cx, field, regions)?;
    prepared.admit_field_independent_boundaries(cx, &axes).map_err(assembly_heat_error)?;
    let layout = MagneticThermalNetwork {
        cells: bulk.cells, links: bulk.links, magnetic_regions: regions, volume_tolerance_m3,
    };
    let mut next = layout.step_nonlinear_eddy_thermal_resolved_candidate(
        thermal, target_s, field, impressed_jz_a_m2, bindings, controls,
        field_controls, assembly, initial_energy_increments_j,
        |trial| resolve_heat(provider, cx, trial),
        |sample: &A::Sample| &sample.as_ref().material, &cancel,
    )?;
    close_assembly_heat_balance(&mut next.report, target_s - thermal.time_s(),
        controls.total_energy_tolerance_j)?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}

impl PreparedHeatNetwork<'_, '_> {
    /// Continue the original planar thermal magnetic history through the same
    /// complete `heat_assembly_candidate` entry.
    #[allow(clippy::too_many_arguments)]
    pub fn step_eddy_thermal_candidate<'t, 'm>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
        field: &PlanarEddy<'m>, impressed_jz_a_m2: &[f64], regions: &[MagneticRegion],
        volume_tolerance_m3: f64, bindings: &[MagneticThermalLaw],
        controls: CircuitEnthalpySolve, field_controls: EddyControls,
        assembly: WeightedAssemblyLimits, initial_energy_increments_j: &[f64],
    ) -> Result<EddyCaloricCandidate<'t, 'm, HeatNetworkSample>, Error> {
        heat_assembly_candidate(self, cx, thermal, target_s, field, impressed_jz_a_m2,
            regions, volume_tolerance_m3, bindings, controls, field_controls,
            assembly, initial_energy_increments_j)
            .map_err(|error| error.map_source(concrete_error))
    }

    /// Continue the original nonlinear thermal magnetic history through the
    /// same complete `nonlinear_heat_assembly_candidate` entry.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_eddy_thermal_candidate<'t, 'm>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
        field: &PlanarNonlinearEddy<'m>, impressed_jz_a_m2: &[f64], regions: &[MagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve, field_controls: NonlinearEddyControls,
        assembly: WeightedAssemblyLimits, initial_energy_increments_j: &[f64],
    ) -> Result<NonlinearEddyCaloricCandidate<'t, 'm, HeatNetworkSample>, Error> {
        nonlinear_heat_assembly_candidate(self, cx, thermal, target_s, field, impressed_jz_a_m2,
            regions, volume_tolerance_m3, bindings, controls, field_controls,
            assembly, initial_energy_increments_j)
            .map_err(|error| error.map_source(concrete_error))
    }
}
