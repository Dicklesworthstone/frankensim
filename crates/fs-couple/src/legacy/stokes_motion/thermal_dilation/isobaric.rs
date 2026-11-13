//! Native material heating drives the same original measured expansion isobar.

use super::{
    CurrentPrismSource, StokesThermalDilationCell, StokesThermalDilationControls,
    StokesThermalDilationError, admit_binding, admit_current_prism, invalid, poll, volume_balance,
};
use crate::stokes_motion::heating::StokesThermalCells;
use fs_flux::ns::boundary::thermal_dilation::isobaric::{
    IsobaricDilationCell, IsobaricDilationError, IsobaricDilationSolution,
};
use fs_flux::{FluxParams, FluxSystem, StokesBoundary};
use fs_material::thermal::transport::isobaric::expansion::liquid::LiquidIsobaricExpansion;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, MaterialId, PowerBalance};
use std::collections::{BTreeMap, BTreeSet};

/// Explicitly liquid original source paths, registered on the thermal timeline
/// as each adapter's `caloric()`, and CURRENT material triangle-prism geometry.
#[derive(Clone, Copy)]
pub struct StokesIsobaricDilationBinding<'a> {
    pub cells: StokesThermalCells<'a>,
    pub depth_m: f64,
    pub caloric_by_material: &'a BTreeMap<MaterialId, &'a LiquidIsobaricExpansion<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum StokesIsobaricDilationError {
    Admission(StokesThermalDilationError),
    Flow(IsobaricDilationError),
}
impl core::fmt::Display for StokesIsobaricDilationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original-source thermal/Stokes dilation: {self:?}")
    }
}
impl std::error::Error for StokesIsobaricDilationError {}
impl From<StokesThermalDilationError> for StokesIsobaricDilationError {
    fn from(error: StokesThermalDilationError) -> Self {
        Self::Admission(error)
    }
}
impl From<IsobaricDilationError> for StokesIsobaricDilationError {
    fn from(error: IsobaricDilationError) -> Self {
        Self::Flow(error)
    }
}
impl From<DepositionError> for StokesIsobaricDilationError {
    fn from(error: DepositionError) -> Self {
        Self::Admission(error.into())
    }
}

/// Original native heat/flow results, with independent source and boundary
/// volume rates. This contains no proposed thermal, mass or geometry timestep.
pub struct StokesIsobaricDilation {
    pub time_s: f64,
    pub heat: PowerBalance,
    pub flow: IsobaricDilationSolution,
    pub cells: Vec<StokesThermalDilationCell>,
    pub expected_volume_rate_m3_s: f64,
    pub boundary_volume_rate_m3_s: f64,
    pub volume_rate_defect_m3_s: f64,
}

/// Solve native outlet flow from the current state's original isobaric laws.
///
/// The thermal owner supplies authored native link, reservoir, ambient and
/// source powers. `Q/m` is the material enthalpy rate; Cp already accounts for
/// maintained-pressure boundary work, so p*dV is never added again. The same
/// original Cp integral/inverse and engineering-length source supply current
/// density and the exact directional dv/dh, with no replacement phase chart.
/// The liquid-path constructor records caller-declared physical applicability;
/// it does not turn solid expansion evidence into liquid phase evidence.
///
/// The original expanding-conduction owner may share this exact caloric law;
/// this entry leaves heat-transfer coefficients with the native heat owner.
/// As with `super::solve_at_thermal_state`, arbitrary/advective observers are
/// excluded, each fixed parcel must currently satisfy m=rho(h)*area*depth,
/// and original source/storage identities must agree. Connected flow cells
/// share one original maintained source pressure. Outward domain rates
/// and incompatible sealed heating refuse. Native stress/pressure and reference
/// density conventions are unchanged. This is an instantaneous one-way response,
/// without finite-time advection, ALE, pressure-work or viscous-heat feedback.
#[allow(clippy::too_many_arguments)]
pub fn solve_at_thermal_state(
    system: &FluxSystem<'_>,
    thermal: &DepositionSimulation<'_>,
    binding: StokesIsobaricDilationBinding<'_>,
    params: FluxParams,
    viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
    force: &dyn Fn([f64; 2]) -> [f64; 2],
    boundaries: &[StokesBoundary<'_>],
    controls: StokesThermalDilationControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<StokesIsobaricDilation, StokesIsobaricDilationError> {
    poll(&mut cancelled)?;
    let depth = binding.depth_m;
    admit_binding(system, thermal, binding.cells, depth, controls)?;
    let heat = thermal.power_balance_observed(&(), || {
        if cancelled() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    })?;
    let mut seen = BTreeSet::new();
    let mut sources = Vec::with_capacity(system.mesh.tris.len());
    let mut cells = Vec::with_capacity(system.mesh.tris.len());
    for (triangle, &id) in binding.cells.triangle_cells.iter().enumerate() {
        if triangle % 256 == 0 {
            poll(&mut cancelled)?;
        }
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(invalid("inactive mapped thermal cell"))?;
        if !seen.insert(id) {
            return Err(invalid("duplicate mapped thermal cell").into());
        }
        let &caloric = binding
            .caloric_by_material
            .get(&cell.material_id)
            .ok_or(invalid("missing original isobaric source for liquid cell"))?;
        let state = caloric
            .dilation_at_specific_enthalpy(
                cell.specific_h_j_kg,
                heat.net_power_w[&id] / cell.mass_kg,
            )
            .map_err(|source| IsobaricDilationError::Material {
                cell: triangle,
                source,
            })?;
        let receipt = admit_current_prism(
            system,
            thermal,
            cell,
            triangle,
            binding.cells,
            depth,
            &heat,
            CurrentPrismSource {
                identity: *state.caloric_identity.as_bytes(),
                temperature_k: state.state.temperature_k(),
                density_kg_m3: state.state.density_kg_m3(),
            },
        )?;
        sources.push(IsobaricDilationCell {
            caloric,
            specific_enthalpy_j_kg: cell.specific_h_j_kg,
            material_specific_enthalpy_rate_w_kg: receipt.material_specific_enthalpy_rate_w_kg,
        });
        cells.push(receipt);
    }
    poll(&mut cancelled)?;
    let assembly =
        system.assemble_stokes_isobaric_dilation(params, viscosity, force, boundaries, &sources)?;
    poll(&mut cancelled)?;
    let flow = assembly.solve(controls.maximum_dilation_defect_per_s)?;
    poll(&mut cancelled)?;
    let (expected, boundary, defect) = volume_balance(
        system,
        depth,
        &mut cells,
        flow.source
            .iter()
            .map(|source| source.specific_volume_derivative_m3_j),
        &flow.flow.field.x,
        controls.maximum_volume_rate_defect_m3_s,
        &mut cancelled,
    )?;
    Ok(StokesIsobaricDilation {
        time_s: thermal.time_s(),
        heat,
        flow,
        cells,
        expected_volume_rate_m3_s: expected,
        boundary_volume_rate_m3_s: boundary,
        volume_rate_defect_m3_s: defect,
    })
}
