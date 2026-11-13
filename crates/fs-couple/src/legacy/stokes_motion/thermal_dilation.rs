//! Instantaneous isobaric expansion from the original native heat balance.
//!
//! The thermal owner supplies actual conduction, reservoir, ambient and authored
//! source watts. Its original fixed parcel mass converts these into MATERIAL
//! `Dh/Dt`; the same caloric chart then supplies `rho * d(1/rho)/dh` to the
//! native Stokes constraint. No separate density, heat capacity or CTE is fitted.
//! Each mapped prism must describe the CURRENT volume of that parcel. Reference
//! heat geometry inconsistent with its current source density refuses.
//!
//! This is a quasistatic response at one accepted thermal state. Neither clock,
//! mass, enthalpy nor geometry advances. A subsequent finite-time composition
//! needs its own material advection and moving-volume owner; this entry does not
//! turn a fixed-mass thermal cell into an Eulerian flow-control volume.
//! Viscous heating and mechanical pressure-work feedback are not added here.

/// The same native coupling with original measured Cp/free-expansion sources.
#[path = "thermal_dilation/isobaric.rs"]
pub mod isobaric;

use std::collections::{BTreeMap, BTreeSet};

use fs_flux::ns::boundary::thermal_dilation::{
    CaloricDilationCell, CaloricDilationError, CaloricDilationSolution,
};
use fs_flux::{FluxParams, FluxSystem, StokesBoundary};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionSimulation, MaterialId, PowerBalance,
    ThermalStorageBasis,
};

use super::heating::StokesThermalCells;
use super::{nonnegative, positive};

/// Original caloric laws and current triangle prisms, all in the mesh's SI frame.
#[derive(Clone, Copy)]
pub struct StokesThermalDilationBinding<'a> {
    pub cells: StokesThermalCells<'a>,
    pub depth_m: f64,
    /// Exact original charts, including their enthalpy datum and density law.
    /// Additional unmapped thermal cells may represent walls or heat sinks.
    pub caloric_by_material: &'a BTreeMap<MaterialId, &'a EquilibriumEnthalpyPhaseCurve>,
}

/// Independent local continuity and whole-boundary volume allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesThermalDilationControls {
    pub maximum_dilation_defect_per_s: f64,
    pub maximum_volume_rate_defect_m3_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum StokesThermalDilationError {
    Binding(&'static str),
    Thermal(DepositionError),
    Flow(CaloricDilationError),
    VolumeBalance {
        expected_m3_s: f64,
        boundary_m3_s: f64,
        defect_m3_s: f64,
        maximum_m3_s: f64,
    },
    Cancelled,
}

impl core::fmt::Display for StokesThermalDilationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native thermal/Stokes dilation: {self:?}")
    }
}
impl std::error::Error for StokesThermalDilationError {}
impl From<DepositionError> for StokesThermalDilationError {
    fn from(error: DepositionError) -> Self {
        match error {
            DepositionError::Cancelled => Self::Cancelled,
            error => Self::Thermal(error),
        }
    }
}
impl From<CaloricDilationError> for StokesThermalDilationError {
    fn from(error: CaloricDilationError) -> Self {
        Self::Flow(error)
    }
}

/// Original thermal identity and physical heat/volume rates, in triangle order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesThermalDilationCell {
    pub cell: CellId,
    pub caloric_identity: [u8; 32],
    pub mass_kg: f64,
    pub heat_power_w: f64,
    pub material_specific_enthalpy_rate_w_kg: f64,
    /// `Q * dv/dh` for this conserved parcel [m³/s], signed outward.
    pub volume_rate_m3_s: f64,
}

/// Both original native results, without publishing a thermal or flow timestep.
pub struct StokesThermalDilation {
    pub time_s: f64,
    /// Complete native thermal balance, including any unmapped walls/sinks.
    pub heat: PowerBalance,
    pub flow: CaloricDilationSolution,
    pub cells: Vec<StokesThermalDilationCell>,
    /// Sum over mapped liquid cells only, from original heat powers and dv/dh.
    pub expected_volume_rate_m3_s: f64,
    /// Independently read from returned BDM moments on the physical boundary.
    pub boundary_volume_rate_m3_s: f64,
    pub volume_rate_defect_m3_s: f64,
}

/// Solve the original Stokes response to this state's native material heating.
///
/// Each distinct mapped cell must use enthalpy storage, the exact supplied
/// caloric identity, a fully liquid current state, and the original triangle
/// prism with `mass = rho(h) * area * depth`. These are arithmetic identity
/// checks, not tunable physical tolerances or a geometry/mass repair.
///
/// Only authored native heat sources, links, ambient and reservoir exchange
/// are evaluated, via `power_balance_observed(&())`. Arbitrary observers are
/// deliberately not accepted: an advective enthalpy rate is an Eulerian
/// transport term, not this parcel's material heating. Pending events are not
/// applied and interval-dependent heat is not estimated from two endpoints.
///
/// The native viscosity/force/traction normalization and effective pressure
/// convention are unchanged; see `FluxSystem::assemble_stokes_caloric_dilation`.
/// Outflow boundaries can release thermal volume. A rigid sealed component
/// with incompatible net dilation refuses; pressure is not an invented EOS.
/// No closed thermal-mechanical energy balance is claimed by this one-way
/// isobaric snapshot. Refining the thermal and flow discretizations remains
/// independent. Preparation polls cancellation and reuses the thermal budget;
/// the synchronous native Stokes solve is polled only before and after it.
#[allow(clippy::too_many_arguments)]
pub fn solve_at_thermal_state(
    system: &FluxSystem<'_>,
    thermal: &DepositionSimulation<'_>,
    binding: StokesThermalDilationBinding<'_>,
    params: FluxParams,
    viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
    force: &dyn Fn([f64; 2]) -> [f64; 2],
    boundaries: &[StokesBoundary<'_>],
    controls: StokesThermalDilationControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<StokesThermalDilation, StokesThermalDilationError> {
    poll(&mut cancelled)?;
    let mesh = system.mesh;
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
    let mut sources = Vec::with_capacity(mesh.tris.len());
    let mut cells = Vec::with_capacity(mesh.tris.len());
    for (triangle, &id) in binding.cells.triangle_cells.iter().enumerate() {
        if triangle % 256 == 0 {
            poll(&mut cancelled)?;
        }
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(invalid("inactive mapped thermal cell"))?;
        if !seen.insert(id) {
            return Err(invalid("duplicate mapped thermal cell"));
        }
        let &caloric = binding
            .caloric_by_material
            .get(&cell.material_id)
            .ok_or(invalid("missing original caloric chart for liquid cell"))?;
        let identity = *caloric.identity().as_bytes();
        let state = caloric
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|source| CaloricDilationError::Material {
                cell: triangle,
                source,
            })?;
        if state.phase() != SolidLiquidPhase::Liquid {
            return Err(CaloricDilationError::NonLiquid { cell: triangle }.into());
        }
        let receipt = admit_current_prism(
            system,
            thermal,
            cell,
            triangle,
            binding.cells,
            depth,
            &heat,
            CurrentPrismSource {
                identity,
                temperature_k: state.temperature_k(),
                density_kg_m3: state.bulk_density_kg_m3(),
            },
        )?;
        sources.push(CaloricDilationCell {
            caloric,
            specific_enthalpy_j_kg: cell.specific_h_j_kg,
            material_specific_enthalpy_rate_w_kg: receipt.material_specific_enthalpy_rate_w_kg,
        });
        cells.push(receipt);
    }
    poll(&mut cancelled)?;
    let assembly =
        system.assemble_stokes_caloric_dilation(params, viscosity, force, boundaries, &sources)?;
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
    Ok(StokesThermalDilation {
        time_s: thermal.time_s(),
        heat,
        flow,
        cells,
        expected_volume_rate_m3_s: expected,
        boundary_volume_rate_m3_s: boundary,
        volume_rate_defect_m3_s: defect,
    })
}

fn admit_binding(
    system: &FluxSystem<'_>,
    thermal: &DepositionSimulation<'_>,
    cells: StokesThermalCells<'_>,
    depth: f64,
    controls: StokesThermalDilationControls,
) -> Result<(), StokesThermalDilationError> {
    let end_z = cells.plane_z_m + depth;
    if cells.triangle_cells.len() != system.mesh.tris.len()
        || system.mesh.tris.len() > thermal.cells().len()
        || !positive(depth)
        || !cells.plane_z_m.is_finite()
        || !same(end_z - cells.plane_z_m, depth)
        || !nonnegative(controls.maximum_dilation_defect_per_s)
        || !nonnegative(controls.maximum_volume_rate_defect_m3_s)
    {
        return Err(invalid(
            "complete current prism map and finite physical allowances required",
        ));
    }
    Ok(())
}

struct CurrentPrismSource {
    identity: [u8; 32],
    temperature_k: f64,
    density_kg_m3: f64,
}

#[allow(clippy::too_many_arguments)]
fn admit_current_prism(
    system: &FluxSystem<'_>,
    thermal: &DepositionSimulation<'_>,
    cell: &DepositedCell,
    triangle: usize,
    mapping: StokesThermalCells<'_>,
    depth: f64,
    heat: &PowerBalance,
    source: CurrentPrismSource,
) -> Result<StokesThermalDilationCell, StokesThermalDilationError> {
    if thermal.material_storage_basis(cell.material_id)? != ThermalStorageBasis::Enthalpy
        || thermal.material_caloric_identity(cell.material_id)? != Some(source.identity)
    {
        return Err(invalid(
            "mapped liquid requires its exact original enthalpy chart",
        ));
    }
    let mesh = system.mesh;
    let [x, y] = mesh.centroids[triangle];
    let start = [x, y, mapping.plane_z_m];
    let end = [x, y, mapping.plane_z_m + depth];
    let volume = mesh.areas[triangle] * depth;
    let mass = source.density_kg_m3 * volume;
    let perimeter: f64 = mesh.tri_edges[triangle]
        .iter()
        .map(|&(edge, _)| mesh.edges[edge].len)
        .sum();
    let surface = 2.0 * mesh.areas[triangle] + perimeter * depth;
    if !positive(volume)
        || !positive(mass)
        || !positive(surface)
        || !same(cell.geometry.volume_m3, volume)
        || !same(cell.geometry.surface_area_m2, surface)
        || !same(cell.mass_kg, mass)
        || !(0..3).all(|i| {
            same(cell.geometry.start_m[i], start[i]) && same(cell.geometry.end_m[i], end[i])
        })
        || !same(heat.temperatures_k[&cell.id], source.temperature_k)
    {
        return Err(invalid(
            "current liquid mass, temperature and geometry must match its original caloric prism",
        ));
    }
    let power = heat.net_power_w[&cell.id];
    let rate = power / cell.mass_kg;
    if !rate.is_finite() || (power != 0.0 && rate == 0.0) {
        return Err(invalid("unrepresentable material enthalpy rate"));
    }
    Ok(StokesThermalDilationCell {
        cell: cell.id,
        caloric_identity: source.identity,
        mass_kg: cell.mass_kg,
        heat_power_w: power,
        material_specific_enthalpy_rate_w_kg: rate,
        volume_rate_m3_s: 0.0,
    })
}

fn volume_balance(
    system: &FluxSystem<'_>,
    depth: f64,
    cells: &mut [StokesThermalDilationCell],
    slopes: impl Iterator<Item = f64>,
    field: &[f64],
    maximum_defect: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64, f64), StokesThermalDilationError> {
    let mesh = system.mesh;
    let mut expected = Dd::from_f64(0.0);
    for (index, (cell, slope)) in cells.iter_mut().zip(slopes).enumerate() {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        let rate = Dd::from_f64(cell.heat_power_w) * Dd::from_f64(slope);
        cell.volume_rate_m3_s = rate.to_f64();
        if !cell.volume_rate_m3_s.is_finite()
            || (cell.heat_power_w != 0.0 && slope != 0.0 && cell.volume_rate_m3_s == 0.0)
        {
            return Err(invalid("unrepresentable original thermal volume rate"));
        }
        expected = expected + rate;
    }
    let mut boundary = Dd::from_f64(0.0);
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        if edge % 256 == 0 {
            poll(cancelled)?;
        }
        if mesh.is_boundary(edge) {
            // Boundary edge normals are outward. The linear BDM moment has
            // zero mean, leaving only the constant moment in the flux integral.
            boundary = boundary
                + Dd::from_f64(depth) * Dd::from_f64(geometry.len) * Dd::from_f64(field[2 * edge]);
        }
    }
    let defect = (boundary - expected).to_f64().abs();
    if !expected.to_f64().is_finite()
        || !boundary.to_f64().is_finite()
        || !defect.is_finite()
        || defect > maximum_defect
    {
        return Err(StokesThermalDilationError::VolumeBalance {
            expected_m3_s: expected.to_f64(),
            boundary_m3_s: boundary.to_f64(),
            defect_m3_s: defect,
            maximum_m3_s: maximum_defect,
        });
    }
    poll(cancelled)?;
    Ok((expected.to_f64(), boundary.to_f64(), defect))
}

fn invalid(reason: &'static str) -> StokesThermalDilationError {
    StokesThermalDilationError::Binding(reason)
}
fn same(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 64.0 * f64::EPSILON * a.abs().max(b.abs())
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), StokesThermalDilationError> {
    if cancelled() {
        Err(StokesThermalDilationError::Cancelled)
    } else {
        Ok(())
    }
}
