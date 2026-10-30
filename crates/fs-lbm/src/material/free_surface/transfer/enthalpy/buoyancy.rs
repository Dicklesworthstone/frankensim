//! The original bounded source projection and actual receiving force frame.
//!
//! Geometric/caloric mass uses rho0. The source rho(T) stays in the material
//! owner; only its gravity response is projected through native Boussinesq.

use super::*;
use crate::material::MaterialLbmError;
use crate::material::free_surface::{MaterialFreeSurface3, MaterialFreeSurfaceInitialCell3};
use crate::material::transport::boussinesq::BoussinesqMaterial;
use fs_material::fluid_transport::FluidTransport;
use fs_material::liquid::ResolvedLiquidState;

pub(super) fn admit_source_transport(
    material: ResolvedLiquidState,
    density: f64,
    law: &LiquidTemperatureLaw<'_>,
    config: BoussinesqConfig<3>,
) -> Result<(), TetEnthalpyTransferError> {
    let reference = law
        .reference_state()
        .transport_state()
        .map_err(MaterialLbmError::from)
        .map_err(FreeSurfaceEnthalpyError::from)?;
    if reference.density_kg_m3().to_bits() != density.to_bits() {
        return Err(invalid(
            "tetrahedral Boussinesq mass needs the original liquid reference density",
        )
        .into());
    }
    // Exactly the shared owner used by native thermal momentum: reference
    // rho/Cp, this actual source's dynamic viscosity and conductivity, and
    // separately bounded unprojected source departures. No second rho law.
    BoussinesqMaterial::try_new(material, &reference, config)
        .map_err(MaterialLbmError::from)
        .map_err(FreeSurfaceEnthalpyError::from)?;
    Ok(())
}

/// Re-read force-corrected native velocity after the original Guo rebind.
/// Reuse the same per-fragment projection/error budget, including its absolute
/// sums, rather than trusting an earlier initializer's unforced report.
pub(super) fn admit_receiver_moments(
    flow: &MaterialFreeSurface3<'_>,
    report: &mut TetLiquidTransferReport,
    dx: f64,
    gravity: [f64; 3],
    limits: TetLiquidTransferLimits,
    cx: &Cx<'_>,
) -> Result<(), TetEnthalpyTransferError> {
    use super::super::{
        LiquidTransferMoments, accumulate_cell, centre, check_projection, components,
    };

    let mut original = Vec::with_capacity(report.components.len());
    let mut actual = Vec::with_capacity(report.components.len());
    let mut initial = Vec::with_capacity(report.wet_cells);
    let mut total = LiquidTransferMoments::default();
    let cell_volume = dx * dx * dx;
    for component in &report.components {
        checkpoint(cx)?;
        let mut target = LiquidTransferMoments::default();
        let mut cell_indices = Vec::with_capacity(component.target_cells.len());
        for &coordinates in &component.target_cells {
            checkpoint(cx)?;
            let [x, y, z] = coordinates;
            let cell = flow
                .sample_at_si(x, y, z)
                .map_err(TetLiquidTransferError::from)?
                .ok_or_else(|| invalid("buoyancy attachment lost an original liquid cell"))?;
            let mass = cell.density_kg_m3 * cell.fill_fraction * cell_volume;
            let position = centre(coordinates, dx);
            accumulate_cell(&mut target, mass, position, cell.velocity_m_s);
            accumulate_cell(&mut total, mass, position, cell.velocity_m_s);
            cell_indices.push(initial.len());
            initial.push(MaterialFreeSurfaceInitialCell3 {
                coordinates,
                fill_fraction: cell.fill_fraction,
                velocity_m_s: cell.velocity_m_s,
            });
        }
        original.push(components::Component {
            source_tetrahedra: component.source_tetrahedra.clone(),
            cell_indices,
            source: component.source,
        });
        actual.push(target);
    }
    if (flow
        .ledger_mass_kg()
        .map_err(TetLiquidTransferError::from)?
        - total.mass_kg)
        .abs()
        > limits.mass_tolerance_kg
    {
        return Err(invalid("buoyancy attachment changed the original native mass ledger").into());
    }
    let (centroid_shift_m, kinetic_energy_change_j, gravitational_energy_change_j) =
        check_projection(&report.source, &total, gravity, limits)?;
    let components = components::reports(original, actual, &initial, gravity, limits, cx)?;
    checkpoint(cx)?;
    report.target = total;
    report.centroid_shift_m = centroid_shift_m;
    report.kinetic_energy_change_j = kinetic_energy_change_j;
    report.gravitational_energy_change_j = gravitational_energy_change_j;
    report.components = components;
    Ok(())
}
