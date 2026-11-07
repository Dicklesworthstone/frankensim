//! Explicit pressure openings inside the existing occupied-voxel pressure graph.

use super::*;
use crate::d3q19::Face3;

/// One liquid face adjoining an inactive voxel, with an authored mechanical
/// correction pressure [Pa]. The inactive neighbour supplies no inertia or
/// transported material. Unlisted inactive faces remain impermeable walls.
///
/// The opening is the complete voxel face. Pressure acts across the occupied
/// half-cell, using its physical density and actual Cartesian width. This is
/// a staircase boundary, not a cut-cell aperture or an inferred VOF surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalPressureFace {
    /// Occupied fluid cell, in the original domain's linear indexing.
    pub cell: usize,
    /// Outward face of that fluid cell.
    pub face: Face3,
    /// Mechanical pressure in the same gauge as exterior Free faces.
    pub pressure_pa: f64,
}

fn admit(
    domain: &VoxelDomain,
    openings: &[InternalPressureFace],
    gate: &CancelGate,
) -> Result<BTreeMap<(usize, usize), f64>, ChtError> {
    poll(gate)?;
    if openings.len() > domain.cell_count().saturating_mul(6) {
        return Err(refuse(
            "flow.internal_pressure_faces",
            "too many pressure faces",
        ));
    }
    let mut pressures = BTreeMap::new();
    for (i, opening) in openings.iter().enumerate() {
        poll_every(gate, i)?;
        let f = opening.face as usize;
        if opening.cell >= domain.cell_count()
            || !domain.is_fluid(opening.cell)
            || !domain
                .neighbor(opening.cell, f)
                .is_some_and(|n| !domain.is_fluid(n))
        {
            return Err(refuse(
                "flow.internal_pressure_role",
                "an internal pressure face must point from fluid to an inactive voxel",
            ));
        }
        finite("flow.internal_pressure_pa", opening.pressure_pa)?;
        if pressures
            .insert((opening.cell, f), opening.pressure_pa)
            .is_some()
        {
            return Err(refuse(
                "flow.internal_pressure_faces",
                "duplicate internal pressure face",
            ));
        }
    }
    Ok(pressures)
}

impl FlowField {
    /// Project physical volume sources with explicit internal pressure openings.
    ///
    /// Reuses the original density-weighted graph, component compatibility,
    /// pressure gauge, Krylov solve and independent flux/work audit. The gas or
    /// other inactive side has no density unknown; the opening's inertia is the
    /// occupied half-cell mass. Positive source is outward volume demand, not
    /// injected physical mass. Each sealed component must balance its own
    /// source and held exterior flux; a pressure opening removes that component's
    /// nullspace and carries its actual signed flux.
    ///
    /// Exterior Free planes prescribe zero correction pressure. Internal face
    /// pressures may vary; source and boundary work use the first open pressure
    /// in each connected component as their shared gauge. Large common pressure
    /// offsets therefore do not enter the small-head solve or its work defect.
    ///
    /// Returned flux across an opening is a kinematic volume flux. It does not
    /// create a gas material state or perform VOF, mass, enthalpy or topology
    /// transport. Ordinary closed-solid transport cannot consume these openings
    /// without its own explicit interface owner. Inputs are never mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn project_density_volume_source_with_internal_pressure_faces(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        openings: &[InternalPressureFace],
        density_kg_m3: &[f64],
        target_outflow_m3_s: &[f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<DensityProjection, ChtError> {
        let pressures = admit(domain, openings, gate)?;
        self.project_density_volume_source_with_reservoirs(
            domain,
            faces,
            BoundaryPressure::Internal(&pressures),
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )
    }

    /// Interpolate an unprojected physical velocity predictor on the same
    /// occupied domain and internal pressure openings. Interior faces use the
    /// two original centre distances; every pressure opening uses only its
    /// adjacent liquid velocity. Inactive velocities are never read.
    ///
    /// This makes no divergence-free correction before the actual volume-source
    /// projection. Normal velocities retain the usual positive-axis convention.
    pub fn from_cell_velocity_predictor_with_internal_pressure_faces(
        domain: &VoxelDomain,
        velocities_m_s: &[[f64; 3]],
        faces: [FlowFace; 6],
        openings: &[InternalPressureFace],
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        let pressures = admit(domain, openings, gate)?;
        if velocities_m_s.len() != domain.cell_count() {
            return Err(refuse(
                "flow.velocities",
                "velocities must match the original domain",
            ));
        }
        let mut field = Self::quiescent(domain);
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            for value in velocities_m_s[c] {
                finite("flow.velocities", value)?;
            }
            for f in 0..6 {
                let axis = f / 2;
                let velocity = match domain.neighbor(c, f) {
                    Some(n) if f % 2 == 1 && domain.is_fluid(n) => {
                        finite("flow.velocities", velocities_m_s[n][axis])?;
                        if domain.is_uniform() {
                            0.5 * (velocities_m_s[c][axis] + velocities_m_s[n][axis])
                        } else {
                            let a = 0.5 * domain.widths(c)[axis];
                            let b = 0.5 * domain.widths(n)[axis];
                            let fraction = a / (a + b);
                            (1.0 - fraction) * velocities_m_s[c][axis]
                                + fraction * velocities_m_s[n][axis]
                        }
                    }
                    None if faces[f] != FlowFace::Wall => velocities_m_s[c][axis],
                    Some(_) if pressures.contains_key(&(c, f)) => velocities_m_s[c][axis],
                    _ => continue,
                };
                let flux = velocity * domain.face_area(c, axis);
                finite("flow.interpolated_flux", flux)?;
                let (slot, _) = field.slot(domain, c, f);
                *field.face_mut(axis, slot) = flux;
            }
        }
        poll(gate)?;
        Ok(field)
    }
}

#[cfg(test)]
mod tests;
