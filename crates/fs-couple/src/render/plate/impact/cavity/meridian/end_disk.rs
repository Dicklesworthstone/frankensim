//! Reciprocal coupling of an actual tensioned film to a chamber end disk.
//! Positive film displacement is explicitly either +z or -z. The overlap is
//! outward from the gas, and the same column is used for pressure work and
//! swept volume by CavityCoupling. No equivalent-volume cylinder is introduced.
use super::{MeridianCavity, ImpactError, invalid};
use crate::vibroacoustic::{CavityModes, StructuralModes, assemble_coupling};
use fs_exec::CancelGate;
use fs_plate::{ModePair, shell::head::TensionedDisk};

/// Which end of the ordered meridian is covered by the film.
#[derive(Clone, Copy, Debug)]
pub enum EndDisk { Lower, Upper }

/// Addresses and direction of the supplied mass-normalized film basis.
#[derive(Clone, Copy, Debug)]
pub struct DiskPort {
    pub end: EndDisk,
    /// Mechanical coordinate of the first film mode, before q/p interleaving.
    pub first_mode: usize,
    /// All structural coordinates, including strikers; not acoustic inertia.
    pub structural_modes: usize,
    /// +1 for displacement toward +z, -1 for displacement toward -z.
    pub positive_z: f64,
}

/// Cold quadrature and storage ceilings. Neither is a convergence certificate.
#[derive(Clone, Copy, Debug)]
pub struct InterfaceBudget {
    pub maximum_points: usize,
    /// Bounds (structural modes + pressure modes) times quadrature points.
    pub maximum_terms: usize,
}

/// Same-basis pressure carrier and structural-by-pressure overlap matrix.
/// Unrelated structural coordinates have exactly zero overlap. Callers may
/// sum overlaps of several moving walls before constructing CavityCoupling.
pub struct EndDiskInterface {
    pub cavity: CavityModes,
    pub overlaps: Vec<f64>,
}

impl MeridianCavity {
    /// Integrate P1 film motion times the actual meridional pressure basis.
    /// Three positive quadrature points per head triangle approximate the
    /// pressure/film overlap. Radial interpolation, Fourier factors and crossing
    /// pressure cells require independent interface refinement. The polygonal film
    /// approaches the circular Fourier boundary as azimuthal resolution grows.
    /// Radius mismatches refuse instead of stretching, snapping or rescaling.
    pub fn end_disk_interface(&self, film: &TensionedDisk, modes: &[ModePair],
        port: DiskPort, budget: InterfaceBudget, gate: &CancelGate)
        -> Result<EndDiskInterface, ImpactError>
    {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        let station = match port.end {
            EndDisk::Lower => self.profile[0],
            EndDisk::Upper => self.profile[self.profile.len()-1],
        };
        let end = port.first_mode.checked_add(modes.len());
        let points = film.mesh.tris.len().checked_mul(3);
        let terms = port.structural_modes.checked_add(self.modes.len())
            .and_then(|n| points.and_then(|p| n.checked_mul(p)));
        if modes.is_empty() || port.structural_modes == 0
            || port.structural_modes > crate::modal_acoustic_time::MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || end.is_none_or(|n| n > port.structural_modes)
            || (port.positive_z != -1.0 && port.positive_z != 1.0)
            || station.radius_m <= 0.0
            || !film.spec.radius_m.is_finite() || film.spec.radius_m <= 0.0
            || (film.spec.radius_m-station.radius_m).abs() > 64.0*f64::EPSILON*station.radius_m
            || points.is_none_or(|n| n == 0 || n > budget.maximum_points)
            || terms.is_none_or(|n| n > budget.maximum_terms)
            || film.mesh.nodes.len().checked_mul(3) != Some(film.model.dof_map.len())
            || film.model.dof_map.iter().flatten().any(|&d| d >= film.model.free)
            || film.mesh.tris.iter().flatten().any(|&i| i >= film.mesh.nodes.len())
            || film.mesh.nodes.iter().any(|&(x,y)| !x.is_finite() || !y.is_finite())
            || modes.iter().any(|m| !m.lambda.is_finite() || m.lambda <= 0.0
                || m.phi.len() != film.model.free || m.phi.iter().any(|v| !v.is_finite()))
        { return Err(invalid("end disk needs matching geometry, complete finite modes, signed addresses and quadrature budgets")); }
        let points = points.expect("admitted quadrature count");
        let outward = port.positive_z * match port.end { EndDisk::Lower => -1.0, EndDisk::Upper => 1.0 };
        let mut positions = Vec::with_capacity(points);
        let mut areas = Vec::with_capacity(points);
        let mut shapes = vec![vec![0.0;points];port.structural_modes];
        for tri in &film.mesh.tris {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            let nodes = tri.map(|i| film.mesh.nodes[i]);
            let area = ((nodes[1].0-nodes[0].0)*(nodes[2].1-nodes[0].1)
                -(nodes[2].0-nodes[0].0)*(nodes[1].1-nodes[0].1))/2.0;
            if !area.is_finite() || area <= 0.0 { return Err(invalid("end disk has a degenerate or reversed triangle")); }
            for bary in [[2.0/3.0,1.0/6.0,1.0/6.0], [1.0/6.0,2.0/3.0,1.0/6.0], [1.0/6.0,1.0/6.0,2.0/3.0]] {
                let point = positions.len();
                positions.push([nodes.iter().zip(bary).map(|(n,b)| n.0*b).sum(),
                    nodes.iter().zip(bary).map(|(n,b)| n.1*b).sum(), station.z_m]);
                areas.push(area/3.0);
                for (j,mode) in modes.iter().enumerate() {
                    shapes[port.first_mode+j][point] = outward * tri.iter().zip(bary)
                        .map(|(&node,b)| b*film.model.dof_map[3*node].map_or(0.0,|d| mode.phi[d])).sum::<f64>();
                }
            }
        }
        let cavity = self.sample(&positions, budget.maximum_terms, gate)?;
        let mut omegas = vec![0.0;port.structural_modes];
        for (j,m) in modes.iter().enumerate() { omegas[port.first_mode+j] = m.lambda.sqrt(); }
        let structure = StructuralModes { omegas, shapes, loss_factor: 0.0 };
        let overlaps = assemble_coupling(&structure, &cavity, &areas)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok(EndDiskInterface { cavity, overlaps })
    }
}

#[cfg(test)]
mod tests;
