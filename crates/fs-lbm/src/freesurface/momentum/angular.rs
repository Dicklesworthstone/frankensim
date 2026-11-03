//! Scalar spin closure on actual local D2Q9 conversion routes.
//!
//! Cell-centre relocation changes first mass moments. This numerical remap
//! closes angular momentum in the declared lattice frame; it is not a
//! Galilean-invariant geometric transfer or a whole-timestep conservation law.
//! Periodic coordinates are refused by the common momentum admission owner.

use super::*;

#[derive(Clone, Copy)]
struct Snapshot {
    mass: f64,
    velocity: [f64; 2],
}

struct Node {
    before: Snapshot,
    parent: usize,
    transferred_mass: f64,
}

#[derive(Default)]
pub(super) struct AngularConversion2 {
    nodes: BTreeMap<usize, Node>,
}

struct Participant {
    cell: usize,
    before: Snapshot,
    after: Snapshot,
    position: [f64; 2],
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}

fn snapshot(surface: &FreeSurface, cell: usize) -> Result<Snapshot, FreeSurfaceMomentumError2> {
    if surface.grid.flags[cell] == Cell::Gas {
        return Ok(Snapshot {
            mass: 0.0,
            velocity: [0.0; 2],
        });
    }
    if !matches!(surface.grid.flags[cell], Cell::Fluid | Cell::Interface) {
        return Err(invalid(
            "angular conversion participant is not liquid or gas",
        ));
    }
    let state = surface.grid.moments(cell);
    let mass = if surface.grid.flags[cell] == Cell::Fluid {
        state.rho
    } else {
        surface.mass[cell]
    };
    if !state.rho.is_finite()
        || state.rho <= 0.0
        || !mass.is_finite()
        || mass < 0.0
        || mass > state.rho
    {
        return Err(invalid(
            "angular conversion needs bounded final cell masses",
        ));
    }
    admit_velocity(state.u)?;
    Ok(Snapshot {
        mass,
        velocity: state.u,
    })
}

// Difference local inventories first, avoiding cancellation between large
// absolute angular totals. The origin remains fixed throughout conversion.
fn difference(
    cells: &[Participant],
    gate: &CancelGate,
) -> Result<(f64, [f64; 2], f64), FreeSurfaceMomentumError2> {
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 2], 0.0);
    for (index, cell) in cells.iter().enumerate() {
        if index % 256 == 0 {
            poll(gate)?;
        }
        let dm = cell.after.mass - cell.before.mass;
        let dp = core::array::from_fn(|d| {
            cell.after.mass.mul_add(
                cell.after.velocity[d] - cell.before.velocity[d],
                dm * cell.before.velocity[d],
            )
        });
        mass += dm;
        for d in 0..2 {
            momentum[d] += dp[d];
        }
        angular += cross(cell.position, dp);
    }
    Ok((mass, momentum, angular))
}

fn offset(coordinate: usize, origin: usize) -> f64 {
    if coordinate >= origin {
        (coordinate - origin) as f64
    } else {
        -((origin - coordinate) as f64)
    }
}

impl AngularConversion2 {
    pub(super) fn capture(&mut self, cell: usize, mass: f64, velocity: [f64; 2]) {
        // Preserve the ORIGINAL signed source, even if it survives and later
        // becomes a recipient or pays a neighboring source's negative debt.
        self.nodes.entry(cell).or_insert(Node {
            before: Snapshot { mass, velocity },
            parent: cell,
            transferred_mass: 0.0,
        });
    }

    fn root(&mut self, mut cell: usize) -> usize {
        while self.nodes[&cell].parent != cell {
            let parent = self.nodes[&cell].parent;
            let grandparent = self.nodes[&parent].parent;
            self.nodes.get_mut(&cell).expect("captured route").parent = grandparent;
            cell = parent;
        }
        cell
    }

    pub(super) fn link(
        &mut self,
        source: usize,
        recipient: usize,
        mass: f64,
        velocity: [f64; 2],
        transferred: f64,
    ) {
        self.capture(recipient, mass, velocity);
        for cell in [source, recipient] {
            self.nodes
                .get_mut(&cell)
                .expect("captured route")
                .transferred_mass += transferred.abs();
        }
        let a = self.root(source);
        let b = self.root(recipient);
        self.nodes
            .get_mut(&a.max(b))
            .expect("captured route")
            .parent = a.min(b);
    }

    pub(super) fn finish(
        mut self,
        surface: &mut FreeSurface,
        report: &mut ConversionMomentumReport2,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceMomentumError2> {
        let mut components = BTreeMap::<usize, Vec<usize>>::new();
        let mut identifiers = Vec::with_capacity(self.nodes.len());
        for (index, &cell) in self.nodes.keys().enumerate() {
            if index % 256 == 0 {
                poll(gate)?;
            }
            identifiers.push(cell);
        }
        for (index, cell) in identifiers.into_iter().enumerate() {
            if index % 256 == 0 {
                poll(gate)?;
            }
            components.entry(self.root(cell)).or_default().push(cell);
        }
        let mut angular_defect = 0.0;
        for (origin, members) in components {
            poll(gate)?;
            let origin = [origin % surface.grid.nx, origin / surface.grid.nx];
            let mut cells = Vec::with_capacity(members.len());
            let (mut mass, mut first, mut mass_scale, mut moment_scale) = (0.0, [0.0; 2], 0.0, 0.0);
            for (index, cell) in members.into_iter().enumerate() {
                if index % 256 == 0 {
                    poll(gate)?;
                }
                let node = &self.nodes[&cell];
                let after = snapshot(surface, cell)?;
                let position = [
                    offset(cell % surface.grid.nx, origin[0]),
                    offset(cell / surface.grid.nx, origin[1]),
                ];
                mass += after.mass;
                for d in 0..2 {
                    first[d] += after.mass * position[d];
                }
                let scale = node.before.mass.abs() + after.mass + node.transferred_mass;
                mass_scale += scale;
                moment_scale += scale * position[0].hypot(position[1]);
                cells.push(Participant {
                    cell,
                    before: node.before,
                    after,
                    position,
                });
            }
            // Same mass/linear roundoff policy as the original owner, scoped
            // to actual connected routes; angular units add local lever arms.
            let mass_roundoff = 128.0 * f64::EPSILON * mass_scale;
            let momentum_roundoff = 256.0 * f64::EPSILON * mass_scale;
            let angular_roundoff = 256.0 * f64::EPSILON * moment_scale;
            let (dm, dp, dl) = difference(&cells, gate)?;
            if !mass.is_finite()
                || !mass_roundoff.is_finite()
                || !momentum_roundoff.is_finite()
                || !angular_roundoff.is_finite()
                || !dm.is_finite()
                || dm.abs() > mass_roundoff
                || dp
                    .iter()
                    .any(|p| !p.is_finite() || p.abs() > momentum_roundoff)
                || !dl.is_finite()
            {
                return Err(invalid(
                    "angular conversion component mass or momentum is unrepresentable",
                ));
            }
            if mass > 0.0 {
                let centroid = first.map(|x| x / mass);
                let translation = dp.map(|p| -p / mass);
                let mut inertia = 0.0;
                for (index, cell) in cells.iter().enumerate() {
                    if index % 256 == 0 {
                        poll(gate)?;
                    }
                    let r = [
                        cell.position[0] - centroid[0],
                        cell.position[1] - centroid[1],
                    ];
                    inertia += cell.after.mass * r[0].mul_add(r[0], r[1] * r[1]);
                }
                // Translating by -dP/M changes orbital L by -c cross dP.
                // The remaining scalar spin therefore supplies c cross dP-dL.
                let required = cross(centroid, dp) - dl;
                let omega = if inertia > 0.0 {
                    required / inertia
                } else if inertia == 0.0 && required.abs() <= 0.5 * angular_roundoff {
                    0.0
                } else {
                    return Err(invalid(
                        "converted cells cannot resolve angular momentum; refine",
                    ));
                };
                if !inertia.is_finite() || !omega.is_finite() {
                    return Err(invalid(
                        "angular conversion inertia or spin is unrepresentable",
                    ));
                }
                for (index, cell) in cells.iter_mut().enumerate() {
                    if index % 256 == 0 {
                        poll(gate)?;
                    }
                    if cell.after.mass == 0.0 {
                        continue;
                    }
                    let r = [
                        cell.position[0] - centroid[0],
                        cell.position[1] - centroid[1],
                    ];
                    let previous = cell.after.velocity;
                    let next = [
                        previous[0] + translation[0] - omega * r[1],
                        previous[1] + translation[1] + omega * r[0],
                    ];
                    admit_velocity(next)?;
                    if next != previous {
                        let rho = surface.grid.f[cell.cell].iter().sum::<f64>();
                        let force = core::array::from_fn(|d| {
                            surface.grid.g[d].mul_add(rho, surface.grid.fext[cell.cell][d])
                        });
                        surface.grid.f[cell.cell] =
                            retarget_velocity(surface.grid.f[cell.cell], previous, next, force)?;
                    }
                    cell.after = snapshot(surface, cell.cell)?;
                    let du = [
                        cell.after.velocity[0] - previous[0],
                        cell.after.velocity[1] - previous[1],
                    ];
                    report.kinetic_energy_change += cell.after.mass
                        * (0..2)
                            .map(|d| (previous[d] + 0.5 * du[d]) * du[d])
                            .sum::<f64>();
                }
            }
            let (dm, dp, dl) = difference(&cells, gate)?;
            if !dm.is_finite()
                || dm.abs() > mass_roundoff
                || dp
                    .iter()
                    .any(|p| !p.is_finite() || p.abs() > momentum_roundoff)
                || !dl.is_finite()
                || dl.abs() > angular_roundoff
            {
                return Err(invalid(
                    "converted populations do not preserve both momenta; refine",
                ));
            }
            // Keep the original local mass/P defects; closure is additional.
            report.absolute_mass_defect += dm.abs();
            for d in 0..2 {
                report.absolute_momentum_defect[d] += dp[d].abs();
            }
            angular_defect += dl.abs();
        }
        poll(gate)?;
        report.absolute_angular_momentum_defect = Some(angular_defect);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
