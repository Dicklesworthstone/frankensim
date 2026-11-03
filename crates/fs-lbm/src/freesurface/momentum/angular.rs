//! Planar angular closure on the actual adjacent conversion graph.
//!
//! With final positive masses m_i, centroid c and I=sum(m_i*|x_i-c|^2),
//! the least mass-weighted velocity change is a translation plus rotation:
//! du_i=-dP/M + omega J(x_i-c), omega=(c cross dP-dL)/I.
//! These close linear and angular momentum in the declared lattice frame.
//! First mass moments are not conserved, so boost-equivariance is not claimed.

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
    a[0] * b[1] - a[1] * b[0]
}

fn snapshot(surface: &FreeSurface, cell: usize) -> Result<Snapshot, FreeSurfaceMomentumError2> {
    if surface.grid.flags[cell] == Cell::Gas {
        return Ok(Snapshot {
            mass: 0.0,
            velocity: [0.0; 2],
        });
    }
    let state = surface.grid.moments(cell);
    let mass = if surface.grid.flags[cell] == Cell::Fluid {
        state.rho
    } else {
        surface.mass[cell]
    };
    if !mass.is_finite() || mass < 0.0 || mass > state.rho {
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

// Subtract each cell's inventory before summing, keeping a local fixed origin
// and retaining the original signed debt of cells that disappear.
fn difference(cells: &[Participant]) -> (f64, [f64; 2], f64) {
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 2], 0.0);
    for cell in cells {
        let dm = cell.after.mass - cell.before.mass;
        let dp = core::array::from_fn(|d| {
            cell.after.mass.mul_add(
                cell.after.velocity[d] - cell.before.velocity[d],
                dm * cell.before.velocity[d],
            )
        });
        mass += dm;
        angular += cross(cell.position, dp);
        for d in 0..2 {
            momentum[d] += dp[d];
        }
    }
    (mass, momentum, angular)
}

fn coordinate_difference(a: usize, b: usize) -> f64 {
    if a >= b {
        (a - b) as f64
    } else {
        -((b - a) as f64)
    }
}

impl AngularConversion2 {
    pub(super) fn capture_sources(
        surface: &FreeSurface,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceMomentumError2> {
        let mut conversion = Self::default();
        for (cell, flag) in surface.grid.flags.iter().enumerate() {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            if *flag != Cell::Interface {
                continue;
            }
            let state = surface.grid.moments(cell);
            let mass = surface.mass[cell];
            if mass < 0.0 || mass > state.rho {
                // Capture BEFORE the class cascade. A positive source can be
                // made interface again by a receding neighbour, then receive
                // debt; its original signed inventory must survive both roles.
                conversion.capture(cell, mass, state.u);
            }
        }
        Ok(conversion)
    }

    fn capture(&mut self, cell: usize, mass: f64, velocity: [f64; 2]) {
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
                .expect("captured source or recipient")
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
        gate: &CancelGate,
        report: &mut ConversionMomentumReport2,
    ) -> Result<(), FreeSurfaceMomentumError2> {
        let mut components = BTreeMap::<usize, Vec<usize>>::new();
        for cell in self.nodes.keys().copied().collect::<Vec<_>>() {
            poll(gate)?;
            components.entry(self.root(cell)).or_default().push(cell);
        }
        let mut angular_defect = 0.0;
        for (origin, members) in components {
            let origin = [origin % surface.grid.nx, origin / surface.grid.nx];
            let mut cells = Vec::with_capacity(members.len());
            let (mut mass, mut first, mut mass_scale, mut moment_scale) = (0.0, [0.0; 2], 0.0, 0.0);
            for cell in members {
                poll(gate)?;
                let node = &self.nodes[&cell];
                let after = snapshot(surface, cell)?;
                let position = [
                    coordinate_difference(cell % surface.grid.nx, origin[0]),
                    coordinate_difference(cell / surface.grid.nx, origin[1]),
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
            let mass_roundoff = 128.0 * f64::EPSILON * mass_scale;
            let momentum_roundoff = 256.0 * f64::EPSILON * mass_scale;
            let angular_roundoff = 256.0 * f64::EPSILON * moment_scale;
            let (dm, dp, dl) = difference(&cells);
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
                    "angular component mass or momentum is unrepresentable",
                ));
            }
            if mass > 0.0 {
                let centroid = first.map(|x| x / mass);
                let translation = dp.map(|p| -p / mass);
                let mut inertia = 0.0;
                for cell in &cells {
                    poll(gate)?;
                    let r = core::array::from_fn(|d| cell.position[d] - centroid[d]);
                    inertia += cell.after.mass * norm2(r);
                }
                let required = cross(centroid, dp) - dl;
                // A single surviving point has zero inertia and cannot carry
                // spin independently of its linear momentum. Never invent an
                // epsilon inertia or borrow a distant cell to hide this loss.
                let omega = if required.is_finite() && required.abs() <= 0.5 * angular_roundoff {
                    0.0
                } else if required.is_finite() && inertia.is_finite() && inertia > 0.0 {
                    required / inertia
                } else {
                    return Err(invalid(
                        "converted cells cannot resolve angular momentum; refine",
                    ));
                };
                if !omega.is_finite() || !centroid.iter().chain(&translation).all(|x| x.is_finite())
                {
                    return Err(invalid("angular correction is unrepresentable; refine"));
                }
                for cell in &mut cells {
                    poll(gate)?;
                    if cell.after.mass == 0.0 {
                        continue;
                    }
                    let r: [f64; 2] = core::array::from_fn(|d| cell.position[d] - centroid[d]);
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
                    let du: [f64; 2] =
                        core::array::from_fn(|d| cell.after.velocity[d] - previous[d]);
                    report.kinetic_energy_change += cell.after.mass
                        * (0..2)
                            .map(|d| (previous[d] + 0.5 * du[d]) * du[d])
                            .sum::<f64>();
                }
            }
            let (dm, dp, dl) = difference(&cells);
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
            report.absolute_mass_defect += dm.abs();
            for d in 0..2 {
                report.absolute_momentum_defect[d] += dp[d].abs();
            }
            angular_defect += dl.abs();
        }
        report.absolute_angular_momentum_defect = Some(angular_defect);
        Ok(())
    }
}
