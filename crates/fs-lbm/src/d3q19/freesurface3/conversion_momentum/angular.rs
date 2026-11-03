//! Angular closure on actual adjacent conversion routes, never detached fronts.
//!
//! This is a numerical remap in the declared lattice frame. Moving mass between
//! cell centres changes its first moment; even uniform boosted flow may need a
//! velocity correction to retain angular momentum. No Galilean-invariant
//! geometric projection or whole-timestep angular conservation is claimed.
//! FreeSurface3 has wall-bounded coordinates, not periodic angular transport.

use super::*;
use crate::angular_momentum::angular_correction;

type CellId = (u64, usize);

#[derive(Clone, Copy)]
struct Snapshot {
    mass: f64,
    velocity: [f64; 3],
}

struct Node {
    before: Snapshot,
    parent: CellId,
    transferred_mass: f64,
}

#[derive(Default)]
pub(super) struct AngularConversion3 {
    nodes: BTreeMap<CellId, Node>,
}

struct Participant {
    cell: CellId,
    before: Snapshot,
    after: Snapshot,
    position: [f64; 3],
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}

fn snapshot(surface: &FreeSurface3, cell: CellId) -> Result<Snapshot, FreeSurfaceError3> {
    let (key, lane) = cell;
    let tile = &surface.tiles[&key];
    if tile.cells[lane] == Cell3::Gas {
        return Ok(Snapshot {
            mass: 0.0,
            velocity: [0.0; 3],
        });
    }
    let (density, velocity) = state(surface, key, lane)?;
    let mass = if tile.cells[lane] == Cell3::Fluid {
        density
    } else {
        tile.mass[lane]
    };
    if !mass.is_finite() || mass < 0.0 || mass > density {
        return Err(refused(
            "angular conversion needs bounded final cell masses",
        ));
    }
    Ok(Snapshot { mass, velocity })
}

// Difference each cell's actual inventory before summing, rather than
// subtracting large global P or L totals. The same fixed local origin is used
// before and after conversion, including sources whose signed mass disappears.
fn difference(cells: &[Participant]) -> (f64, [f64; 3], [f64; 3]) {
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for cell in cells {
        let dm = cell.after.mass - cell.before.mass;
        let dp = core::array::from_fn(|d| {
            cell.after.mass.mul_add(
                cell.after.velocity[d] - cell.before.velocity[d],
                dm * cell.before.velocity[d],
            )
        });
        let dl = cross(cell.position, dp);
        mass += dm;
        for d in 0..3 {
            momentum[d] += dp[d];
            angular[d] += dl[d];
        }
    }
    (mass, momentum, angular)
}

impl AngularConversion3 {
    pub(super) fn capture(&mut self, cell: CellId, mass: f64, velocity: [f64; 3]) {
        // A source may become a recipient later. Keep its ORIGINAL signed
        // inventory, captured before any class change or mass redistribution.
        self.nodes.entry(cell).or_insert(Node {
            before: Snapshot { mass, velocity },
            parent: cell,
            transferred_mass: 0.0,
        });
    }

    fn root(&mut self, mut cell: CellId) -> CellId {
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
        source: CellId,
        recipient: CellId,
        mass: f64,
        velocity: [f64; 3],
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
        surface: &mut FreeSurface3,
        report: &mut ConversionMomentumReport3,
    ) -> Result<(), FreeSurfaceError3> {
        let mut components = BTreeMap::<CellId, Vec<CellId>>::new();
        for cell in self.nodes.keys().copied().collect::<Vec<_>>() {
            components.entry(self.root(cell)).or_default().push(cell);
        }
        let mut angular_defect = [0.0; 3];
        for (origin, members) in components {
            let (ox, oy, oz) = FreeSurface3::global_of(origin.0, origin.1);
            let origin = [ox, oy, oz];
            let mut cells = Vec::with_capacity(members.len());
            let (mut mass, mut first, mut mass_scale, mut moment_scale) = (0.0, [0.0; 3], 0.0, 0.0);
            for cell in members {
                let node = &self.nodes[&cell];
                let after = snapshot(surface, cell)?;
                let (x, y, z) = FreeSurface3::global_of(cell.0, cell.1);
                // Integer differences keep the cell-centre lever arms local.
                let position = core::array::from_fn(|d| ([x, y, z][d] - origin[d]) as f64);
                mass += after.mass;
                for d in 0..3 {
                    first[d] += after.mass * position[d];
                }
                let scale = node.before.mass.abs() + after.mass + node.transferred_mass;
                mass_scale += scale;
                moment_scale += scale * norm(position);
                cells.push(Participant {
                    cell,
                    before: node.before,
                    after,
                    position,
                });
            }
            // Original conversion roundoff policy, summed only over this
            // route-connected set. Angular units add the actual local lever arm.
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
                || dl.iter().any(|l| !l.is_finite())
            {
                return Err(refused(
                    "angular conversion component mass or momentum is unrepresentable",
                ));
            }
            if mass > 0.0 {
                let centroid = first.map(|x| x / mass);
                let translation = dp.map(|p| -p / mass);
                let mut inertia = [0.0; 9];
                for cell in &cells {
                    let r: [f64; 3] = core::array::from_fn(|d| cell.position[d] - centroid[d]);
                    for a in 0..3 {
                        for b in 0..3 {
                            let value = if a == b {
                                (0..3).filter(|&d| d != a).map(|d| r[d] * r[d]).sum()
                            } else {
                                -r[a] * r[b]
                            };
                            inertia[3 * a + b] += cell.after.mass * value;
                        }
                    }
                }
                let orbital = cross(centroid, dp);
                let required = core::array::from_fn(|d| orbital[d] - dl[d]);
                // Reserve half the same allowance for native population readback.
                let omega = angular_correction(&inertia, required, 0.5 * angular_roundoff)
                    .ok_or_else(|| {
                        refused("converted cells cannot resolve angular momentum; refine")
                    })?;
                for cell in &mut cells {
                    if cell.after.mass == 0.0 {
                        continue;
                    }
                    let r = core::array::from_fn(|d| cell.position[d] - centroid[d]);
                    let rotation = cross(omega, r);
                    let previous = cell.after.velocity;
                    let next = core::array::from_fn(|d| previous[d] + translation[d] + rotation[d]);
                    admit_velocity(next)?;
                    if next != previous {
                        let (key, lane) = cell.cell;
                        let slot = surface.grid.slot_of(key).expect("active angular recipient");
                        let populations = forcing::retarget_velocity(
                            surface.grid.populations(slot, lane),
                            previous,
                            next,
                            surface.force_at_lane(key, lane),
                        )?;
                        surface.grid.set_populations(slot, lane, populations);
                    }
                    cell.after = snapshot(surface, cell.cell)?;
                    let du: [f64; 3] =
                        core::array::from_fn(|d| cell.after.velocity[d] - previous[d]);
                    report.kinetic_energy_change += cell.after.mass
                        * (0..3)
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
                || !norm(dl).is_finite()
                || norm(dl) > angular_roundoff
            {
                return Err(refused(
                    "converted populations do not preserve both momenta; refine",
                ));
            }
            report.absolute_mass_defect += dm.abs();
            for d in 0..3 {
                report.absolute_momentum_defect[d] += dp[d].abs();
                angular_defect[d] += dl[d].abs();
            }
        }
        report.absolute_angular_momentum_defect = Some(angular_defect);
        Ok(())
    }
}
