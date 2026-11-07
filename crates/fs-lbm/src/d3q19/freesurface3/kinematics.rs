//! Resolved liquid velocity gradients on the original sparse VOF owner.
//!
//! Public observations use only occupied liquid cell centers. Gas, empty interface
//! cells, solid/inactive tiles and the current piston slab are missing samples,
//! never stationary fluid. One-sided derivatives use contiguous liquid centers
//! and expose their order; an unresolvable axis refuses instead of inventing a
//! zero derivative. The native wall/SGS transaction additionally supplies its
//! explicit halfway face velocity; a free-surface velocity is never invented.

use super::{CS2_3, CancelGate, Cell3, FreeSurface3, PlanarWallVelocityField3, TILE_CELLS};
use crate::MACH_LIMIT;
use crate::kinematics::{FlowKinematics3, KinematicsError, checkpoint};

/// Actual spatial stencil for one column of the liquid velocity gradient.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WetGradientStencil3 {
    CenteredSecondOrder,
    ForwardSecondOrder,
    BackwardSecondOrder,
    ForwardFirstOrder,
    BackwardFirstOrder,
    /// Quadratic reconstruction through the prescribed lower halfway wall
    /// and two occupied liquid centers; neither center is a gas ghost.
    LowerWallSecondOrder,
    /// The same reconstruction at an upper domain wall.
    UpperWallSecondOrder,
}

/// Resolved kinematics in lattice units, with per-axis reconstruction order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceKinematicsSample3 {
    pub kinematics: FlowKinematics3,
    /// x, y, z derivative columns, respectively.
    pub stencils: [WetGradientStencil3; 3],
}

#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceKinematicsError3 {
    Kinematics(KinematicsError),
    /// Neither adjacent liquid cell exists along this spatial axis. In
    /// particular, an unresolved one-cell sheet has no inferred normal curl.
    MissingDerivative {
        cell: [usize; 3],
        axis: usize,
    },
}

impl From<KinematicsError> for FreeSurfaceKinematicsError3 {
    fn from(value: KinematicsError) -> Self {
        Self::Kinematics(value)
    }
}

impl core::fmt::Display for FreeSurfaceKinematicsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-FREE-SURFACE-KINEMATICS3: {self:?}")
    }
}
impl std::error::Error for FreeSurfaceKinematicsError3 {}

/// Reusable sparse velocity cache with a fixed active-tile capacity. Memory
/// scales with active tiles, not the bounding domain, and queries allocate
/// nothing. Newly activated tiles exceeding capacity refuse before rebuilding.
pub struct FreeSurfaceKinematicsWorkspace3 {
    velocity: Vec<Option<[f64; 3]>>,
}

impl FreeSurfaceKinematicsWorkspace3 {
    /// Numerical storage only; excludes allocator bookkeeping. Gas lanes in
    /// active tiles count toward the budget, just as they do in the owner.
    pub fn required_bytes(maximum_active_tiles: usize) -> Result<usize, KinematicsError> {
        maximum_active_tiles
            .checked_mul(TILE_CELLS)
            .and_then(|n| n.checked_mul(core::mem::size_of::<Option<[f64; 3]>>()))
            .ok_or(KinematicsError::Input(
                "sparse kinematics workspace size overflow",
            ))
    }

    /// Reserve once for the chosen active-tile limit. A caller may reserve
    /// more than the current count to accommodate future liquid advance.
    pub fn new(maximum_active_tiles: usize, maximum_bytes: usize) -> Result<Self, KinematicsError> {
        let requested = Self::required_bytes(maximum_active_tiles)?;
        if requested > maximum_bytes {
            return Err(KinematicsError::Budget {
                requested,
                maximum: maximum_bytes,
            });
        }
        let cells = maximum_active_tiles * TILE_CELLS;
        let mut velocity = Vec::new();
        velocity
            .try_reserve_exact(cells)
            .map_err(|_| KinematicsError::Allocation)?;
        velocity.resize(cells, None);
        Ok(Self { velocity })
    }
}

/// Successful read-only snapshot of the owner's current half-force velocity.
/// Borrowing both owner and workspace prevents reuse across geometry changes
/// or sparse slot reshuffles. A new preparation always rebuilds all active lanes.
pub struct FreeSurfaceKinematics3<'a> {
    owner: &'a FreeSurface3,
    gate: &'a CancelGate,
    velocity: &'a [Option<[f64; 3]>],
    // Only the native wall/SGS transaction supplies this boundary epoch.
    // Public liquid-only observations retain their original stencil semantics.
    wall: Option<PlanarWallVelocityField3<'a>>,
}

impl FreeSurface3 {
    /// Cache the native physical velocity once per occupied liquid cell. No
    /// collision, SGS, boundary rule or population is changed. A cancelled or
    /// invalid preparation exposes no partial field and is safe to retry.
    pub fn kinematics<'a>(
        &'a self,
        gate: &'a CancelGate,
        workspace: &'a mut FreeSurfaceKinematicsWorkspace3,
    ) -> Result<FreeSurfaceKinematics3<'a>, FreeSurfaceKinematicsError3> {
        self.kinematics_with_wall(gate, workspace, None)
    }

    pub(crate) fn kinematics_with_wall<'a>(
        &'a self,
        gate: &'a CancelGate,
        workspace: &'a mut FreeSurfaceKinematicsWorkspace3,
        wall: Option<PlanarWallVelocityField3<'a>>,
    ) -> Result<FreeSurfaceKinematics3<'a>, FreeSurfaceKinematicsError3> {
        checkpoint(gate)?;
        if let Some(wall) = wall {
            let cells = self.domain_cells();
            let axis = wall.wall.axis();
            let faces = cells[(axis + 1) % 3].checked_mul(cells[(axis + 2) % 3]);
            if faces != Some(wall.tangential_velocity.len()) {
                return Err(
                    KinematicsError::Input("wall gradient face-field shape differs").into(),
                );
            }
            for (i, tangent) in wall.tangential_velocity.iter().enumerate() {
                if i % 64 == 0 {
                    checkpoint(gate)?;
                }
                let speed = tangent[0].hypot(tangent[1]);
                if !speed.is_finite() || speed / CS2_3.sqrt() >= MACH_LIMIT {
                    return Err(KinematicsError::Input(
                        "wall gradient velocity exceeds the native Mach domain",
                    )
                    .into());
                }
            }
        }
        let active_tiles = self.grid.active_keys().len();
        if active_tiles > workspace.velocity.len() / TILE_CELLS {
            return Err(KinematicsError::Budget {
                requested: FreeSurfaceKinematicsWorkspace3::required_bytes(active_tiles)?,
                maximum: workspace.velocity.len() * core::mem::size_of::<Option<[f64; 3]>>(),
            }
            .into());
        }
        for (slot, &key) in self.grid.active_keys().iter().enumerate() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let index = slot * TILE_CELLS + lane;
                workspace.velocity[index] = None;
                let (x, y, z) = Self::global_of(key, lane);
                if self.tiles[&key].cells[lane] == Cell3::Gas
                    || self.planar_piston_contains([x, y, z])
                {
                    continue;
                }
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                let force = self.force_at_lane(key, lane);
                let velocity = core::array::from_fn(|a| raw[a] + 0.5 * force[a] / rho);
                let fill = self.fill(key, lane);
                if !(rho.is_finite() && rho > 0.0 && fill.is_finite())
                    || velocity.iter().any(|v| !v.is_finite())
                {
                    return Err(KinematicsError::InvalidState {
                        cell: [x as usize, y as usize, z as usize],
                    }
                    .into());
                }
                if fill > 0.0 {
                    workspace.velocity[index] = Some(velocity);
                }
            }
        }
        checkpoint(gate)?;
        Ok(FreeSurfaceKinematics3 {
            owner: self,
            gate,
            velocity: &workspace.velocity[..active_tiles * TILE_CELLS],
            wall,
        })
    }
}

impl FreeSurfaceKinematics3<'_> {
    /// Liquid-cell kinematics. Gas, empty interfaces, solid/piston and
    /// out-of-domain positions return None. Missing liquid neighbors on both
    /// sides of an axis return an explicit spatial-resolution refusal.
    pub fn sample_at(
        &self,
        x: usize,
        y: usize,
        z: usize,
    ) -> Result<Option<FreeSurfaceKinematicsSample3>, FreeSurfaceKinematicsError3> {
        checkpoint(self.gate)?;
        let cell = [x, y, z];
        let (Ok(x), Ok(y), Ok(z)) = (i64::try_from(x), i64::try_from(y), i64::try_from(z)) else {
            return Ok(None);
        };
        let position = [x, y, z];
        let Some(center) = self.velocity_at(position) else {
            return Ok(None);
        };
        let mut gradient = [[0.0; 3]; 3];
        let mut stencils = [WetGradientStencil3::CenteredSecondOrder; 3];
        for axis in 0..3 {
            let mut plus = position;
            let mut minus = position;
            plus[axis] += 1;
            minus[axis] -= 1;
            let left = self.velocity_at(minus);
            let right = self.velocity_at(plus);
            // A prescribed wall value is located halfway from the adjacent
            // liquid center, not at the absent solid/gas lattice node. The
            // quadratic fit therefore has unequal spacings (1/2, 1).
            if let Some(wall) = self.wall_velocity_at(minus, axis) {
                let Some(inward) = right else {
                    return Err(FreeSurfaceKinematicsError3::MissingDerivative { cell, axis });
                };
                for component in 0..3 {
                    gradient[component][axis] = (center[component] - wall[component])
                        + (inward[component] - wall[component]) / 3.0;
                }
                stencils[axis] = WetGradientStencil3::LowerWallSecondOrder;
                continue;
            }
            if let Some(wall) = self.wall_velocity_at(plus, axis) {
                let Some(inward) = left else {
                    return Err(FreeSurfaceKinematicsError3::MissingDerivative { cell, axis });
                };
                for component in 0..3 {
                    gradient[component][axis] = -(center[component] - wall[component])
                        - (inward[component] - wall[component]) / 3.0;
                }
                stencils[axis] = WetGradientStencil3::UpperWallSecondOrder;
                continue;
            }
            let (column, stencil): ([f64; 3], WetGradientStencil3) = match (left, right) {
                (Some(left), Some(right)) => (
                    core::array::from_fn(|a| 0.5 * (right[a] - left[a])),
                    WetGradientStencil3::CenteredSecondOrder,
                ),
                (None, Some(first)) | (Some(first), None) => {
                    let forward = right.is_some();
                    let sign = if forward { 1.0 } else { -1.0 };
                    let mut next = position;
                    next[axis] += if forward { 2 } else { -2 };
                    if let Some(second) = self.velocity_at(next) {
                        (
                            core::array::from_fn(|a| {
                                sign * (2.0 * (first[a] - center[a])
                                    - 0.5 * (second[a] - center[a]))
                            }),
                            if forward {
                                WetGradientStencil3::ForwardSecondOrder
                            } else {
                                WetGradientStencil3::BackwardSecondOrder
                            },
                        )
                    } else if let Some(wall) = self.wall_velocity_at(next, axis) {
                        // The last center in a two-layer liquid film still
                        // has three physical samples: itself, the first
                        // center one cell away, and the wall 3/2 cells away.
                        // This supplies a second-order one-sided derivative
                        // without using the missing free-surface gas velocity.
                        (
                            core::array::from_fn(|a| {
                                sign * (3.0 * (first[a] - center[a])
                                    - (4.0 / 3.0) * (wall[a] - center[a]))
                            }),
                            if forward {
                                WetGradientStencil3::UpperWallSecondOrder
                            } else {
                                WetGradientStencil3::LowerWallSecondOrder
                            },
                        )
                    } else {
                        (
                            core::array::from_fn(|a| sign * (first[a] - center[a])),
                            if forward {
                                WetGradientStencil3::ForwardFirstOrder
                            } else {
                                WetGradientStencil3::BackwardFirstOrder
                            },
                        )
                    }
                }
                (None, None) => {
                    return Err(FreeSurfaceKinematicsError3::MissingDerivative { cell, axis });
                }
            };
            for component in 0..3 {
                gradient[component][axis] = column[component];
            }
            stencils[axis] = stencil;
        }
        Ok(Some(FreeSurfaceKinematicsSample3 {
            kinematics: FlowKinematics3::from_gradient(gradient)?,
            stencils,
        }))
    }

    /// Visit occupied liquid cells in deterministic Morton/lane order. Local
    /// stencil refusals are delivered to the visitor, allowing resolved cells
    /// elsewhere to remain useful. Cancellation stops the sweep, including
    /// across long empty stretches; the visitor should also keep work bounded.
    pub fn visit_samples(
        &self,
        mut visit: impl FnMut(
            [usize; 3],
            Result<FreeSurfaceKinematicsSample3, FreeSurfaceKinematicsError3>,
        ),
    ) -> Result<(), FreeSurfaceKinematicsError3> {
        for (slot, &key) in self.owner.grid.active_keys().iter().enumerate() {
            checkpoint(self.gate)?;
            for lane in 0..TILE_CELLS {
                if self.velocity[slot * TILE_CELLS + lane].is_none() {
                    continue;
                }
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let cell = [x as usize, y as usize, z as usize];
                let sample = self
                    .sample_at(cell[0], cell[1], cell[2])
                    .map(|sample| sample.expect("the borrowed occupied-liquid cache is immutable"));
                checkpoint(self.gate)?;
                visit(cell, sample);
            }
        }
        checkpoint(self.gate)?;
        Ok(())
    }

    fn velocity_at(&self, position: [i64; 3]) -> Option<[f64; 3]> {
        let (slot, lane) = self
            .owner
            .grid
            .resolve_source(position[0], position[1], position[2])?;
        self.velocity[slot * TILE_CELLS + lane]
    }

    // Recognize ONLY the explicitly prescribed domain plane. Inactive tiles,
    // piston slabs, other walls and gas cells retain the liquid-only fallback.
    // The returned velocity belongs to the Cartesian face center used by the
    // original moving-wall reflection, with exactly zero normal component.
    fn wall_velocity_at(&self, position: [i64; 3], axis: usize) -> Option<[f64; 3]> {
        let wall = self.wall?;
        if axis != wall.wall.axis() {
            return None;
        }
        let cells = self.owner.domain_cells();
        let boundary = if wall.wall.upper() {
            cells[axis] as i64
        } else {
            -1
        };
        if position[axis] != boundary {
            return None;
        }
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let ia = usize::try_from(position[a]).ok()?;
        let ib = usize::try_from(position[b]).ok()?;
        if ia >= cells[a] || ib >= cells[b] {
            return None;
        }
        let tangent = wall.tangential_velocity[ib * cells[a] + ia];
        let mut velocity = [0.0; 3];
        velocity[a] = tangent[0];
        velocity[b] = tangent[1];
        Some(velocity)
    }
}
