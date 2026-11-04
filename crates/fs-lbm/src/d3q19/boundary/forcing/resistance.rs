//! Implicit linear porous resistance on the native D3Q19 momentum owner.

use super::*;

/// Prescribed local law `F = body_force + D (skeleton_velocity - u)`.
///
/// All values use lattice units, `dx = dt = 1`. Principal resistances are
/// force density per relative velocity, NOT accelerations or permeabilities.
/// For Darcy/Brinkman drag they are the scaled dynamic-viscosity/permeability
/// ratios. `axes[k]` is principal axis k resolved in grid coordinates, so
/// `D = sum_k principal[k] axes[k] axes[k]^T`. Axes must be orthonormal;
/// finite nonnegative principal values admit anisotropy and free directions.
/// Infinite resistance/sealed directions require an explicit constraint.
///
/// Sample prescribed coefficients, skeleton motion and body force at the
/// step midpoint and hold them for this lattice tick. The skeleton has no
/// evolved inertia here; its equal-and-opposite reaction is returned to the
/// caller. This law replaces BOTH previously held force components.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearResistance3 {
    /// Finite nonnegative drag coefficients in the three principal directions.
    pub principal: [f64; 3],
    /// Unit, mutually orthogonal material axes resolved in grid coordinates.
    pub axes: [[f64; 3]; 3],
    /// Prescribed velocity of the porous skeleton during this step.
    pub skeleton_velocity: [f64; 3],
    /// Additional volume force density, independent of fluid velocity.
    pub body_force: [f64; 3],
}

impl Default for LinearResistance3 {
    fn default() -> Self {
        Self {
            principal: [0.0; 3],
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            skeleton_velocity: [0.0; 3],
            body_force: [0.0; 3],
        }
    }
}

/// One tick's external forcing and reaction on the prescribed skeleton.
/// Angular impulses use cell centres and the requested moment origin.
///
/// In a closed domain with ALL wall links selected, fluid momentum changes,
/// to floating-point precision, by `external_impulse - skeleton_impulse +
/// walls.fluid_population_impulse`. Omitted wall links require their
/// independently measured contribution. Angular momentum uses the corresponding
/// angular terms. Periodic cuts additionally carry angular momentum, as for
/// the native unforced stream.
///
/// Work uses the arithmetic mean of each cell's initial and final physical
/// velocity. `fluid_drag_work + skeleton_work = -relative_dissipation`
/// up to rounding. This is the local force-work identity, NOT the total
/// transported fluid kinetic-energy change: density, pressure, advection and
/// collision remain owned by the original lattice operator. No heat is
/// deposited, and native wall work retains its separate link convention.
/// The skeleton impulse/work is the constitutive midpoint transfer, retained
/// separately even when a larger body force rounds a tiny drag out of the
/// total lattice force. Sub-ULP impulses are not guaranteed to accumulate in
/// the executed fluid momentum; this receipt is not a compensated integrator.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct ResistanceMomentumExchange3 {
    /// Impulse supplied by the independent body force.
    pub external_impulse: [f64; 3],
    /// Moment of the independent body impulse about the requested origin.
    pub external_angular_impulse: [f64; 3],
    /// Opposite of the fluid's trapezoidal drag impulse.
    pub skeleton_impulse: [f64; 3],
    /// Moment of the skeleton reaction about the requested origin.
    pub skeleton_angular_impulse: [f64; 3],
    /// Independent body impulse dotted with midpoint fluid velocity.
    pub external_work: f64,
    /// Fluid drag impulse dotted with midpoint fluid velocity.
    pub fluid_drag_work: f64,
    /// Skeleton impulse dotted with prescribed skeleton velocity.
    pub skeleton_work: f64,
    /// Nonnegative `sum (u_mid-v_s)^T D (u_mid-v_s)` for the admitted law.
    pub relative_dissipation: f64,
    /// Native exchange on the selected wall links, from this same collision.
    pub walls: WallMomentumExchange3,
}

impl LinearResistance3 {
    fn admit(self) -> Result<(), BoundaryStepError3> {
        let invalid = |reason| BoundaryStepError3::EndpointForce { reason };
        if !self.principal.iter().all(|x| x.is_finite() && *x >= 0.0)
            || !self
                .axes
                .iter()
                .flatten()
                .chain(&self.skeleton_velocity)
                .chain(&self.body_force)
                .all(|x| x.is_finite())
        {
            return Err(invalid(
                "linear resistance, axes, skeleton velocity and body force must be finite; principal resistance must be nonnegative",
            ));
        }
        for a in 0..3 {
            for b in 0..=a {
                if (dot(self.axes[a], self.axes[b]) - f64::from(a == b)).abs() > 64.0 * f64::EPSILON
                {
                    return Err(invalid(
                        "linear resistance requires orthonormal material axes",
                    ));
                }
            }
        }
        Ok(())
    }

    fn drag(self, velocity: [f64; 3]) -> [f64; 3] {
        if self.principal == [0.0; 3] {
            return [0.0; 3];
        }
        let relative = core::array::from_fn(|a| self.skeleton_velocity[a] - velocity[a]);
        if self
            .principal
            .iter()
            .all(|&value| value == self.principal[0])
        {
            return relative.map(|v| self.principal[0] * v);
        }
        let principal: [f64; 3] =
            core::array::from_fn(|a| self.principal[a] * dot(self.axes[a], relative));
        core::array::from_fn(|a| (0..3).map(|k| self.axes[k][a] * principal[k]).sum())
    }

    fn force(self, velocity: [f64; 3]) -> [f64; 3] {
        // Preserve the native held-force reduction exactly when drag is zero.
        if self.principal == [0.0; 3] {
            self.body_force
        } else {
            let drag = self.drag(velocity);
            core::array::from_fn(|a| self.body_force[a] + drag[a])
        }
    }

    fn endpoint_force(self, rho: f64, momentum: [f64; 3]) -> Option<([f64; 3], [f64; 3])> {
        if self.principal == [0.0; 3] {
            return Some((self.body_force, [0.0; 3]));
        }
        // (rho I + D/2) (u-v_s) = m_raw + b/2 - rho v_s.
        // Solve in material axes: three divisions, no iteration/factorization.
        let rhs = core::array::from_fn(|a| {
            momentum[a] + 0.5 * self.body_force[a] - rho * self.skeleton_velocity[a]
        });
        let denominator = self.principal.map(|d| rho + 0.5 * d);
        if !rhs.iter().chain(&denominator).all(|v| v.is_finite()) {
            return None;
        }
        let drag = if self
            .principal
            .iter()
            .all(|&value| value == self.principal[0])
        {
            rhs.map(|r| -self.principal[0] * (r / denominator[0]))
        } else {
            let principal: [f64; 3] = core::array::from_fn(|k| {
                -self.principal[k] * (dot(self.axes[k], rhs) / denominator[k])
            });
            core::array::from_fn(|a| (0..3).map(|k| self.axes[k][a] * principal[k]).sum())
        };
        let force = core::array::from_fn(|a| self.body_force[a] + drag[a]);
        force
            .iter()
            .chain(&drag)
            .all(|v| v.is_finite())
            .then_some((force, drag))
    }
}

impl ResistanceMomentumExchange3 {
    fn add(&mut self, other: Self) {
        for a in 0..3 {
            self.external_impulse[a] += other.external_impulse[a];
            self.external_angular_impulse[a] += other.external_angular_impulse[a];
            self.skeleton_impulse[a] += other.skeleton_impulse[a];
            self.skeleton_angular_impulse[a] += other.skeleton_angular_impulse[a];
        }
        self.external_work += other.external_work;
        self.fluid_drag_work += other.fluid_drag_work;
        self.skeleton_work += other.skeleton_work;
        self.relative_dissipation += other.relative_dissipation;
    }

    fn is_finite(self) -> bool {
        [
            self.external_impulse,
            self.external_angular_impulse,
            self.skeleton_impulse,
            self.skeleton_angular_impulse,
        ]
        .into_iter()
        .flatten()
        .chain([
            self.external_work,
            self.fluid_drag_work,
            self.skeleton_work,
            self.relative_dissipation,
        ])
        .all(f64::is_finite)
    }
}

#[derive(Clone, Copy, Default)]
pub(in crate::d3q19::boundary) struct PreparedResistance3 {
    law: LinearResistance3,
    velocity: [f64; 3],
}

impl PreparedResistance3 {
    pub(in crate::d3q19::boundary) fn prepare(
        self,
        populations: &mut [f64; Q3],
        previous: [f64; 3],
    ) -> [f64; 3] {
        let force = self.law.force(self.velocity);
        if force != previous {
            crate::thermal::forcing::reframe3(populations, self.velocity, previous, force);
        }
        force
    }
}

pub(in crate::d3q19::boundary) struct ResistanceGroup3<'a> {
    force: &'a mut [[f64; 3]],
    exchange: ResistanceMomentumExchange3,
}

pub(in crate::d3q19::boundary) struct ResistanceStream3<'a> {
    prepared: &'a [PreparedResistance3],
    pub(in crate::d3q19::boundary) groups: Vec<Mutex<ResistanceGroup3<'a>>>,
    origin: [f64; 3],
}

impl ResistanceStream3<'_> {
    pub(in crate::d3q19::boundary) fn finish_cell(
        &self,
        populations: &[f64; Q3],
        tile: usize,
        lane: usize,
        local_tile: usize,
        cell: [usize; 3],
        output: &mut ResistanceGroup3<'_>,
    ) -> bool {
        let Some((rho, momentum)) = raw_moments(populations) else {
            return false;
        };
        let prepared = self.prepared[tile * TILE_CELLS + lane];
        let law = prepared.law;
        let Some((force, endpoint_drag)) = law.endpoint_force(rho, momentum) else {
            return false;
        };
        let velocity: [f64; 3] = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
        let midpoint = core::array::from_fn(|a| 0.5 * prepared.velocity[a] + 0.5 * velocity[a]);
        let initial_drag = law.drag(prepared.velocity);
        // Preserve a small physical reaction even when the independent body
        // force is large enough to round it away in the total force vector.
        let reaction = core::array::from_fn(|a| -0.5 * initial_drag[a] - 0.5 * endpoint_drag[a]);
        let offset = core::array::from_fn(|a| cell[a] as f64 + 0.5 - self.origin[a]);
        let relative = core::array::from_fn(|a| midpoint[a] - law.skeleton_velocity[a]);
        let relative_dissipation = if law.principal.iter().all(|&d| d == law.principal[0]) {
            law.principal[0] * dot(relative, relative)
        } else {
            (0..3)
                .map(|a| {
                    let r = dot(law.axes[a], relative);
                    (law.principal[a] * r) * r
                })
                .sum()
        };
        let exchange = ResistanceMomentumExchange3 {
            external_impulse: law.body_force,
            external_angular_impulse: cross(offset, law.body_force),
            skeleton_impulse: reaction,
            skeleton_angular_impulse: cross(offset, reaction),
            external_work: dot(law.body_force, midpoint),
            fluid_drag_work: -dot(reaction, midpoint),
            skeleton_work: dot(reaction, law.skeleton_velocity),
            relative_dissipation,
            walls: WallMomentumExchange3::default(),
        };
        if !velocity.iter().all(|v| v.is_finite()) || !exchange.is_finite() {
            return false;
        }
        output.force[local_tile * TILE_CELLS + lane] = force;
        output.exchange.add(exchange);
        true
    }
}

impl BoundaryGrid3 {
    /// Advance one closed/periodic tick with implicit anisotropic resistance.
    ///
    /// Samples fluid cells once in x-fastest order, skipping static solids.
    /// The old state is reframed at fixed physical velocity/stress, then the
    /// native collision uses `F_n = b + D(v_s-u_n)`. After native streaming,
    /// the endpoint half-source relation is solved analytically for `u_(n+1)`
    /// and `F_(n+1)`. Forces and populations publish together; an invalid law,
    /// cancellation or a late refusal leaves both accepted fields unchanged.
    /// Open faces are refused: their reconstruction needs a joint source solve.
    ///
    /// For an isolated uniform source this is implicit midpoint/trapezoidal
    /// drag, with bounded relative kinetic energy for any admitted resistance.
    /// Large resistance can reverse velocity each tick; this method is not
    /// L-stable or exact exponential relaxation. Resolve the relaxation time
    /// for accuracy, and retain the original low-Mach/collision requirements.
    ///
    /// Uses the existing population scratch and grouped collision/pull kernels.
    /// One coefficient/initial-velocity vector and one endpoint-force vector
    /// are sampled; the solve and receipts are fused into streaming. Moments
    /// in this overload use the zero origin; no wall links are measured.
    pub fn step_with_linear_resistance(
        &mut self,
        resistance: impl FnMut([usize; 3]) -> LinearResistance3,
        gate: &CancelGate,
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_linear_resistance_inner::<fs_exec::TilePool>(
            None,
            gate,
            resistance,
            &[],
            [0.0; 3],
        )
    }

    /// Same source step with native selected-wall impulse, torque and work.
    /// Both skeleton and wall moments use `moment_origin`.
    pub fn step_with_linear_resistance_and_wall_momentum(
        &mut self,
        resistance: impl FnMut([usize; 3]) -> LinearResistance3,
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_linear_resistance_inner::<fs_exec::TilePool>(
            None,
            gate,
            resistance,
            links,
            moment_origin,
        )
    }

    /// Pooled counterpart with identical cell arithmetic and ordered receipts.
    pub fn step_pooled_with_linear_resistance<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        resistance: impl FnMut([usize; 3]) -> LinearResistance3,
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_linear_resistance_inner(Some(pool), gate, resistance, &[], [0.0; 3])
    }

    /// Pooled source step returning the same selected native wall exchange.
    pub fn step_pooled_with_linear_resistance_and_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        resistance: impl FnMut([usize; 3]) -> LinearResistance3,
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_linear_resistance_inner(Some(pool), gate, resistance, links, moment_origin)
    }

    fn step_linear_resistance_inner<P: KernelRunner>(
        &mut self,
        pool: Option<&P>,
        gate: &CancelGate,
        mut resistance: impl FnMut([usize; 3]) -> LinearResistance3,
        links: &[BoundaryLink3],
        origin: [f64; 3],
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        endpoint_checkpoint(gate)?;
        self.validate_wall_measurement(links, origin)?;
        if Face3::ALL
            .iter()
            .any(|&face| self.boundaries.face(face).is_open())
        {
            return Err(BoundaryStepError3::EndpointForce {
                reason: "implicit linear resistance requires closed or periodic boundaries",
            });
        }
        let mut prepared = vec![PreparedResistance3::default(); self.f[0].len() * TILE_CELLS];
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if ((z * self.ny + y) * self.nx + x).is_multiple_of(TILE_CELLS) {
                        endpoint_checkpoint(gate)?;
                    }
                    if self.is_solid(x, y, z) {
                        continue;
                    }
                    let cell = [x, y, z];
                    let law = resistance(cell);
                    law.admit()?;
                    let (tile, lane) = self.addr(x, y, z);
                    let populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                    let (rho, momentum) =
                        raw_moments(&populations).ok_or(BoundaryStepError3::Unphysical { cell })?;
                    let previous = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    let velocity =
                        core::array::from_fn(|a| (momentum[a] + 0.5 * previous[a]) / rho);
                    if !velocity
                        .iter()
                        .chain(&law.force(velocity))
                        .all(|v| v.is_finite())
                    {
                        return Err(BoundaryStepError3::Unphysical { cell });
                    }
                    prepared[tile * TILE_CELLS + lane] = PreparedResistance3 { law, velocity };
                }
            }
        }
        {
            let kernel = CollidePass {
                f: &self.f,
                solid: &self.solid,
                post: chunk_groups(&mut self.post),
                model: self.collision_model,
                models: None,
                force: self.force,
                force_field: self.force_field.as_deref(),
                resistance: Some(&prepared),
            };
            let result = if let Some(pool) = pool {
                run_pass(pool, gate, &kernel)?
            } else {
                let mut result = FirstCollision(None);
                for group in 0..kernel.post.len() {
                    endpoint_checkpoint(gate)?;
                    result = kernel.run_group(group);
                    if result.0.is_some() {
                        break;
                    }
                }
                result
            };
            if let FirstCollision(Some((index, source))) = result {
                return Err(BoundaryStepError3::Collision {
                    cell: self.coords(index / TILE_CELLS, index % TILE_CELLS),
                    source,
                });
            }
        }
        let walls = self.wall_momentum(links, origin, Some(gate))?;
        let mut next = self.next.take().unwrap_or_else(|| self.f.clone());
        let mut forces = vec![[0.0; 3]; prepared.len()];
        let outcome = (|| {
            let field = ResistanceStream3 {
                prepared: &prepared,
                origin,
                groups: forces
                    .chunks_mut(POOLED_GROUP_TILES * TILE_CELLS)
                    .map(|force| {
                        Mutex::new(ResistanceGroup3 {
                            force,
                            exchange: ResistanceMomentumExchange3::default(),
                        })
                    })
                    .collect(),
            };
            {
                let kernel = StreamPass {
                    view: StreamView {
                        post: &self.post,
                        pull: &self.pull,
                        link_masks: &self.link_masks,
                        open_link_masks: &self.open_link_masks,
                        stationary_link_masks: &self.stationary_link_masks,
                        boundaries: self.boundaries,
                        dims: [self.nx, self.ny, self.nz],
                    },
                    solid: &self.solid,
                    endpoint_field: None,
                    resistance: Some(&field),
                    next: chunk_groups(&mut next),
                };
                let result = if let Some(pool) = pool {
                    run_pass(pool, gate, &kernel)?
                } else {
                    let mut result = FirstUnphysical(None);
                    for group in 0..kernel.next.len() {
                        endpoint_checkpoint(gate)?;
                        result = kernel.run_group(group);
                        if result.0.is_some() {
                            break;
                        }
                    }
                    result
                };
                if let FirstUnphysical(Some(index)) = result {
                    return Err(BoundaryStepError3::Unphysical {
                        cell: self.coords(index / TILE_CELLS, index % TILE_CELLS),
                    });
                }
            }
            let mut exchange = ResistanceMomentumExchange3::default();
            for group in &field.groups {
                exchange.add(
                    group
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner)
                        .exchange,
                );
            }
            if !exchange.is_finite() {
                return Err(BoundaryStepError3::EndpointForce {
                    reason: "linear-resistance aggregate impulse or work is not representable",
                });
            }
            exchange.walls = walls;
            Ok(exchange)
        })();
        let exchange = match outcome {
            Ok(exchange) => exchange,
            Err(error) => {
                self.next = Some(next);
                return Err(error);
            }
        };
        self.publish_next_with_forcing(
            next,
            gate,
            true,
            Some(EndpointForcing3::Field {
                uniform: [0.0; 3],
                local: forces,
            }),
        )?;
        self.topology_locked = true;
        Ok(exchange)
    }
}

fn raw_moments(populations: &[f64; Q3]) -> Option<(f64, [f64; 3])> {
    let mut rho = 0.0;
    let mut momentum = [0.0; 3];
    for (value, (x, y, z)) in populations.iter().zip(E3) {
        if !value.is_finite() {
            return None;
        }
        rho += value;
        momentum[0] += f64::from(x) * value;
        momentum[1] += f64::from(y) * value;
        momentum[2] += f64::from(z) * value;
    }
    (rho.is_finite() && rho > 0.0 && momentum.iter().all(|v| v.is_finite()))
        .then_some((rho, momentum))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
