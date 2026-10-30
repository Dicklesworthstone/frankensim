//! One simultaneous fluid tick for several additional solid inertias.
//!
//! At a shared cell, `B_a = epsilon_a (tau-1/2)/(1-sum(epsilon)+tau-1/2)`.
//! The population correction is the sum of `B_a (Omega_s,a-Omega_f)`.
//! Each rigid body therefore has its own six-velocity midpoint block; fluid
//! populations and all bodies are published together, never streamed once per
//! body. No dense cross-body inertia is invented for this collision law.

use super::*;

/// One body's response in a simultaneous assembly. Fixed-wall exchange is
/// owned once by [`PartialSaturationBodiesExchange3::walls`].
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PartialSaturationBodyResponse3 {
    /// This body's population correction, impulse, torque and work.
    pub hydrodynamic: PartialSaturationExchange3,
    /// Translational midpoint velocity used in the collision.
    pub velocity_midpoint: [f64; 3],
    /// Angular midpoint velocity used in the collision.
    pub angular_velocity_midpoint: [f64; 3],
    /// Actual accepted body kinetic-energy increment.
    pub kinetic_energy_change: f64,
    /// Work of the external center impulse and couple.
    pub external_work: f64,
    /// Linear and center-angular momentum balance defects.
    pub momentum_residual: [f64; 6],
}

impl From<PartialSaturationBodyExchange3> for PartialSaturationBodyResponse3 {
    fn from(value: PartialSaturationBodyExchange3) -> Self {
        Self {
            hydrodynamic: value.hydrodynamic,
            velocity_midpoint: value.velocity_midpoint,
            angular_velocity_midpoint: value.angular_velocity_midpoint,
            kinetic_energy_change: value.kinetic_energy_change,
            external_work: value.external_work,
            momentum_residual: value.momentum_residual,
        }
    }
}

/// Simultaneous body responses and actual composite-population roundoff.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct PartialSaturationBodiesExchange3<const N: usize> {
    /// Stable caller order; each population correction has exactly one owner.
    pub bodies: [PartialSaturationBodyResponse3; N],
    /// All stationary voxel/exterior wall exchange, counted only once.
    pub walls: WallMomentumExchange3,
    /// Imposed fluid body-force impulse, torque and local midpoint work once
    /// over the whole fluid, including retained interior populations.
    pub forcing: PartialSaturationForcing3,
    /// Unresolved normal sphere-gap resistance, separate from the populations.
    /// Zero for the original uncorrected entry points.
    pub lubrication: SphereLubricationExchange3<N>,
    /// Actual rounded composite correction minus summed per-body corrections.
    pub population_mass_residual: f64,
    /// Composite versus summed corrections, for world linear momentum.
    pub population_impulse_residual: [f64; 3],
    /// Same measurement for angular momentum about the requested origin.
    pub population_angular_impulse_residual: [f64; 3],
}

#[derive(Clone, Copy, Default)]
pub(super) struct Response {
    weight: f64,
    first: [f64; 3],
    rotary: [[f64; 3]; 3],
    increment: [f64; 6],
}

impl Response {
    fn add(&mut self, weight: f64, r: [f64; 3], slip: [f64; 3]) {
        let torque = cross(r, slip);
        self.weight += weight;
        for i in 0..3 {
            self.first[i] = weight.mul_add(r[i], self.first[i]);
            self.increment[i] = weight.mul_add(slip[i], self.increment[i]);
            self.increment[i + 3] = weight.mul_add(torque[i], self.increment[i + 3]);
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            self.rotary[i][i] = weight.mul_add(r[j].mul_add(r[j], r[k] * r[k]), self.rotary[i][i]);
            for j in 0..i {
                self.rotary[i][j] = (-weight * r[i]).mul_add(r[j], self.rotary[i][j]);
            }
        }
    }

    pub(super) fn factor(self, body: PartialSaturationBody3) -> Result<fs_la::factor::Cholesky> {
        let mut matrix = [0.0; 36];
        for i in 0..3 {
            matrix[i * 6 + i] = body.mass + 0.5 * self.weight;
            for j in 0..=i {
                matrix[(i + 3) * 6 + j + 3] = body.inertia_world[i][j] + 0.5 * self.rotary[i][j];
            }
        }
        matrix[3 * 6 + 1] = -0.5 * self.first[2];
        matrix[3 * 6 + 2] = 0.5 * self.first[1];
        matrix[4 * 6] = 0.5 * self.first[2];
        matrix[4 * 6 + 2] = -0.5 * self.first[0];
        matrix[5 * 6] = -0.5 * self.first[1];
        matrix[5 * 6 + 1] = 0.5 * self.first[0];
        for i in 0..6 {
            for j in 0..i {
                matrix[j * 6 + i] = matrix[i * 6 + j];
            }
        }
        if !matrix.into_iter().chain(self.increment).all(f64::is_finite) {
            return Err(invalid("assembly spatial response overflowed"));
        }
        fs_la::factor::cholesky(&matrix, 6)
            .map_err(|_| invalid("assembly coupled inertia is not positive definite"))
    }

    pub(super) fn increment(self) -> [f64; 6] {
        self.increment
    }

    pub(super) fn solve(mut self, body: PartialSaturationBody3) -> Result<PartialSaturationBody3> {
        self.factor(body)?.solve(&mut self.increment);
        Self::apply(body, self.increment)
    }

    pub(super) fn apply(
        body: PartialSaturationBody3,
        increment: [f64; 6],
    ) -> Result<PartialSaturationBody3> {
        let proposed = PartialSaturationBody3 {
            velocity: core::array::from_fn(|i| body.velocity[i] + increment[i]),
            angular_velocity: core::array::from_fn(|i| body.angular_velocity[i] + increment[i + 3]),
            ..body
        };
        if !low_mach(proposed.velocity)
            || !proposed.angular_velocity.into_iter().all(f64::is_finite)
        {
            return Err(invalid(
                "assembly body response exceeds finite low-Mach envelope",
            ));
        }
        Ok(proposed)
    }
}

impl BoundaryGrid3 {
    /// Advance several rigid velocities and the fluid in one atomic tick.
    ///
    /// Fractions are row-major full-grid fields in stable body order. Their
    /// cellwise sum must lie in `[0,1]`; each field describes disjoint occupied
    /// subvolumes, a caller obligation for this geometry-free entry point.
    /// Different bodies may occupy distinct subvolumes of the same cell.
    /// Geometry and world inertia are frozen, and inertia is additional to
    /// retained interior fluid. Each body is solved implicitly at its midpoint;
    /// the fluid collision and streaming occur once. Shared-cell weights use
    /// total occupancy, with no clipping or fraction renormalization.
    ///
    /// The one-body case delegates to the existing bit-compatible owner.
    /// Temporary body work is fixed-size `O(N)`; fraction storage belongs to the
    /// caller and the native fluid `post`/`next` buffers are reused. Native
    /// Cholesky factors allocate two small buffers per body. Reductions follow
    /// body order: permutations agree to roundoff, not necessarily bit for bit.
    ///
    /// # Errors
    /// Inherits single-body BGK, held Guo forcing, finite/SPD, low-Mach, stationary
    /// wall/periodic admission. Also refuses an empty assembly, malformed
    /// fraction fields and aggregate occupancy above one. Every refusal leaves
    /// all bodies, fluid populations and topology admission unchanged.
    pub fn step_with_partial_saturation_bodies<const N: usize>(
        &mut self,
        fractions: [&[f64]; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<PartialSaturationBodiesExchange3<N>> {
        self.step_with_partial_saturation_bodies_checked(
            fractions,
            bodies,
            external,
            origin,
            gate,
            |_, _, _| Ok(()),
        )
        .map(|(receipt, ())| receipt)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::d3q19::boundary) fn step_with_partial_saturation_bodies_checked<
        const N: usize,
        T,
    >(
        &mut self,
        fractions: [&[f64]; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        gate: &CancelGate,
        finish: impl FnOnce(
            &Self,
            &[PartialSaturationBody3; N],
            &PartialSaturationBodiesExchange3<N>,
        ) -> Result<T>,
    ) -> Result<(PartialSaturationBodiesExchange3<N>, T)> {
        self.step_with_partial_saturation_resistance_checked(
            fractions, bodies, external, origin, gate, None, finish,
        )
        .map(|(report, value, _)| (report, value))
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::d3q19::boundary) fn step_with_partial_saturation_resistance_checked<
        const N: usize,
        T,
    >(
        &mut self,
        fractions: [&[f64]; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        gate: &CancelGate,
        lubrication: Option<(SphereLubrication3, [f64; N], &[LubricationPlane3])>,
        finish: impl FnOnce(
            &Self,
            &[PartialSaturationBody3; N],
            &PartialSaturationBodiesExchange3<N>,
        ) -> Result<T>,
    ) -> Result<(
        PartialSaturationBodiesExchange3<N>,
        T,
        Vec<PlaneLubricationExchange3>,
    )> {
        checkpoint(gate)?;
        if N == 0 {
            return Err(invalid("immersed assembly must contain at least one body"));
        }
        if let Some((model, radii, planes)) = lubrication {
            model.admit(radii)?;
            lubrication::planes::admit_planes(planes)?;
        }
        if N == 1 && lubrication.is_none_or(|(_, _, planes)| planes.is_empty()) {
            let mut proposed = *bodies;
            let mut single = proposed[0];
            let (_, (receipt, value)) = self.step_with_partial_saturation_body_checked(
                fractions[0],
                &mut single,
                external[0],
                origin,
                gate,
                |grid, after, exchange| {
                    proposed[0] = *after;
                    let receipt = PartialSaturationBodiesExchange3 {
                        bodies: [(*exchange).into(); N],
                        walls: exchange.walls,
                        forcing: exchange.forcing,
                        lubrication: SphereLubricationExchange3::default(),
                        population_mass_residual: 0.0,
                        population_impulse_residual: [0.0; 3],
                        population_angular_impulse_residual: [0.0; 3],
                    };
                    let value = finish(grid, &proposed, &receipt)?;
                    Ok((receipt, value))
                },
            )?;
            *bodies = proposed;
            return Ok((receipt, value, Vec::new()));
        }
        let CollisionModel3::Bgk { tau } = self.collision_model else {
            return Err(invalid("immersed assembly requires BGK"));
        };
        if !origin.into_iter().all(f64::is_finite)
            || Face3::ALL.into_iter().any(|face| {
                !matches!(self.boundaries.face(face), FaceBoundary3::Periodic)
                && !matches!(self.boundaries.face(face), FaceBoundary3::Wall { velocity } if velocity == [0.0; 3])
            }) {
            return Err(invalid("assembly requires finite origin and periodic/stationary faces"));
        }
        let mut response = [Response::default(); N];
        for a in 0..N {
            checkpoint(gate)?;
            let body = bodies[a];
            if fractions[a].len() != self.nx * self.ny * self.nz
                || !(body.mass.is_finite() && body.mass > 0.0)
                || !body
                    .center
                    .into_iter()
                    .chain(body.angular_velocity)
                    .chain(body.inertia_world.into_iter().flatten())
                    .chain(external[a])
                    .all(f64::is_finite)
                || !low_mach(body.velocity)
            {
                return Err(invalid(
                    "assembly body, fractions or impulse is inadmissible",
                ));
            }
            if (0..3).any(|i| (0..i).any(|j| body.inertia_world[i][j] != body.inertia_world[j][i]))
            {
                return Err(invalid("assembly world inertia must be symmetric"));
            }
            let inertia: [f64; 9] = core::array::from_fn(|i| body.inertia_world[i / 3][i % 3]);
            fs_la::factor::cholesky(&inertia, 3)
                .map_err(|_| invalid("assembly world inertia must be positive definite"))?;
            response[a].increment = external[a];
        }
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = self.coords(tile, lane);
                let index = (cell[2] * self.ny + cell[1]) * self.nx + cell[0];
                let total = admit_fractions(
                    &fractions,
                    index,
                    cell,
                    self.solid[tile] & (1u64 << lane) != 0,
                    gate,
                )?;
                if total == 0.0 {
                    continue;
                }
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                let (rho, u) = moments(self.f.each_ref().map(|f| f[tile].0[lane]), force, cell)?;
                for a in 0..N {
                    if a.is_multiple_of(TILE_CELLS) {
                        checkpoint(gate)?;
                    }
                    let epsilon = fractions[a][index];
                    if epsilon == 0.0 {
                        continue;
                    }
                    let r = core::array::from_fn(|i| cell[i] as f64 + 0.5 - bodies[a].center[i]);
                    let v = point_velocity(bodies[a].velocity, bodies[a].angular_velocity, r);
                    if !low_mach(v) {
                        return Err(PartialSaturationError3::Cell {
                            cell,
                            reason: "assembly initial point exceeds low-Mach envelope",
                        });
                    }
                    let weight = epsilon * (tau - 0.5) / (1.0 - total + tau - 0.5) * rho;
                    response[a].add(weight, r, core::array::from_fn(|i| u[i] - v[i]));
                }
            }
        }
        let mut proposed = *bodies;
        let mut report = PartialSaturationBodiesExchange3 {
            bodies: [PartialSaturationBodyResponse3::default(); N],
            walls: WallMomentumExchange3::default(),
            forcing: PartialSaturationForcing3::default(),
            lubrication: SphereLubricationExchange3::default(),
            population_mass_residual: 0.0,
            population_impulse_residual: [0.0; 3],
            population_angular_impulse_residual: [0.0; 3],
        };
        let plane_reports;
        if let Some((model, radii, planes)) = lubrication {
            (report.lubrication, plane_reports) = lubrication::solve(
                model,
                radii,
                planes,
                bodies,
                &response,
                &mut proposed,
                origin,
                gate,
            )?;
        } else {
            plane_reports = Vec::new();
            for a in 0..N {
                checkpoint(gate)?;
                proposed[a] = response[a].solve(bodies[a])?;
            }
        }
        for a in 0..N {
            report.bodies[a].velocity_midpoint =
                core::array::from_fn(|i| 0.5 * (bodies[a].velocity[i] + proposed[a].velocity[i]));
            report.bodies[a].angular_velocity_midpoint = core::array::from_fn(|i| {
                0.5 * (bodies[a].angular_velocity[i] + proposed[a].angular_velocity[i])
            });
        }
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let cell = self.coords(tile, lane);
                let index = (cell[2] * self.ny + cell[1]) * self.nx + cell[0];
                let f = self.f.each_ref().map(|f| f[tile].0[lane]);
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                let ordinary = collide_cell3(f, self.collision_model, force)
                    .map_err(|source| PartialSaturationError3::Collision { cell, source })?;
                let mut post = ordinary;
                let mut total = 0.0;
                for (a, field) in fractions.iter().enumerate() {
                    if a.is_multiple_of(TILE_CELLS) {
                        checkpoint(gate)?;
                    }
                    total += field[index];
                }
                if total > 0.0 {
                    let (rho, velocity) = moments(f, force, cell)?;
                    let physical = physical_populations(f, velocity, force);
                    let eq = equilibrium3(rho, velocity);
                    let mut summed = [0.0; Q3];
                    for a in 0..N {
                        if a.is_multiple_of(TILE_CELLS) {
                            checkpoint(gate)?;
                        }
                        let epsilon = fractions[a][index];
                        if epsilon == 0.0 {
                            continue;
                        }
                        let response = &mut report.bodies[a];
                        let r =
                            core::array::from_fn(|i| cell[i] as f64 + 0.5 - bodies[a].center[i]);
                        let v = point_velocity(
                            response.velocity_midpoint,
                            response.angular_velocity_midpoint,
                            r,
                        );
                        if !low_mach(v)
                            || !low_mach(point_velocity(
                                proposed[a].velocity,
                                proposed[a].angular_velocity,
                                r,
                            ))
                        {
                            return Err(PartialSaturationError3::Cell {
                                cell,
                                reason: "assembly solved point exceeds low-Mach envelope",
                            });
                        }
                        let solid_eq = equilibrium3(rho, v);
                        let b = epsilon * (tau - 0.5) / (1.0 - total + tau - 0.5);
                        let mut impulse = [0.0; 3];
                        for q in 0..Q3 {
                            let fluid_collision = (eq[q] - physical[q]) / tau;
                            let solid_collision =
                                physical[OPP3[q]] - physical[q] + solid_eq[q] - eq[OPP3[q]];
                            let correction = b * (solid_collision - fluid_collision);
                            summed[q] += correction;
                            response.hydrodynamic.fluid_population_mass_change += correction;
                            impulse[0] += f64::from(E3[q].0) * correction;
                            impulse[1] += f64::from(E3[q].1) * correction;
                            impulse[2] += f64::from(E3[q].2) * correction;
                        }
                        let arm = core::array::from_fn(|i| cell[i] as f64 + 0.5 - origin[i]);
                        let angular = cross(arm, impulse);
                        let spin = cross(r, impulse);
                        for i in 0..3 {
                            response.hydrodynamic.fluid_population_impulse[i] += impulse[i];
                            response.hydrodynamic.body_impulse[i] -= impulse[i];
                            response.hydrodynamic.fluid_population_angular_impulse[i] += angular[i];
                            response.hydrodynamic.body_angular_impulse[i] -= angular[i];
                            response.hydrodynamic.body_angular_impulse_about_center[i] -= spin[i];
                        }
                        response.hydrodynamic.body_work -= dot(impulse, v);
                        response.hydrodynamic.coupled_cells += 1;
                        response.hydrodynamic.coupling_weight_sum += b;
                    }
                    let mut residual = [0.0; 3];
                    for q in 0..Q3 {
                        post[q] = ordinary[q] + summed[q];
                        let difference = (post[q] - ordinary[q]) - summed[q];
                        report.population_mass_residual += difference;
                        residual[0] += f64::from(E3[q].0) * difference;
                        residual[1] += f64::from(E3[q].1) * difference;
                        residual[2] += f64::from(E3[q].2) * difference;
                    }
                    let arm = core::array::from_fn(|i| cell[i] as f64 + 0.5 - origin[i]);
                    let angular = cross(arm, residual);
                    for i in 0..3 {
                        report.population_impulse_residual[i] += residual[i];
                        report.population_angular_impulse_residual[i] += angular[i];
                    }
                }
                if !post.into_iter().all(f64::is_finite) {
                    return Err(PartialSaturationError3::Cell {
                        cell,
                        reason: "assembly blended collision overflowed",
                    });
                }
                report.forcing.add(f, post, force, cell, origin)?;
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
                let mask = self.link_masks[tile].0[lane];
                if mask != 0 {
                    let rho = post.iter().sum();
                    for q in 1..Q3 {
                        if mask & (1u32 << q) != 0 {
                            report.walls.add_link(
                                BoundaryLink3 {
                                    tile,
                                    lane,
                                    cell,
                                    direction: q,
                                },
                                post[OPP3[q]],
                                rho,
                                [0.0; 3],
                                origin,
                            );
                        }
                    }
                }
            }
        }
        report.walls = report
            .walls
            .checked()
            .map_err(|_| invalid("assembly wall exchange overflowed"))?;
        report.forcing = report.forcing.checked()?;
        for a in 0..N {
            checkpoint(gate)?;
            let response = &mut report.bodies[a];
            let dp = core::array::from_fn(|i| {
                bodies[a].mass * (proposed[a].velocity[i] - bodies[a].velocity[i])
            });
            let dw = core::array::from_fn(|i| {
                proposed[a].angular_velocity[i] - bodies[a].angular_velocity[i]
            });
            let dl = bodies[a].inertia_world.map(|row| dot(row, dw));
            response.kinetic_energy_change =
                dot(dp, response.velocity_midpoint) + dot(dl, response.angular_velocity_midpoint);
            response.external_work = dot(
                external[a][..3].try_into().unwrap(),
                response.velocity_midpoint,
            ) + dot(
                external[a][3..].try_into().unwrap(),
                response.angular_velocity_midpoint,
            );
            response.momentum_residual = core::array::from_fn(|i| {
                if i < 3 {
                    dp[i]
                        - external[a][i]
                        - response.hydrodynamic.body_impulse[i]
                        - report.lubrication.body_impulses[a][i]
                } else {
                    dl[i - 3]
                        - external[a][i]
                        - response.hydrodynamic.body_angular_impulse_about_center[i - 3]
                }
            });
            let hydro = response.hydrodynamic;
            if ![
                response.kinetic_energy_change,
                response.external_work,
                hydro.fluid_population_mass_change,
                hydro.body_work,
                hydro.coupling_weight_sum,
            ]
            .into_iter()
            .chain(response.momentum_residual)
            .chain(hydro.fluid_population_impulse)
            .chain(hydro.fluid_population_angular_impulse)
            .chain(hydro.body_angular_impulse_about_center)
            .all(f64::is_finite)
            {
                return Err(invalid("assembly body balance overflowed"));
            }
        }
        if ![report.population_mass_residual]
            .into_iter()
            .chain(report.population_impulse_residual)
            .chain(report.population_angular_impulse_residual)
            .all(f64::is_finite)
        {
            return Err(invalid("assembly population balance overflowed"));
        }
        let mut next = match self.next.take() {
            Some(next) => next,
            None => {
                let mut next = core::array::from_fn(|_| Vec::new());
                for field in &mut next {
                    field
                        .try_reserve_exact(self.f[0].len())
                        .map_err(|_| invalid("assembly streaming scratch allocation refused"))?;
                }
                for tile in 0..self.f[0].len() {
                    checkpoint(gate)?;
                    for (q, field) in next.iter_mut().enumerate() {
                        field.push(self.f[q][tile].clone());
                    }
                }
                next
            }
        };
        let view = StreamView {
            post: &self.post,
            pull: &self.pull,
            link_masks: &self.link_masks,
            open_link_masks: &self.open_link_masks,
            stationary_link_masks: &self.stationary_link_masks,
            boundaries: self.boundaries,
            dims: [self.nx, self.ny, self.nz],
        };
        for tile in 0..self.f[0].len() {
            if let Err(error) = checkpoint(gate) {
                self.next = Some(next);
                return Err(error);
            }
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) == 0 {
                    let post = view.pull_cell(tile, lane);
                    for (field, value) in next.iter_mut().zip(post) {
                        field[tile].0[lane] = value;
                    }
                }
            }
        }
        let value = match checkpoint(gate).and_then(|()| finish(self, &proposed, &report)) {
            Ok(value) => value,
            Err(error) => {
                self.next = Some(next);
                return Err(error);
            }
        };
        if let Err(error) = checkpoint(gate) {
            self.next = Some(next);
            return Err(error);
        }
        core::mem::swap(&mut self.f, &mut next);
        self.next = Some(next);
        self.topology_locked = true;
        *bodies = proposed;
        Ok((report, value, plane_reports))
    }
}

fn admit_fractions<const N: usize>(
    fields: &[&[f64]; N],
    index: usize,
    cell: [usize; 3],
    solid: bool,
    gate: &CancelGate,
) -> Result<f64> {
    let mut total = 0.0;
    for (a, field) in fields.iter().enumerate() {
        if a.is_multiple_of(TILE_CELLS) {
            checkpoint(gate)?;
        }
        let epsilon = field[index];
        if !(epsilon.is_finite() && (0.0..=1.0).contains(&epsilon)) || (solid && epsilon > 0.0) {
            return Err(PartialSaturationError3::Cell {
                cell,
                reason: "assembly fraction must be in [0,1] and positive only on fluid",
            });
        }
        total += epsilon;
    }
    if !(total.is_finite() && total <= 1.0) {
        return Err(PartialSaturationError3::Cell {
            cell,
            reason: "assembly solid fractions exceed cell volume",
        });
    }
    Ok(total)
}
