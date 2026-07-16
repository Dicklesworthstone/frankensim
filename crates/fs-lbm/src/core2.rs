//! The general D2Q9 core (bead tfz.19): cell flags, VECTOR gravity
//! (tilt schedules rotate it), per-cell relaxation time (the
//! non-Newtonian hook), Guo forcing, pull streaming with halfway
//! bounce-back at walls — the substrate the thermal, rheology,
//! refinement, and free-surface extensions all share. Deterministic:
//! fixed row-major cell order, no RNG.

use crate::{CS2, E, OPP, Q, W};

const MAX_REGULARIZED_BOUNDARY_SPEED_SQ: f64 = 0.03;

/// Cell classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cell {
    /// Bulk fluid.
    Fluid,
    /// Solid wall (halfway bounce-back).
    Wall,
    /// Free-surface interface cell (carries partial mass).
    Interface,
    /// Gas cell (no populations).
    Gas,
}

/// The general D2Q9 lattice.
#[derive(Debug, Clone)]
pub struct Grid {
    /// Cells in x.
    pub nx: usize,
    /// Cells in y.
    pub ny: usize,
    /// Cell flags.
    pub flags: Vec<Cell>,
    /// Populations.
    pub f: Vec<[f64; Q]>,
    /// Per-cell relaxation time.
    pub tau: Vec<f64>,
    /// Gravity vector (lattice units).
    pub g: [f64; 2],
    /// Per-cell external force (Boussinesq buoyancy etc.), added to
    /// ρ·g.
    pub fext: Vec<[f64; 2]>,
    /// Periodic in x?
    pub periodic_x: bool,
    /// Periodic in y?
    pub periodic_y: bool,
}

/// Macroscopic moments of one cell.
#[derive(Debug, Clone, Copy)]
pub struct Moments {
    /// Density.
    pub rho: f64,
    /// Velocity (force-corrected).
    pub u: [f64; 2],
}

/// Momentum delivered to a selected set of stationary halfway-bounce-back
/// wall cells during one D2Q9 stream.
///
/// `wall_impulse` uses lattice momentum units. With the lattice time step
/// equal to one it is also the raw lattice force; callers remain responsible
/// for physical-unit and drag/lift normalization.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct MomentumExchange2 {
    /// Net `(x, y)` impulse delivered by the fluid to selected walls.
    pub wall_impulse: [f64; 2],
    /// Number of selected fluid-wall links included in the sum.
    pub measured_links: usize,
}

/// Momentum, torque, and work delivered to selected moving D2Q9 walls.
///
/// `wall_impulse` uses the boundary-relative momentum-exchange expression
/// `(c_out - u_wall) f_out - (c_in - u_wall) f_in`. The separate population
/// and wall-velocity mass terms retain the discrete balance
///
/// `wall_impulse + fluid_population_impulse = wall_velocity_mass_impulse`
///
/// up to floating-point roundoff. All values are raw lattice quantities for
/// one stream step; no physical-unit or reference-area normalization is
/// applied. The relative-velocity force convention follows Wen et al.,
/// *Galilean Invariant Fluid-Solid Interfacial Dynamics in Lattice Boltzmann
/// Simulations* (2014; <https://arxiv.org/abs/1303.0625>).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct MovingWallMomentumExchange2 {
    /// Boundary-relative hydrodynamic impulse delivered to selected walls.
    pub wall_impulse: [f64; 2],
    /// Change in resolved fluid-population momentum across selected links.
    pub fluid_population_impulse: [f64; 2],
    /// Net incoming-minus-outgoing population mass across selected links.
    pub fluid_mass_change: f64,
    /// Sum of `u_wall * (f_in - f_out)` across selected links.
    pub wall_velocity_mass_impulse: [f64; 2],
    /// Scalar 2-D wall angular impulse about the requested origin.
    pub wall_angular_impulse: f64,
    /// Scalar angular impulse of the resolved fluid-population change.
    pub fluid_population_angular_impulse: f64,
    /// Scalar angular impulse of `wall_velocity_mass_impulse`.
    pub wall_velocity_mass_angular_impulse: f64,
    /// Work delivered to selected walls, `sum(wall_impulse_link dot u_wall)`.
    pub wall_work: f64,
    /// Number of selected fluid-wall links included in the receipt.
    pub measured_links: usize,
}

/// Low-Mach regularized velocity inlet at `x = 0` and isothermal
/// pressure/density outlet at `x = nx - 1`.
///
/// The inlet copies density and non-equilibrium stress from `x = 1` while
/// prescribing velocity. The outlet copies velocity and stress from
/// `x = nx - 2` while prescribing density. This boundary pair is intended for
/// x-directed crossflow fixtures with periodic y closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VelocityPressureX2 {
    inlet_velocity: [f64; 2],
    outlet_density: f64,
}

impl VelocityPressureX2 {
    /// Construct a checked low-Mach boundary pair.
    #[must_use]
    pub fn new(inlet_velocity: [f64; 2], outlet_density: f64) -> Self {
        assert!(
            inlet_velocity.iter().all(|component| component.is_finite()),
            "D2Q9 inlet velocity must be finite"
        );
        let speed_sq = inlet_velocity
            .iter()
            .map(|component| component * component)
            .sum::<f64>();
        assert!(
            speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
            "D2Q9 inlet velocity exceeds the low-Mach boundary envelope"
        );
        assert!(
            outlet_density.is_finite() && outlet_density > 0.0,
            "D2Q9 outlet density must be positive and finite"
        );
        Self {
            inlet_velocity,
            outlet_density,
        }
    }

    /// Prescribed inlet velocity in lattice units.
    #[must_use]
    pub const fn inlet_velocity(self) -> [f64; 2] {
        self.inlet_velocity
    }

    /// Prescribed outlet density in lattice units.
    #[must_use]
    pub const fn outlet_density(self) -> f64 {
        self.outlet_density
    }
}

impl Grid {
    /// A grid of fluid at rest (unit density), uniform `tau`.
    #[must_use]
    pub fn uniform(nx: usize, ny: usize, tau: f64) -> Grid {
        assert!(nx > 0 && ny > 0, "grid dimensions must be positive");
        assert!(
            tau.is_finite() && tau > 0.5,
            "relaxation time tau must be finite and greater than 0.5"
        );
        let f0 = crate::equilibrium(1.0, 0.0, 0.0);
        Grid {
            nx,
            ny,
            flags: vec![Cell::Fluid; nx * ny],
            f: vec![f0; nx * ny],
            tau: vec![tau; nx * ny],
            g: [0.0, 0.0],
            fext: vec![[0.0; 2]; nx * ny],
            periodic_x: true,
            periodic_y: true,
        }
    }

    /// Row-major index.
    #[must_use]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Moments of cell `i` (Guo half-force correction).
    #[must_use]
    pub fn moments(&self, i: usize) -> Moments {
        let f = &self.f[i];
        let rho: f64 = f.iter().sum();
        assert!(
            rho.is_finite() && rho > 0.0,
            "moments require positive finite density"
        );
        let mut m = [0.0f64; 2];
        for (q, fi) in f.iter().enumerate() {
            m[0] += f64::from(E[q].0) * fi;
            m[1] += f64::from(E[q].1) * fi;
        }
        let fx = self.g[0].mul_add(rho, self.fext[i][0]);
        let fy = self.g[1].mul_add(rho, self.fext[i][1]);
        Moments {
            rho,
            u: [(m[0] + 0.5 * fx) / rho, (m[1] + 0.5 * fy) / rho],
        }
    }

    /// Total mass over non-gas cells.
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .zip(&self.flags)
            .filter(|&(_, &fl)| fl != Cell::Gas && fl != Cell::Wall)
            .map(|(c, _)| c.iter().sum::<f64>())
            .sum()
    }

    /// Collide (per-cell tau, vector Guo forcing) into `post`.
    pub fn collide_into(&self, post: &mut Vec<[f64; Q]>) {
        post.clear();
        post.resize(self.nx * self.ny, [0.0; Q]);
        for (i, out) in post.iter_mut().enumerate().take(self.nx * self.ny) {
            if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                *out = self.f[i];
                continue;
            }
            let mm = self.moments(i);
            let (rho, ux, uy) = (mm.rho, mm.u[0], mm.u[1]);
            let feq = crate::equilibrium(rho, ux, uy);
            let tau = self.tau[i];
            assert!(
                tau.is_finite() && tau > 0.5,
                "cell relaxation time tau must be finite and greater than 0.5"
            );
            let coef = 1.0 - 0.5 / tau;
            let (gx, gy) = (
                self.g[0].mul_add(rho, self.fext[i][0]),
                self.g[1].mul_add(rho, self.fext[i][1]),
            );
            for q in 0..Q {
                let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
                let eu = ex * ux + ey * uy;
                // Guo forcing, vector form.
                let fx = (ex - ux) / CS2 + eu * ex / (CS2 * CS2);
                let fy = (ey - uy) / CS2 + eu * ey / (CS2 * CS2);
                let force = coef * W[q] * (fx * gx + fy * gy);
                out[q] = self.f[i][q] + (feq[q] - self.f[i][q]) / tau + force;
            }
        }
    }

    /// Source cell for pull-streaming direction `q` into (x, y);
    /// `None` when the pull crosses a non-periodic boundary (treated
    /// as wall bounce-back).
    #[must_use]
    pub fn source(&self, x: usize, y: usize, q: usize) -> Option<usize> {
        let (ex, ey) = E[q];
        let sx = match ex {
            1 => {
                if x == 0 {
                    if self.periodic_x {
                        self.nx - 1
                    } else {
                        return None;
                    }
                } else {
                    x - 1
                }
            }
            -1 => {
                if x + 1 == self.nx {
                    if self.periodic_x {
                        0
                    } else {
                        return None;
                    }
                } else {
                    x + 1
                }
            }
            _ => x,
        };
        let sy = match ey {
            1 => {
                if y == 0 {
                    if self.periodic_y {
                        self.ny - 1
                    } else {
                        return None;
                    }
                } else {
                    y - 1
                }
            }
            -1 => {
                if y + 1 == self.ny {
                    if self.periodic_y {
                        0
                    } else {
                        return None;
                    }
                } else {
                    y + 1
                }
            }
            _ => y,
        };
        Some(self.idx(sx, sy))
    }

    fn validate_stream_input(&self, post: &[[f64; Q]]) {
        assert!(
            post.len() >= self.nx * self.ny,
            "post-collision populations must cover every grid cell"
        );
    }

    fn validate_measured_walls(&self, measured_walls: &[bool]) {
        assert_eq!(
            measured_walls.len(),
            self.nx * self.ny,
            "measured-wall mask length must match the grid"
        );
        for (index, (&measured, &flag)) in measured_walls.iter().zip(&self.flags).enumerate() {
            assert!(
                !measured || flag == Cell::Wall,
                "measured-wall mask selects non-wall cell {index}"
            );
        }
    }

    fn validate_moving_wall_fields(&self, wall_velocities: &[[f64; 2]], moment_origin: [f64; 2]) {
        assert_eq!(
            wall_velocities.len(),
            self.nx * self.ny,
            "moving-wall velocity field length must match the grid"
        );
        assert!(
            moment_origin.into_iter().all(f64::is_finite),
            "moving-wall moment origin must be finite"
        );
        for (index, (&velocity, &flag)) in wall_velocities.iter().zip(&self.flags).enumerate() {
            assert!(
                velocity.into_iter().all(f64::is_finite),
                "moving-wall velocity at cell {index} must be finite"
            );
            let speed_sq = velocity[0].mul_add(velocity[0], velocity[1] * velocity[1]);
            assert!(
                speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
                "moving-wall velocity at cell {index} exceeds the low-Mach admission envelope"
            );
            assert!(
                flag == Cell::Wall || velocity == [0.0; 2],
                "moving-wall velocity field assigns motion to non-wall cell {index}"
            );
        }
    }

    fn validate_moving_wall_post(&self, post: &[[f64; Q]], wall_velocities: &[[f64; 2]]) {
        self.validate_stream_input(post);
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                let mut needs_density = false;
                for q in 0..Q {
                    let Some(source) = self.source(x, y, q) else {
                        continue;
                    };
                    if self.flags[source] != Cell::Wall {
                        continue;
                    }
                    assert!(
                        post[i][OPP[q]].is_finite(),
                        "moving-wall outgoing population must be finite"
                    );
                    needs_density |= wall_velocities[source] != [0.0; 2];
                }
                if needs_density {
                    let rho_post = post[i].iter().sum::<f64>();
                    assert!(
                        rho_post.is_finite() && rho_post > 0.0,
                        "moving-wall bounce-back requires positive finite post-collision density"
                    );
                }
            }
        }
    }

    fn validate_velocity_pressure_x(&self) {
        assert!(
            self.nx >= 3,
            "D2Q9 x-open flow requires at least three columns"
        );
        assert!(
            !self.periodic_x && self.periodic_y,
            "D2Q9 x-open flow requires non-periodic x and periodic y"
        );
        assert!(
            self.g.into_iter().all(|component| component == 0.0)
                && self
                    .fext
                    .iter()
                    .flatten()
                    .all(|component| *component == 0.0),
            "D2Q9 regularized open boundaries currently require zero body force"
        );
        for y in 0..self.ny {
            for x in [0, 1, self.nx - 2, self.nx - 1] {
                assert_eq!(
                    self.flags[self.idx(x, y)],
                    Cell::Fluid,
                    "D2Q9 open faces and first-interior columns must be fluid"
                );
            }
        }
    }

    fn regularized_source(&self, index: usize) -> (f64, [f64; 2], [[f64; 2]; 2]) {
        let moments = self.moments(index);
        let equilibrium = crate::equilibrium(moments.rho, moments.u[0], moments.u[1]);
        let mut stress = [[0.0; 2]; 2];
        for q in 0..Q {
            let nonequilibrium = self.f[index][q] - equilibrium[q];
            let e = [f64::from(E[q].0), f64::from(E[q].1)];
            for row in 0..2 {
                for column in 0..2 {
                    stress[row][column] += e[row] * e[column] * nonequilibrium;
                }
            }
        }
        (moments.rho, moments.u, stress)
    }

    fn apply_velocity_pressure_x(&mut self, boundary: VelocityPressureX2) {
        for y in 0..self.ny {
            let inlet = self.idx(0, y);
            let inlet_source = self.idx(1, y);
            let outlet_source = self.idx(self.nx - 2, y);
            let outlet = self.idx(self.nx - 1, y);
            let (inlet_density, _, inlet_stress) = self.regularized_source(inlet_source);
            let (_, outlet_velocity, outlet_stress) = self.regularized_source(outlet_source);
            self.f[inlet] =
                regularized_populations(inlet_density, boundary.inlet_velocity, inlet_stress);
            self.f[outlet] =
                regularized_populations(boundary.outlet_density, outlet_velocity, outlet_stress);
        }
    }

    fn stream_from_inner(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: Option<&[bool]>,
    ) -> MomentumExchange2 {
        let mut receipt = MomentumExchange2::default();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                for q in 0..Q {
                    let pulled = match self.source(x, y, q) {
                        Some(s) if matches!(self.flags[s], Cell::Wall | Cell::Gas) => {
                            let reflected = post[i][OPP[q]];
                            if self.flags[s] == Cell::Wall
                                && measured_walls.is_some_and(|mask| mask[s])
                            {
                                // Pull direction q points from the wall back into
                                // the fluid. The fluid momentum change is
                                // +2 f_post c_q, so the opposite impulse delivered
                                // to the stationary wall is -2 f_post c_q.
                                receipt.wall_impulse[0] -= 2.0 * reflected * f64::from(E[q].0);
                                receipt.wall_impulse[1] -= 2.0 * reflected * f64::from(E[q].1);
                                receipt.measured_links += 1;
                            }
                            reflected
                        }
                        Some(s) => post[s][q],
                        None => post[i][OPP[q]],
                    };
                    self.f[i][q] = pulled;
                }
            }
        }
        receipt
    }

    fn stream_from_moving_wall_inner(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
        wall_velocities: &[[f64; 2]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        let mut receipt = MovingWallMomentumExchange2::default();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                let mut rho_post = None;
                for q in 0..Q {
                    let pulled = match self.source(x, y, q) {
                        Some(source) if self.flags[source] == Cell::Wall => {
                            let outgoing = post[i][OPP[q]];
                            let wall_velocity = wall_velocities[source];
                            let incoming = if wall_velocity == [0.0; 2] {
                                outgoing
                            } else {
                                let rho =
                                    *rho_post.get_or_insert_with(|| post[i].iter().sum::<f64>());
                                let e = [f64::from(E[q].0), f64::from(E[q].1)];
                                let e_dot_wall =
                                    e[0].mul_add(wall_velocity[0], e[1] * wall_velocity[1]);
                                outgoing + 2.0 * W[q] * rho * e_dot_wall / CS2
                            };

                            if measured_walls[source] {
                                let e = [f64::from(E[q].0), f64::from(E[q].1)];
                                let wall_impulse = if wall_velocity == [0.0; 2] {
                                    [-2.0 * outgoing * e[0], -2.0 * outgoing * e[1]]
                                } else {
                                    [
                                        (-e[0] - wall_velocity[0]) * outgoing
                                            - (e[0] - wall_velocity[0]) * incoming,
                                        (-e[1] - wall_velocity[1]) * outgoing
                                            - (e[1] - wall_velocity[1]) * incoming,
                                    ]
                                };
                                let fluid_population_impulse =
                                    [e[0] * (incoming + outgoing), e[1] * (incoming + outgoing)];
                                let fluid_mass_change = incoming - outgoing;
                                let wall_velocity_mass_impulse = [
                                    wall_velocity[0] * fluid_mass_change,
                                    wall_velocity[1] * fluid_mass_change,
                                ];
                                let link_offset = [
                                    (x as f64 - 0.5 * e[0]) - moment_origin[0],
                                    (y as f64 - 0.5 * e[1]) - moment_origin[1],
                                ];

                                receipt.wall_impulse[0] += wall_impulse[0];
                                receipt.wall_impulse[1] += wall_impulse[1];
                                receipt.fluid_population_impulse[0] += fluid_population_impulse[0];
                                receipt.fluid_population_impulse[1] += fluid_population_impulse[1];
                                receipt.fluid_mass_change += fluid_mass_change;
                                receipt.wall_velocity_mass_impulse[0] +=
                                    wall_velocity_mass_impulse[0];
                                receipt.wall_velocity_mass_impulse[1] +=
                                    wall_velocity_mass_impulse[1];
                                receipt.wall_angular_impulse += link_offset[0] * wall_impulse[1]
                                    - link_offset[1] * wall_impulse[0];
                                receipt.fluid_population_angular_impulse += link_offset[0]
                                    * fluid_population_impulse[1]
                                    - link_offset[1] * fluid_population_impulse[0];
                                receipt.wall_velocity_mass_angular_impulse += link_offset[0]
                                    * wall_velocity_mass_impulse[1]
                                    - link_offset[1] * wall_velocity_mass_impulse[0];
                                receipt.wall_work += wall_impulse[0]
                                    .mul_add(wall_velocity[0], wall_impulse[1] * wall_velocity[1]);
                                receipt.measured_links += 1;
                            }
                            incoming
                        }
                        Some(source) if self.flags[source] == Cell::Gas => post[i][OPP[q]],
                        Some(source) => post[source][q],
                        None => post[i][OPP[q]],
                    };
                    self.f[i][q] = pulled;
                }
            }
        }
        receipt
    }

    /// Stream `post` into `self.f` (fluid pull; wall and out-of-domain
    /// pulls bounce back).
    pub fn stream_from(&mut self, post: &[[f64; Q]]) {
        self.validate_stream_input(post);
        let _ = self.stream_from_inner(post, None);
    }

    /// Stream while measuring momentum delivered to selected wall cells.
    ///
    /// Only bounce-back links whose source cell is both [`Cell::Wall`] and
    /// `true` in `measured_walls` contribute. Gas-boundary and non-periodic
    /// exterior bounces are deliberately excluded.
    pub fn stream_from_with_wall_momentum(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
    ) -> MomentumExchange2 {
        self.validate_stream_input(post);
        self.validate_measured_walls(measured_walls);
        self.stream_from_inner(post, Some(measured_walls))
    }

    /// Stream with per-wall-cell velocities and a moving-wall exchange receipt.
    ///
    /// Moving halfway bounce-back is applied to every wall cell according to
    /// `wall_velocities`; `measured_walls` only selects which links contribute
    /// to the returned receipt. Torque uses each destination-local halfway-link
    /// midpoint about `moment_origin`. All request and post-collision inputs are
    /// validated before populations are mutated.
    pub fn stream_from_with_moving_wall_momentum(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
        wall_velocities: &[[f64; 2]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_moving_wall_fields(wall_velocities, moment_origin);
        self.validate_moving_wall_post(post, wall_velocities);
        self.stream_from_moving_wall_inner(post, measured_walls, wall_velocities, moment_origin)
    }

    /// One plain step (no free-surface bookkeeping).
    pub fn step(&mut self, scratch: &mut Vec<[f64; Q]>) {
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        self.stream_from(&post);
        *scratch = post;
    }

    /// One plain step plus a raw stationary-wall momentum-exchange receipt.
    ///
    /// The mask is validated before collision, so a malformed measurement
    /// request cannot partially advance the grid.
    pub fn step_with_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        measured_walls: &[bool],
    ) -> MomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_inner(&post, Some(measured_walls));
        *scratch = post;
        receipt
    }

    /// One collide-stream step with moving walls and a selected-wall receipt.
    ///
    /// The velocity field, measurement mask, moment origin, and post-collision
    /// state are all admitted before streaming can mutate `self.f`.
    pub fn step_with_moving_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        measured_walls: &[bool],
        wall_velocities: &[[f64; 2]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_moving_wall_fields(wall_velocities, moment_origin);
        self.collide_into(scratch);
        self.validate_moving_wall_post(scratch, wall_velocities);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_moving_wall_inner(
            &post,
            measured_walls,
            wall_velocities,
            moment_origin,
        );
        *scratch = post;
        receipt
    }

    fn step_velocity_pressure_x_inner(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        boundary: VelocityPressureX2,
        measured_walls: Option<&[bool]>,
    ) -> MomentumExchange2 {
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_inner(&post, measured_walls);
        self.apply_velocity_pressure_x(boundary);
        *scratch = post;
        receipt
    }

    /// One collide-stream step followed by regularized x-face reconstruction.
    pub fn step_velocity_pressure_x(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        boundary: VelocityPressureX2,
    ) {
        self.validate_velocity_pressure_x();
        let _ = self.step_velocity_pressure_x_inner(scratch, boundary, None);
    }

    /// Regularized x-flow step plus a selected-wall momentum receipt.
    ///
    /// Grid topology, periodicity, forcing, and the wall mask are validated
    /// before collision, so a refused request cannot partially advance state.
    pub fn step_velocity_pressure_x_with_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        boundary: VelocityPressureX2,
        measured_walls: &[bool],
    ) -> MomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_velocity_pressure_x();
        self.step_velocity_pressure_x_inner(scratch, boundary, Some(measured_walls))
    }
}

fn regularized_populations(rho: f64, velocity: [f64; 2], stress: [[f64; 2]; 2]) -> [f64; Q] {
    let mut populations = crate::equilibrium(rho, velocity[0], velocity[1]);
    let coefficient = 1.0 / (2.0 * CS2 * CS2);
    for q in 0..Q {
        let e = [f64::from(E[q].0), f64::from(E[q].1)];
        let mut contraction = 0.0;
        for row in 0..2 {
            for column in 0..2 {
                let isotropic = if row == column { CS2 } else { 0.0 };
                contraction += (e[row] * e[column] - isotropic) * stress[row][column];
            }
        }
        populations[q] += coefficient * W[q] * contraction;
    }
    populations
}

/// Strain-rate magnitude (sqrt(2 S:S)) of one cell from its
/// non-equilibrium populations — the LOCAL quantity non-Newtonian
/// relaxation adapts to. `feq` must match the cell's moments.
#[must_use]
pub fn shear_rate(f: &[f64; Q], feq: &[f64; Q], rho: f64, tau: f64) -> f64 {
    // S_ab = −(3 / (2 ρ τ)) Σ_q e_qa e_qb (f_q − feq_q)   (c_s² = 1/3)
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    let mut syy = 0.0f64;
    for q in 0..Q {
        let neq = f[q] - feq[q];
        let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
        sxx += ex * ex * neq;
        sxy += ex * ey * neq;
        syy += ey * ey * neq;
    }
    let c = -3.0 / (2.0 * rho * tau);
    let (sxx, sxy, syy) = (c * sxx, c * sxy, c * syy);
    let ss = 2.0f64.mul_add(sxy * sxy, sxx.mul_add(sxx, syy * syy));
    (2.0 * ss).sqrt()
}
