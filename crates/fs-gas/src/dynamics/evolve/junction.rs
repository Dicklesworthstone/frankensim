//! Isentropic two-sided face at an abrupt radius step of a quasi-one-dimensional duct.
//!
//! The quasi-one-dimensional equations assume a smoothly varying area. A radius step inside
//! one cell smears that jump over the cell and, in the finite-volume balance, turns it into
//! a numerical energy sink (about 0.95 Gamma^2 of the acoustic energy per traversal, Gamma
//! the step's reflection coefficient). A junction face keeps the step where it is: the
//! cell on each side sees its own side's area, and the face flux is the solution of a
//! stationary, adiabatic and isentropic area change between the two wave curves of the
//! neighbouring traces. In small signal this is exactly the plane-wave transfer matrix of the
//! step (equal pressure, equal volume flow); at finite amplitude it carries the Bernoulli
//! pressure difference of the source-owned isentrope. It has no friction or mixing loss: a
//! Borda-Carnot loss at an expansion is a different closure with its own printed source.
//!
//! Mass and total enthalpy are conserved across the step, so the mass and energy of the two
//! neighbouring cells telescope; the momentum difference is the axial force of the step
//! annulus on the gas, reported with the geometric area force.

use super::*;
use crate::dynamics::wall::shared::CaloricLaw;
use fs_math::det;

/// The radius steps of a stepped duct, in increasing face order, with their per-stage fluxes.
#[derive(Clone, Debug, Default)]
pub(super) struct StepSet {
    /// Interior face of each step. `Mesh1d::face_areas()[face]` is the area the cell on the
    /// left sees there.
    faces: Vec<usize>,
    /// The area the cell on the right sees.
    area_after: Vec<f64>,
    /// The flux the cell on the right receives, per unit of `area_after`.
    after: Vec<FaceFlux>,
}

impl StepSet {
    pub(super) fn new(mesh: &Mesh1d, steps: &[(usize, f64)]) -> Result<Self, Error> {
        let n = mesh.cell_count();
        let mut set = Self::default();
        let mut previous = 0;
        for &(face, area_after) in steps {
            let area_before = mesh.face_areas()[face.min(n)];
            if face <= previous || face >= n || !positive(area_after) || area_after == area_before {
                return Err(Error::Invalid(
                    "radius steps need increasing interior faces and a different positive area",
                ));
            }
            previous = face;
            set.faces.push(face);
            set.area_after.push(area_after);
            set.after.push(FaceFlux::default());
        }
        Ok(set)
    }
    pub(super) fn is_empty(&self) -> bool {
        self.faces.is_empty()
    }
    /// `(face, area seen by the cell on the right)` of every step.
    pub(super) fn steps(&self) -> Vec<(usize, f64)> {
        self.faces
            .iter()
            .copied()
            .zip(self.area_after.iter().copied())
            .collect()
    }
    fn index(&self, face: usize) -> Option<usize> {
        if self.faces.is_empty() {
            return None;
        }
        self.faces.binary_search(&face).ok()
    }
    /// The factor that turns the velocity of cell `neighbour` into the velocity the same volume
    /// flow has in cell `cell`, when the face between them is a step; 1 otherwise.
    pub(super) fn neighbour_scale(&self, cell: usize, neighbour: usize, mesh: &Mesh1d) -> f64 {
        let face = cell.max(neighbour);
        match self.index(face) {
            Some(k) if neighbour == cell + 1 => self.area_after[k] / mesh.face_areas()[face],
            Some(k) if neighbour + 1 == cell => mesh.face_areas()[face] / self.area_after[k],
            _ => 1.0,
        }
    }
    /// Area and flux the cell whose left face is `face` sees, when that face is a step.
    pub(super) fn entering(&self, face: usize) -> Option<(f64, FaceFlux)> {
        self.index(face)
            .map(|k| (self.area_after[k], self.after[k]))
    }
    /// Area the cell whose left face is `face` sees: `area` when that face is not a step.
    pub(super) fn left_area(&self, face: usize, area: f64) -> f64 {
        self.index(face).map_or(area, |k| self.area_after[k])
    }
    /// The flux the cell on the left of `face` receives when `face` is a step (the other side's
    /// is stored for [`Self::entering`]); `None` for every other face.
    pub(super) fn solve(
        &mut self,
        gas: &FrozenGas,
        face: usize,
        left: Primitive,
        right: Primitive,
        area_before: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<FaceFlux>, Error> {
        let Some(k) = self.index(face) else {
            return Ok(None);
        };
        let [before, after] =
            gas.step_fluxes(left, right, [area_before, self.area_after[k]], cancelled)?;
        self.after[k] = after;
        Ok(Some(before))
    }
    /// Axial force of each step annulus on the gas, `area_after * F_mom(after) - area_before *
    /// F_mom(before)`: the momentum the two neighbouring cells gain that no face flux carries.
    pub(super) fn impulses<'a>(
        &'a self,
        faces: &'a [FaceFlux],
        mesh: &'a Mesh1d,
    ) -> impl Iterator<Item = f64> + 'a {
        self.faces.iter().enumerate().map(move |(k, &face)| {
            self.area_after[k] * self.after[k].momentum
                - mesh.face_areas()[face] * faces[face].momentum
        })
    }
}

/// Signed unknowns of the stationary junction solve, in the frame where gas flows from the
/// upstream trace toward the downstream one.
struct Junction {
    mass_flow: f64,
    /// The mass-flow residual the solve accepted as zero [kg/s].
    mass_tolerance: f64,
    upstream: Primitive,
    downstream: Primitive,
    upstream_signal: f64,
    downstream_signal: f64,
}

impl FrozenGas {
    /// The two face fluxes of a junction face: `[left, right]`, each per unit of the area its
    /// own cell sees (`areas = [area_left, area_right]`). `left` and `right` are the
    /// reconstructed traces of the two cells at the face.
    fn step_fluxes(
        &self,
        left: Primitive,
        right: Primitive,
        areas: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[FaceFlux; 2], Error> {
        poll(cancelled)?;
        let mirror = |p: Primitive| Primitive {
            velocity: -p.velocity,
            ..p
        };
        // The acoustic solution decides which side donates; a root that disagrees with it
        // (only possible within the solver tolerance of rest) is retried the other way.
        let guess = acoustic_guess(left, right, areas);
        let forward_first = areas[0] * guess.1 >= 0.0;
        let mut last = Error::Invalid("radius-step junction");
        for attempt in 0..2 {
            let forward = forward_first == (attempt == 0);
            let solved = if forward {
                self.solve_junction(left, right, areas, cancelled)
            } else {
                self.solve_junction(mirror(right), mirror(left), [areas[1], areas[0]], cancelled)
            };
            match solved {
                // A root of the other sign within the solver tolerance of rest is the same flux.
                Ok(junction) if junction.mass_flow >= -junction.mass_tolerance => {
                    let upstream =
                        self.junction_flux(junction.upstream, junction.upstream_signal)?;
                    let mut downstream =
                        self.junction_flux(junction.downstream, junction.downstream_signal)?;
                    let (up_area, down_area) = if forward {
                        (areas[0], areas[1])
                    } else {
                        (areas[1], areas[0])
                    };
                    // Mass and energy leave one side exactly as they enter the other.
                    downstream.mass = upstream.mass * up_area / down_area;
                    downstream.energy = upstream.energy * up_area / down_area;
                    return Ok(if forward {
                        [upstream, downstream]
                    } else {
                        // Back to the real frame: mass and energy flow along -x.
                        let flip = |f: FaceFlux| FaceFlux {
                            mass: -f.mass,
                            energy: -f.energy,
                            ..f
                        };
                        [flip(downstream), flip(upstream)]
                    });
                }
                Ok(_) => last = Error::Invalid("radius-step junction flow direction"),
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(error) => last = error,
            }
        }
        Err(last)
    }

    fn junction_flux(&self, state: Primitive, signal: f64) -> Result<FaceFlux, Error> {
        let (q, p) = self.state_at(state.density, state.velocity, state.pressure)?;
        let flux = FaceFlux {
            mass: q.momentum,
            momentum: q.momentum * state.velocity + state.pressure,
            energy: (q.energy + state.pressure) * state.velocity,
            signal_speed: signal.max(p.velocity.abs() + p.sound_speed),
            used_hll_fallback: false,
        };
        if ![flux.mass, flux.momentum, flux.energy, flux.signal_speed]
            .iter()
            .all(|x| x.is_finite())
            || !positive(flux.signal_speed)
        {
            return Err(Error::Invalid("radius-step junction flux"));
        }
        Ok(flux)
    }

    /// Gas flows from `up` (area `areas[0]`) to `down` (area `areas[1]`).
    fn solve_junction(
        &self,
        up: Primitive,
        down: Primitive,
        areas: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Junction, Error> {
        let [a1, a2] = areas;
        let r = self.gas_constant();
        let z_down = down.density * down.sound_speed;
        // Unknowns x = (p1, u2): the pressure just upstream of the step and the velocity just
        // downstream of it. The chord matrix is the acoustic Jacobian of the residuals
        // (mass, velocity continuity) in those unknowns; the iteration contracts by O(Mach).
        let (p1_start, u2_start) = acoustic_guess(up, down, areas);
        let evaluate =
            |p1: f64, u2: f64, cancelled: &mut _| -> Result<(f64, f64, Evaluated), Error> {
                let wave_up = self.end_wave(up, 1.0, p1, &mut *cancelled)?;
                let s1 = wave_up.primitive;
                let mass_flow = s1.density * s1.velocity * a1;
                let dh = 0.5 * (s1.velocity * s1.velocity - u2 * u2);
                let t2 = self.shifted_temperature(s1.temperature, dh, &mut *cancelled)?;
                let exponent =
                    self.entropy_change(s1.temperature, p1, t2, p1, &mut *cancelled)? / r;
                let p2 = p1 + p1 * det::expm1(exponent);
                let rho2 = p2 / (r * t2);
                let wave_down = self.end_wave(down, -1.0, p2, &mut *cancelled)?;
                Ok((
                    rho2 * u2 * a2 - mass_flow,
                    u2 - wave_down.primitive.velocity,
                    Evaluated {
                        s1,
                        p2,
                        rho2,
                        signal_up: wave_up.signal_speed,
                        signal_down: wave_down.signal_speed,
                    },
                ))
            };
        let scale = up.velocity.abs() + down.velocity.abs();
        let tolerance_u = 1e-10 * scale + 1e-12;
        let tolerance_m = tolerance_u * up.density * a1.max(a2);
        let (mut p1, mut u2) = (p1_start, u2_start);
        let mut best = f64::INFINITY;
        let mut stalled = 0;
        for _ in 0..100 {
            poll(cancelled)?;
            let (ra, rb, e) = evaluate(p1, u2, cancelled)?;
            let size = (ra.abs() / tolerance_m).max(rb.abs() / tolerance_u);
            if size <= 1.0 {
                let upstream = e.s1;
                let downstream = self.from_rho_u_p(e.rho2, u2, e.p2)?;
                for state in [upstream, downstream] {
                    if state.velocity.abs() >= state.sound_speed * (1.0 - 1e-3) {
                        return Err(Error::Invalid("radius-step junction is choked"));
                    }
                }
                return Ok(Junction {
                    mass_flow: upstream.density * upstream.velocity * a1,
                    mass_tolerance: tolerance_m,
                    upstream,
                    downstream,
                    upstream_signal: e.signal_up,
                    downstream_signal: e.signal_down,
                });
            }
            if size < best {
                best = size;
                stalled = 0;
            } else {
                stalled += 1;
                if stalled >= 6 {
                    break;
                }
            }
            // Chord step: J0 = [[a1/c, rho a2], [-1/z_down, 1]] acting on (dp1, du2).
            let (m00, m01, m10, m11) = (a1 / up.sound_speed, up.density * a2, -1.0 / z_down, 1.0);
            let determinant = m00 * m11 - m01 * m10;
            let dp1 = -(m11 * ra - m01 * rb) / determinant;
            let du2 = -(-m10 * ra + m00 * rb) / determinant;
            p1 += dp1;
            u2 += du2;
            if !positive(p1) || !u2.is_finite() {
                return Err(Error::Invalid(
                    "radius-step junction left the source domain",
                ));
            }
        }
        Err(Error::Invalid("radius-step junction did not converge"))
    }

    /// Temperature T with `h(T) - h(t1) = dh` along the source's own heat capacity.
    fn shifted_temperature(
        &self,
        t1: f64,
        dh: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        if dh == 0.0 {
            return Ok(t1);
        }
        let r = self.gas_constant();
        let domain = self.domain();
        let mut t = t1 + dh / CaloricLaw::heat_capacities(self, t1, 1.0, cancelled)?.cp;
        for _ in 0..16 {
            poll(cancelled)?;
            if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&t) {
                return Err(Error::CaloricPath("radius-step junction temperature"));
            }
            let change = self.internal_energy_change(t1, t, cancelled)? + r * (t - t1);
            let cp = CaloricLaw::heat_capacities(self, t, 1.0, cancelled)?.cp;
            let correction = (dh - change) / cp;
            t += correction;
            if correction.abs() <= 4.0 * f64::EPSILON * t.abs() {
                return Ok(t);
            }
        }
        Err(Error::TemperatureIterations)
    }
}

struct Evaluated {
    s1: Primitive,
    p2: f64,
    rho2: f64,
    signal_up: f64,
    signal_down: f64,
}

/// The linear acoustic junction (equal pressure, equal volume flow) of two traces:
/// `(p, u_downstream)` for the gas flowing from `up` to `down`.
fn acoustic_guess(up: Primitive, down: Primitive, areas: [f64; 2]) -> (f64, f64) {
    let [a1, a2] = areas;
    let z_up = up.density * up.sound_speed;
    let z_down = down.density * down.sound_speed;
    // A1 (u_up - dp/z_up) = A2 (u_down + (dp + p_up - p_down)/z_down), dp = p - p_up.
    let dp = (a1 * up.velocity - a2 * down.velocity + a2 * (down.pressure - up.pressure) / z_down)
        / (a1 / z_up + a2 / z_down);
    (
        up.pressure + dp,
        down.velocity + (up.pressure + dp - down.pressure) / z_down,
    )
}

#[cfg(test)]
mod tests;
