use super::{
    Conserved, Error, FaceFlux, FrozenGas, MAX_EULER_CELLS, Primitive, allocated, poll, positive,
};
use fs_qty::{Length, Time};

/// Physical boundary on one end of a fixed-area interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    /// Zero-gradient ghost state. This is not a nonreflecting subsonic outlet.
    Transmissive,
    /// Stationary impermeable adiabatic wall, with exactly zero mass/energy flux.
    Reflecting,
}

/// End conditions; a periodic face is computed once and shared bit-for-bit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundaries {
    /// Both ends are the same periodic interface.
    Periodic,
    /// Independent physical left and right boundaries.
    Ends {
        /// Boundary at the smallest coordinate.
        left: Boundary,
        /// Boundary at the largest coordinate.
        right: Boundary,
    },
}

/// Spatial reconstruction admitted by this engineering baseline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reconstruction {
    /// Piecewise-constant conservative finite volumes.
    FirstOrder,
    /// Monotonized-central primitive slopes, reduced together to retain the source domain.
    Muscl,
}

/// Explicit per-step stability and retry budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepControl {
    /// Maximum face-speed CFL, in (0, 0.45]. Checked at both SSP RK stages.
    pub cfl: f64,
    /// Spatial order before an admissibility-triggered first-order retry.
    pub reconstruction: Reconstruction,
    /// Maximum halvings after the first attempt, in 0..=32.
    pub maximum_halvings: usize,
}
impl Default for StepControl {
    fn default() -> Self {
        Self {
            cfl: 0.4,
            reconstruction: Reconstruction::Muscl,
            maximum_halvings: 16,
        }
    }
}

/// Extensive integrals over the actual cross-sectional area and cell widths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Totals {
    /// Mass [kg].
    pub mass: f64,
    /// Axial momentum [kg m/s].
    pub momentum: f64,
    /// Referenced total energy [J]. Physical energy adds `gas.energy_reference()*mass`.
    pub energy: f64,
}
impl Totals {
    fn from_array(v: [f64; 3]) -> Self {
        Self {
            mass: v[0],
            momentum: v[1],
            energy: v[2],
        }
    }
    fn array(self) -> [f64; 3] {
        [self.mass, self.momentum, self.energy]
    }
    /// Recover the physical total energy with the original signed source reference.
    pub fn absolute_energy(self, gas: &FrozenGas) -> f64 {
        self.energy + gas.energy_reference() * self.mass
    }
}

/// Accepted single-step result. No stage or failed retry is published.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    /// Actual elapsed time, at most the requested maximum [s].
    pub elapsed: Time,
    /// Maximum face-speed CFL across the two accepted stages.
    pub maximum_cfl: f64,
    /// Number of rejected and halved attempts before acceptance.
    pub halvings: usize,
    /// Number of accepted-stage interfaces that used HLL fallback.
    pub hll_fallback_faces: usize,
    /// Number of accepted-stage cells whose slope was reduced for source admission.
    pub limited_cells: usize,
    /// Whether retry changed the whole step to first-order spatial reconstruction.
    pub first_order_retry: bool,
    /// Original extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Area- and time-integrated left inflow minus right outflow, using both RK stages.
    pub boundary_change: Totals,
    /// `after - before - boundary_change` in the same physical units.
    pub conservation_residual: Totals,
}

struct Workspace {
    primitive: Vec<Primitive>,
    slopes: Vec<[f64; 3]>,
    faces: Vec<FaceFlux>,
    stage: Vec<Conserved>,
    candidate: Vec<Conserved>,
}

/// Uniform one-dimensional finite-volume field of one frozen source-owned mixture.
///
/// Evolution conserves mass, axial momentum and total energy through shared face
/// fluxes. Spatial reconstruction is MUSCL; time integration is SSP RK2. Accepted
/// stages remain inside the declared NASA-9 temperature path, with positive density
/// and pressure. Refusal/cancellation preserves every cell and the original time.
/// Workspace is allocated once at construction, not in a successful step.
pub struct Euler1d {
    gas: FrozenGas,
    states: Vec<Conserved>,
    cell_width: f64,
    cross_section: f64,
    boundaries: Boundaries,
    time: f64,
    work: Workspace,
}

impl Euler1d {
    /// Construct a fixed-area mesh from validated conservative cell averages.
    /// `cross_section` is in m²; coordinates are `i*cell_width` plus a caller-owned origin.
    pub fn new(
        gas: FrozenGas,
        cell_width: Length,
        cross_section: f64,
        states: Vec<Conserved>,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let n = states.len();
        if !(2..=MAX_EULER_CELLS).contains(&n)
            || !positive(cell_width.value())
            || !positive(cross_section)
            || !positive(cell_width.value() * cross_section)
            || !(cell_width.value() * n as f64).is_finite()
        {
            return Err(Error::Invalid("mesh size, cell width or cross section"));
        }
        let mut primitive = allocated(n, Primitive::default())?;
        for (i, &state) in states.iter().enumerate() {
            primitive[i] = gas.decode(state, None, &mut cancelled)?;
        }
        let work = Workspace {
            primitive,
            slopes: allocated(n, [0.0; 3])?,
            faces: allocated(n + 1, FaceFlux::default())?,
            stage: allocated(n, Conserved::default())?,
            candidate: allocated(n, Conserved::default())?,
        };
        poll(&mut cancelled)?;
        Ok(Self {
            gas,
            states,
            cell_width: cell_width.value(),
            cross_section,
            boundaries,
            time: 0.0,
            work,
        })
    }
    /// Original frozen-mixture EOS and source.
    pub const fn gas(&self) -> &FrozenGas {
        &self.gas
    }
    /// Current conservative cell averages, all sharing the owner's energy reference.
    pub fn states(&self) -> &[Conserved] {
        &self.states
    }
    /// Accepted physical time.
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Mesh cell width [m].
    pub fn cell_width(&self) -> Length {
        Length::new(self.cell_width)
    }
    /// Uniform cross-sectional area [m²].
    pub const fn cross_section(&self) -> f64 {
        self.cross_section
    }
    /// Extensive conserved state, accumulated in canonical cell order with compensation.
    pub fn totals(&self) -> Totals {
        totals(&self.states, self.cell_width * self.cross_section)
    }

    /// Advance by one accepted step no greater than `maximum_duration`.
    ///
    /// Initial dt uses cell acoustic speeds; both stages check the actual face
    /// speeds. An inadmissible trial halves dt and uses first-order reconstruction.
    /// No density, pressure, temperature or energy is clipped. Exhausted retries,
    /// invalid source states, and cancellation leave the physical field unchanged.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("time step, CFL or retry controls"));
        }
        let mut speed = 0.0_f64;
        for (i, &state) in self.states.iter().enumerate() {
            let p = self.gas.decode(
                state,
                Some(self.work.primitive[i].temperature),
                &mut cancelled,
            )?;
            self.work.primitive[i] = p;
            speed = speed.max(p.velocity.abs() + p.sound_speed);
        }
        // Reserve ten percent for reconstructed/stage signal-speed growth;
        // the actual face CFL still decides admission at both stages.
        let mut dt = maximum_duration
            .value()
            .min((0.9 * controls.cfl * self.cell_width) / speed);
        let before = self.totals();
        let mut last = Error::Invalid("no admissible time increment");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(Error::Invalid(
                    "time step cannot advance the represented clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            match attempt(
                &self.gas,
                &self.states,
                self.cell_width,
                self.boundaries,
                dt,
                controls.cfl,
                reconstruction,
                &mut self.work,
                &mut cancelled,
            ) {
                Ok((a, b)) => {
                    let after = totals(&self.work.candidate, self.cell_width * self.cross_section);
                    let factor = 0.5 * dt * self.cross_section;
                    let boundary = Totals::from_array(core::array::from_fn(|k| {
                        factor
                            * ((a.left.array()[k] - a.right.array()[k])
                                + (b.left.array()[k] - b.right.array()[k]))
                    }));
                    let residual = Totals::from_array(core::array::from_fn(|k| {
                        after.array()[k] - before.array()[k] - boundary.array()[k]
                    }));
                    if !after
                        .array()
                        .iter()
                        .chain(&boundary.array())
                        .chain(&residual.array())
                        .all(|x| x.is_finite())
                    {
                        return Err(Error::Invalid("non-finite extensive balance"));
                    }
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    self.time += dt;
                    return Ok(StepReport {
                        elapsed: Time::new(dt),
                        maximum_cfl: a.speed.max(b.speed) * dt / self.cell_width,
                        halvings,
                        hll_fallback_faces: a.fallbacks + b.fallbacks,
                        limited_cells: a.limited + b.limited,
                        first_order_retry: first_order_retry
                            && controls.reconstruction == Reconstruction::Muscl,
                        before,
                        after,
                        boundary_change: boundary,
                        conservation_residual: residual,
                    });
                }
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(error) => {
                    // A CFL-only refusal changes dt, not the spatial scheme.
                    // Smooth reconstructed waves must retain MUSCL accuracy.
                    if !matches!(error, Error::CflExceeded { .. }) {
                        first_order_retry = true;
                    }
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(Error::StepRetries {
            attempts: controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}

#[derive(Clone, Copy)]
struct Stage {
    left: FaceFlux,
    right: FaceFlux,
    speed: f64,
    fallbacks: usize,
    limited: usize,
}

#[allow(clippy::too_many_arguments)]
fn attempt(
    gas: &FrozenGas,
    states: &[Conserved],
    dx: f64,
    boundaries: Boundaries,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    work: &mut Workspace,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage), Error> {
    let a = operator(
        gas,
        states,
        boundaries,
        reconstruction,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        cancelled,
    )?;
    let ratio = dt / dx;
    if a.speed * ratio > cfl * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded {
            observed: a.speed * ratio,
            maximum: cfl,
        });
    }
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        work.stage[i] = Conserved::from_array(core::array::from_fn(|k| {
            q.array()[k] - ratio * (work.faces[i + 1].array()[k] - work.faces[i].array()[k])
        }));
        if !gas.admissible_state(work.stage[i]) {
            return Err(Error::Invalid("first RK stage outside source domain"));
        }
    }
    let b = operator(
        gas,
        &work.stage,
        boundaries,
        reconstruction,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        cancelled,
    )?;
    if b.speed * ratio > cfl * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded {
            observed: b.speed * ratio,
            maximum: cfl,
        });
    }
    for (i, &q) in work.stage.iter().enumerate() {
        poll(cancelled)?;
        let second = Conserved::from_array(core::array::from_fn(|k| {
            q.array()[k] - ratio * (work.faces[i + 1].array()[k] - work.faces[i].array()[k])
        }));
        if !gas.admissible_state(second) {
            return Err(Error::Invalid("second RK stage outside source domain"));
        }
        let candidate = Conserved::from_array(core::array::from_fn(|k| {
            0.5 * states[i].array()[k] + 0.5 * second.array()[k]
        }));
        // Final EOS evaluation, including source finite-arithmetic checks, before publication.
        gas.decode(candidate, Some(work.primitive[i].temperature), cancelled)?;
        work.candidate[i] = candidate;
    }
    Ok((a, b))
}

#[allow(clippy::too_many_arguments)]
fn operator(
    gas: &FrozenGas,
    states: &[Conserved],
    boundaries: Boundaries,
    reconstruction: Reconstruction,
    primitive: &mut [Primitive],
    slopes: &mut [[f64; 3]],
    faces: &mut [FaceFlux],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, Error> {
    let n = states.len();
    for (i, &q) in states.iter().enumerate() {
        primitive[i] = gas.decode(q, Some(primitive[i].temperature), cancelled)?;
    }
    let values = |p: Primitive| [p.density, p.velocity, p.pressure];
    let mut limited = 0;
    for i in 0..n {
        poll(cancelled)?;
        slopes[i] = [0.0; 3];
        if reconstruction == Reconstruction::FirstOrder {
            continue;
        }
        let (prev, next) = if boundaries == Boundaries::Periodic {
            ((i + n - 1) % n, (i + 1) % n)
        } else if i == 0 || i + 1 == n {
            continue;
        } else {
            (i - 1, i + 1)
        };
        let (a, b, c) = (
            values(primitive[prev]),
            values(primitive[i]),
            values(primitive[next]),
        );
        let slope = core::array::from_fn(|k| {
            minmod(
                2.0 * (b[k] - a[k]),
                0.5 * (c[k] - a[k]),
                2.0 * (c[k] - b[k]),
            )
        });
        let mut scale = 1.0;
        for reduction in 0..=12 {
            let valid = [-0.5, 0.5].into_iter().all(|side| {
                let rho = b[0] + side * scale * slope[0];
                let u = b[1] + side * scale * slope[1];
                let p = b[2] + side * scale * slope[2];
                let t = (p / rho) / gas.gas_constant();
                positive(rho)
                    && positive(p)
                    && u.is_finite()
                    && (gas.domain().minimum_temperature..=gas.domain().maximum_temperature)
                        .contains(&t)
            });
            if valid {
                slopes[i] = slope.map(|x| scale * x);
                break;
            }
            if reduction == 0 {
                limited += 1;
            }
            scale *= 0.5;
        }
    }
    let reconstruct = |i: usize, sign: f64| {
        if slopes[i] == [0.0; 3] {
            return Ok((states[i], primitive[i]));
        }
        let p = values(primitive[i]);
        gas.state_at(
            p[0] + sign * slopes[i][0],
            p[1] + sign * slopes[i][1],
            p[2] + sign * slopes[i][2],
        )
    };
    let mut speed = 0.0_f64;
    let mut fallbacks = 0;
    for face in 0..=n {
        poll(cancelled)?;
        if face == n && boundaries == Boundaries::Periodic {
            faces[n] = faces[0];
            continue;
        }
        let (left, right, wall) = if face == 0 || face == n {
            match boundaries {
                Boundaries::Periodic => (reconstruct(n - 1, 0.5)?, reconstruct(0, -0.5)?, false),
                Boundaries::Ends { left, right } => {
                    let inner = if face == 0 {
                        reconstruct(0, -0.5)?
                    } else {
                        reconstruct(n - 1, 0.5)?
                    };
                    let boundary = if face == 0 { left } else { right };
                    let mut ghost = inner;
                    if boundary == Boundary::Reflecting {
                        ghost.0.momentum = -ghost.0.momentum;
                        ghost.1.velocity = -ghost.1.velocity;
                    }
                    if face == 0 {
                        (ghost, inner, boundary == Boundary::Reflecting)
                    } else {
                        (inner, ghost, boundary == Boundary::Reflecting)
                    }
                }
            }
        } else {
            (reconstruct(face - 1, 0.5)?, reconstruct(face, -0.5)?, false)
        };
        let mut flux = gas.interface(left.0, left.1, right.0, right.1, cancelled)?;
        if wall {
            flux.mass = 0.0;
            flux.energy = 0.0;
        }
        speed = speed.max(flux.signal_speed);
        fallbacks += usize::from(flux.used_hll_fallback);
        faces[face] = flux;
    }
    Ok(Stage {
        left: faces[0],
        right: faces[n],
        speed,
        fallbacks,
        limited,
    })
}

fn minmod(a: f64, b: f64, c: f64) -> f64 {
    if a > 0.0 && b > 0.0 && c > 0.0 {
        a.min(b).min(c)
    } else if a < 0.0 && b < 0.0 && c < 0.0 {
        a.max(b).max(c)
    } else {
        0.0
    }
}
fn totals(states: &[Conserved], volume: f64) -> Totals {
    let mut sum = [0.0; 3];
    let mut compensation = [0.0; 3];
    for q in states {
        for k in 0..3 {
            let y = q.array()[k] * volume - compensation[k];
            let next = sum[k] + y;
            compensation[k] = (next - sum[k]) - y;
            sum[k] = next;
        }
    }
    Totals::from_array(sum)
}
