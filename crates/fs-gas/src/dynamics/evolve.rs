use super::{Conserved, Error, FaceFlux, FrozenGas, Mesh1d, Primitive, allocated, poll, positive};
use fs_qty::{Length, Time};

/// Physical boundary on one end of a stationary interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Boundary {
    /// Zero-gradient ghost state. This is not a nonreflecting subsonic outlet.
    Transmissive,
    /// Stationary impermeable adiabatic wall, with exactly zero mass/energy flux.
    Reflecting,
    /// Fixed, stationary reservoir state in this gas owner's energy reference.
    /// Its pressure and temperature are stagnation conditions for inflow and
    /// its pressure is the backpressure for outflow. Momentum must be zero.
    /// A source-aware characteristic interface supplies the shared flux;
    /// this is not a zero-gradient or nonreflecting boundary.
    Reservoir(Conserved),
}

/// End conditions; a periodic face is computed once and shared bit-for-bit.
#[derive(Clone, Copy, Debug, PartialEq)]
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
    pub(super) fn from_array(v: [f64; 3]) -> Self {
        Self {
            mass: v[0],
            momentum: v[1],
            energy: v[2],
        }
    }
    pub(super) fn array(self) -> [f64; 3] {
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
    /// Accepted left-face exchange, signed into this interval, using both RK stages.
    pub left_boundary_change: Totals,
    /// Accepted right-face exchange, signed into this interval, using both RK stages.
    pub right_boundary_change: Totals,
    /// Integrated p*dA/dx force on the gas [kg m/s]. Stationary walls do no work.
    pub area_momentum_impulse: f64,
    /// `after - before - boundary_change - area momentum impulse`.
    pub conservation_residual: Totals,
}

/// One SSP RK stage's end cells, offered to caller-owned end ports.
///
/// Boundary cells of a non-periodic interval are never reconstructed, so each
/// trace is exactly its cell average, resolved through the original source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortStage {
    /// 0: the predictor, on the accepted field at the step's start time;
    /// 1: the corrector, on the predicted field at the step's end time.
    pub stage: usize,
    /// Physical time the offered field represents [s].
    pub time: f64,
    /// Trial step duration [s]. A halved retry offers a smaller value.
    pub dt: f64,
    /// Cell at the smallest coordinate and its resolved primitive.
    pub left: (Conserved, Primitive),
    /// Cell at the largest coordinate and its resolved primitive.
    pub right: (Conserved, Primitive),
}

type NoPorts = fn(&PortStage) -> Result<[Option<FaceFlux>; 2], Error>;

pub(super) struct Workspace {
    pub(super) primitive: Vec<Primitive>,
    pub(super) slopes: Vec<[f64; 3]>,
    pub(super) faces: Vec<FaceFlux>,
    pub(super) stage: Vec<Conserved>,
    pub(super) candidate: Vec<Conserved>,
}

/// Quasi-one-dimensional finite-volume field of one frozen source-owned mixture.
///
/// Evolution conserves mass, axial momentum and total energy through shared face
/// fluxes. Spatial reconstruction is MUSCL; time integration is SSP RK2. Accepted
/// stages remain inside the declared NASA-9 temperature path, with positive density
/// and pressure. Refusal/cancellation preserves every cell and the original time.
/// Workspace is allocated once at construction, not in a successful step.
pub struct Euler1d {
    pub(super) gas: FrozenGas,
    pub(super) states: Vec<Conserved>,
    pub(super) mesh: Mesh1d,
    pub(super) boundaries: Boundaries,
    pub(super) time: f64,
    pub(super) work: Workspace,
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
        let mesh = Mesh1d::uniform(states.len(), cell_width, cross_section)?;
        Self::new_on_mesh(gas, mesh, states, boundaries, cancelled)
    }

    /// Construct on a stationary, possibly nonuniform and tapered mesh.
    /// Face areas and integrated cell volumes enter the same conservative
    /// operator; periodic ends must have identical physical area.
    pub fn new_on_mesh(
        gas: FrozenGas,
        mesh: Mesh1d,
        states: Vec<Conserved>,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let n = states.len();
        if n != mesh.cell_count() {
            return Err(Error::Invalid("state count differs from physical mesh"));
        }
        if boundaries == Boundaries::Periodic && mesh.face_areas()[0] != mesh.face_areas()[n] {
            return Err(Error::Invalid(
                "periodic end faces must have identical physical area",
            ));
        }
        validate_reservoirs(&gas, boundaries, &mut cancelled)?;
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
            mesh,
            boundaries,
            time: 0.0,
            work,
        })
    }
    /// Original frozen-mixture EOS and source.
    pub const fn gas(&self) -> &FrozenGas {
        &self.gas
    }
    /// Physical end conditions used by this interval.
    pub const fn boundaries(&self) -> Boundaries {
        self.boundaries
    }
    /// Replace prescribed physical end conditions without advancing the field.
    /// Reservoirs remain fixed throughout each subsequent RK step. A caller
    /// may update them between accepted steps to supply a prescribed drive.
    /// Invalid source states or cancellation retain the previous boundaries.
    pub fn set_boundaries(
        &mut self,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(&mut cancelled)?;
        let areas = self.mesh.face_areas();
        if boundaries == Boundaries::Periodic && areas[0] != areas[self.states.len()] {
            return Err(Error::Invalid(
                "periodic end faces must have identical physical area",
            ));
        }
        validate_reservoirs(&self.gas, boundaries, &mut cancelled)?;
        poll(&mut cancelled)?;
        self.boundaries = boundaries;
        Ok(())
    }
    /// Current conservative cell averages, all sharing the owner's energy reference.
    pub fn states(&self) -> &[Conserved] {
        &self.states
    }
    /// Accepted physical time.
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Complete stationary geometry and integrated physical volumes.
    pub const fn mesh(&self) -> &Mesh1d {
        &self.mesh
    }
    /// Uniform mesh cell width [m], or None on a nonuniform mesh.
    pub fn cell_width(&self) -> Option<Length> {
        self.mesh.uniform_cell_width()
    }
    /// Uniform cross-sectional area [m²], or None in a variable-area mesh.
    pub const fn cross_section(&self) -> Option<f64> {
        self.mesh.uniform_cross_section()
    }
    /// Extensive conserved state, accumulated in canonical cell order with compensation.
    pub fn totals(&self) -> Totals {
        totals(&self.states, self.mesh.volumes(), &mut || false)
            .expect("an infallible cancellation callback cannot refuse")
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        self.advance(maximum_duration, controls, None::<&mut NoPorts>, cancelled)
    }

    /// [`Self::step`] whose end faces may be supplied by caller-owned ports.
    ///
    /// `ports` is called once per RK stage, before that stage's fluxes, with
    /// both end cells. Each `Some` flux replaces that end's declared boundary
    /// for the stage; `None` keeps it. Fluxes are per unit face area, along the
    /// increasing coordinate, like every interior face. The accepted step's end
    /// exchange is the stage average reported in [`StepReport`]. A port error
    /// is a refusal of the trial: it is retried like any other and, when the
    /// retries are exhausted, leaves the field and time unchanged. A port that
    /// keeps its own memory must therefore commit it only after acceptance.
    pub fn step_with_ports(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut ports: impl FnMut(&PortStage) -> Result<[Option<FaceFlux>; 2], Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        if self.boundaries == Boundaries::Periodic {
            return Err(Error::Invalid("end ports need a non-periodic interval"));
        }
        self.advance(maximum_duration, controls, Some(&mut ports), cancelled)
    }

    /// Replace every cell average through a local, caller-owned source update
    /// between accepted steps (wall friction, wall heat exchange).
    ///
    /// `update(cell, average, primitive)` returns the new average. Every
    /// candidate is resolved through the original caloric law before any is
    /// published; a refusal or cancellation leaves the field unchanged. Time
    /// does not advance: the caller owns the duration the sources represent.
    pub fn apply_cell_sources(
        &mut self,
        mut update: impl FnMut(usize, Conserved, Primitive) -> Result<Conserved, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(&mut cancelled)?;
        for (i, &q) in self.states.iter().enumerate() {
            let p = self
                .gas
                .decode(q, Some(self.work.primitive[i].temperature), &mut cancelled)?;
            let next = update(i, q, p)?;
            if !self.gas.admissible_state(next) {
                return Err(Error::Invalid("cell source left the source domain"));
            }
            self.work.primitive[i] = self.gas.decode(next, Some(p.temperature), &mut cancelled)?;
            self.work.candidate[i] = next;
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.states, &mut self.work.candidate);
        Ok(())
    }

    fn advance<P>(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut ports: Option<&mut P>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error>
    where
        P: FnMut(&PortStage) -> Result<[Option<FaceFlux>; 2], Error>,
    {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("time step, CFL or retry controls"));
        }
        let mut signal_rate = 0.0_f64;
        for (i, &state) in self.states.iter().enumerate() {
            let p = self.gas.decode(
                state,
                Some(self.work.primitive[i].temperature),
                &mut cancelled,
            )?;
            self.work.primitive[i] = p;
            let area = 0.5 * self.mesh.face_areas()[i] + 0.5 * self.mesh.face_areas()[i + 1];
            let rate = (p.velocity.abs() + p.sound_speed) * (area / self.mesh.volumes()[i]);
            signal_rate = signal_rate.max(rate);
        }
        // Reserve ten percent for reconstructed/stage signal-speed growth;
        // the actual face CFL still decides admission at both stages.
        let mut dt = maximum_duration
            .value()
            .min((0.9 * controls.cfl) / signal_rate);
        let before = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
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
                &self.mesh,
                self.boundaries,
                [self.time, dt],
                controls.cfl,
                reconstruction,
                &mut self.work,
                ports.as_deref_mut(),
                &mut cancelled,
            ) {
                Ok((a, b)) => {
                    let after = totals(&self.work.candidate, self.mesh.volumes(), &mut cancelled)?;
                    let factor = 0.5 * dt;
                    let left_area = self.mesh.face_areas()[0];
                    let right_area = self.mesh.face_areas()[self.states.len()];
                    let boundary = Totals::from_array(core::array::from_fn(|k| {
                        factor
                            * (left_area * (a.left.array()[k] + b.left.array()[k])
                                - right_area * (a.right.array()[k] + b.right.array()[k]))
                    }));
                    let left_boundary_change = Totals::from_array(core::array::from_fn(|k| {
                        factor * left_area * (a.left.array()[k] + b.left.array()[k])
                    }));
                    let right_boundary_change = Totals::from_array(core::array::from_fn(|k| {
                        -factor * right_area * (a.right.array()[k] + b.right.array()[k])
                    }));
                    let area_momentum_impulse = factor * (a.area_force + b.area_force);
                    let residual = Totals::from_array(core::array::from_fn(|k| {
                        after.array()[k]
                            - before.array()[k]
                            - boundary.array()[k]
                            - if k == 1 { area_momentum_impulse } else { 0.0 }
                    }));
                    if !after
                        .array()
                        .iter()
                        .chain(&boundary.array())
                        .chain(&left_boundary_change.array())
                        .chain(&right_boundary_change.array())
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
                        maximum_cfl: a.signal_rate.max(b.signal_rate) * dt,
                        halvings,
                        hll_fallback_faces: a.fallbacks + b.fallbacks,
                        limited_cells: a.limited + b.limited,
                        first_order_retry: first_order_retry
                            && controls.reconstruction == Reconstruction::Muscl,
                        before,
                        after,
                        boundary_change: boundary,
                        left_boundary_change,
                        right_boundary_change,
                        area_momentum_impulse,
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
pub(super) struct Stage {
    pub(super) left: FaceFlux,
    pub(super) right: FaceFlux,
    pub(super) signal_rate: f64,
    pub(super) area_force: f64,
    pub(super) fallbacks: usize,
    pub(super) limited: usize,
}

#[allow(clippy::too_many_arguments)]
fn attempt<P>(
    gas: &FrozenGas,
    states: &[Conserved],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    [time, dt]: [f64; 2],
    cfl: f64,
    reconstruction: Reconstruction,
    work: &mut Workspace,
    mut ports: Option<&mut P>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage), Error>
where
    P: FnMut(&PortStage) -> Result<[Option<FaceFlux>; 2], Error>,
{
    let ends = offer(
        gas,
        states,
        [time, dt],
        0,
        &work.primitive,
        ports.as_deref_mut(),
        cancelled,
    )?;
    let a = operator(
        gas,
        states,
        mesh,
        boundaries,
        reconstruction,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        ends,
        cancelled,
    )?;
    if a.signal_rate * dt > cfl * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded {
            observed: a.signal_rate * dt,
            maximum: cfl,
        });
    }
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        work.stage[i] = forward(q, i, dt, mesh, &work.primitive, &work.faces);
        if !gas.admissible_state(work.stage[i]) {
            return Err(Error::Invalid("first RK stage outside source domain"));
        }
    }
    let ends = offer(
        gas,
        &work.stage,
        [time + dt, dt],
        1,
        &work.primitive,
        ports,
        cancelled,
    )?;
    let b = operator(
        gas,
        &work.stage,
        mesh,
        boundaries,
        reconstruction,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        ends,
        cancelled,
    )?;
    if b.signal_rate * dt > cfl * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded {
            observed: b.signal_rate * dt,
            maximum: cfl,
        });
    }
    for (i, &q) in work.stage.iter().enumerate() {
        poll(cancelled)?;
        let second = forward(q, i, dt, mesh, &work.primitive, &work.faces);
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

// The same guess the stage operator uses, so an offered trace is the exact
// average and primitive that operator would pair with a declared boundary.
fn offer<P>(
    gas: &FrozenGas,
    states: &[Conserved],
    [time, dt]: [f64; 2],
    stage: usize,
    primitive: &[Primitive],
    ports: Option<&mut P>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Option<FaceFlux>; 2], Error>
where
    P: FnMut(&PortStage) -> Result<[Option<FaceFlux>; 2], Error>,
{
    let Some(port) = ports else {
        return Ok([None; 2]);
    };
    let n = states.len();
    let left = gas.decode(states[0], Some(primitive[0].temperature), cancelled)?;
    let right = gas.decode(states[n - 1], Some(primitive[n - 1].temperature), cancelled)?;
    poll(cancelled)?;
    let ends = port(&PortStage {
        stage,
        time,
        dt,
        left: (states[0], left),
        right: (states[n - 1], right),
    })?;
    for flux in ends.iter().flatten() {
        if ![flux.mass, flux.momentum, flux.energy]
            .iter()
            .all(|x| x.is_finite())
            || !positive(flux.signal_speed)
        {
            return Err(Error::Invalid("end port flux or signal speed"));
        }
    }
    Ok(ends)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn operator(
    gas: &FrozenGas,
    states: &[Conserved],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    reconstruction: Reconstruction,
    primitive: &mut [Primitive],
    slopes: &mut [[f64; 3]],
    faces: &mut [FaceFlux],
    endpoint_overrides: [Option<FaceFlux>; 2],
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
            let x = mesh.center(i);
            let length = mesh.positions()[n] - mesh.positions()[0];
            let xl = mesh.center(prev) - if i == 0 { length } else { 0.0 };
            let xr = mesh.center(next) + if i + 1 == n { length } else { 0.0 };
            let width = mesh.width(i);
            minmod(
                2.0 * (b[k] - a[k]) / (x - xl),
                (c[k] - a[k]) / (xr - xl),
                2.0 * (c[k] - b[k]) / (xr - x),
            ) * width
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
    let mut fallbacks = 0;
    for face in 0..=n {
        poll(cancelled)?;
        if face == 0 || face == n {
            if let Some(flux) = endpoint_overrides[usize::from(face == n)] {
                fallbacks += usize::from(flux.used_hll_fallback);
                faces[face] = flux;
                continue;
            }
            if let Boundaries::Ends { left, right } = boundaries {
                if let Boundary::Reservoir(reservoir) = if face == 0 { left } else { right } {
                    let inner = if face == 0 {
                        reconstruct(0, -0.5)?
                    } else {
                        reconstruct(n - 1, 0.5)?
                    };
                    let exterior = gas.decode(reservoir, None, cancelled)?;
                    faces[face] = gas.plenum_flux(
                        inner.0,
                        inner.1,
                        reservoir,
                        exterior,
                        if face == 0 { -1.0 } else { 1.0 },
                        cancelled,
                    )?;
                    continue;
                }
            }
        }
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
        let flux = if wall {
            let (inside, normal) = if face == 0 {
                (right, -1.0)
            } else {
                (left, 1.0)
            };
            gas.wall_flux(inside.0, inside.1, normal, cancelled)?
        } else {
            gas.interface(left.0, left.1, right.0, right.1, cancelled)?
        };
        fallbacks += usize::from(flux.used_hll_fallback);
        faces[face] = flux;
    }
    let mut signal_rate = 0.0_f64;
    let mut area_force = 0.0;
    let mut compensation = 0.0;
    for i in 0..n {
        poll(cancelled)?;
        let rate = (0.5 * mesh.face_areas()[i] * faces[i].signal_speed
            + 0.5 * mesh.face_areas()[i + 1] * faces[i + 1].signal_speed)
            / mesh.volumes()[i];
        if !positive(rate) {
            return Err(Error::Invalid("unrepresentable geometric face CFL"));
        }
        signal_rate = signal_rate.max(rate);
        let term = primitive[i].pressure * (mesh.face_areas()[i + 1] - mesh.face_areas()[i]);
        let y = term - compensation;
        let next = area_force + y;
        compensation = (next - area_force) - y;
        area_force = next;
    }
    Ok(Stage {
        left: faces[0],
        right: faces[n],
        signal_rate,
        area_force,
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
pub(super) fn forward(
    q: Conserved,
    i: usize,
    dt: f64,
    mesh: &Mesh1d,
    primitive: &[Primitive],
    faces: &[FaceFlux],
) -> Conserved {
    let left = faces[i].array();
    let right = faces[i + 1].array();
    let al = mesh.face_areas()[i];
    let ar = mesh.face_areas()[i + 1];
    let ratio = dt / mesh.volumes()[i];
    Conserved::from_array(core::array::from_fn(|k| {
        // Pair pressure traction with the geometric source before multiplying
        // by area: a uniform-pressure resting state then cancels exactly.
        let pressure = if k == 1 { primitive[i].pressure } else { 0.0 };
        q.array()[k] + ratio * (al * (left[k] - pressure) - ar * (right[k] - pressure))
    }))
}
pub(super) fn totals(
    states: &[Conserved],
    volumes: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Totals, Error> {
    let mut sum = [0.0; 3];
    let mut compensation = [0.0; 3];
    for (i, (q, &volume)) in states.iter().zip(volumes).enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        for k in 0..3 {
            let y = q.array()[k] * volume - compensation[k];
            let next = sum[k] + y;
            compensation[k] = (next - sum[k]) - y;
            sum[k] = next;
        }
    }
    Ok(Totals::from_array(sum))
}

fn validate_reservoirs(
    gas: &FrozenGas,
    boundaries: Boundaries,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    if let Boundaries::Ends { left, right } = boundaries {
        for boundary in [left, right] {
            poll(cancelled)?;
            if let Boundary::Reservoir(state) = boundary {
                if state.momentum != 0.0 {
                    return Err(Error::Invalid("stationary reservoir has nonzero momentum"));
                }
                gas.decode(state, None, cancelled)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "port_tests.rs"]
mod port_tests;
