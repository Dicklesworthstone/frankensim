use super::flux::Face;
use super::{
    Boundaries, Boundary, Cell, Conserved, Error, MAX_COMPONENTS, MAX_SPECIES, Mesh1d, MixtureGas,
    Primitive, Reconstruction, StepControl, allocated, cell, poll, positive, sum,
};
use fs_qty::{Length, Time};

/// Extensive conserved integrals in the canonical source-species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Totals {
    /// Species masses [kg]. Only the first `species_count` entries are active.
    pub species_masses: [f64; MAX_SPECIES],
    /// Active canonical species count; remaining mass entries are exactly zero.
    pub species_count: usize,
    /// Axial momentum [kg m/s].
    pub momentum: f64,
    /// Species-referenced total energy [J].
    pub energy: f64,
}
impl Totals {
    /// Total mass, derived from all species masses.
    pub fn mass(&self) -> f64 {
        sum(self.species_masses[..self.species_count].iter().copied())
    }
    /// Physical total energy with the owner's original formation references.
    pub fn absolute_energy(&self, gas: &MixtureGas) -> Result<f64, Error> {
        if self.species_count != gas.species_count() {
            return Err(Error::Invalid("total/source species axis"));
        }
        let energy = self.energy
            + sum(self.species_masses[..self.species_count]
                .iter()
                .zip(gas.energy_references())
                .map(|(a, b)| a * b));
        if !energy.is_finite() {
            return Err(Error::Invalid("non-finite physical extensive energy"));
        }
        Ok(energy)
    }
    pub(super) fn from_values(v: &[f64], n: usize) -> Self {
        let mut species_masses = [0.0; MAX_SPECIES];
        species_masses[..n].copy_from_slice(&v[..n]);
        Self {
            species_masses,
            species_count: n,
            momentum: v[n],
            energy: v[n + 1],
        }
    }
    pub(super) fn get(&self, k: usize) -> f64 {
        if k < self.species_count {
            self.species_masses[k]
        } else if k == self.species_count {
            self.momentum
        } else {
            self.energy
        }
    }
}

/// One complete accepted SSP RK2 step, including species boundary balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    /// Accepted duration, no greater than the caller's maximum.
    pub elapsed: Time,
    /// Actual maximum stage face-speed CFL.
    pub maximum_cfl: f64,
    /// Number of rejected, halved attempts.
    pub halvings: usize,
    /// Accepted-stage interfaces using the admissible HLL fallback.
    pub hll_fallback_faces: usize,
    /// Accepted-stage cells with reduced reconstruction slopes.
    pub limited_cells: usize,
    /// Whether an admissibility retry switched the step to first-order space.
    pub first_order_retry: bool,
    /// Original extensive conserved quantities.
    pub before: Totals,
    /// Accepted extensive conserved quantities.
    pub after: Totals,
    /// Area/time-integrated left inflow minus right outflow.
    pub boundary_change: Totals,
    /// Accepted exchange into the branch through its left end.
    pub left_boundary_change: Totals,
    /// Accepted exchange into the branch through its right end.
    pub right_boundary_change: Totals,
    /// Stationary area-wall pressure impulse [kg m/s].
    pub area_momentum_impulse: f64,
    /// `after - before - boundary_change - area momentum impulse`.
    pub conservation_residual: Totals,
}

pub(super) struct Workspace {
    pub(super) primitive: Vec<Primitive>,
    pub(super) slopes: Vec<f64>,
    pub(super) faces: Vec<f64>,
    pub(super) face_data: Vec<Face>,
    pub(super) stage: Vec<f64>,
    pub(super) candidate: Vec<f64>,
}

/// Conservative variable-composition finite-volume flow on a stationary mesh.
///
/// A contiguous cell-major species/momentum/energy buffer and preallocated RK
/// and face buffers keep successful numerical steps allocation-free. All
/// species share each face's wave speeds and one conservative energy flux.
/// No species clipping, trace seeding or composition normalization is applied.
/// Periodic, reflecting and transmissive boundaries are supported. The shared
/// frozen-gas `Boundary::Reservoir` variant refuses because it cannot specify
/// the required species masses and their separate formation-energy references.
/// Prescribed species-bearing reservoirs use [`Self::step_with_reservoirs`].
pub struct MixtureEuler1d {
    pub(super) gas: MixtureGas,
    pub(super) states: Vec<f64>,
    pub(super) mesh: Mesh1d,
    pub(super) boundaries: Boundaries,
    pub(super) time: f64,
    pub(super) work: Workspace,
}
impl MixtureEuler1d {
    /// Construct a uniform fixed-area flow from complete conservative cells.
    pub fn new(
        gas: MixtureGas,
        width: Length,
        area: f64,
        states: Vec<Conserved>,
        boundaries: Boundaries,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mesh = Mesh1d::uniform(states.len(), width, area)?;
        Self::new_on_mesh(gas, mesh, states, boundaries, cancelled)
    }
    /// Construct on a stationary nonuniform, possibly variable-area mesh.
    pub fn new_on_mesh(
        gas: MixtureGas,
        mesh: Mesh1d,
        states: Vec<Conserved>,
        boundaries: Boundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        admit_boundaries(boundaries)?;
        let count = mesh.cell_count();
        let n = gas.species_count();
        let m = n + 2;
        if states.len() != count {
            return Err(Error::Invalid("state count differs from mesh"));
        }
        let mut flat = allocated(count.checked_mul(m).ok_or(Error::Allocation)?, 0.0)?;
        for (i, q) in states.iter().enumerate() {
            poll(&mut cancelled)?;
            if q.partial_densities.len() != n {
                return Err(Error::Invalid("cell species axis"));
            }
            flat[i * m..i * m + n].copy_from_slice(&q.partial_densities);
            flat[i * m + n] = q.momentum;
            flat[i * m + n + 1] = q.energy;
        }
        Self::from_flat(gas, mesh, flat, boundaries, 0.0, &mut cancelled)
    }
    fn from_flat(
        gas: MixtureGas,
        mesh: Mesh1d,
        states: Vec<f64>,
        boundaries: Boundaries,
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        admit_boundaries(boundaries)?;
        let count = mesh.cell_count();
        let m = gas.species_count() + 2;
        if states.len() != count * m {
            return Err(Error::Invalid("flattened cell dimension"));
        }
        if boundaries == Boundaries::Periodic && mesh.face_areas()[0] != mesh.face_areas()[count] {
            return Err(Error::Invalid("periodic end areas differ"));
        }
        let mut primitive = allocated(count, Primitive::default())?;
        for (i, q) in states.chunks_exact(m).enumerate() {
            primitive[i] = gas.decode(q, None, cancelled)?;
        }
        // Public totals are infallible only after every extensive quantity has
        // been admitted. Finite densities alone do not prevent volume overflow.
        totals(&states, mesh.volumes(), m - 2, cancelled)?;
        let work = Workspace {
            primitive,
            // Primitive reconstruction carries rho separately from every Yi.
            slopes: allocated(count * (m + 1), 0.0)?,
            faces: allocated((count + 1) * m, 0.0)?,
            face_data: allocated(count + 1, Face::default())?,
            stage: allocated(count * m, 0.0)?,
            candidate: allocated(count * m, 0.0)?,
        };
        poll(cancelled)?;
        Ok(Self {
            gas,
            states,
            mesh,
            boundaries,
            time,
            work,
        })
    }
    /// Source bank, formation-energy gauge and caloric domain.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Number of finite-volume cells.
    pub fn cell_count(&self) -> usize {
        self.mesh.cell_count()
    }
    /// Borrow a cell in canonical species order; no reconstruction or allocation.
    pub fn cell(&self, index: usize) -> Option<Cell<'_>> {
        let m = self.gas.species_count() + 2;
        if index >= self.cell_count() {
            return None;
        }
        Some(cell(&self.states[index * m..(index + 1) * m], m - 2))
    }
    /// Original stationary physical geometry.
    pub const fn mesh(&self) -> &Mesh1d {
        &self.mesh
    }
    /// Declared endpoint conditions of this accepted flow.
    pub fn boundaries(&self) -> Boundaries {
        self.boundaries
    }
    /// Current accepted time.
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Compensated extensive species, momentum and energy integrals.
    pub fn totals(&self) -> Totals {
        totals(
            &self.states,
            self.mesh.volumes(),
            self.gas.species_count(),
            &mut || false,
        )
        .expect("infallible cancellation callback")
    }
    /// Create a same-owner scratch flow once for an atomic coupled integrator.
    pub(super) fn duplicate(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        let mut states = allocated(self.states.len(), 0.0)?;
        states.copy_from_slice(&self.states);
        Self::from_flat(
            self.gas.clone(),
            self.mesh.clone(),
            states,
            self.boundaries,
            self.time,
            cancelled,
        )
    }
    /// Synchronize already admitted physical state without allocating. The
    /// caller holds two instances constructed from the same immutable owner.
    pub(super) fn copy_physical_from(&mut self, other: &Self) -> Result<(), Error> {
        if self.states.len() != other.states.len()
            || self.gas != other.gas
            || self.mesh != other.mesh
            || self.boundaries != other.boundaries
        {
            return Err(Error::Invalid("coupled scratch flow has a different owner"));
        }
        self.states.copy_from_slice(&other.states);
        self.time = other.time;
        Ok(())
    }
    /// Validate a complete in-place scratch field before coupled publication.
    pub(super) fn validate_states(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for q in self.states.chunks_exact(self.gas.species_count() + 2) {
            self.gas.decode(q, None, cancelled)?;
        }
        totals(
            &self.states,
            self.mesh.volumes(),
            self.gas.species_count(),
            cancelled,
        )?;
        Ok(())
    }
    /// Advance one atomic conservative step. Actual face CFL is checked at both
    /// stages; inadmissibility halves dt and may reduce spatial order.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        self.advance(maximum_duration, controls, [None; 2], cancelled)
    }

    /// Advance against prescribed stationary reservoirs at the left/right ends.
    ///
    /// Each `Some` overrides that end's declared physical boundary for this
    /// call; `None` retains it. Both `None` execute the original step exactly.
    /// Periodic intervals cannot have reservoirs. The borrowed cells must use
    /// this gas owner's canonical species order and per-species energy gauges,
    /// exactly as for [`MixtureGas::reservoir_flux`]. Momentum must be zero;
    /// the complete source state is admitted before starting a trial.
    ///
    /// Pressure and temperature are stagnation conditions for inflow; the
    /// reservoir pressure supplies outflow backpressure. Inflow transports the
    /// reservoir composition and stagnation enthalpy, while outflow transports
    /// the evolving interior composition and total enthalpy. The original
    /// characteristic closure handles sonic transitions and original calorics.
    /// Its flux is recomputed from the actual end cell at BOTH RK stages and
    /// every halved retry, with the same face CFL and extensive balance report.
    ///
    /// These reservoirs remain prescribed throughout the accepted step and
    /// are not depleted. Change their source states between calls for an
    /// explicit boundary drive; use the existing finite network owner for
    /// reservoir feedback. No reservoir state or flux is retained, no numerical
    /// workspace is allocated per step, and refusal/cancellation publishes no
    /// cell or time change.
    pub fn step_with_reservoirs(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        if reservoirs.iter().all(Option::is_none) {
            return self.step(maximum_duration, controls, cancelled);
        }
        poll(&mut cancelled)?;
        if self.boundaries == Boundaries::Periodic {
            return Err(Error::Invalid(
                "prescribed mixture reservoirs require physical ends",
            ));
        }
        let mut admitted = [None; 2];
        for (target, reservoir) in admitted.iter_mut().zip(reservoirs) {
            if let Some(state) = reservoir {
                poll(&mut cancelled)?;
                if state.momentum != 0.0 {
                    return Err(Error::Invalid(
                        "stationary mixture reservoir has nonzero momentum",
                    ));
                }
                *target = Some((state, self.gas.primitive(state, &mut cancelled)?));
            }
        }
        self.advance(maximum_duration, controls, admitted, cancelled)
    }

    fn advance(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        reservoirs: ReservoirEnds<'_>,
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
        let n = self.gas.species_count();
        let m = n + 2;
        let mut signal_rate = 0.0_f64;
        for (i, q) in self.states.chunks_exact(m).enumerate() {
            let p = self
                .gas
                .decode(q, Some(self.work.primitive[i].temperature), &mut cancelled)?;
            self.work.primitive[i] = p;
            let area = 0.5 * self.mesh.face_areas()[i] + 0.5 * self.mesh.face_areas()[i + 1];
            signal_rate = signal_rate
                .max((p.velocity.abs() + p.sound_speed) * (area / self.mesh.volumes()[i]));
        }
        let mut dt = maximum_duration
            .value()
            .min(0.9 * controls.cfl / signal_rate);
        let before = totals(&self.states, self.mesh.volumes(), n, &mut cancelled)?;
        let mut last = Error::Invalid("no admissible time increment");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(Error::Invalid("time step cannot advance represented clock"));
            }
            let order = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            match attempt(
                &self.gas,
                &self.states,
                &self.mesh,
                self.boundaries,
                dt,
                controls.cfl,
                order,
                reservoirs,
                &mut self.work,
                &mut cancelled,
            ) {
                Ok((a, b)) => {
                    let after =
                        totals(&self.work.candidate, self.mesh.volumes(), n, &mut cancelled)?;
                    let mut boundary = [0.0; MAX_COMPONENTS];
                    let mut left = [0.0; MAX_COMPONENTS];
                    let mut right = [0.0; MAX_COMPONENTS];
                    let mut residual = [0.0; MAX_COMPONENTS];
                    let area_momentum_impulse = 0.5 * dt * (a.area_force + b.area_force);
                    for k in 0..m {
                        left[k] = 0.5 * dt * self.mesh.face_areas()[0] * (a.left[k] + b.left[k]);
                        right[k] = -0.5
                            * dt
                            * self.mesh.face_areas()[self.cell_count()]
                            * (a.right[k] + b.right[k]);
                        boundary[k] = 0.5
                            * dt
                            * (self.mesh.face_areas()[0] * (a.left[k] + b.left[k])
                                - self.mesh.face_areas()[self.cell_count()]
                                    * (a.right[k] + b.right[k]));
                        residual[k] = after.get(k)
                            - before.get(k)
                            - boundary[k]
                            - if k == n { area_momentum_impulse } else { 0.0 };
                    }
                    if !boundary[..m]
                        .iter()
                        .chain(&left[..m])
                        .chain(&right[..m])
                        .chain(&residual[..m])
                        .all(|x| x.is_finite())
                        || !area_momentum_impulse.is_finite()
                    {
                        return Err(Error::Invalid("non-finite mixture extensive balance"));
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
                        boundary_change: Totals::from_values(&boundary[..m], n),
                        left_boundary_change: Totals::from_values(&left[..m], n),
                        right_boundary_change: Totals::from_values(&right[..m], n),
                        area_momentum_impulse,
                        conservation_residual: Totals::from_values(&residual[..m], n),
                    });
                }
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(e) => {
                    if !matches!(e, Error::CflExceeded { .. }) {
                        first_order_retry = true;
                    }
                    last = e;
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

fn admit_boundaries(boundaries: Boundaries) -> Result<(), Error> {
    if matches!(
        boundaries,
        Boundaries::Ends {
            left: Boundary::Reservoir(_),
            ..
        } | Boundaries::Ends {
            right: Boundary::Reservoir(_),
            ..
        }
    ) {
        return Err(Error::Invalid(
            "mixture reservoir requires species-bearing state and energy references",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy)]
pub(super) struct EndpointFlux<'a> {
    pub(super) values: &'a [f64],
    pub(super) data: Face,
}

type ReservoirEnds<'a> = [Option<(Cell<'a>, Primitive)>; 2];

// The state and primitive refer to the SAME RK stage. Only the prescribed
// outside primitive is reusable; the interior trace and flux are always fresh.
fn reservoir_fluxes<'a>(
    gas: &MixtureGas,
    states: &[f64],
    primitive: &[Primitive],
    reservoirs: ReservoirEnds<'_>,
    storage: &'a mut [[f64; MAX_COMPONENTS]; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Option<EndpointFlux<'a>>; 2], Error> {
    let n = gas.species_count();
    let m = n + 2;
    let mut result = [None; 2];
    for (end, (reservoir, values)) in reservoirs.into_iter().zip(storage.iter_mut()).enumerate() {
        if let Some((outside, outside_primitive)) = reservoir {
            let i = if end == 0 { 0 } else { states.len() / m - 1 };
            let inside = cell(&states[i * m..(i + 1) * m], n);
            let inside_primitive =
                gas.decode_cell(inside, Some(primitive[i].temperature), cancelled)?;
            let data = gas.plenum_flux(
                inside,
                inside_primitive,
                outside,
                outside_primitive,
                if end == 0 { -1.0 } else { 1.0 },
                &mut values[..m],
                cancelled,
            )?;
            result[end] = Some(EndpointFlux {
                values: &values[..m],
                data,
            });
        }
    }
    Ok(result)
}

#[derive(Clone, Copy)]
pub(super) struct Stage {
    pub(super) left: [f64; MAX_COMPONENTS],
    pub(super) right: [f64; MAX_COMPONENTS],
    pub(super) signal_rate: f64,
    pub(super) area_force: f64,
    pub(super) fallbacks: usize,
    pub(super) limited: usize,
}
#[allow(clippy::too_many_arguments)]
fn attempt(
    gas: &MixtureGas,
    states: &[f64],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    reservoirs: ReservoirEnds<'_>,
    work: &mut Workspace,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage), Error> {
    let n = gas.species_count();
    let m = n + 2;
    let mut endpoint_storage = [[0.0; MAX_COMPONENTS]; 2];
    let ends = reservoir_fluxes(
        gas,
        states,
        &work.primitive,
        reservoirs,
        &mut endpoint_storage,
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
        &mut work.face_data,
        ends,
        cancelled,
    )?;
    check_cfl(a.signal_rate * dt, cfl)?;
    for (i, q) in states.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        forward(
            q,
            i,
            dt,
            mesh,
            &work.primitive,
            &work.faces,
            &mut work.stage[i * m..(i + 1) * m],
            n,
        );
        if !gas.admissible(&work.stage[i * m..(i + 1) * m]) {
            return Err(Error::Invalid(
                "first mixture RK stage outside source domain",
            ));
        }
    }
    let ends = reservoir_fluxes(
        gas,
        &work.stage,
        &work.primitive,
        reservoirs,
        &mut endpoint_storage,
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
        &mut work.face_data,
        ends,
        cancelled,
    )?;
    check_cfl(b.signal_rate * dt, cfl)?;
    let mut second = [0.0; MAX_COMPONENTS];
    for (i, q) in work.stage.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        forward(
            q,
            i,
            dt,
            mesh,
            &work.primitive,
            &work.faces,
            &mut second[..m],
            n,
        );
        if !gas.admissible(&second[..m]) {
            return Err(Error::Invalid(
                "second mixture RK stage outside source domain",
            ));
        }
        let candidate = &mut work.candidate[i * m..(i + 1) * m];
        for k in 0..m {
            candidate[k] = 0.5 * states[i * m + k] + 0.5 * second[k];
        }
        gas.decode(candidate, Some(work.primitive[i].temperature), cancelled)?;
    }
    Ok((a, b))
}
fn check_cfl(observed: f64, maximum: f64) -> Result<(), Error> {
    if !observed.is_finite() || observed > maximum * (1.0 + 64.0 * f64::EPSILON) {
        Err(Error::CflExceeded { observed, maximum })
    } else {
        Ok(())
    }
}

// Reconstruct Yi, rho, velocity and pressure, with one common composition
// limiter. The stored cell averages remain the original conserved species.
// A uniform inert or linear elemental fraction must stay uniform through a
// density wave; separately limiting partial densities does not preserve that.
#[allow(clippy::too_many_arguments)]
pub(super) fn operator(
    gas: &MixtureGas,
    states: &[f64],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    reconstruction: Reconstruction,
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    faces: &mut [f64],
    face_data: &mut [Face],
    endpoint_fluxes: [Option<EndpointFlux<'_>>; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, Error> {
    let n = gas.species_count();
    let m = n + 2;
    let primitive_width = n + 3;
    let count = mesh.cell_count();
    for (i, q) in states.chunks_exact(m).enumerate() {
        primitive[i] = gas.decode(q, Some(primitive[i].temperature), cancelled)?;
    }
    slopes.fill(0.0);
    let value = |i: usize, k: usize| {
        if k < n {
            states[i * m + k] / primitive[i].density
        } else if k == n {
            primitive[i].density
        } else if k == n + 1 {
            primitive[i].velocity
        } else {
            primitive[i].pressure
        }
    };
    let center =
        |i: usize| mesh.positions()[i] + 0.5 * (mesh.positions()[i + 1] - mesh.positions()[i]);
    let mut limited = 0;
    let mut trial = [0.0; MAX_COMPONENTS + 1];
    let mut densities = [0.0; MAX_SPECIES];
    let mut out = [0.0; MAX_COMPONENTS];
    for i in 0..count {
        poll(cancelled)?;
        if reconstruction == Reconstruction::FirstOrder {
            continue;
        }
        let (prev, next) = if boundaries == Boundaries::Periodic {
            ((i + count - 1) % count, (i + 1) % count)
        } else if i == 0 || i + 1 == count {
            continue;
        } else {
            (i - 1, i + 1)
        };
        let length = mesh.positions()[count] - mesh.positions()[0];
        let x = center(i);
        let xl = center(prev) - if i == 0 { length } else { 0.0 };
        let xr = center(next) + if i + 1 == count { length } else { 0.0 };
        let width = mesh.positions()[i + 1] - mesh.positions()[i];
        let row = &mut slopes[i * primitive_width..(i + 1) * primitive_width];
        super::reconstruction::slopes_with_geometry(
            &states[prev * m..prev * m + n],
            &states[i * m..i * m + n],
            &states[next * m..next * m + n],
            [width / (x - xl), width / (xr - xl), width / (xr - x)],
            &mut row[..n],
        )?;
        for k in n..primitive_width {
            let (a, b, c) = (value(prev, k), value(i, k), value(next, k));
            row[k] = minmod(
                2.0 * (b - a) / (x - xl),
                (c - a) / (xr - xl),
                2.0 * (c - b) / (xr - x),
            ) * width;
        }
        let mut scale = 1.0;
        let mut valid = false;
        for reduction in 0..=12 {
            valid = true;
            for sign in [-0.5, 0.5] {
                for k in 0..primitive_width {
                    trial[k] = value(i, k) + sign * scale * row[k];
                }
                for k in 0..n {
                    densities[k] = trial[n] * trial[k];
                }
                if !positive(trial[n])
                    || gas
                        .state_at(&densities[..n], trial[n + 1], trial[n + 2], &mut out[..m])
                        .is_err()
                {
                    valid = false;
                    break;
                }
            }
            if valid {
                break;
            }
            if reduction == 0 {
                limited += 1;
            }
            scale *= 0.5;
        }
        for value in row {
            *value = if valid { scale * *value } else { 0.0 };
        }
    }
    let reconstruct = |i: usize, sign: f64, target: &mut [f64]| -> Result<Primitive, Error> {
        let row = &slopes[i * primitive_width..(i + 1) * primitive_width];
        if row.iter().all(|x| *x == 0.0) {
            target.copy_from_slice(&states[i * m..(i + 1) * m]);
            return Ok(primitive[i]);
        }
        let mut densities = [0.0; MAX_SPECIES];
        let density = value(i, n) + sign * row[n];
        if !positive(density) {
            return Err(Error::Invalid("reconstructed mixture density"));
        }
        for k in 0..n {
            densities[k] = density * (value(i, k) + sign * row[k]);
        }
        gas.state_at(
            &densities[..n],
            value(i, n + 1) + sign * row[n + 1],
            value(i, n + 2) + sign * row[n + 2],
            target,
        )
    };
    let mut left = [0.0; MAX_COMPONENTS];
    let mut right = [0.0; MAX_COMPONENTS];
    let mut fallbacks = 0;
    for face in 0..=count {
        poll(cancelled)?;
        if face == 0 || face == count {
            if let Some(flux) = endpoint_fluxes[usize::from(face == count)] {
                if boundaries == Boundaries::Periodic
                    || flux.values.len() != m
                    || !positive(flux.data.signal_speed)
                    || !flux.values.iter().all(|x| x.is_finite())
                {
                    return Err(Error::Invalid("invalid physical mixture endpoint flux"));
                }
                faces[face * m..(face + 1) * m].copy_from_slice(flux.values);
                face_data[face] = flux.data;
                fallbacks += usize::from(flux.data.used_hll_fallback);
                continue;
            }
        }
        if face == count && boundaries == Boundaries::Periodic {
            let (prefix, suffix) = faces.split_at_mut(count * m);
            suffix[..m].copy_from_slice(&prefix[..m]);
            face_data[count] = face_data[0];
            continue;
        }
        let (pl, pr, wall) = if face == 0 || face == count {
            match boundaries {
                Boundaries::Periodic => (
                    reconstruct(count - 1, 0.5, &mut left[..m])?,
                    reconstruct(0, -0.5, &mut right[..m])?,
                    false,
                ),
                Boundaries::Ends {
                    left: lb,
                    right: rb,
                } => {
                    let i = if face == 0 { 0 } else { count - 1 };
                    let sign = if face == 0 { -0.5 } else { 0.5 };
                    let inner = reconstruct(i, sign, &mut left[..m])?;
                    right[..m].copy_from_slice(&left[..m]);
                    let wall = if face == 0 { lb } else { rb } == Boundary::Reflecting;
                    (inner, inner, wall)
                }
            }
        } else {
            (
                reconstruct(face - 1, 0.5, &mut left[..m])?,
                reconstruct(face, -0.5, &mut right[..m])?,
                false,
            )
        };
        let flux = &mut faces[face * m..(face + 1) * m];
        let data = if wall {
            gas.wall_flux(
                cell(&left[..m], n),
                pl,
                if face == 0 { -1.0 } else { 1.0 },
                flux,
                cancelled,
            )?
        } else {
            gas.interface(&left[..m], pl, &right[..m], pr, flux, cancelled)?
        };
        fallbacks += usize::from(data.used_hll_fallback);
        face_data[face] = data;
    }
    let mut signal_rate = 0.0_f64;
    let mut force = 0.0;
    let mut correction = 0.0;
    for i in 0..count {
        poll(cancelled)?;
        let rate = (0.5 * mesh.face_areas()[i] * face_data[i].signal_speed
            + 0.5 * mesh.face_areas()[i + 1] * face_data[i + 1].signal_speed)
            / mesh.volumes()[i];
        if !positive(rate) {
            return Err(Error::Invalid("unrepresentable mixture face CFL"));
        }
        signal_rate = signal_rate.max(rate);
        let term = primitive[i].pressure * (mesh.face_areas()[i + 1] - mesh.face_areas()[i]);
        let y = term - correction;
        let next = force + y;
        correction = (next - force) - y;
        force = next;
    }
    let mut lf = [0.0; MAX_COMPONENTS];
    let mut rf = [0.0; MAX_COMPONENTS];
    lf[..m].copy_from_slice(&faces[..m]);
    rf[..m].copy_from_slice(&faces[count * m..(count + 1) * m]);
    Ok(Stage {
        left: lf,
        right: rf,
        signal_rate,
        area_force: force,
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
    q: &[f64],
    i: usize,
    dt: f64,
    mesh: &Mesh1d,
    primitive: &[Primitive],
    faces: &[f64],
    out: &mut [f64],
    n: usize,
) {
    let m = n + 2;
    let al = mesh.face_areas()[i];
    let ar = mesh.face_areas()[i + 1];
    let ratio = dt / mesh.volumes()[i];
    for k in 0..m {
        let pressure = if k == n { primitive[i].pressure } else { 0.0 };
        out[k] = q[k]
            + ratio
                * (al * (faces[i * m + k] - pressure) - ar * (faces[(i + 1) * m + k] - pressure));
    }
}
pub(super) fn totals(
    states: &[f64],
    volumes: &[f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Totals, Error> {
    let m = n + 2;
    let mut values = [0.0; MAX_COMPONENTS];
    let mut corrections = [0.0; MAX_COMPONENTS];
    for (i, q) in states.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        for k in 0..m {
            let y = q[k] * volumes[i] - corrections[k];
            let next = values[k] + y;
            corrections[k] = (next - values[k]) - y;
            values[k] = next;
        }
    }
    if !values[..m].iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("non-finite extensive mixture state"));
    }
    Ok(Totals::from_values(&values[..m], n))
}

#[cfg(test)]
#[path = "reservoir_boundary_tests.rs"]
mod reservoir_boundary_tests;
