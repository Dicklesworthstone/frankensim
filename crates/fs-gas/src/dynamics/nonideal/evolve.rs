//! Fixed-domain real-gas evolution sharing the native finite-volume machinery.

use super::{RealGas, RealGasFlowError, RealGasWallOperator};
use crate::body_force::ForceChange;
use crate::dynamics::evolve::{Stage, forward, totals};
use crate::dynamics::{
    Conserved, Error, FaceFlux, Mesh1d, Primitive, Reconstruction, StepControl, StepReport, Totals,
    allocated, poll, positive,
};
use crate::multicomponent::wall_exchange::WallControl;
use fs_qty::{Length, Time};

/// One end of a stationary constant-area real-gas interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RealGasBoundary {
    /// Zero-gradient exterior trace; not a nonreflecting subsonic outlet.
    Transmissive,
    /// Mirrored velocity through the same EOS-checked HLLC/HLL kernel. Wall
    /// mass and energy flux are exactly zero; pressure traction is numerical.
    Reflecting,
    /// Prescribed exterior trace in this owner's SAME frozen composition and
    /// energy reference. This is a Riemann ghost state, not a reservoir law.
    Exterior(Conserved),
}

/// The same periodic interface is evaluated once and shared at both ends.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RealGasBoundaries {
    Periodic,
    Ends {
        left: RealGasBoundary,
        right: RealGasBoundary,
    },
}

struct Workspace {
    faces: Vec<FaceFlux>,
    stage: Vec<Conserved>,
    stage_primitive: Vec<Primitive>,
    candidate: Vec<Conserved>,
    candidate_primitive: Vec<Primitive>,
}

/// Conservative density, momentum and total-energy evolution on an original
/// constant-area mesh. Nonuniform cell widths use the retained physical volumes.
///
/// The initial spatial path is piecewise constant. Time stepping uses SSP-RK2
/// and requires `StepControl::reconstruction == Reconstruction::FirstOrder`.
/// Every Euler stage and final average resolves through the original real-gas
/// source. Bounded retries halve dt; refusal/cancellation publishes no state.
/// Numerical scratch is allocated once; native material resolution retains its
/// own allocations. Optional local wall sources use the same gas calorics;
/// no gas phase flash, spatial viscous transport or moving end wall is added.
pub struct RealGasEuler1d {
    gas: RealGas,
    states: Vec<Conserved>,
    primitive: Vec<Primitive>,
    mesh: Mesh1d,
    boundaries: RealGasBoundaries,
    boundary_primitive: [Option<Primitive>; 2],
    time: f64,
    work: Workspace,
}

impl RealGasEuler1d {
    pub fn new(
        gas: RealGas,
        cell_width: Length,
        cross_section_m2: f64,
        states: Vec<Conserved>,
        boundaries: RealGasBoundaries,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        let mesh = Mesh1d::uniform(states.len(), cell_width, cross_section_m2)?;
        Self::new_on_mesh(gas, mesh, states, boundaries, cancelled)
    }

    pub fn new_on_mesh(
        gas: RealGas,
        mesh: Mesh1d,
        states: Vec<Conserved>,
        boundaries: RealGasBoundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        poll(&mut cancelled)?;
        let n = states.len();
        if n != mesh.cell_count() {
            return Err(RealGasFlowError::Invalid(
                "state count differs from original mesh",
            ));
        }
        let area = mesh
            .uniform_cross_section()
            .ok_or(RealGasFlowError::Invalid(
                "homogeneous real-gas field currently requires constant area",
            ))?;
        for (i, &volume) in mesh.volumes().iter().enumerate() {
            poll(&mut cancelled)?;
            // The original uniform constructor retains its declared width;
            // subtracting distant rounded coordinates can lose that precision.
            let width = mesh
                .uniform_cell_width()
                .map_or_else(|| mesh.width(i), Length::value);
            let expected = area * width;
            if !positive(expected)
                || (volume - expected).abs() > 64.0 * f64::EPSILON * volume.max(expected)
            {
                return Err(RealGasFlowError::Invalid(
                    "constant-area cell volume differs from original face geometry",
                ));
            }
        }
        let boundary_primitive = admit_boundaries(&gas, boundaries, &mut cancelled)?;
        let mut primitive = allocated(n, Primitive::default())?;
        for (p, &q) in primitive.iter_mut().zip(&states) {
            *p = gas.decode(q, &mut cancelled)?;
        }
        let work = Workspace {
            faces: allocated(n + 1, FaceFlux::default())?,
            stage: allocated(n, Conserved::default())?,
            stage_primitive: allocated(n, Primitive::default())?,
            candidate: allocated(n, Conserved::default())?,
            candidate_primitive: allocated(n, Primitive::default())?,
        };
        poll(&mut cancelled)?;
        Ok(Self {
            gas,
            states,
            primitive,
            mesh,
            boundaries,
            boundary_primitive,
            time: 0.0,
            work,
        })
    }

    pub const fn gas(&self) -> &RealGas {
        &self.gas
    }
    pub fn states(&self) -> &[Conserved] {
        &self.states
    }
    /// Cached original-source primitives for the complete accepted field.
    pub fn primitives(&self) -> &[Primitive] {
        &self.primitive
    }
    pub const fn mesh(&self) -> &Mesh1d {
        &self.mesh
    }
    pub const fn boundaries(&self) -> RealGasBoundaries {
        self.boundaries
    }
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }

    pub fn totals(&self, mut cancelled: impl FnMut() -> bool) -> Result<Totals, RealGasFlowError> {
        let total = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
        if total.array().iter().any(|value| !value.is_finite()) {
            return Err(RealGasFlowError::Invalid(
                "nonfinite extensive real-gas state",
            ));
        }
        Ok(total)
    }

    /// Admit new prescribed end traces atomically between accepted steps.
    pub fn set_boundaries(
        &mut self,
        boundaries: RealGasBoundaries,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), RealGasFlowError> {
        let primitive = admit_boundaries(&self.gas, boundaries, &mut cancelled)?;
        poll(&mut cancelled)?;
        self.boundaries = boundaries;
        self.boundary_primitive = primitive;
        Ok(())
    }

    /// Advance by at most the supplied duration. CFL and retries use the same
    /// controls/report as the original gas owner. This entry explicitly requires
    /// piecewise-constant spatial reconstruction; no implicit MUSCL downgrade.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        self.step_inner(maximum_duration, controls, None, cancelled)
    }

    /// Advance the original wave field, then implicit wall heat/drag over the
    /// ACTUAL accepted duration. Finite-wall storage, gas states, cached EOS
    /// primitives, reports and clock publish atomically after the final poll.
    ///
    /// This is first-order transport/source splitting. `walls.last_report()`
    /// retains the original finite-storage and reservoir/drive accounts;
    /// `StepReport` includes the actual gas wall impulse, heat and work.
    pub fn step_with_walls(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        walls: &mut RealGasWallOperator<'_>,
        wall_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        self.step_inner(
            maximum_duration,
            controls,
            Some((walls, wall_controls)),
            cancelled,
        )
    }

    #[allow(clippy::too_many_lines)]
    fn step_inner(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut walls: Option<(&mut RealGasWallOperator<'_>, WallControl)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        poll(&mut cancelled)?;
        if let Some((walls, wall_controls)) = &walls {
            walls.validate_binding(self)?;
            wall_controls.validate()?;
        }
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
            || controls.reconstruction != Reconstruction::FirstOrder
        {
            return Err(RealGasFlowError::Invalid(
                "positive duration, CFL <= 0.45, at most 32 halvings and FirstOrder reconstruction required",
            ));
        }
        let area = self.mesh.face_areas()[0];
        let mut rate = 0.0_f64;
        for (p, &volume) in self.primitive.iter().zip(self.mesh.volumes()) {
            poll(&mut cancelled)?;
            let local = (p.velocity.abs() + p.sound_speed) * (area / volume);
            if !positive(local) {
                return Err(RealGasFlowError::Invalid(
                    "unrepresentable real-gas signal rate",
                ));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * controls.cfl / rate);
        let before = self.totals(&mut cancelled)?;
        let mut last = RealGasFlowError::Invalid("no admissible real-gas step");
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || !(self.time + dt).is_finite() || self.time + dt <= self.time {
                return Err(RealGasFlowError::Invalid(
                    "time step cannot advance represented clock",
                ));
            }
            let result = attempt(
                &self.gas,
                &self.states,
                &self.primitive,
                &self.mesh,
                self.boundaries,
                self.boundary_primitive,
                dt,
                controls.cfl,
                &mut self.work,
                &mut cancelled,
            )
            .and_then(|(a, b)| {
                let thermal = match &mut walls {
                    Some((walls, wall_controls)) => Some(walls.prepare(
                        &mut self.work.candidate,
                        &mut self.work.candidate_primitive,
                        Time::new(dt),
                        *wall_controls,
                        &mut cancelled,
                    )?),
                    None => None,
                };
                Ok((a, b, thermal))
            });
            match result {
                Ok((a, b, mut thermal)) => {
                    let after = totals(&self.work.candidate, self.mesh.volumes(), &mut cancelled)?;
                    let factor = 0.5 * dt * area;
                    let left = Totals::from_array(core::array::from_fn(|k| {
                        factor * (a.left.array()[k] + b.left.array()[k])
                    }));
                    let right = Totals::from_array(core::array::from_fn(|k| {
                        -factor * (a.right.array()[k] + b.right.array()[k])
                    }));
                    let boundary = Totals::from_array(core::array::from_fn(|k| {
                        factor
                            * ((a.left.array()[k] + b.left.array()[k])
                                - (a.right.array()[k] + b.right.array()[k]))
                    }));
                    let (wall_impulse, wall_heat, wall_work) =
                        thermal.as_ref().map_or((0.0, 0.0, 0.0), |report| {
                            (
                                -report.wall_impulse_ns,
                                report.heat_to_gas_j,
                                report.mechanical_work_to_gas_j,
                            )
                        });
                    let source = [0.0, wall_impulse, wall_heat + wall_work];
                    let residual = Totals::from_array(core::array::from_fn(|k| {
                        after.array()[k] - before.array()[k] - boundary.array()[k] - source[k]
                    }));
                    if let Some(report) = &mut thermal {
                        report.gas_wall_storage_residual_j = after.energy
                            - before.energy
                            - boundary.energy
                            + report.wall_storage_change_j
                            - report.external_heat_to_gas_j
                            - report.mechanical_work_to_gas_j;
                        if !report.gas_wall_storage_residual_j.is_finite() {
                            return Err(RealGasFlowError::Invalid(
                                "nonfinite real-gas wall storage balance",
                            ));
                        }
                    }
                    if after
                        .array()
                        .iter()
                        .chain(&left.array())
                        .chain(&right.array())
                        .chain(&boundary.array())
                        .chain(&residual.array())
                        .any(|x| !x.is_finite())
                    {
                        return Err(RealGasFlowError::Invalid(
                            "nonfinite extensive real-gas balance",
                        ));
                    }
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.primitive, &mut self.work.candidate_primitive);
                    self.time += dt;
                    if let (Some((walls, _)), Some(report)) = (&mut walls, thermal) {
                        walls.commit(report);
                    }
                    return Ok(StepReport {
                        elapsed: Time::new(dt),
                        maximum_cfl: a.signal_rate.max(b.signal_rate) * dt,
                        halvings,
                        hll_fallback_faces: a.fallbacks + b.fallbacks,
                        limited_cells: 0,
                        first_order_retry: false,
                        before,
                        after,
                        boundary_change: boundary,
                        left_boundary_change: left,
                        right_boundary_change: right,
                        area_momentum_impulse: 0.0,
                        gravity_momentum_impulse: 0.0,
                        gravity_work: 0.0,
                        wall_momentum_impulse: wall_impulse,
                        wall_heat,
                        wall_work,
                        source_change: ForceChange::default(),
                        conservation_residual: residual,
                    });
                }
                Err(RealGasFlowError::Cancelled) => return Err(RealGasFlowError::Cancelled),
                Err(error) => {
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(RealGasFlowError::StepRetries {
            attempts: controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}

fn admit_boundaries(
    gas: &RealGas,
    boundaries: RealGasBoundaries,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Option<Primitive>; 2], RealGasFlowError> {
    poll(cancelled)?;
    let mut primitive = [None; 2];
    if let RealGasBoundaries::Ends { left, right } = boundaries {
        for (i, boundary) in [left, right].into_iter().enumerate() {
            if let RealGasBoundary::Exterior(state) = boundary {
                primitive[i] = Some(gas.decode(state, cancelled)?);
            }
        }
    }
    Ok(primitive)
}

#[allow(clippy::too_many_arguments)]
fn attempt(
    gas: &RealGas,
    states: &[Conserved],
    primitive: &[Primitive],
    mesh: &Mesh1d,
    boundaries: RealGasBoundaries,
    boundary_primitive: [Option<Primitive>; 2],
    dt: f64,
    cfl: f64,
    work: &mut Workspace,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage), RealGasFlowError> {
    let a = operator(
        gas,
        states,
        primitive,
        mesh,
        boundaries,
        boundary_primitive,
        &mut work.faces,
        cancelled,
    )?;
    cfl_gate(a.signal_rate * dt, cfl)?;
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        let next = forward(q, i, dt, mesh, primitive, &work.faces, None);
        work.stage_primitive[i] = if next == q {
            primitive[i]
        } else {
            gas.decode(next, cancelled)?
        };
        work.stage[i] = next;
    }
    let b = operator(
        gas,
        &work.stage,
        &work.stage_primitive,
        mesh,
        boundaries,
        boundary_primitive,
        &mut work.faces,
        cancelled,
    )?;
    cfl_gate(b.signal_rate * dt, cfl)?;
    for (i, &q) in work.stage.iter().enumerate() {
        poll(cancelled)?;
        let second = forward(q, i, dt, mesh, &work.stage_primitive, &work.faces, None);
        let second_primitive = if second == q {
            work.stage_primitive[i]
        } else {
            gas.decode(second, cancelled)?
        };
        let candidate = Conserved::from_array(core::array::from_fn(|k| {
            0.5 * states[i].array()[k] + 0.5 * second.array()[k]
        }));
        work.candidate_primitive[i] = if candidate == states[i] {
            primitive[i]
        } else if candidate == second {
            second_primitive
        } else {
            gas.decode(candidate, cancelled)?
        };
        work.candidate[i] = candidate;
    }
    Ok((a, b))
}

fn cfl_gate(observed: f64, maximum: f64) -> Result<(), RealGasFlowError> {
    if !observed.is_finite() || observed > maximum * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded { observed, maximum }.into());
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn operator(
    gas: &RealGas,
    states: &[Conserved],
    primitive: &[Primitive],
    mesh: &Mesh1d,
    boundaries: RealGasBoundaries,
    boundary_primitive: [Option<Primitive>; 2],
    faces: &mut [FaceFlux],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, RealGasFlowError> {
    let n = states.len();
    let mut fallbacks = 0;
    for face in 0..=n {
        poll(cancelled)?;
        if face == n && boundaries == RealGasBoundaries::Periodic {
            faces[n] = faces[0];
            continue;
        }
        let (left, right, wall) = if face == 0 || face == n {
            match boundaries {
                RealGasBoundaries::Periodic => (
                    (states[n - 1], primitive[n - 1]),
                    (states[0], primitive[0]),
                    false,
                ),
                RealGasBoundaries::Ends { left, right } => {
                    let index = if face == 0 { 0 } else { n - 1 };
                    let interior = (states[index], primitive[index]);
                    let boundary = if face == 0 { left } else { right };
                    let mut exterior = interior;
                    match boundary {
                        RealGasBoundary::Transmissive => (),
                        RealGasBoundary::Reflecting => {
                            exterior.0.momentum = -exterior.0.momentum;
                            exterior.1.velocity = -exterior.1.velocity;
                        }
                        RealGasBoundary::Exterior(q) => {
                            exterior = (
                                q,
                                boundary_primitive[usize::from(face == n)]
                                    .expect("prescribed exterior admitted at construction/update"),
                            );
                        }
                    }
                    if face == 0 {
                        (exterior, interior, boundary == RealGasBoundary::Reflecting)
                    } else {
                        (interior, exterior, boundary == RealGasBoundary::Reflecting)
                    }
                }
            }
        } else {
            (
                (states[face - 1], primitive[face - 1]),
                (states[face], primitive[face]),
                false,
            )
        };
        let mut flux = gas.interface(left.0, left.1, right.0, right.1, cancelled)?;
        if wall {
            // Exact impermeability/work at a stationary mirror boundary. The
            // source-checked numerical Riemann solve supplies pressure traction.
            flux.mass = 0.0;
            flux.energy = 0.0;
        }
        fallbacks += usize::from(flux.used_hll_fallback);
        faces[face] = flux;
    }
    let area = mesh.face_areas()[0];
    let mut signal_rate = 0.0_f64;
    for (i, &volume) in mesh.volumes().iter().enumerate() {
        poll(cancelled)?;
        let rate = area * (0.5 * faces[i].signal_speed + 0.5 * faces[i + 1].signal_speed) / volume;
        if !positive(rate) {
            return Err(RealGasFlowError::Invalid(
                "unrepresentable real-gas face CFL",
            ));
        }
        signal_rate = signal_rate.max(rate);
    }
    Ok(Stage {
        left: faces[0],
        right: faces[n],
        signal_rate,
        area_force: 0.0,
        gravity_force: 0.0,
        gravity_power: 0.0,
        fallbacks,
        limited: 0,
    })
}
