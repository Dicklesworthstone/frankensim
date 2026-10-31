//! A sounding bore on the shared nonlinear gas owner.
//!
//! [`NonlinearBore`] keeps one [`Euler1d`] field of the original NASA-9 source
//! and adds what a played wind bore needs that the inviscid duct lacks:
//!
//! - Radius steps inside cells. Cell volumes are exact; the step's wall force
//!   enters through the same quasi-one-dimensional `p dA` source as a taper.
//!   Like the plane-wave transfer-matrix chain, this resolves no evanescent
//!   step modes.
//! - Zwikker–Kosten wall memory per cell: the `fs_phs::ZkDiffusive` drag and
//!   diffusive arms of the local radius, applied between accepted steps by
//!   their exact propagator (a symmetric negative semidefinite generator, so
//!   the update is contractive). Wall friction turns kinetic into internal
//!   energy; wall heat exchange changes internal energy at fixed density.
//! - A radiating end: a rigid circular piston in an infinite baffle. Its
//!   radiation impedance has the exact impulse response
//!   `z(t) = (rho c / S) [delta(t) - (4/(pi T)) sqrt(1 - (t/T)^2)]`, `0 <= t <= T = 2a/c`
//!   (Stepanishen, JASA 49 (1971) 841), the transform of
//!   `rho c / S [1 - J1(2ka)/ka + i H1(2ka)/ka]`. It is product-integrated
//!   against the piecewise-linear end flow on the fixed step clock and solved
//!   implicitly with the interior's own end wave at every RK stage.
//!
//! The driven end is a caller-owned port. Walls are rigid, the room is at
//! rest, and the bell radiates as a baffled piston: no flare-mode or
//! unflanged-lip radiation is claimed.

use super::{NonlinearDuctError, area, positive, profile_volume, radius, reserved, unsupported};
use crate::{Duct, DuctError, Segment};
use fs_gas::dynamics::{
    Boundaries, Boundary, Conserved, Error as GasError, Euler1d, FaceFlux, FrozenGas, Mesh1d,
    PortStage, Primitive, Reconstruction, StepControl, StepReport,
};
use fs_math::det;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::Time;

/// Zwikker–Kosten wall-memory admission for every cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallMemory {
    /// Viscous and thermal diffusive poles per cell (`fs_phs` cap applies).
    pub poles: usize,
    /// Admitted acoustic band [Hz]; the binding checks both wall operators on it.
    pub band_hz: [f64; 2],
    /// Largest sampled relative error of either wall operator on the band.
    pub maximum_relative_error: f64,
}

/// Transport of the bore gas at the wall, from the same source as the field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallGas {
    /// Density [kg/m³].
    pub density_kg_m3: f64,
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Heat-capacity ratio.
    pub gamma: f64,
    /// Prandtl number.
    pub prandtl: f64,
}

/// Discretization of a sounding bore.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoreSpec {
    /// Largest axial cell [m]. The bore is cut into equal cells; a face that
    /// would fall within a quarter cell of a radius step is moved off it.
    pub maximum_cell_m: f64,
    /// Compact extra volume on the first cell (a lumped cup) [m³], or zero.
    pub inlet_volume_m3: f64,
    /// Fixed accepted step [s]. The radiating end's memory runs on this clock;
    /// a field that cannot take it refuses rather than shortening it.
    pub step_s: f64,
    /// Admitted stage CFL, in (0, 0.45].
    pub cfl: f64,
    /// Wall memory, or `None` for rigid inviscid adiabatic walls.
    pub wall: Option<WallMemory>,
}

/// One accepted bore step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoreStep {
    /// The shared gas owner's accepted-step report, before the wall update.
    pub report: StepReport,
    /// Mean outward mass flow through the radiating end over the step [kg/s].
    pub outlet_mass_flow_kg_s: f64,
    /// That mass flow divided by the room density: the piston's volume
    /// velocity, the source strength the exterior pressure is referenced to [m³/s].
    pub outlet_volume_flow_m3_s: f64,
    /// Radiating-end pressure solved at the step's start [Pa].
    pub outlet_pressure_pa: f64,
}

/// Wall-memory receipt: what the bound arms admit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallReceipt {
    /// Distinct wall radii bound (cells of one radius share one propagator).
    pub distinct_radii: usize,
    /// Largest viscous or thermal arm count in any cell.
    pub maximum_arms: usize,
    /// Worst sampled complex error of the excess wall operators.
    pub maximum_complex_error: f64,
    /// Worst sampled error of their dissipative parts.
    pub maximum_loss_error: f64,
}

/// A bore whose gas evolves nonlinearly between a caller port and a radiating end.
pub struct NonlinearBore {
    solver: Euler1d,
    gas: FrozenGas,
    room: Primitive,
    radiation: PistonRadiation,
    walls: Option<Walls>,
    step_s: f64,
    cfl: f64,
    centers_m: Vec<f64>,
}

impl NonlinearBore {
    /// Mesh `duct`, fill it with room gas at rest and bind its wall memory.
    ///
    /// `room` is the stationary exterior (and initial bore) state: its density
    /// and sound speed set the piston's radiation, and its pressure is the
    /// reference of the thermal wall arms. `wall_gas` is required with
    /// `spec.wall` and is ignored without it.
    pub fn new(
        duct: &Duct,
        gas: FrozenGas,
        room_temperature_k: f64,
        room_pressure_pa: f64,
        wall_gas: Option<WallGas>,
        spec: BoreSpec,
    ) -> Result<Self, NonlinearDuctError> {
        if !positive(spec.step_s) || !(spec.cfl > 0.0 && spec.cfl <= 0.45) {
            return Err(unsupported(
                "bore step must be positive and CFL in (0, 0.45]",
            ));
        }
        let geometry = bore_mesh(duct, spec.maximum_cell_m, spec.inlet_volume_m3)?;
        let room_density = room_pressure_pa / (gas.gas_constant() * room_temperature_k);
        let rest = gas.conserved(room_density, 0.0, room_pressure_pa)?;
        let room = gas.primitive(rest, || false)?;
        let n = geometry.wall_radii.len();
        let mut cells = reserved(n)?;
        cells.resize(n, rest);
        let mouth_radius = duct
            .segments
            .last()
            .map(Segment::outlet_radius)
            .ok_or(DuctError::EmptyDuct)?;
        let radiation = PistonRadiation::new(mouth_radius, room, spec.step_s)?;
        let walls = match (spec.wall, wall_gas) {
            (None, _) => None,
            (Some(_), None) => return Err(unsupported("wall memory needs the wall gas transport")),
            (Some(memory), Some(transport)) => Some(Walls::new(
                &geometry.wall_radii,
                memory,
                transport,
                room.pressure,
                spec.step_s,
            )?),
        };
        let solver = Euler1d::new_on_mesh(
            gas.clone(),
            geometry.mesh,
            cells,
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            || false,
        )?;
        Ok(Self {
            solver,
            gas,
            room,
            radiation,
            walls,
            step_s: spec.step_s,
            cfl: spec.cfl,
            centers_m: geometry.centers_m,
        })
    }

    /// Accept one fixed step. `inlet` supplies the driven end's flux at each
    /// RK stage (`None`: a closed end); the radiating end is this bore's own.
    /// The wall memory is applied after the gas step, over the same duration.
    /// A refused gas step (including one whose CFL bound is shorter than the
    /// fixed step, which the radiating end refuses) leaves the field and both
    /// memories unchanged. A wall refusal after an accepted gas step leaves
    /// the bore mid-step and must end the run.
    pub fn step(
        &mut self,
        mut inlet: impl FnMut(&PortStage) -> Result<Option<FaceFlux>, GasError>,
    ) -> Result<BoreStep, NonlinearDuctError> {
        let radiation = &mut self.radiation;
        let gas = &self.gas;
        let controls = StepControl {
            cfl: self.cfl,
            reconstruction: Reconstruction::Muscl,
            maximum_halvings: 0,
        };
        let report = self.solver.step_with_ports(
            Time::new(self.step_s),
            controls,
            |stage| Ok([inlet(stage)?, Some(radiation.flux(gas, stage)?)]),
            || false,
        )?;
        let outlet_pressure_pa = radiation.commit()?;
        if let Some(walls) = self.walls.as_mut() {
            walls.apply(&mut self.solver)?;
        }
        let outlet_mass_flow_kg_s = -report.right_boundary_change.mass / self.step_s;
        Ok(BoreStep {
            report,
            outlet_mass_flow_kg_s,
            outlet_volume_flow_m3_s: outlet_mass_flow_kg_s / self.room.density,
            outlet_pressure_pa,
        })
    }

    /// The shared gas field, its mesh, source and time.
    #[must_use]
    pub const fn solver(&self) -> &Euler1d {
        &self.solver
    }

    /// The source-owned gas law (the same one the field evolves).
    #[must_use]
    pub const fn gas(&self) -> &FrozenGas {
        &self.gas
    }

    /// The stationary room state outside the radiating end.
    #[must_use]
    pub const fn room(&self) -> Primitive {
        self.room
    }

    /// Cell midpoints [m].
    #[must_use]
    pub fn centers_m(&self) -> &[f64] {
        &self.centers_m
    }

    /// Driven-end face area [m²].
    #[must_use]
    pub fn inlet_area_m2(&self) -> f64 {
        self.solver.mesh().face_areas()[0]
    }

    /// Number of past end flows the radiation memory convolves.
    #[must_use]
    pub fn radiation_taps(&self) -> usize {
        self.radiation.weights.len()
    }

    /// Product-integration weights of the regular piston kernel, newest first
    /// (index 0 multiplies the current end flow). They sum to one: the
    /// radiation impedance vanishes at zero frequency.
    #[must_use]
    pub fn radiation_weights(&self) -> &[f64] {
        &self.radiation.weights
    }

    /// What the wall memory admitted, or `None` for rigid inviscid walls.
    #[must_use]
    pub fn wall_receipt(&self) -> Option<WallReceipt> {
        self.walls.as_ref().map(|walls| walls.receipt)
    }
}

struct BoreMesh {
    mesh: Mesh1d,
    wall_radii: Vec<f64>,
    centers_m: Vec<f64>,
}

// Equal cells over the whole bore, faces moved a quarter cell off radius steps,
// exact piecewise volumes, face areas from the original geometry.
fn bore_mesh(
    duct: &Duct,
    maximum_cell_m: f64,
    inlet_volume_m3: f64,
) -> Result<BoreMesh, NonlinearDuctError> {
    if duct.segments.is_empty() {
        return Err(DuctError::EmptyDuct.into());
    }
    if !positive(maximum_cell_m) || !(inlet_volume_m3 >= 0.0 && inlet_volume_m3.is_finite()) {
        return Err(unsupported("bore cell size and inlet volume"));
    }
    let mut starts = reserved(duct.segments.len() + 1)?;
    starts.push(0.0);
    for segment in &duct.segments {
        segment.validate()?;
        if matches!(segment, Segment::ToneHole { .. }) {
            return Err(unsupported(
                "tone holes require a conservative gas junction",
            ));
        }
        let end = starts.last().copied().unwrap_or(0.0) + segment.axial_length();
        if !positive(end) {
            return Err(unsupported("bore length is not representable"));
        }
        starts.push(end);
    }
    let length = starts[duct.segments.len()];
    let cells = (length / maximum_cell_m).ceil();
    if !(1.0..=super::MAX_NONLINEAR_DUCT_CELLS as f64).contains(&cells) {
        return Err(unsupported("cell count exceeds MAX_NONLINEAR_DUCT_CELLS"));
    }
    let n = cells as usize;
    let h = length / cells;
    let mut faces = reserved(n + 1)?;
    for i in 0..=n {
        faces.push(if i == n {
            length
        } else {
            length * (i as f64 / cells)
        });
    }
    for k in 1..duct.segments.len() {
        let before = duct.segments[k - 1].outlet_radius();
        let after = duct.segments[k].inlet_radius();
        if (before - after).abs() <= 64.0 * f64::EPSILON * before.max(after) {
            continue;
        }
        let step = starts[k];
        let i = (step / h).round() as usize;
        if i > 0 && i < n && (faces[i] - step).abs() < 0.25 * h {
            faces[i] = if faces[i] <= step {
                step - 0.25 * h
            } else {
                step + 0.25 * h
            };
        }
    }
    if faces.windows(2).any(|w| !(w[1] - w[0] >= 0.4 * h)) {
        return Err(unsupported("two radius steps lie within one cell"));
    }
    let locate = |x: f64| -> usize {
        let k = starts.partition_point(|&s| s <= x);
        k.saturating_sub(1).min(duct.segments.len() - 1)
    };
    let mut areas = reserved(n + 1)?;
    for &x in &faces {
        let k = locate(x);
        areas.push(area(radius(
            &duct.segments[k],
            (x - starts[k]).clamp(0.0, duct.segments[k].axial_length()),
        )?)?);
    }
    let mut volumes = reserved(n)?;
    let mut wall_radii = reserved(n)?;
    let mut centers_m = reserved(n)?;
    for i in 0..n {
        let (left, right) = (faces[i], faces[i + 1]);
        let mut volume = 0.0;
        for (k, segment) in duct.segments.iter().enumerate() {
            let a = left.max(starts[k]) - starts[k];
            let b = right.min(starts[k + 1]) - starts[k];
            if b <= a {
                continue;
            }
            let b = b.min(segment.axial_length());
            volume += match segment {
                Segment::Cylinder { radius, .. } => area(*radius)? * (b - a),
                Segment::Cone { .. } => {
                    let (ra, rb) = (radius(segment, a)?, radius(segment, b)?);
                    core::f64::consts::PI * ((b - a) / 3.0) * (ra * ra + ra * rb + rb * rb)
                }
                Segment::Profile { profile, .. } => profile_volume(profile, a, b, &mut || false)?,
                Segment::ToneHole { .. } => unreachable!("topology admitted above"),
            };
        }
        if !positive(volume) {
            return Err(unsupported("cell volume is not positive and finite"));
        }
        // The wall memory sees the cylinder of the same volume and width.
        wall_radii.push(det::sqrt(volume / (core::f64::consts::PI * (right - left))));
        centers_m.push(left.midpoint(right));
        volumes.push(if i == 0 {
            volume + inlet_volume_m3
        } else {
            volume
        });
    }
    Ok(BoreMesh {
        mesh: Mesh1d::new(faces, areas, volumes)?,
        wall_radii,
        centers_m,
    })
}

// The regular part of the piston kernel against hat functions on the step
// grid: closed-form antiderivatives of sqrt(1-s^2) and s sqrt(1-s^2).
fn piston_weights(step_s: f64, transit_s: f64) -> Result<Vec<f64>, NonlinearDuctError> {
    let h = step_s / transit_s;
    if !positive(h) {
        return Err(unsupported("piston transit time and step"));
    }
    let taps = (1.0 / h).ceil() + 1.0;
    if !(taps <= 1_048_576.0) {
        return Err(unsupported(
            "piston radiation memory exceeds its tap budget",
        ));
    }
    let f0 = |s: f64| 0.5 * (s * det::sqrt(1.0 - s * s) + det::asin(s));
    let f1 = |s: f64| -(1.0 - s * s) * det::sqrt(1.0 - s * s) / 3.0;
    let mut weights = reserved(taps as usize)?;
    for j in 0..taps as usize {
        let centre = j as f64 * h;
        let mut weight = 0.0;
        // Rising half on [s_{j-1}, s_j], falling half on [s_j, s_{j+1}], both cut at s = 1.
        if j > 0 {
            let (a, b) = (centre - h, centre.min(1.0));
            if a < b {
                let i0 = f0(b) - f0(a);
                let i1 = f1(b) - f1(a);
                weight += (i1 - a * i0) / h;
            }
        }
        let (a, b) = (centre, (centre + h).min(1.0));
        if a < b {
            let i0 = f0(b) - f0(a);
            let i1 = f1(b) - f1(a);
            weight += ((centre + h) * i0 - i1) / h;
        }
        weights.push(4.0 / core::f64::consts::PI * weight);
    }
    Ok(weights)
}

struct PistonRadiation {
    weights: Vec<f64>,
    history: Vec<f64>,
    newest: usize,
    pending: Option<(f64, f64)>,
    rho_c: f64,
    impedance: f64,
    area: f64,
    room: Primitive,
    step_s: f64,
    guess: f64,
}

impl PistonRadiation {
    fn new(radius_m: f64, room: Primitive, step_s: f64) -> Result<Self, NonlinearDuctError> {
        let area = area(radius_m)?;
        let transit = 2.0 * radius_m / room.sound_speed;
        let weights = piston_weights(step_s, transit)?;
        let mut history = reserved(weights.len())?;
        history.resize(weights.len(), 0.0);
        let rho_c = room.density * room.sound_speed;
        Ok(Self {
            weights,
            history,
            newest: 0,
            pending: None,
            rho_c,
            impedance: rho_c / area,
            area,
            room,
            step_s,
            guess: room.pressure,
        })
    }

    // Past end flows U_{t-j}, j >= 1, convolved with the kernel weights. The
    // corrector's newest past flow is the predictor's uncommitted one.
    fn memory(&self, stage: usize) -> f64 {
        let m = self.history.len();
        let mut sum = 0.0;
        let mut first = 1;
        if stage == 1 {
            sum += self.weights[1] * self.pending.map_or(0.0, |(flow, _)| flow);
            first = 2;
        }
        for j in first..self.weights.len() {
            let back = j - first;
            sum += self.weights[j] * self.history[(self.newest + m - back) % m];
        }
        sum
    }

    fn flux(&mut self, gas: &FrozenGas, stage: &PortStage) -> Result<FaceFlux, GasError> {
        if stage.dt != self.step_s
            || stage.stage > 1
            || (stage.stage == 1 && self.pending.is_none())
        {
            return Err(GasError::Invalid(
                "radiating end runs on its declared fixed step",
            ));
        }
        let interior = stage.right.1;
        let held = self.impedance * self.memory(stage.stage);
        let gain = self.rho_c * (1.0 - self.weights[0]);
        // G(p) = p - p_room - gain u(p) + held is strictly increasing in p.
        let mut pressure = self.guess;
        let mut wave = gas.end_wave(interior, 1.0, pressure, || false)?;
        for _ in 0..64 {
            let residual = pressure - self.room.pressure - gain * wave.primitive.velocity + held;
            let slope = 1.0 - gain * wave.velocity_per_pressure;
            let next = pressure - residual / slope;
            if !positive(next) {
                return Err(GasError::Invalid(
                    "radiating end pressure left the gas domain",
                ));
            }
            let converged = (next - pressure).abs() <= 1e-13 * pressure;
            pressure = next;
            wave = gas.end_wave(interior, 1.0, pressure, || false)?;
            if converged {
                break;
            }
        }
        let velocity = wave.primitive.velocity;
        let residual = pressure - self.room.pressure - gain * velocity + held;
        if !(residual.abs() <= 1e-9 * self.room.pressure) {
            return Err(GasError::Invalid(
                "radiating end characteristic did not converge",
            ));
        }
        // Outflow carries the interior's entropy. Inflow carries the room's: the
        // room brought to the end pressure along this source's wave curve
        // (isentropic in expansion, its Hugoniot in compression).
        let face = if velocity >= 0.0 {
            wave.primitive
        } else {
            let room_side = gas.end_wave(self.room, 1.0, pressure, || false)?;
            Primitive {
                velocity,
                ..room_side.primitive
            }
        };
        if stage.stage == 0 {
            self.pending = Some((self.area * velocity, pressure));
        }
        self.guess = pressure;
        gas.face_flux(face, wave.signal_speed)
    }

    fn commit(&mut self) -> Result<f64, NonlinearDuctError> {
        let (flow, pressure) = self
            .pending
            .take()
            .ok_or_else(|| unsupported("radiating end has no accepted predictor flow"))?;
        let m = self.history.len();
        self.newest = (self.newest + 1) % m;
        self.history[self.newest] = flow;
        self.guess = pressure;
        Ok(pressure)
    }
}

// One cell's exactly propagated wall arms, held in the eigen-coordinates of
// their energy-scaled generator, so a step costs a few passes over the arms.
struct Propagator {
    size: usize,
    // G = Q diag(lambda) Q^T: each eigen-coordinate decays by exp(lambda dt);
    // `head` is Q's first row, the eigen-components of the mean coordinate.
    decay: Vec<f64>,
    head: Vec<f64>,
    #[cfg(test)]
    basis: Vec<f64>,
}

impl Propagator {
    // G = [[-(d + sum g), sqrt(g p)^T], [sqrt(g p), -diag(p)]] is symmetric and
    // negative semidefinite (x^T G x = -d x0^2 - sum (sqrt(g) x0 - sqrt(p) xk)^2),
    // so exp(G dt) = Q diag(exp(lambda dt)) Q^T is a contraction.
    fn new(drag: f64, arms: &[(f64, f64)], dt: f64) -> Result<Self, NonlinearDuctError> {
        let size = arms.len() + 1;
        let mut g = vec![0.0; size * size];
        g[0] = -(drag + arms.iter().map(|&(gain, _)| gain).sum::<f64>());
        for (k, &(gain, rate)) in arms.iter().enumerate() {
            let coupling = det::sqrt(gain * rate);
            g[k + 1] = coupling;
            g[(k + 1) * size] = coupling;
            g[(k + 1) * size + k + 1] = -rate;
        }
        let scale = g.iter().map(|x| x.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || !(dt.is_finite() && dt > 0.0) {
            return Err(unsupported("wall arm generator is not finite"));
        }
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&g, size);
        for (k, &lambda) in values.iter().enumerate() {
            let mut residual = 0.0_f64;
            let mut norm = 0.0;
            for i in 0..size {
                let product: f64 = (0..size)
                    .map(|j| g[i * size + j] * vectors[j * size + k])
                    .sum();
                residual = residual.max((product - lambda * vectors[i * size + k]).abs());
                norm += vectors[i * size + k] * vectors[i * size + k];
            }
            if !lambda.is_finite()
                || lambda > 1e-12 * scale
                || residual > 1e-12 * scale
                || (norm - 1.0).abs() > 1e-12
            {
                return Err(unsupported(
                    "wall arm generator is not resolved as negative semidefinite",
                ));
            }
        }
        Ok(Self {
            size,
            // The exact generator is semidefinite; rounding may leave its
            // (thermal) zero mode a hair above zero, whose decay is one.
            decay: values.iter().map(|&l| det::exp(l.min(0.0) * dt)).collect(),
            head: vectors[..size].to_vec(),
            #[cfg(test)]
            basis: vectors,
        })
    }

    // Set the mean coordinate to `mean` (the arms keep their state), advance
    // one step in place, and return the new mean with |x|^2 before and after.
    fn advance(&self, state: &mut [f64], mean: f64) -> (f64, f64, f64) {
        let current: f64 = self.head.iter().zip(&*state).map(|(h, y)| h * y).sum();
        let shift = mean - current;
        let (mut before, mut after, mut next) = (0.0, 0.0, 0.0);
        for ((y, &h), &d) in state.iter_mut().zip(&self.head).zip(&self.decay) {
            let shifted = *y + shift * h;
            before += shifted * shifted;
            *y = shifted * d;
            after += *y * *y;
            next += h * *y;
        }
        (next, before, after)
    }

    // The dense exp(G dt), for the tests.
    #[cfg(test)]
    fn matrix(&self) -> Vec<f64> {
        let n = self.size;
        let mut out = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                out[r * n + c] = (0..n)
                    .map(|k| self.basis[r * n + k] * self.decay[k] * self.basis[c * n + k])
                    .sum();
            }
        }
        out
    }
}

struct Walls {
    // Per cell: indices of its viscous and thermal propagators, and arm states.
    viscous_of: Vec<usize>,
    thermal_of: Vec<usize>,
    propagators: Vec<Propagator>,
    viscous_state: Vec<f64>,
    thermal_state: Vec<f64>,
    spare_viscous: Vec<f64>,
    spare_thermal: Vec<f64>,
    state_stride: usize,
    reference_pressure: f64,
    receipt: WallReceipt,
}

impl Walls {
    fn new(
        radii: &[f64],
        memory: WallMemory,
        gas: WallGas,
        reference_pressure: f64,
        dt: f64,
    ) -> Result<Self, NonlinearDuctError> {
        let basis = ZkDiffusiveBasis::new(memory.poles)
            .map_err(|_| unsupported("wall memory pole count"))?;
        let band = [
            core::f64::consts::TAU * memory.band_hz[0],
            core::f64::consts::TAU * memory.band_hz[1],
        ];
        let mut bound: Vec<(u64, usize, usize)> = Vec::new();
        let mut propagators = Vec::new();
        let mut viscous_of = reserved(radii.len())?;
        let mut thermal_of = reserved(radii.len())?;
        let mut receipt = WallReceipt {
            distinct_radii: 0,
            maximum_arms: 0,
            maximum_complex_error: 0.0,
            maximum_loss_error: 0.0,
        };
        for &radius in radii {
            let key = radius.to_bits();
            let (v, t) = if let Some(&(_, v, t)) = bound.iter().find(|(k, _, _)| *k == key) {
                (v, t)
            } else {
                let zk = basis
                    .bind(
                        radius,
                        gas.density_kg_m3,
                        gas.dynamic_viscosity_pa_s,
                        gas.gamma,
                        gas.prandtl,
                        band,
                        memory.maximum_relative_error,
                    )
                    .map_err(|_| {
                        unsupported(
                            "wall memory does not resolve this radius on the band: add poles, narrow the band or relax the error",
                        )
                    })?;
                let arms = |terms: &[fs_phs::ZkDiffusiveTerm]| {
                    terms
                        .iter()
                        .map(|t| (t.gain_per_s, t.rate_per_s))
                        .collect::<Vec<_>>()
                };
                let viscous = arms(zk.series_terms());
                let thermal = arms(zk.shunt_terms());
                receipt.maximum_arms = receipt.maximum_arms.max(viscous.len()).max(thermal.len());
                receipt.maximum_complex_error =
                    receipt.maximum_complex_error.max(zk.max_complex_error());
                receipt.maximum_loss_error =
                    receipt.maximum_loss_error.max(zk.max_real_loss_error());
                propagators.push(Propagator::new(zk.series_drag_per_s(), &viscous, dt)?);
                propagators.push(Propagator::new(0.0, &thermal, dt)?);
                let entry = (key, propagators.len() - 2, propagators.len() - 1);
                bound.push(entry);
                (entry.1, entry.2)
            };
            viscous_of.push(v);
            thermal_of.push(t);
        }
        receipt.distinct_radii = bound.len();
        let stride = receipt.maximum_arms + 1;
        Ok(Self {
            viscous_of,
            thermal_of,
            propagators,
            viscous_state: vec![0.0; radii.len() * stride],
            thermal_state: vec![0.0; radii.len() * stride],
            spare_viscous: vec![0.0; radii.len() * stride],
            spare_thermal: vec![0.0; radii.len() * stride],
            state_stride: stride,
            reference_pressure,
            receipt,
        })
    }

    fn apply(&mut self, solver: &mut Euler1d) -> Result<(), NonlinearDuctError> {
        let stride = self.state_stride;
        let reference = self.reference_pressure;
        let Self {
            propagators,
            viscous_of,
            thermal_of,
            viscous_state,
            thermal_state,
            spare_viscous,
            spare_thermal,
            ..
        } = self;
        // Every cell writes its whole next arm state into the spare buffers,
        // which are published only with the accepted cells.
        solver.apply_cell_sources(
            |i, q, p| {
                let at = i * stride;
                let viscous = &propagators[viscous_of[i]];
                let n = viscous.size;
                let arms = &mut spare_viscous[at..at + n];
                arms.copy_from_slice(&viscous_state[at..at + n]);
                let (velocity, before, after) = viscous.advance(arms, p.velocity);
                // The arms hold the boundary layer's extra kinetic energy,
                // rho |z|^2 / 2; only what the propagator removes from the
                // mean flow plus arms is dissipated, and it heats the gas.
                let dissipated = 0.5 * q.density * (before - after);
                let kinetic_change =
                    0.5 * q.density * (velocity * velocity - p.velocity * p.velocity);
                let gamma = p.sound_speed * p.sound_speed * p.density / p.pressure;
                let heated = p.pressure + (gamma - 1.0) * dissipated;
                let moved = q.energy + kinetic_change + dissipated;
                let thermal = &propagators[thermal_of[i]];
                let n = thermal.size;
                let arms = &mut spare_thermal[at..at + n];
                arms.copy_from_slice(&thermal_state[at..at + n]);
                let (exchanged, _, _) = thermal.advance(arms, heated - reference);
                let exchanged = exchanged + reference;
                // Wall heat exchange changes internal energy at fixed density.
                Ok(Conserved {
                    density: q.density,
                    momentum: q.density * velocity,
                    energy: moved + (exchanged - heated) / (gamma - 1.0),
                })
            },
            || false,
        )?;
        core::mem::swap(viscous_state, spare_viscous);
        core::mem::swap(thermal_state, spare_thermal);
        Ok(())
    }
}

#[cfg(test)]
#[path = "bore_tests.rs"]
mod tests;
