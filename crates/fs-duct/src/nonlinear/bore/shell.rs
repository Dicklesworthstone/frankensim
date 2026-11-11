//! The wall of a flaring bore as a breathing shell, one radial oscillator per cell.
//!
//! Brass walls move a little under the bore pressure, most where the bell flares (Kausel, Zietlow and Moore,
//! JASA 128 (2010) 3161: the wall is a ring that breathes). Each cell here owns the radial displacement `s` of
//! its wall band and its velocity `v`:
//!
//! ```text
//! m v' + c v + k s = (p - p_room) A,     A = int 2 pi r dx  (the band's projected lateral area)
//! ```
//!
//! The band's mass is `sigma S` (surface density times the slant area `S = int 2 pi r sqrt(1 + r'^2) dx`) and its
//! stiffness is the classical thin cone membrane's: the hoop resultant of a cone under pressure is `p r / cos phi`,
//! so `s = p r^2 / (E t cos phi)` and `k = (E t / r^2) A cos phi`, with `cos phi = A / S`. The surface density
//! and `E t / r^2` come from a caller-owned [`ShellLaw`] (in the sounding bore, the landed membrane/Lame wall
//! reduction); a cell whose law ran its thick-wall branch carries the cylinder stiffness unchanged. `c` is the
//! viscous damping that gives the cell's breathing mode the stated quality factor at its own resonance,
//! `c = m omega0 / Q`. Bending and axial shell modes, a rim bead, the outside air's added mass and a
//! non-uniform thickness are not carried; the paper's own stiffness (`1 / cos^3 phi` against the membrane's
//! `1 / cos phi`, a factor `1 / cos^2 phi`) is not used.
//!
//! The gas sees the wall through its volume. When the band moves out by `dV = A ds` the cell's gas keeps its mass
//! and fills the larger volume, so the cell average is rescaled,
//! `rho' = rho V / (V + dV)`, `(rho u)' = (rho u) V / (V + dV)`, `E' = (E V - p_mid dV) / (V + dV)`, with `V` the
//! PHYSICAL volume `V_nominal + A s` before the step: the gas does the work `p_mid dV`. Each cell's physical mass
//! and momentum are unchanged by this update, to roundoff. The gas's adiabatic compliance comes from the Euler
//! solver's own equation of state.
//!
//! The coupling is second order in the step. The wall step holds ONE force over the step, and that force is
//! `F = A (p_mid - p_room)`, where `p_mid = p - g ds / 2` is the mean of the pressure before the step and the
//! pressure after it, `p - g ds` (the gas's own stiffness against the wall's motion, `g = rho c^2 A / V`, from the
//! cell's density and frozen sound speed). `ds` depends on `F` linearly, so `F` solves a scalar linear equation
//! (no iteration). Holding the start-of-step pressure instead feeds `g A ds^2 / 2` of energy into the coupled
//! system every step (a very lightly damped wall mode would grow); the midpoint force does not (the lumped-cell model
//! satwall/o1-model2.py: H/H0 stays within 1e-9 over 20000 steps against +9 % for the start-of-step force, and the
//! gap to the continuous solution falls about a thousandfold). The wall step is exact for its held force (the
//! exponential of the oscillator with its input, once, at binding). The update runs once per accepted step after the
//! gas step: Lie splitting of the wall with the rest of the gas dynamics.
//!
//! Energy: for the held force the oscillator obeys `W' - W + D = F ds` exactly (`W` its stored energy,
//! `D >= 0` the work its damper dissipates), and the gas loses `p_mid dV = F ds + p_room dV`, so
//! `gas work = stored energy + dissipated work + work against the room pressure` closes cell by cell and step by
//! step. [`ShellLedger`] reports the four sums. One limit: the gas step conserves the sum of `rho V_nominal` while the
//! stored density is the physical one, so the total physical mass is conserved only to a relative `A s / V`
//! (about `2 s / r`, 1e-7 here) of each flux-driven change.

use super::super::{NonlinearDuctError, positive, radius, unsupported};
use crate::{Duct, Segment};
use fs_gas::dynamics::{Conserved, Error as GasError, Euler1d, FrozenGas};
use fs_math::det;
use fs_phs::linear_exact::exponential;

/// Panels of the composite Simpson rule over one cell of a smooth profile segment.
const PROFILE_PANELS: usize = 64;
/// Relative agreement required between the supplied duct's cell volume and the bore's own.
const VOLUME_AGREEMENT: f64 = 1.0e-5;

/// What a cylinder strip of wall is, at one inner radius: the surface mass `sigma` [kg/m^2] and the stiffness per
/// inner area [Pa/m] of the strip, and which law produced them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellStrip {
    /// Surface density [kg/m^2] of the strip.
    pub surface_density_kg_m2: f64,
    /// Stiffness per inner area `p / s` of the cylinder strip [Pa/m].
    pub stiffness_per_area_pa_m: f64,
    /// True when the thin membrane expression ran (the classical cone factor `cos phi` then applies to the
    /// stiffness); false when the thick-wall law ran (the cylinder stiffness is carried unchanged).
    pub membrane: bool,
}

/// The source of a strip's mass and stiffness. The bore owner holds no material data.
pub trait ShellLaw {
    /// The strip at inner radius `radius_m`, or a named refusal outside the law's support.
    ///
    /// # Errors
    /// Whatever the law refuses (a missing property, an unsupported wall).
    fn strip(&self, radius_m: f64) -> Result<ShellStrip, NonlinearDuctError>;
}

/// What binding the shell admitted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellReceipt {
    /// Cells of the bore.
    pub cells: usize,
    /// Cells whose strip ran the membrane expression.
    pub membrane_cells: usize,
    /// Cells whose strip ran the thick-wall expression.
    pub thick_cells: usize,
    /// Largest half-angle `acos(A / S)` of any cell [rad].
    pub largest_half_angle_rad: f64,
    /// Largest half-angle among the thick-wall cells [rad] (no flare theory is claimed for them); 0 with none.
    pub largest_thick_half_angle_rad: f64,
    /// Lowest breathing resonance of any cell [Hz].
    pub lowest_resonance_hz: f64,
    /// Highest breathing resonance of any cell [Hz].
    pub highest_resonance_hz: f64,
    /// The quality factor of every cell's breathing mode.
    pub quality: f64,
}

/// The four sums of the shell's energy account, in joules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellLedger {
    /// Work the gas did on the walls, `sum p dV`.
    pub gas_work_j: f64,
    /// Work done against the room pressure, `sum p_room dV` (telescopes to `p_room` times the volume change).
    pub room_work_j: f64,
    /// Work dissipated by the walls' dampers, `sum D >= 0`.
    pub dissipated_j: f64,
    /// Energy now stored in the walls, `sum (m v^2 + k s^2) / 2`.
    pub stored_j: f64,
}

impl ShellLedger {
    /// `gas work - (stored + dissipated + room)`: zero when the account closes.
    #[must_use]
    pub fn residual_j(&self) -> f64 {
        self.gas_work_j - (self.stored_j + self.dissipated_j + self.room_work_j)
    }

    /// The residual over the largest of the four sums (0 when all are zero).
    #[must_use]
    pub fn closing_error(&self) -> f64 {
        let largest = [self.gas_work_j, self.room_work_j, self.dissipated_j, self.stored_j]
            .iter()
            .fold(0.0_f64, |m, x| m.max(x.abs()));
        if largest == 0.0 {
            0.0
        } else {
            self.residual_j().abs() / largest
        }
    }
}

/// One cell's band: volume, projected lateral area and slant area.
#[derive(Clone, Copy)]
struct Band {
    volume: f64,
    lateral: f64,
    slant: f64,
}

pub(super) struct ShellCell {
    /// The cell's gas volume [m^3], the bore's own.
    pub(super) volume: f64,
    /// Projected lateral area `A` [m^2]: the volume swept per unit radial displacement.
    pub(super) lateral: f64,
    /// Band mass [kg] and stiffness [N/m].
    pub(super) mass: f64,
    pub(super) stiffness: f64,
    /// The exact held-force step: `s' = a s + b v + e F`, `v' = c s + d v + f F` as `[a, b, c, d, e, f]`.
    pub(super) step: [f64; 6],
}

impl ShellCell {
    /// One held-force wall step from `(s, v)` under `force` [N]: `(s', v')`.
    pub(super) fn advance(&self, state: [f64; 2], force: f64) -> [f64; 2] {
        let [a, b, c, d, e, f] = self.step;
        [
            a * state[0] + b * state[1] + e * force,
            c * state[0] + d * state[1] + f * force,
        ]
    }

    /// The force held over the step, from the cell pressure `pressure` at its start: `A (p_mid - p_room)` with
    /// `p_mid = pressure - g ds / 2`, `g = gas_stiffness A / physical_volume`, `ds` the step's displacement, linear
    /// in the force (module doc).
    pub(super) fn midpoint_force(
        &self,
        state: [f64; 2],
        pressure: f64,
        room_pressure: f64,
        gas_stiffness: f64,
        physical_volume: f64,
    ) -> f64 {
        let [a, b, _, _, e, _] = self.step;
        let coupling = gas_stiffness * self.lateral * self.lateral / physical_volume;
        let free = (a - 1.0) * state[0] + b * state[1];
        (self.lateral * (pressure - room_pressure) - 0.5 * coupling * free) / (1.0 + 0.5 * coupling * e)
    }
}

pub(super) struct Shell {
    pub(super) cells: Vec<ShellCell>,
    pub(super) state: Vec<[f64; 2]>,
    next: Vec<[f64; 2]>,
    gas_work: f64,
    room_work: f64,
    dissipated: f64,
    pub(super) receipt: ShellReceipt,
}

impl Shell {
    /// Bind the shell on the bore whose faces are `positions` and whose cell volumes are `volumes`; `duct` must be
    /// that bore's own duct.
    pub(super) fn new(
        duct: &Duct,
        positions: &[f64],
        volumes: &[f64],
        law: &dyn ShellLaw,
        quality: f64,
        dt: f64,
        inlet_volume_m3: f64,
    ) -> Result<(Self, ShellReceipt), NonlinearDuctError> {
        if !positive(quality) {
            return Err(unsupported("shell quality factor must be positive and finite"));
        }
        if !positive(dt) {
            return Err(unsupported("shell step must be positive and finite"));
        }
        let n = volumes.len();
        if n == 0 || positions.len() != n + 1 {
            return Err(unsupported("shell faces and volumes disagree"));
        }
        let mut starts = Vec::with_capacity(duct.segments.len() + 1);
        starts.push(0.0);
        for segment in &duct.segments {
            let end = starts.last().copied().unwrap_or(0.0) + segment.axial_length();
            if !positive(end) {
                return Err(unsupported("bore length is not representable"));
            }
            starts.push(end);
        }
        let length = starts[duct.segments.len()];
        if (positions[n] - length).abs() > 1.0e-9 * length || positions[0] != 0.0 {
            return Err(unsupported("the duct is not the bore's own duct (length)"));
        }
        let mut cells = Vec::with_capacity(n);
        let mut receipt = ShellReceipt {
            cells: n,
            membrane_cells: 0,
            thick_cells: 0,
            largest_half_angle_rad: 0.0,
            largest_thick_half_angle_rad: 0.0,
            lowest_resonance_hz: f64::INFINITY,
            highest_resonance_hz: 0.0,
            quality,
        };
        for i in 0..n {
            let (left, right) = (positions[i], positions[i + 1]);
            let width = right - left;
            let band = cell_band(duct, &starts, left, right)?;
            // The first cell carries the bore's lumped cup (extra gas volume with no wall).
            let expected = band.volume + if i == 0 { inlet_volume_m3 } else { 0.0 };
            let agrees = (expected - volumes[i]).abs() <= VOLUME_AGREEMENT * volumes[i];
            if !agrees || !positive(band.lateral) || !positive(band.slant) || !positive(width) {
                return Err(unsupported("the duct is not the bore's own duct (cell volume)"));
            }
            let mean_radius = band.lateral / (core::f64::consts::TAU * width);
            let cosine = (band.lateral / band.slant).min(1.0);
            let strip = law.strip(mean_radius)?;
            if !positive(strip.surface_density_kg_m2) || !positive(strip.stiffness_per_area_pa_m) {
                return Err(unsupported(
                    "shell strip needs positive finite surface density and stiffness",
                ));
            }
            let mass = strip.surface_density_kg_m2 * band.slant;
            let flare = if strip.membrane { cosine } else { 1.0 };
            let stiffness = strip.stiffness_per_area_pa_m * band.lateral * flare;
            let omega0 = det::sqrt(stiffness / mass);
            let damping = mass * omega0 / quality;
            // exp([[0, 1, 0], [-k/m, -c/m, 1/m], [0, 0, 0]] dt) acts on (s, v, F) with F held: the exact step.
            let generator = [
                0.0,
                dt,
                0.0,
                -dt * stiffness / mass,
                -dt * damping / mass,
                dt / mass,
                0.0,
                0.0,
                0.0,
            ];
            let e = exponential(&generator, 3)
                .map_err(|_| unsupported("shell step exponential is not representable"))?;
            if e.iter().any(|x| !x.is_finite()) {
                return Err(unsupported("shell step exponential is not representable"));
            }
            let angle = det::acos(cosine);
            receipt.largest_half_angle_rad = receipt.largest_half_angle_rad.max(angle);
            if strip.membrane {
                receipt.membrane_cells += 1;
            } else {
                receipt.thick_cells += 1;
                receipt.largest_thick_half_angle_rad =
                    receipt.largest_thick_half_angle_rad.max(angle);
            }
            let hz = omega0 / core::f64::consts::TAU;
            receipt.lowest_resonance_hz = receipt.lowest_resonance_hz.min(hz);
            receipt.highest_resonance_hz = receipt.highest_resonance_hz.max(hz);
            cells.push(ShellCell {
                volume: volumes[i],
                lateral: band.lateral,
                mass,
                stiffness,
                step: [e[0], e[1], e[3], e[4], e[2], e[5]],
            });
        }
        Ok((
            Self {
                cells,
                state: vec![[0.0, 0.0]; n],
                next: vec![[0.0, 0.0]; n],
                gas_work: 0.0,
                room_work: 0.0,
                dissipated: 0.0,
                receipt,
            },
            receipt,
        ))
    }

    /// Move every wall under the cell pressures and rescale the gas for the volume it now fills. A refusal
    /// leaves the field and the walls unchanged.
    pub(super) fn apply(
        &mut self,
        solver: &mut Euler1d,
        room_pressure: f64,
    ) -> Result<(), NonlinearDuctError> {
        let Self {
            cells,
            state,
            next,
            gas_work,
            room_work,
            dissipated,
            ..
        } = self;
        let (mut work, mut room, mut lost) = (0.0, 0.0, 0.0);
        solver.apply_cell_sources(
            |i, q, p| {
                let cell = &cells[i];
                let [s, v] = state[i];
                let before = cell.volume + cell.lateral * s;
                let gas_stiffness = p.density * p.sound_speed * p.sound_speed;
                let force = cell.midpoint_force(state[i], p.pressure, room_pressure, gas_stiffness, before);
                let [s1, v1] = cell.advance(state[i], force);
                let ds = s1 - s;
                let dv = cell.lateral * ds;
                let after = before + dv;
                if !(force.is_finite() && s1.is_finite() && v1.is_finite() && after.is_finite() && after > 0.0) {
                    return Err(GasError::Invalid("shell displacement left the bore volume"));
                }
                // The pressure the held force stands for: the mean over the step (module doc).
                let mid = room_pressure + force / cell.lateral;
                let stored = |s: f64, v: f64| 0.5 * (cell.mass * v * v + cell.stiffness * s * s);
                lost += force * ds - (stored(s1, v1) - stored(s, v));
                work += mid * dv;
                room += room_pressure * dv;
                next[i] = [s1, v1];
                let ratio = before / after;
                Ok(Conserved {
                    density: q.density * ratio,
                    momentum: q.momentum * ratio,
                    energy: (q.energy * before - mid * dv) / after,
                })
            },
            || false,
        )?;
        core::mem::swap(state, next);
        *gas_work += work;
        *room_work += room;
        *dissipated += lost;
        Ok(())
    }

    pub(super) fn ledger(&self) -> ShellLedger {
        let stored = self
            .cells
            .iter()
            .zip(&self.state)
            .map(|(cell, [s, v])| 0.5 * (cell.mass * v * v + cell.stiffness * s * s))
            .sum();
        ShellLedger {
            gas_work_j: self.gas_work,
            room_work_j: self.room_work,
            dissipated_j: self.dissipated,
            stored_j: stored,
        }
    }
}

// Volume, projected lateral area and slant area of the bore between axial positions `left` and `right`.
fn cell_band(
    duct: &Duct,
    starts: &[f64],
    left: f64,
    right: f64,
) -> Result<Band, NonlinearDuctError> {
    let mut band = Band {
        volume: 0.0,
        lateral: 0.0,
        slant: 0.0,
    };
    for (k, segment) in duct.segments.iter().enumerate() {
        let a = left.max(starts[k]) - starts[k];
        let b = right.min(starts[k + 1]) - starts[k];
        if b <= a {
            continue;
        }
        let b = b.min(segment.axial_length());
        let width = b - a;
        match segment {
            Segment::Cylinder { radius: r, .. } => {
                let side = core::f64::consts::TAU * r * width;
                band.volume += core::f64::consts::PI * r * r * width;
                band.lateral += side;
                band.slant += side;
            }
            Segment::Cone { .. } => {
                let (ra, rb) = (radius(segment, a)?, radius(segment, b)?);
                let mean = core::f64::consts::PI * (ra + rb) * width;
                band.volume +=
                    core::f64::consts::PI * (width / 3.0) * (ra * ra + ra * rb + rb * rb);
                band.lateral += mean;
                band.slant += core::f64::consts::PI
                    * (ra + rb)
                    * det::sqrt(width * width + (rb - ra) * (rb - ra));
            }
            Segment::Profile { .. } => {
                let h = width / PROFILE_PANELS as f64;
                let length = segment.axial_length();
                let probe = h / 8.0;
                let (mut volume, mut lateral, mut slant) = (0.0, 0.0, 0.0);
                for j in 0..=PROFILE_PANELS {
                    let x = if j == PROFILE_PANELS { b } else { a + h * j as f64 };
                    let r = radius(segment, x)?;
                    let (lo, hi) = ((x - probe).max(0.0), (x + probe).min(length));
                    let slope = (radius(segment, hi)? - radius(segment, lo)?) / (hi - lo);
                    let weight = if j == 0 || j == PROFILE_PANELS {
                        1.0
                    } else if j % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    volume += weight * core::f64::consts::PI * r * r;
                    lateral += weight * core::f64::consts::TAU * r;
                    slant += weight * core::f64::consts::TAU * r * det::sqrt(1.0 + slope * slope);
                }
                band.volume += volume * h / 3.0;
                band.lateral += lateral * h / 3.0;
                band.slant += slant * h / 3.0;
            }
            Segment::ToneHole { .. } => {
                return Err(unsupported("tone holes require a conservative gas junction"));
            }
        }
    }
    Ok(band)
}

// ---- a wall that couples the cells: any set of wall modes (the modal bell wall of fs-couple), one `WallDriver` ----

/// What a [`WallDriver`] is given at the start of each accepted step.
#[derive(Debug, Clone, Copy)]
pub struct WallStepInput<'a> {
    /// Gas pressure of every cell [Pa].
    pub pressure: &'a [f64],
    /// `rho c^2` of every cell [Pa]: the gas's stiffness against a volume change (adiabatic, frozen composition).
    pub gas_stiffness: &'a [f64],
    /// Physical volume of every cell before the step [m^3].
    pub volume: &'a [f64],
    /// The rest (room) pressure [Pa].
    pub room_pressure: f64,
}

/// A wall that moves under the pressures of ALL the cells at once and sweeps volume from each. The ring oscillator of
/// [`NonlinearBore::bind_shell`](super::NonlinearBore::bind_shell) is the case where each cell has its own; a driver
/// covers walls whose motion is shared (the modes of a shell). The bore keeps the gas side: the rescale of every
/// cell's average for the volume it now fills, the work `p_mid dV`, and the account.
///
/// A driver stages a step, and the bore commits it only after the gas field accepted the update, so a refused update
/// changes nothing.
pub trait WallDriver {
    /// The number of bore cells it covers.
    fn cells(&self) -> usize;
    /// Stage one step (the bore's fixed step) under `input`: write into `swept[i]` the change of cell `i`'s swept
    /// volume over the step [m^3] and into `mean_pressure[i]` the pressure that cell's work is paid at [Pa] (the
    /// mean of its pressure before and after the step). The staged state is not published.
    ///
    /// # Errors
    /// A named refusal; nothing is changed.
    fn stage(
        &mut self,
        input: &WallStepInput<'_>,
        swept: &mut [f64],
        mean_pressure: &mut [f64],
    ) -> Result<(), NonlinearDuctError>;
    /// Publish the staged step.
    fn commit(&mut self);
    /// Energy stored in the wall now [J].
    fn stored_energy_j(&self) -> f64;
    /// Work its dampers have dissipated so far [J].
    fn dissipated_j(&self) -> f64;
}

pub(super) struct DrivenWall {
    driver: Box<dyn WallDriver>,
    nominal: Vec<f64>,
    swept: Vec<f64>,
    gas_work: f64,
    room_work: f64,
    pressure: Vec<f64>,
    stiffness: Vec<f64>,
    volume: Vec<f64>,
    change: Vec<f64>,
    mean: Vec<f64>,
}

impl DrivenWall {
    pub(super) fn new(driver: Box<dyn WallDriver>, nominal: Vec<f64>) -> Result<Self, NonlinearDuctError> {
        let n = nominal.len();
        if driver.cells() != n {
            return Err(unsupported("the wall driver covers a different number of cells than the bore"));
        }
        Ok(Self {
            driver,
            swept: vec![0.0; n],
            gas_work: 0.0,
            room_work: 0.0,
            pressure: vec![0.0; n],
            stiffness: vec![0.0; n],
            volume: nominal.clone(),
            change: vec![0.0; n],
            mean: vec![0.0; n],
            nominal,
        })
    }

    /// Advance the wall one step under the cell pressures and rescale every cell for the volume it now fills. A refusal
    /// (the driver's, a cell swept to no volume, or the gas owner's) leaves the field, the wall and the account unchanged.
    pub(super) fn apply(
        &mut self,
        solver: &mut Euler1d,
        gas: &FrozenGas,
        room_pressure: f64,
    ) -> Result<(), NonlinearDuctError> {
        let states: Vec<Conserved> = solver.states().to_vec();
        for (i, q) in states.iter().enumerate() {
            let p = gas.primitive(*q, || false)?;
            self.pressure[i] = p.pressure;
            self.stiffness[i] = p.density * p.sound_speed * p.sound_speed;
            self.volume[i] = self.nominal[i] + self.swept[i];
        }
        let input = WallStepInput {
            pressure: &self.pressure,
            gas_stiffness: &self.stiffness,
            volume: &self.volume,
            room_pressure,
        };
        self.driver.stage(&input, &mut self.change, &mut self.mean)?;
        let mut next = Vec::with_capacity(states.len());
        let (mut work, mut room) = (0.0, 0.0);
        for (i, q) in states.iter().enumerate() {
            let (before, dv, mid) = (self.volume[i], self.change[i], self.mean[i]);
            let after = before + dv;
            if !(dv.is_finite() && mid.is_finite() && after.is_finite() && after > 0.0) {
                return Err(GasError::Invalid("wall displacement left the bore volume").into());
            }
            let ratio = before / after;
            next.push(Conserved {
                density: q.density * ratio,
                momentum: q.momentum * ratio,
                energy: (q.energy * before - mid * dv) / after,
            });
            work += mid * dv;
            room += room_pressure * dv;
        }
        solver.apply_cell_sources(|i, _, _| Ok(next[i]), || false)?;
        self.driver.commit();
        for (s, dv) in self.swept.iter_mut().zip(&self.change) {
            *s += dv;
        }
        self.gas_work += work;
        self.room_work += room;
        Ok(())
    }

    pub(super) fn ledger(&self) -> ShellLedger {
        ShellLedger {
            gas_work_j: self.gas_work,
            room_work_j: self.room_work,
            dissipated_j: self.driver.dissipated_j(),
            stored_j: self.driver.stored_energy_j(),
        }
    }
}
