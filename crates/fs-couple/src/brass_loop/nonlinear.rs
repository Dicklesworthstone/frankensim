//! The lip on the shared nonlinear gas bore (bead `frankensim-giuhr`).
//!
//! The same outward-striking lip island as [`super::BrassVoice`] meters the
//! mouth into a bore whose gas obeys the conservative Euler equations of the
//! performance's own NASA-9 mixture ([`NonlinearBore`]). Wave steepening, the
//! Zwikker-Kosten wall memory and the baffled-piston bell all live in that one
//! field; no reflectance line, transfer function or fitted load is used.
//!
//! The lip channel is an isentropic nozzle of area `width * gap` from the
//! mouth (stagnation gas at the blowing pressure and the room temperature) to
//! the pressure on the bore's first face. Its jet mixes out over that face
//! without pressure recovery, the linear loop's Bernoulli law made
//! compressible. When the cup pressure exceeds the mouth's, gas leaves through
//! the same channel. At every Runge-Kutta stage the face pressure solves mass
//! continuity between the metered channel and the interior's own end wave
//! (this source's shock or fan). A closed lip is a wall.
//!
//! The lip steps once per gas step under the cup pressure of the previous
//! step's first stage, and the gap is held through both stages. The radiated
//! observer receives every mechanics sample's mean bell volume flow.
//!
//! Limits, stated rather than hidden: the mouth is a fixed reservoir at the
//! room's temperature and composition (no breath heating or exhaled CO2); the
//! mixed-out inflow enters at the mouth temperature (its kinetic heat,
//! `u^2 / (2 c_p)` with `u` the cup-face speed, is not added); reverse flow
//! discharges the cup's static state; one held bore without valves; the bell
//! radiates as a rigid baffled piston.

use fs_duct::Duct;
use fs_duct::nonlinear::bore::{BoreSpec, NonlinearBore, WallGas, WallMemory, WallReceipt};
use fs_duct::nonlinear::{FaceFlux, FrozenGas, PortStage, Primitive};
use fs_gas::dynamics::Error as GasError;
use fs_phs::linear_exact::{LinearPropagator, LinearSplit, step_split};
use fs_phs::{PortHamiltonian, mass_spring_damper};

use super::{BrassBlockDiag, BrassControl, LipIsland, prepare_lip, zero_crossing_lock_hz};

/// How the nonlinear bore is cut, clocked and walled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearClock {
    /// Largest axial gas cell [m].
    pub maximum_cell_m: f64,
    /// Output sample rate [Hz].
    pub output_rate_hz: u32,
    /// Mechanics (observer) samples per output sample, 1..=16.
    pub mechanics_per_output: u32,
    /// Gas steps per mechanics sample.
    pub gas_steps_per_mechanics: u32,
    /// Wall memory of every cell.
    pub wall: WallMemory,
}

impl NonlinearClock {
    /// Mechanics (observer) clock [Hz].
    #[must_use]
    pub fn mechanics_rate_hz(&self) -> u32 {
        self.output_rate_hz * self.mechanics_per_output
    }

    /// Fixed gas step [s].
    #[must_use]
    pub fn gas_step_s(&self) -> f64 {
        1.0 / (f64::from(self.mechanics_rate_hz()) * f64::from(self.gas_steps_per_mechanics))
    }
}

/// The face the lip channel and the bore agree on at one stage.
#[derive(Debug, Clone, Copy)]
struct LipFace {
    flux: FaceFlux,
    pressure_pa: f64,
    mass_flow_kg_s: f64,
}

/// The lip island driving a nonlinear gas bore.
pub struct NonlinearLip {
    bore: NonlinearBore,
    gas: FrozenGas,
    room: Primitive,
    lip_sys: PortHamiltonian,
    lip_linear: LinearSplit,
    lip_state: Vec<f64>,
    lip: LipIsland,
    blow_pa: f64,
    mouth: Primitive,
    clock: NonlinearClock,
    face_pa: f64,
    block_index: u64,
    diags: Vec<BrassBlockDiag>,
}

impl NonlinearLip {
    /// Fill `duct` with the room gas at rest and attach the lip at its inlet.
    /// `cup_volume_m3` is added to the first gas cell.
    ///
    /// # Errors
    /// A named refusal from the lip island, the gas owner, the mesh, the wall
    /// binding or the clock.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        duct: &Duct,
        gas: FrozenGas,
        room_temperature_k: f64,
        room_pressure_pa: f64,
        wall_gas: WallGas,
        lip: LipIsland,
        cup_volume_m3: f64,
        clock: NonlinearClock,
    ) -> Result<Self, String> {
        if !(1..=16).contains(&clock.mechanics_per_output) || clock.gas_steps_per_mechanics == 0 {
            return Err(
                "nonlinear bore needs 1..=16 mechanics samples per output sample and at least one gas step per mechanics sample"
                    .into(),
            );
        }
        if !(lip.mass_kg > 0.0
            && lip.stiffness_n_m > 0.0
            && lip.damping_n_s_m >= 0.0
            && lip.width_m > 0.0
            && lip.rest_gap_m > 0.0
            && lip.face_area_m2 > 0.0)
        {
            return Err("lip island parameters must be positive".into());
        }
        let bore = NonlinearBore::new(
            duct,
            gas.clone(),
            room_temperature_k,
            room_pressure_pa,
            Some(wall_gas),
            BoreSpec {
                maximum_cell_m: clock.maximum_cell_m,
                inlet_volume_m3: cup_volume_m3,
                step_s: clock.gas_step_s(),
                cfl: 0.45,
                wall: Some(clock.wall),
            },
        )
        .map_err(|e| format!("nonlinear bore refused: {e}"))?;
        let lip_sys = mass_spring_damper(lip.mass_kg, lip.stiffness_n_m, lip.damping_n_s_m)
            .map_err(|_| "lip island refused admission".to_string())?;
        let lip_linear = prepare_lip(&lip_sys).map_err(|e| e.to_string())?;
        let room = bore.room();
        Ok(Self {
            bore,
            gas,
            room,
            lip_sys,
            lip_linear,
            lip_state: vec![0.0, 0.0],
            lip,
            blow_pa: 0.0,
            mouth: room,
            clock,
            face_pa: room.pressure,
            block_index: 0,
            diags: Vec::new(),
        })
    }

    /// Apply a control-rate change between blocks. A held bore has no valves.
    ///
    /// # Errors
    /// An unusable value, a refused lip rebuild, or a valve change.
    pub fn apply(&mut self, control: BrassControl) -> Result<(), String> {
        match control {
            BrassControl::SetBlowingPressure(p) => {
                if !(p.is_finite() && p >= 0.0) {
                    return Err("blowing pressure must be finite non-negative".into());
                }
                let pressure = self.room.pressure + p;
                let density = pressure / (self.gas.gas_constant() * self.room.temperature);
                let state = self
                    .gas
                    .conserved(density, 0.0, pressure)
                    .map_err(|e| format!("mouth reservoir refused: {e}"))?;
                self.mouth = self
                    .gas
                    .primitive(state, || false)
                    .map_err(|e| format!("mouth reservoir refused: {e}"))?;
                self.blow_pa = p;
            }
            BrassControl::SetLipTension(t) => {
                if !(t.is_finite() && t > 0.0) {
                    return Err("lip tension must be finite positive".into());
                }
                // As in the linear loop: stiffness ~ t^2, damping ~ t, same state.
                let lip_sys = mass_spring_damper(
                    self.lip.mass_kg,
                    self.lip.stiffness_n_m * t * t,
                    self.lip.damping_n_s_m * t,
                )
                .map_err(|_| "tensioned lip island refused admission".to_string())?;
                self.lip_linear = prepare_lip(&lip_sys).map_err(|e| e.to_string())?;
                self.lip_sys = lip_sys;
            }
            BrassControl::SetValve { .. } => {
                return Err("the nonlinear bore is one held bore without valves".into());
            }
        }
        Ok(())
    }

    /// Per-block diagnostics so far.
    #[must_use]
    pub fn diagnostics(&self) -> &[BrassBlockDiag] {
        &self.diags
    }

    /// The bore's gas field and its receipts.
    #[must_use]
    pub const fn bore(&self) -> &NonlinearBore {
        &self.bore
    }

    /// The wall memory's admission receipt.
    #[must_use]
    pub fn wall_receipt(&self) -> Option<WallReceipt> {
        self.bore.wall_receipt()
    }

    /// The declared clock.
    #[must_use]
    pub const fn clock(&self) -> NonlinearClock {
        self.clock
    }

    /// Render one block: `mouthpiece` receives the cup-face overpressure [Pa]
    /// at each output sample, `flows` the mean bell volume flow [m^3/s] of
    /// every mechanics sample, `mechanics_per_output` per output sample.
    ///
    /// # Errors
    /// A length mismatch, or any refusal of the lip, the lip face or the gas
    /// bore. A refusal leaves the block partly advanced and must end the run.
    pub fn step_block(&mut self, mouthpiece: &mut [f64], flows: &mut [f64]) -> Result<(), String> {
        let ratio = self.clock.mechanics_per_output as usize;
        if mouthpiece.is_empty() || flows.len() != mouthpiece.len() * ratio {
            return Err("nonlinear lip block needs one flow per mechanics sample".into());
        }
        let steps = self.clock.gas_steps_per_mechanics as usize;
        let dt = self.clock.gas_step_s();
        let mut gap_min = f64::INFINITY;
        let mut gap_max = 0.0f64;
        let mut gap_acc = 0.0f64;
        let mut blow_work = 0.0f64;
        let mut bore_work = 0.0f64;
        let mut worst_defect = 0.0f64;
        for (sample, slot) in mouthpiece.iter_mut().enumerate() {
            for frame in 0..ratio {
                let mut flow = 0.0;
                for _ in 0..steps {
                    let (gap, defect, face, outlet) = self.gas_step(dt)?;
                    gap_min = gap_min.min(gap);
                    gap_max = gap_max.max(gap);
                    gap_acc += gap;
                    worst_defect = worst_defect.max(defect);
                    let mouth_volume = face.mass_flow_kg_s.max(0.0) / self.mouth.density;
                    blow_work += self.blow_pa * mouth_volume * dt;
                    bore_work += (face.pressure_pa - self.room.pressure) * face.mass_flow_kg_s
                        / self.room.density
                        * dt;
                    flow += outlet;
                }
                flows[sample * ratio + frame] = flow / steps as f64;
            }
            *slot = self.face_pa - self.room.pressure;
        }
        let output_dt = 1.0 / f64::from(self.clock.output_rate_hz);
        let n = (mouthpiece.len() * ratio * steps) as f64;
        self.diags.push(BrassBlockDiag {
            block: self.block_index,
            lock_hz: zero_crossing_lock_hz(mouthpiece, output_dt),
            gap_min_m: gap_min,
            gap_mean_m: gap_acc / n,
            gap_max_m: gap_max,
            blow_work_j: blow_work,
            bore_work_j: bore_work,
            supply_defect: worst_defect,
        });
        self.block_index += 1;
        Ok(())
    }

    // One gas step: lip under the last cup pressure, then the bore with the
    // lip face at both stages. Returns the gap, the lip supply defect, the
    // first-stage face and the mean bell volume flow.
    fn gas_step(&mut self, dt: f64) -> Result<(f64, f64, LipFace, f64), String> {
        let force = (self.blow_pa - (self.face_pa - self.room.pressure)) * self.lip.face_area_m2;
        let record = step_split(
            &self.lip_sys,
            &self.lip_linear,
            LinearPropagator::ExactZoh,
            &self.lip_state,
            &[force],
            dt,
        )
        .map_err(|_| "lip island step refused".to_string())?;
        let defect = record.supply_defect().abs();
        self.lip_state = record.x;
        let gap = (self.lip.rest_gap_m + self.lip_state[0]).max(0.0);
        let opening = self.lip.width_m * gap;
        let inlet = self.bore.inlet_area_m2();
        let (gas, mouth) = (&self.gas, self.mouth);
        let mut first = None;
        let step = self
            .bore
            .step(|stage: &PortStage| {
                let face = lip_face(gas, stage.left.1, mouth, opening, inlet)?;
                if stage.stage == 0 {
                    first = Some(face);
                }
                Ok(Some(face.flux))
            })
            .map_err(|e| format!("nonlinear bore step refused: {e}"))?;
        let face = first.ok_or_else(|| "the bore offered no first stage".to_string())?;
        self.face_pa = face.pressure_pa;
        Ok((gap, defect, face, step.outlet_volume_flow_m3_s))
    }
}

// Signed mass flow the lip channel meters at cup-face pressure `p` [kg/s]:
// mouth to cup below the mouth pressure, the cup's static gas back out above it.
fn metered(
    gas: &FrozenGas,
    star: Primitive,
    mouth: Primitive,
    opening: f64,
    p: f64,
) -> Result<f64, GasError> {
    if opening == 0.0 || p == mouth.pressure {
        Ok(0.0)
    } else if p < mouth.pressure {
        let (jet, _) = gas.nozzle_exit(mouth, p, || false)?;
        Ok(jet.density * jet.velocity * opening)
    } else {
        let donor = Primitive {
            velocity: 0.0,
            ..star
        };
        let (jet, _) = gas.nozzle_exit(donor, mouth.pressure, || false)?;
        Ok(-jet.density * jet.velocity * opening)
    }
}

// Mass flux through the face if the interior's end wave reaches `p`: mouth gas
// at the mouth temperature entering, the interior's star state leaving.
fn face_mass_flux(gas: &FrozenGas, star: Primitive, mouth: Primitive, p: f64) -> f64 {
    if star.velocity >= 0.0 {
        p / (gas.gas_constant() * mouth.temperature) * star.velocity
    } else {
        star.density * star.velocity
    }
}

fn lip_face(
    gas: &FrozenGas,
    interior: Primitive,
    mouth: Primitive,
    opening: f64,
    inlet_area: f64,
) -> Result<LipFace, GasError> {
    // The pressure at which the interior's end wave stops the gas: the face
    // pressure of a closed lip, and one end of the bracket for an open one.
    let mut wall = interior.pressure;
    let mut stopped = false;
    for _ in 0..64 {
        let wave = gas.end_wave(interior, -1.0, wall, || false)?;
        let step = wave.primitive.velocity / wave.velocity_per_pressure;
        wall -= step;
        if !(wall > 0.0 && wall.is_finite()) {
            return Err(GasError::Invalid(
                "lip face closed-end pressure left the gas domain",
            ));
        }
        if step.abs() <= 1e-14 * wall {
            stopped = true;
            break;
        }
    }
    if !stopped {
        return Err(GasError::Invalid(
            "lip face closed-end pressure did not converge",
        ));
    }
    // Continuity G(p) = metered(p) / A - face flux(p) falls strictly with p and
    // changes sign between the closed-end and mouth pressures.
    let residual = |p: f64| -> Result<(f64, Primitive, f64), GasError> {
        let wave = gas.end_wave(interior, -1.0, p, || false)?;
        let mass = metered(gas, wave.primitive, mouth, opening, p)?;
        let value = mass / inlet_area - face_mass_flux(gas, wave.primitive, mouth, p);
        Ok((value, wave.primitive, wave.signal_speed))
    };
    let mut pressure = wall;
    let (mut lo, mut hi) = (wall.min(mouth.pressure), wall.max(mouth.pressure));
    if opening > 0.0 && lo < hi {
        // Illinois false position on the bracket [lo, hi], G(lo) >= 0 >= G(hi).
        let (mut g_lo, mut g_hi) = (residual(lo)?.0, residual(hi)?.0);
        let mut side = 0i8;
        let mut converged = false;
        for _ in 0..200 {
            let next = if g_lo > g_hi {
                (lo * g_hi - hi * g_lo) / (g_hi - g_lo)
            } else {
                lo.midpoint(hi)
            };
            let next = if next > lo && next < hi {
                next
            } else {
                lo.midpoint(hi)
            };
            let g = residual(next)?.0;
            pressure = next;
            if g == 0.0 {
                converged = true;
                break;
            }
            if g > 0.0 {
                (lo, g_lo) = (next, g);
                if side == 1 {
                    g_hi *= 0.5;
                }
                side = 1;
            } else {
                (hi, g_hi) = (next, g);
                if side == -1 {
                    g_lo *= 0.5;
                }
                side = -1;
            }
            if hi - lo <= 1e-12 * hi {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(GasError::Invalid("lip face continuity did not converge"));
        }
    }
    let (_, star, signal) = residual(pressure)?;
    let flux = if star.velocity >= 0.0 {
        let density = pressure / (gas.gas_constant() * mouth.temperature);
        let face = gas.primitive(gas.conserved(density, star.velocity, pressure)?, || false)?;
        gas.inflow_flux(face, mouth.temperature, signal)?
    } else {
        gas.face_flux(star, signal)?
    };
    Ok(LipFace {
        flux,
        pressure_pa: pressure,
        mass_flow_kg_s: flux.mass * inlet_area,
    })
}

#[cfg(test)]
#[path = "nonlinear_tests.rs"]
mod tests;
