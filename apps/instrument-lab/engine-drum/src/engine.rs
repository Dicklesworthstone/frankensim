//! The live voice: fs-couple's prepared exact-ZOH linear impact image, the
//! player's stick lifecycle, the energy ledger and Rayleigh radiation.
//!
//! Mechanics per step are EXACTLY `LinearImpactSystem::step` (the offline
//! `drum-modal` image). This module adds no integrator of its own. What it
//! owns is the honest bookkeeping around that owner in a never-ending live
//! session:
//!
//! - **One stick.** A strike launches the stick 0.2 mm above the current
//!   membrane surface with the declared downward speed (kinetic energy, not a
//!   force pulse). After the tip has rebounded and separated by `catch_gap_m`
//!   the player catches it: the stick leaves the mechanical system and its
//!   remaining kinetic energy is booked as recovered by the player.
//! - **Re-strikes.** A strike requested while the tip is in contact waits
//!   until the tip separates (no stored contact energy is ever deleted); a
//!   strike requested while the stick is in free flight replaces it at once.
//! - **Recompilation.** Launch and catch recompile the linear image with the
//!   accepted head states carried over bit-for-bit. The sealed-air spring and
//!   the elastic contact have no hidden history, so nothing is reset.

use crate::instrument::{HERTZ_ALPHA, Instrument, LAUNCH_GAP_M};
use crate::params::{OUTPUT_RATE_HZ, Refusal};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample};
use fs_couple::render::plate::impact::linear::{LinearImpactConfig, LinearImpactSystem, VolumeConnection};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, VolumeSpring};
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

/// A requested stroke: physical station and launch speed only.
#[derive(Debug, Clone, PartialEq)]
pub struct Strike {
    /// Station on the batter head [m], head-centred.
    pub position_m: [f64; 2],
    /// Downward launch speed of the tip [m/s].
    pub speed_m_s: f64,
    /// Batter modal participation at the station [1/sqrt(kg)].
    row: Vec<f64>,
}

/// The stick while it is part of the mechanical system (coordinate 0).
#[derive(Debug, Clone)]
struct Stick {
    position_m: [f64; 2],
    row: Vec<f64>,
}

/// Whole-session energy accounting [J]. Nothing here is corrected or clamped.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Ledger {
    /// Energy stored now: head strain/tension + kinetic, air spring, contact
    /// potential and the stick's kinetic energy while it is in the system.
    pub stored_j: f64,
    /// Cumulative modal-damping plus contact dissipation.
    pub dissipated_j: f64,
    /// Cumulative stick kinetic energy supplied by the player at launches.
    pub supplied_j: f64,
    /// Cumulative stick kinetic energy taken back by the player at catches.
    pub recovered_j: f64,
    /// stored - (supplied - recovered - dissipated): the uncorrected residual.
    pub balance_residual_j: f64,
    /// Largest single-step residual reported by the mechanics owner.
    pub max_step_residual_j: f64,
}

/// What `strike` did with the request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrikeOutcome {
    /// Launches at the next mechanical step.
    Queued,
    /// The tip is in contact; launches the moment it separates.
    DeferredUntilSeparation,
}

/// Read-only view for drawing. All values are accepted physical state.
#[derive(Debug, Clone, Default)]
pub struct Snapshot {
    /// Audio time rendered since the last (re)build [s].
    pub time_s: f64,
    /// Accepted mechanical steps since the last (re)build.
    pub mechanical_steps: u64,
    /// Energy ledger.
    pub ledger: Ledger,
    /// Stick is currently a body of the mechanical system.
    pub stick_present: bool,
    /// Tip is compressing the head (closure > 0).
    pub stick_in_contact: bool,
    /// A strike is waiting for launch.
    pub strike_pending: bool,
    /// Tip station [m] of the current (or last) strike.
    pub stick_xy_m: [f64; 2],
    /// Tip displacement, positive downward, 0 = undeformed head plane [m].
    pub stick_position_m: f64,
    /// Tip velocity, positive downward [m/s].
    pub stick_velocity_m_s: f64,
    /// Tip-into-membrane penetration [m] (0 when separated).
    pub penetration_m: f64,
    /// Elastic Hertz force K p^1.5 at the accepted endpoint [N].
    pub contact_force_n: f64,
    /// Sealed-air overpressure (K/V) * sum(area_i q_i) [Pa].
    pub cavity_pressure_pa: f64,
    /// Strikes launched since the last (re)build.
    pub strikes_launched: u64,
    /// Largest |pressure| at each microphone since the previous snapshot [Pa].
    pub peak_pa: [f64; 2],
    /// Output frames whose |pressure| exceeded full scale since the (re)build.
    pub clipped_frames: u64,
    /// Modal displacements, batter then resonant [m sqrt(kg)].
    pub q: Vec<f64>,
    /// Modal velocities, batter then resonant [m sqrt(kg)/s].
    pub v: Vec<f64>,
}

/// Counters from a benchmark render; wall time is measured by the caller.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BenchReport {
    /// Audio frames rendered.
    pub frames: usize,
    /// Mechanical steps taken.
    pub mechanical_steps: u64,
    /// Steps during which the stick was a body of the system.
    pub steps_with_stick: u64,
    /// Steps during which the tip was in contact.
    pub steps_in_contact: u64,
    /// Strikes launched.
    pub strikes: u64,
    /// Peak |pressure| at the batter microphone [Pa].
    pub peak_pa: f64,
    /// Final ledger.
    pub ledger: Ledger,
}

/// The reference's prepared work envelope (examples/percussion/mechanics.rs
/// `prepared_config` over main.rs `config`), at an explicit integer-Hz clock.
/// Only `max_modes`, `max_setup_terms` and the lifetime step budget differ:
/// they are admission ceilings for a wider optional window and an open-ended
/// session, not physical or numerical parameters.
#[must_use]
pub fn mechanics_config(rate_hz: u32, total_modes: usize) -> LinearImpactConfig {
    let maximum_energy_j = 20.0;
    LinearImpactConfig {
        sample_rate_hz: rate_hz,
        max_steps: 1u64 << 53,
        maximum_generalized_force: 1e6,
        component: ModalAcousticTimeBudget {
            maximum_total_energy_j: maximum_energy_j,
            ..ModalAcousticTimeBudget::audible_reference()
        },
        coupling: ModalCouplingConfig {
            max_modes: total_modes.max(64),
            max_connections: 8,
            max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: maximum_energy_j,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-9,
            energy_relative_tolerance: 1e-6,
        },
        contact: ModalContactConfig {
            max_iterations: 100,
            maximum_force_n: 1e4,
            maximum_penetration_m: 0.01,
            force_absolute_tolerance_n: 1e-9,
            force_relative_tolerance: 1e-9,
        },
        multiple: MultiContactConfig { max_contacts: 32, max_sweeps: 100, max_setup_terms: 100_000 },
    }
}

fn mechanics_refusal(error: impl core::fmt::Display) -> Refusal {
    Refusal::new("mechanics_refused", error.to_string())
}

/// Compile the instrument into fs-couple's prepared linear impact image.
///
/// Coordinate order matches the reference: `[stick?, batter modes, resonant
/// modes]`. `head_states` carries the batter then resonant coordinates.
/// `stick` supplies the tip's mass-normalized state and its batter row; the
/// contact row is `+1/sqrt(m)` on the stick and `-shape` on the batter modes,
/// so closure = tip displacement - membrane displacement, both positive down.
///
/// # Errors
/// Any refusal of the modal, coupling or contact owners.
pub fn compile_system(
    inst: &Instrument,
    rate_hz: u32,
    head_states: &[ModalAcousticState],
    stick: Option<(ModalAcousticState, &[f64])>,
    gate: &CancelGate,
) -> Result<LinearImpactSystem, Refusal> {
    let n0 = inst.heads[0].mode_count();
    let n1 = inst.heads[1].mode_count();
    if head_states.len() != n0 + n1 {
        return Err(Refusal::new("mechanics_refused", "head state does not match the retained modes"));
    }
    let lead = usize::from(stick.is_some());
    let total = lead + n0 + n1;
    let mut bodies = Vec::with_capacity(3);
    let mut contacts = Vec::new();
    if let Some((state, row)) = stick {
        if row.len() != n0 {
            return Err(Refusal::new("mechanics_refused", "strike row does not match the batter modes"));
        }
        let weight = 1.0 / inst.stick.effective_mass_kg.sqrt();
        bodies.push(ImpactBody {
            potential: BodyPotential::Linear(vec![0.0]),
            initial: vec![state],
            damping_per_s: vec![0.0],
        });
        let mut contact = vec![0.0; total];
        contact[0] = weight;
        for (i, shape) in row.iter().enumerate() {
            contact[1 + i] = -shape;
        }
        let obstacle = Obstacle::new(
            contact, 1, total, vec![0.0], vec![1.0], inst.stick.hertz_stiffness, HERTZ_ALPHA,
            "estimated isotropic Hertz tip: E_eff=0.8GPa, R=3mm; not identified hickory-shell contact".into(),
        )
        .and_then(|o| o.with_internal_loss(inst.params.contact_loss_s_m))
        .map_err(mechanics_refusal)?;
        contacts.push(obstacle);
    }
    let mut offset = 0;
    for head in &inst.heads {
        let count = head.mode_count();
        bodies.push(ImpactBody {
            potential: BodyPotential::Linear(head.omegas.clone()),
            initial: head_states[offset..offset + count].to_vec(),
            damping_per_s: head.damping_per_s.clone(),
        });
        offset += count;
    }
    let mut areas = vec![0.0; lead];
    areas.extend(inst.compression_areas());
    let volume = VolumeSpring {
        bulk_modulus_pa: inst.params.bulk_modulus_pa(),
        volume_m3: inst.params.volume_m3(),
        areas,
    };
    let reference_area_m2 = core::f64::consts::PI * inst.params.radius_m * inst.params.radius_m;
    LinearImpactSystem::new_with_dampers(
        bodies,
        contacts,
        vec![VolumeConnection { spring: volume, reference_area_m2 }],
        Vec::new(),
        mechanics_config(rate_hz, total),
        gate,
    )
    .map_err(mechanics_refusal)
}

fn microphone(inst: &Instrument, head: usize, rate_hz: u32) -> Result<BaffledPressure, Refusal> {
    let h = &inst.heads[head];
    let count = h.mode_count();
    // Outward normal velocity: the batter radiates upward (coordinates are
    // positive DOWNWARD, hence -phi); the resonant head radiates downward (+phi).
    let sign = if head == 0 { -1.0 } else { 1.0 };
    let mut surface = Vec::with_capacity(inst.nodes.len());
    for (node, point) in inst.nodes.iter().enumerate() {
        let shapes = &h.node_shapes[node * count..(node + 1) * count];
        // Fixed-rim nodes do not move and radiate nothing.
        if shapes.iter().all(|s| *s == 0.0) {
            continue;
        }
        surface.push(SurfaceSample {
            position_m: [point[0], point[1], 0.0],
            area_m2: inst.nodal_area_m2[node],
            mode_shape: shapes.iter().map(|s| sign * s).collect(),
        });
    }
    let position = if head == 0 { inst.params.mic_batter_m } else { inst.params.mic_resonant_m };
    BaffledPressure::new_multirate(
        &surface,
        count,
        rate_hz,
        OUTPUT_RATE_HZ,
        position,
        RayleighMedium { density: inst.params.gas_density_kg_m3, sound_speed: inst.params.gas_sound_speed_m_s },
    )
    .map_err(|e| Refusal::new("radiation_refused", e))
}

/// One running instrument state at one mechanical rate.
pub struct Voice {
    substeps: usize,
    rate_hz: u32,
    n0: usize,
    n1: usize,
    system: LinearImpactSystem,
    stick: Option<Stick>,
    pending: Option<Strike>,
    last_xy_m: [f64; 2],
    mics: [BaffledPressure; 2],
    traces: [Vec<f64>; 2],
    force: Vec<f64>,
    gate: CancelGate,
    ledger: Ledger,
    frames: u64,
    steps: u64,
    steps_with_stick: u64,
    steps_in_contact: u64,
    strikes: u64,
    peak_pa: [f64; 2],
    clipped: u64,
    halted: Option<Refusal>,
    sqrt_mass: f64,
    areas: Vec<f64>,
    pressure_per_volume: f64,
}

impl Voice {
    /// A silent instrument at rest, compiled at `48 kHz * substeps`.
    ///
    /// # Errors
    /// Clock admission or any owner refusal.
    pub fn new(inst: &Instrument, substeps: usize) -> Result<Self, Refusal> {
        if !matches!(substeps, 1 | 2 | 4 | 8 | 16) {
            return Err(Refusal::new("bad_parameter", "substeps must be 1, 2, 4, 8 or 16"));
        }
        let rate_hz = OUTPUT_RATE_HZ * substeps as u32;
        if inst.params.band_high_hz / f64::from(rate_hz) >= 0.45 {
            return Err(Refusal::new("nyquist_guard", "declared head window crosses the mechanical Nyquist guard"));
        }
        let n0 = inst.heads[0].mode_count();
        let n1 = inst.heads[1].mode_count();
        let gate = CancelGate::new_clock_free();
        let rest = vec![ModalAcousticState::default(); n0 + n1];
        let system = compile_system(inst, rate_hz, &rest, None, &gate)?;
        let mics = [microphone(inst, 0, rate_hz)?, microphone(inst, 1, rate_hz)?];
        Ok(Self {
            substeps,
            rate_hz,
            n0,
            n1,
            force: vec![0.0; system.mode_count()],
            system,
            stick: None,
            pending: None,
            last_xy_m: [0.0; 2],
            mics,
            traces: [vec![0.0; substeps * n0], vec![0.0; substeps * n1]],
            gate,
            ledger: Ledger::default(),
            frames: 0,
            steps: 0,
            steps_with_stick: 0,
            steps_in_contact: 0,
            strikes: 0,
            peak_pa: [0.0; 2],
            clipped: 0,
            halted: None,
            sqrt_mass: inst.stick.effective_mass_kg.sqrt(),
            areas: inst.compression_areas(),
            pressure_per_volume: inst.params.bulk_modulus_pa() / inst.params.volume_m3(),
        })
    }

    /// Mechanical steps per audio frame.
    #[must_use]
    pub const fn substeps(&self) -> usize {
        self.substeps
    }
    /// Mechanical rate [Hz].
    #[must_use]
    pub const fn mechanical_rate_hz(&self) -> u32 {
        self.rate_hz
    }
    /// The refusal that stopped this voice, if any. A halted voice is silent.
    #[must_use]
    pub const fn halted(&self) -> Option<&Refusal> {
        self.halted.as_ref()
    }
    /// Anti-alias (decimator) latency of each microphone, in 48 kHz frames.
    #[must_use]
    pub fn filter_delay_frames(&self) -> [f64; 2] {
        [self.mics[0].filter_delay_samples(), self.mics[1].filter_delay_samples()]
    }
    /// Acoustic travel-time support (nearest, farthest surface point) of each
    /// microphone, in 48 kHz frames. Physical propagation, not processing lag.
    #[must_use]
    pub fn travel_delay_frames(&self) -> [(usize, usize); 2] {
        [self.mics[0].delay_samples, self.mics[1].delay_samples]
    }
    /// Multiply-adds per audio frame spent in both microphones.
    #[must_use]
    pub fn radiation_multiply_adds(&self) -> usize {
        self.mics[0].multiply_adds_per_sample() + self.mics[1].multiply_adds_per_sample()
    }

    /// Request a stroke. Admission (station on the head, speed in range)
    /// happens here; the launch happens inside the mechanical clock.
    ///
    /// # Errors
    /// `strike_outside_head`, `strike_speed`, or `halted`.
    pub fn strike(&mut self, inst: &Instrument, position_m: [f64; 2], speed_m_s: f64) -> Result<StrikeOutcome, Refusal> {
        if let Some(refusal) = &self.halted {
            return Err(Refusal::new("halted", format!("the voice stopped on a refusal ({refusal}); reset or rebuild first")));
        }
        // The reference admits 0..=20 m/s; zero launches no stroke, so the
        // live engine asks for a strictly positive speed.
        if !speed_m_s.is_finite() || speed_m_s <= 0.0 || speed_m_s > 20.0 {
            return Err(Refusal::new("strike_speed", "strike speed must be finite in (0, 20] m/s"));
        }
        let row = inst.batter_shapes_at(position_m)?;
        let outcome = if self.stick.is_some() && self.closure().0 > 0.0 {
            StrikeOutcome::DeferredUntilSeparation
        } else {
            StrikeOutcome::Queued
        };
        self.pending = Some(Strike { position_m, speed_m_s, row });
        Ok(outcome)
    }

    /// (closure [m], tip velocity [m/s]) of the present stick; closure > 0 is
    /// penetration, closure < 0 is the open gap. (0, 0) without a stick.
    fn closure(&self) -> (f64, f64) {
        let Some(stick) = &self.stick else { return (0.0, 0.0) };
        let state = self.system.state();
        let tip = state[0] / self.sqrt_mass;
        let mut membrane = 0.0;
        for (i, shape) in stick.row.iter().enumerate() {
            membrane += shape * state[2 * (1 + i)];
        }
        (tip - membrane, state[1] / self.sqrt_mass)
    }

    fn head_states(&self) -> Vec<ModalAcousticState> {
        let lead = usize::from(self.stick.is_some());
        let state = self.system.state();
        (0..self.n0 + self.n1)
            .map(|i| ModalAcousticState {
                displacement_m_sqrt_kg: state[2 * (lead + i)],
                velocity_m_sqrt_kg_per_s: state[2 * (lead + i) + 1],
            })
            .collect()
    }

    /// The player takes the departing stick out of the system.
    fn catch(&mut self, inst: &Instrument) -> Result<(), Refusal> {
        if self.stick.is_none() {
            return Ok(());
        }
        let velocity = self.system.state()[1];
        let heads = self.head_states();
        let system = compile_system(inst, self.rate_hz, &heads, None, &self.gate)?;
        // Mass-normalized coordinate: kinetic energy is v^2 / 2 exactly.
        self.ledger.recovered_j += 0.5 * velocity * velocity;
        self.system = system;
        self.stick = None;
        self.force = vec![0.0; self.system.mode_count()];
        self.ledger.stored_j = self.system.frame().stored_energy_j;
        Ok(())
    }

    /// Launch the stick 0.2 mm above the CURRENT membrane surface.
    fn launch(&mut self, inst: &Instrument, strike: Strike) -> Result<(), Refusal> {
        let heads = self.head_states();
        let mut membrane = 0.0;
        for (shape, state) in strike.row.iter().zip(&heads) {
            membrane += shape * state.displacement_m_sqrt_kg;
        }
        let tip = ModalAcousticState {
            displacement_m_sqrt_kg: (membrane - LAUNCH_GAP_M) * self.sqrt_mass,
            velocity_m_sqrt_kg_per_s: strike.speed_m_s * self.sqrt_mass,
        };
        let system = compile_system(inst, self.rate_hz, &heads, Some((tip, &strike.row)), &self.gate)?;
        self.ledger.supplied_j += 0.5 * tip.velocity_m_sqrt_kg_per_s * tip.velocity_m_sqrt_kg_per_s;
        self.system = system;
        self.force = vec![0.0; self.system.mode_count()];
        self.last_xy_m = strike.position_m;
        self.stick = Some(Stick { position_m: strike.position_m, row: strike.row });
        self.strikes += 1;
        self.ledger.stored_j = self.system.frame().stored_energy_j;
        Ok(())
    }

    fn mechanical_step(&mut self, inst: &Instrument) -> Result<(), Refusal> {
        if self.stick.is_some() {
            let (closure, tip_velocity) = self.closure();
            if closure <= 0.0 {
                if let Some(strike) = self.pending.take() {
                    self.catch(inst)?;
                    self.launch(inst, strike)?;
                } else if -closure >= inst.params.catch_gap_m && tip_velocity < 0.0 {
                    self.catch(inst)?;
                }
            } else {
                self.steps_in_contact += 1;
            }
        } else if let Some(strike) = self.pending.take() {
            self.launch(inst, strike)?;
        }
        let frame = self.system.step(&self.force, &self.gate).map_err(mechanics_refusal)?;
        self.steps += 1;
        if self.stick.is_some() {
            self.steps_with_stick += 1;
        }
        self.ledger.stored_j = frame.stored_energy_j;
        self.ledger.dissipated_j += frame.dissipated_energy_j;
        // No external generalized force is ever applied; kept for closure.
        self.ledger.supplied_j += frame.supplied_work_j;
        self.ledger.max_step_residual_j = self.ledger.max_step_residual_j.max(frame.balance_residual_j.abs());
        Ok(())
    }

    fn render_frame(&mut self, inst: &Instrument) -> Result<[f64; 2], Refusal> {
        for step in 0..self.substeps {
            self.mechanical_step(inst)?;
            let lead = usize::from(self.stick.is_some());
            let state = self.system.state();
            let (n0, n1) = (self.n0, self.n1);
            for i in 0..n0 {
                self.traces[0][step * n0 + i] = state[2 * (lead + i) + 1];
            }
            for i in 0..n1 {
                self.traces[1][step * n1 + i] = state[2 * (lead + n0 + i) + 1];
            }
        }
        let batter = self.mics[0].step_trace(&self.traces[0]).map_err(|e| Refusal::new("radiation_refused", e))?;
        let resonant = self.mics[1].step_trace(&self.traces[1]).map_err(|e| Refusal::new("radiation_refused", e))?;
        self.frames += 1;
        Ok([batter, resonant])
    }

    /// Render microphone pressure [Pa] for both microphones (batter side,
    /// resonant side). On a refusal the voice halts and the remaining frames
    /// are exact zeros: silence, never a held or faded sample.
    ///
    /// # Errors
    /// The refusal that halted the voice (also retained in `halted()`).
    pub fn render_pa(&mut self, inst: &Instrument, batter: &mut [f64], resonant: &mut [f64]) -> Result<(), Refusal> {
        let frames = batter.len().min(resonant.len());
        if self.halted.is_some() {
            batter.fill(0.0);
            resonant.fill(0.0);
            return Ok(());
        }
        let full_scale = inst.params.full_scale_pa;
        for frame in 0..frames {
            match self.render_frame(inst) {
                Ok([a, b]) => {
                    batter[frame] = a;
                    resonant[frame] = b;
                    self.peak_pa[0] = self.peak_pa[0].max(a.abs());
                    self.peak_pa[1] = self.peak_pa[1].max(b.abs());
                    if a.abs() > full_scale || b.abs() > full_scale {
                        self.clipped += 1;
                    }
                }
                Err(refusal) => {
                    batter[frame..].fill(0.0);
                    resonant[frame..].fill(0.0);
                    self.halted = Some(refusal.clone());
                    return Err(refusal);
                }
            }
        }
        Ok(())
    }

    /// Current ledger with the uncorrected whole-session residual.
    #[must_use]
    pub fn ledger(&self) -> Ledger {
        let mut ledger = self.ledger;
        ledger.balance_residual_j =
            ledger.stored_j - (ledger.supplied_j - ledger.recovered_j - ledger.dissipated_j);
        ledger
    }

    /// Accepted modal (displacement, velocity) of batter then resonant modes.
    #[must_use]
    pub fn modal_state(&self) -> (Vec<f64>, Vec<f64>) {
        let lead = usize::from(self.stick.is_some());
        let state = self.system.state();
        let count = self.n0 + self.n1;
        let mut q = Vec::with_capacity(count);
        let mut v = Vec::with_capacity(count);
        for i in 0..count {
            q.push(state[2 * (lead + i)]);
            v.push(state[2 * (lead + i) + 1]);
        }
        (q, v)
    }

    /// Accepted state for drawing. Resets the per-snapshot pressure peaks.
    pub fn snapshot(&mut self, inst: &Instrument) -> Snapshot {
        let (q, v) = self.modal_state();
        let state = self.system.state();
        let volume: f64 = self.areas.iter().zip(&q).map(|(a, q)| a * q).sum();
        let (closure, tip_velocity) = self.closure();
        let penetration = if self.stick.is_some() { closure.max(0.0) } else { 0.0 };
        let snapshot = Snapshot {
            time_s: self.frames as f64 / f64::from(OUTPUT_RATE_HZ),
            mechanical_steps: self.steps,
            ledger: self.ledger(),
            stick_present: self.stick.is_some(),
            stick_in_contact: penetration > 0.0,
            strike_pending: self.pending.is_some(),
            stick_xy_m: self.stick.as_ref().map_or(self.last_xy_m, |s| s.position_m),
            stick_position_m: if self.stick.is_some() { state[0] / self.sqrt_mass } else { 0.0 },
            stick_velocity_m_s: tip_velocity,
            penetration_m: penetration,
            contact_force_n: inst.stick.hertz_stiffness * penetration.powf(HERTZ_ALPHA),
            cavity_pressure_pa: self.pressure_per_volume * volume,
            strikes_launched: self.strikes,
            peak_pa: self.peak_pa,
            clipped_frames: self.clipped,
            q,
            v,
        };
        self.peak_pa = [0.0; 2];
        snapshot
    }

    /// Render `frames` audio frames with a stroke every `strike_interval`
    /// frames (0 = a single stroke at the start). The caller times this call.
    ///
    /// # Errors
    /// Any refusal while building or rendering.
    pub fn bench(
        inst: &Instrument, substeps: usize, frames: usize, strike_interval: usize,
        position_m: [f64; 2], speed_m_s: f64,
    ) -> Result<BenchReport, Refusal> {
        let mut voice = Self::new(inst, substeps)?;
        let mut peak = 0.0_f64;
        for frame in 0..frames {
            if frame == 0 || (strike_interval > 0 && frame % strike_interval == 0) {
                voice.strike(inst, position_m, speed_m_s)?;
            }
            let [batter, _] = voice.render_frame(inst)?;
            peak = peak.max(batter.abs());
        }
        Ok(BenchReport {
            frames,
            mechanical_steps: voice.steps,
            steps_with_stick: voice.steps_with_stick,
            steps_in_contact: voice.steps_in_contact,
            strikes: voice.strikes,
            peak_pa: peak,
            ledger: voice.ledger(),
        })
    }
}
