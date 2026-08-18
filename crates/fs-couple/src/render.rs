//! Block render API (music bead `frankensim-music-v8-root-3ez8g.2.1`).
//!
//! Every performance image eventually runs in a callback-shaped world:
//! AUDIO-RATE BLOCKS advance composed voices; CONTROL-RATE deltas apply
//! BETWEEN blocks (doctrine D17: state lifts at control boundaries so
//! charts move without clicks); no allocation inside a block on the
//! admitted no-alloc voices; cancellation is polled at block boundaries
//! (request -> complete the current block -> stop; blocks are
//! transactional, so a cancelled render is always a whole number of
//! blocks and resumes bitwise-identically).
//!
//! This module HOSTS the existing steppers — the exact-FIR characteristic
//! line, the scalar reed islands, the exact-ZOH modal runtime — it never
//! adds a fourth integration scheme (consolidation ruling: exactly three
//! steppers exist by design). Determinism contract: same construction
//! inputs + same control schedule => bitwise-identical pascals regardless
//! of block-size choice, because a block boundary performs no arithmetic
//! (the per-sample operation sequence is invariant; only the loop
//! partition moves).
//!
//! Allocation honesty (measured, not aspirational): the massless-reed
//! voice with an empty plate bank is allocation-free per block and is
//! gated by the counting-allocator test. The massive-reed lay path calls
//! `fs_dcontact::Obstacle::dissipative_modal_forces`, which returns a
//! `Vec` per sample, and the pHS plate path allocates inside
//! `fs_phs::step` — both are DISCLOSED allocating voices (fusion
//! candidates for bead 3ez8g.15), not silently admitted to the no-alloc
//! set.

use crate::acoustic_realize::AcousticRealizeError;
use crate::driving_point::characteristic_line;
use crate::modal_acoustic_time::{ModalAcousticTimeError, ModalAcousticTimeModel};
use crate::reed_bore::{blowing_envelope, reed_structural, solve_reed_wave};
use crate::thin_plate::PlateBank;
use fs_dcontact::Obstacle;
use fs_duct::{Duct, Termination};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_scenario::BeatingReed;

/// One reed-on-a-characteristic-line voice: the `realize_reed_bore`
/// physics, restructured so a caller can advance it block by block. The
/// one-shot realizer is now a thin wrapper over this type, so the two
/// paths cannot drift.
pub struct ReedBoreVoice {
    line: fs_vfit::discretize::DelayedFilter,
    plates: PlateBank,
    reed: BeatingReed,
    lay: Option<Obstacle>,
    rho: f64,
    sound_speed: f64,
    zc: f64,
    area_bore: f64,
    listener_m: f64,
    dt: f64,
    reed_y: f64,
    reed_v: f64,
    p_plus_prev: f64,
    f_jet_prev: f64,
    sample_index: u64,
}

impl ReedBoreVoice {
    /// Build the voice exactly as `realize_reed_bore` builds its state
    /// (same admission checks, same line realization, same lay matching,
    /// same `p_plus = 5.0` priming push), then hand stepping to the
    /// caller.
    ///
    /// # Errors
    /// Domain, TMM, or solver refusals — identical to the one-shot path.
    #[allow(clippy::too_many_arguments)] // mirrors the realizer signature
    pub fn new(
        physics: &Duct,
        gas: &GasState,
        reed: BeatingReed,
        termination: Termination,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
    ) -> Result<Self, AcousticRealizeError> {
        if !(reed.rest_opening_m > 0.0
            && reed.width_m > 0.0
            && reed.closing_pressure_pa > 0.0
            && reed.blowing_pressure_pa >= 0.0
            && reed.attack_s >= 0.0
            && reed.mass_kg >= 0.0
            && reed.stiffness_n_m >= 0.0
            && reed.rest_opening_m.is_finite())
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "reed parameters must be physical and finite",
            });
        }
        let inlet_r = physics
            .segments
            .first()
            .ok_or(AcousticRealizeError::InvalidDescription {
                what: "duct has no segments",
            })?
            .outlet_radius();
        let area_bore = core::f64::consts::PI * inlet_r * inlet_r;
        let zc = gas.density * gas.sound_speed / area_bore;
        let mut line = characteristic_line(
            physics,
            gas,
            termination,
            sample_rate_hz,
            line_samples,
            zc,
            wall,
        )
        .map_err(crate::reed_bore::map_drive)?;
        let lay = if reed.mass_kg > 0.0 {
            let k_lay = 1.0e7 * reed.width_m;
            let (_k, r_damp) = reed_structural(reed);
            let chi = r_damp / (k_lay * reed.rest_opening_m * reed.rest_opening_m).max(1.0e-18);
            Some(
                crate::unilateral_contact::slit_lay(k_lay, 2.0)
                    .and_then(|o| o.with_internal_loss(chi))
                    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
            )
        } else {
            None
        };
        let p_plus_prev = 5.0;
        let _ = line
            .push(p_plus_prev)
            .map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?;
        Ok(Self {
            line,
            plates,
            reed,
            lay,
            rho: gas.density,
            sound_speed: gas.sound_speed,
            zc,
            area_bore,
            listener_m,
            dt: 1.0 / f64::from(sample_rate_hz),
            reed_y: reed.rest_opening_m,
            reed_v: 0.0,
            p_plus_prev,
            f_jet_prev: 0.0,
            sample_index: 0,
        })
    }

    /// Total samples rendered since construction.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.sample_index
    }

    /// Mutable access to the plate bank (the one-shot realizer swaps the
    /// caller's bank in and out through this seam so its error semantics
    /// stay byte-identical to the pre-render-API code).
    pub const fn plate_bank_mut(&mut self) -> &mut PlateBank {
        &mut self.plates
    }

    /// Replace the reed's held blowing pressure. Applied between blocks
    /// only (the render context enforces the boundary); the change is a
    /// pure input-parameter move — no stored state depends on it, so the
    /// D17 lift map is empty and the transition is click-bounded by the
    /// reed's own dynamics.
    pub const fn set_blowing_pressure(&mut self, pressure_pa: f64) {
        self.reed.blowing_pressure_pa = pressure_pa;
    }

    /// Advance exactly `out.len()` samples, writing observer pascals.
    ///
    /// The per-sample sequence is byte-for-byte the one-shot realizer's
    /// loop body; a block boundary performs no arithmetic, which is the
    /// whole block-size-invariance contract.
    ///
    /// # Errors
    /// Solver or finite-set refusals; on refusal the voice must be
    /// discarded (mid-sample state is not rewound).
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), AcousticRealizeError> {
        for slot in out.iter_mut() {
            let t = self.sample_index as f64 * self.dt;
            let p_m = blowing_envelope(self.reed, t);
            let p_minus = self.line.incoming();
            let u_body = self.plates.volume_velocity();
            let p_plus = if self.reed.mass_kg > 0.0 {
                let (pp, y, v) = crate::reed_bore::step_massive_reed(
                    self.reed,
                    self.rho,
                    self.zc,
                    p_minus,
                    p_m,
                    self.reed_y,
                    self.reed_v,
                    self.dt,
                    u_body,
                    self.lay.as_ref(),
                )?;
                self.reed_y = y;
                self.reed_v = v;
                pp
            } else {
                solve_reed_wave(
                    self.reed,
                    self.rho,
                    self.zc,
                    0.0,
                    p_minus,
                    p_m,
                    self.p_plus_prev,
                    u_body,
                )?
            };
            self.p_plus_prev = p_plus;
            let p_minus_now = self
                .line
                .push(p_plus)
                .map_err(|_| AcousticRealizeError::Reed {
                    what: "characteristic line left the finite set",
                })?;
            let p_bore = p_plus + p_minus_now;
            let mut p_obs = p_bore;
            p_obs += self.plates.drive_and_radiate(
                p_bore * self.area_bore,
                self.dt,
                self.rho,
                self.listener_m,
            )?;
            let opening = if self.reed.mass_kg > 0.0 {
                self.reed_y.max(0.0)
            } else {
                crate::reed_bore::aperture_of(self.reed)
                    .opening_m(p_m - p_bore)
                    .max(0.0)
            };
            let f_jet = (p_m - p_bore) * self.reed.width_m * opening;
            let dfdt = (f_jet - self.f_jet_prev) / self.dt;
            self.f_jet_prev = f_jet;
            p_obs += self.rho * dfdt
                / (4.0 * core::f64::consts::PI * self.listener_m * self.sound_speed);
            *slot = p_obs;
            self.sample_index += 1;
        }
        Ok(())
    }
}

/// A modal string voice hosting the exact-ZOH runtime: per-sample free
/// decay (or a held generalized force) with the model's own observer
/// pressure. DISCLOSED allocating: `ModalAcousticTimeModel::step` builds a
/// per-sample energy frame `Vec`, so this voice is not in the no-alloc
/// set until a frame-free stepping seam lands (fusion bead 3ez8g.15).
pub struct ModalStringVoice {
    model: ModalAcousticTimeModel,
    held_force: Vec<f64>,
    sample_index: u64,
}

impl ModalStringVoice {
    /// Wrap an admitted model with a held (zero for free decay)
    /// generalized-force vector.
    ///
    /// # Errors
    /// Force-length mismatch against the model's mode count.
    pub fn new(
        model: ModalAcousticTimeModel,
        held_force: Vec<f64>,
    ) -> Result<Self, ModalAcousticTimeError> {
        if held_force.len() != model.modes().len() {
            return Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: model.modes().len(),
                found: held_force.len(),
            });
        }
        Ok(Self {
            model,
            held_force,
            sample_index: 0,
        })
    }

    /// Total samples rendered since construction.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.sample_index
    }

    /// Advance exactly `out.len()` samples of observer pressure.
    ///
    /// # Errors
    /// Model refusals (budget ceilings, non-finite states).
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), ModalAcousticTimeError> {
        for slot in out.iter_mut() {
            let frame = self.model.step(&self.held_force)?;
            *slot = frame.observer_pressure_pa;
            self.sample_index += 1;
        }
        Ok(())
    }
}

/// Typed control deltas applied BETWEEN blocks (control rate). v1 carries
/// the minimal set the current voices support; gesture schedules
/// (bead 3ez8g.2.3) lower into these, and image-swap lifts land with
/// their track beads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ControlDelta {
    /// Replace a reed voice's held blowing pressure [Pa].
    SetBlowingPressure {
        /// Voice index in the context.
        voice: usize,
        /// New held pressure [Pa].
        pressure_pa: f64,
    },
}

/// One applied-control record: what changed at which block boundary, and
/// what state lift it required (empty = pure input-parameter move). This
/// is the D17 lift map made inspectable.
#[derive(Clone, Debug, PartialEq)]
pub struct ControlRecord {
    /// Block index the delta applied BEFORE.
    pub block_index: u64,
    /// The delta, verbatim.
    pub delta: ControlDelta,
    /// Human-readable lift description (empty when no state moved).
    pub lift: &'static str,
}

/// Typed refusal from the render context.
#[derive(Debug)]
pub enum RenderError {
    /// A voice refused mid-block; the context is poisoned.
    Voice(AcousticRealizeError),
    /// The modal voice refused mid-block; the context is poisoned.
    Modal(ModalAcousticTimeError),
    /// A control named a voice that does not exist or cannot accept it.
    Control {
        /// What was wrong.
        what: &'static str,
    },
    /// `block` was called with an empty output slice.
    EmptyBlock,
}

impl core::fmt::Display for RenderError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Voice(e) => write!(f, "voice refusal: {e:?}"),
            Self::Modal(e) => write!(f, "modal voice refusal: {e:?}"),
            Self::Control { what } => write!(f, "control refusal: {what}"),
            Self::EmptyBlock => write!(f, "block output slice is empty"),
        }
    }
}

/// One voice slot: the admitted performance images this context can host.
#[allow(clippy::large_enum_variant)] // a handful of voice slots; boxing would ripple the API
pub enum RenderVoice {
    /// Reed on a characteristic line (wind-reed filling).
    ReedBore(ReedBoreVoice),
    /// Exact-ZOH modal string (string filling).
    ModalString(ModalStringVoice),
}

/// The block render context: a set of voices summed into one observer
/// channel, advanced a block at a time, with controls applied only at
/// block boundaries and cancellation polled only between blocks.
pub struct RenderContext {
    voices: Vec<RenderVoice>,
    /// Scratch buffer reused across blocks so summation allocates nothing
    /// after construction.
    scratch: Vec<f64>,
    blocks_rendered: u64,
    controls_applied: Vec<ControlRecord>,
}

impl RenderContext {
    /// Build a context over admitted voices, pre-sizing the scratch buffer
    /// to `max_block` samples (larger blocks refuse rather than allocate).
    #[must_use]
    pub fn new(voices: Vec<RenderVoice>, max_block: usize) -> Self {
        Self {
            voices,
            scratch: vec![0.0; max_block],
            blocks_rendered: 0,
            controls_applied: Vec::new(),
        }
    }

    /// Blocks rendered so far.
    #[must_use]
    pub const fn blocks_rendered(&self) -> u64 {
        self.blocks_rendered
    }

    /// The applied-control log (the inspectable D17 record).
    #[must_use]
    pub fn control_log(&self) -> &[ControlRecord] {
        &self.controls_applied
    }

    /// Apply control deltas at the CURRENT block boundary (before the next
    /// `block` call). Refusals leave every voice untouched.
    ///
    /// # Errors
    /// Unknown voice index or a delta the voice kind cannot accept.
    pub fn apply_controls(&mut self, deltas: &[ControlDelta]) -> Result<(), RenderError> {
        // Validate everything first: control application is transactional.
        for delta in deltas {
            match delta {
                ControlDelta::SetBlowingPressure { voice, pressure_pa } => {
                    if !pressure_pa.is_finite() || *pressure_pa < 0.0 {
                        return Err(RenderError::Control {
                            what: "blowing pressure must be finite and non-negative",
                        });
                    }
                    match self.voices.get(*voice) {
                        Some(RenderVoice::ReedBore(_)) => {}
                        Some(RenderVoice::ModalString(_)) => {
                            return Err(RenderError::Control {
                                what: "a modal string voice has no blowing pressure",
                            });
                        }
                        None => {
                            return Err(RenderError::Control {
                                what: "control names a voice index that does not exist",
                            });
                        }
                    }
                }
            }
        }
        for delta in deltas {
            match delta {
                ControlDelta::SetBlowingPressure { voice, pressure_pa } => {
                    if let Some(RenderVoice::ReedBore(reed)) = self.voices.get_mut(*voice) {
                        reed.set_blowing_pressure(*pressure_pa);
                    }
                    self.controls_applied.push(ControlRecord {
                        block_index: self.blocks_rendered,
                        delta: *delta,
                        lift: "", // pure input-parameter move; no state lifted
                    });
                }
            }
        }
        Ok(())
    }

    /// Render exactly one block: zero the output, advance every voice
    /// `out.len()` samples, sum. No allocation occurs for blocks within
    /// the pre-sized maximum and no-alloc voices; larger blocks refuse.
    ///
    /// # Errors
    /// Voice refusals poison the context (mid-sample state is not
    /// rewound); an oversized or empty block refuses before any state
    /// moves.
    pub fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if out.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if out.len() > self.scratch.len() {
            return Err(RenderError::Control {
                what: "block exceeds the pre-sized maximum; grow max_block at construction",
            });
        }
        out.fill(0.0);
        for voice in &mut self.voices {
            let scratch = &mut self.scratch[..out.len()];
            match voice {
                RenderVoice::ReedBore(reed) => {
                    reed.step_block(scratch).map_err(RenderError::Voice)?;
                }
                RenderVoice::ModalString(string) => {
                    string.step_block(scratch).map_err(RenderError::Modal)?;
                }
            }
            for (accumulator, sample) in out.iter_mut().zip(scratch.iter()) {
                *accumulator += *sample;
            }
        }
        self.blocks_rendered += 1;
        Ok(())
    }
}

/// Outcome of a gated multi-block render.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GatedRenderOutcome {
    /// Every requested block rendered.
    Completed {
        /// Blocks rendered.
        blocks: u64,
    },
    /// Cancellation was observed at a block boundary; the render stopped
    /// cleanly after a whole number of blocks (request -> drain the
    /// current block -> stop).
    Cancelled {
        /// Whole blocks rendered before stopping.
        blocks: u64,
    },
}

/// Render `blocks` blocks of `block_len` samples into `out` under a
/// cancellation gate, polling ONLY at block boundaries: an in-flight
/// block always completes (drain semantics), so the rendered prefix is a
/// whole number of blocks and a resumed context continues
/// bitwise-identically from the boundary.
///
/// # Errors
/// Voice refusals, oversized blocks, or an output slice smaller than the
/// requested render.
pub fn render_under_gate(
    context: &mut RenderContext,
    gate: &CancelGate,
    out: &mut [f64],
    block_len: usize,
    blocks: usize,
) -> Result<GatedRenderOutcome, RenderError> {
    if block_len == 0 || out.len() < block_len * blocks {
        return Err(RenderError::Control {
            what: "output slice must hold block_len * blocks samples",
        });
    }
    for index in 0..blocks {
        if gate.is_requested() {
            return Ok(GatedRenderOutcome::Cancelled {
                blocks: index as u64,
            });
        }
        let start = index * block_len;
        context.block(&mut out[start..start + block_len])?;
    }
    Ok(GatedRenderOutcome::Completed {
        blocks: blocks as u64,
    })
}
