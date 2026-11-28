//! One-way baffled-pressure Rayleigh microphone observer for modal surfaces (`.11.9.5`).
//! Reference: `fs_couple::pcm_wav::baffled::BaffledPressure::new_multirate`.
//! No radiation backreaction and no shell/room diffraction.

/// Streaming causal Blackman-Harris cascade decimator.
pub mod decimate;

use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

use fs_material::gas::GasState;
use fs_math::{dd::Dd, det};
use fs_modal::time::ModalAcousticState;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftEnvironment, MicrophonePattern, ProbeDecl, ProbeKind,
};

use self::decimate::{Decimator, MAX_MODES};
use crate::lowering::LoweredSystem;
use crate::{NodeId, System};

const MAX_PATCHES: usize = 120_000;
const MAX_HISTORY: usize = 32_768;

/// One free-node planar surface quadrature element (`z = 0`).
#[derive(Debug, Clone, PartialEq)]
pub struct RadiatingPatch {
    /// Planar patch coordinates `[x, y]` [m].
    pub xy_m: [f64; 2],
    /// Lumped P1 nodal area [m^2].
    pub area_m2: f64,
    /// Outward-signed mode shapes `outward_sign * phi_k(x, y)` [1/sqrt(kg)].
    pub shapes: Vec<f64>,
}

/// Lowered radiating modal head metadata attached to [`LoweredSystem`].
#[derive(Debug, Clone, PartialEq)]
pub struct RadiatingModalSurface {
    /// Kernel `ModalStructure` node owning modal velocities.
    pub node: NodeId,
    /// Number of retained modes.
    pub mode_count: usize,
    /// Natural angular frequencies `omega_k` [rad/s] of retained modes.
    pub mode_omegas_rad_s: Vec<f64>,
    /// Outward normal sign (`+1.0` or `-1.0`, from touching `PistonFace :sign`).
    pub outward_sign: f64,
    /// Free-node surface patches (fixed rim omitted).
    pub samples: Vec<RadiatingPatch>,
}

/// Zero-allocation streaming Rayleigh baffled-pressure microphone observer.
#[derive(Clone, Debug)]
pub struct AcousticMicrophone {
    node: NodeId,
    modes: usize,
    output_rate: f64,
    mechanics_dt_s: f64,
    kernels: Vec<(usize, Vec<f64>)>,
    history: Vec<f64>,
    history_frames: usize,
    head: usize,
    decimator: Decimator,
    sub_vel: Vec<f64>,
    sub_idx: usize,
    prev_vel: Vec<f64>,
    accel: Vec<f64>,
    last_pa: f64,
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
fn build_kernels(
    s: &RadiatingModalSurface,
    pos: [f64; 3],
    rate: u32,
    gas: &GasState,
) -> Result<Vec<(usize, Vec<f64>)>, String> {
    let mut bins: BTreeMap<usize, Vec<Dd>> = BTreeMap::new();
    for p in &s.samples {
        if !(p.area_m2.is_finite() && p.area_m2 > 0.0) || p.shapes.len() != s.mode_count {
            return Err("invalid radiating surface patch".into());
        }
        let r =
            det::sqrt((pos[0] - p.xy_m[0]).powi(2) + (pos[1] - p.xy_m[1]).powi(2) + pos[2].powi(2));
        let delay = r * f64::from(rate) / gas.sound_speed;
        if !(delay.is_finite() && delay >= 1.0 && delay < (MAX_HISTORY - 2) as f64) {
            return Err("microphone propagation delay exceeds causal history budget".into());
        }
        let (first, frac, gain) = (
            delay.floor() as usize,
            delay - delay.floor(),
            gas.density * p.area_m2 / (2.0 * PI * r),
        );
        for (d, w) in [(first, 1.0 - frac), (first + 1, frac)] {
            let sc = Dd::from_f64(gain) * Dd::from_f64(w);
            let row = bins
                .entry(d)
                .or_insert_with(|| vec![Dd::ZERO; s.mode_count]);
            for (acc, &sh) in row.iter_mut().zip(&p.shapes) {
                *acc = *acc + sc * Dd::from_f64(sh);
            }
        }
    }
    bins.into_iter()
        .map(|(d, row)| {
            row.into_iter()
                .map(|v| (v.hi.is_finite() && v.lo.is_finite()).then(|| v.to_f64()))
                .collect::<Option<Vec<_>>>()
                .map(|ws| (d, ws))
                .ok_or_else(|| "Rayleigh radiation kernel overflow".into())
        })
        .collect()
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn check_clocks(
    decl: &ProbeDecl,
    s: &RadiatingModalSurface,
    dt_s: f64,
) -> Result<(u32, usize, f64), String> {
    let (sr, r) = (decl.sample_rate.value, decl.sample_rate.value.round());
    if !(sr.is_finite() && (8_000.0..=192_000.0).contains(&r) && (sr - r).abs() <= 1e-9 * r) {
        return Err(format!("microphone sample-rate {sr} Hz outside 8..192 kHz"));
    }
    let out_rate = r as u32;
    if !(dt_s.is_finite() && dt_s > 0.0) {
        return Err("positive finite mechanics dt required".into());
    }
    let rf = 1.0 / (dt_s * f64::from(out_rate));
    let rr = rf.round();
    if !(rf.is_finite() && (1.0..=16.0).contains(&rr) && (rf - rr).abs() <= 1e-9 * rr) {
        return Err(format!(
            "decimation ratio {rf:.6} must be an integer in 1..=16"
        ));
    }
    let (ratio, mech_rate) = (rr as usize, f64::from(out_rate) * rr);
    let max_hz = s.mode_omegas_rad_s.iter().copied().fold(0.0_f64, f64::max) / TAU;
    if max_hz / mech_rate >= 0.45 {
        return Err(format!(
            "max mode {max_hz:.2} Hz crosses 0.45 * {mech_rate:.2} Hz Nyquist guard"
        ));
    }
    Ok((out_rate, ratio, mech_rate))
}

impl AcousticMicrophone {
    /// Construct a cold baffled-pressure microphone from a probe declaration and radiating surface.
    /// # Errors
    /// Pattern, `z <= 0` / `< 0.05 m`, clock/decimation, Nyquist, moist-air, or causal-delay refusal.
    pub fn from_probe(
        decl: &ProbeDecl,
        env: &DraftEnvironment,
        s: &RadiatingModalSurface,
        dt_s: f64,
    ) -> Result<(Self, GasState, usize), String> {
        let ProbeKind::Microphone {
            position, pattern, ..
        } = &decl.kind
        else {
            return Err("expected microphone probe".into());
        };
        if *pattern != MicrophonePattern::Omnidirectional {
            return Err(format!(
                "only omnidirectional microphone is lowered; got {pattern:?}"
            ));
        }
        let pos = [position[0].value, position[1].value, position[2].value];
        if pos.iter().any(|v| !v.is_finite()) {
            return Err("microphone position must be finite".into());
        }
        if pos[2] <= 0.0 {
            return Err("microphone position z must be > 0 along target outward normal".into());
        }
        if pos[2] < 0.05 {
            return Err("microphone position z must be >= 0.05 m above baffled plane".into());
        }
        let (out_rate, ratio, mech_rate) = check_clocks(decl, s, dt_s)?;
        let m = s.mode_count;
        if !(1..=MAX_MODES).contains(&m) || s.samples.is_empty() || s.samples.len() > MAX_PATCHES {
            return Err("invalid radiating surface mode or patch count".into());
        }
        let gas = GasState::try_new_moist_air(
            env.base.ambient_temperature.value,
            env.base.ambient_pressure.value,
            env.relative_humidity,
        )
        .map_err(|e| format!("exterior air state refused: {e}"))?;
        let kernels = build_kernels(s, pos, out_rate, &gas)?;
        let history_frames = kernels.last().ok_or("empty kernel")?.0 + 1;
        let n_hist = history_frames.checked_mul(m).ok_or("history overflow")?;
        let decimator = Decimator::new(ratio, m).map_err(str::to_owned)?;
        let (node, output_rate, mechanics_dt_s) = (s.node, f64::from(out_rate), 1.0 / mech_rate);
        let (history, head, sub_vel, sub_idx) = (vec![0.0; n_hist], 0, vec![0.0; ratio * m], 0);
        let (prev_vel, accel, last_pa) = (vec![0.0; m], vec![0.0; ratio * m], 0.0);
        Ok((
            Self {
                node,
                modes: m,
                output_rate,
                mechanics_dt_s,
                kernels,
                history,
                history_frames,
                head,
                decimator,
                sub_vel,
                sub_idx,
                prev_vel,
                accel,
                last_pa,
            },
            gas,
            ratio,
        ))
    }

    /// Build a microphone for a declared probe from a lowered system and step size `dt_s`.
    /// # Errors
    /// Missing/non-radiating target or [`Self::from_probe`] admission failure.
    pub fn from_draft_and_lowered(
        decl: &ProbeDecl,
        draft: &DraftAssembly,
        lowered: &LoweredSystem,
        dt_s: f64,
    ) -> Result<Self, String> {
        let ProbeKind::Microphone { target, .. } = &decl.kind else {
            return Err("expected microphone probe".into());
        };
        let surf = lowered.radiating_surfaces.get(target).ok_or_else(|| {
            format!(
                "microphone `{}` target `{}` is not a radiating modal surface",
                decl.name,
                target.token()
            )
        })?;
        Ok(Self::from_probe(decl, &draft.input().environment, surf, dt_s)?.0)
    }

    /// Advance by one accepted mechanical step with zero heap allocations.
    /// # Errors
    /// Missing modal state, non-finite velocity, or filter/pressure overflow.
    pub fn advance(&mut self, system: &System) -> Result<Option<f64>, &'static str> {
        let st = system
            .state::<Vec<ModalAcousticState>>(self.node)
            .ok_or("microphone target state missing")?;
        if st.len() != self.modes {
            return Err("microphone target modal state length mismatch");
        }
        let base = self.sub_idx * self.modes;
        for (dst, m) in self.sub_vel[base..base + self.modes].iter_mut().zip(st) {
            if !m.velocity_m_sqrt_kg_per_s.is_finite() {
                return Err("non-finite modal velocity");
            }
            *dst = m.velocity_m_sqrt_kg_per_s;
        }
        self.sub_idx += 1;
        if self.sub_idx < self.decimator.ratio() {
            return Ok(None);
        }
        self.sub_idx = 0;
        let p = self.eval_block()?;
        self.last_pa = p;
        Ok(Some(p))
    }

    /// Directly advance one complete `decimation_ratio() * modes` velocity block.
    /// # Errors
    /// Length mismatch, non-finite velocity, or filter/pressure overflow.
    pub fn step_velocities(&mut self, vel: &[f64]) -> Result<f64, &'static str> {
        if vel.len() != self.sub_vel.len() || vel.iter().any(|v| !v.is_finite()) {
            return Err("microphone needs complete finite substep velocity trace");
        }
        self.sub_vel.copy_from_slice(vel);
        let p = self.eval_block()?;
        self.last_pa = p;
        Ok(p)
    }

    #[allow(clippy::cast_precision_loss)]
    fn eval_block(&mut self) -> Result<f64, &'static str> {
        let (ratio, m) = (self.decimator.ratio(), self.modes);
        let mech_rate = self.output_rate * ratio as f64;
        for s in 0..ratio {
            for i in 0..m {
                let prev = if s == 0 {
                    self.prev_vel[i]
                } else {
                    self.sub_vel[(s - 1) * m + i]
                };
                self.accel[s * m + i] = (self.sub_vel[s * m + i] - prev) * mech_rate;
            }
        }
        self.decimator.preview(&self.accel)?;
        let mut p = 0.0;
        for (delay, ws) in &self.kernels {
            let fr = if self.head >= *delay {
                self.head - *delay
            } else {
                self.head + self.history_frames - *delay
            };
            let off = fr * m;
            p += ws
                .iter()
                .zip(&self.history[off..off + m])
                .map(|(w, a)| w * a)
                .sum::<f64>();
        }
        if !p.is_finite() {
            return Err("radiated pressure overflow");
        }
        let off = self.head * m;
        self.history[off..off + m].copy_from_slice(self.decimator.pending());
        self.prev_vel
            .copy_from_slice(&self.sub_vel[(ratio - 1) * m..]);
        self.decimator.commit();
        self.head = (self.head + 1) % self.history_frames;
        Ok(p)
    }

    /// Most recently emitted pressure [Pa].
    #[must_use]
    pub const fn last_pressure_pa(&self) -> f64 {
        self.last_pa
    }
    /// Mechanics substeps per emitted audio sample (`1..=16`).
    #[must_use]
    pub const fn decimation_ratio(&self) -> usize {
        self.decimator.ratio()
    }
    /// Admitted mechanical substep interval [s].
    #[must_use]
    pub const fn mechanics_dt_s(&self) -> f64 {
        self.mechanics_dt_s
    }
}
