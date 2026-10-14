//! The reduced real-time piano kernel.
//!
//! SAME equations and SAME time discretization as the offline reference bank
//! (`grand_piano/linear.rs`) and engine (`grand_piano/engine.rs`):
//!
//! * every string is a stiff-string modal bank whose partials come from
//!   (T, EI, mu, L); with y(x) = sum(phi_n q_n) + x/L * b the completed-square
//!   coordinate z_n = q_n + beta_n b has unit mass and couples RECIPROCALLY to
//!   the bridge through the potential sum(w_n^2 (z_n - beta_n b)^2)/2;
//! * every diagonal oscillator (string partial, loaded board mode) advances
//!   with fs-couple's exact zero-order-hold map, the cross potential uses
//!   average-displacement forces and a small board Schur solve closes its
//!   work exactly; intrinsic string loss and board loss are the reference's
//!   symmetric velocity flows around that conservative step;
//! * the hammer is the reference's geometry-derived shank with the per-key
//!   WoolFelt card and Prony creep, solved implicitly against the string/board
//!   contact compliance.
//!
//! What is REDUCED (and declared, never hidden):
//!
//! * only the SOUNDING keys' strings are dynamic. Every other string of the
//!   88-key scale stays attached to the bridge quasi-statically (its chord
//!   mass mu L/3 and tension stiffness T/L load the board), the reference's
//!   own image of a string with no retained partials. A key becomes dynamic
//!   at note-on with zero relative motion (energy-continuous) and returns to
//!   the static image once its relative energy is below a declared threshold;
//!   the energy discarded there is reported in the ledger.
//! * string partials are capped per key (mode budget) and by a partial band;
//! * one mechanical step per audio sample (the reference default is four);
//! * dampers are a K-point viscous pad (the reference uses a finer quadrature
//!   of the same estimated pad);
//! * the soundboard keeps the complete FE modal slice below `board_band_hz`.
//!
//! All dense eigen-solves, exact-ZOH coefficient probes and projections are
//! cold (construction, note-on of a new key, voice retirement). The per-sample
//! path allocates nothing.
#![allow(clippy::too_many_lines, clippy::needless_range_loop)]

use crate::refmods::{
    board_geometry, felt,
    geometry::{self, Course},
    linear::{Bank, BoardMode},
    relaxation, shank, steinway_d, steinway_scale,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample};
use fs_material::{Uniaxial, WoolFelt};
use fs_math::{c64::C64, det};
use std::f64::consts::{PI, TAU};

/// Idealized backcheck/rest gap, as in the reference engine.
const CATCH_DISTANCE: f64 = 0.020;
/// Reference launch clearance below the string (engine.rs `arm`).
const LAUNCH_CLEARANCE: f64 = 0.002;
/// Estimated damper pad of the reference (`dampers::Specification::estimated`):
/// centre at 35 % of the speaking length, width min(40 mm, 8 % L), 0.4 N s/m
/// total drag per string, no dampers above MIDI 88.
const DAMPER_CENTRE: f64 = 0.35;
const DAMPER_WIDTH_FRACTION: f64 = 0.08;
const DAMPER_WIDTH_MAX_M: f64 = 0.04;
const DAMPER_DRAG_NS_M: f64 = 0.4;
const LAST_DAMPED_MIDI: u8 = 88;
/// The scale dumped by the offline reference binary (`grand_piano --preset
/// steinway-d --dump-scale`): RT-0425 Appendix A wrapped-string model, 84
/// published courses plus four estimated extensions, tensions adjusted to
/// A4 = 440 Hz first-partial equal temperament.
pub const SCALE_CSV: &str = include_str!("../data/model-d-scale-et.csv");

/// Declared reductions of the live model. Every field is reported by `info`.
#[derive(Clone, Debug)]
pub struct Config {
    /// Audio rate and mechanical rate (one mechanical step per sample).
    pub sample_rate: u32,
    /// Complete FE soundboard modal slice (0, band] Hz is retained.
    pub board_band_hz: f64,
    /// String partials above this frequency are not retained.
    pub partial_band_hz: f64,
    /// Table ceiling of partials per string segment.
    pub max_partials: usize,
    /// Dynamic string modes granted to one sounding key (all unison strings
    /// and duplex segments together).
    pub course_mode_budget: usize,
    /// Sounding keys simulated at once. Beyond it the quietest is retired.
    pub max_courses: usize,
    /// Dynamic string modes simulated at once (the CPU budget). Beyond it the
    /// quietest key is retired; the event is counted.
    pub max_dynamic_modes: usize,
    /// Quadrature points of the estimated viscous damper pad.
    pub damper_points: usize,
    /// Simulate the duplex (aliquot) segments of sounding keys.
    pub duplex: bool,
    /// Keep every non-sounding string attached quasi-statically.
    pub static_attachment: bool,
    /// Rayleigh half-space receiver position [m] in the board frame.
    pub microphone_m: [f64; 3],
    /// Fixed pressure full scale of the float output [Pa].
    pub full_scale_pa: f64,
    /// A sounding, undamped key returns to the static image below this
    /// relative string energy [J].
    pub retire_ringing_j: f64,
    /// Same threshold for a key whose dampers are engaged [J].
    pub retire_damped_j: f64,
}
impl Config {
    pub fn live(sample_rate: u32) -> Self {
        Self {
            sample_rate,
            board_band_hz: 800.0,
            partial_band_hz: 0.45 * f64::from(sample_rate),
            max_partials: 256,
            course_mode_budget: 200,
            max_courses: 24,
            max_dynamic_modes: 2400,
            damper_points: 4,
            duplex: true,
            static_attachment: true,
            microphone_m: [0.675, 1.0, 1.0],
            full_scale_pa: 10.0,
            retire_ringing_j: 1.0e-10,
            retire_damped_j: 1.0e-8,
        }
    }
}

/// The instrument description: real scale, real board, real surface.
pub struct Description {
    pub courses: Vec<Course>,
    pub board: Vec<BoardMode>,
    /// Full viscous matrix [1/s] in the board coordinates (row-major), when the
    /// board is a bridge-informed reduction; `None` = per-mode damping ratios.
    pub damping: Option<Vec<f64>>,
    pub surface: Vec<board_geometry::SurfaceSample>,
    pub board_provenance: String,
    pub board_mass_kg: f64,
    pub board_area_m2: f64,
    /// Leading coordinates that are exact FE eigenmodes (all of them for an
    /// unreduced board; the protected prefix of a bridge-informed reduction).
    pub exact_low_modes: usize,
}

/// The soundboard shipped with the live engine, generated NATIVELY by
/// `examples/gen_board.rs` with the reference's own `steinway_d` +
/// `board_geometry::prepare_reduced` (the offline binary's
/// `--board-band-hz 2400 --board-reduction 32,16,200,400,800,1600,2300`):
/// the certified FE slice of the Model D board below 2400 Hz projected onto
/// 32 coordinates, 16 exact low modes plus 16 bridge-response Ritz vectors,
/// with its full damping matrix. Surface mode shapes are stored in single
/// precision. Solving that 2400 Hz slice inside the browser would take far
/// too long at start-up; nothing else about the board is precomputed.
pub const BOARD_BIN: &[u8] = include_bytes!("../data/model-d-board-2400-r32.bin");
const BOARD_MAGIC: &[u8; 8] = b"FSLPBRD1";

struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}
impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8], String> {
        let end = self.at.checked_add(n).filter(|e| *e <= self.bytes.len()).ok_or("board data is truncated")?;
        let slice = &self.bytes[self.at..end];
        self.at = end;
        Ok(slice)
    }
    fn u32(&mut self) -> Result<usize, String> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
    }
    fn f64(&mut self) -> Result<f64, String> {
        let b = self.take(8)?;
        let value = f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
        if value.is_finite() { Ok(value) } else { Err("board data contains a nonfinite value".into()) }
    }
    fn f32(&mut self) -> Result<f64, String> {
        let b = self.take(4)?;
        let value = f32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        if value.is_finite() { Ok(f64::from(value)) } else { Err("board data contains a nonfinite value".into()) }
    }
    fn text(&mut self) -> Result<String, String> {
        let n = self.u32()?;
        String::from_utf8(self.take(n)?.to_vec()).map_err(|_| "board provenance is not UTF-8".into())
    }
}

impl Description {
    /// The embedded scale and the embedded 32-coordinate board.
    pub fn embedded() -> Result<Self, String> {
        let courses = geometry::read_scale(SCALE_CSV)?;
        if courses.len() != 88 {
            return Err("embedded Model D scale must contain 88 courses".into());
        }
        Self::from_board_bytes(courses, BOARD_BIN)
    }

    /// Serialize the board of this description (`examples/gen_board.rs`).
    /// Layout (little-endian): magic, r, keys, surface points, exact low
    /// modes, mass, area, provenance; per coordinate frequency, damping ratio,
    /// volume; bridge rows [key][coordinate] for MIDI 21..108; damping matrix
    /// (zeros when the board is unreduced); per surface point x, y, area (f64)
    /// then its r mode shapes (f32).
    pub fn board_bytes(&self) -> Vec<u8> {
        let r = self.board.len();
        let mut out = Vec::new();
        out.extend_from_slice(BOARD_MAGIC);
        for value in [r, 88, self.surface.len(), self.exact_low_modes, usize::from(self.damping.is_some())] {
            out.extend_from_slice(&(value as u32).to_le_bytes());
        }
        out.extend_from_slice(&self.board_mass_kg.to_le_bytes());
        out.extend_from_slice(&self.board_area_m2.to_le_bytes());
        out.extend_from_slice(&(self.board_provenance.len() as u32).to_le_bytes());
        out.extend_from_slice(self.board_provenance.as_bytes());
        for b in &self.board {
            for value in [b.frequency_hz, b.damping_ratio, b.volume] {
                out.extend_from_slice(&value.to_le_bytes());
            }
        }
        for key in 0..88 {
            for b in &self.board {
                out.extend_from_slice(&b.bridge[key].to_le_bytes());
            }
        }
        if let Some(c) = &self.damping {
            for value in c {
                out.extend_from_slice(&value.to_le_bytes());
            }
        }
        for p in &self.surface {
            for value in [p.position_m[0], p.position_m[1], p.area_m2] {
                out.extend_from_slice(&value.to_le_bytes());
            }
            for value in &p.mode_shape {
                out.extend_from_slice(&(*value as f32).to_le_bytes());
            }
        }
        out
    }

    /// Read a board written by `board_bytes`. Strict: no value is repaired.
    pub fn from_board_bytes(courses: Vec<Course>, bytes: &[u8]) -> Result<Self, String> {
        let mut input = Reader { bytes, at: 0 };
        if input.take(8)? != BOARD_MAGIC {
            return Err("not a live-piano board file".into());
        }
        let (r, keys, points, exact, has_damping) = (input.u32()?, input.u32()?, input.u32()?, input.u32()?, input.u32()?);
        if !(1..=128).contains(&r) || keys != 88 || points == 0 || points > 120_000 || exact > r || has_damping > 1 {
            return Err("board file header is out of range".into());
        }
        let (mass, area) = (input.f64()?, input.f64()?);
        let provenance = input.text()?;
        let mut board = Vec::with_capacity(r);
        for _ in 0..r {
            let (frequency_hz, damping_ratio, volume) = (input.f64()?, input.f64()?, input.f64()?);
            board.push(BoardMode { frequency_hz, damping_ratio, bridge: [0.0; 88], volume });
        }
        for key in 0..88 {
            for b in board.iter_mut() {
                b.bridge[key] = input.f64()?;
            }
        }
        let damping = if has_damping == 1 {
            let mut c = Vec::with_capacity(r * r);
            for _ in 0..r * r {
                c.push(input.f64()?);
            }
            Some(c)
        } else {
            None
        };
        let mut surface = Vec::with_capacity(points);
        for _ in 0..points {
            let (x, y, area_m2) = (input.f64()?, input.f64()?, input.f64()?);
            let mut mode_shape = Vec::with_capacity(r);
            for _ in 0..r {
                mode_shape.push(input.f32()?);
            }
            surface.push(board_geometry::SurfaceSample { position_m: [x, y, 0.0], area_m2, mode_shape });
        }
        if input.at != bytes.len() {
            return Err("board file has trailing bytes".into());
        }
        Ok(Self { courses, board, damping, surface, board_provenance: provenance, board_mass_kg: mass, board_area_m2: area, exact_low_modes: exact })
    }

    /// Keep only the first `count` coordinates. Legal only inside the exact
    /// low-mode prefix: the result is then the complete FE slice below the
    /// next mode, with the per-mode damping of those exact modes.
    pub fn truncated(mut self, count: usize) -> Result<Self, String> {
        let r = self.board.len();
        if count == 0 || count > self.exact_low_modes {
            return Err(format!("only the {} exact low modes may be kept on their own", self.exact_low_modes));
        }
        self.board.truncate(count);
        for p in &mut self.surface {
            p.mode_shape.truncate(count);
        }
        if let Some(c) = self.damping.take() {
            let mut kept = Vec::with_capacity(count * count);
            for i in 0..count {
                kept.extend_from_slice(&c[i * r..i * r + count]);
            }
            self.damping = Some(kept);
        }
        self.exact_low_modes = count;
        self.board_provenance.push_str(&format!("; truncated to its {count} exact low modes"));
        Ok(self)
    }

    /// The reference's `--preset steinway-d`: its dumped 88-key scale and the
    /// Model D soundboard assembled and eigen-solved by the reference's own
    /// `steinway_d` + `board_geometry` modules (default mesh divisions 8).
    pub fn model_d(board_band_hz: f64) -> Result<Self, String> {
        let courses = geometry::read_scale(SCALE_CSV)?;
        if courses.len() != 88 {
            return Err("embedded Model D scale must contain 88 courses".into());
        }
        Self::with_courses(courses, board_band_hz)
    }
    /// The same board with an explicit (possibly sparse) scale.
    pub fn with_courses(courses: Vec<Course>, board_band_hz: f64) -> Result<Self, String> {
        let preset = steinway_d::build(8)?;
        let keys: Vec<u8> = courses.iter().map(|c| c.midi).collect();
        let prepared = board_geometry::BoardGeometry::read(&preset.geometry)?.prepare(&keys, board_band_hz)?;
        Ok(Self {
            courses,
            board: prepared.modes,
            damping: prepared.physical_damping,
            surface: prepared.surface,
            board_provenance: prepared.provenance,
            board_mass_kg: prepared.mass_kg,
            board_area_m2: prepared.area_m2,
            exact_low_modes: usize::MAX,
        }.with_exact_count())
    }
    fn with_exact_count(mut self) -> Self {
        self.exact_low_modes = self.exact_low_modes.min(self.board.len());
        self
    }
    /// The same board through the reference's explicit bridge-informed
    /// reduction (`--board-reduction max,keep_low,Hz,...` of the offline
    /// binary): a certified FE slice below `source_band_hz` is projected onto
    /// `max` coordinates, `keep_low` exact low modes plus static and harmonic
    /// bridge responses at the listed frequencies, with its full damping matrix.
    pub fn with_reduction(courses: Vec<Course>, source_band_hz: f64, spec: &str) -> Result<Self, String> {
        let options = board_geometry::ritz::RitzOptions::parse(spec)?;
        let preset = steinway_d::build(8)?;
        let keys: Vec<u8> = courses.iter().map(|c| c.midi).collect();
        let prepared = board_geometry::BoardGeometry::read(&preset.geometry)?
            .prepare_reduced(&keys, source_band_hz, false, false, false, false, &options)?;
        if prepared.physical_damping.is_none() {
            return Err("reduced soundboard requires its full projected material damping".into());
        }
        let exact = prepared.reduction.as_ref().map_or(0, |report| report.protected_low_modes);
        Ok(Self {
            exact_low_modes: exact,
            courses,
            board: prepared.modes,
            damping: prepared.physical_damping,
            surface: prepared.surface,
            board_provenance: format!("{}; bridge-informed reduction {spec} of the (0,{source_band_hz}] Hz FE slice", prepared.provenance),
            board_mass_kg: prepared.mass_kg,
            board_area_m2: prepared.area_m2,
        })
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Zoh {
    qq: f64,
    qv: f64,
    vq: f64,
    vv: f64,
    bq: f64,
    bv: f64,
}

/// Cold lowering of fs-couple's public exact-ZOH oscillator to a linear map,
/// exactly as the reference bank does (`linear.rs::transition_chunk`).
fn zoh_probe(rate: u32, omega: &[f64]) -> Result<Vec<Zoh>, String> {
    let mut out = Vec::with_capacity(omega.len());
    for chunk in omega.chunks(4096) {
        let modes: Vec<ModalAcousticMode> = chunk
            .iter()
            .map(|&w| ModalAcousticMode {
                angular_frequency_rad_s: w,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::new(0.0, 0.0),
            })
            .collect();
        let model = ModalAcousticTimeModel::try_new(rate, modes, ModalAcousticTimeBudget::audible_reference())
            .map_err(|e| e.to_string())?;
        for i in 0..chunk.len() {
            let (qq, qv, vq, vv, bq, bv) = model.probe_linear_map(i);
            if !bq.is_finite() || bq <= 0.0 {
                return Err("nonpositive ZOH displacement compliance".into());
            }
            out.push(Zoh { qq, qv, vq, vv, bq, bv });
        }
    }
    Ok(out)
}

fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 1.0;
    }
    a
}

/// Four-lane dot product: lanes `s[0..4]` accumulate the whole 4-blocks, are
/// combined as `(s0 + s1) + (s2 + s3)`, then the tails (`len % 4`) are added in
/// order. Both loop shapes below produce the same blocks, tails and operation
/// order, so the sum is bit-identical on every target; the shape is per target
/// because LLVM vectorizes them differently (measured, optimisation pass 1):
/// - wasm32: `as_chunks` + a single-exit zip; the lanes SLP-vectorize to f64x2.
///   Measured ~5% less CPU for the whole engine in V8 than the two-exit loop.
/// - elsewhere: the two-exit `(&mut ca).zip(&mut cb)` loop. Do NOT "fix" it to a
///   single-exit (trusted-random-access) loop: AArch64's loop vectorizer then
///   applies an ordered FP reduction (lane extracts, serial scalar adds, spills),
///   measured 1.5-1.65x slower per dot on Apple silicon. The `split_at_checked`
///   frame a profiler attributes to this loop is sampling skid from the four
///   add-latency chains, which bound the loop (the sum order fixes them).
#[cfg(target_arch = "wasm32")]
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let (ca, ra) = a.as_chunks::<4>();
    let (cb, rb) = b.as_chunks::<4>();
    let mut s = [0.0_f64; 4];
    for (x, y) in ca.iter().zip(cb) {
        s[0] += x[0] * y[0];
        s[1] += x[1] * y[1];
        s[2] += x[2] * y[2];
        s[3] += x[3] * y[3];
    }
    let mut t = (s[0] + s[1]) + (s[2] + s[3]);
    for (x, y) in ra.iter().zip(rb) {
        t += x * y;
    }
    t
}

/// See the wasm32 `dot` above: same sum, two-exit loop shape for native targets.
///
/// `inline(always)`, not `inline` (optimisation pass 6): with plain `inline`
/// LLVM simplifies `dot` on its own first and inlines it afterwards. At call
/// sites where both slices provably have the same length (`matvec` rows
/// against `x`, `gl . mid`, `cross . v` over one segment) it then folds the
/// two exits into one, and the loop vectorizer emits the slow ordered
/// reduction anyway (ld4.2d, lane extracts, serial scalar adds). Inlined
/// before simplification, the loop keeps the two-exit fadd.2d/fmul.2d form.
/// Same operations in the same order, so the sums are bit-identical.
#[cfg(not(target_arch = "wasm32"))]
#[inline(always)]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = [0.0_f64; 4];
    let mut ca = a.chunks_exact(4);
    let mut cb = b.chunks_exact(4);
    for (x, y) in (&mut ca).zip(&mut cb) {
        s[0] += x[0] * y[0];
        s[1] += x[1] * y[1];
        s[2] += x[2] * y[2];
        s[3] += x[3] * y[3];
    }
    let mut t = (s[0] + s[1]) + (s[2] + s[3]);
    for (x, y) in ca.remainder().iter().zip(cb.remainder()) {
        t += x * y;
    }
    t
}

/// Pad shapes of one segment packed for the fused damper passes (1 to 4
/// pads; otherwise empty and `damp` runs the sequential flows): for each
/// 4-mode block `c`, `g_0[4c..4c+4], ..., g_{np-1}[4c..4c+4]`, then each
/// pad's tail modes (`n % 4`) pad-major. Also the Gram matrix `g_t . g_j`.
fn pack_pads(pad: &[PadPoint], n: usize) -> (Vec<f64>, Vec<f64>) {
    let np = pad.len();
    if !(1..=4).contains(&np) {
        return (Vec::new(), Vec::new());
    }
    let blocks = n / 4;
    let mut pack = Vec::with_capacity(n * np);
    for c in 0..blocks {
        for p in pad {
            pack.extend_from_slice(&p.shape[4 * c..4 * c + 4]);
        }
    }
    for p in pad {
        pack.extend_from_slice(&p.shape[4 * blocks..n]);
    }
    let mut gram = vec![0.0; np * np];
    for t in 0..np {
        for j in 0..np {
            gram[t * np + j] = dot(&pad[t].shape, &pad[j].shape);
        }
    }
    (pack, gram)
}

/// `[g_t . v; N]` for the `N` packed pad shapes in ONE pass over `v`; each
/// sum has `dot`'s four-lane order. Loop shape per target as in `dot`.
#[cfg(target_arch = "wasm32")]
#[inline]
fn pad_dots<const N: usize>(pack: &[f64], v: &[f64]) -> [f64; N] {
    let mut s = [[0.0_f64; 4]; N];
    let (cv, rv) = v.as_chunks::<4>();
    let (cg, _) = pack.as_chunks::<4>();
    for (x, g) in cv.iter().zip(cg.chunks_exact(N)) {
        for t in 0..N {
            s[t][0] += g[t][0] * x[0];
            s[t][1] += g[t][1] * x[1];
            s[t][2] += g[t][2] * x[2];
            s[t][3] += g[t][3] * x[3];
        }
    }
    pad_dots_finish(&s, &pack[cv.len() * 4 * N..], rv)
}

#[cfg(not(target_arch = "wasm32"))]
#[inline]
fn pad_dots<const N: usize>(pack: &[f64], v: &[f64]) -> [f64; N] {
    let mut s = [[0.0_f64; 4]; N];
    let mut cv = v.chunks_exact(4);
    let mut cg = pack.chunks_exact(4 * N);
    for (x, g) in (&mut cv).zip(&mut cg) {
        for t in 0..N {
            s[t][0] += g[4 * t] * x[0];
            s[t][1] += g[4 * t + 1] * x[1];
            s[t][2] += g[4 * t + 2] * x[2];
            s[t][3] += g[4 * t + 3] * x[3];
        }
    }
    let rv = cv.remainder();
    pad_dots_finish(&s, &pack[(v.len() - rv.len()) * N..], rv)
}

#[inline]
fn pad_dots_finish<const N: usize>(s: &[[f64; 4]; N], tail: &[f64], rv: &[f64]) -> [f64; N] {
    let rem = rv.len();
    let mut d = [0.0; N];
    for t in 0..N {
        let mut acc = (s[t][0] + s[t][1]) + (s[t][2] + s[t][3]);
        for (g, x) in tail[t * rem..(t + 1) * rem].iter().zip(rv) {
            acc += g * x;
        }
        d[t] = acc;
    }
    d
}

/// `v += sum_t g_t * c[t]` for the `N` packed pad shapes in ONE pass (each
/// element adds the pads in pack order).
#[inline]
fn pad_axpy<const N: usize>(pack: &[f64], v: &mut [f64], c: &[f64; N]) {
    let body = v.len() / 4 * 4;
    let (head, rv) = v.split_at_mut(body);
    for (x, g) in head.chunks_exact_mut(4).zip(pack.chunks_exact(4 * N)) {
        for t in 0..N {
            x[0] += g[4 * t] * c[t];
            x[1] += g[4 * t + 1] * c[t];
            x[2] += g[4 * t + 2] * c[t];
            x[3] += g[4 * t + 3] * c[t];
        }
    }
    let rem = rv.len();
    let tail = &pack[body * N..];
    for t in 0..N {
        for (x, g) in rv.iter_mut().zip(&tail[t * rem..(t + 1) * rem]) {
            *x += g * c[t];
        }
    }
}

/// The `N` sequential exact rank-one pad flows of one segment (`damp`), fused:
/// with `v0` the entry velocities, pad `t` (in processing order) sees
/// `g_t . v = g_t . v0 + sum_{j before t} (g_t . g_j) c_j` and the bridge
/// speed `beta = gl . bvel`, which each flow moves by `|gl|^2 lift_t c_t`.
/// One multi-dot pass, an N x N recurrence, one multi-axpy pass. `beta` and
/// `push` (= sum lift_t c_t, the bridge velocity change still owed to
/// `bvel` along `gl`) carry across the course's segments. Returns the loss.
#[inline]
fn damp_fused<const N: usize>(s: &ASeg, v: &mut [f64], gl_norm2: f64, beta: &mut f64, push: &mut f64, reverse: bool) -> f64 {
    let d = pad_dots::<N>(&s.pad_pack, v);
    let mut c = [0.0_f64; N];
    let mut loss = 0.0;
    let mut moved = false;
    for k in 0..N {
        let t = if reverse { N - 1 - k } else { k };
        let p = &s.pad[t];
        if p.coef == 0.0 {
            continue;
        }
        let row = &s.pad_gram[t * N..t * N + N];
        let mut string_speed = d[t];
        for kk in 0..k {
            let j = if reverse { N - 1 - kk } else { kk };
            string_speed += row[j] * c[j];
        }
        let speed = string_speed + p.lift * *beta;
        let change = p.coef * speed;
        c[t] = change;
        let lc = p.lift * change;
        *beta += gl_norm2 * lc;
        *push += lc;
        loss += -change * speed - 0.5 * change * change * p.norm;
        moved = true;
    }
    if moved {
        pad_axpy::<N>(&s.pad_pack, v, &c);
    }
    loss
}

/// `out = row^T * a` for a row-major `r x r` matrix `a`.
fn matvec_into_row(row: &[f64], a: &[f64], r: usize, out: &mut [f64]) {
    out.fill(0.0);
    for (b, weight) in row.iter().enumerate() {
        for (o, x) in out.iter_mut().zip(&a[b * r..(b + 1) * r]) {
            *o += weight * x;
        }
    }
}

#[inline]
fn matvec(a: &[f64], x: &[f64], y: &mut [f64]) {
    let r = x.len();
    for (row, out) in a.chunks_exact(r).zip(y.iter_mut()) {
        *out = dot(row, x);
    }
}

#[inline]
fn step_string_loss(
    v: &mut [f64],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    bm: f64,
) -> (f64, f64) {
    let mut l = [0.0_f64; 4];
    let mut a = [0.0_f64; 4];
    let mut cv = v.chunks_exact_mut(4);
    let mut cc1 = c1.chunks_exact(4);
    let mut cc2 = c2.chunks_exact(4);
    let mut cc3 = c3.chunks_exact(4);
    for (((vc, c1c), c2c), c3c) in (&mut cv).zip(&mut cc1).zip(&mut cc2).zip(&mut cc3) {
        for lane in 0..4 {
            let old = vc[lane];
            let dv = c1c[lane] * old + c2c[lane] * bm;
            l[lane] += (old - dv) * (old + dv);
            a[lane] += c3c[lane] * old;
            vc[lane] = dv;
        }
    }
    let mut loss = (l[0] + l[1]) + (l[2] + l[3]);
    let mut acc3 = (a[0] + a[1]) + (a[2] + a[3]);
    for (((v, &c1), &c2), &c3) in cv
        .into_remainder()
        .iter_mut()
        .zip(cc1.remainder())
        .zip(cc2.remainder())
        .zip(cc3.remainder())
    {
        let old = *v;
        let dv = c1 * old + c2 * bm;
        loss += (old - dv) * (old + dv);
        acc3 += c3 * old;
        *v = dv;
    }
    (loss, acc3)
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn step_zoh_free(
    q: &[f64],
    v: &mut [f64],
    fq: &mut [f64],
    qq: &[f64],
    qv: &[f64],
    vq: &[f64],
    vv: &[f64],
    ha: &[f64],
) -> f64 {
    let mut r = [0.0_f64; 4];
    let mut cq = q.chunks_exact(4);
    let mut cv = v.chunks_exact_mut(4);
    let mut cfq = fq.chunks_exact_mut(4);
    let mut cqq = qq.chunks_exact(4);
    let mut cqv = qv.chunks_exact(4);
    let mut cvq = vq.chunks_exact(4);
    let mut cvv = vv.chunks_exact(4);
    let mut cha = ha.chunks_exact(4);
    for (((((((qc, vc), fqc), qqc), qvc), vqc), vvc), hac) in (&mut cq)
        .zip(&mut cv)
        .zip(&mut cfq)
        .zip(&mut cqq)
        .zip(&mut cqv)
        .zip(&mut cvq)
        .zip(&mut cvv)
        .zip(&mut cha)
    {
        for lane in 0..4 {
            let qi = qc[lane];
            let dv = vc[lane];
            let x = qqc[lane] * qi + qvc[lane] * dv;
            fqc[lane] = x;
            vc[lane] = vqc[lane] * qi + vvc[lane] * dv;
            r[lane] += hac[lane] * (qi + x);
        }
    }
    let mut reac = (r[0] + r[1]) + (r[2] + r[3]);
    for (((((((&q, v), fq), &qq), &qv), &vq), &vv), &ha) in cq
        .remainder()
        .iter()
        .zip(cv.into_remainder().iter_mut())
        .zip(cfq.into_remainder().iter_mut())
        .zip(cqq.remainder())
        .zip(cqv.remainder())
        .zip(cvq.remainder())
        .zip(cvv.remainder())
        .zip(cha.remainder())
    {
        let dv = *v;
        let x = qq * q + qv * dv;
        *fq = x;
        *v = vq * q + vv * dv;
        reac += ha * (q + x);
    }
    reac
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn step_zoh_forced(
    q: &mut [f64],
    v: &mut [f64],
    fq: &[f64],
    bqa: &[f64],
    bva: &[f64],
    ha: &[f64],
    cross: &[f64],
    bbar: f64,
) -> (f64, f64) {
    let mut r = [0.0_f64; 4];
    let mut s = [0.0_f64; 4];
    let mut cq = q.chunks_exact_mut(4);
    let mut cv = v.chunks_exact_mut(4);
    let mut cfq = fq.chunks_exact(4);
    let mut cbqa = bqa.chunks_exact(4);
    let mut cbva = bva.chunks_exact(4);
    let mut cha = ha.chunks_exact(4);
    let mut ccross = cross.chunks_exact(4);
    for ((((((qc, vc), fqc), bqac), bvac), hac), crossc) in (&mut cq)
        .zip(&mut cv)
        .zip(&mut cfq)
        .zip(&mut cbqa)
        .zip(&mut cbva)
        .zip(&mut cha)
        .zip(&mut ccross)
    {
        for lane in 0..4 {
            let nq = fqc[lane] + bqac[lane] * bbar;
            let nv = vc[lane] + bvac[lane] * bbar;
            r[lane] += hac[lane] * (qc[lane] + nq);
            s[lane] += crossc[lane] * nv;
            qc[lane] = nq;
            vc[lane] = nv;
        }
    }
    let mut reac = (r[0] + r[1]) + (r[2] + r[3]);
    let mut src = (s[0] + s[1]) + (s[2] + s[3]);
    for ((((((q, v), &fq), &bqa), &bva), &ha), &cross) in cq
        .into_remainder()
        .iter_mut()
        .zip(cv.into_remainder().iter_mut())
        .zip(cfq.remainder())
        .zip(cbqa.remainder())
        .zip(cbva.remainder())
        .zip(cha.remainder())
        .zip(ccross.remainder())
    {
        let nq = fq + bqa * bbar;
        let nv = *v + bva * bbar;
        reac += ha * (*q + nq);
        src += cross * nv;
        *q = nq;
        *v = nv;
    }
    (reac, src)
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn step_zoh_forced_contact(
    q: &mut [f64],
    v: &mut [f64],
    fq: &[f64],
    bqa: &[f64],
    bva: &[f64],
    ha: &[f64],
    cross: &[f64],
    bqhs: &[f64],
    bvhs: &[f64],
    bbar: f64,
    force: f64,
) -> (f64, f64) {
    let mut r = [0.0_f64; 4];
    let mut s = [0.0_f64; 4];
    let mut cq = q.chunks_exact_mut(4);
    let mut cv = v.chunks_exact_mut(4);
    let mut cfq = fq.chunks_exact(4);
    let mut cbqa = bqa.chunks_exact(4);
    let mut cbva = bva.chunks_exact(4);
    let mut cha = ha.chunks_exact(4);
    let mut ccross = cross.chunks_exact(4);
    let mut cbqhs = bqhs.chunks_exact(4);
    let mut cbvhs = bvhs.chunks_exact(4);
    for ((((((((qc, vc), fqc), bqac), bvac), hac), crossc), bqhsc), bvhsc) in (&mut cq)
        .zip(&mut cv)
        .zip(&mut cfq)
        .zip(&mut cbqa)
        .zip(&mut cbva)
        .zip(&mut cha)
        .zip(&mut ccross)
        .zip(&mut cbqhs)
        .zip(&mut cbvhs)
    {
        for lane in 0..4 {
            let nq = fqc[lane] + bqac[lane] * bbar + bqhsc[lane] * force;
            let nv = vc[lane] + bvac[lane] * bbar + bvhsc[lane] * force;
            r[lane] += hac[lane] * (qc[lane] + nq);
            s[lane] += crossc[lane] * nv;
            qc[lane] = nq;
            vc[lane] = nv;
        }
    }
    let mut reac = (r[0] + r[1]) + (r[2] + r[3]);
    let mut src = (s[0] + s[1]) + (s[2] + s[3]);
    for ((((((((q, v), &fq), &bqa), &bva), &ha), &cross), &bqhs), &bvhs) in cq
        .into_remainder()
        .iter_mut()
        .zip(cv.into_remainder().iter_mut())
        .zip(cfq.remainder())
        .zip(cbqa.remainder())
        .zip(cbva.remainder())
        .zip(cha.remainder())
        .zip(ccross.remainder())
        .zip(cbqhs.remainder())
        .zip(cbvhs.remainder())
    {
        let nq = fq + bqa * bbar + bqhs * force;
        let nv = *v + bva * bbar + bvhs * force;
        reac += ha * (*q + nq);
        src += cross * nv;
        *q = nq;
        *v = nv;
    }
    (reac, src)
}

/// One string segment of the scale (speaking unison member or its duplex).
#[derive(Clone, Debug)]
struct Seg {
    course: usize,
    member: usize,
    duplex: bool,
    off: usize,
    n_max: usize,
    /// Chord (static) endpoint mass mu*L/3 [kg].
    endpoint_mass: f64,
    /// Chord stiffness T/L [N/m].
    t_over_l: f64,
    /// Distributed velocity-loss rate 2R [1/s] of this segment's partials.
    mass_rate: f64,
    length_m: f64,
    modal_mass_kg: f64,
    strike_fraction: f64,
}

/// Everything the pre half-flow loop of `step` reads per mode, packed into one
/// 64-byte record (one cache line, one prefetch stream) instead of eight
/// separate arrays. Values are copies of the `Modes` constants of the same
/// name; the ZOH free-motion map `fq = qq*q + qv*v`, `fv = vq*q + vv*v` has no
/// other reader, so it lives only here.
#[derive(Clone, Copy)]
#[repr(C, align(64))]
struct PreFlow {
    c1: f64,
    c2: f64,
    c3: f64,
    qq: f64,
    qv: f64,
    vq: f64,
    vv: f64,
    ha: f64,
}

/// What the finish loop of `step` reads per mode besides the state: the board
/// drive `bqa = bq*a`, `bva = bv*a`, the reaction weight `ha` and the loss
/// cross term `cross`, in one 32-byte record (one stream instead of four).
/// `ha` and `cross` are copies of the `Modes` arrays of the same name (other
/// readers keep those); `bqa`/`bva` live only here.
#[derive(Clone, Copy)]
#[repr(C, align(32))]
struct FinFlow {
    bqa: f64,
    bva: f64,
    ha: f64,
    cross: f64,
}

/// The felt-force drive `bqhs = bq*hs`, `bvhs = bv*hs` the finish loop adds
/// while a segment carries a contact force (16 bytes; read only there).
#[derive(Clone, Copy)]
#[repr(C, align(16))]
struct HitFlow {
    bqhs: f64,
    bvhs: f64,
}

/// Per-mode constants (structure of arrays, bank order).
#[derive(Default)]
struct Modes {
    omega: Vec<f64>,
    w2: Vec<f64>,
    beta: Vec<f64>,
    hs: Vec<f64>,
    /// Pre half-flow record per mode (see `PreFlow`).
    pre: Vec<PreFlow>,
    bq: Vec<f64>,
    /// 0.5*a, 0.5*bq*a with a = w^2*beta.
    ha: Vec<f64>,
    hba: Vec<f64>,
    /// Finish-loop records per mode (see `FinFlow`, `HitFlow`).
    fin: Vec<FinFlow>,
    hit: Vec<HitFlow>,
    /// String-loss half-flow: v' = c1*v + c2*b_mid, bridge reaction -c3*v.
    /// (The post half-flow reads these 3 of the 8 `PreFlow` fields alone, so
    /// it keeps the plain arrays.)
    c1: Vec<f64>,
    c2: Vec<f64>,
    c3: Vec<f64>,
    c4: Vec<f64>,
    cross: Vec<f64>,
    /// Loss Schur weight term of the mode.
    lw: Vec<f64>,
}

struct PadPoint {
    shape: Vec<f64>,
    lift: f64,
    norm_string: f64,
    /// expm1(-drag*norm*dt)/norm for the current engagement; refreshed cold.
    coef: f64,
    norm: f64,
}

struct ASeg {
    seg: usize,
    off: usize,
    n: usize,
    duplex: bool,
    member: usize,
    beta2: f64,
    mek: f64,
    w_schur: f64,
    loss_w: f64,
    sum_c4: f64,
    endpoint_drag: f64,
    hammer_lift: f64,
    lift_total: f64,
    cdiag: f64,
    pad: Vec<PadPoint>,
    /// Pad shapes packed for the fused damper passes (`pack_pads`); empty
    /// unless `1 <= pad.len() <= 4`.
    pad_pack: Vec<f64>,
    /// Pad Gram matrix `g_t . g_j`, row-major `pad.len()` squared (same rule).
    pad_gram: Vec<f64>,
    /// Held contact force of this step [N] (speaking strings only).
    force: f64,
    src: f64,
}

struct ACourse {
    ci: usize,
    segs: Vec<ASeg>,
    /// Bridge row in the current loaded board basis and S^-1 times it.
    gl: Vec<f64>,
    resp: Vec<f64>,
    gl_norm2: f64,
    modes: usize,
    b0: f64,
    damped_now: bool,
    /// Samples since activation (diagnostic).
    age: u64,
    /// Partial count and mode offset of the first speaking string for display
    /// (`0` for duplex-only keys).
    shape_n: usize,
    shape_off: usize,
}

#[derive(Clone, Debug)]
struct Contact {
    state: felt::State,
    memory: relaxation::Memory,
    overlap: f64,
    force: f64,
    /// Recoverable felt energy at the committed overlap [J].
    stored_j: f64,
}

struct Key {
    midi: u8,
    held: bool,
    slot: Option<usize>,
    model: shank::Prepared,
    rest: shank::State,
    rest_energy: f64,
    motion: shank::State,
    active: bool,
    /// Creep memory still relaxing after the hammer was caught.
    relaxing: bool,
    pending_velocity: Option<f64>,
    /// Unforced shank step of the current sample (state, bending loss).
    predicted: Option<(shank::State, f64)>,
    refused: bool,
    law: WoolFelt,
    creep: relaxation::Prepared,
    area: f64,
    thickness: f64,
    contacts: Vec<Contact>,
    has_damper: bool,
    /// Bank string indices of this course: (speaking, duplex) per member.
    segs: Vec<usize>,
    last_force: f64,
}

/// Energy ledger [J]. `residual` is what the books do not explain.
#[derive(Clone, Copy, Debug, Default)]
pub struct Ledger {
    pub input_j: f64,
    pub stored_j: f64,
    pub string_loss_j: f64,
    pub board_loss_j: f64,
    pub damper_loss_j: f64,
    pub felt_loss_j: f64,
    pub shank_loss_j: f64,
    pub catch_loss_j: f64,
    /// Relative string energy discarded when keys returned to the static image.
    pub retired_j: f64,
    pub residual_j: f64,
}
impl Ledger {
    pub fn dissipated_j(&self) -> f64 {
        self.string_loss_j + self.board_loss_j + self.damper_loss_j + self.felt_loss_j + self.shank_loss_j + self.catch_loss_j
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Counters {
    pub samples: u64,
    pub note_ons: u64,
    pub refused_note_ons: u64,
    pub rebuilds: u64,
    pub retirements: u64,
    pub forced_retirements: u64,
    pub contact_failures: u64,
    pub clipped_samples: u64,
    pub peak_pa: f64,
}

pub struct LivePiano {
    cfg: Config,
    r: usize,
    dt: f64,
    courses: Vec<Course>,
    segs: Vec<Seg>,
    m: Modes,
    q: Vec<f64>,
    v: Vec<f64>,
    fq: Vec<f64>,
    fv: Vec<f64>,
    keys: Vec<Key>,
    /// Bridge rows in the STATIC basis (every string attached as a chord),
    /// course-major. In that basis the static pencil is (diag(lam_s), I).
    g_st: Vec<f64>,
    lam_s: Vec<f64>,
    /// Physical board damping matrix in the static basis, row-major.
    c_st: Vec<f64>,
    /// Bare FE modal coordinates = phi_s * static coordinates (display only).
    phi_s: Vec<f64>,
    // Current loaded basis (cold products).
    phi: Vec<f64>,
    k_phys: Vec<f64>,
    flow: Vec<f64>,
    bz: Vec<Zoh>,
    s_inv: Vec<f64>,
    sd_inv: Vec<f64>,
    // Board state in loaded coordinates and scratch.
    bq: Vec<f64>,
    bvel: Vec<f64>,
    bfq: Vec<f64>,
    bfv: Vec<f64>,
    rhs: Vec<f64>,
    end: Vec<f64>,
    mid: Vec<f64>,
    force_b: Vec<f64>,
    tmp: Vec<f64>,
    st_v: Vec<f64>,
    active: Vec<ACourse>,
    /// Keys whose hammer is in flight/contact or whose creep still relaxes.
    busy: Vec<usize>,
    /// Strikes waiting for the next sample: (course, hammer velocity).
    pending_on: Vec<(usize, f64)>,
    sustain: f64,
    engagement_cached: f64,
    mic: BaffledPressure,
    pub ledger: Ledger,
    pub counters: Counters,
    retire_clock: u32,
    board_freq_hz: Vec<f64>,
    board_provenance: String,
    exact_low_modes: usize,
    surface_xy: Vec<f32>,
    surface_shapes: Vec<f32>,
    // contact scratch (bounded by 3 contacts per busy hammer)
    c_key: Vec<usize>,
    c_member: Vec<usize>,
    c_gap: Vec<f64>,
    c_force: Vec<f64>,
    c_h: Vec<f64>,
    c_lift: Vec<f64>,
    c_site: Vec<crate::fastfelt::Site>,
    c_start: Vec<f64>,
    /// Cold caches per course: damper pad projections and display table.
    pad_cache: Vec<Option<(usize, Vec<(Vec<f64>, f64, f64)>)>>,
    shape_cache: Vec<Option<(usize, Vec<f32>)>>,
    last_error: Option<String>,
}

const SHAPE_STATIONS: usize = 48;

impl LivePiano {
    pub fn new(desc: &Description, cfg: Config) -> Result<Self, String> {
        let rate = cfg.sample_rate;
        if !(8_000..=192_000).contains(&rate) {
            return Err("sample rate must lie in 8..192 kHz".into());
        }
        if !(1..=512).contains(&cfg.max_partials) || cfg.damper_points == 0 || cfg.damper_points > 16
            || cfg.max_courses == 0 || cfg.course_mode_budget == 0
            || !cfg.partial_band_hz.is_finite() || cfg.partial_band_hz <= 0.0
            || cfg.partial_band_hz > 0.45 * f64::from(rate)
            || !cfg.full_scale_pa.is_finite() || cfg.full_scale_pa <= 0.0 {
            return Err("invalid live piano configuration".into());
        }
        let courses = desc.courses.clone();
        let board = &desc.board;
        let r = board.len();
        let dt = 1.0 / f64::from(rate);
        // The reference bank is the single source of every per-partial constant
        // (frequency, bridge cross inertia, strike projection, intrinsic loss).
        let bank = Bank::new_with_string_damping(&courses, board, rate, 0.45 * f64::from(rate),
            cfg.max_partials, true, None, false)?;
        let (losses, _) = crate::refmods::bank_string_losses(&bank);
        let n_total = bank.modes.len();
        let omega: Vec<f64> = bank.modes.iter().map(|m| m.omega).collect();
        let zoh = zoh_probe(rate, &omega)?;
        let dur = 0.5 * dt; // each loss half-flow spans half a step
        let half = 0.5 * dur;
        let mut m = Modes::default();
        for (k, mode) in bank.modes.iter().enumerate() {
            let z = zoh[k];
            let a = mode.a;
            let (mr, sr) = (losses[k].mass_rate, losses[k].strain_rate);
            let inv = 1.0 / (1.0 + half * (mr + sr));
            let cross = half * sr * mode.beta * inv;
            m.omega.push(mode.omega);
            m.w2.push(mode.omega * mode.omega);
            m.beta.push(mode.beta);
            m.hs.push(mode.hammer_shape);
            let ha = 0.5 * a;
            let c1 = 1.0 - dur * (mr + sr) * inv;
            let c2 = -dur * ((mr + sr) * cross - sr * mode.beta);
            let c3 = sr * mode.beta * inv;
            m.pre.push(PreFlow { c1, c2, c3, qq: z.qq, qv: z.qv, vq: z.vq, vv: z.vv, ha });
            m.bq.push(z.bq);
            m.ha.push(ha);
            m.hba.push(0.5 * z.bq * a);
            m.fin.push(FinFlow { bqa: z.bq * a, bva: z.bv * a, ha, cross });
            m.hit.push(HitFlow { bqhs: z.bq * mode.hammer_shape, bvhs: z.bv * mode.hammer_shape });
            m.c1.push(c1);
            m.c2.push(c2);
            m.c3.push(c3);
            m.c4.push(sr * mode.beta * (cross - mode.beta));
            m.cross.push(cross);
            m.lw.push(half * sr * mode.beta * mode.beta * (1.0 + half * mr) * inv);
        }
        let mut segs = Vec::with_capacity(bank.strings.len());
        let mut key_segs: Vec<Vec<usize>> = vec![Vec::new(); courses.len()];
        for (si, s) in bank.strings.iter().enumerate() {
            let c = &courses[s.course];
            let cents = (s.member as f64 - 0.5 * (c.unison - 1) as f64) * c.detune_cents;
            let tension = c.tension_at_cents(cents)?;
            let length = if s.duplex { c.duplex_length_m } else { c.length_m };
            let mass_rate = s.modes.clone().next().map_or(0.0, |k| losses[k].mass_rate);
            segs.push(Seg {
                course: s.course,
                member: s.member,
                duplex: s.duplex,
                off: s.modes.start,
                n_max: s.modes.len(),
                endpoint_mass: c.linear_density_kg_m * length / 3.0,
                t_over_l: tension / length,
                mass_rate,
                length_m: length,
                modal_mass_kg: 0.5 * c.linear_density_kg_m * length,
                strike_fraction: c.strike_fraction,
            });
            key_segs[s.course].push(si);
        }
        drop(bank);
        // Bare board (mass-normalized FE modes) with EVERY string of the scale
        // attached as a chord: mass mu*L/3 and stiffness T/L on its bridge row.
        // Its eigenbasis (the fs-modal facility, once, cold) is the fixed
        // "static" basis of the live model: pencil (diag(lam_s), I).
        let mut g_bare = vec![0.0; courses.len() * r];
        for (ci, c) in courses.iter().enumerate() {
            for (j, b) in board.iter().enumerate() {
                g_bare[ci * r + j] = b.bridge[usize::from(c.midi - 21)];
            }
        }
        let mut m_stat = identity(r);
        let mut k_stat = vec![0.0; r * r];
        for (i, b) in board.iter().enumerate() {
            k_stat[i * r + i] = (TAU * b.frequency_hz).powi(2);
        }
        for seg in &segs {
            let g = &g_bare[seg.course * r..(seg.course + 1) * r];
            for i in 0..r {
                for j in 0..r {
                    m_stat[i * r + j] += seg.endpoint_mass * g[i] * g[j];
                    k_stat[i * r + j] += seg.t_over_l * g[i] * g[j];
                }
            }
        }
        let static_modes = fs_modal::eigh_gen_dense(&k_stat, &m_stat, r).map_err(|e| e.to_string())?;
        let mut phi_s = vec![0.0; r * r];
        let mut lam_s = vec![0.0; r];
        for (j, e) in static_modes.iter().enumerate() {
            if !e.lambda.is_finite() || e.lambda <= 0.0 || e.residual > 1e-7 * e.lambda {
                return Err("statically loaded soundboard mode is unresolved".into());
            }
            lam_s[j] = e.lambda;
            for i in 0..r {
                phi_s[i * r + j] = e.phi[i];
            }
        }
        let mut g_st = vec![0.0; courses.len() * r];
        for ci in 0..courses.len() {
            for j in 0..r {
                g_st[ci * r + j] = (0..r).map(|i| g_bare[ci * r + i] * phi_s[i * r + j]).sum();
            }
        }
        // Physical board damping: the supplied full matrix of a reduced board,
        // otherwise diag(2*zeta*omega) of the FE modes; projected once.
        let mut c_bare = vec![0.0; r * r];
        match &desc.damping {
            Some(c) => {
                if c.len() != r * r || c.iter().any(|x| !x.is_finite()) {
                    return Err("board damping must be a complete finite r x r matrix".into());
                }
                c_bare.copy_from_slice(c);
            }
            None => {
                for (a, b) in board.iter().enumerate() {
                    c_bare[a * r + a] = 2.0 * b.damping_ratio * TAU * b.frequency_hz;
                }
            }
        }
        let mut c_phi = vec![0.0; r * r];
        for a in 0..r {
            matvec_into_row(&c_bare[a * r..(a + 1) * r], &phi_s, r, &mut c_phi[a * r..(a + 1) * r]);
        }
        let mut c_st = vec![0.0; r * r];
        for i in 0..r {
            for j in i..r {
                let value: f64 = (0..r).map(|a| phi_s[a * r + i] * c_phi[a * r + j]).sum();
                c_st[i * r + j] = value;
                c_st[j * r + i] = value;
            }
        }
        // Hammers: the reference preset's shank geometry and per-key felt cards.
        let geometry = shank::Geometry::published();
        let mut keys = Vec::with_capacity(courses.len());
        for (ci, c) in courses.iter().enumerate() {
            let (law, prony) = steinway_scale::hammer_material(c)?;
            let area = c.felt_area_m2 / c.unison as f64;
            let thickness = c.felt_thickness_m;
            let creep = relaxation::Spectrum::from_prony(&prony)?.prepare(area, thickness, dt)?;
            let model = shank::Prepared::from_geometry(geometry, c.hammer_mass_kg, rate)?;
            let rest = model.rest(CATCH_DISTANCE);
            let rest_energy = model.energy(&rest, CATCH_DISTANCE);
            let contacts = (0..c.unison)
                .map(|_| Contact { state: law.initial_state(), memory: relaxation::Memory::default(), overlap: -CATCH_DISTANCE, force: 0.0, stored_j: 0.0 })
                .collect();
            keys.push(Key {
                midi: c.midi,
                held: false,
                slot: None,
                model,
                rest,
                rest_energy,
                motion: rest,
                active: false,
                relaxing: false,
                pending_velocity: None,
                predicted: None,
                refused: false,
                law,
                creep,
                area,
                thickness,
                contacts,
                has_damper: c.midi <= LAST_DAMPED_MIDI,
                segs: key_segs[ci].clone(),
                last_force: 0.0,
            });
        }
        // The library's time-domain Rayleigh half-space receiver, fed with the
        // static-basis board velocities; the FE surface shapes are projected
        // into that same basis once, here.
        let gas = fs_material::gas::GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)
            .map_err(|e| e.to_string())?;
        let surface: Vec<SurfaceSample> = desc
            .surface
            .iter()
            .map(|p| SurfaceSample {
                position_m: p.position_m,
                area_m2: p.area_m2,
                mode_shape: (0..r).map(|j| (0..r).map(|i| p.mode_shape[i] * phi_s[i * r + j]).sum()).collect(),
            })
            .collect();
        let mic = BaffledPressure::new_multirate_differentiating(&surface, r, rate, rate, cfg.microphone_m,
            RayleighMedium { density: gas.density, sound_speed: gas.sound_speed })?;
        let surface_xy = desc.surface.iter().flat_map(|p| [p.position_m[0] as f32, p.position_m[1] as f32]).collect();
        let surface_shapes: Vec<f32> = desc.surface.iter().flat_map(|p| p.mode_shape.iter().map(|x| *x as f32)).collect();
        let n_courses = courses.len();
        let mut piano = Self {
            r,
            dt,
            courses,
            segs,
            q: vec![0.0; n_total],
            v: vec![0.0; n_total],
            fq: vec![0.0; n_total],
            fv: vec![0.0; n_total],
            m,
            keys,
            g_st,
            lam_s,
            c_st,
            phi_s,
            phi: identity(r),
            k_phys: vec![0.0; r * r],
            flow: identity(r),
            bz: vec![Zoh::default(); r],
            s_inv: identity(r),
            sd_inv: identity(r),
            bq: vec![0.0; r],
            bvel: vec![0.0; r],
            bfq: vec![0.0; r],
            bfv: vec![0.0; r],
            rhs: vec![0.0; r],
            end: vec![0.0; r],
            mid: vec![0.0; r],
            force_b: vec![0.0; r],
            tmp: vec![0.0; r],
            st_v: vec![0.0; r],
            active: Vec::with_capacity(cfg.max_courses + 1),
            busy: Vec::with_capacity(32),
            pending_on: Vec::with_capacity(32),
            sustain: 0.0,
            engagement_cached: f64::NAN,
            mic,
            ledger: Ledger::default(),
            counters: Counters::default(),
            retire_clock: 0,
            board_freq_hz: board.iter().map(|b| b.frequency_hz).collect(),
            board_provenance: desc.board_provenance.clone(),
            exact_low_modes: desc.exact_low_modes.min(r),
            surface_xy,
            surface_shapes,
            c_key: Vec::with_capacity(96),
            c_member: Vec::with_capacity(96),
            c_gap: Vec::with_capacity(96),
            c_force: Vec::with_capacity(96),
            c_h: Vec::with_capacity(96 * 96),
            c_lift: Vec::with_capacity(96),
            c_site: Vec::with_capacity(96),
            c_start: Vec::with_capacity(96),
            pad_cache: vec![None; n_courses],
            shape_cache: vec![None; n_courses],
            last_error: None,
            cfg,
        };
        piano.rebuild(&[], &[])?;
        Ok(piano)
    }

    pub fn config(&self) -> &Config { &self.cfg }
    pub fn board_modes(&self) -> usize { self.r }
    pub fn board_frequencies_hz(&self) -> &[f64] { &self.board_freq_hz }
    pub fn board_provenance(&self) -> &str { &self.board_provenance }
    /// Leading board coordinates that are exact FE eigenmodes.
    pub fn board_exact_low_modes(&self) -> usize { self.exact_low_modes }
    pub fn surface_shapes(&self) -> &[f32] {
        &self.surface_shapes
    }
    pub fn surface_xy(&self) -> &[f32] { &self.surface_xy }
    pub fn courses(&self) -> &[Course] { &self.courses }
    pub fn sounding_keys(&self) -> usize { self.active.len() }
    pub fn dynamic_modes(&self) -> usize { self.active.iter().map(|a| a.modes).sum() }
    pub fn dynamic_strings(&self) -> usize { self.active.iter().map(|a| a.segs.iter().filter(|s| !s.duplex).count()).sum() }
    pub fn take_error(&mut self) -> Option<String> { self.last_error.take() }
    pub fn microphone_delay_samples(&self) -> (usize, f64) { (self.mic.delay_samples.0, self.mic.filter_delay_samples()) }
    pub fn microphone_multiplies(&self) -> usize { self.mic.multiply_adds_per_sample() }
    pub fn total_strings(&self) -> usize { self.segs.iter().filter(|s| !s.duplex).count() }
    pub fn table_modes(&self) -> usize { self.q.len() }

    fn key_index(&self, midi: u8) -> Option<usize> { self.courses.iter().position(|c| c.midi == midi) }

    /// Change quality limits. Applies to keys that start sounding afterwards.
    pub fn set_quality(&mut self, course_mode_budget: usize, partial_band_hz: f64, max_dynamic_modes: usize) -> Result<(), String> {
        if course_mode_budget == 0 || max_dynamic_modes == 0 || !partial_band_hz.is_finite() || partial_band_hz <= 0.0
            || partial_band_hz > 0.45 * f64::from(self.cfg.sample_rate) {
            return Err("invalid quality limits".into());
        }
        self.cfg.course_mode_budget = course_mode_budget;
        self.cfg.partial_band_hz = partial_band_hz;
        self.cfg.max_dynamic_modes = max_dynamic_modes;
        Ok(())
    }

    pub fn set_sustain(&mut self, value: f64) -> Result<(), String> {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err("sustain must lie in [0,1]".into());
        }
        self.sustain = value;
        Ok(())
    }
    pub fn sustain(&self) -> f64 { self.sustain }

    /// Partial counts (speaking, duplex) granted to a course by the budget.
    /// Both follow the partial band. When the budget (or the table ceiling)
    /// cuts the speaking strings short, the duplex segments stop at the same
    /// frequency, so a key never keeps aliquot partials above its own top.
    fn grant(&self, ci: usize) -> (usize, usize) {
        let key = &self.keys[ci];
        let speaking: Vec<&Seg> = key.segs.iter().map(|&s| &self.segs[s]).filter(|s| !s.duplex).collect();
        let duplex: Vec<&Seg> = key.segs.iter().map(|&s| &self.segs[s]).filter(|s| s.duplex).collect();
        let band = TAU * self.cfg.partial_band_hz;
        let s0 = speaking[0];
        let in_band = |seg: &Seg, limit: f64| self.m.omega[seg.off..seg.off + seg.n_max].iter().take_while(|w| **w <= limit).count();
        let n_band = speaking.iter().map(|s| in_band(s, band)).min().unwrap_or(1).max(1);
        // The table itself is band-limited unless it hit the partial ceiling.
        let table_capped = speaking.iter().any(|s| s.n_max >= self.cfg.max_partials);
        let mut n = n_band;
        loop {
            let cut = n < n_band || (table_capped && n == s0.n_max);
            let limit = if cut { self.m.omega[s0.off + n - 1].min(band) } else { band };
            let d = if self.cfg.duplex { duplex.iter().map(|s| in_band(s, limit)).min().unwrap_or(0) } else { 0 };
            let total = speaking.len() * n + duplex.len() * d;
            if total <= self.cfg.course_mode_budget || n == 1 {
                return (n, d);
            }
            n -= 1;
        }
    }

    fn make_aseg(&mut self, si: usize, n: usize) -> ASeg {
        let seg = &self.segs[si];
        let off = seg.off;
        let m = &self.m;
        let (mut beta2, mut mek, mut w_schur, mut loss_w, mut sum_c4) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut hammer_lift, mut lift_bq, mut cdiag) = (seg.strike_fraction, 0.0, 0.0);
        for k in off..off + n {
            let a = 2.0 * m.ha[k];
            beta2 += m.beta[k] * m.beta[k];
            mek += m.w2[k] * m.beta[k] * m.beta[k];
            w_schur += 0.25 * a * a * m.bq[k];
            loss_w += m.lw[k];
            sum_c4 += m.c4[k];
            hammer_lift -= m.hs[k] * m.beta[k];
            lift_bq += 0.5 * a * m.bq[k] * m.hs[k];
            cdiag += m.bq[k] * m.hs[k] * m.hs[k];
        }
        let residual_mass = seg.endpoint_mass - beta2;
        let endpoint_drag = seg.mass_rate * residual_mass;
        let half = 0.25 * self.dt;
        loss_w += half * endpoint_drag;
        let mut pad = Vec::new();
        if !seg.duplex && self.keys[seg.course].has_damper {
            let ci = seg.course;
            if self.pad_cache[ci].as_ref().is_none_or(|(count, _)| *count != n) {
                // The reference's estimated pad, as `damper_points` midpoint
                // stations: y = phi*z + (x/L - sum(phi*beta))*b at each station.
                let width = (DAMPER_WIDTH_FRACTION * seg.length_m).min(DAMPER_WIDTH_MAX_M);
                let start = DAMPER_CENTRE - 0.5 * width / seg.length_m;
                let span = width / seg.length_m;
                let count = self.cfg.damper_points;
                let root = det::sqrt(seg.modal_mass_kg);
                let mut points = Vec::with_capacity(count);
                for i in 0..count {
                    let fraction = start + span * (i as f64 + 0.5) / count as f64;
                    let shape: Vec<f64> = (1..=n).map(|p| det::sin(p as f64 * PI * fraction) / root).collect();
                    let lift = fraction - shape.iter().zip(&m.beta[off..off + n]).map(|(g, b)| g * b).sum::<f64>();
                    let norm_string = shape.iter().map(|g| g * g).sum::<f64>();
                    points.push((shape, lift, norm_string));
                }
                self.pad_cache[ci] = Some((n, points));
            }
            if let Some((_, points)) = &self.pad_cache[ci] {
                pad = points.iter().map(|(shape, lift, norm_string)| PadPoint {
                    shape: shape.clone(), lift: *lift, norm_string: *norm_string, coef: 0.0, norm: 0.0,
                }).collect();
            }
        }
        let (pad_pack, pad_gram) = pack_pads(&pad, n);
        ASeg {
            seg: si, off, n, duplex: self.segs[si].duplex, member: self.segs[si].member, beta2, mek, w_schur, loss_w, sum_c4,
            endpoint_drag, hammer_lift, lift_total: hammer_lift + lift_bq, cdiag, pad, pad_pack, pad_gram, force: 0.0, src: 0.0,
        }
    }

    /// Cold change of the dynamic string set: retire `(slot, segment)` pairs
    /// (`usize::MAX` = the whole key), give `activate` keys their missing
    /// strings, re-diagonalize the loaded board and carry the physical board
    /// motion across. Energy-continuous for activation; the relative energy of
    /// retired strings goes to `ledger.retired_j`.
    fn rebuild(&mut self, activate: &[usize], retire: &[(usize, usize)]) -> Result<(), String> {
        let r = self.r;
        // Physical board motion (static coordinates) before the basis changes.
        let mut qs = vec![0.0; r];
        let mut vs = vec![0.0; r];
        matvec(&self.phi, &self.bq, &mut qs);
        matvec(&self.phi, &self.bvel, &mut vs);
        for &(slot, which) in retire {
            let ci = self.active[slot].ci;
            let g = &self.g_st[ci * r..(ci + 1) * r];
            let b = dot(g, &qs);
            let bdot = dot(g, &vs);
            let ac = &mut self.active[slot];
            let mut kept = Vec::with_capacity(ac.segs.len());
            for s in ac.segs.drain(..) {
                if which != usize::MAX && s.seg != which {
                    kept.push(s);
                    continue;
                }
                let range = s.off..s.off + s.n;
                let mut e = 0.0;
                for k in range.clone() {
                    let x = self.q[k] - self.m.beta[k] * b;
                    e += 0.5 * self.v[k] * self.v[k] + 0.5 * self.m.w2[k] * x * x;
                }
                // Its chord keeps following the bridge: that kinetic energy stays.
                self.ledger.retired_j += e - 0.5 * s.beta2 * bdot * bdot;
                self.q[range.clone()].fill(0.0);
                self.v[range].fill(0.0);
                self.counters.retirements += 1;
            }
            ac.segs = kept;
        }
        if self.active.iter().any(|ac| ac.segs.is_empty()) {
            for ac in &self.active {
                self.keys[ac.ci].slot = None;
            }
            self.active.retain(|ac| !ac.segs.is_empty());
            for (slot, ac) in self.active.iter().enumerate() {
                self.keys[ac.ci].slot = Some(slot);
            }
        }
        for &ci in activate {
            let (n, d) = self.grant(ci);
            let (b, bdot) = {
                let g = &self.g_st[ci * r..(ci + 1) * r];
                (dot(g, &qs), dot(g, &vs))
            };
            let present: Vec<usize> = match self.keys[ci].slot {
                Some(slot) => self.active[slot].segs.iter().map(|s| s.seg).collect(),
                None => Vec::new(),
            };
            let wanted: Vec<(usize, usize)> = self.keys[ci].segs.iter()
                .map(|&si| (si, if self.segs[si].duplex { d } else { n }))
                .filter(|(si, count)| *count > 0 && !present.contains(si))
                .collect();
            let mut fresh = Vec::with_capacity(wanted.len());
            for (si, count) in wanted {
                fresh.push(self.make_aseg(si, count));
            }
            // A string joins at rest RELATIVE to its chord: q_n = 0, so
            // z_n = beta_n*b. Kinetic and potential energy are continuous.
            for s in &fresh {
                for k in s.off..s.off + s.n {
                    self.q[k] = self.m.beta[k] * b;
                    self.v[k] = self.m.beta[k] * bdot;
                }
            }
            match self.keys[ci].slot {
                Some(slot) => {
                    let ac = &mut self.active[slot];
                    ac.segs.extend(fresh);
                    ac.segs.sort_by_key(|s| s.seg);
                }
                None => {
                    self.keys[ci].slot = Some(self.active.len());
                    self.active.push(ACourse {
                        ci, segs: fresh, gl: vec![0.0; r], resp: vec![0.0; r], gl_norm2: 0.0, modes: 0, b0: 0.0,
                        damped_now: false, age: 0, shape_n: 0, shape_off: 0,
                    });
                }
            }
        }
        // Per-key bookkeeping (display shape table is populated lazily in snapshot()).
        for slot in 0..self.active.len() {
            let ac = &mut self.active[slot];
            ac.modes = ac.segs.iter().map(|s| s.n).sum();
            match ac.segs.iter().find(|s| !s.duplex).map(|s| (s.n, s.off)) {
                Some((n, off)) => {
                    ac.shape_n = n;
                    ac.shape_off = off;
                }
                None => {
                    ac.shape_n = 0;
                    ac.shape_off = 0;
                }
            }
        }
        // Loaded board pencil of the CURRENT dynamic set, in the static basis:
        // a sounding string trades its chord mass for modal mass (-beta^2) and
        // adds its modal endpoint stiffness sum(w^2 beta^2).
        let mut mass = identity(r);
        let mut k_loaded = vec![0.0; r * r];
        for i in 0..r {
            k_loaded[i * r + i] = self.lam_s[i];
        }
        for ac in &self.active {
            let g = &self.g_st[ac.ci * r..(ac.ci + 1) * r];
            let beta2: f64 = ac.segs.iter().map(|s| s.beta2).sum();
            let mek: f64 = ac.segs.iter().map(|s| s.mek).sum();
            for i in 0..r {
                for j in 0..r {
                    mass[i * r + j] -= beta2 * g[i] * g[j];
                    k_loaded[i * r + j] += mek * g[i] * g[j];
                }
            }
        }
        let (lambda, phi) = crate::dense::eigh_gen(&k_loaded, &mass, r)?;
        let band = 0.45 * f64::from(self.cfg.sample_rate);
        if lambda.iter().any(|l| !l.is_finite() || *l <= 0.0 || l.sqrt() / TAU > band) {
            return Err("loaded soundboard mode unresolved or above output band".into());
        }
        self.phi = phi;
        // New loaded coordinates of the same physical motion: Phi^T M x.
        let mut mq = vec![0.0; r];
        let mut mv = vec![0.0; r];
        matvec(&mass, &qs, &mut mq);
        matvec(&mass, &vs, &mut mv);
        for j in 0..r {
            self.bq[j] = (0..r).map(|i| self.phi[i * r + j] * mq[i]).sum();
            self.bvel[j] = (0..r).map(|i| self.phi[i * r + j] * mv[i]).sum();
        }
        // Phi^T K_end Phi for the energy books (K_end = diag(lam_s) excludes the
        // modal endpoint stiffness, which lives in the string potential), and
        // Phi^T C Phi, the physical board damping in the loaded basis.
        let mut c_phi = vec![0.0; r * r];
        for a in 0..r {
            matvec_into_row(&self.c_st[a * r..(a + 1) * r], &self.phi, r, &mut c_phi[a * r..(a + 1) * r]);
        }
        let mut c = vec![0.0; r * r];
        for i in 0..r {
            for j in i..r {
                let stiffness: f64 = (0..r).map(|a| self.phi[a * r + i] * self.lam_s[a] * self.phi[a * r + j]).sum();
                self.k_phys[i * r + j] = stiffness;
                self.k_phys[j * r + i] = stiffness;
                let damping: f64 = (0..r).map(|a| self.phi[a * r + i] * c_phi[a * r + j]).sum();
                c[i * r + j] = 0.5 * self.dt * damping;
                c[j * r + i] = 0.5 * self.dt * damping;
            }
        }
        // Exact half-step flow exp(-C dt/2) of that damping.
        self.flow = crate::dense::exp_neg_small(&c, r)?;
        let board_omega: Vec<f64> = lambda.iter().map(|l| l.sqrt()).collect();
        self.bz = zoh_probe(self.cfg.sample_rate, &board_omega)?;
        // Bridge rows in the loaded basis, coupling Schur and loss Schur.
        let mut schur = vec![0.0; r * r];
        for i in 0..r {
            schur[i * r + i] = 1.0 / self.bz[i].bq;
        }
        let mut loss_schur = identity(r);
        for ac in &mut self.active {
            let g = &self.g_st[ac.ci * r..(ac.ci + 1) * r];
            for j in 0..r {
                ac.gl[j] = (0..r).map(|i| g[i] * self.phi[i * r + j]).sum();
            }
            ac.gl_norm2 = dot(&ac.gl, &ac.gl);
            let w: f64 = ac.segs.iter().map(|s| s.w_schur).sum();
            let lw: f64 = ac.segs.iter().map(|s| s.loss_w).sum();
            for i in 0..r {
                for j in 0..r {
                    schur[i * r + j] -= w * ac.gl[i] * ac.gl[j];
                    loss_schur[i * r + j] += lw * ac.gl[i] * ac.gl[j];
                }
            }
        }
        self.s_inv = crate::dense::spd_inverse(&schur, r, "bridge Schur complement")?;
        self.sd_inv = crate::dense::spd_inverse(&loss_schur, r, "string-loss Schur complement")?;
        for ac in &mut self.active {
            matvec(&self.s_inv, &ac.gl, &mut ac.resp);
        }
        self.engagement_cached = f64::NAN; // pad coefficients depend on gl
        self.counters.rebuilds += 1;
        Ok(())
    }

    /// Relative (chord-following removed) energy of one dynamic string [J].
    fn segment_energy(&self, s: &ASeg, b: f64, bdot: f64) -> f64 {
        let mut e = 0.0;
        for k in s.off..s.off + s.n {
            let x = self.q[k] - self.m.beta[k] * b;
            let u = self.v[k] - self.m.beta[k] * bdot;
            e += 0.5 * u * u + 0.5 * self.m.w2[k] * x * x;
        }
        e
    }
    /// Relative energy of all dynamic strings of a sounding key [J].
    fn relative_energy(&self, ac: &ACourse) -> f64 {
        let b = dot(&ac.gl, &self.bq);
        let bdot = dot(&ac.gl, &self.bvel);
        ac.segs.iter().map(|s| self.segment_energy(s, b, bdot)).sum()
    }

    fn is_damped(&self, ac: &ACourse) -> bool {
        let key = &self.keys[ac.ci];
        key.has_damper && !key.held && self.sustain < 1.0
    }

    /// Post-escapement hammer velocity [m/s]. Returns false when refused.
    pub fn note_on(&mut self, midi: u8, hammer_velocity_m_s: f64) -> Result<(), String> {
        if !hammer_velocity_m_s.is_finite() || hammer_velocity_m_s <= 0.0 || hammer_velocity_m_s > 8.0 {
            self.counters.refused_note_ons += 1;
            return Err("hammer velocity must lie in (0, 8] m/s".into());
        }
        let Some(ci) = self.key_index(midi) else {
            self.counters.refused_note_ons += 1;
            return Err("key is not present in the scale".into());
        };
        self.keys[ci].held = true;
        // Strikes received between two samples start together: one cold
        // voice-set change serves the whole chord (see `apply_pending_notes`).
        self.pending_on.retain(|(key, _)| *key != ci);
        self.pending_on.push((ci, hammer_velocity_m_s));
        self.counters.note_ons += 1;
        Ok(())
    }

    /// Start every strike queued by `note_on` (called at the next sample;
    /// public so offline callers can inspect the voice set immediately).
    pub fn apply_pending_notes(&mut self) {
        if self.pending_on.is_empty() {
            return;
        }
        let pending = std::mem::take(&mut self.pending_on);
        // Which keys need strings that are not dynamic yet, and how many modes?
        let mut activate: Vec<usize> = Vec::new();
        let mut modes = self.dynamic_modes();
        let mut keys = self.active.len();
        for &(ci, _) in &pending {
            let slot = self.keys[ci].slot;
            let speaking = self.keys[ci].segs.iter().filter(|&&s| !self.segs[s].duplex).count();
            let complete = slot.is_some_and(|slot| self.active[slot].segs.iter().filter(|s| !s.duplex).count() == speaking);
            if complete {
                continue;
            }
            let (n, d) = self.grant(ci);
            let present: Vec<usize> = slot.map_or(Vec::new(), |slot| self.active[slot].segs.iter().map(|s| s.seg).collect());
            modes += self.keys[ci].segs.iter().filter(|si| !present.contains(si))
                .map(|&si| if self.segs[si].duplex { d } else { n }).sum::<usize>();
            keys += usize::from(slot.is_none());
            activate.push(ci);
        }
        if !activate.is_empty() {
            // Make room first: the quietest key that is not being struck yields
            // (keys still held down are the last to go). The event is counted.
            let mut retire: Vec<(usize, usize)> = Vec::new();
            if modes > self.cfg.max_dynamic_modes || keys > self.cfg.max_courses {
                let mut victims: Vec<(f64, usize, usize)> = self.active.iter().enumerate()
                    .filter(|(_, ac)| !pending.iter().any(|(ci, _)| *ci == ac.ci)
                        && !self.keys[ac.ci].active && self.keys[ac.ci].pending_velocity.is_none())
                    .map(|(slot, ac)| (self.relative_energy(ac) * if self.keys[ac.ci].held { 1.0e6 } else { 1.0 }, slot, ac.modes))
                    .collect();
                victims.sort_by(|a, b| a.0.total_cmp(&b.0));
                for (_, victim, victim_modes) in victims {
                    if modes <= self.cfg.max_dynamic_modes && keys <= self.cfg.max_courses {
                        break;
                    }
                    retire.push((victim, usize::MAX));
                    modes -= victim_modes;
                    keys -= 1;
                    self.counters.forced_retirements += 1;
                }
                // Still no room: the last strikes of this batch are refused.
                while (modes > self.cfg.max_dynamic_modes || keys > self.cfg.max_courses) && !activate.is_empty() {
                    let ci = activate.pop().unwrap_or(0);
                    let (n, d) = self.grant(ci);
                    let slot = self.keys[ci].slot;
                    let present: Vec<usize> = slot.map_or(Vec::new(), |slot| self.active[slot].segs.iter().map(|s| s.seg).collect());
                    modes -= self.keys[ci].segs.iter().filter(|si| !present.contains(si))
                        .map(|&si| if self.segs[si].duplex { d } else { n }).sum::<usize>();
                    keys -= usize::from(slot.is_none());
                    self.keys[ci].held = false;
                    self.keys[ci].refused = true;
                    self.counters.refused_note_ons += 1;
                    self.last_error = Some(format!("key {}: no room in the dynamic mode budget; strike refused", self.keys[ci].midi));
                }
            }
            // Retiring whole keys shifts slots: do it first, then activate.
            if !retire.is_empty() {
                if let Err(e) = self.rebuild(&[], &retire) {
                    self.last_error = Some(e);
                }
            }
            if !activate.is_empty() {
                if let Err(e) = self.rebuild(&activate, &[]) {
                    self.last_error = Some(e);
                }
            }
        }
        for (ci, velocity) in pending {
            if std::mem::take(&mut self.keys[ci].refused) {
                continue;
            }
            let speaking = self.keys[ci].segs.iter().filter(|&&s| !self.segs[s].duplex).count();
            let ready = self.keys[ci].slot.is_some_and(|slot| self.active[slot].segs.iter().filter(|s| !s.duplex).count() == speaking);
            if !ready {
                self.counters.refused_note_ons += 1;
                continue;
            }
            if self.keys[ci].active {
                // The hammer has not come back yet: relaunch as soon as it is clear.
                self.keys[ci].pending_velocity = Some(velocity);
            } else {
                self.launch(ci, velocity);
            }
        }
    }

    fn contact_position(&self, ac: &ACourse, s: &ASeg) -> f64 {
        let b = dot(&ac.gl, &self.bq);
        dot(&self.m.hs[s.off..s.off + s.n], &self.q[s.off..s.off + s.n]) + s.hammer_lift * b
    }

    fn launch(&mut self, ci: usize, velocity: f64) {
        let Some(slot) = self.keys[ci].slot else { return; };
        let ac = &self.active[slot];
        // engine.rs `arm`: keep 2 mm clearance below the lowest felt face.
        let mut launch = -LAUNCH_CLEARANCE;
        let mut positions = [0.0_f64; 3];
        for s in ac.segs.iter().filter(|s| !s.duplex) {
            let x = self.contact_position(ac, s);
            positions[s.member] = x;
            let key = &self.keys[ci];
            let c = &key.contacts[s.member];
            let free = key.law.eps_residual(&c.state) * key.thickness + key.creep.deformation(&c.memory);
            launch = launch.min(x + free - LAUNCH_CLEARANCE);
        }
        let key = &mut self.keys[ci];
        let before = key.model.energy(&key.motion, CATCH_DISTANCE);
        key.motion = key.model.launch(launch, velocity);
        key.active = true;
        key.relaxing = false;
        let y = key.model.position(&key.motion);
        for (member, c) in key.contacts.iter_mut().enumerate() {
            c.overlap = y - positions[member];
            c.force = 0.0;
            c.stored_j = 0.0;
        }
        self.ledger.input_j += key.model.energy(&key.motion, CATCH_DISTANCE) - before;
        if !self.busy.contains(&ci) {
            self.busy.push(ci);
        }
    }

    pub fn note_off(&mut self, midi: u8) {
        if let Some(ci) = self.key_index(midi) {
            self.keys[ci].held = false;
        }
    }

    pub fn all_notes_off(&mut self) {
        self.pending_on.clear();
        for key in &mut self.keys {
            key.held = false;
            key.pending_velocity = None;
        }
        self.sustain = 0.0;
    }

    fn refresh_pads(&mut self) {
        let engagement = (1.0 - self.sustain).powi(2);
        if engagement == self.engagement_cached {
            return;
        }
        self.engagement_cached = engagement;
        let half_dt = 0.5 * self.dt;
        for ac in &mut self.active {
            for s in &mut ac.segs {
                let count = s.pad.len() as f64;
                for p in &mut s.pad {
                    p.norm = p.norm_string + p.lift * p.lift * ac.gl_norm2;
                    let drag = DAMPER_DRAG_NS_M * engagement / count;
                    p.coef = if p.norm == 0.0 || drag == 0.0 { 0.0 } else { det::expm1(-drag * p.norm * half_dt) / p.norm };
                }
            }
        }
    }

    /// Exact rank-one viscous flows of every engaged pad point for half a
    /// step (`dampers.rs::rank_one`); the post half reverses the order.
    /// The flows stay sequential in exact arithmetic; per segment they are
    /// fused (`damp_fused`), and the bridge velocity change of a whole course
    /// is applied to `bvel` once (it enters later flows of that course only
    /// through the scalar `gl . bvel`, tracked exactly).
    fn damp(&mut self, reverse: bool) -> f64 {
        if self.sustain >= 1.0 {
            return 0.0;
        }
        let mut loss = 0.0;
        let count = self.active.len();
        for index in 0..count {
            let slot = if reverse { count - 1 - index } else { index };
            let ac = &self.active[slot];
            if !ac.damped_now || ac.segs.iter().all(|s| s.pad.is_empty()) {
                continue;
            }
            let mut beta = dot(&ac.gl, &self.bvel);
            let mut push = 0.0;
            let nseg = ac.segs.len();
            for sindex in 0..nseg {
                let s = &ac.segs[if reverse { nseg - 1 - sindex } else { sindex }];
                let np = s.pad.len();
                let v = &mut self.v[s.off..s.off + s.n];
                let g2 = ac.gl_norm2;
                loss += match np {
                    0 => 0.0,
                    1 => damp_fused::<1>(s, v, g2, &mut beta, &mut push, reverse),
                    2 => damp_fused::<2>(s, v, g2, &mut beta, &mut push, reverse),
                    3 => damp_fused::<3>(s, v, g2, &mut beta, &mut push, reverse),
                    4 => damp_fused::<4>(s, v, g2, &mut beta, &mut push, reverse),
                    _ => {
                        // More than four pad points: the plain sequential flows.
                        for (x, g) in self.bvel.iter_mut().zip(&ac.gl) {
                            *x += g * push;
                        }
                        push = 0.0;
                        let mut seg_loss = 0.0;
                        for pindex in 0..np {
                            let p = &s.pad[if reverse { np - 1 - pindex } else { pindex }];
                            if p.coef == 0.0 {
                                continue;
                            }
                            let speed = dot(&p.shape, v) + p.lift * dot(&ac.gl, &self.bvel);
                            let change = p.coef * speed;
                            for (x, g) in v.iter_mut().zip(&p.shape) {
                                *x += g * change;
                            }
                            let lc = p.lift * change;
                            for (x, g) in self.bvel.iter_mut().zip(&ac.gl) {
                                *x += g * lc;
                            }
                            seg_loss += -change * speed - 0.5 * change * change * p.norm;
                        }
                        beta = dot(&ac.gl, &self.bvel);
                        seg_loss
                    }
                };
            }
            if push != 0.0 {
                for (x, g) in self.bvel.iter_mut().zip(&ac.gl) {
                    *x += g * push;
                }
            }
        }
        loss
    }

    /// One mechanical step = one audio sample. Returns pressure [Pa].
    pub fn step(&mut self) -> f64 {
        if !self.pending_on.is_empty() {
            self.apply_pending_notes();
        }
        let r = self.r;
        let dur = 0.5 * self.dt;
        self.refresh_pads();
        for slot in 0..self.active.len() {
            let damped = self.is_damped(&self.active[slot]);
            self.active[slot].damped_now = damped;
            self.active[slot].age += 1;
        }
        let mut damper_loss = self.damp(false);
        let mut string_loss = 0.0;
        let mut board_loss = 0.0;

        // ---- intrinsic string loss, pre half-flow (implicit midpoint) ----
        self.rhs.copy_from_slice(&self.bvel);
        for ac in &self.active {
            let mut source = 0.0;
            for s in &ac.segs {
                source += dot(&self.m.cross[s.off..s.off + s.n], &self.v[s.off..s.off + s.n]);
            }
            for (x, g) in self.rhs.iter_mut().zip(&ac.gl) {
                *x += g * source;
            }
        }
        matvec(&self.sd_inv, &self.rhs, &mut self.mid);
        self.force_b.fill(0.0);
        self.rhs.fill(0.0); // becomes the coupling right-hand side below
        for ac in &mut self.active {
            let bm = dot(&ac.gl, &self.mid);
            let b0 = dot(&ac.gl, &self.bq);
            ac.b0 = b0;
            let mut reaction_loss = 0.0;
            let mut reaction = 0.0;
            for s in &mut ac.segs {
                let range = s.off..s.off + s.n;
                let (mut loss, mut acc3, mut reac) = (0.0, 0.0, 0.0);
                let m = &self.m;
                let it = self.q[range.clone()].iter()
                    .zip(self.v[range.clone()].iter())
                    .zip(self.fq[range.clone()].iter_mut())
                    .zip(self.fv[range.clone()].iter_mut())
                    .zip(m.pre[range.clone()].iter());
                for ((((q, v), fq), fv), p) in it {
                    let old = *v;
                    let dv = p.c1 * old + p.c2 * bm;
                    loss += (old - dv) * (old + dv);
                    acc3 += p.c3 * old;
                    let x = p.qq * q + p.qv * dv;
                    *fq = x;
                    *fv = p.vq * q + p.vv * dv;
                    reac += p.ha * (q + x);
                }
                string_loss += 0.5 * loss;
                reaction_loss += s.endpoint_drag * bm - acc3 - bm * s.sum_c4;
                reaction += reac + b0 * s.w_schur;
            }
            for ((f, x), g) in self.force_b.iter_mut().zip(self.rhs.iter_mut()).zip(&ac.gl) {
                *f += g * reaction_loss;
                *x += g * reaction;
            }
        }
        for (v, f) in self.bvel.iter_mut().zip(&self.force_b) {
            let old = *v;
            *v = old - dur * f;
            string_loss += 0.5 * (old - *v) * (old + *v);
        }
        // ---- board loss, pre half-flow (exact matrix exponential) ----
        matvec(&self.flow, &self.bvel, &mut self.tmp);
        board_loss += 0.5 * self.bvel.iter().zip(&self.tmp).map(|(a, b)| a * a - b * b).sum::<f64>();
        // ---- exact-ZOH free motion and the coupled unforced end position ----
        for j in 0..r {
            let z = self.bz[j];
            let dv = self.tmp[j];
            self.bfq[j] = z.qq * self.bq[j] + z.qv * dv;
            self.bfv[j] = z.vq * self.bq[j] + z.vv * dv;
            self.rhs[j] += self.bfq[j] / z.bq;
        }
        matvec(&self.s_inv, &self.rhs, &mut self.end);

        // ---- hammers: free flight, implicit felt contact ----
        let mut felt_loss = 0.0;
        let mut shank_loss = 0.0;
        let mut catch_loss = 0.0;
        let any_hammer = self.busy.iter().any(|&ci| self.keys[ci].active);
        if any_hammer {
            self.solve_contacts();
        } else {
            for ac in &mut self.active {
                for s in &mut ac.segs {
                    s.force = 0.0;
                }
            }
        }

        // ---- finish with held contact forces ----
        // Board end position: unforced Schur solution plus contact responses.
        self.tmp.copy_from_slice(&self.end);
        if any_hammer {
            for ac in &self.active {
                let push: f64 = ac.segs.iter().map(|s| s.lift_total * s.force).sum();
                if push != 0.0 {
                    for (x, g) in self.tmp.iter_mut().zip(&ac.resp) {
                        *x += g * push;
                    }
                }
            }
        }
        self.rhs.fill(0.0);
        for ac in &mut self.active {
            let bbar = 0.5 * (ac.b0 + dot(&ac.gl, &self.tmp));
            let mut reaction = 0.0;
            for s in &mut ac.segs {
                let range = s.off..s.off + s.n;
                let m = &self.m;
                let (mut reac, mut src) = (0.0, 0.0);
                if s.force == 0.0 {
                    let it = self.q[range.clone()].iter_mut()
                        .zip(self.v[range.clone()].iter_mut())
                        .zip(self.fq[range.clone()].iter())
                        .zip(self.fv[range.clone()].iter())
                        .zip(m.fin[range.clone()].iter());
                    for ((((q, v), fq), fv), p) in it {
                        let nq = fq + p.bqa * bbar;
                        let nv = fv + p.bva * bbar;
                        reac += p.ha * (*q + nq);
                        src += p.cross * nv;
                        *q = nq;
                        *v = nv;
                    }
                } else {
                    let force = s.force;
                    let it = self.q[range.clone()].iter_mut()
                        .zip(self.v[range.clone()].iter_mut())
                        .zip(self.fq[range.clone()].iter())
                        .zip(self.fv[range.clone()].iter())
                        .zip(m.fin[range.clone()].iter())
                        .zip(m.hit[range.clone()].iter());
                    for (((((q, v), fq), fv), p), h) in it {
                        let nq = fq + p.bqa * bbar + h.bqhs * force;
                        let nv = fv + p.bva * bbar + h.bvhs * force;
                        reac += p.ha * (*q + nq);
                        src += p.cross * nv;
                        *q = nq;
                        *v = nv;
                    }
                    reac += s.hammer_lift * force;
                }
                reaction += reac;
                s.src = src;
            }
            for (x, g) in self.rhs.iter_mut().zip(&ac.gl) {
                *x += g * reaction;
            }
        }
        for j in 0..r {
            self.bq[j] = self.tmp[j];
            self.bvel[j] = self.bfv[j] + self.bz[j].bv * self.rhs[j];
        }
        // ---- hammer/felt/creep commit against the accepted string motion ----
        if !self.busy.is_empty() {
            let (f, s, c) = self.advance_hammers();
            felt_loss += f;
            shank_loss += s;
            catch_loss += c;
        }
        // ---- board loss, post half-flow ----
        matvec(&self.flow, &self.bvel, &mut self.tmp);
        board_loss += 0.5 * self.bvel.iter().zip(&self.tmp).map(|(a, b)| a * a - b * b).sum::<f64>();
        self.bvel.copy_from_slice(&self.tmp);
        // ---- intrinsic string loss, post half-flow ----
        self.rhs.copy_from_slice(&self.bvel);
        for ac in &self.active {
            let source: f64 = ac.segs.iter().map(|s| s.src).sum();
            for (x, g) in self.rhs.iter_mut().zip(&ac.gl) {
                *x += g * source;
            }
        }
        matvec(&self.sd_inv, &self.rhs, &mut self.mid);
        self.force_b.fill(0.0);
        for ac in &self.active {
            let bm = dot(&ac.gl, &self.mid);
            let mut reaction_loss = 0.0;
            for s in &ac.segs {
                let range = s.off..s.off + s.n;
                let m = &self.m;
                let (mut loss, mut acc3) = (0.0, 0.0);
                let it = self.v[range.clone()].iter_mut()
                    .zip(m.c1[range.clone()].iter())
                    .zip(m.c2[range.clone()].iter())
                    .zip(m.c3[range.clone()].iter());
                for (((v, c1), c2), c3) in it {
                    let old = *v;
                    let dv = c1 * old + c2 * bm;
                    loss += (old - dv) * (old + dv);
                    acc3 += c3 * old;
                    *v = dv;
                }
                string_loss += 0.5 * loss;
                reaction_loss += s.endpoint_drag * bm - acc3 - bm * s.sum_c4;
            }
            for (f, g) in self.force_b.iter_mut().zip(&ac.gl) {
                *f += g * reaction_loss;
            }
        }
        for (v, f) in self.bvel.iter_mut().zip(&self.force_b) {
            let old = *v;
            *v = old - dur * f;
            string_loss += 0.5 * (old - *v) * (old + *v);
        }
        damper_loss += self.damp(true);

        self.ledger.string_loss_j += string_loss;
        self.ledger.board_loss_j += board_loss;
        self.ledger.damper_loss_j += damper_loss;
        self.ledger.felt_loss_j += felt_loss;
        self.ledger.shank_loss_j += shank_loss;
        self.ledger.catch_loss_j += catch_loss;

        // ---- Rayleigh receiver on the physical (bare) board velocity ----
        matvec(&self.phi, &self.bvel, &mut self.st_v);
        let pressure = match self.mic.step_trace(&self.st_v) {
            Ok(p) => p,
            Err(e) => {
                self.last_error = Some(e);
                0.0
            }
        };
        self.counters.samples += 1;
        let magnitude = pressure.abs();
        if magnitude > self.counters.peak_pa {
            self.counters.peak_pa = magnitude;
        }
        if magnitude >= self.cfg.full_scale_pa {
            self.counters.clipped_samples += 1;
        }
        self.retire_clock += 1;
        if self.retire_clock >= 512 {
            self.retire_clock = 0;
            self.retire_quiet();
        }
        pressure
    }

    /// Predict the free hammer positions and solve the felt forces of every
    /// flying hammer against the reciprocal string/board/shank compliance
    /// (`string_engine.rs::solve_contacts_for_trial` + `scalar_contact`).
    fn solve_contacts(&mut self) {
        self.c_key.clear();
        self.c_member.clear();
        self.c_gap.clear();
        self.c_force.clear();
        self.c_lift.clear();
        for ac in &mut self.active {
            for s in &mut ac.segs {
                s.force = 0.0;
            }
        }
        // Free hammer flight for the step, then the string-side free position.
        let busy = std::mem::take(&mut self.busy);
        for &ci in &busy {
            let key = &mut self.keys[ci];
            if !key.active {
                continue;
            }
            let Some(slot) = key.slot else { continue; };
            // The unforced step of the shank is kept: if no felt force arises
            // it IS the accepted step (same call, same arguments).
            key.predicted = None;
            let free = match key.model.advance(&key.motion, 0.0, 0.0) {
                Ok((state, loss)) => {
                    key.predicted = Some((state, loss));
                    key.model.position(&state)
                }
                Err(e) => {
                    self.last_error = Some(e.into());
                    key.model.position(&key.motion)
                }
            };
            let ac = &self.active[slot];
            let b1 = dot(&ac.gl, &self.end);
            for s in ac.segs.iter().filter(|s| !s.duplex) {
                // sum(hs*(fq + hba*(b0+b1))) + hammer_lift*b1
                let lift_bq = s.lift_total - s.hammer_lift;
                let p_free = dot(&self.m.hs[s.off..s.off + s.n], &self.fq[s.off..s.off + s.n]);
                let position = p_free + (ac.b0 + b1) * lift_bq + s.hammer_lift * b1;
                let old = &key.contacts[s.member];
                // Felt that starts and ends this step clear of its crushed
                // face carries exactly zero force: the scalar solve would
                // return 0 at its first residual. Such sites need no solve.
                let face = key.law.eps_residual(&old.state) * key.thickness;
                let start = old.overlap - key.creep.deformation(&old.memory);
                let end = free - position - key.creep.free_deformation(&old.memory);
                if old.force == 0.0 && start <= face && end <= face && start <= 0.0 && end <= 0.0 {
                    continue;
                }
                self.c_key.push(ci);
                self.c_member.push(s.member);
                self.c_gap.push(free - position);
                self.c_force.push(old.force);
                self.c_lift.push(s.lift_total);
            }
        }
        self.busy = busy;
        let nc = self.c_key.len();
        if nc == 0 {
            return;
        }
        // Reciprocal contact compliance among the active contacts.
        self.c_h.clear();
        self.c_h.resize(nc * nc, 0.0);
        for i in 0..nc {
            let (ki, mi) = (self.c_key[i], self.c_member[i]);
            let ai = &self.active[self.keys[ki].slot.unwrap_or(0)];
            for j in i..nc {
                let kj = self.c_key[j];
                let aj = &self.active[self.keys[kj].slot.unwrap_or(0)];
                let mut h = self.c_lift[i] * self.c_lift[j] * dot(&ai.gl, &aj.resp);
                if ki == kj {
                    h += self.keys[ki].model.compliance();
                    if mi == self.c_member[j] {
                        let s = ai.segs.iter().find(|s| !s.duplex && s.member == mi);
                        h += s.map_or(0.0, |s| s.cdiag);
                    }
                }
                self.c_h[i * nc + j] = h;
                self.c_h[j * nc + i] = h;
            }
        }
        for i in 0..nc {
            for j in 0..nc {
                self.c_gap[i] -= self.c_h[i * nc + j] * self.c_force[j];
            }
        }
        // Prepare each site once: law, committed crush branch, start overlap.
        self.c_site.clear();
        self.c_start.clear();
        for i in 0..nc {
            let key = &self.keys[self.c_key[i]];
            let old = &key.contacts[self.c_member[i]];
            let start = old.overlap - key.creep.deformation(&old.memory);
            self.c_start.push(key.creep.free_deformation(&old.memory));
            self.c_site.push(crate::fastfelt::Site::new(&key.law, &old.state, start, key.thickness, key.area));
        }
        let mut converged = false;
        let mut failed: Option<usize> = None;
        'sweeps: for _ in 0..32 {
            for i in 0..nc {
                let key = &self.keys[self.c_key[i]];
                let diagonal = self.c_h[i * nc + i];
                let free = self.c_gap[i] + diagonal * self.c_force[i] - self.c_start[i];
                match self.c_site[i].solve(free, diagonal + key.creep.compliance(), self.c_force[i]) {
                    Ok(next) => {
                        let change = next - self.c_force[i];
                        self.c_force[i] = next;
                        if change != 0.0 {
                            for j in 0..nc {
                                self.c_gap[j] -= self.c_h[j * nc + i] * change;
                            }
                        }
                    }
                    Err(e) => {
                        self.last_error = Some(format!("key {}: {e}", key.midi));
                        failed = Some(self.c_key[i]);
                        break 'sweeps;
                    }
                }
            }
            converged = true;
            for i in 0..nc {
                let key = &self.keys[self.c_key[i]];
                let end = self.c_gap[i] - self.c_start[i] - key.creep.compliance() * self.c_force[i];
                let expected = self.c_site[i].force(end).0;
                if !expected.is_finite() || (self.c_force[i] - expected).abs() > 1e-5 + 1e-8 * expected.abs() {
                    converged = false;
                }
            }
            if converged {
                break;
            }
        }
        if failed.is_none() && !converged {
            self.last_error = Some("coupled contact solve exhausted 32 sweeps".into());
            failed = Some(self.c_key[0]);
        }
        if let Some(ci) = failed {
            // The strike cannot be resolved at this clock: withdraw that hammer
            // (its energy goes to the catch books) rather than invent a force.
            self.counters.contact_failures += 1;
            for i in 0..nc {
                if self.c_key[i] == ci || !converged {
                    self.c_force[i] = 0.0;
                }
            }
            let key = &mut self.keys[ci];
            let felt_stored: f64 = key.contacts.iter().map(|c| c.stored_j + key.creep.stored(&c.memory)).sum();
            self.ledger.catch_loss_j += key.model.energy(&key.motion, CATCH_DISTANCE) - key.rest_energy + felt_stored;
            key.motion = key.rest;
            key.active = false;
            key.relaxing = false;
            key.pending_velocity = None;
            for c in &mut key.contacts {
                c.memory = relaxation::Memory::default();
                c.overlap = -CATCH_DISTANCE;
                c.force = 0.0;
                c.stored_j = 0.0;
            }
        }
        for i in 0..nc {
            let (ci, member) = (self.c_key[i], self.c_member[i]);
            if let Some(slot) = self.keys[ci].slot {
                if let Some(s) = self.active[slot].segs.iter_mut().find(|s| !s.duplex && s.member == member) {
                    s.force = self.c_force[i];
                }
            }
        }
    }

    /// Advance flying hammers with the accepted forces, commit felt crush and
    /// Prony memory, catch landed hammers (engine.rs `mechanics_step`).
    fn advance_hammers(&mut self) -> (f64, f64, f64) {
        let (mut felt_loss, mut shank_loss, mut catch_loss) = (0.0, 0.0, 0.0);
        let mut busy = std::mem::take(&mut self.busy);
        let mut relaunch: Vec<(usize, f64)> = Vec::new();
        for &ci in &busy {
            let slot = self.keys[ci].slot;
            // String-side contact positions after the accepted step.
            let mut positions = [0.0_f64; 3];
            let mut forces = [0.0_f64; 3];
            if let Some(slot) = slot {
                let ac = &self.active[slot];
                for s in ac.segs.iter().filter(|s| !s.duplex) {
                    positions[s.member] = self.contact_position(ac, s);
                    forces[s.member] = s.force;
                }
            }
            let key = &mut self.keys[ci];
            if key.active {
                let total: f64 = forces.iter().sum();
                key.last_force = total;
                let advanced = match key.predicted.take() {
                    Some(step) if total == 0.0 => Ok(step),
                    _ => key.model.advance(&key.motion, -total, 0.0),
                };
                match advanced {
                    Ok((motion, loss)) => {
                        key.motion = motion;
                        shank_loss += loss;
                    }
                    Err(e) => self.last_error = Some(e.into()),
                }
            } else {
                key.last_force = 0.0;
            }
            let y = key.model.position(&key.motion);
            let mut memory_left = 0.0_f64;
            for (member, old) in key.contacts.iter_mut().enumerate() {
                let force = if key.active { forces[member] } else { 0.0 };
                let (memory, viscous) = key.creep.advance(&old.memory, force);
                felt_loss += viscous;
                if key.active {
                    let end = y - positions[member];
                    let start_elastic = old.overlap - key.creep.deformation(&old.memory);
                    let end_elastic = end - key.creep.deformation(&memory);
                    // fs-material commits the crush history; the recoverable
                    // energy at the committed point is carried from step to step.
                    let state = key.law.update_state(end_elastic / key.thickness, &old.state);
                    let stored = key.area * key.thickness
                        * crate::fastfelt::stored(&key.law, end_elastic / key.thickness, &state);
                    felt_loss += force * (end_elastic - start_elastic) - (stored - old.stored_j);
                    old.stored_j = stored;
                    old.state = state;
                    old.overlap = end;
                }
                old.memory = memory;
                old.force = force;
                memory_left = memory_left.max(key.creep.deformation(&old.memory).abs());
            }
            if key.active {
                let clear = key.contacts.iter().all(|c| c.overlap <= 0.0);
                let landed = y < -CATCH_DISTANCE - 1e-10 && key.model.velocity(&key.motion) < 0.0;
                if clear && (landed || key.pending_velocity.is_some()) {
                    catch_loss += key.model.energy(&key.motion, CATCH_DISTANCE) - key.rest_energy;
                    key.motion = key.rest;
                    key.active = false;
                    key.relaxing = true;
                    key.last_force = 0.0;
                    for c in &mut key.contacts {
                        c.overlap = -CATCH_DISTANCE;
                        c.force = 0.0;
                        c.stored_j = 0.0;
                    }
                    if let Some(velocity) = key.pending_velocity.take() {
                        relaunch.push((ci, velocity));
                    }
                }
            }
            if !key.active && key.relaxing && memory_left < 1.0e-12 {
                // Residual creep energy is below resolution: book it as felt loss.
                for c in &mut key.contacts {
                    felt_loss += key.creep.stored(&c.memory);
                    c.memory = relaxation::Memory::default();
                }
                key.relaxing = false;
            }
        }
        busy.retain(|&ci| self.keys[ci].active || self.keys[ci].relaxing);
        self.busy = busy;
        for (ci, velocity) in relaunch {
            if self.keys[ci].slot.is_some() {
                self.launch(ci, velocity);
            }
        }
        (felt_loss, shank_loss, catch_loss)
    }

    /// Return quiet strings to the static image (cold, at most one rebuild).
    /// Damped speaking strings and the never-damped duplex segments of one
    /// key retire independently.
    fn retire_quiet(&mut self) {
        if self.active.is_empty() || !self.cfg.static_attachment {
            return;
        }
        let mut retire = Vec::new();
        for (slot, ac) in self.active.iter().enumerate().rev() {
            let key = &self.keys[ac.ci];
            if key.active || key.pending_velocity.is_some() || ac.age < 2048 {
                continue;
            }
            let damped = self.is_damped(ac);
            let b = dot(&ac.gl, &self.bq);
            let bdot = dot(&ac.gl, &self.bvel);
            for s in &ac.segs {
                let threshold = if damped && !s.duplex { self.cfg.retire_damped_j } else { self.cfg.retire_ringing_j };
                if self.segment_energy(s, b, bdot) < threshold {
                    retire.push((slot, s.seg));
                }
            }
        }
        if !retire.is_empty() {
            if let Err(e) = self.rebuild(&[], &retire) {
                self.last_error = Some(e);
            }
        }
    }

    /// Mechanical energy now stored in strings, board, hammers and felt [J].
    pub fn stored_energy_j(&self) -> f64 {
        let r = self.r;
        let mut e = 0.5 * self.bvel.iter().map(|x| x * x).sum::<f64>();
        for i in 0..r {
            e += 0.5 * self.bq[i] * dot(&self.k_phys[i * r..(i + 1) * r], &self.bq);
        }
        for ac in &self.active {
            let b = dot(&ac.gl, &self.bq);
            for s in &ac.segs {
                for k in s.off..s.off + s.n {
                    let x = self.q[k] - self.m.beta[k] * b;
                    e += 0.5 * self.v[k] * self.v[k] + 0.5 * self.m.w2[k] * x * x;
                }
            }
        }
        for &ci in &self.busy {
            let key = &self.keys[ci];
            e += key.model.energy(&key.motion, CATCH_DISTANCE) - key.rest_energy;
            for c in &key.contacts {
                e += key.creep.stored(&c.memory);
                if key.active {
                    e += c.stored_j;
                }
            }
        }
        e
    }

    /// Refresh `stored_j` and `residual_j` (input - stored - dissipated - retired).
    pub fn close_ledger(&mut self) -> Ledger {
        self.ledger.stored_j = self.stored_energy_j();
        self.ledger.residual_j = self.ledger.input_j - self.ledger.stored_j - self.ledger.dissipated_j() - self.ledger.retired_j;
        self.ledger
    }

    /// Render `out.len()` samples of pressure divided by the fixed full scale.
    pub fn render(&mut self, out: &mut [f32]) {
        let scale = 1.0 / self.cfg.full_scale_pa;
        for x in out.iter_mut() {
            *x = (self.step() * scale) as f32;
        }
    }

    // ------------------------------------------------------------------
    // Observation (display only; never feeds back into the mechanics).
    // ------------------------------------------------------------------

    /// Physical bare-board (FE) modal displacements, mass-normalized.
    pub fn board_coordinates(&self, out: &mut [f64]) {
        let mut st = vec![0.0; self.r];
        matvec(&self.phi, &self.bq, &mut st);
        matvec(&self.phi_s, &st, out);
    }
    /// Physical bridge velocity under a key [m/s] (sounding or not).
    pub fn bridge_velocity(&self, midi: u8) -> Option<f64> {
        let ci = self.key_index(midi)?;
        let g = &self.g_st[ci * self.r..(ci + 1) * self.r];
        Some(dot(g, &self.st_v))
    }
    /// Total felt force of a key's hammer during the last step [N].
    pub fn hammer_force(&self, midi: u8) -> Option<f64> {
        self.key_index(midi).map(|ci| self.keys[ci].last_force)
    }
    /// Hammer crown position [m] relative to the undeflected string plane.
    pub fn hammer_position(&self, midi: u8) -> Option<f64> {
        self.key_index(midi).map(|ci| self.keys[ci].model.position(&self.keys[ci].motion))
    }
    pub fn stations() -> usize { SHAPE_STATIONS }

    /// Flat snapshot for drawing. Layout (f32):
    /// `[n_keys, stations, board_modes, sustain, ledger x10, counters x4,`
    /// ` board coordinates (bare, m sqrt(kg)) x board_modes,`
    /// ` per key: midi, held, damped, hammer_y_m, hammer_force_n, relative_energy_j,`
    /// `          partials, strings, bridge_m, y(x_j) x stations ]`.
    pub fn snapshot(&mut self, out: &mut Vec<f32>) {
        out.clear();
        let ledger = self.close_ledger();
        out.extend_from_slice(&[self.active.len() as f32, SHAPE_STATIONS as f32, self.r as f32, self.sustain as f32]);
        out.extend_from_slice(&[
            ledger.input_j as f32, ledger.stored_j as f32, ledger.dissipated_j() as f32, ledger.residual_j as f32,
            ledger.string_loss_j as f32, ledger.board_loss_j as f32, ledger.damper_loss_j as f32,
            (ledger.felt_loss_j + ledger.shank_loss_j) as f32, ledger.catch_loss_j as f32, ledger.retired_j as f32,
        ]);
        out.extend_from_slice(&[
            self.dynamic_modes() as f32, self.dynamic_strings() as f32, self.counters.peak_pa as f32,
            self.counters.forced_retirements as f32,
        ]);
        let mut bare = vec![0.0; self.r];
        self.board_coordinates(&mut bare);
        out.extend(bare.iter().map(|x| *x as f32));
        for slot in 0..self.active.len() {
            let (ci, n) = (self.active[slot].ci, self.active[slot].shape_n);
            if n > 0 && self.shape_cache[ci].as_ref().is_none_or(|(count, _)| *count != n) {
                let c = &self.courses[ci];
                let root = (0.5 * c.linear_density_kg_m * c.length_m).sqrt();
                let mut shape = Vec::with_capacity(SHAPE_STATIONS * n);
                for j in 0..SHAPE_STATIONS {
                    let x = (j as f64 + 0.5) / SHAPE_STATIONS as f64;
                    for p in 1..=n {
                        shape.push(((p as f64 * PI * x).sin() / root) as f32);
                    }
                }
                self.shape_cache[ci] = Some((n, shape));
            }
        }
        for ac in &self.active {
            let key = &self.keys[ac.ci];
            let b = dot(&ac.gl, &self.bq);
            out.push(f32::from(key.midi));
            out.push(if key.held { 1.0 } else { 0.0 });
            out.push(if self.is_damped(ac) { 1.0 } else { 0.0 });
            out.push(key.model.position(&key.motion) as f32);
            out.push(key.last_force as f32);
            out.push(self.relative_energy(ac) as f32);
            out.push(ac.shape_n as f32);
            out.push(ac.segs.iter().filter(|s| !s.duplex).count() as f32);
            out.push(b as f32);
            let n = ac.shape_n;
            let shape = if n > 0 {
                self.shape_cache[ac.ci].as_ref().map(|(_, s)| s.as_slice())
            } else {
                None
            };
            for j in 0..SHAPE_STATIONS {
                let x = (j as f64 + 0.5) / SHAPE_STATIONS as f64;
                let mut y = x * b;
                if let Some(table) = shape {
                    let row = &table[j * n..(j + 1) * n];
                    for (p, g) in row.iter().enumerate() {
                        let k = ac.shape_off + p;
                        y += f64::from(*g) * (self.q[k] - self.m.beta[k] * b);
                    }
                }
                out.push(y as f32);
            }
        }
    }

    /// Partial frequencies [Hz] the model would retain for a key (member 0).
    pub fn partials_hz(&self, midi: u8) -> Vec<f64> {
        let Some(ci) = self.key_index(midi) else { return Vec::new(); };
        let (n, _) = self.grant(ci);
        let seg = &self.segs[self.keys[ci].segs[0]];
        self.m.omega[seg.off..seg.off + n].iter().map(|w| w / TAU).collect()
    }
    /// (speaking partials, duplex partials, dynamic modes) granted to a key.
    pub fn grant_for(&self, midi: u8) -> Option<(usize, usize, usize)> {
        let ci = self.key_index(midi)?;
        let (n, d) = self.grant(ci);
        let key = &self.keys[ci];
        let speaking = key.segs.iter().filter(|&&s| !self.segs[s].duplex).count();
        let duplex = key.segs.len() - speaking;
        Some((n, d, speaking * n + if self.cfg.duplex { duplex * d } else { 0 }))
    }
}
