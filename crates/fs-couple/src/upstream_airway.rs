//! The airway upstream of a valve: a passive linear duct between the valve junction and a pressure source at its far
//! end, as the valve sees it sample by sample. The glottis is the first user
//! ([`crate::glottis::GlottalIsland::with_upstream`]): the trachea and the bronchial tree below the folds, the lungs at
//! the far end. A reed or lip junction can take the same seam with a mouth, a vocal tract and the tree behind it
//! ([`AirwayChain::ahead_of`], [`UpstreamPort`]); nothing here names a valve.
//!
//! **Geometry.** An [`AirwayChain`] is a list of [`AirwaySection`]s, valve side first. A section is `count` identical
//! cylinders in parallel, each `length` long and of `radius`, with an optional locally reacting wall pin on each tube's
//! perimeter. A symmetric bifurcating tree, generation `n` of `2^(n-1)` tubes, is exactly such a chain: every tube of a
//! generation has the same pressure at each end and an equal share of the flow, so a generation is one two-port whose
//! series impedance is a tube's over `count` and whose shunt admittance is a tube's times `count`. An asymmetric tree is
//! not a chain and cannot be written here.
//!
//! **One frequency.** Per tube, the gas wavenumber `k` and characteristic impedance `Z_c` come from
//! [`fs_duct::segment_wave`] under the chain's [`LossModel`] ([`LossModel::Bessel`], the Zwikker-Kosten solution at
//! every shear number, which the narrow tubes deep in a tree need; [`LossModel::Lossless`] for closed-form checks). The
//! telegraph pair is rebuilt from them as fs-duct rebuilds it, `z' = conj(-i k Z_c)` and `y' = conj(-i k / Z_c)`, and a
//! wall adds its shunt `conj(2 pi a / Z'_wall)` per metre ([`fs_phs::wall_admittance_per_metre`]) to `y'`. With
//! `count` tubes the pair is `z' / count` and `y' count`; then `gamma = sqrt(z' y')`, `Z = z' / gamma`, and the
//! section's chain matrix is `[[cosh gamma L, Z sinh gamma L], [sinh gamma L / Z, cosh gamma L]]` in the telegraph
//! convention. `M` is the product from the valve to the far end: `(p, U)` at the valve from `(p, U)` at the far end,
//! `U` flowing away from the valve. With the far end held at the source pressure `P_s` (the lungs) and `U_v` the flow
//! INTO the valve, `p = P_s / M22 - Z_in U_v` at the valve with `Z_in = M12 / M22` (`det M = 1`). Every value this
//! module returns is in the crate's `e^{-iwt}` convention, the conjugate of the telegraph one.
//!
//! **Waves at the valve.** Referenced to the lossless characteristic impedance of the first section, `Z_s = rho c /
//! (count pi a^2)`: `p = a + b` and `U_v = (a - b) / Z_s`, so `b` is the wave the valve sends down the airway and `a`
//! the one it receives. Then `a = R b + G P_s`, with the reflectance `R = (Z_in - Z_s) / (Z_in + Z_s)` and the source
//! transfer `G = 1 / (M22 + M12 / Z_s)`. At `w = 0` the matrix is `[[1, R_dc], [0, 1]]` with the Poiseuille sum
//! `R_dc = sum 8 mu L / (pi a^4 count)` (a wall with a stiffness admits nothing at DC), so a steady flow reads `p = P_s
//! - R_dc U_v`.
//!
//! **Time domain.** [`UpstreamAirway::realize`] makes both causal FIRs on the sample clock. `R` goes through the crate's
//! characteristic-line ladder, which sizes the taps until the tail energy and the first impedance peaks converge and
//! refuses a converged line that is not passive on the 4x-oversampled grid; its DC bin is the exact `(R_dc - Z_s) / (R_dc
//! + Z_s)`. `G` goes through the transfer ladder of the mouth-flow line ([`crate::driving_point::outlet`]: tail energy
//! and tap change), with the exact DC bin `Z_s / (Z_s + R_dc)`. A step is split in two so the valve solves its own
//! junction implicitly. [`UpstreamAirway::port`] returns the Thevenin source the valve sees this sample, `p = p_oc - Z
//! U_v`, with `A = sum_{k>=1} r_k b_{n-k} + sum_{k>=0} g_k P_{n-k}`, `p_oc = 2 A / (1 - r_0)` and `Z = Z_s (1 + r_0) /
//! (1 - r_0)`: the same-sample tap `r_0` is folded into the junction, not delayed. [`UpstreamAirway::commit`] takes the
//! flow the valve solved for, sends `b = (A - Z_s U_v) / (1 - r_0)` down the airway and returns the port pressure `p =
//! (1 + r_0) b + A`.
//!
//! **Weibel's tree.** [`WEIBEL_SAAB_2021_TABLE_1_B_MM`] is the symmetric model printed by Saab, Issa, Zeineddine,
//! Tartakovsky and Lakkis (2021), "Impact of airways geometry on transport of gases to blood", arXiv:2111.07627v1,
//! Table 1 column B, which they attribute to Weibel's morphometry of the human lung (their refs 78-80; the book is of
//! 1963): generations `n = 1` (the trachea) to `24`, length and diameter in mm, `2^(n-1)` tubes in generation `n`.
//! [`AirwayChain::weibel_saab_2021_b`] takes the first `generations` of it. How many a model needs is a convergence
//! question its caller answers and states. The far end of the last generation taken is held at the source pressure.
//!
//! **A duct ahead of the tree.** A reed or a pair of lips sits at the mouth, not at the folds: between it and the
//! trachea lie the mouth, the vocal tract and the glottis, held open. [`AirwayChain::ahead_of`] puts one chain in front
//! of another, valve side first, so a duct (each entry of an area function an [`AirwaySection::from_area`] cylinder of
//! the same area, the glottal opening a section of its own) followed by [`AirwayChain::weibel_saab_2021_b`] is one
//! chain, with the source pressure held at the tree's far end and `Z_s` the duct's first section's. Everything above
//! applies to it unchanged; the realization never asks where a section came from. No duct is supplied here: the shape of
//! a player's tract and glottis while playing is a measurement the caller brings and names (a sung vowel is not one).
//!
//! Limits, not hidden: cylinders with no wall or one locally reacting wall pin on every tube (no measured per-generation
//! wall is carried here); a symmetric tree (the real bronchial tree is not); plane waves in every tube, with no junction
//! mass or branching loss, and no end correction or separation loss where the area steps between two sections; a
//! stationary gas (no mean flow and no lung-volume change over a note), so an open glottis ahead of the tree carries no
//! jet loss; and the far end an ideal pressure source (no alveolar compliance).

use crate::driving_point::{
    DrivingPointError, LineDc, LineLadder, LineRealization, PoleState, realize_characteristic_response,
    realize_reflectance_ir_from_positive_spectrum,
};
use fs_duct::{LossModel, segment_wave};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;
use fs_phs::{WallPin, wall_admittance_per_metre};

/// Weibel's symmetric airway tree as Saab et al. (2021) print it in Table 1 column B (module doc): `(length,
/// diameter)` in mm for generations 1 (the trachea) to 24, with `2^(n-1)` tubes in generation `n`.
pub const WEIBEL_SAAB_2021_TABLE_1_B_MM: [(f64, f64); 24] = [
    (120.0, 18.0),
    (47.8, 12.21),
    (19.1, 8.3),
    (7.6, 5.6),
    (12.7, 4.5),
    (10.7, 3.5),
    (9.0, 2.8),
    (7.6, 2.3),
    (6.4, 1.9),
    (5.4, 1.5),
    (4.6, 1.3),
    (3.9, 1.1),
    (3.3, 1.0),
    (2.7, 0.8),
    (2.3, 0.7),
    (2.0, 0.7),
    (1.7, 0.6),
    (1.4, 0.5),
    (1.2, 0.5),
    (1.0, 0.5),
    (0.8, 0.5),
    (0.7, 0.4),
    (0.6, 0.4),
    (0.5, 0.4),
];

/// Largest DFT length the airway ladders accept, the characteristic line's own ceiling.
const MAX_AIRWAY_TAPS: usize = 131_072;

/// One generation of an airway: `count` identical cylinders in parallel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AirwaySection {
    /// Axial length of each tube [m].
    pub length_m: f64,
    /// Radius of each tube [m].
    pub radius_m: f64,
    /// Number of identical tubes in parallel, at least 1.
    pub count: u32,
    /// Locally reacting wall pin on each tube's perimeter; `None` is a rigid wall.
    pub wall: Option<WallPin>,
}

impl AirwaySection {
    /// `count` identical tubes `length_m` long, each of cross-sectional area `area_m2`: one entry of an area function as
    /// the cylinder of the same area, `radius = sqrt(area / pi)` (module doc, a duct ahead of the tree).
    ///
    /// # Errors
    /// A non-positive or non-finite length or area, or a zero count.
    pub fn from_area(length_m: f64, area_m2: f64, count: u32, wall: Option<WallPin>) -> Result<Self, DrivingPointError> {
        let ok = |x: f64| x > 0.0 && x.is_finite();
        let radius_m = (area_m2 / core::f64::consts::PI).sqrt();
        if !(ok(length_m) && ok(area_m2) && ok(radius_m) && count >= 1) {
            return Err(DrivingPointError::Invalid {
                what: "an airway section from an area needs a positive finite length and area and at least one tube",
            });
        }
        Ok(Self { length_m, radius_m, count, wall })
    }
}

/// An airway from the valve (first section) to the source pressure at the last section's far end.
#[derive(Debug, Clone, PartialEq)]
pub struct AirwayChain {
    sections: Vec<AirwaySection>,
    loss: LossModel,
}

/// `exp(z)` from the deterministic real kernels.
fn c_exp(z: C64) -> C64 {
    let m = det::exp(z.re);
    let (s, c) = det::sin_cos(z.im);
    C64::new(m * c, m * s)
}

/// `(cosh z, sinh z)`; a Taylor series below `|z| = 1e-2`, where `(e^z - e^-z) / 2` would cancel.
fn c_cosh_sinh(z: C64) -> (C64, C64) {
    if z.abs() < 1.0e-2 {
        let z2 = z * z;
        let ch = C64::ONE + z2 * (C64::from_re(0.5) + z2 * (C64::from_re(1.0 / 24.0) + z2.scale(1.0 / 720.0)));
        let sh = z * (C64::ONE + z2 * (C64::from_re(1.0 / 6.0) + z2 * (C64::from_re(1.0 / 120.0) + z2.scale(1.0 / 5040.0))));
        (ch, sh)
    } else {
        let e = c_exp(z);
        let inv = e.recip();
        ((e + inv).scale(0.5), (e - inv).scale(0.5))
    }
}

fn mat_mul(a: [[C64; 2]; 2], b: [[C64; 2]; 2]) -> [[C64; 2]; 2] {
    [
        [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
        [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
    ]
}

impl AirwayChain {
    /// A chain of sections, valve side first.
    ///
    /// # Errors
    /// No section, a non-positive or non-finite length or radius, a zero count, or a loss model without an exact DC
    /// limit ([`LossModel::WideTube`] is first order in the shear number and refuses the narrow tubes it would meet).
    pub fn new(sections: Vec<AirwaySection>, loss: LossModel) -> Result<Self, DrivingPointError> {
        if sections.is_empty() {
            return Err(DrivingPointError::Invalid { what: "an airway needs at least one section" });
        }
        if matches!(loss, LossModel::WideTube) {
            return Err(DrivingPointError::Invalid { what: "an airway takes the Bessel or the lossless gas model" });
        }
        let ok = |x: f64| x > 0.0 && x.is_finite();
        if sections.iter().any(|s| !(ok(s.length_m) && ok(s.radius_m) && s.count >= 1)) {
            return Err(DrivingPointError::Invalid {
                what: "airway sections need positive finite lengths and radii and at least one tube",
            });
        }
        Ok(Self { sections, loss })
    }

    /// The first `generations` of [`WEIBEL_SAAB_2021_TABLE_1_B_MM`] with Bessel losses, every tube carrying `wall`.
    ///
    /// # Errors
    /// `generations` outside `1..=24`.
    pub fn weibel_saab_2021_b(generations: usize, wall: Option<WallPin>) -> Result<Self, DrivingPointError> {
        if !(1..=WEIBEL_SAAB_2021_TABLE_1_B_MM.len()).contains(&generations) {
            return Err(DrivingPointError::Invalid { what: "Weibel's tree as printed has generations 1 to 24" });
        }
        let sections = WEIBEL_SAAB_2021_TABLE_1_B_MM[..generations]
            .iter()
            .enumerate()
            .map(|(n, &(length_mm, diameter_mm))| AirwaySection {
                length_m: length_mm * 1.0e-3,
                radius_m: 0.5 * diameter_mm * 1.0e-3,
                count: 1_u32 << n,
                wall,
            })
            .collect();
        Self::new(sections, LossModel::Bessel)
    }

    /// This chain with `far` beyond its far end, valve side first (module doc, a duct ahead of the tree): a duct at the
    /// valve ahead of the airway below it. The source pressure moves to `far`'s far end; `Z_s` stays this chain's first
    /// section's.
    ///
    /// # Errors
    /// The two chains carry different gas loss models.
    pub fn ahead_of(mut self, far: &Self) -> Result<Self, DrivingPointError> {
        if self.loss != far.loss {
            return Err(DrivingPointError::Invalid { what: "an airway joins only chains of one gas loss model" });
        }
        self.sections.extend_from_slice(&far.sections);
        Ok(self)
    }

    /// The sections, valve side first.
    #[must_use]
    pub fn sections(&self) -> &[AirwaySection] {
        &self.sections
    }

    /// The gas loss model of every tube.
    #[must_use]
    pub const fn loss(&self) -> LossModel {
        self.loss
    }

    /// `Z_s = rho c / (count pi a^2)` of the first section [Pa s/m^3], the reference of the port waves.
    #[must_use]
    pub fn port_impedance(&self, gas: &GasState) -> f64 {
        let first = &self.sections[0];
        gas.density * gas.sound_speed
            / (f64::from(first.count) * core::f64::consts::PI * first.radius_m * first.radius_m)
    }

    /// The DC resistance `R_dc = sum 8 mu L / (pi a^4 count)` [Pa s/m^3] (zero on a lossless chain).
    #[must_use]
    pub fn dc_resistance(&self, gas: &GasState) -> f64 {
        if matches!(self.loss, LossModel::Lossless) {
            return 0.0;
        }
        self.sections
            .iter()
            .map(|s| {
                let a2 = s.radius_m * s.radius_m;
                8.0 * gas.dynamic_viscosity * s.length_m / (core::f64::consts::PI * a2 * a2 * f64::from(s.count))
            })
            .sum()
    }

    /// Path length from the valve to the far end [m].
    #[must_use]
    pub fn path_length_m(&self) -> f64 {
        self.sections.iter().map(|s| s.length_m).sum()
    }

    fn telegraph_chain(&self, gas: &GasState, omega: f64) -> Result<[[C64; 2]; 2], DrivingPointError> {
        let minus_i = C64::new(0.0, -1.0);
        let mut m = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
        for s in &self.sections {
            let wave = segment_wave(gas, s.radius_m, omega, self.loss)?;
            let z_tube = (minus_i * wave.wavenumber * wave.characteristic_impedance).conj();
            let mut y_tube = (minus_i * wave.wavenumber * wave.characteristic_impedance.recip()).conj();
            if let Some(pin) = &s.wall {
                let y_phs = wall_admittance_per_metre(pin, s.radius_m, omega)
                    .map_err(|_| DrivingPointError::Invalid { what: "airway wall pin and frequency" })?;
                y_tube = y_tube + y_phs.conj();
            }
            let n = f64::from(s.count);
            let (z, y) = (z_tube.scale(1.0 / n), y_tube.scale(n));
            // `Z = z / gamma` keeps the matrix independent of the root's sign: both odd entries are even in gamma.
            let gamma = (z * y).sqrt();
            let zc = z * gamma.recip();
            let (ch, sh) = c_cosh_sinh(gamma.scale(s.length_m));
            m = mat_mul(m, [[ch, zc * sh], [sh * zc.recip(), ch]]);
        }
        if m.iter().flatten().any(|v| !(v.re.is_finite() && v.im.is_finite())) {
            return Err(DrivingPointError::Invalid { what: "airway chain matrix left the finite set" });
        }
        Ok(m)
    }

    /// The chain matrix `M` from the valve to the far end at `omega` [rad/s] (module doc), `e^{-iwt}` convention.
    ///
    /// # Errors
    /// A gas-wave or wall refusal at this frequency, or a non-finite matrix.
    pub fn chain_matrix(&self, gas: &GasState, omega: f64) -> Result<[[C64; 2]; 2], DrivingPointError> {
        let m = self.telegraph_chain(gas, omega)?;
        Ok([[m[0][0].conj(), m[0][1].conj()], [m[1][0].conj(), m[1][1].conj()]])
    }

    /// `Z_in = M12 / M22` [Pa s/m^3], the airway seen from the valve with the far end held at a fixed pressure.
    ///
    /// # Errors
    /// As [`Self::chain_matrix`], or a lossless pole (`M22 = 0`).
    pub fn input_impedance(&self, gas: &GasState, omega: f64) -> Result<C64, DrivingPointError> {
        let m = self.chain_matrix(gas, omega)?;
        let z = m[0][1] / m[1][1];
        if !(z.re.is_finite() && z.im.is_finite()) {
            return Err(DrivingPointError::Invalid { what: "airway input impedance is unrepresentable at a lossless pole" });
        }
        Ok(z)
    }

    /// `R = (M12 - Z_s M22) / (M12 + Z_s M22)`, which is `(Z_in - Z_s) / (Z_in + Z_s)` without forming `Z_in`.
    ///
    /// # Errors
    /// As [`Self::chain_matrix`].
    pub fn reflectance(&self, gas: &GasState, omega: f64) -> Result<C64, DrivingPointError> {
        let m = self.chain_matrix(gas, omega)?;
        let zs = self.port_impedance(gas);
        let scaled = m[1][1].scale(zs);
        Ok((m[0][1] - scaled) / (m[0][1] + scaled))
    }

    /// `G = Z_s / (Z_s M22 + M12)`, the wave the valve receives per unit source pressure.
    ///
    /// # Errors
    /// As [`Self::chain_matrix`].
    pub fn source_transfer(&self, gas: &GasState, omega: f64) -> Result<C64, DrivingPointError> {
        let m = self.chain_matrix(gas, omega)?;
        let zs = self.port_impedance(gas);
        Ok(C64::from_re(zs) / (m[1][1].scale(zs) + m[0][1]))
    }
}

/// The source the valve's junction sees on its upstream side for one sample: `p = open_circuit_pa - impedance U`, with
/// `U` the flow into the valve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UpstreamPort {
    /// Upstream pressure with no flow into the valve this sample [Pa].
    pub open_circuit_pa: f64,
    /// Series impedance against the flow into the valve [Pa s/m^3]; zero for an ideal source.
    pub impedance: f64,
}

impl UpstreamPort {
    /// An ideal pressure source: the pressure whatever the flow, which is the valve's input before this seam existed.
    #[must_use]
    pub const fn ideal(pressure_pa: f64) -> Self {
        Self {
            open_circuit_pa: pressure_pa,
            impedance: 0.0,
        }
    }

    /// The upstream pressure at the flow `flow_m3_s` into the valve. A zero impedance returns the open-circuit pressure
    /// itself, bit for bit, whatever the flow.
    #[must_use]
    pub fn pressure_pa(&self, flow_m3_s: f64) -> f64 {
        if self.impedance == 0.0 {
            self.open_circuit_pa
        } else {
            self.open_circuit_pa - self.impedance * flow_m3_s
        }
    }

    /// The same port with a pressure source `source_pa` in series on the valve side (a noise source at the junction).
    /// Zero returns the port unchanged, bit for bit.
    #[must_use]
    pub fn in_series(self, source_pa: f64) -> Self {
        if source_pa == 0.0 {
            self
        } else {
            Self {
                open_circuit_pa: self.open_circuit_pa + source_pa,
                ..self
            }
        }
    }
}

/// A plain FIR whose input history is held twice over, so the last `N` inputs are one contiguous slice.
#[derive(Debug, Clone)]
struct PastFir {
    taps: Vec<f64>,
    buf: Vec<f64>,
    pos: usize,
}

impl PastFir {
    fn new(taps: Vec<f64>) -> Self {
        let n = taps.len();
        Self {
            taps,
            buf: vec![0.0; 2 * n],
            pos: 0,
        }
    }

    /// `sum_{k>=1} h_k x_{n-k}`: the response to the inputs before this sample, in tap order.
    fn past(&self) -> f64 {
        let n = self.taps.len();
        // Oldest first: window[n - j] is the input j samples ago.
        let window = &self.buf[self.pos..self.pos + n];
        let mut acc = 0.0;
        for (&h, &x) in self.taps[1..].iter().zip(window[1..].iter().rev()) {
            acc += h * x;
        }
        acc
    }

    fn push(&mut self, x: f64) {
        let n = self.taps.len();
        self.buf[self.pos] = x;
        self.buf[self.pos + n] = x;
        self.pos = (self.pos + 1) % n;
    }
}

/// One rung of the source-transfer ladder.
#[derive(Debug, Clone, PartialEq)]
pub struct TransferRung {
    /// DFT length of this rung.
    pub n_fft: usize,
    /// Relative energy in the last quarter of the taps.
    pub tail_energy: f64,
    /// `||h_n[..m] - h_{n/2}[..m]|| / ||h_n||` over the first `m = n / 4` taps (`+inf` on the first rung).
    pub tap_change: f64,
}

/// What a realized airway reports about itself.
#[derive(Debug, Clone, PartialEq)]
pub struct AirwayReceipt {
    /// Sections (generations) in the chain.
    pub sections: usize,
    /// `Z_s` [Pa s/m^3].
    pub port_impedance: f64,
    /// `R_dc` [Pa s/m^3].
    pub dc_resistance: f64,
    /// The reflectance's same-sample tap `r_0`.
    pub same_sample_reflection: f64,
    /// The series impedance of every sample's port, `Z_s (1 + r_0) / (1 - r_0)` [Pa s/m^3].
    pub series_impedance: f64,
    /// The reflectance line's ladder receipt.
    pub reflectance: LineRealization,
    /// Converged length of the source-transfer FIR.
    pub source_taps: usize,
    /// Every rung of the source-transfer ladder.
    pub source_rungs: Vec<TransferRung>,
}

/// A realized airway on one sample clock (module doc): the reflectance and source-transfer FIRs and their histories.
#[derive(Debug, Clone)]
pub struct UpstreamAirway {
    reflect: PastFir,
    source: PastFir,
    zs: f64,
    r0: f64,
    g0: f64,
    series_impedance: f64,
    sample_rate_hz: u32,
    /// The source pressure and the known part `A` of this sample's received wave, between `port` and `commit`.
    pending: Option<(f64, f64)>,
    receipt: AirwayReceipt,
}

fn check_ladder(ladder: &LineLadder) -> Result<usize, DrivingPointError> {
    if !(ladder.tail_tol > 0.0 && ladder.tail_tol.is_finite() && ladder.peak_tol_cents > 0.0
        && ladder.peak_tol_cents.is_finite())
    {
        return Err(DrivingPointError::Invalid { what: "line ladder tolerances must be positive and finite" });
    }
    if !(256..=MAX_AIRWAY_TAPS).contains(&ladder.max_taps) || !ladder.max_taps.is_power_of_two() {
        return Err(DrivingPointError::Invalid {
            what: "line ladder max_taps must be a power of two in [256, 131072]",
        });
    }
    Ok(ladder.max_taps)
}

/// The source-transfer FIR, on the mouth-flow line's ladder with an exact DC bin.
fn realize_source_transfer(
    chain: &AirwayChain,
    gas: &GasState,
    sample_rate_hz: u32,
    transit: f64,
    dc: f64,
    ladder: &LineLadder,
) -> Result<(Vec<f64>, Vec<TransferRung>), DrivingPointError> {
    let max_taps = check_ladder(ladder)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let demand = (4.0 * transit).ceil();
    if !demand.is_finite() {
        return Err(DrivingPointError::Invalid { what: "airway transit delay is not finite" });
    }
    let mut n_fft = (demand.min(max_taps as f64) as usize).max(256).next_power_of_two().min(max_taps);
    if max_taps >= 512 {
        n_fft = n_fft.min(max_taps / 2).max(256);
    }
    let tap_tol = ladder.tail_tol.sqrt();
    let mut t_pos: Vec<C64> = Vec::new();
    let mut previous: Vec<f64> = Vec::new();
    let mut rungs: Vec<TransferRung> = Vec::new();
    loop {
        let half = n_fft / 2;
        let mut next = vec![C64::ZERO; half + 1];
        if t_pos.is_empty() {
            for (k, slot) in next.iter_mut().enumerate().skip(1) {
                *slot = chain.source_transfer(gas, core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt))?;
            }
        } else {
            // Even bins of this rung are the bins of the previous one.
            for (k_prev, &value) in t_pos.iter().enumerate() {
                next[2 * k_prev] = value;
            }
            for k in (1..half).step_by(2) {
                next[k] = chain.source_transfer(gas, core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt))?;
            }
        }
        next[0] = C64::from_re(dc);
        t_pos = next;
        let ir = realize_reflectance_ir_from_positive_spectrum(&t_pos, n_fft, transit, true);
        let total: f64 = ir.iter().map(|h| h * h).sum();
        let tail: f64 = ir[3 * (n_fft / 4)..].iter().map(|h| h * h).sum();
        let tail_energy = if total > 0.0 { tail / total } else { 0.0 };
        let tap_change = if previous.is_empty() {
            f64::INFINITY
        } else {
            let diff: f64 = ir[..previous.len()].iter().zip(&previous).map(|(a, b)| (a - b) * (a - b)).sum();
            if total > 0.0 { (diff / total).sqrt() } else { 0.0 }
        };
        rungs.push(TransferRung { n_fft, tail_energy, tap_change });
        if (n_fft as f64) > transit + 2.0 && tail_energy < ladder.tail_tol && tap_change < tap_tol {
            return Ok((ir, rungs));
        }
        if n_fft >= max_taps {
            return Err(DrivingPointError::LadderUnconverged {
                max_taps,
                tail_energy,
                tail_tol: ladder.tail_tol,
                // The transfer ladder gates on the tap change, not on peaks.
                peak_shift_cents: f64::NAN,
                peak_tol_cents: ladder.peak_tol_cents,
            });
        }
        previous = ir[..half].to_vec();
        n_fft *= 2;
    }
}

impl UpstreamAirway {
    /// Realize `chain` filled with `gas` on a `sample_rate_hz` clock (module doc), at rest.
    ///
    /// # Errors
    /// A zero rate; either ladder unconverged at `ladder.max_taps`; a reflectance line that is not passive or not a
    /// plain FIR; a gas-wave or wall refusal at any bin.
    pub fn realize(
        chain: &AirwayChain,
        gas: &GasState,
        sample_rate_hz: u32,
        ladder: &LineLadder,
    ) -> Result<Self, DrivingPointError> {
        if sample_rate_hz == 0 {
            return Err(DrivingPointError::Invalid { what: "an airway needs a positive sample rate" });
        }
        let zs = chain.port_impedance(gas);
        let r_dc = chain.dc_resistance(gas);
        if !(zs.is_finite() && zs > 0.0 && r_dc.is_finite() && r_dc >= 0.0) {
            return Err(DrivingPointError::Invalid { what: "airway port impedance and DC resistance must be finite" });
        }
        let transit = chain.path_length_m() / gas.sound_speed * f64::from(sample_rate_hz);
        let round_trip = 2.0 * transit;
        let samples = round_trip.floor().max(0.0) as usize;
        let frac = (round_trip - samples as f64).clamp(0.0, 1.0 - f64::EPSILON);
        let (line, reflectance) = realize_characteristic_response(
            round_trip,
            LineDc::Explicit((r_dc - zs) / (r_dc + zs)),
            sample_rate_hz,
            zs,
            vec![PoleState::BulkDelay { samples, frac }],
            ladder,
            |omega| chain.reflectance(gas, omega),
        )?;
        let r_taps = line
            .quiescent_fir()
            .ok_or(DrivingPointError::Invalid { what: "the airway reflectance line is not a plain FIR" })?
            .1
            .to_vec();
        let (g_taps, source_rungs) =
            realize_source_transfer(chain, gas, sample_rate_hz, transit, zs / (zs + r_dc), ladder)?;
        let (r0, g0) = (r_taps[0], g_taps[0]);
        if !(r0.abs() < 1.0) {
            return Err(DrivingPointError::Invalid { what: "the airway reflects its same-sample wave whole" });
        }
        let series_impedance = zs * (1.0 + r0) / (1.0 - r0);
        let receipt = AirwayReceipt {
            sections: chain.sections().len(),
            port_impedance: zs,
            dc_resistance: r_dc,
            same_sample_reflection: r0,
            series_impedance,
            reflectance,
            source_taps: g_taps.len(),
            source_rungs,
        };
        Ok(Self {
            reflect: PastFir::new(r_taps),
            source: PastFir::new(g_taps),
            zs,
            r0,
            g0,
            series_impedance,
            sample_rate_hz,
            pending: None,
            receipt,
        })
    }

    /// The Thevenin source the valve sees this sample with the source pressure `source_pa` at the far end (module doc).
    /// Nothing moves until [`Self::commit`]; a second call before it replaces the first.
    ///
    /// # Errors
    /// A non-finite source pressure or port.
    pub fn port(&mut self, source_pa: f64) -> Result<UpstreamPort, DrivingPointError> {
        if !source_pa.is_finite() {
            return Err(DrivingPointError::Invalid { what: "airway source pressure must be finite" });
        }
        let known = self.reflect.past() + (self.g0 * source_pa + self.source.past());
        let open_circuit_pa = 2.0 * known / (1.0 - self.r0);
        if !open_circuit_pa.is_finite() {
            return Err(DrivingPointError::Invalid { what: "airway port left the finite set" });
        }
        self.pending = Some((source_pa, known));
        Ok(UpstreamPort {
            open_circuit_pa,
            impedance: self.series_impedance,
        })
    }

    /// Close this sample at the flow `flow_m3_s` into the valve: send the wave it implies down the airway and return the
    /// pressure at the valve [Pa].
    ///
    /// # Errors
    /// No [`Self::port`] this sample, or a non-finite flow or wave; a refusal leaves the airway as it was.
    pub fn commit(&mut self, flow_m3_s: f64) -> Result<f64, DrivingPointError> {
        let (source_pa, known) =
            self.pending.ok_or(DrivingPointError::Invalid { what: "an airway commit needs this sample's port" })?;
        let b = (known - self.zs * flow_m3_s) / (1.0 - self.r0);
        let p = (1.0 + self.r0) * b + known;
        if !(flow_m3_s.is_finite() && b.is_finite() && p.is_finite()) {
            return Err(DrivingPointError::Invalid { what: "airway wave left the finite set" });
        }
        self.pending = None;
        self.reflect.push(b);
        self.source.push(source_pa);
        Ok(p)
    }

    /// `Z_s` [Pa s/m^3], the reference of the port waves.
    #[must_use]
    pub const fn port_impedance(&self) -> f64 {
        self.zs
    }

    /// The clock the FIRs were realized on [Hz].
    #[must_use]
    pub const fn sample_rate_hz(&self) -> u32 {
        self.sample_rate_hz
    }

    /// The reflectance taps, `r_0` first.
    #[must_use]
    pub fn reflectance_taps(&self) -> &[f64] {
        &self.reflect.taps
    }

    /// The source-transfer taps, `g_0` first.
    #[must_use]
    pub fn source_taps(&self) -> &[f64] {
        &self.source.taps
    }

    /// The realization receipt.
    #[must_use]
    pub const fn receipt(&self) -> &AirwayReceipt {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_duct::{Duct, Segment, Termination, input_impedance_load, input_impedance_wall};
    use fs_material::gas::GasSpec;

    const RATE: u32 = 48_000;

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    /// The voice cards' `tract_wall 15 1.48e+06 1.6e4`.
    fn wall() -> WallPin {
        WallPin {
            surface_density: 15.0,
            stiffness_per_area: 1.48e6,
            resistance: 1.6e4,
        }
    }

    fn tube(radius_m: f64, length_m: f64, count: u32, wall: Option<WallPin>) -> AirwaySection {
        AirwaySection { length_m, radius_m, count, wall }
    }

    fn ladder() -> LineLadder {
        LineLadder {
            max_taps: MAX_AIRWAY_TAPS,
            ..LineLadder::default()
        }
    }

    fn rel(a: C64, b: C64) -> f64 {
        (a - b).abs() / b.abs().max(1.0e-300)
    }

    /// The DTFT `sum h_n e^{+i w n dt}` of a realized FIR (`e^{-iwt}` convention).
    fn dtft(taps: &[f64], omega: f64) -> C64 {
        let dt = 1.0 / f64::from(RATE);
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for (n, h) in taps.iter().enumerate() {
            let phase = omega * n as f64 * dt;
            re += h * phase.cos();
            im += h * phase.sin();
        }
        C64::new(re, im)
    }

    /// One cylinder with the far end held at a fixed pressure. Lossless and rigid, `Z_in = -i (rho c / S) tan(w L / c)`
    /// in closed form; Bessel with the cards' wall, it is fs-duct's own transfer-matrix duct with an ideal open end
    /// (`p = 0` there), computed by a different route (Helmholtz basis solutions, not the telegraph matrix here).
    #[test]
    fn ua_001_a_uniform_cylinder_is_its_closed_form_and_fs_ducts_duct() {
        let gas = air();
        let (radius, length) = (0.009, 0.12);
        let lossless = AirwayChain::new(vec![tube(radius, length, 1, None)], LossModel::Lossless).expect("chain");
        let lossy = AirwayChain::new(vec![tube(radius, length, 1, Some(wall()))], LossModel::Bessel).expect("chain");
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        };
        let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let (mut worst_closed, mut worst_duct) = (0.0_f64, 0.0_f64);
        for f in [50.0, 333.0, 1000.0, 2717.0, 5119.0, 9000.0] {
            let omega = core::f64::consts::TAU * f;
            let closed = C64::new(0.0, -zc * (omega * length / gas.sound_speed).tan());
            worst_closed = worst_closed.max(rel(lossless.input_impedance(&gas, omega).expect("z"), closed));
            let reference = input_impedance_wall(&duct, &gas, omega, LossModel::Bessel, Termination::IdealOpen, Some(&wall()))
                .expect("fs-duct")
                .impedance;
            worst_duct = worst_duct.max(rel(lossy.input_impedance(&gas, omega).expect("z"), reference));
        }
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ua-001-uniform-cylinder\",\"worst_rel_closed_form\":{worst_closed:e},\"worst_rel_fs_duct\":{worst_duct:e}}}"
        );
        assert!(worst_closed < 1.0e-10, "lossless cylinder misses -i Zc tan kL by {worst_closed:e}");
        assert!(worst_duct < 1.0e-9, "Bessel cylinder with wall misses fs-duct by {worst_duct:e}");
    }

    /// A generation of two tubes is two tubes in parallel. Input impedance: each daughter alone (fs-duct, `p = 0` at its
    /// far end), the two in parallel as the parent's load (fs-duct's loaded duct), against the chain with `count = 2`.
    /// Open-circuit source gain: the network's nodes solved by hand with one-tube matrices (the junction pressure is
    /// common, the parent carries the sum of the daughters' flows), against `1 / M22`; and `G` from those two.
    #[test]
    fn ua_002_a_count_of_two_is_two_explicit_tubes_in_parallel() {
        let gas = air();
        let parent = tube(0.009, 0.12, 1, Some(wall()));
        let daughter = tube(0.006_105, 0.0478, 1, Some(wall()));
        let chain = AirwayChain::new(vec![parent, AirwaySection { count: 2, ..daughter }], LossModel::Bessel).expect("chain");
        let one = |s: AirwaySection| AirwayChain::new(vec![s], LossModel::Bessel).expect("chain");
        let duct = |s: &AirwaySection| Duct {
            segments: vec![Segment::Cylinder { radius: s.radius_m, length: s.length_m }],
        };
        let zs = chain.port_impedance(&gas);
        let (mut worst_z, mut worst_oc, mut worst_g) = (0.0_f64, 0.0_f64, 0.0_f64);
        for f in [20.0, 180.0, 640.0, 1400.0, 2100.0, 4400.0, 11000.0] {
            let omega = core::f64::consts::TAU * f;
            let z_daughter =
                input_impedance_wall(&duct(&daughter), &gas, omega, LossModel::Bessel, Termination::IdealOpen, Some(&wall()))
                    .expect("daughter")
                    .impedance;
            let z_load = (z_daughter.recip() + z_daughter.recip()).recip();
            let z_explicit = input_impedance_load(&duct(&parent), &gas, omega, LossModel::Bessel, z_load, Some(&wall()))
                .expect("loaded parent")
                .impedance;
            worst_z = worst_z.max(rel(chain.input_impedance(&gas, omega).expect("z"), z_explicit));
            let a = one(parent).chain_matrix(&gas, omega).expect("parent");
            let d = one(daughter).chain_matrix(&gas, omega).expect("daughter");
            // Daughter i: U_far = (p_J - d11 P) / d12 and U_J = d21 P + d22 U_far; the parent at no valve flow:
            // a21 p_J + a22 (U_J1 + U_J2) = 0, so p_J = 2 a22 P / (d12 a21 + 2 a22 d22) and p_valve = p_J / a22.
            let two = C64::from_re(2.0);
            let p_junction = two * a[1][1] / (d[0][1] * a[1][0] + two * a[1][1] * d[1][1]);
            let u_far = (p_junction - d[0][0]) / d[0][1];
            let u_junction = d[1][0] + d[1][1] * u_far;
            let u_parent = u_junction + u_junction;
            let (parent_shunt, parent_series) = (a[1][0] * p_junction, a[1][1] * u_parent);
            let residual = (parent_shunt + parent_series).abs();
            let node_scale = parent_shunt.abs().max(parent_series.abs()).max(1.0e-300);
            assert!(residual < 1.0e-9 * node_scale, "node solve residual {residual:e} of {node_scale:e}");
            let open_circuit = a[0][0] * p_junction + a[0][1] * u_parent;
            let m = chain.chain_matrix(&gas, omega).expect("chain");
            worst_oc = worst_oc.max(rel(m[1][1].recip(), open_circuit));
            let g_explicit = open_circuit * C64::from_re(zs) / (C64::from_re(zs) + z_explicit);
            worst_g = worst_g.max(rel(chain.source_transfer(&gas, omega).expect("g"), g_explicit));
        }
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ua-002-two-generation-network\",\"worst_rel_impedance\":{worst_z:e},\"worst_rel_open_circuit\":{worst_oc:e},\"worst_rel_source_transfer\":{worst_g:e}}}"
        );
        assert!(worst_z < 1.0e-9, "count 2 misses two tubes in parallel by {worst_z:e}");
        assert!(worst_oc < 1.0e-9, "open-circuit gain misses the node solve by {worst_oc:e}");
        assert!(worst_g < 1.0e-9, "source transfer misses the network by {worst_g:e}");
    }

    fn weibel_airway() -> (AirwayChain, UpstreamAirway) {
        let chain = AirwayChain::weibel_saab_2021_b(17, Some(wall())).expect("tree");
        let airway = UpstreamAirway::realize(&chain, &air(), RATE, &ladder()).expect("airway realizes");
        (chain, airway)
    }

    /// Run `steps` samples at a held source pressure and flow; return the last port and pressure.
    fn hold(airway: &mut UpstreamAirway, source_pa: f64, flow: f64, steps: usize) -> (UpstreamPort, f64) {
        let mut last = (UpstreamPort::ideal(0.0), 0.0);
        for _ in 0..steps {
            let port = airway.port(source_pa).expect("port");
            let p = airway.commit(flow).expect("commit");
            last = (port, p);
        }
        last
    }

    /// DC, on Weibel's tree to 17 generations with the cards' wall. After both FIRs' lengths of held input (the loop
    /// that a held flow closes at the valve has rung down long before), no flow reads the source pressure and a steady
    /// flow reads `P_s - R_dc U` with the Poiseuille sum; the Thevenin form the valve solves with and the pressure the
    /// airway returns agree. The realized FIRs reproduce `R` and `G` between their bins, below the first bin included.
    #[test]
    fn ua_003_a_steady_flow_reads_the_source_pressure_less_the_poiseuille_drop() {
        let gas = air();
        let (chain, mut airway) = weibel_airway();
        let settle = airway.reflectance_taps().len() + airway.source_taps().len() + 8;
        let r_dc = chain.dc_resistance(&gas);
        let (p_lung, flow) = (800.0, 2.0e-4);
        let (port, p_open) = hold(&mut airway, p_lung, 0.0, settle);
        let open_miss = (p_open - p_lung).abs() / p_lung;
        let thevenin_open = (port.open_circuit_pa - p_lung).abs() / p_lung;
        let (port, p_flow) = hold(&mut airway, p_lung, flow, settle);
        let drop_miss = (p_flow - (p_lung - r_dc * flow)).abs() / p_lung;
        let thevenin_flow = (port.pressure_pa(flow) - p_flow).abs() / p_lung;
        let (_, p_zero) = hold(&mut airway, 0.0, flow, settle);
        let zero_miss = (p_zero + r_dc * flow).abs() / (r_dc * flow);
        let mut worst_fir = 0.0_f64;
        for f in [1.0, 7.0, 310.0, 1013.0, 2207.0, 5119.0] {
            let omega = core::f64::consts::TAU * f;
            worst_fir = worst_fir
                .max((dtft(airway.reflectance_taps(), omega) - chain.reflectance(&gas, omega).expect("r")).abs())
                .max((dtft(airway.source_taps(), omega) - chain.source_transfer(&gas, omega).expect("g")).abs());
        }
        let receipt = airway.receipt();
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ua-003-steady-flow\",\"r_dc\":{r_dc:e},\"zs\":{:e},\"r0\":{:e},\"reflectance_taps\":{},\"source_taps\":{},\"open_miss\":{open_miss:e},\"drop_miss\":{drop_miss:e},\"zero_source_miss\":{zero_miss:e},\"thevenin_open\":{thevenin_open:e},\"thevenin_flow\":{thevenin_flow:e},\"worst_fir_abs\":{worst_fir:e}}}",
            receipt.port_impedance,
            receipt.same_sample_reflection,
            receipt.reflectance.n_fft,
            receipt.source_taps
        );
        assert!(open_miss < 1.0e-9 && thevenin_open < 1.0e-9, "no flow must read P_s ({open_miss:e}, {thevenin_open:e})");
        assert!(drop_miss < 1.0e-9, "a steady flow must read P_s - R_dc U ({drop_miss:e})");
        assert!(zero_miss < 1.0e-6, "no source must read -R_dc U ({zero_miss:e})");
        assert!(thevenin_flow < 1.0e-12, "the port and the commit disagree ({thevenin_flow:e})");
        assert!(worst_fir < 1.0e-3, "the FIRs depart from R and G by {worst_fir:e}");
    }

    /// One second of seeded white flow into `airway` with no source pressure: the energy absorbed, `sum p (-U) dt` from
    /// rest, its scale `sum |p U| dt`, and the lowest running ratio of the two.
    fn white_flow_work(airway: &mut UpstreamAirway) -> (f64, f64, f64) {
        let dt = 1.0 / f64::from(RATE);
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut uniform = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        };
        let (mut absorbed, mut scale, mut worst) = (0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..RATE {
            let flow = 4.0e-4 * uniform();
            airway.port(0.0).expect("port");
            let p = airway.commit(flow).expect("commit");
            absorbed += -p * flow * dt;
            scale += (p * flow).abs() * dt;
            worst = worst.min(absorbed / scale.max(1.0e-300));
        }
        (absorbed, scale, worst)
    }

    /// Passivity: with no source pressure, the energy the valve's flow puts into the airway, `sum p (-U) dt` from rest,
    /// never falls below zero under a seeded white flow, at every sample.
    #[test]
    fn ua_004_the_airway_never_returns_more_energy_than_it_was_given() {
        let (_, mut airway) = weibel_airway();
        let (absorbed, scale, worst) = white_flow_work(&mut airway);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ua-004-airway-passivity\",\"absorbed_j\":{absorbed:e},\"scale_j\":{scale:e},\"worst_running_ratio\":{worst:e}}}"
        );
        assert!(absorbed > 0.0, "a lossy airway absorbs ({absorbed:e})");
        assert!(worst > -1.0e-9, "the airway returned more than it was given (running ratio {worst:e})");
    }

    /// The ideal port is the plain pressure, bit for bit, at every flow; a zero series source changes nothing.
    #[test]
    fn ua_005_the_ideal_port_is_its_pressure_bit_for_bit() {
        for p in [0.0, -0.0, 812.5, -3.25e-7, 1153.49] {
            let port = UpstreamPort::ideal(p);
            for u in [-1.0e-3, -0.0, 0.0, 2.5e-4, f64::MIN_POSITIVE] {
                assert_eq!(port.pressure_pa(u).to_bits(), f64::to_bits(p), "ideal({p}) at {u}");
            }
            assert_eq!(port.in_series(0.0), port);
            assert_eq!(port.in_series(0.0).open_circuit_pa.to_bits(), f64::to_bits(p));
        }
    }

    /// A manufactured duct ahead of the tree, NOT a measured tract: a 40 mm by 3 cm^2 cavity, a 100 mm by 2 cm^2 tube
    /// and a 5 mm by 0.5 cm^2 opening, each from its area, all with the cards' wall.
    fn fixture_duct() -> AirwayChain {
        let w = Some(wall());
        let sections = [(0.040, 3.0e-4), (0.100, 2.0e-4), (0.005, 0.5e-4)]
            .iter()
            .map(|&(length, area)| AirwaySection::from_area(length, area, 1, w).expect("section"))
            .collect();
        AirwayChain::new(sections, LossModel::Bessel).expect("duct")
    }

    /// A duct ahead of the tree is the duct loaded by the tree. Its input impedance against fs-duct's loaded duct (the
    /// same three cylinders, the tree's own `Z_in` as the load, Helmholtz basis solutions rather than this module's
    /// telegraph matrix); its sections, `Z_s`, path length and DC resistance against the two chains'. An area section
    /// has its area; mixed loss models and bad areas refuse.
    #[test]
    fn ua_006_a_duct_ahead_of_the_tree_is_the_duct_loaded_by_the_tree() {
        let gas = air();
        let duct = fixture_duct();
        let tree = AirwayChain::weibel_saab_2021_b(17, Some(wall())).expect("tree");
        let joined = duct.clone().ahead_of(&tree).expect("joined");
        let cylinders = Duct {
            segments: duct
                .sections()
                .iter()
                .map(|s| Segment::Cylinder { radius: s.radius_m, length: s.length_m })
                .collect(),
        };
        let mut worst_z = 0.0_f64;
        for f in [20.0, 180.0, 640.0, 1400.0, 2100.0, 4400.0, 11000.0] {
            let omega = core::f64::consts::TAU * f;
            let z_tree = tree.input_impedance(&gas, omega).expect("tree z");
            let reference = input_impedance_load(&cylinders, &gas, omega, LossModel::Bessel, z_tree, Some(&wall()))
                .expect("loaded duct")
                .impedance;
            worst_z = worst_z.max(rel(joined.input_impedance(&gas, omega).expect("z"), reference));
        }
        let n = duct.sections().len();
        assert_eq!(&joined.sections()[..n], duct.sections());
        assert_eq!(&joined.sections()[n..], tree.sections());
        assert_eq!(joined.port_impedance(&gas).to_bits(), duct.port_impedance(&gas).to_bits());
        let path_miss = (joined.path_length_m() - (duct.path_length_m() + tree.path_length_m())).abs();
        let r_sum = duct.dc_resistance(&gas) + tree.dc_resistance(&gas);
        let dc_miss = (joined.dc_resistance(&gas) - r_sum).abs() / r_sum;
        let area_miss = duct
            .sections()
            .iter()
            .zip([3.0e-4, 2.0e-4, 0.5e-4])
            .map(|(s, area)| (core::f64::consts::PI * s.radius_m * s.radius_m - area).abs() / area)
            .fold(0.0_f64, f64::max);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ua-006-duct-ahead-of-tree\",\"worst_rel_loaded_duct\":{worst_z:e},\"path_miss_m\":{path_miss:e},\"dc_rel_miss\":{dc_miss:e},\"area_rel_miss\":{area_miss:e}}}"
        );
        assert!(worst_z < 1.0e-9, "the joined chain misses the duct loaded by the tree by {worst_z:e}");
        assert!(path_miss < 1.0e-15 && dc_miss < 1.0e-13, "path {path_miss:e} m, DC {dc_miss:e}");
        assert!(area_miss < 4.0 * f64::EPSILON, "an area section misses its area by {area_miss:e}");
        let lossless = AirwayChain::new(duct.sections().to_vec(), LossModel::Lossless).expect("lossless duct");
        assert!(lossless.ahead_of(&tree).is_err(), "mixed loss models must refuse");
        let refused = [(0.01, 0.0, 1), (0.01, -1.0e-4, 1), (0.0, 1.0e-4, 1), (0.01, f64::NAN, 1), (0.01, 1.0e-4, 0)];
        for (length, area, count) in refused {
            assert!(AirwaySection::from_area(length, area, count, None).is_err(), "({length}, {area}, {count})");
        }
    }

    /// The port behind a duct ahead of the tree is the same port: realized on the same ladders, a steady flow reads
    /// `P_s - R_dc U` with the joined Poiseuille sum, and a seeded white flow never draws energy out of it.
    #[test]
    fn ua_007_a_duct_ahead_of_the_tree_realizes_as_a_passive_port() {
        let gas = air();
        let chain = fixture_duct()
            .ahead_of(&AirwayChain::weibel_saab_2021_b(17, Some(wall())).expect("tree"))
            .expect("joined");
        let mut airway = UpstreamAirway::realize(&chain, &gas, RATE, &ladder()).expect("joined airway realizes");
        let settle = airway.reflectance_taps().len() + airway.source_taps().len() + 8;
        let r_dc = chain.dc_resistance(&gas);
        let (p_lung, flow) = (800.0, 2.0e-4);
        let (_, p_open) = hold(&mut airway, p_lung, 0.0, settle);
        let open_miss = (p_open - p_lung).abs() / p_lung;
        let (port, p_flow) = hold(&mut airway, p_lung, flow, settle);
        let drop_miss = (p_flow - (p_lung - r_dc * flow)).abs() / p_lung;
        let thevenin_flow = (port.pressure_pa(flow) - p_flow).abs() / p_lung;
        let mut fresh = UpstreamAirway::realize(&chain, &gas, RATE, &ladder()).expect("joined airway realizes");
        let (absorbed, scale, worst) = white_flow_work(&mut fresh);
        let receipt = airway.receipt();
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ua-007-duct-ahead-of-tree-port\",\"sections\":{},\"r_dc\":{r_dc:e},\"zs\":{:e},\"r0\":{:e},\"reflectance_taps\":{},\"source_taps\":{},\"open_miss\":{open_miss:e},\"drop_miss\":{drop_miss:e},\"thevenin_flow\":{thevenin_flow:e},\"absorbed_j\":{absorbed:e},\"scale_j\":{scale:e},\"worst_running_ratio\":{worst:e}}}",
            receipt.sections,
            receipt.port_impedance,
            receipt.same_sample_reflection,
            receipt.reflectance.n_fft,
            receipt.source_taps
        );
        assert_eq!(receipt.sections, 3 + 17);
        assert!(open_miss < 1.0e-9, "no flow must read P_s ({open_miss:e})");
        assert!(drop_miss < 1.0e-9, "a steady flow must read P_s - R_dc U ({drop_miss:e})");
        assert!(thevenin_flow < 1.0e-12, "the port and the commit disagree ({thevenin_flow:e})");
        assert!(absorbed > 0.0, "a lossy airway absorbs ({absorbed:e})");
        assert!(worst > -1.0e-9, "the airway returned more than it was given (running ratio {worst:e})");
    }
}
