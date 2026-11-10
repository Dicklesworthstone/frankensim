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
//! Reflection, input impedance and source transfer use the native bounded
//! [`fs_duct::compose::TwoPort`] wave blocks and Redheffer composition. The
//! default source boundary is `a_R + b_R = P_s`, so `R = S11 - S12 S21/(1+S22)` and
//! `G = S12/(1+S22)`. A finite passive source impedance `Z_L(w)` instead imposes
//! `p_R + Z_L U_source = P_s`, with `U_source` INTO the airway, just as
//! [`fs_duct::network::PortExcitation::Pressure`] does. With
//! `r_L = (Z_L-Z_s)/(Z_L+Z_s)` and `g_L = Z_s/(Z_L+Z_s)`, its waves obey
//! `a_R = r_L b_R + g_L P_s`. Thus `R = S11 + S12 r_L S21/(1-r_L S22)` and
//! `G = S12 g_L/(1-r_L S22)`: changing the source changes both the returning
//! reflection and the pressure delivered to the valve. No growing exponential or input-impedance pole enters
//! those wave observables. [`AirwayChain::chain_matrix`] is an explicit ABCD
//! diagnostic only and may refuse when a strongly attenuating line's matrix
//! cannot be represented, even though its reflected waves remain finite.
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
//! [`UpstreamAirway::realize_with_source_impedance`] realizes an explicitly supplied
//! causal source law, caching each frequency's paired `(R,G)` so both FIRs use the
//! identical admitted boundary. Its supplied finite DC resistance adds to the
//! airway's Poiseuille resistance in both exact DC bins. Source impedance is
//! evaluated only during the bounded, fallible preprocessing; the sample loop
//! is unchanged. Nonnegative resistance at the sampled frequencies and the
//! realized reflection's passivity check do not certify an arbitrary law's causality.
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
//! jet loss. The default far end is an ideal pressure source. An explicitly supplied source impedance may include
//! resistance, inertance or compliant relaxation, but its causal law belongs to the caller; no alveolar parameters
//! are inferred. [`UpstreamAirway::realize_with_source_admittance`] also accepts sealed compliance and resonant zeros
//! through U_source=Y_L(P_s-p_R), with exact R(0)=1, G(0)=0 when Y_L(0)=0. This is a zero-state wave response;
//! initialized reservoir storage needs an explicit owner. [`StoredPressureSource`] provides that owner through
//! the original R-L-C storage and native characteristic-network boundary; changing volume remains outside this
//! fixed-coefficient linear model.

/// Initialized source storage at the existing valve or characteristic-network port.
pub mod stored;
pub use stored::{StoredPressureSource, StoredSourceFrame, StoredSourceNetworkFrame};

use crate::driving_point::{
    DrivingPointError, LineDc, LineLadder, LineRealization, PoleState, realize_characteristic_response,
    realize_reflectance_ir_from_positive_spectrum,
};
use fs_duct::{LossModel, SegmentWave, segment_wave};
use fs_duct::compose::{ComposeOptions, TwoPort, compose_chain, compose_chain_reflection};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;
use fs_phs::{WallPin, wall_admittance_per_metre};
use std::collections::BTreeMap;

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

/// Admit a finite passive Thevenin source in the native pressure/volume-flow convention.
fn admit_source_impedance(impedance: C64) -> Result<(), DrivingPointError> {
    if !impedance.re.is_finite() || !impedance.im.is_finite() || impedance.re < 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "airway source impedance must be finite and passive (Re Z >= 0)",
        });
    }
    Ok(())
}

/// Right-port reflection and incident pressure wave per unit open-circuit source pressure.
/// Scale before adding: finite source impedances near f64::MAX still have bounded waves.
fn source_boundary(impedance: C64, reference: f64) -> Result<(C64, C64), DrivingPointError> {
    admit_source_impedance(impedance)?;
    if !reference.is_finite() || reference <= 0.0 {
        return Err(DrivingPointError::Invalid { what: "airway source needs a finite positive wave reference" });
    }
    let scale = impedance.re.max(impedance.im.abs()).max(reference);
    let z = C64::new(impedance.re / scale, impedance.im / scale);
    let r = C64::from_re(reference / scale);
    let denominator = z + r;
    Ok(((z - r) / denominator, r / denominator))
}

/// Admittance form of the same source boundary: U_source = Y_L (P_s - p_R).
/// Y_L=0 is sealed, without forming an infinite impedance. Choose the scaling
/// chart before multiplying so both finite extremes keep bounded wave data.
fn source_admittance_boundary(
    admittance: C64,
    reference: f64,
) -> Result<(C64, C64), DrivingPointError> {
    if !admittance.re.is_finite() || !admittance.im.is_finite() || admittance.re < 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "airway source admittance must be finite and passive (Re Y >= 0)",
        });
    }
    if !reference.is_finite() || reference <= 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "airway source needs a finite positive wave reference",
        });
    }
    let (y, reference_y) = if reference < 1.0 {
        let product = admittance.scale(reference);
        let scale = product.re.max(product.im.abs()).max(1.0);
        (
            C64::new(product.re / scale, product.im / scale),
            1.0 / scale,
        )
    } else {
        let reference_y = 1.0 / reference;
        let scale = admittance.re.max(admittance.im.abs()).max(reference_y);
        (
            C64::new(admittance.re / scale, admittance.im / scale),
            reference_y / scale,
        )
    };
    let reference_y = C64::from_re(reference_y);
    let denominator = reference_y + y;
    Ok(((reference_y - y) / denominator, y / denominator))
}

fn loaded_source_transfer(scattering: TwoPort, impedance: C64) -> Result<C64, DrivingPointError> {
    loaded_source_waves(
        scattering,
        source_boundary(impedance, scattering.reference_impedance)?,
    )
    .map(|response| response.1)
}

/// Contract the authored boundary with the native bounded two-port once.
fn loaded_source_waves(
    scattering: TwoPort,
    boundary: (C64, C64),
) -> Result<(C64, C64), DrivingPointError> {
    let (reflection, drive) = boundary;
    let product = reflection * scattering.scattering[3];
    let denominator = C64::ONE - product;
    let condition = (1.0 + product.abs()) / denominator.abs();
    let gain = scattering.scattering[1] * drive / denominator;
    let returned = scattering.scattering[0]
        + scattering.scattering[1] * reflection * (scattering.scattering[2] / denominator);
    if !condition.is_finite()
        || condition > ComposeOptions::default().maximum_condition
        || !gain.re.is_finite()
        || !gain.im.is_finite()
        || !returned.re.is_finite()
        || !returned.im.is_finite()
    {
        return Err(fs_duct::DuctError::Singular.into());
    }
    Ok((returned, gain))
}

/// Static chain M=[[1,R_dc],[0,1]], directly in the authored admittance chart.
/// This avoids both 1/Y overflow at a sealed limit and Y*(R_dc+Z_s) overflow.
fn admittance_dc_waves(reference: f64, resistance: f64, conductance: f64) -> (f64, f64) {
    if conductance == 0.0 {
        return (1.0, 0.0);
    }
    let scale = reference.max(resistance);
    let r = resistance / scale;
    let z = reference / scale;
    let (one, y) = if conductance <= 1.0 || scale <= 1.0 {
        let product = conductance * scale;
        let normalization = product.max(1.0);
        (1.0 / normalization, product / normalization)
    } else {
        ((1.0 / conductance) / scale, 1.0)
    };
    let denominator = one + y * (r + z);
    ((one + y * (r - z)) / denominator, y * z / denominator)
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
        let mut m = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
        for s in &self.sections {
            let (z_tube, y_tube, _) = self.tube_telegraph(s, gas, omega)?;
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

    // One actual tube's constitutive law, shared by the bounded wave path and
    // the diagnostic matrix. Tube count changes volume flow, never its radius
    // or the viscothermal and locally reacting wall physics.
    fn tube_telegraph(&self, section: &AirwaySection, gas: &GasState, omega: f64)
        -> Result<(C64, C64, f64), DrivingPointError>
    {
        let wave = segment_wave(gas, section.radius_m, omega, self.loss)?;
        let minus_i = C64::new(0.0, -1.0);
        let z = (minus_i * wave.wavenumber * wave.characteristic_impedance).conj();
        let mut y = (minus_i * wave.wavenumber * wave.characteristic_impedance.recip()).conj();
        if let Some(pin) = &section.wall {
            let wall = wall_admittance_per_metre(pin, section.radius_m, omega)
                .map_err(|_| DrivingPointError::Invalid { what: "airway wall pin and frequency" })?;
            y = y + wall.conj();
        }
        Ok((z, y, wave.shear_number))
    }

    fn scattering_chain(&self, gas: &GasState, omega: f64) -> Result<TwoPort, DrivingPointError> {
        let reference = self.port_impedance(gas);
        let mut total = TwoPort::through(reference)?;
        for section in &self.sections {
            let (z, y, shear_number) = self.tube_telegraph(section, gas, omega)?;
            let gamma = (z * y).sqrt();
            let tube_impedance = (z / gamma).conj();
            let tube_area = core::f64::consts::PI * section.radius_m * section.radius_m;
            let count = f64::from(section.count);
            // sqrt gives Re(gamma)>=0; conjugating the telegraph convention
            // gives exp(i k x), k=Im(gamma)+i Re(gamma), with bounded decay.
            let wave = SegmentWave {
                wavenumber: C64::new(gamma.im, gamma.re),
                characteristic_impedance: tube_impedance.scale(1.0 / count),
                specific_impedance: tube_impedance.scale(tube_area),
                shear_number,
            };
            let block = TwoPort::from_uniform_wave(
                &wave, tube_area * count, section.length_m, reference,
            )?;
            total = total.star(block, ComposeOptions::default().maximum_condition)?.0;
        }
        Ok(total)
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
    /// A gas-wave/wall or native scattering-conditioning refusal, or a lossless
    /// input-impedance pole. ABCD growth alone does not refuse this observable.
    pub fn input_impedance(&self, gas: &GasState, omega: f64) -> Result<C64, DrivingPointError> {
        self.input_impedance_with_source_impedance(gas, omega, C64::ZERO)
    }

    /// Valve-side impedance with the source's open-circuit pressure set to zero and its finite
    /// impedance still connected. `source_impedance` is in Pa s/m^3 at this `omega`, in `exp(-iwt)`.
    ///
    /// # Errors
    /// As [`Self::input_impedance`], or a non-finite/negative-real source impedance.
    pub fn input_impedance_with_source_impedance(
        &self, gas: &GasState, omega: f64, source_impedance: C64,
    ) -> Result<C64, DrivingPointError> {
        admit_source_impedance(source_impedance)?;
        let scattering = self.scattering_chain(gas, omega)?;
        Ok(compose_chain(&[scattering], Some(source_impedance), ComposeOptions::default())?.z_in)
    }

    /// `R = (M12 - Z_s M22) / (M12 + Z_s M22)`, which is `(Z_in - Z_s) / (Z_in + Z_s)` without forming `Z_in`.
    ///
    /// # Errors
    /// A gas-wave/wall or native scattering-conditioning refusal. Neither
    /// ABCD growth nor an input-impedance pole refuses a finite reflection.
    pub fn reflectance(&self, gas: &GasState, omega: f64) -> Result<C64, DrivingPointError> {
        self.reflectance_with_source_impedance(gas, omega, C64::ZERO)
    }

    /// Reflected pressure wave per incident valve wave with the far source's impedance connected.
    /// Setting the source pressure to zero does not remove its resistance or reactance.
    ///
    /// # Errors
    /// As [`Self::reflectance`], or a non-finite/negative-real source impedance.
    pub fn reflectance_with_source_impedance(
        &self, gas: &GasState, omega: f64, source_impedance: C64,
    ) -> Result<C64, DrivingPointError> {
        admit_source_impedance(source_impedance)?;
        let scattering = self.scattering_chain(gas, omega)?;
        Ok(compose_chain_reflection(&[scattering], Some(source_impedance), ComposeOptions::default())?.reflection)
    }

    /// `G = Z_s / (Z_s M22 + M12)`, the wave the valve receives per unit source pressure.
    ///
    /// # Errors
    /// A gas-wave/wall or native scattering-conditioning refusal, or an
    /// unrepresentable source transfer. An opaque line transmits zero.
    pub fn source_transfer(&self, gas: &GasState, omega: f64) -> Result<C64, DrivingPointError> {
        self.source_transfer_with_source_impedance(gas, omega, C64::ZERO)
    }

    /// Wave received at the valve per unit open-circuit far-source pressure, including loading by
    /// its finite passive impedance. A source matched to a uniform lossless line launches half
    /// its open-circuit pressure; an ideal source (`Z_L = 0`) launches the whole pressure.
    ///
    /// # Errors
    /// As [`Self::source_transfer`], or a non-finite/negative-real source impedance.
    pub fn source_transfer_with_source_impedance(
        &self, gas: &GasState, omega: f64, source_impedance: C64,
    ) -> Result<C64, DrivingPointError> {
        admit_source_impedance(source_impedance)?;
        loaded_source_transfer(self.scattering_chain(gas, omega)?, source_impedance)
    }

    /// Valve-side impedance when the source is specified as a passive admittance
    /// [m^3/(Pa s)] in exp(-iwt). Zero admittance is a sealed termination.
    ///
    /// # Errors
    /// Non-finite/negative-real admittance, native gas/wall/conditioning refusal,
    /// or an input-impedance pole. Use reflectance at a pressure antinode.
    pub fn input_impedance_with_source_admittance(
        &self,
        gas: &GasState,
        omega: f64,
        source_admittance: C64,
    ) -> Result<C64, DrivingPointError> {
        let reflection = self.reflectance_with_source_admittance(gas, omega, source_admittance)?;
        let denominator = C64::ONE - reflection;
        let condition = (1.0 + reflection.abs()) / denominator.abs();
        let impedance = ((C64::ONE + reflection) / denominator).scale(self.port_impedance(gas));
        if !condition.is_finite()
            || condition > ComposeOptions::default().maximum_condition
            || !impedance.re.is_finite()
            || !impedance.im.is_finite()
        {
            return Err(fs_duct::DuctError::Singular.into());
        }
        Ok(impedance)
    }

    /// Reflection for U_source = Y_L (P_s - p_R), including Y_L=0.
    /// Source pressure zero retains the original admittance. No reciprocal of
    /// Y_L is formed, so a resonant zero or sealed compliance limit is finite.
    ///
    /// # Errors
    /// Non-finite/negative-real admittance or native gas/wall/conditioning refusal.
    pub fn reflectance_with_source_admittance(
        &self,
        gas: &GasState,
        omega: f64,
        source_admittance: C64,
    ) -> Result<C64, DrivingPointError> {
        self.waves_with_source_admittance(gas, omega, source_admittance)
            .map(|response| response.0)
    }

    /// Valve incident wave per unit open-circuit source pressure for the
    /// authored admittance boundary. A sealed source launches no pressure wave;
    /// a compliance Y_L=-i*omega*C retains its frequency-dependent drive.
    ///
    /// # Errors
    /// As [`Self::reflectance_with_source_admittance`].
    pub fn source_transfer_with_source_admittance(
        &self,
        gas: &GasState,
        omega: f64,
        source_admittance: C64,
    ) -> Result<C64, DrivingPointError> {
        self.waves_with_source_admittance(gas, omega, source_admittance)
            .map(|response| response.1)
    }

    fn waves_with_source_admittance(
        &self,
        gas: &GasState,
        omega: f64,
        source_admittance: C64,
    ) -> Result<(C64, C64), DrivingPointError> {
        let boundary = source_admittance_boundary(source_admittance, self.port_impedance(gas))?;
        loaded_source_waves(self.scattering_chain(gas, omega)?, boundary)
    }

    /// Both observables from one physical block and the same source boundary sample.
    fn waves_with_source_impedance(
        &self, gas: &GasState, omega: f64, source_impedance: C64,
    ) -> Result<(C64, C64), DrivingPointError> {
        admit_source_impedance(source_impedance)?;
        let scattering = self.scattering_chain(gas, omega)?;
        let reflection =
            compose_chain_reflection(&[scattering], Some(source_impedance), ComposeOptions::default())?.reflection;
        Ok((reflection, loaded_source_transfer(scattering, source_impedance)?))
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
    /// Airway-only Poiseuille `R_dc` [Pa s/m^3].
    pub dc_resistance: f64,
    /// Source DC resistance [Pa s/m^3]. The impedance route records its exact
    /// authored value; the admittance route records its finite reciprocal when
    /// representable. None means sealed or unrepresentable, NEVER zero.
    /// With Some(R), the steady valve law is p = P_s - (dc_resistance + R) U_v.
    pub source_dc_resistance: Option<f64>,
    /// Source DC conductance [m^3/(Pa s)]. The admittance route records its exact
    /// authored value, including zero for a sealed limit. The impedance route
    /// records its finite reciprocal when representable; None is not zero.
    pub source_dc_admittance: Option<f64>,
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
    /// The source pressure, the known part `A` of this sample's received wave and the shadow's past (zero when
    /// there is no shadow), between `port` and `commit`.
    pending: Option<(f64, f64, f64)>,
    /// The reflectance history of the airway driven by the flow alone, the far-end source at rest: the response a
    /// valve ledger needs to read the work the SOURCE-FREE airway takes. `reflect` cannot give it, because the
    /// waves it holds were launched against the source's own contribution. Off unless asked
    /// ([`Self::track_source_free_response`]); the voice route never asks, and its arithmetic is untouched.
    free: Option<PastFir>,
    /// Samples closed by [`Self::commit`].
    commits: u64,
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
    sample_rate_hz: u32,
    transit: f64,
    dc: f64,
    ladder: &LineLadder,
    mut transfer_at: impl FnMut(f64) -> Result<C64, DrivingPointError>,
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
                *slot = transfer_at(core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt))?;
            }
        } else {
            // Even bins of this rung are the bins of the previous one.
            for (k_prev, &value) in t_pos.iter().enumerate() {
                next[2 * k_prev] = value;
            }
            for k in (1..half).step_by(2) {
                next[k] = transfer_at(core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt))?;
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
    /// plain FIR; a gas-wave or wall refusal at any bin of the ladder (one at a bin only the passivity gain samples, past
    /// the converged rung, keeps the passivity refusal instead).
    pub fn realize(
        chain: &AirwayChain,
        gas: &GasState,
        sample_rate_hz: u32,
        ladder: &LineLadder,
    ) -> Result<Self, DrivingPointError> {
        Self::realize_with_source_impedance(chain, gas, sample_rate_hz, ladder, 0.0, |_| Ok(C64::ZERO))
    }

    /// Realize an airway driven through an explicitly supplied passive source impedance.
    ///
    /// The boundary is `p_R + Z_L U_source = P_s`, with flow INTO the airway. The callback
    /// returns `Z_L(omega)` [Pa s/m^3] on the `exp(-iwt)` clock, at positive frequencies through
    /// Nyquist. `source_dc_resistance` is its exact, finite, nonnegative DC limit; no table is
    /// extrapolated and no stiffness pole is discarded. The supplied law must be causal:
    /// `Re Z_L >= 0` at finitely many samples does not establish causality or broadband passivity.
    ///
    /// Each distinct dyadic-grid frequency is sampled once, and both FIRs reuse its paired
    /// reflected-wave and source-pressure responses. Storage and calls are bounded by
    /// `ladder.max_taps / 2`; callback errors, including a source owner's cancellation, propagate, except a refusal
    /// other than a cancellation while the passivity gain samples its next rung: that is a decline of the gain, and the
    /// line keeps the passivity refusal. A line the gain tries takes that next rung's samples as well, within the same
    /// bound, even when it then keeps the refusal.
    /// After realization neither the callback nor any frequency solve runs in `port`/`commit`.
    /// The existing convergence and 4x-oversampled reflection-passivity gates still apply;
    /// source work is supplied by the open-circuit pressure signal, not clipped to unit gain.
    ///
    /// # Errors
    /// As [`Self::realize`], an invalid DC limit, or a non-finite/negative-real source sample (one at a bin only the
    /// passivity gain samples keeps the passivity refusal instead).
    /// A source with infinite DC impedance can use [`Self::realize_with_source_admittance`]
    /// for its zero-state wave response; this route admits finite DC resistance only.
    pub fn realize_with_source_impedance(
        chain: &AirwayChain,
        gas: &GasState,
        sample_rate_hz: u32,
        ladder: &LineLadder,
        source_dc_resistance: f64,
        mut source_impedance_at: impl FnMut(f64) -> Result<C64, DrivingPointError>,
    ) -> Result<Self, DrivingPointError> {
        let zs = chain.port_impedance(gas);
        let r_dc = chain.dc_resistance(gas);
        if !(zs.is_finite()
            && zs > 0.0
            && r_dc.is_finite()
            && r_dc >= 0.0
            && source_dc_resistance.is_finite()
            && source_dc_resistance >= 0.0)
        {
            return Err(DrivingPointError::Invalid {
                what: "airway port impedance and airway/source DC resistances must be finite and passive",
            });
        }
        // Normalize before adding the two physical DC resistances.
        let scale = zs.max(r_dc).max(source_dc_resistance);
        let resistance = r_dc / scale + source_dc_resistance / scale;
        let reference = zs / scale;
        let dc_reflection = (resistance - reference) / (resistance + reference);
        let dc_source = reference / (resistance + reference);
        let source_dc_admittance = if source_dc_resistance > 0.0 {
            let conductance = 1.0 / source_dc_resistance;
            conductance.is_finite().then_some(conductance)
        } else {
            None
        };
        Self::realize_source_response(
            chain,
            gas,
            sample_rate_hz,
            ladder,
            Some(source_dc_resistance),
            source_dc_admittance,
            (dc_reflection, dc_source),
            |omega| chain.waves_with_source_impedance(gas, omega, source_impedance_at(omega)?),
        )
    }

    /// Realize a causal source through its admittance Y_L(omega) [m^3/(Pa s)].
    ///
    /// The boundary is U_source = Y_L (P_s - p_R), in exp(-iwt). The exact
    /// nonnegative DC conductance is supplied separately; zero admits a sealed
    /// compliance or antiresonance without substituting a large resistance.
    /// In particular, a caller-authored gas compliance uses Y_L=-i*omega*C and
    /// source_dc_admittance=0. No lung volume or material data are inferred.
    ///
    /// Both original FIR ladders share one evaluation per positive dyadic-grid
    /// frequency. Their bounded convergence/passivity checks and the realtime
    /// implicit junction are unchanged. For zero DC conductance the exact bins
    /// are R(0)=1 and G(0)=0; no finite steady-flow resistance is claimed.
    ///
    /// This is a zero-state, fixed-geometry linear wave realization. An
    /// independently initialized reservoir pressure, changing gas volume and
    /// a local stored-energy ledger require a separate state-owning model.
    /// Finite-frequency passivity checks do not certify the caller's causality.
    ///
    /// # Errors
    /// As [`Self::realize`], an invalid DC conductance, a non-finite/negative-real
    /// source admittance, or the source callback's own refusal/cancellation (a refusal other than a cancellation at a
    /// bin only the passivity gain samples keeps the passivity refusal instead).
    pub fn realize_with_source_admittance(
        chain: &AirwayChain,
        gas: &GasState,
        sample_rate_hz: u32,
        ladder: &LineLadder,
        source_dc_admittance: f64,
        mut source_admittance_at: impl FnMut(f64) -> Result<C64, DrivingPointError>,
    ) -> Result<Self, DrivingPointError> {
        let zs = chain.port_impedance(gas);
        let r_dc = chain.dc_resistance(gas);
        if !(zs.is_finite()
            && zs > 0.0
            && r_dc.is_finite()
            && r_dc >= 0.0
            && source_dc_admittance.is_finite()
            && source_dc_admittance >= 0.0)
        {
            return Err(DrivingPointError::Invalid {
                what: "airway port impedance, airway resistance and source DC conductance must be finite and passive",
            });
        }
        let source_dc_resistance = if source_dc_admittance > 0.0 {
            let resistance = 1.0 / source_dc_admittance;
            resistance.is_finite().then_some(resistance)
        } else {
            None
        };
        Self::realize_source_response(
            chain,
            gas,
            sample_rate_hz,
            ladder,
            source_dc_resistance,
            Some(source_dc_admittance),
            admittance_dc_waves(zs, r_dc, source_dc_admittance),
            |omega| chain.waves_with_source_admittance(gas, omega, source_admittance_at(omega)?),
        )
    }

    /// Shared bounded realization; frequency responses and their exact DC
    /// boundary remain paired, independently of the source's coordinate chart.
    fn realize_source_response(
        chain: &AirwayChain,
        gas: &GasState,
        sample_rate_hz: u32,
        ladder: &LineLadder,
        source_dc_resistance: Option<f64>,
        source_dc_admittance: Option<f64>,
        dc: (f64, f64),
        mut waves_at: impl FnMut(f64) -> Result<(C64, C64), DrivingPointError>,
    ) -> Result<Self, DrivingPointError> {
        if sample_rate_hz == 0 {
            return Err(DrivingPointError::Invalid {
                what: "an airway needs a positive sample rate",
            });
        }
        let max_samples = check_ladder(ladder)? / 2;
        let zs = chain.port_impedance(gas);
        let r_dc = chain.dc_resistance(gas);
        let (dc_reflection, dc_source) = dc;
        let mut responses = BTreeMap::new();
        let mut response_at = |omega: f64| -> Result<(C64, C64), DrivingPointError> {
            let key = omega.to_bits();
            if let Some(&response) = responses.get(&key) {
                return Ok(response);
            }
            if responses.len() >= max_samples {
                return Err(DrivingPointError::Invalid {
                    what: "airway source sampling exceeded the declared ladder",
                });
            }
            let response = waves_at(omega)?;
            responses.insert(key, response);
            Ok(response)
        };
        let transit = chain.path_length_m() / gas.sound_speed * f64::from(sample_rate_hz);
        let round_trip = 2.0 * transit;
        let samples = round_trip.floor().max(0.0) as usize;
        let frac = (round_trip - samples as f64).clamp(0.0, 1.0 - f64::EPSILON);
        let (line, reflectance) = realize_characteristic_response(
            round_trip,
            LineDc::Explicit(dc_reflection),
            sample_rate_hz,
            zs,
            vec![PoleState::BulkDelay { samples, frac }],
            ladder,
            // Every receipt writer of an airway prints the gain: the voice, lip, reed and flue provenance, solo and in
            // the ensemble, and voice_rig's stderr.
            crate::driving_point::PassivityFallback::Gain,
            |omega| response_at(omega).map(|response| response.0),
        )?;
        let r_taps = line
            .quiescent_fir()
            .ok_or(DrivingPointError::Invalid {
                what: "the airway reflectance line is not a plain FIR",
            })?
            .1
            .to_vec();
        let (g_taps, source_rungs) =
            realize_source_transfer(sample_rate_hz, transit, dc_source, ladder, |omega| {
                response_at(omega).map(|response| response.1)
            })?;
        let (r0, g0) = (r_taps[0], g_taps[0]);
        if !(r0.abs() < 1.0) {
            return Err(DrivingPointError::Invalid {
                what: "the airway reflects its same-sample wave whole",
            });
        }
        let series_impedance = zs * (1.0 + r0) / (1.0 - r0);
        let receipt = AirwayReceipt {
            sections: chain.sections().len(),
            port_impedance: zs,
            dc_resistance: r_dc,
            source_dc_resistance,
            source_dc_admittance,
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
            free: None,
            commits: 0,
            receipt,
        })
    }

    /// The Thevenin source the valve sees this sample with the open-circuit source pressure `source_pa` at the far end.
    /// A supplied source impedance remains connected when this pressure is zero.
    /// Nothing moves until [`Self::commit`]; a second call before it replaces the first.
    ///
    /// # Errors
    /// A non-finite source pressure or port.
    pub fn port(&mut self, source_pa: f64) -> Result<UpstreamPort, DrivingPointError> {
        if !source_pa.is_finite() {
            return Err(DrivingPointError::Invalid { what: "airway source pressure must be finite" });
        }
        let known = self.reflect.past() + (self.g0 * source_pa + self.source.past());
        let free_past = self.free.as_ref().map_or(0.0, PastFir::past);
        let open_circuit_pa = 2.0 * known / (1.0 - self.r0);
        if !open_circuit_pa.is_finite() {
            return Err(DrivingPointError::Invalid { what: "airway port left the finite set" });
        }
        self.pending = Some((source_pa, known, free_past));
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
        let (source_pa, known, free_past) =
            self.pending.ok_or(DrivingPointError::Invalid { what: "an airway commit needs this sample's port" })?;
        let b = (known - self.zs * flow_m3_s) / (1.0 - self.r0);
        let p = (1.0 + self.r0) * b + known;
        let b_free = (free_past - self.zs * flow_m3_s) / (1.0 - self.r0);
        if !(flow_m3_s.is_finite() && b.is_finite() && p.is_finite() && (self.free.is_none() || b_free.is_finite())) {
            return Err(DrivingPointError::Invalid { what: "airway wave left the finite set" });
        }
        self.pending = None;
        self.reflect.push(b);
        self.source.push(source_pa);
        if let Some(free) = &mut self.free {
            free.push(b_free);
        }
        self.commits += 1;
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

    /// Start the flow-only shadow of the reflectance history, so that [`Self::source_free_open_circuit_pa`] can
    /// answer. By linearity the airway's port is the sum of its response to the flow history with the source at
    /// rest and its response to the source; the shadow carries the first alone: `b_free = (free_past - Z_s U) /
    /// (1 - r_0)` pushed at every commit. It can only start on an airway that has not advanced, since it would
    /// otherwise miss the history already launched. Calling it again is a no-op.
    ///
    /// # Errors
    /// An airway that has already closed a sample.
    pub fn track_source_free_response(&mut self) -> Result<(), DrivingPointError> {
        if self.free.is_some() {
            return Ok(());
        }
        if self.commits != 0 {
            return Err(DrivingPointError::Invalid {
                what: "the source-free response can only start on an airway that has not advanced",
            });
        }
        self.free = Some(PastFir::new(self.reflect.taps.clone()));
        Ok(())
    }

    /// The open-circuit pressure of the airway with its far-end source at rest, between [`Self::port`] and
    /// [`Self::commit`]: `2 free_past / (1 - r_0)` [Pa], the flow history's part of `open_circuit_pa`. `None`
    /// without a pending port or without [`Self::track_source_free_response`]. A valve ledger reads the work the
    /// source-free airway takes from the valve with this.
    #[must_use]
    pub fn source_free_open_circuit_pa(&self) -> Option<f64> {
        let (_, _, free_past) = self.pending?;
        self.free.as_ref()?;
        Some(2.0 * free_past / (1.0 - self.r0))
    }

    /// Samples closed so far.
    #[must_use]
    pub const fn commit_count(&self) -> u64 {
        self.commits
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

    /// G1: wave observables preserve the original independently formed ABCD
    /// law where it is representable, including walls and parallel branches.
    #[test]
    fn ua_008_bounded_waves_match_the_physical_source_boundary() {
        let gas = air();
        for loss in [LossModel::Lossless, LossModel::Bessel] {
            for lining in [None, Some(wall())] {
                let chain = AirwayChain::new(vec![
                    tube(0.009, 0.12, 1, lining),
                    tube(0.004, 0.07, 3, lining),
                    tube(0.002, 0.03, 7, lining),
                ], loss).unwrap();
                let zs = C64::from_re(chain.port_impedance(&gas));
                for hz in [0.01, 51.0, 713.0, 2901.0, 13007.0] {
                    let omega = core::f64::consts::TAU * hz;
                    let m = chain.chain_matrix(&gas, omega).unwrap();
                    let denominator = m[0][1] + zs * m[1][1];
                    let reflection = (m[0][1] - zs * m[1][1]) / denominator;
                    let source = zs / denominator;
                    assert!((chain.reflectance(&gas, omega).unwrap() - reflection).abs() < 2e-11);
                    assert!(rel(chain.source_transfer(&gas, omega).unwrap(), source) < 2e-10);
                    assert!(rel(chain.input_impedance(&gas, omega).unwrap(), m[0][1] / m[1][1]) < 2e-9);
                }
            }
        }
    }

    /// G1/G3: an opaque real duct still reflects from its entrance impedance.
    /// Growth of the inverse transfer matrix cannot destroy that finite field.
    #[test]
    fn ua_009_opaque_airways_retain_reflection_and_attenuated_source() {
        let gas = air();
        let radius = 0.0005;
        let count = 8;
        let omega = core::f64::consts::TAU * 2000.0;
        let wave = segment_wave(&gas, radius, omega, LossModel::Bessel).unwrap();
        let z = wave.characteristic_impedance.scale(1.0 / f64::from(count));
        for attenuation in [40.0, 350.0, 900.0] {
            let length = attenuation / wave.wavenumber.im;
            let chain = AirwayChain::new(vec![tube(radius, length, count, None)], LossModel::Bessel).unwrap();
            let zs = C64::from_re(chain.port_impedance(&gas));
            let entrance = (z - zs) / (z + zs);
            let phase = wave.wavenumber.re * length;
            let travel = C64::new(phase.cos(), phase.sin()).scale((-attenuation).exp());
            // Round-trip corrections are below 2e-35 on this declared range.
            let source = zs.scale(2.0) * travel / (zs + z);
            assert!((chain.reflectance(&gas, omega).unwrap() - entrance).abs() < 2e-12);
            let got_source = chain.source_transfer(&gas, omega).unwrap();
            assert!(rel(got_source, source) < 2e-10, "source {got_source:?}, expected {source:?}");
            assert!(rel(chain.input_impedance(&gas, omega).unwrap(), z) < 2e-12);
            let split = AirwayChain::new(vec![tube(radius, length / 16.0, count, None); 16], LossModel::Bessel).unwrap();
            assert!((split.reflectance(&gas, omega).unwrap() - entrance).abs() < 2e-12);
            assert!(rel(split.source_transfer(&gas, omega).unwrap(), source) < 2e-10);
            if attenuation > 800.0 {
                assert!(chain.chain_matrix(&gas, omega).is_err());
                assert_eq!(got_source, C64::ZERO);
            }
        }
    }

    /// G1: an ideal pressure release returns the negative round-trip wave;
    /// the far source travels once. At a quarter wave the input Z has a pole,
    /// while both physically observable traveling waves stay finite.
    #[test]
    fn ua_010_lossless_reflection_phase_survives_an_impedance_pole() {
        let gas = air();
        let length = 0.13;
        let chain = AirwayChain::new(vec![tube(0.007, length, 3, None)], LossModel::Lossless).unwrap();
        for phase in [0.1, core::f64::consts::FRAC_PI_2, core::f64::consts::PI, 5.0] {
            let omega = phase * gas.sound_speed / length;
            let expected_r = -C64::new((2.0 * phase).cos(), (2.0 * phase).sin());
            let expected_g = C64::new(phase.cos(), phase.sin());
            assert!((chain.reflectance(&gas, omega).unwrap() - expected_r).abs() < 2e-14);
            assert!((chain.source_transfer(&gas, omega).unwrap() - expected_g).abs() < 2e-14);
        }
        let quarter_wave = core::f64::consts::FRAC_PI_2 * gas.sound_speed / length;
        assert!(chain.input_impedance(&gas, quarter_wave).is_err());
        let mut inconsistent = gas;
        inconsistent.sound_speed *= 1.01;
        assert!(chain.reflectance(&inconsistent, quarter_wave).is_err());
        assert!(chain.source_transfer(&inconsistent, quarter_wave).is_err());
        assert!(chain.reflectance(&gas, f64::NAN).is_err());
        assert!(chain.source_transfer(&gas, 0.0).is_err());
    }

    /// A finite generator is still a load when its open-circuit pressure is zero.
    /// The independent ABCD boundary uses p_R = P_s + Z_L U_R with U_R away from the valve.
    #[test]
    fn ua_011_finite_source_matches_the_loaded_abcd_boundary() {
        let gas = air();
        for loss in [LossModel::Lossless, LossModel::Bessel] {
            let chain = AirwayChain::new(vec![
                tube(0.009, 0.12, 1, Some(wall())),
                tube(0.004, 0.07, 3, None),
                tube(0.002, 0.03, 7, Some(wall())),
            ], loss).unwrap();
            let zs = C64::from_re(chain.port_impedance(&gas));
            for hz in [0.01, 51.0, 713.0, 2901.0, 13007.0] {
                let omega = core::f64::consts::TAU * hz;
                let m = chain.chain_matrix(&gas, omega).unwrap();
                for load in [C64::ZERO, zs, zs.scale(3.0), C64::new(0.4, -0.9) * zs, C64::new(0.2, 2.0) * zs] {
                    let numerator = m[0][0] * load + m[0][1];
                    let flow_denominator = m[1][0] * load + m[1][1];
                    let denominator = numerator + zs * flow_denominator;
                    let expected_r = (numerator - zs * flow_denominator) / denominator;
                    let expected_g = zs / denominator;
                    let reflection = chain.reflectance_with_source_impedance(&gas, omega, load).unwrap();
                    let source = chain.source_transfer_with_source_impedance(&gas, omega, load).unwrap();
                    let impedance = chain.input_impedance_with_source_impedance(&gas, omega, load).unwrap();
                    assert!((reflection - expected_r).abs() < 2e-11);
                    assert!(rel(source, expected_g) < 2e-10);
                    assert!(rel(impedance, numerator / flow_denominator) < 2e-9);
                    assert!(reflection.abs() <= 1.0 + 2e-12);
                }
            }
        }
    }

    /// Match, reactive storage and the rigid limit pin phase, source normalization and absorbed power.
    #[test]
    fn ua_012_finite_source_preserves_wave_power_and_matched_launch() {
        let gas = air();
        let length = 0.13;
        let chain = AirwayChain::new(vec![tube(0.007, length, 3, None)], LossModel::Lossless).unwrap();
        let zs = chain.port_impedance(&gas);
        for phase in [0.1, core::f64::consts::FRAC_PI_2, core::f64::consts::PI, 5.0] {
            let omega = phase * gas.sound_speed / length;
            let travel = C64::new(phase.cos(), phase.sin());
            for load in [C64::ZERO, C64::from_re(zs), C64::new(0.3 * zs, -0.7 * zs), C64::new(0.0, 2.0 * zs)] {
                let denominator = load + C64::from_re(zs);
                let local_r = (load - C64::from_re(zs)) / denominator;
                let local_g = C64::from_re(zs) / denominator;
                let reflection = chain.reflectance_with_source_impedance(&gas, omega, load).unwrap();
                let source = chain.source_transfer_with_source_impedance(&gas, omega, load).unwrap();
                assert!((reflection - local_r * travel * travel).abs() < 2e-14);
                assert!((source - local_g * travel).abs() < 2e-14);
                let absorbed = 4.0 * zs * load.re / denominator.norm_sq();
                assert!((1.0 - reflection.norm_sq() - absorbed).abs() < 3e-14);
            }
            let rigid = C64::new(f64::MAX, 0.5 * f64::MAX);
            let reflection = chain.reflectance_with_source_impedance(&gas, omega, rigid).unwrap();
            let source = chain.source_transfer_with_source_impedance(&gas, omega, rigid).unwrap();
            assert!((reflection - travel * travel).abs() < 2e-14);
            assert!(source.re.is_finite() && source.im.is_finite() && source.abs() < 1e-300);
        }
    }

    /// A caller-authored R_s + (R_p || C) source changes both FIRs. The identical sampled
    /// boundary is reused across the two adaptive ladders; the sample loop sees its DC drop.
    #[test]
    fn ua_013_reactive_source_reaches_both_firs_and_the_realtime_port() {
        let gas = air();
        let length = 32.0 * gas.sound_speed / f64::from(RATE);
        let chain = AirwayChain::new(vec![tube(0.007, length, 1, None)], LossModel::Lossless).unwrap();
        let zs = chain.port_impedance(&gas);
        let source_dc = 1.5 * zs;
        let relaxation_s = 4.0 / f64::from(RATE);
        let source_z = |omega: f64| C64::from_re(0.5 * zs)
            + C64::from_re(zs) / C64::new(1.0, -omega * relaxation_s);
        let mut frequencies = std::collections::BTreeSet::new();
        let settings = ladder();
        let mut airway = UpstreamAirway::realize_with_source_impedance(
            &chain, &gas, RATE, &settings, source_dc,
            |omega| {
                assert!(frequencies.insert(omega.to_bits()), "a boundary sample was queried twice");
                Ok(source_z(omega))
            },
        ).unwrap();
        assert_eq!(frequencies.len(), airway.reflectance_taps().len().max(airway.source_taps().len()) / 2);
        assert_eq!(airway.receipt().source_dc_resistance, Some(source_dc));
        assert_eq!(airway.receipt().dc_resistance, 0.0);
        for hz in [7.0, 310.0, 1013.0, 2207.0, 5119.0] {
            let omega = core::f64::consts::TAU * hz;
            let load = source_z(omega);
            let expected_r = chain.reflectance_with_source_impedance(&gas, omega, load).unwrap();
            let expected_g = chain.source_transfer_with_source_impedance(&gas, omega, load).unwrap();
            assert!((dtft(airway.reflectance_taps(), omega) - expected_r).abs() < 1e-3);
            assert!((dtft(airway.source_taps(), omega) - expected_g).abs() < 1e-3);
        }
        let dc_r: f64 = airway.reflectance_taps().iter().sum();
        let dc_g: f64 = airway.source_taps().iter().sum();
        assert!((2.0 * dc_g / (1.0 - dc_r) - 1.0).abs() < 1e-10);
        assert!((zs * (1.0 + dc_r) / (1.0 - dc_r) / source_dc - 1.0).abs() < 1e-10);
        let mut silent_source = airway.clone();
        // Include the FIR's numerical tail in the valve feedback's settling time.
        let settle = 3 * airway.reflectance_taps().len() + 8;
        let source_pressure = 800.0;
        let flow = 2.0e-4;
        let (port, pressure) = hold(&mut airway, source_pressure, flow, settle);
        assert!((pressure - (source_pressure - source_dc * flow)).abs() < 1e-6);
        assert!((port.pressure_pa(flow) - pressure).abs() < 1e-10);
        // Source zero leaves the finite source resistance connected.
        let (_, pressure) = hold(&mut silent_source, 0.0, flow, settle);
        assert!((pressure + source_dc * flow).abs() < 1e-6,
            "residual {}, reflection taps {}, source taps {}", pressure + source_dc * flow,
            airway.reflectance_taps().len(), airway.source_taps().len());
    }

    /// Unphysical source samples, inconsistent finite DC admission and cancellation never produce an audio owner.
    #[test]
    fn ua_014_source_admission_and_cancellation_are_explicit() {
        let gas = air();
        let chain = AirwayChain::new(vec![tube(0.007, 0.13, 1, None)], LossModel::Lossless).unwrap();
        let omega = core::f64::consts::TAU * 713.0;
        for impedance in [C64::new(-1.0, 0.0), C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
            assert!(chain.reflectance_with_source_impedance(&gas, omega, impedance).is_err());
            assert!(chain.source_transfer_with_source_impedance(&gas, omega, impedance).is_err());
            assert!(chain.input_impedance_with_source_impedance(&gas, omega, impedance).is_err());
        }
        for dc in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(UpstreamAirway::realize_with_source_impedance(
                &chain, &gas, RATE, &ladder(), dc,
                |_| panic!("invalid DC must refuse before querying the source"),
            ).is_err());
        }
        let mut calls = 0;
        let stopped = UpstreamAirway::realize_with_source_impedance(
            &chain, &gas, RATE, &ladder(), 0.0,
            |_| {
                calls += 1;
                Err(fs_duct::DuctError::Cancelled { completed_steps: calls }.into())
            },
        );
        assert!(matches!(stopped, Err(DrivingPointError::Duct(fs_duct::DuctError::Cancelled { completed_steps: 1 }))));
        assert_eq!(calls, 1);
    }

    /// G1/G3: use the same physical load in its impedance and admittance charts.
    /// Sealed and zero-admittance resonance limits never manufacture infinity.
    #[test]
    fn ua_015_source_admittance_retains_closed_and_resonant_boundaries() {
        let gas = air();
        for loss in [LossModel::Lossless, LossModel::Bessel] {
            let chain = AirwayChain::new(
                vec![
                    tube(0.009, 0.12, 1, Some(wall())),
                    tube(0.004, 0.07, 3, None),
                ],
                loss,
            )
            .unwrap();
            let zs = chain.port_impedance(&gas);
            for hz in [0.01, 51.0, 713.0, 2901.0, 13007.0] {
                let omega = core::f64::consts::TAU * hz;
                for load in [
                    C64::from_re(zs),
                    C64::new(0.3 * zs, -0.7 * zs),
                    C64::new(0.0, 2.0 * zs),
                ] {
                    let y = load.recip();
                    let r = chain
                        .reflectance_with_source_admittance(&gas, omega, y)
                        .unwrap();
                    let g = chain
                        .source_transfer_with_source_admittance(&gas, omega, y)
                        .unwrap();
                    let z = chain
                        .input_impedance_with_source_admittance(&gas, omega, y)
                        .unwrap();
                    assert!(
                        (r - chain
                            .reflectance_with_source_impedance(&gas, omega, load)
                            .unwrap())
                        .abs()
                            < 2e-12
                    );
                    assert!(
                        (g - chain
                            .source_transfer_with_source_impedance(&gas, omega, load)
                            .unwrap())
                        .abs()
                            < 2e-12
                    );
                    assert!(
                        rel(
                            z,
                            chain
                                .input_impedance_with_source_impedance(&gas, omega, load)
                                .unwrap()
                        ) < 2e-10
                    );
                }
                let m = chain.chain_matrix(&gas, omega).unwrap();
                let denominator = m[0][0] + m[1][0].scale(zs);
                let expected = (m[0][0] - m[1][0].scale(zs)) / denominator;
                assert!(
                    (chain
                        .reflectance_with_source_admittance(&gas, omega, C64::ZERO)
                        .unwrap()
                        - expected)
                        .abs()
                        < 2e-11
                );
                assert_eq!(
                    chain
                        .source_transfer_with_source_admittance(&gas, omega, C64::ZERO)
                        .unwrap(),
                    C64::ZERO
                );
            }
        }
        let length = 0.13;
        let chain =
            AirwayChain::new(vec![tube(0.007, length, 1, None)], LossModel::Lossless).unwrap();
        let zs = chain.port_impedance(&gas);
        let compliance = 1.0 / (1000.0 * zs);
        for phase in [
            0.1,
            core::f64::consts::FRAC_PI_2,
            core::f64::consts::PI,
            5.0,
        ] {
            let omega = phase * gas.sound_speed / length;
            let y = C64::new(0.0, -omega * compliance);
            let travel = C64::new(phase.cos(), phase.sin());
            let local_r = (C64::ONE - y.scale(zs)) / (C64::ONE + y.scale(zs));
            let local_g = y.scale(zs) / (C64::ONE + y.scale(zs));
            let r = chain
                .reflectance_with_source_admittance(&gas, omega, y)
                .unwrap();
            let g = chain
                .source_transfer_with_source_admittance(&gas, omega, y)
                .unwrap();
            assert!((r - local_r * travel * travel).abs() < 2e-14);
            assert!((g - local_g * travel).abs() < 2e-14);
            assert!((r.norm_sq() - 1.0).abs() < 3e-14);
        }
        assert!(
            chain
                .input_impedance_with_source_admittance(
                    &gas,
                    core::f64::consts::PI * gas.sound_speed / length,
                    C64::ZERO,
                )
                .is_err()
        );
    }

    /// G1/G3: the original two ladders admit a closed compliance with exact DC
    /// bins, a single source sample per frequency and its original finite-band phase.
    #[test]
    fn ua_016_compliant_source_reaches_the_original_wave_realization() {
        let gas = air();
        let length = 32.0 * gas.sound_speed / f64::from(RATE);
        let chain =
            AirwayChain::new(vec![tube(0.007, length, 1, None)], LossModel::Bessel).unwrap();
        let compliance = 8.0 / (f64::from(RATE) * chain.port_impedance(&gas));
        let mut frequencies = std::collections::BTreeSet::new();
        let source_y = |omega: f64| C64::new(0.0, -omega * compliance);
        let airway = UpstreamAirway::realize_with_source_admittance(
            &chain,
            &gas,
            RATE,
            &ladder(),
            0.0,
            |omega| {
                assert!(frequencies.insert(omega.to_bits()));
                Ok(source_y(omega))
            },
        )
        .unwrap();
        assert_eq!(airway.receipt().source_dc_resistance, None);
        assert_eq!(airway.receipt().source_dc_admittance, Some(0.0));
        assert_eq!(
            frequencies.len(),
            airway
                .reflectance_taps()
                .len()
                .max(airway.source_taps().len())
                / 2
        );
        assert!((airway.reflectance_taps().iter().sum::<f64>() - 1.0).abs() < 2e-12);
        assert!(airway.source_taps().iter().sum::<f64>().abs() < 2e-12);
        for hz in [7.0, 310.0, 1013.0, 2207.0, 5119.0] {
            let omega = core::f64::consts::TAU * hz;
            let expected_r = chain
                .reflectance_with_source_admittance(&gas, omega, source_y(omega))
                .unwrap();
            let expected_g = chain
                .source_transfer_with_source_admittance(&gas, omega, source_y(omega))
                .unwrap();
            assert!((dtft(airway.reflectance_taps(), omega) - expected_r).abs() < 1e-3);
            assert!((dtft(airway.source_taps(), omega) - expected_g).abs() < 1e-3);
        }
        let sealed = UpstreamAirway::realize_with_source_admittance(
            &chain,
            &gas,
            RATE,
            &ladder(),
            0.0,
            |_| Ok(C64::ZERO),
        )
        .unwrap();
        assert!(sealed.source_taps().iter().all(|&tap| tap == 0.0));
        assert!((sealed.reflectance_taps().iter().sum::<f64>() - 1.0).abs() < 2e-12);
    }

    /// G3/G5: reciprocal source coordinates retain the original finite-DC
    /// realtime junction, including silent-source loading and same-sample flow.
    #[test]
    fn ua_017_admittance_and_impedance_realizations_drive_the_same_valve() {
        let gas = air();
        let length = 32.0 * gas.sound_speed / f64::from(RATE);
        let chain =
            AirwayChain::new(vec![tube(0.007, length, 1, None)], LossModel::Lossless).unwrap();
        let zs = chain.port_impedance(&gas);
        let tau = 4.0 / f64::from(RATE);
        let source_z =
            |omega: f64| C64::from_re(0.5 * zs) + C64::from_re(zs) / C64::new(1.0, -omega * tau);
        let mut z_owner = UpstreamAirway::realize_with_source_impedance(
            &chain,
            &gas,
            RATE,
            &ladder(),
            1.5 * zs,
            |omega| Ok(source_z(omega)),
        )
        .unwrap();
        let mut y_owner = UpstreamAirway::realize_with_source_admittance(
            &chain,
            &gas,
            RATE,
            &ladder(),
            1.0 / (1.5 * zs),
            |omega| Ok(source_z(omega).recip()),
        )
        .unwrap();
        assert_eq!(
            z_owner.reflectance_taps().len(),
            y_owner.reflectance_taps().len()
        );
        assert_eq!(z_owner.source_taps().len(), y_owner.source_taps().len());
        for n in 0..512 {
            let phase = core::f64::consts::TAU * 713.0 * f64::from(n) / f64::from(RATE);
            let pressure = if n < 256 { 800.0 * phase.sin() } else { 0.0 };
            let flow = 2e-4 * (0.7 * phase).cos();
            let a = z_owner.port(pressure).unwrap();
            let b = y_owner.port(pressure).unwrap();
            assert!((a.open_circuit_pa - b.open_circuit_pa).abs() < 2e-9);
            assert!((a.impedance / b.impedance - 1.0).abs() < 2e-14);
            assert!((z_owner.commit(flow).unwrap() - y_owner.commit(flow).unwrap()).abs() < 2e-9);
        }
    }

    /// G0/G4: no impedance reciprocal, infinity arithmetic or callback on invalid
    /// input; both scale charts retain the matched, sealed and release limits.
    #[test]
    fn ua_018_admittance_admission_and_extreme_scales_are_explicit() {
        for reference in [1e-200, 1.0, 1e200] {
            for y in [
                C64::ZERO,
                C64::new(1e-300, -1e-300),
                C64::new(1.0, 2.0),
                C64::new(f64::MAX, -f64::MAX),
            ] {
                let (r, g) = source_admittance_boundary(y, reference).unwrap();
                assert!(r.re.is_finite() && r.im.is_finite() && r.abs() <= 1.0 + 4e-16);
                assert!(g.re.is_finite() && g.im.is_finite());
                assert!((r + g.scale(2.0) - C64::ONE).abs() < 5e-16);
            }
            let (r, g) =
                source_admittance_boundary(C64::from_re(1.0 / reference), reference).unwrap();
            assert!(r.abs() < 2e-16);
            assert!((g.re - 0.5).abs() < 2e-16);
        }
        for (reference, resistance, conductance) in [
            (1.0, 1.0, 0.0),
            (1e200, 1e200, 1e200),
            (1e-200, 1e-200, 1e200),
            (1e200, 1e200, 1e-200),
            (1e-200, 1e-200, 1e-200),
        ] {
            let (r, g) = admittance_dc_waves(reference, resistance, conductance);
            assert!(r.is_finite() && r.abs() <= 1.0);
            assert!(g.is_finite() && (0.0..=1.0).contains(&g));
        }
        let gas = air();
        let chain =
            AirwayChain::new(vec![tube(0.007, 0.13, 1, None)], LossModel::Lossless).unwrap();
        for y in [
            C64::new(-1.0, 0.0),
            C64::new(f64::NAN, 0.0),
            C64::new(0.0, f64::INFINITY),
        ] {
            assert!(
                chain
                    .reflectance_with_source_admittance(&gas, 713.0, y)
                    .is_err()
            );
            assert!(
                chain
                    .source_transfer_with_source_admittance(&gas, 713.0, y)
                    .is_err()
            );
            assert!(
                chain
                    .input_impedance_with_source_admittance(&gas, 713.0, y)
                    .is_err()
            );
        }
        for dc in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                UpstreamAirway::realize_with_source_admittance(
                    &chain,
                    &gas,
                    RATE,
                    &ladder(),
                    dc,
                    |_| panic!("invalid DC must refuse before a source call"),
                )
                .is_err()
            );
        }
        let mut calls = 0;
        let stopped = UpstreamAirway::realize_with_source_admittance(
            &chain,
            &gas,
            RATE,
            &ladder(),
            0.0,
            |_| {
                calls += 1;
                Err(fs_duct::DuctError::Cancelled {
                    completed_steps: calls,
                }
                .into())
            },
        );
        assert!(matches!(
            stopped,
            Err(DrivingPointError::Duct(fs_duct::DuctError::Cancelled {
                completed_steps: 1
            }))
        ));
        assert_eq!(calls, 1);
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
