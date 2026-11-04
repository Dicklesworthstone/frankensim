//! Peng--Robinson real-gas mechanical EOS with classical quadratic mixing.
//!
//! Callers supply critical properties, acentric factors, molar masses,
//! composition and optional binary interactions from their physical source.
//! The evaluator retains every mechanically admissible cubic root and exposes
//! fugacity coefficients. The base evaluator neither fits coefficients nor
//! selects phase amounts; the flash and stability submodules use those same
//! source-bound fugacities for phase-equilibrium calculations.

use core::fmt;

use fs_math::det;
use fs_qty::{Density, MolarMass, Pressure, Temperature};

use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

/// Bounded PT vapor-liquid equilibrium on this same cubic EOS.
pub mod flash;
/// Bounded homogeneous/two-phase selection with explicit stability-screen evidence.
pub mod phase_screen;
/// Bounded tangent-plane searches on the original mixture fugacities.
pub mod stability;
/// Pure-fluid coexistence and conservative phase amounts on this same EOS.
pub mod saturation;

/// Maximum number of components in one bounded mixture.
pub const MAX_PENG_ROBINSON_COMPONENTS_V1: usize = 64;
/// Fixed operation-tree version.
pub const PENG_ROBINSON_EVALUATOR_VERSION_V1: u32 = 4;

const ROOT_TOL: f64 = 128.0 * f64::EPSILON;
const COMPOSITION_TOL: f64 = 256.0 * f64::EPSILON;
const KIJ_TOL: f64 = 256.0 * f64::EPSILON;

/// Pure-component constants for the 1976 Peng--Robinson law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PengRobinsonComponentV1 {
    pub critical_temperature: Temperature,
    pub critical_pressure: Pressure,
    pub acentric_factor: f64,
    pub molar_mass: MolarMass,
}

/// Root selection within one fixed-composition PT evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PengRobinsonPhaseV1 {
    /// Largest admissible Z root.
    Vapor,
    /// Smallest admissible Z root.
    Liquid,
    /// Minimum residual Gibbs root at the retained composition.
    Stable,
}

/// Construction or evaluation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PengRobinsonErrorV1 {
    EmptyMixture,
    TooManyComponents { offered: usize, limit: usize },
    InvalidComponent { component: usize, field: &'static str, bits: u64 },
    CompositionLength { components: usize, fractions: usize },
    InvalidMoleFraction { component: usize, bits: u64 },
    CompositionNotNormalized { sum: f64 },
    BinaryInteractionShape { expected: usize, found: usize },
    InvalidBinaryInteraction { row: usize, column: usize, bits: u64 },
    AsymmetricBinaryInteraction { row: usize, column: usize },
    NonzeroBinaryInteractionDiagonal { component: usize, bits: u64 },
    InvalidTemperature { bits: u64 },
    InvalidPressure { bits: u64 },
    InvalidDensity { bits: u64 },
    OutsideCovolume,
    MechanicallyUnstable,
    Unrepresentable { what: &'static str },
    NoMechanicalRoot,
}

impl fmt::Display for PengRobinsonErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Peng-Robinson EOS refusal: {self:?}")
    }
}
impl std::error::Error for PengRobinsonErrorV1 {}

/// Immutable parameter/composition binding.
#[derive(Debug, Clone, PartialEq)]
pub struct PengRobinsonMixtureV1 {
    components: Vec<PengRobinsonComponentV1>,
    mole_fractions: Vec<f64>,
    binary_interactions: Vec<f64>,
    molar_mass: MolarMass,
}

impl PengRobinsonMixtureV1 {
    /// Admit a fixed-composition mixture. An empty interaction matrix means
    /// k_ij=0. A supplied matrix is row-major, symmetric, finite, has exactly
    /// zero diagonal, and must keep every (1-k_ij) nonnegative.
    pub fn new(
        components: Vec<PengRobinsonComponentV1>,
        mole_fractions: Vec<f64>,
        binary_interactions: Vec<f64>,
    ) -> Result<Self, PengRobinsonErrorV1> {
        let n = components.len();
        if n == 0 {
            return Err(PengRobinsonErrorV1::EmptyMixture);
        }
        if n > MAX_PENG_ROBINSON_COMPONENTS_V1 {
            return Err(PengRobinsonErrorV1::TooManyComponents {
                offered: n,
                limit: MAX_PENG_ROBINSON_COMPONENTS_V1,
            });
        }
        if mole_fractions.len() != n {
            return Err(PengRobinsonErrorV1::CompositionLength {
                components: n,
                fractions: mole_fractions.len(),
            });
        }
        for (i, component) in components.iter().enumerate() {
            for (field, value) in [
                ("critical temperature", component.critical_temperature.value()),
                ("critical pressure", component.critical_pressure.value()),
                ("molar mass", component.molar_mass.value()),
            ] {
                if !value.is_finite() || value <= 0.0 {
                    return Err(PengRobinsonErrorV1::InvalidComponent {
                        component: i,
                        field,
                        bits: value.to_bits(),
                    });
                }
            }
            if !component.acentric_factor.is_finite() {
                return Err(PengRobinsonErrorV1::InvalidComponent {
                    component: i,
                    field: "acentric factor",
                    bits: component.acentric_factor.to_bits(),
                });
            }
        }

        let mut xsum = 0.0;
        let mut molar_mass = 0.0;
        for (i, (&x, component)) in mole_fractions.iter().zip(&components).enumerate() {
            if !x.is_finite() || x < 0.0 {
                return Err(PengRobinsonErrorV1::InvalidMoleFraction {
                    component: i,
                    bits: x.to_bits(),
                });
            }
            xsum += x;
            molar_mass += x * component.molar_mass.value();
        }
        if !xsum.is_finite() || (xsum - 1.0).abs() > COMPOSITION_TOL * n as f64 {
            return Err(PengRobinsonErrorV1::CompositionNotNormalized { sum: xsum });
        }
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "mixture molar mass",
            });
        }

        let kij = if binary_interactions.is_empty() {
            vec![0.0; n * n]
        } else {
            if binary_interactions.len() != n * n {
                return Err(PengRobinsonErrorV1::BinaryInteractionShape {
                    expected: n * n,
                    found: binary_interactions.len(),
                });
            }
            binary_interactions
        };
        for i in 0..n {
            if kij[i * n + i].to_bits() != 0.0f64.to_bits() {
                return Err(PengRobinsonErrorV1::NonzeroBinaryInteractionDiagonal {
                    component: i,
                    bits: kij[i * n + i].to_bits(),
                });
            }
            for j in 0..n {
                let value = kij[i * n + j];
                if !value.is_finite() || value > 1.0 {
                    return Err(PengRobinsonErrorV1::InvalidBinaryInteraction {
                        row: i,
                        column: j,
                        bits: value.to_bits(),
                    });
                }
                if (value - kij[j * n + i]).abs()
                    > KIJ_TOL * value.abs().max(kij[j * n + i].abs()).max(1.0)
                {
                    return Err(PengRobinsonErrorV1::AsymmetricBinaryInteraction {
                        row: i,
                        column: j,
                    });
                }
            }
        }
        Ok(Self {
            components,
            mole_fractions,
            binary_interactions: kij,
            molar_mass: MolarMass::new(molar_mass),
        })
    }

    #[must_use]
    pub fn components(&self) -> &[PengRobinsonComponentV1] {
        &self.components
    }
    #[must_use]
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }
    #[must_use]
    pub fn binary_interactions(&self) -> &[f64] {
        &self.binary_interactions
    }
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }

    fn mixing(&self, t: f64) -> Result<MixedParameters, PengRobinsonErrorV1> {
        if !t.is_finite() || t <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidTemperature { bits: t.to_bits() });
        }

        let n = self.components.len();
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let mut ai = vec![0.0; n];
        let mut dai_dt = vec![0.0; n];
        let mut d2ai_dt2 = vec![0.0; n];
        let mut bi = vec![0.0; n];
        for i in 0..n {
            let component = self.components[i];
            let tc = component.critical_temperature.value();
            let pc = component.critical_pressure.value();
            let omega = component.acentric_factor;
            let kappa = 0.374_64 + 1.542_26 * omega - 0.269_92 * omega * omega;
            let sqrt_reduced_t = det::sqrt(t / tc);
            let alpha_base = 1.0 + kappa * (1.0 - sqrt_reduced_t);
            let alpha = alpha_base * alpha_base;
            let a0 = 0.457_24 * r * r * tc * tc / pc;
            let dbase_dt = -kappa * sqrt_reduced_t / (2.0 * t);
            let d2base_dt2 = kappa * sqrt_reduced_t / (4.0 * t * t);
            ai[i] = a0 * alpha;
            dai_dt[i] = 2.0 * a0 * alpha_base * dbase_dt;
            d2ai_dt2[i] = 2.0 * a0
                * (dbase_dt * dbase_dt + alpha_base * d2base_dt2);
            bi[i] = 0.077_80 * r * tc / pc;
            if !ai[i].is_finite() || ai[i] <= 0.0
                || !dai_dt[i].is_finite() || !d2ai_dt2[i].is_finite()
                || !bi[i].is_finite() || bi[i] <= 0.0 {
                return Err(PengRobinsonErrorV1::Unrepresentable {
                    what: "pure-component attraction/covolume",
                });
            }
        }

        let mut aij = vec![0.0; n * n];
        let mut a_mix = 0.0;
        let mut da_mix_dt = 0.0;
        let mut d2a_mix_dt2 = 0.0;
        let mut b_mix = 0.0;
        for i in 0..n {
            b_mix += self.mole_fractions[i] * bi[i];
            for j in 0..n {
                let value = det::sqrt(ai[i] * ai[j])
                    * (1.0 - self.binary_interactions[i * n + j]);
                if !value.is_finite() || value < 0.0 {
                    return Err(PengRobinsonErrorV1::Unrepresentable {
                        what: "binary attraction coefficient",
                    });
                }
                let li = dai_dt[i] / ai[i];
                let lj = dai_dt[j] / ai[j];
                let dlog = 0.5 * (li + lj);
                let d2log = 0.5
                    * (d2ai_dt2[i] / ai[i] - li * li
                        + d2ai_dt2[j] / ai[j] - lj * lj);
                let dvalue = value * dlog;
                let d2value = value * (dlog * dlog + d2log);
                if !dvalue.is_finite() || !d2value.is_finite() {
                    return Err(PengRobinsonErrorV1::Unrepresentable {
                        what: "binary attraction temperature derivatives",
                    });
                }
                aij[i * n + j] = value;
                let weight = self.mole_fractions[i] * self.mole_fractions[j];
                a_mix += weight * value;
                da_mix_dt += weight * dvalue;
                d2a_mix_dt2 += weight * d2value;
            }
        }
        if !a_mix.is_finite() || a_mix <= 0.0
            || !da_mix_dt.is_finite() || !d2a_mix_dt2.is_finite()
            || !b_mix.is_finite() || b_mix <= 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "mixture attraction/covolume",
            });
        }

        Ok(MixedParameters { a_mix, da_mix_dt, d2a_mix_dt2, b_mix, aij, bi })
    }

    /// Evaluate PT roots and phase-local fugacity/density observables.
    pub fn evaluate_pt(
        &self,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<PengRobinsonEvaluationV1, PengRobinsonErrorV1> {
        let t = temperature.value();
        let p = pressure.value();
        if !t.is_finite() || t <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidTemperature { bits: t.to_bits() });
        }
        if !p.is_finite() || p <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidPressure { bits: p.to_bits() });
        }
        let mixed = self.mixing(t)?;
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let (dimensionless_a, dimensionless_b) = mixed.dimensionless(t, p)?;

        let rt = r * t;

        let c2 = -(1.0 - dimensionless_b);
        let c1 =
            dimensionless_a - 3.0 * dimensionless_b * dimensionless_b - 2.0 * dimensionless_b;
        let c0 = -(dimensionless_a * dimensionless_b
            - dimensionless_b * dimensionless_b
            - dimensionless_b * dimensionless_b * dimensionless_b);
        let mut roots = cubic_real_roots(c2, c1, c0)?;
        roots.retain(|z| {
            z.is_finite() && *z > dimensionless_b * (1.0 + ROOT_TOL)
        });
        roots.sort_by(|a, b| a.total_cmp(b));
        roots.dedup_by(|a, b| (*a - *b).abs() <= ROOT_TOL * a.abs().max(b.abs()));
        if roots.is_empty() {
            return Err(PengRobinsonErrorV1::NoMechanicalRoot);
        }

        let mut phases = Vec::with_capacity(roots.len());
        for z in roots {
            let molar_volume = z * rt / p;
            let density = self.molar_mass.value() / molar_volume;
            phases.push(mixed.phase_root(
                &self.mole_fractions, t, z, molar_volume, density,
                dimensionless_a, dimensionless_b,
            )?);
        }

        Ok(PengRobinsonEvaluationV1 {
            evaluator_version: PENG_ROBINSON_EVALUATOR_VERSION_V1,
            temperature,
            pressure,
            molar_mass: self.molar_mass,
            dimensionless_a,
            dimensionless_b,
            phases,
        })
    }

    /// Evaluate a homogeneous state at conserved mass density and temperature.
    ///
    /// This evaluates p(T,v) directly with v=M/rho, without a cubic root solve
    /// or a density-changing phase selection. It requires p>0, v>b and
    /// (dp/dv)_T<0. Mechanical stability includes metastable homogeneous
    /// states; it does not establish phase equilibrium or exclude coexistence.
    /// Caloric stability additionally requires the caller's ideal-reference Cv.
    pub fn evaluate_density_temperature(
        &self,
        temperature: Temperature,
        density: Density,
    ) -> Result<PengRobinsonDensityEvaluationV1, PengRobinsonErrorV1> {
        let t = temperature.value();
        let rho = density.value();
        if !rho.is_finite() || rho <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidDensity { bits: rho.to_bits() });
        }
        let mixed = self.mixing(t)?;
        let v = self.molar_mass.value() / rho;
        let b = mixed.b_mix;
        if !v.is_finite() || v <= b {
            return Err(PengRobinsonErrorV1::OutsideCovolume);
        }
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t;
        let denominator = v * v + 2.0 * b * v - b * b;
        let p = rt / (v - b) - mixed.a_mix / denominator;
        if !denominator.is_finite() || denominator <= 0.0 || !p.is_finite() || p <= 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "positive homogeneous density-temperature pressure",
            });
        }
        let z = p * v / rt;
        let (a, b) = mixed.dimensionless(t, p)?;
        let root = mixed.phase_root(&self.mole_fractions, t, z, v, rho, a, b)?;
        let (dp_dt, dp_dv) = root.pressure_derivatives();
        if !dp_dt.is_finite() || !dp_dv.is_finite() || dp_dv >= 0.0 {
            return Err(PengRobinsonErrorV1::MechanicallyUnstable);
        }
        Ok(PengRobinsonDensityEvaluationV1 {
            temperature, pressure: Pressure::new(p), molar_mass: self.molar_mass, root,
        })
    }
}

/// Shared temperature-dependent mixing law for both PT and density inputs.
struct MixedParameters {
    a_mix: f64,
    da_mix_dt: f64,
    d2a_mix_dt2: f64,
    b_mix: f64,
    aij: Vec<f64>,
    bi: Vec<f64>,
}
impl MixedParameters {
    fn dimensionless(&self, t: f64, p: f64) -> Result<(f64, f64), PengRobinsonErrorV1> {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t;
        let a = self.a_mix * p / (rt * rt);
        let b = self.b_mix * p / rt;
        if !a.is_finite() || !b.is_finite() || a <= 0.0 || b <= 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "dimensionless Peng-Robinson coefficients",
            });
        }
        Ok((a, b))
    }

    fn phase_root(&self, fractions: &[f64], t: f64, z: f64, volume: f64,
        density: f64, dimensionless_a: f64, dimensionless_b: f64)
        -> Result<PengRobinsonPhaseRootV1, PengRobinsonErrorV1> {
        if !volume.is_finite() || volume <= self.b_mix
            || !density.is_finite() || density <= 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable { what: "phase-root observables" });
        }
        let residual_geometry = ResidualGeometry::new(
            z, dimensionless_a, dimensionless_b, volume, self.b_mix,
        )?;
        let log_phi = fugacity_logs(residual_geometry, t,
            self.a_mix, self.b_mix, fractions, &self.aij, &self.bi)?;
        let residual_gibbs_over_rt = fractions.iter().zip(&log_phi)
            .map(|(x, value)| x * value).sum::<f64>();
        if !residual_gibbs_over_rt.is_finite() {
            return Err(PengRobinsonErrorV1::Unrepresentable { what: "phase-root observables" });
        }
        Ok(PengRobinsonPhaseRootV1 {
            compressibility_factor: z,
            molar_volume_m3_per_mol: volume,
            mass_density_kg_per_m3: density,
            residual_gibbs_over_rt,
            log_fugacity_coefficients: log_phi,
            attraction_parameter: self.a_mix,
            attraction_temperature_derivative: self.da_mix_dt,
            attraction_second_temperature_derivative: self.d2a_mix_dt2,
            covolume_m3_per_mol: self.b_mix,
            temperature_k: t,
            residual_geometry,
        })
    }
}

/// Direct homogeneous state at one exact supplied density; no phase switch.
#[derive(Debug, Clone, PartialEq)]
pub struct PengRobinsonDensityEvaluationV1 {
    temperature: Temperature,
    pressure: Pressure,
    molar_mass: MolarMass,
    root: PengRobinsonPhaseRootV1,
}
impl PengRobinsonDensityEvaluationV1 {
    #[must_use]
    pub const fn temperature(&self) -> Temperature { self.temperature }
    #[must_use]
    pub const fn pressure(&self) -> Pressure { self.pressure }
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass { self.molar_mass }
    #[must_use]
    pub const fn root(&self) -> &PengRobinsonPhaseRootV1 { &self.root }
}

fn cubic_real_roots(
    c2: f64,
    c1: f64,
    c0: f64,
) -> Result<Vec<f64>, PengRobinsonErrorV1> {
    let p = c1 - c2 * c2 / 3.0;
    let q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    let discriminant = 0.25 * q * q + p * p * p / 27.0;
    if !p.is_finite() || !q.is_finite() || !discriminant.is_finite() {
        return Err(PengRobinsonErrorV1::Unrepresentable {
            what: "cubic reduction",
        });
    }

    let scale = (0.25 * q * q)
        .abs()
        .max((p * p * p / 27.0).abs())
        .max(f64::MIN_POSITIVE);
    if discriminant > ROOT_TOL * scale {
        let root_discriminant = det::sqrt(discriminant);
        // Recover the smaller Cardano term from u*v = -p/3. Directly
        // subtracting nearly equal q/2 and sqrt(discriminant) can round its
        // argument to zero; taking a cube root amplifies that loss into large
        // fugacity noise even when this unique mechanical root is well conditioned.
        let u = signed_cuberoot(-0.5 * q - root_discriminant.copysign(q));
        if !u.is_finite() || u == 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "single-root cubic factor",
            });
        }
        let root = polish_cubic_root(u - p / (3.0 * u) - c2 / 3.0, c2, c1, c0);
        return Ok(vec![root]);
    }
    if p >= 0.0 {
        return Ok(vec![-c2 / 3.0]);
    }

    let radius = 2.0 * det::sqrt(-p / 3.0);
    let denominator = det::sqrt(-(p * p * p) / 27.0);
    if !radius.is_finite() || !denominator.is_finite() || denominator <= 0.0 {
        return Err(PengRobinsonErrorV1::Unrepresentable {
            what: "three-root cubic geometry",
        });
    }
    let theta = det::acos((-0.5 * q / denominator).clamp(-1.0, 1.0)) / 3.0;
    let candidates = [
        radius * det::cos(theta) - c2 / 3.0,
        radius * det::cos(theta - 2.0 * core::f64::consts::PI / 3.0) - c2 / 3.0,
        radius * det::cos(theta - 4.0 * core::f64::consts::PI / 3.0) - c2 / 3.0,
    ];
    // At small pressure the gas root is O(1), while both dense roots are O(p).
    // The depressed discriminant subtracts O(1) terms and can lose their sign;
    // the two small trigonometric roots then lose their volumes or even invent
    // a liquid branch above the critical temperature. Use only the best
    // separated candidate (largest |f'|), polish it in the ORIGINAL cubic,
    // and recover the remaining quadratic through Vieta's identities.
    let root = candidates.into_iter().max_by(|a, b| {
        let slope = |z: f64| ((3.0 * z + 2.0 * c2) * z + c1).abs();
        slope(*a).total_cmp(&slope(*b))
    }).expect("three finite cubic candidates");
    let root = polish_cubic_root(root, c2, c1, c0);
    let (sum, product) = if root == 0.0 {
        if c0 != 0.0 {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "cubic deflation root",
            });
        }
        (-c2, c1)
    } else {
        let product = -c0 / root;
        let direct_sum = -c2 - root;
        let sum = if direct_sum.abs() > 0.5 * c2.abs().max(root.abs()) {
            direct_sum
        } else {
            // r*(s+t)+s*t=c1 retains two small roots when r nearly cancels c2.
            (c1 - product) / root
        };
        (sum, product)
    };
    if !sum.is_finite() || !product.is_finite() {
        return Err(PengRobinsonErrorV1::Unrepresentable {
            what: "cubic deflation coefficients",
        });
    }
    let scale = sum.abs().max(det::sqrt(product.abs()));
    if scale == 0.0 {
        return Ok(vec![root, 0.0, 0.0]);
    }
    let scaled_sum = sum / scale;
    let scaled_product = (product / scale) / scale;
    let discriminant = scaled_sum * scaled_sum - 4.0 * scaled_product;
    let allowance = ROOT_TOL * (scaled_sum * scaled_sum).max(4.0 * scaled_product.abs());
    if discriminant < -allowance {
        return Ok(vec![root]);
    }
    // A tiny negative value at a repeated root is roundoff. This allowance is
    // relative to the remaining roots' OWN scale, never to an O(1) gas root.
    let delta = scale * det::sqrt(discriminant.max(0.0));
    let first = 0.5 * (sum + delta.copysign(sum));
    let second = if first == 0.0 { 0.0 } else { product / first };
    Ok(vec![root, first, second])
}

fn polish_cubic_root(mut root: f64, c2: f64, c1: f64, c0: f64) -> f64 {
    // Fixed work; accept only a finite correction that improves the original
    // polynomial residual. No tolerance or constitutive coefficient is changed.
    for _ in 0..2 {
        let residual = ((root + c2) * root + c1) * root + c0;
        let slope = (3.0 * root + 2.0 * c2) * root + c1;
        let next = root - residual / slope;
        if !next.is_finite() || !slope.is_finite() || slope == 0.0
            || (((next + c2) * next + c1) * next + c0).abs() >= residual.abs() {
            break;
        }
        root = next;
    }
    root
}

fn signed_cuberoot(value: f64) -> f64 {
    if value == 0.0 {
        0.0
    } else {
        value.signum() * det::exp(det::ln(value.abs()) / 3.0)
    }
}

/// Shared original-EOS geometry for fugacity and all caloric departures.
/// Retain offsets before adding them to one; a rounded Z=1 is not zero
/// residual energy or fugacity. Immutable roots reuse the logarithms and integral.
#[derive(Debug, Clone, Copy, PartialEq)]
struct ResidualGeometry {
    z_minus_one: f64,
    log_z: f64,
    log_free_volume: f64,
    attraction_integral: f64,
}

impl ResidualGeometry {
    fn new(z: f64, a: f64, b: f64, volume: f64, covolume: f64)
        -> Result<Self, PengRobinsonErrorV1> {
        let mut z_minus_one = z - 1.0;
        if z_minus_one.abs() < 0.5 {
            // The original cubic gives Z-1 = B-A*(Z-B)/(Z²+2BZ-B²).
            // Unlike subtraction from rounded Z, this retains its dilute limit.
            let denominator = z * z + 2.0 * b * z - b * b;
            z_minus_one = b - a * ((z - b) / denominator);
        }
        let log_z = if z_minus_one.abs() < 0.5 {
            log_one_plus(z_minus_one)
        } else {
            // Dense low-pressure roots can have Z so small that Z-1 rounds
            // to -1. Their logarithm must retain the original positive root.
            det::ln(z)
        };
        let eta = covolume / volume;
        let sqrt_2 = det::sqrt(2.0);
        let log_free_volume = log_one_plus(-eta);
        let attraction_integral = log_one_plus(
            2.0 * sqrt_2 * eta / (1.0 + (1.0 - sqrt_2) * eta),
        ) / (2.0 * sqrt_2 * covolume);
        if [z_minus_one, log_z, log_free_volume, attraction_integral]
            .iter().any(|value| !value.is_finite()) {
            return Err(PengRobinsonErrorV1::Unrepresentable { what: "residual EOS geometry" });
        }
        Ok(Self { z_minus_one, log_z, log_free_volume, attraction_integral })
    }
}

/// Kahan's log1p identity using the same deterministic logarithm as the EOS.
fn log_one_plus(offset: f64) -> f64 {
    let rounded = 1.0 + offset;
    if rounded == 1.0 {
        offset
    } else {
        det::ln(rounded) * (offset / (rounded - 1.0))
    }
}

#[allow(clippy::too_many_arguments)]
fn fugacity_logs(
    geometry: ResidualGeometry,
    temperature: f64,
    a_mix: f64,
    b_mix: f64,
    x: &[f64],
    aij: &[f64],
    bi: &[f64],
) -> Result<Vec<f64>, PengRobinsonErrorV1> {
    let attraction_scale = a_mix / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature)
        * geometry.attraction_integral;
    let n = x.len();
    let mut result = vec![0.0; n];
    for i in 0..n {
        let sum_aij = (0..n).map(|j| x[j] * aij[i * n + j]).sum::<f64>();
        let bi_over_b = bi[i] / b_mix;
        let bracket = 2.0 * sum_aij / a_mix - bi_over_b;
        let value = bi_over_b * geometry.z_minus_one
            - geometry.log_z - geometry.log_free_volume
            - attraction_scale * bracket;
        if !value.is_finite() {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "fugacity coefficient",
            });
        }
        result[i] = value;
    }
    Ok(result)
}

/// One admissible real root with phase-local mechanical/fugacity observables.
#[derive(Debug, Clone, PartialEq)]
pub struct PengRobinsonPhaseRootV1 {
    compressibility_factor: f64,
    molar_volume_m3_per_mol: f64,
    mass_density_kg_per_m3: f64,
    residual_gibbs_over_rt: f64,
    log_fugacity_coefficients: Vec<f64>,
    attraction_parameter: f64,
    attraction_temperature_derivative: f64,
    attraction_second_temperature_derivative: f64,
    covolume_m3_per_mol: f64,
    temperature_k: f64,
    residual_geometry: ResidualGeometry,
}

impl PengRobinsonPhaseRootV1 {
    #[must_use]
    pub const fn compressibility_factor(&self) -> f64 {
        self.compressibility_factor
    }
    #[must_use]
    pub const fn molar_volume_m3_per_mol(&self) -> f64 {
        self.molar_volume_m3_per_mol
    }
    #[must_use]
    pub const fn mass_density_kg_per_m3(&self) -> f64 {
        self.mass_density_kg_per_m3
    }
    #[must_use]
    pub const fn residual_gibbs_over_rt(&self) -> f64 {
        self.residual_gibbs_over_rt
    }
    #[must_use]
    pub fn log_fugacity_coefficients(&self) -> &[f64] {
        &self.log_fugacity_coefficients
    }

    /// Residual internal energy U^R [J/mol] relative to the ideal-gas caloric
    /// reference at the same T and composition.
    #[must_use]
    pub fn residual_internal_energy_j_mol(&self) -> f64 {
        (self.temperature_k * self.attraction_temperature_derivative
            - self.attraction_parameter)
            * self.residual_geometry.attraction_integral
    }

    /// Residual enthalpy H^R [J/mol].
    #[must_use]
    pub fn residual_enthalpy_j_mol(&self) -> f64 {
        UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * self.temperature_k
            * self.residual_geometry.z_minus_one
            + self.residual_internal_energy_j_mol()
    }

    /// Entropy departure S^R [J/(mol K)] from an ideal gas at the SAME T,
    /// density and frozen composition. Differentiating the PR residual
    /// Helmholtz energy gives `R ln(1-b/v) + a'(T) I(v)`.
    /// This reference is useful for conservative density/entropy changes;
    /// it differs from the conventional same-pressure departure by `R ln Z`.
    #[must_use]
    pub fn residual_entropy_at_density_j_mol_k(&self) -> f64 {
        UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * self.residual_geometry.log_free_volume
            + self.attraction_temperature_derivative * self.residual_geometry.attraction_integral
    }

    /// Entropy departure S^R [J/(mol K)] from the ideal gas at the SAME T,
    /// pressure and composition, consistent with `H^R/T - R*(G^R/RT)`.
    /// Original ideal-reference and mixing entropies are supplied by the
    /// thermochemical species owner; no phase-equilibrium entropy is inferred.
    #[must_use]
    pub fn residual_entropy_j_mol_k(&self) -> f64 {
        self.residual_entropy_at_density_j_mol_k()
            + UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * self.residual_geometry.log_z
    }

    /// Residual constant-volume heat capacity C_v^R [J/(mol K)].
    #[must_use]
    pub fn residual_cv_j_mol_k(&self) -> f64 {
        self.temperature_k * self.attraction_second_temperature_derivative
            * self.residual_geometry.attraction_integral
    }

    /// (dp/dT)_v [Pa/K] and (dp/dv)_T [Pa mol/m^3] from the original PR law.
    #[must_use]
    pub fn pressure_derivatives(&self) -> (f64, f64) {
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let t = self.temperature_k;
        let v = self.molar_volume_m3_per_mol;
        let b = self.covolume_m3_per_mol;
        let d = v * v + 2.0 * b * v - b * b;
        let dp_dt_v = r / (v - b) - self.attraction_temperature_derivative / d;
        let dp_dv_t = -r * t / ((v - b) * (v - b))
            + self.attraction_parameter * (2.0 * v + 2.0 * b) / (d * d);
        (dp_dt_v, dp_dv_t)
    }

    /// Complete Cp and isentropic sound speed after adding an ideal-reference
    /// molar Cv from the SAME composition and T. Returns None for a mechanically
    /// unstable root or nonpositive complete caloric response.
    #[must_use]
    pub fn caloric_mechanics(
        &self,
        ideal_molar_cv_j_mol_k: f64,
        molar_mass_kg_mol: f64,
    ) -> Option<PengRobinsonCaloricMechanicsV1> {
        let cv = ideal_molar_cv_j_mol_k + self.residual_cv_j_mol_k();
        let (dp_dt_v, dp_dv_t) = self.pressure_derivatives();
        if !cv.is_finite() || cv <= 0.0
            || !molar_mass_kg_mol.is_finite() || molar_mass_kg_mol <= 0.0
            || !dp_dt_v.is_finite() || !dp_dv_t.is_finite() || dp_dv_t >= 0.0
        {
            return None;
        }
        let t = self.temperature_k;
        let v = self.molar_volume_m3_per_mol;
        let cp = cv - t * dp_dt_v * dp_dt_v / dp_dv_t;
        let sound_speed_squared = v * v / molar_mass_kg_mol
            * (-dp_dv_t + t * dp_dt_v * dp_dt_v / cv);
        if !cp.is_finite() || cp <= cv
            || !sound_speed_squared.is_finite() || sound_speed_squared <= 0.0
        {
            return None;
        }
        Some(PengRobinsonCaloricMechanicsV1 {
            molar_cv_j_mol_k: cv,
            molar_cp_j_mol_k: cp,
            sound_speed_m_s: det::sqrt(sound_speed_squared),
            dp_dt_v_pa_k: dp_dt_v,
            dp_dv_t_pa_mol_m3: dp_dv_t,
        })
    }
}

/// Complete local caloric/mechanical derivatives for one stable PR root.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PengRobinsonCaloricMechanicsV1 {
    pub molar_cv_j_mol_k: f64,
    pub molar_cp_j_mol_k: f64,
    pub sound_speed_m_s: f64,
    pub dp_dt_v_pa_k: f64,
    pub dp_dv_t_pa_mol_m3: f64,
}

/// Complete fixed-composition PT evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct PengRobinsonEvaluationV1 {
    evaluator_version: u32,
    temperature: Temperature,
    pressure: Pressure,
    molar_mass: MolarMass,
    dimensionless_a: f64,
    dimensionless_b: f64,
    phases: Vec<PengRobinsonPhaseRootV1>,
}

impl PengRobinsonEvaluationV1 {
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }
    #[must_use]
    pub const fn dimensionless_a(&self) -> f64 {
        self.dimensionless_a
    }
    #[must_use]
    pub const fn dimensionless_b(&self) -> f64 {
        self.dimensionless_b
    }
    #[must_use]
    pub fn phase_roots(&self) -> &[PengRobinsonPhaseRootV1] {
        &self.phases
    }

    /// Select an explicit mechanical branch, or the minimum residual-Gibbs
    /// branch at this same fixed composition. This is not a multiphase flash.
    #[must_use]
    pub fn select(&self, phase: PengRobinsonPhaseV1) -> &PengRobinsonPhaseRootV1 {
        match phase {
            PengRobinsonPhaseV1::Liquid => &self.phases[0],
            PengRobinsonPhaseV1::Vapor => &self.phases[self.phases.len() - 1],
            PengRobinsonPhaseV1::Stable => self
                .phases
                .iter()
                .min_by(|a, b| {
                    a.residual_gibbs_over_rt
                        .total_cmp(&b.residual_gibbs_over_rt)
                })
                .expect("successful evaluation retains a phase root"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn methane() -> PengRobinsonComponentV1 {
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
            molar_mass: MolarMass::new(0.016_042_46),
        }
    }

    fn pure(component: PengRobinsonComponentV1) -> PengRobinsonMixtureV1 {
        PengRobinsonMixtureV1::new(vec![component], vec![1.0], vec![]).unwrap()
    }

    fn methane_ethane(methane_fraction: f64) -> PengRobinsonMixtureV1 {
        let ethane = PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(305.322),
            critical_pressure: Pressure::new(4.8722e6),
            acentric_factor: 0.0995,
            molar_mass: MolarMass::new(0.030_069_04),
        };
        PengRobinsonMixtureV1::new(vec![methane(), ethane],
            vec![methane_fraction, 1.0 - methane_fraction], vec![]).unwrap()
    }

    #[test]
    fn g1_single_root_cardano_retains_small_complement_and_cubic_residual() {
        // This original liquid-composition state has depressed p=5.58e-7,
        // q=0.01389. One direct Cardano argument rounds to zero, giving a
        // cubic residual about1.34e-7 before the product-form correction.
        let state = methane_ethane(0.314_603_967_109_088_94)
            .evaluate_pt(Temperature::new(217.152_272_470_045_6),
                Pressure::new(2_500_001.521_079_547_7)).unwrap();
        assert_eq!(state.phase_roots().len(), 1);
        assert_eq!(state.evaluator_version(), 4);
        let z = state.phase_roots()[0].compressibility_factor();
        // Independent polynomial root, not a cubic-formula round trip.
        assert!((z - 0.076_224_432_883_693_82).abs() < 2e-14);
        let a = state.dimensionless_a();
        let b = state.dimensionless_b();
        let c2 = b - 1.0;
        let c1 = a - 3.0*b*b - 2.0*b;
        let c0 = -(a*b - b*b - b*b*b);
        assert!((((z + c2)*z + c1)*z + c0).abs() < 1e-16);
        for p in [-1e-6, -1e-9, 0.0, 1e-9, 1e-6] {
            let q = 0.013_892_238_324_066_268;
            let roots = cubic_real_roots(0.0, p, q).unwrap();
            assert_eq!(roots.len(), 1);
            let z = roots[0];
            assert!((z*z*z + p*z + q).abs() < 2e-17);
        }
    }

    #[test]
    fn g1_binary_flash_survives_cardano_cancellation_without_weaker_tolerances() {
        let model = methane_ethane(0.5);
        let options = phase_screen::PengRobinsonPhaseScreenOptionsV1 {
            stability: stability::PengRobinsonStabilityOptionsV1 {
                stationarity_tolerance: 1e-12, tangent_plane_tolerance: 1e-12,
                ..Default::default()
            },
            flash: flash::PengRobinsonFlashOptionsV1 {
                log_k_tolerance: 1e-12, ..Default::default()
            },
        };
        let state = model.phase_screen_pt_cancellable(
            Temperature::new(217.152_272_470_045_6),
            Pressure::new(2_500_001.521_079_547_7), options, || false).unwrap();
        let split = state.two_phase().unwrap();
        // Independent simultaneous polynomial/fugacity equations and exact
        // binary lever rule, not production successive substitution.
        assert!((split.liquid_mole_fractions()[0] - 0.314_603_967_109_088_94).abs() < 3e-10);
        assert!((split.vapor_mole_fractions()[0] - 0.803_195_924_585_061_2).abs() < 3e-10);
        assert!((split.vapor_fraction_mol() - 0.379_449_620_596_811_37).abs() < 3e-10);
        assert!(split.log_k_residual() <= options.flash.log_k_tolerance);
        assert!(state.stability().vapor_seed().stationarity_residual() <= 1e-12);
        assert!(state.stability().liquid_seed().stationarity_residual() <= 1e-12);
        assert!(split.iterations() <= options.flash.maximum_iterations);
    }

    #[test]
    fn supercritical_methane_matches_independent_cubic_root() {
        let state = pure(methane())
            .evaluate_pt(Temperature::new(300.0), Pressure::new(5.0e6))
            .unwrap();
        assert_eq!(state.phase_roots().len(), 1);
        let root = state.select(PengRobinsonPhaseV1::Stable);
        assert!((root.compressibility_factor() - 0.901_830_710_952_881_1).abs() < 2e-13);
        assert!((root.mass_density_kg_per_m3() - 35.658_29).abs() < 2e-3);
        assert_eq!(state.select(PengRobinsonPhaseV1::Liquid), root);
        assert_eq!(state.select(PengRobinsonPhaseV1::Vapor), root);
    }

    #[test]
    fn subcritical_methane_retains_all_three_roots() {
        let state = pure(methane())
            .evaluate_pt(Temperature::new(150.0), Pressure::new(1.0e6))
            .unwrap();
        let roots: Vec<_> = state
            .phase_roots()
            .iter()
            .map(PengRobinsonPhaseRootV1::compressibility_factor)
            .collect();
        assert_eq!(roots.len(), 3);
        for (actual, expected) in roots.iter().zip([
            0.033_118_329_926_412_89,
            0.120_348_510_936_658_85,
            0.825_042_640_206_099_7,
        ]) {
            assert!((*actual - expected).abs() < 3e-13);
        }
        assert_eq!(
            state.select(PengRobinsonPhaseV1::Liquid).compressibility_factor(),
            roots[0]
        );
        assert_eq!(
            state.select(PengRobinsonPhaseV1::Vapor).compressibility_factor(),
            roots[2]
        );
        assert!(state
            .phase_roots()
            .iter()
            .all(|root| root.log_fugacity_coefficients().len() == 1));
    }

    #[test]
    fn dilute_limit_recovers_ideal_gas_mechanics() {
        let state = pure(methane())
            .evaluate_pt(Temperature::new(400.0), Pressure::new(1.0))
            .unwrap();
        let root = state.select(PengRobinsonPhaseV1::Stable);
        assert!((root.compressibility_factor() - 1.0).abs() < 1e-7);
        let ideal_density = methane().molar_mass.value()
            / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0);
        assert!((root.mass_density_kg_per_m3() / ideal_density - 1.0).abs() < 1e-7);
        let direct = pure(methane()).evaluate_density_temperature(
            Temperature::new(400.0), Density::new(ideal_density),
        ).unwrap();
        assert!((direct.pressure().value() - 1.0).abs() < 1e-7);
        assert!((direct.root().compressibility_factor() - 1.0).abs() < 1e-7);
    }

    #[test]
    fn g1_low_pressure_cubic_preserves_liquid_volume_and_all_original_rows() {
        let model = pure(methane());
        // Independent 70-digit roots of the original dimensional volume
        // polynomial, with the same declared methane constants at 150 K.
        // These are mechanical/metastable branches, not phase-stability claims.
        let references = [
            (100.0, 0.000_041_760_662_761_526_087),
            (1.0, 0.000_041_760_710_558_101_24),
            (0.01, 0.000_041_760_711_036_069_79),
            (0.0001, 0.000_041_760_711_040_849_47),
            (0.000001, 0.000_041_760_711_040_897_27),
            (0.00000001, 0.000_041_760_711_040_897_75),
        ];
        for (pressure, expected_volume) in references {
            let state = model
                .evaluate_pt(Temperature::new(150.0), Pressure::new(pressure))
                .unwrap();
            assert_eq!(state.phase_roots().len(), 3, "p={pressure} Pa");
            let liquid = state.select(PengRobinsonPhaseV1::Liquid);
            assert!(
                (liquid.molar_volume_m3_per_mol() / expected_volume - 1.0).abs() < 5e-13,
                "p={pressure} Pa, v={}", liquid.molar_volume_m3_per_mol(),
            );
            let a = state.dimensionless_a();
            let b = state.dimensionless_b();
            let c2 = b - 1.0;
            let c1 = a - 3.0 * b * b - 2.0 * b;
            let c0 = -(a * b - b * b - b * b * b);
            for root in state.phase_roots() {
                let z = root.compressibility_factor();
                let residual = ((z + c2) * z + c1) * z + c0;
                let scale = (z * z * z).abs() + (c2 * z * z).abs()
                    + (c1 * z).abs() + c0.abs();
                assert!(residual.abs() < 512.0 * f64::EPSILON * scale,
                    "p={pressure} Pa, z={z}, relative residual={}", residual / scale);
            }
            // The intermediate mechanical branch remains unstable; repairing
            // its root never turns it into an admissible acoustic material.
            assert!(liquid.pressure_derivatives().1 < 0.0);
            assert!(state.phase_roots()[1].pressure_derivatives().1 > 0.0);
            assert!(state.select(PengRobinsonPhaseV1::Vapor).pressure_derivatives().1 < 0.0);
        }
    }

    #[test]
    fn g1_low_pressure_supercritical_gas_has_no_phantom_dense_root() {
        let model = pure(methane());
        for temperature in [200.0, 250.0, 300.0, 400.0] {
            for pressure in [100.0, 1.0, 0.01, 0.0001, 0.000001] {
                let state = model.evaluate_pt(
                    Temperature::new(temperature), Pressure::new(pressure),
                ).unwrap();
                // T exceeds the supplied critical temperature. The original
                // PR polynomial has one real root throughout this pressure set.
                assert_eq!(state.phase_roots().len(), 1,
                    "T={temperature} K, p={pressure} Pa: {:?}", state.phase_roots());
                let vapor = state.select(PengRobinsonPhaseV1::Vapor);
                assert_eq!(state.select(PengRobinsonPhaseV1::Liquid), vapor);
                assert_eq!(state.select(PengRobinsonPhaseV1::Stable), vapor);
                let ideal = methane().molar_mass.value() * pressure
                    / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature);
                assert!((vapor.mass_density_kg_per_m3() / ideal - 1.0).abs() < 1e-5);
                assert!(vapor.caloric_mechanics(20.0, methane().molar_mass.value()).is_some());
            }
        }
    }

    #[test]
    fn g1_cubic_deflation_retains_small_and_repeated_roots() {
        for expected in [
            [0.0, 1.0, 1.0],
            [1e-20, 1.0, 1.0],
            [-2.0, 1.0, 1.0],
            [-1e-8, 1e-8, 1.0],
            [1e-14, 3e-14, 1.0],
            [-1.0, 0.0, 1.0],
        ] {
            let [a, b, c] = expected;
            let mut roots = cubic_real_roots(-(a + b + c), a*b + a*c + b*c, -a*b*c).unwrap();
            roots.sort_by(f64::total_cmp);
            assert_eq!(roots.len(), 3, "expected={expected:?}, actual={roots:?}");
            for (&actual, reference) in roots.iter().zip(expected) {
                if reference == 0.0 {
                    assert_eq!(actual, 0.0);
                } else {
                    // Repeated roots are intrinsically ill-conditioned; the
                    // distinct small roots still retain their own relative scale.
                    let repeated = expected.windows(2).any(|p| p[0] == p[1]);
                    let tolerance = if repeated && reference.abs() >= 1.0 { 1e-7 } else { 2e-13 };
                    assert!((actual / reference - 1.0).abs() < tolerance,
                        "expected={expected:?}, actual={roots:?}");
                }
            }
        }
    }

    #[test]
    fn binary_mixing_is_permutation_invariant_and_component_resolved() {
        let co2 = PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(304.1282),
            critical_pressure: Pressure::new(7.3773e6),
            acentric_factor: 0.22394,
            molar_mass: MolarMass::new(0.044_009_5),
        };
        let a = PengRobinsonMixtureV1::new(
            vec![methane(), co2],
            vec![0.7, 0.3],
            vec![0.0, 0.09, 0.09, 0.0],
        )
        .unwrap();
        let b = PengRobinsonMixtureV1::new(
            vec![co2, methane()],
            vec![0.3, 0.7],
            vec![0.0, 0.09, 0.09, 0.0],
        )
        .unwrap();
        let ea = a
            .evaluate_pt(Temperature::new(280.0), Pressure::new(8.0e6))
            .unwrap();
        let eb = b
            .evaluate_pt(Temperature::new(280.0), Pressure::new(8.0e6))
            .unwrap();
        let ra = ea.select(PengRobinsonPhaseV1::Stable);
        let rb = eb.select(PengRobinsonPhaseV1::Stable);
        assert!((ra.compressibility_factor() - rb.compressibility_factor()).abs() < 1e-14);
        let la = ra.log_fugacity_coefficients();
        let lb = rb.log_fugacity_coefficients();
        assert!((la[0] - lb[1]).abs() < 1e-14);
        assert!((la[1] - lb[0]).abs() < 1e-14);
    }

    #[test]
    fn residual_calorics_and_sound_speed_match_independent_differentiation() {
        let model = pure(methane());
        let t = 300.0;
        let p = 5.0e6;
        let state = model.evaluate_pt(Temperature::new(t), Pressure::new(p)).unwrap();
        let root = state.select(PengRobinsonPhaseV1::Stable);
        let ideal_cv = 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let mechanics = root.caloric_mechanics(ideal_cv, methane().molar_mass.value()).unwrap();

        let dt = 1.0e-3;
        let plus = model.evaluate_pt(Temperature::new(t + dt), Pressure::new(p)).unwrap();
        let minus = model.evaluate_pt(Temperature::new(t - dt), Pressure::new(p)).unwrap();
        let hp = plus.select(PengRobinsonPhaseV1::Stable).residual_enthalpy_j_mol();
        let hm = minus.select(PengRobinsonPhaseV1::Stable).residual_enthalpy_j_mol();
        let cp_departure_numeric = (hp - hm) / (2.0 * dt);
        let cp_departure = mechanics.molar_cp_j_mol_k
            - (ideal_cv + UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K);
        assert!((cp_departure - cp_departure_numeric).abs()
            <= 2.0e-4 * cp_departure_numeric.abs().max(1.0));

        let v = root.molar_volume_m3_per_mol();
        let dv = v * 1.0e-6;
        let pressure_from_v = |volume: f64| {
            let a = root.attraction_parameter;
            let b = root.covolume_m3_per_mol;
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t / (volume - b)
                - a / (volume * volume + 2.0 * b * volume - b * b)
        };
        let dpdv_numeric = (pressure_from_v(v + dv) - pressure_from_v(v - dv)) / (2.0 * dv);
        assert!((mechanics.dp_dv_t_pa_mol_m3 / dpdv_numeric - 1.0).abs() < 2.0e-10);
        assert!(mechanics.sound_speed_m_s.is_finite() && mechanics.sound_speed_m_s > 0.0);
        assert!(root.residual_internal_energy_j_mol().is_finite());
        assert!(root.residual_cv_j_mol_k().is_finite());
    }

    #[test]
    fn g1_density_temperature_round_trips_each_stable_homogeneous_root() {
        let model = pure(methane());
        for (temperature, pressure) in [(150.0, 1.0e6), (300.0, 5.0e6), (500.0, 1.0e5)] {
            let pt = model.evaluate_pt(Temperature::new(temperature), Pressure::new(pressure)).unwrap();
            for root in pt.phase_roots() {
                if root.pressure_derivatives().1 >= 0.0 { continue; }
                let density = root.mass_density_kg_per_m3();
                let direct = model.evaluate_density_temperature(
                    Temperature::new(temperature), Density::new(density),
                ).unwrap();
                assert_eq!(direct.root().mass_density_kg_per_m3().to_bits(), density.to_bits());
                assert!((direct.pressure().value() / pressure - 1.0).abs() < 2.0e-10);
                assert!((direct.root().residual_internal_energy_j_mol()
                    - root.residual_internal_energy_j_mol()).abs() < 1.0e-9);
            }
        }
    }

    #[test]
    fn g1_density_tangent_is_actual_constant_volume_residual_cv() {
        let model = pure(methane());
        let temperature = 300.0;
        let density = Density::new(50.0);
        let center = model.evaluate_density_temperature(Temperature::new(temperature), density).unwrap();
        let h = 1.0e-3;
        let minus = model.evaluate_density_temperature(Temperature::new(temperature - h), density).unwrap();
        let plus = model.evaluate_density_temperature(Temperature::new(temperature + h), density).unwrap();
        let cv = (plus.root().residual_internal_energy_j_mol()
            - minus.root().residual_internal_energy_j_mol()) / (2.0 * h);
        assert!((cv / center.root().residual_cv_j_mol_k() - 1.0).abs() < 2.0e-8);
        let dpdt = (plus.pressure().value() - minus.pressure().value()) / (2.0 * h);
        assert!((dpdt / center.root().pressure_derivatives().0 - 1.0).abs() < 2.0e-9);
    }

    #[test]
    fn g1_entropy_departure_closes_gibbs_caloric_and_maxwell_identities() {
        let model = pure(methane());
        let t = 300.0;
        let rho = 50.0;
        let evaluate = |temperature, density| model.evaluate_density_temperature(
            Temperature::new(temperature), Density::new(density),
        ).unwrap();
        let state = evaluate(t, rho);
        let root = state.root();
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let from_gibbs = root.residual_enthalpy_j_mol() / t
            - r * root.residual_gibbs_over_rt();
        assert!((root.residual_entropy_j_mol_k() - from_gibbs).abs() < 2.0e-12);
        let dt = 1.0e-3;
        let ds_dt = (evaluate(t + dt, rho).root().residual_entropy_at_density_j_mol_k()
            - evaluate(t - dt, rho).root().residual_entropy_at_density_j_mol_k()) / (2.0 * dt);
        assert!((ds_dt * t / root.residual_cv_j_mol_k() - 1.0).abs() < 2.0e-8);
        let v = root.molar_volume_m3_per_mol();
        let dv = v * 1.0e-5;
        let m = model.molar_mass().value();
        let ds_dv = (evaluate(t, m / (v + dv)).root().residual_entropy_at_density_j_mol_k()
            - evaluate(t, m / (v - dv)).root().residual_entropy_at_density_j_mol_k()) / (2.0 * dv)
            + r / v;
        assert!((ds_dv / root.pressure_derivatives().0 - 1.0).abs() < 2.0e-8);
    }

    #[test]
    fn g1_dilute_departures_retain_original_second_virial_response() {
        let model = pure(methane());
        let t = 300.0;
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        for pressure in [1.0e-6, 1.0e-9, 1.0e-12] {
            let state = model.evaluate_pt(Temperature::new(t), Pressure::new(pressure)).unwrap();
            let root = state.select(PengRobinsonPhaseV1::Vapor);
            let v = root.molar_volume_m3_per_mol();
            let a = root.attraction_parameter;
            let da = root.attraction_temperature_derivative;
            let d2a = root.attraction_second_temperature_derivative;
            let b = root.covolume_m3_per_mol;
            // Independently expand the original Helmholtz law to order 1/v.
            // These departures remain representable after rounded Z becomes 1.
            for (name, actual, virial) in [
                ("u", root.residual_internal_energy_j_mol(), (t * da - a) / v),
                ("h", root.residual_enthalpy_j_mol(), (t * da - 2.0 * a + r * t * b) / v),
                ("cv", root.residual_cv_j_mol_k(), t * d2a / v),
                ("s(rho)", root.residual_entropy_at_density_j_mol_k(), (da - r * b) / v),
                ("s(p)", root.residual_entropy_j_mol_k(), (da - a / t) / v),
                ("ln(phi)", root.log_fugacity_coefficients()[0], (b - a / (r * t)) / v),
            ] {
                assert!((actual / virial - 1.0).abs() < 2.0e-10,
                    "p={pressure}, {name}: actual={actual:e}, virial={virial:e}");
            }
            let from_gibbs = root.residual_enthalpy_j_mol() / t
                - r * root.residual_gibbs_over_rt();
            assert!((root.residual_entropy_j_mol_k() / from_gibbs - 1.0).abs() < 1.0e-13);
        }
    }

    #[test]
    fn g1_dilute_binary_fugacity_respects_each_partial_virial() {
        let model = methane_ethane(0.3);
        let t = 300.0;
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t;
        let mixed = model.mixing(t).unwrap();
        let state = model.evaluate_pt(Temperature::new(t), Pressure::new(1.0e-10)).unwrap();
        let root = state.select(PengRobinsonPhaseV1::Vapor);
        let v = root.molar_volume_m3_per_mol();
        for i in 0..2 {
            let partial_a = (0..2).map(|j| model.mole_fractions[j] * mixed.aij[2 * i + j])
                .sum::<f64>();
            let virial = (mixed.bi[i] + (mixed.a_mix - 2.0 * partial_a) / rt) / v;
            assert!((root.log_fugacity_coefficients()[i] / virial - 1.0).abs() < 1.0e-12);
        }
    }

    #[test]
    fn g1_dilute_caloric_and_entropy_tangents_preserve_cv() {
        let model = pure(methane());
        let t = 300.0;
        let density = Density::new(1.0e-15);
        let evaluate = |temperature| model.evaluate_density_temperature(
            Temperature::new(temperature), density,
        ).unwrap();
        let center = evaluate(t);
        let dt = 1.0e-3;
        let plus = evaluate(t + dt);
        let minus = evaluate(t - dt);
        let du_dt = (plus.root().residual_internal_energy_j_mol()
            - minus.root().residual_internal_energy_j_mol()) / (2.0 * dt);
        let ds_dt = (plus.root().residual_entropy_at_density_j_mol_k()
            - minus.root().residual_entropy_at_density_j_mol_k()) / (2.0 * dt);
        let cv = center.root().residual_cv_j_mol_k();
        assert!(cv > 0.0);
        assert!((du_dt / cv - 1.0).abs() < 1.0e-7);
        assert!((t * ds_dt / cv - 1.0).abs() < 1.0e-7);
    }

    #[test]
    fn g0_density_input_refuses_covolume_spinodal_and_invalid_values() {
        let model = pure(methane());
        for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(model.evaluate_density_temperature(Temperature::new(300.0), Density::new(density)),
                Err(PengRobinsonErrorV1::InvalidDensity { .. })));
        }
        assert!(matches!(model.evaluate_density_temperature(Temperature::new(300.0), Density::new(1000.0)),
            Err(PengRobinsonErrorV1::OutsideCovolume)));
        let subcritical = model.evaluate_pt(Temperature::new(150.0), Pressure::new(1.0e6)).unwrap();
        let middle = &subcritical.phase_roots()[1];
        assert!(middle.pressure_derivatives().1 > 0.0);
        assert!(matches!(model.evaluate_density_temperature(Temperature::new(150.0),
            Density::new(middle.mass_density_kg_per_m3())), Err(PengRobinsonErrorV1::MechanicallyUnstable)));
    }

    #[test]
    fn malformed_source_and_state_data_refuse_without_repair() {
        let c = methane();
        assert!(matches!(
            PengRobinsonMixtureV1::new(
                vec![c, c],
                vec![0.5, 0.5],
                vec![0.0, 0.1, 0.2, 0.0]
            ),
            Err(PengRobinsonErrorV1::AsymmetricBinaryInteraction { .. })
        ));
        assert!(matches!(
            PengRobinsonMixtureV1::new(vec![c], vec![0.9], vec![]),
            Err(PengRobinsonErrorV1::CompositionNotNormalized { .. })
        ));
        let model = pure(c);
        assert!(matches!(
            model.evaluate_pt(Temperature::new(0.0), Pressure::new(1.0)),
            Err(PengRobinsonErrorV1::InvalidTemperature { .. })
        ));
        assert!(matches!(
            model.evaluate_pt(Temperature::new(300.0), Pressure::new(f64::NAN)),
            Err(PengRobinsonErrorV1::InvalidPressure { .. })
        ));
    }
}
