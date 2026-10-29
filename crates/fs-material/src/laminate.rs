//! Laminate (layered) plate/shell sections (IRP M1.2-W7).
//!
//! Classical lamination theory on a stack of oriented layers, each carrying its
//! own resolved source material and either an elastic or a viscoelastic
//! (tensor-Prony) plane-stress law. The section integrates, about one declared
//! reference surface, the membrane/coupling/bending stiffness
//! `A`, `B`, `D` (3×3 each, engineering-shear `[ε_x, ε_y, γ_xy]`), the
//! energy-consistent transverse shear stiffness `H` (2×2, Whitney 1973), and the
//! translational/coupling/rotary inertia `I0`, `I1`, `I2`.
//!
//! Each layer is integrated in declaration order between its bottom and top
//! `z` faces, measured from the section's physical bottom; the reference surface
//! then shifts every first moment. A thin 20–50 µm finish on 3 mm wood enters
//! `A`, `B` and `D` exactly: no layer is neglected, and the per-direction output
//! shows its weak-direction contribution.
//!
//! Nothing here is averaged into one effective material. The per-layer rotated
//! stiffness and the source bundles are retained so a consumer (M1.2-W5 loss,
//! the M1.2-W1a MITC element) can recover every layer. Viscoelastic layers add
//! `complex_abd(ω)`, whose imaginary parts are a passive loss operator by
//! construction (each layer's `plane_stress_complex` verifies `Im Q ⪰ 0`).
//!
//! Determinism: layers integrate in declaration order; the only transcendentals
//! are `fs_math::det` rotations. fs-material is D0.

use crate::elastic::OrthotropicElastic;
use crate::plane_stress_viscosity::MaterialPlane;
use crate::state_point::{OrthotropicPlaneElastic, ResolvedMaterialStatePoint};
use crate::visco::ViscoError;
use crate::visco::tensor::TensorProny;
use fs_math::c64::C64;

/// Largest laminate the section admits, matching the rod precedent
/// [`crate::section::MAX_CONCENTRIC_LAYERS`].
pub const MAX_LAMINA_LAYERS: usize = 128;

/// Typed refusals of a laminate section. No partial section is ever returned.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminateError {
    /// The layer count is zero or exceeds [`MAX_LAMINA_LAYERS`].
    LayerCount {
        /// Number of layers supplied.
        count: usize,
    },
    /// A layer thickness is not finite and strictly positive.
    Thickness {
        /// Zero-based layer index.
        layer: usize,
    },
    /// A material constant is not finite/positive, or the rotation or an
    /// aggregate coefficient is unrepresentable.
    Invalid {
        /// What was being formed.
        quantity: &'static str,
    },
    /// The reference-surface offset is not finite.
    Reference,
    /// A viscoelastic layer refused at the queried frequency (out of band or
    /// non-passive reduced operator).
    Layer {
        /// Zero-based layer index.
        layer: usize,
        /// The underlying viscoelastic refusal.
        source: ViscoError,
    },
}

impl core::fmt::Display for LaminateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LaminateError::LayerCount { count } => {
                write!(
                    f,
                    "FS-MATERIAL-LAMINATE-LAYERS: {count} layers, need 1..={MAX_LAMINA_LAYERS}"
                )
            }
            LaminateError::Thickness { layer } => {
                write!(
                    f,
                    "FS-MATERIAL-LAMINATE-THICKNESS: layer {layer} thickness is not finite and positive"
                )
            }
            LaminateError::Invalid { quantity } => {
                write!(f, "FS-MATERIAL-LAMINATE-INVALID: {quantity}")
            }
            LaminateError::Reference => {
                write!(
                    f,
                    "FS-MATERIAL-LAMINATE-REFERENCE: reference-surface offset must be finite"
                )
            }
            LaminateError::Layer { layer, source } => {
                write!(f, "FS-MATERIAL-LAMINATE-LAYER: layer {layer}: {source}")
            }
        }
    }
}
impl std::error::Error for LaminateError {}

/// The plane-stress constitutive law of one layer, in the 1–2 material plane
/// (axis 1 at `orientation_rad` counterclockwise from mesh x).
#[derive(Clone, Debug, PartialEq)]
pub enum LayerLaw {
    /// Fixed-state orthotropic elasticity. The in-plane reduced stiffness comes
    /// from `E1, E2, ν12, G12`; the transverse shear from `G13, G23`.
    Elastic(OrthotropicElastic),
    /// A passive tensor-Prony law (M0.2-W2a). The in-plane reduced stiffness is
    /// its frequency-dependent plane-stress Schur complement; the transverse
    /// shear uses the equilibrium (`ω = 0`) out-of-plane moduli.
    Prony(TensorProny),
    /// Independently characterized reduced plane-stress stiffness and transverse
    /// shear moduli. This does not invent an unmeasured third Young modulus or
    /// out-of-plane Poisson ratios for an existing shell material source.
    PlaneStress {
        /// Symmetric positive-definite engineering-shear matrix [Pa].
        stiffness_pa: [f64; 9],
        /// Principal transverse shear moduli G13, G23 [Pa].
        transverse_shear_pa: [f64; 2],
    },
}

/// One lamina: a thickness, an in-plane orientation, its resolved source
/// material (retained for provenance, never averaged) and its law.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminaLayer {
    /// Layer thickness [m], finite and strictly positive.
    pub thickness_m: f64,
    /// Material axis 1, counterclockwise from mesh x [rad].
    pub orientation_rad: f64,
    /// Density [kg/m³], finite and strictly positive.
    pub density_kg_m3: f64,
    /// Exact source bundle for this layer's constants, retained for provenance
    /// when the caller resolved one. The section math reads only `law`,
    /// `thickness_m`, `orientation_rad` and `density_kg_m3`.
    pub material: Option<ResolvedMaterialStatePoint>,
    /// Plane-stress law.
    pub law: LayerLaw,
}

/// Where the first moments `B`, `I1` are measured from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReferenceSurface {
    /// The geometric midsurface at half the total thickness.
    Midsurface,
    /// A plane `offset` metres above the physical bottom face.
    Offset(f64),
}

/// Per-layer integrals and retained geometry, in declaration (spatial) order.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminaIntegrals {
    /// Bottom and top `z` faces relative to the reference surface [m].
    pub z_bottom_m: f64,
    /// Top face relative to the reference surface [m].
    pub z_top_m: f64,
    /// Rotated elastic (or ω = 0 storage) reduced stiffness `Q̄` [Pa], row-major.
    pub q_bar_pa: [f64; 9],
    /// Rotated transverse shear moduli `[G_xz, G_yz]` [Pa].
    pub transverse_shear_pa: [f64; 2],
}

/// An immutable laminate section with all source bundles, retained layer
/// stiffness and the aggregate CLT coefficients about the reference surface.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateSection {
    layers: Vec<LaminaLayer>,
    integrals: Vec<LaminaIntegrals>,
    reference: ReferenceSurface,
    total_thickness_m: f64,
    a_pa_m: [f64; 9],
    b_pa_m2: [f64; 9],
    d_pa_m3: [f64; 9],
    shear_pa_m: [f64; 4],
    shear_factor: [f64; 2],
    i0_kg_m2: f64,
    i1_kg_m: f64,
    i2_kg: f64,
}

fn finite_positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// `T(θ)ᵀ Q T(θ)`: rotate a reduced plane-stress stiffness (engineering-shear,
/// row-major) by `angle` about the plate normal, exactly as
/// [`crate::sheet::SheetMatrix::rotated`].
fn rotate_reduced(q: &[f64; 9], angle: f64) -> Option<[f64; 9]> {
    if angle == 0.0 {
        return Some(*q);
    }
    let (s, c) = angle.sin_cos();
    let (cc, ss, cs) = (c * c, s * s, c * s);
    let t = [[cc, ss, cs], [ss, cc, -cs], [-2.0 * cs, 2.0 * cs, cc - ss]];
    let mut out = [0.0f64; 9];
    for i in 0..3 {
        for j in i..3 {
            let mut v = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    v += t[a][i] * q[3 * a + b] * t[b][j];
                }
            }
            out[3 * i + j] = v;
            out[3 * j + i] = v;
        }
    }
    out.iter().all(|v| v.is_finite()).then_some(out)
}

/// Rotate the in-plane `[G_xz, G_yz]` diagonal transverse shear by `angle`
/// about the normal: `diag(G_xz, G_yz)` transforms as a symmetric 2-tensor.
fn rotate_transverse(g: [f64; 2], angle: f64) -> [f64; 4] {
    if angle == 0.0 {
        return [g[0], 0.0, 0.0, g[1]];
    }
    let (s, c) = angle.sin_cos();
    let (cc, ss) = (c * c, s * s);
    let cross = (g[0] - g[1]) * c * s;
    [g[0] * cc + g[1] * ss, cross, cross, g[0] * ss + g[1] * cc]
}

/// The reference (band-low storage) reduced in-plane stiffness and the
/// transverse shear moduli of a layer, before rotation. For an elastic layer
/// this is exact; for a Prony layer it is the storage at the band start, the
/// frequency-independent reference for the stored real `A`, `B`, `D` (the
/// frequency dependence lives in [`LaminateSection::complex_abd`]).
fn principal_stiffness(
    law: &LayerLaw,
    layer: usize,
) -> Result<([f64; 9], [f64; 2]), LaminateError> {
    match law {
        LayerLaw::PlaneStress { stiffness_pa, transverse_shear_pa } => {
            crate::sheet::SheetMatrix::from_coefficients(*stiffness_pa).map_err(|_| {
                LaminateError::Invalid { quantity: "reduced layer stiffness must be symmetric positive definite" }
            })?;
            if !transverse_shear_pa.iter().all(|v| finite_positive(*v)) {
                return Err(LaminateError::Invalid { quantity: "reduced layer transverse shear modulus" });
            }
            Ok((*stiffness_pa, *transverse_shear_pa))
        }
        LayerLaw::Elastic(e) => {
            let plane = OrthotropicPlaneElastic {
                plane: MaterialPlane::Axes12,
                e: [e.e[0], e.e[1]],
                nu: e.nu[0],
                g: e.g[0],
                strain_limit: e.strain_limit,
            };
            let q = plane.stiffness();
            if !q.iter().all(|v| v.is_finite()) || q[0] <= 0.0 || q[4] <= 0.0 || q[8] <= 0.0 {
                return Err(LaminateError::Invalid {
                    quantity: "elastic reduced stiffness",
                });
            }
            // G13 (x–z) and G23 (y–z) with the 3-axis normal to the plate.
            // OrthotropicElastic stores [G12, G23, G31].
            let shear = [e.g[2], e.g[1]];
            if !shear.iter().all(|v| finite_positive(*v)) {
                return Err(LaminateError::Invalid {
                    quantity: "transverse shear modulus",
                });
            }
            Ok((q, shear))
        }
        LayerLaw::Prony(p) => {
            let reference = p.band_rad_s.0;
            let (q_re, _q_im) = p
                .plane_stress_complex(MaterialPlane::Axes12, reference)
                .map_err(|source| LaminateError::Layer { layer, source })?;
            // The transverse shear uses the out-of-plane moduli at the same
            // reference, so the correction is a static property of the stack.
            let (storage0, _) = p
                .complex_modulus(reference)
                .map_err(|source| LaminateError::Layer { layer, source })?;
            // Crate Mandel order is [11,22,33,12,13,23]; a shear diagonal holds
            // 2·G. The 3-axis is the plate normal, so G_xz = G13 (index 4) and
            // G_yz = G23 (index 5).
            let shear = [0.5 * storage0[4][4], 0.5 * storage0[5][5]];
            if !shear.iter().all(|v| finite_positive(*v)) {
                return Err(LaminateError::Invalid {
                    quantity: "transverse shear modulus",
                });
            }
            Ok((q_re, shear))
        }
    }
}

impl LaminateSection {
    /// Integrate the classical-lamination-theory coefficients of `layers` about
    /// `reference`. Layer order is spatial and MUST NOT be sorted by material.
    ///
    /// # Errors
    /// [`LaminateError`] for an empty/oversized stack, a nonpositive or
    /// non-finite thickness, a non-finite reference offset, a nonpositive or
    /// unrepresentable material constant, or a viscoelastic layer that refuses
    /// at its own band start.
    pub fn try_new(
        layers: Vec<LaminaLayer>,
        reference: ReferenceSurface,
    ) -> Result<Self, LaminateError> {
        if layers.is_empty() || layers.len() > MAX_LAMINA_LAYERS {
            return Err(LaminateError::LayerCount {
                count: layers.len(),
            });
        }
        let offset = match reference {
            ReferenceSurface::Midsurface => None,
            ReferenceSurface::Offset(o) => {
                if !o.is_finite() {
                    return Err(LaminateError::Reference);
                }
                Some(o)
            }
        };
        let mut total = 0.0f64;
        for (i, layer) in layers.iter().enumerate() {
            if !finite_positive(layer.thickness_m) {
                return Err(LaminateError::Thickness { layer: i });
            }
            if !finite_positive(layer.density_kg_m3) || !layer.orientation_rad.is_finite() {
                return Err(LaminateError::Invalid {
                    quantity: "layer density or orientation",
                });
            }
            total += layer.thickness_m;
        }
        if !finite_positive(total) {
            return Err(LaminateError::Invalid {
                quantity: "total thickness",
            });
        }
        let reference_z = match offset {
            None => 0.5 * total,
            Some(o) => o,
        };
        // Integrate from the physical bottom (z = 0), then subtract the
        // reference height so every first and second moment is about it.
        let mut a = [0.0f64; 9];
        let mut b = [0.0f64; 9];
        let mut d = [0.0f64; 9];
        let mut shear = [0.0f64; 4];
        let (mut i0, mut i1, mut i2) = (0.0f64, 0.0f64, 0.0f64);
        let mut integrals = Vec::with_capacity(layers.len());
        let mut bottom = 0.0f64;
        for (i, layer) in layers.iter().enumerate() {
            let top = bottom + layer.thickness_m;
            let z0 = bottom - reference_z;
            let z1 = top - reference_z;
            let (q_principal, shear_principal) = principal_stiffness(&layer.law, i)?;
            let q = rotate_reduced(&q_principal, layer.orientation_rad).ok_or(
                LaminateError::Invalid {
                    quantity: "rotated reduced stiffness",
                },
            )?;
            let gs = rotate_transverse(shear_principal, layer.orientation_rad);
            // Factored differences keep a thin layer representable: a thin
            // finish's (z1 − z0) stays exact where z1²−z0² would cancel.
            let m0 = z1 - z0;
            let m1 = 0.5 * (z1 - z0) * (z1 + z0);
            let m2 = (z1 * z1 * z1 - z0 * z0 * z0) / 3.0;
            for k in 0..9 {
                a[k] += q[k] * m0;
                b[k] += q[k] * m1;
                d[k] += q[k] * m2;
            }
            // Uncorrected transverse shear stiffness ∫G dz; the energy-consistent
            // Whitney factor multiplies it below.
            for k in 0..4 { shear[k] += gs[k] * m0; }
            let rho = layer.density_kg_m3;
            i0 += rho * m0;
            i1 += rho * m1;
            i2 += rho * m2;
            integrals.push(LaminaIntegrals {
                z_bottom_m: z0,
                z_top_m: z1,
                q_bar_pa: q,
                transverse_shear_pa: [gs[0], gs[3]],
            });
            bottom = top;
        }
        // Energy-consistent Whitney (1973) correction, per in-plane direction,
        // about each direction's own bending-neutral axis (κ = 5/6 for one
        // homogeneous layer). `shear[0]`/`shear[3]` currently hold ∫G dz.
        let kappa_x = transverse_shear_factor(&integrals, 0, 0, shear[0]);
        let kappa_y = transverse_shear_factor(&integrals, 4, 1, shear[3]);
        shear[0] *= kappa_x;
        shear[3] *= kappa_y;
        // A positive diagonal congruence retains the complete transverse
        // tensor and passivity. For off-axis multilayers these directional
        // Whitney factors are an approximation, not a coupled warping solve.
        let cross_correction = kappa_x.sqrt() * kappa_y.sqrt();
        shear[1] *= cross_correction;
        shear[2] = shear[1];
        let finite_9 = |m: &[f64; 9]| m.iter().all(|v| v.is_finite());
        if !(finite_9(&a) && finite_9(&b) && finite_9(&d))
            || !shear.iter().all(|v| v.is_finite())
            || !(i0.is_finite() && i1.is_finite() && i2.is_finite())
            || !(finite_positive(i0) && i2 > 0.0)
            || a[0] <= 0.0
            || a[4] <= 0.0
            || a[8] <= 0.0
        {
            return Err(LaminateError::Invalid {
                quantity: "aggregate section coefficient",
            });
        }
        Ok(Self {
            layers,
            integrals,
            reference,
            total_thickness_m: total,
            a_pa_m: a,
            b_pa_m2: b,
            d_pa_m3: d,
            shear_pa_m: shear,
            shear_factor: [kappa_x, kappa_y],
            i0_kg_m2: i0,
            i1_kg_m: i1,
            i2_kg: i2,
        })
    }

    /// Ordered source layers; no averaged material receipt.
    #[must_use]
    pub fn layers(&self) -> &[LaminaLayer] {
        &self.layers
    }
    /// Per-layer integrals in the same spatial order as [`Self::layers`].
    #[must_use]
    pub fn layer_integrals(&self) -> &[LaminaIntegrals] {
        &self.integrals
    }
    /// The declared reference surface.
    #[must_use]
    pub const fn reference(&self) -> ReferenceSurface {
        self.reference
    }
    /// Total stack thickness [m].
    #[must_use]
    pub const fn total_thickness_m(&self) -> f64 {
        self.total_thickness_m
    }
    /// Membrane stiffness `A` [Pa·m], row-major 3×3.
    #[must_use]
    pub const fn a_pa_m(&self) -> &[f64; 9] {
        &self.a_pa_m
    }
    /// Membrane–bending coupling `B` [Pa·m²], row-major 3×3.
    #[must_use]
    pub const fn b_pa_m2(&self) -> &[f64; 9] {
        &self.b_pa_m2
    }
    /// Bending stiffness `D` [Pa·m³], row-major 3×3.
    #[must_use]
    pub const fn d_pa_m3(&self) -> &[f64; 9] {
        &self.d_pa_m3
    }
    /// Full symmetric transverse shear stiffness `H` [Pa·m], row-major 2×2
    /// in `[xz,yz]`, including off-axis cross coupling.
    #[must_use]
    pub const fn shear_pa_m(&self) -> &[f64; 4] {
        &self.shear_pa_m
    }
    /// Directional Whitney shear correction factors `[κ_x, κ_y]`
    /// (`5/6` for one homogeneous layer). Cross shear uses their geometric
    /// mean; off-axis multilayers do not claim a coupled 3-D warping solution.
    #[must_use]
    pub const fn shear_factor(&self) -> &[f64; 2] {
        &self.shear_factor
    }
    /// Translational inertia `I0 = ∫ρ dz` [kg/m²].
    #[must_use]
    pub const fn i0_kg_m2(&self) -> f64 {
        self.i0_kg_m2
    }
    /// Inertia coupling `I1 = ∫ρz dz` [kg/m]; zero for a symmetric stack.
    #[must_use]
    pub const fn i1_kg_m(&self) -> f64 {
        self.i1_kg_m
    }
    /// Rotary inertia `I2 = ∫ρz² dz` [kg].
    #[must_use]
    pub const fn i2_kg(&self) -> f64 {
        self.i2_kg
    }

    /// Complex `(A, B, D)` [Pa·m, Pa·m², Pa·m³] at `omega` [rad/s]. Each entry is
    /// `C64(storage, loss)`; an elastic layer contributes a real part only, a
    /// Prony layer its frequency-dependent plane-stress Schur complement. The
    /// imaginary parts are a passive (`Im ⪰ 0`) loss operator by construction.
    ///
    /// # Errors
    /// [`LaminateError::Layer`] if any Prony layer is out of band or non-passive
    /// at `omega`.
    pub fn complex_abd(&self, omega: f64) -> Result<([C64; 9], [C64; 9], [C64; 9]), LaminateError> {
        let mut a = [C64::ZERO; 9];
        let mut b = [C64::ZERO; 9];
        let mut d = [C64::ZERO; 9];
        for (i, (layer, integ)) in self.layers.iter().zip(&self.integrals).enumerate() {
            let (q_re, q_im) = match &layer.law {
                LayerLaw::Elastic(_) | LayerLaw::PlaneStress { .. } => (integ.q_bar_pa, [0.0f64; 9]),
                LayerLaw::Prony(p) => {
                    let (re, im) = p
                        .plane_stress_complex(MaterialPlane::Axes12, omega)
                        .map_err(|source| LaminateError::Layer { layer: i, source })?;
                    let re = rotate_reduced(&re, layer.orientation_rad).ok_or(
                        LaminateError::Invalid {
                            quantity: "rotated complex storage",
                        },
                    )?;
                    let im = rotate_reduced(&im, layer.orientation_rad).ok_or(
                        LaminateError::Invalid {
                            quantity: "rotated complex loss",
                        },
                    )?;
                    (re, im)
                }
            };
            let (z0, z1) = (integ.z_bottom_m, integ.z_top_m);
            let m0 = z1 - z0;
            let m1 = 0.5 * (z1 - z0) * (z1 + z0);
            let m2 = (z1 * z1 * z1 - z0 * z0 * z0) / 3.0;
            for k in 0..9 {
                let q = C64::new(q_re[k], q_im[k]);
                a[k] = a[k] + q.scale(m0);
                b[k] = b[k] + q.scale(m1);
                d[k] = d[k] + q.scale(m2);
            }
        }
        Ok((a, b, d))
    }

    /// Integrate one principal-frame resultant matrix per retained layer using
    /// this section's exact geometry and reference surface. For reduced
    /// viscosities [Pa s], the outputs are membrane, coupling and bending
    /// viscosities [N s/m, N s, N m s]. This is the same CLT integration as the
    /// elastic section; no homogeneous h³/12 replacement discards a coating's
    /// offset. The caller admits the constitutive matrices' passivity.
    pub fn integrate_layer_matrices(&self, matrices: &[[f64; 9]])
        -> Result<([f64; 9], [f64; 9], [f64; 9]), LaminateError>
    {
        if matrices.len() != self.layers.len() {
            return Err(LaminateError::Invalid { quantity: "one resultant matrix is required per layer" });
        }
        let (mut a, mut b, mut d) = ([0.0; 9], [0.0; 9], [0.0; 9]);
        for ((matrix, layer), geometry) in matrices.iter().zip(&self.layers).zip(&self.integrals) {
            let q = rotate_reduced(matrix, layer.orientation_rad).ok_or(
                LaminateError::Invalid { quantity: "nonfinite layer resultant" })?;
            let (z0, z1) = (geometry.z_bottom_m, geometry.z_top_m);
            let m0 = z1 - z0;
            let m1 = 0.5 * m0 * (z1 + z0);
            let m2 = m0 * (z0*z0 + z0*z1 + z1*z1) / 3.0;
            for k in 0..9 { a[k] += q[k]*m0; b[k] += q[k]*m1; d[k] += q[k]*m2; }
        }
        if a.iter().chain(&b).chain(&d).any(|v| !v.is_finite()) {
            return Err(LaminateError::Invalid { quantity: "unrepresentable layer resultants" });
        }
        Ok((a,b,d))
    }
}

/// Whitney (1973) energy-consistent transverse shear correction factor in one
/// in-plane direction: it matches the shear strain energy of the parabolic
/// equilibrium stress to the constant first-order-theory one.
///
/// `in_plane` indexes the in-plane modulus in `q_bar` (`0` → Q̄₁₁ for x, `4` →
/// Q̄₂₂ for y); `shear_idx` the transverse modulus (`0` → G_xz, `1` → G_yz);
/// `shear_stiffness` is the uncorrected `∫G dz`.
///
/// Under a transverse shear force the equilibrium stress is
/// `τ(z) = −(V/D_n)·R(z)`, where `R(z) = ∫ Q̄(ζ)(ζ − z_n) dζ` is the first
/// moment of the in-plane modulus about the bending-neutral axis `z_n`
/// (`R` vanishes at both faces). Matching `½V²/(κ·∫G dz)` to `½(V²/D_n²)∫R²/G dz`
/// gives `κ = D_n² / (∫G dz · ∫R²/G dz)`. Each integral is a layerwise exact
/// polynomial. A single homogeneous orthotropic layer yields exactly `5/6`.
fn transverse_shear_factor(
    integrals: &[LaminaIntegrals],
    in_plane: usize,
    shear_idx: usize,
    shear_stiffness: f64,
) -> f64 {
    let fallback = 5.0 / 6.0;
    if !(shear_stiffness.is_finite() && shear_stiffness > 0.0) {
        return fallback;
    }
    // Bending-neutral axis z_n = (∫Q̄ z dz)/(∫Q̄ dz), about the reference.
    let mut a = 0.0f64; // ∫Q̄ dz
    let mut first = 0.0f64; // ∫Q̄ z dz
    for integ in integrals {
        let q = integ.q_bar_pa[in_plane];
        if !(q.is_finite() && q > 0.0) {
            return fallback;
        }
        let (z0, z1) = (integ.z_bottom_m, integ.z_top_m);
        a += q * (z1 - z0);
        first += q * 0.5 * (z1 - z0) * (z1 + z0);
    }
    if !(a.is_finite() && a > 0.0) {
        return fallback;
    }
    let z_n = first / a;
    // D_n = ∫Q̄ (z − z_n)² dz, and energy = ∫ R(z)²/G dz with
    // R(z) = ∫_bottom^z Q̄ (ζ − z_n) dζ (a per-layer quadratic in z).
    let mut d_n = 0.0f64;
    let mut energy = 0.0f64;
    let mut running = 0.0f64; // R at the current layer bottom
    for integ in integrals {
        let q = integ.q_bar_pa[in_plane];
        let g = integ.transverse_shear_pa[shear_idx];
        if !(g.is_finite() && g > 0.0) {
            return fallback;
        }
        let (z0, z1) = (integ.z_bottom_m, integ.z_top_m);
        // Centred moments about z_n (factored to stay representable for thin layers).
        let (u0, u1) = (z0 - z_n, z1 - z_n);
        d_n += q * (u1 * u1 * u1 - u0 * u0 * u0) / 3.0;
        // R(z) = running + q·((z−z_n)² − u0²)/2 = c0 + (q/2)(z−z_n)², with
        // c0 = running − (q/2)u0². Integrate R²/g over [z0, z1] exactly in u.
        let half_q = 0.5 * q;
        let c0 = running - half_q * u0 * u0;
        let int_1 = u1 - u0;
        let int_u2 = (u1.powi(3) - u0.powi(3)) / 3.0;
        let int_u4 = (u1.powi(5) - u0.powi(5)) / 5.0;
        let r2 = c0 * c0 * int_1 + 2.0 * c0 * half_q * int_u2 + half_q * half_q * int_u4;
        energy += r2 / g;
        running += half_q * (u1 * u1 - u0 * u0);
    }
    if !(d_n.is_finite() && d_n > 0.0 && energy.is_finite() && energy > 0.0) {
        return fallback;
    }
    let kappa = d_n * d_n / (shear_stiffness * energy);
    if kappa.is_finite() && kappa > 0.0 {
        kappa
    } else {
        fallback
    }
}
