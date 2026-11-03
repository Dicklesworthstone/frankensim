//! Coupled energy equivalence for declared cylindrical-bending probes.
//!
//! CLT gives sigma(z) = Q(z) (z I - A^-1 B) kappa when the local membrane
//! resultants vanish. For bending in direction n, kappa' = [nx²,ny²,2nx ny]
//! and transverse equilibrium gives tau(z) = -integral(div sigma) dz.
//! Two orthogonal physical probes span a 2x2 resultant matrix V. Their
//! normalized traction T(z) = tau(z) V^-1 integrates to I, so
//! H^-1 = integral(T(z)^T G(z)^-1 T(z)) dz preserves full shear work.
//!
//! This extends the cylindrical energy-equivalence construction in Hadavinia
//! et al., doi:10.3970/sdhm.2006.002.197, equations 9 and 12, to the declared
//! two-probe subspace. It is not a unique, load-independent 3-D warping law:
//! the choice of shear gradients matters for angle-ply laminates (Loredo,
//! <https://arxiv.org/abs/1211.0781>, sections 1 and 3).

use super::{
    LaminaLayer, LaminateError, LaminateSection, LaminateShearModel, LayerLaw, ReferenceSurface,
    principal_stiffness, rotate_transverse,
};
use fs_ad::Real;

fn invalid(quantity: &'static str) -> LaminateError {
    LaminateError::Invalid { quantity }
}

impl LaminateSection {
    /// Build the original elastic A/B/D and inertia with full energy-equivalent
    /// transverse shear for two declared orthogonal cylindrical-bending probes.
    /// `bending_axis_rad` is a physical direction, not an inferred grain axis.
    /// It must co-rotate with the laminate when changing mesh coordinates.
    ///
    /// The probes use local zero membrane resultants, retaining the complete
    /// `A^-1 B` coupling and every layer's rotated transverse shear compliance.
    /// Their continuous traction vanishes at both faces, integrates to the
    /// prescribed shear resultant, and supplies its full cross-energy. A single
    /// homogeneous layer gives `H = (5/6) h G` at every admitted probe angle.
    /// Preparation is O(layers); existing shell consumers use the cached 2x2 H.
    ///
    /// This is a frozen linear-elastic, two-probe reduction. It does not claim
    /// arbitrary loading/warping equivalence, interlaminar normal stress, free
    /// edge accuracy or viscoelastic memory. The thermal owner differentiates
    /// this same kernel using its original source Q/G derivatives at fixed
    /// geometry and physical probe axes. [`Self::try_new`] retains the distinct
    /// directional Whitney model.
    /// The cylindrical energy-equivalence basis is described by
    /// [Hadavinia et al.](https://doi.org/10.3970/sdhm.2006.002.197); the need
    /// to specify shear gradients for angle-ply reductions is discussed by
    /// [Loredo](https://arxiv.org/abs/1211.0781).
    ///
    /// # Errors
    /// As for [`Self::try_new`], plus nonfinite probe angles, Prony layers,
    /// unresolved probe rank, nonpositive shear energy or unrepresentable work.
    pub fn try_new_equilibrated_shear(
        layers: Vec<LaminaLayer>,
        reference: ReferenceSurface,
        bending_axis_rad: f64,
    ) -> Result<Self, LaminateError> {
        if !bending_axis_rad.is_finite() {
            return Err(invalid(
                "equilibrated shear needs a finite physical bending axis",
            ));
        }
        if layers
            .iter()
            .any(|layer| matches!(layer.law, LayerLaw::Prony(_)))
        {
            return Err(invalid(
                "equilibrated shear currently requires elastic layers",
            ));
        }
        let mut section = Self::try_new(layers, reference)?;
        let (stiffness, factors) = prepare(&section, bending_axis_rad)?;
        section.shear_pa_m = stiffness;
        section.shear_factor = factors;
        section.shear_model = LaminateShearModel::EquilibratedCylindrical { bending_axis_rad };
        Ok(section)
    }
}

/// One original ply in section axes. Thickness stays independent of any offset
/// reference surface; Q and G may carry the original thermal derivatives.
pub(super) struct EquilibratedLayer<T> {
    pub thickness_m: f64,
    pub q: [T; 9],
    pub g: [T; 4],
}

fn prepare(section: &LaminateSection, angle: f64) -> Result<([f64; 4], [f64; 2]), LaminateError> {
    let mut layers = Vec::with_capacity(section.layers.len());
    for (i, (layer, geometry)) in section.layers.iter().zip(&section.integrals).enumerate() {
        let (_, principal) = principal_stiffness(&layer.law, i)?;
        layers.push(EquilibratedLayer {
            thickness_m: layer.thickness_m,
            q: geometry.q_bar_pa,
            g: rotate_transverse(principal, layer.orientation_rad),
        });
    }
    equilibrated(&layers, angle, invalid, &mut || Ok(()))
}

/// Columns are the shear resultants/tractions of the two physical probes.
fn divergence<T: Real>(matrix: &[T; 9], axes: [[f64; 2]; 2]) -> [T; 4] {
    let mut out = [T::zero(); 4];
    for (column, n) in axes.iter().enumerate() {
        let curvature = [n[0] * n[0], n[1] * n[1], 2.0 * n[0] * n[1]];
        let stress: [T; 3] = core::array::from_fn(|i| {
            (0..3).fold(T::zero(), |sum, j| {
                sum + matrix[3 * i + j] * T::from_f64(curvature[j])
            })
        });
        out[column] = stress[0] * T::from_f64(n[0]) + stress[2] * T::from_f64(n[1]);
        out[2 + column] = stress[2] * T::from_f64(n[0]) + stress[1] * T::from_f64(n[1]);
    }
    out
}

/// The one coupled shear expression used by mechanical preparation and AD.
/// All numerical scales are derivative-free primal values: each cancels
/// algebraically, so freezing it neither loses source derivatives nor
/// differentiates the irrelevant max/absolute-value choices. Admission checks
/// remain on the primal branch; no derivative is promised across a refusal.
pub(super) fn equilibrated<T: Real, E>(
    layers: &[EquilibratedLayer<T>],
    angle: f64,
    invalid: impl Fn(&'static str) -> E,
    check: &mut impl FnMut() -> Result<(), E>,
) -> Result<([T; 4], [T; 2]), E> {
    check()?;
    // Dimensionless thickness and moduli avoid powers of SI scales and any
    // cancellation from the caller's choice of an offset reference surface.
    let h: f64 = layers.iter().map(|layer| layer.thickness_m).sum();
    let q_scale = layers
        .iter()
        .flat_map(|layer| layer.q)
        .map(|x| x.value().abs())
        .fold(0.0_f64, f64::max);
    let g_scale = layers
        .iter()
        .flat_map(|layer| layer.g)
        .map(|x| x.value().abs())
        .fold(0.0_f64, f64::max);
    if !angle.is_finite()
        || !h.is_finite()
        || h <= 0.0
        || !q_scale.is_finite()
        || q_scale <= 0.0
        || !g_scale.is_finite()
        || g_scale <= 0.0
    {
        return Err(invalid(
            "equilibrated shear scales and physical axes are invalid",
        ));
    }
    let (mut a, mut b, mut d) = ([T::zero(); 9], [T::zero(); 9], [T::zero(); 9]);
    let mut z0 = -0.5;
    for layer in layers {
        check()?;
        let q = layer.q.map(|v| v / T::from_f64(q_scale));
        let dz = layer.thickness_m / h;
        let z1 = z0 + dz;
        let moments = [
            dz,
            0.5 * dz * (z0 + z1),
            dz * (z0 * z0 + z0 * z1 + z1 * z1) / 3.0,
        ];
        for k in 0..9 {
            a[k] = a[k] + q[k] * T::from_f64(moments[0]);
            b[k] = b[k] + q[k] * T::from_f64(moments[1]);
            d[k] = d[k] + q[k] * T::from_f64(moments[2]);
        }
        z0 = z1;
    }
    let neutral = neutral_matrix(a, b, &invalid)?;
    let relaxed: [T; 9] = core::array::from_fn(|ij| {
        d[ij]
            - (0..3).fold(T::zero(), |sum, k| {
                sum + b[3 * (ij / 3) + k] * neutral[3 * k + ij % 3]
            })
    });
    let (sin, cos) = angle.sin_cos();
    let axes = [[cos, sin], [-sin, cos]];
    let v = divergence(&relaxed, axes);
    let v_scale = v.iter().map(|x| x.value().abs()).fold(0.0_f64, f64::max);
    let scaled = v.map(|x| x / T::from_f64(v_scale));
    let determinant = scaled[0] * scaled[3] - scaled[1] * scaled[2];
    if !determinant.value().is_finite() || determinant.value().abs() <= 128.0 * f64::EPSILON {
        return Err(invalid(
            "physical bending probes do not resolve two independent shear resultants",
        ));
    }
    let inverse = [scaled[3], -scaled[1], -scaled[2], scaled[0]]
        .map(|x| x / determinant / T::from_f64(v_scale));

    let (mut primitive, mut compliance, mut raw) = ([T::zero(); 9], [T::zero(); 4], [T::zero(); 4]);
    z0 = -0.5;
    let gauss = (3.0_f64 / 5.0).sqrt();
    for layer in layers {
        check()?;
        let dz = layer.thickness_m / h;
        let q = layer.q.map(|v| v / T::from_f64(q_scale));
        let g = layer.g.map(|x| x / T::from_f64(g_scale));
        for k in 0..4 {
            raw[k] = raw[k] + g[k] * T::from_f64(dz);
        }
        let l00 = g[0].sqrt();
        let l10 = g[2] / l00;
        let pivot = g[3] - l10 * l10;
        if !(l00.value().is_finite()
            && l00.value() > 0.0
            && pivot.value().is_finite()
            && pivot.value() > 0.0)
        {
            return Err(invalid(
                "equilibrated layer shear compliance is not positive definite",
            ));
        }
        let l11 = pivot.sqrt();
        let qn: [T; 9] = core::array::from_fn(|ij| {
            (0..3).fold(T::zero(), |sum, k| {
                sum + q[3 * (ij / 3) + k] * neutral[3 * k + ij % 3]
            })
        });
        let integral = |delta: f64| -> [T; 9] {
            core::array::from_fn(|k| {
                primitive[k] + T::from_f64(delta) * (q[k] * T::from_f64(z0 + 0.5 * delta) - qn[k])
            })
        };
        // The primitive is quadratic, hence positive three-point Gauss
        // integration is exact for its quartic energy within each layer.
        for (node, weight) in [(-gauss, 5.0 / 9.0), (0.0, 8.0 / 9.0), (gauss, 5.0 / 9.0)] {
            let traction = divergence(&integral(0.5 * dz * (1.0 + node)), axes).map(|x| -x);
            let mut whitened = [T::zero(); 4];
            for j in 0..2 {
                let tx = traction[0] * inverse[j] + traction[1] * inverse[2 + j];
                let ty = traction[2] * inverse[j] + traction[3] * inverse[2 + j];
                whitened[j] = tx / l00;
                whitened[2 + j] = (ty - l10 * whitened[j]) / l11;
            }
            for i in 0..2 {
                for j in i..2 {
                    compliance[2 * i + j] = compliance[2 * i + j]
                        + T::from_f64(0.5 * dz * weight)
                            * (whitened[i] * whitened[j] + whitened[2 + i] * whitened[2 + j]);
                }
            }
        }
        primitive = integral(dz);
        z0 += dz;
    }
    if primitive
        .iter()
        .any(|x| !x.value().is_finite() || x.value().abs() > 8192.0 * f64::EPSILON)
    {
        return Err(invalid(
            "equilibrated shear does not close its free-face traction",
        ));
    }
    compliance[2] = compliance[1];
    let c_scale = compliance
        .iter()
        .map(|x| x.value().abs())
        .fold(0.0_f64, f64::max);
    let c = compliance.map(|x| x / T::from_f64(c_scale));
    let det = c[0] * c[3] - c[1] * c[2];
    if !(c[0].value() > 0.0 && c[3].value() > 0.0 && det.value() > 0.0 && det.value().is_finite()) {
        return Err(invalid(
            "equilibrated shear energy is not positive definite",
        ));
    }
    let normalized = [c[3], -c[1], -c[2], c[0]].map(|x| x / det / T::from_f64(c_scale));
    let stiffness = normalized.map(|x| (x * T::from_f64(g_scale)) * T::from_f64(h));
    let factors = [normalized[0] / raw[0], normalized[3] / raw[3]];
    if stiffness
        .iter()
        .chain(&factors)
        .any(|x| !x.value().is_finite())
        || stiffness[0].value() <= 0.0
        || stiffness[3].value() <= 0.0
        || factors[0].value() <= 0.0
        || factors[1].value() <= 0.0
    {
        return Err(invalid("equilibrated shear stiffness is unrepresentable"));
    }
    check()?;
    Ok((stiffness, factors))
}

// Fixed 3x3 Cholesky and three solves. This is fs-la's small-panel operation
// order lifted to Real, so derivatives include the full neutral-axis A^-1 B.
// The arrays stay on the stack; no matrix inverse or generic solver is exposed.
fn neutral_matrix<T: Real, E>(
    mut lower: [T; 9],
    b: [T; 9],
    invalid: &impl Fn(&'static str) -> E,
) -> Result<[T; 9], E> {
    for j in 0..3 {
        let mut diagonal = lower[3 * j + j];
        for k in 0..j {
            diagonal = (-lower[3 * j + k]).mul_add(lower[3 * j + k], diagonal);
        }
        if !diagonal.value().is_finite() || diagonal.value() <= 0.0 {
            return Err(invalid(
                "equilibrated shear membrane stiffness is not positive definite",
            ));
        }
        let pivot = diagonal.sqrt();
        lower[3 * j + j] = pivot;
        for i in j + 1..3 {
            let mut value = lower[3 * i + j];
            for k in 0..j {
                value = (-lower[3 * i + k]).mul_add(lower[3 * j + k], value);
            }
            lower[3 * i + j] = value / pivot;
        }
    }
    let mut neutral = [T::zero(); 9];
    for j in 0..3 {
        let mut rhs = [b[j], b[3 + j], b[6 + j]];
        for i in 0..3 {
            let mut value = rhs[i];
            for (k, &x) in rhs.iter().enumerate().take(i) {
                value = (-lower[3 * i + k]).mul_add(x, value);
            }
            rhs[i] = value / lower[3 * i + i];
        }
        for i in (0..3).rev() {
            let mut value = rhs[i];
            for (k, &x) in rhs.iter().enumerate().skip(i + 1) {
                value = (-lower[3 * k + i]).mul_add(x, value);
            }
            rhs[i] = value / lower[3 * i + i];
        }
        for i in 0..3 {
            neutral[3 * i + j] = rhs[i];
        }
    }
    Ok(neutral)
}
