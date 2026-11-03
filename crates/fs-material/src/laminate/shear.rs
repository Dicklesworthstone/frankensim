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
    /// edge accuracy or viscoelastic/temperature derivatives. The original
    /// [`Self::try_new`] retains its differentiated directional Whitney model.
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
        let (stiffness, factors) = equilibrated(&section, bending_axis_rad)?;
        section.shear_pa_m = stiffness;
        section.shear_factor = factors;
        section.shear_model = LaminateShearModel::EquilibratedCylindrical { bending_axis_rad };
        Ok(section)
    }
}

/// Columns are the shear resultants/tractions of the two physical probes.
fn divergence(matrix: &[f64; 9], axes: [[f64; 2]; 2]) -> [f64; 4] {
    let mut out = [0.0; 4];
    for (column, n) in axes.iter().enumerate() {
        let curvature = [n[0] * n[0], n[1] * n[1], 2.0 * n[0] * n[1]];
        let stress: [f64; 3] =
            core::array::from_fn(|i| (0..3).map(|j| matrix[3 * i + j] * curvature[j]).sum());
        out[column] = stress[0] * n[0] + stress[2] * n[1];
        out[2 + column] = stress[2] * n[0] + stress[1] * n[1];
    }
    out
}

fn equilibrated(
    section: &LaminateSection,
    angle: f64,
) -> Result<([f64; 4], [f64; 2]), LaminateError> {
    // Dimensionless thickness and moduli avoid powers of SI scales and any
    // cancellation from the caller's choice of an offset reference surface.
    let h = section.total_thickness_m;
    let q_scale = section
        .integrals
        .iter()
        .flat_map(|i| i.q_bar_pa)
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let mut q_layers = Vec::with_capacity(section.layers.len());
    let mut g_layers = Vec::with_capacity(section.layers.len());
    let (mut a, mut b, mut d) = ([0.0; 9], [0.0; 9], [0.0; 9]);
    let (mut z0, mut g_scale) = (-0.5, 0.0_f64);
    for (i, (layer, geometry)) in section.layers.iter().zip(&section.integrals).enumerate() {
        let q = geometry.q_bar_pa.map(|v| v / q_scale);
        let dz = layer.thickness_m / h;
        let z1 = z0 + dz;
        let moments = [
            dz,
            0.5 * dz * (z0 + z1),
            dz * (z0 * z0 + z0 * z1 + z1 * z1) / 3.0,
        ];
        for k in 0..9 {
            a[k] += q[k] * moments[0];
            b[k] += q[k] * moments[1];
            d[k] += q[k] * moments[2];
        }
        let (_, principal) = principal_stiffness(&layer.law, i)?;
        let g = rotate_transverse(principal, layer.orientation_rad);
        g_scale = g.iter().map(|v| v.abs()).fold(g_scale, f64::max);
        q_layers.push(q);
        g_layers.push(g);
        z0 = z1;
    }
    let factor = fs_la::factor::cholesky(&a, 3)
        .map_err(|_| invalid("equilibrated shear membrane stiffness is not positive definite"))?;
    let mut neutral = [0.0; 9];
    for j in 0..3 {
        let mut rhs = [b[j], b[3 + j], b[6 + j]];
        factor.solve(&mut rhs);
        for i in 0..3 {
            neutral[3 * i + j] = rhs[i];
        }
    }
    let relaxed: [f64; 9] = core::array::from_fn(|ij| {
        d[ij]
            - (0..3)
                .map(|k| b[3 * (ij / 3) + k] * neutral[3 * k + ij % 3])
                .sum::<f64>()
    });
    let (sin, cos) = angle.sin_cos();
    let axes = [[cos, sin], [-sin, cos]];
    let v = divergence(&relaxed, axes);
    let v_scale = v.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let scaled = v.map(|x| x / v_scale);
    let determinant = scaled[0] * scaled[3] - scaled[1] * scaled[2];
    if !determinant.is_finite() || determinant.abs() <= 128.0 * f64::EPSILON {
        return Err(invalid(
            "physical bending probes do not resolve two independent shear resultants",
        ));
    }
    let inverse = [scaled[3], -scaled[1], -scaled[2], scaled[0]].map(|x| x / determinant / v_scale);

    let (mut primitive, mut compliance, mut raw) = ([0.0; 9], [0.0; 4], [0.0; 4]);
    z0 = -0.5;
    let gauss = (3.0_f64 / 5.0).sqrt();
    for ((layer, q), g) in section.layers.iter().zip(&q_layers).zip(&g_layers) {
        let dz = layer.thickness_m / h;
        let g = g.map(|x| x / g_scale);
        for k in 0..4 {
            raw[k] += g[k] * dz;
        }
        let l00 = g[0].sqrt();
        let l10 = g[2] / l00;
        let pivot = g[3] - l10 * l10;
        if !(l00.is_finite() && l00 > 0.0 && pivot.is_finite() && pivot > 0.0) {
            return Err(invalid(
                "equilibrated layer shear compliance is not positive definite",
            ));
        }
        let l11 = pivot.sqrt();
        let qn: [f64; 9] = core::array::from_fn(|ij| {
            (0..3)
                .map(|k| q[3 * (ij / 3) + k] * neutral[3 * k + ij % 3])
                .sum()
        });
        let integral = |delta: f64| -> [f64; 9] {
            core::array::from_fn(|k| primitive[k] + delta * (q[k] * (z0 + 0.5 * delta) - qn[k]))
        };
        // The primitive is quadratic, hence positive three-point Gauss
        // integration is exact for its quartic energy within each layer.
        for (node, weight) in [(-gauss, 5.0 / 9.0), (0.0, 8.0 / 9.0), (gauss, 5.0 / 9.0)] {
            let traction = divergence(&integral(0.5 * dz * (1.0 + node)), axes).map(|x| -x);
            let mut whitened = [0.0; 4];
            for j in 0..2 {
                let tx = traction[0] * inverse[j] + traction[1] * inverse[2 + j];
                let ty = traction[2] * inverse[j] + traction[3] * inverse[2 + j];
                whitened[j] = tx / l00;
                whitened[2 + j] = (ty - l10 * whitened[j]) / l11;
            }
            for i in 0..2 {
                for j in i..2 {
                    compliance[2 * i + j] += 0.5
                        * dz
                        * weight
                        * (whitened[i] * whitened[j] + whitened[2 + i] * whitened[2 + j]);
                }
            }
        }
        primitive = integral(dz);
        z0 += dz;
    }
    if primitive
        .iter()
        .any(|x| !x.is_finite() || x.abs() > 8192.0 * f64::EPSILON)
    {
        return Err(invalid(
            "equilibrated shear does not close its free-face traction",
        ));
    }
    compliance[2] = compliance[1];
    let c_scale = compliance
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let c = compliance.map(|x| x / c_scale);
    let det = c[0] * c[3] - c[1] * c[2];
    if !(c[0] > 0.0 && c[3] > 0.0 && det > 0.0 && det.is_finite()) {
        return Err(invalid(
            "equilibrated shear energy is not positive definite",
        ));
    }
    let normalized = [c[3], -c[1], -c[2], c[0]].map(|x| x / det / c_scale);
    let stiffness = normalized.map(|x| (x * g_scale) * h);
    let factors = [normalized[0] / raw[0], normalized[3] / raw[3]];
    if stiffness.iter().chain(&factors).any(|x| !x.is_finite())
        || stiffness[0] <= 0.0
        || stiffness[3] <= 0.0
        || factors[0] <= 0.0
        || factors[1] <= 0.0
    {
        return Err(invalid("equilibrated shear stiffness is unrepresentable"));
    }
    Ok((stiffness, factors))
}
