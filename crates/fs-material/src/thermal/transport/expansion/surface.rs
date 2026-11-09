//! Material-face areas from the original free-expansion geometry.
use super::{FreeExpansionGeometrySample, FreeExpansionGeometryTangent, PreparedFreeExpansion};
use crate::thermal::transport::ThermalTransportError;

/// Bounds the work and memory of one oriented surface evaluation. Coplanar
/// patches with the same boundary condition may be combined before admission.
pub const MAX_ORIENTED_EXPANSION_PATCHES: usize = 4096;

/// One explicitly selected set of material faces, by reference principal normal.
/// Areas can combine both signs of the same normal when their boundary law and
/// temperature are shared. Distinct ambient, reservoir and contact patches need
/// separate bindings and an owner that admits their disjoint surface coverage.
///
/// The original expansion path supplies every stretch and source identity.
/// For a face normal to i, A_i(T)=A_i,ref*s_j(T)*s_k(T). No body shape, face area,
/// isotropy, normal rotation, shear, or heat-transfer coefficient is inferred
/// from a resistor's length/section or from the material's density. These are
/// material faces under unconstrained diagonal free expansion. Oblique planar
/// patches use the exact cofactor map of the same deformation gradient; a
/// curved surface needs an explicitly supplied facet/quadrature approximation.
/// This is not a constrained thermoelastic area calculation.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedExpansionSurface {
    source: PreparedFreeExpansion,
    reference_areas_m2: [f64; 3],
    reference_oriented_area_vectors_m2: Vec<[f64; 3]>,
    reference_total_area_m2: f64,
}

/// Areas and exact temperature tangents from that same original source sample.
/// Source receipts and the complete fixed context remain on `source()`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExpansionSurfaceSample {
    geometry: FreeExpansionGeometrySample,
    principal_areas_m2: [f64; 3],
    principal_area_derivatives_m2_k: [f64; 3],
    total_area_m2: f64,
    total_area_derivative_m2_k: f64,
}

impl PreparedExpansionSurface {
    /// Reference areas at this path's exact reference query point. Zero means
    /// no selected faces with that normal; an entirely empty patch is explicit.
    /// No source selection, database query or thermal state is changed.
    pub fn try_new(
        source: PreparedFreeExpansion,
        reference_areas_m2: [f64; 3],
    ) -> Result<Self, ThermalTransportError> {
        if !reference_areas_m2
            .into_iter()
            .all(|a| a.is_finite() && a >= 0.0)
        {
            return Err(invalid(
                "reference surface areas must be finite and nonnegative",
            ));
        }
        let reference_total_area_m2 = reference_areas_m2.into_iter().sum::<f64>();
        if !reference_total_area_m2.is_finite() {
            return Err(invalid("reference total surface area is unrepresentable"));
        }
        Ok(Self {
            source,
            reference_areas_m2,
            reference_oriented_area_vectors_m2: Vec::new(),
            reference_total_area_m2,
        })
    }

    /// Add explicitly authored planar surface patches in the source's fixed
    /// principal frame. Each vector is reference area times its unit normal
    /// [m²]; its norm is the actual area, not the sum of projected face areas.
    /// Opposite signs describe the same area response. Patches must be disjoint
    /// from each other and the principal faces; coverage is a caller obligation.
    ///
    /// Nanson's formula gives `a(T) = cof(F(T)) a_ref`, with
    /// `F = diag(s_1,s_2,s_3)`. Thus an oblique face changes both its area and
    /// normal under anisotropic expansion. No isotropic area approximation or
    /// source query is introduced. The existing thermal-boundary consumers use
    /// the combined total area and tangent without a separate heat source.
    pub fn with_oriented_patches(
        mut self,
        reference_area_vectors_m2: &[[f64; 3]],
    ) -> Result<Self, ThermalTransportError> {
        if reference_area_vectors_m2.len()
            > MAX_ORIENTED_EXPANSION_PATCHES - self.reference_oriented_area_vectors_m2.len()
        {
            return Err(invalid("oriented expansion surface exceeds its patch bound"));
        }
        let mut total = self.reference_total_area_m2;
        for &area in reference_area_vectors_m2 {
            total += vector_area(area)?;
            if !total.is_finite() {
                return Err(invalid("reference total surface area is unrepresentable"));
            }
        }
        self.reference_oriented_area_vectors_m2
            .extend_from_slice(reference_area_vectors_m2);
        self.reference_total_area_m2 = total;
        Ok(self)
    }

    /// Exact original source/card/claim/reference/basis binding.
    pub const fn source(&self) -> &PreparedFreeExpansion {
        &self.source
    }
    /// Only the principal-normal portion; oriented patches contribute to the
    /// total separately and are never counted as projected principal areas.
    pub const fn reference_areas_m2(&self) -> [f64; 3] {
        self.reference_areas_m2
    }
    pub fn reference_oriented_area_vectors_m2(&self) -> &[[f64; 3]] {
        &self.reference_oriented_area_vectors_m2
    }
    pub const fn reference_total_area_m2(&self) -> f64 {
        self.reference_total_area_m2
    }

    /// One original source interpolation, then three principal-face evaluations
    /// and bounded linear work in the explicitly admitted oriented patches.
    /// Successful runtime evaluation allocates nothing and performs no query.
    /// One-sided derivatives at source knots follow the original path exactly.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ExpansionSurfaceSample, ThermalTransportError> {
        self.from_geometry(self.source.at_temperature_with_tangent(temperature_k)?)
    }

    /// Reuse a tangent already evaluated by another consumer of the SAME path.
    /// Numeric agreement from another source is insufficient; its identity must
    /// match, including reference normalization and all fixed query conditions.
    pub fn from_geometry(
        &self,
        tangent: FreeExpansionGeometryTangent,
    ) -> Result<ExpansionSurfaceSample, ThermalTransportError> {
        let geometry = tangent.geometry();
        if geometry.path_identity() != self.source.identity() {
            return Err(invalid(
                "surface and supplied geometry have different expansion sources",
            ));
        }
        let s = geometry.principal_stretches();
        let ds = tangent.principal_stretch_derivatives_per_k();
        let mut areas = [0.0; 3];
        let mut derivatives = [0.0; 3];
        for i in 0..3 {
            let [j, k] = [(i + 1) % 3, (i + 2) % 3];
            if self.reference_areas_m2[i] == 0.0 {
                continue;
            }
            areas[i] = (self.reference_areas_m2[i] * s[j]) * s[k];
            let relative = ds[j] / s[j] + ds[k] / s[k];
            derivatives[i] = areas[i] * relative;
            if !areas[i].is_finite()
                || areas[i] <= 0.0
                || !relative.is_finite()
                || !derivatives[i].is_finite()
                || (relative != 0.0 && derivatives[i] == 0.0)
            {
                return Err(invalid(
                    "expanded surface area or its tangent is unrepresentable",
                ));
            }
        }
        let mut total_area_m2 = areas.into_iter().sum::<f64>();
        let mut total_area_derivative_m2_k = derivatives.into_iter().sum::<f64>();
        let oriented_relative: [f64; 3] = core::array::from_fn(|i| {
            let [j, k] = [(i + 1) % 3, (i + 2) % 3];
            ds[j] / s[j] + ds[k] / s[k]
        });
        for reference in &self.reference_oriented_area_vectors_m2 {
            let mut vector = [0.0; 3];
            let mut vector_derivative = [0.0; 3];
            for i in 0..3 {
                if reference[i] == 0.0 {
                    continue;
                }
                let [j, k] = [(i + 1) % 3, (i + 2) % 3];
                vector[i] = (reference[i] * s[j]) * s[k];
                let relative = oriented_relative[i];
                vector_derivative[i] = vector[i] * relative;
                if !vector[i].is_finite()
                    || vector[i] == 0.0
                    || !relative.is_finite()
                    || !vector_derivative[i].is_finite()
                    || (relative != 0.0 && vector_derivative[i] == 0.0)
                {
                    return Err(invalid("oriented surface cofactor or tangent is unrepresentable"));
                }
            }
            let area = vector_area(vector)?;
            // d|a|/dT = n(T) dot da/dT. Normalize before multiplication
            // to avoid squaring large, but representable, areas.
            let derivative = (0..3)
                .map(|i| (vector[i] / area) * vector_derivative[i])
                .sum::<f64>();
            total_area_m2 += area;
            total_area_derivative_m2_k += derivative;
        }
        if !total_area_m2.is_finite() || !total_area_derivative_m2_k.is_finite() {
            return Err(invalid(
                "total expanded surface area or its tangent is unrepresentable",
            ));
        }
        Ok(ExpansionSurfaceSample {
            geometry,
            principal_areas_m2: areas,
            principal_area_derivatives_m2_k: derivatives,
            total_area_m2,
            total_area_derivative_m2_k,
        })
    }
}
impl ExpansionSurfaceSample {
    pub const fn geometry(self) -> FreeExpansionGeometrySample {
        self.geometry
    }
    /// Areas of the principal-normal portion only. Oriented patches are
    /// included in `total_area_m2`, not projected into these entries.
    pub const fn principal_areas_m2(self) -> [f64; 3] {
        self.principal_areas_m2
    }
    pub const fn principal_area_derivatives_m2_k(self) -> [f64; 3] {
        self.principal_area_derivatives_m2_k
    }
    pub const fn total_area_m2(self) -> f64 {
        self.total_area_m2
    }
    pub const fn total_area_derivative_m2_k(self) -> f64 {
        self.total_area_derivative_m2_k
    }
}
fn vector_area(vector: [f64; 3]) -> Result<f64, ThermalTransportError> {
    if !vector.into_iter().all(f64::is_finite) {
        return Err(invalid("reference area vectors must be finite"));
    }
    let scale = vector.into_iter().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 {
        return Err(invalid("oriented patches require a nonzero reference area vector"));
    }
    let norm = scale * vector.into_iter().map(|v| (v / scale).powi(2)).sum::<f64>().sqrt();
    if !norm.is_finite() || norm <= 0.0 {
        return Err(invalid("oriented surface area is unrepresentable"));
    }
    Ok(norm)
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
