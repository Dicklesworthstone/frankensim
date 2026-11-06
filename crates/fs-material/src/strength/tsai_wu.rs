//! Three-dimensional orthotropic Tsai--Wu strength in fixed material axes.
//!
//! Normalization uses `u_i = sigma_ii / sqrt(Xti Xci)` and
//! `v_ij = tau_ij / Sij`. With explicitly sourced normalized interactions
//! `f12,f23,f31`, the failure polynomial is `A+B`, where
//! `A = u^T C u + sum(v^2)`, `C_ii=1`, `C_ij=fij`, and
//! `B = sum(sigma_ii * (1/Xti - 1/Xci))`.
//!
//! All nine strengths and all three interactions must be supplied. In
//! particular, zero interaction is not inferred from missing biaxial data.
//! The complete normal interaction matrix must be positive definite; checking
//! its three coefficients independently is insufficient in three dimensions.
//! This selects a bounded convex strength surface. Open/degenerate surfaces
//! and determinants unresolved at the documented roundoff margin refuse here.
//! Source coefficients are never clipped or automatically regularized.
//!
//! The positive homogeneous utilization is `(B + sqrt(B^2 + 4 A))/2`.
//! Scaling ALL stresses by its reciprocal reaches the surface; this is not a
//! load factor with residual or thermal stress held fixed. A Cholesky norm,
//! normalized stress scale and rationalized negative-B root avoid squaring
//! dimensional stresses or subtracting nearly equal roots. Prepared point
//! assessment allocates nothing and performs no eigensolve or source lookup.
//!
//! Reference: Tsai and Wu (1971), doi:10.1177/002199837100500106; the 3-D
//! stress polynomial is also given in SIMULIA, Generalized Damage Initiation
//! Criteria for Anisotropic Materials:
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-damageinitgeneral.htm>.
//! This is a phenomenological strength observation, not a failure-mode
//! decomposition, progressive damage law, fracture energy or calibration.

use super::invalid;
use crate::Voigt;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_math::det::{hypot, sqrt};
use fs_qty::{Dims, Pressure, QuantitySpec};

/// Original source strengths and interactions for a full 3-D orthotropic body.
#[derive(Clone, Debug, PartialEq)]
pub struct OrthotropicTsaiWu {
    resolved: ResolvedMaterialStatePoint,
    strength_properties: [String; 9],
    interaction_properties: [String; 3],
    strengths_pa: [f64; 9],
    interactions: [f64; 3],
    scales_pa: [f64; 6],
    linear: [f64; 3],
    // Upper Cholesky factor: C = R^T R, fixed at preparation.
    factor: [[f64; 3]; 3],
}

impl OrthotropicTsaiWu {
    /// Resolve strengths in `[X1t,X1c,X2t,X2c,X3t,X3c,S12,S23,S31]`
    /// order, all positive magnitudes in Pa. Interaction names are in
    /// `[f12,f23,f31]` order, where `fij = Fij / sqrt(Fii Fjj)`.
    /// Dimensional inverse-pressure-squared coefficients are not these inputs.
    /// All twelve distinct claims use the same complete physical query.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        strength_properties: [&str; 9],
        interaction_properties: [&str; 3],
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        let mut requirements = Vec::with_capacity(12);
        for name in strength_properties {
            requirements.push(ScalarPropertyRequirement::try_new(
                name,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
        for name in interaction_properties {
            requirements.push(ScalarPropertyRequirement::try_new(
                name,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: -1.0,
                    upper: 1.0,
                },
            )?);
        }
        Self::try_from_resolved(
            resolve_material_state_point(card, point, &requirements, selection)?,
            strength_properties,
            interaction_properties,
        )
    }

    /// Bind existing source answers without selecting replacements, filling
    /// missing directions or equating compressive and tensile strengths.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        strength_properties: [&str; 9],
        interaction_properties: [&str; 3],
    ) -> Result<Self, MaterialStatePointError> {
        let names: [&str; 12] = core::array::from_fn(|i| {
            if i < 9 {
                strength_properties[i]
            } else {
                interaction_properties[i - 9]
            }
        });
        for i in 0..names.len() {
            if names[..i].contains(&names[i]) {
                return Err(invalid(
                    "3-D Tsai-Wu requires twelve distinct named source claims",
                ));
            }
        }
        let mut values = [0.0; 12];
        for (i, name) in names.iter().enumerate() {
            let source = resolved.property(name).ok_or_else(|| {
                invalid("3-D Tsai-Wu requires every selected strength and interaction")
            })?;
            let dims = if i < 9 { Pressure::DIMS } else { Dims::NONE };
            let value = source.value_si();
            if source.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !source.is_static_at_fixed_external_conditions()
                || !value.is_finite()
                || (i < 9 && value <= 0.0)
                || (i >= 9 && value.abs() >= 1.0)
            {
                return Err(invalid(
                    "3-D Tsai-Wu needs fixed-condition positive Pa strengths and normalized dimensionless interactions in (-1,1)",
                ));
            }
            values[i] = value;
        }
        let strengths_pa = core::array::from_fn(|i| values[i]);
        let interactions = [values[9], values[10], values[11]];
        let [f12, f23, f31] = interactions;
        // A positive rounded Cholesky pivot alone can admit an indefinite
        // source matrix arbitrarily close to singularity. Bound cancellation
        // in the normalized determinant by a conservative arithmetic margin,
        // before factorizing. This rejects unresolved input, not repairs it.
        let triple = 2.0 * f12 * f23 * f31;
        let squares = f12 * f12 + f23 * f23 + f31 * f31;
        let determinant = (1.0 + triple) - squares;
        let roundoff_margin = 64.0 * f64::EPSILON * (1.0 + triple.abs() + squares);
        if determinant <= roundoff_margin {
            return Err(invalid(
                "3-D Tsai-Wu interaction determinant must be positive beyond its roundoff margin",
            ));
        }
        let r11 = sqrt((1.0 - f12) * (1.0 + f12));
        let r12 = (f23 - f12 * f31) / r11;
        let last_pivot = (1.0 - f31) * (1.0 + f31) - r12 * r12;
        if !last_pivot.is_finite() || last_pivot <= 0.0 {
            return Err(invalid(
                "complete 3-D Tsai-Wu interaction matrix must be positive definite",
            ));
        }
        let factor = [
            [1.0, f12, f31],
            [0.0, r11, r12],
            [0.0, 0.0, sqrt(last_pivot)],
        ];
        let scales_pa: [f64; 6] = core::array::from_fn(|i| {
            if i < 3 {
                sqrt(strengths_pa[2 * i]) * sqrt(strengths_pa[2 * i + 1])
            } else {
                strengths_pa[i + 3]
            }
        });
        let linear = core::array::from_fn(|i| {
            scales_pa[i] / strengths_pa[2 * i] - scales_pa[i] / strengths_pa[2 * i + 1]
        });
        if scales_pa.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || linear.iter().any(|v: &f64| !v.is_finite())
        {
            return Err(invalid(
                "3-D Tsai-Wu strength normalization is unrepresentable",
            ));
        }
        Ok(Self {
            resolved,
            strength_properties: strength_properties.map(str::to_owned),
            interaction_properties: interaction_properties.map(str::to_owned),
            strengths_pa,
            interactions,
            scales_pa,
            linear,
            factor,
        })
    }

    /// Refresh exactly the original selections at supported external conditions.
    /// Unsupported queries, a foreign card or loss of convexity refuse atomically.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(
            self.resolved.at_query_point(card, point)?,
            self.strength_properties(),
            self.interaction_properties(),
        )
    }

    /// Complete retained source state, including its typed conditions.
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }
    /// Ordered source names for six normal and three shear strengths.
    pub fn strength_properties(&self) -> [&str; 9] {
        self.strength_properties.each_ref().map(String::as_str)
    }
    /// Ordered names for independently selected `f12,f23,f31`.
    pub fn interaction_properties(&self) -> [&str; 3] {
        self.interaction_properties.each_ref().map(String::as_str)
    }
    /// `[X1t,X1c,X2t,X2c,X3t,X3c,S12,S23,S31]`, all positive Pa.
    pub const fn strengths_pa(&self) -> [f64; 9] {
        self.strengths_pa
    }
    /// Dimensionless `[f12,f23,f31]`, not dimensional polynomial coefficients.
    pub const fn normalized_interactions(&self) -> [f64; 3] {
        self.interactions
    }

    /// Assess tensile-positive tensor-shear `[11,22,33,12,23,31]` stress [Pa]
    /// in the ORIGINAL material axes. A consumer must use the same material
    /// rotation as its stiffness and compatible stress/strength conventions.
    pub fn assess_material_stress(
        &self,
        stress_pa: &Voigt,
    ) -> Result<OrthotropicTsaiWuAssessment, MaterialStatePointError> {
        let normalized: [f64; 6] = core::array::from_fn(|i| stress_pa[i] / self.scales_pa[i]);
        if stress_pa.iter().any(|v| !v.is_finite())
            || normalized
                .iter()
                .enumerate()
                .any(|(i, v)| !v.is_finite() || (*v == 0.0 && stress_pa[i] != 0.0))
        {
            return Err(invalid(
                "3-D Tsai-Wu requires finite representable material-axis stresses",
            ));
        }
        let scale = normalized.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Ok(OrthotropicTsaiWuAssessment {
                material_stress_pa: *stress_pa,
                utilization: 0.0,
            });
        }
        let unit = normalized.map(|v| v / scale);
        if unit
            .iter()
            .zip(normalized)
            .any(|(u, v)| *u == 0.0 && v != 0.0)
        {
            return Err(invalid(
                "3-D Tsai-Wu stress dynamic range is unrepresentable",
            ));
        }
        let mut q = 0.0;
        for row in self.factor {
            q = hypot(q, row[0] * unit[0] + row[1] * unit[1] + row[2] * unit[2]);
        }
        for shear in &unit[3..] {
            q = hypot(q, *shear);
        }
        let b = self.linear[0] * unit[0] + self.linear[1] * unit[1] + self.linear[2] * unit[2];
        let radius = hypot(0.5 * b, q);
        let radial = if b >= 0.0 {
            radius + 0.5 * b
        } else {
            q * (q / (radius - 0.5 * b))
        };
        let utilization = scale * radial;
        if !b.is_finite() || !utilization.is_finite() || utilization <= 0.0 {
            return Err(invalid("3-D Tsai-Wu utilization is unrepresentable"));
        }
        Ok(OrthotropicTsaiWuAssessment {
            material_stress_pa: *stress_pa,
            utilization,
        })
    }
}

/// Full coupled strength observation; no isolated failure mode is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthotropicTsaiWuAssessment {
    material_stress_pa: Voigt,
    utilization: f64,
}
impl OrthotropicTsaiWuAssessment {
    /// Actual tensor-shear stress in the admitted material frame [Pa].
    pub const fn material_stress_pa(&self) -> Voigt {
        self.material_stress_pa
    }
    /// Positive homogeneous distance to the complete strength surface.
    pub const fn utilization(&self) -> f64 {
        self.utilization
    }
    /// Whether the selected interaction criterion reaches its limit.
    pub fn reaches_limit(&self) -> bool {
        self.utilization >= 1.0
    }
}
