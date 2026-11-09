//! Source-bound Tsai-Wu interaction in the original lamina material axes.
//!
//! The five positive strength magnitudes retain the parent module's convention.
//! With `u=s11/sqrt(Xt Xc)`, `v=s22/sqrt(Yt Yc)`, `w=t12/S`, the failure
//! polynomial is `A+B`, where `A=u*u+v*v+2*f*u*v+w*w` and
//! `B=s11*(1/Xt-1/Xc)+s22*(1/Yt-1/Yc)`. The dimensionless interaction `f`
//! is a separately selected source claim from the SAME card and physical state.
//! No biaxial calibration is inferred from uniaxial strengths; even zero
//! interaction must be explicitly supplied. We require `-1 < f < 1` for a
//! strictly convex, bounded surface surrounding zero stress.
//!
//! The reported utilization is `R=(B+sqrt(B*B+4*A))/2`, so `stress/R` lies
//! on the failure surface and `R>=1` indicates a crossing. This positive,
//! homogeneous measure avoids ranking negative polynomial values as safety
//! margins. The reciprocal scales the ENTIRE current stress vector, including
//! residual/thermal stress; it is not a mechanical-load factor with fixed
//! thermal stress. A scaled quadratic norm and stable negative-B branch avoid
//! squared-stress overflow and cancellation. No hot allocation or source query.
//!
//! Native preparation retains both faces of every original layer. The convex
//! utilization of affine through-thickness stress reaches its maximum at a
//! face. This criterion provides no failure-mode decomposition, degradation,
//! delamination, fracture energy or post-failure redistribution.
//!
//! Reference: Abaqus, Plane Stress Orthotropic Failure Measures, Tsai-Wu theory
//! and output failure indices (its compressive inputs are negative; ours are
//! positive magnitudes):
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-failuremeasures.htm>.

use super::{
    LaminaFace, LaminaMaximumStress, LaminateStrengthError, PreparedLaminateStrength, admit_beam,
    invalid, rotate,
};
use crate::laminate::thermal::source::{LaminateThermalSourceError, ResolvedLaminateThermal};
use crate::laminate::thermal::{LayerThermalStress, beam::PreparedLaminateBeam};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_math::det::{hypot, sqrt};
use fs_qty::{Dims, QuantitySpec};

/// Five original strengths plus their independently sourced biaxial interaction.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminaTsaiWu {
    strength: LaminaMaximumStress,
    interaction: ResolvedMaterialStatePoint,
    interaction_property: String,
    coefficient: f64,
    scales_pa: [f64; 3],
    linear: [f64; 2],
    transverse_factor: f64,
}

impl LaminaTsaiWu {
    /// Resolve the interaction at the strengths' exact complete typed query.
    /// The selected source must explicitly describe normalized `F12`, not an
    /// unnormalized coefficient with inverse-square pressure units.
    pub fn resolve(
        card: &MaterialCard,
        strength: LaminaMaximumStress,
        interaction_property: impl Into<String>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        let interaction_property = interaction_property.into();
        let requirement = ScalarPropertyRequirement::try_new(
            &interaction_property,
            Dims::NONE,
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 1.0,
            },
        )?;
        let interaction = resolve_material_state_point(
            card,
            &strength.resolved().query(),
            &[requirement],
            selection,
        )?;
        Self::try_from_resolved(strength, interaction, interaction_property)
    }

    /// Bind existing source answers without losing their original selections.
    /// Foreign cards, different typed conditions and unsupported coefficients
    /// refuse before any criterion is published.
    pub fn try_from_resolved(
        strength: LaminaMaximumStress,
        interaction: ResolvedMaterialStatePoint,
        interaction_property: impl Into<String>,
    ) -> Result<Self, MaterialStatePointError> {
        let interaction_property = interaction_property.into();
        if strength.resolved().card_identity() != interaction.card_identity()
            || strength.resolved().query() != interaction.query()
        {
            return Err(invalid(
                "Tsai-Wu interaction and strengths require the same card and typed physical query",
            ));
        }
        let property = interaction
            .property(&interaction_property)
            .ok_or_else(|| invalid("Tsai-Wu requires its explicitly selected interaction claim"))?;
        let coefficient = property.value_si();
        if property.requirement().quantity() != QuantitySpec::dimensional(Dims::NONE)
            || !property.is_static_at_fixed_external_conditions()
            || !coefficient.is_finite()
            || coefficient.abs() >= 1.0
        {
            return Err(invalid(
                "normalized Tsai-Wu interaction must be static dimensionless and strictly between -1 and 1",
            ));
        }
        let [xt, xc, yt, yc, shear] = strength.limits_pa();
        // Preparation avoids multiplying two strengths, which can overflow
        // even though their geometric mean and all utilizations are finite.
        let scales_pa = [sqrt(xt) * sqrt(xc), sqrt(yt) * sqrt(yc), shear];
        let linear = [
            scales_pa[0] / xt - scales_pa[0] / xc,
            scales_pa[1] / yt - scales_pa[1] / yc,
        ];
        let transverse_factor = sqrt((1.0 - coefficient) * (1.0 + coefficient));
        if scales_pa.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || linear.iter().any(|v| !v.is_finite())
            || !transverse_factor.is_finite()
            || transverse_factor <= 0.0
        {
            return Err(invalid("Tsai-Wu strength normalization is unrepresentable"));
        }
        Ok(Self {
            strength,
            interaction,
            interaction_property,
            coefficient,
            scales_pa,
            linear,
            transverse_factor,
        })
    }

    /// Refresh every original claim together; no new claim is selected.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(
            self.strength.at_query_point(card, point)?,
            self.interaction.at_query_point(card, point)?,
            self.interaction_property.clone(),
        )
    }

    /// Original five directional strength sources and maximum-stress criterion.
    pub const fn strength(&self) -> &LaminaMaximumStress {
        &self.strength
    }
    /// Original interaction claim, material identity and complete typed query.
    pub const fn interaction_source(&self) -> &ResolvedMaterialStatePoint {
        &self.interaction
    }
    /// Explicit source property interpreted as normalized `F12`.
    pub fn interaction_property(&self) -> &str {
        &self.interaction_property
    }
    /// The selected normalized interaction, strictly between -1 and 1.
    pub const fn normalized_interaction(&self) -> f64 {
        self.coefficient
    }

    /// Assess physical `[sigma_11,sigma_22,tau_12]` in material axes.
    /// Unrepresentable inputs or arithmetic refuse instead of reporting zero.
    pub fn assess(
        &self,
        stress_pa: [f64; 3],
    ) -> Result<LaminaTsaiWuAssessment, MaterialStatePointError> {
        let normalized: [f64; 3] = core::array::from_fn(|i| stress_pa[i] / self.scales_pa[i]);
        if stress_pa.iter().any(|v| !v.is_finite())
            || normalized
                .iter()
                .enumerate()
                .any(|(i, v)| !v.is_finite() || (stress_pa[i] != 0.0 && *v == 0.0))
        {
            return Err(invalid(
                "Tsai-Wu requires finite representable material-axis stresses",
            ));
        }
        let scale = normalized.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Ok(LaminaTsaiWuAssessment {
                stress_pa,
                utilization: 0.0,
            });
        }
        let unit = normalized.map(|v| v / scale);
        if unit
            .iter()
            .zip(normalized)
            .any(|(u, v)| *u == 0.0 && v != 0.0)
        {
            return Err(invalid("Tsai-Wu stress dynamic range is unrepresentable"));
        }
        let [u, v, w] = unit;
        // Cholesky form of the positive quadratic term; no cancellation of
        // nearly equal squared normal stresses at an interaction near +/-1.
        let q = hypot(
            hypot(u + self.coefficient * v, self.transverse_factor * v),
            w,
        );
        let b = self.linear[0] * u + self.linear[1] * v;
        let radius = hypot(0.5 * b, q);
        let radial = if b >= 0.0 {
            radius + 0.5 * b
        } else {
            q * (q / (radius - 0.5 * b))
        };
        let utilization = scale * radial;
        if !b.is_finite() || !utilization.is_finite() || utilization <= 0.0 {
            return Err(invalid("Tsai-Wu utilization is unrepresentable"));
        }
        Ok(LaminaTsaiWuAssessment {
            stress_pa,
            utilization,
        })
    }
}

/// Coupled criterion at one original stress sample; no isolated mode is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminaTsaiWuAssessment {
    stress_pa: [f64; 3],
    utilization: f64,
}
impl LaminaTsaiWuAssessment {
    /// Actual physical material-axis stress used by the criterion [Pa].
    pub const fn material_stress_pa(&self) -> [f64; 3] {
        self.stress_pa
    }
    /// Positive radial failure index. `stress/utilization` is on the envelope
    /// for nonzero stress; this scales thermal and mechanical stress together.
    pub const fn utilization(&self) -> f64 {
        self.utilization
    }
    /// Equality or exceedance of the explicitly selected coupled envelope.
    pub fn reaches_limit(&self) -> bool {
        self.utilization >= 1.0
    }
}

/// Explicit interaction selection for one original layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TsaiWuInteractionSelection {
    /// Dimensionless normalized `F12` property name, with no default value.
    pub property: String,
    /// Initial claim selection; all subsequent refreshes pin the same claim.
    pub selection: MaterialPropertySelection,
}

#[derive(Clone, Debug)]
struct PreparedTsaiWuLamina {
    strength: LaminaTsaiWu,
    rotation: [f64; 9],
}

/// Coupled strength on the same original thermal/elastic layers and axes.
#[derive(Clone, Debug)]
pub struct PreparedLaminateTsaiWu {
    thermal: ResolvedLaminateThermal,
    layers: Vec<PreparedTsaiWuLamina>,
}
impl PreparedLaminateTsaiWu {
    /// Extend an admitted laminate by selecting one independent interaction
    /// source per layer. Existing strength sources and rotations are moved,
    /// without another strength query or a duplicate mechanics owner.
    pub fn from_maximum_stress(
        prepared: PreparedLaminateStrength,
        selections: &[TsaiWuInteractionSelection],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateStrengthError> {
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        if selections.len() != prepared.layers.len() {
            return Err(LaminateStrengthError::LayerCount {
                expected: prepared.layers.len(),
                actual: selections.len(),
            });
        }
        let PreparedLaminateStrength {
            thermal,
            layers: original,
        } = prepared;
        let mut layers = Vec::with_capacity(original.len());
        for (layer, (old, selected)) in original.into_iter().zip(selections).enumerate() {
            if cancelled() {
                return Err(LaminateStrengthError::Cancelled);
            }
            let strength = LaminaTsaiWu::resolve(
                &thermal.source().layers()[layer].input().card,
                old.strength,
                selected.property.clone(),
                selected.selection.clone(),
            )
            .map_err(|source| LaminateStrengthError::Material { layer, source })?;
            layers.push(PreparedTsaiWuLamina {
                strength,
                rotation: old.rotation,
            });
        }
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        Ok(Self { thermal, layers })
    }

    /// Refresh the original elasticity, expansion, strengths and interaction
    /// sources together. Failure leaves this immutable snapshot unchanged.
    pub fn at_temperatures(
        &self,
        temperatures_k: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateStrengthError> {
        let thermal = self
            .thermal
            .source()
            .at_temperatures(temperatures_k, &mut cancelled)
            .map_err(|e| match e {
                LaminateThermalSourceError::Cancelled => LaminateStrengthError::Cancelled,
                other => LaminateStrengthError::ThermalSource(other),
            })?;
        let mut layers = Vec::with_capacity(self.layers.len());
        for (layer, old) in self.layers.iter().enumerate() {
            if cancelled() {
                return Err(LaminateStrengthError::Cancelled);
            }
            let strength = old
                .strength
                .at_query_point(
                    &thermal.source().layers()[layer].input().card,
                    &thermal.layers()[layer].elasticity().query(),
                )
                .map_err(|source| LaminateStrengthError::Material { layer, source })?;
            layers.push(PreparedTsaiWuLamina {
                strength,
                rotation: old.rotation,
            });
        }
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        Ok(Self { thermal, layers })
    }

    /// Exact native mechanics, current material states and expansion paths.
    pub const fn thermal(&self) -> &ResolvedLaminateThermal {
        &self.thermal
    }
    /// Criterion and complete sources of one original layer, if in range.
    pub fn layer_strength(&self, layer: usize) -> Option<&LaminaTsaiWu> {
        self.layers.get(layer).map(|p| &p.strength)
    }
    /// Recover both actual layer faces, then rotate into the original material
    /// axes. Kinematics use `[epsilon_xx,epsilon_yy,gamma_xy]` in section axes.
    pub fn assess_layer(
        &self,
        layer: usize,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
    ) -> Result<LayerTsaiWuAssessment, LaminateStrengthError> {
        let stress = self
            .thermal
            .prepared()
            .layer_stress(layer, membrane, gradient_per_m)
            .map_err(LaminateStrengthError::Thermal)?;
        self.assess_stresses(layer, stress)
    }

    /// Assess original material-axis strengths and interaction using the actual
    /// free-sided beam face stresses from this exact prepared thermal snapshot.
    pub fn assess_beam_layer(
        &self,
        beam: &PreparedLaminateBeam,
        layer: usize,
        eps0: f64,
        kappa: f64,
    ) -> Result<LayerTsaiWuAssessment, LaminateStrengthError> {
        admit_beam(&self.thermal, beam)?;
        let stress = beam
            .layer_stress(layer, eps0, kappa)
            .map_err(LaminateStrengthError::Thermal)?;
        self.assess_stresses(layer, stress)
    }

    fn assess_stresses(
        &self,
        layer: usize,
        stress: LayerThermalStress,
    ) -> Result<LayerTsaiWuAssessment, LaminateStrengthError> {
        let prepared = &self.layers[layer];
        let assess = |s| {
            prepared
                .strength
                .assess(rotate(&prepared.rotation, s))
                .map_err(|source| LaminateStrengthError::Material { layer, source })
        };
        Ok(LayerTsaiWuAssessment {
            bottom: assess(stress.bottom_pa)?,
            top: assess(stress.top_pa)?,
        })
    }
    /// Governing face of all original layers, with no hot allocations. The
    /// earlier layer and then bottom face win exact ties. Every layer and the
    /// final publication boundary poll cancellation.
    pub fn assess(
        &self,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateTsaiWuAssessment, LaminateStrengthError> {
        self.assess_from(
            |layer| self.assess_layer(layer, membrane, gradient_per_m),
            cancelled,
        )
    }

    /// Governing beam layer face using the existing Tsai-Wu radial utilization.
    /// All source and temperature identity rules are the same as the original
    /// maximum-stress beam assessment. No failure-mode or degradation law is added.
    pub fn assess_beam(
        &self,
        beam: &PreparedLaminateBeam,
        eps0: f64,
        kappa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateTsaiWuAssessment, LaminateStrengthError> {
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        admit_beam(&self.thermal, beam)?;
        self.assess_from(
            |layer| {
                let stress = beam
                    .layer_stress(layer, eps0, kappa)
                    .map_err(LaminateStrengthError::Thermal)?;
                self.assess_stresses(layer, stress)
            },
            cancelled,
        )
    }

    fn assess_from(
        &self,
        mut recover: impl FnMut(usize) -> Result<LayerTsaiWuAssessment, LaminateStrengthError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateTsaiWuAssessment, LaminateStrengthError> {
        let mut worst: Option<LaminateTsaiWuAssessment> = None;
        for layer in 0..self.layers.len() {
            if cancelled() {
                return Err(LaminateStrengthError::Cancelled);
            }
            let faces = recover(layer)?;
            for (face, assessment) in [
                (LaminaFace::Bottom, faces.bottom),
                (LaminaFace::Top, faces.top),
            ] {
                if worst
                    .as_ref()
                    .is_none_or(|old| assessment.utilization() > old.assessment.utilization())
                {
                    worst = Some(LaminateTsaiWuAssessment {
                        layer,
                        face,
                        assessment,
                    });
                }
            }
        }
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        Ok(worst.expect("source-bound laminate has at least one layer"))
    }
}

/// Coupled indications at both faces of one original layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerTsaiWuAssessment {
    /// The layer's own bottom face, including its thermal stress.
    pub bottom: LaminaTsaiWuAssessment,
    /// The layer's own top face, including its thermal stress.
    pub top: LaminaTsaiWuAssessment,
}
/// Location and coupled indication of the most highly utilized layer face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateTsaiWuAssessment {
    /// Original zero-based layer index.
    pub layer: usize,
    /// Physical face of this layer.
    pub face: LaminaFace,
    /// Stress and radial utilization of the coupled failure envelope.
    pub assessment: LaminaTsaiWuAssessment,
}
