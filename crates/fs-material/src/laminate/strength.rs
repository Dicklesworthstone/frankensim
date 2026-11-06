//! Plane-stress maximum-stress assessment in each lamina's material axes.
//!
//! Five explicit positive strength magnitudes supply independent limits for
//! tension/compression along axes 1 and 2, and for in-plane shear. The index is
//! `max(s11+/Xt, (-s11)+/Xc, s22+/Yt, (-s22)+/Yc, abs(t12)/S)`.
//! No isotropy, interaction coefficient, safety factor or missing strength is
//! inferred. See the maximum-stress theory in the Abaqus documentation:
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-failuremeasures.htm>.
//!
//! The native consumer resolves strength from each thermal snapshot's original
//! material card at its complete current condition point. It recovers both
//! faces, including bending and thermal self-stress, before rotating physical
//! stress into the same axes as the layer's elastic constants. With affine
//! stress and constant strength within each layer, the convex maximum-stress
//! index attains its thickness maximum at a face. Averaged section resultants
//! cannot replace this check. Preparation caches sources and rotations;
//! assessment allocates nothing and performs no source queries or trigonometry.
//!
//! This is an infinitesimal-strain elastic failure indication, not progressive
//! damage or a prediction of ultimate laminate load. It supplies no interaction
//! surface, out-of-plane/interlaminar failure, delamination, buckling, fatigue,
//! fracture energy or post-failure stiffness. Exceedance never removes material.

use super::thermal::LaminateThermalError;
use super::thermal::source::{LaminateThermalSourceError, ResolvedLaminateThermal};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Pressure, QuantitySpec};

/// Coupled plane-stress failure with an explicitly sourced interaction term.
pub mod tsai_wu;

/// Explicit property roles in the original elastic material 1/2 frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LaminaStrengthProperties {
    /// Tensile strength magnitudes along axes `[1, 2]` [Pa].
    pub tensile: [String; 2],
    /// Compressive strength magnitudes along axes `[1, 2]` [Pa], both positive.
    pub compressive: [String; 2],
    /// In-plane shear strength magnitude [Pa], symmetric in shear sign.
    pub shear: String,
}

impl LaminaStrengthProperties {
    fn ordered(&self) -> [&str; 5] {
        [
            &self.tensile[0],
            &self.compressive[0],
            &self.tensile[1],
            &self.compressive[1],
            &self.shear,
        ]
    }
}

/// Independent failure modes, in the deterministic tie-breaking order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaminaFailureMode {
    /// Positive material-axis 1 normal stress.
    LongitudinalTension,
    /// Negative material-axis 1 normal stress.
    LongitudinalCompression,
    /// Positive material-axis 2 normal stress.
    TransverseTension,
    /// Negative material-axis 2 normal stress.
    TransverseCompression,
    /// Absolute physical material-axis 12 shear stress.
    InPlaneShear,
}

const MODES: [LaminaFailureMode; 5] = [
    LaminaFailureMode::LongitudinalTension,
    LaminaFailureMode::LongitudinalCompression,
    LaminaFailureMode::TransverseTension,
    LaminaFailureMode::TransverseCompression,
    LaminaFailureMode::InPlaneShear,
];

/// Original selected strength claims at one complete physical condition point.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminaMaximumStress {
    resolved: ResolvedMaterialStatePoint,
    properties: LaminaStrengthProperties,
    limits_pa: [f64; 5],
}

impl LaminaMaximumStress {
    /// Resolve all five roles together; repeated names explicitly reuse a claim.
    /// A pinned selection must cover the distinct requested properties exactly.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        properties: LaminaStrengthProperties,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        let mut requirements = Vec::with_capacity(5);
        for name in properties.ordered() {
            let requirement = ScalarPropertyRequirement::try_new(
                name,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?;
            if !requirements.contains(&requirement) {
                requirements.push(requirement);
            }
        }
        Self::try_from_resolved(
            resolve_material_state_point(card, point, &requirements, selection)?,
            properties,
        )
    }

    /// Select explicit orthotropic maximum-stress roles from an existing bundle.
    /// Each strength is a positive Pa magnitude at static external conditions;
    /// unowned strain, strain rate and frequency dependence refuse.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        properties: LaminaStrengthProperties,
    ) -> Result<Self, MaterialStatePointError> {
        let mut limits_pa = [0.0; 5];
        for (i, name) in properties.ordered().into_iter().enumerate() {
            let p = resolved.property(name).ok_or_else(|| {
                invalid("lamina strength requires every named tensile, compressive and shear claim")
            })?;
            if p.requirement().quantity() != QuantitySpec::dimensional(Pressure::DIMS)
                || !p.is_static_at_fixed_external_conditions()
                || !p.value_si().is_finite()
                || p.value_si() <= 0.0
            {
                return Err(invalid(
                    "lamina strengths must be positive finite Pa magnitudes at static fixed external conditions",
                ));
            }
            limits_pa[i] = p.value_si();
        }
        Ok(Self {
            resolved,
            properties,
            limits_pa,
        })
    }

    /// Refresh the same original card and exact selected claims together.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(
            self.resolved.at_query_point(card, point)?,
            self.properties.clone(),
        )
    }

    /// Original material, complete condition point and selected-source receipts.
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Explicit assignment of source properties to material-axis strength roles.
    pub const fn properties(&self) -> &LaminaStrengthProperties {
        &self.properties
    }

    /// Positive strength magnitudes `[Xt, Xc, Yt, Yc, S12]` [Pa].
    pub const fn limits_pa(&self) -> [f64; 5] {
        self.limits_pa
    }

    /// Assess physical `[sigma_11, sigma_22, tau_12]` in the original material
    /// axes, tension positive. The caller owns stress/condition compatibility
    /// when using this leaf directly; [`PreparedLaminateStrength`] binds it.
    pub fn assess(
        &self,
        stress_pa: [f64; 3],
    ) -> Result<LaminaStrengthAssessment, MaterialStatePointError> {
        if stress_pa.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "lamina strength assessment requires finite material-axis stress",
            ));
        }
        let [s1, s2, t12] = stress_pa;
        let demand = [
            s1.max(0.0),
            (-s1).max(0.0),
            s2.max(0.0),
            (-s2).max(0.0),
            t12.abs(),
        ];
        let utilizations = core::array::from_fn(|i| demand[i] / self.limits_pa[i]);
        if utilizations
            .iter()
            .enumerate()
            .any(|(i, u)| !u.is_finite() || (demand[i] > 0.0 && *u == 0.0))
        {
            return Err(invalid("lamina strength utilization is unrepresentable"));
        }
        let mut critical = 0;
        for i in 1..5 {
            if utilizations[i] > utilizations[critical] {
                critical = i;
            }
        }
        Ok(LaminaStrengthAssessment {
            stress_pa,
            utilizations,
            critical,
        })
    }
}

/// Local indication with all signed failure modes retained separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminaStrengthAssessment {
    stress_pa: [f64; 3],
    utilizations: [f64; 5],
    critical: usize,
}
impl LaminaStrengthAssessment {
    /// Actual material-axis `[sigma_11, sigma_22, tau_12]` [Pa].
    pub const fn material_stress_pa(&self) -> [f64; 3] {
        self.stress_pa
    }
    /// Ratios in `[1 tension, 1 compression, 2 tension, 2 compression, 12 shear]` order.
    pub const fn mode_utilizations(&self) -> [f64; 5] {
        self.utilizations
    }
    /// Maximum of the five independent demand/strength ratios.
    pub const fn utilization(&self) -> f64 {
        self.utilizations[self.critical]
    }
    /// Governing mode, or `None` at zero stress. Equal indices retain the first
    /// mode in [`LaminaFailureMode`]'s declaration order without a tolerance.
    pub fn critical_mode(&self) -> Option<LaminaFailureMode> {
        (self.utilization() > 0.0).then_some(MODES[self.critical])
    }
    /// Equality/exceedance indicates a limit crossing, not damage evolution.
    pub fn reaches_limit(&self) -> bool {
        self.utilization() >= 1.0
    }
}

/// One explicit strength selection for one original layer, in layer order.
/// Card, temperature and other conditions come from the mechanical snapshot.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LaminaStrengthSelection {
    /// Five source roles in this layer's original elastic material frame.
    pub properties: LaminaStrengthProperties,
    /// Initial selection; subsequent temperature refreshes pin these claims.
    pub selection: MaterialPropertySelection,
}

/// Refusal of a complete preparation or assessment, without partial publication.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminateStrengthError {
    /// Exactly one selection is required per original layer.
    LayerCount { expected: usize, actual: usize },
    /// Strength admission or evaluation failed for the named layer.
    Material {
        layer: usize,
        source: MaterialStatePointError,
    },
    /// Native stress recovery refused the layer or strain inputs.
    Thermal(LaminateThermalError),
    /// Original elastic/thermal sources refused the requested temperatures.
    ThermalSource(LaminateThermalSourceError),
    /// Cancellation occurred before the complete result was available.
    Cancelled,
}
impl core::fmt::Display for LaminateStrengthError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MATERIAL-LAMINATE-STRENGTH: {self:?}")
    }
}
impl std::error::Error for LaminateStrengthError {}

#[derive(Clone, Debug)]
struct PreparedLaminaStrength {
    strength: LaminaMaximumStress,
    // Stress rotation from section axes to the layer's material axes.
    rotation: [f64; 9],
}

/// Immutable strength and native elastic/thermal mechanics at the SAME state.
/// Assessing a new temperature requires a refreshed complete snapshot.
#[derive(Clone, Debug)]
pub struct PreparedLaminateStrength {
    thermal: ResolvedLaminateThermal,
    layers: Vec<PreparedLaminaStrength>,
}
impl PreparedLaminateStrength {
    /// Bind the ordered strength roles to the exact native material cards and
    /// current physical conditions, never to an independently supplied state.
    pub fn resolve(
        thermal: ResolvedLaminateThermal,
        selections: &[LaminaStrengthSelection],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateStrengthError> {
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        if selections.len() != thermal.layers().len() {
            return Err(LaminateStrengthError::LayerCount {
                expected: thermal.layers().len(),
                actual: selections.len(),
            });
        }
        let mut layers = Vec::with_capacity(selections.len());
        for (layer, selection) in selections.iter().enumerate() {
            if cancelled() {
                return Err(LaminateStrengthError::Cancelled);
            }
            let input = thermal.source().layers()[layer].input();
            let strength = LaminaMaximumStress::resolve(
                &input.card,
                &thermal.layers()[layer].elasticity().query(),
                selection.properties.clone(),
                selection.selection.clone(),
            )
            .map_err(|source| LaminateStrengthError::Material { layer, source })?;
            let (s, c) = fs_math::det::sin_cos(input.orientation_rad);
            let (cc, ss, cs) = (c * c, s * s, c * s);
            layers.push(PreparedLaminaStrength {
                strength,
                rotation: [cc, ss, 2.0 * cs, ss, cc, -2.0 * cs, -cs, cs, cc - ss],
            });
        }
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        Ok(Self { thermal, layers })
    }

    /// Refresh pinned strengths, stiffness and integrated expansion together.
    /// Use this returned snapshot's [`Self::thermal`] for native assembly;
    /// failure leaves the old snapshot and its original receipts intact.
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
            layers.push(PreparedLaminaStrength {
                strength,
                rotation: old.rotation,
            });
        }
        if cancelled() {
            return Err(LaminateStrengthError::Cancelled);
        }
        Ok(Self { thermal, layers })
    }

    /// Exact native mechanics, current elasticity and complete expansion paths.
    pub const fn thermal(&self) -> &ResolvedLaminateThermal {
        &self.thermal
    }

    /// Strength source state for one original layer, or `None` for an invalid index.
    pub fn layer_strength(&self, layer: usize) -> Option<&LaminaMaximumStress> {
        self.layers.get(layer).map(|p| &p.strength)
    }

    /// Recover and assess both physical faces, including thermal residual stress.
    /// `membrane` and `gradient_per_m` use `[epsilon_xx, epsilon_yy, gamma_xy]`
    /// in section axes; the latter multiplies the native signed thickness z.
    pub fn assess_layer(
        &self,
        layer: usize,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
    ) -> Result<LayerStrengthAssessment, LaminateStrengthError> {
        let stress = self
            .thermal
            .prepared()
            .layer_stress(layer, membrane, gradient_per_m)
            .map_err(LaminateStrengthError::Thermal)?;
        let prepared = &self.layers[layer];
        let assess = |s| {
            prepared
                .strength
                .assess(rotate(&prepared.rotation, s))
                .map_err(|source| LaminateStrengthError::Material { layer, source })
        };
        Ok(LayerStrengthAssessment {
            bottom: assess(stress.bottom_pa)?,
            top: assess(stress.top_pa)?,
        })
    }

    /// Most critical face across every original layer, without allocation.
    /// Equal indices retain the earlier layer, then its bottom face. Cancellation
    /// is polled per layer and before publication. Per-face details remain
    /// available through [`Self::assess_layer`].
    pub fn assess(
        &self,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateStrengthAssessment, LaminateStrengthError> {
        let mut worst: Option<LaminateStrengthAssessment> = None;
        for layer in 0..self.layers.len() {
            if cancelled() {
                return Err(LaminateStrengthError::Cancelled);
            }
            let faces = self.assess_layer(layer, membrane, gradient_per_m)?;
            for (face, assessment) in [
                (LaminaFace::Bottom, faces.bottom),
                (LaminaFace::Top, faces.top),
            ] {
                if worst
                    .as_ref()
                    .is_none_or(|old| assessment.utilization() > old.assessment.utilization())
                {
                    worst = Some(LaminateStrengthAssessment {
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

/// Face belonging to one original layer, not a shared interlaminar traction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaminaFace {
    /// Face with smaller signed thickness coordinate.
    Bottom,
    /// Face with larger signed thickness coordinate.
    Top,
}

/// Strength indications at both faces of one original layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerStrengthAssessment {
    /// The layer's own bottom-face stress and indication.
    pub bottom: LaminaStrengthAssessment,
    /// The layer's own top-face stress and indication.
    pub top: LaminaStrengthAssessment,
}

/// Location and material-axis indication of the governing layer face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateStrengthAssessment {
    /// Original zero-based layer index.
    pub layer: usize,
    /// Physical face of this layer.
    pub face: LaminaFace,
    /// The stress, all independent mode ratios and the governing mode.
    pub assessment: LaminaStrengthAssessment,
}

fn rotate(q: &[f64; 9], stress: [f64; 3]) -> [f64; 3] {
    let rotated = core::array::from_fn(|i| (0..3).map(|j| q[3 * i + j] * stress[j]).sum::<f64>());
    if rotated.iter().all(|s: &f64| s.is_finite()) {
        return rotated;
    }
    // Only scale when intermediate addition overflowed. Ordinary component
    // arithmetic preserves small independent stresses in widely separated axes.
    let scale = stress.iter().map(|x| x.abs()).fold(0.0, f64::max);
    if scale == 0.0 {
        return stress;
    }
    let normalized = stress.map(|s| s / scale);
    core::array::from_fn(|i| (0..3).map(|j| q[3 * i + j] * normalized[j]).sum::<f64>() * scale)
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
