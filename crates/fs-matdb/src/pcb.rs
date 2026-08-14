//! Evidence-carrying PCB laminate homogenization (bead
//! `frankensim-extreal-program-f85xj.5.6`).
//!
//! A printed-circuit board is a stack of copper/dielectric layers, not an
//! isotropic material name. This module binds every constituent conductivity
//! to an immutable [`MaterialCard`] and [`PropertyUsageReceipt`], keeps copper
//! coverage as a bounded provenance-bearing input, and evaluates the standard
//! laminate first rung:
//!
//! - thickness-weighted parallel conduction in the board plane;
//! - a series rule through the layer stack;
//! - nominal Reuss and Voigt structural bounds;
//! - a shared-source response record for every uncertain coverage fraction.
//!
//! The result is DATA. It does not add vias, spreading resistance, trace-scale
//! hot spots, delamination, or a statistical distribution that the source did
//! not supply.

use std::collections::BTreeSet;
use std::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_qty::Dims;

use crate::{
    MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyUsageReceipt, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};

/// Thermal-conductivity dimensions, W/(m K), in `[m, kg, s, K, A, mol]`.
pub const PCB_THERMAL_CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);

/// Closed schema version of [`PcbHomogenizedConductivity`].
pub const PCB_HOMOGENIZATION_SCHEMA_VERSION: u32 = 1;

/// Content-identity domain of [`PcbHomogenizedConductivity`].
pub const PCB_HOMOGENIZATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-matdb.pcb-homogenization.v1";

/// Typed refusal surface for PCB stackup admission and homogenization.
#[derive(Debug, Clone, PartialEq)]
pub enum PcbHomogenizationError {
    /// A structural field is blank, non-finite, non-positive, out of range, or
    /// otherwise inadmissible.
    InvalidField {
        /// Stable field name.
        field: &'static str,
        /// Actionable diagnosis.
        detail: String,
    },
    /// An upstream material-card query refused.
    MaterialQuery(MatDbError),
    /// The declared feature/thickness scale separation is outside its admitted
    /// homogenization domain.
    ScaleSeparation {
        /// Observed feature-size / board-thickness ratio.
        observed_ratio: f64,
        /// Largest ratio the caller declared admissible.
        maximum_ratio: f64,
    },
}

impl fmt::Display for PcbHomogenizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField { field, detail } => {
                write!(f, "PCB homogenization field '{field}' is invalid: {detail}")
            }
            Self::MaterialQuery(error) => write!(f, "PCB material query refused: {error}"),
            Self::ScaleSeparation {
                observed_ratio,
                maximum_ratio,
            } => write!(
                f,
                "PCB feature/thickness ratio {observed_ratio} exceeds the declared \
                 homogenization limit {maximum_ratio}; resolve the copper geometry or widen the \
                 model only with new evidence"
            ),
        }
    }
}

impl std::error::Error for PcbHomogenizationError {}

impl From<MatDbError> for PcbHomogenizationError {
    fn from(value: MatDbError) -> Self {
        Self::MaterialQuery(value)
    }
}

/// One positive scalar conductivity selected from an immutable material card.
///
/// The nominal value and any source-stated band are copied from the query
/// answer, while the exact usage receipt and material-card identity remain
/// attached. `uncertainty_complete == false` means the source used
/// [`UncertaintyModel::Unstated`]; the nominal value remains usable as an
/// estimate, but no complete material-property band may be inferred.
#[derive(Debug, Clone, PartialEq)]
pub struct PcbConductivityDatum {
    /// Named manufactured material state.
    material_state: MaterialStateId,
    /// Exact immutable card identity.
    material_card: ContentHash,
    /// Exact property-use receipt.
    receipt: PropertyUsageReceipt,
    /// Selected nominal conductivity, W/(m K).
    nominal_w_mk: f64,
    /// Lower source-stated conductivity, or the nominal value when uncertainty
    /// was unstated.
    lower_w_mk: f64,
    /// Upper source-stated conductivity, or the nominal value when uncertainty
    /// was unstated.
    upper_w_mk: f64,
    /// Whether `lower_w_mk..=upper_w_mk` includes a stated uncertainty model.
    uncertainty_complete: bool,
}

impl PcbConductivityDatum {
    /// Query a material card and retain the complete answer as one PCB
    /// constituent datum.
    ///
    /// # Errors
    ///
    /// Returns the upstream material-query refusal, a dimensions refusal, or a
    /// non-positive/overflowed conductivity-band refusal.
    pub fn from_card(
        card: &MaterialCard,
        property: &str,
        point: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<Self, PcbHomogenizationError> {
        let answer = card.claims().query(property, point, policy)?;
        Self::from_answer(card, answer)
    }

    fn from_answer(
        card: &MaterialCard,
        answer: MaterialAnswer,
    ) -> Result<Self, PcbHomogenizationError> {
        let sample = &answer.evidence.value;
        if sample.dims != PCB_THERMAL_CONDUCTIVITY_DIMS {
            return Err(PcbHomogenizationError::InvalidField {
                field: "conductivity-dimensions",
                detail: format!(
                    "expected {:?} for W/(m K), found {:?}",
                    PCB_THERMAL_CONDUCTIVITY_DIMS, sample.dims
                ),
            });
        }
        if !(sample.value.is_finite() && sample.value > 0.0) {
            return Err(PcbHomogenizationError::InvalidField {
                field: "conductivity",
                detail: format!(
                    "material {} answered {}; a constituent conductivity must be finite and \
                     positive",
                    card.id(),
                    sample.value
                ),
            });
        }

        let (lower, upper, complete) = match sample.uncertainty {
            UncertaintyModel::Unstated => (sample.value, sample.value, false),
            UncertaintyModel::HalfWidth { half_width, .. } => {
                (sample.value - half_width, sample.value + half_width, true)
            }
            UncertaintyModel::RelativeHalfWidth { fraction, .. } => (
                sample.value * (1.0 - fraction),
                sample.value * (1.0 + fraction),
                true,
            ),
        };
        if !(lower.is_finite() && upper.is_finite() && lower > 0.0 && lower <= upper) {
            return Err(PcbHomogenizationError::InvalidField {
                field: "conductivity-uncertainty",
                detail: format!(
                    "the source uncertainty around {} W/(m K) produces inadmissible positive \
                     bounds [{lower}, {upper}]",
                    sample.value
                ),
            });
        }

        Ok(Self {
            material_state: card.id().clone(),
            material_card: card.content_hash(),
            receipt: answer.receipt,
            nominal_w_mk: sample.value,
            lower_w_mk: lower,
            upper_w_mk: upper,
            uncertainty_complete: complete,
        })
    }

    /// Named manufactured material state.
    #[must_use]
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }

    /// Exact immutable material-card identity.
    #[must_use]
    pub fn material_card(&self) -> ContentHash {
        self.material_card
    }

    /// Exact property-use receipt retained by this datum.
    #[must_use]
    pub fn receipt(&self) -> &PropertyUsageReceipt {
        &self.receipt
    }

    /// Selected nominal conductivity in W/(m K).
    #[must_use]
    pub fn nominal_w_mk(&self) -> f64 {
        self.nominal_w_mk
    }

    /// Lower source-stated conductivity in W/(m K).
    #[must_use]
    pub fn lower_w_mk(&self) -> f64 {
        self.lower_w_mk
    }

    /// Upper source-stated conductivity in W/(m K).
    #[must_use]
    pub fn upper_w_mk(&self) -> f64 {
        self.upper_w_mk
    }

    /// Whether the lower/upper interval comes from a stated uncertainty model.
    #[must_use]
    pub fn uncertainty_complete(&self) -> bool {
        self.uncertainty_complete
    }
}

/// One bounded copper-area fraction and its load-bearing source.
///
/// Bounds are epistemic/declared bounds, not a probability distribution.
#[derive(Debug, Clone, PartialEq)]
pub struct CopperCoverage {
    /// Stable source identity used by the shared-direction influence record.
    source_id: String,
    /// Nominal copper area fraction.
    nominal: f64,
    /// Lower admitted area fraction.
    lower: f64,
    /// Upper admitted area fraction.
    upper: f64,
    /// Design-data or declared-estimate provenance.
    provenance: Provenance,
}

impl CopperCoverage {
    /// Admit a bounded copper fraction.
    ///
    /// # Errors
    ///
    /// Refuses blank source ids, non-finite or unordered fractions, fractions
    /// outside `[0, 1]`, and incomplete provenance.
    pub fn new(
        source_id: impl Into<String>,
        nominal: f64,
        lower: f64,
        upper: f64,
        provenance: Provenance,
    ) -> Result<Self, PcbHomogenizationError> {
        let source_id = source_id.into();
        if source_id.trim().is_empty() {
            return Err(PcbHomogenizationError::InvalidField {
                field: "coverage-source-id",
                detail: "the coverage source must have a stable nonblank identity".to_string(),
            });
        }
        if !(lower.is_finite()
            && nominal.is_finite()
            && upper.is_finite()
            && 0.0 <= lower
            && lower <= nominal
            && nominal <= upper
            && upper <= 1.0)
        {
            return Err(PcbHomogenizationError::InvalidField {
                field: "copper-coverage",
                detail: format!(
                    "expected 0 <= lower <= nominal <= upper <= 1, got \
                     [{lower}, {nominal}, {upper}]"
                ),
            });
        }
        provenance.validate()?;
        Ok(Self {
            source_id,
            nominal,
            lower,
            upper,
            provenance,
        })
    }

    /// Stable identity of the source that supplied this coverage bound.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Nominal copper area fraction.
    #[must_use]
    pub fn nominal(&self) -> f64 {
        self.nominal
    }

    /// Lower admitted copper area fraction.
    #[must_use]
    pub fn lower(&self) -> f64 {
        self.lower
    }

    /// Upper admitted copper area fraction.
    #[must_use]
    pub fn upper(&self) -> f64 {
        self.upper
    }

    /// Design-data or declared-estimate provenance.
    #[must_use]
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }
}

/// One physical layer in a PCB stackup.
#[derive(Debug, Clone, PartialEq)]
pub struct PcbLayer {
    /// Stable layer name; stack order is identity-bearing.
    name: String,
    /// Positive layer thickness in metres.
    thickness_m: f64,
    /// Copper-state conductivity selected from a material card.
    copper: PcbConductivityDatum,
    /// Matrix/dielectric conductivity selected from a material card.
    matrix: PcbConductivityDatum,
    /// Copper area fraction and its bounded source.
    coverage: CopperCoverage,
}

impl PcbLayer {
    /// Admit one stackup layer.
    ///
    /// # Errors
    ///
    /// Refuses a blank layer name or a non-positive/non-finite thickness.
    pub fn new(
        name: impl Into<String>,
        thickness_m: f64,
        copper: PcbConductivityDatum,
        matrix: PcbConductivityDatum,
        coverage: CopperCoverage,
    ) -> Result<Self, PcbHomogenizationError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(PcbHomogenizationError::InvalidField {
                field: "layer-name",
                detail: "every physical layer needs a stable name".to_string(),
            });
        }
        if !(thickness_m.is_finite() && thickness_m > 0.0) {
            return Err(PcbHomogenizationError::InvalidField {
                field: "layer-thickness",
                detail: format!("layer '{name}' thickness {thickness_m} m is not positive"),
            });
        }
        Ok(Self {
            name,
            thickness_m,
            copper,
            matrix,
            coverage,
        })
    }

    /// Stable layer name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Physical layer thickness in metres.
    #[must_use]
    pub fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Copper-state conductivity datum.
    #[must_use]
    pub fn copper(&self) -> &PcbConductivityDatum {
        &self.copper
    }

    /// Matrix/dielectric conductivity datum.
    #[must_use]
    pub fn matrix(&self) -> &PcbConductivityDatum {
        &self.matrix
    }

    /// Copper area fraction and its bounded source.
    #[must_use]
    pub fn coverage(&self) -> &CopperCoverage {
        &self.coverage
    }
}

/// Right-handed orthonormal principal frame for the laminate.
///
/// Rows are the in-plane x axis, in-plane y axis, and board-normal axis,
/// expressed in the consumer's global coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcbPrincipalFrame {
    /// Principal axes as right-handed orthonormal rows.
    axes: [[f64; 3]; 3],
}

impl PcbPrincipalFrame {
    /// Admit a right-handed orthonormal frame.
    ///
    /// # Errors
    ///
    /// Refuses non-finite, non-unit, non-orthogonal, or reflected frames.
    pub fn new(axes: [[f64; 3]; 3]) -> Result<Self, PcbHomogenizationError> {
        const TOLERANCE: f64 = 1.0e-12;
        for (axis_index, axis) in axes.iter().enumerate() {
            if axis.iter().any(|value| !value.is_finite()) {
                return Err(PcbHomogenizationError::InvalidField {
                    field: "principal-frame",
                    detail: format!("axis {axis_index} contains a non-finite component"),
                });
            }
            let norm2 = axis.iter().map(|value| value * value).sum::<f64>();
            if (norm2 - 1.0).abs() > TOLERANCE {
                return Err(PcbHomogenizationError::InvalidField {
                    field: "principal-frame",
                    detail: format!("axis {axis_index} squared norm is {norm2}, not one"),
                });
            }
        }
        for left in 0..3 {
            for right in (left + 1)..3 {
                let dot = axes[left]
                    .iter()
                    .zip(axes[right])
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                if dot.abs() > TOLERANCE {
                    return Err(PcbHomogenizationError::InvalidField {
                        field: "principal-frame",
                        detail: format!("axes {left} and {right} have dot product {dot}"),
                    });
                }
            }
        }
        let determinant = axes[0][0].mul_add(
            axes[1][1] * axes[2][2] - axes[1][2] * axes[2][1],
            -axes[0][1] * (axes[1][0] * axes[2][2] - axes[1][2] * axes[2][0]),
        ) + axes[0][2] * (axes[1][0] * axes[2][1] - axes[1][1] * axes[2][0]);
        if (determinant - 1.0).abs() > TOLERANCE {
            return Err(PcbHomogenizationError::InvalidField {
                field: "principal-frame",
                detail: format!(
                    "frame determinant is {determinant}; the stackup frame must be right-handed"
                ),
            });
        }
        Ok(Self { axes })
    }

    /// Principal axes as right-handed orthonormal rows.
    #[must_use]
    pub fn axes(&self) -> [[f64; 3]; 3] {
        self.axes
    }
}

impl Default for PcbPrincipalFrame {
    fn default() -> Self {
        Self {
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }
}

/// Declared scale-separation admission rule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcbScaleSeparation {
    /// Representative unresolved copper feature size in metres.
    feature_size_m: f64,
    /// Largest admitted `feature_size / total_board_thickness` ratio.
    maximum_feature_to_thickness_ratio: f64,
}

impl PcbScaleSeparation {
    /// Construct a declared feature/thickness admission rule.
    ///
    /// # Errors
    ///
    /// Refuses non-positive/non-finite feature sizes or ratio limits outside
    /// `(0, 1]`.
    pub fn new(
        feature_size_m: f64,
        maximum_feature_to_thickness_ratio: f64,
    ) -> Result<Self, PcbHomogenizationError> {
        if !(feature_size_m.is_finite() && feature_size_m > 0.0) {
            return Err(PcbHomogenizationError::InvalidField {
                field: "feature-size",
                detail: format!("{feature_size_m} m is not finite and positive"),
            });
        }
        if !(maximum_feature_to_thickness_ratio.is_finite()
            && maximum_feature_to_thickness_ratio > 0.0
            && maximum_feature_to_thickness_ratio <= 1.0)
        {
            return Err(PcbHomogenizationError::InvalidField {
                field: "scale-separation-limit",
                detail: format!(
                    "{maximum_feature_to_thickness_ratio} is outside the admitted interval (0, 1]"
                ),
            });
        }
        Ok(Self {
            feature_size_m,
            maximum_feature_to_thickness_ratio,
        })
    }

    /// Representative unresolved copper feature size in metres.
    #[must_use]
    pub fn feature_size_m(&self) -> f64 {
        self.feature_size_m
    }

    /// Largest admitted feature-size to board-thickness ratio.
    #[must_use]
    pub fn maximum_feature_to_thickness_ratio(&self) -> f64 {
        self.maximum_feature_to_thickness_ratio
    }
}

/// Structural Reuss/Voigt bounds computed from nominal constituent values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcbStructuralBounds {
    /// Reuss (harmonic-mixture) lower structural bound, W/(m K).
    pub reuss_w_mk: f64,
    /// Voigt (arithmetic-mixture) upper structural bound, W/(m K).
    pub voigt_w_mk: f64,
}

/// Nominal and propagated principal conductivities.
///
/// The propagated band includes bounded coverage and every source-stated
/// constituent band. When `material_uncertainty_complete` on the parent result
/// is false, it does not include the missing material uncertainty.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcbPrincipalConductivity {
    /// Nominal `[in-plane x, in-plane y, through-plane]`, W/(m K).
    pub nominal_w_mk: [f64; 3],
    /// Conservative propagated lower values in the principal frame.
    pub lower_w_mk: [f64; 3],
    /// Conservative propagated upper values in the principal frame.
    pub upper_w_mk: [f64; 3],
}

/// One coverage source's coupled response across all principal directions.
///
/// This record preserves dependency: x, y, and z move under the SAME bounded
/// coverage source. It is not a covariance matrix and does not invent a
/// probability law.
#[derive(Debug, Clone, PartialEq)]
pub struct PcbCoverageInfluence {
    /// Layer carrying the coverage source.
    pub layer: String,
    /// Stable source identity.
    pub source_id: String,
    /// Bounded coverage `[lower, nominal, upper]`.
    pub coverage: [f64; 3],
    /// Principal conductivities when only this source is set to its lower
    /// bound and every other input stays nominal.
    pub principal_at_lower_w_mk: [f64; 3],
    /// Principal conductivities when only this source is set to its upper
    /// bound and every other input stays nominal.
    pub principal_at_upper_w_mk: [f64; 3],
    /// Source provenance copied from stackup admission.
    pub provenance: Provenance,
}

/// Explicit status of via-array correction in this first rung.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcbViaCorrection {
    /// No via/thermal-via correction is modeled.
    NotModeled,
}

/// Homogenized, evidence-carrying PCB conductivity.
#[derive(Debug, Clone, PartialEq)]
pub struct PcbHomogenizedConductivity {
    /// Closed output schema version.
    schema_version: u32,
    /// Stable stackup name.
    stackup_id: String,
    /// Right-handed principal frame.
    frame: PcbPrincipalFrame,
    /// Nominal and propagated principal values.
    principal: PcbPrincipalConductivity,
    /// Nominal global symmetric conductivity tensor, W/(m K).
    tensor_w_mk: [[f64; 3]; 3],
    /// Nominal Reuss/Voigt structural bracket.
    structural_bounds: PcbStructuralBounds,
    /// Observed feature-size / board-thickness ratio.
    feature_to_thickness_ratio: f64,
    /// Per-source shared-direction coverage responses.
    coverage_influences: Vec<PcbCoverageInfluence>,
    /// Every material use in stack order, copper then matrix.
    material_uses: Vec<PcbConductivityDatum>,
    /// True only when every material use carried a stated uncertainty model.
    material_uncertainty_complete: bool,
    /// Via-array correction status.
    via_correction: PcbViaCorrection,
    /// Content identity over the admitted stackup, algorithm version, and
    /// result.
    identity: ContentHash,
}

impl PcbHomogenizedConductivity {
    /// Closed output schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Stable stackup name.
    #[must_use]
    pub fn stackup_id(&self) -> &str {
        &self.stackup_id
    }

    /// Right-handed principal frame.
    #[must_use]
    pub const fn frame(&self) -> PcbPrincipalFrame {
        self.frame
    }

    /// Nominal and propagated principal values.
    #[must_use]
    pub const fn principal(&self) -> PcbPrincipalConductivity {
        self.principal
    }

    /// Nominal global symmetric conductivity tensor, W/(m K).
    #[must_use]
    pub const fn tensor_w_mk(&self) -> [[f64; 3]; 3] {
        self.tensor_w_mk
    }

    /// Nominal Reuss/Voigt structural bracket.
    #[must_use]
    pub const fn structural_bounds(&self) -> PcbStructuralBounds {
        self.structural_bounds
    }

    /// Observed feature-size / board-thickness ratio.
    #[must_use]
    pub const fn feature_to_thickness_ratio(&self) -> f64 {
        self.feature_to_thickness_ratio
    }

    /// Per-source shared-direction coverage responses.
    #[must_use]
    pub fn coverage_influences(&self) -> &[PcbCoverageInfluence] {
        &self.coverage_influences
    }

    /// Every material use in stack order, copper then matrix.
    #[must_use]
    pub fn material_uses(&self) -> &[PcbConductivityDatum] {
        &self.material_uses
    }

    /// True only when every material use carried a stated uncertainty model.
    #[must_use]
    pub const fn material_uncertainty_complete(&self) -> bool {
        self.material_uncertainty_complete
    }

    /// Via-array correction status.
    #[must_use]
    pub const fn via_correction(&self) -> PcbViaCorrection {
        self.via_correction
    }

    /// Content identity over the admitted stackup, algorithm version, and
    /// result.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Immutable PCB stackup input.
#[derive(Debug, Clone, PartialEq)]
pub struct PcbStackup {
    id: String,
    layers: Vec<PcbLayer>,
    frame: PcbPrincipalFrame,
    scale_separation: PcbScaleSeparation,
}

impl PcbStackup {
    /// Admit a named stackup.
    ///
    /// Layer order is identity-bearing even though this first effective
    /// through-plane series rule is permutation-invariant.
    ///
    /// # Errors
    ///
    /// Refuses a blank id, an empty stackup, duplicate layer names, or
    /// duplicate coverage source ids (cross-layer coverage correlation is not
    /// represented in v1).
    pub fn new(
        id: impl Into<String>,
        layers: Vec<PcbLayer>,
        frame: PcbPrincipalFrame,
        scale_separation: PcbScaleSeparation,
    ) -> Result<Self, PcbHomogenizationError> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(PcbHomogenizationError::InvalidField {
                field: "stackup-id",
                detail: "the stackup needs a stable nonblank identity".to_string(),
            });
        }
        if layers.is_empty() {
            return Err(PcbHomogenizationError::InvalidField {
                field: "layers",
                detail: "a PCB stackup must contain at least one physical layer".to_string(),
            });
        }
        let mut names = BTreeSet::new();
        let mut sources = BTreeSet::new();
        for layer in &layers {
            if !names.insert(layer.name.clone()) {
                return Err(PcbHomogenizationError::InvalidField {
                    field: "layer-name",
                    detail: format!("layer name '{}' is duplicated", layer.name),
                });
            }
            if !sources.insert(layer.coverage.source_id.clone()) {
                return Err(PcbHomogenizationError::InvalidField {
                    field: "coverage-source-id",
                    detail: format!(
                        "coverage source '{}' is reused across layers; v1 requires one explicit \
                         bounded source per layer rather than implying cross-layer correlation",
                        layer.coverage.source_id
                    ),
                });
            }
        }
        Ok(Self {
            id,
            layers,
            frame,
            scale_separation,
        })
    }

    /// Physical layers in stack order.
    #[must_use]
    pub fn layers(&self) -> &[PcbLayer] {
        &self.layers
    }

    /// Evaluate the evidence-carrying laminate first rung.
    ///
    /// # Errors
    ///
    /// Refuses a total-thickness overflow or a stackup outside its declared
    /// feature/thickness validity domain.
    pub fn homogenize(&self) -> Result<PcbHomogenizedConductivity, PcbHomogenizationError> {
        let total_thickness = self
            .layers
            .iter()
            .try_fold(0.0f64, |sum, layer| {
                let next = sum + layer.thickness_m;
                next.is_finite().then_some(next)
            })
            .ok_or_else(|| PcbHomogenizationError::InvalidField {
                field: "total-thickness",
                detail: "summing layer thicknesses overflowed".to_string(),
            })?;
        let separation_ratio = self.scale_separation.feature_size_m / total_thickness;
        if separation_ratio > self.scale_separation.maximum_feature_to_thickness_ratio {
            return Err(PcbHomogenizationError::ScaleSeparation {
                observed_ratio: separation_ratio,
                maximum_ratio: self.scale_separation.maximum_feature_to_thickness_ratio,
            });
        }

        let (principal, structural_bounds) =
            principal_and_structural(&self.layers, total_thickness);
        let tensor = rotate_principal(self.frame.axes, principal.nominal_w_mk);
        let coverage_influences = self
            .layers
            .iter()
            .enumerate()
            .map(|(index, layer)| PcbCoverageInfluence {
                layer: layer.name.clone(),
                source_id: layer.coverage.source_id.clone(),
                coverage: [
                    layer.coverage.lower,
                    layer.coverage.nominal,
                    layer.coverage.upper,
                ],
                principal_at_lower_w_mk: effective_principal(
                    &self.layers,
                    Some((index, layer.coverage.lower)),
                ),
                principal_at_upper_w_mk: effective_principal(
                    &self.layers,
                    Some((index, layer.coverage.upper)),
                ),
                provenance: layer.coverage.provenance.clone(),
            })
            .collect::<Vec<_>>();
        let material_uses = self
            .layers
            .iter()
            .flat_map(|layer| [layer.copper.clone(), layer.matrix.clone()])
            .collect::<Vec<_>>();
        let material_uncertainty_complete =
            material_uses.iter().all(|datum| datum.uncertainty_complete);
        let identity = homogenization_identity(
            self,
            &principal,
            &structural_bounds,
            separation_ratio,
            &tensor,
        );

        Ok(PcbHomogenizedConductivity {
            schema_version: PCB_HOMOGENIZATION_SCHEMA_VERSION,
            stackup_id: self.id.clone(),
            frame: self.frame,
            principal,
            tensor_w_mk: tensor,
            structural_bounds,
            feature_to_thickness_ratio: separation_ratio,
            coverage_influences,
            material_uses,
            material_uncertainty_complete,
            via_correction: PcbViaCorrection::NotModeled,
            identity,
        })
    }
}

/// Principal interval bounds (series/parallel per axis) and the
/// Reuss/Voigt structural bracket (stage of [`PcbStackup::homogenize`];
/// summation order preserved exactly).
fn principal_and_structural(
    layers: &[PcbLayer],
    total_thickness: f64,
) -> (PcbPrincipalConductivity, PcbStructuralBounds) {
    let nominal = effective_principal(layers, None);
    let mut in_plane_lower_numerator = 0.0f64;
    let mut in_plane_upper_numerator = 0.0f64;
    let mut through_lower_denominator = 0.0f64;
    let mut through_upper_denominator = 0.0f64;
    for layer in layers {
        let (lower, upper) = layer_conductivity_bounds(layer);
        in_plane_lower_numerator = lower.mul_add(layer.thickness_m, in_plane_lower_numerator);
        in_plane_upper_numerator = upper.mul_add(layer.thickness_m, in_plane_upper_numerator);
        through_lower_denominator += layer.thickness_m / lower;
        through_upper_denominator += layer.thickness_m / upper;
    }
    let in_plane_lower = in_plane_lower_numerator / total_thickness;
    let in_plane_upper = in_plane_upper_numerator / total_thickness;
    let through_lower = total_thickness / through_lower_denominator;
    let through_upper = total_thickness / through_upper_denominator;

    let mut reuss_denominator = 0.0f64;
    for layer in layers {
        let weight = layer.thickness_m / total_thickness;
        let coverage = layer.coverage.nominal;
        reuss_denominator += weight
            * (coverage / layer.copper.nominal_w_mk + (1.0 - coverage) / layer.matrix.nominal_w_mk);
    }
    let structural_bounds = PcbStructuralBounds {
        reuss_w_mk: 1.0 / reuss_denominator,
        voigt_w_mk: nominal[0],
    };
    let principal = PcbPrincipalConductivity {
        nominal_w_mk: nominal,
        lower_w_mk: [in_plane_lower, in_plane_lower, through_lower],
        upper_w_mk: [in_plane_upper, in_plane_upper, through_upper],
    };
    (principal, structural_bounds)
}

fn mix(coverage: f64, copper: f64, matrix: f64) -> f64 {
    coverage.mul_add(copper, (1.0 - coverage) * matrix)
}

fn effective_principal(layers: &[PcbLayer], coverage_override: Option<(usize, f64)>) -> [f64; 3] {
    let total_thickness = layers.iter().map(|layer| layer.thickness_m).sum::<f64>();
    let mut in_plane_numerator = 0.0f64;
    let mut through_denominator = 0.0f64;
    for (index, layer) in layers.iter().enumerate() {
        let coverage = coverage_override
            .filter(|(target, _)| *target == index)
            .map_or(layer.coverage.nominal, |(_, value)| value);
        let layer_k = mix(
            coverage,
            layer.copper.nominal_w_mk,
            layer.matrix.nominal_w_mk,
        );
        in_plane_numerator = layer_k.mul_add(layer.thickness_m, in_plane_numerator);
        through_denominator += layer.thickness_m / layer_k;
    }
    let in_plane = in_plane_numerator / total_thickness;
    let through = total_thickness / through_denominator;
    [in_plane, in_plane, through]
}

fn layer_conductivity_bounds(layer: &PcbLayer) -> (f64, f64) {
    let mut lower = f64::INFINITY;
    let mut upper = f64::NEG_INFINITY;
    for coverage in [layer.coverage.lower, layer.coverage.upper] {
        for copper in [layer.copper.lower_w_mk, layer.copper.upper_w_mk] {
            for matrix in [layer.matrix.lower_w_mk, layer.matrix.upper_w_mk] {
                let value = mix(coverage, copper, matrix);
                lower = lower.min(value);
                upper = upper.max(value);
            }
        }
    }
    (lower, upper)
}

fn rotate_principal(axes: [[f64; 3]; 3], principal: [f64; 3]) -> [[f64; 3]; 3] {
    let mut tensor = [[0.0f64; 3]; 3];
    for (axis, conductivity) in axes.iter().zip(principal) {
        for row in 0..3 {
            for column in 0..3 {
                tensor[row][column] =
                    (conductivity * axis[row]).mul_add(axis[column], tensor[row][column]);
            }
        }
    }
    tensor
}

fn homogenization_identity(
    stackup: &PcbStackup,
    principal: &PcbPrincipalConductivity,
    structural: &PcbStructuralBounds,
    separation_ratio: f64,
    tensor: &[[f64; 3]; 3],
) -> ContentHash {
    let mut payload = Vec::new();
    let mut push = |part: &[u8]| {
        payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
        payload.extend_from_slice(part);
    };
    push(&PCB_HOMOGENIZATION_SCHEMA_VERSION.to_le_bytes());
    push(stackup.id.as_bytes());
    push(
        &u64::try_from(stackup.layers.len())
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    for axis in stackup.frame.axes {
        for value in axis {
            push(&value.to_bits().to_le_bytes());
        }
    }
    push(
        &stackup
            .scale_separation
            .feature_size_m
            .to_bits()
            .to_le_bytes(),
    );
    push(
        &stackup
            .scale_separation
            .maximum_feature_to_thickness_ratio
            .to_bits()
            .to_le_bytes(),
    );
    for layer in &stackup.layers {
        push(layer.name.as_bytes());
        push(&layer.thickness_m.to_bits().to_le_bytes());
        for datum in [&layer.copper, &layer.matrix] {
            push(datum.material_state.chemistry.as_bytes());
            push(datum.material_state.phase.as_bytes());
            push(datum.material_state.process.as_bytes());
            push(&datum.material_state.revision.to_le_bytes());
            push(&datum.material_card.0);
            push(&datum.receipt.content_hash().0);
            push(&datum.nominal_w_mk.to_bits().to_le_bytes());
            push(&datum.lower_w_mk.to_bits().to_le_bytes());
            push(&datum.upper_w_mk.to_bits().to_le_bytes());
            push(&[u8::from(datum.uncertainty_complete)]);
        }
        push(layer.coverage.source_id.as_bytes());
        push(&layer.coverage.nominal.to_bits().to_le_bytes());
        push(&layer.coverage.lower.to_bits().to_le_bytes());
        push(&layer.coverage.upper.to_bits().to_le_bytes());
        push(layer.coverage.provenance.source.as_bytes());
        push(layer.coverage.provenance.license.as_bytes());
        match layer.coverage.provenance.artifact {
            None => push(&[0]),
            Some(artifact) => {
                push(&[1]);
                push(&artifact.0);
            }
        }
    }
    for value in principal
        .nominal_w_mk
        .into_iter()
        .chain(principal.lower_w_mk)
        .chain(principal.upper_w_mk)
    {
        push(&value.to_bits().to_le_bytes());
    }
    push(&structural.reuss_w_mk.to_bits().to_le_bytes());
    push(&structural.voigt_w_mk.to_bits().to_le_bytes());
    push(&separation_ratio.to_bits().to_le_bytes());
    for row in tensor {
        for value in row {
            push(&value.to_bits().to_le_bytes());
        }
    }
    hash_domain(PCB_HOMOGENIZATION_IDENTITY_DOMAIN, &payload)
}
