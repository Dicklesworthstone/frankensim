//! Card-backed surface radiation for the steady P1 conduction rung.
//!
//! Two deliberately separate models live here:
//!
//! - [`LinearizedSurfaceRadiation`] produces a Robin row only inside an
//!   explicit small-temperature-departure domain and reports its discrepancy
//!   from the nonlinear Stefan-Boltzmann flux at the evaluated point.
//! - [`GrayDiffuseEnclosure`] solves the opaque, gray, diffuse radiosity
//!   equations for an admitted view-factor matrix.  The outer coupling driver
//!   freezes the resulting face fluxes for one conduction solve and repeats
//!   until the declared surface-temperature criterion closes.
//!
//! General-geometry view-factor generation is not implemented here.  A caller
//! may provide an analytic or externally generated QMC matrix, but admission
//! checks row closure, non-negativity, and area-weighted reciprocity before the
//! matrix can influence a solve.

use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;
use fs_matdb::{MaterialCard, PropertyUsageReceipt, QueryPoint, SelectionPolicy, UncertaintyModel};

use crate::ConductionError;
use crate::assemble::DofMap;
use crate::bc::ThermalBc;
use crate::interface::ThermalInterfaces;
use crate::material::TEMPERATURE_AXIS;
use crate::mesh::{BoundaryFace, ConductionMesh};
use crate::solve::{
    ConductionProblem, ConductionSolution, InitialGuess, SolveConfig, solve, solve_with_interfaces,
};

/// CODATA exact SI value after the 2019 kelvin redefinition, W/(m² K⁴).
pub const STEFAN_BOLTZMANN_W_M2_K4: f64 = 5.670_374_419e-8;

/// Canonical `fs-matdb` property consumed by the radiation models.
pub const SURFACE_EMISSIVITY_PROPERTY: &str = "hemispherical-total-emissivity";

/// Hemispherical emissivity is dimensionless.
pub const EMISSIVITY_DIMS: fs_qty::Dims = fs_qty::Dims::NONE;

const VIEW_FACTOR_IDENTITY_DOMAIN: &str = "org.frankensim.fs-conduction.view-factors.v1";
const RADIATION_MESH_IDENTITY_DOMAIN: &str = "org.frankensim.fs-conduction.radiation-mesh.v1";
const MAX_ENCLOSURE_SURFACES: usize = 256;

fn radiation_error(
    surface: impl Into<String>,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> ConductionError {
    ConductionError::Radiation {
        surface: surface.into(),
        what: what.into(),
        fix: fix.into(),
    }
}

fn require_temperature(
    surface: &str,
    field: &'static str,
    value: f64,
) -> Result<(), ConductionError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(radiation_error(
            surface,
            format!("{field} {value} K must be finite and strictly positive"),
            "supply an absolute thermodynamic temperature in kelvin",
        ))
    }
}

fn fourth_power(value: f64) -> f64 {
    let square = value * value;
    square * square
}

/// One material-card-backed hemispherical total emissivity.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceEmissivity {
    value: f64,
    uncertainty: UncertaintyModel,
    temperature_k: f64,
    temperature_validity_k: Option<(f64, f64)>,
    card_identity: ContentHash,
    material_state: String,
    receipt: PropertyUsageReceipt,
}

impl SurfaceEmissivity {
    /// Resolve emissivity from an immutable as-manufactured material card.
    /// The material-state string, card identity, and exact selection receipt
    /// travel with every downstream radiation report.
    ///
    /// # Errors
    /// A material-query refusal for missing, ambiguous, or out-of-domain data;
    /// [`ConductionError::Dimensions`] for a non-dimensionless claim; or a
    /// typed radiation refusal when emissivity is outside `(0, 1]`.
    pub fn from_card(
        surface_name: &str,
        card: &MaterialCard,
        temperature_k: f64,
        policy: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        if surface_name.trim().is_empty() {
            return Err(radiation_error(
                "<unnamed>",
                "surface name is blank",
                "name the radiating surface so diagnostics and receipts identify it",
            ));
        }
        require_temperature(surface_name, "emissivity query temperature", temperature_k)?;
        let point = QueryPoint::new()
            .with(TEMPERATURE_AXIS, temperature_k)
            .map_err(|error| ConductionError::MaterialQuery {
                property: SURFACE_EMISSIVITY_PROPERTY.to_string(),
                temperature: temperature_k,
                upstream: error.to_string(),
            })?;
        let answer = card
            .claims()
            .query(SURFACE_EMISSIVITY_PROPERTY, &point, policy)
            .map_err(|error| ConductionError::MaterialQuery {
                property: SURFACE_EMISSIVITY_PROPERTY.to_string(),
                temperature: temperature_k,
                upstream: error.to_string(),
            })?;
        let sample = &answer.evidence.value;
        if sample.dims != EMISSIVITY_DIMS {
            return Err(ConductionError::Dimensions {
                context: format!(
                    "surface {surface_name:?} property {SURFACE_EMISSIVITY_PROPERTY:?}"
                ),
                expected: EMISSIVITY_DIMS.0,
                found: sample.dims.0,
            });
        }
        if !(sample.value.is_finite() && sample.value > 0.0 && sample.value <= 1.0) {
            return Err(radiation_error(
                surface_name,
                format!(
                    "hemispherical emissivity {} must lie in (0, 1]",
                    sample.value
                ),
                "attach one in-domain surface-finish claim with physical emissivity",
            ));
        }
        Ok(Self {
            value: sample.value,
            uncertainty: sample.uncertainty.clone(),
            temperature_k,
            temperature_validity_k: answer.evidence.model.validity.bound(TEMPERATURE_AXIS),
            card_identity: card.content_hash(),
            material_state: card.id().to_string(),
            receipt: answer.receipt,
        })
    }

    /// Dimensionless emissivity in `(0, 1]`.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// Temperature at which the property card was queried, K.
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Selected claim's closed temperature-validity interval, K.
    ///
    /// `None` means the selected claim did not constrain temperature.
    #[must_use]
    pub const fn temperature_validity_k(&self) -> Option<(f64, f64)> {
        self.temperature_validity_k
    }

    /// Source uncertainty; `Unstated` remains an explicit unknown.
    #[must_use]
    pub const fn uncertainty(&self) -> &UncertaintyModel {
        &self.uncertainty
    }

    /// Immutable material-card identity, including process/surface state.
    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    /// Human-readable named material state bound into the card identity.
    #[must_use]
    pub fn material_state(&self) -> &str {
        &self.material_state
    }

    /// Exact property-use receipt.
    #[must_use]
    pub const fn receipt(&self) -> &PropertyUsageReceipt {
        &self.receipt
    }

    fn validate_temperature(
        &self,
        surface_name: &str,
        field: &'static str,
        temperature_k: f64,
    ) -> Result<(), ConductionError> {
        require_temperature(surface_name, field, temperature_k)?;
        if let Some((lower_k, upper_k)) = self.temperature_validity_k
            && !(lower_k..=upper_k).contains(&temperature_k)
        {
            return Err(radiation_error(
                surface_name,
                format!(
                    "{field} {temperature_k} K is outside the selected emissivity claim's \
                     retained validity [{lower_k}, {upper_k}] K"
                ),
                "supply an in-domain surface temperature or resolve a different emissivity claim",
            ));
        }
        Ok(())
    }

    fn half_width(&self) -> Option<(f64, f64)> {
        match self.uncertainty {
            UncertaintyModel::Unstated => None,
            UncertaintyModel::HalfWidth {
                half_width,
                confidence,
            } => Some((half_width, confidence)),
            UncertaintyModel::RelativeHalfWidth {
                fraction,
                confidence,
            } => Some((fraction * self.value.abs(), confidence)),
        }
    }
}

/// Evaluation of the linearized surface-to-ambient rung at one surface
/// temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearizedRadiationPoint {
    /// Robin row carrying `h_rad` and the declared ambient temperature.
    pub boundary: ThermalBc,
    /// `4 eps sigma T_m^3`, W/(m² K).
    pub h_rad_w_m2k: f64,
    /// `h_rad (T_s - T_ambient)`, W/m², positive outward.
    pub linearized_outward_flux_w_m2: f64,
    /// `eps sigma (T_s^4 - T_ambient^4)`, W/m², positive outward.
    pub nonlinear_outward_flux_w_m2: f64,
    /// Absolute pointwise difference between the two rungs, W/m².
    pub discrepancy_w_m2: f64,
    /// Propagated half-width on `h_rad`, when emissivity uncertainty is stated.
    pub h_rad_half_width_w_m2k: Option<f64>,
    /// Confidence attached to `h_rad_half_width_w_m2k`.
    pub uncertainty_confidence: Option<f64>,
}

/// Small-departure surface-to-ambient radiation represented as a Robin row.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearizedSurfaceRadiation {
    surface_name: String,
    emissivity: SurfaceEmissivity,
    mean_temperature_k: f64,
    ambient_temperature_k: f64,
    max_abs_departure_k: f64,
}

impl LinearizedSurfaceRadiation {
    /// Construct a card-backed linearization domain.
    ///
    /// Both the ambient and every later surface evaluation must lie within
    /// `max_abs_departure_k` of `mean_temperature_k`; otherwise the model
    /// refuses rather than silently stretching the small-delta assumption.
    pub fn new(
        surface_name: impl Into<String>,
        emissivity: SurfaceEmissivity,
        mean_temperature_k: f64,
        ambient_temperature_k: f64,
        max_abs_departure_k: f64,
    ) -> Result<Self, ConductionError> {
        let surface_name = surface_name.into();
        if surface_name.trim().is_empty() {
            return Err(radiation_error(
                "<unnamed>",
                "surface name is blank",
                "name the linearized radiation model",
            ));
        }
        emissivity.validate_temperature(
            &surface_name,
            "linearization mean temperature",
            mean_temperature_k,
        )?;
        emissivity.validate_temperature(
            &surface_name,
            "ambient temperature",
            ambient_temperature_k,
        )?;
        if !(max_abs_departure_k.is_finite() && max_abs_departure_k > 0.0) {
            return Err(radiation_error(
                &surface_name,
                format!(
                    "linearization half-width {max_abs_departure_k} K must be finite and positive"
                ),
                "declare the small-temperature-departure validity half-width",
            ));
        }
        if (ambient_temperature_k - mean_temperature_k).abs() > max_abs_departure_k {
            return Err(radiation_error(
                &surface_name,
                "ambient temperature lies outside the declared linearization domain",
                "move the expansion point, widen it only with discrepancy evidence, or use full T^4 radiation",
            ));
        }
        Ok(Self {
            surface_name,
            emissivity,
            mean_temperature_k,
            ambient_temperature_k,
            max_abs_departure_k,
        })
    }

    /// Evaluate and lower this rung to a Robin row.
    pub fn evaluate(
        &self,
        surface_temperature_k: f64,
    ) -> Result<LinearizedRadiationPoint, ConductionError> {
        self.emissivity.validate_temperature(
            &self.surface_name,
            "surface temperature",
            surface_temperature_k,
        )?;
        if (surface_temperature_k - self.mean_temperature_k).abs() > self.max_abs_departure_k {
            return Err(radiation_error(
                &self.surface_name,
                format!(
                    "surface temperature {surface_temperature_k} K is outside [{}, {}] K",
                    self.mean_temperature_k - self.max_abs_departure_k,
                    self.mean_temperature_k + self.max_abs_departure_k
                ),
                "use the nonlinear gray-diffuse rung or a justified wider linearization card",
            ));
        }
        let h_rad = 4.0
            * self.emissivity.value
            * STEFAN_BOLTZMANN_W_M2_K4
            * self.mean_temperature_k
            * self.mean_temperature_k
            * self.mean_temperature_k;
        let linearized = h_rad * (surface_temperature_k - self.ambient_temperature_k);
        let nonlinear = self.emissivity.value
            * STEFAN_BOLTZMANN_W_M2_K4
            * (fourth_power(surface_temperature_k) - fourth_power(self.ambient_temperature_k));
        let (h_rad_half_width_w_m2k, uncertainty_confidence) =
            self.emissivity
                .half_width()
                .map_or((None, None), |(half_width, confidence)| {
                    (
                        Some(h_rad * half_width / self.emissivity.value),
                        Some(confidence),
                    )
                });
        Ok(LinearizedRadiationPoint {
            boundary: ThermalBc::robin(h_rad, self.ambient_temperature_k)?,
            h_rad_w_m2k: h_rad,
            linearized_outward_flux_w_m2: linearized,
            nonlinear_outward_flux_w_m2: nonlinear,
            discrepancy_w_m2: (linearized - nonlinear).abs(),
            h_rad_half_width_w_m2k,
            uncertainty_confidence,
        })
    }

    /// Emissivity card and receipt retained by this model.
    #[must_use]
    pub const fn emissivity(&self) -> &SurfaceEmissivity {
        &self.emissivity
    }

    /// Declared linearization center, K.
    #[must_use]
    pub const fn mean_temperature_k(&self) -> f64 {
        self.mean_temperature_k
    }

    /// Declared validity half-width, K.
    #[must_use]
    pub const fn max_abs_departure_k(&self) -> f64 {
        self.max_abs_departure_k
    }
}

/// Evidence attached to an admitted view-factor matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ViewFactorEvidence {
    /// A named analytic limiting or closed-form geometry.
    Analytic {
        /// Stable formula/geometry identifier.
        geometry: String,
    },
    /// An externally generated randomized/QMC estimate.  This crate validates
    /// the matrix and retains the generation coordinates; it does not claim to
    /// implement the geometry sampler.
    ExternalQmc {
        /// Counter-based seed used by the external generator.
        seed: u64,
        /// Number of admitted rays/samples.
        samples: u64,
        /// Stable generator/version identifier.
        generator: String,
    },
}

/// Numerical admission tolerances for a view-factor matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViewFactorTolerance {
    /// Maximum absolute deviation of any row sum from one.
    pub row_sum_abs: f64,
    /// Maximum relative area-weighted reciprocity residual.
    pub reciprocity_rel: f64,
}

impl Default for ViewFactorTolerance {
    fn default() -> Self {
        Self {
            row_sum_abs: 1.0e-12,
            reciprocity_rel: 1.0e-12,
        }
    }
}

/// Admitted enclosure view factors with areas and evidence bound into a
/// stable content identity.
#[derive(Debug, Clone, PartialEq)]
pub struct ViewFactorMatrix {
    areas_m2: Vec<f64>,
    factors: Vec<Vec<f64>>,
    row_sums: Vec<f64>,
    max_reciprocity_residual: f64,
    evidence: ViewFactorEvidence,
    identity: ContentHash,
}

impl ViewFactorMatrix {
    /// Admit a complete matrix after row-closure and reciprocity checks.
    pub fn admit(
        areas_m2: Vec<f64>,
        factors: Vec<Vec<f64>>,
        evidence: ViewFactorEvidence,
        tolerance: ViewFactorTolerance,
    ) -> Result<Self, ConductionError> {
        let n = areas_m2.len();
        if !(2..=MAX_ENCLOSURE_SURFACES).contains(&n) {
            return Err(radiation_error(
                "view-factors",
                format!("surface count {n} is outside 2..={MAX_ENCLOSURE_SURFACES}"),
                "supply a finite, nonempty enclosure within the declared resource cap",
            ));
        }
        match &evidence {
            ViewFactorEvidence::Analytic { geometry } if geometry.trim().is_empty() => {
                return Err(radiation_error(
                    "view-factors",
                    "analytic view-factor evidence has a blank geometry/formula identifier",
                    "name the exact analytic geometry or formula revision",
                ));
            }
            ViewFactorEvidence::ExternalQmc {
                samples, generator, ..
            } if *samples == 0 || generator.trim().is_empty() => {
                return Err(radiation_error(
                    "view-factors",
                    "external-QMC evidence requires a positive sample count and generator identifier",
                    "retain the explicit seed, positive ray count, and generator/version label",
                ));
            }
            _ => {}
        }
        if !(tolerance.row_sum_abs.is_finite()
            && tolerance.row_sum_abs >= 0.0
            && tolerance.row_sum_abs < 1.0
            && tolerance.reciprocity_rel.is_finite()
            && tolerance.reciprocity_rel >= 0.0
            && tolerance.reciprocity_rel < 1.0)
        {
            return Err(radiation_error(
                "view-factors",
                "view-factor tolerances must be finite and lie in [0, 1)",
                "declare bounded row-sum and reciprocity tolerances",
            ));
        }
        if factors.len() != n || factors.iter().any(|row| row.len() != n) {
            return Err(radiation_error(
                "view-factors",
                format!("matrix shape is not {n} by {n}"),
                "supply one factor from every surface to every surface",
            ));
        }
        for (i, &area) in areas_m2.iter().enumerate() {
            if !(area.is_finite() && area > 0.0) {
                return Err(radiation_error(
                    format!("surface-{i}"),
                    format!("area {area} m^2 must be finite and positive"),
                    "supply the same physical surface area used by the conduction trace",
                ));
            }
        }
        let mut row_sums = Vec::with_capacity(n);
        for (i, row) in factors.iter().enumerate() {
            let mut sum = 0.0;
            for (j, &factor) in row.iter().enumerate() {
                if !(factor.is_finite() && factor >= 0.0 && factor <= 1.0) {
                    return Err(radiation_error(
                        format!("view-factor-{i}-{j}"),
                        format!("factor {factor} must lie in [0, 1]"),
                        "regenerate the enclosure matrix without clipping or NaNs",
                    ));
                }
                sum += factor;
            }
            if (sum - 1.0).abs() > tolerance.row_sum_abs {
                return Err(radiation_error(
                    format!("view-factor-row-{i}"),
                    format!(
                        "row sum {sum} differs from one beyond {}",
                        tolerance.row_sum_abs
                    ),
                    "close the enclosure or provide explicit environment surfaces",
                ));
            }
            row_sums.push(sum);
        }
        let mut max_reciprocity_residual = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let lhs = areas_m2[i] * factors[i][j];
                let rhs = areas_m2[j] * factors[j][i];
                let scale = lhs.abs().max(rhs.abs()).max(f64::MIN_POSITIVE);
                let residual = (lhs - rhs).abs() / scale;
                max_reciprocity_residual = max_reciprocity_residual.max(residual);
                if residual > tolerance.reciprocity_rel {
                    return Err(radiation_error(
                        format!("view-factor-{i}-{j}"),
                        format!(
                            "A_i F_ij = {lhs} and A_j F_ji = {rhs} differ by relative {residual}"
                        ),
                        "use an area-consistent reciprocal matrix or retain a larger justified tolerance",
                    ));
                }
            }
        }
        let identity = view_factor_identity(&areas_m2, &factors, &evidence, tolerance);
        Ok(Self {
            areas_m2,
            factors,
            row_sums,
            max_reciprocity_residual,
            evidence,
            identity,
        })
    }

    /// Infinite parallel plates of equal area: `F12 = F21 = 1` exactly.
    pub fn infinite_parallel_plates(area_m2: f64) -> Result<Self, ConductionError> {
        Self::admit(
            vec![area_m2, area_m2],
            vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            ViewFactorEvidence::Analytic {
                geometry: "infinite-parallel-plates".to_string(),
            },
            ViewFactorTolerance {
                row_sum_abs: 0.0,
                reciprocity_rel: 0.0,
            },
        )
    }

    /// Number of enclosure surfaces.
    #[must_use]
    pub fn len(&self) -> usize {
        self.areas_m2.len()
    }

    /// True only for an impossible admitted state; supplied for ordinary
    /// collection ergonomics.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.areas_m2.is_empty()
    }

    /// Areas in matrix order, m².
    #[must_use]
    pub fn areas_m2(&self) -> &[f64] {
        &self.areas_m2
    }

    /// Matrix rows in deterministic surface order.
    #[must_use]
    pub fn factors(&self) -> &[Vec<f64>] {
        &self.factors
    }

    /// Row sums retained from admission.
    #[must_use]
    pub fn row_sums(&self) -> &[f64] {
        &self.row_sums
    }

    /// Largest admitted area-weighted reciprocity residual.
    #[must_use]
    pub const fn max_reciprocity_residual(&self) -> f64 {
        self.max_reciprocity_residual
    }

    /// Formula or external-QMC evidence tag.
    #[must_use]
    pub const fn evidence(&self) -> &ViewFactorEvidence {
        &self.evidence
    }

    /// Stable identity over areas, matrix entries, tolerances, and evidence.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

fn view_factor_identity(
    areas: &[f64],
    factors: &[Vec<f64>],
    evidence: &ViewFactorEvidence,
    tolerance: ViewFactorTolerance,
) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(areas.len() as u64).to_le_bytes());
    for &area in areas {
        bytes.extend_from_slice(&area.to_bits().to_le_bytes());
    }
    for row in factors {
        for &factor in row {
            bytes.extend_from_slice(&factor.to_bits().to_le_bytes());
        }
    }
    bytes.extend_from_slice(&tolerance.row_sum_abs.to_bits().to_le_bytes());
    bytes.extend_from_slice(&tolerance.reciprocity_rel.to_bits().to_le_bytes());
    match evidence {
        ViewFactorEvidence::Analytic { geometry } => {
            bytes.push(0);
            bytes.extend_from_slice(&(geometry.len() as u64).to_le_bytes());
            bytes.extend_from_slice(geometry.as_bytes());
        }
        ViewFactorEvidence::ExternalQmc {
            seed,
            samples,
            generator,
        } => {
            bytes.push(1);
            bytes.extend_from_slice(&seed.to_le_bytes());
            bytes.extend_from_slice(&samples.to_le_bytes());
            bytes.extend_from_slice(&(generator.len() as u64).to_le_bytes());
            bytes.extend_from_slice(generator.as_bytes());
        }
    }
    hash_domain(VIEW_FACTOR_IDENTITY_DOMAIN, &bytes)
}

/// One named P1 boundary trace and its card-backed emissivity.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationSurface {
    name: String,
    face_slots: Vec<usize>,
    area_m2: f64,
    emissivity: SurfaceEmissivity,
    mesh_identity: ContentHash,
}

impl RadiationSurface {
    /// Select a nonempty boundary trace.  Slots are retained in mesh order.
    pub fn new(
        mesh: &ConductionMesh,
        name: impl Into<String>,
        select: impl Fn(&BoundaryFace) -> bool,
        emissivity: SurfaceEmissivity,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(radiation_error(
                "<unnamed>",
                "radiation surface name is blank",
                "name the trace so coupling reports are stable",
            ));
        }
        let face_slots = mesh
            .boundary()
            .iter()
            .enumerate()
            .filter_map(|(slot, face)| select(face).then_some(slot))
            .collect::<Vec<_>>();
        if face_slots.is_empty() {
            return Err(radiation_error(
                &name,
                "surface selector matched no boundary faces",
                "bind the radiation card to a nonempty exterior trace",
            ));
        }
        let area_m2 = face_slots
            .iter()
            .map(|&slot| mesh.boundary()[slot].area)
            .sum::<f64>();
        if !(area_m2.is_finite() && area_m2 > 0.0) {
            return Err(radiation_error(
                &name,
                format!("selected area {area_m2} m^2 must be finite and positive"),
                "repair the boundary geometry",
            ));
        }
        Ok(Self {
            name,
            face_slots,
            area_m2,
            emissivity,
            mesh_identity: radiation_mesh_identity(mesh),
        })
    }

    /// Stable surface name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Boundary-face slots in ascending mesh order.
    #[must_use]
    pub fn face_slots(&self) -> &[usize] {
        &self.face_slots
    }

    /// Integrated P1 trace area, m².
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Card-backed emissivity and receipt.
    #[must_use]
    pub const fn emissivity(&self) -> &SurfaceEmissivity {
        &self.emissivity
    }

    /// Exact area average of the P1 trace: each triangular face average is
    /// the arithmetic mean of its three nodal values.
    pub fn mean_temperature(
        &self,
        mesh: &ConductionMesh,
        temperature: &[f64],
    ) -> Result<f64, ConductionError> {
        self.validate_for(mesh)?;
        if temperature.len() != mesh.vertex_count() {
            return Err(ConductionError::FieldLength {
                field: "radiation surface temperature",
                expected: mesh.vertex_count(),
                found: temperature.len(),
            });
        }
        let mut integral = 0.0;
        for &slot in &self.face_slots {
            let face = &mesh.boundary()[slot];
            let average = face
                .vertices
                .iter()
                .map(|&vertex| temperature[vertex as usize])
                .sum::<f64>()
                / 3.0;
            require_temperature(&self.name, "surface nodal average", average)?;
            integral = face.area.mul_add(average, integral);
        }
        Ok(integral / self.area_m2)
    }

    fn validate_for(&self, mesh: &ConductionMesh) -> Result<(), ConductionError> {
        if self.mesh_identity != radiation_mesh_identity(mesh) {
            return Err(radiation_error(
                &self.name,
                "radiation surface was bound to a different conduction mesh",
                "rebuild RadiationSurface from the exact mesh used by this solve",
            ));
        }
        Ok(())
    }
}

fn radiation_mesh_identity(mesh: &ConductionMesh) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(mesh.vertex_count() as u64).to_le_bytes());
    bytes.extend_from_slice(&(mesh.element_count() as u64).to_le_bytes());
    for position in mesh.positions() {
        for coordinate in position {
            bytes.extend_from_slice(&coordinate.to_bits().to_le_bytes());
        }
    }
    for tet in &mesh.complex().tets {
        for vertex in tet {
            bytes.extend_from_slice(&vertex.to_le_bytes());
        }
    }
    hash_domain(RADIATION_MESH_IDENTITY_DOMAIN, &bytes)
}

/// Nonlinear gray-diffuse enclosure with fixed view factors.
#[derive(Debug, Clone, PartialEq)]
pub struct GrayDiffuseEnclosure {
    surfaces: Vec<RadiationSurface>,
    view_factors: ViewFactorMatrix,
}

impl GrayDiffuseEnclosure {
    /// Bind surface traces to the view-factor order.
    pub fn new(
        surfaces: Vec<RadiationSurface>,
        view_factors: ViewFactorMatrix,
    ) -> Result<Self, ConductionError> {
        if surfaces.len() != view_factors.len() {
            return Err(radiation_error(
                "enclosure",
                format!(
                    "{} surfaces cannot bind a {} by {} view-factor matrix",
                    surfaces.len(),
                    view_factors.len(),
                    view_factors.len()
                ),
                "supply exactly one trace and emissivity per matrix row",
            ));
        }
        let mut names = BTreeSet::new();
        let mut slots = BTreeSet::new();
        let mesh_identity = surfaces.first().map(|surface| surface.mesh_identity);
        for (index, surface) in surfaces.iter().enumerate() {
            if Some(surface.mesh_identity) != mesh_identity {
                return Err(radiation_error(
                    &surface.name,
                    "enclosure surfaces were bound to different conduction meshes",
                    "construct every enclosure surface from one exact mesh",
                ));
            }
            if !names.insert(surface.name.clone()) {
                return Err(radiation_error(
                    &surface.name,
                    "surface name is duplicated",
                    "use stable unique names",
                ));
            }
            for &slot in &surface.face_slots {
                if !slots.insert(slot) {
                    return Err(radiation_error(
                        &surface.name,
                        format!("boundary-face slot {slot} belongs to two radiation surfaces"),
                        "partition the radiating trace without overlap",
                    ));
                }
            }
            let declared = view_factors.areas_m2[index];
            let scale = declared
                .abs()
                .max(surface.area_m2.abs())
                .max(f64::MIN_POSITIVE);
            if (declared - surface.area_m2).abs() > 1.0e-12 * scale {
                return Err(radiation_error(
                    &surface.name,
                    format!(
                        "view-factor area {declared} m^2 differs from mesh trace area {} m^2",
                        surface.area_m2
                    ),
                    "generate the matrix from the same surface partition used by conduction",
                ));
            }
        }
        Ok(Self {
            surfaces,
            view_factors,
        })
    }

    /// Bound surfaces in deterministic matrix order.
    #[must_use]
    pub fn surfaces(&self) -> &[RadiationSurface] {
        &self.surfaces
    }

    /// Admitted view factors and their content identity.
    #[must_use]
    pub const fn view_factors(&self) -> &ViewFactorMatrix {
        &self.view_factors
    }

    /// Solve `J_i - (1-eps_i) sum_j F_ij J_j = eps_i sigma T_i^4`.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        surface_temperatures_k: &[f64],
    ) -> Result<RadiosityReport, ConductionError> {
        let n = self.surfaces.len();
        if surface_temperatures_k.len() != n {
            return Err(radiation_error(
                "enclosure",
                format!(
                    "{} surface temperatures supplied for {n} surfaces",
                    surface_temperatures_k.len()
                ),
                "supply temperatures in view-factor row order",
            ));
        }
        let mut matrix = vec![vec![0.0; n]; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                stage: "radiation-radiosity-assemble",
                at: i,
            })?;
            self.surfaces[i].emissivity.validate_temperature(
                self.surfaces[i].name(),
                "gray-diffuse surface temperature",
                surface_temperatures_k[i],
            )?;
            let emissivity = self.surfaces[i].emissivity.value;
            for j in 0..n {
                matrix[i][j] = if i == j { 1.0 } else { 0.0 }
                    - (1.0 - emissivity) * self.view_factors.factors[i][j];
            }
            rhs[i] =
                emissivity * STEFAN_BOLTZMANN_W_M2_K4 * fourth_power(surface_temperatures_k[i]);
        }
        let radiosity_w_m2 = solve_dense(cx, &matrix, &rhs)?;
        let mut irradiation_w_m2 = vec![0.0; n];
        let mut net_outward_flux_w_m2 = vec![0.0; n];
        let mut net_outward_heat_w = vec![0.0; n];
        let mut heat_sum = 0.0;
        let mut heat_scale = 0.0;
        let mut linear_residual_max = 0.0f64;
        for i in 0..n {
            let mut irradiation = 0.0;
            let mut lhs = 0.0;
            for j in 0..n {
                irradiation =
                    self.view_factors.factors[i][j].mul_add(radiosity_w_m2[j], irradiation);
                lhs = matrix[i][j].mul_add(radiosity_w_m2[j], lhs);
            }
            irradiation_w_m2[i] = irradiation;
            net_outward_flux_w_m2[i] = radiosity_w_m2[i] - irradiation;
            net_outward_heat_w[i] = self.surfaces[i].area_m2 * net_outward_flux_w_m2[i];
            heat_sum += net_outward_heat_w[i];
            heat_scale += net_outward_heat_w[i].abs();
            linear_residual_max = linear_residual_max.max((lhs - rhs[i]).abs());
        }
        Ok(RadiosityReport {
            surface_temperatures_k: surface_temperatures_k.to_vec(),
            radiosity_w_m2,
            irradiation_w_m2,
            net_outward_flux_w_m2,
            net_outward_heat_w,
            enclosure_energy_closure_w: heat_sum,
            enclosure_energy_scale_w: heat_scale.max(f64::MIN_POSITIVE),
            linear_residual_max_w_m2: linear_residual_max,
            view_factor_identity: self.view_factors.identity,
            view_factor_row_sums: self.view_factors.row_sums.clone(),
            max_reciprocity_residual: self.view_factors.max_reciprocity_residual,
            emissivity_uncertainty_complete: self.surfaces.iter().all(|surface| {
                !matches!(surface.emissivity.uncertainty, UncertaintyModel::Unstated)
            }),
        })
    }

    fn surface_temperatures(
        &self,
        mesh: &ConductionMesh,
        temperature: &[f64],
    ) -> Result<Vec<f64>, ConductionError> {
        self.surfaces
            .iter()
            .map(|surface| surface.mean_temperature(mesh, temperature))
            .collect()
    }

    fn validate_for(&self, mesh: &ConductionMesh) -> Result<(), ConductionError> {
        for surface in &self.surfaces {
            surface.validate_for(mesh)?;
        }
        Ok(())
    }
}

/// Evidence from one gray-diffuse radiosity solve.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiosityReport {
    /// Area-averaged P1 surface temperatures used by this solve, K.
    pub surface_temperatures_k: Vec<f64>,
    /// Surface radiosities `J_i`, W/m².
    pub radiosity_w_m2: Vec<f64>,
    /// Surface irradiations `G_i = sum_j F_ij J_j`, W/m².
    pub irradiation_w_m2: Vec<f64>,
    /// `J_i - G_i`, W/m², positive leaving solid `i`.
    pub net_outward_flux_w_m2: Vec<f64>,
    /// `A_i (J_i - G_i)`, W.
    pub net_outward_heat_w: Vec<f64>,
    /// Sum of all enclosure heat rates, W; zero for a closed reciprocal matrix
    /// up to algebraic rounding.
    pub enclosure_energy_closure_w: f64,
    /// Sum of absolute enclosure heat rates, W.
    pub enclosure_energy_scale_w: f64,
    /// Largest absolute residual of the dense radiosity system, W/m².
    pub linear_residual_max_w_m2: f64,
    /// Identity of the exact admitted view-factor matrix.
    pub view_factor_identity: ContentHash,
    /// Row sums retained for structured diagnostics.
    pub view_factor_row_sums: Vec<f64>,
    /// Largest admitted area-weighted reciprocity residual.
    pub max_reciprocity_residual: f64,
    /// True only when every emissivity source stated a quantitative width.
    /// This does not claim that the nonlinear output uncertainty was propagated.
    pub emissivity_uncertainty_complete: bool,
}

impl RadiosityReport {
    /// Relative enclosure energy closure.
    #[must_use]
    pub fn relative_energy_closure(&self) -> f64 {
        self.enclosure_energy_closure_w.abs() / self.enclosure_energy_scale_w
    }
}

/// Outer fixed-point controls for conduction-radiosity coupling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupledRadiationConfig {
    /// Surface-temperature relative tolerance.
    pub surface_temperature_rtol: f64,
    /// Surface-temperature absolute tolerance, K.
    pub surface_temperature_atol_k: f64,
    /// Under-relaxation applied to the temperatures that drive the next
    /// radiosity solve, in `(0, 1]`.
    pub relaxation: f64,
    /// Outer-iteration budget.
    pub max_iterations: usize,
}

impl Default for CoupledRadiationConfig {
    fn default() -> Self {
        Self {
            surface_temperature_rtol: 1.0e-9,
            surface_temperature_atol_k: 1.0e-9,
            relaxation: 0.5,
            max_iterations: 100,
        }
    }
}

impl CoupledRadiationConfig {
    fn validate(self) -> Result<(), ConductionError> {
        if !(self.surface_temperature_rtol.is_finite()
            && self.surface_temperature_rtol >= 0.0
            && self.surface_temperature_atol_k.is_finite()
            && self.surface_temperature_atol_k >= 0.0
            && self.relaxation.is_finite()
            && self.relaxation > 0.0
            && self.relaxation <= 1.0
            && self.max_iterations > 0)
        {
            return Err(radiation_error(
                "coupling",
                "coupling tolerances/relaxation/budget are inadmissible",
                "use finite non-negative tolerances, relaxation in (0, 1], and a positive budget",
            ));
        }
        Ok(())
    }
}

/// Outer-coupling evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledRadiationReport {
    /// Outer fixed-point iterations performed.
    pub iterations: usize,
    /// Maximum raw surface-temperature update per iteration, K.
    pub surface_update_history_k: Vec<f64>,
    /// Threshold used on the accepted iteration, K.
    pub final_threshold_k: f64,
    /// Last radiosity solve whose fluxes were applied to conduction.
    pub radiosity: RadiosityReport,
}

/// Coupled conduction field plus enclosure evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledRadiationSolution {
    /// Accepted inner conduction solution.
    pub conduction: ConductionSolution,
    /// Outer fixed-point and radiosity evidence.
    pub radiation: CoupledRadiationReport,
}

/// Couple a gray-diffuse enclosure to steady conduction with a declared outer
/// fixed point. Radiation faces must be part of the base boundary's explicit
/// adiabatic remainder; this prevents accidental replacement of Dirichlet,
/// Neumann, or convective rows.
pub fn solve_with_gray_diffuse_enclosure(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    enclosure: &GrayDiffuseEnclosure,
    conduction_config: SolveConfig,
    coupling_config: CoupledRadiationConfig,
) -> Result<CoupledRadiationSolution, ConductionError> {
    coupling_config.validate()?;
    enclosure.validate_for(problem.mesh)?;
    let mut conduction = run_conduction(cx, problem, interfaces, conduction_config.clone())?;
    let mut driving_temperatures =
        enclosure.surface_temperatures(problem.mesh, &conduction.temperature)?;
    let mut update_history = Vec::new();

    for iteration in 0..coupling_config.max_iterations {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "radiation-coupling",
            at: iteration,
        })?;
        let radiosity = enclosure.solve(cx, &driving_temperatures)?;
        let rows = enclosure
            .surfaces
            .iter()
            .zip(&radiosity.net_outward_flux_w_m2)
            .map(|(surface, &flux)| (surface.name.clone(), surface.face_slots.clone(), flux))
            .collect::<Vec<_>>();
        let overlaid = problem
            .boundary
            .with_uniform_outward_flux_overlays(problem.mesh, &rows)?;
        let dofs = DofMap::new(&overlaid, problem.mesh.vertex_count())?;
        let mut next_config = conduction_config.clone();
        next_config.initial = InitialGuess::Free(dofs.gather(&conduction.temperature));
        let next_problem = ConductionProblem { element_materials: None, mesh: problem.mesh,
        boundary: &overlaid,
        material: problem.material,
        source: problem.source, };
        let next = run_conduction(cx, next_problem, interfaces, next_config)?;
        let actual_temperatures =
            enclosure.surface_temperatures(problem.mesh, &next.temperature)?;
        for (surface, &temperature_k) in enclosure.surfaces.iter().zip(&actual_temperatures) {
            surface.emissivity.validate_temperature(
                surface.name(),
                "coupled surface temperature",
                temperature_k,
            )?;
        }
        let update = actual_temperatures
            .iter()
            .zip(&driving_temperatures)
            .map(|(next, prior)| (next - prior).abs())
            .fold(0.0f64, f64::max);
        let scale = actual_temperatures
            .iter()
            .copied()
            .map(f64::abs)
            .fold(0.0f64, f64::max);
        let threshold = coupling_config
            .surface_temperature_rtol
            .mul_add(scale, coupling_config.surface_temperature_atol_k);
        update_history.push(update);
        conduction = next;
        if update <= threshold {
            return Ok(CoupledRadiationSolution {
                conduction,
                radiation: CoupledRadiationReport {
                    iterations: iteration + 1,
                    surface_update_history_k: update_history,
                    final_threshold_k: threshold,
                    radiosity,
                },
            });
        }
        for (driving, actual) in driving_temperatures.iter_mut().zip(actual_temperatures) {
            *driving = coupling_config
                .relaxation
                .mul_add(actual, (1.0 - coupling_config.relaxation) * *driving);
        }
    }
    Err(radiation_error(
        "coupling",
        format!(
            "surface-temperature fixed point did not converge in {} iterations; final update was {} K",
            coupling_config.max_iterations,
            update_history.last().copied().unwrap_or(f64::INFINITY)
        ),
        "increase the declared budget, reduce relaxation, or use a monolithic nonlinear radiation operator",
    ))
}

fn run_conduction(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    config: SolveConfig,
) -> Result<ConductionSolution, ConductionError> {
    match interfaces {
        Some(interfaces) => solve_with_interfaces(cx, problem, interfaces, config),
        None => solve(cx, problem, config),
    }
}

fn solve_dense(cx: &Cx<'_>, matrix: &[Vec<f64>], rhs: &[f64]) -> Result<Vec<f64>, ConductionError> {
    let n = rhs.len();
    let mut a = matrix.to_vec();
    let mut b = rhs.to_vec();
    let scale = a
        .iter()
        .flatten()
        .copied()
        .map(f64::abs)
        .fold(0.0f64, f64::max)
        .max(f64::MIN_POSITIVE);
    for column in 0..n {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "radiation-radiosity-factor",
            at: column,
        })?;
        let mut pivot = column;
        let mut pivot_abs = a[column][column].abs();
        for (row, values) in a.iter().enumerate().skip(column + 1) {
            let candidate = values[column].abs();
            if candidate > pivot_abs {
                pivot = row;
                pivot_abs = candidate;
            }
        }
        if !pivot_abs.is_finite() || pivot_abs <= 64.0 * f64::EPSILON * scale {
            return Err(radiation_error(
                "radiosity-system",
                format!("pivot {column} has magnitude {pivot_abs}"),
                "repair the emissivity/view-factor model; the radiosity operator is singular or ill-conditioned",
            ));
        }
        if pivot != column {
            a.swap(column, pivot);
            b.swap(column, pivot);
        }
        let pivot_row = a[column].clone();
        let pivot_rhs = b[column];
        for row in (column + 1)..n {
            let factor = a[row][column] / pivot_row[column];
            a[row][column] = 0.0;
            for entry in (column + 1)..n {
                a[row][entry] = (-factor).mul_add(pivot_row[entry], a[row][entry]);
            }
            b[row] = (-factor).mul_add(pivot_rhs, b[row]);
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut value = b[row];
        for (column, &known) in x.iter().enumerate().skip(row + 1) {
            value = (-a[row][column]).mul_add(known, value);
        }
        x[row] = value / a[row][row];
        if !x[row].is_finite() {
            return Err(radiation_error(
                "radiosity-system",
                format!("solution component {row} is not finite"),
                "repair the enclosure inputs",
            ));
        }
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::{ViewFactorEvidence, ViewFactorMatrix, ViewFactorTolerance};

    #[test]
    fn view_factor_identity_moves_with_evidence() {
        let factors = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let a = ViewFactorMatrix::admit(
            vec![1.0, 1.0],
            factors.clone(),
            ViewFactorEvidence::Analytic {
                geometry: "parallel-plates-a".to_string(),
            },
            ViewFactorTolerance::default(),
        )
        .expect("matrix");
        let b = ViewFactorMatrix::admit(
            vec![1.0, 1.0],
            factors,
            ViewFactorEvidence::Analytic {
                geometry: "parallel-plates-b".to_string(),
            },
            ViewFactorTolerance::default(),
        )
        .expect("matrix");
        assert_ne!(a.identity(), b.identity());
    }
}
