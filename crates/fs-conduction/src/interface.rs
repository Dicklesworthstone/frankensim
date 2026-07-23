//! Thermal contact resistance on explicitly paired, coincident P1 faces.
//!
//! A contact interface is represented by duplicated vertices: each solid owns
//! its temperature trace, so the discrete field may jump.  For area-specific
//! resistance `R''` and conductance `g = 1/R''`, each paired face contributes
//!
//! ```text
//! integral_Gamma g (T_a - T_b) (v_a - v_b) dA
//! ```
//!
//! using the exact P1 triangle mass matrix.  Every geometrically coincident
//! boundary-face pair must be declared.  Missing bindings refuse instead of
//! silently becoming perfect contact or an adiabatic gap.

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::{
    InterfaceSystemCard, PropertyUsageReceipt, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_sparse::Coo;

use crate::ConductionError;
use crate::bc::ThermalBoundary;
use crate::mesh::ConductionMesh;

/// Canonical `fs-matdb` property consumed by the contact operator.
pub const AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY: &str =
    "area-specific-thermal-contact-resistance";

/// SI dimensions of `R''`, m² K/W.
pub const AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS: fs_qty::Dims = fs_qty::Dims([0, -1, 3, 1, 0, 0]);

/// One evidence-bearing area-specific interface resistance.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceResistance {
    value_m2_k_per_w: f64,
    uncertainty: UncertaintyModel,
    card_identity: ContentHash,
    receipt: PropertyUsageReceipt,
}

impl InterfaceResistance {
    /// Resolve one ordered interface card under an explicit query point and
    /// selection policy.  The query receipt and card identity travel with the
    /// resulting operator input.
    ///
    /// # Errors
    /// A typed interface refusal for a missing, ambiguous, out-of-domain, or
    /// non-positive claim; [`ConductionError::Dimensions`] for a claim that is
    /// not in m² K/W.
    pub fn from_card(
        interface_name: &str,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        if interface_name.trim().is_empty() {
            return Err(interface_error(
                "<unnamed>",
                "interface name is blank",
                "name the scenario interface so diagnostics and receipts can identify it",
            ));
        }
        let answer = card
            .claims()
            .query(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, point, policy)
            .map_err(|error| {
                interface_error(
                    interface_name,
                    format!(
                        "fs-matdb refused {AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY:?}: {error}"
                    ),
                    "attach exactly one in-domain thermal-resistance claim to the ordered interface card",
                )
            })?;
        let sample = &answer.evidence.value;
        if sample.dims != AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS {
            return Err(ConductionError::Dimensions {
                context: format!(
                    "interface {interface_name:?} property {AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY:?}"
                ),
                expected: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS.0,
                found: sample.dims.0,
            });
        }
        if !(sample.value.is_finite() && sample.value > 0.0) {
            return Err(interface_error(
                interface_name,
                format!(
                    "area-specific resistance {} m^2 K/W must be finite and strictly positive",
                    sample.value
                ),
                "supply a positive measured or model-card resistance; perfect contact is not an implicit default",
            ));
        }
        Ok(Self {
            value_m2_k_per_w: sample.value,
            uncertainty: sample.uncertainty.clone(),
            card_identity: card.content_hash(),
            receipt: answer.receipt,
        })
    }

    /// Area-specific resistance, m² K/W.
    #[must_use]
    pub const fn value_m2_k_per_w(&self) -> f64 {
        self.value_m2_k_per_w
    }

    /// Stated source uncertainty. [`UncertaintyModel::Unstated`] remains an
    /// explicit unknown, never a zero-width band.
    #[must_use]
    pub const fn uncertainty(&self) -> &UncertaintyModel {
        &self.uncertainty
    }

    /// Ordered interface-system-card identity.
    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    /// Receipt for the exact selected property claim.
    #[must_use]
    pub const fn receipt(&self) -> &PropertyUsageReceipt {
        &self.receipt
    }

    /// Convert the area-specific claim into one K/W series term.
    ///
    /// # Errors
    /// A typed interface refusal when `area_m2` is not finite and positive.
    pub fn term_for_area(
        &self,
        name: impl Into<String>,
        area_m2: f64,
    ) -> Result<ThermalResistanceTerm, ConductionError> {
        let name = name.into();
        if !(area_m2.is_finite() && area_m2 > 0.0) {
            return Err(interface_error(
                &name,
                format!("interface area {area_m2} m^2 must be finite and positive"),
                "supply the physical contact area used by both the operator and resistance budget",
            ));
        }
        let value = self.value_m2_k_per_w / area_m2;
        let uncertainty = match self.uncertainty {
            UncertaintyModel::Unstated => ResistanceUncertainty::Unstated,
            UncertaintyModel::HalfWidth {
                half_width,
                confidence,
            } => ResistanceUncertainty::HalfWidth {
                half_width: half_width / area_m2,
                confidence,
            },
            UncertaintyModel::RelativeHalfWidth {
                fraction,
                confidence,
            } => ResistanceUncertainty::HalfWidth {
                half_width: fraction * value.abs(),
                confidence,
            },
        };
        ThermalResistanceTerm::new(
            name,
            value,
            uncertainty,
            ResistanceOrigin::InterfaceCard {
                card_identity: self.card_identity,
                receipt: self.receipt.clone(),
            },
        )
    }
}

/// Uncertainty on one total-resistance term in K/W.
#[derive(Debug, Clone, PartialEq)]
pub enum ResistanceUncertainty {
    /// No quantitative uncertainty was supplied.  This is an unbounded term,
    /// not a zero-width interval.
    Unstated,
    /// Symmetric half-width at a stated confidence.
    HalfWidth {
        /// Half-width in K/W.
        half_width: f64,
        /// Confidence strictly inside `(0, 1)`.
        confidence: f64,
    },
}

impl ResistanceUncertainty {
    fn validate(&self, name: &str) -> Result<(), ConductionError> {
        match self {
            Self::Unstated => Ok(()),
            Self::HalfWidth {
                half_width,
                confidence,
            } if half_width.is_finite()
                && *half_width >= 0.0
                && confidence.is_finite()
                && *confidence > 0.0
                && *confidence < 1.0 =>
            {
                Ok(())
            }
            Self::HalfWidth {
                half_width,
                confidence,
            } => Err(interface_error(
                name,
                format!(
                    "resistance half-width {half_width} K/W and confidence {confidence} are inadmissible"
                ),
                "use a finite non-negative half-width and confidence strictly inside (0, 1)",
            )),
        }
    }
}

/// Where one K/W term came from.
#[derive(Debug, Clone, PartialEq)]
pub enum ResistanceOrigin {
    /// A caller declaration with a nonblank rationale.
    Declared {
        /// Why this term is authoritative for the current calculation.
        rationale: String,
    },
    /// An ordered interface card plus its selected-property receipt.
    InterfaceCard {
        /// Identity of the ordered surfaces, environment, and history.
        card_identity: ContentHash,
        /// Exact material-property selection receipt.
        receipt: PropertyUsageReceipt,
    },
}

/// One named K/W term in a series thermal-resistance network.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalResistanceTerm {
    name: String,
    value_k_per_w: f64,
    uncertainty: ResistanceUncertainty,
    origin: ResistanceOrigin,
}

impl ThermalResistanceTerm {
    /// Build a caller-declared term.  The declaration is allowed for analytic
    /// fixtures and explicit engineering inputs, but it never masquerades as
    /// material-card provenance.
    pub fn declared(
        name: impl Into<String>,
        value_k_per_w: f64,
        uncertainty: ResistanceUncertainty,
        rationale: impl Into<String>,
    ) -> Result<Self, ConductionError> {
        let rationale = rationale.into();
        if rationale.trim().is_empty() {
            return Err(interface_error(
                "<declared-term>",
                "declared resistance rationale is blank",
                "record why the declared value applies",
            ));
        }
        Self::new(
            name,
            value_k_per_w,
            uncertainty,
            ResistanceOrigin::Declared { rationale },
        )
    }

    /// Create the exact geometric contribution `L/(k A)` as a declared term.
    pub fn slab(
        name: impl Into<String>,
        length_m: f64,
        conductivity_w_per_mk: f64,
        area_m2: f64,
        uncertainty: ResistanceUncertainty,
        rationale: impl Into<String>,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        for (label, value) in [
            ("length", length_m),
            ("conductivity", conductivity_w_per_mk),
            ("area", area_m2),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(interface_error(
                    &name,
                    format!("slab {label} {value} must be finite and positive"),
                    "supply positive SI geometry and conductivity values",
                ));
            }
        }
        Self::declared(
            name,
            length_m / (conductivity_w_per_mk * area_m2),
            uncertainty,
            rationale,
        )
    }

    fn new(
        name: impl Into<String>,
        value_k_per_w: f64,
        uncertainty: ResistanceUncertainty,
        origin: ResistanceOrigin,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(interface_error(
                "<unnamed>",
                "series-resistance term name is blank",
                "name every layer and interface so the budget is actionable",
            ));
        }
        if !(value_k_per_w.is_finite() && value_k_per_w > 0.0) {
            return Err(interface_error(
                &name,
                format!("thermal resistance {value_k_per_w} K/W must be finite and positive"),
                "supply a positive resistance; zero contact resistance is never inferred",
            ));
        }
        uncertainty.validate(&name)?;
        Ok(Self {
            name,
            value_k_per_w,
            uncertainty,
            origin,
        })
    }

    /// Stable term name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Resistance in K/W.
    #[must_use]
    pub const fn value_k_per_w(&self) -> f64 {
        self.value_k_per_w
    }

    /// Quantitative or explicitly unstated uncertainty.
    #[must_use]
    pub const fn uncertainty(&self) -> &ResistanceUncertainty {
        &self.uncertainty
    }

    /// Provenance class and retained receipt, when present.
    #[must_use]
    pub const fn origin(&self) -> &ResistanceOrigin {
        &self.origin
    }
}

/// Conservative uncertainty summary for a series sum.
#[derive(Debug, Clone, PartialEq)]
pub struct SeriesResistanceBudget {
    /// Sum of all term values, K/W.
    pub value_k_per_w: f64,
    /// Sum of every stated symmetric half-width, K/W.
    pub stated_half_width_k_per_w: f64,
    /// Lowest confidence among stated bands, or `None` when none are stated.
    pub confidence_floor: Option<f64>,
    /// Terms whose uncertainty is unbounded/unstated.
    pub unbounded_terms: Vec<String>,
}

impl SeriesResistanceBudget {
    /// Complete conservative half-width, or `None` while any term is
    /// unbounded.  Unknown uncertainty never collapses to zero.
    #[must_use]
    pub fn complete_half_width_k_per_w(&self) -> Option<f64> {
        self.unbounded_terms
            .is_empty()
            .then_some(self.stated_half_width_k_per_w)
    }
}

/// Deterministic series network.  Terms are sorted by name before summation,
/// so input permutation cannot move a floating-point bit.
#[derive(Debug, Clone, PartialEq)]
pub struct SeriesThermalResistance {
    terms: Vec<ThermalResistanceTerm>,
    budget: SeriesResistanceBudget,
}

impl SeriesThermalResistance {
    /// Build a nonempty, uniquely named series network.
    pub fn new(mut terms: Vec<ThermalResistanceTerm>) -> Result<Self, ConductionError> {
        if terms.is_empty() {
            return Err(interface_error(
                "<series-network>",
                "thermal-resistance network has no terms",
                "declare every solid layer and interface in the heat path",
            ));
        }
        terms.sort_by(|a, b| a.name.cmp(&b.name));
        for pair in terms.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(interface_error(
                    &pair[0].name,
                    "series-resistance term name is duplicated",
                    "give each physical layer/interface a unique stable name",
                ));
            }
        }
        let mut value = 0.0f64;
        let mut stated_half_width = 0.0f64;
        let mut confidence_floor: Option<f64> = None;
        let mut unbounded_terms = Vec::new();
        for term in &terms {
            value += term.value_k_per_w;
            if !value.is_finite() {
                return Err(interface_error(
                    &term.name,
                    "series resistance overflowed the finite range",
                    "rescale or correct the resistance inputs",
                ));
            }
            match term.uncertainty {
                ResistanceUncertainty::Unstated => unbounded_terms.push(term.name.clone()),
                ResistanceUncertainty::HalfWidth {
                    half_width,
                    confidence,
                } => {
                    stated_half_width += half_width;
                    if !stated_half_width.is_finite() {
                        return Err(interface_error(
                            &term.name,
                            "series uncertainty half-width overflowed the finite range",
                            "rescale or correct the resistance uncertainty inputs",
                        ));
                    }
                    confidence_floor = Some(
                        confidence_floor.map_or(confidence, |current| current.min(confidence)),
                    );
                }
            }
        }
        Ok(Self {
            terms,
            budget: SeriesResistanceBudget {
                value_k_per_w: value,
                stated_half_width_k_per_w: stated_half_width,
                confidence_floor,
                unbounded_terms,
            },
        })
    }

    /// Canonically sorted terms.
    #[must_use]
    pub fn terms(&self) -> &[ThermalResistanceTerm] {
        &self.terms
    }

    /// Conservative sum and uncertainty completeness.
    #[must_use]
    pub const fn budget(&self) -> &SeriesResistanceBudget {
        &self.budget
    }
}

/// One oriented pair of coincident boundary-face slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct InterfaceFacePair {
    /// Side A boundary-face slot.
    pub side_a: usize,
    /// Side B boundary-face slot.
    pub side_b: usize,
}

impl InterfaceFacePair {
    fn normalized(self) -> (usize, usize) {
        if self.side_a < self.side_b {
            (self.side_a, self.side_b)
        } else {
            (self.side_b, self.side_a)
        }
    }
}

/// One named interface surface, possibly tessellated into many face pairs.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceSurface {
    name: String,
    face_pairs: Vec<InterfaceFacePair>,
    resistance: InterfaceResistance,
}

impl InterfaceSurface {
    /// Bind a named ordered card result to explicit oriented face pairs.
    pub fn new(
        name: impl Into<String>,
        face_pairs: Vec<InterfaceFacePair>,
        resistance: InterfaceResistance,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(interface_error(
                "<unnamed>",
                "interface surface name is blank",
                "use the scenario interface name",
            ));
        }
        if face_pairs.is_empty() {
            return Err(interface_error(
                &name,
                "interface surface has no paired faces",
                "bind at least one coincident face pair",
            ));
        }
        Ok(Self {
            name,
            face_pairs,
            resistance,
        })
    }

    /// Stable interface name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Ordered face pairs.
    #[must_use]
    pub fn face_pairs(&self) -> &[InterfaceFacePair] {
        &self.face_pairs
    }

    /// Evidence-bearing resistance.
    #[must_use]
    pub const fn resistance(&self) -> &InterfaceResistance {
        &self.resistance
    }
}

#[derive(Debug, Clone, PartialEq)]
struct BoundFacePair {
    pair: InterfaceFacePair,
    /// Side-A vertices and corresponding side-B vertices in the same order.
    side_a_vertices: [usize; 3],
    side_b_vertices: [usize; 3],
    area_m2: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct BoundSurface {
    name: String,
    faces: Vec<BoundFacePair>,
    resistance: InterfaceResistance,
}

/// Complete, fail-closed contact-interface set for one mesh/boundary pair.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalInterfaces {
    surfaces: Vec<BoundSurface>,
}

impl ThermalInterfaces {
    /// Validate and bind every coincident boundary-face pair.
    ///
    /// Geometry matching is intentionally exact in coordinate-bit space.  It
    /// supports matching P1 traces only; nonmatching/mortar contact requires a
    /// separately verified projection and is refused here.
    pub fn new(
        mesh: &ConductionMesh,
        boundary: &ThermalBoundary,
        mut surfaces: Vec<InterfaceSurface>,
    ) -> Result<Self, ConductionError> {
        let candidates = coincident_candidates(mesh)?;
        let candidate_keys = candidates
            .iter()
            .map(|pair| pair.normalized())
            .collect::<BTreeSet<_>>();
        surfaces.sort_by(|a, b| a.name.cmp(&b.name));
        for names in surfaces.windows(2) {
            if names[0].name == names[1].name {
                return Err(interface_error(
                    &names[0].name,
                    "interface surface name is duplicated",
                    "bind each named scenario interface exactly once",
                ));
            }
        }

        let mut claimed = BTreeSet::new();
        let mut bound_surfaces = Vec::with_capacity(surfaces.len());
        for mut surface in surfaces {
            surface.face_pairs.sort_by_key(|pair| pair.normalized());
            let mut bound_faces = Vec::with_capacity(surface.face_pairs.len());
            for pair in surface.face_pairs {
                if pair.side_a >= mesh.boundary().len() || pair.side_b >= mesh.boundary().len() {
                    return Err(interface_error(
                        &surface.name,
                        format!(
                            "face pair ({}, {}) is outside the {}-face boundary",
                            pair.side_a,
                            pair.side_b,
                            mesh.boundary().len()
                        ),
                        "use boundary-face slots from this ConductionMesh",
                    ));
                }
                if pair.side_a == pair.side_b {
                    return Err(interface_error(
                        &surface.name,
                        format!("face slot {} is paired with itself", pair.side_a),
                        "pair one face from each duplicated solid trace",
                    ));
                }
                let key = pair.normalized();
                if !candidate_keys.contains(&key) {
                    return Err(interface_error(
                        &surface.name,
                        format!(
                            "face pair ({}, {}) is not exactly coincident matching P1 geometry",
                            pair.side_a, pair.side_b
                        ),
                        "use exact duplicated interface coordinates or a future nonmatching mortar path",
                    ));
                }
                if !claimed.insert(key) {
                    return Err(interface_error(
                        &surface.name,
                        format!("face pair ({}, {}) was bound more than once", key.0, key.1),
                        "assign every coincident pair to exactly one interface card",
                    ));
                }
                for slot in [pair.side_a, pair.side_b] {
                    if boundary.condition_for(slot).is_some() {
                        return Err(interface_error(
                            &surface.name,
                            format!(
                                "interface face slot {slot} also carries an external boundary condition"
                            ),
                            "leave paired faces in the explicit adiabatic remainder; the interface operator owns their transfer",
                        ));
                    }
                }
                let b_for_a = matching_vertices(mesh, pair)?;
                let side_a_vertices = mesh.boundary()[pair.side_a]
                    .vertices
                    .map(|vertex| vertex as usize);
                let side_b_vertices =
                    b_for_a.map(|local| mesh.boundary()[pair.side_b].vertices[local] as usize);
                bound_faces.push(BoundFacePair {
                    pair,
                    side_a_vertices,
                    side_b_vertices,
                    area_m2: mesh.boundary()[pair.side_a].area,
                });
            }
            bound_surfaces.push(BoundSurface {
                name: surface.name,
                faces: bound_faces,
                resistance: surface.resistance,
            });
        }

        if let Some(&(a, b)) = candidate_keys.difference(&claimed).next() {
            return Err(interface_error(
                "<undeclared>",
                format!("coincident boundary faces {a} and {b} have no interface card"),
                "declare the named interface and bind an in-domain fs-matdb resistance; zero resistance is never assumed",
            ));
        }
        Ok(Self {
            surfaces: bound_surfaces,
        })
    }

    /// Coincident matching-face candidates in deterministic boundary-slot
    /// order.  Slot order is diagnostic only; callers must assign semantic A/B
    /// orientation when constructing [`InterfaceFacePair`].
    pub fn coincident_face_pairs(
        mesh: &ConductionMesh,
    ) -> Result<Vec<InterfaceFacePair>, ConductionError> {
        coincident_candidates(mesh)
    }

    pub(crate) fn require_no_undeclared(mesh: &ConductionMesh) -> Result<(), ConductionError> {
        if let Some(pair) = coincident_candidates(mesh)?.first() {
            return Err(interface_error(
                "<undeclared>",
                format!(
                    "coincident boundary faces {} and {} require an interface card",
                    pair.side_a, pair.side_b
                ),
                "construct ThermalInterfaces and pass it to assembly/solve; zero resistance is never assumed",
            ));
        }
        Ok(())
    }

    pub(crate) fn validate_for(
        &self,
        mesh: &ConductionMesh,
        boundary: &ThermalBoundary,
    ) -> Result<(), ConductionError> {
        let candidates = coincident_candidates(mesh)?
            .into_iter()
            .map(InterfaceFacePair::normalized)
            .collect::<BTreeSet<_>>();
        let bound = self
            .surfaces
            .iter()
            .flat_map(|surface| surface.faces.iter())
            .map(|face| face.pair.normalized())
            .collect::<BTreeSet<_>>();
        if candidates != bound {
            return Err(interface_error(
                "<binding>",
                "contact-interface face set does not match the supplied mesh",
                "rebuild ThermalInterfaces from the exact mesh and boundary used by this solve",
            ));
        }
        for surface in &self.surfaces {
            for face in &surface.faces {
                for slot in [face.pair.side_a, face.pair.side_b] {
                    if boundary.condition_for(slot).is_some() {
                        return Err(interface_error(
                            &surface.name,
                            format!(
                                "interface face slot {slot} now carries an external boundary condition"
                            ),
                            "rebuild the boundary with paired faces in its explicit adiabatic remainder",
                        ));
                    }
                }
                let side_a_vertices = mesh.boundary()[face.pair.side_a]
                    .vertices
                    .map(|vertex| vertex as usize);
                let mapping = matching_vertices(mesh, face.pair)?;
                let side_b_vertices =
                    mapping.map(|local| mesh.boundary()[face.pair.side_b].vertices[local] as usize);
                if side_a_vertices != face.side_a_vertices
                    || side_b_vertices != face.side_b_vertices
                    || mesh.boundary()[face.pair.side_a].area.to_bits() != face.area_m2.to_bits()
                {
                    return Err(interface_error(
                        &surface.name,
                        format!(
                            "contact-interface face pair ({}, {}) changed after binding",
                            face.pair.side_a, face.pair.side_b
                        ),
                        "rebuild ThermalInterfaces after changing mesh geometry or topology",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Bound interface-surface count (one receipt per surface).
    #[must_use]
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// Aggregate one flux record per named interface.
    pub fn fluxes(&self, temperature: &[f64]) -> Result<Vec<InterfaceFlux>, ConductionError> {
        let maximum_vertex = self
            .surfaces
            .iter()
            .flat_map(|surface| surface.faces.iter())
            .flat_map(|face| face.side_a_vertices.into_iter().chain(face.side_b_vertices))
            .max()
            .unwrap_or(0);
        if temperature.len() <= maximum_vertex {
            return Err(ConductionError::FieldLength {
                field: "temperature for interface flux",
                expected: maximum_vertex + 1,
                found: temperature.len(),
            });
        }
        let mut out = Vec::with_capacity(self.surfaces.len());
        for surface in &self.surfaces {
            let conductance_density = 1.0 / surface.resistance.value_m2_k_per_w;
            let mut area_m2 = 0.0f64;
            let mut heat_rate = 0.0f64;
            let mut conductance = 0.0f64;
            for face in &surface.faces {
                let mean_jump = face
                    .side_a_vertices
                    .iter()
                    .zip(face.side_b_vertices)
                    .map(|(&a, b)| temperature[a] - temperature[b])
                    .sum::<f64>()
                    / 3.0;
                let face_conductance = face.area_m2 * conductance_density;
                area_m2 += face.area_m2;
                conductance += face_conductance;
                heat_rate = face_conductance.mul_add(mean_jump, heat_rate);
            }
            out.push(InterfaceFlux {
                interface: surface.name.clone(),
                area_m2,
                conductance_w_per_k: conductance,
                mean_jump_k: heat_rate / conductance,
                heat_rate_a_to_b_w: heat_rate,
                card_identity: surface.resistance.card_identity,
                receipt: surface.resistance.receipt.clone(),
            });
        }
        Ok(out)
    }

    pub(crate) fn assemble_into(&self, cx: &Cx<'_>, coo: &mut Coo) -> Result<(), ConductionError> {
        let mut face_index = 0usize;
        for surface in &self.surfaces {
            let conductance_density = 1.0 / surface.resistance.value_m2_k_per_w;
            for face in &surface.faces {
                if face_index % crate::assemble::ASSEMBLY_TILE == 0 {
                    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                        stage: "assemble-interfaces",
                        at: face_index,
                    })?;
                }
                for a in 0..3 {
                    let va = face.side_a_vertices[a];
                    let wa = face.side_b_vertices[a];
                    for b in 0..3 {
                        let vb = face.side_a_vertices[b];
                        let wb = face.side_b_vertices[b];
                        let mass = if a == b {
                            face.area_m2 / 6.0
                        } else {
                            face.area_m2 / 12.0
                        };
                        let value = conductance_density * mass;
                        coo.push(va, vb, value);
                        coo.push(wa, wb, value);
                        coo.push(va, wb, -value);
                        coo.push(wa, vb, -value);
                    }
                }
                face_index += 1;
            }
        }
        Ok(())
    }
}

/// Integrated heat transfer across one named interface.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceFlux {
    /// Stable scenario/interface name.
    pub interface: String,
    /// Total paired area, m².
    pub area_m2: f64,
    /// `area/R''`, W/K.
    pub conductance_w_per_k: f64,
    /// Conductance-weighted `T_a - T_b`, K.
    pub mean_jump_k: f64,
    /// Positive from declared side A to side B, W.
    pub heat_rate_a_to_b_w: f64,
    /// Ordered interface-card identity.
    pub card_identity: ContentHash,
    /// Property selection receipt.
    pub receipt: PropertyUsageReceipt,
}

type PointKey = [u64; 3];
type FaceKey = [PointKey; 3];

fn point_key(point: [f64; 3]) -> PointKey {
    [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()]
}

fn face_key(mesh: &ConductionMesh, slot: usize) -> FaceKey {
    let mut key = mesh.boundary()[slot]
        .vertices
        .map(|vertex| point_key(mesh.positions()[vertex as usize]));
    key.sort_unstable();
    key
}

fn coincident_candidates(mesh: &ConductionMesh) -> Result<Vec<InterfaceFacePair>, ConductionError> {
    let mut groups = BTreeMap::<FaceKey, Vec<usize>>::new();
    for slot in 0..mesh.boundary().len() {
        groups.entry(face_key(mesh, slot)).or_default().push(slot);
    }
    let mut pairs = Vec::new();
    for slots in groups.values() {
        match slots.as_slice() {
            [_] => {}
            [a, b] => {
                let normal_a = mesh.boundary()[*a].outward_normal;
                let normal_b = mesh.boundary()[*b].outward_normal;
                let dot = normal_a[0].mul_add(
                    normal_b[0],
                    normal_a[1].mul_add(normal_b[1], normal_a[2] * normal_b[2]),
                );
                if dot > -1.0 + 64.0 * f64::EPSILON {
                    return Err(interface_error(
                        "<geometry>",
                        format!(
                            "coincident faces {a} and {b} do not have opposing normals (dot={dot})"
                        ),
                        "repair the duplicated solid orientation before applying contact transfer",
                    ));
                }
                pairs.push(InterfaceFacePair {
                    side_a: *a,
                    side_b: *b,
                });
            }
            many => {
                return Err(interface_error(
                    "<geometry>",
                    format!(
                        "{} boundary faces share one exact geometric triangle: {many:?}",
                        many.len()
                    ),
                    "supply a two-sided manifold interface or a separately modeled junction",
                ));
            }
        }
    }
    pairs.sort_by_key(|pair| pair.normalized());
    Ok(pairs)
}

fn matching_vertices(
    mesh: &ConductionMesh,
    pair: InterfaceFacePair,
) -> Result<[usize; 3], ConductionError> {
    let a = &mesh.boundary()[pair.side_a];
    let b = &mesh.boundary()[pair.side_b];
    let mut mapping = [usize::MAX; 3];
    for (a_local, &a_vertex) in a.vertices.iter().enumerate() {
        let key = point_key(mesh.positions()[a_vertex as usize]);
        let matches = b
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, vertex)| point_key(mesh.positions()[**vertex as usize]) == key)
            .map(|(local, _)| local)
            .collect::<Vec<_>>();
        if matches.len() != 1 {
            return Err(interface_error(
                "<geometry>",
                format!(
                    "face pair ({}, {}) does not have a unique exact vertex correspondence",
                    pair.side_a, pair.side_b
                ),
                "use matching P1 traces or a future nonmatching mortar path",
            ));
        }
        mapping[a_local] = matches[0];
    }
    Ok(mapping)
}

fn interface_error(
    interface: impl Into<String>,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> ConductionError {
    ConductionError::Interface {
        interface: interface.into(),
        what: what.into(),
        fix: fix.into(),
    }
}
