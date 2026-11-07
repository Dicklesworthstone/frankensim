//! Material-driven truss layout, catalog sizing, mass and rod re-analysis.
//!
//! The same immutable state-point bundle supplies density, Young's modulus
//! and signed axial strength limits to every stage. Assignment is per member; different cards
//! and temperatures require data changes, never a different equilibrium solver.
//! This is the existing axial, solid-square, pin-jointed equilibrium model. The caller
//! explicitly selects either symmetric yield or separate tension/compression
//! limits; the same choice reaches layout, refitting and catalog sizing.
//! Euler sizing defaults to pinned effective lengths; separately admitted
//! restraint/bracing factors persist through sizing and self-weight rebuilds.
//! It does not infer a code allowable from an ultimate strength, model creep or
//! thermal strain, or add self-weight implicitly to the declared load case.

/// Material mass, nodal body forces and consistent catalog self-weight design.
pub mod body_load;

use fs_exec::Cx;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{Density, Pressure};

use crate::ground::{GroundStructure, TrussConstructionError};
use crate::lp::{LayoutCase, LayoutLimits, LayoutLp, PdhgError, PdhgReport, PdhgSettings};
use crate::sizing::{CatalogAudit, size_and_snap_members};

fn invalid(field: &'static str, requirement: &'static str) -> TrussConstructionError {
    TrussConstructionError::InvalidInput { field, requirement }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TrussConstructionError> {
    cx.checkpoint()
        .map_err(|_| TrussConstructionError::Cancelled {
            stage: "material-driven truss",
        })
}

/// Exactly the material properties used by axial layout and square-member sizing.
/// No elastic, thermal or strength coefficient is inferred from a material name.
#[derive(Clone, Debug, PartialEq)]
pub struct TrussMaterial {
    resolved: ResolvedMaterialStatePoint,
    density: f64,
    youngs: f64,
    tensile_limit: f64,
    compressive_limit: f64,
    tension_property: String,
    compression_property: String,
}

impl TrussMaterial {
    /// Resolve density, Young's modulus and symmetric uniaxial yield stress at
    /// one physical state. The caller explicitly admits that yield model for
    /// both tension and compression; brittle/anisotropic failure is not implied.
    ///
    /// # Errors
    /// Missing, ambiguous, unsupported, mistyped or nonpositive properties refuse
    /// atomically through the shared material resolver. Poisson ratio and thermal
    /// data are not required because this axial sizing model does not use them.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        Self::resolve_axial_limits(
            card,
            point,
            YIELD_STRESS_PROPERTY,
            YIELD_STRESS_PROPERTY,
            selection,
        )
    }

    /// Resolve density, axial modulus and explicitly named tensile/compressive
    /// stress limits at one physical state. Limits need not be equal, and a
    /// property used for both signs is resolved only once.
    ///
    /// These are caller-admitted uniaxial limits, not an inferred yield surface,
    /// allowable stress, safety factor or post-failure law. A directional claim
    /// is appropriate only for a member aligned with that measured direction;
    /// this scalar model does not rotate anisotropic strengths.
    ///
    /// # Errors
    /// Missing, ambiguous, mistyped, nonpositive or dynamic properties refuse
    /// through the shared resolver. No absent compression limit is copied from
    /// tension, and no ultimate strength is silently substituted for yield.
    pub fn resolve_axial_limits(
        card: &MaterialCard,
        point: &QueryPoint,
        tension_property: &str,
        compression_property: &str,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        let mut requirements = vec![
            ScalarPropertyRequirement::try_new(
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            ScalarPropertyRequirement::try_new(
                tension_property,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
        ];
        if compression_property != tension_property {
            requirements.push(ScalarPropertyRequirement::try_new(
                compression_property,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
        let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
        for requirement in &requirements {
            let name = requirement.name();
            if !resolved
                .property(name)
                .expect("required property was resolved")
                .is_static_at_fixed_external_conditions()
            {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: name.into(),
                    reason: "static axial sizing cannot freeze frequency, strain, stress or rate dependence",
                });
            }
        }
        let value = |name| {
            resolved
                .property(name)
                .expect("required property was resolved")
                .value_si()
        };
        Ok(Self {
            density: value(DENSITY_PROPERTY),
            youngs: value(YOUNG_MODULUS_PROPERTY),
            tensile_limit: value(tension_property),
            compressive_limit: value(compression_property),
            tension_property: tension_property.into(),
            compression_property: compression_property.into(),
            resolved,
        })
    }

    /// Original material identity, state coordinates and selected claim receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Current material density [kg/m³].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density
    }

    /// Current axial Young's modulus [Pa].
    #[must_use]
    pub const fn young_modulus_pa(&self) -> f64 {
        self.youngs
    }

    /// Current tensile stress limit [Pa], without an implicit safety factor.
    #[must_use]
    pub const fn tensile_limit_pa(&self) -> f64 {
        self.tensile_limit
    }

    /// Current compressive stress magnitude limit [Pa]. Euler buckling is a
    /// separate limit imposed during catalog sizing.
    #[must_use]
    pub const fn compressive_limit_pa(&self) -> f64 {
        self.compressive_limit
    }

    /// Exact sourced property selected for positive axial force.
    #[must_use]
    pub fn tension_property(&self) -> &str {
        &self.tension_property
    }

    /// Exact sourced property selected for negative axial force.
    #[must_use]
    pub fn compression_property(&self) -> &str {
        &self.compression_property
    }
}

/// Physical quantity minimized by the axial layout before Euler/catalog sizing.
/// For unlike source densities, minimum volume and minimum mass generally select
/// different load paths. This does not optimize the subsequent discrete catalog.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MaterialLayoutObjective {
    /// Minimize `sum length * (q_tension / limit_tension + q_compression / limit_compression)` [m³].
    #[default]
    Volume,
    /// Minimize that same per-member volume multiplied by its source density [kg].
    Mass,
}

/// A material assignment and its admitted layout operator. Borrowing immutable
/// geometry and source states prevents later sizing from accidentally using a
/// different member order, temperature, modulus or signed axial limits.
pub struct MaterialLayout<'a> {
    ground: &'a GroundStructure,
    materials: &'a [TrussMaterial],
    assignments: &'a [usize],
    case: LayoutCase,
    limits: LayoutLimits,
    lp: LayoutLp,
    strength_scale_pa: f64,
    density_scale_kg_m3: f64,
    objective: MaterialLayoutObjective,
    effective_buckling_length_factors: Option<&'a [f64]>,
}

/// Layout iterates and their physical volume and mass. The numerical report uses the
/// normalized objective; its certificates, if separately requested, also do.
#[derive(Debug)]
pub struct MaterialLayoutSolution {
    /// Split member forces `[tension | compression]` [N].
    pub forces: Vec<f64>,
    /// Dual iterate in the normalized objective convention.
    pub dual: Vec<f64>,
    /// Existing PDHG diagnostics. Its legacy `volume` field is the normalized
    /// selected objective: `strength_scale_pa * volume_m3` for Volume, or
    /// `strength_scale_pa * mass_kg / density_scale_kg_m3` for Mass.
    pub numerical: PdhgReport,
    /// Approximate returned-iterate material volume [m³], before catalog sizing.
    pub volume_m3: f64,
    /// Approximate returned-iterate mass [kg], before Euler/catalog sizing.
    pub mass_kg: f64,
    /// Physical objective actually used to select the returned force path.
    pub objective: MaterialLayoutObjective,
}

impl MaterialLayoutSolution {
    /// Selected physical objective [m³ for Volume, kg for Mass].
    #[must_use]
    pub const fn objective_value_si(&self) -> f64 {
        match self.objective {
            MaterialLayoutObjective::Volume => self.volume_m3,
            MaterialLayoutObjective::Mass => self.mass_kg,
        }
    }
}

/// Failure of the common solve or its conversion back to physical material totals.
#[derive(Debug)]
pub enum MaterialLayoutSolveError {
    /// The underlying numerical owner refused its inputs or iterates.
    Solver(PdhgError),
    /// Positive normalized volume overflowed or underflowed in SI units.
    UnrepresentableVolume,
    /// Source-density mass overflowed or underflowed in SI units.
    UnrepresentableMass,
}
impl core::fmt::Display for MaterialLayoutSolveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solver(error) => write!(f, "material layout: {error}"),
            Self::UnrepresentableVolume => {
                f.write_str("material layout volume is unrepresentable in m3")
            }
            Self::UnrepresentableMass => {
                f.write_str("material layout mass is unrepresentable in kg")
            }
        }
    }
}
impl std::error::Error for MaterialLayoutSolveError {}

/// The existing member audit plus physical material totals from those very
/// same selected states. Missing catalog sizes leave totals unavailable.
#[derive(Debug, Clone)]
pub struct MaterialCatalogAudit {
    /// Equilibrium, strength, Euler and catalog checks from the common solver.
    pub audit: CatalogAudit,
    /// Per-survivor mass [kg], in `audit.members` order; no size means `None`.
    pub member_mass_kg: Vec<Option<f64>>,
    /// Total catalog volume [m³], available only when every size is finite.
    pub volume_m3: Option<f64>,
    /// Total mass [kg], available only when every member has a supported size.
    pub mass_kg: Option<f64>,
}

impl<'a> MaterialLayout<'a> {
    /// Assemble the existing LP with one source-material index per member.
    /// Material rows can be reused by any number of members. Costs are
    /// `length * reference_strength / member_limit` separately for each sign:
    /// one common positive scaling
    /// preserves the heterogeneous optimum while avoiding tiny SI objective
    /// coefficients that stall the unscaled PDHG lane.
    ///
    /// # Errors
    /// Refuses incomplete/out-of-range assignments, unrepresentable normalized
    /// properties, ordinary LP admission failures or cancellation.
    pub fn try_new(
        ground: &'a GroundStructure,
        case: &LayoutCase,
        materials: &'a [TrussMaterial],
        assignments: &'a [usize],
        limits: LayoutLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, TrussConstructionError> {
        Self::try_new_with_objective(
            ground,
            case,
            materials,
            assignments,
            MaterialLayoutObjective::Volume,
            limits,
            cx,
        )
    }

    /// Assemble an explicitly selected volume or mass objective on the SAME
    /// equilibrium and signed-strength problem. For mass, each normalized
    /// volume cost is multiplied by `member_density / reference_density`.
    /// The common density scale preserves conditioning and does not change the
    /// optimum; density is never substituted for the actual sizing strength.
    /// A uniform density field therefore gives the original volume operator.
    ///
    /// # Errors
    /// Inherits assignment, source-ratio, resource and cancellation refusals.
    pub fn try_new_with_objective(
        ground: &'a GroundStructure,
        case: &LayoutCase,
        materials: &'a [TrussMaterial],
        assignments: &'a [usize],
        objective: MaterialLayoutObjective,
        limits: LayoutLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, TrussConstructionError> {
        checkpoint(cx)?;
        if assignments.len() != ground.members().len() {
            return Err(TrussConstructionError::VectorLength {
                field: "member material assignments",
                expected: ground.members().len(),
                actual: assignments.len(),
            });
        }
        if assignments.is_empty() || materials.is_empty() {
            return Err(invalid(
                "member materials",
                "must contain at least one assigned material",
            ));
        }
        for (member, &index) in assignments.iter().enumerate() {
            if member % 256 == 0 {
                checkpoint(cx)?;
            }
            if index >= materials.len() {
                return Err(invalid(
                    "member material index",
                    "must refer to an admitted material row",
                ));
            }
        }
        let strength_scale_pa = materials[assignments[0]].tensile_limit;
        let density_scale_kg_m3 = materials[assignments[0]].density;
        let lp = LayoutLp::assemble_member_limits(
            ground,
            case,
            |member| {
                let state = &materials[assignments[member]];
                (
                    state.tensile_limit / strength_scale_pa,
                    state.compressive_limit / strength_scale_pa,
                )
            },
            |member| match objective {
                MaterialLayoutObjective::Volume => 1.0,
                MaterialLayoutObjective::Mass => {
                    materials[assignments[member]].density / density_scale_kg_m3
                }
            },
            limits,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(Self {
            ground,
            materials,
            assignments,
            case: case.clone(),
            limits,
            lp,
            strength_scale_pa,
            density_scale_kg_m3,
            objective,
            effective_buckling_length_factors: None,
        })
    }

    /// Bind one dimensionless Euler effective-length factor `K` per ORIGINAL
    /// ground member. Compression sizing uses `P_cr = pi² E I / (K L)²`;
    /// omission retains pinned `K = 1`. The same borrowed factors reach every
    /// self-weight rebuild and the post-snap audit, including after pruning.
    ///
    /// The caller must establish the restraint/bracing model independently.
    /// This does not infer end fixity from nodal translation supports, change
    /// pin-jointed equilibrium or physical member mass, or establish global
    /// frame stability, inelastic buckling or postbuckling strength. The axial
    /// LP objective and its optional certificates remain unchanged.
    ///
    /// # Errors
    /// Refuses incomplete assignments, nonfinite/nonpositive factors, an
    /// unrepresentable effective length, or cancellation before publication.
    pub fn with_effective_buckling_length_factors(
        mut self,
        factors: &'a [f64],
        cx: &Cx<'_>,
    ) -> Result<Self, TrussConstructionError> {
        checkpoint(cx)?;
        if factors.len() != self.assignments.len() {
            return Err(TrussConstructionError::VectorLength {
                field: "effective buckling length factors",
                expected: self.assignments.len(),
                actual: factors.len(),
            });
        }
        for (member, &factor) in factors.iter().enumerate() {
            if member % 256 == 0 {
                checkpoint(cx)?;
            }
            let effective_length = factor * self.ground.lengths()[member];
            if !factor.is_finite()
                || factor <= 0.0
                || !effective_length.is_finite()
                || effective_length <= 0.0
            {
                return Err(invalid(
                    "effective buckling length factor",
                    "requires finite positive K and representable positive K * member length",
                ));
            }
        }
        checkpoint(cx)?;
        self.effective_buckling_length_factors = Some(factors);
        Ok(self)
    }

    /// Admitted Euler effective-length factor for an original member, or
    /// `None` if the member does not exist. The implicit pinned default is 1.
    #[must_use]
    pub fn effective_buckling_length_factor(&self, member: usize) -> Option<f64> {
        self.assignments.get(member)?;
        Some(
            self.effective_buckling_length_factors
                .map_or(1.0, |factors| factors[member]),
        )
    }

    /// Common equilibrium/PDHG/certificate API. Its objective is uniformly
    /// scaled: divide by [`Self::strength_scale_pa`] for Volume [m³]. For Mass,
    /// additionally multiply by [`Self::density_scale_kg_m3`] to recover kg.
    /// Certificates remain bounds on that explicitly normalized objective.
    #[must_use]
    pub const fn lp(&self) -> &LayoutLp {
        &self.lp
    }

    /// Original applied loads and supports, including loads at supported nodes.
    /// Body-load iterations always add to this immutable case once.
    #[must_use]
    pub const fn load_case(&self) -> &LayoutCase {
        &self.case
    }

    /// Shared positive objective normalization [Pa].
    #[must_use]
    pub const fn strength_scale_pa(&self) -> f64 {
        self.strength_scale_pa
    }

    /// Reference source density [kg/m³] used only to normalize a mass objective.
    #[must_use]
    pub const fn density_scale_kg_m3(&self) -> f64 {
        self.density_scale_kg_m3
    }

    /// Selected quantity minimized by the original LP and any body-load reloads.
    #[must_use]
    pub const fn objective(&self) -> MaterialLayoutObjective {
        self.objective
    }

    /// Original state used for an assigned member, or `None` for an absent member.
    #[must_use]
    pub fn material(&self, member: usize) -> Option<&TrussMaterial> {
        self.assignments
            .get(member)
            .map(|&index| &self.materials[index])
    }

    /// Solve through the existing PDHG owner, retaining both physical totals.
    /// This inherits that solver's iteration bound and synchronous cancellation
    /// limitation; assembly and material publication remain cancellation-polled.
    ///
    /// # Errors
    /// Propagates the common solver's input and numerical refusals.
    pub fn solve(
        &self,
        settings: PdhgSettings,
    ) -> Result<MaterialLayoutSolution, MaterialLayoutSolveError> {
        let (forces, dual, numerical) = self
            .lp
            .solve(None, None, settings)
            .map_err(MaterialLayoutSolveError::Solver)?;
        let mut volume_m3 = 0.0;
        let mut mass_kg = 0.0;
        let member_count = self.assignments.len();
        for member in 0..member_count {
            let state = &self.materials[self.assignments[member]];
            let volume = self.ground.lengths()[member]
                * (forces[member] / state.tensile_limit
                    + forces[member + member_count] / state.compressive_limit);
            volume_m3 += volume;
            if !volume_m3.is_finite()
                || volume < 0.0
                || (volume == 0.0 && (forces[member] > 0.0 || forces[member + member_count] > 0.0))
            {
                return Err(MaterialLayoutSolveError::UnrepresentableVolume);
            }
            let mass = volume * state.density;
            mass_kg += mass;
            if !mass_kg.is_finite() || (volume > 0.0 && mass == 0.0) {
                return Err(MaterialLayoutSolveError::UnrepresentableMass);
            }
        }
        Ok(MaterialLayoutSolution {
            forces,
            dual,
            numerical,
            volume_m3,
            mass_kg,
            objective: self.objective,
        })
    }

    /// Refit equilibrium and size every survivor with the resolved limit for
    /// its REFITTED force sign, its own modulus and its admitted effective
    /// buckling length. Mass comes from its source
    /// density and actual snapped geometry, so swapping a material propagates
    /// through all three stages.
    ///
    /// # Errors
    /// Invalid force shape/domain, catalog ordering, pruning controls, overflow
    /// or cancellation refuse without publishing partial material totals.
    pub fn size_and_snap(
        &self,
        forces: &[f64],
        catalog: &[f64],
        prune_frac: f64,
        cx: &Cx<'_>,
    ) -> Result<MaterialCatalogAudit, TrussConstructionError> {
        checkpoint(cx)?;
        if forces.len() != self.lp.c().len() {
            return Err(TrussConstructionError::VectorLength {
                field: "split member forces",
                expected: self.lp.c().len(),
                actual: forces.len(),
            });
        }
        if forces.iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err(invalid(
                "split member forces",
                "must be finite and nonnegative",
            ));
        }
        if catalog.is_empty()
            || catalog.iter().any(|a| !a.is_finite() || *a <= 0.0)
            || catalog.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "catalog",
                "must contain positive finite areas in strictly increasing order",
            ));
        }
        if !prune_frac.is_finite() || !(0.0..=1.0).contains(&prune_frac) {
            return Err(invalid("prune fraction", "must be finite and in [0,1]"));
        }
        let audit = size_and_snap_members(
            self.ground,
            &self.lp,
            forces,
            |member| {
                let state = &self.materials[self.assignments[member]];
                (
                    state.tensile_limit,
                    state.compressive_limit,
                    state.youngs,
                    self.effective_buckling_length_factors
                        .map_or(1.0, |factors| factors[member]),
                )
            },
            catalog,
            prune_frac,
            fs_constraint::ConstraintKind::Hard,
        );
        let mut member_mass_kg = Vec::with_capacity(audit.members.len());
        let mut volume = 0.0;
        let mut mass = 0.0;
        let mut complete = true;
        for (index, member) in audit.members.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if !member.area_catalog.is_finite() {
                complete = false;
                member_mass_kg.push(None);
                continue;
            }
            let member_volume = member.area_catalog * member.length;
            let member_mass =
                member_volume * self.materials[self.assignments[member.index]].density;
            volume += member_volume;
            mass += member_mass;
            if !volume.is_finite() || !mass.is_finite() || member_mass <= 0.0 {
                return Err(invalid(
                    "catalog mass/volume",
                    "must remain finite and positive",
                ));
            }
            member_mass_kg.push(Some(member_mass));
        }
        checkpoint(cx)?;
        Ok(MaterialCatalogAudit {
            audit,
            member_mass_kg,
            volume_m3: complete.then_some(volume),
            mass_kg: complete.then_some(mass),
        })
    }

    /// Run the existing conservative shear-stiff Cosserat spot check with the
    /// member's geometry and material modulus. This retains that model's fixed
    /// ten-segment clamped-free approximation, not a certified buckling load.
    ///
    /// # Errors
    /// Refuses an absent member or invalid loading/section controls.
    pub fn rod_buckling_check(
        &self,
        member: usize,
        area_m2: f64,
        design_force_n: f64,
        factor: f64,
        imperfection: f64,
    ) -> Result<(bool, f64), TrussConstructionError> {
        let state = self
            .material(member)
            .ok_or_else(|| invalid("member", "must identify an assigned member"))?;
        if [area_m2, factor, imperfection]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || !design_force_n.is_finite()
            || design_force_n == 0.0
        {
            return Err(invalid(
                "rod check",
                "requires finite positive area, factor, imperfection and nonzero force",
            ));
        }
        Ok(crate::rodcheck::rod_buckling_check(
            self.ground.lengths()[member],
            area_m2,
            state.youngs,
            design_force_n,
            factor,
            imperfection,
        ))
    }
}
