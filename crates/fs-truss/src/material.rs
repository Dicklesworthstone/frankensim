//! Material-driven truss layout, catalog sizing, mass and rod re-analysis.
//!
//! The same immutable state-point bundle supplies density, Young's modulus
//! and yield stress to every stage. Assignment is per member; different cards
//! and temperatures require data changes, never a different equilibrium solver.
//! This is the existing axial, symmetric-yield, solid-square, pin-jointed model.
//! It does not infer a code allowable from an ultimate strength, model creep or
//! thermal strain, or add self-weight implicitly to the declared load case.

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
    yield_stress: f64,
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
        let requirements = [
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
                YIELD_STRESS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
        ];
        let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
        for name in [
            DENSITY_PROPERTY,
            YOUNG_MODULUS_PROPERTY,
            YIELD_STRESS_PROPERTY,
        ] {
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
            yield_stress: value(YIELD_STRESS_PROPERTY),
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

    /// Current symmetric yield stress [Pa].
    #[must_use]
    pub const fn yield_stress_pa(&self) -> f64 {
        self.yield_stress
    }
}

/// A material assignment and its admitted layout operator. Borrowing immutable
/// geometry and source states prevents later sizing from accidentally using a
/// different member order, temperature, modulus or yield stress.
pub struct MaterialLayout<'a> {
    ground: &'a GroundStructure,
    materials: &'a [TrussMaterial],
    assignments: &'a [usize],
    lp: LayoutLp,
    yield_scale_pa: f64,
}

/// Layout iterates and their physical volume. The numerical report uses the
/// normalized objective; its certificates, if separately requested, also do.
pub struct MaterialLayoutSolution {
    /// Split member forces `[tension | compression]` [N].
    pub forces: Vec<f64>,
    /// Dual iterate in the normalized objective convention.
    pub dual: Vec<f64>,
    /// Existing PDHG diagnostics; `volume` is `yield_scale_pa * volume_m3`.
    pub numerical: PdhgReport,
    /// Approximate returned-iterate material volume [m³], before catalog sizing.
    pub volume_m3: f64,
}

/// Failure of the common solve or its conversion back to physical volume.
#[derive(Debug)]
pub enum MaterialLayoutSolveError {
    /// The underlying numerical owner refused its inputs or iterates.
    Solver(PdhgError),
    /// Positive normalized volume overflowed or underflowed in SI units.
    UnrepresentableVolume,
}
impl core::fmt::Display for MaterialLayoutSolveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solver(error) => write!(f, "material layout: {error}"),
            Self::UnrepresentableVolume => {
                f.write_str("material layout volume is unrepresentable in m3")
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
    /// `length * reference_yield / member_yield`: one common positive scaling
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
        let yield_scale_pa = materials[assignments[0]].yield_stress;
        let lp = LayoutLp::assemble_member_yields(
            ground,
            case,
            |member| materials[assignments[member]].yield_stress / yield_scale_pa,
            limits,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(Self {
            ground,
            materials,
            assignments,
            lp,
            yield_scale_pa,
        })
    }

    /// Common equilibrium/PDHG/certificate API. Its objective is uniformly
    /// scaled; divide objective values by [`Self::yield_scale_pa`] for m³.
    #[must_use]
    pub const fn lp(&self) -> &LayoutLp {
        &self.lp
    }

    /// Shared positive objective normalization [Pa].
    #[must_use]
    pub const fn yield_scale_pa(&self) -> f64 {
        self.yield_scale_pa
    }

    /// Original state used for an assigned member, or `None` for an absent member.
    #[must_use]
    pub fn material(&self, member: usize) -> Option<&TrussMaterial> {
        self.assignments
            .get(member)
            .map(|&index| &self.materials[index])
    }

    /// Solve through the existing PDHG owner, retaining physical volume units.
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
        let volume_m3 = numerical.volume / self.yield_scale_pa;
        if !volume_m3.is_finite() || volume_m3 < 0.0 || (numerical.volume > 0.0 && volume_m3 == 0.0)
        {
            return Err(MaterialLayoutSolveError::UnrepresentableVolume);
        }
        Ok(MaterialLayoutSolution {
            forces,
            dual,
            numerical,
            volume_m3,
        })
    }

    /// Refit equilibrium and size every survivor with its own resolved yield
    /// stress and modulus. Mass comes from its source density and actual snapped
    /// geometry, so swapping a material propagates through all three stages.
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
                (state.yield_stress, state.youngs)
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
