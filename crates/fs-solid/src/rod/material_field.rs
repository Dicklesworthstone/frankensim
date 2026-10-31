//! Spatially heterogeneous members using the existing Cosserat rod solver.
//!
//! A segment owns its actual concentric material section. Changing one segment
//! changes equilibrium everywhere through shared nodal compatibility, not through
//! an average modulus or separately solved pieces. Reference segments have equal
//! length; cross-section/material jumps are one-dimensional idealizations, not
//! resolved joint/interface stresses. The original clamped-root, dense Newton
//! statics and body-frame tip-moment convention remain in force.
//! An explicitly supplied stress-free centerline and material frames retain
//! intrinsic camber/twist through the SAME energy, loads and thermal refresh.

use fs_exec::CancelGate;
use fs_ga::So3;
use fs_material::section::ConcentricShearSection;
use fs_material::state_point::MaterialStatePointError;
use std::sync::Arc;

/// Sourced nonuniform thermal expansion and compatible section stresses.
pub mod thermal;
/// Prescribed nodal supports and reactions using the same material energy.
pub mod supports;
/// Reference mass, gravity and root-clamp reaction in the same equilibrium.
pub mod gravity;
/// Stress-free curved centerlines and material frames, independent of material.
pub mod reference;
/// Source-bound tension-only links sharing this rod's equilibrium and gravity.
pub mod tension;
/// Frozen reference inertia and analytical force adapter for dynamic consumers.
pub mod motion;
mod derivatives;

use super::{MAX_MATERIAL_ROD_SEGMENTS, Rod, RodSection, TipLoad, TipLoading, WorldTipLoad};
use crate::SolidError;

/// Original solid/material failure, input refusal, or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum SpatialRodError {
    /// The original constitutive/kinematic/equilibrium owner refused.
    Solid(SolidError),
    /// The original thermal section refused a source, strain or derived value.
    Material(MaterialStatePointError),
    /// Invalid spatial assignment or unrepresentable result.
    Invalid(&'static str),
    /// No candidate configuration is published after cancellation.
    Cancelled,
}
impl core::fmt::Display for SpatialRodError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => write!(f, "spatial material rod: {error}"),
            Self::Material(error) => write!(f, "spatial material rod: {error}"),
            Self::Invalid(what) => write!(f, "spatial material rod: {what}"),
            Self::Cancelled => f.write_str("spatial material rod: cancelled"),
        }
    }
}
impl std::error::Error for SpatialRodError {}
impl From<SolidError> for SpatialRodError {
    fn from(error: SolidError) -> Self { Self::Solid(error) }
}
impl From<MaterialStatePointError> for SpatialRodError {
    fn from(error: MaterialStatePointError) -> Self { Self::Material(error) }
}

/// One compatible rod configuration with a source-bound section per segment.
///
/// The uniform `Rod` is kept private: exposing its first section as though it
/// described the complete heterogeneous member would give a wrong energy.
/// Sources, solved coordinates, strains and energy are available explicitly.
/// This is statics, not dynamics, conservative remapping, or material evolution.
/// Density is retained in each source but does not introduce an implicit weight.
#[derive(Clone, Debug)]
pub struct SpatialMaterialRod {
    rod: Rod,
    sections: Vec<ConcentricShearSection>,
    stiffness: Vec<RodSection>,
    reference_length_m: f64,
    thermal: Option<thermal::ThermalField>,
    supports: Vec<supports::RodSupport>,
    free_dofs: Vec<usize>,
    reference_mass: Arc<gravity::RodReferenceMass>,
    gravity: Option<gravity::GravityLoad>,
    reference_geometry: Option<Arc<reference::RodReferenceGeometry>>,
    tension_links: Vec<tension::MaterialTensionLink>,
}

impl SpatialMaterialRod {
    /// Build a straight member; spatial order is preserved, never sorted by card.
    /// Each entry spans `length_m / sections.len()` in the reference coordinate.
    /// Sources must support static, fixed-state elasticity as in `MaterialRod`.
    ///
    /// # Errors
    /// Refuses invalid length, empty/oversized assignments, unrepresentable
    /// segment lengths, and every original material-law admission failure.
    pub fn try_straight(
        length_m: f64,
        sections: Vec<ConcentricShearSection>,
    ) -> Result<Self, SpatialRodError> {
        Self::try_straight_with_mass(length_m, sections, None)
    }

    // Source refreshes retain the reference inventory instead of evaluating
    // current density on frozen geometry. Only this module's descendants use it.
    fn try_straight_with_mass(
        length_m: f64, sections: Vec<ConcentricShearSection>,
        reference_mass: Option<Arc<gravity::RodReferenceMass>>,
    ) -> Result<Self, SpatialRodError> {
        if !length_m.is_finite() || length_m <= 0.0 || sections.is_empty()
            || sections.len() > MAX_MATERIAL_ROD_SEGMENTS
        {
            return Err(SpatialRodError::Invalid("finite positive length and 1..=64 material segments are required"));
        }
        let l0 = length_m / sections.len() as f64;
        if !l0.is_finite() || l0 <= 0.0 {
            return Err(SpatialRodError::Invalid("reference segment length is unrepresentable"));
        }
        let stiffness = sections.iter().map(RodSection::try_from_concentric_material)
            .collect::<Result<Vec<_>, _>>()?;
        let reference_mass = match reference_mass {
            Some(mass) => {
                if mass.segment_mass_kg().len() != sections.len()
                    || mass.segment_length_m().to_bits() != l0.to_bits()
                {
                    return Err(SpatialRodError::Invalid("reference mass and rod discretization must match"));
                }
                mass
            }
            None => Arc::new(gravity::RodReferenceMass::try_new(l0, &sections)?),
        };
        let rod = Rod::straight(length_m, sections.len(), stiffness[0]);
        let free_dofs = (0..rod.ndof()).collect();
        Ok(Self { rod, sections, stiffness, reference_length_m: length_m, thermal: None,
            supports: Vec::new(), free_dofs, reference_mass, gravity: None,
            reference_geometry: None, tension_links: Vec::new() })
    }

    /// Original material assignments and all layer receipts, in spatial order.
    #[must_use]
    pub fn sections(&self) -> &[ConcentricShearSection] { &self.sections }
    /// Accepted nodal positions [m], shared by all incident material segments.
    #[must_use]
    pub fn positions_m(&self) -> &[[f64; 3]] { &self.rod.positions }
    /// Accepted body-to-world nodal frames from the existing SO(3) owner.
    #[must_use]
    pub fn frames(&self) -> &[So3] { &self.rod.frames }
    /// Authored material reference length [m], not the current arc length.
    #[must_use]
    pub const fn reference_length_m(&self) -> f64 { self.reference_length_m }

    /// Elastic Cosserat strain/curvature relative to the stress-free geometry.
    /// Thermal eigenstrain is separate and is applied by the section owner.
    pub fn segment_strains(&self, segment: usize) -> Result<([f64; 3], [f64; 3]), SpatialRodError> {
        if segment >= self.sections.len() {
            return Err(SpatialRodError::Invalid("segment index is outside the material field"));
        }
        self.strains_for(&self.rod, segment)
    }

    /// Local force [N] and moment [N m] resultants at the accepted configuration.
    /// These are section resultants, not three-dimensional interfacial stresses.
    pub fn segment_resultants(&self, segment: usize) -> Result<([f64; 3], [f64; 3]), SpatialRodError> {
        let (g, k) = self.segment_strains(segment)?;
        self.resultants_from_strains(segment, g, k)
    }

    // The same constitutive response supplies observations and analytical
    // assembly. No separately copied thermal-force or intrinsic-strain law.
    fn resultants_from_strains(&self, segment: usize, g: [f64; 3], k: [f64; 3])
        -> Result<([f64; 3], [f64; 3]), SpatialRodError> {
        let s = self.stiffness[segment];
        let axial = match &self.thermal {
            Some(field) => field.summary(segment, g[0])?.0,
            None => s.ea * g[0],
        };
        let force = [axial, s.ga * g[1], s.ga * g[2]];
        let moment = [s.gj * k[0], s.ei * k[1], s.ei * k[2]];
        if !force.iter().chain(&moment).all(|v| v.is_finite()) {
            return Err(SpatialRodError::Invalid("section resultants are unrepresentable"));
        }
        Ok((force, moment))
    }

    /// Stored strain energy [J] from every segment and attached tensile link.
    pub fn energy_j(&self) -> Result<f64, SpatialRodError> {
        self.energy_for(&self.rod)
    }

    // Retained scale for independent numerical oracles and callback fallbacks.
    // Native analytical assembly does not perturb the material strain domain.
    fn probe_translation_scale_m(&self) -> f64 {
        if self.thermal.is_some() || self.reference_geometry.is_some() || !self.tension_links.is_empty() {
            self.rod.l0
        } else { 1.0 }
    }

    fn strains_for(&self, rod: &Rod, segment: usize)
        -> Result<([f64; 3], [f64; 3]), SpatialRodError> {
        self.reference_relative_strains(segment, rod.strains(segment)?)
    }

    fn reference_relative_strains(&self, segment: usize, (mut g, mut k): ([f64; 3], [f64; 3]))
        -> Result<([f64; 3], [f64; 3]), SpatialRodError> {
        if let Some(reference) = &self.reference_geometry {
            let (g0, k0) = reference.strains[segment];
            for c in 0..3 {
                g[c] -= g0[c];
                k[c] -= k0[c];
            }
        }
        if !g.iter().chain(&k).all(|v| v.is_finite()) {
            return Err(SpatialRodError::Invalid("relative material rod strain is unrepresentable"));
        }
        Ok((g, k))
    }

    // One constitutive owner for total storage and element-local derivatives.
    fn segment_energy_for(&self, rod: &Rod, segment: usize) -> Result<f64, SpatialRodError> {
        let section = &self.stiffness[segment];
        let (g, k) = self.strains_for(rod, segment)?;
        let energy = match &self.thermal {
            Some(field) => {
                let axial_energy = field.summary(segment, g[0])?.1;
                rod.l0 * axial_energy + 0.5 * rod.l0
                    * section.twice_energy_per_length([0.0, g[1], g[2]], k)
            }
            None => 0.5 * rod.l0 * section.twice_energy_per_length(g, k),
        };
        if !energy.is_finite() {
            return Err(SpatialRodError::Invalid("material segment strain energy is unrepresentable"));
        }
        Ok(energy)
    }

    fn energy_for(&self, rod: &Rod) -> Result<f64, SpatialRodError> {
        let mut energy = 0.0;
        for segment in 0..self.stiffness.len() {
            energy += self.segment_energy_for(rod, segment)?;
        }
        for link in &self.tension_links {
            energy += link.state_for(rod)?.stored_energy_j;
        }
        if !energy.is_finite() {
            return Err(SpatialRodError::Invalid("material strain energy is unrepresentable"));
        }
        Ok(energy)
    }

    /// Solve ONE compatible equilibrium using the original Newton/Lie-group path.
    /// Analytical segment, nonlocal tensile-link and rotating-gravity forces
    /// and tangents share the canonical SO(3) differential and section owners.
    /// Actual states still pass full energy/source-domain admission. Slack links
    /// remain potential interactions; their taut/slack kink uses the explicitly
    /// selected slack-sided derivative. The original dense linear solve remains.
    /// Cancellation is checked at load/Newton steps, each element/link and
    /// around the bounded dense solve. There is no sparse-solver or dynamics claim.
    /// Only a complete finite candidate replaces the previous geometry.
    ///
    /// `tolerance` retains the legacy norm mixing force and moment residuals;
    /// it is numerical stopping evidence, not a dimensionless error certificate.
    /// Explicit gravity, when present, is a permanent load at every stage;
    /// `steps` ramps only the tip load. Stored and gravitational energy remain
    /// separate observations, not a coupled thermodynamic balance.
    ///
    /// # Errors
    /// Refuses invalid loads/budgets, all source/solver failures and cancellation.
    /// Every failure, including a late load-step failure, preserves the old pose.
    pub fn solve_static(
        &mut self,
        load: &TipLoad,
        steps: usize,
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, SpatialRodError> {
        self.solve_static_loading(load.into(), steps, tolerance, gate)
    }

    /// Equilibrate the original spatial/thermal material rod under a world-held
    /// tip force and couple. Existing supports, intrinsic shape, gravity and
    /// tensile members compose in the same residual and consistent tangent.
    /// Only a complete accepted pose replaces the original state.
    ///
    /// # Errors
    /// Refuses invalid loads/budgets, source/solver failures and cancellation,
    /// with the same atomic state boundary as [`Self::solve_static`].
    pub fn solve_static_world(
        &mut self,
        load: &WorldTipLoad,
        steps: usize,
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, SpatialRodError> {
        self.solve_static_loading(load.into(), steps, tolerance, gate)
    }

    fn solve_static_loading(
        &mut self,
        load: TipLoading,
        steps: usize,
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, SpatialRodError> {
        checkpoint(gate)?;
        if steps == 0 || steps > 256 || !tolerance.is_finite() || tolerance <= 0.0
            || !load.is_finite()
        {
            return Err(SpatialRodError::Invalid("finite loads, positive tolerance and 1..=256 load steps are required"));
        }
        let mut candidate = self.rod.clone();
        let energy = derivatives::SpatialEnergy::new(self);
        let residuals = candidate.solve_static_with_energy(
            load, steps, tolerance, &energy, &|| checkpoint(gate),
            Some(&self.free_dofs), self.probe_translation_scale_m(),
        )?;
        self.equilibrium_energy_for(&candidate)?;
        if !candidate.positions.iter().flatten().all(|v| v.is_finite()) {
            return Err(SpatialRodError::Invalid("solved material rod position is unrepresentable"));
        }
        checkpoint(gate)?;
        self.rod = candidate;
        Ok(residuals)
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), SpatialRodError> {
    if gate.is_requested() { Err(SpatialRodError::Cancelled) } else { Ok(()) }
}
