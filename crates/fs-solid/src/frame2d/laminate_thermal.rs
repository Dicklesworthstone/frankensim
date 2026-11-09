//! Thermal-gradient and bimetal bending in the original force-based frame.
//!
//! The shared material owner supplies a narrow free-sided Euler–Bernoulli
//! section: each layer relaxes transverse normal/shear stress, retaining its
//! original elastic tensor and integrated thermal strain. This is not a
//! cylindrically constrained plate strip. The reference centerline follows the
//! original density centroid, with `epsilon_x = epsilon_0 - y * kappa`.
//! An elastic centroid elsewhere produces genuine axial/bending coupling.
//!
//! Reference thickness, width, mass and rotary inertia are fixed. Heating may
//! change E(T), free curvature, prestress and vibration, but does not replace a
//! closed specimen by its endpoint density. Transverse warping is quasistatic;
//! no transverse-mode inertia, interlaminar stress, torsion, damage or heat
//! evolution is supplied. Prescribed-temperature work belongs to the caller.

use super::{ElementCommit, Frame2d, SectionModel};
use crate::SolidError;
use crate::fiber::SectionState;
use fs_exec::Cx;
use fs_material::laminate::strength::PreparedLaminateStrength;
use fs_material::laminate::thermal::LayerThermalStress;
use fs_material::laminate::thermal::beam::PreparedLaminateBeam;
use fs_material::laminate::thermal::beam::thermoelastic::ResolvedLaminateBeamThermoelastic;
use fs_material::laminate::thermal::source::{LaminateThermalSource, ResolvedLaminateThermal};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use std::sync::Arc;

/// Source-pinned maximum-stress indications at accepted native stations.
pub mod strength;
pub use strength::FrameLaminateStrengthAssessment;

/// Original-source reversible storage and instantaneous thermal feedback.
pub mod thermoelastic;

fn invalid(what: &str) -> SolidError {
    SolidError::MaterialRefused { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}
fn source_error(error: impl core::fmt::Display, cx: &Cx<'_>) -> SolidError {
    if cx.checkpoint().is_err() {
        SolidError::Cancelled
    } else {
        invalid(&error.to_string())
    }
}

/// Original source/reference inventory shared by every station and refresh.
#[derive(Clone, Debug)]
pub(super) struct FrameLaminateBinding {
    source: Arc<LaminateThermalSource>,
    width_m: f64,
    axis_offset_m: f64,
    linear_strain_limit: f64,
    // Original per-station selections. Refresh evaluates these exact claims,
    // even when another claim becomes usable at the new temperature.
    strength: Option<Vec<Arc<PreparedLaminateStrength>>>,
    thermoelastic: Option<Arc<Vec<AffineLaminaThermoelastic>>>,
}

/// Immutable current station with original material cards and reference paths.
#[derive(Clone, Debug)]
pub struct FrameLaminateSection {
    binding: Arc<FrameLaminateBinding>,
    resolved: Arc<ResolvedLaminateThermal>,
    prepared: PreparedLaminateBeam,
    strength: Option<Arc<PreparedLaminateStrength>>,
    thermoelastic: Option<Arc<ResolvedLaminateBeamThermoelastic>>,
}

impl FrameLaminateSection {
    /// Current per-layer elastic and integrated-expansion receipts.
    #[must_use]
    pub fn source(&self) -> &ResolvedLaminateThermal {
        &self.resolved
    }

    /// Shared narrow-beam law, including its original section and axis shift.
    #[must_use]
    pub const fn prepared(&self) -> &PreparedLaminateBeam {
        &self.prepared
    }

    /// Current strength sources bound to this exact mechanical snapshot.
    /// Present only after [`Frame2d::with_laminate_strength`] succeeds.
    #[must_use]
    pub fn strength(&self) -> Option<&PreparedLaminateStrength> {
        self.strength.as_deref()
    }

    /// Reversible potential attached to this exact current mechanical snapshot.
    #[must_use]
    pub fn thermoelastic(&self) -> Option<&ResolvedLaminateBeamThermoelastic> {
        self.thermoelastic.as_deref()
    }

    /// Bound on total, free and elastic in-plane strain tensor norms.
    #[must_use]
    pub fn linear_strain_limit(&self) -> f64 {
        self.binding.linear_strain_limit
    }

    fn admit_strain(&self, eps0: f64, kappa: f64) -> Result<(), SolidError> {
        for layer in 0..self.prepared.source().section().layers().len() {
            let total = self
                .prepared
                .layer_strain_faces(layer, eps0, kappa)
                .map_err(|error| invalid(&error.to_string()))?;
            let free = self
                .prepared
                .source()
                .layer_free_strain_in_section(layer)
                .map_err(|error| invalid(&error.to_string()))?;
            for face in 0..2 {
                let elastic = std::array::from_fn(|axis| total[face][axis] - free[face][axis]);
                for strain in [total[face], free[face], elastic] {
                    // Frobenius norm of the in-plane symmetric strain tensor;
                    // gamma_xy is the engineering shear, twice epsilon_xy.
                    let norm = strain[0]
                        .hypot(strain[1])
                        .hypot(strain[2] / std::f64::consts::SQRT_2);
                    if !norm.is_finite() || norm > self.binding.linear_strain_limit {
                        return Err(invalid(
                            "laminate frame layer strain requires finite-strain or inelastic escalation",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    // Thermal integration and tensor condensation occur once per temperature
    // snapshot. Repeated native residual/tangent evaluations allocate nothing.
    pub(super) fn respond(&self, eps0: f64, kappa: f64) -> Result<SectionState, SolidError> {
        self.admit_strain(eps0, kappa)?;
        let tangent = self.prepared.tangent();
        // Integrate the original layer elastic strain before multiplication;
        // subtracting two large section resultants can erase a small stress.
        let response = self
            .prepared
            .evaluate(eps0, kappa, || false)
            .map_err(|error| invalid(&error.to_string()))?;
        Ok(SectionState {
            n: response.axial_force_n,
            m: response.bending_moment_n_m,
            tangent,
        })
    }
}

/// Physical response at an accepted native station about the mass centerline.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameLaminateSectionResponse {
    /// Total original-centerline axial strain.
    pub axial_strain: f64,
    /// Native curvature [1/m]; positive curvature shortens positive-y fibers.
    pub curvature_per_m: f64,
    /// Physical axial force [N], including thermal restraint and coupling.
    pub axial_force_n: f64,
    /// Work-conjugate signed bending moment [N m].
    pub bending_moment_n_m: f64,
    /// Elastic storage per original member length [J/m], including residual
    /// self-stress even when both free-section resultants vanish.
    pub strain_energy_j_m: f64,
    /// Original layer bottom/top stresses in section axes. No interlaminar
    /// traction is inferred from the free-sided beam approximation.
    pub layer_stress_pa: Vec<LayerThermalStress>,
}

fn prepare_stations(
    binding: &Arc<FrameLaminateBinding>,
    temperatures_k: &[Vec<f64>],
    cx: &Cx<'_>,
) -> Result<Vec<SectionModel>, SolidError> {
    let mut sections = Vec::with_capacity(temperatures_k.len());
    for (station, temperatures) in temperatures_k.iter().enumerate() {
        checkpoint(cx)?;
        let (resolved, strength) = if let Some(original) = &binding.strength {
            let current = original[station]
                .at_temperatures(temperatures, || cx.checkpoint().is_err())
                .map_err(|error| source_error(error, cx))?;
            // The strength owner refreshes mechanics and all original strength
            // claims together. Clone its exact prepared Arc, not a separately
            // queried equivalent section: temperature-dependent strength can
            // change even when every mechanical coefficient remains equal.
            (current.thermal().clone(), Some(Arc::new(current)))
        } else {
            (
                binding
                    .source
                    .at_temperatures(temperatures, || cx.checkpoint().is_err())
                    .map_err(|error| source_error(error, cx))?,
                None,
            )
        };
        let prepared = PreparedLaminateBeam::try_new(
            Arc::clone(resolved.prepared()),
            binding.width_m,
            binding.axis_offset_m,
            || cx.checkpoint().is_err(),
        )
        .map_err(|error| source_error(error, cx))?;
        let thermoelastic = binding
            .thermoelastic
            .as_ref()
            .map(|laws| {
                ResolvedLaminateBeamThermoelastic::try_new(laws, &resolved, &prepared, || {
                    cx.checkpoint().is_err()
                })
                .map(Arc::new)
                .map_err(|error| source_error(error, cx))
            })
            .transpose()?;
        let section = FrameLaminateSection {
            binding: Arc::clone(binding),
            resolved: Arc::new(resolved),
            prepared,
            strength,
            thermoelastic,
        };
        // Admit complete original thermal strain, including opposing strains
        // that would disappear in an averaged section resultant.
        section.admit_strain(0.0, 0.0)?;
        sections.push(SectionModel::LaminateThermal(section));
    }
    checkpoint(cx)?;
    Ok(sections)
}

impl Frame2d {
    /// Add a source-bound narrow laminate with one temperature per original
    /// layer at each of 3..=7 Gauss–Lobatto stations (including both ends).
    /// Temperatures may differ both through thickness and along the member.
    /// Every path starts from the source's ORIGINAL stress-free temperature;
    /// material cards, layer order, orientation and geometry remain attached.
    ///
    /// Frame nodes lie on the original density-centroid axis. The shared
    /// material owner shifts the constitutive law and layer recovery to that
    /// axis; a nonzero elastic first moment remains in the coupled tangent.
    /// Original mass and centered rotary inertia enter the native mass owner,
    /// without discarding a first mass moment about an arbitrary source origin.
    ///
    /// Each layer is free-sided (`sigma_yy = tau_xy = 0`). Mechanical or thermal
    /// shear/twist coupling outside this planar model refuses. This is a narrow
    /// Euler–Bernoulli member, not a constrained plate, Timoshenko shear member,
    /// finite-strain solid or thermal-history law. `linear_strain_limit` must
    /// lie in (0, 0.05]. All failures/cancellation leave the frame unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn add_laminate_thermal_element(
        &mut self,
        i: usize,
        j: usize,
        source: Arc<LaminateThermalSource>,
        width_m: f64,
        temperatures_k: &[Vec<f64>],
        linear_strain_limit: f64,
        cx: &Cx<'_>,
    ) -> Result<usize, SolidError> {
        checkpoint(cx)?;
        if super::lobatto(temperatures_k.len()).is_none()
            || !linear_strain_limit.is_finite()
            || linear_strain_limit <= 0.0
            || linear_strain_limit > 0.05
        {
            return Err(invalid(
                "laminate frame requires 3..=7 stations and a strain bound in (0, 0.05]",
            ));
        }
        let reference = source.reference_section();
        let axis_offset_m = reference.i1_kg_m() / reference.i0_kg_m2();
        let binding = Arc::new(FrameLaminateBinding {
            source,
            width_m,
            axis_offset_m,
            linear_strain_limit,
            strength: None,
            thermoelastic: None,
        });
        let sections = prepare_stations(&binding, temperatures_k, cx)?;
        let SectionModel::LaminateThermal(first) = &sections[0] else {
            unreachable!("prepared stations are laminate sections")
        };
        let [mass, first_moment, spin] = first.prepared.reference_mass_moments();
        let h = first.prepared.source().section().total_thickness_m();
        // A scalar line mass plus rotary inertia cannot represent an omitted
        // first moment. Only floating-point centering residual is admitted.
        let moment_roundoff = 256.0 * f64::EPSILON * mass * h;
        if !first_moment.is_finite()
            || !moment_roundoff.is_finite()
            || first_moment.abs() > moment_roundoff
        {
            return Err(invalid(
                "laminate frame reference mass centroid is unrepresentable",
            ));
        }
        let first = sections[0].clone();
        let mut candidate = self.clone();
        let element =
            candidate.add_element_with_rotary_inertia(i, j, sections.len(), mass, spin, &|| {
                first.clone()
            })?;
        candidate.elements[element].laminate_thermal = Some(binding);
        candidate.elements[element].sections = sections;
        let trial = candidate.element_trial(&candidate.elements[element], &candidate.u)?;
        candidate.elements[element].commit = ElementCommit {
            q: trial.q,
            d: trial.d,
            alpha: trial.alpha,
        };
        checkpoint(cx)?;
        *self = candidate;
        Ok(element)
    }

    /// Refresh prescribed layer temperatures through the same original source.
    /// Source selection, references, geometry, width, strain bound and mass
    /// centroid cannot be replaced or silently rebased. Nodal configuration,
    /// velocities and original inertia are preserved; native compatibility is
    /// solved at that pose before any candidate is returned. The caller then
    /// owns supported equilibrium or frozen-temperature time advancement.
    /// Attached strength claims refresh at the same temperatures and remain
    /// pinned; any unavailable strength refuses the complete thermal update.
    /// Attached reversible potentials retain their original smooth intervals
    /// and share that same snapshot; an unsupported interval also refuses.
    pub fn with_laminate_temperatures(
        &self,
        element: usize,
        temperatures_k: &[Vec<f64>],
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("laminate thermal frame member index is out of range"))?;
        let binding = member.laminate_thermal.as_ref().ok_or_else(|| {
            invalid(
                "laminate temperature refresh requires an original source-bound laminate member",
            )
        })?;
        if temperatures_k.len() != member.sections.len()
            || member.transverse_shear.is_some()
            || member.viscous_section.is_some()
            || member.sections.iter().any(|section| {
                !matches!(section, SectionModel::LaminateThermal(old)
                    if Arc::ptr_eq(&old.binding, binding))
            })
        {
            return Err(invalid(
                "laminate refresh must retain every original elastic station without a separate shear or viscous operator",
            ));
        }
        let sections = prepare_stations(binding, temperatures_k, cx)?;
        for section in &sections {
            let SectionModel::LaminateThermal(section) = section else {
                unreachable!()
            };
            let [mass, _, spin] = section.prepared.reference_mass_moments();
            if mass.to_bits() != member.mass_per_length.to_bits()
                || spin.to_bits() != member.rotary_inertia_per_length.to_bits()
            {
                return Err(invalid(
                    "laminate thermal refresh changed the original mass inventory",
                ));
            }
        }
        let mut candidate = self.clone();
        candidate.elements[element].sections = sections;
        let trial = candidate.element_trial(&candidate.elements[element], &candidate.u)?;
        candidate.elements[element].commit = ElementCommit {
            q: trial.q,
            d: trial.d,
            alpha: trial.alpha,
        };
        if binding.thermoelastic.is_some() {
            // Admit the potential at the newly compatible native strains as
            // part of the same transaction, before returning the hot frame.
            candidate.laminate_thermoelastic_energy(element, cx)?;
        }
        checkpoint(cx)?;
        Ok(candidate)
    }

    /// Current source-bound laminate station, if the member carries one.
    #[must_use]
    pub fn laminate_thermal_section(
        &self,
        element: usize,
        station: usize,
    ) -> Option<&FrameLaminateSection> {
        match self.elements.get(element)?.sections.get(station)? {
            SectionModel::LaminateThermal(section) => Some(section),
            _ => None,
        }
    }

    /// Recover accepted native force, curvature, layer stresses and full stored
    /// energy. Face stresses retain bimetal mismatch hidden by net resultants.
    pub fn laminate_thermal_response(
        &self,
        element: usize,
        station: usize,
    ) -> Result<FrameLaminateSectionResponse, SolidError> {
        let section = self
            .laminate_thermal_section(element, station)
            .ok_or_else(|| {
                invalid("requested frame station has no source-bound thermal laminate")
            })?;
        let (eps0, kappa) = self.elements[element].commit.d[station];
        section.admit_strain(eps0, kappa)?;
        let response = section
            .prepared
            .evaluate(eps0, kappa, || false)
            .map_err(|error| invalid(&error.to_string()))?;
        let layer_stress_pa = (0..section.prepared.source().section().layers().len())
            .map(|layer| {
                section
                    .prepared
                    .layer_stress(layer, eps0, kappa)
                    .map_err(|error| invalid(&error.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(FrameLaminateSectionResponse {
            axial_strain: eps0,
            curvature_per_m: kappa,
            axial_force_n: response.axial_force_n,
            bending_moment_n_m: response.bending_moment_n_m,
            strain_energy_j_m: response.strain_energy_j_m,
            layer_stress_pa,
        })
    }
}
