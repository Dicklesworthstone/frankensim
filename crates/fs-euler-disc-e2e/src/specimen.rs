//! Resolved physical specimens for the Euler-disc campaign.
//!
//! A [`DiscProfileSpec`] is deliberately a small *description* of an
//! axisymmetric solid.  Resolving it produces one `AxisymmetricChart` and then
//! obtains mass, centroid, and principal inertia from that exact same line/arc
//! profile.  A dynamics caller must use the resolved chart for support queries
//! as well; hand-entering cylinder or cone inertia alongside a different
//! contact shape is specifically what this boundary prevents.
//!
//! The geometry and mass integrations are analytic over the represented
//! line/arc profile, but their binary64 evaluations carry `Estimate`/roundoff
//! telemetry only.  This module is therefore an input-consistency layer, not a
//! contact-patch, material-calibration, experimental-validation, or
//! configuration-ranking claim.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_conduction::lumped::{LumpedEnthalpyBody, LumpedEnthalpyMarch};
use fs_exec::Cx;
use fs_material::{
    phase::{EquilibriumPhaseState, SolidLiquidPhase},
    state_point::{
        IsotropicElasticStatePoint, IsotropicSolidStatePoint, OrthotropicElasticStatePoint,
    },
};
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricError, AxisymmetricIdentity, AxisymmetricMassError,
    AxisymmetricMassProperties, AxisymmetricSurfaceAreaError, MeridianPoint, MeridianSegment,
    SquatDiscEdgeTreatment,
};
use fs_solid::{TetElasticError, TetElasticMaterial};

/// Canonical identity domain for the exact retained axisymmetric chart input.
pub const EULER_SPECIMEN_CHART_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.specimen-chart.v1";
/// Canonical identity domain for resolved geometry plus homogeneous density.
pub const EULER_SPECIMEN_PROFILE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.specimen-profile.v1";
/// Canonical identity domain for the resolved production mass properties.
pub const EULER_SPECIMEN_MASS_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.specimen-mass-properties.v1";
/// Canonical identity domain for geometry-derived lumped thermal measures.
pub const EULER_SPECIMEN_THERMAL_GEOMETRY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.specimen-thermal-geometry.v1";
/// Canonical identity domain for a thermal march bound to one exact specimen.
pub const EULER_SPECIMEN_THERMAL_MARCH_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.specimen-thermal-march.v1";
/// Canonical identity domain for the bounded free-isotropic solid expansion law.
pub const EULER_ISOTROPIC_FREE_EXPANSION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.isotropic-free-expansion-law.v1";
/// Canonical identity domain for one thermally evolved solid profile.
pub const EULER_EVOLVED_SOLID_PROFILE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.evolved-solid-profile.v1";
/// Canonical identity domain for geometry plus its resolved material state.
pub const EULER_MATERIAL_SPECIMEN_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.material-specimen.v1";

/// The bounded, user-facing profile families admitted by the Euler campaign.
///
/// Every variant revolves a single simple meridian about its local `z` axis.
/// The local origin is a construction coordinate, not necessarily the center
/// of mass: [`ResolvedDiscProfile::mass_properties`] is authoritative for the
/// latter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiscProfileSpec {
    /// A solid squat cylinder with either sharp or true circular-filleted rims.
    SolidCylinder {
        /// Maximum cylindrical radius (m).
        outer_radius_m: f64,
        /// Distance between the two cap planes (m).
        thickness_m: f64,
        /// Exact outer-rim treatment.
        edge_treatment: SquatDiscEdgeTreatment,
    },
    /// A homogeneous annular cylinder.  Its inner bore is physical geometry,
    /// rather than a mass-only adjustment to a solid-cylinder contact shape.
    AnnularCylinder {
        /// Outer cylindrical radius (m).
        outer_radius_m: f64,
        /// Radius of the through bore (m).
        inner_radius_m: f64,
        /// Distance between cap planes (m).
        thickness_m: f64,
    },
    /// An annular cylinder with true circular fillets at both outer rims.
    ///
    /// The through bore and its cap intersections remain deliberately sharp.
    /// The strictly positive outer fillet is resolved as two circular meridian
    /// arcs, so no sharp ring-contact surrogate or independent inertia factor
    /// is introduced.
    OuterFilletedAnnularCylinder {
        /// Outer cylindrical radius (m).
        outer_radius_m: f64,
        /// Radius of the through bore (m).
        inner_radius_m: f64,
        /// Distance between cap planes (m).
        thickness_m: f64,
        /// Radius of each true outer-rim meridian fillet (m). Must be positive.
        outer_fillet_radius_m: f64,
    },
    /// A symmetric double-conical or double-frustum profile.
    ///
    /// `face_radius_m == 0` is a true bicone whose two points lie on the
    /// revolution axis.  Positive values create equal planar end faces and
    /// true conical flanks.  This parameterization keeps the material
    /// symmetric about `z = 0` while making contact geometry differ from a
    /// cylinder.
    SymmetricTapered {
        /// Maximum radius at the equatorial plane (m).
        outer_radius_m: f64,
        /// Radius of each planar end face (m).
        face_radius_m: f64,
        /// Axial tip-to-tip/end-face separation (m).
        thickness_m: f64,
    },
    /// A solid cylinder with equal straight conical chamfers at both rims.
    ///
    /// Both chamfer distances must be positive.  Use `SolidCylinder::Sharp`
    /// to request no chamfer; accepting a half-zero chamfer would conceal a
    /// caller-unit or topology mistake.
    ChamferedCylinder {
        /// Maximum radius on the retained cylindrical band (m).
        outer_radius_m: f64,
        /// Distance between cap planes (m).
        thickness_m: f64,
        /// Radial inset from the cylindrical band to either cap edge (m).
        chamfer_radial_m: f64,
        /// Axial run of either conical chamfer (m).
        chamfer_axial_m: f64,
    },
}

/// Geometry extents declared by a resolved profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscProfileDimensions {
    /// Largest radial extent of the profile (m).
    pub outer_radius_m: f64,
    /// Difference between maximum and minimum meridian axial coordinates (m).
    pub thickness_m: f64,
}

/// A profile resolved into one validated geometry and its matching properties.
#[derive(Clone, Debug)]
pub struct ResolvedDiscProfile {
    /// Original bounded parameterization for provenance and JSON reporting.
    pub spec: DiscProfileSpec,
    /// Validated line/arc solid of revolution used for every support query.
    pub chart: AxisymmetricChart,
    /// Homogeneous volumetric density used for the mass integration (kg/m³).
    pub density_kg_per_m3: f64,
    /// Deterministic identity of the exact retained meridian input.
    pub identity: AxisymmetricIdentity,
    /// Stable outer-radius and axial-thickness dimensions (m).
    pub dimensions: DiscProfileDimensions,
    /// Analytic line/arc mass, center of mass, and centroidal inertia.
    pub mass_properties: AxisymmetricMassProperties,
}

/// Whole-boundary measures for an isothermal-body thermal rung.
///
/// Both values derive from the same exact line/arc chart and volume already
/// used by mechanics. `characteristic_length_m = volume_m3 / surface_area_m2`
/// is the conventional body length used by a Biot-number admission gate.
/// This whole-boundary value is only appropriate when the declared heat-
/// transfer condition acts over the complete exposed surface; spatially
/// varying or partially insulated boundaries require a partitioned thermal
/// mesh rather than an adjusted caller-entered area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedDiscThermalGeometry {
    /// Exact-formula binary64 area of the complete revolved boundary [m2].
    pub surface_area_m2: f64,
    /// Enclosed volume already bound to the resolved mass properties [m3].
    pub volume_m3: f64,
    /// Whole-body characteristic length `V/A` [m].
    pub characteristic_length_m: f64,
    /// Identity binding the exact chart and all three thermal measures.
    pub identity: ContentHash,
}

/// One axisymmetric specimen whose density and remaining elastic properties
/// come from the same evidence-bearing material state point.
///
/// The contained profile remains the geometry/mass source consumed by legacy
/// mechanics adapters. This wrapper adds the full material-state identity so
/// equal-density but physically different materials cannot alias in an
/// integrated simulation or render.
#[derive(Clone, Debug)]
pub struct ResolvedMaterialDiscProfile {
    /// Geometry and mass properties evaluated from the resolved density.
    pub profile: ResolvedDiscProfile,
    /// Complete isotropic material property/receipt bundle.
    pub material: IsotropicSolidStatePoint,
    /// Identity binding geometry, mass, and complete material state.
    pub identity: ContentHash,
}

/// Axisymmetric geometry and mass bound to a complete linear-elastic tensor.
///
/// This is the material-symmetry-independent specimen consumed by structural
/// modes and vibroacoustics. Contact/plastic solvers may require additional
/// state such as yield and interface properties, but they must reuse this
/// profile and its card/state identities rather than constructing a second
/// material by name.
#[derive(Clone, Debug)]
pub struct ResolvedElasticDiscProfile {
    /// Geometry and mass properties evaluated from the same resolved density.
    pub profile: ResolvedDiscProfile,
    /// Global-frame isotropic or oriented-anisotropic elastic tensor.
    pub elastic_material: TetElasticMaterial,
    /// Immutable material card from which the bulk properties resolved.
    pub material_card_identity: ContentHash,
    /// Complete constitutive-state identity; for anisotropy this includes the
    /// independently identified material-axis orientation.
    pub material_state_identity: ContentHash,
    /// Identity binding geometry, mass, card, and complete constitutive state.
    pub identity: ContentHash,
}

impl ResolvedElasticDiscProfile {
    fn bind(
        profile: ResolvedDiscProfile,
        elastic_material: TetElasticMaterial,
        material_card_identity: ContentHash,
    ) -> Self {
        let material_state_identity = ContentHash(elastic_material.material_state_identity);
        let profile_identities = profile.content_identities();
        let mut identity =
            DomainHasher::new("org.frankensim.fs-euler-disc-e2e.elastic-material-specimen.v1");
        identity.update(profile_identities.profile.as_bytes());
        identity.update(profile_identities.mass_properties.as_bytes());
        identity.update(material_card_identity.as_bytes());
        identity.update(material_state_identity.as_bytes());
        Self {
            profile,
            elastic_material,
            material_card_identity,
            material_state_identity,
            identity: identity.finalize(),
        }
    }
}

impl From<&ResolvedMaterialDiscProfile> for ResolvedElasticDiscProfile {
    fn from(specimen: &ResolvedMaterialDiscProfile) -> Self {
        Self::bind(
            specimen.profile.clone(),
            TetElasticMaterial::from_resolved_state(&specimen.material),
            specimen.material.resolved().card_identity(),
        )
    }
}

/// Refusal while resolving geometry and a spatial elastic constitutive state.
#[derive(Clone, Debug, PartialEq)]
pub enum ElasticDiscProfileError {
    /// Axisymmetric geometry or mass integration refused.
    Profile(DiscProfileError),
    /// Elastic tensor or material orientation refused.
    Elastic(TetElasticError),
}

impl fmt::Display for ElasticDiscProfileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ElasticDiscProfileError {}

impl From<DiscProfileError> for ElasticDiscProfileError {
    fn from(source: DiscProfileError) -> Self {
        Self::Profile(source)
    }
}

impl From<TetElasticError> for ElasticDiscProfileError {
    fn from(source: TetElasticError) -> Self {
        Self::Elastic(source)
    }
}

/// One reference profile resolved at a thermodynamic phase state.
///
/// This carrier remains valid when the state is mushy or liquid, but its
/// `profile` is only the reference axisymmetric material domain evaluated at
/// the supplied equilibrium bulk density. It does not claim the current free
/// surface or deformation. A fixed-solid consumer must call
/// [`Self::try_bind_fixed_solid`] before using rigid or small-strain mechanics.
#[derive(Clone, Debug)]
pub struct ResolvedPhaseDiscProfile {
    /// Reference geometry and mass properties at the phase-state density.
    pub profile: ResolvedDiscProfile,
    /// Equilibrium temperature, density, and solid/liquid fractions.
    pub phase_state: EquilibriumPhaseState,
    /// Identity binding reference geometry and the complete phase state.
    pub identity: ContentHash,
}

/// Geometry response required by a mass-conserving phase state.
///
/// The enum is an admission result, not a deformation model. It prevents an
/// enthalpy solver from silently reusing stale rigid geometry after density or
/// phase changes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscPhaseGeometryRegime {
    /// The exact reference phase state and geometry remain applicable.
    ReferenceGeometry,
    /// The body remains solid, but thermomechanical properties and geometry
    /// must be recomputed at the new state before mechanics continues.
    SolidThermomechanicalUpdateRequired,
    /// A nonzero liquid fraction requires an evolving free-surface/topology
    /// solver rather than fixed-solid or shape-similarity mechanics.
    EvolvingFreeSurfaceRequired,
}

/// One phase state constrained to retain the reference specimen's mass.
///
/// `required_volume_m3 = invariant_mass_kg / bulk_density_kg_m3` is a scalar
/// conservation constraint only. It deliberately does not invent an isotropic
/// scale, thermal strain field, or liquid free surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassConservingDiscPhaseState {
    /// State produced by the material-card-bound phase curve.
    pub phase_state: EquilibriumPhaseState,
    /// Material mass retained from the reference specimen [kg].
    pub invariant_mass_kg: f64,
    /// Volume required by mass conservation at the state's bulk density [m3].
    pub required_volume_m3: f64,
    /// Required volume divided by the reference volume.
    pub volume_ratio: f64,
    /// Solver rung required before mechanics, sound, or rendering continues.
    pub geometry_regime: DiscPhaseGeometryRegime,
    /// Identity binding reference specimen, phase state, and conservation data.
    pub identity: ContentHash,
}

/// Bounded law for a homogeneous, unconstrained, isotropic solid in thermal equilibrium.
///
/// The phase curve supplies equilibrium density, so mass conservation fixes the
/// volume ratio. This law supplies the additional constitutive assertion that
/// the free thermal strain is isotropic, making the unique linear scale the
/// cube root of that ratio. The authority identity must refer to evidence or a
/// caller declaration for that assertion; a chemistry name never selects it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformIsotropicFreeExpansionLaw {
    maximum_absolute_linear_strain: f64,
    authority_identity: ContentHash,
    identity: ContentHash,
}

impl UniformIsotropicFreeExpansionLaw {
    /// Construct a bounded free-expansion law.
    ///
    /// `maximum_absolute_linear_strain` is a dimensionless validity bound on
    /// `abs(linear_scale - 1)`. Zero identities and nonpositive/nonfinite
    /// bounds refuse rather than creating anonymous constitutive authority.
    pub fn try_new(
        maximum_absolute_linear_strain: f64,
        authority_identity: ContentHash,
    ) -> Result<Self, SolidGeometryEvolutionError> {
        if !(maximum_absolute_linear_strain.is_finite() && maximum_absolute_linear_strain > 0.0) {
            return Err(SolidGeometryEvolutionError::InvalidLaw {
                field: "maximum_absolute_linear_strain",
            });
        }
        if authority_identity.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(SolidGeometryEvolutionError::InvalidLaw {
                field: "authority_identity",
            });
        }
        let mut identity = DomainHasher::new(EULER_ISOTROPIC_FREE_EXPANSION_IDENTITY_DOMAIN);
        identity.update(authority_identity.as_bytes());
        identity.update(&maximum_absolute_linear_strain.to_bits().to_le_bytes());
        Ok(Self {
            maximum_absolute_linear_strain,
            authority_identity,
            identity: identity.finalize(),
        })
    }

    /// Maximum admitted magnitude of free linear strain.
    #[must_use]
    pub const fn maximum_absolute_linear_strain(self) -> f64 {
        self.maximum_absolute_linear_strain
    }

    /// Evidence or caller-declaration identity for isotropic free expansion.
    #[must_use]
    pub const fn authority_identity(self) -> ContentHash {
        self.authority_identity
    }

    /// Complete law identity.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }
}

/// A mass-conserving solid profile evolved by one admitted geometry law.
///
/// Consumers must still resolve temperature-dependent elasticity, contact,
/// damping, and optics at `phase_state` before advancing mechanics, sound, or
/// rendering. This value establishes geometry and mass only.
#[derive(Clone, Debug)]
pub struct ResolvedEvolvedSolidDiscProfile {
    /// Geometry and mass properties at the evolved equilibrium density.
    profile: ResolvedDiscProfile,
    /// Solid thermodynamic state that required the geometry update.
    phase_state: EquilibriumPhaseState,
    /// Uniform scale applied to every reference length.
    linear_scale: f64,
    /// Exact free-expansion law identity.
    expansion_law_identity: ContentHash,
    /// Identity binding reference specimen, phase state, law, and evolved profile.
    identity: ContentHash,
}

impl ResolvedEvolvedSolidDiscProfile {
    /// Geometry and mass properties at the evolved equilibrium density.
    #[must_use]
    pub const fn profile(&self) -> &ResolvedDiscProfile {
        &self.profile
    }

    /// Solid thermodynamic state that required the geometry update.
    #[must_use]
    pub const fn phase_state(&self) -> EquilibriumPhaseState {
        self.phase_state
    }

    /// Uniform scale applied to every reference length.
    #[must_use]
    pub const fn linear_scale(&self) -> f64 {
        self.linear_scale
    }

    /// Exact free-expansion law identity.
    #[must_use]
    pub const fn expansion_law_identity(&self) -> ContentHash {
        self.expansion_law_identity
    }

    /// Identity binding the admitted reference, phase state, law, and result.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Bind the evolved geometry to a complete isotropic solid state for contact.
    pub fn try_bind_isotropic_solid(
        &self,
        material: &IsotropicSolidStatePoint,
    ) -> Result<ResolvedMaterialDiscProfile, PhaseDiscBindingError> {
        require_matching_phase_material_state(
            self.phase_state,
            material.resolved().card_identity(),
            material.resolved().query_point(),
            material.density_kg_m3(),
        )?;
        let profile_identities = self.profile.content_identities();
        let mut identity = DomainHasher::new(EULER_MATERIAL_SPECIMEN_IDENTITY_DOMAIN);
        identity.update(profile_identities.profile.as_bytes());
        identity.update(profile_identities.mass_properties.as_bytes());
        identity.update(material.resolved().identity().as_bytes());
        Ok(ResolvedMaterialDiscProfile {
            profile: self.profile.clone(),
            material: material.clone(),
            identity: identity.finalize(),
        })
    }

    /// Bind the evolved geometry to minimal isotropic elasticity for modes and sound.
    pub fn try_bind_isotropic_elastic(
        &self,
        material: &IsotropicElasticStatePoint,
    ) -> Result<ResolvedElasticDiscProfile, PhaseDiscBindingError> {
        require_matching_phase_material_state(
            self.phase_state,
            material.resolved().card_identity(),
            material.resolved().query_point(),
            material.density_kg_m3(),
        )?;
        Ok(ResolvedElasticDiscProfile::bind(
            self.profile.clone(),
            TetElasticMaterial::from_resolved_elastic_state(material),
            material.resolved().card_identity(),
        ))
    }
}

/// One thermal boundary state coupled back to invariant specimen mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscThermalSample {
    /// Physical time [s].
    pub time_s: f64,
    /// Phase, density, required volume, and downstream geometry regime.
    pub mass_conserving_state: MassConservingDiscPhaseState,
    /// Convective power into the body [W].
    pub convection_into_body_w: f64,
    /// Net radiative power into the body [W].
    pub radiation_into_body_w: f64,
    /// Declared volumetric/internal power integrated by the thermal rung [W].
    pub internal_power_w: f64,
    /// Total power into the body [W].
    pub net_power_into_body_w: f64,
    /// Backward-Euler energy residual for the preceding interval [J].
    pub step_energy_residual_j: f64,
}

/// An isothermal enthalpy march proven to use the exact mechanics specimen.
///
/// This adapter does not deform a solid or evolve a liquid surface. Instead,
/// every boundary state names the next solver rung required before mechanics,
/// acoustics, or rendering may consume it.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDiscThermalMarch {
    /// Conservative Biot number that admitted the lumped thermal model.
    pub maximum_biot: f64,
    /// Ordered initial and accepted endpoint states.
    pub samples: Vec<DiscThermalSample>,
    /// Sum of absolute discrete thermal energy residuals [J].
    pub cumulative_absolute_energy_residual_j: f64,
    /// Original generic thermal-march identity.
    pub conduction_march_identity: ContentHash,
    /// Identity binding specimen, exact thermal geometry, and thermal march.
    pub identity: ContentHash,
}

impl ResolvedPhaseDiscProfile {
    /// Bind a generic lumped enthalpy march to this exact specimen.
    ///
    /// The body must use the profile's invariant mass, complete revolved
    /// surface area, `V/A` characteristic length, and phase curve bit for bit.
    /// This prevents a thermally convenient surrogate body from silently
    /// driving the mechanics, sound, or picture of a different specimen.
    pub fn bind_lumped_enthalpy_march(
        &self,
        body: &LumpedEnthalpyBody<'_>,
        march: &LumpedEnthalpyMarch,
        cx: &Cx<'_>,
    ) -> Result<ResolvedDiscThermalMarch, DiscThermalCouplingError> {
        let thermal_geometry = self
            .profile
            .thermal_geometry(cx)
            .map_err(DiscThermalCouplingError::ThermalGeometry)?;
        if march.body_identity() != body.identity() {
            return Err(DiscThermalCouplingError::BodyIdentityMismatch);
        }
        if body.phase_curve().identity() != self.phase_state.phase_curve_identity() {
            return Err(DiscThermalCouplingError::PhaseCurveMismatch);
        }
        for (field, actual, expected) in [
            ("mass_kg", body.mass_kg(), self.profile.mass_properties.mass),
            (
                "surface_area_m2",
                body.surface_area_m2(),
                thermal_geometry.surface_area_m2,
            ),
            (
                "characteristic_length_m",
                body.characteristic_length_m(),
                thermal_geometry.characteristic_length_m,
            ),
        ] {
            if actual.to_bits() != expected.to_bits() {
                return Err(DiscThermalCouplingError::SpecimenQuantityMismatch {
                    field,
                    expected,
                    actual,
                });
            }
        }
        let Some(initial) = march.samples().first() else {
            return Err(DiscThermalCouplingError::EmptyMarch);
        };
        if initial.phase_state.identity() != self.phase_state.identity() {
            return Err(DiscThermalCouplingError::InitialPhaseStateMismatch);
        }
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(march.samples().len())
            .map_err(|_| DiscThermalCouplingError::Capacity {
                requested: march.samples().len(),
            })?;
        for sample in march.samples() {
            let mass_conserving_state = self
                .mass_conserving_state(sample.phase_state)
                .map_err(DiscThermalCouplingError::PhaseBinding)?;
            samples.push(DiscThermalSample {
                time_s: sample.time_s,
                mass_conserving_state,
                convection_into_body_w: sample.convection_into_body_w,
                radiation_into_body_w: sample.radiation_into_body_w,
                internal_power_w: sample.internal_power_w,
                net_power_into_body_w: sample.net_power_into_body_w,
                step_energy_residual_j: sample.step_energy_residual_j,
            });
        }
        let mut identity = DomainHasher::new(EULER_SPECIMEN_THERMAL_MARCH_IDENTITY_DOMAIN);
        identity.update(self.identity.as_bytes());
        identity.update(thermal_geometry.identity.as_bytes());
        identity.update(march.identity().as_bytes());
        for sample in &samples {
            identity.update(sample.mass_conserving_state.identity.as_bytes());
        }
        Ok(ResolvedDiscThermalMarch {
            maximum_biot: march.maximum_biot(),
            samples,
            cumulative_absolute_energy_residual_j: march.cumulative_absolute_energy_residual_j(),
            conduction_march_identity: march.identity(),
            identity: identity.finalize(),
        })
    }

    /// Bind another state on the same phase curve to invariant specimen mass.
    ///
    /// This is the safe handoff from thermal transport to geometry evolution.
    /// It publishes the volume demanded by density and then explicitly names
    /// whether reference, thermomechanical-solid, or free-surface geometry is
    /// required; it never manufactures the missing geometry.
    pub fn mass_conserving_state(
        &self,
        phase_state: EquilibriumPhaseState,
    ) -> Result<MassConservingDiscPhaseState, PhaseDiscBindingError> {
        if phase_state.material_card_identity() != self.phase_state.material_card_identity() {
            return Err(PhaseDiscBindingError::MaterialCardMismatch);
        }
        if phase_state.phase_curve_identity() != self.phase_state.phase_curve_identity() {
            return Err(PhaseDiscBindingError::PhaseCurveMismatch);
        }
        let invariant_mass_kg = self.profile.mass_properties.mass;
        let reference_volume_m3 = self.profile.mass_properties.volume;
        let required_volume_m3 = if phase_state.identity() == self.phase_state.identity() {
            reference_volume_m3
        } else {
            invariant_mass_kg / phase_state.bulk_density_kg_m3()
        };
        let volume_ratio = required_volume_m3 / reference_volume_m3;
        if !(required_volume_m3.is_finite()
            && required_volume_m3 > 0.0
            && volume_ratio.is_finite()
            && volume_ratio > 0.0)
        {
            return Err(PhaseDiscBindingError::InvalidMassConservingVolume {
                invariant_mass_kg,
                bulk_density_kg_m3: phase_state.bulk_density_kg_m3(),
                reference_volume_m3,
            });
        }
        let geometry_regime = if phase_state.identity() == self.phase_state.identity() {
            DiscPhaseGeometryRegime::ReferenceGeometry
        } else if phase_state.phase() == SolidLiquidPhase::Solid {
            DiscPhaseGeometryRegime::SolidThermomechanicalUpdateRequired
        } else {
            DiscPhaseGeometryRegime::EvolvingFreeSurfaceRequired
        };
        let mut identity =
            DomainHasher::new("org.frankensim.fs-euler-disc-e2e.mass-conserving-phase-state.v1");
        identity.update(self.identity.as_bytes());
        identity.update(phase_state.identity().as_bytes());
        for value in [invariant_mass_kg, required_volume_m3, volume_ratio] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        identity.update(&[match geometry_regime {
            DiscPhaseGeometryRegime::ReferenceGeometry => 0,
            DiscPhaseGeometryRegime::SolidThermomechanicalUpdateRequired => 1,
            DiscPhaseGeometryRegime::EvolvingFreeSurfaceRequired => 2,
        }]);
        Ok(MassConservingDiscPhaseState {
            phase_state,
            invariant_mass_kg,
            required_volume_m3,
            volume_ratio,
            geometry_regime,
            identity: identity.finalize(),
        })
    }

    /// Evolve a still-solid body under a bounded free/isotropic/isothermal law.
    ///
    /// All profile lengths, including bores, fillets, chamfers, and taper
    /// dimensions, receive the same scale. The evolved profile is then
    /// re-resolved at the phase state's density, so geometry, support queries,
    /// mass, centroid, and inertia remain one coherent asset. Liquid states and
    /// changes outside the law's validity bound refuse.
    pub fn resolve_uniform_isotropic_free_expansion(
        &self,
        state: MassConservingDiscPhaseState,
        law: UniformIsotropicFreeExpansionLaw,
        cx: &Cx<'_>,
    ) -> Result<ResolvedEvolvedSolidDiscProfile, SolidGeometryEvolutionError> {
        let expected = self
            .mass_conserving_state(state.phase_state)
            .map_err(SolidGeometryEvolutionError::PhaseBinding)?;
        // The carrier is publicly constructible for transport/reporting.  A
        // matching retained hash alone is therefore not an admission token:
        // reject any field mutation that failed to recompute the identity.
        if expected != state {
            return Err(SolidGeometryEvolutionError::PhaseStateMismatch);
        }
        if state.phase_state.phase() != SolidLiquidPhase::Solid {
            return Err(SolidGeometryEvolutionError::EvolvingFreeSurfaceRequired {
                liquid_mass_fraction: state.phase_state.liquid_mass_fraction(),
            });
        }
        let linear_scale = state.volume_ratio.cbrt();
        if !(linear_scale.is_finite() && linear_scale > 0.0) {
            return Err(SolidGeometryEvolutionError::InvalidDerivedScale {
                volume_ratio: state.volume_ratio,
            });
        }
        let absolute_linear_strain = (linear_scale - 1.0).abs();
        if absolute_linear_strain > law.maximum_absolute_linear_strain() {
            return Err(SolidGeometryEvolutionError::LawValidityExceeded {
                absolute_linear_strain,
                maximum_absolute_linear_strain: law.maximum_absolute_linear_strain(),
            });
        }
        let scaled_spec = self
            .profile
            .spec
            .uniformly_scaled(linear_scale)
            .map_err(SolidGeometryEvolutionError::Profile)?;
        let profile = scaled_spec
            .resolve(state.phase_state.bulk_density_kg_m3(), cx)
            .map_err(SolidGeometryEvolutionError::Profile)?;
        let mass_tolerance = 64.0
            * f64::EPSILON
            * state
                .invariant_mass_kg
                .abs()
                .max(profile.mass_properties.mass.abs());
        let volume_tolerance = 64.0
            * f64::EPSILON
            * state
                .required_volume_m3
                .abs()
                .max(profile.mass_properties.volume.abs());
        if (profile.mass_properties.mass - state.invariant_mass_kg).abs() > mass_tolerance
            || (profile.mass_properties.volume - state.required_volume_m3).abs() > volume_tolerance
        {
            return Err(SolidGeometryEvolutionError::MassConservationFailure {
                expected_mass_kg: state.invariant_mass_kg,
                actual_mass_kg: profile.mass_properties.mass,
                expected_volume_m3: state.required_volume_m3,
                actual_volume_m3: profile.mass_properties.volume,
            });
        }
        let identities = profile.content_identities();
        let mut identity = DomainHasher::new(EULER_EVOLVED_SOLID_PROFILE_IDENTITY_DOMAIN);
        identity.update(self.identity.as_bytes());
        identity.update(state.identity.as_bytes());
        identity.update(law.identity().as_bytes());
        identity.update(identities.profile.as_bytes());
        identity.update(identities.mass_properties.as_bytes());
        identity.update(&linear_scale.to_bits().to_le_bytes());
        Ok(ResolvedEvolvedSolidDiscProfile {
            profile,
            phase_state: state.phase_state,
            linear_scale,
            expansion_law_identity: law.identity(),
            identity: identity.finalize(),
        })
    }

    /// Bind a fully solid phase state to the independently resolved elastic
    /// properties consumed by fixed-topology mechanics.
    ///
    /// # Errors
    /// Refuses any liquid fraction, different material card, different
    /// temperature, or different density. No stale solid law can therefore
    /// survive an admitted phase transition.
    pub fn try_bind_fixed_solid(
        &self,
        material: &IsotropicSolidStatePoint,
    ) -> Result<ResolvedMaterialDiscProfile, PhaseDiscBindingError> {
        if self.phase_state.phase() != SolidLiquidPhase::Solid {
            return Err(PhaseDiscBindingError::EvolvingPhaseRequired {
                liquid_mass_fraction: self.phase_state.liquid_mass_fraction(),
            });
        }
        let resolved = material.resolved();
        require_matching_phase_material_state(
            self.phase_state,
            resolved.card_identity(),
            resolved.query_point(),
            material.density_kg_m3(),
        )?;
        let mut identity = DomainHasher::new(EULER_MATERIAL_SPECIMEN_IDENTITY_DOMAIN);
        let profile_identities = self.profile.content_identities();
        identity.update(profile_identities.profile.as_bytes());
        identity.update(profile_identities.mass_properties.as_bytes());
        identity.update(material.resolved().identity().as_bytes());
        Ok(ResolvedMaterialDiscProfile {
            profile: self.profile.clone(),
            material: material.clone(),
            identity: identity.finalize(),
        })
    }
}

fn require_matching_phase_material_state(
    phase_state: EquilibriumPhaseState,
    material_card_identity: ContentHash,
    query_point: &[(String, f64)],
    density_kg_m3: f64,
) -> Result<(), PhaseDiscBindingError> {
    if material_card_identity != phase_state.material_card_identity() {
        return Err(PhaseDiscBindingError::MaterialCardMismatch);
    }
    let temperature_k = query_point
        .binary_search_by(|(axis, _)| axis.as_str().cmp("T"))
        .ok()
        .and_then(|index| query_point.get(index))
        .map(|(_, temperature_k)| *temperature_k)
        .ok_or(PhaseDiscBindingError::MissingTemperatureCoordinate)?;
    if temperature_k.to_bits() != phase_state.temperature_k().to_bits() {
        return Err(PhaseDiscBindingError::TemperatureMismatch {
            phase_temperature_k: phase_state.temperature_k(),
            mechanical_temperature_k: temperature_k,
        });
    }
    if density_kg_m3.to_bits() != phase_state.bulk_density_kg_m3().to_bits() {
        return Err(PhaseDiscBindingError::DensityMismatch {
            phase_density_kg_m3: phase_state.bulk_density_kg_m3(),
            mechanical_density_kg_m3: density_kg_m3,
        });
    }
    Ok(())
}

/// Refusal from joining a phase state to fixed-solid mechanics.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseDiscBindingError {
    /// A nonzero liquid fraction requires phase-aware evolving geometry.
    EvolvingPhaseRequired {
        /// Admitted liquid mass fraction that triggered the refusal.
        liquid_mass_fraction: f64,
    },
    /// Phase and mechanical properties came from different material cards.
    MaterialCardMismatch,
    /// A thermal state came from a different enthalpy/phase curve.
    PhaseCurveMismatch,
    /// Density and invariant mass did not produce a finite positive volume.
    InvalidMassConservingVolume {
        /// Reference specimen mass [kg].
        invariant_mass_kg: f64,
        /// New equilibrium bulk density [kg/m3].
        bulk_density_kg_m3: f64,
        /// Reference specimen volume [m3].
        reference_volume_m3: f64,
    },
    /// The elastic state omitted the canonical absolute-temperature axis.
    MissingTemperatureCoordinate,
    /// Phase and elastic properties were resolved at different temperatures.
    TemperatureMismatch {
        /// Temperature carried by the enthalpy/phase state [K].
        phase_temperature_k: f64,
        /// Temperature carried by the mechanical state [K].
        mechanical_temperature_k: f64,
    },
    /// Phase and elastic properties disagree on homogeneous bulk density.
    DensityMismatch {
        /// Density carried by the enthalpy/phase state [kg/m3].
        phase_density_kg_m3: f64,
        /// Density carried by the mechanical state [kg/m3].
        mechanical_density_kg_m3: f64,
    },
}

/// Refusal while binding a generic thermal march to one exact disc profile.
#[derive(Clone, Debug, PartialEq)]
pub enum DiscThermalCouplingError {
    /// Recomputing whole-boundary thermal geometry refused.
    ThermalGeometry(DiscProfileError),
    /// The march was not produced by the supplied body.
    BodyIdentityMismatch,
    /// The body's phase curve differs from the specimen's admitted curve.
    PhaseCurveMismatch,
    /// A body mass or geometry measure differs from the mechanics specimen.
    SpecimenQuantityMismatch {
        /// Quantity that failed the exact binding.
        field: &'static str,
        /// Exact profile-derived value.
        expected: f64,
        /// Value supplied to the thermal body.
        actual: f64,
    },
    /// A successful march must retain its initial state.
    EmptyMarch,
    /// The march starts from a different thermodynamic state than the profile.
    InitialPhaseStateMismatch,
    /// One phase state could not be conservatively joined to the specimen.
    PhaseBinding(PhaseDiscBindingError),
    /// Retaining the bounded coupled samples failed.
    Capacity {
        /// Requested number of thermal boundary states.
        requested: usize,
    },
}

/// Refusal while evolving a still-solid reference profile.
#[derive(Clone, Debug, PartialEq)]
pub enum SolidGeometryEvolutionError {
    /// Free-expansion law metadata or validity bound was malformed.
    InvalidLaw {
        /// Invalid field.
        field: &'static str,
    },
    /// The state was not derived from this exact reference specimen.
    PhaseStateMismatch,
    /// Rebinding the phase state to the reference specimen refused.
    PhaseBinding(PhaseDiscBindingError),
    /// A liquid fraction requires free-surface evolution instead of similarity scaling.
    EvolvingFreeSurfaceRequired {
        /// Equilibrium liquid mass fraction that triggered escalation.
        liquid_mass_fraction: f64,
    },
    /// Density-derived volume ratio did not produce a usable linear scale.
    InvalidDerivedScale {
        /// Required evolved-to-reference volume ratio.
        volume_ratio: f64,
    },
    /// Required free strain exceeds the supplied constitutive law's validity.
    LawValidityExceeded {
        /// Required magnitude of `linear_scale - 1`.
        absolute_linear_strain: f64,
        /// Maximum magnitude admitted by the law.
        maximum_absolute_linear_strain: f64,
    },
    /// Constructing or resolving the scaled profile refused.
    Profile(DiscProfileError),
    /// The scaled profile failed the independent mass/volume closure check.
    MassConservationFailure {
        /// Invariant reference mass [kg].
        expected_mass_kg: f64,
        /// Reintegrated evolved mass [kg].
        actual_mass_kg: f64,
        /// Density-derived required volume [m3].
        expected_volume_m3: f64,
        /// Reintegrated evolved volume [m3].
        actual_volume_m3: f64,
    },
}

impl fmt::Display for SolidGeometryEvolutionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for SolidGeometryEvolutionError {}

impl fmt::Display for DiscThermalCouplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for DiscThermalCouplingError {}

impl fmt::Display for PhaseDiscBindingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for PhaseDiscBindingError {}

/// Strong content identities used to bind a trajectory to its resolved asset.
///
/// These roots complement the compact diagnostic [`AxisymmetricIdentity`].
/// They hash the complete canonical line/arc input rather than promoting its
/// 64-bit cache fingerprint into a durable content address.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResolvedDiscProfileIdentities {
    /// Exact ordered canonical meridian and construction-schema identity.
    pub chart: ContentHash,
    /// Exact chart plus homogeneous density identity.
    pub profile: ContentHash,
    /// Resolved production mass/centroid/inertia identity.
    pub mass_properties: ContentHash,
}

impl ResolvedDiscProfile {
    /// Derive whole-boundary thermal measures from the mechanics geometry.
    ///
    /// The chart is revalidated and integrated under `cx`; no render mesh,
    /// nominal cylinder formula, or separately supplied area participates.
    pub fn thermal_geometry(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ResolvedDiscThermalGeometry, DiscProfileError> {
        let surface_area_m2 = self
            .chart
            .surface_area(cx)
            .map_err(DiscProfileError::SurfaceArea)?
            .area;
        let volume_m3 = self.mass_properties.volume;
        let characteristic_length_m = volume_m3 / surface_area_m2;
        if !(characteristic_length_m.is_finite() && characteristic_length_m > 0.0) {
            return Err(DiscProfileError::InvalidDerivedThermalGeometry {
                volume_m3,
                surface_area_m2,
            });
        }
        let identities = self.content_identities();
        let mut identity = DomainHasher::new(EULER_SPECIMEN_THERMAL_GEOMETRY_IDENTITY_DOMAIN);
        identity.update(identities.chart.as_bytes());
        for value in [volume_m3, surface_area_m2, characteristic_length_m] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        Ok(ResolvedDiscThermalGeometry {
            surface_area_m2,
            volume_m3,
            characteristic_length_m,
            identity: identity.finalize(),
        })
    }

    /// Compute the versioned strong identities used by trajectory metadata and
    /// visualization asset admission.
    #[must_use]
    pub fn content_identities(&self) -> ResolvedDiscProfileIdentities {
        let mut chart = DomainHasher::new(EULER_SPECIMEN_CHART_IDENTITY_DOMAIN);
        let certificate = self.chart.construction_certificate();
        chart.update(&certificate.schema_version.to_le_bytes());
        chart.update(
            &u64::try_from(self.chart.segments().len())
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        for segment in self.chart.segments() {
            match *segment {
                MeridianSegment::Line { start, end } => {
                    chart.update(&[0]);
                    hash_meridian_point(&mut chart, start);
                    hash_meridian_point(&mut chart, end);
                }
                MeridianSegment::Arc {
                    start,
                    end,
                    center,
                    clockwise,
                } => {
                    chart.update(&[1]);
                    hash_meridian_point(&mut chart, start);
                    hash_meridian_point(&mut chart, end);
                    hash_meridian_point(&mut chart, center);
                    chart.update(&[u8::from(clockwise)]);
                }
            }
        }
        let chart = chart.finalize();

        let mut profile = DomainHasher::new(EULER_SPECIMEN_PROFILE_IDENTITY_DOMAIN);
        profile.update(chart.as_bytes());
        profile.update(&self.density_kg_per_m3.to_bits().to_le_bytes());
        let profile = profile.finalize();

        let mut mass = DomainHasher::new(EULER_SPECIMEN_MASS_IDENTITY_DOMAIN);
        mass.update(profile.as_bytes());
        for value in [
            self.mass_properties.volume,
            self.mass_properties.mass,
            self.mass_properties.center_of_mass.x,
            self.mass_properties.center_of_mass.y,
            self.mass_properties.center_of_mass.z,
            self.mass_properties.principal_inertia.transverse,
            self.mass_properties.principal_inertia.axial,
            self.mass_properties.origin_inertia.transverse,
            self.mass_properties.origin_inertia.axial,
        ] {
            mass.update(&value.to_bits().to_le_bytes());
        }

        ResolvedDiscProfileIdentities {
            chart,
            profile,
            mass_properties: mass.finalize(),
        }
    }
}

fn hash_meridian_point(hasher: &mut DomainHasher, point: MeridianPoint) {
    hasher.update(&point.radius.to_bits().to_le_bytes());
    hasher.update(&point.axial.to_bits().to_le_bytes());
}

/// Refusal from resolving a profile specification.
#[derive(Clone, Debug, PartialEq)]
pub enum DiscProfileError {
    /// A named profile parameter was non-finite or outside its documented domain.
    InvalidParameter {
        /// Stable name of the refused parameter.
        field: &'static str,
        /// Refused numeric value.
        value: f64,
    },
    /// A parameter relationship does not describe the named profile family.
    InvalidRelationship {
        /// Stable description of the violated relationship.
        detail: &'static str,
    },
    /// The generic line/arc chart refused the constructed meridian.
    Geometry(AxisymmetricError),
    /// The exact line/arc mass integration refused to publish properties.
    Mass(AxisymmetricMassError),
    /// The exact line/arc surface integration refused to publish an area.
    SurfaceArea(AxisymmetricSurfaceAreaError),
    /// Finite positive volume and area did not produce a usable `V/A` length.
    InvalidDerivedThermalGeometry {
        /// Resolved enclosed volume [m3].
        volume_m3: f64,
        /// Resolved complete surface area [m2].
        surface_area_m2: f64,
    },
}

impl fmt::Display for DiscProfileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter { field, value } => {
                write!(
                    formatter,
                    "invalid Euler-disc profile parameter {field}={value}"
                )
            }
            Self::InvalidRelationship { detail } => {
                write!(
                    formatter,
                    "invalid Euler-disc profile relationship: {detail}"
                )
            }
            Self::Geometry(error) => {
                write!(formatter, "Euler-disc profile geometry refused: {error}")
            }
            Self::Mass(error) => write!(formatter, "Euler-disc profile mass refused: {error}"),
            Self::SurfaceArea(error) => {
                write!(
                    formatter,
                    "Euler-disc profile surface area refused: {error}"
                )
            }
            Self::InvalidDerivedThermalGeometry {
                volume_m3,
                surface_area_m2,
            } => write!(
                formatter,
                "Euler-disc profile produced invalid thermal geometry V={volume_m3} m3, A={surface_area_m2} m2"
            ),
        }
    }
}

impl std::error::Error for DiscProfileError {}

impl DiscProfileSpec {
    /// Scale every physical length in the profile by one positive factor.
    ///
    /// This operation preserves the exact profile family and all dimensionless
    /// proportions. It does not select a thermomechanical law; callers must
    /// obtain the factor from an admitted constitutive/phase-state coupling.
    pub fn uniformly_scaled(self, linear_scale: f64) -> Result<Self, DiscProfileError> {
        positive("linear_scale", linear_scale)?;
        let scale_edge = |edge| match edge {
            SquatDiscEdgeTreatment::Sharp => SquatDiscEdgeTreatment::Sharp,
            SquatDiscEdgeTreatment::CircularFillet { radius } => {
                SquatDiscEdgeTreatment::CircularFillet {
                    radius: radius * linear_scale,
                }
            }
        };
        let scaled = match self {
            Self::SolidCylinder {
                outer_radius_m,
                thickness_m,
                edge_treatment,
            } => Self::SolidCylinder {
                outer_radius_m: outer_radius_m * linear_scale,
                thickness_m: thickness_m * linear_scale,
                edge_treatment: scale_edge(edge_treatment),
            },
            Self::AnnularCylinder {
                outer_radius_m,
                inner_radius_m,
                thickness_m,
            } => Self::AnnularCylinder {
                outer_radius_m: outer_radius_m * linear_scale,
                inner_radius_m: inner_radius_m * linear_scale,
                thickness_m: thickness_m * linear_scale,
            },
            Self::OuterFilletedAnnularCylinder {
                outer_radius_m,
                inner_radius_m,
                thickness_m,
                outer_fillet_radius_m,
            } => Self::OuterFilletedAnnularCylinder {
                outer_radius_m: outer_radius_m * linear_scale,
                inner_radius_m: inner_radius_m * linear_scale,
                thickness_m: thickness_m * linear_scale,
                outer_fillet_radius_m: outer_fillet_radius_m * linear_scale,
            },
            Self::SymmetricTapered {
                outer_radius_m,
                face_radius_m,
                thickness_m,
            } => Self::SymmetricTapered {
                outer_radius_m: outer_radius_m * linear_scale,
                face_radius_m: face_radius_m * linear_scale,
                thickness_m: thickness_m * linear_scale,
            },
            Self::ChamferedCylinder {
                outer_radius_m,
                thickness_m,
                chamfer_radial_m,
                chamfer_axial_m,
            } => Self::ChamferedCylinder {
                outer_radius_m: outer_radius_m * linear_scale,
                thickness_m: thickness_m * linear_scale,
                chamfer_radial_m: chamfer_radial_m * linear_scale,
                chamfer_axial_m: chamfer_axial_m * linear_scale,
            },
        };
        // Refuse overflow and any relationship lost to binary64 scaling before
        // returning a public specification.
        scaled.chart_and_dimensions()?;
        Ok(scaled)
    }

    /// Resolve reference geometry and mass at one equilibrium phase state.
    ///
    /// This is phase-agnostic ingress: it admits solid, mushy, and liquid
    /// states without pretending the reference shape is the evolved free
    /// surface. Downstream solvers select an appropriate constitutive and
    /// geometry-evolution rung from `phase_state`.
    pub fn resolve_with_phase_state(
        self,
        phase_state: EquilibriumPhaseState,
        cx: &Cx<'_>,
    ) -> Result<ResolvedPhaseDiscProfile, DiscProfileError> {
        let profile = self.resolve(phase_state.bulk_density_kg_m3(), cx)?;
        let mut identity = DomainHasher::new("org.frankensim.fs-euler-disc-e2e.phase-specimen.v1");
        let profile_identities = profile.content_identities();
        identity.update(profile_identities.profile.as_bytes());
        identity.update(profile_identities.mass_properties.as_bytes());
        identity.update(phase_state.identity().as_bytes());
        Ok(ResolvedPhaseDiscProfile {
            profile,
            phase_state,
            identity: identity.finalize(),
        })
    }

    /// Resolve this geometry using the density from one admitted material
    /// state, while retaining that complete state for downstream contact,
    /// structural, acoustic, thermal, optical, and provenance consumers.
    pub fn resolve_with_material_state(
        self,
        material: &IsotropicSolidStatePoint,
        cx: &Cx<'_>,
    ) -> Result<ResolvedMaterialDiscProfile, DiscProfileError> {
        let profile = self.resolve(material.density_kg_m3(), cx)?;
        let profile_identities = profile.content_identities();
        let mut identity = DomainHasher::new(EULER_MATERIAL_SPECIMEN_IDENTITY_DOMAIN);
        identity.update(profile_identities.profile.as_bytes());
        identity.update(profile_identities.mass_properties.as_bytes());
        identity.update(material.resolved().identity().as_bytes());
        Ok(ResolvedMaterialDiscProfile {
            profile,
            material: material.clone(),
            identity: identity.finalize(),
        })
    }

    /// Resolve geometry and mass from the minimal isotropic tangent-elastic
    /// state used by modal vibration and acoustics.
    pub fn resolve_with_isotropic_elastic_state(
        self,
        material: &IsotropicElasticStatePoint,
        cx: &Cx<'_>,
    ) -> Result<ResolvedElasticDiscProfile, DiscProfileError> {
        let profile = self.resolve(material.density_kg_m3(), cx)?;
        Ok(ResolvedElasticDiscProfile::bind(
            profile,
            TetElasticMaterial::from_resolved_elastic_state(material),
            material.resolved().card_identity(),
        ))
    }

    /// Resolve geometry and mass from an oriented orthotropic tangent-elastic
    /// state. `principal_to_world` maps material axes into the specimen-local
    /// frame; its nonzero identity must name the texture/orientation evidence.
    pub fn resolve_with_orthotropic_elastic_state(
        self,
        material: &OrthotropicElasticStatePoint,
        principal_to_world: [[f64; 3]; 3],
        orientation_identity: ContentHash,
        cx: &Cx<'_>,
    ) -> Result<ResolvedElasticDiscProfile, ElasticDiscProfileError> {
        let profile = self.resolve(material.density_kg_m3(), cx)?;
        let elastic_material = TetElasticMaterial::try_from_resolved_orthotropic_state(
            material,
            principal_to_world,
            orientation_identity,
        )?;
        Ok(ResolvedElasticDiscProfile::bind(
            profile,
            elastic_material,
            material.resolved().card_identity(),
        ))
    }

    /// Resolve the specification into a validated profile and matching mass
    /// properties.  The same `Cx` controls chart mass integration and is
    /// retained by callers for later support queries.
    pub fn resolve(
        self,
        density_kg_per_m3: f64,
        cx: &Cx<'_>,
    ) -> Result<ResolvedDiscProfile, DiscProfileError> {
        if !density_kg_per_m3.is_finite() || density_kg_per_m3 <= 0.0 {
            return Err(DiscProfileError::InvalidParameter {
                field: "density_kg_per_m3",
                value: density_kg_per_m3,
            });
        }
        let (chart, dimensions) = self.chart_and_dimensions()?;
        let identity = chart.construction_certificate().identity;
        let mass_properties = chart
            .mass_properties(density_kg_per_m3, cx)
            .map_err(DiscProfileError::Mass)?;
        Ok(ResolvedDiscProfile {
            spec: self,
            chart,
            density_kg_per_m3,
            identity,
            dimensions,
            mass_properties,
        })
    }

    fn chart_and_dimensions(
        self,
    ) -> Result<(AxisymmetricChart, DiscProfileDimensions), DiscProfileError> {
        match self {
            Self::SolidCylinder {
                outer_radius_m,
                thickness_m,
                edge_treatment,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("thickness_m", thickness_m)?;
                let chart =
                    AxisymmetricChart::squat_disc(outer_radius_m, thickness_m, edge_treatment)
                        .map_err(DiscProfileError::Geometry)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::AnnularCylinder {
                outer_radius_m,
                inner_radius_m,
                thickness_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("inner_radius_m", inner_radius_m)?;
                positive("thickness_m", thickness_m)?;
                if inner_radius_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "annular inner_radius_m must be smaller than outer_radius_m",
                    });
                }
                let half = 0.5 * thickness_m;
                let chart = chart_from_segments(vec![
                    line(inner_radius_m, -half, outer_radius_m, -half),
                    line(outer_radius_m, -half, outer_radius_m, half),
                    line(outer_radius_m, half, inner_radius_m, half),
                    line(inner_radius_m, half, inner_radius_m, -half),
                ])?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::OuterFilletedAnnularCylinder {
                outer_radius_m,
                inner_radius_m,
                thickness_m,
                outer_fillet_radius_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("inner_radius_m", inner_radius_m)?;
                positive("thickness_m", thickness_m)?;
                positive("outer_fillet_radius_m", outer_fillet_radius_m)?;
                if inner_radius_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "annular inner_radius_m must be smaller than outer_radius_m",
                    });
                }
                let maximum_fillet_radius =
                    (outer_radius_m - inner_radius_m).min(0.5 * thickness_m);
                if outer_fillet_radius_m > maximum_fillet_radius {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "outer_fillet_radius_m must not exceed both annular cap span and thickness_m / 2",
                    });
                }
                let chart = AxisymmetricChart::annular_disc_outer_fillets(
                    outer_radius_m,
                    inner_radius_m,
                    thickness_m,
                    outer_fillet_radius_m,
                )
                .map_err(DiscProfileError::Geometry)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::SymmetricTapered {
                outer_radius_m,
                face_radius_m,
                thickness_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                nonnegative("face_radius_m", face_radius_m)?;
                positive("thickness_m", thickness_m)?;
                if face_radius_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "symmetric tapered face_radius_m must be smaller than outer_radius_m",
                    });
                }
                let half = 0.5 * thickness_m;
                let segments = if face_radius_m == 0.0 {
                    vec![
                        line(0.0, -half, outer_radius_m, 0.0),
                        line(outer_radius_m, 0.0, 0.0, half),
                        line(0.0, half, 0.0, -half),
                    ]
                } else {
                    vec![
                        line(0.0, -half, face_radius_m, -half),
                        line(face_radius_m, -half, outer_radius_m, 0.0),
                        line(outer_radius_m, 0.0, face_radius_m, half),
                        line(face_radius_m, half, 0.0, half),
                        line(0.0, half, 0.0, -half),
                    ]
                };
                let chart = chart_from_segments(segments)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::ChamferedCylinder {
                outer_radius_m,
                thickness_m,
                chamfer_radial_m,
                chamfer_axial_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("thickness_m", thickness_m)?;
                positive("chamfer_radial_m", chamfer_radial_m)?;
                positive("chamfer_axial_m", chamfer_axial_m)?;
                if chamfer_radial_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "chamfer_radial_m must be smaller than outer_radius_m",
                    });
                }
                let half = 0.5 * thickness_m;
                if chamfer_axial_m > half {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "chamfer_axial_m must not exceed thickness_m / 2",
                    });
                }
                let lower_cap_radius = outer_radius_m - chamfer_radial_m;
                let lower_outer_z = -half + chamfer_axial_m;
                let upper_outer_z = half - chamfer_axial_m;
                let mut segments = vec![
                    line(0.0, -half, lower_cap_radius, -half),
                    line(lower_cap_radius, -half, outer_radius_m, lower_outer_z),
                ];
                if chamfer_axial_m < half {
                    segments.push(line(
                        outer_radius_m,
                        lower_outer_z,
                        outer_radius_m,
                        upper_outer_z,
                    ));
                }
                segments.extend([
                    line(outer_radius_m, upper_outer_z, lower_cap_radius, half),
                    line(lower_cap_radius, half, 0.0, half),
                    line(0.0, half, 0.0, -half),
                ]);
                let chart = chart_from_segments(segments)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
        }
    }
}

fn point(radius: f64, axial: f64) -> MeridianPoint {
    MeridianPoint::new(radius, axial)
}

fn line(start_radius: f64, start_axial: f64, end_radius: f64, end_axial: f64) -> MeridianSegment {
    MeridianSegment::Line {
        start: point(start_radius, start_axial),
        end: point(end_radius, end_axial),
    }
}

fn chart_from_segments(
    segments: Vec<MeridianSegment>,
) -> Result<AxisymmetricChart, DiscProfileError> {
    AxisymmetricChart::try_new(segments).map_err(DiscProfileError::Geometry)
}

fn positive(field: &'static str, value: f64) -> Result<(), DiscProfileError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(DiscProfileError::InvalidParameter { field, value })
    }
}

fn nonnegative(field: &'static str, value: f64) -> Result<(), DiscProfileError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(DiscProfileError::InvalidParameter { field, value })
    }
}
