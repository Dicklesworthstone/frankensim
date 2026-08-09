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
use fs_exec::Cx;
use fs_material::{
    phase::{EquilibriumPhaseState, SolidLiquidPhase},
    state_point::{
        IsotropicElasticStatePoint, IsotropicSolidStatePoint, OrthotropicElasticStatePoint,
    },
};
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricError, AxisymmetricIdentity, AxisymmetricMassError,
    AxisymmetricMassProperties, MeridianPoint, MeridianSegment, SquatDiscEdgeTreatment,
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

impl ResolvedPhaseDiscProfile {
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
        if resolved.card_identity() != self.phase_state.material_card_identity() {
            return Err(PhaseDiscBindingError::MaterialCardMismatch);
        }
        let temperature_k = resolved
            .query_point()
            .binary_search_by(|(axis, _)| axis.as_str().cmp("T"))
            .ok()
            .map(|index| resolved.query_point()[index].1)
            .ok_or(PhaseDiscBindingError::MissingTemperatureCoordinate)?;
        if temperature_k.to_bits() != self.phase_state.temperature_k().to_bits() {
            return Err(PhaseDiscBindingError::TemperatureMismatch {
                phase_temperature_k: self.phase_state.temperature_k(),
                mechanical_temperature_k: temperature_k,
            });
        }
        if material.density_kg_m3().to_bits() != self.phase_state.bulk_density_kg_m3().to_bits() {
            return Err(PhaseDiscBindingError::DensityMismatch {
                phase_density_kg_m3: self.phase_state.bulk_density_kg_m3(),
                mechanical_density_kg_m3: material.density_kg_m3(),
            });
        }
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
        }
    }
}

impl std::error::Error for DiscProfileError {}

impl DiscProfileSpec {
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
