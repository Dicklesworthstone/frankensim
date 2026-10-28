//! Low-Tabor DMT-M sphere contact with an explicit finite-range Dugdale tail.
//!
//! The contact branch retains Hertz geometry and subtracts the constant force
//! 2 pi R W. This is the DMT-M limiting approximation described by Maugis,
//! J. Colloid Interface Sci. 150 (1992), DOI 10.1016/0021-9797(92)90285-T,
//! following Derjaguin, Muller and Toporov (1975),
//! DOI 10.1016/0021-9797(75)90018-1. It is not the full Maugis crossover.
//!
//! The separately selected noncontact continuation integrates a constant
//! attractive traction W/z0 over the rigid parabolic gap g+r^2/(2R).
//! It therefore gives F=-2 pi R W (1-g/z0) for 0<g<z0. This deliberate
//! short-range Derjaguin approximation neglects attraction-induced deformation.
//! Its potential and force join continuously to the contact branch.
//! Tabor <= 0.1 is a model-selection ceiling, not an accuracy certificate.
//! Roughness, viscoelasticity, fluid bridges, mixed-mode fracture and changes
//! of the frozen source conditions require other owners.

#[path = "dmt/energy.rs"]
mod energy;

use core::f64::consts::PI;

use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::{InputAuthority, TriboError, adhesive::AdhesiveContactState};

use super::PointNormalPort;
use crate::interface_binding::{ADHESION_ENERGY_PROPERTY, BoundNormalInterfaceState};

/// Explicit geometric and low-Tabor applicability envelope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtEnvelope {
    /// Maximum positive signed approach, m.
    pub maximum_indentation_m: f64,
    /// Available half-space depth, m.
    pub half_space_depth_m: f64,
    /// Smallest relevant layer thickness, m.
    pub layer_thickness_m: f64,
    /// Upper bound on sqrt(a^2+2 R z0)/R, including the attractive footprint.
    pub maximum_interaction_radius_to_radius: f64,
    /// Upper bound on the same footprint/depth.
    pub maximum_interaction_radius_to_depth: f64,
    /// Upper bound on the same footprint/layer thickness.
    pub maximum_interaction_radius_to_layer: f64,
    /// Upper bound on Hertz nominal strain a/R.
    pub maximum_nominal_strain: f64,
    /// Upper bound on Hertz peak compressive pressure / bulk yield stress.
    pub maximum_repulsive_pressure_to_yield: f64,
    /// Upper bound on the declared cohesive traction (W/z0) / bulk yield stress.
    pub maximum_adhesive_stress_to_yield: f64,
    /// Upper bound on z0/R, no greater than 0.01.
    pub maximum_interaction_range_to_radius: f64,
    /// Upper bound on (R W^2 / (E*^2 z0^3))^(1/3), no greater than 0.1.
    pub maximum_tabor_parameter: f64,
}

/// Geometric regimes; attraction alone is not solid contact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DmtRegime {
    /// Gap is at least the explicit interaction range.
    Separated,
    /// Finite gap with attractive force but zero real contact area.
    AttractiveGap,
    /// Nonnegative indentation, including zero-area geometric touch.
    Contact,
}

/// One conservative local response with original normal-source identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtResponse {
    normal_state_identity: ContentHash,
    /// Geometric regime, independent of the sign of net normal force.
    pub regime: DmtRegime,
    /// Signed approach, m; negative is a finite geometric gap.
    pub indentation_m: f64,
    /// Hertz contact radius, m; zero throughout noncontact.
    pub contact_radius_m: f64,
    /// Actual solid contact area, m2, excluding the attractive annulus.
    pub contact_area_m2: f64,
    /// Hertz repulsive contribution, N.
    pub repulsive_normal_force_n: f64,
    /// Signed resultant force, N; negative denotes attraction.
    pub normal_force_n: f64,
    /// dF/d(delta), N/m, absent at the two piecewise-law corners.
    pub tangent_n_per_m: Option<f64>,
    /// Left and right force derivatives. Equal on every smooth branch.
    pub one_sided_tangents_n_per_m: [f64; 2],
    /// Elastic plus interaction potential, J, zero outside interaction range.
    pub reversible_energy_j: f64,
    /// Model estimate; explicitly synthetic inputs remain synthetic.
    pub authority: InputAuthority,
}

impl DmtResponse {
    /// Exact original normal-state identity for source-bound consumers.
    #[must_use]
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_state_identity
    }

    /// Supply real contact area and repulsive force to the adhesive shear owner.
    /// A noncontact attractive force deliberately supplies zero shear area/load.
    pub fn friction_contact_state(&self) -> Result<AdhesiveContactState, TriboError> {
        AdhesiveContactState::new(self.contact_area_m2, self.repulsive_normal_force_n)
    }

    /// Signed action/reaction at a common origin on a smooth branch.
    /// The unit normal points from the reacting to the acting body.
    /// At touch/cutoff, the consumer must select a one-sided event treatment.
    pub fn point_port(
        &self,
        normal: [f64; 3],
        moment_arm_m: [f64; 3],
    ) -> Result<PointNormalPort, DmtError> {
        let tangent = self.tangent_n_per_m.ok_or(DmtError::CornerTangent)?;
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !length.is_finite()
            || (length - 1.0).abs() > 1e-12
            || moment_arm_m.iter().any(|x| !x.is_finite())
        {
            return Err(DmtError::InvalidInput {
                field: "DMT port frame",
            });
        }
        let f = normal.map(|n| n * self.normal_force_n);
        let r = moment_arm_m;
        let m = [
            r[1] * f[2] - r[2] * f[1],
            r[2] * f[0] - r[0] * f[2],
            r[0] * f[1] - r[1] * f[0],
        ];
        if f.iter().chain(m.iter()).any(|x| !x.is_finite()) {
            return Err(DmtError::Unrepresentable {
                field: "DMT port wrench",
            });
        }
        Ok(PointNormalPort {
            action_force_n: f,
            action_moment_n_m: m,
            reaction_force_n: f.map(|x| -x),
            reaction_moment_n_m: m.map(|x| -x),
            residual_force_n: [0.0; 3],
            tangent_n_per_m: tangent,
            dissipated_power_w: 0.0,
            irreversible_work_j: 0.0,
        })
    }
}

/// Explicit refusal instead of using DMT outside its selected regime.
#[derive(Clone, Debug, PartialEq)]
pub enum DmtError {
    /// Invalid explicit parameter or kinematic input.
    InvalidInput { field: &'static str },
    /// Frozen elastic/adhesion properties have unowned response coordinates.
    NonstaticSource,
    /// An explicit physical applicability limit is exceeded.
    OutsideApplicability {
        ratio: &'static str,
        value: f64,
        limit: f64,
    },
    /// A requested applied load is below the contacting branch's pull-off.
    NoContactEquilibrium { load_n: f64, pull_off_n: f64 },
    /// A smooth tangent was requested at touch or the finite-range cutoff.
    CornerTangent,
    /// A derived number cannot be represented as a finite SI scalar.
    Unrepresentable { field: &'static str },
}
impl core::fmt::Display for DmtError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { field } => write!(f, "invalid DMT input: {field}"),
            Self::NonstaticSource => write!(
                f,
                "DMT requires static sources at fixed external conditions"
            ),
            Self::OutsideApplicability {
                ratio,
                value,
                limit,
            } => write!(f, "DMT {ratio}={value} exceeds limit {limit}"),
            Self::NoContactEquilibrium { load_n, pull_off_n } => write!(
                f,
                "DMT load {load_n} is below contact pull-off {pull_off_n}"
            ),
            Self::CornerTangent => write!(f, "DMT touch/cutoff needs one-sided event treatment"),
            Self::Unrepresentable { field } => write!(f, "unrepresentable DMT quantity: {field}"),
        }
    }
}
impl std::error::Error for DmtError {}

/// Explicit low-Tabor model, sharing the original immutable source owner.
/// Construction screens the operating envelope; evaluation is allocation-free O(1).
#[derive(Clone, Debug, PartialEq)]
pub struct SourceDmtSphereContact {
    normal: BoundNormalInterfaceState,
    radius_m: f64,
    interaction_range_m: f64,
    envelope: DmtEnvelope,
    identity: ContentHash,
    hertz_coefficient: f64,
    pull_off_magnitude_n: f64,
    tail_stiffness_n_per_m: f64,
    tabor_parameter: f64,
}

impl SourceDmtSphereContact {
    /// Select the DMT-M contact branch and explicit constant-traction noncontact
    /// tail. W is read from the existing ordered interface, not guessed from
    /// material names. z0 is an authored or measured cohesive interaction range.
    pub fn new(
        normal: BoundNormalInterfaceState,
        effective_radius_m: f64,
        interaction_range_m: f64,
        envelope: DmtEnvelope,
    ) -> Result<Self, DmtError> {
        let w = normal.adhesion_energy_j_per_m2();
        for (field, value) in [
            ("work of adhesion", w),
            ("effective radius", effective_radius_m),
            ("interaction range", interaction_range_m),
            ("maximum indentation", envelope.maximum_indentation_m),
            ("half-space depth", envelope.half_space_depth_m),
            ("layer thickness", envelope.layer_thickness_m),
        ] {
            positive(value, field)?;
        }
        for (field, value) in [
            (
                "interaction radius/radius",
                envelope.maximum_interaction_radius_to_radius,
            ),
            (
                "interaction radius/depth",
                envelope.maximum_interaction_radius_to_depth,
            ),
            (
                "interaction radius/layer",
                envelope.maximum_interaction_radius_to_layer,
            ),
            ("nominal strain", envelope.maximum_nominal_strain),
            (
                "repulsive pressure/yield",
                envelope.maximum_repulsive_pressure_to_yield,
            ),
            (
                "adhesive stress/yield",
                envelope.maximum_adhesive_stress_to_yield,
            ),
        ] {
            positive(value, field)?;
            if value >= 1.0 {
                return Err(DmtError::InvalidInput { field });
            }
        }
        for (field, value, cap) in [
            (
                "maximum range/radius",
                envelope.maximum_interaction_range_to_radius,
                0.01,
            ),
            (
                "maximum Tabor parameter",
                envelope.maximum_tabor_parameter,
                0.1,
            ),
        ] {
            positive(value, field)?;
            if value > cap {
                return Err(DmtError::InvalidInput { field });
            }
        }
        let elastic = normal.elastic();
        if elastic
            .surface_a_state()
            .resolved()
            .properties()
            .iter()
            .chain(elastic.surface_b_state().resolved().properties())
            .any(|p| !p.is_static_at_fixed_external_conditions())
            || !normal
                .interface_state()
                .property(ADHESION_ENERGY_PROPERTY)
                .is_some_and(|p| p.is_static_at_fixed_external_conditions())
        {
            return Err(DmtError::NonstaticSource);
        }
        let e = elastic.reduced_modulus_pa();
        let yield_stress = elastic.limiting_yield_stress_pa();
        let hertz_coefficient = (4.0 / 3.0) * e * effective_radius_m.sqrt();
        let pull_off_magnitude_n = (2.0 * PI) * (w * effective_radius_m);
        let tail_stiffness_n_per_m = pull_off_magnitude_n / interaction_range_m;
        let tabor_parameter =
            effective_radius_m.cbrt() * (w / e).cbrt().powi(2) / interaction_range_m;
        for (field, value) in [
            ("Hertz coefficient", hertz_coefficient),
            ("pull-off force", pull_off_magnitude_n),
            ("tail stiffness", tail_stiffness_n_per_m),
            ("Tabor parameter", tabor_parameter),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(DmtError::Unrepresentable { field });
            }
        }
        let a_max = effective_radius_m.sqrt() * envelope.maximum_indentation_m.sqrt();
        let footprint =
            a_max.hypot(effective_radius_m.sqrt() * interaction_range_m.sqrt() * 2.0_f64.sqrt());
        for (ratio, value, limit) in [
            (
                "Tabor parameter",
                tabor_parameter,
                envelope.maximum_tabor_parameter,
            ),
            (
                "interaction range/radius",
                interaction_range_m / effective_radius_m,
                envelope.maximum_interaction_range_to_radius,
            ),
            (
                "interaction radius/radius",
                footprint / effective_radius_m,
                envelope.maximum_interaction_radius_to_radius,
            ),
            (
                "interaction radius/depth",
                footprint / envelope.half_space_depth_m,
                envelope.maximum_interaction_radius_to_depth,
            ),
            (
                "interaction radius/layer",
                footprint / envelope.layer_thickness_m,
                envelope.maximum_interaction_radius_to_layer,
            ),
            (
                "nominal strain",
                a_max / effective_radius_m,
                envelope.maximum_nominal_strain,
            ),
            (
                "repulsive pressure/yield",
                (2.0 / PI) * (e / yield_stress) * (a_max / effective_radius_m),
                envelope.maximum_repulsive_pressure_to_yield,
            ),
            (
                "adhesive stress/yield",
                (w / interaction_range_m) / yield_stress,
                envelope.maximum_adhesive_stress_to_yield,
            ),
        ] {
            if !value.is_finite() || value > limit {
                return Err(DmtError::OutsideApplicability {
                    ratio,
                    value,
                    limit,
                });
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.dmt-dugdale-sphere.v1");
        hash.update(normal.identity().as_bytes());
        for value in [
            effective_radius_m,
            interaction_range_m,
            envelope.maximum_indentation_m,
            envelope.half_space_depth_m,
            envelope.layer_thickness_m,
            envelope.maximum_interaction_radius_to_radius,
            envelope.maximum_interaction_radius_to_depth,
            envelope.maximum_interaction_radius_to_layer,
            envelope.maximum_nominal_strain,
            envelope.maximum_repulsive_pressure_to_yield,
            envelope.maximum_adhesive_stress_to_yield,
            envelope.maximum_interaction_range_to_radius,
            envelope.maximum_tabor_parameter,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        let model = Self {
            normal,
            radius_m: effective_radius_m,
            interaction_range_m,
            envelope,
            identity: hash.finalize(),
            hertz_coefficient,
            pull_off_magnitude_n,
            tail_stiffness_n_per_m,
            tabor_parameter,
        };
        // Include the deepest interaction energy and the largest contact response.
        model.evaluate(0.0)?;
        model.evaluate(envelope.maximum_indentation_m)?;
        Ok(model)
    }

    /// Original material and ordered-interface state and complete receipts.
    #[must_use]
    pub const fn normal_state(&self) -> &BoundNormalInterfaceState {
        &self.normal
    }
    /// Model identity includes the explicit tail and every applicability limit.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Effective sphere/plane radius, m.
    #[must_use]
    pub const fn effective_radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Explicit cohesive interaction range, m.
    #[must_use]
    pub const fn interaction_range_m(&self) -> f64 {
        self.interaction_range_m
    }
    /// Selected operating envelope.
    #[must_use]
    pub const fn envelope(&self) -> DmtEnvelope {
        self.envelope
    }
    /// Dimensionless Tabor parameter under the envelope's stated convention.
    #[must_use]
    pub const fn tabor_parameter(&self) -> f64 {
        self.tabor_parameter
    }
    /// Force-controlled contact pull-off, -2 pi R W, N.
    #[must_use]
    pub fn force_controlled_pull_off_n(&self) -> f64 {
        -self.pull_off_magnitude_n
    }

    /// Evaluate signed approach. Force and potential are continuous; force
    /// derivatives jump at -z0 and zero. Negative tail stiffness must not be
    /// replaced with zero: global stability depends on the surrounding system.
    pub fn evaluate(&self, indentation_m: f64) -> Result<DmtResponse, DmtError> {
        if !indentation_m.is_finite() {
            return Err(DmtError::InvalidInput {
                field: "signed indentation",
            });
        }
        if indentation_m > self.envelope.maximum_indentation_m {
            return Err(DmtError::OutsideApplicability {
                ratio: "indentation",
                value: indentation_m,
                limit: self.envelope.maximum_indentation_m,
            });
        }
        let z = self.interaction_range_m;
        let f0 = self.pull_off_magnitude_n;
        let negative_tail = -self.tail_stiffness_n_per_m;
        let (regime, a, repulsive, force, energy, slopes) = if indentation_m <= -z {
            (
                DmtRegime::Separated,
                0.0,
                0.0,
                0.0,
                0.0,
                if indentation_m == -z {
                    [0.0, negative_tail]
                } else {
                    [0.0; 2]
                },
            )
        } else if indentation_m < 0.0 {
            let overlap = indentation_m + z;
            let fraction = overlap / z;
            let force = -f0 * fraction;
            (
                DmtRegime::AttractiveGap,
                0.0,
                0.0,
                force,
                0.5 * force * overlap,
                [negative_tail; 2],
            )
        } else {
            let root_delta = indentation_m.sqrt();
            let a = self.radius_m.sqrt() * root_delta;
            let repulsive = self.hertz_coefficient * indentation_m * root_delta;
            let force = repulsive - f0;
            let energy = 0.4 * repulsive * indentation_m - f0 * (indentation_m + 0.5 * z);
            let tangent = 1.5 * self.hertz_coefficient * root_delta;
            (
                DmtRegime::Contact,
                a,
                repulsive,
                force,
                energy,
                if indentation_m == 0.0 {
                    [negative_tail, 0.0]
                } else {
                    [tangent; 2]
                },
            )
        };
        let area = PI * a * a;
        for (field, value) in [
            ("contact radius", a),
            ("contact area", area),
            ("repulsive force", repulsive),
            ("normal force", force),
            ("reversible energy", energy),
            ("left tangent", slopes[0]),
            ("right tangent", slopes[1]),
        ] {
            if !value.is_finite() {
                return Err(DmtError::Unrepresentable { field });
            }
        }
        if indentation_m > 0.0 && (a <= 0.0 || area <= 0.0 || repulsive <= 0.0) {
            return Err(DmtError::Unrepresentable {
                field: "positive contact response",
            });
        }
        let authority = match self
            .normal
            .elastic()
            .interface()
            .interface()
            .provenance()
            .authority()
        {
            InputAuthority::SyntheticFixture => InputAuthority::SyntheticFixture,
            _ => InputAuthority::Estimated,
        };
        Ok(DmtResponse {
            normal_state_identity: self.normal.identity(),
            regime,
            indentation_m,
            contact_radius_m: a,
            contact_area_m2: area,
            repulsive_normal_force_n: repulsive,
            normal_force_n: force,
            tangent_n_per_m: if slopes[0] == slopes[1] {
                Some(slopes[0])
            } else {
                None
            },
            one_sided_tangents_n_per_m: slopes,
            reversible_energy_j: energy,
            authority,
        })
    }

    /// Solve the stable contacting branch for an applied signed normal load.
    /// The finite-gap negative-stiffness branch is deliberately not a second
    /// equilibrium under pure force control. A coupled machine must solve it.
    pub fn evaluate_contact_load(&self, load_n: f64) -> Result<DmtResponse, DmtError> {
        if !load_n.is_finite() {
            return Err(DmtError::InvalidInput {
                field: "applied normal load",
            });
        }
        let pull_off = self.force_controlled_pull_off_n();
        if load_n < pull_off {
            return Err(DmtError::NoContactEquilibrium {
                load_n,
                pull_off_n: pull_off,
            });
        }
        let repulsive = load_n + self.pull_off_magnitude_n;
        let delta = (repulsive / self.hertz_coefficient).cbrt().powi(2);
        if !delta.is_finite() || (repulsive > 0.0 && delta == 0.0) {
            return Err(DmtError::Unrepresentable {
                field: "load-controlled indentation",
            });
        }
        self.evaluate(delta)
    }
}

fn positive(value: f64, field: &'static str) -> Result<(), DmtError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(DmtError::InvalidInput { field })
    }
}
