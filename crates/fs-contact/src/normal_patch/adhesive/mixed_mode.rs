//! Circular, ideally brittle full-stick JKR contact under normal and shear load.
//!
//! The SAME source supplies E*, W, R and both elastic bodies. Mindlin's shared
//! full-stick compliance gives k_t = 8 G* a. The conservative potential is
//! U = E* a [(delta-a^2/(3R))^2 + 4a^4/(45R^2)]
//!     + k_t |u|^2/2 - pi W a^2.
//! Radius stationarity gives E*(delta-a^2/R)^2 + 4G*|u|^2 = 2pi W a.
//! Thus shear changes the actual area, both forces and the coupled stiffness;
//! adding an independent shear force to an unchanged JKR area is insufficient.
//!
//! This explicitly selects the reversible, mode-independent Griffith model of
//! Savkoor & Briggs (1977), DOI 10.1098/rspa.1977.0123, with circular/perimeter-
//! averaged mode II/III release and decoupled isotropic bulk compliance. It is
//! not a fitted mixed-mode toughness, partial-slip or friction law. This brittle
//! model can overpredict measured area reduction: see Ciavarella & Papangelo
//! (2020), DOI 10.1007/s11249-020-01318-5, Eq. 12 with lambda=1. No empirical
//! loss fraction is introduced. The full-stick branch ends at its control-
//! dependent fold; subsequent sliding, nucleation or dynamic separation needs
//! the caller's physical owner, not an invented zero-area equilibrium or heat.

/// Finite tangential apparatus stiffness and competing sourced slip/peel limits.
pub mod loading;

/// Mixed-mode withdrawal through finite normal and tangential apparatus springs.
pub mod pull_off;

use super::{JkrError, SourceJkrSphereContact};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::{InputAuthority, TriboError, adhesive::AdhesiveContactState};

/// Additional nominal shear screens; original JKR source limits still apply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeEnvelope {
    /// Maximum |u|/a, strictly between 0 and 1.
    pub maximum_tangential_displacement_to_radius: f64,
    /// Maximum mean |T|/(pi a^2) divided by the original limiting yield stress.
    /// This does not bound the singular traction at the contact edge.
    pub maximum_mean_shear_stress_to_yield: f64,
}

/// Which complete set of conjugate coordinates is held during perturbations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrMixedModeControl {
    /// Normal approach and the two tangential displacements are prescribed.
    Displacement,
    /// Normal force and the two tangential forces are prescribed.
    Force,
}

/// No source/geometry substitution or post-instability slip rule is inferred.
#[derive(Clone, Debug, PartialEq)]
pub enum JkrMixedModeError {
    Normal(JkrError),
    Tangential(TriboError),
    Invalid {
        field: &'static str,
    },
    BelowNormalPullOff {
        force_n: f64,
        minimum_force_n: f64,
    },
    /// Limit is metres for displacement control and newtons for force control.
    NoFullStickEquilibrium {
        control: JkrMixedModeControl,
        magnitude: f64,
        limit: f64,
    },
    OutsideApplicability {
        ratio: &'static str,
        value: f64,
        limit: f64,
    },
    NoConvergence,
}
impl From<JkrError> for JkrMixedModeError {
    fn from(e: JkrError) -> Self {
        Self::Normal(e)
    }
}
impl From<TriboError> for JkrMixedModeError {
    fn from(e: TriboError) -> Self {
        Self::Tangential(e)
    }
}
impl core::fmt::Display for JkrMixedModeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "full-stick mixed-mode JKR refusal: {self:?}")
    }
}
impl std::error::Error for JkrMixedModeError {}

/// One attached local minimum, including its limiting marginal point.
/// All efforts enter the contact potential; body reactions have opposite sign.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeResponse {
    normal_identity: ContentHash,
    model_identity: ContentHash,
    pub control: JkrMixedModeControl,
    /// True at loss of radius stability under the declared control.
    pub marginal: bool,
    pub indentation_m: f64,
    pub tangential_displacement_m: [f64; 2],
    pub contact_radius_m: f64,
    pub contact_area_m2: f64,
    /// Original Hertz repulsive component, not the signed resultant.
    pub repulsive_normal_force_n: f64,
    pub normal_force_n: f64,
    pub tangential_force_n: [f64; 2],
    /// Positive square-form normal elastic storage before surface energy.
    pub normal_elastic_energy_j: f64,
    pub tangential_elastic_energy_j: f64,
    /// Negative adhesion potential -W times the actual contact area.
    pub surface_energy_j: f64,
    pub reversible_energy_j: f64,
    pub normal_energy_release_rate_j_m2: f64,
    pub shear_energy_release_rate_j_m2: f64,
    /// G_I + G_II/III - W, without hiding the floating-point equilibrium error.
    pub griffith_residual_j_m2: f64,
    /// Radius second derivative of U at fixed (delta,u1,u2).
    pub displacement_radius_curvature_n_per_m: f64,
    /// Radius second derivative after eliminating all displacement coordinates
    /// at fixed (N,T1,T2). It can be negative on a displacement-stable branch.
    pub force_radius_curvature_n_per_m: f64,
    /// Symmetric d(N,T1,T2)/d(delta,u1,u2), with radius equilibrated. Absent
    /// at the displacement fold, where this condensation is singular.
    pub coupled_stiffness_n_per_m: Option<[[f64; 3]; 3]>,
    /// Symmetric d(delta,u1,u2)/d(N,T1,T2). Absent when the force-controlled
    /// branch is unstable or marginal, even if displacement control is stable.
    pub coupled_compliance_m_per_n: Option<[[f64; 3]; 3]>,
    pub authority: InputAuthority,
}
impl JkrMixedModeResponse {
    #[must_use]
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_identity
    }
    #[must_use]
    pub const fn model_identity(&self) -> ContentHash {
        self.model_identity
    }
    /// Actual reduced area/load for source-bound junction-capacity exchange.
    /// A consumer must solve shear loading consistently: combining this elastic
    /// storage with an independently imposed sliding traction is not automatic.
    pub fn friction_contact_state(&self) -> Result<AdhesiveContactState, TriboError> {
        AdhesiveContactState::new(self.contact_area_m2, self.repulsive_normal_force_n)
    }
}

/// Explicit conservative model with the complete unchanged normal source.
/// Runtime is allocation-free: force control is closed form; displacement
/// control uses two monotone scalar solves of at most 64 iterations each.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrMixedMode {
    normal: SourceJkrSphereContact,
    envelope: JkrMixedModeEnvelope,
    identity: ContentHash,
    // Shared Mindlin stiffness per contact radius, equal to 8 G* [Pa].
    tangential_stiffness_per_radius: f64,
    // E*/(4G*) weights squared shear force in the Griffith equation.
    shear_force_weight: f64,
    tangential_displacement_scale_m: f64,
}
impl SourceJkrMixedMode {
    /// Select circular full stick and constant, mode-independent W explicitly.
    /// No Coulomb coefficient or shear strength is fabricated. Source E and nu
    /// are exactly those already admitted by the original normal owner.
    pub fn new(
        normal: SourceJkrSphereContact,
        envelope: JkrMixedModeEnvelope,
    ) -> Result<Self, JkrMixedModeError> {
        for (value, field) in [
            (
                envelope.maximum_tangential_displacement_to_radius,
                "maximum shear displacement/radius",
            ),
            (
                envelope.maximum_mean_shear_stress_to_yield,
                "maximum mean shear/yield",
            ),
        ] {
            if !(value.is_finite() && value > 0.0 && value < 1.0) {
                return Err(JkrMixedModeError::Invalid { field });
            }
        }
        let elastic = normal.normal_state().elastic();
        let solids = [elastic.surface_a_state(), elastic.surface_b_state()];
        let poisson = solids.map(|s| s.poisson_ratio());
        let shear = solids.map(|s| (0.5 * s.young_modulus_pa()) / (1.0 + s.poisson_ratio()));
        let radius = normal.spinodal_radius_m;
        let k = fs_tribo::mindlin::full_stick_stiffness_n_per_m(radius, shear, poisson)? / radius;
        let weight = 2.0 * elastic.reduced_modulus_pa() / k;
        let shift_scale = normal.displacement_scale_m * weight.sqrt();
        for (value, field) in [
            (k, "tangential stiffness/radius"),
            (weight, "shear force weight"),
            (shift_scale, "tangential displacement scale"),
        ] {
            positive(value, field)?;
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-full-stick-griffith.v1");
        hash.update(normal.identity().as_bytes());
        for value in [
            envelope.maximum_tangential_displacement_to_radius,
            envelope.maximum_mean_shear_stress_to_yield,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            normal,
            envelope,
            identity: hash.finalize(),
            tangential_stiffness_per_radius: k,
            shear_force_weight: weight,
            tangential_displacement_scale_m: shift_scale,
        })
    }
    #[must_use]
    pub const fn normal_source(&self) -> &SourceJkrSphereContact {
        &self.normal
    }
    #[must_use]
    pub const fn envelope(&self) -> JkrMixedModeEnvelope {
        self.envelope
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Shared elastic G*; no friction or mode-mixity parameter is fitted here.
    #[must_use]
    pub fn reduced_shear_modulus_pa(&self) -> f64 {
        self.tangential_stiffness_per_radius / 8.0
    }

    /// Magnitude |u| at the fixed-approach peeling fold. Physical source/shear
    /// limits can end the admitted model earlier; this is not a slip threshold.
    pub fn critical_tangential_displacement_m(
        &self,
        indentation_m: f64,
    ) -> Result<f64, JkrMixedModeError> {
        let (_, critical_square) = self.displacement_fold(indentation_m)?;
        finite(
            self.tangential_displacement_scale_m * critical_square.sqrt(),
            "critical tangential displacement",
        )
    }
    /// Magnitude |T| at the simultaneous force-controlled fold.
    /// N must first lie above the original pure-normal JKR pull-off.
    pub fn critical_tangential_force_n(
        &self,
        normal_force_n: f64,
    ) -> Result<f64, JkrMixedModeError> {
        let discriminant = self.force_discriminant(normal_force_n)?;
        finite(
            (18.0 * self.normal.force_scale_n / self.shear_force_weight.sqrt())
                * discriminant.sqrt(),
            "critical tangential force",
        )
    }

    /// Equilibrate actual radius at fixed signed approach and vector shift.
    /// The larger local minimum is followed to its radius fold, including
    /// displacement-stable states that force control cannot retain.
    pub fn evaluate_displacement_controlled(
        &self,
        indentation_m: f64,
        tangential_displacement_m: [f64; 2],
    ) -> Result<JkrMixedModeResponse, JkrMixedModeError> {
        let magnitude = norm(tangential_displacement_m)?;
        let (fold, critical_square) = self.displacement_fold(indentation_m)?;
        let limit = finite(
            self.tangential_displacement_scale_m * critical_square.sqrt(),
            "critical shear shift",
        )?;
        if magnitude > limit {
            return Err(JkrMixedModeError::NoFullStickEquilibrium {
                control: JkrMixedModeControl::Displacement,
                magnitude,
                limit,
            });
        }
        let (radius, offset) = if magnitude == 0.0 {
            let pure = self.normal.evaluate_attached(indentation_m)?;
            (
                pure.contact_radius_m,
                (pure.contact_radius_m / self.normal.spinodal_radius_m - fold).max(0.0),
            )
        } else {
            // Factor the distance to the fold, retaining adjacent input values.
            let ratio = magnitude / limit;
            let distance = critical_square * (1.0 - ratio) * (1.0 + ratio);
            let a = 4.0 * fold * fold + 8.0 / fold;
            let upper = (distance / a).sqrt().min(distance.sqrt().sqrt());
            let offset = positive_root(distance, upper, |t| {
                (
                    t * t * (t * t + 4.0 * fold * t + a),
                    2.0 * t * (2.0 * t * t + 6.0 * fold * t + a),
                )
            })?;
            (self.normal.spinodal_radius_m * (fold + offset), offset)
        };
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let curvature = (PI * w / 4.0)
            * offset
            * (4.0 * fold * fold + 8.0 / fold + 6.0 * fold * offset + 2.0 * offset * offset);
        self.response(
            radius,
            indentation_m,
            tangential_displacement_m,
            JkrMixedModeControl::Displacement,
            curvature,
            curvature - 4.0 * PI * w,
        )
    }

    /// Stable simultaneous applied normal/tangential force equilibrium. With
    /// H=4E*a^3/(3R), B=3piWR and alpha=E*/(4G*),
    /// H = N+B+sqrt(B^2+2BN-alpha|T|^2). The opposite root is unstable.
    pub fn evaluate_force_controlled(
        &self,
        normal_force_n: f64,
        tangential_force_n: [f64; 2],
    ) -> Result<JkrMixedModeResponse, JkrMixedModeError> {
        let magnitude = norm(tangential_force_n)?;
        let base = self.force_discriminant(normal_force_n)?;
        let limit = self.critical_tangential_force_n(normal_force_n)?;
        if magnitude > limit {
            return Err(JkrMixedModeError::NoFullStickEquilibrium {
                control: JkrMixedModeControl::Force,
                magnitude,
                limit,
            });
        }
        let ratio = if magnitude == 0.0 {
            0.0
        } else {
            magnitude / limit
        };
        let distance = (base * (1.0 - ratio) * (1.0 + ratio)).sqrt();
        let b = 18.0 * self.normal.force_scale_n;
        let radius =
            self.normal.spinodal_radius_m * (18.0 * (normal_force_n / b + 1.0 + distance)).cbrt();
        positive(radius, "force-controlled contact radius")?;
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let indentation =
            normal_force_n / (2.0 * e * radius) + radius * (radius / self.normal.radius_m) / 3.0;
        let stiffness = self.tangential_stiffness_per_radius * radius;
        let shift = tangential_force_n.map(|t| t / stiffness);
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let force_curvature = (2.0 * b / self.normal.radius_m) * distance;
        self.response(
            radius,
            indentation,
            shift,
            JkrMixedModeControl::Force,
            4.0 * PI * w + force_curvature,
            force_curvature,
        )
    }

    fn force_discriminant(&self, normal_force_n: f64) -> Result<f64, JkrMixedModeError> {
        finite(normal_force_n, "normal force")?;
        let minimum = self.normal.force_controlled_pull_off_n();
        if normal_force_n < minimum {
            return Err(JkrMixedModeError::BelowNormalPullOff {
                force_n: normal_force_n,
                minimum_force_n: minimum,
            });
        }
        // Dimensional subtraction preserves the distance to normal pull-off.
        finite(
            (normal_force_n - minimum) / (9.0 * self.normal.force_scale_n),
            "force discriminant",
        )
    }
    fn displacement_fold(&self, delta: f64) -> Result<(f64, f64), JkrMixedModeError> {
        self.normal.validate_indentation(delta)?;
        let minimum = self.normal.minimum_attached_indentation_m();
        if delta < minimum {
            return Err(JkrError::NoAttachedEquilibrium {
                indentation_m: delta,
                minimum_indentation_m: minimum,
            }
            .into());
        }
        let distance = finite(
            (delta - minimum) / self.normal.displacement_scale_m,
            "approach above spinodal",
        )?;
        // For y=a_fold/a_s=1+t, y^3-(delta/h)y-4=0 becomes
        // t*(t+2+4/(1+t))=(delta-delta_s)/h, strictly increasing for t>=0.
        let offset = positive_root(distance, (distance / 5.0).min(distance.sqrt()), |t| {
            let y = 1.0 + t;
            (t * (t + 2.0 + 4.0 / y), 2.0 * y + 4.0 / (y * y))
        })?;
        let fold = 1.0 + offset;
        let critical_square = 16.0 * offset * (fold * fold + fold + 1.0) / (fold * fold);
        Ok((fold, critical_square))
    }

    fn response(
        &self,
        radius: f64,
        delta: f64,
        shift: [f64; 2],
        control: JkrMixedModeControl,
        curvature: f64,
        force_curvature: f64,
    ) -> Result<JkrMixedModeResponse, JkrMixedModeError> {
        self.normal.validate_indentation(delta)?;
        self.normal.check_applicability(radius)?;
        let elastic = self.normal.normal_state().elastic();
        let e = elastic.reduced_modulus_pa();
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let k = self.tangential_stiffness_per_radius;
        let hertz_depth = radius * (radius / self.normal.radius_m);
        let mean_normal_displacement = delta - hertz_depth / 3.0;
        let normal_force_n = 2.0 * e * radius * mean_normal_displacement;
        let stiffness = k * radius;
        let tangential_force_n = shift.map(|u| stiffness * u);
        let area = PI * radius * radius;
        let shear_norm = norm(shift)?;
        for (ratio, value, limit) in [
            (
                "tangential displacement/radius",
                shear_norm / radius,
                self.envelope.maximum_tangential_displacement_to_radius,
            ),
            (
                "mean shear/yield",
                norm(tangential_force_n)? / area / elastic.limiting_yield_stress_pa(),
                self.envelope.maximum_mean_shear_stress_to_yield,
            ),
        ] {
            if !value.is_finite() || value > limit {
                return Err(JkrMixedModeError::OutsideApplicability {
                    ratio,
                    value,
                    limit,
                });
            }
        }
        let normal_energy = (e * radius * mean_normal_displacement) * mean_normal_displacement
            + ((4.0 / 45.0) * e * radius * hertz_depth) * hertz_depth;
        let shear_energy = (0.5 * stiffness * shear_norm) * shear_norm;
        let surface_energy = -w * area;
        let opening = hertz_depth - delta;
        let normal_release = (e * opening) * (opening / radius) / (2.0 * PI);
        let shear_release = shear_energy / (2.0 * area);
        let griffith_residual = normal_release + shear_release - w;
        let diagonal = [2.0 * e * radius, stiffness, stiffness];
        let cross = [-2.0 * e * opening, k * shift[0], k * shift[1]];
        let coupled_stiffness = if curvature > 0.0 {
            let mut matrix = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in i..3 {
                    let value = (if i == j { diagonal[i] } else { 0.0 })
                        - (cross[i] / curvature) * cross[j];
                    matrix[i][j] = value;
                    matrix[j][i] = value;
                }
            }
            Some(matrix)
        } else {
            None
        };
        let coupled_compliance = if force_curvature > 0.0 {
            Some(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    (if i == j { 1.0 / diagonal[i] } else { 0.0 })
                        + (cross[i] / diagonal[i]) * (cross[j] / diagonal[j]) / force_curvature
                })
            }))
        } else {
            None
        };
        let repulsive = (4.0 / 3.0) * e * radius * hertz_depth;
        let energy = normal_energy + shear_energy + surface_energy;
        for value in [
            area,
            normal_force_n,
            repulsive,
            normal_energy,
            shear_energy,
            surface_energy,
            energy,
            normal_release,
            shear_release,
            griffith_residual,
            curvature,
            force_curvature,
        ] {
            finite(value, "mixed-mode response")?;
        }
        for matrix in [coupled_stiffness, coupled_compliance]
            .into_iter()
            .flatten()
        {
            for value in matrix.into_iter().flatten() {
                finite(value, "coupled tangent")?;
            }
        }
        if radius <= 0.0
            || area <= 0.0
            || repulsive <= 0.0
            || curvature < 0.0
            || opening < 0.0
            || griffith_residual.abs() > 512.0 * f64::EPSILON * w
        {
            return Err(JkrMixedModeError::NoConvergence);
        }
        Ok(JkrMixedModeResponse {
            normal_identity: self.normal.normal_state().identity(),
            model_identity: self.identity,
            control,
            marginal: if control == JkrMixedModeControl::Displacement {
                curvature == 0.0
            } else {
                force_curvature == 0.0
            },
            indentation_m: delta,
            tangential_displacement_m: shift,
            contact_radius_m: radius,
            contact_area_m2: area,
            repulsive_normal_force_n: repulsive,
            normal_force_n,
            tangential_force_n,
            normal_elastic_energy_j: normal_energy,
            tangential_elastic_energy_j: shear_energy,
            surface_energy_j: surface_energy,
            reversible_energy_j: energy,
            normal_energy_release_rate_j_m2: normal_release,
            shear_energy_release_rate_j_m2: shear_release,
            griffith_residual_j_m2: griffith_residual,
            displacement_radius_curvature_n_per_m: curvature,
            force_radius_curvature_n_per_m: force_curvature,
            coupled_stiffness_n_per_m: coupled_stiffness,
            coupled_compliance_m_per_n: coupled_compliance,
            authority: self.normal.authority(),
        })
    }
}

/// Monotone, nonnegative root in a known analytic bracket. Both callers use
/// factored offsets, so losing significant digits at a peeling fold is avoided.
fn positive_root(
    target: f64,
    upper: f64,
    f: impl Fn(f64) -> (f64, f64),
) -> Result<f64, JkrMixedModeError> {
    if target == 0.0 {
        return Ok(0.0);
    }
    positive(target, "root target")?;
    let mut hi = upper * (1.0 + 16.0 * f64::EPSILON);
    positive(hi, "root upper bound")?;
    if f(hi).0 < target {
        return Err(JkrMixedModeError::NoConvergence);
    }
    let mut lo = 0.0;
    let mut x = hi;
    for _ in 0..64 {
        let (value, derivative) = f(x);
        let residual = value - target;
        if residual.abs() <= 8.0 * f64::EPSILON * target {
            return Ok(x);
        }
        if residual > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
        let newton = x - residual / derivative;
        let next = if newton > lo && newton < hi {
            newton
        } else {
            f64::midpoint(lo, hi)
        };
        if next == x {
            break;
        }
        x = next;
    }
    Err(JkrMixedModeError::NoConvergence)
}
fn norm(vector: [f64; 2]) -> Result<f64, JkrMixedModeError> {
    finite(vector[0], "tangential component")?;
    finite(vector[1], "tangential component")?;
    finite(vector[0].hypot(vector[1]), "tangential magnitude")
}
fn finite(value: f64, field: &'static str) -> Result<f64, JkrMixedModeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrMixedModeError::Invalid { field })
    }
}
fn positive(value: f64, field: &'static str) -> Result<(), JkrMixedModeError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(JkrMixedModeError::Invalid { field })
    }
}
