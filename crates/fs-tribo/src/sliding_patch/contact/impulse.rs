//! Maximum-dissipation impact impulses on a finite Hertz contact patch.
//!
//! A prescribed compressive normal impulse supplies the total local Coulomb
//! budget. Tangential impulse and normal-axis angular impulse share the exact
//! integrated Hertz surface. Independent point-friction and torsional clamps
//! would overspend that budget during mixed sliding and spinning.
//!
//! This conditional tangent/spin problem does not determine normal restitution
//! or the normal impulse. It assumes the supplied Hertz pressure-impulse shape
//! and radius stay fixed through the impact. No elapsed time, slip history or
//! thermal partition is inferred from this instantaneous mechanical loss.

use super::{
    CircularWrenchLaw, HertzWrenchConfig, HertzWrenchError, HertzWrenchMode, HertzWrenchPhase,
    Matrix, Problem, SPIN_RADIUS, Vector, dot, factor_circular_response, finite_vector, max_abs,
    multiply, project_circular, shape_factors, sub,
};
use crate::TriboError;

/// Finite-patch geometry and one explicitly declared impact-friction set.
/// Static/kinetic phase evolution requires a finite-time contact owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpulsePatch {
    radius_m: f64,
    friction_coefficient: f64,
}

/// Prescribed normal impulse and the complete physical tangent/spin response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpulseInput {
    pub normal_impulse_ns: f64,
    pub velocity_before_m_s: [f64; 2],
    pub spin_before_rad_s: f64,
    /// Maps `[p_x,p_y,L/a]` to `[dv_x,dv_y,a*domega]` [1/kg].
    /// Supply `S*W_physical*S`, `S=diag(1,1,a)`, retaining cross coupling.
    /// The complete matrix must be symmetric and positive definite.
    pub scaled_mobility_per_kg: Matrix,
}

/// Bounded numerical work and independently named physical impulse tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpulseControls {
    pub maximum_static_iterations: usize,
    pub maximum_iterations: usize,
    pub maximum_evaluations: usize,
    pub impulse_tolerance_ns: f64,
    pub spin_impulse_tolerance_nms: f64,
    /// All radius-scaled velocity compatibility coordinates [m/s].
    pub velocity_tolerance_m_s: f64,
    pub energy_tolerance_j: f64,
    pub relative_matrix_tolerance: f64,
    pub maximum_impulse_ns: f64,
    pub maximum_spin_impulse_nms: f64,
}

impl HertzImpulseControls {
    fn projection(self) -> HertzWrenchConfig {
        // The numerical projection pairs conjugate coordinates. These are
        // impulse/velocity scalars, never a call to finite-time contact advance.
        HertzWrenchConfig {
            maximum_static_iterations: self.maximum_static_iterations,
            maximum_iterations: self.maximum_iterations,
            maximum_evaluations: self.maximum_evaluations,
            force_tolerance_n: self.impulse_tolerance_ns,
            moment_tolerance_nm: self.spin_impulse_tolerance_nms,
            displacement_tolerance_m: self.velocity_tolerance_m_s,
            energy_tolerance_j: self.energy_tolerance_j,
            relative_matrix_tolerance: self.relative_matrix_tolerance,
            maximum_force_n: self.maximum_impulse_ns,
            maximum_spin_moment_nm: self.maximum_spin_impulse_nms,
        }
    }

    pub fn validate(self) -> Result<(), HertzImpulseError> {
        self.projection().validate().map_err(Into::into)
    }
}

/// Candidate only; the mechanical owner applies and audits its body impulses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzImpulseReport {
    /// Acting impulse on the side whose relative velocity was supplied [N s].
    pub impulse_ns: [f64; 2],
    /// Acting normal-axis angular impulse about the patch centre [N m s].
    pub spin_impulse_nms: f64,
    /// Actual pre-impact velocity plus the original mobility response.
    /// Small floating arrest residuals are reported, never projected away.
    pub velocity_after_m_s: [f64; 2],
    pub spin_after_rad_s: f64,
    pub mode: HertzWrenchMode,
    /// Nonnegative loss from physical midpoint impulse/velocity work [J].
    pub dissipated_j: f64,
    /// Independently reconstructed in the original inverse-mobility metric [J].
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub velocity_residual_m_s: Vector,
    /// Native Hertz support stationarity at the admitted endpoint [N s].
    pub impulse_residual_ns: f64,
    pub spin_impulse_residual_nms: f64,
    /// Endpoint support optimality gap, distinct from impact loss [J].
    pub duality_gap_j: f64,
    pub iterations: usize,
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HertzImpulseError {
    Invalid(&'static str),
    Constitutive(TriboError),
    Cancelled,
    ImpulseBudget {
        impulse_ns: f64,
        spin_impulse_nms: f64,
    },
    BoundaryUnresolved {
        impulse_bracket_ns: f64,
        spin_impulse_bracket_nms: f64,
        evaluations: usize,
    },
    Unresolved {
        iterations: usize,
        evaluations: usize,
        impulse_residual_ns: f64,
        spin_impulse_residual_nms: f64,
        velocity_residual_m_s: f64,
        energy_residual_j: f64,
    },
    EnergyBalance {
        dissipated_j: f64,
        residual_j: f64,
    },
}

impl From<HertzWrenchError> for HertzImpulseError {
    fn from(error: HertzWrenchError) -> Self {
        match error {
            HertzWrenchError::Invalid(field) => Self::Invalid(field),
            HertzWrenchError::Law(error) => Self::Constitutive(error),
            HertzWrenchError::Cancelled => Self::Cancelled,
            HertzWrenchError::EffortBudget {
                force_n,
                spin_moment_nm,
            } => Self::ImpulseBudget {
                impulse_ns: force_n,
                spin_impulse_nms: spin_moment_nm,
            },
            HertzWrenchError::StaticBoundaryUnresolved {
                force_bracket_n,
                moment_bracket_nm,
                evaluations,
            } => Self::BoundaryUnresolved {
                impulse_bracket_ns: force_bracket_n,
                spin_impulse_bracket_nms: moment_bracket_nm,
                evaluations,
            },
            HertzWrenchError::Unresolved {
                iterations,
                evaluations,
                force_residual_n,
                spin_moment_residual_nm,
                displacement_residual_m,
                energy_residual_j,
            } => Self::Unresolved {
                iterations,
                evaluations,
                impulse_residual_ns: force_residual_n,
                spin_impulse_residual_nms: spin_moment_residual_nm,
                velocity_residual_m_s: displacement_residual_m,
                energy_residual_j,
            },
        }
    }
}

impl core::fmt::Display for HertzImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz contact impulse: {self:?}")
    }
}
impl std::error::Error for HertzImpulseError {}

impl HertzImpulsePatch {
    pub fn new(radius_m: f64, friction_coefficient: f64) -> Result<Self, HertzImpulseError> {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(HertzImpulseError::Invalid("positive_finite_patch_radius"));
        }
        if !friction_coefficient.is_finite() || friction_coefficient < 0.0 {
            return Err(HertzImpulseError::Invalid(
                "nonnegative_finite_friction_coefficient",
            ));
        }
        Ok(Self {
            radius_m,
            friction_coefficient,
        })
    }

    #[must_use]
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }

    #[must_use]
    pub const fn friction_coefficient(&self) -> f64 {
        self.friction_coefficient
    }

    /// Minimize `p dot v_before + p dot W*p/2` over the original integrated
    /// Hertz impulse set, retaining the complete translation/spin response.
    /// The bounded, allocation-free projector shares its exact analytic support
    /// with finite contact. Complete arrest is checked against that joint set.
    pub fn solve(
        &self,
        input: &HertzImpulseInput,
        controls: HertzImpulseControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzImpulseReport, HertzImpulseError> {
        if cancelled() {
            return Err(HertzImpulseError::Cancelled);
        }
        controls.validate()?;
        if !input.normal_impulse_ns.is_finite() || input.normal_impulse_ns < 0.0 {
            return Err(HertzImpulseError::Invalid(
                "nonnegative_finite_normal_impulse",
            ));
        }
        if !input.spin_before_rad_s.is_finite() {
            return Err(HertzImpulseError::Invalid("finite_spin_velocity"));
        }
        let before = [
            input.velocity_before_m_s[0],
            input.velocity_before_m_s[1],
            self.radius_m * input.spin_before_rad_s,
        ];
        finite_vector(before, "finite_contact_velocity")?;
        if input.spin_before_rad_s != 0.0 && before[2] == 0.0 {
            return Err(HertzImpulseError::Invalid("scaled_spin_underflow"));
        }
        let capacity = self.friction_coefficient * input.normal_impulse_ns;
        if !capacity.is_finite()
            || (self.friction_coefficient != 0.0
                && input.normal_impulse_ns != 0.0
                && capacity == 0.0)
        {
            return Err(HertzImpulseError::Invalid("representable_impulse_capacity"));
        }
        let config = controls.projection();
        let factor = factor_circular_response(input.scaled_mobility_per_kg, config)?;
        let law = ImpulseSupport {
            radius_m: self.radius_m,
            capacity,
        };
        let projected = project_circular(
            &law,
            before,
            input.scaled_mobility_per_kg,
            config,
            &mut cancelled,
        )?;
        let impulse_ns = projected.force_n;
        let spin_impulse_nms = projected.spin_moment_nm;
        let impulse = [
            impulse_ns[0],
            impulse_ns[1],
            spin_impulse_nms / self.radius_m,
        ];
        finite_vector(impulse, "representable_impulse")?;
        let response = multiply(input.scaled_mobility_per_kg, impulse);
        let after = core::array::from_fn(|i| before[i] + response[i]);
        finite_vector(after, "finite_endpoint_velocity")?;
        let target = [
            projected.slip_displacement_m[0],
            projected.slip_displacement_m[1],
            self.radius_m * projected.spin_angle_rad,
        ];
        let velocity_residual = sub(after, target);
        let spin_after = after[2] / self.radius_m;
        if !spin_after.is_finite() || (after[2] != 0.0 && spin_after == 0.0) {
            return Err(HertzImpulseError::Invalid("representable_spin_endpoint"));
        }
        if max_abs(velocity_residual) > controls.velocity_tolerance_m_s {
            return Err(HertzImpulseError::Unresolved {
                iterations: projected.iterations,
                evaluations: projected.evaluations,
                impulse_residual_ns: projected.force_residual_n,
                spin_impulse_residual_nms: projected.spin_moment_residual_nm,
                velocity_residual_m_s: max_abs(velocity_residual),
                energy_residual_j: 0.0,
            });
        }
        // Audit the support at the velocity produced by the *returned* physical
        // impulse. Radius conversion and the original (unaveraged) mobility can
        // differ from the internal candidate by floating-point roundoff.
        let mut impulse_residual_ns = projected.force_residual_n;
        let mut spin_impulse_residual_nms = projected.spin_moment_residual_nm;
        let mut duality_gap_j = projected.duality_gap_j;
        let mut evaluations = projected.evaluations;
        if projected.mode == HertzWrenchMode::Sliding {
            let mut audit = Problem {
                law: &law,
                response: input.scaled_mobility_per_kg,
                free: before,
                config,
                evaluations,
                cancelled: &mut cancelled,
            };
            let native = audit.support(after)?;
            evaluations = audit.evaluations;
            let residual = core::array::from_fn(|i| impulse[i] + native[i]);
            impulse_residual_ns = fs_math::det::hypot(residual[0], residual[1]);
            spin_impulse_residual_nms = self.radius_m * residual[2].abs();
            duality_gap_j = dot(residual, after);
            if !impulse_residual_ns.is_finite()
                || !spin_impulse_residual_nms.is_finite()
                || !duality_gap_j.is_finite()
                || impulse_residual_ns > controls.impulse_tolerance_ns
                || spin_impulse_residual_nms > controls.spin_impulse_tolerance_nms
                || duality_gap_j.abs() > controls.energy_tolerance_j
            {
                return Err(HertzImpulseError::Unresolved {
                    iterations: projected.iterations,
                    evaluations,
                    impulse_residual_ns,
                    spin_impulse_residual_nms,
                    velocity_residual_m_s: max_abs(velocity_residual),
                    energy_residual_j: duality_gap_j.abs(),
                });
            }
        }
        // Physical translation and angular ports are evaluated independently
        // of the radius-scaled convex support gap. Arrest can lose energy even
        // though its endpoint support pairing is exactly zero.
        let midpoint = core::array::from_fn(|i| 0.5 * before[i] + 0.5 * after[i]);
        let spin_midpoint = 0.5 * input.spin_before_rad_s + 0.5 * spin_after;
        let dissipated_j = -impulse_ns[0] * midpoint[0]
            - impulse_ns[1] * midpoint[1]
            - spin_impulse_nms * spin_midpoint;
        // Use measured velocity change and the original inverse mobility,
        // avoiding cancellation from subtracting two large kinetic energies.
        let momentum_change =
            factor
                .solve(sub(after, before))
                .ok_or(HertzImpulseError::Invalid(
                    "inverse_mobility_energy_overflow",
                ))?;
        let kinetic_energy_change_j = dot(midpoint, momentum_change);
        let energy_residual_j = kinetic_energy_change_j + dissipated_j;
        if !dissipated_j.is_finite()
            || !kinetic_energy_change_j.is_finite()
            || !energy_residual_j.is_finite()
            || dissipated_j < 0.0
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(HertzImpulseError::EnergyBalance {
                dissipated_j,
                residual_j: energy_residual_j,
            });
        }
        if cancelled() {
            return Err(HertzImpulseError::Cancelled);
        }
        Ok(HertzImpulseReport {
            impulse_ns,
            spin_impulse_nms,
            velocity_after_m_s: [after[0], after[1]],
            spin_after_rad_s: spin_after,
            mode: projected.mode,
            dissipated_j,
            kinetic_energy_change_j,
            energy_residual_j,
            velocity_residual_m_s: velocity_residual,
            impulse_residual_ns,
            spin_impulse_residual_nms,
            duality_gap_j,
            iterations: projected.iterations,
            evaluations,
        })
    }
}

/// The Hertz normalized support is homogeneous. Its capacity here is a real
/// impulse, not a normal force produced from an invented contact duration.
#[derive(Clone, Copy)]
struct ImpulseSupport {
    radius_m: f64,
    capacity: f64,
}

impl CircularWrenchLaw for ImpulseSupport {
    fn radius(&self) -> f64 {
        self.radius_m
    }
    fn static_capacity(&self) -> f64 {
        self.capacity
    }
    fn kinetic_capacity(&self) -> f64 {
        self.capacity
    }
    fn spin_radius_fraction(&self) -> f64 {
        SPIN_RADIUS
    }
    fn shape(&self, speed: f64, edge_speed: f64) -> Result<(f64, f64), TriboError> {
        Ok(shape_factors(speed, edge_speed))
    }
    fn phase_law(&self, _: HertzWrenchPhase) -> Self {
        *self
    }
    fn static_projection(&self) -> Self {
        *self
    }
}
