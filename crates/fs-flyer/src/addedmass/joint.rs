//! One joint factorization for the constant-body and articulated inertia paths.

use super::{AeroGeneralizedLoads, MAX_CONTROL_COORDS, refuse};
use crate::Refusal;
use fs_mbd::articulated::generalized::FreeFloatingTerms;

pub(crate) fn native_index(i: usize) -> usize {
    if i < 3 {
        i + 3
    } else if i < 6 {
        i - 3
    } else {
        i
    }
}

pub(crate) fn admit_fluid(loads: &AeroGeneralizedLoads) -> Result<(), Refusal> {
    let nc = loads.q_control_nonaccel.len();
    if nc > MAX_CONTROL_COORDS
        || loads.m_added_rc.len() != nc
        || loads.m_added_cc.len() != nc
        || loads.m_added_cc.iter().any(|row| row.len() != nc)
        || loads.q_added_bias.len() != 6 + nc
    {
        return Err(refuse(
            "joint-load-shape",
            "inconsistent or oversized generalized fluid blocks".into(),
            "supply matching rigid, control, cross and bias blocks within the control cap",
        ));
    }
    if !loads
        .m_added_rr
        .iter()
        .flatten()
        .chain(loads.m_added_rc.iter().flatten())
        .chain(loads.m_added_cc.iter().flatten())
        .chain(&loads.q_rigid_nonaccel)
        .chain(&loads.q_control_nonaccel)
        .chain(&loads.q_added_bias)
        .all(|x| x.is_finite())
    {
        return Err(refuse(
            "joint-load-invalid",
            "nonfinite fluid mass, load or bias".into(),
            "use finite assembled mass and loads",
        ));
    }
    for i in 0..6 {
        for j in 0..i {
            admit_symmetric(loads.m_added_rr[i][j], loads.m_added_rr[j][i])?;
        }
    }
    for i in 0..nc {
        for j in 0..i {
            admit_symmetric(loads.m_added_cc[i][j], loads.m_added_cc[j][i])?;
        }
    }
    Ok(())
}

fn admit_symmetric(a: f64, b: f64) -> Result<(), Refusal> {
    if (a - b).abs() > 128.0 * f64::EPSILON * a.abs().max(b.abs()) {
        Err(refuse(
            "joint-mass-invalid",
            "nonsymmetric generalized mass".into(),
            "use the original symmetric dry and fluid mass assemblies",
        ))
    } else {
        Ok(())
    }
}

/// Solve the complete moving dry-body and fluid generalized inertia together.
///
/// `dry` is evaluated at the same state as `loads`. The native articulated
/// owner supplies all dry body/control cross blocks, moving COM effects and
/// Euler/configuration bias. This adapter only permutes its `[angular, linear,
/// joints]` coordinates to the Flyer's `[linear, angular, controls]` convention.
/// Applied generalized loads exclude BOTH dry and fluid inertial bias. The
/// dry bias is on the left of the native equation; fluid bias is on the right.
///
/// A diagonal congruence equilibrates the mixed translational/rotational units
/// before the original fixed-order Cholesky. Near-singular normalized pivots
/// refuse; no diagonal regularization or fictitious control inertia is added.
/// Returns body-coordinate velocity derivatives and joint accelerations.
///
/// # Errors
/// Mismatched/oversized blocks, nonfinite or nonsymmetric inputs, nonpositive
/// effective mass, or unresolved scaled pivots.
pub fn solve_articulated_joint(
    dry: &FreeFloatingTerms,
    loads: &AeroGeneralizedLoads,
) -> Result<(Vec<f64>, Vec<f64>), Refusal> {
    solve_articulated_with_diagonal(dry, loads, &[])
}

/// One admitted physical dry-plus-fluid mass, factored once at a configuration.
///
/// All vectors use Flyer `[linear_body, angular_body, joint]` ordering. The
/// same factor serves force solves, every contact-space response column and
/// the final accumulated impulse. Its immutable snapshot must be rebuilt when
/// configuration, dry mass properties, fluid density or strip geometry changes.
/// No applied force, inertial bias, restitution or contact law is stored here.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedArticulatedMass {
    n: usize,
    mass: Vec<f64>,
    scale: Vec<f64>,
    lower: Vec<f64>,
}

impl PreparedArticulatedMass {
    /// Prepare the complete physical inertia at one common dry/fluid state.
    /// The load and bias fields are admitted for shape/finite consistency but
    /// do not enter this mass operator. Impulses never integrate those forces.
    ///
    /// # Errors
    /// Mismatched, nonfinite, nonsymmetric or unresolved positive inertia.
    pub fn new(dry: &FreeFloatingTerms, fluid: &AeroGeneralizedLoads) -> Result<Self, Refusal> {
        Self::with_control_diagonal(dry, fluid, &[])
    }

    /// Effective midpoint operator only: the extra diagonal is eliminated
    /// linear impedance, not physical storage and never used for an impact.
    pub(crate) fn with_control_diagonal(
        dry: &FreeFloatingTerms,
        loads: &AeroGeneralizedLoads,
        control_diagonal: &[f64],
    ) -> Result<Self, Refusal> {
        admit_fluid(loads)?;
        let nc = loads.q_control_nonaccel.len();
        let n = 6 + nc;
        if dry.kinematics.dof_count() != nc
            || dry.mass.len() != n
            || dry.mass.iter().any(|row| row.len() != n)
            || dry.bias.len() != n
            || (!control_diagonal.is_empty() && control_diagonal.len() != nc)
            || !control_diagonal.iter().all(|x| x.is_finite() && *x >= 0.0)
            || !dry
                .mass
                .iter()
                .flatten()
                .chain(&dry.bias)
                .all(|x| x.is_finite())
        {
            return Err(refuse(
                "articulated-mass-invalid",
                "dry generalized terms or implicit diagonal do not match the active controls"
                    .into(),
                "prepare native dry terms and fluid blocks at the same admitted state",
            ));
        }
        let mut mass = vec![0.0; n * n];
        for i in 0..n {
            let ni = native_index(i);
            for j in 0..n {
                let nj = native_index(j);
                admit_symmetric(dry.mass[ni][nj], dry.mass[nj][ni])?;
                let fluid = match (i < 6, j < 6) {
                    (true, true) => loads.m_added_rr[i][j],
                    (true, false) => loads.m_added_rc[j - 6][i],
                    (false, true) => loads.m_added_rc[i - 6][j],
                    (false, false) => loads.m_added_cc[i - 6][j - 6],
                };
                mass[i * n + j] = dry.mass[ni][nj] + fluid;
            }
            if i >= 6 && !control_diagonal.is_empty() {
                mass[i * n + i] += control_diagonal[i - 6];
            }
        }
        let scale: Vec<_> = (0..n).map(|i| mass[i * n + i].sqrt()).collect();
        if !mass.iter().all(|x| x.is_finite()) || !scale.iter().all(|s| s.is_finite() && *s > 0.0) {
            return Err(refuse(
                "joint-mass-not-pd",
                "invalid effective mass diagonal or unrepresentable mass".into(),
                "use representable positive dry and fluid inertia",
            ));
        }
        let mut normalized = mass.clone();
        for i in 0..n {
            for j in 0..n {
                normalized[i * n + j] = (normalized[i * n + j] / scale[i]) / scale[j];
            }
        }
        let lower = factor_positive(normalized, n, 128.0 * f64::EPSILON)?;
        Ok(Self {
            n,
            mass,
            scale,
            lower,
        })
    }

    /// Number of body and joint coordinates.
    #[must_use]
    pub const fn coordinate_count(&self) -> usize {
        self.n
    }

    /// Number of active scalar joints, excluding the six body coordinates.
    #[must_use]
    pub const fn dof_count(&self) -> usize {
        self.n - 6
    }

    fn admit_vector(&self, vector: &[f64]) -> Result<(), Refusal> {
        if vector.len() != self.n || !vector.iter().all(|x| x.is_finite()) {
            Err(refuse(
                "joint-vector-invalid",
                "generalized vector has the wrong shape or a nonfinite coordinate".into(),
                "supply one finite entry per prepared body/joint coordinate",
            ))
        } else {
            Ok(())
        }
    }

    /// Apply `M^-1` to a force or impulse, reusing the immutable factor.
    /// Force input gives a body-coordinate acceleration; impulse gives a
    /// velocity increment. Neither operation adds inertial bias or changes state.
    ///
    /// # Errors
    /// Invalid vector shape, nonfinite input or unrepresentable response.
    pub fn solve(&self, right_hand_side: &[f64]) -> Result<Vec<f64>, Refusal> {
        self.admit_vector(right_hand_side)?;
        let normalized = right_hand_side
            .iter()
            .zip(&self.scale)
            .map(|(x, s)| x / s)
            .collect();
        let mut result = solve_factored(&self.lower, normalized, self.n)?;
        for (value, scale) in result.iter_mut().zip(&self.scale) {
            *value /= scale;
        }
        self.admit_vector(&result)?;
        Ok(result)
    }

    /// Apply the assembled physical mass to a generalized velocity.
    /// The first six entries are body-axis linear momentum and angular
    /// momentum about the declared base origin, including attached fluid.
    ///
    /// # Errors
    /// Invalid vector shape, nonfinite input or unrepresentable momentum.
    pub fn momentum(&self, velocity: &[f64]) -> Result<Vec<f64>, Refusal> {
        self.admit_vector(velocity)?;
        let mut momentum = vec![0.0; self.n];
        for i in 0..self.n {
            for j in 0..self.n {
                momentum[i] += self.mass[i * self.n + j] * velocity[j];
            }
        }
        self.admit_vector(&momentum)?;
        Ok(momentum)
    }

    /// Evaluate `u^T M u/2` through the positive factored kinetic form.
    /// Computing `|L^T D u|^2/2` avoids subtracting large orbital cross terms.
    ///
    /// # Errors
    /// Invalid vector shape, nonfinite input or unrepresentable kinetic energy.
    pub fn kinetic_energy(&self, velocity: &[f64]) -> Result<f64, Refusal> {
        self.admit_vector(velocity)?;
        let mut energy = 0.0;
        for j in 0..self.n {
            let mut mode = 0.0;
            for i in j..self.n {
                mode += self.lower[i * self.n + j] * (self.scale[i] * velocity[i]);
            }
            energy += 0.5 * mode * mode;
        }
        if !energy.is_finite() {
            return Err(refuse(
                "joint-energy-invalid",
                "unrepresentable generalized kinetic energy".into(),
                "inspect mass, geometry and velocity scales",
            ));
        }
        Ok(energy)
    }
}

/// Extra diagonal terms arise only from eliminating implicit linear impedance.
pub(crate) fn solve_articulated_with_diagonal(
    dry: &FreeFloatingTerms,
    loads: &AeroGeneralizedLoads,
    control_diagonal: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), Refusal> {
    let mass = PreparedArticulatedMass::with_control_diagonal(dry, loads, control_diagonal)?;
    let mut rhs = vec![0.0; mass.coordinate_count()];
    for (i, entry) in rhs.iter_mut().enumerate() {
        *entry = if i < 6 {
            loads.q_rigid_nonaccel[i]
        } else {
            loads.q_control_nonaccel[i - 6]
        } + loads.q_added_bias[i]
            - dry.bias[native_index(i)];
    }
    let acceleration = mass.solve(&rhs)?;
    Ok((acceleration[..6].to_vec(), acceleration[6..].to_vec()))
}

pub(crate) fn solve_positive(
    a: Vec<f64>,
    rhs: Vec<f64>,
    n: usize,
    pivot_floor: f64,
) -> Result<Vec<f64>, Refusal> {
    let lower = factor_positive(a, n, pivot_floor)?;
    solve_factored(&lower, rhs, n)
}

fn factor_positive(mut a: Vec<f64>, n: usize, pivot_floor: f64) -> Result<Vec<f64>, Refusal> {
    // The original fixed-order factorization arithmetic.
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                if !sum.is_finite() || sum <= pivot_floor {
                    return Err(refuse(
                        "joint-mass-not-pd",
                        format!("Cholesky pivot {sum:e} at row {i}"),
                        "the total effective mass must be positive and numerically resolved",
                    ));
                }
                a[i * n + j] = sum.sqrt();
            } else {
                a[i * n + j] = sum / a[j * n + j];
            }
        }
    }
    Ok(a)
}

fn solve_factored(a: &[f64], mut rhs: Vec<f64>, n: usize) -> Result<Vec<f64>, Refusal> {
    for i in 0..n {
        let mut sum = rhs[i];
        for k in 0..i {
            sum -= a[i * n + k] * rhs[k];
        }
        rhs[i] = sum / a[i * n + i];
    }
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for k in i + 1..n {
            sum -= a[k * n + i] * rhs[k];
        }
        rhs[i] = sum / a[i * n + i];
    }
    if !rhs.iter().all(|x| x.is_finite()) {
        return Err(refuse(
            "joint-mass-not-pd",
            "nonfinite generalized acceleration".into(),
            "inspect inertia and applied load scales",
        ));
    }
    Ok(rhs)
}
