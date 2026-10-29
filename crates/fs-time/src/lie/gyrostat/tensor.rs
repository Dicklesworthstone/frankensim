//! Fixed body-frame inertia with a cached principal-axis representation.
//!
//! A body need not use its principal axes as geometry/control coordinates.
//! Diagonalization belongs to assembly admission, not the stepping loop. This
//! adapter reuses the native gyrostat map after a fixed change of body basis;
//! it does not introduce another rotational equation or per-step solver.

use super::{HeldTorque, LieStepError, driven_step, momentum_step};
use fs_ga::{Mat3, Quat, So3, Vec3};

const SYMMETRY_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const RECONSTRUCTION_TOLERANCE: f64 = 512.0 * f64::EPSILON;

/// Positive-definite rotational kinetic form in fixed body coordinates.
///
/// Construction checks a measured normalized reconstruction residual. It is
/// ordinary floating-point admission, not an interval eigenvalue certificate
/// or evidence that a particular material geometry realizes the supplied
/// inertia. Reduced rotor inertias are allowed. The tensor is fixed during a
/// step: changing mass, moving rotor axes and deforming bodies need their own
/// momentum transport law.
///
/// [`Self::new`] uses the existing bounded 3-by-3 native eigensolver once;
/// [`Self::from_principal`] needs no eigensolve or allocation. All subsequent
/// transforms and integration calls are allocation-free and constant work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertiaTensor3 {
    body: Mat3,
    principal: Vec3,
    principal_to_body: So3,
    relative_reconstruction_error: f64,
}

impl InertiaTensor3 {
    /// Admit a full body-frame inertia and cache its principal frame.
    ///
    /// Off-diagonal asymmetry beyond 64 machine epsilons of the largest
    /// entry refuses. Roundoff-sized asymmetry is averaged explicitly.
    /// Strictly diagonal input retains its supplied axis order and delegates
    /// to the old diagonal map bit for bit. Otherwise the native eigensolver
    /// supplies ascending principal moments with deterministic tie-breaking.
    ///
    /// # Errors
    /// Refuses non-finite, zero, nonsymmetric or non-positive-definite forms,
    /// unrepresentable principal moments, and failed reconstruction/group
    /// admission. Near-singular forms have numerical, not certified, admission.
    pub fn new(body: Mat3) -> Result<Self, LieStepError> {
        if !body.m.iter().all(|x| x.is_finite()) {
            return Err(invalid("body inertia tensor"));
        }
        let scale = body.max_abs();
        if scale == 0.0 {
            return Err(invalid("body inertia tensor"));
        }
        let mut normalized = body.m.map(|x| x / scale);
        for (i, j) in [(1, 3), (2, 6), (5, 7)] {
            if (normalized[i] - normalized[j]).abs() > SYMMETRY_TOLERANCE {
                return Err(invalid("body inertia symmetry"));
            }
            let average = 0.5 * normalized[i] + 0.5 * normalized[j];
            normalized[i] = average;
            normalized[j] = average;
        }
        if [1, 2, 3, 5, 6, 7].iter().all(|&i| normalized[i] == 0.0) {
            return Self::from_principal(
                Vec3::new(body.m[0], body.m[4], body.m[8]),
                So3::identity(),
            );
        }
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&normalized, 3);
        if !values.iter().all(|x| x.is_finite() && *x > 0.0)
            || !vectors.iter().all(|x| x.is_finite())
        {
            return Err(invalid("positive-definite body inertia"));
        }
        let mut basis = [0.0; 9];
        basis.copy_from_slice(&vectors);
        // Eigenvector signs have no physical meaning. Fix each by its largest
        // component (lowest row wins ties), then choose a proper third axis.
        for column in 0..3 {
            let mut pivot = column;
            for row in 1..3 {
                let candidate = 3 * row + column;
                if basis[candidate].abs() > basis[pivot].abs() {
                    pivot = candidate;
                }
            }
            if basis[pivot] < 0.0 {
                for row in 0..3 {
                    basis[3 * row + column] = -basis[3 * row + column];
                }
            }
        }
        let a = Vec3::new(basis[0], basis[3], basis[6]);
        let b = Vec3::new(basis[1], basis[4], basis[7]);
        let c = Vec3::new(basis[2], basis[5], basis[8]);
        if a.cross(b).dot(c) < 0.0 {
            for row in 0..3 {
                basis[3 * row + 2] = -basis[3 * row + 2];
            }
        }
        let principal_to_body = rotation_from_matrix(Mat3 { m: basis })?;
        let principal = Vec3::new(values[0] * scale, values[1] * scale, values[2] * scale);
        admit_principal(principal)?;
        let reconstructed = reconstruct(
            Vec3::new(values[0], values[1], values[2]),
            principal_to_body,
        );
        let error = normalized
            .iter()
            .zip(reconstructed.m)
            .map(|(expected, actual)| (expected - actual).abs())
            .fold(0.0_f64, f64::max);
        if !error.is_finite() || error > RECONSTRUCTION_TOLERANCE {
            return Err(invalid("body inertia principal-frame reconstruction"));
        }
        Ok(Self {
            body: Mat3 {
                m: normalized.map(|x| x * scale),
            },
            principal,
            principal_to_body,
            relative_reconstruction_error: error,
        })
    }

    /// Construct from already known positive principal moments and a proper
    /// rotation mapping principal coordinates into the body's coordinates.
    ///
    /// # Errors
    /// Refuses non-positive/non-finite moments or an unrepresentable tensor.
    pub fn from_principal(principal: Vec3, principal_to_body: So3) -> Result<Self, LieStepError> {
        admit_principal(principal)?;
        let body = reconstruct(principal, principal_to_body);
        if !body.m.iter().all(|x| x.is_finite()) {
            return Err(invalid("body inertia reconstruction"));
        }
        Ok(Self {
            body,
            principal,
            principal_to_body,
            relative_reconstruction_error: 0.0,
        })
    }

    /// Admitted symmetric tensor in the caller's fixed body coordinates.
    #[must_use]
    pub const fn matrix(self) -> Mat3 {
        self.body
    }

    /// Principal moments used by the native rotation map.
    #[must_use]
    pub const fn principal_moments(self) -> Vec3 {
        self.principal
    }

    /// Fixed coordinate transformation; no time-dependent frame is inferred.
    #[must_use]
    pub const fn principal_to_body(self) -> So3 {
        self.principal_to_body
    }

    /// Measured maximum reconstruction error divided by the largest input
    /// entry. Zero for an explicitly supplied principal representation.
    #[must_use]
    pub const fn relative_reconstruction_error(self) -> f64 {
        self.relative_reconstruction_error
    }

    /// Body angular momentum `I omega`, excluding any separately carried
    /// internal rotor momentum.
    ///
    /// # Errors
    /// Refuses non-finite input or an unrepresentable result.
    pub fn momentum(self, omega_body: Vec3) -> Result<Vec3, LieStepError> {
        finite_vector(omega_body, "body angular velocity")?;
        let momentum = self.body.apply(omega_body);
        finite_vector(momentum, "body angular momentum")?;
        Ok(momentum)
    }

    /// Body angular velocity `I^-1 L`, excluding internal rotor momentum
    /// from `L` before calling. Uses the cached principal representation.
    ///
    /// # Errors
    /// Refuses non-finite input or an unrepresentable result.
    pub fn velocity(self, momentum_body: Vec3) -> Result<Vec3, LieStepError> {
        let p = self.principal_to_body.inverse().rotate(momentum_body)?;
        let omega = Vec3::new(
            p.x / self.principal.x,
            p.y / self.principal.y,
            p.z / self.principal.z,
        );
        let body = self.principal_to_body.rotate(omega)?;
        finite_vector(body, "body angular velocity")?;
        Ok(body)
    }

    /// Free motion with total body momentum `L = I omega + h_rotor`.
    ///
    /// The body pose, total momentum and internal momentum use the caller's
    /// frame. The same fixed basis transforms all three into the original
    /// gyrostat integrator and back. No nonlinear solve or allocation occurs.
    /// This inherits that map's second-order, fixed-step approximation and
    /// roundoff-level spatial angular-momentum balance, not exact energy.
    ///
    /// # Errors
    /// The refusals of [`super::momentum_step`] and group transformations.
    pub fn momentum_step(
        self,
        rotation: So3,
        momentum_body: Vec3,
        rotor_momentum_body: Vec3,
        dt: f64,
    ) -> Result<(So3, Vec3), LieStepError> {
        if self.principal_to_body == So3::identity() {
            return momentum_step(
                rotation,
                momentum_body,
                self.principal,
                rotor_momentum_body,
                dt,
            );
        }
        let basis = self.principal_to_body;
        let inverse = basis.inverse();
        let (increment, _) = momentum_step(
            So3::identity(),
            inverse.rotate(momentum_body)?,
            self.principal,
            inverse.rotate(rotor_momentum_body)?,
            dt,
        )?;
        // Conjugate the *increment*, not the accumulated attitude. Relabeling
        // the quaternion vector avoids two constant-frame quaternion products
        // per step. Transport the original body momentum with this same
        // increment, avoiding repeated full-state basis-roundtrip drift.
        let q = increment.as_quat();
        let v = basis.rotate(Vec3::new(q.x, q.y, q.z))?;
        let body_increment = So3::try_from_quat(Quat {
            w: q.w,
            x: v.x,
            y: v.y,
            z: v.z,
        })?;
        Ok((
            rotation.compose(body_increment)?,
            body_increment.inverse().rotate(momentum_body)?,
        ))
    }

    /// Held external torque and linear internal-momentum transfer in the
    /// caller's body frame. `HeldTorque::Space` remains in the inertial frame.
    ///
    /// # Errors
    /// The refusals of [`super::driven_step`] and group transformations.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_step(
        self,
        rotation: So3,
        omega_body: Vec3,
        rotor_before: Vec3,
        rotor_after: Vec3,
        torque: HeldTorque,
        dt: f64,
    ) -> Result<(So3, Vec3), LieStepError> {
        if self.principal_to_body == So3::identity() {
            return driven_step(
                rotation,
                omega_body,
                self.principal,
                rotor_before,
                rotor_after,
                torque,
                dt,
            );
        }
        super::admit_rotor(rotor_before)?;
        super::admit_rotor(rotor_after)?;
        let midpoint_rotor = rotor_before.scale(0.5) + rotor_after.scale(0.5);
        let torque_at = |orientation: So3| -> Result<Vec3, LieStepError> {
            match torque {
                HeldTorque::Body(value) => Ok(value),
                HeldTorque::Space(value) => Ok(orientation.inverse().rotate(value)?),
            }
        };
        let initial = self.momentum(omega_body)? + rotor_before;
        let kicked = super::super::torque_kick(initial, torque_at(rotation)?, dt * 0.5)?;
        let (next_rotation, drifted) = self.momentum_step(rotation, kicked, midpoint_rotor, dt)?;
        let accepted = super::super::torque_kick(drifted, torque_at(next_rotation)?, dt * 0.5)?;
        Ok((next_rotation, self.velocity(accepted - rotor_after)?))
    }
}

fn invalid(context: &'static str) -> LieStepError {
    LieStepError::InvalidParameter { context }
}

fn finite_vector(value: Vec3, context: &'static str) -> Result<(), LieStepError> {
    if [value.x, value.y, value.z].iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(invalid(context))
    }
}

fn admit_principal(value: Vec3) -> Result<(), LieStepError> {
    if [value.x, value.y, value.z]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
    {
        Ok(())
    } else {
        Err(invalid("positive-definite body inertia"))
    }
}

fn reconstruct(principal: Vec3, basis: So3) -> Mat3 {
    let r = basis.matrix().m;
    let d = [principal.x, principal.y, principal.z];
    let mut m = [0.0; 9];
    for row in 0..3 {
        for column in row..3 {
            let mut value = 0.0;
            for axis in 0..3 {
                value += (r[3 * row + axis] * r[3 * column + axis]) * d[axis];
            }
            m[3 * row + column] = value;
            m[3 * column + row] = value;
        }
    }
    Mat3 { m }
}

fn rotation_from_matrix(matrix: Mat3) -> Result<So3, LieStepError> {
    let m = matrix.m;
    // Select the largest quaternion component to avoid division by a small
    // trace at half-turns. These are direct matrix-to-quaternion formulas;
    // So3 validates the result without normalization or a polar projection.
    let squared = [
        1.0 + m[0] + m[4] + m[8],
        1.0 + m[0] - m[4] - m[8],
        1.0 - m[0] + m[4] - m[8],
        1.0 - m[0] - m[4] + m[8],
    ];
    let mut largest = 0;
    for i in 1..4 {
        if squared[i] > squared[largest] {
            largest = i;
        }
    }
    let component = 0.5 * squared[largest].sqrt();
    let denominator = 4.0 * component;
    let q = match largest {
        0 => Quat {
            w: component,
            x: (m[7] - m[5]) / denominator,
            y: (m[2] - m[6]) / denominator,
            z: (m[3] - m[1]) / denominator,
        },
        1 => Quat {
            w: (m[7] - m[5]) / denominator,
            x: component,
            y: (m[1] + m[3]) / denominator,
            z: (m[2] + m[6]) / denominator,
        },
        2 => Quat {
            w: (m[2] - m[6]) / denominator,
            x: (m[1] + m[3]) / denominator,
            y: component,
            z: (m[5] + m[7]) / denominator,
        },
        _ => Quat {
            w: (m[3] - m[1]) / denominator,
            x: (m[2] + m[6]) / denominator,
            y: (m[5] + m[7]) / denominator,
            z: component,
        },
    };
    Ok(So3::try_from_quat(q)?)
}
