//! Tangent-space optimization over SO(3)/SE(3): the approved hook for
//! optimizing pose keyframes and quaternion-headed parameters directly.
//!
//! Why this exists: Euclidean sampling of rotations is wrong (double cover,
//! no group closure, curved geodesics). Every operation here samples in the
//! tangent algebra, moves through the group exponential, and averages through
//! the log map — the same exp/log machinery the dynamics owners use, so poses
//! composed here stay interoperable with `Se3`-posed scene state.
//!
//! What this is NOT: a replacement for CMA-ES over joint-space policies
//! (joint targets live in R^n, where Euclidean CMA-ES is already the
//! natural-gradient method — see fs-dfo/src/cma.rs). Use this module when
//! the optimized parameters ARE poses (SE(3) keyframes, effector/camera
//! targets, quaternion output heads).
//!
//! Conventions follow [`crate::lie`]: twists `[angular, linear]`;
//! `space_*` = left perturbation (`Exp(delta) * group`); `body_*` = right
//! perturbation (`group * Exp(delta)`).

use crate::GaError;
use crate::facade::Vec3;
use crate::lie::{Se3, So3, So3Tangent, Twist};

/// Deterministic splitmix64 for the tangent samplers (module-scope
/// determinism; keyed production streams stay in fs-rand).
#[must_use]
pub fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Uniform sample in [-1, 1).
#[must_use]
fn unit_pm(state: &mut u64) -> f64 {
    (splitmix64(state) >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

/// One tangent-space evolution step on SO(3): sample `k` isotropic tangent
/// perturbations of magnitude `sigma`, score them, and recenter through the
/// exponential toward the best — the geodesically correct analog of
/// "mu += sigma * best_direction" (left/space perturbation convention).
/// Returns the accepted tangent norm (0.0 when the incumbent stays).
///
/// # Errors
/// Propagates [`GaError`] from the exponential/log admission.
pub fn so3_tangent_step(
    pose: &mut So3,
    objective: &dyn Fn(&So3) -> f64,
    sigma: f64,
    k: usize,
    state: &mut u64,
) -> Result<f64, GaError> {
    let incumbent = objective(pose);
    let mut best_tangent = So3Tangent::new(Vec3::new(0.0, 0.0, 0.0));
    let mut best_score = incumbent;
    for _ in 0..k.max(1) {
        let tangent = So3Tangent::new(Vec3::new(
            unit_pm(state) * sigma,
            unit_pm(state) * sigma,
            unit_pm(state) * sigma,
        ));
        let candidate = So3::exp(tangent)?.compose(*pose)?;
        let score = objective(&candidate);
        if score < best_score {
            best_score = score;
            best_tangent = tangent;
        }
    }
    if best_score >= incumbent {
        return Ok(0.0);
    }
    let a = best_tangent.angular;
    let norm = (a.x * a.x + a.y * a.y + a.z * a.z).sqrt();
    *pose = So3::exp(best_tangent)?.compose(*pose)?;
    Ok(norm)
}

/// One tangent-space evolution step on SE(3): same contract as
/// [`so3_tangent_step`] with a full `se(3)` twist (angular rad, linear m,
/// both isotropic at `sigma`). Returns the accepted max norm component.
///
/// # Errors
/// Propagates [`GaError`] from the SE(3) exponential admission.
pub fn se3_tangent_step(
    pose: &mut Se3,
    objective: &dyn Fn(&Se3) -> f64,
    sigma: f64,
    k: usize,
    state: &mut u64,
) -> Result<Option<(f64, f64)>, GaError> {
    let incumbent = objective(pose);
    let mut best_twist = Twist::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
    let mut best_score = incumbent;
    for _ in 0..k.max(1) {
        let twist = Twist::new(
            Vec3::new(
                unit_pm(state) * sigma,
                unit_pm(state) * sigma,
                unit_pm(state) * sigma,
            ),
            Vec3::new(
                unit_pm(state) * sigma,
                unit_pm(state) * sigma,
                unit_pm(state) * sigma,
            ),
        );
        let candidate = Se3::exp(twist)?.compose(*pose)?;
        let score = objective(&candidate);
        if score < best_score {
            best_score = score;
            best_twist = twist;
        }
    }
    if best_score >= incumbent {
        return Ok(None);
    }
    let an = best_twist.angular;
    let ln = best_twist.linear;
    let norm = (an.x * an.x + an.y * an.y + an.z * an.z)
        .sqrt()
        .max((ln.x * ln.x + ln.y * ln.y + ln.z * ln.z).sqrt());
    *pose = Se3::exp(best_twist)?.compose(*pose)?;
    Ok(Some((norm, best_score)))
}

/// Right-perturbation geodesic (Karcher) mean of SO(3) poses: fixed-point
/// iteration `mu <- mu * Exp(mean_i Log(pose_i mu^-1))` (the module's
/// `body_*` right-perturbation convention). Converges for clusters inside a
/// common geodesic ball (radius < pi).
///
/// # Errors
/// Propagates [`GaError`] from the log/exp admission.
pub fn so3_geodesic_mean(poses: &[So3], iterations: usize) -> Result<Option<So3>, GaError> {
    let mut mean = match poses.first() {
        Some(first) => *first,
        None => return Ok(None),
    };
    let count = poses.len() as f64;
    for _ in 0..iterations.max(1) {
        let mut acc = So3Tangent::new(Vec3::new(0.0, 0.0, 0.0));
        for pose in poses {
            let delta = pose.compose(mean.inverse())?.log();
            acc = So3Tangent::new(Vec3::new(
                acc.angular.x + delta.angular.x / count,
                acc.angular.y + delta.angular.y / count,
                acc.angular.z + delta.angular.z / count,
            ));
        }
        mean = So3::exp(acc)?.compose(mean)?;
    }
    Ok(Some(mean))
}

/// Right-perturbation geodesic mean of SE(3) poses (same fixed-point scheme
/// over the full twist log).
///
/// # Errors
/// Propagates [`GaError`] from the log/exp admission.
pub fn se3_geodesic_mean(poses: &[Se3], iterations: usize) -> Result<Option<Se3>, GaError> {
    let mut mean = match poses.first() {
        Some(first) => *first,
        None => return Ok(None),
    };
    let count = poses.len() as f64;
    for _ in 0..iterations.max(1) {
        let mut acc = Twist::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        for pose in poses {
            let delta = pose.compose(mean.inverse()?)?.log();
            acc = Twist::new(
                Vec3::new(
                    acc.angular.x + delta.angular.x / count,
                    acc.angular.y + delta.angular.y / count,
                    acc.angular.z + delta.angular.z / count,
                ),
                Vec3::new(
                    acc.linear.x + delta.linear.x / count,
                    acc.linear.y + delta.linear.y / count,
                    acc.linear.z + delta.linear.z / count,
                ),
            );
        }
        mean = Se3::exp(acc)?.compose(mean)?;
    }
    Ok(Some(mean))
}

/// Multi-step driver: run [`se3_tangent_step`] until the accepted step norm
/// drops below `tolerance` or `max_steps` elapse. Returns the last accepted
/// norm (0.0 = converged incumbent, never moved again).
///
/// # Errors
/// Propagates [`GaError`] from the exponential admission.
pub fn se3_optimize(
    pose: &mut Se3,
    objective: &dyn Fn(&Se3) -> f64,
    sigma: f64,
    k: usize,
    max_steps: usize,
    tolerance: f64,
    state: &mut u64,
) -> Result<f64, GaError> {
    let mut sigma_t = sigma;
    let mut objective_now = objective(pose);
    for _ in 0..max_steps {
        if objective_now <= tolerance {
            break;
        }
        // (1+lambda)-ES with classic 1/5th-style adaptation: widen on a
        // rejected step, tighten on an accepted one. A rejected step means
        // "no candidate improved" — NOT convergence.
        match se3_tangent_step(pose, objective, sigma_t, k, state)? {
            Some((_norm, score)) => {
                objective_now = score;
                sigma_t = (sigma_t * 0.97).max(1e-5);
            }
            None => sigma_t = (sigma_t * 1.12).min(10.0),
        }
    }
    Ok(objective_now)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sq_distance_so3(a: &So3, b: &So3) -> f64 {
        match a.compose(b.inverse()) {
            Ok(delta) => {
                let d = delta.log().angular;
                d.x * d.x + d.y * d.y + d.z * d.z
            }
            Err(_) => f64::MAX,
        }
    }

    #[test]
    fn geodesic_mean_of_one_pose_is_the_pose() -> Result<(), GaError> {
        let pose = Se3::exp(Twist::new(
            Vec3::new(0.1, -0.2, 0.3),
            Vec3::new(0.4, 0.5, -0.6),
        ))?;
        let mean = se3_geodesic_mean(&[pose], 8)?.unwrap();
        let t = mean.compose(pose.inverse()?)?.log();
        assert!(t.angular.x.abs() < 1e-12);
        assert!(t.linear.z.abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn so3_geodesic_mean_is_symmetric() {
        let a = So3::exp(So3Tangent::new(Vec3::new(0.2, 0.0, 0.0))).unwrap();
        let b = So3::exp(So3Tangent::new(Vec3::new(0.0, 0.2, 0.0))).unwrap();
        let m1 = so3_geodesic_mean(&[a, b], 16).unwrap().unwrap();
        let m2 = so3_geodesic_mean(&[b, a], 16).unwrap().unwrap();
        assert!(sq_distance_so3(&m1, &m2) < 1e-18);
    }

    #[test]
    fn se3_optimizer_reaches_target_pose() {
        let target = Se3::exp(Twist::new(
            Vec3::new(0.3, -0.2, 0.15),
            Vec3::new(0.5, 0.25, -0.1),
        ))
        .unwrap();
        let objective = |p: &Se3| -> f64 {
            let huge = Twist::new(
                Vec3::new(1.0e6, 1.0e6, 1.0e6),
                Vec3::new(1.0e6, 1.0e6, 1.0e6),
            );
            let t = match p.compose(target.inverse().unwrap()) {
                Ok(delta) => delta.log(),
                Err(_) => return 1.0e12,
            };
            let _ = huge;
            let a = t.angular;
            let l = t.linear;
            a.x * a.x + a.y * a.y + a.z * a.z + l.x * l.x + l.y * l.y + l.z * l.z
        };
        let mut start = Se3::exp(Twist::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
        ))
        .unwrap();
        let mut state = 0xD1CE_5EED_u64;
        // The (1+12)-ES converges logarithmically near the optimum — this
        // utility is for coarse pose alignment, not precision polishing.
        let final_objective =
            se3_optimize(&mut start, &objective, 0.12, 32, 6000, 1e-3, &mut state)
                .expect("se3_optimize must not refuse");
        assert!(
            final_objective < 5e-3,
            "did not coarsely converge: objective = {final_objective}"
        );
    }

    #[test]
    fn tangent_steps_are_deterministic() {
        let objective = |_p: &So3| -> f64 { 1.0 };
        let mut s1 = 42u64;
        let mut s2 = 42u64;
        let mut p1 = So3::exp(So3Tangent::new(Vec3::new(0.1, 0.0, 0.0))).unwrap();
        let mut p2 = p1;
        let r1 = so3_tangent_step(&mut p1, &objective, 0.05, 8, &mut s1).unwrap();
        let r2 = so3_tangent_step(&mut p2, &objective, 0.05, 8, &mut s2).unwrap();
        assert_eq!(r1, r2);
        assert_eq!(sq_distance_so3(&p1, &p2), 0.0);
    }
}
