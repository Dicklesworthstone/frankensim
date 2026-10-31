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

/// Extra diagonal terms arise only from eliminating implicit linear impedance.
pub(crate) fn solve_articulated_with_diagonal(
    dry: &FreeFloatingTerms,
    loads: &AeroGeneralizedLoads,
    control_diagonal: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), Refusal> {
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
            "dry generalized terms or implicit diagonal do not match the active controls".into(),
            "prepare native dry terms and fluid blocks at the same admitted state",
        ));
    }
    let mut a = vec![0.0; n * n];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        let ni = native_index(i);
        rhs[i] = if i < 6 {
            loads.q_rigid_nonaccel[i]
        } else {
            loads.q_control_nonaccel[i - 6]
        } + loads.q_added_bias[i]
            - dry.bias[ni];
        for j in 0..n {
            let nj = native_index(j);
            admit_symmetric(dry.mass[ni][nj], dry.mass[nj][ni])?;
            let fluid = match (i < 6, j < 6) {
                (true, true) => loads.m_added_rr[i][j],
                (true, false) => loads.m_added_rc[j - 6][i],
                (false, true) => loads.m_added_rc[i - 6][j],
                (false, false) => loads.m_added_cc[i - 6][j - 6],
            };
            a[i * n + j] = dry.mass[ni][nj] + fluid;
        }
        if i >= 6 && !control_diagonal.is_empty() {
            a[i * n + i] += control_diagonal[i - 6];
        }
    }
    let scale: Vec<_> = (0..n).map(|i| a[i * n + i].sqrt()).collect();
    if !scale.iter().all(|s| s.is_finite() && *s > 0.0) {
        return Err(refuse(
            "joint-mass-not-pd",
            "invalid effective mass diagonal".into(),
            "use representable positive dry and fluid inertia",
        ));
    }
    for i in 0..n {
        rhs[i] /= scale[i];
        for j in 0..n {
            a[i * n + j] = (a[i * n + j] / scale[i]) / scale[j];
        }
    }
    let mut acceleration = solve_positive(a, rhs, n, 128.0 * f64::EPSILON)?;
    for (value, scale) in acceleration.iter_mut().zip(scale) {
        *value /= scale;
    }
    if !acceleration.iter().all(|v| v.is_finite()) {
        return Err(refuse(
            "joint-mass-not-pd",
            "unrepresentable generalized acceleration".into(),
            "inspect inertia and applied load scales",
        ));
    }
    Ok((acceleration[..6].to_vec(), acceleration[6..].to_vec()))
}

pub(crate) fn solve_positive(
    mut a: Vec<f64>,
    mut rhs: Vec<f64>,
    n: usize,
    pivot_floor: f64,
) -> Result<Vec<f64>, Refusal> {
    // The original fixed-order factorization and substitution arithmetic.
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
