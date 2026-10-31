//! Fluid inertia attached to the original articulated body's moving links.
//!
//! A strip is constant in its link frame. The native articulated owner supplies
//! the link twist `V = J u` and acceleration bias `c = Jdot u`; its same map
//! serves dry and fluid inertia. The existing positive strip modes give a
//! constant local metric `G`, so the fluid contribution is `Jᵀ G J` and its
//! right-hand-side bias is `−Jᵀ (G c + ad*_V G V)`. This includes serial moving
//! hinges, sliding/helical joints, all control cross blocks and their analytic
//! configuration terms without differentiating assembled matrices.
//!
//! Only this adapter permutes native `[angular, linear, joints]` columns into
//! Flyer's `[linear, angular, controls]` order. Cross-surface fluid interaction,
//! deforming strips and density time derivatives remain outside this model.

use super::{
    AeroGeneralizedLoads, MAX_CONTROL_COORDS, MAX_STRIPS, Strip, admit_strip, refuse, strip_modes,
};
use crate::Refusal;
use fs_ga::{Vec3, Wrench};
use fs_mbd::articulated::generalized::GeneralizedKinematics;

#[derive(Clone, Copy, Debug, PartialEq)]
struct LocalMode {
    inertia_per_density: f64,
    covector: Wrench,
}

/// A finite strip rigidly attached to one native articulated link.
#[derive(Clone, Debug, PartialEq)]
pub struct LinkedStrip {
    link: usize,
    modes: [LocalMode; 3],
}

impl LinkedStrip {
    /// Bind strip geometry expressed in the link's own coordinates.
    ///
    /// Control fields must be inactive: all control velocities come from the
    /// articulated model, including inherited motion of upstream joints. Link
    /// existence is checked against the shared kinematics at assembly time.
    ///
    /// # Errors
    /// Invalid strip geometry or any preassigned control velocity.
    pub fn new(link: usize, strip_in_link: Strip) -> Result<Self, Refusal> {
        admit_strip(&strip_in_link, 0)?;
        if strip_in_link.control_coord.is_some()
            || strip_in_link.control_gain != 0.0
            || strip_in_link.control_angular_gain != [0.0; 3]
        {
            return Err(refuse(
                "linked-strip-invalid",
                "linked strip has preassigned control velocity".into(),
                "supply inactive control fields; the native link map owns all joint motion",
            ));
        }
        let modes = strip_modes(1.0, &strip_in_link, 0)?.modes.map(|mode| {
            let j = mode.rigid;
            LocalMode {
                inertia_per_density: mode.inertia,
                covector: Wrench::new(Vec3::new(j[3], j[4], j[5]), Vec3::new(j[0], j[1], j[2])),
            }
        });
        Ok(Self { link, modes })
    }

    /// Index in the native model's root-first link order.
    #[must_use]
    pub fn link(&self) -> usize {
        self.link
    }
}

fn native_column(flyer_column: usize) -> usize {
    if flyer_column < 6 {
        (flyer_column + 3) % 6
    } else {
        flyer_column
    }
}

/// Assemble fluid blocks and their full bias using already admitted kinematics.
///
/// Fixed/root links and arbitrarily nested scalar joints use the same strip
/// law. The output's applied-load fields are zero; its bias contains only fluid
/// inertia. Sum it with the complete dry generalized inertia and dry bias from
/// the same native state before solving; adding a constant body and independent
/// scalar control inertias would omit the moving dry linkage's cross terms.
///
/// The native map is prepared once and reused. Each strip's admitted local
/// rank modes are cached at construction, then scaled by density and projected
/// using one bounded stack row. There is no per-strip matrix or numerical
/// configuration derivative. Density is held constant in the bias calculation.
///
/// # Errors
/// Nonpositive/nonfinite density, control or strip caps, missing links, or
/// unrepresentable assembled mass/bias.
pub fn assemble_articulated_strips(
    rho_kg_m3: f64,
    strips: &[LinkedStrip],
    kinematics: &GeneralizedKinematics,
) -> Result<AeroGeneralizedLoads, Refusal> {
    let nc = kinematics.dof_count();
    if nc > MAX_CONTROL_COORDS {
        return Err(refuse(
            "control-coords-outside-domain",
            format!("nc = {nc} exceeds the cap {MAX_CONTROL_COORDS}"),
            "use at most the admitted number of active dynamic coordinates",
        ));
    }
    if strips.len() > MAX_STRIPS {
        return Err(refuse(
            "strip-count-exceeded",
            format!("{} strips exceed {MAX_STRIPS}", strips.len()),
            "aggregate panels into strips",
        ));
    }
    if !rho_kg_m3.is_finite() || rho_kg_m3 <= 0.0 {
        return Err(refuse(
            "non-finite-input",
            "fluid density must be finite and positive".into(),
            "supply the admitted fluid density",
        ));
    }
    if strips.iter().any(|s| s.link >= kinematics.link_count()) {
        return Err(refuse(
            "linked-strip-invalid",
            "strip refers to a missing articulated link".into(),
            "bind strips to links in the same native articulated model",
        ));
    }

    let mut result = AeroGeneralizedLoads {
        q_rigid_nonaccel: [0.0; 6],
        q_control_nonaccel: vec![0.0; nc],
        m_added_rr: [[0.0; 6]; 6],
        m_added_rc: vec![[0.0; 6]; nc],
        m_added_cc: vec![vec![0.0; nc]; nc],
        q_added_bias: vec![0.0; 6 + nc],
    };
    for strip in strips {
        let velocity = kinematics.body_twist(strip.link);
        let acceleration = kinematics.body_bias_acceleration(strip.link);
        let jacobian = kinematics.body_jacobian(strip.link);
        for mode in strip.modes {
            let inertia = rho_kg_m3 * mode.inertia_per_density;
            if !inertia.is_finite() || inertia <= 0.0 {
                return Err(refuse(
                    "linked-strip-invalid",
                    "unrepresentable strip inertia at the supplied density".into(),
                    "use representable density and strip dimensions",
                ));
            }
            let covector = mode.covector;
            let speed = covector.pairing(velocity);
            let bias_speed = covector.pairing(acceleration);
            let mut row = [0.0; 6 + MAX_CONTROL_COORDS];
            for a in 0..6 + nc {
                let column = jacobian[native_column(a)];
                row[a] = covector.pairing(column);
                // <J_a, ad*_V p> = -<p, [V, J_a]>: use the original
                // Lie bracket to include angular AND linear transport once.
                result.q_added_bias[a] += inertia
                    * (speed * covector.pairing(velocity.bracket(column)) - row[a] * bias_speed);
            }
            for a in 0..6 {
                for b in 0..6 {
                    result.m_added_rr[a][b] += inertia * (row[a] * row[b]);
                }
            }
            for a in 0..nc {
                for b in 0..6 {
                    result.m_added_rc[a][b] += inertia * (row[6 + a] * row[b]);
                }
                for b in 0..nc {
                    result.m_added_cc[a][b] += inertia * (row[6 + a] * row[6 + b]);
                }
            }
        }
    }
    if !result
        .m_added_rr
        .iter()
        .flatten()
        .chain(result.m_added_rc.iter().flatten())
        .chain(result.m_added_cc.iter().flatten())
        .chain(&result.q_added_bias)
        .all(|x| x.is_finite())
    {
        return Err(refuse(
            "linked-strip-invalid",
            "unrepresentable fluid inertia or transport bias".into(),
            "use representable geometry, density and joint motion",
        ));
    }
    Ok(result)
}
