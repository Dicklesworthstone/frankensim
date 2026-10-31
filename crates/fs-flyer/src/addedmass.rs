//! Joint (6+nc) generalized added-mass assembly (bead wf-root-guzez.4.3,
//! E3.2a). Plan §5.1.2, Round-3 S-01 correction: the solve is the JOINT
//! generalized system
//!
//!   [ M_rigid + M_rr   M_rc            ] [ ν̇     ]   [ Q_r ]
//!   [ M_rcᵀ            M_ctrl + M_cc   ] [ q̈_ctrl ] = [ Q_c ] + Q_bias
//!
//! — never the rigid six-vector with orphaned cross blocks. One declared
//! deterministic factorization (dense Cholesky of the PD total, fixed
//! order) BEFORE the Lie-group update; cross-coordinate accelerations are
//! never lagged or finite-differenced.
//!
//! `AnalyticStrip` baseline: each surface strip is a flat plate with 2-D
//! added mass m_a = ρ·π·c²/4 per unit span acting along its normal. The
//! strip's centroid normal-velocity Jacobian j = [n; r×n] (plus the
//! control gain g) contributes m_a·[j; g·e_k][j; g·e_k]ᵀ. Finite chord
//! adds centroidal pitch inertia J_a = m_a c²/32 about the span axis t;
//! exact integration along a uniform span b adds m_a b²/12 about t×n.
//! Each contributes a rotational outer product, including the SAME
//! control-coordinate angular velocity. Thus rigid, control and cross
//! blocks all contain the fluid's translational AND rotational inertia.
//!
//! J_a follows from the 2-D ellipse added inertia πρ(a²−b²)²/8, taking
//! semichord a=c/2 and the thin-plate limit b=0, then multiplying by span.
//! Reference: MIT 13.021 Marine Hydrodynamics, supplemental problem A1,
//! <https://web.mit.edu/13.021/13021_2004/suppnotes/supp04.pdf>, page 3.
//!
//! `Q_added_bias` is the Euler–Poincaré gyroscopic term on the added
//! momentum, −ad*_ν(p): WORKLESS for any p (νᵀ·ad*_ν(p) ≡ 0), which the
//! battery checks to machine precision — energy consistency without any
//! runtime finite-differencing. Configuration-independent strips make the
//! control-block bias identically zero (disclosed below). The fixed-strip
//! entry point freezes geometry and control Jacobians. [`hinged`] supplies
//! finite hinge motion and its analytic configuration-inertia terms through
//! the same rank modes, mass blocks and joint solve.
//!
//! No-claims: cross-surface fluid-inertia interference is OMITTED
//! (AnalyticStrip's disclosed boundary; PanelExtracted is the additive
//! referee mode, E4 lane). Propeller disks are EXCLUDED — rotor terms
//! stay in fs-airscrew.

use crate::Refusal;
use crate::spine::RigidBody;

#[path = "addedmass/hinged.rs"]
pub mod hinged;

#[path = "addedmass/joint.rs"]
pub(crate) mod joint;
pub use joint::{PreparedArticulatedMass, solve_articulated_joint};

#[path = "addedmass/articulated.rs"]
pub mod articulated;

/// Active-dynamic control-coordinate cap (refusals at cap AND cap+1).
pub const MAX_CONTROL_COORDS: usize = 8;
/// Strip-count cap.
pub const MAX_STRIPS: usize = 256;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// One flat-plate strip of the AnalyticStrip baseline.
#[derive(Clone, Debug, PartialEq)]
pub struct Strip {
    /// Name (diagnostics).
    pub name: &'static str,
    /// Chord [m].
    pub chord_m: f64,
    /// Span of this strip [m].
    pub span_m: f64,
    /// Centroid position r [m], frd from the mass-solve reference point.
    pub position_m: [f64; 3],
    /// Unit plate normal (body frame).
    pub normal: [f64; 3],
    /// Unit span direction in the body frame, perpendicular to `normal`.
    /// Its sign is arbitrary; it defines the finite strip's pitch axis.
    pub span_axis: [f64; 3],
    /// Dynamic control coordinate driving this strip, if any.
    pub control_coord: Option<usize>,
    /// Normal-velocity gain g [m/rad] of that coordinate
    /// (hinge arm × mode shape at the centroid).
    pub control_gain: f64,
    /// Relative angular velocity per unit dynamic-coordinate rate [rad/rad],
    /// in body axes. A rigid hinge about the span axis uses `span_axis`.
    /// Together with `control_gain`, this defines velocity at every point
    /// of the finite strip; zero means pure centroid-normal translation.
    pub control_angular_gain: [f64; 3],
}

/// The assembled generalized loads (plan §5.1.2 shape).
#[derive(Clone, Debug, PartialEq)]
pub struct AeroGeneralizedLoads {
    /// Rigid inertial right-hand side [6] (force; moment), body frame.
    /// For body-coordinate acceleration, the caller includes the rigid bias
    /// `(-omega cross m*v, -omega cross I*omega)` with the external loads.
    /// Fluid bias is excluded here: `q_added_bias` supplies it separately.
    pub q_rigid_nonaccel: [f64; 6],
    /// Non-acceleration control loads [nc].
    pub q_control_nonaccel: Vec<f64>,
    /// Added-mass rigid block [6×6].
    pub m_added_rr: [[f64; 6]; 6],
    /// Rigid↔control cross block [6×nc].
    pub m_added_rc: Vec<[f64; 6]>,
    /// Control block [nc×nc].
    pub m_added_cc: Vec<Vec<f64>>,
    /// Energy-consistent bias [6+nc].
    pub q_added_bias: Vec<f64>,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[derive(Clone, Copy)]
struct StripMode {
    inertia: f64,
    rigid: [f64; 6],
    control_gain: f64,
}

struct StripModes {
    modes: [StripMode; 3],
    control_coord: Option<usize>,
    // Each moving strip has one rigid, body-fixed hinge. Its control-mode
    // gains are constant; only the six body-velocity coefficients vary.
    derivative: Option<[[f64; 6]; 3]>,
}

fn admit_strip(s: &Strip, nc: usize) -> Result<(), Refusal> {
    let finite = s.chord_m.is_finite()
        && s.span_m.is_finite()
        && s.position_m.iter().all(|v| v.is_finite())
        && s.normal.iter().all(|v| v.is_finite())
        && s.span_axis.iter().all(|v| v.is_finite())
        && s.control_gain.is_finite()
        && s.control_angular_gain.iter().all(|v| v.is_finite());
    let norm2: f64 = s.normal.iter().map(|v| v * v).sum();
    let span_norm2: f64 = s.span_axis.iter().map(|v| v * v).sum();
    let perpendicular: f64 = s.normal.iter().zip(s.span_axis).map(|(n, t)| n * t).sum();
    if !finite
        || s.chord_m <= 0.0
        || s.span_m <= 0.0
        || (norm2 - 1.0).abs() > 1.0e-9
        || (span_norm2 - 1.0).abs() > 1.0e-9
        || perpendicular.abs() > 1.0e-9
    {
        return Err(refuse(
            "strip-invalid",
            format!("strip {}: geometry or control velocity invalid", s.name),
            "positive chord/span; orthonormal normal and span axes; finite control gains",
        ));
    }
    if let Some(k) = s.control_coord {
        if k >= nc {
            return Err(refuse(
                "strip-invalid",
                format!("strip {}: control coordinate {k} outside nc = {nc}", s.name),
                "declare the coordinate in the active-dynamic set",
            ));
        }
    }
    Ok(())
}

fn strip_modes(rho: f64, s: &Strip, nc: usize) -> Result<StripModes, Refusal> {
    admit_strip(s, nc)?;
    let m_a = rho * core::f64::consts::PI * s.chord_m * s.chord_m / 4.0 * s.span_m;
    let pitch_inertia = m_a * (s.chord_m * s.chord_m / 32.0);
    let span_inertia = m_a * (s.span_m * s.span_m / 12.0);
    if ![m_a, pitch_inertia, span_inertia]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
    {
        return Err(refuse(
            "strip-invalid",
            format!(
                "strip {}: unrepresentable fluid mass or rotary inertia",
                s.name
            ),
            "use representable density and strip dimensions",
        ));
    }
    let rxn = cross(s.position_m, s.normal);
    let translation = StripMode {
        inertia: m_a,
        rigid: [
            s.normal[0],
            s.normal[1],
            s.normal[2],
            rxn[0],
            rxn[1],
            rxn[2],
        ],
        control_gain: s.control_gain,
    };
    let rotary = |inertia, axis: [f64; 3]| StripMode {
        inertia,
        rigid: [0.0, 0.0, 0.0, axis[0], axis[1], axis[2]],
        control_gain: axis
            .iter()
            .zip(s.control_angular_gain)
            .map(|(a, h)| a * h)
            .sum(),
    };
    Ok(StripModes {
        modes: [
            translation,
            rotary(pitch_inertia, s.span_axis),
            rotary(span_inertia, cross(s.span_axis, s.normal)),
        ],
        control_coord: s.control_coord,
        derivative: None,
    })
}

/// Density of the working fluid [kg/m³] enters per assembly call (the
/// provenance-bound air state supplies it).
///
/// Assemble the added-mass blocks + bias for the strip set at the current
/// generalized velocity ν = (v_body[3], ω_body[3]) and control rates.
///
/// # Errors
/// `control-coords-outside-domain` (nc cap at cap AND cap+1),
/// `strip-count-exceeded`, `strip-invalid` (non-finite, non-positive
/// chord/span, nonorthonormal normal/span axes, control coordinate out of range),
/// `non-finite-input`.
pub fn assemble_analytic_strip(
    rho_kg_m3: f64,
    strips: &[Strip],
    nc: usize,
    nu: &[f64; 6],
    control_rates: &[f64],
) -> Result<AeroGeneralizedLoads, Refusal> {
    assemble_modes(
        rho_kg_m3,
        strips.len(),
        nc,
        nu,
        control_rates,
        strips.iter().map(|s| strip_modes(rho_kg_m3, s, nc)),
    )
}

fn assemble_modes(
    rho_kg_m3: f64,
    strip_count: usize,
    nc: usize,
    nu: &[f64; 6],
    control_rates: &[f64],
    strips: impl Iterator<Item = Result<StripModes, Refusal>>,
) -> Result<AeroGeneralizedLoads, Refusal> {
    if nc > MAX_CONTROL_COORDS {
        return Err(refuse(
            "control-coords-outside-domain",
            format!("nc = {nc} exceeds the cap {MAX_CONTROL_COORDS}"),
            "the 1903 Flyer has at most a few active dynamic coordinates",
        ));
    }
    if strip_count > MAX_STRIPS {
        return Err(refuse(
            "strip-count-exceeded",
            format!("{strip_count} strips exceed {MAX_STRIPS}"),
            "aggregate panels into strips",
        ));
    }
    if !rho_kg_m3.is_finite()
        || rho_kg_m3 <= 0.0
        || !nu.iter().all(|v| v.is_finite())
        || control_rates.len() != nc
        || !control_rates.iter().all(|v| v.is_finite())
    {
        return Err(refuse(
            "non-finite-input",
            "rho/nu/control_rates".into(),
            "rho > 0, finite state, control_rates.len() == nc",
        ));
    }
    let mut m_rr = [[0.0f64; 6]; 6];
    let mut m_rc: Vec<[f64; 6]> = vec![[0.0; 6]; nc];
    let mut m_cc: Vec<Vec<f64>> = vec![vec![0.0; nc]; nc];
    let mut q_added_bias = vec![0.0f64; 6 + nc];
    for s in strips {
        let s = s?;
        for (index, mode) in s.modes.iter().enumerate() {
            // Keep the fixed-strip accumulation order: translation first,
            // followed by each rotary mode's angular 3x3 block.
            let first = if index == 0 { 0 } else { 3 };
            for a in first..6 {
                for b in first..6 {
                    m_rr[a][b] += mode.inertia * (mode.rigid[a] * mode.rigid[b]);
                }
            }
            if let Some(k) = s.control_coord {
                for a in first..6 {
                    m_rc[k][a] += mode.inertia * mode.rigid[a] * mode.control_gain;
                }
                m_cc[k][k] += mode.inertia * mode.control_gain * mode.control_gain;
                if let Some(derivative) = &s.derivative {
                    // T = a/2 (j(q).nu + g*qdot)^2. A rigid hinge has
                    // constant g. Apply -Mdot*u plus dT/dq to the control
                    // row directly, without allocating derivative matrices.
                    let d = derivative[index];
                    let body_speed: f64 = mode.rigid.iter().zip(nu).map(|(j, v)| j * v).sum();
                    let derivative_speed: f64 = d.iter().zip(nu).map(|(j, v)| j * v).sum();
                    let speed = body_speed + mode.control_gain * control_rates[k];
                    for a in first..6 {
                        q_added_bias[a] -= mode.inertia
                            * control_rates[k]
                            * (d[a] * speed + mode.rigid[a] * derivative_speed);
                    }
                    q_added_bias[6 + k] += mode.inertia * body_speed * derivative_speed;
                }
            }
        }
    }
    // Euler–Poincaré bias on the ADDED momentum p = M_rr·ν + M_rc·q̇:
    // bias_lin = −ω × p_lin, bias_ang = −(ω × p_ang + v × p_lin). Workless.
    let mut p = [0.0f64; 6];
    for a in 0..6 {
        for b in 0..6 {
            p[a] += m_rr[a][b] * nu[b];
        }
        for (k, col) in m_rc.iter().enumerate() {
            p[a] += col[a] * control_rates[k];
        }
    }
    let v = [nu[0], nu[1], nu[2]];
    let w = [nu[3], nu[4], nu[5]];
    let p_lin = [p[0], p[1], p[2]];
    let p_ang = [p[3], p[4], p[5]];
    let bias_lin = cross(w, p_lin);
    let wxpa = cross(w, p_ang);
    let vxpl = cross(v, p_lin);
    for i in 0..3 {
        // Retain the fixed-strip path's signed zeros and determinism digest.
        q_added_bias[i] = if q_added_bias[i] == 0.0 {
            -bias_lin[i]
        } else {
            q_added_bias[i] - bias_lin[i]
        };
        q_added_bias[3 + i] = if q_added_bias[3 + i] == 0.0 {
            -(wxpa[i] + vxpl[i])
        } else {
            q_added_bias[3 + i] - (wxpa[i] + vxpl[i])
        };
    }
    // For fixed strips there are no configuration terms and control bias is
    // exactly zero. Moving-strip terms have already been accumulated above.
    if !m_rr
        .iter()
        .flatten()
        .chain(m_rc.iter().flatten())
        .chain(m_cc.iter().flatten())
        .chain(&q_added_bias)
        .all(|x| x.is_finite())
    {
        return Err(refuse(
            "strip-invalid",
            "non-finite assembled fluid inertia or momentum bias".into(),
            "use representable geometry, density, velocities and control gains",
        ));
    }
    Ok(AeroGeneralizedLoads {
        q_rigid_nonaccel: [0.0; 6],
        q_control_nonaccel: vec![0.0; nc],
        m_added_rr: m_rr,
        m_added_rc: m_rc,
        m_added_cc: m_cc,
        q_added_bias,
    })
}

/// Solve the JOINT (6+nc) system with the declared deterministic dense
/// Cholesky (fixed elimination order; the total effective mass is PD).
/// Returns (ν̇[6], q̈_ctrl[nc]).
///
/// This is an algebraic solve, not a free-body time integrator. It has no
/// velocity input: `q_rigid_nonaccel` must already contain the rigid-body
/// inertial right-hand side, including Euler/Coriolis bias when the result
/// is interpreted as the derivative of body-coordinate velocity. Only the
/// separately assembled fluid bias is added by this function.
///
/// # Errors
/// `control-inertia-invalid` (m_control len/values), `joint-mass-not-pd`
/// (Cholesky pivot failure — the total effective mass must be PD).
pub fn solve_joint(
    body: &RigidBody,
    m_control: &[f64],
    loads: &AeroGeneralizedLoads,
) -> Result<(Vec<f64>, Vec<f64>), Refusal> {
    body.admit()?;
    joint::admit_fluid(loads)?;
    let nc = loads.q_control_nonaccel.len();
    if m_control.len() != nc || !m_control.iter().all(|m| m.is_finite() && *m > 0.0) {
        return Err(refuse(
            "control-inertia-invalid",
            format!("m_control len {} vs nc {nc}", m_control.len()),
            "one positive generalized inertia per active coordinate",
        ));
    }
    let n = 6 + nc;
    let mut a = vec![0.0f64; n * n];
    // Rigid block: mass at the COM and the complete authored-frame inertia.
    // Off-diagonal entries must survive the same joint solve as fluid inertia.
    let inertia = body.inertia_tensor().matrix();
    for i in 0..6 {
        for j in 0..6 {
            a[i * n + j] = loads.m_added_rr[i][j];
            if i >= 3 && j >= 3 {
                a[i * n + j] += inertia.m[3 * (i - 3) + j - 3];
            }
        }
        if i < 3 {
            a[i * n + i] += body.mass_kg;
        }
    }
    for (k, col) in loads.m_added_rc.iter().enumerate() {
        for i in 0..6 {
            a[i * n + (6 + k)] = col[i];
            a[(6 + k) * n + i] = col[i];
        }
    }
    for k in 0..nc {
        for l in 0..nc {
            a[(6 + k) * n + (6 + l)] = loads.m_added_cc[k][l];
        }
        a[(6 + k) * n + (6 + k)] += m_control[k];
    }
    let mut rhs = vec![0.0f64; n];
    for i in 0..6 {
        rhs[i] = loads.q_rigid_nonaccel[i] + loads.q_added_bias[i];
    }
    for k in 0..nc {
        rhs[6 + k] = loads.q_control_nonaccel[k] + loads.q_added_bias[6 + k];
    }
    let rhs = joint::solve_positive(a, rhs, n, 0.0)?;
    Ok((rhs[..6].to_vec(), rhs[6..].to_vec()))
}
