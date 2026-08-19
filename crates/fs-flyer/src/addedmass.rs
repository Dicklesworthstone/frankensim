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
//! strip's normal-velocity Jacobian j = [n; r×n] (plus the control gain g
//! on its dynamic coordinate) contributes m_a·[j; g·e_k][j; g·e_k]ᵀ — so
//! the assembled M_added is symmetric PSD BY CONSTRUCTION and the
//! rigid↔control cross blocks arise from the same outer product (no
//! separate, forgettable code path).
//!
//! `Q_added_bias` is the Euler–Poincaré gyroscopic term on the added
//! momentum, −ad*_ν(p): WORKLESS for any p (νᵀ·ad*_ν(p) ≡ 0), which the
//! battery checks to machine precision — energy consistency without any
//! runtime finite-differencing. Configuration-independent strips make the
//! control-block bias identically zero (disclosed below).
//!
//! No-claims: cross-surface fluid-inertia interference is OMITTED
//! (AnalyticStrip's disclosed boundary; PanelExtracted is the additive
//! referee mode, E4 lane). Propeller disks are EXCLUDED — rotor terms
//! stay in fs-airscrew.

use crate::Refusal;
use crate::spine::RigidBody;

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
    /// Dynamic control coordinate driving this strip, if any.
    pub control_coord: Option<usize>,
    /// Normal-velocity gain g [m/rad·s⁻¹ per rad/s] of that coordinate
    /// (hinge arm × mode shape at the centroid).
    pub control_gain: f64,
}

/// The assembled generalized loads (plan §5.1.2 shape).
#[derive(Clone, Debug, PartialEq)]
pub struct AeroGeneralizedLoads {
    /// Non-acceleration rigid loads [6] (force; moment), body frame.
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

/// Density of the working fluid [kg/m³] enters per assembly call (the
/// provenance-bound air state supplies it).
///
/// Assemble the added-mass blocks + bias for the strip set at the current
/// generalized velocity ν = (v_body[3], ω_body[3]) and control rates.
///
/// # Errors
/// `control-coords-outside-domain` (nc cap at cap AND cap+1),
/// `strip-count-exceeded`, `strip-invalid` (non-finite, non-positive
/// chord/span, non-unit normal, control coordinate out of range),
/// `non-finite-input`.
pub fn assemble_analytic_strip(
    rho_kg_m3: f64,
    strips: &[Strip],
    nc: usize,
    nu: &[f64; 6],
    control_rates: &[f64],
) -> Result<AeroGeneralizedLoads, Refusal> {
    if nc > MAX_CONTROL_COORDS {
        return Err(refuse(
            "control-coords-outside-domain",
            format!("nc = {nc} exceeds the cap {MAX_CONTROL_COORDS}"),
            "the 1903 Flyer has at most a few active dynamic coordinates",
        ));
    }
    if strips.len() > MAX_STRIPS {
        return Err(refuse(
            "strip-count-exceeded",
            format!("{} strips exceed {MAX_STRIPS}", strips.len()),
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
    for s in strips {
        let finite = s.chord_m.is_finite()
            && s.span_m.is_finite()
            && s.position_m.iter().all(|v| v.is_finite())
            && s.normal.iter().all(|v| v.is_finite())
            && s.control_gain.is_finite();
        let norm2: f64 = s.normal.iter().map(|v| v * v).sum();
        if !finite || s.chord_m <= 0.0 || s.span_m <= 0.0 || (norm2 - 1.0).abs() > 1.0e-9 {
            return Err(refuse(
                "strip-invalid",
                format!("strip {}: chord/span/normal invalid", s.name),
                "positive chord/span; unit normal",
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
        // 2-D flat-plate added mass per unit span, times strip span.
        let m_a = rho_kg_m3 * core::f64::consts::PI * s.chord_m * s.chord_m / 4.0 * s.span_m;
        let rxn = cross(s.position_m, s.normal);
        let j = [
            s.normal[0],
            s.normal[1],
            s.normal[2],
            rxn[0],
            rxn[1],
            rxn[2],
        ];
        for a in 0..6 {
            for b in 0..6 {
                m_rr[a][b] += m_a * (j[a] * j[b]);
            }
        }
        if let Some(k) = s.control_coord {
            for a in 0..6 {
                m_rc[k][a] += m_a * j[a] * s.control_gain;
            }
            m_cc[k][k] += m_a * s.control_gain * s.control_gain;
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
    let mut q_added_bias = vec![0.0f64; 6 + nc];
    for i in 0..3 {
        q_added_bias[i] = -bias_lin[i];
        q_added_bias[3 + i] = -(wxpa[i] + vxpl[i]);
    }
    // Control-block bias: identically zero for configuration-independent
    // strips (disclosed in the module no-claims).
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
/// # Errors
/// `control-inertia-invalid` (m_control len/values), `joint-mass-not-pd`
/// (Cholesky pivot failure — the total effective mass must be PD).
pub fn solve_joint(
    body: &RigidBody,
    m_control: &[f64],
    loads: &AeroGeneralizedLoads,
) -> Result<(Vec<f64>, Vec<f64>), Refusal> {
    body.admit()?;
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
    // Rigid block: diag(m, m, m, Ixx, Iyy, Izz) + M_rr.
    let diag = [
        body.mass_kg,
        body.mass_kg,
        body.mass_kg,
        body.inertia_kgm2[0],
        body.inertia_kgm2[1],
        body.inertia_kgm2[2],
    ];
    for i in 0..6 {
        for j in 0..6 {
            a[i * n + j] = loads.m_added_rr[i][j];
        }
        a[i * n + i] += diag[i];
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
    // Fixed-order dense Cholesky (the declared deterministic factorization).
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return Err(refuse(
                        "joint-mass-not-pd",
                        format!("Cholesky pivot {sum:e} at row {i}"),
                        "the total effective mass must be PD; check strip geometry",
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
    Ok((rhs[..6].to_vec(), rhs[6..].to_vec()))
}
