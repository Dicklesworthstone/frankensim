//! Analytic tube constructors: constant-twist screws and the Wankel
//! rotor pose.
//!
//! Blade placement and signs for the rotational and ideal (
//! translational) bivectors are extracted from fs-ga at runtime by
//! probing `Motor::rotor` / `Motor::translator` on the coordinate
//! axes: only the SIGN of a component is read (exact), never a rounded
//! magnitude, so no transcription or rounding enters a polynomial
//! coefficient. Trigonometric models use `cos u = 1 − 2·sin²(u/2)`;
//! no irrational constant is baked into a coefficient.
//!
//! Constructor semantics (contract): the tube encloses its CONSTRUCTED
//! component path. Deviation of that path from the ideal real-number
//! motion (from non-unit axes or constant-pose products rounded in
//! f64) is measured by the reported versor defect, not separately
//! certified.

use crate::MotionError;
use crate::algebra::{SCALAR, TmMv};
use crate::tube::{CertifiedMotorTube, MotorTubeSegment};
use fs_ga::{Motor, Pga};
use fs_ivl::{Interval, TaylorModel1};
use std::sync::OnceLock;

/// Blade index and sign for one coordinate axis of a bivector family.
#[derive(Debug, Clone, Copy)]
struct AxisSlot {
    blade: usize,
    sign: f64,
}

/// Extracted placement of the rotational (Euclidean) and ideal
/// bivector coordinate axes.
#[derive(Debug)]
struct AxisFrame {
    rot: [AxisSlot; 3],
    ideal: [AxisSlot; 3],
}

fn single_grade2_slot(m: &Motor, expected_negative_of: f64) -> AxisSlot {
    let mut found: Option<AxisSlot> = None;
    for (k, &c) in m.0.0.iter().enumerate() {
        if k == SCALAR || c == 0.0 {
            continue;
        }
        assert!(
            found.is_none(),
            "axis probe produced more than one non-scalar component"
        );
        // The probed coefficient is `expected_negative_of · (−sign)`,
        // so only the sign of `c` is consumed — exact.
        let sign = if (c < 0.0) == (expected_negative_of > 0.0) {
            1.0
        } else {
            -1.0
        };
        found = Some(AxisSlot { blade: k, sign });
    }
    found.expect("axis probe produced no non-scalar component")
}

fn extract_axis_frame() -> AxisFrame {
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut rot = [AxisSlot { blade: 0, sign: 1.0 }; 3];
    let mut ideal = [AxisSlot { blade: 0, sign: 1.0 }; 3];
    for (i, axis) in axes.iter().enumerate() {
        // rotor(a, φ) = cos(φ/2) − sin(φ/2)·B_a : coefficient at the
        // B_a blade is −sin(0.5)·sign for φ = 1.
        rot[i] = single_grade2_slot(&Motor::rotor(*axis, 1.0), fs_math::det::sin(0.5));
        // translator(a) = 1 − 0.5·N_a* : coefficient is −0.5·sign.
        ideal[i] = single_grade2_slot(&Motor::translator(axis[0], axis[1], axis[2]), 0.5);
    }
    AxisFrame { rot, ideal }
}

fn axis_frame() -> &'static AxisFrame {
    static FRAME: OnceLock<AxisFrame> = OnceLock::new();
    FRAME.get_or_init(extract_axis_frame)
}

/// The Euclidean bivector `B_a` with `rotor(a, φ) =
/// cos(φ/2) − sin(φ/2)·B_a`.
fn rotational_bivector(axis: [f64; 3]) -> Pga {
    let frame = axis_frame();
    let mut b = Pga::zero();
    for (i, &a) in axis.iter().enumerate() {
        b.0[frame.rot[i].blade] += a * frame.rot[i].sign;
    }
    b
}

/// The ideal bivector `N_a*` with `translator(a·d) = 1 − (d/2)·N_a*`.
fn ideal_bivector(axis: [f64; 3]) -> Pga {
    let frame = axis_frame();
    let mut b = Pga::zero();
    for (i, &a) in axis.iter().enumerate() {
        b.0[frame.ideal[i].blade] += a * frame.ideal[i].sign;
    }
    b
}

fn validate_domain(domain: Interval) -> Result<(), MotionError> {
    if !(domain.lo().is_finite() && domain.hi().is_finite()) || domain.lo() >= domain.hi() {
        return Err(MotionError::EmptyTimeDomain);
    }
    Ok(())
}

fn validate_segments(segments: usize) -> Result<(), MotionError> {
    if segments == 0 {
        return Err(MotionError::InvalidSegments);
    }
    Ok(())
}

fn segment_domains(domain: Interval, segments: usize) -> Vec<Interval> {
    let (lo, hi) = (domain.lo(), domain.hi());
    let n = segments as f64;
    let mut cuts = Vec::with_capacity(segments + 1);
    cuts.push(lo);
    for i in 1..segments {
        cuts.push(lo + (hi - lo) * (i as f64 / n));
    }
    cuts.push(hi);
    cuts.windows(2).map(|w| Interval::new(w[0], w[1])).collect()
}

/// `cos` of a model via `1 − 2·sin²(u/2)` (no irrational coefficient).
fn cos_model(u: &TaylorModel1) -> Result<TaylorModel1, MotionError> {
    let half = u.scale(0.5)?;
    let s = half.sin()?;
    let s2 = s.try_mul(&s)?;
    let one = TaylorModel1::constant(1.0, u.domain(), u.order())?;
    Ok(one.try_sub(&s2.scale(2.0)?)?)
}

/// Accumulate `coeffs · basis(t)` into the multivector's components.
fn accumulate(
    out: &mut TmMv,
    coeffs: &Pga,
    basis: &TaylorModel1,
) -> Result<(), MotionError> {
    for (k, &c) in coeffs.0.iter().enumerate() {
        if c == 0.0 {
            continue;
        }
        let term = basis.scale(c)?;
        let combined = if out.is_nonzero(k) {
            out.component(k).try_add(&term)?
        } else {
            term
        };
        out.set(k, combined)?;
    }
    Ok(())
}

/// Parameters of a constant-twist screw motion: rotation at `omega`
/// (rad/s) about the axis line through `center`, plus translation at
/// `axial_velocity` (m/s) along the axis, all pre-composed with
/// `base_pose` (applied on the left, world side).
#[derive(Debug, Clone, Copy)]
pub struct ScrewParams {
    /// Rotation/translation axis (should be unit; deviation shows up
    /// in the tube's versor defect).
    pub axis: [f64; 3],
    /// A point on the axis.
    pub center: [f64; 3],
    /// Angular rate (rad/s).
    pub omega: f64,
    /// Translation rate along the axis (m/s).
    pub axial_velocity: f64,
    /// World-side constant pose.
    pub base_pose: Motor,
}

/// Build a certified tube for a constant-twist screw over `domain`,
/// split into `segments` pieces of Taylor order `order`.
pub fn screw_tube(
    params: &ScrewParams,
    domain: Interval,
    order: usize,
    segments: usize,
) -> Result<CertifiedMotorTube, MotionError> {
    validate_domain(domain)?;
    validate_segments(segments)?;
    for v in params
        .axis
        .iter()
        .chain(params.center.iter())
        .chain([params.omega, params.axial_velocity].iter())
    {
        if !v.is_finite() {
            return Err(MotionError::NonFiniteInput { what: "screw parameter" });
        }
    }
    let b_axis = rotational_bivector(params.axis);
    let n_axis = ideal_bivector(params.axis);
    // (cosθ − sinθ·B)(1 − (vt/2)·N*) expanded over the four basis
    // functions {cosθ, sinθ, t·cosθ, t·sinθ}:
    let u1 = Pga::scalar(1.0);
    let u2 = b_axis.scale(-1.0);
    let u3 = n_axis.scale(-params.axial_velocity / 2.0);
    let u4 = b_axis.gp(&n_axis).scale(params.axial_velocity / 2.0);
    // World-side and body-side constant motors: A = base ∘ T_center,
    // Bc = T_center⁻¹.
    let t_center = Motor::translator(params.center[0], params.center[1], params.center[2]);
    let t_back = Motor::translator(-params.center[0], -params.center[1], -params.center[2]);
    let a = params.base_pose.compose(&t_center);
    let c1 = a.0.gp(&u1).gp(&t_back.0);
    let c2 = a.0.gp(&u2).gp(&t_back.0);
    let c3 = a.0.gp(&u3).gp(&t_back.0);
    let c4 = a.0.gp(&u4).gp(&t_back.0);

    let mut sealed = Vec::with_capacity(segments);
    for sub in segment_domains(domain, segments) {
        let t = TaylorModel1::variable(sub, order)?;
        let theta = t.scale(params.omega)?.scale(0.5)?;
        let cos_t = cos_model(&theta)?;
        let sin_t = theta.sin()?;
        let t_cos = t.try_mul(&cos_t)?;
        let t_sin = t.try_mul(&sin_t)?;
        let mut mv = TmMv::zero(sub, order)?;
        accumulate(&mut mv, &c1, &cos_t)?;
        accumulate(&mut mv, &c2, &sin_t)?;
        accumulate(&mut mv, &c3, &t_cos)?;
        accumulate(&mut mv, &c4, &t_sin)?;
        sealed.push(MotorTubeSegment::seal(mv)?);
    }
    CertifiedMotorTube::from_segments(sealed)
}

/// Parameters of the Wankel rotor POSE: the rotor center orbits the
/// crank axis (z through the origin) at radius `eccentricity` and
/// crank angle `α(t) = omega·t + crank_phase`, while the rotor spins
/// at one third of the crank rate, `β(t) = α(t)/3 + rotor_phase`. The
/// epitrochoid housing curve is the derived APEX LOCUS — deliberately
/// not constructed here.
#[derive(Debug, Clone, Copy)]
pub struct WankelParams {
    /// Crank eccentricity (m).
    pub eccentricity: f64,
    /// Crank angular rate (rad/s).
    pub omega: f64,
    /// Crank phase at `t = 0` (rad).
    pub crank_phase: f64,
    /// Rotor phase offset at `t = 0` (rad).
    pub rotor_phase: f64,
    /// World-side constant pose.
    pub base_pose: Motor,
}

/// Build a certified tube for the Wankel rotor pose over `domain`.
pub fn wankel_tube(
    params: &WankelParams,
    domain: Interval,
    order: usize,
    segments: usize,
) -> Result<CertifiedMotorTube, MotionError> {
    validate_domain(domain)?;
    validate_segments(segments)?;
    for v in [
        params.eccentricity,
        params.omega,
        params.crank_phase,
        params.rotor_phase,
    ] {
        if !v.is_finite() {
            return Err(MotionError::NonFiniteInput { what: "wankel parameter" });
        }
    }
    let frame_x = ideal_bivector([1.0, 0.0, 0.0]);
    let frame_y = ideal_bivector([0.0, 1.0, 0.0]);
    let b_z = rotational_bivector([0.0, 0.0, 1.0]);
    let mut sealed = Vec::with_capacity(segments);
    for sub in segment_domains(domain, segments) {
        let t = TaylorModel1::variable(sub, order)?;
        // α(t) = ω t + α0 ; β/2 = α/6 + β0/2.
        let alpha = t.scale(params.omega)?
            .try_add(&TaylorModel1::constant(params.crank_phase, sub, order)?)?;
        let cos_a = cos_model(&alpha)?;
        let sin_a = alpha.sin()?;
        let beta_half = alpha.scale(1.0 / 6.0)?.try_add(&TaylorModel1::constant(
            params.rotor_phase * 0.5,
            sub,
            order,
        )?)?;
        let cos_bh = cos_model(&beta_half)?;
        let sin_bh = beta_half.sin()?;
        // T(orbit) = 1 − (e/2)(cosα·N_x* + sinα·N_y*)
        let mut t_mv = TmMv::zero(sub, order)?;
        accumulate(
            &mut t_mv,
            &Pga::scalar(1.0),
            &TaylorModel1::constant(1.0, sub, order)?,
        )?;
        accumulate(
            &mut t_mv,
            &frame_x.scale(-params.eccentricity / 2.0),
            &cos_a,
        )?;
        accumulate(
            &mut t_mv,
            &frame_y.scale(-params.eccentricity / 2.0),
            &sin_a,
        )?;
        // R(β) = cos(β/2) − sin(β/2)·B_z
        let mut r_mv = TmMv::zero(sub, order)?;
        accumulate(&mut r_mv, &Pga::scalar(1.0), &cos_bh)?;
        accumulate(&mut r_mv, &b_z.scale(-1.0), &sin_bh)?;
        let pose = t_mv.gp(&r_mv)?;
        // World-side constant pose, rebuilt on this segment's domain.
        let base_seg = TmMv::constant(&params.base_pose.0, sub, order)?;
        let mv = base_seg.gp(&pose)?;
        sealed.push(MotorTubeSegment::seal(mv)?);
    }
    CertifiedMotorTube::from_segments(sealed)
}
