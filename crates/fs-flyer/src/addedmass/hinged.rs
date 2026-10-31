//! Finite rigid-hinge motion for the existing analytic-strip inertia law.
//!
//! With generalized velocity `u = (nu, qdot)`, fluid storage is
//! `T = 1/2 u^T M(q) u`. The body equations require the original coadjoint
//! bias AND `-Mdot*u`; control row k additionally requires `dT/dq_k`.
//! Dropping those configuration terms creates or destroys kinetic energy as
//! a control moves, even when the instantaneous mass matrix is correct.
//!
//! Each strip rotates around one hinge fixed in the carrier body. Several
//! strips may share a coordinate. Its three positive kinetic modes and their
//! analytic derivatives feed the same assembly and joint solve as fixed
//! strips. Hinges in a serial moving chain, strip deformation, changing fluid
//! density, and cross-surface hydrodynamic interaction are outside this model.

use super::{
    AeroGeneralizedLoads, MAX_CONTROL_COORDS, Strip, StripModes, admit_strip, assemble_modes,
    cross, refuse, strip_modes,
};
use crate::Refusal;
use fs_math::det;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// A finite plate with a rigid hinge fixed in carrier-body coordinates.
/// Geometry is immutable after admission. All velocities are relative to the
/// same body reference as the fixed strips and the joint rigid inertia.
#[derive(Clone, Debug, PartialEq)]
pub struct HingedStrip {
    reference: Strip,
    reference_arm_m: [f64; 3],
    hinge_axis: [f64; 3],
    control_coord: usize,
    control_gains: [f64; 3],
}

impl HingedStrip {
    /// Prepare a strip at zero control angle and its actual hinge geometry.
    ///
    /// `reference` supplies geometry only: its control coordinate must be
    /// `None`, with zero linear/angular control gains. The prepared strip
    /// derives those gains from the hinge; no authored gain is discarded.
    /// `hinge_axis` is an oriented unit vector. Positive angles and rates use
    /// its right-hand sense. A center hinge has zero normal-translation gain
    /// but still retains the plate's pitch and spanwise rotary fluid inertia.
    ///
    /// # Errors
    /// `hinge-invalid` for invalid hinge geometry or preassigned control gains;
    /// `control-coords-outside-domain` for a coordinate at or beyond the cap;
    /// the original strip-geometry refusals also apply.
    pub fn new(
        mut reference: Strip,
        hinge_point_m: [f64; 3],
        hinge_axis: [f64; 3],
        control_coord: usize,
    ) -> Result<Self, Refusal> {
        if control_coord >= MAX_CONTROL_COORDS {
            return Err(refuse(
                "control-coords-outside-domain",
                format!("hinge coordinate {control_coord} outside 0..{MAX_CONTROL_COORDS}"),
                "use an admitted active dynamic coordinate",
            ));
        }
        if reference.control_coord.is_some()
            || reference.control_gain != 0.0
            || reference.control_angular_gain.iter().any(|g| *g != 0.0)
        {
            return Err(refuse(
                "hinge-invalid",
                format!(
                    "strip {} already declares control kinematics",
                    reference.name
                ),
                "provide geometry with inactive control fields; the hinge derives its own gains",
            ));
        }
        admit_strip(&reference, 0)?;
        let norm2 = dot(hinge_axis, hinge_axis);
        if !hinge_point_m
            .iter()
            .chain(&hinge_axis)
            .all(|x| x.is_finite())
            || !norm2.is_finite()
            || (norm2 - 1.0).abs() > 1e-9
        {
            return Err(refuse(
                "hinge-invalid",
                "hinge needs a finite point and a unit axis".into(),
                "supply the physical oriented hinge line in body coordinates",
            ));
        }
        // Remove only admitted roundoff in the nominal unit vector, once.
        let hinge_axis = hinge_axis.map(|x| x / norm2.sqrt());
        let reference_arm_m = core::array::from_fn(|i| reference.position_m[i] - hinge_point_m[i]);
        let normal_gain = dot(reference.normal, cross(hinge_axis, reference_arm_m));
        let control_gains = [
            normal_gain,
            dot(reference.span_axis, hinge_axis),
            dot(cross(reference.span_axis, reference.normal), hinge_axis),
        ];
        if !reference_arm_m
            .iter()
            .chain(&control_gains)
            .all(|x| x.is_finite())
        {
            return Err(refuse(
                "hinge-invalid",
                "unrepresentable hinge arm or control velocity".into(),
                "use representable geometry relative to the body reference",
            ));
        }
        reference.control_coord = Some(control_coord);
        reference.control_gain = normal_gain;
        reference.control_angular_gain = hinge_axis;
        Ok(Self {
            reference,
            reference_arm_m,
            hinge_axis,
            control_coord,
            control_gains,
        })
    }

    /// Active generalized coordinate driven by this hinge.
    #[must_use]
    pub const fn control_coord(&self) -> usize {
        self.control_coord
    }

    fn evaluate(&self, angle_rad: f64) -> Result<(Strip, [f64; 3]), Refusal> {
        if !angle_rad.is_finite() {
            return Err(refuse(
                "non-finite-input",
                "nonfinite hinge angle".into(),
                "supply a finite control state",
            ));
        }
        let mut strip = self.reference.clone();
        let mut arm = self.reference_arm_m;
        if angle_rad != 0.0 {
            // Half-angle Rodrigues coefficients retain the q^2 term near
            // zero. Applying the displacement to r0 avoids pivot + arm
            // cancellation for a small rotation about a distant hinge.
            let (s, c) = (det::sin(0.5 * angle_rad), det::cos(0.5 * angle_rad));
            let sine = 2.0 * s * c;
            let one_minus_cosine = 2.0 * s * s;
            let displacement = |v| {
                let first = cross(self.hinge_axis, v);
                let second = cross(self.hinge_axis, first);
                core::array::from_fn::<_, 3, _>(|i| sine * first[i] + one_minus_cosine * second[i])
            };
            let arm_change = displacement(arm);
            let normal_change = displacement(strip.normal);
            let span_change = displacement(strip.span_axis);
            for i in 0..3 {
                arm[i] += arm_change[i];
                strip.position_m[i] += arm_change[i];
                strip.normal[i] += normal_change[i];
                strip.span_axis[i] += span_change[i];
            }
        }
        admit_strip(&strip, self.control_coord + 1)?;
        let position_derivative = cross(self.hinge_axis, arm);
        if !position_derivative.iter().all(|x| x.is_finite()) {
            return Err(refuse(
                "hinge-invalid",
                "unrepresentable moving hinge velocity".into(),
                "use representable geometry and finite control angles",
            ));
        }
        Ok((strip, position_derivative))
    }

    /// Evaluate actual strip geometry and its derived control Jacobian.
    ///
    /// This is useful for coincident aerodynamic load evaluation. Calling the
    /// fixed-strip assembler on this geometry alone omits its configuration
    /// bias; dynamic callers use [`assemble_hinged_strips`].
    ///
    /// # Errors
    /// Nonfinite angle or unrepresentable transformed geometry.
    pub fn at_angle(&self, angle_rad: f64) -> Result<Strip, Refusal> {
        self.evaluate(angle_rad).map(|(strip, _)| strip)
    }

    fn modes(&self, rho: f64, angle_rad: f64, nc: usize) -> Result<StripModes, Refusal> {
        let (strip, dr) = self.evaluate(angle_rad)?;
        let mut modes = strip_modes(rho, &strip, nc)?;
        let dn = cross(self.hinge_axis, strip.normal);
        let dt = cross(self.hinge_axis, strip.span_axis);
        let first = cross(dr, strip.normal);
        let second = cross(strip.position_m, dn);
        let ds = cross(self.hinge_axis, cross(strip.span_axis, strip.normal));
        modes.derivative = Some([
            [
                dn[0],
                dn[1],
                dn[2],
                first[0] + second[0],
                first[1] + second[1],
                first[2] + second[2],
            ],
            [0.0, 0.0, 0.0, dt[0], dt[1], dt[2]],
            [0.0, 0.0, 0.0, ds[0], ds[1], ds[2]],
        ]);
        // Rigid rotation around the same axis leaves all three projections
        // invariant. Retain their reference values without roundoff drift.
        for (mode, gain) in modes.modes.iter_mut().zip(self.control_gains) {
            mode.control_gain = gain;
        }
        Ok(modes)
    }
}

/// Assemble fixed strips and finite moving hinges in one generalized mass law.
///
/// `control_angles` and `control_rates` are ordered by active coordinate. The
/// returned bias obeys `u.dot(bias) = -1/2 u^T Mdot u`, including nonzero
/// control-row bias; its body coadjoint part remains workless. Fluid density is
/// held fixed in this mechanical derivative. All strips share the original
/// normal-translation and centroidal rotary coefficients and existing caps.
/// No configuration derivative matrices or temporary geometry vectors are built.
///
/// # Errors
/// The original assembly refusals, angle/rate count mismatch, nonfinite angles,
/// or a hinge coordinate missing from the active control state.
pub fn assemble_hinged_strips(
    rho: f64,
    fixed: &[Strip],
    moving: &[HingedStrip],
    control_angles: &[f64],
    nu: &[f64; 6],
    control_rates: &[f64],
) -> Result<AeroGeneralizedLoads, Refusal> {
    let nc = control_angles.len();
    if control_rates.len() != nc || !control_angles.iter().all(|q| q.is_finite()) {
        return Err(refuse(
            "non-finite-input",
            "hinge angles/rates have inconsistent counts or nonfinite angles".into(),
            "one finite angle and rate per active coordinate",
        ));
    }
    if moving.iter().any(|s| s.control_coord >= nc) {
        return Err(refuse(
            "control-coords-outside-domain",
            "hinge coordinate is missing from the active state".into(),
            "include every hinge's angle and rate in coordinate order",
        ));
    }
    let count = fixed.len().checked_add(moving.len()).ok_or_else(|| {
        refuse(
            "strip-count-exceeded",
            "combined fixed/moving strip count overflow".into(),
            "aggregate panels into strips",
        )
    })?;
    assemble_modes(
        rho,
        count,
        nc,
        nu,
        control_rates,
        fixed.iter().map(|s| strip_modes(rho, s, nc)).chain(
            moving
                .iter()
                .map(|s| s.modes(rho, control_angles[s.control_coord], nc)),
        ),
    )
}
