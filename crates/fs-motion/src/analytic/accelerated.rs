//! Prescribed accelerating motion with rigorously paired time derivatives.

use super::{
    MotorRateTube, accumulate, chained_flip, cos_model, ideal_bivector, rotational_bivector,
    segment_domains, validate_domain, validate_segments,
};
use crate::MotionError;
use crate::algebra::{
    BLADES, POINT_BLADE_NX, POINT_BLADE_NZ, POINT_BLADE_W, POINT_BLADE_Y, TmMv, homogeneous_point,
    point_to_mv,
};
use crate::tube::{CertifiedMotorTube, EnclosureClass, MotorTubeSegment};
use fs_exec::Cx;
use fs_ga::{Motor, Pga};
use fs_geom::Point3;
use fs_ivl::{Interval, TaylorModel1};

/// A prescribed fixed-axis spin ramp and uniformly accelerating translation.
/// All vectors use the fixed reference frame before `base_pose` is applied;
/// translation is independent of the rotating frame. At `reference_time_s`
/// the pose equals `base_pose`.
#[derive(Debug, Clone, Copy)]
pub struct AcceleratedRigidParams {
    /// Unit rotation axis; squared-norm roundoff up to `64 * EPSILON` is admitted.
    pub axis: [f64; 3],
    /// Rotation center in the reference frame (m).
    pub center: [f64; 3],
    /// Angular velocity at the reference time (rad/s), along `axis`.
    pub angular_velocity: f64,
    /// Constant angular acceleration along `axis` (rad/s²).
    pub angular_acceleration: f64,
    /// Translational velocity at the reference time (m/s).
    pub linear_velocity: [f64; 3],
    /// Constant reference-frame translational acceleration (m/s²).
    pub linear_acceleration: [f64; 3],
    /// Absolute time at which the initial pose and velocities apply (s).
    pub reference_time_s: f64,
    /// Constant world-side pose, applied on the left.
    pub base_pose: Motor,
}

/// World-coordinate motion of a body-fixed point, enclosed over a time span.
#[derive(Debug, Clone, Copy)]
pub struct PointMotionEnclosure {
    /// Position (m).
    pub position: [Interval; 3],
    /// Inertial velocity (m/s).
    pub velocity: [Interval; 3],
    /// Kinematic acceleration (m/s²), including tangential and centripetal terms.
    /// This is not accelerometer specific force; no gravity is subtracted.
    pub acceleration: [Interval; 3],
    /// Enclosure strength for the constructed component path.
    pub class: EnclosureClass,
    /// Maximum versor defect of the queried segments.
    pub defect: f64,
}

/// A motor tube bound to its first and second component derivatives.
/// It plugs into native event guards and swept geometry without time sampling.
#[derive(Debug, Clone)]
pub struct AcceleratedRigidTube {
    motion: CertifiedMotorTube,
    rate: MotorRateTube,
    acceleration: Vec<TmMv>,
}

impl AcceleratedRigidTube {
    /// The certified pose path.
    #[must_use]
    pub fn motion(&self) -> &CertifiedMotorTube {
        &self.motion
    }

    /// First component derivative, sign-locked to [`Self::motion`].
    #[must_use]
    pub fn rate(&self) -> &MotorRateTube {
        &self.rate
    }

    /// Transfer the pose path into a spacetime or swept chart.
    #[must_use]
    pub fn into_motion(self) -> CertifiedMotorTube {
        self.motion
    }

    /// Enclose a body-fixed point's world motion for every absolute time in
    /// `span`. Derivatives belong to the same constructed path as the pose;
    /// a homogeneous weight containing zero refuses instead of dividing through.
    pub fn point_motion_over(
        &self,
        point: Point3,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<PointMotionEnclosure, MotionError> {
        if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
            return Err(MotionError::NonFiniteInput { what: "point" });
        }
        let domain = self.motion.domain();
        if !domain.encloses(span) {
            return Err(MotionError::OutOfDomain {
                lo: span.lo(),
                hi: span.hi(),
                domain_lo: domain.lo(),
                domain_hi: domain.hi(),
            });
        }
        let mut result: Option<PointMotionEnclosure> = None;
        for (i, segment) in self.motion.segments().iter().enumerate() {
            let lo = span.lo().max(segment.domain().lo());
            let hi = span.hi().min(segment.domain().hi());
            if lo > hi {
                continue;
            }
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let sub = Interval::new(lo, hi);
            let m = segment.components();
            let d = &self.rate.segments[i];
            let e = &self.acceleration[i];
            let p = TmMv::constant(
                &point_to_mv(point.x, point.y, point.z),
                m.domain(),
                m.order(),
            )?;
            let (mp, dp, ep) = (m.gp(&p)?, d.gp(&p)?, e.gp(&p)?);
            let (rm, rd, re) = (m.reverse()?, d.reverse()?, e.reverse()?);
            let x = mp.gp(&rm)?.eval_all(sub)?;
            let dx = dp.gp(&rm)?.add_componentwise(&mp.gp(&rd)?)?.eval_all(sub)?;
            let middle = dp.gp(&rd)?;
            let ddx = ep
                .gp(&rm)?
                .add_componentwise(&middle)?
                .add_componentwise(&middle)?
                .add_componentwise(&mp.gp(&re)?)?
                .eval_all(sub)?;
            let position = homogeneous_point(&x)?;
            let (w, dw, ddw) = (x[POINT_BLADE_W], dx[POINT_BLADE_W], ddx[POINT_BLADE_W]);
            let (dn, ddn) = (point_numerators(&dx), point_numerators(&ddx));
            let velocity = std::array::from_fn(|k| (dn[k] - position[k] * dw) / w);
            let acceleration = std::array::from_fn(|k| {
                (ddn[k] - position[k] * ddw - Interval::point(2.0) * velocity[k] * dw) / w
            });
            match &mut result {
                Some(prev) => {
                    for k in 0..3 {
                        prev.position[k] = prev.position[k].hull(position[k]);
                        prev.velocity[k] = prev.velocity[k].hull(velocity[k]);
                        prev.acceleration[k] = prev.acceleration[k].hull(acceleration[k]);
                    }
                    prev.defect = prev.defect.max(segment.defect());
                }
                None => {
                    result = Some(PointMotionEnclosure {
                        position,
                        velocity,
                        acceleration,
                        class: EnclosureClass::Certified,
                        defect: segment.defect(),
                    })
                }
            }
        }
        cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
        result.ok_or(MotionError::EmptyTimeDomain)
    }
}

fn point_numerators(x: &[Interval; BLADES]) -> [Interval; 3] {
    [
        Interval::point(0.0) - x[POINT_BLADE_NX],
        x[POINT_BLADE_Y],
        Interval::point(0.0) - x[POINT_BLADE_NZ],
    ]
}

/// Build a rigorous accelerated path and both time derivatives. For
/// `tau = t - reference_time_s`, the path is
/// `base * T(center) * T(v*tau + a*tau²/2) * R(axis, omega*tau + alpha*tau²/2)
/// * T(-center)`. All time-dependent rounding/truncation stays in Taylor
/// remainders. Rounded constant motor products define the constructed path;
/// its versor defect is reported, and no force-integration error is claimed.
pub fn accelerated_rigid_tube(
    params: &AcceleratedRigidParams,
    domain: Interval,
    order: usize,
    segments: usize,
) -> Result<AcceleratedRigidTube, MotionError> {
    validate_domain(domain)?;
    validate_segments(segments)?;
    if params
        .axis
        .iter()
        .chain(params.center.iter())
        .chain(params.linear_velocity.iter())
        .chain(params.linear_acceleration.iter())
        .chain(params.base_pose.0.0.iter())
        .chain(
            [
                params.angular_velocity,
                params.angular_acceleration,
                params.reference_time_s,
            ]
            .iter(),
        )
        .any(|v| !v.is_finite())
    {
        return Err(MotionError::NonFiniteInput {
            what: "accelerated rigid parameter",
        });
    }
    let norm2: f64 = params.axis.iter().map(|v| v * v).sum();
    if (norm2 - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(MotionError::InvalidConfiguration {
            what: "rotation axis must be unit",
        });
    }
    let b = rotational_bivector(params.axis);
    let nv = ideal_bivector(params.linear_velocity);
    let na = ideal_bivector(params.linear_acceleration);
    let tc = Motor::translator(params.center[0], params.center[1], params.center[2]);
    let tb = Motor::translator(-params.center[0], -params.center[1], -params.center[2]);
    let left = params.base_pose.compose(&tc);
    // Translation is on the left of rotation. N*B (not B*N) handles translation
    // perpendicular to the spin axis. Share these exact rounded constants
    // between the primal and both derivatives.
    let coefficients = [
        Pga::scalar(1.0),
        b.scale(-1.0),
        nv.scale(-1.0),
        nv.gp(&b),
        na.scale(-1.0),
        na.gp(&b),
    ]
    .map(|u| left.0.gp(&u).gp(&tb.0));
    let mut sealed = Vec::with_capacity(segments);
    let mut rates = Vec::with_capacity(segments);
    let mut accelerations = Vec::with_capacity(segments);
    for sub in segment_domains(domain, segments) {
        validate_domain(sub)?;
        let constant = |value| TaylorModel1::constant(value, sub, order);
        let tau =
            TaylorModel1::variable(sub, order)?.try_sub(&constant(params.reference_time_s)?)?;
        let tau2 = tau.try_mul(&tau)?;
        let theta = tau
            .scale(params.angular_velocity)?
            .scale(0.5)?
            .try_add(&tau2.scale(params.angular_acceleration)?.scale(0.25)?)?;
        let p = constant(params.angular_velocity)?
            .scale(0.5)?
            .try_add(&tau.scale(params.angular_acceleration)?.scale(0.5)?)?;
        let q = constant(params.angular_acceleration)?.scale(0.5)?;
        let p2 = p.try_mul(&p)?;
        let c = cos_model(&theta)?;
        let s = theta.sin()?;
        let dc = p.try_mul(&s)?.scale(-1.0)?;
        let ds = p.try_mul(&c)?;
        let ddc = q.try_mul(&s)?.try_add(&p2.try_mul(&c)?)?.scale(-1.0)?;
        let dds = q.try_mul(&c)?.try_sub(&p2.try_mul(&s)?)?;
        let mut mv = TmMv::zero(sub, order)?;
        let mut rate = TmMv::zero(sub, order)?;
        let mut acceleration = TmMv::zero(sub, order)?;
        for (j, (f, df, ddf)) in [(&c, &dc, &ddc), (&s, &ds, &dds)].into_iter().enumerate() {
            accumulate(&mut mv, &coefficients[j], f)?;
            accumulate(&mut rate, &coefficients[j], df)?;
            accumulate(&mut acceleration, &coefficients[j], ddf)?;
            // d(tau*f), d²(tau*f), d(tau²*f), d²(tau²*f).
            let tf = tau.try_mul(f)?;
            let dtf = f.try_add(&tau.try_mul(df)?)?;
            let ddtf = df.scale(2.0)?.try_add(&tau.try_mul(ddf)?)?;
            let t2f = tau2.try_mul(f)?;
            let dt2f = tf.scale(2.0)?.try_add(&tau2.try_mul(df)?)?;
            let ddt2f = f
                .scale(2.0)?
                .try_add(&tau.try_mul(df)?.scale(4.0)?)?
                .try_add(&tau2.try_mul(ddf)?)?;
            for (k, factor, f, df, ddf) in
                [(j + 2, 0.5, tf, dtf, ddtf), (j + 4, 0.25, t2f, dt2f, ddt2f)]
            {
                accumulate(&mut mv, &coefficients[k], &f.scale(factor)?)?;
                accumulate(&mut rate, &coefficients[k], &df.scale(factor)?)?;
                accumulate(&mut acceleration, &coefficients[k], &ddf.scale(factor)?)?;
            }
        }
        let flip = chained_flip(&sealed, &mv)?;
        let (segment, _) = MotorTubeSegment::seal_with_sign(mv, flip)?;
        sealed.push(segment);
        rates.push(if flip { rate.negate()? } else { rate });
        accelerations.push(if flip {
            acceleration.negate()?
        } else {
            acceleration
        });
    }
    Ok(AcceleratedRigidTube {
        motion: CertifiedMotorTube::from_segments(sealed)?,
        rate: MotorRateTube { segments: rates },
        acceleration: accelerations,
    })
}
