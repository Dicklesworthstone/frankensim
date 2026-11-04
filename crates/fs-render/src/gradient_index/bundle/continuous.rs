//! Fixed-plane-coordinate variations through a continuous-index gradient jump.
//!
//! Let N be the original plane normal, u the unit tangent, k=u', and tau
//! the launch-angle derivative of arc length at fixed h=N.r. The ordinary
//! fixed-arc columns X,V become Q=X+u*tau and W=V+k*tau, with
//! tau=-N.X/(N.u). Along the central ray, q=N.W/(N.u) gives
//! Q'=W-u*q, W'=A*Q+B*W-k*q, tau'=-q, where A,B are the original
//! smooth variational operator. Q is tangent to each constant-h plane.
//!
//! These variables remain continuous across the original face: the
//! distributional normal Hessian has zero action on Q. Transforming back
//! gives X=Q-u*tau and V=W-k*tau, hence the fixed-arc saltation
//! V+ = V- + (k+ - k-)*(N.X-)/(N.u), with X+ = X-.
//! The native central trajectory actually crosses the face inside its guarded
//! collar; no jump is applied early at the approximate arrival endpoint.
//!
//! The general crossing-time correction is the identity-reset case in Kong
//! et al., *Saltation Matrices* (2023), <https://arxiv.org/abs/2306.06862>.
//! This coordinate form is admitted only for a planar, continuous phase law
//! with one-sided analytic gradients/Hessians, and a transverse central ray.
//! It does not supply a dielectric-index jump or curved-interface mapping.

use super::*;
use crate::gradient_index::event::{BoundarySide, GradientDeparture, GradientEventOptions};

#[derive(Clone, Copy)]
pub(super) struct PlaneCoordinates {
    normal: Vec3,
    origin: Point3,
    min_cosine: f64,
    clear_distance_m: f64,
}

#[derive(Clone, Copy)]
pub(super) struct PlaneClockTrial {
    pub(super) value: [f64; 2],
    error: [f64; 2],
    projection: [f64; 2],
}

impl PlaneCoordinates {
    pub(super) fn new(
        departure: GradientDeparture,
        options: GradientEventOptions,
    ) -> Result<Self, GradientIndexError> {
        let sign = match departure.side {
            BoundarySide::Positive => 1.0,
            BoundarySide::Negative => -1.0,
            BoundarySide::OnSurface => {
                return Err(GradientIndexError::UnsupportedBoundaryDeparture);
            }
        };
        let normal = unit(departure.witness.normal)
            .ok_or(GradientIndexError::InvalidBoundaryQuery)?
            .scale(sign);
        Ok(Self {
            normal,
            origin: departure.witness.surface_point,
            min_cosine: options.min_approach_cosine,
            clear_distance_m: 2.0 * options.surface_tolerance_m,
        })
    }

    fn cosine(self, direction: Vec3) -> Result<f64, GradientIndexError> {
        let cosine = self.normal.dot(direction);
        if !cosine.is_finite() || cosine < self.min_cosine {
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        Ok(cosine)
    }

    pub(super) fn has_left_collar(self, position: Point3) -> bool {
        self.normal.dot(position.delta_from(self.origin)) > self.clear_distance_m
    }

    pub(super) fn enter(
        self,
        state: &mut BundleState,
        jet: IndexJet,
    ) -> Result<[f64; 2], GradientIndexError> {
        let (u, k) = motion(state, jet)?;
        let cosine = self.cosine(u)?;
        let mut clock = [0.0; 2];
        for (i, tau) in clock.iter_mut().enumerate() {
            let x = read_vec(state, 12 + 6 * i);
            let v = read_vec(state, 15 + 6 * i);
            *tau = -self.normal.dot(x) / cosine;
            write_vec(state, 12 + 6 * i, sub(x, u.scale(-*tau)));
            write_vec(state, 15 + 6 * i, sub(v, k.scale(-*tau)));
        }
        if state.iter().chain(clock.iter()).any(|x| !x.is_finite()) {
            return Err(GradientIndexError::NumericalRange);
        }
        Ok(clock)
    }

    pub(super) fn rhs(
        self,
        state: &BundleState,
        jet: IndexJet,
        track_group: bool,
    ) -> Result<(BundleState, [f64; 2]), GradientIndexError> {
        // Reuse the native optical equations and analytic A*Q+B*W. Only
        // the change of independent coordinate modifies the variations.
        let mut derivative = bundle_rhs(state, jet, track_group)?;
        let u = read_vec(&derivative, 0);
        let k = read_vec(&derivative, 3);
        let cosine = self.cosine(u)?;
        let mut clock_derivative = [0.0; 2];
        for (i, tau_prime) in clock_derivative.iter_mut().enumerate() {
            let w = read_vec(&derivative, 12 + 6 * i);
            let q = self.normal.dot(w) / cosine;
            let v_prime = read_vec(&derivative, 15 + 6 * i);
            write_vec(&mut derivative, 12 + 6 * i, sub(w, u.scale(q)));
            write_vec(&mut derivative, 15 + 6 * i, sub(v_prime, k.scale(q)));
            *tau_prime = -q;
        }
        if derivative
            .iter()
            .chain(clock_derivative.iter())
            .any(|x| !x.is_finite())
        {
            return Err(GradientIndexError::NumericalRange);
        }
        Ok((derivative, clock_derivative))
    }

    pub(super) fn clock_trial(
        self,
        state: &mut BundleState,
        jet: IndexJet,
        initial: [f64; 2],
        stages: [[f64; 2]; 7],
        step: f64,
    ) -> Result<PlaneClockTrial, GradientIndexError> {
        let (u, k) = motion(state, jet)?;
        let cosine = self.cosine(u)?;
        let mut trial = PlaneClockTrial {
            value: initial,
            error: [0.0; 2],
            projection: [0.0; 2],
        };
        for i in 0..2 {
            let mut increment = 0.0;
            for j in 0..6 {
                increment += A[5][j] * stages[j][i];
            }
            trial.value[i] += step * increment;
            for j in 0..7 {
                trial.error[i] += ERROR[j] * stages[j][i];
            }
            trial.error[i] *= step;
            // Correct floating-point drift in N.Q=0 by a coordinate change
            // that leaves the represented fixed-arc X,V unchanged.
            let q = read_vec(state, 12 + 6 * i);
            let w = read_vec(state, 15 + 6 * i);
            let correction = self.normal.dot(q) / cosine;
            write_vec(state, 12 + 6 * i, sub(q, u.scale(correction)));
            write_vec(state, 15 + 6 * i, sub(w, k.scale(correction)));
            trial.value[i] -= correction;
            trial.projection[i] = correction.abs();
        }
        if trial
            .value
            .iter()
            .chain(&trial.error)
            .chain(&trial.projection)
            .any(|x| !x.is_finite())
        {
            return Err(GradientIndexError::NumericalRange);
        }
        Ok(trial)
    }

    pub(super) fn fixed_arc(
        self,
        state: &BundleState,
        jet: IndexJet,
        clock: [f64; 2],
    ) -> Result<[AngularDerivative; 2], GradientIndexError> {
        let (u, k) = motion(state, jet)?;
        self.cosine(u)?;
        let columns = core::array::from_fn(|i| {
            let x = sub(read_vec(state, 12 + 6 * i), u.scale(clock[i]));
            let v = sub(read_vec(state, 15 + 6 * i), k.scale(clock[i]));
            AngularDerivative {
                position_m_per_rad: x,
                direction_per_rad: sub(v, u.scale(u.dot(v))),
            }
        });
        if columns
            .iter()
            .any(|c| !finite_vec(c.position_m_per_rad) || !finite_vec(c.direction_per_rad))
        {
            return Err(GradientIndexError::NumericalRange);
        }
        Ok(columns)
    }

    pub(super) fn error(
        self,
        before: &BundleState,
        before_jet: IndexJet,
        before_clock: [f64; 2],
        trial: &BundleTrial,
        options: BundleOptions,
    ) -> Result<f64, GradientIndexError> {
        let clock = trial.plane_clock.ok_or(GradientIndexError::InvalidSample)?;
        let a = self.fixed_arc(before, before_jet, before_clock)?;
        let b = self.fixed_arc(&trial.state, trial.jet, clock.value)?;
        let (u, k) = motion(&trial.state, trial.jet)?;
        let g = trial
            .jet
            .optical
            .phase_gradient_per_m
            .scale(1.0 / trial.jet.optical.phase_index);
        let dx = read_vec(&trial.error, 0);
        let du = read_vec(&trial.error, 3);
        let du = sub(du, u.scale(u.dot(du)));
        let kg = sub(
            trial
                .jet
                .hessian
                .apply(dx)
                .scale(1.0 / trial.jet.optical.phase_index),
            g.scale(g.dot(dx)),
        );
        let dk = sub(
            sub(sub(kg, u.scale(u.dot(kg))), du.scale(u.dot(g))),
            u.scale(g.dot(du)),
        );
        let mut error = 0.0_f64;
        for i in 0..2 {
            // tau has units metres/radian and participates in its OWN embedded
            // estimate. Also control the conversion back to physical columns,
            // including central tangent/curvature errors multiplied by tau.
            let scale = options.position_derivative_tolerance_m
                + options.integration.relative_tolerance
                    * before_clock[i].abs().max(clock.value[i].abs());
            error = error.max(ratio(clock.error[i].abs().max(clock.projection[i]), scale)?);
            let ex = sub(
                sub(read_vec(&trial.error, 12 + 6 * i), u.scale(clock.error[i])),
                du.scale(clock.value[i]),
            );
            let ev = sub(
                sub(read_vec(&trial.error, 15 + 6 * i), k.scale(clock.error[i])),
                dk.scale(clock.value[i]),
            );
            for (estimate, old, new, absolute) in [
                (
                    ex,
                    a[i].position_m_per_rad,
                    b[i].position_m_per_rad,
                    options.position_derivative_tolerance_m,
                ),
                (
                    ev,
                    a[i].direction_per_rad,
                    b[i].direction_per_rad,
                    options.direction_derivative_tolerance,
                ),
            ] {
                for (e, x, y) in [
                    (estimate.x, old.x, new.x),
                    (estimate.y, old.y, new.y),
                    (estimate.z, old.z, new.z),
                ] {
                    let scale =
                        absolute + options.integration.relative_tolerance * x.abs().max(y.abs());
                    error = error.max(ratio(e.abs(), scale)?);
                }
            }
        }
        Ok(error)
    }
}

fn motion(state: &BundleState, jet: IndexJet) -> Result<(Vec3, Vec3), GradientIndexError> {
    let derivative = rhs(&optical_state(state), jet.optical, false)?;
    Ok((get_vec(&derivative, 0), get_vec(&derivative, 3)))
}
fn ratio(error: f64, scale: f64) -> Result<f64, GradientIndexError> {
    let result = error / scale;
    if !scale.is_finite() || scale <= 0.0 || !result.is_finite() {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(result)
}
