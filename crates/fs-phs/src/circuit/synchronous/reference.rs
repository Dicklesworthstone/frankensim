//! Torque references from the original machine law at an authored d current.
//!
//! For fixed i_d, torque is p*(lambda_d*i_q-lambda_q*i_d). The inverse of each
//! original piecewise-affine q curve is also piecewise affine, so every feasible
//! torque branch can be searched exactly, without a fitted torque map or Newton
//! iteration. Current-circle and source-support limits apply to the reference;
//! they do not certify the transient current of a controlled plant.

use super::{SynchronousMotor, SynchronousMotorParameters, invalid, poll, saturation::MagneticLaw};
use crate::{PhsError, PreparedStepError};
use fs_math::{dd::Dd, det};

/// A feasible constitutive current target, before any voltage/control dynamics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorqueCurrentReference {
    pub requested_torque_nm: f64,
    /// Torque of the ORIGINAL magnetic law at the returned current target.
    pub achievable_torque_nm: f64,
    pub current_dq_a: [f64; 2],
    /// Power-invariant d/q vector magnitude limit [A], not phase-peak current.
    pub current_limit_a: f64,
    /// True when current/source limits prevent the requested constitutive torque.
    /// A feasible affine root can retain a floating-point torque remainder.
    pub torque_limited: bool,
}

#[derive(Clone, Copy)]
struct Point {
    current: f64,
    torque: f64,
}

struct Search {
    requested: f64,
    d_current: f64,
    d_flux: f64,
    poles: f64,
    linear_coefficient: Option<f64>,
    best: Point,
    best_error: f64,
    has_root: bool,
}

fn finite(value: f64) -> Result<f64, PhsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "synchronous torque reference arithmetic is unrepresentable",
        ))
    }
}

impl Search {
    fn point(&self, current: f64, flux: f64) -> Result<Point, PhsError> {
        if let Some(coefficient) = self.linear_coefficient {
            return Ok(Point {
                current,
                torque: finite(coefficient * current)?,
            });
        }
        let torque = Dd::from_f64(self.poles)
            * (Dd::from_f64(self.d_flux) * Dd::from_f64(current)
                - Dd::from_f64(flux) * Dd::from_f64(self.d_current));
        Ok(Point {
            current,
            torque: finite(torque.to_f64())?,
        })
    }

    fn consider(&mut self, point: Point, root: bool) -> Result<(), PhsError> {
        let error =
            finite((Dd::from_f64(point.torque) - Dd::from_f64(self.requested)).to_f64())?.abs();
        let lower_current = point.current.abs() < self.best.current.abs()
            || point.current.abs() == self.best.current.abs() && point.current < self.best.current;
        let replace = if root {
            !self.has_root || lower_current
        } else {
            !self.has_root && (error < self.best_error || error == self.best_error && lower_current)
        };
        if replace {
            self.best = point;
            self.best_error = error;
        }
        self.has_root |= root;
        Ok(())
    }

    fn segment(
        &mut self,
        lo: Point,
        hi: Point,
        flux_at: impl Fn(f64) -> Result<f64, PhsError>,
    ) -> Result<(), PhsError> {
        self.consider(lo, lo.torque == self.requested)?;
        self.consider(hi, hi.torque == self.requested)?;
        if self.requested < lo.torque.min(hi.torque) || self.requested > lo.torque.max(hi.torque) {
            return Ok(());
        }
        let current = if lo.torque == hi.torque {
            // An entire branch can have zero torque authority. Its least
            // current representative is deterministic, including zero torque.
            0.0_f64.clamp(lo.current, hi.current)
        } else {
            let fraction = finite(
                ((Dd::from_f64(self.requested) - Dd::from_f64(lo.torque))
                    / (Dd::from_f64(hi.torque) - Dd::from_f64(lo.torque)))
                .to_f64(),
            )?;
            if !(0.0..=1.0).contains(&fraction) {
                return Err(invalid(
                    "synchronous torque root left its original affine segment",
                ));
            }
            let span = Dd::from_f64(hi.current) - Dd::from_f64(lo.current);
            // Anchor at the closer end, retaining small targets next to the
            // original zero-current knot without subtracting large currents.
            finite(if fraction <= 0.5 {
                (Dd::from_f64(lo.current) + Dd::from_f64(fraction) * span).to_f64()
            } else {
                (Dd::from_f64(hi.current) - Dd::from_f64(1.0 - fraction) * span).to_f64()
            })?
        };
        if current < lo.current || current > hi.current {
            return Err(invalid(
                "synchronous torque root is not representable inside its source segment",
            ));
        }
        let point = self.point(current, flux_at(current)?)?;
        let residual = finite(point.torque - self.requested)?.abs();
        let scale = lo
            .torque
            .abs()
            .max(hi.torque.abs())
            .max(self.requested.abs());
        if residual > 512.0 * f64::EPSILON * scale {
            return Err(invalid(
                "synchronous torque reference does not reproduce its original law",
            ));
        }
        self.consider(point, true)
    }
}

impl SynchronousMotor {
    /// Allocate q current for the requested rotor torque at a fixed authored d
    /// current, using this machine's original linear or bounded nonlinear law.
    ///
    /// All q-curve branches inside the current circle are considered. If several
    /// roots exist, choose minimum |i_q|, then the smaller signed current. If the
    /// torque is unreachable, return the closest achievable torque, with the
    /// same tie rule and `torque_limited=true`; a flat zero-authority machine
    /// therefore returns zero q current and explicitly limited nonzero torque.
    /// The authored d current itself is never silently changed or extrapolated.
    ///
    /// This is not MTPA (d current is fixed), field-weakening optimization, a
    /// voltage-feasibility test, or a transient overcurrent guarantee. A sampled
    /// controller and actual inverter still own voltage limiting and tracking.
    /// Search is allocation-free and bounded by the original q-curve knot count,
    /// with cancellation before work, at every source segment and before return.
    pub fn current_reference_for_torque_controlled(
        &self,
        requested_torque_nm: f64,
        d_current_a: f64,
        current_limit_a: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TorqueCurrentReference, PreparedStepError> {
        allocate(
            self.parameters,
            &self.magnetic_law,
            requested_torque_nm,
            d_current_a,
            current_limit_a,
            cancelled,
        )
    }
}

/// Shared by the cold machine and the prepared inverter that owns its law.
pub(super) fn allocate(
    parameters: SynchronousMotorParameters,
    magnetic_law: &MagneticLaw,
    requested_torque_nm: f64,
    d_current_a: f64,
    current_limit_a: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<TorqueCurrentReference, PreparedStepError> {
    poll(&mut cancelled)?;
    if [requested_torque_nm, d_current_a, current_limit_a]
        .iter()
        .any(|v| !v.is_finite())
        || current_limit_a < 0.0
        || d_current_a.abs() > current_limit_a
    {
        return Err(invalid(
            "finite torque/d-current and compatible nonnegative current magnitude required",
        )
        .into());
    }
    let q_limit = if current_limit_a == 0.0 {
        0.0
    } else {
        let ratio = d_current_a.abs() / current_limit_a;
        finite(current_limit_a * det::sqrt((1.0 - ratio) * (1.0 + ratio)))?
    };
    let d_flux = magnetic_law.flux_dq(parameters, [d_current_a, 0.0])?[0];
    let mut search = Search {
        requested: requested_torque_nm,
        d_current: d_current_a,
        d_flux,
        poles: f64::from(parameters.pole_pairs),
        linear_coefficient: if matches!(magnetic_law, MagneticLaw::Linear) {
            // Preserve the original exact zero-saliency/zero-PM nullspace;
            // rounded inverse flux products must not manufacture torque.
            Some(finite(
                f64::from(parameters.pole_pairs)
                    * (parameters.d_inductance_h - parameters.q_inductance_h)
                        .mul_add(d_current_a, parameters.pm_flux_linkage_wb),
            )?)
        } else {
            None
        },
        best: Point {
            current: 0.0,
            torque: 0.0,
        },
        best_error: requested_torque_nm.abs(),
        has_root: requested_torque_nm == 0.0,
    };
    match magnetic_law {
        MagneticLaw::Linear => {
            let flux_at = |current: f64| finite(parameters.q_inductance_h * current);
            let zero = search.point(0.0, 0.0)?;
            let lo = search.point(-q_limit, flux_at(-q_limit)?)?;
            let hi = search.point(q_limit, flux_at(q_limit)?)?;
            // Retain zero as an anchor even for the unbounded linear law.
            search.segment(lo, zero, flux_at)?;
            poll(&mut cancelled)?;
            search.segment(zero, hi, flux_at)?;
        }
        MagneticLaw::Curves { q, .. } => {
            let knots = q.knots();
            let lo = (-q_limit).max(knots[0].1);
            let hi = q_limit.min(knots[knots.len() - 1].1);
            let flux_at = |current: f64| {
                q.flux_wb(current)
                    .map_err(|_| invalid("synchronous q-current reference left its original curve"))
            };
            let mut previous = search.point(lo, flux_at(lo)?)?;
            for &(flux, current) in knots {
                poll(&mut cancelled)?;
                if current > lo && current < hi {
                    let next = search.point(current, flux)?;
                    search.segment(previous, next, flux_at)?;
                    previous = next;
                }
            }
            let last = search.point(hi, flux_at(hi)?)?;
            search.segment(previous, last, flux_at)?;
        }
    }
    poll(&mut cancelled)?;
    Ok(TorqueCurrentReference {
        requested_torque_nm,
        achievable_torque_nm: search.best.torque,
        current_dq_a: [d_current_a, search.best.current],
        current_limit_a,
        torque_limited: !search.has_root,
    })
}
