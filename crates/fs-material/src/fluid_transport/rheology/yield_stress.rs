//! Original-card Herschel–Bulkley pressure/Couette flow between planar walls.
//!
//! The fully developed, homogeneous, incompressible laminar reduction uses
//! `tau(y)=c+dp*(h/2-y)/L` and the original slurry card's set-valued yield law.
//! One bounded scalar solve enforces both wall velocities; yielded-layer
//! integrals are analytic, including the integrable differential singularity
//! of shear-thickening cards at a yield surface. No regularization, apparent
//! plug viscosity, superposition of nonlinear flows, or fitted surrogate is
//! introduced. The temperature and solids mass fraction select support of the
//! original constant-parameter card; they do not imply an extra thermal law.
//!
//! This is a wide planar passage: side walls, entrance/inertial effects,
//! cavitation, thermal evolution, particle migration and wall slip are outside
//! this reduction. The caller explicitly declares a laminar shear-rate ceiling.

use crate::slurry::rheology::{HerschelBulkleyCard, ShearResponse, SlurryConditions, SlurryError};
use fs_math::det;
use fs_qty::{Frequency, Length, Pressure, Stress, Velocity};

mod moments;
use moments::{Moments, distribution, product};

fn require(ok: bool, message: &'static str) -> Result<(), SlurryError> {
    if ok {
        Ok(())
    } else {
        Err(SlurryError::Invalid(message))
    }
}

fn outside(axis: &'static str, value: f64, maximum: f64) -> SlurryError {
    SlurryError::OutsideValidity {
        axis,
        value,
        minimum: 0.0,
        maximum,
    }
}

/// Numerical closure of the relative wall speed, independent of common motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressSlitControls {
    /// Relative tolerance on the nonzero imposed difference of wall velocities.
    pub relative_wall_velocity_tolerance: f64,
    /// Safeguarded stress-offset iterations; permitted range is 1 through 128.
    pub maximum_iterations: u32,
}

impl Default for YieldStressSlitControls {
    fn default() -> Self {
        Self {
            relative_wall_velocity_tolerance: 1.0e-11,
            maximum_iterations: 96,
        }
    }
}

impl YieldStressSlitControls {
    fn validate(self) -> Result<(), SlurryError> {
        require(
            self.relative_wall_velocity_tolerance.is_finite()
                && (16.0 * f64::EPSILON..=1.0e-3).contains(&self.relative_wall_velocity_tolerance)
                && (1..=128).contains(&self.maximum_iterations),
            "invalid yield-stress slit solve controls",
        )
    }
}

/// Differential response at fixed temperature and composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressSlitTangent {
    /// dQ/d(dp) with both wall velocities fixed [m^3/(s Pa)].
    pub pressure_conductance_m3_s_pa: f64,
    /// dQ/d(U_upper-U_lower), holding mean wall velocity and dp fixed [m^2].
    pub flow_relative_wall_velocity_derivative_m2: f64,
    /// d(w L c)/d(U_upper-U_lower) at fixed dp [N s/m].
    pub relative_wall_resistance_n_s_m: f64,
}

/// Coupled hydraulic, wall-work and heat accounts for one admitted solution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressSlitFlow {
    /// Inlet minus outlet pressure [Pa].
    pub pressure_drop_pa: f64,
    /// Signed volume flow, including common-wall entrainment [m^3/s].
    pub volume_flow_m3_s: f64,
    /// Flow after removing `width*gap*mean_wall_velocity` [m^3/s].
    pub relative_volume_flow_m3_s: f64,
    /// Requested lower/upper wall velocities [m/s].
    pub requested_wall_velocities_m_s: [f64; 2],
    /// Wall velocities actually closed by the returned profile [m/s].
    pub accepted_wall_velocities_m_s: [f64; 2],
    /// Accepted minus requested relative wall velocity [m/s].
    pub relative_wall_velocity_residual_m_s: f64,
    /// Stress offset c at the midplane [Pa]. For a static plug this is a
    /// symmetric representative of `static_midplane_stress_interval_pa`.
    pub midplane_stress_pa: f64,
    /// Lower/upper shear tractions exerted by the walls on the fluid [Pa].
    pub wall_tractions_pa: [f64; 2],
    /// Inclusive traction intervals [Pa], distinct from a selected static
    /// representative. Flowing states have equal interval endpoints.
    pub wall_traction_bounds_pa: [[f64; 2]; 2],
    /// Whole-gap static stress offsets, when the relative wall speed is zero.
    /// The lower wall traction is `-c-p` and the upper is `c-p`, p=dp*h/(2L).
    pub static_midplane_stress_interval_pa: Option<[f64; 2]>,
    /// Physical lower/upper boundaries of the unyielded interval [m], if it
    /// has nonzero width. An isolated zero-rate surface is not a finite plug.
    pub unyielded_interval_m: Option<[f64; 2]>,
    /// Pressure work `dp*Q` [W], including work on a translating rigid plug.
    pub pressure_work_w: f64,
    /// Work of wall tractions at the requested wall speeds [W].
    pub wall_work_w: f64,
    /// Nonnegative original-law volume integral of `tau*du/dy` [W].
    pub dissipation_w: f64,
    /// Pressure work plus requested wall work minus dissipation [W], evaluated
    /// in relative coordinates to avoid cancelling common translation.
    pub work_balance_residual_w: f64,
    /// Finite differentiable response, or None for a static yield interval or
    /// a zero-yield power-law origin with a singular/degenerate differential.
    pub tangent: Option<YieldStressSlitTangent>,
    /// Number of bounded nonlinear stress-offset iterations performed.
    pub stress_iterations: u32,
}

/// Physical sample: a rigid plug has exactly zero gradient and zero heating.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressSlitSample {
    /// Distance above the lower wall [m].
    pub position_m: f64,
    /// Axial fluid velocity [m/s].
    pub velocity_m_s: f64,
    /// Signed shear stress [Pa].
    pub shear_stress_pa: f64,
    /// Signed physical velocity gradient du/dy [1/s].
    pub shear_rate_per_s: f64,
    /// Local irreversible heating [W/m^3].
    pub dissipation_w_m3: f64,
    /// True exactly where the original stress lies inside the yield interval.
    pub unyielded: bool,
}

/// Prepared geometry and original source, with explicit source/laminar bounds.
#[derive(Clone, Copy, Debug)]
pub struct LaminarHerschelBulkleySlit<'a> {
    source: &'a HerschelBulkleyCard,
    conditions: SlurryConditions,
    width: f64,
    gap: f64,
    length: f64,
    yield_stress: f64,
    inverse_index: f64,
    maximum_rate: f64,
    maximum_stress: f64,
}

#[derive(Clone, Copy, Debug)]
struct Segment {
    wall: usize,
    sign: f64,
    fraction: f64,
    high_excess: f64,
    low_excess: f64,
    stress_span: f64,
    high_rate: f64,
}

/// Immutable solved profile; queries retain the original card and conditions.
#[derive(Clone, Debug)]
pub struct YieldStressSlitProfile<'a> {
    slit: LaminarHerschelBulkleySlit<'a>,
    flow: YieldStressSlitFlow,
    half_pressure_stress: f64,
    wall_stresses: [f64; 2],
    segments: [Option<Segment>; 2],
    uniform_rate: Option<f64>,
}

impl<'a> LaminarHerschelBulkleySlit<'a> {
    /// Prepare a no-slip planar passage using one original slurry card.
    /// Source support must include zero rate, because plugs and yield surfaces
    /// cannot be filled by extrapolating a positive-rate-only calibration.
    ///
    /// # Errors
    /// Refuses unsupported conditions, invalid dimensions/rate ceilings and
    /// unrepresentable source responses. No material or laminar limit is inferred.
    pub fn new(
        source: &'a HerschelBulkleyCard,
        conditions: SlurryConditions,
        width: Length,
        gap: Length,
        length: Length,
        maximum_laminar_shear_rate: Frequency,
    ) -> Result<Self, SlurryError> {
        let (width, gap, length) = (width.value(), gap.value(), length.value());
        require(
            [width, gap, length, maximum_laminar_shear_rate.value()]
                .into_iter()
                .all(|x| x.is_finite() && x > 0.0),
            "slit dimensions and laminar shear-rate ceiling must be positive",
        )?;
        let support = source.domain().shear_rate_bounds();
        require(
            support[0].value() == 0.0,
            "yield-stress slit requires source support down to rest",
        )?;
        source.stress_at_shear_rate(Frequency::new(0.0), conditions)?;
        let maximum_rate = maximum_laminar_shear_rate.value().min(support[1].value());
        let maximum_stress =
            match source.stress_at_shear_rate(Frequency::new(maximum_rate), conditions)? {
                ShearResponse::Yielded { stress, .. } => stress.value(),
                ShearResponse::Unyielded { .. } => {
                    return Err(SlurryError::Invalid("positive slit ceiling did not yield"));
                }
            };
        let parameters = source.parameters();
        Ok(Self {
            source,
            conditions,
            width,
            gap,
            length,
            yield_stress: parameters.yield_stress.value(),
            inverse_index: 1.0 / parameters.flow_index.value(),
            maximum_rate,
            maximum_stress,
        })
    }

    /// Effective inclusive source/laminar shear-rate ceiling [1/s].
    #[must_use]
    pub const fn maximum_shear_rate_per_s(&self) -> f64 {
        self.maximum_rate
    }

    /// Solve simultaneous pressure and lower/upper wall motion.
    ///
    /// # Errors
    /// Refuses source/laminar violations, invalid inputs, unresolved numerical
    /// responses and exhausted offset iterations; no partial profile is returned.
    pub fn profile_at_pressure(
        &self,
        pressure_drop: Pressure,
        wall_velocities: [Velocity; 2],
        controls: YieldStressSlitControls,
    ) -> Result<YieldStressSlitProfile<'a>, SlurryError> {
        self.profile_at_pressure_with_check(pressure_drop, wall_velocities, controls, || Ok(()))
    }

    /// Bounded, cancellable solve. The caller's check is polled before work,
    /// each offset iteration, each yielded layer, and before publishing output.
    /// The source and prepared geometry are immutable on every refusal.
    ///
    /// # Errors
    /// Propagates cancellation unchanged in addition to `profile_at_pressure` errors.
    pub fn profile_at_pressure_with_check(
        &self,
        pressure_drop: Pressure,
        wall_velocities: [Velocity; 2],
        controls: YieldStressSlitControls,
        mut check: impl FnMut() -> Result<(), SlurryError>,
    ) -> Result<YieldStressSlitProfile<'a>, SlurryError> {
        check()?;
        controls.validate()?;
        let dp = pressure_drop.value();
        let walls = wall_velocities.map(|u| u.value());
        require(
            dp.is_finite() && walls.iter().all(|u| u.is_finite()),
            "slit loads must be finite",
        )?;
        let relative = walls[1] - walls[0];
        require(
            relative.is_finite(),
            "relative wall speed is unrepresentable",
        )?;
        let p = dp.signum() * product(&[dp.abs(), self.gap], &[2.0, self.length])?;
        if p.abs() > self.maximum_stress {
            return Err(outside(
                "slit pressure half-stress",
                p.abs(),
                self.maximum_stress,
            ));
        }
        let target_rate = product(&[relative.abs()], &[self.gap])?;
        if target_rate > self.maximum_rate {
            return Err(outside(
                "slit mean shear-rate magnitude",
                target_rate,
                self.maximum_rate,
            ));
        }
        let (c, accepted_relative, uniform_rate, iterations) = if p == 0.0 {
            let rate = relative.signum() * target_rate;
            let c = if relative == 0.0 {
                0.0
            } else {
                self.stress_at_rate(rate)?
            };
            (c, relative, Some(rate), 0)
        } else if relative == 0.0 {
            (0.0, 0.0, None, 0)
        } else {
            let (magnitude, actual_rate, iterations) =
                self.offset(p.abs(), target_rate, controls, &mut check)?;
            (
                relative.signum() * magnitude,
                relative.signum() * product(&[actual_rate, self.gap], &[])?,
                None,
                iterations,
            )
        };
        let mut wall_stresses = [c + p, c - p];
        // Preserve the admitted source/caller endpoint instead of sending a
        // rounded forward/inverse round trip just beyond the rate ceiling.
        if c != 0.0 && c.abs() == self.maximum_stress - p.abs() {
            let wall = usize::from(c.signum() != p.signum());
            wall_stresses[wall] = c.signum() * self.maximum_stress;
        }
        for stress in wall_stresses {
            self.rate_at_stress(stress)?;
        }
        let mut segments = [None, None];
        let moments = if let Some(rate) = uniform_rate {
            self.uniform_moments(rate)?
        } else {
            let mut moments = Moments::default();
            for wall in 0..2 {
                check()?;
                if let Some(segment) = self.segment(wall, wall_stresses[wall], p)? {
                    self.add_segment_moments(segment, &mut moments)?;
                    segments[wall] = Some(segment);
                }
            }
            moments
        };
        let flow = self.accounts(
            dp,
            walls,
            relative,
            accepted_relative,
            c,
            p,
            wall_stresses,
            &segments,
            moments,
            iterations,
            controls.relative_wall_velocity_tolerance,
        )?;
        check()?;
        Ok(YieldStressSlitProfile {
            slit: *self,
            flow,
            half_pressure_stress: p,
            wall_stresses,
            segments,
            uniform_rate,
        })
    }

    fn stress_at_rate(&self, rate: f64) -> Result<f64, SlurryError> {
        if rate.abs() > self.maximum_rate {
            return Err(outside(
                "slit shear-rate magnitude",
                rate.abs(),
                self.maximum_rate,
            ));
        }
        match self
            .source
            .stress_at_shear_rate(Frequency::new(rate), self.conditions)?
        {
            ShearResponse::Yielded { stress, .. } => Ok(stress.value()),
            ShearResponse::Unyielded { .. } => Ok(0.0),
        }
    }

    fn rate_at_stress(&self, stress: f64) -> Result<f64, SlurryError> {
        if !stress.is_finite() || stress.abs() > self.maximum_stress {
            return Err(outside(
                "slit shear-stress magnitude",
                stress.abs(),
                self.maximum_stress,
            ));
        }
        if stress.abs() == self.maximum_stress {
            return Ok(stress.signum() * self.maximum_rate);
        }
        let rate = self
            .source
            .shear_rate_at_stress(Stress::new(stress), self.conditions)?
            .value();
        if rate.abs() > self.maximum_rate {
            return Err(outside(
                "slit shear-rate magnitude",
                rate.abs(),
                self.maximum_rate,
            ));
        }
        Ok(rate)
    }

    // Positive-stress interval average. The known span survives subtraction of
    // nearly equal endpoints, and the true yield interval contributes zero.
    fn average(
        &self,
        low: f64,
        high: f64,
        span: f64,
        differential: bool,
    ) -> Result<f64, SlurryError> {
        if high <= self.yield_stress {
            return Ok(0.0);
        }
        let high_excess = high - self.yield_stress;
        let rate = self.rate_at_stress(high)?;
        if differential && low <= self.yield_stress {
            return product(&[rate], &[span]);
        }
        let low_excess = (low - self.yield_stress).max(0.0);
        let yielded_span = if low_excess == 0.0 { high_excess } else { span };
        let a = self.inverse_index + if differential { 0.0 } else { 1.0 };
        let d = distribution(a, low_excess, high_excess, yielded_span)?;
        if differential {
            product(&[self.inverse_index, rate, d.average], &[high_excess])
        } else {
            product(&[yielded_span, rate, d.average], &[span])
        }
    }

    fn mean_rate(&self, c: f64, p: f64) -> Result<(f64, f64), SlurryError> {
        if c == 0.0 {
            return Ok((0.0, product(&[self.rate_at_stress(p)?], &[p])?));
        }
        let high = if c == self.maximum_stress - p {
            self.maximum_stress
        } else {
            c + p
        };
        let low = (c - p).abs();
        let span = product(&[2.0, c.min(p)], &[])?;
        if high <= self.yield_stress {
            return Ok((0.0, 0.0));
        }
        let average = self.average(low, high, span, false)?;
        let mean = if c <= p {
            product(&[average, c], &[p])?
        } else {
            average
        };
        let differential = if c <= p {
            product(
                &[0.5 * self.rate_at_stress(high)? + 0.5 * self.rate_at_stress(low)?],
                &[p],
            )?
        } else {
            self.average(low, high, span, true)?
        };
        Ok((mean, differential))
    }

    fn offset(
        &self,
        p: f64,
        target: f64,
        controls: YieldStressSlitControls,
        check: &mut impl FnMut() -> Result<(), SlurryError>,
    ) -> Result<(f64, f64, u32), SlurryError> {
        let mut upper = self.maximum_stress - p;
        let (maximum, _) = self.mean_rate(upper, p)?;
        if target > maximum {
            return Err(outside(
                "source-bounded slit mean shear rate",
                target,
                maximum,
            ));
        }
        if target == maximum {
            return Ok((upper, maximum, 0));
        }
        let mut lower = (self.yield_stress - p).max(0.0);
        let pure = self.stress_at_rate(target)?;
        let mut c = pure.min(upper).max(lower);
        for iteration in 1..=controls.maximum_iterations {
            check()?;
            let (mean, differential) = self.mean_rate(c, p)?;
            let residual = mean - target;
            if residual.abs() <= controls.relative_wall_velocity_tolerance * target {
                return Ok((c, mean, iteration));
            }
            if residual > 0.0 {
                upper = c;
            } else {
                lower = c;
            }
            let newton = c - residual / differential;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            require(
                next != c,
                "slit stress offset stalled before relative wall closure",
            )?;
            c = next;
        }
        Err(SlurryError::NoConvergence {
            iterations: controls.maximum_iterations,
        })
    }

    fn segment(&self, wall: usize, stress: f64, p: f64) -> Result<Option<Segment>, SlurryError> {
        let descending_sign = if wall == 0 { p.signum() } else { -p.signum() };
        if stress.abs() <= self.yield_stress || stress.signum() != descending_sign {
            return Ok(None);
        }
        let high_excess = stress.abs() - self.yield_stress;
        let full_span = product(&[2.0, p.abs()], &[])?;
        let stress_span = high_excess.min(full_span);
        let fraction = if stress_span == full_span {
            1.0
        } else {
            product(&[stress_span], &[full_span])?
        };
        let low_excess = if stress_span == high_excess {
            0.0
        } else {
            high_excess - stress_span
        };
        Ok(Some(Segment {
            wall,
            sign: stress.signum(),
            fraction,
            high_excess,
            low_excess,
            stress_span,
            high_rate: self.rate_at_stress(stress)?.abs(),
        }))
    }

    fn add_segment_moments(
        &self,
        segment: Segment,
        moments: &mut Moments,
    ) -> Result<(), SlurryError> {
        let d = distribution(
            self.inverse_index,
            segment.low_excess,
            segment.high_excess,
            segment.stress_span,
        )?;
        let weight = if segment.low_excess == 0.0 {
            product(
                &[segment.fraction, segment.high_rate],
                &[segment.high_excess],
            )?
        } else {
            product(
                &[
                    segment.fraction,
                    self.inverse_index,
                    segment.high_rate,
                    d.average,
                ],
                &[segment.high_excess],
            )?
        };
        let mean = (if segment.wall == 0 { 1.0 } else { -1.0 }) * (0.5 - segment.fraction * d.mean);
        let variance = product(&[segment.fraction, segment.fraction, d.variance], &[])?;
        let wall_moment = product(
            &[
                segment.fraction,
                (1.0 - segment.fraction) * d.mean + segment.fraction * d.wall_moment,
            ],
            &[],
        )?;
        moments.add(weight, mean, variance, wall_moment)
    }

    fn uniform_moments(&self, rate: f64) -> Result<Moments, SlurryError> {
        let parameters = self.source.parameters();
        let differential = if rate == 0.0 {
            if self.yield_stress != 0.0 || parameters.flow_index.value() != 1.0 {
                return Ok(Moments::default());
            }
            product(
                &[parameters.reference_shear_rate.value()],
                &[parameters.reference_excess_stress.value()],
            )?
        } else {
            let excess = product(
                &[
                    parameters.reference_excess_stress.value(),
                    det::pow(
                        rate.abs() / parameters.reference_shear_rate.value(),
                        parameters.flow_index.value(),
                    ),
                ],
                &[],
            )?;
            product(&[self.inverse_index, rate.abs()], &[excess])?
        };
        Ok(Moments {
            weight: differential,
            mean: 0.0,
            variance: product(&[differential], &[12.0])?,
            wall_weight: product(&[differential], &[6.0])?,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn accounts(
        &self,
        dp: f64,
        walls: [f64; 2],
        requested_relative: f64,
        accepted_relative: f64,
        c: f64,
        p: f64,
        wall_stresses: [f64; 2],
        segments: &[Option<Segment>; 2],
        moments: Moments,
        iterations: u32,
        relative_tolerance: f64,
    ) -> Result<YieldStressSlitFlow, SlurryError> {
        let mean_wall = 0.5 * walls[0] + 0.5 * walls[1];
        let accepted_walls = if accepted_relative == requested_relative {
            walls
        } else {
            [
                mean_wall - 0.5 * accepted_relative,
                mean_wall + 0.5 * accepted_relative,
            ]
        };
        let represented_relative = accepted_walls[1] - accepted_walls[0];
        require(
            represented_relative.is_finite()
                && (represented_relative - accepted_relative).abs()
                    <= relative_tolerance * accepted_relative.abs().max(requested_relative.abs()),
            "common wall motion hides the accepted relative slit velocity",
        )?;
        let relative_flow = p.signum()
            * product(
                &[self.width, self.gap, self.gap, p.abs(), moments.wall_weight],
                &[],
            )?;
        let entrainment =
            mean_wall.signum() * product(&[self.width, self.gap, mean_wall.abs()], &[])?;
        let volume_flow = entrainment + relative_flow;
        let relative_residual = accepted_relative - requested_relative;
        let shear_heat = product(
            &[self.width, self.length, c.abs(), accepted_relative.abs()],
            &[],
        )?;
        let pressure_heat = product(&[dp.abs(), relative_flow.abs()], &[])?;
        let dissipation = shear_heat + pressure_heat;
        let pressure_work = dp * volume_flow;
        let wall_work = self.width * self.length * (c * requested_relative - 2.0 * p * mean_wall);
        let work_residual = -c.signum()
            * relative_residual.signum()
            * product(
                &[self.width, self.length, c.abs(), relative_residual.abs()],
                &[],
            )?;
        let is_static = requested_relative == 0.0 && p.abs() <= self.yield_stress;
        let static_interval = is_static.then(|| {
            let extent = self.yield_stress - p.abs();
            [-extent, extent]
        });
        let tractions = [-wall_stresses[0], wall_stresses[1]];
        let traction_bounds = if let Some([lo, hi]) = static_interval {
            [[-hi - p, -lo - p], [lo - p, hi - p]]
        } else {
            tractions.map(|t| [t, t])
        };
        let plug = if is_static {
            Some([0.0, self.gap])
        } else if p != 0.0 {
            let lower = segments[0].map_or(0.0, |s| s.fraction);
            let upper = 1.0 - segments[1].map_or(0.0, |s| s.fraction);
            (lower < upper && self.yield_stress > 0.0)
                .then_some([self.gap * lower, self.gap * upper])
        } else {
            None
        };
        let tangent = if moments.weight > 0.0 {
            Some(YieldStressSlitTangent {
                pressure_conductance_m3_s_pa: product(
                    &[self.width, self.gap, self.gap, self.gap, moments.variance],
                    &[self.length],
                )?,
                flow_relative_wall_velocity_derivative_m2: moments.mean.signum()
                    * product(&[self.width, self.gap, moments.mean.abs()], &[])?,
                relative_wall_resistance_n_s_m: product(
                    &[self.width, self.length],
                    &[self.gap, moments.weight],
                )?,
            })
        } else {
            None
        };
        require(
            [
                volume_flow,
                relative_flow,
                dissipation,
                pressure_work,
                wall_work,
                work_residual,
                accepted_walls[0],
                accepted_walls[1],
            ]
            .iter()
            .all(|v| v.is_finite()),
            "slit work or flow account overflow",
        )?;
        Ok(YieldStressSlitFlow {
            pressure_drop_pa: dp,
            volume_flow_m3_s: volume_flow,
            relative_volume_flow_m3_s: relative_flow,
            requested_wall_velocities_m_s: walls,
            accepted_wall_velocities_m_s: accepted_walls,
            relative_wall_velocity_residual_m_s: relative_residual,
            midplane_stress_pa: c,
            wall_tractions_pa: tractions,
            wall_traction_bounds_pa: traction_bounds,
            static_midplane_stress_interval_pa: static_interval,
            unyielded_interval_m: plug,
            pressure_work_w: pressure_work,
            wall_work_w: wall_work,
            dissipation_w: dissipation,
            work_balance_residual_w: work_residual,
            tangent,
            stress_iterations: iterations,
        })
    }
}

impl YieldStressSlitProfile<'_> {
    /// Complete flow/work accounts and differential response.
    #[must_use]
    pub const fn flow(&self) -> &YieldStressSlitFlow {
        &self.flow
    }

    /// Sample the physical no-slip profile at a distance above the lower wall.
    ///
    /// # Errors
    /// Refuses positions outside the gap and unrepresentable original-law samples.
    pub fn sample(&self, position: Length) -> Result<YieldStressSlitSample, SlurryError> {
        self.sample_with_check(position, || Ok(()))
    }

    /// Cancellable immutable profile query.
    ///
    /// # Errors
    /// Propagates the caller's check unchanged in addition to `sample` errors.
    pub fn sample_with_check(
        &self,
        position: Length,
        mut check: impl FnMut() -> Result<(), SlurryError>,
    ) -> Result<YieldStressSlitSample, SlurryError> {
        check()?;
        let y = position.value();
        require(
            y.is_finite() && (0.0..=self.slit.gap).contains(&y),
            "slit sample outside physical gap",
        )?;
        let x = y / self.slit.gap;
        let stress = if y == 0.0 {
            self.wall_stresses[0]
        } else if y == self.slit.gap {
            self.wall_stresses[1]
        } else {
            self.flow.midplane_stress_pa + self.half_pressure_stress * (1.0 - 2.0 * x)
        };
        let rate = if let Some(rate) = self.uniform_rate {
            rate
        } else {
            self.slit.rate_at_stress(stress)?
        };
        let walls = self.flow.accepted_wall_velocities_m_s;
        let velocity = if y == 0.0 {
            walls[0]
        } else if y == self.slit.gap {
            walls[1]
        } else if let Some(rate) = self.uniform_rate {
            walls[0] + rate * y
        } else if x <= 0.5 {
            walls[0] + self.slit.gap * self.integral(0.0, x, &mut check)?
        } else {
            walls[1] - self.slit.gap * self.integral(x, 1.0, &mut check)?
        };
        let dissipation = product(&[stress.abs(), rate.abs()], &[])?;
        require(
            velocity.is_finite() && stress * rate >= 0.0,
            "invalid physical slit sample",
        )?;
        check()?;
        Ok(YieldStressSlitSample {
            position_m: y,
            velocity_m_s: velocity,
            shear_stress_pa: stress,
            shear_rate_per_s: rate,
            dissipation_w_m3: dissipation,
            unyielded: stress.abs() <= self.slit.yield_stress,
        })
    }

    fn integral(
        &self,
        from: f64,
        to: f64,
        check: &mut impl FnMut() -> Result<(), SlurryError>,
    ) -> Result<f64, SlurryError> {
        let mut result = 0.0;
        for segment in self.segments.iter().flatten() {
            check()?;
            let (from, to) = if segment.wall == 0 {
                (from, to)
            } else {
                (1.0 - to, 1.0 - from)
            };
            let lower = from.max(0.0);
            let upper = to.min(segment.fraction);
            if lower >= upper {
                continue;
            }
            let fraction = (upper - lower) / segment.fraction;
            let drop = segment.stress_span / segment.high_excess;
            let remaining = 1.0 - drop * (lower / segment.fraction);
            let high = product(&[segment.high_excess, remaining], &[])?;
            let low = if upper == segment.fraction {
                segment.low_excess
            } else {
                segment.high_excess * (1.0 - drop * (upper / segment.fraction))
            };
            let span = if low == 0.0 {
                high
            } else {
                product(&[segment.stress_span, fraction], &[])?
            };
            let d = distribution(self.slit.inverse_index + 1.0, low, high, span)?;
            let rate = product(
                &[
                    segment.high_rate,
                    det::pow(remaining, self.slit.inverse_index),
                ],
                &[],
            )?;
            result += segment.sign * product(&[upper - lower, rate, d.average], &[])?;
        }
        require(result.is_finite(), "slit velocity integral overflow")?;
        Ok(result)
    }
}
