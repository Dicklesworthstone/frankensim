//! Circular Couette flow of an original Herschel--Bulkley source, including plugs.
//!
//! Momentum fixes tau(r)=tau_i*(a/r)^2 and physical shear is r*domega/dr.
//! Stresses below the source yield stress produce a rigid outer plug, with zero
//! heat. Equal wall speeds leave a static torque INTERVAL, not a measured zero
//! torque or a finite apparent viscosity. An imposed torque resolves that stress.
//!
//! The original source's stress inverse supplies the inner rate. With n the
//! source index, gamma=gamma_i*z^p and p=max(1,5/(n+1)), positive quadrature
//! integrates dOmega=gamma*d(tau)/(2*tau) without an endpoint singularity.
//! The exact mobility is dOmega/dtau_i=(gamma_i-gamma_o)/(2*tau_i), including
//! movement of the yield surface. Zero-yield power laws are analytic.
//! No solve or sample allocates; Newton, quadrature and cancellation are bounded.
//!
//! Formula context: Chatzimina, Georgiou & Alexandrou, Applied Rheology 19
//! (2009) 34288, DOI:10.3933/ApplRheol-19-34288. This is a steady, isothermal,
//! fully developed, concentric no-slip laminar reduction. The source's stated
//! temperature/composition apply unchanged. End effects, pressure load capacity,
//! inertia, Taylor instability, settling, slip and thermal evolution are absent.

use super::JournalGeometry;
use fs_material::slurry::rheology::{
    HerschelBulkleyCard, ShearResponse, SlurryConditions, SlurryError,
};
use fs_math::det;
use fs_qty::{Frequency, Stress};

#[derive(Clone, Debug, PartialEq)]
pub enum YieldStressJournalError {
    Source(SlurryError),
    Invalid(&'static str),
    Cancelled,
    QuadratureBudget,
    QuadratureResolution,
    StressSolveBudget,
    UnsupportedMotion,
}
impl From<SlurryError> for YieldStressJournalError {
    fn from(value: SlurryError) -> Self {
        Self::Source(value)
    }
}
impl core::fmt::Display for YieldStressJournalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "yield-stress journal: {self:?}")
    }
}
impl std::error::Error for YieldStressJournalError {}
type Result<T> = core::result::Result<T, YieldStressJournalError>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressJournalControls {
    pub maximum_stress_iterations: u32,
    pub maximum_quadrature_refinements: u32,
    pub maximum_quadrature_nodes: u64,
    pub relative_angular_tolerance: f64,
    pub absolute_angular_tolerance_rad_s: f64,
    pub relative_quadrature_tolerance: f64,
}
impl Default for YieldStressJournalControls {
    fn default() -> Self {
        Self {
            maximum_stress_iterations: 64,
            maximum_quadrature_refinements: 20,
            maximum_quadrature_nodes: 262_144,
            relative_angular_tolerance: 1e-10,
            absolute_angular_tolerance_rad_s: 1e-12,
            relative_quadrature_tolerance: 1e-11,
        }
    }
}
impl YieldStressJournalControls {
    fn admit(self) -> Result<()> {
        if !(1..=128).contains(&self.maximum_stress_iterations)
            || !(3..=24).contains(&self.maximum_quadrature_refinements)
            || self.maximum_quadrature_nodes == 0
            || !positive(self.absolute_angular_tolerance_rad_s)
            || !positive(self.relative_angular_tolerance)
            || self.relative_angular_tolerance >= 1.0
            || !positive(self.relative_quadrature_tolerance)
            || self.relative_quadrature_tolerance >= 1.0
        {
            return Err(YieldStressJournalError::Invalid(
                "invalid work or accuracy controls",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JournalYieldRegime {
    Unyielded,
    PartiallyYielded,
    FullyYielded,
}

/// Original source and fixed geometry/conditions; zero shear must be supported.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressJournal<'a> {
    source: &'a HerschelBulkleyCard,
    conditions: SlurryConditions,
    geometry: JournalGeometry,
    maximum_rate: f64,
    maximum_stress: f64,
    torque_factor: f64,
    static_torque: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressJournalFlow {
    pub regime: JournalYieldRegime,
    /// Actual outer-minus-inner angular speed, before common rotation is added.
    pub relative_angular_velocity_rad_s: f64,
    pub accepted_wall_angular_velocities_rad_s: [f64; 2],
    pub requested_wall_angular_velocities_rad_s: Option<[f64; 2]>,
    pub angular_velocity_residual_rad_s: f64,
    pub requested_outer_torque_nm: Option<f64>,
    pub outer_torque_residual_nm: Option<f64>,
    /// Torques ON THE FLUID. None means equal wall speeds did not determine
    /// which static stress is present; see the complete admissible interval.
    pub wall_torques_nm: Option<[f64; 2]>,
    pub static_outer_torque_bounds_nm: [f64; 2],
    /// Yielded material occupies [inner_radius, yielded_outer_radius].
    pub yielded_outer_radius_m: f64,
    pub wall_shear_rates_per_s: [f64; 2],
    pub wall_power_w: f64,
    pub dissipated_power_w: f64,
    pub work_residual_w: f64,
    /// d(outer torque)/d(Omega_outer-Omega_inner), for flowing states.
    /// Static states have no inferred finite apparent resistance.
    pub differential_resistance_nm_s: Option<f64>,
    pub quadrature_relative_change: f64,
    pub stress_iterations: u32,
    pub quadrature_nodes: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressJournalProfile<'a> {
    journal: YieldStressJournal<'a>,
    flow: YieldStressJournalFlow,
    inner_stress: Option<f64>,
    controls: YieldStressJournalControls,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YieldStressJournalSample {
    pub angular_velocity_rad_s: f64,
    pub azimuthal_velocity_m_s: f64,
    /// r*domega/dr; rigid rotation is excluded.
    pub shear_rate_per_s: f64,
    /// Unknown in the static velocity-controlled case, never invented as zero.
    pub radial_shear_stress_pa: Option<f64>,
    pub dissipation_w_m3: f64,
    pub yielded: bool,
}

#[derive(Clone, Copy)]
struct Evaluation {
    speed: f64,
    rates: [f64; 2],
    rate_drop_fraction: f64,
    change: f64,
}
struct Budget {
    used: u64,
    maximum: u64,
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(YieldStressJournalError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(YieldStressJournalError::Invalid("nonfinite response"))
    }
}
fn product_ratio(numerator: &[f64], denominator: &[f64]) -> Result<f64> {
    if numerator.iter().any(|x| !x.is_finite() || *x < 0.0)
        || denominator.iter().any(|x| !positive(*x))
    {
        return Err(YieldStressJournalError::Invalid("invalid product factors"));
    }
    if numerator.contains(&0.0) {
        return Ok(0.0);
    }
    let direct = numerator.iter().product::<f64>() / denominator.iter().product::<f64>();
    if positive(direct) {
        return Ok(direct);
    }
    let log = numerator.iter().map(|x| det::ln(*x)).sum::<f64>()
        - denominator.iter().map(|x| det::ln(*x)).sum::<f64>();
    let value = det::exp(log);
    if positive(value) {
        Ok(value)
    } else {
        Err(YieldStressJournalError::Invalid(
            "positive response is unrepresentable",
        ))
    }
}
fn log1p(x: f64) -> f64 {
    let rounded = 1.0 + x;
    det::ln(rounded) + (x - (rounded - 1.0)) / rounded
}
fn radial_drop(inner: f64, outer: f64) -> f64 {
    ((outer - inner) / outer) * (1.0 + inner / outer)
}

impl<'a> YieldStressJournal<'a> {
    pub fn new(
        source: &'a HerschelBulkleyCard,
        conditions: SlurryConditions,
        geometry: JournalGeometry,
        maximum_laminar_shear_rate_per_s: f64,
    ) -> Result<Self> {
        if !positive(geometry.inner_radius_m)
            || !positive(geometry.outer_radius_m)
            || !positive(geometry.length_m)
            || geometry.outer_radius_m <= geometry.inner_radius_m
            || !positive(maximum_laminar_shear_rate_per_s)
        {
            return Err(YieldStressJournalError::Invalid(
                "invalid geometry or physical shear ceiling",
            ));
        }
        source.shear_rate_at_stress(Stress::new(0.0), conditions)?;
        let ShearResponse::Yielded { stress, .. } = source
            .stress_at_shear_rate(Frequency::new(maximum_laminar_shear_rate_per_s), conditions)?
        else {
            return Err(YieldStressJournalError::Invalid(
                "positive ceiling must be yielded",
            ));
        };
        let factor = product_ratio(
            &[
                core::f64::consts::TAU,
                geometry.length_m,
                geometry.inner_radius_m,
                geometry.inner_radius_m,
            ],
            &[],
        )?;
        product_ratio(&[factor, stress.value()], &[])?;
        let static_torque =
            product_ratio(&[factor, source.parameters().yield_stress.value()], &[])?;
        Ok(Self {
            source,
            conditions,
            geometry,
            maximum_rate: maximum_laminar_shear_rate_per_s,
            maximum_stress: stress.value(),
            torque_factor: factor,
            static_torque,
        })
    }
    pub const fn source(&self) -> &'a HerschelBulkleyCard {
        self.source
    }
    pub const fn conditions(&self) -> SlurryConditions {
        self.conditions
    }
    pub const fn geometry(&self) -> JournalGeometry {
        self.geometry
    }
    pub const fn maximum_laminar_shear_rate_per_s(&self) -> f64 {
        self.maximum_rate
    }
    pub fn static_outer_torque_bounds_nm(&self) -> [f64; 2] {
        [-self.static_torque, self.static_torque]
    }

    pub fn at_wall_speeds(
        &self,
        walls: [f64; 2],
        controls: YieldStressJournalControls,
    ) -> Result<YieldStressJournalProfile<'a>> {
        self.with_cancel(walls, controls, || false)
    }
    pub fn with_cancel(
        &self,
        walls: [f64; 2],
        controls: YieldStressJournalControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<YieldStressJournalProfile<'a>> {
        controls.admit()?;
        poll(&mut cancelled)?;
        if walls.iter().any(|x| !x.is_finite()) {
            return Err(YieldStressJournalError::Invalid("nonfinite wall speed"));
        }
        let relative = finite(walls[1] - walls[0])?;
        let mut budget = Budget {
            used: 0,
            maximum: controls.maximum_quadrature_nodes,
        };
        if relative == 0.0 {
            return self.finish(
                (self.static_torque == 0.0).then_some(0.0),
                walls[0],
                Some(walls),
                None,
                Evaluation {
                    speed: 0.0,
                    rates: [0.0; 2],
                    rate_drop_fraction: 0.0,
                    change: 0.0,
                },
                controls,
                0,
                budget.used,
                &mut cancelled,
            );
        }
        let target = relative.abs();
        let maximum = self.evaluate(
            self.maximum_stress,
            self.geometry.outer_radius_m,
            controls,
            &mut budget,
            &mut cancelled,
        )?;
        if target > maximum.speed {
            return Err(YieldStressJournalError::UnsupportedMotion);
        }
        if target == maximum.speed {
            return self.finish(
                Some(self.maximum_stress.copysign(relative)),
                walls[0],
                Some(walls),
                None,
                maximum,
                controls,
                0,
                budget.used,
                &mut cancelled,
            );
        }
        let yield_stress = self.source.parameters().yield_stress.value();
        let n = self.source.parameters().flow_index.value();
        let (mut low, mut high) = (yield_stress, self.maximum_stress);
        let fraction = det::pow(target / maximum.speed, n / (n + 1.0));
        let mut stress = (high - low).mul_add(fraction, low);
        if stress <= low || stress >= high {
            stress = 0.5 * low + 0.5 * high;
        }
        for iteration in 1..=controls.maximum_stress_iterations {
            poll(&mut cancelled)?;
            let value = self.evaluate(
                stress,
                self.geometry.outer_radius_m,
                controls,
                &mut budget,
                &mut cancelled,
            )?;
            let defect = value.speed - target;
            let mobility =
                product_ratio(&[0.5, value.rates[0], value.rate_drop_fraction], &[stress])?;
            if !positive(mobility) {
                return Err(YieldStressJournalError::UnsupportedMotion);
            }
            let correction = defect / mobility;
            let tolerance = controls.absolute_angular_tolerance_rad_s
                + controls.relative_angular_tolerance * target;
            if defect.abs() <= tolerance
                && (defect == 0.0
                    || correction.abs()
                        <= controls.relative_angular_tolerance * (stress - yield_stress))
            {
                return self.finish(
                    Some(stress.copysign(relative)),
                    walls[0],
                    Some(walls),
                    None,
                    value,
                    controls,
                    iteration,
                    budget.used,
                    &mut cancelled,
                );
            }
            if defect < 0.0 {
                low = stress;
            } else {
                high = stress;
            }
            let proposed = stress - correction;
            let next = if proposed.is_finite() && proposed > low && proposed < high {
                proposed
            } else {
                0.5 * low + 0.5 * high
            };
            if next <= low || next >= high || next == stress {
                return Err(YieldStressJournalError::UnsupportedMotion);
            }
            stress = next;
        }
        Err(YieldStressJournalError::StressSolveBudget)
    }

    /// Signed torque ON THE FLUID at the outer wall resolves the static stress
    /// as well as flowing states. The inner wall speed anchors common rotation.
    pub fn at_outer_torque(
        &self,
        torque_nm: f64,
        inner_speed_rad_s: f64,
        controls: YieldStressJournalControls,
    ) -> Result<YieldStressJournalProfile<'a>> {
        self.torque_with_cancel(torque_nm, inner_speed_rad_s, controls, || false)
    }
    pub fn torque_with_cancel(
        &self,
        torque_nm: f64,
        inner_speed_rad_s: f64,
        controls: YieldStressJournalControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<YieldStressJournalProfile<'a>> {
        controls.admit()?;
        poll(&mut cancelled)?;
        if !torque_nm.is_finite() || !inner_speed_rad_s.is_finite() {
            return Err(YieldStressJournalError::Invalid(
                "nonfinite torque or wall speed",
            ));
        }
        let maximum_torque = product_ratio(&[self.torque_factor, self.maximum_stress], &[])?;
        if torque_nm.abs() > maximum_torque {
            return Err(YieldStressJournalError::UnsupportedMotion);
        }
        let mut stress =
            product_ratio(&[torque_nm.abs()], &[self.torque_factor])?.min(self.maximum_stress);
        if torque_nm.abs() <= self.static_outer_torque_bounds_nm()[1] {
            stress = stress.min(self.source.parameters().yield_stress.value());
        }
        let mut budget = Budget {
            used: 0,
            maximum: controls.maximum_quadrature_nodes,
        };
        let value = self.evaluate(
            stress,
            self.geometry.outer_radius_m,
            controls,
            &mut budget,
            &mut cancelled,
        )?;
        self.finish(
            Some(stress.copysign(torque_nm)),
            inner_speed_rad_s,
            None,
            Some(torque_nm),
            value,
            controls,
            0,
            budget.used,
            &mut cancelled,
        )
    }

    fn evaluate(
        &self,
        stress: f64,
        outer: f64,
        controls: YieldStressJournalControls,
        budget: &mut Budget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        poll(cancelled)?;
        if stress > self.maximum_stress {
            return Err(YieldStressJournalError::UnsupportedMotion);
        }
        let source_rate = self
            .source
            .shear_rate_at_stress(Stress::new(stress), self.conditions)?
            .value();
        // The declared ceiling was evaluated through the same original law at
        // construction. Preserve that exact endpoint across forward/inverse
        // rounding; interior stresses remain bounded by its monotone source.
        let rate = if stress == self.maximum_stress {
            self.maximum_rate
        } else {
            source_rate.min(self.maximum_rate)
        };
        if rate == 0.0 {
            return Ok(Evaluation {
                speed: 0.0,
                rates: [0.0; 2],
                rate_drop_fraction: 0.0,
                change: 0.0,
            });
        }
        let yield_stress = self.source.parameters().yield_stress.value();
        let excess = stress - yield_stress;
        let drop = product_ratio(
            &[stress, radial_drop(self.geometry.inner_radius_m, outer)],
            &[],
        )?;
        let n = self.source.parameters().flow_index.value();
        let log_rate_ratio = if yield_stress == 0.0 {
            let a = self.geometry.inner_radius_m;
            let logarithm = if (outer - a) < a {
                log1p((outer - a) / a)
            } else {
                det::ln(outer) - det::ln(a)
            };
            -2.0 * logarithm / n
        } else if drop >= excess {
            f64::NEG_INFINITY
        } else {
            log1p(-drop / excess) / n
        };
        let rate_drop_fraction = -det::expm1(log_rate_ratio);
        let outer_rate = if log_rate_ratio == f64::NEG_INFINITY {
            0.0
        } else {
            product_ratio(&[rate, det::exp(log_rate_ratio)], &[])?
        };
        if yield_stress == 0.0 {
            return Ok(Evaluation {
                speed: product_ratio(&[0.5, n, rate, rate_drop_fraction], &[])?,
                rates: [rate, outer_rate],
                rate_drop_fraction,
                change: 0.0,
            });
        }
        let p = (5.0 / (n + 1.0)).max(1.0);
        let width = -det::expm1(log_rate_ratio / p);
        let prefactor = product_ratio(&[0.5, n, p, rate, width], &[])?;
        let e = excess / stress;
        let y = yield_stress / stress;
        let mut node = |t: f64| -> Result<f64> {
            poll(cancelled)?;
            if budget.used == budget.maximum {
                return Err(YieldStressJournalError::QuadratureBudget);
            }
            budget.used += 1;
            let z = (-width).mul_add(t, 1.0);
            if z == 0.0 {
                return Ok(0.0);
            }
            let numerator = e * det::pow(z, p * n);
            finite(det::pow(z, p - 1.0) * (numerator / (y + numerator)))
        };
        let (integral, change) = integrate(&mut node, controls)?;
        Ok(Evaluation {
            speed: product_ratio(&[prefactor, integral], &[])?,
            rates: [rate, outer_rate],
            rate_drop_fraction,
            change,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        stress: Option<f64>,
        inner_speed: f64,
        requested_walls: Option<[f64; 2]>,
        requested_torque: Option<f64>,
        value: Evaluation,
        controls: YieldStressJournalControls,
        iterations: u32,
        nodes: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<YieldStressJournalProfile<'a>> {
        let sign = stress.unwrap_or(0.0).signum();
        let relative = value.speed * sign;
        let walls = [inner_speed, finite(inner_speed + relative)?];
        let torque = stress
            .map(|s| product_ratio(&[self.torque_factor, s.abs()], &[]).map(|t| t.copysign(s)))
            .transpose()?;
        let heat = product_ratio(&[torque.unwrap_or(0.0).abs(), value.speed], &[])?;
        let requested_relative = requested_walls.map_or(relative, |w| w[1] - w[0]);
        let power = finite(requested_torque.or(torque).unwrap_or(0.0) * requested_relative)?;
        let yielded_outer = if value.rates[0] == 0.0 {
            self.geometry.inner_radius_m
        } else if value.rates[1] > 0.0 || self.source.parameters().yield_stress.value() == 0.0 {
            self.geometry.outer_radius_m
        } else {
            let s = stress.unwrap().abs();
            let y = self.source.parameters().yield_stress.value();
            let logarithm = if (s - y) / y < 1.0 {
                log1p((s - y) / y)
            } else {
                det::ln(s) - det::ln(y)
            };
            product_ratio(
                &[self.geometry.inner_radius_m, det::exp(0.5 * logarithm)],
                &[],
            )?
            .min(self.geometry.outer_radius_m)
        };
        let regime = if value.rates[0] == 0.0 {
            JournalYieldRegime::Unyielded
        } else if yielded_outer < self.geometry.outer_radius_m {
            JournalYieldRegime::PartiallyYielded
        } else {
            JournalYieldRegime::FullyYielded
        };
        let resistance = if value.speed == 0.0 {
            None
        } else {
            Some(product_ratio(
                &[2.0, self.torque_factor, stress.unwrap().abs()],
                &[value.rates[0], value.rate_drop_fraction],
            )?)
        };
        poll(cancelled)?;
        Ok(YieldStressJournalProfile {
            journal: *self,
            controls,
            inner_stress: stress,
            flow: YieldStressJournalFlow {
                regime,
                relative_angular_velocity_rad_s: relative,
                accepted_wall_angular_velocities_rad_s: walls,
                requested_wall_angular_velocities_rad_s: requested_walls,
                angular_velocity_residual_rad_s: finite(relative - requested_relative)?,
                requested_outer_torque_nm: requested_torque,
                outer_torque_residual_nm: requested_torque.map(|t| torque.unwrap() - t),
                wall_torques_nm: torque.map(|t| [-t, t]),
                static_outer_torque_bounds_nm: self.static_outer_torque_bounds_nm(),
                yielded_outer_radius_m: yielded_outer,
                wall_shear_rates_per_s: value.rates.map(|r| sign * r),
                wall_power_w: power,
                dissipated_power_w: heat,
                work_residual_w: finite(power - heat)?,
                differential_resistance_nm_s: resistance,
                quadrature_relative_change: value.change,
                stress_iterations: iterations,
                quadrature_nodes: nodes,
            },
        })
    }
}

impl<'a> YieldStressJournalProfile<'a> {
    pub const fn journal(&self) -> &YieldStressJournal<'a> {
        &self.journal
    }
    pub const fn flow(&self) -> &YieldStressJournalFlow {
        &self.flow
    }
    pub fn sample(&self, radius_m: f64) -> Result<YieldStressJournalSample> {
        self.sample_with_cancel(radius_m, || false)
    }
    pub fn sample_with_cancel(
        &self,
        radius: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<YieldStressJournalSample> {
        poll(&mut cancelled)?;
        let j = &self.journal;
        if !radius.is_finite()
            || radius < j.geometry.inner_radius_m
            || radius > j.geometry.outer_radius_m
        {
            return Err(YieldStressJournalError::Invalid(
                "sample is outside the annular gap",
            ));
        }
        let mut rate = 0.0;
        let mut speed = self.flow.accepted_wall_angular_velocities_rad_s[0];
        let stress = if let Some(inner) = self.inner_stress {
            let drop = radial_drop(j.geometry.inner_radius_m, radius);
            let signed = if j.geometry.inner_radius_m / radius < 0.5 {
                inner * (j.geometry.inner_radius_m / radius).powi(2)
            } else {
                (-inner).mul_add(drop, inner)
            };
            if self.flow.wall_shear_rates_per_s[0] != 0.0 {
                if radius == j.geometry.inner_radius_m {
                    rate = self.flow.wall_shear_rates_per_s[0];
                } else if radius == j.geometry.outer_radius_m {
                    rate = self.flow.wall_shear_rates_per_s[1];
                    speed = self.flow.accepted_wall_angular_velocities_rad_s[1];
                } else if signed.abs() <= j.source.parameters().yield_stress.value() {
                    // A rigid plug has the outer wall's angular speed. Avoid
                    // repeating any constitutive inversion or quadrature there.
                    speed = self.flow.accepted_wall_angular_velocities_rad_s[1];
                } else {
                    let mut budget = Budget {
                        used: 0,
                        maximum: self.controls.maximum_quadrature_nodes,
                    };
                    let partial = j.evaluate(
                        inner.abs(),
                        radius,
                        self.controls,
                        &mut budget,
                        &mut cancelled,
                    )?;
                    rate = partial.rates[1].copysign(inner);
                    speed = if rate == 0.0 {
                        self.flow.accepted_wall_angular_velocities_rad_s[1]
                    } else {
                        finite(speed + partial.speed.copysign(inner))?
                    };
                }
            }
            Some(signed)
        } else {
            None
        };
        let heat = product_ratio(&[stress.unwrap_or(0.0).abs(), rate.abs()], &[])?;
        poll(&mut cancelled)?;
        Ok(YieldStressJournalSample {
            angular_velocity_rad_s: speed,
            azimuthal_velocity_m_s: finite(radius * speed)?,
            shear_rate_per_s: rate,
            radial_shear_stress_pa: stress,
            dissipation_w_m3: heat,
            yielded: rate != 0.0,
        })
    }
}

fn integrate(
    node: &mut impl FnMut(f64) -> Result<f64>,
    controls: YieldStressJournalControls,
) -> Result<(f64, f64)> {
    let mut trapezoid = 0.5 * node(0.0)? + 0.5 * node(1.0)?;
    let mut previous: Option<f64> = None;
    let mut agreements = 0;
    for level in 1..=controls.maximum_quadrature_refinements {
        let panels = 1_u32 << level;
        let h = 1.0 / f64::from(panels);
        let (mut added, mut compensation) = (0.0, 0.0);
        for i in (1..panels).step_by(2) {
            let adjusted = node(f64::from(i) * h)? - compensation;
            let next = added + adjusted;
            compensation = (next - added) - adjusted;
            added = next;
        }
        let value = finite((trapezoid + 4.0 * h * added) / 3.0)?;
        if let Some(old) = previous {
            let change = if value == old {
                0.0
            } else {
                (value - old).abs() / value
            };
            if change <= controls.relative_quadrature_tolerance {
                agreements += 1;
                if agreements == 2 {
                    return Ok((value, change));
                }
            } else {
                agreements = 0;
            }
        }
        previous = Some(value);
        trapezoid = 0.5 * trapezoid + h * added;
    }
    Err(YieldStressJournalError::QuadratureResolution)
}
