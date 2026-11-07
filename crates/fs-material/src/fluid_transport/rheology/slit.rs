//! No-slip planar Couette--Poiseuille flow of the original thermal Carreau law.
//!
//! Momentum fixes `tau(y) = c + delta_p*(h/2-y)/L`. The ONE stress offset c
//! enforces the relative wall speed through the canonical constitutive inverse.
//! Sliding and pressure therefore share a velocity profile, wall tractions and
//! heat. A frozen apparent viscosity, or a sum of separately solved nonlinear
//! Couette and Poiseuille profiles, would give a different physical problem.
//!
//! This is an instantaneous, fully developed, incompressible planar reduction.
//! Width must be large enough that side walls can be neglected. Temperature is
//! uniform and frozen; the caller owns thin-gap/laminar admission, pressure
//! limits, entrance effects, material provenance and thermal time integration.
//! No cavitation, slip, yield surface, elastic memory or thermal feedback is
//! inferred. A declared shear ceiling is not a turbulence certificate.

use super::pipe::{CarreauPipeControls, PreparedCarreau, RateState, integrate, product_ratio};
use super::{MaterialError, ThermalCarreauYasuda, parameters, positive, state_error};

#[cfg(test)]
mod tests;

/// Explicit numerical work and residual budgets; not source uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlitControls {
    /// Original shared constitutive inverse and positive quadrature controls.
    pub local: CarreauPipeControls,
    /// Relative no-slip wall-speed residual AND local stress correction.
    pub relative_wall_velocity_tolerance: f64,
    /// At most this many offset corrections, from one through 128.
    pub maximum_wall_stress_iterations: u32,
}

impl Default for CarreauSlitControls {
    fn default() -> Self {
        Self {
            local: CarreauPipeControls::default(),
            relative_wall_velocity_tolerance: 1e-9,
            maximum_wall_stress_iterations: 96,
        }
    }
}

impl CarreauSlitControls {
    fn admit(self) -> Result<(), MaterialError> {
        self.local.admit()?;
        if !positive(self.relative_wall_velocity_tolerance)
            || self.relative_wall_velocity_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_wall_stress_iterations)
        {
            return Err(parameters("invalid Carreau slit wall-speed solve budgets"));
        }
        Ok(())
    }
}

/// Original SI law, uniform temperature, physical geometry and rate envelope.
/// Both zero shear and the entire interval to the ceiling must be supported.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarCarreauSlit {
    prepared: PreparedCarreau,
    width_m: f64,
    gap_m: f64,
    length_m: f64,
}

/// Physical response of one accepted profile. Powers are instantaneous [W].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlitFlow {
    pub pressure_drop_pa: f64,
    /// Original no-slip data, in lower/upper order [m/s].
    pub requested_wall_velocities_m_s: [f64; 2],
    /// Same mean translation, with the solved relative speed. This explicit
    /// constitutive residual is never silently substituted into requested work.
    pub accepted_wall_velocities_m_s: [f64; 2],
    /// Accepted minus requested upper-minus-lower speed [m/s].
    pub wall_velocity_residual_m_s: f64,
    pub volumetric_flow_m3_s: f64,
    /// Excess over mean-wall translation; kept even under large entrainment.
    pub pressure_flow_m3_s: f64,
    pub entrained_flow_m3_s: f64,
    /// Resisting lower/upper tractions [Pa], conjugate to the wall velocities.
    pub wall_tractions_pa: [f64; 2],
    pub pressure_power_w: f64,
    /// Work at the REQUESTED wall velocities; either input port may be negative.
    pub wall_power_w: f64,
    /// Positive pressure-profile heat plus positive relative-sliding heat.
    pub dissipated_power_w: f64,
    /// Requested total input minus actual-profile heat, evaluated without
    /// cancelling the potentially large common-translation power exchange.
    pub work_residual_w: f64,
    /// dQ/d(delta_p) at fixed wall velocities [m^3/(s Pa)], strictly positive.
    pub pressure_conductance_m3_s_pa: f64,
    /// dQ/d(Uupper-Ulower) at fixed mean wall speed and pressure [m^2].
    /// The derivative of relative resisting force with respect to pressure
    /// is its negative, from the same profile and momentum balance.
    pub flow_relative_wall_velocity_derivative_m2: f64,
    /// d(width*length*c)/d(Uupper-Ulower) at fixed pressure [N s/m].
    pub relative_wall_resistance_n_s_m: f64,
    /// Largest observed relative refinement change; an estimate, not a bound.
    pub quadrature_relative_change: f64,
    pub wall_stress_iterations: u32,
}

/// Immutable accepted profile, retaining its original source and solve controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlitProfile {
    slit: LaminarCarreauSlit,
    controls: CarreauSlitControls,
    flow: CarreauSlitFlow,
    offset_pa: f64,
    half_pressure_stress_pa: f64,
}

/// Prescribed net flow and moving walls resolved through one accepted profile.
/// Pressure, tractions and heat always belong to that profile; requested flow
/// is never substituted into its physical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlitFlowSolution {
    profile: CarreauSlitProfile,
    requested_volumetric_flow_m3_s: f64,
    flow_residual_m3_s: f64,
    flow_iterations: u32,
}

impl CarreauSlitFlowSolution {
    pub const fn profile(&self) -> &CarreauSlitProfile {
        &self.profile
    }
    pub const fn requested_volumetric_flow_m3_s(&self) -> f64 {
        self.requested_volumetric_flow_m3_s
    }
    /// Accepted minus requested flow [m³/s], evaluated after subtracting common
    /// wall entrainment so a large translation does not hide the residual.
    pub const fn flow_residual_m3_s(&self) -> f64 {
        self.flow_residual_m3_s
    }
    pub const fn flow_iterations(&self) -> u32 {
        self.flow_iterations
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlitSample {
    pub axial_velocity_m_s: f64,
    pub shear_stress_pa: f64,
    pub velocity_gradient_per_s: f64,
    pub dynamic_viscosity_pa_s: f64,
    pub dissipation_w_m3: f64,
    pub quadrature_relative_change: f64,
}

impl LaminarCarreauSlit {
    pub fn new(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        width_m: f64,
        gap_m: f64,
        length_m: f64,
        maximum_laminar_shear_rate_per_s: f64,
    ) -> Result<Self, MaterialError> {
        if ![width_m, gap_m, length_m].into_iter().all(positive) {
            return Err(parameters(
                "Carreau slit geometry must be positive and finite",
            ));
        }
        let prepared = PreparedCarreau::new(law, temperature_k, maximum_laminar_shear_rate_per_s)?;
        if !positive(prepared.ceiling.stress()) {
            return Err(parameters(
                "Carreau slit needs a representable supported stress ceiling",
            ));
        }
        Ok(Self {
            prepared,
            width_m,
            gap_m,
            length_m,
        })
    }

    pub const fn law(&self) -> &ThermalCarreauYasuda {
        &self.prepared.law
    }
    pub const fn temperature_k(&self) -> f64 {
        self.prepared.temperature.temperature_k()
    }
    pub const fn width_m(&self) -> f64 {
        self.width_m
    }
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    /// Rebind the actual positive hydraulic aperture while retaining the
    /// original prepared temperature, law, support, width and passage length.
    /// Reynolds pressure iterations need not repeat thermal preparation.
    pub fn with_gap(&self, gap_m: f64) -> Result<Self, MaterialError> {
        if !positive(gap_m) {
            return Err(parameters("Carreau slit gap must be positive and finite"));
        }
        Ok(Self { gap_m, ..*self })
    }
    /// Original thermally shifted zero-shear viscosity [Pa s], useful for a
    /// source-consistent Newtonian limit or pressure-solver preconditioner.
    pub const fn zero_shear_viscosity_pa_s(&self) -> f64 {
        self.prepared.rest.viscosity
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    pub const fn maximum_laminar_shear_rate_per_s(&self) -> f64 {
        self.prepared.ceiling.rate
    }

    /// Solve simultaneous signed pressure drop and lower/upper wall velocities.
    /// There is no heap allocation, hidden material extension or mutable history.
    pub fn profile_at_pressure_and_walls(
        &self,
        pressure_drop_pa: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauSlitControls,
    ) -> Result<CarreauSlitProfile, MaterialError> {
        self.profile_with_check(pressure_drop_pa, wall_velocities_m_s, controls, || Ok(()))
    }

    /// Resolve the signed pressure needed for prescribed net volume flow and
    /// lower/upper wall velocities, including blockage or reversed pumping.
    ///
    /// Newton corrections use the accepted profile's positive pressure
    /// conductance. A bracket in physical wall stress admits both wall states;
    /// sliding and pressure cannot independently spend the same rate envelope.
    /// Constant-viscosity laws take their exact Couette--Poiseuille branch.
    /// All profiles retain the original coupled constitutive law and heat.
    ///
    /// `relative_flow_tolerance` gates the flow residual relative to the flow
    /// in excess of wall entrainment AND the relative pressure correction.
    /// No pressure is needed when that excess is exactly zero. Wall residuals
    /// and quadrature estimates retain their independent `controls` budgets.
    /// `maximum_flow_iterations` is in 1..=128 and bounds nonzero-pressure
    /// trials, in addition to one allocation-free pure-Couette evaluation.
    pub fn profile_at_flow_and_walls(
        &self,
        volumetric_flow_m3_s: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauSlitControls,
        relative_flow_tolerance: f64,
        maximum_flow_iterations: u32,
    ) -> Result<CarreauSlitFlowSolution, MaterialError> {
        self.profile_at_flow_and_walls_with_check(
            volumetric_flow_m3_s,
            wall_velocities_m_s,
            controls,
            relative_flow_tolerance,
            maximum_flow_iterations,
            || Ok(()),
        )
    }

    /// Same inverse with caller-owned checks throughout every local profile,
    /// source-span query and pressure correction, including before publication.
    /// Unsupported forcing, exhausted work and cancellation publish no result.
    pub fn profile_at_flow_and_walls_with_check(
        &self,
        volumetric_flow_m3_s: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauSlitControls,
        relative_flow_tolerance: f64,
        maximum_flow_iterations: u32,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauSlitFlowSolution, MaterialError> {
        controls.admit()?;
        if !volumetric_flow_m3_s.is_finite()
            || !positive(relative_flow_tolerance)
            || relative_flow_tolerance >= 1.0
            || !(1..=128).contains(&maximum_flow_iterations)
        {
            return Err(parameters(
                "invalid prescribed Carreau slit flow or inverse controls",
            ));
        }
        let rest = self.profile_with_check(0.0, wall_velocities_m_s, controls, &mut check)?;
        let excess = volumetric_flow_m3_s - rest.flow.entrained_flow_m3_s;
        if !excess.is_finite() {
            return Err(state_error(
                "prescribed Carreau slit flow relative to wall entrainment is unrepresentable",
            ));
        }
        let accepted = |profile: CarreauSlitProfile, flow_iterations| CarreauSlitFlowSolution {
            flow_residual_m3_s: profile.flow.pressure_flow_m3_s - excess,
            profile,
            requested_volumetric_flow_m3_s: volumetric_flow_m3_s,
            flow_iterations,
        };
        if excess == 0.0 {
            check()?;
            return Ok(accepted(rest, 0));
        }
        let target = excess.abs();
        if self.prepared.constant() {
            let pressure = product_ratio(&[target], &[rest.flow.pressure_conductance_m3_s_pa])?
                .copysign(excess);
            let profile =
                self.profile_with_check(pressure, wall_velocities_m_s, controls, &mut check)?;
            let difference = profile.flow.pressure_flow_m3_s - excess;
            if (difference / target).abs() > relative_flow_tolerance
                || (difference / profile.flow.pressure_flow_m3_s).abs() > relative_flow_tolerance
            {
                return Err(state_error(
                    "Carreau slit exact flow inverse cannot resolve the requested tolerance",
                ));
            }
            check()?;
            return Ok(accepted(profile, 1));
        }
        let ceiling = self.prepared.ceiling.stress();
        let uniform_stress = rest.offset_pa.abs();
        let relative_rate = (wall_velocities_m_s[1] - wall_velocities_m_s[0]).abs() / self.gap_m;
        if uniform_stress >= ceiling {
            return Err(state_error(
                "prescribed Carreau slit flow leaves no supported pressure-stress span",
            ));
        }
        let mut low = 0.0;
        let mut high = ceiling;
        // This is only a Newton predictor. If its linear extrapolation is
        // unrepresentable, the finite original stress ceiling still brackets
        // every supported nonlinear solution.
        let mut stress = product_ratio(
            &[target, self.gap_m],
            &[2.0, self.length_m, rest.flow.pressure_conductance_m3_s_pa],
        )
        .map_or(high, |value| value.min(high));
        for iterations in 1..=maximum_flow_iterations {
            check()?;
            // Concavity of the admitted shear-thinning stress law makes the
            // Couette offset an upper bound on the actual offset. Only trials
            // beyond this inexpensive safe span need the full source-limit
            // query. At the limiting offset, one wall reaches the ceiling.
            let supported = if relative_rate == 0.0 || stress <= ceiling - uniform_stress {
                true
            } else if stress == ceiling {
                false
            } else {
                let maximum = self
                    .relative_rate(ceiling - stress, stress, controls.local, &mut check)?
                    .0;
                relative_rate <= maximum
            };
            let candidate = if supported {
                let pressure =
                    product_ratio(&[2.0, self.length_m, stress], &[self.gap_m])?.copysign(excess);
                let profile =
                    self.profile_with_check(pressure, wall_velocities_m_s, controls, &mut check)?;
                let difference = profile.flow.pressure_flow_m3_s.abs() - target;
                // Express the correction in the same bounded stress coordinate
                // as the bracket; no finite-difference material calls are used.
                let correction = if difference == 0.0 {
                    0.0
                } else {
                    product_ratio(
                        &[difference.abs(), self.gap_m],
                        &[
                            2.0,
                            self.length_m,
                            profile.flow.pressure_conductance_m3_s_pa,
                        ],
                    )?
                    .copysign(difference)
                };
                if (difference / target).abs() <= relative_flow_tolerance
                    && (correction / stress).abs() <= relative_flow_tolerance
                {
                    check()?;
                    return Ok(accepted(profile, iterations));
                }
                if difference < 0.0 {
                    if stress == ceiling {
                        return Err(state_error(
                            "prescribed flow exceeds the Carreau slit shear-rate envelope",
                        ));
                    }
                    low = stress;
                } else {
                    high = stress;
                }
                Some(stress - correction)
            } else {
                high = stress;
                None
            };
            let next = candidate
                .filter(|value| value.is_finite() && *value > low && *value < high)
                .unwrap_or_else(|| low + 0.5 * (high - low));
            if next == stress || next <= low || next >= high {
                return Err(state_error(
                    "Carreau slit flow inverse cannot resolve the requested flow inside its stress envelope",
                ));
            }
            stress = next;
        }
        Err(state_error(
            "Carreau slit flow inverse exhausted its iteration budget",
        ))
    }

    /// Same solve with bounded caller-owned checks between quadrature samples
    /// and offset corrections. At most one capped constitutive inverse occurs
    /// between checks. Refusal publishes no partially accepted profile.
    pub fn profile_with_check(
        &self,
        pressure_drop_pa: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauSlitControls,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauSlitProfile, MaterialError> {
        controls.admit()?;
        check()?;
        let [lower, upper] = wall_velocities_m_s;
        let relative = upper - lower;
        let mean = if lower.abs() <= upper.abs() {
            0.5_f64.mul_add(relative, lower)
        } else {
            (-0.5_f64).mul_add(relative, upper)
        };
        if [pressure_drop_pa, lower, upper, relative, mean]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(parameters(
                "Carreau slit pressure or wall velocities are nonfinite",
            ));
        }
        let half_stress = if pressure_drop_pa == 0.0 {
            0.0
        } else {
            product_ratio(&[pressure_drop_pa.abs(), self.gap_m], &[2.0, self.length_m])?
                .copysign(pressure_drop_pa)
        };
        let target_rate = relative.abs() / self.gap_m;
        if !target_rate.is_finite()
            || relative != 0.0 && target_rate == 0.0
            || target_rate > self.prepared.ceiling.rate
            || half_stress.abs() > self.prepared.ceiling.stress()
        {
            return Err(state_error(
                "Carreau slit forcing exceeds its supported rate/stress envelope",
            ));
        }
        let uniform = self.prepared.sample_rate(target_rate)?;
        let (offset, actual_rate, iterations, offset_change) =
            if half_stress == 0.0 || self.prepared.constant() {
                (uniform.stress().copysign(relative), target_rate, 0, 0.0)
            } else if relative == 0.0 {
                (0.0, 0.0, 0, 0.0)
            } else {
                let (c, rate, iterations, change) = self.solve_offset(
                    half_stress.abs(),
                    target_rate,
                    uniform.stress(),
                    controls,
                    &mut check,
                )?;
                (c.copysign(relative), rate, iterations, change)
            };
        // Both actual wall states are admitted, including on exact Newtonian
        // and pure-Couette branches; their shear may exceed the mean rate.
        for stress in [offset + half_stress, offset - half_stress] {
            check()?;
            self.prepared
                .inverse(stress.abs(), self.prepared.ceiling, controls.local)?;
        }
        let (moments, quadrature_change, viscosity_scale) =
            if half_stress == 0.0 || self.prepared.constant() {
                (
                    Moments {
                        weight: 1.0,
                        mean: 0.5,
                        variance: 1.0 / 12.0,
                        pressure: 1.0 / 6.0,
                    },
                    0.0,
                    uniform.tangent,
                )
            } else {
                let scale = self.prepared.rest.viscosity;
                let (moments, change) =
                    self.moments(offset, half_stress, scale, controls.local, &mut check)?;
                (moments, change, scale)
            };
        let conductance = product_ratio(
            &[
                self.width_m,
                self.gap_m,
                self.gap_m,
                self.gap_m,
                moments.variance,
            ],
            &[self.length_m, viscosity_scale],
        )?;
        let resistance = product_ratio(
            &[self.width_m, self.length_m, viscosity_scale],
            &[self.gap_m, moments.weight],
        )?;
        let pressure_flow = if pressure_drop_pa == 0.0 {
            0.0
        } else {
            product_ratio(
                &[
                    self.width_m,
                    self.gap_m,
                    self.gap_m,
                    self.gap_m,
                    pressure_drop_pa.abs(),
                    moments.pressure,
                ],
                &[2.0, self.length_m, viscosity_scale],
            )?
            .copysign(pressure_drop_pa)
        };
        let entrained = self.width_m * self.gap_m * mean;
        let total_flow = entrained + pressure_flow;
        let actual_relative = (self.gap_m * actual_rate).copysign(relative);
        let residual = actual_relative - relative;
        let accepted_walls = [mean - 0.5 * actual_relative, mean + 0.5 * actual_relative];
        let area = self.width_m * self.length_m;
        let tractions = [-offset - half_stress, offset - half_stress];
        let pressure_heat = pressure_drop_pa * pressure_flow;
        let shear_heat = area * offset * actual_relative;
        let heat = pressure_heat + shear_heat;
        let pressure_power = pressure_drop_pa * total_flow;
        let wall_power = area * offset * relative - pressure_drop_pa * entrained;
        let work_residual = -area * offset * residual;
        let cross = self.width_m * self.gap_m * (0.5 - moments.mean);
        if [
            entrained,
            total_flow,
            actual_relative,
            residual,
            area,
            pressure_heat,
            shear_heat,
            heat,
            pressure_power,
            wall_power,
            work_residual,
            cross,
            accepted_walls[0],
            accepted_walls[1],
            tractions[0],
            tractions[1],
        ]
        .iter()
        .any(|x| !x.is_finite())
            || pressure_heat < 0.0
            || shear_heat < 0.0
        {
            return Err(state_error(
                "Carreau slit flow, traction or physical work is unrepresentable",
            ));
        }
        check()?;
        Ok(CarreauSlitProfile {
            slit: *self,
            controls,
            offset_pa: offset,
            half_pressure_stress_pa: half_stress,
            flow: CarreauSlitFlow {
                pressure_drop_pa,
                requested_wall_velocities_m_s: wall_velocities_m_s,
                accepted_wall_velocities_m_s: accepted_walls,
                wall_velocity_residual_m_s: residual,
                volumetric_flow_m3_s: total_flow,
                pressure_flow_m3_s: pressure_flow,
                entrained_flow_m3_s: entrained,
                wall_tractions_pa: tractions,
                pressure_power_w: pressure_power,
                wall_power_w: wall_power,
                dissipated_power_w: heat,
                work_residual_w: work_residual,
                pressure_conductance_m3_s_pa: conductance,
                flow_relative_wall_velocity_derivative_m2: cross,
                relative_wall_resistance_n_s_m: resistance,
                quadrature_relative_change: offset_change.max(quadrature_change),
                wall_stress_iterations: iterations,
            },
        })
    }

    // Symmetry removes cancellation in nearly pure pressure flow. The mean
    // signed rate is min(c/s,1) times the positive inverse-stress average on
    // [abs(c-s),c+s], never the difference of two nearly equal primitives.
    fn relative_rate(
        &self,
        c: f64,
        s: f64,
        controls: CarreauPipeControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<(f64, f64, f64), MaterialError> {
        if c == 0.0 {
            return Ok((0.0, 0.0, 0.0));
        }
        let lower = (c - s).abs();
        let width = 2.0 * c.min(s);
        check()?;
        let (high, _) = self
            .prepared
            .inverse(c + s, self.prepared.ceiling, controls)?;
        check()?;
        let (low, _) = self
            .prepared
            .inverse(lower, self.prepared.ceiling, controls)?;
        let split = self
            .prepared
            .transition_rate(high.rate)
            .map(|rate| {
                self.prepared
                    .sample_rate(rate)
                    .map(|p| (p.stress() - lower) / width)
            })
            .transpose()?;
        let (integral, change) = integrate(
            &mut |t| {
                check()?;
                if t == 0.0 {
                    return Ok(low.rate / high.rate);
                }
                if t == 1.0 {
                    return Ok(1.0);
                }
                let stress = lower + width * t;
                Ok(self
                    .prepared
                    .inverse(stress, self.prepared.ceiling, controls)?
                    .0
                    .rate
                    / high.rate)
            },
            split,
            controls,
        )?;
        let factor = if c < s { c / s } else { 1.0 };
        let rate = product_ratio(&[factor, high.rate, integral], &[])?;
        let slope = if c <= s {
            (0.5 * high.rate + 0.5 * low.rate) / s
        } else {
            (0.5 * (high.rate - low.rate)) / s
        };
        // This derivative proposes a root correction only. When neighbouring
        // endpoint rates round together, the local tangent supplies a positive
        // proposal; the independently integrated no-slip residual still gates it.
        let slope = if positive(slope) {
            slope
        } else {
            1.0 / high.tangent
        };
        if !positive(slope) {
            return Err(state_error(
                "Carreau slit offset derivative is unrepresentable",
            ));
        }
        Ok((rate, slope, change))
    }

    fn solve_offset(
        &self,
        s: f64,
        target: f64,
        seed: f64,
        controls: CarreauSlitControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<(f64, f64, u32, f64), MaterialError> {
        let mut low = 0.0;
        let mut high = self.prepared.ceiling.stress() - s;
        if !positive(high) {
            return Err(state_error(
                "Carreau slit has no supported wall-stress offset",
            ));
        }
        let maximum = self.relative_rate(high, s, controls.local, check)?.0;
        if target > maximum {
            return Err(state_error(
                "Carreau slit sliding exceeds the supported stress span",
            ));
        }
        let mut c = seed.min(high);
        let mut change = 0.0_f64;
        for iteration in 1..=controls.maximum_wall_stress_iterations {
            check()?;
            let (actual, derivative, local_change) =
                self.relative_rate(c, s, controls.local, check)?;
            change = change.max(local_change);
            let residual = actual - target;
            let correction = residual / derivative;
            if (residual / target).abs() <= controls.relative_wall_velocity_tolerance
                && (correction / c).abs() <= controls.relative_wall_velocity_tolerance
            {
                return Ok((c, actual, iteration, change));
            }
            if residual < 0.0 {
                low = c;
            } else {
                high = c;
            }
            let candidate = c - correction;
            let next = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                low + 0.5 * (high - low)
            };
            if next == c || next <= low || next >= high {
                return Err(state_error(
                    "Carreau slit offset solve lost floating-point resolution",
                ));
            }
            c = next;
        }
        Err(state_error(
            "Carreau slit offset solve exhausted its iteration budget",
        ))
    }

    fn moments(
        &self,
        c: f64,
        s: f64,
        scale: f64,
        controls: CarreauPipeControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<(Moments, f64), MaterialError> {
        let mut cuts = [0.0, 1.0, 1.0, 1.0, 1.0];
        let mut count = 2;
        let transition = self
            .prepared
            .transition_rate(self.prepared.ceiling.rate)
            .map(|g| self.prepared.sample_rate(g).map(RateState::stress))
            .transpose()?;
        for stress in [Some(0.0), transition, transition.map(|v| -v)]
            .into_iter()
            .flatten()
        {
            let t = 0.5 + 0.5 * ((c - stress) / s);
            if t > 0.0 && t < 1.0 {
                cuts[count] = t;
                count += 1;
            }
        }
        cuts[..count].sort_unstable_by(f64::total_cmp);
        let mut total = Moments::default();
        let mut change = 0.0_f64;
        for endpoints in cuts[..count].windows(2) {
            if endpoints[0] == endpoints[1] {
                continue;
            }
            let (part, local_change) = integrate_moments(
                &mut |t| {
                    check()?;
                    let stress = c + s * (1.0 - 2.0 * t);
                    let local = self
                        .prepared
                        .inverse(stress.abs(), self.prepared.ceiling, controls)?
                        .0;
                    let weight = scale / local.tangent;
                    if !positive(weight) {
                        return Err(state_error(
                            "Carreau slit inverse tangent is unrepresentable",
                        ));
                    }
                    Ok(Moments {
                        weight,
                        mean: t,
                        variance: 0.0,
                        pressure: t * (1.0 - t) * weight,
                    })
                },
                endpoints[0],
                endpoints[1],
                controls,
            )?;
            total = total.combine(part);
            change = change.max(local_change);
        }
        if !total.valid() {
            return Err(state_error(
                "Carreau slit positive mobility moments are unrepresentable",
            ));
        }
        Ok((total, change))
    }
}

// Positive weighted central moments avoid subtracting M2-M1^2/M0 when a
// strongly varying viscosity concentrates the inverse tangent in a thin band.
// Nested Simpson refinement reuses every previous sample, with four scalars
// per interval and no per-query allocation or constitutive resampling table.
#[derive(Clone, Copy, Default)]
struct Moments {
    weight: f64,
    mean: f64,
    variance: f64,
    pressure: f64,
}

impl Moments {
    fn scaled(self, factor: f64) -> Self {
        Self {
            weight: self.weight * factor,
            mean: self.mean,
            variance: self.variance * factor,
            pressure: self.pressure * factor,
        }
    }
    fn combine(self, other: Self) -> Self {
        if self.weight == 0.0 {
            return other;
        }
        if other.weight == 0.0 {
            return self;
        }
        let weight = self.weight + other.weight;
        let delta = other.mean - self.mean;
        Self {
            weight,
            mean: self.mean + delta * (other.weight / weight),
            variance: self.variance
                + other.variance
                + delta * delta * (self.weight / weight) * other.weight,
            pressure: self.pressure + other.pressure,
        }
    }
    fn valid(self) -> bool {
        [self.weight, self.variance, self.pressure]
            .into_iter()
            .all(positive)
            && self.mean.is_finite()
            && (0.0..=1.0).contains(&self.mean)
    }
    fn change(self, previous: Self) -> f64 {
        let mut change = 0.0_f64;
        for (value, old) in [
            (self.weight, previous.weight),
            (self.weight * self.mean, previous.weight * previous.mean),
            (self.variance, previous.variance),
            (self.pressure, previous.pressure),
        ] {
            change = change.max(((value - old) / value).abs());
        }
        change
    }
}

fn integrate_moments(
    sample: &mut impl FnMut(f64) -> Result<Moments, MaterialError>,
    lower: f64,
    upper: f64,
    controls: CarreauPipeControls,
) -> Result<(Moments, f64), MaterialError> {
    let width = upper - lower;
    let mut trapezoid = sample(lower)?
        .scaled(0.5)
        .combine(sample(upper)?.scaled(0.5));
    let mut previous = None;
    let mut agreements = 0;
    for level in 1..=controls.maximum_quadrature_refinements {
        let panels = 1_u32 << level;
        let h = 1.0 / f64::from(panels);
        let mut added = Moments::default();
        for i in (1..panels).step_by(2) {
            added = added.combine(sample(lower + width * f64::from(i) * h)?);
        }
        let simpson = trapezoid
            .scaled(1.0 / 3.0)
            .combine(added.scaled(4.0 * h / 3.0));
        if !simpson.valid() {
            return Err(state_error(
                "Carreau slit quadrature moments are unrepresentable",
            ));
        }
        if let Some(old) = previous {
            let change = simpson.change(old);
            if change <= controls.relative_quadrature_tolerance {
                agreements += 1;
                if agreements == 2 {
                    return Ok((simpson.scaled(width), change));
                }
            } else {
                agreements = 0;
            }
        }
        previous = Some(simpson);
        trapezoid = trapezoid.scaled(0.5).combine(added.scaled(h));
    }
    Err(state_error(
        "Carreau slit quadrature exhausted its refinement budget",
    ))
}

impl CarreauSlitProfile {
    pub const fn slit(&self) -> &LaminarCarreauSlit {
        &self.slit
    }
    pub const fn controls(&self) -> CarreauSlitControls {
        self.controls
    }
    pub const fn flow(&self) -> &CarreauSlitFlow {
        &self.flow
    }
    /// Solved midpoint shear stress c [Pa]. Retaining this small sliding
    /// component avoids subtracting nearly equal pressure-driven wall stresses.
    pub const fn stress_offset_pa(&self) -> f64 {
        self.offset_pa
    }

    /// Sample physical y in [0,h] from the lower wall. Endpoints use the
    /// explicitly accepted no-slip velocities; any requested-speed discrepancy
    /// remains in the flow report. The same source gives local stress and heat.
    pub fn sample(
        &self,
        distance_from_lower_wall_m: f64,
    ) -> Result<CarreauSlitSample, MaterialError> {
        self.sample_with_check(distance_from_lower_wall_m, || Ok(()))
    }

    /// Same physical sample with a caller-owned check before each capped
    /// constitutive inverse and before publishing the accepted sample.
    pub fn sample_with_check(
        &self,
        distance_from_lower_wall_m: f64,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauSlitSample, MaterialError> {
        check()?;
        let h = self.slit.gap_m;
        let y = distance_from_lower_wall_m;
        if !y.is_finite() || y < 0.0 || y > h {
            return Err(parameters(
                "Carreau slit sample must lie between its physical walls",
            ));
        }
        let t = y / h;
        let stress = self.offset_pa + self.half_pressure_stress_pa * (1.0 - 2.0 * t);
        let local = self
            .slit
            .prepared
            .inverse(
                stress.abs(),
                self.slit.prepared.ceiling,
                self.controls.local,
            )?
            .0;
        let rate = local.rate.copysign(stress);
        let (velocity, change) = if y == 0.0 {
            (self.flow.accepted_wall_velocities_m_s[0], 0.0)
        } else if y == h {
            (self.flow.accepted_wall_velocities_m_s[1], 0.0)
        } else if self.slit.prepared.constant() || self.half_pressure_stress_pa == 0.0 {
            let lower = self.flow.accepted_wall_velocities_m_s[0];
            let delta = self.flow.accepted_wall_velocities_m_s[1] - lower;
            let pressure = self.flow.pressure_drop_pa * y * (h - y)
                / (2.0 * self.slit.length_m * local.viscosity);
            (lower + t * delta + pressure, 0.0)
        } else {
            let (lower, upper, wall, sign) = if t <= 0.5 {
                (0.0, t, self.flow.accepted_wall_velocities_m_s[0], 1.0)
            } else {
                (t, 1.0, self.flow.accepted_wall_velocities_m_s[1], -1.0)
            };
            let zero = 0.5 + 0.5 * self.offset_pa / self.half_pressure_stress_pa;
            let split = if zero > lower && zero < upper {
                zero
            } else {
                upper
            };
            let mut sum = 0.0;
            let mut change = 0.0_f64;
            for (a, b) in [(lower, split), (split, upper)] {
                if a == b {
                    continue;
                }
                let scale = self.slit.prepared.ceiling.rate;
                let (integral, local_change) = integrate(
                    &mut |q| {
                        check()?;
                        let x = a + (b - a) * q;
                        let stress =
                            self.offset_pa + self.half_pressure_stress_pa * (1.0 - 2.0 * x);
                        Ok(self
                            .slit
                            .prepared
                            .inverse(
                                stress.abs(),
                                self.slit.prepared.ceiling,
                                self.controls.local,
                            )?
                            .0
                            .rate
                            / scale)
                    },
                    None,
                    self.controls.local,
                )?;
                let midpoint_stress = self.offset_pa + self.half_pressure_stress_pa * (1.0 - a - b);
                sum += ((b - a) * scale * integral).copysign(midpoint_stress);
                change = change.max(local_change);
            }
            (wall + sign * h * sum, change)
        };
        let heat = stress * rate;
        if !velocity.is_finite() || !heat.is_finite() || heat < 0.0 {
            return Err(state_error(
                "Carreau slit profile velocity or local heat is unrepresentable",
            ));
        }
        check()?;
        Ok(CarreauSlitSample {
            axial_velocity_m_s: velocity,
            shear_stress_pa: stress,
            velocity_gradient_per_s: rate,
            dynamic_viscosity_pa_s: local.viscosity,
            dissipation_w_m3: heat,
            quadrature_relative_change: change,
        })
    }
}
