//! Fully developed circular-pipe flow of the original thermal Carreau-Yasuda law.
//!
//! Radial momentum fixes `tau(r) = tau_w r/R`. The constitutive inverse and
//! positive quadrature resolve the varying viscosity; a single apparent
//! viscosity is never substituted into Poiseuille's law. Temperature is frozen,
//! and the caller supplies a laminar wall-rate ceiling. No transition, entrance,
//! slip, elastic memory, pressure-dependent material law or thermal feedback is
//! inferred. Reported heat is available to a thermal owner, not deposited here.
//! A declared rate ceiling alone is not a Reynolds-number or turbulence certificate.

use super::{MaterialError, TemperatureShiftState, ThermalCarreauYasuda, ViscosityShiftMode};
use super::{parameters, positive, state_error};
use fs_math::det;

/// Radially resolved nonlinear viscous startup, relaxation and reversal.
pub mod transient;

/// The original steady pipe as a pressure/flow/heat Component.
#[cfg(feature = "component")]
pub mod component;

#[cfg(test)]
mod tests;

/// Numerical tolerances and hard work bounds, separate from source uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauPipeControls {
    /// Relative stress residual AND local Newton rate correction required by
    /// each monotone constitutive inverse. Must lie strictly between zero and one.
    pub relative_constitutive_tolerance: f64,
    /// Maximum iterations per constitutive inverse, from 1 through 128.
    pub maximum_constitutive_iterations: u32,
    /// Two consecutive Simpson refinements must change each positive integral
    /// by no more than this fraction. This is an estimate, not a rigorous bound.
    pub relative_quadrature_tolerance: f64,
    /// At most 2^N panels on each of at most two source-transition intervals;
    /// N must be 1..=20. Fewer than three refinements cannot establish agreement.
    pub maximum_quadrature_refinements: u32,
}

impl Default for CarreauPipeControls {
    fn default() -> Self {
        Self {
            relative_constitutive_tolerance: 1e-11,
            maximum_constitutive_iterations: 96,
            relative_quadrature_tolerance: 1e-8,
            maximum_quadrature_refinements: 16,
        }
    }
}

impl CarreauPipeControls {
    pub(super) fn admit(self) -> Result<(), MaterialError> {
        if ![
            self.relative_constitutive_tolerance,
            self.relative_quadrature_tolerance,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0 && v < 1.0)
            || !(1..=128).contains(&self.maximum_constitutive_iterations)
            || !(1..=20).contains(&self.maximum_quadrature_refinements)
        {
            return Err(parameters(
                "invalid Carreau pipe tolerances or finite work budgets",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct RateState {
    pub(super) rate: f64,
    pub(super) viscosity: f64,
    pub(super) tangent: f64,
}

impl RateState {
    pub(super) fn stress(self) -> f64 {
        self.rate * self.viscosity
    }
}

/// Canonical immutable constitutive preparation shared by pipes, slits and
/// resolved transient films. Geometry never enters this source-law owner.
/// The original temperature shift is evaluated once; local queries retain the
/// original reference curve, positive stress tangent and complete rate support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedCarreau {
    pub(super) law: ThermalCarreauYasuda,
    pub(super) temperature: TemperatureShiftState,
    pub(super) rest: RateState,
    pub(super) ceiling: RateState,
}

impl PreparedCarreau {
    pub fn new(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        maximum_rate_per_s: f64,
    ) -> Result<Self, MaterialError> {
        if !positive(maximum_rate_per_s) {
            return Err(parameters(
                "Carreau constitutive rate ceiling must be positive and finite",
            ));
        }
        let rest = law.evaluate(temperature_k, 0.0)?;
        let ceiling = law.evaluate(temperature_k, maximum_rate_per_s)?;
        Ok(Self {
            law,
            temperature: *rest.temperature(),
            rest: RateState {
                rate: 0.0,
                viscosity: rest.dynamic_viscosity_pa_s(),
                tangent: rest.stress_tangent_pa_s(),
            },
            ceiling: RateState {
                rate: maximum_rate_per_s,
                viscosity: ceiling.dynamic_viscosity_pa_s(),
                tangent: ceiling.stress_tangent_pa_s(),
            },
        })
    }

    /// Original material calibration; no apparent-viscosity fit is introduced.
    pub const fn law(&self) -> &ThermalCarreauYasuda {
        &self.law
    }

    /// Uniform held temperature of this constitutive preparation.
    pub const fn temperature_k(&self) -> f64 {
        self.temperature.temperature_k()
    }

    /// Physical shear-rate ceiling, unchanged by subsequent queries.
    pub const fn maximum_shear_rate_per_s(&self) -> f64 {
        self.ceiling.rate
    }

    /// Original zero-shear viscosity, for a diffusion scale or preconditioner.
    pub const fn zero_shear_viscosity_pa_s(&self) -> f64 {
        self.rest.viscosity
    }

    /// Apparent viscosity and d(mu*g)/dg [Pa s] at a nonnegative physical rate.
    /// Reuses the prepared shift and the original reference evaluator. The
    /// positive tangent remains defined when the scalar viscosity derivative
    /// is singular at rest; neither quantity is clipped or extended.
    pub fn viscosity_and_tangent(&self, rate_per_s: f64) -> Result<[f64; 2], MaterialError> {
        let sample = self.sample_rate(rate_per_s)?;
        Ok([sample.viscosity, sample.tangent])
    }

    /// Invert a signed shear stress through the SAME bounded constitutive
    /// solve used by pipe and slit profiles. Returns the signed rate and the
    /// number of local source evaluations; cached zero/ceiling queries use zero.
    /// The caller may reduce the iteration budget to its remaining work budget.
    pub fn signed_shear_rate_at_stress(
        &self,
        stress_pa: f64,
        controls: CarreauPipeControls,
    ) -> Result<(f64, u32), MaterialError> {
        controls.admit()?;
        if !stress_pa.is_finite() {
            return Err(state_error("Carreau shear stress must be finite"));
        }
        let (sample, iterations) = self.inverse(stress_pa.abs(), self.ceiling, controls)?;
        Ok((sample.rate.copysign(stress_pa), iterations))
    }

    pub(super) fn constant(&self) -> bool {
        let [zero, infinity, time, n, _] = self.law.curve().parameters();
        time == 0.0 || n == 1.0 || zero == infinity
    }

    // Reuse the canonical reference evaluator and frozen, original shift. Only
    // stress and its positive tangent are needed; no singular scalar mu' is used.
    pub(super) fn sample_rate(&self, rate: f64) -> Result<RateState, MaterialError> {
        if rate == 0.0 {
            return Ok(self.rest);
        }
        if rate == self.ceiling.rate {
            return Ok(self.ceiling);
        }
        if !positive(rate) || rate > self.ceiling.rate {
            return Err(state_error(
                "Carreau constitutive rate leaves the declared laminar envelope",
            ));
        }
        let reduced = match self.law.mode() {
            ViscosityShiftMode::Vertical => rate,
            ViscosityShiftMode::VerticalAndRate => self.temperature.factor() * rate,
        };
        if !positive(reduced) {
            return Err(state_error(
                "Carreau constitutive reduced rate is unrepresentable",
            ));
        }
        let reference = self.law.curve().sample(reduced)?;
        let viscosity = self.temperature.factor() * reference.viscosity;
        let tangent = self.temperature.factor() * reference.stress_tangent;
        if !positive(viscosity) || !positive(tangent) {
            return Err(state_error(
                "Carreau constitutive viscosity or stress tangent is unrepresentable",
            ));
        }
        Ok(RateState {
            rate,
            viscosity,
            tangent,
        })
    }

    pub(super) fn inverse(
        &self,
        target: f64,
        ceiling: RateState,
        controls: CarreauPipeControls,
    ) -> Result<(RateState, u32), MaterialError> {
        if target == 0.0 {
            return Ok((self.rest, 0));
        }
        if !positive(target) || target > ceiling.stress() {
            return Err(state_error(
                "Carreau constitutive stress exceeds the source or laminar wall-rate ceiling",
            ));
        }
        if target == ceiling.stress() {
            return Ok((ceiling, 0));
        }
        let mut low = 0.0;
        let mut high = ceiling.rate;
        let seed = target / self.rest.viscosity;
        // The smallest representable rate is a numerical initial iterate, never
        // a material cutoff. Its stress still has to meet the requested state.
        let mut rate = if seed > 0.0 { seed } else { f64::from_bits(1) };
        if rate >= high || !rate.is_finite() {
            rate = 0.5 * high;
        }
        for iterations in 1..=controls.maximum_constitutive_iterations {
            let sample = self.sample_rate(rate)?;
            let stress = sample.stress();
            let difference = stress - target;
            let correction = difference / sample.tangent;
            if positive(stress)
                && (difference / target).abs() <= controls.relative_constitutive_tolerance
                && (correction / rate).abs() <= controls.relative_constitutive_tolerance
            {
                return Ok((sample, iterations));
            }
            if stress < target {
                low = rate;
            } else {
                high = rate;
            }
            let candidate = rate - correction;
            let next = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                low + 0.5 * (high - low)
            };
            if next == rate || next <= low || next >= high {
                return Err(state_error(
                    "Carreau constitutive inverse lost floating-point resolution",
                ));
            }
            rate = next;
        }
        Err(state_error(
            "Carreau constitutive inverse exhausted its iteration budget",
        ))
    }

    // Split quadrature at the original curve's transition, when it is interior.
    // This is a numerical integration site, not an extension or a fitted law.
    pub(super) fn transition_rate(&self, wall_rate: f64) -> Option<f64> {
        let time = self.law.curve().parameters()[2];
        if self.constant() || time == 0.0 {
            return None;
        }
        let log_rate = -det::ln(time)
            - match self.law.mode() {
                ViscosityShiftMode::Vertical => 0.0,
                ViscosityShiftMode::VerticalAndRate => self.temperature.log_factor(),
            };
        let rate = det::exp(log_rate);
        (rate.is_finite() && rate > 0.0 && rate < wall_rate).then_some(rate)
    }
}

/// One original law, frozen temperature, no-slip geometry and explicit laminar
/// rate envelope. The entire rate interval from rest to the ceiling must be
/// source-supported: a pipe cannot invent an uncharacterized zero-shear core.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarCarreauPipe {
    prepared: PreparedCarreau,
    radius_m: f64,
    length_m: f64,
}

/// Integrated flow of the accepted constitutive state. Pressure and volume
/// flow have the same sign. Its pressure differs from the requested pressure
/// only by the explicitly reported constitutive inversion residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauPipeFlow {
    /// Original signed driving pressure, excluding elevation and entrance losses.
    pub requested_pressure_drop_pa: f64,
    /// Accepted signed pressure, obtained from `2 L mu(g_w) g_w / R`.
    pub pressure_drop_pa: f64,
    /// Accepted minus requested pressure [Pa]; never hidden in the heat account.
    pub pressure_residual_pa: f64,
    /// Signed integrated volume flow [m³/s].
    pub volumetric_flow_m3_s: f64,
    /// Signed driving wall stress [Pa]. Physical sigma_rz has the opposite sign.
    pub wall_stress_pa: f64,
    /// Signed wall shear-rate magnitude [1/s], oriented with the pressure drop.
    pub wall_shear_rate_per_s: f64,
    /// Nonnegative accepted pressure times volume flow [W].
    pub dissipated_power_w: f64,
    /// Last normalized change of the converged positive flow integral.
    pub quadrature_relative_change: f64,
    /// Actual constitutive iterations for the wall state.
    pub constitutive_iterations: u32,
}

/// Original pipe/law and accepted wall solution. Sampling performs bounded
/// constitutive solves and quadrature with these same controls; no external
/// material queries, source refits or mutable history are introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauPipeProfile {
    pipe: LaminarCarreauPipe,
    controls: CarreauPipeControls,
    flow: CarreauPipeFlow,
    wall: RateState,
    transition_stress_fraction: Option<f64>,
    normalized_flow_integral: f64,
}

/// A prescribed signed flow resolved through the original radial constitutive
/// law. The accepted profile owns pressure and heat; the flow residual is never
/// hidden by substituting the requested flow into either quantity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauPipeFlowSolution {
    profile: CarreauPipeProfile,
    requested_volumetric_flow_m3_s: f64,
    flow_residual_m3_s: f64,
    flow_iterations: u32,
}

impl CarreauPipeFlowSolution {
    pub const fn profile(&self) -> &CarreauPipeProfile {
        &self.profile
    }
    pub const fn requested_volumetric_flow_m3_s(&self) -> f64 {
        self.requested_volumetric_flow_m3_s
    }
    /// Accepted minus requested flow, in m³/s.
    pub const fn flow_residual_m3_s(&self) -> f64 {
        self.flow_residual_m3_s
    }
    pub const fn flow_iterations(&self) -> u32 {
        self.flow_iterations
    }
}

/// Local physical axial velocity, radial traction and irreversible heating.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauPipeSample {
    /// Signed axial velocity [m/s], exactly zero at the no-slip wall.
    pub axial_velocity_m_s: f64,
    /// Physical sigma_rz on the outward radial face [Pa].
    pub radial_shear_stress_pa: f64,
    /// Physical du_z/dr [1/s], with the same sign as sigma_rz.
    pub radial_velocity_gradient_per_s: f64,
    /// Original law's dynamic viscosity at this actual rate and temperature [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Nonnegative sigma_rz * du_z/dr [W/m³].
    pub dissipation_w_m3: f64,
    /// Last normalized change of the velocity integral; zero for exact limits.
    pub quadrature_relative_change: f64,
}

impl LaminarCarreauPipe {
    /// Admit finite positive radius, length and laminar physical wall-rate ceiling.
    /// Both the rest state and ceiling are checked against the original law's
    /// temperature and reduced-rate support. The thermal shift is prepared once.
    pub fn new(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        radius_m: f64,
        length_m: f64,
        maximum_laminar_wall_rate_per_s: f64,
    ) -> Result<Self, MaterialError> {
        if ![radius_m, length_m, maximum_laminar_wall_rate_per_s]
            .into_iter()
            .all(positive)
        {
            return Err(parameters(
                "Carreau pipe requires positive finite geometry and laminar rate ceiling",
            ));
        }
        Ok(Self {
            prepared: PreparedCarreau::new(law, temperature_k, maximum_laminar_wall_rate_per_s)?,
            radius_m,
            length_m,
        })
    }

    /// Exact authored law, including its original temperature-shift declaration.
    pub const fn law(&self) -> &ThermalCarreauYasuda {
        &self.prepared.law
    }
    pub const fn temperature_k(&self) -> f64 {
        self.prepared.temperature.temperature_k()
    }
    /// Re-evaluate the same declared source at an explicit uniform temperature.
    /// Geometry and the authored physical rate ceiling remain unchanged. The
    /// entire interval from rest to that ceiling must still lie inside the
    /// original reduced-rate support; a thermal shift never clips the ceiling.
    /// This is a frozen constitutive preparation, not a heat-evolution model.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<Self, MaterialError> {
        Ok(Self {
            prepared: PreparedCarreau::new(
                self.prepared.law,
                temperature_k,
                self.prepared.ceiling.rate,
            )?,
            ..*self
        })
    }
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Rebind only the positive physical length, retaining the prepared
    /// original temperature shift, constitutive support and radius. Moving
    /// menisci can reuse the same law without resolving its thermal state.
    pub fn with_length(&self, length_m: f64) -> Result<Self, MaterialError> {
        if !positive(length_m) {
            return Err(parameters(
                "Carreau pipe length must be positive and finite",
            ));
        }
        Ok(Self { length_m, ..*self })
    }
    pub const fn maximum_laminar_wall_rate_per_s(&self) -> f64 {
        self.prepared.ceiling.rate
    }

    /// Resolve a signed frictional pressure into a complete frozen radial state.
    ///
    /// With t=g/g_w and tau=mu(g)g, the volume-flow integral is
    /// `Q = pi R^3 g_w integral_0^1 (tau/tau_w)^3 (tau'/mu) dt`.
    /// Both factors in its integrand are nonnegative and bounded by one. This
    /// change of coordinate needs only ONE constitutive inverse for total flow.
    /// The constant-viscosity limit uses the exact Poiseuille integral.
    ///
    /// The returned pressure is obtained from the accepted wall source response;
    /// its difference from the request is retained. Thus the profile, flux and
    /// `pressure * flow` heat describe the same accepted constitutive state.
    /// Unsupported, unresolved or nonconverged states return no partial profile.
    pub fn profile_at_pressure(
        &self,
        pressure_drop_pa: f64,
        controls: CarreauPipeControls,
    ) -> Result<CarreauPipeProfile, MaterialError> {
        self.profile_at_pressure_poll(pressure_drop_pa, controls, || Ok(()))
    }

    // The Component uses the same arithmetic and source admissions, with
    // cancellation at every quadrature sample. The constitutive inverse has
    // its existing hard limit of 128 evaluations between polls.
    fn profile_at_pressure_poll(
        &self,
        pressure_drop_pa: f64,
        controls: CarreauPipeControls,
        mut poll: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauPipeProfile, MaterialError> {
        poll()?;
        controls.admit()?;
        if !pressure_drop_pa.is_finite() {
            return Err(parameters("Carreau pipe pressure drop must be finite"));
        }
        let target = if pressure_drop_pa == 0.0 {
            0.0
        } else {
            product_ratio(
                &[pressure_drop_pa.abs(), self.radius_m, 0.5],
                &[self.length_m],
            )?
        };
        let (wall, iterations) = self
            .prepared
            .inverse(target, self.prepared.ceiling, controls)?;
        self.profile_from_wall_poll(wall, pressure_drop_pa, iterations, controls, &mut poll)
    }

    /// Resolve a steady profile at an authored signed wall shear rate.
    ///
    /// The rate is admitted directly against the original prepared source.
    /// In particular, the exact declared ceiling can be queried without an
    /// external pressure reconstruction rounding beyond its supported stress.
    /// The original positive radial quadrature still owns flow and heat.
    pub fn profile_at_wall_rate(
        &self,
        signed_wall_rate_per_s: f64,
        controls: CarreauPipeControls,
    ) -> Result<CarreauPipeProfile, MaterialError> {
        controls.admit()?;
        if !signed_wall_rate_per_s.is_finite() {
            return Err(parameters("Carreau pipe wall rate must be finite"));
        }
        let wall = self.prepared.sample_rate(signed_wall_rate_per_s.abs())?;
        let pressure = if wall.rate == 0.0 {
            0.0
        } else {
            product_ratio(&[2.0, self.length_m, wall.stress()], &[self.radius_m])?
                .copysign(signed_wall_rate_per_s)
        };
        self.profile_from_wall(wall, pressure, 0, controls)
    }

    /// Resolve the pressure required by a prescribed signed volume flow.
    ///
    /// The bounded Newton solve works directly in wall shear rate, reusing the
    /// original positive flow quadrature. It does not nest pressure-to-rate
    /// inversions or use finite-difference derivatives. Its derivative follows
    /// from the Rabinowitsch relation: if Q = pi R³ g_w I, then
    /// dQ/dg_w = pi R³ (tau'_w/mu_w) (1 - 3 I).
    /// For the admitted shear-thinning laws, I <= 1/4, so this expression does
    /// not subtract nearly equal terms. The zero-flow pressure and heat are
    /// exactly zero; the zero-flow conductance remains finite.
    ///
    /// `relative_flow_tolerance` applies to both the computed flow residual and
    /// the Newton wall-rate correction. Quadrature accuracy is independently
    /// controlled by `controls`, and remains an estimate, not a certified error.
    /// `maximum_flow_iterations` is in 1..=128 and bounds all profile queries.
    /// Unsupported flow, exhausted work or unrepresentable responses refuse.
    pub fn profile_at_flow(
        &self,
        volumetric_flow_m3_s: f64,
        controls: CarreauPipeControls,
        relative_flow_tolerance: f64,
        maximum_flow_iterations: u32,
    ) -> Result<CarreauPipeFlowSolution, MaterialError> {
        controls.admit()?;
        if !volumetric_flow_m3_s.is_finite()
            || !relative_flow_tolerance.is_finite()
            || !(0.0..1.0).contains(&relative_flow_tolerance)
            || relative_flow_tolerance == 0.0
            || !(1..=128).contains(&maximum_flow_iterations)
        {
            return Err(parameters(
                "invalid prescribed Carreau pipe flow or inverse controls",
            ));
        }
        let accepted = |profile: CarreauPipeProfile, flow_iterations| {
            let residual = profile.flow.volumetric_flow_m3_s - volumetric_flow_m3_s;
            CarreauPipeFlowSolution {
                profile,
                requested_volumetric_flow_m3_s: volumetric_flow_m3_s,
                flow_residual_m3_s: residual,
                flow_iterations,
            }
        };
        if volumetric_flow_m3_s == 0.0 {
            return Ok(accepted(self.profile_at_pressure(0.0, controls)?, 0));
        }
        let at_rate = |wall: RateState| {
            let pressure = product_ratio(&[2.0, self.length_m, wall.stress()], &[self.radius_m])?
                .copysign(volumetric_flow_m3_s);
            self.profile_from_wall(wall, pressure, 0, controls)
        };
        let target = volumetric_flow_m3_s.abs();
        // The Newtonian relation is a lower wall-rate estimate for shear
        // thinning. The original declared ceiling supplies the upper bracket.
        let seed = product_ratio(
            &[4.0, target],
            &[
                core::f64::consts::PI,
                self.radius_m,
                self.radius_m,
                self.radius_m,
            ],
        )?;
        let mut low = 0.0;
        let mut high = self.prepared.ceiling.rate;
        let mut bracketed = false;
        let mut rate = seed.min(high);
        for iterations in 1..=maximum_flow_iterations {
            let profile = at_rate(self.prepared.sample_rate(rate)?)?;
            let difference = profile.flow.volumetric_flow_m3_s.abs() - target;
            if rate == self.prepared.ceiling.rate && difference < 0.0 {
                return Err(state_error(
                    "prescribed flow exceeds the Carreau pipe wall-rate envelope",
                ));
            }
            let derivative = profile.wall_rate_flow_derivative()?;
            let correction = difference / derivative;
            if (difference / target).abs() <= relative_flow_tolerance
                && (correction / rate).abs() <= relative_flow_tolerance
            {
                return Ok(accepted(profile, iterations));
            }
            if difference < 0.0 {
                low = rate;
            } else {
                high = rate;
                bracketed = true;
            }
            let candidate = rate - correction;
            let next = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else if !bracketed {
                // Evaluate the source ceiling only if needed to bound the
                // requested flow. A distant ceiling may have unrepresentable
                // pressure/heat even when the requested solution is ordinary.
                high
            } else {
                low + 0.5 * (high - low)
            };
            if next == rate || next <= low || next > high || (bracketed && next == high) {
                return Err(state_error(
                    "Carreau pipe flow inverse lost floating-point resolution",
                ));
            }
            rate = next;
        }
        Err(state_error(
            "Carreau pipe flow inverse exhausted its iteration budget",
        ))
    }

    fn profile_from_wall(
        &self,
        wall: RateState,
        pressure_drop_pa: f64,
        iterations: u32,
        controls: CarreauPipeControls,
    ) -> Result<CarreauPipeProfile, MaterialError> {
        self.profile_from_wall_poll(wall, pressure_drop_pa, iterations, controls, &mut || Ok(()))
    }

    fn profile_from_wall_poll(
        &self,
        wall: RateState,
        pressure_drop_pa: f64,
        iterations: u32,
        controls: CarreauPipeControls,
        poll: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauPipeProfile, MaterialError> {
        poll()?;
        let transition = self.prepared.transition_rate(wall.rate);
        let (integral, change) = if wall.rate == 0.0 {
            (0.25, 0.0)
        } else if self.prepared.constant() {
            (0.25, 0.0)
        } else {
            integrate(
                &mut |t| {
                    poll()?;
                    if t == 0.0 {
                        return Ok(0.0);
                    }
                    let local = if t == 1.0 {
                        wall
                    } else {
                        self.prepared.sample_rate(wall.rate * t)?
                    };
                    let fraction = product_ratio(&[local.rate, local.viscosity], &[wall.stress()])?;
                    Ok(fraction * fraction * fraction * (local.tangent / local.viscosity))
                },
                transition.map(|rate| rate / wall.rate),
                controls,
            )?
        };
        poll()?;
        let (pressure, volume, heat) = if wall.rate == 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let pressure = product_ratio(&[2.0, self.length_m, wall.stress()], &[self.radius_m])?;
            let volume = product_ratio(
                &[
                    core::f64::consts::PI,
                    self.radius_m,
                    self.radius_m,
                    self.radius_m,
                    wall.rate,
                    integral,
                ],
                &[],
            )?;
            let heat = product_ratio(&[pressure, volume], &[])?;
            (
                pressure.copysign(pressure_drop_pa),
                volume.copysign(pressure_drop_pa),
                heat,
            )
        };
        let residual = pressure - pressure_drop_pa;
        if !residual.is_finite()
            || (pressure_drop_pa != 0.0
                && (residual / pressure_drop_pa).abs()
                    > controls.relative_constitutive_tolerance + 16.0 * f64::EPSILON)
        {
            return Err(state_error(
                "Carreau pipe pressure reconstruction failed its constitutive tolerance",
            ));
        }
        let transition_stress_fraction = transition
            .map(|rate| {
                let sample = self.prepared.sample_rate(rate)?;
                product_ratio(&[sample.rate, sample.viscosity], &[wall.stress()])
            })
            .transpose()?;
        Ok(CarreauPipeProfile {
            pipe: *self,
            controls,
            wall,
            transition_stress_fraction,
            normalized_flow_integral: integral,
            flow: CarreauPipeFlow {
                requested_pressure_drop_pa: pressure_drop_pa,
                pressure_drop_pa: pressure,
                pressure_residual_pa: residual,
                volumetric_flow_m3_s: volume,
                wall_stress_pa: wall.stress().copysign(pressure_drop_pa),
                wall_shear_rate_per_s: wall.rate.copysign(pressure_drop_pa),
                dissipated_power_w: heat,
                quadrature_relative_change: change,
                constitutive_iterations: iterations,
            },
        })
    }
}

impl CarreauPipeProfile {
    pub const fn pipe(&self) -> &LaminarCarreauPipe {
        &self.pipe
    }
    pub const fn controls(&self) -> CarreauPipeControls {
        self.controls
    }
    pub const fn flow(&self) -> &CarreauPipeFlow {
        &self.flow
    }

    /// Positive differential hydraulic conductance dQ/d(delta_p) [m³/(s Pa)].
    /// It is even under flow reversal and has the exact Poiseuille limit at
    /// rest. This uses the already accepted radial integral without another
    /// constitutive solve or quadrature; it is not the secant Q/delta_p.
    pub fn pressure_conductance_m3_s_pa(&self) -> Result<f64, MaterialError> {
        product_ratio(
            &[
                core::f64::consts::PI,
                self.pipe.radius_m,
                self.pipe.radius_m,
                self.pipe.radius_m,
                self.pipe.radius_m,
                1.0 - 3.0 * self.normalized_flow_integral,
            ],
            &[2.0, self.pipe.length_m, self.wall.viscosity],
        )
    }

    /// Analytic dQ/dT [m³/(s K)] at fixed frictional pressure and geometry.
    ///
    /// This retains the original shift calibration and accepted radial
    /// profile. For a vertical viscosity shift, changing temperature is
    /// equivalent to dividing pressure by `a_T`, so the derivative is
    /// `-d(log a_T)/dT * pressure * dQ/dpressure`. For a vertical-and-rate
    /// shift, the inverse stress profile scales as `1/a_T`, giving
    /// `-d(log a_T)/dT * Q`. No perturbed source solve or quadrature is needed.
    /// The derivative is odd under pressure reversal and exactly zero at rest.
    pub fn temperature_flow_derivative_m3_s_k(&self) -> Result<f64, MaterialError> {
        let pressure = self.flow.pressure_drop_pa;
        if pressure == 0.0 {
            return Ok(0.0);
        }
        let rate = -self.pipe.prepared.temperature.log_derivative_per_k();
        let magnitude = match self.pipe.prepared.law.mode() {
            ViscosityShiftMode::Vertical => product_ratio(
                &[rate, pressure.abs(), self.pressure_conductance_m3_s_pa()?],
                &[],
            )?,
            ViscosityShiftMode::VerticalAndRate => {
                product_ratio(&[rate, self.flow.volumetric_flow_m3_s.abs()], &[])?
            }
        };
        Ok(magnitude.copysign(pressure))
    }

    fn wall_rate_flow_derivative(&self) -> Result<f64, MaterialError> {
        product_ratio(
            &[
                core::f64::consts::PI,
                self.pipe.radius_m,
                self.pipe.radius_m,
                self.pipe.radius_m,
                self.wall.tangent,
                1.0 - 3.0 * self.normalized_flow_integral,
            ],
            &[self.wall.viscosity],
        )
    }

    /// Sample the inclusive physical radius [0,R]. Velocity integrates the same
    /// inverse law in stress coordinates. Factoring (R-r)/R retains a small
    /// positive near-wall velocity even when r/R rounds to one. Local shear and
    /// heat follow tau(r); no shear-rate or viscosity cutoff is introduced.
    pub fn sample(&self, radius_m: f64) -> Result<CarreauPipeSample, MaterialError> {
        if !radius_m.is_finite() || radius_m < 0.0 || radius_m > self.pipe.radius_m {
            return Err(parameters(
                "Carreau pipe sample radius must lie between centre and wall",
            ));
        }
        let x = radius_m / self.pipe.radius_m;
        let distance = (self.pipe.radius_m - radius_m) / self.pipe.radius_m;
        let stress = if radius_m == 0.0 || self.wall.rate == 0.0 {
            0.0
        } else {
            product_ratio(&[self.wall.stress(), radius_m], &[self.pipe.radius_m])?
        };
        let (local, _) = self
            .pipe
            .prepared
            .inverse(stress, self.wall, self.controls)?;
        let (speed, change) = if self.wall.rate == 0.0 || radius_m == self.pipe.radius_m {
            (0.0, 0.0)
        } else if self.pipe.prepared.constant() {
            (
                product_ratio(
                    &[0.5, self.pipe.radius_m, self.wall.rate, distance, 1.0 + x],
                    &[],
                )?,
                0.0,
            )
        } else {
            let split = self.transition_stress_fraction.map(|s| (s - x) / distance);
            let (integral, change) = integrate(
                &mut |q| {
                    if q == 0.0 {
                        return Ok(local.rate / self.wall.rate);
                    }
                    if q == 1.0 {
                        return Ok(1.0);
                    }
                    let fraction = if x <= 0.5 {
                        x + distance * q
                    } else {
                        1.0 - distance * (1.0 - q)
                    };
                    let (state, _) = self.pipe.prepared.inverse(
                        self.wall.stress() * fraction,
                        self.wall,
                        self.controls,
                    )?;
                    Ok(state.rate / self.wall.rate)
                },
                split,
                self.controls,
            )?;
            (
                product_ratio(
                    &[self.pipe.radius_m, self.wall.rate, distance, integral],
                    &[],
                )?,
                change,
            )
        };
        let heat = if stress == 0.0 || local.rate == 0.0 {
            0.0
        } else {
            product_ratio(&[stress, local.rate], &[])?
        };
        let direction = self.flow.pressure_drop_pa;
        Ok(CarreauPipeSample {
            axial_velocity_m_s: speed.copysign(direction),
            radial_shear_stress_pa: (-stress).copysign(-direction),
            radial_velocity_gradient_per_s: (-local.rate).copysign(-direction),
            dynamic_viscosity_pa_s: local.viscosity,
            dissipation_w_m3: heat,
            quadrature_relative_change: change,
        })
    }
}

// Ordinary SI arithmetic stays on the direct path. Log accumulation rescues
// representable products whose intermediate geometry/constitutive scale would
// overflow or underflow; a truly unrepresentable nonzero response still refuses.
pub(super) fn product_ratio(
    numerators: &[f64],
    denominators: &[f64],
) -> Result<f64, MaterialError> {
    let mut value = 1.0;
    let mut ordinary = true;
    for &factor in numerators {
        if !positive(factor) {
            return Err(state_error(
                "Carreau pipe nonzero response is unrepresentable",
            ));
        }
        value *= factor;
        ordinary &= value.is_normal();
    }
    for &factor in denominators {
        if !positive(factor) {
            return Err(state_error(
                "Carreau pipe response scale is unrepresentable",
            ));
        }
        value /= factor;
        ordinary &= value.is_normal();
    }
    if !ordinary {
        let logarithm: f64 = numerators.iter().map(|&v| det::ln(v)).sum::<f64>()
            - denominators.iter().map(|&v| det::ln(v)).sum::<f64>();
        value = det::exp(logarithm);
    }
    if !positive(value) {
        return Err(state_error(
            "Carreau pipe nonzero response is unrepresentable",
        ));
    }
    Ok(value)
}

pub(super) fn integrate(
    f: &mut impl FnMut(f64) -> Result<f64, MaterialError>,
    split: Option<f64>,
    controls: CarreauPipeControls,
) -> Result<(f64, f64), MaterialError> {
    let (value, change) = if let Some(mid) = split.filter(|&x| x > 0.0 && x < 1.0) {
        let (left, left_change) = integrate_interval(f, 0.0, mid, controls)?;
        let (right, right_change) = integrate_interval(f, mid, 1.0, controls)?;
        (left + right, left_change + right_change)
    } else {
        integrate_interval(f, 0.0, 1.0, controls)?
    };
    if !positive(value) || !change.is_finite() || change < 0.0 {
        return Err(state_error(
            "Carreau pipe positive integral is unrepresentable",
        ));
    }
    Ok((value, change / value))
}

// Nested positive Simpson sums reuse all previous samples. No extrapolated
// negative weights or numerical heat are introduced. Two consecutive observed
// agreements are required; exhaustion returns no accepted flow/profile sample.
fn integrate_interval(
    f: &mut impl FnMut(f64) -> Result<f64, MaterialError>,
    lower: f64,
    upper: f64,
    controls: CarreauPipeControls,
) -> Result<(f64, f64), MaterialError> {
    let width = upper - lower;
    let mut checked = |x| {
        let value = f(x)?;
        if !value.is_finite() || value < 0.0 {
            return Err(state_error(
                "Carreau pipe quadrature received a nonfinite or negative integrand",
            ));
        }
        Ok(value)
    };
    let mut trapezoid = 0.5 * (checked(lower)? + checked(upper)?);
    let mut previous = None;
    let mut agreements = 0;
    for level in 1..=controls.maximum_quadrature_refinements {
        let panels = 1_u32 << level;
        let h = 1.0 / f64::from(panels);
        let mut sum = 0.0;
        let mut compensation = 0.0;
        for i in (1..panels).step_by(2) {
            let value = checked(lower + width * (f64::from(i) * h))?;
            let corrected = value - compensation;
            let next = sum + corrected;
            compensation = (next - sum) - corrected;
            sum = next;
        }
        let simpson = (trapezoid + 4.0 * h * sum) / 3.0;
        if let Some(old) = previous {
            let change = f64::abs(simpson - old);
            if simpson > 0.0 && change / simpson <= controls.relative_quadrature_tolerance {
                agreements += 1;
                if agreements == 2 {
                    return Ok((width * simpson, width * change));
                }
            } else {
                agreements = 0;
            }
        }
        previous = Some(simpson);
        trapezoid = 0.5 * trapezoid + h * sum;
    }
    Err(state_error(
        "Carreau pipe quadrature exhausted its refinement budget",
    ))
}
