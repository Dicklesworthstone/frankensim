//! First-order liquid compressibility and acoustic radiation back-reaction.
//!
//! With `U=R'` and `Delta=p_g+p_v-p_inf-2 sigma/R-4 mu U/R`,
//! `(1-U/c) R R'' + 3/2 (1-U/(3c)) U^2 = (1+U/c) Delta/rho
//! + R/(rho c) dDelta/dt`. This is the Newtonian specialization of
//! Warnez and Johnsen (2015), equations 1--2:
//! <https://pmc.ncbi.nlm.nih.gov/articles/PMC4474959/>.
//!
//! Differentiation includes both the imposed pressure and viscous wall stress.
//! The latter contributes `4 mu/(rho c)` to the acceleration denominator.
//! Temperature, saturation pressure and other liquid coefficients remain
//! frozen. This is a low-Mach spherical approximation, not shock-resolving
//! dynamics; strong compression, thermal transport, changing vapor inventory
//! and nonspherical collapse remain outside its domain. The same conservative
//! `|U|/c <= 0.1` ceiling as RP applies. The caller's imposed field must also be
//! spatially uniform over the small bubble; local time derivatives cannot
//! certify this assumption or detect unresolved drive bandwidth.
//!
//! Acoustic observation remains passive. Radiation back-reaction is already
//! in this equation; subtracting a second monopole power loss double-counts it.
//! The common incompressible kinetic-energy diagnostic does not include the
//! compressible near field and outgoing wave energy.

use super::{
    BubbleGas, BubbleLimits, BubbleLinearization, BubbleNumerics, BubbleObservation, BubbleState,
    BubbleStep, CavitationError, CavitationState, RadialDynamics, RayleighPlessetBubble,
    SphericalBubbleModel, finite, observation, positive, reference_linearization, stepper,
    validate_state,
};

/// A smooth imposed absolute pressure and its analytic time derivative.
///
/// Both values must refer to the same time and undisturbed local field. Finite
/// tensile pressures are admitted. A discontinuous pressure jump is not a
/// differentiable forcing and must not be disguised as a zero derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblePressureSample {
    pub absolute_pressure_pa: f64,
    pub time_derivative_pa_s: f64,
}

impl BubblePressureSample {
    pub fn try_new(
        absolute_pressure_pa: f64,
        time_derivative_pa_s: f64,
    ) -> Result<Self, CavitationError> {
        let sample = Self {
            absolute_pressure_pa,
            time_derivative_pa_s,
        };
        sample.validate()?;
        Ok(sample)
    }

    fn validate(self) -> Result<(), CavitationError> {
        if !self.absolute_pressure_pa.is_finite() || !self.time_derivative_pa_s.is_finite() {
            return Err(CavitationError::InvalidInput(
                "imposed absolute pressure and its time derivative",
            ));
        }
        Ok(())
    }
}

/// Small-perturbation liquid bounds in addition to the shared geometric domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KellerMiksisLimits {
    pub bubble: BubbleLimits,
    /// Largest `|p_wall-p_inf|/(rho c^2)`, in `(0, 0.1]`. Checking this even
    /// at zero wall velocity prevents a large compression from passing only
    /// because its initial wall Mach number happens to be zero.
    pub maximum_pressure_ratio: f64,
}

/// Material-driven Keller--Miksis bubble using the original shared DP5(4).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KellerMiksisBubble {
    liquid: CavitationState,
    gas: BubbleGas,
    limits: KellerMiksisLimits,
}

impl KellerMiksisBubble {
    pub fn try_new(
        liquid: CavitationState,
        gas: BubbleGas,
        limits: KellerMiksisLimits,
    ) -> Result<Self, CavitationError> {
        // Reuse the original material, vapor-density and spherical-domain
        // admission; compressibility does not waive those assumptions.
        RayleighPlessetBubble::try_new(liquid, gas, limits.bubble)?;
        if !positive(limits.maximum_pressure_ratio) || limits.maximum_pressure_ratio > 0.1 {
            return Err(CavitationError::InvalidInput(
                "Keller-Miksis pressure/bulk-modulus limit",
            ));
        }
        Ok(Self {
            liquid,
            gas,
            limits,
        })
    }

    pub const fn liquid(self) -> CavitationState {
        self.liquid
    }
    pub const fn gas(self) -> BubbleGas {
        self.gas
    }
    pub const fn limits(self) -> BubbleLimits {
        self.limits.bubble
    }
    pub const fn compressibility_limits(self) -> KellerMiksisLimits {
        self.limits
    }

    /// Exact linearization of this implemented equation at the gas reference
    /// equilibrium. It includes viscous inertia, radiation damping AND the
    /// pressure-derivative forcing; it adds no thermal damping surrogate.
    pub fn linearization(&self) -> Result<BubbleLinearization, CavitationError> {
        let mut response = reference_linearization(self.liquid, self.gas, self.limits.bubble)?;
        let tau = response.equilibrium_radius_m / self.liquid.liquid_sound_speed_m_s();
        if !positive(tau) {
            return Err(CavitationError::NonFiniteDerived(
                "bubble acoustic crossing time",
            ));
        }
        response.inertia_multiplier = finite(
            1.0 + response.damping_coefficient_per_s * tau,
            "linear compressible inertia",
        )?;
        response.damping_coefficient_per_s = finite(
            response.damping_coefficient_per_s + response.stiffness_coefficient_per_s2 * tau,
            "linear viscous and radiation damping",
        )?;
        response.pressure_derivative_time_s = tau;
        Ok(response)
    }

    pub fn observe(
        &self,
        state: &BubbleState,
        pressure: BubblePressureSample,
    ) -> Result<BubbleObservation, CavitationError> {
        let (gas_pressure, acceleration) = self.dynamics(state, pressure)?;
        observation(
            self.liquid,
            state,
            pressure.absolute_pressure_pa,
            gas_pressure,
            acceleration,
        )
    }

    /// Advance a smooth prescribed drive with its exact matching derivative.
    pub fn advance<Pressure: Fn(f64) -> BubblePressureSample>(
        &self,
        state: &mut BubbleState,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
    ) -> Result<BubbleStep, CavitationError> {
        self.advance_controlled(state, end_time_s, numerics, pressure, &mut || false)
    }

    /// Fixed-stack adaptive integration with full-interval atomic publication.
    /// The pressure derivative is never estimated by probing a callback.
    pub fn advance_controlled<Pressure, Cancel>(
        &self,
        state: &mut BubbleState,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<BubbleStep, CavitationError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Cancel: FnMut() -> bool,
    {
        stepper::advance(self, state, end_time_s, numerics, pressure, cancelled)
    }

    fn dynamics(
        &self,
        state: &BubbleState,
        pressure: BubblePressureSample,
    ) -> Result<(f64, f64), CavitationError> {
        validate_state(state, self.liquid, self.limits.bubble)?;
        pressure.validate()?;
        let radius = state.radius_m;
        let velocity = state.wall_velocity_m_s;
        let density = self.liquid.liquid_density_kg_m3();
        let sound_speed = self.liquid.liquid_sound_speed_m_s();
        let viscosity = self.liquid.dynamic_viscosity_pa_s();
        let surface_tension = self.liquid.surface_tension_n_m();
        let gas_pressure = self.gas.pressure(radius)?;
        let capillary_pressure = 2.0 * surface_tension / radius;
        // Keep the original equilibrium-preserving pressure ordering.
        let pressure_jump = finite(
            gas_pressure
                - (pressure.absolute_pressure_pa - self.liquid.saturation_pressure_pa()
                    + capillary_pressure)
                - 4.0 * viscosity * velocity / radius,
            "liquid wall pressure jump",
        )?;
        let pressure_ratio = finite(
            ((pressure_jump / density) / sound_speed).abs() / sound_speed,
            "liquid pressure/bulk-modulus ratio",
        )?;
        if pressure_ratio > self.limits.maximum_pressure_ratio {
            return Err(CavitationError::PressureRatioExceeded {
                pressure_ratio,
                maximum_pressure_ratio: self.limits.maximum_pressure_ratio,
            });
        }
        let mach = velocity / sound_speed;
        let denominator = radius * (1.0 - mach) + (4.0 * viscosity / density) / sound_speed;
        if !positive(denominator) {
            return Err(CavitationError::NonFiniteDerived(
                "Keller-Miksis acceleration denominator",
            ));
        }
        let radial_strain_rate = velocity / radius;
        let pressure_rate_without_acceleration =
            (-3.0 * self.gas.polytropic_exponent * gas_pressure + capillary_pressure)
                * radial_strain_rate
                + 4.0 * viscosity * radial_strain_rate * radial_strain_rate
                - pressure.time_derivative_pa_s;
        let acceleration = finite(
            ((1.0 + mach) * pressure_jump / density - (1.5 - 0.5 * mach) * velocity * velocity
                + (radius / sound_speed) * pressure_rate_without_acceleration / density)
                / denominator,
            "Keller-Miksis wall acceleration",
        )?;
        Ok((gas_pressure, acceleration))
    }
}

impl SphericalBubbleModel for KellerMiksisBubble {
    type Pressure = BubblePressureSample;

    fn liquid(&self) -> CavitationState {
        self.liquid
    }
    fn gas(&self) -> BubbleGas {
        self.gas
    }
    fn limits(&self) -> BubbleLimits {
        self.limits.bubble
    }
    fn linearization(&self) -> Result<BubbleLinearization, CavitationError> {
        KellerMiksisBubble::linearization(self)
    }
    fn observe(
        &self,
        state: &BubbleState,
        pressure: BubblePressureSample,
    ) -> Result<BubbleObservation, CavitationError> {
        KellerMiksisBubble::observe(self, state, pressure)
    }
    fn advance_controlled<Pressure, Cancel>(
        &self,
        state: &mut BubbleState,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<BubbleStep, CavitationError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Cancel: FnMut() -> bool,
    {
        KellerMiksisBubble::advance_controlled(
            self, state, end_time_s, numerics, pressure, cancelled,
        )
    }
}

impl RadialDynamics for KellerMiksisBubble {
    fn dynamics(
        &self,
        state: &BubbleState,
        pressure: BubblePressureSample,
    ) -> Result<(f64, f64), CavitationError> {
        KellerMiksisBubble::dynamics(self, state, pressure)
    }
}
