//! Nonlinear spherical bubbles with the ORIGINAL finite gas calorics.
//!
//! The liquid uses the same Rayleigh--Plesset or Keller--Miksis momentum
//! kernel as the polytropic owners. The noncondensable inventory instead obeys
//! `U_g' = Q_heat - p_g V'`, `p_g = m R_g T(U_g)/V`. The original NASA-9
//! source integrates Cv and inverts its energy increment, so formation-energy
//! offsets cannot affect motion. For KM the SAME stage contributes
//! `p_g' = -3 gamma(T) p_g R'/R + R_g Q_heat/(Cv(T) V)`.
//!
//! The gas is a lumped, spatially uniform ideal mixture with fixed mass and
//! composition. Gas/liquid inertia and gas wall-Mach bounds are checked at
//! every stage; incident bandwidth and the uniform-temperature approximation
//! remain caller conditions. Heat is an explicit signed power port, not a
//! derived conduction closure. Liquid temperature, viscosity, surface tension,
//! saturation pressure and sound speed stay frozen. The vapor pressure
//! fraction must remain small: vapor heat, phase mass, latent heat, chemistry,
//! dissolution, gas diffusion and shocks are not supplied by this owner.
//! The KM equation already owns acoustic back-reaction; an observer adds none.
//!
//! Available with `caloric-cavitation`; execution of its new regression tests
//! is pending because the implementation environment became unavailable.
//!
//! See Brennen, *Cavitation and Bubble Dynamics*, sections 2.2--2.4:
//! <https://media.library.caltech.edu/CaltechBOOK:1995.001/chap2.htm>.

use super::{
    BubbleLimits, BubbleNumerics, BubbleObservation, BubblePressureSample, BubbleState,
    CavitationError, CavitationState, KellerMiksisLimits, finite, keller_miksis, observation,
    positive, rayleigh_plesset_acceleration, stepper, validate_state,
};
use fs_material::thermochemical_gas::energy::increment::GasCaloricBranch;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use std::cell::RefCell;

/// Original liquid momentum equation and its original applicability limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CaloricBubbleEquation {
    RayleighPlesset(BubbleLimits),
    KellerMiksis(KellerMiksisLimits),
}

impl CaloricBubbleEquation {
    pub const fn bubble_limits(self) -> BubbleLimits {
        match self {
            Self::RayleighPlesset(limits) => limits,
            Self::KellerMiksis(limits) => limits.bubble,
        }
    }

    fn validate(self, liquid: CavitationState) -> Result<(), CavitationError> {
        let limits = self.bubble_limits();
        limits.validate()?;
        if let Some(vapor) = liquid.saturated_vapor() {
            let ratio = vapor.density_kg_m3() / liquid.liquid_density_kg_m3();
            if ratio > limits.maximum_vapor_density_ratio {
                return Err(CavitationError::VaporDensityRatioExceeded {
                    vapor_liquid_density_ratio: ratio,
                    maximum_vapor_density_ratio: limits.maximum_vapor_density_ratio,
                });
            }
        }
        if let Self::KellerMiksis(limits) = self {
            if !positive(limits.maximum_pressure_ratio) || limits.maximum_pressure_ratio > 0.1 {
                return Err(CavitationError::InvalidInput(
                    "Keller-Miksis pressure/bulk-modulus limit",
                ));
            }
        }
        Ok(())
    }
}

/// Bounds for the noncondensable, uniform-pressure gas approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricBubbleLimits {
    /// Largest rho_g/rho_l, in (0, 0.1], neglecting gas inertia.
    pub maximum_gas_liquid_density_ratio: f64,
    /// Largest p_v/(p_g+p_v), in (0, 0.1], neglecting vapor calorics.
    pub maximum_vapor_pressure_fraction: f64,
    /// Largest |R'|/sqrt(gamma(T) R_g T), in (0, 0.1].
    pub maximum_gas_wall_mach_number: f64,
}

impl CaloricBubbleLimits {
    fn validate(self) -> Result<(), CaloricBubbleError> {
        for value in [
            self.maximum_gas_liquid_density_ratio,
            self.maximum_vapor_pressure_fraction,
            self.maximum_gas_wall_mach_number,
        ] {
            if !positive(value) || value > 0.1 {
                return Err(CaloricBubbleError::InvalidInput(
                    "caloric bubble limits must lie in (0,0.1]",
                ));
            }
        }
        Ok(())
    }
}

/// Shared radial integrator, gas/heat/work error scale and whole-window balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricBubbleNumerics {
    pub radial: BubbleNumerics,
    /// Local absolute error scale for delta-U, integrated heat and gas work [J].
    pub energy_absolute_tolerance_j: f64,
    /// Absolute whole-window gas energy residual budget [J]. Never multiplied
    /// by step count and never scaled by the absolute formation energy.
    pub energy_balance_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CaloricBubbleError {
    Dynamics(CavitationError),
    Gas(GasEnergyError),
    InvalidInput(&'static str),
    OutsideModel {
        quantity: &'static str,
        value: f64,
        maximum: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for CaloricBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "caloric bubble calculation refused: {self:?}")
    }
}
impl std::error::Error for CaloricBubbleError {}

impl From<CavitationError> for CaloricBubbleError {
    fn from(error: CavitationError) -> Self {
        if error == CavitationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Dynamics(error)
        }
    }
}

impl From<GasEnergyError> for CaloricBubbleError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(error)
        }
    }
}

/// Gas state at a native radial stage; available to a deterministic heat law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricGasObservation {
    pub time_s: f64,
    pub radius_m: f64,
    pub wall_velocity_m_s: f64,
    pub volume_m3: f64,
    pub mass_kg: f64,
    pub temperature_k: f64,
    pub gas_pressure_pa: f64,
    pub cv_j_kg_k: f64,
}

/// Original liquid/source observation and coherent finite gas state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricBubbleObservation {
    pub bubble: BubbleObservation,
    pub gas: CaloricGasObservation,
    pub heat_power_w: f64,
}

/// Complete accepted window and independent original-source energy check.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricBubbleStep {
    pub observation: CaloricBubbleObservation,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub rhs_evaluations: usize,
    pub maximum_normalized_local_error: f64,
    /// Original-law m integral Cv(T) dT, not subtraction of formation energies.
    pub gas_energy_change_j: f64,
    pub supplied_heat_j: f64,
    pub work_by_gas_j: f64,
    /// Delta U_g + integral p_g dV - integral Q_heat dt.
    pub gas_energy_residual_j: f64,
    /// Original-source Delta U_g minus the integrated conservative coordinate.
    pub source_energy_residual_j: f64,
}

/// One finite material inventory and one accepted radial clock.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricBubble {
    liquid: CavitationState,
    gas: ClosedGasParcel,
    equation: CaloricBubbleEquation,
    state: BubbleState,
    limits: CaloricBubbleLimits,
    solve: GasEnergySolve,
}

/// Borrowing proposal: downstream observation can refuse without cloning the
/// current source model or publishing any of the candidate physical state.
pub struct CaloricBubbleProposal<'a> {
    target: &'a mut CaloricBubble,
    gas: ClosedGasParcel,
    state: BubbleState,
    step: CaloricBubbleStep,
}

impl CaloricBubbleProposal<'_> {
    pub const fn step(&self) -> &CaloricBubbleStep {
        &self.step
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub fn commit(self) -> CaloricBubbleStep {
        self.target.gas = self.gas;
        self.target.state = self.state;
        self.step
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CaloricBubbleError> {
    if cancelled() {
        Err(CaloricBubbleError::Cancelled)
    } else {
        Ok(())
    }
}

fn bound(quantity: &'static str, value: f64, maximum: f64) -> Result<(), CaloricBubbleError> {
    if !value.is_finite() || value < 0.0 || value > maximum {
        return Err(CaloricBubbleError::OutsideModel {
            quantity,
            value,
            maximum,
        });
    }
    Ok(())
}

fn volume(radius: f64) -> Result<f64, CaloricBubbleError> {
    let value = 4.0 * core::f64::consts::PI * radius * radius * (radius / 3.0);
    if !positive(value) {
        return Err(CaloricBubbleError::InvalidInput(
            "bubble volume is not representable",
        ));
    }
    Ok(value)
}

impl CaloricBubble {
    /// The gas volume must agree with the supplied spherical radius to 64
    /// epsilons. Inventory is never recreated from an equilibrium pressure.
    /// The declared water component is refused; other species phase validity
    /// remains the original material owner's/caller's responsibility.
    pub fn new(
        liquid: CavitationState,
        gas: ClosedGasParcel,
        equation: CaloricBubbleEquation,
        state: BubbleState,
        initial_pressure: BubblePressureSample,
        limits: CaloricBubbleLimits,
        solve: GasEnergySolve,
    ) -> Result<Self, CaloricBubbleError> {
        equation.validate(liquid)?;
        limits.validate()?;
        solve.validate()?;
        let spherical_volume = volume(state.radius_m)?;
        if (spherical_volume - gas.volume_m3()).abs()
            > 64.0 * f64::EPSILON * spherical_volume.max(gas.volume_m3())
            || gas.model().water_species().is_some()
        {
            return Err(CaloricBubbleError::InvalidInput(
                "matching spherical volume and noncondensable gas inventory required",
            ));
        }
        let temperature = gas.state().state().temperature;
        if temperature < solve.minimum_temperature_k || temperature > solve.maximum_temperature_k {
            return Err(CaloricBubbleError::InvalidInput(
                "initial gas outside caloric bracket",
            ));
        }
        gas.model().prepare_caloric_branch(solve, || false)?;
        let owner = Self {
            liquid,
            gas,
            equation,
            state,
            limits,
            solve,
        };
        owner.observe(initial_pressure, 0.0)?;
        Ok(owner)
    }

    pub const fn liquid(&self) -> CavitationState {
        self.liquid
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn state(&self) -> BubbleState {
        self.state
    }
    pub const fn equation(&self) -> CaloricBubbleEquation {
        self.equation
    }
    pub const fn limits(&self) -> CaloricBubbleLimits {
        self.limits
    }
    pub const fn solve(&self) -> GasEnergySolve {
        self.solve
    }

    /// Accepted instantaneous state. KM needs the instantaneous heat rate as
    /// well as the imposed pressure derivative to compute its acceleration.
    pub fn observe(
        &self,
        pressure: BubblePressureSample,
        heat_power_w: f64,
    ) -> Result<CaloricBubbleObservation, CaloricBubbleError> {
        let gas = self.gas.state().state();
        let cv = self
            .gas
            .state()
            .thermodynamics()
            .properties()
            .mass_specific_cv()
            .value();
        let gas = self.gas_observation(self.state, gas.temperature, gas.pressure, cv)?;
        self.observation(self.state, pressure, gas, heat_power_w)
    }

    fn gas_observation(
        &self,
        state: BubbleState,
        temperature_k: f64,
        gas_pressure_pa: f64,
        cv_j_kg_k: f64,
    ) -> Result<CaloricGasObservation, CaloricBubbleError> {
        validate_state(&state, self.liquid, self.equation.bubble_limits())?;
        let volume_m3 = volume(state.radius_m)?;
        let ratio = self.gas.state().state().specific_gas_constant / cv_j_kg_k;
        let gamma = 1.0 + ratio;
        if !positive(temperature_k)
            || !positive(gas_pressure_pa)
            || !positive(cv_j_kg_k)
            || !positive(ratio)
            || gamma <= 1.0
            || gamma > 5.0 / 3.0 + 1e-9
        {
            return Err(CaloricBubbleError::InvalidInput(
                "invalid original gas caloric sample",
            ));
        }
        let density = self.gas.mass_kg() / volume_m3;
        let sound_speed = fs_math::det::sqrt(
            gamma * self.gas.state().state().specific_gas_constant * temperature_k,
        );
        if !positive(density) || !positive(sound_speed) {
            return Err(CaloricBubbleError::InvalidInput(
                "unrepresentable gas inertia or sound speed",
            ));
        }
        bound(
            "gas/liquid density ratio",
            density / self.liquid.liquid_density_kg_m3(),
            self.limits.maximum_gas_liquid_density_ratio,
        )?;
        bound(
            "gas wall Mach number",
            state.wall_velocity_m_s.abs() / sound_speed,
            self.limits.maximum_gas_wall_mach_number,
        )?;
        let vapor = self.liquid.saturation_pressure_pa();
        let total = finite(gas_pressure_pa + vapor, "total bubble pressure")?;
        bound(
            "vapor pressure fraction",
            vapor / total,
            self.limits.maximum_vapor_pressure_fraction,
        )?;
        Ok(CaloricGasObservation {
            time_s: state.time_s,
            radius_m: state.radius_m,
            wall_velocity_m_s: state.wall_velocity_m_s,
            volume_m3,
            mass_kg: self.gas.mass_kg(),
            temperature_k,
            gas_pressure_pa,
            cv_j_kg_k,
        })
    }

    fn observation(
        &self,
        state: BubbleState,
        pressure: BubblePressureSample,
        gas: CaloricGasObservation,
        heat_power_w: f64,
    ) -> Result<CaloricBubbleObservation, CaloricBubbleError> {
        pressure.validate()?;
        if !heat_power_w.is_finite() {
            return Err(CaloricBubbleError::InvalidInput(
                "finite gas heat power required",
            ));
        }
        let acceleration = match self.equation {
            CaloricBubbleEquation::RayleighPlesset(_) => rayleigh_plesset_acceleration(
                self.liquid,
                &state,
                pressure.absolute_pressure_pa,
                gas.gas_pressure_pa,
            )?,
            CaloricBubbleEquation::KellerMiksis(limits) => {
                let ratio = self.gas.state().state().specific_gas_constant / gas.cv_j_kg_k;
                let heat_rate = finite(
                    (ratio / gas.volume_m3) * heat_power_w,
                    "gas heat contribution to pressure rate",
                )?;
                if heat_power_w != 0.0 && heat_rate == 0.0 {
                    return Err(CaloricBubbleError::InvalidInput(
                        "gas heat pressure rate underflows",
                    ));
                }
                keller_miksis::acceleration_with_gas(
                    self.liquid,
                    &state,
                    pressure,
                    limits.maximum_pressure_ratio,
                    gas.gas_pressure_pa,
                    -3.0 * (1.0 + ratio) * gas.gas_pressure_pa,
                    heat_rate,
                )?
            }
        };
        Ok(CaloricBubbleObservation {
            bubble: observation(
                self.liquid,
                &state,
                pressure.absolute_pressure_pa,
                gas.gas_pressure_pa,
                acceleration,
            )?,
            gas,
            heat_power_w,
        })
    }

    fn stage(
        &self,
        branch: &GasCaloricBranch<'_>,
        state: BubbleState,
        energy_change_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CaloricGasObservation, CaloricBubbleError> {
        validate_state(&state, self.liquid, self.equation.bubble_limits())?;
        let density = self.gas.mass_kg() / volume(state.radius_m)?;
        let increment = energy_change_j / self.gas.mass_kg();
        if !positive(density)
            || !increment.is_finite()
            || (energy_change_j != 0.0 && increment == 0.0)
        {
            return Err(CaloricBubbleError::InvalidInput(
                "unrepresentable gas energy coordinate",
            ));
        }
        let sample = branch.sample_density_internal_energy_increment(
            self.gas.point(),
            density,
            increment,
            cancelled,
        )?;
        self.gas_observation(
            state,
            sample.temperature_k,
            sample.pressure_pa,
            sample.cv_j_kg_k,
        )
    }

    /// Integrate radial state and actual caloric energy on one adaptive clock.
    /// Heat is watts INTO the gas. Callbacks must be deterministic and bounded:
    /// rejected trials can revisit times, and callback side effects are not
    /// transactional. The shared maximum step must resolve both drives.
    pub fn advance_controlled<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &CaloricBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<CaloricBubbleStep, CaloricBubbleError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, CaloricBubbleError>,
        Cancel: FnMut() -> bool,
    {
        Ok(self
            .prepare_advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?
            .commit())
    }

    /// Prepare all physical state and the source/heat/work balance before any
    /// mutation. Dropping a proposal preserves the original gas and radial state.
    pub fn prepare_advance_controlled<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &CaloricBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<CaloricBubbleProposal<'_>, CaloricBubbleError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, CaloricBubbleError>,
        Cancel: FnMut() -> bool,
    {
        poll(cancelled)?;
        numerics.radial.validate()?;
        let inner_tolerance = self.gas.mass_kg() * self.solve.energy_tolerance_j_kg;
        if !positive(numerics.energy_absolute_tolerance_j)
            || !positive(numerics.energy_balance_tolerance_j)
            || !positive(inner_tolerance)
            || inner_tolerance > numerics.energy_absolute_tolerance_j
            || inner_tolerance > numerics.energy_balance_tolerance_j
        {
            return Err(CaloricBubbleError::InvalidInput(
                "positive energy scales must resolve the original gas inversion tolerance",
            ));
        }
        let result = {
            let branch = self
                .gas
                .model()
                .prepare_caloric_branch(self.solve, &mut *cancelled)?;
            // One source is polled both at RK stages and inside caloric roots.
            // The temporary RefCell borrows never overlap.
            let cancellation = RefCell::new(&mut *cancelled);
            let mut rhs =
                |time_s: f64, values: &[f64; 5]| -> Result<[f64; 5], CaloricBubbleError> {
                    let state = BubbleState {
                        time_s,
                        radius_m: values[0],
                        wall_velocity_m_s: values[1],
                    };
                    let gas = self.stage(&branch, state, values[2], &mut || {
                        (**cancellation.borrow_mut())()
                    })?;
                    let power = heat(time_s, &gas)?;
                    let frame = self.observation(state, pressure(time_s), gas, power)?;
                    let work = finite(
                        gas.gas_pressure_pa * frame.bubble.volume_rate_m3_s,
                        "gas pressure work rate",
                    )?;
                    let energy = finite(power - work, "gas energy rate")?;
                    Ok([
                        values[1],
                        frame.bubble.wall_acceleration_m_s2,
                        energy,
                        power,
                        work,
                    ])
                };
            stepper::advance_array(
                self.state.time_s,
                [
                    self.state.radius_m,
                    self.state.wall_velocity_m_s,
                    0.0,
                    0.0,
                    0.0,
                ],
                end_time_s,
                &numerics.radial,
                [
                    numerics.radial.radius_absolute_tolerance_m,
                    numerics.radial.velocity_absolute_tolerance_m_s,
                    numerics.energy_absolute_tolerance_j,
                    numerics.energy_absolute_tolerance_j,
                    numerics.energy_absolute_tolerance_j,
                ],
                &mut rhs,
                &mut || (**cancellation.borrow_mut())(),
                &|error| match error {
                    CaloricBubbleError::Dynamics(error) => stepper::refinable(error),
                    CaloricBubbleError::OutsideModel { .. }
                    | CaloricBubbleError::Gas(GasEnergyError::OutsideBracket { .. }) => true,
                    _ => false,
                },
            )?
        };
        poll(cancelled)?;
        let state = BubbleState::try_new(result.time_s, result.values[0], result.values[1])?;
        let (gas, actual) = if end_time_s == self.state.time_s {
            (self.gas.clone(), 0.0)
        } else {
            self.gas.at_volume_internal_energy_increment(
                volume(state.radius_m)?,
                result.values[2],
                self.solve,
                &mut *cancelled,
            )?
        };
        let source_residual = finite(actual - result.values[2], "source gas energy residual")?;
        let balance = finite(
            actual + result.values[4] - result.values[3],
            "gas energy balance",
        )?;
        for residual in [source_residual, balance] {
            if residual.abs() > numerics.energy_balance_tolerance_j {
                return Err(CaloricBubbleError::EnergyBalance {
                    residual_j: residual,
                    tolerance_j: numerics.energy_balance_tolerance_j,
                });
            }
        }
        let cv = gas
            .state()
            .thermodynamics()
            .properties()
            .mass_specific_cv()
            .value();
        let gas_frame = self.gas_observation(
            state,
            gas.state().state().temperature,
            gas.state().state().pressure,
            cv,
        )?;
        poll(cancelled)?;
        let observation = self.observation(
            state,
            pressure(state.time_s),
            gas_frame,
            heat(state.time_s, &gas_frame)?,
        )?;
        poll(cancelled)?;
        let step = CaloricBubbleStep {
            observation,
            accepted_steps: result.accepted_steps,
            rejected_steps: result.rejected_steps,
            rhs_evaluations: result.rhs_evaluations,
            maximum_normalized_local_error: result.maximum_normalized_local_error,
            gas_energy_change_j: actual,
            supplied_heat_j: result.values[3],
            work_by_gas_j: result.values[4],
            gas_energy_residual_j: balance,
            source_energy_residual_j: source_residual,
        };
        Ok(CaloricBubbleProposal {
            target: self,
            gas,
            state,
            step,
        })
    }
}

#[cfg(test)]
#[path = "caloric/tests.rs"]
mod tests;
