//! Homogeneous real-gas normal-shock compression on the original NASA/PR laws.
//!
//! In the shock frame, mass and momentum give a Rayleigh pressure at each
//! candidate density. The original material owner reconstructs the corresponding
//! Hugoniot internal-energy increment, including BOTH density departures. A
//! bounded scalar solve then matches that Rayleigh pressure to the same EOS.
//! No constant gamma, ideal acoustic speed or absolute formation-energy
//! subtraction enters the jump. The final state must also close stagnation
//! enthalpy and satisfy compressive, entropy and classical sonic admission.
//!
//! Conservation equations: NASA TN D-8274, *Real-gas effects associated with
//! one-dimensional transonic flow of cryogenic nitrogen*, section "Flow
//! Through Normal Shocks", https://ntrs.nasa.gov/citations/19770008402.
//! The thermal/caloric closure here is the caller's original Peng--Robinson
//! and NASA source; that report does not calibrate this material model.
//!
//! The caller declares one continuous homogeneous mechanically/calorically
//! stable branch over the density and temperature brackets. Native local
//! stability checks do not establish global phase stability. Phase fronts,
//! chemistry, vibrational nonequilibrium and finite shock thickness are not
//! supplied. This boundary-state calculation does not evolve a spatial field.

use fs_matdb::QueryPoint;
use fs_material::MaterialError;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel,
};
use fs_qty::{Density, Velocity};

/// Explicit source domain, compression branch, physical residuals and work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasShockControl {
    /// Native temperature/energy inversion at every candidate density.
    pub energy: GasEnergySolve,
    /// Strictly increasing rho2/rho1 bounds, BOTH strictly greater than one.
    /// Excluding identity is essential: it solves the jump equations at every
    /// upstream Mach number but is not a shock.
    pub compression_ratio_bracket: [f64; 2],
    /// Absolute downstream-minus-upstream mass-flux allowance [kg/(m2 s)].
    pub mass_flux_tolerance_kg_m2_s: f64,
    /// Absolute momentum-flux allowance [Pa].
    pub momentum_tolerance_pa: f64,
    /// Absolute Hugoniot AND stagnation-enthalpy residual allowance [J/kg].
    pub jump_energy_tolerance_j_kg: f64,
    /// Allowed numerical shortfall from nonnegative Delta s [J/(kg K)].
    /// The actual signed change is reported without clamping.
    pub entropy_tolerance_j_kg_k: f64,
    /// At most 128 safeguarded density iterations, plus both endpoints.
    pub maximum_density_iterations: usize,
}

impl RealGasShockControl {
    pub fn validate(self) -> Result<(), RealGasShockError> {
        self.energy.validate()?;
        let [lo, hi] = self.compression_ratio_bracket;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo <= 1.0
            || hi <= lo
            || !positive(self.mass_flux_tolerance_kg_m2_s)
            || !positive(self.momentum_tolerance_pa)
            || !positive(self.jump_energy_tolerance_j_kg)
            || !positive(self.entropy_tolerance_j_kg_k)
            || !(1..=128).contains(&self.maximum_density_iterations)
        {
            return Err(RealGasShockError::Invalid(
                "strict compression bracket, positive physical tolerances and bounded density work required",
            ));
        }
        Ok(())
    }
}

/// No partially admitted downstream state is returned on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum RealGasShockError {
    Invalid(&'static str),
    Material(MaterialError),
    Energy(GasEnergyError),
    Cancelled,
    /// The explicit compression interval does not bracket the Rayleigh/EOS root.
    OutsideBracket {
        lower_momentum_residual_pa: f64,
        upper_momentum_residual_pa: f64,
    },
    NotConverged {
        density_iterations: usize,
        momentum_residual_pa: f64,
    },
    /// A converged jump fails the classical compressive/entropy conditions.
    Admissibility {
        upstream_mach: f64,
        downstream_mach: f64,
        entropy_change_j_kg_k: f64,
    },
    Conservation {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
}

impl core::fmt::Display for RealGasShockError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "real-gas normal shock: {self:?}")
    }
}
impl std::error::Error for RealGasShockError {}
impl From<MaterialError> for RealGasShockError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<GasEnergyError> for RealGasShockError {
    fn from(error: GasEnergyError) -> Self {
        match error {
            GasEnergyError::Cancelled => Self::Cancelled,
            other => Self::Energy(other),
        }
    }
}

/// Complete source-bound shock endpoints and physical jump defects.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNormalShock {
    upstream: ResolvedThermochemicalRealGasState,
    downstream: ResolvedThermochemicalRealGasState,
    pub upstream_speed: Velocity,
    pub downstream_speed: Velocity,
    pub mass_flux_kg_m2_s: f64,
    pub upstream_mach: f64,
    pub downstream_mach: f64,
    pub internal_energy_increase_j_kg: f64,
    pub entropy_change_j_kg_k: f64,
    pub mass_flux_residual_kg_m2_s: f64,
    pub momentum_flux_residual_pa: f64,
    pub hugoniot_residual_j_kg: f64,
    pub stagnation_enthalpy_residual_j_kg: f64,
    /// Interior outer iterations, excluding the two bracket endpoints.
    pub density_iterations: usize,
    /// Sum of native caloric iterations over all attempted density samples.
    /// Each native inversion additionally evaluates its two thermal endpoints.
    pub caloric_iterations: usize,
}

impl RealGasNormalShock {
    pub const fn upstream(&self) -> &ResolvedThermochemicalRealGasState {
        &self.upstream
    }
    pub const fn downstream(&self) -> &ResolvedThermochemicalRealGasState {
        &self.downstream
    }
    pub fn into_downstream(self) -> ResolvedThermochemicalRealGasState {
        self.downstream
    }
}

struct Jump<'a> {
    model: &'a ThermochemicalPengRobinsonModel,
    context: &'a QueryPoint,
    upstream: ResolvedThermochemicalRealGasState,
    velocity: f64,
    mass_flux: f64,
    dynamic_pressure: f64,
    energy: GasEnergySolve,
}

struct Sample {
    ratio: f64,
    state: ResolvedThermochemicalRealGasState,
    residual: f64,
    slope: f64,
    velocity: f64,
    caloric_iterations: usize,
}

impl Jump<'_> {
    fn sample(
        &self,
        ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sample, RealGasShockError> {
        poll(cancelled)?;
        let rho1 = self.upstream.density_kg_m3();
        let density = rho1 * ratio;
        // Keep the finite compression fraction before reciprocal subtraction.
        let fraction = (ratio - 1.0) / ratio;
        let delta_volume = fraction / rho1;
        let pressure_jump = self.dynamic_pressure * fraction;
        let rayleigh_pressure = self.upstream.pressure_pa() + pressure_jump;
        let target = (self.upstream.pressure_pa() + 0.5 * pressure_jump) * delta_volume;
        if !positive(density)
            || !positive(delta_volume)
            || !positive(rayleigh_pressure)
            || !positive(target)
        {
            return Err(RealGasShockError::Invalid(
                "compressed density, Rayleigh pressure and Hugoniot increment must be representable",
            ));
        }
        let solved = self
            .model
            .reconstruct_density_internal_energy_increment_between_densities(
                self.context,
                rho1,
                density,
                target,
                self.energy,
                &mut *cancelled,
            )?;
        let caloric_iterations = solved.iterations();
        let state = solved.into_state();
        let velocity = self.mass_flux / density;
        let residual = state.pressure_pa() - rayleigh_pressure;
        // d(p_EOS-p_R)/d(rho2/rho1) along the target Hugoniot energy:
        // rho1*[c_s^2-u2^2-p_T*(p_EOS-p_R)/(rho2^2 Cv)]. The sound speed,
        // thermal derivative and Cv all belong to the same native PR state.
        let thermal = state.root().pressure_derivatives().0;
        let slope = rho1
            * (state.sound_speed_m_s().powi(2)
                - velocity * velocity
                - (residual / density) * (thermal / density) / state.specific_heat_cv_j_kg_k());
        if !positive(velocity) || !residual.is_finite() {
            return Err(RealGasShockError::Invalid(
                "finite shock pressure residual and positive downstream velocity required",
            ));
        }
        // A nonfinite/zero tangent only disables Newton; bisection remains valid.
        Ok(Sample {
            ratio,
            state,
            residual,
            slope,
            velocity,
            caloric_iterations,
        })
    }
}

/// Resolve a stationary normal shock in the caller's shock-relative frame.
///
/// `upstream_speed` is positive toward the shock and must exceed the original
/// adiabatic sound speed. `context` supplies T and retained material axes; its
/// pressure is replaced by the native EOS pressure at `upstream_density`.
/// Source composition, critical properties and pair interactions stay fixed.
/// Every candidate uses direct conserved-density calorics; there is no fitted
/// gamma or temperature surrogate. Newton uses the original thermodynamic
/// tangent with a bounded bracket safeguard; failure/cancellation mutates no
/// input. All temperatures/densities in the brackets must belong to the
/// declared continuous homogeneous branch.
#[allow(clippy::too_many_lines)]
pub fn normal_shock(
    model: &ThermochemicalPengRobinsonModel,
    context: &QueryPoint,
    upstream_density: Density,
    upstream_speed: Velocity,
    controls: RealGasShockControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasNormalShock, RealGasShockError> {
    controls.validate()?;
    poll(&mut cancelled)?;
    let rho1 = upstream_density.value();
    let velocity = upstream_speed.value();
    if !positive(rho1) || !positive(velocity) {
        return Err(RealGasShockError::Invalid(
            "positive finite upstream density and shock-relative speed required",
        ));
    }
    let upstream = model.evaluate_density_temperature(context, rho1)?;
    if velocity <= upstream.sound_speed_m_s() {
        return Err(RealGasShockError::Invalid(
            "classical compression shock requires supersonic upstream flow",
        ));
    }
    let mass_flux = rho1 * velocity;
    let dynamic_pressure = mass_flux * velocity;
    let [lo, hi] = controls.compression_ratio_bracket;
    // Budget the coupled residuals before solving: Delta h_0 = energy_error
    // + pressure_error/rho2; Hugoniot_error = energy_error - dp*Delta v/2.
    let volume_scale = (1.0 / lo).max(0.5 * ((hi - 1.0) / hi)) / rho1;
    let pressure_tolerance = controls
        .momentum_tolerance_pa
        .min(0.25 * controls.jump_energy_tolerance_j_kg / volume_scale);
    let energy = GasEnergySolve {
        energy_tolerance_j_kg: controls
            .energy
            .energy_tolerance_j_kg
            .min(0.25 * controls.jump_energy_tolerance_j_kg),
        ..controls.energy
    };
    if !positive(mass_flux)
        || !positive(dynamic_pressure)
        || !positive(volume_scale)
        || !positive(pressure_tolerance)
        || !positive(energy.energy_tolerance_j_kg)
    {
        return Err(RealGasShockError::Invalid(
            "shock fluxes and coupled numerical tolerances are unrepresentable",
        ));
    }
    let jump = Jump {
        model,
        context,
        upstream,
        velocity,
        mass_flux,
        dynamic_pressure,
        energy,
    };
    let mut low = jump.sample(lo, &mut cancelled)?;
    let mut high = jump.sample(hi, &mut cancelled)?;
    let mut caloric_iterations = low.caloric_iterations + high.caloric_iterations;
    if low.residual > 0.0 || high.residual < 0.0 {
        return Err(RealGasShockError::OutsideBracket {
            lower_momentum_residual_pa: low.residual,
            upper_momentum_residual_pa: high.residual,
        });
    }
    for density_iterations in 0..=controls.maximum_density_iterations {
        poll(&mut cancelled)?;
        // Near the trivial root, a small residual is insufficient. Require
        // the downstream sonic branch before attempting complete admission.
        if low.residual.abs() <= pressure_tolerance && low.velocity < low.state.sound_speed_m_s() {
            return finish(
                jump,
                low,
                controls,
                density_iterations,
                caloric_iterations,
                &mut cancelled,
            );
        }
        if high.residual.abs() <= pressure_tolerance && high.velocity < high.state.sound_speed_m_s()
        {
            return finish(
                jump,
                high,
                controls,
                density_iterations,
                caloric_iterations,
                &mut cancelled,
            );
        }
        let best = if low.residual.abs() < high.residual.abs() {
            &low
        } else {
            &high
        };
        if density_iterations == controls.maximum_density_iterations {
            return Err(RealGasShockError::NotConverged {
                density_iterations,
                momentum_residual_pa: best.residual,
            });
        }
        let width = high.ratio - low.ratio;
        let newton = best.ratio - best.residual / best.slope;
        let ratio = if newton.is_finite()
            && newton > low.ratio + 0.25 * width
            && newton < high.ratio - 0.25 * width
        {
            newton
        } else {
            low.ratio + 0.5 * width
        };
        if !(ratio > low.ratio && ratio < high.ratio) {
            return Err(RealGasShockError::NotConverged {
                density_iterations,
                momentum_residual_pa: best.residual,
            });
        }
        let sample = jump.sample(ratio, &mut cancelled)?;
        caloric_iterations += sample.caloric_iterations;
        if sample.residual <= 0.0 {
            low = sample;
        } else {
            high = sample;
        }
    }
    unreachable!("bounded density loop returns a state or work refusal")
}

fn finish(
    jump: Jump<'_>,
    sample: Sample,
    controls: RealGasShockControl,
    density_iterations: usize,
    caloric_iterations: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealGasNormalShock, RealGasShockError> {
    poll(cancelled)?;
    let upstream = jump.upstream;
    let downstream = sample.state;
    let rho1 = upstream.density_kg_m3();
    let rho2 = downstream.density_kg_m3();
    let t1 = upstream.temperature_k();
    let t2 = downstream.temperature_k();
    let delta_volume = ((rho2 - rho1) / rho2) / rho1;
    let du = jump
        .model
        .specific_internal_energy_change_between_densities(rho1, t1, rho2, t2)?;
    let dh = jump
        .model
        .specific_enthalpy_change_between_densities(rho1, t1, rho2, t2)?;
    let entropy = jump
        .model
        .specific_entropy_change_between_densities(rho1, t1, rho2, t2)?;
    let upstream_mach = jump.velocity / upstream.sound_speed_m_s();
    let downstream_mach = sample.velocity / downstream.sound_speed_m_s();
    if rho2 <= rho1
        || downstream.pressure_pa() <= upstream.pressure_pa()
        || !upstream_mach.is_finite()
        || upstream_mach <= 1.0
        || !positive(downstream_mach)
        || downstream_mach >= 1.0
        || !entropy.is_finite()
        || entropy < -controls.entropy_tolerance_j_kg_k
    {
        return Err(RealGasShockError::Admissibility {
            upstream_mach,
            downstream_mach,
            entropy_change_j_kg_k: entropy,
        });
    }
    let mass_residual = rho2 * sample.velocity - jump.mass_flux;
    let momentum_residual = (downstream.pressure_pa() - upstream.pressure_pa())
        + jump.mass_flux * (sample.velocity - jump.velocity);
    let hugoniot_residual =
        du - (0.5 * downstream.pressure_pa() + 0.5 * upstream.pressure_pa()) * delta_volume;
    let enthalpy_residual =
        dh + 0.5 * (sample.velocity - jump.velocity) * (sample.velocity + jump.velocity);
    for (quantity, residual, tolerance) in [
        (
            "mass flux [kg/(m2 s)]",
            mass_residual,
            controls.mass_flux_tolerance_kg_m2_s,
        ),
        (
            "momentum flux [Pa]",
            momentum_residual,
            controls.momentum_tolerance_pa,
        ),
        (
            "Hugoniot energy [J/kg]",
            hugoniot_residual,
            controls.jump_energy_tolerance_j_kg,
        ),
        (
            "stagnation enthalpy [J/kg]",
            enthalpy_residual,
            controls.jump_energy_tolerance_j_kg,
        ),
    ] {
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(RealGasShockError::Conservation {
                quantity,
                residual,
                tolerance,
            });
        }
    }
    poll(cancelled)?;
    Ok(RealGasNormalShock {
        upstream,
        downstream,
        upstream_speed: Velocity::new(jump.velocity),
        downstream_speed: Velocity::new(sample.velocity),
        mass_flux_kg_m2_s: jump.mass_flux,
        upstream_mach,
        downstream_mach,
        internal_energy_increase_j_kg: du,
        entropy_change_j_kg_k: entropy,
        mass_flux_residual_kg_m2_s: mass_residual,
        momentum_flux_residual_pa: momentum_residual,
        hugoniot_residual_j_kg: hugoniot_residual,
        stagnation_enthalpy_residual_j_kg: enthalpy_residual,
        density_iterations,
        caloric_iterations,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RealGasShockError> {
    if cancelled() {
        Err(RealGasShockError::Cancelled)
    } else {
        Ok(())
    }
}
