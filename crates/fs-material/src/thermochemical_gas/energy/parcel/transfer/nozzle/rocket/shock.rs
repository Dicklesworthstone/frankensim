//! Frozen, thermally perfect normal shocks from the original NASA-9 gas.
//!
//! Mass flux j and momentum flux p+j*u determine T(u)=u*(B-u)/R,
//! with B=u1+R*T1/u1. The nontrivial energy root is solved after factoring out
//! u-u1, so the reversible upstream root cannot masquerade as a shock.
//! Enthalpy differences cancel formation references before rounding.
//!
//! This is a stationary, planar, adiabatic normal shock in the nozzle frame.
//! It supplies no reaction relaxation, radiation, viscous shock thickness,
//! oblique-shock geometry or boundary-layer separation model.
use super::offdesign::{Station, root};
use super::*;

/// Accepted source-resolved conservative normal-shock jump.
#[derive(Clone, Debug, PartialEq)]
pub struct NormalShockFlow {
    upstream: ResolvedThermochemicalGasState,
    downstream: ResolvedThermochemicalGasState,
    /// Incoming velocity normal to the stationary shock [m/s].
    pub upstream_velocity_m_s: f64,
    /// Outgoing velocity normal to the stationary shock [m/s].
    pub downstream_velocity_m_s: f64,
    /// Supersonic incoming frozen Mach number.
    pub upstream_mach_number: f64,
    /// Subsonic outgoing frozen Mach number.
    pub downstream_mach_number: f64,
    /// Downstream mass flux divided by upstream mass flux, minus one.
    pub relative_mass_flux_residual: f64,
    /// Actual downstream-minus-upstream normal momentum flux [Pa].
    pub momentum_residual_pa: f64,
    /// Actual stagnation enthalpy change across the shock [J/kg].
    pub energy_residual_j_kg: f64,
    /// Actual entropy increase; tiny signed residuals retain their value [J/(kg K)].
    pub entropy_increase_j_kg_k: f64,
    /// Ratio of downstream/upstream stagnation pressures, exp(-delta_s/R).
    pub stagnation_pressure_ratio: f64,
    /// Interior iterations needed for the conservative shock root.
    pub iterations: usize,
}
impl NormalShockFlow {
    /// Original upstream material/source state.
    #[must_use]
    pub const fn upstream(&self) -> &ResolvedThermochemicalGasState {
        &self.upstream
    }
    /// Original downstream material/source state.
    #[must_use]
    pub const fn downstream(&self) -> &ResolvedThermochemicalGasState {
        &self.downstream
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct RawShock {
    pub upstream: Station,
    pub downstream: Station,
    pub iterations: usize,
}

pub(super) fn calculate(
    parcel: &ClosedGasParcel,
    upstream: Station,
    stagnation_temperature: f64,
    controls: GasEnergySolve,
    total_iterations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RawShock, GasEnergyError> {
    let (t1, p1, u1) = (upstream.temperature, upstream.pressure, upstream.velocity);
    if !positive(t1) || !positive(p1) || !positive(u1) || upstream.mach() <= 1.0 {
        return Err(invalid(
            "normal shock requires finite positive supersonic upstream flow",
        ));
    }
    let thermo = parcel.model().thermodynamics();
    let r = parcel.state().state().specific_gas_constant;
    let properties = thermo
        .evaluate_properties(Temperature::new(t1), Pressure::new(p1))
        .map_err(source_error)?;
    let cp1 = properties.mass_specific_cp().value();
    let lower = r * t1 / u1;
    let b = u1 + lower;
    let peak_temperature = (0.25 * b / r) * b;
    if !positive(b) || !peak_temperature.is_finite() {
        return Err(invalid(
            "normal shock momentum invariant is unrepresentable",
        ));
    }
    // If the momentum parabola rises above total temperature, stop at its
    // first crossing of T0. No trial asks the caloric source above T0 merely
    // because an algebraic continuation of the momentum relation would.
    let upper = if peak_temperature > stagnation_temperature {
        let discriminant = b * b - 4.0 * r * stagnation_temperature;
        2.0 * r * stagnation_temperature / (b + det::sqrt(discriminant))
    } else {
        u1
    };
    if !(upper > lower && upper <= u1) {
        return Err(invalid(
            "normal shock has no representable nontrivial energy bracket",
        ));
    }
    let before = *total_iterations;
    let u2 = root(
        lower,
        upper,
        controls.energy_tolerance_j_kg / u1.max(1.0),
        controls.maximum_iterations,
        total_iterations,
        cancelled,
        |u| {
            let delta_t = (u - u1) * (lower - u) / r;
            let t = t1 + delta_t;
            let actual_delta_t = t - t1;
            let secant_cp = if actual_delta_t == 0.0 {
                cp1
            } else {
                thermo
                    .specific_enthalpy_change(Temperature::new(t1), Temperature::new(t))
                    .map_err(source_error)?
                    .value()
                    / actual_delta_t
            };
            if !positive(secant_cp) {
                return Err(invalid("normal shock has nonpositive secant heat capacity"));
            }
            // [h(T(u))-h(T1)+(u²-u1²)/2] / (u-u1).
            // The exact upstream limit uses cp(T1), eliminating the trivial root.
            Ok(secant_cp * (lower - u) / r + 0.5 * (u + u1))
        },
    )?;
    let temperature = t1 + (u2 - u1) * (lower - u2) / r;
    let pressure = p1 * (b - u2) / lower;
    let properties = thermo
        .evaluate_properties(Temperature::new(temperature), Pressure::new(pressure))
        .map_err(source_error)?;
    let gamma = properties.mass_specific_cp().value() / properties.mass_specific_cv().value();
    if !positive(temperature) || !positive(pressure) || !gamma.is_finite() || gamma <= 1.0 {
        return Err(invalid("normal shock downstream state is unrepresentable"));
    }
    Ok(RawShock {
        upstream,
        downstream: Station {
            temperature,
            pressure,
            velocity: u2,
            sound_speed: det::sqrt(gamma * r * temperature),
        },
        iterations: *total_iterations - before,
    })
}

pub(super) fn resolve(
    parcel: &ClosedGasParcel,
    shock: RawShock,
    controls: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NormalShockFlow, GasEnergyError> {
    poll(cancelled)?;
    let a = shock.upstream;
    let b = shock.downstream;
    let upstream =
        parcel
            .model()
            .evaluate(&updated_point(parcel.point(), a.temperature, a.pressure)?)?;
    poll(cancelled)?;
    let downstream =
        parcel
            .model()
            .evaluate(&updated_point(parcel.point(), b.temperature, b.pressure)?)?;
    let thermo = parcel.model().thermodynamics();
    let ja = upstream.state().density * a.velocity;
    let jb = downstream.state().density * b.velocity;
    let mass = jb / ja - 1.0;
    let momentum = (b.pressure - a.pressure) + (jb * b.velocity - ja * a.velocity);
    let energy = thermo
        .specific_enthalpy_change(
            Temperature::new(a.temperature),
            Temperature::new(b.temperature),
        )
        .map_err(source_error)?
        .value()
        + 0.5 * (b.velocity - a.velocity) * (b.velocity + a.velocity);
    let entropy = thermo
        .specific_entropy_change(
            Temperature::new(a.temperature),
            Pressure::new(a.pressure),
            Temperature::new(b.temperature),
            Pressure::new(b.pressure),
        )
        .map_err(source_error)?
        .value();
    let r = upstream.state().specific_gas_constant;
    let ratio = det::exp(-entropy / r);
    let (mach1, mach2) = (
        a.velocity / upstream.state().sound_speed,
        b.velocity / downstream.state().sound_speed,
    );
    let momentum_scale = a.pressure + ja * a.velocity;
    let mach_tolerance = 4.0 * controls.energy_tolerance_j_kg
        / (downstream.state().sound_speed * downstream.state().sound_speed);
    if !mass.is_finite()
        || mass.abs() > 256.0 * f64::EPSILON
        || !momentum.is_finite()
        || momentum.abs() > 256.0 * f64::EPSILON * momentum_scale
        || !energy.is_finite()
        || energy.abs() > 4.0 * controls.energy_tolerance_j_kg
        || !entropy.is_finite()
        || entropy < -2.0 * controls.energy_tolerance_j_kg / a.temperature
        || !positive(ratio)
        || !positive(mach1)
        || mach1 <= 1.0
        || !positive(mach2)
        || mach2 > 1.0 + mach_tolerance
        || b.temperature < a.temperature
        || b.pressure < a.pressure
    {
        return Err(invalid(
            "normal shock failed mass, momentum, energy or entropy admission",
        ));
    }
    poll(cancelled)?;
    Ok(NormalShockFlow {
        upstream,
        downstream,
        upstream_velocity_m_s: a.velocity,
        downstream_velocity_m_s: b.velocity,
        upstream_mach_number: mach1,
        downstream_mach_number: mach2,
        relative_mass_flux_residual: mass,
        momentum_residual_pa: momentum,
        energy_residual_j_kg: energy,
        entropy_increase_j_kg_k: entropy,
        stagnation_pressure_ratio: ratio,
        iterations: shock.iterations,
    })
}

impl ClosedGasParcel {
    /// Resolve a stationary normal shock with this parcel as the incoming static
    /// gas. Its model, composition, context and source domains are retained.
    /// The positive normal speed must be supersonic. The admitted temperature
    /// range must include the corresponding stagnation temperature. This is a
    /// pure constitutive query: the parcel inventory and state remain unchanged.
    pub fn normal_shock_flow(
        &self,
        upstream_velocity_m_s: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalShockFlow, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let state = self.state().state();
        if !positive(upstream_velocity_m_s)
            || upstream_velocity_m_s <= state.sound_speed
            || state.temperature < controls.minimum_temperature_k
            || state.temperature > controls.maximum_temperature_k
        {
            return Err(invalid(
                "normal shock requires an admitted supersonic upstream state",
            ));
        }
        let mut iterations = 0;
        let kinetic = 0.5 * upstream_velocity_m_s * upstream_velocity_m_s;
        if !positive(kinetic) {
            return Err(invalid("normal shock kinetic energy is unrepresentable"));
        }
        let thermo = self.model().thermodynamics();
        let t0 = root(
            state.temperature,
            controls.maximum_temperature_k,
            controls.energy_tolerance_j_kg,
            controls.maximum_iterations,
            &mut iterations,
            &mut cancelled,
            |t| {
                Ok(thermo
                    .specific_enthalpy_change(
                        Temperature::new(state.temperature),
                        Temperature::new(t),
                    )
                    .map_err(source_error)?
                    .value()
                    - kinetic)
            },
        )?;
        let upstream = Station {
            temperature: state.temperature,
            pressure: state.pressure,
            velocity: upstream_velocity_m_s,
            sound_speed: state.sound_speed,
        };
        let raw = calculate(
            self,
            upstream,
            t0,
            controls,
            &mut iterations,
            &mut cancelled,
        )?;
        let mut result = resolve(self, raw, controls, &mut cancelled)?;
        result.iterations = iterations;
        Ok(result)
    }
}
