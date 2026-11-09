//! Pressure-controlled discharge on the original homogeneous real-gas branch.
//!
//! The incoming acoustic family keeps the outgoing entropy characteristic.
//! Expansion integrates du = -a d(rho)/rho on the native isentrope;
//! compression closes the native Rankine--Hugoniot jump. A wave moving wholly
//! outward cannot impose its remote state on the face. A transonic expansion
//! is sampled at u=a instead of allowing the backpressure to overdrive flow.

use super::wall::Isentrope;
use super::{RealGas, RealGasFlowError, RealGasWallControl};
use crate::dynamics::{Conserved, FaceFlux, Primitive, poll, positive};
use crate::material::nonideal_shock::pressure_compression;
use fs_material::MaterialError;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_math::det;
use fs_qty::{Density, Pressure};

impl RealGas {
    /// Physical flux at a stationary pressure-controlled outlet.
    ///
    /// `outward_normal` is -1 at a left end and +1 at a right end; returned
    /// fluxes retain the global coordinate direction. The original NASA/PR
    /// source supplies pressure, sound speed, entropy and transported energy.
    /// `control.compression.momentum_tolerance_pa` bounds the pressure root.
    /// Subsonic expansion follows the interior isentrope and saturates at the
    /// sonic state. Compression uses a source-admitted shock, including an
    /// incoming shock against supersonic discharge. Backflow requires an
    /// explicit exterior thermodynamic state and is refused here.
    ///
    /// A weak compression below floating-point shock resolution uses the
    /// source acoustic tangent and an admitted isentropic state; its pressure
    /// error must fit the same authored pressure tolerance. Its front uses
    /// the mean native incoming-family speeds, retaining the first-order
    /// nonlinear shock motion near sonic discharge. The caller must supply
    /// a continuous classical homogeneous branch, as for the physical wall.
    /// This does not infer phase stability or a reservoir stagnation law.
    /// Zero backpressure denotes exterior vacuum and samples the same
    /// positive-density sonic face. The original source must cover that face;
    /// no internal vacuum or exterior rarefied expansion is supplied.
    pub fn pressure_outlet_flux(
        &self,
        inside: Conserved,
        outward_normal: f64,
        backpressure: Pressure,
        control: RealGasWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        control.validate(self)?;
        let primitive = self.decode(inside, &mut cancelled)?;
        self.outlet_flux(
            (inside, primitive),
            outward_normal,
            backpressure.value(),
            control,
            &mut cancelled,
        )
    }

    pub(super) fn outlet_flux(
        &self,
        inside: (Conserved, Primitive),
        normal: f64,
        pressure: f64,
        control: RealGasWallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        poll(cancelled)?;
        if !matches!(normal, -1.0 | 1.0) || !pressure.is_finite() || pressure < 0.0 {
            return Err(RealGasFlowError::Invalid(
                "invalid real-gas pressure outlet",
            ));
        }
        let p = Primitive {
            velocity: normal * inside.1.velocity,
            ..inside.1
        };
        if p.velocity < 0.0 {
            return Err(backflow());
        }
        let signal = p.velocity + p.sound_speed;
        if pressure == p.pressure || (p.velocity >= p.sound_speed && pressure <= p.pressure) {
            return physical_flux(inside, signal);
        }
        let context = self
            .context
            .clone()
            .with("temperature", p.temperature)
            .map_err(|error| MaterialError::State {
                what: format!("real-gas outlet temperature coordinate: {error}"),
            })?;
        let curve = Isentrope {
            gas: self,
            context,
            inside: p,
            control,
        };
        let (state, velocity, wave_signal) = if pressure > p.pressure {
            let increment = pressure - p.pressure;
            let acoustic_velocity = increment / (p.density * p.sound_speed);
            let acoustic_limit = control
                .velocity_tolerance_m_s
                .max(32.0 * det::sqrt(f64::EPSILON) * p.sound_speed);
            if acoustic_velocity <= acoustic_limit {
                // Classical weak compression tends continuously to its native
                // acoustic tangent. Never subtract two large shock speeds.
                let z = -acoustic_velocity / p.sound_speed;
                let density = p.density * det::exp(-z);
                if !self.domain.admits_density(density)
                    || density > p.density * control.compression.compression_ratio_bracket[1]
                {
                    return Err(RealGasFlowError::Invalid(
                        "outlet acoustic compression leaves density path",
                    ));
                }
                let energy = control.compression.energy;
                for (temperature, lower) in [
                    (energy.minimum_temperature_k, true),
                    (energy.maximum_temperature_k, false),
                ] {
                    let entropy = self.model.specific_entropy_change_between_densities(
                        p.density,
                        p.temperature,
                        density,
                        temperature,
                    )?;
                    if (lower && entropy > 0.0) || (!lower && entropy < 0.0) {
                        return Err(RealGasFlowError::Invalid(
                            "outlet acoustic compression leaves temperature path",
                        ));
                    }
                }
                let state = curve.state(z, cancelled)?;
                if (state.pressure_pa() - pressure).abs()
                    > control.compression.momentum_tolerance_pa
                {
                    return Err(RealGasFlowError::Invalid(
                        "outlet acoustic pressure exceeds tolerance",
                    ));
                }
                let velocity = p.velocity - acoustic_velocity;
                // A weak shock travels at the mean of the incoming-family
                // speeds to first order in its strength. Both sound speeds
                // come from the original EOS; no fitted gamma is needed.
                // Testing only u_inside >= a_inside loses incoming shocks at
                // exactly sonic and slightly supersonic discharge.
                let shock_speed =
                    0.5 * (p.velocity - p.sound_speed) + 0.5 * (velocity - state.sound_speed_m_s());
                if shock_speed >= 0.0 {
                    return physical_flux(inside, signal.max(shock_speed));
                }
                (state, velocity, signal.max(shock_speed.abs()))
            } else {
                let mut compression = control.compression;
                compression.compression_ratio_bracket[0] = compression.compression_ratio_bracket[0]
                    .min(1.0 + 0.25 * acoustic_velocity / (p.sound_speed + acoustic_velocity));
                if !self
                    .domain
                    .admits_density(p.density * compression.compression_ratio_bracket[1])
                {
                    return Err(RealGasFlowError::Invalid(
                        "outlet compression bracket exceeds density domain",
                    ));
                }
                // Weak waves need enough caloric precision to distinguish the
                // two sonic branches, as in the source physical-wall solve.
                let weak_energy = (acoustic_velocity / 64.0) * acoustic_velocity;
                compression.energy.energy_tolerance_j_kg =
                    compression.energy.energy_tolerance_j_kg.min(weak_energy);
                compression.momentum_tolerance_pa = compression
                    .momentum_tolerance_pa
                    .min(p.density * weak_energy);
                let wave = pressure_compression(
                    &self.model,
                    &curve.context,
                    Density::new(p.density),
                    Pressure::new(pressure),
                    compression,
                    &mut *cancelled,
                )?;
                let shock_speed = p.velocity - wave.upstream_speed.value();
                if shock_speed >= 0.0 {
                    return physical_flux(inside, signal.max(shock_speed));
                }
                // j*Delta(v) avoids cancellation of nearly equal velocities.
                let density = wave.downstream().density_kg_m3();
                let drop = wave.mass_flux_kg_m2_s * ((density - p.density) / density) / p.density;
                (
                    wave.into_downstream(),
                    p.velocity - drop,
                    signal.max(shock_speed.abs()),
                )
            }
        } else {
            expansion(&curve, pressure, cancelled)?
        };
        if !velocity.is_finite() || velocity < 0.0 {
            return Err(backflow());
        }
        let face = self.state_at(
            state.density_kg_m3(),
            normal * velocity,
            state.temperature_k(),
            cancelled,
        )?;
        physical_flux(
            face,
            wave_signal.max(signal).max(velocity + face.1.sound_speed),
        )
    }
}

fn expansion(
    curve: &Isentrope<'_>,
    pressure: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ResolvedThermochemicalRealGasState, f64, f64), RealGasFlowError> {
    let p = curve.inside;
    let control = curve.control;
    if control.minimum_rarefaction_density_kg_m3 >= p.density {
        return Err(RealGasFlowError::Invalid(
            "outlet expansion requires a lower admitted density",
        ));
    }
    let maximum = det::ln(p.density / control.minimum_rarefaction_density_kg_m3);
    let mut upper = ((p.pressure - pressure) / (p.density * p.sound_speed.powi(2)))
        .min(0.5 * (p.sound_speed - p.velocity) / p.sound_speed)
        .min(maximum);
    if !positive(upper) {
        return Err(RealGasFlowError::Invalid(
            "unrepresentable outlet expansion",
        ));
    }
    // Grow a local bracket only until pressure or sonic discharge is reached.
    // Do not traverse a cold remote tail after choking has isolated the face.
    let mut bracket = None;
    let mut admitted = 0.0;
    let mut limit = maximum;
    let mut source_refusal = None;
    for _ in 0..control.compression.maximum_density_iterations {
        poll(cancelled)?;
        let state = match curve.state(upper, cancelled) {
            Ok(state) => state,
            Err(RealGasFlowError::Cancelled) => return Err(RealGasFlowError::Cancelled),
            Err(error) => {
                // A trial may overshoot a source endpoint even though sonic
                // discharge lies before it. Search only the connected segment
                // below that refused trial; never extrapolate or cross it.
                source_refusal = Some(error);
                limit = upper;
                let next = admitted + 0.5 * (upper - admitted);
                if next <= admitted || next >= upper {
                    break;
                }
                upper = next;
                continue;
            }
        };
        let velocity = p.velocity + curve.integral(upper, cancelled)?;
        if state.pressure_pa() <= pressure || velocity >= state.sound_speed_m_s() {
            bracket = Some((state, velocity));
            break;
        }
        if upper == maximum {
            break;
        }
        admitted = upper;
        let next = if 2.0 * upper < limit {
            2.0 * upper
        } else if source_refusal.is_some() {
            upper + 0.5 * (limit - upper)
        } else {
            limit
        };
        if next <= upper {
            break;
        }
        upper = next;
    }
    let Some((edge, edge_velocity)) = bracket else {
        return Err(source_refusal.unwrap_or(RealGasFlowError::Invalid(
            "outlet backpressure leaves declared isentropic path",
        )));
    };
    if edge_velocity >= edge.sound_speed_m_s() {
        let (z, sonic) = sonic_state(curve, upper, cancelled)?;
        if pressure <= sonic.pressure_pa() {
            let velocity = sonic.sound_speed_m_s();
            return Ok((
                sonic,
                velocity,
                (p.velocity + p.sound_speed).max(2.0 * velocity),
            ));
        }
        upper = z;
    }
    let state = pressure_state(curve, pressure, upper, cancelled)?;
    let z = det::ln(p.density / state.density_kg_m3());
    let velocity = p.velocity + curve.integral(z, cancelled)?;
    if velocity > state.sound_speed_m_s() + control.velocity_tolerance_m_s {
        return Err(RealGasFlowError::Invalid(
            "nonclassical real-gas outlet expansion",
        ));
    }
    let signal = (p.velocity + p.sound_speed).max(velocity + state.sound_speed_m_s());
    Ok((state, velocity, signal))
}

fn pressure_state(
    curve: &Isentrope<'_>,
    pressure: f64,
    mut high: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ResolvedThermochemicalRealGasState, RealGasFlowError> {
    let mut low = 0.0;
    let mut trial = 0.5 * high;
    for _ in 0..curve.control.compression.maximum_density_iterations {
        poll(cancelled)?;
        let state = curve.state(trial, cancelled)?;
        let residual = state.pressure_pa() - pressure;
        if residual.abs() <= curve.control.compression.momentum_tolerance_pa {
            return Ok(state);
        }
        if residual > 0.0 {
            low = trial;
        } else {
            high = trial;
        }
        // dp/dln(rho)|s = rho*a² from the SAME original EOS.
        let newton = trial + residual / (state.density_kg_m3() * state.sound_speed_m_s().powi(2));
        let next = if newton > low && newton < high {
            newton
        } else {
            low + 0.5 * (high - low)
        };
        if next == trial {
            break;
        }
        trial = next;
    }
    Err(RealGasFlowError::Invalid(
        "real-gas outlet pressure root exhausted work bound",
    ))
}

fn sonic_state(
    curve: &Isentrope<'_>,
    mut high: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, ResolvedThermochemicalRealGasState), RealGasFlowError> {
    let mut low = 0.0;
    for _ in 0..curve.control.compression.maximum_density_iterations {
        poll(cancelled)?;
        let trial = low + 0.5 * (high - low);
        if trial <= low || trial >= high {
            break;
        }
        let state = curve.state(trial, cancelled)?;
        let residual =
            curve.inside.velocity + curve.integral(trial, cancelled)? - state.sound_speed_m_s();
        if residual.abs() + curve.control.integral_tolerance_m_s
            <= curve.control.velocity_tolerance_m_s
        {
            return Ok((trial, state));
        }
        if residual < 0.0 {
            low = trial;
        } else {
            high = trial;
        }
    }
    Err(RealGasFlowError::Invalid(
        "real-gas outlet sonic root exhausted work bound",
    ))
}

fn physical_flux(
    (q, p): (Conserved, Primitive),
    signal_speed: f64,
) -> Result<FaceFlux, RealGasFlowError> {
    let flux = FaceFlux {
        mass: q.momentum,
        momentum: q.momentum * p.velocity + p.pressure,
        energy: (q.energy + p.pressure) * p.velocity,
        signal_speed,
        used_hll_fallback: false,
    };
    if !positive(signal_speed) || flux.array().iter().any(|value| !value.is_finite()) {
        return Err(RealGasFlowError::Invalid(
            "unrepresentable real-gas outlet flux",
        ));
    }
    Ok(flux)
}

fn backflow() -> RealGasFlowError {
    RealGasFlowError::Invalid("real-gas outlet backflow requires an exterior thermodynamic state")
}
