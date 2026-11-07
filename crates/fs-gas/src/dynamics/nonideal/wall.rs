//! Original-EOS pressure traction at a stationary impermeable wall.
//!
//! Compression uses the existing Rankine--Hugoniot piston solve. Expansion
//! integrates the native isentrope, du = -a d(rho)/rho. Neither branch fits a
//! constant gamma or takes tensile pressure from a mirrored numerical flux.

use super::{RealGas, RealGasFlowError};
use crate::dynamics::{Conserved, FaceFlux, Primitive, poll, positive};
use crate::material::nonideal_shock::{RealGasShockControl, piston_compression};
use fs_matdb::QueryPoint;
use fs_material::MaterialError;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_math::det;
use fs_qty::{Density, Velocity};

/// Bounded physical wall-wave controls on the owner's homogeneous branch.
///
/// The compression bracket is a downstream/upstream density ratio. The
/// rarefaction endpoint is an absolute density on the SAME isentrope as the
/// current boundary cell. Both entire paths must remain inside the owner's
/// density and temperature domain; unsupported paths are refused, not clipped.
/// The compression search may extend its lower ratio toward one to include
/// a weak shock connected continuously to the admitted upstream state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasWallControl {
    /// Original native compression solve, including its density-ratio bracket.
    pub compression: RealGasShockControl,
    /// Lowest density [kg/m³] admitted for the expansion integral.
    pub minimum_rarefaction_density_kg_m3: f64,
    /// Absolute wall characteristic-velocity residual [m/s]. Waves smaller
    /// than this, or below the floating-point weak-shock resolution, use the
    /// source acoustic tangent with its O(velocity²) error. Every nonzero
    /// continuation still admits its original isentropic density/temperature
    /// path. The compression pressure tolerance gates consistency with that
    /// isentrope at the linearized density, not a certified wall-wave error.
    pub velocity_tolerance_m_s: f64,
    /// Absolute quadrature error estimate budget [m/s], smaller than the
    /// characteristic tolerance. The embedded Gauss estimate is not a proof.
    pub integral_tolerance_m_s: f64,
    /// At most this many bisections of one quadrature panel, in 0..=12.
    /// Numerical scratch is a fixed stack; original material solves may allocate.
    pub maximum_quadrature_depth: usize,
}

impl RealGasWallControl {
    pub(super) fn validate(self, gas: &RealGas) -> Result<(), RealGasFlowError> {
        self.compression.validate()?;
        if !positive(self.minimum_rarefaction_density_kg_m3)
            || !gas
                .domain
                .admits_density(self.minimum_rarefaction_density_kg_m3)
            || !positive(self.velocity_tolerance_m_s)
            || !positive(self.integral_tolerance_m_s)
            || self.integral_tolerance_m_s > 0.25 * self.velocity_tolerance_m_s
            || self.maximum_quadrature_depth > 12
            || self.compression.energy.minimum_temperature_k
                < gas.domain.energy.minimum_temperature_k
            || self.compression.energy.maximum_temperature_k
                > gas.domain.energy.maximum_temperature_k
        {
            return Err(RealGasFlowError::Invalid(
                "physical wall controls exceed original real-gas domain or work/tolerance bounds",
            ));
        }
        Ok(())
    }
}

impl RealGas {
    /// Physical stationary-wall flux from the original real-gas source.
    ///
    /// `outward_normal` is -1 at a left wall and +1 at a right wall. Positive
    /// outward velocity generates a compression shock; negative velocity a
    /// rarefaction. The result has exactly zero mass and energy flux, native
    /// positive pressure traction, and a wave-speed bound for CFL on the
    /// caller-declared classical branch. In particular, the rarefaction fan
    /// must have monotone characteristic speeds so its endpoints bound the
    /// fan; material stability alone does not establish this precondition.
    /// This method does not certify that a source excludes nonclassical waves,
    /// nor does it substitute a phase transition or a vacuum state.
    pub fn reflecting_wall_flux(
        &self,
        inside: Conserved,
        outward_normal: f64,
        control: RealGasWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        control.validate(self)?;
        let primitive = self.decode(inside, &mut cancelled)?;
        self.physical_wall_flux(primitive, outward_normal, control, &mut cancelled)
    }

    pub(super) fn physical_wall_flux(
        &self,
        inside: Primitive,
        outward_normal: f64,
        control: RealGasWallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        poll(cancelled)?;
        if !matches!(outward_normal, -1.0 | 1.0) {
            return Err(RealGasFlowError::Invalid(
                "wall outward normal must be -1 or +1",
            ));
        }
        let speed = outward_normal * inside.velocity;
        let signal = inside.sound_speed + speed.abs();
        if speed == 0.0 {
            return wall_flux(inside.pressure, signal);
        }
        let context = self
            .context
            .clone()
            .with("temperature", inside.temperature)
            .map_err(|error| MaterialError::State {
                what: format!("real-gas wall temperature coordinate: {error}"),
            })?;
        let curve = Isentrope {
            gas: self,
            context,
            inside,
            control,
        };
        // A weak shock's upstream/downstream sonic interval in density ratio
        // is O(M²). Below this roundoff scale a ratio-based shock solve cannot
        // reliably distinguish its two sonic endpoints. Bridge that interval
        // with the SAME source-admitted acoustic continuation on both signs.
        let acoustic_limit = control
            .velocity_tolerance_m_s
            .max(32.0 * det::sqrt(f64::EPSILON) * inside.sound_speed);
        if speed.abs() <= acoustic_limit {
            return curve.acoustic_wall(speed, signal, cancelled);
        }
        if speed > 0.0 {
            let mut compression = control.compression;
            let weak_lower = 1.0 + (0.5 * speed) / (inside.sound_speed + speed);
            compression.compression_ratio_bracket[0] =
                compression.compression_ratio_bracket[0].min(weak_lower);
            // Resolve the O(speed²) sonic separation instead of letting a
            // coarse caloric/pressure residual obscure a weak compression.
            // Strong waves retain the caller's original tolerances exactly.
            let weak_energy_tolerance = (speed / 64.0) * speed;
            compression.energy.energy_tolerance_j_kg = compression
                .energy
                .energy_tolerance_j_kg
                .min(weak_energy_tolerance);
            compression.momentum_tolerance_pa = compression
                .momentum_tolerance_pa
                .min(inside.density * weak_energy_tolerance);
            let maximum_density = inside.density * compression.compression_ratio_bracket[1];
            if !self.domain.admits_density(maximum_density) {
                return Err(RealGasFlowError::Invalid(
                    "wall compression bracket exceeds the original density domain",
                ));
            }
            // Piston velocities are signed FROM the wall into gas. A stopped
            // wall sees upstream normal velocity -speed, and piston speed 0.
            let wave = piston_compression(
                &self.model,
                &curve.context,
                Density::new(inside.density),
                Velocity::new(-speed),
                Velocity::new(0.0),
                compression,
                &mut *cancelled,
            )?;
            let downstream = wave.shock.downstream();
            return wall_flux(
                downstream.pressure_pa(),
                signal
                    .max(downstream.sound_speed_m_s())
                    .max(wave.shock_normal_velocity.value()),
            );
        }
        curve.expanded_wall(-speed, signal, cancelled)
    }
}

fn wall_flux(pressure: f64, signal_speed: f64) -> Result<FaceFlux, RealGasFlowError> {
    if !positive(pressure) || !positive(signal_speed) {
        return Err(RealGasFlowError::Invalid(
            "unrepresentable physical real-gas wall flux",
        ));
    }
    Ok(FaceFlux {
        mass: 0.0,
        momentum: pressure,
        energy: 0.0,
        signal_speed,
        used_hll_fallback: false,
    })
}

struct Isentrope<'a> {
    gas: &'a RealGas,
    context: QueryPoint,
    inside: Primitive,
    control: RealGasWallControl,
}

impl Isentrope<'_> {
    fn acoustic_wall(
        &self,
        speed: f64,
        signal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        poll(cancelled)?;
        let expansion = -speed / self.inside.sound_speed;
        let density = self.inside.density * det::exp(-expansion);
        if !self.gas.domain.admits_density(density)
            || (speed < 0.0 && density < self.control.minimum_rarefaction_density_kg_m3)
            || (speed > 0.0
                && density
                    > self.inside.density * self.control.compression.compression_ratio_bracket[1])
        {
            return Err(RealGasFlowError::Invalid(
                "acoustic wall continuation leaves the declared density path",
            ));
        }
        let energy = self.control.compression.energy;
        let upstream = self
            .gas
            .model
            .evaluate_density_temperature(&self.context, self.inside.density)?;
        // dT/dln(rho)|s = T * (dp/dT)|rho / (rho Cv). Use the original
        // thermal derivative's sign even if an extremely small density
        // increment rounds to zero. A nonzero outward thermal continuation
        // at a source endpoint must not disappear inside a root tolerance.
        let thermal = upstream.root().pressure_derivatives().0;
        let heats = (thermal > 0.0 && speed > 0.0) || (thermal < 0.0 && speed < 0.0);
        let cools = (thermal > 0.0 && speed < 0.0) || (thermal < 0.0 && speed > 0.0);
        if (heats && self.inside.temperature >= energy.maximum_temperature_k)
            || (cools && self.inside.temperature <= energy.minimum_temperature_k)
        {
            return Err(RealGasFlowError::Invalid(
                "acoustic wall continuation leaves the declared temperature path",
            ));
        }
        // Check the actual entropy bracket before native inversion; its
        // numerical residual tolerance must not widen a physical source bound.
        let entropy_at = |temperature| {
            self.gas.model.specific_entropy_change_between_densities(
                self.inside.density,
                self.inside.temperature,
                density,
                temperature,
            )
        };
        let lower = entropy_at(energy.minimum_temperature_k)?;
        poll(cancelled)?;
        let upper = entropy_at(energy.maximum_temperature_k)?;
        if lower > 0.0 || upper < 0.0 {
            return Err(RealGasFlowError::Invalid(
                "acoustic wall isentrope is outside the original temperature bracket",
            ));
        }
        let state = self.state(expansion, cancelled)?;
        let pressure = self.inside.pressure + self.inside.density * self.inside.sound_speed * speed;
        // This compares two states at a linearized density displacement; it
        // does not certify the error against the fully resolved wall wave.
        if (state.pressure_pa() - pressure).abs() > self.control.compression.momentum_tolerance_pa {
            return Err(RealGasFlowError::Invalid(
                "acoustic wall continuation exceeds the authored pressure tolerance",
            ));
        }
        wall_flux(pressure, signal.max(state.sound_speed_m_s()))
    }

    fn state(
        &self,
        expansion: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalRealGasState, RealGasFlowError> {
        poll(cancelled)?;
        // Integrate in z = ln(rho_inside / rho). This remains well scaled over
        // many density decades and removes the explicit 1/rho integrand.
        let density = self.inside.density * det::exp(-expansion);
        let energy = self.control.compression.energy;
        let solution = self.gas.model.reconstruct_density_entropy_change(
            &self.context,
            self.inside.density,
            density,
            0.0,
            GasEntropySolve {
                minimum_temperature_k: energy.minimum_temperature_k,
                maximum_temperature_k: energy.maximum_temperature_k,
                entropy_tolerance_j_kg_k: energy.energy_tolerance_j_kg
                    / energy.maximum_temperature_k,
                maximum_iterations: energy.maximum_iterations,
            },
            &mut *cancelled,
        )?;
        Ok(solution.into_state())
    }

    fn expanded_wall(
        &self,
        speed: f64,
        signal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        if self.control.minimum_rarefaction_density_kg_m3 >= self.inside.density {
            return Err(RealGasFlowError::Invalid(
                "wall expansion requires an admitted density below the boundary cell",
            ));
        }
        let mut lo = 0.0;
        let mut hi = det::ln(self.inside.density / self.control.minimum_rarefaction_density_kg_m3);
        // Check the endpoint separately: the Gaussian nodes are interior and
        // must not silently accept an unsupported endpoint temperature.
        self.state(hi, cancelled)?;
        let available = self.integral(hi, cancelled)?;
        if available + self.control.integral_tolerance_m_s < speed {
            return Err(RealGasFlowError::Invalid(
                "wall retreat exceeds the declared homogeneous isentropic path",
            ));
        }
        let mut trial = (speed / self.inside.sound_speed).min(hi);
        for _ in 0..self.control.compression.maximum_density_iterations {
            poll(cancelled)?;
            let state = self.state(trial, cancelled)?;
            let residual = self.integral(trial, cancelled)? - speed;
            if residual.abs() + self.control.integral_tolerance_m_s
                <= self.control.velocity_tolerance_m_s
            {
                return wall_flux(state.pressure_pa(), signal.max(state.sound_speed_m_s()));
            }
            if residual < 0.0 {
                lo = trial;
            } else {
                hi = trial;
            }
            let newton = trial - residual / state.sound_speed_m_s();
            let next = if newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == trial {
                break;
            }
            trial = next;
        }
        Err(RealGasFlowError::Invalid(
            "real-gas wall characteristic root exhausted its work bound",
        ))
    }

    fn integral(
        &self,
        end: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, RealGasFlowError> {
        if end == 0.0 {
            return Ok(0.0);
        }
        #[derive(Clone, Copy, Default)]
        struct Panel {
            lo: f64,
            hi: f64,
            tolerance: f64,
            depth: usize,
        }
        // Depth-first subdivision needs at most depth+1 pending panels.
        let mut stack = [Panel::default(); 13];
        stack[0] = Panel {
            lo: 0.0,
            hi: end,
            tolerance: self.control.integral_tolerance_m_s,
            depth: 0,
        };
        let mut pending = 1;
        let mut total = 0.0;
        while pending > 0 {
            poll(cancelled)?;
            pending -= 1;
            let panel = stack[pending];
            let center = 0.5 * (panel.lo + panel.hi);
            let half = 0.5 * (panel.hi - panel.lo);
            let middle = self.state(center, cancelled)?.sound_speed_m_s();
            let sample_pair = |node: f64, cancelled: &mut _| -> Result<f64, RealGasFlowError> {
                Ok(self
                    .state(center - half * node, cancelled)?
                    .sound_speed_m_s()
                    + self
                        .state(center + half * node, cancelled)?
                        .sound_speed_m_s())
            };
            let coarse = half
                * ((8.0 / 9.0) * middle
                    + (5.0 / 9.0) * sample_pair(0.774_596_669_241_483_4, cancelled)?);
            let fine = half
                * (0.568_888_888_888_888_9 * middle
                    + 0.478_628_670_499_366_5 * sample_pair(0.538_469_310_105_683_1, cancelled)?
                    + 0.236_926_885_056_189_1 * sample_pair(0.906_179_845_938_664, cancelled)?);
            if !positive(fine) || !coarse.is_finite() {
                return Err(RealGasFlowError::Invalid(
                    "nonfinite real-gas wall characteristic integral",
                ));
            }
            if (fine - coarse).abs() <= panel.tolerance {
                total += fine;
            } else {
                if panel.depth == self.control.maximum_quadrature_depth {
                    return Err(RealGasFlowError::Invalid(
                        "real-gas wall quadrature exhausted its work bound",
                    ));
                }
                stack[pending] = Panel {
                    lo: center,
                    hi: panel.hi,
                    tolerance: 0.5 * panel.tolerance,
                    depth: panel.depth + 1,
                };
                stack[pending + 1] = Panel {
                    lo: panel.lo,
                    hi: center,
                    tolerance: 0.5 * panel.tolerance,
                    depth: panel.depth + 1,
                };
                pending += 2;
            }
        }
        Ok(total)
    }
}
