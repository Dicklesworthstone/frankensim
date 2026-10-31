//! Quasi-steady converging-nozzle flow for the source-bound Peng--Robinson parcel.
//!
//! The upstream frozen-composition real gas follows its original reversible
//! isentrope by density. Back-pressure matching and sonic choking are solved as
//! two independent scalar conditions; the first state encountered during
//! expansion sets the throat. The same NASA/PR state supplies enthalpy, density
//! and sound speed. No ideal-gas gamma, pressure law or caloric surrogate enters.
//!
//! This is a homogeneous single-phase engineering nozzle model. It assumes one
//! continuous mechanically/calorically stable PR branch over the caller's
//! density/temperature bounds. It does not perform a flash, resolve shocks,
//! model a diverging nozzle, wall friction or jet momentum in the receiving tank.

use super::{
    IsochoricRealGasParcel, RealGasMassTransferCandidate, RealGasMassTransferSolve,
};
use crate::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use crate::thermochemical_gas::energy::{
    GasEnergyError, poll,
    entropy::{GasEntropyError, GasEntropySolve},
    parcel::transfer::nozzle::GasNozzle,
};
use crate::thermochemical_gas::positive;

/// Bounded isentropic/nozzle and finite-volume transfer controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasNozzleStepSolve {
    /// Temperature/entropy bracket for every state on the homogeneous isentrope.
    pub entropy: GasEntropySolve,
    /// Endpoint energy solves and extensive conservation tolerances.
    pub transfer: RealGasMassTransferSolve,
    /// Lowest explored density divided by the accepted upstream density.
    /// Must lie strictly in (0,1).
    pub minimum_density_fraction: f64,
    /// Absolute pressure residual for an unchoked exit [Pa].
    pub pressure_tolerance_pa: f64,
    /// Absolute residual of c^2-v^2 at a sonic throat [m2/s2].
    pub sonic_tolerance_m2_s2: f64,
    /// Bisection work for each outer pressure/sonic density solve.
    pub maximum_iterations: usize,
    /// Largest admitted fraction of donor inventory transferred in one step.
    pub maximum_donor_fraction: f64,
}
impl RealGasNozzleStepSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.entropy.validate().map_err(map_entropy)?;
        self.transfer.validate()?;
        if !positive(self.minimum_density_fraction)
            || self.minimum_density_fraction >= 1.0
            || !positive(self.pressure_tolerance_pa)
            || !positive(self.sonic_tolerance_m2_s2)
            || !(1..=128).contains(&self.maximum_iterations)
            || !positive(self.maximum_donor_fraction)
            || self.maximum_donor_fraction >= 1.0
        {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle requires bounded density, pressure/sonic tolerances, work and donor fraction",
            ));
        }
        Ok(())
    }
}

/// One source-resolved real-gas throat before inventory transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleFlow {
    throat: ResolvedThermochemicalRealGasState,
    /// Positive mass rate from donor toward receiver [kg/s].
    pub mass_flow_kg_s: f64,
    /// Isentropic exit velocity from upstream stagnation enthalpy [m/s].
    pub velocity_m_s: f64,
    /// Velocity divided by the same PR throat sound speed.
    pub mach_number: f64,
    /// True when the sonic state is encountered before the back-pressure state.
    pub choked: bool,
    /// Original NASA/PR upstream total enthalpy [J/kg].
    pub upstream_enthalpy_j_kg: f64,
    /// Actual source entropy defect of the selected throat [J/(kg K)].
    pub entropy_residual_j_kg_k: f64,
    /// Sum of inner entropy and outer density iterations.
    pub iterations: usize,
}
impl RealGasNozzleFlow {
    #[must_use]
    pub const fn throat(&self) -> &ResolvedThermochemicalRealGasState { &self.throat }
}

#[derive(Clone)]
struct IsentropicSample {
    state: ResolvedThermochemicalRealGasState,
    entropy_residual: f64,
    inner_iterations: usize,
    velocity_squared: f64,
}

fn map_entropy(error: GasEntropyError) -> GasEnergyError {
    match error {
        GasEntropyError::Material(error) => GasEnergyError::Material(error),
        GasEntropyError::Cancelled => GasEnergyError::Cancelled,
        GasEntropyError::Invalid(message) => GasEnergyError::Invalid(message),
        GasEntropyError::OutsideBracket { .. } => GasEnergyError::Invalid(
            "real-gas nozzle isentrope leaves the admitted temperature bracket",
        ),
        GasEntropyError::NotConverged { .. } => GasEnergyError::Invalid(
            "real-gas nozzle entropy reconstruction exhausted its work budget",
        ),
    }
}

impl IsochoricRealGasParcel {
    fn nozzle_isentropic_sample(
        &self,
        density: f64,
        controls: GasEntropySolve,
        h0: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IsentropicSample, GasEnergyError> {
        let solved = self.model().reconstruct_density_entropy_change(
            self.point(), self.density_kg_m3(), density, 0.0, controls, &mut *cancelled,
        ).map_err(map_entropy)?;
        let inner_iterations = solved.iterations();
        let entropy_residual = solved.residual_j_kg_k();
        let state = solved.into_state();
        let velocity_squared = 2.0 * (h0 - state.specific_enthalpy_j_kg());
        if !velocity_squared.is_finite() || velocity_squared < 0.0 {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle isentrope produced negative or nonfinite kinetic energy",
            ));
        }
        Ok(IsentropicSample {
            state, entropy_residual, inner_iterations, velocity_squared,
        })
    }

    fn nozzle_density_root(
        &self,
        mut lo: f64,
        mut hi: f64,
        tolerance: f64,
        controls: RealGasNozzleStepSolve,
        h0: f64,
        mut residual: impl FnMut(&IsentropicSample) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(IsentropicSample, usize), GasEnergyError> {
        let mut low = self.nozzle_isentropic_sample(lo, controls.entropy, h0, cancelled)?;
        let mut high = self.nozzle_isentropic_sample(hi, controls.entropy, h0, cancelled)?;
        let mut fl = residual(&low);
        let fh = residual(&high);
        if !fl.is_finite() || !fh.is_finite() || fl > 0.0 || fh < 0.0 {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle density root is not bracketed",
            ));
        }
        let mut outer = 0;
        loop {
            if fl.abs() <= tolerance {
                return Ok((low, outer));
            }
            if outer == controls.maximum_iterations {
                return Err(GasEnergyError::NotConverged {
                    iterations: outer,
                    residual_j_kg: fl,
                });
            }
            let mid = lo + 0.5 * (hi - lo);
            if !(mid > lo && mid < hi) {
                return Err(GasEnergyError::NotConverged {
                    iterations: outer,
                    residual_j_kg: fl,
                });
            }
            poll(cancelled)?;
            let sample = self.nozzle_isentropic_sample(mid, controls.entropy, h0, cancelled)?;
            let fm = residual(&sample);
            if !fm.is_finite() {
                return Err(GasEnergyError::Invalid(
                    "real-gas nozzle density residual is nonfinite",
                ));
            }
            outer += 1;
            if fm <= 0.0 {
                lo = mid;
                low = sample;
                fl = fm;
            } else {
                hi = mid;
                high = sample;
            }
            if fm.abs() <= tolerance {
                return Ok((if fm <= 0.0 { low } else { high }, outer));
            }
        }
    }

    /// Resolve a converging-nozzle throat against a lower back pressure.
    ///
    /// Density decreases monotonically from the accepted donor. The pressure
    /// root enforces p*=p_back. The independent sonic root enforces v*=c* with
    /// v* from the original total enthalpy. The larger-density root is reached
    /// first and therefore owns the exit. Vacuum can only select the sonic root.
    pub fn nozzle_flow_to(
        &self,
        back_pressure_pa: f64,
        nozzle: GasNozzle,
        controls: RealGasNozzleStepSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasNozzleFlow, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !nozzle.area_m2.is_finite() || nozzle.area_m2 < 0.0
            || !nozzle.discharge_coefficient.is_finite()
            || !(0.0..=1.0).contains(&nozzle.discharge_coefficient)
        {
            return Err(GasEnergyError::Invalid(
                "finite nonnegative nozzle area and discharge coefficient in [0,1] required",
            ));
        }
        let p0 = self.state().pressure_pa();
        let rho0 = self.density_kg_m3();
        let h0 = self.state().specific_enthalpy_j_kg();
        if !back_pressure_pa.is_finite() || back_pressure_pa < 0.0 || back_pressure_pa > p0 {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle requires finite back pressure between vacuum and upstream pressure",
            ));
        }
        if back_pressure_pa == p0 || nozzle.area_m2 == 0.0 || nozzle.discharge_coefficient == 0.0 {
            poll(&mut cancelled)?;
            return Ok(RealGasNozzleFlow {
                throat: self.state().clone(), mass_flow_kg_s: 0.0,
                velocity_m_s: 0.0, mach_number: 0.0, choked: false,
                upstream_enthalpy_j_kg: h0, entropy_residual_j_kg_k: 0.0,
                iterations: 0,
            });
        }
        let rho_min = rho0 * controls.minimum_density_fraction;
        if !positive(rho_min) || rho_min == rho0 {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle minimum density is unrepresentable",
            ));
        }
        let low = self.nozzle_isentropic_sample(
            rho_min, controls.entropy, h0, &mut cancelled,
        )?;
        let pressure_root = if back_pressure_pa > 0.0
            && low.state.pressure_pa() <= back_pressure_pa
        {
            Some(self.nozzle_density_root(
                rho_min, rho0, controls.pressure_tolerance_pa, controls, h0,
                |sample| sample.state.pressure_pa() - back_pressure_pa,
                &mut cancelled,
            )?)
        } else {
            None
        };
        let sonic_low = low.state.sound_speed_m_s().powi(2) - low.velocity_squared;
        let sonic_root = if sonic_low <= 0.0 {
            Some(self.nozzle_density_root(
                rho_min, rho0, controls.sonic_tolerance_m2_s2, controls, h0,
                |sample| sample.state.sound_speed_m_s().powi(2) - sample.velocity_squared,
                &mut cancelled,
            )?)
        } else {
            None
        };
        let (selected, outer, choked) = match (pressure_root, sonic_root) {
            (Some((pressure, pi)), Some((sonic, si))) => {
                if sonic.state.density_kg_m3() >= pressure.state.density_kg_m3() {
                    (sonic, si, true)
                } else {
                    (pressure, pi, false)
                }
            }
            (Some((pressure, pi)), None) => (pressure, pi, false),
            (None, Some((sonic, si))) => (sonic, si, true),
            (None, None) => {
                return Err(GasEnergyError::Invalid(
                    "real-gas nozzle density bracket reaches neither back pressure nor sonic choking",
                ));
            }
        };
        let velocity = selected.velocity_squared.sqrt();
        let c = selected.state.sound_speed_m_s();
        let mass_flow = nozzle.discharge_coefficient * nozzle.area_m2
            * selected.state.density_kg_m3() * velocity;
        let mach = velocity / c;
        if !positive(mass_flow) || !positive(velocity) || !positive(mach)
            || (choked && (mach - 1.0).abs() > 1e-5)
            || (!choked && mach > 1.0 + 1e-5)
        {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle throat failed positive-flow or Mach admission",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RealGasNozzleFlow {
            throat: selected.state,
            mass_flow_kg_s: mass_flow,
            velocity_m_s: velocity,
            mach_number: mach,
            choked,
            upstream_enthalpy_j_kg: h0,
            entropy_residual_j_kg_k: selected.entropy_residual,
            iterations: selected.inner_iterations + outer,
        })
    }

    /// Resolve pressure-driven flow and conservatively transfer its mass over one
    /// explicit physical interval. Direction is upstream-to-downstream; reverse
    /// flow requires swapping the parcels. A pressure crossing after the fixed
    /// upstream-rate step refuses rather than clipping the transferred mass.
    pub fn through_nozzle_to(
        &self,
        receiver: &Self,
        nozzle: GasNozzle,
        duration_s: f64,
        controls: RealGasNozzleStepSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(RealGasMassTransferCandidate, RealGasNozzleFlow), GasEnergyError> {
        controls.validate()?;
        if !positive(duration_s) {
            return Err(GasEnergyError::Invalid(
                "positive real-gas nozzle duration required",
            ));
        }
        let flow = self.nozzle_flow_to(
            receiver.state().pressure_pa(), nozzle, controls, &mut cancelled,
        )?;
        let mass = flow.mass_flow_kg_s * duration_s;
        if !mass.is_finite()
            || (flow.mass_flow_kg_s > 0.0 && mass == 0.0)
            || mass / self.mass_kg() > controls.maximum_donor_fraction
        {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle timestep exceeds donor inventory or mass resolution",
            ));
        }
        let candidate = self.transfer_mass_to(receiver, mass, controls.transfer, &mut cancelled)?;
        if candidate.donor().state().pressure_pa() < candidate.receiver().state().pressure_pa() {
            return Err(GasEnergyError::Invalid(
                "real-gas nozzle timestep crosses pressure equilibrium; reduce the interval",
            ));
        }
        poll(&mut cancelled)?;
        Ok((candidate, flow))
    }
}

#[cfg(test)]
mod tests;
