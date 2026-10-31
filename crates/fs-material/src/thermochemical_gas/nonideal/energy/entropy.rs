//! Fixed-composition entropy changes from the original NASA-9 and PR owners.
//!
//! The NASA entropy retains its pressure and mixing conventions; the PR
//! departure is derived from the same Helmholtz potential as pressure, energy
//! and Cv. At prescribed density, the exact thermal tangent is Cv/T. A zero
//! entropy increment therefore follows the homogeneous reversible adiabat
//! without imposing an ideal-gas gamma or subtracting formation references.
//!
//! The caller admits a continuous mechanically and calorically stable branch
//! over the temperature bracket. Local stability checks do not perform a
//! two-phase flash, select equilibrium phase fractions or supply latent heat.

use fs_matdb::QueryPoint;
use fs_qty::{Density, Pressure, Temperature};

use super::{ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel};
use crate::MaterialError;
use crate::thermochemical_gas::energy::entropy::{
    GasEntropyError, GasEntropySolve, map_context, map_root,
};
use crate::thermochemical_gas::energy::{poll, root, updated_point};
use crate::thermochemical_gas::{coordinates, failure, invalid, positive};

/// Complete accepted real-fluid state with its specific entropy residual.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasEntropySolution {
    state: ResolvedThermochemicalRealGasState,
    residual_j_kg_k: f64,
    iterations: usize,
}
impl RealGasEntropySolution {
    /// The same resolved material for thermal, mechanical and acoustic use.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalRealGasState {
        &self.state
    }
    /// Move the accepted state into an atomic consumer update without cloning.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalRealGasState {
        self.state
    }
    /// Actual entropy change minus the requested value [J/(kg K)].
    #[must_use]
    pub const fn residual_j_kg_k(&self) -> f64 {
        self.residual_j_kg_k
    }
    /// Bounded interior iterations, excluding endpoint and reference queries.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

struct EntropySample {
    pressure_pa: f64,
    residual_entropy_j_kg_k: f64,
    cv_j_kg_k: f64,
}

impl ResolvedThermochemicalRealGasState {
    /// Total specific entropy [J/(kg K)] in the original species reference.
    /// The ideal-reference source is evaluated at the actual EOS pressure, so
    /// its departure must use the same-T,p convention, including R ln Z.
    pub fn specific_entropy_j_kg_k(&self) -> Result<f64, MaterialError> {
        let molar_mass = self
            .source()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
        let entropy = self.source().specific_entropy_j_kg_k()?
            + self.root().residual_entropy_j_mol_k() / molar_mass;
        if !entropy.is_finite() {
            return Err(invalid("mass-specific real-gas entropy is unrepresentable"));
        }
        Ok(entropy)
    }
}

impl ThermochemicalPengRobinsonModel {
    fn entropy_at_density(&self, density: f64, t: f64) -> Result<EntropySample, MaterialError> {
        let direct = self
            .mechanics()
            .evaluate_density_temperature(Temperature::new(t), Density::new(density))
            .map_err(|error| failure("real-gas entropy density mechanics", error))?;
        let properties = self
            .source()
            .thermodynamics()
            .evaluate_properties(Temperature::new(t), direct.pressure())
            .map_err(|error| failure("real-gas reference entropy calorics", error))?;
        let molar_mass = properties.molar_mass().value();
        let residual_entropy_j_kg_k = direct.root().residual_entropy_j_mol_k() / molar_mass;
        let cv_j_kg_k = properties.mass_specific_cv().value()
            + direct.root().residual_cv_j_mol_k() / molar_mass;
        if !residual_entropy_j_kg_k.is_finite() || !positive(cv_j_kg_k) {
            return Err(invalid(
                "finite real-gas entropy departure and positive complete Cv required",
            ));
        }
        Ok(EntropySample {
            pressure_pa: direct.pressure().value(),
            residual_entropy_j_kg_k,
            cv_j_kg_k,
        })
    }

    fn entropy_increment(
        &self,
        t0: f64,
        p0: f64,
        residual0: f64,
        t1: f64,
        p1: f64,
        residual1: f64,
    ) -> Result<f64, MaterialError> {
        // The original-law difference cancels species entropy references and
        // ideal mixing BEFORE rounding, even when absolute entropy is huge.
        let ideal = self
            .source()
            .thermodynamics()
            .specific_entropy_change(
                Temperature::new(t0),
                Pressure::new(p0),
                Temperature::new(t1),
                Pressure::new(p1),
            )
            .map_err(|error| failure("real-gas reference entropy increment", error))?
            .value();
        let change = ideal + (residual1 - residual0);
        if !change.is_finite() {
            return Err(invalid("real-gas entropy increment is unrepresentable"));
        }
        Ok(change)
    }

    /// Specific entropy change [J/(kg K)] between explicit density/temperature
    /// endpoints at the retained composition. Uses the original source laws
    /// and both PR departures. This constitutive query supplies neither an
    /// entropy production model nor an assertion about the intervening path.
    /// Transport is not queried; its domain is checked on final resolution.
    pub fn specific_entropy_change_between_densities(
        &self,
        initial_density_kg_m3: f64,
        initial_temperature_k: f64,
        final_density_kg_m3: f64,
        final_temperature_k: f64,
    ) -> Result<f64, MaterialError> {
        let initial = self.entropy_at_density(initial_density_kg_m3, initial_temperature_k)?;
        let final_state = self.entropy_at_density(final_density_kg_m3, final_temperature_k)?;
        self.entropy_increment(
            initial_temperature_k,
            initial.pressure_pa,
            initial.residual_entropy_j_kg_k,
            final_temperature_k,
            final_state.pressure_pa,
            final_state.residual_entropy_j_kg_k,
        )
    }

    /// Resolve a new density and entropy change from the original state.
    /// `context` supplies the reference temperature and retained material axes;
    /// its pressure is replaced by the EOS pressure at the initial density.
    /// Explicit densities avoid choosing a different branch from a PT cubic.
    /// Zero change at different densities computes an isentropic endpoint.
    ///
    /// Each bounded root evaluation uses direct density mechanics and original
    /// NASA calorics. Transport and the PT root inventory resolve only for the
    /// initial and final full states. Final density is retained bit for bit;
    /// cancellation or any source refusal returns no candidate state.
    pub fn reconstruct_density_entropy_change(
        &self,
        context: &QueryPoint,
        initial_density_kg_m3: f64,
        density_kg_m3: f64,
        entropy_change_j_kg_k: f64,
        controls: GasEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasEntropySolution, GasEntropyError> {
        controls.validate()?;
        poll(&mut cancelled).map_err(map_context)?;
        let (initial_temperature, _) = coordinates(context)?;
        let initial_t = initial_temperature.value();
        if !positive(initial_density_kg_m3)
            || !positive(density_kg_m3)
            || !entropy_change_j_kg_k.is_finite()
            || initial_t < controls.minimum_temperature_k
            || initial_t > controls.maximum_temperature_k
        {
            return Err(GasEntropyError::Invalid(
                "real-gas entropy closure needs positive densities, finite entropy change and bracketed reference temperature",
            ));
        }
        let initial = self.evaluate_density_temperature(context, initial_density_kg_m3)?;
        let initial_p = initial.pressure_pa();
        let residual_initial =
            initial.root().residual_entropy_j_mol_k() / self.mechanics().molar_mass().value();
        if entropy_change_j_kg_k == 0.0 && density_kg_m3 == initial_density_kg_m3 {
            poll(&mut cancelled).map_err(map_context)?;
            return Ok(RealGasEntropySolution {
                state: initial,
                residual_j_kg_k: 0.0,
                iterations: 0,
            });
        }
        let tolerance = if entropy_change_j_kg_k == 0.0 {
            controls.entropy_tolerance_j_kg_k
        } else {
            controls
                .entropy_tolerance_j_kg_k
                .min(0.5 * entropy_change_j_kg_k.abs())
        };
        if !positive(tolerance) {
            return Err(GasEntropyError::Invalid(
                "real-gas entropy increment is below represented solve resolution",
            ));
        }
        let (lo, hi) = if density_kg_m3 != initial_density_kg_m3 {
            (
                controls.minimum_temperature_k,
                controls.maximum_temperature_k,
            )
        } else if entropy_change_j_kg_k > 0.0 {
            (initial_t, controls.maximum_temperature_k)
        } else {
            (controls.minimum_temperature_k, initial_t)
        };
        let solved = root::solve(
            lo,
            hi,
            tolerance,
            controls.maximum_iterations,
            |t| {
                let sample = self.entropy_at_density(density_kg_m3, t)?;
                Ok::<_, GasEntropyError>(root::Sample {
                    residual: self.entropy_increment(
                        initial_t,
                        initial_p,
                        residual_initial,
                        t,
                        sample.pressure_pa,
                        sample.residual_entropy_j_kg_k,
                    )? - entropy_change_j_kg_k,
                    tangent: sample.cv_j_kg_k / t,
                    pressure_pa: sample.pressure_pa,
                })
            },
            &mut cancelled,
        )
        .map_err(map_root)?;
        poll(&mut cancelled).map_err(map_context)?;
        let point = updated_point(context, solved.temperature_k, solved.sample.pressure_pa)
            .map_err(map_context)?;
        let state = self.evaluate_density_temperature(&point, density_kg_m3)?;
        let change = self.entropy_increment(
            initial_t,
            initial_p,
            residual_initial,
            state.temperature_k(),
            state.pressure_pa(),
            state.root().residual_entropy_j_mol_k() / self.mechanics().molar_mass().value(),
        )?;
        let residual = change - entropy_change_j_kg_k;
        if !residual.is_finite()
            || residual.abs() > tolerance
            || (entropy_change_j_kg_k != 0.0
                && (change == 0.0
                    || change.is_sign_positive() != entropy_change_j_kg_k.is_sign_positive()
                    || (density_kg_m3 == initial_density_kg_m3
                        && state.temperature_k() == initial_t)))
        {
            return Err(GasEntropyError::NotConverged {
                iterations: solved.iterations,
                residual_j_kg_k: residual,
            });
        }
        if state.density_kg_m3().to_bits() != density_kg_m3.to_bits() {
            return Err(GasEntropyError::Invalid(
                "real-gas entropy reconstruction changed conserved density",
            ));
        }
        poll(&mut cancelled).map_err(map_context)?;
        Ok(RealGasEntropySolution {
            state,
            residual_j_kg_k: residual,
            iterations: solved.iterations,
        })
    }
}

#[cfg(test)]
mod tests;
