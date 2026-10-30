//! Quasistatic adiabatic motion of one perfectly mixed gas boundary.
//!
//! This owner changes a compact plenum's physical volume while preserving every
//! original species mass. The endpoint lies on the original fixed-composition
//! ideal-gas isentrope, using NASA-9 calorics through the same source owner as
//! shock/characteristic flow. Boundary work is the resulting change in gas
//! internal energy; it is never injected as heat or approximated by p*Delta V.
//!
//! The operation is intentionally clock-free. It is a thermodynamic endpoint
//! for a piston/cylinder or moving-network driver, not a finite-speed moving-wall
//! Riemann solver. The caller owns piston inertia, prescribed motion, heat leaks,
//! and splitting with branch flow. A refused endpoint publishes nothing.

use super::Plenum;
use crate::dynamics::wall::shared::CaloricLaw;
use crate::multicomponent::wall::LocalCaloric;
use crate::multicomponent::{Error, MixtureGas, allocated, poll, positive};

/// Conservative receipt for one accepted reversible adiabatic volume change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingVolumeReport {
    /// Initial physical volume [m^3].
    pub initial_volume_m3: f64,
    /// Accepted physical volume [m^3].
    pub final_volume_m3: f64,
    /// Initial source-owned pressure [Pa].
    pub initial_pressure_pa: f64,
    /// Final source-owned pressure [Pa].
    pub final_pressure_pa: f64,
    /// Initial source-owned temperature [K].
    pub initial_temperature_k: f64,
    /// Final source-owned temperature [K].
    pub final_temperature_k: f64,
    /// Mechanical work done ON the gas [J]. Positive compression; negative expansion.
    pub boundary_work_on_gas_j: f64,
    /// Change in physical internal energy [J], evaluated as an original-source
    /// energy difference so large formation references cannot erase small work.
    pub internal_energy_change_j: f64,
    /// Actual endpoint entropy change of the closed gas [J/K].
    pub entropy_change_j_k: f64,
    /// Largest canonical species-mass defect between the two volumes [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Referenced-energy change minus source internal-energy change [J]. Species
    /// references cancel exactly in the continuum because every species mass is fixed.
    pub energy_gauge_residual_j: f64,
    /// Safeguarded temperature iterations used by the endpoint solve.
    pub iterations: usize,
}

impl Plenum {
    /// Return a reversible adiabatic candidate at \`final_volume_m3\`.
    ///
    /// Composition and every species mass are held fixed. The source-owned
    /// caloric law solves \`Delta s(T,p)=0\` with \`p = rho R_mix T\` at the new
    /// density. Safeguarded Newton/bisection is bounded by the gas's existing
    /// caloric-domain iteration budget and cannot leave its temperature interval.
    ///
    /// No clock is advanced and \`self\` is never modified. This is therefore a
    /// suitable unpublished thermodynamic trial for a later piston transaction.
    pub fn quasistatic_adiabatic_to(
        &self,
        gas: &MixtureGas,
        final_volume_m3: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, MovingVolumeReport), Error> {
        poll(&mut cancelled)?;
        if !positive(final_volume_m3)
            || self.state.partial_densities.len() != gas.species_count()
            || self.state.momentum != 0.0
        {
            return Err(Error::Invalid(
                "moving plenum needs positive volume, matching species axis and zero bulk momentum",
            ));
        }
        let initial = gas.primitive(self.state(), &mut cancelled)?;
        let initial_volume = self.volume_m3;
        let law = LocalCaloric::new(gas, &self.state.partial_densities)?;
        let total_mass = initial.density * initial_volume;
        if !positive(total_mass) {
            return Err(Error::Invalid("moving plenum mass is unrepresentable"));
        }

        if final_volume_m3.to_bits() == initial_volume.to_bits() {
            poll(&mut cancelled)?;
            return Ok((
                self.clone(),
                MovingVolumeReport {
                    initial_volume_m3: initial_volume,
                    final_volume_m3,
                    initial_pressure_pa: initial.pressure,
                    final_pressure_pa: initial.pressure,
                    initial_temperature_k: initial.temperature,
                    final_temperature_k: initial.temperature,
                    boundary_work_on_gas_j: 0.0,
                    internal_energy_change_j: 0.0,
                    entropy_change_j_k: 0.0,
                    maximum_species_mass_residual_kg: 0.0,
                    energy_gauge_residual_j: 0.0,
                    iterations: 0,
                },
            ));
        }

        let scale = initial_volume / final_volume_m3;
        if !positive(scale) {
            return Err(Error::Invalid("moving plenum density scale"));
        }
        let mut final_densities = allocated(gas.species_count(), 0.0)?;
        let mut maximum_mass_residual = 0.0_f64;
        for (i, density) in self.state.partial_densities.iter().copied().enumerate() {
            poll(&mut cancelled)?;
            let next = density * scale;
            if !next.is_finite() || next < 0.0 {
                return Err(Error::Invalid("moving plenum species density overflow"));
            }
            final_densities[i] = next;
            let before = density * initial_volume;
            let after = next * final_volume_m3;
            if !before.is_finite() || !after.is_finite() {
                return Err(Error::Invalid("moving plenum species mass overflow"));
            }
            maximum_mass_residual = maximum_mass_residual.max((after - before).abs());
        }
        let final_density = initial.density * scale;
        let gas_constant = law.gas_constant();
        if !positive(final_density) || !positive(gas_constant) {
            return Err(Error::Invalid("moving plenum final density or gas constant"));
        }

        let domain = gas.domain();
        let entropy_tolerance = domain.absolute_energy_tolerance / domain.minimum_temperature
            + domain.relative_energy_tolerance * gas_constant
            + 128.0 * f64::EPSILON * gas_constant;
        if !positive(entropy_tolerance) {
            return Err(Error::Invalid("moving plenum entropy tolerance"));
        }
        let residual = |temperature: f64, cancelled: &mut _|
         -> Result<(f64, f64, f64), Error> {
            let pressure = final_density * gas_constant * temperature;
            if !positive(pressure) {
                return Err(Error::Invalid("moving plenum endpoint pressure"));
            }
            let entropy = law.entropy_change(
                initial.temperature,
                initial.pressure,
                temperature,
                pressure,
                cancelled,
            )?;
            // Along the target isochore, ds/dT = Cv/T because p=rho*R*T.
            let derivative = law
                .heat_capacities(temperature, pressure, cancelled)?
                .cv
                / temperature;
            if !entropy.is_finite() || !positive(derivative) {
                return Err(Error::Invalid("moving plenum entropy derivative"));
            }
            Ok((entropy, derivative, pressure))
        };

        let (mut lo, mut hi) = (domain.minimum_temperature, domain.maximum_temperature);
        let (flo, _, _) = residual(lo, &mut cancelled)?;
        let (fhi, _, _) = residual(hi, &mut cancelled)?;
        if flo > entropy_tolerance || fhi < -entropy_tolerance {
            return Err(Error::EnergyDomain);
        }

        let mut temperature = initial.temperature.clamp(lo, hi);
        let mut iterations = 0;
        let (pressure, specific_entropy_change) = loop {
            poll(&mut cancelled)?;
            let (f, derivative, pressure) = residual(temperature, &mut cancelled)?;
            if f.abs() <= entropy_tolerance {
                break (pressure, f);
            }
            if iterations >= domain.maximum_iterations {
                return Err(Error::TemperatureIterations);
            }
            iterations += 1;
            if f < 0.0 {
                lo = temperature;
            } else {
                hi = temperature;
            }
            if !(lo < hi) {
                return Err(Error::TemperatureIterations);
            }
            let trial = temperature - f / derivative;
            temperature = if trial > lo && trial < hi && trial != temperature {
                trial
            } else {
                lo.midpoint(hi)
            };
            if !(temperature > lo && temperature < hi) {
                return Err(Error::TemperatureIterations);
            }
        };

        let state = gas.conserved(&final_densities, 0.0, pressure)?;
        let candidate = Plenum::new(final_volume_m3, state)?;
        // Re-decode the rounded conservative endpoint before any publication.
        let final_primitive = gas.primitive(candidate.state(), &mut cancelled)?;

        // Use the original source's reference-invariant caloric difference for
        // physical work. Subtracting two huge absolute formation energies here
        // would destroy the very small piston work this owner must preserve.
        let specific_du = law.internal_energy_change(
            initial.temperature,
            final_primitive.temperature,
            &mut cancelled,
        )?;
        let internal_energy_change = total_mass * specific_du;
        let referenced_change =
            candidate.state.energy * final_volume_m3 - self.state.energy * initial_volume;
        let gauge_residual = referenced_change - internal_energy_change;
        let entropy_change = total_mass * specific_entropy_change;
        if ![
            internal_energy_change,
            referenced_change,
            gauge_residual,
            entropy_change,
            final_primitive.temperature,
            final_primitive.pressure,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(Error::Invalid("moving plenum endpoint balance"));
        }
        poll(&mut cancelled)?;
        Ok((
            candidate,
            MovingVolumeReport {
                initial_volume_m3: initial_volume,
                final_volume_m3,
                initial_pressure_pa: initial.pressure,
                final_pressure_pa: final_primitive.pressure,
                initial_temperature_k: initial.temperature,
                final_temperature_k: final_primitive.temperature,
                boundary_work_on_gas_j: internal_energy_change,
                internal_energy_change_j: internal_energy_change,
                entropy_change_j_k: entropy_change,
                maximum_species_mass_residual_kg: maximum_mass_residual,
                energy_gauge_residual_j: gauge_residual,
                iterations,
            },
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multicomponent::characteristic_tests::caloric_gas;

    fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * scale.max(1.0),
            "actual={actual:.15e}, expected={expected:.15e}, scale={scale:.9e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn constant_cp_compression_matches_independent_isentrope_and_round_trips() {
        let gas = caloric_gas(false, 0.0);
        let initial = Plenum::new(
            1.0,
            gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap(),
        )
        .unwrap();
        let initial_primitive = gas.primitive(initial.state(), || false).unwrap();
        let (compressed, report) = initial
            .quasistatic_adiabatic_to(&gas, 0.5, || false)
            .unwrap();
        let final_primitive = gas.primitive(compressed.state(), || false).unwrap();

        // Pure A: R=1, Cp=3.5, Cv=2.5, gamma=1.4.
        let gamma = 1.4_f64;
        let temperature = initial_primitive.temperature * 2.0_f64.powf(gamma - 1.0);
        let pressure = initial_primitive.pressure * 2.0_f64.powf(gamma);
        let work = (pressure * 0.5 - initial_primitive.pressure) / (gamma - 1.0);
        close(final_primitive.temperature, temperature, temperature, 2e-12);
        close(final_primitive.pressure, pressure, pressure, 2e-12);
        close(report.boundary_work_on_gas_j, work, work.abs(), 3e-12);
        assert_eq!(report.boundary_work_on_gas_j, report.internal_energy_change_j);
        assert!(report.boundary_work_on_gas_j > 0.0);
        assert!(report.entropy_change_j_k.abs() < 2e-12);
        assert!(report.maximum_species_mass_residual_kg < 2e-16);
        assert!(report.energy_gauge_residual_j.abs() < 2e-12);

        let (restored, expansion) = compressed
            .quasistatic_adiabatic_to(&gas, 1.0, || false)
            .unwrap();
        let restored_primitive = gas.primitive(restored.state(), || false).unwrap();
        close(
            restored_primitive.temperature,
            initial_primitive.temperature,
            initial_primitive.temperature,
            3e-12,
        );
        close(
            restored_primitive.pressure,
            initial_primitive.pressure,
            initial_primitive.pressure,
            3e-12,
        );
        for (&actual, &expected) in restored
            .state()
            .partial_densities
            .iter()
            .zip(initial.state().partial_densities)
        {
            close(actual, expected, expected.abs(), 3e-15);
        }
        close(
            report.boundary_work_on_gas_j + expansion.boundary_work_on_gas_j,
            0.0,
            work.abs(),
            5e-12,
        );
        assert!(expansion.boundary_work_on_gas_j < 0.0);
    }

    #[test]
    fn variable_calorics_and_huge_formation_gauges_do_not_change_piston_physics() {
        let mut reference: Option<(f64, f64, f64)> = None;
        for formation in [0.0, 1.0e18, -1.0e18] {
            let gas = caloric_gas(true, formation);
            let initial = Plenum::new(
                0.8,
                gas.conserved(&[0.7, 0.3], 0.0, 1.2).unwrap(),
            )
            .unwrap();
            let initial_temperature = gas.primitive(initial.state(), || false).unwrap().temperature;
            let (next, report) = initial
                .quasistatic_adiabatic_to(&gas, 0.5, || false)
                .unwrap();
            let primitive = gas.primitive(next.state(), || false).unwrap();
            assert!(primitive.temperature > initial_temperature);
            assert!(report.boundary_work_on_gas_j > 0.0);
            assert!(report.entropy_change_j_k.abs() < 5e-12);
            assert!(report.maximum_species_mass_residual_kg < 5e-16);
            assert!(report.energy_gauge_residual_j.abs() < 5e-10);
            let tuple = (
                primitive.temperature,
                primitive.pressure,
                report.boundary_work_on_gas_j,
            );
            if let Some(expected) = reference {
                close(tuple.0, expected.0, tuple.0, 2e-11);
                close(tuple.1, expected.1, tuple.1, 2e-11);
                close(tuple.2, expected.2, tuple.2.abs(), 2e-10);
            } else {
                reference = Some(tuple);
            }
        }
    }

    #[test]
    fn unsupported_compression_and_cancellation_leave_the_original_plenum_untouched() {
        let gas = caloric_gas(true, 0.0);
        let initial = Plenum::new(
            1.0,
            gas.conserved(&[0.6, 0.4], 0.0, 1.0).unwrap(),
        )
        .unwrap();
        let before = initial.clone();
        assert!(matches!(
            initial.quasistatic_adiabatic_to(&gas, 1e-12, || false),
            Err(Error::EnergyDomain)
        ));
        assert!(matches!(
            initial.quasistatic_adiabatic_to(&gas, 0.0, || false),
            Err(Error::Invalid(_))
        ));
        let mut calls = 0;
        assert_eq!(
            initial
                .quasistatic_adiabatic_to(&gas, 0.6, || {
                    calls += 1;
                    calls >= 4
                })
                .unwrap_err(),
            Error::Cancelled
        );
        assert_eq!(initial, before);
    }
}
