//! Quasistatic adiabatic motion of one perfectly mixed gas boundary.
//!
//! This owner changes a compact plenum's physical volume while preserving every
//! original species mass. The endpoint lies on the original fixed-composition
//! ideal-gas isentrope, using NASA-9 calorics through the same source owner as
//! shock/characteristic flow. Boundary work is the resulting change in gas
//! internal energy; it is never injected as heat or approximated by p*Delta V.
//!
//! The endpoint operation is intentionally clock-free. A running network can
//! accept it between transport steps, keeping its branches and clock unchanged.
//! For simultaneous finite-time port transport and volume motion, use
//! [`MixtureNetwork1d::step_with_volume_rates`], which integrates pressure work
//! and extensive storage in the transport owner's RK stages. Both paths assume
//! a perfectly mixed chamber; neither resolves spatial piston shocks or piston
//! inertia. Refusal publishes neither a volume change nor a new report.

use super::{Accumulator, MixtureNetwork1d, NetworkTotals, Plenum};
use crate::dynamics::wall::shared::CaloricLaw;
use crate::multicomponent::evolve::totals;
use crate::multicomponent::wall::LocalCaloric;
use crate::multicomponent::{Error, MixtureGas, allocated, poll, positive};
use fs_qty::Time;

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

/// Accepted compression/expansion of one connected network chamber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkVolumeReport {
    /// Index in the network's retained plenum order.
    pub plenum: usize,
    /// Accepted network clock; this thermodynamic operation advances no time.
    pub time: Time,
    /// Original-source thermodynamics and mechanical work of the changed volume.
    pub volume: MovingVolumeReport,
    /// Combined branch and plenum storage before the operation.
    pub before: NetworkTotals,
    /// Combined storage after the operation, in the same species energy gauges.
    pub after: NetworkTotals,
    /// Species inventory changes, and `after.energy - before.energy - work`.
    /// Energy subtraction has the rounding scale of the complete network;
    /// `volume.energy_gauge_residual_j` retains the local work-balance diagnostic.
    pub conservation_residual: NetworkTotals,
}

impl MixtureNetwork1d {
    /// Compress or expand one chamber reversibly between accepted flow steps.
    ///
    /// The chamber keeps every species mass and uses the original caloric
    /// source's isentrope. All branch fields, other plenums, connections and
    /// clocks are unchanged. The next transport step uses the accepted volume,
    /// pressure and temperature, including the updated plenum CFL constraint.
    ///
    /// Mechanical work is returned separately from flow and heat exchanges.
    /// After an actual volume change, the preceding transport report and port
    /// exchanges are invalidated: they no longer describe current storage.
    /// An exact no-op, refusal or cancellation preserves those reports.
    ///
    /// This clock-free source split assumes a uniform, quasistatically evolving
    /// chamber while ports are held closed for the volume change. The caller
    /// owns scheduling with port transport; simultaneous finite-speed piston
    /// motion or a temporal order for that composition is not claimed.
    pub fn quasistatic_adiabatic_plenum_to(
        &mut self,
        plenum: usize,
        final_volume_m3: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkVolumeReport, Error> {
        poll(&mut cancelled)?;
        let original = self
            .plenums
            .get(plenum)
            .ok_or(Error::Invalid("moving plenum index outside network"))?;
        let (candidate, volume) =
            original.quasistatic_adiabatic_to(self.gas(), final_volume_m3, &mut cancelled)?;
        let n = self.gas().species_count();
        let mut before = Accumulator::new(n);
        let mut after = Accumulator::new(n);
        // One branch-field pass admits both extensive totals. Checking only the
        // local chamber would allow a finite candidate to overflow the complete
        // network's previously admitted storage invariant.
        for branch in &self.branches {
            let value = totals(&branch.states, branch.mesh.volumes(), n, &mut cancelled)?;
            before.add(&value.species_masses[..n], value.energy)?;
            after.add(&value.species_masses[..n], value.energy)?;
        }
        for (index, existing) in self.plenums.iter().enumerate() {
            poll(&mut cancelled)?;
            before.add_cell(existing.state(), existing.volume_m3)?;
            let next = if index == plenum {
                &candidate
            } else {
                existing
            };
            after.add_cell(next.state(), next.volume_m3)?;
        }
        let before = before.totals();
        let after = after.totals();
        let mut conservation_residual = after;
        for k in 0..n {
            conservation_residual.species_masses[k] -= before.species_masses[k];
        }
        conservation_residual.energy =
            (after.energy - before.energy) - volume.boundary_work_on_gas_j;
        if !conservation_residual.energy.is_finite() {
            return Err(Error::Invalid("moving network energy balance overflow"));
        }
        let report = NetworkVolumeReport {
            plenum,
            time: self.time(),
            volume,
            before,
            after,
            conservation_residual,
        };
        // Only this publication point changes accepted state or report validity.
        poll(&mut cancelled)?;
        if candidate.volume_m3.to_bits() != original.volume_m3.to_bits() {
            self.plenums[plenum] = candidate;
            self.last_report = None;
        }
        Ok(report)
    }
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
            return Err(Error::Invalid(
                "moving plenum final density or gas constant",
            ));
        }

        let domain = gas.domain();
        let entropy_tolerance = domain.absolute_energy_tolerance / domain.minimum_temperature
            + domain.relative_energy_tolerance * gas_constant
            + 128.0 * f64::EPSILON * gas_constant;
        if !positive(entropy_tolerance) {
            return Err(Error::Invalid("moving plenum entropy tolerance"));
        }
        let residual = |temperature: f64, cancelled: &mut _| -> Result<(f64, f64, f64), Error> {
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
            let derivative =
                law.heat_capacities(temperature, pressure, cancelled)?.cv / temperature;
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
    use crate::multicomponent::network::{BranchEnd, PlenumConnection};
    use crate::multicomponent::{Boundaries, Boundary, Conserved, MixtureEuler1d, StepControl};
    use fs_qty::Length;

    fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * scale.max(1.0),
            "actual={actual:.15e}, expected={expected:.15e}, scale={scale:.9e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn constant_cp_compression_matches_independent_isentrope_and_round_trips() {
        let gas = caloric_gas(false, 0.0);
        let initial = Plenum::new(1.0, gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap()).unwrap();
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
        assert_eq!(
            report.boundary_work_on_gas_j,
            report.internal_energy_change_j
        );
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
            let initial = Plenum::new(0.8, gas.conserved(&[0.7, 0.3], 0.0, 1.2).unwrap()).unwrap();
            let initial_temperature = gas
                .primitive(initial.state(), || false)
                .unwrap()
                .temperature;
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
        let initial = Plenum::new(1.0, gas.conserved(&[0.6, 0.4], 0.0, 1.0).unwrap()).unwrap();
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

    fn connected_chamber(gas: &MixtureGas) -> MixtureNetwork1d {
        let branch = MixtureEuler1d::new(
            gas.clone(),
            Length::new(0.05),
            0.01,
            vec![gas.conserved(&[0.2, 0.8], 0.0, 1.0).unwrap(); 8],
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            || false,
        )
        .unwrap();
        MixtureNetwork1d::new(
            vec![branch],
            vec![
                Plenum::new(0.004, gas.conserved(&[0.8, 0.2], 0.0, 1.0).unwrap()).unwrap(),
                Plenum::new(0.006, gas.conserved(&[0.0, 1.0], 0.0, 1.0).unwrap()).unwrap(),
            ],
            vec![PlenumConnection {
                branch: 0,
                end: BranchEnd::Left,
                plenum: 0,
            }],
            || false,
        )
        .unwrap()
    }

    fn branch_states(network: &MixtureNetwork1d) -> Vec<Conserved> {
        network
            .branches()
            .iter()
            .flat_map(|branch| {
                (0..branch.cell_count()).map(move |i| {
                    let q = branch.cell(i).unwrap();
                    Conserved {
                        partial_densities: q.partial_densities.to_vec(),
                        momentum: q.momentum,
                        energy: q.energy,
                    }
                })
            })
            .collect()
    }

    // Independent physical storage integration, including signed formation
    // energies. This does not use either network totals or report residuals.
    fn physical_storage(network: &MixtureNetwork1d) -> ([f64; 2], f64) {
        let mut masses = [0.0; 2];
        let mut energy = 0.0;
        let mut add = |q: crate::multicomponent::Cell<'_>, volume: f64| {
            for k in 0..2 {
                masses[k] += volume * q.partial_densities[k];
            }
            energy += volume * network.gas().absolute_energy_density(q).unwrap();
        };
        for branch in network.branches() {
            for (i, &volume) in branch.mesh().volumes().iter().enumerate() {
                add(branch.cell(i).unwrap(), volume);
            }
        }
        for plenum in network.plenums() {
            add(plenum.state(), plenum.volume_m3());
        }
        (masses, energy)
    }

    #[test]
    fn accepted_chamber_compression_drives_connected_flow_and_closes_work_balance() {
        let gas = caloric_gas(true, 5.0);
        let mut network = connected_chamber(&gas);
        network
            .step(Time::new(1e-4), StepControl::default(), || false)
            .unwrap();
        let before = physical_storage(&network);
        let branches = branch_states(&network);
        let other = network.plenums()[1].clone();
        let clock = network.time();
        assert!(network.last_report().is_some());

        let compression = network
            .quasistatic_adiabatic_plenum_to(0, 0.002, || false)
            .unwrap();
        assert!(compression.volume.boundary_work_on_gas_j > 0.0);
        assert!(compression.volume.final_temperature_k > compression.volume.initial_temperature_k);
        assert!(compression.volume.final_pressure_pa > 2.0);
        assert_eq!(compression.plenum, 0);
        assert_eq!(compression.time, clock);
        assert_eq!(network.time(), clock);
        assert_eq!(branch_states(&network), branches);
        assert_eq!(network.plenums()[1], other);
        assert_eq!(network.plenums()[0].volume_m3(), 0.002);
        assert_eq!(network.last_report(), None);
        assert_eq!(network.accepted_port_exchanges(), None);
        assert_eq!(network.accepted_branch_reports(), None);
        assert_eq!(compression.after, network.totals());
        close(compression.conservation_residual.energy, 0.0, 1.0, 2e-12);

        let mut mass_from_chamber = 0.0;
        for _ in 0..30 {
            let report = network
                .step(Time::new(0.001), StepControl::default(), || false)
                .unwrap();
            assert_eq!(report.boundary_change.energy, 0.0);
            assert_eq!(&report.boundary_change.species_masses[..2], &[0.0; 2]);
            mass_from_chamber -= network.accepted_port_exchanges().unwrap()[0].mass_to_plenum();
        }
        assert!(
            mass_from_chamber > 0.0,
            "accepted compression did not reach the branch port"
        );
        let inlet = gas
            .primitive(network.branches()[0].cell(0).unwrap(), || false)
            .unwrap();
        assert!(inlet.pressure > 1.0);
        assert!(inlet.velocity > 0.0);
        assert!(network.time().value() > clock.value());
        let after = physical_storage(&network);
        for k in 0..2 {
            close(after.0[k], before.0[k], before.0[k], 2e-12);
        }
        close(
            after.1 - before.1,
            compression.volume.boundary_work_on_gas_j,
            1.0,
            3e-12,
        );
        assert_eq!(network.plenums()[1], other);
    }

    #[test]
    fn network_compression_expansion_keeps_closed_inventory_and_cancels_work() {
        let gas = caloric_gas(true, 0.0);
        let mut network = connected_chamber(&gas);
        let before = physical_storage(&network);
        let initial = gas
            .primitive(network.plenums()[0].state(), || false)
            .unwrap();
        let mut work = 0.0;
        for volume in [0.003, 0.002, 0.0025, 0.004] {
            let report = network
                .quasistatic_adiabatic_plenum_to(0, volume, || false)
                .unwrap();
            work += report.volume.boundary_work_on_gas_j;
            assert_eq!(report.time, Time::new(0.0));
            for residual in &report.conservation_residual.species_masses[..2] {
                close(*residual, 0.0, 1.0, 2e-16);
            }
        }
        let final_state = gas
            .primitive(network.plenums()[0].state(), || false)
            .unwrap();
        close(
            final_state.temperature,
            initial.temperature,
            initial.temperature,
            2e-11,
        );
        close(
            final_state.pressure,
            initial.pressure,
            initial.pressure,
            2e-11,
        );
        close(work, 0.0, 1.0, 2e-12);
        let after = physical_storage(&network);
        close(after.1 - before.1, work, 1.0, 2e-12);
        for k in 0..2 {
            close(after.0[k], before.0[k], before.0[k], 2e-15);
        }
    }

    #[test]
    fn network_volume_refusal_cancellation_and_noop_preserve_fields_clock_and_reports() {
        let gas = caloric_gas(true, 0.0);
        let mut network = connected_chamber(&gas);
        network
            .step(Time::new(1e-4), StepControl::default(), || false)
            .unwrap();
        let plenums = network.plenums().to_vec();
        let branches = branch_states(&network);
        let clock = network.time();
        let report = network.last_report();
        let ports = network.accepted_port_exchanges().unwrap().to_vec();
        let branch_reports = network.accepted_branch_reports().unwrap().to_vec();

        assert!(matches!(
            network.quasistatic_adiabatic_plenum_to(2, 0.002, || false),
            Err(Error::Invalid(_))
        ));
        assert!(matches!(
            network.quasistatic_adiabatic_plenum_to(0, 1e-15, || false),
            Err(Error::EnergyDomain)
        ));
        let mut calls = 0;
        assert_eq!(
            network
                .quasistatic_adiabatic_plenum_to(0, 0.002, || {
                    calls += 1;
                    calls >= 8
                })
                .unwrap_err(),
            Error::Cancelled,
        );
        let mut polls = 0;
        // A successful trial on an independent owner counts through the final
        // cancellation boundary, so this also exercises late rollback.
        let mut probe = connected_chamber(&gas);
        probe
            .step(Time::new(1e-4), StepControl::default(), || false)
            .unwrap();
        probe
            .quasistatic_adiabatic_plenum_to(0, 0.002, || {
                polls += 1;
                false
            })
            .unwrap();
        let mut seen = 0;
        assert_eq!(
            network
                .quasistatic_adiabatic_plenum_to(0, 0.002, || {
                    seen += 1;
                    seen >= polls
                })
                .unwrap_err(),
            Error::Cancelled
        );
        let noop = network
            .quasistatic_adiabatic_plenum_to(0, 0.004, || false)
            .unwrap();
        assert_eq!(noop.volume.boundary_work_on_gas_j, 0.0);
        assert_eq!(noop.before, noop.after);
        assert_eq!(network.plenums(), plenums);
        assert_eq!(branch_states(&network), branches);
        assert_eq!(network.time(), clock);
        assert_eq!(network.last_report(), report);
        assert_eq!(network.accepted_port_exchanges(), Some(ports.as_slice()));
        assert_eq!(
            network.accepted_branch_reports(),
            Some(branch_reports.as_slice())
        );
    }

    fn isolated_chamber(gas: &MixtureGas, pressure: f64) -> MixtureNetwork1d {
        let branch = MixtureEuler1d::new(
            gas.clone(),
            Length::new(10.0),
            0.01,
            vec![gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap(); 2],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        MixtureNetwork1d::new(
            vec![branch],
            vec![Plenum::new(1.0, gas.conserved(&[1.0, 0.0], 0.0, pressure).unwrap()).unwrap()],
            vec![],
            || false,
        )
        .unwrap()
    }

    #[test]
    fn g1_finite_time_chamber_compression_and_expansion_have_second_order_isentropic_limits() {
        for variable in [false, true] {
            let gas = caloric_gas(variable, 0.0);
            for rate in [-2.0_f64, 2.0] {
                let mut errors = Vec::new();
                for steps in [20, 40, 80] {
                    let mut network = isolated_chamber(&gas, 1.0);
                    let before = physical_storage(&network);
                    let dt = 0.25 / steps as f64;
                    let mut work = 0.0;
                    for _ in 0..steps {
                        let report = network
                            .step_with_volume_rates(
                                Time::new(dt),
                                StepControl::default(),
                                &[rate],
                                || false,
                            )
                            .unwrap();
                        assert_eq!(report.elapsed.value(), dt);
                        assert_eq!(report.halvings, 0);
                        assert_eq!(report.boundary_change.energy, 0.0);
                        work += report.boundary_work_on_gas_j;
                        close(report.conservation_residual.energy, 0.0, 1.0, 3e-14);
                    }
                    let chamber = &network.plenums()[0];
                    let volume = chamber.volume_m3();
                    close(volume, 1.0 + rate * 0.25, 1.0, 2e-14);
                    let actual = gas.primitive(chamber.state(), || false).unwrap();
                    // Independent source integral for pure A, R=1:
                    // Cv(T)=2.5+b*T; 2.5*ln(T)+b*(T-1)+ln(V)=0.
                    let b = if variable { 0.25 } else { 0.0 };
                    let (mut lo, mut hi) = (0.5_f64, 2.0_f64);
                    for _ in 0..80 {
                        let t = lo.midpoint(hi);
                        if 2.5 * t.ln() + b * (t - 1.0) + volume.ln() > 0.0 {
                            hi = t;
                        } else {
                            lo = t;
                        }
                    }
                    let exact_temperature = lo.midpoint(hi);
                    let exact_pressure = exact_temperature / volume;
                    errors.push((actual.pressure - exact_pressure).abs());
                    let after = physical_storage(&network);
                    close(after.0[0], before.0[0], before.0[0], 3e-14);
                    assert_eq!(after.0[1], 0.0);
                    close(after.1 - before.1, work, work.abs(), 8e-14);
                    assert_eq!(work.is_sign_positive(), rate.is_sign_negative());
                }
                for pair in errors.windows(2) {
                    assert!(
                        pair[1] > 1e-10,
                        "refinement error reached roundoff: {errors:?}"
                    );
                    let order = (pair[0] / pair[1]).log2();
                    assert!(
                        (1.8..2.2).contains(&order),
                        "observed order={order}, errors={errors:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn finite_time_connected_compression_closes_each_chambers_work_and_port_enthalpy() {
        let gas = caloric_gas(true, 5.0);
        let mut network = connected_chamber(&gas);
        let before = physical_storage(&network);
        let mut total_work = 0.0;
        let mut transported_mass = 0.0;
        for _ in 0..60 {
            let original = network.plenums().to_vec();
            let report = network
                .step_with_volume_rates(
                    Time::new(0.001),
                    StepControl::default(),
                    &[-0.02, 0.005],
                    || false,
                )
                .unwrap();
            let port = &network.accepted_port_exchanges().unwrap()[0];
            let mechanical = network.accepted_plenum_boundary_work().unwrap();
            for i in 0..2 {
                let old = &original[i];
                let next = &network.plenums()[i];
                let work = &mechanical[i];
                assert_eq!(work.initial_volume_m3, old.volume_m3());
                assert_eq!(work.final_volume_m3, next.volume_m3());
                let energy_change = next.volume_m3()
                    * gas.absolute_energy_density(next.state()).unwrap()
                    - old.volume_m3() * gas.absolute_energy_density(old.state()).unwrap();
                let port_energy = if i == 0 {
                    port.absolute_energy_to_plenum(&gas).unwrap()
                } else {
                    0.0
                };
                close(energy_change, port_energy + work.work_on_gas_j, 1.0, 2e-14);
                for k in 0..2 {
                    let mass_change = next.volume_m3() * next.state().partial_densities[k]
                        - old.volume_m3() * old.state().partial_densities[k];
                    let port_mass = if i == 0 {
                        port.species_masses_to_plenum[k]
                    } else {
                        0.0
                    };
                    close(mass_change, port_mass, 1.0, 4e-17);
                }
            }
            assert!(mechanical[0].work_on_gas_j > 0.0);
            assert!(mechanical[1].work_on_gas_j < 0.0);
            close(
                report.boundary_work_on_gas_j,
                mechanical[0].work_on_gas_j + mechanical[1].work_on_gas_j,
                1.0,
                2e-18,
            );
            close(report.conservation_residual.energy, 0.0, 1.0, 3e-14);
            total_work += report.boundary_work_on_gas_j;
            transported_mass += port.mass_to_plenum();
        }
        assert!(
            transported_mass < -1e-6,
            "compression did not drive finite-port outflow"
        );
        let after = physical_storage(&network);
        for k in 0..2 {
            close(after.0[k], before.0[k], 1.0, 2e-16);
        }
        close(after.1 - before.1, total_work, 1.0, 2e-13);
    }

    #[test]
    fn zero_volume_rates_preserve_rigid_step_fields_and_reports_exactly() {
        let gas = caloric_gas(true, 0.0);
        let mut rigid = connected_chamber(&gas);
        let mut moving = connected_chamber(&gas);
        for network in [&mut rigid, &mut moving] {
            network
                .quasistatic_adiabatic_plenum_to(0, 0.003, || false)
                .unwrap();
        }
        for _ in 0..4 {
            let expected = rigid
                .step(Time::new(0.001), StepControl::default(), || false)
                .unwrap();
            let actual = moving
                .step_with_volume_rates(
                    Time::new(0.001),
                    StepControl::default(),
                    &[0.0, -0.0],
                    || false,
                )
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(actual.boundary_work_on_gas_j, 0.0);
            assert_eq!(moving.plenums(), rigid.plenums());
            assert_eq!(branch_states(&moving), branch_states(&rigid));
            assert_eq!(moving.time(), rigid.time());
            assert_eq!(
                moving.accepted_port_exchanges(),
                rigid.accepted_port_exchanges()
            );
            assert_eq!(
                moving.accepted_branch_reports(),
                rigid.accepted_branch_reports()
            );
            assert_eq!(
                moving.accepted_plenum_boundary_work(),
                rigid.accepted_plenum_boundary_work()
            );
        }
    }

    #[test]
    fn finite_time_volume_rate_limits_admit_only_positive_geometry_and_source_states() {
        let gas = caloric_gas(false, 0.0);
        let mut network = isolated_chamber(&gas, 1.0);
        let report = network
            .step_with_volume_rates(Time::new(1.0), StepControl::default(), &[-100.0], || false)
            .unwrap();
        assert!(report.elapsed.value() < 0.01);
        assert!(report.halvings > 0);
        assert!(report.maximum_cfl <= StepControl::default().cfl);
        let work = network.accepted_plenum_boundary_work().unwrap()[0];
        assert_eq!(work.final_volume_m3, 1.0 - 100.0 * report.elapsed.value());
        assert!(work.final_volume_m3 > 0.0 && work.final_volume_m3 < 1.0);
        assert!(work.work_on_gas_j > 0.0);

        let mut hot = isolated_chamber(&gas, 49.9);
        let original = hot.plenums().to_vec();
        let controls = StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        };
        assert!(matches!(
            hot.step_with_volume_rates(Time::new(0.01), controls, &[-1.0], || false),
            Err(Error::StepRetries { .. })
        ));
        assert_eq!(hot.plenums(), original);
        assert_eq!(hot.time().value(), 0.0);
        assert!(hot.last_report().is_none());
        assert!(hot.accepted_plenum_boundary_work().is_none());
    }

    #[test]
    fn finite_time_volume_refusal_and_early_late_cancellation_preserve_accepted_transaction() {
        let gas = caloric_gas(true, 0.0);
        let prepare = || {
            let mut network = connected_chamber(&gas);
            network
                .step_with_volume_rates(
                    Time::new(0.001),
                    StepControl::default(),
                    &[-0.02, 0.005],
                    || false,
                )
                .unwrap();
            network
        };
        let snapshot = |network: &MixtureNetwork1d| {
            (
                branch_states(network),
                network.plenums().to_vec(),
                network.time(),
                network.last_report(),
                network.accepted_port_exchanges().unwrap().to_vec(),
                network.accepted_branch_reports().unwrap().to_vec(),
                network.accepted_plenum_boundary_work().unwrap().to_vec(),
            )
        };
        let mut network = prepare();
        let before = snapshot(&network);
        for rates in [
            &[][..],
            &[0.0][..],
            &[f64::NAN, 0.0][..],
            &[f64::MAX, 0.0][..],
            &[f64::MIN_POSITIVE, 0.0][..],
        ] {
            assert!(
                network
                    .step_with_volume_rates(Time::new(0.001), StepControl::default(), rates, || {
                        false
                    })
                    .is_err()
            );
            assert_eq!(snapshot(&network), before);
        }
        let mut polls = 0;
        prepare()
            .step_with_volume_rates(
                Time::new(0.001),
                StepControl::default(),
                &[-0.02, 0.005],
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        for stop in [1, polls / 2, polls] {
            let mut network = prepare();
            let before = snapshot(&network);
            let mut seen = 0;
            assert_eq!(
                network
                    .step_with_volume_rates(
                        Time::new(0.001),
                        StepControl::default(),
                        &[-0.02, 0.005],
                        || {
                            seen += 1;
                            seen >= stop
                        }
                    )
                    .unwrap_err(),
                Error::Cancelled
            );
            assert_eq!(snapshot(&network), before);
        }
    }
}
