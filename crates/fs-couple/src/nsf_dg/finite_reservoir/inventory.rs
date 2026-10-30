//! Retained material inventories for the same conservative reservoir exchange.
//! NASA9 keeps its original parcel and incremental solve. The named near-ambient
//! model owns a constant Cv, so its inverse is analytic; no transport/EOS fit or
//! thermochemical source pack is synthesized from a numerical gas snapshot.

use super::{NsfDgError, NsfFiniteReservoirActuationFrame, NsfFiniteReservoirSource, material, poll};
use fs_exec::CancelGate;
use fs_lbm::nsf3::{NamedMixtureNsfGas, NsfGasMedium, NsfOpenState};
use fs_material::gas::GasState;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_qty::Temperature;

#[derive(Clone)]
pub(super) enum Inventory {
    Thermochemical {
        gas: ClosedGasParcel,
        energy_solve: GasEnergySolve,
        energy_reference_j_kg: f64,
    },
    NamedMixture {
        source: NamedMixtureNsfGas,
        gas: GasState,
        volume_m3: f64,
        temperature_window_k: [f64; 2],
    },
}

impl Inventory {
    pub(super) fn new(
        medium: &NsfGasMedium,
        opening: NsfOpenState,
        volume_m3: f64,
        temperature_window_k: [f64; 2],
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        match medium {
            NsfGasMedium::Thermochemical(source) => {
                let point = source.context().clone().with("temperature", opening.temperature_k)
                    .map_err(|_| NsfDgError::Input("finite reservoir temperature coordinate is invalid"))?;
                // Preserve the original zero-increment pressure reconstruction.
                let solved = source.model().reconstruct_density_internal_energy_increment(
                    &point, opening.density_kg_m3, 0.0, source.energy_solve(),
                    || gate.is_requested(),
                ).map_err(material)?;
                let point = point.with("temperature", solved.state().state().temperature)
                    .and_then(|p| p.with("pressure", solved.state().state().pressure))
                    .map_err(|_| NsfDgError::Input("finite reservoir source coordinates are invalid"))?;
                let gas = ClosedGasParcel::new(source.model().clone(), point, volume_m3)
                    .map_err(material)?;
                Ok(Self::Thermochemical {
                    gas, energy_solve: source.energy_solve(),
                    energy_reference_j_kg: source.energy_reference_j_kg(),
                })
            }
            NsfGasMedium::NamedMixture(source) => {
                let pressure = opening.density_kg_m3
                    * source.reference_state().specific_gas_constant * opening.temperature_k;
                let gas = source.at(opening.temperature_k, pressure).map_err(NsfDgError::Flow)?;
                if !(temperature_window_k[0]..=temperature_window_k[1]).contains(&gas.temperature) {
                    return Err(NsfDgError::Input("finite reservoir temperature lies outside the NSF run window"));
                }
                Ok(Self::NamedMixture {
                    source: source.clone(), gas, volume_m3, temperature_window_k,
                })
            }
            NsfGasMedium::CaloricallyPerfect(_) => Err(NsfDgError::Input(
                "finite reservoir requires the original named-mixture or thermochemical NSF medium",
            )),
        }
    }

    pub(super) fn source_kind(&self) -> NsfFiniteReservoirSource {
        match self {
            Self::Thermochemical { .. } => NsfFiniteReservoirSource::Thermochemical,
            Self::NamedMixture { .. } => NsfFiniteReservoirSource::NamedMixture,
        }
    }

    pub(super) fn thermochemical(&self) -> Option<&ClosedGasParcel> {
        match self { Self::Thermochemical { gas, .. } => Some(gas), _ => None }
    }

    pub(super) fn state(&self) -> &GasState {
        match self {
            Self::Thermochemical { gas, .. } => gas.state().state(),
            Self::NamedMixture { gas, .. } => gas,
        }
    }

    pub(super) fn volume_m3(&self) -> f64 {
        match self {
            Self::Thermochemical { gas, .. } => gas.volume_m3(),
            Self::NamedMixture { volume_m3, .. } => *volume_m3,
        }
    }

    pub(super) fn mass_kg(&self) -> f64 {
        match self {
            Self::Thermochemical { gas, .. } => gas.mass_kg(),
            Self::NamedMixture { gas, volume_m3, .. } => gas.density * volume_m3,
        }
    }

    pub(super) fn energy_reference_j_kg(&self) -> f64 {
        match self {
            Self::Thermochemical { energy_reference_j_kg, .. } => *energy_reference_j_kg,
            Self::NamedMixture { .. } => 0.0,
        }
    }

    pub(super) fn sensible_specific_energy_j_kg(&self) -> f64 {
        match self {
            Self::Thermochemical { gas, energy_reference_j_kg, .. } =>
                gas.state().specific_internal_energy_j_kg() - energy_reference_j_kg,
            Self::NamedMixture { gas, .. } =>
                (gas.specific_heat_cp - gas.specific_gas_constant) * gas.temperature,
        }
    }

    pub(super) fn temperature_window_k(&self) -> [f64; 2] {
        match self {
            Self::Thermochemical { energy_solve, .. } =>
                [energy_solve.minimum_temperature_k, energy_solve.maximum_temperature_k],
            Self::NamedMixture { temperature_window_k, .. } => *temperature_window_k,
        }
    }

    pub(super) fn matches_source(&self, medium: &NsfGasMedium, window: [f64; 2]) -> bool {
        match (self, medium) {
            (Self::Thermochemical { gas, energy_solve, energy_reference_j_kg },
             NsfGasMedium::Thermochemical(source)) => gas.model() == source.model()
                && energy_reference_j_kg.to_bits() == source.energy_reference_j_kg().to_bits()
                && *energy_solve == source.energy_solve(),
            (Self::NamedMixture { source: retained, temperature_window_k, .. },
             NsfGasMedium::NamedMixture(source)) => retained.species() == source.species()
                && retained.composition() == source.composition()
                && retained.reference_state() == source.reference_state()
                && *temperature_window_k == window,
            _ => false,
        }
    }

    pub(super) fn specific_energy_change(
        &self, from_k: f64, to_k: f64, refusal: &'static str,
    ) -> Result<f64, NsfDgError> {
        match self {
            Self::Thermochemical { gas, .. } => gas.model().thermodynamics()
                .specific_internal_energy_change(Temperature::new(from_k), Temperature::new(to_k))
                .map(|value| value.value()).map_err(|_| NsfDgError::Input(refusal)),
            Self::NamedMixture { source, .. } => {
                let state = source.reference_state();
                let change = (state.specific_heat_cp - state.specific_gas_constant) * (to_k - from_k);
                if !change.is_finite() { return Err(NsfDgError::Input(refusal)); }
                Ok(change)
            }
        }
    }

    pub(super) fn reconstruct_energy_change(
        &self, target_mass: f64, change: f64, gate: &CancelGate,
    ) -> Result<(Self, usize), NsfDgError> {
        match self {
            Self::Thermochemical { gas, energy_solve, energy_reference_j_kg } => {
                let solved = gas.model().reconstruct_density_internal_energy_increment(
                    gas.point(), target_mass / gas.volume_m3(), change, *energy_solve,
                    || gate.is_requested(),
                ).map_err(material)?;
                let point = gas.point().clone().with("temperature", solved.state().state().temperature)
                    .and_then(|p| p.with("pressure", solved.state().state().pressure))
                    .map_err(|_| NsfDgError::Input("finite reservoir endpoint coordinates are invalid"))?;
                let next = ClosedGasParcel::new(gas.model().clone(), point, gas.volume_m3())
                    .map_err(material)?;
                Ok((Self::Thermochemical {
                    gas: next, energy_solve: *energy_solve,
                    energy_reference_j_kg: *energy_reference_j_kg,
                }, solved.iterations()))
            }
            Self::NamedMixture { source, gas, volume_m3, temperature_window_k } => {
                poll(gate)?;
                let reference = source.reference_state();
                let cv = reference.specific_heat_cp - reference.specific_gas_constant;
                // Stable source-owned energy increment, retaining the original
                // total-energy debit and the existing subsequent storage audit.
                let temperature = gas.temperature + change / cv;
                let pressure = target_mass / volume_m3 * reference.specific_gas_constant * temperature;
                let next = source.at(temperature, pressure).map_err(NsfDgError::Flow)?;
                if !(temperature_window_k[0]..=temperature_window_k[1]).contains(&next.temperature) {
                    return Err(NsfDgError::Input("finite reservoir endpoint leaves the NSF run window"));
                }
                poll(gate)?;
                Ok((Self::NamedMixture {
                    source: source.clone(), gas: next, volume_m3: *volume_m3,
                    temperature_window_k: *temperature_window_k,
                }, 0))
            }
        }
    }

    /// Preserve the original parcel's simultaneous heat/volume update, source
    /// laws, fixed mass and midpoint work. Named mixtures remain usable for
    /// passive finite storage; they do not have this original mechanical owner.
    pub(super) fn midpoint_heat_volume(
        &self, volume_m3: f64, heat_j: f64, gate: &CancelGate,
    ) -> Result<(Self, NsfFiniteReservoirActuationFrame), NsfDgError> {
        let Self::Thermochemical { gas, energy_solve, energy_reference_j_kg } = self else {
            return Err(NsfDgError::Input(
                "finite reservoir actuation requires its original thermochemical parcel",
            ));
        };
        let candidate = gas.midpoint_heat_volume(
            volume_m3, heat_j, *energy_solve, || gate.is_requested(),
        ).map_err(material)?;
        let midpoint_pressure_pa = candidate.midpoint().state().pressure;
        let (gas, step) = candidate.into_parts();
        Ok((Self::Thermochemical {
            gas, energy_solve: *energy_solve,
            energy_reference_j_kg: *energy_reference_j_kg,
        }, NsfFiniteReservoirActuationFrame { midpoint_pressure_pa, step }))
    }

    /// Named Cv and R are constant and positive: T = E/(m Cv), p = E R/(V Cv).
    /// Outward rounding bounds all chronological prefixes of the original
    /// exchange variations. The original Buck saturation pressure increases
    /// throughout this named model's admitted 273.15..323.15 K window, so
    /// Tmin/Pmax is a conservative phase bound even when those extrema cannot
    /// occur together. No successful endpoint can hide a condensing prefix.
    pub(super) fn admit_domain_envelope(
        &self, low_mass: f64, high_mass: f64, low_energy: f64, high_energy: f64,
        gate: &CancelGate,
    ) -> Result<(), NsfDgError> {
        let Self::NamedMixture { source, volume_m3, temperature_window_k, .. } = self else {
            return Ok(());
        };
        poll(gate)?;
        let reference = source.reference_state();
        let cv = reference.specific_heat_cp - reference.specific_gas_constant;
        let r = reference.specific_gas_constant;
        let t_low = ((low_energy / high_mass).next_down() / cv).next_down();
        let t_high = ((high_energy / low_mass).next_up() / cv).next_up();
        let p_low = (((low_energy / volume_m3).next_down() / cv).next_down() * r).next_down();
        let p_high = (((high_energy / volume_m3).next_up() / cv).next_up() * r).next_up();
        if ![t_low, t_high, p_low, p_high].into_iter().all(|v| v.is_finite() && v > 0.0)
            || t_low < temperature_window_k[0] || t_high > temperature_window_k[1]
        {
            return Err(NsfDgError::Input(
                "finite reservoir named-mixture prefix envelope leaves the NSF run window; refine the interval",
            ));
        }
        source.at(t_low, p_high).map_err(NsfDgError::Flow)?;
        source.at(t_high, p_low).map_err(NsfDgError::Flow)?;
        poll(gate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas_species::GasSpecies;
    use fs_qty::semantic::{Composition, CompositionBasis};

    fn named(water: f64, pressure: f64) -> Inventory {
        let (species, fractions) = if water == 0.0 {
            (vec![GasSpecies::Nitrogen], vec![1.0])
        } else {
            (vec![GasSpecies::Nitrogen, GasSpecies::WaterVapor], vec![1.0 - water, water])
        };
        let source = NamedMixtureNsfGas::new(species,
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(), 300.0, 100000.0).unwrap();
        // The whole NSF reference-pressure window admits this mixture. A
        // separate pressurized opening can still approach its phase boundary.
        source.at(290.0, 100000.0).unwrap();
        let density = source.at(300.0, pressure).unwrap().density;
        Inventory::new(&source.into(), NsfOpenState {
            density_kg_m3: density, velocity_m_s: [0.0; 3], temperature_k: 300.0,
        }, 1e-5, [290.0, 320.0], &CancelGate::new_clock_free()).unwrap()
    }

    #[test]
    fn g0_named_prefix_domain_checks_pressure_and_phase_before_publication() {
        let gate = CancelGate::new_clock_free();
        for (water, pressure) in [(0.018, 190000.0), (0.0, 199000.0)] {
            let gas = named(water, pressure);
            let mass = gas.mass_kg();
            let energy = mass * gas.sensible_specific_energy_j_kg();
            // The unchanged endpoint is admitted, and a small surrounding
            // interval remains usable. Larger opposing transfers can cross
            // phase/pressure support even when their net changes are zero.
            gas.admit_domain_envelope(mass * 0.9999, mass * 1.0001,
                energy * 0.9999, energy * 1.0001, &gate).unwrap();
            assert!(matches!(gas.admit_domain_envelope(mass * 0.99, mass * 1.01,
                energy * 0.99, energy * 1.01, &gate), Err(NsfDgError::Flow(_))));
        }
        let gas = named(0.0, 100000.0);
        let mass = gas.mass_kg();
        let energy = mass * gas.sensible_specific_energy_j_kg();
        assert!(matches!(gas.admit_domain_envelope(mass * 0.95, mass * 1.05,
            energy * 0.95, energy * 1.05, &gate), Err(NsfDgError::Input(_))));
    }
}
