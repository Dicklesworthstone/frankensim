//! Retained material inventories for the same conservative reservoir exchange.
//! NASA9 keeps its original parcel and incremental solve. The named near-ambient
//! model owns a constant Cv, so its inverse is analytic; no transport/EOS fit or
//! thermochemical source pack is synthesized from a numerical gas snapshot.

use super::{NsfDgError, NsfFiniteReservoirActuationFrame, NsfFiniteReservoirLimits,
    NsfFiniteReservoirSource, material, poll};
use fs_exec::CancelGate;
use fs_lbm::nsf3::{NamedMixtureNsfGas, NsfGasMedium, NsfOpenState};
use fs_material::gas::GasState;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_qty::Temperature;

mod actuation;

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
        // Flow reconstruction retains its original rho*V arithmetic. Closed
        // heat/volume motion then preserves this mass exactly, as NASA does.
        mass_kg: f64,
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
                    mass_kg: gas.density * volume_m3,
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
            Self::NamedMixture { mass_kg, .. } => *mass_kg,
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
            Self::NamedMixture { source, gas, volume_m3, temperature_window_k, .. } => {
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
                    mass_kg: next.density * volume_m3,
                    temperature_window_k: *temperature_window_k,
                }, 0))
            }
        }
    }

    /// Preserve the original parcel's simultaneous heat/volume update and
    /// receipt. Named mixtures solve the same conservative midpoint with their
    /// original constant Cv. Return its actual material midpoint for NSF
    /// admission, without changing the published receipt shape.
    pub(super) fn midpoint_heat_volume(
        &self, volume_m3: f64, heat_j: f64,
        limits: &NsfFiniteReservoirLimits, gate: &CancelGate,
    ) -> Result<(Self, NsfFiniteReservoirActuationFrame, GasState), NsfDgError> {
        actuation::candidate(self, volume_m3, heat_j, limits, gate)
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
#[path = "../../../../fs-lbm/tests/nsf3_thermochemical/source.rs"]
mod source;

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


    fn actuation_limits() -> NsfFiniteReservoirLimits {
        NsfFiniteReservoirLimits {
            maximum_relative_mass_exchange: 0.1,
            maximum_relative_energy_exchange: 0.1,
            mass_absolute_tolerance_kg: 1e-16,
            energy_absolute_tolerance_j: 1e-10,
            relative_tolerance: 1e-8,
        }
    }

    #[test]
    fn named_heat_and_motion_preserve_mass_and_original_source_energy() {
        let old = named(0.005, 100000.0);
        let gate = CancelGate::new_clock_free();
        for (volume_ratio, heat) in [(0.98, 0.001), (1.02, -0.001), (1.0, 0.001), (1.0, -0.001)] {
            let volume = old.volume_m3() * volume_ratio;
            let (next, frame, midpoint) = old.midpoint_heat_volume(volume, heat, &actuation_limits(), &gate).unwrap();
            let step = frame.step;
            assert_eq!(next.mass_kg().to_bits(), old.mass_kg().to_bits());
            assert_eq!(next.volume_m3(), volume);
            assert_eq!(step.supplied_heat_j, heat);
            assert_eq!(step.iterations, 0);
            assert!(step.energy_residual_j.abs() < 1e-10);
            assert!((next.state().density * volume - old.mass_kg()).abs() < 1e-16);
            assert!((midpoint.density * (0.5 * volume + 0.5 * old.volume_m3()) - old.mass_kg()).abs() < 1e-16);
            assert!((step.stored_energy_change_j + step.work_by_gas_j - heat).abs() < 1e-10);
            assert_eq!(step.work_by_gas_j, frame.midpoint_pressure_pa * (volume - old.volume_m3()));
            assert_eq!(midpoint.pressure, frame.midpoint_pressure_pa);
            assert!(step.entropy_change_j_k.is_finite());
            if volume_ratio < 1.0 {
                assert!(step.work_by_gas_j < 0.0);
                assert!(next.state().pressure > old.state().pressure);
            } else if volume_ratio > 1.0 {
                assert!(step.work_by_gas_j > 0.0);
                assert!(next.state().pressure < old.state().pressure);
            } else {
                assert_eq!(step.work_by_gas_j, 0.0);
                let cv = old.state().specific_heat_cp - old.state().specific_gas_constant;
                let expected = old.state().temperature + heat / (old.mass_kg() * cv);
                assert!((next.state().temperature - expected).abs() < 1e-12);
                assert_eq!(step.entropy_change_j_k.signum(), heat.signum());
            }
            let Inventory::NamedMixture { source, .. } = &next else { panic!("lost named source") };
            assert_eq!(*next.state(), source.at(next.state().temperature, next.state().pressure).unwrap());
            assert_eq!(midpoint, source.at(midpoint.temperature, midpoint.pressure).unwrap());
            let Inventory::NamedMixture { source: before, .. } = &old else { unreachable!() };
            assert_eq!(source.species(), before.species());
            assert_eq!(source.composition(), before.composition());
            assert_eq!(source.reference_state(), before.reference_state());
        }
    }

    #[test]
    fn named_receipt_derivatives_follow_the_actual_midpoint_response() {
        let old = named(0.005, 100000.0);
        let gate = CancelGate::new_clock_free();
        let volume = old.volume_m3() * 0.98;
        let heat = 0.001;
        let response = |v, q| old.midpoint_heat_volume(v, q, &actuation_limits(), &gate).unwrap().1;
        let receipt = response(volume, heat);
        let dv = old.volume_m3() * 1e-5;
        let dq = 1e-6;
        let volume_derivative = (response(volume + dv, heat).midpoint_pressure_pa
            - response(volume - dv, heat).midpoint_pressure_pa) / (2.0 * dv);
        let heat_derivative = (response(volume, heat + dq).midpoint_pressure_pa
            - response(volume, heat - dq).midpoint_pressure_pa) / (2.0 * dq);
        for (actual, expected) in [
            (volume_derivative, receipt.step.pressure_final_volume_derivative_pa_m3),
            (heat_derivative, receipt.step.pressure_heat_derivative_pa_j),
        ] {
            assert!((actual - expected).abs() <= 2e-8 * expected.abs());
        }
    }

    #[test]
    fn unchanged_actuation_keeps_passive_storage_and_source_bits() {
        let old = named(0.005, 100000.0);
        let (next, frame, midpoint) = old.midpoint_heat_volume(old.volume_m3(), 0.0, &actuation_limits(),
            &CancelGate::new_clock_free()).unwrap();
        assert_eq!(next.state(), old.state());
        assert_eq!(next.mass_kg().to_bits(), old.mass_kg().to_bits());
        assert_eq!(next.volume_m3().to_bits(), old.volume_m3().to_bits());
        assert_eq!(next.sensible_specific_energy_j_kg().to_bits(), old.sensible_specific_energy_j_kg().to_bits());
        assert_eq!(frame.midpoint_pressure_pa, old.state().pressure);
        assert_eq!(midpoint, *old.state());
        assert_eq!(frame.step.stored_energy_change_j, 0.0);
        assert_eq!(frame.step.work_by_gas_j, 0.0);
        assert_eq!(frame.step.energy_residual_j, 0.0);
        assert_eq!(frame.step.entropy_change_j_k, 0.0);
    }

    #[test]
    fn named_motion_refuses_original_phase_pressure_window_and_cancellation() {
        let gate = CancelGate::new_clock_free();
        for (water, pressure, ratio) in [(0.018, 190000.0, 1.02), (0.0, 199000.0, 0.98)] {
            let old = named(water, pressure);
            let before = (*old.state(), old.mass_kg(), old.volume_m3());
            assert!(matches!(old.midpoint_heat_volume(old.volume_m3() * ratio, 0.0, &actuation_limits(), &gate),
                Err(NsfDgError::Flow(_))));
            assert_eq!((*old.state(), old.mass_kg(), old.volume_m3()), before);
        }
        let old = named(0.0, 100000.0);
        assert!(old.midpoint_heat_volume(old.volume_m3(), 1.0, &actuation_limits(), &gate).is_err());
        gate.request();
        assert!(matches!(old.midpoint_heat_volume(old.volume_m3() * 0.99, 0.0, &actuation_limits(), &gate),
            Err(NsfDgError::Cancelled)));
    }

    #[test]
    fn nasa_actuation_reuses_original_midpoint_and_complete_receipt() {
        use fs_lbm::nsf3::ThermochemicalNsfGas;
        let gate = CancelGate::new_clock_free();
        let controls = GasEnergySolve { minimum_temperature_k: 290.0, maximum_temperature_k: 320.0,
            energy_tolerance_j_kg: 1e-6, maximum_iterations: 64 };
        for offset in [0.0, -8e8] {
            let model = source::model(offset, 250000.0, 1.0);
            let point = source::point(300.0, 100000.0);
            let reference = *model.evaluate(&point).unwrap().state();
            let source = ThermochemicalNsfGas::new(model, point, controls).unwrap();
            let old = Inventory::new(&source.into(), NsfOpenState {
                density_kg_m3: reference.density, velocity_m_s: [0.0; 3], temperature_k: reference.temperature,
            }, 1e-5, [290.0, 320.0], &gate).unwrap();
            for (ratio, heat) in [(0.98, 0.001), (1.0, 0.0)] {
                let volume = old.volume_m3() * ratio;
                let expected = old.thermochemical().unwrap()
                    .midpoint_heat_volume(volume, heat, controls, || false).unwrap();
                let (next, frame, midpoint) = old.midpoint_heat_volume(volume, heat, &actuation_limits(), &gate).unwrap();
                assert_eq!(next.thermochemical().unwrap(), expected.endpoint());
                assert_eq!(next.energy_reference_j_kg().to_bits(), old.energy_reference_j_kg().to_bits());
                assert_eq!(next.mass_kg().to_bits(), old.mass_kg().to_bits());
                assert_eq!(frame.midpoint_pressure_pa, expected.midpoint().state().pressure);
                assert_eq!(midpoint, *expected.midpoint().state());
                assert_eq!(frame.step, *expected.step());
                assert!(frame.step.energy_residual_j.abs() < 1e-10);
            }
        }
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
