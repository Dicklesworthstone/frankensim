//! Compensated original-cell and original-beam physical accounts.

use super::*;

pub(super) struct BeamTotals<const B: usize> {
    transmitted: [Sum; B],
    absorbed: [Sum; B],
    passive: [Sum; B],
}

impl<const B: usize> BeamTotals<B> {
    pub(super) fn new() -> Self {
        Self {
            transmitted: [Sum::default(); B],
            absorbed: [Sum::default(); B],
            passive: [Sum::default(); B],
        }
    }

    pub(super) fn add_packet(
        &mut self,
        lane: usize,
        packet: cell::ChannelStep,
    ) -> Result<(), SingletTripletError> {
        self.absorbed[lane].add(packet.absorbed)?;
        self.passive[lane].add(packet.passive)?;
        Ok(())
    }

    pub(super) fn add_transmitted(
        &mut self,
        rates: [f64; B],
        duration: f64,
    ) -> Result<(), SingletTripletError> {
        for (sum, rate) in self.transmitted.iter_mut().zip(rates) {
            sum.add(positive_product(rate, duration, "transmitted photons")?)?;
        }
        Ok(())
    }

    pub(super) fn finish(
        self,
        rates: [f64; B],
        duration: f64,
        quantum: [f64; B],
        tolerance: f64,
    ) -> Result<PumpBeamStep<B>, SingletTripletError> {
        let mut incident_photons = [0.0; B];
        for lane in 0..B {
            incident_photons[lane] = positive_product(rates[lane], duration, "incident photons")?;
        }
        let transmitted_photons = values(self.transmitted)?;
        let absorbed_photons = values(self.absorbed)?;
        let passive_absorbed_photons = values(self.passive)?;
        let photon_residual = photon_balance(
            incident_photons, transmitted_photons, absorbed_photons, passive_absorbed_photons,
        )?;
        admit_photons(
            incident_photons, transmitted_photons, absorbed_photons,
            passive_absorbed_photons, photon_residual, tolerance,
        )?;
        Ok(PumpBeamStep {
            incident_photons, transmitted_photons, absorbed_photons, passive_absorbed_photons,
            incident_energy_j: energies(incident_photons, quantum)?,
            transmitted_energy_j: energies(transmitted_photons, quantum)?,
            absorbed_energy_j: energies(absorbed_photons, quantum)?,
            passive_absorbed_energy_j: energies(passive_absorbed_photons, quantum)?,
            photon_residual,
        })
    }
}

impl<const B: usize> SaturablePumpBeams<B> {
    pub(super) fn energy_balance(
        &self,
        initial: &[[f64; 3]],
        cells: &[SaturableSingletTripletStep],
        passive: &[f64],
        beams: &[PumpBeamStep<B>],
        tolerance: f64,
    ) -> Result<f64, SingletTripletError> {
        let mut residual = Sum::default();
        let mut scale = Sum::default();
        let mut spectral_minus_target = Sum::default();
        let mut pump_scale = Sum::default();
        for beam in beams {
            for lane in 0..B {
                residual.add(beam.incident_energy_j[lane])?;
                residual.add(-beam.transmitted_energy_j[lane])?;
                scale.add(beam.incident_energy_j[lane])?;
                scale.add(beam.transmitted_energy_j[lane])?;
                spectral_minus_target.add(beam.absorbed_energy_j[lane])?;
                pump_scale.add(beam.absorbed_energy_j[lane])?;
            }
        }
        let mut heat = Sum::default();
        for (index, cell) in cells.iter().enumerate() {
            for target in 0..2 {
                let old = positive_product(
                    initial[index][target + 1],
                    self.models[index].excitation_energies_j()[target],
                    "initial excitation",
                )?;
                residual.add(old)?;
                scale.add(old)?;
                for value in [cell.stored_energy_j[target], cell.radiative_decay_energy_j[target]] {
                    residual.add(-value)?;
                    scale.add(value)?;
                }
                spectral_minus_target.add(-cell.absorbed_pump_energy_j[target])?;
                pump_scale.add(cell.absorbed_pump_energy_j[target])?;
                heat.add(cell.pump_relaxation_heat_j[target])?;
                heat.add(cell.nonradiative_heat_j[target])?;
            }
            heat.add(cell.crossing_heat_j)?;
            heat.add(-cell.reverse_crossing_heat_absorbed_j)?;
            heat.add(passive[index])?;
        }
        admit_energy(spectral_minus_target.value()?, pump_scale.value()?, tolerance)?;
        let heat = heat.value()?;
        residual.add(-heat)?;
        scale.add(heat.abs())?;
        let residual = residual.value()?;
        admit_energy(residual, scale.value()?, tolerance)?;
        Ok(residual)
    }

    pub(super) fn components(
        &self,
        step: &PumpBeamsStep<B>,
        duration: f64,
    ) -> Result<Vec<f64>, SingletTripletError> {
        let mut result = Vec::new();
        for (index, (model, k)) in self.models.iter().zip(&step.cells).enumerate() {
            result.extend(k.stored_energy_j);
            result.extend(k.absorbed_pump_energy_j);
            result.extend(k.radiative_decay_energy_j);
            result.extend(k.nonradiative_heat_j);
            result.extend(k.pump_relaxation_heat_j);
            result.push(k.crossing_heat_j);
            result.push(k.reverse_crossing_heat_absorbed_j);
            result.push(step.passive_heat_j[index]);
            let quantum = model.excitation_energies_j();
            result.push(positive_product(
                k.populations[0], quantum[0], "ground convergence weight",
            )?);
            for target in 0..2 {
                result.push(product3(
                    k.stored_energy_j[target], model.radiative_rates_s()[target], duration,
                )?);
                for events in [
                    k.pump_events[target], k.radiative_events[target],
                    k.nonradiative_events[target], k.crossing_events[target],
                ] {
                    result.push(positive_product(
                        events, quantum[target], "event convergence weight",
                    )?);
                }
            }
        }
        for beam in &step.beams {
            result.extend(beam.transmitted_energy_j);
            result.extend(beam.absorbed_energy_j);
            result.extend(beam.passive_absorbed_energy_j);
        }
        Ok(result)
    }
}
