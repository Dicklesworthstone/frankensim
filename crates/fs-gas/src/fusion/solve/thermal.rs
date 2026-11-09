//! Original-source thermal operators within the finite-fuel midpoint solve.

use super::super::thermal::{
    admit_confinement_resolution, confined_energy, evaluate, exchange_energy, radiated_energy,
};
use super::*;

pub(super) use super::super::thermal::ThermalSources as Sources;

impl Problem<'_> {
    pub(super) fn sources(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sources, DtFusionError> {
        let temperatures = self.admit(point, cancelled)?;
        let amounts = self.amounts(point);
        let volume = self.initial.volume();
        evaluate(self.transport, &amounts, volume, temperatures, cancelled)
    }

    /// Freeze coefficients at the predictor midpoint, then compose thermal
    /// half steps around the exact finite-fuel/external-power map. The thermal
    /// map itself symmetrically composes confinement, analytic free-free cooling
    /// and exact reciprocal exchange. Every stage uses actual heat capacities.
    pub(super) fn frozen_coupled_step(
        &self,
        old: Point,
        dt: f64,
        sources: Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        if self.transport.energy_confinement.is_some() {
            let reacting = old.limiting_fuel > 0.0;
            admit_confinement_resolution(
                self.transport,
                sources,
                dt,
                [
                    reacting || self.heating.heavy_power.value() != 0.0,
                    (reacting && self.heating.alpha_electron_fraction > 0.0)
                        || self.heating.electron_power.value() != 0.0,
                ],
            )?;
        }
        let supplied = self.exchange_particles(old, 0.5 * dt, cancelled)?;
        let first = self.thermal_step(supplied, 0.5 * dt, sources, cancelled)?;
        let burned = self.frozen_step(first, dt, sources.temperatures[0])?;
        self.admit(burned, cancelled)?;
        let next = self.thermal_step(burned, 0.5 * dt, sources, cancelled)?;
        self.admit(next, cancelled)?;
        self.exchange_particles(next, 0.5 * dt, cancelled)
    }

    fn thermal_step(
        &self,
        old: Point,
        dt: f64,
        sources: Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        poll(cancelled)?;
        if sources.conductance_w_k == 0.0
            && sources.radiation_w_per_sqrt_k == 0.0
            && self.transport.energy_confinement.is_none()
        {
            return Ok(old);
        }
        let confined = self.confine(old, 0.5 * dt, cancelled)?;
        let first = self.radiate(confined, 0.5 * dt, sources.radiation_w_per_sqrt_k)?;
        self.admit(first, cancelled)?;
        let exchanged = self.exchange(first, dt, sources.conductance_w_k)?;
        self.admit(exchanged, cancelled)?;
        let radiated = self.radiate(exchanged, 0.5 * dt, sources.radiation_w_per_sqrt_k)?;
        self.admit(radiated, cancelled)?;
        self.confine(radiated, 0.5 * dt, cancelled)
    }

    /// Fixed-composition Cv values from the admitted original monatomic bank.
    fn capacities(&self, point: Point) -> [f64; 2] {
        let amounts = self.amounts(point);
        [
            CV_MOLAR
                * self
                    .plasma
                    .indices()
                    .iter()
                    .map(|&i| amounts[i].value())
                    .sum::<f64>(),
            CV_MOLAR * amounts[self.plasma.electron_index()].value(),
        ]
    }

    /// Exact constant-G relaxation: the temperature difference decays as
    /// exp[-G(1/Ch+1/Ce)dt]. Reciprocal energy transfer cannot overshoot the
    /// shared equilibrium, including steps much longer than the exchange time.
    fn exchange(&self, old: Point, dt: f64, conductance: f64) -> Result<Point, DtFusionError> {
        if conductance == 0.0 || dt == 0.0 {
            return Ok(old);
        }
        let transfer = exchange_energy(
            self.temperatures(old)?,
            self.capacities(old),
            dt,
            conductance,
        )?;
        let mut next = old;
        next.coulomb_transfer.add(transfer);
        self.add_pool_energy(&mut next, [transfer, -transfer]);
        Ok(next)
    }

    /// With frozen charge populations and Gaunt factor, P=C*sqrt(Te), hence
    /// sqrt(Te)'=-C/(2Ce). Compute the outgoing energy without subtracting
    /// nearly identical squared temperatures. Crossing zero is refused, never
    /// clipped into an unphysical cold plasma or an extrapolated source domain.
    fn radiate(&self, old: Point, dt: f64, coefficient: f64) -> Result<Point, DtFusionError> {
        if coefficient == 0.0 || dt == 0.0 {
            return Ok(old);
        }
        let te = self.temperatures(old)?[1];
        let ce = self.capacities(old)[1];
        let loss = radiated_energy(te, ce, dt, coefficient)?;
        let mut next = old;
        next.photon_energy.add(loss);
        self.add_pool_energy(&mut next, [0.0, -loss]);
        Ok(next)
    }

    fn confine(
        &self,
        old: Point,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        let Some(confinement) = self.transport.energy_confinement else {
            return Ok(old);
        };
        let energies = self.pool_energy(old);
        let (remaining, outgoing) =
            confined_energy(energies, self.capacities(old), dt, confinement)?;
        let mut next = old;
        for i in 0..2 {
            next.transport_energy_out[i].add(outgoing[i]);
            if outgoing[i] > 0.5 * energies[i] {
                next.thermal_energy[i] = Sum::new(remaining[i]);
            } else {
                next.thermal_energy[i].add(-outgoing[i]);
            }
        }
        self.admit(next, cancelled)?;
        Ok(next)
    }
}
