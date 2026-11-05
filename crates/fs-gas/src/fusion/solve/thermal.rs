//! Original-source thermal operators within the finite-fuel midpoint solve.

use super::super::thermal::{evaluate, exchange_energy, radiated_energy};
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
    /// map itself symmetrically composes analytic free-free cooling and exact
    /// reciprocal exchange. Each thermal stage uses its actual heat capacities.
    pub(super) fn frozen_coupled_step(
        &self,
        old: Point,
        dt: f64,
        sources: Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        let first = self.thermal_step(old, 0.5 * dt, sources, cancelled)?;
        let burned = self.frozen_step(first, dt, sources.temperatures[0])?;
        self.admit(burned, cancelled)?;
        let next = self.thermal_step(burned, 0.5 * dt, sources, cancelled)?;
        self.admit(next, cancelled)?;
        Ok(next)
    }

    fn thermal_step(
        &self,
        old: Point,
        dt: f64,
        sources: Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        poll(cancelled)?;
        if sources.conductance_w_k == 0.0 && sources.radiation_w_per_sqrt_k == 0.0 {
            return Ok(old);
        }
        let first = self.radiate(old, 0.5 * dt, sources.radiation_w_per_sqrt_k)?;
        self.admit(first, cancelled)?;
        let exchanged = self.exchange(first, dt, sources.conductance_w_k)?;
        self.admit(exchanged, cancelled)?;
        let next = self.radiate(exchanged, 0.5 * dt, sources.radiation_w_per_sqrt_k)?;
        self.admit(next, cancelled)?;
        Ok(next)
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
        Ok(next)
    }
}
