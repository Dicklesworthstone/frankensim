//! Original-source thermal operators within the finite-fuel midpoint solve.

use super::*;
use fs_math::det::sqrt;

#[derive(Clone, Copy)]
pub(super) struct Sources {
    pub temperatures: [f64; 2],
    conductance_w_k: f64,
    radiation_w_per_sqrt_k: f64,
}

impl Problem<'_> {
    pub(super) fn sources(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sources, DtFusionError> {
        let temperatures = self.admit(point, cancelled)?;
        let [th, te] = temperatures.map(Temperature::new);
        let amounts = self.amounts(point);
        let volume = self.initial.volume();
        let conductance_w_k = if let Some(model) = &self.transport.coulomb {
            model
                .evaluate_point(&amounts, volume, th, te, &mut *cancelled)
                .map_err(coulomb_error)?
                .electron_heavy_thermal_conductance_w_m3_k()
                * volume.value()
        } else {
            0.0
        };
        let radiation_w_per_sqrt_k = if let Some(model) = &self.transport.bremsstrahlung {
            model
                .evaluate_point(&amounts, volume, th, te, &mut *cancelled)
                .map_err(radiation_error)?
                .total_power_over_sqrt_temperature_w_k_sqrt()
        } else {
            0.0
        };
        if [conductance_w_k, radiation_w_per_sqrt_k]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(DtFusionError::Arithmetic(
                "finite thermal source coefficients",
            ));
        }
        poll(cancelled)?;
        Ok(Sources {
            temperatures,
            conductance_w_k,
            radiation_w_per_sqrt_k,
        })
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
        let [th, te] = self.temperatures(old)?;
        let [ch, ce] = self.capacities(old);
        let exposure = (conductance / ch + conductance / ce) * dt;
        if exposure.is_nan() || exposure < 0.0 {
            return Err(DtFusionError::Arithmetic("Coulomb thermal exposure"));
        }
        // ch*ce/(ch+ce), without multiplying two extensive capacities.
        let reduced_capacity = ch.min(ce) / (1.0 + ch.min(ce) / ch.max(ce));
        let transfer = reduced_capacity * (te - th) * (-expm1(-exposure));
        if !transfer.is_finite() {
            return Err(DtFusionError::Arithmetic(
                "reciprocal Coulomb energy transfer",
            ));
        }
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
        let root = sqrt(te);
        let decrement = (0.5 * coefficient / ce) * dt;
        if !decrement.is_finite() || decrement >= root {
            return Err(DtFusionError::OutsideDomain);
        }
        let loss = (coefficient * dt) * (root - 0.5 * decrement);
        if !loss.is_finite() || loss < 0.0 {
            return Err(DtFusionError::Arithmetic(
                "escaping free-free photon energy",
            ));
        }
        let mut next = old;
        next.photon_energy.add(loss);
        Ok(next)
    }
}

fn coulomb_error(error: CoulombErrorV1) -> DtFusionError {
    match error {
        CoulombErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::Coulomb(error),
    }
}

fn radiation_error(error: BremsstrahlungErrorV1) -> DtFusionError {
    match error {
        BremsstrahlungErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::Radiation(error),
    }
}
