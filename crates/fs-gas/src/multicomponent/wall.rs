//! A stationary wall with the incident trace's actual species composition.
//!
//! Composition is frozen only across this local inviscid wave; every later
//! face solve borrows its own current partial densities. Original NASA-9
//! sources and their energy/entropy differences remain authoritative.

use super::{CaloricDomain, Cell, Error, FaceFlux, MixtureGas, Primitive, poll, positive};
use crate::dynamics::wall::shared::{CaloricLaw, HeatCapacities, WallError, WallSolver};
use fs_math::det;
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

impl MixtureGas {
    pub(super) fn wall_flux(
        &self,
        inside: Cell<'_>,
        primitive: Primitive,
        outward_normal: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        if out.len() != self.species_count() + 2
            || !self.admissible_cell(inside)
            || !matches!(outward_normal, -1.0 | 1.0)
            || !positive(primitive.density)
            || !positive(primitive.pressure)
            || !positive(primitive.temperature)
            || !positive(primitive.sound_speed)
            || !primitive.velocity.is_finite()
        {
            return Err(Error::Invalid(
                "wall normal, interior trace or output dimension",
            ));
        }
        let (density, concentration) = self.admit_densities(inside.partial_densities)?;
        let law = LocalCaloric {
            gas: self,
            densities: inside.partial_densities,
            density,
            gas_constant: (concentration / density) * R,
        };
        if !positive(law.gas_constant) {
            return Err(Error::Invalid("unrepresentable local mixture gas constant"));
        }
        let wall = WallSolver::new(&law)
            .stopped_wall(primitive, outward_normal * primitive.velocity, cancelled)
            .map_err(wall_error)?;
        if !positive(wall.density) || !positive(wall.temperature) {
            return Err(Error::Invalid(
                "unrepresentable wall density or temperature",
            ));
        }
        poll(cancelled)?;
        out.fill(0.0);
        out[self.species_count()] = wall.pressure;
        Ok(FaceFlux {
            signal_speed: wall.signal_speed,
            used_hll_fallback: false,
        })
    }
}

/// Borrows an already admitted owner and current trace; construction and all
/// successful evaluations allocate nothing and never rebuild a frozen owner.
struct LocalCaloric<'a> {
    gas: &'a MixtureGas,
    densities: &'a [f64],
    density: f64,
    gas_constant: f64,
}
impl LocalCaloric<'_> {
    fn temperature(&self, t: f64) -> Result<Temperature, Error> {
        let domain = self.gas.domain();
        if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&t) {
            return Err(Error::EnergyDomain);
        }
        Ok(Temperature::new(t))
    }
    fn weight(&self, i: usize) -> f64 {
        (self.densities[i] / self.density) * self.gas.inverse_molar_masses()[i]
    }
}
impl CaloricLaw for LocalCaloric<'_> {
    type Error = Error;
    fn domain(&self) -> CaloricDomain {
        self.gas.domain()
    }
    fn gas_constant(&self) -> f64 {
        self.gas_constant
    }
    fn heat_capacities(
        &self,
        t: f64,
        _pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeatCapacities, Error> {
        let temperature = self.temperature(t)?;
        let (mut cp, mut cp_correction, mut cv, mut cv_correction) = (0.0, 0.0, 0.0, 0.0);
        for (i, model) in self.gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            if self.densities[i] == 0.0 {
                continue;
            }
            let value = model.evaluate_properties(temperature)?.cp().value();
            let weight = self.weight(i);
            accumulate(&mut cp, &mut cp_correction, weight * value);
            accumulate(&mut cv, &mut cv_correction, weight * (value - R));
        }
        if !positive(cp) || !positive(cv) {
            return Err(Error::Invalid(
                "unrepresentable local mixture heat capacity",
            ));
        }
        Ok(HeatCapacities { cp, cv })
    }
    fn internal_energy_change(
        &self,
        initial_t: f64,
        final_t: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let initial_t = self.temperature(initial_t)?;
        let final_t = self.temperature(final_t)?;
        let (mut value, mut correction) = (0.0, 0.0);
        for (i, model) in self.gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            if self.densities[i] != 0.0 {
                accumulate(
                    &mut value,
                    &mut correction,
                    self.weight(i)
                        * model
                            .molar_internal_energy_change(initial_t, final_t)?
                            .value(),
                );
            }
        }
        if !value.is_finite() {
            return Err(Error::Invalid(
                "unrepresentable local mixture energy change",
            ));
        }
        Ok(value)
    }
    fn entropy_change(
        &self,
        initial_t: f64,
        initial_p: f64,
        final_t: f64,
        final_p: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let initial_t = self.temperature(initial_t)?;
        let final_t = self.temperature(final_t)?;
        if !positive(initial_p) || !positive(final_p) {
            return Err(Error::Invalid("invalid local mixture entropy pressure"));
        }
        let (mut value, mut correction) = (0.0, 0.0);
        for (i, model) in self.gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            if self.densities[i] != 0.0 {
                accumulate(
                    &mut value,
                    &mut correction,
                    self.weight(i) * model.molar_entropy_change(initial_t, final_t)?.value(),
                );
            }
        }
        // The ideal mixing entropy cancels at this fixed local composition.
        accumulate(
            &mut value,
            &mut correction,
            -self.gas_constant * log_ratio(final_p, initial_p),
        );
        if !value.is_finite() {
            return Err(Error::Invalid(
                "unrepresentable local mixture entropy change",
            ));
        }
        Ok(value)
    }
    fn sound_speed(
        &self,
        density: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let t = (pressure / density) / self.gas_constant;
        let capacities = self.heat_capacities(t, pressure, cancelled)?;
        let square = (capacities.cp / capacities.cv) * (pressure / density);
        if !positive(square) {
            return Err(Error::Invalid("unrepresentable local mixture sound speed"));
        }
        Ok(det::sqrt(square))
    }
    fn next_seam(
        &self,
        start: f64,
        mut end: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        for (i, model) in self.gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            if self.densities[i] == 0.0 {
                continue;
            }
            for region in model.regions() {
                let seam = region.temperature_min().value();
                if seam > start && seam < end {
                    end = seam;
                }
            }
        }
        Ok(end)
    }
}

fn wall_error(error: WallError<Error>) -> Error {
    match error {
        WallError::Source(error) => error,
        WallError::Invalid(reason) => Error::Invalid(reason),
        WallError::CaloricPath(reason) => Error::WallCaloricPath(reason),
        WallError::Iterations => Error::TemperatureIterations,
        WallError::Cancelled => Error::Cancelled,
    }
}
fn accumulate(total: &mut f64, correction: &mut f64, term: f64) {
    let y = term - *correction;
    let next = *total + y;
    *correction = (next - *total) - y;
    *total = next;
}
// Preserve weak-wave pressure changes and avoid overflow of a/b.
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    sign * if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        det::ln(rounded) * (delta / (rounded - 1.0))
    } else {
        det::ln(hi) - det::ln(lo)
    }
}
