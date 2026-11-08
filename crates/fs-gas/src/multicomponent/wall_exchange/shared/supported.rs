//! Separate original gas and wall caloric domains for constant/lagged films.

use super::*;

impl Group<'_> {
    /// A hot/cold finite wall can have a supported coupled endpoint even when
    /// holding its initial temperature would drive a gas outside its chart.
    /// Build a supported interval before any such nested reservoir trial.
    /// Ordinary in-domain and inactive contacts spend no extra source queries.
    pub(super) fn supported_heat_interval<G: WallGasCaloric>(
        &self,
        gas: &impl Fn(SharedWallGas, f64) -> G,
        wall_initial: f64,
        duration: Time,
        control: WallControl,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<[f64; 2]>, WallExchangeError> {
        let dt = duration.value();
        if !positive(dt) {
            return Ok(None);
        }
        #[cfg(not(feature = "material"))]
        let state_dependent = false;
        #[cfg(feature = "material")]
        let state_dependent = self
            .conduction
            .as_ref()
            .is_some_and(|c| c.implicit || c.integrated);
        let mut outside = false;
        for i in 0..self.accepted.ports.len() {
            let port = &self.accepted.ports[i];
            poll(cancelled)?;
            let conductance = port.conductance_w_per_k;
            #[cfg(feature = "material")]
            let conductance = match &self.conduction {
                Some(binding) if !state_dependent => binding.conductance(
                    i,
                    conductance,
                    self.initial_wall,
                    self.accepted.storage.heat_storage.as_ref(),
                    work,
                )?,
                _ => conductance,
            };
            if dt * conductance == 0.0 {
                continue;
            }
            let domain = gas(port.gas, wall_initial).domain();
            if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&wall_initial) {
                outside = true;
                break;
            }
        }
        if !outside {
            return Ok(None);
        }
        if state_dependent {
            return Err(WallExchangeError::Invalid(
                "shared wall outside gas interval requires constant or lagged conductance",
            ));
        }
        let (mut minimum, mut maximum) = (wall_initial, wall_initial);
        let (mut gas_low, mut gas_high) = (f64::NEG_INFINITY, f64::INFINITY);
        for (i, port) in self.accepted.ports.iter().enumerate() {
            poll(cancelled)?;
            let conductance = port.conductance_w_per_k;
            #[cfg(feature = "material")]
            let conductance = match &self.conduction {
                Some(binding) => binding.conductance(
                    i,
                    conductance,
                    self.initial_wall,
                    self.accepted.storage.heat_storage.as_ref(),
                    work,
                )?,
                None => conductance,
            };
            // InitialState reuses the already-prepared series coefficient;
            // authored films and source storage are never overwritten.
            let integrated = dt * conductance;
            if !nonnegative(integrated) {
                return Err(WallExchangeError::Invalid(
                    "shared wall integrated conductance",
                ));
            }
            if integrated == 0.0 {
                continue;
            }
            let original = MeteredGas {
                original: gas(port.gas, wall_initial),
                work,
            };
            let state = self.initial[i];
            let primitive = original.primitive(state[0], state[1], cancelled)?;
            let target = state[1] - 0.5 * state[0] * primitive.velocity;
            let temperature = refine_caloric_temperature(
                &original,
                primitive.temperature,
                target,
                self.volumes[i],
                integrated,
                wall_initial,
                0.0,
                control,
                cancelled,
            )?;
            minimum = minimum.min(temperature);
            maximum = maximum.max(temperature);
            let domain = original.domain();
            // For positive original Cv and constant G, the gas root lies in
            // [Tmin,Tmax] iff Tw lies between T + Delta E(T)/(dt G) at those
            // endpoints. Delta E uses the same represented kinetic subtraction
            // and volume as the actual contact equation, never a fitted Cv.
            let limit = |temperature| -> Result<f64, WallExchangeError> {
                let change = (original.caloric(temperature)?.0 - target) * self.volumes[i];
                if !change.is_finite() {
                    return Err(WallExchangeError::Invalid("shared wall gas energy bound"));
                }
                Ok(temperature + change / integrated)
            };
            gas_low = gas_low.max(limit(domain.minimum_temperature)?);
            gas_high = gas_high.min(limit(domain.maximum_temperature)?);
        }
        // Tiny G can make the mathematical bound infinite. Intersect with the
        // finite physical passive hull before admission; no gas is queried at
        // that wall-temperature bound or beyond its own original chart.
        let lower_temperature = minimum.max(gas_low);
        let upper_temperature = maximum.min(gas_high);
        if !lower_temperature.is_finite()
            || !upper_temperature.is_finite()
            || lower_temperature > upper_temperature
        {
            return Err(FlowError::EnergyDomain.into());
        }
        let (mut low, mut high) = if let Some(source) = &self.accepted.storage.heat_storage {
            let contact = WallContact {
                drag: WallDrag::None,
                thermal: self.initial_wall,
            };
            work.charge()?;
            let cold = source.heat_limits(contact, minimum)?;
            work.charge()?;
            let hot = source.heat_limits(contact, maximum)?;
            (-cold[1].max(hot[1]), -cold[0].min(hot[0]))
        } else {
            let WallThermal::Lumped {
                heat_capacity_j_per_k,
                ..
            } = self.initial_wall
            else {
                unreachable!("finite storage admitted");
            };
            (
                heat_capacity_j_per_k * (minimum - wall_initial),
                heat_capacity_j_per_k * (maximum - wall_initial),
            )
        };
        if low < 0.0 {
            low = low.next_up();
        }
        if high > 0.0 {
            high = high.next_down();
        }
        low = self.bound(
            self.initial_wall,
            low,
            minimum,
            true,
            control,
            work,
            cancelled,
        )?;
        high = self.bound(
            self.initial_wall,
            high,
            maximum,
            false,
            control,
            work,
            cancelled,
        )?;
        if !low.is_finite() || !high.is_finite() || low > high {
            return Err(WallExchangeError::Invalid(
                "shared wall supported heat interval",
            ));
        }
        // Both endpoints are in the original wall chart. Verify the opposite
        // endpoint before using it as an inside anchor: zero change may be
        // excluded from BOTH gas-supported temperature inequalities.
        if self.after(high, work)?.temperature < lower_temperature
            || self.after(low, work)?.temperature > upper_temperature
        {
            return Err(FlowError::EnergyDomain.into());
        }
        low = self.bound_from(
            self.initial_wall,
            low,
            high,
            lower_temperature,
            true,
            control,
            work,
            cancelled,
        )?;
        if self.after(low, work)?.temperature > upper_temperature {
            return Err(FlowError::EnergyDomain.into());
        }
        high = self.bound_from(
            self.initial_wall,
            high,
            low,
            upper_temperature,
            false,
            control,
            work,
            cancelled,
        )?;
        if low > high {
            return Err(FlowError::EnergyDomain.into());
        }
        Ok(Some([low, high]))
    }
}
