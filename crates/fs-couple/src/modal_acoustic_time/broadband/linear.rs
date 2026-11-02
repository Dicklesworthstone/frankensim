//! Continuous receiver filtering along the original first-order-held force
//! trajectory. Only an extra force-moment column is retained after cold
//! preparation; the existing held-force operator and accepted histories stay
//! the owners of the common part of every step.

use super::*;

impl BroadbandModalObserver {
    /// Prepare exact linear-force observation at the original sample clock.
    ///
    /// The continuous filter is augmented with the physical force and its
    /// endpoint increment. In normalized time u=t/dt, dF/du=F_end-F_start.
    /// Thus preparation never divides a force difference by a tiny duration.
    /// All original filter states, modes, fit errors and propagation delays
    /// remain active. This cold operation allocates; [`Self::step_linear`]
    /// subsequently reuses the existing workspace without allocation.
    ///
    /// A constant real velocity-to-pressure gain is an ordinary feedthrough.
    /// A constant complex narrow-band quadrature is not a causal broadband
    /// transfer for an arbitrary ramp, and refuses this new path explicitly.
    pub fn prepare_linear_forcing(&mut self) -> Result<(), BroadbandModalError> {
        if self.linear_force_gains.is_some() {
            return Ok(());
        }
        if self
            .constant
            .as_ref()
            .is_some_and(|h| h.iter().any(|x| x.im != 0.0))
        {
            return Err(BroadbandModalError::Invalid(
                "linear forcing requires a continuous receiver or a real constant velocity gain",
            ));
        }
        let entries = self
            .filters
            .iter()
            .try_fold(0usize, |sum, filter| {
                let n = filter.n.checked_add(1)?;
                sum.checked_add(n.checked_mul(n)?)
            })
            .ok_or(BroadbandModalError::Invalid(
                "linear observer shape overflow",
            ))?;
        if entries > MAX_MODAL_AUGMENTED_MATRIX_ENTRIES {
            return Err(BroadbandModalError::Cap {
                what: "linear-force observer matrix entries",
                found: entries,
                cap: MAX_MODAL_AUGMENTED_MATRIX_ENTRIES,
            });
        }
        let dt = self.model.sample_period_s();
        let mut gains = Vec::with_capacity(self.filters.len());
        for (mode, filter) in self.model.modes().iter().zip(&self.filters) {
            let ss = filter.model.state_space();
            let n = ss.n + 4;
            let force = n - 2;
            let increment = n - 1;
            let mut a = vec![0.0; n * n];
            let w = mode.angular_frequency_rad_s;
            a[1] = filter.q_scale * dt;
            a[n] = -w * w / filter.q_scale * dt;
            a[n + 1] = -2.0 * mode.damping_ratio * w * dt;
            a[n + force] = dt;
            a[force * n + increment] = 1.0;
            for i in 0..ss.n {
                a[(i + 2) * n + 1] = ss.b[i] * dt;
                for j in 0..ss.n {
                    a[(i + 2) * n + j + 2] = ss.a[i * ss.n + j] * dt;
                }
            }
            let transition =
                fs_phs::linear_exact::exponential(&a, n).map_err(BroadbandModalError::Exact)?;
            gains.push(
                (0..ss.n)
                    .map(|i| transition[(i + 2) * n + increment])
                    .collect(),
            );
        }
        self.model
            .prepare_linear_forcing()
            .map_err(BroadbandModalError::Modal)?;
        self.linear_force_gains = Some(gains);
        Ok(())
    }

    /// Advance mechanics and continuous receiver under a linearly interpolated
    /// physical force record. Preparation is explicit; this warm method does
    /// not allocate, fit a transfer, or evaluate a matrix exponential.
    ///
    /// Pressure, mechanical, force and energy gates precede publication of
    /// the mechanical state, receiver memory and delay history together.
    /// The receiver remains one-way: it contributes no mechanical heat/load.
    /// Delay interpolation retains its original sampled, band-limited scope.
    pub fn step_linear(
        &mut self,
        force_start: &[f64],
        force_end: &[f64],
    ) -> Result<f64, BroadbandModalError> {
        if force_start.len() != self.model.modes().len()
            || force_end.len() != force_start.len()
            || force_start.iter().chain(force_end).any(|f| !f.is_finite())
        {
            return Err(BroadbandModalError::Invalid(
                "finite endpoint forces required for every retained mode",
            ));
        }
        let gains = self
            .linear_force_gains
            .as_ref()
            .ok_or(BroadbandModalError::Invalid(
                "prepare linear-force observation before the warm step",
            ))?;
        if force_start == force_end && self.constant.is_none() {
            return self.step(force_start);
        }
        let mut raw = 0.0;
        if let Some(h) = &self.constant {
            for (i, &gain) in h.iter().enumerate() {
                let next = self
                    .model
                    .advance_exact_linear_step(
                        i,
                        self.model.states()[i],
                        force_start[i],
                        force_end[i],
                    )
                    .map_err(BroadbandModalError::Modal)?;
                raw += gain.re * next.velocity_m_sqrt_kg_per_s;
            }
        } else {
            for (i, (filter, ramp)) in self.filters.iter_mut().zip(gains).enumerate() {
                let old = self.model.states()[i];
                let next = self
                    .model
                    .advance_exact_linear_step(i, old, force_start[i], force_end[i])
                    .map_err(BroadbandModalError::Modal)?;
                let difference = force_end[i] - force_start[i];
                for (row, &gain) in ramp.iter().enumerate() {
                    let r = (row + 2) * filter.n;
                    let mut value =
                        filter.transition[r] * filter.q_scale * old.displacement_m_sqrt_kg
                            + filter.transition[r + 1] * old.velocity_m_sqrt_kg_per_s;
                    for (j, &state) in filter.state.iter().enumerate() {
                        value += filter.transition[r + j + 2] * state;
                    }
                    value += filter.transition[r + filter.n - 1] * force_start[i];
                    value += if difference.is_finite() {
                        gain * difference
                    } else {
                        gain * force_end[i] - gain * force_start[i]
                    };
                    if !value.is_finite() {
                        return Err(BroadbandModalError::NonFinite {
                            what: "linear-force continuous filter state",
                        });
                    }
                    filter.candidate[row] = value;
                }
                let mut pressure = filter.d * next.velocity_m_sqrt_kg_per_s;
                for (&c, &state) in filter.c.iter().zip(&filter.candidate) {
                    pressure += c * state;
                }
                raw += pressure;
            }
        }
        if !raw.is_finite() {
            return Err(BroadbandModalError::NonFinite {
                what: "linear-force pressure superposition",
            });
        }
        let pressure = self.delay.preview(raw)?;
        if pressure.abs() > self.model.budget().maximum_abs_pressure_pa {
            return Err(BroadbandModalError::Invalid(
                "linear-force observer exceeds the declared pressure ceiling",
            ));
        }
        self.model
            .step_linear_into(force_start, force_end, &mut self.workspace)
            .map_err(BroadbandModalError::Modal)?;
        for filter in &mut self.filters {
            filter.state.copy_from_slice(&filter.candidate);
        }
        self.delay.commit(raw);
        Ok(pressure)
    }
}

#[cfg(test)]
mod tests;
