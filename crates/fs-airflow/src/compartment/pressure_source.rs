//! Source operators inside the existing coupled opening-pressure equation.
//!
//! Every residual, Jacobian and line-search evaluation starts from the accepted
//! room inventories, applies transport, then evaluates the source on that actual
//! transported candidate. Opening pressure therefore closes after chemistry or
//! another finite source. Donor temperature, composition and interface geometry
//! in the vent flux remain frozen at the start of the interval. This is still a
//! first-order estimated composition, not a fully implicit thermochemical flux,
//! second-order integration method or discrete entropy theorem.

use super::{
    LayerRates, MAX_ROOMS, NetworkControl, NetworkError, NetworkStepReport, Room, Time,
    ZoneNetwork, poll, positive,
};

#[cfg(test)]
mod tests;

impl ZoneNetwork {
    /// Close opening pressure against transport followed by the supplied source.
    ///
    /// `prepare` runs once per attempted duration, from accepted rooms, and
    /// supplies the explicit transport sources plus immutable preparation data.
    /// `apply` runs for **every** nonlinear pressure candidate, including finite
    /// differences and line searches. It receives the actual transported rooms
    /// and must preserve their original source, geometry and physical clocks.
    /// The original reaction owner's `apply_split` satisfies that clock rule.
    ///
    /// Add the source's actual aggregate species and referenced-energy rates to
    /// `additional` exactly once. Do not repeat the transport or prepared rates.
    /// Pressure residuals are then recomputed from the source-updated rooms.
    ///
    /// Both callbacks must leave external physical owners unchanged. Return
    /// candidate wall states or diagnostics as `D`; only the value associated
    /// with the converged pressure candidate is returned. No `Clone` bound or
    /// speculative external publication is required. A failed trial discards
    /// its diagnostics, and a whole-step retry rebuilds the preparation data.
    ///
    /// The result is `(balance, preparation, accepted_source_result)`. The
    /// pressure residual includes the supplied source, and the reported
    /// post-pressure-solve split increment is zero because no source is applied
    /// after the converged candidate. Cancellation or exhaustion publishes no
    /// partial room state or payload.
    pub fn step_with_pressure_sources<P, D>(
        &mut self,
        maximum_duration: Time,
        controls: NetworkControl,
        mut prepare: impl FnMut(
            &[Room],
            Time,
            &mut [[LayerRates; 2]],
            &mut dyn FnMut() -> bool,
        ) -> Result<P, NetworkError>,
        mut apply: impl FnMut(
            &mut [Room],
            Time,
            &P,
            &mut LayerRates,
            &mut dyn FnMut() -> bool,
        ) -> Result<D, NetworkError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkStepReport, P, D), NetworkError> {
        poll(&mut cancelled)?;
        let n = self.rooms.len();
        let mut dt = maximum_duration.value();
        if !positive(dt)
            || !(1..=128).contains(&controls.maximum_iterations)
            || controls.maximum_halvings > 32
            || !positive(controls.minimum_step_s)
            || !positive(controls.pressure_tolerance_pa)
            || !positive(controls.pressure_difference_pa)
        {
            return Err(NetworkError::Invalid("network duration or controls"));
        }
        let mut last = NetworkError::PressureIterations;
        let mut attempts = 0;
        let mut source_storage = [[LayerRates::default(); 2]; MAX_ROOMS];
        for halving in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if dt < controls.minimum_step_s {
                break;
            }
            attempts += 1;
            let sources = &mut source_storage[..n];
            sources.fill([LayerRates::default(); 2]);
            for i in 0..n {
                self.pressures[i] = self.rooms[i].gas.pressure().value();
            }
            let candidate: Result<_, NetworkError> = (|| {
                let duration = Time::new(dt);
                let preparation = prepare(&self.rooms, duration, sources, &mut cancelled)?;
                let mut diagnostics = None;
                let (iterations, residual) = self.attempt_with_evaluation(
                    dt,
                    sources,
                    controls,
                    &mut cancelled,
                    |network, dt, sources, zone, cancelled| {
                        // A failed evaluation cannot leave diagnostics from a
                        // previous candidate eligible for publication.
                        diagnostics = None;
                        network.evaluate(dt, sources, zone, cancelled)?;
                        let mut additional = LayerRates::default();
                        let result = apply(
                            &mut network.trial,
                            Time::new(dt),
                            &preparation,
                            &mut additional,
                            cancelled,
                        )?;
                        network.admit_source_trial(dt, &additional)?;
                        network.external.add(&additional)?;
                        for (i, room) in network.trial.iter().enumerate() {
                            network.evaluated[i] =
                                room.gas.pressure().value() - network.perturbation[i];
                        }
                        poll(cancelled)?;
                        diagnostics = Some(result);
                        Ok(())
                    },
                )?;
                let result = diagnostics.ok_or(NetworkError::Invalid(
                    "converged pressure source has no candidate result",
                ))?;
                Ok((iterations, residual, preparation, result))
            })();
            match candidate {
                Ok((iterations, residual, preparation, result)) => {
                    let report =
                        self.publish_trial(dt, halving, iterations, residual, 0.0, &mut cancelled)?;
                    return Ok((report, preparation, result));
                }
                Err(NetworkError::Cancelled) => return Err(NetworkError::Cancelled),
                Err(error) => {
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(NetworkError::RetryLimit {
            attempts,
            last: Box::new(last),
        })
    }
}
