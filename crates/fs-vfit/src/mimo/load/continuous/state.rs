//! Source-preserving adoption of an original rational trajectory.
use super::MimoContinuous;
use crate::mimo::load::{MimoLoadError, PreparedMimoLoad};
use crate::model::PoleTerm;
use fs_math::c64::C64;

impl MimoContinuous<'_> {
    /// Accepted coordinates in the original `MimoModel::state_space` chart.
    /// Pair blocks contain all Re(z), then all MINUS Im(z); the source ZOH
    /// owner stores the same poles as complex z in term/input-port order.
    pub fn real_coordinates(&self) -> Vec<f64> {
        let p = self.load.ports();
        let mut values = Vec::with_capacity(self.dimension());
        for (term, states) in self
            .load
            .source_model()
            .poles
            .iter()
            .zip(self.load.pole_states().chunks_exact(p))
        {
            values.extend(states.iter().map(|z| z.re));
            if matches!(term, PoleTerm::Pair { .. }) {
                values.extend(states.iter().map(|z| -z.im));
            }
        }
        values
    }

    /// Return a candidate of THIS admitted owner, preserving the source, clock,
    /// reports and every original pole. Publication belongs to the enclosing
    /// transaction. Finite endpoints are not a passivity or storage certificate;
    /// the enclosing operator must independently account its actual port work.
    pub fn candidate_from_coordinates(
        &self,
        endpoint: &[f64],
        endpoint_input: &[f64],
    ) -> Result<PreparedMimoLoad, MimoLoadError> {
        let p = self.load.ports();
        if endpoint.len() != self.dimension()
            || endpoint_input.len() != p
            || endpoint
                .iter()
                .chain(endpoint_input)
                .any(|v| !v.is_finite())
        {
            return Err(MimoLoadError(
                "continuous MIMO endpoint must retain finite complete original coordinates and ports",
            ));
        }
        let mut candidate = self.load.clone();
        let mut at = 0;
        for (index, term) in self.load.source_model().poles.iter().enumerate() {
            let pair = matches!(term, PoleTerm::Pair { .. });
            for port in 0..p {
                candidate.runtime.state[index * p + port] = C64::new(
                    endpoint[at + port],
                    if pair { -endpoint[at + p + port] } else { 0.0 },
                );
            }
            at += if pair { 2 * p } else { p };
        }
        candidate
            .runtime
            .candidate_state
            .clone_from(&candidate.runtime.state);
        candidate.runtime.last_input.copy_from_slice(endpoint_input);
        Ok(candidate)
    }
}
