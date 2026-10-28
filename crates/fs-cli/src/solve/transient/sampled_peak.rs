//! The selected region's maximum over accepted fine-grid states.
//!
//! The original initial state is sample zero. This is a forward observation:
//! it neither requests an adjoint nor supplies a continuous-time bound.

use super::{SolveRefusal, bad, json_string, number};

pub(super) struct SampledPeak {
    value_k: f64,
    time_s: f64,
    state: usize,
    samples: usize,
    last_time_s: f64,
}

impl SampledPeak {
    pub(super) fn new(initial_k: f64) -> Result<Self, SolveRefusal> {
        if !initial_k.is_finite() || initial_k <= 0.0 {
            return Err(bad("sampled thermal peak needs a positive finite initial temperature"));
        }
        Ok(Self { value_k: initial_k, time_s: 0.0, state: 0, samples: 1, last_time_s: 0.0 })
    }

    pub(super) fn record(&mut self, time_s: f64, value_k: f64) -> Result<(), SolveRefusal> {
        if !time_s.is_finite() || time_s <= self.last_time_s
            || !value_k.is_finite() || value_k <= 0.0
        {
            return Err(bad("sampled thermal peak lost its increasing accepted time grid"));
        }
        let next = self.samples.checked_add(1)
            .ok_or_else(|| bad("sampled thermal peak count overflows"))?;
        // Strict comparison retains the earliest sample for exact ties.
        if value_k > self.value_k {
            self.value_k = value_k;
            self.time_s = time_s;
            self.state = self.samples;
        }
        self.samples = next;
        self.last_time_s = time_s;
        Ok(())
    }

    pub(super) fn receipt(
        &self, region: &str, fine_steps: usize, final_time_s: f64,
    ) -> Result<String, SolveRefusal> {
        if region.is_empty() || fine_steps.checked_add(1) != Some(self.samples)
            || self.last_time_s != final_time_s
        {
            return Err(bad("sampled thermal peak does not cover the selected complete trajectory"));
        }
        Ok(format!(
            "{{\"region\":{},\"value_k\":{},\"time_s\":{},\"state\":{},\"samples\":{},\"includes_initial\":true,\"time_scope\":\"accepted-fine-grid\",\"continuous_time_bound\":false,\"authority\":\"Estimated\"}}",
            json_string(region), number(self.value_k)?, number(self.time_s)?,
            self.state, self.samples,
        ))
    }
}
