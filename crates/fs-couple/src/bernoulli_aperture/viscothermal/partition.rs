//! Geometry-only partition for the existing distributed-loss builder.
//! Every quarter uses the SAME section resolver and allpass phase as execution.

use super::{AcousticRealizeError, CancelGate, GasState, TubeSection,
    ViscothermalSectionSpec, checkpoint, invalid, wave_error};
use crate::bernoulli_aperture::fractional::{FractionalSection, FractionalTubePolicy};

pub(super) struct FractionalPartition {
    quarters: Vec<FractionalSection>,
    tolerance_m: f64,
    represented_m: f64,
    integer_samples: usize,
}
impl FractionalPartition {
    pub(super) fn new(section: &TubeSection, gas: &GasState, dt: f64,
        selection: &ViscothermalSectionSpec, policy: FractionalTubePolicy, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        let omega = policy.angular_frequency(dt)?;
        if selection.maximum_frequency_hz > policy.maximum_frequency_hz {
            return Err(invalid("fractional propagation band must cover every selected viscothermal band"));
        }
        // The caller has admitted 1..=128 cells and aggregate scratch storage.
        let count = 4*selection.cells;
        let tolerance_m = section.max_length_error_m/count as f64;
        let mut quarters = Vec::new();
        quarters.try_reserve_exact(count)
            .map_err(|_| invalid("fractional viscothermal partition allocation failed"))?;
        let mut represented_m = 0.0;
        let mut integer_samples = 0usize;
        let mut phase_error = 0.0;
        let mut left = 0.0;
        for k in 0..count {
            checkpoint(gate)?;
            // Split physical coordinates, NOT rounded L/(c dt). The final
            // boundary is exactly the original endpoint. All loads below use
            // these physical lengths, not a sample count multiplied by c dt.
            let right = if k+1 == count { section.length_m }
                else { section.length_m*((k+1) as f64/count as f64) };
            let quarter = TubeSection { nodes: section.nodes, length_m: right-left,
                radius_m: section.radius_m, max_length_error_m: tolerance_m };
            let (_, detail) = quarter.realize_fractional(gas, dt, omega)?;
            represented_m += detail.realized_dc_length_m;
            integer_samples = integer_samples.checked_add(detail.integer_samples)
                .ok_or_else(|| invalid("fractional viscothermal transit size overflow"))?;
            phase_error += detail.phase_error_bound_rad;
            quarters.push(detail);
            left = right;
        }
        if !represented_m.is_finite() || !section.max_length_error_m.is_finite()
            || section.max_length_error_m < 0.0
            || (represented_m-section.length_m).abs() > section.max_length_error_m {
            return Err(invalid("fractional viscothermal DC length exceeds the original geometry tolerance"));
        }
        if !phase_error.is_finite() || phase_error > policy.maximum_summed_phase_error_rad {
            return Err(invalid("fractional viscothermal quarters exceed the declared phase-error budget"));
        }
        Ok(Self { quarters, tolerance_m, represented_m, integer_samples })
    }
    pub(super) fn quarter_length_m(&self, k: usize) -> f64 {
        self.quarters[k].requested_length_m
    }
    pub(super) fn cell_length_m(&self, cell: usize) -> f64 {
        self.quarters[4*cell..4*cell+4].iter().map(|q| q.requested_length_m).sum()
    }
    pub(super) const fn represented_length_m(&self) -> f64 { self.represented_m }
    pub(super) const fn integer_samples(&self) -> usize { self.integer_samples }
    pub(super) const fn quarter_tolerance_m(&self) -> f64 { self.tolerance_m }
    pub(super) fn quarter_phase_rad(&self, k: usize, omega: f64)
        -> Result<f64, AcousticRealizeError>
    {
        let quarter = &self.quarters[k];
        Ok(omega*quarter.integer_samples as f64
            + quarter.transit.phase_lag_rad(omega).map_err(wave_error)?)
    }
    pub(super) fn cell_phase_rad(&self, cell: usize, omega: f64)
        -> Result<f64, AcousticRealizeError>
    {
        let mut phase = 0.0;
        for k in 4*cell..4*cell+4 { phase += self.quarter_phase_rad(k, omega)?; }
        Ok(phase)
    }
}
