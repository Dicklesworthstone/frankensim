//! Exact compression histories of the same finite-rate regional water law.
//!
//! Positive compression e means V/V0 = 1-e for a fixed parcel mass. The
//! pressure perturbation is K_eq e + delta_K q, with
//! qdot = edot-q/tau. Storage per reference volume is
//! (K_eq e^2 + delta_K q^2)/2 and loss is delta_K q^2/tau.
//! The original Maxwell ramp evaluator integrates that storage, work and heat.
//!
//! Density is reconstructed from the parcel volume, rho=rho0/(1-e), rather
//! than advanced by a separate or clipped density equation. This geometric
//! identity does not turn the constitutive tangent into a nonlinear IF97 EOS.
//! The caller supplies a compression bound appropriate to the retained
//! coexistence state and declared kinetics. Thermal state, phase inventory,
//! reference pressure and coefficients remain fixed. Work and energy below
//! are perturbation quantities per ORIGINAL volume; background p0 work and
//! deposited heat must be handled by the hosting thermodynamic system.

use super::RelaxingRegionalWater;
use crate::MaterialError;
use crate::visco::PronyRamp;

fn invalid(what: &'static str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

/// Accepted parcel geometry, compression memory and acoustic storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCompressionState {
    /// Positive relative volume decrease, (V0-V)/V0.
    pub compression: f64,
    /// Elastic strain of the positive Maxwell compression arm.
    pub branch_compression: f64,
    /// Increment above the retained thermal state's pressure [Pa].
    pub pressure_perturbation_pa: f64,
    pub density_kg_m3: f64,
    /// Stable rho-rho0, retained even when absolute rho rounds to rho0.
    pub density_perturbation_kg_m3: f64,
    pub specific_volume_m3_kg: f64,
    pub stored_energy_j_m3: f64,
}

/// One accepted exact linear-compression ramp, including physical heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCompressionStep {
    pub initial: WaterCompressionState,
    pub state: WaterCompressionState,
    pub duration_s: f64,
    pub compression_work_j_m3: f64,
    pub dissipated_energy_j_m3: f64,
    /// Work minus the actual change in storage minus independent arm heat.
    pub energy_residual_j_m3: f64,
}

/// Prepared exponential coefficients, reusable by histories with the same tau.
/// Preparation has no source query or state change; application is scalar work.
#[derive(Clone, Copy, Debug)]
pub struct PreparedWaterCompressionRamp {
    duration_s: f64,
    relaxation_time_s: f64,
    ramp: PronyRamp,
}

impl PreparedWaterCompressionRamp {
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
}

/// Local, source-retaining constitutive history with atomic accepted updates.
/// No wave solver, spatial mass transport or second time integrator is added.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCompressionHistory {
    source: RelaxingRegionalWater,
    maximum_abs_compression: f64,
    state: WaterCompressionState,
}

impl WaterCompressionHistory {
    /// Start at the original thermal state with relaxed compression memory.
    ///
    /// The finite positive bound must be below one, ensuring positive parcel
    /// volume for every admitted ramp. Choosing it small enough to preserve
    /// the original phase/linearization regime remains the caller's duty.
    pub fn try_new(
        source: RelaxingRegionalWater,
        maximum_abs_compression: f64,
    ) -> Result<Self, MaterialError> {
        if !maximum_abs_compression.is_finite()
            || maximum_abs_compression <= 0.0
            || maximum_abs_compression >= 1.0
        {
            return Err(invalid(
                "water compression requires an explicit bound in (0,1)",
            ));
        }
        let state = Self::sample(&source, 0.0, 0.0)?;
        Ok(Self {
            source,
            maximum_abs_compression,
            state,
        })
    }

    pub const fn source(&self) -> &RelaxingRegionalWater {
        &self.source
    }

    pub const fn maximum_abs_compression(&self) -> f64 {
        self.maximum_abs_compression
    }

    pub const fn state(&self) -> WaterCompressionState {
        self.state
    }

    /// Prepare once for a shared duration and compression relaxation time.
    /// A zero duration denotes an instantaneous frozen elastic jump.
    pub fn prepare_ramp(
        &self,
        duration_s: f64,
    ) -> Result<PreparedWaterCompressionRamp, MaterialError> {
        let tau = self.source.arm().relaxation_time_s;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(invalid(
                "water compression duration must be finite and nonnegative",
            ));
        }
        let ratio = duration_s / tau;
        if !ratio.is_finite() || (duration_s > 0.0 && ratio == 0.0) {
            return Err(invalid(
                "water compression duration is unrepresentable relative to tau",
            ));
        }
        Ok(PreparedWaterCompressionRamp {
            duration_s,
            relaxation_time_s: tau,
            ramp: PronyRamp::new(ratio),
        })
    }

    /// Advance one prescribed linear compression ramp, atomically.
    pub fn advance(
        &mut self,
        compression: f64,
        duration_s: f64,
    ) -> Result<WaterCompressionStep, MaterialError> {
        let prepared = self.prepare_ramp(duration_s)?;
        self.advance_prepared(compression, &prepared)
    }

    /// Apply prepared exact-ramp coefficients without allocation or exponential
    /// evaluation. Every validation precedes publication of the next state.
    ///
    /// Positive heat is integrated independently from the work/storage balance.
    /// Zero-duration jumps add frozen elastic energy and exactly zero heat.
    /// Reversals, holds and recovery all use the same physical memory equation.
    pub fn advance_prepared(
        &mut self,
        compression: f64,
        prepared: &PreparedWaterCompressionRamp,
    ) -> Result<WaterCompressionStep, MaterialError> {
        if !compression.is_finite() || compression.abs() > self.maximum_abs_compression {
            return Err(invalid(
                "water compression exceeds the declared linearization bound",
            ));
        }
        let arm = self.source.arm();
        if prepared.relaxation_time_s != arm.relaxation_time_s {
            return Err(invalid(
                "water compression ramp belongs to a different relaxation time",
            ));
        }
        let initial = self.state;
        let delta = compression - initial.compression;
        let (memory, mean_memory, heat_per_modulus) =
            prepared.ramp.apply(initial.branch_compression, delta);
        let next = Self::sample(&self.source, compression, memory)?;
        let equilibrium = self
            .source
            .equilibrium()
            .state()
            .adiabatic_bulk_modulus_pa();
        let mean_pressure = equilibrium * f64::midpoint(initial.compression, compression)
            + arm.delta_bulk_modulus_pa * mean_memory;
        let work = mean_pressure * delta;
        let heat = arm.delta_bulk_modulus_pa * heat_per_modulus;
        let residual = work - (next.stored_energy_j_m3 - initial.stored_energy_j_m3) - heat;
        let scale = initial.stored_energy_j_m3 + next.stored_energy_j_m3 + work.abs() + heat;
        if ![work, heat, residual, scale]
            .into_iter()
            .all(f64::is_finite)
            || heat < 0.0
            || residual.abs() > 128.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "water compression work, storage or heat is unrepresentable",
            ));
        }
        let step = WaterCompressionStep {
            initial,
            state: next,
            duration_s: prepared.duration_s,
            compression_work_j_m3: work,
            dissipated_energy_j_m3: heat,
            energy_residual_j_m3: residual,
        };
        self.state = next;
        Ok(step)
    }

    fn sample(
        source: &RelaxingRegionalWater,
        compression: f64,
        memory: f64,
    ) -> Result<WaterCompressionState, MaterialError> {
        let acoustic = source.equilibrium().state();
        let equilibrium = acoustic.adiabatic_bulk_modulus_pa();
        let strength = source.arm().delta_bulk_modulus_pa;
        let volume_ratio = 1.0 - compression;
        let density = acoustic.density_kg_m3() / volume_ratio;
        let density_change = acoustic.density_kg_m3() * (compression / volume_ratio);
        let specific_volume = source.source().specific_volume_m3_kg() * volume_ratio;
        let pressure = equilibrium * compression + strength * memory;
        let stored =
            0.5 * (equilibrium * compression) * compression + 0.5 * (strength * memory) * memory;
        if ![
            compression,
            memory,
            density,
            density_change,
            specific_volume,
            pressure,
            stored,
        ]
        .into_iter()
        .all(f64::is_finite)
            || density <= 0.0
            || specific_volume <= 0.0
            || stored < 0.0
        {
            return Err(invalid(
                "water compression state or fixed-mass geometry is unrepresentable",
            ));
        }
        Ok(WaterCompressionState {
            compression,
            branch_compression: memory,
            pressure_perturbation_pa: pressure,
            density_kg_m3: density,
            density_perturbation_kg_m3: density_change,
            specific_volume_m3_kg: specific_volume,
            stored_energy_j_m3: stored,
        })
    }
}

#[cfg(test)]
#[path = "transient/tests.rs"]
mod tests;
