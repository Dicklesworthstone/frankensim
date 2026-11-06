//! Original volume-stress arrays inside a finite layer, radiating into an exterior.
//!
//! Each physical source caches the native inside-source 4-by-6 pressure/velocity
//! response. The shared stress-array owner retains all cell, tensor and internal-
//! path cross terms before coherent or statistically specified field reduction.

use super::PreparedInteriorSourceFullWaveQuadrupole;
use crate::AeroacError;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::curle3d::refraction::full_wave::layer::{
    FullWaveLayerQuadrupoleDiagnostics, PlanarFluidLayer,
};
use crate::curle3d::refraction::full_wave::quadrupole::array::{
    ArrayField, FullWaveQuadrupoleSourceLimits, Transfer,
};
use crate::curle3d::refraction::full_wave::{QuadratureLimits, poll, positive};
use crate::curle3d::statistics::invalid;
use crate::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use crate::observation::HarmonicFieldLine;
use fs_fft::welch::WelchCrossSpectrum;
use fs_math::c64::C64;

mod recorded;

/// Applicability of the prescribed inside-layer stress array and its exterior field.
pub const FULL_WAVE_INTERIOR_SOURCE_QUADRUPOLE_ARRAY_SCOPE: &str = "Original prescribed \
stationary symmetric volume-stress moments [N m] strictly inside one infinite planar \
homogeneous lossless fluid layer; coherent full-wave pressure and actual particle \
velocity strictly in the original incident exterior. All tensor components, internal \
returns, critical, reactive and evanescent fields retain the original source clock. \
Source statistics are explicit; slower-than-both trapping layers refuse. Per-source \
quadrature estimates and spectral-tail bounds give no relative array guarantee at \
interference nulls. No source-size/SGS validation, source loading, zero-thickness \
interior, face trace, absorption, moving fluids, finite edges or interval certificate.";

/// Prepared original inside-layer stress cells observed before the first face.
/// One bounded transfer table shares a preparation work cap; warm fields and
/// recorded factors use the common allocation-free fixed-order DD reductions.
#[derive(Clone, Debug)]
pub struct PreparedInteriorSourceFullWaveQuadrupoleObserver {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<FullWaveLayerQuadrupoleDiagnostics>>,
    evaluations: usize,
}

impl PreparedInteriorSourceFullWaveQuadrupoleObserver {
    /// Every source lies strictly between the original faces; the receiver lies
    /// strictly before the first face. Empty arrays still admit a positive layer
    /// thickness, exterior receiver, frequency and numerical controls. To observe
    /// the other exterior, explicitly reverse the original faces and exterior
    /// states while retaining physical source positions and Cartesian moments.
    ///
    /// An absent second interface reduces to the original source-fluid interface.
    /// Equal speeds alone do not remove density-layer returns. Absolute native
    /// error estimates multiply each physical stress Frobenius norm and add;
    /// they imply no relative guarantee after destructive cell interference.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits.quadrupole).admit()?;
        if source_positions_m.len() > limits.max_sources || limits.maximum_total_evaluations == 0 {
            return Err(invalid("inside-layer stress array source or work budget"));
        }
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|value| !value.is_finite())
        {
            return Err(invalid("inside-layer stress array frequency or observer"));
        }
        if !positive(layer.thickness_m())
            || !positive((-layer.exterior.distance(observer_m)).to_f64())
        {
            return Err(invalid(
                "inside-layer stress array requires positive thickness and an exterior observer before the first face",
            ));
        }
        let exact = layer.layer_state() == layer.receiving_state()
            && layer.layer_state().sound_speed_m_s() == layer.incident_state().sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("inside-layer stress array allocation"))?;
        let mut evaluations = 0usize;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            let remaining = limits.maximum_total_evaluations - evaluations;
            if !exact && remaining == 0 {
                return Err(invalid(
                    "inside-layer stress array exhausted aggregate work",
                ));
            }
            let mut source_limits = limits.quadrupole;
            if !exact {
                source_limits.maximum_evaluations =
                    source_limits.maximum_evaluations.min(remaining);
            }
            let original = PreparedInteriorSourceFullWaveQuadrupole::new(
                layer,
                frequency_hz,
                source_m,
                observer_m,
                source_limits,
                cancelled,
            )?;
            // Physical unit moments extract the already prepared native map.
            // No quadrature, differentiation or alternative shear convention.
            let mut field = [[C64::ZERO; 6]; 4];
            for column in 0..6 {
                let mut unit = QuadrupolePhasor::default();
                unit.stress_moment_n_m[column] = C64::ONE;
                let line = original.evaluate_harmonic_line(unit)?;
                field[0][column] = line.pressure_pa;
                for axis in 0..3 {
                    field[axis + 1][column] = line.particle_velocity_m_s[axis];
                }
            }
            let diagnostics = *original.diagnostics();
            evaluations = evaluations
                .checked_add(diagnostics.evaluations)
                .filter(|&work| work <= limits.maximum_total_evaluations)
                .ok_or_else(|| invalid("inside-layer stress array exceeded aggregate work"))?;
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            layer,
            frequency_hz,
            observer_m,
            transfers,
            evaluations,
        })
    }

    pub const fn layer(&self) -> PlanarFluidLayer {
        self.layer
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|entry| entry.source_m)
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    pub const fn scope(&self) -> &'static str {
        FULL_WAVE_INTERIOR_SOURCE_QUADRUPOLE_ARRAY_SCOPE
    }
    /// Original inside-source transfer estimates in their physical stress units.
    pub fn source_diagnostics(&self, index: usize) -> Option<&FullWaveLayerQuadrupoleDiagnostics> {
        self.transfers.get(index).map(|entry| &entry.diagnostics)
    }

    pub fn evaluate(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }
    pub fn evaluate_cancellable(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate_cancellable(sources, cancelled)
    }
    /// Original peak-source clock, with the complete native propagation phase.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }
    /// Explicit cell independence retains coherent stress and returns per cell.
    pub fn evaluate_independent(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_independent(sources, cancelled)
    }
    /// Each supplied peak-moment mode retains all cell and tensor cross terms.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[QuadrupolePhasor]],
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view()
            .evaluate_coherence_modes(modes, budget, cancelled)
    }
    /// Joint source-major six-channel Welch factors, with native bin normalization.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view()
            .evaluate_welch_bin(spectrum, bin, budget, cancelled)
    }

    fn view(&self) -> ArrayField<'_, FullWaveLayerQuadrupoleDiagnostics> {
        ArrayField {
            transfers: &self.transfers,
            frequency_hz: self.frequency_hz,
            scope: self.scope(),
        }
    }
}
