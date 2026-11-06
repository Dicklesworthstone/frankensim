//! Original volume-stress arrays in all supported finite-layer observation domains.
//!
//! Each cell caches the native complete 4-by-6 pressure/velocity response. The
//! original stress-array reducer owns coherent, independent, low-rank and joint
//! Welch observations, including every cell, tensor and internal-path cross term.

use super::super::super::quadrupole::array::{
    ArrayField, FullWaveQuadrupoleSourceLimits, Transfer,
};
use super::super::super::{QuadratureLimits, poll, positive};
use super::super::interior::{PreparedInteriorFullWaveQuadrupole, reduces_to_interface};
use super::super::reflection::{
    PreparedLayeredReflectedFullWaveQuadrupole, reduces_to_first_interface,
};
use super::{
    FullWaveLayerQuadrupoleDiagnostics, PlanarFluidLayer, PreparedLayeredFullWaveQuadrupole,
};
use crate::AeroacError;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use crate::observation::HarmonicFieldLine;
use fs_fft::welch::WelchCrossSpectrum;
use fs_math::{c64::C64, dd::Dd};

mod recorded;

/// The original receiver domain and the direct-path policy are explicit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LayerQuadrupoleField {
    /// Receiver strictly beyond the second face.
    Transmitted,
    /// Receiver before the first face; only boundary returns are retained.
    Reflected,
    /// Receiver before the first face; direct and returned fields are coherent.
    IncidentTotal,
    /// Receiver strictly between both faces, retaining both standing-wave directions.
    Interior,
}

impl LayerQuadrupoleField {
    const fn scope(self) -> &'static str {
        match self {
            Self::Transmitted => {
                "Original prescribed stationary volume-stress cells [N m]; coherent full-wave transmission through an infinite planar lossless nontrapping fluid layer. All tensor components, near/evanescent fields and internal returns remain coherent. Source statistics are explicit. Per-cell absolute estimates give no relative array guarantee at nulls. No source-size/SGS validation, loading, moving media, absorption, finite edges or interval certificate."
            }
            Self::Reflected => {
                "Original prescribed stationary volume-stress cells [N m]; reflected-only full-wave observation before an infinite planar lossless nontrapping fluid layer, including every coherent return and finite reflected self-field. Source statistics are explicit. Per-cell absolute estimates give no relative array guarantee at nulls. No direct field, source-size/SGS validation, loading, moving media, absorption, finite edges or interval certificate."
            }
            Self::IncidentTotal => {
                "Original prescribed stationary volume-stress cells [N m]; coherent direct-plus-reflected full-wave observation before an infinite planar lossless nontrapping fluid layer. Tensor, cell and path interference precede energy statistics. Error estimates control reflected transfer only; no relative total-field guarantee at nulls. No singular total self-field, source-size/SGS validation, loading, moving media, absorption, finite edges or interval certificate."
            }
            Self::Interior => {
                "Original prescribed stationary volume-stress cells [N m] before an infinite planar lossless nontrapping fluid layer; full-wave observation strictly inside it, retaining both standing-wave directions, all tensor components and near/evanescent fields. Source statistics are explicit. Per-cell absolute estimates give no relative array guarantee at nulls. No interior source, source-size/SGS validation, loading, moving media, absorption, finite edges or interval certificate."
            }
        }
    }
}

/// Complete prepared array on its original source clock and geometry.
/// Preparation allocates one bounded transfer table and shares one quadrature
/// cap. Warm fields and recorded factors reuse allocation-free DD reductions.
#[derive(Clone, Debug)]
pub struct PreparedLayeredFullWaveQuadrupoleObserver {
    layer: PlanarFluidLayer,
    field_kind: LayerQuadrupoleField,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<FullWaveLayerQuadrupoleDiagnostics>>,
    evaluations: usize,
}

impl PreparedLayeredFullWaveQuadrupoleObserver {
    /// Every original source lies strictly before the first face. The explicit
    /// field kind admits the receiver's physical fluid, including empty arrays.
    /// Equal speeds alone do not remove density-layer internal returns.
    /// Per-source estimates multiply the physical stress Frobenius norm; their
    /// absolute contributions add without a relative guarantee at array nulls.
    pub fn new(
        layer: PlanarFluidLayer,
        field_kind: LayerQuadrupoleField,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits.quadrupole).admit()?;
        if source_positions_m.len() > limits.max_sources || limits.maximum_total_evaluations == 0 {
            return Err(invalid("layered stress array source or work budget"));
        }
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("layered stress array frequency or observer"));
        }
        let distance = layer.exterior.distance(observer_m);
        let thickness = Dd::from_f64(layer.thickness_m());
        let admitted = match field_kind {
            LayerQuadrupoleField::Transmitted => positive((distance - thickness).to_f64()),
            LayerQuadrupoleField::Reflected | LayerQuadrupoleField::IncidentTotal => {
                positive((-distance).to_f64())
            }
            LayerQuadrupoleField::Interior => {
                positive(distance.to_f64()) && positive((thickness - distance).to_f64())
            }
        };
        if !admitted {
            return Err(invalid(
                "layered stress array receiver is outside its declared fluid domain",
            ));
        }
        let reduced = match field_kind {
            LayerQuadrupoleField::Transmitted => layer.reduces_to_interface(),
            LayerQuadrupoleField::Reflected | LayerQuadrupoleField::IncidentTotal => {
                reduces_to_first_interface(layer)
            }
            LayerQuadrupoleField::Interior => reduces_to_interface(layer),
        };
        let exact = reduced
            && layer.incident_state().sound_speed_m_s()
                == layer.receiving_state().sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("layered stress array allocation"))?;
        let mut evaluations = 0usize;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if field_kind == LayerQuadrupoleField::IncidentTotal && source_m == observer_m {
                return Err(invalid(
                    "layered stress array has a singular total self-field",
                ));
            }
            let remaining = limits.maximum_total_evaluations - evaluations;
            if !exact && remaining == 0 {
                return Err(invalid(
                    "layered stress array exhausted its aggregate work budget",
                ));
            }
            let mut source_limits = limits.quadrupole;
            if !exact {
                source_limits.maximum_evaluations =
                    source_limits.maximum_evaluations.min(remaining);
            }
            let mut field = [[C64::ZERO; 6]; 4];
            let diagnostics = match field_kind {
                LayerQuadrupoleField::Transmitted => {
                    let original = PreparedLayeredFullWaveQuadrupole::new(
                        layer,
                        frequency_hz,
                        source_m,
                        observer_m,
                        source_limits,
                        cancelled,
                    )?;
                    extract_columns(&mut field, |source| original.evaluate_harmonic_line(source))?;
                    *original.diagnostics()
                }
                LayerQuadrupoleField::Interior => {
                    let original = PreparedInteriorFullWaveQuadrupole::new(
                        layer,
                        frequency_hz,
                        source_m,
                        observer_m,
                        source_limits,
                        cancelled,
                    )?;
                    extract_columns(&mut field, |source| original.evaluate_harmonic_line(source))?;
                    *original.diagnostics()
                }
                LayerQuadrupoleField::Reflected | LayerQuadrupoleField::IncidentTotal => {
                    let original = PreparedLayeredReflectedFullWaveQuadrupole::new(
                        layer,
                        frequency_hz,
                        source_m,
                        observer_m,
                        source_limits,
                        cancelled,
                    )?;
                    extract_columns(&mut field, |source| {
                        if field_kind == LayerQuadrupoleField::Reflected {
                            original.evaluate_reflected_harmonic_line(source)
                        } else {
                            original.evaluate_total_harmonic_line(source)
                        }
                    })?;
                    *original.diagnostics()
                }
            };
            evaluations = evaluations
                .checked_add(diagnostics.evaluations)
                .filter(|&work| work <= limits.maximum_total_evaluations)
                .ok_or_else(|| invalid("layered stress array exceeded aggregate work budget"))?;
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            layer,
            field_kind,
            frequency_hz,
            observer_m,
            transfers,
            evaluations,
        })
    }

    pub const fn layer(&self) -> PlanarFluidLayer {
        self.layer
    }
    pub const fn field_kind(&self) -> LayerQuadrupoleField {
        self.field_kind
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
        self.field_kind.scope()
    }
    /// Reflected-transfer diagnostics for both incident-side views; complete
    /// transmitted/interior transfer diagnostics for the other views.
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
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }
    /// Independence applies between original cells; each cell's tensor and
    /// every direct/returned/standing-wave path remain coherent internally.
    pub fn evaluate_independent(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_independent(sources, cancelled)
    }
    /// Explicit peak-moment factors, reduced coherently before their powers add.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[QuadrupolePhasor]],
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view()
            .evaluate_coherence_modes(modes, budget, cancelled)
    }
    /// Native joint six-channel-per-cell factors and one-sided bin normalization.
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

// Six native unit queries copy the already prepared map, without quadrature,
// source differentiation or an alternative tensor/shear convention.
fn extract_columns(
    field: &mut [[C64; 6]; 4],
    mut evaluate: impl FnMut(QuadrupolePhasor) -> Result<HarmonicFieldLine<3>, AeroacError>,
) -> Result<(), AeroacError> {
    for column in 0..6 {
        let mut unit = QuadrupolePhasor::default();
        unit.stress_moment_n_m[column] = C64::ONE;
        let line = evaluate(unit)?;
        field[0][column] = line.pressure_pa;
        for axis in 0..3 {
            field[axis + 1][column] = line.particle_velocity_m_s[axis];
        }
    }
    Ok(())
}
