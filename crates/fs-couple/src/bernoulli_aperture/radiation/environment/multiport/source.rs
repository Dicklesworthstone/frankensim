//! One source selection for the existing multiport reaction and observation.
//!
//! Imported operators keep their independent source gas. Computed operators
//! instead sample the actual declared apertures in the receiver's resolved gas.
//! Neither route modifies the live mechanics or adds another terminal runtime.

use super::{MimoRadiationLimits, MimoRadiationModel};
use super::bake::{BaffledMimoBakeOptions, BaffledMimoBakeReport, bake_baffled_mimo};
use super::bake::mixed::{DrivenMouthSurface, bake_driven_mouth, bake_driven_mouths};
use super::input::{MimoRadiationInput, MimoRadiationInputError, MimoRadiationRecord};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::baffled::openings::duct::{DuctPressure, DuctRadiationSpec, DuctView};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources};
use fs_vfit::causality::CausalityOptions;

fn bad(what: &'static str) -> MimoRadiationInputError {
    MimoRadiationInputError::Input { line: 0, what }
}

/// Explicit cold-build controls. The upper physical band is always the actual
/// mechanical Nyquist frequency; a narrower microphone band cannot truncate it.
#[derive(Clone, Copy, Debug)]
pub struct ComputedBaffledInput {
    minimum_frequency_hz: f64,
    options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits,
}
impl ComputedBaffledInput {
    /// Number of exact whitespace-delimited native fields after the selector.
    pub const FIELDS: usize = 14;

    /// Read `FMIN_HZ SAMPLES ORDER1 ORDER2 ORDER3 ITERATIONS MAX_REL_ERROR
    /// FIT_MEMORY_BYTES SUBTRACTIONS TAIL_BOUND_RATIO DATA_TOLERANCE MAX_STATES
    /// RUNTIME_MEMORY_BYTES PASSIVITY_POINTS`. No file or material name implies
    /// a numerical choice. The original bake owner validates its full domain.
    pub fn parse(text: &str) -> Result<Self, MimoRadiationInputError> {
        let fields: Vec<_> = text.split_ascii_whitespace().take(Self::FIELDS + 1).collect();
        if fields.len() != Self::FIELDS {
            return Err(bad("computed-baffled radiation requires exactly fourteen numerical fields"));
        }
        let scalar = |index: usize| -> Result<f64, MimoRadiationInputError> {
            let value: f64 = fields[index].parse().map_err(|_| bad("invalid computed-baffled scalar"))?;
            if value.is_finite() { Ok(value) } else { Err(bad("computed-baffled scalars must be finite")) }
        };
        let count = |index: usize, minimum: usize, maximum: usize| -> Result<usize, MimoRadiationInputError> {
            let value = fields[index].parse().map_err(|_| bad("invalid computed-baffled unsigned count"))?;
            if (minimum..=maximum).contains(&value) { Ok(value) }
            else { Err(bad("computed-baffled count exceeds the bounded source domain")) }
        };
        let minimum_frequency_hz = scalar(0)?;
        let options = BaffledMimoBakeOptions {
            frequency_samples: count(1, 32, 513)?,
            orders: [count(2, 1, 64)?, count(3, 1, 64)?, count(4, 1, 64)?],
            fit_iterations: count(5, 1, 32)?,
            maximum_relative_error: scalar(6)?,
            max_fit_memory_bytes: count(7, 1, 512 * 1024 * 1024)?,
            causality: CausalityOptions {
                subtractions: count(8, 1, 8)?,
                tail_bound_ratio: scalar(9)?,
                relative_data_tolerance: scalar(10)?,
            },
        };
        let limits = MimoRadiationLimits {
            max_states: count(11, 1, 256)?,
            max_memory_bytes: count(12, 1, 64 * 1024 * 1024)?,
            sampled_passivity_points: count(13, 17, 65_537)?,
        };
        if minimum_frequency_hz <= 0.0
            || !options.orders.windows(2).all(|pair| pair[0] < pair[1])
            || !(0.0 < options.maximum_relative_error && options.maximum_relative_error <= 1.0)
            || options.causality.tail_bound_ratio < 1.0
            || !(0.0..=1.0).contains(&options.causality.relative_data_tolerance)
        {
            return Err(bad("computed-baffled radiation needs a positive lower band, three increasing orders, valid error bounds and causal-tail assumptions"));
        }
        Ok(Self { minimum_frequency_hz, options, limits })
    }
    /// Original positive lower band edge in Hz.
    pub const fn minimum_frequency_hz(&self) -> f64 { self.minimum_frequency_hz }
    /// Original numerical refinement and source-data assumptions.
    pub const fn options(&self) -> BaffledMimoBakeOptions { self.options }
    /// Unchanged terminal runtime limits, shared by imported and computed sources.
    pub const fn limits(&self) -> MimoRadiationLimits { self.limits }

    fn build(
        &self,
        graph: &TubeNetworkSpec,
        gases: &[GasState],
        radiation: &DuctRadiationSpec,
        exterior: &GasRecord,
        gate: &CancelGate,
    ) -> Result<(MimoRadiationModel, ComputedBaffledRecord), MimoRadiationInputError> {
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        // Resolve the observer before expensive sampling, using its existing
        // geometric, gas, physical-clock and outlet-order admissions.
        let observer = DuctPressure::new(&[DuctView { graph, gases }], radiation.clone())
            .map_err(|source| MimoRadiationInputError::Physics(AcousticRealizeError::Nonlinear(source)))?;
        if observer.exterior_gas() != exterior.state() {
            return Err(bad("computed radiation must retain the actual receiver's original gas record"));
        }
        let mut apertures = Vec::with_capacity(radiation.openings.len());
        for opening in &radiation.openings {
            if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
            if opening.network != 0 || !matches!(graph.nodes.get(opening.node), Some(NetworkNode::Termination { .. })) {
                return Err(bad("computed radiation requires explicitly observed memoryless physical endpoints"));
            }
            let mut adjacent = graph.sections.iter().filter(|section| section.nodes.contains(&opening.node));
            let section = adjacent.next().ok_or_else(|| bad("computed radiation aperture has no physical duct section"))?;
            if adjacent.next().is_some() {
                return Err(bad("computed radiation aperture must have exactly one physical duct section"));
            }
            apertures.push(CircularOpening { center_m: opening.center_m, radius_m: section.radius_m });
        }
        self.build_apertures(exterior, &apertures,
            radiation.reconstruction.mechanics_rate_hz(), gate)
    }

    /// Build the same physical baffled operator from original circular openings,
    /// independent of the propagation owner. Geometry and exterior gas are
    /// explicit; the original BEM/VF bake retains its full-band accuracy and
    /// passivity admissions. No waveguide surrogate or observer is constructed.
    pub fn build_apertures(
        &self,
        exterior: &GasRecord,
        apertures: &[CircularOpening],
        mechanics_rate_hz: u32,
        gate: &CancelGate,
    ) -> Result<(MimoRadiationModel, ComputedBaffledRecord), MimoRadiationInputError> {
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        let band_hz = [self.minimum_frequency_hz, 0.5 * f64::from(mechanics_rate_hz)];
        let (model, report) = bake_baffled_mimo(*exterior.state(), apertures, mechanics_rate_hz,
            band_hz, self.options, self.limits, gate).map_err(MimoRadiationInputError::Bake)?;
        let record = ComputedBaffledRecord {
            declaration_identity: model.declaration_identity,
            realization_identity: model.realization_identity,
            gas_record: exterior.clone(), apertures: apertures.to_vec(), band_hz,
            options: self.options, limits: self.limits, report, driven_mouths: Vec::new(),
        };
        Ok((model, record))
    }
}

/// Original physical source and measured cold-build results, with no duplicate
/// copy of the rational states or live instrument mechanics.
#[derive(Clone, Debug)]
pub struct ComputedBaffledRecord {
    declaration_identity: [u8; 32],
    realization_identity: [u8; 32],
    gas_record: GasRecord,
    apertures: Vec<CircularOpening>,
    band_hz: [f64; 2],
    options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits,
    report: BaffledMimoBakeReport,
    driven_mouths: Vec<DrivenMouthSurface>,
}
impl ComputedBaffledRecord {
    /// Identity of actual source geometry, atmosphere and numerical declaration.
    pub const fn declaration_identity(&self) -> [u8; 32] { self.declaration_identity }
    /// Identity of the accepted original fitted impedance realization.
    pub const fn realization_identity(&self) -> [u8; 32] { self.realization_identity }
    /// Original resolved exterior, including thermochemical source receipts.
    pub const fn gas_record(&self) -> &GasRecord { &self.gas_record }
    /// Actual circular terminals in their original observer/terminal order.
    /// Driven mouths lead the joint matrix in their declared order, followed
    /// by these same circular channels without reordering.
    pub fn apertures(&self) -> &[CircularOpening] { &self.apertures }
    /// Original singleton view; multiple independent mouths have no one-mouth
    /// representation. Use driven_mouths for the complete leading channel group.
    pub const fn driven_mouth(&self) -> Option<DrivenMouthSurface> {
        match self.driven_mouths.as_slice() { [mouth] => Some(*mouth), _ => None }
    }
    pub fn driven_mouths(&self) -> &[DrivenMouthSurface] { &self.driven_mouths }

    /// Replace this cold circular image with one joint rectangle/terminal fit.
    /// The original full gas record, complete mechanical band and numerical
    /// allowances are retained. Identities and reports describe the new joint
    /// model; no independent old load remains in parallel or masquerades as its
    /// source. The returned model has mouth-first pressure/volume-flow channels.
    pub fn with_driven_mouth(
        &self, mouth: DrivenMouthSurface, mechanics_rate_hz: u32, gate: &CancelGate,
    ) -> Result<(fs_vfit::mimo::MimoModel, Self), MimoRadiationInputError> {
        self.with_driven_mouths(&[mouth], mechanics_rate_hz, gate)
    }

    /// One original reciprocal image of every supplied rectangle and circular
    /// outlet. Leading channels follow the slice order. Self and mutual terms
    /// share one physical source, fit and causal radiation state.
    pub fn with_driven_mouths(
        &self, mouths: &[DrivenMouthSurface], mechanics_rate_hz: u32, gate: &CancelGate,
    ) -> Result<(fs_vfit::mimo::MimoModel, Self), MimoRadiationInputError> {
        if !self.driven_mouths.is_empty() || mouths.is_empty()
            || self.band_hz[1].to_bits() != (0.5 * f64::from(mechanics_rate_hz)).to_bits()
        {
            return Err(bad("joint mouth replacement requires the original clock and one unmodified circular source"));
        }
        let baked = if let [mouth] = mouths {
            bake_driven_mouth(self.gas_record.state(), *mouth, &self.apertures,
                mechanics_rate_hz, self.band_hz, self.options, self.limits, gate)
        } else {
            bake_driven_mouths(self.gas_record.state(), mouths, &self.apertures,
                mechanics_rate_hz, self.band_hz, self.options, self.limits, gate)
        };
        let (model, report, declaration_identity, realization_identity) =
            baked.map_err(MimoRadiationInputError::Bake)?;
        Ok((model, Self { declaration_identity, realization_identity,
            gas_record: self.gas_record.clone(), apertures: self.apertures.clone(),
            band_hz: self.band_hz, options: self.options, limits: self.limits,
            report, driven_mouths: mouths.to_vec() }))
    }
    /// Full admitted physical band, ending at mechanical Nyquist.
    pub const fn band_hz(&self) -> [f64; 2] { self.band_hz }
    /// Original fitting and source-data assumptions.
    pub const fn options(&self) -> BaffledMimoBakeOptions { self.options }
    /// Runtime storage and sampled-passivity limits.
    pub const fn limits(&self) -> MimoRadiationLimits { self.limits }
    /// Accuracy and passivity outcomes from the original cold builder.
    pub const fn report(&self) -> &BaffledMimoBakeReport { &self.report }
}

/// Either original imported coefficients or a computed operator of the actual
/// baffled geometry and material gas. Both attach through one terminal owner.
#[derive(Clone, Debug)]
pub enum MimoSourceInput {
    /// Independently authored operator with its own source declaration.
    Artifact(MimoRadiationInput),
    /// Actual physical opening array, recomputed at source construction.
    ComputedBaffled(ComputedBaffledInput),
}
/// Source-specific records; neither variant owns an additional filter.
#[derive(Clone, Debug)]
pub enum MimoSourceRecord {
    /// Exact imported-byte identity and independent source receipts.
    Artifact(MimoRadiationRecord),
    /// Actual geometry/gas and computed approximation diagnostics.
    ComputedBaffled(ComputedBaffledRecord),
}
impl MimoSourceInput {
    /// Runtime limits shared by either source route.
    pub const fn limits(&self) -> MimoRadiationLimits {
        match self { Self::Artifact(input) => input.limits(), Self::ComputedBaffled(input) => input.limits() }
    }
    /// Build the source once against the original graph and observer context.
    /// Only imported selection calls the external artifact loader. The actual
    /// attach operation still belongs to the unchanged instrument/runtime owner.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        &self,
        graph: &TubeNetworkSpec,
        gases: &[GasState],
        radiation: &DuctRadiationSpec,
        exterior: &GasRecord,
        gas_sources: Option<&ThermochemicalSources<'_>>,
        loader: &mut impl FnMut(&str) -> Result<Vec<u8>, String>,
        gate: &CancelGate,
    ) -> Result<(MimoRadiationModel, MimoSourceRecord), MimoRadiationInputError> {
        match self {
            Self::Artifact(input) => input.load(gas_sources, loader, gate)
                .map(|(model, record)| (model, MimoSourceRecord::Artifact(record))),
            Self::ComputedBaffled(input) => input.build(graph, gases, radiation, exterior, gate)
                .map(|(model, record)| (model, MimoSourceRecord::ComputedBaffled(record))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "1 129 8 12 16 12 0.01 67108864 2 2 0.0001 64 1048576 257";

    #[test]
    fn numerical_source_choices_are_explicit_bounded_and_retained() {
        let source = ComputedBaffledInput::parse(INPUT).unwrap();
        assert_eq!(source.minimum_frequency_hz(), 1.0);
        assert_eq!(source.options().orders, [8, 12, 16]);
        assert_eq!(source.options().causality.subtractions, 2);
        assert_eq!(source.limits().max_states, 64);
        for text in [
            INPUT.replacen("1 ", "0 ", 1), INPUT.replace("129", "514"),
            INPUT.replace("8 12 16", "8 8 16"), INPUT.replace("0.01", "NaN"),
            INPUT.replace("67108864", "536870913"), INPUT.replace("2 2 0.0001", "2 0.5 0.0001"),
            INPUT.replace("1048576", "0"), INPUT.replace("257", "16"),
            format!("{INPUT} extra"), INPUT.rsplit_once(' ').unwrap().0.into(),
        ] { assert!(ComputedBaffledInput::parse(&text).is_err(), "admitted {text}"); }
    }
}
