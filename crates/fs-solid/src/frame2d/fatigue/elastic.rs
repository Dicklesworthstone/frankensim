//! Signed axial/bending stress-life on original material-defined frame layers.
//!
//! Every gauge is a physical point in one named annulus at an existing Lobatto
//! station. The common section's original layer modulus acts on the accepted
//! native `epsilon_0 - y*kappa`; no effective homogeneous fatigue material is
//! constructed. This is a local normal-stress estimate, not a shear/torsion or
//! multiaxial criterion. The stress-free reference and source temperature stay
//! fixed. Members with material viscosity refuse because their total local
//! stress would also require the separately resolved viscous contribution.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fatigue::strain_life::StrainStressSample;
use fs_material::fatigue::{
    FatigueScope, FatigueSummary, MeanStressCorrection, RainflowLimits, RainflowMiner, SnCurve,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{Frame2d, FrameFatigueError, SectionModel, gauge_error, poll};

/// A local signed normal-stress history in an ORIGINAL concentric material layer.
#[derive(Debug, Clone)]
pub struct FrameElasticFatigueGauge {
    pub element: usize,
    pub section: usize,
    pub layer: usize,
    /// `[y,z]` in metres in the member cross-section. `y` is the native
    /// in-plane bending arm. The radius must lie in the explicitly named layer;
    /// a material interface may be observed separately from either side.
    pub cross_section_position_m: [f64; 2],
    /// Exact retained absolute-temperature coordinate; no alias is guessed.
    pub temperature_axis: String,
    pub curve: Arc<SnCurve>,
    pub mean_stress: MeanStressCorrection,
}

#[derive(Debug, Clone)]
pub(super) struct ElasticGaugeState {
    binding: FrameElasticFatigueGauge,
    point: ElasticPoint,
    index: usize,
    counter: RainflowMiner,
    sample: StrainStressSample,
}

/// Shared original-layer admission and axial/bending stress recovery. Tensor
/// observers add their independently owned shear field to this same normal law.
#[derive(Debug, Clone)]
pub(super) struct ElasticPoint {
    pub(super) material: Arc<ResolvedMaterialStatePoint>,
    element: usize,
    section: usize,
    coordinates: [f64; 2],
    outer_radius: f64,
    youngs_modulus_pa: f64,
    maximum_strain: f64,
}

impl ElasticPoint {
    pub(super) fn new(
        frame: &Frame2d,
        indices: [usize; 3],
        coordinates: [f64; 2],
        maximum_strain: f64,
    ) -> Result<Self, FrameFatigueError> {
        let [element_index, station, layer_index] = indices;
        let invalid = FrameFatigueError::Invalid;
        let element = frame
            .elements
            .get(element_index)
            .ok_or_else(|| invalid("elastic fatigue member is outside the original frame"))?;
        if element.viscous_section.is_some() {
            return Err(invalid(
                "elastic-only fatigue cannot omit a member's material viscous stress",
            ));
        }
        let section = element.material_section.as_ref().ok_or_else(|| {
            invalid("elastic fatigue requires an original material-defined member")
        })?;
        if !matches!(
            element.sections.get(station),
            Some(SectionModel::Elastic { .. })
        ) {
            return Err(invalid(
                "elastic fatigue station must be an original elastic Lobatto section",
            ));
        }
        let layer = section
            .layers()
            .get(layer_index)
            .ok_or_else(|| invalid("fatigue layer is outside the original concentric section"))?;
        let [y, z] = coordinates;
        let radius = y.hypot(z);
        let inner = if layer_index == 0 {
            section.inner_radius_m()
        } else {
            section.layers()[layer_index - 1].outer_radius_m
        };
        if !y.is_finite()
            || !z.is_finite()
            || !radius.is_finite()
            || radius < inner
            || radius > layer.outer_radius_m
        {
            return Err(invalid(
                "fatigue point must lie within its named material annulus",
            ));
        }
        Ok(Self {
            material: Arc::new(layer.material.clone()),
            element: element_index,
            section: station,
            coordinates,
            outer_radius: layer.outer_radius_m,
            youngs_modulus_pa: section.layer_integrals()[layer_index].young_modulus_pa(),
            maximum_strain,
        })
    }

    pub(super) fn sample(
        &self,
        frame: &Frame2d,
        whole_annulus: bool,
    ) -> Result<StrainStressSample, FrameFatigueError> {
        // The path exposes only immutable mechanics; original source bindings
        // and these accepted Lobatto station indices cannot change mid-record.
        let (epsilon_0, curvature) = frame.elements[self.element].commit.d[self.section];
        let strain = epsilon_0 - self.coordinates[0] * curvature;
        let bound = if whole_annulus {
            epsilon_0.abs() + self.outer_radius * curvature.abs()
        } else {
            strain.abs()
        };
        if !strain.is_finite() || !bound.is_finite() || bound > self.maximum_strain {
            return Err(FrameFatigueError::Invalid(
                "elastic fatigue strain exceeds the declared small-strain domain",
            ));
        }
        let stress = self.youngs_modulus_pa * strain;
        if !stress.is_finite() {
            return Err(FrameFatigueError::Invalid(
                "elastic fatigue point stress is unrepresentable",
            ));
        }
        Ok(StrainStressSample {
            mechanical_strain: strain,
            normal_stress_pa: stress,
        })
    }
}

/// `scope_material` is the name the calibration's source gives this layer's
/// material: an S-N source writes "chemistry; phase" (fs-material
/// `fatigue::source`, as fs-plate's S-N gauge reads it); a Findley source writes
/// the chemistry.
pub(super) fn admit_scope(
    material: &ResolvedMaterialStatePoint,
    temperature_axis: &str,
    scope: &FatigueScope,
    scope_material: &str,
    index: usize,
) -> Result<(), FrameFatigueError> {
    let mismatch = |what| FrameFatigueError::SourceMismatch { gauge: index, what };
    let absolute_temperature = material
        .axis_quantities()
        .get(temperature_axis)
        .is_some_and(|quantity| {
            *quantity == QuantitySpec::dimensional(Temperature::DIMS)
                || quantity.semantic_type().is_some_and(|semantic| {
                    semantic.kind() == QuantityKind::AbsoluteTemperature
                        && semantic.form() == ValueForm::Static
                })
        });
    if temperature_axis.trim().is_empty()
        || material.state_coordinate(temperature_axis) != Some(scope.temperature_k)
        || !absolute_temperature
    {
        return Err(mismatch(
            "fatigue needs the layer's actual dimensioned absolute-temperature coordinate",
        ));
    }
    if scope.material != scope_material || scope.process != material.material().process {
        return Err(mismatch(
            "fatigue material/process differs from the original layer",
        ));
    }
    Ok(())
}

impl ElasticGaugeState {
    pub(super) fn new(
        frame: &Frame2d,
        gauge: FrameElasticFatigueGauge,
        index: usize,
        limits: RainflowLimits,
        maximum_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, FrameFatigueError> {
        poll(cx)?;
        let mismatch = |what| FrameFatigueError::SourceMismatch { gauge: index, what };
        let point = ElasticPoint::new(
            frame,
            [gauge.element, gauge.section, gauge.layer],
            gauge.cross_section_position_m,
            maximum_strain,
        )?;
        let material = point.material.as_ref();
        let named = material.material();
        admit_scope(
            material,
            &gauge.temperature_axis,
            gauge.curve.scope(),
            &format!("{}; {}", named.chemistry, named.phase),
            index,
        )?;
        if let Some(source) = gauge.curve.source() {
            if gauge.temperature_axis != source.request().temperature_axis
                || !source.matches_material_state(
                    material.card_identity(),
                    material.material(),
                    &material.query(),
                )
            {
                return Err(mismatch(
                    "S-N source card, state or operating conditions differ from the native layer",
                ));
            }
        }
        let counter = RainflowMiner::new(gauge.curve.clone(), gauge.mean_stress, limits)
            .map_err(|error| gauge_error(index, error))?;
        let mut state = Self {
            point,
            binding: gauge,
            index,
            counter,
            sample: StrainStressSample {
                mechanical_strain: 0.0,
                normal_stress_pa: 0.0,
            },
        };
        state.observe(frame, cx)?;
        Ok(state)
    }

    pub(super) fn observe(
        &mut self,
        frame: &Frame2d,
        cx: &Cx<'_>,
    ) -> Result<(), FrameFatigueError> {
        poll(cx)?;
        let sample = self.point.sample(frame, false)?;
        self.counter
            .push_with_cancel(sample.normal_stress_pa, &mut || cx.checkpoint().is_err())
            .map_err(|error| gauge_error(self.index, error))?;
        self.sample = sample;
        Ok(())
    }

    pub(super) fn sample(&self) -> StrainStressSample {
        self.sample
    }
    pub(super) fn summary(&self) -> FatigueSummary {
        self.counter.summary()
    }
    pub(super) fn curve(&self) -> &SnCurve {
        &self.binding.curve
    }
    pub(super) fn material_source(&self) -> &ResolvedMaterialStatePoint {
        &self.point.material
    }

    pub(super) fn preview_finish(&self, cx: &Cx<'_>) -> Result<FatigueSummary, FrameFatigueError> {
        self.counter
            .preview_finish_with_cancel(&mut || cx.checkpoint().is_err())
            .map_err(|error| gauge_error(self.index, error))
    }

    pub(super) fn finish(&mut self, cx: &Cx<'_>) -> Result<FatigueSummary, FrameFatigueError> {
        self.counter
            .finish_with_cancel(&mut || cx.checkpoint().is_err())
            .map_err(|error| gauge_error(self.index, error))
    }
}
