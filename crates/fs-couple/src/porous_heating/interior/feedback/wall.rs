//! Original material walls on the same accepted thermal/acoustic timeline.

use super::InteriorWallHeatCell;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::wall::material::LocalMaterialWall;
use crate::cavity_modal::interior::heterogeneous::mechanical::{
    CylindricalMaterialBoundary, MaterialMechanicalError, MaterialMechanicalHarmonicError,
};
use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::MechanicalBoundary;
use fs_matdb::MaterialCard;
use std::sync::Arc;

/// One complete original membrane patch and its thermally lumped receiver.
///
/// The boundary retains the thin-cylinder source's geometry, restraint,
/// material axes and applicability. The original card and explicit absolute
/// temperature coordinate select the same property claims at each accepted
/// native temperature. The source reference need not start at the receiver's
/// temperature: cold preparation evaluates it there before the first solve.
///
/// Density, thickness, mass and geometry remain fixed. Source Young-viscosity,
/// source-pair Maxwell and admitted principal-aligned material laws use the
/// existing continuation owner; unsupported source/phase/expansion models
/// refuse. This quasisteady law carries no transient wall memory between
/// harmonic solves and does not infer parametric acoustic stored-energy work.
#[derive(Clone)]
pub struct TemperatureInteriorWall<'t, 'm> {
    pub destination: InteriorWallHeatCell<'t>,
    pub boundary: CylindricalMaterialBoundary<'m>,
    pub card: &'m MaterialCard,
    pub temperature_axis: &'m str,
}

/// Clones share complete immutable source receipts, not just numeric moduli.
#[derive(Clone)]
pub(super) struct State {
    pub(super) temperature_k: f64,
    pub(super) material: Arc<LocalMaterialWall>,
}

/// Evaluate once at the requested temperature, then lower through the same
/// geometric and constitutive owner as MaterialMechanicalHarmonic.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    binding: &TemperatureInteriorWall<'_, '_>,
    index: usize,
    current: Option<&State>,
    temperature_k: f64,
    geometry: &HelmholtzFe<'_>,
    frequency_hz: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(State, MechanicalBoundary), MaterialMechanicalHarmonicError> {
    let source = current.map_or(binding.boundary.source, |state| state.material.as_ref());
    let source_error = |source| MaterialMechanicalError::Source {
        boundary: index,
        source,
    };
    let previous_temperature = current
        .map(|state| state.temperature_k)
        .or_else(|| {
            source
                .material()
                .query_point()
                .iter()
                .find(|(name, _)| name == binding.temperature_axis)
                .map(|(_, value)| *value)
        })
        .ok_or_else(|| {
            source_error(AcousticRealizeError::Nonlinear(
                "material wall source has no declared temperature coordinate".into(),
            ))
        })?;
    let material = source
        .at_temperature_pinned(
            binding.card,
            binding.temperature_axis,
            previous_temperature,
            temperature_k,
            cancelled,
        )
        .map_err(source_error)?;
    let boundary = CylindricalMaterialBoundary {
        faces: binding.boundary.faces.clone(),
        source: &material,
        axis: binding.boundary.axis,
        radius_m: binding.boundary.radius_m,
        relative_geometry_tolerance: binding.boundary.relative_geometry_tolerance,
    };
    let numerical = boundary.prepare_harmonic(index, geometry, frequency_hz, cancelled)?;
    Ok((
        State {
            temperature_k,
            material: Arc::new(material),
        },
        numerical,
    ))
}
