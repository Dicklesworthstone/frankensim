//! Original distributed injection in complete pressure coordinates.
//!
//! This adapter reuses the heterogeneous owner's exact volume/surface
//! projections, coherent source load and signed work observation. Only the
//! receiving pressure operator's rho scaling differs: its unscaled weak
//! equation has i omega rho b on the right-hand side.

use super::*;
use super::super::heterogeneous::source::PressureSourceOwner;
use super::super::heterogeneous::{
    BoundaryVelocitySource, HeterogeneousError, PreparedHeterogeneous, VolumeVelocitySource,
};
use std::collections::BTreeSet;
use std::sync::Arc;

/// Retained only by constructors that copy every original FE coordinate.
pub(super) struct SourceGeometry {
    identity: Arc<()>,
    excluded_faces: BTreeSet<[u32; 3]>,
}

impl SourceGeometry {
    pub(super) fn new(
        geometry: &HelmholtzFe<'_>,
        ports: &[BoundaryMode],
        walls: Option<&[[u32; 3]]>,
    ) -> Self {
        let excluded_faces = ports
            .iter()
            .flat_map(|port| port.traces.iter().map(|&(face, _)| face))
            .chain(walls.into_iter().flatten().copied())
            .map(|mut face| {
                face.sort_unstable();
                face
            })
            .collect();
        Self {
            identity: Arc::clone(&geometry.source_identity),
            excluded_faces,
        }
    }
}

impl InteriorMultiport {
    fn source_owner<'m, 'f, 'c>(
        &'m self,
        geometry: &'f HelmholtzFe<'c>,
    ) -> Result<PressureSourceOwner<'m, 'f, 'c>, HeterogeneousError> {
        let original = self.source_geometry.as_ref().ok_or_else(|| {
            invalid("distributed sources require the complete original pressure coordinates")
        })?;
        if !Arc::ptr_eq(&original.identity, &geometry.source_identity) {
            return Err(invalid("source geometry is not this full pressure operator's owner").into());
        }
        Ok(PressureSourceOwner::Interior { model: self, geometry })
    }

    /// Cache an original signed volume-injection shape [m^-3].
    ///
    /// Supply the same immutable FE owner used by the full-coordinate
    /// constructor, even if either owner has moved since construction.
    /// Geometrically identical substitutions and reduced models refuse.
    /// Amplitudes in subsequent responses are peak volume flow [m^3/s].
    pub fn volume_velocity_source<'m, 'f, 'c, E: From<HeterogeneousError>>(
        &'m self,
        geometry: &'f HelmholtzFe<'c>,
        quadrature_order: usize,
        density_per_m3: &mut impl FnMut(usize, [f64; 3]) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VolumeVelocitySource<'m, 'f, 'c>, E> {
        VolumeVelocitySource::prepare(
            self.source_owner(geometry)?,
            quadrature_order,
            density_per_m3,
            cancelled,
        )
    }

    /// Cache six-node inward boundary motion on the original curved/affine
    /// faces. Values are dimensionless nodal samples, not hierarchical
    /// coefficients; subsequent peak amplitudes are m/s. Signed and zero-net
    /// shapes retain all pressure coordinates without area normalization.
    /// Supports must avoid every port/interface and boundary-layer face.
    pub fn quadratic_velocity_source<'m, 'f, 'c>(
        &'m self,
        geometry: &'f HelmholtzFe<'c>,
        traces: &[([u32; 3], [f64; 6])],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BoundaryVelocitySource<'m, 'f, 'c>, HeterogeneousError> {
        self.quadratic_velocity_source_with_quadrature(
            geometry,
            traces,
            geometry.face_quadrature_order(geometry.order + 2),
            cancelled,
        )
    }

    /// Same source with an explicit independent physical face-quadrature order.
    pub fn quadratic_velocity_source_with_quadrature<'m, 'f, 'c>(
        &'m self,
        geometry: &'f HelmholtzFe<'c>,
        traces: &[([u32; 3], [f64; 6])],
        quadrature_order: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BoundaryVelocitySource<'m, 'f, 'c>, HeterogeneousError> {
        let source = BoundaryVelocitySource::prepare(
            self.source_owner(geometry)?,
            traces,
            quadrature_order,
            cancelled,
        )?;
        let excluded = &self.source_geometry.as_ref().expect("source owner admitted").excluded_faces;
        for face in source.support_faces() {
            if cancelled() {
                return Err(HeterogeneousError::Cancelled);
            }
            if excluded.contains(face) {
                return Err(invalid("boundary source overlaps a port or boundary-layer face").into());
            }
        }
        Ok(source)
    }

    fn source_load(
        &self,
        boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Option<Vec<C64>>, HeterogeneousError> {
        let Some(owner) = boundary.first().map(|(source, _)| source.owner)
            .or_else(|| volume.first().map(|(source, _)| source.owner))
        else {
            return Ok(None);
        };
        match owner {
            PressureSourceOwner::Interior { model, .. } if core::ptr::eq(model, self) => {}
            _ => return Err(invalid("source belongs to another original pressure operator").into()),
        }
        let original = self.source_geometry.as_ref().ok_or_else(|| {
            invalid("distributed sources require the complete original pressure coordinates")
        })?;
        Ok(Some(PreparedHeterogeneous::project_sources_with_boundary_exclusions(
            owner,
            &original.excluded_faces,
            &[],
            boundary,
            volume,
            &mut || cancelled(),
        )?))
    }
}

/// Retain the old port RHS arithmetic; append one coherently combined source.
pub(super) fn add_load(
    rhs: &mut [C64],
    load: Option<&[C64]>,
    omega_rho: f64,
    cancelled: &dyn Fn() -> bool,
) -> Result<(), AcousticFeError> {
    if let Some(load) = load {
        for (entry, &source) in rhs.iter_mut().zip(load) {
            poll(cancelled)?;
            *entry = *entry + C64::new(-source.im, source.re).scale(omega_rho);
            if !finite(*entry) {
                return Err(invalid("interior distributed source load overflow"));
            }
        }
    }
    Ok(())
}

fn record_work(
    mut field: InteriorPortResponse,
    boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
    volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
    cancelled: &dyn Fn() -> bool,
) -> Result<InteriorPortResponse, HeterogeneousError> {
    if !boundary.is_empty() || !volume.is_empty() {
        let pressure = field.amplitudes.iter().map(|p| C64::new(p[0], p[1])).collect::<Vec<_>>();
        let (boundary_work, volume_work) = PreparedHeterogeneous::source_power(
            &pressure, boundary, volume, &mut || cancelled(),
        )?;
        field.boundary_source_power_w = boundary_work;
        field.volume_source_power_w = volume_work;
    }
    Ok(field)
}

impl PreparedInteriorMultiport<'_> {
    /// Coherent original modal flows, quadratic boundary motion [m/s] and
    /// volume injection [m^3/s], using the retained unloaded factor.
    pub fn response_with_sources(
        &self,
        inward_flow_m3_s: &[C64],
        boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<InteriorPortResponse, HeterogeneousError> {
        let load = self.model.source_load(boundary, volume, cancelled)?;
        let field = self.response_with_load(inward_flow_m3_s, load.as_deref(), cancelled)?;
        record_work(field, boundary, volume, cancelled)
    }
}

impl PreparedInteriorScattering<'_> {
    /// Coherent incident power waves and original pressure-region sources.
    /// The same real-reference loaded factor remains usable at an unloaded
    /// cavity pole. Signed source work is separate from actual port work,
    /// boundary-layer loss and molecular volume dissipation.
    pub fn response_with_sources(
        &self,
        incident: &[C64],
        boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<InteriorPortResponse, HeterogeneousError> {
        let load = self.model.source_load(boundary, volume, cancelled)?;
        let field = self.response_with_load(incident, load.as_deref(), cancelled)?;
        record_work(field, boundary, volume, cancelled)
    }

    /// Original affine outgoing source waves for reciprocal network closure.
    /// No per-source solve or passive-scattering modification is introduced.
    pub fn source_outgoing(
        &self,
        boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        if boundary.is_empty() && volume.is_empty() {
            poll(cancelled)?;
            return Ok(vec![C64::ZERO; self.reference.len()]);
        }
        let load = self.model.source_load(boundary, volume, cancelled)?;
        let field = self.response_with_load(
            &vec![C64::ZERO; self.reference.len()], load.as_deref(), cancelled,
        )?;
        Ok(field.pressure_pa.iter().zip(&self.roots)
            .map(|(p, &root)| p.scale(1.0 / root)).collect())
    }
}
