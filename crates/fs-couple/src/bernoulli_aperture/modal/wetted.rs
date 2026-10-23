//! Actual fluid-facing geometry in the ORIGINAL modal acoustic coupling.
//!
//! Each selected midsurface triangle names which physical side faces upstream;
//! the opposite side faces downstream. The shared shell pressure owner computes
//! A_i = integral(phi_i dot n_upstream dA) on those original facets. The existing
//! valve then uses -A_i dp for force AND -sum(A_i v_i) for swept volume. Swapping
//! the upstream physical face reverses pressure forcing relative to gap motion
//! without supplying a signed effective area or changing the material pencil.
//!
//! This is a fixed P1 normal trace of the retained flat DKT displacement. It
//! does not infer changing wetting, channel pressure, finite-thickness moments,
//! curved-shell motion or a measured operating threshold. Excluded triangles
//! retain their mechanical/contact participation but do not border this paired
//! pressure port. The independent prescribed-pressure ports remain separate.

use super::{AcousticRealizeError, ModalAperture, invalid, nonlinear};
use super::material::MaterialModalPlate;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_plate::shell::{ShellMesh, pressure::{OrientedShellFacet, ShellPressurePort,
    projection_terms}, profile::ProfileBudget};

fn geometry_error(error: fs_plate::PlateError) -> AcousticRealizeError {
    match error {
        fs_plate::PlateError::Cancelled => AcousticRealizeError::Cancelled,
        other => nonlinear(other),
    }
}

/// Source-bound original triangle sides facing the upstream fluid.
/// A facet's opposite physical side faces the downstream fluid. Absence means
/// no participation in this pressure-difference port, not removal of material.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaterialWettedSurface {
    specimen: ContentHash,
    upstream_facets: Vec<OrientedShellFacet>,
    provenance: String,
}
impl MaterialWettedSurface {
    /// Bind original triangle addresses and winding-side geometry to a specimen.
    /// No signed pressure area or modal coefficient is accepted from the caller.
    /// The receiving valve checks identity, uniqueness, bounds and projection.
    pub fn from_plate(plate: &MaterialModalPlate, upstream_facets: Vec<OrientedShellFacet>,
        provenance: String) -> Self
    {
        Self { specimen: plate.source().specimen_identity(), upstream_facets, provenance }
    }
    /// Original geometry/material identity; a changed specimen needs a new binding.
    pub const fn specimen_identity(&self) -> ContentHash { self.specimen }
    /// Actual upstream-facing sides of the original triangles, in source order.
    pub fn upstream_facets(&self) -> &[OrientedShellFacet] { &self.upstream_facets }
    /// Independent geometry source or explicit manufactured-model declaration.
    pub fn provenance(&self) -> &str { &self.provenance }
}

/// One shared physical projection for acoustic force and displaced volume.
#[derive(Clone, Debug)]
pub struct MaterialWettedPort {
    surface: MaterialWettedSurface,
    area_m2: f64,
    modal_area_m2: Vec<f64>,
}

/// Cold description data used by the material constructor before any acoustic
/// clock admission. Source identity is bound to the actual reduced blade there.
pub(crate) struct MaterialWettedInput {
    pub blade: usize,
    pub facets: Vec<OrientedShellFacet>,
    pub provenance: String,
    pub max_projection_terms: usize,
}
impl MaterialWettedInput {
    pub(super) fn bind(self, plates: &mut [MaterialModalPlate], gate: &CancelGate)
        -> Result<(), AcousticRealizeError>
    {
        let plate = plates.get_mut(self.blade)
            .ok_or_else(|| invalid("wetted surface names an absent original material blade"))?;
        let surface = MaterialWettedSurface::from_plate(plate, self.facets, self.provenance);
        let port = MaterialWettedPort::project(plate, surface, self.max_projection_terms, gate)?;
        plate.bind_wetted_surface(port);
        Ok(())
    }
}

impl MaterialWettedPort {
    /// Exact source geometry and its material-specimen binding.
    pub const fn surface(&self) -> &MaterialWettedSurface { &self.surface }
    /// Selected physical midsurface area [m²], independent of modal cancellation.
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Signed overlaps [m²] in this blade's original local mode order.
    /// These same values appear in ModalApertureSpec at the blade's mode_range.
    pub fn modal_area_m2(&self) -> &[f64] { &self.modal_area_m2 }

    fn project(plate: &MaterialModalPlate, surface: MaterialWettedSurface,
        max_projection_terms: usize, gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        let poll = || if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) };
        poll()?;
        if plate.wetted_surface().is_some() || surface.specimen != plate.source().specimen_identity()
            || surface.provenance.trim().is_empty() || surface.provenance.len() > 4096 {
            return Err(invalid("wetted geometry requires one cold binding to the original specimen and a source"));
        }
        let original = &plate.source().chart().mesh;
        let n = plate.mode_range().len();
        if projection_terms(original.nodes.len(), surface.upstream_facets.len(), n)
            .is_none_or(|terms| terms > max_projection_terms) {
            return Err(invalid("wetted geometry projection exceeds its declared work allowance"));
        }
        // Embed the ORIGINAL planar chart without fitting or remeshing it. The
        // shared oriented-facet owner supplies normal signs and area quadrature.
        let mesh = ShellMesh::new(original.nodes.iter().map(|&(x,y)| [x,y,0.0]).collect(),
            original.tris.clone()).map_err(geometry_error)?;
        let geometry = ShellPressurePort::new(&mesh, surface.upstream_facets.clone(),
            ProfileBudget { max_nodes: plate.options().max_nodes,
                max_triangles: plate.options().max_triangles, max_feature_evaluations: 0 },
            || gate.is_requested()).map_err(geometry_error)?;
        let mut traces = vec![vec![0.0; 6 * original.nodes.len()]; n];
        for (node, shape) in plate.displacement_shapes().iter().enumerate() {
            poll()?;
            for (mode, trace) in traces.iter_mut().enumerate() {
                trace[6 * node + 2] = shape[mode];
            }
        }
        let area = geometry.project_modes(&traces, max_projection_terms,
            || gate.is_requested()).map_err(geometry_error)?;
        poll()?;
        Ok(Self { surface, area_m2: geometry.area_m2(), modal_area_m2: area })
    }
}

impl ModalAperture {
    /// Replace one original blade's whole-face acoustic projection with explicit
    /// wetted geometry before stepping. Pressure traction and swept volume change
    /// together in the existing same-sample nonlinear junction. Mass, stiffness,
    /// force ports, gap/contact geometry and every material history are retained.
    ///
    /// Full-C midpoint and exact/Maxwell pressure responses are rebound through
    /// their already factored mechanical operators, without another eigensolve
    /// or resetting the initial internal strains. An explicit empty selection
    /// shelters this blade from the coupled pressure while external forces and
    /// contact can still move it and change the physical flow opening.
    ///
    /// Omitting this attachment preserves the original arithmetic. Only flat
    /// original material plates are admitted; this is not a shell-wall adapter.
    pub fn with_material_wetted_surface(mut self, blade: usize, surface: MaterialWettedSurface,
        max_projection_terms: usize, gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        if self.samples != 0 {
            return Err(invalid("wetted geometry must be attached before the first accepted sample"));
        }
        let plate = self.material_plates.get(blade)
            .ok_or_else(|| invalid("wetted surface names an absent original material blade"))?;
        let port = MaterialWettedPort::project(plate, surface, max_projection_terms, gate)?;
        let range = plate.mode_range();
        for (mode, projection) in self.spec.modes[range].iter_mut().zip(port.modal_area_m2()) {
            mode.pressure_area_m2 = *projection;
        }
        if let Some(viscosity) = &mut self.viscosity {
            viscosity.rebind_pressure(&self.spec)?;
        }
        if let Some(exact) = &mut self.exact {
            exact.rebind_pressure(&self.spec)?;
        }
        self.admit_linear_clock()?;
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        self.material_plates[blade].bind_wetted_surface(port);
        Ok(self)
    }
}
