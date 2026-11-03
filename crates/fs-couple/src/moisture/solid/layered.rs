//! Heterogeneous moisture layers bonded into one compatible native solid.
//!
//! Water-transfer contact resistance and mechanical bonding are distinct laws:
//! the original moisture stack owns water flux; this consumer shares interface
//! mesh nodes and therefore assumes perfect mechanical bonding. Each original
//! layer retains its own dry skeleton, source datum, elasticity and swelling
//! directions. Ordinary steps retain reference-volume transport. The explicit
//! `advance_spatial_diffusion` path couples the axial spatial transport metric
//! back to the actual compatible deformation while conserving the dry skeleton.

use super::*;
#[path = "layered/deforming.rs"]
mod deforming;
use crate::moisture::{
    adaptive::{LayeredMoistureAdvanceReport, MoistureAdvanceError, MoistureTimeControls},
    atmosphere::{
        AtmosphereMoistureError, AtmosphereMoistureSolve, AtmosphereMoistureStep,
        FiniteMoistureAtmosphere,
    },
    layered::{
        MoistureLayeredFace, MoistureLayeredReport, MoistureLayeredSlab,
        atmosphere::LayeredAtmosphereMoistureCandidate,
    },
};
pub use deforming::MoistureSpatialDiffusionReport;

/// Mechanical source choices for exactly one layer in original stack order.
#[derive(Clone, Debug)]
pub struct MoistureLayerSolidSource {
    pub stiffness: SolidStiffnessModel,
    pub elastic_selection: MaterialPropertySelection,
    pub expansion_selection: MaterialPropertySelection,
    pub axes: HygroscopicAxes,
}

/// One conforming solid with distinct original material/datum bindings.
#[derive(Clone, Debug)]
pub struct MoistureLayeredSolid {
    definition: ContentHash,
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    sources: Vec<MoistureSolidSource>,
    /// Physical cell ranges in global source/expansion output order.
    layer_cells: Vec<std::ops::Range<usize>>,
    reference_cell_volumes_m3: Vec<f64>,
    fixed_dofs: Vec<usize>,
    mesh_frame: ContentHash,
    controls: MoistureSolidControls,
}

/// Transport and all compatible mechanical outputs at the same accepted time.
pub struct MoistureLayeredSolidStep<Report = MoistureLayeredReport> {
    pub moisture: MoistureLayeredSlab,
    pub transport: Report,
    pub solid: MoistureSolidState,
}

impl MoistureLayeredSolid {
    /// Bond the actual layer cells end to end along x. The common yz section is
    /// `width = area / thickness`; each original cell gets six tetrahedra, with
    /// four shared corner nodes at every cell and material interface. No element
    /// straddles a moisture cell. Source arrays follow layer order and then each
    /// layer's own cell order; zero prescribed DOFs use global node-major xyz.
    ///
    /// The original native equilibrium admits the support/geometry/material
    /// definition immediately. All source frames must map into the same mesh
    /// frame, but each layer may use different principal directions and laws.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        moisture: &MoistureLayeredSlab,
        thickness_m: f64,
        sources: Vec<MoistureLayerSolidSource>,
        fixed_dofs: Vec<usize>,
        controls: MoistureSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        if sources.len() != moisture.layers().len() || sources.is_empty() {
            return Err(MoistureError::Invalid(
                "one mechanical source per original moisture layer required",
            ));
        }
        let count = moisture.activity().len();
        let nodes = count
            .checked_add(1)
            .and_then(|n| n.checked_mul(4))
            .ok_or(MoistureError::Invalid("layered solid node count overflow"))?;
        let tets = count.checked_mul(6).ok_or(MoistureError::Invalid(
            "layered solid element count overflow",
        ))?;
        let mesh_frame = sources[0].axes.mesh_frame;
        if !positive(thickness_m)
            || nodes > controls.assembly.maximum_nodes
            || tets > controls.assembly.maximum_tetrahedra
            || !positive(controls.maximum_absolute_mandel_strain)
            || !positive(controls.maximum_nodal_displacement_m)
            || !positive(controls.mass_tolerance_kg)
            || sources
                .iter()
                .any(|source| source.axes.mesh_frame != mesh_frame)
        {
            return Err(MoistureError::Invalid(
                "layered solid requires bounded geometry, positive controls and one mesh frame",
            ));
        }
        let width = moisture.layers()[0].slab().geometry.cross_section_m2 / thickness_m;
        if !positive(width) {
            return Err(MoistureError::Invalid(
                "layered solid width is unrepresentable",
            ));
        }
        let mut nodes_m = Vec::with_capacity(nodes);
        let mut tetrahedra = Vec::with_capacity(tets);
        let mut reference_cell_volumes_m3 = Vec::with_capacity(count);
        let mut layer_cells = Vec::with_capacity(sources.len());
        let mut bindings = Vec::with_capacity(sources.len());
        let mut offset = 0.0;
        nodes_m.extend([
            [0.0, 0.0, 0.0],
            [0.0, width, 0.0],
            [0.0, width, thickness_m],
            [0.0, 0.0, thickness_m],
        ]);
        for (layer, source) in moisture.layers().iter().zip(sources) {
            poll(cx)?;
            let slab = layer.slab();
            let first = reference_cell_volumes_m3.len();
            let cell_count = slab.moisture().len();
            for i in 1..=cell_count {
                poll(cx)?;
                let x = offset + slab.geometry.length_m * (i as f64 / cell_count as f64);
                let a = nodes_m.len() - 4;
                let b = nodes_m.len();
                let volume = (x - nodes_m[a][0]) * width * thickness_m;
                if !positive(volume) || !x.is_finite() {
                    return Err(MoistureError::Invalid(
                        "layered solid cell geometry is unrepresentable",
                    ));
                }
                nodes_m.extend([
                    [x, 0.0, 0.0],
                    [x, width, 0.0],
                    [x, width, thickness_m],
                    [x, 0.0, thickness_m],
                ]);
                tetrahedra.extend([
                    [a, a + 1, a + 2, b + 2],
                    [a, a + 2, a + 3, b + 2],
                    [a, a + 3, b + 3, b + 2],
                    [a, b + 3, b, b + 2],
                    [a, b, b + 1, b + 2],
                    [a, b + 1, a + 1, b + 2],
                ]);
                reference_cell_volumes_m3.push(volume);
            }
            offset = nodes_m[nodes_m.len() - 4][0];
            layer_cells.push(first..reference_cell_volumes_m3.len());
            bindings.push(MoistureSolidSource::bind(
                slab.field(),
                source.stiffness,
                source.elastic_selection,
                source.expansion_selection,
                source.axes,
            )?);
        }
        let result = Self {
            definition: moisture.definition_identity(),
            nodes_m,
            tetrahedra,
            sources: bindings,
            layer_cells,
            reference_cell_volumes_m3,
            fixed_dofs,
            mesh_frame,
            controls,
        };
        result.resolve(moisture, cx)?;
        Ok(result)
    }

    /// Four yz corners at every global cell boundary, including material joins.
    pub fn reference_nodes_m(&self) -> &[[f64; 3]] {
        &self.nodes_m
    }
    /// Six consecutive tetrahedra per original cell in physical layer order.
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        &self.tetrahedra
    }
    /// Ranges into the accepted elastic and integrated-expansion output arrays.
    pub fn layer_cells(&self) -> &[std::ops::Range<usize>] {
        &self.layer_cells
    }

    /// Resolve the original layers at their own accepted local moisture and
    /// solve ONE globally compatible equilibrium. Swelling mismatch therefore
    /// creates physical stress and support reactions through the native owner.
    /// Wet mass is the exact sum of original dry skeletons plus transported
    /// water. Every tensor, strain, deformation and inertia gate remains native.
    pub fn resolve(
        &self,
        moisture: &MoistureLayeredSlab,
        cx: &Cx<'_>,
    ) -> Result<MoistureSolidState, MoistureError> {
        poll(cx)?;
        if moisture.definition_identity() != self.definition
            || moisture.layers().len() != self.sources.len()
        {
            return Err(MoistureError::Invalid(
                "layered solid requires its original ordered moisture stack",
            ));
        }
        let mut cells = MoistureSolidCells::default();
        for ((layer, source), range) in moisture
            .layers()
            .iter()
            .zip(&self.sources)
            .zip(&self.layer_cells)
        {
            source.append_cells(
                layer.slab().field(),
                &self.reference_cell_volumes_m3[range.clone()],
                self.controls.maximum_absolute_mandel_strain,
                &mut cells,
                cx,
            )?;
        }
        solve_native_moisture_solid(
            &self.nodes_m,
            &self.tetrahedra,
            &self.fixed_dofs,
            self.mesh_frame,
            self.controls,
            moisture.time_s(),
            moisture.wet_mass_kg(),
            cells,
            cx,
        )
    }

    /// Original conservative heterogeneous transport and native mechanics are
    /// published together; late source, strain or geometry refusal accepts none.
    pub fn advance(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        faces: &[MoistureLayeredFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<MoistureLayeredSolidStep, MoistureError> {
        let (moisture, report) = moisture.advance(dt_s, faces, transport, nonlinear, cx)?;
        let solid = self.resolve(&moisture, cx)?;
        poll(cx)?;
        Ok(MoistureLayeredSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }

    /// Finite vapor, heterogeneous moisture and compatible solid mechanics form
    /// one transactional endpoint. Gas pressure is not an invented traction;
    /// gas internal-energy change is not a heat-of-sorption account.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_atmosphere(
        &self,
        atmosphere: &FiniteMoistureAtmosphere,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<(LayeredAtmosphereMoistureCandidate, MoistureSolidState), AtmosphereMoistureError>
    {
        let candidate = atmosphere
            .advance_layered(moisture, dt_s, exposed, transport, nonlinear, controls, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }

    /// Advance the original adaptive water trajectory, then resolve compatible
    /// mechanics at its exact final time. This is an endpoint projection: the
    /// native water controller estimates moisture/activity error, not mechanical
    /// interpolation error. A late mechanical refusal discards the whole window.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_to(
        &self,
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        faces: &[MoistureLayeredFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        MoistureLayeredSolidStep<LayeredMoistureAdvanceReport>,
        MoistureAdvanceError<MoistureError>,
    > {
        let (moisture, report) =
            moisture.advance_adaptive_to(target_s, faces, transport, nonlinear, time, cx)?;
        let solid = self.resolve(&moisture, cx).map_err(|error| match error {
            MoistureError::Cancelled => MoistureAdvanceError::Cancelled,
            error => MoistureAdvanceError::FinalInventory(error),
        })?;
        Ok(MoistureLayeredSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }

    /// Full finite-vapor adaptive window plus its native final solid response.
    /// Reservoir, all water layers and mechanics remain candidates until the
    /// final source/equilibrium/geometry/mass gates pass. There is no mechanical
    /// heat source, stress feedback to diffusion or mode-state remapping.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_finite_atmosphere_adaptive_to(
        &self,
        atmosphere: &FiniteMoistureAtmosphere,
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            FiniteMoistureAtmosphere,
            MoistureLayeredSolidStep<
                LayeredMoistureAdvanceReport<AtmosphereMoistureStep<MoistureLayeredReport>>,
            >,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        let (atmosphere, moisture, report) = atmosphere.advance_layered_adaptive_to(
            moisture, target_s, exposed, transport, nonlinear, physical, time, cx,
        )?;
        let solid = self.resolve(&moisture, cx).map_err(|error| match error {
            MoistureError::Cancelled => MoistureAdvanceError::Cancelled,
            error => MoistureAdvanceError::FinalInventory(error.into()),
        })?;
        Ok((
            atmosphere,
            MoistureLayeredSolidStep {
                moisture,
                transport: report,
                solid,
            },
        ))
    }
}
