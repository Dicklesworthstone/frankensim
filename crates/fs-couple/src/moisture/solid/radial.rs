//! Original radial moisture drives a compatible native cylindrical solid.
//!
//! A conforming polygonal cylinder retains every original radial interface.
//! Each shell's analytic dry inventory is distributed by its actual tetrahedral
//! reference volumes. The resulting density correction is exactly the admitted
//! curved/polygonal volume discrepancy; it vanishes with angular refinement.
//! Geometry, free strain, local constitutive state and wet inertia stay distinct.

use super::*;
use crate::moisture::adaptive::{
    MoistureAdvanceError, MoistureAdvanceReport, MoistureTimeControls,
};
use crate::moisture::radial::{MoistureRadial, RadialMoistureSymmetry};

#[path = "radial/finite.rs"]
mod finite;
#[path = "radial/mesh.rs"]
mod mesh;
#[path = "radial/solution.rs"]
mod solution;

/// Resolution of the actual cylindrical reference geometry. The native solid's
/// assembly controls separately bound total nodes, tetrahedra and free DOFs.
#[derive(Clone, Copy, Debug)]
pub struct MoistureRadialSolidMesh {
    /// At least three circumferential sectors, with a shared periodic seam.
    pub angular_divisions: usize,
    /// Positive number of axial intervals; diffusion remains radially uniform.
    pub axial_divisions: usize,
    /// Relative difference of each polygonal and analytic curved shell volume.
    /// Finite in (0,1). This is a geometric admission limit, not an error bound
    /// on displacement, stress or vibration. Refine all relevant mesh directions.
    pub maximum_relative_shell_volume_error: f64,
}

/// Radial source binding, actual three-dimensional geometry and native supports.
/// Fixed material axes are in the cylinder's xyz mesh frame; its axis is z.
/// A radially varying moisture field need not produce an axisymmetric elastic
/// response when the authored material or supports are anisotropic.
#[derive(Clone, Debug)]
pub struct MoistureRadialSolid {
    source: MoistureSolidSource,
    mesh: mesh::CylinderMesh,
    fixed_dofs: Vec<usize>,
    controls: MoistureSolidControls,
}

/// Transport and its dependent native mechanical response at one accepted time.
pub struct MoistureRadialSolidStep<Report = MoistureStepReport> {
    pub moisture: MoistureRadial,
    pub transport: Report,
    pub solid: MoistureSolidState,
}

impl MoistureRadialSolid {
    /// Build a conforming native tetrahedral solid on a cylinder or annulus.
    /// The supplied material frame is held in xyz; rotating radial/tangential
    /// grain frames require a spatial constitutive owner and are not inferred.
    ///
    /// The diffusion grid and dry skeleton stay on their original chart. Small
    /// compatible strain changes physical geometry and density, not the water
    /// diffusion metric. Stress-assisted diffusion, sorption heat and finite
    /// deformation are separate couplings. Spherical geometry refuses explicitly.
    ///
    /// Nodes are z-plane major, then increasing radius, then angle from +x.
    /// Each plane starts with a center node for a solid cylinder; annuli have no
    /// center node. Supports are zero displacements in node-major xyz order.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        moisture: &MoistureRadial,
        resolution: MoistureRadialSolidMesh,
        stiffness: SolidStiffnessModel,
        elastic_selection: MaterialPropertySelection,
        expansion_selection: MaterialPropertySelection,
        axes: HygroscopicAxes,
        fixed_dofs: Vec<usize>,
        controls: MoistureSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let RadialMoistureSymmetry::Cylinder { length_m } = moisture.geometry().symmetry else {
            return Err(MoistureError::Invalid(
                "radial solid requires a cylinder or annulus",
            ));
        };
        if !positive(controls.maximum_absolute_mandel_strain)
            || !positive(controls.maximum_nodal_displacement_m)
            || !positive(controls.mass_tolerance_kg)
        {
            return Err(MoistureError::Invalid(
                "positive native mechanical admission limits required",
            ));
        }
        let mesh = mesh::cylinder(
            &moisture.geometry().faces_m,
            length_m,
            resolution.angular_divisions,
            resolution.axial_divisions,
            controls.assembly.maximum_nodes,
            controls.assembly.maximum_tetrahedra,
            resolution.maximum_relative_shell_volume_error,
            || poll(cx),
        )
        .map_err(|error| match error {
            mesh::MeshError::Invalid(message) => MoistureError::Invalid(message),
            mesh::MeshError::Capacity => MoistureError::Invalid("radial solid allocation refused"),
            mesh::MeshError::Checkpoint(error) => error,
        })?;
        let source = MoistureSolidSource::bind_original(
            moisture.material_card(),
            moisture.reference(),
            moisture.identity(),
            stiffness,
            elastic_selection,
            expansion_selection,
            axes,
        )?;
        let result = Self {
            source,
            mesh,
            fixed_dofs,
            controls,
        };
        result.resolve(moisture, cx)?;
        Ok(result)
    }

    pub fn reference_nodes_m(&self) -> &[[f64; 3]] {
        &self.mesh.nodes_m
    }
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        &self.mesh.tetrahedra
    }
    /// Original radial cell of every native element, including across the seam.
    pub fn shell_for_element(&self) -> &[usize] {
        &self.mesh.shell_for_element
    }
    /// Actual relative geometric discrepancy used in conserving shell inventory.
    pub const fn maximum_relative_shell_volume_error(&self) -> f64 {
        self.mesh.maximum_relative_volume_error
    }

    /// Re-evaluate the original local elastic and swelling sources, then use the
    /// existing solid owner for equilibrium, stress, deformation and K_e+K_g/M.
    /// Nonuniform swelling produces compatible stress and support reactions;
    /// it is never applied as independent, incompatible shell displacements.
    pub fn resolve(
        &self,
        moisture: &MoistureRadial,
        cx: &Cx<'_>,
    ) -> Result<MoistureSolidState, MoistureError> {
        poll(cx)?;
        if self.source.definition != moisture.identity()
            || moisture.moisture().len() != self.mesh.shell_volumes_m3.len()
        {
            return Err(MoistureError::Invalid(
                "radial solid requires its original material, radial mesh and datum",
            ));
        }
        let count = moisture.moisture().len();
        let elements = self.mesh.tetrahedra.len();
        let mut cells = MoistureSolidCells::default();
        let mut materials = Vec::with_capacity(count);
        let mut strains = Vec::with_capacity(count);
        cells.elastic.reserve(count);
        cells.expansion.reserve(count);
        cells.materials.reserve(elements);
        cells.strains.reserve(elements);
        cells.element_mass_kg.reserve(elements);
        for (shell, &m) in moisture.moisture().iter().enumerate() {
            poll(cx)?;
            let point = moisture
                .transport_context()
                .query()
                .with("moisture_content", m)
                .map_err(|_| MoistureError::Invalid("unrepresentable radial mechanical state"))?;
            let state = self
                .source
                .elastic
                .at_query_point(moisture.material_card(), &point)?;
            let expansion = integrate_orthotropic_hygroscopic_expansion(
                moisture.material_card(),
                &moisture.reference().query(),
                &point,
                self.source.expansion_selection.clone(),
            )?;
            let wet_mass = moisture.dry_masses_kg()[shell] * (1.0 + m);
            let volume = self.mesh.shell_volumes_m3[shell];
            let density = wet_mass / volume;
            if !positive(wet_mass) || !positive(density) {
                return Err(MoistureError::Invalid(
                    "radial wet element density is unrepresentable",
                ));
            }
            let mut hash = DomainHasher::new("org.frankensim.fs-couple.radial-moisture-solid.v1");
            hash.update(moisture.identity().as_bytes());
            hash.update(state.identity().as_bytes());
            hash.update(&(shell as u64).to_le_bytes());
            hash.update(&wet_mass.to_bits().to_le_bytes());
            hash.update(&volume.to_bits().to_le_bytes());
            let material = self.source.stiffness.material(
                &state,
                density,
                hash.finalize(),
                self.controls.maximum_absolute_mandel_strain,
            )?;
            let strain = TetEigenstrainState::try_from_principal_strains(
                &material,
                *expansion.free_principal_strains(),
                self.source.axes.principal_to_mesh,
                self.source.axes.principal_frame,
                self.source.axes.mesh_frame,
                expansion.identity(),
                self.controls.maximum_absolute_mandel_strain,
            )?;
            materials.push(material);
            strains.push(strain);
            cells.elastic.push(state);
            cells.expansion.push(expansion);
        }
        for (element, (&shell, &volume)) in self
            .mesh
            .shell_for_element
            .iter()
            .zip(&self.mesh.element_volumes_m3)
            .enumerate()
        {
            if element % 64 == 0 {
                poll(cx)?;
            }
            let mass = moisture.dry_masses_kg()[shell]
                * (1.0 + moisture.moisture()[shell])
                * (volume / self.mesh.shell_volumes_m3[shell]);
            if !positive(mass) {
                return Err(MoistureError::Invalid(
                    "radial element inventory is unrepresentable",
                ));
            }
            cells.materials.push(materials[shell].clone());
            cells.strains.push(strains[shell]);
            cells.element_mass_kg.push(mass);
        }
        solve_native_moisture_solid(
            &self.mesh.nodes_m,
            &self.mesh.tetrahedra,
            &self.fixed_dofs,
            self.source.axes.mesh_frame,
            self.controls,
            moisture.time_s(),
            moisture.wet_mass_kg(),
            cells,
            cx,
        )
    }

    /// Publish one conservative radial water step and all dependent native
    /// mechanics together. A late constitutive or mechanical refusal preserves
    /// the entire original moisture checkpoint.
    pub fn advance(
        &self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<MoistureRadialSolidStep, MoistureError> {
        let (moisture, report) = moisture.advance(dt_s, faces, transport, nonlinear, cx)?;
        let solid = self.resolve(&moisture, cx)?;
        poll(cx)?;
        Ok(MoistureRadialSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }

    /// Resolve native mechanics only after both conservative SDIRK stages are
    /// admitted. The original state survives any stage or mechanical refusal.
    pub fn advance_sdirk2(
        &self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<MoistureRadialSolidStep<MoistureSdirk2Report>, MoistureError> {
        let (moisture, report) = moisture.advance_sdirk2(dt_s, faces, transport, nonlinear, cx)?;
        let solid = self.resolve(&moisture, cx)?;
        poll(cx)?;
        Ok(MoistureRadialSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }

    /// Resolve compatible mechanics after an accuracy-controlled complete
    /// backward-Euler duration. Any late refusal discards every private prefix.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_adaptive_to(
        &self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[MoistureFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        MoistureRadialSolidStep<MoistureAdvanceReport<MoistureStepReport, 2>>,
        MoistureAdvanceError<MoistureError>,
    > {
        let (moisture, report) =
            moisture.advance_adaptive_to(target_s, faces, transport, nonlinear, time, cx)?;
        self.resolve_adaptive(moisture, report, cx)
    }

    /// The same atomic solid endpoint after second-order time-error control.
    /// The original reference-chart diffusion and small-strain mechanical
    /// assumptions remain; no evolving mechanical trajectory is reconstructed.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[MoistureFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        MoistureRadialSolidStep<MoistureAdvanceReport<MoistureSdirk2Report, 2>>,
        MoistureAdvanceError<MoistureError>,
    > {
        let (moisture, report) =
            moisture.advance_sdirk2_adaptive_to(target_s, faces, transport, nonlinear, time, cx)?;
        self.resolve_adaptive(moisture, report, cx)
    }

    fn resolve_adaptive<Report>(
        &self,
        moisture: MoistureRadial,
        report: Report,
        cx: &Cx<'_>,
    ) -> Result<MoistureRadialSolidStep<Report>, MoistureAdvanceError<MoistureError>> {
        let map = |error| match error {
            MoistureError::Cancelled => MoistureAdvanceError::Cancelled,
            error => MoistureAdvanceError::FinalInventory(error),
        };
        let solid = self.resolve(&moisture, cx).map_err(map)?;
        poll(cx).map_err(map)?;
        Ok(MoistureRadialSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }

    /// Resolve the native mechanical endpoint against the actual finite-air
    /// humidity solution. A late material or mechanical refusal publishes
    /// neither the radial water field nor the atmosphere's changed inventory.
    /// Gas pressure is not automatically a surface traction; sorption heat and
    /// mechanical feedback into the reference diffusion metric remain separate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_atmosphere(
        &self,
        atmosphere: &atmosphere::FiniteMoistureAtmosphere,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: atmosphere::AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        (
            crate::moisture::radial::atmosphere::RadialAtmosphereMoistureCandidate,
            MoistureSolidState,
        ),
        atmosphere::AtmosphereMoistureError,
    > {
        let candidate =
            atmosphere.advance_radial(moisture, dt_s, faces, transport, nonlinear, controls, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }

    /// Keep both chambers, their signed face exchanges and the dependent native
    /// solid response in one accepted candidate, including steady through-flow.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired_atmospheres(
        &self,
        first: &atmosphere::FiniteMoistureAtmosphere,
        second: &atmosphere::FiniteMoistureAtmosphere,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<atmosphere::pair::AtmosphereMoistureFace>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: atmosphere::pair::PairedAtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        (
            crate::moisture::radial::atmosphere::PairedRadialAtmosphereMoistureCandidate,
            MoistureSolidState,
        ),
        atmosphere::AtmosphereMoistureError,
    > {
        let candidate = first.advance_paired_radial(
            second, moisture, dt_s, faces, transport, nonlinear, controls, cx,
        )?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }
}
