//! The original hygroscopic solid owner on actual three-dimensional water cells.

use super::*;
use crate::moisture::cartesian::{MoistureCartesian, MoistureCartesianReport};

/// Transport and compatible mechanical response publish as one endpoint.
pub struct MoistureCartesianSolidStep {
    pub moisture: MoistureCartesian,
    pub transport: MoistureCartesianReport,
    pub solid: MoistureSolidState,
}

impl MoistureSolid {
    /// Use the complete Cartesian reference box with six conforming tetrahedra
    /// per moisture cell. This preserves transverse and through-thickness
    /// moisture gradients in source elasticity, free swelling, stress and wet
    /// inertia; no strip averaging or separate material evaluator is involved.
    /// Node order is x-fastest on the (nx+1)*(ny+1)*(nz+1) corner lattice.
    /// Tetrahedra are consecutive in the original moisture-cell order.
    #[allow(clippy::too_many_arguments)]
    pub fn from_cartesian(
        moisture: &MoistureCartesian,
        stiffness: SolidStiffnessModel,
        elastic_selection: MaterialPropertySelection,
        expansion_selection: MaterialPropertySelection,
        axes: HygroscopicAxes,
        fixed_dofs: Vec<usize>,
        controls: MoistureSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let geometry = moisture.geometry();
        let cells = geometry.cells;
        let counts = cells.map(|n| n + 1);
        let nodes = counts
            .iter()
            .try_fold(1usize, |n, &m| n.checked_mul(m))
            .ok_or(MoistureError::Invalid(
                "Cartesian solid node count overflow",
            ))?;
        let tets = moisture
            .moisture()
            .len()
            .checked_mul(6)
            .ok_or(MoistureError::Invalid("Cartesian solid tet count overflow"))?;
        if nodes > controls.assembly.maximum_nodes
            || tets > controls.assembly.maximum_tetrahedra
            || !positive(controls.maximum_absolute_mandel_strain)
            || !positive(controls.maximum_nodal_displacement_m)
            || !positive(controls.mass_tolerance_kg)
            || axes.mesh_frame != geometry.frame
        {
            return Err(MoistureError::Invalid(
                "Cartesian solid needs bounded mesh, positive admission limits and the original grid frame",
            ));
        }
        let node = |x: usize, y: usize, z: usize| x + counts[0] * (y + counts[1] * z);
        let mut nodes_m = Vec::with_capacity(nodes);
        for z in 0..counts[2] {
            for y in 0..counts[1] {
                for x in 0..counts[0] {
                    poll(cx)?;
                    nodes_m.push(std::array::from_fn(|k| {
                        geometry.lengths_m[k] * ([x, y, z][k] as f64 / cells[k] as f64)
                    }));
                }
            }
        }
        let mut tetrahedra = Vec::with_capacity(tets);
        let mut volumes = Vec::with_capacity(moisture.moisture().len());
        for z in 0..cells[2] {
            for y in 0..cells[1] {
                for x in 0..cells[0] {
                    poll(cx)?;
                    let a = [
                        node(x, y, z),
                        node(x, y + 1, z),
                        node(x, y + 1, z + 1),
                        node(x, y, z + 1),
                    ];
                    let b = [
                        node(x + 1, y, z),
                        node(x + 1, y + 1, z),
                        node(x + 1, y + 1, z + 1),
                        node(x + 1, y, z + 1),
                    ];
                    tetrahedra.extend([
                        [a[0], a[1], a[2], b[2]],
                        [a[0], a[2], a[3], b[2]],
                        [a[0], a[3], b[3], b[2]],
                        [a[0], b[3], b[0], b[2]],
                        [a[0], b[0], b[1], b[2]],
                        [a[0], b[1], a[1], b[2]],
                    ]);
                    let volume = (0..3)
                        .map(|k| nodes_m[b[2]][k] - nodes_m[a[0]][k])
                        .product();
                    if !positive(volume) {
                        return Err(MoistureError::Invalid(
                            "Cartesian solid cell volume is unrepresentable",
                        ));
                    }
                    volumes.push(volume);
                }
            }
        }
        Self::bind_mesh(
            moisture.field(),
            nodes_m,
            tetrahedra,
            volumes,
            stiffness,
            elastic_selection,
            expansion_selection,
            axes,
            fixed_dofs,
            controls,
            cx,
        )
    }

    /// Re-evaluate every original 3D cell's elasticity and swelling at its own
    /// accepted moisture; the shared native solid computes compatible shape.
    pub fn resolve_cartesian(
        &self,
        moisture: &MoistureCartesian,
        cx: &Cx<'_>,
    ) -> Result<MoistureSolidState, MoistureError> {
        self.resolve_field(moisture.field(), cx)
    }

    /// Accept all six water-face transfers and their native mechanical response
    /// together. A late source or mechanical refusal publishes neither candidate.
    pub fn advance_cartesian(
        &self,
        moisture: &MoistureCartesian,
        dt_s: f64,
        faces: &[MoistureFace; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<MoistureCartesianSolidStep, MoistureError> {
        let (moisture, report) = moisture.advance(dt_s, faces, transport, nonlinear, cx)?;
        let solid = self.resolve_cartesian(&moisture, cx)?;
        poll(cx)?;
        Ok(MoistureCartesianSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }
}
