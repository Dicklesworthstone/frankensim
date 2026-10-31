//! One physical shell point for reciprocal axial, normal and other force ports.
use super::*;
use crate::thin_plate::PlateMaterialModel;
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Geometry in the original shell mesh's global XYZ frame.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellPointGeometry {
    /// Original assembly part; zero for a single material shell.
    pub shell_component: usize,
    /// Original physical mesh node within that part.
    pub node: usize,
    /// Positive displacement/force direction, a finite unit vector.
    pub direction: [f64; 3],
    /// Vector [m] from the shell node to the physical contact point.
    pub arm_m: [f64; 3],
}

impl ShellPointGeometry {
    /// Admit explicit geometry without normalizing or inventing an orientation.
    pub fn validate(&self) -> Result<(), ReciprocalStringError> {
        self.nodal_row().map(|_| ())
    }

    /// Work-conjugate force and moment [d, r cross d]. Its dot product with
    /// [U, Theta] equals d dot (U + Theta cross r).
    pub(crate) fn nodal_row(&self) -> Result<[f64; 6], ReciprocalStringError> {
        fs_plate::shell::point::nodal_point_row(self.direction, self.arm_m)
            .map_err(|e| invalid(e.reason()))
    }
}

/// Actual original-body map and physical material at one retained shell node.
/// This does not certify a particular contact law or infer a surface coating.
#[derive(Clone, Debug)]
pub struct MaterialShellPoint<'a> {
    /// Original supplied-body component zero, before the string prefix.
    pub body: ModalAttachment,
    /// Same complete source state on EVERY element incident to the node.
    pub material: &'a ResolvedMaterialStatePoint,
    /// Minimum physical section thickness over those incident elements [m].
    pub thickness_m: f64,
    /// Exact elastic source selection when every incident section agrees.
    /// `None` permits a generic force port across differing selections, but
    /// cannot supply one unambiguous constitutive source to a contact law.
    pub elastic_model: Option<PlateMaterialModel>,
    /// Every incident section was actually compiled with isotropic elasticity.
    /// Extra isotropic answers on an orthotropic card do not establish this.
    pub isotropic_elasticity: bool,
}

impl MaterialStringBankSource {
    /// The original material component and its physical GLOBAL modal rows.
    /// Point-force and torque queries share this same source/basis selection.
    pub(super) fn physical_shell_component(&self, component: usize)
        -> Result<(&ResolvedShearShell, &NodalDirectors, &[Vec<f64>]), ReciprocalStringError>
    {
        if let Some(body) = &self.material_shell_assembly {
            let part = body.specimen.components.get(component)
                .ok_or_else(|| invalid("shell port names an absent assembly part"))?;
            let image = self.material_shell_assembly_image()
                .ok_or_else(|| invalid("shell port lost the assembled body basis"))?;
            let modes = image.physical_modes().get(component)
                .ok_or_else(|| invalid("shell port lost its original part basis"))?;
            Ok((&part.specimen, &part.directors, modes.as_slice()))
        } else if let Some(body) = &self.material_shell {
            if component != 0 {
                return Err(invalid("single-shell port requires component zero"));
            }
            let image = self.material_shell_image()
                .ok_or_else(|| invalid("shell port lost the original body basis"))?;
            Ok((&body.specimen, &body.directors, image.physical_modes()))
        } else {
            Err(invalid("physical shell port requires a sourced shell or assembled-shell body"))
        }
    }

    /// Project a force at an actual shell/assembly node into the original common
    /// body basis. Reactions use its transpose in the existing coupled owner.
    /// Different incident material states, absent nodes and non-shell receivers
    /// refuse; an assembly seam never borrows another part's material identity.
    pub fn shell_point(
        &self,
        geometry: &ShellPointGeometry,
        gate: &CancelGate,
    ) -> Result<MaterialShellPoint<'_>, ReciprocalStringError> {
        poll(gate)?;
        let row = geometry.nodal_row()?;
        let (specimen, _, modes) = self.physical_shell_component(geometry.shell_component)?;
        let mesh = specimen.mesh();
        if geometry.node >= mesh.nodes.len() {
            return Err(invalid("shell point names an absent physical node"));
        }
        let size = mesh
            .nodes
            .len()
            .checked_mul(6)
            .ok_or_else(|| invalid("shell point footprint overflow"))?;
        let start = geometry.node * 6;
        let mut footprint_row = vec![0.0; size];
        footprint_row[start..start + 6].copy_from_slice(&row);
        let mut material = None;
        super::super::interfaces::footprint(
            mesh.nodes.len(),
            &mesh.tris,
            &footprint_row,
            6,
            |element| specimen.material_at_element(element),
            &mut material,
            Some(gate),
        )?;
        let material =
            material.ok_or_else(|| invalid("shell point has no actual material surface"))?;
        let mut thickness_m = f64::INFINITY;
        let mut elastic_model = None;
        let mut common_elastic_model = true;
        let mut isotropic_elasticity = true;
        for region in specimen.regions() {
            poll(gate)?;
            if region
                .input()
                .region
                .region
                .triangle_indices
                .iter()
                .any(|&element| mesh.tris[element].contains(&geometry.node))
            {
                thickness_m = thickness_m.min(region.thickness_m());
                let selected = region.input().region.model;
                if elastic_model.is_some_and(|previous| previous != selected) {
                    common_elastic_model = false;
                }
                elastic_model = Some(selected);
                isotropic_elasticity &= selected.is_isotropic();
            }
        }
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(invalid(
                "shell point has no positive physical section thickness",
            ));
        }
        if !common_elastic_model { elastic_model = None; }
        let shapes = modes
            .iter()
            .map(|mode| {
                poll(gate)?;
                let nodal = mode
                    .get(start..start + 6)
                    .ok_or_else(|| invalid("shell point basis lost its physical node"))?;
                let shape = nodal.iter().zip(row).map(|(q, f)| q * f).sum::<f64>();
                if !shape.is_finite() {
                    return Err(invalid("shell point body projection overflow"));
                }
                Ok(shape)
            })
            .collect::<Result<Vec<_>, ReciprocalStringError>>()?;
        Ok(MaterialShellPoint {
            body: ModalAttachment {
                component: 0,
                shapes,
            },
            material,
            thickness_m,
            elastic_model,
            isotropic_elasticity,
        })
    }
}
