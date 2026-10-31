//! Instantaneous nonlinear electric loads and the original P1 shape derivative.
//!
//! The frozen isotropic constitutive potential defines stress
//! `E tensor D - w*(E) I`. Here w* is the ORIGINAL material coenergy density,
//! including every crossed source segment. It is not generally E.D/2.
//! Constitutive strain dependence, electrostriction and ambient mechanical
//! pressure require their own coupled potential and are not inferred here.

use super::*;
use crate::{
    force::{
        MaxwellForceLimits,
        contour::{self, ContourError},
    },
    transfer::PlanarFrame,
};
use fs_ga::Wrench;

fn isotropic(material: &DielectricMaterial) -> bool {
    matches!(
        material,
        DielectricMaterial::Reversible(_)
            | DielectricMaterial::Linear {
                law: PlanarDielectricLaw::Isotropic,
                ..
            }
    )
}

/// Apply the isotropic stress through its parallel/perpendicular eigenvalues.
/// Original w and w* are integrated independently, so the longitudinal force
/// retains small storage under a large ED/coenergy bias without subtraction.
fn stress_action(
    e: [f64; 2],
    energy: f64,
    coenergy: f64,
    normal: [f64; 2],
) -> Result<[f64; 2], Error> {
    let magnitude = finite(e[0].hypot(e[1]))?;
    if magnitude == 0.0 {
        return Ok([0.0; 2]);
    }
    let n = [e[0] / magnitude, e[1] / magnitude];
    let parallel = finite(energy * n[0].mul_add(normal[0], n[1] * normal[1]))?;
    let perpendicular = finite(coenergy * (-n[1]).mul_add(normal[0], n[0] * normal[1]))?;
    Ok([
        finite(n[0].mul_add(parallel, n[1] * perpendicular))?,
        finite(n[1].mul_add(parallel, -n[0] * perpendicular))?,
    ])
}

/// Discrete electric forces conjugate to original vertex positions and the
/// actual extrusion thickness. These are derivatives of the accepted mixed
/// electrical functional with its declared voltage/charge controls held fixed.
///
/// Material cells carry their original total free charge during the virtual
/// deformation, so rho scales inversely with cell volume. Constitutive contexts
/// remain fixed and isotropic; this supplies no strain/material-axis law.
/// A mesh-coordinate force is not a claim that every mesh node is a solid node.
pub struct ElectrostaticShapeForces<'a, 'm> {
    field: &'a NonlinearElectricField<'m>,
    nodal_force_n: Vec<[f64; 2]>,
    thickness_force_n: f64,
}

impl<'a, 'm> ElectrostaticShapeForces<'a, 'm> {
    pub fn field(&self) -> &'a NonlinearElectricField<'m> {
        self.field
    }
    /// Forces in the original mesh's planar coordinate frame [N].
    pub fn nodal_force_n(&self) -> &[[f64; 2]] {
        &self.nodal_force_n
    }
    /// Generalized force conjugate to the extrusion thickness [N].
    pub fn thickness_force_n(&self) -> f64 {
        self.thickness_force_n
    }
}

impl<'m> NonlinearElectricField<'m> {
    /// Integrate the nonlinear electric stress on an interior closed contour.
    /// The original geometric owner supplies exterior-side P1 fields, depth,
    /// orientation, first moments and compensated force/torque sums.
    ///
    /// Both incident cells on every contour edge must carry ONE identical
    /// complete isotropic material source, with zero original free charge and
    /// no electrode vertex. The enclosing region may contain charges, material
    /// interfaces and nonlinear response. Mesh-boundary and folded traces
    /// refuse. The returned native Wrench is a net load about the supplied
    /// global origin, not a map of tractions onto an unrepresented solid surface.
    ///
    /// The exterior constitutive potential is fixed with respect to strain;
    /// electrostriction and mechanical/fluid pressure are excluded. This
    /// instantaneous load does not authorize advancing a moving body without
    /// re-solving the field and closing source/field/mechanical work. Do not add
    /// another rho E load for the same enclosed matter. Spatial/contour
    /// refinement is required for accuracy; method agreement is not a bound.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, Error> {
        poll(cancelled)?;
        let Metric2::Planar { thickness } = self.problem.mesh.metric() else {
            return Err(invalid(
                "electric force requires the original planar metric",
            ));
        };
        let materials = &self.problem.materials;
        let mut exterior = None;
        contour::integrate_traction(
            &self.problem.mesh,
            thickness,
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            |inner, outer, vertices, normal_area| {
                if !isotropic(&materials[inner]) || !isotropic(&materials[outer]) {
                    return Err(ContourError::Invalid(
                        "nonlinear electric force needs an explicitly isotropic exterior potential",
                    ));
                }
                let source = *exterior.get_or_insert(&materials[inner]);
                if materials[inner] != *source || materials[outer] != *source {
                    return Err(ContourError::Invalid(
                        "electric contour needs one homogeneous material source on both sides",
                    ));
                }
                if self.free_charge_density_c_m3[inner] != 0.0
                    || self.free_charge_density_c_m3[outer] != 0.0
                    || vertices
                        .iter()
                        .any(|&v| self.nodal_boundary_charge_c[v as usize].is_some())
                {
                    return Err(ContourError::Invalid(
                        "electric contour must avoid free-charge sources and electrode vertices",
                    ));
                }
                let volume = self.problem.elements[outer].volume_m3;
                stress_action(
                    self.evaluation.electric[outer],
                    self.evaluation.cell_energy[outer] / volume,
                    self.evaluation.cell_coenergy[outer] / volume,
                    normal_area,
                )
                .map_err(|_| ContourError::Invalid("electric stress is not representable"))
            },
            cancelled,
        )
        .map_err(|error| match error {
            ContourError::Invalid(reason) => invalid(reason),
            ContourError::Cancelled => ElectrostaticError::Cancelled.into(),
        })
    }

    /// Local shape forces on vertices 1 and 2; the force on vertex 0 is their
    /// negative sum. This preserves cell translation invariance exactly.
    fn shape_cell(&self, cell: usize) -> Result<[[f64; 2]; 2], Error> {
        if !isotropic(&self.problem.materials[cell]) {
            return Err(invalid(
                "shape work requires a frozen isotropic material potential",
            ));
        }
        let element = &self.problem.elements[cell];
        let e = self.evaluation.electric[cell];
        let mut forces = [[0.0; 2]; 2];
        for (out, g) in forces.iter_mut().zip(element.gradient) {
            // d(volume*w*) = volume*(w* div(u) - E.(grad(u) D)).
            // Integrated principal stresses retain the same original discrete
            // potential without losing longitudinal pressure to cancellation.
            let traction = stress_action(
                e,
                self.evaluation.cell_energy[cell],
                self.evaluation.cell_coenergy[cell],
                g,
            )?;
            *out = [-traction[0], -traction[1]];
        }
        Ok(forces)
    }

    /// Differentiate the ORIGINAL accepted mixed-boundary P1 functional with
    /// respect to every vertex and the extrusion thickness, without an extra
    /// field solve, tangent assembly or finite-difference probe.
    ///
    /// Original voltage controls remain fixed; each floating total charge also
    /// remains fixed. Each cell transports its original total volume charge
    /// with its material, preserving the consistent q/3 nodal source. Source
    /// contexts are held fixed: no electrostriction, material axis transport,
    /// charge injection, thermal entropy or mechanical pressure is inferred.
    /// Anisotropic cells refuse until their deformation convention is supplied.
    ///
    /// The envelope derivative uses the admitted primal stationarity, so its
    /// accuracy inherits the field residual and spatial discretization. For
    /// fixed voltages it differentiates coenergy and source work; for fixed
    /// charges it equals minus the stored-energy derivative. The accepted
    /// field retains those controls. Failure/cancellation returns no partial
    /// force vector and never changes the accepted field.
    pub fn shape_forces(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ElectrostaticShapeForces<'_, 'm>, Error> {
        poll(cancelled)?;
        let n = self.problem.mesh.vertices().len();
        let mut nodal_force_n = vec![[0.0; 2]; n];
        let mut corrections = vec![[0.0; 2]; n];
        for (cell, element) in self.problem.elements.iter().enumerate() {
            poll(cancelled)?;
            let [a, b] = self.shape_cell(cell)?;
            let local = [[-finite(a[0] + b[0])?, -finite(a[1] + b[1])?], a, b];
            for (vertex, force) in element.vertices.into_iter().zip(local) {
                for axis in 0..2 {
                    add(
                        &mut nodal_force_n[vertex][axis],
                        &mut corrections[vertex][axis],
                        force[axis],
                    )?;
                }
            }
        }
        for (force, correction) in nodal_force_n.iter_mut().zip(corrections) {
            for axis in 0..2 {
                force[axis] = finite(force[axis] + correction[axis])?;
            }
        }
        let Metric2::Planar { thickness } = self.problem.mesh.metric() else {
            unreachable!()
        };
        let thickness_force_n = finite(self.evaluation.coenergy / thickness)?;
        poll(cancelled)?;
        Ok(ElectrostaticShapeForces {
            field: self,
            nodal_force_n,
            thickness_force_n,
        })
    }

    /// First-order electric work [J] along `x(t)=x+t*displacement_m` and
    /// `depth(t)=depth+t*thickness_change_m`, where t is dimensionless.
    /// The mixed electrode controls, transported cell total charges and frozen
    /// isotropic contexts are EXACTLY those declared by `shape_forces`.
    ///
    /// Direct cell contractions avoid assembling a nodal force vector. Local
    /// displacement differences preserve exact rigid-translation invariance,
    /// including a large common displacement offset. This is an infinitesimal
    /// derivative, not finite mechanical work or a moving-mesh time step.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 2]],
        thickness_change_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        if displacement_m.len() != self.problem.mesh.vertices().len()
            || displacement_m.iter().flatten().any(|v| !v.is_finite())
            || !thickness_change_m.is_finite()
        {
            return Err(invalid(
                "complete finite virtual displacement and depth change required",
            ));
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for (cell, element) in self.problem.elements.iter().enumerate() {
            poll(cancelled)?;
            let forces = self.shape_cell(cell)?;
            let u = element.vertices.map(|v| displacement_m[v]);
            for local in 0..2 {
                for axis in 0..2 {
                    add(
                        &mut work,
                        &mut correction,
                        forces[local][axis] * finite(u[local + 1][axis] - u[0][axis])?,
                    )?;
                }
            }
        }
        let Metric2::Planar { thickness } = self.problem.mesh.metric() else {
            unreachable!()
        };
        add(
            &mut work,
            &mut correction,
            (self.evaluation.coenergy / thickness) * thickness_change_m,
        )?;
        poll(cancelled)?;
        Ok(finite(work + correction)?)
    }
}
