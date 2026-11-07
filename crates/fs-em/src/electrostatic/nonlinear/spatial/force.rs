//! Instantaneous shape forces from the accepted nonlinear tetrahedral field.
//!
//! The frozen isotropic constitutive stress is `E tensor D - w*(E) I`.
//! Longitudinal stress uses the independently integrated original storage,
//! not a subtraction of potentially almost equal `E.D` and coenergy terms.
//! These are derivatives with respect to original mesh coordinates, not an
//! inferred finite motion, electrostriction law or mechanical-body mapping.

use super::*;

fn isotropic(material: &SpatialDielectricMaterial) -> bool {
    matches!(
        material,
        SpatialDielectricMaterial::Reversible(_)
            | SpatialDielectricMaterial::Linear {
                law: SpatialDielectricLaw::Isotropic,
                ..
            }
    )
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

/// Original nodal electric loads conjugate to the tetrahedral vertex positions.
/// The source field retains the admitted mixed voltage/total-charge controls.
/// Material cells carry their original total volume charge during the virtual
/// deformation, so rho changes inversely with cell volume. Constitutive source
/// contexts remain fixed and explicitly isotropic.
///
/// A mesh-coordinate force does not imply that every mesh vertex is a solid
/// material point. The caller must supply that mechanical interpretation.
pub struct SpatialElectrostaticShapeForces<'a, 'm> {
    field: &'a SpatialNonlinearElectricField<'m>,
    nodal_force_n: Vec<[f64; 3]>,
}

impl<'a, 'm> SpatialElectrostaticShapeForces<'a, 'm> {
    pub fn field(&self) -> &'a SpatialNonlinearElectricField<'m> {
        self.field
    }

    /// Cartesian forces [N], in the original geometry's coordinate frame.
    pub fn nodal_force_n(&self) -> &[[f64; 3]] {
        &self.nodal_force_n
    }
}

impl<'m> SpatialNonlinearElectricField<'m> {
    fn require_shape_transport(&self, cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
        if self.has_additional_sources() {
            return Err(invalid(
                "shape work with impressed polarization or surface sources requires their transport law",
            ));
        }
        for law in &self.problem.materials {
            poll(cancelled)?;
            if !isotropic(law) {
                return Err(invalid(
                    "shape work requires a frozen isotropic material potential; tensor transport is undeclared",
                ));
            }
        }
        Ok(())
    }

    /// Local forces on vertices 1, 2 and 3. Vertex 0 receives their negative
    /// sum, preserving the original constant-displacement nullspace.
    fn shape_cell(&self, cell: usize) -> Result<[[f64; 3]; 3], Error> {
        let element = self.problem.element(cell);
        let e = self.evaluation.electric[cell];
        let magnitude = finite(e[0].hypot(e[1]).hypot(e[2]))?;
        if magnitude == 0.0 {
            return Ok([[0.0; 3]; 3]);
        }
        let n = e.map(|v| v / magnitude);
        // Choose the least parallel coordinate direction, so the transverse
        // normalization is well conditioned even for an axis-aligned field.
        let mut least = 0;
        for axis in 1..3 {
            if n[axis].abs() < n[least].abs() {
                least = axis;
            }
        }
        let mut seed = [0.0; 3];
        seed[least] = 1.0;
        let transverse = cross(n, seed);
        let norm = transverse[0].hypot(transverse[1]).hypot(transverse[2]);
        let t = transverse.map(|v| v / norm);
        let s = cross(n, t);
        let energy = self.evaluation.cell_energy[cell];
        let coenergy = self.evaluation.cell_coenergy[cell];
        let mut forces = [[0.0; 3]; 3];
        for (force, &gradient) in forces.iter_mut().zip(&element.gradients[1..]) {
            // Negative stress times grad(N): parallel -U, transverse +W*.
            // Integrated U and W* already include the original native volume.
            // A transverse basis avoids recovering small U through ED-W*.
            let parallel = finite(-energy * dot(n, gradient))?;
            let first = finite(coenergy * dot(t, gradient))?;
            let second = finite(coenergy * dot(s, gradient))?;
            for axis in 0..3 {
                force[axis] =
                    finite(n[axis].mul_add(parallel, t[axis].mul_add(first, s[axis] * second)))?;
            }
        }
        Ok(forces)
    }

    /// Differentiate the accepted original mixed electrical functional with
    /// respect to all tetrahedral vertex positions. No additional field solve,
    /// tangent matrix, secant-permittivity surrogate or geometric probe is used.
    ///
    /// Voltage controls and each floating electrode's total charge remain
    /// fixed. Every cell transports its original total volume free charge,
    /// keeping its exact P1 nodal source q/4 fixed. The stationary functional is
    /// `W* - sum(phi*b) - sum(floating_voltage*prescribed_charge)`; voltage-source
    /// work is thereby included, and fixed-charge work equals minus the storage
    /// derivative. Constitutive coefficients and source contexts remain fixed.
    ///
    /// Tensor materials refuse until their strain/frame transport convention
    /// is declared. Electrostriction, external pressure, thermal work and finite
    /// trajectories require separate coupled physics. This derivative inherits
    /// the original field residual and discretization error. Cancellation or
    /// failure leaves the accepted field unchanged and returns no partial load.
    pub fn shape_forces(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectrostaticShapeForces<'_, 'm>, Error> {
        poll(cancelled)?;
        self.require_shape_transport(cancelled)?;
        let n = self.problem.mesh().vertex_count;
        let mut nodal_force_n = vec![[0.0; 3]; n];
        let mut corrections = vec![[0.0; 3]; n];
        for cell in 0..self.problem.mesh().tets.len() {
            poll(cancelled)?;
            let [a, b, c] = self.shape_cell(cell)?;
            let mut zero = [0.0; 3];
            for axis in 0..3 {
                zero[axis] = -finite(a[axis] + b[axis] + c[axis])?;
            }
            let element = self.problem.element(cell);
            for (&vertex, force) in element.vertices.iter().zip([zero, a, b, c]) {
                for axis in 0..3 {
                    add(
                        &mut nodal_force_n[vertex as usize][axis],
                        &mut corrections[vertex as usize][axis],
                        force[axis],
                    )?;
                }
            }
        }
        for (force, correction) in nodal_force_n.iter_mut().zip(corrections) {
            poll(cancelled)?;
            for axis in 0..3 {
                force[axis] = finite(force[axis] + correction[axis])?;
            }
        }
        poll(cancelled)?;
        Ok(SpatialElectrostaticShapeForces {
            field: self,
            nodal_force_n,
        })
    }

    /// First-order electric work [J] along `x(t)=x+t*displacement_m`, with
    /// dimensionless t and the transport convention of [`Self::shape_forces`].
    /// Native cell contractions avoid assembling a global force vector. Each
    /// local displacement is taken relative to vertex 0 before multiplication,
    /// preserving exact zero work for rigid translation, even at large common
    /// displacement offsets. This is an instantaneous derivative, not work
    /// along a finite deformation or a moving-mesh time step.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        self.require_shape_transport(cancelled)?;
        if displacement_m.len() != self.problem.mesh().vertex_count {
            return Err(invalid(
                "one finite virtual displacement per original vertex required",
            ));
        }
        for displacement in displacement_m {
            poll(cancelled)?;
            if displacement.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "one finite virtual displacement per original vertex required",
                ));
            }
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for cell in 0..self.problem.mesh().tets.len() {
            poll(cancelled)?;
            let forces = self.shape_cell(cell)?;
            let element = self.problem.element(cell);
            let u = element.vertices.map(|v| displacement_m[v as usize]);
            for local in 0..3 {
                for axis in 0..3 {
                    add(
                        &mut work,
                        &mut correction,
                        forces[local][axis] * finite(u[local + 1][axis] - u[0][axis])?,
                    )?;
                }
            }
        }
        poll(cancelled)?;
        Ok(finite(work + correction)?)
    }
}
