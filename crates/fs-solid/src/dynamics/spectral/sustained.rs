//! Complete-spectrum finite contact forces through the shared physical owner.

use super::*;
use crate::dynamics::sustained::{self as contact, Native, State};
use contact::{
    SolidSustainedContact, SolidSustainedControls, SolidSustainedError, SolidSustainedPath,
    SolidSustainedStep,
};

pub(crate) struct Workspace {
    force: Vec<f64>,
    modal_displacement: Vec<f64>,
    modal_velocity: Vec<f64>,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
}

impl CompleteSpectralSolid {
    /// Advance simultaneous held normal and two-tangent Coulomb forces using
    /// every mode of the original consistent-mass pencil, including rigid modes.
    ///
    /// The complete cached oscillators and basis kernels supply read-only
    /// responses with five reusable O(n) vectors. The existing shared contact
    /// solve, full-coordinate mass/energy audit and whole-path admission are
    /// identical to the dense exact owner's. No modal cutoff, penalty stiffness,
    /// factorization, matrix clone or dissipative surrogate is introduced.
    ///
    /// Original interval/cumulative work gates and the final cancellation poll
    /// still precede acceptance. A basis admitted with loose residuals cannot
    /// bypass those physical gates. Contact loss remains distinct from heat;
    /// this original complete-spectrum owner has no constitutive viscosity.
    pub fn advance_sustained_contacts<P, E>(
        &mut self,
        contacts: &[SolidSustainedContact<'_>],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedPath<'_, '_, Self>) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        SolidSustainedStep<P>,
        SolidSustainedError<E, CompleteSpectralSolidError>,
    > {
        contact::advance(
            self,
            contacts,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }
}

impl Native for CompleteSpectralSolid {
    type Error = CompleteSpectralSolidError;
    type Workspace = Workspace;

    fn state(&self) -> State<'_> {
        State {
            displacement: &self.displacement,
            velocity: &self.velocity,
            mass: &self.mass,
            stiffness: &self.stiffness,
            dt_s: self.dt_s,
        }
    }
    fn admit(&self, cx: &Cx<'_>) -> Result<()> {
        poll(cx)?;
        self.next_interval().map(|_| ())
    }
    fn cancelled() -> Self::Error {
        CompleteSpectralSolidError::Cancelled
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        matches!(error, CompleteSpectralSolidError::Cancelled)
    }
    fn invalid(what: &'static str) -> Self::Error {
        super::invalid(what)
    }
    fn workspace_bytes(&self) -> u128 {
        // Five response vectors plus the final original affine acceptance's
        // scalar and vector staging. Prepared O(n²) matrices are only borrowed.
        8 * (16 * self.free_dofs.len() as u128 + 4)
    }
    fn response_multiply_terms(&self) -> u128 {
        let n = self.free_dofs.len() as u128;
        // Original full native interval bound also covers the final replay;
        // adapter contact reductions and full-M/K audit are charged separately.
        16 * n * n + 512 * n
    }
    fn workspace(&self) -> Result<Workspace> {
        let n = self.free_dofs.len();
        Ok(Workspace {
            force: zeros(n)?,
            modal_displacement: zeros(n)?,
            modal_velocity: zeros(n)?,
            displacement: zeros(n)?,
            velocity: zeros(n)?,
        })
    }
    fn evaluate(&self, force: &[f64], scratch: &mut Workspace, cx: &Cx<'_>) -> Result<()> {
        project_into(&self.basis, force, &mut scratch.force, cx)?;
        for i in 0..self.free_dofs.len() {
            poll(cx)?;
            let [q, v, _] = self.affine_mode(i, scratch.force[i], 0.0);
            scratch.modal_displacement[i] = q;
            scratch.modal_velocity[i] = v;
        }
        reconstruct_into(
            &self.basis,
            &scratch.modal_displacement,
            &mut scratch.displacement,
            cx,
        )?;
        reconstruct_into(
            &self.basis,
            &scratch.modal_velocity,
            &mut scratch.velocity,
            cx,
        )
    }
    fn displacement<'a>(&self, scratch: &'a Workspace) -> &'a [f64] {
        &scratch.displacement
    }
    fn velocity<'a>(&self, scratch: &'a Workspace) -> &'a [f64] {
        &scratch.velocity
    }
    fn work(&self, force: &[f64], scratch: &Workspace, cx: &Cx<'_>) -> Result<(f64, f64)> {
        let (mut supplied, mut correction) = (0.0, 0.0);
        for (i, &force) in force.iter().enumerate() {
            poll(cx)?;
            add(
                &mut supplied,
                &mut correction,
                force.mul_add(scratch.displacement[i] - self.displacement[i], 0.0),
            );
        }
        if !supplied.is_finite() {
            return Err(invalid("unrepresentable held source work"));
        }
        Ok((supplied, 0.0))
    }
    fn accept(
        &mut self,
        scratch: Workspace,
        force: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        // Replays the same original basis and oscillator kernels once, with
        // unchanged accepted state and identical held inputs. This preserves
        // the native source/energy/cumulative gates and sole publication point.
        // Drop response scratch before its O(n) native staging is allocated.
        drop(scratch);
        self.advance(force, cx)
    }
}
