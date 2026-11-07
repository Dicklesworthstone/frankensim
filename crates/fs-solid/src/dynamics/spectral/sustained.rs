//! Complete-spectrum finite contact forces through the shared physical owner.

use super::*;
use crate::dynamics::sustained::{self as contact, Native, State, Work};
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
    /// Resolve distinct static/kinetic principal strengths on the original
    /// complete spectrum. Each kinetic pair uses its static contact's material
    /// axes; coupled breakaway re-solves share budgets and one native commit.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_static_kinetic_elliptic_contacts<P, E>(
        &mut self,
        contacts: &[contact::elliptic::SolidSustainedEllipticContact<'_>],
        kinetic_coefficients_mu: &[[f64; 2]],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(
            &contact::elliptic::SolidSustainedEllipticPath<'_, '_, Self>,
        ) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        contact::elliptic::SolidSustainedEllipticStep<P, ExactLinearSolidStep>,
        SolidSustainedError<E, CompleteSpectralSolidError>,
    > {
        contact::elliptic::advance_static_kinetic(
            self,
            contacts,
            kinetic_coefficients_mu,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

    /// Locate a stop, static recapture or reversal along one rigid contact
    /// trace, retaining every native elastic mode and the original full energy
    /// and cumulative-work gates. See the dense owner's event method for the
    /// compressive, fixed-frame, collinear and constant-load admission domain.
    pub fn advance_collinear_coulomb_events<P, E>(
        &mut self,
        law: &contact::events::CollinearCoulombContact<'_>,
        initial_force_n: [f64; 3],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(
            &contact::events::SolidCoulombEventPath<'_, '_, Self>,
        ) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        contact::events::SolidCoulombEventStep<P>,
        SolidSustainedError<E, CompleteSpectralSolidError>,
    > {
        contact::events::advance(
            self,
            law,
            initial_force_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

    /// Advance an oriented Coulomb ellipse through the same original complete
    /// response, work gates and atomic acceptance as sustained disk contact.
    /// Physical contacts/forces reach the mandatory whole-path callback; the
    /// returned nonlinear diagnostics explicitly retain solver coordinates.
    pub fn advance_sustained_elliptic_contacts<P, E>(
        &mut self,
        contacts: &[contact::elliptic::SolidSustainedEllipticContact<'_>],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(
            &contact::elliptic::SolidSustainedEllipticPath<'_, '_, Self>,
        ) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        contact::elliptic::SolidSustainedEllipticStep<P, ExactLinearSolidStep>,
        SolidSustainedError<E, CompleteSpectralSolidError>,
    > {
        contact::elliptic::advance(
            self,
            contacts,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

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

    /// Use separate static and kinetic bounds on the complete modal response.
    /// Contact coefficients are static; `kinetic_coefficients` use the same
    /// input order and must lie between zero and their static bounds. Incoming
    /// slip selects kinetic friction; resting contacts that break away trigger
    /// a coupled re-solve under the shared response/work budgets. The original
    /// physical audits and mandatory whole-path admission still gate acceptance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_static_kinetic_contacts<P, E>(
        &mut self,
        contacts: &[SolidSustainedContact<'_>],
        kinetic_coefficients: &[f64],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedPath<'_, '_, Self>) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        SolidSustainedStep<P>,
        SolidSustainedError<E, CompleteSpectralSolidError>,
    > {
        contact::advance_with_kinetic(
            self,
            contacts,
            Some(kinetic_coefficients),
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }
}

impl Native for CompleteSpectralSolid {
    type Owner = Self;
    type Step = ExactLinearSolidStep;
    type Error = CompleteSpectralSolidError;
    type Workspace = Workspace;

    fn owner(&self) -> &Self {
        self
    }
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
    fn work(&self, force: &[f64], scratch: &Workspace, cx: &Cx<'_>) -> Result<Work> {
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
        Ok(Work {
            applied_j: supplied,
            viscosity_j: 0.0,
            support_j: 0.0,
        })
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

impl contact::events::EventNative for CompleteSpectralSolid {
    fn trace_response(
        &self,
        weights: [&[f64]; 3],
        cx: &Cx<'_>,
    ) -> Result<contact::events::TraceResponse> {
        let mut modal_columns: [Vec<f64>; 3] = core::array::from_fn(|_| Vec::new());
        let mut columns: [Vec<f64>; 3] = core::array::from_fn(|_| Vec::new());
        for a in 0..3 {
            modal_columns[a] = project(&self.basis, weights[a], cx)?;
            if modal_columns[a]
                .iter()
                .zip(&self.oscillators)
                .any(|(force, oscillator)| *force != 0.0 && oscillator.lambda != 0.0)
            {
                return Err(invalid(
                    "event contact must excite only original zero eigenvalue modes",
                ));
            }
            columns[a] = reconstruct(&self.basis, &modal_columns[a], cx)?;
        }
        Ok(contact::events::TraceResponse {
            columns,
            modal_columns: Some(modal_columns),
            support_reaction_columns: None,
            acceleration_offset: [0.0; 3],
        })
    }

    fn add_contact_moments(
        &self,
        scratch: &mut Workspace,
        response: &contact::events::TraceResponse,
        path: &contact::events::ContactForcePath<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<f64> {
        use fs_math::dd::Dd;
        let columns = response
            .modal_columns
            .as_ref()
            .ok_or_else(|| invalid("complete modal contact response required"))?;
        for i in 0..self.free_dofs.len() {
            poll(cx)?;
            let mut q = Dd::from_f64(scratch.modal_displacement[i]);
            let mut v = Dd::from_f64(scratch.modal_velocity[i]);
            for a in 0..3 {
                let column = Dd::from_f64(columns[a][i]);
                q = q + column * path.displacement_moment[a];
                v = v + column * path.impulse[a];
            }
            scratch.modal_displacement[i] = contact::finite(q).map_err(invalid)?;
            scratch.modal_velocity[i] = contact::finite(v).map_err(invalid)?;
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
        )?;
        Ok(0.0)
    }

    fn accept_force_path(
        &mut self,
        scratch: Workspace,
        supplied_j: f64,
        viscosity_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        if viscosity_j != 0.0 {
            return Err(invalid(
                "complete-spectrum owner has no constitutive viscosity",
            ));
        }
        let (step, time_s) = self.next_interval()?;
        CompleteSpectralSolid::accept_force_path(
            self,
            scratch.modal_displacement,
            scratch.modal_velocity,
            scratch.displacement,
            scratch.velocity,
            supplied_j,
            step,
            time_s,
            cx,
        )
    }
}
