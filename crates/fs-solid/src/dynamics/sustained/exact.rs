//! Existing dense exact owner behind the shared sustained-contact admission.

use super::*;

pub(crate) struct Workspace {
    lifted: Vec<f64>,
    balanced: Vec<f64>,
    candidate: Vec<f64>,
    integrals: Vec<f64>,
}

impl Native for ExactLinearSolid {
    type Error = ExactLinearSolidError;
    type Workspace = Workspace;

    fn state(&self) -> State<'_> {
        State {
            displacement: self.displacement_m(),
            velocity: self.velocity_mps(),
            mass: &self.mass,
            stiffness: &self.stiffness,
            dt_s: self.dt_s,
        }
    }
    fn admit(&self, cx: &Cx<'_>) -> Result<(), Self::Error> {
        super::super::poll(cx)?;
        if self.affine_force_scale_n.is_some() {
            return Err(super::super::invalid(
                "held preparation required for sustained contact",
            ));
        }
        self.next_interval().map(|_| ())
    }
    fn cancelled() -> Self::Error {
        ExactLinearSolidError::Cancelled
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        matches!(error, ExactLinearSolidError::Cancelled)
    }
    fn invalid(what: &'static str) -> Self::Error {
        super::super::invalid(what)
    }
    fn workspace_bytes(&self) -> u128 {
        8 * (16 * self.free_dofs.len() as u128 + 4)
    }
    fn response_multiply_terms(&self) -> u128 {
        let n = self.free_dofs.len() as u128;
        32 * n * n + 4 * (3 * n + 1).pow(2)
    }
    fn workspace(&self) -> Result<Workspace, Self::Error> {
        let n = self.free_dofs.len();
        Ok(Workspace {
            lifted: vec![0.0; 3 * n + 1],
            balanced: vec![0.0; 2 * n],
            candidate: vec![0.0; 2 * n],
            integrals: vec![0.0; usize::from(self.viscous)],
        })
    }
    fn evaluate(
        &self,
        force: &[f64],
        scratch: &mut Workspace,
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        super::super::poll(cx)?;
        self.flow.evaluate_into(
            &self.state,
            force,
            &mut scratch.lifted,
            &mut scratch.balanced,
            &mut scratch.candidate,
            &mut scratch.integrals,
        )?;
        super::super::poll(cx)
    }
    fn displacement<'a>(&self, scratch: &'a Workspace) -> &'a [f64] {
        &scratch.candidate[..self.free_dofs.len()]
    }
    fn velocity<'a>(&self, scratch: &'a Workspace) -> &'a [f64] {
        &scratch.candidate[self.free_dofs.len()..]
    }
    fn work(
        &self,
        force: &[f64],
        scratch: &Workspace,
        _: &Cx<'_>,
    ) -> Result<(f64, f64), Self::Error> {
        self.held_force_work(force, &scratch.candidate, &scratch.integrals)
    }
    fn accept(
        &mut self,
        scratch: Workspace,
        force: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep, Self::Error> {
        let (work, loss) = self.held_force_work(force, &scratch.candidate, &scratch.integrals)?;
        let (step, time_s) = self.next_interval()?;
        self.accept(scratch.candidate, work, loss, step, time_s, cx)
    }
}
