//! Existing dense exact owner behind the shared sustained-contact admission.

use super::*;

pub(crate) struct Workspace {
    lifted: Vec<f64>,
    balanced: Vec<f64>,
    candidate: Vec<f64>,
    integrals: Vec<f64>,
}

impl Native for ExactLinearSolid {
    type Owner = Self;
    type Step = ExactLinearSolidStep;
    type Error = ExactLinearSolidError;
    type Workspace = Workspace;

    fn owner(&self) -> &Self {
        self
    }
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
    fn work(&self, force: &[f64], scratch: &Workspace, _: &Cx<'_>) -> Result<Work, Self::Error> {
        self.held_force_work(force, &scratch.candidate, &scratch.integrals)
            .map(|(applied_j, viscosity_j)| Work {
                applied_j,
                viscosity_j,
                support_j: 0.0,
            })
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

impl events::EventNative for ExactLinearSolid {
    fn trace_response(
        &self,
        weights: [&[f64]; 3],
        cx: &Cx<'_>,
    ) -> Result<events::TraceResponse, Self::Error> {
        let n = self.free_dofs.len();
        let mut columns: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; n]);
        for a in 0..3 {
            for i in 0..n {
                super::super::poll(cx)?;
                columns[a][i] = finite(dot(&self.inverse_mass[i * n..(i + 1) * n], weights[a]))
                    .map_err(Self::invalid)?;
            }
        }
        Ok(events::TraceResponse {
            columns,
            modal_columns: None,
            support_reaction_columns: None,
            acceleration_offset: [0.0; 3],
        })
    }

    fn damping(&self) -> Option<&[f64]> {
        self.damping.as_deref()
    }

    fn add_contact_moments(
        &self,
        scratch: &mut Workspace,
        response: &events::TraceResponse,
        path: &events::ContactForcePath<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<f64, Self::Error> {
        let n = self.free_dofs.len();
        for i in 0..n {
            super::super::poll(cx)?;
            let mut q = Dd::from_f64(scratch.candidate[i]);
            let mut v = Dd::from_f64(scratch.candidate[n + i]);
            for a in 0..3 {
                let column = Dd::from_f64(response.columns[a][i]);
                q = q + column * path.displacement_moment[a];
                v = v + column * path.impulse[a];
            }
            scratch.candidate[i] = finite(q).map_err(Self::invalid)?;
            scratch.candidate[n + i] = finite(v).map_err(Self::invalid)?;
        }
        Ok(0.0)
    }

    fn accept_force_path(
        &mut self,
        scratch: Workspace,
        supplied_j: f64,
        viscosity_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep, Self::Error> {
        let (step, time_s) = self.next_interval()?;
        // Contact corrections lie in C's nullspace, so the original observed
        // viscosity integral is unchanged, including cross terms. Every native
        // original-coordinate and cumulative energy gate still runs here.
        self.accept(scratch.candidate, supplied_j, viscosity_j, step, time_s, cx)
    }
}
