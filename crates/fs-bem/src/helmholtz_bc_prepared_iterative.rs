//! Reusable Robin source solves through the original FMM/GMRES owner.
//! Preparation is shared; every RHS retains its own admitted original-equation
//! residual and resource receipt. No full boundary matrix or alternate kernel
//! is built; the original bounded near-field preconditioner blocks remain.
use super::*;
use fs_solver::krylov::SolveReport;

pub struct PreparedIterativeBoundary<'a> {
    operator: IterativeRadiation<'a>,
    surface: &'a SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'a>,
    law_stamp: u64,
    maximum_iterations: usize,
}

impl<'a> PreparedIterativeBoundary<'a> {
    pub fn new(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        law: Admittance<'a>,
        formulation: Formulation,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        let operator =
            IterativeRadiation::prepare_bc(surface, k, medium, formulation, options, law, gate)?;
        Ok(Self {
            operator,
            surface,
            k,
            medium,
            law,
            law_stamp: law.stamp(),
            maximum_iterations: options.max_iters,
        })
    }

    /// One prescribed boundary drive in m/s. Original Robin work accounting
    /// and convergence refusal are retained. The receipt includes the shared
    /// setup time, so receipts must not be summed as disjoint setup costs.
    pub fn solve(
        &self,
        drive_velocity: &[C64],
        gate: &CancelGate,
    ) -> Result<(BoundarySolution, SolveReport, ResourceReceipt), IterativeError> {
        if gate.is_requested() {
            return Err(IterativeError::Cancelled { matvecs: 0 });
        }
        self.check_law()?;
        let mut state = self.operator.start(drive_velocity, gate)?;
        self.operator
            .advance(&mut state, self.maximum_iterations.saturating_add(1), gate)?;
        let (solution, report, receipt) = self.operator.finish(&state)?;
        if gate.is_requested() {
            return Err(IterativeError::Cancelled {
                matvecs: receipt.matvecs,
            });
        }
        let bc = SurfaceBc {
            drive_velocity,
            admittance: self.law,
        };
        let solution = finish_boundary(self.surface, solution, bc)?;
        self.check_law()?;
        Ok((solution, report, receipt))
    }

    /// Incident pressure AND normal velocity enter the same Robin drive.
    /// Absorbed work uses the total field through the original finalizer.
    pub fn scatter(
        &self,
        incident: &IncidentTraces,
        gate: &CancelGate,
    ) -> Result<(ScatteringSolution, SolveReport, ResourceReceipt), IterativeError> {
        if gate.is_requested() {
            return Err(IterativeError::Cancelled { matvecs: 0 });
        }
        incident.validate_binding(self.surface, self.k, self.medium)?;
        let drive = incident.drive(self.law)?;
        let (solution, report, receipt) = self.solve(&drive, gate)?;
        let solution = finish_scattering(
            self.surface,
            solution,
            incident.pressure_pa().to_vec(),
            incident.normal_velocity_m_s().to_vec(),
            self.law,
        )?;
        self.check_law()?;
        Ok((solution, report, receipt))
    }

    fn check_law(&self) -> Result<(), IterativeError> {
        if self.law.stamp() != self.law_stamp {
            return Err(
                bad("prepared iterative boundary admittance changed after preparation").into(),
            );
        }
        Ok(())
    }
}
