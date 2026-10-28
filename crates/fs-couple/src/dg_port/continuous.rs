//! Continuous single-duct composition keeps the original source and transaction.
use super::*;
use fs_duct::dg::continuous::{DgContinuous, state::DgTrajectoryPort};

impl DgPortLoad {
    pub(crate) fn continuous_single_duct(&self) -> Result<DgContinuous<'_>, AcousticRealizeError> {
        if self.terminal.is_some()
            || !self.branches.is_empty()
            || !self.terminal_channels.is_empty()
            || self.has_driven_mouth()
            || self.duct.continuous_operator().ports().len() != 1
        {
            return Err(invalid(
                "continuous DG composition requires one original duct and inlet, with no rational terminal, branch, gap or driven mouth",
            ));
        }
        Ok(self.duct.continuous_operator())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn preview_continuous_trajectory(
        &mut self,
        endpoint: Vec<f64>,
        dt: f64,
        pressure_impulse: f64,
        inlet_work_j: f64,
        dissipated_j: f64,
        absolute_tolerance_j: f64,
        relative_tolerance: f64,
    ) -> Result<DgPortTrial<'_>, AcousticRealizeError> {
        let original = self.continuous_single_duct()?;
        let (state, output) = original
            .candidate_from_trajectory(
                &self.state,
                endpoint,
                dt,
                &[DgTrajectoryPort {
                    id: self.inlet,
                    pressure_impulse_pa_s: pressure_impulse,
                    work_into_duct_j: inlet_work_j,
                }],
                dissipated_j,
                absolute_tolerance_j,
                relative_tolerance,
            )
            .map_err(AcousticRealizeError::Duct)?;
        Ok(DgPortTrial {
            owner: self,
            state,
            output,
            terminal: None,
            terminal_frame: None,
            terminal_frames: Vec::new(),
            inlet_work_j,
            source_work_j: inlet_work_j,
            driven_mouth_frame: None,
            branch_states: Vec::new(),
            branch_frames: Vec::new(),
            terminal_channels: Vec::new(),
            terminal_channel_frames: Vec::new(),
        })
    }
}
