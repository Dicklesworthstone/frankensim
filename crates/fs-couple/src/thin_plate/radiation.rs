//! Reciprocal exterior loading of the original nonlinear plate and its histories.
use super::{AcousticRealizeError, VkBody, prepare_render_linear, render_step_error};
use crate::render::plate::impact::radiation::{MAX_POLES, Model};
use fs_exec::CancelGate;

impl VkBody {
    /// Attach the complete shared passive radiation model in this modal basis.
    ///
    /// Model column k is the SAME mass-normalized velocity p_k as this body's
    /// original q/p pair. Every retained mode must have a column, including
    /// zero columns. Signed mutual terms act on all of these physical ports;
    /// the independent point-force input and surface observation keep their
    /// original projections. Neither model nor area is collapsed to that drive.
    ///
    /// The existing modal-bank owner supplies acoustic storage and damping, and
    /// the existing lossless interconnection supplies the reaction -L^T p_air.
    /// The complete material, nonlinear and acoustic dynamics then use the SAME
    /// exact/Gonzalez owner, including `into_prepared`. No split-force update or
    /// second integrator is introduced. The work record includes acoustic loss
    /// once; that loss must not also be deposited as solid material heat.
    ///
    /// Every accepted host coordinate and history stays unchanged. New acoustic
    /// q/p pairs follow the current state, in model pole order, initially at rest.
    /// Attaching during motion explicitly introduces a quiescent load at that
    /// instant; it does not infer its unobserved prehistory. A zero coupling is
    /// an exact no-op. Multiple calls add distinct declared loads; the caller
    /// must supply each physical exterior only once.
    ///
    /// The model producer owns geometry, material/gas state, modal-basis identity
    /// and fitted-band validity. Passivity of this composition is not evidence
    /// that an authored pole model represents a physical surface. The compact
    /// pressure readout remains the body's original observation; it is not a
    /// finite-wavelength receiver transfer supplied by this impedance.
    ///
    /// # Errors
    /// Invalid complete model/basis/state, exceeded explicit pole or original
    /// dense render capacity, cancellation, or refusal by the existing owners.
    pub fn with_radiation_load(
        self,
        model: &Model,
        maximum_poles: usize,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        let poll = || {
            if gate.is_requested() {
                Err(AcousticRealizeError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        model.validate().map_err(AcousticRealizeError::Nonlinear)?;
        let modes = self.n_modes();
        let host_states = self.sys.state_dim();
        if modes == 0
            || model.ports != modes
            || maximum_poles > MAX_POLES
            || model.poles.len() > maximum_poles
            || modes.checked_mul(2).is_none_or(|n| n > host_states)
            || self.x.len() != host_states
            || self.x.iter().chain(&self.areas).any(|v| !v.is_finite())
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "radiation needs every original mode, complete finite state and its declared pole capacity",
            });
        }
        if model
            .poles
            .iter()
            .all(|p| p.coupling.iter().all(|v| *v == 0.0))
        {
            poll()?;
            return Ok(self);
        }
        let extent = model
            .poles
            .len()
            .checked_mul(2)
            .and_then(|n| host_states.checked_add(n));
        let block = extent
            .and_then(|n| n.checked_add(self.sys.port_dim()))
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(2));
        if block.is_none_or(|n| n > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "complete acoustic, material and mechanical state exceeds the original dense render capacity",
            });
        }
        let omega: Vec<_> = model.poles.iter().map(|p| p.omega).collect();
        let zeta: Vec<_> = model.poles.iter().map(|p| p.zeta).collect();
        let drives: Vec<Vec<_>> = (0..modes)
            .map(|mode| model.poles.iter().map(|p| p.coupling[mode]).collect())
            .collect();
        let drive_refs: Vec<_> = drives.iter().map(Vec::as_slice).collect();
        poll()?;
        let load =
            fs_phs::modal_bank_ports(&omega, &zeta, &drive_refs).map_err(render_step_error)?;
        let mut columns = vec![0.0; host_states * modes];
        for mode in 0..modes {
            columns[(2 * mode + 1) * modes + mode] = 1.0;
        }
        poll()?;
        let sys = self
            .sys
            .with_parallel_load_columns(load, &columns)
            .map_err(render_step_error)?;
        let linear = prepare_render_linear(&sys)?;
        poll()?;
        let mut x = self.x;
        x.resize(sys.state_dim(), 0.0);
        Ok(Self {
            sys,
            linear,
            x,
            areas: self.areas,
        })
    }
}
