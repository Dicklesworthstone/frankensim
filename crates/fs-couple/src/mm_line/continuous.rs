//! Causal coupled reflection from an original smooth lossless duct.
//!
//! The spatial source is fs-duct's full finite-mode Galerkin image of its
//! continuous Riccati operator. Every state and every inlet mode enters ONE
//! shared exact multiport scattering flow. For each pressure wave a held
//! throughout an interval, U=(2a-P)/Z0 and b=P-a. No fitted damping, scalar
//! substitution, discarded higher modes, extra latency or midpoint map.
//!
//! Warm push uses only prepared buffers. It returns interval MEAN b;
//! independently integrated b(t)^2/Z0 remains in the port accounts. The
//! square of the mean is not the reflected energy. This linear ideal-end,
//! lossless spatial rung does not qualify a nonlinear junction, real exterior
//! radiation or a spatial/clock accuracy band without separate refinement.

use fs_duct::RadiusProfile;
use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::{
    ModalEnd, ModalGalerkinSource, ModalTimeError, assemble_lossless,
};
use fs_material::gas::GasState;
use fs_phs::linear_exact::{
    LinearSplitBudget, MultiportScatteringStepRecord, PreparedMultiportScatteringZoh,
    ScatteringPortRecord,
};

/// Prepared full-state modal region and its accepted physical history.
///
/// Output pressures and waves use the original mode normalization. The
/// common positive real reference is rho*c/S at the original inlet; this
/// is a wave-coordinate choice, not a matched evanescent characteristic.
pub struct ContinuousModalLine {
    source: ModalGalerkinSource,
    flow: PreparedMultiportScatteringZoh,
    state: Vec<f64>,
    next: Vec<f64>,
    mean_pressures: Vec<f64>,
    reflected_means: Vec<f64>,
    ports: Vec<ScatteringPortRecord>,
}

impl ContinuousModalLine {
    /// Prepare original geometry, every retained coordinate and all inlet
    /// modes. The shared exact owner admits the complete dense work; no
    /// resolution is reduced to fit. Cancellation is polled during spatial
    /// assembly and after bounded exact-flow preparation.
    ///
    /// # Errors
    /// Original spatial, state-capacity, cancellation and exact-flow refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        profile: &RadiusProfile,
        gas: &GasState,
        modes: ModeCount,
        elements: usize,
        end: ModalEnd,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless(
            profile,
            gas,
            modes,
            elements,
            end,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, system) = image.into_parts();
        let states = 2 * source.coordinates();
        let count = source.modes().0;
        let reference = source.inlet_reference_impedance();
        let selected: Vec<_> = (0..count).map(|mode| (mode, reference)).collect();
        let flow = PreparedMultiportScatteringZoh::new_owned(system, &selected, duration, budget)?;
        if cancelled() {
            return Err(fs_duct::DuctError::Cancelled {
                completed_steps: elements,
            }
            .into());
        }
        Ok(Self {
            source,
            flow,
            state: vec![0.0; states],
            next: vec![0.0; states],
            mean_pressures: vec![0.0; count],
            reflected_means: vec![0.0; count],
            ports: vec![ScatteringPortRecord::default(); count],
        })
    }

    /// Original continuous source and the actual full spatial operators.
    pub const fn source(&self) -> &ModalGalerkinSource {
        &self.source
    }
    /// Accepted canonical physical state [potential; M*pressure].
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    /// Actual exact-flow interval [s].
    pub fn duration(&self) -> f64 {
        self.flow.duration()
    }
    /// Stored physical acoustic energy [J], with original spatial quadratics.
    pub fn stored_energy_j(&self) -> f64 {
        self.flow.system().hamiltonian(&self.state)
    }
    /// Accepted interval-mean reflected pressure waves [Pa], in mode order.
    pub fn incoming(&self) -> &[f64] {
        &self.reflected_means
    }
    /// Accepted interval-mean inlet modal pressure [Pa].
    pub fn mean_pressures(&self) -> &[f64] {
        &self.mean_pressures
    }
    /// Independently integrated incident/outgoing energy and boundary work.
    /// Zero entries before the first accepted push represent no interval.
    pub fn port_accounts(&self) -> &[ScatteringPortRecord] {
        &self.ports
    }

    /// Same-interval affine relation b_bar=free+direct*a, using accepted
    /// physical history. The caller supplies N and N*N output buffers.
    /// There is no extra sample of delay or history commit during this query.
    ///
    /// # Errors
    /// Shared exact owner's shape/nonfinite/refusal; buffers stay unchanged.
    pub fn reflection_relation_into(
        &self,
        free: &mut [f64],
        direct: &mut [f64],
    ) -> Result<(), ModalTimeError> {
        // All original inputs are selected waves, so their supplied values
        // are ignored by the shared affine relation. Reuse an existing N slice.
        self.flow
            .average_relation_into(&self.state, &self.reflected_means, free, direct)?;
        let count = self.source.modes().0;
        for mode in 0..count {
            direct[mode * count + mode] -= 1.0;
        }
        Ok(())
    }

    /// Consume one held incident pressure wave per ORIGINAL inlet mode.
    /// Commit physical history only after the complete coupled candidate
    /// and every independent energy account succeed. No warm allocation.
    ///
    /// # Errors
    /// Shared shape, finite-arithmetic, exact-flow or physical-work refusal.
    /// On any refusal accepted state, pressure, reflection and accounts remain
    /// unchanged; prepared numerical scratch is not physical history.
    pub fn push(
        &mut self,
        incident: &[f64],
    ) -> Result<MultiportScatteringStepRecord, ModalTimeError> {
        let record = self.flow.step_into(
            &self.state,
            incident,
            &mut self.next,
            &mut self.mean_pressures,
            &mut self.ports,
        )?;
        core::mem::swap(&mut self.state, &mut self.next);
        for (wave, port) in self.reflected_means.iter_mut().zip(&self.ports) {
            *wave = port.outgoing_wave;
        }
        Ok(record)
    }
}
