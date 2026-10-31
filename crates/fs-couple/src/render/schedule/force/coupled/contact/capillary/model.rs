//! Original constitutive adapters for the shared native mechanical transaction.
//! Geometry, source identity, work and numerical applicability remain with the
//! existing bridge owners. No replacement capillary law is introduced here.

use super::{CapillaryModalError, CapillaryProposal, CapillaryStepView};
use fs_contact::capillary_bridge::axisymmetric::load_line::{
    AxisymmetricBridgeLoadLine, AxisymmetricLoadLineConfig,
};
use fs_contact::capillary_bridge::axisymmetric::{
    AxisymmetricState, AxisymmetricStep, SourceAxisymmetricCapillaryBridge,
};
use fs_contact::capillary_bridge::load_line::{BridgeLoadLineConfig, PlanarBridgeLoadLine};
use fs_contact::capillary_bridge::{BridgeState, BridgeStep, SourcePlanarCapillaryBridge};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;

mod asymmetric_circular_contact_line;
mod asymmetric_contact_line;
mod axisymmetric_contact_line;
mod contact_line;

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::PlanarBridgeLoadLine<'_> {}
    impl Sealed for super::AxisymmetricBridgeLoadLine<'_> {}
}

/// The original capillary models admitted by the native mechanics owner.
/// This sealed contract keeps one mechanical transaction for planar bridges
/// with topology changes and axisymmetric bridges of fixed topology.
pub trait NativeCapillaryModel<'a>: sealed::Sealed + Sized {
    type Source: 'a;
    type State: Copy;
    type Step: Copy + core::fmt::Debug + PartialEq;
    type Config: Copy;

    #[doc(hidden)]
    fn maximum_evaluations(config: Self::Config) -> usize;
    #[doc(hidden)]
    fn create(
        source: &'a Self::Source,
        compliance_m_per_n: f64,
        config: Self::Config,
    ) -> Result<Self, CapillaryModalError>;
    #[doc(hidden)]
    fn substrate_materials(source: &Self::Source) -> [&MaterialStateId; 2];
    #[doc(hidden)]
    fn gap_energy(&self, state: Self::State) -> Result<(f64, f64), CapillaryModalError>;
    #[doc(hidden)]
    fn compliance(&self) -> f64;
    #[doc(hidden)]
    fn maximum_gap_increment(&self) -> f64;
    #[doc(hidden)]
    fn propose(
        &self,
        previous: Self::State,
        free_gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<CapillaryProposal<Self::Step>, CapillaryModalError>;
    #[doc(hidden)]
    fn advance_actual(
        &self,
        previous: Self::State,
        actual_gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<Self::Step, CapillaryModalError>;
    #[doc(hidden)]
    fn view(step: &Self::Step) -> CapillaryStepView<Self::State>;
    #[doc(hidden)]
    fn force_error(
        &self,
        previous: Self::State,
        step: &Self::Step,
    ) -> Result<f64, CapillaryModalError>;
}

impl<'a> NativeCapillaryModel<'a> for PlanarBridgeLoadLine<'a> {
    type Source = SourcePlanarCapillaryBridge;
    type State = BridgeState;
    type Step = BridgeStep;
    type Config = BridgeLoadLineConfig;

    fn maximum_evaluations(config: Self::Config) -> usize {
        config.max_evaluations
    }
    fn create(
        source: &'a Self::Source,
        compliance_m_per_n: f64,
        mut config: Self::Config,
    ) -> Result<Self, CapillaryModalError> {
        // The native constructor checks >=2 before reserving this evaluation.
        config.max_evaluations -= 1;
        Ok(Self::new(source, compliance_m_per_n, config)?)
    }
    fn substrate_materials(source: &Self::Source) -> [&MaterialStateId; 2] {
        [
            &source.wetting()[0].surfaces()[0].material,
            &source.wetting()[1].surfaces()[0].material,
        ]
    }
    fn gap_energy(&self, state: Self::State) -> Result<(f64, f64), CapillaryModalError> {
        let profile = self.source().evaluate(state)?;
        Ok((state.gap_m(), profile.surface_free_energy_j))
    }
    fn compliance(&self) -> f64 {
        self.compliance_m_per_n()
    }
    fn maximum_gap_increment(&self) -> f64 {
        self.config().maximum_gap_increment_m
    }
    fn propose(
        &self,
        previous: Self::State,
        free_gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<CapillaryProposal<Self::Step>, CapillaryModalError> {
        let result = self.trial(previous, free_gap_m, dt_s, gate)?;
        Ok(CapillaryProposal {
            bridge: result.bridge,
            compatibility_tolerance_m: result.compatibility_tolerance_m,
            estimated_compatibility_error_m: 0.0,
            iterations: result.iterations,
            evaluations: result.evaluations,
        })
    }
    fn advance_actual(
        &self,
        previous: Self::State,
        actual_gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<Self::Step, CapillaryModalError> {
        Ok(self.source().advance(previous, actual_gap_m, dt_s, gate)?)
    }
    fn view(step: &Self::Step) -> CapillaryStepView<Self::State> {
        CapillaryStepView {
            state: step.state(),
            gap_m: step.state().gap_m(),
            surface_energy_j: step.profile.surface_free_energy_j,
            mean_force_n: step.mean_normal_force_n,
            contact_line_dissipation_j: 0.0,
            topology_release_j: step.irreversible_release_j,
            estimated_energy_error_j: 0.0,
        }
    }
    fn force_error(
        &self,
        _previous: Self::State,
        _step: &Self::Step,
    ) -> Result<f64, CapillaryModalError> {
        // The planar path uses closed-form geometry/work, without quadrature.
        Ok(0.0)
    }
}

impl<'a> NativeCapillaryModel<'a> for AxisymmetricBridgeLoadLine<'a> {
    type Source = SourceAxisymmetricCapillaryBridge;
    type State = AxisymmetricState;
    type Step = AxisymmetricStep;
    type Config = AxisymmetricLoadLineConfig;

    fn maximum_evaluations(config: Self::Config) -> usize {
        config.mechanical.max_evaluations
    }
    fn create(
        source: &'a Self::Source,
        compliance_m_per_n: f64,
        mut config: Self::Config,
    ) -> Result<Self, CapillaryModalError> {
        config.mechanical.max_evaluations -= 1;
        Ok(Self::new(source, compliance_m_per_n, config)?)
    }
    fn substrate_materials(source: &Self::Source) -> [&MaterialStateId; 2] {
        [
            &source.wetting()[0].surfaces()[0].material,
            &source.wetting()[1].surfaces()[0].material,
        ]
    }
    fn gap_energy(&self, state: Self::State) -> Result<(f64, f64), CapillaryModalError> {
        let profile = self.source().evaluate(state)?;
        Ok((profile.gap_m, profile.surface_free_energy_j))
    }
    fn compliance(&self) -> f64 {
        self.compliance_m_per_n()
    }
    fn maximum_gap_increment(&self) -> f64 {
        self.config().mechanical.maximum_gap_increment_m
    }
    fn propose(
        &self,
        previous: Self::State,
        free_gap_m: f64,
        _dt_s: f64,
        gate: &CancelGate,
    ) -> Result<CapillaryProposal<Self::Step>, CapillaryModalError> {
        let result = self.trial(previous, free_gap_m, gate)?;
        Ok(CapillaryProposal {
            bridge: result.bridge,
            compatibility_tolerance_m: result.compatibility_tolerance_m,
            estimated_compatibility_error_m: result.estimated_compatibility_error_m,
            iterations: result.iterations,
            evaluations: result.evaluations,
        })
    }
    fn advance_actual(
        &self,
        previous: Self::State,
        actual_gap_m: f64,
        _dt_s: f64,
        gate: &CancelGate,
    ) -> Result<Self::Step, CapillaryModalError> {
        Ok(self.advance_actual(previous, actual_gap_m, gate)?)
    }
    fn view(step: &Self::Step) -> CapillaryStepView<Self::State> {
        let state = step.state();
        CapillaryStepView {
            state,
            gap_m: state.profile().gap_m,
            surface_energy_j: state.profile().surface_free_energy_j,
            mean_force_n: step.mean_normal_force_n,
            contact_line_dissipation_j: 0.0,
            // The original axisymmetric owner has no topology transition.
            topology_release_j: 0.0,
            estimated_energy_error_j: step.estimated_work_error_j
                + step.surface_energy_balance_residual_j.abs(),
        }
    }
    fn force_error(
        &self,
        previous: Self::State,
        step: &Self::Step,
    ) -> Result<f64, CapillaryModalError> {
        Ok(self.estimated_force_error_n(previous, step)?)
    }
}
