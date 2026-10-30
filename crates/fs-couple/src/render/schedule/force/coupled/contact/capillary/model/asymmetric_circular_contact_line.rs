//! Unlike circular wetting surfaces through the original native mechanics.
use super::super::{CapillaryModalError, CapillaryProposal, CapillaryStepView};
use super::{NativeCapillaryModel, sealed};
use fs_contact::capillary_bridge::axisymmetric::asymmetric_contact_line::{
    AsymmetricCircularContactLineLoadLine, AsymmetricCircularContactLineState,
    AsymmetricCircularContactLineStep, SourceAsymmetricCircularContactLineBridge,
};
use fs_contact::capillary_bridge::contact_line::load_line::ContactLineLoadLineConfig;
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;

impl sealed::Sealed for AsymmetricCircularContactLineLoadLine<'_> {}
impl<'a> NativeCapillaryModel<'a> for AsymmetricCircularContactLineLoadLine<'a> {
    type Source = SourceAsymmetricCircularContactLineBridge;
    type State = AsymmetricCircularContactLineState;
    type Step = AsymmetricCircularContactLineStep;
    type Config = ContactLineLoadLineConfig;

    fn maximum_evaluations(config: Self::Config) -> usize {
        config.mechanical.max_evaluations
    }
    fn create(
        source: &'a Self::Source,
        compliance_m_per_n: f64,
        mut config: Self::Config,
    ) -> Result<Self, CapillaryModalError> {
        // Reserve the final original-source replay on actual body motion.
        config.mechanical.max_evaluations -= 1;
        Ok(Self::new(source, compliance_m_per_n, config)?)
    }
    fn substrate_materials(source: &Self::Source) -> [&MaterialStateId; 2] {
        [
            &source.original().wetting()[0].surfaces()[0].material,
            &source.original().wetting()[1].surfaces()[0].material,
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
        self.config().mechanical.maximum_gap_increment_m
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
            estimated_compatibility_error_m: result.estimated_compatibility_error_m,
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
        Ok(self.advance_actual(previous, actual_gap_m, dt_s, gate)?)
    }
    fn view(step: &Self::Step) -> CapillaryStepView<Self::State> {
        CapillaryStepView {
            state: step.state(),
            gap_m: step.state().gap_m(),
            surface_energy_j: step.profile.surface_free_energy_j,
            mean_force_n: step.mean_normal_force_n,
            contact_line_dissipation_j: step.contact_line_dissipation_j,
            topology_release_j: 0.0,
            estimated_energy_error_j: step.estimated_energy_error_j,
        }
    }
    fn force_error(
        &self,
        _previous: Self::State,
        step: &Self::Step,
    ) -> Result<f64, CapillaryModalError> {
        Ok(step.estimated_force_error_n)
    }
}
