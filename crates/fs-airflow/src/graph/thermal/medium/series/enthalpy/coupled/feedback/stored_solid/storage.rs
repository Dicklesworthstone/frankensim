//! The existing uniform and heterogeneous FEM owners behind one cooling state.
//!
//! This is dispatch, not another constitutive law or spatial solver. A chart
//! always belongs to its original vertex; the heterogeneous owner's private
//! uniform implementation is never exposed as a single-material model.

use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionProblem, ThermalInterfaces,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    heterogeneous::HeterogeneousEnthalpyBackwardEuler,
    radiation::EnthalpyRadiationStepSolution,
};
use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::EnthalpyBackwardEuler<'_, '_> {}
    impl Sealed for super::HeterogeneousEnthalpyBackwardEuler<'_, '_> {}
}

/// A prepared spatial enthalpy owner accepted by the shared cooling transaction.
///
/// Implemented only for the two existing conduction owners. They retain their
/// original mesh, invariant reference masses, phase/conductivity assignments,
/// solver controls and typed failures. Conductivity and interface bindings stay
/// on the caller's `StoredSolidProblem`; no material averaging is introduced.
pub trait SolidEnthalpyStorage: sealed::Sealed + core::fmt::Debug {
    /// Invariant reference mass per mesh vertex, kg.
    fn reference_nodal_masses_kg(&self) -> &[f64];

    /// Original chart for a mesh vertex, or `None` outside this storage model.
    fn phase_curve_at(&self, vertex: usize) -> Option<&EquilibriumEnthalpyPhaseCurve>;

    /// Delegate one physical interval to the existing enthalpy FEM owner.
    #[allow(clippy::too_many_arguments)]
    fn advance(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, EnthalpyError>;

    /// Delegate radiation to that SAME owner and immutable old enthalpy field.
    #[allow(clippy::too_many_arguments)]
    fn advance_with_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, patches: &[AmbientRadiationPatch],
        radiation: AmbientRadiationConfig,
    ) -> Result<EnthalpyRadiationStepSolution, EnthalpyError>;
}

impl SolidEnthalpyStorage for EnthalpyBackwardEuler<'_, '_> {
    fn reference_nodal_masses_kg(&self) -> &[f64] {
        self.reference_nodal_masses_kg()
    }

    fn phase_curve_at(&self, vertex: usize) -> Option<&EquilibriumEnthalpyPhaseCurve> {
        (vertex < self.reference_nodal_masses_kg().len()).then(|| self.phase_curve())
    }

    fn advance(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        EnthalpyBackwardEuler::advance(self, cx, problem, interfaces, old_h, dt_s, config)
    }

    fn advance_with_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, patches: &[AmbientRadiationPatch],
        radiation: AmbientRadiationConfig,
    ) -> Result<EnthalpyRadiationStepSolution, EnthalpyError> {
        EnthalpyBackwardEuler::advance_with_ambient_radiation(
            self, cx, problem, interfaces, old_h, dt_s, config, patches, radiation,
        )
    }
}

impl SolidEnthalpyStorage for HeterogeneousEnthalpyBackwardEuler<'_, '_> {
    fn reference_nodal_masses_kg(&self) -> &[f64] {
        self.reference_nodal_masses_kg()
    }

    fn phase_curve_at(&self, vertex: usize) -> Option<&EquilibriumEnthalpyPhaseCurve> {
        self.phase_curve_at(vertex)
    }

    fn advance(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        HeterogeneousEnthalpyBackwardEuler::advance(
            self, cx, problem, interfaces, old_h, dt_s, config,
        )
    }

    fn advance_with_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, patches: &[AmbientRadiationPatch],
        radiation: AmbientRadiationConfig,
    ) -> Result<EnthalpyRadiationStepSolution, EnthalpyError> {
        HeterogeneousEnthalpyBackwardEuler::advance_with_ambient_radiation(
            self, cx, problem, interfaces, old_h, dt_s, config, patches, radiation,
        )
    }
}
