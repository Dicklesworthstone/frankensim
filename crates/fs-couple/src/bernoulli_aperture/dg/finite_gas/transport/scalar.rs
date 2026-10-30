//! The original scalar moving aperture on the shared finite-gas transaction.
//!
//! This adapter changes only native source dispatch. Species, partial enthalpy,
//! endpoint work, driving-force checks and joint publication remain in the same
//! transport owner used by material walls and independent-blade apertures.

use super::{
    FiniteGasApertureDgStep, FiniteGasApertureTransportInput, FiniteGasApertureTransportStep,
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgTransportCandidate,
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    FiniteGasTransportNode, GasCompositionStep, GasTransportPaths, LoadedApertureDg,
    advance_finite_gas_transport_with_paths, finite_gas_transport_candidate_with_paths,
};
use crate::bernoulli_aperture::dg::ApertureDg;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use std::collections::BTreeMap;

/// Complete private endpoints of the scalar source and all finite inventories.
pub type FiniteGasScalarApertureTransportCandidate =
    FiniteGasDgTransportCandidate<ApertureDg, FiniteGasApertureDgStep<GasCompositionStep>>;

impl FiniteGasTransportNode for ApertureDg {
    type Step<S> = FiniteGasApertureDgStep<S>;
    type Input = FiniteGasApertureTransportInput;

    fn transport_input(
        input: &Self::Input,
        bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str> {
        <LoadedApertureDg as FiniteGasTransportNode>::transport_input(input, bound)
    }

    fn dt(&self) -> f64 {
        self.load().propagator().dt()
    }

    fn step_count(&self) -> u64 {
        self.duct_state().step_count()
    }

    fn stored_energy_j(&self) -> f64 {
        ApertureDg::<crate::bernoulli_aperture::dynamic::DynamicAperture>::stored_energy_j(self)
    }

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &Self::Input,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>> {
        let source = input.source.ok_or(FiniteGasDgError::Invalid(
            "scalar aperture transport is missing its original physical drive",
        ))?;
        if self.aperture().accepted_steps() != self.duct_state().step_count()
            || self.aperture().spec().time_step_s.to_bits() != self.dt().to_bits()
            || !source.drive.upstream_pressure_pa.is_finite()
            || !source.drive.body_flow_m3_s.is_finite()
            || !source.forces_n[0].is_finite()
            || source.forces_n[1] != 0.0
            || source.pressures_pa != [0.0; 2]
        {
            return Err(FiniteGasDgError::Invalid(
                "scalar aperture transport requires one clock and its finite scalar force; paired forces and modal pressure are not scalar inputs",
            ));
        }
        super::admit(
            self,
            gas,
            input.thermal.acoustic.ok_or(FiniteGasDgError::Invalid(
                "scalar aperture transport is missing acoustic work controls",
            ))?,
            cancelled,
        )
    }

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &Self::Input,
        molecular_heat_j: f64,
        evolve: impl FnMut(
            &ClosedGasParcel, f64, f64, &mut C,
        ) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<Self::Step<S>, FiniteGasDgError<E>> {
        if molecular_heat_j != 0.0 {
            return Err(FiniteGasDgError::Invalid(
                "scalar apertures do not own a molecular heat integral",
            ));
        }
        let source = input.source.expect("admitted scalar source drive");
        self.step_with_finite_gas_controlled(
            gas, source.drive, source.forces_n[0],
            input.thermal.acoustic.expect("admitted acoustic controls"),
            evolve, cancelled,
        )
    }

    fn accounts<S>(step: &Self::Step<S>) -> [f64; 4] {
        <LoadedApertureDg as FiniteGasTransportNode>::accounts(step)
    }

    fn source<S>(step: &Self::Step<S>) -> &S {
        &step.gas.background
    }
}

/// Advance the original scalar mechanics and finite reservoir network together.
/// All path laws, fixed-volume/first-order limits and rollback semantics are
/// those of the existing moving-aperture transport owner. Acoustic flow never
/// silently becomes thermodynamic mass transport.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_scalar_aperture_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, ApertureDg>,
    inputs: &BTreeMap<u64, FiniteGasApertureTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasApertureTransportStep, FiniteGasDgTransportError> {
    advance_finite_gas_transport_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}

/// Preview the same complete native transaction without advancing any owner.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_scalar_aperture_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, ApertureDg>,
    inputs: &BTreeMap<u64, FiniteGasApertureTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasScalarApertureTransportCandidate, FiniteGasDgTransportError> {
    finite_gas_transport_candidate_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}
