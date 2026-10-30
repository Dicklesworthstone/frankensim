//! One finite source gas and the original loaded thermoplastic solid.
//!
//! The shared gas/enthalpy Newton residual calls the SAME native tetrahedral
//! return-map/energy proposal at every endpoint. Plastic and viscous heat,
//! reversible exchange, original hardening storage and numerical loss retain
//! their separate meanings. Accepted displacement and history seed every trial;
//! neither Newton probes nor thermal proposals advance the constitutive clock.
//!
//! Gas mass, volume and composition and solid reference geometry remain fixed.
//! This is quasistatic small-strain mechanics over the original fully solid
//! source domain. Gas pressure traction, moving gas boundaries and chemistry
//! are not inferred from this thermal connection.

use super::*;
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkSample, PreparedHeatNetwork};
use crate::circuit_enthalpy::material::{MaterialNetworkSample, PreparedMaterialNetwork};
use crate::solid_enthalpy::gas::{
    GasHeatTransport, SolidGasContact, SolidGasError, SolidGasReport, SolidGasSolve,
    StructuralEndpoint, StructuralExchange, StructuralSource, step_structural_gas,
};
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;

/// The original gas/thermal audit plus the native mechanical subaccounts.
#[derive(Clone, Debug)]
pub struct ThermoplasticGasReport<S = MaterialNetworkSample> {
    /// Its combined residual includes BOTH caloric inventories, mechanical
    /// internal energy, numerical loss and actual boundary mechanical work.
    pub thermal_gas: SolidGasReport<S>,
    /// Original endpoint-coefficient J2 accounts, before coefficient correction.
    pub element_energy: Vec<EnergyIncrement>,
    pub physical_heat_j: f64,
    pub reversible_heat_j: f64,
    pub coefficient_storage_change_j: f64,
    pub thermoelastic_storage_change_j: f64,
    /// Actual free-energy change, including coefficient work exactly once.
    pub stored_mechanical_change_j: f64,
    /// Constitutive backward-Euler loss; never deposited as heat.
    pub numerical_loss_j: f64,
    pub mechanical_work_j: f64,
    pub mechanical_residual_j: f64,
}

/// One atomic candidate containing the actual gas, enthalpies and plastic history.
pub struct ThermoplasticGasCandidate<'m, 't, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    gas: IsochoricGasParcel,
    solid: ThermoplasticSolid<'m>,
    report: ThermoplasticGasReport<S>,
}
impl<'m, 't, S> ThermoplasticGasCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    pub const fn solid(&self) -> &ThermoplasticSolid<'m> {
        &self.solid
    }
    pub const fn report(&self) -> &ThermoplasticGasReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        IsochoricGasParcel,
        ThermoplasticSolid<'m>,
        ThermoplasticGasReport<S>,
    ) {
        (self.thermal, self.gas, self.solid, self.report)
    }
}

impl<'m> ThermoplasticSolid<'m> {
    /// Solve gas heat, every thermal cell and loaded plastic equilibrium at
    /// the SAME endpoint, preserving each original accepted history.
    ///
    /// Source viscosity uses this actual physical interval when overstress
    /// is selected. Material moduli, yield, expansion and all
    /// interface transport resolve at the actual endpoint temperatures.
    /// Optional prepared internal conduction uses the original material owner;
    /// electrical work bindings refuse because this entry point advances no circuit.
    ///
    /// Guesses are sorted CellId energy increments followed by total gas heat
    /// in joules. The original gas controls bound trials, Newton work, source
    /// excursions and each independent energy balance. Their total-energy
    /// tolerance also bounds the native mechanical work account. The separate
    /// yield tolerance is in Pa. All failures leave the three inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_isochoric_gas<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        yield_tolerance_pa: f64,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticGasCandidate<'m, 't>, SolidGasError> {
        self.advance_gas_transport(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            prescribed_m,
            nodal_forces_n,
            settings,
            yield_tolerance_pa,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    /// The same coupled endpoint with the complete prepared heat assembly.
    ///
    /// Original thermal contacts, integrated conduction layers and closed
    /// radiation all feed the SAME cell/gas/mechanical residual. Their complete
    /// source sample remains in the accepted report. Gas patches borrow its
    /// original material samples; contact and radiation heat are never repeated
    /// as external supply. Electrical bindings refuse as on the bulk entry point.
    /// All interval, initial-guess and atomicity rules of
    /// [`Self::advance_with_isochoric_gas`] also apply.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_isochoric_gas_heat_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: &PreparedHeatNetwork<'_, '_>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        yield_tolerance_pa: f64,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticGasCandidate<'m, 't, HeatNetworkSample>, SolidGasError> {
        self.advance_gas_transport(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            prescribed_m,
            nodal_forces_n,
            settings,
            yield_tolerance_pa,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_gas_transport<'t, T: GasHeatTransport>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: T,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        yield_tolerance_pa: f64,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticGasCandidate<'m, 't, T::Sample>, SolidGasError> {
        poll::<Infallible>(cx)?;
        let duration_s = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !positive(duration_s)
            || !target_s.is_finite()
            || !positive(yield_tolerance_pa)
            || !positive(controls.total_energy_tolerance_j)
        {
            return Err(SolidGasError::Invalid(
                "matching accepted clocks, positive interval and mechanical energy/yield tolerances required",
            ));
        }
        self.enthalpies::<Infallible>(thermal, true, cx)?;
        for (i, material) in self.materials.iter().enumerate() {
            poll::<Infallible>(cx)?;
            let bound = &self.cells[i];
            let phase = bound
                .source
                .caloric_source()
                .curve()
                .state_at_specific_enthalpy(bound.accepted.specific_h_j_kg)
                .map_err(super::super::SolidEnthalpyError::from)?;
            let strain = bound
                .source
                .mechanical_strain(phase.temperature_k(), &self.solution.element_strain[i])
                .map_err(ThermoplasticError::<Infallible>::Material)?;
            bound
                .source
                .validate_committed_history(
                    material,
                    &strain,
                    &self.solution.element_states[i],
                    yield_tolerance_pa,
                )
                .map_err(ThermoplasticError::<Infallible>::Energy)?;
        }
        let sources: Vec<_> = self
            .cells
            .iter()
            .map(|bound| StructuralSource {
                cell: bound.accepted.id,
                curve: bound.source.caloric_source().curve(),
                point: bound.source.thermal_material_point(),
                reference_h_j_kg: bound.reference_h_j_kg,
            })
            .collect();
        // Only the native mechanical energy/yield fields are used here. The
        // common gas Newton owner supplies all thermal nonlinear iteration.
        let mechanical_controls = ThermoplasticConfig {
            maximum_iterations: 1,
            temperature_tolerance_k: f64::EPSILON,
            energy_tolerance_j: controls.total_energy_tolerance_j,
            yield_tolerance_pa,
            enthalpy_relaxation: 1.0,
        };
        let candidate = step_structural_gas(
            &sources,
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
            |endpoint| {
                let enthalpies = self.enthalpies::<Infallible>(endpoint, false, cx)?;
                let mechanical = self.mechanics_with_motion::<Infallible>(
                    &enthalpies,
                    prescribed_m,
                    nodal_forces_n,
                    settings,
                    mechanical_controls,
                    duration_s,
                    None,
                    cx,
                )?;
                let exchange = StructuralExchange {
                    heat_j: mechanical.heat_j.clone(),
                    stored_change_j: mechanical.stored_change_j
                        + mechanical.thermoelastic_storage_change_j,
                    numerical_loss_j: mechanical.numerical_loss_j,
                    external_work_j: mechanical.external_work_j,
                };
                Ok(StructuralEndpoint {
                    state: mechanical,
                    exchange: Some(exchange),
                })
            },
        )?;
        let mut gas_reports = candidate.report.gases;
        let mut gases = candidate.gases;
        let gas = gases
            .pop()
            .ok_or(SolidGasError::Invalid("missing gas endpoint"))?
            .gas;
        let single = gas_reports
            .pop()
            .ok_or(SolidGasError::Invalid("missing gas account"))?;
        let mechanical = candidate.state;
        let thermal = candidate.thermal;
        let mut solid = self.clone();
        solid.time_s = target_s;
        solid.materials = mechanical.materials;
        solid.solution = mechanical.solution;
        for cell in &mut solid.cells {
            cell.accepted = thermal.cells()[&cell.accepted.id].clone();
        }
        let report = ThermoplasticGasReport {
            thermal_gas: SolidGasReport {
                thermal: candidate.report.thermal,
                gas: single.gas,
                contacts: single.contacts,
                transport: candidate.report.transport,
                closed_radiation_net_j: candidate.report.closed_radiation_net_j,
                external_thermal_heat_j: candidate.report.external_thermal_heat_j,
                external_gas_heat_j: single.external_heat_j,
                solid_to_gas_heat_j: single.solid_to_gas_heat_j,
                gas_energy_residual_j: single.energy_residual_j,
                combined_energy_residual_j: candidate.report.combined_energy_residual_j,
                maximum_cell_residual_j: candidate.report.maximum_cell_residual_j,
                nonlinear: candidate.report.nonlinear,
                trials: candidate.report.trials,
            },
            element_energy: mechanical.energy,
            physical_heat_j: mechanical.physical_heat_j,
            reversible_heat_j: mechanical.reversible_heat_j,
            coefficient_storage_change_j: mechanical.coefficient_storage_change_j,
            thermoelastic_storage_change_j: mechanical.thermoelastic_storage_change_j,
            stored_mechanical_change_j: mechanical.stored_change_j,
            numerical_loss_j: mechanical.numerical_loss_j,
            mechanical_work_j: mechanical.external_work_j,
            mechanical_residual_j: mechanical.mechanical_residual_j,
        };
        poll::<Infallible>(cx)?;
        Ok(ThermoplasticGasCandidate {
            thermal,
            gas,
            solid,
            report,
        })
    }
}
