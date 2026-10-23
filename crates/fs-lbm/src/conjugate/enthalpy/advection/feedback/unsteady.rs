//! Material feedback on the retained native unsteady momentum owner.

use super::{
    DepositedCell, DepositionError, DepositionObserver, DepositionSimulation,
    FlowTransportSnapshot, FluidProperties, MaterialError, MaterialFlowError, MaterialFlowFeedback,
    MaterialMomentum, check_gate, invalid,
};
use crate::conjugate::{FlowStepRecord, FvFlow};
use fs_conduction::transient::deposition::CellId;
use std::collections::BTreeMap;

impl<F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'_, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Retain native face momentum instead of equilibrating it at every epoch.
    ///
    /// The closed laminar domain starts at rest at the bound thermal time.
    /// Every accepted substep first transports total h with the accepted old
    /// velocity, then takes exactly one native backward-Euler momentum step over
    /// that actual (possibly CFL-reduced) duration. Both stages use the same
    /// starting material epoch: local viscosity, Boussinesq forces and declared
    /// h-dependent mushy drag. Conductivity and caloric/latent storage keep their
    /// existing owners. This first-order split requires timestep refinement.
    ///
    /// Native face velocities and pressure survive subsequent advances. The
    /// thermal owner's clone-and-publish transaction also owns this history, so
    /// refusal/cancellation of either stage publishes neither. A source callback
    /// must not mutate external shared history. Resume with the returned thermal
    /// state and observer together; mixing clocks or enthalpy states refuses.
    ///
    /// Only stationary-normal walls/symmetry and no fans are admitted. Moving
    /// tangential walls and source buoyancy can supply mechanical work, but this
    /// mode does not claim a closed thermal-plus-kinetic energy account. Molecular,
    /// Darcy and pressure work are not deposited as heat. The steady creeping
    /// heat APIs therefore refuse in either attachment order. Fixed mass/density,
    /// phase mobility and free-surface restrictions of this owner still apply.
    /// Exact-solid continuation from `try_new_phase` refuses: changing native
    /// occupancy needs a solid momentum/kinetic-energy handoff. The original
    /// mobile-only mushy closure remains supported on its fixed native domain.
    /// SIMPLEC's existing max_iterations/tolerance bound each momentum step.
    pub fn with_unsteady_momentum(mut self) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if self.exact_phase {
            return Err(invalid(
                "retained momentum cannot change exact-solid phase occupancy without a solid momentum/kinetic-energy handoff; use quasisteady phase flow or a fixed-domain mobile closure",
            ).into());
        }
        if self.unsteady.is_some() || self.applied.is_some() {
            return Err(invalid(
                "attach retained momentum before the first advance; an accepted flow history cannot be replaced by rest",
            ).into());
        }
        if self.viscous_heating.is_some() || self.external_force_n_m3.is_some() {
            return Err(invalid(
                "steady viscous/body-force heat omits kinetic storage and cannot be combined with retained momentum",
            ).into());
        }
        let fluid = self.reference_fluid()?;
        self.unsteady = Some(MaterialMomentum::at_rest(
            self.domain,
            &fluid,
            self.config,
            self.bound_time_s,
            self.gate,
        )?);
        Ok(self)
    }

    /// Accepted endpoint momentum, distinct from applied().flow used by heat
    /// advection at the start. Before the first step this is the rest state with
    /// zero solver iterations, not a steady solution of the driven domain.
    pub fn unsteady_flow(&self) -> Option<&FvFlow> {
        self.unsteady.as_ref().map(|owner| &owner.flow)
    }

    /// Exact accepted thermal clock of the retained native momentum state.
    pub fn unsteady_time_s(&self) -> Option<f64> {
        self.unsteady.as_ref().map(|owner| owner.time_s)
    }

    /// Last accepted native transient solve, including diagnostic kinetic energy.
    /// None before the first step. This receipt is not an additional heat input.
    pub fn unsteady_step(&self) -> Option<&FlowStepRecord> {
        self.unsteady
            .as_ref()
            .and_then(|owner| owner.record.as_ref())
    }

    fn reference_fluid(&self) -> Result<FluidProperties, MaterialFlowError<O::Error>> {
        let reference = self.reference.ok_or(invalid(
            "retained momentum requires the admitted mobile reference; no solid transport is fabricated",
        ))?;
        Ok(FluidProperties {
            density_kg_m3: self.density_kg_m3,
            specific_heat_j_kg_k: reference.specific_heat_j_kg_k(),
            conductivity_w_m_k: reference.thermal_conductivity_w_m_k(),
            kinematic_viscosity_m2_s: reference.kinematic_viscosity_m2_s(),
        })
    }

    pub(super) fn check_unsteady_start(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        let Some(owner) = &self.unsteady else {
            return Ok(());
        };
        owner.check_time(thermal.time_s())?;
        for sample in &self.endpoint {
            check_gate(self.gate, check)?;
            let cell = thermal.cells().get(&sample.cell).ok_or(invalid(
                "retained momentum requires its accepted thermal cells",
            ))?;
            if cell.specific_h_j_kg.to_bits() != sample.specific_enthalpy_j_kg.to_bits() {
                return Err(invalid(
                    "retained momentum must resume with its accepted enthalpy state",
                )
                .into());
            }
        }
        Ok(())
    }

    pub(super) fn advance_unsteady(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        if self.unsteady.is_none() {
            return Ok(());
        }
        check_gate(self.gate, check)?;
        let fluid = self.reference_fluid()?;
        let applied = self.applied.as_ref().ok_or(invalid(
            "retained momentum needs the applied material epoch",
        ))?;
        if applied.time_s.to_bits() != start_s.to_bits() {
            return Err(
                invalid("retained momentum source epoch differs from the thermal start").into(),
            );
        }
        let reference = self.reference.ok_or(invalid(
            "retained momentum lacks its admitted mobile reference",
        ))?;
        let mut viscosity = vec![reference.dynamic_viscosity_pa_s(); self.domain.cell_count()];
        for sample in &applied.samples {
            check_gate(self.gate, check)?;
            if old.get(&sample.cell).map(|t| t.to_bits()) != Some(sample.temperature_k.to_bits()) {
                return Err(invalid(
                    "retained momentum must use the actual starting material samples",
                )
                .into());
            }
            let transport = sample.transport.mobile_transport().ok_or(invalid(
                "retained momentum cannot restore motion onto an exact solid constraint",
            ))?;
            viscosity[sample.voxel] = transport.dynamic_viscosity_pa_s();
        }
        self.unsteady
            .as_mut()
            .expect("checked retained momentum")
            .advance(
                self.domain,
                &fluid,
                self.config,
                &viscosity,
                applied.buoyancy_force_n_m3.as_deref(),
                applied.mushy_drag_pa_s_m2.as_deref(),
                start_s,
                dt_s,
                thermal.time_s(),
                self.gate,
            )?;
        check_gate(self.gate, check)?;
        Ok(())
    }
}
