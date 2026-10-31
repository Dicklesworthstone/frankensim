//! Material feedback on the retained native unsteady momentum owner.

use super::{
    DepositedCell, DepositionError, DepositionObserver, DepositionSimulation,
    FlowTransportSnapshot, FluidProperties, MaterialError, MaterialFlowError, MaterialFlowFeedback,
    MaterialMomentum, check_gate, invalid,
};
use crate::conjugate::{
    FlowStepRecord, FvFlow, UnsteadyHeatingAudit, UnsteadyMechanicalStep, ViscousHeatControls,
};
use fs_conduction::transient::deposition::CellId;
use std::borrow::Cow;
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
    /// mode alone does not close thermal-plus-kinetic energy. Its molecular and
    /// Darcy losses are not deposited as heat; select
    /// Self::with_unsteady_viscous_heating for the explicit unsteady Stokes
    /// energy closure. Steady heat APIs still refuse in either attachment order.
    /// Fixed mass/density,
    /// phase mobility and free-surface restrictions of this owner still apply.
    /// Exact-solid phase cells retain their caloric/material identities. A
    /// melting face opens at its existing zero speed while every other native
    /// velocity survives. A newly frozen face must ALREADY have zero native
    /// speed, otherwise the step refuses an unavailable momentum/kinetic-energy
    /// handoff. No tolerance clipping or projection energy becomes heat.
    /// All-solid intervals carry the same clock without a momentum solve or
    /// fabricated liquid transport. This is a fixed stationary phase skeleton,
    /// not moving solid mechanics, shrinkage or stress-reference restoration.
    /// SIMPLEC's existing max_iterations/tolerance bound each momentum step.
    pub fn with_unsteady_momentum(mut self) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
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
        let active = self.momentum_domain(&self.endpoint, &mut || Ok(()))?;
        let domain = active.as_ref().unwrap_or(self.domain);
        self.check_phase_drives(domain, &mut || Ok(()))?;
        self.reference = self.reference.or_else(|| self.endpoint.iter()
            .find_map(|sample| sample.transport.mobile_transport()));
        self.unsteady = Some(if domain.fluid_count() == 0 {
            MaterialMomentum::at_rest_immobile(domain, self.bound_time_s, self.gate)?
        } else {
            let fluid = self.reference_fluid()?;
            MaterialMomentum::at_rest(
                domain,
                &fluid,
                self.config,
                self.bound_time_s,
                self.gate,
            )?
        });
        Ok(self)
    }

    /// Evolve native inertia and deposit physical molecular and porous heat.
    ///
    /// This explicitly selects unsteady Stokes: native backward-Euler inertia,
    /// full source-resolved stress and stationary Darcy/Forchheimer drag remain, but
    /// momentum convection is omitted under the declared Reynolds ceiling.
    /// Old accepted velocity still advects h. Once its combined CFL/event step
    /// is known, ONE native momentum solve supplies endpoint physical heat for
    /// that exact interval. No guessed dt or second thermal integrator is used.
    ///
    /// The native face mass norm stores kinetic energy. Backward-Euler loss
    /// has its own nonnegative numerical account and never heats material.
    /// Wall/body work and negligible pressure-constraint work close against
    /// kinetic change, physical heat and that numerical loss to the declared
    /// power allowance. Refine time and space for accuracy.
    ///
    /// Material Boussinesq force uses density from the same starting epoch and
    /// adds to any external force in the native momentum RHS. Its signed work
    /// is a separately reported subset of total body work, never a second heat
    /// source. Molecular, Darcy and explicitly sourced Forchheimer dissipation
    /// enter the thermal owner with separate physical-loss subaccounts.
    /// The reference-density inertia and fixed cell masses remain unchanged.
    ///
    /// Phase locking of a moving face, openings, fans, configured resistance
    /// and turbulence refuse. A stationary Darcy skeleton must share the
    /// mapped bulk enthalpy at local thermal equilibrium. This mode does not
    /// debit a finite actuator or gravitational-potential reservoir, nor model
    /// pressure compressibility, moving solid work or phase-density changes.
    ///
    /// Attach before the first momentum advance. Reattaching controls to an
    /// already enabled history at an accepted endpoint retains all momentum;
    /// it never restarts the fluid. Existing steady heating remains separate.
    pub fn with_unsteady_viscous_heating(
        mut self,
        controls: ViscousHeatControls,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if self.viscous_heating.is_some() {
            return Err(invalid(
                "unsteady mechanical heating cannot be combined with steady heating",
            ).into());
        }
        crate::conjugate::simple::viscous_heat::admit_heating(self.domain, self.config, controls)?;
        if self.unsteady.is_none() {
            self = self.with_unsteady_momentum()?;
        }
        self.unsteady
            .as_mut()
            .expect("attached native momentum")
            .enable_heating(controls)?;
        Ok(self)
    }

    /// Supply signed external actuator force to the same unsteady heat owner.
    ///
    /// One vector [N/m3] per native voxel; fixed solid entries must be zero.
    /// Work is evaluated on the SAME endpoint faces as inertia/dissipation;
    /// negative work is permitted. Only physical losses heat material.
    /// This force adds to material buoyancy. The mechanical receipt keeps
    /// total body work and its signed buoyancy subset; their difference is
    /// actuator work. Neither work account is an additional thermal source.
    /// Updating this drive between accepted advances preserves face momentum,
    /// so turning it off permits actual kinetic decay rather than a restart.
    pub fn with_unsteady_body_force_heating(
        self,
        force_density_n_m3: Vec<[f64; 3]>,
        controls: ViscousHeatControls,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if force_density_n_m3.len() != self.domain.cell_count() {
            return Err(invalid("unsteady heating needs one body-force vector per voxel").into());
        }
        for (c, force) in force_density_n_m3.iter().enumerate() {
            check_gate(self.gate, &mut || Ok(()))?;
            if force
                .iter()
                .any(|f| !f.is_finite() || (!self.domain.is_fluid(c) && *f != 0.0))
            {
                return Err(invalid(
                    "unsteady body forces must be finite and zero on fixed solids",
                )
                .into());
            }
        }
        let mut coupled = self.with_unsteady_viscous_heating(controls)?;
        coupled.external_force_n_m3 = Some(force_density_n_m3);
        let active = coupled.momentum_domain(&coupled.endpoint, &mut || Ok(()))?;
        coupled.check_phase_drives(active.as_ref().unwrap_or(coupled.domain), &mut || Ok(()))?;
        Ok(coupled)
    }

    /// Latest accepted exact-interval native kinetic/work/physical-loss receipt.
    pub fn unsteady_mechanical_step(&self) -> Option<&UnsteadyMechanicalStep> {
        self.unsteady
            .as_ref()
            .and_then(|owner| owner.mechanical_step.as_ref())
    }

    /// This accepted advance's supply, stored kinetic energy and physical heat.
    /// Generated heat is already counted by the original thermal EnergyAudit.
    pub fn unsteady_heating_audit(&self) -> Option<UnsteadyHeatingAudit> {
        self.unsteady
            .as_ref()
            .and_then(MaterialMomentum::heating_audit)
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

    pub(super) fn reference_fluid(&self) -> Result<FluidProperties, MaterialFlowError<O::Error>> {
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
        let Some(owner) = &mut self.unsteady else {
            return Ok(());
        };
        if owner.has_heating() {
            // Prepared before thermal integration with the accepted interval.
            // Do not integrate momentum a second time at the thermal endpoint.
            owner.accept_heating_step(start_s, dt_s, thermal.time_s())?;
            return check_gate(self.gate, check).map_err(Into::into);
        }
        self.advance_native_momentum(start_s, dt_s, thermal.time_s(), old, check)
    }

    pub(super) fn prepare_unsteady_heat(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, MaterialFlowError<O::Error>> {
        if !self
            .unsteady
            .as_ref()
            .is_some_and(MaterialMomentum::has_heating)
        {
            return Ok(BTreeMap::new());
        }
        self.check_unsteady_start(thermal, check)?;
        if thermal.time_s().to_bits() != start_s.to_bits() {
            return Err(
                invalid("transient heat must prepare from the actual thermal start").into(),
            );
        }
        let applied = self.applied.as_ref().ok_or(invalid(
            "transient heat lacks its source-start material epoch",
        ))?;
        let mut old = BTreeMap::new();
        for sample in &applied.samples {
            check_gate(self.gate, check)?;
            let cell = thermal
                .cells()
                .get(&sample.cell)
                .ok_or(invalid("transient heat source cell is inactive"))?;
            if cell.specific_h_j_kg.to_bits() != sample.specific_enthalpy_j_kg.to_bits() {
                return Err(invalid(
                    "transient heat must retain its source-start enthalpy samples",
                )
                .into());
            }
            old.insert(sample.cell, thermal.temperature_k(sample.cell)?);
        }
        self.advance_native_momentum(start_s, dt_s, end_s, &old, check)?;
        let step = self
            .unsteady_mechanical_step()
            .ok_or(invalid("native momentum did not supply physical heat"))?;
        let applied = self.applied.as_ref().expect("checked source epoch");
        let mut sources = BTreeMap::new();
        for sample in &applied.samples {
            check_gate(self.gate, check)?;
            let power = step.power.cell_heat_w[sample.voxel];
            if power != 0.0 {
                sources.insert(sample.cell, power);
            }
        }
        Ok(sources)
    }

    fn advance_native_momentum(
        &mut self,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        old: &BTreeMap<CellId, f64>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        check_gate(self.gate, check)?;
        let applied = self.applied.as_ref().ok_or(invalid(
            "retained momentum needs the applied material epoch",
        ))?;
        if applied.time_s.to_bits() != start_s.to_bits() {
            return Err(
                invalid("retained momentum source epoch differs from the thermal start").into(),
            );
        }
        let active = self.momentum_domain(&applied.samples, check)?;
        let domain = active.as_ref().unwrap_or(self.domain);
        self.check_phase_drives(domain, check)?;
        for sample in &applied.samples {
            check_gate(self.gate, check)?;
            if old.get(&sample.cell).map(|t| t.to_bits()) != Some(sample.temperature_k.to_bits()) {
                return Err(invalid(
                    "retained momentum must use the actual starting material samples",
                )
                .into());
            }
        }
        if domain.fluid_count() == 0 {
            self.unsteady.as_mut().expect("checked retained momentum").advance_immobile(
                domain,
                start_s,
                dt_s,
                end_s,
                applied.mushy_drag_pa_s_m2.is_some() || applied.darcy_resistance_pa_s_m2.is_some(),
                applied.external_force_n_m3.is_some(),
                applied.buoyancy_force_n_m3.is_some(),
                self.gate,
            )?;
            return check_gate(self.gate, check).map_err(Into::into);
        }
        let fluid = self.reference_fluid()?;
        let reference = self.reference.ok_or(invalid(
            "retained momentum lacks its admitted mobile reference",
        ))?;
        let mut viscosity = vec![reference.dynamic_viscosity_pa_s(); self.domain.cell_count()];
        for sample in &applied.samples {
            check_gate(self.gate, check)?;
            if let Some(transport) = sample.transport.mobile_transport() {
                viscosity[sample.voxel] = transport.dynamic_viscosity_pa_s();
            }
        }
        // Preserve both signed drives. The source fields remain separate in
        // the applied epoch, while the native RHS receives their actual sum.
        let force_n_m3: Option<Cow<'_, [[f64; 3]]>> = match (
            applied.external_force_n_m3.as_deref(),
            applied.buoyancy_force_n_m3.as_deref(),
        ) {
            (Some(external), Some(buoyancy)) => {
                let mut total = Vec::with_capacity(domain.cell_count());
                for (external, buoyancy) in external.iter().zip(buoyancy) {
                    check_gate(self.gate, check)?;
                    let force = std::array::from_fn(|a| external[a] + buoyancy[a]);
                    if force.iter().any(|f| !f.is_finite()) {
                        return Err(invalid("combined material body force must remain finite").into());
                    }
                    total.push(force);
                }
                Some(Cow::Owned(total))
            }
            (Some(force), None) | (None, Some(force)) => Some(Cow::Borrowed(force)),
            (None, None) => None,
        };
        self.unsteady
            .as_mut()
            .expect("checked retained momentum")
            .advance(
                domain,
                &fluid,
                self.config,
                &viscosity,
                force_n_m3.as_deref(),
                applied.buoyancy_force_n_m3.as_deref(),
                if applied.darcy_resistance_pa_s_m2.is_some() { None }
                    else { applied.mushy_drag_pa_s_m2.as_deref() },
                applied.darcy_resistance_pa_s_m2.as_deref(),
                applied.forchheimer_resistance.as_deref(),
                start_s,
                dt_s,
                end_s,
                self.gate,
            )?;
        let flow = &self.unsteady.as_ref().expect("advanced native momentum").flow;
        self.admit_forchheimer_speed(&applied.porous_samples, flow, check)?;
        check_gate(self.gate, check)?;
        Ok(())
    }
}
