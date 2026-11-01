//! Bounded constitutive fixed point around the ORIGINAL native flow stage.
//! Iterations query one immutable source at the current caloric temperature;
//! they never call a wrapped thermal observer or deposit candidate heat.
use super::*;

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Close source viscosity against the actual quasisteady native strain.
    ///
    /// Supply one explicit finite nonnegative initial rate per voxel; zero is
    /// not presumed to lie in a material's support. The immutable source admits
    /// those seeds immediately. Before advancing, select the original closed
    /// viscous/body-force or reservoir pressure-work heating owner. Retained
    /// inertia and buoyancy remain separate modes.
    ///
    /// Each trial uses the original source/porous/native Stokes stage. Actual
    /// symmetric-strain RMS rates, including half-cell walls and reservoir
    /// mouth edges, drive the next source query. Rate AND transport closure are
    /// required; the exact queried coordinates remain on applied(). Thermal
    /// conductivity and viscosity-dependent porous resistance therefore share
    /// that same accepted source epoch. Every native trial consumes the existing
    /// max_solves budget. Only the converged flow, heat and advection graph reach
    /// the original thermal observer, and every failure rolls back its candidate.
    ///
    /// This fixed-point iteration has no constitutive memory or uniqueness
    /// certificate. Source support, relaxation/iteration and mesh refinement
    /// remain explicit. Thermal endpoint T is re-admitted at the last solved
    /// physical rate; a fresh equilibrium is solved at the NEXT source epoch.
    pub fn with_quasisteady_rheology(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'a dyn FlowRheologySource<F>,
        initial_shear_rate_per_s: Vec<f64>,
        controls: QuasisteadyRheologyControls,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        controls.admit()?;
        self.check_map(thermal, &mut || Ok(()))?;
        if self.unsteady.is_some()
            || self.buoyancy.is_some()
            || self.rheology.is_some()
            || self.applied.is_some()
        {
            return Err(invalid("attach one quasisteady source before the first advance, without retained inertia or buoyancy").into());
        }
        if initial_shear_rate_per_s.len() != self.domain.cell_count()
            || initial_shear_rate_per_s
                .iter()
                .enumerate()
                .any(|(c, rate)| {
                    !rate.is_finite() || *rate < 0.0 || (!self.domain.is_fluid(c) && *rate != 0.0)
                })
        {
            return Err(invalid("quasisteady rheology needs one finite nonnegative rate seed per voxel, zero in fixed solid geometry").into());
        }
        self.rheology = Some(RheologyBinding {
            source,
            rates: initial_shear_rate_per_s,
            quasisteady: Some(controls),
            report: None,
        });
        self.endpoint = self.sample_materials(thermal, &mut || Ok(()))?;
        admit_density(&self.endpoint, self.density_kg_m3)?;
        self.mushy_drag(&self.endpoint, &mut || Ok(()))?;
        self.check_porous_endpoint(&mut || Ok(()))?;
        self.check_phase_conduction(&self.endpoint, &mut || Ok(()))?;
        Ok(self)
    }

    pub fn quasisteady_rheology_report(&self) -> Option<QuasisteadyRheologyReport> {
        self.rheology
            .as_ref()
            .filter(|binding| binding.quasisteady.is_some())
            .and_then(|binding| binding.report)
    }

    pub(in crate::conjugate::enthalpy::advection::feedback) fn solve_rheology_epoch(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(MaterialFlowEpoch<F, T>, Option<FluidTransportState>), MaterialFlowError<O::Error>>
    {
        let Some(controls) = self
            .rheology
            .as_ref()
            .and_then(|binding| binding.quasisteady)
        else {
            return self.solve_material_flow_epoch(thermal, check);
        };
        let heat = self.viscous_heating.as_ref().ok_or(invalid(
            "quasisteady rheology needs the original closed or pressure-work Stokes heating owner",
        ))?;
        let (heat_controls, pressure_work) = (heat.controls, heat.pressure_work);
        for iteration in 1..=controls.max_iterations {
            check_gate(self.gate, check)?;
            let (epoch, reference) = self.solve_material_flow_epoch(thermal, check)?;
            let active = self.momentum_domain(&epoch.samples, check)?;
            let domain = active.as_ref().unwrap_or(self.domain);
            let rates = if domain.fluid_count() == 0 {
                vec![0.0; domain.cell_count()]
            } else {
                let source = reference.ok_or(invalid(
                    "quasisteady rate observation lacks mobile source transport",
                ))?;
                let fluid = FluidProperties {
                    density_kg_m3: self.density_kg_m3,
                    specific_heat_j_kg_k: source.specific_heat_j_kg_k(),
                    conductivity_w_m_k: source.thermal_conductivity_w_m_k(),
                    kinematic_viscosity_m2_s: source.kinematic_viscosity_m2_s(),
                };
                crate::conjugate::simple::viscous_heat::physical_flow_shear_rates(
                    domain,
                    &fluid,
                    &epoch.flow.field,
                    self.config,
                    heat_controls,
                    pressure_work,
                    self.gate,
                )?
            };
            let (_, closed, _) =
                self.observe_rheology(thermal, &epoch, rates, iteration, controls, check)?;
            let queried = epoch.shear_rate_per_s.as_ref().expect("bound source epoch");
            let binding = self.rheology.as_mut().expect("binding retained");
            if closed {
                return Ok((epoch, reference));
            }
            // Relax source coordinates, never velocities, forces or generated
            // energy. The next native trial rebuilds all coupled coefficients.
            if controls.relaxation != 1.0 {
                for (actual, &query) in binding.rates.iter_mut().zip(queried) {
                    *actual = controls
                        .relaxation
                        .mul_add(*actual, (1.0 - controls.relaxation) * query);
                }
            }
        }
        Err(DepositionError::Budget("quasisteady constitutive iterations").into())
    }
}
