//! Source-bound liquid enthalpy on the native stationary-skeleton Stokes flow.

use super::*;
use fs_conduction::transient::deposition::advection::liquid::LiquidCarrierKnot;

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Carry liquid enthalpy through the stationary solid while the original
    /// thermal cells retain the complete sensible/latent mixture energy.
    ///
    /// Requires the already bound original pure-substance phase chart, equal
    /// constant phase density, and creeping flow with physical viscous heating
    /// (`with_viscous_heating`, `with_body_force_heating`, or the existing
    /// `with_pressure_work_heating` reservoir-mouth model). The native
    /// velocity is interpreted as liquid volume flux per total cell-face area.
    /// The stationary solid carries no mass flux or reaction work. Every native
    /// liquid mass rate remains locally balanced, so changes in liquid inventory
    /// are exactly the melting/freezing supplied by the original h chart.
    ///
    /// At coexistence a donor carries its liquid-side enthalpy, independently
    /// of how much stationary solid remains; fully liquid donors retain their
    /// original h. Isothermal circulation therefore cannot homogenize solid
    /// fractions. Superheat, conduction, Joule and viscous heat still melt the
    /// same original cells, and latent energy is never supplied twice.
    ///
    /// Open flow uses the ordinary mobile mushy entry, explicit pure-liquid
    /// inlet enthalpy on every actual inflow (including backflow), and the
    /// original pressure-work heating owner. Outlet enthalpy is the current
    /// liquid's h. Stored and incoming h must use the same fixed thermodynamic
    /// reference pressure. Do not add hydrodynamic p/rho to inlet h: the native
    /// mechanical pressure supply already enters once through physical heat.
    /// Its pressure gauge remains independent of the caloric enthalpy datum.
    ///
    /// Attach before the first advance. Exact-solid topology remains closed;
    /// open flow refuses below the source mobility's admitted liquid fraction.
    /// Transient homogeneous inertia, non-isothermal alloy coexistence, unequal
    /// density and moving skeletons are unsupported. No liquid-occupancy inertia or separate
    /// phase momentum is inferred. Explicit and implicit Stokes thermal paths
    /// use the same original `Advected` owner and source-bound descriptor.
    pub fn with_liquid_carrier(mut self) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if self.unsteady.is_some() || self.applied.is_some() || self.viscous_heating.is_none() {
            return Err(invalid(
                "liquid carriage requires physical Stokes heating before the first advance",
            ).into());
        }
        if self.viscous_heating.as_ref().is_some_and(|heat| heat.pressure_work) {
            if self.exact_phase {
                return Err(invalid("open liquid carriage has no exact-solid pressure-mouth closure").into());
            }
            // The existing pressure-work attachment already admits the native
            // opposed normal-traction mouths and their physical power account.
        } else {
            phase::admit_phase_config(self.config)?;
        }
        let phase = self.mushy.ok_or(invalid(
            "liquid carriage requires the original stationary phase-mobility source",
        ))?;
        if phase.fixed_density_kg_m3() != Some(self.density_kg_m3) {
            return Err(invalid("liquid carriage requires equal constant phase density").into());
        }
        let mut knots = Vec::with_capacity(phase.curve().knots().len());
        for (index, knot) in phase.curve().knots().iter().enumerate() {
            if index % 256 == 0 { check_gate(self.gate, &mut || Ok(()))?; }
            knots.push(LiquidCarrierKnot {
                specific_enthalpy_j_kg: knot.specific_enthalpy_j_kg,
                temperature_k: knot.temperature_k,
                liquid_mass_fraction: knot.liquid_mass_fraction,
                bulk_density_kg_m3: knot.bulk_density_kg_m3,
            });
        }
        let carrier = LiquidEnthalpyCarrier::from_source_chart(
            self.input.material_id, self.input.caloric, &knots,
            || check_gate(self.gate, &mut || Ok(())),
        )?;
        self.advected = Some(self.advected.take()
            .ok_or(invalid("unfinished material-flow source epoch"))?
            .with_liquid_carrier(carrier)?);
        self.liquid_carrier = Some(carrier);
        Ok(self)
    }
}

impl<F, R, O, T> MaterialFlowFeedback<'_, F, R, O, T> {
    /// The original liquid-carrier chart, absent for homogeneous bulk-h transport.
    pub fn liquid_carrier(&self) -> Option<&LiquidEnthalpyCarrier> {
        self.liquid_carrier.as_ref()
    }

    pub(super) fn bind_liquid_carrier<P>(
        &self,
        advected: Advected<P>,
    ) -> Result<Advected<P>, DepositionError> {
        match self.liquid_carrier {
            Some(carrier) => advected.with_liquid_carrier(carrier),
            None => Ok(advected),
        }
    }
}
