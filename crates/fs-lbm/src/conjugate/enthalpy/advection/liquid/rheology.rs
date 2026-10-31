//! Same pinned liquid/caloric bundle as ordinary thermal feedback, with an
//! explicit physical shear-rate axis owned by native staggered momentum.

use super::super::feedback::rheology::{FlowRheologySource, QuasisteadyRheologyControls};
use super::*;

/// Immutable source binding; no independent material table or caloric chart.
#[derive(Clone)]
pub struct LiquidFlowRheology<'a> {
    input: LiquidFlowEnthalpyInput<'a>,
    rate_axis: String,
}

impl<'a> LiquidFlowRheology<'a> {
    /// Admit the original Cp/rho pins, material identity, fixed density and pure
    /// liquid span. Only T and the explicitly typed physical 1/s rate evolve;
    /// Cp is rate-independent because the retained caloric law is h(T).
    pub fn try_new(
        input: LiquidFlowEnthalpyInput<'a>,
        rate_axis: &str,
    ) -> Result<Self, LiquidEnthalpyError> {
        admit_law_with_rheology(input, false, true, Some(rate_axis))?;
        Ok(Self {
            input,
            rate_axis: rate_axis.to_owned(),
        })
    }

    /// Bind this original liquid/caloric bundle to quasisteady constitutive
    /// closure. Supply one declared initial physical rate per voxel, then
    /// select `with_pressure_work_heating`, `with_viscous_heating` or
    /// `with_body_force_heating` on the returned owner before advancing.
    /// Every trial consumes the original max_solves budget; wrapped source,
    /// conduction, advection and heat hooks receive only the converged epoch.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn quasisteady_feedback<'s, O: DepositionObserver>(
        &'s self,
        thermal: &DepositionSimulation<'_>,
        domain: &'s VoxelDomain,
        config: &'s SimpleConfig,
        inner: O,
        gate: &'s CancelGate,
        max_solves: usize,
        initial_shear_rate_per_s: Vec<f64>,
        controls: QuasisteadyRheologyControls,
    ) -> Result<
        MaterialFlowFeedback<
            's,
            ResolvedLiquidState,
            impl Clone
            + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedLiquidState, MaterialError>
            + use<'s, 'a, O>,
            O,
        >,
        MaterialFlowError<O::Error>,
    > {
        // The already admitted source reference initializes the common owner;
        // the explicit physical seeds below immediately replace that query.
        let rate = self.input.liquid.reference_point().axes()[&self.rate_axis];
        let resolve = move |voxel, cell: &DepositedCell, temperature| {
            self.resolve(voxel, cell, temperature, rate)
        };
        MaterialFlowFeedback::try_new(
            thermal, domain, self.input.native(), config, inner, gate, max_solves, resolve,
        )?
        .with_quasisteady_rheology(thermal, self, initial_shear_rate_per_s, controls)
    }

    /// Create the original retained MAC/enthalpy owner at rest and immediately
    /// admit its actual shear-rate field, including moving tangential walls.
    /// Source k(T,rate) can drive `with_fluid_conduction`; molecular heat is
    /// explicitly selected by `with_unsteady_viscous_heating` or body-force
    /// heating. Finite-atmosphere and existing electrical wrappers retain this
    /// same binding automatically. No new fluid or thermal stepper is used.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn feedback<'s, O: DepositionObserver>(
        &'s self,
        thermal: &DepositionSimulation<'_>,
        domain: &'s VoxelDomain,
        config: &'s SimpleConfig,
        inner: O,
        gate: &'s CancelGate,
        max_solves: usize,
    ) -> Result<
        MaterialFlowFeedback<
            's,
            ResolvedLiquidState,
            impl Clone
            + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedLiquidState, MaterialError>
            + use<'s, 'a, O>,
            O,
        >,
        MaterialFlowError<O::Error>,
    > {
        // The reference rate admits only initialization. Immediately below,
        // native faces replace it with the real cold-start strain invariant.
        let rate = self.input.liquid.reference_point().axes()[&self.rate_axis];
        let resolve = move |voxel, cell: &DepositedCell, temperature| {
            self.resolve(voxel, cell, temperature, rate)
        };
        MaterialFlowFeedback::try_new(
            thermal,
            domain,
            self.input.native(),
            config,
            inner,
            gate,
            max_solves,
            resolve,
        )?
        .with_unsteady_momentum()?
        .with_native_rheology(thermal, self)
    }
}

impl FlowRheologySource<ResolvedLiquidState> for LiquidFlowRheology<'_> {
    fn resolve(
        &self,
        _: usize,
        cell: &DepositedCell,
        temperature: f64,
        rate: f64,
    ) -> Result<ResolvedLiquidState, MaterialError> {
        let invalid = |what: &str| MaterialError::State { what: what.into() };
        let phase = self
            .input
            .caloric
            .curve()
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|error| invalid(&error.to_string()))?;
        if phase.phase() != SolidLiquidPhase::Liquid
            || phase.liquid_mass_fraction() != 1.0
            || phase.temperature_k().to_bits() != temperature.to_bits()
            || !rate.is_finite()
            || rate < 0.0
        {
            return Err(invalid(
                "liquid rheology requires the original pure-liquid T(h) and finite physical rate",
            ));
        }
        let law = self.input.liquid;
        let point = law
            .reference_point()
            .clone()
            .with(law.temperature_axis(), temperature)
            .and_then(|point| point.with(&self.rate_axis, rate))
            .map_err(|_| invalid("invalid native liquid temperature/shear coordinates"))?;
        let material = law
            .reference_state()
            .at_query_point(law.card(), &point)
            .map_err(|error| invalid(&error.to_string()))?;
        if material.density_kg_m3().to_bits() != law.reference_state().density_kg_m3().to_bits()
            || !same_positive(phase.bulk_density_kg_m3(), material.density_kg_m3())
        {
            return Err(invalid(
                "native liquid rheology cannot change fixed mass/volume density",
            ));
        }
        Ok(material)
    }
}
