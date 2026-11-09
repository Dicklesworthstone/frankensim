//! Original material enthalpy in the native immersed-body heat balance.
//!
//! With the original post-collision traces frozen for one native tick, solve
//! `m (h_new-h_old)/E_ref + K (theta(h_new)-theta_trace) = Q`.
//! The monotone material chart includes its exact sensible Cp(T) primitive and
//! latent plateaus. A safeguarded Newton/bisection solve acts on enthalpy, never
//! on an invented large apparent heat capacity. Actual reflected populations
//! and actual m*delta_h still pass the original native first-law gates.
//! Caloric convention: <https://wtt-pro.nist.gov/wtt-pro/help/properties/enthalpy_cond.html>.

use super::{ThermalOpenError3, ThermalOpenField3, ThermalOpenLimits3, invalid};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use std::sync::Arc;

mod balance;
use balance::solve_balance;

/// Physical interpretation of the SAME native heat and temperature coordinate.
/// Fixed body mass is independent of the chart's temperature-dependent density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImmersedCaloricBinding3 {
    pub mass_kg: f64,
    pub initial_specific_enthalpy_j_kg: f64,
    /// Explicit supported operating span; no material extrapolation is allowed.
    /// Its endpoint temperatures must fit the native thermal limits.
    pub specific_enthalpy_range_j_kg: [f64; 2],
    /// `T_K = temperature_origin_k + temperature_scale_k * theta_native`.
    pub temperature_origin_k: f64,
    /// Positive so positive native heat retains its physical sign.
    pub temperature_scale_k: f64,
    /// The reference fluid-cell capacity already used to interpret D3Q7 heat.
    /// One native energy unit equals this capacity times temperature_scale_k.
    pub fluid_cell_heat_capacity_j_k: f64,
}

/// Bounded material queries and a physical root-solve allowance. The original
/// local and whole-grid heat-balance limits remain independently authoritative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImmersedCaloricControls3 {
    pub maximum_iterations: usize,
    /// Bounds the retained immutable chart's original knot count. Its external
    /// Arc-owned source storage is not copied into the native body allocation.
    pub maximum_source_knots: usize,
    pub absolute_energy_tolerance_j: f64,
}

/// The body's ONE accepted material inventory. Sensible and latent heat are
/// already included in this enthalpy and are never credited separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImmersedCaloricState3 {
    pub mass_kg: f64,
    pub phase: EquilibriumPhaseState,
}

/// Original source error, domain limit, unrepresentable balance or exhausted
/// bounded solve. Residuals use the original native energy coordinate.
pub type ImmersedCaloricError3 = balance::CaloricBalanceError<PhaseStateError>;

pub(super) struct CaloricBody3 {
    curve: Arc<EquilibriumEnthalpyPhaseCurve>,
    binding: ImmersedCaloricBinding3,
    controls: ImmersedCaloricControls3,
    energy_scale_j: f64,
    mass_per_native_energy: f64,
    trial: ImmersedCaloricState3,
}

impl ThermalOpenField3 {
    /// Replace fixed-C storage of one newly installed independent immersed
    /// body with its original immutable material chart, including latent heat.
    ///
    /// The supplied h must reproduce the body's installed native temperature
    /// EXACTLY. Binding after an accepted body exchange, rebinding, and binding
    /// any conducting body/voxel graph are refused. Thus this operation neither
    /// resets a physical trajectory nor silently reuses a fixed-C graph factor.
    /// Original resistance, held heat source, links and fluid populations stay
    /// in place. The initial fixed heat capacity becomes inactive declaration
    /// data; all future storage and acceptance use the single m*h inventory.
    /// Every caloric body bound to this field must share identical temperature
    /// origin/scale and reference fluid-cell capacity, so shared native traces
    /// and summed heat retain one physical interpretation.
    ///
    /// The body remains isothermal and stationary, even when its caloric chart
    /// contains liquid. This models a lumped, retained phase-change store; it
    /// does not melt/remove the native obstacle, move liquid, infer geometry
    /// from density, resolve internal temperature gradients, or perform
    /// pressure/expansion work. The caller must justify those approximations
    /// and supply the SAME physical normalization as the native fluid field.
    /// Preparation is allocation-free apart from the caller-owned shared chart;
    /// accepted body state is changed only after all admission succeeds.
    pub fn set_immersed_body_caloric(
        &mut self,
        body: usize,
        curve: Arc<EquilibriumEnthalpyPhaseCurve>,
        binding: ImmersedCaloricBinding3,
        controls: ImmersedCaloricControls3,
    ) -> Result<(), ThermalOpenError3> {
        let owner = self
            .immersed
            .as_mut()
            .ok_or_else(|| invalid("caloric binding requires an installed immersed body"))?;
        if owner.conduction.is_some() {
            return Err(invalid(
                "caloric binding does not support a fixed-C conducting graph",
            ));
        }
        if owner
            .bodies
            .iter()
            .filter_map(|node| node.caloric.as_ref())
            .any(|caloric| {
                caloric.binding.temperature_origin_k != binding.temperature_origin_k
                    || caloric.binding.temperature_scale_k != binding.temperature_scale_k
                    || caloric.binding.fluid_cell_heat_capacity_j_k
                        != binding.fluid_cell_heat_capacity_j_k
            })
        {
            return Err(invalid(
                "caloric bodies sharing native fluid must share temperature and energy normalization",
            ));
        }
        let original = owner
            .bodies
            .get_mut(body)
            .ok_or_else(|| invalid("unknown immersed thermal body"))?;
        if original.caloric.is_some() || original.accepted.last_exchange.is_some() {
            return Err(invalid(
                "caloric storage must be bound once before the body's first step",
            ));
        }
        let [lo, hi] = binding.specific_enthalpy_range_j_kg;
        let energy_scale_j = binding.fluid_cell_heat_capacity_j_k * binding.temperature_scale_k;
        let coefficient = binding.mass_kg / energy_scale_j;
        if !binding.mass_kg.is_finite()
            || binding.mass_kg <= 0.0
            || !binding.temperature_origin_k.is_finite()
            || !binding.temperature_scale_k.is_finite()
            || binding.temperature_scale_k <= 0.0
            || !binding.fluid_cell_heat_capacity_j_k.is_finite()
            || binding.fluid_cell_heat_capacity_j_k <= 0.0
            || !energy_scale_j.is_finite()
            || energy_scale_j <= 0.0
            || !coefficient.is_finite()
            || coefficient <= 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || hi <= lo
            || !binding.initial_specific_enthalpy_j_kg.is_finite()
            || binding.initial_specific_enthalpy_j_kg < lo
            || binding.initial_specific_enthalpy_j_kg > hi
            || !(1..=256).contains(&controls.maximum_iterations)
            || !(2..=65_536).contains(&controls.maximum_source_knots)
            || curve.knots().len() > controls.maximum_source_knots
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(controls.absolute_energy_tolerance_j / energy_scale_j).is_finite()
        {
            return Err(invalid(
                "invalid immersed caloric mass, scaling, domain or work controls",
            ));
        }
        let source_error = |source| ThermalOpenError3::Caloric {
            body,
            source: ImmersedCaloricError3::Material(source),
        };
        let phase = curve
            .state_at_specific_enthalpy(binding.initial_specific_enthalpy_j_kg)
            .map_err(source_error)?;
        for h in [lo, hi] {
            let t = curve
                .temperature_at_specific_enthalpy(h)
                .map_err(source_error)?;
            let native = (t - binding.temperature_origin_k) / binding.temperature_scale_k;
            super::check_boundary_temperature(6, native, self.limits)?;
        }
        let native =
            (phase.temperature_k() - binding.temperature_origin_k) / binding.temperature_scale_k;
        if native != original.accepted.solid.temperature {
            return Err(invalid(
                "caloric initial state must reproduce the original body temperature",
            ));
        }
        let state = ImmersedCaloricState3 {
            mass_kg: binding.mass_kg,
            phase,
        };
        original.caloric = Some(CaloricBody3 {
            curve,
            binding,
            controls,
            energy_scale_j,
            mass_per_native_energy: coefficient,
            trial: state,
        });
        original.accepted.caloric = Some(state);
        Ok(())
    }

    /// Original pinned material curve, including its source receipts. Reading
    /// it does not recreate the body, reselect material or advance any heat.
    pub fn immersed_body_caloric_curve(
        &self,
        body: usize,
    ) -> Option<&EquilibriumEnthalpyPhaseCurve> {
        Some(
            &self
                .immersed
                .as_ref()?
                .bodies
                .get(body)?
                .caloric
                .as_ref()?
                .curve,
        )
    }
}

impl CaloricBody3 {
    pub(super) fn solve<E: From<ThermalOpenError3>>(
        &mut self,
        body: usize,
        old: ImmersedCaloricState3,
        source: f64,
        conductance: f64,
        trace: f64,
        limits: ThermalOpenLimits3,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<f64, E> {
        let map = |source| E::from(ThermalOpenError3::Caloric { body, source });
        let old_h = old.phase.specific_enthalpy_j_kg();
        let tolerance = limits
            .balance_tolerance
            .min(self.controls.absolute_energy_tolerance_j / self.energy_scale_j);
        let mut evaluate = |h| {
            let temperature = self
                .curve
                .temperature_at_specific_enthalpy(h)
                .map_err(ImmersedCaloricError3::Material)?;
            let slope = self
                .curve
                .temperature_derivative_at_specific_enthalpy(h)
                .map_err(ImmersedCaloricError3::Material)?;
            Ok((
                (temperature - self.binding.temperature_origin_k)
                    / self.binding.temperature_scale_k,
                slope / self.binding.temperature_scale_k,
            ))
        };
        let solution = solve_balance(
            old_h,
            self.binding.specific_enthalpy_range_j_kg,
            self.mass_per_native_energy,
            source,
            conductance,
            trace,
            tolerance,
            self.controls.maximum_iterations,
            &mut evaluate,
            &mut poll,
            &map,
        )?;
        let phase = self
            .curve
            .state_at_specific_enthalpy(solution.0)
            .map_err(|error| map(ImmersedCaloricError3::Material(error)))?;
        super::check_boundary_temperature(6, solution.1, limits)?;
        self.trial = ImmersedCaloricState3 {
            mass_kg: old.mass_kg,
            phase,
        };
        Ok(solution.1)
    }

    pub(super) fn stored_change(&self, old: ImmersedCaloricState3) -> f64 {
        self.mass_per_native_energy
            * (self.trial.phase.specific_enthalpy_j_kg() - old.phase.specific_enthalpy_j_kg())
    }

    pub(super) fn trial_state(&self) -> ImmersedCaloricState3 {
        self.trial
    }

    pub(super) fn balance_tolerance(&self) -> f64 {
        self.controls.absolute_energy_tolerance_j / self.energy_scale_j
    }
}

#[cfg(test)]
mod tests;
