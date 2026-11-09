//! Positive Maxwell heat from the original exponential impact trajectory.
//!
//! The complete mechanical solve supplies the source integrals. Finite material
//! enthalpy is staged before its existing acceptance gate; every history, clock,
//! and thermal store stays unchanged on refusal. Material coefficients remain
//! frozen over each sample and are not inferred from these temperatures.
use super::{
    CaloricReservoir, HeatedImpactFrame, MemoryRoute, RelaxationHeatRoute, caloric::CaloricDeposit,
    invalid, positive,
};
use crate::render::plate::impact::{
    ImpactError, audio::ImpactSource, exponential::ExponentialImpactSystem,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_phs::linear_exact::DenseExponentialBudget;

/// Exact-trajectory material loss into finite source-bound caloric stores.
///
/// Every selected original Maxwell branch has one explicit destination. Several
/// branches may heat the same specimen; unselected resistance, radiation and
/// felt conditioning stay in the original mechanical account. No branch storage
/// is converted to heat and no dissipation is charged twice.
///
/// The exponential discrete-gradient owner's trajectory is exact for the
/// quadratic held-input limit. Nonlinear material/contact heat and motion still
/// require time refinement. This wrapper supplies one-way thermalization, not
/// thermal expansion, changing moduli, or temperature feedback.
pub struct HeatedExponentialImpactSystem<'a> {
    mechanics: ExponentialImpactSystem,
    stores: Vec<CaloricReservoir<'a>>,
    routes: Vec<MemoryRoute>,
    candidate: Vec<CaloricDeposit>,
    heat: Vec<Dd>,
    heat_tolerance_j: f64,
}

impl<'a> HeatedExponentialImpactSystem<'a> {
    /// Attach at the current accepted mechanical time without replaying old loss.
    ///
    /// Routes use the original relaxation insertion order, including material
    /// body histories. Preparation retains complete acoustic/contact/body state
    /// and the current Newton limit. Explicit dense allowances include all source
    /// integral caches; no physical mode or relaxation branch is discarded.
    pub fn new(
        mechanics: ExponentialImpactSystem,
        stores: Vec<CaloricReservoir<'a>>,
        routes: Vec<RelaxationHeatRoute>,
        heat_tolerance_j: f64,
        budget: DenseExponentialBudget,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if stores.is_empty()
            || stores.len() > 64
            || routes.is_empty()
            || routes.len() > 256
            || !positive(heat_tolerance_j)
        {
            return Err(invalid(
                "exponential heating needs 1..64 stores, 1..256 original arms and a positive joule allowance",
            ));
        }
        for (i, store) in stores.iter().enumerate() {
            if stores[..i].iter().any(|old| old.key == store.key) {
                return Err(invalid("duplicate exponential heat store"));
            }
        }
        let mut compiled: Vec<MemoryRoute> = Vec::with_capacity(routes.len());
        for route in routes {
            if route.branch >= mechanics.relaxation_branch_count()
                || compiled.iter().any(|old| old.branch == route.branch)
            {
                return Err(invalid(
                    "exponential heat must name an existing material arm exactly once",
                ));
            }
            let to = stores
                .iter()
                .position(|store| store.key == route.to)
                .ok_or_else(|| invalid("unknown exponential material heat store"))?;
            compiled.push(MemoryRoute {
                branch: route.branch,
                to,
            });
        }
        compiled.sort_by_key(|route| route.branch);
        let dimension = mechanics.state().len();
        // Admit even the transient source witnesses before allocating them.
        // The native owner separately admits all flow/integral caches and work.
        let witness_bytes =
            16_u128 * compiled.len() as u128 * dimension as u128 * dimension as u128;
        if witness_bytes > budget.max_dense_bytes as u128 {
            return Err(invalid(
                "exponential heat source factors exceed the byte allowance",
            ));
        }
        let memory = mechanics
            .relaxation
            .as_ref()
            .ok_or_else(|| invalid("exponential material heat needs the original memory owner"))?;
        let mut channels = Vec::with_capacity(compiled.len());
        for route in &compiled {
            if gate.is_requested() {
                return Err(ImpactError::Cancelled);
            }
            channels.push(memory.branch_resistance_witness(route.branch, dimension)?);
        }
        let mechanics = mechanics.with_heat_channels(&channels, budget, gate)?;
        let candidate = stores
            .iter()
            .map(|store| store.stage_deposit(0.0, heat_tolerance_j))
            .collect::<Result<Vec<_>, _>>()?;
        let heat = vec![Dd::from_f64(0.0); stores.len()];
        Ok(Self {
            mechanics,
            stores,
            routes: compiled,
            candidate,
            heat,
            heat_tolerance_j,
        })
    }

    /// Original read-only accepted state, material sources and motion observations.
    pub fn mechanics(&self) -> &ExponentialImpactSystem {
        &self.mechanics
    }

    /// Original finite charts and accepted thermal inventory.
    pub fn reservoirs(&self) -> &[CaloricReservoir<'a>] {
        &self.stores
    }

    /// Bound Newton work without changing any physical coefficient.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.mechanics.set_iteration_limit(limit)
    }

    /// Separate the current owners without changing state, clock or cumulative heat.
    pub fn into_parts(self) -> (ExponentialImpactSystem, Vec<CaloricReservoir<'a>>) {
        (self.mechanics, self.stores)
    }

    /// Advance one held-force sample and publish mechanics and heat together.
    ///
    /// Warm work reuses all host buffers and native positive integral factors.
    /// Constitutive callbacks retain their own allocation and refinement limits.
    pub fn step(
        &mut self,
        forces: &[f64],
        gate: &CancelGate,
    ) -> Result<HeatedImpactFrame, ImpactError> {
        let Self {
            mechanics,
            stores,
            routes,
            candidate,
            heat,
            heat_tolerance_j,
        } = self;
        heat.fill(Dd::from_f64(0.0));
        let mut deposited = 0.0;
        let mut thermal_change = 0.0;
        let mechanical = mechanics.step_checked(forces, gate, |record, channels| {
            if channels.len() != routes.len() {
                return Err(invalid(
                    "exponential heat lost its original channel mapping",
                ));
            }
            let mut total = Dd::from_f64(0.0);
            for (route, &value) in routes.iter().zip(channels) {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                if !value.is_finite() || value < 0.0 {
                    return Err(invalid(
                        "original exponential material heat is unrepresentable",
                    ));
                }
                let value = Dd::from_f64(value);
                heat[route.to] = heat[route.to] + value;
                total = total + value;
            }
            deposited = total.to_f64();
            let roundoff = 256.0_f64.max(routes.len() as f64)
                * f64::EPSILON
                * (record.dissipated.abs() + deposited.abs());
            if !deposited.is_finite() || deposited > record.dissipated + roundoff {
                return Err(invalid(
                    "exponential material heat exceeds original total dissipation",
                ));
            }
            let mut change = Dd::from_f64(0.0);
            for (i, store) in stores.iter().enumerate() {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                candidate[i] = store.stage_deposit(heat[i].to_f64(), *heat_tolerance_j)?;
                change = change + Dd::from_f64(candidate[i].energy_change_j);
            }
            thermal_change = change.to_f64();
            if !thermal_change.is_finite() || (thermal_change - deposited).abs() > *heat_tolerance_j
            {
                return Err(invalid(
                    "exponential caloric transfer exceeds its joule allowance",
                ));
            }
            Ok(())
        })?;
        // Mechanics performs the final cancellation and physical-state checks.
        // This publication is infallible and invokes no callbacks.
        for (store, &deposit) in stores.iter_mut().zip(candidate.iter()) {
            store.accept_deposit(deposit);
        }
        Ok(HeatedImpactFrame {
            deposited_heat_j: deposited,
            thermal_energy_change_j: thermal_change,
            remaining_mechanical_loss_j: mechanical.dissipated_energy_j - deposited,
            combined_energy_residual_j: mechanical.balance_residual_j + thermal_change - deposited,
            mechanical,
        })
    }
}

impl ImpactSource for HeatedExponentialImpactSystem<'_> {
    type Frame = HeatedImpactFrame;
    fn state(&self) -> &[f64] {
        self.mechanics.state()
    }
    fn mode_count(&self) -> usize {
        self.mechanics.mode_count()
    }
    fn samples_rendered(&self) -> u64 {
        self.mechanics.samples_rendered()
    }
    fn sample_period_s(&self) -> f64 {
        self.mechanics.sample_period_s()
    }
    fn remaining_steps(&self) -> u64 {
        self.mechanics.remaining_steps()
    }
    fn maximum_generalized_force(&self) -> f64 {
        self.mechanics.maximum_generalized_force()
    }
    fn advance(&mut self, forces: &[f64], gate: &CancelGate) -> Result<Self::Frame, ImpactError> {
        self.step(forces, gate)
    }
}
