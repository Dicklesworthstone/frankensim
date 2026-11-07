//! Adaptive physical time for the native heat-before-motion cavity split.
//!
//! A trial compares one complete interval with two sequential half intervals.
//! Only the fine path is retained; no extrapolated plastic state, gas state or
//! caloric history is invented. All three solves consume the same global heat
//! trial budget. Native source, nonlinear and observer refusals remain explicit.

use super::*;

/// Local error controls for a first-order physical compression stroke.
///
/// Relative scales use changes over this interval, never absolute enthalpy,
/// formation energy, displacement or plastic-metric reference values. The step
/// bounds refer to a comparison interval; its accepted physical steps are its
/// two halves. A short final remainder may be below `minimum_step_s`.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveCavityHeatSolve {
    pub initial_step_s: f64,
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Includes rejected comparison intervals.
    pub maximum_intervals: usize,
    /// Includes every coarse, fine, Jacobian and final native heat evaluation.
    pub maximum_heat_trials: usize,
    pub absolute_displacement_tolerance_m: f64,
    pub absolute_velocity_tolerance_mps: f64,
    pub absolute_temperature_tolerance_k: f64,
    pub absolute_specific_enthalpy_tolerance_j_kg: f64,
    pub absolute_gas_volume_tolerance_m3: f64,
    /// Porous inverse plastic metric, accumulated strain and plastic Jacobian.
    pub absolute_plastic_state_tolerance: f64,
    /// Each gas's direct storage and heat, and each native work/energy channel.
    pub absolute_energy_tolerance_j: f64,
    pub relative_tolerance: f64,
    /// Whole-stroke absolute residual sum, interface defect and storage closure.
    pub window_energy_tolerance_j: f64,
}

/// Raw full-step versus two-half-step differences, without Richardson scaling.
#[derive(Clone, Copy, Debug, Default)]
pub struct CavityHeatLocalError {
    pub displacement_m: f64,
    pub velocity_mps: f64,
    pub temperature_k: f64,
    pub specific_enthalpy_j_kg: f64,
    pub gas_volume_m3: f64,
    pub plastic_state: f64,
    pub energy_j: f64,
    pub normalized: f64,
}

/// A refusal publishes no part of any provisional fine or coarse trajectory.
#[derive(Debug)]
pub enum CavityHeatMarchError<E = Infallible> {
    Invalid(&'static str),
    Native {
        start_s: f64,
        target_s: f64,
        source: CavityHeatError<E>,
    },
    Thermal(DepositionError),
    Gas(GasEnergyError),
    Budget {
        time_s: f64,
        intervals: usize,
        heat_trials: usize,
    },
    Accuracy {
        time_s: f64,
        step_s: f64,
        error: CavityHeatLocalError,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for CavityHeatMarchError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive thermal porous cavity stroke refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CavityHeatMarchError<E> {}

/// Physical first-law terms summed only along accepted fine steps.
///
/// Native endpoint work defects and backward-Euler kinetic loss remain named
/// numerical accounts; neither is deposited as material or gas heat.
#[derive(Clone, Debug, Default)]
pub struct CavityHeatWindowEnergy {
    pub gas_energy_change_j: f64,
    pub thermal_storage_change_j: f64,
    pub solid_internal_energy_change_j: f64,
    pub kinetic_change_j: f64,
    pub kinetic_numerical_loss_j: f64,
    pub external_work_j: f64,
    pub mechanical_work_defect_j: f64,
    pub external_gas_heat_j: f64,
    pub external_solid_heat_j: f64,
    pub incoming_thermal_energy_j: f64,
    pub gas_to_solid_heat_j: f64,
    pub material_heat_j: f64,
    pub plastic_dissipation_j: f64,
    pub interface_heat_defect_j: f64,
    pub combined_residual_j: f64,
    pub absolute_residual_sum_j: f64,
    /// Sum of absolute direct-original-gas versus accumulated storage differences.
    pub gas_storage_defect_j: f64,
}

/// Compact stroke account; full native receipts remain on the last endpoint.
#[derive(Clone, Debug)]
pub struct CavityHeatMarchReport {
    pub start_s: f64,
    pub end_s: f64,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub advance_calls: usize,
    pub heat_trials: usize,
    pub accepted_thermal_substeps: usize,
    pub next_step_s: f64,
    /// Componentwise maxima over accepted comparison intervals.
    pub maximum_accepted_error: CavityHeatLocalError,
    pub energy: CavityHeatWindowEnergy,
}

/// The endpoint owns the actual final gas/solid/thermal/observer histories.
/// Its detailed balance describes the last physical half step, while `report`
/// accounts for the complete accepted stroke. Memory does not grow with the
/// number of time steps or store a duplicate mesh for each historical endpoint.
pub struct AdaptiveCavityHeatCandidate<'a, 't, O, G = ClosedGasParcel> {
    pub endpoint: CavityHeatCandidate<'a, 't, O, G>,
    pub report: CavityHeatMarchReport,
}

fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), CavityHeatMarchError<E>> {
    cx.checkpoint().map_err(|_| CavityHeatMarchError::Cancelled)
}
fn finite<E>(x: f64) -> Result<f64, CavityHeatMarchError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(CavityHeatMarchError::Invalid(
            "unrepresentable adaptive time or energy account",
        ))
    }
}

impl AdaptiveCavityHeatSolve {
    fn validate<E>(self, heat: CavityHeatSolve) -> Result<(), CavityHeatMarchError<E>> {
        let positive = [
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
            self.absolute_displacement_tolerance_m,
            self.absolute_velocity_tolerance_mps,
            self.absolute_temperature_tolerance_k,
            self.absolute_specific_enthalpy_tolerance_j_kg,
            self.absolute_gas_volume_tolerance_m3,
            self.absolute_plastic_state_tolerance,
            self.absolute_energy_tolerance_j,
            self.window_energy_tolerance_j,
        ];
        if positive.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 0.0
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || self.maximum_intervals == 0
            || self.maximum_heat_trials == 0
            || !heat.energy_tolerance_j.is_finite()
            || heat.energy_tolerance_j <= 0.0
            || heat.energy_tolerance_j > self.absolute_energy_tolerance_j
            || heat.energy_tolerance_j > self.window_energy_tolerance_j
        {
            return Err(CavityHeatMarchError::Invalid(
                "positive bounded time/error controls and compatible native energy gates required",
            ));
        }
        Ok(())
    }
}

// One common list serves both local work/heat comparison and compensated
// accepted-path accounting. Every entry is an interval increment, not an
// absolute formation energy. Internal film heat is deliberately separate.
fn energy_channels<G>(b: &CavityHeatBalance<G>) -> [f64; 16] {
    let n = &b.native;
    let m = &n.material_energy;
    [
        b.total_gas_energy_change_j,
        n.thermal_audit.stored_energy_change_j,
        m.material.internal_energy_change_j,
        n.inertia.as_ref().map_or(0.0, |x| x.kinetic_change_j),
        n.inertia.as_ref().map_or(0.0, |x| x.numerical_loss_j),
        (Dd::from_f64(m.material.dead_load_endpoint_work_j)
            + Dd::from_f64(m.material.support_endpoint_work_j)
            + Dd::from_f64(m.pressure.as_ref().map_or(0.0, |x| x.work_j))
            + Dd::from_f64(m.exterior_work_j))
        .to_f64(),
        m.external_work_defect_j,
        b.gases
            .iter()
            .fold(Dd::ZERO, |sum, gas| sum + Dd::from_f64(gas.external_heat_j))
            .to_f64(),
        b.external_solid_heat_j,
        n.thermal_audit.incoming_energy_j,
        b.gas_to_solid_heat_j,
        n.applied_heat_j,
        m.material.plastic_dissipation_j,
        b.interface_heat_defect_j,
        b.combined_residual_j,
        b.combined_residual_j.abs(),
    ]
}

struct Totals {
    energy: [Dd; 16],
    gas_storage: Vec<Dd>,
    gas_external: Vec<Dd>,
    thermal_substeps: usize,
}
impl Totals {
    fn new(gases: usize) -> Self {
        Self {
            energy: [Dd::ZERO; 16],
            gas_storage: vec![Dd::ZERO; gases],
            gas_external: vec![Dd::ZERO; gases],
            thermal_substeps: 0,
        }
    }
    fn add<E, G>(&mut self, b: &CavityHeatBalance<G>) -> Result<(), CavityHeatMarchError<E>> {
        for (sum, value) in self.energy.iter_mut().zip(energy_channels(b)) {
            *sum = *sum + Dd::from_f64(finite(value)?);
            finite::<E>(sum.to_f64())?;
        }
        for (i, gas) in b.gases.iter().enumerate() {
            self.gas_storage[i] = self.gas_storage[i] + Dd::from_f64(gas.total_energy_change_j);
            self.gas_external[i] = self.gas_external[i] + Dd::from_f64(gas.external_heat_j);
        }
        self.thermal_substeps = self
            .thermal_substeps
            .checked_add(b.native.thermal_audit.substeps)
            .ok_or(CavityHeatMarchError::Invalid(
                "thermal substep count overflow",
            ))?;
        Ok(())
    }
    fn finish<E>(
        &self,
        storage_defect: f64,
    ) -> Result<CavityHeatWindowEnergy, CavityHeatMarchError<E>> {
        let mut v = [0.0; 16];
        for (out, sum) in v.iter_mut().zip(self.energy) {
            *out = finite(sum.to_f64())?;
        }
        Ok(CavityHeatWindowEnergy {
            gas_energy_change_j: v[0],
            thermal_storage_change_j: v[1],
            solid_internal_energy_change_j: v[2],
            kinetic_change_j: v[3],
            kinetic_numerical_loss_j: v[4],
            external_work_j: v[5],
            mechanical_work_defect_j: v[6],
            external_gas_heat_j: v[7],
            external_solid_heat_j: v[8],
            incoming_thermal_energy_j: v[9],
            gas_to_solid_heat_j: v[10],
            material_heat_j: v[11],
            plastic_dissipation_j: v[12],
            interface_heat_defect_j: v[13],
            combined_residual_j: v[14],
            absolute_residual_sum_j: v[15],
            gas_storage_defect_j: storage_defect,
        })
    }
}

fn compare<E>(
    before: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    maximum: &mut f64,
    normalized: &mut f64,
) -> Result<(), CavityHeatMarchError<E>> {
    let error = finite((Dd::from_f64(fine) - Dd::from_f64(coarse)).to_f64())?.abs();
    let scale = finite((coarse - before).abs().max((fine - before).abs()))?;
    let tolerance = finite(absolute + relative * scale)?;
    *maximum = maximum.max(error);
    *normalized = normalized.max(finite(error / tolerance)?);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn estimate<E, O, G: ThermalCavityGas>(
    before: &ThermalPorousGasSystem<'_, G>,
    thermal: &DepositionSimulation<'_>,
    coarse: &CavityHeatCandidate<'_, '_, O, G>,
    half: &CavityHeatCandidate<'_, '_, O, G>,
    fine: &CavityHeatCandidate<'_, '_, O, G>,
    controls: AdaptiveCavityHeatSolve,
    cx: &Cx<'_>,
) -> Result<CavityHeatLocalError, CavityHeatMarchError<E>> {
    let mut error = CavityHeatLocalError::default();
    let rel = controls.relative_tolerance;
    let old = before.solid().mechanics();
    let a = coarse.system.solid().mechanics();
    let b = fine.system.solid().mechanics();
    for ((old, a), b) in old
        .displacement_m()
        .as_flattened()
        .iter()
        .zip(a.displacement_m().as_flattened())
        .zip(b.displacement_m().as_flattened())
    {
        checkpoint(cx)?;
        compare(
            *old,
            *a,
            *b,
            controls.absolute_displacement_tolerance_m,
            rel,
            &mut error.displacement_m,
            &mut error.normalized,
        )?;
    }
    match (old.velocity_mps(), a.velocity_mps(), b.velocity_mps()) {
        (Some(old), Some(a), Some(b)) => {
            for ((old, a), b) in old
                .as_flattened()
                .iter()
                .zip(a.as_flattened())
                .zip(b.as_flattened())
            {
                checkpoint(cx)?;
                compare(
                    *old,
                    *a,
                    *b,
                    controls.absolute_velocity_tolerance_mps,
                    rel,
                    &mut error.velocity_mps,
                    &mut error.normalized,
                )?;
            }
        }
        (None, None, None) => {}
        _ => {
            return Err(CavityHeatMarchError::Invalid(
                "inertial ownership changed during adaptive stroke",
            ));
        }
    }
    for ((old, a), b) in old
        .element_states()
        .iter()
        .zip(a.element_states())
        .zip(b.element_states())
    {
        checkpoint(cx)?;
        for ((old, a), b) in old
            .inverse_plastic_cauchy_green
            .iter()
            .zip(&a.inverse_plastic_cauchy_green)
            .zip(&b.inverse_plastic_cauchy_green)
        {
            compare(
                *old,
                *a,
                *b,
                controls.absolute_plastic_state_tolerance,
                rel,
                &mut error.plastic_state,
                &mut error.normalized,
            )?;
        }
        compare(
            old.alpha,
            a.alpha,
            b.alpha,
            controls.absolute_plastic_state_tolerance,
            rel,
            &mut error.plastic_state,
            &mut error.normalized,
        )?;
    }
    for ((old, a), b) in old
        .element_plastic_jacobians()
        .iter()
        .zip(a.element_plastic_jacobians())
        .zip(b.element_plastic_jacobians())
    {
        checkpoint(cx)?;
        compare(
            *old,
            *a,
            *b,
            controls.absolute_plastic_state_tolerance,
            rel,
            &mut error.plastic_state,
            &mut error.normalized,
        )?;
    }
    if thermal.cells().len() != coarse.thermal.cells().len()
        || thermal.cells().len() != fine.thermal.cells().len()
    {
        return Err(CavityHeatMarchError::Invalid(
            "adaptive stroke requires unchanged thermal cell inventory",
        ));
    }
    for (&id, old) in thermal.cells() {
        checkpoint(cx)?;
        let a = coarse
            .thermal
            .cells()
            .get(&id)
            .ok_or(CavityHeatMarchError::Invalid(
                "adaptive thermal cell disappeared",
            ))?;
        let b = fine
            .thermal
            .cells()
            .get(&id)
            .ok_or(CavityHeatMarchError::Invalid(
                "adaptive thermal cell disappeared",
            ))?;
        if old.mass_kg.to_bits() != a.mass_kg.to_bits()
            || old.mass_kg.to_bits() != b.mass_kg.to_bits()
            || old.material_id != a.material_id
            || old.material_id != b.material_id
        {
            return Err(CavityHeatMarchError::Invalid(
                "adaptive thermal material inventory changed",
            ));
        }
        compare(
            old.specific_h_j_kg,
            a.specific_h_j_kg,
            b.specific_h_j_kg,
            controls.absolute_specific_enthalpy_tolerance_j_kg,
            rel,
            &mut error.specific_enthalpy_j_kg,
            &mut error.normalized,
        )?;
        compare(
            thermal
                .temperature_k(id)
                .map_err(CavityHeatMarchError::Thermal)?,
            coarse
                .thermal
                .temperature_k(id)
                .map_err(CavityHeatMarchError::Thermal)?,
            fine.thermal
                .temperature_k(id)
                .map_err(CavityHeatMarchError::Thermal)?,
            controls.absolute_temperature_tolerance_k,
            rel,
            &mut error.temperature_k,
            &mut error.normalized,
        )?;
    }
    for i in 0..before.cavities().len() {
        checkpoint(cx)?;
        let old = &before.cavities()[i].gas;
        let a = &coarse.system.cavities()[i].gas;
        let b = &fine.system.cavities()[i].gas;
        compare(
            old.temperature_k(),
            a.temperature_k(),
            b.temperature_k(),
            controls.absolute_temperature_tolerance_k,
            rel,
            &mut error.temperature_k,
            &mut error.normalized,
        )?;
        compare(
            old.volume_m3(),
            a.volume_m3(),
            b.volume_m3(),
            controls.absolute_gas_volume_tolerance_m3,
            rel,
            &mut error.gas_volume_m3,
            &mut error.normalized,
        )?;
        // Direct native caloric subtraction also retains PR density departures.
        let delta = a
            .stored_energy_change_to(b, || cx.checkpoint().is_err())
            .map_err(CavityHeatMarchError::Gas)?;
        let c = &coarse.balance.gases[i];
        let h = &half.balance.gases[i];
        let f = &fine.balance.gases[i];
        let fine_storage = (Dd::from_f64(h.total_energy_change_j)
            + Dd::from_f64(f.total_energy_change_j))
        .to_f64();
        let scale = c.total_energy_change_j.abs().max(fine_storage.abs());
        error.energy_j = error.energy_j.max(finite(delta)?.abs());
        error.normalized = error.normalized.max(finite(
            delta.abs() / finite(controls.absolute_energy_tolerance_j + rel * scale)?,
        )?);
        for (coarse, first, second) in [
            (
                c.heat_into_solid_j,
                h.heat_into_solid_j,
                f.heat_into_solid_j,
            ),
            (
                c.rigid_heat_stage_energy_change_j,
                h.rigid_heat_stage_energy_change_j,
                f.rigid_heat_stage_energy_change_j,
            ),
        ] {
            compare(
                0.0,
                coarse,
                (Dd::from_f64(first) + Dd::from_f64(second)).to_f64(),
                controls.absolute_energy_tolerance_j,
                rel,
                &mut error.energy_j,
                &mut error.normalized,
            )?;
        }
    }
    for ((coarse, half), fine) in energy_channels(&coarse.balance)
        .into_iter()
        .zip(energy_channels(&half.balance))
        .zip(energy_channels(&fine.balance))
        .take(13)
    {
        compare(
            0.0,
            coarse,
            (Dd::from_f64(half) + Dd::from_f64(fine)).to_f64(),
            controls.absolute_energy_tolerance_j,
            rel,
            &mut error.energy_j,
            &mut error.normalized,
        )?;
    }
    Ok(error)
}

impl<'a, G: ThermalCavityGas> ThermalPorousGasSystem<'a, G> {
    /// Resolve a complete heat/compression stroke in physical time.
    ///
    /// Prescribed coordinates move linearly from their ACTUAL accepted values
    /// to `prescribed_m` at `target_s`. Nodal/body loads and optional pressure
    /// are held throughout the stroke. Each external gas heat is a whole-stroke
    /// integral supplied at constant power; `initial_heat_j` is a whole-stroke
    /// Newton guess, scaled to each trial interval. These are explicit loading
    /// assumptions, not an interpolation of accepted plastic or gas states.
    ///
    /// Original fixed cell inventories, transport, observer transitions, source
    /// conductivity/convection, retained inertia and phase admission remain
    /// native. Cell birth/remap requires a separate admitted stroke boundary.
    /// Observers must clone independent history; their native begin/finish hooks
    /// run for each physical solve. Coarse and rejected histories are discarded.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        feedback: FiniteThermoplasticConfig,
        method: CouplingMethod,
        external_gas_heat_j: &[f64],
        initial_heat_j: &[f64],
        heat: CavityHeatSolve,
        controls: AdaptiveCavityHeatSolve,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveCavityHeatCandidate<'a, 't, O, G>, CavityHeatMarchError<O::Error>> {
        checkpoint(cx)?;
        controls.validate(heat)?;
        let start_s = self.time_s();
        let duration = finite(target_s - start_s)?;
        let count = self.cavities().len();
        if duration <= 0.0
            || thermal.time_s().to_bits() != start_s.to_bits()
            || external_gas_heat_j.len() != count
            || initial_heat_j.len() != count
            || external_gas_heat_j
                .iter()
                .chain(initial_heat_j)
                .any(|x| !x.is_finite())
        {
            return Err(CavityHeatMarchError::Invalid(
                "forward stroke and one finite heat/guess per gas required",
            ));
        }
        let initial_displacement = self.solid().mechanics().displacement_m().as_flattened();
        let mut used = BTreeSet::new();
        for &(dof, value) in prescribed_m {
            checkpoint(cx)?;
            if dof >= initial_displacement.len() || !value.is_finite() || !used.insert(dof) {
                return Err(CavityHeatMarchError::Invalid(
                    "finite unique original prescribed coordinates required",
                ));
            }
        }
        let mut accepted: Option<CavityHeatCandidate<'a, 't, O, G>> = None;
        let mut totals = Totals::new(count);
        let mut accepted_intervals = 0usize;
        let mut rejected_intervals = 0usize;
        let mut advance_calls = 0usize;
        let mut heat_trials = 0usize;
        let mut step_s = controls.initial_step_s;
        let mut maximum_error = CavityHeatLocalError::default();
        let mut time_s = start_s;
        while time_s < target_s {
            checkpoint(cx)?;
            let intervals = accepted_intervals + rejected_intervals;
            if intervals >= controls.maximum_intervals {
                return Err(CavityHeatMarchError::Budget {
                    time_s,
                    intervals,
                    heat_trials,
                });
            }
            let remaining = finite(target_s - time_s)?;
            let end_s = if step_s >= remaining {
                target_s
            } else {
                finite(time_s + step_s)?
            };
            let actual_step = finite(end_s - time_s)?;
            let mid_s = finite(time_s + 0.5 * actual_step)?;
            if !(time_s < mid_s && mid_s < end_s) {
                return Err(CavityHeatMarchError::Invalid(
                    "adaptive half interval is not representable",
                ));
            }
            let (old_system, old_thermal, old_observer) =
                accepted.as_ref().map_or((self, thermal, observer), |x| {
                    (&x.system, &x.thermal, &x.observer)
                });
            let mut coarse_heat = Vec::with_capacity(count);
            let mut first_heat = Vec::with_capacity(count);
            let mut second_heat = Vec::with_capacity(count);
            for (i, &total) in external_gas_heat_j.iter().enumerate() {
                let q = if end_s.to_bits() == target_s.to_bits() {
                    (Dd::from_f64(total) - totals.gas_external[i]).to_f64()
                } else {
                    total * (actual_step / duration)
                };
                let first = finite(q * ((mid_s - time_s) / actual_step))?;
                coarse_heat.push(finite(q)?);
                first_heat.push(first);
                second_heat.push(finite((Dd::from_f64(q) - Dd::from_f64(first)).to_f64())?);
            }
            let mut advance = |system: &Self,
                               thermal: &DepositionSimulation<'t>,
                               observer: &O,
                               end: f64,
                               q: &[f64]| {
                checkpoint(cx)?;
                let left = controls.maximum_heat_trials - heat_trials;
                if left == 0 {
                    return Err(CavityHeatMarchError::Budget {
                        time_s,
                        intervals,
                        heat_trials,
                    });
                }
                let fraction = finite((end - start_s) / duration)?;
                let mut prescribed = Vec::with_capacity(prescribed_m.len());
                for &(dof, target) in prescribed_m {
                    let value = if end.to_bits() == target_s.to_bits() {
                        target
                    } else {
                        (Dd::from_f64(initial_displacement[dof]) * Dd::from_f64(1.0 - fraction)
                            + Dd::from_f64(target) * Dd::from_f64(fraction))
                        .to_f64()
                    };
                    prescribed.push((dof, finite(value)?));
                }
                let fraction = finite((end - system.time_s()) / duration)?;
                let guess = initial_heat_j
                    .iter()
                    .map(|q| finite(q * fraction))
                    .collect::<Result<Vec<_>, CavityHeatMarchError<O::Error>>>()?;
                advance_calls =
                    advance_calls
                        .checked_add(1)
                        .ok_or(CavityHeatMarchError::Invalid(
                            "native advance count overflow",
                        ))?;
                let next = system
                    .advance(
                        thermal,
                        observer,
                        end,
                        &prescribed,
                        nodal_forces_n,
                        body_force_n_m3,
                        pressure,
                        settings,
                        mechanical_controls,
                        thermal_config,
                        feedback,
                        method,
                        q,
                        &guess,
                        CavityHeatSolve {
                            maximum_trials: heat.maximum_trials.min(left),
                            ..heat
                        },
                        cx,
                    )
                    .map_err(|source| CavityHeatMarchError::Native {
                        start_s: system.time_s(),
                        target_s: end,
                        source,
                    })?;
                heat_trials += next.trials;
                Ok(next)
            };
            let coarse = advance(old_system, old_thermal, old_observer, end_s, &coarse_heat)?;
            let half = advance(old_system, old_thermal, old_observer, mid_s, &first_heat)?;
            let fine = advance(
                &half.system,
                &half.thermal,
                &half.observer,
                end_s,
                &second_heat,
            )?;
            let error = estimate(old_system, old_thermal, &coarse, &half, &fine, controls, cx)?;
            if error.normalized > 1.0 {
                rejected_intervals += 1;
                if actual_step <= controls.minimum_step_s {
                    return Err(CavityHeatMarchError::Accuracy {
                        time_s,
                        step_s: actual_step,
                        error,
                    });
                }
                step_s = (actual_step * (0.8 / error.normalized.sqrt()).clamp(0.1, 0.5))
                    .max(controls.minimum_step_s);
                continue;
            }
            totals.add(&half.balance)?;
            totals.add(&fine.balance)?;
            if totals.energy[13].to_f64() > controls.window_energy_tolerance_j
                || totals.energy[15].to_f64() > controls.window_energy_tolerance_j
            {
                return Err(CavityHeatMarchError::Invalid(
                    "accepted stroke exceeds whole-window energy gate",
                ));
            }
            accepted_intervals += 1;
            for (maximum, value) in [
                (&mut maximum_error.displacement_m, error.displacement_m),
                (&mut maximum_error.velocity_mps, error.velocity_mps),
                (&mut maximum_error.temperature_k, error.temperature_k),
                (
                    &mut maximum_error.specific_enthalpy_j_kg,
                    error.specific_enthalpy_j_kg,
                ),
                (&mut maximum_error.gas_volume_m3, error.gas_volume_m3),
                (&mut maximum_error.plastic_state, error.plastic_state),
                (&mut maximum_error.energy_j, error.energy_j),
                (&mut maximum_error.normalized, error.normalized),
            ] {
                *maximum = maximum.max(value);
            }
            let factor = if error.normalized == 0.0 {
                2.0
            } else {
                (0.9 / error.normalized.sqrt()).clamp(0.5, 2.0)
            };
            step_s = (actual_step * factor).clamp(controls.minimum_step_s, controls.maximum_step_s);
            time_s = end_s;
            accepted = Some(fine);
        }
        let endpoint = accepted.ok_or(CavityHeatMarchError::Invalid("empty adaptive stroke"))?;
        let mut storage_defect = Dd::ZERO;
        for (i, entry) in self.cavities().iter().enumerate() {
            checkpoint(cx)?;
            let direct = entry
                .gas
                .stored_energy_change_to(&endpoint.system.cavities()[i].gas, || {
                    cx.checkpoint().is_err()
                })
                .map_err(CavityHeatMarchError::Gas)?;
            storage_defect = storage_defect
                + Dd::from_f64(
                    (Dd::from_f64(direct) - totals.gas_storage[i])
                        .to_f64()
                        .abs(),
                );
            let heat_defect =
                finite((totals.gas_external[i] - Dd::from_f64(external_gas_heat_j[i])).to_f64())?;
            if heat_defect.abs() > controls.window_energy_tolerance_j {
                return Err(CavityHeatMarchError::Invalid(
                    "whole-stroke external gas heat did not close",
                ));
            }
        }
        let storage_defect = finite(storage_defect.to_f64())?;
        if storage_defect > controls.window_energy_tolerance_j {
            return Err(CavityHeatMarchError::Invalid(
                "whole-stroke original gas storage did not close",
            ));
        }
        let energy = totals.finish(storage_defect)?;
        checkpoint(cx)?;
        Ok(AdaptiveCavityHeatCandidate {
            endpoint,
            report: CavityHeatMarchReport {
                start_s,
                end_s: target_s,
                accepted_intervals,
                rejected_intervals,
                advance_calls,
                heat_trials,
                accepted_thermal_substeps: totals.thermal_substeps,
                next_step_s: step_s,
                maximum_accepted_error: maximum_error,
                energy,
            },
        })
    }
}
