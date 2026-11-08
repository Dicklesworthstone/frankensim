//! Original wall sources surrounding the same pressure/inertia transaction.

use super::*;
use crate::body_force::ForceChange;
use crate::dynamics::nonideal::RealGasWallReport;

impl RealGasEuler1d {
    /// Advance wall heat/drag, real-gas pressure work and a finite right piston.
    ///
    /// A half wall-source interval precedes the original pressure/inertia solve,
    /// so heating can drive the piston during this step. The second half uses
    /// the actual transported density and moved cell volumes. Each finite wall
    /// carries its first-half inventory into the second half; reservoir heat,
    /// prescribed sliding-wall work and internal slip heating remain distinct.
    /// Finite material walls retain their original H/U and pressure-work convention.
    /// Shared regions and their conduction links use each single original wall
    /// inventory through both halves, within one shared source-work allowance.
    ///
    /// The first source endpoint is computed once per attempted duration and
    /// reused by every midpoint-speed proposal. A shortened or refused attempt
    /// restarts both sources and pressure evolution from accepted gas and walls.
    /// Gas, thermodynamic caches, geometry, piston, wall inventories and reports
    /// publish together. The wall binding follows only accepted mesh motion.
    /// Source and geometry buffers allocate once and are reused.
    ///
    /// `PistonStep::gas` covers the complete interval, including its wall-source
    /// accounts. The mechanical pressure/kinetic check remains distinct from
    /// the complete gas/piston/finite-wall balance. Backward-Euler local heat
    /// makes this split composition first order, despite its symmetric ordering.
    /// Conductances, finite wall amounts and sliding velocities are authored;
    /// neither changing contact area nor a thermal piston mass is inferred.
    /// The original closed constant-area, frozen-composition branch restrictions
    /// remain. Chemistry, distributed forcing and phase flash are not composed.
    ///
    /// # Errors
    /// Invalid bindings/controls, source or full-interval CFL refusal, failed
    /// mechanical/storage balances, exhausted retries and cancellation leave
    /// every accepted gas, piston and wall owner unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston_and_walls(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        controls: PistonControl,
        walls: &mut RealGasWallOperator<'_>,
        wall_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        self.step_with_piston_suspension(
            piston,
            maximum_duration,
            gas_controls,
            load,
            None,
            controls,
            Some((walls, wall_controls)),
            None,
            cancelled,
        )
        .map(|step| step.piston)
    }

    /// Retain spring storage and viscous loss in the same thermal piston step.
    ///
    /// This uses [`Self::step_with_piston_and_walls`]'s original heat/pressure/
    /// heat composition. Suspension damping remains a separately reported loss;
    /// it is not silently returned as gas or wall heat. Authored wall drag has
    /// its own support impulse, sliding work and internal slip heating.
    ///
    /// # Errors
    /// The same complete-transaction refusals apply as for the free piston.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_suspended_piston_and_walls(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        spring: PistonSuspension,
        controls: PistonControl,
        walls: &mut RealGasWallOperator<'_>,
        wall_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep, PistonError> {
        self.step_with_piston_suspension(
            piston,
            maximum_duration,
            gas_controls,
            load,
            Some(spring),
            controls,
            Some((walls, wall_controls)),
            None,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn finish(
    step: &mut SuspendedPistonStep,
    before: Totals,
    after: Totals,
    stages: Option<[RealGasWallReport; 2]>,
    body: ForceChange,
    piston: Piston,
    load: PistonLoad,
    area: f64,
    controls: PistonControl,
) -> Result<Option<RealGasWallReport>, PistonError> {
    let report = &mut step.piston;
    let mut thermal = if let Some(stages) = stages {
        let mut thermal = stages[0];
        thermal.elapsed = report.gas.elapsed;
        for (total, second) in [
            (&mut thermal.wall_impulse_ns, stages[1].wall_impulse_ns),
            (&mut thermal.heat_to_gas_j, stages[1].heat_to_gas_j),
            (
                &mut thermal.wall_storage_change_j,
                stages[1].wall_storage_change_j,
            ),
            (
                &mut thermal.external_heat_to_gas_j,
                stages[1].external_heat_to_gas_j,
            ),
            (
                &mut thermal.mechanical_work_to_gas_j,
                stages[1].mechanical_work_to_gas_j,
            ),
            (&mut thermal.drag_heating_j, stages[1].drag_heating_j),
            (
                &mut thermal.source_storage_residual_j,
                stages[1].source_storage_residual_j,
            ),
        ] {
            *total += second;
            if !total.is_finite() {
                return Err(PistonError::Invalid(
                    "thermal piston source account overflow",
                ));
            }
        }
        thermal.maximum_thermal_residual_j = thermal
            .maximum_thermal_residual_j
            .max(stages[1].maximum_thermal_residual_j);
        thermal.maximum_iterations = thermal.maximum_iterations.max(stages[1].maximum_iterations);
        Some(thermal)
    } else {
        None
    };
    let (wall_impulse, heat, wall_storage, external_heat, wall_work) =
        thermal
            .as_ref()
            .map_or((0.0, 0.0, 0.0, 0.0, 0.0), |report| {
                (
                    report.wall_impulse_ns,
                    report.heat_to_gas_j,
                    report.wall_storage_change_j,
                    report.external_heat_to_gas_j,
                    report.mechanical_work_to_gas_j,
                )
            });
    report.gas.before = before;
    report.gas.after = after;
    report.gas.wall_momentum_impulse = -wall_impulse;
    report.gas.wall_heat = heat;
    report.gas.wall_work = wall_work;
    report.gas.source_change = body;
    let source = [
        0.0,
        sum([-wall_impulse, body.impulse]),
        sum([heat, wall_work, body.work]),
    ];
    report.gas.conservation_residual = Totals::from_array(core::array::from_fn(|i| {
        sum([
            after.array()[i],
            -before.array()[i],
            -report.gas.boundary_change.array()[i],
            -source[i],
        ])
    }));
    if let Some(thermal) = &mut thermal {
        thermal.gas_wall_storage_residual_j = sum([
            after.energy,
            -before.energy,
            -report.gas.boundary_change.energy,
            wall_storage,
            -external_heat,
            -wall_work,
            -body.work,
        ]);
    }
    report.combined_momentum_residual_n_s = sum([
        after.momentum,
        -before.momentum,
        piston.mass_kg() * (report.velocity_m_s[1] - report.velocity_m_s[0]),
        report.fixed_wall_impulse_n_s,
        wall_impulse,
        -report.gas.elapsed.value() * (load.external_force_n - area * load.exterior_pressure_pa),
        -step.suspension_impulse_n_s,
        -body.impulse,
    ]);
    report.total_energy_residual_j = sum([
        after.energy,
        -before.energy,
        report.kinetic_energy_change_j,
        report.exterior_pressure_work_j,
        -report.external_force_work_j,
        wall_storage,
        -external_heat,
        -wall_work,
        -body.work,
    ]);
    step.total_energy_residual_j = sum([
        after.energy,
        -before.energy,
        report.kinetic_energy_change_j,
        step.spring_energy_change_j,
        step.damping_energy_j,
        report.exterior_pressure_work_j,
        -step.applied_force_work_j,
        wall_storage,
        -external_heat,
        -wall_work,
        -body.work,
    ]);
    if report
        .gas
        .conservation_residual
        .array()
        .iter()
        .any(|value| !value.is_finite())
        || thermal
            .as_ref()
            .is_some_and(|report| !report.gas_wall_storage_residual_j.is_finite())
        || !report.combined_momentum_residual_n_s.is_finite()
        || !report.total_energy_residual_j.is_finite()
        || !step.total_energy_residual_j.is_finite()
    {
        return Err(PistonError::Invalid(
            "sourced piston complete balance overflow",
        ));
    }
    if !report.admitted(controls)
        || step.total_energy_residual_j.abs() > controls.energy_tolerance_j
    {
        return Err(PistonError::Balance {
            momentum_residual_n_s: report
                .momentum_residual_n_s
                .abs()
                .max(report.combined_momentum_residual_n_s.abs()),
            energy_residual_j: report
                .mechanical_energy_residual_j
                .abs()
                .max(report.total_energy_residual_j.abs())
                .max(step.total_energy_residual_j.abs()),
        });
    }
    Ok(thermal)
}

fn sum<const N: usize>(values: [f64; N]) -> f64 {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    total + correction
}
