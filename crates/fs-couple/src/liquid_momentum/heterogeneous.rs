//! Heterogeneous directional liquid/body momentum exchange.
//!
//! Each original wet cell has three nonnegative authored relaxation rates and
//! a proper principal frame in world axes. The law is
//! du_i/dt = Gamma_i (V + omega cross r_i - u_i), with equal and opposite
//! force/couple on the original finite-inertia body. These rates are not a
//! guessed viscosity, impermeable-wall condition or displaced-fluid correction.
//!
//! Every local cell/body exchange is an exact kinetic-metric exponential,
//! including when directional resistance and off-center mobility do not commute.
//! A forward-half/reverse-half composition is second order for the simultaneous
//! frozen source under subdivision. It is not the exact heterogeneous global
//! exponential: callers choose and refine the explicit substep count. Exact
//! local contraction prevents stiff-step slip reversal and energy creation,
//! but does not remove the composition error of a coarse multi-cell step.
//!
//! All cell positions, source axes, body pose and clocks remain fixed. One
//! compact native transaction publishes only after final linear/angular
//! momentum, native work and independently integrated source dissipation agree.
//! The positive loss is reported; no heat partition or enthalpy update is made.

use super::*;
use fs_mbd::UnitQuaternion;

mod pairwise;

/// Authored directional relaxation on one selected ORIGINAL liquid cell.
/// Supply one entry per selected index, in the same canonical index order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrincipalLiquidDrag3 {
    /// Principal inverse times [1/s]. Exact zeros retain free-slip directions.
    pub rates_per_s: [f64; 3],
    /// Proper rotation from the resistance principal frame to world axes.
    /// The frame stays fixed for this frozen-position source operation.
    pub axes_to_world: UnitQuaternion,
}

/// Native physical gates plus an explicit bounded source-composition schedule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeterogeneousRigidDragControls {
    pub exchange: LiquidRigidDragControls,
    /// Positive number of symmetric composition substeps. Refinement controls
    /// temporal composition error independently of momentum/work closure.
    pub substeps: usize,
    /// At most (2 * selected_cells - 1) * substeps exact local exchanges.
    /// Checked before either the native stencil or response cache is allocated.
    pub maximum_pair_exchanges: usize,
}

/// Both accepted original owners and the physical loss of the composed source.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct HeterogeneousRigidDragReport {
    pub source_duration_s: f64,
    pub substeps: usize,
    /// Actual active local exchanges; zero-rate/no-duration sources do no work.
    pub pair_exchanges: usize,
    /// Uses independently integrated nonnegative local source dissipation.
    /// Its sum is checked against both original owners' actual kinetic change.
    pub exchange: LiquidRigidExchangeReport,
}

fn kernel_error(error: pairwise::Error) -> LiquidMomentumError {
    match error {
        pairwise::Error::Invalid(reason) => invalid(reason),
        pairwise::Error::Cancelled => LiquidMomentumError::Cancelled,
    }
}

impl HeterogeneousRigidDragControls {
    pub(super) fn admit(
        self,
        selected_cells: usize,
        laws: &[PrincipalLiquidDrag3],
        source_duration_s: f64,
    ) -> Result<(), LiquidMomentumError> {
        self.exchange.admit(selected_cells)?;
        if laws.len() != selected_cells
            || !source_duration_s.is_finite()
            || source_duration_s < 0.0
            || laws
                .iter()
                .any(|law| law.rates_per_s.iter().any(|v| !v.is_finite() || *v < 0.0))
        {
            return Err(invalid(
                "heterogeneous drag needs matching nonnegative rates and finite source duration",
            ));
        }
        pairwise::work_bound(selected_cells, self.substeps, self.maximum_pair_exchanges)
            .map_err(kernel_error)?;
        Ok(())
    }
}

/// Apply heterogeneous principal drag to original native liquid populations.
///
/// Partial cells use their retained liquid mass. The finite body recoils with
/// its complete original inertia tensor, through the same native momentum
/// owner and final conservation/work gates as the uniform exact source.
/// Source coefficients and small eigensystems are prepared once per cell.
/// The symmetric schedule costs O(stencil_cells * substeps), stores O(stencil)
/// data, and performs no full-grid clone or dense global exponential.
///
/// Cancellation, coefficient/domain refusal, unresolved modes or a failed
/// actual native balance publishes neither owner. This is a source substep:
/// pressure, transport, rigid drift/free rotation and heat remain separate.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_drag_native(
    liquid: &mut FreeSurface3,
    body: &mut RigidBodyState,
    properties: MassProperties,
    indices: &[usize],
    laws: &[PrincipalLiquidDrag3],
    units: LiquidMomentumUnits3,
    source_duration_s: f64,
    controls: HeterogeneousRigidDragControls,
    gate: &CancelGate,
) -> Result<HeterogeneousRigidDragReport, LiquidMomentumError> {
    poll(gate)?;
    controls.admit(indices.len(), laws, source_duration_s)?;
    let scales = units.admit()?;
    let stage =
        liquid.prepare_impulse_stencil(indices, controls.exchange.native(units, scales)?, gate)?;
    apply_staged(
        stage,
        body,
        properties,
        laws,
        units,
        source_duration_s,
        controls,
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn apply_staged(
    stage: FreeSurfaceImpulseStencil3<'_>,
    body: &mut RigidBodyState,
    properties: MassProperties,
    laws: &[PrincipalLiquidDrag3],
    units: LiquidMomentumUnits3,
    source_duration_s: f64,
    controls: HeterogeneousRigidDragControls,
    gate: &CancelGate,
) -> Result<HeterogeneousRigidDragReport, LiquidMomentumError> {
    poll(gate)?;
    controls.admit(stage.cells().len(), laws, source_duration_s)?;
    let (cells, basis) = prepare_source(&stage, body, properties, units, controls.exchange, gate)?;
    let mut source = Vec::new();
    source
        .try_reserve_exact(cells.len())
        .map_err(|_| invalid("heterogeneous liquid/body source allocation refused"))?;
    for (index, (cell, law)) in cells.iter().zip(laws).enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let rotational_velocity = cross(basis.omega, cell.arm);
        let slip = core::array::from_fn(|axis| {
            (dd(a(cell.velocity)[axis]) - dd(a(basis.velocity)[axis]) - rotational_velocity[axis])
                .to_f64()
        });
        source.push(pairwise::Cell {
            mass: cell.mass,
            slip,
            body_columns: basis.columns(cell.arm).map(a),
            principal_rates_per_s: law.rates_per_s,
            principal_axes: [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ]
            .map(|axis| a(law.axes_to_world.rotate_body_to_world(axis))),
        });
    }
    let result = pairwise::solve(
        &source,
        source_duration_s,
        controls.substeps,
        controls.maximum_pair_exchanges,
        || gate.is_requested(),
    )
    .map_err(kernel_error)?;
    let (body_linear, body_angular) = basis.physical_impulse(result.body_increment);
    let mut impulses = Vec::new();
    impulses
        .try_reserve_exact(result.impulses.len())
        .map_err(|_| invalid("heterogeneous liquid impulse allocation refused"))?;
    for impulse in result.impulses {
        impulses.push(v(impulse));
    }
    let exchange = finish_source(
        stage,
        body,
        properties,
        units,
        controls.exchange,
        &cells,
        basis,
        relaxation::Plan {
            impulses,
            body_linear,
            body_angular,
        },
        Some(result.dissipation_j),
        gate,
    )?;
    Ok(HeterogeneousRigidDragReport {
        source_duration_s,
        substeps: controls.substeps,
        pair_exchanges: result.pair_exchanges,
        exchange,
    })
}
