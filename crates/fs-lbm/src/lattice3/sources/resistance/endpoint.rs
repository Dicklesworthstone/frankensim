//! Physical endpoint laws and explicit temporal force-work quadratures.

use super::*;
use crate::d3q19::PeriodicMomentumTransport3;

/// Two second-order work quadratures for one force/velocity pair over dt=1.
/// Neither is a complete transported fluid energy balance.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ResistanceWork3 {
    /// Endpoint-averaged force impulse dotted with mean endpoint velocity.
    pub impulse_work: f64,
    /// Trapezoidal endpoint power, `(f0.v0 + f1.v1)/2`.
    pub endpoint_work: f64,
    /// Signed difference `impulse_work - endpoint_work`, independently
    /// evaluated as `-(f1-f0).(v1-v0)/4`. This is not heat or dissipation.
    pub quadrature_remainder: f64,
}

/// One accepted step between independently prescribed porous-law endpoints.
///
/// Drag `d_i` is the force ON fluid, and relative velocity is `v_i=u_i-v_s,i`.
/// External and skeleton impulses are `(b0+b1)/2` and `-(d0+d1)/2`. Native wall
/// and periodic transport close the same momentum balance as held-law steps.
///
/// Each work entry distinguishes midpoint impulse work from trapezoidal
/// physical endpoint power. In particular, endpoint drag and skeleton work
/// sum to `-endpoint_relative_dissipation`, while their impulse works sum to
/// `-endpoint_relative_dissipation + relative_quadrature_remainder`, up to
/// floating-point precision. The signed remainder can have EITHER sign when
/// resistance or prescribed skeleton motion evolves; it is never deposited
/// as heat. Body-force work and its own temporal remainder remain separate.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct ResistanceEndpointExchange3 {
    /// Independent body-force impulse over the elapsed tick.
    pub external_impulse: [f64; 3],
    pub external_angular_impulse: [f64; 3],
    /// Opposite of the actual retained endpoint-averaged fluid drag.
    pub skeleton_impulse: [f64; 3],
    pub skeleton_angular_impulse: [f64; 3],
    /// Independent force paired with physical fluid velocity.
    pub external_work: ResistanceWork3,
    /// Fluid drag paired with physical fluid velocity.
    pub fluid_drag_work: ResistanceWork3,
    /// Opposite drag paired with each prescribed skeleton velocity.
    pub skeleton_work: ResistanceWork3,
    /// Nonnegative `.5*(v0.D0.v0 + beta0*|v0|^3 + v1.D1.v1 + beta1*|v1|^3)`.
    /// This endpoint power quadrature is NOT held-law midpoint dissipation.
    pub endpoint_relative_dissipation: f64,
    /// `-(d1-d0).(v1-v0)/4`, evaluated directly from the relative velocities.
    /// Equals the sum of the fluid-drag and skeleton work remainders.
    pub relative_quadrature_remainder: f64,
    /// Sum of ABSOLUTE per-cell relative quadrature remainders, in lattice
    /// energy units. Spatial cancellation cannot reduce the acceptance budget.
    pub absolute_relative_quadrature_remainder: f64,
    /// Every native stationary exterior and solid wall link.
    pub walls: WallMomentumExchange3,
    /// Actual periodic coordinate jumps from this step's native fluid pulls.
    pub periodic: PeriodicMomentumTransport3,
}

pub(super) struct EndpointLaws3 {
    pub(super) laws: Vec<DarcyForchheimer3>,
    pub(super) maximum_absolute_work_remainder: f64,
}

impl<L: ForcedLattice3<Q>, const Q: usize> ResistanceGrid3<L, Q> {
    /// Advance to a complete law field sampled at the NEXT physical endpoint.
    /// Collision uses the already accepted `F_n` and its retained drag. The
    /// new tensor, skeleton velocity and body force do not replace that history.
    /// The native pull solves `F_(n+1)` from the supplied law and endpoint raw
    /// moments. Laws, populations, forces, drag and clock publish together.
    ///
    /// `next_laws` has one admitted entry per cell, including unused solids,
    /// in x-fastest order. The law's coefficients and forces use lattice units.
    /// Smooth histories use trapezoidal source timing; resolve
    /// changes in time. This does not interpolate material geometry, infer a
    /// temperature/damage model, advance a skeleton body, or remap populations.
    /// `maximum_absolute_work_remainder` is a finite nonnegative bound, in
    /// lattice energy units, on the SUM of absolute per-cell relative temporal
    /// work remainders. It is checked before publication, without a unit floor
    /// or signed cancellation. Large coefficient drops can inject impulse work
    /// despite positive endpoint dissipation; resolve such histories in time.
    /// Current collision coefficients retain `step_serial`'s meaning.
    /// Any invalid law, unresolved endpoint, collision or cancellation refusal
    /// leaves the complete accepted state AND law field unchanged.
    pub fn advance_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
    ) -> Result<ResistanceEndpointExchange3, GridError3<L::Error>> {
        let endpoint = EndpointLaws3 {
            laws: next_laws,
            maximum_absolute_work_remainder,
        };
        self.step_inner::<fs_exec::TilePool>(None, collisions, Some(endpoint), || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
        .map(|(_, endpoint)| endpoint)
    }

    /// Same physical endpoint timing, native kernels and atomic publication,
    /// with fixed64-cell arithmetic independent of the worker count.
    pub fn advance_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
    ) -> Result<ResistanceEndpointExchange3, GridError3<L::Error>> {
        let endpoint = EndpointLaws3 {
            laws: next_laws,
            maximum_absolute_work_remainder,
        };
        self.step_inner(Some((pool, gate)), collisions, Some(endpoint), || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })
        .map(|(_, endpoint)| endpoint)
    }
}

pub(super) fn admit_laws<L: ForcedLattice3<Q>, const Q: usize>(
    laws: &[DarcyForchheimer3],
    count: usize,
    mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<(), GridError3<L::Error>> {
    checkpoint()?;
    if laws.len() != count {
        return Err(GridError3::Input("resistance endpoint field length"));
    }
    for (cell, &law) in laws.iter().enumerate() {
        if cell.is_multiple_of(64) {
            checkpoint()?;
        }
        law.admit().map_err(|error| GridError3::MaterialField {
            cell,
            reason: match error {
                BoundaryStepError3::EndpointForce { reason } => reason,
                _ => "invalid endpoint resistance law",
            },
        })?;
    }
    checkpoint()
}

impl ResistanceEndpointExchange3 {
    pub(super) fn add(&mut self, other: Self) {
        for a in 0..3 {
            self.external_impulse[a] += other.external_impulse[a];
            self.external_angular_impulse[a] += other.external_angular_impulse[a];
            self.skeleton_impulse[a] += other.skeleton_impulse[a];
            self.skeleton_angular_impulse[a] += other.skeleton_angular_impulse[a];
        }
        for (target, source) in [
            (&mut self.external_work, other.external_work),
            (&mut self.fluid_drag_work, other.fluid_drag_work),
            (&mut self.skeleton_work, other.skeleton_work),
        ] {
            target.impulse_work += source.impulse_work;
            target.endpoint_work += source.endpoint_work;
            target.quadrature_remainder += source.quadrature_remainder;
        }
        self.endpoint_relative_dissipation += other.endpoint_relative_dissipation;
        self.relative_quadrature_remainder += other.relative_quadrature_remainder;
        self.absolute_relative_quadrature_remainder += other.absolute_relative_quadrature_remainder;
    }

    pub(super) fn is_finite(self) -> bool {
        [
            self.external_impulse,
            self.external_angular_impulse,
            self.skeleton_impulse,
            self.skeleton_angular_impulse,
        ]
        .into_iter()
        .flatten()
        .chain(
            [self.external_work, self.fluid_drag_work, self.skeleton_work]
                .into_iter()
                .flat_map(|w| [w.impulse_work, w.endpoint_work, w.quadrature_remainder]),
        )
        .chain([
            self.endpoint_relative_dissipation,
            self.relative_quadrature_remainder,
            self.absolute_relative_quadrature_remainder,
        ])
        .all(f64::is_finite)
    }
}

pub(super) fn exchange<L: ForcedLattice3<Q>, const Q: usize>(
    initial: SourceEndpoint3,
    end: SourceEndpoint3,
    cell: usize,
    position: [usize; 3],
    origin: [f64; 3],
) -> Result<ResistanceEndpointExchange3, GridError3<L::Error>> {
    // Exact unchanged laws retain the held step's impulse/work arithmetic,
    // including its subnormal and signed-zero behavior and admission guard.
    let held = if initial.law == end.law {
        Some(held_exchange::<L, Q>(initial, end, cell, position, origin)?)
    } else {
        None
    };
    let b0 = initial.law.linear.body_force;
    let b1 = end.law.linear.body_force;
    let s0 = initial.law.linear.skeleton_velocity;
    let s1 = end.law.linear.skeleton_velocity;
    let external_impulse = held.map_or_else(|| mean(b0, b1), |r| r.external_impulse);
    let skeleton_impulse = core::array::from_fn(|a| -0.5 * initial.drag[a] - 0.5 * end.drag[a]);
    let u_mid = mean(initial.velocity, end.velocity);
    let arm = core::array::from_fn(|a| position[a] as f64 + 0.5 - origin[a]);
    let external_impulse_work =
        held.map_or_else(|| dot(external_impulse, u_mid), |r| r.external_work);
    let drag_impulse_work =
        held.map_or_else(|| -dot(skeleton_impulse, u_mid), |r| r.fluid_drag_work);
    let skeleton_impulse_work =
        held.map_or_else(|| dot(skeleton_impulse, mean(s0, s1)), |r| r.skeleton_work);
    let failure = || GridError3::MaterialField {
        cell,
        reason: "endpoint resistance work is not representable",
    };
    let (external_work, _) = work(
        b0,
        b1,
        initial.velocity,
        end.velocity,
        external_impulse_work,
    )
    .ok_or_else(failure)?;
    let (fluid_drag_work, fluid_scale) = work(
        initial.drag,
        end.drag,
        initial.velocity,
        end.velocity,
        drag_impulse_work,
    )
    .ok_or_else(failure)?;
    let (skeleton_work, skeleton_scale) = work(
        initial.drag.map(|f| -f),
        end.drag.map(|f| -f),
        s0,
        s1,
        skeleton_impulse_work,
    )
    .ok_or_else(failure)?;
    let endpoint_relative_dissipation = 0.5 * initial.law.dissipation_power(initial.velocity)
        + 0.5 * end.law.dissipation_power(end.velocity);
    let relative0 = core::array::from_fn(|a| initial.velocity[a] - s0[a]);
    let relative1 = core::array::from_fn(|a| end.velocity[a] - s1[a]);
    let relative_quadrature_remainder = -dot(
        half_difference(end.drag, initial.drag),
        half_difference(relative1, relative0),
    );
    let report = ResistanceEndpointExchange3 {
        external_impulse,
        external_angular_impulse: cross(arm, external_impulse),
        skeleton_impulse,
        skeleton_angular_impulse: cross(arm, skeleton_impulse),
        external_work,
        fluid_drag_work,
        skeleton_work,
        endpoint_relative_dissipation,
        relative_quadrature_remainder,
        absolute_relative_quadrature_remainder: relative_quadrature_remainder.abs(),
        ..ResistanceEndpointExchange3::default()
    };
    let endpoint_defect =
        fluid_drag_work.endpoint_work + skeleton_work.endpoint_work + endpoint_relative_dissipation;
    let impulse_defect =
        fluid_drag_work.impulse_work + skeleton_work.impulse_work + endpoint_relative_dissipation
            - relative_quadrature_remainder;
    let scale = fluid_scale
        + skeleton_scale
        + endpoint_relative_dissipation
        + relative_quadrature_remainder.abs();
    if !report.is_finite()
        || endpoint_relative_dissipation < 0.0
        || !scale.is_finite()
        || endpoint_defect.abs().max(impulse_defect.abs()) > 4096.0 * f64::EPSILON * scale
    {
        return Err(failure());
    }
    Ok(report)
}

fn work(
    f0: [f64; 3],
    f1: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    impulse_work: f64,
) -> Option<(ResistanceWork3, f64)> {
    let endpoint_work = dot(f0.map(|f| 0.5 * f), v0) + dot(f1.map(|f| 0.5 * f), v1);
    let quadrature_remainder = -dot(half_difference(f1, f0), half_difference(v1, v0));
    let scale = (0..3)
        .map(|a| 0.5 * f0[a].abs() * v0[a].abs() + 0.5 * f1[a].abs() * v1[a].abs())
        .sum::<f64>()
        + impulse_work.abs()
        + quadrature_remainder.abs();
    let defect = impulse_work - endpoint_work - quadrature_remainder;
    ([impulse_work, endpoint_work, quadrature_remainder, scale]
        .into_iter()
        .all(f64::is_finite)
        && defect.abs() <= 4096.0 * f64::EPSILON * scale)
        .then_some((
            ResistanceWork3 {
                impulse_work,
                endpoint_work,
                quadrature_remainder,
            },
            scale,
        ))
}

fn mean(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| 0.5 * a[i] + 0.5 * b[i])
}
fn half_difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| 0.5 * a[i] - 0.5 * b[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
