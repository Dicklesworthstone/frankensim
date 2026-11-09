//! Evolving porous laws on the original open-flow endpoint transaction.

use super::*;
use crate::lattice3::sources::resistance::ResistanceEndpointExchange3;

/// The original endpoint force/work account plus actual reservoir transport.
///
/// `resistance` retains body and skeleton impulses, independent endpoint-power
/// and impulse-work quadratures, every wall link and periodic angular transport.
/// `open` measures the original reservoir replacement, not an inferred residual.
/// Changing material or skeleton motion can make the temporal work remainder
/// signed. Only endpoint relative dissipation is nonnegative; the remainder is
/// never heat. Geometry and prescribed skeleton inertia remain outside this law.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[must_use]
pub struct ResistanceOpenEndpointExchange3 {
    pub resistance: ResistanceEndpointExchange3,
    pub open: OpenBoundaryMomentumTransport3,
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ResistanceOpenGrid3<L, Q> {
    /// Evolving-law receipt from the last accepted interval, or None when that
    /// interval held its law with `step_*`. Failures retain the accepted receipt.
    pub const fn last_endpoint_exchange(&self) -> Option<ResistanceOpenEndpointExchange3> {
        self.last_endpoint_exchange
    }

    /// Advance to a complete law field at the NEXT physical endpoint.
    ///
    /// The accepted force still drives collision. The original native pull and
    /// reservoir reconstruction resolve each next force from its own raw moments
    /// and law. Complete face traces use their prescribed velocity and next law;
    /// their interior donors also use their own next law. No same-time rebinding,
    /// population reset, frozen-force endpoint or cloned flow step is used.
    ///
    /// Supply one admitted lattice-unit law per cell, including unused solids.
    /// Active resistance on a mixed wall/open seam still refuses. An unchanged
    /// law field produces the same native populations and force history as
    /// `step_serial`, with the distinct endpoint-work receipt explicitly named.
    /// `maximum_absolute_work_remainder` bounds the sum of absolute per-cell
    /// relative temporal work remainders in lattice energy units. Resolve rapid
    /// changes in time; a coefficient drop can inject impulse work even though
    /// the endpoint drag powers are dissipative. Nothing is clipped to pass.
    /// All laws, populations, forces, drag, receipts and time publish together.
    pub fn advance_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
    ) -> Result<ResistanceOpenEndpointExchange3, GridError3<L::Error>> {
        self.advance_with_poll::<fs_exec::TilePool>(
            None,
            collisions,
            next_laws,
            maximum_absolute_work_remainder,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Identical constitutive and reservoir ownership with the native executor's
    /// fixed tile arithmetic and canonical force/work reduction order.
    pub fn advance_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
    ) -> Result<ResistanceOpenEndpointExchange3, GridError3<L::Error>> {
        self.advance_with_poll(
            Some((pool, gate)),
            collisions,
            next_laws,
            maximum_absolute_work_remainder,
            &mut || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
    }

    fn advance_with_poll<P: KernelRunner>(
        &mut self,
        pool: Option<(&P, &CancelGate)>,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<ResistanceOpenEndpointExchange3, GridError3<L::Error>> {
        let (_, receipt, step) = self.stage_with_laws(
            pool,
            collisions,
            Some((&next_laws, maximum_absolute_work_remainder)),
            poll,
        )?;
        poll()?;
        self.publish_endpoint_step(step, receipt, next_laws);
        Ok(receipt)
    }

    pub(super) fn admit_next_laws(
        &self,
        laws: &[DarcyForchheimer3],
        maximum_absolute_work_remainder: f64,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(), GridError3<L::Error>> {
        if !maximum_absolute_work_remainder.is_finite() || maximum_absolute_work_remainder < 0.0 {
            return Err(GridError3::Input("invalid resistance temporal work budget"));
        }
        source_endpoint::admit_laws::<L, Q>(laws, self.laws.len(), &mut *poll)?;
        let grid = &self.flow.inner.grid;
        let data = StepData::<L, Q> {
            dims: grid.dims,
            periodic: grid.periodic,
            solid: &grid.solid,
            collision: grid.collision,
            collisions: None,
            force: grid.force,
            forces: None,
            budget: grid.budget,
            marker: PhantomData,
        };
        let plan = self.flow.inner.plan;
        for (cell, law) in laws.iter().enumerate() {
            if cell.is_multiple_of(64) {
                poll()?;
            }
            if !grid.solid[cell]
                && (law.linear.principal != [0.0; 3] || law.quadratic != 0.0)
                && plan.locate(cell).is_some()
            {
                let (wall, open) = plan.masks(&data, cell);
                if wall != 0 && open != 0 {
                    return Err(GridError3::MaterialField {
                        cell,
                        reason: "active resistance requires a complete open-face trace; wall/open seams are unsupported",
                    });
                }
            }
        }
        poll()
    }

    pub(super) fn publish_endpoint_step(
        &mut self,
        step: usize,
        receipt: ResistanceOpenEndpointExchange3,
        next_laws: Vec<DarcyForchheimer3>,
    ) {
        self.publish_fields(step);
        self.laws = next_laws;
        self.last_exchange = None;
        self.last_endpoint_exchange = Some(receipt);
    }
}
