//! Conducting finite wall pairs on the actual native open-flow heat carrier.
//!
//! A pair is selected by its normal Cartesian axis. The graph borrows the
//! original lower/upper wall arrays; there is no copied solid temperature
//! field and no second heat deposition. Different axis pairs can be installed
//! independently. Connections across pairs or interior solid depth require a
//! different physical discretization and are not inferred from coordinates.

use super::{D3Q7, ThermalOpenError3, ThermalOpenFace3, ThermalOpenField3, ThermalWall, invalid};
use crate::thermal::boundary::conduction::{
    WallConduction, WallConductionControls, WallConductionLink,
};

/// Inverse of the carrier's x-fast surface index, for any normal axis.
pub(super) fn boundary_cell(dims: [usize; 3], axis: usize, side: usize, link: usize) -> usize {
    let stride: usize = dims[..axis].iter().product();
    (link / stride) * (stride * dims[axis]) + link % stride + side * (dims[axis] - 1) * stride
}

impl ThermalOpenField3 {
    /// Prepare implicit conduction between existing finite cells on one wall
    /// pair. `axis` is 0/1/2 for X/Y/Z; graph side 0 is the minimum face and
    /// side 1 the maximum. Link indices use the unchanged x-fast surface order.
    ///
    /// Every endpoint must face an active fluid cell. Blocked slots in the
    /// indexed wall arrays cannot consume sources or act as hidden solid nodes.
    /// Sparse-factor and work budgets are explicit for EACH installed pair,
    /// additional to the carrier's original thermal storage budget. Preparing
    /// a replacement briefly retains its previous factor until success.
    ///
    /// The existing wall cells remain the sole physical temperature owner.
    /// A refused graph leaves the previous plan and all physical state intact;
    /// successful preparation changes neither temperature nor accepted time.
    pub fn set_wall_conduction(
        &mut self,
        axis: usize,
        links: &[WallConductionLink],
        controls: WallConductionControls,
    ) -> Result<(), ThermalOpenError3> {
        if axis >= 3
            || (0..2).any(|side| {
                !matches!(
                    self.boundaries.faces[2 * axis + side],
                    ThermalOpenFace3::Wall(_)
                )
            })
        {
            return Err(invalid(
                "wall conduction requires one original stationary wall pair",
            ));
        }
        if links.len() > controls.max_links {
            return Err(invalid(
                "wall conduction link count exceeds its explicit graph limit",
            ));
        }
        let plane = self.populations.len() / self.dims[axis];
        for link in links {
            for endpoint in [link.a, link.b] {
                if endpoint.side >= 2
                    || endpoint.link >= plane
                    || self.solid[boundary_cell(self.dims, axis, endpoint.side, endpoint.link)]
                {
                    return Err(invalid(
                        "wall conduction endpoint must face an active fluid cell",
                    ));
                }
            }
        }
        let walls: &[Vec<ThermalWall>; 2] = self.walls[2 * axis..2 * axis + 2]
            .try_into()
            .expect("admitted wall pair");
        let (positive, _) = D3Q7.direction_pair(axis);
        let prepared =
            WallConduction::prepare(walls, plane, D3Q7.weights[positive], links, controls)
                .map_err(ThermalOpenError3::Conduction)?;
        self.conduction[axis] = Some(prepared);
        Ok(())
    }

    /// Inspect the immutable graph preparation; it owns no physical wall data.
    pub fn wall_conduction(&self, axis: usize) -> Option<&WallConduction> {
        self.conduction.get(axis).and_then(Option::as_ref)
    }

    /// Return this pair to its original independent finite-wall update.
    /// Current fluid and solid temperatures, histories and time remain intact.
    pub fn clear_wall_conduction(&mut self, axis: usize) -> Result<(), ThermalOpenError3> {
        let graph = self
            .conduction
            .get_mut(axis)
            .ok_or_else(|| invalid("invalid wall axis"))?;
        *graph = None;
        Ok(())
    }

    /// Prepare scratch once per installed pair, after the scalar collision.
    /// Failed/cancelled outer steps need no graph rollback: every later attempt
    /// rebuilds these reflections from the accepted wall copies and new post g.
    pub(super) fn prepare_wall_conduction(&mut self) -> Result<(), ThermalOpenError3> {
        let dims = self.dims;
        let post = &self.post;
        for axis in 0..3 {
            let Some(conduction) = &mut self.conduction[axis] else {
                continue;
            };
            let walls: &[Vec<ThermalWall>; 2] = self.trial_walls[2 * axis..2 * axis + 2]
                .try_into()
                .expect("admitted wall pair");
            let (positive, negative) = D3Q7.direction_pair(axis);
            conduction
                .prepare_reflections(walls, |side, link| {
                    let cell = boundary_cell(dims, axis, side, link);
                    post[7 * cell + if side == 0 { negative } else { positive }]
                })
                .map_err(ThermalOpenError3::Conduction)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
