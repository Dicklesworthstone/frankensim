//! Complete source-owned volumetric relaxation on the original sparse liquid.

use super::transport::MassLocation3;
use super::{Cell3, FreeSurface3, FreeSurfaceError3, TILE_CELLS, check_cancel};
use crate::d3q19::CollisionModel3;
use fs_exec::CancelGate;
use std::collections::BTreeMap;

impl FreeSurface3 {
    /// Install source-resolved trace rates on this candidate. The field must
    /// cover all and only current wet collision cells, including empty
    /// interface carriers. No inferred rate is assigned to gas or carry.
    pub(crate) fn rebind_bulk_relaxations(
        &mut self,
        field: BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        if !matches!(self.collision_model(), CollisionModel3::CentralMomentBulk { .. }) {
            return Err(FreeSurfaceError3::Transport(
                "local bulk viscosity requires independent central-moment collision",
            ));
        }
        self.admit_bulk_relaxations(&field, gate)?;
        self.local_bulk_relaxation = Some(field);
        Ok(())
    }

    /// Select a source field atomically, including its reference collision.
    /// Every allocation and fallible field check precedes native publication.
    pub(crate) fn select_bulk_relaxations(
        &mut self,
        reference_rate: f64,
        higher_order_rate: f64,
        field: BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.admit_bulk_relaxations(&field, gate)?;
        self.set_central_moment_bulk(reference_rate, higher_order_rate).map_err(|_| {
            FreeSurfaceError3::Transport("invalid source bulk reference collision")
        })?;
        self.local_bulk_relaxation = Some(field);
        Ok(())
    }

    fn admit_bulk_relaxations(
        &self,
        field: &BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let mut wet = 0usize;
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let rate = field.get(&(key, lane)).copied().ok_or(
                    FreeSurfaceError3::Transport("missing source bulk relaxation at a wet cell"),
                )?;
                if !rate.is_finite() || rate <= 0.0 || rate >= 2.0 {
                    return Err(FreeSurfaceError3::Transport("invalid source bulk relaxation"));
                }
                wet += 1;
            }
        }
        if field.len() != wet {
            return Err(FreeSurfaceError3::Transport("source bulk field contains non-liquid cells"));
        }
        check_cancel(gate)?;
        Ok(())
    }

    pub(super) fn collision_with_bulk(
        &self,
        key: u64,
        lane: usize,
        model: CollisionModel3,
    ) -> Result<CollisionModel3, FreeSurfaceError3> {
        let Some(field) = &self.local_bulk_relaxation else { return Ok(model) };
        let bulk_rate = field.get(&(key, lane)).copied().ok_or(
            FreeSurfaceError3::Transport("missing source bulk relaxation at a wet cell"),
        )?;
        match model {
            CollisionModel3::CentralMomentBulk { shear_rate, higher_order_rate, .. } => {
                Ok(CollisionModel3::CentralMomentBulk { shear_rate, bulk_rate, higher_order_rate })
            }
            _ => Err(FreeSurfaceError3::Transport(
                "local bulk field lost its independent collision owner",
            )),
        }
    }

    /// Actual configured trace rate at an occupied collision cell. A newly wet
    /// cell needs source admission before a local volumetric reconstruction;
    /// only its raw shear/trace numerators can be read before that admission.
    pub(crate) fn collision_model_at(
        &self,
        coordinates: [usize; 3],
    ) -> Result<CollisionModel3, FreeSurfaceError3> {
        let Some(MassLocation3::Cell(key, lane)) = Self::transport_location(coordinates) else {
            return Err(FreeSurfaceError3::Transport("bulk relaxation requires a located cell"));
        };
        if self.tiles.get(&key).is_none_or(|tile| tile.cells[lane] == Cell3::Gas) {
            return Err(FreeSurfaceError3::Transport("bulk relaxation requires an active wet cell"));
        }
        self.collision_with_bulk(key, lane, self.collision_model())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::{E3, Q3, collide_cell3};
    use crate::d3q19::freesurface3::ContactModel3;

    #[test]
    fn source_bulk_field_controls_actual_local_collision_and_requires_complete_wet_coverage() {
        let gate = CancelGate::new_clock_free();
        let force = [0.0002, -0.0001, 0.0003];
        let mut flow = FreeSurface3::new(
            4, 4, 4, 0.9, force, 0.0, ContactModel3::Neutral, |_, _, _| true,
        ).unwrap();
        flow.set_central_moment_bulk(1.1, 1.2).unwrap();
        let mut bulk = BTreeMap::new();
        let mut shear = BTreeMap::new();
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..4 {
                    let Some(MassLocation3::Cell(key, lane)) = FreeSurface3::transport_location([x, y, z]) else { unreachable!() };
                    bulk.insert((key, lane), 0.7 + 0.2 * y as f64);
                    shear.insert((key, lane), 0.75 + 0.05 * x as f64);
                }
            }
        }
        flow.rebind_bulk_relaxations(bulk.clone(), &gate).unwrap();
        for (&(key, lane), &tau) in &shear {
            let c = FreeSurface3::transport_coordinates(MassLocation3::Cell(key, lane)).unwrap();
            flow.initialize_strain_at(c, [0.003, 0.003, 0.003, 0.002, 0.0, 0.0], tau).unwrap();
            let rates = flow.transport_strain_rates(MassLocation3::Cell(key, lane), tau).unwrap();
            assert!((rates[0] - 0.004).abs() < 2e-13);
            assert!((rates[1] - 0.009 * (2.0_f64 / 3.0).sqrt()).abs() < 2e-13);
        }
        let location = FreeSurface3::transport_location([2; 3]).unwrap();
        let MassLocation3::Cell(key, lane) = location else { unreachable!() };
        let mut incomplete = bulk.clone();
        incomplete.remove(&(key, lane));
        assert!(flow.rebind_bulk_relaxations(incomplete.clone(), &gate).is_err());
        assert_eq!(flow.local_bulk_relaxation.as_ref(), Some(&bulk));
        // Emulate a newly wet cell before its material source is resolved.
        // Raw shear can seed the source solve; actual trace/collision cannot.
        flow.local_bulk_relaxation = Some(incomplete);
        assert!(flow.transport_shear_rate(location, 1.0).unwrap() > 0.0);
        assert!(flow.transport_strain_rates(location, 1.0).is_err());
        let before = flow.grid.state_bits();
        assert!(flow.step_with_gate(&gate).is_err());
        assert_eq!(flow.grid.state_bits(), before);
        flow.rebind_bulk_relaxations(bulk.clone(), &gate).unwrap();
        let mut expected = [0.0; Q3];
        for (q, e) in E3.iter().enumerate() {
            let c = [2 - e.0, 2 - e.1, 2 - e.2];
            let (slot, lane) = flow.grid.resolve_source(i64::from(c[0]), i64::from(c[1]), i64::from(c[2])).unwrap();
            let key = flow.grid.active_keys()[slot];
            let model = CollisionModel3::CentralMomentBulk {
                shear_rate: 1.0 / shear[&(key, lane)],
                bulk_rate: bulk[&(key, lane)], higher_order_rate: 1.2,
            };
            expected[q] = collide_cell3(flow.grid.populations(slot, lane), model, force).unwrap()[q];
        }
        flow.step_inner(&gate, None, Some(&shear), None).unwrap();
        let (slot, lane) = flow.grid.resolve_source(2, 2, 2).unwrap();
        for (actual, expected) in flow.grid.populations(slot, lane).iter().zip(expected) {
            assert!((actual - expected).abs() < 3e-14);
        }
    }
}
