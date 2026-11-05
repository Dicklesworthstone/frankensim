//! Additional cuts of an already split mesh, without collapsing its old fields.
//!
//! Every new nodal star must be a subset of one old star. Its displacement and
//! velocity therefore have ONE original receiving value even when opposite
//! crack faces already move independently. Looking up an unsplit node ID would
//! weld those histories back together. The element/corner correspondence below
//! constructs and checks the correct new-node to previous-node map instead.
//!
//! Existing seam identities, A/B roles and original vertex order remain the
//! prefix of the new seam list. P1 consumers must retain that reference order
//! for quadrature and local frames; sorting new DOF IDs can permute old history.
//! The public facade can verify retained P1 traces before any history is reused.
//! This operation admits kinematics and inventory, not a fracture criterion or
//! a traction/energy initialization for newly introduced surface laws.

use super::{
    FractureMotionTransfer, FractureTopologyBudget, FractureTopologyError, InteriorSeam, SplitFace,
    TetFractureTopology, finite, inventory, poll,
};

/// Immutable refinement and its unambiguous motion-transfer map. The previous
/// topology stays borrowed, so the map cannot silently target another mesh.
#[derive(Debug)]
pub struct TetFractureRefinement<'p, 'g> {
    previous: &'p TetFractureTopology<'g>,
    topology: TetFractureTopology<'g>,
    previous_nodes: Vec<usize>,
}

impl<'g> TetFractureTopology<'g> {
    /// Add original-reference faces to the existing cut set. Existing cuts
    /// cannot be removed, merged, reordered, reassigned or silently recreated.
    /// The supplied budget applies to the COMPLETE resulting topology.
    ///
    /// The original admitted geometry and node-star splitter remain the only
    /// geometric owner. A refused candidate leaves both old cuts and motion
    /// reusable. No mass, force, damage or time is changed by preparation.
    pub fn refine_controlled<'p, C: FnMut() -> bool>(
        &'p self,
        additional_seams: &[InteriorSeam],
        budget: FractureTopologyBudget,
        cancelled: &mut C,
    ) -> Result<TetFractureRefinement<'p, 'g>, FractureTopologyError> {
        poll(cancelled)?;
        let count = self
            .seams
            .len()
            .checked_add(additional_seams.len())
            .ok_or(FractureTopologyError::Budget("refined seam count overflow"))?;
        if count > budget.maximum_seams {
            return Err(FractureTopologyError::Budget("complete refined seam count"));
        }
        let mut seams = Vec::with_capacity(count);
        for seam in &self.seams {
            poll(cancelled)?;
            seams.push(InteriorSeam {
                id: seam.id,
                nodes: seam.original_nodes,
                side_a_element: seam.elements[0],
            });
        }
        for &seam in additional_seams {
            poll(cancelled)?;
            seams.push(seam);
        }
        let topology = Self::prepare(
            self.original_nodes_m,
            self.original_tetrahedra,
            self.volumes_m3,
            &seams,
            budget,
            cancelled,
        )?;
        let mut previous_nodes = vec![usize::MAX; topology.nodes_m.len()];
        for (old_tet, new_tet) in self.tetrahedra.iter().zip(&topology.tetrahedra) {
            poll(cancelled)?;
            for local in 0..4 {
                let old = old_tet[local];
                let new = new_tet[local];
                let mapped = &mut previous_nodes[new];
                if (*mapped != usize::MAX && *mapped != old)
                    || self.original_nodes[old] != topology.original_nodes[new]
                    || self.nodes_m[old] != topology.nodes_m[new]
                {
                    return Err(FractureTopologyError::Invalid(
                        "refinement would weld distinct previous nodal stars",
                    ));
                }
                *mapped = old;
            }
        }
        for chunk in previous_nodes.chunks(256) {
            poll(cancelled)?;
            if chunk.iter().any(|&node| node == usize::MAX) {
                return Err(FractureTopologyError::Invalid("unmapped refined node"));
            }
        }
        poll(cancelled)?;
        Ok(TetFractureRefinement {
            previous: self,
            topology,
            previous_nodes,
        })
    }
}

impl<'p, 'g> TetFractureRefinement<'p, 'g> {
    /// Exact previous split mesh to which input nodal fields must belong.
    #[must_use]
    pub const fn previous(&self) -> &'p TetFractureTopology<'g> {
        self.previous
    }

    /// Prepared mesh containing every old seam followed by the supplied cuts.
    #[must_use]
    pub const fn topology(&self) -> &TetFractureTopology<'g> {
        &self.topology
    }

    /// NEW node -> PREVIOUS split-mesh node, not an unsplit-reference node ID.
    /// This map copies intensive fields. External forces, prescribed supports
    /// and matter sources require their own physical boundary/load assignment.
    #[must_use]
    pub fn previous_nodes(&self) -> &[usize] {
        &self.previous_nodes
    }

    /// Original seams remain this unchanged prefix in reference identity/order.
    /// Retaining IDs alone does not establish identical quadrature/local frames.
    #[must_use]
    pub fn retained_seam_count(&self) -> usize {
        self.previous.seams.len()
    }

    /// Newly introduced surfaces; this method initializes no cohesive history.
    #[must_use]
    pub fn additional_seams(&self) -> &[SplitFace] {
        &self.topology.seams[self.retained_seam_count()..]
    }

    /// Adopt the prepared topology after all required state transfers succeed.
    /// Material/source and element-history arrays retain their original order.
    #[must_use]
    pub fn into_topology(self) -> TetFractureTopology<'g> {
        self.topology
    }

    /// Transfer fields on the PREVIOUS split mesh, preserving independent
    /// motion on old crack faces. Consistent P1 inventories are independently
    /// reconstructed before/after using original element order and densities.
    /// The exact same represented integrals must agree before publication.
    /// No element is deleted, no inertia is lumped and no velocity is averaged.
    pub fn transfer_motion_controlled<C: FnMut() -> bool>(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        density_kg_m3: &[f64],
        cancelled: &mut C,
    ) -> Result<FractureMotionTransfer, FractureTopologyError> {
        poll(cancelled)?;
        if displacement_m.len() != 3 * self.previous.nodes_m.len()
            || velocity_m_s.len() != displacement_m.len()
            || density_kg_m3.len() != self.previous.tetrahedra.len()
        {
            return Err(FractureTopologyError::Invalid(
                "motion must belong to the previous split topology",
            ));
        }
        for values in [displacement_m, velocity_m_s, density_kg_m3] {
            for chunk in values.chunks(256) {
                poll(cancelled)?;
                for &value in chunk {
                    finite(value)?;
                }
            }
        }
        if density_kg_m3.iter().any(|&density| density <= 0.0) {
            return Err(FractureTopologyError::Invalid(
                "positive original reference mass density required",
            ));
        }
        let before = inventory(
            &self.previous.nodes_m,
            &self.previous.tetrahedra,
            self.previous.volumes_m3,
            density_kg_m3,
            displacement_m,
            velocity_m_s,
            None,
            cancelled,
        )?
        .0;
        let mut displacement = Vec::with_capacity(3 * self.previous_nodes.len());
        let mut velocity = Vec::with_capacity(3 * self.previous_nodes.len());
        for (index, &node) in self.previous_nodes.iter().enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            displacement.extend_from_slice(&displacement_m[3 * node..3 * node + 3]);
            velocity.extend_from_slice(&velocity_m_s[3 * node..3 * node + 3]);
        }
        let (after, fragments, row_sum_mass_kg) = inventory(
            &self.topology.nodes_m,
            &self.topology.tetrahedra,
            self.topology.volumes_m3,
            density_kg_m3,
            &displacement,
            &velocity,
            Some((
                &self.topology.element_fragments,
                self.topology.fragment_count,
            )),
            cancelled,
        )?;
        if before != after {
            return Err(FractureTopologyError::Conservation);
        }
        poll(cancelled)?;
        Ok(FractureMotionTransfer {
            displacement_m: displacement,
            velocity_m_s: velocity,
            row_sum_mass_kg,
            before,
            after,
            fragments,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::{P1PhysicalInventory, cross};
    use super::*;

    const NODES: [[f64; 3]; 5] = [
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [-0.5, 0.9, 0.0],
        [-0.5, -0.9, 0.0],
    ];
    const TETS: [[usize; 4]; 3] = [[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 2]];
    const VOLUMES: [f64; 3] = [0.3; 3];
    const FIRST: InteriorSeam = InteriorSeam {
        id: 19,
        nodes: [3, 0, 1],
        side_a_element: 1,
    };
    const MORE: [InteriorSeam; 2] = [
        InteriorSeam {
            id: 27,
            nodes: [0, 1, 4],
            side_a_element: 2,
        },
        InteriorSeam {
            id: 38,
            nodes: [2, 1, 0],
            side_a_element: 0,
        },
    ];
    fn budget() -> FractureTopologyBudget {
        FractureTopologyBudget {
            maximum_output_nodes: 32,
            maximum_face_entries: 12,
            maximum_seams: 3,
        }
    }
    fn first() -> TetFractureTopology<'static> {
        TetFractureTopology::prepare(&NODES, &TETS, &VOLUMES, &[FIRST], budget(), &mut || false)
            .unwrap()
    }
    fn fields(topology: &TetFractureTopology<'_>) -> (Vec<f64>, Vec<f64>) {
        // Old split vertices intentionally differ despite identical reference
        // coordinates: looking up original_nodes() would erase this motion.
        let count = 3 * topology.nodes_m.len();
        (
            (0..count)
                .map(|i| 0.001 * ((i * 11 % 19) as f64 - 9.0))
                .collect(),
            (0..count)
                .map(|i| 0.2 * ((i * 7 % 23) as f64 - 11.0))
                .collect(),
        )
    }
    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 4e-13 * a.abs().max(b.abs()).max(1.0),
            "{a:e} != {b:e}"
        );
    }
    fn compare(a: P1PhysicalInventory, b: P1PhysicalInventory) {
        near(a.mass_kg, b.mass_kg);
        near(a.kinetic_energy_j, b.kinetic_energy_j);
        for c in 0..3 {
            near(a.first_moment_kg_m[c], b.first_moment_kg_m[c]);
            near(a.linear_momentum_kg_m_s[c], b.linear_momentum_kg_m_s[c]);
            near(a.angular_momentum_kg_m2_s[c], b.angular_momentum_kg_m2_s[c]);
        }
    }
    // Independent positive degree-two tetrahedral quadrature, including the
    // cross terms of x cross v and v dot v. No production mass formula is used.
    fn quadrature(
        topology: &TetFractureTopology<'_>,
        u: &[f64],
        v: &[f64],
        rho: &[f64],
    ) -> P1PhysicalInventory {
        let a = 0.5854101966249685;
        let b = 0.1381966011250105;
        let mut result = P1PhysicalInventory::default();
        for (element, tet) in topology.tetrahedra.iter().enumerate() {
            for q in 0..4 {
                let mut x = [0.0; 3];
                let mut velocity = [0.0; 3];
                for i in 0..4 {
                    let shape = if i == q { a } else { b };
                    for c in 0..3 {
                        x[c] += shape * (topology.nodes_m[tet[i]][c] + u[3 * tet[i] + c]);
                        velocity[c] += shape * v[3 * tet[i] + c];
                    }
                }
                let mass = rho[element] * VOLUMES[element] / 4.0;
                result.mass_kg += mass;
                let moment = cross(x, velocity);
                for c in 0..3 {
                    result.first_moment_kg_m[c] += mass * x[c];
                    result.linear_momentum_kg_m_s[c] += mass * velocity[c];
                    result.angular_momentum_kg_m2_s[c] += mass * moment[c];
                    result.kinetic_energy_j += 0.5 * mass * velocity[c] * velocity[c];
                }
            }
        }
        result
    }

    #[test]
    fn growing_cuts_preserve_already_discontinuous_motion_and_consistent_integrals() {
        let previous = first();
        let (u, v) = fields(&previous);
        let rho = [6.0, 12.0, 9.0];
        let original_tets = previous.tetrahedra.clone();
        let refinement = previous
            .refine_controlled(&MORE, budget(), &mut || false)
            .unwrap();
        assert_eq!(refinement.retained_seam_count(), 1);
        assert_eq!(
            refinement
                .additional_seams()
                .iter()
                .map(|s| s.id)
                .collect::<Vec<_>>(),
            vec![27, 38]
        );
        assert_eq!(refinement.topology().element_fragments(), &[0, 1, 2]);
        assert_eq!(
            refinement.topology().seams()[0].elements,
            previous.seams()[0].elements
        );
        let transfer = refinement
            .transfer_motion_controlled(&u, &v, &rho, &mut || false)
            .unwrap();
        assert_eq!(transfer.before, transfer.after);
        compare(transfer.before, quadrature(&previous, &u, &v, &rho));
        compare(
            transfer.after,
            quadrature(
                refinement.topology(),
                &transfer.displacement_m,
                &transfer.velocity_m_s,
                &rho,
            ),
        );
        for (new, &old) in refinement.previous_nodes().iter().enumerate() {
            for c in 0..3 {
                assert_eq!(
                    transfer.displacement_m[3 * new + c].to_bits(),
                    u[3 * old + c].to_bits()
                );
                assert_eq!(
                    transfer.velocity_m_s[3 * new + c].to_bits(),
                    v[3 * old + c].to_bits()
                );
            }
        }
        let new_face = refinement.topology().seams()[0];
        let old_face = previous.seams()[0];
        for i in 0..3 {
            assert_eq!(
                refinement.previous_nodes()[new_face.side_a[i]],
                old_face.side_a[i]
            );
            assert_eq!(
                refinement.previous_nodes()[new_face.side_b[i]],
                old_face.side_b[i]
            );
            for c in 0..3 {
                near(
                    transfer.displacement_m[3 * new_face.side_b[i] + c]
                        - transfer.displacement_m[3 * new_face.side_a[i] + c],
                    u[3 * old_face.side_b[i] + c] - u[3 * old_face.side_a[i] + c],
                );
            }
        }
        struct History {
            loss_j: f64,
        }
        let history = [
            History { loss_j: 2.0 },
            History { loss_j: 7.0 },
            History { loss_j: 11.0 },
        ];
        let retained = refinement
            .topology()
            .retain_element_values(&history)
            .unwrap();
        assert!(core::ptr::eq(retained.as_ptr(), history.as_ptr()));
        assert_eq!(retained.iter().map(|h| h.loss_j).sum::<f64>(), 20.0);
        assert_eq!(previous.tetrahedra, original_tets);
        near(
            transfer.fragments.iter().map(|x| x.mass_kg).sum(),
            transfer.after.mass_kg,
        );
    }

    #[test]
    fn sequential_refinements_compose_the_same_motion_map_as_one_combined_cut() {
        let old = first();
        let (u, v) = fields(&old);
        let rho = [6.0, 12.0, 9.0];
        let first = old
            .refine_controlled(&MORE[..1], budget(), &mut || false)
            .unwrap();
        let moved = first
            .transfer_motion_controlled(&u, &v, &rho, &mut || false)
            .unwrap();
        let second = first
            .topology()
            .refine_controlled(&MORE[1..], budget(), &mut || false)
            .unwrap();
        let sequential = second
            .transfer_motion_controlled(
                &moved.displacement_m,
                &moved.velocity_m_s,
                &rho,
                &mut || false,
            )
            .unwrap();
        let direct = old
            .refine_controlled(&MORE, budget(), &mut || false)
            .unwrap();
        let combined = direct
            .transfer_motion_controlled(&u, &v, &rho, &mut || false)
            .unwrap();
        assert_eq!(
            second.topology().tetrahedra(),
            direct.topology().tetrahedra()
        );
        assert_eq!(sequential.displacement_m, combined.displacement_m);
        assert_eq!(sequential.velocity_m_s, combined.velocity_m_s);
        assert_eq!(sequential.after, combined.after);
        for (new, &intermediate) in second.previous_nodes().iter().enumerate() {
            assert_eq!(
                first.previous_nodes()[intermediate],
                direct.previous_nodes()[new]
            );
        }
        // Identity refinement keeps independently moving old faces intact too.
        let identity = old.refine_controlled(&[], budget(), &mut || false).unwrap();
        let retained = identity
            .transfer_motion_controlled(&u, &v, &rho, &mut || false)
            .unwrap();
        assert_eq!(retained.displacement_m, u);
        assert_eq!(retained.velocity_m_s, v);
    }

    #[test]
    fn old_cuts_cannot_be_redefined_and_complete_output_budgets_apply() {
        let old = first();
        for repeated in [
            FIRST,
            InteriorSeam { id: 80, ..FIRST },
            InteriorSeam {
                id: FIRST.id,
                ..MORE[0]
            },
        ] {
            assert!(
                old.refine_controlled(&[repeated], budget(), &mut || false)
                    .is_err()
            );
        }
        let full = old
            .refine_controlled(&MORE, budget(), &mut || false)
            .unwrap();
        let mut cap = budget();
        cap.maximum_output_nodes = full.topology().nodes_m.len();
        assert!(old.refine_controlled(&MORE, cap, &mut || false).is_ok());
        cap.maximum_output_nodes -= 1;
        assert!(matches!(
            old.refine_controlled(&MORE, cap, &mut || false),
            Err(FractureTopologyError::Budget(_))
        ));
        cap = budget();
        cap.maximum_seams = 2;
        assert!(matches!(
            old.refine_controlled(&MORE, cap, &mut || false),
            Err(FractureTopologyError::Budget(_))
        ));
        let (u, v) = fields(&old);
        assert!(
            full.transfer_motion_controlled(&u[..15], &v[..15], &[6.0; 3], &mut || false)
                .is_err()
        );
        assert!(
            full.transfer_motion_controlled(&u, &v, &[6.0; 2], &mut || false)
                .is_err()
        );
        assert!(
            full.transfer_motion_controlled(&u, &v, &[6.0, 0.0, 6.0], &mut || false)
                .is_err()
        );
        let mut invalid = v.clone();
        invalid[0] = f64::NAN;
        assert!(
            full.transfer_motion_controlled(&u, &invalid, &[6.0; 3], &mut || false)
                .is_err()
        );
    }

    #[test]
    fn cancellation_during_refinement_and_motion_leaves_accepted_fields_reusable() {
        let old = first();
        let (u, v) = fields(&old);
        let saved = (u.clone(), v.clone());
        let mut calls = 0;
        let prepared = old
            .refine_controlled(&MORE, budget(), &mut || {
                calls += 1;
                false
            })
            .unwrap();
        for stop in [1, calls / 2, calls] {
            let mut observed = 0;
            assert!(matches!(
                old.refine_controlled(&MORE, budget(), &mut || {
                    observed += 1;
                    observed == stop
                }),
                Err(FractureTopologyError::Cancelled)
            ));
        }
        let mut calls = 0;
        let result = prepared
            .transfer_motion_controlled(&u, &v, &[6.0; 3], &mut || {
                calls += 1;
                false
            })
            .unwrap();
        for stop in [1, calls / 2, calls] {
            let mut observed = 0;
            assert!(matches!(
                prepared.transfer_motion_controlled(&u, &v, &[6.0; 3], &mut || {
                    observed += 1;
                    observed == stop
                }),
                Err(FractureTopologyError::Cancelled)
            ));
        }
        assert_eq!((u.clone(), v.clone()), saved);
        assert_eq!(
            result,
            prepared
                .transfer_motion_controlled(&u, &v, &[6.0; 3], &mut || false)
                .unwrap()
        );
    }
}
