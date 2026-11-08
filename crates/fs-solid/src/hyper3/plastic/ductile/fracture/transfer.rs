//! Original consistent mass and P1 matter inventory through a nodal-star split.
//! Graded density stays nodal; its exact cubic density/shape mass is not replaced
//! by reporting cell means. Every physical integral is evaluated in unchanged
//! element/local-corner order, so the copied fields must give identical receipts.

use std::sync::Arc;

use super::{
    DifferentiableDuctileInitiationLaw, FiniteDuctileFractureError, FiniteDuctileFractureInventory,
    FiniteDuctileTetMotion, TetFractureTopology, TetLoadGeometry3, checkpoint, cross, energy_gate,
};
use crate::inertia::{P1ReferenceMass, p1_tet_nodal_mass};
use crate::plastic3::SmallStrainTetError;
use fs_exec::Cx;

pub(super) struct TransferredMotion {
    pub(super) mass: P1ReferenceMass,
    pub(super) nodal_density: Option<Arc<[f64]>>,
    pub(super) before: FiniteDuctileFractureInventory,
    pub(super) after: FiniteDuctileFractureInventory,
    pub(super) fragments: Vec<FiniteDuctileFractureInventory>,
    pub(super) kinetic_before_j: f64,
    pub(super) kinetic_after_j: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn motion<L>(
    old: &FiniteDuctileTetMotion<'_, '_, L>,
    topology: &TetFractureTopology<'_>,
    target: &TetLoadGeometry3<'_>,
    displacement: &[[f64; 3]],
    velocity: &[[f64; 3]],
    energy_tolerance_j: f64,
    cx: &Cx<'_>,
) -> Result<TransferredMotion, FiniteDuctileFractureError>
where
    L: DifferentiableDuctileInitiationLaw + Clone + core::fmt::Debug + PartialEq,
{
    use FiniteDuctileFractureError as E;
    checkpoint(cx)?;
    let map = topology.original_nodes();
    let nodal_density: Option<Arc<[f64]>> =
        old.reference_nodal_density_kg_m3.as_ref().map(|density| {
            map.iter()
                .map(|&node| density[node])
                .collect::<Vec<_>>()
                .into()
        });
    let geometry = || {
        target
            .tetrahedra()
            .iter()
            .zip(target.element_volumes_m3())
            .map(|(&nodes, &volume)| (nodes, volume))
    };
    let mass = if let Some(density) = &nodal_density {
        let (mass, means) = P1ReferenceMass::from_nodal_density(
            geometry(),
            density,
            3 * target.nodes_m().len(),
            cx,
        )
        .map_err(SmallStrainTetError::from)?;
        if means.as_slice() != old.reference_density_kg_m3.as_ref() {
            return Err(E::Invalid(
                "a split cannot change retained graded-density cell means",
            ));
        }
        mass
    } else {
        P1ReferenceMass::new(
            geometry(),
            &old.reference_density_kg_m3,
            3 * target.nodes_m().len(),
            cx,
        )
        .map_err(SmallStrainTetError::from)?
    };
    if mass.reference_mass_kg() != old.mass.reference_mass_kg() {
        return Err(E::Invalid(
            "a split cannot change original consistent reference mass",
        ));
    }
    let (before, _) = inventory(
        old.solid.reference_nodes_m(),
        old.solid.reference_tetrahedra(),
        target.element_volumes_m3(),
        &old.reference_density_kg_m3,
        old.reference_nodal_density_kg_m3.as_deref(),
        old.solid.displacement_m(),
        &old.velocity_mps,
        None,
        cx,
    )?;
    let (after, fragments) = inventory(
        target.nodes_m(),
        target.tetrahedra(),
        target.element_volumes_m3(),
        &old.reference_density_kg_m3,
        nodal_density.as_deref(),
        displacement,
        velocity,
        Some((topology.element_fragments(), topology.fragment_count())),
        cx,
    )?;
    if before != after {
        return Err(E::Invalid(
            "native nodal-star transfer changed matter or motion inventory",
        ));
    }
    let kinetic_before_j = kinetic(&old.mass, &old.velocity_mps, cx)?;
    let kinetic_after_j = kinetic(&mass, velocity, cx)?;
    for (quantity, defect) in [
        (
            "original consistent kinetic inventory",
            kinetic_before_j - before.kinetic_energy_j,
        ),
        (
            "split consistent kinetic inventory",
            kinetic_after_j - after.kinetic_energy_j,
        ),
        (
            "zero-duration kinetic continuity",
            kinetic_after_j - kinetic_before_j,
        ),
    ] {
        energy_gate(quantity, defect, energy_tolerance_j)?;
    }
    checkpoint(cx)?;
    Ok(TransferredMotion {
        mass,
        nodal_density,
        before,
        after,
        fragments,
        kinetic_before_j,
        kinetic_after_j,
    })
}

fn kinetic(
    mass: &P1ReferenceMass,
    velocity: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<f64, FiniteDuctileFractureError> {
    let momentum = mass
        .action(velocity.as_flattened(), cx)
        .map_err(SmallStrainTetError::from)?;
    let mut value = 0.0;
    for (dof, (&v, &p)) in velocity.as_flattened().iter().zip(&momentum).enumerate() {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        value = v.mul_add(p, value);
    }
    value *= 0.5;
    if !value.is_finite() || value < 0.0 {
        return Err(FiniteDuctileFractureError::Invalid(
            "finite nonnegative native kinetic storage required",
        ));
    }
    Ok(value)
}

#[allow(clippy::too_many_arguments)]
fn inventory(
    nodes: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    volumes: &[f64],
    means: &[f64],
    nodal_density: Option<&[f64]>,
    displacement: &[[f64; 3]],
    velocity: &[[f64; 3]],
    fragments: Option<(&[usize], usize)>,
    cx: &Cx<'_>,
) -> Result<
    (
        FiniteDuctileFractureInventory,
        Vec<FiniteDuctileFractureInventory>,
    ),
    FiniteDuctileFractureError,
> {
    let mut total = InventorySum::default();
    let mut pieces = vec![InventorySum::default(); fragments.map_or(0, |(_, count)| count)];
    for (element, &tet) in tetrahedra.iter().enumerate() {
        checkpoint(cx)?;
        let rho =
            nodal_density.map_or([means[element]; 4], |density| tet.map(|node| density[node]));
        let local = p1_tet_nodal_mass(volumes[element], rho).map_err(SmallStrainTetError::from)?;
        if local.mean_density_kg_m3 != means[element] {
            return Err(FiniteDuctileFractureError::Invalid(
                "retained native density inventory mismatch",
            ));
        }
        let x: [[f64; 3]; 4] =
            tet.map(|node| core::array::from_fn(|i| nodes[node][i] + displacement[node][i]));
        let v = tet.map(|node| velocity[node]);
        let mut one = InventorySum::default();
        one.add_value(0, local.mass_kg)?;
        for a in 0..4 {
            let row_mass: f64 = local.coefficients[a].iter().sum();
            let momentum: [f64; 3] = core::array::from_fn(|i| {
                (0..4).fold(0.0, |sum, b| local.coefficients[a][b].mul_add(v[b][i], sum))
            });
            let moment = cross(x[a], momentum);
            for i in 0..3 {
                one.add_value(1 + i, row_mass * x[a][i])?;
                one.add_value(4 + i, momentum[i])?;
                one.add_value(7 + i, moment[i])?;
            }
        }
        // Positive exact-polynomial form for integral rho*|v|^2/2:
        // each positive nodal-density term has Dirichlet weights (2,1,1,1).
        // This retains consistent cross terms without cancellation of signed
        // nodal work. Multiplying the mass weight first avoids velocity squares
        // overflowing when the physical kinetic energy remains representable.
        let density_scale = rho.into_iter().fold(0.0_f64, f64::max);
        let normalized = rho.map(|r| r / density_scale);
        let density_sum: f64 = normalized.iter().sum();
        let sum_v: [f64; 3] = core::array::from_fn(|i| v.iter().map(|q| q[i]).sum());
        for k in 0..4 {
            let weight = (local.mass_kg / 60.0) * (normalized[k] / density_sum);
            for i in 0..3 {
                let shifted = sum_v[i] + v[k][i];
                one.add_value(10, shifted * (weight * shifted))?;
                one.add_value(10, v[k][i] * (weight * v[k][i]))?;
                for va in &v {
                    one.add_value(10, va[i] * (weight * va[i]))?;
                }
            }
        }
        total.add(one)?;
        if let Some((ids, _)) = fragments {
            pieces[ids[element]].add(one)?;
        }
    }
    checkpoint(cx)?;
    Ok((
        total.finish(),
        pieces.into_iter().map(InventorySum::finish).collect(),
    ))
}

#[derive(Clone, Copy, Default)]
struct InventorySum {
    value: [f64; 11],
    correction: [f64; 11],
}
impl InventorySum {
    fn add_value(&mut self, i: usize, value: f64) -> Result<(), FiniteDuctileFractureError> {
        let adjusted = value - self.correction[i];
        let next = self.value[i] + adjusted;
        self.correction[i] = (next - self.value[i]) - adjusted;
        self.value[i] = next;
        if !value.is_finite() || !next.is_finite() || !self.correction[i].is_finite() {
            return Err(FiniteDuctileFractureError::Invalid(
                "finite represented matter inventory required",
            ));
        }
        Ok(())
    }
    fn add(&mut self, other: Self) -> Result<(), FiniteDuctileFractureError> {
        for i in 0..11 {
            self.add_value(i, other.value[i])?;
        }
        Ok(())
    }
    fn finish(self) -> FiniteDuctileFractureInventory {
        let v = self.value;
        FiniteDuctileFractureInventory {
            mass_kg: v[0],
            first_moment_kg_m: [v[1], v[2], v[3]],
            linear_momentum_kg_m_s: [v[4], v[5], v[6]],
            angular_momentum_kg_m2_s: [v[7], v[8], v[9]],
            kinetic_energy_j: v[10],
        }
    }
}
