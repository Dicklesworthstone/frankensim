//! Observe the original native mass operator; never rebuild or lump its inertia.

use super::*;
use fs_lbm::material::free_surface::transfer::LiquidKineticFrame;

pub(super) struct Native {
    pub moments: LiquidTransferMoments,
    pub frame: LiquidKineticFrame,
    pub maximum_inertial_force_n: f64,
}

pub(super) fn native(
    body: &FiniteSolidBody<'_>,
    positions: &[[f64; 3]],
    velocity: &[[f64; 3]],
    origin: [f64; 3],
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<Native, FiniteSolidLiquidReleaseError> {
    if origin.iter().any(|x| !x.is_finite())
        || body.free.len() != 3 * positions.len()
        || body.state.a.len() != body.free.len()
    {
        return Err(invalid(
            "unconstrained original mass coordinates and finite moment origin required",
        ));
    }
    let reference = *velocity
        .first()
        .ok_or_else(|| invalid("empty native velocity field"))?;
    let mut frame = LiquidKineticFrame {
        reference_velocity_m_s: reference,
        mass_kg: Dd::ZERO,
        relative_momentum_kg_m_s: [Dd::ZERO; 3],
        relative_kinetic_energy_j: [Dd::ZERO; 3],
    };
    let mut first = [Dd::ZERO; 3];
    let mut momentum = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    let mut maximum_inertial_force_n = 0.0_f64;
    for (node, position) in positions.iter().enumerate() {
        check(cx, gate)?;
        // The native 3-D P1 reference matrix contains the SAME scalar block
        // on all three displacement axes. Read its existing x rows, retaining
        // every off-diagonal mass entry and its original assembly roundoff.
        let (columns, values) = body.mass.row(3 * node);
        let mut row_mass = Dd::ZERO;
        let mut relative_momentum = [Dd::ZERO; 3];
        let mut inertial_force = [Dd::ZERO; 3];
        for (&column, &coefficient) in columns.iter().zip(values) {
            if column % 3 != 0 {
                if coefficient != 0.0 {
                    return Err(invalid(
                        "native reference inertia no longer has its original scalar Cartesian blocks",
                    ));
                }
                continue;
            }
            let other = column / 3;
            let mass = Dd::from_f64(coefficient);
            row_mass = row_mass + mass;
            for axis in 0..3 {
                relative_momentum[axis] = relative_momentum[axis]
                    + mass * (Dd::from_f64(velocity[other][axis]) - Dd::from_f64(reference[axis]));
                inertial_force[axis] =
                    inertial_force[axis] + mass * Dd::from_f64(body.state.a[3 * other + axis]);
            }
        }
        let relative_position = core::array::from_fn::<_, 3, _>(|axis| {
            Dd::from_f64(position[axis]) - Dd::from_f64(origin[axis])
        });
        let physical_momentum = core::array::from_fn::<_, 3, _>(|axis| {
            relative_momentum[axis] + row_mass * Dd::from_f64(reference[axis])
        });
        frame.mass_kg = frame.mass_kg + row_mass;
        for axis in 0..3 {
            first[axis] = first[axis] + row_mass * relative_position[axis];
            momentum[axis] = momentum[axis] + physical_momentum[axis];
            let j = (axis + 1) % 3;
            let k = (axis + 2) % 3;
            angular[axis] = angular[axis] + relative_position[j] * physical_momentum[k]
                - relative_position[k] * physical_momentum[j];
            frame.relative_momentum_kg_m_s[axis] =
                frame.relative_momentum_kg_m_s[axis] + relative_momentum[axis];
            frame.relative_kinetic_energy_j[axis] = frame.relative_kinetic_energy_j[axis]
                + Dd::from_f64(0.5)
                    * (Dd::from_f64(velocity[node][axis]) - Dd::from_f64(reference[axis]))
                    * relative_momentum[axis];
        }
        if !inertial_force.into_iter().all(finite_dd) {
            return Err(invalid(
                "native inertial-force observation is unrepresentable",
            ));
        }
        maximum_inertial_force_n = maximum_inertial_force_n.max(
            inertial_force
                .into_iter()
                .fold(0.0_f64, |n, f| n.hypot(f.to_f64())),
        );
    }
    let moments = LiquidTransferMoments {
        mass_kg: frame.mass_kg.to_f64(),
        first_mass_moment_kg_m: first.map(Dd::to_f64),
        momentum_kg_m_s: momentum.map(Dd::to_f64),
        angular_momentum_kg_m2_s: angular.map(Dd::to_f64),
        kinetic_energy_j: body.kinetic_energy_j()?,
    };
    if !finite_frame(frame) || !finite_moments(moments) || !maximum_inertial_force_n.is_finite() {
        return Err(invalid(
            "original consistent-mass inventory is unrepresentable",
        ));
    }
    Ok(Native {
        moments,
        frame,
        maximum_inertial_force_n,
    })
}

pub(super) fn kinetic_change(
    before: LiquidKineticFrame,
    after: LiquidKineticFrame,
) -> Result<Dd, FiniteSolidLiquidReleaseError> {
    if before.reference_velocity_m_s != after.reference_velocity_m_s
        || !finite_frame(before)
        || !finite_frame(after)
    {
        return Err(invalid(
            "native and tetrahedral kinetic inventories require one finite velocity reference",
        ));
    }
    let mass_change = after.mass_kg - before.mass_kg;
    let mut relative = Dd::ZERO;
    let mut momentum = Dd::ZERO;
    let mut mass = Dd::ZERO;
    for axis in 0..3 {
        let reference = Dd::from_f64(before.reference_velocity_m_s[axis]);
        relative = relative
            + (after.relative_kinetic_energy_j[axis] - before.relative_kinetic_energy_j[axis]);
        momentum = momentum
            + reference
                * (after.relative_momentum_kg_m_s[axis] - before.relative_momentum_kg_m_s[axis]);
        mass = mass + Dd::from_f64(0.5) * reference * reference * mass_change;
    }
    let change = relative + momentum + mass;
    finite_dd(change)
        .then_some(change)
        .ok_or_else(|| invalid("original-to-tetrahedral kinetic identity is unrepresentable"))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn audit(
    original: LiquidTransferMoments,
    original_h: f64,
    actual: LiquidTransferMoments,
    actual_h: f64,
    kinetic: Dd,
    controls: FiniteSolidLiquidReleaseControls,
    stage: FiniteSolidLiquidBalanceStage,
) -> Result<FiniteSolidLiquidDifference, FiniteSolidLiquidReleaseError> {
    let difference = FiniteSolidLiquidDifference {
        mass_kg: actual.mass_kg - original.mass_kg,
        momentum_kg_m_s: core::array::from_fn(|i| {
            actual.momentum_kg_m_s[i] - original.momentum_kg_m_s[i]
        }),
        angular_momentum_kg_m2_s: core::array::from_fn(|i| {
            actual.angular_momentum_kg_m2_s[i] - original.angular_momentum_kg_m2_s[i]
        }),
        centroid_m: core::array::from_fn(|i| {
            actual.first_mass_moment_kg_m[i] / actual.mass_kg
                - original.first_mass_moment_kg_m[i] / original.mass_kg
        }),
        kinetic_j: kinetic.to_f64(),
        enthalpy_j: actual_h - original_h,
    };
    let (centroid_limit, kinetic_limit) = match stage {
        FiniteSolidLiquidBalanceStage::TetrahedralSource => (
            controls.source_centroid_tolerance_m,
            controls.source_kinetic_tolerance_j,
        ),
        FiniteSolidLiquidBalanceStage::NativeLiquid => (
            controls.transfer.mechanical.maximum_centroid_shift_m,
            controls
                .transfer
                .mechanical
                .maximum_projection_energy_change_j,
        ),
    };
    let norm = |a: [f64; 3]| a[0].hypot(a[1]).hypot(a[2]);
    if !finite_moments(actual)
        || !finite_dd(kinetic)
        || !difference.enthalpy_j.is_finite()
        || !difference.mass_kg.is_finite()
        || !norm(difference.centroid_m).is_finite()
        || difference.mass_kg.abs() > controls.transfer.mechanical.mass_tolerance_kg
        || norm(difference.momentum_kg_m_s) > controls.transfer.mechanical.momentum_tolerance_kg_m_s
        || norm(difference.angular_momentum_kg_m2_s)
            > controls
                .transfer
                .mechanical
                .angular_momentum_tolerance_kg_m2_s
        || norm(difference.centroid_m) > centroid_limit
        || difference.enthalpy_j.abs() > controls.transfer.enthalpy_tolerance_j
        || Dd::from_f64(kinetic_limit).lt(kinetic.abs())
    {
        return Err(FiniteSolidLiquidReleaseError::Balance { stage, difference });
    }
    Ok(difference)
}

fn finite_frame(frame: LiquidKineticFrame) -> bool {
    frame.reference_velocity_m_s.iter().all(|v| v.is_finite())
        && Dd::ZERO.lt(frame.mass_kg)
        && core::iter::once(frame.mass_kg)
            .chain(frame.relative_momentum_kg_m_s)
            .chain(frame.relative_kinetic_energy_j)
            .all(finite_dd)
}

fn finite_moments(value: LiquidTransferMoments) -> bool {
    value.mass_kg.is_finite()
        && value.mass_kg > 0.0
        && value.kinetic_energy_j.is_finite()
        && value.kinetic_energy_j >= 0.0
        && value
            .first_mass_moment_kg_m
            .iter()
            .chain(&value.momentum_kg_m_s)
            .chain(&value.angular_momentum_kg_m2_s)
            .all(|x| x.is_finite())
}
