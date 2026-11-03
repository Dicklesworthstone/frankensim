//! Actual finite-pair exchanges and original native spatial-source reductions.

use super::{
    ArticulatedSustainedControls, ArticulatedSustainedError, ContactGeometryError, Dd, Frame,
    Impulses, LinkContactPoint, NativeError, Prepared, SpatialSum, Vec3, components, difference,
    dot_dd, finite_vector, point_parts, prepare_contact_frame, vector,
};

/// Signed exchange between two original finite material witnesses.
///
/// Forces are held in world axes and rebuilt at each original native stage.
/// Angular impulses use the actual midpoint arms about the step's declared
/// angular origin. Residuals retain relative geometry before rounding the
/// individual side receipts. These are internal exchanges, not source work.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedSustainedPairExchange {
    pub contact_id: u64,
    pub first: LinkContactPoint,
    pub second: LinkContactPoint,
    pub first_linear_impulse_world: Vec3,
    pub second_linear_impulse_world: Vec3,
    pub first_angular_impulse_about_origin_world: Vec3,
    pub second_angular_impulse_about_origin_world: Vec3,
    pub linear_impulse_residual_world: Vec3,
    pub angular_impulse_residual_world: Vec3,
    /// Actual applied world force paired with each original material displacement.
    pub first_work_j: f64,
    pub second_work_j: f64,
    /// Sum retained before rounding; closes normal work minus frictional loss.
    pub relative_work_j: f64,
}

pub(super) struct PairResponse {
    pub(super) exchange: ArticulatedSustainedPairExchange,
    geometry_error: Option<ContactGeometryError>,
    stage_moments: [Vec3; 3],
}

impl PairResponse {
    pub(super) fn admit<E>(
        &self,
        controls: ArticulatedSustainedControls,
    ) -> Result<(), ArticulatedSustainedError<E>> {
        if let Some(error) = &self.geometry_error {
            return Err(ArticulatedSustainedError::Geometry(error.clone()));
        }
        let e = &self.exchange;
        if [
            e.first_linear_impulse_world,
            e.second_linear_impulse_world,
            e.first_angular_impulse_about_origin_world,
            e.second_angular_impulse_about_origin_world,
            e.linear_impulse_residual_world,
            e.angular_impulse_residual_world,
        ]
        .into_iter()
        .chain(self.stage_moments)
        .any(|v| !finite_vector(v))
            || [e.first_work_j, e.second_work_j, e.relative_work_j]
                .into_iter()
                .any(|value| !value.is_finite())
        {
            return Err(ArticulatedSustainedError::Invalid(
                "nonfinite reciprocal contact exchange",
            ));
        }
        for &moment in &self.stage_moments {
            if moment.norm() > controls.angular_momentum_tolerance_nms {
                return Err(ArticulatedSustainedError::Reciprocity {
                    contact_id: e.contact_id,
                    angular_impulse_world: moment,
                });
            }
        }
        Ok(())
    }
}

/// Keep the original reducer's low components through source/origin reduction.
pub(super) fn add_spatial(sum: &mut Impulses, value: SpatialSum, origin: Vec3) {
    sum.add(
        vector(value.linear.map(|x| Dd::from_f64(x.hi))),
        vector(value.angular.map(|x| Dd::from_f64(x.hi))),
        origin,
    );
    sum.add(
        vector(value.linear.map(|x| Dd::from_f64(x.lo))),
        vector(value.angular.map(|x| Dd::from_f64(x.lo))),
        origin,
    );
}

#[allow(clippy::too_many_arguments)]
pub(super) fn record(
    prepared: &Prepared,
    second: LinkContactPoint,
    contact_index: usize,
    force_world: Vec3,
    middle: &[Frame],
    after: &[Frame],
    first_delta: [Dd; 3],
    second_delta: [Dd; 3],
    controls: ArticulatedSustainedControls,
) -> Result<PairResponse, NativeError> {
    let c = prepared.contact;
    let first_middle = point_parts(&middle[c.point.system_index], c.point)
        .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
    let second_middle = point_parts(&middle[second.system_index], second)
        .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
    let first_end = point_parts(&after[c.point.system_index], c.point)
        .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
    let second_end = point_parts(&after[second.system_index], second)
        .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
    let force = components(force_world).map(Dd::from_f64);
    let opposite = force.map(|x| -x);
    let mut first = SpatialSum::default();
    first.add_force(force, &first_middle);
    let mut second_sum = SpatialSum::default();
    second_sum.add_force(opposite, &second_middle);
    let mut first_receipt = Impulses::default();
    let mut second_receipt = Impulses::default();
    add_spatial(
        &mut first_receipt,
        first.scaled(controls.duration),
        controls.angular_origin_world,
    );
    add_spatial(
        &mut second_receipt,
        second_sum.scaled(controls.duration),
        controls.angular_origin_world,
    );

    // Relative geometry is reduced BEFORE crossing with the shared force. Two
    // huge rounded world moments must not hide a small artificial pair couple.
    let mut stage_moments = [Vec3::default(); 3];
    let mut paired_force = SpatialSum::default();
    paired_force.add_force(force, &[]);
    paired_force.add_force(opposite, &[]);
    let linear_residual = vector(paired_force.scaled(controls.duration).linear);
    for (slot, (a, b)) in stage_moments.iter_mut().zip([
        (prepared.initial_parts, prepared.initial_other_parts),
        (first_middle, second_middle),
        (first_end, second_end),
    ]) {
        let delta = difference(a, b);
        let parts = [
            vector(delta.map(|x| Dd::from_f64(x.hi))),
            vector(delta.map(|x| Dd::from_f64(x.lo))),
        ];
        let mut pair = SpatialSum::default();
        pair.add_force(force, &parts);
        pair.add_force(opposite, &[]);
        *slot = vector(pair.scaled(controls.duration).angular);
    }

    // Numerical trial points may leave the contact manifold. Retain geometry
    // failures for FINAL adoption; never clip the force or repair trial states.
    let mut geometry_error = None;
    if force_world != Vec3::default() {
        for stage in [middle, after] {
            if let Err(error) = prepare_contact_frame(
                c.point,
                c.other,
                c.normal_world,
                c.tangent_world,
                contact_index,
                stage,
            ) {
                geometry_error = Some(error);
                break;
            }
        }
    }
    let first_work = dot_dd(force_world, first_delta);
    let second_work = -dot_dd(force_world, second_delta);
    Ok(PairResponse {
        exchange: ArticulatedSustainedPairExchange {
            contact_id: c.id,
            first: c.point,
            second,
            first_linear_impulse_world: first_receipt.linear,
            second_linear_impulse_world: second_receipt.linear,
            first_angular_impulse_about_origin_world: first_receipt.angular,
            second_angular_impulse_about_origin_world: second_receipt.angular,
            linear_impulse_residual_world: linear_residual,
            angular_impulse_residual_world: stage_moments[1],
            first_work_j: first_work.to_f64(),
            second_work_j: second_work.to_f64(),
            relative_work_j: (first_work + second_work).to_f64(),
        },
        geometry_error,
        stage_moments,
    })
}
