//! Shared bilateral geometry and constrained contact mobility.
//!
//! Geometry is independent of the mass provider. The same normalized Schur
//! construction serves native articulated and complete dry/fluid responses;
//! constitutive contact laws and physical endpoint audits remain their owners.

use super::ConstrainedImpactConfig;
use super::bilateral::{
    ImpactBilateralConstraint, ImpactBilateralRow, ImpactBilateralRowKind,
    MAX_IMPACT_BILATERAL_ROWS,
};
use crate::articulated::contact_geometry::{
    LinkContactPoint, PointDirection, prepare_point, unit_normal,
};
use crate::articulated::{Kinematics, solve_floating_base_system};
use crate::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, Prepared, dynamics, finite, poll,
};
use fs_ga::{Vec3, Wrench};
use fs_la::factor::cholesky;
use std::borrow::Borrow;

/// One objective bilateral row expressed as two original link covectors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BilateralImpactGeometry {
    /// Positive-side point force or pure angular couple.
    pub first: PointDirection,
    /// Negative-side covector, including the point-closure geometric couple.
    pub second: PointDirection,
    /// Row axis in world coordinates; angular for weld angular rows.
    pub axis_world: Vec3,
    /// Original positive-side material witness.
    pub world_point: Vec3,
    /// Original constraint index and physical row interpretation.
    pub row: ImpactBilateralRow,
    /// Supplied distance/point closure error; angular weld rows are zero.
    pub position_residual_m: f64,
    /// Objective supplied rate: m/s for point/range, rad/s for angular rows.
    pub velocity_before: f64,
}

impl BilateralImpactGeometry {
    pub(super) fn prepared(self) -> Prepared {
        Prepared {
            first: self.first,
            second: Some(self.second),
            axis_world: self.axis_world,
            world_point: self.world_point,
            plane_speed: 0.,
            minimum_speed: None,
        }
    }
}

/// Admit the original native constraint controls before doing physical work.
/// The native velocity budget applies in the units of each declared row.
pub fn validate_constraint_config(
    cfg: ConstrainedImpactConfig,
) -> Result<(), ArticulatedContactError> {
    if [
        cfg.position_tolerance_m,
        cfg.velocity_tolerance_mps,
        cfg.minimum_scaled_pivot,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.)
        || cfg.minimum_scaled_pivot < 4096. * f64::EPSILON
        || cfg.minimum_scaled_pivot >= 1.
    {
        return Err(ArticulatedContactError::Invalid(
            "invalid constrained-impact controls",
        ));
    }
    Ok(())
}

/// Prepare distance, spherical and weld covectors from the same native maps.
///
/// Initial geometry and objective row rates must already satisfy the supplied
/// tolerances. No position or velocity is repaired. Point closures retain the
/// A-side geometric couple of `R_A^T(x_B-x_A)`, even at an admitted nonzero
/// gap, so their internal impulse has zero net world spatial moment.
/// Borrowed maps avoid cloning whole articulated trees.
pub fn prepare_bilateral_rows(
    constraints: &[ImpactBilateralConstraint],
    cfg: ConstrainedImpactConfig,
    kinematics: &[impl Borrow<Kinematics>],
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<BilateralImpactGeometry>, ArticulatedContactError> {
    poll(&mut cancelled)?;
    validate_constraint_config(cfg)?;
    if constraints.len() > MAX_IMPACT_BILATERAL_ROWS {
        return Err(ArticulatedContactError::Invalid(
            "bilateral-impact work bound exceeded",
        ));
    }
    let count: usize = constraints
        .iter()
        .map(|c| match c {
            ImpactBilateralConstraint::Distance(_) => 1,
            ImpactBilateralConstraint::PointCoincidence(_) => 3,
            ImpactBilateralConstraint::Weld(_) => 6,
        })
        .sum();
    if count > MAX_IMPACT_BILATERAL_ROWS {
        return Err(ArticulatedContactError::Invalid(
            "bilateral-impact work bound exceeded",
        ));
    }
    let mut rows = Vec::with_capacity(count);
    for (constraint, c) in constraints.iter().copied().enumerate() {
        poll(&mut cancelled)?;
        match c {
            ImpactBilateralConstraint::Distance(c) => {
                if !c.length_m.is_finite() || c.length_m <= 0. || same_link(c.a, c.b) {
                    return Err(ArticulatedContactError::Invalid(
                        "invalid impact distance endpoints/length",
                    ));
                }
                let (_, a) = prepare_point(c.a, Vec3::new(1., 0., 0.), kinematics)?;
                let (_, b) = prepare_point(c.b, Vec3::new(1., 0., 0.), kinematics)?;
                let delta = b - a;
                let gap = finite(finite(delta.norm())? - c.length_m)?;
                if gap.abs() > cfg.position_tolerance_m {
                    return Err(ArticulatedContactError::Invalid(
                        "initial rigid distance is inconsistent",
                    ));
                }
                let axis_world = unit_normal(delta)?;
                let (first, world_point) = prepare_point(c.b, axis_world, kinematics)?;
                let (second, _) = prepare_point(c.a, axis_world.scale(-1.), kinematics)?;
                let rate = row_rate(first, second, kinematics)?;
                if rate.abs() > cfg.velocity_tolerance_mps {
                    return Err(ArticulatedContactError::Invalid(
                        "initial rigid distance velocity is not tangent",
                    ));
                }
                rows.push(BilateralImpactGeometry {
                    first,
                    second,
                    axis_world,
                    world_point,
                    row: ImpactBilateralRow {
                        constraint,
                        kind: ImpactBilateralRowKind::Distance,
                    },
                    position_residual_m: gap,
                    velocity_before: rate,
                });
            }
            ImpactBilateralConstraint::PointCoincidence(c) => {
                append_point(c.a, c.b, constraint, false, cfg, kinematics, &mut rows)?;
            }
            ImpactBilateralConstraint::Weld(c) => {
                append_point(c.a, c.b, constraint, true, cfg, kinematics, &mut rows)?;
                append_angular(c.a, c.b, constraint, cfg, kinematics, &mut rows)?;
            }
        }
    }
    poll(&mut cancelled)?;
    Ok(rows)
}

fn row_rate(
    first: PointDirection,
    second: PointDirection,
    kinematics: &[impl Borrow<Kinematics>],
) -> Result<f64, ArticulatedContactError> {
    finite(
        [first, second]
            .iter()
            .map(|term| {
                term.direction.pairing(
                    kinematics[term.point.system_index].borrow().body_twist[term.point.link_index],
                )
            })
            .sum(),
    )
}

const BASIS: [Vec3; 3] = [
    Vec3 {
        x: 1.,
        y: 0.,
        z: 0.,
    },
    Vec3 {
        x: 0.,
        y: 1.,
        z: 0.,
    },
    Vec3 {
        x: 0.,
        y: 0.,
        z: 1.,
    },
];

#[allow(clippy::too_many_arguments)]
fn append_point(
    a: LinkContactPoint,
    b: LinkContactPoint,
    constraint: usize,
    weld: bool,
    cfg: ConstrainedImpactConfig,
    kinematics: &[impl Borrow<Kinematics>],
    rows: &mut Vec<BilateralImpactGeometry>,
) -> Result<(), ArticulatedContactError> {
    if same_link(a, b) {
        return Err(ArticulatedContactError::Invalid(
            "point constraint requires distinct rigid links",
        ));
    }
    let (_, a_world) = prepare_point(a, BASIS[0], kinematics)?;
    let (_, b_world) = prepare_point(b, BASIS[0], kinematics)?;
    let pose_a = kinematics[a.system_index].borrow().world_from_link[a.link_index];
    let delta_world = b_world - a_world;
    let delta_a = pose_a.rotation().inverse().rotate(delta_world)?;
    if delta_a.norm() > cfg.position_tolerance_m {
        return Err(ArticulatedContactError::Invalid(
            "initial point closure is inconsistent",
        ));
    }
    for (component, axis_a) in BASIS.into_iter().enumerate() {
        let axis_world = pose_a.rotation().rotate(axis_a)?;
        let (first, world_point) = prepare_point(b, axis_world, kinematics)?;
        let (mut second, _) = prepare_point(a, axis_world.scale(-1.), kinematics)?;
        let couple_world = delta_world.cross(axis_world).scale(-1.);
        second.direction.torque =
            second.direction.torque + pose_a.rotation().inverse().rotate(couple_world)?;
        let rate = row_rate(first, second, kinematics)?;
        if rate.abs() > cfg.velocity_tolerance_mps {
            return Err(ArticulatedContactError::Invalid(
                "initial point-closure velocity is not tangent",
            ));
        }
        rows.push(BilateralImpactGeometry {
            first,
            second,
            axis_world,
            world_point,
            row: ImpactBilateralRow {
                constraint,
                kind: if weld {
                    ImpactBilateralRowKind::WeldPoint { component }
                } else {
                    ImpactBilateralRowKind::Point { component }
                },
            },
            position_residual_m: [delta_a.x, delta_a.y, delta_a.z][component],
            velocity_before: rate,
        });
    }
    Ok(())
}

fn append_angular(
    a: LinkContactPoint,
    b: LinkContactPoint,
    constraint: usize,
    cfg: ConstrainedImpactConfig,
    kinematics: &[impl Borrow<Kinematics>],
    rows: &mut Vec<BilateralImpactGeometry>,
) -> Result<(), ArticulatedContactError> {
    // append_point has already admitted both material endpoints and maps.
    let pose_a = kinematics[a.system_index].borrow().world_from_link[a.link_index];
    let pose_b = kinematics[b.system_index].borrow().world_from_link[b.link_index];
    for (component, axis_a) in BASIS.into_iter().enumerate() {
        let axis_world = pose_a.rotation().rotate(axis_a)?;
        let first = PointDirection {
            point: b,
            direction: Wrench::new(
                pose_b.rotation().inverse().rotate(axis_world)?,
                Vec3::default(),
            ),
        };
        let second = PointDirection {
            point: a,
            direction: Wrench::new(
                pose_a.rotation().inverse().rotate(axis_world.scale(-1.))?,
                Vec3::default(),
            ),
        };
        let rate = row_rate(first, second, kinematics)?;
        if rate.abs() > cfg.velocity_tolerance_mps {
            return Err(ArticulatedContactError::Invalid(
                "initial weld angular velocity is not tangent",
            ));
        }
        rows.push(BilateralImpactGeometry {
            first,
            second,
            axis_world,
            world_point: pose_b.transform_point(b.point_body)?,
            row: ImpactBilateralRow {
                constraint,
                kind: ImpactBilateralRowKind::WeldAngular { component },
            },
            position_residual_m: 0.,
            velocity_before: rate,
        });
    }
    Ok(())
}

fn same_link(a: LinkContactPoint, b: LinkContactPoint) -> bool {
    a.system_index == b.system_index && a.link_index == b.link_index
}

/// The existing independent-endpoint dry mobility used by native adapters.
pub(super) fn native_reference_mobility(
    systems: &[ArticulatedContactSystem<'_>],
    rows: &[Prepared],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, ArticulatedContactError> {
    let mut reference = Vec::with_capacity(rows.len());
    for row in rows {
        poll(cancelled)?;
        let mut free = 0.;
        for term in row.terms() {
            let system = &systems[term.point.system_index];
            let inertia = system.model.links()[term.point.link_index].inertia();
            let matrix = dynamics(term.point.system_index, inertia.matrix())?;
            let response = dynamics(
                term.point.system_index,
                solve_floating_base_system(&matrix, term.direction),
            )?;
            free += term.direction.pairing(response);
        }
        reference.push(free);
    }
    Ok(reference)
}

/// One normalized bilateral factorization and its contact-reaction recovery.
///
/// The provider supplies `[D B; B^T A]` in contact-first order. The physical
/// contact mobility is `D - B A^-1 B^T`; recovered reactions are
/// `-A^-1 B^T p`. No penalty stiffness, inverse generalized mass, regularizing
/// diagonal or velocity correction is introduced.
#[derive(Debug, Clone)]
pub struct ConstrainedMobility {
    /// Condensed contact matrix, row-major in the original contact row order.
    pub contact_mobility: Vec<f64>,
    /// Smallest squared Cholesky pivot of the normalized bilateral block.
    pub minimum_scaled_pivot: f64,
    contact_rows: usize,
    bilateral_rows: usize,
    recovery: Vec<f64>,
}

impl ConstrainedMobility {
    /// Factor up to 32 independent bilateral rows and condense up to 63 contact
    /// rows. `reference_mobility` gives a positive, uncancelled physical response
    /// scale per bilateral row. It must use the provider's actual inertia: a
    /// dry-only scale must not reject a legitimately large attached-fluid mass.
    pub fn new(
        gram: &[f64],
        contact_rows: usize,
        reference_mobility: &[f64],
        minimum_scaled_pivot: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ArticulatedContactError> {
        poll(&mut cancelled)?;
        let n = contact_rows;
        let k = reference_mobility.len();
        if n > 63
            || k > MAX_IMPACT_BILATERAL_ROWS
            || !minimum_scaled_pivot.is_finite()
            || !(4096. * f64::EPSILON..1.).contains(&minimum_scaled_pivot)
        {
            return Err(ArticulatedContactError::Invalid(
                "invalid bilateral mobility bounds",
            ));
        }
        let dimension = n + k;
        if gram.len() != dimension * dimension || !gram.iter().all(|x| x.is_finite()) {
            return Err(ArticulatedContactError::Invalid(
                "invalid augmented impact mobility",
            ));
        }
        if k == 0 {
            return Ok(Self {
                contact_mobility: gram.to_vec(),
                minimum_scaled_pivot: 1.,
                contact_rows: n,
                bilateral_rows: 0,
                recovery: Vec::new(),
            });
        }
        let mut scale = Vec::with_capacity(k);
        for i in 0..k {
            poll(&mut cancelled)?;
            let d = gram[(n + i) * dimension + n + i];
            let free = reference_mobility[i];
            if !free.is_finite() || free <= 0. || d <= 0. || d / free < minimum_scaled_pivot {
                return Err(ArticulatedContactError::Invalid(
                    "bilateral row has unresolved native mobility",
                ));
            }
            scale.push(d.sqrt());
        }
        let mut a = vec![0.; k * k];
        for i in 0..k {
            a[i * k + i] = 1.;
            for j in 0..i {
                let x = (gram[(n + i) * dimension + n + j] / scale[i]) / scale[j];
                let y = (gram[(n + j) * dimension + n + i] / scale[j]) / scale[i];
                if !x.is_finite() || !y.is_finite() || (x - y).abs() > 1e-10 {
                    return Err(ArticulatedContactError::Unresolved(
                        "bilateral mobility reciprocity",
                    ));
                }
                a[i * k + j] = 0.5 * (x + y);
                a[j * k + i] = a[i * k + j];
            }
        }
        poll(&mut cancelled)?;
        let factor = cholesky(&a, k).map_err(|_| {
            ArticulatedContactError::Invalid("dependent impact bilateral constraints")
        })?;
        let mut minimum_pivot = 1.0_f64;
        for i in 0..k {
            let pivot = factor.l(i, i).powi(2);
            if !pivot.is_finite() || pivot < minimum_scaled_pivot {
                return Err(ArticulatedContactError::Invalid(
                    "ill-conditioned impact bilateral constraints",
                ));
            }
            minimum_pivot = minimum_pivot.min(pivot);
        }
        let mut recovery = vec![0.; n * k];
        for column in 0..n {
            poll(&mut cancelled)?;
            let rhs = &mut recovery[column * k..(column + 1) * k];
            for i in 0..k {
                rhs[i] = gram[(n + i) * dimension + column] / scale[i];
            }
            factor.solve(rhs);
            for i in 0..k {
                rhs[i] = finite(rhs[i] / scale[i])?;
            }
        }
        let mut mobility = vec![0.; n * n];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                let correction: f64 = (0..k)
                    .map(|r| gram[i * dimension + n + r] * recovery[j * k + r])
                    .sum();
                mobility[i * n + j] = finite(gram[i * dimension + j] - correction)?;
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            contact_mobility: mobility,
            minimum_scaled_pivot: minimum_pivot,
            contact_rows: n,
            bilateral_rows: k,
            recovery,
        })
    }

    /// Recover bilateral impulses from the same factor used for condensation.
    /// A contact velocity-only witness does not confer unique reaction loads.
    pub fn recover(
        &self,
        contact_impulses: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, ArticulatedContactError> {
        poll(&mut cancelled)?;
        if contact_impulses.len() != self.contact_rows
            || !contact_impulses.iter().all(|x| x.is_finite())
        {
            return Err(ArticulatedContactError::Invalid(
                "invalid contact impulses for bilateral recovery",
            ));
        }
        let mut result = Vec::with_capacity(self.bilateral_rows);
        for r in 0..self.bilateral_rows {
            poll(&mut cancelled)?;
            result.push(finite(
                -(0..self.contact_rows)
                    .map(|i| self.recovery[i * self.bilateral_rows + r] * contact_impulses[i])
                    .sum::<f64>(),
            )?);
        }
        Ok(result)
    }
}
