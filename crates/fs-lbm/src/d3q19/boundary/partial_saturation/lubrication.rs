//! Missing normal Stokes resistance between separated spherical surfaces.
//!
//! This is a subgrid correction, not a resolved film or a contact law. The
//! gap geometry is frozen at the same pose as the native fluid/body exchange.

use super::{multiple::Response, *};

pub(in crate::d3q19::boundary) mod planes;
pub use planes::{LubricatedPlanesExchange3, LubricationPlane3, PlaneLubricationExchange3};

/// Constant Newtonian normal lubrication correction in lattice units.
///
/// For actual gap `h > 0`, effective radius `R = Ra Rb/(Ra+Rb)` and
/// `h_eff = max(h, roughness_gap)`, the extra resistance is
/// `zeta = 6 pi mu R^2 (1/h_eff - 1/cutoff_gap)` below the cutoff, zero above.
/// This subtracts the resistance assigned to the resolved-fluid cutoff.
/// Dynamic viscosity `mu` has lattice mass/(length*time) units; it is not the
/// kinematic viscosity `(tau-1/2)/3`. The caller owns matching it to the fluid
/// density and calibrating the cutoff for the native partial-saturation grid.
/// The positive roughness scale explicitly regularizes the Stokes singularity.
///
/// All active pairs are solved simultaneously at midpoint, together with the
/// native six-velocity body responses. Midpoint damping is passive but not
/// L-stable: an unresolved stiff relaxation can reverse its endpoint slip.
/// Prescribed planes can join the same solve through the explicit plane entry
/// points. This model does not supply contact, temporal accuracy near contact,
/// tangential lubrication, cavitation, or a high-volume-fraction closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereLubrication3 {
    dynamic_viscosity: f64,
    cutoff_gap: f64,
    roughness_gap: f64,
    maximum_pairs: usize,
}

impl SphereLubrication3 {
    /// Admit finite positive viscosity and `0 < roughness_gap <= cutoff_gap`.
    /// At most 64 bodies and 64 active pairs are supported by this small-
    /// assembly solver; `maximum_pairs` must be in `1..=64` and bounds the TOTAL
    /// number of sphere-pair and sphere-plane films when planes are supplied.
    ///
    /// # Errors
    /// Refuses unsupported scales or pair limits.
    pub fn new(
        dynamic_viscosity: f64,
        cutoff_gap: f64,
        roughness_gap: f64,
        maximum_pairs: usize,
    ) -> Result<Self> {
        if ![dynamic_viscosity, cutoff_gap, roughness_gap]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
            || roughness_gap > cutoff_gap
            || !(1..=64).contains(&maximum_pairs)
        {
            return Err(invalid("invalid sphere lubrication scales or pair limit"));
        }
        Ok(Self {
            dynamic_viscosity,
            cutoff_gap,
            roughness_gap,
            maximum_pairs,
        })
    }

    pub(in crate::d3q19::boundary) fn admit<const N: usize>(self, radii: [f64; N]) -> Result<()> {
        if N == 0 || N > 64 || !radii.into_iter().all(|r| r.is_finite() && r > 0.0) {
            return Err(invalid(
                "lubrication requires 1 through 64 positive sphere radii",
            ));
        }
        Ok(())
    }
}

/// Actual equal/opposite gap impulses and their midpoint work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereLubricationExchange3<const N: usize> {
    /// Center-normal impulse on each body; no spin impulse is added.
    pub body_impulses: [[f64; 3]; N],
    /// Orbital angular impulse about the caller's fixed world origin.
    pub body_angular_impulses: [[f64; 3]; N],
    /// Work using the ACTUAL returned midpoint velocities and net body impulse.
    pub body_work: [f64; N],
    /// Sum of missing-film reaction impulses received by prescribed planes.
    /// These are separate from the original resolved-fluid wall exchange.
    pub plane_impulse: [f64; 3],
    /// Total plane reaction angular impulse about the same fixed world origin.
    pub plane_angular_impulse: [f64; 3],
    /// Actual work received by prescribed translating planes.
    pub plane_work: f64,
    /// Sum of `zeta * relative_midpoint_normal_velocity^2`, nonnegative.
    /// This is viscous model loss; no thermal destination is implied.
    pub dissipated_energy: f64,
    /// Actual body work plus plane work plus viscous loss, measuring roundoff.
    pub work_residual: f64,
    /// Maximum absolute pair residual `impulse + zeta * midpoint_slip`.
    pub constitutive_residual: f64,
    /// Number of positive-resistance sphere pairs in the simultaneous solve.
    pub active_pairs: usize,
    /// Number of positive-resistance sphere-plane films in that SAME solve.
    pub active_plane_pairs: usize,
}

impl<const N: usize> Default for SphereLubricationExchange3<N> {
    fn default() -> Self {
        Self {
            body_impulses: [[0.0; 3]; N],
            body_angular_impulses: [[0.0; 3]; N],
            body_work: [0.0; N],
            plane_impulse: [0.0; 3],
            plane_angular_impulse: [0.0; 3],
            plane_work: 0.0,
            dissipated_energy: 0.0,
            work_residual: 0.0,
            constitutive_residual: 0.0,
            active_pairs: 0,
            active_plane_pairs: 0,
        }
    }
}

impl BoundaryGrid3 {
    /// Advance spherical frozen geometries with missing normal film resistance.
    ///
    /// Fraction fields must describe the disjoint spheres of the supplied radii
    /// about their body centers. The geometry-free caller owns that obligation;
    /// [`Self::step_with_lubricated_spheres`] samples the actual spheres instead.
    /// Positive actual gaps are required even below the admitted roughness scale.
    /// Guo fluid forces and solid external impulses retain their original owners.
    ///
    /// # Errors
    /// Inherits native assembly admission and rejects intersecting spheres,
    /// invalid radii, exhausted pair limits or unrepresentable spatial solves.
    /// No population or body state changes on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_partial_saturation_bodies<const N: usize>(
        &mut self,
        fractions: [&[f64]; N],
        radii: [f64; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<PartialSaturationBodiesExchange3<N>> {
        self.step_with_partial_saturation_resistance_checked(
            fractions,
            bodies,
            external,
            origin,
            gate,
            Some((model, radii, &[])),
            |_, _, _| Ok(()),
        )
        .map(|(report, (), _)| report)
    }
}

#[derive(Clone, Copy)]
struct Pair {
    a: usize,
    counterpart: Counterpart,
    normal: [f64; 3],
    resistance: f64,
}

#[derive(Clone, Copy)]
enum Counterpart {
    Sphere(usize),
    Plane { index: usize, point: [f64; 3] },
}

impl Pair {
    fn bodies(self) -> [Option<(usize, f64)>; 2] {
        [
            Some((self.a, 1.0)),
            match self.counterpart {
                Counterpart::Sphere(b) => Some((b, -1.0)),
                Counterpart::Plane { .. } => None,
            },
        ]
    }
}

#[derive(Clone, Copy, Default)]
struct FreeResponse {
    increment: [f64; 6],
    // Only the three linear-impulse columns of the native six-velocity map.
    linear_columns: [[f64; 3]; 6],
}

pub(super) fn solve<const N: usize>(
    model: SphereLubrication3,
    radii: [f64; N],
    planes: &[LubricationPlane3],
    bodies: &[PartialSaturationBody3; N],
    response: &[Response; N],
    proposed: &mut [PartialSaturationBody3; N],
    origin: [f64; 3],
    gate: &CancelGate,
) -> Result<(
    SphereLubricationExchange3<N>,
    Vec<PlaneLubricationExchange3>,
)> {
    let mut pairs = Vec::new();
    pairs
        .try_reserve_exact(model.maximum_pairs)
        .map_err(|_| invalid("lubrication pair allocation refused"))?;
    for a in 0..N {
        checkpoint(gate)?;
        for b in a + 1..N {
            let delta = core::array::from_fn(|i| bodies[b].center[i] - bodies[a].center[i]);
            let distance = delta[0].hypot(delta[1]).hypot(delta[2]);
            let gap = (distance - radii[a]) - radii[b];
            if !(gap.is_finite() && gap > 0.0) {
                return Err(invalid(
                    "lubrication requires strictly positive actual sphere gaps",
                ));
            }
            let effective_gap = gap.max(model.roughness_gap);
            if effective_gap >= model.cutoff_gap {
                continue;
            }
            if pairs.len() == model.maximum_pairs {
                return Err(invalid("lubrication active-pair limit exhausted"));
            }
            let small = radii[a].min(radii[b]);
            let large = radii[a].max(radii[b]);
            let radius = small / (1.0 + small / large);
            let resistance = 6.0
                * core::f64::consts::PI
                * model.dynamic_viscosity
                * radius
                * radius
                * (1.0 / effective_gap - 1.0 / model.cutoff_gap);
            if !(resistance.is_finite() && resistance > 0.0) {
                return Err(invalid("lubrication resistance is unrepresentable"));
            }
            pairs.push(Pair {
                a,
                counterpart: Counterpart::Sphere(b),
                normal: delta.map(|x| x / distance),
                resistance,
            });
        }
    }
    let sphere_pairs = pairs.len();
    let mut plane_reports = Vec::new();
    plane_reports
        .try_reserve_exact(planes.len())
        .map_err(|_| invalid("lubrication plane report allocation refused"))?;
    for (index, plane) in planes.iter().copied().enumerate() {
        checkpoint(gate)?;
        plane_reports.push(PlaneLubricationExchange3::default());
        for a in 0..N {
            let distance = plane.distance(bodies[a].center);
            let gap = distance - radii[a];
            if !(gap.is_finite() && gap > 0.0) {
                return Err(invalid(
                    "lubrication requires strictly positive actual plane gaps",
                ));
            }
            let effective_gap = gap.max(plane.roughness_gap);
            if effective_gap >= plane.cutoff_gap {
                continue;
            }
            if pairs.len() == model.maximum_pairs {
                return Err(invalid("lubrication active-pair limit exhausted"));
            }
            let resistance = 6.0
                * core::f64::consts::PI
                * model.dynamic_viscosity
                * radii[a]
                * radii[a]
                * (1.0 / effective_gap - 1.0 / plane.cutoff_gap);
            let point = core::array::from_fn(|i| bodies[a].center[i] - distance * plane.normal[i]);
            if !(resistance.is_finite() && resistance > 0.0)
                || !point.into_iter().all(f64::is_finite)
            {
                return Err(invalid(
                    "lubrication plane resistance or footpoint overflowed",
                ));
            }
            pairs.push(Pair {
                a,
                counterpart: Counterpart::Plane { index, point },
                normal: plane.normal,
                resistance,
            });
            plane_reports[index].active_spheres += 1;
        }
    }
    if pairs.is_empty() {
        for a in 0..N {
            checkpoint(gate)?;
            proposed[a] = response[a].solve(bodies[a])?;
        }
        return Ok((SphereLubricationExchange3::default(), plane_reports));
    }
    let mut free = [FreeResponse::default(); N];
    for a in 0..N {
        checkpoint(gate)?;
        let factor = response[a].factor(bodies[a])?;
        free[a].increment = response[a].increment();
        factor.solve(&mut free[a].increment);
        for axis in 0..3 {
            let mut rhs = [0.0; 6];
            rhs[axis] = 1.0;
            factor.solve(&mut rhs);
            for (row, value) in rhs.into_iter().enumerate() {
                free[a].linear_columns[row][axis] = value;
            }
        }
    }
    let k = pairs.len();
    let mut matrix = Vec::new();
    matrix
        .try_reserve_exact(k * k)
        .map_err(|_| invalid("lubrication response allocation refused"))?;
    matrix.resize(k * k, 0.0);
    let mut impulse = [0.0; 64];
    for (i, pair) in pairs.iter().copied().enumerate() {
        checkpoint(gate)?;
        let slip = core::array::from_fn(|axis| match pair.counterpart {
            Counterpart::Sphere(b) => {
                bodies[pair.a].velocity[axis] - bodies[b].velocity[axis]
                    + 0.5 * (free[pair.a].increment[axis] - free[b].increment[axis])
            }
            Counterpart::Plane { index, .. } => {
                bodies[pair.a].velocity[axis] - planes[index].velocity[axis]
                    + 0.5 * free[pair.a].increment[axis]
            }
        });
        impulse[i] = -dot(pair.normal, slip);
        for (j, other) in pairs[..=i].iter().copied().enumerate() {
            let mut entry = if i == j { 1.0 / pair.resistance } else { 0.0 };
            for (a, sign_a) in pair.bodies().into_iter().flatten() {
                for (b, sign_b) in other.bodies().into_iter().flatten() {
                    if a == b {
                        let velocity = core::array::from_fn(|axis| {
                            dot(free[a].linear_columns[axis], other.normal)
                        });
                        entry += 0.5 * sign_a * sign_b * dot(pair.normal, velocity);
                    }
                }
            }
            matrix[i * k + j] = entry;
        }
    }
    if !matrix
        .iter()
        .chain(impulse[..k].iter())
        .all(|x| x.is_finite())
    {
        return Err(invalid("lubrication coupled response overflowed"));
    }
    if k > 0 {
        fs_la::factor::cholesky(&matrix, k)
            .map_err(|_| invalid("lubrication response is not numerically positive definite"))?
            .solve(&mut impulse[..k]);
    }
    let mut report = SphereLubricationExchange3 {
        active_pairs: sphere_pairs,
        active_plane_pairs: k - sphere_pairs,
        ..Default::default()
    };
    for (pair, lambda) in pairs.iter().zip(impulse) {
        for axis in 0..3 {
            let value = lambda * pair.normal[axis];
            report.body_impulses[pair.a][axis] += value;
            match pair.counterpart {
                Counterpart::Sphere(b) => report.body_impulses[b][axis] -= value,
                Counterpart::Plane { index, .. } => plane_reports[index].impulse[axis] -= value,
            }
        }
        if let Counterpart::Plane { index, point } = pair.counterpart {
            let reaction = pair.normal.map(|x| -lambda * x);
            let angular = cross(core::array::from_fn(|i| point[i] - origin[i]), reaction);
            for (total, value) in plane_reports[index].angular_impulse.iter_mut().zip(angular) {
                *total += value;
            }
        }
    }
    for a in 0..N {
        checkpoint(gate)?;
        let increment = core::array::from_fn(|axis| {
            free[a].increment[axis] + dot(free[a].linear_columns[axis], report.body_impulses[a])
        });
        proposed[a] = Response::apply(bodies[a], increment)?;
        let midpoint =
            core::array::from_fn(|i| 0.5 * (bodies[a].velocity[i] + proposed[a].velocity[i]));
        report.body_work[a] = dot(report.body_impulses[a], midpoint);
        report.body_angular_impulses[a] = cross(
            core::array::from_fn(|i| bodies[a].center[i] - origin[i]),
            report.body_impulses[a],
        );
    }
    for (pair, lambda) in pairs.iter().zip(impulse) {
        let slip = dot(
            pair.normal,
            core::array::from_fn(|i| match pair.counterpart {
                Counterpart::Sphere(b) => {
                    0.5 * ((bodies[pair.a].velocity[i] + proposed[pair.a].velocity[i])
                        - (bodies[b].velocity[i] + proposed[b].velocity[i]))
                }
                Counterpart::Plane { index, .. } => {
                    0.5 * ((bodies[pair.a].velocity[i] - planes[index].velocity[i])
                        + (proposed[pair.a].velocity[i] - planes[index].velocity[i]))
                }
            }),
        );
        report.dissipated_energy += pair.resistance * slip * slip;
        report.constitutive_residual = report
            .constitutive_residual
            .max((lambda + pair.resistance * slip).abs());
    }
    for (plane, exchange) in planes.iter().zip(&mut plane_reports) {
        exchange.work = dot(exchange.impulse, plane.velocity);
        report.plane_work += exchange.work;
        for i in 0..3 {
            report.plane_impulse[i] += exchange.impulse[i];
            report.plane_angular_impulse[i] += exchange.angular_impulse[i];
        }
    }
    report.work_residual =
        report.body_work.into_iter().sum::<f64>() + report.plane_work + report.dissipated_energy;
    if !report
        .body_impulses
        .into_iter()
        .flatten()
        .chain(report.body_angular_impulses.into_iter().flatten())
        .chain(report.body_work)
        .chain(report.plane_impulse)
        .chain(report.plane_angular_impulse)
        .chain([
            report.plane_work,
            report.dissipated_energy,
            report.work_residual,
            report.constitutive_residual,
        ])
        .all(f64::is_finite)
    {
        return Err(invalid("lubrication momentum/work account overflowed"));
    }
    Ok((report, plane_reports))
}
