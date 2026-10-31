//! Missing Stokes resistance between separated spherical surfaces and planes.
//!
//! This is a subgrid correction, not a resolved film or a contact law. The
//! gap geometry is frozen at the same pose as the native fluid/body exchange.

use super::{multiple::Response, *};

mod pairs;
pub(in crate::d3q19::boundary) mod planes;
pub use planes::{LubricatedPlanesExchange3, LubricationPlane3, PlaneLubricationExchange3};

/// Constant Newtonian missing-film resistance in lattice units.
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
/// points. Both pairs and planes can opt into reciprocal tangent/spin terms.
/// This model does not supply contact, temporal accuracy near contact,
/// complete finite-gap hydrodynamics, cavitation, or a dense-suspension closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereLubrication3 {
    dynamic_viscosity: f64,
    cutoff_gap: f64,
    roughness_gap: f64,
    maximum_pairs: usize,
    maximum_modes: usize,
    tangential_pairs: Option<(f64, f64)>,
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
            maximum_modes: 64,
            tangential_pairs: None,
        })
    }

    /// Enable unequal-sphere leading-logarithmic tangent/parallel-spin films.
    ///
    /// With radii `a,b`, actual center distance `d`, tangent `t` and
    /// `s = n cross t`, set `u = (v_a-v_b).t`, `z_a=a*(omega_a.s+u/d)`
    /// and `z_b=b*(omega_b.s+u/d)`. The resistance in these two velocities is
    /// `k * [[4,1],[1,4]]`, `k=4*pi*mu*(a*b/(a+b))/5*ln(cutoff/h_eff)`.
    /// Both positive modes are retained: zero surface slip alone does not remove
    /// the counter-rotation resistance. Rotation about `n` has no divergent
    /// leading term and receives no extra resistance.
    ///
    /// This is the objective finite-gap completion of the Jeffrey--Onishi
    /// leading logarithmic `YA/YB/YC` block, with corrected normalization from
    /// Townsend, Phys. Fluids 35, 127126 (2023), equations 19-25 and 33-34a:
    /// <https://arxiv.org/abs/1802.08226>. Actual `d`, rather than `a+b`, makes
    /// the frame follow the line of centers. The resulting power-conjugate
    /// forces/torques conserve orbital-plus-spin momentum and vanish under a
    /// common rigid motion. This changes only omitted-order `h*ln(h)` terms;
    /// it does NOT reproduce the complete higher-order resistance functions.
    ///
    /// Cutoff/roughness are separate lattice lengths from the normal film. The
    /// cutoff must be <= 0.02 times every supplied radius, a declared small-gap
    /// domain, not a certified error bound. Actual gaps must remain positive.
    /// Four modes per active tangent film join the same native response solve.
    ///
    /// # Errors
    /// Refuses nonfinite/nonpositive scales or roughness above the cutoff.
    pub fn with_tangential_pairs(mut self, cutoff_gap: f64, roughness_gap: f64) -> Result<Self> {
        if !(cutoff_gap.is_finite() && cutoff_gap > 0.0)
            || !(roughness_gap.is_finite() && roughness_gap > 0.0)
            || roughness_gap > cutoff_gap
        {
            return Err(invalid("invalid tangential sphere-pair lubrication scales"));
        }
        self.tangential_pairs = Some((cutoff_gap, roughness_gap));
        Ok(self)
    }

    /// Bound the simultaneous scalar resistance solve independently of the
    /// physical-film count. The default is 64; an explicit limit in `1..=320`
    /// admits at most five modes per each of the 64 supported films. A normal
    /// film uses one mode and an active tangent pair/wall block uses four more.
    /// Dense storage is O(modes squared), factorization O(modes cubed).
    ///
    /// # Errors
    /// Refuses zero or a mode limit above 320.
    pub fn with_resistance_mode_limit(mut self, maximum_modes: usize) -> Result<Self> {
        if !(1..=320).contains(&maximum_modes) {
            return Err(invalid("invalid lubrication resistance-mode limit"));
        }
        self.maximum_modes = maximum_modes;
        Ok(self)
    }

    /// Refresh the same fluid's viscosity while preserving film geometry,
    /// regularization, tangential terms and work limits.
    ///
    /// This changes no accepted body or fluid history. A material adapter must
    /// update resolved-fluid viscosity from the same physical snapshot before
    /// the next coupled step; this method does not infer that correspondence.
    ///
    /// # Errors
    /// Refuses nonpositive or nonfinite dynamic viscosity before mutation.
    pub fn with_dynamic_viscosity(mut self, dynamic_viscosity: f64) -> Result<Self> {
        if !dynamic_viscosity.is_finite() || dynamic_viscosity <= 0.0 {
            return Err(invalid("invalid sphere lubrication viscosity"));
        }
        self.dynamic_viscosity = dynamic_viscosity;
        Ok(self)
    }

    pub(in crate::d3q19::boundary) fn admit<const N: usize>(self, radii: [f64; N]) -> Result<()> {
        if N == 0 || N > 64 || !radii.into_iter().all(|r| r.is_finite() && r > 0.0) {
            return Err(invalid(
                "lubrication requires 1 through 64 positive sphere radii",
            ));
        }
        if let Some((cutoff, _)) = self.tangential_pairs
            && radii.into_iter().any(|r| cutoff > 0.02 * r)
        {
            return Err(invalid(
                "tangential sphere-pair cutoff exceeds 0.02 sphere radius",
            ));
        }
        Ok(())
    }
}

/// Actual equal/opposite gap impulses and their midpoint work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereLubricationExchange3<const N: usize> {
    /// Net linear film impulse on each body.
    pub body_impulses: [[f64; 3]; N],
    /// Net film spin impulse about each body's center.
    pub body_spin_impulses: [[f64; 3]; N],
    /// Orbital PLUS spin angular impulse about the caller's fixed world origin.
    pub body_angular_impulses: [[f64; 3]; N],
    /// Full wrench work using ACTUAL returned midpoint translation and spin.
    pub body_work: [f64; N],
    /// Sum of missing-film reaction impulses received by prescribed planes.
    /// These are separate from the original resolved-fluid wall exchange.
    pub plane_impulse: [f64; 3],
    /// Total plane reaction angular impulse about the same fixed world origin.
    pub plane_angular_impulse: [f64; 3],
    /// Actual work received by prescribed translating planes.
    pub plane_work: f64,
    /// Sum of `resistance * midpoint_mode_speed^2`, nonnegative. Tangential
    /// modes use radius-scaled spin so their speeds also have velocity units.
    /// This is viscous model loss; no thermal destination is implied.
    pub dissipated_energy: f64,
    /// Actual body work plus plane work plus viscous loss, measuring roundoff.
    pub work_residual: f64,
    /// Maximum absolute scalar mode residual `impulse + resistance * speed`.
    pub constitutive_residual: f64,
    /// Number of positive-resistance sphere pairs in the simultaneous solve.
    pub active_pairs: usize,
    /// Active sphere pairs with reciprocal tangent/parallel-spin resistance.
    pub active_tangential_pairs: usize,
    /// Number of positive-resistance sphere-plane films in that SAME solve.
    pub active_plane_pairs: usize,
    /// Number of active plane films with translation/parallel-spin correction.
    pub active_tangential_plane_pairs: usize,
    /// Actual scalar dimension of the shared resistance solve.
    pub active_resistance_modes: usize,
}

impl<const N: usize> Default for SphereLubricationExchange3<N> {
    fn default() -> Self {
        Self {
            body_impulses: [[0.0; 3]; N],
            body_spin_impulses: [[0.0; 3]; N],
            body_angular_impulses: [[0.0; 3]; N],
            body_work: [0.0; N],
            plane_impulse: [0.0; 3],
            plane_angular_impulse: [0.0; 3],
            plane_work: 0.0,
            dissipated_energy: 0.0,
            work_residual: 0.0,
            constitutive_residual: 0.0,
            active_pairs: 0,
            active_tangential_pairs: 0,
            active_plane_pairs: 0,
            active_tangential_plane_pairs: 0,
            active_resistance_modes: 0,
        }
    }
}

impl BoundaryGrid3 {
    /// Advance spherical frozen geometries with missing film resistance.
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
    // Zero for the original center-normal modes. A tangent mode has linear
    // covector `normal` and angular covector `angular` about the sphere center.
    angular: Option<[f64; 3]>,
    resistance: f64,
}

#[derive(Clone, Copy)]
enum Counterpart {
    // The second body's entire wrench has sign -1. Its angular covector need
    // not equal the first body's: pair rolling modes require independent spins.
    Sphere {
        index: usize,
        angular: Option<[f64; 3]>,
    },
    Plane {
        index: usize,
        point: [f64; 3],
    },
}

#[derive(Clone, Copy)]
struct Direction {
    linear: [f64; 3],
    angular: Option<[f64; 3]>,
}

fn tangent_basis(normal: [f64; 3]) -> [[f64; 3]; 2] {
    // Isotropic tangent planes share a deterministic, well-conditioned gauge.
    let axis = (0..3)
        .min_by(|&a, &b| normal[a].abs().total_cmp(&normal[b].abs()))
        .expect("three axes");
    let mut hint = [0.0; 3];
    hint[axis] = 1.0;
    let tangent = cross(normal, hint);
    let length = tangent[0].hypot(tangent[1]).hypot(tangent[2]);
    let tangent = tangent.map(|x| x / length);
    [tangent, cross(normal, tangent)]
}

impl Pair {
    fn bodies(self) -> [Option<(usize, f64, Direction)>; 2] {
        [
            Some((
                self.a,
                1.0,
                Direction {
                    linear: self.normal,
                    angular: self.angular,
                },
            )),
            match self.counterpart {
                Counterpart::Sphere { index, angular } => Some((
                    index,
                    -1.0,
                    Direction {
                        linear: self.normal,
                        angular,
                    },
                )),
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
    // Only bodies touched by active tangential modes need these three columns.
    angular_columns: Option<[[f64; 3]; 6]>,
}

impl FreeResponse {
    fn action(&self, row: usize, mode: Direction) -> f64 {
        let linear = dot(self.linear_columns[row], mode.linear);
        match mode.angular {
            Some(angular) => {
                linear
                    + dot(
                        self.angular_columns.as_ref().expect("active spin response")[row],
                        angular,
                    )
            }
            None => linear,
        }
    }

    fn mobility(&self, test: Direction, trial: Direction) -> f64 {
        let linear = dot(test.linear, core::array::from_fn(|i| self.action(i, trial)));
        match test.angular {
            Some(angular) => {
                linear + dot(angular, core::array::from_fn(|i| self.action(i + 3, trial)))
            }
            None => linear,
        }
    }
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
    let modes_per_film = if model.tangential_pairs.is_some()
        || planes.iter().any(|plane| plane.tangential.is_some())
    {
        5
    } else {
        1
    };
    pairs
        .try_reserve_exact(
            model
                .maximum_modes
                .min(model.maximum_pairs * modes_per_film),
        )
        .map_err(|_| invalid("lubrication pair allocation refused"))?;
    let mut sphere_pairs = 0;
    let mut tangential_sphere_pairs = 0;
    let mut spin_response = [false; N];
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
            let normal_active = effective_gap < model.cutoff_gap;
            let normal = delta.map(|x| x / distance);
            let tangential = pairs::tangential_modes(
                [radii[a], radii[b]],
                distance,
                normal,
                gap,
                model.tangential_pairs,
                model.dynamic_viscosity,
            )?;
            if !normal_active && tangential.is_none() {
                continue;
            }
            if sphere_pairs == model.maximum_pairs {
                return Err(invalid("lubrication active-pair limit exhausted"));
            }
            let added_modes = usize::from(normal_active) + if tangential.is_some() { 4 } else { 0 };
            if pairs.len() + added_modes > model.maximum_modes {
                return Err(invalid("lubrication resistance-mode limit exhausted"));
            }
            if normal_active {
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
                    counterpart: Counterpart::Sphere {
                        index: b,
                        angular: None,
                    },
                    normal,
                    angular: None,
                    resistance,
                });
            }
            if let Some(tangential) = tangential {
                for mode in tangential {
                    pairs.push(Pair {
                        a,
                        counterpart: Counterpart::Sphere {
                            index: b,
                            angular: Some(mode.angular_b.map(|x| -x)),
                        },
                        normal: mode.linear,
                        angular: Some(mode.angular_a),
                        resistance: mode.resistance,
                    });
                }
                spin_response[a] = true;
                spin_response[b] = true;
                tangential_sphere_pairs += 1;
            }
            sphere_pairs += 1;
        }
    }
    let mut plane_pairs = 0;
    let mut tangential_plane_pairs = 0;
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
            let normal_active = effective_gap < plane.cutoff_gap;
            let tangential = plane.tangential_modes(radii[a], gap, model.dynamic_viscosity)?;
            if !normal_active && tangential.is_none() {
                continue;
            }
            if sphere_pairs + plane_pairs == model.maximum_pairs {
                return Err(invalid("lubrication active-pair limit exhausted"));
            }
            let added_modes = usize::from(normal_active) + if tangential.is_some() { 4 } else { 0 };
            if pairs.len() + added_modes > model.maximum_modes {
                return Err(invalid("lubrication resistance-mode limit exhausted"));
            }
            let point = core::array::from_fn(|i| bodies[a].center[i] - distance * plane.normal[i]);
            if !point.into_iter().all(f64::is_finite) {
                return Err(invalid("lubrication plane footpoint overflowed"));
            }
            if normal_active {
                let resistance = 6.0
                    * core::f64::consts::PI
                    * model.dynamic_viscosity
                    * radii[a]
                    * radii[a]
                    * (1.0 / effective_gap - 1.0 / plane.cutoff_gap);
                if !(resistance.is_finite() && resistance > 0.0) {
                    return Err(invalid("lubrication plane resistance overflowed"));
                }
                pairs.push(Pair {
                    a,
                    counterpart: Counterpart::Plane { index, point },
                    normal: plane.normal,
                    angular: None,
                    resistance,
                });
            }
            if let Some(tangential) = tangential {
                for mode in tangential {
                    pairs.push(Pair {
                        a,
                        counterpart: Counterpart::Plane { index, point },
                        normal: mode.linear,
                        angular: Some(mode.angular),
                        resistance: mode.resistance,
                    });
                }
                spin_response[a] = true;
                tangential_plane_pairs += 1;
                plane_reports[index].tangential_spheres += 1;
            }
            plane_pairs += 1;
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
        if spin_response[a] {
            let mut columns = [[0.0; 3]; 6];
            for axis in 0..3 {
                let mut rhs = [0.0; 6];
                rhs[axis + 3] = 1.0;
                factor.solve(&mut rhs);
                for (row, value) in rhs.into_iter().enumerate() {
                    columns[row][axis] = value;
                }
            }
            free[a].angular_columns = Some(columns);
        }
    }
    let k = pairs.len();
    let mut matrix = Vec::new();
    matrix
        .try_reserve_exact(k * k)
        .map_err(|_| invalid("lubrication response allocation refused"))?;
    matrix.resize(k * k, 0.0);
    let mut impulse = [0.0; 320];
    for (i, pair) in pairs.iter().copied().enumerate() {
        checkpoint(gate)?;
        let slip = core::array::from_fn(|axis| match pair.counterpart {
            Counterpart::Sphere { index: b, .. } => {
                bodies[pair.a].velocity[axis] - bodies[b].velocity[axis]
                    + 0.5 * (free[pair.a].increment[axis] - free[b].increment[axis])
            }
            Counterpart::Plane { index, .. } => {
                bodies[pair.a].velocity[axis] - planes[index].velocity[axis]
                    + 0.5 * free[pair.a].increment[axis]
            }
        });
        impulse[i] = -dot(pair.normal, slip);
        if let Some(angular) = pair.angular {
            let omega = core::array::from_fn(|axis| {
                bodies[pair.a].angular_velocity[axis] + 0.5 * free[pair.a].increment[axis + 3]
            });
            impulse[i] -= dot(angular, omega);
        }
        if let Counterpart::Sphere {
            index: b,
            angular: Some(angular),
        } = pair.counterpart
        {
            let omega = core::array::from_fn(|axis| {
                bodies[b].angular_velocity[axis] + 0.5 * free[b].increment[axis + 3]
            });
            impulse[i] += dot(angular, omega);
        }
        for (j, other) in pairs[..=i].iter().copied().enumerate() {
            let mut entry = if i == j { 1.0 / pair.resistance } else { 0.0 };
            for (a, sign_a, direction_a) in pair.bodies().into_iter().flatten() {
                for (b, sign_b, direction_b) in other.bodies().into_iter().flatten() {
                    if a == b {
                        entry += 0.5 * sign_a * sign_b * free[a].mobility(direction_a, direction_b);
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
        active_tangential_pairs: tangential_sphere_pairs,
        active_plane_pairs: plane_pairs,
        active_tangential_plane_pairs: tangential_plane_pairs,
        active_resistance_modes: k,
        ..Default::default()
    };
    for (pair, lambda) in pairs.iter().zip(impulse) {
        for axis in 0..3 {
            let value = lambda * pair.normal[axis];
            report.body_impulses[pair.a][axis] += value;
            match pair.counterpart {
                Counterpart::Sphere { index: b, .. } => report.body_impulses[b][axis] -= value,
                Counterpart::Plane { index, .. } => plane_reports[index].impulse[axis] -= value,
            }
        }
        if let Some(angular) = pair.angular {
            for i in 0..3 {
                report.body_spin_impulses[pair.a][i] += lambda * angular[i];
            }
        }
        if let Counterpart::Sphere {
            index: b,
            angular: Some(angular),
        } = pair.counterpart
        {
            for i in 0..3 {
                report.body_spin_impulses[b][i] -= lambda * angular[i];
            }
        }
        if let Counterpart::Plane { index, point } = pair.counterpart {
            let reaction = pair.normal.map(|x| -lambda * x);
            let mut angular = cross(core::array::from_fn(|i| point[i] - origin[i]), reaction);
            if let Some(spin_axis) = pair.angular {
                let spin = spin_axis.map(|x| lambda * x);
                let lever = cross(
                    core::array::from_fn(|i| bodies[pair.a].center[i] - point[i]),
                    reaction,
                );
                for i in 0..3 {
                    // Opposite total body moment, expressed as force at the
                    // real wall point PLUS a free wall couple. Treating this as
                    // contact-point friction alone would lose the pressure moment.
                    let couple = lever[i] - spin[i];
                    plane_reports[index].couple_impulse[i] += couple;
                    angular[i] += couple;
                }
            }
            for (total, value) in plane_reports[index].angular_impulse.iter_mut().zip(angular) {
                *total += value;
            }
        }
    }
    for a in 0..N {
        checkpoint(gate)?;
        let increment = core::array::from_fn(|axis| {
            let linear = free[a].increment[axis]
                + dot(free[a].linear_columns[axis], report.body_impulses[a]);
            match &free[a].angular_columns {
                Some(columns) => linear + dot(columns[axis], report.body_spin_impulses[a]),
                None => linear,
            }
        });
        proposed[a] = Response::apply(bodies[a], increment)?;
        let midpoint =
            core::array::from_fn(|i| 0.5 * (bodies[a].velocity[i] + proposed[a].velocity[i]));
        report.body_work[a] = dot(report.body_impulses[a], midpoint);
        report.body_angular_impulses[a] = cross(
            core::array::from_fn(|i| bodies[a].center[i] - origin[i]),
            report.body_impulses[a],
        );
        if spin_response[a] {
            let omega = core::array::from_fn(|i| {
                0.5 * (bodies[a].angular_velocity[i] + proposed[a].angular_velocity[i])
            });
            report.body_work[a] += dot(report.body_spin_impulses[a], omega);
            for i in 0..3 {
                report.body_angular_impulses[a][i] += report.body_spin_impulses[a][i];
            }
        }
    }
    for (pair, lambda) in pairs.iter().zip(impulse) {
        let mut slip = dot(
            pair.normal,
            core::array::from_fn(|i| match pair.counterpart {
                Counterpart::Sphere { index: b, .. } => {
                    0.5 * ((bodies[pair.a].velocity[i] + proposed[pair.a].velocity[i])
                        - (bodies[b].velocity[i] + proposed[b].velocity[i]))
                }
                Counterpart::Plane { index, .. } => {
                    0.5 * ((bodies[pair.a].velocity[i] - planes[index].velocity[i])
                        + (proposed[pair.a].velocity[i] - planes[index].velocity[i]))
                }
            }),
        );
        if let Some(angular) = pair.angular {
            slip += dot(
                angular,
                core::array::from_fn(|i| {
                    0.5 * (bodies[pair.a].angular_velocity[i]
                        + proposed[pair.a].angular_velocity[i])
                }),
            );
        }
        if let Counterpart::Sphere {
            index: b,
            angular: Some(angular),
        } = pair.counterpart
        {
            slip -= dot(
                angular,
                core::array::from_fn(|i| {
                    0.5 * (bodies[b].angular_velocity[i] + proposed[b].angular_velocity[i])
                }),
            );
        }
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
        .chain(report.body_spin_impulses.into_iter().flatten())
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
