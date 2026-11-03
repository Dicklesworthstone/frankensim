//! Finite surface recoil inside the original simultaneous film response.

use super::*;

/// A locally planar film surface carried by finite spatial inertia.
///
/// Geometry, center and world inertia are frozen during this exchange, just
/// as for [`PartialSaturationBody3`]. The complete center/spin velocity is
/// solved with the spheres and fluid. The plane supplies only the existing
/// missing-film correction: no resolved wall geometry, contact law, thermal
/// destination, pose drift or gyroscopic rotation is introduced.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteLubricationPlane3 {
    geometry: LubricationPlane3,
    body: PartialSaturationBody3,
}

impl FiniteLubricationPlane3 {
    /// Bind a plane to a positive mass and full world-frame inertia. The
    /// plane point need not coincide with the body's center of mass.
    /// All quantities have the native lattice units (`dx = dt = 1`).
    ///
    /// # Errors
    /// Refuses invalid plane scales, nonfinite state, nonpositive mass,
    /// nonsymmetric/non-SPD inertia or inadmissible center velocity.
    pub fn new(
        point: [f64; 3],
        normal: [f64; 3],
        body: PartialSaturationBody3,
        cutoff_gap: f64,
        roughness_gap: f64,
    ) -> Result<Self> {
        body_factor(body)?;
        let geometry =
            LubricationPlane3::new(point, normal, body.velocity, cutoff_gap, roughness_gap)?;
        Ok(Self { geometry, body })
    }

    /// Enable the original leading-logarithmic tangent/parallel-spin law.
    /// Both surface and sphere spin enter relative rigid-frame velocities;
    /// a common translation or rigid rotation produces no film loss.
    ///
    /// # Errors
    /// Inherits [`LubricationPlane3::with_tangential`] scale admission.
    pub fn with_tangential(mut self, cutoff_gap: f64, roughness_gap: f64) -> Result<Self> {
        self.geometry = self.geometry.with_tangential(cutoff_gap, roughness_gap)?;
        Ok(self)
    }

    // Keep finite state packaging shared by the mixed-surface entry point.
    pub(in crate::d3q19::boundary) fn driven(
        self,
        external: [f64; 6],
    ) -> Result<LubricationPlane3> {
        if !external.into_iter().all(f64::is_finite) {
            return Err(invalid("nonfinite finite-surface external impulse"));
        }
        Ok(LubricationPlane3 {
            recoil: Some(FinitePlaneInput {
                body: self.body,
                external,
                carrier_anchor: None,
            }),
            ..self.geometry
        })
    }

    // Called only after the common owner has checked and published its candidate.
    pub(super) fn accept(&mut self, after: PartialSaturationBody3) {
        self.body = after;
        self.geometry.velocity = after.velocity;
    }

    /// Accepted spatial body state; geometry and inertia are unchanged by a tick.
    pub fn body(self) -> PartialSaturationBody3 {
        self.body
    }

    /// Point on the frozen planar surface.
    pub fn point(self) -> [f64; 3] {
        self.geometry.point
    }

    /// Unit normal into the spheres' permitted fluid half-space.
    pub fn normal(self) -> [f64; 3] {
        self.geometry.normal
    }
}

/// Full finite-plane response to its film wrench and authored external impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinitePlaneResponse3 {
    /// Accepted center/spin velocities at the unchanged exchange geometry.
    pub body_after: PartialSaturationBody3,
    /// Actual translational midpoint velocity.
    pub velocity_midpoint: [f64; 3],
    /// Actual angular midpoint velocity.
    pub angular_velocity_midpoint: [f64; 3],
    /// Film angular impulse about the plane body's center, accumulated directly
    /// from local arms rather than subtracting distant world-origin moments.
    pub spin_impulse: [f64; 3],
    /// Actual change of translational plus rotational kinetic energy.
    pub kinetic_energy_change: f64,
    /// Work of the external center impulse and center couple.
    pub external_work: f64,
    /// Linear and center-angular impulse balance roundoff.
    pub momentum_residual: [f64; 6],
    /// Kinetic increment minus external work minus received film work.
    pub work_residual: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct FinitePlaneInput {
    pub body: PartialSaturationBody3,
    pub external: [f64; 6],
    // Private shared-carrier bindings name their FIRST face. Equal body values
    // alone never merge inertia. Legacy surfaces remain independent with None.
    pub carrier_anchor: Option<usize>,
}

impl BoundaryGrid3 {
    /// Exchange resolved fluid and unresolved film momentum with finite surfaces.
    ///
    /// Sphere bodies and plane bodies participate in ONE implicit midpoint
    /// resistance solve. Plane mobility changes all shared film impulses; no
    /// reaction is applied as a later kick. For a mode `(n, a)` on a sphere,
    /// the opposite plane wrench uses `(n, (c_s-c_p) cross n + a)`. This retains
    /// world linear/angular momentum and vanishes under common rigid motion.
    /// The original roughness/cutoff laws and simultaneous mode budgets apply.
    ///
    /// `plane_external` contains center impulses followed by center couples,
    /// separate from film loads. Returned plane work enters finite-body kinetic
    /// energy; it is not prescribed-wall reservoir work. Viscous loss remains
    /// in `coupling.lubrication.dissipated_energy`, without a thermal destination.
    /// Fluid is collided and streamed once. All sphere velocities, plane
    /// velocities and populations publish together after admission completes.
    ///
    /// This is a frozen-geometry exchange. The caller owns subsequent pose,
    /// world-inertia and resolved-wall updates. A plane's extent and mass are
    /// not inferred from its infinite local geometric approximation. Resolved
    /// stationary-wall reaction remains a separate external port.
    ///
    /// # Errors
    /// Inherits native BGK/TRT, sphere, film and distinct-plane admission.
    /// Also refuses nonfinite plane loads, non-SPD inertia and inadmissible
    /// before/after plane velocities at coupled footpoints. There are at most
    /// 64 planes. Any failure or cancellation preserves all accepted state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_partial_saturation_finite_planes<const N: usize, const P: usize>(
        &mut self,
        fractions: [&[f64]; N],
        radii: [f64; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        planes: &mut [FiniteLubricationPlane3; P],
        plane_external: [[f64; 6]; P],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<LubricatedPlanesExchange3<N>> {
        checkpoint(gate)?;
        if P > 64 || !plane_external.into_iter().flatten().all(f64::is_finite) {
            return Err(invalid("invalid finite-plane count or external impulses"));
        }
        let native_planes: [LubricationPlane3; P] = core::array::from_fn(|p| LubricationPlane3 {
            recoil: Some(FinitePlaneInput {
                body: planes[p].body,
                external: plane_external[p],
                carrier_anchor: None,
            }),
            ..planes[p].geometry
        });
        let (coupling, (), exchanges) = self.step_with_partial_saturation_resistance_checked(
            fractions,
            bodies,
            external,
            origin,
            gate,
            Some((model, radii, &native_planes)),
            |_, _, _| Ok(()),
        )?;
        for (plane, exchange) in planes.iter_mut().zip(&exchanges) {
            let after = exchange
                .finite_body
                .expect("every admitted finite plane has a checked response")
                .body_after;
            plane.body = after;
            plane.geometry.velocity = after.velocity;
        }
        Ok(LubricatedPlanesExchange3 {
            coupling,
            planes: exchanges,
        })
    }
}

fn body_factor(body: PartialSaturationBody3) -> Result<fs_la::factor::Cholesky> {
    if !(body.mass.is_finite() && body.mass > 0.0)
        || !body
            .center
            .into_iter()
            .chain(body.velocity)
            .chain(body.angular_velocity)
            .chain(body.inertia_world.into_iter().flatten())
            .all(f64::is_finite)
        || !low_mach(body.velocity)
        || (0..3).any(|i| (0..i).any(|j| body.inertia_world[i][j] != body.inertia_world[j][i]))
    {
        return Err(invalid(
            "finite plane requires finite state and positive symmetric inertia",
        ));
    }
    Response::default().factor(body)
}

impl FreeResponse {
    pub(super) fn finite_plane(input: FinitePlaneInput) -> Result<Self> {
        let factor = body_factor(input.body)?;
        let mut response = Self {
            increment: input.external,
            ..Self::default()
        };
        factor.solve(&mut response.increment);
        let mut angular = [[0.0; 3]; 6];
        for axis in 0..6 {
            let mut rhs = [0.0; 6];
            rhs[axis] = 1.0;
            factor.solve(&mut rhs);
            for row in 0..6 {
                if axis < 3 {
                    response.linear_columns[row][axis] = rhs[row];
                } else {
                    angular[row][axis - 3] = rhs[row];
                }
            }
        }
        response.angular_columns = Some(angular);
        Ok(response)
    }
}

impl FinitePlaneResponse3 {
    pub(super) fn initial(body: PartialSaturationBody3) -> Self {
        Self {
            body_after: body,
            velocity_midpoint: body.velocity,
            angular_velocity_midpoint: body.angular_velocity,
            spin_impulse: [0.0; 3],
            kinetic_energy_change: 0.0,
            external_work: 0.0,
            momentum_residual: [0.0; 6],
            work_residual: 0.0,
        }
    }
}

pub(super) fn finish_response(
    input: FinitePlaneInput,
    free: &FreeResponse,
    exchange: &mut PlaneLubricationExchange3,
) -> Result<()> {
    let finite = exchange
        .finite_body
        .as_mut()
        .expect("finite plane response");
    let angular = free
        .angular_columns
        .as_ref()
        .expect("finite plane mobility");
    let increment = core::array::from_fn(|i| {
        free.increment[i]
            + dot(free.linear_columns[i], exchange.impulse)
            + dot(angular[i], finite.spin_impulse)
    });
    let after = Response::apply(input.body, increment)?;
    let dp =
        core::array::from_fn(|i| input.body.mass * (after.velocity[i] - input.body.velocity[i]));
    let dw = core::array::from_fn(|i| after.angular_velocity[i] - input.body.angular_velocity[i]);
    let dl = input.body.inertia_world.map(|row| dot(row, dw));
    finite.body_after = after;
    finite.velocity_midpoint =
        core::array::from_fn(|i| 0.5 * (input.body.velocity[i] + after.velocity[i]));
    finite.angular_velocity_midpoint = core::array::from_fn(|i| {
        0.5 * (input.body.angular_velocity[i] + after.angular_velocity[i])
    });
    finite.kinetic_energy_change =
        dot(dp, finite.velocity_midpoint) + dot(dl, finite.angular_velocity_midpoint);
    finite.external_work = dot(
        core::array::from_fn(|i| input.external[i]),
        finite.velocity_midpoint,
    ) + dot(
        core::array::from_fn(|i| input.external[i + 3]),
        finite.angular_velocity_midpoint,
    );
    finite.momentum_residual = core::array::from_fn(|i| {
        if i < 3 {
            dp[i] - input.external[i] - exchange.impulse[i]
        } else {
            dl[i - 3] - input.external[i] - finite.spin_impulse[i - 3]
        }
    });
    exchange.work = dot(exchange.impulse, finite.velocity_midpoint)
        + dot(finite.spin_impulse, finite.angular_velocity_midpoint);
    finite.work_residual = finite.kinetic_energy_change - finite.external_work - exchange.work;
    if !finite
        .spin_impulse
        .into_iter()
        .chain(finite.momentum_residual)
        .chain([
            finite.kinetic_energy_change,
            finite.external_work,
            finite.work_residual,
            exchange.work,
        ])
        .all(f64::is_finite)
    {
        return Err(invalid("finite-plane momentum/work account overflowed"));
    }
    Ok(())
}

pub(super) fn admit_point_speed(body: PartialSaturationBody3, point: [f64; 3]) -> Result<()> {
    let arm = core::array::from_fn(|i| point[i] - body.center[i]);
    if !low_mach(point_velocity(body.velocity, body.angular_velocity, arm)) {
        return Err(invalid(
            "finite-plane footpoint velocity exceeds low-Mach envelope",
        ));
    }
    Ok(())
}

// Shared identity is an explicit canonical first-face index, not floating
// equality of geometry or inertia. Validate it before any matrix assembly.
pub(super) fn carrier_owners(
    planes: &[LubricationPlane3],
    gate: &CancelGate,
) -> Result<[usize; 64]> {
    if planes.len() > 64 {
        return Err(invalid("lubrication supports at most 64 carrier faces"));
    }
    let mut owners = core::array::from_fn(|i| i);
    for (index, plane) in planes.iter().enumerate() {
        checkpoint(gate)?;
        let Some(input) = plane.recoil else { continue };
        let Some(anchor) = input.carrier_anchor else {
            continue;
        };
        if anchor > index {
            return Err(invalid("shared carrier must name its first face"));
        }
        let canonical = planes[anchor]
            .recoil
            .ok_or_else(|| invalid("shared carrier anchor must own finite inertia"))?;
        if canonical.carrier_anchor != Some(anchor)
            || canonical.body != input.body
            || canonical.external != input.external
        {
            return Err(invalid(
                "shared carrier faces disagree on inertia or external load",
            ));
        }
        owners[index] = anchor;
    }
    Ok(owners)
}

// Sum all face wrenches before solving ONE endpoint per carrier. Per-face
// spin impulses and work remain local; the private shared-surface wrapper
// projects the repeated endpoint metadata into a single public carrier receipt.
// No old public plane entry point can construct a carrier_anchor.
pub(super) fn finish_shared_responses(
    planes: &[LubricationPlane3],
    owners: &[usize; 64],
    free: &[FreeResponse],
    reports: &mut [PlaneLubricationExchange3],
    gate: &CancelGate,
) -> Result<()> {
    let mut linear = [[Dd::ZERO; 3]; 64];
    let mut spin = [[Dd::ZERO; 3]; 64];
    for (index, report) in reports.iter().enumerate() {
        checkpoint(gate)?;
        if let Some(finite) = report.finite_body {
            let owner = owners[index];
            for i in 0..3 {
                linear[owner][i] = linear[owner][i] + Dd::from_f64(report.impulse[i]);
                spin[owner][i] = spin[owner][i] + Dd::from_f64(finite.spin_impulse[i]);
            }
        }
    }
    let mut complete = [None; 64];
    for (index, plane) in planes.iter().enumerate() {
        checkpoint(gate)?;
        if let Some(input) = plane.recoil
            && owners[index] == index
        {
            let mut response = PlaneLubricationExchange3 {
                impulse: linear[index].map(Dd::to_f64),
                finite_body: Some(FinitePlaneResponse3 {
                    spin_impulse: spin[index].map(Dd::to_f64),
                    ..FinitePlaneResponse3::initial(input.body)
                }),
                ..PlaneLubricationExchange3::default()
            };
            finish_response(input, &free[index], &mut response)?;
            complete[index] = response.finite_body;
        }
    }
    let mut work = [Dd::ZERO; 64];
    for (index, report) in reports.iter_mut().enumerate() {
        checkpoint(gate)?;
        if let Some(local) = report.finite_body {
            let owner = owners[index];
            let full = complete[owner].expect("admitted shared carrier response");
            report.work = dot(report.impulse, full.velocity_midpoint)
                + dot(local.spin_impulse, full.angular_velocity_midpoint);
            if !report.work.is_finite() {
                return Err(invalid("shared carrier face work overflowed"));
            }
            work[owner] = work[owner] + Dd::from_f64(report.work);
            report.finite_body = Some(FinitePlaneResponse3 {
                spin_impulse: local.spin_impulse,
                ..full
            });
        }
    }
    for (index, report) in reports.iter_mut().enumerate() {
        checkpoint(gate)?;
        if let Some(finite) = &mut report.finite_body {
            finite.work_residual =
                finite.kinetic_energy_change - finite.external_work - work[owners[index]].to_f64();
            if !finite.work_residual.is_finite() {
                return Err(invalid("shared carrier work balance overflowed"));
            }
        }
    }
    Ok(())
}
