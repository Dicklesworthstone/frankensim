//! Existing missing-film resistance inside the simultaneous physical-dry solve.

use super::*;

/// Accepted dry spheres/fluid and the prescribed surfaces' separate reactions.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct DryLubricatedPlanesExchange3<const N: usize> {
    /// Physical dry-body correction and native simultaneous film exchange.
    /// Pair and aggregate plane impulses/work remain in
    /// `motion.native.coupling.lubrication`, separate from resolved fluid loads.
    pub motion: DryBodiesExchange3<N>,
    /// Film reactions RECEIVED by the prescribed surfaces, in input order.
    /// These planes have no finite-inertia recoil. Their driving reservoirs
    /// supply the opposite impulse, angular impulse and work to this system.
    pub planes: Vec<PlaneLubricationExchange3>,
}

impl BoundaryGrid3 {
    /// Advance 1 through 8 physical dry spheres with simultaneous missing films.
    ///
    /// The existing sphere-pair normal and optional reciprocal tangent/spin
    /// resistance are evaluated INSIDE every full coupled dry-body candidate.
    /// Each original positive dry mass and complete inertia tensor is retained;
    /// there is no excess-density mass or sequential post-step lubrication kick.
    /// The full 6N retained-fluid correction sees the film-modified velocities
    /// and the same single native collision/streaming result for all bodies.
    ///
    /// Actual supplied radii define the sampled spherical geometry. The model's
    /// independent roughness/cutoff scales, positive gaps, film/mode caps and
    /// viscosity admission are unchanged. Interior Guo momentum and midpoint
    /// fluid-source subtraction retain the ordinary dry-body convention.
    /// Film impulses/torques/work and nonnegative modeled loss belong solely to
    /// `native.coupling.lubrication`; they are not applied to populations again.
    /// A singleton and zero active resistance retain the original dry arithmetic.
    ///
    /// # Errors
    /// Inherits both the bounded dry-assembly and native sphere-film admission.
    /// Any geometry, solve, budget or cancellation refusal leaves all body and
    /// fluid state intact. Contact, a GCL, exact no-slip, unconditional added-mass
    /// stability and conversion of film loss into heat remain outside this API.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_lubricated_spheres<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryBodiesExchange3<N>> {
        self.step_with_dry_lubricated_spheres_and_planes(
            samplers,
            radii,
            bodies,
            external,
            model,
            &[],
            origin,
            config,
            gate,
        )
        .map(|report| report.motion)
    }

    /// The same coupled dry-sphere/film solve with endpoint-sampled fluid
    /// forcing. Each body's fictitious-interior source uses the shared force
    /// trapezoid; pair-film exchange stays internal and is counted once.
    ///
    /// # Errors
    /// Inherits the dry sphere/film and endpoint-force admission. Any refusal
    /// preserves all bodies, fluid populations and the accepted source scale.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_lubricated_spheres_and_endpoint_force_scale<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        next_scale: f64,
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryBodiesExchange3<N>> {
        self.step_with_dry_lubricated_spheres_and_planes_and_endpoint_force_scale(
            samplers, radii, bodies, external, model, &[], next_scale, origin, config, gate,
        )
        .map(|report| report.motion)
    }

    /// Advance physical dry spheres with shared sphere-pair/prescribed-plane films.
    ///
    /// The native resistance owner retains all shared-body and tangent/spin
    /// coupling. Plane points are supplied at the tick START; their constant
    /// prescribed translations determine midpoint gaps and both admitted affine
    /// half-drifts. Planes provide only missing-film resistance and its external
    /// reaction/work: they do not create or move resolved fluid boundaries, and
    /// they have no finite mass or rotational recoil. The caller supplies any
    /// matching resolved-wall model and advances prescribed plane positions
    /// between accepted ticks with `LubricationPlane3::translated`.
    ///
    /// Individual dry-body momentum includes its film wrench once. The combined
    /// exterior-fluid-plus-dry momentum account includes the opposite received
    /// plane reaction about the SAME fixed world origin, including each reaction
    /// couple. Internal sphere-pair wrenches cancel from that common account.
    /// Authored dry loads, fictitious-fluid corrections, imposed fluid sources,
    /// resolved wall exchange and prescribed film ports remain separate.
    ///
    /// The original positive film law owns body/plane midpoint work and modeled
    /// viscous loss. A driven plane may inject body work; its reservoir supplies
    /// the negative received plane work. No thermodynamic heat destination or
    /// complete finite-gap/contact law is inferred. Empty planes reproduce the
    /// sphere-only entry point exactly; a singleton WITHOUT planes delegates to
    /// the original single-dry owner.
    ///
    /// # Errors
    /// Inherits the dry-assembly candidate cap, conditioning and volume checks,
    /// plus the existing distinct-plane, positive-gap, film/mode caps and tangent
    /// asymptotic limits. All final geometry, population and physical momentum
    /// accounts complete before the sole body/fluid publication. Refused probes
    /// and cancellation change only reusable scratch; all plane inputs are values.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_lubricated_spheres_and_planes<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        planes: &[LubricationPlane3],
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedPlanesExchange3<N>> {
        self.dry_lubricated_spheres_and_planes_transition(
            samplers, radii, bodies, external, model, planes, origin, config, gate, None,
        )
    }

    /// Endpoint fluid-force history through the original coupled sphere-pair
    /// and prescribed-plane resistance owner. Plane reaction/work remain the
    /// separate external film port; no source is applied twice to populations.
    ///
    /// # Errors
    /// Inherits all held-force dry/film/plane checks and endpoint-scale admission.
    /// Body, fluid and source state publish together after the final accounts.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_lubricated_spheres_and_planes_and_endpoint_force_scale<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        planes: &[LubricationPlane3],
        next_scale: f64,
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedPlanesExchange3<N>> {
        self.validate_immersed_endpoint_force_scale(next_scale, gate)?;
        self.dry_lubricated_spheres_and_planes_transition(
            samplers, radii, bodies, external, model, planes, origin, config, gate,
            Some(next_scale),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn dry_lubricated_spheres_and_planes_transition<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        planes: &[LubricationPlane3],
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
        endpoint_force_scale: Option<f64>,
    ) -> Result<DryLubricatedPlanesExchange3<N>> {
        model.admit(radii)?;
        // Validate even the singleton fast path before it can publish. Plane
        // geometry/rate invariants are already enforced by its constructor.
        crate::d3q19::boundary::partial_saturation::lubrication::planes::admit_planes(planes)?;
        let spheres = radii
            .map(|r| BodyFrameEllipsoid3::sphere(r).expect("lubrication admitted every radius"));
        self.step_with_dry_partial_saturation_resistance_with_endpoint_force(
            samplers,
            core::array::from_fn(|a| &spheres[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            config,
            gate,
            Some((model, planes)),
            None,
            endpoint_force_scale,
        )
        .map(|(report, _)| report)
    }
}
