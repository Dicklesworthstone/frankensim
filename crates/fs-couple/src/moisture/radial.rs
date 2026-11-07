//! Source-integrated diffusion in cylinders, annuli, spheres and spherical shells.
//!
//! The radial conservation law is `rho_dry dM/dt = div(rho_dry grad K(M))`,
//! where the original material owner supplies `K'(M)=D(M)`. Each shell stores
//! its actual dry mass. A face uses the exact radial resistance `integral dr/A(r)`
//! between cell centers, so steady cylindrical logarithms and spherical inverse
//! radii are retained even on nonuniform meshes. Replacing this by a slab loses
//! both the curved-face flux and the radial storage measure.
//!
//! This is isothermal bound-water diffusion in a fixed, homogeneous dry skeleton.
//! The supplied property must describe isotropic or radial-principal transport;
//! arbitrary anisotropic tensors, axial/angular gradients, swelling feedback,
//! hysteresis remain outside this reduction. The atmosphere adapter supplies
//! finite exterior vapor inventories through the same original radial fluxes.

use super::*;

/// Finite vapor inventories coupled to native radial moisture diffusion.
pub mod atmosphere;
/// Source-partitioned diffusion through concentric coatings and material layers.
pub mod layered;
mod sdirk;
mod temperature;

/// Physical area and storage measure for a radial reduction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RadialMoistureSymmetry {
    /// Uniform along a cylinder's axis; both axial ends are sealed.
    Cylinder { length_m: f64 },
    /// Uniform in both angular coordinates.
    Sphere,
}

impl RadialMoistureSymmetry {
    fn factor(self) -> f64 {
        match self {
            Self::Cylinder { length_m } => 2.0 * core::f64::consts::PI * length_m,
            Self::Sphere => 4.0 * core::f64::consts::PI,
        }
    }

    fn area(self, radius: f64) -> f64 {
        match self {
            Self::Cylinder { .. } => self.factor() * radius,
            Self::Sphere => self.factor() * radius * radius,
        }
    }

    fn volume(self, a: f64, b: f64) -> f64 {
        // Factored differences retain thin shells without subtracting powers.
        match self {
            Self::Cylinder { .. } => self.factor() * (b - a) * (0.5 * a + 0.5 * b),
            Self::Sphere => self.factor() * (b - a) * (a * a + a * b + b * b) / 3.0,
        }
    }

    fn resistance(self, a: f64, b: f64) -> f64 {
        match self {
            Self::Cylinder { .. } => ((b - a) / a).ln_1p() / self.factor(),
            Self::Sphere => ((b - a) / b) / a / self.factor(),
        }
    }
}

/// Fixed radial shell mesh and source applicability; radii are ordered in meters.
#[derive(Clone, Debug)]
pub struct MoistureRadialGeometry {
    /// Radial symmetry, including cylinder length for physical inventories.
    pub symmetry: RadialMoistureSymmetry,
    /// N+1 strictly increasing faces for N shells; the first may be zero.
    pub faces_m: Vec<f64>,
    /// Admitted dry-basis moisture upper bound, kg water/kg dry solid.
    pub maximum_moisture: f64,
    /// Caller work bound, additionally limited by `MAX_MOISTURE_CELLS`.
    pub maximum_cells: usize,
}

/// A complete accepted radial moisture checkpoint. Preparation caches geometry;
/// each Picard iteration uses O(N) storage/work and the native tridiagonal LU.
#[derive(Clone, Debug)]
pub struct MoistureRadial {
    card: MaterialCard,
    reference: ResolvedMaterialStatePoint,
    transport: ResolvedMaterialStatePoint,
    property: String,
    geometry: MoistureRadialGeometry,
    centers_m: Vec<f64>,
    dry_masses_kg: Vec<f64>,
    mass_scale_kg: f64,
    // rho_dry / (mass_scale * integral dr/A), excluding material D and time.
    interior_metric: Vec<f64>,
    boundary_metric: [f64; 2],
    // rho_dry * A / mass_scale, excluding source surface speed and time.
    surface_metric: [f64; 2],
    moisture: Vec<f64>,
    time_s: f64,
    identity: ContentHash,
}

impl MoistureRadial {
    /// Bind one original material and a possibly nonuniform radial mesh.
    /// Zero diffusivity, including a dry D(0)=0 front, is admitted through the
    /// original Kirchhoff source integral. Density is fixed at the reference
    /// state and converted to dry-skeleton density on its explicit moisture basis.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        card: MaterialCard,
        reference_point: QueryPoint,
        diffusivity_property: &str,
        selection: MaterialPropertySelection,
        geometry: MoistureRadialGeometry,
        initial_moisture: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let n = initial_moisture.len();
        if n == 0
            || n > geometry.maximum_cells
            || n > MAX_MOISTURE_CELLS
            || geometry.faces_m.len() != n + 1
            || !positive(geometry.symmetry.factor())
            || !positive(geometry.maximum_moisture)
            || geometry.faces_m.iter().any(|r| !r.is_finite() || *r < 0.0)
            || geometry.faces_m.windows(2).any(|r| r[0] >= r[1])
            || initial_moisture
                .iter()
                .any(|m| !m.is_finite() || *m < 0.0 || *m > geometry.maximum_moisture)
        {
            return Err(MoistureError::Invalid(
                "bounded radial mesh, positive symmetry measure and supported moisture required",
            ));
        }
        let (reference, transport, density) = binding::bind_source(
            &card,
            &reference_point,
            &[diffusivity_property],
            selection,
            geometry.maximum_moisture,
            true,
        )?;
        let centers_m: Vec<_> = geometry
            .faces_m
            .windows(2)
            .map(|r| r[0] + 0.5 * (r[1] - r[0]))
            .collect();
        let dry_masses_kg: Vec<_> = geometry
            .faces_m
            .windows(2)
            .map(|r| density * geometry.symmetry.volume(r[0], r[1]))
            .collect();
        let mass_scale_kg = dry_masses_kg.iter().copied().fold(0.0, f64::max);
        let metric = |a, b| density / geometry.symmetry.resistance(a, b) / mass_scale_kg;
        let interior_metric: Vec<_> = centers_m.windows(2).map(|r| metric(r[0], r[1])).collect();
        let inner = geometry.faces_m[0];
        let outer = geometry.faces_m[n];
        let boundary_metric = [
            if inner == 0.0 {
                0.0
            } else {
                metric(inner, centers_m[0])
            },
            metric(centers_m[n - 1], outer),
        ];
        let surface_metric = [
            density * geometry.symmetry.area(inner) / mass_scale_kg,
            density * geometry.symmetry.area(outer) / mass_scale_kg,
        ];
        if dry_masses_kg
            .iter()
            .any(|m| !positive(*m) || !positive(*m / mass_scale_kg))
            || !positive(dry_masses_kg.iter().sum())
            || interior_metric.iter().any(|m| !positive(*m))
            || !positive(boundary_metric[1])
            || !positive(surface_metric[1])
            || (inner > 0.0 && (!positive(boundary_metric[0]) || !positive(surface_metric[0])))
            || centers_m
                .iter()
                .zip(geometry.faces_m.windows(2))
                .any(|(r, f)| *r <= f[0] || *r >= f[1])
        {
            return Err(MoistureError::Invalid(
                "radial storage or geometric conductance is unrepresentable",
            ));
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-couple.moisture-radial.v1");
        identity.update(reference.identity().as_bytes());
        match geometry.symmetry {
            RadialMoistureSymmetry::Cylinder { length_m } => {
                identity.update(&[0]);
                identity.update(&length_m.to_bits().to_le_bytes());
            }
            RadialMoistureSymmetry::Sphere => identity.update(&[1]),
        }
        identity.update(&geometry.maximum_moisture.to_bits().to_le_bytes());
        for radius in &geometry.faces_m {
            identity.update(&radius.to_bits().to_le_bytes());
        }
        let result = Self {
            card,
            reference,
            transport,
            property: diffusivity_property.into(),
            geometry,
            centers_m,
            dry_masses_kg,
            mass_scale_kg,
            interior_metric,
            boundary_metric,
            surface_metric,
            moisture: initial_moisture,
            time_s: 0.0,
            identity: identity.finalize(),
        };
        if !result.water_mass_kg().is_finite() || !positive(result.wet_mass_kg()) {
            return Err(MoistureError::Invalid(
                "radial water inventory is unrepresentable",
            ));
        }
        result.transport_states(&result.moisture, cx)?;
        Ok(result)
    }

    /// Accepted shell contents, kg water/kg dry solid, ordered radially outward.
    pub fn moisture(&self) -> &[f64] {
        &self.moisture
    }
    /// Arithmetic shell centers where the moisture field is represented.
    pub fn centers_m(&self) -> &[f64] {
        &self.centers_m
    }
    /// Fixed dry-skeleton inventory of each physical shell.
    pub fn dry_masses_kg(&self) -> &[f64] {
        &self.dry_masses_kg
    }
    /// Immutable admitted mesh and physical symmetry.
    pub const fn geometry(&self) -> &MoistureRadialGeometry {
        &self.geometry
    }
    /// Elapsed time at this accepted checkpoint.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Original source/reference and geometry identity, independent of history.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Total dry-skeleton inventory.
    pub fn dry_mass_kg(&self) -> f64 {
        self.dry_masses_kg.iter().sum()
    }
    /// Actual volume-weighted water inventory.
    pub fn water_mass_kg(&self) -> f64 {
        self.water(&self.moisture)
    }
    /// Total dry-skeleton and bound-water mass.
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_kg() + self.water_mass_kg()
    }
    /// Original pinned source state that supplies the dry reference density.
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }

    pub(super) const fn material_card(&self) -> &MaterialCard {
        &self.card
    }

    /// Originally pinned diffusivity source in the current prescribed environment.
    /// The dry-density and stress-free reference remain available in `reference`.
    pub const fn transport_context(&self) -> &ResolvedMaterialStatePoint {
        &self.transport
    }

    pub(super) fn field(&self) -> binding::MoistureField<'_> {
        binding::MoistureField {
            card: &self.card,
            reference: &self.reference,
            transport: &self.transport,
            moisture: &self.moisture,
            dry_masses: binding::DryMasses::Cells(&self.dry_masses_kg),
            time_s: self.time_s,
            identity: self.identity,
        }
    }

    fn water(&self, moisture: &[f64]) -> f64 {
        moisture
            .iter()
            .zip(&self.dry_masses_kg)
            .map(|(m, mass)| m * mass)
            .sum()
    }

    fn transport_states(
        &self,
        moisture: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<ResolvedMaterialStatePoint>, MoistureError> {
        let mut result = Vec::with_capacity(moisture.len());
        for (i, &m) in moisture.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            let point = self
                .transport
                .query()
                .with("moisture_content", m)
                .map_err(|_| MoistureError::Invalid("invalid radial moisture coordinate"))?;
            result.push(self.transport.at_query_point(&self.card, &point)?);
        }
        Ok(result)
    }

    /// Conservative backward Euler with the complete endpoint-dependent source
    /// law. Faces are inner/outer; a solid center has exactly zero area and must
    /// be sealed. Finite source surface exchange uses its real curved area.
    /// Every rejected solve leaves this accepted checkpoint unchanged.
    pub fn advance(
        &self,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureStepReport), MoistureError> {
        poll(cx)?;
        let temperature = self.temperature_kelvin();
        let boundaries = [
            faces[0].resolve(temperature)?,
            faces[1].resolve(temperature)?,
        ];
        self.advance_stage(dt_s, dt_s, boundaries, step, nonlinear, None, cx)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::moisture) fn advance_stage(
        &self,
        dt_s: f64,
        implicit_dt_s: f64,
        boundaries: [ResolvedMoistureFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        known: Option<&crate::moisture::sdirk::StageLoad>,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureStepReport), MoistureError> {
        poll(cx)?;
        nonlinear.validate()?;
        let end = self.time_s + dt_s;
        if !positive(dt_s)
            || !positive(implicit_dt_s)
            || !end.is_finite()
            || end <= self.time_s
            || !positive(step.water_balance_tolerance_kg)
            || !positive(step.maximum_content_increment)
            || (self.geometry.faces_m[0] == 0.0 && boundaries[0].equilibrium.is_some())
        {
            return Err(MoistureError::Invalid(
                "forward radial interval, positive limits and sealed zero-area center required",
            ));
        }
        let reservoir = [boundaries[0].equilibrium, boundaries[1].equilibrium];
        if reservoir
            .into_iter()
            .flatten()
            .any(|m| m > self.geometry.maximum_moisture)
        {
            return Err(MoistureError::Invalid(
                "radial surface exceeds moisture support",
            ));
        }
        let lower = diffusion_reference(&self.moisture, &boundaries);
        let upper = self
            .moisture
            .iter()
            .copied()
            .chain(
                boundaries
                    .iter()
                    .filter(|b| b.coefficient_m_s != Some(0.0))
                    .filter_map(|b| b.equilibrium),
            )
            .fold(0.0, f64::max);
        // A finite reservoir can move between RK stages. The actual retained
        // first-stage load may produce content outside the old/current-probe
        // envelope, while still inside the physical admitted stage range.
        // Extend only source support; preserve the algebraic reference shift.
        let source_bounds = known.map_or([lower, upper], |load| {
            let [a, b] = load.source_bounds();
            [lower.min(a), upper.max(b)]
        });
        let profile = KirchhoffMoistureDiffusivity::from_binding(
            &self.card,
            &self.transport,
            &self.property,
            source_bounds,
        )?;
        let start = self.transport_states(&self.moisture, cx)?;
        let mut iterate = self.moisture.clone();
        let mut workspace = fs_la::tridiagonal::TridiagonalWorkspace::new(iterate.len())
            .map_err(MoistureError::Tridiagonal)?;
        let mut system = self.system(implicit_dt_s, &boundaries, &profile, &iterate, lower, cx)?;
        if let Some(known) = known {
            known.add_to(&mut system.1, cx)?;
        }
        let mut solves = 0;
        for iteration in 1..=nonlinear.maximum_iterations {
            poll(cx)?;
            let solution =
                system
                    .0
                    .solve_tridiagonal(&system.1, step.linear, &mut workspace, cx)?;
            solves += solution.report.iters;
            if !solution.report.converged_euclidean() {
                return Err(MoistureError::NotConverged(solution.report));
            }
            let candidate: Vec<_> = solution.x.iter().map(|m| m + lower).collect();
            if let Some(known) = known {
                known.check_range_with(
                    &candidate,
                    self.geometry.maximum_moisture,
                    system
                        .2
                        .iter()
                        .enumerate()
                        .filter(|(_, coupling)| coupling.0 != 0.0)
                        .map(|(side, _)| reservoir[side].expect("active radial face")),
                    cx,
                )?;
            }
            let stored = [self.water_mass_kg(), self.water(&candidate)];
            if !stored[1].is_finite() || !positive(self.dry_mass_kg() + stored[1]) {
                return Err(MoistureError::Invalid(
                    "radial candidate water or wet inventory is unrepresentable",
                ));
            }
            let mut endpoint =
                self.system(implicit_dt_s, &boundaries, &profile, &candidate, lower, cx)?;
            if let Some(known) = known {
                known.add_to(&mut endpoint.1, cx)?;
            }
            let (_, residual) = endpoint.0.defect(&solution.x, &endpoint.1, cx)?;
            let mut boundary_water_kg =
                self.boundary_transfers(&candidate, &boundaries, &endpoint.2)?;
            if let Some(known) = known {
                for (total, prior) in boundary_water_kg.iter_mut().zip(known.boundary_water_kg) {
                    *total += prior;
                }
            }
            if boundary_water_kg.iter().any(|water| !water.is_finite()) {
                return Err(MoistureError::Invalid(
                    "radial boundary water inventory is unrepresentable",
                ));
            }
            let balance = stored[1] - stored[0] - boundary_water_kg[0] - boundary_water_kg[1];
            let report = MoistureNonlinearReport {
                iterations: iteration,
                tridiagonal_solves: solves,
                relative_residual: residual,
                water_balance_residual_kg: balance,
            };
            if residual < nonlinear.relative_tolerance
                && balance.is_finite()
                && balance.abs() <= step.water_balance_tolerance_kg
            {
                if let Some(known) = known {
                    known.check_range_with(
                        &candidate,
                        self.geometry.maximum_moisture,
                        endpoint
                            .2
                            .iter()
                            .enumerate()
                            .filter(|(_, coupling)| coupling.0 != 0.0)
                            .map(|(side, _)| reservoir[side].expect("active radial face")),
                        cx,
                    )?;
                }
                for (i, (&old, &m)) in self.moisture.iter().zip(&candidate).enumerate() {
                    if !m.is_finite()
                        || m < 0.0
                        || m > self.geometry.maximum_moisture
                        || (known.is_none() && (m < lower || m > upper))
                    {
                        return Err(MoistureError::Invalid(
                            "radial candidate violates moisture maximum principle",
                        ));
                    }
                    if (m - old).abs() > step.maximum_content_increment {
                        return Err(MoistureError::ContentIncrement {
                            cell: i,
                            increment: (m - old).abs(),
                            maximum: step.maximum_content_increment,
                        });
                    }
                }
                let transport_end = self.transport_states(&candidate, cx)?;
                let surface_moisture = std::array::from_fn(|side| {
                    reservoir[side].map(|m| {
                        let i = if side == 0 { 0 } else { candidate.len() - 1 };
                        endpoint.2[side].1.mul_add(m - candidate[i], candidate[i])
                    })
                });
                let [left, right] = boundaries;
                let report = MoistureStepReport {
                    interval_s: [self.time_s, end],
                    boundary_water_kg,
                    surface_moisture,
                    reservoir_moisture: reservoir,
                    surface_sorption: [left.sorption, right.sorption],
                    surface_transfer: [left.transfer, right.transfer],
                    stored_water_kg: stored,
                    water_balance_residual_kg: balance,
                    transport_start: start,
                    transport_end,
                    kirchhoff_source: Some(profile.identity()),
                    linear: solution.report,
                    nonlinear: Some(report),
                };
                let mut next = self.clone();
                next.moisture = candidate;
                next.time_s = end;
                poll(cx)?;
                return Ok((next, report));
            }
            if iteration == nonlinear.maximum_iterations {
                return Err(MoistureError::NonlinearNotConverged(report));
            }
            if nonlinear.relaxation == 1.0 {
                system = endpoint;
            } else {
                for (old, new) in iterate.iter_mut().zip(candidate) {
                    *old = nonlinear.relaxation.mul_add(new - *old, *old);
                }
                system = self.system(implicit_dt_s, &boundaries, &profile, &iterate, lower, cx)?;
                if let Some(known) = known {
                    known.add_to(&mut system.1, cx)?;
                }
            }
        }
        unreachable!("positive nonlinear iteration budget")
    }

    fn boundary_transfers(
        &self,
        moisture: &[f64],
        faces: &[ResolvedMoistureFace; 2],
        coupling: &[(f64, f64); 2],
    ) -> Result<[f64; 2], MoistureError> {
        let water = std::array::from_fn(|side| {
            let i = if side == 0 { 0 } else { moisture.len() - 1 };
            let difference = faces[side].equilibrium.unwrap_or(moisture[i]) - moisture[i];
            // Preserve the original range-safe physical transfer, shared by
            // complete backward-Euler steps and the weighted SDIRK face load.
            let mut factors = [self.mass_scale_kg, coupling[side].0, difference.abs()];
            factors.sort_unstable_by(f64::total_cmp);
            ((factors[0] * factors[2]) * factors[1]).copysign(difference)
        });
        if water.iter().any(|value| !value.is_finite()) {
            return Err(MoistureError::Invalid(
                "radial boundary water inventory is unrepresentable",
            ));
        }
        Ok(water)
    }

    fn system(
        &self,
        dt: f64,
        faces: &[ResolvedMoistureFace; 2],
        profile: &KirchhoffMoistureDiffusivity,
        moisture: &[f64],
        reference: f64,
        cx: &Cx<'_>,
    ) -> Result<(DiffusionStep, Vec<f64>, [(f64, f64); 2]), MoistureError> {
        let mass: Vec<_> = self
            .dry_masses_kg
            .iter()
            .map(|m| m / self.mass_scale_kg)
            .collect();
        let mut operator = DiffusionStep {
            diagonal: mass.clone(),
            storage_boundary: mass.clone(),
            off: Vec::with_capacity(mass.len() - 1),
        };
        let mut rhs: Vec<_> = mass
            .iter()
            .zip(&self.moisture)
            .map(|(mass, m)| mass * (m - reference))
            .collect();
        for (i, &metric) in self.interior_metric.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            let a = kirchhoff::scaled_diffusivity(
                dt * metric,
                profile.mean_diffusivity(moisture[i], moisture[i + 1])?,
            )?;
            operator.off.push(a);
            operator.diagonal[i] += a;
            operator.diagonal[i + 1] += a;
        }
        let mut coupling = [(0.0, 0.0); 2];
        for (side, i) in [0, mass.len() - 1].into_iter().enumerate() {
            if let Some(equilibrium) = faces[side].equilibrium {
                coupling[side] = self.surface(side, dt, &faces[side], profile, moisture[i], cx)?;
                operator.diagonal[i] += coupling[side].0;
                operator.storage_boundary[i] += coupling[side].0;
                rhs[i] += coupling[side].0 * (equilibrium - reference);
            }
        }
        if operator.diagonal.iter().any(|a| !positive(*a)) || rhs.iter().any(|b| !b.is_finite()) {
            return Err(MoistureError::Invalid("radial diffusion system overflow"));
        }
        Ok((operator, rhs, coupling))
    }

    fn surface(
        &self,
        side: usize,
        dt: f64,
        face: &ResolvedMoistureFace,
        profile: &KirchhoffMoistureDiffusivity,
        cell: f64,
        cx: &Cx<'_>,
    ) -> Result<(f64, f64), MoistureError> {
        let equilibrium = face.equilibrium.expect("active radial surface");
        if face.coefficient_m_s == Some(0.0) {
            return Ok((0.0, 0.0));
        }
        let evaluate = |fraction: f64| -> Result<(f64, f64), MoistureError> {
            let surface = fraction.mul_add(equilibrium - cell, cell);
            let mean = profile.mean_diffusivity(cell, surface)?;
            let Some(h) = face.coefficient_m_s else {
                return Ok((
                    kirchhoff::scaled_diffusivity(dt * self.boundary_metric[side], mean)?,
                    1.0,
                ));
            };
            if mean == 0.0 {
                return Ok((0.0, 1.0));
            }
            // Combine transfer speeds before scaling by time and area. A very
            // large h must approach surface equilibrium without overflowing an
            // exterior conductance that is never the limiting resistance.
            let speed = self.boundary_metric[side] / self.surface_metric[side] * mean;
            if !positive(speed) {
                return Err(MoistureError::Invalid(
                    "radial half-shell transfer speed is unrepresentable",
                ));
            }
            let fraction = if h <= speed {
                let ratio = h / speed;
                ratio / (1.0 + ratio)
            } else {
                1.0 / (1.0 + speed / h)
            };
            let effective = h.min(speed) / (1.0 + h.min(speed) / h.max(speed));
            if !positive(fraction) || !positive(effective) {
                return Err(MoistureError::Invalid("radial surface transfer underflow"));
            }
            let conductance =
                kirchhoff::scaled_diffusivity(dt * self.surface_metric[side], effective)?;
            Ok((conductance, fraction))
        };
        if face.coefficient_m_s.is_none() || cell == equilibrium {
            return evaluate(1.0);
        }
        let mut bracket = [0.0, 1.0];
        let mut iterations = 0;
        for iteration in 1..=128 {
            iterations = iteration;
            poll(cx)?;
            let fraction = 0.5 * (bracket[0] + bracket[1]);
            let coupling = evaluate(fraction)?;
            let checked = evaluate(coupling.1)?;
            let close = |a: f64, b: f64| {
                a == b || (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
            };
            if close(coupling.0, checked.0) && close(coupling.1, checked.1) {
                return Ok(checked);
            }
            if fraction == bracket[0] || fraction == bracket[1] {
                break;
            }
            if fraction > coupling.1 {
                bracket[1] = fraction;
            } else {
                bracket[0] = fraction;
            }
        }
        Err(MoistureError::KirchhoffSurfaceNotConverged {
            iterations,
            fraction_bracket: bracket,
        })
    }
}
