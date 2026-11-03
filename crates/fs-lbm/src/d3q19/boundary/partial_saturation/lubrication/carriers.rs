//! Several film faces carried by one finite spatial inertia.

use super::*;

/// Geometry and missing-resistance scales, without an independent rigid state.
///
/// Points and normals use world axes in the frozen entry point and body axes
/// when attached by a moving owner. A patch does not infer a finite face extent,
/// resolved wall, carrier shape, displaced volume or thermal destination.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LubricationPlanePatch3 {
    geometry: LubricationPlane3,
}

impl LubricationPlanePatch3 {
    /// Admit a point, nonzero normal and original positive film scales.
    ///
    /// # Errors
    /// Inherits the original plane geometry and normal-gap admission.
    pub fn new(
        point: [f64; 3],
        normal: [f64; 3],
        cutoff_gap: f64,
        roughness_gap: f64,
    ) -> Result<Self> {
        Ok(Self {
            geometry: LubricationPlane3::new(point, normal, [0.0; 3], cutoff_gap, roughness_gap)?,
        })
    }

    /// Enable the original reciprocal tangent/parallel-spin missing resistance.
    ///
    /// # Errors
    /// Inherits the original tangent-gap scale admission.
    pub fn with_tangential(mut self, cutoff_gap: f64, roughness_gap: f64) -> Result<Self> {
        self.geometry = self.geometry.with_tangential(cutoff_gap, roughness_gap)?;
        Ok(self)
    }

    /// Point in the coordinate frame declared by the enclosing owner.
    pub fn point(self) -> [f64; 3] {
        self.geometry.point()
    }
    /// Unit normal into the spheres' permitted half-space.
    pub fn normal(self) -> [f64; 3] {
        self.geometry.normal()
    }
    /// Original normal-film cutoff in lattice lengths.
    pub fn cutoff_gap(self) -> f64 {
        self.geometry.cutoff_gap
    }
    /// Original positive normal-gap regularization in lattice lengths.
    pub fn roughness_gap(self) -> f64 {
        self.geometry.roughness_gap
    }
    /// Optional tangent/parallel-spin cutoff and roughness in lattice lengths.
    pub fn tangential_gaps(self) -> Option<(f64, f64)> {
        self.geometry.tangential
    }

    pub(in crate::d3q19::boundary) fn from_plane(plane: LubricationPlane3) -> Self {
        Self {
            geometry: LubricationPlane3 {
                velocity: [0.0; 3],
                axial_spin: 0.0,
                recoil: None,
                ..plane
            },
        }
    }

    pub(in crate::d3q19::boundary) fn at(self, point: [f64; 3], normal: [f64; 3]) -> Result<Self> {
        let mut patch = Self::new(point, normal, self.cutoff_gap(), self.roughness_gap())?;
        if let Some((cutoff, roughness)) = self.tangential_gaps() {
            patch = patch.with_tangential(cutoff, roughness)?;
        }
        Ok(patch)
    }

    pub(in crate::d3q19::boundary) fn carried(self, carrier: usize) -> SharedLubricationSurface3 {
        SharedLubricationSurface3::Carried(CarrierLubricationFace3 {
            carrier,
            geometry: self.geometry,
        })
    }
}

/// A frozen planar film face referring to one entry of a carrier array.
/// It owns geometry and film scales, but no additional mass or inertia.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierLubricationFace3 {
    carrier: usize,
    geometry: LubricationPlane3,
}

impl CarrierLubricationFace3 {
    /// Bind world-frame geometry to a carrier index. The supplied carrier owns
    /// translation and spin; this face does not prescribe an additional rate.
    ///
    /// # Errors
    /// Refuses invalid plane geometry or film scales. Carrier-index admission
    /// occurs when the complete surface/carrier assembly is supplied.
    pub fn new(
        carrier: usize,
        point: [f64; 3],
        normal: [f64; 3],
        cutoff_gap: f64,
        roughness_gap: f64,
    ) -> Result<Self> {
        Ok(Self {
            carrier,
            geometry: LubricationPlane3::new(point, normal, [0.0; 3], cutoff_gap, roughness_gap)?,
        })
    }

    /// Enable the original reciprocal tangent/parallel-spin film correction.
    ///
    /// # Errors
    /// Inherits the original positive-scale and small-gap admission.
    pub fn with_tangential(mut self, cutoff_gap: f64, roughness_gap: f64) -> Result<Self> {
        self.geometry = self.geometry.with_tangential(cutoff_gap, roughness_gap)?;
        Ok(self)
    }

    /// Stable index of the one mechanical carrier.
    pub const fn carrier(self) -> usize {
        self.carrier
    }

    /// Point on the frozen surface in world lattice coordinates.
    pub fn point(self) -> [f64; 3] {
        self.geometry.point()
    }

    /// Unit normal into the spheres' permitted half-space.
    pub fn normal(self) -> [f64; 3] {
        self.geometry.normal()
    }
}

/// A prescribed film surface or a face of a shared finite carrier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SharedLubricationSurface3 {
    /// Original prescribed translation/axial spin and its actuator-work port.
    Prescribed(LubricationPlane3),
    /// Geometry attached to a separately supplied finite carrier.
    Carried(CarrierLubricationFace3),
}

/// One face's film reaction, without a duplicated carrier energy inventory.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierFaceExchange3 {
    /// Carrier index, or `None` for a prescribed actuator surface.
    pub carrier: Option<usize>,
    /// Received linear film impulse.
    pub impulse: [f64; 3],
    /// Received angular impulse about the requested world origin.
    pub angular_impulse: [f64; 3],
    /// Free reaction couples at this face's film footpoints.
    pub couple_impulse: [f64; 3],
    /// Work received at the actual face/carrier midpoint motion.
    pub work: f64,
    /// Number of active sphere films on this face.
    pub active_spheres: usize,
    /// Number of active tangent/parallel-spin film blocks on this face.
    pub tangential_spheres: usize,
}

/// One mechanical carrier's aggregate response, counted exactly once.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierLubricationExchange3 {
    /// Sum of the carrier's received face impulses.
    pub impulse: [f64; 3],
    /// Sum of received face moments about the requested world origin.
    pub angular_impulse: [f64; 3],
    /// Sum of face work, entering this carrier's kinetic state.
    pub work: f64,
    /// Complete carrier endpoint, center-spin impulse, kinetic/external work
    /// and measured momentum/work residuals. Spin impulse sums all its faces.
    pub motion: FinitePlaneResponse3,
    /// Active sphere/face films, counting distinct faces separately.
    pub active_films: usize,
    /// Active tangent/parallel-spin blocks, counting faces separately.
    pub tangential_films: usize,
}

/// One native fluid/body tick, face reactions and unique carrier inventories.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct SharedCarrierLubricationExchange3<const N: usize> {
    /// Original single simultaneous immersed-body/fluid response.
    pub coupling: PartialSaturationBodiesExchange3<N>,
    /// Face-only reactions in the original surface order.
    pub surfaces: Vec<CarrierFaceExchange3>,
    /// Unique finite responses in the original carrier order.
    pub carriers: Vec<CarrierLubricationExchange3>,
}

impl BoundaryGrid3 {
    /// Couple spheres, fluid, prescribed surfaces and multi-face carriers once.
    ///
    /// Every face sharing a carrier uses its ONE mass, full world inertia,
    /// center velocity and spin. Cross-face mobility enters the original
    /// simultaneous film matrix. External center impulses/couples are supplied
    /// once per carrier, independent of its number of faces.
    ///
    /// Face reactions and carrier inventories have separate result types:
    /// summing faces cannot silently count a carrier's kinetic energy twice.
    /// Prescribed face work belongs to its actuator; carried face work enters
    /// the corresponding finite carrier. Film dissipation retains the original
    /// missing-film account and has no inferred thermal destination.
    ///
    /// The complete result is prepared before the original fluid/sphere
    /// publication. Carrier assignment afterward is infallible. Every carrier
    /// must own at least one face; empty surfaces require empty carriers.
    /// Original positive-gap, low-Mach, film-count and mode budgets apply.
    ///
    /// Geometry, carrier centers and world inertia are frozen for this tick.
    /// These local infinite planes do not create resolved walls, face extents,
    /// pose trajectories, gyroscopic evolution or solid contact.
    ///
    /// # Errors
    /// Refuses more than 64 surfaces/carriers, invalid carrier references,
    /// unused carriers, mismatched/nonfinite external loads and all original
    /// geometry, inertia, response, allocation and cancellation failures.
    /// Failure preserves fluid populations, spheres and every carrier.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_partial_saturation_shared_carriers<const N: usize>(
        &mut self,
        fractions: [&[f64]; N],
        radii: [f64; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        surfaces: &[SharedLubricationSurface3],
        carriers: &mut [PartialSaturationBody3],
        carrier_external: &[[f64; 6]],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<SharedCarrierLubricationExchange3<N>> {
        let native = prepare(surfaces, carriers, carrier_external, gate)?;
        let (coupling, (surface_reports, carrier_reports), _) = self
            .partial_saturation_resistance_transition(
                fractions,
                bodies,
                external,
                origin,
                gate,
                Some((model, radii, &native)),
                true,
                |_, _, _, planes| collect(surfaces, carriers, planes, gate),
            )?;
        // Every fallible operation and callback completed before publication.
        for (body, report) in carriers.iter_mut().zip(&carrier_reports) {
            *body = report.motion.body_after;
        }
        Ok(SharedCarrierLubricationExchange3 {
            coupling,
            surfaces: surface_reports,
            carriers: carrier_reports,
        })
    }
}

pub(in crate::d3q19::boundary) fn prepare(
    surfaces: &[SharedLubricationSurface3],
    carriers: &[PartialSaturationBody3],
    external: &[[f64; 6]],
    gate: &CancelGate,
) -> Result<Vec<LubricationPlane3>> {
    checkpoint(gate)?;
    if surfaces.len() > 64
        || carriers.len() > 64
        || external.len() != carriers.len()
        || !external.iter().flatten().all(|x| x.is_finite())
    {
        return Err(invalid(
            "invalid shared carrier surface count or external loads",
        ));
    }
    let mut anchors = [None; 64];
    let mut native = Vec::new();
    native
        .try_reserve_exact(surfaces.len())
        .map_err(|_| invalid("shared carrier descriptor allocation refused"))?;
    for (index, surface) in surfaces.iter().copied().enumerate() {
        checkpoint(gate)?;
        native.push(match surface {
            SharedLubricationSurface3::Prescribed(plane) => {
                if plane.recoil.is_some() {
                    return Err(invalid(
                        "prescribed surface contains finite carrier metadata",
                    ));
                }
                plane
            }
            SharedLubricationSurface3::Carried(face) => {
                let body = *carriers
                    .get(face.carrier)
                    .ok_or_else(|| invalid("shared face references a missing carrier"))?;
                let anchor = *anchors[face.carrier].get_or_insert(index);
                // Geometry was admitted at face construction. Preserve its
                // normalization; the common solver admits/factors this one
                // carrier inertia once, regardless of its number of faces.
                LubricationPlane3 {
                    velocity: body.velocity,
                    recoil: Some(finite_planes::FinitePlaneInput {
                        body,
                        external: external[face.carrier],
                        carrier_anchor: Some(anchor),
                    }),
                    ..face.geometry
                }
            }
        });
    }
    if anchors[..carriers.len()].iter().any(Option::is_none) {
        return Err(invalid("every finite carrier must own a film face"));
    }
    checkpoint(gate)?;
    Ok(native)
}

pub(in crate::d3q19::boundary) fn collect(
    surfaces: &[SharedLubricationSurface3],
    carriers: &[PartialSaturationBody3],
    reports: &[PlaneLubricationExchange3],
    gate: &CancelGate,
) -> Result<(Vec<CarrierFaceExchange3>, Vec<CarrierLubricationExchange3>)> {
    checkpoint(gate)?;
    if surfaces.len() != reports.len() || carriers.len() > 64 {
        return Err(invalid("shared carrier response shape changed"));
    }
    let mut faces = Vec::new();
    faces
        .try_reserve_exact(surfaces.len())
        .map_err(|_| invalid("shared face receipt allocation refused"))?;
    let mut groups = Vec::new();
    groups
        .try_reserve_exact(carriers.len())
        .map_err(|_| invalid("shared carrier receipt allocation refused"))?;
    for &body in carriers {
        groups.push(CarrierLubricationExchange3 {
            impulse: [0.0; 3],
            angular_impulse: [0.0; 3],
            work: 0.0,
            motion: FinitePlaneResponse3::initial(body),
            active_films: 0,
            tangential_films: 0,
        });
    }
    let mut seen = [false; 64];
    let mut impulses = [[Dd::ZERO; 3]; 64];
    let mut angular = [[Dd::ZERO; 3]; 64];
    let mut spin = [[Dd::ZERO; 3]; 64];
    let mut work = [Dd::ZERO; 64];
    for (&surface, report) in surfaces.iter().zip(reports) {
        checkpoint(gate)?;
        let carrier = match surface {
            SharedLubricationSurface3::Prescribed(_) => {
                if report.finite_body.is_some() {
                    return Err(invalid("prescribed face returned finite carrier storage"));
                }
                None
            }
            SharedLubricationSurface3::Carried(face) => {
                let motion = report
                    .finite_body
                    .ok_or_else(|| invalid("shared face lost its finite response"))?;
                let group = groups
                    .get_mut(face.carrier)
                    .ok_or_else(|| invalid("shared response references a missing carrier"))?;
                if !seen[face.carrier] {
                    group.motion = motion;
                    seen[face.carrier] = true;
                }
                for i in 0..3 {
                    impulses[face.carrier][i] =
                        impulses[face.carrier][i] + Dd::from_f64(report.impulse[i]);
                    angular[face.carrier][i] =
                        angular[face.carrier][i] + Dd::from_f64(report.angular_impulse[i]);
                    spin[face.carrier][i] =
                        spin[face.carrier][i] + Dd::from_f64(motion.spin_impulse[i]);
                }
                work[face.carrier] = work[face.carrier] + Dd::from_f64(report.work);
                group.active_films += report.active_spheres;
                group.tangential_films += report.tangential_spheres;
                Some(face.carrier)
            }
        };
        faces.push(CarrierFaceExchange3 {
            carrier,
            impulse: report.impulse,
            angular_impulse: report.angular_impulse,
            couple_impulse: report.couple_impulse,
            work: report.work,
            active_spheres: report.active_spheres,
            tangential_spheres: report.tangential_spheres,
        });
    }
    for (index, group) in groups.iter_mut().enumerate() {
        checkpoint(gate)?;
        if !seen[index] {
            return Err(invalid("shared carrier response has no owning face"));
        }
        group.impulse = impulses[index].map(Dd::to_f64);
        group.angular_impulse = angular[index].map(Dd::to_f64);
        group.motion.spin_impulse = spin[index].map(Dd::to_f64);
        group.work = work[index].to_f64();
        group.motion.work_residual =
            group.motion.kinetic_energy_change - group.motion.external_work - group.work;
        if !group
            .impulse
            .into_iter()
            .chain(group.angular_impulse)
            .chain(group.motion.spin_impulse)
            .chain([group.work, group.motion.work_residual])
            .all(f64::is_finite)
        {
            return Err(invalid("shared carrier receipt overflowed"));
        }
    }
    checkpoint(gate)?;
    Ok((faces, groups))
}

#[cfg(test)]
mod tests;
