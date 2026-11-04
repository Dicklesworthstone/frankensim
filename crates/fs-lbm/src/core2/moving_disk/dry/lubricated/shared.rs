//! Several missing-film patches on one translating finite inertia.

use super::translating::{TranslationResponse2, admitted_velocity, normal_exchange};
use super::*;

const MAX_CARRIERS: usize = 64;

/// One finite mass shared by any number of attached planar film patches.
///
/// Motion is translation in world lattice coordinates, per unit out-of-plane
/// depth. Orientation is held by a zero-work guide; no rotary inertia,
/// displaced volume, buoyancy or resolved-fluid wall is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TranslatingFilmCarrier2 {
    mass: f64,
    center: [f64; 2],
    velocity: [f64; 2],
}
impl TranslatingFilmCarrier2 {
    /// Admit positive mass, finite center and a complete low-Mach velocity vector.
    pub fn new(mass: f64, center: [f64; 2], velocity: [f64; 2]) -> Result<Self, DryDiskError2> {
        if !(mass.is_finite() && mass > 0.0)
            || !center.into_iter().all(f64::is_finite)
            || !admitted_velocity(velocity)
        {
            return Err(DryDiskError2::Invalid(
                "invalid shared film carrier mass or motion",
            ));
        }
        Ok(Self {
            mass,
            center,
            velocity,
        })
    }
    /// Physical mass per unit out-of-plane depth.
    pub const fn mass(self) -> f64 {
        self.mass
    }
    /// Actual accepted world center of mass.
    pub const fn center(self) -> [f64; 2] {
        self.center
    }
    /// Actual accepted world velocity, including normal and tangential recoil.
    pub const fn velocity(self) -> [f64; 2] {
        self.velocity
    }
}

/// A planar film patch attached to a stable index in a finite carrier registry.
///
/// The point is relative to its carrier COM; normal and relative offset stay
/// fixed in world axes because the carrier's orientation is constrained.
/// Equal carrier values at different indices still mean distinct physical masses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierLubricationPatch2 {
    carrier: usize,
    local: LubricationPlane2,
}
impl CarrierLubricationPatch2 {
    /// Bind the original normal-film geometry and scales to a carrier index.
    /// Index admission occurs against the complete registry supplied to a step.
    pub fn new(
        carrier: usize,
        point_relative: [f64; 2],
        normal: [f64; 2],
        cutoff_gap: f64,
        minimum_gap: f64,
    ) -> Result<Self, DryDiskError2> {
        Ok(Self {
            carrier,
            local: LubricationPlane2::new(
                point_relative,
                normal,
                [0.0; 2],
                cutoff_gap,
                minimum_gap,
            )?,
        })
    }
    /// Enable the original independent cylinder sliding and spin resistance.
    pub fn with_tangential(
        mut self,
        cutoff_gap: f64,
        minimum_gap: f64,
    ) -> Result<Self, DryDiskError2> {
        self.local = self.local.with_tangential(cutoff_gap, minimum_gap)?;
        Ok(self)
    }
    /// Index of the unique mechanical carrier.
    pub const fn carrier(self) -> usize {
        self.carrier
    }
    /// Fixed point relative to that carrier's center of mass.
    pub const fn point_relative(self) -> [f64; 2] {
        self.local.point
    }
    /// Actual world geometry and velocity at the supplied accepted registry.
    pub fn geometry(
        self,
        carriers: &[TranslatingFilmCarrier2],
    ) -> Result<LubricationPlane2, DryDiskError2> {
        let carrier = carriers.get(self.carrier).ok_or(DryDiskError2::Invalid(
            "film patch carrier index is out of range",
        ))?;
        self.at(*carrier)
    }
    fn at(self, carrier: TranslatingFilmCarrier2) -> Result<LubricationPlane2, DryDiskError2> {
        let point = core::array::from_fn(|a| carrier.center[a] + self.local.point[a]);
        if !point.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid(
                "shared film patch geometry overflowed",
            ));
        }
        Ok(LubricationPlane2 {
            point,
            velocity: carrier.velocity,
            ..self.local
        })
    }
}

/// One unique carrier's aggregate film transfer, trajectory and storage.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct CarrierLubricationExchange2 {
    /// Complete endpoint assigned once to the registry.
    pub state_after: TranslatingFilmCarrier2,
    /// Common carrier COM at the disk/fluid/film impulse.
    pub collision_center: [f64; 2],
    /// One solved world midpoint velocity used by every attached patch.
    pub velocity_midpoint: [f64; 2],
    /// Sum of received patch impulses [Px, Py, Lz] about the requested origin.
    pub film_impulse: [f64; 3],
    /// Sum of received patch work, entering this carrier's kinetic storage.
    pub film_work: f64,
    /// Actual [Px, Py, Lz] change along both free half-drifts.
    pub momentum_change: [f64; 3],
    /// Authored center impulse and its moment about the same fixed origin.
    pub external_impulse: [f64; 3],
    /// Pure guide couple DELIVERED to this carrier; no supported linear force.
    /// Its work is zero because orientation remains fixed.
    pub guide_couple: f64,
    /// Work of the once-per-carrier authored impulse at the common midpoint.
    pub external_work: f64,
    /// Actual change in the one carrier's kinetic energy.
    pub kinetic_energy_change: f64,
    /// Actual momentum change minus received film, authored and guide inputs.
    pub momentum_residual: [f64; 3],
    /// Kinetic increment minus authored work minus aggregate received film work.
    pub work_residual: f64,
    /// Prepared 1-norm condition estimate of the scaled two-direction response.
    pub response_condition_number: f64,
}

/// One dry-disk/fluid step with distinct patch and carrier accounts.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DryLubricatedCarriersExchange2 {
    /// Original disk/fluid result. Its planes hold patch-only wrench and work
    /// in input patch order; they do not repeat carrier kinetic storage.
    pub motion: DryLubricatedDiskExchange2,
    /// One aggregate physical trajectory and storage account per registry entry.
    pub carriers: Vec<CarrierLubricationExchange2>,
    /// Exterior fluid, dry disk and all carriers after subtracting authored
    /// loads, guide couples and the original fluid/boundary sources.
    pub combined_momentum_residual: [f64; 3],
}

impl DryDisk2 {
    /// Advance the native disk/fluid and shared translating film carriers once.
    ///
    /// Each carrier owns ONE mass and authored WORLD center impulse. For all
    /// its attached patches, assemble R = sum(k_n n n^T + k_t t t^T).
    /// The shared midpoint solves (2m I+R) v_mid = 2m v_old+J_ext+R v_disk,mid.
    /// Its bounded two-by-two response is prepared once; every original
    /// dry-fluid candidate uses that same response for all attached patches.
    /// Independent per-patch mass elimination would omit cross-patch recoil.
    ///
    /// Original normal, sliding and cylinder-spin laws, cutoffs and floors
    /// remain authoritative. The orientation guide opposes the sum of all
    /// face moments about the actual common carrier COM. Patches follow its
    /// two free half-drifts; actual solved endpoints admit every relative gap.
    /// Unrepresentable COM/local cancellation that prevents the patch from
    /// following its common translation refuses at the local lattice scale.
    ///
    /// At most 64 patches and 64 carriers are allowed, with at least one patch
    /// per carrier. The existing maximum_condition_number also bounds the
    /// prepared response. Probes allocate no carrier vectors or factorizations.
    /// All finite endpoints and accounts finish before native publication.
    /// Any refusal/cancellation preserves all accepted fluid, disk and carriers.
    ///
    /// These local infinite patches add no resolved moving walls, face extents,
    /// contact, GCL, thermal destination or rotation. The native midpoint
    /// method remains non-L-stable; stiff relaxation needs time resolution.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubrication_carriers(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
        moment_origin: [f64; 2],
        dynamic_viscosity: f64,
        carriers: &mut [TranslatingFilmCarrier2],
        carrier_external_impulses: &[[f64; 2]],
        patches: &[CarrierLubricationPatch2],
        config: DryDiskCouplingConfig2,
        cancel: &CancelGate,
    ) -> Result<DryLubricatedCarriersExchange2, DryDiskError2> {
        checkpoint(cancel)?;
        config.admit()?;
        if !(dynamic_viscosity.is_finite() && dynamic_viscosity > 0.0)
            || !moment_origin.into_iter().all(f64::is_finite)
            || carriers.len() > MAX_CARRIERS
            || patches.len() > 64
            || carrier_external_impulses.len() != carriers.len()
        {
            return Err(DryDiskError2::Invalid(
                "invalid shared carrier viscosity, origin or count",
            ));
        }
        let descriptors = patches
            .iter()
            .map(|patch| patch.geometry(carriers))
            .collect::<Result<Vec<_>, _>>()?;
        let body = self.body();
        let disk_center = core::array::from_fn(|a| body.center[a] + 0.5 * body.velocity[a]);
        let radius = self.radius();
        let mut prepared = PreparedPlanes2 {
            planes: Vec::with_capacity(patches.len()),
            shared: None,
        };
        for (i, &patch) in patches.iter().enumerate() {
            checkpoint(cancel)?;
            PreparedPlanes2::admit_initial_plane(&descriptors, i, body.center, radius)?;
            let initial = carriers[patch.carrier];
            let midpoint = TranslatingFilmCarrier2 {
                center: core::array::from_fn(|a| initial.center[a] + 0.5 * initial.velocity[a]),
                ..initial
            };
            let geometry = patch.at(midpoint)?;
            admit_patch_drift(
                descriptors[i].point,
                geometry.point,
                initial.velocity,
                disk_center,
                radius,
            )?;
            prepared.planes.push(PreparedPlane2::at_midpoint(
                geometry,
                disk_center,
                radius,
                dynamic_viscosity,
                None,
            )?);
        }
        prepared.shared = Some(PreparedCarriers2::new(
            carriers,
            carrier_external_impulses,
            patches,
            &mut prepared.planes,
            config.maximum_condition_number,
            cancel,
        )?);
        let report = self.step_inner(
            grid,
            external_impulse,
            external_couple,
            moment_origin,
            config,
            cancel,
            Some(&prepared),
            |motion| prepared.finish_carriers(motion, moment_origin, cancel),
        )?;
        for (carrier, receipt) in carriers.iter_mut().zip(&report.carriers) {
            *carrier = receipt.state_after;
        }
        Ok(report)
    }
}

pub(super) struct PreparedCarrierPatch2 {
    patch: CarrierLubricationPatch2,
    normal_reference: UnderResolvedLubricationExchange2,
}
pub(super) struct PreparedCarriers2 {
    carriers: Vec<PreparedCarrier2>,
}
struct PreparedCarrier2 {
    initial: TranslatingFilmCarrier2,
    center: [f64; 2],
    external: [f64; 2],
    free_midpoint: [f64; 2],
    scale: f64,
    mass_scaled: f64,
    inverse: [f64; 4],
    retained: [f64; 4],
    recoil: [f64; 4],
    active: bool,
    condition: f64,
}
#[derive(Clone, Copy, Default)]
pub(super) struct CarrierResponse2 {
    midpoint: [f64; 2],
    solution: [f64; 2],
}

impl PreparedCarriers2 {
    fn new(
        carriers: &[TranslatingFilmCarrier2],
        external: &[[f64; 2]],
        patches: &[CarrierLubricationPatch2],
        planes: &mut [PreparedPlane2],
        maximum_condition: f64,
        cancel: &CancelGate,
    ) -> Result<Self, DryDiskError2> {
        let mut scales = [0.0_f64; MAX_CARRIERS];
        let mut used = [false; MAX_CARRIERS];
        let mut active = [false; MAX_CARRIERS];
        for (i, carrier) in carriers.iter().enumerate() {
            scales[i] = carrier.mass;
        }
        for (plane, &patch) in planes.iter_mut().zip(patches) {
            checkpoint(cancel)?;
            let mut pair = plane.pair;
            pair.surface_velocities = [[0.0; 2]; 2];
            let reference = pair.correction([0.0; 2]);
            let owner = patch.carrier;
            used[owner] = true;
            active[owner] |=
                reference.resistance_increment > 0.0 || plane.tangential_resistance > 0.0;
            scales[owner] = scales[owner]
                .max(reference.resistance_increment)
                .max(plane.tangential_resistance);
            plane.carrier = Some(PreparedCarrierPatch2 {
                patch,
                normal_reference: reference,
            });
        }
        if used[..carriers.len()].iter().any(|&value| !value) {
            return Err(DryDiskError2::Invalid(
                "every finite carrier must own at least one patch",
            ));
        }
        let mut gram = [[Dd::ZERO; 3]; MAX_CARRIERS];
        for plane in planes.iter() {
            checkpoint(cancel)?;
            let patch = plane.carrier.as_ref().expect("all shared patches prepared");
            let owner = patch.patch.carrier;
            for (k, direction) in [
                (
                    patch.normal_reference.resistance_increment,
                    patch.patch.local.normal,
                ),
                (plane.tangential_resistance, plane.tangent),
            ] {
                let weight = 0.5 * (k / scales[owner]);
                let values = [
                    weight * direction[0] * direction[0],
                    weight * direction[0] * direction[1],
                    weight * direction[1] * direction[1],
                ];
                for (sum, value) in gram[owner].iter_mut().zip(values) {
                    *sum = *sum + Dd::from_f64(value);
                }
            }
        }
        let mut prepared = Vec::with_capacity(carriers.len());
        for (i, &initial) in carriers.iter().enumerate() {
            checkpoint(cancel)?;
            let scale = scales[i];
            let mass_scaled = initial.mass / scale;
            let g = gram[i].map(Dd::to_f64);
            let resistance = [g[0], g[1], g[1], g[2]];
            let matrix = [mass_scaled + g[0], g[1], g[1], mass_scaled + g[2]];
            let factor = fs_la::factor::cholesky(&matrix, 2).map_err(|_| {
                DryDiskError2::Invalid("shared carrier response is not positive definite")
            })?;
            let mut first = [1.0, 0.0];
            let mut second = [0.0, 1.0];
            factor.solve(&mut first);
            factor.solve(&mut second);
            let off_diagonal = 0.5 * first[1] + 0.5 * second[0];
            let inverse = [first[0], off_diagonal, off_diagonal, second[1]];
            let condition = matrix_norm(matrix) * matrix_norm(inverse);
            if !condition.is_finite() || condition > maximum_condition {
                return Err(DryDiskError2::Invalid(
                    "shared carrier response is ill-conditioned",
                ));
            }
            let retained = inverse.map(|x| mass_scaled * x);
            let recoil = matrix_product(inverse, resistance);
            let free_midpoint = core::array::from_fn(|a| {
                initial.velocity[a] + (0.5 * external[i][a]) / initial.mass
            });
            let center = core::array::from_fn(|a| initial.center[a] + 0.5 * initial.velocity[a]);
            if !external[i]
                .into_iter()
                .chain(free_midpoint)
                .chain(center)
                .chain(retained)
                .chain(recoil)
                .all(f64::is_finite)
            {
                return Err(DryDiskError2::Invalid("shared carrier response overflowed"));
            }
            prepared.push(PreparedCarrier2 {
                initial,
                center,
                external: external[i],
                free_midpoint,
                scale,
                mass_scaled,
                inverse,
                retained,
                recoil,
                active: active[i],
                condition,
            });
        }
        Ok(Self { carriers: prepared })
    }

    fn responses(
        &self,
        velocity: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<[CarrierResponse2; MAX_CARRIERS], DryDiskError2> {
        let mut out = [CarrierResponse2::default(); MAX_CARRIERS];
        for (i, carrier) in self.carriers.iter().enumerate() {
            checkpoint(cancel)?;
            let delta = subtract(velocity, carrier.free_midpoint);
            let solution = matrix_vector(carrier.inverse, delta);
            let midpoint = if !carrier.active || delta == [0.0; 2] {
                carrier.free_midpoint
            } else {
                let free = matrix_vector(carrier.retained, carrier.free_midpoint);
                let disk = matrix_vector(carrier.recoil, velocity);
                core::array::from_fn(|a| (Dd::from_f64(free[a]) + Dd::from_f64(disk[a])).to_f64())
            };
            if !midpoint.into_iter().chain(solution).all(f64::is_finite) {
                return Err(DryDiskError2::Invalid("shared carrier midpoint overflowed"));
            }
            out[i] = CarrierResponse2 { midpoint, solution };
        }
        Ok(out)
    }
}

impl PreparedCarrier2 {
    fn endpoint(
        &self,
        response: CarrierResponse2,
    ) -> Result<TranslatingFilmCarrier2, DryDiskError2> {
        if !admitted_velocity(response.midpoint) {
            return Err(DryDiskError2::Invalid(
                "shared carrier midpoint exceeds low-Mach limit",
            ));
        }
        let velocity =
            core::array::from_fn(|a| 2.0 * response.midpoint[a] - self.initial.velocity[a]);
        let center = core::array::from_fn(|a| self.center[a] + 0.5 * velocity[a]);
        TranslatingFilmCarrier2::new(self.initial.mass, center, velocity)
    }

    fn mode(&self, response: CarrierResponse2, k: f64, direction: [f64; 2]) -> (f64, f64) {
        let projected = dot(response.solution, direction);
        let rate = self.mass_scaled * projected;
        // Keep the finite force when k/m is large without multiplying a huge
        // original resistance by an already-underflowed relative velocity.
        let impulse = if k <= self.initial.mass {
            -k * rate
        } else {
            -self.initial.mass * ((k / self.scale) * projected)
        };
        (impulse, rate)
    }

    fn finish(
        &self,
        response: CarrierResponse2,
        totals: [f64; 5],
        origin: [f64; 2],
    ) -> Result<CarrierLubricationExchange2, DryDiskError2> {
        let after = self.endpoint(response)?;
        let delta_v = subtract(after.velocity, self.initial.velocity);
        let delta = delta_v.map(|v| self.initial.mass * v);
        let old_momentum = self.initial.velocity.map(|v| self.initial.mass * v);
        let momentum_change = [
            delta[0],
            delta[1],
            cross(subtract(after.center, origin), delta)
                + cross(subtract(after.center, self.initial.center), old_momentum),
        ];
        let arm = subtract(self.center, origin);
        let external_impulse = [
            self.external[0],
            self.external[1],
            cross(arm, self.external),
        ];
        let film_impulse = [totals[0], totals[1], totals[2]];
        let film_work = totals[3];
        let guide_couple = totals[4];
        let momentum_residual = core::array::from_fn(|a| {
            momentum_change[a]
                - film_impulse[a]
                - external_impulse[a]
                - if a == 2 { guide_couple } else { 0.0 }
        });
        let external_work = dot(self.external, response.midpoint);
        let sum_v = core::array::from_fn(|a| after.velocity[a] + self.initial.velocity[a]);
        let kinetic_energy_change = 0.5 * self.initial.mass * dot(delta_v, sum_v);
        let work_residual = kinetic_energy_change - external_work - film_work;
        let scale = norm(momentum_change)
            .max(norm(film_impulse))
            .max(norm(external_impulse))
            .max(guide_couple.abs())
            .max(norm([
                old_momentum[0],
                old_momentum[1],
                cross(arm, old_momentum),
            ]));
        let energy_scale = kinetic_energy_change
            .abs()
            .max(external_work.abs())
            .max(film_work.abs())
            .max(self.initial.mass * dot(self.initial.velocity, self.initial.velocity))
            .max(self.initial.mass * dot(after.velocity, after.velocity));
        if !momentum_change
            .into_iter()
            .chain(external_impulse)
            .chain(momentum_residual)
            .chain(totals)
            .chain([
                external_work,
                kinetic_energy_change,
                work_residual,
                scale,
                energy_scale,
            ])
            .all(f64::is_finite)
            || norm(momentum_residual) > 16384.0 * f64::EPSILON * scale
            || work_residual.abs() > 16384.0 * f64::EPSILON * energy_scale
        {
            return Err(DryDiskError2::Invalid(
                "shared carrier momentum/work balance is unresolved",
            ));
        }
        Ok(CarrierLubricationExchange2 {
            state_after: after,
            collision_center: self.center,
            velocity_midpoint: response.midpoint,
            film_impulse,
            film_work,
            momentum_change,
            external_impulse,
            guide_couple,
            external_work,
            kinetic_energy_change,
            momentum_residual,
            work_residual,
            response_condition_number: self.condition,
        })
    }
}

impl PreparedPlanes2 {
    pub(super) fn shared_responses(
        &self,
        velocity: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<Option<[CarrierResponse2; MAX_CARRIERS]>, DryDiskError2> {
        self.shared
            .as_ref()
            .map(|shared| shared.responses(velocity, cancel))
            .transpose()
    }

    pub(super) fn exchange_at(
        &self,
        index: usize,
        velocity: [f64; 2],
        angular_velocity: f64,
        origin: [f64; 2],
        responses: Option<&[CarrierResponse2; MAX_CARRIERS]>,
    ) -> Result<PlaneLubricationExchange2, DryDiskError2> {
        let plane = &self.planes[index];
        let Some(shared) = &self.shared else {
            return plane.exchange(velocity, angular_velocity, origin);
        };
        let patch = plane
            .carrier
            .as_ref()
            .ok_or(DryDiskError2::Invalid("missing shared film patch"))?;
        let states = responses.ok_or(DryDiskError2::Invalid("missing shared carrier response"))?;
        let owner = patch.patch.carrier;
        let carrier = &shared.carriers[owner];
        let state = states[owner];
        let (normal_impulse, normal_rate) = carrier.mode(
            state,
            patch.normal_reference.resistance_increment,
            patch.patch.local.normal,
        );
        let (tangent_impulse, tangent_rate) =
            carrier.mode(state, plane.tangential_resistance, plane.tangent);
        let response = TranslationResponse2 {
            velocity_midpoint: state.midpoint,
            normal_impulse,
            normal_rate,
            tangent_impulse,
            tangent_rate,
            tangent_midpoint: dot(state.midpoint, plane.tangent),
        };
        let normal = normal_exchange(
            patch.patch.local.normal,
            plane.pair.line_of_action_point,
            patch.normal_reference,
            response,
            origin,
        )?;
        plane.exchange_with_response(velocity, angular_velocity, origin, Some((response, normal)))
    }

    pub(super) fn shared_endpoint(
        &self,
        plane: &PreparedPlane2,
        responses: &[CarrierResponse2; MAX_CARRIERS],
    ) -> Result<LubricationPlane2, DryDiskError2> {
        let shared = self
            .shared
            .as_ref()
            .ok_or(DryDiskError2::Invalid("missing carrier registry"))?;
        let patch = plane
            .carrier
            .as_ref()
            .ok_or(DryDiskError2::Invalid("missing shared film patch"))?;
        let owner = patch.patch.carrier;
        let carrier = shared.carriers[owner].endpoint(responses[owner])?;
        let geometry = patch.patch.at(carrier)?;
        admit_patch_drift(
            plane.midpoint_point,
            geometry.point,
            carrier.velocity,
            plane.pair.line_of_action_point,
            plane.pair.reduced_radius,
        )?;
        Ok(geometry)
    }

    fn finish_carriers(
        &self,
        motion: DryLubricatedDiskExchange2,
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<DryLubricatedCarriersExchange2, DryDiskError2> {
        let shared = self
            .shared
            .as_ref()
            .ok_or(DryDiskError2::Invalid("missing carrier registry"))?;
        let responses = shared.responses(motion.motion.native.velocity_midpoint, cancel)?;
        let mut sums = [[Dd::ZERO; 5]; MAX_CARRIERS];
        for (plane, film) in self.planes.iter().zip(&motion.planes) {
            checkpoint(cancel)?;
            let patch = plane
                .carrier
                .as_ref()
                .ok_or(DryDiskError2::Invalid("missing shared film patch"))?;
            let owner = patch.patch.carrier;
            let force = [film.plane_impulse[0], film.plane_impulse[1]];
            let guide_couple = -film.plane_couple
                - cross(
                    subtract(film.plane_footpoint, shared.carriers[owner].center),
                    force,
                );
            let values = [
                force[0],
                force[1],
                film.plane_impulse[2],
                film.plane_work,
                guide_couple,
            ];
            for (sum, value) in sums[owner].iter_mut().zip(values) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let mut carriers = Vec::with_capacity(shared.carriers.len());
        for (i, carrier) in shared.carriers.iter().enumerate() {
            checkpoint(cancel)?;
            carriers.push(carrier.finish(responses[i], sums[i].map(Dd::to_f64), origin)?);
        }
        let mut combined = motion.motion.combined_momentum_residual;
        for a in 0..3 {
            let mut sum = Dd::from_f64(combined[a]) + Dd::from_f64(motion.body_impulse[a]);
            for carrier in &carriers {
                sum = sum + Dd::from_f64(carrier.momentum_change[a])
                    - Dd::from_f64(carrier.external_impulse[a])
                    - Dd::from_f64(if a == 2 { carrier.guide_couple } else { 0.0 });
            }
            combined[a] = sum.to_f64();
        }
        if !combined.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid(
                "shared carrier aggregate balance overflowed",
            ));
        }
        Ok(DryLubricatedCarriersExchange2 {
            motion,
            carriers,
            combined_momentum_residual: combined,
        })
    }
}

fn matrix_vector(a: [f64; 4], x: [f64; 2]) -> [f64; 2] {
    [
        a[0].mul_add(x[0], a[1] * x[1]),
        a[2].mul_add(x[0], a[3] * x[1]),
    ]
}
fn matrix_product(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    core::array::from_fn(|i| {
        let row = i / 2;
        let column = i % 2;
        (Dd::from_f64(a[2 * row]) * Dd::from_f64(b[column])
            + Dd::from_f64(a[2 * row + 1]) * Dd::from_f64(b[2 + column]))
        .to_f64()
    })
}
fn matrix_norm(a: [f64; 4]) -> f64 {
    (a[0].abs() + a[2].abs()).max(a[1].abs() + a[3].abs())
}

fn admit_patch_drift(
    before: [f64; 2],
    after: [f64; 2],
    velocity: [f64; 2],
    disk_center: [f64; 2],
    radius: f64,
) -> Result<(), DryDiskError2> {
    // COM/local cancellation must not freeze a nominally moving face. Use
    // the local lattice-interface scale, not the possibly enormous COM.
    let scale = radius
        .max(disk_center[0].abs())
        .max(disk_center[1].abs())
        .max(1.0);
    if (0..2)
        .any(|a| ((after[a] - before[a]) - 0.5 * velocity[a]).abs() > 128.0 * f64::EPSILON * scale)
    {
        return Err(DryDiskError2::Invalid(
            "shared patch translation is unrepresentable",
        ));
    }
    Ok(())
}
