//! Source-backed capillary adhesion between parallel, identically wetting plates.
//!
//! This is the exact two-dimensional, gravity-free circular-meniscus model:
//! an invariant out-of-plane direction (reported over a supplied width), fixed
//! liquid volume, equal acute Young angles, and quasistatic mobile contact lines.
//! The actual source owners remain immutable. With A=V/width, c=cos(theta),
//! s=sin(theta), alpha=pi/2-theta, the wetted span is A/h+k h, where
//! k=(alpha-sc)/(2c²). The same surface free energy gives
//! E=gamma width[-2c A/h+(alpha/c+s)h] and F_normal=-dE/dh.
//!
//! The chosen topology rule keeps a connected profile until its geometric neck
//! vanishes, then relaxes into two equal-volume circular caps. On closing, caps
//! coalesce when they touch. These are explicit symmetric quasistatic topology
//! rules, NOT a certified stability threshold or a resolved pinch-off mechanism.
//! Released free energy is exported, not silently deposited as heat. Gravity,
//! inertia and rate screens are applicability limits, not accuracy certificates.
//! No finite-width end curvature, three-dimensional instability, unequal liquid
//! partition, evaporation, cavitation, roughness, pinning or dynamic contact angle is inferred.
//!
//! Geometry: Teixeira & Teixeira, J. Phys. Condens. Matter 32 (2020) 034002,
//! https://doi.org/10.1088/1361-648X/ab48b7; the zero-gravity circular profile and
//! zero-neck condition also appear in Yang & Hong, arXiv:2411.15021, App. B/C.
//! The mobile-Young-angle energy and two-cap transitions here are explicit
//! derivations within that restricted geometry, not empirical rupture fits.

use crate::normal_patch::PointNormalPort;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_material::liquid::surface::adhesion::YoungDupreLaw;
use fs_material::liquid::surface::{ResolvedLiquidSurfaceState, ResolvedLiquidWettingState};
use fs_material::state_point::MaterialStatePointError;

mod geometry;

/// Physical geometry, finite-work and quasistatic admission supplied by caller.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeLimits {
    pub minimum_gap_m: f64,
    pub maximum_gap_m: f64,
    /// Largest wet footprint along either plate, including detached caps.
    pub maximum_wetted_span_m: f64,
    pub maximum_abs_force_n: f64,
    pub maximum_abs_surface_energy_j: f64,
    pub maximum_abs_step_work_j: f64,
    pub maximum_transition_release_j: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
    /// Actual gravity magnitude; the profile neglects it subject to the Bo screen.
    pub gravity_m_s2: f64,
    /// Each dimensionless ceiling must be in (0,0.01]. No accuracy is certified.
    pub maximum_bond_number: f64,
    pub maximum_capillary_number: f64,
    pub maximum_weber_number: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BridgeError {
    Material(MaterialStatePointError),
    Invalid {
        field: &'static str,
    },
    Limit {
        field: &'static str,
        value: f64,
        limit: f64,
    },
    StateMismatch,
    NoConnectedProfile,
    OverlappingCaps,
    /// The proposed instantaneous topology switch would require outside energy.
    EndothermicTransition {
        energy_increase_j: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl From<MaterialStatePointError> for BridgeError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for BridgeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar capillary bridge refused: {self:?}")
    }
}
impl std::error::Error for BridgeError {}

/// Immutable accepted branch, bound to the original model and source snapshots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeState {
    identity: ContentHash,
    gap_m: f64,
    attached: bool,
}
impl BridgeState {
    #[must_use]
    pub const fn gap_m(self) -> f64 {
        self.gap_m
    }
    #[must_use]
    pub const fn attached(self) -> bool {
        self.attached
    }
}

/// Actual geometry and signed normal force usable by a mechanical owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeProfile {
    pub gap_m: f64,
    pub attached: bool,
    pub liquid_volume_m3: f64,
    /// Rupture conserves all volume in two explicitly retained equal caps.
    pub detached_cap_volumes_m3: Option<[f64; 2]>,
    pub wetted_span_m: f64,
    pub neck_width_m: Option<f64>,
    pub liquid_gas_area_m2: f64,
    /// Total wetted area on BOTH plates.
    pub wetted_area_m2: f64,
    pub liquid_minus_gas_pressure_pa: f64,
    /// Young--Dupre surface energy relative to the initially dry substrates.
    pub surface_free_energy_j: f64,
    /// Force on the upper plate along increasing gap; negative attracts.
    pub normal_force_n: f64,
    /// dF_normal/d(gap), at fixed volume and on this smooth branch.
    pub force_gap_derivative_n_per_m: f64,
}
impl BridgeProfile {
    /// Equal/opposite resultant at one common origin. The normal points from
    /// lower to upper plate. This does not resolve spatial end effects. The
    /// existing normal port's tangent is with respect to APPROACH=-gap.
    /// For a finite path, use `BridgeStep::mean_normal_force_n`; this endpoint
    /// force must not replace the work-conjugate force across a topology event.
    pub fn point_port(
        &self,
        normal: [f64; 3],
        arm_m: [f64; 3],
    ) -> Result<PointNormalPort, BridgeError> {
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !length.is_finite()
            || (length - 1.0).abs() > 1e-12
            || arm_m.iter().any(|x| !x.is_finite())
            || self.neck_width_m == Some(0.0)
        {
            return Err(BridgeError::Invalid {
                field: "smooth bridge port frame",
            });
        }
        let force = normal.map(|n| n * self.normal_force_n);
        let moment = [
            arm_m[1] * force[2] - arm_m[2] * force[1],
            arm_m[2] * force[0] - arm_m[0] * force[2],
            arm_m[0] * force[1] - arm_m[1] * force[0],
        ];
        for value in force.iter().chain(moment.iter()) {
            finite(*value, "bridge wrench")?;
        }
        Ok(PointNormalPort {
            action_force_n: force,
            action_moment_n_m: moment,
            reaction_force_n: force.map(|x| -x),
            reaction_moment_n_m: moment.map(|x| -x),
            residual_force_n: [0.0; 3],
            tangent_n_per_m: -self.force_gap_derivative_n_per_m,
            dissipated_power_w: 0.0,
            irreversible_work_j: 0.0,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BridgeEvent {
    None,
    Ruptured,
    Reformed,
}

/// One monotone-gap proposal. Publish its state only with accepted mechanics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeStep {
    state: BridgeState,
    pub profile: BridgeProfile,
    pub event: BridgeEvent,
    /// Integral of dE/dgap along connected portions, excluding topology drops.
    /// Positive when mechanical work is supplied TO the liquid surface store.
    pub work_on_liquid_j: f64,
    /// Work-conjugate force over the entire gap increment: -work/Delta_gap.
    pub mean_normal_force_n: f64,
    /// Nonnegative unresolved free-energy release during the chosen relaxation.
    /// Its mechanical/thermal destination is external; it is NOT automatic heat.
    pub irreversible_release_j: f64,
    pub surface_energy_change_j: f64,
    /// Delta(surface energy) + release - work_on_liquid.
    pub energy_balance_residual_j: f64,
}
impl BridgeStep {
    #[must_use]
    pub const fn state(self) -> BridgeState {
        self.state
    }
}

/// Original liquid and both ordered substrate/liquid states, never refreshed or
/// copied per sample. All geometric evaluations and transitions allocate nothing.
#[derive(Clone, Debug)]
pub struct SourcePlanarCapillaryBridge {
    liquid: ResolvedLiquidSurfaceState,
    wetting: [ResolvedLiquidWettingState; 2],
    law: YoungDupreLaw,
    geometry: geometry::Geometry,
    limits: BridgeLimits,
    identity: ContentHash,
}
impl SourcePlanarCapillaryBridge {
    /// Bind the original liquid and both ordered surface systems at one frozen
    /// source point. Volume belongs to an invariant/periodic segment of the
    /// supplied width; it does not describe a finite bridge with curved ends.
    /// The caller's minimum gap must admit continuum capillarity without
    /// cavitation or disjoining-pressure effects. No source refresh is inferred.
    pub fn new(
        liquid: ResolvedLiquidSurfaceState,
        wetting: [ResolvedLiquidWettingState; 2],
        volume_m3: f64,
        invariant_width_m: f64,
        limits: BridgeLimits,
    ) -> Result<Self, BridgeError> {
        limits.validate()?;
        if liquid.rheology_axes().is_some() {
            return Err(BridgeError::Invalid {
                field: "frozen Newtonian liquid source",
            });
        }
        for boundary in &wetting {
            boundary.require_liquid_state(&liquid)?;
        }
        if wetting[0].contact_angle_radians() != wetting[1].contact_angle_radians()
            || wetting[0].context() != wetting[1].context()
        {
            return Err(BridgeError::Invalid {
                field: "symmetric wetting and common surface environment",
            });
        }
        let sigma = liquid.surface_tension_n_m().ok_or(BridgeError::Invalid {
            field: "active sourced liquid/gas tension",
        })?;
        // A directly measured Young angle also selects the ideal Young energy
        // relation explicitly; its derived work is never fabricated as a claim.
        let laws = [wetting[0].young_dupre_law(), wetting[1].young_dupre_law()];
        if let [Some(a), Some(b)] = laws {
            if a != b {
                return Err(BridgeError::Invalid {
                    field: "equal original wetting free energies",
                });
            }
        }
        let law = if let Some(law) = laws[0].or(laws[1]) {
            law
        } else {
            let cosine = fs_math::det::cos(wetting[0].contact_angle_radians());
            YoungDupreLaw::new(sigma * (1.0 + cosine), sigma).map_err(|_| BridgeError::Invalid {
                field: "representable ideal partial wetting",
            })?
        };
        let geometry = geometry::Geometry::new(volume_m3, invariant_width_m, law)?;
        positive(liquid.density_kg_m3() * volume_m3, "retained liquid mass")?;
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.planar-capillary-bridge.v1");
        hash.update(liquid.material().identity().as_bytes());
        for boundary in &wetting {
            hash.update(boundary.resolved().identity().as_bytes());
        }
        for value in [
            volume_m3,
            invariant_width_m,
            law.work_of_adhesion_j_m2(),
            limits.minimum_gap_m,
            limits.maximum_gap_m,
            limits.maximum_wetted_span_m,
            limits.maximum_abs_force_n,
            limits.maximum_abs_surface_energy_j,
            limits.maximum_abs_step_work_j,
            limits.maximum_transition_release_j,
            limits.energy_absolute_tolerance_j,
            limits.energy_relative_tolerance,
            limits.gravity_m_s2,
            limits.maximum_bond_number,
            limits.maximum_capillary_number,
            limits.maximum_weber_number,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        let result = Self {
            liquid,
            wetting,
            law,
            geometry,
            limits,
            identity: hash.finalize(),
        };
        // Neither topology switch may silently inject missing energy.
        result.transition_release(
            geometry
                .attached(geometry.rupture_gap)?
                .surface_free_energy_j,
            geometry.cap_energy,
        )?;
        result.transition_release(
            geometry.cap_energy,
            geometry
                .attached(geometry.formation_gap)?
                .surface_free_energy_j,
        )?;
        Ok(result)
    }
    #[must_use]
    pub const fn liquid(&self) -> &ResolvedLiquidSurfaceState {
        &self.liquid
    }
    #[must_use]
    pub const fn wetting(&self) -> &[ResolvedLiquidWettingState; 2] {
        &self.wetting
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn limits(&self) -> BridgeLimits {
        self.limits
    }
    #[must_use]
    pub const fn surface_law(&self) -> YoungDupreLaw {
        self.law
    }
    #[must_use]
    pub const fn rupture_gap_m(&self) -> f64 {
        self.geometry.rupture_gap
    }
    #[must_use]
    pub const fn formation_gap_m(&self) -> f64 {
        self.geometry.formation_gap
    }

    /// Select an explicit initial branch. Its energy is an initial condition,
    /// not an invented previous formation work. Detached caps cannot overlap.
    pub fn initialize(&self, gap_m: f64, attached: bool) -> Result<BridgeState, BridgeError> {
        let state = BridgeState {
            identity: self.identity,
            gap_m,
            attached,
        };
        if attached && gap_m >= self.geometry.rupture_gap {
            return Err(BridgeError::NoConnectedProfile);
        }
        self.evaluate(state)?;
        Ok(state)
    }
    pub fn evaluate(&self, state: BridgeState) -> Result<BridgeProfile, BridgeError> {
        if state.identity != self.identity {
            return Err(BridgeError::StateMismatch);
        }
        self.check_gap(state.gap_m)?;
        let profile = if state.attached {
            self.geometry.attached(state.gap_m)?
        } else {
            self.geometry.detached(state.gap_m)?
        };
        self.check_profile(&profile)?;
        Ok(profile)
    }

    /// Monotone, constant-rate gap motion over dt. The proposed state and all
    /// source owners remain untouched on source/domain/limit/cancellation refusal.
    pub fn advance(
        &self,
        previous: BridgeState,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<BridgeStep, BridgeError> {
        poll(gate)?;
        positive(dt_s, "bridge step duration")?;
        let before = self.evaluate(previous)?;
        self.check_gap(gap_m)?;
        let g = self.geometry;
        let (attached, event, from, to, release) = if previous.attached && gap_m >= g.rupture_gap {
            let end = g.attached(g.rupture_gap)?;
            self.check_profile(&end)?;
            (
                false,
                BridgeEvent::Ruptured,
                previous.gap_m,
                g.rupture_gap,
                self.transition_release(end.surface_free_energy_j, g.cap_energy)?,
            )
        } else if !previous.attached && gap_m <= g.formation_gap {
            let start = g.attached(g.formation_gap)?;
            self.check_profile(&start)?;
            (
                true,
                BridgeEvent::Reformed,
                g.formation_gap,
                gap_m,
                self.transition_release(g.cap_energy, start.surface_free_energy_j)?,
            )
        } else {
            (
                previous.attached,
                BridgeEvent::None,
                previous.gap_m,
                gap_m,
                0.0,
            )
        };
        let state = BridgeState {
            identity: self.identity,
            gap_m,
            attached,
        };
        let profile = self.evaluate(state)?;
        let connected = previous.attached || attached;
        let work = if connected { g.work(from, to)? } else { 0.0 };
        limit(
            work.abs(),
            self.limits.maximum_abs_step_work_j,
            "absolute bridge step work",
        )?;
        let speed = finite((gap_m - previous.gap_m).abs() / dt_s, "plate speed")?;
        if connected {
            self.check_rate(from, to, speed)?;
        }
        let mean_force = if gap_m == previous.gap_m {
            profile.normal_force_n
        } else {
            finite(
                -work / (gap_m - previous.gap_m),
                "work-conjugate normal force",
            )?
        };
        limit(
            mean_force.abs(),
            self.limits.maximum_abs_force_n,
            "mean bridge force",
        )?;
        let change = finite(
            profile.surface_free_energy_j - before.surface_free_energy_j,
            "surface energy change",
        )?;
        let residual = finite(change + release - work, "surface energy closure")?;
        let scale = before
            .surface_free_energy_j
            .abs()
            .max(profile.surface_free_energy_j.abs())
            .max(work.abs())
            .max(release);
        let tolerance = finite(
            self.limits.energy_absolute_tolerance_j + self.limits.energy_relative_tolerance * scale,
            "energy tolerance",
        )?;
        if residual.abs() > tolerance {
            return Err(BridgeError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        poll(gate)?;
        Ok(BridgeStep {
            state,
            profile,
            event,
            work_on_liquid_j: work,
            mean_normal_force_n: mean_force,
            irreversible_release_j: release,
            surface_energy_change_j: change,
            energy_balance_residual_j: residual,
        })
    }

    fn check_gap(&self, gap: f64) -> Result<(), BridgeError> {
        positive(gap, "plate gap")?;
        if gap < self.limits.minimum_gap_m {
            return Err(BridgeError::Invalid {
                field: "minimum continuum gap",
            });
        }
        limit(gap, self.limits.maximum_gap_m, "maximum plate gap")
    }
    fn check_profile(&self, p: &BridgeProfile) -> Result<(), BridgeError> {
        limit(
            p.wetted_span_m,
            self.limits.maximum_wetted_span_m,
            "wet footprint",
        )?;
        limit(
            p.normal_force_n.abs(),
            self.limits.maximum_abs_force_n,
            "bridge force",
        )?;
        limit(
            p.surface_free_energy_j.abs(),
            self.limits.maximum_abs_surface_energy_j,
            "absolute surface energy",
        )?;
        let scale = p.gap_m.max(p.wetted_span_m);
        let bond = if self.limits.gravity_m_s2 == 0.0 {
            0.0
        } else {
            finite(
                self.liquid.density_kg_m3() * self.limits.gravity_m_s2 * scale * scale
                    / self.law.surface_tension_n_m(),
                "Bond number",
            )?
        };
        limit(bond, self.limits.maximum_bond_number, "Bond number")
    }
    fn check_rate(&self, from: f64, to: f64, gap_speed: f64) -> Result<(), BridgeError> {
        let g = self.geometry;
        let line_speed = [from, to]
            .into_iter()
            .map(|h| 0.5 * (g.segment - (g.area / h) / h).abs() * gap_speed)
            .fold(gap_speed, f64::max);
        finite(line_speed, "maximum contact-line speed")?;
        let sigma = self.law.surface_tension_n_m();
        let ca = finite(
            self.liquid.dynamic_viscosity_pa_s() * line_speed / sigma,
            "capillary number",
        )?;
        let scale = from
            .max(to)
            .max(g.attached(from)?.wetted_span_m)
            .max(g.attached(to)?.wetted_span_m);
        let we = finite(
            self.liquid.density_kg_m3() * line_speed * line_speed * scale / sigma,
            "Weber number",
        )?;
        limit(ca, self.limits.maximum_capillary_number, "capillary number")?;
        limit(we, self.limits.maximum_weber_number, "Weber number")
    }
    fn transition_release(&self, before: f64, after: f64) -> Result<f64, BridgeError> {
        let released = finite(before - after, "topology free-energy release")?;
        if released < 0.0 {
            return Err(BridgeError::EndothermicTransition {
                energy_increase_j: -released,
            });
        }
        limit(
            released,
            self.limits.maximum_transition_release_j,
            "topology free-energy release",
        )?;
        Ok(released)
    }
}

impl BridgeLimits {
    fn validate(self) -> Result<(), BridgeError> {
        for value in [
            self.minimum_gap_m,
            self.maximum_gap_m,
            self.maximum_wetted_span_m,
            self.maximum_abs_force_n,
            self.maximum_abs_surface_energy_j,
            self.maximum_abs_step_work_j,
            self.maximum_transition_release_j,
            self.energy_absolute_tolerance_j,
            self.energy_relative_tolerance,
        ] {
            positive(value, "bridge limit")?;
        }
        if self.energy_relative_tolerance >= 1.0 {
            return Err(BridgeError::Invalid {
                field: "relative energy tolerance",
            });
        }
        if self.minimum_gap_m >= self.maximum_gap_m
            || !self.gravity_m_s2.is_finite()
            || self.gravity_m_s2 < 0.0
        {
            return Err(BridgeError::Invalid {
                field: "gap range or gravity",
            });
        }
        for value in [
            self.maximum_bond_number,
            self.maximum_capillary_number,
            self.maximum_weber_number,
        ] {
            positive(value, "quasistatic screen")?;
            limit(value, 0.01, "quasistatic screen ceiling")?;
        }
        Ok(())
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, BridgeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BridgeError::Invalid { field })
    }
}
fn positive(value: f64, field: &'static str) -> Result<f64, BridgeError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(BridgeError::Invalid { field })
    }
}
fn limit(value: f64, maximum: f64, field: &'static str) -> Result<(), BridgeError> {
    if !value.is_finite() || value > maximum {
        Err(BridgeError::Limit {
            field,
            value,
            limit: maximum,
        })
    } else {
        Ok(())
    }
}
fn poll(gate: &CancelGate) -> Result<(), BridgeError> {
    if gate.is_requested() {
        Err(BridgeError::Cancelled)
    } else {
        Ok(())
    }
}
