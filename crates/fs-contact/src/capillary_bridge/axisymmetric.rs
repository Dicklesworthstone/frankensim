//! Actual three-dimensional, axisymmetric capillary adhesion between plates.
//!
//! Each mechanical continuation has fixed liquid volume; an explicit source
//! endpoint can change it under a closed-mass constraint. The two mobile contact
//! lines obey the original (possibly different) acute Young angles.
//! Gravity-free Young--Laplace profiles
//! have one neck and nonpositive liquid-minus-gas pressure. With neck radius a,
//! q=-Delta_p*a/gamma >= 0, and tangent angle beta measured from the axis,
//! r cos(beta)+q r²/(2a)=a(1+q/2). This first integral reduces the shape/volume
//! boundary-value problem to smooth quadrature and a bracketed scalar solve.
//! Pressure and the line resultant give the same force at either plate or neck.
//! No toroidal approximation, sphere fit or two-dimensional extrusion is used.
//!
//! The caller selects a finite q interval. Admission additionally requires a
//! sufficient local stability bound for AXISYMMETRIC, first-order volume-
//! preserving radial graph perturbations with mobile contact lines and fixed
//! plates. For a perturbation u, the Young boundary terms cancel and
//! delta²(E-pV)/(2 pi gamma)=integral[r cos³(beta) u'²-cos(beta)/r u²] dz.
//! The volume constraint integral(r u dz)=0 supplies a zero of u; a one-sided
//! Poincare inequality gives positivity when
//! (pi a/(2h))² min(sin(theta_i))³ > 1. This is a conservative sufficient
//! criterion, not the exact stability boundary. Lateral translations are neutral;
//! no full three-dimensional, nonlinear or dynamic stability is certified.
//!
//! Continuation is quasistatic and conservative, without a time/rate model.
//! It returns the surface-energy work and its work-conjugate mean force. Tiny
//! steps integrate the same force directly when endpoint-energy subtraction
//! loses accuracy; the stored-energy roundoff residual remains explicit.
//! Numerical quadrature/work errors are estimates, not interval certificates.
//! Refusal does not mean physical rupture. Positive-pressure branches, obtuse
//! angles, pinning, finite gravity, cavitation, disjoining pressure, evaporation,
//! inertia and pinch-off need different owners. The caller supplies continuum
//! gap/suction limits and must establish quasistatic applicability in mechanics.
//!
//! First integral/force: Fortes, J. Colloid Interface Sci. (1982),
//! https://doi.org/10.1016/0021-9797(82)90263-6; Millet et al., IJNAMG (2026),
//! https://doi.org/10.1002/nag.70363. Sufficient energy bounding follows the
//! variational approach of Bostwick & Steen, IMA J. Appl. Math. 80 (2015),
//! https://doi.org/10.1093/imamat/hxv016; the explicit conservative bound above
//! is derived here for radial graphs with the stated boundary/volume constraints.

use super::{BridgeError, BridgeVolumeConstraint, finite, limit, poll, positive, refresh_sources};
use crate::normal_patch::PointNormalPort;
use core::f64::consts::{FRAC_PI_2, PI};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_matdb::{InterfaceSystemCard, MaterialCard, QueryPoint};
use fs_material::liquid::surface::adhesion::YoungDupreLaw;
use fs_material::liquid::surface::{ResolvedLiquidSurfaceState, ResolvedLiquidWettingState};
use fs_math::det;

mod numerics;

/// Work-conjugate coupling to the original mechanical held-force compliance.
pub mod load_line;

/// Finite-rate circular contact-line motion on the original Young--Laplace shape.
pub mod contact_line;

/// Independently moving circular contact lines on unlike original surfaces.
pub mod asymmetric_contact_line;

/// Explicit suction branch; q=-Delta_p*a/gamma, a the actual neck radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricBranch {
    pub minimum_suction_parameter: f64,
    pub maximum_suction_parameter: f64,
    /// Positive margin above the conservative axisymmetric stability criterion.
    pub minimum_stability_margin: f64,
}

/// Physical continuum and load envelope. No time/rate assumption is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricLimits {
    pub minimum_gap_m: f64,
    pub maximum_gap_m: f64,
    pub maximum_wetted_radius_m: f64,
    /// Caller must place this below cavitation or other missing physics.
    pub maximum_suction_pa: f64,
    pub maximum_abs_force_n: f64,
    pub maximum_abs_surface_energy_j: f64,
    pub maximum_abs_step_work_j: f64,
    pub gravity_m_s2: f64,
    /// Neglected-gravity applicability screen in (0,0.01], not an error bound.
    pub maximum_bond_number: f64,
}

/// Deterministic work/accuracy limits; at most two angular rules per shape call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricNumerics {
    pub quadrature_relative_tolerance: f64,
    pub gap_relative_tolerance: f64,
    /// Power of two, at least eight; dyadic Simpson reuses all existing samples.
    pub maximum_quadrature_panels: u32,
    pub maximum_root_iterations: u32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricWorkBudget {
    /// Also bounds each dyadic force-work rule by `maximum_quadrature_panels`.
    pub equilibrium: AxisymmetricNumerics,
    pub absolute_tolerance_j: f64,
    pub relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AxisymmetricError {
    Physical(BridgeError),
    QuadratureBudget,
    RootBudget,
    /// Direct finite-work quadrature exhausted its declared panel ceiling.
    WorkQuadratureBudget,
    /// Endpoint of the declared q interval fails the sufficient criterion.
    StabilityUnestablished {
        margin: f64,
        required: f64,
    },
    TurningBranch,
    GapOutsideBranch {
        minimum_m: f64,
        maximum_m: f64,
    },
    StateMismatch,
    /// The finite gap or force-work integral cannot support the requested accuracy.
    WorkUnresolved {
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl From<BridgeError> for AxisymmetricError {
    fn from(error: BridgeError) -> Self {
        Self::Physical(error)
    }
}
impl core::fmt::Display for AxisymmetricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "axisymmetric capillary equilibrium refused: {self:?}")
    }
}
impl std::error::Error for AxisymmetricError {}

/// One actual circular-symmetric bridge, with signed forces along increasing gap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricProfile {
    /// Converged physical gap, within the requested numerical tolerance.
    pub gap_m: f64,
    pub liquid_volume_m3: f64,
    pub neck_radius_m: f64,
    pub neck_height_m: f64,
    pub contact_radii_m: [f64; 2],
    pub contact_angles_radians: [f64; 2],
    pub liquid_gas_area_m2: f64,
    pub wetted_areas_m2: [f64; 2],
    pub liquid_minus_gas_pressure_pa: f64,
    pub normal_force_n: f64,
    pub force_gap_derivative_n_per_m: f64,
    pub contact_radius_gap_derivatives: [f64; 2],
    pub surface_free_energy_j: f64,
    /// Sufficient positivity margin for the perturbation family in module docs.
    pub axisymmetric_stability_margin: f64,
}
impl AxisymmetricProfile {
    /// Endpoint force port; a finite continuation uses its work-conjugate mean.
    /// Native normal-port tangent is with respect to approach=-gap.
    pub fn point_port(
        &self,
        normal: [f64; 3],
        arm_m: [f64; 3],
    ) -> Result<PointNormalPort, AxisymmetricError> {
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !length.is_finite()
            || (length - 1.0).abs() > 1e-12
            || arm_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("axisymmetric port frame"));
        }
        let force = normal.map(|n| n * self.normal_force_n);
        let moment = [
            arm_m[1] * force[2] - arm_m[2] * force[1],
            arm_m[2] * force[0] - arm_m[0] * force[2],
            arm_m[0] * force[1] - arm_m[1] * force[0],
        ];
        for value in force.iter().chain(moment.iter()) {
            finite(*value, "axisymmetric wrench")?;
        }
        Ok(PointNormalPort {
            action_force_n: force,
            action_moment_n_m: moment,
            reaction_force_n: force.map(|v| -v),
            reaction_moment_n_m: moment.map(|v| -v),
            residual_force_n: [0.0; 3],
            tangent_n_per_m: -self.force_gap_derivative_n_per_m,
            dissipated_power_w: 0.0,
            irreversible_work_j: 0.0,
        })
    }
}

/// Refinement estimates and actual solver work, not a certified enclosure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricSolveReport {
    pub gap_residual_m: f64,
    pub quadrature_relative_error_estimate: f64,
    pub surface_energy_error_estimate_j: f64,
    pub root_iterations: u32,
    pub quadrature_evaluations: u64,
}

/// Immutable accepted equilibrium; full source owners are never copied per solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricState {
    identity: ContentHash,
    q: f64,
    profile: AxisymmetricProfile,
    report: AxisymmetricSolveReport,
}
impl AxisymmetricState {
    #[must_use]
    pub const fn profile(&self) -> &AxisymmetricProfile {
        &self.profile
    }
    #[must_use]
    pub const fn report(&self) -> &AxisymmetricSolveReport {
        &self.report
    }
    #[must_use]
    pub const fn suction_parameter(&self) -> f64 {
        self.q
    }
}

/// Conservative proposal. The consumer publishes the successor with mechanics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricStep {
    state: AxisymmetricState,
    pub requested_gap_m: f64,
    /// Positive work is supplied to the liquid surface free-energy store.
    pub work_on_liquid_j: f64,
    /// Force-work quadrature/roundoff uncertainty; independent of subtraction
    /// of the two much larger endpoint surface-energy stores.
    pub estimated_work_error_j: f64,
    /// Delta(stored surface energy) - work. This is numerical roundoff/error,
    /// never physical dissipation or an unreported thermal exchange.
    pub surface_energy_balance_residual_j: f64,
    /// -work/Delta(actual converged gap), never an endpoint-force shortcut.
    pub mean_normal_force_n: f64,
    pub mean_force_gap_derivative_n_per_m: f64,
}
impl AxisymmetricStep {
    #[must_use]
    pub const fn state(&self) -> AxisymmetricState {
        self.state
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeridianPoint {
    pub height_m: f64,
    pub radius_m: f64,
    pub radius_axial_derivative: f64,
    /// Signed principal curvatures sum to Delta_p/gamma.
    pub meridional_curvature_per_m: f64,
    pub azimuthal_curvature_per_m: f64,
}

#[derive(Clone, Debug)]
pub struct SourceAxisymmetricCapillaryBridge {
    liquid: ResolvedLiquidSurfaceState,
    wetting: [ResolvedLiquidWettingState; 2],
    laws: [YoungDupreLaw; 2],
    volume: f64,
    sigma: f64,
    alpha: [f64; 2],
    cosine: [f64; 2],
    minimum_sine: f64,
    branch: AxisymmetricBranch,
    limits: AxisymmetricLimits,
    identity: ContentHash,
}

/// Original-source material endpoint on the same admitted Young--Laplace branch.
/// Publish this candidate only with the receiving mechanics' acceptance.
#[derive(Clone, Debug)]
pub struct AxisymmetricBridgeSourceRefresh {
    source: SourceAxisymmetricCapillaryBridge,
    previous: AxisymmetricState,
    state: AxisymmetricState,
    gap_m: f64,
    parameter_work_j: f64,
    estimated_parameter_work_error_j: f64,
}
impl AxisymmetricBridgeSourceRefresh {
    #[must_use]
    pub const fn source(&self) -> &SourceAxisymmetricCapillaryBridge {
        &self.source
    }
    #[must_use]
    pub const fn previous_state(&self) -> AxisymmetricState {
        self.previous
    }
    #[must_use]
    pub const fn state(&self) -> AxisymmetricState {
        self.state
    }
    /// The predecessor's requested gap, recovered from its explicit root residual.
    #[must_use]
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    /// Signed change of surface free energy, separate from plate work and heat.
    #[must_use]
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }
    /// Original quadrature and gap-root work error estimates, not an enclosure.
    #[must_use]
    pub const fn estimated_parameter_work_error_j(&self) -> f64 {
        self.estimated_parameter_work_error_j
    }
}
impl SourceAxisymmetricCapillaryBridge {
    pub fn new(
        liquid: ResolvedLiquidSurfaceState,
        wetting: [ResolvedLiquidWettingState; 2],
        volume_m3: f64,
        branch: AxisymmetricBranch,
        limits: AxisymmetricLimits,
    ) -> Result<Self, AxisymmetricError> {
        positive(volume_m3, "axisymmetric volume")?;
        positive(
            liquid.density_kg_m3() * volume_m3,
            "axisymmetric retained mass",
        )?;
        limits.validate()?;
        if !branch.minimum_suction_parameter.is_finite()
            || branch.minimum_suction_parameter < 0.0
            || !branch.maximum_suction_parameter.is_finite()
            || branch.maximum_suction_parameter <= branch.minimum_suction_parameter
            || !(branch.maximum_suction_parameter * (2.0 + branch.maximum_suction_parameter))
                .is_finite()
        {
            return Err(invalid("finite ordered suction branch"));
        }
        positive(
            branch.minimum_stability_margin,
            "axisymmetric stability margin",
        )?;
        if liquid.rheology_axes().is_some() {
            return Err(invalid("frozen liquid source"));
        }
        for boundary in &wetting {
            boundary
                .require_liquid_state(&liquid)
                .map_err(BridgeError::from)?;
        }
        // The liquid and surrounding environment are common; each original
        // solid surface may have its own preparation/contamination history.
        // Those histories remain in the ordered source identities below.
        let first = wetting[0].context();
        let second = wetting[1].context();
        if first.medium != second.medium
            || first.third_body != second.third_body
            || first.environment != second.environment
        {
            return Err(invalid("common original surface environment"));
        }
        let sigma = liquid
            .surface_tension_n_m()
            .ok_or_else(|| invalid("active sourced tension"))?;
        let mut laws = [YoungDupreLaw::new(sigma, sigma).map_err(|_| invalid("surface law"))?; 2];
        let mut alpha = [0.0; 2];
        let mut cosine = [0.0; 2];
        let mut minimum_sine: f64 = 1.0;
        for i in 0..2 {
            // An original direct angle selects the ideal Young relation; the
            // resulting W is numerical, never fabricated as a material claim.
            laws[i] = wetting[i]
                .young_dupre_law()
                .map(Ok)
                .unwrap_or_else(|| {
                    YoungDupreLaw::new(
                        sigma * (1.0 + det::cos(wetting[i].contact_angle_radians())),
                        sigma,
                    )
                })
                .map_err(|_| invalid("representable ideal partial wetting"))?;
            cosine[i] = (laws[i].work_of_adhesion_j_m2() - sigma) / sigma;
            if !(cosine[i] > 0.0 && cosine[i] < 1.0) {
                return Err(invalid("strictly acute sourced Young angles"));
            }
            let sine = ((1.0 - cosine[i]) * (1.0 + cosine[i])).sqrt();
            minimum_sine = minimum_sine.min(sine);
            alpha[i] = det::atan2(cosine[i], sine);
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-contact.axisymmetric-capillary-bridge.v1");
        hash.update(liquid.material().identity().as_bytes());
        for w in &wetting {
            hash.update(w.resolved().identity().as_bytes());
        }
        for value in [
            volume_m3,
            laws[0].work_of_adhesion_j_m2(),
            laws[1].work_of_adhesion_j_m2(),
            branch.minimum_suction_parameter,
            branch.maximum_suction_parameter,
            branch.minimum_stability_margin,
            limits.minimum_gap_m,
            limits.maximum_gap_m,
            limits.maximum_wetted_radius_m,
            limits.maximum_suction_pa,
            limits.maximum_abs_force_n,
            limits.maximum_abs_surface_energy_j,
            limits.maximum_abs_step_work_j,
            limits.gravity_m_s2,
            limits.maximum_bond_number,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            liquid,
            wetting,
            laws,
            volume: volume_m3,
            sigma,
            alpha,
            cosine,
            minimum_sine,
            branch,
            limits,
            identity: hash.finalize(),
        })
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
    pub const fn surface_laws(&self) -> &[YoungDupreLaw; 2] {
        &self.laws
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Continue the same original liquid and ordered surfaces at a new endpoint.
    ///
    /// `FixedMass` resolves thermal expansion from the refreshed density; no
    /// material, source claim, phase, branch, plate gap or physical limit is
    /// substituted. The existing shape solver must establish the same selected
    /// suction/stability branch at the predecessor's requested physical gap.
    /// Quadrature and root-gap contributions must fit the supplied work budget.
    ///
    /// The signed surface parameter work is not deposited heat, bulk pressure
    /// work, latent heat or a complete thermodynamic balance. This operation
    /// establishes endpoint support only, not a transient thermal/meniscus path.
    pub fn refresh_at_query_point(
        &self,
        previous: AxisymmetricState,
        material: &MaterialCard,
        interfaces: [&InterfaceSystemCard; 2],
        point: &QueryPoint,
        constraint: BridgeVolumeConstraint,
        budget: AxisymmetricWorkBudget,
        gate: &CancelGate,
    ) -> Result<AxisymmetricBridgeSourceRefresh, AxisymmetricError> {
        poll(gate)?;
        let before = self.evaluate(previous)?;
        budget.equilibrium.validate()?;
        positive(
            budget.absolute_tolerance_j,
            "source parameter-work absolute tolerance",
        )?;
        if !budget.relative_tolerance.is_finite()
            || budget.relative_tolerance < 0.0
            || budget.relative_tolerance >= 1.0
        {
            return Err(invalid("source parameter-work relative tolerance"));
        }
        let (liquid, wetting, volume) = refresh_sources(
            &self.liquid,
            &self.wetting,
            self.volume,
            material,
            interfaces,
            point,
            constraint,
        )?;
        let source = Self::new(liquid, wetting, volume, self.branch, self.limits)?;
        let gap_m = positive(
            before.gap_m - previous.report.gap_residual_m,
            "retained physical source-refresh gap",
        )?;
        let state = source.solve_gap(gap_m, Some(previous.q), budget.equilibrium, gate)?;
        let after = state.profile;
        let parameter_work_j = finite(
            after.surface_free_energy_j - before.surface_free_energy_j,
            "axisymmetric surface parameter work",
        )?;
        let error = finite(
            previous.report.surface_energy_error_estimate_j
                + state.report.surface_energy_error_estimate_j
                + before.normal_force_n.abs() * previous.report.gap_residual_m.abs()
                + after.normal_force_n.abs() * state.report.gap_residual_m.abs(),
            "axisymmetric source parameter-work error estimate",
        )?;
        let tolerance = finite(
            budget.absolute_tolerance_j + budget.relative_tolerance * parameter_work_j.abs(),
            "axisymmetric source parameter-work tolerance",
        )?;
        if error > tolerance {
            return Err(AxisymmetricError::WorkUnresolved {
                estimated_error_j: error,
                tolerance_j: tolerance,
            });
        }
        limit(
            parameter_work_j.abs(),
            self.limits.maximum_abs_step_work_j,
            "axisymmetric source parameter work",
        )?;
        poll(gate)?;
        Ok(AxisymmetricBridgeSourceRefresh {
            source,
            previous,
            state,
            gap_m,
            parameter_work_j,
            estimated_parameter_work_error_j: error,
        })
    }

    /// Authenticate and read an accepted equilibrium without solving it again.
    pub fn evaluate(
        &self,
        state: AxisymmetricState,
    ) -> Result<AxisymmetricProfile, AxisymmetricError> {
        self.require_state(state)?;
        Ok(state.profile)
    }

    /// Direct branch-parameter initialization, also useful for continuation decks.
    pub fn solve_at_suction(
        &self,
        q: f64,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<AxisymmetricState, AxisymmetricError> {
        budget.validate()?;
        if !q.is_finite()
            || q < self.branch.minimum_suction_parameter
            || q > self.branch.maximum_suction_parameter
        {
            return Err(invalid("q within selected suction branch"));
        }
        let shape = self.shape(q, budget, gate)?;
        self.accept(shape, 0.0, 0)
    }

    /// Solve the fixed-volume shape at the requested plate spacing. Actual gap
    /// and residual are explicit; callers use the converged gap in path work.
    pub fn solve_at_gap(
        &self,
        gap_m: f64,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<AxisymmetricState, AxisymmetricError> {
        self.solve_gap(gap_m, None, budget, gate)
    }

    /// Quasistatic continuation on the same original-source branch. No topology
    /// release or heat is invented. A cancellation-prone energy secant uses
    /// bounded direct force quadrature on the same branch. A geometrically
    /// unresolved increment still refuses. This is not a
    /// finite-speed contact-line or liquid-momentum integrator.
    pub fn continue_to_gap(
        &self,
        previous: AxisymmetricState,
        gap_m: f64,
        budget: AxisymmetricWorkBudget,
        gate: &CancelGate,
    ) -> Result<AxisymmetricStep, AxisymmetricError> {
        poll(gate)?;
        self.require_state(previous)?;
        budget.equilibrium.validate()?;
        positive(
            budget.absolute_tolerance_j,
            "axisymmetric work absolute tolerance",
        )?;
        if !budget.relative_tolerance.is_finite()
            || budget.relative_tolerance < 0.0
            || budget.relative_tolerance >= 1.0
        {
            return Err(invalid("axisymmetric work relative tolerance"));
        }
        if gap_m == previous.profile.gap_m {
            return Ok(AxisymmetricStep {
                state: previous,
                requested_gap_m: gap_m,
                work_on_liquid_j: 0.0,
                estimated_work_error_j: 0.0,
                surface_energy_balance_residual_j: 0.0,
                mean_normal_force_n: previous.profile.normal_force_n,
                mean_force_gap_derivative_n_per_m: 0.5
                    * previous.profile.force_gap_derivative_n_per_m,
            });
        }
        let state = self.solve_gap(gap_m, Some(previous.q), budget.equilibrium, gate)?;
        let dh = state.profile.gap_m - previous.profile.gap_m;
        let energy_difference = finite(
            state.profile.surface_free_energy_j - previous.profile.surface_free_energy_j,
            "axisymmetric branch work",
        )?;
        let subtraction_error = finite(
            previous.report.surface_energy_error_estimate_j
                + state.report.surface_energy_error_estimate_j,
            "axisymmetric work error estimate",
        )?;
        // A few ulps of gap are not a resolved mechanical trajectory, even if
        // the absolute energy allowance would admit it. Keep rest distinct.
        if dh.abs() <= 64.0 * f64::EPSILON * gap_m.abs().max(previous.profile.gap_m.abs()) {
            return Err(AxisymmetricError::WorkUnresolved {
                estimated_error_j: subtraction_error,
                tolerance_j: budget.absolute_tolerance_j,
            });
        }
        // An absolute energy tolerance can hide an inaccurate force after
        // division by a very small dh. Retain the inexpensive energy secant
        // only when its relative accuracy also resolves that mean force.
        let relative = budget.relative_tolerance.max(32.0 * f64::EPSILON);
        let (work, error, mean, tangent) = if energy_difference.signum() == dh.signum()
            && subtraction_error <= relative * energy_difference.abs()
        {
            let mean = finite(-energy_difference / dh, "axisymmetric mean force")?;
            (
                energy_difference,
                subtraction_error,
                mean,
                finite(
                    (state.profile.normal_force_n - mean) / dh,
                    "axisymmetric mean force derivative",
                )?,
            )
        } else {
            self.integrate_force_work(previous, state, budget, gate)?
        };
        let residual = finite(
            energy_difference - work,
            "axisymmetric surface-energy residual",
        )?;
        let tolerance = finite(
            budget.absolute_tolerance_j + budget.relative_tolerance * work.abs(),
            "axisymmetric work tolerance",
        )?;
        if error > tolerance
            || work == 0.0
            || work.signum() != dh.signum()
            || error >= work.abs()
            || residual.abs() > subtraction_error + error
        {
            return Err(AxisymmetricError::WorkUnresolved {
                estimated_error_j: error.max(residual.abs() - subtraction_error),
                tolerance_j: tolerance,
            });
        }
        limit(
            work.abs(),
            self.limits.maximum_abs_step_work_j,
            "axisymmetric step work",
        )?;
        limit(
            mean.abs(),
            self.limits.maximum_abs_force_n,
            "axisymmetric mean force",
        )?;
        poll(gate)?;
        Ok(AxisymmetricStep {
            state,
            requested_gap_m: gap_m,
            work_on_liquid_j: work,
            estimated_work_error_j: error,
            surface_energy_balance_residual_j: residual,
            mean_normal_force_n: mean,
            mean_force_gap_derivative_n_per_m: tangent,
        })
    }

    /// Sample the true meridian without allocating. Parameter 0 is the lower
    /// plate and 1 the upper; it varies the tangent angle, not the axial height.
    pub fn meridian_point(
        &self,
        state: AxisymmetricState,
        fraction: f64,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<MeridianPoint, AxisymmetricError> {
        self.require_state(state)?;
        budget.validate()?;
        poll(gate)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(invalid("meridian parameter"));
        }
        let beta = -self.alpha[0] + fraction * (self.alpha[0] + self.alpha[1]);
        let angle = beta.abs();
        let (x, _, _) = numerics::radius(angle, state.q);
        let height = if fraction == 0.0 {
            0.0
        } else if fraction == 1.0 {
            state.profile.gap_m
        } else if angle == 0.0 {
            state.profile.neck_height_m
        } else {
            state.profile.neck_height_m
                + beta.signum()
                    * state.profile.neck_radius_m
                    * numerics::integrate(angle, state.q, budget, gate)?.values[0]
        };
        let radius = positive(
            state.profile.neck_radius_m * x,
            "axisymmetric sampled radius",
        )?;
        let (s, c) = det::sin_cos(beta);
        let azimuthal = finite(c / radius, "azimuthal curvature")?;
        Ok(MeridianPoint {
            height_m: finite(height, "meridian height")?,
            radius_m: radius,
            radius_axial_derivative: finite(s / c, "meridian slope")?,
            azimuthal_curvature_per_m: azimuthal,
            meridional_curvature_per_m: finite(
                -state.q / state.profile.neck_radius_m - azimuthal,
                "meridional curvature",
            )?,
        })
    }

    fn require_state(&self, state: AxisymmetricState) -> Result<(), AxisymmetricError> {
        if state.identity != self.identity {
            Err(AxisymmetricError::StateMismatch)
        } else {
            Ok(())
        }
    }
    fn check_gap(&self, gap: f64) -> Result<(), AxisymmetricError> {
        positive(gap, "axisymmetric gap")?;
        if gap < self.limits.minimum_gap_m {
            return Err(invalid("axisymmetric minimum continuum gap"));
        }
        limit(gap, self.limits.maximum_gap_m, "axisymmetric gap")?;
        Ok(())
    }
    fn accept(
        &self,
        shape: numerics::Shape,
        gap_residual: f64,
        iterations: u32,
    ) -> Result<AxisymmetricState, AxisymmetricError> {
        let p = shape.profile;
        self.check_gap(p.gap_m)?;
        let radius = p.contact_radii_m[0].max(p.contact_radii_m[1]);
        limit(
            radius,
            self.limits.maximum_wetted_radius_m,
            "axisymmetric wetted radius",
        )?;
        limit(
            -p.liquid_minus_gas_pressure_pa,
            self.limits.maximum_suction_pa,
            "axisymmetric suction",
        )?;
        limit(
            p.normal_force_n.abs(),
            self.limits.maximum_abs_force_n,
            "axisymmetric force",
        )?;
        limit(
            p.surface_free_energy_j.abs(),
            self.limits.maximum_abs_surface_energy_j,
            "axisymmetric surface energy",
        )?;
        let length = p.gap_m.max(2.0 * radius);
        let bond = finite(
            self.liquid.density_kg_m3() * self.limits.gravity_m_s2 / self.sigma * length * length,
            "axisymmetric Bond number",
        )?;
        limit(
            bond,
            self.limits.maximum_bond_number,
            "axisymmetric Bond number",
        )?;
        Ok(AxisymmetricState {
            identity: self.identity,
            q: shape.q,
            profile: p,
            report: AxisymmetricSolveReport {
                gap_residual_m: gap_residual,
                quadrature_relative_error_estimate: shape.relative_error,
                surface_energy_error_estimate_j: shape.energy_error,
                root_iterations: iterations,
                quadrature_evaluations: shape.evaluations,
            },
        })
    }
}

fn invalid(field: &'static str) -> AxisymmetricError {
    BridgeError::Invalid { field }.into()
}
impl AxisymmetricLimits {
    fn validate(self) -> Result<(), AxisymmetricError> {
        for (v, n) in [
            (self.minimum_gap_m, "minimum gap"),
            (self.maximum_gap_m, "maximum gap"),
            (self.maximum_wetted_radius_m, "maximum wetted radius"),
            (self.maximum_suction_pa, "maximum suction"),
            (self.maximum_abs_force_n, "maximum force"),
            (self.maximum_abs_surface_energy_j, "maximum energy"),
            (self.maximum_abs_step_work_j, "maximum work"),
            (self.maximum_bond_number, "maximum Bond number"),
        ] {
            positive(v, n)?;
        }
        if self.maximum_gap_m <= self.minimum_gap_m
            || !self.gravity_m_s2.is_finite()
            || self.gravity_m_s2 < 0.0
            || self.maximum_bond_number > 0.01
        {
            return Err(invalid("axisymmetric physical limits"));
        }
        Ok(())
    }
}
impl AxisymmetricNumerics {
    fn validate(self) -> Result<(), AxisymmetricError> {
        if !self.quadrature_relative_tolerance.is_finite()
            || !self.gap_relative_tolerance.is_finite()
            || self.quadrature_relative_tolerance < 64.0 * f64::EPSILON
            || self.quadrature_relative_tolerance >= 0.01
            || self.gap_relative_tolerance < 64.0 * f64::EPSILON
            || self.gap_relative_tolerance >= 0.01
            || self.maximum_quadrature_panels < 8
            || !self.maximum_quadrature_panels.is_power_of_two()
            || self.maximum_root_iterations == 0
        {
            return Err(invalid("axisymmetric numerical budget"));
        }
        Ok(())
    }
}
