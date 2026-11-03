//! Steady annular mass/angular-momentum transport to downstream lifting surfaces.
//!
//! This explicitly selected engineering model transports each accepted native
//! blade torque Q_i through a fixed cylindrical annulus: L_i=Q_i/(rho U A_i),
//! v_theta=L_i/r. Thus integral(rho U r v_theta dA)=Q_i without inserting or
//! removing a Prandtl factor a second time. Only the no-loss, matching-flow,
//! narrow-annulus limit gives v_theta=2u; u is the native AT-DISK induction.
//! See the angular balance in NREL/TP-500-36881, equations (9),(10),
//! <https://www.nlr.gov/docs/libraries/wind-docs/ad-theory.pdf>.
//!
//! U is the existing axial wash model at that downstream section. Fixed area
//! and varying U imply varying mass flux: no contraction/entrainment solution,
//! radial pressure balance, resolved helical wake or full fluid energy theorem
//! is claimed. These are unobstructed source fields: momentum transferred to a
//! wing is not advected into a reconstructed post-wing exhaust profile.
//! The analytic resolved kinetic-power defect is bounded against
//! the native swirl power by caller-admitted tolerances; it is NOT heat. With
//! Prandtl losses, unresolved azimuthal structure can carry some of that power.
//!
//! The observation window is strictly downstream and shorter than 2R. Inside
//! the cylinder, samples on the disk or beyond that window REFUSE, rather than
//! deleting angular momentum at a cutoff. Upstream/lateral/hub samples are zero.
//! Receiving disks intersecting another wake refuse: the native axial rotor
//! API cannot consume the incident angular profile or predict swirl recovery.
//!
//! The opt-in MaterialStreamtubes mode replaces downstream axial and angular
//! wash with per-native-annulus mass/torque transport and continuity-driven
//! radial velocity. It retains the upstream axial model. See streamtube for
//! its explicit fixed hub, pressure/exterior-flow and source-field limits.

use super::{
    Controls, CoupledSwirlStep, DiskAxes, IncidentInflow, InfluenceOperator, Panel,
    PropCouplingSolverSpec, PropDisk, Refusal, Rotor, SectionClosure, SectionFrame, Solution,
    StripSpec, WarmStart, coupled_impl,
};
use fs_blake3::hash_domain;
use fs_ga::Vec3;
use fs_math::det;

#[path = "wake/streamtube.rs"]
mod streamtube;
pub use streamtube::MaterialAnnulus;

/// Explicit downstream source-field closure. Historical presets keep their
/// original axial model; even this angular-wake API is opt-in.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WakeTransport {
    #[default]
    FixedCylinder,
    /// Preserve each native annulus's disk mass flux and signed torque while
    /// resolving axial development and radial contraction/expansion.
    MaterialStreamtubes,
}

fn vec(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn reject(code: &'static str, detail: &str) -> Refusal {
    Refusal {
        code,
        message: detail.into(),
        ranked_repairs: vec![
            "use separated disks and downstream surfaces inside the declared wake window; inspect the physical energy budget".into()
        ],
    }
}

/// Bounded field/wing/rotor fixed point and explicit projection-error budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeControls {
    /// Complete nested wing/rotor evaluations, in 1..=128.
    pub maximum_field_evaluations: u32,
    /// Fixed field relaxation, in (0,1].
    pub relaxation: f64,
    /// Maximum Euclidean change at ANY sampled strip [m/s].
    pub velocity_tolerance_mps: f64,
    /// Absolute allowance on sum_i |P_resolved_i-P_native_i| [W].
    pub absolute_power_tolerance_w: f64,
    /// Additional allowance relative to native total swirl power, in [0,1].
    /// This is a model-error admission, not a numerical convergence tolerance.
    pub relative_power_tolerance: f64,
}
impl Default for WakeControls {
    fn default() -> Self {
        Self {
            maximum_field_evaluations: 48,
            relaxation: 0.5,
            velocity_tolerance_mps: 1e-7,
            absolute_power_tolerance_w: 0.0,
            relative_power_tolerance: 0.5,
        }
    }
}
impl WakeControls {
    fn admit(self) -> Result<(), Refusal> {
        if !(1..=128).contains(&self.maximum_field_evaluations)
            || ![
                self.relaxation,
                self.velocity_tolerance_mps,
                self.absolute_power_tolerance_w,
                self.relative_power_tolerance,
            ]
            .into_iter()
            .all(f64::is_finite)
            || self.relaxation <= 0.0
            || self.relaxation > 1.0
            || self.velocity_tolerance_mps <= 0.0
            || self.absolute_power_tolerance_w < 0.0
            || !(0.0..=1.0).contains(&self.relative_power_tolerance)
        {
            return Err(reject(
                "swirl-wake-invalid",
                "invalid field iteration or power-error controls",
            ));
        }
        Ok(())
    }
}

/// Original radial quadrature cell, with no post-hoc rescaling of blade torque.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Annulus {
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    /// Signed blade torque in the positive blade-spin convention [N m].
    pub torque_nm: f64,
    /// Native (u/r) Q for this cell [W].
    pub native_swirl_power_w: f64,
}
impl Annulus {
    fn area(self) -> f64 {
        core::f64::consts::PI
            * (self.outer_radius_m - self.inner_radius_m)
            * (self.outer_radius_m + self.inner_radius_m)
    }
    fn resolved_power(self, density: f64, axial: f64) -> f64 {
        let specific = self.torque_nm / (density * axial * self.area());
        core::f64::consts::PI
            * density
            * axial
            * specific
            * specific
            * det::ln(self.outer_radius_m / self.inner_radius_m)
    }
}

/// Aggregate fluxes of the explicitly selected annular source-field closure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrossSection {
    /// Uniform speed in FixedCylinder; area-weighted mean in MaterialStreamtubes.
    /// Use material_annuli for the actual radial speed distribution.
    pub axial_speed_mps: f64,
    /// Total annular mass flux [kg/s]. Constant in MaterialStreamtubes; varying
    /// with distance is the disclosed FixedCylinder continuity defect.
    pub mass_flux_kg_s: f64,
    /// Signed angular momentum per second about the positive thrust axis.
    pub angular_momentum_flux_nm: f64,
    pub resolved_swirl_power_w: f64,
    pub native_swirl_power_w: f64,
    /// Resolved minus native; neither erased nor classified as heat [W].
    pub power_defect_w: f64,
    /// Sum of absolute cell defects; opposite cell errors cannot cancel [W].
    pub absolute_cell_power_defect_w: f64,
}

/// Immutable admitted angular field from one actual native rotor solution.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnularWake {
    center_m: [f64; 3],
    axis: [f64; 3],
    spin_sign: f64,
    radius_m: f64,
    density: f64,
    incident_axial_mps: f64,
    disk_induction_mps: f64,
    cells: Vec<Annulus>,
    maximum_power_defect_w: f64,
    streamtubes: Option<streamtube::Streamtubes>,
}
impl AnnularWake {
    #[allow(clippy::too_many_arguments)]
    fn from_solution(
        rotor: &Rotor,
        solution: &Solution,
        center_m: [f64; 3],
        axis: [f64; 3],
        spin_sign: f64,
        density: f64,
        axial: f64,
        induction: f64,
        controls: WakeControls,
        transport: WakeTransport,
    ) -> Result<Self, Refusal> {
        let mut cells = Vec::with_capacity(rotor.stations.len());
        for (i, station) in solution.stations.iter().enumerate() {
            let low = if i == 0 {
                0.8 * station.r_over_r
            } else {
                0.5 * (rotor.stations[i - 1].r_over_r + station.r_over_r)
            };
            let high = if i + 1 == rotor.stations.len() {
                1.0
            } else {
                0.5 * (station.r_over_r + rotor.stations[i + 1].r_over_r)
            };
            cells.push(Annulus {
                inner_radius_m: low * rotor.radius_m,
                outer_radius_m: high * rotor.radius_m,
                torque_nm: station.torque_nm,
                native_swirl_power_w: station.induction.tangential_mps * station.torque_nm
                    / (station.r_over_r * rotor.radius_m),
            });
        }
        let mut wake = Self {
            center_m,
            axis,
            spin_sign,
            radius_m: rotor.radius_m,
            density,
            incident_axial_mps: axial,
            disk_induction_mps: induction,
            cells,
            maximum_power_defect_w: 0.0,
            streamtubes: None,
        };
        if transport == WakeTransport::MaterialStreamtubes {
            let tubes = streamtube::Streamtubes::new(&wake, solution, controls)?;
            wake.maximum_power_defect_w = tubes.maximum_power_defect_w();
            wake.streamtubes = Some(tubes);
            return Ok(wake);
        }
        // Each resolved cell power is C_i/U. Sum |C_i/U-P_i| is convex in
        // 1/U, so checking both monotone window endpoints bounds its interior.
        let near = wake.cross_section(0.0)?;
        let far = wake.cross_section(2.0 * rotor.radius_m)?;
        wake.maximum_power_defect_w = near
            .absolute_cell_power_defect_w
            .max(far.absolute_cell_power_defect_w);
        let budget = controls.absolute_power_tolerance_w
            + controls.relative_power_tolerance * near.native_swirl_power_w;
        if !budget.is_finite() || wake.maximum_power_defect_w > budget {
            return Err(reject(
                "swirl-wake-power-budget",
                &format!(
                    "annular projection defect {} W exceeds admitted {} W",
                    wake.maximum_power_defect_w, budget
                ),
            ));
        }
        Ok(wake)
    }

    #[must_use]
    pub fn annuli(&self) -> &[Annulus] {
        &self.cells
    }
    /// The selected immutable source-field physics.
    #[must_use]
    pub fn transport(&self) -> WakeTransport {
        if self.streamtubes.is_some() {
            WakeTransport::MaterialStreamtubes
        } else {
            WakeTransport::FixedCylinder
        }
    }

    /// Actual material annuli, radial interface slopes and kinetic transports.
    /// Original disk annuli remain available through annuli().
    ///
    /// # Errors
    /// Refuses for FixedCylinder, invalid distance or unrepresentable geometry.
    pub fn material_annuli(&self, distance_m: f64) -> Result<Vec<MaterialAnnulus>, Refusal> {
        self.streamtubes
            .as_ref()
            .ok_or_else(|| {
                reject(
                    "swirl-wake-transport-mismatch",
                    "select MaterialStreamtubes",
                )
            })?
            .material_annuli(self, distance_m)
    }

    /// Upper bound over the complete admitted observation window [W].
    #[must_use]
    pub fn maximum_power_defect_w(&self) -> f64 {
        self.maximum_power_defect_w
    }

    /// Fluxes at downstream distance s in [0,2R], including endpoint limits.
    /// These limits do not admit a lifting surface on the singular disk plane.
    ///
    /// # Errors
    /// Nonfinite/out-of-window distance, nonpositive throughflow or overflow.
    pub fn cross_section(&self, distance_m: f64) -> Result<CrossSection, Refusal> {
        if let Some(tubes) = &self.streamtubes {
            return tubes.cross_section(self, distance_m);
        }
        if !distance_m.is_finite() || distance_m < 0.0 || distance_m > 2.0 * self.radius_m {
            return Err(reject(
                "swirl-wake-outside-domain",
                "cross section outside 0..=2R",
            ));
        }
        let axial = self.axial_speed(distance_m)?;
        let mut flux = CrossSection {
            axial_speed_mps: axial,
            mass_flux_kg_s: 0.0,
            angular_momentum_flux_nm: 0.0,
            resolved_swirl_power_w: 0.0,
            native_swirl_power_w: 0.0,
            power_defect_w: 0.0,
            absolute_cell_power_defect_w: 0.0,
        };
        for &cell in &self.cells {
            let power = cell.resolved_power(self.density, axial);
            flux.mass_flux_kg_s += self.density * axial * cell.area();
            flux.angular_momentum_flux_nm += self.spin_sign * cell.torque_nm;
            flux.resolved_swirl_power_w += power;
            flux.native_swirl_power_w += cell.native_swirl_power_w;
            flux.absolute_cell_power_defect_w += (power - cell.native_swirl_power_w).abs();
        }
        flux.power_defect_w = flux.resolved_swirl_power_w - flux.native_swirl_power_w;
        if ![
            flux.mass_flux_kg_s,
            flux.angular_momentum_flux_nm,
            flux.resolved_swirl_power_w,
            flux.native_swirl_power_w,
            flux.absolute_cell_power_defect_w,
        ]
        .into_iter()
        .all(f64::is_finite)
            || flux.native_swirl_power_w < 0.0
        {
            return Err(reject("swirl-wake-invalid", "unrepresentable annular flux"));
        }
        Ok(flux)
    }

    fn axial_speed(&self, distance_m: f64) -> Result<f64, Refusal> {
        let x = distance_m / self.radius_m;
        let factor = 1.0 + x / det::sqrt(1.0 + x * x);
        let axial = self.incident_axial_mps + factor * self.disk_induction_mps;
        if !axial.is_finite() || axial <= 0.0 {
            return Err(reject(
                "swirl-wake-outside-domain",
                "nonpositive downstream throughflow",
            ));
        }
        Ok(axial)
    }

    /// Increment to MATERIAL-RELATIVE wing velocity. FixedCylinder supplies
    /// swirl; MaterialStreamtubes supplies its complete downstream induced
    /// axial, radial and swirl velocity. Radial/swirl signs oppose the fluid.
    /// Positive blade spin sends fluid along axis cross radius. Handedness and
    /// signed windmilling torque are applied once each.
    ///
    /// # Errors
    /// Nonfinite geometry; an in-annulus point on the disk or at/beyond 2R.
    pub fn relative_velocity_at(&self, point_m: [f64; 3]) -> Result<[f64; 3], Refusal> {
        if let Some(tubes) = &self.streamtubes {
            return tubes.relative_velocity_at(self, point_m);
        }
        if !point_m.into_iter().all(f64::is_finite) {
            return Err(reject("swirl-wake-invalid", "nonfinite sample point"));
        }
        let offset = vec(point_m) - vec(self.center_m);
        let axial = offset.dot(vec(self.axis));
        let radial = offset - vec(self.axis).scale(axial);
        let radius = radial.norm();
        if !axial.is_finite() || !radius.is_finite() {
            return Err(reject(
                "swirl-wake-invalid",
                "unrepresentable sample geometry",
            ));
        }
        let Some(cell) = self
            .cells
            .iter()
            .find(|c| radius >= c.inner_radius_m && radius < c.outer_radius_m)
        else {
            return Ok([0.0; 3]);
        };
        if axial > 0.0 {
            return Ok([0.0; 3]);
        }
        let distance = -axial;
        if distance <= 0.0 || distance >= 2.0 * self.radius_m {
            return Err(reject(
                "swirl-wake-outside-domain",
                "in-wake surface outside 0<s<2R",
            ));
        }
        let axial_speed = self.axial_speed(distance)?;
        let specific = cell.torque_nm / (self.density * axial_speed * cell.area());
        let velocity = vec(self.axis)
            .cross(radial)
            .scale(-self.spin_sign * specific / (radius * radius));
        let result = array(velocity);
        if !result.into_iter().all(f64::is_finite) {
            return Err(reject(
                "swirl-wake-invalid",
                "unrepresentable swirl velocity",
            ));
        }
        Ok(result)
    }
}

/// Accepted actual wing/rotor solve plus the angular field used in its forces.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledWakeStep {
    pub solution: CoupledSwirlStep,
    /// Fields from the SAME accepted native rotor states. The reported field
    /// residual compares their desired contribution with the actual accepted
    /// force velocity, including material mode's legacy-wash replacement.
    pub wakes: [AnnularWake; 2],
    /// Exact numerical velocity correction passed to the inner wing solve.
    /// FixedCylinder adds only swirl. MaterialStreamtubes also replaces the
    /// actual downstream legacy wash with axial/radial streamtube transport.
    /// The historical field name is retained; the correction accessor below
    /// names its complete meaning. No second authoritative velocity is stored.
    pub used_strip_swirl_velocity_mps: Vec<[f64; 3]>,
    pub field_evaluations: u32,
    /// Maximum Euclidean used-versus-predicted field discrepancy [m/s].
    pub field_residual_mps: f64,
}

/// Explicit numerical continuation; no extra physical angular momentum state.
#[derive(Clone, Debug, PartialEq)]
pub struct WakeWarmStart {
    pub rotor: WarmStart,
    /// Numerical correction seed; see CoupledWakeStep::strip_velocity_correction_mps.
    pub strip_swirl_velocity_mps: Vec<[f64; 3]>,
}
impl CoupledWakeStep {
    /// The exact all-component correction used by the accepted inner solve.
    /// Complete force velocities remain in solution.coupled.strip_velocity_mps.
    #[must_use]
    pub fn strip_velocity_correction_mps(&self) -> &[[f64; 3]] {
        &self.used_strip_swirl_velocity_mps
    }

    #[must_use]
    pub fn warm_start(&self) -> WakeWarmStart {
        WakeWarmStart {
            rotor: self.solution.warm_start(),
            strip_swirl_velocity_mps: self.used_strip_swirl_velocity_mps.clone(),
        }
    }
}

fn separated_disks(axes: DiskAxes, disks: &[PropDisk; 2], radius: f64) -> Result<(), Refusal> {
    let directions = axes.directions();
    for k in 0..2 {
        let other = 1 - k;
        let offset = vec(disks[other].center_m) - vec(disks[k].center_m);
        let axial = offset.dot(vec(directions[k]));
        let radial = (offset - vec(directions[k]).scale(axial)).norm();
        let half_depth = radius * vec(directions[k]).cross(vec(directions[other])).norm();
        // A conservative disk bounding volume deliberately refuses marginal
        // overlap, rather than averaging a nonzero incident swirl to zero.
        if ![axial, radial, half_depth].into_iter().all(f64::is_finite) {
            return Err(reject("swirl-wake-invalid", "nonfinite disk geometry"));
        }
        if radial < 2.0 * radius && -axial + half_depth >= 0.0 {
            return Err(reject(
                "swirl-wake-disk-interaction",
                "receiving disk may intersect another downstream swirl cylinder",
            ));
        }
    }
    Ok(())
}

/// Jointly solve field -> wing -> native rotors -> torque-transport field.
///
/// A bounded relaxed outer field iteration wraps the unchanged inner axial
/// two-way solve. Every outer trial uses the caller's fixed native branch seeds.
/// The inner solve returns its ORIGINAL disk axial samples, without recovering
/// them from thrust/power ratios. Prepared wing geometry/operator are reused.
/// spin_sign is +/-1 about each positive-thrust axis, not a gear-ratio guess.
///
/// # Errors
/// Existing native refusals; invalid controls/seeds/layout; unsupported
/// receiving-disk swirl; power-budget violation; outside-window wing samples;
/// "swirl-wake-nonconvergence" on exhaustion. No partial loads are returned.
#[allow(clippy::too_many_arguments)]
pub fn coupled_step(
    axes: DiskAxes,
    spin_sign: [f64; 2],
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    frames: Option<&[SectionFrame]>,
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    native: Controls,
    controls: WakeControls,
    warm: Option<&WakeWarmStart>,
    incident: IncidentInflow<'_>,
) -> Result<CoupledWakeStep, Refusal> {
    coupled_step_with_transport(
        WakeTransport::FixedCylinder,
        axes,
        spin_sign,
        op,
        panels,
        strips,
        closure,
        frames,
        rotor,
        disks,
        freestream,
        rho,
        spec,
        native,
        controls,
        warm,
        incident,
    )
}

/// Select source-field physics while sharing the same native nested iteration.
///
/// MaterialStreamtubes removes the ACTUALLY evaluated downstream legacy wash
/// before adding its full native-annulus field. Original upstream wash remains.
/// Rejected outer trials never become new native branch anchors.
///
/// # Errors
/// The same admission, source-window, energy-budget and convergence refusals as
/// coupled_step, including expanded wake/receiver-disk overlap.
#[allow(clippy::too_many_arguments)]
pub fn coupled_step_with_transport(
    transport: WakeTransport,
    axes: DiskAxes,
    spin_sign: [f64; 2],
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    frames: Option<&[SectionFrame]>,
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    native: Controls,
    controls: WakeControls,
    warm: Option<&WakeWarmStart>,
    incident: IncidentInflow<'_>,
) -> Result<CoupledWakeStep, Refusal> {
    controls.admit()?;
    if !spin_sign.into_iter().all(|v| v == -1.0 || v == 1.0)
        || incident.strip_velocity_mps.len() != strips.len()
        || strips.is_empty()
        || strips.iter().any(|s| {
            s.panel_indices.is_empty() || s.panel_indices.iter().any(|&p| p >= panels.len())
        })
        || warm.is_some_and(|w| {
            w.strip_swirl_velocity_mps.len() != strips.len()
                || !w
                    .strip_swirl_velocity_mps
                    .iter()
                    .flatten()
                    .all(|v| v.is_finite())
        })
    {
        return Err(reject(
            "swirl-wake-invalid",
            "invalid spin, layout or field seed",
        ));
    }
    rotor.admit().map_err(|e| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    })?;
    separated_disks(axes, disks, rotor.radius_m)?;
    let points: Vec<_> = strips
        .iter()
        .map(|s| {
            let panel = panels[s.panel_indices[0]];
            array((vec(panel.a) + vec(panel.b)).scale(0.5))
        })
        .collect();
    // Upstream axial wash retains its original owner. Prepare these factors
    // once; the private callback returns the actual accepted x, never G(x).
    let upstream = if transport == WakeTransport::MaterialStreamtubes {
        let centers = [disks[0].center_m, disks[1].center_m];
        let directions = axes.directions();
        points
            .iter()
            .map(|&point| {
                let mut factors = axes.wash_factors(point, centers, rotor.radius_m)?;
                for k in 0..2 {
                    if (vec(point) - vec(centers[k])).dot(vec(directions[k])) <= 0.0 {
                        factors[k] = 0.0;
                    }
                }
                Ok(factors)
            })
            .collect::<Result<Vec<_>, Refusal>>()?
    } else {
        Vec::new()
    };
    let mut used = warm.map_or_else(
        || vec![[0.0; 3]; strips.len()],
        |w| w.strip_swirl_velocity_mps.clone(),
    );
    let mut last_residual = f64::INFINITY;
    for evaluation in 1..=controls.maximum_field_evaluations {
        let velocity: Vec<_> = incident
            .strip_velocity_mps
            .iter()
            .zip(&used)
            .map(|(&v, &u)| array(vec(v) + vec(u)))
            .collect();
        let (mut solution, axial, applied_wash) = coupled_impl(
            Some(axes),
            op,
            panels,
            strips,
            closure,
            frames,
            rotor,
            disks,
            freestream,
            rho,
            spec,
            native,
            warm.map(|w| &w.rotor),
            IncidentInflow {
                strip_velocity_mps: &velocity,
                disk_velocity_mps: incident.disk_velocity_mps,
            },
        )?;
        let directions = axes.directions();
        let make = |k: usize| {
            AnnularWake::from_solution(
                rotor,
                &solution.rotors[k],
                disks[k].center_m,
                directions[k],
                spin_sign[k],
                rho,
                axial[k],
                solution.coupled.w_slip[k],
                controls,
                transport,
            )
        };
        let wakes = [make(0)?, make(1)?];
        if transport == WakeTransport::MaterialStreamtubes {
            for k in 0..2 {
                let offset = vec(disks[1 - k].center_m) - vec(disks[k].center_m);
                let axial = offset.dot(vec(directions[k]));
                let radial = (offset - vec(directions[k]).scale(axial)).norm();
                let depth =
                    rotor.radius_m * vec(directions[k]).cross(vec(directions[1 - k])).norm();
                let radius = wakes[k]
                    .streamtubes
                    .as_ref()
                    .expect("selected material mode")
                    .maximum_radius_m();
                if radial < radius + rotor.radius_m && -axial + depth >= 0.0 {
                    return Err(reject(
                        "swirl-wake-disk-interaction",
                        "receiving disk may intersect an expanded material streamtube",
                    ));
                }
            }
        }
        let mut predicted = Vec::with_capacity(points.len());
        last_residual = 0.0;
        for (i, (&point, &current)) in points.iter().zip(&used).enumerate() {
            let left = vec(wakes[0].relative_velocity_at(point)?);
            let right = vec(wakes[1].relative_velocity_at(point)?);
            if left.norm() > 0.0 && right.norm() > 0.0 {
                return Err(reject(
                    "swirl-wake-interaction",
                    "overlapping angular source fields require a wake-mixing model",
                ));
            }
            let mut next = left + right;
            if transport == WakeTransport::MaterialStreamtubes {
                let actual = solution
                    .coupled
                    .strip_velocity_mps
                    .as_ref()
                    .expect("native coupling retains accepted force velocities")[i];
                let legacy = vec(actual) - vec(velocity[i]);
                let retained_upstream = axes.velocity([0.0; 3], upstream[i], applied_wash);
                next = next + vec(retained_upstream) - legacy;
            }
            last_residual = last_residual.max((next - vec(current)).norm());
            predicted.push(array(next));
        }
        if last_residual <= controls.velocity_tolerance_mps {
            let mut identity = solution.coupled.spec_digest.as_bytes().to_vec();
            identity.extend_from_slice(&controls.maximum_field_evaluations.to_le_bytes());
            for value in [
                spin_sign[0],
                spin_sign[1],
                controls.relaxation,
                controls.velocity_tolerance_mps,
                controls.absolute_power_tolerance_w,
                controls.relative_power_tolerance,
            ] {
                identity.extend_from_slice(&value.to_bits().to_le_bytes());
            }
            solution.coupled.spec_digest = hash_domain(
                match transport {
                    WakeTransport::FixedCylinder => {
                        "org.frankensim.fs-flyer.torque-transport-wake.v1"
                    }
                    WakeTransport::MaterialStreamtubes => {
                        "org.frankensim.fs-flyer.material-streamtube-wake.v1"
                    }
                },
                &identity,
            )
            .to_hex();
            return Ok(CoupledWakeStep {
                solution,
                wakes,
                used_strip_swirl_velocity_mps: used,
                field_evaluations: evaluation,
                field_residual_mps: last_residual,
            });
        }
        for (current, next) in used.iter_mut().zip(predicted) {
            *current =
                array(vec(*current) + (vec(next) - vec(*current)).scale(controls.relaxation));
        }
    }
    Err(reject(
        "swirl-wake-nonconvergence",
        &format!(
            "field residual {last_residual} m/s after {} complete evaluations",
            controls.maximum_field_evaluations
        ),
    ))
}

#[cfg(test)]
mod projection_tests {
    use super::*;

    // Analytic lossless annular angular actuator with constant disk u,
    // constant axial transport U and equal radial cells. This isolates
    // projection error from the disclosed Prandtl/fixed-cylinder errors.
    fn projection(count: usize, u: f64) -> AnnularWake {
        let density = 1.225;
        let axial = 10.0;
        let width = 0.6 / count as f64;
        let cells = (0..count)
            .map(|i| {
                let inner = 0.4 + i as f64 * width;
                let outer = 0.4 + (i + 1) as f64 * width;
                let radius = 0.5 * (inner + outer);
                let torque =
                    4.0 * core::f64::consts::PI * density * axial * radius * radius * u * width;
                Annulus {
                    inner_radius_m: inner,
                    outer_radius_m: outer,
                    torque_nm: torque,
                    native_swirl_power_w: u / radius * torque,
                }
            })
            .collect();
        AnnularWake {
            center_m: [0.0; 3],
            axis: [1.0, 0.0, 0.0],
            spin_sign: 1.0,
            radius_m: 1.0,
            density,
            incident_axial_mps: axial,
            disk_induction_mps: 0.0,
            cells,
            maximum_power_defect_w: 0.0,
            streamtubes: None,
        }
    }

    #[test]
    fn no_loss_matching_flow_recovers_twice_disk_u_and_zero_limit() {
        for u in [0.0, 1e-9, 0.2, -0.2] {
            let wake = projection(16, u);
            for cell in wake.annuli() {
                let radius = 0.5 * (cell.inner_radius_m + cell.outer_radius_m);
                let velocity = wake.relative_velocity_at([-0.5, radius, 0.0]).unwrap();
                assert!((velocity[2] + 2.0 * u).abs() <= 2e-14 * u.abs() + 1e-30);
                assert_eq!(velocity[0], 0.0);
                assert_eq!(velocity[1], 0.0);
            }
            let flux = wake.cross_section(0.5).unwrap();
            if u == 0.0 {
                assert_eq!(flux.angular_momentum_flux_nm, 0.0);
                assert_eq!(flux.resolved_swirl_power_w, 0.0);
            } else {
                assert!(flux.angular_momentum_flux_nm * u > 0.0);
                assert!(flux.resolved_swirl_power_w > 0.0);
            }
        }
    }

    #[test]
    fn finite_annulus_resolved_energy_error_converges_quadratically() {
        let mut errors = Vec::new();
        for count in [4, 8, 16, 32] {
            let wake = projection(count, 0.2);
            let flux = wake.cross_section(0.5).unwrap();
            let exact = 0.5 * 1.225 * 10.0 * 0.4 * 0.4 * core::f64::consts::PI * (1.0 - 0.4 * 0.4);
            assert!((flux.native_swirl_power_w - exact).abs() < 2e-14);
            let error = flux.resolved_swirl_power_w - exact;
            assert!(error > 0.0);
            errors.push(error);
        }
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!((3.9..4.1).contains(&ratio), "{ratio}");
        }
    }
}
