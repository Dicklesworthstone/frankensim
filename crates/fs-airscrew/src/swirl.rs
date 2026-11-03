//! Axial and angular annular momentum coupled to the original blade law.
//!
//! With disk axial induction `w` and tangential induction `u`, the blade sees
//! `(V+w, Omega*r-u)`. The native Prandtl factor F enters BOTH balances:
//! `dT/dr = 4 pi rho F r (V+w) w` and
//! `dQ/dr = 4 pi rho F r² (V+w) u`.
//! These are the propeller-sign form of Moriarty and Hansen, NREL/TP-500-36881,
//! equations (9), (10), <https://www.nlr.gov/docs/fy05osti/36881.pdf>.
//! Propeller kinematics also follow NASA's 2016 high-lift design presentation,
//! <https://ntrs.nasa.gov/citations/20160007767>, slides 11 and 13.
//!
//! The two inductions change the actual section forces; angular momentum is
//! not inferred afterward from unchanged no-swirl thrust/torque. A bounded
//! local Newton solve retains the original section law, radial quadrature and
//! Prandtl convention. Explicit seeds allow branch continuation; no global
//! root uniqueness, resolved helical wake, radial pressure, contraction,
//! transient swirl, compressibility, skew or historical validation is claimed.

use super::{MAX_BRAKING_INDUCTION, Refusal, Rotor, StationProblem, StationSample, refuse};
use fs_math::det;

/// Hard count of original blade-force evaluations per station, including
/// Jacobian probes and rejected line-search trials.
pub const MAX_SWIRL_EVALUATIONS: u32 = 256;

/// Actual disk-induced velocity components in the positive blade-spin frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Induction {
    /// Axial induction [m/s]; negative for admitted braking annuli.
    pub axial_mps: f64,
    /// Tangential induction along positive blade rotation [m/s].
    /// Negative values transfer angular momentum from air to the shaft.
    pub tangential_mps: f64,
}

/// Bounded local solve; final acceptance checks both original momentum laws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Controls {
    /// Original force evaluations admitted for each station, in 3..=256.
    pub maximum_station_evaluations: u32,
    /// Relative tolerance applied separately to actual thrust and torque
    /// versus their corresponding annular momentum fluxes, in [1e-13,1e-6].
    /// A 16-epsilon section-force roundoff floor admits genuine zero crossings.
    pub relative_tolerance: f64,
}
impl Default for Controls {
    fn default() -> Self {
        Self {
            maximum_station_evaluations: MAX_SWIRL_EVALUATIONS,
            relative_tolerance: 1e-10,
        }
    }
}

/// One complete accepted annulus, with SI forces and unrepaired residuals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Station {
    /// Original radial fraction.
    pub r_over_r: f64,
    /// Actual disk induction used by the blade law and both momentum fluxes.
    pub induction: Induction,
    /// Original blade-force evaluations, including unsuccessful probes.
    pub evaluations: u32,
    /// Original section incidence [rad].
    pub alpha_rad: f64,
    /// Native combined tip/root Prandtl factor.
    pub prandtl_f: f64,
    /// Blade thrust over this radial quadrature interval [N].
    pub thrust_n: f64,
    /// Required blade torque over the same interval [N m].
    pub torque_nm: f64,
    /// Prandtl-corrected axial momentum flux over the same interval [N].
    pub axial_momentum_flux_n: f64,
    /// Prandtl-corrected angular momentum flux over the same interval [N m].
    pub angular_momentum_flux_nm: f64,
    /// Actual blade thrust minus original axial momentum flux [N].
    pub axial_residual_n: f64,
    /// Actual blade torque minus original angular momentum flux [N m].
    pub angular_residual_nm: f64,
}

/// Blade power including the angular-momentum channel, all in watts.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Power {
    /// Signed drive power `Omega Q`.
    pub shaft_w: f64,
    /// Signed axial thrust power `V T`.
    pub thrust_w: f64,
    /// Actual sum of `w dT` over the accepted annuli.
    pub axial_induced_w: f64,
    /// Actual sum of `(u/r) dQ`; a wake transfer, not blade heating.
    pub swirl_w: f64,
    /// Original section profile drag times its relative air speed.
    pub profile_drag_w: f64,
    /// Shaft minus the four original component powers; numerical error.
    pub residual_w: f64,
}

/// Steady swirling rotor result in the same CT/CQ/J conventions as BEMT.
#[derive(Clone, Debug, PartialEq)]
pub struct Solution {
    /// One-rotor axial thrust [N].
    pub thrust_n: f64,
    /// Required one-rotor torque in the positive spin convention [N m].
    pub torque_nm: f64,
    /// `T/(rho n² D⁴)`, with n in revolutions per second.
    pub ct: f64,
    /// `Q/(rho n² D⁵)`.
    pub cq: f64,
    /// `V/(n D)`.
    pub j: f64,
    /// Actual section-power decomposition, including swirl transfer.
    pub power: Power,
    /// Ordered original station results; suitable for explicit continuation.
    pub stations: Vec<Station>,
}

fn error(code: &'static str, station: usize, detail: &str) -> Refusal {
    refuse(
        code,
        format!("swirl station {station}: {detail}"),
        "supply normal axial flow and forward blade-relative rotation; inspect the seed or refine the operating-point continuation",
    )
}

struct Equation<'a> {
    rotor: &'a Rotor,
    index: usize,
    density: f64,
    axial: f64,
    rotation: f64,
    radius: f64,
    width: f64,
    controls: Controls,
    evaluations: u32,
}

#[derive(Clone, Copy)]
struct Sample {
    state: Induction,
    blade: StationSample,
    /// Original momentum residuals divided only by the common geometric
    /// coefficients; both coordinates have units m²/s² for Newton's solve.
    residual: [f64; 2],
    thrust: f64,
    torque: f64,
    momentum: [f64; 2],
}
impl Sample {
    fn norm(self) -> f64 {
        self.residual[0].abs().max(self.residual[1].abs())
    }
    fn converged(self, tolerance: f64, radius: f64) -> bool {
        // Lift and drag can cancel one projected component while the blade
        // still carries a nonzero resultant. Use only the floating-point
        // roundoff of that original force for the zero-crossing floor, never
        // a fixed engineering tolerance or repaired momentum residual.
        let roundoff = 16.0 * f64::EPSILON * (self.thrust.abs() + self.torque.abs() / radius);
        (self.thrust - self.momentum[0]).abs()
            <= tolerance * (self.thrust.abs() + self.momentum[0].abs()) + roundoff
            && (self.torque - self.momentum[1]).abs()
                <= tolerance * (self.torque.abs() + self.momentum[1].abs()) + roundoff * radius
    }
}
impl Equation<'_> {
    fn in_domain(&self, state: Induction) -> bool {
        state.axial_mps.is_finite()
            && state.tangential_mps.is_finite()
            && state.axial_mps >= -0.5 * self.axial
            && state.tangential_mps < self.rotation
    }

    fn evaluate(&mut self, state: Induction) -> Result<Sample, Refusal> {
        if self.evaluations == self.controls.maximum_station_evaluations {
            return Err(error(
                "swirl-nonconvergence",
                self.index,
                "force-evaluation budget exhausted",
            ));
        }
        self.evaluations += 1;
        if !self.in_domain(state) {
            return Err(error(
                "swirl-wake-domain",
                self.index,
                "reversed throughflow or blade-relative rotation",
            ));
        }
        let native = StationProblem {
            rotor: self.rotor,
            station: &self.rotor.stations[self.index],
            index: self.index,
            rho: self.density,
            axial: self.axial,
            tangential: self.rotation - state.tangential_mps,
            radius: self.radius,
            root_radius: 0.8 * self.rotor.stations[0].r_over_r * self.rotor.radius_m,
        };
        let blade = native.evaluate(state.axial_mps)?;
        let coefficient = 4.0 * core::f64::consts::PI * self.density * blade.prandtl * self.radius;
        let through = self.axial + state.axial_mps;
        let rotational = through * state.tangential_mps;
        let momentum = [
            coefficient * through * state.axial_mps * self.width,
            coefficient * self.radius * rotational * self.width,
        ];
        let residual = [
            blade.balance,
            blade.torque / (coefficient * self.radius) - rotational,
        ];
        let thrust = blade.thrust * self.width;
        let torque = blade.torque * self.width;
        if ![
            coefficient,
            momentum[0],
            momentum[1],
            residual[0],
            residual[1],
            thrust,
            torque,
        ]
        .into_iter()
        .all(f64::is_finite)
            || coefficient <= 0.0
        {
            return Err(error(
                "swirl-numerical-invalid",
                self.index,
                "unrepresentable momentum equation",
            ));
        }
        Ok(Sample {
            state,
            blade,
            residual,
            thrust,
            torque,
            momentum,
        })
    }

    fn solve(&mut self, initial: Induction) -> Result<Sample, Refusal> {
        let mut current = self.evaluate(initial)?;
        loop {
            if current.converged(self.controls.relative_tolerance, self.radius) {
                if current.state.axial_mps < -MAX_BRAKING_INDUCTION * self.axial
                    || current.blade.alpha.abs() > 1.3
                    || self.axial + current.state.axial_mps <= 0.0
                {
                    return Err(error(
                        "swirl-wake-domain",
                        self.index,
                        "converged root exceeds ordinary braking, section-angle or positive-throughflow domain",
                    ));
                }
                return Ok(current);
            }
            // Forward differences use the physical velocity scale; the final
            // gate always checks the original force and torque, never an update.
            let h = 1e-6 * self.rotation.abs().max(self.axial).max(1.0);
            let a = self.evaluate(Induction {
                axial_mps: current.state.axial_mps + h,
                ..current.state
            })?;
            let signed_h = if current.state.tangential_mps + h < self.rotation {
                h
            } else {
                -h
            };
            let b = self.evaluate(Induction {
                tangential_mps: current.state.tangential_mps + signed_h,
                ..current.state
            })?;
            let j00 = (a.residual[0] - current.residual[0]) / h;
            let j10 = (a.residual[1] - current.residual[1]) / h;
            let j01 = (b.residual[0] - current.residual[0]) / signed_h;
            let j11 = (b.residual[1] - current.residual[1]) / signed_h;
            let determinant = j00 * j11 - j01 * j10;
            let change = [
                (j01 * current.residual[1] - j11 * current.residual[0]) / determinant,
                (j10 * current.residual[0] - j00 * current.residual[1]) / determinant,
            ];
            if determinant == 0.0 || !change.into_iter().all(f64::is_finite) {
                return Err(error(
                    "swirl-nonconvergence",
                    self.index,
                    "singular local momentum Jacobian",
                ));
            }
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..20 {
                let state = Induction {
                    axial_mps: current.state.axial_mps + fraction * change[0],
                    tangential_mps: current.state.tangential_mps + fraction * change[1],
                };
                if self.in_domain(state) {
                    let trial = self.evaluate(state)?;
                    if trial.norm() < current.norm()
                        || trial.converged(self.controls.relative_tolerance, self.radius)
                    {
                        accepted = Some(trial);
                        break;
                    }
                }
                fraction *= 0.5;
            }
            current = accepted.ok_or_else(|| {
                error(
                    "swirl-nonconvergence",
                    self.index,
                    "bounded local line search could not reduce the original momentum residual",
                )
            })?;
        }
    }
}

/// Solve steady axial and angular momentum with the native rotor/section law.
///
/// `warm` supplies one original induction per station for local continuation.
/// Without it, the first station starts at `(0.1 max(V,Omega*r),0)`; subsequent
/// stations inherit axial induction and angular velocity from their neighbor.
/// This is a local root policy, not a guarantee of unique stalled-blade roots.
/// Positive density/shaft speed, V >= 0 and the original rotor domain apply.
/// Moderate windmilling is signed; braking beyond `a=0.4`, reversed relative
/// tangential flow, and static zero-throughflow refuse. No result is partially
/// published when a later station fails.
pub fn bemt_solve(
    rotor: &Rotor,
    density_kg_m3: f64,
    axial_mps: f64,
    omega_rad_s: f64,
    controls: Controls,
    warm: Option<&[Induction]>,
) -> Result<Solution, Refusal> {
    rotor.admit()?;
    if ![
        density_kg_m3,
        axial_mps,
        omega_rad_s,
        controls.relative_tolerance,
    ]
    .into_iter()
    .all(f64::is_finite)
        || density_kg_m3 <= 0.0
        || axial_mps < 0.0
        || omega_rad_s <= 0.0
        || !(3..=MAX_SWIRL_EVALUATIONS).contains(&controls.maximum_station_evaluations)
        || !(1e-13..=1e-6).contains(&controls.relative_tolerance)
        || warm.is_some_and(|s| s.len() != rotor.stations.len())
    {
        return Err(error(
            "swirl-input-invalid",
            0,
            "invalid operating point, controls or seed count",
        ));
    }
    let mut solution = Solution {
        thrust_n: 0.0,
        torque_nm: 0.0,
        ct: 0.0,
        cq: 0.0,
        j: 0.0,
        power: Power::default(),
        stations: Vec::with_capacity(rotor.stations.len()),
    };
    let mut preceding = Induction {
        axial_mps: 0.1 * axial_mps.max(omega_rad_s * rotor.radius_m * rotor.stations[0].r_over_r),
        tangential_mps: 0.0,
    };
    let mut preceding_radius = rotor.stations[0].r_over_r * rotor.radius_m;
    for (index, station) in rotor.stations.iter().enumerate() {
        let radius = rotor.radius_m * station.r_over_r;
        let rotation = omega_rad_s * radius;
        let width = rotor.station_width_m(index);
        if ![radius, rotation, width]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
        {
            return Err(error(
                "swirl-input-invalid",
                index,
                "unrepresentable station geometry or speed",
            ));
        }
        let initial = warm.map_or(
            Induction {
                axial_mps: preceding.axial_mps,
                tangential_mps: preceding.tangential_mps * (radius / preceding_radius),
            },
            |s| s[index],
        );
        let mut equation = Equation {
            rotor,
            index,
            density: density_kg_m3,
            axial: axial_mps,
            rotation,
            radius,
            width,
            controls,
            evaluations: 0,
        };
        let sample = equation.solve(initial)?;
        let state = sample.state;
        solution.thrust_n += sample.thrust;
        solution.torque_nm += sample.torque;
        solution.power.axial_induced_w += state.axial_mps * sample.thrust;
        solution.power.swirl_w += (state.tangential_mps / radius) * sample.torque;
        let relative_axial = axial_mps + state.axial_mps;
        let relative_tangent = rotation - state.tangential_mps;
        let speed =
            det::sqrt(relative_axial * relative_axial + relative_tangent * relative_tangent);
        solution.power.profile_drag_w += sample.blade.profile_drag * width * speed;
        solution.stations.push(Station {
            r_over_r: station.r_over_r,
            induction: state,
            evaluations: equation.evaluations,
            alpha_rad: sample.blade.alpha,
            prandtl_f: sample.blade.prandtl,
            thrust_n: sample.thrust,
            torque_nm: sample.torque,
            axial_momentum_flux_n: sample.momentum[0],
            angular_momentum_flux_nm: sample.momentum[1],
            axial_residual_n: sample.thrust - sample.momentum[0],
            angular_residual_nm: sample.torque - sample.momentum[1],
        });
        preceding = state;
        preceding_radius = radius;
    }
    let n = omega_rad_s / core::f64::consts::TAU;
    let diameter = 2.0 * rotor.radius_m;
    solution.ct = solution.thrust_n / (density_kg_m3 * n * n * diameter.powi(4));
    solution.cq = solution.torque_nm / (density_kg_m3 * n * n * diameter.powi(5));
    solution.j = axial_mps / (n * diameter);
    solution.power.shaft_w = solution.torque_nm * omega_rad_s;
    solution.power.thrust_w = solution.thrust_n * axial_mps;
    solution.power.residual_w = solution.power.shaft_w
        - solution.power.thrust_w
        - solution.power.axial_induced_w
        - solution.power.swirl_w
        - solution.power.profile_drag_w;
    if ![
        solution.ct,
        solution.cq,
        solution.j,
        solution.power.shaft_w,
        solution.power.thrust_w,
        solution.power.axial_induced_w,
        solution.power.swirl_w,
        solution.power.profile_drag_w,
        solution.power.residual_w,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(error(
            "swirl-numerical-invalid",
            rotor.stations.len() - 1,
            "unrepresentable rotor loads or power",
        ));
    }
    Ok(solution)
}
