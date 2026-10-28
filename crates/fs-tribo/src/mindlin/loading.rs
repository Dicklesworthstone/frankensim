//! Hertz--Mindlin memory under explicit normal loading and unloading.
//!
//! The MDR representation has tangential stiffness density `2 G_MDR` on the
//! positive half-line, where `G_MDR = 4 / sum((2-nu_i)/G_i)`. For a paraboloid
//! of equivalent radius R, its elastic spring displacement is bounded by
//! `|e(x)| <= mu E_star/G_MDR * (a*a-x*x)/R`. Normal growth adds stress-free
//! springs; normal unloading clips existing strains and releases their stored
//! energy as physical microslip heat. The retained strain is piecewise affine
//! in x squared, so all history changes and energy integrals are analytic.
//!
//! `normal_step` is a monotonic normal-load segment at fixed macroscopic shear
//! displacement. `TangentialLaw::advance` is a shear segment at the retained
//! normal load, using the original fixed-load owner. Callers explicitly choose
//! that path ordering and refine it when approximating simultaneous motion.
//! There is no implied instantaneous rescaling, thermal coefficient evolution,
//! two-dimensional shear, or normal/tangential elastic coupling. Normal work
//! and normal Hertz storage are separate from the tangential work port.
//!
//! This is the decoupled Cattaneo--Mindlin/MDR approximation described by
//! Popov, Hess and Willert (2019), section 4.4, equations 4.21--4.25:
//! <https://doi.org/10.1007/978-3-662-58709-6_4>. In that convention
//! `k_t = 2 a G_MDR = 8 a / sum((2-nu_i)/G_i)`.

use fs_math::det;

use super::{
    MAX_HISTORY_INTERVALS, MemoryInterval, MindlinDeresiewiczPoint, MindlinPointParameters,
    MindlinState, Moments,
};
use crate::{
    FrictionRegime, TangentialLaw, TangentialStep, TangentialStepInput, TriboError, checked_add,
    checked_mul, finite, nonnegative_finite, positive_finite,
};

/// Immutable geometry and elastic/interface parameters for a Hertz contact.
/// The caller owns admission of the small-strain, half-space approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzMindlinParameters {
    /// Equivalent positive sphere/paraboloid radius of curvature [m].
    pub effective_radius_m: f64,
    /// Shear moduli of the two ordered elastic bodies [Pa].
    pub shear_moduli_pa: [f64; 2],
    /// Their Poisson ratios, with `-1 < nu <= 1/2`.
    pub poisson_ratios: [f64; 2],
    /// Constant local Coulomb coefficient; zero admits frictionless contact.
    pub friction_coefficient: f64,
    /// Explicit analytic-memory ceiling in `1..=MAX_HISTORY_INTERVALS`.
    pub maximum_history_intervals: usize,
}

/// Analytic Hertz contact with evolving normal load and retained scalar shear.
#[derive(Clone, Debug, PartialEq)]
pub struct HertzMindlinContact {
    parameters: HertzMindlinParameters,
    normal_modulus_pa: f64,
}

/// Accepted geometry and complete shear history, including the open state.
#[derive(Clone, Debug, PartialEq)]
pub struct HertzMindlinState {
    owner: HertzMindlinContact,
    fixed: Option<MindlinState>,
    slip_m: f64,
    indentation_m: f64,
    normal_energy_j: f64,
}

impl HertzMindlinState {
    /// Retained compressive normal load [N]; zero means open/unloaded contact.
    #[must_use]
    pub fn normal_force_n(&self) -> f64 {
        self.fixed
            .as_ref()
            .map_or(0.0, |s| s.owner.parameters.normal_force_n)
    }
    /// Retained circular contact radius [m], zero when open.
    #[must_use]
    pub fn contact_radius_m(&self) -> f64 {
        self.fixed
            .as_ref()
            .map_or(0.0, |s| s.owner.parameters.contact_radius_m)
    }
    /// Retained Hertz indentation [m].
    #[must_use]
    pub const fn indentation_m(&self) -> f64 {
        self.indentation_m
    }
    /// Normal elastic storage [J], excluded from `TangentialLaw::stored_energy_j`.
    #[must_use]
    pub const fn normal_energy_j(&self) -> f64 {
        self.normal_energy_j
    }
    /// Tangential endpoint effort entering the interface [N].
    #[must_use]
    pub fn endpoint_force_n(&self) -> f64 {
        self.fixed
            .as_ref()
            .map_or(0.0, MindlinState::endpoint_force_n)
    }
    /// Accumulated solved scalar shear, including motion while open [m].
    #[must_use]
    pub const fn slip_m(&self) -> f64 {
        self.slip_m
    }
    /// Number of retained analytic intervals; zero while open.
    #[must_use]
    pub fn history_intervals(&self) -> usize {
        self.fixed
            .as_ref()
            .map_or(0, MindlinState::history_intervals)
    }
    /// Current original fixed-load state, exposed without mutation.
    #[must_use]
    pub const fn fixed_load_state(&self) -> Option<&MindlinState> {
        self.fixed.as_ref()
    }
}

/// One monotonic normal change at fixed macroscopic tangential displacement.
#[derive(Clone, Debug, PartialEq)]
pub struct MindlinNormalStep {
    /// Candidate geometry and history; errors leave the input untouched.
    pub state: HertzMindlinState,
    /// Signed external normal work, equal to the Hertz storage change [J].
    pub normal_work_j: f64,
    /// Normal Hertz elastic-storage change [J].
    pub normal_stored_change_j: f64,
    /// Shear storage change at zero macroscopic shear work [J].
    pub tangential_stored_change_j: f64,
    /// Physical microslip heat released by unloading, never normal damping [J].
    pub dissipated_j: f64,
    /// Both storage changes plus heat minus external normal work [J].
    pub closure_residual_j: f64,
}

impl HertzMindlinContact {
    /// Admit the immutable elastic Hertz geometry and Coulomb interface.
    pub fn new(parameters: HertzMindlinParameters) -> Result<Self, TriboError> {
        positive_finite(parameters.effective_radius_m, "mindlin_effective_radius_m")?;
        nonnegative_finite(
            parameters.friction_coefficient,
            "mindlin_friction_coefficient",
        )?;
        if !(1..=MAX_HISTORY_INTERVALS).contains(&parameters.maximum_history_intervals) {
            return Err(invalid("mindlin_maximum_history_intervals"));
        }
        let mut compliance = 0.0;
        for (g, nu) in parameters
            .shear_moduli_pa
            .into_iter()
            .zip(parameters.poisson_ratios)
        {
            positive_finite(g, "mindlin_shear_modulus_pa")?;
            if !(nu > -1.0 && nu <= 0.5) {
                return Err(invalid("mindlin_poisson_ratio"));
            }
            compliance = checked_add(
                compliance,
                0.5 * (1.0 - nu) / g,
                "mindlin_normal_compliance",
            )?;
        }
        let normal_modulus_pa = 1.0 / compliance;
        positive_finite(normal_modulus_pa, "mindlin_effective_normal_modulus_pa")?;
        Ok(Self {
            parameters,
            normal_modulus_pa,
        })
    }

    /// Immutable curvature, ordered elastic values, coefficient and memory cap.
    #[must_use]
    pub const fn parameters(&self) -> HertzMindlinParameters {
        self.parameters
    }

    /// Effective Hertz modulus `1 / sum((1-nu_i)/(2 G_i))` [Pa].
    #[must_use]
    pub const fn normal_modulus_pa(&self) -> f64 {
        self.normal_modulus_pa
    }

    /// Change normal load monotonically while holding macroscopic shear fixed.
    /// Newly contacting regions start with zero shear; detached regions lose
    /// their memory only after its entire recoverable energy is exported.
    /// Recontact therefore cannot resurrect shear from a previous encounter.
    pub fn normal_step(
        &self,
        state: &HertzMindlinState,
        normal_force_n: f64,
    ) -> Result<MindlinNormalStep, TriboError> {
        self.check_owner(state)?;
        nonnegative_finite(normal_force_n, "normal_force_n")?;
        if normal_force_n == state.normal_force_n() {
            return Ok(MindlinNormalStep {
                state: state.clone(),
                normal_work_j: 0.0,
                normal_stored_change_j: 0.0,
                tangential_stored_change_j: 0.0,
                dissipated_j: 0.0,
                closure_residual_j: 0.0,
            });
        }
        let before = self.stored_energy_j(state);
        let (fixed, indentation_m, normal_energy_j, heat) = if normal_force_n == 0.0 {
            (None, 0.0, 0.0, before)
        } else {
            let radius_cube = checked_mul(
                0.75 * (normal_force_n / self.normal_modulus_pa),
                self.parameters.effective_radius_m,
                "mindlin_hertz_radius_cube",
            )?;
            positive_finite(radius_cube, "mindlin_hertz_radius_cube")?;
            let radius = det::pow(radius_cube, 1.0 / 3.0);
            let indentation = radius * (radius / self.parameters.effective_radius_m);
            positive_finite(indentation, "mindlin_hertz_indentation_m")?;
            let normal_energy =
                checked_mul(0.4 * normal_force_n, indentation, "mindlin_normal_energy_j")?;
            let law = MindlinDeresiewiczPoint::new(MindlinPointParameters {
                normal_force_n,
                contact_radius_m: radius,
                shear_moduli_pa: self.parameters.shear_moduli_pa,
                poisson_ratios: self.parameters.poisson_ratios,
                friction_coefficient: self.parameters.friction_coefficient,
                maximum_history_intervals: self.parameters.maximum_history_intervals,
            })?;
            let (fixed, released) = if let Some(old) = state.fixed.as_ref()
                && self.parameters.friction_coefficient > 0.0
            {
                transfer(old, law)?
            } else {
                let mut fresh = law.initial_state();
                fresh.slip_m = state.slip_m;
                (fresh, 0.0)
            };
            (Some(fixed), indentation, normal_energy, released)
        };
        let normal_change = normal_energy_j - state.normal_energy_j;
        finite(normal_change, "mindlin_normal_storage_change_j")?;
        let after = fixed.as_ref().map_or(0.0, |s| s.stored_energy_j);
        let tangential_change = after - before;
        let closure = tangential_change + heat;
        finite(closure, "mindlin_normal_step_closure_j")?;
        Ok(MindlinNormalStep {
            state: HertzMindlinState {
                owner: self.clone(),
                fixed,
                slip_m: state.slip_m,
                indentation_m,
                normal_energy_j,
            },
            normal_work_j: normal_change,
            normal_stored_change_j: normal_change,
            tangential_stored_change_j: tangential_change,
            dissipated_j: heat,
            closure_residual_j: closure,
        })
    }

    fn check_owner(&self, state: &HertzMindlinState) -> Result<(), TriboError> {
        if state.owner != *self {
            return Err(invalid("mindlin_loading_checkpoint_owner"));
        }
        Ok(())
    }
}

impl TangentialLaw for HertzMindlinContact {
    type State = HertzMindlinState;

    fn model_id(&self) -> &str {
        "hertz-mindlin-explicit-normal-path-v1"
    }

    fn initial_state(&self) -> Self::State {
        HertzMindlinState {
            owner: self.clone(),
            fixed: None,
            slip_m: 0.0,
            indentation_m: 0.0,
            normal_energy_j: 0.0,
        }
    }

    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        if state.owner != *self {
            f64::NAN
        } else {
            state.fixed.as_ref().map_or(0.0, |s| s.stored_energy_j)
        }
    }

    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.check_owner(state)?;
        if input.normal_n != state.normal_force_n() {
            return Err(invalid("mindlin_loading_requires_explicit_normal_step"));
        }
        if let Some(fixed) = &state.fixed {
            let step = fixed.owner.advance(fixed, input)?;
            let next = HertzMindlinState {
                owner: self.clone(),
                slip_m: step.state.slip_m,
                fixed: Some(step.state),
                indentation_m: state.indentation_m,
                normal_energy_j: state.normal_energy_j,
            };
            return Ok(TangentialStep {
                traction_n: step.traction_n,
                slip_increment_m: step.slip_increment_m,
                state: next,
                regime: step.regime,
                stored_change_j: step.stored_change_j,
                dissipated_j: step.dissipated_j,
                closure_residual_j: step.closure_residual_j,
                compatibility_residual_m: step.compatibility_residual_m,
                iterations: step.iterations,
            });
        }
        input.solve.validate()?;
        positive_finite(input.dt_s, "dt_s")?;
        nonnegative_finite(input.compliance_m_per_n, "compliance_m_per_n")?;
        for value in input.slip_increment_m {
            finite(value, "slip_increment_m")?;
        }
        if input.slip_increment_m[1] != 0.0 {
            return Err(invalid("mindlin_requires_proportional_first_axis_motion"));
        }
        if let Some(t) = input.contact_temperature_k {
            positive_finite(t, "temperature_kelvin")?;
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "contact_pressure_pa")?;
        }
        let mut next = state.clone();
        next.slip_m = checked_add(
            state.slip_m,
            input.slip_increment_m[0],
            "mindlin_accumulated_slip_m",
        )?;
        Ok(TangentialStep {
            traction_n: [0.0; 2],
            slip_increment_m: input.slip_increment_m,
            state: next,
            regime: FrictionRegime::Sliding,
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            closure_residual_j: 0.0,
            compatibility_residual_m: [0.0; 2],
            iterations: 0,
        })
    }
}

fn invalid(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}

// Work in the old r=1-(x/a_old)^2 coordinate. If the contact shrinks, its new
// local bound is ±(r-cutoff) in old displacement units. Every branch is affine;
// the two possible intersections therefore give the exact history partition.
fn transfer(
    old: &MindlinState,
    law: MindlinDeresiewiczPoint,
) -> Result<(MindlinState, f64), TriboError> {
    let load_ratio = law.parameters.normal_force_n / old.owner.parameters.normal_force_n;
    // For a small load change, subtracting two independently computed Hertz
    // radii destroys the relative accuracy of the released heat. Hertz gives
    // a_new^2/a_old^2 = (P_new/P_old)^(2/3); expm1 retains its small difference.
    let cutoff = if (0.5..=2.0).contains(&load_ratio) {
        -det::expm1((2.0 / 3.0) * det::ln(load_ratio))
    } else {
        let ratio = law.parameters.contact_radius_m / old.owner.parameters.contact_radius_m;
        1.0 - ratio * ratio
    };
    let area_ratio = 1.0 - cutoff;
    positive_finite(area_ratio, "mindlin_patch_area_ratio")?;
    let mut memory =
        Vec::with_capacity((old.memory.len() + 2).min(law.parameters.maximum_history_intervals));
    let mut heat = 0.0;
    if area_ratio > 1.0 {
        law.append_memory(
            &mut memory,
            MemoryInterval {
                upper: -cutoff / area_ratio,
                slope: 0.0,
                offset: 0.0,
            },
        )?;
    }
    let mut lower = 0.0;
    for interval in &old.memory {
        let upper = interval.upper;
        if cutoff > lower {
            let detached_upper = upper.min(cutoff);
            if detached_upper > lower {
                let moments = Moments::new(lower, detached_upper);
                let mean = interval.slope * moments.mean + interval.offset;
                heat += 0.5
                    * old.owner.energy_scale_j
                    * moments.weight
                    * (mean * mean + interval.slope * interval.slope * moments.variance);
            }
        }
        let retained_lower = lower.max(cutoff);
        if retained_lower < upper {
            let mut cuts = [retained_lower, upper, upper, upper];
            let mut count = 2;
            if area_ratio < 1.0 {
                for direction in [-1.0, 1.0] {
                    let denominator = interval.slope - direction;
                    if denominator != 0.0 {
                        let root = (-direction * cutoff - interval.offset) / denominator;
                        if root > retained_lower && root < upper {
                            cuts[count] = root;
                            count += 1;
                        }
                    }
                }
            }
            cuts[..count].sort_by(f64::total_cmp);
            for pair in cuts[..count].windows(2) {
                let lo = pair[0];
                let hi = pair[1];
                if hi <= lo {
                    continue;
                }
                let mid = f64::midpoint(lo, hi);
                let old_strain = interval.slope * mid + interval.offset;
                let bound = mid - cutoff;
                let direction = if area_ratio < 1.0 && old_strain > bound {
                    1.0
                } else if area_ratio < 1.0 && old_strain < -bound {
                    -1.0
                } else {
                    0.0
                };
                let (slope, offset) = if direction != 0.0 {
                    heat += release_energy(old, interval, lo, hi, cutoff, direction);
                    (direction, 0.0)
                } else {
                    (
                        interval.slope,
                        (interval.offset + interval.slope * cutoff) / area_ratio,
                    )
                };
                let mapped_upper = if hi == 1.0 {
                    1.0
                } else {
                    ((hi - cutoff) / area_ratio).clamp(0.0, 1.0)
                };
                if mapped_upper <= memory.last().map_or(0.0, |piece| piece.upper) {
                    return Err(invalid("mindlin_normal_history_below_resolution"));
                }
                law.append_memory(
                    &mut memory,
                    MemoryInterval {
                        upper: mapped_upper,
                        slope,
                        offset,
                    },
                )?;
            }
        }
        lower = upper;
    }
    nonnegative_finite(heat, "mindlin_normal_release_heat_j")?;
    let allowance = 256.0 * f64::EPSILON * old.stored_energy_j;
    if heat > old.stored_energy_j + allowance {
        return Err(invalid("mindlin_normal_release_exceeds_storage"));
    }
    // Only the admitted rounding excess can reach this min; no physical heat
    // is discarded or guessed when an entire contact opens.
    heat = heat.min(old.stored_energy_j);
    let energy = old.stored_energy_j - heat;
    let mut endpoint = 0.0;
    let mut integrated_energy = 0.0;
    let mut lower = 0.0;
    for interval in &memory {
        let moments = Moments::new(lower, interval.upper);
        let mean = interval.slope * moments.mean + interval.offset;
        endpoint += moments.weight * mean;
        integrated_energy += 0.5
            * moments.weight
            * (mean * mean + interval.slope * interval.slope * moments.variance);
        lower = interval.upper;
    }
    let represented_energy = law.energy_scale_j * integrated_energy;
    nonnegative_finite(represented_energy, "mindlin_normal_transfer_energy_j")?;
    let scale = old.stored_energy_j.max(represented_energy);
    if (represented_energy - energy).abs() > 1024.0 * f64::EPSILON * scale {
        return Err(invalid("mindlin_normal_history_transfer_energy"));
    }
    let endpoint_force_n = law.force_scale_n * endpoint;
    finite(endpoint_force_n, "mindlin_normal_transfer_endpoint_n")?;
    Ok((
        MindlinState {
            owner: law,
            memory,
            slip_m: old.slip_m,
            endpoint_force_n,
            stored_energy_j: energy,
        },
        heat,
    ))
}

// Positive Bernstein moments integrate (|old strain|-new bound) times
// (|old strain|+new bound), avoiding subtraction of two nearly equal total
// elastic energies for a small normal decrement.
fn release_energy(
    old: &MindlinState,
    interval: &MemoryInterval,
    lower: f64,
    upper: f64,
    cutoff: f64,
    direction: f64,
) -> f64 {
    let a = (1.0 - lower).sqrt();
    let b = (1.0 - upper).sqrt();
    let w = (upper - lower) / (a + b);
    let denominator = (a + b) * (a + b);
    let variance = a * b / 3.0 + 4.0 * w * w / 45.0;
    let left = (variance + (a + 2.0 * b).powi(2) / 9.0) / denominator;
    let right = (variance + (2.0 * a + b).powi(2) / 9.0) / denominator;
    let mixed = (2.0 / 3.0) * (a * b + w * w / 5.0) / denominator;
    let values = [lower, upper].map(|r| {
        let old_abs = direction * (interval.slope * r + interval.offset);
        let bound = r - cutoff;
        // At an analytically computed branch crossing the difference is zero;
        // a negative last-bit residual has no physical negative dissipation.
        ((old_abs - bound).max(0.0), (old_abs + bound).max(0.0))
    });
    let integral = left * values[0].0 * values[0].1
        + mixed * (values[0].0 * values[1].1 + values[1].0 * values[0].1)
        + right * values[1].0 * values[1].1;
    0.5 * old.owner.energy_scale_j * w * integral
}
