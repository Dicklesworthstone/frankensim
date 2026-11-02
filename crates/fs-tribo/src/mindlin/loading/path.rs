//! Analytic simultaneous Hertz--Mindlin indentation and scalar shear paths.
//!
//! Both indentation and macroscopic shear are affine in time. On a fixed
//! reference patch, the MDR Coulomb bound is `B(r,t) = r-c0+(c0-c1)t` and
//! the imposed shear rate is constant. Each elastic spring either remains
//! elastic or reaches one moving bound; newly contacting springs enter with
//! zero strain, and departing springs reach zero capacity before removal.
//! Yield/entry/exit times are affine in the radial coordinate. Their force,
//! storage and positive plastic-work integrals therefore require only the
//! original analytic radial moments, not a spring grid or time substeps.
//!
//! This integrates the same decoupled MDR model as the parent module. It is
//! not a trajectory solver, an elastic-wave solution or a full 3-D friction
//! solution. The normal FORCE is generally not affine: Hertz gives P~d^(3/2).
//! Reference: Popov, Hess and Willert (2019), section 4.4, Eqs. 4.21--4.25,
//! <https://doi.org/10.1007/978-3-662-58709-6_4>.

use super::super::{MemoryInterval, MindlinDeresiewiczPoint, MindlinState};
use super::{HertzMindlinContact, HertzMindlinState, invalid};
use crate::{
    FrictionRegime, TangentialLaw, TangentialStep, TangentialStepInput, TriboError, checked_add,
    checked_mul, finite, nonnegative_finite, positive_finite,
};
use fs_math::det;

mod compliance;

/// Complete path response; shear work uses the time-average force
/// in `tangential`, and normal work is audited separately.
#[derive(Clone, Debug, PartialEq)]
pub struct MindlinPathStep {
    /// Candidate history and time-average shear effort; includes actual heat.
    pub tangential: TangentialStep<HertzMindlinState>,
    /// Signed external normal work [J], independent of the shear port.
    pub normal_work_j: f64,
    /// Conservative Hertz storage change [J], equal to normal work.
    pub normal_stored_change_j: f64,
}

/// The original constitutive refusal or bounded caller cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum MindlinPathError {
    /// Original input, numerical representation or memory/force refusal.
    Constitutive(TriboError),
    /// Caller cancelled before publishing the candidate history.
    Cancelled,
}
impl From<TriboError> for MindlinPathError {
    fn from(error: TriboError) -> Self {
        Self::Constitutive(error)
    }
}
impl core::fmt::Display for MindlinPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz--Mindlin path: {self:?}")
    }
}
impl std::error::Error for MindlinPathError {}

impl HertzMindlinContact {
    /// Advance one simultaneous straight path in indentation and scalar shear.
    ///
    /// `input.normal_n` is the END load; its Hertz indentation and the accepted
    /// indentation bound the affine normal path. The first slip component is
    /// prescribed macroscopic travel. Scalar body compliance must be zero;
    /// deriving/accepting body trajectories belongs to the mechanical owner.
    /// `maximum_force_n` bounds the absolute returned average shear effort, not
    /// an unexamined peak force. All original memory/thermal/input checks remain.
    ///
    /// Fixed-normal paths delegate exactly to the original fixed-load law.
    /// General paths take O(retained intervals) time and bounded memory. A
    /// cancellation poll surrounds each original interval and each bounded
    /// event partition. No partial history is published on error/cancellation.
    pub fn advance_prescribed_path(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinPathStep, MindlinPathError> {
        self.advance_path(state, input, true, cancelled)
    }

    fn advance_path(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        prescribed: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinPathStep, MindlinPathError> {
        poll(cancelled)?;
        self.check_owner(state)?;
        input.solve.validate()?;
        positive_finite(input.dt_s, "dt_s")?;
        nonnegative_finite(input.normal_n, "normal_force_n")?;
        finite(input.slip_increment_m[0], "slip_increment_m")?;
        if (prescribed && input.compliance_m_per_n != 0.0) || input.slip_increment_m[1] != 0.0 {
            return Err(invalid("mindlin_path_requires_prescribed_first_axis_motion").into());
        }
        nonnegative_finite(input.compliance_m_per_n, "compliance_m_per_n")?;
        if let Some(t) = input.contact_temperature_k {
            positive_finite(t, "temperature_kelvin")?;
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "contact_pressure_pa")?;
        }
        if input.normal_n == state.normal_force_n() {
            let tangential = self.advance(state, input)?;
            poll(cancelled)?;
            return Ok(MindlinPathStep {
                tangential,
                normal_work_j: 0.0,
                normal_stored_change_j: 0.0,
            });
        }
        if self.parameters.friction_coefficient == 0.0 {
            let normal = self.normal_step(state, input.normal_n)?;
            let tangential = self.advance(&normal.state, input)?;
            poll(cancelled)?;
            return Ok(MindlinPathStep {
                tangential,
                normal_work_j: normal.normal_work_j,
                normal_stored_change_j: normal.normal_stored_change_j,
            });
        }
        let prepared = PreparedPath::new(self, state, input.normal_n)?;
        let (increment_m, iterations) = prepared.solve_increment(state, input, cancelled)?;
        let mut memory = Vec::with_capacity(
            (state.history_intervals() + 3).min(self.parameters.maximum_history_intervals),
        );
        let Response {
            average_n: average,
            endpoint_n: endpoint,
            energy_j: energy,
            heat_j: heat,
            ..
        } = prepared.response(state, increment_m, Some(&mut memory), cancelled)?;
        if average.abs() > input.solve.maximum_force_n {
            return Err(TriboError::OutsideApplicability {
                field: "maximum_force_n",
                value: average.abs(),
                minimum: 0.0,
                maximum: input.solve.maximum_force_n,
            }
            .into());
        }
        let change = energy - self.stored_energy_j(state);
        let work = checked_mul(average, increment_m, "mindlin_path_shear_work_j")?;
        let residual = change + heat - work;
        finite(residual, "mindlin_path_closure_j")?;
        let normal_change = prepared.normal_energy_j - state.normal_energy_j;
        finite(normal_change, "mindlin_path_normal_work_j")?;
        let slip_m = checked_add(state.slip_m, increment_m, "mindlin_accumulated_slip_m")?;
        let compatibility = increment_m
            - (input.slip_increment_m[0]
                - checked_mul(input.compliance_m_per_n, average, "mindlin_body_response_m")?);
        finite(compatibility, "mindlin_compatibility_residual_m")?;
        let fixed = prepared.final_law.map(|owner| MindlinState {
            owner,
            memory,
            slip_m,
            endpoint_force_n: endpoint,
            stored_energy_j: energy,
        });
        let next = HertzMindlinState {
            owner: self.clone(),
            fixed,
            slip_m,
            indentation_m: prepared.indentation_m,
            normal_energy_j: prepared.normal_energy_j,
        };
        poll(cancelled)?;
        Ok(MindlinPathStep {
            tangential: TangentialStep {
                traction_n: [average, 0.0],
                slip_increment_m: [increment_m, 0.0],
                state: next,
                regime: if heat > 0.0 {
                    FrictionRegime::Sliding
                } else {
                    FrictionRegime::Sticking
                },
                stored_change_j: change,
                dissipated_j: heat,
                closure_residual_j: residual,
                compatibility_residual_m: [compatibility, 0.0],
                iterations,
            },
            normal_work_j: normal_change,
            normal_stored_change_j: normal_change,
        })
    }
}

// Normal geometry and scales are held fixed throughout the mechanical solve.
// Trial responses never allocate or modify history; only the accepted response
// receives an output memory, as in the original fixed-load constitutive owner.
struct PreparedPath {
    c0: f64,
    c1: f64,
    final_law: Option<MindlinDeresiewiczPoint>,
    indentation_m: f64,
    normal_energy_j: f64,
    force_scale_n: f64,
    energy_scale_j: f64,
    critical_slip_m: f64,
    stiffness_n_per_m: f64,
    capacity_n: f64,
}

struct Response {
    average_n: f64,
    endpoint_n: f64,
    energy_j: f64,
    heat_j: f64,
    derivative_n_per_m: f64,
}

impl PreparedPath {
    fn new(
        owner: &HertzMindlinContact,
        state: &HertzMindlinState,
        normal_n: f64,
    ) -> Result<Self, MindlinPathError> {
        let (final_law, indentation_m, normal_energy_j) = if normal_n == 0.0 {
            (None, 0.0, 0.0)
        } else {
            let (law, d, energy) = owner.normal_geometry(normal_n)?;
            (Some(law), d, energy)
        };
        let growing = normal_n > state.normal_force_n();
        let reference = if growing {
            final_law.as_ref().expect("positive growing contact")
        } else {
            &state
                .fixed
                .as_ref()
                .expect("positive shrinking contact")
                .owner
        };
        let small_load = normal_n.min(state.normal_force_n());
        let ratio = small_load / reference.parameters.normal_force_n;
        let cutoff = if ratio == 0.0 {
            1.0
        } else {
            -det::expm1((2.0 / 3.0) * det::ln(ratio))
        };
        if !cutoff.is_finite() || cutoff <= 0.0 || (small_load > 0.0 && cutoff >= 1.0) {
            return Err(invalid("mindlin_path_patch_change_below_resolution").into());
        }
        let (c0, c1) = if growing {
            (cutoff, 0.0)
        } else {
            (0.0, cutoff)
        };
        let (force_scale_n, energy_scale_j, critical_slip_m, stiffness_n_per_m, capacity_n) = (
            reference.force_scale_n,
            reference.energy_scale_j,
            reference.critical_slip_m,
            reference.stiffness_n_per_m,
            reference.capacity_n,
        );
        Ok(Self {
            c0,
            c1,
            final_law,
            indentation_m,
            normal_energy_j,
            force_scale_n,
            energy_scale_j,
            critical_slip_m,
            stiffness_n_per_m,
            capacity_n,
        })
    }

    fn response(
        &self,
        state: &HertzMindlinState,
        increment_m: f64,
        memory: Option<&mut Vec<MemoryInterval>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Response, MindlinPathError> {
        let du = increment_m / self.critical_slip_m;
        finite(du, "mindlin_path_normalized_increment")?;
        let mut path = Path {
            c0: self.c0,
            c1: self.c1,
            du,
            final_law: self.final_law.as_ref(),
            memory,
            average: 0.0,
            endpoint: 0.0,
            energy: 0.0,
            heat: 0.0,
            derivative: 0.0,
        };
        if self.c0 > 0.0 {
            path.interval(0.0, self.c0, [0.0; 2], true, false, cancelled)?;
        }
        if let Some(old) = state.fixed.as_ref() {
            let mut lower = self.c0;
            for interval in &old.memory {
                poll(cancelled)?;
                let upper = if interval.upper == 1.0 {
                    1.0
                } else {
                    self.c0 + (1.0 - self.c0) * interval.upper
                };
                if upper <= lower {
                    return Err(invalid("mindlin_path_history_below_resolution").into());
                }
                let strain = [
                    interval.slope,
                    (1.0 - self.c0) * interval.offset - interval.slope * self.c0,
                ];
                if lower < self.c1 {
                    path.interval(lower, upper.min(self.c1), strain, false, true, cancelled)?;
                }
                if upper > self.c1 {
                    path.interval(lower.max(self.c1), upper, strain, false, false, cancelled)?;
                }
                lower = upper;
            }
        }
        let energy_j = checked_mul(self.energy_scale_j, path.energy, "mindlin_path_storage_j")?;
        let heat_j = checked_mul(self.energy_scale_j, path.heat, "mindlin_path_heat_j")?;
        let derivative_n_per_m = checked_mul(
            self.stiffness_n_per_m,
            path.derivative,
            "mindlin_path_tangent_n_per_m",
        )?;
        nonnegative_finite(energy_j, "mindlin_path_storage_j")?;
        nonnegative_finite(heat_j, "mindlin_path_heat_j")?;
        nonnegative_finite(derivative_n_per_m, "mindlin_path_tangent_n_per_m")?;
        Ok(Response {
            average_n: checked_mul(self.force_scale_n, path.average, "mindlin_path_average_n")?,
            endpoint_n: checked_mul(self.force_scale_n, path.endpoint, "mindlin_path_endpoint_n")?,
            energy_j,
            heat_j,
            derivative_n_per_m,
        })
    }
}

struct Path<'a> {
    c0: f64,
    c1: f64,
    du: f64,
    final_law: Option<&'a MindlinDeresiewiczPoint>,
    memory: Option<&'a mut Vec<MemoryInterval>>,
    average: f64,
    endpoint: f64,
    energy: f64,
    heat: f64,
    derivative: f64,
}

impl Path<'_> {
    // Contact entry/exit, initial strain and final elastic trial are affine in r.
    fn local(&self, r: f64, strain: [f64; 2], new: bool, detached: bool) -> Local {
        let db = self.c0 - self.c1;
        let enter = if new { (self.c0 - r) / db } else { 0.0 };
        let exit = if detached { (r - self.c0) / -db } else { 1.0 };
        let duration = (exit - enter).max(0.0);
        let initial = if new { 0.0 } else { strain[0] * r + strain[1] };
        Local {
            duration,
            initial,
            bound_before: if new { 0.0 } else { r - self.c0 },
            bound_after: if detached { 0.0 } else { r - self.c1 },
            trial: initial + self.du * duration,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn interval(
        &mut self,
        lower: f64,
        upper: f64,
        strain: [f64; 2],
        new: bool,
        detached: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MindlinPathError> {
        if lower == upper {
            return Ok(());
        }
        let mut cuts = [lower, upper, upper, upper];
        let mut count = 2;
        let ends = [lower, upper].map(|r| self.local(r, strain, new, detached));
        for sign in [-1.0, 1.0].into_iter().filter(|_| !new) {
            let a = ends[0].trial - sign * ends[0].bound_after;
            let b = ends[1].trial - sign * ends[1].bound_after;
            if (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) {
                let root = lower + (upper - lower) * (a / (a - b));
                if lower < root && root < upper {
                    cuts[count] = root;
                    count += 1;
                }
            }
        }
        cuts[..count].sort_by(f64::total_cmp);
        for pair in cuts[..count].windows(2) {
            poll(cancelled)?;
            let [lo, hi] = [pair[0], pair[1]];
            if lo == hi {
                continue;
            }
            let middle = self.local(f64::midpoint(lo, hi), strain, new, detached);
            let db = self.c0 - self.c1;
            // Newly entering springs all share the same ratio of shear to
            // growing capacity; decide it without roundoff-dependent radial
            // crossings. The resulting elastic slope stays in [-1, 1].
            let sign = if new {
                if self.du.abs() > db {
                    self.du.signum()
                } else {
                    0.0
                }
            } else if middle.trial > middle.bound_after {
                1.0
            } else if middle.trial < -middle.bound_after {
                -1.0
            } else {
                0.0
            };
            let rate = sign * self.du - db;
            if sign != 0.0 {
                positive_finite(rate, "mindlin_path_plastic_rate")?;
            }
            let values = [lo, hi].map(|r| {
                let p = self.local(r, strain, new, detached);
                let elastic = if sign == 0.0 {
                    p.duration
                } else {
                    ((p.bound_before - sign * p.initial) / rate).clamp(0.0, p.duration)
                };
                let plastic = p.duration - elastic;
                let at_yield = if sign == 0.0 {
                    p.trial
                } else {
                    sign * (p.bound_after - db * plastic).max(0.0)
                };
                [
                    p.initial + at_yield,
                    elastic,
                    if sign == 0.0 {
                        0.0
                    } else {
                        sign * at_yield + p.bound_after
                    },
                    plastic,
                    if detached {
                        0.0
                    } else if sign == 0.0 {
                        p.trial
                    } else {
                        sign * p.bound_after
                    },
                ]
            });
            let moments = ProductMoments::new(lo, hi);
            let product = |i: usize, j: usize| {
                moments.integral([values[0][i], values[1][i]], [values[0][j], values[1][j]])
            };
            let sliding = 0.5 * product(2, 3);
            self.average += 0.5 * product(0, 1) + sign * sliding;
            // Before yield, d(strain)/d(du) is elapsed elastic time. After
            // yield the strain follows the prescribed normal bound, whose
            // derivative is zero. Continuity cancels the moving event term.
            self.derivative += 0.5 * product(1, 1);
            self.heat += if sign == 0.0 { 0.0 } else { rate * sliding };
            self.energy += 0.5 * product(4, 4);
            self.endpoint += moments.integral([values[0][4], values[1][4]], [1.0; 2]);
            if !detached && let Some(memory) = self.memory.as_deref_mut() {
                let (slope, offset) = if sign != 0.0 {
                    (sign, -sign * self.c1)
                } else if new {
                    (self.du / db, -self.du / db * self.c1)
                } else {
                    (strain[0], strain[1] + self.du)
                };
                let mapped_upper = if hi == 1.0 {
                    1.0
                } else {
                    (hi - self.c1) / (1.0 - self.c1)
                };
                if mapped_upper <= memory.last().map_or(0.0, |p| p.upper) {
                    return Err(invalid("mindlin_path_history_below_resolution").into());
                }
                self.final_law
                    .expect("retained contact has positive final load")
                    .append_memory(
                        memory,
                        MemoryInterval {
                            upper: mapped_upper,
                            slope,
                            offset: (offset + slope * self.c1) / (1.0 - self.c1),
                        },
                    )?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Local {
    duration: f64,
    initial: f64,
    bound_before: f64,
    bound_after: f64,
    trial: f64,
}

// Exact positive Bernstein moments of the product of two affine functions.
// Unlike subtracting accumulated work and storage, this resolves tiny heat.
struct ProductMoments {
    weight: f64,
    left: f64,
    right: f64,
    mixed: f64,
}
impl ProductMoments {
    fn new(lo: f64, hi: f64) -> Self {
        let a = (1.0 - lo).sqrt();
        let b = (1.0 - hi).sqrt();
        let weight = (hi - lo) / (a + b);
        let denominator = (a + b) * (a + b);
        let variance = a * b / 3.0 + 4.0 * weight * weight / 45.0;
        Self {
            weight,
            left: (variance + (a + 2.0 * b).powi(2) / 9.0) / denominator,
            right: (variance + (2.0 * a + b).powi(2) / 9.0) / denominator,
            mixed: (2.0 / 3.0) * (a * b + weight * weight / 5.0) / denominator,
        }
    }
    fn integral(&self, f: [f64; 2], g: [f64; 2]) -> f64 {
        self.weight
            * (self.left * f[0] * g[0]
                + self.mixed * (f[0] * g[1] + f[1] * g[0])
                + self.right * f[1] * g[1])
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MindlinPathError> {
    if cancelled() {
        Err(MindlinPathError::Cancelled)
    } else {
        Ok(())
    }
}
