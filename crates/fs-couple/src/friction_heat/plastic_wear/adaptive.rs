//! Adaptive refinement of the original hardness/contact/wear/heat split.
//!
//! A full interval and two half intervals start from the same accepted history.
//! Only the two native half steps may be accepted; no pressure, indentation,
//! temperature or source state is extrapolated. The defect estimates local
//! first-order splitting error, not a global bound or material validation.
//!
//! Motion is a prescribed, constant ACTUAL velocity at each original grid cell.
//! Tangential body compliance is zero. Dividing a released displacement while
//! retaining a nonzero body compliance would change the physical problem, so
//! this adapter does not accept that input. The original normal boundary is held
//! throughout the span. Thermal events and material export remain explicit
//! endpoint operations; neither is consumed by this march.

use super::*;
use fs_thermal::enthalpy::backward_euler::BackwardEulerStep;

/// Componentwise physical scales for the full-step/two-half-step difference.
#[derive(Clone, Copy, Debug)]
pub struct PlasticSlidingMarchTolerance {
    pub cell_energy_j: f64,
    /// Absolute interface-temperature scale; an initial solver guess is not a
    /// physical temperature history from which to infer a relative scale.
    pub interface_temperature_k: f64,
    pub pressure_pa: f64,
    /// Applies separately to approach, plastic depth, each body's recession and
    /// each component of the original accumulated tangential slip history.
    pub length_m: f64,
    /// Applies separately to friction, plastic, boundary, parameter and heat work.
    pub work_j: f64,
    /// Relative scales use changes from the accepted interval origin, not large
    /// background pressures, indentation histories or formation energies.
    pub relative: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct PlasticSlidingMarchControls {
    pub initial_step_s: f64,
    /// Smallest attempted comparison interval; a final remainder may be smaller.
    /// Each accepted comparison publishes two native half steps.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Includes rejected comparisons, each containing three bounded native solves.
    pub maximum_attempts: usize,
    /// Bounds retained source/mechanics receipts and accepted native half steps.
    pub maximum_accepted_substeps: usize,
    pub tolerance: PlasticSlidingMarchTolerance,
    /// Independent whole-span contact, thermal and combined account gates [J].
    /// This bound is not multiplied by the number of accepted substeps.
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PlasticSlidingMarchDefect {
    pub cell_energy_j: f64,
    pub interface_temperature_k: f64,
    pub pressure_pa: f64,
    pub length_m: f64,
    pub work_j: f64,
    pub normalized: f64,
}

pub struct AcceptedPlasticSlidingHeatStep<P> {
    pub start_s: f64,
    pub end_s: f64,
    /// Original accepted source payload, mechanics, thermal and work receipts.
    pub report: PlasticSlidingHeatReport<P>,
}

pub struct PlasticSlidingMarchReport<P> {
    pub steps: Vec<AcceptedPlasticSlidingHeatStep<P>>,
    pub attempts: usize,
    pub rejected_attempts: usize,
    /// Includes coarse and rejected native solves, not just accepted work.
    pub thermal_trials: usize,
    pub maximum_accepted_normalized_defect: f64,
    pub frictional_work_j: f64,
    pub normal_plastic_work_j: f64,
    pub normal_plastic_heat_j: [f64; 2],
    pub nonthermal_plastic_work_j: f64,
    pub normal_work_j: f64,
    pub temperature_parameter_work_j: f64,
    pub wear_parameter_work_j: f64,
    pub other_heat_j: f64,
    /// Signed interface heat delivered to each original body, excluding direct
    /// normal plastic heat (reported separately above).
    pub body_heat_j: [f64; 2],
    pub contact_elastic_energy_change_j: f64,
    pub thermal_energy_change_j: f64,
    /// Accepted native contact identity, excluding tangential body storage.
    pub contact_energy_residual_j: f64,
    pub thermal_energy_residual_j: f64,
    /// Contact elastic + thermal storage + unconverted plastic work minus
    /// normal, parameter, friction and other-heat supplies. This is not a
    /// complete body thermodynamic or passivity certificate.
    pub combined_energy_residual_j: f64,
}

pub struct PlasticSlidingMarchCandidate<'a, 't, P> {
    contact: PlasticSlidingHeatContact<'a>,
    thermal: DepositionSimulation<'t>,
    report: PlasticSlidingMarchReport<P>,
}
impl<'a, 't, P> PlasticSlidingMarchCandidate<'a, 't, P> {
    pub const fn contact(&self) -> &PlasticSlidingHeatContact<'a> {
        &self.contact
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &PlasticSlidingMarchReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PlasticSlidingHeatContact<'a>,
        DepositionSimulation<'t>,
        PlasticSlidingMarchReport<P>,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

#[derive(Debug)]
pub enum PlasticSlidingMarchError<E = Infallible> {
    Invalid(&'static str),
    Step(PlasticSlidingHeatError<E>),
    AttemptBudget,
    AcceptedSubstepBudget,
    Unresolved {
        step_s: f64,
        defect: PlasticSlidingMarchDefect,
    },
    Energy {
        contact_j: f64,
        thermal_j: f64,
        combined_j: f64,
    },
}
impl<E> From<PlasticSlidingHeatError<E>> for PlasticSlidingMarchError<E> {
    fn from(error: PlasticSlidingHeatError<E>) -> Self {
        Self::Step(error)
    }
}
impl<E> From<DepositionError> for PlasticSlidingMarchError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Step(PlasticSlidingHeatError::Thermal(error))
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for PlasticSlidingMarchError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive hardness/sliding/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for PlasticSlidingMarchError<E> {}

impl<'a> PlasticSlidingHeatContact<'a> {
    /// Propose an entire constant-velocity span from immutable accepted owners.
    ///
    /// The callback has the same candidate-only contract as `step_candidate`.
    /// Only accepted half-step payloads are retained. Any refusal, including one
    /// after private accepted substeps, publishes no partial span. Native source
    /// domain, footprint, removal, parameter-work and nonlinear refusals propagate
    /// directly: refinement never hides missing physical coverage. Each native
    /// solve retains its own original cell/trial/Newton budgets and residual gates.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn march_prescribed_velocity_candidate<'t, P, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        velocity_m_s: &[[f64; 2]],
        boundary: PlasticSlidingHeatBoundary,
        step_controls: PlasticSlidingHeatControls,
        controls: PlasticSlidingMarchControls,
        interface_temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<PlasticSlidingHeatSource<P>, E>,
        cx: &Cx<'_>,
    ) -> Result<PlasticSlidingMarchCandidate<'a, 't, P>, PlasticSlidingMarchError<E>> {
        checkpoint(cx)?;
        self.admit_thermal(thermal)?;
        let t = controls.tolerance;
        if [
            controls.initial_step_s,
            controls.minimum_step_s,
            controls.maximum_step_s,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
            || controls.minimum_step_s > controls.initial_step_s
            || controls.initial_step_s > controls.maximum_step_s
            || controls.maximum_attempts == 0
            || controls.maximum_accepted_substeps < 2
            || [
                t.cell_energy_j,
                t.interface_temperature_k,
                t.pressure_pa,
                t.length_m,
                t.work_j,
                t.relative,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || !interface_temperature_guess_k.is_finite()
            || interface_temperature_guess_k <= 0.0
            || velocity_m_s.len() != self.state.tangential_states().len()
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
            || thermal.cells().len() > step_controls.maximum_cells
        {
            return Err(PlasticSlidingMarchError::Invalid(
                "finite prescribed velocities, ordered positive steps and bounded controls required",
            ));
        }
        // Admit the WHOLE interval before allocating trials or querying sources.
        // A pending endpoint event remains pending; an interior/start event refuses.
        let _ = BackwardEulerStep::new(thermal, target_s, || checkpoint(cx))?;
        let zero_energy = vec![0.0; thermal.cells().len()];
        let mut accepted = self.clone();
        let mut accepted_thermal = thermal.clone();
        let mut interface_k = interface_temperature_guess_k;
        let mut step_s = controls.initial_step_s;
        let mut steps = Vec::new();
        let mut attempts = 0usize;
        let mut rejected = 0usize;
        let mut trials = 0usize;
        let mut maximum_defect: f64 = 0.0;
        let mut totals = [Sum::default(); 13];
        let mut residuals = [Sum::default(); 3];
        while accepted.time_s < target_s {
            checkpoint(cx)?;
            if attempts >= controls.maximum_attempts {
                return Err(PlasticSlidingMarchError::AttemptBudget);
            }
            if steps.len() > controls.maximum_accepted_substeps - 2 {
                return Err(PlasticSlidingMarchError::AcceptedSubstepBudget);
            }
            let start_s = accepted.time_s;
            let remaining_s = target_s - start_s;
            let planned_s = step_s.min(remaining_s);
            let end_s = if planned_s == remaining_s {
                target_s
            } else {
                start_s + planned_s
            };
            let half_s = start_s + 0.5 * (end_s - start_s);
            if half_s <= start_s || half_s >= end_s {
                return Err(PlasticSlidingMarchError::Invalid(
                    "time refinement is unrepresentable",
                ));
            }
            attempts += 1;
            let run = |contact: &PlasticSlidingHeatContact<'a>,
                       state: &DepositionSimulation<'t>,
                       end: f64,
                       guess: f64| {
                let dt = end - contact.time_s;
                let mut motion = Vec::with_capacity(velocity_m_s.len());
                for (i, velocity) in velocity_m_s.iter().enumerate() {
                    if i % 64 == 0 {
                        checkpoint(cx)?;
                    }
                    let released_increment_m =
                        [finite(velocity[0] * dt)?, finite(velocity[1] * dt)?];
                    motion.push(PlasticSlidingInput {
                        released_increment_m,
                        compliance_m_per_n: 0.0,
                    });
                }
                contact.step_candidate(
                    state,
                    end,
                    &motion,
                    boundary,
                    step_controls,
                    &zero_energy,
                    guess,
                    &resolve,
                    cx,
                )
            };
            let coarse = run(&accepted, &accepted_thermal, end_s, interface_k)?;
            let first = run(&accepted, &accepted_thermal, half_s, interface_k)?;
            let fine = run(
                first.contact(),
                first.thermal(),
                end_s,
                first.report().heat.interface_temperature_k,
            )?;
            for count in [
                coarse.report().trials,
                first.report().trials,
                fine.report().trials,
            ] {
                trials = trials
                    .checked_add(count)
                    .ok_or(PlasticSlidingMarchError::AttemptBudget)?;
            }
            let defect = compare(&accepted, &accepted_thermal, &coarse, &first, &fine, t, cx)?;
            if defect.normalized > 1.0 {
                rejected += 1;
                let reduced_s = 0.5 * (end_s - start_s);
                if reduced_s < controls.minimum_step_s {
                    return Err(PlasticSlidingMarchError::Unresolved {
                        step_s: end_s - start_s,
                        defect,
                    });
                }
                step_s = reduced_s;
                continue;
            }
            maximum_defect = maximum_defect.max(defect.normalized);
            let (_, _, first_report) = first.into_parts();
            let (next, next_thermal, fine_report) = fine.into_parts();
            interface_k = fine_report.heat.interface_temperature_k;
            for (begin, end, report) in [
                (start_s, half_s, first_report),
                (half_s, end_s, fine_report),
            ] {
                let values = accounts(&report);
                for (sum, value) in totals.iter_mut().zip(values) {
                    sum.add(value)?;
                }
                // Incremental native storage avoids subtracting large baseline
                // caloric energies or cumulative indentation work.
                for value in [values[11], values[1], -values[5], -values[6], -values[7]] {
                    residuals[0].add(value)?;
                }
                for value in [values[12], -values[0], -values[2], -values[3], -values[8]] {
                    residuals[1].add(value)?;
                }
                for value in [
                    values[11], values[12], values[4], -values[5], -values[6], -values[7],
                    -values[0], -values[8],
                ] {
                    residuals[2].add(value)?;
                }
                steps.push(AcceptedPlasticSlidingHeatStep {
                    start_s: begin,
                    end_s: end,
                    report,
                });
            }
            accepted = next;
            accepted_thermal = next_thermal;
            if defect.normalized <= 0.25 {
                step_s = (2.0 * (end_s - start_s)).min(controls.maximum_step_s);
            } else {
                step_s = end_s - start_s;
            }
        }
        let [contact_j, thermal_j, combined_j] = residuals.map(|sum| sum.value);
        if [contact_j, thermal_j, combined_j]
            .iter()
            .any(|v| v.abs() > controls.energy_tolerance_j)
        {
            return Err(PlasticSlidingMarchError::Energy {
                contact_j,
                thermal_j,
                combined_j,
            });
        }
        let v = totals.map(|sum| sum.value);
        checkpoint(cx)?;
        Ok(PlasticSlidingMarchCandidate {
            contact: accepted,
            thermal: accepted_thermal,
            report: PlasticSlidingMarchReport {
                steps,
                attempts,
                rejected_attempts: rejected,
                thermal_trials: trials,
                maximum_accepted_normalized_defect: maximum_defect,
                frictional_work_j: v[0],
                normal_plastic_work_j: v[1],
                normal_plastic_heat_j: [v[2], v[3]],
                nonthermal_plastic_work_j: v[4],
                normal_work_j: v[5],
                temperature_parameter_work_j: v[6],
                wear_parameter_work_j: v[7],
                other_heat_j: v[8],
                body_heat_j: [v[9], v[10]],
                contact_elastic_energy_change_j: v[11],
                thermal_energy_change_j: v[12],
                contact_energy_residual_j: contact_j,
                thermal_energy_residual_j: thermal_j,
                combined_energy_residual_j: combined_j,
            },
        })
    }
}

fn accounts<P>(r: &PlasticSlidingHeatReport<P>) -> [f64; 13] {
    [
        r.frictional_work_j,
        r.normal_plastic_work_j,
        r.normal_plastic_heat_j[0],
        r.normal_plastic_heat_j[1],
        r.nonthermal_plastic_work_j,
        r.normal_work_j,
        r.temperature_parameter_work_j,
        r.wear_parameter_work_j,
        r.other_heat_j,
        r.heat.body_heat_j[0],
        r.heat.body_heat_j[1],
        r.mechanics.temperature_normal().elastic_energy_change_j
            + r.mechanics.wear_normal().elastic_energy_change_j,
        r.heat.thermal.stored_energy_change_j,
    ]
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), PlasticSlidingMarchError<E>> {
        let y = value - self.correction;
        let next = checked(self.value + y)?;
        self.correction = checked((next - self.value) - y)?;
        self.value = next;
        Ok(())
    }
}
fn checked<E>(value: f64) -> Result<f64, PlasticSlidingMarchError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PlasticSlidingMarchError::Invalid(
            "unrepresentable march coordinate or account",
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn compare<P, E>(
    old: &PlasticSlidingHeatContact<'_>,
    old_thermal: &DepositionSimulation<'_>,
    coarse: &PlasticSlidingHeatCandidate<'_, '_, P>,
    first: &PlasticSlidingHeatCandidate<'_, '_, P>,
    fine: &PlasticSlidingHeatCandidate<'_, '_, P>,
    tolerance: PlasticSlidingMarchTolerance,
    cx: &Cx<'_>,
) -> Result<PlasticSlidingMarchDefect, PlasticSlidingMarchError<E>> {
    let mut out = PlasticSlidingMarchDefect::default();
    let mut compare_value = |a: f64, b: f64, origin: f64, absolute: f64, relative: f64| {
        let error = checked((a - b).abs())?;
        let scale = checked(absolute + relative * (a - origin).abs().max((b - origin).abs()))?;
        let normalized = if error == 0.0 {
            0.0
        } else if scale == 0.0 {
            f64::INFINITY
        } else {
            error / scale
        };
        out.normalized = out.normalized.max(normalized);
        Ok::<_, PlasticSlidingMarchError<E>>(error)
    };
    for (i, (id, cell)) in old_thermal.cells().iter().enumerate() {
        if i % 64 == 0 {
            checkpoint(cx)?;
        }
        let a = &coarse.thermal().cells()[id];
        let b = &fine.thermal().cells()[id];
        out.cell_energy_j = out.cell_energy_j.max(compare_value(
            cell.mass_kg * (a.specific_h_j_kg - cell.specific_h_j_kg),
            cell.mass_kg * (b.specific_h_j_kg - cell.specific_h_j_kg),
            0.0,
            tolerance.cell_energy_j,
            tolerance.relative,
        )?);
    }
    out.interface_temperature_k = compare_value(
        coarse.report().heat.interface_temperature_k,
        fine.report().heat.interface_temperature_k,
        0.0,
        tolerance.interface_temperature_k,
        0.0,
    )?;
    let o = old.state.wear_state().contact_state().contact_state();
    let a = coarse
        .contact
        .state
        .wear_state()
        .contact_state()
        .contact_state();
    let b = fine
        .contact
        .state
        .wear_state()
        .contact_state()
        .contact_state();
    out.length_m = compare_value(
        a.approach_m(),
        b.approach_m(),
        o.approach_m(),
        tolerance.length_m,
        tolerance.relative,
    )?;
    for i in 0..o.pressure_pa().len() {
        if i % 64 == 0 {
            checkpoint(cx)?;
        }
        out.pressure_pa = out.pressure_pa.max(compare_value(
            a.pressure_pa()[i],
            b.pressure_pa()[i],
            o.pressure_pa()[i],
            tolerance.pressure_pa,
            tolerance.relative,
        )?);
        out.length_m = out.length_m.max(compare_value(
            a.plastic_depth_m()[i],
            b.plastic_depth_m()[i],
            o.plastic_depth_m()[i],
            tolerance.length_m,
            tolerance.relative,
        )?);
        for side in 0..2 {
            out.length_m = out.length_m.max(compare_value(
                coarse.contact.state.wear_state().recession_m()[side][i],
                fine.contact.state.wear_state().recession_m()[side][i],
                old.state.wear_state().recession_m()[side][i],
                tolerance.length_m,
                tolerance.relative,
            )?);
            out.length_m = out.length_m.max(compare_value(
                coarse.contact.state.tangential_states()[i].slip_m()[side],
                fine.contact.state.tangential_states()[i].slip_m()[side],
                old.state.tangential_states()[i].slip_m()[side],
                tolerance.length_m,
                tolerance.relative,
            )?);
        }
    }
    for ((a, b), c) in accounts(coarse.report())
        .into_iter()
        .zip(accounts(first.report()))
        .zip(accounts(fine.report()))
    {
        out.work_j = out.work_j.max(compare_value(
            a,
            checked(b + c)?,
            0.0,
            tolerance.work_j,
            tolerance.relative,
        )?);
    }
    checkpoint(cx)?;
    Ok(out)
}
