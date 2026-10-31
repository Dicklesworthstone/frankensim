//! One finite solution root for any admitted collection of native specimens.

use super::*;
use crate::moisture::reservoir::water_sum;

/// One common solution debit, with all signed specimen transfers retained in
/// their original native reports. Solute amounts are unchanged, even when
/// opposing adsorption/desorption leaves the net solvent inventory unchanged.
#[derive(Clone, Debug)]
pub struct SharedSolutionMoistureStep {
    pub interval_s: [f64; 2],
    pub water_activity: [f64; 2],
    /// Diagnostic only; this is not an applied mechanical pressure.
    pub osmotic_pressure_pa: [f64; 2],
    pub solvent_mass_kg: [f64; 2],
    pub solid_water_kg: [f64; 2],
    pub water_to_solids_kg: f64,
    pub water_balance_residual_kg: f64,
    /// Independent balance using original/final represented stored inventories.
    pub direct_water_balance_residual_kg: f64,
    /// Interior scalar trials; initial/bracket-endpoint trials are additional.
    pub iterations: usize,
}

/// One actual finite solution and every original field admitted atomically.
/// Specimens are in ascending stable-ID order, independent of enumeration.
#[derive(Clone, Debug)]
pub struct SharedSolutionMoistureCandidate<Field = MoistureSlab, Report = MoistureStepReport> {
    reservoir: FiniteMoistureSolution,
    solids: Vec<SharedMoistureSolid<Field, Report>>,
    step: SharedSolutionMoistureStep,
}
impl<Field, Report> SharedSolutionMoistureCandidate<Field, Report> {
    pub const fn reservoir(&self) -> &FiniteMoistureSolution {
        &self.reservoir
    }
    pub fn solids(&self) -> &[SharedMoistureSolid<Field, Report>] {
        &self.solids
    }
    pub const fn step(&self) -> &SharedSolutionMoistureStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        FiniteMoistureSolution,
        Vec<SharedMoistureSolid<Field, Report>>,
        SharedSolutionMoistureStep,
    ) {
        (self.reservoir, self.solids, self.step)
    }
}

/// One finite aqueous inventory and complete, independently sourced 3D fields.
pub type SharedCartesianSolutionMoistureCandidate =
    SharedSolutionMoistureCandidate<MoistureCartesian, MoistureCartesianReport>;

struct Trial<Field, Report> {
    solution: AqueousSolution,
    solids: Vec<SharedMoistureSolid<Field, Report>>,
    water_to_solids_kg: f64,
    residual_kg: f64,
}

fn admit_count(count: usize, maximum: usize) -> Result<(), SolutionMoistureError> {
    if count == 0
        || maximum == 0
        || maximum > crate::moisture::MAX_MOISTURE_CELLS
        || count > maximum
    {
        return Err(SolutionMoistureError::Invalid(
            "a nonempty collection of solution-coupled specimens within the aggregate cell budget is required",
        ));
    }
    Ok(())
}

impl FiniteMoistureSolution {
    /// Simultaneously exchange water with independently sourced slabs through
    /// one common endpoint solution activity. Every trial restarts all original
    /// fields; a specimen's uptake can be supplied by another's desorption.
    ///
    /// Stable IDs must be unique, accepted clocks and temperatures must agree,
    /// and the total cell budget includes sealed specimens whose internal
    /// diffusion still advances. The original scalar and native solver work
    /// bounds apply to the entire coupled proposal and each trial respectively.
    /// All participants and the one reservoir publish together or not at all.
    pub fn advance_shared(
        &self,
        participants: &[MoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSolutionMoistureCandidate, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
        let fields: Vec<_> = participants
            .iter()
            .map(|body| ReservoirParticipant {
                id: body.id,
                field: body.slab,
                faces: body.faces,
                transport: body.transport,
            })
            .collect();
        self.advance_fields(&fields, dt_s, maximum_total_cells, controls, cx)
    }

    /// Simultaneous six-face exchange for genuine Cartesian fields. Different
    /// geometries, source diffusivities and local surface laws share one solvent
    /// inventory and one water-activity root, without sequential reservoir debits.
    pub fn advance_shared_cartesian(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedCartesianSolutionMoistureCandidate, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
        let fields: Vec<_> = participants
            .iter()
            .map(|body| ReservoirParticipant {
                id: body.id,
                field: body.moisture,
                faces: body.faces,
                transport: (body.transport, body.nonlinear),
            })
            .collect();
        self.advance_fields(&fields, dt_s, maximum_total_cells, controls, cx)
    }

    pub(super) fn advance_fields<F: ReservoirField<N>, const N: usize>(
        &self,
        participants: &[ReservoirParticipant<'_, F, N>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSolutionMoistureCandidate<F, F::Report>, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(SolutionMoistureError::Invalid(
                "shared-solution specimen IDs must be unique",
            ));
        }
        let time_s = ordered[0].field.inventory().time_s;
        let temperature = self.solution.temperature_k();
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            let inventory = body.field.inventory();
            cells = cells
                .checked_add(inventory.moisture.len())
                .filter(|count| *count <= maximum_total_cells)
                .ok_or(SolutionMoistureError::Invalid(
                    "shared-solution total cell budget exhausted",
                ))?;
            if inventory.time_s.to_bits() != time_s.to_bits()
                || inventory.reference.query().axes()["temperature"].to_bits()
                    != temperature.to_bits()
            {
                return Err(SolutionMoistureError::Invalid(
                    "shared solution specimens need the same accepted clock and solution temperature",
                ));
            }
        }
        let initial_volume = self.solution.solvent_volume_m3();
        let [mut lower, mut upper] = controls.solvent_volume_bracket_m3;
        if !positive(lower)
            || !positive(upper)
            || lower >= upper
            || initial_volume < lower
            || initial_volume > upper
            || !positive(controls.mass_tolerance_kg)
            || controls.maximum_iterations == 0
            || controls.maximum_iterations > 256
        {
            return Err(SolutionMoistureError::Invalid(
                "positive source-supported volume bracket containing the start, finite mass tolerance and bounded iterations required",
            ));
        }
        let evaluate = |volume: f64| -> Result<Trial<F, F::Report>, SolutionMoistureError> {
            poll(cx)?;
            let solution = self.solution.with_solvent_volume_m3(volume)?;
            let activity = solution.water_activity()?;
            let mut solids = Vec::with_capacity(ordered.len());
            for body in &ordered {
                poll(cx)?;
                let (moisture, transport) = body.field.advance_activity(
                    dt_s,
                    body.faces,
                    temperature,
                    activity,
                    body.transport,
                    cx,
                )?;
                solids.push(SharedMoistureSolid {
                    id: body.id,
                    moisture,
                    transport,
                });
            }
            // Each dry inventory appears once. Opposing signed changes are
            // compensated in stable ID order, not averaged into equilibrium.
            let water_to_solids_kg = water_sum(
                ordered
                    .iter()
                    .zip(&solids)
                    .map(|(body, next)| body.field.water_change_to(&next.moisture)),
            );
            let residual_kg =
                (volume - initial_volume) * self.solvent_density_kg_m3 + water_to_solids_kg;
            if !residual_kg.is_finite() || !water_to_solids_kg.is_finite() {
                return Err(SolutionMoistureError::Invalid(
                    "finite solution water balance is unrepresentable",
                ));
            }
            Ok(Trial {
                solution,
                solids,
                water_to_solids_kg,
                residual_kg,
            })
        };

        let mut iterations = 0;
        let start = evaluate(initial_volume)?;
        let accepted = if start.residual_kg.abs() <= controls.mass_tolerance_kg {
            start
        } else if ordered
            .iter()
            .all(|body| body.faces.iter().all(Option::is_none))
        {
            return Err(SolutionMoistureError::MassBalance {
                residual_kg: start.residual_kg,
                tolerance_kg: controls.mass_tolerance_kg,
            });
        } else if self
            .solution
            .solute_moles()
            .iter()
            .all(|&amount| amount == 0.0)
        {
            // Pure-solvent activity is exactly one at every positive volume.
            // Every native field is already correct; only one finite external
            // inventory moves by the sum of actual local exchanges.
            let volume = initial_volume - start.water_to_solids_kg / self.solvent_density_kg_m3;
            if !positive(volume) || volume < lower || volume > upper {
                return Err(SolutionMoistureError::Unbracketed {
                    residual_kg: [
                        (lower - initial_volume) * self.solvent_density_kg_m3
                            + start.water_to_solids_kg,
                        (upper - initial_volume) * self.solvent_density_kg_m3
                            + start.water_to_solids_kg,
                    ],
                });
            }
            Trial {
                solution: self.solution.with_solvent_volume_m3(volume)?,
                residual_kg: (volume - initial_volume) * self.solvent_density_kg_m3
                    + start.water_to_solids_kg,
                ..start
            }
        } else {
            drop(start);
            let low = evaluate(lower)?;
            let high = evaluate(upper)?;
            if low.residual_kg.abs() <= controls.mass_tolerance_kg {
                low
            } else if high.residual_kg.abs() <= controls.mass_tolerance_kg {
                high
            } else {
                if low.residual_kg.is_sign_positive() == high.residual_kg.is_sign_positive() {
                    return Err(SolutionMoistureError::Unbracketed {
                        residual_kg: [low.residual_kg, high.residual_kg],
                    });
                }
                let mut low_residual = low.residual_kg;
                let mut high_residual = high.residual_kg;
                // Root refinement retains signs, not two extra complete sets
                // of source fields and specimens throughout the solve.
                drop(low);
                drop(high);
                loop {
                    poll(cx)?;
                    let scale = low_residual.abs().max(high_residual.abs());
                    let low_weight = low_residual.abs() / scale;
                    let high_weight = high_residual.abs() / scale;
                    let fraction = (low_weight / (low_weight + high_weight)).clamp(0.05, 0.95);
                    let volume = lower + fraction * (upper - lower);
                    if volume <= lower || volume >= upper {
                        return Err(SolutionMoistureError::NotConverged {
                            residual_kg: low_residual.abs().min(high_residual.abs()),
                            iterations,
                        });
                    }
                    let trial = evaluate(volume)?;
                    iterations += 1;
                    if trial.residual_kg.abs() <= controls.mass_tolerance_kg {
                        break trial;
                    }
                    if iterations == controls.maximum_iterations {
                        return Err(SolutionMoistureError::NotConverged {
                            residual_kg: trial.residual_kg,
                            iterations,
                        });
                    }
                    if trial.residual_kg.is_sign_positive() == low_residual.is_sign_positive() {
                        lower = volume;
                        low_residual = trial.residual_kg;
                    } else {
                        upper = volume;
                        high_residual = trial.residual_kg;
                    }
                }
            }
        };

        let reservoir = Self::new(accepted.solution, self.solvent_molar_mass_kg_mol)?;
        let solvent_mass_kg = [self.solvent_mass_kg(), reservoir.solvent_mass_kg()];
        let residual = solvent_mass_kg[1] - solvent_mass_kg[0] + accepted.water_to_solids_kg;
        let mut originals = Vec::with_capacity(ordered.len());
        let mut endpoints = Vec::with_capacity(ordered.len());
        for (body, next) in ordered.iter().zip(&accepted.solids) {
            poll(cx)?;
            let old = body.field.inventory();
            let end = next.moisture.inventory();
            originals.push(water_sum(
                old.moisture.iter().map(|m| m * old.dry_mass_per_cell_kg),
            ));
            endpoints.push(water_sum(
                end.moisture.iter().map(|m| m * end.dry_mass_per_cell_kg),
            ));
        }
        let solid_water_kg = [
            water_sum(originals.iter().copied()),
            water_sum(endpoints.iter().copied()),
        ];
        let direct_residual = water_sum(
            [solvent_mass_kg[1], -solvent_mass_kg[0]].into_iter().chain(
                endpoints
                    .iter()
                    .zip(&originals)
                    .flat_map(|(&new, &old)| [new, -old]),
            ),
        );
        if solvent_mass_kg
            .iter()
            .chain(&solid_water_kg)
            .any(|m| !m.is_finite())
        {
            return Err(SolutionMoistureError::Invalid(
                "shared solution original/final water inventory is unrepresentable",
            ));
        }
        for residual_kg in [residual, direct_residual] {
            if !residual_kg.is_finite() || residual_kg.abs() > controls.mass_tolerance_kg {
                return Err(SolutionMoistureError::MassBalance {
                    residual_kg,
                    tolerance_kg: controls.mass_tolerance_kg,
                });
            }
        }
        let step = SharedSolutionMoistureStep {
            interval_s: [time_s, time_s + dt_s],
            water_activity: [
                self.solution.water_activity()?,
                reservoir.solution.water_activity()?,
            ],
            osmotic_pressure_pa: [
                self.solution.osmotic_pressure_pa()?,
                reservoir.solution.osmotic_pressure_pa()?,
            ],
            solvent_mass_kg,
            solid_water_kg,
            water_to_solids_kg: accepted.water_to_solids_kg,
            water_balance_residual_kg: residual,
            direct_water_balance_residual_kg: direct_residual,
            iterations,
        };
        poll(cx)?;
        Ok(SharedSolutionMoistureCandidate {
            reservoir,
            solids: accepted.solids,
            step,
        })
    }
}
