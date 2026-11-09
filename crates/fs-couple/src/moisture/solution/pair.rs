//! Native bound-water diffusion between two independently evolving finite baths.
//!
//! A first-bath volume root is solved inside each second-bath volume trial.
//! Every native field proposal starts from the same original specimen and both
//! original bath stocks. Only that bath's assigned signed face transfers enter
//! its balance; through-flow is retained even with zero net specimen uptake.
//! No intermediate bath or field is published.
//!
//! Both source-supported brackets must remain admissible along the nested probe
//! path. An unbracketed inner root can refuse even when a coupled root exists
//! elsewhere. Brackets are never widened and no solvent or moisture is clipped.
//! This retains the original isothermal, equal-pressure vapor-equilibrium
//! bridge. It adds neither immersed-liquid permeability nor solute permeation,
//! pressure-driven osmosis, thermal feedback or higher-order time accuracy.

use super::*;
use crate::moisture::{
    radial::MoistureRadial, reservoir::water_sum, surface::ResolvedMoistureFace,
};

/// Position of one independent original solution inventory in the pair.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoistureSolutionPort {
    First,
    Second,
}
impl MoistureSolutionPort {
    pub(in crate::moisture) const fn index(self) -> usize {
        match self {
            Self::First => 0,
            Self::Second => 1,
        }
    }
}

/// One face communicates with one explicitly named bath through the original
/// sourced vapor-path sorption and surface-transfer law. `None` seals a face.
#[derive(Clone, Debug, PartialEq)]
pub struct SolutionMoistureFace {
    pub solution: MoistureSolutionPort,
    pub sorption: DryBasisSorption,
}

/// Independent finite-bath roots and one cumulative coupled-work allowance.
#[derive(Clone, Copy, Debug)]
pub struct PairedSolutionMoistureSolve {
    /// First (inner), then second (outer) original solvent-volume controls.
    pub solutions: [SolutionMoistureSolve; 2],
    /// All native field proposals across both roots, in 1..=100_000.
    /// Counts before surface-source or field admission, including failed trials.
    pub maximum_evaluations: usize,
    /// One whole-system AND field/face water tolerance [kg], independently
    /// checked in addition to each bath's own solvent-balance tolerance.
    pub mass_tolerance_kg: f64,
}

/// Actual native field and two separately balanced original solution accounts.
#[derive(Clone, Debug)]
pub struct PairedSolutionMoistureStep<Report = MoistureStepReport, const FACES: usize = 2> {
    pub transport: Report,
    pub face_solutions: [Option<MoistureSolutionPort>; FACES],
    pub solutions: [SolutionMoistureExchange; 2],
    pub solid_water_kg: [f64; 2],
    /// Compensated sum of represented local solid-water changes.
    pub water_change_kg: f64,
    /// Local solid-water change minus the two assigned native face totals.
    pub field_water_balance_residual_kg: f64,
    /// Local solid-water changes plus actual represented solvent-mass changes.
    pub water_balance_residual_kg: f64,
    /// Independent balance from each represented cell and bath endpoint/start.
    pub direct_water_balance_residual_kg: f64,
    pub field_evaluations: usize,
}

/// Two finite solutions and the original field accepted atomically.
#[derive(Clone, Debug)]
pub struct PairedSolutionMoistureCandidate<
    Field = MoistureSlab,
    Report = MoistureStepReport,
    const FACES: usize = 2,
> {
    solutions: [FiniteMoistureSolution; 2],
    moisture: Field,
    step: PairedSolutionMoistureStep<Report, FACES>,
}
impl<Field, Report, const FACES: usize> PairedSolutionMoistureCandidate<Field, Report, FACES> {
    pub(in crate::moisture) fn from_accepted(
        solutions: [FiniteMoistureSolution; 2],
        moisture: Field,
        step: PairedSolutionMoistureStep<Report, FACES>,
    ) -> Self {
        Self {
            solutions,
            moisture,
            step,
        }
    }

    pub const fn solutions(&self) -> &[FiniteMoistureSolution; 2] {
        &self.solutions
    }
    pub const fn moisture(&self) -> &Field {
        &self.moisture
    }
    pub const fn step(&self) -> &PairedSolutionMoistureStep<Report, FACES> {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        [FiniteMoistureSolution; 2],
        Field,
        PairedSolutionMoistureStep<Report, FACES>,
    ) {
        (self.solutions, self.moisture, self.step)
    }
}

pub type PairedCartesianSolutionMoistureCandidate =
    PairedSolutionMoistureCandidate<MoistureCartesian, MoistureCartesianReport, 6>;
pub type PairedRadialSolutionMoistureCandidate = PairedSolutionMoistureCandidate<MoistureRadial>;

impl FiniteMoistureSolution {
    /// Advance the original slab between two separate finite solutions.
    /// The receiver is First; `second` is Second, even when their values match.
    /// Each retains its own nonvolatile stocks, optional phases and activity law.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired(
        &self,
        second: &Self,
        moisture: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<SolutionMoistureFace>; 2],
        transport: MoistureStepConfig,
        controls: PairedSolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<PairedSolutionMoistureCandidate, SolutionMoistureError> {
        advance_pair(
            [self, second],
            moisture,
            dt_s,
            faces,
            transport,
            controls,
            cx,
            |report| report.boundary_water_kg,
            || Ok(()),
        )
    }

    /// Bind each actual -x,+x,-y,+y,-z,+z face to its named finite bath.
    /// Original anisotropic/nonlinear diffusion and all six signed receipts stay
    /// with the Cartesian owner; no homogeneous replacement field is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired_cartesian(
        &self,
        second: &Self,
        moisture: &MoistureCartesian,
        dt_s: f64,
        faces: &[Option<SolutionMoistureFace>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: PairedSolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<PairedCartesianSolutionMoistureCandidate, SolutionMoistureError> {
        advance_pair(
            [self, second],
            moisture,
            dt_s,
            faces,
            (transport, nonlinear),
            controls,
            cx,
            |report| report.boundary_water_kg,
            || Ok(()),
        )
    }

    /// Bind inner/outer cylindrical or spherical faces to independent finite
    /// baths. Native shell volumes and curved resistances remain authoritative;
    /// a zero-area center cannot be exposed as a physical exchange surface.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired_radial(
        &self,
        second: &Self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<SolutionMoistureFace>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: PairedSolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<PairedRadialSolutionMoistureCandidate, SolutionMoistureError> {
        advance_pair(
            [self, second],
            moisture,
            dt_s,
            faces,
            (transport, nonlinear),
            controls,
            cx,
            |report| report.boundary_water_kg,
            || Ok(()),
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(in crate::moisture) fn advance_pair<F: ReservoirField<N>, const N: usize>(
    solutions: [&FiniteMoistureSolution; 2],
    field: &F,
    dt_s: f64,
    faces: &[Option<SolutionMoistureFace>; N],
    transport: F::Controls,
    controls: PairedSolutionMoistureSolve,
    cx: &Cx<'_>,
    face_water: impl Fn(&F::Report) -> [f64; N],
    charge_field: impl FnMut() -> Result<(), SolutionMoistureError>,
) -> Result<PairedSolutionMoistureCandidate<F, F::Report, N>, SolutionMoistureError> {
    poll(cx)?;
    let [first, second] = solutions;
    let temperature = first.solution().temperature_k();
    let assignments = std::array::from_fn(|i| faces[i].as_ref().map(|face| face.solution));
    let sealed = [
        !assignments.contains(&Some(MoistureSolutionPort::First)),
        !assignments.contains(&Some(MoistureSolutionPort::Second)),
    ];
    let ([first_next, second_next], (moisture, report), [first_step, second_step], evaluations) =
        advance_baths(
            solutions,
            field.inventory().temperature_kelvin(),
            dt_s,
            sealed,
            controls,
            cx,
            |activities| {
                let mut boundary = std::array::from_fn(|_| ResolvedMoistureFace::default());
                for (side, face) in faces.iter().enumerate() {
                    poll(cx)?;
                    if let Some(face) = face {
                        boundary[side] = face
                            .sorption
                            .resolve(temperature, activities[face.solution.index()])?;
                    }
                }
                let (next, report) = field.advance_resolved(dt_s, boundary, transport, cx)?;
                let transfers = face_water(&report);
                let received = std::array::from_fn(|port| {
                    water_sum(transfers.iter().zip(&assignments).filter_map(
                        |(&water, assignment)| match *assignment {
                            Some(p) if p.index() == port => Some(water),
                            _ => None,
                        },
                    ))
                });
                Ok(((next, report), received))
            },
            charge_field,
        )?;
    let old = field.inventory();
    let next = moisture.inventory();
    let water_change_kg = field.water_change_to(&moisture);
    let solid_water_kg = [
        water_sum(
            old.moisture
                .iter()
                .enumerate()
                .map(|(i, m)| m * old.dry_mass_at(i)),
        ),
        water_sum(
            next.moisture
                .iter()
                .enumerate()
                .map(|(i, m)| m * next.dry_mass_at(i)),
        ),
    ];
    let field_water_balance_residual_kg = water_sum([
        water_change_kg,
        -first_step.water_to_solid_kg,
        -second_step.water_to_solid_kg,
    ]);
    let bath_changes = [
        first_next.solvent_mass_kg(),
        -first.solvent_mass_kg(),
        second_next.solvent_mass_kg(),
        -second.solvent_mass_kg(),
    ];
    let water_balance_residual_kg = water_sum(bath_changes.into_iter().chain([water_change_kg]));
    let direct_water_balance_residual_kg = water_sum(
        bath_changes.into_iter().chain(
            old.moisture
                .iter()
                .zip(next.moisture)
                .enumerate()
                .flat_map(|(i, (&a, &b))| [b * next.dry_mass_at(i), -a * old.dry_mass_at(i)]),
        ),
    );
    for residual_kg in [
        field_water_balance_residual_kg,
        water_balance_residual_kg,
        direct_water_balance_residual_kg,
    ] {
        if !residual_kg.is_finite() || residual_kg.abs() > controls.mass_tolerance_kg {
            return Err(SolutionMoistureError::MassBalance {
                residual_kg,
                tolerance_kg: controls.mass_tolerance_kg,
            });
        }
    }
    if solid_water_kg.iter().any(|mass| !mass.is_finite()) {
        return Err(SolutionMoistureError::Invalid(
            "paired solution specimen water inventory is unrepresentable",
        ));
    }
    poll(cx)?;
    Ok(PairedSolutionMoistureCandidate {
        solutions: [first_next, second_next],
        moisture,
        step: PairedSolutionMoistureStep {
            transport: report,
            face_solutions: assignments,
            solutions: [first_step, second_step],
            solid_water_kg,
            water_change_kg,
            field_water_balance_residual_kg,
            water_balance_residual_kg,
            direct_water_balance_residual_kg,
            field_evaluations: evaluations,
        },
    })
}

/// Original nested solvent balance for any native field/solid transaction.
/// The callback returns each bath's actual assigned signed face transfer.
/// No trial may restart from a previously rejected field or phase inventory.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(in crate::moisture) fn advance_baths<T>(
    solutions: [&FiniteMoistureSolution; 2],
    field_temperature_k: f64,
    dt_s: f64,
    sealed: [bool; 2],
    controls: PairedSolutionMoistureSolve,
    cx: &Cx<'_>,
    mut propose: impl FnMut([f64; 2]) -> Result<(T, [f64; 2]), SolutionMoistureError>,
    mut charge_field: impl FnMut() -> Result<(), SolutionMoistureError>,
) -> Result<
    (
        [FiniteMoistureSolution; 2],
        T,
        [SolutionMoistureExchange; 2],
        usize,
    ),
    SolutionMoistureError,
> {
    poll(cx)?;
    if controls.maximum_evaluations == 0
        || controls.maximum_evaluations > 100_000
        || !positive(controls.mass_tolerance_kg)
        || !positive(dt_s)
    {
        return Err(SolutionMoistureError::Invalid(
            "paired solutions need positive duration and water tolerance, with bounded cumulative evaluations",
        ));
    }
    let [first, second] = solutions;
    let temperature = first.solution().temperature_k();
    if second.solution().temperature_k().to_bits() != temperature.to_bits()
        || field_temperature_k.to_bits() != temperature.to_bits()
        || first.solvent_molar_mass_kg_mol().to_bits()
            != second.solvent_molar_mass_kg_mol().to_bits()
    {
        return Err(SolutionMoistureError::Invalid(
            "paired solutions require the same declared water molar mass and prescribed field temperature",
        ));
    }
    if let (Some(a), Some(b)) = (
        first.electrolyte_solubility(),
        second.electrolyte_solubility(),
    ) && a.solvent() != b.solvent()
    {
        return Err(SolutionMoistureError::Invalid(
            "paired ionic solutions must identify the same original solvent species",
        ));
    }
    let mut evaluations = 0_usize;
    let (second_next, accepted, second_step) = second.advance_exchange(
        dt_s,
        sealed[1],
        controls.solutions[1],
        cx,
        |second_activity| {
            let (first_next, accepted, first_step) = first.advance_exchange(
                dt_s,
                sealed[0],
                controls.solutions[0],
                cx,
                |first_activity| {
                    poll(cx)?;
                    if evaluations == controls.maximum_evaluations {
                        return Err(SolutionMoistureError::EvaluationBudget {
                            evaluations,
                            maximum: controls.maximum_evaluations,
                        });
                    }
                    charge_field()?;
                    evaluations += 1;
                    let (value, received) = propose([first_activity, second_activity])?;
                    Ok(((value, received[1]), received[0]))
                },
            )?;
            let (value, received) = accepted;
            Ok(((first_next, value, first_step), received))
        },
    )?;
    let (first_next, value, first_step) = accepted;
    Ok((
        [first_next, second_next],
        value,
        [first_step, second_step],
        evaluations,
    ))
}
