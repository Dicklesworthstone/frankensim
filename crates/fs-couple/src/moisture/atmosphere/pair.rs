//! One native moisture field between two independent, finite vapor inventories.
//!
//! Each exposed face names its chamber explicitly. The two original atmosphere
//! balances use only the signed native face transfers assigned to that chamber;
//! a zero net change of solid water can still transfer water between chambers.
//! Every nested trial restarts the original field. The cumulative evaluation
//! limit bounds all inner and outer work, and only a fully balanced endpoint is
//! published.
//!
//! The first chamber is solved inside each trial of the second chamber. Its
//! source-supported RH root must remain bracketed at every outer probe. This
//! bounded nested solve can refuse even when a coupled root exists outside those
//! admitted probe paths; it does not widen brackets, clip water, or guarantee
//! convergence for arbitrary nonlinear sorption laws.
//!
//! Both gases and the specimen share one prescribed temperature. Original gas
//! volumes and solid reference geometry stay fixed. Carrier permeation,
//! condensation and heat feedback are excluded. The native gas storage changes
//! in the reports are not total sorption heat or thermostat duty.

use super::{
    AtmosphereMoistureError, AtmosphereMoistureSolve, DryBasisSorption,
    FiniteMoistureAtmosphere, SharedAtmosphereMoistureStep, invalid,
};
use crate::moisture::{
    MoistureNonlinearConfig, MoistureSlab, MoistureStepConfig, MoistureStepReport,
    cartesian::{MoistureCartesian, MoistureCartesianReport},
    poll, positive,
    reservoir::{ReservoirField, water_sum},
    surface::ResolvedMoistureFace,
};
use fs_exec::Cx;

/// Position of one physical gas inventory in the accepted pair.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoistureAtmospherePort {
    First,
    Second,
}
impl MoistureAtmospherePort {
    const fn index(self) -> usize {
        match self {
            Self::First => 0,
            Self::Second => 1,
        }
    }
}

/// One face contacts exactly one finite atmosphere through its original source
/// law. An absent face binding is sealed; finite surface resistance, when
/// present in the sorption binding, remains part of the native diffusion solve.
#[derive(Clone, Debug, PartialEq)]
pub struct AtmosphereMoistureFace {
    pub atmosphere: MoistureAtmospherePort,
    pub sorption: DryBasisSorption,
}

/// Individual reservoir gates plus one cumulative work and conservation budget.
#[derive(Clone, Copy, Debug)]
pub struct PairedAtmosphereMoistureSolve {
    /// First (inner), then second (outer) original finite-atmosphere controls.
    pub atmospheres: [AtmosphereMoistureSolve; 2],
    /// Total field proposals across every nested trial, in 1..=100_000.
    /// A proposal counts before its source and native transport admission.
    pub maximum_evaluations: usize,
    /// One absolute whole-system water AND total-mass tolerance [kg].
    /// This is independently checked; per-chamber tolerances are not added.
    pub mass_tolerance_kg: f64,
}

/// Native field and separately debited gas accounts at one accepted endpoint.
#[derive(Clone, Debug)]
pub struct PairedAtmosphereMoistureStep<Report = MoistureStepReport, const FACES: usize = 2> {
    pub transport: Report,
    /// Same ordering as the native face-transfer array in the transport report.
    pub face_atmospheres: [Option<MoistureAtmospherePort>; FACES],
    /// First, second. Each water_to_solids_kg contains only that port's signed
    /// face transfers, including through-flow that does not change total storage.
    pub atmospheres: [SharedAtmosphereMoistureStep; 2],
    /// Represented endpoint minus accepted-start solid water [kg].
    pub water_change_kg: f64,
    /// Actual vapor changes in BOTH native gas owners plus solid-water change.
    pub water_balance_residual_kg: f64,
    /// Actual total gas-mass changes plus solid-water change [kg].
    pub total_mass_residual_kg: f64,
    pub field_evaluations: usize,
}

/// Two original gas owners and one original spatial field, accepted atomically.
#[derive(Clone, Debug)]
pub struct PairedAtmosphereMoistureCandidate<
    Field = MoistureSlab,
    Report = MoistureStepReport,
    const FACES: usize = 2,
> {
    atmospheres: [FiniteMoistureAtmosphere; 2],
    moisture: Field,
    step: PairedAtmosphereMoistureStep<Report, FACES>,
}
impl<Field, Report, const FACES: usize>
    PairedAtmosphereMoistureCandidate<Field, Report, FACES>
{
    pub const fn atmospheres(&self) -> &[FiniteMoistureAtmosphere; 2] {
        &self.atmospheres
    }
    pub const fn moisture(&self) -> &Field {
        &self.moisture
    }
    pub const fn step(&self) -> &PairedAtmosphereMoistureStep<Report, FACES> {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        [FiniteMoistureAtmosphere; 2],
        Field,
        PairedAtmosphereMoistureStep<Report, FACES>,
    ) {
        (self.atmospheres, self.moisture, self.step)
    }
}

pub type PairedCartesianAtmosphereMoistureCandidate =
    PairedAtmosphereMoistureCandidate<MoistureCartesian, MoistureCartesianReport, 6>;

impl FiniteMoistureAtmosphere {
    /// Diffuse through the original slab while its two faces may contact two
    /// separate finite chambers. The receiver is First; second is Second.
    /// The two arguments denote separate physical inventories even when equal.
    /// Neither accepted gas nor the field changes on any refusal or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired(
        &self,
        second: &Self,
        moisture: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<AtmosphereMoistureFace>; 2],
        transport: MoistureStepConfig,
        controls: PairedAtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<PairedAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        advance_pair(
            [self, second], moisture, dt_s, faces, transport, controls, cx,
            |report| report.boundary_water_kg,
        )
    }

    /// Advance an actual 3D field with explicit chamber assignment in native
    /// -x,+x,-y,+y,-z,+z face order. All six face accounts remain separate until
    /// each is assigned to its one original gas inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired_cartesian(
        &self,
        second: &Self,
        moisture: &MoistureCartesian,
        dt_s: f64,
        faces: &[Option<AtmosphereMoistureFace>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: PairedAtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<PairedCartesianAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        advance_pair(
            [self, second], moisture, dt_s, faces, (transport, nonlinear), controls, cx,
            |report| report.boundary_water_kg,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_pair<F: ReservoirField<N>, const N: usize>(
    atmospheres: [&FiniteMoistureAtmosphere; 2],
    field: &F,
    dt_s: f64,
    faces: &[Option<AtmosphereMoistureFace>; N],
    transport: F::Controls,
    controls: PairedAtmosphereMoistureSolve,
    cx: &Cx<'_>,
    face_water: impl Fn(&F::Report) -> [f64; N],
) -> Result<PairedAtmosphereMoistureCandidate<F, F::Report, N>, AtmosphereMoistureError> {
    poll(cx)?;
    let temperature = atmospheres[0].gas().state().state().temperature;
    if field.inventory().temperature_kelvin().to_bits()
        != temperature.to_bits()
    {
        return Err(invalid(
            "paired atmospheres and specimen need the same prescribed temperature",
        ));
    }
    let face_atmospheres =
        std::array::from_fn(|i| faces[i].as_ref().map(|face| face.atmosphere));
    advance_pair_with(
        atmospheres, field.inventory().time_s, face_atmospheres, controls, cx,
        |rh| {
            let mut boundary: [ResolvedMoistureFace; N] =
                std::array::from_fn(|_| ResolvedMoistureFace::default());
            for (side, face) in faces.iter().enumerate() {
                poll(cx)?;
                if let Some(face) = face {
                    boundary[side] =
                        face.sorption.resolve(temperature, rh[face.atmosphere.index()])?;
                }
            }
            let (next, report) = field.advance_resolved(dt_s, boundary, transport, cx)?;
            let transfers = face_water(&report);
            let end_time = next.inventory().time_s;
            Ok((next, report, transfers, end_time))
        },
        |next| field.water_change_to(next),
    )
}

/// One shared nested gas balance for original homogeneous or heterogeneous
/// field owners. Each proposal starts from its accepted field, retains every
/// native face transfer, and returns the actual accepted clock. The separate
/// water-change callback measures original-to-final storage, not summed flux.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(in crate::moisture) fn advance_pair_with<Field, Report, const N: usize>(
    atmospheres: [&FiniteMoistureAtmosphere; 2],
    time_s: f64,
    face_atmospheres: [Option<MoistureAtmospherePort>; N],
    controls: PairedAtmosphereMoistureSolve,
    cx: &Cx<'_>,
    mut propose: impl FnMut([f64; 2])
        -> Result<(Field, Report, [f64; N], f64), AtmosphereMoistureError>,
    water_change: impl Fn(&Field) -> f64,
) -> Result<PairedAtmosphereMoistureCandidate<Field, Report, N>, AtmosphereMoistureError> {
    poll(cx)?;
    if controls.maximum_evaluations == 0
        || controls.maximum_evaluations > 100_000
        || !positive(controls.mass_tolerance_kg)
    {
        return Err(invalid(
            "paired atmospheres need positive whole-system mass tolerance and bounded cumulative evaluations",
        ));
    }
    let [first, second] = atmospheres;
    let temperature = first.gas().state().state().temperature;
    if second.gas().state().state().temperature.to_bits() != temperature.to_bits()
    {
        return Err(invalid(
            "paired atmospheres and specimen need the same prescribed temperature",
        ));
    }
    let first_water = &first.catalog.thermodynamics().components()[first.water];
    let second_water = &second.catalog.thermodynamics().components()[second.water];
    if first_water.species() != second_water.species()
        || first_water.molar_mass().value().to_bits()
            != second_water.molar_mass().value().to_bits()
    {
        return Err(invalid(
            "paired atmospheres must identify the same water species and molar mass",
        ));
    }
    let sealed = [
        !face_atmospheres.contains(&Some(MoistureAtmospherePort::First)),
        !face_atmospheres.contains(&Some(MoistureAtmospherePort::Second)),
    ];
    let mut evaluations = 0usize;
    let (second_next, accepted, second_step) = second.advance_exchange(
        time_s, sealed[1], controls.atmospheres[1], cx, |second_rh| {
            let (first_next, accepted, first_step) = first.advance_exchange(
                time_s, sealed[0], controls.atmospheres[0], cx, |first_rh| {
                    poll(cx)?;
                    if evaluations == controls.maximum_evaluations {
                        return Err(AtmosphereMoistureError::EvaluationBudget {
                            evaluations,
                            maximum: controls.maximum_evaluations,
                        });
                    }
                    evaluations += 1;
                    let (next, report, transfers, end_time) = propose([first_rh, second_rh])?;
                    let received: [f64; 2] = std::array::from_fn(|port| {
                        water_sum(transfers.iter().zip(&face_atmospheres).filter_map(
                            |(&water, assignment)| match *assignment {
                                Some(p) if p.index() == port => Some(water),
                                _ => None,
                            },
                        ))
                    });
                    Ok(((next, report, received, end_time), received[0], end_time))
                },
            )?;
            let (next, report, received, end_time) = accepted;
            Ok(((first_next, next, report, first_step), received[1], end_time))
        },
    )?;
    let (first_next, moisture, report, first_step) = accepted;
    let water_change_kg = water_change(&moisture);
    let water_balance_residual_kg = water_sum([
        water_change_kg,
        first_next.vapor_mass_kg()? - first.vapor_mass_kg()?,
        second_next.vapor_mass_kg()? - second.vapor_mass_kg()?,
    ]);
    let total_mass_residual_kg = water_sum([
        water_change_kg,
        first_next.gas().mass_kg() - first.gas().mass_kg(),
        second_next.gas().mass_kg() - second.gas().mass_kg(),
    ]);
    for residual in [water_balance_residual_kg, total_mass_residual_kg] {
        if !residual.is_finite() || residual.abs() > controls.mass_tolerance_kg {
            return Err(AtmosphereMoistureError::MassBalance {
                residual_kg: residual,
                tolerance_kg: controls.mass_tolerance_kg,
            });
        }
    }
    poll(cx)?;
    Ok(PairedAtmosphereMoistureCandidate {
        atmospheres: [first_next, second_next],
        moisture,
        step: PairedAtmosphereMoistureStep {
            transport: report,
            face_atmospheres,
            atmospheres: [first_step, second_step],
            water_change_kg,
            water_balance_residual_kg,
            total_mass_residual_kg,
            field_evaluations: evaluations,
        },
    })
}
