//! Two independently evolving finite chambers under the original time driver.
//!
//! Total vapor and total specimen water can both stay constant while water
//! crosses from one chamber to the other. Each chamber's RH and vapor mass
//! therefore enter the temporal comparison separately. Only accepted native
//! half-step face transfers contribute to its whole-window inventory account.

use super::*;
use crate::moisture::{
    atmosphere::{
        AtmosphereMoistureError, FiniteMoistureAtmosphere,
        pair::{
            AtmosphereMoistureFace, MoistureAtmospherePort, PairedAtmosphereMoistureSolve,
            PairedAtmosphereMoistureStep,
        },
    },
    layered::{MoistureLayeredReport, MoistureLayeredSlab},
};

/// One original gas inventory over the complete requested duration.
#[derive(Clone, Debug)]
pub struct AtmosphereMoistureAdvanceAccount {
    pub relative_humidity: [f64; 2],
    pub vapor_mass_kg: [f64; 2],
    /// Signed sum of accepted transfers at this chamber's assigned faces.
    pub water_to_solids_kg: f64,
    pub water_balance_residual_kg: f64,
    pub maximum_carrier_mass_residual_kg: f64,
    pub total_mass_residual_kg: f64,
    /// Native gas m*u change at prescribed T, not sorption heat or heat duty.
    pub gas_internal_energy_change_j: f64,
}

/// Native whole-duration transport plus separately accountable chamber changes.
#[derive(Clone, Debug)]
pub struct PairedAtmosphereMoistureAdvanceReport<
    Transport = MoistureAdvanceReport<PairedAtmosphereMoistureStep, 2>,
> {
    /// Original field report, including each material layer for layered stacks.
    /// Its reservoir-water endpoints are the sum of both native vapor stores;
    /// no scalar RH is assigned to two independent chambers.
    pub transport: Transport,
    /// First and second, using the same port order as the native endpoint API.
    pub atmospheres: [AtmosphereMoistureAdvanceAccount; 2],
}

type PairState<F> = ([FiniteMoistureAtmosphere; 2], F);
type PairReport<R, const N: usize> = MoistureAdvanceReport<PairedAtmosphereMoistureStep<R, N>, N>;

impl<R: WaterStep<N>, const N: usize> WaterStep<N> for PairedAtmosphereMoistureStep<R, N> {
    fn boundary_water(&self) -> [f64; N] {
        self.transport.boundary_water()
    }
}

fn chamber_error<F>(
    old: &PairState<F>,
    coarse: &PairState<F>,
    fine: &PairState<F>,
    controls: MoistureTimeControls,
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<AtmosphereMoistureError>> {
    let mut ratio = 0.0_f64;
    for port in 0..2 {
        checkpoint(cx)?;
        let humidity = [
            old.0[port].relative_humidity(),
            coarse.0[port].relative_humidity(),
            fine.0[port].relative_humidity(),
        ];
        let [a, b, c] = humidity.map(|x| x.map_err(MoistureAdvanceError::FinalInventory));
        ratio = ratio.max(scaled_error(
            a?,
            b?,
            c?,
            controls.activity_absolute_tolerance,
            controls.relative_tolerance,
        )?);
        let water = [
            old.0[port].vapor_mass_kg(),
            coarse.0[port].vapor_mass_kg(),
            fine.0[port].vapor_mass_kg(),
        ];
        let [a, b, c] = water.map(|x| x.map_err(MoistureAdvanceError::FinalInventory));
        ratio = ratio.max(scaled_error(
            a?,
            b?,
            c?,
            controls.reservoir_water_absolute_tolerance_kg,
            controls.relative_tolerance,
        )?);
    }
    Ok(ratio)
}

fn chamber_account(
    original: &FiniteMoistureAtmosphere,
    next: &FiniteMoistureAtmosphere,
    received: f64,
    tolerance: f64,
) -> Result<AtmosphereMoistureAdvanceAccount, MoistureAdvanceError<AtmosphereMoistureError>> {
    let [water, carrier, mass] = original
        .exchange_residuals(next, received)
        .map_err(MoistureAdvanceError::FinalInventory)?;
    for residual in [water, carrier, mass] {
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(MoistureAdvanceError::FinalInventory(
                AtmosphereMoistureError::MassBalance {
                    residual_kg: residual,
                    tolerance_kg: tolerance,
                },
            ));
        }
    }
    let gas_internal_energy_change_j = next.gas().mass_kg()
        * next.gas().state().specific_internal_energy_j_kg()
        - original.gas().mass_kg() * original.gas().state().specific_internal_energy_j_kg();
    if !gas_internal_energy_change_j.is_finite() {
        return Err(MoistureAdvanceError::Invalid(
            "paired whole-window gas storage change is unrepresentable",
        ));
    }
    Ok(AtmosphereMoistureAdvanceAccount {
        relative_humidity: [
            original
                .relative_humidity()
                .map_err(MoistureAdvanceError::FinalInventory)?,
            next.relative_humidity()
                .map_err(MoistureAdvanceError::FinalInventory)?,
        ],
        vapor_mass_kg: [
            original
                .vapor_mass_kg()
                .map_err(MoistureAdvanceError::FinalInventory)?,
            next.vapor_mass_kg()
                .map_err(MoistureAdvanceError::FinalInventory)?,
        ],
        water_to_solids_kg: received,
        water_balance_residual_kg: water,
        maximum_carrier_mass_residual_kg: carrier,
        total_mass_residual_kg: mass,
        gas_internal_energy_change_j,
    })
}

fn finish<R, const N: usize>(
    original: [&FiniteMoistureAtmosphere; 2],
    next: &[FiniteMoistureAtmosphere; 2],
    report: &mut PairReport<R, N>,
    tolerance: f64,
    cx: &Cx<'_>,
) -> Result<[AtmosphereMoistureAdvanceAccount; 2], MoistureAdvanceError<AtmosphereMoistureError>> {
    let mut received = [WaterSum::default(); 2];
    for (side, assignment) in report.last_substep.face_atmospheres.iter().enumerate() {
        if let Some(port) = assignment {
            let index = usize::from(*port == MoistureAtmospherePort::Second);
            received[index].add(report.boundary_water_kg[side]);
        }
    }
    let accounts = [
        chamber_account(original[0], &next[0], received[0].total(), tolerance)?,
        chamber_account(original[1], &next[1], received[1].total(), tolerance)?,
    ];
    let mut vapor = [WaterSum::default(); 2];
    let mut water_balance = WaterSum::default();
    let mut mass_balance = WaterSum::default();
    water_balance.add(report.water_change_kg);
    mass_balance.add(report.water_change_kg);
    for (port, account) in accounts.iter().enumerate() {
        checkpoint(cx)?;
        for (sum, value) in vapor.iter_mut().zip(account.vapor_mass_kg) {
            sum.add(value);
        }
        water_balance.add(account.vapor_mass_kg[1] - account.vapor_mass_kg[0]);
        mass_balance.add(next[port].gas().mass_kg() - original[port].gas().mass_kg());
    }
    let water_balance = water_balance.total();
    let mass_balance = mass_balance.total();
    for residual in [water_balance, mass_balance] {
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(MoistureAdvanceError::FinalInventory(
                AtmosphereMoistureError::MassBalance {
                    residual_kg: residual,
                    tolerance_kg: tolerance,
                },
            ));
        }
    }
    let vapor = vapor.map(WaterSum::total);
    if vapor.iter().any(|m| !m.is_finite()) {
        return Err(MoistureAdvanceError::Invalid(
            "paired whole-window vapor inventory is unrepresentable",
        ));
    }
    report.reservoir_water_kg = Some(vapor);
    report.water_balance_residual_kg = water_balance;
    report.total_mass_residual_kg = Some(mass_balance);
    report.maximum_carrier_mass_residual_kg = Some(
        accounts[0]
            .maximum_carrier_mass_residual_kg
            .max(accounts[1].maximum_carrier_mass_residual_kg),
    );
    checkpoint(cx)?;
    Ok(accounts)
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<F: ReservoirField<N> + Clone, const N: usize>(
    original: [&FiniteMoistureAtmosphere; 2],
    field: &F,
    target_s: f64,
    time: MoistureTimeControls,
    mut propose: impl FnMut(
        &PairState<F>,
        f64,
    ) -> Result<
        (PairState<F>, PairedAtmosphereMoistureStep<F::Report, N>),
        AtmosphereMoistureError,
    >,
    cx: &Cx<'_>,
) -> Result<
    (
        [FiniteMoistureAtmosphere; 2],
        F,
        PairedAtmosphereMoistureAdvanceReport<PairReport<F::Report, N>>,
    ),
    MoistureAdvanceError<AtmosphereMoistureError>,
>
where
    F::Report: WaterStep<N>,
{
    checkpoint(cx)?;
    let initial = ([original[0].clone(), original[1].clone()], field.clone());
    let ((atmospheres, field), report) = evolve_fields_with_error(
        &initial,
        target_s,
        time,
        |state, dt| propose(state, dt).map(|(state, report)| (state, OneReport(report))),
        |state| field_snapshot::<F, N>(&state.1).map_err(AtmosphereMoistureError::from),
        |old, a, b| chamber_error(old, a, b, time, cx),
        cx,
    )?;
    let mut transport = report.summary.map_step(|report| report.0);
    let accounts = finish(
        original,
        &atmospheres,
        &mut transport,
        time.water_balance_tolerance_kg,
        cx,
    )?;
    Ok((
        atmospheres,
        field,
        PairedAtmosphereMoistureAdvanceReport {
            transport,
            atmospheres: accounts,
        },
    ))
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureAtmosphere {
    /// Complete slab diffusion between two evolving finite chambers to target.
    /// Each chamber's RH and vapor mass have independent temporal scales, even
    /// when opposite transfers leave total vapor and specimen water unchanged.
    /// Native face laws, fixed temperature/volume and physical root controls
    /// are unchanged. All accepted inputs survive any late refusal unchanged.
    pub fn advance_paired_adaptive_to(
        &self,
        second: &Self,
        moisture: &MoistureSlab,
        target_s: f64,
        faces: &[Option<AtmosphereMoistureFace>; 2],
        transport: MoistureStepConfig,
        physical: PairedAtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            [Self; 2],
            MoistureSlab,
            PairedAtmosphereMoistureAdvanceReport,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        advance(
            [self, second],
            moisture,
            target_s,
            time,
            |state, dt| {
                let (air, field, report) = state.0[0]
                    .advance_paired(&state.0[1], &state.1, dt, faces, transport, physical, cx)?
                    .into_parts();
                Ok(((air, field), report))
            },
            cx,
        )
    }

    /// Complete Cartesian diffusion with all six original chamber assignments
    /// and source-dependent principal diffusivities. Separate accepted face
    /// accounts survive aggregation into the two final gas inventories.
    pub fn advance_paired_cartesian_adaptive_to(
        &self,
        second: &Self,
        moisture: &MoistureCartesian,
        target_s: f64,
        faces: &[Option<AtmosphereMoistureFace>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: PairedAtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            [Self; 2],
            MoistureCartesian,
            PairedAtmosphereMoistureAdvanceReport<PairReport<MoistureCartesianReport, 6>>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        advance(
            [self, second],
            moisture,
            target_s,
            time,
            |state, dt| {
                let (air, field, report) = state.0[0]
                    .advance_paired_cartesian(
                        &state.0[1],
                        &state.1,
                        dt,
                        faces,
                        transport,
                        nonlinear,
                        physical,
                        cx,
                    )?
                    .into_parts();
                Ok(((air, field), report))
            },
            cx,
        )
    }

    /// Complete heterogeneous barrier diffusion and both finite chamber histories.
    /// Every material layer retains its original dry capacity, source activity,
    /// sorption, nonlinear diffusivity and finite interface law. Per-layer
    /// accuracy and signed interface accounts supplement the two chamber norms.
    /// Only the fully completed target returns the original native owners.
    pub fn advance_paired_layered_adaptive_to(
        &self,
        second: &Self,
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        faces: &[Option<MoistureAtmospherePort>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: PairedAtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            [Self; 2],
            MoistureLayeredSlab,
            PairedAtmosphereMoistureAdvanceReport<
                LayeredMoistureAdvanceReport<PairedAtmosphereMoistureStep<MoistureLayeredReport>>,
            >,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = ([self.clone(), second.clone()], moisture.clone());
        let ((atmospheres, field), report) = evolve_fields_with_error(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (air, field, report) = state.0[0]
                    .advance_paired_layered(
                        &state.0[1],
                        &state.1,
                        dt,
                        faces,
                        transport,
                        nonlinear,
                        physical,
                        cx,
                    )?
                    .into_parts();
                Ok(((air, field), layered::LayeredStep(report)))
            },
            |state| Ok(layered::observe(&state.1)),
            |old, a, b| chamber_error(old, a, b, time, cx),
            cx,
        )?;
        let mut transport = layered::finish(report, time.water_balance_tolerance_kg, cx)?;
        let accounts = finish(
            [self, second],
            &atmospheres,
            &mut transport.summary,
            time.water_balance_tolerance_kg,
            cx,
        )?;
        Ok((
            atmospheres,
            field,
            PairedAtmosphereMoistureAdvanceReport {
                transport,
                atmospheres: accounts,
            },
        ))
    }
}
