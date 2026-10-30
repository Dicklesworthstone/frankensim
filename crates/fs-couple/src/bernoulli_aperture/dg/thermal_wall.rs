//! Finite-body wall heat on the original moving-aperture/DG transaction.
//!
//! Mechanics and the complete load march once with the initial material.
//! Original solved main-port flows supply optional Kelvin/Maxwell heat. One
//! retained caloric chart then pays the native wall's signed parameter work;
//! every endpoint trial starts from that same marched state and initial body.

/// Reciprocal finite gas and material-wall heat on the continuing acoustic step.
pub mod gas;

use super::{ApertureDg, DgApertureFrame, accept_acoustic_candidate_observed, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::ApertureDrive;
use crate::bernoulli_aperture::modal::pressure::MaterialPressureFrame;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::bernoulli_aperture::wall::material::LocalMaterialWall;
use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, feedback};
use crate::dg_port::finite_gas::{FiniteGasDgError, refresh_output};
use crate::dg_port::{
    DgPortLoad, DgPortTrial, DgPortWallSource, DgPortWallUpdate, sum_terminal_work,
};
use fs_matdb::MaterialCard;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::state_point::DENSITY_PROPERTY;
use fs_phs::linear_exact::DenseExponentialBudget;

/// Explicit finite-body root, work and original source-integration budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteWallThermalSolve {
    /// The bracket is signed work OUT of the body and INTO the wall [J].
    /// Every queried body endpoint must lie in the original material/chart
    /// support. Source and energy tolerances remain distinct.
    pub exchange: HeatExchangeSolve,
    /// Absolute native wall parameter-work allowance [J].
    pub maximum_absolute_parameter_work_j: f64,
    /// Optional original solid-only heat integral on the solved DG interval.
    /// None leaves all native loss exported from this accounting boundary.
    pub dissipation_heat: Option<DenseExponentialBudget>,
}

/// One accepted interval of the original body, mechanics and complete load.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWallApertureDgStep {
    /// Actual interval pressures, impulses and work, with the final material's
    /// endpoint pressure/storage. Its native loss includes wall heat; that loss
    /// and `exported_dissipation_j` describe different accounting boundaries.
    pub coupled: DgApertureFrame,
    pub body_state: EquilibriumPhaseState,
    pub body_storage_change_j: f64,
    /// Prescribed signed energy into the finite body [J].
    pub external_heat_j: f64,
    /// Independently integrated native main-wall Kelvin/Maxwell heat [J].
    pub wall_heat_j: f64,
    /// Original wall sources and their native state-transfer/work report.
    pub wall: DgPortWallUpdate,
    /// Original mechanical/acoustic loss less wall heat deposited in this body.
    /// Terminal exchange stays separate; no rational storage is fabricated.
    pub exported_dissipation_j: f64,
    /// Body storage plus wall work less external heat and wall heat [J].
    pub body_energy_residual_j: f64,
    /// Body plus complete moving-load storage and exported losses, less actual
    /// mechanical/pressure inputs and external heat, plus terminal exchange [J].
    pub combined_energy_residual_j: f64,
    /// Proposed work out of the body less the actual native wall work [J].
    pub source_residual_j: f64,
    pub body_trials: usize,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticRealizeError> {
    if cancelled() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn same_density(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a / b - 1.0).abs() <= 64.0 * f64::EPSILON
}

fn contains_inventory(available: f64, amount: f64) -> bool {
    positive(amount) && (amount <= available || amount / available - 1.0 <= 64.0 * f64::EPSILON)
}

fn phase_error(error: fs_material::phase::PhaseStateError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}

fn duct_error(error: fs_duct::DuctError) -> AcousticRealizeError {
    match error {
        fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
        error => AcousticRealizeError::Duct(error),
    }
}

fn admit_phase(state: EquilibriumPhaseState, density: f64) -> Result<(), AcousticRealizeError> {
    if !state.is_fully_solid() || !same_density(state.bulk_density_kg_m3(), density) {
        return Err(invalid(
            "finite wall feedback requires an all-solid fixed-density body",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn admit_wall(
    source: &LocalMaterialWall,
    length: f64,
    radius: f64,
    body: UniformEnthalpyStepInput<'_>,
    material: &MaterialCard,
    temperature_axis: &str,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 2], AcousticRealizeError> {
    let source = source.at_temperature_pinned(
        material,
        temperature_axis,
        body.initial.temperature_k(),
        body.initial.temperature_k(),
        cancelled,
    )?;
    let density = source
        .material()
        .property(DENSITY_PROPERTY)
        .ok_or_else(|| invalid("finite wall requires its original material density"))?
        .value_si();
    if !same_density(density, body.mass_kg / body.volume_m3) {
        return Err(invalid(
            "finite body density must match each original wall material",
        ));
    }
    let volume = core::f64::consts::TAU * radius * length * source.thickness_m();
    let mass = volume * density;
    if !positive(volume) || !positive(mass) {
        return Err(invalid("finite wall inventory must be representable"));
    }
    Ok([volume, mass])
}

#[allow(clippy::too_many_arguments)]
fn admit<A>(
    system: &ApertureDg<A>,
    body: UniformEnthalpyStepInput<'_>,
    material: &MaterialCard,
    temperature_axis: &str,
    controls: FiniteWallThermalSolve,
    prepare_heat: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AcousticRealizeError> {
    poll(cancelled)?;
    if system.upstream.is_some() {
        return Err(invalid(
            "finite wall feedback needs an explicit joint inventory for an upstream tract",
        ));
    }
    if !positive(body.mass_kg)
        || !positive(body.volume_m3)
        || !positive(body.surface_area_m2)
        || !positive(body.duration_s)
        || body.duration_s.to_bits() != system.load.propagator().dt().to_bits()
        || !body.internal_heat_j.is_finite()
        || body.curve.material_card_identity() != material.content_hash()
        || !controls.maximum_absolute_parameter_work_j.is_finite()
        || controls.maximum_absolute_parameter_work_j < 0.0
    {
        return Err(invalid(
            "finite wall feedback requires one original body/card, native clock and finite work budget",
        ));
    }
    feedback::admit::<core::convert::Infallible>(controls.exchange).map_err(|_| {
        invalid("finite wall feedback requires a valid signed work bracket and energy allowances")
    })?;
    if body
        .curve
        .state_at_specific_enthalpy(body.initial.specific_enthalpy_j_kg())
        .map_err(phase_error)?
        != body.initial
    {
        return Err(invalid(
            "finite wall body must retain its canonical original enthalpy state",
        ));
    }
    admit_phase(body.initial, body.mass_kg / body.volume_m3)?;
    let duct = system.load.duct();
    let radius = duct
        .uniform_circular_radius_m()
        .ok_or_else(|| invalid("finite wall feedback requires original fixed circular geometry"))?;
    let source = system.load.material_wall_source().ok_or_else(|| {
        invalid("finite wall feedback requires retained original material wall sources")
    })?;
    let [volume, mass] = match source {
        DgPortWallSource::Uniform(wall) => admit_wall(
            wall,
            duct.total_length(),
            radius,
            body,
            material,
            temperature_axis,
            cancelled,
        )?,
        DgPortWallSource::Regional(walls) => {
            if walls.is_empty() || walls.len() != duct.mesh().elements {
                return Err(invalid(
                    "finite wall regions must retain every original main-bore element",
                ));
            }
            let positions = duct.node_positions();
            let local = duct.mesh().order + 1;
            let mut total = [0.0; 2];
            for (element, wall) in walls.iter().enumerate() {
                poll(cancelled)?;
                let base = element * local;
                let length = positions[base + local - 1] - positions[base];
                let inventory = admit_wall(
                    wall,
                    length,
                    radius,
                    body,
                    material,
                    temperature_axis,
                    cancelled,
                )?;
                for i in 0..2 {
                    let next = total[i] + inventory[i];
                    if !next.is_finite() || next <= total[i] {
                        return Err(invalid("regional wall inventory cannot be accumulated"));
                    }
                    total[i] = next;
                }
            }
            total
        }
    };
    if !contains_inventory(body.volume_m3, volume) || !contains_inventory(body.mass_kg, mass) {
        return Err(invalid(
            "finite body cannot contain less material than its complete main-bore wall",
        ));
    }
    if prepare_heat {
        if let Some(budget) = controls.dissipation_heat {
            duct.prepare_compliant_wall_heat(system.load.propagator(), budget, cancelled)
                .map_err(duct_error)?;
        }
    }
    poll(cancelled)
}

/// Cold renderer admission does not advance either physical participant.
pub(crate) fn admit_initial(
    system: &ApertureDg,
    body: UniformEnthalpyStepInput<'_>,
    material: &MaterialCard,
    temperature_axis: &str,
    controls: FiniteWallThermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AcousticRealizeError> {
    admit(
        system,
        body,
        material,
        temperature_axis,
        controls,
        true,
        cancelled,
    )
}

fn interval_heat(
    trial: &DgPortTrial<'_>,
    controls: FiniteWallThermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, AcousticRealizeError> {
    poll(cancelled)?;
    let Some(budget) = controls.dissipation_heat else {
        return Ok(0.0);
    };
    let heat = trial.wall_dissipated_heat_j(budget, cancelled)?;
    let remaining = trial.main_dissipated_j() - heat;
    if !remaining.is_finite() || remaining < -controls.exchange.energy_tolerance_j {
        return Err(invalid(
            "main-wall heat exceeds its original main-bore native loss allowance",
        ));
    }
    poll(cancelled)?;
    Ok(heat)
}

struct Candidate {
    load: DgPortLoad,
    state: EquilibriumPhaseState,
    wall: DgPortWallUpdate,
    storage_change_j: f64,
    source_residual_j: f64,
    body_energy_residual_j: f64,
    trials: usize,
}

fn feedback_error(error: feedback::ScalarError<AcousticRealizeError>) -> AcousticRealizeError {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Cancelled => AcousticRealizeError::Cancelled,
        feedback::ScalarError::Invalid(what) => invalid(what),
        feedback::ScalarError::OutsideBracket { lower, upper } => AcousticRealizeError::Nonlinear(
            format!("finite wall work lies outside its bracket: residuals {lower} J, {upper} J"),
        ),
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => AcousticRealizeError::Nonlinear(format!(
            "finite wall work did not converge after {iterations} trials: residual {residual} J"
        )),
    }
}

/// Reuse the same original enthalpy/native-wall work solve in joint transactions.
/// The observer bounds nested material work before evaluating each endpoint.
#[allow(clippy::too_many_arguments)]
fn body_candidate(
    marched: &DgPortLoad,
    body: UniformEnthalpyStepInput<'_>,
    material: &MaterialCard,
    temperature_axis: &str,
    controls: FiniteWallThermalSolve,
    wall_heat_j: f64,
    observed: &mut impl FnMut() -> Result<(), AcousticRealizeError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Candidate, AcousticRealizeError> {
    let heat = body.internal_heat_j + wall_heat_j;
    if !heat.is_finite()
        || (wall_heat_j != 0.0 && heat == body.internal_heat_j)
        || (body.internal_heat_j != 0.0 && heat == wall_heat_j)
    {
        return Err(invalid(
            "external or native wall heat disappears in the finite-body input",
        ));
    }
    let h0 = body.initial.specific_enthalpy_j_kg();
    let density = body.mass_kg / body.volume_m3;
    feedback::solve_scalar(
        controls.exchange.heat_bracket_j,
        controls.exchange.source_tolerance_j,
        controls.exchange.maximum_iterations,
        |work, trials, cancelled| {
            observed()?;
            let net = heat - work;
            let dh = net / body.mass_kg;
            let h1 = h0 + dh;
            if ![net, dh, h1].iter().all(|x| x.is_finite())
                || (net != 0.0 && (dh == 0.0 || h1 == h0))
                || (work != 0.0 && net == heat)
                || (heat != 0.0 && net == -work)
            {
                return Err(invalid(
                    "finite wall net energy disappears in its conserved enthalpy coordinate",
                ));
            }
            let state = body
                .curve
                .state_at_specific_enthalpy(h1)
                .map_err(phase_error)?;
            admit_phase(state, density)?;
            // Interior chart knots also belong to this fixed-volume path:
            // equal endpoint density cannot conceal intermediate expansion.
            for knot in body.curve.knots() {
                poll(cancelled)?;
                if knot.specific_enthalpy_j_kg >= h0.min(h1)
                    && knot.specific_enthalpy_j_kg <= h0.max(h1)
                    && (knot.liquid_mass_fraction != 0.0
                        || !same_density(knot.bulk_density_kg_m3, density))
                {
                    return Err(invalid(
                        "finite wall heating crossed unsupported phase or density change",
                    ));
                }
            }
            let mut load = marched.clone();
            let wall = load.rebind_wall_temperature_controlled(
                material,
                temperature_axis,
                body.initial.temperature_k(),
                state.temperature_k(),
                controls.maximum_absolute_parameter_work_j,
                cancelled,
            )?;
            let stored = body.mass_kg * (state.specific_enthalpy_j_kg() - h0);
            let residual = work - wall.parameter_work_j;
            let body_residual = (stored + wall.parameter_work_j) - heat;
            if ![stored, residual, body_residual]
                .iter()
                .all(|x| x.is_finite())
            {
                return Err(invalid("finite wall energy account must remain finite"));
            }
            Ok((
                Candidate {
                    load,
                    state,
                    wall,
                    storage_change_j: stored,
                    source_residual_j: residual,
                    body_energy_residual_j: body_residual,
                    trials,
                },
                residual,
            ))
        },
        cancelled,
    )
    .map_err(feedback_error)
}

#[allow(clippy::too_many_arguments)]
fn finish(
    marched: &DgPortLoad,
    mut frame: DgApertureFrame,
    initial_energy_j: f64,
    body: UniformEnthalpyStepInput<'_>,
    material: &MaterialCard,
    temperature_axis: &str,
    controls: FiniteWallThermalSolve,
    wall_heat_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(DgPortLoad, FiniteWallApertureDgStep), AcousticRealizeError> {
    let mut candidate = body_candidate(
        marched, body, material, temperature_axis, controls, wall_heat_j,
        &mut || Ok(()), cancelled,
    )?;
    refresh_output::<core::convert::Infallible>(&mut candidate.load, &mut frame.duct, cancelled)
        .map_err(|error| match error {
            FiniteGasDgError::Cancelled => AcousticRealizeError::Cancelled,
            FiniteGasDgError::Acoustic(error) => error,
            error => {
                AcousticRealizeError::Nonlinear(format!("finite wall endpoint refresh: {error:?}"))
            }
        })?;
    frame.parameter_work_j = candidate.wall.parameter_work_j;
    frame.stored_energy_j = frame.aperture.stored_energy_j + candidate.load.stored_energy_j();
    frame.storage_change_j = frame.stored_energy_j - initial_energy_j;
    let exported_dissipation_j = frame.inferred_total_dissipation_j - wall_heat_j;
    let combined = candidate.storage_change_j
        + frame.storage_change_j
        + exported_dissipation_j
        + sum_terminal_work(&frame.terminals)
        - frame.upstream_work_j
        - frame.body_work_j
        - frame.aperture.mechanical_work_j
        - body.internal_heat_j;
    let acoustic_residual = frame.balance_residual_j();
    let tolerance = controls.exchange.energy_tolerance_j;
    if ![
        frame.stored_energy_j,
        frame.storage_change_j,
        exported_dissipation_j,
        combined,
        acoustic_residual,
    ]
    .iter()
    .all(|x| x.is_finite())
        || (wall_heat_j != 0.0 && exported_dissipation_j == frame.inferred_total_dissipation_j)
        || (wall_heat_j > 0.0 && exported_dissipation_j < -tolerance)
        || candidate.body_energy_residual_j.abs() > tolerance
        || acoustic_residual.abs() > tolerance
        || combined.abs() > tolerance
    {
        return Err(invalid(
            "finite body, moving aperture and complete acoustic load do not close their energy balance",
        ));
    }
    poll(cancelled)?;
    let report = FiniteWallApertureDgStep {
        coupled: frame,
        body_state: candidate.state,
        body_storage_change_j: candidate.storage_change_j,
        external_heat_j: body.internal_heat_j,
        wall_heat_j,
        wall: candidate.wall,
        exported_dissipation_j,
        body_energy_residual_j: candidate.body_energy_residual_j,
        combined_energy_residual_j: combined,
        source_residual_j: candidate.source_residual_j,
        body_trials: candidate.trials,
    };
    Ok((candidate.load, report))
}

impl ApertureDg {
    /// Advance the original scalar aperture, complete load and one finite wall
    /// body atomically. The body contains the complete main-bore wall inventory;
    /// finite branches retain independent wall sources. Gas remains fixed.
    ///
    /// The native wall source integral receives the ACTUAL solved inlet,
    /// terminal and branch-junction inputs, once before endpoint trials. Body
    /// enthalpy satisfies ΔH + W_wall = Q_external + H_wall. The native acoustic
    /// frame retains its original loss; the enclosing body account exports
    /// that loss less H_wall. No loss remainder is relabeled physical heat.
    ///
    /// This is a first-order frozen-material partition. The original all-solid
    /// caloric chart, card, source conditions, geometry and density must remain
    /// compatible; expansion, phase transfer and gas/body conduction require
    /// their own physical coupling. An upstream tract is explicitly refused.
    /// Errors and cancellation preserve every accepted participant and clock.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_wall_heat_controlled(
        &mut self,
        body_state: &mut EquilibriumPhaseState,
        body: UniformEnthalpyStepInput<'_>,
        material: &MaterialCard,
        temperature_axis: &str,
        drive: TubeDrive,
        force_n: f64,
        controls: FiniteWallThermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteWallApertureDgStep, AcousticRealizeError> {
        if *body_state != body.initial {
            return Err(invalid(
                "finite wall input must begin at the accepted body state",
            ));
        }
        admit(
            self,
            body,
            material,
            temperature_axis,
            controls,
            false,
            cancelled,
        )?;
        let initial_energy_j = self.stored_energy_j();
        let relation = self.load.port_response()?;
        let mechanical = self.aperture.preview_step_with_force(
            ApertureDrive {
                upstream_pressure_pa: drive.upstream_pressure_pa,
                incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                body_flow_m3_s: drive.body_flow_m3_s,
            },
            force_n,
        )?;
        let mut marched = self.load.clone();
        let (frame, heat) = accept_acoustic_candidate_observed(
            &mut marched,
            relation,
            mechanical.frame,
            [MaterialPressureFrame::default(); 2],
            drive,
            |trial| interval_heat(trial, controls, cancelled),
        )?;
        let (load, report) = finish(
            &marched,
            frame,
            initial_energy_j,
            body,
            material,
            temperature_axis,
            controls,
            heat,
            cancelled,
        )?;
        // This last source update validates before its first mutation. No
        // fallible operation or cancellation poll follows publication.
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())?;
        self.aperture.accept_frame(mechanical);
        self.load = load;
        *body_state = report.body_state;
        Ok(report)
    }
}

impl ApertureDg<LoadedAperture> {
    /// The same finite-wall transaction for original single, paired or modal
    /// mechanics. Independent force and spatial-pressure work remain in their
    /// existing owners and enter the coupled account exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_wall_heat_controlled(
        &mut self,
        body_state: &mut EquilibriumPhaseState,
        body: UniformEnthalpyStepInput<'_>,
        material: &MaterialCard,
        temperature_axis: &str,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        controls: FiniteWallThermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteWallApertureDgStep, AcousticRealizeError> {
        if *body_state != body.initial {
            return Err(invalid(
                "finite wall input must begin at the accepted body state",
            ));
        }
        admit(
            self,
            body,
            material,
            temperature_axis,
            controls,
            false,
            cancelled,
        )?;
        let initial_energy_j = self.stored_energy_j();
        let relation = self.load.port_response()?;
        let mechanical = self.aperture.preview(
            ApertureDrive {
                upstream_pressure_pa: drive.upstream_pressure_pa,
                incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                body_flow_m3_s: drive.body_flow_m3_s,
            },
            forces_n,
            pressures_pa,
        )?;
        let mut marched = self.load.clone();
        let (frame, heat) = accept_acoustic_candidate_observed(
            &mut marched,
            relation,
            mechanical.frame(),
            mechanical.spatial_pressure(),
            drive,
            |trial| interval_heat(trial, controls, cancelled),
        )?;
        let (load, report) = finish(
            &marched,
            frame,
            initial_energy_j,
            body,
            material,
            temperature_axis,
            controls,
            heat,
            cancelled,
        )?;
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())?;
        self.aperture.accept(mechanical);
        self.load = load;
        *body_state = report.body_state;
        Ok(report)
    }
}
