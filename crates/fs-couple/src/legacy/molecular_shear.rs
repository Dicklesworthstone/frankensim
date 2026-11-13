//! Finite gas and finite wall inertia under sealed Maxwell molecular shear.
//!
//! Both original owners move tangentially to one plane, without wall spin or
//! normal motion. The gas retains its original finite mass, volume, species,
//! calorics and bulk momentum. The wall retains its native rigid momentum and
//! original enthalpy/phase chart. No species is emitted, captured or removed.
//!
//! Explicit Maxwell diffuse fractions give the SAME weighted incident masses
//! and source-caloric reflected heat as `ThermalAccommodation`. With their sum
//! C and reduced mass mu, slip decays by exp(-C/mu). This exact frozen-coefficient
//! mechanical update cannot acquire the stiff slip reversal of an explicit
//! impulse. Its relative kinetic loss heats each original thermal store by D/2;
//! this split follows the diffuse half-range energy moments. Specular reflection
//! supplies normal pressure but no tangential drag or shear heating.
//!
//! Temperatures and encounter rates are frozen for this interval. Endpoint
//! calorics are fully nonlinear, including constant-density fusion plateaus;
//! this is not an exact integration of changing temperatures or a continuum
//! boundary layer. Explicit duration/temperature limits and time refinement
//! remain necessary. A chart segment that changes density refuses: rigid
//! geometry supplies neither expansion work nor moving normal boundaries.
//!
//! Normal pressure is COMPLETE molecular traction, not an extra pressure term
//! to add again. Separate stationary normal supports confine gas and wall, and
//! a nonrotating wall guide absorbs the patch moment. The wall's normal support
//! acts at its COM; the separate guide cancels the molecular COM moment. Their
//! work is zero. Gas
//! spin, spatial gas motion, and unconstrained three-dimensional angular closure
//! are not inferred from this homogeneous, tangential bulk-velocity owner.

use crate::liquid_gas::evaporation::{EvaporatingLiquidGasError, momentum::MovingGasParcel};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use fs_material::sublimation::gas::{
    GasSurfaceError, accommodation::ReflectedHeatStep, momentum::MaxwellReflection,
};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve, GasHeatStep};
use fs_math::{dd::Dd, det};
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};
use fs_thermal::sublimation::molecular::momentum::maxwell_reflected_normal_pressure;

/// One material patch on a plane fixed in its normal direction. The wall may
/// translate anywhere in its tangent plane; no cover/vacate geometry is modeled.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularShearSurface {
    pub area_m2: f64,
    /// Unit world normal directed from the wall into the gas.
    pub outward_normal: [f64; 3],
    /// Material patch centroid relative to the wall COM, in body coordinates.
    /// It translates with the nonrotating wall, retaining one constant world arm.
    pub centroid_arm_body_m: [f64; 3],
}

/// The existing rigid and caloric states bound to the same physical wall mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularShearWall<'a> {
    body: RigidBodyState,
    properties: MassProperties,
    thermal_mass_kg: f64,
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    phase: EquilibriumPhaseState,
}
impl<'a> MolecularShearWall<'a> {
    pub fn new(
        body: RigidBodyState,
        properties: MassProperties,
        thermal_mass_kg: f64,
        chart: &'a EquilibriumEnthalpyPhaseCurve,
        phase: EquilibriumPhaseState,
    ) -> Result<Self, MolecularShearError> {
        if !positive(thermal_mass_kg)
            || thermal_mass_kg != properties.mass()
            || body.angular_momentum_body() != Vec3::ZERO
            || chart.identity() != phase.phase_curve_identity()
        {
            return Err(invalid(
                "same rigid/thermal mass, original chart and zero wall spin required",
            ));
        }
        body.center_of_mass_velocity_world(properties)
            .map_err(MolecularShearError::Mechanics)?;
        Ok(Self {
            body,
            properties,
            thermal_mass_kg,
            chart,
            phase,
        })
    }
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    pub const fn thermal_mass_kg(&self) -> f64 {
        self.thermal_mass_kg
    }
    pub const fn chart(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.chart
    }
    pub const fn phase(&self) -> EquilibriumPhaseState {
        self.phase
    }
}

/// Physical temporal limits, bounded source work and independent SI budgets.
#[derive(Clone, Copy, Debug)]
pub struct MolecularShearControls {
    pub gas_energy: GasEnergySolve,
    pub maximum_duration_s: f64,
    /// Maximum absolute endpoint temperature change for EACH original owner.
    /// This does not bound thermal-rate truncation error on an isothermal plateau.
    pub maximum_temperature_change_k: f64,
    pub maximum_species: usize,
    pub maximum_phase_knots: usize,
    pub momentum_absolute_n_s: f64,
    pub energy_absolute_j: f64,
    pub position_absolute_m: f64,
    /// Relative allowance on exchanged momentum, heat/loss and displacement;
    /// never on large common velocity, total energy, position or pressure load.
    pub relative_tolerance: f64,
}

#[derive(Debug)]
pub enum MolecularShearError {
    Invalid(&'static str),
    Reflection(GasSurfaceError),
    Gas(GasEnergyError),
    MovingGas(EvaporatingLiquidGasError),
    Phase(PhaseStateError),
    Mechanics(DynamicsError),
    Balance {
        account: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MolecularShearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite molecular shear refused: {self:?}")
    }
}
impl std::error::Error for MolecularShearError {}

/// Original receipts and separately tested actual finite-state increments.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularShearStep {
    pub duration_s: f64,
    /// Original zero-drift reflected heat Q0 and per-species diffuse encounters.
    /// The net gas heat below additionally includes half the relative kinetic loss.
    pub reflection: ReflectedHeatStep,
    pub diffuse_encounter_mass_kg: f64,
    pub reduced_mass_kg: f64,
    /// 1-exp(-C/mu), in [0,1], with a stable small-exposure limit.
    pub relaxed_slip_fraction: f64,
    pub slip_before_m_s: [f64; 3],
    pub slip_after_m_s: [f64; 3],
    pub tangential_impulse_on_wall_n_s: [f64; 3],
    pub normal_molecular_impulse_n_s: f64,
    pub molecular_impulse_on_wall_n_s: [f64; 3],
    pub molecular_impulse_on_gas_n_s: [f64; 3],
    pub normal_support_impulse_on_wall_n_s: [f64; 3],
    pub normal_support_impulse_on_gas_n_s: [f64; 3],
    pub molecular_angular_impulse_on_wall_n_m_s: [f64; 3],
    pub guide_angular_impulse_on_wall_n_m_s: [f64; 3],
    /// Exactly zero under the admitted normal and rotational constraints.
    pub support_work_j: f64,
    /// Native receipt for NET molecular plus guide impulse. Pose drift is
    /// integrated separately over the same interval in the returned wall.
    pub wall_impulse: MomentumImpulseReceipt,
    pub gas_momentum_increment_n_s: [f64; 3],
    pub momentum_residual_n_s: [f64; 3],
    pub gas_kinetic_energy_change_j: f64,
    pub wall_kinetic_energy_change_j: f64,
    pub dissipated_relative_kinetic_energy_j: f64,
    pub mechanical_energy_residual_j: f64,
    pub gas_heat: GasHeatStep,
    pub requested_wall_heat_j: f64,
    pub stored_wall_heat_j: f64,
    pub wall_heat_residual_j: f64,
    /// Sum of original thermal storage minus D, separately bounded on D even
    /// when a much larger Q0 passes between the two stores.
    pub thermal_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    /// Exact frozen-coefficient COM displacement before position rounding.
    pub wall_displacement_m: [f64; 3],
    pub wall_position_residual_m: [f64; 3],
}

/// Both complete endpoints escape together; every refusal leaves inputs intact.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularShearCandidate<'a> {
    gas: MovingGasParcel,
    wall: MolecularShearWall<'a>,
    step: MolecularShearStep,
}
impl<'a> MolecularShearCandidate<'a> {
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn wall(&self) -> &MolecularShearWall<'a> {
        &self.wall
    }
    pub const fn step(&self) -> &MolecularShearStep {
        &self.step
    }
    pub fn into_parts(self) -> (MovingGasParcel, MolecularShearWall<'a>, MolecularShearStep) {
        (self.gas, self.wall, self.step)
    }
}

/// One sealed, normally confined interval using original gas and wall calorics.
/// There is no normal-velocity projection: nonzero normal velocity or spin
/// refuses. The single species exposure report drives both heat and shear;
/// only one original gas heat reconstruction is attempted. Allocation is the
/// bounded original species report and one original gas candidate, not a new
/// gas model, trial history or alternate source registry.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'a>(
    gas: &MovingGasParcel,
    wall: MolecularShearWall<'a>,
    surface: MolecularShearSurface,
    reflection: &MaxwellReflection,
    duration_s: f64,
    controls: MolecularShearControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MolecularShearCandidate<'a>, MolecularShearError> {
    poll(&mut cancelled)?;
    controls
        .gas_energy
        .validate()
        .map_err(MolecularShearError::Gas)?;
    if !positive(duration_s)
        || !positive(controls.maximum_duration_s)
        || duration_s > controls.maximum_duration_s
        || !positive(controls.maximum_temperature_change_k)
        || !positive(controls.momentum_absolute_n_s)
        || !positive(controls.energy_absolute_j)
        || !positive(controls.position_absolute_m)
        || !controls.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_tolerance)
        || controls.maximum_species == 0
        || controls.maximum_phase_knots < 2
        || gas.parcel().model().thermodynamics().components().len() > controls.maximum_species
        || wall.chart.knots().len() > controls.maximum_phase_knots
        || !nonnegative(surface.area_m2)
    {
        return Err(invalid(
            "finite physical interval, source-work bounds and positive SI budgets required",
        ));
    }
    let normal = vector(surface.outward_normal);
    let arm_body = vector(surface.centroid_arm_body_m);
    let gas_before = vector(gas.velocity_m_s());
    let wall_before = wall
        .body
        .center_of_mass_velocity_world(wall.properties)
        .map_err(MolecularShearError::Mechanics)?;
    if !normal.is_finite()
        || !arm_body.is_finite()
        || (normal.norm_squared() - 1.0).abs() > 64.0 * f64::EPSILON
        || !exactly_orthogonal(normal, gas_before)
        || !exactly_orthogonal(normal, wall_before)
    {
        return Err(invalid(
            "unit surface normal, finite arm and exactly tangential velocities required",
        ));
    }
    let gas_temperature = gas.parcel().state().state().temperature;
    let initial_pressure = gas.parcel().state().state().pressure;
    let wall_temperature = wall.phase.temperature_k();
    if !(controls.gas_energy.minimum_temperature_k..=controls.gas_energy.maximum_temperature_k)
        .contains(&gas_temperature)
    {
        return Err(invalid(
            "initial gas temperature outside declared energy bracket",
        ));
    }
    let thermal = reflection.thermal_accommodation();
    let reflected = thermal
        .evaluate_sealed(
            gas.parcel(),
            wall_temperature,
            surface.area_m2,
            duration_s,
            &mut cancelled,
        )
        .map_err(MolecularShearError::Reflection)?;
    let mut encounters = Dd::ZERO;
    for row in &reflected.species {
        poll(&mut cancelled)?;
        encounters = encounters + dd(row.accommodated_incident_mass_kg);
    }
    let encounters = encounters.to_f64();
    let mg = gas.parcel().mass_kg();
    let mw = wall.thermal_mass_kg;
    let reduced = mg.min(mw) / (1.0 + mg.min(mw) / mg.max(mw));
    if !nonnegative(encounters) || !positive(reduced) {
        return Err(invalid(
            "finite diffuse encounters and represented reduced mass required",
        ));
    }
    let gamma = encounters / reduced;
    if gamma.is_nan() || (encounters != 0.0 && gamma == 0.0) {
        return Err(invalid("molecular relaxation scale is unrepresentable"));
    }
    let relaxed = -det::expm1(-gamma);
    let slip = gas_before.sub(wall_before);
    let impulse = slip.scale(reduced * relaxed);
    let loss_root = det::sqrt(reduced)
        * core::f64::consts::FRAC_1_SQRT_2
        * det::sqrt(relaxed * (2.0 - relaxed))
        * norm(slip);
    let loss = loss_root * loss_root;
    if !slip.is_finite()
        || !impulse.is_finite()
        || !nonnegative(loss)
        || !(0.0..=1.0).contains(&relaxed)
        || (relaxed != 0.0 && norm(slip) != 0.0 && (norm(impulse) == 0.0 || loss == 0.0))
    {
        return Err(invalid(
            "molecular impulse or relative kinetic loss is unrepresentable",
        ));
    }
    poll(&mut cancelled)?;
    let wall_impulse = wall
        .body
        .apply_momentum_impulse(wall.properties, impulse, Vec3::ZERO)
        .map_err(MolecularShearError::Mechanics)?;
    let gas_momentum_before = vector(gas.momentum_kg_m_s());
    let gas_momentum_after = gas_momentum_before.sub(impulse);
    let gas_increment = gas_momentum_after.sub(gas_momentum_before);
    let wall_increment = wall_impulse.linear_increment_world;
    let momentum_residual = gas_increment.add(wall_increment);
    let impulse_scale = norm(impulse);
    for (account, residual) in [
        ("gas impulse", norm(gas_increment.add(impulse))),
        ("wall impulse", norm(wall_impulse.linear_residual_world)),
        ("joint momentum", norm(momentum_residual)),
    ] {
        gate(
            account,
            residual,
            controls.momentum_absolute_n_s,
            impulse_scale,
            controls,
        )?;
    }
    let gas_after = vector(array(gas_momentum_after).map(|p| p / mg));
    let wall_after = wall_impulse
        .state_after
        .center_of_mass_velocity_world(wall.properties)
        .map_err(MolecularShearError::Mechanics)?;
    if !gas_after.is_finite()
        || !exactly_orthogonal(normal, gas_after)
        || !exactly_orthogonal(normal, wall_after)
    {
        return Err(invalid(
            "represented endpoints must retain normal confinement",
        ));
    }
    let gas_kinetic = kinetic_increment(gas_before, gas_after, gas_increment);
    let wall_kinetic = kinetic_increment(wall_before, wall_after, wall_increment);
    let mechanical_residual = (gas_kinetic + wall_kinetic + dd(loss)).to_f64();
    gate(
        "mechanical loss",
        mechanical_residual,
        controls.energy_absolute_j,
        loss,
        controls,
    )?;

    let gas_heat_requested = reflected.heat_into_gas_j + 0.5 * loss;
    let wall_heat_requested = -reflected.heat_into_gas_j + 0.5 * loss;
    let dh = wall_heat_requested / mw;
    if !gas_heat_requested.is_finite()
        || !wall_heat_requested.is_finite()
        || !dh.is_finite()
        || (wall_heat_requested != 0.0 && dh == 0.0)
    {
        return Err(invalid("requested heat increments are unrepresentable"));
    }
    let phase = wall
        .chart
        .advance_specific_energy(wall.phase, dh)
        .map_err(MolecularShearError::Phase)?;
    admit_rigid_density(wall.chart, wall.phase, phase, &mut cancelled)?;
    let wall_stored = mw * (phase.specific_enthalpy_j_kg() - wall.phase.specific_enthalpy_j_kg());
    let wall_heat_residual = wall_stored - wall_heat_requested;
    if wall_heat_requested != 0.0 && wall_stored == 0.0 {
        return Err(invalid("wall heat was lost at its enthalpy datum"));
    }
    gate(
        "wall heat",
        wall_heat_residual,
        controls.energy_absolute_j,
        wall_heat_requested.abs(),
        controls,
    )?;
    let mut parcel = gas.parcel().clone();
    let gas_heat = parcel
        .add_heat(gas_heat_requested, controls.gas_energy, &mut cancelled)
        .map_err(MolecularShearError::Gas)?;
    if gas_heat_requested != 0.0 && gas_heat.stored_energy_change_j == 0.0 {
        return Err(invalid(
            "gas heat was not represented by its original calorics",
        ));
    }
    gate(
        "gas heat",
        gas_heat.residual_j,
        controls.energy_absolute_j,
        gas_heat_requested.abs(),
        controls,
    )?;
    let thermal_residual =
        (dd(gas_heat.stored_energy_change_j) + dd(wall_stored) - dd(loss)).to_f64();
    gate(
        "dissipation storage",
        thermal_residual,
        controls.energy_absolute_j,
        loss,
        controls,
    )?;
    let next_temperature = parcel.state().state().temperature;
    let [tmin, tmax] = thermal.temperature_range_k();
    if (next_temperature - gas_temperature).abs() > controls.maximum_temperature_change_k
        || (phase.temperature_k() - wall_temperature).abs() > controls.maximum_temperature_change_k
        || !(tmin..=tmax).contains(&next_temperature)
        || !(tmin..=tmax).contains(&phase.temperature_k())
    {
        return Err(invalid(
            "endpoint temperature exceeds interval or reflection source limits",
        ));
    }
    let total_residual =
        (gas_kinetic + wall_kinetic + dd(gas_heat.stored_energy_change_j) + dd(wall_stored))
            .to_f64();
    gate(
        "total energy",
        total_residual,
        controls.energy_absolute_j,
        loss.max(reflected.heat_into_gas_j.abs()),
        controls,
    )?;
    let gas = MovingGasParcel::from_momentum(parcel, array(gas_momentum_after))
        .map_err(MolecularShearError::MovingGas)?;

    // Integrate the same exponential wall velocity, including the gamma->0
    // limit. A midpoint drift would not reproduce the finite-relaxation path.
    let displacement = wall_before
        .add(slip.scale((reduced / mw) * mean_relaxation(gamma, relaxed)))
        .scale(duration_s);
    let position_before = wall.body.pose().position_world();
    let position_after = position_before.add(displacement);
    let position_residual = position_after.sub(position_before).sub(displacement);
    if !exactly_orthogonal(normal, displacement)
        || !exactly_orthogonal(normal, position_after.sub(position_before))
    {
        return Err(invalid(
            "represented wall displacement must retain its fixed plane",
        ));
    }
    gate(
        "wall position",
        norm(position_residual),
        controls.position_absolute_m,
        norm(displacement),
        controls,
    )?;
    let body = RigidBodyState::new(
        Pose::new(position_after, wall.body.pose().orientation())
            .map_err(MolecularShearError::Mechanics)?,
        wall_impulse.state_after.linear_momentum_world(),
        Vec3::ZERO,
    )
    .map_err(MolecularShearError::Mechanics)?;

    let mut normal_impulse = Dd::ZERO;
    let model = gas.parcel().model().thermodynamics();
    // Pressure and rates belong to the ORIGINAL snapshot, before either heat.
    for (component, &fraction) in model
        .components()
        .iter()
        .zip(model.mole_composition().fractions())
    {
        poll(&mut cancelled)?;
        let chi = *thermal
            .coefficients()
            .get(component.species())
            .ok_or(invalid("Maxwell reflection missing an active species"))?;
        let pressure = maxwell_reflected_normal_pressure(
            fraction * initial_pressure,
            wall_temperature,
            gas_temperature,
            0.0,
            chi,
        )
        .map_err(|error| MolecularShearError::Reflection(GasSurfaceError::Surface(error)))?;
        let load = pressure * (surface.area_m2 * duration_s);
        if !nonnegative(load) || (surface.area_m2 != 0.0 && pressure != 0.0 && load == 0.0) {
            return Err(invalid("normal molecular impulse is unrepresentable"));
        }
        normal_impulse = normal_impulse + dd(load);
    }
    let normal_impulse = normal_impulse.to_f64();
    let normal_load = normal.scale(normal_impulse);
    let molecular_wall = impulse.sub(normal_load);
    let arm_world = wall
        .body
        .pose()
        .orientation()
        .rotate_body_to_world(arm_body);
    let molecular_moment = arm_world.cross(molecular_wall);
    if !normal_load.is_finite() || !molecular_wall.is_finite() || !molecular_moment.is_finite() {
        return Err(invalid(
            "complete molecular load or guide moment is unrepresentable",
        ));
    }
    let step = MolecularShearStep {
        duration_s,
        reflection: reflected,
        diffuse_encounter_mass_kg: encounters,
        reduced_mass_kg: reduced,
        relaxed_slip_fraction: relaxed,
        slip_before_m_s: array(slip),
        slip_after_m_s: array(gas_after.sub(wall_after)),
        tangential_impulse_on_wall_n_s: array(impulse),
        normal_molecular_impulse_n_s: normal_impulse,
        molecular_impulse_on_wall_n_s: array(molecular_wall),
        molecular_impulse_on_gas_n_s: array(molecular_wall.scale(-1.0)),
        normal_support_impulse_on_wall_n_s: array(normal_load),
        normal_support_impulse_on_gas_n_s: array(normal_load.scale(-1.0)),
        molecular_angular_impulse_on_wall_n_m_s: array(molecular_moment),
        guide_angular_impulse_on_wall_n_m_s: array(molecular_moment.scale(-1.0)),
        support_work_j: 0.0,
        wall_impulse,
        gas_momentum_increment_n_s: array(gas_increment),
        momentum_residual_n_s: array(momentum_residual),
        gas_kinetic_energy_change_j: gas_kinetic.to_f64(),
        wall_kinetic_energy_change_j: wall_kinetic.to_f64(),
        dissipated_relative_kinetic_energy_j: loss,
        mechanical_energy_residual_j: mechanical_residual,
        gas_heat,
        requested_wall_heat_j: wall_heat_requested,
        stored_wall_heat_j: wall_stored,
        wall_heat_residual_j: wall_heat_residual,
        thermal_energy_residual_j: thermal_residual,
        total_energy_residual_j: total_residual,
        wall_displacement_m: array(displacement),
        wall_position_residual_m: array(position_residual),
    };
    poll(&mut cancelled)?;
    Ok(MolecularShearCandidate {
        gas,
        wall: MolecularShearWall {
            body,
            phase,
            ..wall
        },
        step,
    })
}

fn admit_rigid_density(
    chart: &EquilibriumEnthalpyPhaseCurve,
    before: EquilibriumPhaseState,
    after: EquilibriumPhaseState,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MolecularShearError> {
    let rho = before.bulk_density_kg_m3();
    if rho != after.bulk_density_kg_m3() {
        return Err(invalid(
            "rigid wall cannot accept a density-changing caloric step",
        ));
    }
    let lo = before
        .specific_enthalpy_j_kg()
        .min(after.specific_enthalpy_j_kg());
    let hi = before
        .specific_enthalpy_j_kg()
        .max(after.specific_enthalpy_j_kg());
    if lo != hi {
        for pair in chart.knots().windows(2) {
            poll(cancelled)?;
            if pair[0].specific_enthalpy_j_kg < hi
                && pair[1].specific_enthalpy_j_kg > lo
                && (pair[0].bulk_density_kg_m3 != rho || pair[1].bulk_density_kg_m3 != rho)
            {
                return Err(invalid(
                    "crossed wall source segment requires unsupported expansion work",
                ));
            }
        }
    }
    Ok(())
}

fn mean_relaxation(gamma: f64, relaxed: f64) -> f64 {
    if gamma < 1e-3 {
        // 1-(1-exp(-gamma))/gamma, with a finite, accurate zero limit.
        gamma
            * (0.5
                + gamma
                    * (-1.0 / 6.0 + gamma * (1.0 / 24.0 + gamma * (-1.0 / 120.0 + gamma / 720.0))))
    } else {
        1.0 - relaxed / gamma
    }
}
fn kinetic_increment(before: Vec3, after: Vec3, increment: Vec3) -> Dd {
    array(before)
        .into_iter()
        .zip(array(after))
        .zip(array(increment))
        .fold(Dd::ZERO, |sum, ((old, new), p)| {
            sum + (dd(0.5 * old) + dd(0.5 * new)) * dd(p)
        })
}
fn exactly_orthogonal(a: Vec3, b: Vec3) -> bool {
    let dot = array(a)
        .into_iter()
        .zip(array(b))
        .fold(Dd::ZERO, |sum, (x, y)| sum + dd(x) * dd(y));
    dot.hi == 0.0 && dot.lo == 0.0
}
fn gate(
    account: &'static str,
    residual: f64,
    absolute: f64,
    scale: f64,
    controls: MolecularShearControls,
) -> Result<(), MolecularShearError> {
    let allowance = absolute + controls.relative_tolerance * scale;
    if !residual.is_finite() || !allowance.is_finite() || residual.abs() > allowance {
        Err(MolecularShearError::Balance {
            account,
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MolecularShearError> {
    if cancelled() {
        Err(MolecularShearError::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(why: &'static str) -> MolecularShearError {
    MolecularShearError::Invalid(why)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn norm(a: Vec3) -> f64 {
    a.x.hypot(a.y).hypot(a.z)
}
