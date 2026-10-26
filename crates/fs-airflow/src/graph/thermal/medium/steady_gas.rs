//! Source-resolved convection and hydraulic losses on the shared steady gas field.
//!
//! NASA-9 enthalpy, species mixing and thermal transport remain owned by
//! `SteadyGasNetwork`, also used by `SteadyAxialGas`. This adapter closes wall conductances against
//! the actual local mixture and temperature, then evaluates pressure losses
//! from those same states. Flow is prescribed mass flow, never frozen volume
//! flow. Local pressure losses must be small relative to the common EOS pressure.

use fs_convection::{CorrelationId, HeatTransferCoefficient, NusseltEvaluation, ThermalDirection};
use fs_evidence::Evidence;
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas_field::steady::network::{
    GasNetworkCell, SteadyGasNetwork, SteadyGasNetworkSolution,
};
use fs_material::gas_field::steady::{
    MAX_STEADY_GAS_CELLS, SteadyAxialGas, SteadyGasError, SteadyGasSolution, SteadyGasSolve,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_qty::{Length, MassFlowRate};

use super::convection::{
    ForcedConvectionSettings, MaterialConvectionError, evaluate_forced_from_material,
};
use super::series::hydraulic::{MaterialFlowResistance, ResistanceError, ResistanceEvaluation};

/// Explicit per-cell correlation and hydraulic geometry. Material properties
/// are resolved from the same downstream gas state used by native donor-cell
/// transport; no independent density, heat capacity or volume flow is supplied.
#[derive(Clone, Copy, Debug)]
pub struct GasCellTransport {
    pub characteristic_length: Length,
    pub correlation: CorrelationId,
    pub settings: ForcedConvectionSettings,
    pub resistance: MaterialFlowResistance,
}

/// Bounds on the material/thermal fixed point, separate from native thermal work.
#[derive(Clone, Copy, Debug)]
pub struct GasTransportFeedback {
    /// At most 64 calls to the existing steady thermal/species solver.
    pub maximum_iterations: usize,
    /// Total native thermal substeps over ALL feedback iterations.
    pub maximum_substeps: usize,
    /// Fraction of a newly evaluated conductance applied at the next iteration.
    /// In (0, 1]; relaxation affects iteration, not the accepted physical law.
    pub conductance_relaxation: f64,
    /// Maximum total pressure loss/common EOS pressure, in (0, 0.1].
    /// This is an explicit low-pressure-drop model limit, not an error certificate.
    pub maximum_pressure_drop_fraction: f64,
}

/// The original material state and both transport responses at one accepted cell.
#[derive(Debug)]
pub struct GasCellTransportSample {
    pub cell: usize,
    pub gas: ResolvedThermochemicalGasState,
    pub transport: FluidTransportState,
    pub velocity_m_s: f64,
    pub nusselt: NusseltEvaluation,
    pub coefficient: Evidence<HeatTransferCoefficient>,
    pub conductance_w_k: f64,
    /// Signed heat received by the gas from its prescribed wall [W].
    pub wall_heat_w: f64,
    pub hydraulic: ResistanceEvaluation,
}

struct ThermalSample {
    gas: ResolvedThermochemicalGasState,
    transport: FluidTransportState,
    velocity_m_s: f64,
    nusselt: NusseltEvaluation,
    coefficient: Evidence<HeatTransferCoefficient>,
    conductance_w_k: f64,
    wall_heat_w: f64,
}

/// A closed shared field, with the actual laws and total work used to obtain it.
#[derive(Debug)]
pub struct MaterialSteadyGasSolution<T = SteadyGasSolution> {
    pub thermal: T,
    pub cells: Vec<GasCellTransportSample>,
    pub iterations: usize,
    pub substeps: usize,
    /// Native caloric/steady residual bound plus the accepted wall-law defect [W].
    pub maximum_power_residual_bound_w: f64,
    /// Sum of the nonnegative cell pressure-loss estimates [Pa]. This is the
    /// series pressure drop for an axial path. On a branched network it is a
    /// conservative sum over all cells, NOT a solved inlet/outlet pressure
    /// difference or a proof that the prescribed flows satisfy momentum balance.
    pub total_pressure_drop_pa: f64,
    /// Hydraulic dp*Q [W], reported once and NOT added to wall heat or enthalpy.
    /// This reduced thermal model neglects mechanical dissipation heating.
    pub hydraulic_dissipation_w: f64,
}

/// The same material feedback result with a general thermal-network solution.
pub type MaterialSteadyGasNetworkSolution = MaterialSteadyGasSolution<SteadyGasNetworkSolution>;

// Only a view of the two original result types. No state, material law or
// thermal balance is reconstructed by this common feedback adapter.
trait SolvedGas {
    fn substeps(&self) -> usize;
    fn maximum_power_residual_bound_w(&self) -> f64;
    fn mass_flow_kg_s(&self) -> &[f64];
    fn cell_gas(&self, cell: usize) -> Result<ResolvedThermochemicalGasState, MaterialError>;
}
impl SolvedGas for SteadyGasSolution {
    fn substeps(&self) -> usize {
        SteadyGasSolution::substeps(self)
    }
    fn maximum_power_residual_bound_w(&self) -> f64 {
        SteadyGasSolution::maximum_power_residual_bound_w(self)
    }
    fn mass_flow_kg_s(&self) -> &[f64] {
        SteadyGasSolution::mass_flow_kg_s(self)
    }
    fn cell_gas(&self, cell: usize) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        self.field().state_at(self.field().stations_m()[cell + 1])
    }
}
impl SolvedGas for SteadyGasNetworkSolution {
    fn substeps(&self) -> usize {
        SteadyGasNetworkSolution::substeps(self)
    }
    fn maximum_power_residual_bound_w(&self) -> f64 {
        SteadyGasNetworkSolution::maximum_power_residual_bound_w(self)
    }
    fn mass_flow_kg_s(&self) -> &[f64] {
        SteadyGasNetworkSolution::mass_flow_kg_s(self)
    }
    fn cell_gas(&self, cell: usize) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        Ok(self.states()[cell].clone())
    }
}

/// No failure publishes a partial field or changes the input description.
#[derive(Debug)]
pub enum MaterialSteadyGasError {
    Invalid(&'static str),
    Material(MaterialError),
    Convection(MaterialConvectionError),
    Resistance(ResistanceError),
    Thermal(SteadyGasError),
    NotConverged {
        iterations: usize,
        substeps: usize,
        residual_bound_w: f64,
    },
    PressureDrop {
        fraction: f64,
        maximum: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MaterialSteadyGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "steady material gas transport: {self:?}")
    }
}
impl std::error::Error for MaterialSteadyGasError {}
impl From<MaterialError> for MaterialSteadyGasError {
    fn from(e: MaterialError) -> Self {
        Self::Material(e)
    }
}
impl From<MaterialConvectionError> for MaterialSteadyGasError {
    fn from(e: MaterialConvectionError) -> Self {
        Self::Convection(e)
    }
}
impl From<ResistanceError> for MaterialSteadyGasError {
    fn from(e: ResistanceError) -> Self {
        Self::Resistance(e)
    }
}
impl From<SteadyGasError> for MaterialSteadyGasError {
    fn from(e: SteadyGasError) -> Self {
        Self::Thermal(e)
    }
}
fn check(cx: &Cx<'_>) -> Result<(), MaterialSteadyGasError> {
    cx.checkpoint()
        .map_err(|_| MaterialSteadyGasError::Cancelled)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> MaterialSteadyGasError {
    MaterialSteadyGasError::Invalid(reason)
}

/// Close gas enthalpy/species transport, material convection and hydraulic loss.
///
/// Each cell must start with `wall_conductance_w_k = 0`: this adapter is the
/// sole owner of that coefficient. Initial temperatures are numerical guesses;
/// the first native solve obtains the mixture before assigning its convection.
/// All other source, geometry, inlet, injection and phase data remain unchanged.
/// The original stationary-species solve determines each local composition.
///
/// Each subsequent solve uses the previous field's material conductances and
/// recomputes them at the new field. Publication requires the native steady
/// residual PLUS the change in actual wall power to fit `steady_power_tolerance_w`.
/// Work caps cover all solves; cancellation is checked between cells and by the
/// existing native solver. Refusal leaves `problem` and previous results intact.
///
/// The chosen correlations must admit a downstream mixed-cell bulk property
/// convention and the supplied geometry. This is a first-order steady reduced
/// model, not a resolved boundary layer. Pressure loss is evaluated at common
/// static pressure and bounded; it does not become an invented axial EOS field.
/// Hydraulic dissipation is reported separately and neglected in the thermal
/// equation. No pressure-work, reaction, compressible acceleration, finite
/// reservoir, phase-change or experimental-accuracy claim is added.
pub fn solve_material_gas_transport(
    cx: &Cx<'_>,
    problem: &SteadyAxialGas,
    transport: &[GasCellTransport],
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
) -> Result<MaterialSteadyGasSolution, MaterialSteadyGasError> {
    check(cx)?;
    let geometry = problem.to_network()?.cells;
    let mut trial = problem.clone();
    solve_feedback(
        cx,
        geometry,
        problem.static_pressure_pa,
        transport,
        controls,
        feedback,
        |cells, native| {
            for (cell, next) in trial.cells.iter_mut().zip(cells) {
                cell.initial_temperature_k = next.initial_temperature_k;
                cell.wall_conductance_w_k = next.wall_conductance_w_k;
            }
            trial.solve(native, || cx.checkpoint().is_err())
        },
    )
}

/// Close material convection on the SAME branched/recirculating gas transport.
///
/// Each cell is a well-mixed flow section whose declared area carries its total
/// outgoing mass flow. The supplied correlation and resistance must describe
/// that section; an unmodelled branch throat or mixing chamber is not inferred.
/// Mixtures, source-specific enthalpies and circulation come from the original
/// network owner. Its wall conductances must be zero before this sole material
/// convection owner is attached.
///
/// This uses the identical wall-power fixed point and cumulative work budget as
/// `solve_material_gas_transport`. Correlation and hydraulic evaluations retain
/// the complete accepted local material states. Pressure losses are reported per
/// cell; their sum conservatively gates the common-pressure approximation. No
/// pressure graph, pump, junction-loss model, or prescribed-flow feasibility is
/// invented. Hydraulic dissipation remains separate from the thermal balance.
pub fn solve_material_gas_network(
    cx: &Cx<'_>,
    problem: &SteadyGasNetwork,
    transport: &[GasCellTransport],
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
) -> Result<MaterialSteadyGasNetworkSolution, MaterialSteadyGasError> {
    check(cx)?;
    let mut trial = problem.clone();
    solve_feedback(
        cx,
        problem.cells.clone(),
        problem.static_pressure_pa,
        transport,
        controls,
        feedback,
        |cells, native| {
            trial.cells.clone_from_slice(cells);
            trial.solve(native, || cx.checkpoint().is_err())
        },
    )
}

fn solve_feedback<S: SolvedGas>(
    cx: &Cx<'_>,
    mut trial: Vec<GasNetworkCell>,
    static_pressure_pa: f64,
    transport: &[GasCellTransport],
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
    mut solve: impl FnMut(&[GasNetworkCell], SteadyGasSolve) -> Result<S, SteadyGasError>,
) -> Result<MaterialSteadyGasSolution<S>, MaterialSteadyGasError> {
    check(cx)?;
    if trial.is_empty()
        || trial.len() > MAX_STEADY_GAS_CELLS
        || transport.len() != trial.len()
        || !(1..=64).contains(&feedback.maximum_iterations)
        || feedback.maximum_substeps == 0
        || feedback.maximum_substeps > 1_000_000
        || !positive(feedback.conductance_relaxation)
        || feedback.conductance_relaxation > 1.0
        || !positive(feedback.maximum_pressure_drop_fraction)
        || feedback.maximum_pressure_drop_fraction > 0.1
        || !positive(controls.steady_power_tolerance_w)
        || !positive(static_pressure_pa)
    {
        return Err(invalid(
            "matched bounded cells/work, positive residual tolerance and feedback limits required",
        ));
    }
    let mut areas = Vec::with_capacity(trial.len());
    for (cell, law) in trial.iter().zip(transport) {
        check(cx)?;
        if cell.wall_conductance_w_k != 0.0 {
            return Err(invalid(
                "material convection cannot duplicate an authored wall conductance",
            ));
        }
        let area = cell.flow_area_m2;
        if !positive(area)
            || !positive(law.characteristic_length.value())
            || (area - law.resistance.flow_area().value()).abs()
                > 64.0 * f64::EPSILON * area.max(law.resistance.flow_area().value())
        {
            return Err(invalid(
                "thermal and hydraulic geometry must use the same positive flow area",
            ));
        }
        if matches!(
            law.correlation,
            CorrelationId::CircularDuctLaminarChf | CorrelationId::RectangularDuctLaminarChf
        ) {
            return Err(invalid(
                "prescribed wall temperatures cannot use a uniform-heat-flux correlation",
            ));
        }
        areas.push(area);
    }
    let mut substeps = 0;
    let mut residual_bound_w = f64::INFINITY;
    for iteration in 1..=feedback.maximum_iterations {
        check(cx)?;
        if substeps >= feedback.maximum_substeps {
            return Err(MaterialSteadyGasError::NotConverged {
                iterations: iteration - 1,
                substeps,
                residual_bound_w,
            });
        }
        let native = SteadyGasSolve {
            maximum_substeps: controls
                .maximum_substeps
                .min(feedback.maximum_substeps - substeps),
            steady_power_tolerance_w: controls.steady_power_tolerance_w * 0.5,
            ..controls
        };
        let solution = solve(&trial, native);
        check(cx)?;
        let solution = solution?;
        substeps += solution.substeps();
        let mut samples = Vec::with_capacity(trial.len());
        let mut wall_defect: f64 = 0.0;
        for (index, ((cell, law), &area)) in trial.iter().zip(transport).zip(&areas).enumerate() {
            check(cx)?;
            let gas = solution.cell_gas(index)?;
            let mass = solution.mass_flow_kg_s()[index];
            let velocity = (mass / gas.state().density) / area;
            if !positive(velocity) {
                return Err(invalid("local gas velocity is not representable"));
            }
            let difference = cell.wall_temperature_k - gas.state().temperature;
            let mut settings = law.settings;
            if difference > 0.0 {
                settings.direction = ThermalDirection::HeatingFluid;
            } else if difference < 0.0 {
                settings.direction = ThermalDirection::CoolingFluid;
            }
            let evaluated = evaluate_forced_from_material(
                law.correlation,
                &gas,
                velocity,
                law.characteristic_length,
                settings,
            )?;
            let (snapshot, nusselt, coefficient) = evaluated.into_parts();
            let conductance = coefficient.value.value() * cell.geometry.surface_area_m2;
            let power = conductance * difference;
            let defect = (conductance - cell.wall_conductance_w_k) * difference;
            if !positive(conductance) || !power.is_finite() || !defect.is_finite() {
                return Err(invalid(
                    "convective conductance or wall power is not representable",
                ));
            }
            wall_defect = wall_defect.max(defect.abs());
            samples.push(ThermalSample {
                gas,
                transport: snapshot,
                velocity_m_s: velocity,
                nusselt,
                coefficient,
                conductance_w_k: conductance,
                wall_heat_w: power,
            });
        }
        residual_bound_w = solution.maximum_power_residual_bound_w() + wall_defect;
        if !residual_bound_w.is_finite() {
            return Err(invalid("gas transport residual is not representable"));
        }
        check(cx)?;
        if residual_bound_w <= controls.steady_power_tolerance_w {
            // Loss laws do not drive this prescribed-mass thermal solve. Admit
            // them at the physical fixed point, never at an artificial thermal
            // guess that may have a different Reynolds number or pressure loss.
            let mut cells = Vec::with_capacity(samples.len());
            let mut pressure_drop = 0.0;
            let mut dissipation = 0.0;
            for (index, (sample, law)) in samples.into_iter().zip(transport).enumerate() {
                check(cx)?;
                let hydraulic = law.resistance.evaluate(
                    &sample.gas,
                    MassFlowRate::new(solution.mass_flow_kg_s()[index]),
                )?;
                pressure_drop += hydraulic.pressure_drop_pa;
                dissipation += hydraulic.dissipated_power_w;
                cells.push(GasCellTransportSample {
                    cell: index,
                    gas: sample.gas,
                    transport: sample.transport,
                    velocity_m_s: sample.velocity_m_s,
                    nusselt: sample.nusselt,
                    coefficient: sample.coefficient,
                    conductance_w_k: sample.conductance_w_k,
                    wall_heat_w: sample.wall_heat_w,
                    hydraulic,
                });
            }
            if !pressure_drop.is_finite() || !dissipation.is_finite() {
                return Err(invalid("total hydraulic response is not representable"));
            }
            let fraction = pressure_drop / static_pressure_pa;
            if !fraction.is_finite() || fraction > feedback.maximum_pressure_drop_fraction {
                return Err(MaterialSteadyGasError::PressureDrop {
                    fraction,
                    maximum: feedback.maximum_pressure_drop_fraction,
                });
            }
            check(cx)?;
            return Ok(MaterialSteadyGasSolution {
                thermal: solution,
                cells,
                iterations: iteration,
                substeps,
                maximum_power_residual_bound_w: residual_bound_w,
                total_pressure_drop_pa: pressure_drop,
                hydraulic_dissipation_w: dissipation,
            });
        }
        for (cell, sample) in trial.iter_mut().zip(samples) {
            check(cx)?;
            cell.initial_temperature_k = sample.gas.state().temperature;
            cell.wall_conductance_w_k = (1.0 - feedback.conductance_relaxation)
                * cell.wall_conductance_w_k
                + feedback.conductance_relaxation * sample.conductance_w_k;
        }
    }
    Err(MaterialSteadyGasError::NotConverged {
        iterations: feedback.maximum_iterations,
        substeps,
        residual_bound_w,
    })
}
