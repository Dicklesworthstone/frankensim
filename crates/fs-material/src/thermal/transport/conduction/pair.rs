//! Reciprocal finite-body heat exchange on original native caloric storage.
//!
//! The same bounded source-interval driver serves the original phase owner.
//! This adapter retains h for Cp, u for Cv, and the declared storage basis for
//! volumetric/amount-bound laws. It constructs neither a phase chart nor entropy.

use super::{CaloricConductionLayer, OriginalConduction};
use super::super::capacity::{CapacityCurve, CapacityKnot};
use crate::phase::exchange::pair::{
    self as phase_pair, FinitePhaseBody, FinitePhaseExchangeError,
    kernel::{NumericPath, Pair, Segment, Solution, Work},
    nonnegative, positive, unresolved,
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{DepositionError, ThermalStorageBasis};

/// Same physical trajectory tolerances and whole-call work ceilings as the
/// shared original-phase driver. Native exchange makes no extra entropy query.
pub use phase_pair::FinitePhaseExchangeControls as FiniteCaloricExchangeControls;

/// Fixed mass and the accepted energy on the SAME original material law.
///
/// The caller retains mechanical/volume constraints and spatial lumping.
/// The amount-bound storage variant additionally admits this exact reference
/// mass here; a cell consumer must also call the original validate_cell.
#[derive(Clone, Copy, Debug)]
pub struct FiniteCaloricBody<'l, 'a, 'm> {
    pub layer: &'l CaloricConductionLayer<'a, 'm>,
    pub initial_specific_energy_j_kg: f64,
    pub mass_kg: f64,
}

/// One atomic proposal, in the supplied first/second order.
///
/// One Q debits and credits the two original native accounts. A physical
/// interpretation as total internal energy still follows each original law's
/// mechanical convention; this transport adds no pressure-work or entropy law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteCaloricExchange {
    pub caloric_identities: [Option<[u8; 32]>; 2],
    pub storage_bases: [ThermalStorageBasis; 2],
    pub initial_specific_energies_j_kg: [f64; 2],
    pub final_specific_energies_j_kg: [f64; 2],
    pub initial_temperatures_k: [f64; 2],
    /// Original inverse of each returned native energy, after quantization.
    pub temperatures_k: [f64; 2],
    pub conductance_time_j_k: f64,
    /// Positive from the first body to the second; the only exchanged heat.
    pub heat_first_to_second_j: f64,
    pub specific_heat_j_kg: [f64; 2],
    pub represented_energy_change_j: [f64; 2],
    /// Represented stored energy minus each body's signed supplied heat [J].
    /// These native-state residuals are separate from flow integration error.
    pub energy_residual_j: [f64; 2],
    pub exposure_residual_j_k: f64,
    pub exposure_error_allowance_j_k: f64,
    pub heat_error_bound_j: f64,
    pub chart_intervals: usize,
    pub quadrature_panels: usize,
    pub root_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteCaloricExchangeError {
    Material(DepositionError),
    /// Original shared numerical/source-support or cancellation refusal.
    Integration(FinitePhaseExchangeError),
}
impl From<DepositionError> for FiniteCaloricExchangeError {
    fn from(error: DepositionError) -> Self { Self::Material(error) }
}
impl From<FinitePhaseExchangeError> for FiniteCaloricExchangeError {
    fn from(error: FinitePhaseExchangeError) -> Self { Self::Integration(error) }
}
impl core::fmt::Display for FiniteCaloricExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite native caloric exchange: {self:?}")
    }
}
impl std::error::Error for FiniteCaloricExchangeError {}
type Result<T> = core::result::Result<T, FiniteCaloricExchangeError>;

/// Integrate dQ/dX=T_first(e_first-Q/m_first)-T_second(e_second+Q/m_second)
/// with held X=G*dt [J/K]. Original capacity intervals and native source limits
/// stay authoritative. Constant capacities use analytic attenuation; nonlinear
/// capacities share the source-split quadrature and bounded scalar inverse.
///
/// Preparation copies only original interval endpoints/capacities, once per
/// call, and never queries a material card. No timeline or accepted state changes.
/// Relative Q remains available even if a native energy increment rounds away;
/// both explicit residuals require the caller's physical energy gate.
pub fn exchange(
    first: FiniteCaloricBody<'_, '_, '_>,
    second: FiniteCaloricBody<'_, '_, '_>,
    conductance_time_j_k: f64,
    controls: FiniteCaloricExchangeControls,
) -> Result<FiniteCaloricExchange> {
    exchange_with_checkpoint(first, second, conductance_time_j_k, controls, || false)
}

/// Cancellation covers both original source preparations, the shared numeric
/// solve and original inverse admission of both proposed native endpoints.
pub fn exchange_with_checkpoint(
    first: FiniteCaloricBody<'_, '_, '_>,
    second: FiniteCaloricBody<'_, '_, '_>,
    conductance_time_j_k: f64,
    controls: FiniteCaloricExchangeControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FiniteCaloricExchange> {
    let mut work = Work { cancelled: &mut cancelled, controls, charts: 0, panels: 0 };
    work.poll()?;
    controls.validate()?;
    nonnegative(conductance_time_j_k, "conductance-time exposure")?;
    let bodies = [first, second];
    let mut temperatures = [0.0; 2];
    for (side, body) in bodies.iter().enumerate() {
        work.poll()?;
        positive(body.mass_kg, "fixed body mass")?;
        temperatures[side] = body.layer.state_at_specific_energy(
            body.initial_specific_energy_j_kg,
        )?.temperature_k;
        if let OriginalConduction::HeatStorage(law) = body.layer.original {
            let reference = law.caloric().reference_mass_kg();
            if (body.mass_kg - reference).abs()
                > 64.0 * f64::EPSILON * body.mass_kg.max(reference)
            {
                return Err(DepositionError::InvalidInput(
                    "finite body mass differs from its original heat-storage assembly",
                ).into());
            }
        }
    }
    if conductance_time_j_k == 0.0 || temperatures[0] == temperatures[1] {
        return finish(bodies, temperatures, None, conductance_time_j_k,
            0.0, Solution::default(), &mut work);
    }
    let first_hot = temperatures[0] > temperatures[1];
    let directions = if first_hot { [-1.0, 1.0] } else { [1.0, -1.0] };
    let first_path = prepare(first, temperatures[0], directions[0], &mut work)?;
    let second_path = prepare(second, temperatures[1], directions[1], &mut work)?;
    let pair = if first_hot {
        Pair::new(&first_path, &second_path)
    } else { Pair::new(&second_path, &first_path) };
    let solution = pair.solve(conductance_time_j_k, &mut work)?;
    let signed_heat = if first_hot { solution.heat } else { -solution.heat };
    finish(bodies, temperatures, Some([&first_path, &second_path]),
        conductance_time_j_k, signed_heat, solution, &mut work)
}

fn prepare(
    body: FiniteCaloricBody<'_, '_, '_>,
    initial_temperature: f64,
    direction: f64,
    work: &mut Work<'_>,
) -> Result<NumericPath> {
    let bounds = body.layer.energy_bounds_j_kg();
    let end_energy = bounds[usize::from(direction > 0.0)];
    if direction * (end_energy - body.initial_specific_energy_j_kg) <= 0.0 {
        return Err(FinitePhaseExchangeError::OutsideSourceSupport.into());
    }
    let end_temperature = body.layer.state_at_specific_energy(end_energy)?.temperature_k;
    match body.layer.original {
        OriginalConduction::EquilibriumPhase(law) => {
            let curve = law.caloric_curve();
            let initial = curve.state_at_specific_enthalpy(body.initial_specific_energy_j_kg)
                .map_err(FinitePhaseExchangeError::from)?;
            Ok(phase_pair::prepare_path(
                FinitePhaseBody { curve, initial, mass_kg: body.mass_kg },
                direction, Some(end_energy), work,
            )?)
        }
        OriginalConduction::Isobaric(law) => prepare_capacity(
            body, &law.caloric().source_line().knots, 1.0,
            initial_temperature, end_temperature, end_energy, direction, work,
        ),
        OriginalConduction::ExpandingIsobaric(law) => prepare_capacity(
            body, &law.caloric().caloric().source_line().knots, 1.0,
            initial_temperature, end_temperature, end_energy, direction, work,
        ),
        OriginalConduction::Isochoric(law) => prepare_capacity(
            body, &law.caloric().source_line().knots, 1.0,
            initial_temperature, end_temperature, end_energy, direction, work,
        ),
        OriginalConduction::FixedVolume(law) => prepare_capacity(
            body, &law.path().source_lines()[0].knots, law.reference_density_kg_m3(),
            initial_temperature, end_temperature, end_energy, direction, work,
        ),
        OriginalConduction::HeatStorage(law) => prepare_capacity(
            body, law.caloric().storage().capacity_knots(), law.caloric().reference_mass_kg(),
            initial_temperature, end_temperature, end_energy, direction, work,
        ),
    }
}

// Keep C in its original units. The sole factor is m for Cp/Cv, V=m/rho_ref
// for volumetric capacity, or m/m_ref for a correctly admitted amount assembly.
// Source Cp is never divided and multiplied back through a surrogate table.
#[allow(clippy::too_many_arguments)]
fn prepare_capacity<K: CapacityKnot>(
    body: FiniteCaloricBody<'_, '_, '_>,
    knots: &[K],
    source_energy_per_specific: f64,
    initial_temperature: f64,
    end_temperature: f64,
    end_energy: f64,
    direction: f64,
    work: &mut Work<'_>,
) -> Result<NumericPath> {
    let factor = positive(body.mass_kg / source_energy_per_specific, "native source energy factor")?;
    let native_limit = positive(
        body.mass_kg * (end_energy - body.initial_specific_energy_j_kg).abs(),
        "native source energy limit",
    )?;
    let upper = if direction > 0.0 {
        knots.partition_point(|k| k.temperature() <= initial_temperature)
    } else {
        knots.partition_point(|k| k.temperature() < initial_temperature)
    };
    if upper == 0 || upper == knots.len()
        || direction * (end_temperature - initial_temperature) <= 0.0
    {
        return Err(FinitePhaseExchangeError::OutsideSourceSupport.into());
    }
    let mut index = upper - 1;
    let mut temperature = initial_temperature;
    let mut total = Dd::ZERO;
    let mut segments = Vec::new();
    let mut admitted_limit = native_limit;
    loop {
        work.poll()?;
        if work.charts == work.controls.maximum_chart_intervals {
            return Err(FinitePhaseExchangeError::WorkLimit.into());
        }
        work.charts += 1;
        let lower = &knots[index];
        let upper = &knots[index + 1];
        let original = [
            [lower.temperature(), lower.capacity()],
            [upper.temperature(), upper.capacity()],
        ];
        let curve = CapacityCurve::new(&original);
        let (start_t, end_t, capacities) = if direction > 0.0 {
            (lower.temperature(), upper.temperature(), [lower.capacity(), upper.capacity()])
        } else {
            (upper.temperature(), lower.temperature(), [upper.capacity(), lower.capacity()])
        };
        positive(factor * capacities[0], "lower native total capacity")?;
        positive(factor * capacities[1], "upper native total capacity")?;
        let source_energy_width = curve.energy_change(start_t, end_t)
            .map_err(|_| DepositionError::MaterialDomain)?.abs();
        let start_q = if segments.is_empty() {
            -factor * curve.energy_change(start_t, initial_temperature)
                .map_err(|_| DepositionError::MaterialDomain)?.abs()
        } else { total.to_f64() };
        let reaches_temperature_bound = direction * (end_t - end_temperature) >= 0.0;
        if reaches_temperature_bound {
            let supported = factor * curve.energy_change(temperature, end_temperature)
                .map_err(|_| DepositionError::MaterialDomain)?.abs();
            admitted_limit = admitted_limit.min((total + Dd::from_f64(supported)).to_f64());
        }
        let increment = positive(
            factor * curve.energy_change(temperature, end_t)
                .map_err(|_| DepositionError::MaterialDomain)?.abs(),
            "original native source heat interval",
        )?;
        total = total + Dd::from_f64(increment);
        let end_q = total.to_f64();
        if !start_q.is_finite() || !end_q.is_finite() || end_q <= start_q {
            return Err(unresolved("original native source breakpoint").into());
        }
        segments.push(Segment {
            start_q, end_q, start_t, end_t, end_energy: None,
            capacities, energy_factor: factor, source_energy_width,
        });
        // The independent Q limit honors BOTH native energy and prepared
        // temperature support without altering original inverse endpoints.
        if reaches_temperature_bound || end_q >= admitted_limit { break; }
        temperature = end_t;
        if direction > 0.0 {
            if index + 2 == knots.len() { break; }
            index += 1;
        } else {
            if index == 0 { break; }
            index -= 1;
        }
    }
    Ok(NumericPath {
        initial_temperature, segments,
        limit: positive(total.to_f64().min(admitted_limit), "admitted native heat interval")?,
    })
}

#[allow(clippy::too_many_arguments)]
fn finish(
    bodies: [FiniteCaloricBody<'_, '_, '_>; 2],
    initial_temperatures: [f64; 2],
    paths: Option<[&NumericPath; 2]>,
    exposure: f64,
    signed_heat: f64,
    solution: Solution,
    work: &mut Work<'_>,
) -> Result<FiniteCaloricExchange> {
    let mut final_energies = [0.0; 2];
    let mut temperatures = initial_temperatures;
    let mut specific = [0.0; 2];
    let mut stored = [0.0; 2];
    let mut residuals = [0.0; 2];
    for (side, heat) in [-signed_heat, signed_heat].into_iter().enumerate() {
        work.poll()?;
        let body = bodies[side];
        specific[side] = heat / body.mass_kg;
        if !specific[side].is_finite() || (heat != 0.0 && specific[side] == 0.0) {
            return Err(DepositionError::Arithmetic.into());
        }
        final_energies[side] = if heat == 0.0 {
            body.initial_specific_energy_j_kg
        } else if paths.is_some_and(|p| heat.abs() == p[side].limit) {
            body.layer.energy_bounds_j_kg()[usize::from(heat > 0.0)]
        } else { body.initial_specific_energy_j_kg + specific[side] };
        temperatures[side] = body.layer.state_at_specific_energy(final_energies[side])?.temperature_k;
        stored[side] = body.mass_kg
            * (final_energies[side] - body.initial_specific_energy_j_kg);
        residuals[side] = stored[side] - heat;
        if !stored[side].is_finite() || !residuals[side].is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
    }
    if (signed_heat > 0.0 && temperatures[0] < temperatures[1])
        || (signed_heat < 0.0 && temperatures[1] < temperatures[0])
    {
        return Err(unresolved("represented native pair state crosses thermal equilibrium").into());
    }
    work.poll()?;
    Ok(FiniteCaloricExchange {
        caloric_identities: bodies.map(|b| b.layer.caloric_law().caloric_identity()),
        storage_bases: bodies.map(|b| b.layer.caloric_law().storage_basis()),
        initial_specific_energies_j_kg: bodies.map(|b| b.initial_specific_energy_j_kg),
        final_specific_energies_j_kg: final_energies,
        initial_temperatures_k: initial_temperatures, temperatures_k: temperatures,
        conductance_time_j_k: exposure, heat_first_to_second_j: signed_heat,
        specific_heat_j_kg: specific, represented_energy_change_j: stored,
        energy_residual_j: residuals,
        exposure_residual_j_k: solution.residual,
        exposure_error_allowance_j_k: solution.allowance,
        heat_error_bound_j: solution.heat_error,
        chart_intervals: work.charts, quadrature_panels: work.panels,
        root_iterations: solution.iterations,
    })
}
