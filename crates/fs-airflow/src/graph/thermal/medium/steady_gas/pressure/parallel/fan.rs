//! Density-aware native fans on the original heated parallel gas branches.
//!
//! The retained fan's nominal STATIC pressure matches the native manifold
//! pressure. Source density affects both fan affinity and every gas branch.
//! The inlet thermodynamic state is frozen while omitted fan work is bounded;
//! this is not a fan heating model, total-pressure curve or uncertainty enclosure.

use fs_exec::Cx;
use fs_material::gas_field::steady::SteadyGasSolve;
use fs_qty::{Pressure, VolumetricFlowRate};

use super::{
    GasTransportFeedback, ParallelGasBranch, ParallelGasBudget, ParallelGasError,
    ParallelGasSolution, checkpoint, positive, solve_parallel_pressure_driven_material_gas,
};
use crate::graph::thermal::medium::{
    FrozenFlowMedium,
    fan::{FanGasDomain, MaterialFanError},
};
use crate::{AirflowError, FanBank};

/// Original equipment and explicit physical applicability, with no efficiency default.
#[derive(Clone, Copy, Debug)]
pub struct GasFanSource<'a> {
    /// Original unbound pressure/volume curve, arrangement and speed.
    pub reference_bank: &'a FanBank,
    /// Original material at which this curve was rated.
    pub reference_medium: FrozenFlowMedium<'a>,
    /// Native density/Reynolds similarity and common-pressure applicability.
    pub gas_domain: FanGasDomain,
    /// Independently authored lower bound on BANK static efficiency
    /// eta_static = (static pressure rise * actual inlet volume flow)/shaft power.
    /// In (0,1], valid at the declared speed, medium and EVERY evaluated flow.
    /// A static pressure curve alone does not establish this bound.
    pub minimum_static_efficiency: f64,
}

/// Bounded nominal pressure matching and explicit omitted fan-work admission.
#[derive(Clone, Copy, Debug)]
pub struct GasFanControls {
    /// Positive increasing pressure bracket [Pa]; no automatic expansion.
    /// Every trial must be inside every branch and fan source domain.
    pub pressure_bracket_pa: (f64, f64),
    pub absolute_pressure_tolerance_pa: f64,
    pub relative_pressure_tolerance: f64,
    /// Complete network solves including endpoints, in 2..=128.
    pub maximum_network_evaluations: usize,
    /// Maximum declared upper bound on omitted bank shaft power [W].
    /// It is separate from branch/manifold mechanical omission allowances.
    pub maximum_neglected_shaft_power_w: f64,
}

/// The exact accepted thermal field and the native fan evaluation at its inlet.
#[derive(Debug)]
pub struct GasFanOperatingPoint<'a> {
    pub source: GasFanSource<'a>,
    /// Density-bound native bank; source, stall range, tolerance and speed survive.
    pub fan: FanBank,
    pub network: ParallelGasSolution,
    pub fan_pressure: Pressure,
    /// Native manifold pressure minus nominal fan static rise [Pa].
    pub pressure_residual: Pressure,
    /// Nominal static pressure-volume supply [W], not heat or shaft power.
    pub static_pressure_power_w: f64,
    /// Q * nominal static rise * (1 + original pressure tolerance) / eta_min.
    /// Conditional on the declared efficiency bound; not a measured heat input.
    pub declared_shaft_power_upper_w: f64,
    pub network_evaluations: usize,
    /// Actual native work over ALL network/pressure/material trials.
    pub total_flow_evaluations: usize,
    pub total_substeps: usize,
}

/// No failed, unsupported or unconverged trial publishes a candidate.
#[derive(Debug)]
pub enum GasFanError {
    Invalid(&'static str),
    Network(ParallelGasError),
    Binding(MaterialFanError),
    Fan(AirflowError),
    WorkBudget {
        network_evaluations: usize,
        flow_evaluations: usize,
        substeps: usize,
    },
    NeglectedShaftPower {
        upper_w: f64,
        allowance_w: f64,
    },
    NotBracketed {
        lower_residual_pa: f64,
        upper_residual_pa: f64,
    },
    NotConverged {
        evaluations: usize,
        last_residual_pa: f64,
    },
}
impl core::fmt::Display for GasFanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source gas fan operating point: {self:?}")
    }
}
impl std::error::Error for GasFanError {}
impl From<ParallelGasError> for GasFanError {
    fn from(error: ParallelGasError) -> Self {
        Self::Network(error)
    }
}
impl From<MaterialFanError> for GasFanError {
    fn from(error: MaterialFanError) -> Self {
        Self::Binding(error)
    }
}
impl From<AirflowError> for GasFanError {
    fn from(error: AirflowError) -> Self {
        Self::Fan(error)
    }
}

/// Match an original source-rated fan to the actual heated gas flow split.
///
/// Each positive pressure trial invokes the existing parallel solve with only
/// the remaining shared native work. Its EXACT source-resolved inlet binds the
/// original fan through FrozenFlowMedium::bind_fan_bank; total actual inlet
/// volume flow, not mixed-outlet or reference-density flow, selects the fan point.
/// Native source, stall, density/Reynolds and pressure limits remain active.
///
/// The common thermodynamic inlet stays frozen across omitted fan heating and
/// compression. The authored efficiency bound and retained pressure tolerance
/// bound omitted shaft power, which is reported separately and never added as
/// heat or subtracted again from branch dissipation. Motor electrical losses,
/// transient fan dynamics and a complete energy closure are not inferred.
///
/// Heating/buoyancy can make the system nonmonotone. The explicit sign bracket
/// selects a candidate, without uniqueness/stability or propagated uncertainty
/// certification. Any unsupported trial refuses; the accepted native network is
/// returned directly, without a second thermal or flow solve.
pub fn solve_parallel_material_gas_with_fan<'a>(
    cx: &Cx<'_>,
    branches: &[ParallelGasBranch<'_>],
    source: GasFanSource<'a>,
    thermal: SteadyGasSolve,
    feedback: GasTransportFeedback,
    work: ParallelGasBudget,
    controls: GasFanControls,
) -> Result<GasFanOperatingPoint<'a>, GasFanError> {
    checkpoint(cx)?;
    work.admit()?;
    let (mut low, mut high) = controls.pressure_bracket_pa;
    if !positive(low)
        || !positive(high)
        || low >= high
        || !positive(controls.absolute_pressure_tolerance_pa)
        || !(0.0..1.0).contains(&controls.relative_pressure_tolerance)
        || !(2..=128).contains(&controls.maximum_network_evaluations)
        || !positive(source.minimum_static_efficiency)
        || source.minimum_static_efficiency > 1.0
        || !controls.maximum_neglected_shaft_power_w.is_finite()
        || controls.maximum_neglected_shaft_power_w < 0.0
    {
        return Err(GasFanError::Invalid(
            "positive fan pressure bracket, bounded work, efficiency and shaft-power allowance required",
        ));
    }
    let (mut total_flow, mut total_substeps) = (0, 0);
    let mut probe = |pressure: f64, calls: usize| -> Result<_, GasFanError> {
        checkpoint(cx)?;
        let remaining_flow = work.maximum_flow_evaluations.saturating_sub(total_flow);
        let remaining_substeps = work.maximum_substeps.saturating_sub(total_substeps);
        if remaining_flow < 2 || remaining_substeps == 0 {
            return Err(GasFanError::WorkBudget {
                network_evaluations: calls - 1,
                flow_evaluations: total_flow,
                substeps: total_substeps,
            });
        }
        let network = solve_parallel_pressure_driven_material_gas(
            cx,
            branches,
            Pressure::new(pressure),
            thermal,
            feedback,
            ParallelGasBudget {
                maximum_flow_evaluations: remaining_flow,
                maximum_substeps: remaining_substeps,
                ..work
            },
        )?;
        total_flow += network.flow_evaluations;
        total_substeps += network.substeps;
        let bound = FrozenFlowMedium::ThermochemicalGas(&network.inlet).bind_fan_bank(
            source.reference_bank,
            source.reference_medium,
            source.gas_domain,
        )?;
        let volume = network.inlet_volume_flow_m3_per_s;
        let fan_pressure = bound.pressure_at(VolumetricFlowRate::new(volume))?;
        // Positive static efficiency cannot bound a running fan at zero static
        // rise. Such an endpoint needs separate shaft-power data and refuses.
        let static_power = fan_pressure.value() * volume;
        let shaft_upper = static_power * (1.0 + bound.curve().pressure_tolerance_rel())
            / source.minimum_static_efficiency;
        let residual = pressure - fan_pressure.value();
        if !positive(fan_pressure.value())
            || !positive(static_power)
            || !positive(shaft_upper)
            || !residual.is_finite()
            || network.driving_pressure_pa != pressure
        {
            return Err(GasFanError::Invalid(
                "fan static supply, shaft bound or native pressure is inconsistent/unrepresentable",
            ));
        }
        if shaft_upper > controls.maximum_neglected_shaft_power_w {
            return Err(GasFanError::NeglectedShaftPower {
                upper_w: shaft_upper,
                allowance_w: controls.maximum_neglected_shaft_power_w,
            });
        }
        let fan = bound.bank().clone();
        drop(bound);
        checkpoint(cx)?;
        Ok(GasFanOperatingPoint {
            source,
            fan,
            network,
            fan_pressure,
            pressure_residual: Pressure::new(residual),
            static_pressure_power_w: static_power,
            declared_shaft_power_upper_w: shaft_upper,
            network_evaluations: calls,
            total_flow_evaluations: total_flow,
            total_substeps,
        })
    };
    let closed = |point: &GasFanOperatingPoint<'_>| {
        point.pressure_residual.value().abs()
            <= controls.absolute_pressure_tolerance_pa.max(
                controls.relative_pressure_tolerance
                    * point
                        .network
                        .driving_pressure_pa
                        .max(point.fan_pressure.value()),
            )
    };
    let lower = probe(low, 1)?;
    if closed(&lower) {
        return Ok(lower);
    }
    let mut low_residual = lower.pressure_residual.value();
    drop(lower);
    let upper = probe(high, 2)?;
    if closed(&upper) {
        return Ok(upper);
    }
    let mut last_residual = upper.pressure_residual.value();
    if low_residual.is_sign_positive() == last_residual.is_sign_positive() {
        return Err(GasFanError::NotBracketed {
            lower_residual_pa: low_residual,
            upper_residual_pa: last_residual,
        });
    }
    drop(upper);
    let mut evaluations = 2;
    for calls in 3..=controls.maximum_network_evaluations {
        let middle = 0.5 * low + 0.5 * high;
        if middle == low || middle == high {
            break;
        }
        let point = probe(middle, calls)?;
        evaluations = calls;
        if closed(&point) {
            return Ok(point);
        }
        last_residual = point.pressure_residual.value();
        if last_residual.is_sign_positive() == low_residual.is_sign_positive() {
            low = middle;
            low_residual = last_residual;
        } else {
            high = middle;
        }
    }
    Err(GasFanError::NotConverged {
        evaluations,
        last_residual_pa: last_residual,
    })
}
