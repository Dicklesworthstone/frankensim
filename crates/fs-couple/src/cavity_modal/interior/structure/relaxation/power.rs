//! Actual original heat and work extrema of the complete causal cavity.

use super::{
    Error, PreparedStructuralInteriorRelaxation, PreparedStructuralInteriorRelaxationDirection,
    StructuralInteriorRelaxationCurvatureResponse,
    StructuralInteriorRelaxationFrequencyParameterResponse,
    StructuralInteriorRelaxationFrequencyResponse,
    StructuralInteriorRelaxationSourceOrder as Order, StructuralInteriorRelaxationSourceSpectrum,
    fields,
};
use core::cell::RefCell;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    self as query, HalfPowerBandwidth, HalfPowerBandwidthParameterResponse,
    HarmonicPowerCurvatureSample, HarmonicPowerParameterJet, HarmonicPowerSample,
    HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet, PowerExtremum,
    PowerExtremumParameterResponse, PowerParameterSurvey, PowerSurvey, ResonanceControls,
    ResonanceError,
};

/// One original signed work account, or explicitly identified constitutive heat.
///
/// All channels use watts. Source/interface/branch work can be negative; no
/// absolute-value transform is applied. Total dissipation sums structural,
/// memoryless-wall and original Maxwell heat exactly once. A branch is indexed
/// in the caller's original order, including equal-time and dark arms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StructuralInteriorRelaxationPowerChannel {
    StructuralInput,
    AcousticInput,
    FluidToStructure,
    StructuralDissipation,
    WallDissipation,
    RelaxationDissipation,
    StructuralToRelaxation,
    AcousticToRelaxation,
    TotalDissipation,
    BranchDissipation(usize),
    BranchStructuralWork(usize),
    BranchAcousticWork(usize),
}
type Channel = StructuralInteriorRelaxationPowerChannel;

impl PreparedStructuralInteriorRelaxation<'_, '_, '_, '_> {
    /// Resolve one isolated maximum or minimum of actual signed physical work.
    /// The bracket must resolve one simple slope reversal. Only two-state
    /// frequency samples are solved while refining; one final three-state
    /// sample supplies the original curvature. All material states survive.
    pub fn power_extremum<S: StructuralInteriorRelaxationSourceSpectrum + ?Sized>(
        &self,
        channel: Channel,
        sources: &mut S,
        bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<PowerExtremum, ResonanceError<Error>> {
        let source = RefCell::new(sources);
        query::refine_power_extremum(
            bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|p| p.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Curvature,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(curvature)
            },
        )
    }

    /// Follow the same simple work extremum under an authored original
    /// M/C/K/B, material k/tau/strain and source direction. The direction is
    /// prepared once; only the reserved final sample performs five state solves.
    /// This is local differentiation, not continuation across feature switches.
    pub fn power_extremum_with_parameter_derivative<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
    >(
        &self,
        channel: Channel,
        sources: &mut S,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<PowerExtremumParameterResponse, ResonanceError<Error>> {
        let source = RefCell::new(sources);
        query::refine_power_extremum_with_parameter_derivative(
            bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|p| p.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Parameter,
                    Some(&direction.inner),
                    &mut || cx.checkpoint().is_err(),
                )
            },
        )
    }

    /// Resolve both actual P_peak/2 crossings around the SAME positive maximum.
    /// Neither power squared nor half amplitude defines this linewidth. The
    /// reported Q is measured peak frequency/width, not a pole decay rate.
    /// The supplied peak, source spectrum and channel must come from one query.
    pub fn power_bandwidth<S: StructuralInteriorRelaxationSourceSpectrum + ?Sized>(
        &self,
        channel: Channel,
        sources: &mut S,
        peak: &PowerExtremum,
        lower_bracket_rad_s: (f64, f64),
        upper_bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<HalfPowerBandwidth, ResonanceError<Error>> {
        let source = RefCell::new(sources);
        query::refine_power_bandwidth(
            peak,
            lower_bracket_rad_s,
            upper_bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|p| p.sample)
            },
        )
    }

    /// Differentiate both physical half-power crossings, their separation and
    /// observed Q. Each edge uses the same source/material jet as the supplied
    /// peak. Only the two reserved final edge calls request all five orders.
    pub fn power_bandwidth_with_parameter_derivative<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
    >(
        &self,
        channel: Channel,
        sources: &mut S,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        peak: &PowerExtremumParameterResponse,
        lower_bracket_rad_s: (f64, f64),
        upper_bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<Error>> {
        let source = RefCell::new(sources);
        query::refine_power_bandwidth_with_parameter_derivative(
            peak,
            lower_bracket_rad_s,
            upper_bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|p| p.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Parameter,
                    Some(&direction.inner),
                    &mut || cx.checkpoint().is_err(),
                )
            },
        )
    }

    /// Resolve the signed maxima and minima visible on a positive rad/s grid.
    /// One bounded shared budget covers the cached grid and all refinements;
    /// one final curvature call is reserved per candidate. This is not a
    /// complete free-mode census, and unresolved peaks remain caller-owned.
    pub fn survey_power_extrema<S: StructuralInteriorRelaxationSourceSpectrum + ?Sized>(
        &self,
        channel: Channel,
        sources: &mut S,
        frequencies_rad_s: &[f64],
        maximum_extrema: usize,
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<PowerSurvey, ResonanceError<Error>> {
        let source = RefCell::new(sources);
        query::survey_power_extrema(
            frequencies_rad_s,
            controls,
            maximum_extrema,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|p| p.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Curvature,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(curvature)
            },
        )
    }

    /// Local material/source sensitivity of every grid-resolved work extremum.
    /// The original physical direction is reused across frequencies; all
    /// complete state, field and work jets are evaluated only at final extrema.
    pub fn survey_power_extrema_with_parameter_derivative<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
    >(
        &self,
        channel: Channel,
        sources: &mut S,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        frequencies_rad_s: &[f64],
        maximum_extrema: usize,
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<PowerParameterSurvey, ResonanceError<Error>> {
        let source = RefCell::new(sources);
        query::survey_power_extrema_with_parameter_derivative(
            frequencies_rad_s,
            controls,
            maximum_extrema,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|p| p.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    channel,
                    &source,
                    w,
                    Order::Parameter,
                    Some(&direction.inner),
                    &mut || cx.checkpoint().is_err(),
                )
            },
        )
    }
}

/// Shared physical evaluation used by the public query and native regressions.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate<S: StructuralInteriorRelaxationSourceSpectrum + ?Sized>(
    prepared: &fields::Prepared<'_>,
    channel: Channel,
    sources: &RefCell<&mut S>,
    omega: f64,
    order: Order,
    direction: Option<&fields::PreparedDirection<'_, '_>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerParameterJet, Error> {
    if cancelled() {
        return Err(Error::Cancelled);
    }
    let mut source = sources.borrow_mut();
    let jet = source.jet(omega, order)?;
    if cancelled() {
        return Err(Error::Cancelled);
    }
    let factor = prepared.factor(omega, cancelled)?;
    let result = match order {
        Order::Frequency => {
            let response = factor.solve_with_frequency_derivative(
                jet.values.structural_force,
                jet.values.inward_flow_load,
                jet.angular_derivative
                    .map(|s| (s.structural_force, s.inward_flow_load)),
                cancelled,
            )?;
            lift(channel.frequency(&response)?)
        }
        Order::Curvature => {
            let response = factor.solve_with_frequency_curvature(
                jet.values.structural_force,
                jet.values.inward_flow_load,
                jet.angular_derivative
                    .map(|s| (s.structural_force, s.inward_flow_load)),
                jet.angular_second_derivative
                    .map(|s| (s.structural_force, s.inward_flow_load)),
                cancelled,
            )?;
            let sample = channel.curvature(&response)?;
            HarmonicPowerParameterJet {
                angular_second_derivative: sample.angular_second_derivative,
                ..lift(sample.sample)
            }
        }
        Order::Parameter => {
            let direction = direction.ok_or(Error::Invalid(
                "original prepared physical direction required",
            ))?;
            let response =
                factor.solve_with_frequency_parameter_derivative(jet, direction, cancelled)?;
            channel.parameter(&response)?
        }
    };
    if cancelled() {
        return Err(Error::Cancelled);
    }
    Ok(result)
}

impl Channel {
    fn frequency(
        self,
        r: &StructuralInteriorRelaxationFrequencyResponse,
    ) -> Result<HarmonicPowerSample, Error> {
        let p = r.power;
        Ok(match self {
            Self::StructuralInput => p.structural_input,
            Self::AcousticInput => p.acoustic_input,
            Self::FluidToStructure => p.fluid_to_structure,
            Self::StructuralDissipation => p.structural_dissipation,
            Self::WallDissipation => p.wall_dissipation,
            Self::RelaxationDissipation => p.relaxation_dissipation,
            Self::StructuralToRelaxation => p.structural_to_relaxation,
            Self::AcousticToRelaxation => p.acoustic_to_relaxation,
            Self::TotalDissipation => fields::power::sum(&[
                p.structural_dissipation,
                p.wall_dissipation,
                p.relaxation_dissipation,
            ])?,
            Self::BranchDissipation(i) => r.memory.get(i).ok_or(branch_error())?.dissipation,
            Self::BranchStructuralWork(i) => r.memory.get(i).ok_or(branch_error())?.structural_work,
            Self::BranchAcousticWork(i) => r.memory.get(i).ok_or(branch_error())?.acoustic_work,
        })
    }
    fn curvature(
        self,
        r: &StructuralInteriorRelaxationCurvatureResponse,
    ) -> Result<HarmonicPowerCurvatureSample, Error> {
        let p = r.power;
        Ok(match self {
            Self::StructuralInput => p.structural_input,
            Self::AcousticInput => p.acoustic_input,
            Self::FluidToStructure => p.fluid_to_structure,
            Self::StructuralDissipation => p.structural_dissipation,
            Self::WallDissipation => p.wall_dissipation,
            Self::RelaxationDissipation => p.relaxation_dissipation,
            Self::StructuralToRelaxation => p.structural_to_relaxation,
            Self::AcousticToRelaxation => p.acoustic_to_relaxation,
            Self::TotalDissipation => curvature(sum([
                p.structural_dissipation,
                p.wall_dissipation,
                p.relaxation_dissipation,
            ]
            .map(|v| HarmonicPowerParameterJet {
                angular_second_derivative: v.angular_second_derivative,
                ..lift(v.sample)
            }))?),
            Self::BranchDissipation(i) => r.memory.get(i).ok_or(branch_error())?.dissipation,
            Self::BranchStructuralWork(i) => r.memory.get(i).ok_or(branch_error())?.structural_work,
            Self::BranchAcousticWork(i) => r.memory.get(i).ok_or(branch_error())?.acoustic_work,
        })
    }
    fn parameter(
        self,
        r: &StructuralInteriorRelaxationFrequencyParameterResponse,
    ) -> Result<HarmonicPowerParameterJet, Error> {
        let p = r.power;
        Ok(match self {
            Self::StructuralInput => p.structural_input,
            Self::AcousticInput => p.acoustic_input,
            Self::FluidToStructure => p.fluid_to_structure,
            Self::StructuralDissipation => p.structural_dissipation,
            Self::WallDissipation => p.wall_dissipation,
            Self::RelaxationDissipation => p.relaxation_dissipation,
            Self::StructuralToRelaxation => p.structural_to_relaxation,
            Self::AcousticToRelaxation => p.acoustic_to_relaxation,
            Self::TotalDissipation => sum([
                p.structural_dissipation,
                p.wall_dissipation,
                p.relaxation_dissipation,
            ])?,
            Self::BranchDissipation(i) => r.memory.get(i).ok_or(branch_error())?.dissipation,
            Self::BranchStructuralWork(i) => r.memory.get(i).ok_or(branch_error())?.structural_work,
            Self::BranchAcousticWork(i) => r.memory.get(i).ok_or(branch_error())?.acoustic_work,
        })
    }
}
fn branch_error() -> Error {
    Error::Invalid("original memory branch index is out of range")
}
fn lift(sample: HarmonicPowerSample) -> HarmonicPowerParameterJet {
    HarmonicPowerParameterJet {
        sample,
        angular_second_derivative: 0.0,
        parameter_derivative: 0.0,
        angular_parameter_derivative: 0.0,
    }
}
fn curvature(p: HarmonicPowerParameterJet) -> HarmonicPowerCurvatureSample {
    HarmonicPowerCurvatureSample {
        sample: p.sample,
        angular_second_derivative: p.angular_second_derivative,
    }
}
fn sum(values: [HarmonicPowerParameterJet; 3]) -> Result<HarmonicPowerParameterJet, Error> {
    let real = |f: fn(HarmonicPowerParameterJet) -> f64| values.map(|v| C64::from_re(f(v)));
    let total = HarmonicTransferParameterJet::project_real(
        &[1.0; 3],
        HarmonicTransferFieldParameterJet {
            values: &real(|p| p.sample.average_power_w),
            angular_derivative: &real(|p| p.sample.angular_derivative),
            angular_second_derivative: &real(|p| p.angular_second_derivative),
            parameter_derivative: &real(|p| p.parameter_derivative),
            angular_parameter_derivative: &real(|p| p.angular_parameter_derivative),
        },
    )
    .map_err(|_| Error::Invalid("total original heat jet is unrepresentable"))?;
    Ok(HarmonicPowerParameterJet {
        sample: HarmonicPowerSample {
            average_power_w: total.sample.amplitude.re,
            angular_derivative: total.sample.angular_derivative.re,
        },
        angular_second_derivative: total.angular_second_derivative.re,
        parameter_derivative: total.parameter_derivative.re,
        angular_parameter_derivative: total.angular_parameter_derivative.re,
    })
}
