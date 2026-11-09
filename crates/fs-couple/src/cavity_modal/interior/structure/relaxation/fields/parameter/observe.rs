//! Independent original material, source and interface work parameter jets.
//! This recovery consumes one native solve and owns no additional factorization.

use super::{Error, Prepared, PreparedDirection, State, poll};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicPowerParameterJet, HarmonicPowerSample, HarmonicReceiverParameterJet,
    HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet,
};
use fs_modal::harmonic::{
    HarmonicResponseError, harmonic_mean_power_parameter_jet_with_cancel,
    harmonic_quadratic_power_parameter_jet_with_cancel,
};

/// One original Maxwell arm's local material/source parameter accounts.
/// Signed derivatives are per authored parameter unit and per rad/s.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationMemoryParameterJet {
    /// Parameter derivative of the original cycle-mean spring storage [J].
    pub mean_stored_energy_parameter_derivative: f64,
    /// Mixed frequency/parameter derivative of that same spring storage.
    pub mean_stored_energy_mixed_derivative: f64,
    /// Independent original |stress|²/(2*k*tau) heat, at all five orders.
    pub dissipation: HarmonicPowerParameterJet,
    /// Signed original structural motion into this branch, including b_theta.
    pub structural_work: HarmonicPowerParameterJet,
    /// Signed original pressure into this branch, including b_theta.
    pub acoustic_work: HarmonicPowerParameterJet,
}

/// Original physical work in watts and its full frequency/parameter jet.
/// Interface exchange and signed branch work remain separate from material heat.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationPowerParameterJet {
    pub structural_input: HarmonicPowerParameterJet,
    pub acoustic_input: HarmonicPowerParameterJet,
    /// Positive from fluid into solid through the original B/B^T interface.
    pub fluid_to_structure: HarmonicPowerParameterJet,
    pub structural_dissipation: HarmonicPowerParameterJet,
    pub wall_dissipation: HarmonicPowerParameterJet,
    pub relaxation_dissipation: HarmonicPowerParameterJet,
    pub structural_to_relaxation: HarmonicPowerParameterJet,
    pub acoustic_to_relaxation: HarmonicPowerParameterJet,
}
impl StructuralInteriorRelaxationPowerParameterJet {
    pub fn structural_balance(&self) -> Result<HarmonicPowerParameterJet, Error> {
        signed_sum(
            [
                self.structural_input,
                self.fluid_to_structure,
                self.structural_dissipation,
                self.structural_to_relaxation,
            ],
            [1., 1., -1., -1.],
        )
    }
    pub fn acoustic_balance(&self) -> Result<HarmonicPowerParameterJet, Error> {
        signed_sum(
            [
                self.acoustic_input,
                self.fluid_to_structure,
                self.wall_dissipation,
                self.acoustic_to_relaxation,
            ],
            [1., -1., -1., -1.],
        )
    }
    pub fn relaxation_balance(&self) -> Result<HarmonicPowerParameterJet, Error> {
        signed_sum(
            [
                self.structural_to_relaxation,
                self.acoustic_to_relaxation,
                self.relaxation_dissipation,
            ],
            [1., 1., -1.],
        )
    }
    pub fn total_balance(&self) -> Result<HarmonicPowerParameterJet, Error> {
        signed_sum(
            [
                self.structural_input,
                self.acoustic_input,
                self.structural_dissipation,
                self.wall_dissipation,
                self.relaxation_dissipation,
            ],
            [1., 1., -1., -1., -1.],
        )
    }
}

pub(super) struct Observed {
    pub(super) memory: Vec<StructuralInteriorRelaxationMemoryParameterJet>,
    pub(super) power: StructuralInteriorRelaxationPowerParameterJet,
    pub(super) parameter_force: Vec<C64>,
    pub(super) parameter_flow: Vec<C64>,
    pub(super) mixed_force: Vec<C64>,
    pub(super) mixed_flow: Vec<C64>,
}

/// Observe all original physical channels without inferring heat from a balance.
/// States and sources follow value, omega, omega-omega, theta, omega-theta order.
pub(super) fn observe(
    p: &Prepared<'_>,
    direction: &PreparedDirection<'_, '_>,
    native: &fs_phs::harmonic::HarmonicFrequencyParameterResponse,
    states: [&State; 5],
    sources: HarmonicTransferFieldParameterJet<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Observed, Error> {
    poll(cancelled)?;
    let (n, ns) = (p.source.coordinates, p.source.structural);
    let na = n - ns;
    let velocity = field(states.map(|state| state.velocity.as_slice()));
    let solid = slice(velocity, 0..ns);
    let pressure = slice(velocity, ns..n);
    let mut projection = Projection::new(n);
    let mut force: [Vec<C64>; 5] = core::array::from_fn(|_| Vec::with_capacity(ns));
    for s in 0..ns {
        let projected = projection.row(
            pressure,
            |a| {
                (
                    p.source.coupling[a * ns + s],
                    direction.coupling[a * ns + s],
                )
            },
            cancelled,
        )?;
        for (out, value) in force.iter_mut().zip(projected) {
            out.push(value);
        }
    }
    let mut flow: [Vec<C64>; 5] = core::array::from_fn(|_| Vec::with_capacity(na));
    for a in 0..na {
        let projected = projection.row(
            solid,
            |s| {
                (
                    p.source.coupling[a * ns + s],
                    direction.coupling[a * ns + s],
                )
            },
            cancelled,
        )?;
        for (out, value) in flow.iter_mut().zip(projected) {
            out.push(value);
        }
    }
    let structural_input = work(solid, slice(sources, 0..ns), cancelled)?;
    let acoustic_input = work(pressure, slice(sources, ns..n), cancelled)?;
    let fluid_to_structure = work(solid, field(force.each_ref().map(Vec::as_slice)), cancelled)?;
    let structural_dissipation = heat(p, direction, solid, 0, cancelled)?;
    let wall_dissipation = heat(p, direction, pressure, ns, cancelled)?;

    let curvature = &native.frequency_response;
    let first = &curvature.frequency_response;
    let effort = [
        &first.response.effort,
        &first.effort_derivative,
        &curvature.effort_second_derivative,
        &native.parameter_response.effort_derivative,
        &native.effort_mixed_derivative,
    ];
    let mut memory = Vec::with_capacity(p.branches.len());
    for (j, branch) in p.branches.iter().enumerate() {
        poll(cancelled)?;
        // Canonical effort is -stress/sqrt(k). Its full native parameter
        // map already contains k_theta; the spring coefficient is fixed1/2.
        let e = effort.map(|values| [values[2 * n + j]]);
        let e = field(e.each_ref().map(|v| v.as_slice()));
        let stored = quadratic(e, &[0.5], None, cancelled)?;
        let rate = branch.relaxation_time_s.recip();
        let rate_derivative = direction.native.branch_rate_derivatives()[j];
        let dissipation = quadratic(e, &[rate], Some(&[rate_derivative]), cancelled)?;
        let structural_flow = projection
            .row(
                solid,
                |i| (branch.strain[i], direction.branches[j].projection[i]),
                cancelled,
            )?
            .map(|value| [value]);
        let acoustic_flow = projection
            .row(
                pressure,
                |a| {
                    (
                        branch.strain[ns + a],
                        direction.branches[j].projection[ns + a],
                    )
                },
                cancelled,
            )?
            .map(|value| [value]);
        let stress = states.map(|state| [state.stress[j]]);
        let stress = field(stress.each_ref().map(|v| v.as_slice()));
        let structural_work = work(
            field(structural_flow.each_ref().map(|v| v.as_slice())),
            stress,
            cancelled,
        )?;
        let acoustic_work = work(
            field(acoustic_flow.each_ref().map(|v| v.as_slice())),
            stress,
            cancelled,
        )?;
        memory.push(StructuralInteriorRelaxationMemoryParameterJet {
            mean_stored_energy_parameter_derivative: stored.parameter_derivative,
            mean_stored_energy_mixed_derivative: stored.angular_parameter_derivative,
            dissipation,
            structural_work,
            acoustic_work,
        });
    }
    let power = StructuralInteriorRelaxationPowerParameterJet {
        structural_input,
        acoustic_input,
        fluid_to_structure,
        structural_dissipation,
        wall_dissipation,
        relaxation_dissipation: memory_sum(&memory, |m| m.dissipation, cancelled)?,
        structural_to_relaxation: memory_sum(&memory, |m| m.structural_work, cancelled)?,
        acoustic_to_relaxation: memory_sum(&memory, |m| m.acoustic_work, cancelled)?,
    };
    let [_, _, _, parameter_force, mixed_force] = force;
    let [_, _, _, parameter_flow, mixed_flow] = flow;
    poll(cancelled)?;
    Ok(Observed {
        memory,
        power,
        parameter_force,
        parameter_flow,
        mixed_force,
        mixed_flow,
    })
}

fn heat(
    p: &Prepared<'_>,
    direction: &PreparedDirection<'_, '_>,
    velocity: HarmonicTransferFieldParameterJet<'_>,
    offset: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerParameterJet, Error> {
    let (n, size) = (p.source.coordinates, velocity.values.len());
    let mut matrix = Vec::with_capacity(size * size);
    let mut tangent = Vec::with_capacity(size * size);
    for i in 0..size {
        poll(cancelled)?;
        let range = (offset + i) * n + offset..(offset + i) * n + offset + size;
        matrix.extend_from_slice(&p.source.damping[range.clone()]);
        tangent.extend_from_slice(&direction.blocks[1][range]);
    }
    // The shared quadratic owner symmetrizes original pairs before work,
    // preserving signed directions, subnormal coefficients and gyro cancellation.
    quadratic(velocity, &matrix, Some(&tangent), cancelled)
}

/// Reuse two original rows for all reciprocal/interface and branch projections.
/// Each changing-map product is normalized jointly with its changing field.
struct Projection {
    weights: Vec<C64>,
    direction: Vec<C64>,
}
impl Projection {
    fn new(maximum: usize) -> Self {
        Self {
            weights: Vec::with_capacity(maximum),
            direction: Vec::with_capacity(maximum),
        }
    }
    fn row(
        &mut self,
        values: HarmonicTransferFieldParameterJet<'_>,
        mut coefficient: impl FnMut(usize) -> (f64, f64),
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[C64; 5], Error> {
        poll(cancelled)?;
        let n = values.values.len();
        self.weights.resize(n, C64::ZERO);
        self.direction.resize(n, C64::ZERO);
        for j in 0..n {
            let (value, derivative) = coefficient(j);
            self.weights[j] = C64::from_re(value);
            self.direction[j] = C64::from_re(derivative);
        }
        let projected = HarmonicTransferParameterJet::project_complex_with_derivatives(
            HarmonicReceiverParameterJet {
                weights: &self.weights,
                angular_derivative: None,
                angular_second_derivative: None,
                parameter_derivative: Some(&self.direction),
                angular_parameter_derivative: None,
            },
            values,
        )
        .map_err(observation)?;
        poll(cancelled)?;
        Ok([
            projected.sample.amplitude,
            projected.sample.angular_derivative,
            projected.angular_second_derivative,
            projected.parameter_derivative,
            projected.angular_parameter_derivative,
        ])
    }
}

fn field(orders: [&[C64]; 5]) -> HarmonicTransferFieldParameterJet<'_> {
    HarmonicTransferFieldParameterJet {
        values: orders[0],
        angular_derivative: orders[1],
        angular_second_derivative: orders[2],
        parameter_derivative: orders[3],
        angular_parameter_derivative: orders[4],
    }
}
fn slice(
    input: HarmonicTransferFieldParameterJet<'_>,
    range: core::ops::Range<usize>,
) -> HarmonicTransferFieldParameterJet<'_> {
    field(
        [
            input.values,
            input.angular_derivative,
            input.angular_second_derivative,
            input.parameter_derivative,
            input.angular_parameter_derivative,
        ]
        .map(|v| &v[range.clone()]),
    )
}
fn quadratic(
    values: HarmonicTransferFieldParameterJet<'_>,
    matrix: &[f64],
    direction: Option<&[f64]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerParameterJet, Error> {
    harmonic_quadratic_power_parameter_jet_with_cancel(values, matrix, direction, cancelled)
        .map_err(observation)
}
fn work(
    flow: HarmonicTransferFieldParameterJet<'_>,
    effort: HarmonicTransferFieldParameterJet<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerParameterJet, Error> {
    harmonic_mean_power_parameter_jet_with_cancel(flow, effort, cancelled).map_err(observation)
}
fn observation(error: HarmonicResponseError) -> Error {
    match error {
        HarmonicResponseError::Cancelled => Error::Cancelled,
        _ => Error::Invalid("original physical parameter observation is unrepresentable"),
    }
}
fn from_projection(value: HarmonicTransferParameterJet) -> HarmonicPowerParameterJet {
    HarmonicPowerParameterJet {
        sample: HarmonicPowerSample {
            average_power_w: value.sample.amplitude.re,
            angular_derivative: value.sample.angular_derivative.re,
        },
        angular_second_derivative: value.angular_second_derivative.re,
        parameter_derivative: value.parameter_derivative.re,
        angular_parameter_derivative: value.angular_parameter_derivative.re,
    }
}
fn signed_sum<const N: usize>(
    powers: [HarmonicPowerParameterJet; N],
    signs: [f64; N],
) -> Result<HarmonicPowerParameterJet, Error> {
    let values = powers.map(|p| C64::from_re(p.sample.average_power_w));
    let first = powers.map(|p| C64::from_re(p.sample.angular_derivative));
    let second = powers.map(|p| C64::from_re(p.angular_second_derivative));
    let parameter = powers.map(|p| C64::from_re(p.parameter_derivative));
    let mixed = powers.map(|p| C64::from_re(p.angular_parameter_derivative));
    HarmonicTransferParameterJet::project_real(
        &signs,
        field([&values, &first, &second, &parameter, &mixed]),
    )
    .map(from_projection)
    .map_err(observation)
}
fn memory_sum(
    memory: &[StructuralInteriorRelaxationMemoryParameterJet],
    select: fn(&StructuralInteriorRelaxationMemoryParameterJet) -> HarmonicPowerParameterJet,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerParameterJet, Error> {
    if memory.is_empty() {
        return signed_sum(
            [HarmonicPowerParameterJet {
                sample: HarmonicPowerSample {
                    average_power_w: 0.,
                    angular_derivative: 0.,
                },
                angular_second_derivative: 0.,
                parameter_derivative: 0.,
                angular_parameter_derivative: 0.,
            }],
            [1.],
        );
    }
    let mut columns: [Vec<C64>; 5] = core::array::from_fn(|_| Vec::with_capacity(memory.len()));
    for account in memory {
        poll(cancelled)?;
        let value = select(account);
        for (column, value) in columns.iter_mut().zip([
            value.sample.average_power_w,
            value.sample.angular_derivative,
            value.angular_second_derivative,
            value.parameter_derivative,
            value.angular_parameter_derivative,
        ]) {
            column.push(C64::from_re(value));
        }
    }
    let signs = vec![1.; memory.len()];
    HarmonicTransferParameterJet::project_real(&signs, field(columns.each_ref().map(Vec::as_slice)))
        .map(from_projection)
        .map_err(observation)
}
