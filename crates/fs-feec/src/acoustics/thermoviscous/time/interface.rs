//! Causal normal-traction connection of two COMPLETE resolved fluid owners.
//!
//! The original P1 traction traces determine C with psi_left C=psi_right.
//! p_left=C lambda, p_right=lambda and C^T Q_left+Q_right=0 preserve work.
//! The shared pHS common-flow constraint evolves all original fields together;
//! no frequency response, fitted pole, modal truncation or partitioned lag enters.
//!
//! `connect` retains its Galerkin normal-flow, insulated interface. The separate
//! `connect_conforming` binding supplies the COMPLETE original P2 vector and
//! thermal trace to this SAME trajectory owner. Equilibrium geometry, pressure
//! and temperature are fixed; composition is not advected or mixed.
use super::*;
use fs_phs::linear_exact::CommonFlowProjection;
use std::collections::BTreeMap;

mod conforming;
pub use conforming::{
    ConformingThermoviscousInterface, ConformingThermoviscousTimePairStep,
    PreparedConformingThermoviscousTimePair,
};

const QUADRATIC_FORMS: usize = 2 * HEAT_CHANNELS + 4;

// Original load columns and their physical kind, not another dynamical owner.
struct InterfaceColumns {
    count: usize,
    mechanical: usize,
    loads: Vec<f64>,
}

struct PairEvaluation {
    regions: [ThermoviscousTimeStep; 2],
    mechanical_work_j: [f64; 2],
    conductive_work_j: [f64; 2],
    external_pressure_work_j: f64,
    mean_reaction: Vec<f64>,
    endpoint_reaction: Vec<f64>,
}

/// Both admitted original field endpoints and same-trajectory physical work.
#[derive(Clone, Debug)]
pub struct ThermoviscousTimePairStep {
    /// Left then right. Each region's pressure work includes its SIGNED
    /// interface work; the two internal transfers cancel in the enclosing pair.
    pub regions: [ThermoviscousTimeStep; 2],
    /// Independently integrated lambda(t)*Q(t), never an energy remainder.
    pub interface_work_j: [f64; 2],
    /// Actual external pressure sources only, excluding both internal transfers.
    pub external_pressure_work_j: f64,
    /// Multiplier coefficients in the right region's selected trace basis.
    pub mean_interface_traction_pa: Vec<f64>,
    pub endpoint_interface_traction_pa: Vec<f64>,
}
impl ThermoviscousTimePairStep {
    pub fn balance_residual_j(&self) -> f64 {
        self.regions
            .iter()
            .map(|r| r.work.delta_energy_j + r.work.dissipated_j() - r.work.structural_work_j)
            .sum::<f64>()
            - self.external_pressure_work_j
    }
}

/// One fixed-clock trajectory, borrowing both immutable source owners. All
/// external pressure coordinates and original structural forces remain inputs.
/// Warm calls allocate candidate vectors; this is a bounded dense local owner.
pub struct PreparedThermoviscousTimePair<'p, 'a, 'c> {
    regions: [&'p PreparedThermoviscousTime<'a, 'c>; 2],
    selected: [Vec<usize>; 2],
    traction_map: Vec<f64>,
    external_pressure: Vec<(usize, usize)>,
    interface_count: usize,
    mechanical_count: usize,
    system: PortHamiltonian,
    projection: CommonFlowProjection,
    flow: PreparedQuadraticFlow,
    maximum_relative_energy_residual: f64,
}

impl<'a, 'c> PreparedThermoviscousTime<'a, 'c> {
    /// Connect original physical traces to another original resolved time owner.
    /// Both preparations must use the same clock and reference pressure/temperature.
    /// Unselected pressures stay external in left-then-right original mode order.
    pub fn connect<'p>(
        &'p self,
        right: &'p Self,
        left_modes: &[usize],
        right_modes: &[usize],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousTimePair<'p, 'a, 'c>, ThermoviscousTimeError> {
        PreparedThermoviscousTimePair::new(self, right, left_modes, right_modes, limits, cancelled)
    }
}

impl<'p, 'a, 'c> PreparedThermoviscousTimePair<'p, 'a, 'c> {
    fn new(
        left: &'p PreparedThermoviscousTime<'a, 'c>,
        right: &'p PreparedThermoviscousTime<'a, 'c>,
        left_modes: &[usize],
        right_modes: &[usize],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousTimeError> {
        poll(cancelled)?;
        let regions = [left, right];
        let selected = [left_modes.to_vec(), right_modes.to_vec()];
        let count = left_modes.len();
        for side in 0..2 {
            if count == 0
                || selected[side].len() != count
                || selected[side]
                    .iter()
                    .any(|&p| p >= regions[side].model.loads.len())
                || selected[side]
                    .iter()
                    .copied()
                    .collect::<BTreeSet<_>>()
                    .len()
                    != count
            {
                return Err(invalid(
                    "unique equal-rank original interface pressure modes required",
                ));
            }
        }
        let traction_map = trace_map(left.model, right.model, left_modes, right_modes, cancelled)?;
        Self::prepare(regions, selected, traction_map, None, limits, cancelled)
    }

    fn prepare(
        regions: [&'p PreparedThermoviscousTime<'a, 'c>; 2],
        selected: [Vec<usize>; 2],
        traction_map: Vec<f64>,
        columns: Option<InterfaceColumns>,
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousTimeError> {
        poll(cancelled)?;
        let [left, right] = regions;
        let left_modes = &selected[0];
        let right_modes = &selected[1];
        let count = columns.as_ref().map_or(left_modes.len(), |c| c.count);
        let mechanical_count = columns.as_ref().map_or(count, |c| c.mechanical);
        let same = |a: f64, b: f64| (a - b).abs() <= 64.0 * f64::EPSILON * a.abs().max(b.abs());
        if left.duration_s.to_bits() != right.duration_s.to_bits()
            || !same(
                left.model.medium.pressure_pa,
                right.model.medium.pressure_pa,
            )
            || !same(
                left.model.medium.temperature_k,
                right.model.medium.temperature_k,
            )
        {
            return Err(invalid(
                "coupled original fields require a common clock and equilibrium",
            ));
        }
        let sizes = [left.system.state_dim(), right.system.state_dim()];
        let starts = [0, sizes[0]];
        let n = sizes.iter().sum::<usize>();
        let original_ports = [left.system.port_dim(), right.system.port_dim()];
        let port_starts = [0, original_ports[0]];
        let all_ports = original_ports.iter().sum::<usize>();
        let external_pressure: Vec<_> = (0..2)
            .flat_map(|side| {
                (0..regions[side].model.loads.len()).filter_map({
                    let selected = &selected;
                    move |p| (!selected[side].contains(&p)).then_some((side, p))
                })
            })
            .collect();
        let external: Vec<_> = external_pressure
            .iter()
            .copied()
            .chain((0..2).flat_map(|side| {
                (regions[side].model.loads.len()..original_ports[side]).map(move |p| (side, p))
            }))
            .collect();
        let ports = external.len();
        if n > limits.max_states
            || 2 * (n + ports + 1) > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || !limits.maximum_relative_energy_residual.is_finite()
            || limits.maximum_relative_energy_residual <= 0.0
            || limits.maximum_relative_energy_residual > 1e-6
        {
            return Err(invalid(
                "complete coupled state/port capacity and bounded tolerance required",
            ));
        }
        let mut q = vec![0.0; n * n];
        let mut j = vec![0.0; n * n];
        let mut r = vec![0.0; n * n];
        let mut g = vec![0.0; n * ports];
        let custom = columns.is_some();
        let mut b = columns.map_or_else(|| vec![0.0; n * count], |c| c.loads);
        if count == 0 || count > n || mechanical_count > count || b.len() != n * count {
            return Err(invalid("complete independent physical interface columns required"));
        }
        let mut forms = vec![vec![0.0; n * n]; QUADRATIC_FORMS];
        let mut observations = vec![vec![0.0; n]; all_ports];
        for side in 0..2 {
            let source = regions[side];
            let size = sizes[side];
            let start = starts[side];
            let (sj, sr, sg) = source.system.structure();
            for col in 0..size {
                poll(cancelled)?;
                let mut unit = vec![0.0; size];
                unit[col] = 1.0;
                let effort = source.system.effort(&unit);
                for row in 0..size {
                    let at = (start + row) * n + start + col;
                    q[at] = effort[row];
                    j[at] = sj[row * size + col];
                    r[at] = sr[row * size + col];
                    for h in 0..HEAT_CHANNELS {
                        forms[side * HEAT_CHANNELS + h][at] = source.heat[h][row * size + col];
                    }
                }
                for p in 0..original_ports[side] {
                    observations[port_starts[side] + p][start + col] = (0..size)
                        .map(|row| sg[row * original_ports[side] + p] * effort[row])
                        .sum();
                }
            }
            for row in 0..size {
                for (p, &(owner, source_port)) in external.iter().enumerate() {
                    if owner == side {
                        g[(start + row) * ports + p] = sg[row * original_ports[side] + source_port];
                    }
                }
                for a in 0..if custom { 0 } else { count } {
                    b[(start + row) * count + a] = if side == 0 {
                        (0..count)
                            .map(|i| {
                                sg[row * original_ports[side] + left_modes[i]]
                                    * traction_map[i * count + a]
                            })
                            .sum()
                    } else {
                        sg[row * original_ports[side] + right_modes[a]]
                    };
                }
            }
        }
        let original = PortHamiltonian::new(
            n,
            ports,
            j,
            r,
            g,
            Box::new(QuadraticStorage::new(q.clone(), n)?),
        )?;
        let projection =
            CommonFlowProjection::new(&original, &b, count, limits.max_states, cancelled)?;
        let (j, r, g) = projection.structure();
        let system = PortHamiltonian::new(
            n,
            ports,
            j.to_vec(),
            r.to_vec(),
            g.to_vec(),
            Box::new(QuadraticStorage::new(q.clone(), n)?),
        )?;
        let (reaction_state, _) = projection.reaction_maps();
        let mut interface_flows: [Vec<Vec<f64>>; 2] =
            std::array::from_fn(|_| vec![vec![0.0; n]; count]);
        for side in 0..2 {
            for a in 0..count {
                for col in 0..n {
                    interface_flows[side][a][col] = (starts[side]..starts[side] + sizes[side])
                        .map(|row| b[row * count + a] * q[row * n + col])
                        .sum();
                }
                let kind = usize::from(a >= mechanical_count);
                for row in 0..n {
                    poll(cancelled)?;
                    for col in 0..n {
                        forms[2 * HEAT_CHANNELS + 2 * side + kind][row * n + col] +=
                            reaction_state[a * n + row] * interface_flows[side][a][col];
                    }
                }
            }
        }
        observations.extend(interface_flows.into_iter().flatten());
        observations.extend((0..count).map(|a| reaction_state[a * n..(a + 1) * n].to_vec()));
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: limits.max_states,
                max_omega_dt: 1.0,
            },
        )?;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            left.duration_s,
            &forms,
            &observations,
            limits.preparation,
            &mut *cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            regions,
            selected,
            traction_map,
            external_pressure,
            interface_count: count,
            mechanical_count,
            system,
            projection,
            flow,
            maximum_relative_energy_residual: limits
                .maximum_relative_energy_residual
                .min(left.maximum_relative_energy_residual)
                .min(right.maximum_relative_energy_residual),
        })
    }

    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn duration_s(&self) -> f64 {
        self.regions[0].duration_s
    }
    /// P_left=C P_right, row-major selected-left by selected-right modes.
    pub fn traction_map(&self) -> &[f64] {
        &self.traction_map
    }
    pub fn interface_modes(&self) -> [&[usize]; 2] {
        [&self.selected[0], &self.selected[1]]
    }
    /// (region 0/1, original pressure mode), preserving original ordering.
    pub fn external_pressure_modes(&self) -> &[(usize, usize)] {
        &self.external_pressure
    }
    pub fn zero_states(&self) -> [ThermoviscousTimeState; 2] {
        [self.regions[0].zero_state(), self.regions[1].zero_state()]
    }

    /// Evolve both original fields under held external pressures and original
    /// structural forces. No state or clock is published on any refusal.
    pub fn step(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousTimePairStep, ThermoviscousTimeError> {
        let frame = self.evaluate(initial, external_traction_pa, wall_forces_n_per_sqrt_kg, cancelled)?;
        Ok(ThermoviscousTimePairStep {
            regions: frame.regions,
            interface_work_j: frame.mechanical_work_j,
            external_pressure_work_j: frame.external_pressure_work_j,
            mean_interface_traction_pa: frame.mean_reaction,
            endpoint_interface_traction_pa: frame.endpoint_reaction,
        })
    }

    fn evaluate(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PairEvaluation, ThermoviscousTimeError> {
        poll(cancelled)?;
        if external_traction_pa.len() != self.external_pressure.len()
            || (0..2)
                .any(|s| wall_forces_n_per_sqrt_kg[s].len() != self.regions[s].wall_scale.len())
        {
            return Err(invalid(
                "complete external pressure and original structural drives required",
            ));
        }
        let input: Vec<_> = external_traction_pa
            .iter()
            .chain(wall_forces_n_per_sqrt_kg[0])
            .chain(wall_forces_n_per_sqrt_kg[1])
            .copied()
            .collect();
        if input.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite coupled held inputs required"));
        }
        let mut x = self.regions[0].encode(initial[0])?;
        x.extend(self.regions[1].encode(initial[1])?);
        self.projection.admit_state(&x)?;
        let candidate = self.flow.evaluate(&x, &input)?;
        poll(cancelled)?;
        self.projection.admit_state(&candidate.state)?;
        let endpoint_reaction = self.projection.reaction(&candidate.state, &input)?;
        let count = self.interface_count;
        let all_ports = self
            .regions
            .iter()
            .map(|s| s.system.port_dim())
            .sum::<usize>();
        let outputs = &candidate.integrals[QUADRATIC_FORMS..QUADRATIC_FORMS + all_ports];
        let (_, reaction_input) = self.projection.reaction_maps();
        let feed: Vec<_> = (0..count)
            .map(|a| {
                input
                    .iter()
                    .enumerate()
                    .map(|(p, u)| reaction_input[a * input.len() + p] * u)
                    .sum::<f64>()
            })
            .collect();
        let mean_reaction: Vec<_> = (0..count)
            .map(|a| {
                candidate.integrals[QUADRATIC_FORMS + all_ports + 2 * count + a]
                    / self.duration_s() + feed[a]
            })
            .collect();
        let mut interface_work_j = [0.0; 2];
        let mut conductive_work_j = [0.0; 2];
        let mut external_work = [0.0; 2];
        let port_starts = [0, self.regions[0].system.port_dim()];
        for (p, &(side, original)) in self.external_pressure.iter().enumerate() {
            external_work[side] += external_traction_pa[p] * outputs[port_starts[side] + original];
        }
        for side in 0..2 {
            interface_work_j[side] = candidate.integrals[2 * HEAT_CHANNELS + 2 * side];
            conductive_work_j[side] = candidate.integrals[2 * HEAT_CHANNELS + 2 * side + 1];
            for a in 0..count {
                let integrated_flow = candidate.integrals[QUADRATIC_FORMS + all_ports + side * count + a];
                if a < self.mechanical_count {
                    interface_work_j[side] += feed[a] * integrated_flow;
                } else {
                    conductive_work_j[side] += feed[a] * integrated_flow;
                }
            }
        }
        let mut frames = Vec::with_capacity(2);
        let mut start = 0;
        for side in 0..2 {
            let source = self.regions[side];
            let end = start + source.system.state_dim();
            let heat = &candidate.integrals[side * HEAT_CHANNELS..(side + 1) * HEAT_CHANNELS];
            let flows = &outputs[port_starts[side]..port_starts[side] + source.system.port_dim()];
            let work = ThermoviscousTimeWork {
                initial_energy_j: source.system.hamiltonian(&x[start..end]),
                final_energy_j: source.system.hamiltonian(&candidate.state[start..end]),
                delta_energy_j: energy_change_j(
                    &source.system,
                    &x[start..end],
                    &candidate.state[start..end],
                ),
                pressure_work_j: external_work[side] + interface_work_j[side],
                conductive_boundary_work_j: conductive_work_j[side],
                structural_work_j: wall_forces_n_per_sqrt_kg[side]
                    .iter()
                    .enumerate()
                    .map(|(i, u)| u * flows[source.model.loads.len() + i])
                    .sum(),
                fluid_viscous_heat_j: heat[0],
                fluid_conductive_heat_j: heat[1],
                structural_damping_heat_j: heat[2],
                molecular_relaxation_heat_j: heat[3],
                solid_conductive_heat_j: heat[4],
            };
            admit_work(work, heat, self.maximum_relative_energy_residual)?;
            let state = source.decode(&candidate.state[start..end]);
            source.encode(&state)?;
            frames.push(ThermoviscousTimeStep {
                state,
                mean_inward_flow_m3_s: flows[..source.model.loads.len()]
                    .iter()
                    .map(|v| v / self.duration_s())
                    .collect(),
                mean_wall_velocity_m_sqrt_kg_per_s: flows[source.model.loads.len()..]
                    .iter()
                    .map(|v| v / self.duration_s())
                    .collect(),
                work,
            });
            start = end;
        }
        let result = PairEvaluation {
            regions: frames
                .try_into()
                .map_err(|_| invalid("complete pair candidate required"))?,
            mechanical_work_j: interface_work_j,
            conductive_work_j,
            external_pressure_work_j: external_work.iter().sum(),
            mean_reaction,
            endpoint_reaction,
        };
        let scale = result
            .regions
            .iter()
            .map(|r| {
                r.work.initial_energy_j.abs()
                    + r.work.final_energy_j.abs()
                    + r.work.supplied_j().abs()
                    + r.work.dissipated_j().abs()
            })
            .sum::<f64>()
            .max(f64::MIN_POSITIVE);
        if candidate.integrals.iter().any(|v| !v.is_finite())
            || !scale.is_finite()
            || result
                .mean_reaction
                .iter()
                .any(|v| !v.is_finite())
            || result.mechanical_work_j.iter().sum::<f64>().abs()
                > self.maximum_relative_energy_residual * scale
            || result.conductive_work_j.iter().sum::<f64>().abs()
                > self.maximum_relative_energy_residual * scale
            || (result.regions.iter().map(|r| r.work.delta_energy_j + r.work.dissipated_j()
                    - r.work.structural_work_j).sum::<f64>() - result.external_pressure_work_j).abs()
                > self.maximum_relative_energy_residual * scale
        {
            return Err(invalid(
                "original coupled interface work or enclosing energy balance refused",
            ));
        }
        poll(cancelled)?;
        Ok(result)
    }
}

type FaceKey = [[u64; 3]; 3];
type FaceModes<'a> = BTreeMap<FaceKey, Vec<(usize, &'a PhysicalTrace)>>;
fn faces<'a>(model: &'a ThermoviscousFe<'_>, selected: &[usize]) -> FaceModes<'a> {
    let mut result = BTreeMap::new();
    for (local, &global) in selected.iter().enumerate() {
        for trace in &model.physical_modes[global] {
            result
                .entry(trace.vertices)
                .or_insert_with(Vec::new)
                .push((local, trace));
        }
    }
    result
}
fn trace_map(
    left: &ThermoviscousFe<'_>,
    right: &ThermoviscousFe<'_>,
    left_modes: &[usize],
    right_modes: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, ThermoviscousTimeError> {
    let count = left_modes.len();
    let faces = [faces(left, left_modes), faces(right, right_modes)];
    if faces[0].keys().ne(faces[1].keys()) {
        return Err(invalid("original interface triangulations must coincide"));
    }
    let mut gram = [vec![0.0; count * count], vec![0.0; count * count]];
    let mut cross = vec![0.0; count * count];
    let mut values = Vec::new();
    for (face, left_traces) in &faces[0] {
        poll(cancelled)?;
        let right_traces = &faces[1][face];
        let a = left_traces[0].1;
        let b = right_traces[0].1;
        let dot: f64 = a
            .outward_normal
            .iter()
            .zip(b.outward_normal)
            .map(|(a, b)| a * b)
            .sum();
        if (dot + 1.0).abs() > 1e-12 || (a.area_m2 - b.area_m2).abs() > 1e-12 * a.area_m2 {
            return Err(invalid(
                "opposite original interface normals and equal face areas required",
            ));
        }
        let mut shape = [vec![[0.0; 3]; count], vec![[0.0; 3]; count]];
        for (side, traces) in [left_traces, right_traces].iter().enumerate() {
            for &(mode, trace) in *traces {
                shape[side][mode] = trace.values;
            }
        }
        let moment = |u: [f64; 3], v: [f64; 3]| {
            a.area_m2 / 12.0
                * (u.iter().sum::<f64>() * v.iter().sum::<f64>()
                    + u.iter().zip(v).map(|(u, v)| u * v).sum::<f64>())
        };
        for i in 0..count {
            for j in 0..count {
                for side in 0..2 {
                    gram[side][i * count + j] += moment(shape[side][i], shape[side][j]);
                }
                cross[i * count + j] += moment(shape[0][i], shape[1][j]);
            }
        }
        values.push(shape);
    }
    let mut roots = [Vec::new(), Vec::new()];
    let mut factors = Vec::new();
    for side in 0..2 {
        roots[side] = (0..count)
            .map(|i| det::sqrt(gram[side][i * count + i]))
            .collect();
        if roots[side].iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(invalid("nonzero original interface trace rank required"));
        }
        let scaled: Vec<_> = (0..count * count)
            .map(|i| gram[side][i] / roots[side][i / count] / roots[side][i % count])
            .collect();
        let factor = cholesky(&scaled, count)?;
        if (0..count).any(|i| factor.l(i, i).powi(2) <= 512.0 * f64::EPSILON * count as f64) {
            return Err(invalid("original interface trace rank is unresolved at working precision"));
        }
        factors.push(factor);
    }
    let mut map = vec![0.0; count * count];
    for j in 0..count {
        let mut rhs: Vec<_> = (0..count)
            .map(|i| cross[i * count + j] / roots[0][i])
            .collect();
        factors[0].solve(&mut rhs);
        for i in 0..count {
            map[i * count + j] = rhs[i] / roots[0][i];
        }
    }
    // Full-rank overlap alone does not prove matching spaces. Verify the
    // actual P1 values at every original triangle vertex, not just face area.
    for shape in values {
        poll(cancelled)?;
        for j in 0..count {
            for v in 0..3 {
                let value: f64 = (0..count)
                    .map(|i| shape[0][i][v] * map[i * count + j])
                    .sum();
                let bound = shape[1][j][v].abs()
                    + (0..count)
                        .map(|i| (shape[0][i][v] * map[i * count + j]).abs())
                        .sum::<f64>();
                if !value.is_finite()
                    || !bound.is_finite()
                    || (value - shape[1][j][v]).abs() > 512.0 * f64::EPSILON * count as f64 * bound
                {
                    return Err(invalid(
                        "original interface traction spaces must match without truncation",
                    ));
                }
            }
        }
    }
    Ok(map)
}

#[cfg(test)]
mod tests;
