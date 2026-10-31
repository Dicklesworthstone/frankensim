//! Causal evolution of the COMPLETE original resolved fluid/structure equations.
//!
//! The original real mixed mass, viscous, conductive, pressure-work and wall
//! operators are used without modal truncation or a rational approximation.
//! Every original structural displacement is retained. A diagonal change of
//! units balances the energy; the original mass matrix remains the quadratic
//! Hamiltonian. Its inverse is computed by fs-la and checked against EVERY
//! original entry before the admitted fs-phs exact held-input owner is called.
//!
//! Fluid viscosity (shear plus bulk), fluid conduction, molecular heat exchange
//! finite-wall conduction and original structural damping are integrated
//! along that SAME trajectory.
//! Pressure and structural-force work use their actual interval-integrated dual flows.
//! No response curve, fitted decay or endpoint-loss approximation is involved.
//!
//! Original caloric reservoirs contribute their full P2 temperature fields,
//! frozen heat capacity and positive translational/molecular heat exchange. Insulated
//! and isothermal thermal boundaries constrain only the translational field.
//! Finite solid slabs may retain their actual P1 normal heat fields, sharing
//! the original P2 fluid interface; resolution is explicit, with no fitted poles.
//! Geometry and equilibrium remain fixed; no mean flow, ALE,
//! phase change, nonlinear mechanics or acoustic-resolution claim is inferred.

use super::*;
use fs_la::factor::{FactorError, cholesky};
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    prepared::PreparedQuadraticFlow, split_quadratic,
};
use fs_phs::{PhsError, PortHamiltonian, PreparedStepError, QuadraticStorage};

mod caloric;
pub use caloric::ThermoviscousCaloricReservoir;
mod wall;
pub use wall::{ConductiveWallTimeLayout, TangentialWallConductivity, ThermoviscousWallDiscretization};
mod interface;
pub use interface::{
    ConformingThermoviscousInterface, ConformingThermoviscousTimePairStep,
    PreparedConformingThermoviscousTimePair, PreparedThermoviscousTimePair,
    ThermoviscousTimePairStep,
};
const HEAT_CHANNELS: usize = 5;

/// Original-owner refusals; no state is published on failure or cancellation.
#[derive(Debug)]
pub enum ThermoviscousTimeError {
    Fluid(ThermoviscousError),
    Factor(FactorError),
    Propagation(PhsError),
}
impl std::fmt::Display for ThermoviscousTimeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FS-FEEC-LNS-TIME: {self:?}")
    }
}
impl std::error::Error for ThermoviscousTimeError {}
impl From<ThermoviscousError> for ThermoviscousTimeError {
    fn from(e: ThermoviscousError) -> Self {
        Self::Fluid(e)
    }
}
impl From<FactorError> for ThermoviscousTimeError {
    fn from(e: FactorError) -> Self {
        Self::Factor(e)
    }
}
impl From<PhsError> for ThermoviscousTimeError {
    fn from(e: PhsError) -> Self {
        Self::Propagation(e)
    }
}
impl From<PreparedStepError> for ThermoviscousTimeError {
    fn from(e: PreparedStepError) -> Self {
        match e {
            PreparedStepError::Cancelled => ThermoviscousError::Cancelled.into(),
            PreparedStepError::Solver(e) => e.into(),
        }
    }
}
fn invalid(what: &'static str) -> ThermoviscousTimeError {
    ThermoviscousError::InvalidInput(what).into()
}

/// Explicit small-domain cold preparation and physical accounting admission.
#[derive(Clone, Copy, Debug)]
pub struct ThermoviscousTimeLimits {
    /// Complete mixed, molecular and finite-wall fields plus every structural displacement.
    /// Exceeding the cap is a refusal; no original field is truncated.
    pub max_states: usize,
    /// Existing pHS flow-integral byte/work limits. Coordinate maps additionally
    /// use O(max_states^2) storage; no per-step exponential is computed.
    pub preparation: DenseExponentialBudget,
    /// Maximum energy-balance defect relative to endpoint energy/work scale.
    /// This is numerical admission, not a material or mesh accuracy certificate.
    pub maximum_relative_energy_residual: f64,
}
impl Default for ThermoviscousTimeLimits {
    fn default() -> Self {
        Self {
            max_states: 127,
            preparation: DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: 64 * 1024 * 1024,
                max_multiply_terms: 20_000_000_000,
            },
            maximum_relative_energy_residual: 1e-10,
        }
    }
}

/// Real instantaneous fields in their original physical coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoviscousTimeState {
    /// P1 pressure perturbation [Pa], in vertex order.
    pub pressure_pa: Vec<f64>,
    /// Original hierarchical P2 Cartesian velocity coefficients [m/s].
    pub velocity_m_s: Vec<[f64; 3]>,
    /// Original hierarchical P2 temperature PERTURBATION coefficients [K].
    pub temperature_k: Vec<f64>,
    /// Original full P2 molecular-temperature perturbations [K], in the
    /// supplied reservoir order. A zero-capacity reservoir has an empty vector.
    /// Molecular fields have no independent isothermal boundary constraint.
    pub molecular_temperature_k: Vec<Vec<f64>>,
    /// Full finite-solid temperature perturbations [K], in source-region order.
    /// Each vector is depth-major, then the original P2 trace order from
    /// `conductive_wall_layouts()`. The fluid interface and any zero-temperature
    /// backing are shared/constrained coordinates and are not duplicated.
    pub conductive_wall_temperature_k: Vec<Vec<f64>>,
    /// Original mass-normalized structural displacements [m sqrt(kg)].
    pub wall_displacement_m_sqrt_kg: Vec<f64>,
    /// Original mass-normalized structural velocities [m sqrt(kg)/s].
    pub wall_velocity_m_sqrt_kg_per_s: Vec<f64>,
}

/// Same-window energy/work [J]. No loss or balance residual is clipped.
#[derive(Clone, Copy, Debug)]
pub struct ThermoviscousTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    pub pressure_work_j: f64,
    pub structural_work_j: f64,
    /// Signed conductive boundary transfer of acoustic available energy,
    /// integral((T'/T0) q_in dA dt) [J]. This is not the first-order caloric
    /// heat integral(q_in dA dt), nor a volume dissipation channel.
    pub conductive_boundary_work_j: f64,
    /// Combined deviatoric and dilatational viscosity in the original fluid.
    pub fluid_viscous_heat_j: f64,
    pub fluid_conductive_heat_j: f64,
    pub structural_damping_heat_j: f64,
    /// Heat exchange between the translational and molecular temperatures.
    /// This is independent of viscosity and of spatial thermal conduction.
    pub molecular_relaxation_heat_j: f64,
    /// Original finite-solid Fourier conduction, along the same trajectory.
    /// Includes declared in-plane conduction on the spatial wall path.
    pub solid_conductive_heat_j: f64,
}
impl ThermoviscousTimeWork {
    pub fn dissipated_j(self) -> f64 {
        self.fluid_viscous_heat_j
            + self.fluid_conductive_heat_j
            + self.structural_damping_heat_j
            + self.molecular_relaxation_heat_j
            + self.solid_conductive_heat_j
    }
    pub fn supplied_j(self) -> f64 {
        self.pressure_work_j + self.structural_work_j + self.conductive_boundary_work_j
    }
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j + self.dissipated_j() - self.supplied_j()
    }
}

/// One fully admitted candidate. The caller owns the accepted state and clock.
#[derive(Clone, Debug)]
pub struct ThermoviscousTimeStep {
    pub state: ThermoviscousTimeState,
    /// Actual trajectory means, not endpoint flow or sampled pressure.
    pub mean_inward_flow_m3_s: Vec<f64>,
    pub mean_wall_velocity_m_sqrt_kg_per_s: Vec<f64>,
    pub work: ThermoviscousTimeWork,
}

/// One immutable fixed-clock cache of the original complete pHS trajectory.
/// Warm work allocates result vectors and performs bounded matrix contractions.
/// This is not an allocation-free or large-domain sparse time solver.
pub struct PreparedThermoviscousTime<'a, 'c> {
    model: &'a ThermoviscousFe<'c>,
    duration_s: f64,
    scale: Vec<f64>,
    wall_scale: Vec<f64>,
    frozen_cp_j_kg_k: f64,
    reservoirs: Vec<ThermoviscousCaloricReservoir>,
    reservoir_offsets: Vec<Option<usize>>,
    conductive_wall_fields: wall::WallFields,
    // Original physical forms, retained for causal multi-region composition.
    // Never replace these with projected resistance or an energy remainder.
    heat: [Vec<f64>; HEAT_CHANNELS],
    // The SAME checked normalized mixed-mass inverse used to construct J/R/G.
    // Retained so physical P2 interface loads use the original field owner.
    inverse_mass: Vec<f64>,
    system: PortHamiltonian,
    flow: PreparedQuadraticFlow,
    mean_state: Vec<f64>,
    mean_input: Vec<f64>,
    maximum_relative_energy_residual: f64,
}

impl<'c> ThermoviscousFe<'c> {
    /// Prepare all original fields and wall coordinates for real held forces.
    /// Finite slabs require explicit meshes through `prepare_time_with_thermal_fields`.
    pub fn prepare_time<'a>(
        &'a self,
        duration_s: f64,
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousTime<'a, 'c>, ThermoviscousTimeError> {
        self.prepare_time_with_caloric_reservoirs(
            duration_s,
            self.medium.cp_j_kg_k,
            &[],
            limits,
            cancelled,
        )
    }

    /// Evolve the original positive caloric reservoirs on the SAME complete
    /// velocity/pressure/temperature and structural trajectory. Each positive
    /// capacity retains all P2 molecular coefficients, including those at an
    /// isothermal translational wall. No fitted pole or bulk-loss overlay enters.
    ///
    /// The original equilibrium Cp must equal frozen Cp plus the supplied
    /// capacities to roundoff. Source adapters retain species/process identity;
    /// this material-independent owner preserves their order and coefficients.
    pub fn prepare_time_with_caloric_reservoirs<'a>(
        &'a self,
        duration_s: f64,
        frozen_cp_j_kg_k: f64,
        reservoirs: &[ThermoviscousCaloricReservoir],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousTime<'a, 'c>, ThermoviscousTimeError> {
        PreparedThermoviscousTime::new(
            self,
            duration_s,
            frozen_cp_j_kg_k,
            reservoirs,
            None,
            None,
            limits,
            cancelled,
        )
    }

    /// Evolve original fluid, molecular, finite solid and structural coordinates
    /// together. Wall normal elements discretize the ORIGINAL heat equation,
    /// not its frequency response. Refine the physical depth mesh to converge
    /// to the original harmonic slab admittance; no finite count is exact.
    pub fn prepare_time_with_thermal_fields<'a>(
        &'a self,
        duration_s: f64,
        frozen_cp_j_kg_k: f64,
        reservoirs: &[ThermoviscousCaloricReservoir],
        walls: ThermoviscousWallDiscretization<'_>,
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousTime<'a, 'c>, ThermoviscousTimeError> {
        PreparedThermoviscousTime::new(
            self,
            duration_s,
            frozen_cp_j_kg_k,
            reservoirs,
            Some(walls),
            None,
            limits,
            cancelled,
        )
    }
}

impl<'a, 'c> PreparedThermoviscousTime<'a, 'c> {
    fn new(
        model: &'a ThermoviscousFe<'c>,
        duration_s: f64,
        frozen_cp_j_kg_k: f64,
        reservoirs: &[ThermoviscousCaloricReservoir],
        wall_resolution: Option<ThermoviscousWallDiscretization<'_>>,
        tangential: Option<&[TangentialWallConductivity]>,
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousTimeError> {
        poll(cancelled)?;
        let (caloric_end, reservoir_offsets) =
            caloric::layout(model, frozen_cp_j_kg_k, reservoirs)?;
        let (conductive_wall_fields, base) =
            wall::WallFields::new(model, duration_s, wall_resolution, tangential, caloric_end, cancelled)?;
        let walls = model.pressure_start - 3 * model.velocity_count;
        let n = base + walls;
        let ports = model.loads.len() + walls;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !limits.maximum_relative_energy_residual.is_finite()
            || limits.maximum_relative_energy_residual <= 0.0
            || limits.maximum_relative_energy_residual > 1e-6
            || n > limits.max_states
            || 2 * (n + ports + 1) > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        {
            return Err(invalid(
                "complete transient state/port capacity, positive clock and bounded energy tolerance required",
            ));
        }
        let (mut mass, molecular_exchange) = caloric::assemble(
            model,
            frozen_cp_j_kg_k,
            reservoirs,
            &reservoir_offsets,
            base,
            cancelled,
        )?;
        let solid_conduction =
            conductive_wall_fields.assemble(model, &mut mass, base, cancelled)?;
        let mut wall_scale = Vec::with_capacity(walls);
        for i in 0..walls {
            let row = 3 * model.velocity_count + i;
            mass[row * base + row] += 1.0;
            let omega = model
                .moving_wall
                .as_ref()
                .expect("wall source")
                .structure
                .modes()[i]
                .angular_frequency_rad_s;
            wall_scale.push(if omega > 0.0 { omega } else { 1.0 });
        }
        let mut scale = Vec::with_capacity(base);
        for i in 0..base {
            let value = det::sqrt(mass[i * base + i]);
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid("positive representable mixed storage required"));
            }
            scale.push(value);
        }
        for i in 0..base {
            poll(cancelled)?;
            for j in 0..base {
                mass[i * base + j] /= scale[i] * scale[j];
            }
        }
        // Retain the full original normalized storage. Do not replace it with
        // an identity, symmetrize it, or silently discard small thermal terms.
        let mut q = vec![0.0; n * n];
        for i in 0..base {
            q[i * n..i * n + base].copy_from_slice(&mass[i * base..(i + 1) * base]);
        }
        let factor = cholesky(&mass, base)?;
        poll(cancelled)?;
        let mut inverse = vec![0.0; n * n];
        for col in 0..base {
            poll(cancelled)?;
            let mut rhs = vec![0.0; base];
            rhs[col] = 1.0;
            factor.solve(&mut rhs);
            for row in 0..base {
                inverse[row * n + col] = rhs[row];
            }
        }
        // Cholesky reads the lower triangle. Checking the full inverse against
        // ALL original entries makes an inconsistent upper triangle a refusal.
        for i in 0..base {
            poll(cancelled)?;
            for j in 0..base {
                let mut value = 0.0;
                let mut bound = 0.0;
                for k in 0..base {
                    let term = mass[i * base + k] * inverse[k * n + j];
                    value += term;
                    bound += term.abs();
                }
                let expected: f64 = if i == j { 1.0 } else { 0.0 };
                if !value.is_finite()
                    || (value - expected).abs()
                        > 512.0 * f64::EPSILON * base as f64 * bound.max(expected)
                {
                    return Err(invalid(
                        "original full mixed mass inverse failed its residual",
                    ));
                }
            }
        }
        let mut physical_j = vec![0.0; n * n];
        let mut heat: [Vec<f64>; HEAT_CHANNELS] = std::array::from_fn(|_| vec![0.0; n * n]);
        let mut input = vec![0.0; n * ports];
        for i in 0..model.n {
            poll(cancelled)?;
            for j in 0..model.n {
                let value = model.constant[i * model.n + j] / (scale[i] * scale[j]);
                if i < model.pressure_start && j < model.pressure_start {
                    heat[0][i * n + j] = value;
                } else if i >= model.temperature_start && j >= model.temperature_start {
                    heat[1][i * n + j] = value;
                } else {
                    physical_j[i * n + j] = -value;
                }
            }
            for (p, load) in model.loads.iter().enumerate() {
                input[i * ports + p] = load[i] / scale[i];
            }
        }
        for i in 0..base {
            poll(cancelled)?;
            for j in 0..base {
                heat[3][i * n + j] = molecular_exchange[i * base + j] / (scale[i] * scale[j]);
                heat[4][i * n + j] = solid_conduction[i * base + j] / (scale[i] * scale[j]);
            }
        }
        for i in 0..walls {
            let row = 3 * model.velocity_count + i;
            let mode = &model
                .moving_wall
                .as_ref()
                .expect("wall source")
                .structure
                .modes()[i];
            let displacement = base + i;
            inverse[displacement * n + displacement] = 1.0;
            q[displacement * n + displacement] = if mode.angular_frequency_rad_s > 0.0 {
                1.0
            } else {
                0.0
            };
            let coupling = wall_scale[i] / scale[row];
            physical_j[row * n + displacement] = -coupling;
            physical_j[displacement * n + row] = coupling;
            heat[2][row * n + row] =
                2.0 * mode.damping_ratio * mode.angular_frequency_rad_s / (scale[row] * scale[row]);
            input[row * ports + model.loads.len() + i] = 1.0 / scale[row];
        }
        let j = congruence(&physical_j, &inverse, n, cancelled)?;
        let total: Vec<_> = (0..n * n)
            .map(|i| heat.iter().map(|channel| channel[i]).sum())
            .collect();
        let r = congruence(&total, &inverse, n, cancelled)?;
        let mut g = vec![0.0; n * ports];
        for i in 0..n {
            poll(cancelled)?;
            for p in 0..ports {
                g[i * ports + p] = (0..n)
                    .map(|k| inverse[k * n + i] * input[k * ports + p])
                    .sum();
            }
        }
        let storage = QuadraticStorage::new(q, n)?;
        let system = PortHamiltonian::new(n, ports, j, r, g, Box::new(storage))?;
        poll(cancelled)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: limits.max_states,
                max_omega_dt: 1.0,
            },
        )?;
        let observations: Vec<_> = (0..ports)
            .map(|p| (0..n).map(|i| input[i * ports + p]).collect::<Vec<_>>())
            .collect();
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &heat,
            &observations,
            limits.preparation,
            &mut *cancelled,
        )?;
        poll(cancelled)?;
        // Exact original observation maps from this SAME cached trajectory.
        // Keeping the complete state columns includes caloric, conductive-wall
        // and structural memory; no probe solve or reduced response is used.
        let mut mean_state = Vec::with_capacity(ports * n);
        let mut mean_input = Vec::with_capacity(ports * ports);
        for port in 0..ports {
            poll(cancelled)?;
            let (state, input) = flow.linear_integral_map(HEAT_CHANNELS + port)?;
            mean_state.extend(state.into_iter().map(|value| value / duration_s));
            mean_input.extend(input.into_iter().map(|value| value / duration_s));
        }
        if mean_state.iter().chain(&mean_input).any(|value| !value.is_finite()) {
            return Err(invalid("original transient mean-port map is nonfinite"));
        }
        Ok(Self {
            model,
            duration_s,
            scale,
            wall_scale,
            frozen_cp_j_kg_k,
            reservoirs: reservoirs.to_vec(),
            reservoir_offsets,
            conductive_wall_fields,
            heat,
            inverse_mass: inverse,
            system,
            flow,
            mean_state,
            mean_input,
            maximum_relative_energy_residual: limits.maximum_relative_energy_residual,
        })
    }

    pub fn inner(&self) -> &ThermoviscousFe<'c> {
        self.model
    }
    /// Original admitted system; no source coefficients can be edited in place.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    pub const fn frozen_cp_j_kg_k(&self) -> f64 {
        self.frozen_cp_j_kg_k
    }

    /// Unchanged capacities and relaxation times in the original source order,
    /// including zero-capacity entries (which own no dynamic coordinates).
    pub fn caloric_reservoirs(&self) -> &[ThermoviscousCaloricReservoir] {
        &self.reservoirs
    }

    /// Full original surface and normal-node ordering, in source-region order.
    pub fn conductive_wall_layouts(&self) -> &[ConductiveWallTimeLayout] {
        &self.conductive_wall_fields.layouts
    }

    pub fn maximum_wall_frequency_hz(&self) -> Option<f64> {
        self.conductive_wall_fields.maximum_frequency_hz
    }

    /// Explicit all-rest initial condition, including structural coordinates.
    /// The structural source's separate time state is not silently advanced.
    pub fn zero_state(&self) -> ThermoviscousTimeState {
        self.decode(&vec![0.0; self.system.state_dim()])
    }

    pub fn stored_energy_j(
        &self,
        state: &ThermoviscousTimeState,
    ) -> Result<f64, ThermoviscousTimeError> {
        let x = self.encode(state)?;
        Ok(self.system.hamiltonian(&x))
    }

    /// Exact same-window relation `inward_flow = free + mobility * traction`
    /// for one original pressure port. All other held pressures, wall forces
    /// and complete initial fields contribute to `free`. The selected traction
    /// slot is ignored, so an owning junction can solve it before propagation.
    /// This contracts cached original-coordinate maps; it runs no trial step,
    /// exponential, constitutive fit or physical state mutation.
    pub fn average_relation(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        selected_port: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ThermoviscousTimeError> {
        poll(cancelled)?;
        let pressure_ports = self.model.loads.len();
        if selected_port >= pressure_ports
            || traction_pa.len() != pressure_ports
            || wall_forces_n_per_sqrt_kg.len() != self.wall_scale.len()
            || traction_pa.iter().enumerate().any(|(p, value)| p != selected_port && !value.is_finite())
            || wall_forces_n_per_sqrt_kg.iter().any(|value| !value.is_finite())
        {
            return Err(invalid("complete original fields and finite other held drives required"));
        }
        let state = self.encode(initial)?;
        let n = self.system.state_dim();
        let ports = pressure_ports + self.wall_scale.len();
        let mut free: f64 = self.mean_state[selected_port * n..(selected_port + 1) * n]
            .iter().zip(state).map(|(a, x)| a * x).sum();
        let row = &self.mean_input[selected_port * ports..(selected_port + 1) * ports];
        free += row.iter().zip(traction_pa.iter().chain(wall_forces_n_per_sqrt_kg))
            .enumerate().filter(|(p, _)| *p != selected_port)
            .map(|(_, (a, u))| a * u).sum::<f64>();
        if !free.is_finite() {
            return Err(invalid("original transient free mean flow is nonfinite"));
        }
        poll(cancelled)?;
        Ok((free, row[selected_port]))
    }

    /// Complete same-window relation `Q = free + mobility * pressure` for the
    /// selected original pressure ports, in the supplied order. `mobility` is
    /// row-major and retains EVERY cross-port term from the original trajectory.
    /// The selected traction slots are ignored; all remaining tractions,
    /// structural forces and complete initial fields contribute to `free`.
    ///
    /// The initial field is encoded once. These are cached original-coordinate
    /// contractions, without probe steps, factorization, response fitting or
    /// mutation. In particular, a connected duct must not replace this matrix
    /// by its diagonal when another duct loads the same spatial field.
    pub fn average_relations(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        selected_ports: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), ThermoviscousTimeError> {
        poll(cancelled)?;
        let pressure_ports = self.model.loads.len();
        if selected_ports.is_empty() || selected_ports.len() > pressure_ports
            || traction_pa.len() != pressure_ports
            || wall_forces_n_per_sqrt_kg.len() != self.wall_scale.len()
        {
            return Err(invalid("nonempty original pressure selection and complete held drives required"));
        }
        let mut selected = vec![false; pressure_ports];
        for &port in selected_ports {
            if port >= pressure_ports || selected[port] {
                return Err(invalid("selected original pressure ports must be valid and unique"));
            }
            selected[port] = true;
        }
        if traction_pa.iter().enumerate().any(|(port, value)| !selected[port] && !value.is_finite())
            || wall_forces_n_per_sqrt_kg.iter().any(|value| !value.is_finite())
        {
            return Err(invalid("complete original fields and finite other held drives required"));
        }
        let state = self.encode(initial)?;
        let n = self.system.state_dim();
        let ports = pressure_ports + self.wall_scale.len();
        let entries = selected_ports.len().checked_mul(selected_ports.len())
            .ok_or_else(|| invalid("original pressure response storage overflow"))?;
        let mut free = Vec::with_capacity(selected_ports.len());
        let mut mobility = Vec::with_capacity(entries);
        for &port in selected_ports {
            poll(cancelled)?;
            let mut value: f64 = self.mean_state[port * n..(port + 1) * n]
                .iter().zip(&state).map(|(a, x)| a * x).sum();
            let row = &self.mean_input[port * ports..(port + 1) * ports];
            value += row.iter().zip(traction_pa.iter().chain(wall_forces_n_per_sqrt_kg))
                .enumerate()
                .filter(|(input, _)| *input >= pressure_ports || !selected[*input])
                .map(|(_, (a, u))| a * u).sum::<f64>();
            if !value.is_finite() {
                return Err(invalid("original transient free mean flow is nonfinite"));
            }
            free.push(value);
            mobility.extend(selected_ports.iter().map(|&other| row[other]));
        }
        poll(cancelled)?;
        Ok((free, mobility))
    }

    fn encode(&self, state: &ThermoviscousTimeState) -> Result<Vec<f64>, ThermoviscousTimeError> {
        let model = self.model;
        let walls = self.wall_scale.len();
        if state.pressure_pa.len() != model.temperature_start - model.pressure_start
            || state.velocity_m_s.len() != model.space.ndof
            || state.temperature_k.len() != model.space.ndof
            || state.molecular_temperature_k.len() != self.reservoir_offsets.len()
            || state.wall_displacement_m_sqrt_kg.len() != walls
            || state.wall_velocity_m_sqrt_kg_per_s.len() != walls
            || state
                .pressure_pa
                .iter()
                .chain(&state.temperature_k)
                .chain(&state.wall_displacement_m_sqrt_kg)
                .chain(&state.wall_velocity_m_sqrt_kg_per_s)
                .chain(state.velocity_m_s.iter().flatten())
                .chain(state.molecular_temperature_k.iter().flatten())
                .chain(state.conductive_wall_temperature_k.iter().flatten())
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite original real field and wall state required",
            ));
        }
        let mut x = vec![0.0; self.system.state_dim()];
        for (i, map) in model.velocity_map.iter().enumerate() {
            for c in 0..3 {
                if let Some(j) = map {
                    x[c * model.velocity_count + j] = state.velocity_m_s[i][c];
                } else {
                    let mut expected = 0.0;
                    let mut scale = state.velocity_m_s[i][c].abs();
                    if let Some(wall) = &model.moving_wall {
                        for (shape, &u) in wall
                            .velocity_shapes
                            .iter()
                            .zip(&state.wall_velocity_m_sqrt_kg_per_s)
                        {
                            let value = shape[i][c] * u;
                            expected += value;
                            scale += value.abs();
                        }
                    }
                    if (state.velocity_m_s[i][c] - expected).abs() > 64.0 * f64::EPSILON * scale {
                        return Err(invalid(
                            "initial real velocity violates the original no-slip wall trace",
                        ));
                    }
                }
            }
        }
        x[model.pressure_start..model.temperature_start].copy_from_slice(&state.pressure_pa);
        for (i, map) in model.temperature_map.iter().enumerate() {
            if let Some(j) = map {
                x[model.temperature_start + j] =
                    state.temperature_k[i] / model.medium.temperature_k;
            } else if state.temperature_k[i] != 0.0 {
                return Err(invalid(
                    "initial temperature violates an isothermal constraint",
                ));
            }
        }
        for i in 0..walls {
            x[3 * model.velocity_count + i] = state.wall_velocity_m_sqrt_kg_per_s[i];
            x[self.scale.len() + i] = state.wall_displacement_m_sqrt_kg[i] * self.wall_scale[i];
        }
        for (field, offset) in state
            .molecular_temperature_k
            .iter()
            .zip(&self.reservoir_offsets)
        {
            match offset {
                Some(start) if field.len() == model.space.ndof => {
                    for (i, temperature) in field.iter().enumerate() {
                        x[start + i] = temperature / model.medium.temperature_k;
                    }
                }
                None if field.is_empty() => {}
                _ => {
                    return Err(invalid(
                        "complete unconstrained molecular temperature fields required",
                    ));
                }
            }
        }
        self.conductive_wall_fields
            .encode(model, &state.conductive_wall_temperature_k, &mut x)?;
        for (value, scale) in x.iter_mut().zip(&self.scale) {
            *value *= scale;
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(invalid("normalized initial state overflow"));
        }
        Ok(x)
    }

    fn decode(&self, state: &[f64]) -> ThermoviscousTimeState {
        let model = self.model;
        let mixed: Vec<_> = state.iter().zip(&self.scale).map(|(x, s)| x / s).collect();
        let wall_velocity = mixed[3 * model.velocity_count..model.pressure_start].to_vec();
        let velocity = model
            .velocity_map
            .iter()
            .enumerate()
            .map(|(i, map)| {
                std::array::from_fn(|c| {
                    map.map_or(0.0, |j| mixed[c * model.velocity_count + j])
                        + model.moving_wall.as_ref().map_or(0.0, |wall| {
                            wall.velocity_shapes
                                .iter()
                                .zip(&wall_velocity)
                                .map(|(phi, u)| phi[i][c] * u)
                                .sum::<f64>()
                        })
                })
            })
            .collect();
        ThermoviscousTimeState {
            pressure_pa: mixed[model.pressure_start..model.temperature_start].to_vec(),
            velocity_m_s: velocity,
            temperature_k: model
                .temperature_map
                .iter()
                .map(|map| {
                    map.map_or(0.0, |j| {
                        mixed[model.temperature_start + j] * model.medium.temperature_k
                    })
                })
                .collect(),
            molecular_temperature_k: self
                .reservoir_offsets
                .iter()
                .map(|offset| {
                    offset.map_or_else(Vec::new, |start| {
                        mixed[start..start + model.space.ndof]
                            .iter()
                            .map(|theta| theta * model.medium.temperature_k)
                            .collect()
                    })
                })
                .collect(),
            conductive_wall_temperature_k: self.conductive_wall_fields.decode(model, &mixed),
            wall_displacement_m_sqrt_kg: self
                .wall_scale
                .iter()
                .enumerate()
                .map(|(i, s)| state[self.scale.len() + i] / s)
                .collect(),
            wall_velocity_m_sqrt_kg_per_s: wall_velocity,
        }
    }

    /// One exact held-input candidate through the existing pHS propagator.
    /// Pressure tractions [Pa] precede generalized wall forces [N/sqrt(kg)].
    /// All source fields, initial state and accepted clock remain unchanged
    /// on any refusal, including cancellation after candidate propagation.
    pub fn step(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousTimeStep, ThermoviscousTimeError> {
        poll(cancelled)?;
        if traction_pa.len() != self.model.loads.len()
            || wall_forces_n_per_sqrt_kg.len() != self.wall_scale.len()
            || traction_pa
                .iter()
                .chain(wall_forces_n_per_sqrt_kg)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite held traction and structural drives required",
            ));
        }
        let input: Vec<_> = traction_pa
            .iter()
            .chain(wall_forces_n_per_sqrt_kg)
            .copied()
            .collect();
        let x = self.encode(initial)?;
        let candidate = self.flow.evaluate(&x, &input)?;
        poll(cancelled)?;
        let initial_energy_j = self.system.hamiltonian(&x);
        let final_energy_j = self.system.hamiltonian(&candidate.state);
        let delta_energy_j = energy_change_j(&self.system, &x, &candidate.state);
        let pressure_work_j = traction_pa
            .iter()
            .enumerate()
            .map(|(p, u)| u * candidate.integrals[HEAT_CHANNELS + p])
            .sum();
        let structural_work_j = wall_forces_n_per_sqrt_kg
            .iter()
            .enumerate()
            .map(|(p, u)| u * candidate.integrals[HEAT_CHANNELS + self.model.loads.len() + p])
            .sum();
        let work = ThermoviscousTimeWork {
            initial_energy_j,
            final_energy_j,
            delta_energy_j,
            pressure_work_j,
            structural_work_j,
            conductive_boundary_work_j: 0.0,
            fluid_viscous_heat_j: candidate.integrals[0],
            fluid_conductive_heat_j: candidate.integrals[1],
            structural_damping_heat_j: candidate.integrals[2],
            molecular_relaxation_heat_j: candidate.integrals[3],
            solid_conductive_heat_j: candidate.integrals[4],
        };
        admit_work(work, &candidate.integrals, self.maximum_relative_energy_residual)?;
        let means: Vec<_> = candidate.integrals[HEAT_CHANNELS..]
            .iter()
            .map(|v| v / self.duration_s)
            .collect();
        let state = self.decode(&candidate.state);
        // Recheck physical finiteness and the exact retained constraints before
        // returning anything the caller could publish as an accepted endpoint.
        self.encode(&state)?;
        poll(cancelled)?;
        Ok(ThermoviscousTimeStep {
            state,
            mean_inward_flow_m3_s: means[..self.model.loads.len()].to_vec(),
            mean_wall_velocity_m_sqrt_kg_per_s: means[self.model.loads.len()..].to_vec(),
            work,
        })
    }
}

// Symmetric bilinear energy difference, without cancellation of large endpoints.
fn energy_change_j(system: &PortHamiltonian, initial: &[f64], next: &[f64]) -> f64 {
    let mid: Vec<_> = initial.iter().zip(next).map(|(a, b)| 0.5 * a + 0.5 * b).collect();
    let delta: Vec<_> = initial.iter().zip(next).map(|(a, b)| b - a).collect();
    let effort = system.effort(&mid);
    let delta_effort = system.effort(&delta);
    delta.iter().zip(&effort).zip(mid.iter().zip(&delta_effort))
        .map(|((d, e), (m, de))| 0.5 * d * e + 0.5 * m * de).sum()
}

fn admit_work(work: ThermoviscousTimeWork, integrals: &[f64], relative: f64)
    -> Result<(), ThermoviscousTimeError>
{
    let scale = work.initial_energy_j.abs().max(work.final_energy_j.abs())
        .max(work.supplied_j().abs()).max(work.dissipated_j().abs()).max(f64::MIN_POSITIVE);
    let tolerance = relative * scale;
    if [work.initial_energy_j, work.final_energy_j, work.balance_residual_j(), scale]
        .iter().any(|v| !v.is_finite())
        || integrals.iter().any(|v| !v.is_finite())
        || integrals[..HEAT_CHANNELS].iter().any(|v| *v < -tolerance)
        || work.balance_residual_j().abs() > tolerance
    {
        return Err(invalid("original transient field/structure energy balance refused"));
    }
    Ok(())
}

// Standard congruence T^T A T; never average or repair the source matrix.
fn congruence(
    a: &[f64],
    t: &[f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, ThermoviscousTimeError> {
    let mut right = vec![0.0; n * n];
    for i in 0..n {
        poll(cancelled)?;
        for k in 0..n {
            if a[i * n + k] != 0.0 {
                for j in 0..n {
                    right[i * n + j] += a[i * n + k] * t[k * n + j];
                }
            }
        }
    }
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        poll(cancelled)?;
        for k in 0..n {
            if t[k * n + i] != 0.0 {
                for j in 0..n {
                    result[i * n + j] += t[k * n + i] * right[k * n + j];
                }
            }
        }
    }
    Ok(result)
}
