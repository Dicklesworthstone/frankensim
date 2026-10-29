//! Actual source tensor memory in the existing exact pressure/contact owner.
use super::super::exact::{CoupledArm, ExactModes};
use super::super::relaxation::{InitialModalMemory, ModalRelaxation};
use super::super::{ApertureMode, MAX_MODES, ModalAperture, ModalApertureSpec, invalid};
use super::{MaterialSolidOptions, poll, scaled, source_error};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::paired::{BladePropagator, BladeState};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;
use fs_solid::reduce::visco::HereditaryModalValve;

/// Runtime clock/ports and an explicit retained-frequency sampling ceiling.
#[derive(Clone, Copy, Debug)]
pub struct HereditarySolidOptions {
    pub aperture: MaterialSolidOptions,
    /// Maximum instantaneous retained angular frequency times dt, in (0, pi).
    /// Exact linear propagation still requires resolved sampled port histories.
    pub maximum_angular_step: f64,
}

/// Original equilibrium modes, tensor laws and physical nodal observation.
/// Accepted history is owned by the aperture's SAME ModalRelaxation, not here.
pub struct HereditaryModalSolid {
    source: HereditaryModalValve,
    root_mass: f64,
}
impl HereditaryModalSolid {
    pub fn source(&self) -> &HereditaryModalValve {
        &self.source
    }
    pub fn coordinate_mass_kg(&self) -> f64 {
        self.source.total_mass_kg()
    }
    pub fn mass_normalized_coordinates(
        &self,
        coordinates: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        poll(gate)?;
        if coordinates.len() != self.source.modes().len()
            || coordinates.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "hereditary observation needs all finite physical coordinates",
            ));
        }
        let result = coordinates
            .iter()
            .map(|x| scaled(*x, self.root_mass))
            .collect::<Result<Vec<_>, _>>()?;
        poll(gate)?;
        Ok(result)
    }
    pub fn expand_nodal(
        &self,
        coordinates: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, AcousticRealizeError> {
        let q = self.mass_normalized_coordinates(coordinates, gate)?;
        self.source
            .expand_nodal(&q, gate)
            .map_err(|e| source_error(e, gate))
    }
    pub(in crate::bernoulli_aperture::modal) fn validate_state(
        &self,
        states: &[BladeState],
    ) -> Result<(), AcousticRealizeError> {
        let mut displacement = [0.0; MAX_MODES];
        for (i, state) in states.iter().enumerate() {
            displacement[i] = state.displacement_m;
        }
        self.source
            .validate_displacement_m(&displacement[..states.len()])
            .map_err(super::super::nonlinear)
    }
}
impl ModalAperture {
    /// Consume all equilibrium modes and every original tensor-Prony/Kelvin
    /// term in one exact structural block. q=sqrt(total_mass)*x; pressure,
    /// force, swept-volume, gap and nodal observations use that same map.
    ///
    /// Relaxed/unrelaxed specifies actual initial viscous memory. Equal-time
    /// source channels have an energy-preserving orthogonal compression; no
    /// source time, weak damping direction or coupled mode is discarded.
    /// Kelvin and Maxwell heat are independently integrated on the trajectory
    /// by the original exact owner and appear in the existing energy ledger.
    ///
    /// Requires dt/tau<=2 for every active source time and the explicit
    /// instantaneous frequency ceiling. Original source strain limits are
    /// checked at initial and candidate endpoints; within-sample strain,
    /// sampled pressure/contact, modal and spatial refinement remain caller
    /// obligations. This does not reconstruct per-element heat deposition.
    #[allow(clippy::too_many_arguments)]
    pub fn from_hereditary_material_solid(
        source: HereditaryModalValve,
        force_patterns_per_n: [Option<&[[f64; 3]]>; 2],
        gas: GasState,
        options: HereditarySolidOptions,
        state: Vec<BladeState>,
        initial: InitialModalMemory,
        contact: Obstacle,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        gas.validate_coherence().map_err(super::super::nonlinear)?;
        let n = source.modes().len();
        let channels = source.memory_channel_count();
        let dt = options.aperture.time_step_s;
        if n == 0
            || n > MAX_MODES
            || state.len() != n
            || 3 * n + channels > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2
            || 5 * n + channels > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        {
            return Err(invalid(
                "complete hereditary source exceeds the exact aperture state capacity",
            ));
        }
        if !dt.is_finite()
            || dt <= 0.0
            || !options.maximum_angular_step.is_finite()
            || options.maximum_angular_step <= 0.0
            || options.maximum_angular_step >= core::f64::consts::PI
            || source.instantaneous_frequency_bound_rad_s() * dt > options.maximum_angular_step
            || source
                .branches()
                .iter()
                .any(|b| dt / b.relaxation_time_s() > 2.0)
        {
            return Err(invalid(
                "hereditary source is unresolved by its original sampled clock",
            ));
        }
        let mass = source.total_mass_kg();
        let root_mass = mass.sqrt();
        if !root_mass.is_finite() || root_mass <= 0.0 {
            return Err(invalid("hereditary coordinate mass is unrepresentable"));
        }
        let mut force = [vec![0.0; n], vec![0.0; n]];
        for (port, pattern) in force_patterns_per_n.into_iter().enumerate() {
            if let Some(pattern) = pattern {
                force[port] = source
                    .project_nodal_force(pattern, gate)
                    .map_err(|e| source_error(e, gate))?;
            }
        }
        let mut modes = Vec::with_capacity(n);
        for (i, mode) in source.modes().iter().enumerate() {
            poll(gate)?;
            modes.push(ApertureMode {
                mass_kg: mass,
                stiffness_n_m: scaled(mode.lambda, mass)?,
                damping_n_s_m: 0.0,
                pressure_area_m2: scaled(
                    source.pressure_projection_m2_per_sqrt_kg()[i],
                    root_mass,
                )?,
                gap_projection: scaled(source.opening_projection_per_sqrt_kg()[i], root_mass)?,
                force_projection: [
                    scaled(force[0][i], root_mass)?,
                    scaled(force[1][i], root_mass)?,
                ],
            });
        }
        let viscosity = source.viscosity_n_s_m().to_vec();
        let spec = ModalApertureSpec {
            rest_opening_m: source.rest_opening_m(),
            width_m: source.width_m(),
            modes,
            density_kg_m3: gas.density,
            impedance_pa_s_m3: options.aperture.impedance_pa_s_m3,
            time_step_s: dt,
            max_steps: options.aperture.max_steps,
            limits: options.aperture.limits,
        };
        let mut valve = Self::new_with_propagator(spec, state, contact, BladePropagator::ExactZoh)?
            .with_viscosity(viscosity.clone(), gate)?;
        let mut arms = Vec::with_capacity(channels);
        for branch in source.branches() {
            poll(gate)?;
            for row in branch.rows() {
                arms.push(CoupledArm {
                    coupling: row.clone(),
                    relaxation_time_s: branch.relaxation_time_s(),
                });
            }
        }
        let exact = ExactModes::with_coupled_memory(&valve.spec, &viscosity, arms.clone())?;
        let relaxation = ModalRelaxation::from_coupled_tensor(arms, &valve.state, initial)?;
        let energy =
            valve.energy(&valve.state) + relaxation.energy(&valve.state, &relaxation.states);
        if !energy.is_finite() || energy > valve.spec.limits.maximum_energy_j {
            return Err(invalid(
                "initial hereditary material memory exceeds the original energy budget",
            ));
        }
        valve.hereditary_solid = Some(HereditaryModalSolid { source, root_mass });
        valve.relaxation = Some(relaxation);
        valve.exact = Some(exact);
        valve.material_gas = Some(gas);
        valve.validate_state(&valve.state)?;
        poll(gate)?;
        Ok(valve)
    }

    /// Source tensor laws, original modes and work-conjugate physical mappings.
    /// Accepted internal coordinates are in material_relaxation().memory_sqrt_j().
    pub fn hereditary_solid(&self) -> Option<&HereditaryModalSolid> {
        self.hereditary_solid.as_ref()
    }
}
