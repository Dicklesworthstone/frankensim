//! Sourced adhesive Coulomb disks driving the complete original modal network.
//!
//! Each ordered original DMT/shear pair supplies its real contact area and one
//! isotropic junction strength. Two complete work columns span its tangent
//! plane. All disks enter one full-network solve; independently clipping two
//! scalar forces would describe a different, non-objective friction law.
//!
//! This explicitly selects equal static and kinetic source strengths. The
//! backward-Euler variational owner retains true set-valued sticking and
//! direction-changing sliding without velocity smoothing. Its bounded internal
//! cell size is a physical accuracy choice: the scheme is first order, and its
//! independently computed numerical loss is neither friction nor heat. Native
//! bodies, spring/damper work ports, pressure, source bounds, numerical loss and
//! work uncertainty all pass before state, diagnostics and clock publish.
//!
//! Indentation, source conditions and the two authored tangent columns remain
//! fixed. This does not resolve normal impacts, rotating frames, velocity
//! weakening, dependent contact reactions or exact continuous event times.

use super::super::adhesive_events::linear_network::vector::{
    MAX_VECTOR_NETWORK_CONTACTS, VectorLinearNetworkConfig, VectorLinearNetworkError,
    VectorLinearNetworkFlow, VectorLinearNetworkStep,
};
use super::super::adhesive_events::native::linear_network::{
    NativeLinearEventError, NativeLinearEventTrajectory, NativeLinearNetworkWorkspace,
    network_parameters, require_linear_network,
};
use super::super::adhesive_events::native::{projected_position, projected_velocity};
use super::super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::{DmtError, DmtResponse, SourceDmtSphereContact};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::det;

/// Work coordinates are physical metres in an orthonormal tangent chart. The
/// two columns map the full mass-normalized native velocity to that chart; the
/// columns themselves need not be orthogonal in mechanical mobility.
pub struct VectorAdhesiveNetworkJunction<'a> {
    pub normal: &'a SourceDmtSphereContact,
    pub shear: &'a SourceAdhesiveShear,
    pub surface_materials: [&'a MaterialStateId; 2],
    pub indentation_m: f64,
    pub tangent_columns: [Vec<f64>; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VectorAdhesiveNetworkModalConfig {
    pub flow: VectorLinearNetworkConfig,
    pub maximum_contact_force_n: f64,
    pub maximum_tangent_speed_m_s: f64,
    /// Sum of physical tangent travel over this sample's internal cells.
    pub maximum_accumulated_slip_m: f64,
}

#[derive(Debug)]
pub enum VectorAdhesiveNetworkModalError {
    Network(ModalCouplingError),
    Normal(DmtError),
    Flow(VectorLinearNetworkError),
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for VectorAdhesiveNetworkModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native vector adhesive friction refused: {self:?}")
    }
}
impl std::error::Error for VectorAdhesiveNetworkModalError {}
impl From<ModalCouplingError> for VectorAdhesiveNetworkModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<DmtError> for VectorAdhesiveNetworkModalError {
    fn from(error: DmtError) -> Self {
        Self::Normal(error)
    }
}
impl From<VectorLinearNetworkError> for VectorAdhesiveNetworkModalError {
    fn from(error: VectorLinearNetworkError) -> Self {
        Self::Flow(error)
    }
}
impl From<NativeLinearEventError> for VectorAdhesiveNetworkModalError {
    fn from(error: NativeLinearEventError) -> Self {
        match error {
            NativeLinearEventError::Network(error) => Self::Network(error),
            NativeLinearEventError::EnergyAccuracy {
                residual_j,
                estimated_error_j,
                tolerance_j,
            } => Self::EnergyAccuracy {
                residual_j,
                estimated_error_j,
                tolerance_j,
            },
        }
    }
}
type Error = VectorAdhesiveNetworkModalError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VectorAdhesiveNetworkModalFrame {
    pub sample: u64,
    pub contact_count: usize,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub friction_dissipation_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    /// Independently computed integrator loss, never a physical heat source.
    pub numerical_dissipation_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub tangent_positions_m: [[f64; 2]; MAX_VECTOR_NETWORK_CONTACTS],
    pub tangent_velocities_m_s: [[f64; 2]; MAX_VECTOR_NETWORK_CONTACTS],
    pub capacities_n: [f64; MAX_VECTOR_NETWORK_CONTACTS],
    pub path: VectorLinearNetworkStep,
}

pub struct VectorAdhesiveNetworkModalSystem<'a> {
    network: CoupledModalSystem,
    junctions: Vec<VectorAdhesiveNetworkJunction<'a>>,
    contacts: Vec<DmtResponse>,
    capacities_n: [f64; MAX_VECTOR_NETWORK_CONTACTS],
    surface_energy_j: f64,
    flow: VectorLinearNetworkFlow,
    workspace: NativeLinearNetworkWorkspace,
    config: VectorAdhesiveNetworkModalConfig,
    last: Option<VectorAdhesiveNetworkModalFrame>,
}

impl<'a> VectorAdhesiveNetworkModalSystem<'a> {
    pub fn new(
        network: CoupledModalSystem,
        junctions: Vec<VectorAdhesiveNetworkJunction<'a>>,
        config: VectorAdhesiveNetworkModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(Some(gate))?;
        require_linear_network(
            &network,
            config.flow.maximum_modes,
            config.flow.maximum_connections,
        )?;
        let count = junctions.len();
        if count == 0 || count > MAX_VECTOR_NETWORK_CONTACTS || count > config.flow.maximum_contacts
        {
            return Err(invalid("bounded nonempty vector adhesive junction set required").into());
        }
        for value in [
            config.maximum_contact_force_n,
            config.maximum_tangent_speed_m_s,
            config.maximum_accumulated_slip_m,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid("positive finite vector-junction limits required").into());
            }
        }
        let mut contacts = Vec::with_capacity(count);
        let mut capacities_n = [0.0; MAX_VECTOR_NETWORK_CONTACTS];
        let mut surface_energy_j = 0.0;
        for (i, site) in junctions.iter().enumerate() {
            poll(Some(gate))?;
            for column in &site.tangent_columns {
                if column.len() != network.mode_count() || column.iter().any(|x| !x.is_finite()) {
                    return Err(
                        invalid("both tangent columns must cover all original modes").into(),
                    );
                }
            }
            let materials = site
                .normal
                .normal_state()
                .elastic()
                .interface()
                .surface_materials();
            if site.surface_materials[0] != &materials[0]
                || site.surface_materials[1] != &materials[1]
                || site.shear.normal_state_identity() != site.normal.normal_state().identity()
            {
                return Err(invalid(
                    "vector junction must retain its original ordered normal/shear sources",
                )
                .into());
            }
            let law = site.shear.law();
            if law.static_shear_strength_pa() != law.kinetic_shear_strength_pa() {
                return Err(invalid("vector variational friction requires equal sourced static and kinetic strengths").into());
            }
            let contact = site.normal.evaluate(site.indentation_m)?;
            capacities_n[i] = finite(contact.contact_area_m2 * law.static_shear_strength_pa())?;
            limit(
                "original vector adhesive capacity",
                capacities_n[i],
                config.maximum_contact_force_n,
            )?;
            surface_energy_j = finite(surface_energy_j + contact.reversible_energy_j)?;
            contacts.push(contact);
        }
        let (stiffness, damping, links) = network_parameters(&network)?;
        let flow = VectorLinearNetworkFlow::new(
            stiffness,
            damping,
            links,
            junctions
                .iter()
                .map(|site| site.tangent_columns.clone())
                .collect(),
            config.flow,
        )?;
        let workspace = NativeLinearNetworkWorkspace::new(&network);
        limit(
            "initial vector adhesive energy magnitude",
            finite(network.total_energy_j()? + surface_energy_j)?.abs(),
            network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(Self {
            network,
            junctions,
            contacts,
            capacities_n,
            surface_energy_j,
            flow,
            workspace,
            config,
            last: None,
        })
    }

    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    pub fn contact_count(&self) -> usize {
        self.junctions.len()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub fn junctions(&self) -> &[VectorAdhesiveNetworkJunction<'a>] {
        &self.junctions
    }
    pub fn normal_responses(&self) -> &[DmtResponse] {
        &self.contacts
    }
    pub const fn config(&self) -> VectorAdhesiveNetworkModalConfig {
        self.config
    }
    pub const fn last_frame(&self) -> Option<&VectorAdhesiveNetworkModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, Error> {
        Ok(finite(
            self.network.total_energy_j()? + self.surface_energy_j,
        )?)
    }
    pub fn step(&mut self, external: &[f64]) -> Result<&VectorAdhesiveNetworkModalFrame, Error> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&VectorAdhesiveNetworkModalFrame, Error> {
        poll(Some(gate))?;
        if external.len() != self.mode_count() || external.iter().any(|f| !f.is_finite()) {
            return Err(
                invalid("vector junction force input must cover every original mode").into(),
            );
        }
        self.workspace.capture(&self.network);
        let count = self.contact_count();
        let path = self.flow.advance(
            self.workspace.positions(),
            self.workspace.velocities(),
            external,
            self.network.dt,
            &self.capacities_n[..count],
            gate,
        )?;
        for i in 0..count {
            limit(
                "individual vector junction travel",
                path.accumulated_slip_m[i],
                self.config.maximum_accumulated_slip_m,
            )?;
            limit(
                "individual vector junction force",
                path.maximum_contact_force_n[i],
                self.config.maximum_contact_force_n,
            )?;
            limit(
                "individual vector junction speed",
                path.maximum_contact_speed_m_s[i],
                self.config.maximum_tangent_speed_m_s,
            )?;
        }
        let native = self.workspace.stage_with_numerical_loss(
            &mut self.network,
            external,
            self.surface_energy_j,
            NativeLinearEventTrajectory {
                positions: self.flow.positions(),
                velocities: self.flow.velocities(),
                connection_dissipation_j: self.flow.connection_dissipation_j(),
                initial_energy_j: path.initial_energy_j,
                modal_energy_j: path.modal_energy_j,
                connection_energy_j: path.connection_energy_j,
                external_work_j: path.external_work_j,
                modal_damping_dissipation_j: path.modal_damping_dissipation_j,
                connection_damping_dissipation_j: path.connection_damping_dissipation_j,
                friction_dissipation_j: path.friction_dissipation_j,
                numerical_work_uncertainty_j: path.numerical_work_uncertainty_j,
            },
            path.numerical_dissipation_j,
            Some(self.flow.connection_numerical_dissipation_j()),
            gate,
        )?;
        let mut tangent_positions_m = [[0.0; 2]; MAX_VECTOR_NETWORK_CONTACTS];
        let mut tangent_velocities_m_s = [[0.0; 2]; MAX_VECTOR_NETWORK_CONTACTS];
        let mut friction_sum = 0.0;
        for (i, site) in self.junctions.iter().enumerate() {
            poll(Some(gate))?;
            let norm_squared = finite(
                dot(&site.tangent_columns[0], &site.tangent_columns[0])?
                    + dot(&site.tangent_columns[1], &site.tangent_columns[1])?,
            )?;
            limit(
                "vector junction passive-energy speed bound",
                finite(det::sqrt(norm_squared) * native.root_energy_bound)?,
                self.config.maximum_tangent_speed_m_s,
            )?;
            for axis in 0..2 {
                tangent_positions_m[i][axis] =
                    projected_position(&self.network.candidates, &site.tangent_columns[axis])?;
                tangent_velocities_m_s[i][axis] =
                    projected_velocity(&self.network.candidates, &site.tangent_columns[axis])?;
            }
            let loss = finite(path.contact_friction_dissipation_j[i])?;
            if loss < 0.0 {
                return Err(invalid("negative vector junction dissipation").into());
            }
            friction_sum = finite(friction_sum + loss)?;
        }
        let residual = finite(friction_sum - path.friction_dissipation_j)?;
        if finite(residual.abs() + path.numerical_work_uncertainty_j)? > native.energy_tolerance_j {
            return Err(Error::EnergyAccuracy {
                residual_j: residual,
                estimated_error_j: path.numerical_work_uncertainty_j,
                tolerance_j: native.energy_tolerance_j,
            });
        }
        let accepted = VectorAdhesiveNetworkModalFrame {
            sample: native.sample,
            contact_count: count,
            observer_pressure_pa: native.observer_pressure_pa,
            network_energy_j: native.network_energy_j,
            surface_energy_j: self.surface_energy_j,
            external_work_j: native.external_work_j,
            friction_dissipation_j: path.friction_dissipation_j,
            modal_damping_dissipation_j: path.modal_damping_dissipation_j,
            connection_damping_dissipation_j: native.connection_damping_dissipation_j,
            numerical_dissipation_j: path.numerical_dissipation_j,
            numerical_work_uncertainty_j: path.numerical_work_uncertainty_j,
            energy_residual_j: native.energy_residual_j,
            energy_tolerance_j: native.energy_tolerance_j,
            tangent_positions_m,
            tangent_velocities_m_s,
            capacities_n: self.capacities_n,
            path,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        Ok(self.last.insert(accepted))
    }
}
