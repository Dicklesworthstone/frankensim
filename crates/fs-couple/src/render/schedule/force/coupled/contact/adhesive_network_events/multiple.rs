//! Simultaneous sourced adhesive junctions on the complete native network.
//!
//! Every junction retains its original ordered normal/shear sources and a
//! complete scalar work column. The joint event owner resolves coupled holding
//! reactions, stops and releases; sequential pairwise mechanics is not used.
//! Original component admission, pressure, connection ports and energy gates
//! use the same transaction as the single-junction owner. Accepted bodies,
//! junction phases, diagnostics and clock publish together after the last gate.
//!
//! This admits a bounded independent set of fixed scalar columns, held normal
//! indentations and fixed source conditions. It does not infer vector Coulomb
//! disks from two scalar contacts, update normal geometry, or select reactions
//! in a rank-deficient constraint set. Such input is explicitly refused.

use super::super::adhesive_events::linear_network::LinearNetworkPhase;
use super::super::adhesive_events::linear_network::multiple::{
    MAX_LINEAR_NETWORK_CONTACTS, MultipleLinearNetworkEventConfig, MultipleLinearNetworkEventFlow,
    MultipleLinearNetworkEventStep,
};
use super::super::adhesive_events::native::linear_network::{
    NativeLinearEventTrajectory, NativeLinearNetworkWorkspace, network_parameters,
    require_linear_network,
};
use super::super::adhesive_events::native::{projected_position, projected_velocity};
use super::super::{CoupledModalSystem, dot, finite, invalid, limit, poll};
use super::AdhesiveNetworkEventModalError;
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::{DmtResponse, SourceDmtSphereContact};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::det;

type Error = AdhesiveNetworkEventModalError;

/// One original, ordered physical junction. The work column maps all native
/// mass-normalized velocities to its scalar relative tangential velocity.
pub struct AdhesiveNetworkEventJunction<'a> {
    pub normal: &'a SourceDmtSphereContact,
    pub shear: &'a SourceAdhesiveShear,
    pub surface_materials: [&'a MaterialStateId; 2],
    pub indentation_m: f64,
    pub column: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleAdhesiveNetworkEventModalConfig {
    pub events: MultipleLinearNetworkEventConfig,
    /// Each original junction must satisfy these physical limits separately.
    pub maximum_contact_force_n: f64,
    pub maximum_tangent_speed_m_s: f64,
    pub maximum_accumulated_slip_m: f64,
}

/// Per-junction entries are meaningful only up to `contact_count`. Work and
/// losses remain global physical sums, with per-contact values in `path`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleAdhesiveNetworkEventModalFrame {
    pub sample: u64,
    pub contact_count: usize,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub friction_dissipation_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub tangent_positions_m: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub tangent_velocities_m_s: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub static_capacities_n: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub kinetic_capacities_n: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub path: MultipleLinearNetworkEventStep,
}

pub struct MultipleAdhesiveNetworkEventModalSystem<'a> {
    network: CoupledModalSystem,
    junctions: Vec<AdhesiveNetworkEventJunction<'a>>,
    contacts: Vec<DmtResponse>,
    surface_energy_j: f64,
    static_capacities_n: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    kinetic_capacities_n: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    phases: [LinearNetworkPhase; MAX_LINEAR_NETWORK_CONTACTS],
    flow: MultipleLinearNetworkEventFlow,
    workspace: NativeLinearNetworkWorkspace,
    config: MultipleAdhesiveNetworkEventModalConfig,
    last: Option<MultipleAdhesiveNetworkEventModalFrame>,
}

impl<'a> MultipleAdhesiveNetworkEventModalSystem<'a> {
    /// Attach every junction to the unchanged sample-zero native network.
    /// Admission uses each source's original real area and ordered identities;
    /// later mechanical solves cannot replace a source or silently merge sites.
    pub fn new(
        network: CoupledModalSystem,
        junctions: Vec<AdhesiveNetworkEventJunction<'a>>,
        config: MultipleAdhesiveNetworkEventModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(Some(gate))?;
        require_linear_network(
            &network,
            config.events.network.maximum_modes,
            config.events.network.maximum_connections,
        )?;
        let count = junctions.len();
        if count == 0
            || count > MAX_LINEAR_NETWORK_CONTACTS
            || count > config.events.maximum_contacts
        {
            return Err(
                invalid("bounded nonempty set of original adhesive junctions required").into(),
            );
        }
        for value in [
            config.maximum_contact_force_n,
            config.maximum_tangent_speed_m_s,
            config.maximum_accumulated_slip_m,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid("positive finite per-junction limits required").into());
            }
        }
        let mut contacts = Vec::with_capacity(count);
        let mut static_capacities_n = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut kinetic_capacities_n = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut surface_energy_j = 0.0;
        for (i, junction) in junctions.iter().enumerate() {
            poll(Some(gate))?;
            if junction.column.len() != network.mode_count()
                || junction.column.iter().any(|b| !b.is_finite())
            {
                return Err(invalid("each adhesive column must cover every original mode").into());
            }
            let materials = junction
                .normal
                .normal_state()
                .elastic()
                .interface()
                .surface_materials();
            if junction.surface_materials[0] != &materials[0]
                || junction.surface_materials[1] != &materials[1]
                || junction.shear.normal_state_identity()
                    != junction.normal.normal_state().identity()
            {
                return Err(invalid(
                    "each adhesive junction must retain its original ordered normal/shear source",
                )
                .into());
            }
            let contact = junction.normal.evaluate(junction.indentation_m)?;
            static_capacities_n[i] =
                finite(contact.contact_area_m2 * junction.shear.law().static_shear_strength_pa())?;
            kinetic_capacities_n[i] =
                finite(contact.contact_area_m2 * junction.shear.law().kinetic_shear_strength_pa())?;
            limit(
                "original adhesive static capacity",
                static_capacities_n[i],
                config.maximum_contact_force_n,
            )?;
            surface_energy_j = finite(surface_energy_j + contact.reversible_energy_j)?;
            contacts.push(contact);
        }
        let (stiffness, damping, links) = network_parameters(&network)?;
        let flow = MultipleLinearNetworkEventFlow::new(
            stiffness,
            damping,
            links,
            junctions.iter().map(|site| site.column.clone()).collect(),
            config.events,
        )?;
        let workspace = NativeLinearNetworkWorkspace::new(&network);
        limit(
            "initial total adhesive network energy magnitude",
            finite(network.total_energy_j()? + surface_energy_j)?.abs(),
            network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(Self {
            network,
            junctions,
            contacts,
            surface_energy_j,
            static_capacities_n,
            kinetic_capacities_n,
            phases: [LinearNetworkPhase::Fresh; MAX_LINEAR_NETWORK_CONTACTS],
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
    pub fn junctions(&self) -> &[AdhesiveNetworkEventJunction<'a>] {
        &self.junctions
    }
    pub fn normal_responses(&self) -> &[DmtResponse] {
        &self.contacts
    }
    pub fn phases(&self) -> &[LinearNetworkPhase] {
        &self.phases[..self.contact_count()]
    }
    pub const fn config(&self) -> MultipleAdhesiveNetworkEventModalConfig {
        self.config
    }
    pub const fn last_frame(&self) -> Option<&MultipleAdhesiveNetworkEventModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, Error> {
        Ok(finite(
            self.network.total_energy_j()? + self.surface_energy_j,
        )?)
    }
    pub fn step(
        &mut self,
        external: &[f64],
    ) -> Result<&MultipleAdhesiveNetworkEventModalFrame, Error> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&MultipleAdhesiveNetworkEventModalFrame, Error> {
        poll(Some(gate))?;
        if external.len() != self.mode_count() || external.iter().any(|f| !f.is_finite()) {
            return Err(invalid("external event forces must cover every original mode").into());
        }
        self.workspace.capture(&self.network);
        let count = self.contact_count();
        let path = self.flow.advance(
            self.workspace.positions(),
            self.workspace.velocities(),
            external,
            &self.phases[..count],
            self.network.dt,
            &self.static_capacities_n[..count],
            &self.kinetic_capacities_n[..count],
            gate,
        )?;
        for slip in &path.accumulated_slip_m[..count] {
            limit(
                "individual junction accumulated physical slip",
                *slip,
                self.config.maximum_accumulated_slip_m,
            )?;
        }
        let native = self.workspace.stage(
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
            gate,
        )?;
        let mut tangent_positions_m = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut tangent_velocities_m_s = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut total_friction = 0.0;
        for (i, junction) in self.junctions.iter().enumerate() {
            poll(Some(gate))?;
            limit(
                "individual junction speed over complete trajectory",
                finite(
                    det::sqrt(dot(&junction.column, &junction.column)?) * native.root_energy_bound,
                )?,
                self.config.maximum_tangent_speed_m_s,
            )?;
            tangent_positions_m[i] =
                projected_position(&self.network.candidates, &junction.column)?;
            tangent_velocities_m_s[i] =
                projected_velocity(&self.network.candidates, &junction.column)?;
            let loss = finite(path.contact_friction_dissipation_j[i])?;
            if loss < 0.0 {
                return Err(
                    invalid("individual junction friction loss must be nonnegative").into(),
                );
            }
            total_friction = finite(total_friction + loss)?;
        }
        let friction_residual = finite(total_friction - path.friction_dissipation_j)?;
        if finite(friction_residual.abs() + path.numerical_work_uncertainty_j)?
            > native.energy_tolerance_j
        {
            return Err(Error::EnergyAccuracy {
                residual_j: friction_residual,
                estimated_error_j: path.numerical_work_uncertainty_j,
                tolerance_j: native.energy_tolerance_j,
            });
        }
        let accepted = MultipleAdhesiveNetworkEventModalFrame {
            sample: native.sample,
            contact_count: count,
            observer_pressure_pa: native.observer_pressure_pa,
            network_energy_j: native.network_energy_j,
            surface_energy_j: self.surface_energy_j,
            external_work_j: native.external_work_j,
            friction_dissipation_j: path.friction_dissipation_j,
            modal_damping_dissipation_j: path.modal_damping_dissipation_j,
            connection_damping_dissipation_j: native.connection_damping_dissipation_j,
            energy_residual_j: native.energy_residual_j,
            energy_tolerance_j: native.energy_tolerance_j,
            numerical_work_uncertainty_j: path.numerical_work_uncertainty_j,
            tangent_positions_m,
            tangent_velocities_m_s,
            static_capacities_n: self.static_capacities_n,
            kinetic_capacities_n: self.kinetic_capacities_n,
            path,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        self.phases = path.phases;
        Ok(self.last.insert(accepted))
    }
}
