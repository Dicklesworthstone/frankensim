//! Existing material/contact apertures at explicit fixed graph endpoints.
//!
//! Upstream is the exterior reservoir and downstream is the chimney. Positive
//! chimney pressure opens the boundary; an applied negative opening force closes
//! it. The original aperture supplies both jet flow and swept volume, and every
//! endpoint reads old waves before any state is accepted. Closing a gap therefore
//! retains chimney pressure, all propagation waves, contact and material
//! histories. This is a moving slit model, not a resolved circular key-pad domain.
use super::{TubeNetworkSpec, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::{
    ApertureDrive, ApertureFrame, ApertureTrial, DynamicAperture,
};
use fs_material::gas::GasState;
use fs_vfit::waveguide::network::{DrivenWave, NetworkFrame, WaveguideNetwork};

/// Original mechanics to attach before any participant is advanced.
pub struct BoundaryAperture {
    /// Original degree-one driven graph node.
    pub node: usize,
    /// Existing source-backed or explicitly authored aperture owner.
    pub aperture: DynamicAperture,
}

struct BoundAperture {
    node: usize,
    gas: GasState,
    aperture: DynamicAperture,
}

/// Independent reservoir pressure and work-conjugate mechanical force.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryDrive {
    /// Original node address, in the same order as the admitted boundaries.
    pub node: usize,
    /// Reservoir acoustic pressure relative to the same equilibrium [Pa].
    pub exterior_pressure_pa: f64,
    /// Force increasing the opening coordinate [N]; a closing key uses negative.
    pub opening_force_n: f64,
}

/// One shared accepted step including every moving boundary and duct load.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BoundaryNetworkFrame {
    /// Original wave/load result, including separate MIMO terminal exchange.
    pub network: NetworkFrame,
    /// Wave, local duct load, and aperture mechanical/contact/material storage [J].
    pub stored_energy_j: f64,
    /// Change of independently computed complete storage [J].
    pub storage_change_j: f64,
    /// Aperture dissipation and local terminal/interior loss [J].
    pub dissipated_energy_j: f64,
    /// Sum of exterior reservoir pressure times actual boundary flow times dt [J].
    pub exterior_work_j: f64,
    /// Sum of physical boundary actuation work [J].
    pub mechanical_work_j: f64,
}
impl BoundaryNetworkFrame {
    /// Complete balance; the external rational radiation operator has its own
    /// storage owner, so its signed work remains an outgoing exchange here.
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j + self.network.coupled_terminal_work_j
            - self.network.inlet_work_j
            - self.exterior_work_j
            - self.mechanical_work_j
    }
}

/// Mechanical boundary owners and their fixed physical graph addresses.
/// The caller retains the network, but any independent advance is detected
/// before the next coupled step. Trials allocate through the existing aperture.
pub struct BoundaryApertures {
    boundaries: Vec<BoundAperture>,
    last_apertures: Vec<ApertureFrame>,
    last: Option<BoundaryNetworkFrame>,
    accepted_steps: u64,
    time_step_s: f64,
}

/// Internal complete candidate for a larger physical source/network owner.
/// Only the original owner can publish it; no accepted state is changed by its
/// construction. The primary inlet reed can add its work before acceptance.
pub(super) struct BoundaryTrial {
    pub frame: BoundaryNetworkFrame,
    accepted_steps: u64,
    inlet_outgoing_pa: f64,
    waves: Vec<DrivenWave>,
    apertures: Vec<ApertureTrial>,
}
fn wave(error: fs_vfit::waveguide::WaveguideError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}

impl BoundaryApertures {
    /// Bind every declared driven endpoint once, in ascending node order, to
    /// its actual local gas, impedance and clock. No inlet geometry is reused
    /// for another hole. A moving boundary requires a represented chimney.
    pub fn new(
        graph: &TubeNetworkSpec,
        gases: &[GasState],
        network: &WaveguideNetwork,
        inputs: Vec<BoundaryAperture>,
        max_apertures: usize,
    ) -> Result<Self, AcousticRealizeError> {
        if !(1..=64).contains(&max_apertures)
            || inputs.is_empty()
            || inputs.len() > max_apertures
            || network.accepted_steps() != 0
        {
            return Err(invalid(
                "moving boundaries require a cold graph and 1..=64 explicitly budgeted apertures",
            ));
        }
        graph.validate_section_gases(gases)?;
        if !network
            .driven_nodes()
            .eq(inputs.iter().map(|input| input.node))
        {
            return Err(invalid(
                "every driven graph endpoint needs one aperture in ascending node order",
            ));
        }
        let mut boundaries = Vec::with_capacity(inputs.len());
        for input in inputs {
            let boundary = network.driven_boundary(input.node).map_err(wave)?;
            let sections: Vec<_> = graph
                .sections
                .iter()
                .enumerate()
                .filter(|(_, section)| section.nodes.contains(&input.node))
                .collect();
            if sections.len() != 1
                || sections[0].0 != boundary.section
                || graph.nodes.get(input.node) != network.node(input.node).as_ref()
            {
                return Err(invalid(
                    "moving boundary must preserve its original degree-one graph address",
                ));
            }
            let gas = gases[sections[0].0];
            // Re-admit the original local geometry through its actual integer
            // or fractional owner, including travel time as well as impedance.
            // The zero-frequency phase query adds no new bandwidth policy.
            let realized = if let Some(transit) = boundary.fractional_transit {
                let (realized, detail) =
                    sections[0]
                        .1
                        .realize_fractional(&gas, network.time_step_s(), 0.0)?;
                if transit != detail.transit {
                    return Err(invalid(
                        "moving boundary must retain its original fractional chimney transit",
                    ));
                }
                realized
            } else {
                sections[0]
                    .1
                    .realize(gas.sound_speed, gas.density, network.time_step_s())?
            };
            let spec = input.aperture.spec();
            if input.aperture.accepted_steps() != 0
                || spec.time_step_s.to_bits() != network.time_step_s().to_bits()
                || realized.one_way_samples != boundary.one_way_samples
                || realized.impedance_pa_s_m3.to_bits() != boundary.impedance_pa_s_m3.to_bits()
                || spec.impedance_pa_s_m3.to_bits() != boundary.impedance_pa_s_m3.to_bits()
                || spec.density_kg_m3.to_bits() != gas.density.to_bits()
                || input.aperture.viscous_slit().is_some_and(|law|
                    law.gas().is_some_and(|source| *source != gas))
            {
                return Err(invalid(
                    "moving boundary clock, state, gas and characteristic load must match its own chimney",
                ));
            }
            boundaries.push(BoundAperture {
                node: input.node,
                gas,
                aperture: input.aperture,
            });
        }
        Ok(Self {
            last_apertures: vec![ApertureFrame::default(); boundaries.len()],
            boundaries,
            last: None,
            accepted_steps: 0,
            time_step_s: network.time_step_s(),
        })
    }

    /// Physical nodes in original graph order.
    pub fn nodes(&self) -> impl Iterator<Item = usize> + '_ {
        self.boundaries.iter().map(|b| b.node)
    }
    /// Original mechanical/material/contact owner and its accepted history.
    pub fn aperture(&self, node: usize) -> Option<&DynamicAperture> {
        self.boundaries
            .iter()
            .find(|b| b.node == node)
            .map(|b| &b.aperture)
    }
    /// Complete local source gas used to admit this aperture's flow law.
    pub fn gas(&self, node: usize) -> Option<&GasState> {
        self.boundaries
            .iter()
            .find(|b| b.node == node)
            .map(|b| &b.gas)
    }
    /// Last accepted boundary pressure, flow, motion and material loss.
    pub fn aperture_frame(&self, node: usize) -> Option<&ApertureFrame> {
        self.last.as_ref()?;
        self.boundaries
            .iter()
            .position(|b| b.node == node)
            .map(|i| &self.last_apertures[i])
    }
    /// Complete last accepted work/storage balance.
    pub const fn last_frame(&self) -> Option<&BoundaryNetworkFrame> {
        self.last.as_ref()
    }
    /// Boundary-only storage [J]; the caller adds the network's actual storage.
    pub fn stored_energy_j(&self) -> f64 {
        self.boundaries
            .iter()
            .map(|b| b.aperture.stored_energy_j())
            .sum()
    }

    pub(super) fn validate_step_budget(&self, total: u64) -> Result<(), AcousticRealizeError> {
        for boundary in &self.boundaries {
            boundary.aperture.validate_step_budget(total)?;
        }
        Ok(())
    }

    // The enclosing source validates this complete set before changing its own
    // horizon. Accepted physical state and material history remain untouched.
    pub(super) fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        self.validate_step_budget(total)?;
        for boundary in &mut self.boundaries {
            boundary.aperture.extend_step_budget(total)?;
        }
        Ok(())
    }

    /// Solve all apertures against old waves, preview the shared duct graph and
    /// verify finite combined work before publishing any participant. Input or
    /// solver refusal preserves every accepted physical coordinate.
    pub fn step(
        &mut self,
        network: &mut WaveguideNetwork,
        inlet_outgoing_pa: f64,
        drives: &[BoundaryDrive],
    ) -> Result<BoundaryNetworkFrame, AcousticRealizeError> {
        let trial = self.preview_step(network, inlet_outgoing_pa, drives)?;
        self.accept_trial(network, trial)
    }

    pub(super) fn preview_step(
        &self,
        network: &mut WaveguideNetwork,
        inlet_outgoing_pa: f64,
        drives: &[BoundaryDrive],
    ) -> Result<BoundaryTrial, AcousticRealizeError> {
        if drives.len() != self.boundaries.len()
            || network.accepted_steps() != self.accepted_steps
            || network.time_step_s().to_bits() != self.time_step_s.to_bits()
            || !network.driven_nodes().eq(self.nodes())
        {
            return Err(invalid(
                "moving boundaries must retain the original shared graph clock and complete assignment",
            ));
        }
        let mut trials = Vec::with_capacity(drives.len());
        let mut waves = Vec::with_capacity(drives.len());
        let mut frame = BoundaryNetworkFrame::default();
        for (boundary, drive) in self.boundaries.iter().zip(drives) {
            let port = network.driven_boundary(boundary.node).map_err(wave)?;
            if drive.node != boundary.node
                || !drive.exterior_pressure_pa.is_finite()
                || !drive.opening_force_n.is_finite()
                || boundary.aperture.accepted_steps() != self.accepted_steps
                || boundary.aperture.spec().impedance_pa_s_m3.to_bits()
                    != port.impedance_pa_s_m3.to_bits()
            {
                return Err(invalid(
                    "moving boundary requires its original local impedance and finite ordered drive",
                ));
            }
            let trial = boundary.aperture.preview_step_with_force(
                ApertureDrive {
                    upstream_pressure_pa: drive.exterior_pressure_pa,
                    incoming_pressure_pa: port.incoming_pressure_pa,
                    body_flow_m3_s: 0.0,
                },
                drive.opening_force_n,
            )?;
            let aperture = trial.frame;
            frame.stored_energy_j += aperture.stored_energy_j;
            frame.storage_change_j += aperture.storage_change_j;
            frame.dissipated_energy_j += aperture.dissipated_energy_j;
            frame.exterior_work_j +=
                drive.exterior_pressure_pa * (aperture.bore_flow_m3_s * self.time_step_s);
            frame.mechanical_work_j += aperture.mechanical_work_j;
            waves.push(DrivenWave {
                node: boundary.node,
                outgoing_pressure_pa: aperture.outgoing_pressure_pa,
            });
            trials.push(trial);
        }
        frame.network = network
            .preview_step_with_driven(inlet_outgoing_pa, &waves)
            .map_err(wave)?;
        frame.stored_energy_j += frame.network.stored_energy_j;
        frame.storage_change_j += frame.network.storage_change_j;
        frame.dissipated_energy_j += frame.network.terminal_loss_j + frame.network.interior_loss_j;
        if ![
            frame.stored_energy_j,
            frame.storage_change_j,
            frame.dissipated_energy_j,
            frame.exterior_work_j,
            frame.mechanical_work_j,
            frame.balance_residual_j(),
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(invalid(
                "moving-boundary network storage or work left the finite set",
            ));
        }
        Ok(BoundaryTrial {
            frame,
            accepted_steps: self.accepted_steps,
            inlet_outgoing_pa,
            waves,
            apertures: trials,
        })
    }

    // The owning composite must not mutate geometry, load controls or mechanics
    // between this trial's preview and acceptance. Existing aperture acceptance
    // is infallible after the original graph has validated and committed.
    pub(super) fn accept_trial(
        &mut self,
        network: &mut WaveguideNetwork,
        trial: BoundaryTrial,
    ) -> Result<BoundaryNetworkFrame, AcousticRealizeError> {
        if self.accepted_steps != trial.accepted_steps
            || network.accepted_steps() != trial.accepted_steps
        {
            return Err(invalid(
                "moving boundary trial must retain its original shared accepted-step clock",
            ));
        }
        network
            .step_with_driven(trial.inlet_outgoing_pa, &trial.waves)
            .map_err(wave)?;
        for ((boundary, last), trial) in self
            .boundaries
            .iter_mut()
            .zip(&mut self.last_apertures)
            .zip(trial.apertures)
        {
            *last = trial.frame;
            boundary.aperture.accept_frame(trial);
        }
        self.accepted_steps = network.accepted_steps();
        self.last = Some(trial.frame);
        Ok(trial.frame)
    }
}
