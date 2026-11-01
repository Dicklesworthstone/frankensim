//! Physical modal admittance source, before any discretization in time.

use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::ModalTimeError;
use fs_duct::{LossModel, segment_wave};
use fs_material::gas::GasState;
use fs_phs::{PortHamiltonian, QuadraticStorage};

/// Original mouth basis and the load's own coherent gas. Mode order and
/// normalization are exactly those of `fs_duct::modal`'s axisymmetric basis:
/// pressure is conjugate to modal volume flow, with plane mode first.
///
/// The producing physical source must declare its actual coordinates here.
/// This declaration does not turn an arbitrary transfer fit into a physical
/// storage law or establish the source's spatial or frequency accuracy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalLoadSource {
    radius_m: f64,
    modes: ModeCount,
    gas: GasState,
}
impl ModalLoadSource {
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    pub const fn modes(&self) -> ModeCount {
        self.modes
    }
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
}

/// Complete original pressure-driven load. Its N inputs are mouth modal
/// pressures [Pa], and its N outputs are volume flows INTO the load [m³/s].
/// There are no additional unbound source ports. The original system may
/// contain arbitrary coupling between modes and linear material-memory states.
///
/// The runtime composes this actual system before exact propagation. It never
/// converts interval-average feedback into an instantaneous boundary law and
/// never discards load states. Source units/basis remain the producer's contract;
/// the constructor checks their declared geometry, gas and port dimensions.
pub struct ModalAdmittanceLoad {
    source: ModalLoadSource,
    system: Option<PortHamiltonian>,
    direct: Vec<f64>,
}
impl ModalAdmittanceLoad {
    /// Bind the original source system to its actual circular modal mouth.
    /// Linear/quadratic storage and complete-state capacity are admitted when
    /// the load and duct are prepared together; no reduction is performed.
    pub fn new(
        system: PortHamiltonian,
        radius_m: f64,
        modes: ModeCount,
        gas: GasState,
    ) -> Result<Self, ModalTimeError> {
        modes.validate()?;
        segment_wave(&gas, radius_m, 1.0, LossModel::Lossless)?;
        if system.port_dim() != modes.0 || system.state_dim() == 0 {
            return Err(ModalTimeError::Invalid(
                "modal load requires one pressure-input/flow-output port per retained mouth mode",
            ));
        }
        Ok(Self {
            source: ModalLoadSource {
                radius_m,
                modes,
                gas,
            },
            system: Some(system),
            direct: vec![0.0; modes.0 * modes.0],
        })
    }
    /// Original memoryless relation U=D*P, with symmetric positive semidefinite
    /// D [m³/(Pa s)]. All self and mutual entries are retained. Zero D is the
    /// exact rigid boundary; no artificial memory coordinate is introduced.
    pub fn from_conductance(
        conductance: Vec<f64>,
        radius_m: f64,
        modes: ModeCount,
        gas: GasState,
    ) -> Result<Self, ModalTimeError> {
        modes.validate()?;
        segment_wave(&gas, radius_m, 1.0, LossModel::Lossless)?;
        let admitted = QuadraticStorage::new(conductance, modes.0)?;
        Ok(Self {
            source: ModalLoadSource {
                radius_m,
                modes,
                gas,
            },
            system: None,
            direct: admitted.matrix().to_vec(),
        })
    }
    /// Set the actual direct conductance of a dynamic source:
    /// U=Gᵀ grad(H_load)+D*P. This is physical loss, independently integrated
    /// with the original load loss; it is not reference-wave loading.
    pub fn with_direct_conductance(
        mut self,
        conductance: Vec<f64>,
    ) -> Result<Self, ModalTimeError> {
        let admitted = QuadraticStorage::new(conductance, self.source.modes.0)?;
        self.direct = admitted.matrix().to_vec();
        Ok(self)
    }
    /// Original physical port declaration, retained after attachment.
    pub const fn source(&self) -> &ModalLoadSource {
        &self.source
    }
    /// Complete original source system, without a fitted replacement.
    pub const fn system(&self) -> Option<&PortHamiltonian> {
        self.system.as_ref()
    }
    /// Original row-major physical direct conductance in declared mode order.
    pub fn direct_conductance(&self) -> &[f64] {
        &self.direct
    }
    pub(super) fn into_parts(self) -> (ModalLoadSource, Option<PortHamiltonian>, Vec<f64>) {
        (self.source, self.system, self.direct)
    }
}
