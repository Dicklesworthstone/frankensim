//! Public material-wall cavity physics through the ordinary integration seam.
use fs_couple::cavity_modal::shared::{
    CavitySurface, CylindricalMaterialWallCompression, CylindricalMaterialWallOptions, SharedCavity,
};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, ModalCouplingError};
use fs_duct::LossModel;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_material::relaxation::RelaxationState;
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSample};
use fs_math::{c64::C64, det};
use fs_scenario::HelmholtzCavity;

#[path = "../src/cavity_modal/shared/thermal/material_wall/tests.rs"]
mod tests;
