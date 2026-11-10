//! fs-nlmodal's `ModalStructure` refuses a modal model whose modes carry an observer transfer: the
//! Component books its modes' work and heat only, so a model that also radiates to an observer
//! would lose that energy from the ledger. The same model without its transfers is admitted.
//! Run with `cargo test -p fs-nlmodal --test modal_structure`.

use fs_math::c64::C64;
use fs_nlmodal::component::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalStructure,
    ModalStructureRefusal,
};
use fs_ports::component::Component;

#[test]
fn a_model_with_an_observer_transfer_is_refused() {
    let budget = ModalAcousticTimeBudget::audible_reference();
    let mode = |hz: f64, transfer: C64| ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * hz,
        damping_ratio: 1.0e-4,
        pressure_per_modal_velocity: transfer,
    };
    let model = ModalAcousticTimeModel::try_new(
        44_100,
        vec![mode(300.0, C64::ZERO), mode(900.0, C64::new(0.0, 1.0e-3))],
        budget,
    )
    .expect("model");
    assert_eq!(
        ModalStructure::new(model.clone(), Vec::new()).unwrap_err(),
        ModalStructureRefusal::ObserverTransfer { mode: 1 }
    );
    let silent =
        ModalStructure::new(model.without_pressure_transfers(), Vec::new()).expect("no transfer");
    assert_eq!((silent.modes(), silent.ports().len()), (2, 1));
}
