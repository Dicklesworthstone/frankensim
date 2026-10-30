//! G0/G3: borrowed material owners retain exact admission, not raw-data fallbacks.
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState,
    ThermalFluidTransport, ThermalFluidTransportState};
use fs_material::gas::{GasSpec, GasState};

fn transport<F: FluidTransport>(source: F) -> Result<FluidTransportState, MaterialError> {
    source.transport_state()
}
fn thermal<F: ThermalFluidTransport>(source: F) -> Result<ThermalFluidTransportState, MaterialError> {
    source.thermal_transport_state()
}

#[test]
fn heterogeneous_borrows_reuse_original_transport_without_a_material_copy() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, 120_000.0).unwrap();
    let explicit = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let sources: [&dyn FluidTransport; 2] = [&gas, &explicit];
    assert_eq!(transport(&gas).unwrap(), gas.transport_state().unwrap());
    assert_eq!(transport(&explicit).unwrap(), explicit);
    assert_eq!(transport(sources[0]).unwrap(), gas.transport_state().unwrap());
    assert_eq!(transport(sources[1]).unwrap(), explicit);
    assert_eq!(transport(&sources[0]).unwrap(), transport(sources[0]).unwrap());
}

#[test]
fn borrowed_thermal_sources_preserve_unknown_expansion_and_the_same_temperature() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 321.0, 100_000.0).unwrap();
    let liquid = ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap(), 301.0, None).unwrap();
    let sources: [&dyn ThermalFluidTransport; 2] = [&gas, &liquid];
    let a = thermal(sources[0]).unwrap();
    let b = thermal(sources[1]).unwrap();
    assert_eq!(a, gas.thermal_transport_state().unwrap());
    assert_eq!(b, liquid);
    assert_eq!(a.transport(), transport(sources[0]).unwrap());
    assert_eq!(b.transport(), transport(sources[1]).unwrap());
    assert_eq!(a.temperature_k(), 321.0);
    assert_eq!(b.temperature_k(), 301.0);
    assert_eq!(b.volumetric_expansion_per_k(), None);
}

#[test]
fn borrowed_source_errors_are_not_reinterpreted_as_authored_transport() {
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    gas.prandtl = 100.0;
    let momentum: &dyn FluidTransport = &gas;
    let heat: &dyn ThermalFluidTransport = &gas;
    assert_eq!(transport(momentum).unwrap_err().to_string(), gas.transport_state().unwrap_err().to_string());
    assert_eq!(thermal(heat).unwrap_err().to_string(), gas.thermal_transport_state().unwrap_err().to_string());
    assert!(transport(&gas).is_err());
    assert!(thermal(&gas).is_err());
}
