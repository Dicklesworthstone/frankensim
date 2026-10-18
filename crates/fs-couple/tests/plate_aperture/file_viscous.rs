//! G3: source geometry and inlet material reach the same-sample valve solve.
use super::*;
use fs_couple::bernoulli_aperture::performance::CoupledAperture;
use fs_couple::bernoulli_aperture::viscous::ViscousSlit;

fn with_slit(text: &str, length: &str) -> String {
    text.replace("\ninitial ", &format!("\nslit_length {length}\ninitial "))
}

#[test]
fn file_channel_matches_direct_spatial_flow_and_changes_retained_mechanics() {
    let source = load(&with_slit(INPUT, "0.001"));
    let seed = load(INPUT);
    let CoupledAperture::Tube(tube) = seed.renderer().system() else {
        panic!("fixture must retain its original tube");
    };
    let a = tube.aperture();
    let spec = a.spec();
    let mut valve = DynamicAperture::from_plate_with_closure_and_propagator(
        a.plate_reduction().unwrap().clone(),
        a.plate_closure().unwrap().spec().clone(),
        spec.density_kg_m3,
        tube.spec().characteristic_impedance(spec.density_kg_m3).unwrap(),
        spec.time_step_s,
        spec.max_steps,
        a.state(),
        a.blade_propagator(),
    ).unwrap().with_plate_relaxation(
        a.relaxation().unwrap().spec().clone(), InitialApertureMemory::Relaxed,
    ).unwrap();
    valve.set_viscous_slit(ViscousSlit::from_gas(seed.ambient_gas(), 0.001).unwrap()).unwrap();
    let mut direct = ApertureTube::new(valve, *tube.spec()).unwrap();
    let schedule = seed.renderer().schedule().clone();
    let mut baseline = seed.into_renderer();
    let mut actual = source.into_renderer();
    let mut changed = false;
    for n in 0..2048 {
        let frame = direct.step(TubeDrive {
            upstream_pressure_pa: schedule.sample("mouth", n * 700 / 96000).unwrap(),
            body_flow_m3_s: 0.0,
        }).unwrap();
        let mut out = [0.0];
        actual.block(&mut out).unwrap();
        baseline.block(&mut [0.0]).unwrap();
        assert_eq!(out[0].to_bits(), frame.aperture.bore_pressure_pa.to_bits());
        changed |= actual.system().aperture().state() != baseline.system().aperture().state();
    }
    assert!(changed, "channel loss must change valve motion, not only reported material data");
    assert_eq!(actual.system().aperture().state(), direct.aperture().state());
}

#[test]
fn regional_inlet_composition_supplies_both_channel_viscosity_and_density() {
    const REGIONAL: &str = include_str!("../../examples/plate-valve-regional-gas.performance");
    let input = with_slit(REGIONAL, "0.001").replace(
        "duct_section 0 1 0.125 0.007 0.002 viscothermal",
        "duct_section 0 1 0.125 0.007 0.002 gas 303.15 0.2 viscothermal",
    );
    let source = load(&input);
    let inlet = GasState::try_new_moist_air(303.15, 101325.0, 0.2).unwrap();
    let valve = source.renderer().system().aperture();
    let law = valve.viscous_slit().unwrap();
    assert_eq!(law.gas(), Some(&inlet));
    assert_eq!(law.dynamic_viscosity_pa_s(), inlet.dynamic_viscosity);
    assert_eq!(valve.spec().density_kg_m3, inlet.density);
    assert_ne!(law.dynamic_viscosity_pa_s(), source.ambient_gas().dynamic_viscosity);
    // A different surrounding gas changes radiation, not the channel medium.
    let exterior = input.replace("ambient 293.15 101325 0", "ambient 293.15 101325 0\nexterior 283.15 101325 0");
    let outside = load(&exterior);
    assert_eq!(outside.renderer().system().aperture().viscous_slit().unwrap().gas(), Some(&inlet));
}

#[test]
fn invalid_or_duplicate_channel_geometry_refuses_instead_of_falling_back() {
    for length in ["0", "-0.001", "NaN", "inf", "0.001 extra", "0.001\nslit_length 0.002"] {
        let text = with_slit(INPUT, length);
        assert!(PlateValvePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new()).is_err(),
            "accepted {length:?}");
    }
    assert!(load(INPUT).renderer().system().aperture().viscous_slit().is_none());
    let symmetric = with_slit(INPUT, "0.001").replace(
        "aperture 0.0002 0.25 0.1", "aperture 0.0002 0.25 0.1 symmetric-opposed-pair",
    );
    let pair = load(&symmetric);
    assert_eq!(pair.info().moving_blades, 2);
    assert_eq!(pair.renderer().system().aperture().viscous_slit().unwrap().length_m(), 0.001);
}
