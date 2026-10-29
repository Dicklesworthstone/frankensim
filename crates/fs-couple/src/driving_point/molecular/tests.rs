//! G1/G3: one constitutive law reaches analytic waves and actual reed feedback.
use super::*;
use crate::render::ReedBoreVoice;
use crate::thin_plate::PlateBank;
use fs_duct::{LossModel, Segment};
use fs_material::gas::GasSpec;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_scenario::BeatingReed;
use fs_thermochem::SpeciesId;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
// Synthetic internal reservoirs test material-to-consumer causality. They do
// not assert measured molecular constants for the dry-air numerical baseline.
fn state(count: usize, rate_scale: f64) -> RelaxationState {
    let id = SpeciesId::new("synthetic").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: (0..count)
            .map(|j| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * j as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (2e-5 + j as f64 * 3e-4) * rate_scale,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas(), vec![(id, 1.0)]).unwrap())
    .unwrap()
}
fn geometry() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.15,
        }],
    }
}
fn point(count: usize, scale: f64, wall: bool) -> MolecularDrivingPoint {
    MolecularDrivingPoint::new(
        &geometry(),
        state(count, scale),
        wall.then_some(MolecularWallResolution {
            poles_per_field: 16,
            band_rad_s: [1e-4, 30000.0],
            maximum_relative_error: 0.05,
        }),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        gas(),
    )
    .unwrap()
}
fn zc(point: &MolecularDrivingPoint) -> f64 {
    let gas = point.inlet_gas();
    gas.density * gas.sound_speed / (core::f64::consts::PI * point.inlet_radius_m().powi(2))
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn zero_reservoirs_recover_existing_scattering_including_lossless_poles() {
    let model = point(0, 1.0, false);
    let radius = model.inlet_radius_m();
    for termination in [
        Termination::IdealOpen,
        Termination::Closed,
        Termination::UnflangedOpen,
    ] {
        for omega in [100.0, 1000.0, 6000.0] {
            let expected = fs_duct::compose::input_reflection_wall(
                &geometry(),
                &gas(),
                omega,
                LossModel::Lossless,
                termination,
                None,
            )
            .unwrap()
            .referenced_to(zc(&model))
            .unwrap();
            near(
                model.reflectance(omega, termination, zc(&model)).unwrap(),
                expected,
                2e-12,
            );
        }
    }
    let omega = core::f64::consts::PI * gas().sound_speed / (2.0 * 0.15);
    let pole = model
        .reflectance(omega, Termination::IdealOpen, zc(&model))
        .unwrap();
    near(pole, C64::ONE, 2e-12);
    assert_eq!(radius, 0.001);
}

#[test]
fn source_relaxation_changes_the_actual_port_and_keeps_wavefront_causality() {
    let slow = point(2, 1.0, false);
    let fast = point(2, 0.1, false);
    assert_eq!(slow.inlet_gas(), fast.inlet_gas());
    assert_eq!(slow.one_way_wavefront_s(), fast.one_way_wavefront_s());
    assert!(slow.one_way_wavefront_s() < 0.15 / gas().sound_speed);
    let w = 5000.0;
    let phase = slow
        .material()
        .state()
        .relaxing_wavenumber(w)
        .unwrap()
        .scale(0.15);
    let attenuation = (-2.0 * phase.im).exp();
    let travelling = C64::new((2.0 * phase.re).cos(), (2.0 * phase.re).sin()).scale(attenuation);
    let k = slow.material().state().relaxing_wavenumber(w).unwrap();
    let wave_z = C64::from_re(w * gas().density / (core::f64::consts::PI * 0.001_f64.powi(2))) / k;
    let input_z = wave_z * (C64::ONE - travelling) / (C64::ONE + travelling);
    let expected = fs_vfit::discretize::reflectance(input_z, zc(&slow));
    let actual = slow
        .reflectance(w, Termination::IdealOpen, zc(&slow))
        .unwrap();
    near(actual, expected, 2e-12);
    let changed = fast
        .reflectance(w, Termination::IdealOpen, zc(&fast))
        .unwrap();
    assert!((actual - changed).abs() > 1e-3);
}

fn reed() -> BeatingReed {
    BeatingReed {
        rest_opening_m: 4e-4,
        width_m: 0.003,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa: 2800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.0,
    }
}
fn voice(scale: f64) -> ReedBoreVoice {
    ReedBoreVoice::new_with_molecular_driving_point(
        point(2, scale, true),
        Termination::UnflangedOpen,
        reed(),
        PlateBank::default(),
        1.0,
        8000,
        &LineLadder::default(),
    )
    .unwrap()
}

#[test]
fn changed_reservoir_data_reaches_reed_pressure_and_block_replay() {
    let mut slow = voice(1.0);
    let mut fast = voice(0.1);
    let mut replay = voice(1.0);
    assert_eq!(
        slow.molecular_driving_point().unwrap().material().state(),
        &state(2, 1.0)
    );
    let (mut a, mut b, mut c) = (vec![0.0; 1024], vec![0.0; 1024], vec![0.0; 1024]);
    slow.step_block(&mut a).unwrap();
    fast.step_block(&mut b).unwrap();
    for block in c.chunks_mut(37) {
        replay.step_block(block).unwrap();
    }
    assert_eq!(a, c);
    assert_eq!(slow.samples_rendered(), replay.samples_rendered());
    let peak = a.iter().map(|p| p.abs()).fold(0.0, f64::max);
    let difference = a
        .iter()
        .zip(&b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    assert!(
        peak > 1e-3 && difference > 1e-5 * peak,
        "peak={peak}, difference={difference}"
    );
    // This is internal pressure too: the difference precedes any observer law.
    assert_ne!(
        slow.last_observer_parts().internal_pa,
        fast.last_observer_parts().internal_pa
    );
}

#[test]
fn molecular_material_band_refuses_without_an_equilibrium_fallback() {
    let model = MolecularDrivingPoint::new(
        &geometry(),
        state(2, 1.0),
        Some(MolecularWallResolution {
            poles_per_field: 16,
            band_rad_s: [1.0, 30000.0],
            maximum_relative_error: 0.05,
        }),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        gas(),
    )
    .unwrap();
    assert!(
        model
            .realize_line(
                Termination::UnflangedOpen,
                8000,
                zc(&model),
                &LineLadder::default()
            )
            .is_err()
    );
    let mut stale = gas();
    stale.pressure *= 2.0;
    assert!(
        MolecularDrivingPoint::new(
            &geometry(),
            state(2, 1.0),
            None,
            ConeSlices::Fixed(8),
            ComposeOptions::default(),
            stale
        )
        .is_err()
    );
    let unequal_pressure =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 2.0 * gas().pressure).unwrap();
    assert!(
        MolecularDrivingPoint::new(
            &geometry(),
            state(2, 1.0),
            None,
            ConeSlices::Fixed(8),
            ComposeOptions::default(),
            unequal_pressure,
        )
        .is_err()
    );
}
