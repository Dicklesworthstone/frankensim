//! The retained field reaches the same pressure feedback and mouth observation.
use super::*;
use crate::driving_point::regional::field::tests::{cylinder, field};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_duct::HoleState;
use fs_material::relaxation::{
    InternalReservoir, RelaxationProcess, RelaxationRateLaw, RelaxingMedium, SpeciesRelaxation,
    VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn law(time_scale: f64) -> RelaxationAbsorption {
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: SpeciesId::new("analytic-field-gas").unwrap(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k: 1300.0,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 20.0 * time_scale,
                t_ref_k: 300.0,
                exponent: 1.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }],
    }])
    .unwrap()
}

fn point(geometry: &Duct, source: &AxialGasField, time_scale: f64) -> MolecularDrivingPoint {
    MolecularDrivingPoint::from_gas_field(
        geometry,
        source,
        &law(time_scale),
        None,
        RiccatiSteps(24),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        *source.state_at(source.length_m()).unwrap().state(),
    )
    .unwrap()
}

fn reference(model: &MolecularDrivingPoint) -> f64 {
    model.inlet_gas().characteristic_impedance
        / (core::f64::consts::PI * model.inlet_radius_m().powi(2))
}

#[test]
fn uniform_field_preserves_original_feedback_and_flow_bits() {
    let source = field(300.0);
    let geometry = cylinder();
    let material = law(1.0);
    let state = material
        .evaluate(&RelaxingMedium::from_thermochemical(&source.state_at(0.0).unwrap()).unwrap())
        .unwrap();
    let direct = MolecularDrivingPoint::new(
        &geometry,
        state,
        None,
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        *source.state_at(0.0).unwrap().state(),
    )
    .unwrap();
    let continuous = point(&geometry, &source, 1.0);
    assert_eq!(continuous.material().gas_field(), Some(&source));
    assert_eq!(
        continuous.one_way_wavefront_s(),
        direct.one_way_wavefront_s()
    );
    for omega in [100.0, 1000.0, 5000.0] {
        for termination in [
            Termination::IdealOpen,
            Termination::Closed,
            Termination::FlangedOpen,
        ] {
            assert_eq!(
                continuous
                    .reflectance(omega, termination, reference(&continuous))
                    .unwrap(),
                direct
                    .reflectance(omega, termination, reference(&direct))
                    .unwrap(),
            );
        }
        assert_eq!(
            continuous.flanged_outlet_flow(omega).unwrap(),
            direct.flanged_outlet_flow(omega).unwrap()
        );
    }
}

#[test]
fn continuous_molecular_field_reaches_feedback_and_each_causal_mouth() {
    let source = field(420.0);
    let radius = cylinder().segments[0].inlet_radius();
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius,
                length: 0.21,
            },
            Segment::ToneHole {
                hole_radius: 0.002,
                bore_radius: radius,
                chimney_height: 0.004,
                state: HoleState::Vent(0.4),
            },
            Segment::Cylinder {
                radius,
                length: source.length_m() - 0.21,
            },
        ],
    };
    let model = point(&geometry, &source, 1.0);
    assert_eq!(model.inlet_gas(), source.state_at(0.0).unwrap().state());
    assert_eq!(model.material().outlet_state().gas(), model.exterior_gas());
    let omega = 2800.0;
    let response = model
        .material()
        .response_with_holes(
            omega,
            model
                .terminal_load(omega, Termination::FlangedOpen)
                .unwrap(),
            &model.hole_mouth_loads(omega).unwrap(),
            ComposeOptions::default(),
        )
        .unwrap();
    let z0 = reference(&model);
    let expected =
        (response.duct.z_in - C64::from_re(z0)) / (response.duct.z_in + C64::from_re(z0));
    assert!(
        (model
            .reflectance(omega, Termination::FlangedOpen, z0)
            .unwrap()
            - expected)
            .abs()
            < 2e-12
    );
    let uniform = point(&geometry, &field(300.0), 1.0);
    assert!(
        (model
            .reflectance(omega, Termination::FlangedOpen, z0)
            .unwrap()
            - uniform
                .reflectance(omega, Termination::FlangedOpen, z0)
                .unwrap())
        .abs()
            > 1e-3
    );
    let faster = point(&geometry, &source, 0.1);
    assert_eq!(model.one_way_wavefront_s(), faster.one_way_wavefront_s());
    assert!(
        (model.flanged_outlet_flow(omega).unwrap() - faster.flanged_outlet_flow(omega).unwrap())
            .abs()
            > 1e-4 * model.flanged_outlet_flow(omega).unwrap().abs()
    );

    let rate = 8000;
    let receiver = CircularOutletReceiver {
        position_m: [0.0, 0.0, 1.0],
        radial_rings: 2,
        angular_points: 8,
        maximum_frequency_hz: 800.0,
    };
    let (terminal, mut holes) = model
        .baffled_outlet_and_open_holes_delayed_at(
            rate,
            512,
            receiver,
            model.exterior_gas().sound_speed,
        )
        .unwrap();
    assert_eq!(holes.len(), 1);
    assert_eq!(
        terminal.info().axial_delay_samples,
        (model.one_way_wavefront_s() * rate as f64).floor() as usize
    );
    assert_eq!(
        holes[0].info().axial_delay_samples,
        (model.open_hole_wavefront_times_s()[0] * rate as f64).floor() as usize
    );
    let mut observers = [terminal, holes.pop().unwrap()];
    let mut energy = [0.0; 2];
    for n in 0..1024 {
        for (observer, energy) in observers.iter_mut().zip(&mut energy) {
            let first =
                observer.info().axial_delay_samples + observer.info().exterior_delay_samples.0;
            let pressure = observer.step(if n == 0 { 1.0 } else { 0.0 }).unwrap();
            assert!(pressure.is_finite());
            if n < first {
                assert_eq!(pressure, 0.0);
            }
            *energy += pressure * pressure;
        }
    }
    assert!(energy.iter().all(|e| *e > 0.0));
}

#[test]
fn field_consumer_keeps_source_and_static_pressure_refusals() {
    let source = field(420.0);
    let mut exterior = *source.state_at(0.0).unwrap().state();
    exterior.pressure *= 1.1;
    assert!(
        MolecularDrivingPoint::from_gas_field(
            &cylinder(),
            &source,
            &law(1.0),
            None,
            RiccatiSteps(16),
            ConeSlices::Fixed(8),
            ComposeOptions::default(),
            exterior,
        )
        .is_err()
    );
    let unknown = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: SpeciesId::new("missing-source-gas").unwrap(),
        processes: vec![],
    }])
    .unwrap();
    assert!(
        MolecularDrivingPoint::from_gas_field(
            &cylinder(),
            &source,
            &unknown,
            None,
            RiccatiSteps(16),
            ConeSlices::Fixed(8),
            ComposeOptions::default(),
            *source.state_at(source.length_m()).unwrap().state(),
        )
        .is_err()
    );
}
