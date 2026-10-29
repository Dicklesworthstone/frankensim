use super::*;
use crate::dg::{DgBoundary, DgDuct, DgMesh, PortId};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
// Synthetic reservoir data isolate the constitutive equations, not empirical
// validation of a measured gas. p*tau is held when pressure changes.
fn law(count: usize, scale: f64, gas: GasState) -> RelaxationState {
    let id = SpeciesId::new("synthetic").unwrap();
    let data = SpeciesRelaxation {
        species: id.clone(),
        processes: (0..count)
            .map(|j| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * j as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (2e-5 + j as f64 * 3e-4) * scale,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    };
    RelaxationAbsorption::new(vec![data])
        .unwrap()
        .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
        .unwrap()
}
fn geometry(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}
fn wall() -> MolecularWallResolution {
    MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [0.1, 7000.0],
        maximum_relative_error: 0.02,
    }
}
fn prepare(
    geometry: &Duct,
    state: RelaxationState,
    wall: Option<MolecularWallResolution>,
) -> MolecularDuct {
    MolecularDuct::new(geometry, state, wall, ConeSlices::Fixed(16), 4096).unwrap()
}
fn sin_cos(z: C64) -> (C64, C64) {
    let (s, c) = det::sin_cos(z.re);
    let e = det::exp(z.im);
    let ch = (e + 1.0 / e) * 0.5;
    let sh = (e - 1.0 / e) * 0.5;
    (C64::new(s * ch, c * sh), C64::new(c * ch, -s * sh))
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}
fn local_wave(model: &MolecularDuct, omega: f64) -> SegmentWave {
    model
        .wave(
            0,
            omega,
            model.state.caloric_response(omega).unwrap().cp_j_kg_k,
            model.state.relaxing_wavenumber(omega).unwrap(),
        )
        .unwrap()
}

#[test]
fn molecular_cylinder_matches_material_and_closed_form_terminated_wave() {
    let state = law(2, 1.0, gas());
    let radius = 0.001;
    let length = 0.05;
    let model = prepare(&geometry(radius, length), state.clone(), None);
    let swapped = prepare(&geometry(radius, length), law(2, 0.1, gas()), None);
    let mut rate_effect = 0.0_f64;
    for omega in [100.0, 1000.0, 7000.0, 30000.0] {
        let k = state.relaxing_wavenumber(omega).unwrap();
        let z = C64::from_re(omega * gas().density / (core::f64::consts::PI * radius * radius)) / k;
        let (sin, cos) = sin_cos(k.scale(length));
        near(local_wave(&model, omega).wavenumber, k, 1e-14);
        let open = model
            .response(omega, Some(C64::ZERO), ComposeOptions::default())
            .unwrap();
        let closed = model
            .response(omega, None, ComposeOptions::default())
            .unwrap();
        near(open.z_in, C64::new(0.0, -1.0) * z * sin / cos, 2e-11);
        near(closed.z_in, C64::new(0.0, 1.0) * z * cos / sin, 2e-11);
        near(open.u_out_over_u_in, cos.recip(), 2e-11);
        assert!(open.z_in.re > 0.0 && closed.z_in.re > 0.0);
        let changed = swapped
            .response(omega, None, ComposeOptions::default())
            .unwrap();
        rate_effect = rate_effect.max((closed.z_in - changed.z_in).abs() / closed.z_in.abs());
    }
    assert!(
        rate_effect > 1e-3,
        "actual solved impedance must respond to reservoir times"
    );
}

#[test]
fn harmonic_wall_and_molecular_response_matches_original_time_domain_owner() {
    let shape = geometry(0.001, 0.05);
    let state = law(2, 1.0, gas());
    let direct = prepare(&shape, state.clone(), Some(wall()));
    let mut errors = Vec::new();
    for elements in [1, 2] {
        let dg = DgDuct::new(
            &shape,
            &gas(),
            DgMesh { elements, order: 2 },
            DgBoundary::Port(PortId(7)),
            DgBoundary::IdealOpen,
        )
        .unwrap()
        .with_viscothermal(
            wall().poles_per_field,
            wall()
                .band_rad_s
                .map(|omega| omega / core::f64::consts::TAU),
            wall().maximum_relative_error,
        )
        .unwrap()
        .with_molecular_relaxation(state.clone())
        .unwrap();
        let mut worst = 0.0_f64;
        for frequency in [100.0, 500.0, 1000.0] {
            let omega = core::f64::consts::TAU * frequency;
            let expected = direct
                .response(omega, Some(C64::ZERO), ComposeOptions::default())
                .unwrap();
            let actual = dg
                .harmonic_response(
                    omega,
                    &[(PortId(7), C64::ONE)],
                    1_000_000_000,
                    1e-10,
                    &mut || false,
                )
                .unwrap();
            worst = worst
                .max((actual.ports[0].pressure_pa - expected.z_in).abs() / expected.z_in.abs());
            assert!(actual.molecular_dissipation_w > 0.0 && actual.wall_dissipation_w > 0.0);
            assert!(actual.power_residual_w.abs() < 1e-7 * actual.input_power_w.abs());
        }
        errors.push(worst);
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[1] < 2e-3,
        "DG refinement to harmonic solution: {errors:?}"
    );
}

#[test]
fn empty_reservoirs_recover_existing_lossless_and_bessel_physics() {
    let shape = geometry(0.001, 0.05);
    let lossless = prepare(&shape, law(0, 1.0, gas()), None);
    let diffusive = prepare(&shape, law(0, 1.0, gas()), Some(wall()));
    for omega in [0.1, 10.0, 100.0, 1000.0, 6000.0] {
        let exact = crate::segment_wave(&gas(), 0.001, omega, LossModel::Lossless).unwrap();
        let retained = local_wave(&lossless, omega);
        assert_eq!(retained.wavenumber, exact.wavenumber);
        assert_eq!(
            retained.characteristic_impedance,
            exact.characteristic_impedance
        );
        let bessel = crate::segment_wave(&gas(), 0.001, omega, LossModel::Bessel).unwrap();
        let actual = local_wave(&diffusive, omega);
        near(actual.wavenumber, bessel.wavenumber, 0.01);
        near(
            actual.characteristic_impedance,
            bessel.characteristic_impedance,
            0.01,
        );
    }
}

#[test]
fn molecular_and_wall_storage_reach_one_isothermal_limit() {
    let state = law(2, 1.0, gas());
    let model = prepare(&geometry(0.001, 0.05), state, Some(wall()));
    let omega = 0.1;
    let wave = local_wave(&model, omega);
    let compliance = wave.wavenumber / wave.specific_impedance.scale(omega);
    let isothermal = 1.0 / (gas().density * gas().specific_gas_constant * gas().temperature);
    assert!((compliance.re / isothermal - 1.0).abs() < 1e-5);
    assert!(compliance.im > 0.0);
    // Changing pressure refreshes BOTH collision times and wall diffusion.
    let high_pressure = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202650.0).unwrap();
    let denser = prepare(
        &geometry(0.001, 0.05),
        law(2, 1.0, high_pressure),
        Some(wall()),
    );
    let a = local_wave(&model, 5000.0);
    let b = local_wave(&denser, 5000.0);
    assert!((a.wavenumber - b.wavenumber).abs() > 0.01);
    assert!((a.specific_impedance - b.specific_impedance).abs() > 100.0);
}

#[test]
fn repeated_cylinders_cache_diffusion_and_long_attenuating_lines_stay_finite() {
    let state = law(2, 1.0, gas());
    let single = prepare(&geometry(0.001, 0.1), state.clone(), Some(wall()));
    let split = prepare(
        &Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: 0.001,
                    length: 0.05,
                },
                Segment::Cylinder {
                    radius: 0.001,
                    length: 0.05,
                },
            ],
        },
        state.clone(),
        Some(wall()),
    );
    assert_eq!(split.blocks(), 2);
    assert_eq!(split.material_points(), 1);
    for omega in [100.0, 1000.0, 6000.0] {
        near(
            split
                .response(omega, None, ComposeOptions::default())
                .unwrap()
                .z_in,
            single
                .response(omega, None, ComposeOptions::default())
                .unwrap()
                .z_in,
            1e-11,
        );
    }
    let long = prepare(&geometry(0.001, 1e6), state, Some(wall()));
    let response = long
        .response(1000.0, Some(C64::ZERO), ComposeOptions::default())
        .unwrap();
    near(
        response.z_in,
        local_wave(&long, 1000.0).characteristic_impedance,
        1e-12,
    );
    assert_eq!(response.u_out_over_u_in, C64::ZERO);
    assert!(response.condition_estimate.is_finite());
}

#[test]
fn cones_and_profiles_use_original_geometry_and_refine_to_the_same_bore() {
    let a = 0.001;
    let b = 0.002;
    let length = 0.05;
    let state = law(2, 1.0, gas());
    let cone = prepare(
        &Duct {
            segments: vec![Segment::Cone {
                inlet_radius: a,
                outlet_radius: b,
                length,
            }],
        },
        state.clone(),
        None,
    );
    // Linear-radius Webster profile has the same exact spherical solutions.
    let profile = std::sync::Arc::new(
        crate::RadiusProfile::spline(vec![(0.0, a), (length, b)], Some([(b - a) / length; 2]))
            .unwrap(),
    );
    let omega = 5000.0;
    let expected = cone
        .response(omega, Some(C64::ZERO), ComposeOptions::default())
        .unwrap()
        .z_in;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let smooth = prepare(
            &Duct {
                segments: vec![Segment::Profile {
                    profile: profile.clone(),
                    steps: crate::RiccatiSteps(steps),
                }],
            },
            state.clone(),
            None,
        );
        let actual = smooth
            .response(omega, Some(C64::ZERO), ComposeOptions::default())
            .unwrap()
            .z_in;
        errors.push((actual - expected).abs() / expected.abs());
    }
    assert!(
        errors[1] < errors[0] / 8.0 && errors[2] < errors[1] / 8.0,
        "fourth-order Magnus refinement: {errors:?}"
    );
    assert!(errors[2] < 1e-5);
}

#[test]
fn unresolved_bands_branch_models_work_and_cancellation_refuse() {
    let shape = geometry(0.001, 0.05);
    let state = law(2, 1.0, gas());
    assert!(MolecularDuct::new(&shape, state.clone(), None, ConeSlices::Auto, 0).is_err());
    let holes = Duct {
        segments: vec![Segment::ToneHole {
            hole_radius: 0.0005,
            chimney_height: 0.001,
            bore_radius: 0.001,
            state: crate::HoleState::Closed,
        }],
    };
    assert!(MolecularDuct::new(&holes, state.clone(), None, ConeSlices::Auto, 4).is_err());
    let unresolved = MolecularWallResolution {
        poles_per_field: 1,
        maximum_relative_error: 1e-12,
        ..wall()
    };
    assert!(
        MolecularDuct::new(&shape, state.clone(), Some(unresolved), ConeSlices::Auto, 4).is_err()
    );
    assert!(
        MolecularDuct::new_controlled(
            &shape,
            state.clone(),
            Some(wall()),
            ConeSlices::Auto,
            4,
            &mut || true
        )
        .is_err()
    );
    let model = prepare(&shape, state, Some(wall()));
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, 0.01, 8000.0] {
        assert!(
            model
                .response(omega, None, ComposeOptions::default())
                .is_err()
        );
    }
    assert!(
        model
            .response_controlled(1000.0, None, ComposeOptions::default(), &mut || true)
            .is_err()
    );
    assert!(
        model
            .response(1000.0, None, ComposeOptions::default())
            .is_ok()
    );
}
