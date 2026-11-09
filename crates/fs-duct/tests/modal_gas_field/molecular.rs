//! G1/G3/G4: original molecular law in the continuous multimodal field.
use super::*;
use fs_duct::compose::ComposeOptions;
use fs_duct::modal::j1_roots;
use fs_duct::molecular::{MolecularDuct, MolecularWallResolution};
use fs_duct::ConeSlices;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    RelaxationState, RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::det;
use fs_thermochem::SpeciesId;

fn law(active: bool, maximum_temperature: f64) -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        [("A", 900.0), ("B", 1300.0)]
            .into_iter()
            .map(|(id, theta_k)| SpeciesRelaxation {
                species: SpeciesId::new(id).unwrap(),
                processes: if active {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 101_325.0 * 3e-5,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: maximum_temperature,
                        },
                    }]
                } else {
                    vec![]
                },
            })
            .collect(),
    )
    .unwrap()
}

fn original_state(field: &AxialGasField, law: &RelaxationAbsorption, x: f64) -> RelaxationState {
    law.evaluate(&RelaxingMedium::from_thermochemical(&field.state_at(x).unwrap()).unwrap())
        .unwrap()
}

// Independent elimination of the original heat-capacity reservoirs. No native
// wave, diffusion spectrum, modal propagator or matched-boundary map is used.
fn wave(state: &RelaxationState, omega: f64, radius: f64, gamma: f64) -> (C64, C64, C64) {
    let gas = state.gas();
    let cv = state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k()),
        |sum, p| {
            sum + C64::from_re(p.heat_capacity_j_kg_k)
                / C64::new(1.0, -omega * p.relaxation_time_s)
        },
    );
    let bulk = ((cv + C64::from_re(gas.specific_gas_constant)) / cv).scale(gas.pressure);
    let transverse = gamma / radius;
    let k = (bulk.recip().scale(gas.density * omega * omega)
        - C64::from_re(transverse * transverse))
    .sqrt();
    assert!(k.re > 0.0 && k.im > 0.0);
    let area = core::f64::consts::PI * radius * radius;
    let z = C64::from_re(gas.density * omega / area) / k;
    (k, z, bulk)
}

fn phase(k: C64, x: f64) -> C64 {
    let (s, c) = det::sin_cos(k.re * x);
    C64::new(c, s).scale(det::exp(-k.im * x))
}

fn backward(k: C64, z: C64, length: f64, endpoint: [C64; 2]) -> [C64; 2] {
    let forward = phase(k, length);
    let reverse = phase(k, -length);
    let cos = (forward + reverse).scale(0.5);
    let minus_i_sin = (reverse - forward).scale(0.5);
    [
        cos * endpoint[0] + minus_i_sin * z * endpoint[1],
        minus_i_sin / z * endpoint[0] + cos * endpoint[1],
    ]
}

#[test]
fn matched_modal_phase_cutoff_and_heat_follow_the_original_caloric_reservoirs() {
    let length = 0.04;
    let radius = 0.03;
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let field = field(vec![0.0, length], &[300.0, 300.0], &[0.8, 0.8], false);
    let law = law(true, 1000.0);
    let source = original_state(&field, &law, 0.0);
    let roots = j1_roots(2).unwrap();
    let gammas = [0.0, roots[0], roots[1]];
    // The first radial mode lies near its equilibrium cutoff; the second is
    // evanescent. Molecular dispersion must alter their phase and induction.
    let omega = source.gas().sound_speed * (1.02 * roots[0] / radius);
    let binding = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(64))
        .unwrap()
        .with_molecular_relaxation(law.clone(), None, &mut || false)
        .unwrap();
    assert_eq!(binding.relaxation_law(), &law);
    assert_eq!(binding.material_states(), &[source.clone()]);
    assert_eq!(binding.outlet_relaxation(), &source);
    let prepared = binding
        .prepare(omega, MouthLoad::Matched, ModeCount(3), &mut || false)
        .unwrap();
    assert!(std::ptr::eq(
        binding.material_states().as_ptr(),
        prepared.molecular_field().unwrap().material_states().as_ptr(),
    ));
    let zref = prepared.reference_impedance_pa_s_m3();
    let drives = [
        PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::new(2.0, -0.5),
            impedance_pa_s_m3: C64::new(0.7, -0.2).scale(zref),
        },
        PortExcitation::VolumeFlow(C64::new(1e-6, -2e-6)),
        PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::new(-0.4, 0.9),
            impedance_pa_s_m3: C64::new(0.3, 0.8).scale(zref),
        },
    ];
    let result = prepared.drive_sources(&drives, &mut || false).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let mut heat = 0.0;
    for (mode, gamma) in gammas.into_iter().enumerate() {
        let (k, z, bulk) = wave(&source, omega, radius, gamma);
        let u0 = match drives[mode] {
            PortExcitation::VolumeFlow(u) => u,
            PortExcitation::Pressure {
                open_circuit_pressure_pa: p,
                impedance_pa_s_m3: zs,
            } => p / (z + zs),
        };
        let p0 = z * u0;
        for (station, &x) in result.field.x.iter().enumerate() {
            close(result.field.flow[station][mode], u0 * phase(k, x), u0.abs(), 3e-10);
            close(result.field.pressure[station][mode], p0 * phase(k, x), p0.abs(), 3e-10);
        }
        // Orthogonality gives int |p|^2 dA = A sum |P_n|^2. The
        // original local heat is omega Im(1/K)|p|^2/2, including evanescence.
        heat += 0.5 * omega * bulk.recip().im * area * p0.abs().powi(2)
            * (1.0 - det::exp(-2.0 * k.im * length)) / (2.0 * k.im);
    }
    let outward = power(result.field.pressure.last().unwrap(), result.field.flow.last().unwrap());
    let absorbed = result.inlet.inlet_power_w - outward;
    assert!(heat > 0.0 && outward > 0.0);
    assert!((absorbed - heat).abs() < 2e-9 * heat);
    assert!((result.inlet.generator_power_w - result.inlet.source_dissipation_w
        - result.inlet.inlet_power_w).abs() < 2e-13);
    let equilibrium = binding.original()
        .prepare(omega, LossModel::Lossless, MouthLoad::Matched, ModeCount(3), &mut || false)
        .unwrap()
        .drive_sources(&drives, &mut || false)
        .unwrap();
    assert!((result.field.pressure[0][1] - equilibrium.field.pressure[0][1]).abs()
        > 1e-3 * result.field.pressure[0][1].abs());
}

#[test]
fn original_material_jumps_refract_each_mode_against_independent_layer_fields() {
    let length = 0.06;
    let radius = 0.03;
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let field = field(
        vec![0.0, length / 2.0, length],
        &[300.0, 300.0, 540.0],
        &[0.8, 0.8, 0.2],
        true,
    );
    let law = law(true, 1000.0);
    let original = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(48)).unwrap();
    let count = original.steps_used();
    let binding = original
        .with_molecular_relaxation(law.clone(), None, &mut || false)
        .unwrap();
    assert_eq!(binding.original().steps_used(), count);
    for step in 0..count {
        for node in 0..2 {
            let x = binding.quadrature_position_m(step, node).unwrap();
            assert_eq!(
                binding.relaxation_at_quadrature(step, node).unwrap(),
                &original_state(&field, &law, x),
            );
        }
    }
    let states = [
        original_state(&field, &law, length / 4.0),
        original_state(&field, &law, 3.0 * length / 4.0),
    ];
    assert_eq!(binding.outlet_relaxation(), &states[1]);
    let omega = 55_000.0;
    let prepared = binding
        .prepare(omega, MouthLoad::Matched, ModeCount(2), &mut || false)
        .unwrap();
    let drives = flows(2);
    let actual = prepared.drive_sources(&drives, &mut || false).unwrap();
    assert!(actual.field.x.contains(&(length / 2.0)));
    for (mode, gamma) in [0.0, j1_roots(1).unwrap()[0]].into_iter().enumerate() {
        let (k1, z1, _) = wave(&states[0], omega, radius, gamma);
        let (k2, z2, _) = wave(&states[1], omega, radius, gamma);
        let mouth = [z2, C64::ONE];
        let interface = backward(k2, z2, length / 2.0, mouth);
        let inlet = backward(k1, z1, length / 2.0, interface);
        let PortExcitation::VolumeFlow(flow) = drives[mode] else { unreachable!() };
        let scale = flow / inlet[1];
        for (station, &x) in actual.field.x.iter().enumerate() {
            let state = if x >= length / 2.0 {
                backward(k2, z2, length - x, mouth)
            } else {
                backward(k1, z1, length / 2.0 - x, interface)
            };
            close(actual.field.pressure[station][mode], state[0] * scale, 1.0, 4e-10);
            close(actual.field.flow[station][mode], state[1] * scale, flow.abs(), 4e-10);
        }
    }
    let outward = power(actual.field.pressure.last().unwrap(), actual.field.flow.last().unwrap());
    assert!(actual.inlet.inlet_power_w > outward && outward > 0.0);
}

#[test]
fn coupled_wall_heat_uses_the_original_plane_owner_and_full_modal_mouth_matrix() {
    let length = 0.04;
    let profile = RadiusProfile::exponential(0.001, 2.0, length).unwrap();
    let field = field(vec![0.0, length], &[290.0, 350.0], &[0.8, 0.4], false);
    let wall = MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [100.0, 4000.0],
        maximum_relative_error: 0.5,
    };
    let resolution = RiccatiSteps(48);
    for active in [false, true] {
        let law = law(active, 1000.0);
        let binding = ModalGasFieldProfile::new(&profile, &field, resolution)
            .unwrap()
            .with_molecular_relaxation(law.clone(), Some(wall), &mut || false)
            .unwrap();
        assert_eq!(binding.wall_resolution(), Some(wall));
        let duct = Duct {
            segments: vec![Segment::Profile {
                profile: Arc::new(profile.clone()),
                steps: resolution,
            }],
        };
        let plane = MolecularDuct::from_gas_field(
            &duct, &field, &law, Some(wall), resolution, ConeSlices::Fixed(8), 4096,
        ).unwrap();
        for omega in [200.0, 1000.0, 3000.0] {
            let expected = plane.response(omega, Some(C64::ZERO), ComposeOptions::default()).unwrap();
            let actual = binding.prepare(
                omega, MouthLoad::Analytic(Termination::IdealOpen), ModeCount(1), &mut || false,
            ).unwrap().drive_sources(&flows(1), &mut || false).unwrap();
            let input_flow = actual.field.flow[0][0];
            close(actual.field.pressure[0][0] / input_flow, expected.z_in, 1.0, 4e-9);
            close(actual.field.flow.last().unwrap()[0] / input_flow, expected.u_out_over_u_in, 1.0, 4e-9);
            assert!(actual.inlet.inlet_power_w > 0.0);
        }
        assert!(binding.prepare(8000.0, MouthLoad::Matched, ModeCount(1), &mut || false).is_err());
    }

    // A wider flared source drives several coupled modes and the complete
    // exterior matrix. The prepared field keeps every signed mouth cross term.
    let profile = RadiusProfile::exponential(0.025, 6.0, length).unwrap();
    let binding = ModalGasFieldProfile::new(&profile, &field, resolution).unwrap()
        .with_molecular_relaxation(law(true, 1000.0), None, &mut || false).unwrap();
    let omega = 60_000.0;
    let modes = ModeCount(3);
    let exterior = *field.state_at(length).unwrap().state();
    let matrix = fs_duct::flanged::flanged_modal_impedance(
        &exterior, profile.outlet_radius(), omega, modes.0,
    ).unwrap();
    let load = ModalMouthImpedance::new(
        omega, profile.outlet_radius(), exterior, modes, matrix.clone(),
    ).unwrap();
    let full = binding.prepare_with_flanged_mouth(&exterior, omega, modes, &mut || false).unwrap();
    let declared = binding.prepare_with_mouth_matrix(omega, &load, modes, &mut || false).unwrap();
    let actual = full.drive_sources(&flows(3), &mut || false).unwrap();
    same_fields(&actual, &declared.drive_sources(&flows(3), &mut || false).unwrap());
    let (p, u) = (actual.field.pressure.last().unwrap(), actual.field.flow.last().unwrap());
    for i in 0..3 {
        let expected = (0..3).fold(C64::ZERO, |sum, j| sum + matrix[3 * i + j] * u[j]);
        close(p[i], expected, 1.0, 3e-10);
    }
    let outward = power(p, u);
    assert!(actual.inlet.inlet_power_w > outward && outward > 0.0);
}

#[test]
fn absent_reservoirs_recover_original_fields_and_refusals_leave_the_source_reusable() {
    let length = 0.06;
    let profile = RadiusProfile::exponential(0.025, 5.0, length).unwrap();
    for temperatures in [[300.0, 300.0], [300.0, 540.0]] {
        let field = field(vec![0.0, length], &temperatures, &[0.8, 0.8], false);
        let original = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(48)).unwrap();
        let binding = original.clone()
            .with_molecular_relaxation(law(false, 1000.0), None, &mut || false).unwrap();
        for load in [
            MouthLoad::Matched,
            MouthLoad::Analytic(Termination::Closed),
            MouthLoad::Analytic(Termination::IdealOpen),
        ] {
            let plain = original.prepare(50_000.0, LossModel::Lossless, load, ModeCount(3), &mut || false).unwrap();
            let molecular = binding.prepare(50_000.0, load, ModeCount(3), &mut || false).unwrap();
            let expected = plain.drive_sources(&flows(3), &mut || false).unwrap();
            same_fields(&expected, &molecular.drive_sources(&flows(3), &mut || false).unwrap());
            let refs = [plain.reference_impedance_pa_s_m3(); 3];
            exact_complex(&plain.at_references(&refs).unwrap(), &molecular.at_references(&refs).unwrap());
        }
        assert!(original.clone().with_molecular_relaxation(law(true, 1000.0), None, &mut || true).is_err());
        assert!(binding.prepare(50_000.0, MouthLoad::Matched, ModeCount(3), &mut || true).is_err());
        assert!(binding.prepare(0.0, MouthLoad::Matched, ModeCount(3), &mut || false).is_err());
        assert!(binding.prepare(50_000.0, MouthLoad::Matched, ModeCount(0), &mut || false).is_err());
        let source_temperatures: Vec<_> = binding.material_states().iter().map(|s| s.gas().temperature).collect();
        assert!(binding.prepare(50_000.0, MouthLoad::Matched, ModeCount(3), &mut || false).is_ok());
        assert_eq!(
            source_temperatures,
            binding.material_states().iter().map(|s| s.gas().temperature).collect::<Vec<_>>(),
        );
        if temperatures[1] > 500.0 {
            assert!(original.with_molecular_relaxation(law(true, 500.0), None, &mut || false).is_err());
        }
    }
}
