//! Exterior material data determines the same physical load and pressure medium.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::radiation::{BaffledRadiationLoad,
    environment::bind_baffled_environment};
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctRadiationSpec, DuctReconstruction,
};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 48_000;
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101_325.0).unwrap()
}
fn graph(g: GasState, radius: f64) -> TubeNetworkSpec {
    TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.7 }],
        sections: vec![
            TubeSection { nodes: [0,1], length_m: 0.12, radius_m: radius, max_length_error_m: 0.005 },
            TubeSection { nodes: [1,2], length_m: 0.18, radius_m: radius, max_length_error_m: 0.005 }],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1<<20 }
}
fn radiation(exterior: Option<GasState>, count: usize) -> DuctRadiationSpec {
    DuctRadiationSpec {
        openings: (0..count).map(|network| DuctOpening { network, node: 2,
            center_m: [0.5*network as f64, 0.0] }).collect(), exterior,
        receiver: CircularOutletReceiver { position_m: [0.0,0.0,0.5], radial_rings: 16,
            angular_points: 64, maximum_frequency_hz: 400.0 },
        reconstruction: DuctReconstruction::Differentiating {
            mechanics_rate_hz: RATE, output_rate_hz: RATE/2 },
    }
}

#[test]
fn actual_radius_and_one_exterior_bind_identical_owner_loads_without_changing_the_bore() {
    let local = [gas(293.15), gas(315.15)]; let exterior = gas(333.15);
    let mut graphs = [graph(local[0],0.007), graph(local[1],0.005)];
    let before = graphs.clone(); let gases = [[local[0];2], [local[1];2]];
    let input = radiation(Some(exterior),2);
    let bound = bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],input.clone(),
        &CancelGate::new_clock_free()).unwrap();
    assert_eq!(bound.radiation(),&input);
    for (i,binding) in bound.bindings().iter().enumerate() {
        let expected = BaffledRadiationLoad::new(before[i].sections[1].radius_m,
            exterior.density,exterior.sound_speed,1.0/f64::from(RATE),400.0,
            &CancelGate::new_clock_free()).unwrap();
        assert_eq!(binding.load(),&expected);
        assert_eq!(graphs[i].nodes[2],expected.termination());
        assert_eq!(binding.previous_reflection(),-0.7);
        assert_eq!(*binding.interior_gas(),local[i]);
        assert_eq!(binding.opening(),input.openings[i]);
        assert_eq!(graphs[i].nodes[..2],before[i].nodes[..2]);
        assert_eq!(format!("{:?}",graphs[i].sections),format!("{:?}",before[i].sections));
        assert_eq!(graphs[i].sound_speed_m_s,before[i].sound_speed_m_s);
        assert_eq!(graphs[i].max_wave_memory_bytes,before[i].max_wave_memory_bytes);
    }
}

#[test]
fn implicit_single_port_exterior_is_resolved_from_its_local_gas_not_the_inlet() {
    let inlet = gas(293.15); let outlet = gas(315.15);
    let mut graphs = [graph(inlet,0.007)];
    let bound = bind_baffled_environment(&mut graphs,&[&[inlet,outlet]],radiation(None,1),
        &CancelGate::new_clock_free()).unwrap();
    assert_eq!(bound.radiation().exterior,Some(outlet));
    assert_eq!(bound.bindings()[0].load().density_kg_m3(),outlet.density);
    assert_eq!(bound.bindings()[0].load().sound_speed_m_s(),outlet.sound_speed);
}

#[test]
fn later_boundary_failure_reactive_replacement_and_cancel_leave_all_graphs_unchanged() {
    let g=gas(293.15); let gases=[[g;2];2];
    let mut graphs=[graph(g,0.007),graph(g,0.08)];
    let before=format!("{graphs:?}");
    // The receiver spatial rule admits both disks, but the second load violates
    // ka<=0.5. The first successful load must not already have changed its node.
    assert!(bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],radiation(Some(g),2),
        &CancelGate::new_clock_free()).is_err());
    assert_eq!(format!("{graphs:?}"),before);
    graphs[1]=graph(g,0.005);
    let existing=BaffledRadiationLoad::new(0.005,g.density,g.sound_speed,
        1.0/f64::from(RATE),400.0,&CancelGate::new_clock_free()).unwrap();
    graphs[1].nodes[2]=existing.termination(); let before=format!("{graphs:?}");
    assert!(bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],radiation(Some(g),2),
        &CancelGate::new_clock_free()).is_err());
    assert_eq!(format!("{graphs:?}"),before);
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],radiation(Some(g),2),&gate).is_err());
    assert_eq!(format!("{graphs:?}"),before);
}

#[test]
fn absent_mappings_duplicate_ports_and_pressure_mismatch_do_not_become_a_boundary_model() {
    let g=gas(293.15);let gases=[g;2];let mut graphs=[graph(g,0.007)];
    let original=format!("{graphs:?}");
    let gate=CancelGate::new_clock_free();
    assert!(bind_baffled_environment(&mut graphs,&[],radiation(None,1),&gate).is_err());
    let mut duplicate=radiation(Some(g),1);duplicate.openings.push(duplicate.openings[0]);
    assert!(bind_baffled_environment(&mut graphs,&[&gases],duplicate,&gate).is_err());
    let wrong=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,90000.0).unwrap();
    assert!(bind_baffled_environment(&mut graphs,&[&gases],radiation(Some(wrong),1),&gate).is_err());
    assert_eq!(format!("{graphs:?}"),original);
}

#[path = "baffled_environment/instruments.rs"]
mod instruments;

mod harmonic {
    use fs_bem::near_field::FirstOrder;
    use fs_couple::air_path::GasAbsorptionModel;
    use fs_couple::pcm_wav::baffled::gas::{GasPropagation, GasPropagationLimits};
    use fs_couple::pcm_wav::baffled::harmonic::{Field, Options, Prepared};
    use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample};
    use fs_material::gas::{ConductivityModel, GasSpec, GasState};
    use fs_material::relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxationState, RelaxingMedium, SpeciesRelaxation,
    };
    use fs_math::c64::C64;
    use fs_thermochem::SpeciesId;
    use std::f64::consts::TAU;

    // Manufactured caloric material, not measured species data: cv0=4R,
    // cvinf=3R and one rotational reservoir C=R with tau=1 ms.
    fn material(reservoir: bool) -> RelaxationState {
        let spec = GasSpec {
            gamma: 1.25,
            conductivity: ConductivityModel::Eucken,
            ..GasSpec::dry_air_ussa1976()
        };
        let gas = GasState::try_new(&spec, 293.15, 101_325.0).unwrap();
        let species = SpeciesId::new("A").unwrap();
        let processes = if reservoir {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Rotational {
                    degrees_of_freedom: 2,
                },
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 0.001 * gas.pressure,
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 250.0,
                    t_max_k: 350.0,
                },
            }]
        } else {
            Vec::new()
        };
        RelaxationAbsorption::new(vec![SpeciesRelaxation {
            species: species.clone(),
            processes,
        }])
        .unwrap()
        .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap()
    }

    // Independent real arithmetic for k² = w²/(RT) * (4-3ix)/(5-4ix)
    // / (1-iy), x=w*tau, y=w*delta. No production complex square root or
    // wavenumber evaluator is used for the oracle.
    fn wave(gas: &GasState, omega: f64, reservoir: bool) -> C64 {
        let x = omega * 0.001;
        let (a, b) = if reservoir {
            (
                (20.0 + 12.0 * x * x) / (25.0 + 16.0 * x * x),
                x / (25.0 + 16.0 * x * x),
            )
        } else {
            (0.8, 0.0)
        };
        let delta = (4.0 * gas.dynamic_viscosity / 3.0
            + (gas.gamma - 1.0) * gas.thermal_conductivity / gas.specific_heat_cp)
            / (gas.density * gas.sound_speed * gas.sound_speed);
        let y = omega * delta;
        let factor = omega * omega / (gas.specific_gas_constant * gas.temperature);
        let re = factor * (a - b * y) / (1.0 + y * y);
        let im = factor * (b + a * y) / (1.0 + y * y);
        let real = (0.5 * (re.hypot(im) + re)).sqrt();
        C64::new(real, im / (2.0 * real))
    }

    fn patch(x: f64, shape: Vec<f64>) -> SurfaceSample {
        SurfaceSample {
            position_m: [x, 0.0, 0.0],
            area_m2: 0.125,
            mode_shape: shape,
        }
    }
    fn close(actual: C64, expected: C64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1e-25),
            "actual={actual:?}, expected={expected:?}"
        );
    }
    fn green(k: C64, radius: f64) -> C64 {
        C64::new((k.re * radius).cos(), (k.re * radius).sin())
            .scale((-k.im * radius).exp() / (TAU * radius))
    }

    #[test]
    fn molecular_harmonics_keep_dispersion_gradient_and_absorption_exactly_once() {
        let surface = [patch(0.0, vec![8.0])]; // one unit of volume flow
        let radii = [0.5, 3.0];
        for reservoir in [false, true] {
            let source = material(reservoir);
            let gas = source.gas();
            for frequency in [20.0, 160.0, 2000.0, 16_000.0] {
                let omega = TAU * frequency;
                let k = wave(gas, omega, reservoir);
                let prepared = Prepared::new_relaxing(
                    &surface,
                    1,
                    &radii.map(|r| [0.0, 0.0, r]),
                    frequency,
                    &source,
                    Options::default(),
                    || false,
                )
                .unwrap();
                assert_eq!(prepared.relaxation(), Some(&source));
                assert_eq!(prepared.medium(), RayleighMedium::from(gas));
                close(prepared.wavenumber_per_m(), k, 2e-14);
                let mut power = [0.0; 2];
                for (i, radius) in radii.into_iter().enumerate() {
                    let field = prepared.transfer(i, 0).unwrap();
                    let pressure = green(k, radius) * C64::new(0.0, -gas.density * omega);
                    let velocity = pressure
                        * C64::new(k.re, k.im + 1.0 / radius).scale(1.0 / (gas.density * omega));
                    close(field.pressure_pa, pressure, 2e-12);
                    close(field.particle_velocity_m_s[2], velocity, 2e-12);
                    close(
                        field.pressure_gradient_pa_m[2],
                        pressure * C64::new(-k.im - 1.0 / radius, k.re),
                        2e-12,
                    );
                    power[i] = radius
                        * radius
                        * (field.pressure_pa * field.particle_velocity_m_s[2].conj()).re;
                    assert!(power[i] > 0.0);
                }
                let attenuation = (-2.0 * k.im * (radii[1] - radii[0])).exp();
                assert!((power[1] / power[0] - attenuation).abs() < 2e-12);
            }
        }
    }

    #[test]
    fn molecular_signed_paths_preserve_spatial_gradient_and_microphone_reference() {
        let source = material(true);
        let surface = [patch(-0.3, vec![1.0, -0.2]), patch(0.25, vec![-0.6, 0.8])];
        let center = [0.11, 0.13, 0.4];
        let step = 1e-6;
        let mut positions = vec![center];
        for axis in 0..3 {
            let mut plus = center;
            let mut minus = center;
            plus[axis] += step;
            minus[axis] -= step;
            positions.extend([plus, minus]);
        }
        let prepared = Prepared::new_relaxing(
            &surface,
            2,
            &positions,
            1600.0,
            &source,
            Options::default(),
            || false,
        )
        .unwrap();
        let inputs = [C64::new(0.3, -0.7), C64::new(-0.2, 0.4)];
        let mut fields = vec![Field::default(); positions.len()];
        prepared.evaluate_into(&inputs, &mut fields).unwrap();
        let omega = TAU * 1600.0;
        let k = wave(source.gas(), omega, true);
        let mut expected = C64::ZERO;
        for sample in &surface {
            let radius = (0..3)
                .map(|axis| (center[axis] - sample.position_m[axis]).powi(2))
                .sum::<f64>()
                .sqrt();
            let flow = (inputs[0].scale(sample.mode_shape[0])
                + inputs[1].scale(sample.mode_shape[1]))
            .scale(sample.area_m2);
            expected =
                expected + flow * green(k, radius) * C64::new(0.0, -source.gas().density * omega);
        }
        close(fields[0].pressure_pa, expected, 2e-13);
        for axis in 0..3 {
            let difference = (fields[1 + 2 * axis].pressure_pa - fields[2 + 2 * axis].pressure_pa)
                .scale(0.5 / step);
            close(fields[0].pressure_gradient_pa_m[axis], difference, 2e-8);
        }
        let pattern = FirstOrder::new(0.5, [0.0, 0.0, -1.0]).unwrap();
        let mut heard = vec![C64::ZERO; positions.len()];
        prepared
            .observe_into(&inputs, &vec![pattern; positions.len()], &mut heard)
            .unwrap();
        close(
            heard[0],
            fields[0].pressure_pa.scale(0.5)
                + fields[0].particle_velocity_m_s[2]
                    .scale(0.5 * source.gas().density * source.gas().sound_speed),
            1e-14,
        );
        let before = heard.clone();
        assert!(
            prepared
                .observe_into(
                    &[C64::new(f64::NAN, 0.0); 2],
                    &vec![pattern; positions.len()],
                    &mut heard
                )
                .is_err()
        );
        assert_eq!(heard, before);
        assert!(
            Prepared::new_relaxing(
                &surface,
                2,
                &positions,
                1600.0,
                &source,
                Options {
                    maximum_kernel_terms: 27,
                    ..Options::default()
                },
                || false
            )
            .is_err()
        );
        assert!(
            Prepared::new_relaxing(
                &surface,
                2,
                &positions,
                1600.0,
                &source,
                Options::default(),
                || true
            )
            .is_err()
        );
    }

    #[test]
    fn molecular_harmonics_agree_with_the_original_causal_gas_pressure_observer() {
        let source = material(true);
        let surface = [patch(-0.15, vec![1.0, -0.4]), patch(0.2, vec![-0.3, 0.7])];
        let position = [0.12, 0.0, 0.5];
        let inputs = [C64::new(0.003, 0.002), C64::new(-0.002, 0.001)];
        let propagation = GasPropagation {
            law: GasAbsorptionModel::MolecularRelaxation,
            relaxation: Some(source.clone()),
            limits: GasPropagationLimits::surface_reference(),
        };
        let (base, _) = BaffledPressure::new_multirate_gas_differentiating(
            &surface,
            2,
            48_000,
            48_000,
            position,
            source.gas(),
            &propagation,
            &mut || false,
        )
        .unwrap();
        for frequency in [200.0, 1600.0, 8000.0] {
            let prepared = Prepared::new_relaxing(
                &surface,
                2,
                &[position],
                frequency,
                &source,
                Options::default(),
                || false,
            )
            .unwrap();
            let mut field = [Field::default()];
            prepared.evaluate_into(&inputs, &mut field).unwrap();
            let mut causal = base.clone();
            let latency = causal.filter_delay_samples();
            let step = TAU * frequency / 48_000.0;
            let (mut error, mut signal) = (0.0, 0.0);
            for n in 0..4800 {
                let phase = step * n as f64;
                let velocity = inputs.map(|v| v.re * phase.cos() + v.im * phase.sin());
                let heard = causal.step(&velocity).unwrap();
                if n >= 2400 {
                    let phase = step * (n as f64 - latency);
                    let expected = field[0].pressure_pa.re * phase.cos()
                        + field[0].pressure_pa.im * phase.sin();
                    error += (heard - expected).powi(2);
                    signal += expected * expected;
                }
            }
            // The original causal filter certifies retained response energy,
            // not machine-precision agreement with a continuous phasor.
            assert!(
                (error / signal).sqrt() < 3e-3,
                "frequency={frequency}, relative RMS={}",
                (error / signal).sqrt()
            );
        }
    }
}
