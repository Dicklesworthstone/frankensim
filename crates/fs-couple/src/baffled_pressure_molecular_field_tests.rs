//! Independent Rayleigh/one-reservoir fields, coherent nulls and physical onset.
use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

const DURATION: f64 = 0.04;
const RATE: f64 = 800.0;
const RELAXATION_TIME: f64 = 0.0008;

// A declared caloric oracle, not a measured species/air attenuation fixture.
fn material(active: bool) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            gamma: 1.25,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let species = SpeciesId::new("rayleigh-reservoir-oracle").unwrap();
    let processes = if active {
        vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * RELAXATION_TIME,
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }]
    } else {
        vec![]
    };
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes,
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap()
}

fn config(samples: usize) -> MolecularTransientConfig {
    MolecularTransientConfig {
        timestep_s: DURATION / (samples - 1) as f64,
        samples,
        transform_len: 4 * samples.next_power_of_two(),
        contour_exponent: 24.0,
    }
}

fn clock(owner: &PreparedMolecularBaffledField, epoch: f64) -> MolecularSampleClock {
    MolecularSampleClock::new(epoch, owner.config().timestep_s, 0.0).unwrap()
}

fn prepared(
    surface: &[SurfaceSample],
    inputs: usize,
    receiver: [f64; 3],
    active: bool,
    samples: usize,
) -> PreparedMolecularBaffledField {
    PreparedMolecularBaffledField::new(
        surface,
        inputs,
        receiver,
        material(active),
        config(samples),
        MolecularTransientBudget::default(),
        Options::default(),
        || false,
    )
    .unwrap()
}

// Smooth compatible onset. The velocity transform is RATE^6/(s+RATE)^7.
fn pulse(time: f64) -> InputMotion {
    if time <= 0.0 {
        return InputMotion::default();
    }
    let x = RATE * time;
    let exponential = (-x).exp() / 720.0;
    InputMotion {
        velocity: x.powi(6) * exponential,
        acceleration: RATE * (6.0 * x.powi(5) - x.powi(6)) * exponential,
    }
}

fn fields(field: &RetardedSourceField) -> [f64; 4] {
    let [x, y, z] = field.particle_velocity_m_s;
    [field.pressure_pa, x, y, z]
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-24),
        "{actual:e} != {expected:e}; scale={scale:e}, tolerance={tolerance:e}"
    );
}

fn surface() -> Vec<SurfaceSample> {
    vec![
        SurfaceSample {
            position_m: [-0.2, 0.1, 0.0],
            area_m2: 0.03,
            mode_shape: vec![1.2, -0.4],
        },
        SurfaceSample {
            position_m: [0.35, -0.2, 0.0],
            area_m2: 0.07,
            mode_shape: vec![-0.25, 0.8],
        },
    ]
}

fn modal_pulse(input: usize, time: f64) -> InputMotion {
    let value = pulse(time - if input == 0 { 0.0 } else { 0.002 });
    let amplitude = if input == 0 { 0.7 } else { -0.3 };
    InputMotion {
        velocity: amplitude * value.velocity,
        acceleration: amplitude * value.acceleration,
    }
}

#[test]
fn lossless_surface_recovers_exact_half_space_field_and_each_retarded_path() {
    let surface = surface();
    let receiver = [0.1, 0.2, 0.55];
    let owner = prepared(&surface, 2, receiver, false, 129);
    let epoch = 1e9;
    let clock = clock(&owner, epoch);
    let record = owner
        .observe_initially_resting(
            clock,
            |input, time| Ok(modal_pulse(input, (time - Dd::from_f64(epoch)).to_f64())),
            || false,
        )
        .unwrap();
    assert_eq!(record.clock(), clock);
    assert_eq!(owner.input_count(), 2);
    assert_eq!(owner.surface_sample_count(), surface.len());
    assert_ne!(owner.source_delay_s(0), owner.source_delay_s(1));
    let gas = owner.molecular_state().gas();
    let c = owner.molecular_state().frozen_sound_speed_m_s();
    let impedance = gas.density * c;
    let peak = record
        .fields()
        .iter()
        .map(|f| f.pressure_pa.abs())
        .fold(0.0, f64::max);
    assert!(peak > 0.01);
    let mut cross_terms: f64 = 0.0;
    for (index, actual) in record.fields().iter().enumerate() {
        let elapsed = Dd::from_f64(clock.timestep_s()) * Dd::from_f64(index as f64);
        let mut expected = [0.0; 4];
        let mut separate_intensity = [0.0; 3];
        for (source, sample) in surface.iter().enumerate() {
            assert_eq!(owner.source_position_m(source), Some(sample.position_m));
            let delta: [f64; 3] =
                core::array::from_fn(|axis| receiver[axis] - sample.position_m[axis]);
            let radius = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
            let time = (elapsed - Dd::from_f64(radius / c)).to_f64();
            let mut velocity = 0.0;
            let mut acceleration = 0.0;
            for (input, &shape) in sample.mode_shape.iter().enumerate() {
                let motion = modal_pulse(input, time);
                velocity += shape * motion.velocity;
                acceleration += shape * motion.acceleration;
            }
            // Independent retarded Rayleigh formulas: original half-space 2 pi,
            // with BOTH volume flow and its derivative in particle velocity.
            let green = sample.area_m2 / (core::f64::consts::TAU * radius);
            let pressure = gas.density * acceleration * green;
            expected[0] += pressure;
            for axis in 0..3 {
                let u = delta[axis] / radius * green * (velocity / radius + acceleration / c);
                expected[axis + 1] += u;
                separate_intensity[axis] += pressure * u;
            }
        }
        for (field, value) in fields(actual).into_iter().enumerate() {
            close(
                value,
                expected[field],
                peak / if field == 0 { 1.0 } else { impedance },
                8e-13,
            );
        }
        for axis in 0..3 {
            assert_eq!(
                actual.instantaneous_intensity_w_m2[axis],
                actual.pressure_pa * actual.particle_velocity_m_s[axis]
            );
            cross_terms = cross_terms
                .max((actual.instantaneous_intensity_w_m2[axis] - separate_intensity[axis]).abs());
        }
        if elapsed.to_f64()
            < owner
                .source_delay_s(0)
                .unwrap()
                .min(owner.source_delay_s(1).unwrap())
        {
            assert_eq!(fields(actual), [0.0; 4]);
        }
        assert_eq!(
            record.time_at_clock_dd(index).unwrap(),
            Dd::from_f64(epoch) + elapsed
        );
    }
    assert!(cross_terms > 1e-4 * peak * peak / impedance);
    assert!(record.time_at_clock_dd(record.fields().len()).is_none());
}

#[test]
fn relaxing_pressure_null_retains_velocity_and_original_microphone_calibration() {
    let surface = [
        SurfaceSample {
            position_m: [-0.2, 0.0, 0.0],
            area_m2: 0.04,
            mode_shape: vec![1.0],
        },
        SurfaceSample {
            position_m: [0.2, 0.0, 0.0],
            area_m2: 0.04,
            mode_shape: vec![-1.0],
        },
    ];
    let owner = prepared(&surface, 1, [0.0, 0.0, 0.5], true, 257);
    let record = owner
        .observe_initially_resting(clock(&owner, 0.0), |_, t| Ok(pulse(t.to_f64())), || false)
        .unwrap();
    let gas = owner.molecular_state().gas();
    assert_eq!(record.microphone_medium(), RayleighMedium::from(gas));
    let impedance = gas.density * gas.sound_speed;
    let frozen_impedance = gas.density * owner.molecular_state().frozen_sound_speed_m_s();
    assert!((frozen_impedance / impedance - 1.0).abs() > 0.01);
    let x = FirstOrder::new(0.0, [-1.0, 0.0, 0.0]).unwrap();
    let z = FirstOrder::new(0.0, [0.0, 0.0, -1.0]).unwrap();
    let omni = FirstOrder::new(1.0, [0.0, 0.0, -1.0]).unwrap();
    let peak = record
        .fields()
        .iter()
        .map(|f| f.particle_velocity_m_s[0].abs())
        .fold(0.0, f64::max);
    assert!(peak > 1e-6);
    for (index, field) in record.fields().iter().enumerate() {
        assert_eq!(field.pressure_pa, 0.0);
        assert_eq!(field.particle_velocity_m_s[1..], [0.0; 2]);
        assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
        close(
            record.microphone_pressure_pa(index, x).unwrap(),
            impedance * field.particle_velocity_m_s[0],
            impedance * peak,
            3e-15,
        );
        assert_eq!(record.microphone_pressure_pa(index, z).unwrap(), 0.0);
        assert_eq!(record.microphone_pressure_pa(index, omni).unwrap(), 0.0);
    }
    assert!(
        record
            .microphone_pressure_pa(record.fields().len(), x)
            .is_err()
    );
}

// Independent single-reservoir energy elimination and Cartesian Rayleigh
// gradient; does not call production wavenumbers, projections or CQ kernels.
fn laplace_field(
    material: &RelaxationState,
    surface: &[SurfaceSample],
    receiver: [f64; 3],
    s: C64,
) -> [C64; 4] {
    let gas = material.gas();
    let cv0 = gas.specific_heat_cp - gas.specific_gas_constant;
    let cvi = cv0 - gas.specific_gas_constant;
    let cp0 = cv0 + gas.specific_gas_constant;
    let cpi = cvi + gas.specific_gas_constant;
    let bulk = ((C64::from_re(cp0) + s.scale(RELAXATION_TIME * cpi))
        / (C64::from_re(cv0) + s.scale(RELAXATION_TIME * cvi)))
    .scale(gas.pressure);
    let kappa = s * bulk.recip().scale(gas.density).sqrt();
    let mut denominator = C64::ONE;
    for _ in 0..7 {
        denominator = denominator * (s + C64::from_re(RATE));
    }
    let velocity = denominator.recip().scale(RATE.powi(6));
    let mut field = [C64::ZERO; 4];
    for sample in surface {
        let delta: [f64; 3] = core::array::from_fn(|axis| receiver[axis] - sample.position_m[axis]);
        let radius = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
        let e = kappa.scale(-radius);
        let green = C64::new(e.im.cos(), e.im.sin()).scale(
            e.re.exp() * sample.area_m2 * sample.mode_shape[0] / (core::f64::consts::TAU * radius),
        );
        let source_green = velocity * green;
        field[0] = field[0] + source_green * s.scale(gas.density);
        for axis in 0..3 {
            field[axis + 1] = field[axis + 1]
                + source_green * (kappa + C64::from_re(1.0 / radius)).scale(delta[axis] / radius);
        }
    }
    field
}

fn bromwich(
    material: &RelaxationState,
    surface: &[SurfaceSample],
    receiver: [f64; 3],
    times: &[f64],
    extent: usize,
) -> Vec<[f64; 4]> {
    let sigma = 100.0;
    let dw = 20.0;
    let mut output = vec![[0.0; 4]; times.len()];
    for bin in 0..=extent {
        let omega = bin as f64 * dw;
        let values = laplace_field(material, surface, receiver, C64::new(sigma, omega));
        let weight = if bin == 0 || bin == extent { 0.5 } else { 1.0 };
        for (&time, row) in times.iter().zip(&mut output) {
            let phase = C64::new((omega * time).cos(), (omega * time).sin());
            for field in 0..4 {
                row[field] += weight * (values[field] * phase).re;
            }
        }
    }
    for (&time, row) in times.iter().zip(&mut output) {
        for value in row {
            *value *= dw * (sigma * time).exp() / core::f64::consts::PI;
        }
    }
    output
}

#[test]
fn molecular_surface_fields_refine_against_independent_material_green_inversion() {
    let surface = [
        SurfaceSample {
            position_m: [-0.2, 0.1, 0.0],
            area_m2: 0.03,
            mode_shape: vec![1.2],
        },
        SurfaceSample {
            position_m: [0.35, -0.2, 0.0],
            area_m2: 0.07,
            mode_shape: vec![-0.35],
        },
    ];
    let receiver = [0.1, 0.2, 0.55];
    let material = material(true);
    let times = [0.01, 0.015, 0.02, 0.025, 0.03];
    let reference = bromwich(&material, &surface, receiver, &times, 12000);
    let half_extent = bromwich(&material, &surface, receiver, &times, 6000);
    let impedance = material.gas().density * material.frozen_sound_speed_m_s();
    let peak = reference.iter().map(|f| f[0].abs()).fold(0.0, f64::max);
    assert!(peak > 0.01);
    for (a, b) in reference.iter().zip(half_extent) {
        for field in 0..4 {
            close(
                a[field],
                b[field],
                peak / if field == 0 { 1.0 } else { impedance },
                2e-7,
            );
        }
    }
    let mut errors = Vec::new();
    for samples in [257, 513, 1025] {
        let owner = prepared(&surface, 1, receiver, true, samples);
        let record = owner
            .observe_initially_resting(clock(&owner, 0.0), |_, t| Ok(pulse(t.to_f64())), || false)
            .unwrap();
        let mut error: f64 = 0.0;
        for (&time, expected) in times.iter().zip(&reference) {
            let index = (time / owner.config().timestep_s).round() as usize;
            for (field, actual) in fields(&record.fields()[index]).into_iter().enumerate() {
                error = error.max(
                    (actual - expected[field]).abs() * if field == 0 { 1.0 } else { impedance },
                );
            }
        }
        errors.push(error / peak);
    }
    assert!(errors[0] > 1e-7, "nonvacuous material response: {errors:?}");
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.5 * pair[0],
            "Rayleigh CQ timestep refinement: {errors:?}"
        );
    }
    assert!(
        errors[2] < 4e-4,
        "resolved complete Rayleigh field: {errors:?}"
    );
}

#[test]
fn signed_modal_projection_preserves_quiet_surface_residual() {
    let large = 2.0_f64.powi(54);
    let full = [SurfaceSample {
        position_m: [0.0; 3],
        area_m2: 0.0625,
        mode_shape: vec![large, 1.0, -large],
    }];
    let weak = [SurfaceSample {
        mode_shape: vec![1.0],
        ..full[0].clone()
    }];
    let a = prepared(&full, 3, [0.1, 0.2, 0.5], false, 65);
    let b = prepared(&weak, 1, [0.1, 0.2, 0.5], false, 65);
    let a = a
        .observe_initially_resting(clock(&a, 0.0), |_, t| Ok(pulse(t.to_f64())), || false)
        .unwrap();
    let b = b
        .observe_initially_resting(clock(&b, 0.0), |_, t| Ok(pulse(t.to_f64())), || false)
        .unwrap();
    let peak = b
        .fields()
        .iter()
        .map(|f| f.pressure_pa.abs())
        .fold(0.0, f64::max);
    assert!(peak > 0.01);
    let impedance = b.microphone_medium().density * b.microphone_medium().sound_speed;
    for (actual, expected) in a.fields().iter().zip(b.fields()) {
        for field in 0..4 {
            close(
                fields(actual)[field],
                fields(expected)[field],
                peak / if field == 0 { 1.0 } else { impedance },
                2e-14,
            );
        }
    }
}

#[test]
fn physical_onset_mapping_budgets_geometry_and_cancellation_refuse_whole_records() {
    let surface = surface();
    let receiver = [0.1, 0.2, 0.55];
    let owner = prepared(&surface, 2, receiver, false, 65);
    let make = |budget, options| {
        PreparedMolecularBaffledField::new(
            &surface,
            2,
            receiver,
            material(false),
            owner.config(),
            budget,
            options,
            || false,
        )
    };
    let exact_budget = MolecularTransientBudget {
        maximum_workspace_bytes: owner.workspace_bytes(),
        ..MolecularTransientBudget::default()
    };
    let exact_options = Options {
        maximum_projection_terms: owner.projection_terms(),
        ..Options::default()
    };
    assert!(make(exact_budget, exact_options).is_ok());
    assert!(
        make(
            MolecularTransientBudget {
                maximum_workspace_bytes: owner.workspace_bytes() - 1,
                ..exact_budget
            },
            exact_options
        )
        .is_err()
    );
    assert!(
        make(
            exact_budget,
            Options {
                maximum_projection_terms: owner.projection_terms() - 1,
                ..exact_options
            }
        )
        .is_err()
    );
    assert_eq!(
        owner.workspace_bytes() - owner.propagation.workspace_bytes(),
        surface.len() * (2 * core::mem::size_of::<Dd>() + core::mem::size_of::<[f64; 3]>())
    );
    let clock = clock(&owner, 0.0);
    for invalid_initial in [
        InputMotion {
            velocity: 1e-30,
            acceleration: 0.0,
        },
        InputMotion {
            velocity: 0.0,
            acceleration: f64::NAN,
        },
    ] {
        assert!(
            owner
                .observe_initially_resting(clock, |_, _| Ok(invalid_initial), || false)
                .is_err()
        );
    }
    // Nonzero cancelling input velocities cannot masquerade as physical rest.
    let cancelled_shapes = [SurfaceSample {
        position_m: [0.0; 3],
        area_m2: 1.0,
        mode_shape: vec![1.0, -1.0],
    }];
    let cancelled_owner = prepared(&cancelled_shapes, 2, receiver, false, 65);
    assert!(
        cancelled_owner
            .observe_initially_resting(
                clock,
                |_, _| Ok(InputMotion {
                    velocity: 1.0,
                    acceleration: 0.0
                }),
                || false
            )
            .is_err()
    );
    let jump = owner
        .observe_initially_resting(
            clock,
            |_, t| {
                Ok(InputMotion {
                    velocity: t.to_f64(),
                    acceleration: 1.0,
                })
            },
            || false,
        )
        .unwrap();
    assert!(jump.fields().iter().any(|field| field.pressure_pa != 0.0));
    let mut polls = 0;
    assert!(matches!(
        owner.observe_initially_resting(
            clock,
            |_, t| Ok(pulse(t.to_f64())),
            || {
                polls += 1;
                polls == 3
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        owner
            .observe_initially_resting(clock, |_, t| Ok(pulse(t.to_f64())), || false)
            .is_ok()
    );
    assert!(matches!(
        PreparedMolecularBaffledField::new(
            &surface,
            2,
            receiver,
            material(false),
            owner.config(),
            exact_budget,
            exact_options,
            || true
        ),
        Err(MaterialAeroacError::Aeroac(AeroacError::Cancelled { .. }))
    ));
    let wrong_clock = MolecularSampleClock::new(0.0, clock.timestep_s() * 2.0, 0.0).unwrap();
    assert!(
        owner
            .observe_initially_resting(
                wrong_clock,
                |_, _| panic!("mismatched clock must refuse before reading motion"),
                || false
            )
            .is_err()
    );
    let mut nonflat = surface.clone();
    nonflat[0].position_m[2] = 0.01;
    assert!(
        PreparedMolecularBaffledField::new(
            &nonflat,
            2,
            receiver,
            material(false),
            owner.config(),
            exact_budget,
            exact_options,
            || false
        )
        .is_err()
    );
    assert!(
        PreparedMolecularBaffledField::new(
            &surface,
            2,
            [0.0, 0.0, 0.01],
            material(false),
            owner.config(),
            exact_budget,
            exact_options,
            || false
        )
        .is_err()
    );
}
