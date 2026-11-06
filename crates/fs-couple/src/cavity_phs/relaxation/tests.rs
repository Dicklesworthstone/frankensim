use super::*;
use fs_material::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use fs_qty::Pressure;
use fs_scenario::NeckEnd;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;

fn material(pressure: f64, quality: f64) -> RelaxingRegionalWater {
    let (lo, hi) = if pressure < 1e6 {
        (300.0, 500.0)
    } else {
        (625.0, 645.0)
    };
    let law = RegionalWater::new(
        Pressure::new(pressure),
        RegionalWaterControls {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap();
    RelaxingRegionalWater::try_new(law.state_at_quality(quality, || false).unwrap(), 0.003).unwrap()
}

fn geometry() -> HelmholtzCavity {
    HelmholtzCavity {
        volume_m3: 0.003,
        neck_radius_m: 0.003,
        neck_length_m: 0.1,
        inner_end: NeckEnd::Flanged,
        outer_end: NeckEnd::Unflanged,
    }
}

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 8,
        max_omega_dt: 2.5,
    }
}

fn near(actual: f64, expected: f64, scale: f64, relative: f64) {
    let error = (actual - expected).abs() / scale.abs().max(expected.abs()).max(1e-30);
    assert!(
        actual.is_finite() && error <= relative,
        "{actual:e} != {expected:e}; scaled error {error:e}"
    );
}

fn independent_inertance(rho: f64, geometry: HelmholtzCavity) -> f64 {
    // Independent declared end-correction reference, not the native constructor.
    let length = geometry.neck_length_m + (0.8216 + 0.6133) * geometry.neck_radius_m;
    rho * length / (core::f64::consts::PI * geometry.neck_radius_m.powi(2))
}

#[test]
fn g1_source_density_sets_native_inertia_and_frozen_relaxed_compression() {
    for pressure in [1e5, 20e6] {
        for quality in [0.2, 0.6] {
            let material = material(pressure, quality);
            let original = *material.source();
            let geometry = geometry();
            let rho = original.bulk_density_kg_m3();
            let cavity = RelaxingWaterCavity::try_new(material, geometry, 0.0).unwrap();
            assert_eq!(cavity.material(), &material);
            assert_eq!(cavity.material().source(), &original);
            assert_eq!(cavity.geometry(), &geometry);
            assert_eq!(
                (cavity.system().state_dim(), cavity.system().port_dim()),
                (3, 1)
            );
            near(
                cavity.neck_inertance_pa_s2_m3(),
                independent_inertance(rho, geometry),
                1.0,
                2e-14,
            );

            let q = 1e-7 * geometry.volume_m3;
            let frozen = cavity.initial_state(q, 0.0, Some(0.0)).unwrap();
            let relaxed = cavity.initial_state(q, 0.0, None).unwrap();
            let a = cavity.observe(&frozen).unwrap();
            let b = cavity.observe(&relaxed).unwrap();
            near(
                a.pressure_increment_pa,
                material.frozen().state().adiabatic_bulk_modulus_pa() * q / geometry.volume_m3,
                1.0,
                3e-14,
            );
            near(
                b.pressure_increment_pa,
                material.equilibrium().state().adiabatic_bulk_modulus_pa() * q / geometry.volume_m3,
                1.0,
                3e-14,
            );
            assert!(a.pressure_increment_pa > b.pressure_increment_pa);
            near(
                a.density_increment_kg_m3,
                rho * q / geometry.volume_m3,
                rho * 1e-7,
                3e-14,
            );
            near(
                a.excess_mass_kg,
                geometry.volume_m3 * a.density_increment_kg_m3,
                rho * q,
                3e-14,
            );
            assert_eq!(a.density_increment_kg_m3, b.density_increment_kg_m3);
            assert_eq!(a.pressure_increment_pa, cavity.system().output(&frozen)[0]);
            assert_eq!(b.compression_dissipation_w, 0.0);
            assert!(a.compression_dissipation_w > 0.0);

            let flow = 3e-7;
            let moving = cavity.initial_state(0.0, flow, None).unwrap();
            let observation = cavity.observe(&moving).unwrap();
            near(observation.neck_outflow_m3_s, flow, flow, 2e-14);
            near(
                observation.stored_energy_j,
                0.5 * independent_inertance(rho, geometry) * flow * flow,
                1e-20,
                3e-14,
            );
            assert_eq!(
                observation.linearized_density_kg_m3.to_bits(),
                rho.to_bits()
            );
        }
    }
}

// Independent physical integration in excess mass, outward neck flow and
// relaxed volume, with heat/work integrated as physical scalars. No pHS
// matrix, energy-normalized memory or native time propagator enters this oracle.
fn reference(
    mut state: [f64; 5],
    duration: f64,
    steps: usize,
    material: RelaxingRegionalWater,
    geometry: HelmholtzCavity,
    resistance: f64,
    input: f64,
) -> [f64; 5] {
    let rho = material.source().bulk_density_kg_m3();
    let inertia = independent_inertance(rho, geometry);
    let equilibrium = material.equilibrium().state().adiabatic_bulk_modulus_pa();
    let strength = material.arm().delta_bulk_modulus_pa;
    let tau = material.arm().relaxation_time_s;
    let rhs = |x: &[f64; 5]| {
        let q = x[0] / rho;
        let slip = q - x[2];
        let pressure = (equilibrium * q + strength * slip) / geometry.volume_m3;
        [
            rho * (input - x[1]),
            (pressure - resistance * x[1]) / inertia,
            slip / tau,
            resistance * x[1] * x[1] + strength * slip * slip / (geometry.volume_m3 * tau),
            input * pressure,
        ]
    };
    let h = duration / steps as f64;
    for _ in 0..steps {
        let a = rhs(&state);
        let b = rhs(&core::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = rhs(&core::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = rhs(&core::array::from_fn(|i| state[i] + h * c[i]));
        for i in 0..5 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn g1_native_prepared_flow_matches_independent_mass_momentum_heat_and_work() {
    for pressure in [1e5, 20e6] {
        let material = material(pressure, 0.2);
        let geometry = geometry();
        let rho = material.source().bulk_density_kg_m3();
        let inertia = independent_inertance(rho, geometry);
        let stiffness =
            material.equilibrium().state().adiabatic_bulk_modulus_pa() / geometry.volume_m3;
        let resistance = 0.2 * (inertia * stiffness).sqrt();
        let cavity = RelaxingWaterCavity::try_new(material, geometry, resistance).unwrap();
        let q0 = geometry.volume_m3 * 1e-6;
        let input = 1e-7;
        let duration = material.arm().relaxation_time_s;
        let dt = duration / 64.0;
        let mut prepared = cavity.prepare(dt, budget()).unwrap();
        let mut state = cavity.initial_state(q0, 0.0, Some(0.0)).unwrap();
        let initial_energy = cavity.system().hamiltonian(&state);
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..64 {
            let mut next = [0.0; 3];
            let mut mean_pressure = [0.0];
            let record = prepared
                .step_into(&state, &[input], &mut next, &mut mean_pressure)
                .unwrap();
            near(
                record.supplied,
                dt * input * mean_pressure[0],
                initial_energy,
                3e-12,
            );
            heat += record.dissipated;
            work += record.supplied;
            state = next;
        }
        let observed = cavity.observe(&state).unwrap();
        let initial = [rho * q0, 0.0, 0.0, 0.0, 0.0];
        let coarse = reference(
            initial, duration, 2048, material, geometry, resistance, input,
        );
        let fine = reference(
            initial, duration, 4096, material, geometry, resistance, input,
        );
        let scales = [rho * q0, q0 / duration, q0, initial_energy, initial_energy];
        let actual = [
            observed.excess_mass_kg,
            observed.neck_outflow_m3_s,
            state[2] / cavity.memory_scale,
            heat,
            work,
        ];
        for i in 0..5 {
            near(coarse[i], fine[i], scales[i], 2e-9);
            near(actual[i], fine[i], scales[i], 2e-8);
        }
        near(
            observed.stored_energy_j - initial_energy + heat,
            work,
            initial_energy,
            3e-9,
        );
        assert!(heat > 0.0);
        assert_eq!(cavity.material(), &material);
    }
}

#[test]
fn g3_native_mechanical_transformer_preserves_pressure_flow_power_and_total_energy() {
    let material = material(1e5, 0.3);
    let geometry = geometry();
    let resistance = 1000.0;
    let cavity = RelaxingWaterCavity::try_new(material, geometry, resistance).unwrap();
    let inertia = cavity.neck_inertance_pa_s2_m3();
    let root = cavity.memory_scale;
    let (mass, stiffness, damping, area) = (0.02, 2000.0, 0.2, 1e-4);
    let (displacement, velocity, q, flow) = (3e-7, 2e-4, 2e-7 * geometry.volume_m3, 1e-8);
    let initial_cavity = cavity
        .initial_state(q, flow, Some(0.25 * root * q))
        .unwrap();
    let observed = cavity.observe(&initial_cavity).unwrap();
    let (retained, network) = cavity.into_parts();
    assert_eq!(retained, material);
    let mechanical = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let system = fs_phs::transformer(mechanical, network, 0, 0, area).unwrap();
    assert_eq!((system.state_dim(), system.port_dim()), (5, 0));
    let mut state = [
        displacement,
        mass * velocity,
        q,
        initial_cavity[1],
        initial_cavity[2],
    ];
    let effort = system.effort(&state);
    let mut derivative = [0.0; 5];
    system
        .vector_field_from_effort_into(&effort, &[], &mut derivative)
        .unwrap();
    let pressure = observed.pressure_increment_pa;
    near(derivative[0], velocity, velocity, 3e-14);
    near(
        derivative[1],
        area * pressure - stiffness * displacement - damping * velocity,
        (area * pressure).abs() + stiffness * displacement + damping * velocity,
        3e-14,
    );
    near(
        derivative[2],
        -area * velocity - flow,
        area * velocity + flow,
        3e-14,
    );
    near(
        derivative[3],
        -pressure + resistance * flow,
        pressure.abs(),
        3e-14,
    );
    near(
        derivative[4],
        (root * q - state[4]) / material.arm().relaxation_time_s,
        root * q / material.arm().relaxation_time_s,
        3e-14,
    );
    let exchanged_power = (area * pressure) * velocity + pressure * (-area * velocity);
    near(
        exchanged_power,
        0.0,
        (area * pressure * velocity).abs(),
        3e-14,
    );
    let native_power: f64 = effort.iter().zip(derivative).map(|(a, b)| a * b).sum();
    let physical_loss = damping * velocity * velocity
        + observed.neck_dissipation_w
        + observed.compression_dissipation_w;
    near(native_power, -physical_loss, physical_loss, 2e-12);
    let initial_energy = system.hamiltonian(&state);
    near(
        initial_energy,
        0.5 * mass * velocity * velocity
            + 0.5 * stiffness * displacement * displacement
            + observed.stored_energy_j,
        initial_energy,
        3e-14,
    );
    let mut prepared =
        PreparedExactZoh::new(&system, material.arm().relaxation_time_s / 100.0, budget()).unwrap();
    let mut heat = 0.0;
    for _ in 0..100 {
        let mut next = [0.0; 5];
        let record = prepared.step_into(&state, &[], &mut next, &mut []).unwrap();
        assert_eq!(record.supplied, 0.0);
        heat += record.dissipated;
        state = next;
    }
    near(
        system.hamiltonian(&state) + heat,
        initial_energy,
        initial_energy,
        3e-9,
    );
    assert!(heat > 0.0 && system.hamiltonian(&state) < initial_energy);
    assert!((state[2] - q).abs() > 1e-4 * q);
    assert!(inertia > 0.0);
}

#[test]
fn g0_native_admission_and_small_density_perturbations_keep_original_source() {
    let material = material(20e6, 0.2);
    let geometry = geometry();
    for invalid_geometry in [
        HelmholtzCavity {
            volume_m3: 0.0,
            ..geometry
        },
        HelmholtzCavity {
            neck_radius_m: 0.0,
            ..geometry
        },
        HelmholtzCavity {
            neck_length_m: -0.1,
            ..geometry
        },
    ] {
        assert!(RelaxingWaterCavity::try_new(material, invalid_geometry, 0.0).is_err());
    }
    assert!(RelaxingWaterCavity::try_new(material, geometry, -1.0).is_err());
    let cavity = RelaxingWaterCavity::try_new(material, geometry, 0.0).unwrap();
    assert!(cavity.prepare(-1.0, budget()).is_err());
    assert!(
        cavity
            .prepare(
                0.001,
                LinearSplitBudget {
                    max_modes: 2,
                    ..budget()
                }
            )
            .is_err()
    );
    assert!(cavity.initial_state(f64::NAN, 0.0, None).is_err());
    assert!(cavity.initial_state(0.0, f64::INFINITY, None).is_err());
    assert!(cavity.initial_state(0.0, 0.0, Some(f64::NAN)).is_err());
    assert!(cavity.observe(&[0.0, 0.0]).is_err());
    assert!(cavity.observe(&[-geometry.volume_m3, 0.0, 0.0]).is_err());
    let tiny = cavity
        .initial_state(geometry.volume_m3 * 1e-20, 0.0, None)
        .unwrap();
    let observed = cavity.observe(&tiny).unwrap();
    assert!(observed.excess_mass_kg > 0.0 && observed.density_increment_kg_m3 > 0.0);
    assert_eq!(
        observed.linearized_density_kg_m3.to_bits(),
        material.source().bulk_density_kg_m3().to_bits()
    );
    let frozen = cavity
        .initial_state(geometry.volume_m3 * 1e-7, 0.0, Some(0.0))
        .unwrap();
    let record = cavity
        .prepare(0.0, budget())
        .unwrap()
        .step(&frozen, &[0.0])
        .unwrap();
    assert_eq!(record.x, frozen);
    assert_eq!(record.supplied, 0.0);
    assert_eq!(record.dissipated, 0.0);
    assert_eq!(cavity.material(), &material);
}
