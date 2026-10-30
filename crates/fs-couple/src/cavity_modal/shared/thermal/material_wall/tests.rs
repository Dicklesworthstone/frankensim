//! Synthetic source data test conjugate heat equations, not measured materials.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::gas::GasSpec;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_material::thermal::transport::ThermalTransportSourceQuery;
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::SpeciesId;

fn wall(rho: f64, cp: f64, k: f64) -> PreparedThermalTransport {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let keys = [
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0])),
    ];
    let mut claims = ClaimSet::new();
    for (key, value) in keys.iter().zip([rho, cp, k]) {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value,
                    dims: key.quantity().dims(),
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic fixed-geometry radial thermal source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
            phase: "solid".into(),
            process: "authored numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = QueryPoint::new()
        .with_quantity("T", temperature, 280.0)
        .unwrap();
    let upper = QueryPoint::new()
        .with_quantity("T", temperature, 360.0)
        .unwrap();
    PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        },
    )
    .unwrap()
}

fn source(gas: GasState) -> HelmholtzModalCavity {
    HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: 1e-7,
            neck_radius_m: 1e-4,
            neck_length_m: 0.1,
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
        },
        gas,
        CavityNeckLoss::Resistance(0.0),
    )
    .unwrap()
}
fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn molecular(gas: &GasState) -> RelaxationState {
    let species = SpeciesId::new("synthetic-rotor").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 101325.0 * 0.025,
                t_ref_k: 293.15,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 400.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas.clone(), vec![(species, 1.0)]).unwrap())
    .unwrap()
}
fn options() -> CylindricalMaterialWallOptions {
    CylindricalMaterialWallOptions {
        radius_m: 0.001,
        thickness_m: 0.0004,
        gas_cells: 1,
        wall_cells: 1,
        band_rad_s: [0.01, 10000.0],
        maximum_coordinates: 128,
    }
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}
fn complex_close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}

// Independent two-temperature Fourier circuit from radial resistances, before
// capacity normalization or pressure elimination. The molecular caloric response
// remains frequency dependent in this direct physical reference.
fn two_cell_parameters(
    source: &HelmholtzModalCavity,
    options: CylindricalMaterialWallOptions,
    wall: ThermalTransportSample,
) -> (f64, f64, f64) {
    let a = options.radius_m;
    let b = a + options.thickness_m;
    let length = source.geometry().volume_m3 / (core::f64::consts::PI * a * a);
    let rg = a / det::sqrt(2.0);
    let rw = det::sqrt((a * a + b * b) / 2.0);
    let g = 2.0 * core::f64::consts::PI * length
        / (det::ln(a / rg) / source.gas().thermal_conductivity
            + det::ln(rw / a) / wall.conductivity_w_m_k());
    let exterior =
        2.0 * core::f64::consts::PI * length * wall.conductivity_w_m_k() / det::ln(b / rw);
    let capacity =
        wall.volumetric_heat_capacity_j_m3_k() * core::f64::consts::PI * length * (b * b - a * a);
    (g, exterior, capacity)
}
fn two_cell_bulk(
    source: &HelmholtzModalCavity,
    wall: ThermalTransportSample,
    options: CylindricalMaterialWallOptions,
    relaxation: Option<&RelaxationState>,
    omega: f64,
) -> C64 {
    let gas = source.gas();
    let volume = source.geometry().volume_m3;
    let (g, exterior, cw) = two_cell_parameters(source, options, wall);
    let s = C64::new(0.0, omega);
    let cp = relaxation
        .map_or(Ok(C64::from_re(gas.specific_heat_cp)), |law| {
            law.caloric_response(omega).map(|c| c.cp_j_kg_k.conj())
        })
        .unwrap();
    let gas_capacity = cp.scale(gas.density * volume);
    let heat = C64::from_re(g) - C64::from_re(g * g) / (s.scale(cw) + C64::from_re(g + exterior));
    let temperature_per_pressure = s.scale(volume) / (s * gas_capacity + heat);
    (C64::ONE - temperature_per_pressure.scale(gas.pressure / gas.temperature))
        .recip()
        .scale(gas.pressure)
}

#[test]
fn wall_and_gas_sources_drive_one_caloric_law_against_direct_heat_equations() {
    let gate = CancelGate::new_clock_free();
    let mut heavy = GasSpec::dry_air_ussa1976();
    heavy.molar_mass = 0.06;
    heavy.gamma = 1.25;
    heavy.conductivity = fs_material::gas::ConductivityModel::Eucken;
    for gas in [
        air(),
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202650.0).unwrap(),
        GasState::try_new(&heavy, 293.15, 101325.0).unwrap(),
    ] {
        let source = source(gas);
        for material in [wall(10.0, 50.0, 0.03), wall(40.0, 120.0, 0.15)] {
            for enabled in [false, true] {
                let molecular = enabled.then(|| molecular(source.gas()));
                let law = CylindricalMaterialWallCompression::new(
                    &source,
                    options(),
                    material.clone(),
                    molecular.clone(),
                    &gate,
                )
                .unwrap();
                assert_eq!(law.wall().identity(), material.identity());
                assert_eq!(law.wall_sample().temperature_k(), source.gas().temperature);
                assert_eq!(law.thermal_coordinates(), if enabled { 3 } else { 2 });
                let scale = source.neck_area_m2().powi(2) / source.geometry().volume_m3;
                for omega in [0.01, 0.3, 4.0, 40.0, 200.0, 1000.0, 10000.0] {
                    let reference = two_cell_bulk(
                        &source,
                        law.wall_sample(),
                        options(),
                        molecular.as_ref(),
                        omega,
                    )
                    .scale(scale);
                    let actual = law.stiffness_at(omega).unwrap();
                    complex_close(actual, reference, 1e-11);
                    close(actual.im, reference.im, 1e-18, 2e-10);
                    assert!(actual.re > 0.0 && actual.im > 0.0);
                }
                let cv_inf = molecular.as_ref().map_or(
                    source.gas().specific_heat_cp - source.gas().specific_gas_constant,
                    RelaxationState::frozen_cv_j_kg_k,
                );
                close(
                    law.instantaneous_stiffness_n_m(),
                    scale
                        * source.gas().pressure
                        * (1.0 + source.gas().specific_gas_constant / cv_inf),
                    1e-18,
                    2e-12,
                );
                close(
                    law.equilibrium_stiffness_n_m(),
                    scale * source.gas().pressure,
                    1e-18,
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn radial_material_wall_converges_and_recovers_isothermal_and_insulating_limits() {
    let gate = CancelGate::new_clock_free();
    let source = source(air());
    let gas = source.gas();
    let base = options();
    let conductive = wall(1000.0, 100.0, 100.0);
    let mut previous = f64::INFINITY;
    for gas_cells in [8, 16, 32] {
        let law = CylindricalMaterialWallCompression::new(
            &source,
            CylindricalMaterialWallOptions {
                gas_cells,
                wall_cells: 2,
                thickness_m: 0.0001,
                ..base
            },
            conductive.clone(),
            None,
            &gate,
        )
        .unwrap();
        let mut error = 0.0_f64;
        for omega in [10.0, 40.0, 150.0, 500.0, 1000.0] {
            let wave = fs_duct::segment_wave(gas, base.radius_m, omega, LossModel::Bessel).unwrap();
            let compliance = (wave.wavenumber / wave.characteristic_impedance).scale(
                gas.density * gas.sound_speed.powi(2)
                    / (omega * core::f64::consts::PI * base.radius_m.powi(2)),
            );
            let reference = C64::from_re(source.compression_stiffness_n_m()) / compliance.conj();
            error =
                error.max((law.stiffness_at(omega).unwrap() - reference).abs() / reference.abs());
        }
        assert!(
            error < previous,
            "gas cells={gas_cells}, error={error:e}, previous={previous:e}"
        );
        previous = error;
    }
    assert!(previous < 0.003, "Bessel limit error={previous:e}");

    // Independent wall refinement with the gas grid held fixed. Neither axis
    // is hidden inside a fitted heat-transfer coefficient.
    let material = wall(10.0, 50.0, 0.003);
    let prepare = |wall_cells| {
        CylindricalMaterialWallCompression::new(
            &source,
            CylindricalMaterialWallOptions {
                gas_cells: 24,
                wall_cells,
                ..base
            },
            material.clone(),
            None,
            &gate,
        )
        .unwrap()
    };
    let reference = prepare(48);
    let mut previous = f64::INFINITY;
    for wall_cells in [2, 4, 8] {
        let law = prepare(wall_cells);
        let mut error = 0.0_f64;
        for omega in [3.0, 15.0, 80.0, 300.0] {
            let expected = reference.stiffness_at(omega).unwrap();
            error = error.max((law.stiffness_at(omega).unwrap() - expected).abs() / expected.abs());
        }
        assert!(
            error < previous,
            "wall cells={wall_cells}, error={error:e}, previous={previous:e}"
        );
        previous = error;
    }
    assert!(previous < 0.001);

    // A vanishing radial conductance suppresses wall heat escape at finite
    // frequency. Its zero-frequency limit remains isothermal, as it must.
    let insulated = CylindricalMaterialWallCompression::new(
        &source,
        CylindricalMaterialWallOptions {
            gas_cells: 8,
            wall_cells: 2,
            ..base
        },
        wall(10.0, 50.0, 1e-8),
        None,
        &gate,
    )
    .unwrap();
    complex_close(
        insulated.stiffness_at(100.0).unwrap(),
        C64::from_re(source.compression_stiffness_n_m()),
        2e-5,
    );
    close(
        insulated.equilibrium_stiffness_n_m(),
        source.compression_stiffness_n_m() / gas.gamma,
        1e-18,
        2e-12,
    );

    let make = |thickness_m, cp| {
        CylindricalMaterialWallCompression::new(
            &source,
            CylindricalMaterialWallOptions {
                thickness_m,
                gas_cells: 4,
                wall_cells: 4,
                ..base
            },
            wall(10.0, cp, 0.003),
            None,
            &gate,
        )
        .unwrap()
    };
    let thin = make(0.0001, 50.0);
    let thick = make(0.001, 50.0);
    let massive = make(0.001, 500.0);
    let at = |law: &CylindricalMaterialWallCompression| law.stiffness_at(40.0).unwrap();
    assert!((at(&thin) - at(&thick)).abs() > 0.005 * at(&thin).abs());
    assert!((at(&massive) - at(&thick)).abs() > 0.005 * at(&thick).abs());
}

fn coupled(wall: PreparedThermalTransport) -> SharedCavity {
    let source = source(air());
    let molecular = molecular(source.gas());
    let rate = 48000;
    let surface = CavitySurface {
        model: ModalAcousticTimeModel::try_new(
            rate,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 500.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap(),
        outward_volume_shapes: vec![1e-4],
    };
    let config = ModalCouplingConfig {
        max_modes: 8,
        max_connections: 6,
        max_setup_terms: 20000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8,
    };
    SharedCavity::new_with_cylindrical_material_wall(
        source,
        vec![surface],
        rate,
        ModalAcousticState {
            displacement_m_sqrt_kg: 2e-7,
            velocity_m_sqrt_kg_per_s: 0.0,
        },
        ModalAcousticTimeBudget::audible_reference(),
        config,
        None,
        options(),
        wall,
        Some(molecular),
        None,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn wall_temperatures_molecular_heat_and_pressure_share_the_original_runtime_energy() {
    let gate = CancelGate::new_clock_free();
    let mut actual = coupled(wall(10.0, 50.0, 0.03));
    let source = actual.source().clone();
    let gas = source.gas();
    let material = actual.cavity_material_wall_compression().unwrap();
    assert!(actual.cavity_thermal_compression().is_none());
    assert!(actual.cavity_molecular_compression().is_none());
    assert_eq!(
        actual.network().material_memory().len(),
        material.arms().len()
    );
    let law = material.molecular_relaxation().unwrap();
    let process = &law.processes()[0];
    let tau = process.relaxation_time_s;
    let volume = source.geometry().volume_m3;
    let cv = gas.density * volume * law.frozen_cv_j_kg_k();
    let cr = gas.density * volume * process.heat_capacity_j_kg_k;
    let (g, exterior, cw) = two_cell_parameters(&source, options(), material.wall_sample());
    let inertia = source.neck_mass_kg() / source.neck_area_m2().powi(2);
    let shapes = [1e-4, source.neck_volume_shape_m2_per_sqrt_kg()];
    for omega in [50.0, 300.0, 800.0] {
        let bulk =
            two_cell_bulk(&source, material.wall_sample(), options(), Some(law), omega).conj();
        let gain = bulk.scale(volume.recip());
        let a = C64::from_re(500.0_f64.powi(2) - omega * omega) + gain.scale(shapes[0] * shapes[0]);
        let b = gain.scale(shapes[0] * shapes[1]);
        let d = C64::from_re(-omega * omega) + gain.scale(shapes[1] * shapes[1]);
        let force = C64::new(1e-4, 2e-5);
        let determinant = a * d - b * b;
        let expected = [force * d / determinant, -force * b / determinant];
        let response = actual
            .harmonic_response(omega, &[force], C64::ZERO, 1_000_000, &gate)
            .unwrap();
        for (got, want) in response.network.displacement_m_sqrt_kg.iter().zip(expected) {
            complex_close(*got, want, 2e-9);
        }
        let delta_v = expected[0].scale(shapes[0]) + expected[1].scale(shapes[1]);
        complex_close(response.pressure_increment_pa, -gain * delta_v, 2e-9);
        let heat = -0.5 * omega * bulk.im / volume * delta_v.norm_sq();
        close(response.network.material_dissipation_w, heat, 1e-22, 2e-9);
        close(response.network.power_residual_w, 0.0, heat * 2e-9, 0.0);
    }
    // deltaV, neck flow, surface q/v, gas temperature, reservoir temperature,
    // wall temperature. No production spectral coordinates enter this oracle.
    let mut x = [0.0; 7];
    x[0] = actual.observe().unwrap().volume_change_m3;
    let pressure =
        |x: &[f64; 7]| -gas.pressure * x[0] / volume + gas.pressure * x[4] / gas.temperature;
    let derivative = |x: &[f64; 7], incident: f64| {
        let p = pressure(x);
        let flow = x[1] + 1e-4 * x[3];
        [
            flow,
            (p - incident) / inertia,
            x[3],
            -500.0_f64.powi(2) * x[2] + 1e-4 * p,
            (-gas.pressure * flow - g * (x[4] - x[6]) - cr / tau * (x[4] - x[5])) / cv,
            (x[4] - x[5]) / tau,
            (g * (x[4] - x[6]) - exterior * x[6]) / cw,
        ]
    };
    let initial = actual.network().total_energy_j().unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for tick in 0..160 {
        let incident = if tick < 20 { 0.05 } else { 0.0 };
        let dt = 1.0 / (48000.0 * 16.0);
        for _ in 0..16 {
            let stage =
                |k: &[f64; 7], factor: f64| std::array::from_fn(|i| x[i] + factor * dt * k[i]);
            let a = derivative(&x, incident);
            let b = derivative(&stage(&a, 0.5), incident);
            let c = derivative(&stage(&b, 0.5), incident);
            let d = derivative(&stage(&c, 1.0), incident);
            for i in 0..7 {
                x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
        let frame = actual.step(&[0.0], incident, &gate).unwrap().clone();
        let observed = actual.observe().unwrap();
        close(observed.pressure_increment_pa, pressure(&x), 1e-8, 1e-8);
        close(observed.neck_outflow_m3_s, x[1], 1e-18, 1e-8);
        close(observed.volume_change_m3, x[0], 1e-21, 1e-8);
        let availability = 0.5 * gas.pressure * x[0] * x[0] / volume
            + 0.5 * (cv * x[4] * x[4] + cr * x[5] * x[5] + cw * x[6] * x[6]) / gas.temperature;
        close(observed.compression_energy_j, availability, 1e-19, 2e-8);
        let total = availability
            + 0.5 * inertia * x[1] * x[1]
            + 0.5 * (x[3] * x[3] + 500.0_f64.powi(2) * x[2] * x[2]);
        close(frame.total_energy_j(), total, 1e-19, 2e-8);
        work += frame.external_work_j;
        loss += frame.total_dissipation_j();
        assert!(frame.material_dissipation_j >= 0.0);
    }
    assert!(loss > 0.0);
    close(
        actual.network().total_energy_j().unwrap() + loss,
        initial + work,
        1e-19,
        2e-8,
    );

    let before = actual.observe().unwrap();
    let memory = actual.network().material_memory().to_vec();
    let states: Vec<_> = actual
        .network()
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect();
    assert!(actual.step(&[0.0], f64::NAN, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        actual.step(&[0.0], 0.0, &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert_eq!(actual.observe().unwrap(), before);
    assert_eq!(actual.network().material_memory(), memory);
    assert_eq!(
        actual
            .network()
            .components()
            .iter()
            .map(|m| m.states().to_vec())
            .collect::<Vec<_>>(),
        states
    );
}

#[test]
fn finite_wall_admission_refuses_stale_sources_and_unbounded_or_unresolved_geometry() {
    let source = source(air());
    let gate = CancelGate::new_clock_free();
    let material = wall(10.0, 50.0, 0.03);
    let stale = GasState::try_new(&GasSpec::dry_air_ussa1976(), 303.15, 101325.0).unwrap();
    assert!(
        CylindricalMaterialWallCompression::new(
            &source,
            options(),
            material.clone(),
            Some(molecular(&stale)),
            &gate
        )
        .is_err()
    );
    for invalid in [
        CylindricalMaterialWallOptions {
            thickness_m: 0.0,
            ..options()
        },
        CylindricalMaterialWallOptions {
            gas_cells: 0,
            ..options()
        },
        CylindricalMaterialWallOptions {
            gas_cells: 129,
            ..options()
        },
        CylindricalMaterialWallOptions {
            maximum_coordinates: 1,
            ..options()
        },
    ] {
        assert!(
            CylindricalMaterialWallCompression::new(
                &source,
                invalid,
                material.clone(),
                None,
                &gate
            )
            .is_err()
        );
    }
    assert!(
        CylindricalMaterialWallCompression::new(
            &source,
            CylindricalMaterialWallOptions {
                maximum_coordinates: 2,
                ..options()
            },
            material.clone(),
            Some(molecular(source.gas())),
            &gate
        )
        .is_err()
    );
    let hot =
        self::source(GasState::try_new(&GasSpec::dry_air_ussa1976(), 370.0, 101325.0).unwrap());
    assert!(
        CylindricalMaterialWallCompression::new(&hot, options(), material.clone(), None, &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        CylindricalMaterialWallCompression::new(&source, options(), material, None, &gate),
        Err(ModalCouplingError::Cancelled)
    ));
}
