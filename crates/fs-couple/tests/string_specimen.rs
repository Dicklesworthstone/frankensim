//! G0/G1/G3: material resolution → circular specimen → existing pressure solver.
//! Synthetic material data check the implementation, not measured wire fidelity.

use fs_conduction::lumped::LumpedThermalEnvironment as ThermalStringEnvironment;
use fs_couple::acoustic_realize::{
    AcousticRealizeError, LinearMaterialStringRuntime, ThermalMaterialStringRuntime,
    ThermalStringError, realize_assembly, string_mode_omega,
};
use fs_couple::string_specimen::{
    BENDING_RELAXATION_TIME_PROPERTY, EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, RELAXING_BENDING_MODULUS_PROPERTY,
    StringGeometryConstraint, StringPrestress, with_uniform_circular_material_and_constraints,
    with_uniform_circular_material_and_prestress, with_uniform_circular_material_state,
    with_uniform_circular_thermal_extension,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, IntegratedIsotropicThermalExpansion,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_qty::{
    Density, Dims, DynViscosity, Pressure, QuantitySpec, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_scenario::acoustic::BendingRelaxationProperties;
use fs_scenario::{
    AcousticAssembly, AmbientGas, Listener, Pluck, PrestressedString, RayleighParams,
};

fn state(name: &str, properties: &[(&str, Dims, f64)]) -> ResolvedMaterialStatePoint {
    let quantities: Vec<_> = properties
        .iter()
        .map(|&(key, dims, value)| (key, QuantitySpec::dimensional(dims), value))
        .collect();
    state_with_quantities(name, &quantities)
}

fn state_with_quantities(
    name: &str,
    properties: &[(&str, QuantitySpec, f64)],
) -> ResolvedMaterialStatePoint {
    state_with_domain(
        name,
        properties,
        ValidityDomain::unconstrained().with("T", 290.0, 300.0),
        QueryPoint::new().with("T", 293.15).unwrap(),
        None,
    )
}

fn state_with_domain(
    name: &str,
    properties: &[(&str, QuantitySpec, f64)],
    validity: ValidityDomain,
    point: QueryPoint,
    curve_property: Option<&str>,
) -> ResolvedMaterialStatePoint {
    state_with_property_domains(
        name,
        properties,
        |_| validity.clone(),
        point,
        curve_property,
    )
}

fn state_with_property_domains(
    name: &str,
    properties: &[(&str, QuantitySpec, f64)],
    validity: impl Fn(&str) -> ValidityDomain,
    point: QueryPoint,
    curve_property: Option<&str>,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(key, quantity, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(key, quantity),
                value: if curve_property == Some(key) {
                    PropertyValue::Curve {
                        abscissa: "omega".into(),
                        abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
                        knots: vec![(1.0, value), (20_000.0, 2.0 * value)],
                        dims: quantity.dims(),
                    }
                } else {
                    PropertyValue::Scalar {
                        value,
                        dims: quantity.dims(),
                    }
                },
                validity: validity(key),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve_property == Some(key) {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic uniform string".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key,
                quantity,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn elastic(name: &str, rho: f64, young: f64) -> ResolvedMaterialStatePoint {
    state(
        name,
        &[
            ("density", Density::DIMS, rho),
            ("young_modulus", Pressure::DIMS, young),
        ],
    )
}

fn thermal_state(
    temperature_k: f64,
) -> (
    ResolvedMaterialStatePoint,
    IntegratedIsotropicThermalExpansion,
) {
    thermal_state_with_density(
        temperature_k,
        PropertyValue::Scalar {
            value: 7800.0,
            dims: Density::DIMS,
        },
    )
}

fn thermal_state_with_density(
    temperature_k: f64,
    density: PropertyValue,
) -> (
    ResolvedMaterialStatePoint,
    IntegratedIsotropicThermalExpansion,
) {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, density),
        (
            "young_modulus",
            Pressure::DIMS,
            PropertyValue::Scalar {
                value: 200e9,
                dims: Pressure::DIMS,
            },
        ),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                knots: vec![
                    (250.0, 5e-6),
                    (300.0, 10e-6),
                    (350.0, 20e-6),
                    (400.0, 30e-6),
                    (450.0, 40e-6),
                ],
                dims: INVERSE_TEMPERATURE_DIMS,
            },
        ),
    ] {
        let interpolation = if matches!(&value, PropertyValue::Curve { .. }) {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value,
                validity: ValidityDomain::unconstrained().with("T", 250.0, 450.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic thermal string".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal wire".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", temperature_k).unwrap();
    let requirements = [
        ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "young_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    let state = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let expansion = integrate_isotropic_thermal_expansion(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (state, expansion)
}

#[test]
fn g1_thermal_eigenstrain_changes_fixed_support_tension_with_retained_receipts() {
    let mut tensions = Vec::new();
    for (temperature, thermal_strain) in [(300.0, 0.0), (400.0, 0.002), (250.0, -0.000375)] {
        let (state, expansion) = thermal_state(temperature);
        let specimen = with_uniform_circular_thermal_extension(
            PrestressedString {
                length_m: 0.502,
                ..template()
            },
            StringGeometryConstraint::FixedMass(0.001),
            &state,
            &expansion,
            0.5,
            0.01,
        )
        .unwrap();
        assert_eq!(specimen.material(), &state);
        assert_eq!(specimen.thermal_expansion(), Some(&expansion));
        assert_eq!(
            specimen.material().properties().len(),
            2,
            "no unrelated property requirements"
        );
        close(specimen.mass_kg(), 0.001);
        let area = 0.001 / (7800.0 * 0.502);
        let expected = 200e9 * area * ((0.502 - 0.5) / 0.5 - thermal_strain);
        close(specimen.string().tension_n, expected);
        tensions.push(expected);
        let force_controlled = with_uniform_circular_material_and_constraints(
            specimen.string(),
            StringGeometryConstraint::FixedMass(0.001),
            &state,
            StringPrestress::FixedTension(20.0),
        )
        .unwrap();
        close(force_controlled.string().tension_n, 20.0);
        assert!(force_controlled.thermal_expansion().is_none());
    }
    assert!(tensions[1] < tensions[0] && tensions[2] > tensions[0]);
}

#[test]
fn g3_thermal_expansion_changes_actual_plucked_pressure_pitch() {
    let mut pitches = Vec::new();
    for temperature in [300.0, 400.0] {
        let (state, expansion) = thermal_state(temperature);
        let specimen = with_uniform_circular_thermal_extension(
            PrestressedString {
                length_m: 0.502,
                ..template()
            },
            StringGeometryConstraint::FixedMass(0.001),
            &state,
            &expansion,
            0.5,
            0.01,
        )
        .unwrap();
        let string = specimen.string();
        let expected = string_mode_omega(&string, 1) / core::f64::consts::TAU;
        let pitch = measured_hz(&pressure(string));
        assert!(
            (pitch / expected - 1.0).abs() < 0.01,
            "physical pressure pitch {pitch} vs mode {expected}"
        );
        pitches.push(pitch);
    }
    assert!(
        (pitches[1] / pitches[0] - 0.5_f64.sqrt()).abs() < 0.01,
        "heating reduces elastic prestress: {pitches:?}"
    );
    eprintln!("G3 thermal string 300/400 K pressure pitches: {pitches:?} Hz");
}

#[test]
fn g0_thermal_string_refuses_mismatched_state_slack_and_excess_strain() {
    let (hot, expansion) = thermal_state(400.0);
    let (cold, _) = thermal_state(300.0);
    let other = state_with_domain(
        "another wire",
        &[
            ("density", QuantitySpec::dimensional(Density::DIMS), 7800.0),
            (
                "young_modulus",
                QuantitySpec::dimensional(Pressure::DIMS),
                200e9,
            ),
        ],
        ValidityDomain::unconstrained().with("T", 250.0, 450.0),
        QueryPoint::new().with("T", 400.0).unwrap(),
        None,
    );
    let string = PrestressedString {
        length_m: 0.502,
        ..template()
    };
    for state in [&cold, &other] {
        assert!(
            with_uniform_circular_thermal_extension(
                string.clone(),
                StringGeometryConstraint::FixedMass(0.001),
                state,
                &expansion,
                0.5,
                0.01
            )
            .is_err()
        );
    }
    for (length, reference, limit, moving) in [
        (0.5005, 0.5, 0.01, false), // expansion exceeds initial stretch: compression
        (0.51, 0.5, 0.01, false),   // total reference strain too large
        (0.502, 0.5, 0.001, false), // thermal strain too large
        (0.502, 0.0, 0.01, false),
        (0.502, 0.5, f64::NAN, false),
        (0.502, 0.5, 0.01, true),
    ] {
        assert!(
            with_uniform_circular_thermal_extension(
                PrestressedString {
                    length_m: length,
                    moving_end: moving,
                    ..string.clone()
                },
                StringGeometryConstraint::FixedMass(0.001),
                &hot,
                &expansion,
                reference,
                limit
            )
            .is_err()
        );
    }
}

fn template() -> PrestressedString {
    PrestressedString {
        length_m: 0.5,
        tension_n: 20.0,
        // These independent legacy inputs must be replaced by the binding.
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 1.0,
        bending_stiffness_n_m2: 1.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        width_m: 1.0,
        n_modes: 1,
        damping_ratio: 0.0,
        rayleigh: Some(RayleighParams {
            alpha_per_s: 2.0,
            beta_s: 0.0,
        }),
        polarization_detune: 0.0,
        moving_end: false,
    }
}

fn with_string_cx<R>(f: impl FnOnce(&fs_exec::Cx<'_>, &fs_exec::CancelGate) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 5,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn incremental_specimen(
    temperature_k: f64,
    n_modes: usize,
) -> fs_couple::string_specimen::ResolvedStringSpecimen {
    let (state, expansion) = thermal_state(temperature_k);
    with_uniform_circular_thermal_extension(
        PrestressedString {
            length_m: 0.502,
            n_modes,
            ..template()
        },
        StringGeometryConstraint::FixedMass(0.001),
        &state,
        &expansion,
        0.5,
        0.01,
    )
    .unwrap()
}

fn incremental_ambient() -> AmbientGas {
    AmbientGas {
        temperature_k: 300.0,
        pressure_pa: 101_325.0,
        relative_humidity: 0.0,
    }
}

#[test]
fn g1_incremental_material_hot_rebind_preserves_motion_and_accounts_parameter_work() {
    with_string_cx(|cx, _| {
        let cold = incremental_specimen(300.0, 1);
        let hot = incremental_specimen(400.0, 1);
        let mut runtime = LinearMaterialStringRuntime::try_new(
            cx,
            cold,
            Some(Pluck {
                station_frac: 0.4,
                height_m: 1e-5,
            }),
            incremental_ambient(),
            1.0,
            48_000,
        )
        .unwrap();
        for _ in 0..8 {
            runtime.step(cx, &[0.0]).unwrap();
        }
        let mut unchanged = runtime.clone();
        let before = runtime.states()[0];
        let old_omega = runtime.modes()[0].angular_frequency_rad_s;
        let update = runtime.rebind(cx, hot.clone(), 1.0).unwrap();
        assert_eq!(runtime.states(), &[before]);
        assert_eq!(runtime.accepted_samples(), 8);
        assert_eq!(runtime.epoch(), 9);
        assert_eq!(runtime.specimen().material(), hot.material());
        let omega = runtime.modes()[0].angular_frequency_rad_s;
        assert!(omega < old_omega);
        let expected_work =
            0.5 * (omega.powi(2) - old_omega.powi(2)) * before.displacement_m_sqrt_kg.powi(2);
        close(update.parameter_work_j, expected_work);
        close(
            update.energy_after_j - update.energy_before_j,
            expected_work,
        );

        // Independent damped-oscillator solution, Rayleigh alpha=2 => decay=1/s.
        let dt: f64 = 1.0 / 48_000.0;
        let wd = (omega * omega - 1.0).sqrt();
        let q = before.displacement_m_sqrt_kg;
        let v = before.velocity_m_sqrt_kg_per_s;
        let expected_q = (-dt).exp() * (q * (wd * dt).cos() + (v + q) / wd * (wd * dt).sin());
        let expected_v =
            (-dt).exp() * (v * (wd * dt).cos() - (v + omega * omega * q) / wd * (wd * dt).sin());
        let frame = runtime.step(cx, &[0.0]).unwrap();
        close(runtime.states()[0].displacement_m_sqrt_kg, expected_q);
        close(runtime.states()[0].velocity_m_sqrt_kg_per_s, expected_v);
        let string = hot.string();
        let gas = fs_material::gas::GasState::try_new_moist_air(300.0, 101_325.0, 0.0).unwrap();
        let scale = (string.lin_density_kg_m * string.length_m / 2.0).sqrt();
        let area_weight = gas.density * string.width_m * string.length_m
            / (2.0 * core::f64::consts::PI.powi(2) * scale);
        let expected_pressure = area_weight * (-omega * omega * expected_q - 2.0 * expected_v);
        close(frame.acoustic.observer_pressure_pa, expected_pressure);
        let cold_frame = unchanged.step(cx, &[0.0]).unwrap();
        assert!(
            (frame.acoustic.observer_pressure_pa - cold_frame.acoustic.observer_pressure_pa).abs()
                > cold_frame.acoustic.observer_pressure_pa.abs() * 0.01
        );
        assert_eq!(runtime.epoch(), 10);
    });
}

#[test]
fn g3_incremental_fixed_mass_density_cycle_preserves_basis_and_changes_pressure() {
    // These are prescribed current-state density/strain data, not a closed
    // heat/Poisson-contraction model. The same immutable card supplies both states.
    let bind = |temperature_k| {
        let (state, expansion) = thermal_state_with_density(
            temperature_k,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                knots: vec![(300.0, 7800.0), (400.0, 7784.4)],
                dims: Density::DIMS,
            },
        );
        with_uniform_circular_thermal_extension(
            PrestressedString {
                length_m: 0.502,
                n_modes: 1,
                ..template()
            },
            StringGeometryConstraint::FixedMass(0.001),
            &state,
            &expansion,
            0.5,
            0.01,
        )
        .unwrap()
    };
    with_string_cx(|cx, _| {
        let cold = bind(300.0);
        let hot = bind(400.0);
        let mut runtime = LinearMaterialStringRuntime::try_new(
            cx,
            cold.clone(),
            Some(Pluck {
                station_frac: 0.4,
                height_m: 1e-5,
            }),
            incremental_ambient(),
            1.0,
            48_000,
        )
        .unwrap();
        for _ in 0..8 {
            runtime.step(cx, &[0.0]).unwrap();
        }
        let mut cold_continuation = runtime.clone();
        let before = runtime.states().to_vec();
        let cold_modes = runtime.modes().to_vec();
        let update = runtime.rebind(cx, hot.clone(), 1.0).unwrap();
        assert_eq!(runtime.states(), before);
        assert_eq!(runtime.accepted_samples(), 8);
        assert_eq!(runtime.epoch(), 9);
        assert!(update.parameter_work_j < 0.0);
        for specimen in [&cold, &hot] {
            assert_eq!(specimen.mass_kg().to_bits(), 0.001_f64.to_bits());
            assert_eq!(
                specimen.string().lin_density_kg_m.to_bits(),
                (0.001_f64 / 0.502).to_bits()
            );
        }
        assert_ne!(cold.specimen_identity(), hot.specimen_identity());
        close(
            hot.string().width_m / cold.string().width_m,
            (7800.0_f64 / 7784.4).sqrt(),
        );
        close(
            hot.string().bending_stiffness_n_m2 / cold.string().bending_stiffness_n_m2,
            (7800.0_f64 / 7784.4).powi(2),
        );

        let frame = runtime.step(cx, &[0.0]).unwrap();
        let mode = runtime.modes()[0];
        let state = runtime.states()[0];
        let acceleration = -mode.angular_frequency_rad_s.powi(2) * state.displacement_m_sqrt_kg
            - 2.0
                * mode.damping_ratio
                * mode.angular_frequency_rad_s
                * state.velocity_m_sqrt_kg_per_s;
        let gas = fs_material::gas::GasState::try_new_moist_air(300.0, 101_325.0, 0.0).unwrap();
        let area_weight = gas.density * hot.string().width_m * 0.502
            / (2.0 * core::f64::consts::PI.powi(2) * (0.001_f64 / 2.0).sqrt());
        close(
            frame.acoustic.observer_pressure_pa,
            area_weight * acceleration,
        );
        let cold_frame = cold_continuation.step(cx, &[0.0]).unwrap();
        assert!(
            (frame.acoustic.observer_pressure_pa - cold_frame.acoustic.observer_pressure_pa).abs()
                > 0.01 * cold_frame.acoustic.observer_pressure_pa.abs()
        );

        let hot_motion = runtime.states().to_vec();
        let reverse = runtime.rebind(cx, cold.clone(), 1.0).unwrap();
        assert!(reverse.parameter_work_j > 0.0);
        assert_eq!(runtime.states(), hot_motion);
        assert_eq!(runtime.modes(), cold_modes);
        assert_eq!(runtime.specimen(), &cold);
        assert_eq!(runtime.accepted_samples(), 9);
        assert_eq!(runtime.epoch(), 11);
    });
}

#[test]
fn g0_fixed_mass_identity_distinguishes_masses_with_the_same_rounded_radius() {
    let material = elastic("mass-rounding", 7800.0, 200e9);
    let bind = |mass| {
        with_uniform_circular_material_and_constraints(
            PrestressedString {
                length_m: 0.502,
                ..template()
            },
            StringGeometryConstraint::FixedMass(mass),
            &material,
            StringPrestress::FixedTension(20.0),
        )
        .unwrap()
    };
    let first = bind(0.001);
    let second = bind(f64::from_bits(0.001_f64.to_bits() + 1));
    assert_eq!(
        first.string().width_m.to_bits(),
        second.string().width_m.to_bits()
    );
    assert_ne!(first.mass_kg().to_bits(), second.mass_kg().to_bits());
    assert_ne!(
        first.string().lin_density_kg_m.to_bits(),
        second.string().lin_density_kg_m.to_bits()
    );
    assert_ne!(first.specimen_identity(), second.specimen_identity());
}

#[test]
fn g4_incremental_material_refusals_preserve_epoch_and_replay_suffix() {
    with_string_cx(|cx, _| {
        let mut runtime = LinearMaterialStringRuntime::try_new(
            cx,
            incremental_specimen(300.0, 2),
            Some(Pluck {
                station_frac: 0.4,
                height_m: 1e-5,
            }),
            incremental_ambient(),
            1.0,
            48_000,
        )
        .unwrap();
        for _ in 0..8 {
            runtime.step(cx, &[0.0, 0.0]).unwrap();
        }
        let mut reference = runtime.clone();
        let hot = incremental_specimen(400.0, 2);
        assert!(matches!(
            runtime.rebind(cx, hot.clone(), 0.0),
            Err(AcousticRealizeError::InvalidDescription {
                what: "material update exceeds its parameter-work budget"
            })
        ));
        assert!(
            runtime
                .rebind(cx, incremental_specimen(400.0, 1), 1.0)
                .is_err()
        );
        for (length_m, mass_kg) in [(0.503, 0.001), (0.502, 0.002)] {
            let changed = with_uniform_circular_material_and_constraints(
                PrestressedString {
                    length_m,
                    ..hot.string()
                },
                StringGeometryConstraint::FixedMass(mass_kg),
                hot.material(),
                StringPrestress::FixedTension(20.0),
            )
            .unwrap();
            assert!(runtime.rebind(cx, changed, 1.0).is_err());
        }
        assert!(runtime.step(cx, &[f64::NAN, 0.0]).is_err());
        with_string_cx(|cancelled, gate| {
            gate.request();
            assert_eq!(
                runtime.step(cancelled, &[0.0, 0.0]),
                Err(AcousticRealizeError::Cancelled)
            );
            assert_eq!(
                runtime.rebind(cancelled, hot.clone(), 1.0),
                Err(AcousticRealizeError::Cancelled)
            );
        });
        assert_eq!(runtime.states(), reference.states());
        assert_eq!(runtime.modes(), reference.modes());
        assert_eq!(
            runtime.specimen().material(),
            reference.specimen().material()
        );
        assert_eq!(runtime.epoch(), reference.epoch());
        assert_eq!(runtime.accepted_samples(), reference.accepted_samples());
        let before_acceleration: Vec<f64> = runtime
            .modes()
            .iter()
            .zip(runtime.states())
            .map(|(mode, state)| {
                -mode.angular_frequency_rad_s.powi(2) * state.displacement_m_sqrt_kg
                    - 2.0
                        * mode.damping_ratio
                        * mode.angular_frequency_rad_s
                        * state.velocity_m_sqrt_kg_per_s
            })
            .collect();
        assert_eq!(
            runtime.rebind(cx, hot.clone(), 1.0).unwrap(),
            reference.rebind(cx, hot, 1.0).unwrap()
        );
        for sample in 0..32 {
            let frame = runtime.step(cx, &[0.0, 0.0]).unwrap();
            assert_eq!(frame, reference.step(cx, &[0.0, 0.0]).unwrap());
            if sample == 0 {
                let string = runtime.specimen().string();
                let gas =
                    fs_material::gas::GasState::try_new_moist_air(300.0, 101_325.0, 0.0).unwrap();
                let pi = core::f64::consts::PI;
                let factor = gas.density * string.width_m
                    / (4.0 * pi * (string.lin_density_kg_m * string.length_m / 2.0).sqrt());
                let acceleration: Vec<f64> = runtime
                    .modes()
                    .iter()
                    .zip(runtime.states())
                    .map(|(mode, state)| {
                        -mode.angular_frequency_rad_s.powi(2) * state.displacement_m_sqrt_kg
                            - 2.0
                                * mode.damping_ratio
                                * mode.angular_frequency_rad_s
                                * state.velocity_m_sqrt_kg_per_s
                    })
                    .collect();
                // The even mode's jerk must span the material update, retaining
                // the last accepted cold acceleration rather than reinitializing it.
                let expected = factor
                    * (2.0 * string.length_m / pi * acceleration[0]
                        - string.length_m.powi(2) / (2.0 * pi * gas.sound_speed)
                            * (acceleration[1] - before_acceleration[1])
                            * 48_000.0);
                close(frame.acoustic.observer_pressure_pa, expected);
            }
        }
        assert_eq!(runtime.states(), reference.states());
    });
}

#[test]
fn g4_incremental_material_late_pressure_refusal_preserves_mechanics() {
    with_string_cx(|cx, _| {
        // Deliberately pathological observation geometry forces a late pressure
        // refusal while the candidate oscillator remains inside its state budget.
        let mut runtime = LinearMaterialStringRuntime::try_new(
            cx,
            incremental_specimen(300.0, 1),
            None,
            incremental_ambient(),
            1e-12,
            48_000,
        )
        .unwrap();
        let mut reference = runtime.clone();
        assert!(matches!(
            runtime.step(cx, &[1.0]),
            Err(AcousticRealizeError::InvalidDescription {
                what: "compact string observer pressure is nonfinite or exceeds its pressure budget"
            })
        ));
        assert_eq!(runtime.states(), reference.states());
        assert_eq!(runtime.epoch(), 0);
        assert_eq!(runtime.accepted_samples(), 0);
        assert_eq!(
            runtime.step(cx, &[0.0]).unwrap(),
            reference.step(cx, &[0.0]).unwrap()
        );
    });
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 2.0e-13,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_circular_mass_and_stiffness_match_diameter_formulas_and_beam_modes() {
    let material = elastic("specimen", 7800.0, 200.0e9);
    let resolved = with_uniform_circular_material_state(template(), 0.0005, &material).unwrap();
    let string = resolved.string();
    let diameter: f64 = 0.001;
    let area = core::f64::consts::PI * diameter.powi(2) / 4.0;
    let moment = core::f64::consts::PI * diameter.powi(4) / 64.0;
    close(resolved.area_m2(), area);
    close(resolved.second_moment_m4(), moment);
    close(resolved.mass_kg(), 7800.0 * area * 0.5);
    close(string.lin_density_kg_m, 7800.0 * area);
    close(string.axial_stiffness_n, 200.0e9 * area);
    close(string.bending_stiffness_n_m2, 200.0e9 * moment);
    assert_eq!(string.width_m.to_bits(), diameter.to_bits());
    assert_eq!(string.tension_n.to_bits(), template().tension_n.to_bits());
    assert_eq!(string.rayleigh, template().rayleigh);
    assert_eq!(resolved.material(), &material);
    for n in 1..=4 {
        let k = f64::from(n) * core::f64::consts::PI / 0.5;
        let omega = ((20.0 * k.powi(2) + 200.0e9 * moment * k.powi(4)) / (7800.0 * area)).sqrt();
        close(string_mode_omega(&string, n as usize), omega);
    }
}

#[test]
fn g3_geometric_scaling_and_material_replacement_move_all_derived_inputs() {
    let material = elastic("first", 1000.0, 2.0e9);
    let first = with_uniform_circular_material_state(template(), 0.0005, &material).unwrap();
    let twice = with_uniform_circular_material_state(
        PrestressedString {
            length_m: 1.0,
            ..template()
        },
        0.001,
        &material,
    )
    .unwrap();
    close(twice.mass_kg() / first.mass_kg(), 8.0);
    close(
        twice.string().lin_density_kg_m / first.string().lin_density_kg_m,
        4.0,
    );
    close(
        twice.string().axial_stiffness_n / first.string().axial_stiffness_n,
        4.0,
    );
    close(
        twice.string().bending_stiffness_n_m2 / first.string().bending_stiffness_n_m2,
        16.0,
    );
    let changed = with_uniform_circular_material_state(
        template(),
        0.0005,
        &elastic("replacement", 4000.0, 6.0e9),
    )
    .unwrap();
    close(changed.mass_kg() / first.mass_kg(), 4.0);
    close(
        changed.string().axial_stiffness_n / first.string().axial_stiffness_n,
        3.0,
    );
    close(
        changed.string().bending_stiffness_n_m2 / first.string().bending_stiffness_n_m2,
        3.0,
    );
    assert_ne!(changed.specimen_identity(), first.specimen_identity());
    let same_values = with_uniform_circular_material_state(
        template(),
        0.0005,
        &elastic("different-lot", 1000.0, 2.0e9),
    )
    .unwrap();
    assert_eq!(same_values.string(), first.string());
    assert_ne!(same_values.specimen_identity(), first.specimen_identity());
}

#[test]
fn g1_fixed_mass_derives_geometry_stiffness_and_independent_prestress() {
    let mass = 0.001;
    let length = template().length_m;
    let young = 200.0e9;
    let constraint = StringGeometryConstraint::FixedMass(mass);
    for density in [1000.0, 4000.0] {
        let material = elastic("fixed-mass", density, young);
        for prestress in [
            StringPrestress::FixedTension(20.0),
            StringPrestress::FixedExtension {
                stress_free_length_m: 0.498,
                linear_strain_limit: 0.005,
            },
            StringPrestress::FixedThermalExtension {
                reference_stress_free_length_m: 0.498,
                free_thermal_strain: 0.001,
                linear_strain_limit: 0.005,
            },
            StringPrestress::TargetFundamentalHz(160.0),
        ] {
            let specimen = with_uniform_circular_material_and_constraints(
                template(),
                constraint,
                &material,
                prestress,
            )
            .unwrap();
            let string = specimen.string();
            let area = mass / (density * length);
            let moment = area.powi(2) / (4.0 * core::f64::consts::PI);
            close(specimen.mass_kg(), mass);
            close(specimen.area_m2(), area);
            close(specimen.second_moment_m4(), moment);
            close(string.lin_density_kg_m, mass / length);
            close(string.axial_stiffness_n, young * area);
            close(string.bending_stiffness_n_m2, young * moment);
            close(string.width_m, (4.0 * area / core::f64::consts::PI).sqrt());
            let expected_tension = match prestress {
                StringPrestress::FixedTension(force) => force,
                StringPrestress::FixedExtension {
                    stress_free_length_m,
                    ..
                } => young * area * (length / stress_free_length_m - 1.0),
                StringPrestress::FixedThermalExtension {
                    reference_stress_free_length_m,
                    free_thermal_strain,
                    ..
                } => {
                    young
                        * area
                        * (length / reference_stress_free_length_m - 1.0 - free_thermal_strain)
                }
                StringPrestress::TargetFundamentalHz(hz) => {
                    4.0 * mass * length * hz.powi(2)
                        - young * moment * (core::f64::consts::PI / length).powi(2)
                }
            };
            close(string.tension_n, expected_tension);
            assert_eq!(specimen.geometry_constraint(), constraint);
            assert_eq!(specimen.prestress(), prestress);
            assert_eq!(specimen.material(), &material);
            assert_eq!(string.rayleigh, template().rayleigh);
            let same_geometry = with_uniform_circular_material_and_prestress(
                template(),
                string.width_m / 2.0,
                &material,
                prestress,
            )
            .unwrap();
            assert_eq!(same_geometry.string(), string);
            assert_eq!(
                same_geometry.specimen_identity(),
                specimen.specimen_identity()
            );
            assert_ne!(same_geometry.geometry_constraint(), constraint);
        }
    }
}

#[test]
fn g3_fixed_mass_material_swap_reaches_pressure_with_bending_change() {
    let bind = |density| {
        with_uniform_circular_material_and_constraints(
            template(),
            StringGeometryConstraint::FixedMass(0.001),
            &elastic("fixed-mass-pressure", density, 200.0e9),
            StringPrestress::FixedTension(20.0),
        )
        .unwrap()
    };
    let light = bind(1000.0);
    let heavy = bind(4000.0);
    close(heavy.mass_kg(), light.mass_kg());
    close(heavy.string().width_m / light.string().width_m, 0.5);
    close(
        heavy.string().axial_stiffness_n / light.string().axial_stiffness_n,
        0.25,
    );
    close(
        heavy.string().bending_stiffness_n_m2 / light.string().bending_stiffness_n_m2,
        1.0 / 16.0,
    );
    let mut pitches = Vec::new();
    for (specimen, density) in [(&light, 1000.0_f64), (&heavy, 4000.0_f64)] {
        // Eliminate A and I independently: f_n² = n² T/(4mL)
        // + n⁴ E pi m/(16 rho² L⁵). Fixed mass preserves only the tension term.
        let frequency = |n: f64| {
            (n.powi(2) * 20.0 / (4.0 * 0.001 * 0.5)
                + n.powi(4) * 200.0e9 * core::f64::consts::PI * 0.001
                    / (16.0 * density.powi(2) * 0.5_f64.powi(5)))
            .sqrt()
        };
        for n in 1..=4 {
            close(
                string_mode_omega(&specimen.string(), n) / core::f64::consts::TAU,
                frequency(n as f64),
            );
        }
        let waveform = pressure(specimen.string());
        assert!(waveform.iter().all(|p| p.is_finite()));
        assert_eq!(waveform, pressure(specimen.string()));
        let measured = measured_hz(&waveform);
        assert!((measured / frequency(1.0) - 1.0).abs() < 0.01);
        pitches.push(measured);
    }
    assert!(
        pitches[0] > 1.04 * pitches[1],
        "bending must change pitch at fixed mass"
    );
    eprintln!(
        "G3 fixed mass, density x4: pressure pitch {:.6} -> {:.6} Hz",
        pitches[0], pitches[1]
    );
}

fn pressure(string: PrestressedString) -> Vec<f64> {
    realize_assembly(&assembly(string)).unwrap().pressure_pa
}

fn assembly(string: PrestressedString) -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas::sea_level(),
        string: Some(string),
        duct: None,
        pluck: Some(Pluck {
            station_frac: 0.5,
            height_m: 1.0e-6,
        }),
        bow: None,
        blow: None,
        reed: None,
        soundboard: None,
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 8000,
        duration_s: 0.12,
    }
}

fn measured_hz(pressure: &[f64]) -> f64 {
    let crossings: Vec<_> = pressure
        .windows(2)
        .enumerate()
        .skip(160)
        .filter(|(_, p)| p[0] > 0.0 && p[1] <= 0.0)
        .map(|(i, p)| i as f64 + p[0] / (p[0] - p[1]))
        .collect();
    assert!(
        crossings.len() >= 5,
        "a live decaying pressure wave is required"
    );
    8000.0 * (crossings.len() - 1) as f64 / (crossings.last().unwrap() - crossings[0])
}

fn viscous_material(
    eta: f64,
    quantity: QuantitySpec,
    band: (f64, f64),
) -> ResolvedMaterialStatePoint {
    state_with_domain(
        "synthetic Kelvin-Voigt solid",
        &[
            ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
            (
                "young_modulus",
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0e9,
            ),
            (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, quantity, eta),
        ],
        ValidityDomain::unconstrained()
            .with("T", 290.0, 300.0)
            .with("omega", band.0, band.1),
        QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", band.0)
            .unwrap(),
        None,
    )
}

fn loss_template() -> PrestressedString {
    PrestressedString {
        rayleigh: None,
        ..template()
    }
}

fn heated_string_parts(
    cx: &fs_exec::Cx<'_>,
    rate: u32,
    elastic_curve: bool,
    viscosity_policy: InterpolationPolicy,
) -> (
    LinearMaterialStringRuntime,
    MaterialCard,
    fs_material::phase::EquilibriumEnthalpyPhaseCurve,
) {
    heated_string_parts_with_conductor(cx, rate, elastic_curve, viscosity_policy, false)
}

fn heated_string_parts_with_conductor(
    cx: &fs_exec::Cx<'_>,
    rate: u32,
    elastic_curve: bool,
    viscosity_policy: InterpolationPolicy,
    electrical: bool,
) -> (
    LinearMaterialStringRuntime,
    MaterialCard,
    fs_material::phase::EquilibriumEnthalpyPhaseCurve,
) {
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    let mut properties = vec![
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 2e9),
        (
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            1e7,
        ),
    ];
    if electrical {
        properties.push((
            fs_material::conductor::ELECTRICAL_RESISTIVITY_PROPERTY,
            fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS,
            1e-6,
        ));
    }
    for (name, dims, value) in properties {
        let is_curve = name == KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY
            || name == fs_material::conductor::ELECTRICAL_RESISTIVITY_PROPERTY
            || elastic_curve && name == "young_modulus";
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if is_curve {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(300.0, value), (340.0, 9.0 * value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 340.0)
                    .with("omega", 1.0, 20000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if is_curve {
                    viscosity_policy
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic dissipative thermal feedback".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal string".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 1.0)
        .unwrap();
    let state = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let specimen = with_uniform_circular_material_state(loss_template(), 0.0015, &state)
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap();
    let runtime = LinearMaterialStringRuntime::try_new(
        cx,
        specimen,
        Some(Pluck {
            station_frac: 0.4,
            height_m: 1e-3,
        }),
        incremental_ambient(),
        1.0,
        rate,
    )
    .unwrap();
    // Deliberately small synthetic cp=0.001 J/(kg K) makes the feedback
    // measurable in a short test. This is not a physical material dataset.
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.1,
                temperature_k: 400.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.2,
                temperature_k: 400.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    (runtime, card, curve)
}

#[test]
fn g1_current_heating_refines_to_independent_ohmic_temperature_solution() {
    with_string_cx(|cx, _| {
        let mut errors = Vec::new();
        for rate in [4000, 8000, 16000] {
            let (mechanical, card, curve) = heated_string_parts_with_conductor(
                cx,
                rate,
                false,
                InterpolationPolicy::LinearInside,
                true,
            );
            let mass = mechanical.specimen().mass_kg();
            let length = mechanical.specimen().string().length_m;
            let area = mechanical.specimen().area_m2();
            let quiet = LinearMaterialStringRuntime::try_new(
                cx,
                mechanical.specimen().clone(),
                None,
                incremental_ambient(),
                1.0,
                rate,
            )
            .unwrap();
            let mut runtime =
                ThermalMaterialStringRuntime::try_new(cx, quiet, card, curve, 0.0, 5.0).unwrap();
            let current = 0.05_f64;
            let r0 = 1e-6 * length / area;
            let steps = rate / 100;
            let mut joules = 0.0;
            for epoch in 0..u64::from(steps) {
                let frame = runtime
                    .step_with_current(cx, epoch, &[0.0], current)
                    .unwrap();
                let before = frame.conductor_before.resistance_ohm();
                close(
                    frame.coupled.external_heat_j,
                    current.powi(2) * before / f64::from(rate),
                );
                close(
                    frame.conductor_after.resistance_ohm(),
                    r0 * (1.0 + 0.2 * (frame.coupled.thermal.temperature_k() - 300.0)),
                );
                assert_eq!(frame.coupled.vibration.epoch, epoch + 1);
                assert_eq!(frame.coupled.vibration.dissipation.material_heat_j, 0.0);
                joules += frame.coupled.external_heat_j;
            }
            // cp=0.001 and rho_e(T)=rho_e0[1+0.2(T-300)] are synthetic.
            // Solve m cp dT/dt = I² R0[1+alpha(T-300)] independently.
            let exponent = current.powi(2) * r0 * 0.2 * 0.01 / (mass * 0.001);
            let exact = 300.0 + exponent.exp_m1() / 0.2;
            errors.push((runtime.thermal().temperature_k() - exact).abs());
            assert!((mass * runtime.thermal().specific_enthalpy_j_kg() - joules).abs() < 1e-12);
        }
        eprintln!("G1 current-heating temperature errors at4/8/16kHz: {errors:?}");
        assert!(errors[0] > 1e-6 && errors[0] < 0.01);
        assert!((0.45..0.55).contains(&(errors[1] / errors[0])));
        assert!((0.45..0.55).contains(&(errors[2] / errors[1])));
    });
}

#[test]
fn g3_current_sign_preserves_heating_while_resistance_damping_and_sound_respond() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) = heated_string_parts_with_conductor(
            cx,
            8000,
            false,
            InterpolationPolicy::LinearInside,
            true,
        );
        let mut positive =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 5.0).unwrap();
        let mut negative = positive.clone();
        let mut unpowered = positive.clone();
        let mut pressure_difference = 0.0;
        for epoch in 0..80 {
            let a = positive.step_with_current(cx, epoch, &[0.0], 0.05).unwrap();
            let b = negative
                .step_with_current(cx, epoch, &[0.0], -0.05)
                .unwrap();
            let control = unpowered.step_with_current(cx, epoch, &[0.0], 0.0).unwrap();
            assert_eq!(a.coupled, b.coupled);
            assert_eq!(a.conductor_after, b.conductor_after);
            assert!(a.conductor_after.resistance_ohm() > a.conductor_before.resistance_ohm());
            assert_eq!(
                a.conductor_after.material().card_identity(),
                positive.mechanical().specimen().material().card_identity()
            );
            assert_eq!(
                a.conductor_after.material().query_point(),
                positive.mechanical().specimen().material().query_point()
            );
            pressure_difference += (a.coupled.vibration.acoustic.observer_pressure_pa
                - control.coupled.vibration.acoustic.observer_pressure_pa)
                .abs();
        }
        assert!(positive.thermal().temperature_k() > unpowered.thermal().temperature_k() + 0.1);
        assert!(
            positive.mechanical().modes()[0].damping_ratio
                > unpowered.mechanical().modes()[0].damping_ratio
        );
        assert!(pressure_difference > 1e-9);
    });
}

#[test]
fn g4_current_refusal_preserves_electrothermal_state_and_retry() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mut missing =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 200.0).unwrap();
        assert!(matches!(
            missing.step_with_current(cx, 0, &[0.0], 0.05),
            Err(ThermalStringError::Admission(_))
        ));
        let (mechanical, card, curve) = heated_string_parts_with_conductor(
            cx,
            8000,
            false,
            InterpolationPolicy::LinearInside,
            true,
        );
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 200.0).unwrap();
        let mut retry = runtime.clone();
        assert!(matches!(
            runtime.step_with_current(cx, 0, &[0.0], f64::NAN),
            Err(ThermalStringError::Electrical(_))
        ));
        // Joule heating stays inside the thermal chart but exits the pinned
        // mechanical/electrical claim domain, after the vibration substep.
        assert!(matches!(
            runtime.step_with_current(cx, 0, &[0.0], 5.0),
            Err(ThermalStringError::Material(_))
        ));
        assert_eq!(runtime.thermal(), retry.thermal());
        assert_eq!(
            runtime.mechanical().specimen(),
            retry.mechanical().specimen()
        );
        assert_eq!(runtime.mechanical().states(), retry.mechanical().states());
        assert_eq!(runtime.mechanical().accepted_samples(), 0);
        for epoch in 0..16 {
            assert_eq!(
                runtime.step_with_current(cx, epoch, &[0.0], 0.05).unwrap(),
                retry.step_with_current(cx, epoch, &[0.0], 0.05).unwrap()
            );
            assert!(matches!(
                runtime.step_with_current(cx, epoch, &[0.0], 0.05),
                Err(ThermalStringError::Epoch { .. })
            ));
        }
    });
}

#[test]
fn g1_thermal_string_heating_updates_sourced_damping_and_pressure_at_one_epoch() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mut frozen = mechanical.clone();
        let mass = mechanical.specimen().mass_kg();
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 5.0).unwrap();
        let mut supplied = 0.0;
        let mut material_heat = 0.0;
        let mut pressure_difference = 0.0;
        for epoch in 0..160 {
            let external = if epoch % 2 == 0 { 1e-8 } else { -0.5e-8 };
            let frame = runtime.step(cx, epoch, &[0.0], external).unwrap();
            let reference = frozen.step(cx, &[0.0]).unwrap();
            assert_eq!(frame.vibration.epoch, epoch + 1);
            assert_eq!(runtime.mechanical().epoch(), epoch + 1);
            assert_eq!(runtime.mechanical().accepted_samples(), epoch + 1);
            assert_eq!(frame.thermal, runtime.thermal());
            let state = runtime.mechanical().specimen().material();
            assert_eq!(
                state
                    .query_point()
                    .iter()
                    .find(|(axis, _)| axis == "T")
                    .unwrap()
                    .1,
                frame.thermal.temperature_k()
            );
            close(
                state
                    .property(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
                    .unwrap()
                    .value_si(),
                1e7 * (1.0 + 0.2 * (frame.thermal.temperature_k() - 300.0)),
            );
            assert!(frame.energy_balance_residual_j.abs() <= frame.energy_roundoff_tolerance_j);
            assert_eq!(frame.vibration.dissipation.authored_loss_j, 0.0);
            supplied += external;
            material_heat += frame.vibration.dissipation.material_heat_j;
            pressure_difference += (frame.vibration.acoustic.observer_pressure_pa
                - reference.acoustic.observer_pressure_pa)
                .abs();
        }
        close(
            mass * runtime.thermal().specific_enthalpy_j_kg(),
            supplied + material_heat,
        );
        assert!(runtime.thermal().temperature_k() > 300.5);
        assert!(
            runtime.mechanical().modes()[0].damping_ratio > 1.05 * frozen.modes()[0].damping_ratio
        );
        assert!(pressure_difference > 1e-9);
    });
}

fn heated_mass_string_parts(
    cx: &fs_exec::Cx<'_>,
    mass_kg: f64,
) -> (
    LinearMaterialStringRuntime,
    MaterialCard,
    fs_material::phase::EquilibriumEnthalpyPhaseCurve,
) {
    let (base, card, curve) =
        heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
    let specimen = with_uniform_circular_material_and_constraints(
        base.specimen().string(),
        StringGeometryConstraint::FixedMass(mass_kg),
        base.specimen().material(),
        StringPrestress::FixedTension(20.0),
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    let mechanical = LinearMaterialStringRuntime::try_new(
        cx,
        specimen,
        Some(Pluck {
            station_frac: 0.4,
            height_m: 1e-3,
        }),
        incremental_ambient(),
        1.0,
        8000,
    )
    .unwrap();
    (mechanical, card, curve)
}

#[test]
fn g3_fixed_mass_thermal_string_matches_the_same_resolved_radius() {
    with_string_cx(|cx, _| {
        for mass in [0.002, 0.005] {
            let (mechanical, card, curve) = heated_mass_string_parts(cx, mass);
            let original = mechanical.specimen().clone();
            let radius_specimen = with_uniform_circular_material_state(
                original.string(),
                original.radius_m(),
                original.material(),
            )
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .unwrap();
            let radius_mechanical = LinearMaterialStringRuntime::try_new(
                cx,
                radius_specimen,
                Some(Pluck {
                    station_frac: 0.4,
                    height_m: 1e-3,
                }),
                incremental_ambient(),
                1.0,
                8000,
            )
            .unwrap();
            let initial_damping = mechanical.modes()[0].damping_ratio;
            let mut runtime = ThermalMaterialStringRuntime::try_new(
                cx,
                mechanical,
                card.clone(),
                curve.clone(),
                0.0,
                5.0,
            )
            .unwrap();
            let mut reference =
                ThermalMaterialStringRuntime::try_new(cx, radius_mechanical, card, curve, 0.0, 5.0)
                    .unwrap();
            let environment = string_environment(0.03, 0.01);
            let mut pressure_error = 0.0;
            let mut pressure_scale = 0.0;
            let mut heat = 0.0;
            for epoch in 0..80 {
                let (frame, control) = if epoch % 2 == 0 {
                    let actual = runtime
                        .step_with_thermal_transport(cx, epoch, &[0.0], |input| {
                            assert_eq!(input.mass_kg.to_bits(), mass.to_bits());
                            assert_eq!(
                                input.volume_m3,
                                original.area_m2() * original.string().length_m
                            );
                            let r = original.radius_m();
                            assert_eq!(
                                input.surface_area_m2,
                                2.0 * core::f64::consts::PI * r * (original.string().length_m + r)
                            );
                            environment.advance(cx, input)
                        })
                        .unwrap();
                    let expected = reference
                        .step_with_thermal_transport(cx, epoch, &[0.0], |input| {
                            environment.advance(cx, input)
                        })
                        .unwrap();
                    (actual.coupled, expected.coupled)
                } else {
                    (
                        runtime.step(cx, epoch, &[0.0], 1e-8).unwrap(),
                        reference.step(cx, epoch, &[0.0], 1e-8).unwrap(),
                    )
                };
                let current = runtime.mechanical().specimen();
                assert_eq!(
                    current.geometry_constraint(),
                    StringGeometryConstraint::FixedMass(mass)
                );
                assert_eq!(current.mass_kg().to_bits(), mass.to_bits());
                assert_eq!(
                    current.string().lin_density_kg_m.to_bits(),
                    (mass / original.string().length_m).to_bits()
                );
                assert_eq!(current.radius_m().to_bits(), original.radius_m().to_bits());
                assert_eq!(current.area_m2(), original.area_m2());
                assert_eq!(
                    current.string().bending_stiffness_n_m2,
                    original.string().bending_stiffness_n_m2
                );
                assert_eq!(frame.vibration.epoch, epoch + 1);
                // The radius description reconstructs mass with floating-point
                // roundoff. Compare physical trajectories, not receipt identities.
                assert!(
                    (frame.thermal.temperature_k() - control.thermal.temperature_k()).abs() < 2e-8
                );
                assert!(
                    frame.energy_balance_residual_j.abs()
                        <= frame.energy_roundoff_tolerance_j + frame.thermal_solve_tolerance_j
                );
                let p = frame.vibration.acoustic.observer_pressure_pa;
                let p_ref = control.vibration.acoustic.observer_pressure_pa;
                pressure_error += (p - p_ref).powi(2);
                pressure_scale += p_ref.powi(2);
                heat += frame.external_heat_j + frame.vibration.dissipation.material_heat_j;
            }
            assert!(pressure_scale > 0.0);
            assert!(pressure_error < 1e-14 * pressure_scale); // relative RMS < 1e-7
            assert!((mass * runtime.thermal().specific_enthalpy_j_kg() - heat).abs() < 1e-12);
            assert!(runtime.thermal().temperature_k() > 300.5);
            assert!(runtime.mechanical().modes()[0].damping_ratio > initial_damping);
        }
    });
}

#[test]
fn g4_fixed_mass_thermal_string_late_refusal_preserves_geometry_and_retry() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) = heated_mass_string_parts(cx, 0.003777001);
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 200.0).unwrap();
        let mut reference = runtime.clone();
        let mut hot = string_environment(100.0, 0.0);
        hot.temperature_k = 370.0;
        hot.radiation_temperature_k = 370.0;
        let mut thermal_solved = false;
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 0, &[0.0], |input| {
                let proposed = hot.advance(cx, input)?;
                thermal_solved = true;
                Ok::<_, fs_conduction::ConductionError>(proposed)
            }),
            Err(ThermalStringError::Material(_))
        ));
        assert!(thermal_solved, "refusal must follow the real thermal solve");
        assert_eq!(runtime.thermal(), reference.thermal());
        assert_eq!(
            runtime.mechanical().specimen(),
            reference.mechanical().specimen()
        );
        assert_eq!(
            runtime.mechanical().states(),
            reference.mechanical().states()
        );
        assert_eq!(runtime.mechanical().accepted_samples(), 0);
        let environment = string_environment(0.03, 0.01);
        for epoch in 0..16 {
            assert_eq!(
                runtime
                    .step_with_thermal_transport(cx, epoch, &[0.0], |input| environment
                        .advance(cx, input))
                    .unwrap(),
                reference
                    .step_with_thermal_transport(cx, epoch, &[0.0], |input| environment
                        .advance(cx, input))
                    .unwrap()
            );
        }
    });
}

fn string_environment(convection: f64, emissivity: f64) -> ThermalStringEnvironment {
    ThermalStringEnvironment {
        temperature_k: 320.0,
        radiation_temperature_k: 320.0,
        convection_w_per_m2_k: convection,
        transport: fs_conduction::lumped::LumpedThermalTransport::try_declared(100.0, emissivity)
            .unwrap(),
        enthalpy_tolerance_j_kg: 1e-13,
        maximum_thermal_residual_j: 1e-14,
    }
}

#[test]
fn g1_ambient_string_convection_matches_analytical_heating_cooling_and_refines() {
    with_string_cx(|cx, _| {
        // Independently integrated lumped Newton cooling: T = Ta + (T0-Ta) exp(-hAt/mc).
        // Two radii exercise the actual area/mass and V/A geometry binding.
        for radius in [0.001, 0.002] {
            for initial_temperature in [300.0, 330.0] {
                let mut errors = Vec::new();
                for rate in [4000, 8000, 16000] {
                    let (base, card, curve) =
                        heated_string_parts(cx, rate, false, InterpolationPolicy::LinearInside);
                    let specimen = with_uniform_circular_material_state(
                        base.specimen().string().clone(),
                        radius,
                        base.specimen().material(),
                    )
                    .unwrap()
                    .with_kelvin_voigt_bending_loss()
                    .unwrap();
                    let length = specimen.string().length_m;
                    let mass = specimen.mass_kg();
                    let mechanical = LinearMaterialStringRuntime::try_new(
                        cx,
                        specimen,
                        None,
                        incremental_ambient(),
                        1.0,
                        rate,
                    )
                    .unwrap();
                    let mut runtime = ThermalMaterialStringRuntime::try_new(
                        cx, mechanical, card, curve, 0.0, 40.0,
                    )
                    .unwrap();
                    runtime
                        .step(cx, 0, &[0.0], mass * 0.001 * (initial_temperature - 300.0))
                        .unwrap();
                    let environment = string_environment(0.03, 0.0);
                    let area = 2.0 * core::f64::consts::PI * radius * (length + radius);
                    let volume = core::f64::consts::PI * radius * radius * length;
                    let decay = environment.convection_w_per_m2_k * area / (mass * 0.001);
                    let duration = 0.02;
                    let mut external_heat = 0.0;
                    for sample in 0..rate / 50 {
                        let frame = runtime
                            .step_with_thermal_transport(
                                cx,
                                u64::from(sample) + 1,
                                &[0.0],
                                |input| environment.advance(cx, input),
                            )
                            .unwrap();
                        close(frame.transport.maximum_biot(), 0.03 * volume / area / 100.0);
                        let end = frame.transport.samples().last().unwrap();
                        assert_eq!(end.phase_state, frame.coupled.thermal);
                        assert_eq!(end.time_s, 1.0 / f64::from(rate));
                        assert_eq!(end.internal_power_w, 0.0);
                        assert_eq!(end.radiation_into_body_w, 0.0);
                        assert_eq!(frame.coupled.vibration.acoustic.observer_pressure_pa, 0.0);
                        assert_eq!(frame.coupled.vibration.epoch, u64::from(sample) + 2);
                        assert_eq!(
                            end.convection_into_body_w.is_sign_positive(),
                            initial_temperature < 320.0
                        );
                        external_heat += frame.coupled.external_heat_j;
                    }
                    let exact = 320.0 + (initial_temperature - 320.0) * (-decay * duration).exp();
                    errors.push((runtime.thermal().temperature_k() - exact).abs());
                    assert!(
                        (mass * 0.001 * (runtime.thermal().temperature_k() - initial_temperature)
                            - external_heat)
                            .abs()
                            < 1e-12
                    );
                }
                eprintln!(
                    "ambient convection radius={radius} initial_T={initial_temperature} errors={errors:?}"
                );
                assert!(errors[0] > 1e-4 && errors[0] < 0.1);
                assert!((0.45..0.55).contains(&(errors[1] / errors[0])));
                assert!((0.45..0.55).contains(&(errors[2] / errors[1])));
            }
        }
    });
}

#[test]
fn g1_ambient_string_fluxes_change_material_damping_and_pressure_without_double_heat() {
    with_string_cx(|cx, _| {
        for (convection, emissivity) in [(0.03, 0.0), (0.0, 0.01), (0.03, 0.01)] {
            let (mechanical, card, curve) =
                heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
            let mass = mechanical.specimen().mass_kg();
            let length = mechanical.specimen().string().length_m;
            let radius = 0.0015;
            let area = 2.0 * core::f64::consts::PI * radius * (length + radius);
            let mut runtime =
                ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 5.0)
                    .unwrap();
            let mut insulated = runtime.clone();
            let environment = string_environment(convection, emissivity);
            let mut total_heat = 0.0;
            let mut pressure_difference = 0.0;
            for epoch in 0..80 {
                let frame = runtime
                    .step_with_thermal_transport(cx, epoch, &[0.0], |input| {
                        environment.advance(cx, input)
                    })
                    .unwrap();
                let control = insulated.step(cx, epoch, &[0.0], 0.0).unwrap();
                let end = frame.transport.samples().last().unwrap();
                let temperature = end.phase_state.temperature_k();
                let convective_power = convection * area * (320.0 - temperature);
                let radiative_power = emissivity
                    * 5.670_374_419e-8
                    * area
                    * (320.0_f64.powi(4) - temperature.powi(4));
                assert!(
                    (end.convection_into_body_w - convective_power).abs()
                        <= 1e-12 * convective_power.abs().max(1e-20)
                );
                assert!(
                    (end.radiation_into_body_w - radiative_power).abs()
                        <= 1e-12 * radiative_power.abs().max(1e-20)
                );
                close(
                    frame.coupled.external_heat_j,
                    (convective_power + radiative_power) / 8000.0,
                );
                close(
                    end.internal_power_w / 8000.0,
                    frame.coupled.vibration.dissipation.material_heat_j,
                );
                assert!(
                    frame.coupled.energy_balance_residual_j.abs()
                        <= frame.coupled.energy_roundoff_tolerance_j
                            + frame.coupled.thermal_solve_tolerance_j
                );
                assert!(
                    (frame.coupled.energy_balance_residual_j + end.step_energy_residual_j).abs()
                        <= frame.coupled.energy_roundoff_tolerance_j
                );
                total_heat += frame.coupled.external_heat_j
                    + frame.coupled.vibration.dissipation.material_heat_j;
                pressure_difference += (frame.coupled.vibration.acoustic.observer_pressure_pa
                    - control.vibration.acoustic.observer_pressure_pa)
                    .abs();
            }
            assert!((mass * runtime.thermal().specific_enthalpy_j_kg() - total_heat).abs() < 1e-12);
            assert!(runtime.thermal().temperature_k() > insulated.thermal().temperature_k() + 1.0);
            assert!(
                runtime.mechanical().modes()[0].damping_ratio
                    > insulated.mechanical().modes()[0].damping_ratio
            );
            assert!(pressure_difference > 1e-9);
        }
    });
}

#[test]
fn g1_hot_enclosure_and_cool_fluid_drive_shared_string_temperature_and_sound() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mass = mechanical.specimen().mass_kg();
        let radius = mechanical.specimen().radius_m();
        let length = mechanical.specimen().string().length_m;
        let area = 2.0 * core::f64::consts::PI * radius * (length + radius);
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 5.0).unwrap();
        let mut control = runtime.clone();
        let mut environment = string_environment(0.03, 0.01);
        // Both runs start at the chart's 300 K lower boundary. Keep the fluid
        // there so the control does not demand out-of-domain cooling; the hot
        // enclosure raises the test string above the same fluid temperature.
        environment.temperature_k = 300.0;
        environment.radiation_temperature_k = 330.0;
        let mut control_environment = environment.clone();
        control_environment.radiation_temperature_k = 300.0;
        let mut total_heat = 0.0;
        let mut pressure_difference = 0.0;
        for epoch in 0..80 {
            let frame = runtime
                .step_with_thermal_transport(cx, epoch, &[0.0], |input| {
                    environment.advance(cx, input)
                })
                .unwrap();
            let comparison = control
                .step_with_thermal_transport(cx, epoch, &[0.0], |input| {
                    control_environment.advance(cx, input)
                })
                .unwrap();
            let end = frame.transport.samples().last().unwrap();
            let temperature = frame.coupled.thermal.temperature_k();
            close(
                end.convection_into_body_w,
                0.03 * area * (300.0 - temperature),
            );
            close(
                end.radiation_into_body_w,
                0.01 * 5.670_374_419e-8 * area * (330.0_f64.powi(4) - temperature.powi(4)),
            );
            assert!(end.convection_into_body_w < 0.0);
            assert!(end.radiation_into_body_w > 0.0);
            assert_eq!(frame.coupled.vibration.epoch, epoch + 1);
            assert!(
                frame.coupled.energy_balance_residual_j.abs()
                    <= frame.coupled.energy_roundoff_tolerance_j
                        + frame.coupled.thermal_solve_tolerance_j
            );
            total_heat +=
                frame.coupled.external_heat_j + frame.coupled.vibration.dissipation.material_heat_j;
            pressure_difference += (frame.coupled.vibration.acoustic.observer_pressure_pa
                - comparison.coupled.vibration.acoustic.observer_pressure_pa)
                .abs();
        }
        assert!((mass * runtime.thermal().specific_enthalpy_j_kg() - total_heat).abs() < 1e-12);
        assert!(runtime.thermal().temperature_k() > control.thermal().temperature_k() + 1.0);
        assert!(
            runtime.mechanical().modes()[0].damping_ratio
                > control.mechanical().modes()[0].damping_ratio
        );
        assert!(pressure_difference > 1e-9);

        let mut retry = runtime.clone();
        let mut invalid = environment.clone();
        invalid.radiation_temperature_k = f64::NAN;
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 80, &[0.0], |input| {
                invalid.advance(cx, input)
            }),
            Err(ThermalStringError::Transport(_))
        ));
        assert_eq!(runtime.thermal(), retry.thermal());
        assert_eq!(
            runtime.mechanical().specimen(),
            retry.mechanical().specimen()
        );
        assert_eq!(runtime.mechanical().states(), retry.mechanical().states());
        assert_eq!(
            runtime.mechanical().accepted_samples(),
            retry.mechanical().accepted_samples()
        );
        assert_eq!(
            runtime
                .step_with_thermal_transport(cx, 80, &[0.0], |input| environment.advance(cx, input))
                .unwrap(),
            retry
                .step_with_thermal_transport(cx, 80, &[0.0], |input| environment.advance(cx, input))
                .unwrap()
        );
    });
}

#[test]
fn g4_ambient_string_refusal_preserves_temperature_material_sound_and_retry() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 200.0).unwrap();
        let mut reference = runtime.clone();
        let environment = string_environment(0.03, 0.01);
        let mut bad = environment.clone();
        bad.transport =
            fs_conduction::lumped::LumpedThermalTransport::try_declared(1e-8, 0.01).unwrap();
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 0, &[0.0], |input| bad.advance(cx, input)),
            Err(ThermalStringError::Transport(_))
        ));
        bad = environment.clone();
        bad.enthalpy_tolerance_j_kg = 0.01;
        bad.maximum_thermal_residual_j = 0.0;
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 0, &[0.0], |input| bad.advance(cx, input)),
            Err(ThermalStringError::Transport(message)) if message.contains("ambient thermal solve exceeds its energy-residual budget")
        ));
        // Conduction succeeds inside its solid chart but T exits the selected
        // mechanical property's 340 K domain. Both owners must roll back.
        bad = string_environment(100.0, 0.0);
        bad.temperature_k = 370.0;
        bad.radiation_temperature_k = 370.0;
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 0, &[0.0], |input| bad.advance(cx, input)),
            Err(ThermalStringError::Material(_))
        ));
        // Thermal solver also supports latent heat, but this string has no
        // liquid mechanics/state transfer and must not publish a melted string.
        bad.temperature_k = 400.8;
        bad.radiation_temperature_k = 400.8;
        bad.convection_w_per_m2_k = 1000.0;
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 0, &[0.0], |input| bad.advance(cx, input)),
            Err(ThermalStringError::Admission(_))
        ));
        bad = environment.clone();
        bad.temperature_k = f64::NAN;
        assert!(
            runtime
                .step_with_thermal_transport(cx, 0, &[0.0], |input| bad.advance(cx, input))
                .is_err()
        );
        // Reject a valid transport result with a substituted phase chart.
        assert!(matches!(
            runtime.step_with_thermal_transport(cx, 0, &[0.0], |input| {
                let mut knots = input.curve.knots().to_vec();
                knots.last_mut().unwrap().specific_enthalpy_j_kg += 0.1;
                let foreign = fs_material::phase::EquilibriumEnthalpyPhaseCurve::try_new(
                    input.curve.material_card_identity(),
                    knots,
                )
                .unwrap();
                let mut proposed = environment.advance(cx, input)?;
                proposed.state = foreign
                    .state_at_specific_enthalpy(proposed.state.specific_enthalpy_j_kg())
                    .unwrap();
                Ok::<_, fs_conduction::ConductionError>(proposed)
            }),
            Err(ThermalStringError::Admission(_))
        ));
        with_string_cx(|cancelled, gate| {
            gate.request();
            assert!(matches!(
                runtime.step_with_thermal_transport(cancelled, 0, &[0.0], |input| environment
                    .advance(cancelled, input)),
                Err(ThermalStringError::Acoustic(
                    AcousticRealizeError::Cancelled
                ))
            ));
        });
        with_string_cx(|late_cx, gate| {
            let mut transport_solved = false;
            assert!(matches!(
                runtime.step_with_thermal_transport(late_cx, 0, &[0.0], |input| {
                    let proposed = environment.advance(late_cx, input)?;
                    transport_solved = true;
                    gate.request();
                    Ok::<_, fs_conduction::ConductionError>(proposed)
                }),
                Err(ThermalStringError::Acoustic(
                    AcousticRealizeError::Cancelled
                ))
            ));
            assert!(
                transport_solved,
                "late cancellation must reach the real thermal solve"
            );
        });
        assert_eq!(runtime.thermal(), reference.thermal());
        assert_eq!(
            runtime.mechanical().states(),
            reference.mechanical().states()
        );
        assert_eq!(
            runtime.mechanical().specimen(),
            reference.mechanical().specimen()
        );
        assert_eq!(runtime.mechanical().accepted_samples(), 0);
        for epoch in 0..16 {
            assert_eq!(
                runtime
                    .step_with_thermal_transport(cx, epoch, &[0.01], |input| environment
                        .advance(cx, input))
                    .unwrap(),
                reference
                    .step_with_thermal_transport(cx, epoch, &[0.01], |input| environment
                        .advance(cx, input))
                    .unwrap()
            );
            assert!(matches!(
                runtime.step_with_thermal_transport(cx, epoch, &[0.01], |input| environment
                    .advance(cx, input)),
                Err(ThermalStringError::Epoch { .. })
            ));
        }
        let mut resumed = runtime.clone();
        assert_eq!(
            runtime
                .step_with_thermal_transport(cx, 16, &[0.0], |input| environment.advance(cx, input))
                .unwrap(),
            resumed
                .step_with_thermal_transport(cx, 16, &[0.0], |input| environment.advance(cx, input))
                .unwrap()
        );
    });
}

#[test]
fn g4_thermal_string_refusals_and_duplicate_epoch_preserve_both_owners() {
    with_string_cx(|cx, _| {
        let (mechanical, card, curve) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mass = mechanical.specimen().mass_kg();
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 200.0).unwrap();
        let mut reference = runtime.clone();
        // Mechanics succeeds, then T=350 K lies outside the pinned viscosity
        // claim while still inside the solid thermal chart.
        assert!(matches!(
            runtime.step(cx, 0, &[0.0], mass * 0.05),
            Err(ThermalStringError::Material(_))
        ));
        assert!(matches!(
            runtime.step(cx, 0, &[0.0], mass * 0.15),
            Err(ThermalStringError::Admission(_))
        ));
        assert!(matches!(
            runtime.step(cx, 0, &[0.0], mass),
            Err(ThermalStringError::Phase(_))
        ));
        assert!(runtime.step(cx, 0, &[0.0], f64::NAN).is_err());
        with_string_cx(|cancelled, gate| {
            gate.request();
            assert!(matches!(
                runtime.step(cancelled, 0, &[0.0], 0.0),
                Err(ThermalStringError::Acoustic(
                    AcousticRealizeError::Cancelled
                ))
            ));
        });
        assert_eq!(runtime.thermal(), reference.thermal());
        assert_eq!(
            runtime.mechanical().states(),
            reference.mechanical().states()
        );
        assert_eq!(
            runtime.mechanical().specimen(),
            reference.mechanical().specimen()
        );
        for epoch in 0..32 {
            assert_eq!(
                runtime.step(cx, epoch, &[0.01], 1e-8).unwrap(),
                reference.step(cx, epoch, &[0.01], 1e-8).unwrap()
            );
            let accepted = runtime.thermal();
            assert!(matches!(
                runtime.step(cx, epoch, &[0.01], 1e-8),
                Err(ThermalStringError::Epoch { .. })
            ));
            assert_eq!(runtime.thermal(), accepted);
        }
        let mut resumed = runtime.clone();
        for epoch in 32..64 {
            assert_eq!(
                runtime.step(cx, epoch, &[0.0], 0.0).unwrap(),
                resumed.step(cx, epoch, &[0.0], 0.0).unwrap()
            );
        }
    });
}

#[test]
fn g0_thermal_string_admission_preserves_model_and_temperature_budgets() {
    with_string_cx(|cx, _| {
        for (elastic_curve, policy) in [
            (true, InterpolationPolicy::LinearInside),
            (false, InterpolationPolicy::TabulatedOnly),
        ] {
            let (mechanical, card, curve) = heated_string_parts(cx, 8000, elastic_curve, policy);
            assert!(matches!(
                ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 1.0),
                Err(ThermalStringError::Admission(_))
            ));
        }
        let (mechanical, card, curve) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mass = mechanical.specimen().mass_kg();
        let unbound = with_uniform_circular_material_state(
            loss_template(),
            0.0015,
            mechanical.specimen().material(),
        )
        .unwrap();
        let unbound = LinearMaterialStringRuntime::try_new(
            cx,
            unbound,
            None,
            incremental_ambient(),
            1.0,
            8000,
        )
        .unwrap();
        assert!(matches!(
            ThermalMaterialStringRuntime::try_new(
                cx,
                unbound,
                card.clone(),
                curve.clone(),
                0.0,
                1.0,
            ),
            Err(ThermalStringError::Admission(
                "initial string must already use its sourced Kelvin-Voigt bending law"
            ))
        ));
        for (h, limit) in [(0.001, 1.0), (0.0, 0.0), (0.0, f64::NAN)] {
            assert!(
                ThermalMaterialStringRuntime::try_new(
                    cx,
                    mechanical.clone(),
                    card.clone(),
                    curve.clone(),
                    h,
                    limit
                )
                .is_err()
            );
        }
        let mut runtime =
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 0.1).unwrap();
        let before = runtime.thermal();
        assert!(matches!(
            runtime.step(cx, 0, &[0.0], mass * 0.001),
            Err(ThermalStringError::Admission(_))
        ));
        assert_eq!(runtime.thermal(), before);
        assert_eq!(runtime.mechanical().accepted_samples(), 0);
    });
}

#[test]
fn g1_thermal_string_partition_converges_to_independent_coupled_ode() {
    with_string_cx(|cx, _| {
        let (mechanical, _, _) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let string = mechanical.specimen().string();
        let mass = mechanical.specimen().mass_kg();
        let k = core::f64::consts::PI / string.length_m;
        let moment = core::f64::consts::PI * 0.0015_f64.powi(4) / 4.0;
        let omega2 = (20.0 * k.powi(2) + 2e9 * moment * k.powi(4)) / string.lin_density_kg_m;
        let gas = fs_material::gas::GasState::try_new_moist_air(300.0, 101325.0, 0.0).unwrap();
        let air_c = fs_couple::air_path::oscillating_cylinder_air_resistance_per_length(
            0.0015,
            omega2.sqrt(),
            &gas,
        )
        .unwrap()
            / string.lin_density_kg_m;
        let initial = mechanical.states()[0];
        let rhs = |y: [f64; 3]| {
            let eta = 1e7 * (1.0 + 0.2 * y[2] / 0.001);
            let c = eta * moment * k.powi(4) / string.lin_density_kg_m;
            [
                y[1],
                -omega2 * y[0] - (c + air_c) * y[1],
                c * y[1] * y[1] / mass,
            ]
        };
        // Independent RK4 integrates q, v and h continuously. The production
        // algorithm instead uses exact frozen-temperature mechanical substeps.
        let reference = |steps: usize| {
            let dt = 0.02 / steps as f64;
            let mut y = [
                initial.displacement_m_sqrt_kg,
                initial.velocity_m_sqrt_kg_per_s,
                0.0,
            ];
            let shifted =
                |y: [f64; 3], k: [f64; 3], scale: f64| std::array::from_fn(|i| y[i] + scale * k[i]);
            for _ in 0..steps {
                let a = rhs(y);
                let b = rhs(shifted(y, a, dt / 2.0));
                let c = rhs(shifted(y, b, dt / 2.0));
                let d = rhs(shifted(y, c, dt));
                y = std::array::from_fn(|i| {
                    y[i] + dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i])
                });
            }
            y
        };
        let exact = reference(16384);
        for (a, b) in exact.iter().zip(reference(8192)) {
            assert!((a - b).abs() < 1e-11);
        }
        let mut errors = Vec::new();
        for rate in [4000, 8000, 16000] {
            let (mechanical, card, curve) =
                heated_string_parts(cx, rate, false, InterpolationPolicy::LinearInside);
            let mut runtime =
                ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 5.0)
                    .unwrap();
            for epoch in 0..u64::from(rate / 50) {
                runtime.step(cx, epoch, &[0.0], 0.0).unwrap();
            }
            let state = runtime.mechanical().states()[0];
            let error =
                ((state.displacement_m_sqrt_kg - exact[0]) / initial.displacement_m_sqrt_kg).abs()
                    + ((state.velocity_m_sqrt_kg_per_s - exact[1])
                        / (omega2.sqrt() * initial.displacement_m_sqrt_kg))
                        .abs()
                    + ((runtime.thermal().specific_enthalpy_j_kg() - exact[2]) / exact[2]).abs();
            errors.push(error);
        }
        eprintln!("G1 thermal string first-order errors: {errors:?}");
        assert!(errors[0] < 0.02 && errors[2] > 1e-8);
        for pair in errors.windows(2) {
            assert!(pair[1] < 0.6 * pair[0] && pair[1] > 0.4 * pair[0]);
        }
    });
}

#[test]
fn g0_thermal_string_refuses_a_negative_viscosity_between_positive_endpoints() {
    with_string_cx(|cx, _| {
        let (mechanical, card, _) =
            heated_string_parts(cx, 8000, false, InterpolationPolicy::LinearInside);
        let mut claims = ClaimSet::new();
        for (_, original) in card.claims().claims_ordered() {
            let mut claim = original.clone();
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                knots.insert(1, (320.0, -1e5));
            }
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap();
        let requirements: Vec<_> = mechanical
            .specimen()
            .material()
            .properties()
            .iter()
            .map(|property| property.requirement().clone())
            .collect();
        let point = QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("omega", 1.0)
            .unwrap();
        let state = resolve_material_state_point(
            &card,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let specimen = with_uniform_circular_material_state(loss_template(), 0.0015, &state)
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .unwrap();
        let mechanical = LinearMaterialStringRuntime::try_new(
            cx,
            specimen,
            None,
            incremental_ambient(),
            1.0,
            8000,
        )
        .unwrap();
        use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            card.content_hash(),
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 300.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 1.0,
                    temperature_k: 400.0,
                    liquid_mass_fraction: 1.0,
                    bulk_density_kg_m3: 1000.0,
                },
            ],
        )
        .unwrap();
        assert!(matches!(
            ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 100.0),
            Err(ThermalStringError::Admission(
                "thermal string viscosity must be nonnegative over its whole source curve"
            ))
        ));
    });
}

#[test]
fn g1_material_loss_drives_enthalpy_without_crediting_air_or_pressure_twice() {
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

    let eta = 1.0e5;
    let material = viscous_material(
        eta,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 1e6),
    );
    let specimen = with_uniform_circular_material_state(
        PrestressedString {
            n_modes: 2,
            ..loss_template()
        },
        0.0015,
        &material,
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    let string = specimen.string();
    let mass = specimen.mass_kg();
    let moment = specimen.second_moment_m4();
    // Synthetic constant-cp solid absorber using the existing enthalpy owner.
    // The final knot is required by that owner's full phase-curve contract;
    // this test never leaves the solid branch or exercises latent heat.
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        material.card_identity(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 293.15,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1000.0,
                temperature_k: 294.15,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 2000.0,
                temperature_k: 294.15,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    with_string_cx(|cx, _| {
        let mut thermal = curve.state_at_specific_enthalpy(0.0).unwrap();
        let mut runtime = LinearMaterialStringRuntime::try_new(
            cx,
            specimen,
            Some(Pluck {
                station_frac: 0.4,
                height_m: 1e-3,
            }),
            incremental_ambient(),
            1.0,
            48_000,
        )
        .unwrap();
        let gas = fs_material::gas::GasState::try_new_moist_air(300.0, 101_325.0, 0.0).unwrap();
        let energy = |runtime: &LinearMaterialStringRuntime| {
            runtime
                .modes()
                .iter()
                .zip(runtime.states())
                .map(|(m, s)| {
                    0.5 * (s.velocity_m_sqrt_kg_per_s.powi(2)
                        + (m.angular_frequency_rad_s * s.displacement_m_sqrt_kg).powi(2))
                })
                .sum::<f64>()
        };
        let initial = energy(&runtime);
        let mut heat = 0.0;
        let mut air = 0.0;
        let mut work = 0.0;
        let mut numerical = 0.0;
        let mut allowance = 0.0;
        for step in 0..64 {
            let forces = [0.1, -0.05];
            let before = runtime.states().to_vec();
            let modes = runtime.modes().to_vec();
            let frame = runtime.step(cx, &forces).unwrap();
            assert_eq!(frame.epoch, step + 1);
            assert_eq!(frame.epoch, runtime.epoch());
            let mut expected_heat = 0.0;
            let mut expected_air = 0.0;
            for i in 0..2 {
                let k = (i + 1) as f64 * core::f64::consts::PI / string.length_m;
                let omega = ((string.tension_n * k * k
                    + string.bending_stiffness_n_m2 * k.powi(4))
                    / string.lin_density_kg_m)
                    .sqrt();
                let solid_c = eta * moment * k.powi(4) / string.lin_density_kg_m;
                let air_c = fs_couple::air_path::oscillating_cylinder_air_resistance_per_length(
                    string.width_m / 2.0,
                    omega,
                    &gas,
                )
                .unwrap()
                    / string.lin_density_kg_m;
                close(modes[i].angular_frequency_rad_s, omega);
                close(2.0 * modes[i].damping_ratio * omega, solid_c + air_c);
                // Independent analytic velocity and Simpson integration of v^2.
                // Coefficients are constant throughout this one ZOH interval.
                let a = (solid_c + air_c) / 2.0;
                let b = (omega * omega - a * a).sqrt();
                let q = before[i].displacement_m_sqrt_kg;
                let v = before[i].velocity_m_sqrt_kg_per_s;
                let dt = 1.0 / 48_000.0;
                let mut integral = 0.0;
                for j in 0..=128 {
                    let t = dt * j as f64 / 128.0;
                    let velocity = (-a * t).exp()
                        * (v * (b * t).cos()
                            - (a * v + omega * omega * q - forces[i]) / b * (b * t).sin());
                    let weight = if j == 0 || j == 128 {
                        1.0
                    } else if j % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    integral += weight * velocity * velocity;
                }
                integral *= dt / (3.0 * 128.0);
                expected_heat += solid_c * integral;
                expected_air += air_c * integral;
            }
            let loss = frame.dissipation;
            assert_eq!(loss.authored_loss_j, 0.0);
            assert!(loss.material_heat_j > 0.0 && loss.air_loss_j > 0.0);
            assert!(
                (loss.material_heat_j - expected_heat).abs()
                    <= loss.roundoff_tolerance_j + 1e-10 * expected_heat
            );
            assert!(
                (loss.air_loss_j - expected_air).abs()
                    <= loss.roundoff_tolerance_j + 1e-10 * expected_air
            );
            thermal = curve
                .advance_specific_energy(thermal, loss.material_heat_j / mass)
                .unwrap();
            heat += loss.material_heat_j;
            air += loss.air_loss_j;
            work += frame.acoustic.input_work_j;
            numerical += loss.roundoff_residual_j;
            allowance += loss.roundoff_tolerance_j;
        }
        let delta_mechanical = energy(&runtime) - initial;
        let absorbed = mass * thermal.specific_enthalpy_j_kg();
        close(absorbed, heat);
        assert!(thermal.temperature_k() > 293.15 && thermal.solid_mass_fraction() == 1.0);
        let residual = delta_mechanical + absorbed + air + numerical - work;
        let tolerance = allowance + 256.0 * f64::EPSILON * (initial + work.abs());
        assert!(
            residual.abs() <= tolerance,
            "balance {residual} vs {tolerance}"
        );
        // A second credit (including crediting all damping to the solid), or
        // omitting the mechanical side, cannot hide inside numerical allowance.
        assert!((residual + absorbed).abs() > tolerance);
        assert!((residual + air).abs() > tolerance);
        assert!((residual - delta_mechanical).abs() > tolerance);
    });
}

#[test]
fn g3_rebinding_loss_law_replaces_destinations_and_keeps_roundoff_out_of_heat() {
    let material = viscous_material(
        1e5,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 1e6),
    );
    let bind =
        |template| with_uniform_circular_material_state(template, 0.0015, &material).unwrap();
    with_string_cx(|cx, _| {
        let mut runtime = LinearMaterialStringRuntime::try_new(
            cx,
            bind(loss_template())
                .with_kelvin_voigt_bending_loss()
                .unwrap(),
            Some(Pluck {
                station_frac: 0.4,
                height_m: 1e-3,
            }),
            incremental_ambient(),
            1.0,
            48_000,
        )
        .unwrap();
        let material_frame = runtime.step(cx, &[0.0]).unwrap();
        assert!(material_frame.dissipation.material_heat_j > 0.0);
        assert!(material_frame.dissipation.air_loss_j > 0.0);
        assert_eq!(material_frame.dissipation.authored_loss_j, 0.0);
        for alpha_per_s in [2.0, 0.0] {
            let next = bind(PrestressedString {
                rayleigh: Some(RayleighParams {
                    alpha_per_s,
                    beta_s: 0.0,
                }),
                ..template()
            });
            let before = runtime.states().to_vec();
            assert_eq!(runtime.rebind(cx, next, 0.0).unwrap().parameter_work_j, 0.0);
            assert_eq!(runtime.states(), before);
            let frame = runtime.step(cx, &[0.0]).unwrap();
            assert_eq!(frame.dissipation.material_heat_j, 0.0);
            // Existing Rayleigh semantics replace the entire damping model.
            assert_eq!(frame.dissipation.air_loss_j, 0.0);
            if alpha_per_s > 0.0 {
                assert!(frame.dissipation.authored_loss_j > 0.0);
            } else {
                assert_eq!(frame.dissipation.authored_loss_j, 0.0);
                assert_eq!(
                    frame.dissipation.roundoff_residual_j,
                    frame.acoustic.viscous_dissipation_j
                );
            }
            assert!(
                frame.dissipation.roundoff_residual_j.abs()
                    <= frame.dissipation.roundoff_tolerance_j
            );
        }
    });
}

fn relaxing_material(
    delta: f64,
    tau: f64,
    band: (f64, f64),
    curve: Option<&str>,
) -> ResolvedMaterialStatePoint {
    state_with_domain(
        "synthetic standard-linear-solid",
        &[
            ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
            (
                "young_modulus",
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0e9,
            ),
            (
                EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0e9,
            ),
            (
                RELAXING_BENDING_MODULUS_PROPERTY,
                QuantitySpec::dimensional(Pressure::DIMS),
                delta,
            ),
            (
                BENDING_RELAXATION_TIME_PROPERTY,
                QuantitySpec::dimensional(Time::DIMS),
                tau,
            ),
        ],
        ValidityDomain::unconstrained().with("omega", band.0, band.1),
        QueryPoint::new().with("omega", band.0).unwrap(),
        curve,
    )
}

fn prony_material(
    terms: &[(f64, f64)],
    curve: Option<&str>,
) -> (ResolvedMaterialStatePoint, Vec<BendingRelaxationProperties>) {
    let pairs: Vec<_> = (0..terms.len())
        .map(|j| BendingRelaxationProperties {
            modulus: format!("relaxing_bending_modulus_{j}"),
            relaxation_time: format!("bending_relaxation_time_{j}"),
        })
        .collect();
    let mut properties = vec![
        ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
        (
            "young_modulus",
            QuantitySpec::dimensional(Pressure::DIMS),
            2.0e9,
        ),
        (
            EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
            QuantitySpec::dimensional(Pressure::DIMS),
            2.0e9,
        ),
    ];
    for (pair, &(delta, tau)) in pairs.iter().zip(terms) {
        properties.push((
            pair.modulus.as_str(),
            QuantitySpec::dimensional(Pressure::DIMS),
            delta,
        ));
        properties.push((
            pair.relaxation_time.as_str(),
            QuantitySpec::dimensional(Time::DIMS),
            tau,
        ));
    }
    let material = state_with_domain(
        "synthetic Prony spectrum",
        &properties,
        ValidityDomain::unconstrained().with("omega", 1.0, 20_000.0),
        QueryPoint::new().with("omega", 1.0).unwrap(),
        curve,
    );
    (material, pairs)
}

#[test]
fn g3_frequency_units_preserve_string_loss_bands_and_pressure() {
    use fs_qty::semantic::FrequencyConvention::{Angular, Cyclic};
    let frequency = |kind| {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(kind),
            ValueForm::Static,
        ))
    };
    let (_, pairs) = prony_material(&[(1e10, 0.005), (5e9, 0.02)], None);
    let modulus = QuantitySpec::dimensional(Pressure::DIMS);
    let time = QuantitySpec::dimensional(Time::DIMS);
    let rows = [
        (
            "density",
            QuantitySpec::dimensional(Density::DIMS),
            1000.0,
            (1.0, 140.0),
        ),
        ("young_modulus", modulus, 2e9, (2.0, 130.0)),
        (
            EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
            modulus,
            2e9,
            (3.0, 120.0),
        ),
        (
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(DynViscosity::DIMS),
            3e7,
            (4.0, 90.0),
        ),
        (pairs[0].modulus.as_str(), modulus, 1e10, (5.0, 110.0)),
        (pairs[0].relaxation_time.as_str(), time, 0.005, (6.0, 100.0)),
        (pairs[1].modulus.as_str(), modulus, 5e9, (7.0, 80.0)),
        (pairs[1].relaxation_time.as_str(), time, 0.02, (8.0, 60.0)),
    ];
    let properties: Vec<_> = rows
        .iter()
        .map(|&(key, quantity, value, _)| (key, quantity, value))
        .collect();
    let mut reference: [Option<Vec<f64>>; 2] = [None, None];
    for (axis, quantity, scale) in [
        ("omega", None, core::f64::consts::TAU),
        ("response", Some(frequency(Angular)), core::f64::consts::TAU),
        ("frequency", Some(frequency(Cyclic)), 1.0),
        ("omega", Some(frequency(Cyclic)), 1.0),
    ] {
        let point = match quantity {
            Some(q) => QueryPoint::new()
                .with_quantity(axis, q, 50.0 * scale)
                .unwrap(),
            None => QueryPoint::new().with(axis, 50.0 * scale).unwrap(),
        };
        let state = state_with_property_domains(
            "frequency-unit string",
            &properties,
            |key| {
                let (_, _, _, (lo, hi)) = rows.iter().find(|row| row.0 == key).unwrap();
                match quantity {
                    Some(q) => ValidityDomain::unconstrained().with_quantity(
                        axis,
                        q,
                        lo * scale,
                        hi * scale,
                    ),
                    None => ValidityDomain::unconstrained().with(axis, lo * scale, hi * scale),
                }
            },
            point,
            None,
        );
        for (law, reference) in reference.iter_mut().enumerate() {
            let specimen =
                with_uniform_circular_material_state(loss_template(), 0.002, &state).unwrap();
            let (string, expected_hz) = if law == 0 {
                (
                    specimen.with_kelvin_voigt_bending_loss().unwrap().string(),
                    (4.0, 90.0),
                )
            } else {
                (
                    specimen.with_prony_bending_loss(&pairs).unwrap().string(),
                    (8.0, 60.0),
                )
            };
            let band = if law == 0 {
                string
                    .kelvin_voigt_bending
                    .as_ref()
                    .unwrap()
                    .omega_band_rad_s
            } else {
                string.relaxation_bending.as_ref().unwrap().omega_band_rad_s
            };
            assert_eq!(
                band,
                (
                    expected_hz.0 * core::f64::consts::TAU,
                    expected_hz.1 * core::f64::consts::TAU
                )
            );
            let actual = pressure(string);
            assert!(actual.iter().any(|p| p.abs() > 0.0));
            match reference {
                Some(expected) => assert_eq!(
                    actual.as_slice(),
                    expected.as_slice(),
                    "{axis}/{quantity:?}, law={law}"
                ),
                None => *reference = Some(actual),
            }
        }
    }
    // Narrow only the density source: it must constrain both loss laws, not
    // just eta or the relaxation-time claims. The query remains in every band.
    let state = state_with_property_domains(
        "restricted density",
        &properties,
        |key| {
            let band = if key == "density" {
                (20.0, 30.0)
            } else {
                (1.0, 140.0)
            };
            ValidityDomain::unconstrained().with_quantity(
                "frequency",
                frequency(Cyclic),
                band.0,
                band.1,
            )
        },
        QueryPoint::new()
            .with_quantity("frequency", frequency(Cyclic), 25.0)
            .unwrap(),
        None,
    );
    for prony in [false, true] {
        let specimen =
            with_uniform_circular_material_state(loss_template(), 0.002, &state).unwrap();
        let specimen = if prony {
            specimen.with_prony_bending_loss(&pairs)
        } else {
            specimen.with_kelvin_voigt_bending_loss()
        }
        .unwrap();
        let error = realize_assembly(&assembly(specimen.string())).unwrap_err();
        assert!(error.to_string().contains("frequenc"), "{error}");
    }
}

#[test]
fn g1_temperature_curves_drive_every_prony_branch_and_emitted_pressure() {
    let (_, pairs) = prony_material(&[(1e10, 0.004), (5e9, 0.02)], None);
    let pressure_quantity = QuantitySpec::dimensional(Pressure::DIMS);
    let time_quantity = QuantitySpec::dimensional(Time::DIMS);
    let properties = [
        (
            "density",
            QuantitySpec::dimensional(Density::DIMS),
            1000.0,
            1200.0,
        ),
        ("young_modulus", pressure_quantity, 2e9, 1e9),
        (
            EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
            pressure_quantity,
            2e9,
            1e9,
        ),
        (pairs[0].modulus.as_str(), pressure_quantity, 1e10, 5e9),
        (
            pairs[0].relaxation_time.as_str(),
            time_quantity,
            0.004,
            0.008,
        ),
        (pairs[1].modulus.as_str(), pressure_quantity, 5e9, 7e9),
        (pairs[1].relaxation_time.as_str(), time_quantity, 0.02, 0.01),
    ];
    let domain = ValidityDomain::unconstrained()
        .with("T", 300.0, 400.0)
        .with("omega", 1.0, 20_000.0);
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(key, quantity, cold, hot) in &properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(key, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(300.0, cold), (400.0, hot)],
                    dims: quantity.dims(),
                },
                validity: domain.clone(),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic thermal Prony coefficients".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key,
                quantity,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal Prony solid".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut responses = Vec::new();
    for temperature in [300.0, 325.0, 400.0] {
        let point = QueryPoint::new()
            .with("T", temperature)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap();
        let material = resolve_material_state_point(
            &card,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let fraction = (temperature - 300.0) / 100.0;
        let expected: Vec<_> = properties
            .iter()
            .map(|&(key, quantity, cold, hot)| (key, quantity, cold + fraction * (hot - cold)))
            .collect();
        // Independent scalar source exercises the same solver with explicitly
        // computed coefficients; the curve receipts must remain interpolated.
        let reference = state_with_domain(
            "thermal Prony scalar reference",
            &expected,
            domain.clone(),
            point,
            None,
        );
        let bound = with_uniform_circular_material_state(loss_template(), 0.002, &material)
            .unwrap()
            .with_prony_bending_loss(&pairs)
            .unwrap();
        let law = bound.string().relaxation_bending.unwrap();
        assert_eq!(law.material_state_identity, Some(material.identity()));
        let second_moment = core::f64::consts::PI * 0.004_f64.powi(4) / 64.0;
        for (j, branch) in law.branches.iter().enumerate() {
            close(
                branch.relaxing_stiffness_n_m2,
                expected[3 + 2 * j].2 * second_moment,
            );
            close(branch.relaxation_time_s, expected[4 + 2 * j].2);
        }
        for &(key, _, value) in &expected {
            let property = material.property(key).unwrap();
            close(property.value_si(), value);
            if temperature == 325.0 {
                assert!(matches!(
                    property.answer().receipt.decision,
                    fs_matdb::EvaluationDecision::LinearInside { .. }
                ));
            }
        }
        let reference = with_uniform_circular_material_state(loss_template(), 0.002, &reference)
            .unwrap()
            .with_prony_bending_loss(&pairs)
            .unwrap();
        let run = |string| {
            let mut input = assembly(string);
            input.sample_rate_hz = 48_000;
            input.duration_s = 0.06;
            realize_assembly(&input).unwrap().pressure_pa
        };
        let actual = run(bound.string());
        let expected = run(reference.string());
        let peak = expected.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        assert!(peak > 0.0, "the pressure path must be live");
        assert_eq!(actual.len(), expected.len());
        assert!(
            actual
                .iter()
                .zip(&expected)
                .all(|(a, e)| (a - e).abs() < 1e-10 * peak)
        );
        responses.push(actual);
    }
    assert_ne!(responses[0], responses[1]);
    assert_ne!(responses[1], responses[2]);
}

#[test]
fn g1_prony_material_binding_rebinds_every_selected_branch() {
    let terms = [(1.0e10, 0.005), (0.0, 1.0e-9), (5.0e9, 0.02)];
    let (material, pairs) = prony_material(&terms, None);
    let specimen = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_prony_bending_loss(&pairs)
        .unwrap();
    let string = specimen.string();
    let law = string.relaxation_bending.as_ref().unwrap();
    assert_eq!(law.source_properties.as_ref(), Some(&pairs));
    assert_eq!(law.material_state_identity, Some(material.identity()));
    assert_eq!(specimen.material(), &material);
    assert_eq!(law.branches.len(), 3);
    for (branch, &(delta, tau)) in law.branches.iter().zip(&terms) {
        if delta == 0.0 {
            assert_eq!(branch.relaxing_stiffness_n_m2.to_bits(), 0.0_f64.to_bits());
        } else {
            close(
                branch.relaxing_stiffness_n_m2,
                delta * core::f64::consts::PI * 0.004_f64.powi(4) / 64.0,
            );
        }
        close(branch.relaxation_time_s, tau);
    }
    let replacement_terms: Vec<_> = terms.iter().map(|&(e, t)| (2.0 * e, 3.0 * t)).collect();
    let (replacement, _) = prony_material(&replacement_terms, None);
    let rebound =
        with_uniform_circular_material_state(string.clone(), 0.004, &replacement).unwrap();
    let rebound_string = rebound.string();
    let next = rebound_string.relaxation_bending.as_ref().unwrap();
    for (before, after) in law.branches.iter().zip(&next.branches) {
        if before.relaxing_stiffness_n_m2 == 0.0 {
            assert_eq!(after.relaxing_stiffness_n_m2.to_bits(), 0.0_f64.to_bits());
        } else {
            close(
                after.relaxing_stiffness_n_m2,
                32.0 * before.relaxing_stiffness_n_m2,
            );
        }
        close(after.relaxation_time_s, 3.0 * before.relaxation_time_s);
    }
    assert_eq!(next.material_state_identity, Some(replacement.identity()));
    assert_eq!(next.source_properties.as_ref(), Some(&pairs));
    for selectors in [None, Some(vec![])] {
        let mut authored = string.clone();
        authored
            .relaxation_bending
            .as_mut()
            .unwrap()
            .source_properties = selectors;
        assert!(with_uniform_circular_material_state(authored, 0.004, &replacement).is_err());
    }
    let (missing_later_branch, _) = prony_material(&terms[..2], None);
    assert!(with_uniform_circular_material_state(string, 0.002, &missing_later_branch).is_err());
}

#[test]
fn g0_prony_material_refuses_bad_later_sources_and_duplicate_pairs() {
    let terms = [(1.0e10, 0.005), (5.0e9, 0.02)];
    let (material, pairs) = prony_material(&terms, None);
    let base = with_uniform_circular_material_state(loss_template(), 0.002, &material).unwrap();
    assert!(
        base.clone()
            .with_prony_bending_loss(&[pairs[0].clone(), pairs[0].clone()])
            .is_err()
    );
    let mut missing = pairs.clone();
    missing[1].relaxation_time = "absent_time".into();
    assert!(base.clone().with_prony_bending_loss(&missing).is_err());
    missing[1] = pairs[1].clone();
    missing[1].modulus = "absent_modulus".into();
    assert!(base.with_prony_bending_loss(&missing).is_err());
    for key in [&pairs[1].modulus, &pairs[1].relaxation_time] {
        let (curve, _) = prony_material(&terms, Some(key));
        assert!(
            with_uniform_circular_material_state(loss_template(), 0.002, &curve)
                .unwrap()
                .with_prony_bending_loss(&pairs)
                .is_err()
        );
    }
    let (empty, _) = prony_material(&[], None);
    let elastic = with_uniform_circular_material_state(loss_template(), 0.002, &empty)
        .unwrap()
        .with_prony_bending_loss(&[])
        .unwrap();
    let (zero, pair) = prony_material(&[(0.0, 1.0e-9)], None);
    let zero = with_uniform_circular_material_state(loss_template(), 0.002, &zero)
        .unwrap()
        .with_prony_bending_loss(&pair)
        .unwrap();
    assert_eq!(
        pressure(elastic.string()),
        pressure(zero.string()),
        "zero branches add no memory state or loss"
    );
}

#[test]
fn g0_prony_binding_intersects_all_sources_and_requires_declared_time_bands() {
    let properties = [
        ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
        (
            "young_modulus",
            QuantitySpec::dimensional(Pressure::DIMS),
            2.0e9,
        ),
        (
            EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
            QuantitySpec::dimensional(Pressure::DIMS),
            2.0e9,
        ),
        ("delta_a", QuantitySpec::dimensional(Pressure::DIMS), 1.0e10),
        ("tau_a", QuantitySpec::dimensional(Time::DIMS), 0.005),
        ("delta_b", QuantitySpec::dimensional(Pressure::DIMS), 5.0e9),
        ("tau_b", QuantitySpec::dimensional(Time::DIMS), 0.02),
    ];
    let pairs = [
        BendingRelaxationProperties {
            modulus: "delta_a".into(),
            relaxation_time: "tau_a".into(),
        },
        BendingRelaxationProperties {
            modulus: "delta_b".into(),
            relaxation_time: "tau_b".into(),
        },
    ];
    for omitted in [
        None,
        Some("tau_b"),
        Some(EQUILIBRIUM_YOUNG_MODULUS_PROPERTY),
    ] {
        let material = state_with_property_domains(
            "individual Prony domains",
            &properties,
            |key| {
                if Some(key) == omitted {
                    return ValidityDomain::unconstrained();
                }
                let (lo, hi) = match key {
                    "density" => (200.0, 800.0),
                    "delta_b" => (100.0, 900.0),
                    "tau_b" => (50.0, 1000.0),
                    _ => (1.0, 1200.0),
                };
                ValidityDomain::unconstrained().with("omega", lo, hi)
            },
            QueryPoint::new().with("omega", 250.0).unwrap(),
            None,
        );
        let base = with_uniform_circular_material_state(loss_template(), 0.002, &material).unwrap();
        let bound = base.clone().with_prony_bending_loss(&pairs);
        if omitted == Some("tau_b") {
            assert!(bound.is_err(), "the second arm needs its own band");
        } else {
            assert_eq!(
                bound
                    .unwrap()
                    .string()
                    .relaxation_bending
                    .unwrap()
                    .omega_band_rad_s,
                (200.0, 800.0)
            );
        }
        if omitted == Some(EQUILIBRIUM_YOUNG_MODULUS_PROPERTY) {
            assert!(
                base.with_prony_bending_loss(&[]).is_err(),
                "empty spectrum cannot borrow a density-only band"
            );
        }
    }
}

#[test]
fn g0_prony_admission_checks_total_stiffness_and_each_active_rate() {
    let (material, pairs) = prony_material(&[(1.0e10, 0.005); 3], None);
    let mut string = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_prony_bending_loss(&pairs)
        .unwrap()
        .string();
    let w = string_mode_omega(&string, 1);
    string.relaxation_bending.as_mut().unwrap().omega_band_rad_s = (0.9 * w, 1.2 * w);
    assert!(
        realize_assembly(&assembly(string.clone())).is_err(),
        "sum exceeds band"
    );
    let mut single = string.clone();
    single
        .relaxation_bending
        .as_mut()
        .unwrap()
        .branches
        .truncate(1);
    assert!(
        realize_assembly(&assembly(single)).is_ok(),
        "each branch alone fits"
    );
    let law = string.relaxation_bending.as_mut().unwrap();
    law.omega_band_rad_s = (1.0, 20_000.0);
    law.branches[2].relaxation_time_s = 1.0e-9;
    assert!(
        realize_assembly(&assembly(string.clone()))
            .unwrap_err()
            .to_string()
            .contains("dt/tau")
    );
    string.relaxation_bending.as_mut().unwrap().branches[2].relaxing_stiffness_n_m2 = 0.0;
    assert!(
        realize_assembly(&assembly(string)).is_ok(),
        "zero stiffness has no active pole"
    );
}

#[test]
fn g1_prony_pressure_converges_to_independent_hereditary_dynamics() {
    let terms = [(1.0e10, 0.005), (5.0e9, 0.02)];
    let (material, pairs) = prony_material(&terms, None);
    let mut string = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_prony_bending_loss(&pairs)
        .unwrap()
        .string();
    string.axial_stiffness_n = 0.0;
    let mut errors = Vec::new();
    for rate in [8_000, 16_000] {
        let mut scene = assembly(string.clone());
        scene.duration_s = 0.12;
        scene.sample_rate_hz = rate;
        let output = realize_assembly(&scene).unwrap();
        // Independent physical displacement + hereditary strain coordinates,
        // integrated by RK4; no pHS structure, memory scaling or stepper reused.
        let pi = core::f64::consts::PI;
        let k = pi / string.length_m;
        let bend_factor = string.width_m.powi(2) * k.powi(4) / (16.0 * 1000.0);
        let w2 = string.tension_n * k * k / string.lin_density_kg_m + 2.0e9 * bend_factor;
        let a = [terms[0].0 * bend_factor, terms[1].0 * bend_factor];
        let gas = output.gas;
        let drag = (2.0 * pi * gas.dynamic_viscosity
            + 2.0
                * pi
                * string.width_m
                * (gas.dynamic_viscosity * gas.density * w2.sqrt() / 2.0).sqrt())
            / string.lin_density_kg_m;
        let rhs = |x: [f64; 4]| {
            [
                x[1],
                -w2 * x[0] - drag * x[1] - a[0] * (x[0] - x[2]) - a[1] * (x[0] - x[3]),
                (x[0] - x[2]) / terms[0].1,
                (x[0] - x[3]) / terms[1].1,
            ]
        };
        let q0 = 8.0 * scene.pluck.unwrap().height_m / pi.powi(2);
        let mut x = [q0, 0.0, q0, q0]; // held until fully relaxed
        let h = 1.0 / (f64::from(rate) * 16.0);
        let weight = gas.density * string.width_m * string.length_m
            / (2.0 * pi.powi(2) * scene.listener.distance_m);
        let mut reference = Vec::with_capacity(output.pressure_pa.len());
        for _ in &output.pressure_pa {
            for _ in 0..16 {
                let k1 = rhs(x);
                let k2 = rhs(core::array::from_fn(|i| x[i] + 0.5 * h * k1[i]));
                let k3 = rhs(core::array::from_fn(|i| x[i] + 0.5 * h * k2[i]));
                let k4 = rhs(core::array::from_fn(|i| x[i] + h * k3[i]));
                x = core::array::from_fn(|i| {
                    x[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
                });
            }
            reference.push(weight * rhs(x)[1]);
        }
        // The shared propagation operator is outside this constitutive check;
        // apply it to both histories, after independent mechanical integration.
        fs_couple::air_path::absorb_pressure_history(
            &mut reference,
            1.0 / f64::from(rate),
            scene.listener.distance_m,
            &gas,
            scene.ambient.relative_humidity,
        );
        let (mut error, mut norm) = (0.0, 0.0);
        for (pressure, expected) in output.pressure_pa.iter().zip(reference) {
            error += (pressure - expected).powi(2);
            norm += expected.powi(2);
        }
        errors.push((error / norm).sqrt());
    }
    assert!(errors[1] < 5.0e-4, "pressure error {errors:?}");
    assert!(
        errors[1] < 0.3 * errors[0],
        "second-order refinement {errors:?}"
    );
    eprintln!("G1 Prony relative RMS pressure error at 8/16 kHz: {errors:?}");
}

#[test]
fn g3_prony_branch_permutation_and_equivalent_splitting_preserve_pressure() {
    let run = |terms: &[(f64, f64)]| {
        let (material, pairs) = prony_material(terms, None);
        let mut string = with_uniform_circular_material_state(loss_template(), 0.002, &material)
            .unwrap()
            .with_prony_bending_loss(&pairs)
            .unwrap()
            .string();
        string.polarization_detune = 0.01;
        pressure(string)
    };
    let base = run(&[(1.0e10, 0.005), (5.0e9, 0.02)]);
    for equivalent in [
        run(&[(5.0e9, 0.02), (1.0e10, 0.005)]),
        run(&[(4.0e9, 0.005), (0.0, 1.0e-9), (5.0e9, 0.02), (6.0e9, 0.005)]),
    ] {
        let error: f64 = base
            .iter()
            .zip(&equivalent)
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        let norm: f64 = base.iter().map(|a| a * a).sum();
        assert!(
            (error / norm).sqrt() < 1.0e-7,
            "equivalent spectra changed pressure"
        );
    }
}

#[test]
fn g1_sls_material_binding_recomputes_geometry_and_preserves_sources() {
    let material = relaxing_material(1.0e10, 0.005, (1.0, 20_000.0), None);
    let specimen = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_standard_linear_solid_bending_loss()
        .unwrap();
    let law = specimen.string().relaxation_bending.unwrap();
    close(
        law.branches[0].relaxing_stiffness_n_m2,
        1.0e10 * core::f64::consts::PI * 0.004_f64.powi(4) / 64.0,
    );
    close(law.branches[0].relaxation_time_s, 0.005);
    assert_eq!(law.material_state_identity, Some(material.identity()));
    assert_eq!(specimen.material(), &material);
    assert!(
        !material
            .property(RELAXING_BENDING_MODULUS_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .observation_backed
    );
    let replacement = relaxing_material(2.0e10, 0.01, (1.0, 30_000.0), None);
    let rebound =
        with_uniform_circular_material_state(specimen.string(), 0.004, &replacement).unwrap();
    let new_law = rebound.string().relaxation_bending.unwrap();
    close(
        new_law.branches[0].relaxing_stiffness_n_m2 / law.branches[0].relaxing_stiffness_n_m2,
        32.0,
    );
    close(new_law.branches[0].relaxation_time_s, 0.01);
    assert_eq!(
        new_law.material_state_identity,
        Some(replacement.identity())
    );
    assert_eq!(new_law.omega_band_rad_s, (1.0, 30_000.0));
    assert!(
        with_uniform_circular_material_state(
            specimen.string(),
            0.004,
            &elastic("missing relaxation", 1000.0, 2.0e9)
        )
        .is_err()
    );
}

#[test]
fn g0_sls_material_refuses_missing_equilibrium_wrong_time_curves_and_duplicate_loss() {
    let bind = |material: &ResolvedMaterialStatePoint| {
        with_uniform_circular_material_state(loss_template(), 0.002, material)
            .unwrap()
            .with_standard_linear_solid_bending_loss()
    };
    assert!(bind(&elastic("missing", 1000.0, 2.0e9)).is_err());
    for (equilibrium, tau_dims, band) in [
        (3.0e9, Time::DIMS, true),
        (2.0e9, Pressure::DIMS, true),
        (2.0e9, Time::DIMS, false),
    ] {
        let material = state_with_domain(
            "invalid SLS",
            &[
                ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
                (
                    "young_modulus",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    2.0e9,
                ),
                (
                    EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
                    QuantitySpec::dimensional(Pressure::DIMS),
                    equilibrium,
                ),
                (
                    RELAXING_BENDING_MODULUS_PROPERTY,
                    QuantitySpec::dimensional(Pressure::DIMS),
                    1.0e10,
                ),
                (
                    BENDING_RELAXATION_TIME_PROPERTY,
                    QuantitySpec::dimensional(tau_dims),
                    0.005,
                ),
            ],
            if band {
                ValidityDomain::unconstrained().with("omega", 1.0, 20_000.0)
            } else {
                ValidityDomain::unconstrained()
            },
            QueryPoint::new().with("omega", 1.0).unwrap(),
            None,
        );
        assert!(bind(&material).is_err());
    }
    for key in [
        "density",
        "young_modulus",
        EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
        RELAXING_BENDING_MODULUS_PROPERTY,
        BENDING_RELAXATION_TIME_PROPERTY,
    ] {
        assert!(
            bind(&relaxing_material(
                1.0e10,
                0.005,
                (1.0, 20_000.0),
                Some(key)
            ))
            .is_err()
        );
    }
    for (delta, tau) in [(-1.0, 0.005), (1.0e10, 0.0)] {
        assert!(bind(&relaxing_material(delta, tau, (1.0, 20_000.0), None)).is_err());
    }
    let material = relaxing_material(1.0e10, 0.005, (1.0, 20_000.0), None);
    assert!(
        with_uniform_circular_material_state(template(), 0.002, &material)
            .unwrap()
            .with_standard_linear_solid_bending_loss()
            .is_err()
    );
    let specimen = bind(&material).unwrap();
    assert!(specimen.clone().with_kelvin_voigt_bending_loss().is_err());
    let mut duplicate = specimen.string();
    duplicate.damping_ratio = 0.01;
    assert!(realize_assembly(&assembly(duplicate)).is_err());
}

#[test]
fn g0_sls_realization_checks_instantaneous_band_nyquist_and_polarizations() {
    let material = relaxing_material(1.0e10, 0.005, (1.0, 20_000.0), None);
    let string = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_standard_linear_solid_bending_loss()
        .unwrap()
        .string();
    let w = string_mode_omega(&string, 1);
    let mut narrow = string.clone();
    narrow
        .relaxation_bending
        .as_mut()
        .unwrap()
        .omega_band_rad_s
        .1 = 1.01 * w;
    assert!(
        realize_assembly(&assembly(narrow)).is_err(),
        "instantaneous stiffness exceeds admitted band"
    );
    let mut valid = string.clone();
    valid.relaxation_bending.as_mut().unwrap().omega_band_rad_s = (0.9 * w, 1.2 * w);
    assert!(realize_assembly(&assembly(valid.clone())).is_ok());
    for modified in [
        PrestressedString {
            n_modes: 2,
            ..valid.clone()
        },
        PrestressedString {
            polarization_detune: 0.5,
            ..valid.clone()
        },
        PrestressedString {
            moving_end: true,
            ..valid
        },
        PrestressedString {
            moving_end: true,
            polarization_detune: 0.01,
            ..string.clone()
        },
    ] {
        assert!(realize_assembly(&assembly(modified)).is_err());
    }
    let mut too_stiff = string.clone();
    too_stiff.relaxation_bending.as_mut().unwrap().branches[0].relaxing_stiffness_n_m2 *= 1.0e6;
    too_stiff
        .relaxation_bending
        .as_mut()
        .unwrap()
        .omega_band_rad_s
        .1 = 1.0e9;
    assert!(
        realize_assembly(&assembly(too_stiff)).is_err(),
        "material band cannot license aliased modes"
    );
    let mut fast = string.clone();
    fast.relaxation_bending.as_mut().unwrap().branches[0].relaxation_time_s = 1.0e-6;
    assert!(
        realize_assembly(&assembly(fast))
            .unwrap_err()
            .to_string()
            .contains("dt/tau")
    );
    let mut coarse = assembly(string);
    coarse.sample_rate_hz = 1000; // below Nyquist, but not an accurate midpoint phase
    coarse
        .string
        .as_mut()
        .unwrap()
        .relaxation_bending
        .as_mut()
        .unwrap()
        .branches[0]
        .relaxation_time_s = 0.05;
    assert!(
        realize_assembly(&coarse)
            .unwrap_err()
            .to_string()
            .contains("reference phase")
    );
}

#[test]
fn g1_sls_pressure_pitch_and_decay_match_independent_characteristic_roots() {
    for (nonlinear, moving_end) in [(false, false), (true, false), (true, true)] {
        let mut pitches = Vec::new();
        for delta in [0.0, 1.0e10] {
            let material = relaxing_material(delta, 0.005, (1.0, 20_000.0), None);
            let mut string =
                with_uniform_circular_material_state(loss_template(), 0.002, &material)
                    .unwrap()
                    .with_standard_linear_solid_bending_loss()
                    .unwrap()
                    .string();
            string.moving_end = moving_end;
            if !nonlinear {
                string.axial_stiffness_n = 0.0;
            }
            let mut scene = assembly(string.clone());
            scene.sample_rate_hz = 16_000;
            scene.duration_s = 0.4;
            let output = realize_assembly(&scene).unwrap();
            let pi = core::f64::consts::PI;
            let k = if moving_end { 0.5 } else { 1.0 } * pi / string.length_m;
            // Independent diameter-form reduction, not the runtime branch builder.
            let w2 = string.tension_n * k * k / string.lin_density_kg_m
                + 2.0e9 * string.width_m.powi(2) * k.powi(4) / (16.0 * 1000.0);
            let a = delta * string.width_m.powi(2) * k.powi(4) / (16.0 * 1000.0);
            let gas = output.gas;
            let c = (2.0 * pi * gas.dynamic_viscosity
                + 2.0
                    * pi
                    * string.width_m
                    * (gas.dynamic_viscosity * gas.density * w2.sqrt() / 2.0).sqrt())
                / string.lin_density_kg_m;
            // q'' + c q' + w0² q + a(q-v)=0, v'=(q-v)/tau.
            // The real root lies in [-1/tau,0]; divide the cubic to get the pair.
            let rate = 200.0;
            let (mut lo, mut hi) = (-rate, 0.0);
            for _ in 0..80 {
                let s = 0.5 * (lo + hi);
                if (s * s + c * s + w2) * (s + rate) + a * s > 0.0 {
                    hi = s;
                } else {
                    lo = s;
                }
            }
            let real = 0.5 * (lo + hi);
            let decay = 0.5 * (c + rate + real);
            let hz = (-w2 * rate / real - decay * decay).sqrt() / (2.0 * pi);
            let peaks: Vec<_> = output
                .pressure_pa
                .windows(3)
                .enumerate()
                .filter(|(i, p)| {
                    (960..5600).contains(i) && p[1] > 0.0 && p[1] > p[0] && p[1] >= p[2]
                })
                .map(|(i, p)| ((i + 1) as f64 / 16_000.0, p[1]))
                .collect();
            assert!(peaks.len() >= 4, "live oscillation required");
            let (t0, p0) = peaks[0];
            let (t1, p1) = *peaks.last().unwrap();
            let measured_hz = (peaks.len() - 1) as f64 / (t1 - t0);
            let measured_decay = (p0 / p1).ln() / (t1 - t0);
            assert!(
                (measured_hz / hz - 1.0).abs() < 0.005,
                "{measured_hz} vs {hz} Hz"
            );
            if delta > 0.0 {
                assert!(
                    (measured_decay / decay - 1.0).abs() < 0.03,
                    "nonlinear={nonlinear}, moving={moving_end}: {measured_decay} vs {decay} /s"
                );
            }
            pitches.push(measured_hz);
            eprintln!(
                "G1 SLS deltaE={delta}, nonlinear={nonlinear}, moving={moving_end}: {measured_hz} Hz / {measured_decay} per s; oracle {hz} / {decay}"
            );
        }
        assert!(
            pitches[1] > 1.005 * pitches[0],
            "relaxation must change storage, not just attenuate PCM"
        );
    }
}

#[test]
fn g1_material_bending_viscosity_follows_geometry_and_retains_sources() {
    let material = viscous_material(
        3.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 20_000.0),
    );
    let specimen = with_uniform_circular_material_state(loss_template(), 0.0005, &material)
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap();
    let loss = specimen.string().kelvin_voigt_bending.unwrap();
    close(
        loss.viscous_stiffness_n_m2_s,
        3.0e7 * core::f64::consts::PI * 0.001_f64.powi(4) / 64.0,
    );
    assert_eq!(loss.material_state_identity, Some(material.identity()));
    assert_eq!(specimen.material(), &material);
    assert!(
        !specimen
            .material()
            .property(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .observation_backed
    );
    let replacement = viscous_material(
        6.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 30_000.0),
    );
    let rebound =
        with_uniform_circular_material_state(specimen.string(), 0.001, &replacement).unwrap();
    let new_loss = rebound.string().kelvin_voigt_bending.unwrap();
    close(
        new_loss.viscous_stiffness_n_m2_s / loss.viscous_stiffness_n_m2_s,
        32.0,
    );
    assert_eq!(
        new_loss.material_state_identity,
        Some(replacement.identity())
    );
    assert_eq!(new_loss.omega_band_rad_s.0.to_bits(), 1.0_f64.to_bits());
    assert_eq!(
        new_loss.omega_band_rad_s.1.to_bits(),
        30_000.0_f64.to_bits()
    );
    assert_ne!(rebound.specimen_identity(), specimen.specimen_identity());
    assert!(
        with_uniform_circular_material_state(
            specimen.string(),
            0.001,
            &elastic("missing viscosity", 1000.0, 2.0e9)
        )
        .is_err()
    );
}

#[test]
fn g3_material_viscosity_changes_pressure_decay_without_changing_elastic_modes() {
    let pi = core::f64::consts::PI;
    // The same material loss enters linear, nonlinear KC, and moving-end modes.
    for (nonlinear, moving_end) in [(false, false), (true, false), (true, true)] {
        let mut decays = Vec::new();
        let mut reference_frequency = None;
        for eta in [0.0, 3.0e7, 1.2e8] {
            let material = viscous_material(
                eta,
                QuantitySpec::dimensional(DynViscosity::DIMS),
                (1.0, 20_000.0),
            );
            let specimen = with_uniform_circular_material_state(loss_template(), 0.0005, &material)
                .unwrap()
                .with_kelvin_voigt_bending_loss()
                .unwrap();
            let mut string = specimen.string();
            string.moving_end = moving_end;
            if !nonlinear {
                string.axial_stiffness_n = 0.0;
            }
            let wave_number = if moving_end { 0.5 } else { 1.0 } * pi / string.length_m;
            let omega = ((string.tension_n * wave_number.powi(2)
                + string.bending_stiffness_n_m2 * wave_number.powi(4))
                / string.lin_density_kg_m)
                .sqrt();
            if let Some(first) = reference_frequency {
                assert_eq!(omega.to_bits(), first);
            }
            reference_frequency = Some(omega.to_bits());
            let mut scene = assembly(string.clone());
            scene.sample_rate_hz = 16_000;
            scene.duration_s = 0.4;
            let output = realize_assembly(&scene).unwrap();
            assert!(output.pressure_pa.iter().all(|p| p.is_finite()));
            // Independent diameter-form projection: I/A = d²/16. The static
            // tensile energy must not be counted as dissipative bending energy.
            let bending_decay =
                eta * string.width_m.powi(2) * wave_number.powi(4) / (32.0 * 1000.0);
            let gas = output.gas;
            let resistance = 2.0 * pi * gas.dynamic_viscosity
                + 2.0
                    * pi
                    * string.width_m
                    * (gas.dynamic_viscosity * gas.density * omega / 2.0).sqrt();
            let expected = bending_decay + resistance / (2.0 * string.lin_density_kg_m);
            let peaks: Vec<_> = output
                .pressure_pa
                .windows(3)
                .enumerate()
                .filter(|(i, p)| {
                    (800..5600).contains(i) && p[1] > 0.0 && p[1] > p[0] && p[1] >= p[2]
                })
                .map(|(i, p)| ((i + 1) as f64 / 16_000.0, p[1]))
                .collect();
            assert!(peaks.len() > 20, "live pressure oscillation required");
            let (t0, p0) = peaks[0];
            let (t1, p1) = *peaks.last().unwrap();
            let measured = (p0 / p1).ln() / (t1 - t0);
            assert!(
                (measured / expected - 1.0).abs() < 0.02,
                "nonlinear={nonlinear}, moving={moving_end}, eta={eta}: {measured} vs {expected} /s"
            );
            decays.push(measured);
        }
        let ratio = (decays[2] - decays[0]) / (decays[1] - decays[0]);
        assert!(
            (ratio - 4.0).abs() < 0.12,
            "material-only decay ratio {ratio}"
        );
        eprintln!(
            "G3 Kelvin-Voigt viscosity 0, 3e7, 1.2e8 Pa s; nonlinear={nonlinear}, moving={moving_end}: {decays:?} /s"
        );
    }
}

#[test]
fn g0_material_bending_loss_refuses_missing_data_aliases_and_duplicate_losses() {
    let bind = |material: &ResolvedMaterialStatePoint| {
        with_uniform_circular_material_state(loss_template(), 0.0005, material)
            .unwrap()
            .with_kelvin_voigt_bending_loss()
    };
    assert!(bind(&elastic("missing", 1000.0, 2.0e9)).is_err());
    assert!(
        bind(&state(
            "no band",
            &[
                ("density", Density::DIMS, 1000.0),
                ("young_modulus", Pressure::DIMS, 2.0e9),
                (
                    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                    DynViscosity::DIMS,
                    3.0e7
                )
            ]
        ))
        .is_err()
    );
    for (eta, quantity) in [
        (-1.0, QuantitySpec::dimensional(DynViscosity::DIMS)),
        (3.0e7, QuantitySpec::dimensional(Pressure::DIMS)),
    ] {
        assert!(bind(&viscous_material(eta, quantity, (1.0, 20_000.0))).is_err());
    }
    let material = viscous_material(
        3.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 20_000.0),
    );
    assert!(
        with_uniform_circular_material_state(template(), 0.0005, &material)
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .is_err()
    );
    let specimen = bind(&material).unwrap();
    let mut duplicate = specimen.string();
    duplicate.damping_ratio = 0.01;
    assert!(realize_assembly(&assembly(duplicate)).is_err());
    duplicate = specimen.string();
    duplicate.rayleigh = template().rayleigh;
    assert!(realize_assembly(&assembly(duplicate)).is_err());
}

#[test]
fn g0_material_bending_loss_does_not_freeze_a_sampled_frequency_curve() {
    for key in [
        "density",
        "young_modulus",
        KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
    ] {
        let material = state_with_domain(
            "frequency-dependent synthetic solid",
            &[
                ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
                (
                    "young_modulus",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    2.0e9,
                ),
                (
                    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                    QuantitySpec::dimensional(DynViscosity::DIMS),
                    3.0e7,
                ),
            ],
            ValidityDomain::unconstrained().with("omega", 1.0, 20_000.0),
            QueryPoint::new().with("omega", 1000.0).unwrap(),
            Some(key),
        );
        assert!(matches!(
            material.property(key).unwrap().answer().receipt.decision,
            fs_matdb::EvaluationDecision::LinearInside { .. }
        ));
        let sampled =
            with_uniform_circular_material_state(loss_template(), 0.0005, &material).unwrap();
        let error = sampled.with_kelvin_voigt_bending_loss().unwrap_err();
        assert!(
            error.to_string().contains("validity-wide scalar constants"),
            "{error}"
        );
    }
}

#[test]
fn g0_material_bending_loss_checks_every_retained_frequency_and_polarization() {
    let omega = string_mode_omega(
        &with_uniform_circular_material_state(
            loss_template(),
            0.0005,
            &elastic("base", 1000.0, 2.0e9),
        )
        .unwrap()
        .string(),
        1,
    );
    let material = viscous_material(
        3.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (omega * 0.9, omega * 1.1),
    );
    let string = with_uniform_circular_material_state(loss_template(), 0.0005, &material)
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap()
        .string();
    assert!(realize_assembly(&assembly(string.clone())).is_ok());
    for modified in [
        PrestressedString {
            n_modes: 2,
            ..string.clone()
        },
        PrestressedString {
            polarization_detune: 0.5,
            ..string.clone()
        },
        PrestressedString {
            moving_end: true,
            ..string
        },
    ] {
        assert!(realize_assembly(&assembly(modified)).is_err());
    }
}

#[test]
fn g3_material_density_changes_pitch_in_the_existing_pressure_simulation() {
    let light =
        with_uniform_circular_material_state(template(), 0.0005, &elastic("light", 1000.0, 2.0e9))
            .unwrap();
    let heavy =
        with_uniform_circular_material_state(template(), 0.0005, &elastic("heavy", 4000.0, 2.0e9))
            .unwrap();
    let a = pressure(light.string());
    let b = pressure(heavy.string());
    assert!(a.iter().chain(&b).all(|p| p.is_finite()));
    assert_eq!(a, pressure(light.string()));
    let f_a = measured_hz(&a);
    let f_b = measured_hz(&b);
    eprintln!("G3 material density 1000 -> 4000 kg/m3: pressure pitch {f_a:.6} -> {f_b:.6} Hz");
    assert!(
        (f_a / f_b - 2.0).abs() < 0.015,
        "density ×4 must halve pitch: {f_a}, {f_b}"
    );
    let diameter: f64 = 0.001;
    let area = core::f64::consts::PI * diameter.powi(2) / 4.0;
    let moment = core::f64::consts::PI * diameter.powi(4) / 64.0;
    let k = core::f64::consts::PI / 0.5;
    let expected = ((20.0 * k.powi(2) + 2.0e9 * moment * k.powi(4)) / (1000.0 * area)).sqrt()
        / (2.0 * core::f64::consts::PI);
    assert!(
        (f_a / expected - 1.0).abs() < 0.01,
        "{f_a} vs beam reference {expected}"
    );
}

#[test]
fn g0_binding_refuses_missing_wrong_dimension_and_unrepresentable_inputs() {
    let material = elastic("valid", 1000.0, 2.0e9);
    for radius in [
        0.0,
        -0.001,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(with_uniform_circular_material_state(template(), radius, &material).is_err());
    }
    for bad in [
        state("missing", &[("density", Density::DIMS, 1000.0)]),
        state(
            "wrong-dimension",
            &[
                ("density", Dims::NONE, 1000.0),
                ("young_modulus", Pressure::DIMS, 2.0e9),
            ],
        ),
        elastic("negative", -1000.0, 2.0e9),
    ] {
        assert!(with_uniform_circular_material_state(template(), 0.0005, &bad).is_err());
    }
    let original = template();
    assert!(
        with_uniform_circular_material_state(
            PrestressedString {
                tension_n: f64::INFINITY,
                ..original.clone()
            },
            0.0005,
            &material,
        )
        .is_err()
    );
    assert_eq!(original, template());
}

#[test]
fn g0_binding_preserves_quantity_kinds_and_value_forms_at_the_solver_boundary() {
    let density = QuantitySpec::dimensional(Density::DIMS);
    let modulus = QuantitySpec::dimensional(Pressure::DIMS);
    for (property, quantity) in [
        (
            "density",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::MassConcentration,
                ValueForm::Static,
            )),
        ),
        (
            "young_modulus",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AcousticPressure,
                ValueForm::Rms,
            )),
        ),
    ] {
        let mut properties = [
            ("density", density, 1000.0),
            ("young_modulus", modulus, 2.0e9),
        ];
        let entry = properties.iter_mut().find(|p| p.0 == property).unwrap();
        assert_eq!(entry.1.dims(), quantity.dims());
        entry.1 = quantity;
        let material = state_with_quantities("explicit-semantic-alias", &properties);
        assert_eq!(
            material
                .property(property)
                .unwrap()
                .requirement()
                .quantity(),
            quantity
        );
        assert!(with_uniform_circular_material_state(template(), 0.0005, &material).is_err());
        assert!(
            with_uniform_circular_material_and_constraints(
                template(),
                StringGeometryConstraint::FixedMass(0.001),
                &material,
                StringPrestress::FixedTension(20.0),
            )
            .is_err()
        );
    }
}

#[test]
fn g0_fixed_mass_refuses_invalid_and_unrepresentable_geometry() {
    let material = elastic("mass-admission", 1000.0, 2.0e9);
    for mass in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(
            with_uniform_circular_material_and_constraints(
                template(),
                StringGeometryConstraint::FixedMass(mass),
                &material,
                StringPrestress::FixedTension(20.0),
            )
            .is_err()
        );
    }
    for length in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(
            with_uniform_circular_material_and_constraints(
                PrestressedString {
                    length_m: length,
                    ..template()
                },
                StringGeometryConstraint::FixedMass(0.001),
                &material,
                StringPrestress::FixedTension(20.0),
            )
            .is_err()
        );
    }
}

#[test]
fn g1_prestress_prescriptions_preserve_material_authority_and_solve_beam_tension() {
    let material = elastic("same-specimen", 1000.0, 2.0e9);
    let fixed = with_uniform_circular_material_state(template(), 0.0005, &material).unwrap();
    let extension = StringPrestress::FixedExtension {
        stress_free_length_m: 0.498,
        linear_strain_limit: 0.005,
    };
    let extended = with_uniform_circular_material_and_prestress(
        PrestressedString {
            tension_n: f64::NAN,
            ..template()
        },
        0.0005,
        &material,
        extension,
    )
    .unwrap();
    let area = core::f64::consts::PI * 0.001_f64.powi(2) / 4.0;
    close(extended.string().tension_n, 2.0e9 * area * 0.002 / 0.498);
    assert_eq!(extended.prestress(), extension);
    assert_eq!(extended.material(), &material);
    assert_eq!(extended.specimen_identity(), fixed.specimen_identity());
    assert_eq!(extended.mass_kg().to_bits(), fixed.mass_kg().to_bits());
    assert_eq!(extended.string().rayleigh, fixed.string().rayleigh);
    assert_ne!(
        extended.string().tension_n.to_bits(),
        fixed.string().tension_n.to_bits()
    );

    let target = StringPrestress::TargetFundamentalHz(160.0);
    let tuned = with_uniform_circular_material_and_prestress(template(), 0.0005, &material, target)
        .unwrap();
    let moment = core::f64::consts::PI * 0.001_f64.powi(4) / 64.0;
    let expected_tension = 4.0 * 1000.0 * area * 0.5_f64.powi(2) * 160.0_f64.powi(2)
        - core::f64::consts::PI.powi(2) * 2.0e9 * moment / 0.5_f64.powi(2);
    close(tuned.string().tension_n, expected_tension);
    close(
        string_mode_omega(&tuned.string(), 1) / core::f64::consts::TAU,
        160.0,
    );
    assert_eq!(tuned.prestress(), target);
    assert_eq!(tuned.specimen_identity(), fixed.specimen_identity());
    assert_eq!(fixed.prestress(), StringPrestress::FixedTension(20.0));
}

#[test]
fn g3_constraint_choice_changes_material_swap_pressure_response() {
    let extension = StringPrestress::FixedExtension {
        stress_free_length_m: 0.498,
        linear_strain_limit: 0.005,
    };
    let compliant = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("compliant", 1000.0, 2.0e9),
        extension,
    )
    .unwrap();
    let stiff = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("stiff", 1000.0, 8.0e9),
        extension,
    )
    .unwrap();
    close(stiff.string().tension_n / compliant.string().tension_n, 4.0);
    let compliant_hz = measured_hz(&pressure(compliant.string()));
    let stiff_hz = measured_hz(&pressure(stiff.string()));
    assert!(
        (stiff_hz / compliant_hz - 2.0).abs() < 0.015,
        "fixed extension and E x4 must double small-amplitude pitch: {compliant_hz}, {stiff_hz}"
    );

    let target = StringPrestress::TargetFundamentalHz(160.0);
    let light = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("light-tuned", 1000.0, 2.0e9),
        target,
    )
    .unwrap();
    let heavy = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("heavy-tuned", 4000.0, 2.0e9),
        target,
    )
    .unwrap();
    assert!(heavy.string().tension_n > 4.0 * light.string().tension_n);
    let light_hz = measured_hz(&pressure(light.string()));
    let heavy_hz = measured_hz(&pressure(heavy.string()));
    for actual in [light_hz, heavy_hz] {
        assert!(
            (actual / 160.0 - 1.0).abs() < 0.01,
            "target pressure pitch: {actual}"
        );
    }
    eprintln!(
        "G3 fixed extension E x4: {compliant_hz:.6} -> {stiff_hz:.6} Hz; target pitch density x4: {light_hz:.6} -> {heavy_hz:.6} Hz"
    );
}

#[test]
fn g0_prestress_refuses_slack_overstrain_invalid_targets_and_moving_end_mismatch() {
    let material = elastic("valid", 1000.0, 2.0e9);
    let original = template();
    let bind = |prestress| {
        with_uniform_circular_material_and_prestress(original.clone(), 0.0005, &material, prestress)
    };
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(bind(StringPrestress::FixedTension(invalid)).is_err());
        assert!(bind(StringPrestress::TargetFundamentalHz(invalid)).is_err());
        assert!(
            bind(StringPrestress::FixedExtension {
                stress_free_length_m: invalid,
                linear_strain_limit: 0.005,
            })
            .is_err()
        );
        assert!(
            bind(StringPrestress::FixedExtension {
                stress_free_length_m: 0.498,
                linear_strain_limit: invalid,
            })
            .is_err()
        );
    }
    for stress_free_length_m in [0.5, 0.6, f64::MIN_POSITIVE] {
        assert!(
            bind(StringPrestress::FixedExtension {
                stress_free_length_m,
                linear_strain_limit: 0.005,
            })
            .is_err()
        );
    }
    assert!(
        bind(StringPrestress::FixedExtension {
            stress_free_length_m: 0.498,
            linear_strain_limit: 0.001,
        })
        .is_err()
    );
    for hz in [1.0e-10, f64::MAX] {
        assert!(bind(StringPrestress::TargetFundamentalHz(hz)).is_err());
    }
    let moving = PrestressedString {
        moving_end: true,
        ..original.clone()
    };
    assert!(
        with_uniform_circular_material_and_prestress(
            moving.clone(),
            0.0005,
            &material,
            StringPrestress::FixedTension(20.0),
        )
        .is_ok()
    );
    for prestress in [
        StringPrestress::TargetFundamentalHz(160.0),
        StringPrestress::FixedExtension {
            stress_free_length_m: 0.498,
            linear_strain_limit: 0.005,
        },
    ] {
        assert!(
            with_uniform_circular_material_and_prestress(
                moving.clone(),
                0.0005,
                &material,
                prestress,
            )
            .is_err()
        );
    }
    assert_eq!(original, template());
}
