//! G0/G1/G3/G4: sourced hereditary vibration and physical heat, not measured fidelity.
use core::convert::Infallible;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_couple::thermovisco_dynamics::{Controls, Error, Geometry, Oscillator};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermovisco::source::{ResolvedThermoMaxwell, ThermoMaxwellSourceQuery};
use fs_material::thermovisco::{
    Limits,
    thermal::{FeedbackConfig, FeedbackError},
};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

struct Sources {
    mechanical: ResolvedThermoMaxwell,
    caloric: ResolvedHeatCapacityCurve,
    transport: LumpedThermalTransport,
}
fn temperature(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn sources(cp: f64, density: f64, branch: bool, thermal_shift: bool) -> Sources {
    sources_at_pressure(cp, density, branch, thermal_shift, None)
}
fn sources_at_pressure(
    cp: f64,
    density: f64,
    branch: bool,
    thermal_shift: bool,
    pressure_pa: Option<f64>,
) -> Sources {
    let mut claims = ClaimSet::new();
    let quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let pressure_quantity =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    for (name, dims, values) in [
        ("E_inf", Pressure::DIMS, [1e6; 4]),
        ("E0", Pressure::DIMS, [3e6; 4]),
        ("tau0", Time::DIMS, [0.005; 4]),
        ("scale", Dims::NONE, [1.0; 4]),
        (
            "shift",
            Dims::NONE,
            if thermal_shift {
                [1.0, 1.0, 0.05, 0.05]
            } else {
                [1.0; 4]
            },
        ),
        ("free", Dims::NONE, [0.0; 4]),
        ("density", Density::DIMS, [density; 4]),
        ("cp", SPECIFIC_HEAT_CAPACITY_DIMS, [cp; 4]),
        ("conductivity", THERMAL_CONDUCTIVITY_DIMS, [10.0; 4]),
    ] {
        let mut validity =
            ValidityDomain::unconstrained().with_quantity("T", quantity, 280.0, 340.0);
        let value = if pressure_pa.is_some() && name == "shift" {
            validity = validity.with_quantity("p", pressure_quantity, 1e5, 3e5);
            PropertyValue::Surface2D {
                axes: ["T".into(), "p".into()],
                axis_dims: [Temperature::DIMS, Pressure::DIMS],
                samples: [280.0, 300.0, 320.0, 340.0]
                    .into_iter()
                    .zip(values)
                    .flat_map(|(t, value)| [([t, 1e5], value), ([t, 3e5], 0.1 * value)])
                    .collect(),
                triangles: (0_u32..3)
                    .flat_map(|index| {
                        let left = 2 * index;
                        [[left, left + 2, left + 1], [left + 2, left + 3, left + 1]]
                    })
                    .collect(),
                dims,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [280.0, 300.0, 320.0, 340.0]
                    .into_iter()
                    .zip(values)
                    .collect(),
                dims,
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured inertial heat test".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured axial source".into(),
            phase: "solid".into(),
            process: "numerical test only".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let at = |t| match pressure_pa {
        Some(pressure) => temperature(t)
            .with_quantity("p", pressure_quantity, pressure)
            .unwrap(),
        None => temperature(t),
    };
    let (low, high) = (at(280.0), at(340.0));
    let branches = if branch {
        vec![(
            PropertyKey::new("E0", Pressure::DIMS),
            PropertyKey::new("tau0", Time::DIMS),
        )]
    } else {
        vec![]
    };
    let mechanical = ResolvedThermoMaxwell::try_new(
        &card,
        ThermoMaxwellSourceQuery {
            equilibrium_modulus: &PropertyKey::new("E_inf", Pressure::DIMS),
            branches: &branches,
            modulus_scale: &PropertyKey::new("scale", Dims::NONE),
            relaxation_scale: &PropertyKey::new("shift", Dims::NONE),
            free_strain: &PropertyKey::new("free", Dims::NONE),
            temperature_axis: "T",
            lower: &low,
            upper: &high,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        Limits {
            max_strain: 0.04,
            max_strain_rate_s: 10.0,
        },
        || false,
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("cp", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &low,
            upper: &high,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * cp,
        1e-8,
    )
    .unwrap();
    let transport = LumpedThermalTransport::from_material_card_without_radiation(
        &card,
        &PropertyKey::new("conductivity", THERMAL_CONDUCTIVITY_DIMS),
        "T",
        &[low, high],
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    Sources {
        mechanical,
        caloric,
        transport,
    }
}
fn geometry(payload: f64) -> Geometry {
    Geometry {
        length_m: 0.1,
        cross_section_m2: 1e-5,
        thermal_surface_m2: 0.0013,
        end_payload_kg: payload,
    }
}
fn config() -> Controls {
    Controls {
        feedback: FeedbackConfig {
            max_iterations: 100,
            temperature_tolerance_k: 1e-9,
            energy_tolerance_j: 1e-10,
            relaxation: 1.0,
        },
        force_tolerance_n: 1e-8,
        energy_tolerance_j: 2e-10,
    }
}
fn adiabatic(input: UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<()>, Infallible> {
    let h = input.initial.specific_enthalpy_j_kg() + input.internal_heat_j / input.mass_kg;
    Ok(UniformEnthalpyStep {
        state: input.curve.state_at_specific_enthalpy(h).unwrap(),
        external_heat_j: 0.0,
        energy_residual_tolerance_j: 1e-11,
        report: (),
    })
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_same_density_gives_galerkin_inertia_and_the_independent_maxwell_endpoint() {
    for rho in [1000.0, 2000.0] {
        let source = sources(1000.0, rho, true, false);
        let g = geometry(0.001);
        let old =
            Oscillator::try_new(&source.mechanical, &source.caloric, g, 300.0, 0.0, 0.2).unwrap();
        let volume = g.length_m * g.cross_section_m2;
        let mass = rho * volume / 3.0 + g.end_payload_kg;
        close(old.point().mass_kg(), rho * volume, 1e-16);
        close(old.effective_mass_kg(), mass, 1e-16);
        let (dt, force) = (0.0001, 0.03);
        let branch_factor = 1.0 / (1.0 + dt / 0.005);
        let tangent = 1e6 + 3e6 * branch_factor;
        let velocity = (mass * 0.2 + dt * force)
            / (mass + dt * dt * g.cross_section_m2 * tangent / g.length_m);
        let strain = dt * velocity / g.length_m;
        let branch_strain = branch_factor * strain;
        let heat = volume * (3e6 / 0.005) * branch_strain * branch_strain * dt;
        let material_numerical =
            volume * 0.5 * (1e6 * strain * strain + 3e6 * branch_strain * branch_strain);
        let inertial_numerical = 0.5 * mass * (velocity - 0.2).powi(2);
        let step = old
            .step(force, dt, config(), adiabatic, &CancelGate::new())
            .unwrap();
        close(step.next.velocity_m_s(), velocity, 1e-13);
        close(
            step.next.point().mechanical().branch_strains()[0],
            branch_strain,
            1e-15,
        );
        close(step.material_heat_j, heat, 1e-14);
        close(
            step.numerical_loss_j,
            material_numerical + inertial_numerical,
            1e-14,
        );
        close(step.external_work_j, force * dt * velocity, 1e-14);
        close(
            old.point().mass_kg()
                * (step.next.point().caloric_state().specific_enthalpy_j_kg()
                    - old.point().caloric_state().specific_enthalpy_j_kg()),
            heat,
            1e-12,
        );
        assert!(step.combined_residual_j.abs() < 1e-12);
        assert_eq!(old.point().time_s(), 0.0);
        assert_eq!(old.point().mechanical().branch_strains(), &[0.0]);
    }
}

#[test]
fn g1_elastic_mode_converges_and_never_turns_backward_euler_loss_into_heat() {
    let source = sources(10.0, 1000.0, false, false);
    let g = geometry(0.0);
    let initial =
        Oscillator::try_new(&source.mechanical, &source.caloric, g, 300.0, 0.0, 0.2).unwrap();
    let omega = (g.cross_section_m2 * 1e6 / g.length_m / initial.effective_mass_kg()).sqrt();
    let end = std::f64::consts::FRAC_PI_2 / omega;
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut state = initial.clone();
        let mut numerical = 0.0;
        for _ in 0..count {
            let step = state
                .step(
                    0.0,
                    end / count as f64,
                    config(),
                    adiabatic,
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(step.material_heat_j, 0.0);
            assert_eq!(step.applied_material_heat_j, 0.0);
            assert_eq!(
                step.next.point().caloric_state(),
                initial.point().caloric_state()
            );
            numerical += step.numerical_loss_j;
            state = step.next;
        }
        let energy = state.kinetic_energy_j()
            + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3();
        assert!(numerical > 0.0);
        close(energy + numerical, initial.kinetic_energy_j(), 1e-12);
        errors.push((state.extension_m() - 0.2 / omega).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            (1.8..2.2).contains(&(pair[0] / pair[1])),
            "first-order inertial scheme: {errors:?}"
        );
    }
}

#[test]
fn g3_caloric_data_changes_actual_vibration_through_the_existing_thermal_owner() {
    let low = sources(0.02, 1000.0, true, true);
    let high = sources(2.0, 1000.0, true, true);
    assert_eq!(
        low.mechanical.law(),
        high.mechanical.law(),
        "only caloric data changed"
    );
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x544d44,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let environment = |source: &Sources| LumpedThermalEnvironment {
            temperature_k: 300.0,
            radiation_temperature_k: 300.0,
            convection_w_per_m2_k: 0.0,
            transport: source.transport.clone(),
            enthalpy_tolerance_j_kg: 1e-12,
            maximum_thermal_residual_j: 1e-11,
        };
        let mut a = Oscillator::try_new(
            &low.mechanical,
            &low.caloric,
            geometry(0.0),
            300.0,
            0.0,
            0.5,
        )
        .unwrap();
        let mut b = Oscillator::try_new(
            &high.mechanical,
            &high.caloric,
            geometry(0.0),
            300.0,
            0.0,
            0.5,
        )
        .unwrap();
        let initial_energy = a.kinetic_energy_j();
        let initial_enthalpies = [
            a.point().caloric_state().specific_enthalpy_j_kg(),
            b.point().caloric_state().specific_enthalpy_j_kg(),
        ];
        let (ea, eb) = (environment(&low), environment(&high));
        let mut numerical = [0.0; 2];
        let mut max_velocity_difference = 0.0_f64;
        for index in 0..800 {
            let sa = a
                .step(
                    0.0,
                    0.00005,
                    config(),
                    |input| ea.advance(&cx, input),
                    &gate,
                )
                .unwrap();
            let sb = b
                .step(
                    0.0,
                    0.00005,
                    config(),
                    |input| eb.advance(&cx, input),
                    &gate,
                )
                .unwrap();
            if index == 0 {
                // The returned inertial solution uses the original branch
                // owner at its ACTUAL thermal endpoint, not lagged damping.
                let expected = a
                    .point()
                    .mechanical()
                    .advance(
                        sa.next.point().caloric_state().temperature_k(),
                        sa.next.point().mechanical().total_strain(),
                        0.00005,
                        config().feedback.energy_tolerance_j / a.point().volume_m3(),
                    )
                    .unwrap();
                assert_eq!(sa.next.point().mechanical(), &expected.point);
            }
            assert!(sa.material_heat_j >= 0.0 && sb.material_heat_j >= 0.0);
            numerical[0] += sa.numerical_loss_j;
            numerical[1] += sb.numerical_loss_j;
            a = sa.next;
            b = sb.next;
            max_velocity_difference =
                max_velocity_difference.max((a.velocity_m_s() - b.velocity_m_s()).abs());
        }
        assert!(
            a.point().caloric_state().temperature_k()
                > b.point().caloric_state().temperature_k() + 0.01
        );
        assert!(max_velocity_difference > 1e-6);
        for (i, state) in [&a, &b].into_iter().enumerate() {
            let energy = state.kinetic_energy_j()
                + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3()
                + state.point().mass_kg()
                    * (state.point().caloric_state().specific_enthalpy_j_kg()
                        - initial_enthalpies[i]);
            close(energy + numerical[i], initial_energy, 1e-8);
            assert_eq!(
                state.point().mechanical().temperature_k(),
                state.point().caloric_state().temperature_k()
            );
        }
    });
}

#[test]
fn g3_pressure_surface_changes_motion_and_physical_heat_through_the_existing_owners() {
    let sources = [
        sources_at_pressure(0.02, 1000.0, true, true, Some(1.5e5)),
        sources_at_pressure(0.02, 1000.0, true, true, Some(2.5e5)),
    ];
    assert_eq!(
        sources[0].mechanical.samples()[0].card_identity(),
        sources[1].mechanical.samples()[0].card_identity(),
        "only the fixed pressure query changes"
    );
    assert_ne!(sources[0].mechanical.law(), sources[1].mechanical.law());
    // Each source-cell crossing, including the one inside the active 300–320 K
    // interval, must survive lowering into the actual hereditary law.
    for (source, crossed_temperature) in sources.iter().zip([315.0, 305.0]) {
        assert!(source.mechanical.samples().iter().any(|sample| {
            sample
                .query_point()
                .contains(&("T".into(), crossed_temperature))
        }));
    }
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x544d50,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let environments = sources.each_ref().map(|source| LumpedThermalEnvironment {
            temperature_k: 300.0,
            radiation_temperature_k: 300.0,
            convection_w_per_m2_k: 0.0,
            transport: source.transport.clone(),
            enthalpy_tolerance_j_kg: 1e-12,
            maximum_thermal_residual_j: 1e-11,
        });
        let mut states = sources.each_ref().map(|source| {
            Oscillator::try_new(
                &source.mechanical,
                &source.caloric,
                geometry(0.0),
                300.0,
                0.0,
                0.5,
            )
            .unwrap()
        });
        let initial_energy = states[0].kinetic_energy_j();
        let initial_enthalpies = states
            .each_ref()
            .map(|state| state.point().caloric_state().specific_enthalpy_j_kg());
        let mut numerical = [0.0; 2];
        let mut heat = [0.0; 2];
        let mut max_velocity_difference = 0.0_f64;
        let mut max_temperature_difference = 0.0_f64;
        for _ in 0..200 {
            for index in 0..2 {
                let old = &states[index];
                let step = old
                    .step(
                        0.0,
                        0.00005,
                        config(),
                        |input| environments[index].advance(&cx, input),
                        &gate,
                    )
                    .unwrap();
                let expected = old
                    .point()
                    .mechanical()
                    .advance(
                        step.next.point().caloric_state().temperature_k(),
                        step.next.point().mechanical().total_strain(),
                        0.00005,
                        config().feedback.energy_tolerance_j / old.point().volume_m3(),
                    )
                    .unwrap();
                assert_eq!(step.next.point().mechanical(), &expected.point);
                assert!(step.material_heat_j >= 0.0);
                assert!(step.combined_residual_j.abs() <= config().energy_tolerance_j);
                heat[index] += step.material_heat_j;
                numerical[index] += step.numerical_loss_j;
                states[index] = step.next;
            }
            max_velocity_difference = max_velocity_difference
                .max((states[0].velocity_m_s() - states[1].velocity_m_s()).abs());
            max_temperature_difference = max_temperature_difference.max(
                (states[0].point().caloric_state().temperature_k()
                    - states[1].point().caloric_state().temperature_k())
                .abs(),
            );
        }
        assert!(max_velocity_difference > 1e-4);
        assert!(max_temperature_difference > 1e-4);
        for (index, state) in states.iter().enumerate() {
            let caloric_gain = state.point().mass_kg()
                * (state.point().caloric_state().specific_enthalpy_j_kg()
                    - initial_enthalpies[index]);
            close(caloric_gain, heat[index], 1e-9);
            let energy = state.kinetic_energy_j()
                + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3()
                + caloric_gain;
            close(energy + numerical[index], initial_energy, 1e-8);
            assert_eq!(
                state.point().mechanical().temperature_k(),
                state.point().caloric_state().temperature_k()
            );
        }
    });
}

#[test]
fn g4_failed_thermal_trials_and_cancellation_preserve_the_whole_vibrating_history() {
    let source = sources(0.02, 1000.0, true, true);
    let initial = Oscillator::try_new(
        &source.mechanical,
        &source.caloric,
        geometry(0.0),
        300.0,
        0.0,
        0.5,
    )
    .unwrap();
    let gate = CancelGate::new();
    let old = initial
        .step(0.0, 0.0001, config(), adiabatic, &gate)
        .unwrap()
        .next;
    let saved = old.clone();
    let result = old.step(
        0.0,
        0.0001,
        config(),
        |_: UniformEnthalpyStepInput<'_>| {
            Err::<UniformEnthalpyStep<()>, _>("finite reservoir refused")
        },
        &gate,
    );
    assert!(matches!(
        result,
        Err(Error::Feedback(FeedbackError::Transport(
            "finite reservoir refused"
        )))
    ));
    let cancelled = CancelGate::new();
    let result = old.step(
        0.0,
        0.0001,
        config(),
        |input| {
            let candidate = adiabatic(input);
            cancelled.request();
            candidate
        },
        &cancelled,
    );
    assert!(matches!(
        result,
        Err(Error::Feedback(FeedbackError::Cancelled))
    ));
    assert_eq!(old.point(), saved.point());
    assert_eq!(old.velocity_m_s().to_bits(), saved.velocity_m_s().to_bits());
    let mut a = old;
    let mut b = saved;
    for force in [0.03, -0.02, 0.0, 0.01, 0.0] {
        let sa = a.step(force, 0.0001, config(), adiabatic, &gate).unwrap();
        let sb = b.step(force, 0.0001, config(), adiabatic, &gate).unwrap();
        assert_eq!(sa.external_work_j, sb.external_work_j);
        assert_eq!(sa.material_heat_j, sb.material_heat_j);
        a = sa.next;
        b = sb.next;
        assert_eq!(a.point(), b.point());
        assert_eq!(a.velocity_m_s().to_bits(), b.velocity_m_s().to_bits());
    }
    let other = sources(1.0, 1000.0, true, true);
    assert!(
        Oscillator::try_new(
            &source.mechanical,
            &other.caloric,
            geometry(0.0),
            300.0,
            0.0,
            0.0
        )
        .is_err()
    );
    assert!(
        Oscillator::try_new(
            &source.mechanical,
            &source.caloric,
            geometry(-1.0),
            300.0,
            0.0,
            0.0
        )
        .is_err()
    );
    assert!(
        Oscillator::try_new(
            &source.mechanical,
            &source.caloric,
            geometry(0.0),
            300.0,
            0.0,
            2.0
        )
        .is_err()
    );
    assert!(a.step(1e9, 0.0001, config(), adiabatic, &gate).is_err());
}

#[path = "thermovisco_dynamics/finite_gas.rs"]
mod finite_gas;
