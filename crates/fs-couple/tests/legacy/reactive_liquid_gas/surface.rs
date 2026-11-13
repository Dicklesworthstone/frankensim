//! G1/G3/G4: reacting finite gas, radiation and source surface stress share one step.
//! All coefficients are manufactured; these are coupling checks, not calibration.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceInterfaceFace, FreeSurfaceInterfaceSolve,
};
use fs_lbm::material::free_surface::enthalpy::thermocapillary::PlanarThermocapillarity;
use fs_material::liquid::surface::{
    ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

const DX: f64 = 0.001;

fn surface(surroundings_k: f64, emissivity: f64) -> FreeSurfaceInterfaceFace {
    FreeSurfaceInterfaceFace {
        face: face(),
        emissivity: emissivity.into(),
        surroundings_k,
    }
}

// Independent zero-storage Fourier/Stefan--Boltzmann balance in watts.
fn reference_surface(tl: f64, tg: f64, tr: f64, emissivity: f64, kg: f64) -> (f64, f64, f64) {
    let (mut lo, mut hi) = (tl.min(tg).min(tr), tl.max(tg).max(tr));
    for _ in 0..80 {
        let ts = 0.5 * (lo + hi);
        let balance = 0.5 * 1e-6 / 0.0005 * (tl - ts)
            + kg * 1e-6 / 0.0005 * (tg - ts)
            + emissivity * 5.670374419e-8 * 1e-6 * (tr.powi(4) - ts.powi(4));
        if balance > 0.0 {
            lo = ts;
        } else {
            hi = ts;
        }
    }
    let ts = 0.5 * (lo + hi);
    (
        ts,
        kg * 1e-6 / 0.0005 * (tg - ts),
        emissivity * 5.670374419e-8 * 1e-6 * (tr.powi(4) - ts.powi(4)),
    )
}

#[test]
fn g1_reacting_surface_closes_independent_radiation_and_finite_capacity_balance() {
    for (surroundings, emissivity) in [(600.0, 0.8), (350.0, 0.8), (400.0, 0.0)] {
        with_liquid(false, |liquid, gate| {
            let law = reaction(500.0, 20000.0, 8.0);
            let original = gas(&law);
            let next = law
                .liquid_surface_candidate(
                    liquid,
                    &original,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &[surface(surroundings, emissivity)],
                    FreeSurfaceInterfaceSolve::default(),
                    controls(),
                    gate,
                )
                .unwrap();
            let extent = DT * VOLUME * 500.0;
            let amount = 100000.0 * VOLUME / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0);
            let capacity = 2.5 * 100000.0 * VOLUME / 400.0;
            let kg = 0.03 * (1.0 + 7.0 * extent / amount);
            let (mut tl, mut tg) = (400.0, 400.0 + 20000.0 * extent / capacity);
            let (mut gas_heat, mut radiant_heat) = (0.0, 0.0);
            for _ in 0..4 {
                let (_, qg, qr) = reference_surface(tl, tg, surroundings, emissivity, kg);
                tl += DT / 4.0 * (qg + qr) / 0.002;
                tg -= DT / 4.0 * qg / capacity;
                gas_heat += DT / 4.0 * qg;
                radiant_heat += DT / 4.0 * qr;
            }
            let report = next.report();
            let exchange = report.thermal.conduction.surface.as_ref().unwrap();
            assert!(report.thermal.conduction.environment.is_none());
            assert!(report.thermal.conduction.finite_gas.is_some());
            assert_eq!(report.thermal.conduction.thermal.substeps, 4);
            let endpoint = &exchange.thermal_endpoint[0];
            near(endpoint.liquid_temperature_k, tl, 2e-9);
            near(
                endpoint.surface_temperature_k,
                reference_surface(tl, tg, surroundings, emissivity, kg).0,
                2e-8,
            );
            near(next.gas().state().state().temperature, tg, 2e-8);
            near(
                next.gas().state().state().pressure,
                100000.0 * tg / 400.0,
                1e-5,
            );
            near(exchange.heat_from_gas_j, gas_heat, 2e-12);
            near(exchange.heat_from_radiation_j, radiant_heat, 2e-12);
            near(report.gas_internal_energy_change_j, -gas_heat, 2e-12);
            near(
                report.liquid_enthalpy_change_j,
                gas_heat + radiant_heat,
                2e-11,
            );
            near(report.supplied_thermal_energy_j, radiant_heat, 2e-12);
            assert!(report.combined_energy_residual_j.abs() < 2e-9);
            near(report.reaction.extent_mol, extent, 1e-20);
            if emissivity == 0.0 {
                let old = law
                    .liquid_candidate(
                        liquid,
                        &original,
                        &[],
                        &[],
                        &[],
                        &[],
                        &[],
                        &[face()],
                        controls(),
                        gate,
                    )
                    .unwrap();
                near(
                    next.gas().state().state().temperature,
                    old.gas().state().state().temperature,
                    2e-8,
                );
                near(
                    report.liquid_enthalpy_change_j,
                    old.report().liquid_enthalpy_change_j,
                    2e-11,
                );
            }
        });
    }
}

fn planar_surfaces() -> Vec<FreeSurfaceInterfaceFace> {
    (0..8)
        .flat_map(|y| {
            (0..8).map(move |x| FreeSurfaceInterfaceFace {
                face: FreeSurfaceGasFace {
                    liquid: [x, y, 3],
                    gas: [x, y, 4],
                    area_m2: DX * DX,
                    distances_m: [0.25 * DX, 0.5 * DX],
                },
                emissivity: 0.8.into(),
                surroundings_k: 350.0 + 20.0 * x as f64,
            })
        })
        .collect()
}

fn with_planar<R>(
    slope: f64,
    run: impl FnOnce(&MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let base = liquid_card();
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.01 - 100.0 * slope), (500.0, 0.01 + 100.0 * slope)],
                dims: SURFACE_TENSION_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 500.0)
                .with("pressure", 100000.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured surface-tension law".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let point = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap();
    let mut requirements = [
        ("density", Dims([-3, 1, 0, 0, 0, 0])),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0])),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect::<Vec<_>>();
    requirements.push(
        ScalarPropertyRequirement::try_with_quantity(
            "thermal_conductivity",
            quantity(QuantityKind::ThermalConductivity),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = ResolvedLiquidSurfaceState::try_from_material(
        material,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    )
    .unwrap();
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        [300.0, 600.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: (temperature_k - 300.0) * 2000.0,
                temperature_k,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let cells: Vec<_> = (0..4)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (0..8).map(move |x| MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if z == 3 { 0.5 } else { 1.0 },
                    velocity_m_s: [0.0; 3],
                })
            })
        })
        .collect();
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let liquid = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_thermocapillarity(
                MaterialFreeSurfaceDomain3 {
                    cells: [8; 3],
                    gravity_m_s2: [0.0; 3],
                    contact: ContactModel3::Neutral,
                    solid_tiles: &[],
                    max_cells: 512,
                },
                &cells,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                PlanarThermocapillarity {
                    normal_axis: 2,
                    outward_sign: 1,
                    surface_position_m: 3.5 * DX,
                    maximum_departure_m: 0.01 * DX,
                    minimum_surface_fill: 0.25,
                },
                None,
                None,
                |_| 200000.0,
                &cx,
                &gate,
            )
            .unwrap();
        run(&liquid, &gate)
    })
}

#[test]
fn g3_chemical_and_surface_data_reach_actual_native_pressure_and_marangoni_motion() {
    let run = |rate, slope| {
        with_planar(slope, |liquid, gate| {
            let law = reaction(rate, 20000.0, 8.0);
            let original = gas(&law);
            let mut limits = controls();
            limits.thermal.maximum_faces = 64;
            let next = law
                .liquid_surface_candidate(
                    liquid,
                    &original,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &planar_surfaces(),
                    FreeSurfaceInterfaceSolve::default(),
                    limits,
                    gate,
                )
                .unwrap();
            let report = next.report();
            let exchange = report.thermal.conduction.surface.as_ref().unwrap();
            assert_eq!(exchange.native_endpoint.len(), 64);
            let sample = |x| {
                exchange
                    .native_endpoint
                    .iter()
                    .find(|s| s.surface.face.liquid == [x, 3, 3])
                    .unwrap()
            };
            let sigma = |x| {
                let state = sample(x);
                let value = 0.01 + slope * (state.surface_temperature_k - 400.0);
                near(
                    state.surface_tension.as_ref().unwrap().value_si(),
                    value,
                    2e-12,
                );
                value
            };
            let native = next.liquid().flow().lattice();
            let fill = native.sample_at(3, 3, 3).unwrap().fill_fraction;
            let expected_force =
                (sigma(4) - sigma(2)) / (2.0 * DX) / (fill * DX) / 1000.0 * DT * DT / DX;
            near(
                native.force_density_at(3, 3, 3).unwrap()[0],
                expected_force,
                2e-10,
            );
            near(
                native.surface_tension_at(3, 3, 3).unwrap(),
                sigma(3) * DT * DT / (1000.0 * DX.powi(3)),
                2e-12,
            );
            let speed = next
                .liquid()
                .flow()
                .sample_at_si(3, 3, 3)
                .unwrap()
                .unwrap()
                .velocity_m_s[0];
            if slope < 0.0 {
                assert!(speed < 0.0);
            } else {
                near(speed, 0.0, 1e-15);
            }
            near(
                report.supplied_thermal_energy_j,
                exchange.heat_from_radiation_j,
                2e-11,
            );
            assert!(report.combined_energy_residual_j.abs() < limits.total_energy_tolerance_j);
            (
                native.grid().state_bits(),
                sample(3).surface_temperature_k.to_bits(),
                speed.to_bits(),
            )
        })
    };
    let base = run(500.0, -0.00001);
    assert_eq!(base, run(500.0, -0.00001));
    let faster = run(1000.0, -0.00001);
    assert_ne!(base.0, faster.0);
    assert_ne!(base.1, faster.1);
    assert_ne!(base.0, run(500.0, 0.0).0);
}

#[test]
fn g4_reacting_surface_continuation_and_late_refusal_keep_original_owners() {
    with_liquid(false, |liquid, gate| {
        let law = reaction(500.0, 20000.0, 4.0);
        let original = gas(&law);
        let boundary = [surface(600.0, 0.8)];
        let first = law
            .liquid_surface_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &boundary,
                FreeSurfaceInterfaceSolve::default(),
                controls(),
                gate,
            )
            .unwrap();
        let expected = law
            .parcel_candidate(
                &ClosedGasParcel::from_isochoric(first.gas().clone()),
                DT,
                controls().reaction,
                || false,
            )
            .unwrap();
        let next = law
            .liquid_surface_candidate(
                first.liquid(),
                first.gas(),
                first.bodies(),
                &[],
                &[],
                &[],
                &[],
                &boundary,
                FreeSurfaceInterfaceSolve::default(),
                controls(),
                gate,
            )
            .unwrap();
        assert_eq!(next.report().reaction, *expected.step());
        assert!(next.report().reaction.species[1].initial_mol > 0.0);
        let bits = first.liquid().flow().lattice().grid().state_bits();
        let enthalpy = first.liquid().total_enthalpy_j().unwrap();
        let gas_before = first.gas().clone();
        let mut missing = boundary;
        missing[0].face.gas = [2, 3, 3];
        let mut too_tight = controls();
        too_tight.total_energy_tolerance_j = 0.0;
        for (surfaces, limits) in [(&missing[..], controls()), (&boundary[..], too_tight)] {
            assert!(
                law.liquid_surface_candidate(
                    first.liquid(),
                    first.gas(),
                    first.bodies(),
                    &[],
                    &[],
                    &[],
                    &[],
                    surfaces,
                    FreeSurfaceInterfaceSolve::default(),
                    limits,
                    gate
                )
                .is_err()
            );
            assert_eq!(first.liquid().flow().lattice().grid().state_bits(), bits);
            assert_eq!(first.liquid().total_enthalpy_j().unwrap(), enthalpy);
            assert_eq!(first.gas(), &gas_before);
        }
        gate.request();
        assert!(matches!(
            law.liquid_surface_candidate(
                first.liquid(),
                first.gas(),
                first.bodies(),
                &[],
                &[],
                &[],
                &[],
                &boundary,
                FreeSurfaceInterfaceSolve::default(),
                controls(),
                gate
            ),
            Err(ReactiveLiquidGasError::Cancelled)
        ));
        assert_eq!(first.liquid().flow().lattice().grid().state_bits(), bits);
        assert_eq!(first.gas(), &gas_before);
    });
}
