//! Independent NASA-polynomial and finite-volume references for gas thermosiphons.
//! Synthetic coefficients exercise implementation, not calibrated gas accuracy.
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::transient::{
    FiniteGasLoop, FiniteLoopControls, FiniteLoopError, FiniteLoopLimits, GasLoopLimits, PipeCell,
};
use fs_airflow::graph::thermal::medium::series::hydraulic::{
    MaterialFlowResistance, pressure::SeriesGravity,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::transient::deposition::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
    backward_euler::{BoundaryStepCoefficients, StepCoefficients},
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::{
    SpeciesTransportMaterial, ThermochemicalGasModel,
    energy::{GasEnergySolve, enthalpy::GasEnthalpyLaw},
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Area, DynViscosity, Length, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, mixture::FrozenIdealGasMixtureModelV1,
    required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;
use std::convert::Infallible;

const AREA: f64 = 1e-3;
const PRESSURE: f64 = 140000.0;
const MOLAR_MASS: f64 = 0.75 * 0.028 + 0.25 * 0.004;
const VISCOUS_GEOMETRY: f64 = 5e7;
const SLOPE: f64 = 0.01;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), PRESSURE)
        .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic finite gas circuit".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn species(name: &str, mass: f64, cp_r: f64, slope: f64, datum: f64) -> Nasa9StandardStateModelV1 {
    let mut a = [0.0; 9];
    a[2] = cp_r;
    a[3] = slope;
    a[7] = mass * datum / R;
    let mut parameters = BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 250.0, 450.0),
        sources: vec![],
        provenance: provenance(),
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("fixture-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn transport(name: &str, mu: f64, k: f64) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        (
            PropertyKey::new("dynamic_viscosity", DynViscosity::DIMS),
            mu,
            DynViscosity::DIMS,
        ),
        (
            PropertyKey::with_quantity(
                "thermal_conductivity",
                typed(QuantityKind::ThermalConductivity),
            ),
            k,
            fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        250.0,
                        450.0,
                    )
                    .with_quantity(
                        "pressure",
                        typed(QuantityKind::Pressure),
                        100000.0,
                        200000.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "synthetic".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn gas_law() -> GasEnthalpyLaw {
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap();
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", 0.028, 3.5, SLOPE, -100000.0),
                species("B", 0.004, 2.5, 0.0, 300000.0),
            ],
            composition.clone(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap();
    GasEnthalpyLaw::from_mole(
        &model,
        composition,
        point(300.0),
        GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 360.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 64,
        },
    )
    .unwrap()
}
fn vertices() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
}
fn sections() -> Vec<PipeCell> {
    let law = MaterialFlowResistance::try_new(
        VISCOUS_GEOMETRY,
        0.0,
        Area::new(AREA),
        Length::new(0.03),
        (0.0, 2300.0),
    )
    .unwrap();
    (0..4)
        .map(|i| PipeCell {
            cell: CellId(i),
            resistance: law,
        })
        .collect()
}
fn gravity() -> SeriesGravity {
    SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap()
}
fn limits(maximum_mach: f64) -> GasLoopLimits {
    GasLoopLimits {
        inventory: FiniteLoopLimits {
            maximum_cells: 8,
            maximum_relative_density_departure: 0.1,
        },
        maximum_mach,
    }
}
fn controls() -> FiniteLoopControls {
    FiniteLoopControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_newton_iterations: 32,
        maximum_residual_evaluations: 2048,
        energy_scale_j: 1.0,
        flow_scale_m3_s: 1e-4,
        finite_difference_step: 1e-3,
        cell_tolerance_j: 1e-6,
        total_thermal_tolerance_j: 1e-6,
        momentum_tolerance_pa_s: 1e-6,
        mechanical_energy_tolerance_j: 1e-9,
        maximum_neglected_power_w: 1.0,
        maximum_neglected_energy_j: 1.0,
        maximum_drag_step: 2.0,
    }
}
fn heat(powers: [f64; 4]) -> BoundaryStepCoefficients<()> {
    StepCoefficients {
        sources_w: powers
            .into_iter()
            .enumerate()
            .map(|(i, p)| (CellId(i as u64), p))
            .collect(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into()
}
fn initial(law: &GasEnthalpyLaw) -> DepositionSimulation<'_> {
    let reference = law.model().evaluate(law.context()).unwrap();
    let mass = reference.state().density * AREA;
    let path = vertices();
    let h = law.enthalpy_at_temperature(300.0).unwrap();
    let schedule = (0..4)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: mass,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: path[i],
                    end_m: path[i + 1],
                    volume_m3: AREA,
                    surface_area_m2: 0.12,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links: 16,
            max_events: 8,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn context(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 95,
                kernel_id: 745,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

// Independent NASA fixture h(T)=datum+A*T+B*T^2, including both formation datums.
fn coefficients() -> (f64, f64, f64) {
    (
        (0.75 * 0.028 * (-100000.0) + 0.25 * 0.004 * 300000.0) / MOLAR_MASS,
        R / MOLAR_MASS * 3.25,
        R / MOLAR_MASS * 0.75 * SLOPE / 2.0,
    )
}
fn h(t: f64) -> f64 {
    let (c, a, b) = coefficients();
    c + a * t + b * t * t
}
fn t_from_h(h: f64) -> f64 {
    let (c, a, b) = coefficients();
    2.0 * (h - c) / (a + (a * a + 4.0 * b * (h - c)).sqrt())
}
fn finite_volume(mass: f64, rho0: f64, q: f64, powers: [f64; 4], dt: f64) -> [f64; 4] {
    let c = rho0 * q.abs() * dt / mass;
    let ratio = c / (1.0 + c);
    let b: [f64; 4] = core::array::from_fn(|i| (h(300.0) + powers[i] * dt / mass) / (1.0 + c));
    core::array::from_fn(|i| {
        let mut sum = 0.0;
        let mut weight = 1.0;
        for k in 0..4 {
            let upstream = if q >= 0.0 {
                (i + 4 - k) % 4
            } else {
                (i + k) % 4
            };
            sum += weight * b[upstream];
            weight *= ratio;
        }
        sum / (1.0 - ratio.powi(4))
    })
}

#[test]
fn nasa_enthalpy_inventory_and_donor_transport_drive_real_inertial_gas_buoyancy() {
    context(|cx| {
        let source = gas_law();
        let mut thermal = initial(&source);
        let mut owner = FiniteGasLoop::try_new(
            cx,
            &source,
            &thermal,
            &sections(),
            gravity(),
            0.0,
            limits(0.1),
        )
        .unwrap();
        assert!(core::ptr::eq(owner.source(), &source));
        let rho0 = PRESSURE * MOLAR_MASS / (R * 300.0);
        let mass = rho0 * AREA;
        let inertia = 4.0 * rho0 / AREA;
        let reference = source.model().evaluate(source.context()).unwrap();
        // Pure transport fixtures are constant; original Wilke mixture is fixed.
        let resistance = 4.0 * VISCOUS_GEOMETRY * reference.state().dynamic_viscosity;
        let powers = [5.0, 10.0, 100.0, 15.0];
        let dt = 0.5;
        let residual = |q: f64| {
            let hs = finite_volume(mass, rho0, q, powers, dt);
            let ts = hs.map(t_from_h);
            let head = 9.81 * PRESSURE * MOLAR_MASS / R * (1.0 / ts[0] - 1.0 / ts[2]);
            q - dt * head / (2.0 * inertia + dt * resistance)
        };
        let (mut lo, mut hi) = (0.0, 1e-3);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            if residual(q) > 0.0 {
                hi = q;
            } else {
                lo = q;
            }
        }
        let q = 0.5 * (lo + hi);
        let expected = finite_volume(mass, rho0, q, powers, dt);
        let report = owner
            .advance(
                cx,
                &mut thermal,
                dt,
                0.0,
                |_, _| Ok::<_, Infallible>(heat(powers)),
                controls(),
                &[0.0; 4],
            )
            .unwrap();
        assert!((report.hydraulics.average_volume_flow_m3_s - q).abs() < 1e-10);
        assert!((owner.volume_flow_m3_s() - 2.0 * q).abs() < 2e-10);
        let mut measured = 0.0;
        let mut endpoint_cp_shortcut = 0.0;
        for i in 0..4 {
            let cell = &thermal.cells()[&CellId(i as u64)];
            let temperature = thermal.temperature_k(cell.id).unwrap();
            assert!((temperature - t_from_h(expected[i])).abs() < 2e-6);
            assert!((cell.specific_h_j_kg - expected[i]).abs() < 2e-3);
            measured += cell.mass_kg * (cell.specific_h_j_kg - h(300.0));
            let state = &report.hydraulics.sources[i];
            assert_eq!(state.state().pressure, PRESSURE);
            assert_eq!(
                state.thermodynamics().receipt().declared_fraction_bits(),
                reference
                    .thermodynamics()
                    .receipt()
                    .declared_fraction_bits()
            );
            for (actual, original) in state.components().iter().zip(reference.components()) {
                assert_eq!(actual.card_identity(), original.card_identity());
            }
            let (_, a, b) = coefficients();
            endpoint_cp_shortcut += mass * (a + 2.0 * b * temperature) * (temperature - 300.0);
            let upstream = thermal.cells()[&CellId(((i + 3) % 4) as u64)].specific_h_j_kg;
            let defect = cell.mass_kg * (cell.specific_h_j_kg - h(300.0))
                - dt * (powers[i]
                    + rho0
                        * report.hydraulics.average_volume_flow_m3_s
                        * (upstream - cell.specific_h_j_kg));
            assert!(defect.abs() < 1e-6);
            assert_eq!(
                report.hydraulics.losses[i].transport.density_kg_m3(),
                owner.reference_density_kg_m3()
            );
        }
        assert!((measured - 65.0).abs() < 1e-6);
        assert!((report.thermal.stored_energy_change_j - 65.0).abs() < 1e-6);
        assert!((report.thermal.external_heat_j - 65.0).abs() < 1e-12);
        assert!((endpoint_cp_shortcut - measured).abs() > 0.01);
        assert!(report.hydraulics.sources[2].state().density < rho0);
        assert!(report.hydraulics.buoyancy_work_j > 0.0);
        assert!(report.hydraulics.mechanical_balance_residual_j.abs() < 1e-9);
        assert_eq!(owner.accepted_time_s(), thermal.time_s());
    });
}

#[test]
fn gas_reversal_checks_endpoint_mach_and_preserves_both_histories_on_refusal() {
    context(|cx| {
        let source = gas_law();
        let mut thermal = initial(&source);
        let reference = source.model().evaluate(source.context()).unwrap();
        let initial_mach = 1e-5;
        let initial_q = AREA * reference.state().sound_speed * initial_mach;
        let mut owner = FiniteGasLoop::try_new(
            cx,
            &source,
            &thermal,
            &sections(),
            gravity(),
            initial_q,
            limits(1.5 * initial_mach),
        )
        .unwrap();
        let before = thermal.cells().clone();
        let retained_flow = owner.volume_flow_m3_s();
        let inertia = owner.inertance_pa_s2_m3();
        let resistance = 4.0 * VISCOUS_GEOMETRY * reference.state().dynamic_viscosity;
        let dt = 0.1;
        let rejected_final = -2.0 * initial_q;
        let drive = inertia * (rejected_final - initial_q) / dt
            + resistance * 0.5 * (initial_q + rejected_final);
        let error = owner
            .advance(
                cx,
                &mut thermal,
                dt,
                drive,
                |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                controls(),
                &[0.0; 4],
            )
            .unwrap_err();
        match error {
            FiniteLoopError::GasMach { value, maximum, .. } => {
                assert!(value > maximum);
                assert!((value - 2.0 * initial_mach).abs() < 1e-10);
            }
            other => panic!("endpoint Mach must refuse despite a smaller average: {other:?}"),
        }
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(owner.accepted_time_s(), 0.0);
        assert_eq!(owner.volume_flow_m3_s(), retained_flow);

        let accepted_final = -1.25 * initial_q;
        let accepted_average = 0.5 * (initial_q + accepted_final);
        let drive = inertia * (accepted_final - initial_q) / dt + resistance * accepted_average;
        let report = owner
            .advance(
                cx,
                &mut thermal,
                dt,
                drive,
                |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                controls(),
                &[0.0; 4],
            )
            .unwrap();
        assert!(report.hydraulics.average_volume_flow_m3_s < 0.0);
        assert!((report.hydraulics.average_volume_flow_m3_s - accepted_average).abs() < 1e-10);
        assert!((owner.volume_flow_m3_s() - accepted_final).abs() < 1e-10);
        for (id, retained) in &before {
            let cell = &thermal.cells()[id];
            assert_eq!(cell.mass_kg, retained.mass_kg);
            assert_eq!(cell.geometry, retained.geometry);
            let storage_change =
                cell.mass_kg * (cell.specific_h_j_kg - retained.specific_h_j_kg);
            assert!(storage_change.abs() < controls().cell_tolerance_j);
        }
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert_eq!(thermal.time_s(), dt);
    });
}
