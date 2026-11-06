//! G1/G3/G4: original spatial polarization, switching heat and finite enthalpy.
//! Scaled analytical SI fixtures exercise the equations, not measured materials.
#![cfg(feature = "experimental-dielectric-hysteresis")]

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        dielectric::hysteretic::{
            DielectricRegion, HystereticDielectricThermalError, ThermalHystereticElectrostatics,
        },
        material::{ConductionLink, MaterialCell, MaterialNetwork},
    },
    material_circuit::MaterialCircuitSolve,
};
use fs_em::electrostatic::nonlinear::{
    DielectricMaterial, ElectrodeControl, NonlinearElectrostaticControls,
    hysteretic::{
        PlanarHystereticElectrostatics,
        path::{HystereticPathControls, HystereticPathError},
    },
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{
            hysteresis::{
                source::{
                    DielectricHysteresisBranchSources, DielectricHysteresisSources,
                    ResolvedDielectricHysteresis,
                },
                vector::{ResolvedVectorDielectricHysteresis, VectorDielectricHysteresisControls},
            },
            reversible::{ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS},
        },
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};
use fs_thermal::enthalpy::*;
use std::sync::Arc;

const DEPTH: f64 = 0.3;
const EPS: f64 = 0.01;
const CHI: f64 = 0.2;
const EC: f64 = 2.0;
const PS: f64 = 1.0;
const DENSITY: f64 = 2.0;
const CP: f64 = 3.0;
const CAPACITY: f64 = DEPTH * DENSITY * CP;

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-8 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn roles(switching: bool) -> DielectricHysteresisSources {
    DielectricHysteresisSources {
        background_permittivity: PropertyKey::new("epsilon_background", ELECTRIC_PERMITTIVITY_DIMS),
        field_limit: PropertyKey::new("field_limit", ELECTRIC_FIELD_DIMS),
        branches: if switching {
            vec![DielectricHysteresisBranchSources {
                susceptibility: PropertyKey::new("switching_slope", ELECTRIC_PERMITTIVITY_DIMS),
                pinning_field: PropertyKey::new("pinning", ELECTRIC_FIELD_DIMS),
                saturation_polarization: PropertyKey::new(
                    "polarization_limit",
                    ELECTRIC_DISPLACEMENT_DIMS,
                ),
            }]
        } else {
            vec![]
        },
    }
}
fn conductivity() -> PropertyKey {
    PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS)
}
fn card(explicit_temperature: bool) -> MaterialCard {
    let domain = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        500.0,
    );
    let r = roles(true);
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value, validity, interpolation| {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "analytical spatial dielectric/enthalpy fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (key, value) in [
        (r.background_permittivity, EPS),
        (r.field_limit, 30.0),
        (r.branches[0].susceptibility.clone(), CHI),
        (r.branches[0].pinning_field.clone(), EC),
        (r.branches[0].saturation_polarization.clone(), PS),
    ] {
        let dims = key.dims();
        insert(
            key,
            PropertyValue::Scalar { value, dims },
            if explicit_temperature {
                domain.clone()
            } else {
                ValidityDomain::unconstrained()
            },
            InterpolationPolicy::ConstantWithinValidity,
        );
    }
    for (key, value) in [
        (PropertyKey::new("density", Density::DIMS), DENSITY),
        (
            PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            CP,
        ),
    ] {
        let dims = key.dims();
        insert(
            key,
            PropertyValue::Scalar { value, dims },
            domain.clone(),
            InterpolationPolicy::ConstantWithinValidity,
        );
    }
    insert(
        conductivity(),
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 1.0), (500.0, 4.0)],
            dims: THERMAL_CONDUCTIVITY_DIMS,
        },
        domain,
        InterpolationPolicy::LinearInside,
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "scaled spatial switching and heat fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn local_controls() -> VectorDielectricHysteresisControls {
    VectorDielectricHysteresisControls {
        relative_tolerance: 1e-12,
        absolute_field_tolerance_v_m: 1e-12,
        absolute_displacement_tolerance_c_m2: 1e-13,
        work_tolerance_j_m3: 1e-10,
        maximum_branch_iterations: 96,
        maximum_inverse_iterations: 96,
        maximum_backtracks: 32,
    }
}
fn field_controls() -> NonlinearElectrostaticControls {
    NonlinearElectrostaticControls {
        absolute_residual_tolerance_c: 1e-11,
        relative_residual_tolerance: 1e-12,
        charge_tolerance_c: 1e-10,
        energy_tolerance_j: 2e-9,
        max_newton_iterations: 48,
        max_backtracks: 32,
        max_operator_applications: 2048,
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 128,
            max_corrections: 2,
        },
    }
}
fn heat_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 2e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-13,
            linear_restart: 4,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn field(
    card: &MaterialCard,
    switching: bool,
    t: f64,
    rotation: [f64; 2],
) -> PlanarHystereticElectrostatics {
    let source = ResolvedVectorDielectricHysteresis::from_isotropic_source(
        ResolvedDielectricHysteresis::resolve(
            Arc::new(card.clone()),
            &point(t),
            roles(switching),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
    .unwrap();
    let [c, s] = rotation;
    let vertices = [
        [0.0, 0.0],
        [0.5, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 1.0],
        [1.0, 1.0],
    ];
    PlanarHystereticElectrostatics::new(
        TriComplex2::from_indexed_triangles(
            "dielectric finite thermal slab",
            vertices
                .map(|[x, y]| [c * x - s * y, s * x + c * y])
                .to_vec(),
            vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
            Metric2::planar(DEPTH).unwrap(),
        )
        .unwrap(),
        vec![source; 4],
        &[&[0, 3], &[2, 5]],
        local_controls(),
        WeightedAssemblyLimits {
            max_cells: 32,
            max_dofs: 64,
            max_triplets: 288,
        },
        &mut || false,
    )
    .unwrap()
}
fn voltage(v: f64) -> [ElectrodeControl; 2] {
    [ElectrodeControl::Voltage(v), ElectrodeControl::Voltage(0.0)]
}
fn thermal(
    curve: &ResolvedHeatCapacityCurve,
    count: u64,
    linked: bool,
    t: f64,
) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = (0..count)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i),
                segment_id: i,
                body_id: i,
                material_id: MaterialId(0),
                mass_kg: DENSITY * DEPTH,
                specific_h_j_kg: CP * (t - 200.0),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: DEPTH,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: t,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(0),
                b: CellId(1),
                area_m2: DEPTH,
                conductance_w_k: 999.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
    count: usize,
    conductivity: Option<&'a PropertyKey>,
) -> Vec<MaterialCell<'a>> {
    (0..count)
        .map(|i| MaterialCell {
            cell: CellId(i as u64),
            material_id: MaterialId(0),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: conductivity.map(|key| (key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect()
}
fn regions() -> Vec<DielectricRegion> {
    (0..4)
        .map(|triangle| DielectricRegion::Cell {
            triangle,
            cell: CellId(0),
        })
        .collect()
}
fn play(old: f64, electric: f64) -> f64 {
    old.clamp(CHI * (electric - EC), CHI * (electric + EC))
        .clamp(-PS, PS)
}

#[test]
fn shared_receiver_minor_loops_match_independent_play_and_only_physical_heat() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let mut thermal = thermal(&curve, 1, false, 300.0);
    let mut coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [0.6, 0.8]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let (mut p, mut energy, mut d, mut expected_t) = (0.0, 0.0, 0.0, 300.0);
    for electric in [12.0, 0.0, 6.0, 1.0, 5.5, 1.0, 6.0, 0.0, -12.0, 0.0] {
        let next_p = play(p, electric);
        let next_d = EPS * electric + next_p;
        let next_energy = DEPTH * (0.5 * EPS * electric * electric + next_p * next_p / (2.0 * CHI));
        let heat = DEPTH * EC * (next_p - p).abs();
        let work = DEPTH * electric * (next_d - d);
        expected_t += heat / CAPACITY;
        let next = coupled
            .step_candidate(
                &thermal,
                thermal.time_s() + 0.25,
                &[0.0; 4],
                &voltage(electric),
                heat_controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        near(next.deposited_heat_j(), heat);
        near(next.thermal().temperature_k(CellId(0)).unwrap(), expected_t);
        near(
            next.report().electrical.stored_change_j,
            next_energy - energy,
        );
        near(next.report().electrical.endpoint_source_work_j, work);
        near(
            next.report().electrical.numerical_dissipation_j,
            work - (next_energy - energy) - heat,
        );
        near(next.report().total_energy_residual_j, 0.0);
        assert!(next.report().circuit_trials > 1);
        for triangle in 0..4 {
            let field = next.field().field();
            let polarization = field.cell_history(triangle).unwrap().polarizations_c_m2()[0];
            near(polarization[0], 0.6 * next_p);
            near(polarization[1], 0.8 * next_p);
            let DielectricMaterial::Hysteretic(source) = &field.problem().materials()[triangle]
            else {
                panic!("missing history")
            };
            near(source.source().material().query().axes()["T"], expected_t);
        }
        assert_eq!(coupled.time_s(), thermal.time_s());
        (thermal, coupled, _, _) = next.into_parts();
        (p, energy, d) = (next_p, next_energy, next_d);
    }
}

#[test]
fn complete_heat_assembly_matches_independent_two_cell_temperature_dependent_conduction() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let key = conductivity();
    let bindings = cells(&card, &curve, &context, 2, Some(&key));
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    }];
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &links,
    };
    let thermal = thermal(&curve, 2, true, 300.0);
    let coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [1.0, 0.0]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let electrical = ElectricalContactNetwork {
        bulk: &network,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let dt = 0.4;
    let next = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 17,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat.prepare(&cx, &thermal).unwrap();
        coupled
            .step_candidate_with_heat_assembly(
                &cx,
                &prepared,
                &thermal,
                dt,
                &[0.0; 4],
                &voltage(12.0),
                heat_controls(),
                field_controls(),
                &[0.0; 2],
            )
            .unwrap()
    });
    // Saturated P=1 gives Q=volume*Ec. Total caloric energy fixes theta0+theta1;
    // the second cell obeys C*theta1 = dt*G(T0,T1)*(theta0-theta1).
    // G is the independently derived series conductance of two half-cell paths.
    let q = DEPTH * EC * PS;
    let total_rise = q / CAPACITY;
    let (mut low, mut high) = (0.0, total_rise / 2.0);
    for _ in 0..64 {
        let theta1 = 0.5 * (low + high);
        let theta0 = total_rise - theta1;
        let (k0, k1) = (2.0 + 0.01 * theta0, 2.0 + 0.01 * theta1);
        let conductance = DEPTH / (0.5 / k0 + 0.5 / k1);
        if CAPACITY * theta1 < dt * conductance * (theta0 - theta1) {
            low = theta1;
        } else {
            high = theta1;
        }
    }
    let theta1 = 0.5 * (low + high);
    near(
        next.thermal().temperature_k(CellId(0)).unwrap(),
        300.0 + total_rise - theta1,
    );
    near(
        next.thermal().temperature_k(CellId(1)).unwrap(),
        300.0 + theta1,
    );
    near(next.deposited_heat_j(), q);
    near(next.report().thermal.stored_energy_change_j, q);
    near(next.report().non_electrical_heat_j, 0.0);
    near(next.report().total_energy_residual_j, 0.0);
    assert!(next.report().electrical.numerical_dissipation_j > q);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        coupled
            .field()
            .cell_history(0)
            .unwrap()
            .polarizations_c_m2(),
        &[[0.0; 3]]
    );
    assert!(next.report().source.material.links[0].conductance_w_k > 0.6);
}

#[test]
fn fixed_zero_loss_background_is_valid_and_unmapped_switching_heat_refuses() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = thermal(&curve, 1, false, 300.0);
    let fixed: Vec<_> = (0..4)
        .map(|triangle| DielectricRegion::Fixed { triangle })
        .collect();
    let baseline = ThermalHystereticElectrostatics::bind(
        field(&card, false, 300.0, [1.0, 0.0]),
        &thermal,
        &network,
        &fixed,
        1e-12,
        || false,
    )
    .unwrap();
    let next = baseline
        .step_candidate(
            &thermal,
            1.0,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(next.deposited_heat_j(), 0.0);
    near(next.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
    assert!(next.report().electrical.numerical_dissipation_j > 0.0);

    let switching = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [1.0, 0.0]),
        &thermal,
        &network,
        &fixed,
        1e-12,
        || false,
    )
    .unwrap();
    assert!(
        switching
            .step_candidate(
                &thermal,
                1.0,
                &[0.0; 4],
                &voltage(12.0),
                heat_controls(),
                field_controls(),
                &[0.0],
                || false
            )
            .is_err()
    );
    assert_eq!(switching.time_s(), 0.0);
    assert_eq!(switching.field().terminal_charge_c(), &[0.0; 2]);
}

#[test]
fn physical_volume_and_original_temperature_coverage_are_required_before_heating() {
    for explicit in [false, true] {
        let card = card(explicit);
        let curve = caloric(&card);
        let context = point(300.0);
        let bindings = cells(&card, &curve, &context, 1, None);
        let network = MaterialNetwork {
            cells: &bindings,
            resistors: &[],
            links: &[],
        };
        let thermal = thermal(&curve, 1, false, 300.0);
        let mut mapping = regions();
        if explicit {
            mapping[3] = DielectricRegion::Fixed { triangle: 3 };
        }
        assert!(
            ThermalHystereticElectrostatics::bind(
                field(&card, true, 300.0, [1.0, 0.0]),
                &thermal,
                &network,
                &mapping,
                1e-12,
                || false
            )
            .is_err()
        );
    }
    let original = card(true);
    let foreign = card(false);
    let curve = caloric(&original);
    let context = point(300.0);
    let bindings = cells(&original, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = thermal(&curve, 1, false, 300.0);
    assert!(
        ThermalHystereticElectrostatics::bind(
            field(&foreign, true, 300.0, [1.0, 0.0]),
            &thermal,
            &network,
            &regions(),
            1e-12,
            || false
        )
        .is_err()
    );
}

#[test]
fn late_cancellation_budgets_and_caloric_domain_leave_both_histories_replayable() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = thermal(&curve, 1, false, 300.0);
    let coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [1.0, 0.0]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let old_cells = thermal.cells().clone();
    let old_field = coupled.field().cell_history(0).unwrap().clone();
    let mut polls = 0;
    let expected = coupled
        .step_candidate(
            &thermal,
            0.25,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            field_controls(),
            &[0.0],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut stopped = 0;
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.25,
                &[0.0; 4],
                &voltage(12.0),
                heat_controls(),
                field_controls(),
                &[0.0],
                || {
                    stopped += 1;
                    stopped >= polls - 1
                }
            )
            .is_err()
    );
    let mut budget = heat_controls();
    budget.maximum_circuit_trials = 1;
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.25,
                &[0.0; 4],
                &voltage(12.0),
                budget,
                field_controls(),
                &[0.0],
                || false
            )
            .is_err()
    );
    let mut field_budget = field_controls();
    field_budget.max_operator_applications = 0;
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.25,
                &[0.0; 4],
                &voltage(12.0),
                heat_controls(),
                field_budget,
                &[0.0],
                || false
            )
            .is_err()
    );
    let retried = coupled
        .step_candidate(
            &thermal,
            0.25,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(coupled.field().cell_history(0).unwrap(), &old_field);
    assert_eq!(coupled.time_s(), 0.0);
    assert_eq!(retried.thermal().cells(), expected.thermal().cells());
    assert_eq!(
        retried.report().electrical.cells,
        expected.report().electrical.cells
    );

    // This dyadic initial temperature and its enthalpy are exactly represented.
    // Q/C=1/3 K cannot fit in the remaining 1/4 K of admitted caloric support.
    let hot_context = point(499.75);
    let hot_bindings = cells(&card, &curve, &hot_context, 1, None);
    let hot_network = MaterialNetwork {
        cells: &hot_bindings,
        resistors: &[],
        links: &[],
    };
    let hot = self::thermal(&curve, 1, false, 499.75);
    let hot_coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 499.75, [1.0, 0.0]),
        &hot,
        &hot_network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    assert!(
        hot_coupled
            .step_candidate(
                &hot,
                0.25,
                &[0.0; 4],
                &voltage(12.0),
                heat_controls(),
                field_controls(),
                &[0.0],
                || false
            )
            .is_err()
    );
    assert_eq!(hot.time_s(), 0.0);
    assert_eq!(hot_coupled.field().terminal_charge_c(), &[0.0; 2]);
}

#[test]
fn refined_loading_deposits_cumulative_physical_heat_and_keeps_global_work_receipts() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = thermal(&curve, 1, false, 300.0);
    let coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [1.0, 0.0]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let path = HystereticPathControls {
        field: field_controls(),
        relative_tolerance: 1e-4,
        absolute_polarization_tolerance_c_m2: 1e-8,
        absolute_heat_tolerance_j: 1e-8,
        absolute_work_tolerance_j: 0.015,
        energy_tolerance_j: 2e-8,
        maximum_field_solves: 8192,
        maximum_operator_applications: 200_000,
        maximum_refinement_depth: 20,
    };
    let endpoint = coupled
        .step_candidate(
            &thermal,
            0.25,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    let refined = coupled
        .step_candidate(
            &thermal,
            0.25,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            path,
            &[0.0],
            || false,
        )
        .unwrap();
    // This monotone collinear path has exact total Q=Ec*volume*Ps independent
    // of subdivision. Refinement resolves electrical work; it cannot turn
    // the endpoint scheme's much larger numerical loss into extra heating.
    let q = DEPTH * EC * PS;
    near(refined.deposited_heat_j(), q);
    near(
        refined.thermal().temperature_k(CellId(0)).unwrap(),
        300.0 + q / CAPACITY,
    );
    near(refined.deposited_heat_j(), endpoint.deposited_heat_j());
    assert!(
        refined.report().electrical.numerical_dissipation_j
            < 0.1 * endpoint.report().electrical.numerical_dissipation_j
    );
    near(refined.report().total_energy_residual_j, 0.0);
    assert!(endpoint.path_summary().is_none());
    let summary = *refined.path_summary().unwrap();
    assert!(summary.field_solves > 3 && summary.accepted_substeps > 2);
    assert!(summary.rejected_intervals > 0);
    assert!(summary.field_solves <= path.maximum_field_solves);
    assert!(summary.operator_applications <= path.maximum_operator_applications);
    let (_, _, _, retained_summary) = refined.into_parts();
    assert_eq!(retained_summary, Some(summary));

    let mut exhausted = path;
    exhausted.maximum_field_solves = 3;
    assert!(matches!(
        coupled.step_candidate(
            &thermal,
            0.25,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            exhausted,
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::Source(
            HystereticDielectricThermalError::Path(HystereticPathError::FieldSolveBudget)
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(coupled.field().terminal_charge_c(), &[0.0; 2]);
}
