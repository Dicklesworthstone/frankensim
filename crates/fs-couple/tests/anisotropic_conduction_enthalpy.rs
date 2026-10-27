//! A rotated, temperature-dependent source tensor closes local heat and field.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::config;
use fs_blake3::ContentHash;
use fs_couple::circuit_enthalpy::conduction::circuit::{
    ConductionCircuit, ConductionCircuitTerminal,
};
use fs_couple::circuit_enthalpy::conduction::{ConductionRegion, ConductionThermalNetwork};
use fs_couple::circuit_enthalpy::material::MaterialCell;
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::conductor::planar::{
    PRINCIPAL_RESISTIVITY_PROPERTIES, PlanarElectricalFrame, ResolvedOrthotropicConductor2,
    resolve_orthotropic_conductor2,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use fs_thermal::enthalpy::*;

const VOLTAGES: [(usize, f64); 4] = [(0, 0.), (1, 1.), (2, 3.), (3, 2.)];
const ROUTES: [ConductionRegion; 4] = [
    ConductionRegion {
        triangle: 0,
        cell: CellId(7),
    },
    ConductionRegion {
        triangle: 1,
        cell: CellId(7),
    },
    ConductionRegion {
        triangle: 2,
        cell: CellId(7),
    },
    ConductionRegion {
        triangle: 3,
        cell: CellId(7),
    },
];
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn card(alpha: f64, electrical_upper: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let rho = ElectricalPropertyBasis::Resistivity;
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, values, upper) in [
        (heat.clone(), [2., 2.], 500.),
        (density.clone(), [1., 1.], 500.),
        (
            PropertyKey::new(PRINCIPAL_RESISTIVITY_PROPERTIES[0], rho.dimensions()),
            [1. - 10. * alpha, 1. + (electrical_upper - 300.) * alpha],
            electrical_upper,
        ),
        (
            PropertyKey::new(PRINCIPAL_RESISTIVITY_PROPERTIES[1], rho.dimensions()),
            [4., 4.],
            500.,
        ),
    ] {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(290., values[0]), (upper, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 290., upper),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic tensor thermal fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic conductor".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(290.),
            upper: &point(500.),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.,
        1e-8,
    )
    .unwrap();
    (card, curve)
}
fn initial(curve: &ResolvedHeatCapacityCurve) -> DepositionSimulation<'_> {
    let cell = DepositedCell {
        id: CellId(7),
        segment_id: 0,
        body_id: 0,
        material_id: MaterialId(3),
        mass_kg: 1.,
        specific_h_j_kg: 20.,
        geometry: CellGeometry {
            start_m: [0., 0., 0.],
            end_m: [1., 0., 0.],
            volume_m3: 1.,
            surface_area_m2: 6.,
        },
        boundary: Boundary {
            ambient_k: 300.,
            htc_w_m2_k: 0.,
            reservoir: None,
            source_w: 0.,
        },
    };
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(3), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 8,
        max_dofs: 32,
        max_triplets: 128,
    }
}
fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 100,
            max_corrections: 2,
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn problem(card: &MaterialCard, angle: f64) -> PlanarConduction<ResolvedOrthotropicConductor2> {
    rotated_problem(card, angle, 0.)
}
fn rotated_problem(
    card: &MaterialCard,
    angle: f64,
    geometry_angle: f64,
) -> PlanarConduction<ResolvedOrthotropicConductor2> {
    let frame =
        PlanarElectricalFrame::new(ContentHash([1; 32]), ContentHash([2; 32]), angle).unwrap();
    let source = resolve_orthotropic_conductor2(
        card,
        &point(300.),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
        frame,
    )
    .unwrap();
    let (sin, cos) = geometry_angle.sin_cos();
    let vertices = [[0., 0.], [1., 0.], [1., 1.], [0., 1.], [0.5, 0.5]]
        .into_iter()
        .map(|[x, y]| [cos * x - sin * y, sin * x + cos * y])
        .collect();
    let mesh = TriComplex2::from_indexed_triangles(
        "anisotropic-heated-square",
        vertices,
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(1.).unwrap(),
    )
    .unwrap();
    PlanarConduction::new(mesh, vec![source; 4], limits(), &mut || false).unwrap()
}
fn binding<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(3),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }
}
fn network<'a>(cells: &'a [MaterialCell<'a>]) -> ConductionThermalNetwork<'a> {
    ConductionThermalNetwork {
        cells,
        links: &[],
        regions: &ROUTES,
        volume_tolerance_m3: 1e-12,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-6, "{a} != {b}");
}
fn endpoint(previous: f64, dt: f64, a: f64, b: f64, alpha: f64) -> f64 {
    // C(theta-theta_old)=dt*(a/(1+alpha*theta)+b), C=2.
    let k = 2. * previous + dt * b;
    let linear = 2. - alpha * k;
    2. * (k + dt * a) / (linear + (linear * linear + 8. * alpha * (k + dt * a)).sqrt())
}

#[test]
fn oriented_temperature_sources_close_the_same_spatial_nonlinear_heat_balance() {
    let alpha = 0.02;
    let (card, curve) = card(alpha, 500.);
    let thermal = initial(&curve);
    let context = point(300.);
    let cells = [binding(&card, &curve, &context)];
    let mut rises = Vec::new();
    for (angle, a, b) in [(0., 1., 1.), (core::f64::consts::FRAC_PI_4, 4.5, 0.125)] {
        let old = problem(&card, angle);
        let expected = endpoint(0., 0.5, a, b, alpha);
        let accepted = network(&cells)
            .step_candidate(
                &thermal,
                0.5,
                &old,
                &VOLTAGES,
                limits(),
                field_controls(),
                config(),
                &[0.],
                || false,
            )
            .unwrap();
        let temperature = accepted.thermal().temperature_k(CellId(7)).unwrap();
        near(temperature, 300. + expected);
        let report = accepted.report();
        near(report.electrical.heat_j, 2. * expected);
        near(report.electrical.supplied_work_j, 2. * expected);
        for power in &report.electrical.cell_power_w {
            near(*power, (a / (1. + alpha * expected) + b) / 4.);
        }
        for (old, new) in old.conductors().iter().zip(accepted.problem().conductors()) {
            near(
                new.principal_conductivity_s_m()[0],
                1. / (1. + alpha * expected),
            );
            near(new.principal_conductivity_s_m()[1], 0.25);
            assert_eq!(new.frame(), old.frame());
            assert_ne!(new.identity(), old.identity());
            for (a, b) in old
                .material()
                .properties()
                .iter()
                .zip(new.material().properties())
            {
                assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
            }
        }
        // The accepted field genuinely has transverse tensor transport; use
        // the endpoint matrix to verify J, independently of its heat account.
        let sigma = accepted.problem().conductors()[0].conductivity_s_m();
        for j in &report.electrical.current_density_a_m2 {
            near(j[0], -sigma[0][0] - 2. * sigma[0][1]);
            near(j[1], -sigma[1][0] - 2. * sigma[1][1]);
        }
        let next_expected = endpoint(expected, 0.5, a, b, alpha);
        let continued = network(&cells)
            .step_candidate(
                accepted.thermal(),
                1.,
                accepted.problem(),
                &VOLTAGES,
                limits(),
                field_controls(),
                config(),
                &[0.],
                || false,
            )
            .unwrap();
        near(
            continued.thermal().temperature_k(CellId(7)).unwrap(),
            300. + next_expected,
        );
        near(
            continued.report().electrical.heat_j,
            2. * (next_expected - expected),
        );
        rises.push(expected);
    }
    assert!(rises[1] > rises[0]);
    near(thermal.temperature_k(CellId(7)).unwrap(), 300.);
}

#[test]
fn tensor_source_domain_and_late_cancellation_publish_no_partial_heat_or_field() {
    for upper in [300.01, 500.] {
        let (card, curve) = card(0.02, upper);
        let thermal = initial(&curve);
        let original = problem(&card, 0.4);
        let context = point(300.);
        let cells = [binding(&card, &curve, &context)];
        let sources = original.conductors().to_vec();
        if upper < 301. {
            assert!(
                network(&cells)
                    .step_candidate(
                        &thermal,
                        0.5,
                        &original,
                        &VOLTAGES,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.],
                        || false
                    )
                    .is_err()
            );
        } else {
            let mut calls = 0;
            network(&cells)
                .step_candidate(
                    &thermal,
                    0.5,
                    &original,
                    &VOLTAGES,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.],
                    || {
                        calls += 1;
                        false
                    },
                )
                .unwrap();
            let mut count = 0;
            assert!(
                network(&cells)
                    .step_candidate(
                        &thermal,
                        0.5,
                        &original,
                        &VOLTAGES,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.],
                        || {
                            count += 1;
                            count == calls - 1
                        }
                    )
                    .is_err()
            );
        }
        assert_eq!(original.conductors(), sources);
        near(thermal.time_s(), 0.);
        near(thermal.temperature_k(CellId(7)).unwrap(), 300.);
    }
}

#[test]
fn rotated_principal_resistivity_controls_native_circuit_storage_and_local_heat() {
    let alpha = 0.02;
    let capacitance = 0.5;
    let dt = 0.2;
    let geometry_angle = core::f64::consts::FRAC_PI_4;
    let (sin, cos) = geometry_angle.sin_cos();
    let (card, curve) = card(alpha, 500.);
    let context = point(300.);
    let cells = [binding(&card, &curve, &context)];
    let mut outcomes = Vec::new();
    for parallel in [true, false] {
        let angle = geometry_angle
            + if parallel {
                0.
            } else {
                core::f64::consts::FRAC_PI_2
            };
        let field = rotated_problem(&card, angle, geometry_angle);
        let sources = field.conductors().to_vec();
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: 999. }),
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: capacitance,
                    },
                ),
            ],
            transformers: vec![],
        };
        let mut circuit = ConductionCircuit::new(
            field,
            graph,
            ConductionCircuitTerminal {
                branch: 0,
                positive: vec![1, 2],
                negative: vec![0, 3],
            },
            0.,
            vec![],
            vec![2.],
            &[],
            field_controls(),
            || false,
        )
        .unwrap();
        let mut thermal = initial(&curve);
        let initial_heat = thermal.energy_j();
        let mut q0 = 2.0;
        let mut t0 = 300.;
        for step in 1..=3 {
            // Independent discrete RC/heat equations for a unit square whose
            // principal axis follows its electrodes and insulated side walls:
            // q1/q0=(2 R(T1) C-dt)/(2 R(T1) C+dt),
            // 2 (T1-T0)=(q0^2-q1^2)/(2 C).
            // The mesh is rotated in the geometry plane, so the actual field
            // still uses off-diagonal tensor entries in both orientations.
            let mut lo = 0.;
            let mut hi = q0;
            for _ in 0..90 {
                let q = 0.5 * (lo + hi);
                let t = t0 + (q0 * q0 - q * q) / (4. * capacitance);
                let r = if parallel {
                    1. + alpha * (t - 300.)
                } else {
                    4.
                };
                let rhs = q0 * (2. * r * capacitance - dt) / (2. * r * capacitance + dt);
                if q < rhs {
                    lo = q;
                } else {
                    hi = q;
                }
            }
            let q1 = 0.5 * (lo + hi);
            let t1 = t0 + (q0 * q0 - q1 * q1) / (4. * capacitance);
            let r = if parallel {
                1. + alpha * (t1 - 300.)
            } else {
                4.
            };
            let current = (q0 - q1) / dt;
            let heat = (q0 * q0 - q1 * q1) / (2. * capacitance);
            let next = network(&cells)
                .step_circuit_candidate(
                    &thermal,
                    step as f64 * dt,
                    &circuit,
                    &[],
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.],
                    || false,
                )
                .unwrap();
            let report = next.report();
            let electrical = &report.electrical;
            near(next.circuit().charges_c()[0], q1);
            near(next.thermal().temperature_k(CellId(7)).unwrap(), t1);
            near(electrical.spatial_current_a, current);
            near(electrical.spatial_voltage_v, current * r);
            near(electrical.native.supplied, 0.);
            near(electrical.heat_j, heat);
            near(report.deposited_joule_j, heat);
            near(next.thermal().energy_j() - thermal.energy_j(), heat);
            near(
                next.thermal().energy_j() + electrical.stored_after_j,
                initial_heat + 4.,
            );
            assert!(report.total_energy_residual_j.abs() < config().total_energy_tolerance_j);
            for triangle in 0..4 {
                let e = electrical.spatial.electric_field_v_m[triangle];
                let j = electrical.spatial.current_density_a_m2[triangle];
                near(e[0], -current * r * cos);
                near(e[1], -current * r * sin);
                near(j[0], -current * cos);
                near(j[1], -current * sin);
                near(
                    electrical.spatial.cell_power_w[triangle],
                    current * current * r / 4.,
                );
                let source = &next.circuit().problem().conductors()[triangle];
                near(
                    source.principal_conductivity_s_m()[0],
                    1. / (1. + alpha * (t1 - 300.)),
                );
                assert!(source.conductivity_s_m()[0][1].abs() > 0.1);
                assert_eq!(source.frame(), sources[triangle].frame());
                for (before, after) in sources[triangle]
                    .material()
                    .properties()
                    .iter()
                    .zip(source.material().properties())
                {
                    assert_eq!(
                        before.answer().receipt.selected,
                        after.answer().receipt.selected
                    );
                }
            }
            let (next_thermal, next_circuit, _) = next.into_parts();
            thermal = next_thermal;
            circuit = next_circuit;
            q0 = q1;
            t0 = t1;
        }
        outcomes.push((q0, t0));
    }
    assert!(outcomes[0].0 < outcomes[1].0);
    assert!(outcomes[0].1 > outcomes[1].1);
}
