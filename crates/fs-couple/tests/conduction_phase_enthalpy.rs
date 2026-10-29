//! Confined phase change uses one sourced field, circuit and total enthalpy.
//! All numbers below are synthetic constitutive data, not empirical metals.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "conduction_phase_enthalpy/dielectric.rs"]
mod dielectric;
#[path = "conduction_phase_enthalpy/response.rs"]
mod response;

use common::{config, near};
use fs_blake3::ContentHash;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    conduction::circuit::{ConductionCircuit, ConductionCircuitTerminal},
    conduction::{ConductionDrive, ConductionRegion, ConductionThermalNetwork},
    material::MaterialCell,
};
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::planar::{
        PRINCIPAL_RESISTIVITY_PROPERTIES, PlanarConductorMaterial, PlanarElectricalFrame,
        ResolvedOrthotropicConductor2, ResolvedPlanarConductor,
    },
    conductor::{ElectricalPropertyBasis, resolve_conductor_state_with_basis},
    phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve},
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use fs_thermal::enthalpy::*;

const FRACTION: &str = "liquid_mass_fraction";
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
fn point(t: f64, f: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity(FRACTION, QuantitySpec::dimensional(Dims([0; 6])), f)
        .unwrap()
}
fn card(slope: f64, upper_fraction: f64) -> MaterialCard {
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut claims = ClaimSet::new();
    for (name, base, slope) in [
        (basis.property_name(), 1., slope),
        (PRINCIPAL_RESISTIVITY_PROPERTIES[0], 1., slope),
        (PRINCIPAL_RESISTIVITY_PROPERTIES[1], 4., 2.),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, basis.dimensions()),
                value: PropertyValue::Curve {
                    abscissa: FRACTION.into(),
                    abscissa_dims: Dims([0; 6]),
                    knots: vec![(0., base), (upper_fraction, base + slope * upper_fraction)],
                    dims: basis.dimensions(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.,
                        350.,
                    )
                    .with_quantity(
                        FRACTION,
                        QuantitySpec::dimensional(Dims([0; 6])),
                        0.,
                        upper_fraction,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic equal-density phase conductor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid-liquid chart".into(),
            process: "authored test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(card: &MaterialCard, liquid_density: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0., 300., 0., 1.),
            (20., 310., 0., 1.),
            (120., 310., 1., liquid_density),
            (200., 350., 1., liquid_density),
        ]
        .into_iter()
        .map(|(h, t, f, density)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: density,
        })
        .collect(),
    )
    .unwrap()
}
fn initial(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.,
            event: Event::Birth(DepositedCell {
                id: CellId(7),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(3),
                mass_kg: 1.,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [0.; 3],
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
            }),
        }],
        &[(MaterialId(3), curve as &dyn ThermalLaw)],
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
fn controls() -> ConductionControls {
    ConductionControls {
        linear: CheckedCgConfig {
            tolerance: 1e-11,
            max_iterations: 100,
            max_corrections: 2,
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn mesh() -> TriComplex2 {
    let (sin, cos) = core::f64::consts::FRAC_PI_4.sin_cos();
    TriComplex2::from_indexed_triangles(
        "confined rotated phase conductor",
        [[0., 0.], [1., 0.], [1., 1.], [0., 1.], [0.5, 0.5]]
            .into_iter()
            .map(|[x, y]| [cos * x - sin * y, sin * x + cos * y])
            .collect(),
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(1.).unwrap(),
    )
    .unwrap()
}
fn source(card: &MaterialCard, context: &QueryPoint, tensor: bool) -> ResolvedPlanarConductor {
    let basis = ElectricalPropertyBasis::Resistivity;
    if !tensor {
        return resolve_conductor_state_with_basis(
            card,
            context,
            MaterialPropertySelection::SingleClaimOnly,
            basis,
        )
        .unwrap()
        .with_liquid_mass_fraction_axis(FRACTION)
        .unwrap()
        .into();
    }
    let requirements = PRINCIPAL_RESISTIVITY_PROPERTIES.map(|name| {
        ScalarPropertyRequirement::try_new(
            name,
            basis.dimensions(),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()
    });
    let material = resolve_material_state_point(
        card,
        context,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let frame = PlanarElectricalFrame::new(
        ContentHash([1; 32]),
        ContentHash([2; 32]),
        core::f64::consts::FRAC_PI_4,
    )
    .unwrap();
    ResolvedOrthotropicConductor2::try_from_phase_material(material, basis, frame, FRACTION)
        .unwrap()
        .into()
}
fn problem(
    card: &MaterialCard,
    context: &QueryPoint,
    tensor: bool,
) -> PlanarConduction<ResolvedPlanarConductor> {
    PlanarConduction::new(
        mesh(),
        vec![source(card, context, tensor); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn binding<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(3),
        curve,
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: Some(FRACTION),
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
fn current() -> ConductionDrive<'static> {
    ConductionDrive::Current {
        positive: &[1, 2],
        negative: &[0, 3],
        current_a: 2.,
    }
}
fn check_sources(
    old: &PlanarConduction<ResolvedPlanarConductor>,
    next: &PlanarConduction<ResolvedPlanarConductor>,
    t: f64,
    f: f64,
    slope: f64,
) {
    for (a, b) in old.conductors().iter().zip(next.conductors()) {
        assert_eq!(b.liquid_mass_fraction_axis(), Some(FRACTION));
        near(b.material().query().axes()["T"], t);
        near(b.material().query().axes()[FRACTION], f);
        assert_eq!(a.frame(), b.frame());
        for (a, b) in a.source_properties().zip(b.source_properties()) {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
        match b {
            ResolvedPlanarConductor::Isotropic(s) => near(s.resistivity_ohm_m(), 1. + slope * f),
            ResolvedPlanarConductor::Orthotropic(s) => {
                near(s.principal_conductivity_s_m()[0], 1. / (1. + slope * f));
                near(s.principal_conductivity_s_m()[1], 1. / (4. + 2. * f));
                assert!(s.conductivity_s_m()[0][1].abs() > 0.01);
            }
        }
    }
}

#[test]
fn current_and_voltage_follow_liquid_fraction_on_an_isothermal_plateau() {
    for tensor in [false, true] {
        for slope in [1., 3.] {
            let card = card(slope, 1.);
            let curve = curve(&card, 1.);
            let thermal = initial(&curve, 10.);
            let context = point(305., 0.);
            let old = problem(&card, &context, tensor);
            let cells = [binding(&card, &curve, &context)];
            let layout = network(&cells);
            // Independent BE balance: h1=10+5*2²*(1+slope*f),
            // h1=20+100*f on the T=310 K latent plateau (mass=1 kg).
            let f = 10. / (100. - 20. * slope);
            let h = 20. + 100. * f;
            let resistance = 1. + slope * f;
            let next = layout
                .step_driven_candidate(
                    &thermal,
                    5.,
                    &old,
                    current(),
                    limits(),
                    controls(),
                    config(),
                    &[20.],
                    || false,
                )
                .unwrap();
            near(next.thermal().energy_j(), h);
            near(next.thermal().temperature_k(CellId(7)).unwrap(), 310.);
            near(
                next.report().source.cells[0].phase.liquid_mass_fraction(),
                f,
            );
            near(next.report().electrical.heat_j, h - 10.);
            near(next.report().electrical.supplied_work_j, h - 10.);
            near(next.report().deposited_joule_j, h - 10.);
            near(next.report().electrical.voltage_v[1], 2. * resistance);
            let direction = core::f64::consts::FRAC_1_SQRT_2;
            for i in 0..4 {
                near(
                    next.report().electrical.current_density_a_m2[i][0],
                    -2. * direction,
                );
                near(
                    next.report().electrical.current_density_a_m2[i][1],
                    -2. * direction,
                );
                near(next.report().electrical.cell_power_w[i], resistance);
            }
            check_sources(&old, next.problem(), 310., f, slope);
            assert!(
                next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
            );
            // The independently chosen voltage V=I*R(h1) has the same endpoint;
            // an initial-temperature or temperature-only refresh cannot match it.
            let voltage = [(0, 0.), (3, 0.), (1, 2. * resistance), (2, 2. * resistance)];
            let voltage_next = layout
                .step_candidate(
                    &thermal,
                    5.,
                    &old,
                    &voltage,
                    limits(),
                    controls(),
                    config(),
                    &[20.],
                    || false,
                )
                .unwrap();
            near(voltage_next.thermal().energy_j(), h);
            check_sources(&old, voltage_next.problem(), 310., f, slope);
            if slope == 3. {
                let liquid = layout
                    .step_driven_candidate(
                        next.thermal(),
                        10.,
                        next.problem(),
                        current(),
                        limits(),
                        controls(),
                        config(),
                        &[80.],
                        || false,
                    )
                    .unwrap();
                // All 100 latent joules are already in h, followed by Cp=2.
                near(liquid.thermal().energy_j(), 125.);
                near(liquid.thermal().temperature_k(CellId(7)).unwrap(), 312.5);
                near(liquid.report().electrical.heat_j, 80.);
                check_sources(next.problem(), liquid.problem(), 312.5, 1., slope);
            }
            near(thermal.energy_j(), 10.);
            assert_eq!(thermal.time_s(), 0.);
        }
    }
}

fn capacitor(
    problem: PlanarConduction<ResolvedPlanarConductor>,
) -> ConductionCircuit<ResolvedPlanarConductor> {
    ConductionCircuit::new(
        problem,
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: 999. }),
                (1, 0, Branch::Capacitor { farads: 1. }),
            ],
            transformers: vec![],
        },
        ConductionCircuitTerminal {
            branch: 0,
            positive: vec![1, 2],
            negative: vec![0, 3],
        },
        0.,
        vec![],
        vec![12.],
        &[],
        controls(),
        || false,
    )
    .unwrap()
}
fn capacitor_reference(q0: f64, h0: f64, slope: f64) -> (f64, f64) {
    // C=1 F, dt=1 s. Independent scalar root combines exact storage loss
    // with the native midpoint RC relation at the endpoint phase resistance.
    let mut low = 0.;
    let mut high = q0;
    for _ in 0..90 {
        let q = 0.5 * (low + high);
        let h = h0 + 0.5 * (q0 * q0 - q * q);
        let r = 1. + slope * (h - 20.) / 100.;
        let rhs = q0 * (2. * r - 1.) / (2. * r + 1.);
        if q < rhs {
            low = q;
        } else {
            high = q;
        }
    }
    let q = 0.5 * (low + high);
    (q, h0 + 0.5 * (q0 * q0 - q * q))
}

#[test]
fn capacitor_storage_melts_the_same_spatial_material_and_changes_charge_decay() {
    let mut final_charges = Vec::new();
    for tensor in [false, true] {
        for slope in [0., 3.] {
            let card = card(slope, 1.);
            let curve = curve(&card, 1.);
            let thermal = initial(&curve, 20.);
            let context = point(310., 0.);
            let cells = [binding(&card, &curve, &context)];
            let layout = network(&cells);
            let circuit = capacitor(problem(&card, &context, tensor));
            let next = layout
                .step_circuit_candidate(
                    &thermal,
                    1.,
                    &circuit,
                    &[],
                    &[],
                    limits(),
                    controls(),
                    config(),
                    &[30.],
                    || false,
                )
                .unwrap();
            let (q1, h1) = capacitor_reference(12., 20., slope);
            near(next.circuit().charges_c()[0], q1);
            near(next.thermal().energy_j(), h1);
            near(next.thermal().temperature_k(CellId(7)).unwrap(), 310.);
            near(next.report().electrical.spatial_current_a, 12. - q1);
            near(next.report().electrical.native.supplied, 0.);
            near(next.report().deposited_joule_j, 72. - q1 * q1 / 2.);
            near(
                next.thermal().energy_j() + next.report().electrical.stored_after_j,
                92.,
            );
            check_sources(
                circuit.problem(),
                next.circuit().problem(),
                310.,
                (h1 - 20.) / 100.,
                slope,
            );
            let second = layout
                .step_circuit_candidate(
                    next.thermal(),
                    2.,
                    next.circuit(),
                    &[],
                    &[],
                    limits(),
                    controls(),
                    config(),
                    &[10.],
                    || false,
                )
                .unwrap();
            let (q2, h2) = capacitor_reference(q1, h1, slope);
            near(second.circuit().charges_c()[0], q2);
            near(second.thermal().energy_j(), h2);
            near(
                second.thermal().energy_j() + second.report().electrical.stored_after_j,
                92.,
            );
            check_sources(
                next.circuit().problem(),
                second.circuit().problem(),
                310.,
                (h2 - 20.) / 100.,
                slope,
            );
            assert_eq!(second.circuit().time_s(), second.thermal().time_s());
            assert_eq!(circuit.charges_c(), &[12.]);
            near(thermal.energy_j(), 20.);
            final_charges.push(q2);
        }
    }
    // Only source data changed; the phase-dependent higher resistance retains
    // more capacitor charge. Scalar and rotated tensor use the same owner.
    assert!(final_charges[1] > final_charges[0]);
    near(final_charges[0], final_charges[2]);
    near(final_charges[1], final_charges[3]);
}

#[test]
fn phase_admission_source_refusal_and_late_cancellation_preserve_all_inputs() {
    let full_card = card(3., 1.);
    let caloric = curve(&full_card, 1.);
    let thermal = initial(&caloric, 10.);
    let context = point(305., 0.);
    let raw = resolve_conductor_state_with_basis(
        &full_card,
        &context,
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    assert!(PlanarConduction::new(mesh(), vec![raw.clone(); 4], limits(), &mut || false).is_err());
    assert!(raw.clone().with_liquid_mass_fraction_axis("T").is_err());
    assert!(raw.with_liquid_mass_fraction_axis("missing").is_err());
    let old = problem(&full_card, &context, true);
    let original_sources = old.conductors().to_vec();
    let ResolvedPlanarConductor::Orthotropic(tensor) = &original_sources[0] else {
        panic!("tensor fixture required");
    };
    assert!(
        ResolvedOrthotropicConductor2::try_from_material(
            tensor.material().clone(),
            tensor.source_basis(),
            tensor.frame()
        )
        .is_err()
    );
    let cells = [binding(&full_card, &caloric, &context)];
    let layout = network(&cells);
    let circuit = capacitor(problem(&full_card, &context, true));
    let mut polls = 0;
    let good = layout
        .step_circuit_candidate(
            &thermal,
            1.,
            &circuit,
            &[],
            &[],
            limits(),
            controls(),
            config(),
            &[30.],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(good.thermal().energy_j() > 20.);
    let mut calls = 0;
    assert!(matches!(
        layout.step_circuit_candidate(
            &thermal,
            1.,
            &circuit,
            &[],
            &[],
            limits(),
            controls(),
            config(),
            &[30.],
            || {
                calls += 1;
                calls >= polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    for axis in [None, Some("another_phase")] {
        let mut cell = binding(&full_card, &caloric, &context);
        cell.point.liquid_fraction_axis = axis;
        assert!(
            network(&[cell])
                .step_driven_candidate(
                    &thermal,
                    5.,
                    &old,
                    current(),
                    limits(),
                    controls(),
                    config(),
                    &[20.],
                    || false
                )
                .is_err()
        );
    }
    let unequal_density = curve(&full_card, 2.);
    let unequal_thermal = initial(&unequal_density, 10.);
    let unequal_cell = binding(&full_card, &unequal_density, &context);
    assert!(
        network(&[unequal_cell])
            .step_driven_candidate(
                &unequal_thermal,
                5.,
                &old,
                current(),
                limits(),
                controls(),
                config(),
                &[20.],
                || false
            )
            .is_err()
    );
    let restricted = card(3., 0.2);
    let restricted_curve = curve(&restricted, 1.);
    let restricted_thermal = initial(&restricted_curve, 10.);
    let restricted_problem = problem(&restricted, &context, true);
    let restricted_cell = binding(&restricted, &restricted_curve, &context);
    assert!(
        network(&[restricted_cell])
            .step_driven_candidate(
                &restricted_thermal,
                5.,
                &restricted_problem,
                current(),
                limits(),
                controls(),
                config(),
                &[20.],
                || false
            )
            .is_err()
    );
    assert_eq!(old.conductors(), original_sources.as_slice());
    assert_eq!(circuit.charges_c(), &[12.]);
    assert_eq!(circuit.time_s(), 0.);
    assert_eq!(thermal.time_s(), 0.);
    near(thermal.energy_j(), 10.);
    near(restricted_thermal.energy_j(), 10.);
}
