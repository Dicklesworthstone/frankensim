//! Independent reversible-dielectric force, circuit and motion regressions.
//! Constitutive sources are authored fixtures, not measured dielectric data.
#[path = "electrical_contacts/source.rs"]
#[allow(dead_code)]
mod source;

use fs_couple::material_circuit::motion::dielectric::{
    ReversibleDielectricGap, ReversibleDielectricGapMotion, ReversibleDielectricMotionInitial,
};
use fs_couple::material_circuit::motion::support::AxialSupport;
use fs_couple::material_circuit::motion::{FieldMotionInput, MaterialFieldMotion};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor};
use fs_material::field_storage::dielectric::reversible::{
    ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
    ReversibleDielectricMaterial,
};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor,
};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::capacitance::{CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};

fn card(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored dielectric actuator".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "authored reversible actuator fixture".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn field_material(value: PropertyValue, maximum: f64) -> ReversibleDielectricMaterial {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
                ELECTRIC_DISPLACEMENT_DIMS,
            ),
            value,
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "E",
                    QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
                    0.0,
                    maximum,
                )
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                    450.0,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let point = QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    ReversibleDielectricMaterial::resolve(
        &card(claims),
        &point,
        "E",
        maximum,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn curve(points: &[(f64, f64)]) -> ReversibleDielectricMaterial {
    field_material(
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots: points.to_vec(),
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        points.last().unwrap().0,
    )
}
fn nonlinear() -> ReversibleDielectricMaterial {
    curve(&[(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)])
}
// Closed forms derived directly from the authored D(E) segments.
fn density(d: f64) -> f64 {
    let d = d.abs();
    if d <= 2.0 {
        d * d / 4.0
    } else if d <= 3.0 {
        let x = d - 2.0;
        1.0 + x + x * x
    } else {
        let x = d - 3.0;
        3.0 + 3.0 * x + x * x / 8.0
    }
}
fn field(d: f64) -> f64 {
    let a = d.abs();
    d.signum()
        * if a <= 2.0 {
            a / 2.0
        } else if a <= 3.0 {
            1.0 + 2.0 * (a - 2.0)
        } else {
            3.0 + (a - 3.0) / 4.0
        }
}
fn body(mass: f64) -> ResolvedBody {
    let card = source::scalar_card("density", Density::DIMS, mass);
    resolve_body(
        &[RegionBinding {
            key: "guided electrode",
            material: &card,
            point: &QueryPoint::new(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        }],
        1,
        || false,
    )
    .unwrap()
}
fn resistor(branch: usize, resistance: f64) -> MaterialResistor {
    let card = source::scalar_card(
        "electrical_resistivity",
        ELECTRICAL_RESISTIVITY_DIMS,
        resistance,
    );
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
fn graph(resistive: bool) -> CircuitGraph {
    let mut branches = vec![(1, 0, Branch::Capacitor { farads: 123.0 })];
    if resistive {
        branches.push((1, 0, Branch::Resistor { ohms: f64::NAN }));
    }
    CircuitGraph {
        node_count: 2,
        branches,
        transformers: vec![],
    }
}
fn initial(q: f64, gap: f64, velocity: f64) -> ReversibleDielectricMotionInitial {
    ReversibleDielectricMotionInitial {
        fluxes_wb: vec![],
        charges_c: vec![q],
        gap_m: gap,
        velocity_m_s: velocity,
        electrical_inputs: vec![],
    }
}
fn isolated(
    material: ReversibleDielectricMaterial,
    q: f64,
    gap: f64,
    velocity: f64,
    mass: f64,
    support: Option<AxialSupport>,
) -> ReversibleDielectricGapMotion {
    ReversibleDielectricGapMotion::new(
        graph(false),
        vec![],
        vec![],
        vec![],
        vec![],
        ReversibleDielectricGap::new(material, 0, 0.5, [0.2, 1.8]).unwrap(),
        body(mass),
        support,
        initial(q, gap, velocity),
    )
    .unwrap()
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 2e-10,
        power_tolerance_w: 2e-10,
    }
}
fn input(dt: f64, force: f64) -> FieldMotionInput<'static> {
    FieldMotionInput {
        electrical_ports: &[],
        external_force_n: force,
        duration_s: dt,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn fixed_charge_force_is_the_energy_derivative_for_both_polarities_and_source_masses() {
    for q in [-1.3, 0.0, 1.3] {
        for mass in [2.0, 7.0] {
            let moving = isolated(nonlinear(), q, 1.0, 0.02, mass, None);
            let before = moving.clone();
            let force = -0.5 * density(q / 0.5);
            let energy = 0.5 * density(q / 0.5);
            near(moving.field_force_n().unwrap(), force, 2e-15);
            near(moving.voltage_v().unwrap(), field(q / 0.5), 2e-15);
            near(moving.dielectric_energy_j().unwrap(), energy, 2e-15);
            if q != 0.0 {
                assert!(
                    (energy - 0.5 * q * moving.voltage_v().unwrap()).abs() > 0.4,
                    "nonlinear storage must not use half qV"
                );
            }
            let h = 1e-5;
            let plus = isolated(nonlinear(), q, 1.0 + h, 0.0, mass, None)
                .dielectric_energy_j()
                .unwrap();
            let minus = isolated(nonlinear(), q, 1.0 - h, 0.0, mass, None)
                .dielectric_energy_j()
                .unwrap();
            near(-(plus - minus) / (2.0 * h), force, 2e-11);
            let (dt, external) = (0.1, 0.2);
            let acceleration = (force + external) / mass;
            let gap = 1.0 + 0.02 * dt + 0.5 * acceleration * dt * dt;
            let velocity = 0.02 + acceleration * dt;
            let proposal = moving
                .step_candidate(vec![], input(dt, external), controls(), || false)
                .unwrap();
            let next = proposal.model();
            let step = proposal.step();
            near(next.gap_m(), gap, 3e-12);
            near(next.velocity_m_s(), velocity, 3e-12);
            near(next.charge_c(), q, 1e-14);
            near(step.electrical_energy_change_j, -force * (gap - 1.0), 3e-12);
            near(
                step.kinetic_energy_change_j,
                0.5 * mass * (velocity * velocity - 0.02 * 0.02),
                3e-12,
            );
            near(step.field_work_on_body_j, force * (gap - 1.0), 3e-12);
            near(step.mechanical_supply_j, external * (gap - 1.0), 3e-12);
            near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
            assert_eq!(step.heat_j, 0.0);
            assert_eq!(next.dielectric().material(), moving.dielectric().material());
            assert_eq!(next.body(), moving.body());
            assert_eq!(moving, before);
        }
    }
}

// Independent continuous ODE: q'=-d E(q/A)/R, d'=v, m v'=F-A w(q/A).
fn rhs(x: [f64; 3]) -> [f64; 3] {
    [
        -x[1] * field(x[0] / 0.5) / 2.0,
        x[2],
        (0.3 - 0.5 * density(x[0] / 0.5)) / 2.0,
    ]
}
fn rk4() -> [f64; 3] {
    let mut x = [1.3, 1.0, 0.04];
    let dt = 0.8 / 16000.0;
    for _ in 0..16000 {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + dt * c[i]));
        x = std::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}
#[test]
fn nonlinear_rc_motion_crosses_source_knots_and_refines_to_the_independent_ode() {
    let exact = rk4();
    let run = |steps: u32| {
        let mut moving = ReversibleDielectricGapMotion::new(
            graph(true),
            vec![resistor(1, 2.0)],
            vec![],
            vec![],
            vec![],
            ReversibleDielectricGap::new(nonlinear(), 0, 0.5, [0.2, 1.8]).unwrap(),
            body(2.0),
            None,
            initial(1.3, 1.0, 0.04),
        )
        .unwrap();
        let energy0 = 0.5 * density(2.6) + 0.5 * 2.0 * 0.04 * 0.04;
        let mut heat = 0.0;
        for _ in 0..steps {
            let proposal = moving
                .step_candidate(
                    vec![resistor(1, 2.0)],
                    input(0.8 / f64::from(steps), 0.3),
                    controls(),
                    || false,
                )
                .unwrap();
            let (next, step) = proposal.into_parts();
            assert!(step.heat_j > 0.0);
            near(step.heat_j, step.resistors[0].heat_j, 1e-14);
            heat += step.heat_j;
            near(
                next.dielectric_energy_j().unwrap() + next.velocity_m_s().powi(2) + heat,
                energy0 + 0.3 * (next.gap_m() - 1.0),
                2e-9,
            );
            moving = next;
        }
        assert!(
            moving.charge_c() < 1.0,
            "the source D=2 knot must be crossed"
        );
        let answer = [moving.charge_c(), moving.gap_m(), moving.velocity_m_s()];
        answer
            .into_iter()
            .zip(exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
    };
    let coarse = run(20);
    let medium = run(40);
    let fine = run(80);
    assert!(
        medium < coarse / 2.8,
        "first step refinement: coarse {coarse:e}, medium {medium:e}"
    );
    assert!(
        fine < medium / 2.8,
        "second step refinement: medium {medium:e}, fine {fine:e}"
    );
    assert!(fine < 2e-5, "fine independent state error {fine:e}");
}

#[test]
fn linear_material_limit_matches_the_existing_reciprocal_motion_owner() {
    let scalar = source::scalar_card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        2.0,
    );
    let cap = ResolvedCapacitor::resolve(
        &scalar,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.5,
        1.0,
    )
    .unwrap();
    let old = MaterialCircuit::new_with_storage(
        graph(true),
        vec![resistor(1, 2.0)],
        vec![MaterialStorage {
            branch: 0,
            element: MaterialStorageElement::Capacitor(cap),
        }],
        vec![],
        vec![1.3],
        &[],
    )
    .unwrap();
    let mut old = MaterialFieldMotion::new(old, 0, body(2.0), 0.03, [0.2, 1.8]).unwrap();
    let mut new = ReversibleDielectricGapMotion::new(
        graph(true),
        vec![resistor(1, 2.0)],
        vec![],
        vec![],
        vec![],
        ReversibleDielectricGap::new(curve(&[(0.0, 0.0), (4.0, 8.0)]), 0, 0.5, [0.2, 1.8]).unwrap(),
        body(2.0),
        None,
        initial(1.3, 1.0, 0.03),
    )
    .unwrap();
    for _ in 0..8 {
        let a = old
            .step_candidate(vec![resistor(1, 2.0)], input(0.02, 0.1), controls(), || {
                false
            })
            .unwrap();
        let b = new
            .step_candidate(vec![resistor(1, 2.0)], input(0.02, 0.1), controls(), || {
                false
            })
            .unwrap();
        near(a.model().gap_m(), b.model().gap_m(), 2e-10);
        near(a.model().velocity_m_s(), b.model().velocity_m_s(), 2e-10);
        near(
            a.model().circuit().charges_c()[0],
            b.model().charge_c(),
            2e-10,
        );
        near(a.step().heat_j, b.step().heat_j, 2e-10);
        old = a.into_parts().0;
        new = b.into_parts().0;
    }
}

#[test]
fn unrelated_nonlinear_capacitor_mutual_inductance_and_contact_heat_survive_motion() {
    let graph = CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 123.0 }),
            (1, 0, Branch::Resistor { ohms: 1.0 }),
            (2, 0, Branch::Capacitor { farads: 99.0 }),
            (2, 0, Branch::Resistor { ohms: 1.0 }),
            (3, 0, Branch::Inductor { henries: 2.0 }),
            (4, 0, Branch::Inductor { henries: 3.0 }),
            (3, 0, Branch::Resistor { ohms: 1.0 }),
            (4, 0, Branch::Resistor { ohms: 1.0 }),
        ],
        transformers: vec![],
    };
    let extra = CapacitorChargeCurve::try_new(vec![
        (-2.0, -3.0),
        (-0.5, -0.5),
        (0.0, 0.0),
        (0.5, 0.5),
        (2.0, 3.0),
    ])
    .unwrap();
    let bindings = vec![CapacitorCurve {
        branch: 2,
        curve: extra.clone(),
    }];
    let mutual = vec![MutualInductance {
        branches: (4, 5),
        henries: 0.4,
    }];
    let bulk = || {
        vec![
            resistor(1, 2.0),
            resistor(3, 1.0),
            resistor(6, 1.2),
            resistor(7, 0.7),
        ]
    };
    let a = source::scalar_card("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 2.0);
    let b = source::scalar_card("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 1.0);
    let interface = source::card(&a, &b, 0.4, 0.0);
    let contacts = vec![source::contact(11, 1, &interface, 1.0)];
    let mut state = initial(1.3, 1.0, 0.0);
    state.charges_c.push(-0.7);
    state.fluxes_wb = vec![0.4, -0.3];
    let mut moving = ReversibleDielectricGapMotion::new(
        graph,
        bulk(),
        contacts.clone(),
        mutual.clone(),
        bindings.clone(),
        ReversibleDielectricGap::new(nonlinear(), 0, 0.5, [0.2, 1.8]).unwrap(),
        body(2.0),
        None,
        state,
    )
    .unwrap();
    let magnetic = |x: &[f64]| {
        0.5 * (3.0 * x[0] * x[0] - 0.8 * x[0] * x[1] + 2.0 * x[1] * x[1]) / (6.0 - 0.16)
    };
    for _ in 0..12 {
        let proposal = moving
            .step_with_contacts_candidate(
                bulk(),
                contacts.clone(),
                input(0.03, 0.2),
                controls(),
                || false,
            )
            .unwrap();
        let (next, step) = proposal.into_parts();
        let routing = next.routing_circuit();
        let physical = next.dielectric_energy_j().unwrap()
            + extra.stored_energy_j(routing.charges_c()[1]).unwrap()
            + magnetic(routing.fluxes_wb());
        near(step.electrical_energy_j[1], physical, 2e-12);
        let old = moving.routing_circuit();
        let change = next.dielectric_energy_j().unwrap() - moving.dielectric_energy_j().unwrap()
            + extra.stored_energy_j(routing.charges_c()[1]).unwrap()
            - extra.stored_energy_j(old.charges_c()[1]).unwrap()
            + magnetic(routing.fluxes_wb())
            - magnetic(old.fluxes_wb());
        near(step.electrical_energy_change_j, change, 2e-12);
        assert_eq!(routing.capacitor_curves(), bindings.as_slice());
        assert_eq!(routing.mutual_inductance(), mutual.as_slice());
        assert_eq!(routing.contacts(), contacts.as_slice());
        let current = step
            .resistors
            .iter()
            .find(|heat| heat.material.branch == 1)
            .unwrap()
            .current_a;
        near(
            step.contacts[0].heat_j,
            0.03 * 0.4 * current * current,
            2e-12,
        );
        near(step.contacts[0].current_a, current, 2e-12);
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert_eq!(
            routing.graph().branches[0].2,
            Branch::Capacitor { farads: 123.0 }
        );
        moving = next;
    }
    assert!(
        moving.routing_circuit().charges_c()[1] > -0.5,
        "the unrelated capacitor must cross its own knot"
    );
    assert!((moving.routing_circuit().fluxes_wb()[0] - 0.4).abs() > 1e-3);
}

fn support() -> AxialSupport {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 3.0),
        ("young_modulus", Pressure::DIMS, 100.0),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            0.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 250.0, 450.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let card = card(claims);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &[
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
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let expansion = integrate_isotropic_thermal_expansion(
        &card,
        &point,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: (0.1 / core::f64::consts::PI).sqrt(),
            material,
        }],
    )
    .unwrap();
    AxialSupport::try_new(
        ConcentricThermoelasticSection::try_new(section, vec![expansion]).unwrap(),
        2.0,
        1.0,
        0.2,
    )
    .unwrap()
}
#[test]
fn original_material_support_adds_reciprocal_restoring_force_energy_and_modal_inertia() {
    let support = support();
    let moving = isolated(nonlinear(), 1.3, 1.02, 0.02, 2.0, Some(support.clone()));
    let mass = 2.0 + 3.0 * 0.1 * 2.0 / 3.0;
    let stiffness = 100.0 * 0.1 / 2.0;
    near(moving.generalized_mass_kg(), mass, 1e-14);
    let dt = 0.05;
    let external = 0.2;
    let a = dt * dt * stiffness / (4.0 * mass);
    let force = external - 0.5 * density(2.6);
    let v = ((1.0 - a) * 0.02 + dt * (force - stiffness * 0.02) / mass) / (1.0 + a);
    let gap = 1.02 + 0.5 * dt * (0.02 + v);
    let proposal = moving
        .step_candidate(vec![], input(dt, external), controls(), || false)
        .unwrap();
    near(proposal.model().gap_m(), gap, 2e-11);
    near(proposal.model().velocity_m_s(), v, 2e-11);
    near(
        proposal.step().elastic_energy_j[1],
        0.5 * stiffness * (gap - 1.0).powi(2),
        2e-12,
    );
    near(
        proposal.step().elastic_energy_change_j,
        0.5 * stiffness * ((gap - 1.0).powi(2) - 0.02 * 0.02),
        2e-12,
    );
    assert_eq!(proposal.model().support(), Some(&support));
    near(
        proposal.step().energy_residual_j,
        0.0,
        controls().energy_tolerance_j,
    );
}

#[test]
fn small_si_dielectric_motion_keeps_work_below_an_unrelated_large_stored_bias() {
    let epsilon = 8.8541878128e-12;
    let material = curve(&[(0.0, 0.0), (1e6, epsilon * 1e6), (3e6, epsilon * 2e6)]);
    let (area, gap, dt) = (0.01, 1e-4, 1e-3);
    let q = area * epsilon * 1.5e6;
    let w = 0.5 * epsilon * 1e12 + 1.5e6 * (0.5 * epsilon * 1e6);
    let force = -area * w;
    let moving = ReversibleDielectricGapMotion::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (2, 0, Branch::Capacitor { farads: 1.0 }),
            ],
            transformers: vec![],
        },
        vec![],
        vec![],
        vec![],
        vec![],
        ReversibleDielectricGap::new(material, 0, area, [5e-5, 2e-4]).unwrap(),
        body(1.0),
        None,
        ReversibleDielectricMotionInitial {
            charges_c: vec![q, 1e8],
            fluxes_wb: vec![],
            gap_m: gap,
            velocity_m_s: 0.0,
            electrical_inputs: vec![],
        },
    )
    .unwrap();
    let strict = MaterialCircuitSolve {
        energy_tolerance_j: 1e-16,
        power_tolerance_w: 1e-13,
    };
    let proposal = moving
        .step_candidate(vec![], input(dt, 0.0), strict, || false)
        .unwrap();
    let displacement = 0.5 * force * dt * dt;
    near(proposal.model().gap_m(), gap + displacement, 2e-17);
    near(proposal.model().velocity_m_s(), force * dt, 2e-14);
    near(
        proposal.step().electrical_energy_change_j,
        -force * displacement,
        2e-17,
    );
    near(
        proposal.step().kinetic_energy_change_j,
        0.5 * (force * dt).powi(2),
        2e-17,
    );
    near(
        proposal.step().energy_residual_j,
        0.0,
        strict.energy_tolerance_j,
    );
    assert_eq!(
        proposal.model().routing_circuit().charges_c()[1].to_bits(),
        1e8_f64.to_bits()
    );
    // Absolute stored-energy diagnostics cannot represent this tiny transfer;
    // the published balance must use the independently evaluated increments.
    assert_eq!(
        proposal.step().electrical_energy_j[0],
        proposal.step().electrical_energy_j[1]
    );
    assert_ne!(proposal.step().electrical_energy_change_j, 0.0);
}

#[test]
fn unsupported_source_axes_domains_constraints_and_cancelled_proposals_are_atomic() {
    let moving = isolated(nonlinear(), 1.3, 1.0, 0.02, 2.0, None);
    let accepted = moving.clone();
    assert!(matches!(
        moving.step_candidate(vec![], input(0.1, 0.0), controls(), || true),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert!(
        moving
            .step_candidate(vec![], input(1.0, 1e6), controls(), || false)
            .is_err()
    );
    assert!(
        moving
            .step_candidate(vec![], input(0.0, 0.0), controls(), || false)
            .is_err()
    );
    assert_eq!(moving, accepted);
    let build = |state: ReversibleDielectricMotionInitial, graph: CircuitGraph| {
        ReversibleDielectricGapMotion::new(
            graph,
            vec![],
            vec![],
            vec![],
            vec![],
            ReversibleDielectricGap::new(nonlinear(), 0, 0.5, [0.2, 1.8]).unwrap(),
            body(2.0),
            None,
            state,
        )
    };
    assert!(build(initial(3.5001, 1.0, 0.0), graph(false)).is_err());
    assert!(build(initial(1.3, 0.2, 0.0), graph(false)).is_err());
    assert!(
        ReversibleDielectricGapMotion::new(
            graph(false),
            vec![],
            vec![],
            vec![],
            vec![],
            ReversibleDielectricGap::new(nonlinear(), 0, 4.0, [0.2, 1.8]).unwrap(),
            body(2.0),
            None,
            initial(f64::from_bits(1), 1.0, 0.0)
        )
        .is_err(),
        "nonzero charge must not underflow into an uncharged material response"
    );
    let mut clamped = graph(false);
    clamped
        .branches
        .push((1, 0, Branch::VoltageSource { port: 0 }));
    let mut wrong = initial(1.3, 1.0, 0.0);
    wrong.electrical_inputs = vec![1.0];
    assert!(
        build(wrong, clamped.clone()).is_err(),
        "ideal clamp must not rewrite charge"
    );
    let mut compatible = initial(1.3, 1.0, 0.0);
    compatible.electrical_inputs = vec![2.2];
    let compatible = build(compatible, clamped).unwrap();
    assert!(
        compatible
            .step_candidate(
                vec![],
                FieldMotionInput {
                    electrical_ports: &[2.2],
                    external_force_n: 0.0,
                    duration_s: 0.1
                },
                controls(),
                || false
            )
            .is_err(),
        "a finite step must not conceal a nonlinear endpoint clamp inconsistency"
    );
    let dependent = field_material(
        PropertyValue::Surface2D {
            axes: ["E".into(), "T".into()],
            axis_dims: [ELECTRIC_FIELD_DIMS, Temperature::DIMS],
            samples: vec![
                ([0.0, 250.0], 0.0),
                ([1.0, 250.0], 2.0),
                ([4.0, 250.0], 5.0),
                ([0.0, 450.0], 0.0),
                ([1.0, 450.0], 3.0),
                ([4.0, 450.0], 8.0),
            ],
            triangles: vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        4.0,
    );
    assert!(
        ReversibleDielectricGap::new(dependent, 0, 0.5, [0.2, 1.8]).is_err(),
        "frozen T-dependent field potential cannot supply absent reversible heat"
    );
}
