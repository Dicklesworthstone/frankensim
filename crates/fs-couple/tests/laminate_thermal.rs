//! Native heat -> original material sources -> compatible MITC deformation.
//! Authored numerical cards and independent analytical oracles, not measurements.

use fs_couple::laminate_thermal::{
    LaminateShellDiscretization, LaminateShellError, LaminateShellLimits, LaminateShellSnapshot,
    LayeredThermalShell,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY as K_PROPERTY,
};
use fs_material::thermal::conduction::LayerConductionConfig;
use fs_material::thermal::transport::stack::LayerStackSourceQuery;
use fs_material::thermal::transport::stack::transient::{
    BoundLayeredThermalBody, LayeredThermalIds, PreparedLayeredThermalBody,
    TransientLayerSourceQuery,
    implicit::{LayeredThermalSolveControls, LayeredThermalWorkspace},
};
use fs_material::thermal::transport::storage::{HeatStorageAmount, HeatStorageComponentQuery};
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{DepositionBudget, DepositionSimulation, StepConfig};
use std::sync::Arc;

const AREA: f64 = 0.6;
const H: f64 = 0.02;
const RHO: f64 = 1000.0;
const CP: f64 = 100.0;
const K: f64 = 10.0;
const E: f64 = 2e8;

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn at(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn young(t: f64) -> f64 {
    E * (1.0 - 0.002 * (t - 300.0))
}
fn free(t: f64) -> f64 {
    let dt = t - 300.0;
    1e-5 * dt + 0.5e-7 * dt * dt
}
fn card(name: &str) -> Arc<MaterialCard> {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            vec![(280.0, 1040.0), (300.0, RHO), (360.0, 880.0)],
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            vec![(280.0, young(280.0)), (300.0, E), (360.0, young(360.0))],
        ),
        ("poisson_ratio", Dims::NONE, vec![(300.0, 0.0)]),
        (
            CTE,
            INVERSE_TEMPERATURE_DIMS,
            vec![(280.0, 8e-6), (320.0, 12e-6), (360.0, 16e-6)],
        ),
        (K_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, vec![(300.0, K)]),
        (CP_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, vec![(300.0, CP)]),
    ] {
        let constant = values.len() == 1;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if constant {
                    PropertyValue::Scalar {
                        value: values[0].1,
                        dims,
                    }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: values,
                        dims,
                    }
                },
                interpolation: if constant {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored coupled laminate oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "solid".into(),
                process: "numerical thermal shell".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
fn mechanical(cards: &[Arc<MaterialCard>]) -> Vec<LaminaSourceInput> {
    cards
        .iter()
        .map(|card| LaminaSourceInput {
            card: Arc::clone(card),
            reference: at(300.0),
            thickness_m: H / cards.len() as f64,
            orientation_rad: 0.0,
            elasticity: LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson {
                young_modulus: ScalarPropertyRequirement::try_new(
                    "young_modulus",
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
                poisson_ratio: ScalarPropertyRequirement::try_new(
                    "poisson_ratio",
                    Dims::NONE,
                    ScalarAdmissibility::Finite,
                )
                .unwrap(),
            }),
            expansion: LaminaExpansionSource::Isotropic,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect()
}
fn prepare_body(
    cards: &[Arc<MaterialCard>],
    counts: &[usize],
    mass_scale: f64,
) -> PreparedLayeredThermalBody {
    let lower = at(280.0);
    let upper = at(360.0);
    let conductivity = PropertyKey::new(K_PROPERTY, THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new(CP_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let h = H / cards.len() as f64;
    let mass = mass_scale * AREA * h * RHO;
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let queries: Vec<_> = cards
        .iter()
        .zip(counts)
        .enumerate()
        .map(|(i, (card, &count))| TransientLayerSourceQuery {
            conduction: LayerStackSourceQuery {
                id: i as u64,
                card,
                conductivity: &conductivity,
                thickness_m: h,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            capacity: HeatStorageComponentQuery {
                id: i as u64,
                card,
                capacity: &capacity,
                amount: HeatStorageAmount::MassKg(mass),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            mass_kg: mass,
            cell_count: count,
        })
        .collect();
    PreparedLayeredThermalBody::try_new(&queries, &vec![0.0; cards.len() - 1], AREA, 300.0).unwrap()
}
fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 10,
        first_material: 30,
        first_link: 50,
        body_id: 70,
    }
}
fn conduction() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    }
}
fn initial<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperatures: &[f64],
) -> DepositionSimulation<'a> {
    let events = bound
        .birth_events(0.0, temperatures, conduction(), || Ok(()))
        .unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: temperatures.len(),
            max_links: temperatures.len() - 1,
            max_events: 2 * temperatures.len() - 1,
            max_substeps: 64,
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
fn discretization() -> LaminateShellDiscretization {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.6, 0.0],
            [1.0, 0.6, 0.0],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 4], 0.9).unwrap();
    LaminateShellDiscretization {
        mesh,
        directors,
        constraints: vec![ShellConstraint {
            node: 0,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        }],
        section_angles_rad: vec![0.0; 2],
        mass: ShellMass::Consistent,
        budget: MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        drilling_alpha: 1e-3,
    }
}
fn limits() -> LaminateShellLimits {
    LaminateShellLimits {
        maximum_translation_m: 0.2,
        maximum_rotation_rad: 0.5,
        maximum_reduced_elastic_strain: 0.01,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}

// Independent scalar CLT: nu=0 decouples the two identical normal directions.
// Integrate each original cell's constant E and free strain analytically.
fn clt(thicknesses: &[f64], temperatures: &[f64]) -> (f64, f64, f64) {
    let (mut a, mut b, mut d, mut n, mut m) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut z = -0.5 * thicknesses.iter().sum::<f64>();
    for (&h, &t) in thicknesses.iter().zip(temperatures) {
        let top = z + h;
        let e = young(t);
        let moment = (top * top - z * z) / 2.0;
        a += e * h;
        b += e * moment;
        d += e * (top.powi(3) - z.powi(3)) / 3.0;
        n += e * free(t) * h;
        m += e * free(t) * moment;
        z = top;
    }
    let determinant = a * d - b * b;
    let eps = (d * n - b * m) / determinant;
    let kappa = (a * m - b * n) / determinant;
    let mut energy = 0.0;
    z = -0.5 * thicknesses.iter().sum::<f64>();
    for (&h, &t) in thicknesses.iter().zip(temperatures) {
        let top = z + h;
        let r = eps - free(t);
        // Two normal directions cancel the 1/2 in each quadratic energy.
        energy += young(t)
            * (r * r * h
                + r * kappa * (top * top - z * z)
                + kappa * kappa * (top.powi(3) - z.powi(3)) / 3.0);
        z = top;
    }
    (eps, kappa, AREA * energy)
}
fn assert_motion(state: &LaminateShellSnapshot, thicknesses: &[f64], temperatures: &[f64]) {
    let (eps, kappa, energy) = clt(thicknesses, temperatures);
    let mesh = discretization().mesh;
    for (node, p) in mesh.nodes.iter().enumerate() {
        let expected = [
            eps * p[0],
            eps * p[1],
            -0.5 * kappa * (p[0] * p[0] + p[1] * p[1]),
            -kappa * p[1],
            kappa * p[0],
            0.0,
        ];
        for (i, want) in expected.into_iter().enumerate() {
            close(
                state.response().displacement_6n[6 * node + i],
                want,
                2e-8 * want.abs().max(0.001),
            );
        }
    }
    close(
        state.response().strain_energy_j,
        energy,
        2e-8 * energy.abs().max(1e-8),
    );
}

#[test]
fn g1_implicit_heat_exchange_drives_original_temperature_dependent_bimetal() {
    let cards = [card("lower"), card("upper")];
    let body = prepare_body(&cards, &[1, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[300.0, 340.0]);
    let old_cells = initial.cells().clone();
    let old_links = initial.links().clone();
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let dt = 0.5;
    let step = owner
        .advance_to(
            &initial,
            dt,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    // Two half-layer distances in series: G=A*k/(h/2), C=rho*A*(h/2)*cp.
    let g = AREA * K / (H / 2.0);
    let c = RHO * AREA * (H / 2.0) * CP;
    let difference = 40.0 / (1.0 + 2.0 * g * dt / c);
    let expected = [320.0 - difference / 2.0, 320.0 + difference / 2.0];
    for (i, &t) in expected.iter().enumerate() {
        close(step.shell.temperatures_k()[i], t, 5e-10);
        let layer = &step.shell.material().layers()[i];
        close(layer.elasticity().state_coordinate("T").unwrap(), t, 5e-10);
        close(
            layer
                .elasticity()
                .property("young_modulus")
                .unwrap()
                .value_si(),
            young(t),
            1e-5,
        );
        close(
            step.shell.material().prepared().profiles()[i].faces()[0][0],
            free(t),
            1e-14,
        );
        assert_eq!(layer.elasticity().card_identity(), cards[i].content_hash());
    }
    assert_motion(&step.shell, &[H / 2.0; 2], &expected);
    assert!(step.shell.response().strain_energy_j > 0.0);
    close(step.thermal.mass_kg(), RHO * H * AREA, 1e-12);
    let section = step.shell.material().prepared().section();
    let reference = owner.source().reference_section();
    assert_eq!(section.i0_kg_m2().to_bits(), reference.i0_kg_m2().to_bits());
    assert_eq!(section.i1_kg_m().to_bits(), reference.i1_kg_m().to_bits());
    assert_eq!(section.i2_kg().to_bits(), reference.i2_kg().to_bits());
    // Heat conservation remains the caloric owner's separate account; stored
    // mechanical energy is measured above and is not inserted into this audit.
    close(step.thermal.energy_j(), initial.energy_j(), 1e-8);
    close(step.heat_report.energy.stored_energy_change_j, 0.0, 1e-8);
    close(step.heat_report.energy.external_heat_j, 0.0, 0.0);
    close(step.heat_report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(step.heat_report.energy.substeps, 1);
    assert_eq!(step.shell.time_s(), dt);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &old_cells);
    assert_eq!(initial.links(), &old_links);
}

#[test]
fn g1_physical_layers_retain_every_unequal_thermal_cell_without_averaging() {
    let cards = [card("two-cell material"), card("three-cell material")];
    let body = prepare_body(&cards, &[2, 3], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let temperatures = [300.0, 340.0, 305.0, 325.0, 350.0];
    let initial = initial(&bound, &temperatures);
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let state = owner.resolve(&initial, &[0.0; 24], || false).unwrap();
    let thicknesses = [H / 4.0, H / 4.0, H / 6.0, H / 6.0, H / 6.0];
    assert_eq!(owner.source().layers().len(), 5);
    for (i, source) in owner.source().layers().iter().enumerate() {
        let physical = usize::from(i >= 2);
        assert!(Arc::ptr_eq(&source.input().card, &cards[physical]));
        close(source.input().thickness_m, thicknesses[i], 1e-18);
        close(state.temperatures_k()[i], temperatures[i], 5e-10);
        close(
            state.material().prepared().profiles()[i].faces()[0][0],
            free(temperatures[i]),
            1e-14,
        );
    }
    assert_motion(&state, &thicknesses, &temperatures);
    let exact = clt(&thicknesses, &temperatures);
    let averaged = clt(&[H / 2.0; 2], &[320.0, (305.0 + 325.0 + 350.0) / 3.0]);
    assert!(
        (exact.1 - averaged.1).abs() > 1e-3,
        "fixture must distinguish actual sublayers from temperature averages"
    );
    close(
        AREA * state.material().prepared().section().i0_kg_m2(),
        initial.mass_kg(),
        1e-12,
    );
}

#[test]
fn g1_uniform_free_expansion_has_near_zero_elastic_strain_bound_and_storage() {
    let original = card("uniform");
    let cards = [Arc::clone(&original), original];
    let body = prepare_body(&cards, &[2, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[330.0; 3]);
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let state = owner.resolve(&initial, &[0.0; 24], || false).unwrap();
    assert_motion(&state, &[H / 4.0, H / 4.0, H / 2.0], &[330.0; 3]);
    assert!(state.reduced_elastic_strain_upper_bound() < 1e-10);
    assert!(state.response().strain_energy_j < 1e-17);
    close(state.response().displacement_6n[6], free(330.0), 1e-11);
    // The card's density falls with temperature, but the original inventory
    // remains exactly rho(300 K) times reference thickness and patch area.
    close(
        AREA * state.material().prepared().section().i0_kg_m2(),
        RHO * H * AREA,
        1e-12,
    );
}

#[test]
fn g0_binding_envelope_and_late_cancellation_leave_original_timeline_unchanged() {
    let cards = [card("binding lower"), card("binding upper")];
    let body = prepare_body(&cards, &[1, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[300.0, 340.0]);
    let old_cells = initial.cells().clone();
    let old_links = initial.links().clone();
    let mut wrong_card = mechanical(&cards);
    wrong_card[0].card = card("different immutable source");
    assert!(matches!(
        LayeredThermalShell::try_new(
            &bound,
            wrong_card,
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false
        ),
        Err(LaminateShellError::Invalid(_))
    ));
    let mut wrong_thickness = mechanical(&cards);
    wrong_thickness[0].thickness_m *= 1.1;
    assert!(matches!(
        LayeredThermalShell::try_new(
            &bound,
            wrong_thickness,
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false
        ),
        Err(LaminateShellError::Invalid(_))
    ));
    let wrong_mass_body = prepare_body(&cards, &[1, 1], 1.1);
    let wrong_mass_laws = wrong_mass_body.thermal_laws(false).unwrap();
    let wrong_mass_bound = wrong_mass_body.bind(&wrong_mass_laws, ids()).unwrap();
    assert!(matches!(
        LayeredThermalShell::try_new(
            &wrong_mass_bound,
            mechanical(&cards),
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false
        ),
        Err(LaminateShellError::Invalid(
            "thermal and mechanical reference cell mass differ"
        ))
    ));
    let restrictive = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        LaminateShellLimits {
            maximum_translation_m: 1e-8,
            ..limits()
        },
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    assert!(matches!(
        restrictive.advance_to(
            &initial,
            0.5,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false
        ),
        Err(LaminateShellError::Envelope {
            quantity: "translation",
            ..
        })
    ));
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let mut checkpoints = 0;
    let accepted = owner
        .advance_to(
            &initial,
            0.5,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(accepted.thermal.time_s(), 0.5);
    let mut retry = 0;
    // Cancel at the final checkpoint, after both physical candidates were
    // computed: neither candidate may alter the input owner's accepted state.
    assert!(matches!(
        owner.advance_to(
            &initial,
            0.5,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || {
                retry += 1;
                retry == checkpoints
            }
        ),
        Err(LaminateShellError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &old_cells);
    assert_eq!(initial.links(), &old_links);
}
