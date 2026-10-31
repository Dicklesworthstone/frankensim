//! Manufactured finite coolant/J2 feedback on the original shared thermal graph.
#[path = "solid/electrical.rs"]
mod electrical;
use super::*;
use fs_couple::solid_enthalpy::plastic::{
    ThermoplasticCell, ThermoplasticConfig, ThermoplasticSolid,
};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::thermal::SourceQuery;
use fs_material::state_point::{
    POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [0.01, 0.0, 0.0],
    [0.0, 0.02, 0.0],
    [0.0, 0.0, 0.03],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const VOLUME: f64 = 1e-6;
const SHEAR: f64 = 0.01;
const TRACTION: f64 = 0.5e6;
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;

fn solid_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let validity =
        || ValidityDomain::unconstrained().with_quantity("T", temperature_quantity(), 280.0, 320.0);
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        // This separate manufactured mechanical specimen declares T-only
        // applicability. The parent fixture's untyped lot applicability is
        // not invented into a pressure/composition law to bypass admission.
        claim.validity = validity();
        claims.insert_claim(claim).unwrap();
    }
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2),
        (BASIS.property_name(), Pressure::DIMS, 1e8),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6),
    ] {
        let strength = name == YIELD_STRESS_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if strength {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(280.0, 5e6), (300.0, 3e6), (320.0, 1e6)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: validity(),
                interpolation: if strength {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite-coolant thermoplastic specimen; not measured data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &body_point(280.0),
            upper: &body_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, caloric)
}

fn supports(shear: f64) -> Vec<(usize, f64)> {
    // One actual equilibrium DOF remains: the apex's x displacement.
    (0..12)
        .filter(|&d| d != 9)
        .map(|d| {
            (
                d,
                if d % 3 == 0 {
                    2.0 * shear * NODES[d / 3][1]
                } else {
                    0.0
                },
            )
        })
        .collect()
}

fn solid_controls() -> ThermoplasticConfig {
    ThermoplasticConfig {
        maximum_iterations: 80,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        enthalpy_relaxation: 1.0,
    }
}

// Independent two-shear radial-return reduction: prescribed epsilon_xy,
// actual apex traction sigma_zx and scalar energetic hardening.
fn apex_reference(temperature_k: f64) -> (f64, f64) {
    let fy = 3e6 - 1e5 * (temperature_k - 300.0);
    let g = 5e8;
    let h = 1e8;
    let q = |gamma: f64| 2.0 * g * 3.0_f64.sqrt() * (SHEAR * SHEAR + gamma * gamma).sqrt();
    let shear_stress = |gamma: f64| 2.0 * g * gamma * (h + 3.0 * g * fy / q(gamma)) / (3.0 * g + h);
    let (mut low, mut high) = (0.0, 0.03);
    assert!(shear_stress(high) > TRACTION);
    for _ in 0..80 {
        let gamma = 0.5 * (low + high);
        if shear_stress(gamma) < TRACTION {
            low = gamma;
        } else {
            high = gamma;
        }
    }
    let gamma = 0.5 * (low + high);
    (2.0 * NODES[3][2] * gamma, (q(gamma) - fy) / (3.0 * g + h))
}

#[test]
fn finite_coolant_data_reaches_native_solid_equilibrium_heat_and_retained_history() {
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 173,
                kernel_id: 0x534f4c49,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let mut observations = Vec::new();
        let mut body_ids = Vec::new();
        for fluid_cp in [2000.0, 8000.0] {
            let (body_card, body_curve) = solid_material();
            let (fluid_card, fluid_curve) = fluid_material(1.0, fluid_cp);
            let thermal = initial_pair(&body_curve, &fluid_curve, fluid_cp, false, 0.0);
            let original_cells = thermal.cells().clone();
            let original_links = thermal.links().clone();
            let fluid_context = fluid_point(290.0);
            let fluid = liquid_law(&fluid_card);
            let solid = ThermoplasticSolid::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![ThermoplasticCell {
                    cell: BODY,
                    query: SourceQuery {
                        card: &body_card,
                        caloric: &body_curve,
                        context: body_point(310.0),
                        temperature_axis: "T".into(),
                        selection: MaterialPropertySelection::SingleClaimOnly,
                        hardening_basis: BASIS,
                        thermal_expansion: None,
                        strain_limit: 0.04,
                    },
                }],
                Default::default(),
                &cx,
            )
            .unwrap();
            let original_solid = solid.solution().clone();
            let roles = [None];
            let mut cells = solid.thermal_material_cells(&roles).unwrap();
            cells.push(MaterialCell {
                cell: FLUID,
                material_id: MaterialId(1),
                curve: fluid_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &fluid_card,
                    context: &fluid_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            });
            let bulk = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &[],
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let links = [binding()];
            let initial = evaluate_source_link_convection(
                &thermal,
                &links,
                1,
                |_, _, t| fluid.at_temperature(t),
                || false,
            )
            .unwrap();
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let assembly = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let thermal_config = StepConfig {
                max_dt_s: DT / 10.0,
                energy_tolerance_j: 1e-9,
            };
            let (loaded, audit) = solid
                .advance_heat_assembly_candidate(
                    &thermal,
                    &assembly,
                    &(),
                    DT,
                    &supports(SHEAR),
                    &[(9, TRACTION * VOLUME / NODES[3][2])],
                    Default::default(),
                    thermal_config,
                    solid_controls(),
                    &cx,
                )
                .unwrap();
            let wall = loaded.thermal().temperature_k(BODY).unwrap();
            let coolant = loaded.thermal().temperature_k(FLUID).unwrap();
            let expected = apex_reference(wall);
            let solution = loaded.solid().solution();
            let report = loaded.report();
            assert!(wall < 310.0 && wall > coolant && coolant > 290.0);
            near(solution.displacement_m[3][0], expected.0, 2e-11);
            near(solution.element_states[0].alpha, expected.1, 2e-10);
            near(solution.element_stress[0][5], TRACTION, 1e-2);
            assert!(solution.element_states[0].alpha > 0.0);
            assert!(report.physical_heat_j > 0.0);
            assert!(report.numerical_loss_j > 0.0);
            near(report.reversible_heat_j, 0.0, 1e-12);
            near(audit.external_heat_j, 0.0, 1e-10);
            near(
                audit.thermal_residual_j,
                0.0,
                solid_controls().energy_tolerance_j,
            );
            near(
                audit.combined_residual_j,
                0.0,
                solid_controls().energy_tolerance_j,
            );
            let stored = loaded.thermal().energy_j() - thermal.energy_j();
            near(
                stored,
                report.physical_heat_j,
                solid_controls().energy_tolerance_j,
            );
            near(
                MASS * CP * (wall - 310.0) + MASS * fluid_cp * (coolant - 290.0),
                report.physical_heat_j,
                solid_controls().energy_tolerance_j,
            );
            near(
                stored + report.stored_mechanical_change_j + report.numerical_loss_j,
                report.mechanical_work_j,
                solid_controls().energy_tolerance_j,
            );

            let sample = loaded.observer().first.sample().unwrap();
            let finite = sample.finite_convection.as_ref().unwrap();
            let link = &finite.samples()[0];
            assert_eq!((link.wall_cell, link.fluid_cell), (BODY, FLUID));
            near(link.wall_temperature_k, wall, 1e-12);
            near(link.fluid_temperature_k, coolant, 1e-12);
            near(link.property_temperature_k, 0.5 * (wall + coolant), 1e-12);
            assert_eq!(
                link.material.material().card_identity(),
                fluid_card.content_hash()
            );
            near(
                link.material
                    .material()
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .unwrap()
                    .value_si(),
                fluid_cp,
                1e-12,
            );
            for property in link.material.material().properties() {
                assert_eq!(
                    property.answer().receipt.selected,
                    fluid
                        .reference_state()
                        .material()
                        .property(property.requirement().name())
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
            assert_eq!(link.material.material().query().axes()["p"], 100000.0);
            assert_eq!(
                link.material.material().query().axes()["solute_mass_fraction"],
                0.2
            );
            let structural = loaded.solid().materials()[0].material();
            assert_eq!(structural.card_identity(), body_card.content_hash());
            for property in structural.properties() {
                assert!(
                    property
                        .answer()
                        .receipt
                        .query_point
                        .contains(&("T".into(), wall))
                );
                assert_eq!(
                    property.answer().receipt.selected,
                    solid.materials()[0]
                        .material()
                        .property(property.requirement().name())
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
            for body in &sample.heat.material.cells {
                near(
                    body.phase.temperature_k(),
                    loaded.thermal().temperature_k(body.cell).unwrap(),
                    1e-12,
                );
            }
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(thermal.links(), &original_links);
            assert_eq!(solid.solution(), &original_solid);
            near(loaded.solid().time_s(), DT, 1e-16);
            near(loaded.thermal().time_s(), DT, 1e-16);
            body_ids.push(body_card.content_hash());
            observations.push([
                wall,
                coolant,
                solution.displacement_m[3][0],
                solution.element_states[0].alpha,
            ]);

            // Continue the original finite coolant and loaded solid. Elastic
            // unloading retains the accepted plastic strain/alpha, while the
            // same source assembly keeps exchanging internal thermal energy.
            let history = solution.element_states.clone();
            let (unloaded, audit) = loaded
                .solid()
                .advance_heat_assembly_candidate(
                    loaded.thermal(),
                    &assembly,
                    &loaded.observer().second,
                    2.0 * DT,
                    &supports(SHEAR - 0.0005),
                    &[],
                    Default::default(),
                    thermal_config,
                    solid_controls(),
                    &cx,
                )
                .unwrap();
            assert_eq!(unloaded.solid().solution().element_states, history);
            near(
                unloaded.solid().solution().displacement_m[3][0],
                2.0 * NODES[3][2] * history[0].plastic_strain[5],
                2e-11,
            );
            near(unloaded.report().physical_heat_j, 0.0, 1e-10);
            near(
                unloaded.thermal().energy_j() - loaded.thermal().energy_j(),
                0.0,
                1e-8,
            );
            near(audit.external_heat_j, 0.0, 1e-10);
            near(
                audit.combined_residual_j,
                0.0,
                solid_controls().energy_tolerance_j,
            );
            assert!(unloaded.thermal().temperature_k(BODY).unwrap() < wall);
            assert!(unloaded.thermal().temperature_k(FLUID).unwrap() > coolant);
            near(unloaded.solid().time_s(), 2.0 * DT, 1e-16);
            assert_eq!(loaded.solid().solution().element_states, history);
            assert_eq!(solid.solution(), &original_solid);
        }
        assert_eq!(body_ids[0], body_ids[1]);
        // Only the original coolant Cp data changed: it changes BOTH finite
        // coolant storage and transport, then actual mechanical equilibrium.
        assert!(observations[1][0] < observations[0][0] - 0.1);
        assert!(observations[1][1] < observations[0][1]);
        assert!(observations[1][2] < observations[0][2]);
        assert!(observations[1][3] < observations[0][3]);
    });
}
