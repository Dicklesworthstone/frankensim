//! Same-card original thermal hardening observations on the reciprocal caloric owner.
//! The independent coaxial oracle includes the table's exact storage and entropy.

#[path = "tabulated/pressure.rs"]
mod pressure;

use super::*;
use fs_matdb::ClaimSelection;
use fs_material::plastic::finite::thermal::source::tabulated::THERMAL_FINITE_J2_TABULATED_LAW_ID;
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS, KIRCHHOFF_YIELD_CURVE_PROPERTY,
};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const MU: f64 = 80.0e6;
const KNOTS: [(f64, f64); 5] = [
    (0.0, 2.0e6),
    (0.04, 5.0e6),
    (0.12, 6.0e6),
    (0.3, 6.0e6),
    (0.8, 20.0e6),
];
const SLOPES: [f64; 5] = [-4_000.0, -10_000.0, -14_000.0, -14_000.0, -35_000.0];
const SAMPLE_TEMPERATURES: [f64; 3] = [280.0, 300.0, 420.0];
const SAMPLE_LIMIT: usize = KNOTS.len() * SAMPLE_TEMPERATURES.len();

fn surface() -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(KIRCHHOFF_YIELD_CURVE_PROPERTY, PA),
        value: PropertyValue::Surface2D {
            axes: [ACCUMULATED_PLASTIC_STRAIN_AXIS.into(), "temperature".into()],
            axis_dims: [Dims::NONE, K],
            samples: SAMPLE_TEMPERATURES
                .into_iter()
                .flat_map(|temperature| {
                    KNOTS
                        .into_iter()
                        .zip(SLOPES)
                        .map(move |((alpha, y), slope)| {
                            ([alpha, temperature], y + slope * (temperature - 300.0))
                        })
                })
                .collect(),
            triangles: vec![],
            dims: PA,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                ACCUMULATED_PLASTIC_STRAIN_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                0.8,
            )
            .with_quantity("temperature", temperature_quantity(), 280.0, 420.0),
        interpolation: InterpolationPolicy::TabulatedOnly,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    }
}

fn sourced_fixture(identity: usize) -> Fixture {
    let mut constitutive = model();
    constitutive.law = LawId(THERMAL_FINITE_J2_TABULATED_LAW_ID.into());
    for name in [
        "yield_stress_pa",
        "hardening_modulus_pa",
        "yield_stress_slope_pa_k",
        "hardening_modulus_slope_pa_k",
    ] {
        assert!(constitutive.parameters.remove(name).is_some());
    }
    // Constant elastic coefficients isolate entropy from the measured hardening
    // surface. Thermal expansion and elastic coefficient derivatives are covered
    // by the parent regression; this test has a fully independent closed form.
    for name in [
        "shear_modulus_slope_pa_k",
        "bulk_modulus_slope_pa_k",
        "logarithmic_expansion_per_k",
        "logarithmic_expansion_slope_per_k2",
    ] {
        constitutive.parameters.get_mut(name).unwrap().value = 0.0;
    }
    constitutive
        .parameters
        .get_mut("maximum_accumulated_plastic_strain")
        .unwrap()
        .value = 0.8;
    let mut claims = ClaimSet::new();
    claims.insert_claim(scalar(CP_NAME, CP_DIMS, CP)).unwrap();
    claims
        .insert_claim(scalar(DENSITY_PROPERTY, Density::DIMS, RHO))
        .unwrap();
    claims.insert_claim(surface()).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic thermal tabulated coupled fixture {identity}"),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        vec![constitutive],
    )
    .unwrap();
    // Resolve only after adding the surface, so all three source owners retain
    // the identity of the same complete original MaterialCard.
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "temperature",
            lower: &source_context(280.0),
            upper: &source_context(420.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1.0e-6,
    )
    .unwrap();
    Fixture { card, caloric }
}

fn ordinates(temperature: f64) -> [f64; 5] {
    core::array::from_fn(|i| KNOTS[i].1 + SLOPES[i] * (temperature - 300.0))
}

fn interpolate(alpha: f64, values: [f64; 5]) -> f64 {
    assert!((0.0..=0.8).contains(&alpha));
    for i in 0..KNOTS.len() - 1 {
        if alpha <= KNOTS[i + 1].0 {
            return values[i]
                + (values[i + 1] - values[i]) * (alpha - KNOTS[i].0)
                    / (KNOTS[i + 1].0 - KNOTS[i].0);
        }
    }
    unreachable!()
}

/// Exact integral of the piecewise-linear ordinates relative to initial yield.
fn primitive(alpha: f64, values: [f64; 5]) -> f64 {
    assert!((0.0..=0.8).contains(&alpha));
    let mut integral = 0.0;
    for i in 0..KNOTS.len() - 1 {
        if alpha <= KNOTS[i].0 {
            break;
        }
        let dx = alpha.min(KNOTS[i + 1].0) - KNOTS[i].0;
        integral += dx
            * (values[i] - values[0]
                + 0.5 * dx * (values[i + 1] - values[i]) / (KNOTS[i + 1].0 - KNOTS[i].0));
    }
    integral
}

fn coaxial(axial_log: f64) -> [f64; 9] {
    [
        axial_log.exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * axial_log).exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * axial_log).exp(),
    ]
}

fn displacement(f: &[f64; 9]) -> [[f64; 3]; 4] {
    core::array::from_fn(|node| {
        core::array::from_fn(|axis| (f[3 * axis + axis] - 1.0) * NODES[node][axis])
    })
}

struct Expected {
    time_s: f64,
    temperature: f64,
    old_temperature: f64,
    alpha: f64,
    deformation: [f64; 9],
    displacement: [[f64; 3]; 4],
    old_displacement: [[f64; 3]; 4],
    first_piola: [f64; 9],
    free_energy: f64,
    entropy: f64,
    internal_energy: f64,
    internal_change: f64,
    plastic_heat: f64,
    reversible_heat: f64,
    work_defect: f64,
    boundary_heat_j: f64,
}

impl Expected {
    fn prescribed(&self) -> Vec<(usize, f64)> {
        self.displacement
            .iter()
            .flatten()
            .enumerate()
            .filter(|(dof, _)| *dof != 3)
            .map(|(dof, value)| (dof, *value))
            .collect()
    }
    fn free_force_n(&self) -> f64 {
        V / L * self.first_piola[0]
    }
    fn caloric_heat(&self) -> f64 {
        self.plastic_heat + self.reversible_heat
    }
}

fn oracle_path() -> Vec<Expected> {
    let mut old_temperature = 300.0;
    let mut old_alpha = 0.0;
    let mut old_f = IDENTITY;
    let mut old_free = 0.0;
    let mut old_entropy = 0.0;
    let mut path = Vec::new();
    for (index, temperature) in [310.0, 400.0, 300.0].into_iter().enumerate() {
        // Hold total F while heating across the .12 hardening knot, then remove
        // elastic strain at the cold reference temperature, preserving Cp.
        let axial_log = if index == 2 { old_alpha } else { 0.142 };
        let y = ordinates(temperature);
        let mut alpha = old_alpha;
        if 3.0 * MU * (axial_log - old_alpha) > interpolate(old_alpha, y) {
            let mut lower = old_alpha;
            let mut upper = axial_log;
            for _ in 0..100 {
                let middle = 0.5 * (lower + upper);
                if interpolate(middle, y) < 3.0 * MU * (axial_log - middle) {
                    lower = middle;
                } else {
                    upper = middle;
                }
            }
            alpha = 0.5 * (lower + upper);
        }
        let q = 3.0 * MU * (axial_log - alpha);
        let f = coaxial(axial_log);
        let mut piola = [0.0; 9];
        piola[0] = (2.0 / 3.0) * q / f[0];
        piola[4] = (-1.0 / 3.0) * q / f[4];
        piola[8] = (-1.0 / 3.0) * q / f[8];
        let free_energy = q * q / (6.0 * MU) + primitive(alpha, y);
        let entropy = -primitive(alpha, SLOPES);
        let internal_energy = free_energy + temperature * entropy;
        let coefficient_change = (temperature - old_temperature) * primitive(old_alpha, SLOPES);
        let entropy_storage_change = temperature * entropy - old_temperature * old_entropy;
        let plastic_heat = y[0] * (alpha - old_alpha);
        let reversible_heat = -coefficient_change - entropy_storage_change;
        let work: f64 = (0..9).map(|i| piola[i] * (f[i] - old_f[i])).sum();
        let internal_change = free_energy - old_free + entropy_storage_change;
        let work_defect = work - (free_energy - old_free - coefficient_change) - plastic_heat;
        close(
            work - internal_change - plastic_heat - reversible_heat - work_defect,
            0.0,
            1e-8,
        );
        path.push(Expected {
            time_s: (index + 1) as f64,
            temperature,
            old_temperature,
            alpha,
            deformation: f,
            displacement: displacement(&f),
            old_displacement: displacement(&old_f),
            first_piola: piola,
            free_energy,
            entropy,
            internal_energy,
            internal_change,
            plastic_heat,
            reversible_heat,
            work_defect,
            boundary_heat_j: RHO * V * CP * (temperature - old_temperature)
                - V * (plastic_heat + reversible_heat),
        });
        old_temperature = temperature;
        old_alpha = alpha;
        old_f = f;
        old_free = free_energy;
        old_entropy = entropy;
    }
    assert!(path[0].alpha > 0.04 && path[0].alpha < 0.12);
    assert!(path[1].alpha > 0.12);
    assert_eq!(path[2].alpha, path[1].alpha);
    path
}

fn initial<'a>(fixture: &'a Fixture, path: &[Expected]) -> DepositionSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CELL,
            material_id: MaterialId(0),
            segment_id: 1,
            body_id: 7,
            mass_kg: RHO * V,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [L; 3],
                volume_m3: V,
                surface_area_m2: 3.0 * L * L,
            },
            boundary: boundary(path[0].boundary_heat_j),
        }),
    }];
    for (i, endpoint) in path.iter().enumerate().skip(1) {
        events.push(TimedEvent {
            time_s: i as f64,
            event: Event::Boundary(CELL, boundary(endpoint.boundary_heat_j)),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn bind<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    cx: &Cx<'_>,
) -> FiniteThermoplasticSolid<'a> {
    FiniteThermoplasticSolid::try_new_with_tabulated_hardening(
        thermal,
        &NODES,
        &TETS,
        vec![FiniteThermoplasticCell {
            cell: CELL,
            query: fixture.query(),
        }],
        SINGLE,
        SAMPLE_LIMIT,
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

#[test]
fn original_thermal_surface_solves_heat_driven_knot_crossing_and_retained_unloaded_shape() {
    context(|cx| {
        let fixture = sourced_fixture(0);
        let path = oracle_path();
        let mut thermal = initial(&fixture, &path);
        let mut solid = bind(&fixture, &thermal, cx);
        for (index, expected) in path.iter().enumerate() {
            let candidate = solid
                .advance(
                    &thermal,
                    &(),
                    expected.time_s,
                    &expected.prescribed(),
                    &[(3, expected.free_force_n())],
                    [0.0; 3],
                    settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let actual = candidate.solid().mechanics();
            let report = candidate.report();
            close(
                candidate.thermal().temperature_k(CELL).unwrap(),
                expected.temperature,
                2e-7,
            );
            close(actual.temperatures_k()[0], expected.temperature, 2e-7);
            for (u, want) in actual
                .displacement_m()
                .iter()
                .flatten()
                .zip(expected.displacement.iter().flatten())
            {
                close(*u, *want, 2e-11);
            }
            for (f, want) in actual.element_deformation()[0]
                .iter()
                .zip(expected.deformation)
            {
                close(*f, want, 2e-9);
            }
            close(actual.element_states()[0].alpha, expected.alpha, 3e-9);
            let inverse_cp = [
                (-2.0 * expected.alpha).exp(),
                expected.alpha.exp(),
                expected.alpha.exp(),
                0.0,
                0.0,
                0.0,
            ];
            for (metric, want) in actual.element_states()[0]
                .inverse_plastic_cauchy_green
                .iter()
                .zip(inverse_cp)
            {
                close(*metric, want, 4e-9);
            }
            for (piola, want) in actual.element_first_piola_pa()[0]
                .iter()
                .zip(expected.first_piola)
            {
                close(*piola, want, 0.05);
            }
            close(actual.free_energy_j(), V * expected.free_energy, 2e-9);
            close(
                actual.internal_energy_j(),
                V * expected.internal_energy,
                2e-9,
            );
            close(
                report.material_energy.plastic_dissipation_j,
                V * expected.plastic_heat,
                2e-9,
            );
            close(
                report.material_energy.reversible_heat_j,
                V * expected.reversible_heat,
                2e-9,
            );
            close(
                report.material_energy.caloric_heat_j,
                V * expected.caloric_heat(),
                2e-9,
            );
            close(
                report.material_energy.internal_energy_change_j,
                V * expected.internal_change,
                2e-9,
            );
            close(
                report.material_energy.endpoint_work_defect_j,
                V * expected.work_defect,
                2e-9,
            );
            close(report.external_heat_j, expected.boundary_heat_j, 2e-9);
            close(
                report.thermal_audit.stored_energy_change_j,
                RHO * V * CP * (expected.temperature - expected.old_temperature),
                2e-8,
            );
            close(
                report.thermal_audit.stored_energy_change_j - report.external_heat_j,
                report.material_energy.caloric_heat_j,
                2e-9,
            );
            let mut total_work = expected.free_force_n()
                * (expected.displacement[1][0] - expected.old_displacement[1][0]);
            for &(dof, reaction) in &report.equilibrium.reactions_n {
                total_work += reaction
                    * (expected.displacement[dof / 3][dof % 3]
                        - expected.old_displacement[dof / 3][dof % 3]);
            }
            close(
                total_work + report.external_heat_j,
                report.material_energy.internal_energy_change_j
                    + report.thermal_audit.stored_energy_change_j
                    + report.material_energy.endpoint_work_defect_j,
                2e-8,
            );
            close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            assert!(report.iterations > 1);
            assert!(report.inertia.is_none());

            let material = &actual.families()[0];
            assert!(std::ptr::eq(material.caloric_source(), &fixture.caloric));
            assert_eq!(
                material.query().card.content_hash(),
                fixture.card.content_hash()
            );
            assert_eq!(material.family().model_card(), &fixture.card.models()[0]);
            let source = material
                .family()
                .tabulated_hardening_source()
                .unwrap()
                .source();
            assert_eq!(source.card_identity(), fixture.card.content_hash());
            assert_eq!(source.claim(), &surface());
            assert_eq!(source.context(), &source_context(300.0));
            assert_eq!(source.selection(), SINGLE);
            assert_eq!(source.observations().len(), SAMPLE_LIMIT);
            for (_, answer) in source.observations() {
                fixture
                    .card
                    .claims()
                    .verify_receipt(&answer.receipt)
                    .unwrap();
            }
            let inspection = material
                .family()
                .law()
                .thermodynamics(
                    &actual.element_deformation()[0],
                    &actual.element_states()[0],
                    actual.temperatures_k()[0],
                )
                .unwrap();
            close(inspection.entropy_j_m3_k, expected.entropy, 3e-5);
            // The model's interpolation does not promote original observations
            // into an interpolated source claim at an unmeasured coordinate.
            let intermediate = source_context(400.0)
                .with_quantity(
                    ACCUMULATED_PLASTIC_STRAIN_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    expected.alpha,
                )
                .unwrap();
            assert!(
                fixture
                    .card
                    .claims()
                    .query_typed(
                        &source.claim().key,
                        &intermediate,
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .is_err()
            );
            if index == 1 {
                close(report.material_energy.internal_endpoint_work_j, 0.0, 2e-8);
                assert!(report.material_energy.plastic_dissipation_j > 1e-5);
                assert!(expected.reversible_heat.abs() > expected.plastic_heat);
            }
            if index == 2 {
                close(report.material_energy.plastic_dissipation_j, 0.0, 2e-9);
                assert!(actual.displacement_m()[1][0] > 1e-3);
                assert!(actual.displacement_m()[2][1] < -5e-4);
            }
            assert_eq!(solid.time_s(), index as f64);
            assert_eq!(thermal.time_s(), index as f64);
            let (next_thermal, next_solid, (), _) = candidate.into_parts();
            thermal = next_thermal;
            solid = next_solid;
        }
    });
}

#[derive(Clone)]
struct ExcessHeat;
impl DepositionObserver for ExcessHeat {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok([(CELL, 1_000.0)].into())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn foreign_caloric_evidence_and_unsupported_heat_or_history_preserve_the_accepted_pair() {
    context(|cx| {
        let fixture = sourced_fixture(0);
        let foreign = sourced_fixture(1);
        let mut wrong = fixture.query();
        wrong.caloric = &foreign.caloric;
        assert!(matches!(
            FiniteThermoplasticMaterial::try_new_with_tabulated_hardening(
                wrong,
                0.0,
                SINGLE,
                SAMPLE_LIMIT,
            ),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        assert!(FiniteThermoplasticMaterial::try_new(fixture.query(), 0.0).is_err());
        assert!(
            FiniteThermoplasticMaterial::try_new_with_tabulated_hardening(
                fixture.query(),
                0.0,
                SINGLE,
                SAMPLE_LIMIT - 1,
            )
            .is_err()
        );
        let path = oracle_path();
        let thermal = initial(&fixture, &path);
        let solid = bind(&fixture, &thermal, cx);
        let first = &path[0];
        let accepted = solid
            .advance(
                &thermal,
                &(),
                first.time_s,
                &first.prescribed(),
                &[(3, first.free_force_n())],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, solid, (), _) = accepted.into_parts();
        let cells = thermal.cells().clone();
        let history = solid.mechanics().element_states().to_vec();
        let deformation = solid.mechanics().element_deformation().to_vec();
        let motion = solid.mechanics().displacement_m().to_vec();
        let stress = solid.mechanics().element_first_piola_pa().to_vec();
        let temperature = solid.mechanics().temperatures_k().to_vec();
        let free = solid.mechanics().free_energy_j();
        let internal = solid.mechanics().internal_energy_j();
        let next = &path[1];
        assert!(
            solid.mechanics().families()[0]
                .family()
                .law()
                .at_temperature(420.0001)
                .is_err()
        );
        // Thermal source power cannot extrapolate the same-card caloric chart.
        assert!(
            solid
                .advance(
                    &thermal,
                    &ExcessHeat,
                    next.time_s,
                    &next.prescribed(),
                    &[(3, next.free_force_n())],
                    [0.0; 3],
                    settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .is_err()
        );
        // Complete prescribed deformation demands alpha > .8. Elastic and
        // tensor-metric bounds permit the state; the measured flow support does not.
        let outside: Vec<_> = displacement(&coaxial(1.0))
            .into_iter()
            .flatten()
            .enumerate()
            .collect();
        assert!(
            solid
                .advance(
                    &thermal,
                    &(),
                    next.time_s,
                    &outside,
                    &[],
                    [0.0; 3],
                    settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(solid.mechanics().element_states(), history);
        assert_eq!(solid.mechanics().element_deformation(), deformation);
        assert_eq!(solid.mechanics().displacement_m(), motion);
        assert_eq!(solid.mechanics().element_first_piola_pa(), stress);
        assert_eq!(solid.mechanics().temperatures_k(), temperature);
        assert_eq!(solid.mechanics().free_energy_j(), free);
        assert_eq!(solid.mechanics().internal_energy_j(), internal);
        assert_eq!(thermal.time_s(), 1.0);
        assert_eq!(solid.time_s(), 1.0);
        let retry = solid
            .advance(
                &thermal,
                &(),
                next.time_s,
                &next.prescribed(),
                &[(3, next.free_force_n())],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(
            retry.thermal().temperature_k(CELL).unwrap(),
            next.temperature,
            2e-7,
        );
        close(
            retry.solid().mechanics().element_states()[0].alpha,
            next.alpha,
            3e-9,
        );
        assert!(retry.solid().mechanics().element_states()[0].alpha > 0.12);
        close(
            retry.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(solid.mechanics().element_states(), history);
    });
}
