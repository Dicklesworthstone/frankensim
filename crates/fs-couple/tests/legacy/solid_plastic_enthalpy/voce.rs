//! G1/G3 heated loaded spatial member with independently sourced Voce softening.
//! The constitutive and caloric numbers are synthetic numerical fixtures.
use super::*;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
};

const Q0: f64 = 3e6;
const B_VALUE: f64 = 120.0;
const CP_VALUE: f64 = 100.0;
const DURATION: f64 = 0.001;
const HEAT_PER_VOLUME: f64 = RHO * CP_VALUE * 4.0;
const STRAIN: f64 = 0.006;

fn fixture(slope: f64) -> Fixture {
    fixture_with_rate(slope, 1e9, 1.0)
}

fn fixture_with_rate(slope: f64, viscosity_pa_s: f64, exponent: f64) -> Fixture {
    let original = Fixture::new(CP_VALUE);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.card.claims().claims_ordered() {
        // Hold the initial yield offset fixed so only the selected saturation
        // data cause further plastic flow during the heated strain hold.
        if claim.key.name() != YIELD_STRESS_PROPERTY {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    let mut saturation = scalar(Q, Pressure::DIMS, Q0);
    saturation.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: [280.0, 300.0, 320.0]
            .map(|temperature| (temperature, Q0 + slope * (temperature - 300.0)))
            .to_vec(),
        dims: Pressure::DIMS,
    };
    saturation.interpolation = InterpolationPolicy::LinearInside;
    saturation.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    for claim in [
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD),
        saturation,
        scalar(B, Dims([0; 6]), B_VALUE),
        scalar(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            viscosity_pa_s,
        ),
        scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 1e6),
        scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims([0; 6]), exponent),
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let mut identity = original.card.id().clone();
    identity.chemistry = "synthetic temperature-dependent Voce solid".into();
    let card = MaterialCard::assemble(identity, claims, vec![]).unwrap();
    let lower = QueryPoint::new().with("T", 280.0).unwrap();
    let upper = QueryPoint::new().with("T", 320.0).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP_VALUE,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

fn member() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let (mut nodes, tetrahedra) = cube();
    for node in &mut nodes {
        node[0] *= 3.0;
    }
    (nodes, tetrahedra)
}

fn initial_member<'a>(fixture: &'a Fixture) -> DepositionSimulation<'a> {
    // The six connected tetrahedra partition a 3L by L by L prismatic member.
    // Each native tet therefore has volume 3V. The structural binding checks it.
    let (nodes, tetrahedra) = member();
    let events = tetrahedra
        .iter()
        .enumerate()
        .map(|(index, tet)| {
            let area = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
                .into_iter()
                .map(|face| {
                    let a: [f64; 3] =
                        core::array::from_fn(|i| nodes[tet[face[1]]][i] - nodes[tet[face[0]]][i]);
                    let b: [f64; 3] =
                        core::array::from_fn(|i| nodes[tet[face[2]]][i] - nodes[tet[face[0]]][i]);
                    let cross = [
                        a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0],
                    ];
                    0.5 * cross
                        .into_iter()
                        .map(|value| value * value)
                        .sum::<f64>()
                        .sqrt()
                })
                .sum();
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(100 + index as u64),
                    material_id: MaterialId(0),
                    segment_id: index as u64,
                    body_id: 7,
                    mass_kg: RHO * 3.0 * V,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [0.0; 3],
                        end_m: [3.0 * L, L, L],
                        volume_m3: 3.0 * V,
                        // No boundary exchange is selected; only source watts act.
                        surface_area_m2: area,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: HEAT_PER_VOLUME * 3.0 * V / DURATION,
                    },
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 6,
            max_links: 0,
            max_events: 6,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

#[derive(Clone, Copy)]
struct Reference {
    temperature: f64,
    strain: f64,
    alpha: f64,
}

fn primitive(alpha: f64) -> f64 {
    alpha + (-B_VALUE * alpha).exp_m1() / B_VALUE
}

impl Reference {
    fn stress(self) -> f64 {
        E * (self.strain - self.alpha)
    }
    fn q(self, slope: f64) -> f64 {
        Q0 + slope * (self.temperature - 300.0)
    }
    fn entropy(self, slope: f64) -> f64 {
        -slope * primitive(self.alpha)
    }
    fn free_energy(self, slope: f64) -> f64 {
        0.5 * self.stress().powi(2) / E
            + 0.5 * H * self.alpha.powi(2)
            + self.q(slope) * primitive(self.alpha)
    }
    fn at(self, temperature: f64, strain: f64, slope: f64) -> Self {
        let q = Q0 + slope * (temperature - 300.0);
        // Independent free-sided monotone axial equilibrium, with no native
        // return-map call: E(epsilon-alpha)=Y+H alpha+Q(T)(1-exp(-b alpha)).
        let residual =
            |alpha: f64| E * (strain - alpha) - YIELD - H * alpha + q * (-B_VALUE * alpha).exp_m1();
        let alpha = if residual(self.alpha) <= 0.0 {
            self.alpha
        } else {
            let (mut lower, mut upper) = (self.alpha, strain);
            assert!(residual(upper) < 0.0);
            for _ in 0..80 {
                let middle = 0.5 * (lower + upper);
                if residual(middle) > 0.0 {
                    lower = middle;
                } else {
                    upper = middle;
                }
            }
            0.5 * (lower + upper)
        };
        Self {
            temperature,
            strain,
            alpha,
        }
    }
}

fn coefficient(old: Reference, new: Reference, slope: f64) -> f64 {
    (new.q(slope) - old.q(slope)) * primitive(old.alpha)
}
fn entropy_storage(old: Reference, new: Reference, slope: f64) -> f64 {
    new.temperature * new.entropy(slope) - old.temperature * old.entropy(slope)
}
fn reversible(old: Reference, new: Reference, slope: f64) -> f64 {
    -coefficient(old, new, slope) - entropy_storage(old, new, slope)
}
fn endpoint(old: Reference, slope: f64) -> Reference {
    let residual = |temperature| {
        let new = old.at(temperature, STRAIN, slope);
        RHO * CP_VALUE * (temperature - old.temperature)
            - HEAT_PER_VOLUME
            - YIELD * (new.alpha - old.alpha)
            - reversible(old, new, slope)
    };
    let (mut lower, mut upper) = (old.temperature, 320.0);
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (lower + upper);
        if residual(middle) > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    old.at(0.5 * (lower + upper), STRAIN, slope)
}

fn viscous_overstress(delta_alpha: f64, viscosity: f64, exponent: f64) -> f64 {
    1e6 * (viscosity * delta_alpha / (DURATION * 1e6)).powf(1.0 / exponent)
}

fn viscous_endpoint(old: Reference, slope: f64, viscosity: f64, exponent: f64) -> Reference {
    // Independent axial flow equation, solved directly in accumulated plastic
    // strain. The production law instead uses its tensor return and secants.
    let at = |temperature| {
        let q = Q0 + slope * (temperature - 300.0);
        let residual = |alpha: f64| {
            E * (STRAIN - alpha) - YIELD - H * alpha + q * (-B_VALUE * alpha).exp_m1()
                - viscous_overstress(alpha - old.alpha, viscosity, exponent)
        };
        let (mut lower, mut upper) = (old.alpha, STRAIN);
        assert!(residual(lower) > 0.0 && residual(upper) < 0.0);
        for _ in 0..80 {
            let middle = 0.5 * (lower + upper);
            if residual(middle) > 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        Reference {
            temperature,
            strain: STRAIN,
            alpha: 0.5 * (lower + upper),
        }
    };
    let residual = |temperature| {
        let new = at(temperature);
        let da = new.alpha - old.alpha;
        RHO * CP_VALUE * (temperature - old.temperature)
            - HEAT_PER_VOLUME
            - (YIELD + viscous_overstress(da, viscosity, exponent)) * da
            - reversible(old, new, slope)
    };
    let (mut lower, mut upper) = (old.temperature, 320.0);
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (lower + upper);
        if residual(middle) > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    at(0.5 * (lower + upper))
}

#[test]
fn heated_member_relaxes_with_joint_voce_kinetics_and_one_physical_clock() {
    context(|cx| {
        let (nodes, tetrahedra) = member();
        let prescribed = grips(&nodes, STRAIN);
        let selections = vec![Some(MaterialPropertySelection::SingleClaimOnly); 6];
        let volume = 18.0 * V;
        let viscosity = 1e6;
        for exponent in [1.0, 3.0] {
            let slope = -8e4;
            let fixture = fixture_with_rate(slope, viscosity, exponent);
            let mut final_states = Vec::new();
            // Different thermal substeps and feedback iteration counts must
            // converge to the same physical interval. They do not age history.
            for refined in [false, true] {
                let mut thermal = initial_member(&fixture);
                let plain = bind(&thermal, &[&fixture; 6], &nodes, &tetrahedra, cx);
                let mut solid = if refined {
                    plain
                        .with_power_law_overstress(&selections, cx)
                        .unwrap()
                        .with_voce_hardening(&selections, cx)
                        .unwrap()
                } else if exponent == 1.0 {
                    plain
                        .with_voce_hardening(&selections, cx)
                        .unwrap()
                        .with_linear_overstress(&selections, cx)
                        .unwrap()
                } else {
                    plain
                        .with_voce_hardening(&selections, cx)
                        .unwrap()
                        .with_power_law_overstress(&selections, cx)
                        .unwrap()
                };
                let mut old = Reference {
                    temperature: 300.0,
                    strain: 0.0,
                    alpha: 0.0,
                };
                let mut previous_force = f64::INFINITY;
                for step in 1..=3 {
                    let mut thermal_settings = thermal_config();
                    let mut coupling = config();
                    if refined {
                        thermal_settings.max_dt_s *= 0.5;
                        coupling.enthalpy_relaxation = 0.4;
                    }
                    let accepted = solid.solution().clone();
                    let next = solid
                        .advance(
                            &thermal,
                            &(),
                            step as f64 * DURATION,
                            &prescribed,
                            &[],
                            IncrementSettings::default(),
                            thermal_settings,
                            coupling,
                            cx,
                        )
                        .unwrap();
                    let expected = viscous_endpoint(old, slope, viscosity, exponent);
                    let da = expected.alpha - old.alpha;
                    let overstress = viscous_overstress(da, viscosity, exponent);
                    let force = tension(next.solid(), &prescribed);
                    close(force, L * L * expected.stress(), 2e-6);
                    assert!(force < previous_force);
                    assert!(overstress > 0.0 && da > 0.0);
                    let solution = next.solid().solution();
                    let lateral = -0.2 * (STRAIN - expected.alpha) - 0.5 * expected.alpha;
                    close(solution.displacement_m[7][1], L * lateral, 3e-11);
                    for element in 0..6 {
                        close(
                            next.thermal()
                                .temperature_k(CellId(100 + element as u64))
                                .unwrap(),
                            expected.temperature,
                            5e-8,
                        );
                        close(solution.element_states[element].alpha, expected.alpha, 3e-9);
                        close(solution.element_stress[element][0], expected.stress(), 3e-3);
                    }
                    let report = next.report();
                    close(
                        report.physical_heat_j,
                        volume * (YIELD + overstress) * da,
                        1e-8,
                    );
                    close(
                        report.coefficient_storage_change_j,
                        volume * coefficient(old, expected, slope),
                        1e-8,
                    );
                    close(
                        report.thermoelastic_storage_change_j,
                        volume * entropy_storage(old, expected, slope),
                        1e-8,
                    );
                    close(
                        report.reversible_heat_j,
                        volume * reversible(old, expected, slope),
                        1e-8,
                    );
                    close(
                        report.stored_mechanical_change_j,
                        volume * (expected.free_energy(slope) - old.free_energy(slope)),
                        1e-8,
                    );
                    close(
                        report.mechanical_work_j,
                        volume * expected.stress() * (expected.strain - old.strain),
                        1e-8,
                    );
                    let de = expected.strain - expected.alpha - old.strain + old.alpha;
                    let voce_loss = expected.q(slope)
                        * (-(-B_VALUE * expected.alpha).exp_m1() * da - primitive(expected.alpha)
                            + primitive(old.alpha));
                    close(
                        report.numerical_loss_j,
                        volume * (0.5 * E * de * de + 0.5 * H * da * da + voce_loss),
                        1e-8,
                    );
                    close(
                        report.thermal_audit.stored_energy_change_j,
                        report.external_heat_j + report.physical_heat_j + report.reversible_heat_j,
                        1e-8,
                    );
                    close(report.combined_residual_j, 0.0, coupling.energy_tolerance_j);
                    assert!(report.physical_heat_j > 0.0 && report.reversible_heat_j < 0.0);
                    if step > 1 {
                        assert!(report.coefficient_storage_change_j < 0.0);
                        close(report.mechanical_work_j, 0.0, 1e-8);
                    }
                    assert_eq!(solid.solution(), &accepted);
                    let (next_thermal, next_solid, _, _) = next.into_parts();
                    thermal = next_thermal;
                    solid = next_solid;
                    old = expected;
                    previous_force = force;
                }
                final_states.push((
                    thermal.temperature_k(CellId(100)).unwrap(),
                    solid.solution().element_states[0].alpha,
                    previous_force,
                ));
            }
            close(final_states[0].0, final_states[1].0, 5e-8);
            close(final_states[0].1, final_states[1].1, 3e-9);
            close(final_states[0].2, final_states[1].2, 2e-6);
        }
    });
}

#[test]
fn heated_loaded_member_softens_moves_and_balances_voce_storage_without_invented_heat() {
    context(|cx| {
        let (nodes, tetrahedra) = member();
        let prescribed = grips(&nodes, STRAIN);
        let volume = 18.0 * V;
        let mut outcomes = Vec::new();
        for slope in [0.0, -8e4] {
            let fixture = fixture(slope);
            let mut thermal = initial_member(&fixture);
            let mut solid = bind(&thermal, &[&fixture; 6], &nodes, &tetrahedra, cx)
                .with_voce_hardening(
                    &vec![Some(MaterialPropertySelection::SingleClaimOnly); 6],
                    cx,
                )
                .unwrap();
            let mut old = Reference {
                temperature: 300.0,
                strain: 0.0,
                alpha: 0.0,
            };
            let mut first_force = 0.0;
            let mut first_transverse = 0.0;
            for step in 1..=3 {
                let original_solution = solid.solution().clone();
                let original_cells = thermal.cells().clone();
                let next = solid
                    .advance(
                        &thermal,
                        &(),
                        step as f64 * DURATION,
                        &prescribed,
                        &[],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let expected = endpoint(old, slope);
                let solution = next.solid().solution();
                close(
                    tension(next.solid(), &prescribed),
                    expected.stress() * L * L,
                    2e-6,
                );
                close(solution.displacement_m[1][0], 3.0 * L * STRAIN, 1e-15);
                let lateral = -0.2 * (STRAIN - expected.alpha) - 0.5 * expected.alpha;
                close(solution.displacement_m[7][1], L * lateral, 3e-11);
                close(solution.displacement_m[7][2], L * lateral, 3e-11);
                for element in 0..6 {
                    close(
                        next.thermal()
                            .temperature_k(CellId(100 + element as u64))
                            .unwrap(),
                        expected.temperature,
                        5e-8,
                    );
                    close(solution.element_states[element].alpha, expected.alpha, 3e-9);
                    close(solution.element_stress[element][0], expected.stress(), 3e-3);
                    for stress in &solution.element_stress[element][1..] {
                        close(*stress, 0.0, 3e-3);
                    }
                    let material = next.solid().materials()[element].material();
                    close(
                        material.property(Q).unwrap().value_si(),
                        expected.q(slope),
                        0.02,
                    );
                    assert_eq!(material.property(B).unwrap().value_si(), B_VALUE);
                    for name in [Q, B] {
                        let selected = material.property(name).unwrap();
                        fixture
                            .card
                            .claims()
                            .verify_receipt(&selected.answer().receipt)
                            .unwrap();
                        assert_eq!(
                            selected.answer().receipt.selected,
                            solid.materials()[element]
                                .material()
                                .property(name)
                                .unwrap()
                                .answer()
                                .receipt
                                .selected
                        );
                    }
                    assert_eq!(
                        next.thermal().cells()[&CellId(100 + element as u64)].mass_kg,
                        RHO * 3.0 * V
                    );
                }
                let report = next.report();
                close(
                    report.coefficient_storage_change_j,
                    volume * coefficient(old, expected, slope),
                    1e-8,
                );
                close(
                    report.thermoelastic_storage_change_j,
                    volume * entropy_storage(old, expected, slope),
                    1e-8,
                );
                close(
                    report.reversible_heat_j,
                    volume * reversible(old, expected, slope),
                    1e-8,
                );
                close(
                    report.stored_mechanical_change_j,
                    volume * (expected.free_energy(slope) - old.free_energy(slope)),
                    1e-8,
                );
                close(
                    report.physical_heat_j,
                    volume * YIELD * (expected.alpha - old.alpha),
                    1e-8,
                );
                close(
                    report.mechanical_work_j,
                    volume * expected.stress() * (expected.strain - old.strain),
                    1e-8,
                );
                close(report.external_heat_j, volume * HEAT_PER_VOLUME, 1e-10);
                close(
                    report.thermal_audit.stored_energy_change_j,
                    volume * RHO * CP_VALUE * (expected.temperature - old.temperature),
                    1e-8,
                );
                close(
                    report.thermal_audit.stored_energy_change_j,
                    report.external_heat_j + report.physical_heat_j + report.reversible_heat_j,
                    1e-8,
                );
                let da = expected.alpha - old.alpha;
                let de = (expected.strain - expected.alpha) - (old.strain - old.alpha);
                let voce_loss = expected.q(slope)
                    * (-(-B_VALUE * expected.alpha).exp_m1() * da - primitive(expected.alpha)
                        + primitive(old.alpha));
                close(
                    report.numerical_loss_j,
                    volume * (0.5 * E * de * de + 0.5 * H * da * da + voce_loss),
                    1e-8,
                );
                close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                if step == 1 {
                    first_force = tension(next.solid(), &prescribed);
                    first_transverse = solution.displacement_m[7][1];
                    assert!(
                        next.solid()
                            .clone()
                            .with_voce_hardening(
                                &vec![Some(MaterialPropertySelection::SingleClaimOnly); 6],
                                cx,
                            )
                            .is_err()
                    );
                } else if slope < 0.0 {
                    assert!(tension(next.solid(), &prescribed) < first_force);
                    assert!(solution.displacement_m[7][1] < first_transverse);
                    assert!(solution.element_states[0].alpha > old.alpha);
                    assert!(report.coefficient_storage_change_j < 0.0);
                    assert!(report.physical_heat_j > 0.0);
                }
                assert_eq!(solid.solution(), &original_solution);
                assert_eq!(thermal.cells(), &original_cells);
                let (next_thermal, next_solid, _, _) = next.into_parts();
                thermal = next_thermal;
                solid = next_solid;
                old = expected;
            }
            outcomes.push((
                tension(&solid, &prescribed),
                solid.solution().displacement_m[7][1],
                old.alpha,
            ));
        }
        assert!(outcomes[1].0 < outcomes[0].0 - 1.0);
        assert!(outcomes[1].1 < outcomes[0].1 - 1e-7);
        assert!(outcomes[1].2 > outcomes[0].2 + 1e-5);
    });
}

#[test]
fn per_element_voce_selection_keeps_other_laws_and_combines_overstress_in_either_order() {
    context(|cx| {
        let (nodes, tetrahedra) = member();
        let fixture = fixture(-8e4);
        let thermal = initial_member(&fixture);
        let original = bind(&thermal, &[&fixture; 6], &nodes, &tetrahedra, cx);
        let mut selections = vec![None; 6];
        selections[2] = Some(MaterialPropertySelection::SingleClaimOnly);
        let selected = original
            .clone()
            .with_voce_hardening(&selections, cx)
            .unwrap();
        for (element, material) in selected.materials().iter().enumerate() {
            assert_eq!(material.material().property(Q).is_some(), element == 2);
            assert_eq!(material.material().property(B).is_some(), element == 2);
        }
        assert!(original.clone().with_voce_hardening(&[], cx).is_err());
        let voce_then_rate = selected.with_linear_overstress(&selections, cx).unwrap();
        let rate_then_voce = original
            .clone()
            .with_linear_overstress(&selections, cx)
            .unwrap()
            .with_voce_hardening(&selections, cx)
            .unwrap();
        assert_eq!(voce_then_rate.materials(), rate_then_voce.materials());
        for (element, material) in voce_then_rate.materials().iter().enumerate() {
            assert_eq!(
                material
                    .material()
                    .property(OVERSTRESS_VISCOSITY_PROPERTY)
                    .is_some(),
                element == 2
            );
        }
        assert!(
            original
                .materials()
                .iter()
                .all(|m| m.material().property(Q).is_none())
        );
        let missing = Fixture::new(CP_VALUE);
        let thermal = initial(&[&missing], &[], &[0.0]);
        let plain = bind(&thermal, &[&missing], &NODES, &TETS, cx);
        assert!(
            plain
                .with_voce_hardening(&[Some(MaterialPropertySelection::SingleClaimOnly)], cx)
                .is_err()
        );
    });
}
