//! G0/G1: complete cyclic memories and mixed original caloric inventories.
//! Coefficients are manufactured; no experimental alloy calibration is claimed.

use super::*;
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::thermal::chaboche::ThermalChabocheSource;
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_solid::fiber::FiberLaw;
use fs_solid::frame2d::creep::FrameThermalFiberLaw;

const H: f64 = 1e8;
const Q: f64 = 1e6;
const B: f64 = 200.0;
const C: [f64; 2] = [4e8, 1e8];
const DC: [f64; 2] = [-4e6, 1e6];
const GAMMA: [f64; 2] = [120.0, 15.0];

fn branches() -> [BackstressProperties; 2] {
    [
        BackstressProperties::new("fast_C", "fast_gamma"),
        BackstressProperties::new("slow_C", "slow_gamma"),
    ]
}

fn cyclic_fixture(revision: u32) -> Fixture {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, DE),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2, 0.0),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6, -10e3),
        (BASIS.property_name(), Pressure::DIMS, H, 0.0),
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK, DHK),
        (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, Q, 0.0),
        (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, B, 0.0),
        ("specific_heat_capacity", CP_DIMS, CP, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, RHO, 0.0),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            EXPANSION,
            0.0,
        ),
        (OVERSTRESS_VISCOSITY_PROPERTY, OVERSTRESS_VISCOSITY_DIMS, 1e7, -50e3),
        (OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 2e6, 0.0),
        (OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0, 0.0),
    ] {
        claims.insert_claim(claim(name, dims, value, slope)).unwrap();
    }
    for (i, names) in branches().iter().enumerate() {
        claims
            .insert_claim(claim(&names.modulus, Pressure::DIMS, C[i], DC[i]))
            .unwrap();
        claims
            .insert_claim(claim(&names.recovery, Dims::NONE, GAMMA[i], 0.0))
            .unwrap();
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "cyclic frame fixture".into(),
            phase: "solid".into(),
            process: "manufactured two-branch calibration".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    // A lineage starts at revision 0; a later revision is its successor.
    for _ in 0..revision {
        card =
            MaterialCard::supersede(&card, card.claims().clone(), card.models().to_vec()).unwrap();
    }
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(340.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1e-6,
    )
    .unwrap();
    Fixture { card, caloric }
}

fn cyclic_source(fixture: &Fixture, rate: bool) -> ThermalChabocheSource<'_, 2> {
    let base = ThermoplasticSource::try_new(
        SourceQuery {
            card: &fixture.card,
            caloric: &fixture.caloric,
            context: point(300.0),
            temperature_axis: "T".into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            hardening_basis: BASIS,
            thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
            strain_limit: 0.04,
        },
        0.0,
    )
    .unwrap()
    .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
    .unwrap();
    let base = if rate {
        base.with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
            .unwrap()
    } else {
        base
    };
    base.with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, branches())
        .unwrap()
}

fn cyclic_law(fixture: &Fixture, rate: bool) -> FiberLaw {
    FiberMaterial::thermal_chaboche(
        cyclic_source(fixture, rate)
            .at_uniaxial_temperatures(300.0, 300.0, 1.0)
            .unwrap(),
    )
    .unwrap()
    .law()
    .clone()
}

fn frame_with_laws(laws: [FiberLaw; 2]) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame
        .add_element(0, 1, 3, RHO * V, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .into_iter()
                    .enumerate()
                    .map(|(i, y)| Fiber {
                        y,
                        area: AREA,
                        law: laws[i].clone(),
                    })
                    .collect(),
            ))
        })
        .unwrap();
    frame
}

#[test]
fn heated_cyclic_frame_retains_each_branch_through_reversal_and_release() {
    let fixture = cyclic_fixture(0);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut thermal = fixture.thermal(1, 1.0);
        let law = cyclic_law(&fixture, false);
        let mut frame = ThermoplasticFrame::try_new(
            frame_with_laws([law.clone(), law]),
            &thermal,
            vec![FrameThermalCell {
                cell: CellId(0),
                source: cyclic_source(&fixture, false).into(),
            }],
            fixture.mappings(1),
            6,
            cx,
        )
        .unwrap();
        let mut observer = Heating::default();
        let mut old = [300.0, 0.0, 0.0, 0.0, 0.0];
        let mut target = 0.0;
        for (index, (dt, stress, external_j)) in [
            (0.11, 6e6, 0.0),
            (0.21, 6e6, RHO * V * CP * 15.0),
            (0.13, -2.5e6, 0.0),
            (0.17, 0.0, -RHO * V * CP * 10.0),
        ]
        .into_iter()
        .enumerate()
        {
            target += dt;
            observer = observer.with_heat(&[(CellId(0), external_j)], dt);
            let accepted_u = frame.frame().displacements().to_vec();
            let step = frame
                .advance(
                    &thermal,
                    &observer,
                    target,
                    &axial_load(stress),
                    FrameThermalMechanics::Quasistatic,
                    mechanical(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert_eq!(frame.frame().displacements(), accepted_u);
            let temperature = fixture
                .caloric
                .curve()
                .state_at_specific_enthalpy(
                    step.thermal().cells()[&CellId(0)].specific_h_j_kg,
                )
                .unwrap()
                .temperature_k();
            let points = step
                .frame()
                .frame()
                .thermal_inelastic_fibers(6, || false)
                .unwrap();
            let mut next = old;
            let mut recovery_j = 0.0;
            for point in points {
                let FrameThermalFiberLaw::Chaboche { state, .. } = point.law else {
                    panic!("cyclic history changed material family");
                };
                assert_eq!(state.branch_count(), 2);
                assert_eq!(state.temperature_k(), temperature);
                let p = state.plastic_strain();
                let alpha = state.equivalent_plastic_strain();
                let dp = p - old[1];
                let da = alpha - old[2];
                near(da, dp.abs(), 2e-10);
                let mut beta = [0.0; 2];
                for i in 0..2 {
                    let c0 = C[i] + DC[i] * (old[0] - 300.0);
                    let c1 = C[i] + DC[i] * (temperature - 300.0);
                    beta[i] = state.backstress_pa(i).unwrap();
                    // Independent axial Armstrong-Frederick recurrence and
                    // fixed strain-like memory under changing temperature.
                    let expected =
                        (old[3 + i] * c1 / c0 + c1 * dp) / (1.0 + GAMMA[i] * da);
                    near(beta[i], expected, 0.1);
                    recovery_j += point.reference_volume_m3
                        * GAMMA[i] * beta[i] * beta[i] / c1 * da;
                }
                near(
                    youngs(temperature)
                        * (point.axial_strain - p - EXPANSION * (temperature - 300.0)),
                    stress,
                    0.2,
                );
                let radius = 3e6 - 10e3 * (temperature - 300.0)
                    + H * alpha + Q * (1.0 - (-B * alpha).exp());
                let effective = stress - hk(temperature) * p - beta.iter().sum::<f64>();
                if da > 1e-9 {
                    near(effective.abs(), radius, 0.2);
                } else {
                    assert!(effective.abs() <= radius + 0.2);
                }
                if index == 2 {
                    assert!(dp < -1e-6, "reverse flow must retain prior directional memory");
                }
                if index == 3 {
                    near(alpha, old[2], 1e-9);
                }
                next = [temperature, p, alpha, beta[0], beta[1]];
            }
            let report = step.report();
            near(report.mechanical_energy.recovery_dissipation_j, recovery_j, 1e-7);
            if index < 3 {
                assert!(recovery_j > 0.0);
            }
            near(report.external_heat_j, external_j, 1e-7);
            near(
                report.cell_heat_j[&CellId(0)],
                RHO * V * CP * (temperature - old[0]) - external_j,
                1e-7,
            );
            near(report.thermal_residual_j, 0.0, 1e-7);
            near(report.combined_residual_j, 0.0, 1e-7);
            let (next_thermal, next_frame, next_observer, _) = step.into_parts();
            thermal = next_thermal;
            frame = next_frame;
            observer = next_observer;
            old = next;
        }
        assert!(old[2] > old[1].abs());
        assert!(old[1] > 1e-4);
        near(
            frame.frame().displacements()[3],
            old[1] + EXPANSION * (old[0] - 300.0),
            2e-10,
        );
    });
}

fn mixed_thermal<'a>(left: &'a Fixture, right: &'a Fixture) -> DepositionSimulation<'a> {
    let events = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i),
                segment_id: i,
                body_id: 1,
                material_id: MaterialId(i + 1),
                mass_kg: RHO * AREA,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: AREA,
                    surface_area_m2: 0.02,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(1), left.caloric.curve() as &dyn ThermalLaw),
            (MaterialId(2), right.caloric.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 4,
            max_events: 4,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, heat_config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn mixed_tabulated_and_rate_cyclic_fibers_use_their_original_cells_atomically() {
    let tabulated = Fixture::new(0);
    let cyclic = cyclic_fixture(1);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let thermal = mixed_thermal(&tabulated, &cyclic);
        let table = FiberMaterial::thermal_tabulated(
            tabulated.source(false)
                .at_uniaxial_temperatures(300.0, 300.0, 1.0)
                .unwrap(),
        )
        .unwrap();
        let raw = frame_with_laws([table.law().clone(), cyclic_law(&cyclic, true)]);
        let sources = vec![
            FrameThermalCell { cell: CellId(0), source: tabulated.source(false).into() },
            FrameThermalCell { cell: CellId(1), source: cyclic_source(&cyclic, true).into() },
        ];
        let mut swapped = sources.clone();
        swapped.swap(0, 1);
        swapped[0].cell = CellId(0);
        swapped[1].cell = CellId(1);
        assert!(ThermoplasticFrame::try_new(
            raw.clone(), &thermal, swapped, tabulated.mappings(2), 6, cx,
        )
        .is_err());
        let frame = ThermoplasticFrame::try_new(
            raw, &thermal, sources, tabulated.mappings(2), 6, cx,
        )
        .unwrap();
        let observer = Heating::default().with_heat(
            &[(CellId(0), RHO * AREA * CP * 5.0), (CellId(1), RHO * AREA * CP * 15.0)],
            0.17,
        );
        let before_u = frame.frame().displacements().to_vec();
        let mut failing = observer.clone();
        failing.fail = true;
        assert!(frame.advance(
            &thermal, &failing, 0.17, &axial_load(5e6),
            FrameThermalMechanics::Quasistatic,
            mechanical(), heat_config(), config(), cx,
        ).is_err());
        assert_eq!(frame.frame().displacements(), before_u);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(failing.delivered_j, 0.0);
        let step = frame.advance(
            &thermal, &observer, 0.17, &axial_load(5e6),
            FrameThermalMechanics::Quasistatic,
            mechanical(), heat_config(), config(), cx,
        ).unwrap();
        let mut counts = [0; 2];
        for point in step.frame().frame().thermal_inelastic_fibers(6, || false).unwrap() {
            let accepted_temperature = point.law.temperature_k();
            let (family, caloric) = match point.law {
                FrameThermalFiberLaw::Tabulated { state, .. } => {
                    assert!(state.mechanical.equivalent_plastic_strain() > 0.0);
                    (0, &tabulated.caloric)
                }
                FrameThermalFiberLaw::Chaboche { state, .. } => {
                    assert_eq!(state.branch_count(), 2);
                    assert!(state.equivalent_plastic_strain() > 0.0);
                    (1, &cyclic.caloric)
                }
            };
            counts[family] += 1;
            let actual_temperature = caloric.curve()
                .state_at_specific_enthalpy(
                    step.thermal().cells()[&CellId(family as u64)].specific_h_j_kg,
                ).unwrap().temperature_k();
            assert_eq!(accepted_temperature, actual_temperature);
        }
        assert_eq!(counts, [3, 3]);
        assert!(step.report().mechanical_energy.recovery_dissipation_j > 0.0);
        assert!(step.report().mechanical_energy.overstress_dissipation_j > 0.0);
        near(step.report().combined_residual_j, 0.0, 1e-7);
        near(step.report().thermal_residual_j, 0.0, 1e-7);
    });
}
