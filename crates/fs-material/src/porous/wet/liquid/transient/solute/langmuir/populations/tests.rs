//! Original-source fixtures and independent multiple-population balances.
//! Cards below are authored arithmetic fixtures, not experimental validation.
use super::super::super::{dispersion, transport::DiluteSoluteBoundary};
use super::*;
use crate::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, CAPACITY_PARAMETER, HALF_SATURATION_PARAMETER,
    flow::{LOCAL_LDF_LAW, RATE_PARAMETER},
    kinetics::{
        AdsorptionBatch, AdsorptionNumerics, ELAPSED_TIME_AXIS, LDF_LAW, LangmuirKinetics,
        SorbentPopulation,
    },
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{
    Amount, Dims, Frequency, Length, Mass, Qty, Time,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn card(
    law: &str,
    params: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "porous-multiple-population-fixture",
        b"authored arithmetic source",
    );
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: params
            .iter()
            .map(|&(key, value, dims)| (key.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "Authored independent-site arithmetic, not measured material behavior".into(),
            license: "project fixture".into(),
            artifact: Some(artifact),
        },
    }
}

fn equilibrium(capacity: f64, half: f64, support: [f64; 2], element: &str) -> AqueousLangmuir {
    AqueousLangmuir::from_card(
        ElementId::new(element).unwrap(),
        card(
            AQUEOUS_LANGMUIR_LAW,
            &[
                (CAPACITY_PARAMETER, capacity, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
            ],
            ValidityDomain::unconstrained()
                .with_quantity(
                    CONCENTRATION_AXIS,
                    QuantitySpec::dimensional(AmountConcentration::DIMS),
                    support[0],
                    support[1],
                )
                .with_quantity(
                    TEMPERATURE_AXIS,
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    290.0,
                    310.0,
                )
                .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0, 12.0),
        ),
    )
    .unwrap()
}

fn kinetic(capacity: f64, half: f64, rate: f64) -> LocalLangmuirTransfer {
    kinetic_from(equilibrium(capacity, half, [0.0, 1000.0], "Au"), rate)
}
fn kinetic_from(eq: AqueousLangmuir, rate: f64) -> LocalLangmuirTransfer {
    let validity = eq
        .source()
        .validity
        .clone()
        .with_quantity(
            LOADING_AXIS,
            QuantitySpec::dimensional(AdsorbentLoading::DIMS),
            0.0,
            eq.capacity().value(),
        )
        .with_quantity(
            VOLUME_AXIS,
            QuantitySpec::dimensional(Volume::DIMS),
            1e-6,
            1e6,
        );
    LocalLangmuirTransfer::from_card(
        eq,
        card(
            LOCAL_LDF_LAW,
            &[(RATE_PARAMETER, rate, Frequency::DIMS)],
            validity,
        ),
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)
        .unwrap()
}
fn population<'a>(
    id: &'a str,
    cell: usize,
    law: LangmuirSoluteLaw<'a>,
    mass: f64,
    bound: f64,
    point: &'a QueryPoint,
) -> LangmuirSolutePopulation<'a> {
    LangmuirSolutePopulation {
        id,
        cell,
        species: 0,
        sorption: LangmuirSoluteSorption {
            law,
            dry_mass_kg_m2: mass,
            temperature: Temperature::new(298.15),
            conditions: point,
        },
        sorbed_amount_mol_m2: bound,
        vacant_amount_mol_m2: mass * law.equilibrium().capacity().value() - bound,
    }
}
fn controls<const N: usize>(maximum_populations: usize) -> MultiLangmuirSoluteControls<N> {
    MultiLangmuirSoluteControls {
        maximum_populations,
        langmuir: LangmuirSoluteControls {
            transport: DiluteSoluteControls {
                absolute_equation_tolerance_mol_m2: [1e-11; N],
                relative_equation_tolerance: 1e-12,
                amount_balance_tolerance_mol_m2: [1e-11; N],
                water_balance_tolerance_m: 1e-13,
                maximum_cell_species: 131072,
                maximum_storage_bytes: 64 * 1024 * 1024,
            },
            maximum_iterations: 128,
            maximum_line_searches: 32,
            maximum_linear_solves: 4096,
            maximum_source_evaluations: 1_000_000,
        },
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

#[allow(clippy::too_many_arguments)]
fn solve<'a, const N: usize>(
    dt: f64,
    water: &[f64],
    widths: &[f64],
    faces: &[f64],
    mobile: &[[f64; N]],
    populations: &'a [LangmuirSolutePopulation<'a>],
    area: f64,
    d: &[[f64; N]],
    tau: &[f64],
    boundaries: [DiluteSoluteBoundary<N>; 2],
    mechanical: Option<&[[f64; 2]]>,
    ctl: MultiLangmuirSoluteControls<N>,
    mut cancelled: impl FnMut() -> bool,
) -> core::result::Result<MultiLangmuirSoluteStep<N>, MultiLangmuirSoluteError> {
    let mut storage =
        MultiLangmuirStorage::new(mobile.len(), dt, populations, area, ctl, &mut cancelled)?;
    let result = transport::transport_with_dispersion(
        [0.0, dt],
        dt,
        water,
        water,
        widths,
        faces,
        mobile,
        d,
        tau,
        boundaries,
        ctl.langmuir.transport,
        &mut storage,
        mechanical,
        &mut cancelled,
    );
    let step = storage.recover(result)?;
    storage.finish(mobile, step, &mut cancelled)
}
fn closed<'a>(
    dt: f64,
    mobile: f64,
    populations: &'a [LangmuirSolutePopulation<'a>],
) -> MultiLangmuirSoluteStep<1> {
    solve(
        dt,
        &[0.4],
        &[1.0],
        &[0.0; 2],
        &[[mobile]],
        populations,
        1.0,
        &[[0.0]],
        &[1.0],
        [DiluteSoluteBoundary::Outflow; 2],
        None,
        controls(128),
        || false,
    )
    .unwrap()
}

// Independent positive root of A*c + Q*c/(H+c) = total.
fn root(a: f64, q: f64, h: f64, total: f64) -> f64 {
    let b = a * h + q - total;
    let discriminant = (b * b + 4.0 * a * h * total).sqrt();
    if b >= 0.0 {
        2.0 * h * total / (discriminant + b)
    } else {
        (discriminant - b) / (2.0 * a)
    }
}

#[test]
fn g1_two_rate_closed_cell_matches_independent_eliminated_quadratic() {
    let fast = kinetic(1.0, 0.5, 4.0);
    let slow = kinetic(2.0, 0.5, 0.03);
    let context = point();
    let populations = [
        population(
            "fast",
            0,
            LangmuirSoluteLaw::Kinetic(&fast),
            1.0,
            0.0,
            &context,
        ),
        population(
            "slow",
            0,
            LangmuirSoluteLaw::Kinetic(&slow),
            1.0,
            1.6,
            &context,
        ),
    ];
    for dt in [1e-6, 0.5, 1e8] {
        let answer = closed(dt, 0.4, &populations);
        let weights = [4.0 * dt / (1.0 + 4.0 * dt), 0.03 * dt / (1.0 + 0.03 * dt)];
        let c = root(
            0.4,
            weights[0] + 2.0 * weights[1],
            0.5,
            0.4 + 1.6 * weights[1],
        );
        near(answer.sorption.concentration_mol_m3[0][0], c, 1e-11);
        for (p, (capacity, old)) in [(1.0, 0.0), (2.0, 1.6)].into_iter().enumerate() {
            let expected = (1.0 - weights[p]) * old + weights[p] * capacity * c / (0.5 + c);
            near(answer.populations[p].sorbed_amount_mol_m2, expected, 1e-11);
            near(
                answer.populations[p].vacant_amount_mol_m2,
                capacity - expected,
                1e-11,
            );
        }
        near(
            answer.sorption.report.total.amount_before_mol_m2[0],
            2.0,
            1e-14,
        );
        near(
            answer.sorption.report.total.amount_after_mol_m2[0],
            2.0,
            1e-11,
        );
        near(
            answer.sorption.report.total.amount_change_mol_m2[0],
            0.0,
            1e-11,
        );
        assert!(answer.populations[0].exchange_into_sorbed_mol_m2 > 0.0);
        assert!(answer.populations[1].exchange_into_sorbed_mol_m2 < 0.0);
        if dt > 1e7 {
            near(c, root(0.4, 3.0, 0.5, 2.0), 3e-7);
        }
    }
}

#[test]
fn g1_equal_total_loading_retains_distinct_fast_and_slow_desorption_memory() {
    let fast = kinetic(1.0, 0.5, 4.0);
    let slow = kinetic(1.0, 0.5, 0.02);
    let context = point();
    let a = [
        population(
            "fast",
            0,
            LangmuirSoluteLaw::Kinetic(&fast),
            1.0,
            0.8,
            &context,
        ),
        population(
            "slow",
            0,
            LangmuirSoluteLaw::Kinetic(&slow),
            1.0,
            0.0,
            &context,
        ),
    ];
    let mut b = a;
    b[0].sorbed_amount_mol_m2 = 0.0;
    b[0].vacant_amount_mol_m2 = 1.0;
    b[1].sorbed_amount_mol_m2 = 0.8;
    b[1].vacant_amount_mol_m2 = 0.2;
    let first = closed(1.0, 0.0, &a);
    let second = closed(1.0, 0.0, &b);
    assert!(
        first.sorption.aqueous_amount_mol_m2[0][0]
            > 10.0 * second.sorption.aqueous_amount_mol_m2[0][0]
    );
    for step in [&first, &second] {
        near(
            step.sorption.report.total.amount_before_mol_m2[0],
            0.8,
            1e-14,
        );
        near(
            step.sorption.report.total.amount_after_mol_m2[0],
            0.8,
            1e-11,
        );
        assert!(
            step.populations
                .iter()
                .all(|p| p.sorbed_amount_mol_m2 >= 0.0 && p.vacant_amount_mol_m2 >= 0.0)
        );
    }
}

#[test]
fn g3_single_population_matches_original_and_frozen_populations_keep_stocks() {
    let model = kinetic(2.0, 0.5, 0.4);
    let frozen = kinetic(2.0, 0.5, 0.0);
    let slow = kinetic(2.0, 0.5, 1e-12);
    let context = point();
    let populations = [population(
        "one",
        0,
        LangmuirSoluteLaw::Kinetic(&model),
        0.5,
        0.2,
        &context,
    )];
    let answer = closed(0.7, 0.3, &populations);
    let bound = [[0.2]];
    let vacant = [[0.8]];
    let laws = [[Some(populations[0].sorption)]];
    let ctl = controls::<1>(128).langmuir;
    let mut old =
        LangmuirStorage::new(1, 0.7, &bound, &vacant, &laws, 1.0, ctl, &mut || false).unwrap();
    let result = transport::transport_with_storage(
        [0.0, 0.7],
        0.7,
        &[0.4],
        &[0.4],
        &[1.0],
        &[0.0; 2],
        &[[0.3]],
        &[[0.0]],
        &[1.0],
        [DiluteSoluteBoundary::Outflow; 2],
        ctl.transport,
        &mut old,
        || false,
    );
    let result = old.recover(result).unwrap();
    let original = old.finish(&[[0.3]], result, &mut || false).unwrap();
    assert_eq!(
        answer.sorption.aqueous_amount_mol_m2,
        original.sorption.aqueous_amount_mol_m2
    );
    assert_eq!(
        answer.sorption.sorbed_amount_mol_m2,
        original.sorption.sorbed_amount_mol_m2
    );
    assert_eq!(
        answer.sorption.exchange_into_sorbed_mol_m2,
        original.sorption.exchange_into_sorbed_mol_m2
    );
    assert_eq!(answer.source_evaluations, original.source_evaluations);
    for law in [&frozen, &slow] {
        let sites = [
            population("a", 0, LangmuirSoluteLaw::Kinetic(law), 0.5, 0.2, &context),
            population("b", 0, LangmuirSoluteLaw::Kinetic(law), 0.3, 0.4, &context),
        ];
        let result = closed(0.7, 0.3, &sites);
        near(result.sorption.aqueous_amount_mol_m2[0][0], 0.3, 1e-11);
        for (old, new) in sites.iter().zip(&result.populations) {
            if law.rate().value() == 0.0 {
                assert_eq!(
                    old.sorbed_amount_mol_m2.to_bits(),
                    new.sorbed_amount_mol_m2.to_bits()
                );
                assert_eq!(
                    old.vacant_amount_mol_m2.to_bits(),
                    new.vacant_amount_mol_m2.to_bits()
                );
            }
            near(old.sorbed_amount_mol_m2, new.sorbed_amount_mol_m2, 1e-11);
        }
    }
}

fn dispersion_model(alpha: f64) -> dispersion::LongitudinalDispersivity {
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            dispersion::PORE_SPEED_AXIS,
            QuantitySpec::dimensional(Qty::<1, 0, -1, 0, 0, 0>::DIMS),
            0.0,
            10.0,
        )
        .with_quantity(
            dispersion::LIQUID_SATURATION_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            1.0,
        )
        .with_quantity(
            dispersion::POROSITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.1,
            1.0,
        )
        .with_quantity(
            dispersion::TRANSPORT_SCALE_AXIS,
            QuantitySpec::dimensional(Length::DIMS),
            0.1,
            10.0,
        );
    dispersion::LongitudinalDispersivity::from_card(card(
        dispersion::LONGITUDINAL_DISPERSION_LAW,
        &[(dispersion::DISPERSIVITY_PARAMETER, alpha, Length::DIMS)],
        validity,
    ))
    .unwrap()
}

#[test]
fn g1_sparse_heterogeneous_sites_share_actual_advection_dispersion_and_face_oracle() {
    let fast = kinetic(1.0, 0.5, 3.0);
    let slow = kinetic(1.0, 0.5, 0.07);
    let eq = equilibrium(1.0, 0.5, [0.0, 1000.0], "Au");
    let context = point();
    let populations = [
        population(
            "slow",
            1,
            LangmuirSoluteLaw::Kinetic(&slow),
            0.3,
            0.2,
            &context,
        ),
        population(
            "fast",
            0,
            LangmuirSoluteLaw::Kinetic(&fast),
            0.3,
            0.1,
            &context,
        ),
        population(
            "instant",
            1,
            LangmuirSoluteLaw::Equilibrium(&eq),
            0.4,
            0.4 * 0.2 / 0.7,
            &context,
        ),
        population(
            "slow",
            0,
            LangmuirSoluteLaw::Kinetic(&slow),
            0.4,
            0.2,
            &context,
        ),
    ];
    let models = [dispersion_model(0.03), dispersion_model(0.07)];
    let extra = QueryPoint::new();
    let cells = core::array::from_fn::<_, 2, _>(|i| dispersion::SoluteDispersionCell {
        model: &models[i],
        transport_scale: Length::new(1.25),
        conditions: &extra,
    });
    let prepared = dispersion::prepare(
        [0.0, 0.3],
        0.3,
        &cells,
        dispersion::SoluteDispersionControls {
            maximum_source_points: 4,
            maximum_storage_bytes: 1_000_000,
        },
        |_| dispersion::DispersionCellState {
            porosity: 0.5,
            liquid_saturation: 0.8,
            liquid_volume_fraction: 0.4,
        },
        |_| Ok(0.03),
        &mut || false,
    )
    .unwrap();
    let original_ctl = controls::<1>(4);
    let mut ctl = original_ctl;
    ctl.langmuir.transport = prepared.transport_controls(ctl.langmuir.transport).unwrap();
    let boundaries = [
        DiluteSoluteBoundary::Reservoir {
            concentration_mol_m3: [0.8],
            film_resistance_s_m: [0.0],
        },
        DiluteSoluteBoundary::Outflow,
    ];
    let run = |sites: &[LangmuirSolutePopulation<'_>], ctl: MultiLangmuirSoluteControls<1>| {
        solve(
            0.3,
            &[0.2, 0.3],
            &[0.5, 0.75],
            &[0.03; 3],
            &[[0.3], [0.06]],
            sites,
            1.0,
            &[[0.02], [0.05]],
            &[2.0, 3.0],
            boundaries,
            Some(&prepared.report.mechanical_bulk_m2_s),
            ctl,
            || false,
        )
    };
    let answer = run(&populations, ctl).unwrap();
    let mechanical = [0.003, 0.007];
    let bulk = [
        0.4 * 0.02 / 2.0 + mechanical[0],
        0.4 * 0.05 / 3.0 + mechanical[1],
    ];
    let g = 0.3 / (0.25 / bulk[0] + 0.375 / bulk[1]);
    let inlet = 0.3 * bulk[0] / 0.25;
    let weight_fast = 0.9 / 1.9;
    let weight_slow = 0.021 / 1.021;
    let a = [0.2 + 0.03 + g + inlet, 0.3 + 0.03 + g];
    let q = [
        0.3 * weight_fast + 0.4 * weight_slow,
        0.4 + 0.3 * weight_slow,
    ];
    let rhs = [
        0.3 + (0.03 + inlet) * 0.8 + weight_fast * 0.1 + weight_slow * 0.2,
        0.06 + 0.4 * 0.2 / 0.7 + weight_slow * 0.2,
    ];
    let second = |first: f64| root(a[1], q[1], 0.5, rhs[1] + (0.03 + g) * first);
    let (mut lo, mut hi) = (0.0, 10.0);
    for _ in 0..100 {
        let c = 0.5 * (lo + hi);
        let f = a[0] * c + q[0] * c / (0.5 + c) - g * second(c) - rhs[0];
        if f > 0.0 {
            hi = c;
        } else {
            lo = c;
        }
    }
    let c = [0.5 * (lo + hi), second(0.5 * (lo + hi))];
    for i in 0..2 {
        near(answer.sorption.concentration_mol_m3[i][0], c[i], 1e-11);
    }
    near(
        answer.sorption.report.total.diffusive_face_amount_mol_m2[1][0],
        g * (c[0] - c[1]),
        1e-12,
    );
    near(
        answer.sorption.report.total.advective_face_amount_mol_m2[1][0],
        0.03 * c[0],
        1e-12,
    );
    near(
        answer.sorption.report.total.amount_balance_defect_mol_m2[0],
        0.0,
        1e-11,
    );
    for (p, state) in populations.iter().zip(&answer.populations) {
        let rate = p.sorption.law.rate();
        let weight = rate.map_or(1.0, |k| k * 0.3 / (1.0 + k * 0.3));
        let target = p.sorption.dry_mass_kg_m2 * c[p.cell] / (0.5 + c[p.cell]);
        near(
            state.sorbed_amount_mol_m2,
            (1.0 - weight) * p.sorbed_amount_mol_m2 + weight * target,
            1e-11,
        );
    }
    let reversed = [
        populations[3],
        populations[2],
        populations[1],
        populations[0],
    ];
    let permutation = run(&reversed, ctl).unwrap();
    for i in 0..2 {
        near(permutation.sorption.concentration_mol_m3[i][0], c[i], 1e-11);
    }
    for i in 0..4 {
        near(
            permutation.populations[3 - i].sorbed_amount_mol_m2,
            answer.populations[i].sorbed_amount_mol_m2,
            1e-13,
        );
    }
    let required =
        answer.sorption.report.total.workspace_bytes + prepared.report.preparation_storage_bytes;
    let mut exact = original_ctl;
    exact.langmuir.transport.maximum_storage_bytes = required;
    exact.langmuir.transport = prepared
        .transport_controls(exact.langmuir.transport)
        .unwrap();
    assert!(run(&populations, exact).is_ok());
    exact.langmuir.transport.maximum_storage_bytes -= 1;
    let error = run(&populations, exact).unwrap_err();
    assert_eq!(
        prepared.multi_langmuir_error(error),
        dispersion::SoluteDispersionError::MultiLangmuir(MultiLangmuirSoluteError::Transport(
            DiluteSoluteError::Workspace {
                required_bytes: required,
                maximum_bytes: required - 1
            }
        ))
    );
}

#[test]
fn g3_no_spatial_flux_matches_original_multi_population_batch_at_same_be_steps() {
    let models = [kinetic(2.0, 0.5, 3.0), kinetic(2.0, 0.5, 0.3)];
    let context = point();
    let area = 2.5;
    let mut populations = [
        population(
            "fast",
            0,
            LangmuirSoluteLaw::Kinetic(&models[0]),
            0.6,
            0.1,
            &context,
        ),
        population(
            "slow",
            0,
            LangmuirSoluteLaw::Kinetic(&models[1]),
            0.8,
            0.6,
            &context,
        ),
    ];
    let batch_populations = populations
        .iter()
        .zip(&models)
        .map(|(p, model)| {
            let source = card(
                LDF_LAW,
                &[(RATE_PARAMETER, model.rate().value(), Frequency::DIMS)],
                model.source().validity.clone().with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    0.0,
                    1000.0,
                ),
            );
            let kinetic = LangmuirKinetics::from_card(model.equilibrium().clone(), source).unwrap();
            SorbentPopulation::new(
                p.id,
                kinetic,
                Mass::new(area * p.sorption.dry_mass_kg_m2),
                Amount::new(area * p.sorbed_amount_mol_m2),
                context.clone(),
            )
            .unwrap()
        })
        .collect();
    let mut batch = AdsorptionBatch::new(
        batch_populations,
        Amount::new(area * 0.16),
        Volume::new(area * 0.4),
        Temperature::new(298.15),
    )
    .unwrap();
    let batch_step = batch
        .advance(
            Time::new(0.4),
            AdsorptionNumerics {
                absolute_amount_error: Amount::new(10.0),
                maximum_step: Time::new(0.4),
                ..AdsorptionNumerics::default()
            },
            || false,
        )
        .unwrap();
    assert_eq!(batch_step.accepted_steps, 1);
    assert_eq!(batch_step.rejected_steps, 0);
    let mut mobile = 0.16;
    for _ in 0..2 {
        let result = solve(
            0.2,
            &[0.4],
            &[1.0],
            &[0.0; 2],
            &[[mobile]],
            &populations,
            area,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            None,
            controls(2),
            || false,
        )
        .unwrap();
        mobile = result.sorption.aqueous_amount_mol_m2[0][0];
        for (input, state) in populations.iter_mut().zip(&result.populations) {
            input.sorbed_amount_mol_m2 = state.sorbed_amount_mol_m2;
            input.vacant_amount_mol_m2 = state.vacant_amount_mol_m2;
        }
    }
    near(area * mobile, batch.dissolved().value(), 1e-10);
    for (p, state) in populations.iter().zip(batch.populations()) {
        assert_eq!(p.id, state.id);
        near(area * p.sorbed_amount_mol_m2, state.adsorbed.value(), 1e-10);
        near(
            area * p.vacant_amount_mol_m2,
            state.vacant_capacity.value(),
            1e-10,
        );
    }
}

#[test]
fn g0_each_population_keeps_own_source_support_identity_and_initial_equilibrium() {
    let first = kinetic_from(equilibrium(1.0, 0.5, [0.0, 0.2], "Au"), 0.4);
    let second = kinetic_from(equilibrium(1.0, 0.5, [0.3, 1.0], "Au"), 0.2);
    let context = point();
    let sites = [
        population(
            "first",
            0,
            LangmuirSoluteLaw::Kinetic(&first),
            1.0,
            0.0,
            &context,
        ),
        population(
            "second",
            0,
            LangmuirSoluteLaw::Kinetic(&second),
            1.0,
            0.0,
            &context,
        ),
    ];
    assert!(matches!(
        MultiLangmuirStorage::new(1, 1.0, &sites, 1.0, controls::<1>(2), &mut || false),
        Err(MultiLangmuirSoluteError::Population { population: 1, .. })
    ));
    let broad = kinetic(1.0, 0.5, 0.5);
    let eq = equilibrium(1.0, 0.5, [0.0, 1000.0], "Au");
    let good = population(
        "valid",
        0,
        LangmuirSoluteLaw::Kinetic(&broad),
        1.0,
        0.0,
        &context,
    );
    let wrong = population(
        "instant",
        0,
        LangmuirSoluteLaw::Equilibrium(&eq),
        1.0,
        0.0,
        &context,
    );
    let trial = |items: &[LangmuirSolutePopulation<'_>]| {
        solve(
            1.0,
            &[0.4],
            &[1.0],
            &[0.0; 2],
            &[[0.4]],
            items,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            None,
            controls(128),
            || false,
        )
    };
    assert!(matches!(
        trial(&[good, wrong]),
        Err(MultiLangmuirSoluteError::Population {
            population: 1,
            error: LangmuirSoluteError::Transport(DiluteSoluteError::Sorption { .. })
        })
    ));
    assert!(matches!(
        trial(&[good, good]),
        Err(MultiLangmuirSoluteError::Population { population: 1, .. })
    ));
    let bad_context = QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0)
        .unwrap();
    let bad = population(
        "bad-pH",
        0,
        LangmuirSoluteLaw::Kinetic(&broad),
        1.0,
        0.0,
        &bad_context,
    );
    assert!(matches!(
        trial(&[good, bad]),
        Err(MultiLangmuirSoluteError::Population {
            population: 1,
            error: LangmuirSoluteError::KineticSource {
                cell: 0,
                species: 0,
                ..
            }
        })
    ));
    let forged = context
        .clone()
        .with_quantity(VOLUME_AXIS, QuantitySpec::dimensional(Volume::DIMS), 0.4)
        .unwrap();
    let forged = population(
        "forged-volume",
        0,
        LangmuirSoluteLaw::Kinetic(&broad),
        1.0,
        0.0,
        &forged,
    );
    assert!(matches!(
        trial(&[good, forged]),
        Err(MultiLangmuirSoluteError::Population { population: 1, .. })
    ));
}

#[test]
fn g3_sparse_constituent_slots_preserve_passive_diffusion_and_element_binding() {
    let gold = kinetic(1.0, 0.5, 0.5);
    let silver = kinetic_from(equilibrium(1.0, 0.5, [0.0, 1000.0], "Ag"), 0.2);
    let context = point();
    let a = population(
        "gold",
        0,
        LangmuirSoluteLaw::Kinetic(&gold),
        0.5,
        0.0,
        &context,
    );
    let mut b = population(
        "silver",
        1,
        LangmuirSoluteLaw::Kinetic(&silver),
        0.5,
        0.0,
        &context,
    );
    b.species = 1;
    let sites = [b, a];
    let mobile = [[0.4, 0.0, 0.2], [0.0, 0.4, 0.1]];
    let boundaries = [DiluteSoluteBoundary::Outflow; 2];
    let ctl = controls::<3>(2);
    let answer = solve(
        0.3,
        &[0.4; 2],
        &[1.0; 2],
        &[0.0; 3],
        &mobile,
        &sites,
        1.0,
        &[[0.05; 3]; 2],
        &[1.0; 2],
        boundaries,
        None,
        ctl,
        || false,
    )
    .unwrap();
    let passive = transport::transport(
        [0.0, 0.3],
        0.3,
        &[0.4; 2],
        &[0.4; 2],
        &[1.0; 2],
        &[0.0; 3],
        &mobile,
        &[[0.05; 3]; 2],
        &[1.0; 2],
        boundaries,
        ctl.langmuir.transport,
        || false,
    )
    .unwrap();
    for i in 0..2 {
        assert_eq!(
            answer.sorption.aqueous_amount_mol_m2[i][2].to_bits(),
            passive.amount_mol_m2[i][2].to_bits()
        );
        assert_eq!(answer.sorption.sorbed_amount_mol_m2[i][2], 0.0);
    }
    assert!(
        answer.populations[0].sorbed_amount_mol_m2 > 0.0
            && answer.populations[1].sorbed_amount_mol_m2 > 0.0
    );
    let mut wrong = sites;
    wrong[0].sorption.law = LangmuirSoluteLaw::Kinetic(&gold);
    assert!(matches!(
        MultiLangmuirStorage::new(2, 0.3, &wrong, 1.0, ctl, &mut || false),
        Err(MultiLangmuirSoluteError::Transport(
            DiluteSoluteError::Invalid(_)
        ))
    ));
    let mut wrong_temperature = sites;
    wrong_temperature[0].cell = 0;
    wrong_temperature[0].sorption.temperature = Temperature::new(299.15);
    assert!(matches!(
        MultiLangmuirStorage::new(2, 0.3, &wrong_temperature, 1.0, ctl, &mut || false),
        Err(MultiLangmuirSoluteError::Population { population: 0, .. })
    ));
}

#[test]
fn g4_shared_source_population_memory_and_cancellation_limits_are_atomic() {
    let model = kinetic(1.0, 0.5, 0.5);
    let context = point();
    let sites = [
        population(
            "a",
            0,
            LangmuirSoluteLaw::Kinetic(&model),
            1.0,
            0.0,
            &context,
        ),
        population(
            "b",
            0,
            LangmuirSoluteLaw::Kinetic(&model),
            0.5,
            0.1,
            &context,
        ),
    ];
    let run = |ctl: MultiLangmuirSoluteControls<1>| {
        solve(
            1.0,
            &[0.4],
            &[1.0],
            &[0.0; 2],
            &[[0.4]],
            &sites,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            None,
            ctl,
            || false,
        )
    };
    let mut ctl = controls(2);
    let answer = run(ctl).unwrap();
    ctl.langmuir.transport.maximum_storage_bytes = answer.sorption.report.total.workspace_bytes;
    ctl.langmuir.maximum_source_evaluations = answer.source_evaluations;
    assert!(run(ctl).is_ok());
    let mut insufficient = ctl;
    insufficient.langmuir.maximum_source_evaluations -= 1;
    assert!(matches!(
        run(insufficient),
        Err(MultiLangmuirSoluteError::Population {
            error: LangmuirSoluteError::SourceBudget { .. },
            ..
        })
    ));
    insufficient = ctl;
    insufficient.maximum_populations = 1;
    assert_eq!(
        run(insufficient).unwrap_err(),
        MultiLangmuirSoluteError::PopulationWork {
            required_populations: 2,
            maximum_populations: 1
        }
    );
    insufficient = ctl;
    insufficient.langmuir.transport.maximum_storage_bytes -= 1;
    assert_eq!(
        run(insufficient).unwrap_err(),
        MultiLangmuirSoluteError::Transport(DiluteSoluteError::Workspace {
            required_bytes: ctl.langmuir.transport.maximum_storage_bytes,
            maximum_bytes: ctl.langmuir.transport.maximum_storage_bytes - 1
        })
    );
    let mut polls = 0;
    let result = solve(
        1.0,
        &[0.4],
        &[1.0],
        &[0.0; 2],
        &[[0.4]],
        &sites,
        1.0,
        &[[0.0]],
        &[1.0],
        [DiluteSoluteBoundary::Outflow; 2],
        None,
        ctl,
        || {
            polls += 1;
            polls == 9
        },
    );
    assert!(matches!(
        result,
        Err(MultiLangmuirSoluteError::Transport(
            DiluteSoluteError::Cancelled
        ))
    ));
    assert_eq!(sites[0].sorbed_amount_mol_m2, 0.0);
    assert_eq!(sites[1].sorbed_amount_mol_m2, 0.1);
}
