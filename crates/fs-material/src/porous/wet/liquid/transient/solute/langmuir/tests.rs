//! G1/G3 source-backed spatial Langmuir storage and original face-flux limits.
//! All cards are authored arithmetic fixtures, not measured adsorption data.
use super::super::transport::{DiluteSoluteBoundary, transport, transport_with_storage};
use super::*;
use crate::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, CAPACITY_PARAMETER, HALF_SATURATION_PARAMETER,
    flow::{LOCAL_LDF_LAW, RATE_PARAMETER},
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{
    Dims, Frequency,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn source(
    law: &str,
    params: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "porous-langmuir-test",
        b"authored finite-capacity transport fixture",
    );
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: params
            .iter()
            .map(|&(n, value, dims)| (n.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "Authored spatial adsorption arithmetic, not measured material validation"
                .into(),
            license: "project fixture".into(),
            artifact: Some(artifact),
        },
    }
}
fn equilibrium(capacity: f64, half: f64, upper: f64) -> AqueousLangmuir {
    AqueousLangmuir::from_card(
        ElementId::new("Au").unwrap(),
        source(
            AQUEOUS_LANGMUIR_LAW,
            &[
                (CAPACITY_PARAMETER, capacity, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
            ],
            ValidityDomain::unconstrained()
                .with_quantity(
                    CONCENTRATION_AXIS,
                    QuantitySpec::dimensional(AmountConcentration::DIMS),
                    0.0,
                    upper,
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
    let eq = equilibrium(capacity, half, 1000.0);
    let validity = eq
        .source()
        .validity
        .clone()
        .with_quantity(
            LOADING_AXIS,
            QuantitySpec::dimensional(AdsorbentLoading::DIMS),
            0.0,
            capacity,
        )
        .with_quantity(
            VOLUME_AXIS,
            QuantitySpec::dimensional(Volume::DIMS),
            1e-6,
            1e6,
        );
    LocalLangmuirTransfer::from_card(
        eq,
        source(
            LOCAL_LDF_LAW,
            &[(RATE_PARAMETER, rate, Frequency::DIMS)],
            validity,
        ),
    )
    .unwrap()
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)
        .unwrap()
}
fn binding<'a>(
    law: LangmuirSoluteLaw<'a>,
    dry_mass_kg_m2: f64,
    conditions: &'a QueryPoint,
) -> Option<LangmuirSoluteSorption<'a>> {
    Some(LangmuirSoluteSorption {
        law,
        dry_mass_kg_m2,
        conditions,
        temperature: Temperature::new(298.15),
    })
}
fn controls<const N: usize>() -> LangmuirSoluteControls<N> {
    LangmuirSoluteControls {
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
        maximum_source_evaluations: 100_000,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}

#[allow(clippy::too_many_arguments)]
fn solve<'a, const N: usize>(
    dt: f64,
    before: &[f64],
    after: &[f64],
    widths: &[f64],
    water: &[f64],
    mobile: &[[f64; N]],
    bound: &'a [[f64; N]],
    vacant: &'a [[f64; N]],
    laws: &'a [[Option<LangmuirSoluteSorption<'a>>; N]],
    area: f64,
    d: &[[f64; N]],
    tau: &[f64],
    boundaries: [DiluteSoluteBoundary<N>; 2],
    ctl: LangmuirSoluteControls<N>,
    mut cancelled: impl FnMut() -> bool,
) -> core::result::Result<LangmuirSoluteStep<N>, LangmuirSoluteError> {
    let mut storage = LangmuirStorage::new(
        mobile.len(),
        dt,
        bound,
        vacant,
        laws,
        area,
        ctl,
        &mut cancelled,
    )?;
    let result = transport_with_storage(
        [0.0, dt],
        dt,
        before,
        after,
        widths,
        water,
        mobile,
        d,
        tau,
        boundaries,
        ctl.transport,
        &mut storage,
        &mut cancelled,
    );
    let step = storage.recover(result)?;
    storage.finish(mobile, step, &mut cancelled)
}

// Independent positive quadratic root of V*c + Q*c/(H+c) = total.
fn concentration(v: f64, q: f64, h: f64, total: f64) -> f64 {
    let b = v * h + q - total;
    let disc = (b * b + 4.0 * v * total * h).sqrt();
    if b >= 0.0 {
        2.0 * total * h / (disc + b)
    } else {
        (disc - b) / (2.0 * v)
    }
}

#[test]
fn g1_closed_kinetic_stock_matches_independent_two_phase_quadratic() {
    let model = kinetic(2.0, 0.3, 0.7);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Kinetic(&model), 2.0, &point)]];
    for dt in [1e-4, 0.5, 1e6] {
        let result = solve(
            dt,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &[[0.2]],
            &[[0.1]],
            &[[3.9]],
            &laws,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false,
        )
        .unwrap();
        let weight = 0.7 * dt / (1.0 + 0.7 * dt);
        let c = concentration(0.25, weight * 4.0, 0.3, 0.2 + weight * 0.1);
        let b = (0.1 + 0.7 * dt * 4.0 * c / (0.3 + c)) / (1.0 + 0.7 * dt);
        near(result.sorption.concentration_mol_m3[0][0], c, 1e-11);
        near(result.sorption.sorbed_amount_mol_m2[0][0], b, 1e-11);
        near(result.vacant_amount_mol_m2[0][0], 4.0 - b, 1e-11);
        near(
            result.sorption.report.total.amount_change_mol_m2[0],
            0.0,
            1e-11,
        );
        assert!(result.sorption.report.total.linear_solves > 0);
    }
}

#[test]
fn g1_heterogeneous_equilibrium_diffusion_matches_independent_face_root() {
    let a = equilibrium(2.0, 0.5, 100.0);
    let b = equilibrium(2.0, 1.5, 100.0);
    let point = conditions();
    let mass = [0.4, 0.8];
    let v = [0.2, 0.4];
    let initial = [2.0, 0.25];
    let half = [0.5, 1.5];
    let q = [0.8, 1.6];
    let aqueous = core::array::from_fn::<_, 2, _>(|i| [v[i] * initial[i]]);
    let bound = core::array::from_fn::<_, 2, _>(|i| [q[i] * initial[i] / (half[i] + initial[i])]);
    let vacant = core::array::from_fn::<_, 2, _>(|i| [q[i] * half[i] / (half[i] + initial[i])]);
    let laws = [
        [binding(LangmuirSoluteLaw::Equilibrium(&a), mass[0], &point)],
        [binding(LangmuirSoluteLaw::Equilibrium(&b), mass[1], &point)],
    ];
    let dt = 2.5;
    let g = 1.0 / (0.4 / (0.2 / 0.8 * 0.04 / 2.0) + 0.6 / (0.4 / 1.2 * 0.09 / 3.0));
    let total = [aqueous[0][0] + bound[0][0], aqueous[1][0] + bound[1][0]];
    let (mut low, mut high) = (0.0, total[0]);
    for _ in 0..100 {
        let j = (low + high) * 0.5;
        let c0 = concentration(v[0], q[0], half[0], total[0] - j);
        let c1 = concentration(v[1], q[1], half[1], total[1] + j);
        if j - dt * g * (c0 - c1) > 0.0 {
            high = j;
        } else {
            low = j;
        }
    }
    let j = (low + high) * 0.5;
    let result = solve(
        dt,
        &v,
        &v,
        &[0.8, 1.2],
        &[0.0; 3],
        &aqueous,
        &bound,
        &vacant,
        &laws,
        0.2,
        &[[0.04], [0.09]],
        &[2.0, 3.0],
        [DiluteSoluteBoundary::Outflow; 2],
        controls(),
        || false,
    )
    .unwrap();
    near(
        result.sorption.report.total.diffusive_face_amount_mol_m2[1][0],
        j,
        1e-11,
    );
    for i in 0..2 {
        let c = concentration(v[i], q[i], half[i], total[i] + if i == 0 { -j } else { j });
        near(result.sorption.concentration_mol_m3[i][0], c, 1e-11);
        near(
            result.sorption.sorbed_amount_mol_m2[i][0] + result.vacant_amount_mol_m2[i][0],
            q[i],
            1e-12,
        );
    }
}

#[test]
fn g1_saturating_breakthrough_and_outlet_integrals_match_stepwise_quadratics() {
    let model = equilibrium(2.0, 0.5, 100.0);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Equilibrium(&model), 0.5, &point)]];
    let (mut aqueous, mut bound, mut vacant) = ([[0.0]], [[0.0]], [[1.0]]);
    let (dt, v, flow, inlet) = (0.5, 0.4, 0.2, 4.0);
    let mut outlet = 0.0;
    let mut oracle_total = 0.0;
    let mut early = 0.0;
    for step in 0..80 {
        let c = concentration(v + dt * flow, 1.0, 0.5, oracle_total + dt * flow * inlet);
        oracle_total = v * c + c / (0.5 + c);
        let result = solve(
            dt,
            &[v],
            &[v],
            &[1.0],
            &[dt * flow; 2],
            &aqueous,
            &bound,
            &vacant,
            &laws,
            1.0,
            &[[0.0]],
            &[1.0],
            [
                DiluteSoluteBoundary::InflowConcentration([inlet]),
                DiluteSoluteBoundary::Outflow,
            ],
            controls(),
            || false,
        )
        .unwrap();
        near(result.sorption.concentration_mol_m3[0][0], c, 2e-10);
        assert!(result.sorption.sorbed_amount_mol_m2[0][0] < 1.0);
        outlet += result.sorption.report.total.advective_face_amount_mol_m2[1][0];
        aqueous[0] = result.sorption.aqueous_amount_mol_m2[0];
        bound[0] = result.sorption.sorbed_amount_mol_m2[0];
        vacant[0] = result.vacant_amount_mol_m2[0];
        near(
            aqueous[0][0] + bound[0][0] + outlet,
            (step + 1) as f64 * dt * flow * inlet,
            2e-9,
        );
        if step == 0 {
            early = c;
        }
    }
    assert!(early < inlet * 0.2);
    near(aqueous[0][0] / v, inlet, 1e-6);
    near(bound[0][0], inlet / (0.5 + inlet), 1e-6);
}

#[test]
fn g3_zero_rate_and_inert_slots_preserve_original_aqueous_arithmetic() {
    let model = kinetic(2.0, 0.5, 0.0);
    let point = conditions();
    let laws = [[
        binding(LangmuirSoluteLaw::Kinetic(&model), 1.0, &point),
        None,
    ]; 2];
    let mobile = [[0.3, 0.1], [0.0, 0.6]];
    let boundary = [
        DiluteSoluteBoundary::Reservoir {
            concentration_mol_m3: [0.5, 0.2],
            film_resistance_s_m: [0.0, 3.0],
        },
        DiluteSoluteBoundary::Outflow,
    ];
    let plain = transport(
        [0.0, 1.0],
        1.0,
        &[0.3; 2],
        &[0.3; 2],
        &[1.0; 2],
        &[0.0; 3],
        &mobile,
        &[[0.1, 0.2]; 2],
        &[1.0; 2],
        boundary,
        controls().transport,
        || false,
    )
    .unwrap();
    let result = solve(
        1.0,
        &[0.3; 2],
        &[0.3; 2],
        &[1.0; 2],
        &[0.0; 3],
        &mobile,
        &[[0.7, 0.0]; 2],
        &[[1.3, 0.0]; 2],
        &laws,
        1.0,
        &[[0.1, 0.2]; 2],
        &[1.0; 2],
        boundary,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.sorption.aqueous_amount_mol_m2, plain.amount_mol_m2);
    assert_eq!(
        result.sorption.concentration_mol_m3,
        plain.concentration_mol_m3
    );
    assert_eq!(
        result.sorption.report.total.diffusive_face_amount_mol_m2,
        plain.report.diffusive_face_amount_mol_m2
    );
    assert_eq!(result.sorption.sorbed_amount_mol_m2, vec![[0.7, 0.0]; 2]);
    assert_eq!(result.vacant_amount_mol_m2, vec![[1.3, 0.0]; 2]);
}

#[test]
fn g3_uniform_concentration_survives_actual_changing_water_storage() {
    let model = equilibrium(2.0, 0.5, 100.0);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Equilibrium(&model), 0.5, &point)]; 3];
    let before = [0.25, 0.6, 0.3];
    let water = [0.02, 0.01, 0.015, 0.015];
    let c = 3.0;
    let after = core::array::from_fn::<_, 3, _>(|i| before[i] + water[i] - water[i + 1]);
    let b = c / (0.5 + c);
    let vacant = 0.5 / (0.5 + c);
    let result = solve(
        1.0,
        &before,
        &after,
        &[1.0; 3],
        &water,
        &before.map(|v| [v * c]),
        &[[b]; 3],
        &[[vacant]; 3],
        &laws,
        1.0,
        &[[0.02]; 3],
        &[1.0, 2.0, 3.0],
        [
            DiluteSoluteBoundary::InflowConcentration([c]),
            DiluteSoluteBoundary::Outflow,
        ],
        controls(),
        || false,
    )
    .unwrap();
    for value in result.sorption.concentration_mol_m3 {
        near(value[0], c, 1e-13);
    }
    near(
        result.sorption.report.total.amount_change_mol_m2[0],
        c * (water[0] - water[3]),
        1e-13,
    );
}

#[test]
fn g1_finite_rate_loading_has_first_order_time_refinement() {
    let model = kinetic(2.0, 0.5, 1.0);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Kinetic(&model), 0.5, &point)]];
    // Independent fixed-step RK4 of the continuous closed-bath scalar ODE.
    let derivative = |b: f64| {
        let c = (0.6 - b) / 0.5;
        c / (0.5 + c) - b
    };
    let mut exact = 0.1;
    let h = 1.0 / 4096.0;
    for _ in 0..4096 {
        let a = derivative(exact);
        let b = derivative(exact + h * a / 2.0);
        let c = derivative(exact + h * b / 2.0);
        let d = derivative(exact + h * c);
        exact += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let errors = [16, 32, 64].map(|steps| {
        let (mut mobile, mut bound, mut vacant) = ([[0.5]], [[0.1]], [[0.9]]);
        for _ in 0..steps {
            let r = solve(
                1.0 / steps as f64,
                &[0.5],
                &[0.5],
                &[1.0],
                &[0.0; 2],
                &mobile,
                &bound,
                &vacant,
                &laws,
                1.0,
                &[[0.0]],
                &[1.0],
                [DiluteSoluteBoundary::Outflow; 2],
                controls(),
                || false,
            )
            .unwrap();
            mobile[0] = r.sorption.aqueous_amount_mol_m2[0];
            bound[0] = r.sorption.sorbed_amount_mol_m2[0];
            vacant[0] = r.vacant_amount_mol_m2[0];
        }
        (bound[0][0] - exact).abs()
    });
    assert!(
        errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8,
        "{errors:?}"
    );
}

#[test]
fn g3_saturated_occupied_stock_keeps_independent_vacant_tail() {
    let model = equilibrium(2.0, 1.0, 1e25);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Equilibrium(&model), 1.0, &point)]];
    let result = solve(
        1.0,
        &[0.5],
        &[0.5],
        &[1.0],
        &[0.0; 2],
        &[[5e19]],
        &[[2.0]],
        &[[2e-20]],
        &laws,
        1.0,
        &[[0.0]],
        &[1.0],
        [DiluteSoluteBoundary::Outflow; 2],
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.sorption.sorbed_amount_mol_m2[0][0], 2.0);
    near(result.vacant_amount_mol_m2[0][0] / 2e-20, 1.0, 1e-15);
    assert_eq!(result.sorption.report.total.linear_solves, 0);
}

#[test]
fn g0_original_sources_actual_volume_equilibrium_and_site_stocks_refuse() {
    let model = kinetic(2.0, 0.5, 1.0);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Kinetic(&model), 1.0, &point)]];
    let call = |area, b, v| {
        solve(
            1.0,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &[[0.2]],
            &[[b]],
            &[[v]],
            &laws,
            area,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false,
        )
    };
    assert!(
        matches!(call(1e-8,0.1,1.9),Err(LangmuirSoluteError::KineticSource {error:AdsorptionFlowError::Source(AqueousAdsorptionError::OutsideValidity{axis,..}),..}) if axis==VOLUME_AXIS)
    );
    assert!(matches!(
        call(1.0, 0.1, 1.8),
        Err(LangmuirSoluteError::Transport(DiluteSoluteError::Invalid(
            _
        )))
    ));
    assert_eq!(site_balance(1e16, 1.0, 1e16).unwrap(), 1.0);
    let large_laws = [[binding(LangmuirSoluteLaw::Kinetic(&model), 5e15, &point)]];
    assert!(matches!(
        solve(
            1.0,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &[[0.2]],
            &[[1e16]],
            &[[1.0]],
            &large_laws,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false
        ),
        Err(LangmuirSoluteError::Transport(DiluteSoluteError::Invalid(
            _
        )))
    ));
    let equilibrium_laws = [[binding(
        LangmuirSoluteLaw::Equilibrium(model.equilibrium()),
        1.0,
        &point,
    )]];
    assert!(matches!(
        solve(
            1.0,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &[[0.2]],
            &[[0.1]],
            &[[1.9]],
            &equilibrium_laws,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false
        ),
        Err(LangmuirSoluteError::Transport(
            DiluteSoluteError::Sorption { .. }
        ))
    ));
    let missing = QueryPoint::new();
    let missing_laws = [[binding(LangmuirSoluteLaw::Kinetic(&model), 1.0, &missing)]];
    assert!(matches!(
        solve(
            1.0,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &[[0.2]],
            &[[0.1]],
            &[[1.9]],
            &missing_laws,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false
        ),
        Err(LangmuirSoluteError::KineticSource { .. })
    ));
}

#[test]
fn g4_nonlinear_work_source_payload_and_cancellation_leave_inputs_intact() {
    let model = kinetic(2.0, 0.3, 1.0);
    let point = conditions();
    let laws = [[binding(LangmuirSoluteLaw::Kinetic(&model), 1.0, &point)]];
    let mobile = [[0.2]];
    let bound = [[0.1]];
    let vacant = [[1.9]];
    let call = |ctl, stop| {
        let mut polls = 0;
        solve(
            1.0,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &mobile,
            &bound,
            &vacant,
            &laws,
            1.0,
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            ctl,
            || {
                polls += 1;
                polls == stop
            },
        )
    };
    let mut ctl = controls();
    ctl.maximum_source_evaluations = 1;
    assert!(matches!(
        call(ctl, usize::MAX),
        Err(LangmuirSoluteError::SourceBudget { evaluations: 1, .. })
    ));
    let mut ctl = controls();
    ctl.maximum_linear_solves = 1;
    assert!(matches!(
        call(ctl, usize::MAX),
        Err(LangmuirSoluteError::Transport(
            DiluteSoluteError::NonlinearWork { .. }
        ))
    ));
    let mut ctl = controls();
    ctl.transport.maximum_storage_bytes = transport::admit_extent(1, ctl.transport).unwrap().1;
    assert!(matches!(
        call(ctl, usize::MAX),
        Err(LangmuirSoluteError::Transport(
            DiluteSoluteError::Workspace { .. }
        ))
    ));
    assert_eq!(
        call(controls(), 8),
        Err(LangmuirSoluteError::Transport(DiluteSoluteError::Cancelled))
    );
    assert_eq!(mobile, [[0.2]]);
    assert_eq!(bound, [[0.1]]);
    assert_eq!(vacant, [[1.9]]);
}
