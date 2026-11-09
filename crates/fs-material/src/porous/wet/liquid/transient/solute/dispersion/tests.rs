//! G1/G3 longitudinal transport, source coordinates and conservative phase coupling.
//! Parameters are authored numerical fixtures, not measured porous calibration.
use super::super::{
    langmuir::{
        LangmuirSoluteControls, LangmuirSoluteLaw, LangmuirSoluteSorption, LangmuirStorage,
    },
    sorption::{LinearSoluteSorption, SorptionStorage},
    transport::{DiluteSoluteBoundary, DiluteSoluteStep, transport, transport_with_dispersion},
};
use super::*;
use crate::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
    CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER, TEMPERATURE_AXIS,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{
    AmountConcentration, Temperature,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn source(
    law: &str,
    params: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let hash = hash_domain(
        "porous-dispersion-test",
        b"authored source-qualified advection-dispersion fixture",
    );
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: params
            .iter()
            .map(|&(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![hash],
        provenance: Provenance {
            source: "Authored transport arithmetic, not measured material calibration".into(),
            license: "project fixture".into(),
            artifact: Some(hash),
        },
    }
}
fn model(alpha: f64) -> LongitudinalDispersivity {
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            PORE_SPEED_AXIS,
            QuantitySpec::dimensional(PoreSpeed::DIMS),
            0.0,
            10.0,
        )
        .with_quantity(
            LIQUID_SATURATION_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            1.0,
        )
        .with_quantity(
            POROSITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.1,
            1.0,
        )
        .with_quantity(
            TRANSPORT_SCALE_AXIS,
            QuantitySpec::dimensional(Length::DIMS),
            1.0,
            100.0,
        )
        .with_quantity(
            "background-pH",
            QuantitySpec::dimensional(Dims::NONE),
            6.0,
            8.0,
        );
    LongitudinalDispersivity::from_card(source(
        LONGITUDINAL_DISPERSION_LAW,
        &[(DISPERSIVITY_PARAMETER, alpha, Length::DIMS)],
        validity,
    ))
    .unwrap()
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("background-pH", QuantitySpec::dimensional(Dims::NONE), 7.0)
        .unwrap()
}
fn binding<'a>(
    model: &'a LongitudinalDispersivity,
    conditions: &'a QueryPoint,
) -> SoluteDispersionCell<'a> {
    SoluteDispersionCell {
        model,
        transport_scale: Length::new(10.0),
        conditions,
    }
}
fn controls<const N: usize>() -> DiluteSoluteControls<N> {
    DiluteSoluteControls {
        absolute_equation_tolerance_mol_m2: [1e-10; N],
        relative_equation_tolerance: 1e-12,
        amount_balance_tolerance_mol_m2: [1e-10; N],
        water_balance_tolerance_m: 1e-12,
        maximum_cell_species: 131072,
        maximum_storage_bytes: 64 * 1024 * 1024,
    }
}
fn preparation_controls() -> SoluteDispersionControls {
    SoluteDispersionControls {
        maximum_source_points: 8192,
        maximum_storage_bytes: 1024 * 1024,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn prepared<'a>(
    dt: f64,
    water: &[f64],
    widths: &[f64],
    faces: &[f64],
    bindings: &'a [SoluteDispersionCell<'a>],
) -> PreparedDispersion<'a> {
    prepare(
        [0.0, dt],
        dt,
        bindings,
        preparation_controls(),
        |i| DispersionCellState {
            porosity: 0.5,
            liquid_saturation: 2.0 * water[i] / widths[i],
            liquid_volume_fraction: water[i] / widths[i],
        },
        |i| Ok(faces[i]),
        || false,
    )
    .unwrap()
}
#[allow(clippy::too_many_arguments)]
fn solve<'a, const N: usize>(
    dt: f64,
    before: &[f64],
    after: &[f64],
    widths: &[f64],
    faces: &[f64],
    mobile: &[[f64; N]],
    d: &[[f64; N]],
    tau: &[f64],
    boundaries: [DiluteSoluteBoundary<N>; 2],
    prepared: PreparedDispersion<'a>,
    ctl: DiluteSoluteControls<N>,
) -> Result<DispersiveSoluteStep<'a, DiluteSoluteStep<N>>> {
    let reduced = prepared.transport_controls(ctl)?;
    let mut step = transport_with_dispersion(
        [0.0, dt],
        dt,
        before,
        after,
        widths,
        faces,
        mobile,
        d,
        tau,
        boundaries,
        reduced,
        &mut (),
        Some(&prepared.report.mechanical_bulk_m2_s),
        || false,
    )
    .map_err(|e| prepared.transport_error(e))?;
    step.report.workspace_bytes += prepared.report.preparation_storage_bytes;
    Ok(prepared.finish(step))
}

#[test]
fn g3_zero_dispersivity_and_stagnation_keep_original_molecular_arithmetic() {
    let point = context();
    let d = [[0.1, 0.2]; 2];
    let mobile = [[0.3, 0.1], [0.0, 0.6]];
    let bounds = [
        DiluteSoluteBoundary::Reservoir {
            concentration_mol_m3: [0.5, 0.2],
            film_resistance_s_m: [0.0, 3.0],
        },
        DiluteSoluteBoundary::Outflow,
    ];
    for (alpha, q) in [(0.0, 0.02), (0.4, 0.0)] {
        let m = model(alpha);
        let bindings = [binding(&m, &point); 2];
        let faces = [q; 3];
        let plain = transport(
            [0.0, 1.0],
            1.0,
            &[0.3; 2],
            &[0.3; 2],
            &[1.0; 2],
            &faces,
            &mobile,
            &d,
            &[1.0, 3.0],
            bounds,
            controls(),
            || false,
        )
        .unwrap();
        let mixed = solve(
            1.0,
            &[0.3; 2],
            &[0.3; 2],
            &[1.0; 2],
            &faces,
            &mobile,
            &d,
            &[1.0, 3.0],
            bounds,
            prepared(1.0, &[0.3; 2], &[1.0; 2], &faces, &bindings),
            controls(),
        )
        .unwrap();
        assert_eq!(mixed.transport.amount_mol_m2, plain.amount_mol_m2);
        assert_eq!(
            mixed.transport.concentration_mol_m3,
            plain.concentration_mol_m3
        );
        assert_eq!(
            mixed.transport.report.advective_face_amount_mol_m2,
            plain.report.advective_face_amount_mol_m2
        );
        assert_eq!(
            mixed.transport.report.diffusive_face_amount_mol_m2,
            plain.report.diffusive_face_amount_mol_m2
        );
        assert_eq!(mixed.dispersion.mechanical_bulk_m2_s, vec![[0.0; 2]; 2]);
        assert_eq!(mixed.models()[0].model.dispersivity().value(), alpha);
    }
}

#[test]
fn g1_opposed_face_inflows_mix_independently_without_molecular_tortuosity() {
    let m = model(0.3);
    let point = context();
    let bindings = [binding(&m, &point)];
    let faces = [0.1, -0.1];
    let boundaries = [
        DiluteSoluteBoundary::Reservoir {
            concentration_mol_m3: [1.0],
            film_resistance_s_m: [0.0],
        },
        DiluteSoluteBoundary::Reservoir {
            concentration_mol_m3: [0.0],
            film_resistance_s_m: [0.0],
        },
    ];
    let mut values = Vec::new();
    for tau in [1.0, 1000.0] {
        let p = prepared(1.0, &[0.4], &[1.0], &faces, &bindings);
        assert_eq!(p.report.pore_speed_m_s, vec![[0.25; 2]]);
        assert_eq!(p.report.mechanical_bulk_m2_s, vec![[0.03; 2]]);
        let r = solve(
            1.0,
            &[0.2],
            &[0.4],
            &[1.0],
            &faces,
            &[[0.0]],
            &[[0.0]],
            &[tau],
            boundaries,
            p,
            controls(),
        )
        .unwrap();
        let c = (0.1 + 0.06) / (0.4 + 0.12);
        near(r.transport.concentration_mol_m3[0][0], c, 1e-15);
        near(
            r.transport.report.diffusive_face_amount_mol_m2[0][0],
            0.06 * (1.0 - c),
            1e-15,
        );
        near(
            r.transport.report.diffusive_face_amount_mol_m2[1][0],
            0.06 * c,
            1e-15,
        );
        near(
            r.transport.report.amount_balance_defect_mol_m2[0],
            0.0,
            1e-15,
        );
        values.push(r.transport);
    }
    assert_eq!(values[0], values[1]);
}

#[test]
fn g1_g3_heterogeneous_half_resistances_and_reflected_flow_match_two_cell_oracle() {
    let a = model(0.2);
    let b = model(0.6);
    let point = context();
    let widths = [0.5, 1.5];
    let water = [0.1, 0.6];
    let initial = [[0.05], [1.2]];
    let d = [[0.01], [0.03]];
    let tau = [2.0, 3.0];
    let dt = 0.7;
    let q = 0.02;
    let faces = [dt * q; 3];
    let g = 1.0 / (0.25 / (0.2 * 0.01 / 2.0 + 0.2 * q) + 0.75 / (0.4 * 0.03 / 3.0 + 0.6 * q));
    let left = 1.0 / (0.25 / (0.2 * 0.01 / 2.0 + 0.2 * q) + 20.0);
    let matrix = [
        water[0] + dt * (q + left + g),
        -dt * g,
        -dt * (q + g),
        water[1] + dt * (q + g),
    ];
    let rhs = [initial[0][0] + dt * (q + left), initial[1][0]];
    let det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    let expected = [
        (rhs[0] * matrix[3] - matrix[1] * rhs[1]) / det,
        (matrix[0] * rhs[1] - matrix[2] * rhs[0]) / det,
    ];
    let reservoir = DiluteSoluteBoundary::Reservoir {
        concentration_mol_m3: [1.0],
        film_resistance_s_m: [20.0],
    };
    let bindings = [binding(&a, &point), binding(&b, &point)];
    let r = solve(
        dt,
        &water,
        &water,
        &widths,
        &faces,
        &initial,
        &d,
        &tau,
        [reservoir, DiluteSoluteBoundary::Outflow],
        prepared(dt, &water, &widths, &faces, &bindings),
        controls(),
    )
    .unwrap();
    for i in 0..2 {
        near(r.transport.concentration_mol_m3[i][0], expected[i], 2e-15);
    }
    near(
        r.transport.report.diffusive_face_amount_mol_m2[1][0],
        dt * g * (expected[0] - expected[1]),
        2e-16,
    );
    near(r.dispersion.pore_speed_m_s[0][1], 0.1, 1e-16);
    near(r.dispersion.pore_speed_m_s[1][0], 0.05, 1e-16);
    let bindings = [binding(&b, &point), binding(&a, &point)];
    let faces = [-dt * q; 3];
    let reflected = solve(
        dt,
        &[water[1], water[0]],
        &[water[1], water[0]],
        &[widths[1], widths[0]],
        &faces,
        &[initial[1], initial[0]],
        &[d[1], d[0]],
        &[tau[1], tau[0]],
        [DiluteSoluteBoundary::Outflow, reservoir],
        prepared(
            dt,
            &[water[1], water[0]],
            &[widths[1], widths[0]],
            &faces,
            &bindings,
        ),
        controls(),
    )
    .unwrap();
    for i in 0..2 {
        near(
            reflected.transport.concentration_mol_m3[1 - i][0],
            expected[i],
            2e-15,
        );
    }
    for face in 0..3 {
        near(
            reflected.transport.report.diffusive_face_amount_mol_m2[2 - face][0],
            -r.transport.report.diffusive_face_amount_mol_m2[face][0],
            2e-16,
        );
    }
}

#[test]
fn g1_advecting_gaussian_resolves_physical_dispersion_and_first_order_error() {
    let m = model(0.1);
    let point = context();
    let (length, theta, velocity, time, sigma, x0) = (10.0, 0.3, 0.5, 0.8, 0.25, 3.0);
    let physical_d = 0.02 / 2.0 + 0.1 * velocity;
    let mut errors = Vec::new();
    for (n, steps) in [(100, 16), (200, 32), (400, 64)] {
        let dx = length / n as f64;
        let dt = time / steps as f64;
        let water = vec![theta * dx; n];
        let widths = vec![dx; n];
        let faces = vec![theta * velocity * dt; n + 1];
        let d = vec![[0.02]; n];
        let tau = vec![2.0; n];
        let bindings = vec![binding(&m, &point); n];
        let mut mobile: Vec<_> = (0..n)
            .map(|i| {
                let x = (i as f64 + 0.5) * dx;
                [theta * dx * (-0.5 * ((x - x0) / sigma).powi(2)).exp()]
            })
            .collect();
        let original_total: f64 = mobile.iter().map(|x| x[0]).sum();
        for _ in 0..steps {
            let r = solve(
                dt,
                &water,
                &water,
                &widths,
                &faces,
                &mobile,
                &d,
                &tau,
                [
                    DiluteSoluteBoundary::InflowConcentration([0.0]),
                    DiluteSoluteBoundary::Outflow,
                ],
                prepared(dt, &water, &widths, &faces, &bindings),
                controls(),
            )
            .unwrap();
            mobile = r.transport.amount_mol_m2;
        }
        let total: f64 = mobile.iter().map(|x| x[0]).sum();
        let mean: f64 = mobile
            .iter()
            .enumerate()
            .map(|(i, a)| (i as f64 + 0.5) * dx * a[0])
            .sum::<f64>()
            / total;
        let variance: f64 = mobile
            .iter()
            .enumerate()
            .map(|(i, a)| (((i as f64 + 0.5) * dx) - mean).powi(2) * a[0])
            .sum::<f64>()
            / total;
        near(total, original_total, 2e-10);
        near(mean, x0 + velocity * time, 2e-8);
        // Exact infinite-grid BE moments include separately identified donor-
        // cell and time-discretization broadening; neither is physical alpha.
        let discrete_variance = sigma * sigma
            + 2.0 * physical_d * time
            + velocity * dx * time
            + velocity * velocity * time * dt;
        near(variance, discrete_variance, 2e-7);
        let exact_variance = sigma * sigma + 2.0 * physical_d * time;
        let error = (0..n)
            .map(|i| {
                let x = (i as f64 + 0.5) * dx;
                let exact = sigma / exact_variance.sqrt()
                    * (-0.5 * (x - x0 - velocity * time).powi(2) / exact_variance).exp();
                dx * (mobile[i][0] / (theta * dx) - exact).powi(2)
            })
            .sum::<f64>()
            .sqrt();
        errors.push(error);
    }
    assert!(
        errors[0] / errors[1] > 1.7 && errors[1] / errors[2] > 1.7,
        "{errors:?}"
    );
    assert!(errors[2] < 0.02, "{errors:?}");
}

#[test]
fn g1_linear_sorbing_dispersion_keeps_physical_flux_and_combined_inventory() {
    let m = model(0.2);
    let point = context();
    let bindings = [binding(&m, &point); 2];
    let dt = 0.5;
    let faces = [0.02; 3];
    let water = [0.3; 2];
    let p = prepared(dt, &water, &[1.0; 2], &faces, &bindings);
    let ctl = p.transport_controls(controls()).unwrap();
    let initial = [[0.3], [0.075]];
    let bound = [[0.2], [0.075]];
    let laws = [
        [LinearSoluteSorption::Equilibrium { capacity_m: 0.2 }],
        [LinearSoluteSorption::Equilibrium { capacity_m: 0.3 }],
    ];
    let mut storage = SorptionStorage::new(2, dt, &bound, &laws, ctl, &mut || false).unwrap();
    let aqueous = transport_with_dispersion(
        [0.0, dt],
        dt,
        &water,
        &water,
        &[1.0; 2],
        &faces,
        &initial,
        &[[0.0]; 2],
        &[1.0; 2],
        [
            DiluteSoluteBoundary::InflowConcentration([0.6]),
            DiluteSoluteBoundary::Outflow,
        ],
        ctl,
        &mut storage,
        Some(&p.report.mechanical_bulk_m2_s),
        || false,
    )
    .unwrap();
    let result = storage.finish(&initial, aqueous, &mut || false).unwrap();
    let (a, b, c, e) = (0.524, -0.004, -0.024, 0.624);
    let (r0, r1) = (0.512, 0.15);
    let determinant = a * e - b * c;
    let expected = [
        (r0 * e - b * r1) / determinant,
        (a * r1 - c * r0) / determinant,
    ];
    for i in 0..2 {
        near(result.concentration_mol_m3[i][0], expected[i], 1e-14);
    }
    near(
        result.report.total.diffusive_face_amount_mol_m2[1][0],
        0.004 * (expected[0] - expected[1]),
        1e-16,
    );
    near(
        result.report.total.amount_change_mol_m2[0],
        0.02 * (0.6 - expected[1]),
        1e-15,
    );
    near(
        result.report.total.amount_balance_defect_mol_m2[0],
        0.0,
        1e-15,
    );
}

#[test]
fn g1_source_langmuir_dispersion_uses_one_nonlinear_total_balance() {
    let m = model(0.3);
    let point = context();
    let bindings = [binding(&m, &point)];
    let dt = 0.4;
    let faces = [0.008; 2];
    let p = prepared(dt, &[0.3], &[1.0], &faces, &bindings);
    let equilibrium = AqueousLangmuir::from_card(
        ElementId::new("Au").unwrap(),
        source(
            AQUEOUS_LANGMUIR_LAW,
            &[
                (CAPACITY_PARAMETER, 2.0, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, 0.5, AmountConcentration::DIMS),
            ],
            ValidityDomain::unconstrained()
                .with_quantity(
                    CONCENTRATION_AXIS,
                    QuantitySpec::dimensional(AmountConcentration::DIMS),
                    0.0,
                    100.0,
                )
                .with_quantity(
                    TEMPERATURE_AXIS,
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    290.0,
                    310.0,
                ),
        ),
    )
    .unwrap();
    let context = QueryPoint::new();
    let laws = [[Some(LangmuirSoluteSorption {
        law: LangmuirSoluteLaw::Equilibrium(&equilibrium),
        dry_mass_kg_m2: 0.25,
        temperature: Temperature::new(298.15),
        conditions: &context,
    })]];
    let ctl = LangmuirSoluteControls {
        transport: p.transport_controls(controls()).unwrap(),
        maximum_iterations: 128,
        maximum_line_searches: 32,
        maximum_linear_solves: 256,
        maximum_source_evaluations: 10000,
    };
    let initial = [[0.15]];
    let bound = [[0.25]];
    let vacant = [[0.25]];
    let mut storage =
        LangmuirStorage::new(1, dt, &bound, &vacant, &laws, 1.0, ctl, &mut || false).unwrap();
    let trial = transport_with_dispersion(
        [0.0, dt],
        dt,
        &[0.3],
        &[0.3],
        &[1.0],
        &faces,
        &initial,
        &[[0.0]],
        &[100.0],
        [
            DiluteSoluteBoundary::Reservoir {
                concentration_mol_m3: [2.0],
                film_resistance_s_m: [0.0],
            },
            DiluteSoluteBoundary::Reservoir {
                concentration_mol_m3: [0.0],
                film_resistance_s_m: [0.0],
            },
        ],
        ctl.transport,
        &mut storage,
        Some(&p.report.mechanical_bulk_m2_s),
        || false,
    );
    let aqueous = storage
        .recover(trial)
        .map_err(|e| p.langmuir_error(e))
        .unwrap();
    let result = storage.finish(&initial, aqueous, &mut || false).unwrap();
    let (a, q, h, total) = (0.3176_f64, 0.5, 0.5, 0.4256);
    let b = a * h + q - total;
    let expected = 2.0 * total * h / ((b * b + 4.0 * a * total * h).sqrt() + b);
    near(result.sorption.concentration_mol_m3[0][0], expected, 1e-10);
    near(
        result.sorption.sorbed_amount_mol_m2[0][0],
        0.5 * expected / (0.5 + expected),
        1e-10,
    );
    near(
        result.sorption.report.total.amount_balance_defect_mol_m2[0],
        0.0,
        1e-10,
    );
    near(
        result.sorption.sorbed_amount_mol_m2[0][0] + result.vacant_amount_mol_m2[0][0],
        0.5,
        1e-14,
    );
}

#[test]
fn g0_original_speed_scale_units_and_context_are_enforced() {
    let m = model(0.3);
    let point = context();
    let mut card = m.source().clone();
    card.parameters
        .get_mut(DISPERSIVITY_PARAMETER)
        .unwrap()
        .dims = Dims::NONE;
    assert!(LongitudinalDispersivity::from_card(card).is_err());
    let mut card = m.source().clone();
    card.validity = card.validity.with_quantity(
        PORE_SPEED_AXIS,
        QuantitySpec::dimensional(PoreSpeed::DIMS),
        0.0,
        0.2,
    );
    let limited = LongitudinalDispersivity::from_card(card).unwrap();
    let bindings = [binding(&limited, &point)];
    let geometry = DispersionCellState {
        porosity: 0.5,
        liquid_saturation: 0.8,
        liquid_volume_fraction: 0.4,
    };
    let call = |bindings: &[SoluteDispersionCell<'_>]| {
        prepare(
            [0.0, 1.0],
            1.0,
            bindings,
            preparation_controls(),
            |_| geometry,
            |i| Ok(if i == 0 { 0.1 } else { -0.1 }),
            || false,
        )
        .map(|_| ())
    };
    assert!(
        matches!(call(&bindings),Err(SoluteDispersionError::OutsideValidity {axis,value,..}) if axis==PORE_SPEED_AXIS && value==0.25)
    );
    let forged = point
        .clone()
        .with_quantity(
            PORE_SPEED_AXIS,
            QuantitySpec::dimensional(PoreSpeed::DIMS),
            0.0,
        )
        .unwrap();
    assert!(matches!(
        call(&[binding(&m, &forged)]),
        Err(SoluteDispersionError::Invalid(_))
    ));
    let missing = QueryPoint::new();
    assert!(matches!(
        call(&[binding(&m, &missing)]),
        Err(SoluteDispersionError::Invalid(_))
    ));
    let mut too_small = binding(&m, &point);
    too_small.transport_scale = Length::new(0.1);
    assert!(
        matches!(call(&[too_small]),Err(SoluteDispersionError::OutsideValidity {axis,..}) if axis==TRANSPORT_SCALE_AXIS)
    );
}

#[test]
fn g4_combined_payload_source_work_and_cancellation_are_bounded() {
    let m = model(0.2);
    let point = context();
    let bindings = [binding(&m, &point)];
    let water = [0.3];
    let widths = [1.0];
    let faces = [0.01; 2];
    let mobile = [[0.2]];
    let p = prepared(1.0, &water, &widths, &faces, &bindings);
    let preparation = p.report.preparation_storage_bytes;
    let base = transport::admit_extent(1, controls::<1>()).unwrap().1;
    let mut ctl = controls();
    ctl.maximum_storage_bytes = preparation + base + core::mem::size_of::<f64>();
    let bounds = [
        DiluteSoluteBoundary::InflowConcentration([0.0]),
        DiluteSoluteBoundary::Outflow,
    ];
    let good = solve(
        1.0,
        &water,
        &water,
        &widths,
        &faces,
        &mobile,
        &[[0.0]],
        &[1.0],
        bounds,
        p,
        ctl,
    )
    .unwrap();
    assert_eq!(
        good.transport.report.workspace_bytes,
        ctl.maximum_storage_bytes
    );
    ctl.maximum_storage_bytes -= 1;
    let fail = solve(
        1.0,
        &water,
        &water,
        &widths,
        &faces,
        &mobile,
        &[[0.0]],
        &[1.0],
        bounds,
        prepared(1.0, &water, &widths, &faces, &bindings),
        ctl,
    );
    assert!(
        matches!(fail,Err(SoluteDispersionError::Transport(DiluteSoluteError::Workspace {required_bytes,maximum_bytes})) if required_bytes==ctl.maximum_storage_bytes+1 && maximum_bytes==ctl.maximum_storage_bytes)
    );
    let geometry = DispersionCellState {
        porosity: 0.5,
        liquid_saturation: 0.6,
        liquid_volume_fraction: 0.3,
    };
    let mut small = preparation_controls();
    small.maximum_source_points = 1;
    assert!(matches!(
        prepare(
            [0.0, 1.0],
            1.0,
            &bindings,
            small,
            |_| geometry,
            |i| Ok(faces[i]),
            || false
        ),
        Err(SoluteDispersionError::SourceWork {
            required_points: 2,
            ..
        })
    ));
    let mut polls = 0;
    assert!(matches!(
        prepare(
            [0.0, 1.0],
            1.0,
            &bindings,
            preparation_controls(),
            |_| geometry,
            |i| Ok(faces[i]),
            || {
                polls += 1;
                polls == 3
            }
        ),
        Err(SoluteDispersionError::Transport(
            DiluteSoluteError::Cancelled
        ))
    ));
    assert_eq!(mobile, [[0.2]]);
    assert_eq!(faces, [0.01; 2]);
}
