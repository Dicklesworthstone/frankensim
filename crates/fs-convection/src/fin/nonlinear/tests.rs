use super::*;
use crate::fin::{FinTip, StraightFinResponse};
use std::cell::Cell;

fn geometry() -> StraightFinGeometry {
    StraightFinGeometry::try_new(0.04, 0.0002, 0.042).unwrap()
}
fn controls(segments: usize) -> NonlinearFinControls {
    NonlinearFinControls {
        segments,
        maximum_segments: 512,
        maximum_iterations: 150,
        maximum_source_evaluations: 200_000,
        absolute_heat_tolerance_w: 1e-8,
        relative_heat_tolerance: 0.0,
        absolute_temperature_tolerance_k: 1e-9,
        relaxation: 0.7,
        maximum_transverse_biot: 1.0,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance={tolerance}"
    );
}
fn constant(
    n: usize,
    base: f64,
    fluid: f64,
    k: f64,
    h: f64,
    ht: Option<f64>,
) -> NonlinearFinResponse<FinPoint> {
    let mut c = controls(n);
    c.relaxation = 1.0;
    solve(
        geometry(),
        base,
        fluid,
        if ht.is_some() {
            NonlinearFinTip::Convection
        } else {
            NonlinearFinTip::Adiabatic
        },
        c,
        |point| {
            Ok::<_, ()>(FinProperties {
                conductivity_w_m_k: k,
                coefficient_w_m2_k: if point.surface == FinSurface::Tip {
                    ht.unwrap()
                } else {
                    h
                },
                source: point,
            })
        },
        |_, _| Ok(k),
        || false,
    )
    .unwrap()
}

#[test]
fn constant_laws_refine_to_original_closed_form_with_separate_tip_and_signed_heat() {
    let analytic = StraightFinResponse::try_new(
        geometry(),
        7.0,
        80.0,
        FinTip::Convection {
            coefficient_w_m2_k: 120.0,
        },
    )
    .unwrap();
    let exact_q = analytic.base_heat_w(350.0, 300.0).unwrap();
    let exact_tip = analytic.temperature_k_at(0.04, 350.0, 300.0).unwrap();
    let mut errors = Vec::new();
    let mut tips = Vec::new();
    for n in [32, 64, 128] {
        let result = constant(n, 350.0, 300.0, 7.0, 80.0, Some(120.0));
        assert_eq!(result.iterations, 2);
        assert_eq!(result.source_evaluations, (2 * n + 2) * 2);
        near(
            result.base_heat_w,
            result.side_heat_w + result.tip_heat_w,
            1e-9,
        );
        assert!(result.sum_absolute_node_residual_w < 1e-8);
        assert_eq!(result.nodes[0].point.position_m, 0.0);
        assert_eq!(result.nodes[n].point.position_m, 0.04);
        assert_eq!(result.tip.as_ref().unwrap().point.surface, FinSurface::Tip);
        for node in &result.nodes {
            assert_eq!(node.point, node.properties.source);
        }
        errors.push((result.base_heat_w - exact_q).abs());
        tips.push((result.nodes[n].point.temperature_k - exact_tip).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
    assert!(
        tips[0] / tips[1] > 3.9 && tips[1] / tips[2] > 3.9,
        "{tips:?}"
    );
    let heating = constant(128, 300.0, 350.0, 7.0, 80.0, Some(120.0));
    let cooling = constant(128, 350.0, 300.0, 7.0, 80.0, Some(120.0));
    assert_eq!(heating.base_heat_w, -cooling.base_heat_w);
    assert_eq!(heating.side_heat_w, -cooling.side_heat_w);
    assert_eq!(heating.tip_heat_w, -cooling.tip_heat_w);
    assert_eq!(
        heating.effective_conductance_w_k,
        cooling.effective_conductance_w_k
    );
    assert_eq!(cooling, constant(128, 350.0, 300.0, 7.0, 80.0, Some(120.0)));
}

#[test]
fn nonlinear_kirchhoff_field_refines_with_actual_local_material_and_film_laws() {
    // Independent continuum solution: w/k0 = theta + a theta^2/2.
    // k=k0(1+a theta), h=h0(1+a theta/2) give A w''=P h0 w/k0.
    // Its insulated-tip w profile is the elementary cosh solution; invert
    // the quadratic analytically instead of sharing the production FV solve.
    let g = geometry();
    let k0 = 7.0;
    let h0 = 80.0;
    let a = 0.02;
    let base = 350.0;
    let fluid = 300.0;
    let theta = base - fluid;
    let w_base = k0 * (theta + 0.5 * a * theta * theta);
    let m: f64 = (h0 * g.exposed_perimeter_m() / (k0 * g.cross_section_m2())).sqrt();
    let exact_q = g.cross_section_m2() * m * w_base * (m * g.length_m()).tanh();
    let w_tip = w_base / (m * g.length_m()).cosh();
    let y = w_tip / k0;
    let exact_tip = fluid + 2.0 * y / (1.0 + (1.0 + 2.0 * a * y).sqrt());
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let result = solve(
            g,
            base,
            fluid,
            NonlinearFinTip::Adiabatic,
            controls(n),
            |p| {
                Ok::<_, ()>(FinProperties {
                    conductivity_w_m_k: k0 * (1.0 + a * (p.temperature_k - fluid)),
                    coefficient_w_m2_k: h0 * (1.0 + 0.5 * a * (p.temperature_k - fluid)),
                    source: p,
                })
            },
            |left, right| Ok(k0 * (1.0 + a * ((left - fluid) + (right - fluid)) * 0.5)),
            || false,
        )
        .unwrap();
        near(result.base_heat_w, result.side_heat_w, 1e-8);
        near(result.nodes[n].point.temperature_k, exact_tip, 0.03);
        assert!(
            result.nodes[n].properties.conductivity_w_m_k
                < result.nodes[0].properties.conductivity_w_m_k
        );
        assert!(
            result.nodes[n].properties.coefficient_w_m2_k
                < result.nodes[0].properties.coefficient_w_m2_k
        );
        errors.push((result.base_heat_w - exact_q).abs());
        assert!(
            result
                .nodes
                .windows(2)
                .all(|w| w[0].point.temperature_k >= w[1].point.temperature_k)
        );
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
    let frozen = StraightFinResponse::try_new(
        g,
        k0 * (1.0 + a * theta),
        h0 * (1.0 + 0.5 * a * theta),
        FinTip::Adiabatic,
    )
    .unwrap()
    .base_heat_w(base, fluid)
    .unwrap();
    assert!((frozen - exact_q).abs() / exact_q > 0.01);
}

#[test]
fn spatial_film_and_nonlinear_tip_only_kirchhoff_balance_have_independent_oracles() {
    let g = geometry();
    let c = 0.4;
    let k = 10.0;
    let base = 350.0;
    let fluid = 300.0;
    // Manufactured theta/50=c+(1-c)(1-x/L)^2. Its k A theta''/P
    // divided by theta supplies a spatial film, never a temperature output.
    let mut control = controls(20);
    control.relaxation = 1.0;
    let exact = solve(
        g,
        base,
        fluid,
        NonlinearFinTip::Adiabatic,
        control,
        |p| {
            let ratio = c + (1.0 - c) * (1.0 - p.position_m / g.length_m()).powi(2);
            Ok::<_, ()>(FinProperties {
                conductivity_w_m_k: k,
                coefficient_w_m2_k: 2.0 * k * g.cross_section_m2() * (1.0 - c)
                    / (g.exposed_perimeter_m() * g.length_m().powi(2) * ratio),
                source: p,
            })
        },
        |_, _| Ok(k),
        || false,
    )
    .unwrap();
    near(
        exact.base_heat_w,
        2.0 * k * g.cross_section_m2() * (base - fluid) * (1.0 - c) / g.length_m(),
        1e-10,
    );
    for node in exact.nodes {
        let expected = fluid
            + (base - fluid)
                * (c + (1.0 - c) * (1.0 - node.point.position_m / g.length_m()).powi(2));
        near(node.point.temperature_k, expected, 2e-11);
    }

    // No lateral exchange: exact integrated Fourier flux equals one tip loss
    // for k(T)=k0(1+a(T-fluid)). This check is independent of axial resolution.
    let a = 0.02;
    let ht = 100.0;
    let evaluate = |t: f64| {
        k / g.length_m() * ((base - t) + 0.5 * a * ((base - fluid).powi(2) - (t - fluid).powi(2)))
            - ht * (t - fluid)
    };
    let (mut lower, mut upper) = (fluid, base);
    for _ in 0..100 {
        let mid = 0.5 * (lower + upper);
        if evaluate(mid) > 0.0 {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    let tip = 0.5 * (lower + upper);
    for n in [1, 16, 64] {
        let solved = solve(
            g,
            base,
            fluid,
            NonlinearFinTip::Convection,
            controls(n),
            |p| {
                Ok::<_, ()>(FinProperties {
                    conductivity_w_m_k: k * (1.0 + a * (p.temperature_k - fluid)),
                    coefficient_w_m2_k: if p.surface == FinSurface::Tip {
                        ht
                    } else {
                        0.0
                    },
                    source: p,
                })
            },
            |left, right| Ok(k * (1.0 + 0.5 * a * ((left - fluid) + (right - fluid)))),
            || false,
        )
        .unwrap();
        near(solved.nodes[n].point.temperature_k, tip, 2e-8);
        near(
            solved.base_heat_w,
            ht * g.cross_section_m2() * (tip - fluid),
            1e-8,
        );
        assert_eq!(solved.side_heat_w, 0.0);
        near(solved.base_heat_w, solved.tip_heat_w, 1e-8);
    }
}

#[test]
fn equilibrium_and_insulated_limits_preserve_sources_without_querying_unused_bulk() {
    let hot = constant(32, 350.0, 300.0, 10.0, 50.0, Some(60.0));
    let equilibrium = constant(32, 300.0, 300.0, 10.0, 50.0, Some(60.0));
    assert_eq!(equilibrium.iterations, 1);
    assert_eq!(equilibrium.base_heat_w, 0.0);
    assert_eq!(equilibrium.side_heat_w, 0.0);
    assert_eq!(equilibrium.tip_heat_w, 0.0);
    near(
        equilibrium.effective_conductance_w_k,
        hot.effective_conductance_w_k,
        1e-14,
    );
    assert!(
        equilibrium
            .nodes
            .iter()
            .all(|n| n.point.temperature_k == 300.0)
    );
    let insulated = solve(
        geometry(),
        350.0,
        100.0,
        NonlinearFinTip::Adiabatic,
        controls(10),
        |p| {
            assert_eq!(p.temperature_k, 350.0);
            Ok::<_, ()>(FinProperties {
                conductivity_w_m_k: 10.0,
                coefficient_w_m2_k: 0.0,
                source: p,
            })
        },
        |a, b| {
            assert_eq!((a, b), (350.0, 350.0));
            Ok(10.0)
        },
        || false,
    )
    .unwrap();
    assert_eq!(insulated.iterations, 1);
    assert_eq!(insulated.effective_conductance_w_k, 0.0);
    assert!(
        insulated
            .nodes
            .iter()
            .all(|n| n.point.temperature_k == 350.0)
    );
}

#[test]
fn unrepresentable_kelvin_gradient_cannot_claim_conducted_heat() {
    // The ratio solve can retain a formal unit-response conductance when k is
    // enormous, but its tiny temperature drop cannot be recorded in Kelvin.
    // The accepted original-source flux must come from that public state.
    let mut control = controls(8);
    control.maximum_iterations = 2;
    control.relaxation = 1.0;
    let outcome = solve(
        geometry(),
        350.0,
        300.0,
        NonlinearFinTip::Adiabatic,
        control,
        |point| {
            Ok::<_, ()>(FinProperties {
                conductivity_w_m_k: 1e20,
                coefficient_w_m2_k: 50.0,
                source: point,
            })
        },
        |_, _| Ok(1e20),
        || false,
    );
    match outcome {
        Err(NonlinearFinError::NotConverged {
            sum_absolute_node_residual_w,
            energy_residual_w,
            ..
        }) => {
            assert!(sum_absolute_node_residual_w > 1.0);
            assert!(energy_residual_w.abs() > 1.0);
        }
        other => panic!("unresolved physical Kelvin drop admitted: {other:?}"),
    }
}

#[test]
fn source_work_cancellation_biot_and_nonlinear_limits_refuse_without_partial_output() {
    let c = controls(8);
    let calls = Cell::new(0usize);
    let local = |p: FinPoint| {
        calls.set(calls.get() + 1);
        Ok::<_, &'static str>(FinProperties {
            conductivity_w_m_k: 10.0,
            coefficient_w_m2_k: 30.0,
            source: p,
        })
    };
    let mean = |_: f64, _: f64| {
        calls.set(calls.get() + 1);
        Ok(10.0)
    };
    let mut work = c;
    work.maximum_source_evaluations = 16;
    assert!(matches!(
        solve(
            geometry(),
            350.0,
            300.0,
            NonlinearFinTip::Adiabatic,
            work,
            local,
            mean,
            || false
        ),
        Err(NonlinearFinError::WorkLimit {
            source_evaluations: 16
        })
    ));
    assert_eq!(calls.get(), 16);
    calls.set(0);
    assert!(matches!(
        solve(
            geometry(),
            350.0,
            300.0,
            NonlinearFinTip::Adiabatic,
            c,
            local,
            mean,
            || calls.get() == 12
        ),
        Err(NonlinearFinError::Cancelled)
    ));
    assert_eq!(calls.get(), 12);
    let mut biot = c;
    biot.maximum_transverse_biot = 0.0;
    assert!(matches!(
        solve(
            geometry(),
            350.0,
            300.0,
            NonlinearFinTip::Adiabatic,
            biot,
            local,
            mean,
            || false
        ),
        Err(NonlinearFinError::TransverseBiot { .. })
    ));
    let mut short = c;
    short.maximum_iterations = 1;
    assert!(matches!(
        solve(
            geometry(),
            350.0,
            300.0,
            NonlinearFinTip::Adiabatic,
            short,
            local,
            mean,
            || false
        ),
        Err(NonlinearFinError::NotConverged { iterations: 1, .. })
    ));
    assert!(matches!(
        solve(
            geometry(),
            350.0,
            300.0,
            NonlinearFinTip::Adiabatic,
            c,
            |p| if p.position_m == geometry().length_m() {
                Err("original tip domain")
            } else {
                local(p)
            },
            mean,
            || false
        ),
        Err(NonlinearFinError::Source {
            error: "original tip domain",
            ..
        })
    ));
    let complete = solve(
        geometry(),
        350.0,
        300.0,
        NonlinearFinTip::Adiabatic,
        c,
        local,
        mean,
        || false,
    )
    .unwrap();
    let retry = solve(
        geometry(),
        350.0,
        300.0,
        NonlinearFinTip::Adiabatic,
        c,
        local,
        mean,
        || false,
    )
    .unwrap();
    assert_eq!(complete, retry);
    let mut invalid = c;
    invalid.segments = invalid.maximum_segments + 1;
    assert!(matches!(
        solve(
            geometry(),
            350.0,
            300.0,
            NonlinearFinTip::Adiabatic,
            invalid,
            local,
            mean,
            || false
        ),
        Err(NonlinearFinError::Controls(_))
    ));
}
