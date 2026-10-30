//! Original P1 motion, source rheology and heat survive one native transfer.
//! Manufactured source tables isolate the equations, not measured materials.

use super::*;
use crate::material::free_surface::enthalpy::FreeSurfaceEnthalpyError;
use crate::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_qty::Frequency;

fn source_card(rate_slope: f64, maximum_rate: f64) -> MaterialCard {
    let original = thermal_card(1000.0, 2e-5);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        let viscosity = claim.key.name() == "dynamic_viscosity";
        if viscosity || claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [
                    [300.0, 0.0],
                    [500.0, 0.0],
                    [300.0, maximum_rate],
                    [500.0, maximum_rate],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        if viscosity {
                            viscosity_at(p[0], p[1], rate_slope)
                        } else {
                            0.5 + 0.001 * (p[0] - 400.0) + 0.01 * p[1]
                        },
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claim.validity = claim.validity.with_quantity(
            "rate",
            QuantitySpec::dimensional(Frequency::DIMS),
            0.0,
            maximum_rate,
        );
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn viscosity_at(temperature: f64, rate: f64, slope: f64) -> f64 {
    0.1 - 0.0004 * (temperature - 300.0) + slope * rate
}

fn shear_point(rate: f64) -> QueryPoint {
    point()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), rate)
        .unwrap()
}

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

fn bounds() -> TetEnthalpyTransferLimits {
    TetEnthalpyTransferLimits {
        max_source_queries: 512,
        ..thermal_limits()
    }
}

#[test]
fn original_p1_shear_law_drives_receiving_motion_and_molecular_heat() {
    let (p, t, _) = cube([0.0; 3]);
    let v: Vec<_> = p
        .iter()
        .map(|x| [0.01 + 2.0 * (x[1] - 0.004), 0.0, 0.0])
        .collect();
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let before = (p.clone(), t.clone(), v.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    let mut outcomes = Vec::new();
    with_cx(&gate, |cx| {
        for slope in [0.0, 0.02] {
            let card = source_card(slope, 100.0);
            // This authored reference differs from both the actual P1 rate (2)
            // and the initially equilibrated receiving native rate (near zero).
            let law = LiquidTemperatureLaw::try_new(
                &card,
                &shear_point(1.5),
                "T",
                LiquidViscosityBasis::Dynamic,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                law.reference_state().material().clone(),
                LiquidViscosityBasis::Dynamic,
                "T",
                "rate",
            )
            .unwrap();
            let caloric = chart(&card, 1000.0);
            let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1).unwrap();
            assert!(
                scaling
                    .free_surface3_from_tetrahedra_with_enthalpy(
                        domain(),
                        [0.0; 3],
                        source(&p, &t, &v, &h),
                        &law,
                        &caloric,
                        bounds(),
                        Default::default(),
                        cx,
                        &gate,
                    )
                    .is_err()
            );
            let (thermal, report) = scaling
                .free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
                    domain(),
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    Default::default(),
                    solve(),
                    cx,
                    &gate,
                )
                .unwrap();
            let mut thermal = thermal.with_time_origin(1.25).unwrap();
            assert_eq!(report.source_shear_rates_per_s.len(), t.len());
            for rate in &report.source_shear_rates_per_s {
                close(*rate, 2.0, 1e-13);
            }
            assert!(report.source_queries > 4 * t.len());
            assert_eq!(report.mechanical.components.len(), 1);
            assert_eq!(report.components.len(), 1);
            let mass = 1000.0 * 0.004_f64.powi(3);
            let enthalpy = mass * affine_h([0.004; 3]);
            close(report.mechanical.source.mass_kg, mass, 1e-17);
            close(report.mechanical.target.mass_kg, mass, 1e-17);
            close(report.source_enthalpy_j, enthalpy, 1e-12);
            close(report.target_enthalpy_j, enthalpy, 1e-11);
            close(report.components[0].source_enthalpy_j, enthalpy, 1e-12);
            close(report.components[0].target_enthalpy_j, enthalpy, 1e-11);
            for axis in 0..3 {
                close(
                    report.mechanical.target.momentum_kg_m_s[axis],
                    report.mechanical.source.momentum_kg_m_s[axis],
                    1e-16,
                );
                close(
                    report.mechanical.target.angular_momentum_kg_m2_s[axis],
                    report.mechanical.source.angular_momentum_kg_m2_s[axis],
                    1e-18,
                );
            }
            assert!(thermal.rheology_at([3, 3, 3]).unwrap().shear_rate_per_s < 1e-8);
            let mut heat = 0.0;
            let mut maximum_native_rate = 0.0_f64;
            for _ in 0..4 {
                let step = thermal.step_with_viscous_heating(&[], &gate).unwrap();
                let mut expected_heat = 0.0;
                for sample in &step.applied {
                    let point = sample.material.material().query();
                    let temperature = point.axes()["T"];
                    let query_rate = point.axes()["rate"];
                    close(query_rate, sample.shear_rate_per_s, 1e-7);
                    let mu = viscosity_at(temperature, query_rate, slope);
                    close(sample.material.dynamic_viscosity_pa_s(), mu, 1e-11);
                    expected_heat +=
                        mu * sample.shear_rate_per_s.powi(2) * sample.liquid_volume_m3 * 0.0005;
                    maximum_native_rate = maximum_native_rate.max(sample.shear_rate_per_s);
                    assert_eq!(
                        sample.material.material().card_identity(),
                        card.content_hash()
                    );
                }
                close(step.flow.viscous_heat_j, expected_heat, 1e-12);
                close(step.flow.external_heat_j, 0.0, 0.0);
                close(step.flow.balance_defect_j, 0.0, 1e-9);
                heat += step.flow.viscous_heat_j;
            }
            assert!(maximum_native_rate > 1e-3);
            assert!(heat > 0.0);
            close(thermal.total_enthalpy_j().unwrap(), enthalpy + heat, 1e-8);
            close(thermal.physical_time_s(), 1.252, 1e-14);
            let native = thermal.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
            outcomes.push((native.velocity_m_s, heat));
        }
    });
    assert!(
        outcomes[0]
            .0
            .iter()
            .zip(outcomes[1].0)
            .any(|(a, b)| (a - b).abs() > 1e-10)
    );
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 1e-13);
    assert_eq!((p, t, v, h), before);
}

#[test]
fn actual_p1_rate_not_reference_rate_controls_source_admission() {
    let card = source_card(0.02, 1.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &shear_point(0.5),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(),
        LiquidViscosityBasis::Dynamic,
        "T",
        "rate",
    )
    .unwrap();
    let caloric = chart(&card, 1000.0);
    let (p, t, rotation) = cube([0.0; 3]);
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let shear: Vec<_> = p
        .iter()
        .map(|x| [0.01 + 2.0 * (x[1] - 0.004), 0.0, 0.0])
        .collect();
    let before = (p.clone(), rotation.clone(), shear.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1).unwrap();
        let call = |v: &[[f64; 3]], limits, controls, rheology| {
            scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
                domain(),
                [0.0; 3],
                source(&p, &t, v, &h),
                &law,
                &caloric,
                limits,
                controls,
                rheology,
                cx,
                &gate,
            )
        };
        assert!(matches!(
            call(&shear, bounds(), Default::default(), solve()),
            Err(TetEnthalpyTransferError::Enthalpy(
                FreeSurfaceEnthalpyError::Material(_)
            ))
        ));
        let (_, report) = call(&rotation, bounds(), Default::default(), solve()).unwrap();
        assert!(
            report
                .source_shear_rates_per_s
                .iter()
                .all(|rate| *rate < 1e-12)
        );
        let mut short = bounds();
        short.max_source_queries = 2;
        assert!(call(&rotation, short, Default::default(), solve()).is_err());
        assert!(
            call(
                &rotation,
                bounds(),
                Default::default(),
                RheologySolve {
                    maximum_iterations: 0,
                    ..solve()
                }
            )
            .is_err()
        );
        // Admission of actual P1 sources has passed, then the receiving native
        // thermal owner refuses atomically without dropping its shear law.
        let mut late = FreeSurfaceEnthalpyControls::default();
        late.max_cells = 1;
        assert!(call(&rotation, bounds(), late, solve()).is_err());
        gate.request();
        assert!(matches!(
            call(&rotation, bounds(), Default::default(), solve()),
            Err(TetEnthalpyTransferError::Transfer(
                TetLiquidTransferError::Cancelled
            ))
        ));
    });
    assert_eq!((p, rotation, shear, h), before);
}
