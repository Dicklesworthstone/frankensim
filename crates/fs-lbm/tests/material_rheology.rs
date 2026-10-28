//! G1/G3: the real thermal/flow kernels consume a source-backed viscosity
//! curve at a self-consistent physical shear rate. Data below is synthetic;
//! the independent moment and quadratic references are numerical evidence.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_lbm::material::transport::rheology::{RheologyPoint, RheologySolve};
use fs_lbm::material::{
    MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3,
    plan_scaling_from_material,
};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::{MaterialError, fluid_transport::FluidTransportState};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Frequency, QuantitySpec, Temperature};

const DX: f64 = 0.01;
const DT: f64 = 0.01;
const E2: [[f64; 2]; 9] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];
const W2: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1000.0, mu).unwrap()
}
fn controls() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 100,
        relaxation: 0.75,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-12,
    }
}
fn limits() -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: 16,
        temperature_range_k: [290.0, 330.0],
        balance_tolerance_j: 1e-8,
    }
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&fluid(1.0), DX, DT, 0.1)
        .unwrap()
        .thermal_slab(MaterialThermalSlab {
            nx: 2,
            ny: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 16,
        })
        .unwrap()
}
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&fluid(1.0), DX, DT, 0.1)
        .unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2,
            ny: 1,
            nz: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 16,
        })
        .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

// Independently prescribe kinetic first/second moments with the Hermite
// projection, not by invoking the production force or strain-rate helpers.
// D_xy=D_yx=gamma/2, so sqrt(2 D:D)=|gamma|. The transformed distribution's
// first moment is rho*u-F/2, and its nonequilibrium second moment is
// -2 rho cs^2 tau D -(u F^T+F u^T)/2.
fn shear_populations<const D: usize, const Q: usize>(
    equilibrium: [f64; Q],
    e: [[f64; D]; Q],
    w: [f64; Q],
    rho: f64,
    u: [f64; D],
    force: [f64; D],
    gamma: f64,
    tau: f64,
) -> [f64; Q] {
    std::array::from_fn(|q| {
        let momentum: f64 = (0..D).map(|d| -0.5 * e[q][d] * force[d]).sum();
        let mut contraction = 0.0;
        for a in 0..D {
            for b in 0..D {
                let strain = if (a == 0 && b == 1) || (a == 1 && b == 0) {
                    gamma / 2.0
                } else {
                    0.0
                };
                let moment =
                    -2.0 * rho * tau * strain / 3.0 - (u[a] * force[b] + force[a] * u[b]) / 2.0;
                contraction += (e[q][a] * e[q][b] - if a == b { 1.0 / 3.0 } else { 0.0 }) * moment;
            }
        }
        equilibrium[q] + w[q] * (3.0 * momentum + 4.5 * contraction)
    })
}
fn set_shear2(slab: &mut MaterialSlab, gamma: f64, force: [f64; 2]) {
    let (rho, u, tau) = (1.0, [0.03, -0.02], 0.8);
    slab.lattice.grid.f[6] = shear_populations(
        fs_lbm::equilibrium(rho, u[0], u[1]),
        E2,
        W2,
        rho,
        u,
        force,
        gamma,
        tau,
    );
    slab.lattice.grid.fext[6] = force;
}
fn set_shear3(slab: &mut MaterialSlab3, gamma: f64, force: [f64; 3]) {
    let (rho, u, tau) = (1.0, [0.03, -0.02, 0.01], 0.8);
    slab.lattice.grid.f[6] = shear_populations(
        fs_lbm::d3q19::equilibrium3(rho, u),
        fs_lbm::d3q19::E3.map(|e| [f64::from(e.0), f64::from(e.1), f64::from(e.2)]),
        fs_lbm::d3q19::W3,
        rho,
        u,
        force,
        gamma,
        tau,
    );
    slab.lattice.grid.fext[6] = force;
}

#[test]
fn g1_constant_viscosity_is_the_existing_step_and_guo_force_does_not_fake_shear() {
    let gamma = 0.02;
    for force in [[0.0, 0.0, 0.0], [0.001, -0.0007, 0.0003]] {
        let mut a = slab3();
        let mut b = slab3();
        set_shear3(&mut a, gamma, force);
        set_shear3(&mut b, gamma, force);
        let old = a
            .step_with_material_transport(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                &CancelGate::new(),
                |_, _, _| Ok(fluid(1.0)),
            )
            .unwrap();
        let new = b
            .step_with_material_rheology(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                &CancelGate::new(),
                |_, _| Ok(fluid(1.0)),
            )
            .unwrap();
        assert_eq!(new.heat, old.heat);
        assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
        assert_eq!(a.lattice.gpop, b.lattice.gpop);
        near(
            new.applied
                .iter()
                .find(|s| s.cell == 6)
                .unwrap()
                .material
                .point()
                .shear_rate_per_s,
            gamma / DT,
            1e-12,
        );

        let mut a = slab2();
        let mut b = slab2();
        let force = [force[0], force[1]];
        set_shear2(&mut a, gamma, force);
        set_shear2(&mut b, gamma, force);
        let old = a
            .step_with_material_transport(
                &[0.0; 16],
                &[[0.0; 2]; 16],
                DX,
                limits(),
                &CancelGate::new(),
                |_, _, _| Ok(fluid(1.0)),
            )
            .unwrap();
        let new = b
            .step_with_material_rheology(
                &[0.0; 16],
                &[[0.0; 2]; 16],
                DX,
                limits(),
                controls(),
                &CancelGate::new(),
                |_, _| Ok(fluid(1.0)),
            )
            .unwrap();
        assert_eq!(new.heat, old.heat);
        assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
        assert_eq!(a.lattice.gpop, b.lattice.gpop);
        near(
            new.applied
                .iter()
                .find(|s| s.cell == 6)
                .unwrap()
                .material
                .point()
                .shear_rate_per_s,
            gamma / DT,
            1e-12,
        );
    }
}

fn query(temperature: f64, shear_rate: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "shear_rate",
            QuantitySpec::dimensional(Frequency::DIMS),
            shear_rate,
        )
        .unwrap()
}
fn card(slope: f64) -> MaterialCard {
    card_curve(
        "shear_rate",
        Frequency::DIMS,
        vec![(0.0, 1.0), (100.0, 1.0 + 100.0 * slope)],
        1.0,
    )
}
fn card_curve(
    axis: &str,
    axis_dims: Dims,
    knots: Vec<(f64, f64)>,
    thermal_scale: f64,
) -> MaterialCard {
    card_value(
        PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: axis_dims,
            knots,
            dims: Dims([-1, 1, -1, 0, 0, 0]),
        },
        thermal_scale,
    )
}
fn card_value(viscosity_value: PropertyValue, thermal_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), 1000.0),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            1000.0 * thermal_scale,
        ),
        (
            "thermal_conductivity",
            Dims([1, 1, -3, -1, 0, 0]),
            1000.0 * thermal_scale,
        ),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), 1.0),
    ] {
        let viscosity = name == "dynamic_viscosity";
        claims
            .insert_claim(PropertyClaim {
                key: if name == "thermal_conductivity" {
                    PropertyKey::with_quantity(
                        name,
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::ThermalConductivity,
                            ValueForm::Static,
                        )),
                    )
                } else {
                    PropertyKey::new(name, dims)
                },
                value: if viscosity {
                    viscosity_value.clone()
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        290.0,
                        330.0,
                    )
                    .with_quantity(
                        "shear_rate",
                        QuantitySpec::dimensional(Frequency::DIMS),
                        0.0,
                        100.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if viscosity {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic rheology fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(
    card: &MaterialCard,
    reference: &ResolvedLiquidState,
    point: RheologyPoint,
) -> Result<ResolvedLiquidState, MaterialError> {
    reference
        .at_query_point(card, &query(point.temperature_k, point.shear_rate_per_s))
        .map_err(|error| MaterialError::State {
            what: error.to_string(),
        })
}

fn card_surface(temperature_slope: f64, shear_slope: f64, thermal_scale: f64) -> MaterialCard {
    card_value(
        PropertyValue::Surface2D {
            axes: ["T".into(), "shear_rate".into()],
            axis_dims: [Temperature::DIMS, Frequency::DIMS],
            samples: [[290.0, 0.0], [330.0, 0.0], [290.0, 100.0], [330.0, 100.0]]
                .into_iter()
                .map(|point| {
                    (
                        point,
                        1.0 + temperature_slope * (point[0] - 300.0) + shear_slope * point[1],
                    )
                })
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: Dims([-1, 1, -1, 0, 0, 0]),
        },
        thermal_scale,
    )
}

#[test]
fn g1_g3_two_axis_source_changes_actual_flow_at_the_constitutive_root() {
    let gamma = 0.02;
    let mut populations = Vec::new();
    for (temperature_slope, shear_slope) in [(0.0, 0.0), (0.02, 0.0), (0.02, 0.2)] {
        let source = card_surface(temperature_slope, shear_slope, 1.0);
        let reference = resolve_liquid_state(
            &source,
            &query(305.0, 0.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut slab = slab3();
        set_shear3(&mut slab, gamma, [0.0; 3]);
        slab.lattice.grid.tau[6] = 1.7;
        let step = slab
            .step_with_material_rheology(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                &CancelGate::new(),
                |point, _| resolve(&source, &reference, point),
            )
            .unwrap();
        let sample = step.applied.iter().find(|s| s.cell == 6).unwrap();
        // Independent moment law: rate=C/(tau*dt), C=0.8*gamma. For the
        // manufactured mu=1+a(T-300)+b*rate, tau solves tau²-A*tau-B=0.
        let a = 0.5 + 0.3 * (1.0 + temperature_slope * (sample.temperature_k - 300.0));
        let b = 0.3 * shear_slope * 0.8 * gamma / DT;
        let expected_tau = 0.5 * (a + (a * a + 4.0 * b).sqrt());
        near(sample.momentum_relaxation, expected_tau, 2e-12);
        near(
            sample.material.point().shear_rate_per_s,
            0.8 * gamma / (DT * expected_tau),
            2e-11,
        );
        for sample in step.applied.iter().chain(&step.endpoint) {
            let material = sample.material.source().material();
            assert_eq!(material.card_identity(), source.content_hash());
            assert_eq!(material.properties().len(), 4);
            assert_eq!(
                material.query(),
                query(
                    sample.temperature_k,
                    sample.material.point().shear_rate_per_s
                )
            );
            let viscosity = material.property("dynamic_viscosity").unwrap();
            assert_eq!(viscosity.source_axes().len(), 2);
            assert!(!viscosity.is_constant_at_fixed_temperature());
            assert_eq!(
                viscosity.answer().receipt.selected,
                reference
                    .material()
                    .property("dynamic_viscosity")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            near(
                sample.transport.dynamic_viscosity_pa_s(),
                1.0 + temperature_slope * (sample.temperature_k - 300.0)
                    + shear_slope * sample.material.point().shear_rate_per_s,
                2e-12,
            );
        }
        populations.push(slab.lattice.grid.f);
    }
    assert_ne!(
        populations[0], populations[1],
        "temperature law must change momentum"
    );
    assert_ne!(
        populations[1], populations[2],
        "shear law must change momentum"
    );
}

#[test]
fn g3_two_axis_viscosity_and_viscous_heat_share_the_same_evolving_material_state() {
    let source = card_surface(0.02, 0.2, 0.001);
    let reference = resolve_liquid_state(
        &source,
        &query(305.0, 0.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let make = || {
        plan_scaling_from_material(&reference, DX, DT, 0.1)
            .unwrap()
            .thermal_slab3(MaterialThermalSlab3 {
                nx: 2,
                ny: 1,
                nz: 6,
                bottom_temperature_k: 310.0,
                top_temperature_k: 300.0,
                gravity_m_s2: 0.0,
                expansion_per_k: 0.0,
                max_cells: 16,
            })
            .unwrap()
    };
    let mut hot = make();
    let mut cold = make();
    set_shear3(&mut hot, 0.02, [0.0; 3]);
    set_shear3(&mut cold, 0.02, [0.0; 3]);
    let gate = CancelGate::new();
    let heated = hot
        .step_with_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |point, _| resolve(&source, &reference, point),
        )
        .unwrap();
    let unheated = cold
        .step_with_material_rheology(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |point, _| resolve(&source, &reference, point),
        )
        .unwrap();
    assert_eq!(hot.lattice.grid.f, cold.lattice.grid.f);
    let warm = heated
        .transport
        .endpoint
        .iter()
        .find(|s| s.cell == 6)
        .unwrap();
    let cool = unheated.endpoint.iter().find(|s| s.cell == 6).unwrap();
    assert!(warm.temperature_k > cool.temperature_k);
    near(
        warm.transport.dynamic_viscosity_pa_s() - cool.transport.dynamic_viscosity_pa_s(),
        0.02 * (warm.temperature_k - cool.temperature_k)
            + 0.2
                * (warm.material.point().shear_rate_per_s - cool.material.point().shear_rate_per_s),
        2e-12,
    );
    let applied = heated
        .transport
        .applied
        .iter()
        .find(|s| s.cell == 6)
        .unwrap();
    near(
        heated.viscous.power_w_m3[6],
        applied.transport.dynamic_viscosity_pa_s()
            * applied.material.point().shear_rate_per_s.powi(2),
        2e-12,
    );
    hot.step_with_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &gate,
        |point, _| resolve(&source, &reference, point),
    )
    .unwrap();
    cold.step_with_material_rheology(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &gate,
        |point, _| resolve(&source, &reference, point),
    )
    .unwrap();
    assert_ne!(hot.lattice.grid.f, cold.lattice.grid.f);
}

#[test]
fn g1_g3_card_curve_changes_flow_and_solves_the_constitutive_quadratic() {
    let gamma = 0.02;
    let mut results = Vec::new();
    for slope in [0.0, 0.2] {
        let source = card(slope);
        let reference = resolve_liquid_state(
            &source,
            &query(305.0, 0.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut slab = slab3();
        set_shear3(&mut slab, gamma, [0.0; 3]);
        // Deliberately unrelated background: it is an initial iterate, not
        // the viscosity that may be used to interpret the evolving stress.
        slab.lattice.grid.tau[6] = 1.7;
        let step = slab
            .step_with_material_rheology(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                &CancelGate::new(),
                |point, _| resolve(&source, &reference, point),
            )
            .unwrap();
        let sample = step.applied.iter().find(|s| s.cell == 6).unwrap();
        // tau = 0.8 + 0.3*slope*C/(tau*dt), C=0.8*gamma.
        let expected_tau = 0.5 * (0.8 + (0.64 + 4.0 * 0.3 * slope * 0.8 * gamma / DT).sqrt());
        near(sample.momentum_relaxation, expected_tau, 2e-12);
        near(
            sample.material.point().shear_rate_per_s,
            0.8 * gamma / (expected_tau * DT),
            2e-11,
        );
        for sample in step.applied.iter().chain(&step.endpoint) {
            assert_eq!(
                sample.material.source().material().card_identity(),
                source.content_hash()
            );
            assert_eq!(sample.material.source().material().properties().len(), 4);
            assert_eq!(
                sample.material.source().material().query(),
                query(
                    sample.temperature_k,
                    sample.material.point().shear_rate_per_s
                )
            );
            near(
                sample.transport.dynamic_viscosity_pa_s(),
                1.0 + slope * sample.material.point().shear_rate_per_s,
                1e-12,
            );
        }
        results.push(slab.lattice.grid.f);
    }
    assert_ne!(
        results[0], results[1],
        "changing only constitutive data must change actual momentum evolution"
    );
}

#[test]
fn g1_g3_strong_source_thickening_matches_constitutive_root_and_collision() {
    let source = card_curve(
        "shear_rate",
        Frequency::DIMS,
        vec![(0.0, 1.0), (1.2, 1.0), (1.4, 3.0), (100.0, 3.0)],
        1.0,
    );
    let reference = resolve_liquid_state(
        &source,
        &query(305.0, 0.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut slab = slab3();
    let mut prescribed = slab3();
    set_shear3(&mut slab, 0.02, [0.0; 3]);
    set_shear3(&mut prescribed, 0.02, [0.0; 3]);

    // The independently imposed stress has gamma*tau=1.6 /s. In the
    // source transition mu=10*gamma-11 Pa s, so tau^2+2.8*tau-4.8=0.
    // The unique positive root is tau=1.2, gamma=4/3 /s and mu=7/3 Pa s.
    // The former damped fixed point has slope -2 at this admissible root.
    prescribed
        .step_with_material_transport(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            &CancelGate::new(),
            |cell, _, _| Ok(fluid(if cell == 6 { 7.0 / 3.0 } else { 1.0 })),
        )
        .unwrap();
    let step = slab
        .step_with_material_rheology(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &CancelGate::new(),
            |point, _| resolve(&source, &reference, point),
        )
        .unwrap();
    let sample = step.applied.iter().find(|sample| sample.cell == 6).unwrap();
    near(sample.momentum_relaxation, 1.2, 2e-12);
    near(sample.material.point().shear_rate_per_s, 4.0 / 3.0, 2e-12);
    near(sample.transport.dynamic_viscosity_pa_s(), 7.0 / 3.0, 2e-11);
    assert!(sample.material.iterations() <= 20);
    for sample in step.applied.iter().chain(&step.endpoint) {
        assert_eq!(
            sample.material.source().material().card_identity(),
            source.content_hash()
        );
        assert_eq!(
            sample.material.source().material().query(),
            query(
                sample.temperature_k,
                sample.material.point().shear_rate_per_s
            )
        );
        assert!(sample.material.relaxation_defect() <= 2e-12);
    }
    for (actual, expected) in slab.lattice.grid.f.iter().zip(&prescribed.lattice.grid.f) {
        for (&actual, &expected) in actual.iter().zip(expected) {
            near(actual, expected, 2e-13);
        }
    }
    for (actual, expected) in slab.lattice.gpop.iter().zip(&prescribed.lattice.gpop) {
        for (&actual, &expected) in actual.iter().zip(expected) {
            near(actual, expected, 2e-13);
        }
    }
}

#[test]
fn nonconvergence_source_refusal_and_cancelled_local_solve_are_atomic() {
    let mut slab = slab3();
    set_shear3(&mut slab, 0.02, [0.0; 3]);
    let (f, g, tau, loads) = (
        slab.lattice.grid.f.clone(),
        slab.lattice.gpop.clone(),
        slab.lattice.grid.tau.clone(),
        slab.lattice.grid.fext.clone(),
    );
    let result = slab.step_with_material_rheology(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        RheologySolve {
            maximum_iterations: 1,
            ..controls()
        },
        &CancelGate::new(),
        |point, _| Ok(fluid(1.0 + point.shear_rate_per_s)),
    );
    assert!(matches!(
        result,
        Err(MaterialTransportError::RheologyConvergence {
            cell: 6,
            iterations: 1,
            ..
        })
    ));
    let mut calls = 0;
    let result = slab.step_with_material_rheology(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &CancelGate::new(),
        |_, _| {
            calls += 1;
            if calls == 13 {
                Err(MaterialError::State {
                    what: "endpoint outside source support".into(),
                })
            } else {
                Ok(fluid(1.0))
            }
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::State { .. })));
    let gate = CancelGate::new();
    let result = slab.step_with_material_rheology(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &gate,
        |point, _| {
            if point.cell == 6 {
                gate.request();
            }
            Ok(fluid(2.0))
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::Cancelled)));
    assert_eq!(slab.lattice.grid.f, f);
    assert_eq!(slab.lattice.gpop, g);
    assert_eq!(slab.lattice.grid.tau, tau);
    assert_eq!(slab.lattice.grid.fext, loads);
}

#[test]
fn g1_viscous_heat_matches_independent_shear_power_and_existing_heat_transport() {
    let gamma = 0.02;
    let expected_power = (gamma / DT) * (gamma / DT); // mu = 1 Pa s
    let mut authored = [0.0; 16];
    authored[8] = 17.0;
    authored[10] = -5.0;
    let mut prescribed = authored;
    prescribed[6] += expected_power;
    let force = [0.001, -0.0007, 0.0003];
    let gate = CancelGate::new();

    let mut hot = slab3();
    let mut manual = slab3();
    set_shear3(&mut hot, gamma, force);
    set_shear3(&mut manual, gamma, force);
    let heated = hot
        .step_with_material_rheology_and_viscous_heating(
            &authored,
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    let expected = manual
        .step_with_material_rheology(
            &prescribed,
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    let check = |heat: &fs_lbm::material::transport::rheology::heating::ViscousHeating,
                 balance: &fs_lbm::thermal::sources::HeatBalance,
                 reference: &fs_lbm::thermal::sources::HeatBalance,
                 volume: f64| {
        near(heat.power_w_m3[6], expected_power, 1e-11);
        assert!(
            heat.power_w_m3
                .iter()
                .enumerate()
                .all(|(cell, &p)| p >= 0.0 && (cell == 6 || p < 1e-20))
        );
        near(heat.viscous_heat_j, expected_power * volume * DT, 1e-18);
        near(heat.authored_heat_j, 12.0 * volume * DT, 1e-18);
        near(heat.source_partition_defect_j, 0.0, 1e-18);
        near(
            balance.source_j,
            (12.0 + expected_power) * volume * DT,
            1e-18,
        );
        near(balance.source_j, reference.source_j, 1e-18);
        near(balance.stored_change_j, reference.stored_change_j, 1e-12);
        near(balance.bottom_wall_j, reference.bottom_wall_j, 1e-12);
        near(balance.top_wall_j, reference.top_wall_j, 1e-12);
    };
    check(
        &heated.viscous,
        &heated.transport.heat,
        &expected.heat,
        DX * DX * DX,
    );
    assert_eq!(hot.lattice.grid.f, manual.lattice.grid.f);
    for (a, b) in hot
        .lattice
        .gpop
        .iter()
        .flatten()
        .zip(manual.lattice.gpop.iter().flatten())
    {
        near(*a, *b, 2e-15);
    }

    let mut hot = slab2();
    let mut manual = slab2();
    set_shear2(&mut hot, gamma, [force[0], force[1]]);
    set_shear2(&mut manual, gamma, [force[0], force[1]]);
    let depth = 2.0 * DX;
    let heated = hot
        .step_with_material_rheology_and_viscous_heating(
            &authored,
            &[[0.0; 2]; 16],
            depth,
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    let expected = manual
        .step_with_material_rheology(
            &prescribed,
            &[[0.0; 2]; 16],
            depth,
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    check(
        &heated.viscous,
        &heated.transport.heat,
        &expected.heat,
        DX * DX * depth,
    );
    assert_eq!(hot.lattice.grid.f, manual.lattice.grid.f);
    for (a, b) in hot
        .lattice
        .gpop
        .iter()
        .flatten()
        .zip(manual.lattice.gpop.iter().flatten())
    {
        near(*a, *b, 2e-15);
    }
}

#[test]
fn zero_stress_adds_no_heat_and_preserves_both_existing_paths() {
    // rho=9/8 makes every zero-velocity quadrature population exactly binary
    // representable: this zero-stress oracle excludes equilibrium roundoff.
    let gate = CancelGate::new();
    let mut hot = slab3();
    let mut cold = slab3();
    for state in [&mut hot, &mut cold] {
        state
            .lattice
            .grid
            .f
            .fill(fs_lbm::d3q19::equilibrium3(1.125, [0.0; 3]));
    }
    let heated = hot
        .step_with_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    let expected = cold
        .step_with_material_rheology(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    assert!(heated.viscous.power_w_m3.iter().all(|&p| p == 0.0));
    assert_eq!(heated.viscous.viscous_heat_j, 0.0);
    assert_eq!(heated.transport, expected);
    assert_eq!(hot.lattice.grid.f, cold.lattice.grid.f);
    assert_eq!(hot.lattice.gpop, cold.lattice.gpop);

    let mut hot = slab2();
    let mut cold = slab2();
    for state in [&mut hot, &mut cold] {
        state
            .lattice
            .grid
            .f
            .fill(fs_lbm::equilibrium(1.125, 0.0, 0.0));
    }
    let heated = hot
        .step_with_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.0; 2]; 16],
            DX,
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    let expected = cold
        .step_with_material_rheology(
            &[0.0; 16],
            &[[0.0; 2]; 16],
            DX,
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    assert!(heated.viscous.power_w_m3.iter().all(|&p| p == 0.0));
    assert_eq!(heated.viscous.viscous_heat_j, 0.0);
    assert_eq!(heated.transport, expected);
    assert_eq!(hot.lattice.grid.f, cold.lattice.grid.f);
    assert_eq!(hot.lattice.gpop, cold.lattice.gpop);
}

#[test]
fn g3_source_temperature_curve_closes_viscous_heat_back_into_actual_flow() {
    for slope in [0.0, 0.02] {
        // Explicitly synthetic low caloric storage makes the feedback readily
        // measurable; conductivity scales with Cp to keep the thermal stencil.
        let source = card_curve(
            "T",
            Temperature::DIMS,
            vec![(290.0, 1.0 - 10.0 * slope), (330.0, 1.0 + 30.0 * slope)],
            0.001,
        );
        let reference = resolve_liquid_state(
            &source,
            &query(305.0, 0.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let make = || {
            plan_scaling_from_material(&reference, DX, DT, 0.1)
                .unwrap()
                .thermal_slab3(MaterialThermalSlab3 {
                    nx: 2,
                    ny: 1,
                    nz: 6,
                    bottom_temperature_k: 310.0,
                    top_temperature_k: 300.0,
                    gravity_m_s2: 0.0,
                    expansion_per_k: 0.0,
                    max_cells: 16,
                })
                .unwrap()
        };
        let mut hot = make();
        let mut cold = make();
        set_shear3(&mut hot, 0.02, [0.0; 3]);
        set_shear3(&mut cold, 0.02, [0.0; 3]);
        let gate = CancelGate::new();
        let heated = hot
            .step_with_material_rheology_and_viscous_heating(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                &gate,
                |point, _| resolve(&source, &reference, point),
            )
            .unwrap();
        let unheated = cold
            .step_with_material_rheology(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                &gate,
                |point, _| resolve(&source, &reference, point),
            )
            .unwrap();
        assert_eq!(
            hot.lattice.grid.f, cold.lattice.grid.f,
            "initial-state coefficients give the same first momentum step without buoyancy"
        );
        let warm = heated
            .transport
            .endpoint
            .iter()
            .find(|s| s.cell == 6)
            .unwrap();
        let cool = unheated.endpoint.iter().find(|s| s.cell == 6).unwrap();
        assert!(warm.temperature_k > cool.temperature_k + 1e-6);
        near(
            warm.transport.dynamic_viscosity_pa_s() - cool.transport.dynamic_viscosity_pa_s(),
            slope * (warm.temperature_k - cool.temperature_k),
            1e-14,
        );
        for sample in heated
            .transport
            .applied
            .iter()
            .chain(&heated.transport.endpoint)
        {
            assert_eq!(
                sample.material.source().material().card_identity(),
                source.content_hash()
            );
            assert_eq!(
                sample.material.source().material().query(),
                query(
                    sample.temperature_k,
                    sample.material.point().shear_rate_per_s
                )
            );
        }
        hot.step_with_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |point, _| resolve(&source, &reference, point),
        )
        .unwrap();
        cold.step_with_material_rheology(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &gate,
            |point, _| resolve(&source, &reference, point),
        )
        .unwrap();
        if slope == 0.0 {
            assert_eq!(hot.lattice.grid.f, cold.lattice.grid.f);
        } else {
            assert_ne!(
                hot.lattice.grid.f, cold.lattice.grid.f,
                "the changed material temperature must change the next momentum step"
            );
        }
    }
}

#[test]
fn viscous_overheat_endpoint_refusal_and_cancellation_do_not_publish() {
    let gate = CancelGate::new();
    let mut slab = slab3();
    set_shear3(&mut slab, 0.02, [0.0; 3]);
    let wall = slab.lattice.gpop[0];
    slab.lattice.gpop[2..14].fill(wall);
    let (f, g, tau, forces) = (
        slab.lattice.grid.f.clone(),
        slab.lattice.gpop.clone(),
        slab.lattice.grid.tau.clone(),
        slab.lattice.grid.fext.clone(),
    );
    let result = slab.step_with_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        HeatSourceLimits {
            temperature_range_k: [290.0, 310.0],
            ..limits()
        },
        controls(),
        &gate,
        |_, _| Ok(fluid(1.0)),
    );
    assert!(matches!(
        result,
        Err(MaterialTransportError::Thermal(
            fs_lbm::thermal::sources::ThermalSourceError::Temperature { .. }
        ))
    ));
    let mut calls = 0;
    let result = slab.step_with_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &gate,
        |_, _| {
            calls += 1;
            if calls == 13 {
                Err(MaterialError::State {
                    what: "heated endpoint outside source support".into(),
                })
            } else {
                Ok(fluid(1.0))
            }
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::State { .. })));
    let result = slab.step_with_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &gate,
        |point, _| {
            if point.cell == 6 {
                gate.request();
            }
            Ok(fluid(1.0))
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::Cancelled)));
    assert_eq!(slab.lattice.grid.f, f);
    assert_eq!(slab.lattice.gpop, g);
    assert_eq!(slab.lattice.grid.tau, tau);
    assert_eq!(slab.lattice.grid.fext, forces);
}

#[test]
fn unrepresentable_viscous_power_refuses_before_mutating_the_slab() {
    // A representable lattice/thermal map can still have an unrepresentable
    // dimensional stress power. This is a conversion boundary, not a material.
    let dt = 1e-160;
    let dx = (0.01_f64 * dt).sqrt();
    let mut slab = plan_scaling_from_material(&fluid(1.0), dx, dt, 0.0)
        .unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2,
            ny: 1,
            nz: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 16,
        })
        .unwrap();
    set_shear3(&mut slab, 0.02, [0.0; 3]);
    let (f, g, tau, forces) = (
        slab.lattice.grid.f.clone(),
        slab.lattice.gpop.clone(),
        slab.lattice.grid.tau.clone(),
        slab.lattice.grid.fext.clone(),
    );
    let result = slab.step_with_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &CancelGate::new(),
        |_, _| Ok(fluid(1.0)),
    );
    assert!(matches!(
        result,
        Err(MaterialTransportError::Invalid(
            "physical viscous heating is unrepresentable"
        ))
    ));
    assert_eq!(slab.lattice.grid.f, f);
    assert_eq!(slab.lattice.gpop, g);
    assert_eq!(slab.lattice.grid.tau, tau);
    assert_eq!(slab.lattice.grid.fext, forces);
}

#[test]
fn g1_combined_buoyancy_and_viscous_heat_match_independent_force_and_power_in_both_dimensions() {
    let gamma = 0.02;
    let expected_power = (gamma / DT) * (gamma / DT); // mu = 1 Pa s
    let material = FluidTransportState::try_new(990.0, 1000.0, 1000.0, 1.0).unwrap();
    let mut prescribed = [0.0; 16];
    prescribed[6] = expected_power;
    let gate = CancelGate::new();
    let config = BoussinesqConfig {
        gravity_m_s2: [1.0, 2.0, -9.0],
        max_relative_density_departure: 0.02,
        max_relative_specific_heat_departure: 0.0,
    };
    let mut coupled = slab3();
    let mut manual = slab3();
    set_shear3(&mut coupled, gamma, [0.0; 3]);
    set_shear3(&mut manual, gamma, [0.0; 3]);
    // An independent background buoyancy must be replaced, not counted twice.
    coupled.lattice.gbeta = 0.07;
    let report = coupled
        .step_with_boussinesq_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.5, -0.25, 0.75]; 16],
            limits(),
            controls(),
            config,
            &gate,
            |_, _| Ok(material),
        )
        .unwrap();
    // Independently calculate F = F_external + (990 - 1000) * gravity.
    let expected = manual
        .step_with_material_rheology(
            &prescribed,
            &[[-9.5, -20.25, 90.75]; 16],
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    assert_eq!(coupled.lattice.grid.f, manual.lattice.grid.f);
    assert_eq!(coupled.lattice.grid.fext, manual.lattice.grid.fext);
    for (a, b) in coupled
        .lattice
        .gpop
        .iter()
        .flatten()
        .zip(manual.lattice.gpop.iter().flatten())
    {
        near(*a, *b, 2e-15);
    }
    near(
        report.viscous.viscous_heat_j,
        expected_power * DX * DX * DX * DT,
        1e-18,
    );
    near(
        report.transport.heat.stored_change_j,
        expected.heat.stored_change_j,
        1e-12,
    );
    near(
        report.transport.heat.source_j,
        expected.heat.source_j,
        1e-18,
    );
    assert_eq!(coupled.lattice.gbeta, 0.07);

    let config = BoussinesqConfig {
        gravity_m_s2: [1.0, -9.0],
        max_relative_density_departure: 0.02,
        max_relative_specific_heat_departure: 0.0,
    };
    let mut coupled = slab2();
    let mut manual = slab2();
    set_shear2(&mut coupled, gamma, [0.0; 2]);
    set_shear2(&mut manual, gamma, [0.0; 2]);
    coupled.lattice.gbeta = 0.07;
    let depth = 2.0 * DX;
    let report = coupled
        .step_with_boussinesq_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.5, -0.25]; 16],
            depth,
            limits(),
            controls(),
            config,
            &gate,
            |_, _| Ok(material),
        )
        .unwrap();
    let expected = manual
        .step_with_material_rheology(
            &prescribed,
            &[[-9.5, 89.75]; 16],
            depth,
            limits(),
            controls(),
            &gate,
            |_, _| Ok(fluid(1.0)),
        )
        .unwrap();
    assert_eq!(coupled.lattice.grid.f, manual.lattice.grid.f);
    assert_eq!(coupled.lattice.grid.fext, manual.lattice.grid.fext);
    for (a, b) in coupled
        .lattice
        .gpop
        .iter()
        .flatten()
        .zip(manual.lattice.gpop.iter().flatten())
    {
        near(*a, *b, 2e-15);
    }
    near(
        report.viscous.viscous_heat_j,
        expected_power * DX * DX * depth * DT,
        1e-18,
    );
    near(
        report.transport.heat.stored_change_j,
        expected.heat.stored_change_j,
        1e-12,
    );
    near(
        report.transport.heat.source_j,
        expected.heat.source_j,
        1e-18,
    );
    assert_eq!(coupled.lattice.gbeta, 0.07);
}

fn card_with_density(density_slope: f64, shear_slope: f64) -> MaterialCard {
    let original = card_surface(0.02, shear_slope, 0.001);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (290.0, 1000.0 + 10.0 * density_slope),
                    (330.0, 1000.0 - 30.0 * density_slope),
                ],
                dims: Dims([-3, 1, 0, 0, 0, 0]),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

#[test]
fn g3_one_material_card_drives_buoyancy_rheology_and_heat_without_rebinding_mass() {
    let config = BoussinesqConfig {
        gravity_m_s2: [0.0, 0.0, -9.0],
        max_relative_density_departure: 0.03,
        max_relative_specific_heat_departure: 0.0,
    };
    let mut flows = Vec::new();
    let mut heats = Vec::new();
    for (density_slope, shear_slope) in [(0.0, 0.0), (0.5, 0.0), (0.5, 0.2)] {
        let source = card_with_density(density_slope, shear_slope);
        let reference = resolve_liquid_state(
            &source,
            &query(300.0, 0.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut slab = plan_scaling_from_material(&reference, DX, DT, 0.1)
            .unwrap()
            .thermal_slab3(MaterialThermalSlab3 {
                nx: 2,
                ny: 1,
                nz: 6,
                bottom_temperature_k: 310.0,
                top_temperature_k: 300.0,
                gravity_m_s2: 0.0,
                expansion_per_k: 0.0,
                max_cells: 16,
            })
            .unwrap();
        set_shear3(&mut slab, 0.02, [0.0; 3]);
        let report = slab
            .step_with_boussinesq_material_rheology_and_viscous_heating(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                config,
                &CancelGate::new(),
                |point, _| resolve(&source, &reference, point),
            )
            .unwrap();
        for sample in report
            .transport
            .applied
            .iter()
            .chain(&report.transport.endpoint)
        {
            let point = sample.material.point();
            let original = sample.material.source().material();
            assert_eq!(original.material().card_identity(), source.content_hash());
            assert_eq!(
                original.material().query(),
                query(point.temperature_k, point.shear_rate_per_s)
            );
            near(
                original.density_kg_m3(),
                1000.0 - density_slope * (point.temperature_k - 300.0),
                2e-12,
            );
            near(
                sample.transport.dynamic_viscosity_pa_s(),
                1.0 + 0.02 * (point.temperature_k - 300.0) + shear_slope * point.shear_rate_per_s,
                2e-12,
            );
            // Inertia and sensible storage still belong to the same reference.
            assert_eq!(sample.transport.density_kg_m3(), 1000.0);
            assert_eq!(sample.transport.specific_heat_j_kg_k(), 1.0);
        }
        let mut expected_heat_j = 0.0;
        for sample in &report.transport.applied {
            let original = sample.material.source().material();
            let rate = sample.material.point().shear_rate_per_s;
            let power = original.dynamic_viscosity_pa_s() * rate * rate;
            near(report.viscous.power_w_m3[sample.cell], power, 2e-12);
            expected_heat_j += power * DX * DX * DX * DT;
            near(
                slab.lattice.grid.fext[sample.cell][2],
                (original.density_kg_m3() - 1000.0) * -9.0 / 1000.0 * DT * DT / DX,
                1e-16,
            );
        }
        near(report.viscous.viscous_heat_j, expected_heat_j, 1e-18);
        near(report.viscous.source_partition_defect_j, 0.0, 1e-18);
        near(report.transport.heat.source_j, expected_heat_j, 1e-18);
        // Repeated use samples the actually accepted heated state, rather than
        // restarting temperature or copying the initial material snapshot.
        let next = slab
            .step_with_boussinesq_material_rheology_and_viscous_heating(
                &[0.0; 16],
                &[[0.0; 3]; 16],
                limits(),
                controls(),
                config,
                &CancelGate::new(),
                |point, _| resolve(&source, &reference, point),
            )
            .unwrap();
        for (previous, next) in report
            .transport
            .endpoint
            .iter()
            .zip(&next.transport.applied)
        {
            assert_eq!(previous.temperature_k, next.temperature_k);
            assert_eq!(previous.material.source(), next.material.source());
        }
        flows.push(slab.lattice.grid.f);
        heats.push(report.viscous.viscous_heat_j);
    }
    assert_ne!(
        flows[0], flows[1],
        "density data must alter the actual evolving flow"
    );
    assert_ne!(
        flows[1], flows[2],
        "shear-response data must alter the same flow"
    );
    assert_ne!(
        heats[1], heats[2],
        "the same viscosity change must alter physical heat"
    );
}

#[test]
fn combined_material_refusals_and_cancellation_preserve_all_accepted_state() {
    let config = BoussinesqConfig {
        gravity_m_s2: [0.0, 0.0, -9.0],
        max_relative_density_departure: 0.02,
        max_relative_specific_heat_departure: 0.0,
    };
    let mut slab = slab3();
    set_shear3(&mut slab, 0.02, [0.0; 3]);
    slab.lattice.gbeta = 0.17;
    let (f, g, tau, force) = (
        slab.lattice.grid.f.clone(),
        slab.lattice.gpop.clone(),
        slab.lattice.grid.tau.clone(),
        slab.lattice.grid.fext.clone(),
    );
    for endpoint in [false, true] {
        let mut calls = 0;
        let result = slab.step_with_boussinesq_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            config,
            &CancelGate::new(),
            |_, _| {
                calls += 1;
                let rho = if !endpoint || calls > 12 {
                    900.0
                } else {
                    990.0
                };
                FluidTransportState::try_new(rho, 1000.0, 1000.0, 1.0)
            },
        );
        assert!(matches!(result, Err(MaterialTransportError::State { .. })));
        if endpoint {
            assert!(calls > 12);
        }
        assert_eq!(slab.lattice.grid.f, f);
        assert_eq!(slab.lattice.gpop, g);
        assert_eq!(slab.lattice.grid.tau, tau);
        assert_eq!(slab.lattice.grid.fext, force);
        assert_eq!(slab.lattice.gbeta, 0.17);
    }
    let gate = CancelGate::new();
    let result = slab.step_with_boussinesq_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        config,
        &gate,
        |point, _| {
            if point.cell == 6 {
                gate.request();
            }
            FluidTransportState::try_new(990.0, 1000.0, 1000.0, 1.0)
        },
    );
    assert_eq!(result, Err(MaterialTransportError::Cancelled));
    assert_eq!(slab.lattice.grid.f, f);
    assert_eq!(slab.lattice.gpop, g);
    assert_eq!(slab.lattice.grid.tau, tau);
    assert_eq!(slab.lattice.grid.fext, force);
    assert_eq!(slab.lattice.gbeta, 0.17);
}
