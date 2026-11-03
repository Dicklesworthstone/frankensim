//! Two original fusion sources, separate enthalpy datums and physical contact.
#[path = "heterogeneous/history.rs"]
mod history;
use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::{
    EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler as Engine, HeterogeneousEnthalpyError,
};
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA as G};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RDIM,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RKEY, ConductionProblem, InterfaceResistance,
    InterfaceSurface, ThermalBoundary, ThermalInterfaces,
};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};

const INITIAL: [f64; 8] = [60.0, 180.0, 240.0, 420.0, 1060.0, 1220.0, 1380.0, 1660.0];
const LOAD: [f64; 8] = [8.0, -3.0, 12.0, 5.0, -2.0, 7.0, 1.0, 4.0];
const WEIGHT: [f64; 8] = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, -0.2, 0.8];
const DIRECTION: [f64; 8] = [0.4, -0.7, 0.2, 0.6, -0.5, 0.3, 0.8, -0.1];

fn second_source() -> MaterialMeltingTransport {
    // Second card: Tm=400 K, L=400 J/kg, rho_s=3, rho_l=1.5,
    // twice the conductivity and a disjoint h datum offset of 1000 J/kg.
    let base = card(true, false);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        let name = claim.key.name();
        let factor = if name == "rho" {
            1.5
        } else if name == "k" || name == "fusion" {
            2.0
        } else {
            1.0
        };
        match &mut claim.value {
            PropertyValue::Curve { knots, .. } => {
                for (t, value) in knots {
                    *t += 50.0;
                    *value *= factor;
                }
            }
            PropertyValue::Scalar { value, .. } => {
                *value = if name == "Tm" {
                    *value + 50.0
                } else {
                    *value * factor
                };
            }
            _ => unreachable!("synthetic scalar/curve card"),
        }
        if let Some((lo, hi)) = claim.validity.bounds().get("T").copied() {
            claim.validity =
                claim
                    .validity
                    .with_quantity("T", absolute_temperature(), lo + 50.0, hi + 50.0);
        }
        claims.insert_claim(claim).unwrap();
    }
    let mut id = base.id().clone();
    id.chemistry = "second synthetic fusion source".into();
    let card = MaterialCard::assemble(id, claims, vec![]).unwrap();
    let cp = PropertyKey::new("cp", CP);
    let rho = PropertyKey::new("rho", RHO);
    let k = PropertyKey::new("k", K);
    let points = [
        point(350.0, Some(0.0), 100_000.0),
        point(400.0, Some(0.0), 100_000.0),
        point(400.0, Some(1.0), 100_000.0),
        point(500.0, Some(1.0), 100_000.0),
        point(400.0, None, 100_000.0),
    ];
    let branch = |i| HeatCapacitySourceQuery {
        heat_capacity: &cp,
        density: &rho,
        temperature_axis: "T",
        lower: &points[i],
        upper: &points[i + 1],
        selection: SelectionPolicy::SingleClaimOnly,
    };
    MaterialMeltingTransport::from_isotropic_card(
        &card,
        MaterialMeltingSource {
            caloric: MeltingSourceQuery {
                solid: branch(0),
                liquid: branch(2),
                melting_point: &PropertyKey::with_quantity("Tm", absolute_temperature()),
                latent_heat_fusion: &PropertyKey::new("fusion", ENERGY),
                transition: &points[4],
                transition_selection: SelectionPolicy::SingleClaimOnly,
                phase_selector: Some("phase"),
            },
            conductivity: [&k, &k],
            reference_phase: SolidLiquidPhase::Solid,
            reference: &points[0],
            lower_specific_enthalpy_j_kg: 1000.0,
            maximum_interpolation_error_j_kg: 1e-6,
            coexistence: CoexistenceConductivity::Series,
        },
    )
    .unwrap()
}
fn coupon() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., -1.],
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
        ],
    )
    .unwrap()
}
fn allowance() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 8,
        max_elements: 2,
    }
}
fn policy() -> EnthalpyStepConfig {
    let mut c = config();
    c.newton.linear_restart = 8;
    c
}
fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 64,
        restart: 8,
    }
}
fn contact(mesh: &ConductionMesh, boundary: &ThermalBoundary) -> ThermalInterfaces {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RKEY, RDIM),
            value: PropertyValue::Scalar {
                value: 2.0,
                dims: RDIM,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured coupon contact".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid/liquid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        texture_frame: "reference".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("first"),
        surface("second"),
        SystemContext {
            medium: "declared finite interface".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let resistance = InterfaceResistance::from_card(
        "joint",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    let mut pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
    for pair in &mut pairs {
        if mesh.boundary()[pair.side_a].element != 0 {
            std::mem::swap(&mut pair.side_a, &mut pair.side_b);
        }
    }
    ThermalInterfaces::new(
        mesh,
        boundary,
        vec![InterfaceSurface::new("joint", pairs, resistance).unwrap()],
    )
    .unwrap()
}
fn state(h: f64, material: usize) -> (f64, f64) {
    if material == 0 {
        let s = reference(h, CoexistenceConductivity::Parallel);
        return (s.0, s.2);
    }
    let h = h - 1000.0;
    if h < 100.0 {
        (350.0 + h / 2.0, 4.0 + h / 25.0)
    } else if h < 500.0 {
        (400.0, 8.0 * (1.0 + (h - 100.0) / 400.0))
    } else {
        (400.0 + (h - 500.0) / 4.0, 16.0 + (h - 500.0) / 50.0)
    }
}
fn balance(h: &[f64], old: &[f64], dt: f64) -> [f64; 8] {
    let t: [f64; 8] = std::array::from_fn(|i| state(h[i], i / 4).0);
    let mut result = [0.0; 8];
    for material in 0..2 {
        let offset = 4 * material;
        let k = (0..4)
            .map(|i| state(h[offset + i], material).1 / 4.0)
            .sum::<f64>();
        let local = &t[offset..offset + 4];
        let laplace = [
            3.0 * local[0] - local[1] - local[2] - local[3],
            local[1] - local[0],
            local[2] - local[0],
            local[3] - local[0],
        ];
        for i in 0..4 {
            result[offset + i] = (2.0 + material as f64) / 24.0 * (h[offset + i] - old[offset + i])
                + dt * (k * laplace[i] / 6.0
                    - (LOAD[offset..offset + 4].iter().sum::<f64>() + LOAD[offset + i]) / 120.0);
        }
    }
    let difference = [t[0] - t[4], t[1] - t[5], t[2] - t[6]];
    for i in 0..3 {
        let q = (difference.iter().sum::<f64>() + difference[i]) / 48.0;
        result[i] += dt * q;
        result[i + 4] -= dt * q;
    }
    result
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn original_regional_fusion_sources_drive_flux_balance_tangent_and_both_integrators() {
    let mesh = coupon();
    let materials = [model(CoexistenceConductivity::Parallel), second_source()];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let contact = contact(&mesh, &boundary);
    let source = ScalarField::nodal("heating", 8, LOAD.to_vec()).unwrap();
    with_cx(|cx| {
        let engine =
            Engine::from_melting_materials(cx, &mesh, &materials, &[0, 1], allowance()).unwrap();
        let problem = engine
            .melting_transport_problem(&boundary, &source)
            .unwrap();
        assert_eq!(
            materials[0].caloric().plateau_enthalpy_bounds_j_kg(),
            [100.0, 300.0]
        );
        assert_eq!(
            materials[1].caloric().plateau_enthalpy_bounds_j_kg(),
            [1100.0, 1500.0]
        );
        assert!(
            materials[0]
                .conductivity_at_specific_enthalpy(INITIAL[4])
                .is_err()
        );
        assert!(
            materials[1]
                .conductivity_at_specific_enthalpy(INITIAL[0])
                .is_err()
        );
        for i in 0..8 {
            assert!(std::ptr::eq(
                engine.phase_curve_at(i).unwrap(),
                materials[i / 4].caloric().curve()
            ));
            near(
                engine.reference_nodal_masses_kg()[i],
                (2.0 + (i / 4) as f64) / 24.0,
                1e-15,
            );
        }
        assert!(std::ptr::eq(
            engine.melting_transport_source(1).unwrap(),
            &materials[1]
        ));
        assert!(engine.melting_transport_source(2).is_none());
        let response = engine
            .linearize_step(cx, problem, Some(&contact), &INITIAL, DT, policy())
            .unwrap();
        let thermal = response.transport();
        let h = &thermal.primal().specific_enthalpy_j_kg;
        for r in balance(h, &INITIAL, DT) {
            near(r, 0.0, 2e-10);
        }
        let flux = engine.melting_transport_heat_flux(cx, h).unwrap();
        for m in 0..2 {
            let t: [f64; 4] = std::array::from_fn(|i| state(h[4 * m + i], m).0);
            let k = (0..4).map(|i| state(h[4 * m + i], m).1 / 4.0).sum::<f64>();
            for axis in 0..3 {
                near(
                    flux[m][axis],
                    -k * (t[axis + 1] - t[0]) * if axis == 2 && m == 0 { -1.0 } else { 1.0 },
                    2e-10,
                );
            }
        }
        for column in 0..8 {
            let mut direction = [0.0; 8];
            direction[column] = 1.0;
            let action = thermal.apply_jacobian(cx, &direction).unwrap();
            let mut plus = h.clone();
            let mut minus = h.clone();
            plus[column] += 0.001;
            minus[column] -= 0.001;
            let plus = balance(&plus, &INITIAL, DT);
            let minus = balance(&minus, &INITIAL, DT);
            for row in 0..8 {
                near(action[row], (plus[row] - minus[row]) / 0.002, 5e-9);
            }
        }
        let inputs = [EnthalpySdirk2Input {
            problem,
            interfaces: Some(&contact),
            radiation: EnthalpyRadiation::None,
            contacts: None,
        }; 2];
        let sdirk = engine
            .linearize_sdirk2(cx, inputs, &INITIAL, DT, policy(), 0)
            .unwrap();
        let stages = sdirk.primal();
        let first = &stages.first_stage.specific_enthalpy_j_kg;
        for r in balance(first, &INITIAL, G * DT) {
            near(r, 0.0, 2e-10);
        }
        let residual = balance(
            &stages.second_stage.specific_enthalpy_j_kg,
            &INITIAL,
            G * DT,
        );
        let carry = balance(first, first, (1.0 - G) * DT);
        for i in 0..8 {
            near(residual[i] + carry[i], 0.0, 2e-10);
        }
        let gradient = sdirk.pullback(cx, &WEIGHT, linear()).unwrap();
        let evaluate = |delta| {
            let old: [f64; 8] = std::array::from_fn(|i| INITIAL[i] + delta * DIRECTION[i]);
            dot(
                &engine
                    .advance_sdirk2(cx, inputs, &old, DT, policy())
                    .unwrap()
                    .second_stage
                    .specific_enthalpy_j_kg,
                &WEIGHT,
            )
        };
        let e = 0.02;
        near(
            dot(&gradient.previous_specific_enthalpy, &DIRECTION),
            (-evaluate(2.0 * e) + 8.0 * evaluate(e) - 8.0 * evaluate(-e) + evaluate(-2.0 * e))
                / (12.0 * e),
            2e-8,
        );
        let other = ConductivityModel::isotropic_declared(1.0).unwrap();
        assert!(
            engine
                .advance(
                    cx,
                    ConductionProblem {
                        material: &other,
                        ..problem
                    },
                    Some(&contact),
                    &INITIAL,
                    DT,
                    policy()
                )
                .is_err()
        );
    });
}

#[test]
fn unequal_fusion_plateaus_exchange_internal_heat_with_independent_reference_inventory() {
    let mesh = coupon();
    let materials = [model(CoexistenceConductivity::Parallel), second_source()];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let contact = contact(&mesh, &boundary);
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let density = [1.0, 2.0, 3.0, 4.0, 3.0, 4.0, 5.0, 6.0];
        let engine = Engine::from_melting_materials_with_nodal_reference_density(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            &density,
            allowance(),
        )
        .unwrap();
        let original_masses = engine.reference_nodal_masses_kg().to_vec();
        let problem = engine
            .melting_transport_problem(&boundary, &source)
            .unwrap();
        let old = [200.0, 200.0, 200.0, 200.0, 1300.0, 1300.0, 1300.0, 1300.0];
        let masses: [f64; 8] = std::array::from_fn(|i| {
            (density[4 * (i / 4)..4 * (i / 4) + 4].iter().sum::<f64>() + density[i]) / 120.0
        });
        for i in 0..8 {
            near(engine.reference_nodal_masses_kg()[i], masses[i], 1e-15);
        }
        assert!(engine.reference_material(0).is_none());
        assert!(engine.melting_transport_source(1).is_some());
        let input = EnthalpySdirk2Input {
            problem,
            interfaces: Some(&contact),
            radiation: EnthalpyRadiation::None,
            contacts: None,
        };
        let step = engine
            .advance_sdirk2(cx, [input; 2], &old, 0.1, policy())
            .unwrap();
        // R''=2 m2K/W, area=1/2: Q=12.5 W from the 400 K plateau
        // into the 350 K plateau, independently of either latent fraction.
        for (stage, weight) in [(&step.first_stage, G), (&step.second_stage, 1.0)] {
            for i in 0..8 {
                let q = if i % 4 == 3 {
                    0.0
                } else if i < 4 {
                    12.5 / 3.0
                } else {
                    -12.5 / 3.0
                };
                near(
                    stage.specific_enthalpy_j_kg[i],
                    old[i] + 0.1 * weight * q / masses[i],
                    2e-9,
                );
                near(stage.temperature[i], if i < 4 { 350.0 } else { 400.0 }, 0.0);
                let (start, latent) = if i < 4 {
                    (100.0, 200.0)
                } else {
                    (1100.0, 400.0)
                };
                near(
                    stage.liquid_mass_fraction[i],
                    (stage.specific_enthalpy_j_kg[i] - start) / latent,
                    1e-14,
                );
            }
        }
        near(step.stored_energy_change_j, 0.0, 2e-10);
        near(step.source_j, 0.0, 0.0);
        near(step.energy_residual_j, 0.0, 2e-10);
        for flux in engine
            .melting_transport_heat_flux(cx, &step.second_stage.specific_enthalpy_j_kg)
            .unwrap()
        {
            for value in flux {
                near(value, 0.0, 1e-12);
            }
        }
        assert_eq!(engine.reference_nodal_masses_kg(), original_masses);
    });
}

#[test]
fn assigned_melting_retains_uniform_results_and_refuses_invalid_material_ownership() {
    let mesh = mesh();
    let material = [model(CoexistenceConductivity::Parallel)];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let source = ScalarField::nodal("heating", 4, Q.to_vec()).unwrap();
    with_cx(|cx| {
        let uniform = material[0].prepare_uniform(cx, &mesh, budget()).unwrap();
        let assigned =
            Engine::from_melting_materials(cx, &mesh, &material, &[0], budget()).unwrap();
        let actual = assigned
            .advance(
                cx,
                assigned
                    .melting_transport_problem(&boundary, &source)
                    .unwrap(),
                None,
                &OLD,
                DT,
                config(),
            )
            .unwrap();
        let expected = uniform
            .advance(
                cx,
                material[0].problem(&mesh, &boundary, &source),
                None,
                &OLD,
                DT,
                config(),
            )
            .unwrap();
        assert_eq!(
            actual.specific_enthalpy_j_kg,
            expected.specific_enthalpy_j_kg
        );
        assert_eq!(
            assigned.melting_transport_heat_flux(cx, &OLD).unwrap(),
            material[0].element_heat_flux(cx, &mesh, &OLD).unwrap()
        );
        assert!(assigned.melting_transport_heat_flux(cx, &[0.0; 3]).is_err());
        assert!(
            assigned
                .melting_transport_heat_flux(cx, &[-1.0; 4])
                .is_err()
        );
        assert!(Engine::from_melting_materials(cx, &mesh, &material, &[1], budget()).is_err());
        assert!(Engine::from_melting_materials(cx, &mesh, &material, &[], budget()).is_err());
        assert!(Engine::from_melting_materials(cx, &mesh, &[], &[0], budget()).is_err());
        assert!(
            assigned
                .with_nodal_reference_density(cx, &[2.0; 4])
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                Engine::from_melting_materials(cancelled, &mesh, &material, &[0], budget())
                    .is_err()
            );
            assert!(
                material[0]
                    .element_heat_flux(cancelled, &mesh, &OLD)
                    .is_err()
            );
        });
        let shared = ConductionMesh::new(
            TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]),
            vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [0., 0., -1.],
            ],
        )
        .unwrap();
        let materials = [model(CoexistenceConductivity::Parallel), second_source()];
        assert!(matches!(
            Engine::from_melting_materials(cx, &shared, &materials, &[0, 1], allowance()),
            Err(HeterogeneousEnthalpyError::SharedVertex { .. })
        ));
    });
}
