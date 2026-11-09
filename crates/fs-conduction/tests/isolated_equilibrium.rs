//! G1/G3 isolated endpoints from independent caloric and heat-balance equations.
mod support;

use fs_blake3::ContentHash;
use fs_conduction::material::{ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::transient::enthalpy::equilibrium::{
    IsolatedEquilibriumConfig, IsolatedEquilibriumError,
};
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, SurfaceHeatStorage,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    InterfaceResistance, InterfaceSurface, ScalarField, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder, ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::thermal::transport::storage::TabulatedHeatStorage;
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}, allowance {tolerance:e}");
}

fn tet() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ).unwrap()
}

fn pair(touching: bool) -> ConductionMesh {
    let shift = if touching { 0.0 } else { 3.0 };
    ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 6, 5, 7]]),
        vec![
            [0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [shift, 0.0, 0.0], [shift + 1.0, 0.0, 0.0],
            [shift, 1.0, 0.0], [shift, 0.0, -1.0],
        ],
    ).unwrap()
}

fn adiabatic(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh).adiabatic_remainder().finish().unwrap()
}
fn budget(mesh: &ConductionMesh) -> EnthalpyBudget {
    EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() }
}
fn config() -> IsolatedEquilibriumConfig {
    IsolatedEquilibriumConfig {
        max_temperature_evaluations: 96,
        energy_tolerance_j: 1e-8,
        temperature_tolerance_k: 1e-9,
        heat_residual_tolerance_w: 1e-9,
        entropy_tolerance_j_k: 1e-10,
    }
}
fn step_config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 4,
        energy_tolerance_j: 1e-8,
    }
}
fn problem<'a>(
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    material: &'a ConductivityModel,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem { mesh, boundary, material, source, element_materials: None }
}
fn sensible(knots: &[(f64, f64)], datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x62; 32]), SolidLiquidPhase::Solid,
        knots.iter().map(|&(t, h)| EnthalpyPhaseKnot {
            temperature_k: t, specific_enthalpy_j_kg: h + datum,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 37.0,
        }).collect(),
    ).unwrap()
}
fn latent() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x63; 32]),
        [(300.0, 0.0, 0.0), (347.0, 100.0, 0.0),
         (347.0, 300.0, 1.0), (500.0, 600.0, 1.0)]
            .into_iter().map(|(t, h, f)| EnthalpyPhaseKnot {
                temperature_k: t, specific_enthalpy_j_kg: h,
                liquid_mass_fraction: f, bulk_density_kg_m3: 37.0 - 20.0 * f,
            }).collect(),
    ).unwrap()
}

fn contact(mesh: &ConductionMesh, boundary: &ThermalBoundary, value: f64) -> ThermalInterfaces {
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(), phase: "solid".into(),
            process: "G1 fixture".into(), revision: 0,
        },
        texture_frame: "declared isotropic".into(),
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),
        value: PropertyValue::Scalar { value, dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: "G1 isolated calorimetric contact fixture".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let card = InterfaceSystemCard::assemble(surface("a"), surface("b"), SystemContext {
        medium: "vacuum".into(), third_body: Some("declared contact".into()),
        environment: "closed calorimeter".into(), history: "fixed".into(),
    }, claims, Vec::new()).unwrap();
    let resistance = InterfaceResistance::from_card(
        "interface", &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly,
    ).unwrap();
    let faces = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
    assert_eq!(faces.len(), 1);
    ThermalInterfaces::new(mesh, boundary, vec![
        InterfaceSurface::new("interface", faces, resistance).unwrap(),
    ]).unwrap()
}

#[test]
fn g1_nonlinear_caloric_equilibrium_preserves_graded_reference_mass_and_holds_in_native_transport() {
    with_cx(|cx| {
        let mesh = tet();
        let boundary = adiabatic(&mesh);
        let law = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(300.0, 0.3), (450.0, 2.0)],
        ).unwrap());
        let source = ScalarField::Uniform(0.0);
        // Cp=2 below 350 K and Cp=4 above; density is the independent P1
        // inventory [1,2,3,4], not the chart's unrelated phase-state density.
        let curve = sensible(&[(300.0, 0.0), (350.0, 100.0), (450.0, 500.0)], 0.0);
        let storage = EnthalpyBackwardEuler::from_nodal_reference_density(
            cx, &mesh, &curve, &[1.0, 2.0, 3.0, 4.0], budget(&mesh),
        ).unwrap();
        let old = [0.0, 40.0, 300.0, 500.0];
        let p = problem(&mesh, &boundary, &law, &source);
        let family = storage.isolated_equilibrium(cx, p, None, &old, config()).unwrap();
        assert_eq!(family.component_of_vertex(), &[0; 4]);
        assert!(!family.requires_allocation());
        let state = family.unique_state(cx).unwrap();
        // Exact P1 masses are [11,12,13,14]/120 kg; h_eq=227.6 J/kg.
        let expected_t = 381.9;
        for (i, (&h, &t)) in state.specific_enthalpy_j_kg.iter().zip(&state.temperature).enumerate() {
            close(storage.reference_nodal_masses_kg()[i], (11 + i) as f64 / 120.0, 1e-15);
            close(h, 227.6, 3e-8);
            close(t, expected_t, 1e-8);
        }
        let entropy_primitive = |t: f64| {
            if t <= 350.0 { 2.0 * (t / 300.0).ln() }
            else { 2.0 * (350.0_f64 / 300.0).ln() + 4.0 * (t / 350.0).ln() }
        };
        let entropy: f64 = [300.0, 320.0, 400.0, 450.0].into_iter().enumerate()
            .map(|(i, t)| (11 + i) as f64 / 120.0
                * (entropy_primitive(expected_t) - entropy_primitive(t))).sum();
        close(state.component_entropy_change_j_k[0], entropy, 1e-10);
        assert!(entropy > 0.0);
        close(state.component_energy_residual_j[0], 0.0, config().energy_tolerance_j);
        let held = storage.advance(cx, p, None, &state.specific_enthalpy_j_kg, 1.0, step_config()).unwrap();
        assert_eq!(held.specific_enthalpy_j_kg, state.specific_enthalpy_j_kg);
        assert_eq!(old, [0.0, 40.0, 300.0, 500.0]);
    });
}

#[test]
fn g1_original_nonlinear_film_energy_changes_equilibrium_and_survives_caloric_datum_shifts() {
    with_cx(|cx| {
        let mesh = tet();
        let face = mesh.boundary().iter().position(|f| !f.vertices.contains(&3)).unwrap();
        let vertices = mesh.boundary()[face].vertices;
        let film = TabulatedHeatStorage::try_new(
            &[[280.0, 0.6], [300.0, 1.0], [450.0, 4.0]],
        ).unwrap();
        let profile = [0.0, 1.0, 3.0];
        let fractions = profile.map(|v| (4.0 + v) / 16.0);
        let face_storage = SurfaceHeatStorage::tabulated(face, &film)
            .with_inventory_profile(profile).unwrap();
        let boundary = adiabatic(&mesh);
        let law = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old_t = [290.0, 310.0, 400.0, 430.0];
        // The independent total balance in x=T-300 is
        // .01*x^2 + 3*x = 115 + sum(face_share*(x_old+.01*x_old^2)).
        let old_film: f64 = vertices.iter().zip(fractions)
            .map(|(&v, w)| { let x = old_t[v as usize] - 300.0; w * (x + 0.01 * x * x) }).sum();
        let x = (9.0 + 0.04 * (115.0 + old_film)).sqrt();
        let expected_t = 300.0 + 2.0 * (115.0 + old_film) / (3.0 + x);
        let mut temperatures = Vec::new();
        for datum in [0.0, 1e6] {
            let curve = sensible(&[(280.0, 0.0), (450.0, 340.0)], datum);
            let storage = EnthalpyBackwardEuler::uniform(
                cx, &mesh, &curve, 6.0, budget(&mesh),
            ).unwrap().with_surface_heat_storage(cx, &[face_storage]).unwrap();
            let old = old_t.map(|t| 2.0 * (t - 280.0) + datum);
            let family = storage.isolated_equilibrium(
                cx, problem(&mesh, &boundary, &law, &source), None, &old, config(),
            ).unwrap();
            let state = family.unique_state(cx).unwrap();
            for &t in &state.temperature { close(t, expected_t, 2e-8); }
            let bulk_entropy: f64 = old_t.iter().map(|&t| 0.5 * (expected_t / t).ln()).sum();
            let film_entropy: f64 = vertices.iter().zip(fractions).map(|(&v, w)| {
                let old = old_t[v as usize];
                w * (-5.0 * (expected_t / old).ln() + 0.02 * (expected_t - old))
            }).sum();
            close(state.component_entropy_change_j_k[0], bulk_entropy + film_entropy, 2e-10);
            close(state.component_energy_residual_j[0], 0.0, config().energy_tolerance_j);
            temperatures.push(state.temperature[0]);
        }
        close(temperatures[0], temperatures[1], 2e-8);
        assert!((expected_t - 357.5).abs() > 0.1, "film inventory must affect the answer");
    });
}

#[test]
fn g1_sharp_latent_equilibrium_keeps_allocation_ambiguity_and_original_dynamic_rest() {
    with_cx(|cx| {
        let mesh = tet();
        let boundary = adiabatic(&mesh);
        let law = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = latent();
        let storage = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 6.0, budget(&mesh)).unwrap();
        let old = [0.0, 50.0, 500.0, 250.0]; // Total mass 1 kg, h_bar=200 J/kg.
        let p = problem(&mesh, &boundary, &law, &source);
        let controls = IsolatedEquilibriumConfig { max_temperature_evaluations: 3, ..config() };
        let family = storage.isolated_equilibrium(cx, p, None, &old, controls).unwrap();
        assert_eq!(family.components()[0].temperature_evaluations, 3);
        assert_eq!(family.components()[0].temperature_k, 347.0);
        close(family.components()[0].latent_energy_above_lower_j, 100.0, 1e-12);
        assert_eq!(family.specific_enthalpy_intervals_j_kg(), &[[100.0, 300.0]; 4]);
        assert!(matches!(family.unique_state(cx), Err(IsolatedEquilibriumError::AllocationRequired)));
        let first = family.materialize(cx, &[100.0, 100.0, 300.0, 300.0]).unwrap();
        let second = family.materialize(cx, &[200.0; 4]).unwrap();
        assert_eq!(first.temperature, vec![347.0; 4]);
        assert_eq!(first.liquid_mass_fraction, vec![0.0, 0.0, 1.0, 1.0]);
        assert_eq!(second.liquid_mass_fraction, vec![0.5; 4]);
        close(first.component_entropy_change_j_k[0], second.component_entropy_change_j_k[0], 1e-12);
        assert!(first.component_entropy_change_j_k[0] > 0.0);
        assert!(matches!(
            family.materialize(cx, &[100.0; 4]), Err(IsolatedEquilibriumError::EnergyBalance { .. })
        ));
        assert!(matches!(
            family.materialize(cx, &[99.0, 101.0, 300.0, 300.0]),
            Err(IsolatedEquilibriumError::AllocationOutsideInterval { .. })
        ));
        let held = storage.advance(cx, p, None, &first.specific_enthalpy_j_kg, 3.0, step_config()).unwrap();
        assert_eq!(held.specific_enthalpy_j_kg, first.specific_enthalpy_j_kg);
        assert!(matches!(
            storage.isolated_equilibrium(cx, p, None, &old,
                IsolatedEquilibriumConfig { max_temperature_evaluations: 2, ..config() }),
            Err(IsolatedEquilibriumError::WorkLimit { evaluations: 2, .. })
        ));
        // At an exact extreme, positive masses force every coordinate to that
        // extreme; there is no remaining latent allocation to invent.
        let endpoint_old = [100.0; 4];
        let endpoint = storage.isolated_equilibrium(cx, p, None, &endpoint_old, config()).unwrap();
        assert!(!endpoint.requires_allocation());
        assert_eq!(endpoint.unique_state(cx).unwrap().specific_enthalpy_j_kg, endpoint_old);
    });
}

#[test]
fn g1_disconnected_inventories_stay_separate_and_finite_contacts_establish_common_equilibrium() {
    with_cx(|cx| {
        let first = sensible(&[(300.0, 0.0), (500.0, 400.0)], 0.0);
        let second = sensible(&[(300.0, 0.0), (500.0, 800.0)], 0.0);
        let materials = [
            ReferenceEnthalpyMaterial { curve: &first, reference_density_kg_m3: 6.0 },
            ReferenceEnthalpyMaterial { curve: &second, reference_density_kg_m3: 18.0 },
        ];
        let source = ScalarField::Uniform(0.0);
        let law = ConductivityModel::isotropic_declared(1.0).unwrap();
        let old = [40.0, 40.0, 40.0, 40.0, 400.0, 400.0, 400.0, 400.0];
        for touching in [false, true] {
            let mesh = pair(touching);
            let boundary = adiabatic(&mesh);
            let storage = HeterogeneousEnthalpyBackwardEuler::new(
                cx, &mesh, &materials, &[0, 1], budget(&mesh),
            ).unwrap();
            let p = problem(&mesh, &boundary, &law, &source);
            for resistance in [0.1, 1e100] {
                let interface = touching.then(|| contact(&mesh, &boundary, resistance));
                let family = storage.isolated_equilibrium(
                    cx, p, interface.as_ref(), &old, config(),
                ).unwrap();
                let state = family.unique_state(cx).unwrap();
                assert_eq!(family.components().len(), if touching { 1 } else { 2 });
                for (v, &t) in state.temperature.iter().enumerate() {
                    // Total Cp is 2 J/K and 12 J/K, from actual reference
                    // masses. A finite resistance changes rates, not this endpoint.
                    let expected = if touching { 5440.0 / 14.0 }
                        else if v < 4 { 320.0 } else { 400.0 };
                    close(t, expected, 1e-8);
                }
                if touching {
                    let t: f64 = 5440.0 / 14.0;
                    let entropy = 2.0 * (t / 320.0).ln() + 12.0 * (t / 400.0).ln();
                    close(state.component_entropy_change_j_k[0], entropy, 1e-10);
                    assert!(entropy > 0.0);
                } else {
                    assert_eq!(family.component_of_vertex(), &[0, 0, 0, 0, 1, 1, 1, 1]);
                    for e in state.component_energy_residual_j { close(e, 0.0, 1e-9); }
                }
            }
        }
    });
}

#[test]
fn g3_disconnected_latent_energy_defects_cannot_cancel_globally() {
    with_cx(|cx| {
        let mesh = pair(false);
        let curve = latent();
        let boundary = adiabatic(&mesh);
        let law = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [200.0; 8];
        let storage = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 6.0, budget(&mesh)).unwrap();
        let family = storage.isolated_equilibrium(
            cx, problem(&mesh, &boundary, &law, &source), None, &old, config(),
        ).unwrap();
        assert_eq!(family.components().len(), 2);
        // Total change is exactly zero, but each disconnected body changes by
        // ten joules. Global conservation alone would accept impossible heat.
        let wrong = [190.0, 190.0, 190.0, 190.0, 210.0, 210.0, 210.0, 210.0];
        assert!(matches!(family.materialize(cx, &wrong),
            Err(IsolatedEquilibriumError::EnergyBalance { component: 0, .. })));
        let accepted = family.materialize(cx, &old).unwrap();
        assert_eq!(accepted.specific_enthalpy_j_kg, old);
        assert_eq!(accepted.component_entropy_change_j_k, vec![0.0; 2]);
    });
}

#[test]
fn g0_open_systems_original_source_limits_and_cancellation_refuse_without_history_changes() {
    with_cx(|cx| {
        let mesh = tet();
        let curve = sensible(&[(280.0, 0.0), (450.0, 340.0)], 0.0);
        let boundary = adiabatic(&mesh);
        let law = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [40.0; 4];
        let storage = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 6.0, budget(&mesh)).unwrap();
        let balanced = ScalarField::Nodal(vec![1.0, -1.0, 1.0, -1.0]);
        assert!(storage.isolated_equilibrium(
            cx, problem(&mesh, &boundary, &law, &balanced), None, &old, config(),
        ).is_err());
        for condition in [
            ThermalBc::neumann(1.0).unwrap(),
            ThermalBc::dirichlet(300.0).unwrap(),
            ThermalBc::robin(1.0, 300.0).unwrap(),
        ] {
            let open = ThermalBoundaryBuilder::new(&mesh).remainder("bath", condition).unwrap()
                .finish().unwrap();
            assert!(storage.isolated_equilibrium(
                cx, problem(&mesh, &open, &law, &source), None, &old, config(),
            ).is_err());
        }
        let family = storage.isolated_equilibrium(
            cx, problem(&mesh, &boundary, &law, &source), None, &old, config(),
        ).unwrap();
        with_cancelled_cx(|cancelled| {
            assert!(matches!(family.unique_state(cancelled),
                Err(IsolatedEquilibriumError::Original(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. })))));
            assert!(storage.isolated_equilibrium(
                cancelled, problem(&mesh, &boundary, &law, &source), None, &old, config(),
            ).is_err());
        });
        assert_eq!(family.unique_state(cx).unwrap().specific_enthalpy_j_kg, old);
        assert_eq!(old, [40.0; 4]);

        // A narrow film source can support the initial local field but not the
        // final shared temperature. No extrapolation expands that source.
        let face = mesh.boundary().iter().position(|f| !f.vertices.contains(&3)).unwrap();
        let narrow = TabulatedHeatStorage::try_new(&[[300.0, 1.0], [310.0, 1.2]]).unwrap();
        let coated = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 6.0, budget(&mesh)).unwrap()
            .with_surface_heat_storage(cx, &[SurfaceHeatStorage::tabulated(face, &narrow)]).unwrap();
        assert!(matches!(coated.isolated_equilibrium(
            cx, problem(&mesh, &boundary, &law, &source), None, &[40.0, 40.0, 40.0, 300.0],
            config(),
        ), Err(IsolatedEquilibriumError::EnergyOutsideDomain { .. })));
    });
}

#[test]
fn g0_final_materialization_retains_heterogeneous_conductivity_source_domains() {
    with_cx(|cx| {
        let mesh = pair(true);
        let boundary = adiabatic(&mesh);
        let interface = contact(&mesh, &boundary, 0.1);
        let curve = sensible(&[(300.0, 0.0), (500.0, 400.0)], 0.0);
        let calorics = [
            ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: 6.0 },
            ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: 18.0 },
        ];
        let source = ScalarField::Uniform(0.0);
        let fallback = ConductivityModel::isotropic_declared(1.0).unwrap();
        let models = MaterialTable::new([
            (MaterialId(1), ConductivityModel::isotropic(ConductivityTable::declared_curve(
                vec![(300.0, 1.0), (340.0, 1.0)],
            ).unwrap())),
            (MaterialId(2), ConductivityModel::isotropic(ConductivityTable::declared_curve(
                vec![(350.0, 1.0), (500.0, 1.0)],
            ).unwrap())),
        ]).unwrap();
        let assignment = ElementMaterials::new(models, vec![MaterialId(1), MaterialId(2)]).unwrap();
        let p = ConductionProblem {
            element_materials: Some(&assignment),
            ..problem(&mesh, &boundary, &fallback, &source)
        };
        let storage = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh, &calorics, &[0, 1], budget(&mesh),
        ).unwrap();
        let old = [0.0, 0.0, 0.0, 0.0, 200.0, 200.0, 200.0, 200.0];
        let family = storage.isolated_equilibrium(cx, p, Some(&interface), &old, config()).unwrap();
        close(family.components()[0].temperature_k, 375.0, 1e-8);
        assert!(matches!(family.unique_state(cx),
            Err(IsolatedEquilibriumError::Original(EnthalpyError::Conduction(
                ConductionError::OutsideTemperatureSpan { .. })))));
        assert_eq!(old, [0.0, 0.0, 0.0, 0.0, 200.0, 200.0, 200.0, 200.0]);
    });
}
