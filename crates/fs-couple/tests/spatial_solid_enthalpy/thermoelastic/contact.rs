//! Live ordered contact heat changes native prestress and P1 inertial motion.
//! The matching interface traces are mechanically fixed. Only the two outer
//! apex DOFs move, so these tests do not infer a changing gap or pressure.
use super::*;
use fs_conduction::ThermalInterfaces;
use fs_conduction::interface::temperature::{
    TemperatureContact, TemperatureContacts, enthalpy::ContactEnthalpyNewtonConfig,
};
use fs_conduction::transient::enthalpy::nodal_heat::NodalHeatEndpoint;
use fs_couple::solid_enthalpy::spatial::thermoelastic::dynamic::SpatialRayleighDestination;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_solid::loads3::transient::TransientConfig3;

const CP: [f64; 2] = [10.0, 20.0];
const ALPHA: [f64; 2] = [1e-3, -5e-4];
const INITIAL: [f64; 8] = [310.0, 310.0, 310.0, 310.0, 290.0, 290.0, 290.0, 290.0];
const V: f64 = 1.0 / 6.0;
const DT: f64 = 0.2;

fn joined_mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0],
        ],
    ).unwrap()
}
fn fixed() -> Vec<usize> { (0..24).filter(|&d| d != 11 && d != 23).collect() }
fn old_h() -> [f64; 8] { std::array::from_fn(|i| CP[i / 4] * (INITIAL[i] - LOW)) }
fn resistance(a: f64, b: f64) -> f64 { 0.25 + 0.003 * (a - T0) + 0.001 * (b - T0) }

fn interface_card() -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Surface2D {
            axes: ["TA".into(), "TB".into()], axis_dims: [Temperature::DIMS; 2],
            samples: [[LOW, LOW], [310.25, LOW], [LOW, HIGH], [310.25, HIGH]]
                .map(|p| (p, resistance(p[0], p[1]))).to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2]], dims,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("TA", absolute(), LOW, 310.25)
            .with_quantity("TB", absolute(), LOW, HIGH)
            .with_quantity("normal_pressure", QuantitySpec::dimensional(Pressure::DIMS), 1e6, 1e6),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside, observations: vec![],
        provenance: Provenance {
            source: "manufactured ordered live contact; not measured data".into(),
            license: "test-only".into(), artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: "manufactured nodal solid".into(), phase: "solid".into(),
            process: "analytical source".into(), revision: 0,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        side("held upper trace"), side("held lower trace"),
        SystemContext {
            medium: "dry".into(), third_body: None,
            environment: "vacuum".into(), history: "unaged".into(),
        }, claims, vec![],
    ).unwrap()
}
fn contacts<'m, 'c>(
    cx: &Cx<'_>, mesh: &'m ConductionMesh, bc: &ThermalBoundary, card: &'c InterfaceSystemCard,
) -> TemperatureContacts<'m, 'c> {
    let mut pair = ThermalInterfaces::coincident_face_pairs(mesh).unwrap()[0];
    if mesh.boundary()[pair.side_a].vertices[0] >= 4 {
        std::mem::swap(&mut pair.side_a, &mut pair.side_b);
    }
    let point = QueryPoint::new().with_quantity("TA", absolute(), T0).unwrap()
        .with_quantity("TB", absolute(), T0).unwrap()
        .with_quantity("normal_pressure", QuantitySpec::dimensional(Pressure::DIMS), 1e6).unwrap();
    let law = TemperatureContact::new_two_sided(
        "joint", vec![pair], card, &point, ["TA", "TB"], SelectionPolicy::SingleClaimOnly,
    ).unwrap();
    TemperatureContacts::new(cx, mesh, bc, vec![], vec![law], 1).unwrap()
}
fn contact_controls() -> ContactEnthalpyNewtonConfig {
    let mut step = thermal_config();
    step.newton.relative_tolerance = 0.0;
    step.newton.absolute_tolerance = 1e-12;
    step.newton.linear_restart = 8;
    step.newton.max_linear_cycles = 16;
    step.max_newton_iterations = 64;
    ContactEnthalpyNewtonConfig {
        step, max_residual_evaluations: 1024, max_face_pairs: 1, residual_tolerance_j: 2e-10,
    }
}
fn coupling() -> ThermoelasticSolve {
    // The fresh-contact Euclidean target must still hold when the outer heat
    // iteration target and original max-row joule target are both looser.
    ThermoelasticSolve { heat_tolerance_j: 1e-4, ..physical_controls() }
}

struct Endpoint { q: f64, v: f64, work: f64, internal: f64, heat: f64, kinetic: f64 }
fn axial(t: f64, body: usize, dynamic: bool) -> Endpoint {
    let theta = t - T0;
    let theta0 = INITIAL[4 * body] - T0;
    let stiffness = |theta| V * (E0 + ET * theta) * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let force = |theta| V * (E0 + ET * theta) / (1.0 - 2.0 * NU) * free_strain(ALPHA[body], theta);
    let (k0, k1, f0, f1) = (stiffness(theta0), stiffness(theta), force(theta0), force(theta));
    let q0 = f0 / k0;
    // Integral rho*N_apex^2 = rho*V/10 = .1 kg, with all other DOFs fixed.
    let mass = RHO * V / 10.0;
    let q = if dynamic {
        ((2.0 * mass / DT.powi(2) - 0.5 * k0) * q0 + 0.5 * (f0 + f1))
            / (2.0 * mass / DT.powi(2) + 0.5 * k1)
    } else { f1 / k1 };
    let v = if dynamic { 2.0 * (q - q0) / DT } else { 0.0 };
    let work = 0.5 * (k0 * q0 - f0 + k1 * q - f1) * (q - q0);
    let internal = energy(theta, ALPHA[body], [q, 0.0, 0.0], V)
        - energy(theta0, ALPHA[body], [q0, 0.0, 0.0], V);
    Endpoint { q, v, work, internal, heat: work - internal, kinetic: 0.5 * mass * v * v }
}

// Independent eight-node equations: lumped original calorics, exact P1 tetra
// gradients, matching-triangle mass matrix, and analytical thermoelastic work.
// This tiny dense reference is test-only; production uses the original native
// contact Newton and native static/Newmark owners.
fn residual(t: &[f64; 8], dynamic: bool) -> [f64; 8] {
    let gram = [[3.0, -1.0, -1.0, -1.0], [-1.0, 1.0, 0.0, 0.0],
        [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]];
    let mut r = std::array::from_fn(|i| 0.25 * CP[i / 4] * (t[i] - INITIAL[i]));
    for body in 0..2 {
        let offset = 4 * body;
        let mean = t[offset..offset + 4].iter().sum::<f64>() / 4.0;
        let heat = axial(mean, body, dynamic).heat;
        for i in 0..4 {
            r[offset + i] += DT / 3.0 * (0..4).map(|j| gram[i][j] * t[offset + j]).sum::<f64>()
                - heat / 4.0;
        }
    }
    let lower = [4, 6, 5]; // Actual coincident vertices, not lower local order.
    let ta = t[..3].iter().sum::<f64>() / 3.0;
    let tb = t[4..7].iter().sum::<f64>() / 3.0;
    let jumps: [f64; 3] = std::array::from_fn(|i| t[i] - t[lower[i]]);
    for i in 0..3 {
        let transfer = DT * (jumps.iter().sum::<f64>() + jumps[i]) / (24.0 * resistance(ta, tb));
        r[i] += transfer;
        r[lower[i]] -= transfer;
    }
    r
}
fn reference(dynamic: bool) -> [f64; 8] {
    let mut t = INITIAL;
    for _ in 0..20 {
        let r = residual(&t, dynamic);
        if r.iter().all(|v| v.abs() < 1e-11) { return t; }
        let mut matrix = [[0.0; 8]; 8];
        for column in 0..8 {
            let mut plus = t;
            let mut minus = t;
            plus[column] += 1e-4;
            minus[column] -= 1e-4;
            let a = residual(&plus, dynamic);
            let b = residual(&minus, dynamic);
            for row in 0..8 { matrix[row][column] = (a[row] - b[row]) / 2e-4; }
        }
        let mut rhs = r.map(|v| -v);
        for column in 0..8 {
            let pivot = (column..8).max_by(|&a, &b|
                matrix[a][column].abs().total_cmp(&matrix[b][column].abs())).unwrap();
            matrix.swap(column, pivot);
            rhs.swap(column, pivot);
            assert!(matrix[column][column].abs() > 1.0);
            for row in column + 1..8 {
                let factor = matrix[row][column] / matrix[column][column];
                for j in column..8 { matrix[row][j] -= factor * matrix[column][j]; }
                rhs[row] -= factor * rhs[column];
            }
        }
        let mut delta = [0.0; 8];
        for row in (0..8).rev() {
            delta[row] = (rhs[row] - (row + 1..8).map(|j| matrix[row][j] * delta[j]).sum::<f64>())
                / matrix[row][row];
        }
        for i in 0..8 { t[i] += delta[i]; }
    }
    panic!("independent eight-node contact reference did not converge");
}
fn audit_contact(
    endpoint: &NodalHeatEndpoint, expected: &[f64; 8],
    contacts: &TemperatureContacts<'_, '_>, bc: &ThermalBoundary,
    card: &InterfaceSystemCard, cx: &Cx<'_>,
) -> f64 {
    for (&actual, &expected) in endpoint.temperature.iter().zip(expected) { near(actual, expected, 3e-7); }
    let actual = endpoint.temperature_contacts.as_ref().unwrap();
    let state = &actual.states()[0];
    let before = contacts.evaluate(cx, bc, &INITIAL).unwrap();
    let ta = endpoint.temperature[..3].iter().sum::<f64>() / 3.0;
    let tb = endpoint.temperature[4..7].iter().sum::<f64>() / 3.0;
    near(state.side_a_temperature_k, ta, 2e-12);
    near(state.side_b_temperature_k, tb, 2e-12);
    near(state.resistance.value_m2_k_per_w(), resistance(ta, tb), 2e-14);
    assert!((state.resistance.value_m2_k_per_w() - resistance(310.0, 290.0)).abs() > 1e-4);
    assert_eq!(state.resistance.card_identity(), card.content_hash());
    assert_eq!(state.resistance.receipt().selected, before.states()[0].resistance.receipt().selected);
    assert!(state.resistance.receipt().query_point.contains(&("normal_pressure".into(), 1e6)));
    assert!(endpoint.residual_norm_j <= contact_controls().residual_tolerance_j);
    assert!(endpoint.maximum_free_residual_j <= contact_controls().step.energy_tolerance_j);
    let transfer = endpoint.contact_fluxes[0].heat_rate_a_to_b_w;
    near(transfer, 0.5 * (ta - tb) / resistance(ta, tb), 2e-10);
    assert!(transfer > 0.0);
    transfer
}

#[test]
fn live_ordered_contact_drives_native_prestress_with_local_reciprocal_heat() {
    with_cx(|cx| {
        let (card_a, caloric_a) = chart(ALPHA[0], CP[0]);
        let (card_b, caloric_b) = chart(ALPHA[1], CP[1]);
        let mesh = joined_mesh();
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: caloric_a.curve(), reference_density_kg_m3: RHO },
              ReferenceEnthalpyMaterial { curve: caloric_b.curve(), reference_density_kg_m3: RHO }],
            &[0, 1], EnthalpyBudget { max_vertices: 8, max_elements: 2 },
        ).unwrap();
        let original_mass = native.reference_nodal_masses_kg().to_vec();
        let context = typed_point(T0);
        let source = SpatialSolidEnthalpy::try_new(
            &native,
            vec![original_source(&card_a, caloric_a.curve(), &context, CP[0]),
                 original_source(&card_b, caloric_b.curve(), &context, CP[1])],
            TetAssemblyBudget::standard(), solid_config().mesh_frame, 0.1, cx,
        ).unwrap();
        let geometry = source.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let bc = boundary(&mesh);
        let card = interface_card();
        let contacts = contacts(cx, &mesh, &bc, &card);
        let model = ReciprocalSpatialSolid::try_new(
            &source, &old_h(), 1.0, &loads, &fixed(), solid_config(), coupling(),
            SpatialCaloricReference::StressFreeSolid, cx,
        ).unwrap();
        let old_u = model.equilibrium().solution.displacement.displacement_m().to_vec();
        let result = model.advance_with_temperature_contacts(
            DT, &bc, &contacts, &ScalarField::Uniform(0.0), &loads,
            contact_controls(), coupling(), EnthalpyRadiation::None, cx,
        ).unwrap();
        let expected = reference(false);
        let transfer = audit_contact(&result.thermal.endpoint, &expected, &contacts, &bc, &card, cx);
        let u = result.solid.equilibrium().solution.displacement.displacement_m();
        for body in 0..2 {
            let e = axial(expected[4 * body..4 * body + 4].iter().sum::<f64>() / 4.0, body, false);
            near(if body == 0 { u[3][2] } else { -u[7][2] }, e.q, 3e-10);
            near(result.report.elements[body].reversible_heat_j, e.heat, 2e-8);
            let stress = result.solid.equilibrium().solution.stress.elements()[body].stress_mandel_pa();
            let young = E0 + ET * (expected[4 * body..4 * body + 4].iter().sum::<f64>() / 4.0 - T0);
            near(stress[0], -young * e.q / (1.0 + NU), 2e-8);
            near(stress[1], stress[0], 2e-10);
            near(stress[2], 0.0, 2e-9);
            let caloric = 0.25 * CP[body] * (0..4).map(|i|
                result.thermal.endpoint.temperature[4 * body + i] - INITIAL[4 * body + i]).sum::<f64>();
            near(caloric + DT * transfer * if body == 0 { 1.0 } else { -1.0 } - e.heat, 0.0, 2e-8);
            assert_eq!(result.solid.state().samples()[body].elastic.card_identity(),
                if body == 0 { card_a.content_hash() } else { card_b.content_hash() });
        }
        assert!(result.report.elements.iter().any(|e| e.reversible_heat_j.abs() > 1e-6));
        near(result.report.external_thermal_heat_j, 0.0, 1e-12);
        near(result.report.external_work_j, 0.0, 1e-12);
        near(result.report.total_energy_residual_j, 0.0, 5e-8);
        assert_eq!(original_mass, native.reference_nodal_masses_kg());
        assert_eq!(model.specific_enthalpy_j_kg(), &old_h());
        assert_eq!(model.equilibrium().solution.displacement.displacement_m(), old_u);
        assert_eq!(model.time_s(), 1.0);
        assert_eq!(result.solid.time_s(), 1.0 + DT);
    });
}

#[test]
fn live_contact_drives_inertia_and_late_source_exit_preserves_the_accepted_pair() {
    with_cx(|cx| {
        let (card_a, caloric_a) = chart(ALPHA[0], CP[0]);
        let (card_b, caloric_b) = chart(ALPHA[1], CP[1]);
        let mesh = joined_mesh();
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: caloric_a.curve(), reference_density_kg_m3: RHO },
              ReferenceEnthalpyMaterial { curve: caloric_b.curve(), reference_density_kg_m3: RHO }],
            &[0, 1], EnthalpyBudget { max_vertices: 8, max_elements: 2 },
        ).unwrap();
        let original_mass = native.reference_nodal_masses_kg().to_vec();
        let context = typed_point(T0);
        let source = SpatialSolidEnthalpy::try_new(
            &native,
            vec![original_source(&card_a, caloric_a.curve(), &context, CP[0]),
                 original_source(&card_b, caloric_b.curve(), &context, CP[1])],
            TetAssemblyBudget::standard(), solid_config().mesh_frame, 0.1, cx,
        ).unwrap();
        let geometry = source.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let bc = boundary(&mesh);
        let card = interface_card();
        let contacts = contacts(cx, &mesh, &bc, &card);
        let model = ReciprocalSpatialSolid::try_new(
            &source, &old_h(), 1.0, &loads, &fixed(), solid_config(), coupling(),
            SpatialCaloricReference::StressFreeSolid, cx,
        ).unwrap().begin_transient(
            &loads, &[[0.0; 3]; 8],
            TransientConfig3 {
                maximum_steps: 8, maximum_displacement_m: 0.1,
                absolute_force_tolerance_n: 1e-10, absolute_energy_tolerance_j: 1e-10,
                relative_energy_tolerance: 1e-9, ..TransientConfig3::default()
            }, SpatialRayleighDestination::ExternalSink, cx,
        ).unwrap();
        let old_motion = model.nodal_motion(cx).unwrap();
        let result = model.advance_with_temperature_contacts(
            DT, &bc, &contacts, &ScalarField::Uniform(0.0), &loads,
            contact_controls(), coupling(), EnthalpyRadiation::None, cx,
        ).unwrap();
        let expected = reference(true);
        let transfer = audit_contact(&result.thermal.endpoint, &expected, &contacts, &bc, &card, cx);
        let motion = result.solid.nodal_motion(cx).unwrap();
        let mut totals = [0.0; 4];
        for body in 0..2 {
            let e = axial(expected[4 * body..4 * body + 4].iter().sum::<f64>() / 4.0, body, true);
            let (node, sign) = if body == 0 { (3, 1.0) } else { (7, -1.0) };
            near(sign * motion[0][node][2], e.q, 3e-10);
            near(sign * motion[1][node][2], e.v, 3e-9);
            let caloric = 0.25 * CP[body] * (0..4).map(|i|
                result.thermal.endpoint.temperature[4 * body + i] - INITIAL[4 * body + i]).sum::<f64>();
            near(caloric + e.internal + e.kinetic + DT * transfer * sign, 0.0, 2e-8);
            for (sum, value) in totals.iter_mut().zip([e.work, e.internal, e.heat, e.kinetic]) { *sum += value; }
            assert_eq!(result.solid.state().samples()[body].elastic.card_identity(),
                if body == 0 { card_a.content_hash() } else { card_b.content_hash() });
        }
        near(result.report.stress_work_j, totals[0], 2e-9);
        near(result.report.internal_energy_change_j, totals[1], 2e-8);
        near(result.report.reversible_heat_j, totals[2], 2e-8);
        near(result.report.kinetic_energy_change_j, totals[3], 2e-9);
        assert!(result.report.kinetic_energy_change_j > 1e-8);
        near(result.report.damping_export_j, 0.0, 0.0);
        near(result.report.external_work_j, 0.0, 1e-12);
        near(result.report.external_thermal_heat_j, 0.0, 1e-12);
        near(result.report.total_energy_residual_j, 0.0, 5e-8);
        assert_eq!(result.report.native.native.step_size_s.to_bits(), DT.to_bits());
        assert_eq!(original_mass, native.reference_nodal_masses_kg());
        assert_eq!(model.nodal_motion(cx).unwrap(), old_motion);
        assert_eq!(model.specific_enthalpy_j_kg(), &old_h());
        assert_eq!(model.step(), 0);

        let accepted = result.solid;
        let accepted_h = accepted.specific_enthalpy_j_kg().to_vec();
        let accepted_motion = accepted.nodal_motion(cx).unwrap();
        // Bulk support extends to 330 K, but this ORIGINAL contact claim ends
        // at TA=310.25 K. A failed hot continuation cannot replace that claim
        // or publish only the heat or only the motion candidate.
        let outside_contact = [311.0, 311.0, 311.0, 311.0, 295.0, 295.0, 295.0, 295.0];
        let supported_h: [f64; 8] = std::array::from_fn(|i| CP[i / 4] * (outside_contact[i] - LOW));
        assert!(source.resolve_at(&supported_h, accepted.time_s() + DT, cx).is_ok());
        assert!(contacts.evaluate(cx, &bc, &outside_contact).is_err());
        let hot = ScalarField::Nodal(vec![2400.0, 2400.0, 2400.0, 2400.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(accepted.advance_with_temperature_contacts(
            DT, &bc, &contacts, &hot, &loads,
            contact_controls(), coupling(), EnthalpyRadiation::None, cx,
        ).is_err());
        let limited = ContactEnthalpyNewtonConfig { max_residual_evaluations: 1, ..contact_controls() };
        assert!(accepted.advance_with_temperature_contacts(
            DT, &bc, &contacts, &ScalarField::Uniform(0.0), &loads,
            limited, coupling(), EnthalpyRadiation::None, cx,
        ).is_err());
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate, arena, StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic,
            );
            assert!(matches!(accepted.advance_with_temperature_contacts(
                DT, &bc, &contacts, &ScalarField::Uniform(0.0), &loads,
                contact_controls(), coupling(), EnthalpyRadiation::None, &cancelled,
            ), Err(SpatialThermoelasticError::Cancelled)));
        });
        assert_eq!(accepted.specific_enthalpy_j_kg(), accepted_h);
        assert_eq!(accepted.nodal_motion(cx).unwrap(), accepted_motion);
        assert_eq!(accepted.time_s(), 1.0 + DT);
        assert_eq!(accepted.step(), 1);
    });
}
