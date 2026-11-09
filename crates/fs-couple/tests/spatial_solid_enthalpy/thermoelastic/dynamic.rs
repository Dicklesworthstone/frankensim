//! G1/G3/G4: real inertia, source thermoelastic heat and original P1 calorics.
use super::*;
use fs_couple::solid_enthalpy::spatial::thermoelastic::dynamic::SpatialRayleighDestination;
use fs_solid::loads3::transient::TransientConfig3;

const VOLUME: f64 = 1.0 / 6.0;
const MASS: f64 = RHO * VOLUME;

fn dynamics() -> TransientConfig3 {
    TransientConfig3 {
        maximum_steps: 128,
        maximum_displacement_m: 0.1,
        absolute_force_tolerance_n: 1e-10,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-9,
        ..TransientConfig3::default()
    }
}
fn fixed() -> Vec<usize> { (0..12).filter(|&dof| dof != 3).collect() }
fn velocity(speed: f64) -> [[f64; 3]; 4] {
    [[0.0; 3], [speed, 0.0, 0.0], [0.0; 3], [0.0; 3]]
}
fn stiffness(theta: f64) -> f64 {
    VOLUME * (E0 + ET * theta) * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU))
}
fn thermal_force(theta: f64, alpha: f64) -> f64 {
    VOLUME * (E0 + ET * theta) / (1.0 - 2.0 * NU) * free_strain(alpha, theta)
}

struct Endpoint {
    theta: f64,
    u: f64,
    v: f64,
    work: f64,
    internal: f64,
    heat: f64,
    kinetic: f64,
    damping: f64,
}

// Independent one-free-DOF Newmark equation with actual P1 consistent mass.
// The thermal equation below uses analytical Psi and entropy from the original
// E(T), nu and alpha(T), never the native exchange implementation.
#[allow(clippy::too_many_arguments)]
fn at_temperature(
    theta: f64, old: [f64; 3], forces: [f64; 2], alpha: f64,
    mass: f64, dt: f64, eta: f64, zeta: f64,
) -> Endpoint {
    let [theta0, u0, v0] = old;
    let k0 = stiffness(theta0);
    let k1 = stiffness(theta);
    let f0 = thermal_force(theta0, alpha);
    let f1 = thermal_force(theta, alpha);
    let c = eta * mass + 0.5 * zeta * (k0 + k1);
    let u = ((2.0 * mass / dt.powi(2) + c / dt - 0.5 * k0) * u0
        + 2.0 * mass / dt * v0 + 0.5 * (forces[0] + forces[1] + f0 + f1))
        / (2.0 * mass / dt.powi(2) + c / dt + 0.5 * k1);
    let v = 2.0 * (u - u0) / dt - v0;
    let work = 0.5 * (k0 * u0 - f0 + k1 * u - f1) * (u - u0);
    let internal = energy(theta, alpha, [u, 0.0, 0.0], VOLUME)
        - energy(theta0, alpha, [u0, 0.0, 0.0], VOLUME);
    Endpoint {
        theta, u, v, work, internal, heat: work - internal,
        kinetic: 0.5 * mass * (v * v - v0 * v0),
        damping: c * (u - u0).powi(2) / dt,
    }
}

#[test]
fn finite_heating_changes_actual_inertia_and_exports_only_declared_damping() {
    with_cx(|cx| {
        for alpha in [1e-3, -1e-3] {
            for (eta, zeta) in [(0.0, 0.0), (0.3, 0.01)] {
                let cp = 10.0;
                let (card, caloric) = chart(alpha, cp);
                let mesh = mesh();
                let native = native(&mesh, caloric.curve(), RHO, cx);
                let context = typed_point(T0);
                let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
                let geometry = binding.geometry(cx).unwrap();
                let initial = geometry.assemble(&[], VolumeLoad3::None, &[(3, 0.02)], cx).unwrap();
                let applied = geometry.assemble(&[], VolumeLoad3::None, &[(3, 0.03)], cx).unwrap();
                let old_h = [cp * (T0 - LOW); 4];
                let initial_time = 2.0;
                let model = ReciprocalSpatialSolid::try_new(
                    &binding, &old_h, initial_time, &initial, &fixed(), solid_config(),
                    physical_controls(), SpatialCaloricReference::StressFreeSolid, cx,
                ).unwrap().begin_transient(
                    &initial, &velocity(0.02),
                    TransientConfig3 { mass_damping_per_s: eta, stiffness_damping_s: zeta, ..dynamics() },
                    SpatialRayleighDestination::ExternalSink, cx,
                ).unwrap();
                let old_motion = model.nodal_motion(cx).unwrap();
                let dt = 0.03;
                assert_ne!(((initial_time + dt) - initial_time).to_bits(), dt.to_bits());
                let old = [0.0, 0.02 / stiffness(0.0), 0.02];
                let theta = root(-3.0, 3.0, |theta| {
                    let expected = at_temperature(theta, old, [0.02, 0.03], alpha, 0.1, dt, eta, zeta);
                    MASS * cp * theta - 2.0 * dt - expected.heat
                });
                let expected = at_temperature(theta, old, [0.02, 0.03], alpha, 0.1, dt, eta, zeta);
                let result = model.advance(
                    dt, &boundary(&mesh), None, &ScalarField::Uniform(12.0), &applied,
                    thermal_config(), physical_controls(), cx,
                ).unwrap();
                let motion = result.solid.nodal_motion(cx).unwrap();
                for &temperature in &result.thermal.endpoint.temperature {
                    near(temperature, T0 + expected.theta, 3e-8);
                }
                near(motion[0][1][0], expected.u, 3e-10);
                near(motion[1][1][0], expected.v, 2e-8);
                near(result.report.stress_work_j, expected.work, 2e-9);
                near(result.report.kinetic_energy_change_j, expected.kinetic, 2e-9);
                near(result.report.internal_energy_change_j, expected.internal, 2e-8);
                near(result.report.reversible_heat_j, expected.heat, 2e-8);
                near(result.report.damping_export_j, expected.damping, 2e-10);
                near(result.report.external_work_j, 0.025 * (expected.u - old[1]), 2e-10);
                near(result.report.external_thermal_heat_j, 2.0 * dt, 2e-10);
                near(result.report.mechanical_residual_j, 0.0, 2e-9);
                near(result.report.total_energy_residual_j, 0.0, 5e-8);
                for q in &result.report.nodal_reversible_heat_j {
                    near(q.to_f64(), expected.heat / 4.0, 1e-8);
                }
                assert_eq!(result.report.native.native.step_size_s.to_bits(), dt.to_bits());
                assert_eq!(result.solid.time_s(), initial_time + dt);
                assert_eq!(result.solid.step(), 1);
                assert_eq!(model.time_s(), initial_time);
                assert_eq!(model.step(), 0);
                assert_eq!(model.specific_enthalpy_j_kg(), &old_h);
                assert_eq!(model.nodal_motion(cx).unwrap(), old_motion);
                assert!(result.thermal.endpoint.maximum_free_residual_j <= thermal_config().energy_tolerance_j);
                let static_u = (0.03 + thermal_force(theta, alpha)) / stiffness(theta);
                assert!((expected.u - static_u).abs() > 1e-4);
            }
        }
    });
}

#[test]
fn already_heated_launch_keeps_clock_shape_and_exact_graded_consistent_mass() {
    with_cx(|cx| {
        let (alpha, cp) = (1e-3, 10.0);
        let (card, caloric) = chart(alpha, cp);
        let mesh = mesh();
        let density = [3.0, 6.0, 9.0, 12.0];
        let native = native(&mesh, caloric.curve(), RHO, cx)
            .with_nodal_reference_density(cx, &density).unwrap();
        let original_masses = native.reference_nodal_masses_kg().to_vec();
        let context = typed_point(T0);
        let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
        let geometry = binding.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let model = ReciprocalSpatialSolid::try_new(
            &binding, &[cp * (T0 - LOW); 4], 2.0, &loads, &fixed(), solid_config(),
            physical_controls(), SpatialCaloricReference::StressFreeSolid, cx,
        ).unwrap();
        let warm_boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("warm", ThermalBc::dirichlet(302.0).unwrap()).unwrap().finish().unwrap();
        let warm = model.advance(
            0.01, &warm_boundary, None, &ScalarField::Uniform(0.0), &loads,
            thermal_config(), physical_controls(), cx,
        ).unwrap();
        let warm_u = warm.solid.equilibrium().solution.displacement.displacement_m().to_vec();
        let warm_time = warm.solid.time_s();
        let moving = warm.solid.begin_transient(
            &loads, &velocity(0.01), dynamics(), SpatialRayleighDestination::ExternalSink, cx,
        ).unwrap();
        assert_eq!(moving.time_s(), warm_time);
        assert_eq!(moving.nodal_motion(cx).unwrap()[0], warm_u);
        // Integral rho_0*N_1^2 = V*(sum rho + 2 rho_1)/60, not lumped mass.
        let mass = VOLUME * (density.iter().sum::<f64>() + 2.0 * density[1]) / 60.0;
        near(moving.kinetic_energy_j(), 0.5 * mass * 0.01_f64.powi(2), 1e-15);
        let dt = 0.03;
        let old = [2.0, thermal_force(2.0, alpha) / stiffness(2.0), 0.01];
        let expected = at_temperature(10.0, old, [0.0; 2], alpha, mass, dt, 0.0, 0.0);
        let hot_boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("hot", ThermalBc::dirichlet(310.0).unwrap()).unwrap().finish().unwrap();
        let result = moving.advance(
            dt, &hot_boundary, None, &ScalarField::Uniform(0.0), &loads,
            thermal_config(), physical_controls(), cx,
        ).unwrap();
        let motion = result.solid.nodal_motion(cx).unwrap();
        near(motion[0][1][0], expected.u, 3e-10);
        near(motion[1][1][0], expected.v, 2e-8);
        near(result.report.kinetic_energy_change_j, expected.kinetic, 2e-9);
        near(result.report.reversible_heat_j, expected.heat, 2e-8);
        let total_mass = VOLUME * density.iter().sum::<f64>() / 4.0;
        near(result.report.external_thermal_heat_j, total_mass * cp * 8.0 - expected.heat, 5e-8);
        near(result.report.total_energy_residual_j, 0.0, 5e-8);
        assert_eq!(result.solid.time_s(), warm_time + dt);
        assert_eq!(result.report.native.native.step_size_s.to_bits(), dt.to_bits());
        assert_eq!(original_masses.as_slice(), native.reference_nodal_masses_kg());
        assert_eq!(result.solid.state().materials()[0].density_kg_m3, RHO);
        assert_eq!(moving.time_s(), warm_time);
        assert_eq!(moving.step(), 0);
    });
}

// Independent continuous one-DOF thermoelastic ODE. This reference is test-only.
// Cp*Tdot + C_el*Tdot = -T*S_u*udot for this closed, undamped specimen.
fn continuous_reference(alpha: f64, cp: f64, speed: f64, duration: f64) -> [f64; 3] {
    let derivative = |[theta, u, v]: [f64; 3]| {
        let e = E0 + ET * theta;
        let bulk = e / (3.0 * (1.0 - 2.0 * NU));
        let shear = e / (2.0 * (1.0 + NU));
        let bulk_t = ET / (3.0 * (1.0 - 2.0 * NU));
        let shear_t = ET / (2.0 * (1.0 + NU));
        let a = alpha * (1.0 + theta / 100.0);
        let dilation = u - 3.0 * free_strain(alpha, theta);
        let entropy_u = VOLUME * (3.0 * bulk * a - bulk_t * dilation - 4.0 * shear_t * u / 3.0);
        let elastic_capacity = (T0 + theta) * VOLUME
            * ((6.0 * bulk_t * a + 3.0 * bulk * alpha / 100.0) * dilation - 9.0 * bulk * a * a);
        [
            -(T0 + theta) * entropy_u * v / (MASS * cp + elastic_capacity),
            v,
            -VOLUME * (bulk * dilation + 4.0 * shear * u / 3.0) / 0.1,
        ]
    };
    let mut state = [0.0, 0.0, speed];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(state);
        let b = derivative(std::array::from_fn(|i| state[i] + 0.5 * dt * a[i]));
        let c = derivative(std::array::from_fn(|i| state[i] + 0.5 * dt * b[i]));
        let d = derivative(std::array::from_fn(|i| state[i] + dt * c[i]));
        state = std::array::from_fn(|i| state[i] + dt * (a[i] + 2.0*b[i] + 2.0*c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn closed_oscillation_refines_and_failed_continuations_preserve_the_accepted_prefix() {
    with_cx(|cx| {
        let (alpha, cp, speed) = (1e-3, 10.0, 0.02);
        let (card, caloric) = chart(alpha, cp);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), RHO, cx);
        let context = typed_point(T0);
        let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
        let geometry = binding.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let boundary = boundary(&mesh);
        let zero = ScalarField::Uniform(0.0);
        let omega = (stiffness(0.0) / 0.1).sqrt();
        let duration = 1.1 / omega;
        let exact = continuous_reference(alpha, cp, speed, duration);
        let displacement_scale = speed / omega;
        let temperature_scale = T0 * VOLUME * E0 / (1.0 - 2.0 * NU) * alpha
            * displacement_scale / (MASS * cp);
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let mut model = ReciprocalSpatialSolid::try_new(
                &binding, &[cp * (T0 - LOW); 4], 0.0, &loads, &fixed(), solid_config(),
                physical_controls(), SpatialCaloricReference::StressFreeSolid, cx,
            ).unwrap().begin_transient(
                &loads, &velocity(speed), TransientConfig3 { maximum_steps: count, ..dynamics() },
                SpatialRayleighDestination::ExternalSink, cx,
            ).unwrap();
            let dt = duration / count as f64;
            let mut bounded = physical_controls();
            bounded.maximum_iterations = 1;
            assert!(model.advance(dt, &boundary, None, &zero, &loads, thermal_config(), bounded, cx).is_err());
            assert_eq!(model.step(), 0);
            let mut heat_seen = false;
            for step in 1..=count {
                let result = model.advance(
                    dt, &boundary, None, &zero, &loads, thermal_config(), physical_controls(), cx,
                ).unwrap();
                heat_seen |= result.report.reversible_heat_j.abs() > 1e-7;
                near(result.report.external_thermal_heat_j, 0.0, 1e-14);
                near(result.report.external_work_j, 0.0, 1e-14);
                near(result.report.damping_export_j, 0.0, 0.0);
                near(result.report.total_energy_residual_j, 0.0, 5e-8);
                assert_eq!(result.solid.step(), step);
                model = result.solid;
            }
            assert!(heat_seen);
            let motion = model.nodal_motion(cx).unwrap();
            let theta = model.state().samples()[0].phase.temperature_k() - T0;
            let error = ((theta - exact[0]) / temperature_scale).abs()
                .max(((motion[0][1][0] - exact[1]) / displacement_scale).abs())
                .max(((motion[1][1][0] - exact[2]) / speed).abs());
            errors.push(error);
            let accepted_h = model.specific_enthalpy_j_kg().to_vec();
            let accepted_time = model.time_s();
            // Native step budget fails after a heat proposal, without resetting
            // its accepted count or publishing the proposed caloric endpoint.
            assert!(model.advance(
                dt, &boundary, None, &zero, &loads, thermal_config(), physical_controls(), cx,
            ).is_err());
            let gate = CancelGate::new();
            gate.request();
            ArenaPool::new(ArenaConfig::default()).scope(|arena| {
                let cancelled = Cx::new(
                    &gate, arena, StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
                    Budget::INFINITE, ExecMode::Deterministic,
                );
                assert!(matches!(model.advance(
                    dt, &boundary, None, &zero, &loads, thermal_config(), physical_controls(), &cancelled,
                ), Err(SpatialThermoelasticError::Cancelled)));
            });
            assert_eq!(model.specific_enthalpy_j_kg(), accepted_h.as_slice());
            assert_eq!(model.nodal_motion(cx).unwrap(), motion);
            assert_eq!(model.time_s(), accepted_time);
            assert_eq!(model.step(), count);
        }
        assert!(errors[0] > 1e-4);
        assert!(errors[1] < 0.6 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.6 * errors[1], "{errors:?}");
    });
}

fn gray_surface() -> fs_conduction::SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_conduction::EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: 0.8, dims: fs_conduction::EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", LOW, HIGH),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "analytical constant gray surface".into(), license: "test-only".into(), artifact: None,
        },
    }).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray".into(), phase: "solid".into(),
            process: "analytical radiation".into(), revision: 0,
        },
        claims, vec![],
    ).unwrap();
    fs_conduction::SurfaceEmissivity::from_card(
        "skin", &card, T0, SelectionPolicy::SingleClaimOnly,
    ).unwrap()
}

#[test]
fn native_pointwise_radiation_drives_equilibrium_and_real_inertial_heat_feedback() {
    use fs_conduction::{AmbientRadiationPatch, PointwiseAmbientRadiation, STEFAN_BOLTZMANN_W_M2_K4};
    with_cx(|cx| {
        let (alpha, cp, dt) = (1e-3, 10.0, 0.01);
        let (card, caloric) = chart(alpha, cp);
        let vertices = [
            [0.0; 3], [1.0, 0.0, 0.0], [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
        ];
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]), vertices.to_vec(),
        ).unwrap();
        let native = native(&mesh, caloric.curve(), RHO, cx);
        let context = typed_point(T0);
        let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
        let geometry = binding.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("skin", |_| true, ThermalBc::neumann(0.0).unwrap()).unwrap().finish().unwrap();
        let patches = [AmbientRadiationPatch::new("skin", gray_surface(), 320.0).unwrap()];
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let make = || ReciprocalSpatialSolid::try_new(
            &binding, &[cp * (T0 - LOW); 4], 0.0, &loads, &FREE, solid_config(),
            physical_controls(), SpatialCaloricReference::StressFreeSolid, cx,
        ).unwrap();
        let volume = 2.0_f64.sqrt() / 12.0;
        let capacity = RHO * cp * volume;
        let outward = |theta: f64| dt * 3.0_f64.sqrt() * 0.8 * STEFAN_BOLTZMANN_W_M2_K4
            * ((T0 + theta).powi(4) - 320.0_f64.powi(4));
        let expected = root(0.0, 2.0, |theta| capacity * theta + outward(theta));
        let equilibrium = make().advance_with_radiation(
            dt, &boundary, None, &ScalarField::Uniform(0.0), &loads,
            thermal_config(), physical_controls(), EnthalpyRadiation::PointwiseAmbient(&radiation), cx,
        ).unwrap();
        for &t in &equilibrium.thermal.endpoint.temperature { near(t, T0 + expected, 3e-8); }
        near(equilibrium.report.reversible_heat_j, 0.0, 2e-9);
        near(equilibrium.report.external_thermal_heat_j, -outward(expected), 2e-8);
        near(equilibrium.thermal.endpoint.energy.radiation_out_j.to_f64(), outward(expected), 2e-8);
        let moving = make().begin_transient(
            &loads, &[[0.0; 3]; 4], dynamics(), SpatialRayleighDestination::ExternalSink, cx,
        ).unwrap();
        let result = moving.advance_with_radiation(
            dt, &boundary, None, &ScalarField::Uniform(0.0), &loads,
            thermal_config(), physical_controls(), EnthalpyRadiation::PointwiseAmbient(&radiation), cx,
        ).unwrap();
        let theta = result.solid.state().samples()[0].phase.temperature_k() - T0;
        for &t in &result.thermal.endpoint.temperature { near(t, T0 + theta, 3e-8); }
        let strain = result.solid.response().stress.elements()[0].total_strain_mandel();
        // Extend the independent principal-strain energy to all Mandel shear
        // modes: their U coefficient is V*(G - T*G_T).
        let shear_storage = volume * (E0 - T0 * ET) / (2.0 * (1.0 + NU))
            * strain[3..].iter().map(|s| s * s).sum::<f64>();
        let internal = energy(theta, alpha, [strain[0], strain[1], strain[2]], volume) + shear_storage;
        near(result.report.internal_energy_change_j, internal, 2e-9);
        near(capacity * theta + internal + result.solid.kinetic_energy_j() + outward(theta), 0.0, 5e-8);
        near(result.report.external_thermal_heat_j, -outward(theta), 2e-8);
        near(result.thermal.endpoint.energy.radiation_out_j.to_f64(), outward(theta), 2e-8);
        near(result.report.total_energy_residual_j, 0.0, 5e-8);
        assert!(result.solid.kinetic_energy_j() > 0.0);
        assert!(result.report.reversible_heat_j.abs() > 1e-7);
        assert!(result.solid.nodal_motion(cx).unwrap()[0].iter().flatten().any(|u| u.abs() > 1e-10));
        assert_eq!(moving.step(), 0);
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(moving.specific_enthalpy_j_kg(), &[cp * (T0 - LOW); 4]);
    });
}
