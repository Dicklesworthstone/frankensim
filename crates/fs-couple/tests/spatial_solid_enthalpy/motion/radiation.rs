//! G1/G3/G4: source radiation changes the continuing physical oscillator.
use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::{
    AmbientRadiationPatch, EMISSIVITY_DIMS, PointwiseAmbientRadiation, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ThermalBc,
};

const RESERVOIR_K: f64 = 450.0;

fn surface(emissivity: f64) -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: emissivity,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic gray coating for native motion integration".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray surface".into(),
            phase: "solid".into(),
            process: "manufactured integration reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    AmbientRadiationPatch::new(
        "radiating",
        SurfaceEmissivity::from_card("radiating", &card, 300.0, SelectionPolicy::SingleClaimOnly)
            .unwrap(),
        RESERVOIR_K,
    )
    .unwrap()
}

fn with_radiative_run<R>(
    emissivity: f64,
    f: impl FnOnce(
        &mut SpatialHeatedValveMotion<'_, '_, '_, '_, '_>,
        &ThermalBoundary,
        &PointwiseAmbientRadiation<'_>,
        &AmbientRadiationPatch,
        &Cx<'_>,
    ) -> R,
) -> R {
    with_cx(|cx| {
        let (card, caloric) = fixture(2.0, 1000.0, 1e-5, 10.0, 500.0, 750.0);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6.0, cx);
        let context = point(300.0);
        let source = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("radiating", ThermalBc::adiabatic())
            .unwrap()
            .finish()
            .unwrap();
        let patches = [surface(emissivity)];
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let fixed: Vec<_> = (0..11).collect();
        let pressure = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 1.]];
        let reduction = HeatedValveReduction {
            fixed_dofs: &fixed,
            orifice: OrificeSpec {
                face_nodes: &[1, 2, 3],
                opening_axis: [0., 0., 1.],
                width_axis: [1., 0., 0.],
                opposing_plane_offset_m: -0.5,
            },
            equilibrium_nodal_forces_n: &[(11, DEAD_FORCE)],
            nodal_load_per_pa_m2: &pressure,
            window_hz: (0.01, 20.0),
            retained_compliance_floor: 0.999,
            loss_factor: 0.0,
            source_id: "manufactured spatial-heat tetrahedron",
            static_solve: TetStaticSolveConfig::default(),
            maximum_nodal_displacement_m: 0.1,
        };
        let mut run = SpatialHeatedValveMotion::try_new(
            &source,
            vec![0.0; 4],
            reduction,
            RATE,
            &[ModalAcousticState {
                displacement_m_sqrt_kg: 3e-4,
                velocity_m_sqrt_kg_per_s: 4e-3,
            }],
            observation,
            ModalAcousticTimeBudget::audible_reference(),
            cx,
        )
        .unwrap();
        f(&mut run, &boundary, &radiation, &patches[0], cx)
    })
}

// Analytic P1 triangle mean of T^4: each degree-four homogeneous monomial
// has coefficient 1/15 after exact barycentric integration. This does not
// call the radiation owner's quadrature or substitute a mean temperature.
fn outward_heat(temperature: &[f64], emissivity: f64) -> f64 {
    let mut heat = 0.0;
    for (vertices, area) in [
        ([0, 1, 2], 0.5),
        ([0, 1, 3], 0.5),
        ([0, 2, 3], 0.5),
        ([1, 2, 3], 3.0_f64.sqrt() / 2.0),
    ] {
        let [a, b, c] = vertices.map(|i| temperature[i]);
        let mut quartic = 0.0;
        for i in 0..=4 {
            for j in 0..=4 - i {
                quartic += a.powi(i) * b.powi(j) * c.powi(4 - i - j);
            }
        }
        heat += area * (quartic / 15.0 - RESERVOIR_K.powi(4));
    }
    emissivity * STEFAN_BOLTZMANN_W_M2_K4 * heat
}

fn simulate(
    emissivity: f64,
    inject_refusals: bool,
) -> (Vec<f64>, Vec<ModalAcousticState>, Vec<f64>) {
    with_radiative_run(emissivity, |run, boundary, radiation, patch, cx| {
        let no_generation = ScalarField::uniform("no prescribed heating", 0.0).unwrap();
        let mut temperatures = vec![300.0; 4];
        let mut pressure_trace = Vec::new();
        for step in 1..=8 {
            let before_h = run.specific_enthalpy_j_kg().to_vec();
            let before = run.motion().clone();
            let current = before.nodal_motion(cx).unwrap();
            let (u, v) = (current.displacement_m[3][2], current.velocity_m_s[3][2]);
            let (equilibrium, k) = pencil(temperatures.iter().sum::<f64>() / 4.0);
            let omega = (k / 0.1).sqrt();
            let (sin, cos) = (omega * DT).sin_cos();
            let driven = equilibrium + 0.04 / k;
            let expected_u = driven + (u - driven) * cos + v * sin / omega;
            let expected_v = -(u - driven) * omega * sin + v * cos;
            if inject_refusals && step == 3 {
                let unchanged = |run: &SpatialHeatedValveMotion<'_, '_, '_, '_, '_>| {
                    assert_eq!(run.specific_enthalpy_j_kg(), before_h);
                    assert_eq!(run.motion().model().states(), before.model().states());
                    assert_eq!(run.motion().valve().identity(), before.valve().identity());
                    assert_eq!(run.motion().accepted_steps(), before.accepted_steps());
                };
                let failed = run.step_with_radiation(
                    0.04,
                    boundary,
                    None,
                    &no_generation,
                    thermal_config(),
                    EnthalpyRadiation::PointwiseAmbient(radiation),
                    controls(),
                    |_| Err::<Vec<C64>, _>("new observer refused"),
                    cx,
                );
                assert!(matches!(
                    failed,
                    Err(SpatialValveMotionError::Observation("new observer refused"))
                ));
                unchanged(run);
                let mut no_work = controls();
                no_work.maximum_abs_parameter_work_j = 0.0;
                assert!(
                    run.step_with_radiation(
                        0.04,
                        boundary,
                        None,
                        &no_generation,
                        thermal_config(),
                        EnthalpyRadiation::PointwiseAmbient(radiation),
                        no_work,
                        observation,
                        cx,
                    )
                    .is_err()
                );
                unchanged(run);
            }
            let result = run
                .step_with_radiation(
                    0.04,
                    boundary,
                    None,
                    &no_generation,
                    thermal_config(),
                    EnthalpyRadiation::PointwiseAmbient(radiation),
                    controls(),
                    observation,
                    cx,
                )
                .unwrap();
            let thermal = result.thermal.conduction();
            let report = thermal.pointwise_radiation.as_ref().unwrap();
            assert_eq!(report.patches.len(), 1);
            assert_eq!(&report.patches[0].patch, patch);
            let radiated_w = outward_heat(&thermal.temperature, emissivity);
            near(report.outward_heat_w, radiated_w, 1e-8);
            assert!(radiated_w < 0.0);
            assert_eq!(thermal.source_w, 0.0);
            assert_eq!(thermal.neumann_out_w, 0.0);
            assert_eq!(thermal.robin_out_w, 0.0);
            // Original reference inventory: each nodal mass is 1/4 kg.
            let inventory_change = thermal
                .specific_enthalpy_j_kg
                .iter()
                .zip(&before_h)
                .map(|(after, before)| 0.25 * (after - before))
                .sum::<f64>();
            near(thermal.stored_energy_change_j, inventory_change, 1e-10);
            near(inventory_change + DT * radiated_w, 0.0, 1e-8);
            assert!(result.thermal.physical_energy_residual_j().abs() < 1e-8);
            temperatures = thermal.temperature.clone();
            let (next_equilibrium, next_k) = pencil(temperatures.iter().sum::<f64>() / 4.0);
            let old_energy =
                0.05 * expected_v.powi(2) + 0.5 * k * (expected_u - equilibrium).powi(2);
            let new_energy =
                0.05 * expected_v.powi(2) + 0.5 * next_k * (expected_u - next_equilibrium).powi(2);
            let next = run.motion().nodal_motion(cx).unwrap();
            near(next.displacement_m[3][2], expected_u, 2e-12);
            near(next.velocity_m_s[3][2], expected_v, 2e-11);
            near(
                result.endpoint_observer_pressure_pa,
                2.0 * expected_v,
                4e-11,
            );
            near(
                result.transfer.parameter_work_j,
                new_energy - old_energy,
                1e-11,
            );
            assert_eq!(run.motion().accepted_steps(), step);
            near(
                run.motion().valve().equilibrium().pencil().total_mass_kg,
                1.0,
                1e-12,
            );
            pressure_trace.push(result.endpoint_observer_pressure_pa);
        }
        let spread = temperatures
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
            - temperatures.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            spread > 1e-4,
            "radiative hot spots must remain a real nodal field"
        );
        (
            run.specific_enthalpy_j_kg().to_vec(),
            run.motion().model().states().to_vec(),
            pressure_trace,
        )
    })
}

#[test]
fn g1_g3_g4_radiating_surface_changes_real_motion_and_pressure_with_atomic_retry() {
    let dark = simulate(0.8, false);
    let replay = simulate(0.8, true);
    assert_eq!(
        replay, dark,
        "late refusal must not consume heat, motion or observation history"
    );
    let light = simulate(0.2, false);
    assert!(dark.0.iter().sum::<f64>() > light.0.iter().sum::<f64>());
    assert!((dark.1[0].displacement_m_sqrt_kg - light.1[0].displacement_m_sqrt_kg).abs() > 1e-8);
    assert!((dark.2.last().unwrap() - light.2.last().unwrap()).abs() > 1e-8);
}
