//! G1/G3/G4: native field/material -> original string motion -> Faraday voltage.
//! Material cards are manufactured equations, not calibrated pickup measurements.
use fs_couple::{
    acoustic_realize::LinearMaterialStringRuntime,
    pickup::{
        PickupPose,
        field::{MaterialPickupMap, PickupFieldConfig, PickupFieldError, PlanarPickupSpec},
    },
    string_specimen::{ResolvedStringSpecimen, with_uniform_circular_material_state},
};
use fs_em::{
    planar::permanent::{PlanarMagneticMaterial, PlanarRecoilMaterial},
    winding::WindingConfig,
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    gas::{GasSpec, GasState},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_scenario::{Pluck, PrestressedString};
use fs_solver::CheckedCgConfig;
use std::sync::Arc;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 54,
                    kernel_id: 12,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn material(mu: f64, rho: f64, br: Option<f64>) -> ResolvedMaterialStatePoint {
    let mut values = vec![
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, 2e11),
        ("kelvin_voigt_bending_viscosity", DynViscosity::DIMS, 10.0),
    ];
    if let Some(br) = br {
        values.extend([
            (RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
            (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, br),
            (
                RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                -1e7,
            ),
            (
                RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                1e7,
            ),
            (
                RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                1e7,
            ),
        ]);
    } else {
        values.push((
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            mu,
        ));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in values {
        let mut validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        if name == "kelvin_voigt_bending_viscosity" {
            validity = validity.with("omega", 0.0, 100_000.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured magnetic/string regression".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "same manufactured material label".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn induced(mu: f64) -> PlanarMagneticMaterial {
    PlanarMagneticMaterial::Linear(material(mu, 7800.0, None))
}
fn magnet(br: f64) -> PlanarMagneticMaterial {
    PlanarMagneticMaterial::Recoil(
        PlanarRecoilMaterial::new(
            ResolvedRecoilMaterial::try_from_material(material(2e-6, 7800.0, Some(br))).unwrap(),
            [0.0, 1.0],
        )
        .unwrap(),
    )
}
fn config() -> PickupFieldConfig {
    PickupFieldConfig {
        max_field_evaluations: 128,
        derivative_step_m: 1e-5,
        derivative_tolerance_wb_m: 1e-8,
        interpolation_tolerance_wb: 1e-10,
        interpolation_slope_tolerance_wb_m: 1e-6,
        assembly: WeightedAssemblyLimits {
            max_cells: 1000,
            max_dofs: 2000,
            max_triplets: 9000,
        },
        winding: WindingConfig {
            max_windings: 1,
            max_basis_values: 10_000,
            pairing_tolerance_h: 1e-10,
            current_balance_tolerance_turns: 1e-9,
            inverse_residual_tolerance: 1e-10,
        },
        bias_stationarity_tolerance_wb: 1e-9,
        linear: CheckedCgConfig {
            tolerance: 1e-11,
            max_iterations: 1000,
            max_corrections: 2,
        },
    }
}
fn pose(height: f64) -> PickupPose {
    PickupPose {
        station_frac: 0.25,
        height_m: height,
        aperture_frac: 0.1,
    }
}
fn layered(mu: f64, br: f64) -> PlanarPickupSpec {
    let mut vertices = Vec::new();
    for y in [0.0, 0.02] {
        for x in [0.0, 0.01, 0.02, 0.03, 0.04] {
            vertices.push([x, y]);
        }
    }
    let mut faces = Vec::new();
    for i in 0..4_u32 {
        faces.extend([[i, i + 1, i + 6], [i, i + 6, i + 5]]);
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "material-pickup-layers",
        vertices,
        faces,
        Metric2::planar(0.03).unwrap(),
    )
    .unwrap();
    let mut motion = vec![[0.0; 2]; 10];
    motion[2] = [1.0, 0.0];
    motion[7] = [1.0, 0.0];
    // Moving leading interface changes the gap while the opposite face is fixed.
    // This is an exact layered field test, not a rigid translating string model.
    PlanarPickupSpec {
        mesh,
        reference_gap_m: 0.01,
        vertex_gap_directions: motion,
        materials: vec![
            magnet(br),
            magnet(br),
            induced(1e-6),
            induced(1e-6),
            induced(mu),
            induced(mu),
            induced(1e-6),
            induced(1e-6),
        ],
        winding_profile_per_m2: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3e4, -3e4],
        prescribed_vertices: vec![0, 5, 4, 9],
        end_inductance_h: 0.0,
        sensed_cells: vec![4, 5],
    }
}
fn expected_linkage(gap: f64, mu: f64, br: f64) -> f64 {
    // Tangential H is common to the layers. Integral B_y dx=0 because both
    // end potentials are fixed, hence H=-Br*l_magnet / sum(mu_i*l_i).
    // Exact P1 triangle integrals of the signed winding profile give A.J dV.
    let delta = gap - 0.01;
    let h = -br * 0.01 / (2e-6 * 0.01 + 1e-6 * (0.02 + delta) + mu * (0.01 - delta));
    -3e4 * 0.03 * 0.01_f64.powi(2) * 0.02 / 6.0 * 1e-6 * h
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn original_material_fields_supply_linkage_faraday_and_reciprocity() {
    with_cx(|cx, _| {
        let map =
            MaterialPickupMap::compile(cx, layered(3e-6, 0.1), &[0.009, 0.01, 0.011], config())
                .unwrap();
        assert_eq!(map.samples().len(), 25);
        for sample in map.samples() {
            near(
                sample.permanent_linkage_wb,
                expected_linkage(sample.gap_m, 3e-6, 0.1),
                1e-14,
            );
            assert!(sample.bias_report.as_ref().unwrap().converged_euclidean());
            assert!(
                sample
                    .unit_current_report
                    .as_ref()
                    .unwrap()
                    .converged_euclidean()
            );
        }
        let pickup = map.bind(pose(0.01), 3).unwrap();
        let q = [2e-5, -1e-5, 0.0];
        let v = [0.1, -0.03, 0.02];
        let gap = 0.01
            - q.iter()
                .zip(pickup.mode_gains())
                .map(|(q, g)| q * g)
                .sum::<f64>();
        near(
            pickup.flux_linkage_wb(&q).unwrap(),
            expected_linkage(gap, 3e-6, 0.1),
            1e-10,
        );
        let emf = pickup.emf_at_v(&q, &v).unwrap();
        assert!(emf.abs() > 1e-7);
        let dt = 1e-7;
        let q1 = std::array::from_fn::<_, 3, _>(|i| q[i] + dt * v[i]);
        near(pickup.emf_over_step_v(&q, &q1, dt).unwrap(), emf, 1e-10);
        let mut f = [0.0; 3];
        pickup.modal_step_force_n(&q, &q1, 0.002, &mut f).unwrap();
        near(
            f.iter()
                .enumerate()
                .map(|(i, f)| f * (q1[i] - q[i]))
                .sum::<f64>()
                + pickup.emf_over_step_v(&q, &q1, dt).unwrap() * 0.002 * dt,
            0.0,
            1e-24,
        );
        assert_eq!(map.source().materials[4], layered(3e-6, 0.1).materials[4]);
    });
}

#[test]
fn changing_permeability_or_remanence_changes_voltage_without_a_material_name_switch() {
    with_cx(|cx, _| {
        let build = |mu, br| {
            MaterialPickupMap::compile(cx, layered(mu, br), &[0.009, 0.01, 0.011], config())
                .unwrap()
        };
        let a = build(3e-6, 0.1);
        let b = build(3e-6, 0.2);
        let neutral = build(1e-6, 0.1);
        let voltage = |map: &MaterialPickupMap| {
            map.bind(pose(0.01), 2)
                .unwrap()
                .emf_at_v(&[0.0; 2], &[0.1, 0.0])
                .unwrap()
        };
        near(voltage(&b), 2.0 * voltage(&a), 1e-13);
        assert!(voltage(&a).abs() > 1e-7);
        assert!(
            voltage(&neutral).abs() < 1e-12,
            "no magnetic contrast means no interface-motion signal"
        );
        // Same physical names; the immutable property data alone changes.
        let PlanarMagneticMaterial::Linear(x) = &a.source().materials[4] else {
            panic!()
        };
        let PlanarMagneticMaterial::Linear(y) = &neutral.source().materials[4] else {
            panic!()
        };
        assert_eq!(x.material(), y.material());
        assert_ne!(x.card_identity(), y.card_identity());
    });
}

fn translating(body: ResolvedMaterialStatePoint) -> PlanarPickupSpec {
    let (nx, ny) = (6usize, 4usize);
    let vertex = |i: usize, j: usize| j * (nx + 1) + i;
    let mut vertices = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([i as f64 * 0.01, j as f64 * 0.01]);
        }
    }
    let mut faces = Vec::new();
    let mut materials = Vec::new();
    let mut profile = Vec::new();
    let mut sensed_cells = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let (a, b, c, d) = (
                vertex(i, j) as u32,
                vertex(i + 1, j) as u32,
                vertex(i + 1, j + 1) as u32,
                vertex(i, j + 1) as u32,
            );
            let model = if i == 1 && (j == 1 || j == 2) {
                magnet(0.1)
            } else if i == 3 && j == 1 {
                sensed_cells.extend([faces.len(), faces.len() + 1]);
                PlanarMagneticMaterial::Linear(body.clone())
            } else {
                induced(1e-6)
            };
            faces.extend([[a, b, c], [a, c, d]]);
            materials.extend([model.clone(), model]);
            let density = if i == 0 && (j == 2 || j == 3) {
                3e4
            } else if i == 5 && (j == 0 || j == 1) {
                -3e4
            } else {
                0.0
            };
            profile.extend([density, density]);
        }
    }
    let mut motion = vec![[0.0; 2]; vertices.len()];
    for j in 1..=2 {
        for i in 3..=4 {
            motion[vertex(i, j)] = [0.0, 1.0];
        }
    }
    let prescribed_vertices = (0..vertices.len())
        .filter(|&v| {
            v % (nx + 1) == 0 || v % (nx + 1) == nx || v / (nx + 1) == 0 || v / (nx + 1) == ny
        })
        .collect();
    let mesh = TriComplex2::from_indexed_triangles(
        "material-pickup-translating",
        vertices,
        faces,
        Metric2::planar(0.03).unwrap(),
    )
    .unwrap();
    PlanarPickupSpec {
        mesh,
        reference_gap_m: 0.005,
        vertex_gap_directions: motion,
        materials,
        winding_profile_per_m2: profile,
        prescribed_vertices,
        end_inductance_h: 0.0,
        sensed_cells,
    }
}
fn specimen(source: &ResolvedMaterialStatePoint) -> ResolvedStringSpecimen {
    with_uniform_circular_material_state(
        PrestressedString {
            length_m: 0.5,
            tension_n: 20.0,
            lin_density_kg_m: 1.0,
            axial_stiffness_n: 1.0,
            width_m: 1.0,
            n_modes: 4,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor: None,
            bending_loss_factor: None,
        },
        0.0003,
        source,
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap()
}
fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

#[test]
fn field_observer_reads_existing_material_string_and_rejects_stale_material() {
    with_cx(|cx, _| {
        let source = material(3e-6, 7800.0, None);
        let mut c = config();
        c.derivative_step_m = 1e-6;
        let map = Arc::new(
            MaterialPickupMap::compile(
                cx,
                translating(source.clone()),
                &[0.0049, 0.005, 0.0051],
                c,
            )
            .unwrap(),
        );
        let mut string = LinearMaterialStringRuntime::try_new_with_gas(
            cx,
            specimen(&source),
            Some(Pluck {
                station_frac: 0.2,
                height_m: 1e-5,
            }),
            &gas(101325.0),
            1.0,
            48000,
        )
        .unwrap();
        let observer = map.bind_material_string(pose(0.005), &string).unwrap();
        let pickup = map.bind(pose(0.005), string.states().len()).unwrap();
        let mass_root = (string.specimen().string().lin_density_kg_m * 0.5 * 0.5).sqrt();
        let mut maximum_voltage = 0.0_f64;
        for _ in 0..64 {
            string.step(cx, &[0.0; 4]).unwrap();
            let before = string.states().to_vec();
            let epoch = string.epoch();
            let observation = observer.observe(cx, &string).unwrap();
            let q = before
                .iter()
                .map(|s| s.displacement_m_sqrt_kg / mass_root)
                .collect::<Vec<_>>();
            let v = before
                .iter()
                .map(|s| s.velocity_m_sqrt_kg_per_s / mass_root)
                .collect::<Vec<_>>();
            near(observation.emf_v, pickup.emf_at_v(&q, &v).unwrap(), 1e-18);
            maximum_voltage = maximum_voltage.max(observation.emf_v.abs());
            assert_eq!(observation.epoch, epoch);
            assert_eq!(string.states(), before);
            assert_eq!(string.epoch(), epoch);
        }
        assert!(maximum_voltage > 1e-12);
        string.rebind_gas(cx, &gas(70000.0)).unwrap();
        observer.observe(cx, &string).unwrap();
        let replacement = material(4e-6, 9000.0, None);
        string
            .substitute_material(cx, specimen(&replacement), 1.0)
            .unwrap();
        assert!(observer.observe(cx, &string).is_err());
        assert!(map.bind_material_string(pose(0.005), &string).is_err());
        let rebuilt = Arc::new(
            MaterialPickupMap::compile(cx, translating(replacement), &[0.0049, 0.005, 0.0051], c)
                .unwrap(),
        );
        let new_observer = rebuilt.bind_material_string(pose(0.005), &string).unwrap();
        assert!(new_observer.observe(cx, &string).unwrap().emf_v.is_finite());
        let new_pickup = rebuilt.bind(pose(0.005), 4).unwrap();
        let q = [0.0; 4];
        let v = [0.1, 0.0, 0.0, 0.0];
        assert!(
            (new_pickup.emf_at_v(&q, &v).unwrap() - pickup.emf_at_v(&q, &v).unwrap()).abs() > 1e-12,
            "recompiling only the source material data must change the electrical response"
        );
    });
}

#[test]
fn field_resolution_source_correspondence_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let gaps = [0.009, 0.01, 0.011];
        let mut c = config();
        c.max_field_evaluations = 24;
        assert!(MaterialPickupMap::compile(cx, layered(3e-6, 0.1), &gaps, c).is_err());
        let mut c = config();
        c.derivative_tolerance_wb_m = 1e-16;
        assert!(matches!(
            MaterialPickupMap::compile(cx, layered(3e-6, 0.1), &gaps, c),
            Err(PickupFieldError::Derivative { .. })
        ));
        let mut c = config();
        c.interpolation_tolerance_wb = 1e-16;
        assert!(matches!(
            MaterialPickupMap::compile(cx, layered(3e-6, 0.1), &[0.005, 0.015], c),
            Err(PickupFieldError::Interpolation { .. })
        ));
        let mut wrong = layered(3e-6, 0.1);
        wrong.vertex_gap_directions[3] = [1.0, 0.0];
        assert!(MaterialPickupMap::compile(cx, wrong, &gaps, config()).is_err());
        let mut wrong = layered(3e-6, 0.1);
        wrong.materials.pop();
        assert!(MaterialPickupMap::compile(cx, wrong, &gaps, config()).is_err());
        assert!(
            MaterialPickupMap::compile(cx, layered(3e-6, 0.1), &[0.009, 0.025], config()).is_err()
        );
        gate.request();
        assert!(matches!(
            MaterialPickupMap::compile(cx, layered(3e-6, 0.1), &gaps, config()),
            Err(PickupFieldError::Cancelled)
        ));
    });
}
