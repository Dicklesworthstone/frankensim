//! A constrained synthetic closed material shell, not a measured instrument.
use super::*;
use fs_bem::panel3d::SpherePanels;
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::bem::{
    ShellBemError, ShellBemOptions,
};
use std::collections::BTreeMap;

fn soft_material(rho: f64) -> ResolvedMaterialStatePoint {
    soft_material_scaled(rho, 1.0)
}

fn soft_material_scaled(rho: f64, stiffness_scale: f64) -> ResolvedMaterialStatePoint {
    let mut entries = vec![
        ("density", Density::DIMS, rho),
        ("young_modulus_1", Pressure::DIMS, 1.0e6 * stiffness_scale),
        ("young_modulus_2", Pressure::DIMS, 4.0e5 * stiffness_scale),
        ("poisson_ratio_12", Dims::NONE, 0.25),
        ("shear_modulus_12", Pressure::DIMS, 1.8e5 * stiffness_scale),
        ("shear_modulus_13", Pressure::DIMS, 8.0e4 * stiffness_scale),
        ("shear_modulus_23", Pressure::DIMS, 3.0e4 * stiffness_scale),
    ];
    entries.extend(
        PLANE_STRESS_VISCOSITY_PROPERTIES
            .into_iter()
            .zip([1.0, 0.7, 0.3, 0.2, 0.04, -0.03])
            .map(|(key, value)| (key, DynViscosity::DIMS, value)),
    );
    entries.extend(
        TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES
            .into_iter()
            .zip([0.4, 0.8, 0.05])
            .map(|(key, value)| (key, DynViscosity::DIMS, value)),
    );
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in entries {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", 0.0, 100_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: "synthetic closed soft shell".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shell".into(),
            phase: "solid".into(),
            process: "fixture".into(),
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

fn closed_body(rho: f64, gate: &CancelGate) -> SharedMaterialShellBody {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let mut vertices = BTreeMap::new();
    let mut nodes = Vec::new();
    let tris = surface
        .triangles()
        .unwrap()
        .iter()
        .map(|tri| {
            tri.map(|p| {
                let key = p.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
                *vertices.entry(key).or_insert_with(|| {
                    nodes.push(p);
                    nodes.len() - 1
                })
            })
        })
        .collect();
    let node = (0..nodes.len())
        .max_by(|&a, &b| nodes[a][0].total_cmp(&nodes[b][0]))
        .unwrap();
    let vectors = nodes
        .iter()
        .map(|p| {
            let r = p[0].hypot(p[1]).hypot(p[2]);
            p.map(|x| x / r)
        })
        .collect();
    let mesh = ShellMesh::new(nodes, tris).unwrap();
    let directors = NodalDirectors::new(&mesh, vectors, 0.05).unwrap();
    // Only one physical nodal translation is free. This makes the source mode
    // independently identifiable and removes rigid modes without frequency tuning.
    let constraints = (0..mesh.nodes.len())
        .map(|i| {
            let mut mask = [true; 6];
            if i == node {
                mask[0] = false;
            }
            ShellConstraint {
                node: i,
                mask,
                frame: ConstraintFrame::Global,
            }
        })
        .collect();
    let budget = MitcBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
    };
    let region = ShellRegionMaterial {
        region: PlateRegionMaterial {
            region: PlateRegion {
                name: "skin".into(),
                triangle_indices: (0..mesh.tris.len()).collect(),
            },
            material: soft_material(rho),
            model: PlateMaterialModel::Orthotropic12,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.003),
            material_angle_rad: 0.17,
        },
        shear_correction: 5.0 / 6.0,
        axis1_global: [1.0, 0.173, 0.519],
    };
    let specimen = compile_material_shear_shell(mesh, vec![region], budget, gate).unwrap();
    let mut port = vec![0.0; 6 * specimen.mesh().nodes.len()];
    port[6 * node] = 1.0;
    SharedMaterialShellBody {
        specimen,
        directors,
        constraints,
        mass: ShellMass::Consistent,
        drilling_alpha: 1e-3,
        eigenvalue_window: (0.0, 1e9),
        n_modes: 1,
        force_ports: vec![port],
        bridges: vec![SharedPlateBridge {
            station_fraction: 0.8,
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }],
        assembly_budget: budget,
        projection_budget: ShellProjectionBudget {
            max_modes: 1,
            max_projection_terms: 128_000_000,
            reference_resistance_n_s_m: 1.0,
        },
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    }
}

fn closed_bank(cx: &Cx<'_>, gate: &CancelGate, rho: f64, pressure: f64) -> ReciprocalStringBank {
    closed_bank_in_gas(cx, gate, rho, gas(pressure))
}

fn closed_bank_in_gas(
    cx: &Cx<'_>,
    gate: &CancelGate,
    rho: f64,
    atmosphere: GasState,
) -> ReciprocalStringBank {
    let mut settings = options(0);
    settings.network.max_setup_terms = 131_072;
    ReciprocalStringBank::from_material_shear_shell(
        cx,
        vec![member(7800.0).source],
        atmosphere,
        closed_body(rho, gate),
        None,
        law(),
        settings,
        gate,
    )
    .unwrap()
}

fn resolution() -> ShellBemOptions {
    ShellBemOptions {
        minimum_hz: 10.0,
        maximum_hz: 750.0,
        frequency_samples: 65,
        maximum_surface_panels: 80,
        // Actual P1/P0 source work is about 10.54 billion for this 65-row sweep.
        maximum_sweep_terms: 12_000_000_000,
        maximum_poles: 64,
    }
}

#[test]
fn original_curved_material_modes_receive_reciprocal_air_reaction() {
    with_cx(|cx, gate| {
        let mut bare = closed_bank(cx, gate, 700.0, 101_325.0);
        let (mut loaded, report) = closed_bank(cx, gate, 700.0, 101_325.0)
            .with_material_shell_bem(resolution(), gate)
            .unwrap();
        assert_eq!(report.surface_panels, 80);
        assert_eq!(report.body_modes, 1);
        assert_eq!(report.reciprocal.unknowns,122);
        assert_eq!(report.sweep_terms,65*report.reciprocal.frequency_work);
        assert!(report.sweep_terms<=resolution().maximum_sweep_terms);
        assert!(report.acoustic_poles > 0 && report.acoustic_poles <= 64);
        assert!(report.minimum_panels_per_wavelength >= 6.0);
        assert!(
            report.maximum_condition_lower_bound.is_finite()
                && report.maximum_condition_lower_bound >= 1.0
        );
        assert!(report.reciprocity_relative_residual <= 1e-12);
        assert!(report.peak_error <= 0.05 && report.rms_error <= 0.02);
        assert!(report.resistance_peak_error <= 0.10 && report.resistance_rms_error <= 0.05);
        assert_eq!(
            report.basis_identity,
            loaded.source().body_radiation().unwrap().basis_identity()
        );
        assert_eq!(&report.gas, loaded.source().body_radiation().unwrap().gas());
        let binding = loaded.source().body_radiation().unwrap();
        let retained = binding.material_shell_bem().unwrap();
        assert_eq!(retained.loading_identity, report.loading_identity);
        assert_eq!(retained.basis_identity, report.basis_identity);
        assert_eq!(retained.gas, report.gas);
        assert_eq!(retained.options, report.options);
        // Even the same numerical model, supplied through the raw API, does not
        // establish that this source's closed-shell producer attached the load.
        let raw = closed_bank(cx, gate, 700.0, 101_325.0)
            .with_body_radiation_load(binding.model().clone(), 64, gate)
            .unwrap();
        assert!(raw.source().body_radiation().unwrap().material_shell_bem().is_none());
        assert_eq!(
            loaded.system().components()[1].modes(),
            bare.system().components()[1].modes()
        );
        let mut force = vec![0.0; loaded.system().mode_count()];
        force[4] = 0.02;
        let released = FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        };
        let mut pressure_difference = 0.0;
        for _ in 0..32 {
            let a = *bare.step(&force, released).unwrap();
            let b = *loaded.step(&force, released).unwrap();
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            pressure_difference += (a.observer_pressure_pa - b.observer_pressure_pa).abs();
        }
        assert!(pressure_difference > 1e-12);
        assert_ne!(
            loaded.system().components()[0].states(),
            bare.system().components()[0].states()
        );
        assert!(loaded.system().radiation_energy_j().unwrap() > 0.0);
        assert!(
            loaded
                .system()
                .coupled_frame()
                .unwrap()
                .radiation_dissipation_j
                > 0.0
        );
    });
}

#[test]
fn atmosphere_pressure_and_solid_density_propagate_without_reauthored_acoustics() {
    with_cx(|cx, gate| {
        let (ordinary, a) = closed_bank(cx, gate, 700.0, 101_325.0)
            .prepare_material_shell_bem(resolution(), gate)
            .unwrap();
        let (pressurized, b) = closed_bank(cx, gate, 700.0, 202_650.0)
            .prepare_material_shell_bem(resolution(), gate)
            .unwrap();
        let (dense, d) = closed_bank(cx, gate, 1400.0, 101_325.0)
            .prepare_material_shell_bem(resolution(), gate)
            .unwrap();
        // An additional named species changes the coherent thermodynamics and
        // transport together; no acoustic coefficient is supplied separately.
        use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
        use fs_qty::semantic::{Composition, CompositionBasis};
        let composition =
            Composition::new(CompositionBasis::MoleFraction, vec![0.4, 0.1, 0.5]).unwrap();
        let mixture = resolve_gas_mixture(
            &[
                GasSpecies::Nitrogen,
                GasSpecies::Oxygen,
                GasSpecies::CarbonDioxide,
            ],
            &composition,
            293.15,
            101_325.0,
        )
        .unwrap();
        let (mixed, m) = closed_bank_in_gas(cx, gate, 700.0, mixture)
            .prepare_material_shell_bem(resolution(), gate)
            .unwrap();
        assert_eq!(a.basis_identity, m.basis_identity);
        assert_ne!(a.loading_identity, m.loading_identity);
        assert_eq!(m.gas, mixture);
        assert_ne!(a.gas.sound_speed, m.gas.sound_speed);
        assert_ne!(a.gas.density, m.gas.density);
        assert_eq!(a.basis_identity, b.basis_identity);
        assert_ne!(a.basis_identity, d.basis_identity);
        assert_ne!(a.loading_identity, b.loading_identity);
        assert_ne!(a.loading_identity, d.loading_identity);
        assert_eq!(a.gas.sound_speed, b.gas.sound_speed);
        assert_eq!(b.gas.density, 2.0 * a.gas.density);
        for hz in [75.0, 250.0, 750.0] {
            let w = core::f64::consts::TAU * hz;
            let z = ordinary.impedance(w).unwrap()[0];
            let changed = mixed.impedance(w).unwrap()[0];
            assert!((changed - z).abs() > 1e-3 * z.abs());
            let p = pressurized.impedance(w).unwrap()[0];
            let h = dense.impedance(w).unwrap()[0];
            assert!((p - z.scale(2.0)).abs() < 1e-7 * z.abs());
            assert!((h - z.scale(0.5)).abs() < 1e-6 * z.abs());
        }
    });
}

#[test]
fn open_skin_incomplete_bands_work_limits_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let open = build_shell(cx, gate, 700.0, 0.02, 0.17);
        assert!(matches!(
            open.prepare_material_shell_bem(resolution(), gate),
            Err(ShellBemError::Solve(_))
        ));
        let closed = closed_bank(cx, gate, 700.0, 101_325.0);
        for case in 0..4 {
            let mut o = resolution();
            match case {
                0 => o.maximum_surface_panels = 79,
                1 => o.maximum_sweep_terms = 1,
                2 => {
                    o.minimum_hz = 1.0;
                    o.maximum_hz = 2.0;
                }
                _ => o.frequency_samples = 34,
            }
            assert!(matches!(
                closed.prepare_material_shell_bem(o, gate),
                Err(ShellBemError::Invalid(_))
            ));
        }
        gate.request();
        assert!(matches!(
            closed.prepare_material_shell_bem(resolution(), gate),
            Err(ShellBemError::Cancelled)
        ));
    });
}

#[path = "bem/openings.rs"]
mod openings;

#[path = "bem/skin.rs"]
mod skin;
