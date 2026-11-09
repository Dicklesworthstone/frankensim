use super::*;
use crate::vibroacoustic::{
    CavityOpeningRadiation, DistributedCavityOpeningRadiation, ModalRadiationLoading,
    ROOM_TEMPERATURE_K, STANDARD_ATMOSPHERE_PA, StructuralModes, radiation_sample,
    with_radiation_cx,
};
use fs_material::gas::GasSpec;

fn gas() -> GasState {
    GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        ROOM_TEMPERATURE_K,
        STANDARD_ATMOSPHERE_PA,
    )
    .unwrap()
}
fn basis(gas: &GasState) -> CavityCompartmentBasis {
    let local = |volume| CavityModes {
        omegas: vec![0.0],
        lambdas: vec![volume],
        interface: vec![vec![1.0]],
        loss_factor: 0.0,
        rho0: gas.density,
        c0: gas.sound_speed,
    };
    CavityCompartmentBasis::try_new(&[local(0.004), local(0.006)], 4).unwrap()
}
fn model(basis: &CavityCompartmentBasis, coupling: Vec<f64>) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![400.0],
            shapes: vec![vec![1.0, -0.8]],
            loss_factor: 0.02,
        },
        basis.modes(),
        coupling,
        None,
    )
    .unwrap()
}
fn controls() -> DistributedOpeningControls {
    DistributedOpeningControls {
        maximum_isothermal_ka: 0.5,
        maximum_composition_condition: 1e10,
        boundary_relative_tolerance: 2e-10,
    }
}
fn link() -> CavityDuctConnection {
    CavityDuctConnection {
        left_compartment: 0,
        right_compartment: 1,
        left_pressure_projection: vec![1.0],
        right_pressure_projection: vec![1.0],
        radius_m: 0.003,
        length_m: 0.4,
        end_corrections: [
            Some(fs_phs::MouthFlange::Flanged),
            Some(fs_phs::MouthFlange::Unflanged),
        ],
    }
}
fn close(a: C64, b: C64) {
    assert!((a - b).abs() < 8e-10 * b.abs().max(1e-30), "{a:?} != {b:?}");
}
fn state(gas: GasState, scale: f64) -> RelaxationState {
    use fs_material::relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxingMedium, SpeciesRelaxation, VibrationalMode,
    };
    let id = fs_thermochem::SpeciesId::new("synthetic").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k: 1300.0,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 101325.0 * 0.002 * scale,
                t_ref_k: 293.15,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap()
}

#[test]
fn two_compartments_keep_native_internal_ducts_and_both_joint_exterior_routes() {
    use fs_bem::helmholtz::rayleigh_modal::AreaSample;
    use fs_bem::modal_loading::{LoadingRoute, SurfaceModeTraces, modal_radiation_loading};
    let gas = gas();
    let basis = basis(&gas);
    let frequencies = [250.0, 700.0, 1200.0];
    let outer = [
        CylindricalCavityNeck {
            radius_m: 0.004,
            length_m: 0.35,
            inner_end: None,
        },
        CylindricalCavityNeck {
            radius_m: 0.005,
            length_m: 0.55,
            inner_end: None,
        },
    ];
    let a = core::f64::consts::PI * outer[0].radius_m * outer[0].radius_m;
    let b = core::f64::consts::PI * outer[1].radius_m * outer[1].radius_m;
    let samples = [
        AreaSample {
            position_m: [0.0, 0.0],
            area_m2: 0.01,
            mode_shape: vec![1.0, 0.0, 0.0],
        },
        AreaSample {
            position_m: [0.12, 0.02],
            area_m2: a,
            mode_shape: vec![0.0, 1.0 / a, 0.0],
        },
        AreaSample {
            position_m: [0.2, 0.04],
            area_m2: b,
            mode_shape: vec![0.0, 0.0, 1.0 / b],
        },
    ];
    let bake = modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples: &samples },
        LoadingRoute::PlanarBaffledRayleigh,
        &frequencies,
        fs_bem::helmholtz::Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        },
        fs_blake3::hash_domain(
            "connected-cavity-fixture",
            b"two compartments and one joint exterior",
        ),
        12,
        || false,
    )
    .unwrap();
    let (inside, external, compact) = with_radiation_cx(|cx, _| {
        let inside = InternalCavityDucts::from_cylindrical_ducts(
            basis.clone(),
            vec![link()],
            &gas,
            &frequencies,
            controls(),
            cx,
        )
        .unwrap();
        let loading = ModalRadiationLoading::from_bem_bake(&bake, cx).unwrap();
        let projection = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let external = DistributedCavityOpeningRadiation::from_cylindrical_ducts(
            1,
            projection.clone(),
            &outer,
            &gas,
            loading.clone(),
            controls(),
            cx,
        )
        .unwrap();
        let mut short = outer;
        short[0].length_m = 0.01;
        short[1].length_m = 0.01;
        let compact = CavityOpeningRadiation::from_cylindrical_necks(
            1, projection, &short, &gas, loading, 0.25, cx,
        )
        .unwrap();
        (inside, external, compact)
    });
    let forward = model(&basis, vec![0.01, -0.008])
        .with_internal_cavity_ducts(inside.clone())
        .unwrap()
        .with_distributed_cavity_openings(external.clone())
        .unwrap();
    let reverse = model(&basis, vec![0.01, -0.008])
        .with_distributed_cavity_openings(external)
        .unwrap()
        .with_internal_cavity_ducts(inside.clone())
        .unwrap();
    let compact = model(&basis, vec![0.01, -0.008])
        .with_cavity_openings(compact)
        .unwrap()
        .with_internal_cavity_ducts(inside)
        .unwrap();
    let force = [C64::new(0.7, 0.3)];
    let source = [C64::new(1e-6, -2e-7), C64::new(-2e-7, 1e-7)];
    for omega in frequencies {
        let response = forward.frf_with_sources(omega, &force, &source).unwrap();
        let other = reverse.frf_with_sources(omega, &force, &source).unwrap();
        assert_eq!(response.b, other.b);
        assert_eq!(response.a, other.a);
        let inside = response.internal_ducts.as_ref().unwrap();
        assert!(inside.dissipation_w > 0.0);
        let physical = 0.5
            * (inside.left_pressure_pa[0] * inside.left_volume_flow_m3_s[0].conj()
                - inside.right_pressure_pa[0] * inside.right_volume_flow_m3_s[0].conj())
            .re;
        assert!(
            (physical - inside.dissipation_w).abs() < 2e-10 * inside.duct_power[0].incident_power_w
        );
        assert_eq!(response.power.internal_ducts, inside.dissipation_w);
        assert!(response.power.balance_residual.abs() < 2e-9);
        assert!(response.power.interface_residual.abs() < 2e-9);
        let observed = forward
            .frf_with_sources_and_convergence(omega, &force, &source, 1e-11)
            .unwrap();
        assert_eq!(observed.mode_counts, (1, 2, 1, 2));
        assert_eq!(observed.internal_duct_flow_delta, 0.0);
        assert!(
            forward
                .frf_truncated_with_sources(omega, &force, &source, 1, 1)
                .is_err()
        );
        let short = compact.frf_with_sources(omega, &force, &source).unwrap();
        assert!(short.power.balance_residual.abs() < 2e-9);
        assert!(short.power.interface_residual.abs() < 2e-9);
    }
    assert!(forward.frf(701.0, &force).is_err());
}

#[test]
fn molecular_internal_graph_preserves_the_original_cavity_law_in_both_attachment_orders() {
    let gas = gas();
    let basis = basis(&gas);
    let material = state(gas, 1.0);
    let different = state(gas, 2.0);
    let omega = 600.0;
    let mut geometry = link();
    geometry.radius_m = 0.0008;
    geometry.length_m = 0.3;
    geometry.end_corrections = [None, None];
    let wall = MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [100.0, 1200.0],
        maximum_relative_error: 0.02,
    };
    let graph = with_radiation_cx(|cx, _| {
        InternalCavityDucts::from_molecular_cylindrical_ducts(
            basis.clone(),
            vec![geometry.clone()],
            material.clone(),
            Some(wall),
            &[omega],
            controls(),
            cx,
        )
        .unwrap()
    });
    let attached = model(&basis, vec![0.0, 0.0])
        .with_internal_cavity_ducts(graph.clone())
        .unwrap();
    assert_eq!(attached.cavity_relaxation.as_ref(), Some(&material));
    assert!(
        attached
            .clone()
            .with_cavity_relaxation(material.clone())
            .is_ok()
    );
    assert!(
        attached
            .clone()
            .with_cavity_relaxation(different.clone())
            .is_err()
    );
    assert!(
        model(&basis, vec![0.0, 0.0])
            .with_cavity_relaxation(different)
            .unwrap()
            .with_internal_cavity_ducts(graph.clone())
            .is_err()
    );
    let before = model(&basis, vec![0.0, 0.0])
        .with_cavity_relaxation(material.clone())
        .unwrap()
        .with_internal_cavity_ducts(graph.clone())
        .unwrap();
    let source = [C64::new(1e-7, 3e-8), C64::new(-2e-8, 1e-8)];
    let response = attached
        .frf_with_sources(omega, &[C64::ZERO], &source)
        .unwrap();
    assert_eq!(
        response.a,
        before
            .frf_with_sources(omega, &[C64::ZERO], &source)
            .unwrap()
            .a
    );
    let original = graph.two_ports(omega).unwrap()[0].abcd().unwrap();
    let beta = material
        .caloric_response(omega)
        .unwrap()
        .bulk_modulus_pa
        .recip();
    let (cl, cr) = (beta.scale(0.004), beta.scale(0.006));
    let iw = C64::new(0.0, omega);
    let aa = original[0] - iw * cr * original[1];
    let bb = original[2] - iw * cr * original[3];
    let pr = (source[0] + (original[3] - iw * cl * original[1]) * source[1]) / (bb - iw * cl * aa);
    let qr = -source[1] - iw * cr * pr;
    close(response.a[1], pr);
    close(response.a[0], original[0] * pr + original[1] * qr);
    close(
        response
            .internal_ducts
            .as_ref()
            .unwrap()
            .right_volume_flow_m3_s[0],
        qr,
    );
    assert!(response.power.cavity > 0.0 && response.power.internal_ducts > 0.0);
    assert!(response.power.balance_residual.abs() < 3e-9);

    let (equilibrium, molecular_exterior, equilibrium_exterior) = with_radiation_cx(|cx, _| {
        let equilibrium = InternalCavityDucts::from_cylindrical_ducts(
            basis.clone(),
            vec![geometry.clone()],
            &gas,
            &[omega],
            controls(),
            cx,
        )
        .unwrap();
        let loading = ModalRadiationLoading::try_new(
            omega,
            2,
            vec![C64::ONE, C64::ZERO, C64::ZERO, C64::new(100.0, -200.0)],
            cx,
        )
        .unwrap();
        let mouth = [CylindricalCavityNeck {
            radius_m: 0.0008,
            length_m: 0.2,
            inner_end: None,
        }];
        let molecular = DistributedCavityOpeningRadiation::from_molecular_cylindrical_ducts(
            1,
            vec![vec![1.0], vec![0.0]],
            &mouth,
            material.clone(),
            Some(wall),
            loading.clone(),
            controls(),
            cx,
        )
        .unwrap();
        let ordinary = DistributedCavityOpeningRadiation::from_cylindrical_ducts(
            1,
            vec![vec![1.0], vec![0.0]],
            &mouth,
            &gas,
            loading,
            controls(),
            cx,
        )
        .unwrap();
        (equilibrium, molecular, ordinary)
    });
    assert!(
        attached
            .clone()
            .with_distributed_cavity_openings(molecular_exterior.clone())
            .is_ok()
    );
    assert!(
        model(&basis, vec![0.0, 0.0])
            .with_distributed_cavity_openings(molecular_exterior.clone())
            .unwrap()
            .with_internal_cavity_ducts(graph.clone())
            .is_ok()
    );
    assert!(
        attached
            .with_distributed_cavity_openings(equilibrium_exterior.clone())
            .is_err()
    );
    assert!(
        model(&basis, vec![0.0, 0.0])
            .with_distributed_cavity_openings(equilibrium_exterior)
            .unwrap()
            .with_internal_cavity_ducts(graph)
            .is_err()
    );
    assert!(
        model(&basis, vec![0.0, 0.0])
            .with_cavity_relaxation(material.clone())
            .unwrap()
            .with_internal_cavity_ducts(equilibrium.clone())
            .is_err()
    );
    assert!(
        model(&basis, vec![0.0, 0.0])
            .with_internal_cavity_ducts(equilibrium.clone())
            .unwrap()
            .with_cavity_relaxation(material)
            .is_err()
    );
    assert!(
        model(&basis, vec![0.0, 0.0])
            .with_internal_cavity_ducts(equilibrium)
            .unwrap()
            .with_distributed_cavity_openings(molecular_exterior)
            .is_err()
    );
}

#[test]
fn original_face_support_geometry_and_frequency_bounds_are_not_silently_replaced() {
    let gas = gas();
    let basis = basis(&gas);
    with_radiation_cx(|cx, gate| {
        let make = |links, frequencies: &[f64]| {
            InternalCavityDucts::from_cylindrical_ducts(
                basis.clone(),
                links,
                &gas,
                frequencies,
                controls(),
                cx,
            )
        };
        let mut invalid = link();
        invalid.right_compartment = 0;
        assert!(make(vec![invalid], &[600.0]).is_err());
        let mut wrong = link();
        wrong.left_pressure_projection[0] = 0.0;
        assert!(make(vec![wrong], &[600.0]).is_err());
        assert!(make(vec![link()], &[600.0, 500.0]).is_err());
        let original = make(vec![link()], &[600.0]).unwrap();
        let mut changed = model(&basis, vec![0.0, 0.0]);
        changed.lambdas[1] *= 2.0;
        assert!(changed.with_internal_cavity_ducts(original).is_err());
        gate.request();
        assert!(matches!(
            make(vec![link()], &[600.0]),
            Err(VibroError::Cancelled)
        ));
    });
}
