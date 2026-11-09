use super::*;
use crate::vibroacoustic::{
    CavityModes, ModalStructuralRelaxation, ROOM_TEMPERATURE_K, STANDARD_ATMOSPHERE_PA,
    StructuralModes, radiation_sample, with_radiation_cx,
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
fn neck() -> CylindricalCavityNeck {
    CylindricalCavityNeck {
        radius_m: 0.006,
        length_m: 1.1,
        inner_end: None,
    }
}
fn model(gas: &GasState, coupling: f64) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![400.0],
            shapes: vec![vec![1.0]],
            loss_factor: 0.0,
        },
        &CavityModes {
            omegas: vec![0.0],
            lambdas: vec![0.004],
            interface: vec![vec![1.0]],
            loss_factor: 0.0,
            rho0: gas.density,
            c0: gas.sound_speed,
        },
        vec![coupling],
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
fn close(a: C64, b: C64) {
    assert!((a - b).abs() < 4e-11 * b.abs().max(1e-30), "{a:?} != {b:?}");
}
fn molecular_state(gas: GasState, scale: f64) -> RelaxationState {
    use fs_material::relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxingMedium, SpeciesRelaxation, VibrationalMode,
    };
    let id = fs_thermochem::SpeciesId::new("opening-fixture").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
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
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
        .unwrap()
}

#[test]
fn finite_bessel_ducts_drive_joint_rayleigh_radiation_and_relaxing_structure() {
    use fs_bem::helmholtz::rayleigh_modal::AreaSample;
    use fs_bem::modal_loading::{LoadingRoute, SurfaceModeTraces, modal_radiation_loading};
    let gas = gas();
    let neck = neck();
    let area = core::f64::consts::PI * neck.radius_m * neck.radius_m;
    // Disjoint physical wall and mouth; mouth source is UNIT VOLUME FLOW.
    let cells = vec![
        AreaSample {
            position_m: [0.0, 0.0],
            area_m2: 0.01,
            mode_shape: vec![1.0, 0.0],
        },
        AreaSample {
            position_m: [0.12, 0.02],
            area_m2: area,
            mode_shape: vec![0.0, 1.0 / area],
        },
    ];
    let omegas = [300.0, 1000.0, 1800.0];
    let bake = modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples: &cells },
        LoadingRoute::PlanarBaffledRayleigh,
        &omegas,
        fs_bem::helmholtz::Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        },
        fs_blake3::hash_domain(
            "joint-cavity-opening-fixture",
            b"one real wall and one uniform mouth",
        ),
        12,
        || false,
    )
    .unwrap();
    let openings = with_radiation_cx(|cx, _| {
        DistributedCavityOpeningRadiation::from_cylindrical_ducts(
            1,
            vec![vec![1.0]],
            &[neck],
            &gas,
            ModalRadiationLoading::from_bem_bake(&bake, cx).unwrap(),
            controls(),
            cx,
        )
        .unwrap()
    });
    assert_eq!(
        openings.exterior_loading().source_declaration(),
        Some(bake.declaration_hash())
    );
    assert_ne!(
        openings.two_ports(omegas[0]).unwrap()[0].scattering,
        openings.two_ports(omegas[2]).unwrap()[0].scattering
    );
    let material = ModalStructuralRelaxation::try_new(
        1,
        &[fs_phs::RelaxationBranch {
            projection: vec![1.0],
            stiffness: 40_000.0,
            relaxation_time_s: 0.003,
        }],
        (100.0, 2000.0),
        1,
        &fs_exec::CancelGate::new_clock_free(),
    )
    .unwrap();
    let model = model(&gas, 0.01)
        .with_structural_relaxation(material)
        .unwrap()
        .with_distributed_cavity_openings(openings)
        .unwrap();
    for omega in omegas {
        let response = model
            .frf_with_sources(omega, &[C64::new(1.0, 0.2)], &[C64::new(5e-7, -2e-7)])
            .unwrap();
        let fields = response.openings.as_ref().unwrap();
        close(
            fields.radiation_port_velocities[0],
            C64::new(0.0, -omega) * response.b[0],
        );
        assert_eq!(
            fields.radiation_port_velocities[1],
            fields.volume_flow_m3_s[0]
        );
        assert_eq!(fields.distributed_duct_power.len(), 1);
        assert!(
            (fields.cavity_volume_flow_m3_s[0] - fields.volume_flow_m3_s[0]).abs()
                > 0.01 * fields.volume_flow_m3_s[0].abs()
        );
        close(
            fields.distributed_duct_power[0].inward_volume_flow_m3_s[0],
            fields.cavity_volume_flow_m3_s[0],
        );
        close(
            fields.distributed_duct_power[0].inward_volume_flow_m3_s[1],
            -fields.volume_flow_m3_s[0],
        );
        let z = bake.at_omega(omega).unwrap().impedance();
        let v = &fields.radiation_port_velocities;
        let effort = [z[0] * v[0] + z[1] * v[1], z[2] * v[0] + z[3] * v[1]];
        let independent_power = 0.5 * (effort[0] * v[0].conj() + effort[1] * v[1].conj()).re;
        assert!(
            (response.power.radiated - independent_power).abs() < 2e-12 * independent_power.abs()
        );
        assert!(response.power.structural_relaxation > 0.0);
        assert!(response.power.opening_necks > 0.0 && response.power.radiated > 0.0);
        assert!(response.power.balance_residual.abs() < 2e-11);
        assert!(response.power.interface_residual.abs() < 2e-11);
        let converged = model
            .frf_with_convergence(omega, &[C64::ONE], 1e-12)
            .unwrap();
        assert_eq!(converged.opening_flow_delta, 0.0);
    }
}

#[test]
fn original_molecular_duct_and_cavity_compression_share_one_physical_material() {
    let gas = gas();
    let state = molecular_state(gas, 1.0);
    let different = molecular_state(gas, 2.0);
    let omega = 600.0;
    let geometry = CylindricalCavityNeck {
        radius_m: 0.0008,
        length_m: 0.3,
        inner_end: None,
    };
    let wall = MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [100.0, 1200.0],
        maximum_relative_error: 0.02,
    };
    let (r, x) = fs_phs::compact_radiation_impedance(
        gas.density,
        gas.sound_speed,
        geometry.radius_m,
        omega,
        fs_phs::MouthFlange::Flanged,
    )
    .unwrap();
    let outer = C64::new(r, x);
    let exterior = radiation_sample(omega, 2, vec![C64::ZERO, C64::ZERO, C64::ZERO, outer]);
    let opening = with_radiation_cx(|cx, _| {
        DistributedCavityOpeningRadiation::from_molecular_cylindrical_ducts(
            1,
            vec![vec![1.0]],
            &[geometry],
            state.clone(),
            Some(wall),
            exterior.clone(),
            controls(),
            cx,
        )
        .unwrap()
    });
    let original = MolecularDuct::new_controlled(
        &fs_duct::Duct {
            segments: vec![fs_duct::Segment::Cylinder {
                radius: geometry.radius_m,
                length: geometry.length_m,
            }],
        },
        state.clone(),
        Some(wall),
        fs_duct::ConeSlices::Fixed(1),
        1,
        &mut || false,
    )
    .unwrap()
    .prepared_chain_controlled(omega, &mut || false)
    .unwrap()
    .blocks[0];
    assert_eq!(
        opening.two_ports(omega).unwrap()[0].scattering,
        original.scattering
    );
    assert_eq!(opening.molecular_state(), Some(&state));
    assert_eq!(
        opening.molecular_wall_resolution().unwrap().poles_per_field,
        8
    );
    assert!(opening.two_ports(omega + 1.0).is_err());
    let attached = model(&gas, 0.0)
        .with_distributed_cavity_openings(opening.clone())
        .unwrap();
    assert_eq!(attached.cavity_relaxation.as_ref(), Some(&state));
    assert!(
        attached
            .clone()
            .with_cavity_relaxation(state.clone())
            .is_ok()
    );
    assert!(
        attached
            .clone()
            .with_cavity_relaxation(different.clone())
            .is_err()
    );
    assert!(
        model(&gas, 0.0)
            .with_cavity_relaxation(different)
            .unwrap()
            .with_distributed_cavity_openings(opening.clone())
            .is_err()
    );

    // Independent finite-length pressure/flow transfer and original complex
    // cavity compliance: no equilibrium Cp or compact one-flow substitution.
    let m = original.abcd().unwrap();
    let a = m[0] * outer + m[1];
    let b = m[2] * outer + m[3];
    let compliance = state
        .caloric_response(omega)
        .unwrap()
        .bulk_modulus_pa
        .recip()
        .scale(0.004);
    let source = C64::new(1e-7, 3e-8);
    let qo = source / (b - C64::new(0.0, omega) * compliance * a);
    let response = attached
        .frf_with_sources(omega, &[C64::ZERO], &[source])
        .unwrap();
    let fields = response.openings.as_ref().unwrap();
    close(fields.volume_flow_m3_s[0], qo);
    close(fields.cavity_volume_flow_m3_s[0], b * qo);
    close(response.a[0], a * qo);
    assert!(response.power.cavity > 0.0 && response.power.opening_necks > 0.0);
    assert!(response.power.balance_residual.abs() < 5e-10);
    assert!(response.power.interface_residual.abs() < 5e-10);

    let equilibrium = with_radiation_cx(|cx, _| {
        DistributedCavityOpeningRadiation::from_cylindrical_ducts(
            1,
            vec![vec![1.0]],
            &[geometry],
            &gas,
            exterior,
            controls(),
            cx,
        )
        .unwrap()
    });
    assert!(
        model(&gas, 0.0)
            .with_cavity_relaxation(state.clone())
            .unwrap()
            .with_distributed_cavity_openings(equilibrium.clone())
            .is_err()
    );
    assert!(
        model(&gas, 0.0)
            .with_distributed_cavity_openings(equilibrium)
            .unwrap()
            .with_cavity_relaxation(state)
            .is_err()
    );
}

#[test]
fn physical_cavity_basis_transverse_domain_and_source_admission_remain_explicit() {
    let gas = gas();
    let omega = 600.0;
    let exterior = radiation_sample(
        omega,
        2,
        vec![C64::ONE, C64::ZERO, C64::ZERO, C64::new(100.0, -200.0)],
    );
    with_radiation_cx(|cx, gate| {
        let make = |geometry, projection, limits| {
            DistributedCavityOpeningRadiation::from_cylindrical_ducts(
                1,
                projection,
                &[geometry],
                &gas,
                exterior.clone(),
                limits,
                cx,
            )
        };
        let original = make(neck(), vec![vec![1.0]], controls()).unwrap();
        assert!(original.two_ports(omega + 1.0).is_err());
        let attached = model(&gas, 0.01)
            .with_distributed_cavity_openings(original.clone())
            .unwrap();
        assert!(
            attached
                .with_distributed_cavity_openings(original.clone())
                .is_err()
        );
        let mut resonator = model(&gas, 0.0);
        resonator.omegas_a[0] = 100.0;
        assert!(
            resonator
                .with_distributed_cavity_openings(original.clone())
                .is_err()
        );
        let mut duplicate = model(&gas, 0.0);
        duplicate.z_modal = Some(radiation_sample(omega, 1, vec![C64::ONE]));
        assert!(
            duplicate
                .with_distributed_cavity_openings(original)
                .is_err()
        );
        let mut entrance = neck();
        entrance.inner_end = Some(fs_phs::MouthFlange::Flanged);
        assert!(make(entrance, vec![vec![1.0], vec![0.2]], controls()).is_err());
        // A legitimate uniform-mode entrance is the native inertive two-port
        // composed with the distributed cylinder; it is not an outer end.
        assert!(make(entrance, vec![vec![1.0]], controls()).is_ok());
        let mut broad = entrance;
        broad.radius_m = 0.3 * (gas.specific_gas_constant * gas.temperature).sqrt() / omega;
        assert!(make(broad, vec![vec![1.0]], controls()).is_err());
        broad.inner_end = None;
        assert!(make(broad, vec![vec![1.0]], controls()).is_ok());
        let mut outside = controls();
        outside.maximum_isothermal_ka = 0.001;
        assert!(make(neck(), vec![vec![1.0]], outside).is_err());
        gate.request();
        assert!(matches!(
            make(neck(), vec![vec![1.0]], controls()),
            Err(VibroError::Cancelled)
        ));
    });
}
