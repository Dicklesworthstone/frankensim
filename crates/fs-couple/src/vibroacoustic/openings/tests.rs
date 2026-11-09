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
        radius_m: 0.012,
        length_m: 0.025,
        inner_end: Some(fs_phs::MouthFlange::Flanged),
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
fn opening(gas: &GasState, omega: f64, z: Vec<C64>) -> CavityOpeningRadiation {
    with_radiation_cx(|cx, _| {
        CavityOpeningRadiation::from_cylindrical_necks(
            1,
            vec![vec![1.0]],
            &[neck()],
            gas,
            ModalRadiationLoading::try_new(omega, 2, z, cx).unwrap(),
            0.25,
            cx,
        )
        .unwrap()
    })
}
fn close(a: C64, b: C64) {
    assert!((a - b).abs() < 4e-12 * b.abs().max(1e-30), "{a:?} != {b:?}");
}

#[test]
fn native_bessel_neck_and_outer_radiation_recover_the_physical_helmholtz_transfer() {
    let gas = gas();
    let neck = neck();
    let volume = 0.004;
    let compliance = volume / (gas.density * gas.sound_speed * gas.sound_speed);
    for omega in [90.0, 450.0, 900.0] {
        let (r, x) = fs_phs::compact_radiation_impedance(
            gas.density,
            gas.sound_speed,
            neck.radius_m,
            omega,
            fs_phs::MouthFlange::Flanged,
        )
        .unwrap();
        let outer = C64::new(r, x);
        let loading = opening(&gas, omega, vec![C64::ZERO, C64::ZERO, C64::ZERO, outer]);
        let internal = loading.neck_impedances(omega).unwrap()[0];
        let wave =
            fs_duct::segment_wave(&gas, neck.radius_m, omega, fs_duct::LossModel::Bessel).unwrap();
        let delta = fs_phs::neck_end_correction(neck.radius_m, neck.inner_end.unwrap()).unwrap();
        let area = core::f64::consts::PI * neck.radius_m * neck.radius_m;
        let original = (C64::new(0.0, -1.0) * wave.wavenumber * wave.characteristic_impedance)
            .scale(neck.length_m)
            + C64::new(0.0, -omega * gas.density * delta / area);
        close(internal, original);
        let model = model(&gas, 0.0).with_cavity_openings(loading).unwrap();
        let source = C64::new(1e-6, 2e-7);
        let response = model
            .frf_with_sources(omega, &[C64::ZERO], &[source])
            .unwrap();
        let total = internal + outer;
        let exact_q = source / (C64::ONE - C64::new(0.0, omega * compliance) * total);
        let openings = response.openings.as_ref().unwrap();
        close(openings.volume_flow_m3_s[0], exact_q);
        close(response.a[0], total * exact_q);
        close(
            openings.cavity_pressure_pa[0] - openings.exterior_pressure_pa[0],
            internal * exact_q,
        );
        assert!(response.power.opening_necks > 0.0 && response.power.radiated > 0.0);
        assert!(response.power.balance_residual.abs() < 2e-12);
        assert!(response.power.interface_residual.abs() < 2e-12);
        assert!(model.frf(omega + 1.0, &[C64::ZERO]).is_err());
    }
}

#[test]
fn native_joint_rayleigh_bake_drives_the_relaxing_structure_and_cavity_on_its_exact_grid() {
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
    let omegas = [150.0, 420.0, 700.0];
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
        CavityOpeningRadiation::from_cylindrical_necks(
            1,
            vec![vec![1.0]],
            &[neck],
            &gas,
            ModalRadiationLoading::from_bem_bake(&bake, cx).unwrap(),
            0.25,
            cx,
        )
        .unwrap()
    });
    assert_eq!(
        openings.exterior_loading().source_declaration(),
        Some(bake.declaration_hash())
    );
    assert_ne!(
        openings.neck_impedances(omegas[0]).unwrap(),
        openings.neck_impedances(omegas[2]).unwrap()
    );
    let material = ModalStructuralRelaxation::try_new(
        1,
        &[fs_phs::RelaxationBranch {
            projection: vec![1.0],
            stiffness: 40_000.0,
            relaxation_time_s: 0.003,
        }],
        (100.0, 800.0),
        1,
        &fs_exec::CancelGate::new_clock_free(),
    )
    .unwrap();
    let model = model(&gas, 0.01)
        .with_structural_relaxation(material)
        .unwrap()
        .with_cavity_openings(openings)
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
fn openings_refuse_double_neck_double_exterior_wrong_medium_and_noncompact_geometry() {
    let gas = gas();
    let omega = 300.0;
    let z = vec![C64::ONE, C64::ZERO, C64::ZERO, C64::new(10.0, -100.0)];
    let original = opening(&gas, omega, z.clone());
    let mut resonator = model(&gas, 0.01);
    resonator.omegas_a[0] = 200.0; // An already neck-loaded Helmholtz carrier.
    assert!(resonator.with_cavity_openings(original.clone()).is_err());
    let mut duplicated = model(&gas, 0.01);
    duplicated.z_modal = Some(radiation_sample(omega, 1, vec![C64::ONE]));
    assert!(duplicated.with_cavity_openings(original.clone()).is_err());
    let attached = model(&gas, 0.01)
        .with_cavity_openings(original.clone())
        .unwrap();
    assert!(attached.with_cavity_openings(original.clone()).is_err());
    let mut other = model(&gas, 0.01);
    other.medium_a.rho0 *= 2.0;
    assert!(other.with_cavity_openings(original).is_err());
    with_radiation_cx(|cx, gate| {
        let exterior = ModalRadiationLoading::try_new(omega, 2, z.clone(), cx).unwrap();
        let mut long = neck();
        long.length_m = 1.0;
        assert!(
            CavityOpeningRadiation::from_cylindrical_necks(
                1,
                vec![vec![1.0]],
                &[long],
                &gas,
                exterior.clone(),
                0.25,
                cx,
            )
            .is_err()
        );
        assert!(
            CavityOpeningRadiation::from_cylindrical_necks(
                1,
                vec![vec![f64::NAN]],
                &[neck()],
                &gas,
                exterior.clone(),
                0.25,
                cx,
            )
            .is_err()
        );
        assert!(
            CavityOpeningRadiation::from_cylindrical_necks(
                1,
                vec![vec![1.0], vec![-0.2]],
                &[neck()],
                &gas,
                exterior.clone(),
                0.25,
                cx,
            )
            .is_err()
        );
        let mut resolved = neck();
        resolved.inner_end = None;
        assert!(
            CavityOpeningRadiation::from_cylindrical_necks(
                1,
                vec![vec![1.0], vec![-0.2]],
                &[resolved],
                &gas,
                exterior.clone(),
                0.25,
                cx,
            )
            .is_ok()
        );
        gate.request();
        assert!(matches!(
            CavityOpeningRadiation::from_cylindrical_necks(
                1,
                vec![vec![1.0]],
                &[neck()],
                &gas,
                exterior,
                0.25,
                cx,
            ),
            Err(VibroError::Cancelled)
        ));
    });
}
