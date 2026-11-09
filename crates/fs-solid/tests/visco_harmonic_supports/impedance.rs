//! Original intrinsic heat stays separate from moving-boundary port work.

use super::*;
use fs_modal::harmonic::HarmonicPortControls;

fn port_controls() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 2,
        maximum_projection_terms: 10_000,
        maximum_payload_bytes: 1 << 20,
    }
}

#[test]
fn material_base_excitation_retains_loaded_reaction_and_original_quadrature_heat() {
    scope(|cx, _| {
        let source = solid(cx);
        let map = [0.7, 1.2];
        let impedance = C64::new(0.8, -0.3);
        let force = [C64::new(0.0002, -0.0001), C64::new(0.0001, 0.0003)];
        let support = C64::new(0.00015, -0.00006);
        for omega in [0.2, 1.0, 5.0, 17.0] {
            let prepared = source
                .prepare_supported_harmonic_response_with_port_impedance(
                    &[0],
                    omega,
                    &map,
                    1,
                    &[impedance],
                    controls(),
                    port_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(prepared.represented_dofs(), &[0, 3]);
            let [a00, a01, a11] = pencil(omega);
            let coefficient = C64::new(impedance.im, -impedance.re).scale(omega);
            let loaded00 = a00 + coefficient.scale(map[0] * map[0]);
            let loaded01 = a01 + coefficient.scale(map[0] * map[1]);
            let loaded11 = a11 + coefficient.scale(map[1] * map[1]);
            let q = (force[1] - loaded01 * support) / loaded11;
            let reaction = loaded00 * support + loaded01 * q - force[0];
            let flow = C64::new(
                (support.scale(map[0]) + q.scale(map[1])).im,
                -(support.scale(map[0]) + q.scale(map[1])).re,
            )
            .scale(omega);
            let result = prepared.solve(&force, &[support], cx).unwrap();
            complex_near(result.displacement_m[1], q);
            complex_near(result.support_reaction_n[0], reaction);
            let load = result.port_load.as_ref().unwrap();
            complex_near(load.flow[0], flow);
            complex_near(load.effort[0], impedance * flow);
            let intrinsic = 0.5 * omega * (-a11.im) * (q - support).norm_sq();
            let absorbed = 0.5 * impedance.re * flow.norm_sq();
            near(result.average_material_heat_w, intrinsic, intrinsic);
            near(load.average_absorbed_power_w, absorbed, absorbed);
            near(result.power_balance_residual_w, 0.0, intrinsic + absorbed);
            let observed = source
                .prepare_harmonic_material_observer(omega, observation_budget(), cx)
                .unwrap()
                .observe(&result.displacement_m, cx)
                .unwrap();
            near(observed.mean_material_heat_w, intrinsic, intrinsic);
            assert!(absorbed > 0.0);
        }
    });
}

#[test]
fn loaded_material_source_strain_and_port_admission_remain_transactional() {
    scope(|cx, gate| {
        let source = solid(cx);
        let prepare = |omega, z| {
            source.prepare_supported_harmonic_response_with_port_impedance(
                &[0, 3],
                omega,
                &[0.7, 1.2],
                1,
                &[z],
                controls(),
                port_controls(),
                cx,
            )
        };
        assert!(prepare(0.05, C64::ONE).is_err());
        assert!(prepare(1.0, C64::from_re(-1.0)).is_err());
        let prepared = prepare(1.0, C64::new(0.8, -0.3)).unwrap();
        let force = [C64::ZERO; 2];
        assert!(
            prepared
                .solve(&force, &[C64::new(0.1, 0.0), C64::ZERO], cx)
                .is_err()
        );
        let motion = [C64::new(0.0001, 0.0001), C64::new(-0.0001, 0.0002)];
        let expected = prepared.solve(&force, &motion, cx).unwrap();
        assert!(
            expected
                .port_load
                .as_ref()
                .unwrap()
                .average_absorbed_power_w
                > 0.0
        );
        let retry = prepared.solve(&force, &motion, cx).unwrap();
        assert_eq!(expected.displacement_m, retry.displacement_m);
        gate.request();
        assert!(prepared.solve(&force, &motion, cx).is_err());
    });
}
