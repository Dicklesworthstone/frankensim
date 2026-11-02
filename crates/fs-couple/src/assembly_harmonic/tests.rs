use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;
use fs_vfit::passivity::CertificateClass;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 65,
                kernel_id: 65,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 2,
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}

#[test]
fn g1_original_mimo_axis_mutual_entries_and_source_admission_survive_harmonic_access() {
    let model = MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real {
            pole: -400.0,
            residue: 0.0,
        }],
        residues: vec![vec![
            C64::from_re(40.0),
            C64::from_re(10.0),
            C64::from_re(10.0),
            C64::from_re(30.0),
        ]],
        d: vec![2.0, 0.2, 0.2, 1.5],
        e: vec![0.0; 4],
        reciprocal: true,
    };
    let mut source = PreparedMimoLoad::new(
        &model,
        1.0 / 48_000.0,
        MimoLoadLimits {
            band_rad_s: [1.0, 160_000.0],
            sampled_passivity_points: 33,
            max_states: 2,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    source.step_mean_into(&[0.7, -0.2], &mut [0.0; 2]).unwrap();
    let history = source.pole_states().to_vec();
    let continuous_class = source.reports().0.class;
    assert_eq!(source.reports().1.class, CertificateClass::GridOnly);
    with_cx(&CancelGate::new(), |cx| {
        for omega in [1.0, 400.0, 3000.0, 159_000.0] {
            let observed = source_impedance(&source, 2, omega, controls(), cx).unwrap();
            for j in 0..4 {
                let expected =
                    C64::from_re(model.d[j]) + model.residues[0][j] / C64::new(400.0, -omega);
                assert!((observed[j] - expected).abs() <= 1e-13 * expected.abs());
            }
            assert!(observed[0].im > 0.0);
        }
        for omega in [0.0, 0.5, 160_001.0, f64::NAN] {
            assert!(matches!(
                source_impedance(&source, 2, omega, controls(), cx),
                Err(AssemblyMimoHarmonicError::Source(_))
            ));
        }
        assert!(source_impedance(&source, 1, 400.0, controls(), cx).is_err());
        assert!(
            source_impedance(
                &source,
                2,
                400.0,
                HarmonicPortControls {
                    maximum_ports: 1,
                    ..controls()
                },
                cx
            )
            .is_err()
        );
        assert!(
            source_impedance(
                &source,
                2,
                400.0,
                HarmonicPortControls {
                    maximum_payload_bytes: 1,
                    ..controls()
                },
                cx
            )
            .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                source_impedance(&source, 2, 400.0, controls(), cancelled),
                Err(AssemblyMimoHarmonicError::Cancelled)
            ))
        });
    });
    assert_eq!(source.pole_states(), history);
    assert_eq!(source.source_model(), &model);
    assert_eq!(source.reports().0.class, continuous_class);
    assert_eq!(source.reports().1.class, CertificateClass::GridOnly);
}
