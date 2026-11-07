//! G0/G1/G3/G4: signed DKT surface stress, native source binding and rainflow.
#[path = "fatigue/fixture.rs"]
#[allow(dead_code)]
mod fixture;
use fixture::*;

use fs_material::fatigue::MeanStressCorrection;
use fs_plate::fatigue::{FatiguePlate, PlateFatigueError};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateSection, assemble};

fn options(tension: f64) -> AssemblyOptions {
    AssemblyOptions {
        pretension: tension,
        support: EdgeSupport::Clamped,
    }
}
// Independent quadratic Kirchhoff field: native beta = -grad(w).
fn field(mesh: &PlateMesh, k: [f64; 3]) -> Vec<f64> {
    mesh.nodes
        .iter()
        .flat_map(|&(x, y)| {
            [
                -0.5 * (k[0] * x * x + k[1] * y * y + k[2] * x * y),
                -k[0] * x - 0.5 * k[2] * y,
                -k[1] * y - 0.5 * k[2] * x,
            ]
        })
        .collect()
}

#[test]
fn g1_signed_faces_twist_and_pretension_match_plane_stress_equations() {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh {
        nodes: vec![(0.001, -0.002), (0.024, 0.001), (0.004, 0.019)],
        tris: vec![[0, 1, 2]],
    };
    let mut gauges = vec![];
    for point in 0..3 {
        for depth in [-H / 2.0, 0.0, H / 2.0] {
            for normal in [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]] {
                let mut g = gauge(&f, 0, depth, normal);
                g.surface_point = point;
                gauges.push(g);
            }
        }
    }
    let tension = 2e3;
    let plate = FatiguePlate::assemble(
        &mesh,
        vec![section(&f)],
        &[],
        &options(tension),
        gauges,
        config(100),
        &mut || false,
    )
    .unwrap();
    let k = [0.03, -0.02, 0.04];
    let stresses = plate.stresses(&field(&mesh, k), &mut || false).unwrap();
    for (i, actual) in stresses.iter().enumerate() {
        let g = plate.gauge(i).unwrap();
        let [nx, ny] = g.normal_in_mesh;
        let norm2 = nx * nx + ny * ny;
        let sx = E / (1.0 - NU * NU) * g.depth_m * (k[0] + NU * k[1]) + tension / H;
        let sy = E / (1.0 - NU * NU) * g.depth_m * (k[1] + NU * k[0]) + tension / H;
        let tau = E / (2.0 * (1.0 + NU)) * g.depth_m * k[2];
        close(
            *actual,
            (nx * nx * sx + ny * ny * sy + 2.0 * nx * ny * tau) / norm2,
            2e-13,
        );
    }
    // Exactly the original material section and support operator are used.
    let native = assemble(
        &mesh,
        &PlateSection::isotropic(E, NU, H, 7800.0).unwrap(),
        &[],
        &[],
        &options(tension),
    )
    .unwrap();
    assert_eq!(plate.model().k.to_dense(), native.k.to_dense());
    assert_eq!(plate.model().m.to_dense(), native.m.to_dense());
    let rigid: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|&(x, y)| [0.001 + 0.002 * x - 0.003 * y, 0.002, -0.003])
        .collect();
    for stress in plate.stresses(&rigid, &mut || false).unwrap() {
        close(stress, tension / H, 2e-12);
    }
}

#[test]
fn g1_repeated_bending_counts_signed_cycles_and_goodman_pretension() {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let mut g = gauge(&f, 0, H / 2.0, [1.0, 0.0]);
    let ultimate = 100e6;
    g.mean_stress = MeanStressCorrection::GoodmanTensionOnly {
        ultimate_strength_pa: ultimate,
    };
    let tension = 20e3;
    let plate = FatiguePlate::assemble(
        &mesh,
        vec![section(&f); 2],
        &[],
        &options(tension),
        vec![g],
        config(100),
        &mut || false,
    )
    .unwrap();
    let amplitude = 2e6;
    let curvature = amplitude * (1.0 - NU * NU) / (E * H / 2.0);
    let positive = field(&mesh, [curvature, 0.0, 0.0]);
    let negative = field(&mesh, [-curvature, 0.0, 0.0]);
    let mut observed = plate.observe(&negative, &mut || false).unwrap();
    let mut checkpoint = observed.clone();
    for _ in 0..3 {
        observed.push(&positive, &mut || false).unwrap();
        observed.push(&negative, &mut || false).unwrap();
    }
    for _ in 0..3 {
        checkpoint.push(&positive, &mut || false).unwrap();
        checkpoint.push(&negative, &mut || false).unwrap();
    }
    let before = observed.summaries();
    let preview = observed.preview_finish(&mut || false).unwrap();
    assert_eq!(observed.summaries(), before);
    let expected = 3.0 * amplitude / (1.0 - tension / H / ultimate) / LIFE_STRESS;
    close(preview[0].miner_damage, expected, 3e-13);
    let finished = observed.finish(&mut || false).unwrap();
    assert_eq!(preview, finished);
    assert_eq!(observed.finish(&mut || false).unwrap(), finished);
    assert_eq!(checkpoint.finish(&mut || false).unwrap(), finished);
    assert!(matches!(
        observed.push(&positive, &mut || false),
        Err(PlateFatigueError::Finished)
    ));
}

#[test]
fn g3_piecewise_material_stress_uses_its_own_section_and_source() {
    let a = fixture(E, 0.2);
    let b = fixture(E * 0.5, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let plate = FatiguePlate::assemble(
        &mesh,
        vec![section(&a), section(&b)],
        &[],
        &options(0.0),
        vec![
            gauge(&a, 0, H / 2.0, [1.0, 0.0]),
            gauge(&b, 1, H / 2.0, [1.0, 0.0]),
        ],
        config(10),
        &mut || false,
    )
    .unwrap();
    let stresses = plate
        .stresses(&field(&mesh, [0.01, 0.0, 0.0]), &mut || false)
        .unwrap();
    close(stresses[0], 2.0 * stresses[1], 1e-13);
    assert_eq!(
        plate.sections()[1].material.resolved().card_identity(),
        b.material.resolved().card_identity()
    );
    let wrong = fixture(E, 0.7);
    assert!(matches!(
        FatiguePlate::assemble(
            &mesh,
            vec![section(&a); 2],
            &[],
            &options(0.0),
            vec![gauge(&wrong, 0, 0.0, [1.0, 0.0])],
            config(10),
            &mut || false
        ),
        Err(PlateFatigueError::SourceMismatch { .. })
    ));
    assert!(matches!(
        FatiguePlate::assemble(
            &mesh,
            vec![section(&a); 2],
            &[],
            &options(0.0),
            vec![gauge(&b, 0, 0.0, [1.0, 0.0])],
            config(10),
            &mut || false
        ),
        Err(PlateFatigueError::SourceMismatch { .. })
    ));
}

#[test]
fn g4_late_gauge_refusal_and_cancellation_preserve_all_history() {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let mut weak = gauge(&f, 1, H / 2.0, [1.0, 0.0]);
    weak.mean_stress = MeanStressCorrection::GoodmanTensionOnly {
        ultimate_strength_pa: 1e6,
    };
    let plate = FatiguePlate::assemble(
        &mesh,
        vec![section(&f); 2],
        &[],
        &options(0.0),
        vec![gauge(&f, 0, H / 2.0, [1.0, 0.0]), weak],
        config(10),
        &mut || false,
    )
    .unwrap();
    let zero = vec![0.0; plate.model().free];
    let mut observed = plate.observe(&zero, &mut || false).unwrap();
    let old = observed.summaries();
    let q = field(&mesh, [0.02, 0.0, 0.0]);
    assert!(matches!(
        observed.push(&q, &mut || false),
        Err(PlateFatigueError::Gauge { gauge: 1, .. })
    ));
    assert_eq!(observed.summaries(), old);
    assert_eq!(observed.stresses_pa(), &[0.0, 0.0]);
    assert!(matches!(
        observed.push(&zero, &mut || true),
        Err(PlateFatigueError::Cancelled)
    ));
    assert_eq!(observed.summaries(), old);
    let huge = field(&mesh, [100.0, 0.0, 0.0]);
    assert!(observed.push(&huge, &mut || false).is_err());
    assert_eq!(observed.summaries(), old);
    assert!(observed.finish(&mut || true).is_err());
    assert!(!observed.is_finished());
    observed.push(&zero, &mut || false).unwrap();
    assert_eq!(observed.finish(&mut || false).unwrap()[0].miner_damage, 0.0);
}
