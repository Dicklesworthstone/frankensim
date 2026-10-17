//! IRP M6.3-W1 (bead pexfp.182): representation-aware import. Three generated
//! fixtures (closed box, midsurface plate, box with one side hole) and one
//! typed refusal per declared fact. Under `--nocapture` each admitted fixture
//! prints its receipt as one JSON line.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::quarantine::{Quarantined, quarantine};
use fs_io::{
    DeclaredLoop, Frame, LengthUnit, MAX_DECLARED_LOOPS, MaterialRef, PartMap, PartSelector,
    ReferenceState, Representation, RepresentationDeclaration, RepresentationReceipt,
    RepresentationRefusal, ThicknessField, admit_representation, promote,
};
use fs_rep_mesh::Soup;
use std::collections::BTreeMap;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x0182,
                kernel_id: 63,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn quarantined(soup: &Soup) -> Quarantined<Soup> {
    quarantine(soup.clone(), "generated", format!("{soup:?}").as_bytes())
}

/// The admitted soup, compared bitwise (`Soup` itself has no `PartialEq`).
#[derive(Debug)]
struct Admitted(Soup);

impl PartialEq for Admitted {
    fn eq(&self, other: &Self) -> bool {
        let bits = |p: &Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        self.0.triangles == other.0.triangles
            && self.0.positions.len() == other.0.positions.len()
            && self
                .0
                .positions
                .iter()
                .zip(&other.0.positions)
                .all(|(a, b)| bits(a) == bits(b))
    }
}

impl std::ops::Deref for Admitted {
    type Target = Soup;
    fn deref(&self) -> &Soup {
        &self.0
    }
}

fn admit(
    soup: &Soup,
    d: &RepresentationDeclaration,
) -> Result<(Admitted, RepresentationReceipt), RepresentationRefusal> {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        admit_representation(quarantined(soup), d, cx).map(|(e, r)| (Admitted(e.value), r))
    })
}

/// Surface of an axis-aligned cube `[0, side]^3`, each face a `k x k` grid,
/// outward-oriented; `hole` drops the centre cell of the `+x` face (`k` odd).
fn grid_box(k: u32, side: f64, hole: bool) -> Soup {
    let mut index: BTreeMap<[u32; 3], u32> = BTreeMap::new();
    let mut positions = Vec::new();
    let mut vertex = |p: [u32; 3]| -> u32 {
        *index.entry(p).or_insert_with(|| {
            positions.push(Point3::new(
                side * f64::from(p[0]) / f64::from(k),
                side * f64::from(p[1]) / f64::from(k),
                side * f64::from(p[2]) / f64::from(k),
            ));
            u32::try_from(positions.len() - 1).expect("fixture size")
        })
    };
    let mut quads = Vec::new();
    for axis in 0..3usize {
        for (level, outward) in [(0, -1.0), (k, 1.0)] {
            let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
            for i in 0..k {
                for j in 0..k {
                    if hole && axis == 0 && level == k && i == k / 2 && j == k / 2 {
                        continue;
                    }
                    let corner = |di: u32, dj: u32| {
                        let mut p = [0u32; 3];
                        p[axis] = level;
                        p[u] = i + di;
                        p[w] = j + dj;
                        p
                    };
                    quads.push((
                        [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)],
                        axis,
                        outward,
                    ));
                }
            }
        }
    }
    let indexed: Vec<_> = quads
        .into_iter()
        .map(|(quad, axis, outward)| (quad.map(&mut vertex), axis, outward))
        .collect();
    let mut triangles = Vec::new();
    for (ids, axis, outward) in indexed {
        for mut t in [[ids[0], ids[1], ids[2]], [ids[0], ids[2], ids[3]]] {
            // Orient outward along `axis`.
            if normal(&positions, t)[axis] * outward < 0.0 {
                t.swap(1, 2);
            }
            triangles.push(t);
        }
    }
    Soup {
        positions,
        triangles,
    }
}

/// Edges used by exactly one triangle.
fn open_edge_count(soup: &Soup) -> usize {
    let mut uses: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for t in &soup.triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *uses.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    uses.values().filter(|&&n| n == 1).count()
}

fn normal(positions: &[Point3], t: [u32; 3]) -> [f64; 3] {
    let [a, b, c] = t.map(|i| positions[i as usize]);
    let (u, v) = (
        [b.x - a.x, b.y - a.y, b.z - a.z],
        [c.x - a.x, c.y - a.y, c.z - a.z],
    );
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

/// A flat `n x n`-cell plate on `z = 0`, side `size`, oriented `+z`; `hole`
/// drops the centre cell (`n` odd).
fn plate(n: u32, size: f64, hole: bool) -> Soup {
    let id = |i: u32, j: u32| i * (n + 1) + j;
    let mut positions = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            positions.push(Point3::new(
                size * f64::from(i) / f64::from(n),
                size * f64::from(j) / f64::from(n),
                0.0,
            ));
        }
    }
    let mut triangles = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if hole && i == n / 2 && j == n / 2 {
                continue;
            }
            triangles.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            triangles.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    Soup {
        positions,
        triangles,
    }
}

fn whole(soup: &Soup, material: MaterialRef) -> Vec<PartMap> {
    vec![PartMap {
        selector: PartSelector::Range {
            start: 0,
            end: u32::try_from(soup.triangles.len()).expect("fixture size"),
        },
        part_token: "body".into(),
        material,
    }]
}

fn spruce() -> MaterialRef {
    MaterialRef::Card {
        pack_id: "fixture-pack".into(),
        card_id: "spruce".into(),
    }
}

/// M8-W1 ledger rows of the three generated fixtures (`data/instruments/ledger.tsv`).
const FIXTURE_ROWS: [&str; 3] = [
    "repr-fixture-closed-box",
    "repr-fixture-midsurface-plate",
    "repr-fixture-opened-box",
];

fn declaration(soup: &Soup, representation: Representation) -> RepresentationDeclaration {
    let ledger_row = match representation {
        Representation::ClosedSolid => FIXTURE_ROWS[0],
        Representation::Midsurface { .. } => FIXTURE_ROWS[1],
        Representation::OpenedSolid { .. } => FIXTURE_ROWS[2],
    };
    RepresentationDeclaration {
        representation,
        unit: Some(LengthUnit {
            label: "m".into(),
            metres_per_unit: 1.0,
        }),
        frame: Frame::IDENTITY,
        reference_state: Some(ReferenceState::Unloaded),
        parts: whole(soup, spruce()),
        ledger_row: ledger_row.into(),
        tolerance_m: 1.0e-9,
    }
}

/// The committed M8-W1 row for `id`, split into its 14 fields.
fn ledger_row(id: &str) -> Option<Vec<&'static str>> {
    include_str!("../../../data/instruments/ledger.tsv")
        .lines()
        .skip(2)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|fields| fields.first() == Some(&id))
}

#[test]
fn geometry_rows_resolve_in_the_instrument_ledger() {
    for id in FIXTURE_ROWS {
        let fields = ledger_row(id).unwrap_or_else(|| panic!("{id} is not in the M8-W1 ledger"));
        assert_eq!(fields.len(), 14, "{id}");
        assert_eq!(
            (fields[1], fields[2], fields[10]),
            ("M6.3-W1", "geometry", "fixture-only")
        );
    }
    assert!(ledger_row("repr-fixture-unledgered").is_none());
    let soup = grid_box(2, 0.1, false);
    let (_, receipt) = admit(&soup, &declaration(&soup, Representation::ClosedSolid))
        .expect("ledgered box admits");
    assert!(
        ledger_row(&receipt.ledger_row).is_some(),
        "receipt row resolves"
    );
}

fn loop_at(label: &str, probe: [f64; 3]) -> DeclaredLoop {
    DeclaredLoop {
        label: label.into(),
        probe,
    }
}

fn midsurface(soup: &Soup, edges: Vec<DeclaredLoop>) -> RepresentationDeclaration {
    declaration(
        soup,
        Representation::Midsurface {
            thickness: ThicknessField::PerVertex(vec![0.003; soup.positions.len()]),
            edges,
        },
    )
}

fn side_hole() -> DeclaredLoop {
    // A point on the hole's lower edge, centre cell of the +x face of a 3x3 box.
    loop_at("side-hole", [0.3, 0.15, 0.1])
}

fn report(fixture: &str, receipt: &RepresentationReceipt) {
    println!(
        "{{\"fixture\":\"{fixture}\",\"receipt\":{}}}",
        receipt.to_json()
    );
}

#[test]
fn g0_closed_solid_admits_watertight_manifold() {
    let soup = grid_box(2, 0.1, false);
    let (_, receipt) = admit(&soup, &declaration(&soup, Representation::ClosedSolid))
        .expect("watertight outward box admits");
    report("closed-box", &receipt);
    assert_eq!(receipt.representation, "closed-solid");
    assert_eq!(
        (
            receipt.open_edges,
            receipt.nonmanifold_edges,
            receipt.boundary_loops
        ),
        (0, 0, 0)
    );
    assert_eq!(receipt.repairs, 0);
    let volume = receipt
        .signed_volume_m3
        .expect("closed solids record volume");
    assert!(
        (volume - 1.0e-3).abs() <= 1.0e-12 * 1.0e-3,
        "volume {volume}"
    );
}

#[test]
fn g0_closed_solid_with_hole_refuses_without_repair() {
    let soup = grid_box(3, 0.3, true);
    assert_eq!(
        admit(&soup, &declaration(&soup, Representation::ClosedSolid)),
        Err(RepresentationRefusal::NotClosed { open_edges: 4 })
    );
}

#[test]
fn g0_flipped_triangle_and_inside_out_refuse() {
    let mut flipped = grid_box(2, 0.1, false);
    flipped.triangles[0].swap(1, 2);
    assert!(matches!(
        admit(
            &flipped,
            &declaration(&flipped, Representation::ClosedSolid)
        ),
        Err(RepresentationRefusal::Orientation {
            inconsistent_edges: 3
        })
    ));
    let mut inside_out = grid_box(2, 0.1, false);
    for t in &mut inside_out.triangles {
        t.swap(1, 2);
    }
    assert!(matches!(
        admit(&inside_out, &declaration(&inside_out, Representation::ClosedSolid)),
        Err(RepresentationRefusal::InsideOut { signed_volume_m3 }) if signed_volume_m3 < 0.0
    ));
}

#[test]
fn g0_midsurface_declared_edges_and_thickness_admit() {
    let soup = plate(3, 0.3, false);
    let d = midsurface(&soup, vec![loop_at("rim", [0.15, 0.0, 0.0])]);
    let (_, receipt) = admit(&soup, &d).expect("declared rim and thickness admit");
    report("midsurface-plate", &receipt);
    assert_eq!(receipt.representation, "midsurface");
    assert_eq!((receipt.boundary_loops, receipt.nonmanifold_edges), (1, 0));
    assert_eq!(receipt.declared_loops.len(), 1);
    assert_eq!(receipt.declared_loops[0].label, "rim");
    assert_eq!(receipt.declared_loops[0].vertices.len(), 12);
    assert_eq!(receipt.open_edges, 12);
    assert_eq!(receipt.thickness_layout, Some("per-vertex"));
    assert_eq!(receipt.thickness_range_m, Some((0.003, 0.003)));
}

#[test]
fn g0_midsurface_undeclared_loop_refuses() {
    let soup = plate(3, 0.3, true);
    let d = midsurface(&soup, vec![loop_at("rim", [0.15, 0.0, 0.0])]);
    assert!(matches!(
        admit(&soup, &d),
        Err(RepresentationRefusal::UndeclaredLoop { nearest_probe: Some(label), .. }) if label == "rim"
    ));
}

#[test]
fn g0_midsurface_missing_declared_loop_refuses() {
    let soup = plate(3, 0.3, false);
    let d = midsurface(
        &soup,
        vec![
            loop_at("rim", [0.15, 0.0, 0.0]),
            loop_at("ghost", [5.0, 5.0, 5.0]),
        ],
    );
    assert_eq!(
        admit(&soup, &d),
        Err(RepresentationRefusal::MissingDeclaredLoop {
            label: "ghost".into()
        })
    );
}

#[test]
fn g0_midsurface_nonmanifold_edge_refuses() {
    let mut soup = plate(3, 0.3, false);
    // A fin on the interior edge (5, 6) makes that edge used three times.
    soup.positions.push(Point3::new(0.15, 0.15, 0.1));
    let apex = u32::try_from(soup.positions.len() - 1).expect("fixture size");
    soup.triangles.push([5, 6, apex]);
    let d = midsurface(&soup, vec![loop_at("rim", [0.15, 0.0, 0.0])]);
    assert_eq!(
        admit(&soup, &d),
        Err(RepresentationRefusal::NonManifold { edges: 1 })
    );
}

#[test]
fn g0_thickness_length_cap_and_cap_plus_one() {
    let soup = plate(3, 0.3, false);
    let n = soup.positions.len();
    for len in [n - 1, n, n + 1] {
        let d = declaration(
            &soup,
            Representation::Midsurface {
                thickness: ThicknessField::PerVertex(vec![0.003; len]),
                edges: vec![loop_at("rim", [0.15, 0.0, 0.0])],
            },
        );
        let result = admit(&soup, &d);
        if len == n {
            assert!(result.is_ok(), "n admits: {result:?}");
        } else {
            assert_eq!(
                result,
                Err(RepresentationRefusal::ThicknessLength {
                    expected: n,
                    got: len
                })
            );
        }
    }
}

#[test]
fn g0_nonpositive_or_nonfinite_thickness_refuses() {
    let soup = plate(3, 0.3, false);
    for bad in [0.0, -0.001, f64::NAN, f64::INFINITY] {
        let mut values = vec![0.003; soup.triangles.len()];
        values[7] = bad;
        let d = declaration(
            &soup,
            Representation::Midsurface {
                thickness: ThicknessField::PerTriangle(values),
                edges: vec![loop_at("rim", [0.15, 0.0, 0.0])],
            },
        );
        assert!(matches!(
            admit(&soup, &d),
            Err(RepresentationRefusal::Thickness { index: 7, .. })
        ));
    }
}

#[test]
fn g0_real_opening_preserved_bitwise() {
    let soup = grid_box(3, 0.3, true);
    let d = declaration(
        &soup,
        Representation::OpenedSolid {
            openings: vec![side_hole()],
        },
    );
    let (admitted, receipt) = admit(&soup, &d).expect("declared opening admits");
    report("opened-box", &receipt);
    assert_eq!(
        (receipt.boundary_loops, receipt.open_edges, receipt.repairs),
        (1, 4, 0)
    );
    let hole = &receipt.declared_loops[0];
    assert_eq!(hole.label, "side-hole");
    assert_eq!(hole.vertices.len(), 4);
    assert_eq!(admitted.triangles, soup.triangles);
    for &v in &hole.vertices {
        let (before, after) = (soup.positions[v as usize], admitted.positions[v as usize]);
        assert_eq!(
            [before.x.to_bits(), before.y.to_bits(), before.z.to_bits()],
            [after.x.to_bits(), after.y.to_bits(), after.z.to_bits()],
            "opening vertex {v}"
        );
    }
}

#[test]
fn g0_legacy_promote_would_fill_the_opening() {
    // Liveness: the repairing path closes the same opening, so the guard matters.
    let soup = grid_box(3, 0.3, true);
    assert_eq!(open_edge_count(&soup), 4);
    let (promoted, _) = promote(quarantined(&soup), 8).expect("legacy promote repairs");
    assert!(promoted.value.triangles.len() > soup.triangles.len());
    assert_eq!(
        open_edge_count(&promoted.value),
        0,
        "the opening was fan-filled"
    );
}

#[test]
fn g0_length_unit_missing_refuses() {
    let soup = grid_box(1, 0.1, false);
    for unit in [
        None,
        Some(LengthUnit {
            label: String::new(),
            metres_per_unit: 1.0,
        }),
        Some(LengthUnit {
            label: "mm".into(),
            metres_per_unit: 0.0,
        }),
    ] {
        let mut d = declaration(&soup, Representation::ClosedSolid);
        d.unit = unit;
        assert_eq!(admit(&soup, &d), Err(RepresentationRefusal::UnitMissing));
    }
}

#[test]
fn g0_mm_scaling_recorded_once() {
    let soup = grid_box(2, 100.0, false);
    let mut d = declaration(&soup, Representation::ClosedSolid);
    d.unit = Some(LengthUnit {
        label: "mm".into(),
        metres_per_unit: 1.0e-3,
    });
    let (admitted, receipt) = admit(&soup, &d).expect("mm box admits");
    report("closed-box-mm", &receipt);
    assert_eq!(
        (receipt.unit.as_str(), receipt.metres_per_unit),
        ("mm", 1.0e-3)
    );
    let volume = receipt.signed_volume_m3.expect("volume");
    assert!(
        (volume - 1.0e-3).abs() <= 1.0e-12 * 1.0e-3,
        "volume {volume}"
    );
    for (src, out) in soup.positions.iter().zip(&admitted.positions) {
        assert_eq!(out.x.to_bits(), (src.x * 1.0e-3).to_bits());
        assert_eq!(out.z.to_bits(), (src.z * 1.0e-3).to_bits());
    }
}

#[test]
fn g0_left_handed_frame_refuses() {
    let soup = grid_box(1, 0.1, false);
    let mut d = declaration(&soup, Representation::ClosedSolid);
    d.frame.axes[2] = [0.0, 0.0, -1.0];
    assert!(matches!(
        admit(&soup, &d),
        Err(RepresentationRefusal::Frame { det, .. }) if (det + 1.0).abs() < 1e-15
    ));
    // A proper rotation (90 degrees about z) admits and keeps the volume.
    d.frame.axes = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let (_, receipt) = admit(&soup, &d).expect("rotation admits");
    let volume = receipt.signed_volume_m3.expect("volume");
    assert!((volume - 1.0e-3).abs() <= 1.0e-12 * 1.0e-3);
}

#[test]
fn g0_non_orthonormal_frame_refuses() {
    let soup = grid_box(1, 0.1, false);
    let mut d = declaration(&soup, Representation::ClosedSolid);
    d.frame.axes[0][0] = 1.0 + 1.0e-9;
    assert!(matches!(
        admit(&soup, &d),
        Err(RepresentationRefusal::Frame { orthogonality, .. }) if orthogonality > 1e-12
    ));
}

#[test]
fn g0_unmapped_and_doubly_mapped_triangles_refuse() {
    let soup = grid_box(1, 0.1, false);
    let n = u32::try_from(soup.triangles.len()).expect("fixture size");
    let part = |start, end, token: &str| PartMap {
        selector: PartSelector::Range { start, end },
        part_token: token.into(),
        material: spruce(),
    };
    let mut d = declaration(&soup, Representation::ClosedSolid);
    d.parts = vec![part(0, n - 1, "most")];
    assert_eq!(
        admit(&soup, &d),
        Err(RepresentationRefusal::Unmapped { triangles: 1 })
    );
    d.parts = vec![part(0, n, "all"), part(n - 2, n, "overlap")];
    assert_eq!(
        admit(&soup, &d),
        Err(RepresentationRefusal::DoublyMapped { triangles: 2 })
    );
    d.parts = vec![PartMap {
        selector: PartSelector::Triangles(vec![n]),
        part_token: "beyond".into(),
        material: spruce(),
    }];
    assert_eq!(
        admit(&soup, &d),
        Err(RepresentationRefusal::SelectorOutOfRange {
            part: "beyond".into()
        })
    );
}

#[test]
fn g0_usemtl_label_is_not_a_material() {
    let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\ng top\nusemtl spruce\n\
               f 1 3 2\nf 1 2 4\nf 2 3 4\nf 3 1 4\n";
    let document = fs_io::obj::read_obj_document(obj).expect("fixture OBJ parses");
    // A converter that forwards the region's usemtl label as the material is refused.
    let parts: Vec<PartMap> = document
        .regions
        .iter()
        .map(|region| PartMap {
            selector: PartSelector::Range {
                start: u32::try_from(region.triangles.start).expect("fixture size"),
                end: u32::try_from(region.triangles.end).expect("fixture size"),
            },
            part_token: format!("region-{}", region.triangles.start),
            material: region
                .material
                .clone()
                .map_or_else(spruce, MaterialRef::SourceLabel),
        })
        .collect();
    assert!(
        parts
            .iter()
            .any(|part| matches!(part.material, MaterialRef::SourceLabel(_))),
        "fixture OBJ must carry a usemtl region"
    );
    let mut d = declaration(&document.soup, Representation::ClosedSolid);
    d.parts = parts;
    assert_eq!(
        admit(&document.soup, &d),
        Err(RepresentationRefusal::MaterialFromLabel {
            label: "spruce".into()
        })
    );
    // The same regions with a declared card admit.
    for part in &mut d.parts {
        part.material = spruce();
    }
    assert!(admit(&document.soup, &d).is_ok());
}

#[test]
fn g0_reference_state_required() {
    let soup = grid_box(1, 0.1, false);
    let mut d = declaration(&soup, Representation::ClosedSolid);
    d.reference_state = None;
    assert_eq!(
        admit(&soup, &d),
        Err(RepresentationRefusal::ReferenceStateMissing)
    );
}

#[test]
fn g0_declared_loop_cap_and_cap_plus_one() {
    // MAX_DECLARED_LOOPS disjoint triangles: one free-edge loop each.
    let mut positions = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..MAX_DECLARED_LOOPS {
        let x = 2.0 * f64::from(u32::try_from(i).expect("fixture size"));
        let base = u32::try_from(positions.len()).expect("fixture size");
        positions.extend([
            Point3::new(x, 0.0, 0.0),
            Point3::new(x + 1.0, 0.0, 0.0),
            Point3::new(x, 1.0, 0.0),
        ]);
        triangles.push([base, base + 1, base + 2]);
    }
    let soup = Soup {
        positions,
        triangles,
    };
    let edges = |count: usize| -> Vec<DeclaredLoop> {
        (0..count)
            .map(|i| {
                let x = 2.0 * f64::from(u32::try_from(i).expect("fixture size"));
                loop_at(&format!("edge-{i}"), [x, 0.0, 0.0])
            })
            .collect()
    };
    let at_cap = declaration(
        &soup,
        Representation::Midsurface {
            thickness: ThicknessField::PerTriangle(vec![0.001; MAX_DECLARED_LOOPS]),
            edges: edges(MAX_DECLARED_LOOPS),
        },
    );
    let (_, receipt) = admit(&soup, &at_cap).expect("cap admits");
    assert_eq!(receipt.boundary_loops, MAX_DECLARED_LOOPS);
    let over = declaration(
        &soup,
        Representation::Midsurface {
            thickness: ThicknessField::PerTriangle(vec![0.001; MAX_DECLARED_LOOPS]),
            edges: edges(MAX_DECLARED_LOOPS + 1),
        },
    );
    assert_eq!(
        admit(&soup, &over),
        Err(RepresentationRefusal::Cap {
            what: "declared_loops",
            cap: MAX_DECLARED_LOOPS,
            got: MAX_DECLARED_LOOPS + 1
        })
    );
}

#[test]
fn g4_cancellation_mid_census_publishes_nothing() {
    let soup = grid_box(3, 0.3, true);
    let d = declaration(
        &soup,
        Representation::OpenedSolid {
            openings: vec![side_hole()],
        },
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    let result = with_cx(&gate, |cx| admit_representation(quarantined(&soup), &d, cx));
    assert!(matches!(
        result,
        Err(RepresentationRefusal::Census(refusal)) if refusal.cancelled
    ));
}

#[test]
fn g5_receipt_bytes_stable_across_runs() {
    let soup = grid_box(3, 0.3, true);
    let d = declaration(
        &soup,
        Representation::OpenedSolid {
            openings: vec![side_hole()],
        },
    );
    let (a, ra) = admit(&soup, &d).expect("first run");
    let (b, rb) = admit(&soup, &d).expect("second run");
    assert_eq!(ra.to_json(), rb.to_json());
    assert_eq!(a.triangles, b.triangles);
    for (p, q) in a.positions.iter().zip(&b.positions) {
        assert_eq!(
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()],
            [q.x.to_bits(), q.y.to_bits(), q.z.to_bits()]
        );
    }
}
