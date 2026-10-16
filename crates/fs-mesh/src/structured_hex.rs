//! Deterministic structured quadratic hexes and piecewise-linear profile
//! extrusions. This is a mapped-block generator, not the frontier hex mesher.
use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};

/// Displacement/geometry nodes retained by the structured mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HexFamily {
    Serendipity20,
    Lagrange27,
}
/// Corners, bottom/top perimeter edges, vertical edges, then face centres
/// (z-,z+,y-,x+,y+,x-) and body centre. Coordinates are in [-1,1]^3.
pub const HEX27_REFERENCE_NODES: [[i8; 3]; 27] = [
    [-1, -1, -1],
    [1, -1, -1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, 1],
    [-1, 1, 1],
    [0, -1, -1],
    [1, 0, -1],
    [0, 1, -1],
    [-1, 0, -1],
    [0, -1, 1],
    [1, 0, 1],
    [0, 1, 1],
    [-1, 0, 1],
    [-1, -1, 0],
    [1, -1, 0],
    [1, 1, 0],
    [-1, 1, 0],
    [0, 0, -1],
    [0, 0, 1],
    [0, -1, 0],
    [1, 0, 0],
    [0, 1, 0],
    [-1, 0, 0],
    [0, 0, 0],
];
pub const HEX20_REFERENCE_NODES: [[i8; 3]; 20] = {
    let mut a = [[0; 3]; 20];
    let mut i = 0;
    while i < 20 {
        a[i] = HEX27_REFERENCE_NODES[i];
        i += 1;
    }
    a
};
pub const HEX_EDGES: [[usize; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];
/// Gmsh ids 17/12 -> canonical local indices. The importer owns actual parsing.
/// Node numbering verified against Gmsh src/geo/MHexahedron.h edges_hexa /
/// faces_hexa (https://gmsh.info/doc/texinfo/#Node-ordering).
pub const GMSH20_TO_CANONICAL: [usize; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14,
];
pub const GMSH27_TO_CANONICAL: [usize; 27] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14, 20, 22, 25, 23, 24, 21,
    26,
];

#[derive(Debug, Clone, PartialEq)]
pub enum HexElements {
    Hex20(Vec<[usize; 20]>),
    Hex27(Vec<[usize; 27]>),
}
impl HexElements {
    pub fn len(&self) -> usize {
        match self {
            Self::Hex20(v) => v.len(),
            Self::Hex27(v) => v.len(),
        }
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn family(&self) -> HexFamily {
        match self {
            Self::Hex20(_) => HexFamily::Serendipity20,
            Self::Hex27(_) => HexFamily::Lagrange27,
        }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredHexMesh {
    pub nodes_m: Vec<[f64; 3]>,
    pub elements: HexElements,
}
#[derive(Debug, Clone, Copy)]
pub struct HexMeshBudget {
    pub maximum_nodes: usize,
    pub maximum_elements: usize,
}
/// Cross sections are rectangles normal to x. All five coordinates/lengths
/// are required SI data; no depth or centre is inferred from an instrument.
#[derive(Debug, Clone, Copy)]
pub struct ExtrusionStation {
    pub x_m: f64,
    pub y_center_m: f64,
    pub z_center_m: f64,
    pub width_m: f64,
    pub depth_m: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub enum StructuredHexError {
    InvalidInput(&'static str),
    Budget {
        what: &'static str,
        requested: usize,
        maximum: usize,
    },
    SizeOverflow,
    Cancelled,
}
impl std::fmt::Display for StructuredHexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "structured hex: {self:?}")
    }
}
impl std::error::Error for StructuredHexError {}
fn poll(cx: &Cx<'_>) -> Result<(), StructuredHexError> {
    cx.checkpoint().map_err(|_| StructuredHexError::Cancelled)
}
fn counts(n: [usize; 3], family: HexFamily, b: HexMeshBudget) -> Result<(), StructuredHexError> {
    if n.contains(&0) {
        return Err(StructuredHexError::InvalidInput(
            "positive cell counts required",
        ));
    }
    let mul = |a: usize, b: usize, c: usize| {
        a.checked_mul(b)
            .and_then(|v| v.checked_mul(c))
            .ok_or(StructuredHexError::SizeOverflow)
    };
    let plus = |v: usize| v.checked_add(1).ok_or(StructuredHexError::SizeOverflow);
    let [x, y, z] = n;
    let [xx, yy, zz] = [plus(x)?, plus(y)?, plus(z)?];
    let nodes = match family {
        HexFamily::Lagrange27 => mul(
            plus(x.checked_mul(2).ok_or(StructuredHexError::SizeOverflow)?)?,
            plus(y.checked_mul(2).ok_or(StructuredHexError::SizeOverflow)?)?,
            plus(z.checked_mul(2).ok_or(StructuredHexError::SizeOverflow)?)?,
        )?,
        HexFamily::Serendipity20 => mul(xx, yy, zz)?
            .checked_add(mul(x, yy, zz)?)
            .and_then(|v| v.checked_add(mul(xx, y, zz).ok()?))
            .and_then(|v| v.checked_add(mul(xx, yy, z).ok()?))
            .ok_or(StructuredHexError::SizeOverflow)?,
    };
    for (what, requested, maximum) in [
        ("nodes", nodes, b.maximum_nodes),
        ("elements", mul(x, y, z)?, b.maximum_elements),
    ] {
        if requested > maximum {
            return Err(StructuredHexError::Budget {
                what,
                requested,
                maximum,
            });
        }
    }
    Ok(())
}
/// A uniform box [0,Lx] x [0,Ly] x [0,Lz]. Node IDs follow sorted logical
/// half-grid coordinates; adjacent elements share every retained face node.
pub fn extruded_box_hex(
    cells: [usize; 3],
    extent_m: [f64; 3],
    family: HexFamily,
    budget: HexMeshBudget,
    cx: &Cx<'_>,
) -> Result<StructuredHexMesh, StructuredHexError> {
    poll(cx)?;
    counts(cells, family, budget)?;
    if extent_m.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(StructuredHexError::InvalidInput(
            "positive finite box lengths required",
        ));
    }
    let stations = (0..=cells[0])
        .map(|i| {
            poll(cx)?;
            Ok(ExtrusionStation {
                x_m: extent_m[0] * i as f64 / cells[0] as f64,
                y_center_m: extent_m[1] / 2.0,
                z_center_m: extent_m[2] / 2.0,
                width_m: extent_m[1],
                depth_m: extent_m[2],
            })
        })
        .collect::<Result<Vec<_>, StructuredHexError>>()?;
    profile_extruded_hex(&stations, [cells[1], cells[2]], family, budget, cx)
}
/// Interpolate the supplied profile linearly between every retained station.
/// Midside sections use that same interpolation, yielding trilinear cells.
pub fn profile_extruded_hex(
    stations: &[ExtrusionStation],
    cross_cells: [usize; 2],
    family: HexFamily,
    budget: HexMeshBudget,
    cx: &Cx<'_>,
) -> Result<StructuredHexMesh, StructuredHexError> {
    poll(cx)?;
    if stations.len() < 2 {
        return Err(StructuredHexError::InvalidInput(
            "two profile stations required",
        ));
    }
    let cells = [stations.len() - 1, cross_cells[0], cross_cells[1]];
    counts(cells, family, budget)?;
    for (i, s) in stations.iter().enumerate() {
        poll(cx)?;
        if [s.x_m, s.y_center_m, s.z_center_m, s.width_m, s.depth_m]
            .iter()
            .any(|x| !x.is_finite())
            || s.width_m <= 0.0
            || s.depth_m <= 0.0
            || (i > 0 && s.x_m <= stations[i - 1].x_m)
        {
            return Err(StructuredHexError::InvalidInput(
                "finite increasing stations and positive section sizes required",
            ));
        }
    }
    let local: &[[i8; 3]] = match family {
        HexFamily::Serendipity20 => &HEX20_REFERENCE_NODES,
        HexFamily::Lagrange27 => &HEX27_REFERENCE_NODES,
    };
    let mut keys = BTreeSet::new();
    let mut topology = Vec::new();
    for i in 0..cells[0] {
        for j in 0..cells[1] {
            for k in 0..cells[2] {
                poll(cx)?;
                let row = local
                    .iter()
                    .map(|p| {
                        [
                            2 * i + (p[0] + 1) as usize,
                            2 * j + (p[1] + 1) as usize,
                            2 * k + (p[2] + 1) as usize,
                        ]
                    })
                    .collect::<Vec<_>>();
                keys.extend(row.iter().copied());
                topology.push(row);
            }
        }
    }
    let mut ids = BTreeMap::new();
    let mut nodes = Vec::with_capacity(keys.len());
    for key in keys {
        poll(cx)?;
        let index = key[0] / 2;
        let s = if key[0] % 2 == 0 {
            stations[index]
        } else {
            let a = stations[index];
            let b = stations[index + 1];
            ExtrusionStation {
                x_m: f64::midpoint(a.x_m, b.x_m),
                y_center_m: f64::midpoint(a.y_center_m, b.y_center_m),
                z_center_m: f64::midpoint(a.z_center_m, b.z_center_m),
                width_m: f64::midpoint(a.width_m, b.width_m),
                depth_m: f64::midpoint(a.depth_m, b.depth_m),
            }
        };
        let p = [
            s.x_m,
            s.y_center_m + s.width_m * (key[1] as f64 / (2 * cells[1]) as f64 - 0.5),
            s.z_center_m + s.depth_m * (key[2] as f64 / (2 * cells[2]) as f64 - 0.5),
        ];
        if p.iter().any(|v| !v.is_finite()) {
            return Err(StructuredHexError::InvalidInput(
                "unrepresentable mapped coordinate",
            ));
        }
        ids.insert(key, nodes.len());
        nodes.push(p);
    }
    let elements = match family {
        HexFamily::Serendipity20 => HexElements::Hex20(
            topology
                .iter()
                .map(|v| {
                    poll(cx)?;
                    Ok(core::array::from_fn(|i| ids[&v[i]]))
                })
                .collect::<Result<Vec<_>, StructuredHexError>>()?,
        ),
        HexFamily::Lagrange27 => HexElements::Hex27(
            topology
                .iter()
                .map(|v| {
                    poll(cx)?;
                    Ok(core::array::from_fn(|i| ids[&v[i]]))
                })
                .collect::<Result<Vec<_>, StructuredHexError>>()?,
        ),
    };
    poll(cx)?;
    Ok(StructuredHexMesh {
        nodes_m: nodes,
        elements,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn with_cx(f: impl FnOnce(&Cx<'_>)) {
        let gate = fs_exec::CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|a| {
            f(&Cx::new(
                &gate,
                a,
                fs_exec::StreamKey {
                    seed: 59,
                    kernel_id: 59,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ))
        });
    }
    #[test]
    fn structured_hex_counts_caps_replay_and_shared_faces() {
        with_cx(|cx| {
            for (family, count) in [(HexFamily::Serendipity20, 32), (HexFamily::Lagrange27, 45)] {
                let budget = HexMeshBudget {
                    maximum_nodes: count,
                    maximum_elements: 2,
                };
                let a = extruded_box_hex([2, 1, 1], [2.0, 1.0, 1.0], family, budget, cx).unwrap();
                assert_eq!(a.nodes_m.len(), count);
                assert_eq!(
                    a,
                    extruded_box_hex([2, 1, 1], [2.0, 1.0, 1.0], family, budget, cx).unwrap()
                );
                assert!(
                    extruded_box_hex(
                        [2, 1, 1],
                        [2.0, 1.0, 1.0],
                        family,
                        HexMeshBudget {
                            maximum_nodes: count - 1,
                            ..budget
                        },
                        cx
                    )
                    .is_err()
                );
                assert!(
                    extruded_box_hex(
                        [2, 1, 1],
                        [2.0, 1.0, 1.0],
                        family,
                        HexMeshBudget {
                            maximum_elements: 1,
                            ..budget
                        },
                        cx
                    )
                    .is_err()
                );
                let shared = match &a.elements {
                    HexElements::Hex20(e) => e[0].iter().filter(|v| e[1].contains(v)).count(),
                    HexElements::Hex27(e) => e[0].iter().filter(|v| e[1].contains(v)).count(),
                };
                assert_eq!(
                    shared,
                    if family == HexFamily::Serendipity20 {
                        8
                    } else {
                        9
                    }
                );
            }
        });
    }
    #[test]
    fn structured_hex_profile_stations_and_refusals() {
        with_cx(|cx| {
            let budget = HexMeshBudget {
                maximum_nodes: 100,
                maximum_elements: 10,
            };
            let p = [
                ExtrusionStation {
                    x_m: 0.0,
                    y_center_m: 0.0,
                    z_center_m: 0.0,
                    width_m: 2.0,
                    depth_m: 1.0,
                },
                ExtrusionStation {
                    x_m: 0.5,
                    y_center_m: 0.1,
                    z_center_m: 0.05,
                    width_m: 1.0,
                    depth_m: 0.5,
                },
                ExtrusionStation {
                    x_m: 1.0,
                    y_center_m: 0.0,
                    z_center_m: 0.0,
                    width_m: 2.0,
                    depth_m: 1.0,
                },
            ];
            let m = profile_extruded_hex(&p, [1, 1], HexFamily::Lagrange27, budget, cx).unwrap();
            assert!(m.nodes_m.contains(&[0.5, 0.1, 0.05]));
            assert_eq!(m.nodes_m.len(), 45);
            for s in p {
                for y in [-0.5, 0.5] {
                    for z in [-0.5, 0.5] {
                        assert!(m.nodes_m.contains(&[
                            s.x_m,
                            s.y_center_m + y * s.width_m,
                            s.z_center_m + z * s.depth_m
                        ]));
                    }
                }
            }
            let mut bad = p;
            bad[1].depth_m = 0.0;
            assert!(profile_extruded_hex(&bad, [1, 1], HexFamily::Lagrange27, budget, cx).is_err());
            assert!(
                extruded_box_hex([0, 1, 1], [1.0; 3], HexFamily::Lagrange27, budget, cx).is_err()
            );
        });
    }
}
