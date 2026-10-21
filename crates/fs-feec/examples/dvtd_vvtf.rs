//! Volume-velocity transfer function of a 3-D vocal tract from its tagged Gmsh volume mesh (music bead pexfp.211 / pexfp.93),
//! for comparison with the measured and FE-calculated transfer functions of the Dresden Vocal Tract Dataset (DVTD, CC0).
//!
//! usage: `dvtd_vvtf MSH ORDER F_MAX_HZ STATIC SCAN_MAX_HZ [nd]`
//! - MSH: an MSH 2.2 file whose triangles carry the physical tags 1 = glottis, 2 = wall, 3 = lips (DVTD's convention) and whose
//!   coordinates are millimetres;
//! - ORDER: element order 1..=4; F_MAX_HZ: top of the rigid-mode window; STATIC: 1 adds one static vector per port, 0 does not;
//!   SCAN_MAX_HZ: top of the response scan; `nd` selects nested-dissection ordering for the sparse factorizations.
//!
//! Boundary laws are INPUT laws, declared and cited: the glottis is a prescribed inward particle velocity; the wall has the real
//! specific admittance 1/500 (DVTD's FE used Z = 500 rho c, Fleischer et al. 2018); the lips have the radiation impedance
//! Z = rho c ((kr)^2 / (1 + (kr)^2) + j kr / (1 + kr)) with r = sqrt(A_lips / (2 pi)) (DVTD's FE, from Fleischer et al.); the gas is
//! c = 345 m/s, rho = 1.18 kg/m^3 (22 C). One JSON line is printed: counts, timings and the peaks of |H| above 3.

use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium, ReducedHelmholtz};
use fs_io::gmsh::{GmshElementType, parse_gmsh};
use fs_modal::SliceOptions;
use fs_rep_mesh::TetComplex;
use fs_sparse::direct::DirectOrdering;
use std::collections::BTreeMap;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 6 {
        eprintln!("usage: dvtd_vvtf MSH ORDER F_MAX_HZ STATIC SCAN_MAX_HZ [nd]");
        std::process::exit(64);
    }
    let order: usize = args[2].parse().expect("order");
    let f_max: f64 = args[3].parse().expect("f_max");
    let with_static = args[4] == "1";
    let scan_max: f64 = args[5].parse().expect("scan_max");
    let nested = args.get(6).is_some_and(|s| s == "nd");
    let t0 = Instant::now();
    let text = std::fs::read_to_string(&args[1]).expect("read msh");
    let (mesh, _receipt) = parse_gmsh(&text).expect("parse msh");
    let mesh = mesh.into_inner();
    let index: BTreeMap<u64, u32> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, n)| (n.id, u32::try_from(i).expect("node index")))
        .collect();
    let positions: Vec<[f64; 3]> = mesh
        .nodes
        .iter()
        .map(|n| [n.coords[0] * 1e-3, n.coords[1] * 1e-3, n.coords[2] * 1e-3])
        .collect();
    let mut tets = Vec::new();
    let mut faces: BTreeMap<i32, Vec<[u32; 3]>> = BTreeMap::new();
    for e in &mesh.elements {
        let ids: Vec<u32> = e.node_ids.iter().map(|id| index[id]).collect();
        match e.element_type {
            GmshElementType::Tetrahedron => tets.push([ids[0], ids[1], ids[2], ids[3]]),
            GmshElementType::Triangle => faces
                .entry(*e.tags.first().expect("physical tag"))
                .or_default()
                .push([ids[0], ids[1], ids[2]]),
            _ => {}
        }
    }
    let complex = TetComplex::from_tets(positions.len(), tets);
    let medium = Medium {
        density_kg_m3: 1.18,
        sound_speed_m_s: 345.0,
    };
    let fe = HelmholtzFe::new(&complex, &positions, order, medium).expect("fe");
    let t_asm = t0.elapsed().as_secs_f64();
    let glottis = &faces[&1];
    let wall = &faces[&2];
    let lips = &faces[&3];
    let opts = SliceOptions {
        ordering: if nested {
            DirectOrdering::NestedDissection
        } else {
            DirectOrdering::Amd
        },
        ..SliceOptions::default()
    };
    let t1 = Instant::now();
    let modes = fe.modes((-1.0, f_max), &opts).expect("modes");
    let t_modes = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let rom = ReducedHelmholtz::new(&fe, &modes, &[glottis, wall, lips], with_static).expect("rom");
    let t_rom = t2.elapsed().as_secs_f64();
    let (_, a_lips) = fe.boundary_load(lips).expect("lips area");
    let r_lips = (a_lips / (2.0 * core::f64::consts::PI)).sqrt();
    let t3 = Instant::now();
    let df = 2.0;
    let mut grid = Vec::new();
    let mut f = 100.0;
    while f <= scan_max {
        let k = 2.0 * core::f64::consts::PI * f / medium.sound_speed_m_s;
        let kr = k * r_lips;
        let (re, im) = (kr * kr / (1.0 + kr * kr), kr / (1.0 + kr));
        let d = re * re + im * im;
        let y_lips = [re / d, -im / d];
        let r = rom
            .response(f, 0, 1.0, &[[0.0, 0.0], [1.0 / 500.0, 0.0], y_lips])
            .expect("response");
        let h = r.transfer_to(2);
        grid.push((f, h[0].hypot(h[1])));
        f += df;
    }
    let t_scan = t3.elapsed().as_secs_f64();
    let mut peaks: Vec<(f64, f64)> = Vec::new();
    for i in 1..grid.len() - 1 {
        let (f0, m0) = grid[i];
        if m0 >= 3.0 && m0 > grid[i - 1].1 && m0 >= grid[i + 1].1 {
            let (a, b, c) = (grid[i - 1].1, m0, grid[i + 1].1);
            let denom = a - 2.0 * b + c;
            let shift = if denom.abs() > 0.0 {
                0.5 * (a - c) / denom
            } else {
                0.0
            };
            let fp = f0 + shift * df;
            if peaks.last().is_none_or(|p| (fp - p.0).abs() > 150.0) {
                peaks.push((fp, b));
            } else if peaks.last().is_some_and(|p| b > p.1) {
                *peaks.last_mut().expect("peak") = (fp, b);
            }
        }
    }
    let peak_json: Vec<String> = peaks
        .iter()
        .map(|(f, m)| format!("[{f:.2},{m:.3}]"))
        .collect();
    println!(
        "{{\"suite\":\"fs-feec/acoustics\",\"case\":\"dvtd_vvtf\",\"msh\":\"{}\",\"order\":{order},\"dofs\":{},\"vertices\":{},\"tets\":{},\
         \"modes\":{},\"basis\":{},\"static_vectors\":{with_static},\"f_max_hz\":{f_max},\"lips_area_mm2\":{:.3},\
         \"t_assemble_s\":{t_asm:.2},\"t_modes_s\":{t_modes:.2},\"t_rom_s\":{t_rom:.2},\"t_scan_s\":{t_scan:.2},\"peaks\":[{}]}}",
        args[1],
        fe.ndof(),
        positions.len(),
        complex.tets.len(),
        modes.lambda.len(),
        rom.dim(),
        a_lips * 1e6,
        peak_json.join(",")
    );
}
