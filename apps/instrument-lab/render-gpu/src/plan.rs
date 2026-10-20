//! An exported instrument (fs_vkplate::InstrumentImage) flattened into the arrays a GPU step reads: one
//! global node numbering over all bodies, work items of at most 256 nodes inside one body, contacts and
//! frictions grouped into node-disjoint batches in list order, microphone taps grouped by delay, and the
//! dense inverse of each plate's Airy operator. Pure data; no device code. A batch of simulations shares
//! everything here except the state and the strike lists.

use fs_vkplate::{GasState, Instrument, InstrumentImage, Law};
use std::collections::BTreeSet;
use std::sync::Arc;
use fs_vkplate::materials::KitMaterials;

pub const TPB: usize = 256;
/// Microphone taps per segment (two per lane of a warp).
const SEG: usize = 64;

#[derive(Default)]
pub struct Plan {
    pub nt: usize,
    pub nb: usize,
    pub nc: usize,
    pub rate: f64,
    pub k: f64,
    pub gas_density: f64,
    pub gas_bulk_modulus: f64,
    // Every batch member shares gas-dependent operators and microphone taps.
    // Keep the complete admitted state, not just the two kernel coefficients.
    gas: Option<GasState>,
    // Batches also share solid operators, contact laws and striker inertia.
    // Retain the full source binding, even if two cards resolve equal numbers.
    materials: Option<Arc<KitMaterials>>,
    pub body_off: Vec<usize>,
    // per node
    pub bias: Vec<f64>,
    pub ia: Vec<f64>,
    pub cw: Vec<f64>,
    pub cwp: Vec<f64>,
    pub cs1: Vec<f64>,
    pub ms0: Vec<f64>,
    pub mass: Vec<f64>,
    pub z0: Vec<f64>,
    pub lap_ptr: Vec<i32>,
    pub lap_col: Vec<i32>,
    pub lap_val: Vec<f64>,
    pub stf_ptr: Vec<i32>,
    pub stf_col: Vec<i32>,
    pub stf_val: Vec<f64>,
    pub gt_ptr: Vec<i32>,
    pub gt_core: Vec<i32>,
    pub gt_slot: Vec<i32>,
    // the same three operators as fixed-width column-major tables (entry q of row i at q * nt + i), rows
    // padded with zero weights on their own node (the CSR order kept, so each row sums as on the CPU)
    pub kl: usize,
    pub lap_ecol: Vec<i32>,
    pub lap_eval: Vec<f64>,
    pub ks: usize,
    pub stf_ecol: Vec<i32>,
    pub stf_eval: Vec<f64>,
    pub kg: usize,
    pub gt_ecore: Vec<i32>,
    pub gt_eslot: Vec<i32>,
    pub bstiff: Vec<i32>,
    // per body
    pub bh: Vec<f64>,
    pub bs1: Vec<f64>,
    pub bten: Vec<f64>,
    pub bber: Vec<f64>,
    pub bairm: Vec<f64>,
    pub barea: Vec<f64>,
    pub beh: Vec<f64>,
    pub bvk: Vec<i32>,
    pub bcoff: Vec<i32>,
    pub bnc: Vec<i32>,
    pub bitem0: Vec<i32>,
    pub bitem1: Vec<i32>,
    pub bcitem0: Vec<i32>,
    pub bcitem1: Vec<i32>,
    pub bcav: Vec<i32>,
    pub brad: Vec<i32>,
    // core nodes of the von Karman plates
    pub cnb: Vec<i32>,
    pub cl0: Vec<f64>,
    pub vk_body: Vec<i32>,
    pub airy_off: Vec<i64>,
    pub airy: Vec<f64>,
    // work items
    pub item_body: Vec<i32>,
    pub item_start: Vec<i32>,
    pub item_len: Vec<i32>,
    pub citem_body: Vec<i32>,
    pub citem_start: Vec<i32>,
    pub citem_len: Vec<i32>,
    // cavities, radiation damping, links, drives
    pub cav_hptr: Vec<i32>,
    pub cav_head: Vec<i32>,
    pub cav_vol: Vec<f64>,
    pub cav_port: Vec<i32>,
    pub port_area: Vec<f64>,
    pub port_len: Vec<f64>,
    pub port_res: Vec<f64>,
    pub rad_ptr: Vec<i32>,
    pub rad_body: Vec<i32>,
    pub rad_res: Vec<f64>,
    pub lk_a: Vec<i32>,
    pub lk_b: Vec<i32>,
    pub lk_k: Vec<f64>,
    pub lk_c: Vec<f64>,
    pub dv_node: Vec<i32>,
    pub dv_ptr: Vec<i32>,
    /// Per node, the links and drives acting on it (in list order).
    pub lnk_ptr: Vec<i32>,
    pub lnk_idx: Vec<i32>,
    pub drv_ptr: Vec<i32>,
    pub drv_idx: Vec<i32>,
    /// Nodes a contact, friction or strike touches (corrected by the contact block, not their own).
    pub cnode: Vec<i32>,
    pub cn_node: Vec<i32>,
    pub cn_body: Vec<i32>,
    pub dv_t: Vec<f64>,
    pub dv_v: Vec<f64>,
    // contacts and frictions
    pub ct_a0: Vec<i32>,
    pub ct_b0: Vec<i32>,
    pub ct_end: Vec<i32>,
    pub ct_node: Vec<i32>,
    pub ct_w: Vec<f64>,
    pub ct_gap: Vec<f64>,
    pub ct_kk: Vec<f64>,
    pub ct_alpha: Vec<f64>,
    pub ct_mu: Vec<f64>,
    pub ct_comp: Vec<f64>,
    pub ct_fx: Vec<i32>,
    pub cb_ptr: Vec<i32>,
    pub cb_idx: Vec<i32>,
    pub bc_ptr: Vec<i32>,
    pub bc_idx: Vec<i32>,
    pub fr_x0: Vec<i32>,
    pub fr_y0: Vec<i32>,
    pub fr_end: Vec<i32>,
    pub fr_node: Vec<i32>,
    pub fr_w: Vec<f64>,
    pub fr_contact: Vec<i32>,
    pub fr_mu: Vec<f64>,
    pub fr_cxx: Vec<f64>,
    pub fr_cyy: Vec<f64>,
    pub fr_cxy: Vec<f64>,
    pub fb_ptr: Vec<i32>,
    pub fb_idx: Vec<i32>,
    // microphones
    pub md_ptr: Vec<i32>,
    /// Each delay's taps in segments of at most SEG (boundaries into md_node), and each delay's segments.
    pub seg_q0: Vec<i32>,
    pub md_seg: Vec<i32>,
    pub md_node: Vec<i32>,
    pub md_w: Vec<f64>,
    pub md_mic: Vec<i32>,
    pub md_delay: Vec<i32>,
    pub mic_roff: Vec<i32>,
    pub mic_rmask: Vec<i32>,
    pub ringtot: usize,
    pub pt_ptr: Vec<i32>,
    pub pt_cav: Vec<i32>,
    pub pt_delay: Vec<i32>,
    pub pt_lo: Vec<f64>,
    pub pt_hi: Vec<f64>,
    pub mic_tau: Vec<f64>,
    pub nm: usize,
}

/// Per-simulation state and strike list (the part of a batch that differs).
#[derive(Default, Clone)]
pub struct SimState {
    pub w: Vec<f64>,
    pub wp: Vec<f64>,
    pub lwp: Vec<f64>,
    pub psi: Vec<f64>,
    pub rprev: Vec<f64>,
    pub budget: Vec<f64>,
    pub ct_active: Vec<i32>,
    pub ct_force: Vec<f64>,
    pub ct_dnow: Vec<f64>,
    pub ct_peak: Vec<f64>,
    pub port_x: Vec<f64>,
    pub port_u: Vec<f64>,
    pub port_du: Vec<f64>,
    pub mic_hp: Vec<f64>,
    /// (time, striker body, contact, speed, clear), in firing order
    pub strikes: Vec<(f64, i32, i32, f64, f64)>,
}

fn i32s(v: usize) -> i32 {
    i32::try_from(v).expect("index fits i32")
}

/// Node-disjoint batches of items in list order: an item joins the open batch unless it shares a node
/// with it. Items sharing a node keep their list order, so solving each batch in parallel gives what
/// solving them one at a time gives.
fn disjoint_batches(items: &[Vec<usize>]) -> (Vec<i32>, Vec<i32>) {
    let mut ptr = vec![0];
    let mut idx = vec![];
    let mut used: BTreeSet<usize> = BTreeSet::new();
    for (j, nodes) in items.iter().enumerate() {
        if nodes.iter().any(|n| used.contains(n)) {
            ptr.push(i32s(idx.len()));
            used.clear();
        }
        used.extend(nodes.iter().copied());
        idx.push(i32s(j));
    }
    if !idx.is_empty() {
        ptr.push(i32s(idx.len()));
    }
    (ptr, idx)
}

/// Node-disjoint groups by greedy colouring in list order: each item takes the lowest colour no earlier
/// item it shares a node with holds. Groups run in colour order, so every item still sees the latest state
/// of its nodes (a Gauss-Seidel sweep), but items sharing a node no longer keep their list order. For the
/// rim frictions of a closed hi-hat (each shares slope nodes with its neighbours) this is 6 groups where
/// list order needs 603.
fn coloured_batches(items: &[Vec<usize>]) -> (Vec<i32>, Vec<i32>) {
    let mut colour_of_node: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    let mut colours: Vec<usize> = vec![];
    for nodes in items {
        let taken: BTreeSet<usize> = nodes.iter().flat_map(|n| colour_of_node.get(n).into_iter().flatten().copied()).collect();
        let c = (0..).find(|c| !taken.contains(c)).unwrap();
        for &n in nodes {
            colour_of_node.entry(n).or_default().push(c);
        }
        colours.push(c);
    }
    let n = colours.iter().max().map_or(0, |&m| m + 1);
    let mut ptr = vec![0];
    let mut idx = vec![];
    for c in 0..n {
        idx.extend(colours.iter().enumerate().filter(|&(_, &k)| k == c).map(|(j, _)| i32s(j)));
        ptr.push(i32s(idx.len()));
    }
    (ptr, idx)
}

impl Plan {
    /// Flatten `img`; the dense Airy inverses come from `ins`'s own factors (the same instrument).
    /// `item`: nodes per work item (at most TPB; smaller spreads a small instrument over more blocks).
    /// `exact_friction`: frictions in the CPU's list order (one group per run of node-disjoint frictions:
    /// 603 groups on a closed hi-hat) instead of a coloured sweep (6 groups); for checks against the CPU.
    pub fn new(img: &InstrumentImage, ins: &Instrument, threads: usize, item: usize, exact_friction: bool) -> Result<Plan, String> {
        let item = item.clamp(32, TPB);
        if img.step_index != 0 {
            return Err("export the instrument before its first step (the microphone rings are not exported)".into());
        }
        let k = 1.0 / img.rate;
        if img.gas != *ins.gas() {
            return Err("exported and live instrument gas states differ".into());
        }
        if img.materials.as_deref() != ins.materials() {
            return Err("exported and live instrument solid material bindings differ".into());
        }
        let mut p = Plan {
            rate: img.rate, k, nb: img.bodies.len(),
            gas_density: img.gas.density,
            gas_bulk_modulus: img.gas.density * img.gas.sound_speed.powi(2),
            gas: Some(img.gas),
            materials: img.materials.clone(),
            ..Default::default()
        };
        let mut off = 0;
        for b in &img.bodies {
            p.body_off.push(off);
            off += b.n;
        }
        p.nt = off;
        let bo = p.body_off.clone();
        let g = |b: usize, i: usize| bo[b] + i;
        // nodes, operators, work items
        let (mut lap_ptr, mut stf_ptr) = (vec![0i32], vec![0i32]);
        let mut coff = 0usize;
        for (bi, b) in img.bodies.iter().enumerate() {
            let o = p.body_off[bi];
            if b.cf.len() != b.n {
                return Err(format!("{}: no update coefficients", b.name));
            }
            for i in 0..b.n {
                p.bias.push(b.bias[i]);
                let [ia, cw, cwp, cs1] = b.cf[i];
                p.ia.push(ia);
                p.cw.push(cw);
                p.cwp.push(cwp);
                p.cs1.push(cs1);
                p.ms0.push(b.mass[i] * b.s0[i]);
                p.mass.push(b.mass[i]);
                for (csr, ptr, col, val) in
                    [(&b.lap, &mut lap_ptr, &mut p.lap_col, &mut p.lap_val), (&b.stiff, &mut stf_ptr, &mut p.stf_col, &mut p.stf_val)]
                {
                    for q in csr.ptr[i]..csr.ptr[i + 1] {
                        col.push(i32s(o + csr.col[q]));
                        val.push(csr.val[q]);
                    }
                    ptr.push(i32s(col.len()));
                }
            }
            p.z0.extend(b.vk.as_ref().map_or(vec![0.0; b.n], |v| v.z0.clone()));
            p.bh.push(b.h);
            p.bs1.push(b.s1);
            p.bten.push(b.tension);
            p.bber.push(b.berger);
            p.bairm.push(b.air_mass);
            p.barea.push(b.area);
            p.bitem0.push(i32s(p.item_body.len()));
            for s0 in (0..b.n).step_by(item) {
                p.item_body.push(i32s(bi));
                p.item_start.push(i32s(o + s0));
                p.item_len.push(i32s(item.min(b.n - s0)));
            }
            p.bitem1.push(i32s(p.item_body.len()));
            p.bcoff.push(i32s(coff));
            p.bcitem0.push(i32s(p.citem_body.len()));
            if let Some(v) = &b.vk {
                let nc = v.nb.len();
                p.beh.push(v.eh);
                p.bvk.push(1);
                p.bnc.push(i32s(nc));
                p.vk_body.push(i32s(bi));
                for c in 0..nc {
                    for s in v.nb[c] {
                        p.cnb.push(s.map_or(-1, |i| i32s(o + i)));
                    }
                    p.cl0.push(v.l0[c]);
                }
                for s0 in (0..nc).step_by(item) {
                    p.citem_body.push(i32s(bi));
                    p.citem_start.push(i32s(coff + s0));
                    p.citem_len.push(i32s(item.min(nc - s0)));
                }
                coff += nc;
            } else {
                p.beh.push(0.0);
                p.bvk.push(0);
                p.bnc.push(0);
            }
            p.bcitem1.push(i32s(p.citem_body.len()));
        }
        p.nc = coff;
        p.lap_ptr = lap_ptr;
        p.stf_ptr = stf_ptr;
        // the nonlinear gradient as a gather: node k sums its core contributions in core order, the order
        // the CPU's scatter adds them in (each core node touches a node through at most one slot)
        let mut gather: Vec<Vec<(i32, i32)>> = vec![vec![]; p.nt];
        for (bi, b) in img.bodies.iter().enumerate() {
            if let Some(v) = &b.vk {
                let c0 = p.bcoff[bi] as usize;
                for (c, s9) in v.nb.iter().enumerate() {
                    for slot in [1usize, 2, 0, 3, 4, 5, 6, 7, 8] {
                        if let Some(node) = s9[slot] {
                            gather[g(bi, node)].push((i32s(c0 + c), slot as i32));
                        }
                    }
                }
            }
        }
        p.gt_ptr.push(0);
        for list in &gather {
            for &(c, s) in list {
                p.gt_core.push(c);
                p.gt_slot.push(s);
            }
            p.gt_ptr.push(i32s(p.gt_core.len()));
        }
        let nt = p.nt;
        let ell = |ptr: &[i32], col: &[i32], val: &[f64]| -> (usize, Vec<i32>, Vec<f64>) {
            let kmax = (0..nt).map(|i| (ptr[i + 1] - ptr[i]) as usize).max().unwrap_or(0).max(1);
            let (mut ec, mut ev) = (vec![0i32; kmax * nt], vec![0.0; kmax * nt]);
            for i in 0..nt {
                for q in 0..kmax {
                    let e = ptr[i] as usize + q;
                    let (c, v) = if e < ptr[i + 1] as usize { (col[e], val[e]) } else { (i32s(i), 0.0) };
                    ec[q * nt + i] = c;
                    ev[q * nt + i] = v;
                }
            }
            (kmax, ec, ev)
        };
        (p.kl, p.lap_ecol, p.lap_eval) = ell(&p.lap_ptr, &p.lap_col, &p.lap_val);
        (p.ks, p.stf_ecol, p.stf_eval) = ell(&p.stf_ptr, &p.stf_col, &p.stf_val);
        p.kg = gather.iter().map(|l| l.len()).max().unwrap_or(0).max(1);
        (p.gt_ecore, p.gt_eslot) = (vec![-1; p.kg * nt], vec![0; p.kg * nt]);
        for (i, list) in gather.iter().enumerate() {
            for (q, &(c, sl)) in list.iter().enumerate() {
                p.gt_ecore[q * nt + i] = c;
                p.gt_eslot[q * nt + i] = sl;
            }
        }
        p.bstiff = img.bodies.iter().map(|b| (!b.stiff.col.is_empty()) as i32).collect();
        // dense Airy inverses (row r of A^-1: column r by symmetry of A, taken as the solve of e_r)
        for (bi, b) in img.bodies.iter().enumerate() {
            if b.vk.is_none() {
                continue;
            }
            let nc = p.bnc[bi] as usize;
            p.airy_off.push(p.airy.len() as i64);
            let base = p.airy.len();
            p.airy.resize(base + nc * nc, 0.0);
            let body = &ins.bodies[bi];
            let rows: Vec<&mut [f64]> = p.airy[base..].chunks_mut(nc).collect();
            let per = nc.div_ceil(threads.max(1));
            std::thread::scope(|sc| {
                let mut rows = rows;
                let mut r0 = 0;
                while !rows.is_empty() {
                    let take = per.min(rows.len());
                    let chunk: Vec<&mut [f64]> = rows.drain(..take).collect();
                    let start = r0;
                    r0 += take;
                    sc.spawn(move || {
                        for (t, row) in chunk.into_iter().enumerate() {
                            row.fill(0.0);
                            row[start + t] = 1.0;
                            assert!(body.airy_solve(row));
                        }
                    });
                }
            });
        }
        // cavities, ports, radiation damping
        p.bcav = vec![-1; p.nb];
        p.brad = vec![-1; p.nb];
        p.cav_hptr.push(0);
        for (ci, c) in img.cavities.iter().enumerate() {
            for &h in &c.heads {
                if p.bcav[h] >= 0 {
                    return Err("a head of two cavities".into());
                }
                p.bcav[h] = i32s(ci);
                p.cav_head.push(i32s(h));
            }
            p.cav_hptr.push(i32s(p.cav_head.len()));
            p.cav_vol.push(c.volume);
            p.cav_port.push(c.port.is_some() as i32);
            let q = c.port.as_ref();
            p.port_area.push(q.map_or(0.0, |q| q.area));
            p.port_len.push(q.map_or(1.0, |q| q.length));
            p.port_res.push(q.map_or(0.0, |q| q.resistance));
        }
        p.rad_ptr.push(0);
        for (ri, r) in img.rads.iter().enumerate() {
            for &b in &r.bodies {
                if p.brad[b] >= 0 {
                    return Err("a body in two radiation groups".into());
                }
                p.brad[b] = i32s(ri);
                p.rad_body.push(i32s(b));
            }
            p.rad_ptr.push(i32s(p.rad_body.len()));
            p.rad_res.push(r.resistance);
        }
        for l in &img.links {
            p.lk_a.push(i32s(g(l.a.0, l.a.1)));
            p.lk_b.push(l.b.map_or(-1, |(b, i)| i32s(g(b, i))));
            p.lk_k.push(l.k);
            p.lk_c.push(l.c);
        }
        p.dv_ptr.push(0);
        for d in &img.drives {
            p.dv_node.push(i32s(g(d.body, d.node)));
            for &(t, v) in &d.knots {
                p.dv_t.push(t);
                p.dv_v.push(v);
            }
            p.dv_ptr.push(i32s(p.dv_t.len()));
        }
        let per_node = |pairs: Vec<(usize, i32)>| {
            let mut lists: Vec<Vec<i32>> = vec![vec![]; p.nt];
            for (node, idx) in pairs {
                lists[node].push(idx);
            }
            let mut ptr = vec![0i32];
            let mut idx = vec![];
            for l in lists {
                idx.extend(l);
                ptr.push(i32s(idx.len()));
            }
            (ptr, idx)
        };
        let mut lpairs = vec![];
        for (li, l) in img.links.iter().enumerate() {
            lpairs.push((g(l.a.0, l.a.1), i32s(li)));
            if let Some((b, i)) = l.b {
                lpairs.push((g(b, i), i32s(li)));
            }
        }
        (p.lnk_ptr, p.lnk_idx) = per_node(lpairs);
        (p.drv_ptr, p.drv_idx) = per_node(img.drives.iter().enumerate().map(|(di, d)| (g(d.body, d.node), i32s(di))).collect());
        // contacts: (node, weight) lists, compliance sums in the CPU's order, disjoint batches
        let cm = |b: usize, i: usize| img.bodies[b].cf[i][0];
        let mut ct_nodes = vec![];
        let mut body_contacts: Vec<Vec<i32>> = vec![vec![]; p.nb];
        for (j, c) in img.contacts.iter().enumerate() {
            p.ct_a0.push(i32s(p.ct_node.len()));
            let mut comp = 0.0;
            for &(b, i, w) in &c.a {
                p.ct_node.push(i32s(g(b, i)));
                p.ct_w.push(w);
                comp += w * w * cm(b, i);
            }
            p.ct_b0.push(i32s(p.ct_node.len()));
            for &(b, i, w) in &c.b {
                p.ct_node.push(i32s(g(b, i)));
                p.ct_w.push(w);
                comp += w * w * cm(b, i);
            }
            p.ct_end.push(i32s(p.ct_node.len()));
            let Law::Power { k: kk, alpha, mu } = c.law;
            p.ct_gap.push(c.gap);
            p.ct_kk.push(kk);
            p.ct_alpha.push(alpha);
            p.ct_mu.push(mu);
            p.ct_comp.push(comp);
            p.ct_fx.push(c.fixed_point_exit as i32);
            ct_nodes.push(c.a.iter().chain(&c.b).map(|&(b, i, _)| g(b, i)).collect::<Vec<_>>());
            let bodies: BTreeSet<usize> = c.a.iter().chain(&c.b).map(|&(b, _, _)| b).collect();
            for b in bodies {
                body_contacts[b].push(i32s(j));
            }
        }
        (p.cb_ptr, p.cb_idx) = disjoint_batches(&ct_nodes);
        p.bc_ptr.push(0);
        for l in &body_contacts {
            p.bc_idx.extend(l);
            p.bc_ptr.push(i32s(p.bc_idx.len()));
        }
        let mut fr_nodes = vec![];
        for fr in &img.frictions {
            p.fr_x0.push(i32s(p.fr_node.len()));
            let (mut cxx, mut cyy, mut cxy) = (0.0, 0.0, 0.0);
            for &(b, i, w) in &fr.sx {
                p.fr_node.push(i32s(g(b, i)));
                p.fr_w.push(w);
                let c = cm(b, i);
                cxx += w * w * c;
                for &(b2, i2, w2) in &fr.sy {
                    if b2 == b && i2 == i {
                        cxy += w * w2 * c;
                    }
                }
            }
            p.fr_y0.push(i32s(p.fr_node.len()));
            for &(b, i, w) in &fr.sy {
                p.fr_node.push(i32s(g(b, i)));
                p.fr_w.push(w);
                cyy += w * w * cm(b, i);
            }
            p.fr_end.push(i32s(p.fr_node.len()));
            p.fr_contact.push(i32s(fr.contact));
            p.fr_mu.push(fr.mu);
            p.fr_cxx.push(cxx);
            p.fr_cyy.push(cyy);
            p.fr_cxy.push(cxy);
            fr_nodes.push(fr.sx.iter().chain(&fr.sy).map(|&(b, i, _)| g(b, i)).collect::<Vec<_>>());
        }
        // frictions as a coloured sweep (each projection only removes energy given the state it sees, so
        // any Gauss-Seidel order is as valid as the CPU's list order; the hi-hat is chaotic either way)
        (p.fb_ptr, p.fb_idx) = if exact_friction { disjoint_batches(&fr_nodes) } else { coloured_batches(&fr_nodes) };
        p.cnode = vec![0; p.nt];
        for &n in ct_nodes.iter().chain(&fr_nodes).flatten() {
            p.cnode[n] = 1;
        }
        let node_body: Vec<usize> = (0..p.nb).flat_map(|b| std::iter::repeat_n(b, img.bodies[b].n)).collect();
        for (n, &c) in p.cnode.iter().enumerate() {
            if c != 0 {
                p.cn_node.push(i32s(n));
                p.cn_body.push(i32s(node_body[n]));
            }
        }
        if p.nb > 64 {
            return Err("more than 64 bodies".into());
        }
        // microphones: taps regrouped by the ring slot they land in (delay d: the low weight of taps at d,
        // the high weight of taps at d - 1)
        p.nm = img.mics.len();
        p.md_ptr.push(0);
        p.pt_ptr.push(0);
        let mut roff = 0;
        for (mi, m) in img.mics.iter().enumerate() {
            let dmax = m.taps.iter().map(|t| t.2 + 1).max().unwrap_or(0);
            let mut by: Vec<Vec<(i32, f64)>> = vec![vec![]; dmax + 1];
            for &(b, i, di, lo, hi) in &m.taps {
                by[di].push((i32s(g(b, i)), lo));
                by[di + 1].push((i32s(g(b, i)), hi));
            }
            for (d, list) in by.iter().enumerate() {
                if list.is_empty() {
                    continue;
                }
                if p.md_seg.is_empty() {
                    p.md_seg.push(0);
                    p.seg_q0.push(0);
                }
                for chunk in list.chunks(SEG) {
                    for &(node, w) in chunk {
                        p.md_node.push(node);
                        p.md_w.push(w);
                    }
                    p.seg_q0.push(i32s(p.md_node.len()));
                }
                p.md_seg.push(i32s(p.seg_q0.len() - 1));
                p.md_ptr.push(i32s(p.md_node.len()));
                p.md_mic.push(i32s(mi));
                p.md_delay.push(i32s(d));
            }
            if !m.ring_len.is_power_of_two() {
                return Err("microphone ring is not a power of two".into());
            }
            p.mic_roff.push(i32s(roff));
            p.mic_rmask.push(i32s(m.ring_len - 1));
            roff += m.ring_len;
            for &(ci, di, lo, hi) in &m.ptaps {
                p.pt_cav.push(i32s(ci));
                p.pt_delay.push(i32s(di));
                p.pt_lo.push(lo);
                p.pt_hi.push(hi);
            }
            p.pt_ptr.push(i32s(p.pt_cav.len()));
            p.mic_tau.push(m.unbaffled.map_or(-1.0, |a| 4.0 * a / (3.0 * std::f64::consts::PI * img.gas.sound_speed)));
        }
        p.ringtot = roff;
        if p.md_seg.is_empty() {
            p.md_seg.push(0);
            p.seg_q0.push(0);
        }
        Ok(p)
    }

    /// The state and strikes of one simulation of this plan's instrument.
    pub fn state(&self, img: &InstrumentImage) -> Result<SimState, String> {
        if self.gas != Some(img.gas) {
            return Err("a batch member's atmosphere differs from the first member's".into());
        }
        if self.materials.as_deref() != img.materials.as_deref() {
            return Err("a batch member's solid material bindings differ from the first member's".into());
        }
        if img.bodies.len() != self.nb || img.contacts.len() != self.ct_gap.len() || img.mics.len() != self.nm {
            return Err("a batch member is a different instrument (bodies, contacts or microphones differ)".into());
        }
        // the batch shares every operator: the same bodies, and every contact touching the same nodes with the
        // same weights, gap and law
        if img.bodies.iter().zip(&self.body_off).enumerate().any(|(b, (bi, &o))| o + bi.n != self.body_off.get(b + 1).copied().unwrap_or(self.nt)) {
            return Err("a batch member's bodies differ from the first member's".into());
        }
        let g = |b: usize, i: usize| self.body_off[b] + i;
        for (j, c) in img.contacts.iter().enumerate() {
            let nodes: Vec<i32> = c.a.iter().chain(&c.b).map(|&(b, i, _)| g(b, i) as i32).collect();
            let weights: Vec<u64> = c.a.iter().chain(&c.b).map(|&(_, _, w)| w.to_bits()).collect();
            let (a0, end) = (self.ct_a0[j] as usize, self.ct_end[j] as usize);
            let Law::Power { k: kk, alpha, mu } = c.law;
            if nodes != self.ct_node[a0..end]
                || weights.iter().ne(self.ct_w[a0..end].iter().map(|w| w.to_bits()).collect::<Vec<_>>().iter())
                || [c.gap, kk, alpha, mu].map(f64::to_bits) != [self.ct_gap[j], self.ct_kk[j], self.ct_alpha[j], self.ct_mu[j]].map(f64::to_bits)
            {
                return Err(format!("a batch member's contact {j} differs from the first member's"));
            }
        }
        let mut s = SimState::default();
        for b in &img.bodies {
            s.w.extend(&b.w);
            s.wp.extend(&b.wp);
            s.lwp.extend(&b.lwp);
            s.psi.push(b.psi);
            s.rprev.push(b.r_prev);
            s.budget.push(b.budget);
        }
        for c in &img.contacts {
            s.ct_active.push(c.active as i32);
            s.ct_force.push(c.force);
            s.ct_dnow.push(c.d_now);
            s.ct_peak.push(c.peak_force);
        }
        for c in &img.cavities {
            let q = c.port.as_ref();
            s.port_x.push(q.map_or(0.0, |q| q.x));
            s.port_u.push(q.map_or(0.0, |q| q.u));
            s.port_du.push(q.map_or(0.0, |q| q.du));
        }
        for m in &img.mics {
            s.mic_hp.extend([m.hp.0, m.hp.1]);
        }
        let mut last = f64::NEG_INFINITY;
        for st in img.strikes.iter().filter(|st| !st.done) {
            if st.time < last {
                return Err("strikes are not in time order".into());
            }
            last = st.time;
            s.strikes.push((st.time, i32s(st.striker), i32s(st.contact), st.speed, st.clear));
            if img.bodies[st.striker].n != 1 {
                return Err("a striker is not a point mass".into());
            }
        }
        Ok(s)
    }
}
