//! The program as flat device buffers for the GPU kernel (`kernel.h`): real constants, damper reals,
//! integer maps, damper flags and the state, with the offsets the kernel reads from its parameter block.
use crate::plan::{Mat, Program};

/// Parameter-block slots, in order; `kernel.h` sees them as `#define P_<NAME> <index>`.
pub const PARAMS: &[&str] = &[
    "N",
    "R",
    "NG",
    "NP",
    "NC",
    "STEPS",
    "DAMP",
    "MAXG",
    "MAXS",
    "GROUP_BLOCKS",
    "BOARD_BLOCKS",
    "RPW",
    "SKIP",
    // constants
    "QQ",
    "QV",
    "VQ",
    "VV",
    "BQ",
    "BV",
    "OM2",
    "A",
    "BETA",
    "HAM",
    "INV",
    "CROSS",
    "MASS",
    "STRAIN",
    "DUR",
    "GM",
    "R1",
    "R2",
    "R3",
    "R4",
    "R5",
    "R6",
    "K1",
    "K2",
    "K3",
    "K4",
    "K5",
    "K6",
    "K7",
    "E2",
    "E3",
    "GT",
    "CLIFT",
    // damper reals
    "ULOSS0",
    "ULOSS1",
    "SHAPE",
    "UC0",
    "UC1",
    "TM0",
    "TM1",
    "GTM0",
    "GTM1",
    "LAM0",
    "LAM1",
    "RHO0",
    "RHO1",
    "MU0",
    "MU1",
    "GUC0",
    "GUC1",
    "P0",
    "Q0",
    "V0",
    "AM0",
    "P1",
    "Q1",
    "V1",
    "AM1",
    // integers
    "GSTART",
    "GPAD",
    "GCON",
    "CONLIST",
    "CSTART",
    "CEND",
    "PSTART",
    "PEND",
    "PST",
    "PSH",
    "TOFF",
    "STPAD",
    "MPAD",
    "GBLK",
    // state
    "SQ",
    "SV",
    "SFQ",
    "SFV",
    "SVBT",
    "SSRC",
    "SDD",
    "SREA",
    "SREA2",
    "SSRC2",
    "SDD2",
    "SA0",
    "SA1",
];

pub fn slot(name: &str) -> usize {
    PARAMS
        .iter()
        .position(|p| *p == name)
        .expect("parameter slot")
}

/// `#define P_<NAME> <index>` for every slot.
pub fn defines() -> String {
    PARAMS
        .iter()
        .enumerate()
        .map(|(i, p)| format!("#define P_{p} {i}\n"))
        .collect()
}

#[derive(Default)]
struct Builder<T> {
    data: Vec<T>,
}
impl<T: Copy + Default> Builder<T> {
    fn push(&mut self, params: &mut [i32], name: &str, values: &[T]) {
        params[slot(name)] = self.data.len() as i32;
        self.data.extend_from_slice(values);
        // keep every array 4-aligned
        while self.data.len() % 4 != 0 {
            self.data.push(T::default());
        }
    }
}

pub struct Packed {
    pub params: Vec<i32>,
    pub consts: Vec<f64>,
    pub ints: Vec<i32>,
    pub damper: Vec<f64>,
    pub state_len: usize,
    pub group_blocks: usize,
    pub board_blocks: usize,
}

/// Warps per block (FS_PIANO_GPU_WARPS, default 8).
pub fn warps() -> usize {
    std::env::var("FS_PIANO_GPU_WARPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .filter(|w| [4, 8, 16, 32].contains(w))
        .unwrap_or(8)
}
pub fn threads() -> usize {
    32 * warps()
}

/// Group blocks and board blocks within `budget` co-resident blocks: one bridge group per block when
/// they fit, otherwise the lightest adjacent pairs share a block; board warps take one or two rows.
/// Returns the group block boundaries (in groups), the board block count and rows per warp.
fn assign(p: &Program, budget: usize) -> Result<(Vec<usize>, usize, usize), String> {
    let (ng, r, w) = (p.ng, p.r, warps());
    let mut bounds: Vec<usize> = (0..=ng).collect();
    if ng + r.div_ceil(w) <= budget {
        return Ok((bounds, r.div_ceil(w), 1));
    }
    let bb = r.div_ceil(2 * w);
    let gb = budget.saturating_sub(bb);
    if gb == 0 || ng > 2 * gb {
        return Err(format!(
            "{ng} bridge groups and {r} board rows do not fit {budget} co-resident blocks"
        ));
    }
    let weight = |g: usize| {
        let (s, e) = p.group_modes[g];
        let (pb, pe) = p.group_pads[g];
        (e - s)
            + p.pads[pb..pe]
                .iter()
                .map(|x| x.stations() * x.len())
                .sum::<usize>()
                / 4
    };
    let mut pairs: Vec<(usize, usize)> = (0..ng.saturating_sub(1))
        .map(|g| (weight(g) + weight(g + 1), g))
        .collect();
    pairs.sort();
    let mut merged = vec![false; ng];
    let mut need = ng - gb;
    for (_, g) in pairs {
        if need == 0 {
            break;
        }
        if merged[g] || merged[g + 1] {
            continue;
        }
        merged[g] = true;
        merged[g + 1] = true;
        need -= 1;
        bounds.retain(|&b| b != g + 1);
    }
    if need > 0 {
        return Err("could not pair enough bridge groups".into());
    }
    Ok((bounds, bb, 2))
}

pub fn pack(p: &Program, budget: usize) -> Result<Packed, String> {
    let (n, r, ng, np, nc) = (p.n, p.r, p.ng, p.np(), p.contacts.len());
    let (bounds, board_blocks, rpw) = assign(p, budget)?;
    let group_blocks = bounds.len() - 1;
    let mut params = vec![0i32; PARAMS.len()];
    params[slot("N")] = n as i32;
    params[slot("R")] = r as i32;
    params[slot("NG")] = ng as i32;
    params[slot("NP")] = np as i32;
    params[slot("NC")] = nc as i32;
    params[slot("GROUP_BLOCKS")] = group_blocks as i32;
    params[slot("BOARD_BLOCKS")] = board_blocks as i32;
    params[slot("RPW")] = rpw as i32;
    let block_pads = |b: usize| p.group_pads[bounds[b + 1] - 1].1 - p.group_pads[bounds[b]].0;
    let block_stations = |b: usize| {
        p.pads[p.group_pads[bounds[b]].0..p.group_pads[bounds[b + 1] - 1].1]
            .iter()
            .map(|x| x.stations())
            .sum::<usize>()
    };
    for b in 0..group_blocks {
        if block_pads(b) > warps() {
            return Err(format!(
                "a block has {} damper pads, more than its {} warps",
                block_pads(b),
                warps()
            ));
        }
        if block_stations(b) > threads() {
            return Err("a block has more damper stations than threads".into());
        }
    }
    if let Some(pad) = p.pads.iter().position(|x| x.stations() > 32) {
        return Err(format!("damper pad {pad} has more than 32 stations"));
    }
    let mut c = Builder::<f64>::default();
    for (name, v) in [
        ("QQ", &p.qq),
        ("QV", &p.qv),
        ("VQ", &p.vq),
        ("VV", &p.vv),
        ("BQ", &p.bq),
        ("BV", &p.bv),
        ("OM2", &p.omega2),
        ("A", &p.a),
        ("BETA", &p.beta),
        ("HAM", &p.hammer),
        ("INV", &p.inverse),
        ("CROSS", &p.cross),
        ("MASS", &p.mass),
        ("STRAIN", &p.strain),
    ] {
        c.push(&mut params, name, v);
    }
    c.push(&mut params, "DUR", &[p.duration]);
    let mats: [(&str, &Mat); 17] = [
        ("GM", &p.gm),
        ("R1", &p.r1),
        ("R2", &p.r2),
        ("R3", &p.r3),
        ("R4", &p.r4),
        ("R5", &p.r5),
        ("R6", &p.r6),
        ("K1", &p.k1),
        ("K2", &p.k2),
        ("K3", &p.k3),
        ("K4", &p.k4),
        ("K5", &p.k5),
        ("K6", &p.k6),
        ("K7", &p.k7),
        ("E2", &p.e2),
        ("E3", &p.e3),
        ("GT", &p.gt),
    ];
    for (name, m) in mats {
        c.push(&mut params, name, &m.data);
    }
    c.push(
        &mut params,
        "CLIFT",
        &p.contacts.iter().map(|x| x.lift).collect::<Vec<_>>(),
    );
    let mut i = Builder::<i32>::default();
    let gstart: Vec<i32> = p
        .group_modes
        .iter()
        .map(|g| g.0 as i32)
        .chain([n as i32])
        .collect();
    i.push(&mut params, "GSTART", &gstart);
    let gpad: Vec<i32> = p
        .group_pads
        .iter()
        .map(|g| g.0 as i32)
        .chain([np as i32])
        .collect();
    i.push(&mut params, "GPAD", &gpad);
    let mut gcon = vec![0i32];
    let mut conlist = Vec::new();
    for list in &p.group_contacts {
        conlist.extend(list.iter().map(|&x| x as i32));
        gcon.push(conlist.len() as i32);
    }
    i.push(&mut params, "GCON", &gcon);
    i.push(&mut params, "CONLIST", &conlist);
    i.push(
        &mut params,
        "CSTART",
        &p.contacts
            .iter()
            .map(|x| x.modes.0 as i32)
            .collect::<Vec<_>>(),
    );
    i.push(
        &mut params,
        "CEND",
        &p.contacts
            .iter()
            .map(|x| x.modes.1 as i32)
            .collect::<Vec<_>>(),
    );
    i.push(
        &mut params,
        "PSTART",
        &p.pads.iter().map(|x| x.modes.0 as i32).collect::<Vec<_>>(),
    );
    i.push(
        &mut params,
        "PEND",
        &p.pads.iter().map(|x| x.modes.1 as i32).collect::<Vec<_>>(),
    );
    let mut pst = vec![0i32];
    let mut psh = Vec::new();
    let mut shape_len = 0usize;
    for pad in &p.pads {
        pst.push(pst.last().unwrap() + pad.stations() as i32);
        psh.push(shape_len as i32);
        shape_len += pad.shape.len();
    }
    i.push(&mut params, "PST", &pst);
    i.push(&mut params, "PSH", &psh);
    let mut toff = Vec::new();
    let mut acc = 0usize;
    for pad in &p.pads {
        toff.push(acc as i32);
        acc += pad.stations() * pad.stations();
    }
    i.push(&mut params, "TOFF", &toff);
    let stpad: Vec<i32> = p
        .pads
        .iter()
        .enumerate()
        .flat_map(|(k, pad)| std::iter::repeat_n(k as i32, pad.stations()))
        .collect();
    i.push(&mut params, "STPAD", &stpad);
    let mut mpad = vec![-1i32; n];
    for (k, pad) in p.pads.iter().enumerate() {
        for x in &mut mpad[pad.modes.0..pad.modes.1] {
            *x = k as i32;
        }
    }
    i.push(&mut params, "MPAD", &mpad);
    i.push(
        &mut params,
        "GBLK",
        &bounds.iter().map(|&b| b as i32).collect::<Vec<_>>(),
    );
    let maxg = (0..group_blocks)
        .map(|b| p.group_modes[bounds[b + 1] - 1].1 - p.group_modes[bounds[b]].0)
        .max()
        .unwrap_or(0);
    let maxs = (0..group_blocks).map(block_stations).max().unwrap_or(0);
    params[slot("MAXG")] = maxg.max(1) as i32;
    params[slot("MAXS")] = maxs.max(1) as i32;
    // state
    let mut state_len = 0usize;
    for (name, len) in [
        ("SQ", n + r),
        ("SV", n + r),
        ("SFQ", n + r),
        ("SFV", n + r),
        ("SVBT", r),
        ("SSRC", ng),
        ("SDD", ng),
        ("SREA", ng),
        ("SREA2", ng),
        ("SSRC2", ng),
        ("SDD2", ng),
        ("SA0", np),
        ("SA1", np),
    ] {
        params[slot(name)] = state_len as i32;
        state_len += len.div_ceil(4) * 4;
    }
    let mut packed = Packed {
        params,
        consts: c.data,
        ints: i.data,
        damper: Vec::new(),
        state_len,
        group_blocks,
        board_blocks,
    };
    pack_dampers(p, &mut packed);
    Ok(packed)
}

/// The damper reals (they change with the pedal and the lifted courses); also sets the DAMP flag.
pub fn pack_dampers(p: &Program, packed: &mut Packed) {
    let np = p.np();
    let params = &mut packed.params;
    let mut d = Builder::<f64>::default();
    d.push(
        params,
        "ULOSS0",
        &p.pads.iter().map(|x| x.unit_loss[0]).collect::<Vec<_>>(),
    );
    d.push(
        params,
        "ULOSS1",
        &p.pads.iter().map(|x| x.unit_loss[1]).collect::<Vec<_>>(),
    );
    let cat =
        |f: &dyn Fn(&crate::plan::Pad) -> Vec<f64>| p.pads.iter().flat_map(f).collect::<Vec<f64>>();
    d.push(params, "SHAPE", &cat(&|x| x.shape.clone()));
    d.push(params, "UC0", &cat(&|x| x.uc[0].clone()));
    d.push(params, "UC1", &cat(&|x| x.uc[1].clone()));
    // the chains as matrices (plan::ChainOp): T and Lambda per pad (stations^2), rho and mu per station
    for dir in 0..2 {
        let name = |s: &str| format!("{s}{dir}");
        d.push(params, &name("TM"), &cat(&|x| x.op[dir].t.clone()));
        d.push(params, &name("GTM"), &cat(&|x| x.op[dir].gt.clone()));
        d.push(params, &name("GUC"), &cat(&|x| x.op[dir].guc.clone()));
        d.push(params, &name("LAM"), &cat(&|x| x.op[dir].lambda.clone()));
        d.push(params, &name("RHO"), &cat(&|x| x.op[dir].rho.clone()));
        d.push(params, &name("MU"), &cat(&|x| x.op[dir].mu.clone()));
    }
    // the same sizes whether or not any pad is engaged, so a reload never reallocates
    let r = p.r;
    for (dir, (pn, qn, vn, an)) in [("P0", "Q0", "V0", "AM0"), ("P1", "Q1", "V1", "AM1")]
        .into_iter()
        .enumerate()
    {
        let dm = &p.dampers;
        let get = |m: &[Mat; 2], len: usize| {
            if dm.any {
                m[dir].data.clone()
            } else {
                vec![0.0; len]
            }
        };
        d.push(params, pn, &get(&dm.p, np * r));
        d.push(params, qn, &get(&dm.q, np * np));
        d.push(params, vn, &get(&dm.v, r * r));
        d.push(params, an, &get(&dm.a, r * np));
    }
    params[slot("DAMP")] = p.dampers.any as i32;
    packed.damper = d.data;
}

/// Reals of threadgroup memory per block: the group's velocities, its stations' u and c, 4 NW scratch.
pub fn threadgroup_reals(p: &Packed) -> usize {
    let at = |s: &str| p.params[slot(s)] as usize;
    (at("MAXG") + at("MAXS") + threads() + 8 * warps()).div_ceil(2) * 2
}

/// The pads' active flags (an inactive pad's chain and superposition are skipped).
pub fn active_flags(p: &Program) -> Vec<i32> {
    let mut v: Vec<i32> = (0..p.np())
        .map(|i| p.dampers.active.get(i).copied().unwrap_or(false) as i32)
        .collect();
    while v.len() % 4 != 0 || v.is_empty() {
        v.push(0);
    }
    v
}

/// f64 -> double-float (hi, lo) pairs for a GPU without FP64.
pub fn to_df64(x: &[f64]) -> Vec<[f32; 2]> {
    x.iter()
        .map(|&v| {
            let hi = v as f32;
            [hi, (v - f64::from(hi)) as f32]
        })
        .collect()
}
pub fn from_df64(x: &[[f32; 2]], out: &mut [f64]) {
    for (o, v) in out.iter_mut().zip(x) {
        *o = f64::from(v[0]) + f64::from(v[1]);
    }
}
