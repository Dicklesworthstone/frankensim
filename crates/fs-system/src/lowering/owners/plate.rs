//! Modal circular plate owner (`PartLoweringKey(Shell, {Mechanics}, Modal)`, `.11.9.2`..`.11.9.4`).
//! Lowers `Shell` + `TensionedDisk` + `FidelityRung::Modal` onto `ModalStructure` via
//! `fs_modal::slice_window`. Damping is film viscoelastic bending loss only:
//! `zeta_k = 0.5 * eta(f_k) * (phi_k^T K_bend phi_k) / lambda_k` on the zero-prestress DKT
//! bending pencil (no air loading, radiation, or rim loss).

use core::f64::consts::{PI, TAU};
use core::num::NonZeroUsize;

use fs_blake3::ContentHash as Hash;
use fs_matdb::{MaterialCard, QueryPoint, SelectionPolicy};
use fs_modal::{ModePair, SliceOptions, slice_window};
use fs_nlmodal::component::{ModalAcousticTimeBudget, ModalAttachment, ModalStructure};
use fs_plate::shell::head::{TensionedDisk, TensionedDiskSpec};
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::{AssemblyOptions, EdgeSupport, assemble};
use fs_ports::{PortKind, StableId};
use fs_qty::semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QtyAny, QuantitySpec};
use fs_scenario::assembly::draft::{
    DraftGeometrySource, DraftInterfaceDecl, DraftInterfaceKind, DraftPartKind, FidelityRung,
    PhysicsDomain, SourceTarget,
};
use fs_scenario::entity::EntityId;

use super::{Scoped, owning_part};
use crate::RouteTemperature;
use crate::lowering::{
    AREA_DIMS, FREQUENCY_DIMS, FidelityRungKind, LENGTH_DIMS, LoweredPartOutput, LoweringError,
    LoweringGraph, LoweringRegistry, PRESSURE_DIMS, PartLoweringContext, PartLoweringKey,
    TEMPERATURE_DIMS, ambient_heat_sink,
};

/// Relative eigenvalue gap below which a cut between `lambda_N` and `lambda_{N+1}` refuses.
pub const DEGENERATE_GAP_TOL: f64 = 1.0e-6;
/// Maximum Weyl-window doublings while certifying at least `min(N + 1, free)` modes.
pub const MAX_WINDOW_DOUBLINGS: usize = 24;
const J01: f64 = 2.404_825_557_695_773;
const LOSS_PROP: &str = "loss_tangent_tensile";
const DENSITY_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const TENSION_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const DISK_BUDGET: ProfileBudget = ProfileBudget {
    max_nodes: 10_000,
    max_triangles: 20_000,
    max_feature_evaluations: 100_000,
};

type Disk<'a> = (TensionedDisk, Vec<ModePair>, &'a MaterialCard, Hash, f64);
type Strike = (EntityId, [f64; 2]);
type Ports<'a> = (&'a [EntityId], &'a [Strike], bool);

/// Register `PartLoweringKey(Shell, {Mechanics}, Modal)`.
pub(crate) fn register(registry: &mut LoweringRegistry) {
    let k = PartLoweringKey::new(
        DraftPartKind::Shell,
        [PhysicsDomain::Mechanics],
        FidelityRungKind::Modal,
    );
    registry.register_part(k, lower_modal_plate);
}

fn is_part_face_power_port(
    ctx: &PartLoweringContext<'_>,
    idx: usize,
    iface: &DraftInterfaceDecl,
) -> bool {
    let DraftInterfaceKind::PowerPort { schema } = &iface.kind else {
        return false;
    };
    let inp = ctx.draft.input();
    (owning_part(inp, iface.feature_a) == Some(idx)
        || owning_part(inp, iface.feature_b) == Some(idx))
        && ctx
            .cards
            .power_port_schema(schema)
            // Matches fs-runtime's RUNTIME_INPUT_FAMILY live-input convention.
            .is_some_and(|c| c.family == "runtime-input" && c.description == "face")
}

/// Distinct surface features of `part` touched by a `PistonFace` interface or a
/// surface-targeted `mechanics` source, in draft order, plus `HertzContact` strikes
/// and whether `(part, "face")` is targeted.
fn collect_ports(
    ctx: &PartLoweringContext<'_>,
    it: &Scoped,
) -> Result<(Vec<EntityId>, Vec<Strike>, bool), LoweringError> {
    let (inp, part) = (ctx.draft.input(), ctx.part.entity);
    let Some(idx) = inp.parts.iter().position(|p| p.entity == part) else {
        return Ok((Vec::new(), Vec::new(), false));
    };
    let owns = |e| owning_part(inp, e) == Some(idx);
    let (mut out, mut strikes, mut part_face) = (Vec::new(), Vec::new(), false);
    let mut push = |e: EntityId| {
        if !out.contains(&e) {
            out.push(e);
        }
    };
    for iface in &inp.interfaces {
        if !matches!(iface.kind, DraftInterfaceKind::PistonFace { .. }) {
            continue;
        }
        for feat in [iface.feature_a, iface.feature_b] {
            if owns(feat) {
                if feat == part {
                    return Err(it.invalid(
                        "piston-face on modal shell must target a surface feature, not the part entity",
                    ));
                }
                push(feat);
            }
        }
    }
    for src in &inp.sources {
        if let SourceTarget::Port { entity, port } = &src.target {
            if *entity == part && port == "face" {
                part_face = true;
            } else if port == "mechanics" && *entity != part && owns(*entity) {
                push(*entity);
            }
        }
    }
    part_face |= inp
        .interfaces
        .iter()
        .any(|i| is_part_face_power_port(ctx, idx, i));
    for iface in &inp.interfaces {
        let DraftInterfaceKind::HertzContact { state } = &iface.kind else {
            continue;
        };
        let xy = |n| {
            state
                .iter()
                .find(|(k, _)| k.name == n)
                .map_or(0.0, |(_, q)| q.value)
        };
        for feat in [iface.feature_a, iface.feature_b] {
            if !owns(feat) {
                continue;
            }
            if feat == part {
                return Err(
                    it.invalid("hertz-contact on modal shell must target a surface feature")
                );
            }
            if out.contains(&feat) {
                return Err(
                    it.invalid("surface feature cannot combine piston-face and hertz-contact")
                );
            }
            if strikes.iter().any(|(f, _)| *f == feat) {
                return Err(it.invalid(
                    "surface feature cannot be shared by two hertz-contact interfaces; use a second feature",
                ));
            }
            strikes.push((feat, [xy("station_x"), xy("station_y")]));
        }
    }
    Ok((out, strikes, part_face))
}

fn has_center_mech_connection(
    ctx: &PartLoweringContext<'_>,
    faces: &[EntityId],
    strikes: &[Strike],
) -> bool {
    let (inp, part) = (ctx.draft.input(), ctx.part.entity);
    let Some(idx) = inp.parts.iter().position(|p| p.entity == part) else {
        return false;
    };
    let owns_non_face = |e| {
        owning_part(inp, e) == Some(idx)
            && !faces.contains(&e)
            && !strikes.iter().any(|(f, _)| *f == e)
    };
    inp.sources.iter().any(|s| {
        matches!(&s.target, SourceTarget::Port { entity, port } if *entity == part && port == "mechanics")
    }) || inp.joints.iter().any(|j| j.parent == part || j.child == part)
        || inp.interfaces.iter().any(|i| {
            !is_part_face_power_port(ctx, idx, i)
                && (owns_non_face(i.feature_a) || owns_non_face(i.feature_b))
        })
}

fn select_non_splitting_modes(
    film: &TensionedDisk,
    n: usize,
    b0: f64,
    it: &Scoped,
) -> Result<Vec<ModePair>, LoweringError> {
    let free = film.model.free;
    if n == 0 || n > free {
        return Err(it.invalid(format!("mode_count {n} out of range 1..={free}")));
    }
    let certify = |need: usize| -> Result<Vec<ModePair>, LoweringError> {
        let (mut b, opts) = (b0, SliceOptions::default());
        for _ in 0..=MAX_WINDOW_DOUBLINGS {
            if !(b.is_finite() && b > 0.0) {
                break;
            }
            let rep = slice_window(&film.model.k, &film.model.m, (0.0, b), &opts)
                .map_err(|e| it.invalid(format!("modal slice_window refused: {e}")))?;
            if rep.modes.len() >= need {
                return Ok(rep.modes);
            }
            b *= 2.0;
        }
        Err(it.invalid(format!(
            "modal window doubling limit ({MAX_WINDOW_DOUBLINGS}) exhausted before {need} modes"
        )))
    };
    let rel = |ms: &[ModePair], k: usize| {
        let (p, c) = (ms[k - 1].lambda, ms[k].lambda);
        (p.is_finite() && p > 0.0 && c.is_finite() && c >= p)
            .then_some((c - p) / p)
            .ok_or_else(|| it.invalid("non-positive or non-monotone modal eigenvalues"))
    };
    let mut modes = certify((n + 1).min(free))?;
    if n < free && rel(&modes, n)? <= DEGENERATE_GAP_TOL {
        let (gap, mut lo, mut hi) = (rel(&modes, n)?, n - 1, n + 1);
        while lo > 0 && rel(&modes, lo)? <= DEGENERATE_GAP_TOL {
            lo -= 1;
        }
        while hi < free {
            if modes.len() <= hi {
                modes = certify((hi + 1).min(free))?;
            }
            if rel(&modes, hi)? > DEGENERATE_GAP_TOL {
                break;
            }
            hi += 1;
        }
        let hint = match lo {
            0 => format!("{hi}"),
            _ => format!("{lo} or {hi}"),
        };
        return Err(it.invalid(format!(
            "mode_count {n} splits a degenerate pair (gap {gap:.3e} <= {DEGENERATE_GAP_TOL:.0e}); choose {hint}"
        )));
    }
    modes.truncate(n);
    Ok(modes)
}

fn bend_shares(f: &TensionedDisk, ms: &[ModePair], it: &Scoped) -> Result<Vec<f64>, LoweringError> {
    let opts = AssemblyOptions {
        pretension: 0.0,
        support: EdgeSupport::SimplySupported,
    };
    let bend = assemble(&f.mesh, &f.section, f.outer_nodes(), &[], &opts)
        .map_err(|e| it.invalid(format!("zero-prestress bending assembly refused: {e}")))?;
    if bend.free != f.model.free || bend.dof_map != f.model.dof_map {
        return Err(it.invalid("bending-only assembly changed reduced coordinates"));
    }
    let (mut prod, mut out) = (vec![0.0; bend.free], Vec::with_capacity(ms.len()));
    for m in ms {
        if m.phi.len() != bend.free || !m.lambda.is_finite() || m.lambda <= 0.0 {
            return Err(it.invalid("non-positive or mismatched modal eigenpair"));
        }
        bend.k.spmv(&m.phi, &mut prod);
        let sh = m.phi.iter().zip(&prod).map(|(a, b)| a * b).sum::<f64>() / m.lambda;
        if !(sh.is_finite() && (0.0..=1.0).contains(&sh)) {
            return Err(it.invalid(format!("bending share {sh} outside [0, 1]")));
        }
        out.push(sh);
    }
    Ok(out)
}

fn loss_tan(
    (card, hash, pol): (&MaterialCard, Hash, SelectionPolicy),
    (temp_k, freq_hz): (f64, f64),
    it: &mut Scoped,
) -> Result<f64, LoweringError> {
    let qs = |k| QuantitySpec::semantic(SemanticType::new(k, ValueForm::Static));
    let pt = QueryPoint::new()
        .with_quantity("T", qs(QuantityKind::AbsoluteTemperature), temp_k)
        .and_then(|qp| {
            qp.with_quantity(
                "frequency",
                qs(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
                freq_hz,
            )
        })
        .map_err(|e| it.invalid(format!("loss query point refused: {e}")))?;
    let ans = card
        .claims()
        .query(LOSS_PROP, &pt, pol)
        .map_err(|s| LoweringError::MatDb {
            item: it.item.clone(),
            source: Box::new(s),
        })?;
    let v = ans.evidence.value.value;
    if !(ans.evidence.value.dims == Dims::NONE && v.is_finite() && v > 0.0) {
        return Err(it.invalid(format!("{LOSS_PROP} at {freq_hz:.3} Hz not positive [1]")));
    }
    it.r.record_matdb_claim(LOSS_PROP, &ans.evidence.value, hash, &ans.receipt);
    Ok(v)
}

fn load_elastic<'a>(
    ctx: &'a PartLoweringContext<'_>,
    it: &mut Scoped,
) -> Result<(&'a MaterialCard, Hash, [f64; 4]), LoweringError> {
    let [mat] = ctx.part.materials.as_slice() else {
        return Err(it.invalid("modal shell requires exactly one material card"));
    };
    if mat.require_known(LOSS_PROP).is_err() {
        return Err(LoweringError::UnknownPropertyDefaulted {
            item: it.item.clone(),
            property: LOSS_PROP.into(),
        });
    }
    let env = ctx.query_env();
    let mut chk = |p, dm, lo, hi| {
        let (s, _) =
            it.r.query_material_property(&it.item, ctx.cards, mat, p, p, env)?;
        (s.dims == dm && s.value.is_finite() && s.value > lo && s.value < hi)
            .then_some(s.value)
            .ok_or_else(|| it.invalid(format!("{p} out of physical domain")))
    };
    let rho = chk("density", DENSITY_DIMS, 0.0, f64::INFINITY)?;
    let e = chk("young_modulus", PRESSURE_DIMS, 0.0, f64::INFINITY)?;
    let nu = chk("poisson_ratio", Dims::NONE, -1.0, 0.5)?;
    let card = ctx
        .cards
        .material(&mat.card)
        .ok_or_else(|| it.invalid("missing card"))?;
    let amb = ctx.draft.input().environment.base.ambient_temperature;
    let t_q = mat
        .query_point
        .iter()
        .find(|(k, _)| k.name == "T" || k.name == "temperature")
        .map_or(amb, |(_, q)| *q);
    if !(t_q.dims == TEMPERATURE_DIMS && t_q.value.is_finite() && t_q.value > 0.0) {
        return Err(it.invalid("material query temperature must be positive finite [K]"));
    }
    Ok((card, mat.card, [rho, e, nu, t_q.value]))
}

fn build_disk<'a>(
    ctx: &'a PartLoweringContext<'_>,
    it: &mut Scoped,
) -> Result<Disk<'a>, LoweringError> {
    let FidelityRung::Modal { mode_count } = *ctx.chosen_rung else {
        return Err(it.invalid("expected Modal fidelity rung"));
    };
    if let Some(i) = &ctx.part.initial {
        let z = |s: &[QtyAny]| s.iter().all(|q| q.value == 0.0);
        let id_q = (i.orientation[0].abs() - 1.0).abs() <= 1.0e-12
            && i.orientation[1..].iter().all(|&c| c == 0.0);
        if !(z(&i.position) && z(&i.linear_velocity) && z(&i.angular_velocity) && id_q) {
            return Err(it.invalid("modal shell does not admit non-zero rigid-body initial state"));
        }
    }
    let DraftGeometrySource::TensionedDisk(d) = &ctx.part.geometry else {
        return Err(it.invalid("modal shell requires TensionedDisk geometry"));
    };
    let mut drv = |k, v, dm| it.derive(ctx.part, k, v, dm, k);
    drv("radius", d.radius.value, LENGTH_DIMS);
    drv("thickness", d.thickness.value, LENGTH_DIMS);
    drv("tension", d.tension.value, TENSION_DIMS);
    drv(
        "radial_intervals",
        f64::from(d.radial_intervals),
        Dims::NONE,
    );
    drv("azimuths", f64::from(d.azimuths), Dims::NONE);
    let (card, hash, [rho, e, nu, temp_k]) = load_elastic(ctx, it)?;
    let spec = TensionedDiskSpec {
        radius_m: d.radius.value,
        thickness_m: d.thickness.value,
        young_pa: e,
        poisson: nu,
        density_kg_m3: rho,
        tension_n_m: d.tension.value,
        radial_intervals: d.radial_intervals as usize,
        azimuths: d.azimuths as usize,
    };
    let film = TensionedDisk::new(spec, DISK_BUDGET)
        .map_err(|er| it.invalid(format!("TensionedDisk refused: {er}")))?;
    let (n, k0) = (mode_count as usize, J01 / spec.radius_m);
    let d_rig = e * spec.thickness_m.powi(3) / (12.0 * (1.0 - nu * nu));
    let b0 = ((2 * (n + 1)) as f64).max(4.0) * (spec.tension_n_m * k0 * k0 + d_rig * k0.powi(4))
        / (rho * spec.thickness_m);
    let modes = select_non_splitting_modes(&film, n, b0, it)?;
    Ok((film, modes, card, hash, temp_k))
}

fn build_attachments(
    film: &TensionedDisk,
    modes: &[ModePair],
    (faces, strikes, part_face): Ports<'_>,
    ctx: &PartLoweringContext<'_>,
    it: &mut Scoped,
) -> Result<Vec<ModalAttachment>, LoweringError> {
    let ctr = film
        .mesh
        .nodes
        .iter()
        .position(|&(x, y)| x == 0.0 && y == 0.0)
        .and_then(|i| film.model.dof_map[3 * i])
        .ok_or_else(|| it.invalid("no free center DOF"))?;
    let disk_area = PI * film.spec.radius_m * film.spec.radius_m;
    if part_face || !faces.is_empty() {
        it.derive(ctx.part, "disk_area", disk_area, AREA_DIMS, "pi*r^2");
    }
    let mk = |name: &str, weights| -> Result<ModalAttachment, LoweringError> {
        let id = StableId::new(name).map_err(|e| it.invalid(format!("port id: {e:?}")))?;
        Ok(ModalAttachment {
            id,
            kind: PortKind::MechanicalForceVelocity,
            width: NonZeroUsize::MIN,
            weights,
        })
    };
    let mut atts = Vec::with_capacity(1 + usize::from(part_face) + faces.len() + strikes.len());
    atts.push(mk("mechanics", modes.iter().map(|m| m.phi[ctr]).collect())?);
    if part_face || !faces.is_empty() {
        let face_w: Vec<f64> = modes
            .iter()
            .map(|m| {
                let a = film
                    .modal_area(&m.phi)
                    .map_err(|e| it.invalid(format!("modal_area refused: {e}")))?;
                a.is_finite()
                    .then_some(a / disk_area)
                    .ok_or_else(|| it.invalid("non-finite modal_area"))
            })
            .collect::<Result<_, _>>()?;
        if part_face {
            atts.push(mk("face", face_w.clone())?);
        }
        for idx in 0..faces.len() {
            atts.push(mk(&format!("face_{idx}"), face_w.clone())?);
        }
    }
    for (idx, &(_, [x, y])) in strikes.iter().enumerate() {
        let sw: Vec<f64> = modes
            .iter()
            .map(|m| {
                film.sample_mode_at(&m.phi, x, y)
                    .map_err(|e| it.invalid(format!("strike station ({x}, {y}) refused: {e}")))
            })
            .collect::<Result<_, _>>()?;
        if sw.iter().all(|w| w.abs() <= 1.0e-12) {
            return Err(it.invalid("strike station has zero modal participation (fixed rim)"));
        }
        atts.push(mk(&format!("strike_{idx}"), sw)?);
    }
    Ok(atts)
}

fn lower_modal_plate(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let mut it = Scoped::new("part", &ctx.part.name);
    let (film, modes, card, hash, temp_k) = build_disk(ctx, &mut it)?;
    let shares = bend_shares(&film, &modes, &it)?;
    let q = ctx.draft.input().time_plan.default_sample_rate;
    it.declare_on("time_plan", "default_sample_rate", q);
    let (v, r) = (q.value, q.value.round());
    let ok_r = q.dims == FREQUENCY_DIMS && v.is_finite() && v > 0.0 && r <= f64::from(u32::MAX);
    if !ok_r || (v - r).abs() > 1.0e-9 * r.max(1.0) {
        return Err(it.invalid(format!("default_sample_rate {v} Hz invalid")));
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let rate = r as u32;
    let (pol, mut pairs) = (
        (card, hash, ctx.selection_policy),
        Vec::with_capacity(modes.len()),
    );
    let s_note =
        "phi^T K_bend phi / lambda (film material bending loss only; no air/radiation/rim loss)";
    let z_note = "0.5 * loss_tangent_tensile(f_k) * bending_share[k]";
    for (k, (mode, share)) in modes.iter().zip(shares).enumerate() {
        let omega = mode.lambda.sqrt();
        let eta = loss_tan(pol, (temp_k, omega / TAU), &mut it)?;
        let zeta = 0.5 * eta * share;
        let mut rec = |n, v, f| it.r.record_declared_field(n, v, Dims::NONE, &it.item, f);
        rec(format!("bending_share[{k}]"), share, s_note);
        rec(format!("modal_damping_ratio[{k}]"), zeta, z_note);
        pairs.push((omega, zeta));
    }
    let (faces, strikes, part_face) = collect_ports(ctx, &it)?;
    let atts = build_attachments(&film, &modes, (&faces, &strikes, part_face), ctx, &mut it)?;
    let bdg = ModalAcousticTimeBudget::audible_reference();
    let st = ModalStructure::from_modes(rate, &pairs, bdg, atts)
        .map_err(|e| it.invalid(format!("ModalStructure refused: {e:?}")))?;
    let rest = st.at_rest();
    let node = graph.add_with_tangent(&ctx.part.name, Some(ctx.part.entity), st, rest)?;
    let (mech, heat) = (graph.port(node, "mechanics")?, graph.port(node, "heat")?);
    graph.expose_port(ctx.part.entity, "mechanics", mech);
    if part_face {
        graph.expose_port(ctx.part.entity, "face", graph.port(node, "face")?);
    }
    for (idx, &feat) in faces.iter().enumerate() {
        let fp = graph.port(node, &format!("face_{idx}"))?;
        graph.expose_port(feat, "mechanics", fp);
    }
    for (idx, &(feat, _)) in strikes.iter().enumerate() {
        let sp = graph.port(node, &format!("strike_{idx}"))?;
        graph.expose_port(feat, "mechanics", sp);
    }
    if !has_center_mech_connection(ctx, &faces, &strikes) {
        graph.bind_zero_boundary(mech);
    }
    let t_amb = ctx.ambient_temperature_k();
    let sink = ambient_heat_sink(&it.item, t_amb)?;
    graph.add_heat_route(heat, vec![], Some(sink), RouteTemperature::Fixed(t_amb));
    let surf = radiating_surface(ctx, &film, &modes, node, &it)?;
    graph.register_radiating_surface(ctx.part.entity, surf);
    Ok(it.part(&ctx.part.name))
}

fn radiating_surface(
    ctx: &PartLoweringContext<'_>,
    film: &TensionedDisk,
    modes: &[ModePair],
    node: crate::NodeId,
    it: &Scoped,
) -> Result<crate::acoustic::RadiatingModalSurface, LoweringError> {
    let inp = ctx.draft.input();
    let idx = inp.parts.iter().position(|p| p.entity == ctx.part.entity);
    let mut sign: Option<f64> = None;
    for i in &inp.interfaces {
        let DraftInterfaceKind::PistonFace { state } = &i.kind else {
            continue;
        };
        if idx.is_none_or(|k| {
            owning_part(inp, i.feature_a) != Some(k) && owning_part(inp, i.feature_b) != Some(k)
        }) {
            continue;
        }
        let s = if state
            .iter()
            .find(|(k, _)| k.name == "sign")
            .map_or(1.0, |(_, q)| q.value)
            < 0.0
        {
            -1.0
        } else {
            1.0
        };
        if sign.is_some_and(|p| (p - s).abs() > 0.0) {
            return Err(it.invalid("conflicting piston-face signs on modal shell"));
        }
        sign = Some(s);
    }
    let outward_sign = sign.unwrap_or(1.0);
    let mut nodal_area = vec![0.0; film.mesh.nodes.len()];
    for tri in &film.mesh.tris {
        let [a, b, c] = tri.map(|i| film.mesh.nodes[i]);
        let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1));
        if !(area.is_finite() && area > 0.0) {
            return Err(it.invalid("disk mesh has non-positive facet"));
        }
        for &ni in tri {
            nodal_area[ni] += area / 3.0;
        }
    }
    let samples = film
        .mesh
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(ni, &(x, y))| {
            let dof = film.model.dof_map[3 * ni]?;
            Some(crate::acoustic::RadiatingPatch {
                xy_m: [x, y],
                area_m2: nodal_area[ni],
                shapes: modes.iter().map(|m| outward_sign * m.phi[dof]).collect(),
            })
        })
        .collect();
    let mode_omegas_rad_s = modes.iter().map(|m| m.lambda.sqrt()).collect();
    Ok(crate::acoustic::RadiatingModalSurface {
        node,
        mode_count: modes.len(),
        mode_omegas_rad_s,
        outward_sign,
        samples,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn doublet_cut_and_bending_share_bounds() {
        let s = TensionedDiskSpec {
            radius_m: 0.1778,
            thickness_m: 0.00025,
            young_pa: 4.903_325e9,
            poisson: 0.38,
            density_kg_m3: 1390.0,
            tension_n_m: 2500.0,
            radial_intervals: 6,
            azimuths: 24,
        };
        let film = TensionedDisk::new(s, DISK_BUDGET).unwrap();
        let it = Scoped::new("part", "head");
        let e = format!(
            "{}",
            select_non_splitting_modes(&film, 2, 1.0e7, &it).unwrap_err()
        );
        assert!(e.contains("splits a degenerate pair") && e.contains("choose 1 or 3"));
        let ms = select_non_splitting_modes(&film, 3, 1.0e7, &it).unwrap();
        let sh = bend_shares(&film, &ms, &it).unwrap();
        assert!(sh.len() == 3 && sh[0] > 0.0 && sh[0] < 1.0e-3);
        let a0 = film.modal_area(&ms[0].phi).unwrap();
        let a1 = film.modal_area(&ms[1].phi).unwrap();
        assert!(a0.abs() > 0.1 && a1.abs() < 1.0e-12);
        let c0 = film.sample_mode_at(&ms[0].phi, 0.0, 0.0).unwrap();
        let c1 = film.sample_mode_at(&ms[1].phi, 0.0, 0.0).unwrap();
        assert!(c0.abs() > 1.0 && c1.abs() < 1.0e-12);
        assert!(film.sample_mode_at(&ms[0].phi, 0.2, 0.0).is_err());
    }
}
