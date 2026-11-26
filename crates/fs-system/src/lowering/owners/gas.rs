//! Closed gas cavities lowered onto `fs-gas::component::GasCavity`.
//!
//! Lowers `Cavity` parts (`{Fluid}` insulated, or `{Fluid, Thermal}`) from a
//! `Cylinder { radius, height }` and one material card (`molar_mass`, NASA-9,
//! transport claims). Attaches no hidden ambient pressure force.

use core::f64::consts::PI;

use fs_blake3::ContentHash;
use fs_gas::component::{GasCavity, GasEnergySolve, parcel_from_material_card};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QtyAny, QuantitySpec};
use fs_scenario::assembly::draft::{
    DraftGeometrySource, DraftInterfaceKind, DraftPartKind, FidelityErrorBudget, PhysicsDomain,
    PrimitiveShape,
};

use super::{Scoped, owning_part};
use crate::lowering::{
    AREA_DIMS, ConstitutiveModelRef, FidelityRungKind, InterfaceKindTag, InterfaceLoweringContext,
    InterfaceLoweringKey, LoweredInterfaceOutput, LoweredPartOutput, LoweringError, LoweringGraph,
    LoweringRegistry, PISTON_FACE_FAMILY, PartLoweringContext, PartLoweringKey, TEMPERATURE_DIMS,
    VOLUME_DIMS,
};

const P_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const E_DIMS: Dims = Dims([2, 0, -2, 0, 0, 0]);
const NASA9_LAW: &str = "nasa9-standard-state";
/// Minimum inversion tolerance [J/kg], above the ~5e-11 J/kg ULP of `c_v * T`.
const MIN_TOL_J_KG: f64 = 1.0e-9;

/// Register the gas-cavity part constructors and the piston-face interface constructor.
pub(crate) fn register(registry: &mut LoweringRegistry) {
    let f = PhysicsDomain::Fluid;
    for (d, iso) in [(vec![f], true), (vec![f, PhysicsDomain::Thermal], false)] {
        let k = PartLoweringKey::new(DraftPartKind::Cavity, d, FidelityRungKind::Lumped);
        registry.register_part(k, move |ctx, g| lower_cavity(ctx, g, iso));
    }
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::PistonFace, PISTON_FACE_FAMILY),
        lower_piston_face,
    );
}

fn lower_cavity(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
    insulated: bool,
) -> Result<LoweredPartOutput, LoweringError> {
    let mut it = Scoped::new("part", &ctx.part.name);
    let DraftGeometrySource::Primitive(PrimitiveShape::Cylinder { radius, height }) =
        &ctx.part.geometry
    else {
        return Err(it.invalid("gas cavity requires Cylinder { radius, height } geometry"));
    };
    let (r, h) = (radius.value, height.value);
    if !(r.is_finite() && r > 0.0 && h.is_finite() && h > 0.0) {
        return Err(it.invalid("cylinder radius and height must be positive finite"));
    }
    let v0 = PI * r * r * h;
    it.derive(ctx.part, "bore_area", PI * r * r, AREA_DIMS, "pi*r^2");
    it.derive(ctx.part, "volume", v0, VOLUME_DIMS, "pi*r^2*h");

    let [mat] = ctx.part.materials.as_slice() else {
        return Err(it.invalid("gas cavity requires exactly one material card"));
    };
    let card = ctx
        .cards
        .material(&mat.card)
        .ok_or(LoweringError::MissingCard {
            item: it.item.clone(),
            card_kind: "material",
            card_hash: mat.card,
        })?;
    let env = &ctx.draft.input().environment.base;
    let mut pick = |decl: Option<QtyAny>, eq: QtyAny, (f, ef, d)| {
        let q = decl.unwrap_or(eq);
        if decl.is_some() {
            it.declare(f, q);
        } else {
            it.declare_on("environment", ef, q);
        }
        (q.dims == d && q.value.is_finite() && q.value > 0.0)
            .then_some(q.value)
            .ok_or_else(|| it.invalid(format!("{f} must be positive finite")))
    };
    let t_d = ctx.part.initial.as_ref().and_then(|i| i.temperature);
    let tk = (
        "initial.temperature",
        "ambient_temperature",
        TEMPERATURE_DIMS,
    );
    let t0 = pick(t_d, env.ambient_temperature, tk)?;
    let p_q = mat
        .query_point
        .iter()
        .find(|(k, _)| k.name == "pressure" || k.name == "P");
    let pk = ("material.query_point.pressure", "ambient_pressure", P_DIMS);
    let p0 = pick(p_q.map(|(_, p)| *p), env.ambient_pressure, pk)?;

    let band = nasa9_band(card, mat.card, t0, &mut it)?;
    let qs = |k| QuantitySpec::semantic(SemanticType::new(k, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity("temperature", qs(QuantityKind::AbsoluteTemperature), t0)
        .and_then(|qp| qp.with_quantity("pressure", qs(QuantityKind::Pressure), p0))
        .map_err(|e| it.invalid(format!("gas query point refused: {e}")))?;
    for prop in ["molar_mass", "dynamic_viscosity", "thermal_conductivity"] {
        if mat.require_known(prop).is_err() {
            return Err(LoweringError::UnknownPropertyDefaulted {
                item: it.item,
                property: prop.into(),
            });
        }
        let ans = card
            .claims()
            .query(prop, &point, ctx.selection_policy)
            .map_err(|s| LoweringError::MatDb {
                item: it.item.clone(),
                source: Box::new(s),
            })?;
        it.r.record_matdb_claim(prop, &ans.evidence.value, mat.card, &ans.receipt);
    }

    let parcel = parcel_from_material_card(card, point, v0)
        .map_err(|e| it.invalid(format!("gas parcel refused: {e}")))?;
    let controls = solve_controls(&ctx.part.fidelity.error_budget, band, &mut it)?;
    let cavity = GasCavity::new(parcel.clone(), controls)
        .and_then(|c| if insulated { c.insulated() } else { Ok(c) })
        .map_err(|e| it.invalid(format!("GasCavity refused: {e}")))?;
    let node = graph.add_with_tangent(&ctx.part.name, Some(ctx.part.entity), cavity, parcel)?;
    graph.expose_port(ctx.part.entity, "pressure", graph.port(node, "p")?);
    if !insulated {
        graph.expose_port(ctx.part.entity, "thermal", graph.port(node, "heat")?);
    }
    Ok(it.part(&ctx.part.name))
}

fn solve_controls(
    b: &FidelityErrorBudget,
    (lo, hi): (f64, f64),
    it: &mut Scoped,
) -> Result<GasEnergySolve, LoweringError> {
    let q = b.absolute_tol;
    if !((q.dims == Dims::NONE || q.dims == E_DIMS) && q.value.is_finite() && q.value > 0.0) {
        return Err(it.invalid("error_budget.absolute_tol must be positive [1] or [J/kg]"));
    }
    let tol = q.value.max(MIN_TOL_J_KG);
    let clamp = |u| usize::try_from(u).unwrap_or(128).clamp(1, 128);
    let iters = b.max_work_units.map_or(128, clamp);
    let tol_q = QtyAny::new(tol, E_DIMS);
    it.declare("fidelity.error_budget.absolute_tol", tol_q);
    #[allow(clippy::cast_precision_loss)]
    let max_q = QtyAny::dimensionless(iters as f64);
    it.declare("fidelity.error_budget.max_work_units", max_q);
    Ok(GasEnergySolve {
        minimum_temperature_k: lo,
        maximum_temperature_k: hi,
        energy_tolerance_j_kg: tol,
        maximum_iterations: iters,
    })
}

fn nasa9_band(
    card: &MaterialCard,
    hash: ContentHash,
    t0: f64,
    it: &mut Scoped,
) -> Result<(f64, f64), LoweringError> {
    let mut bands = Vec::new();
    for (i, m) in card
        .models()
        .iter()
        .filter(|m| m.law.0 == NASA9_LAW)
        .enumerate()
    {
        let mref = |p| ConstitutiveModelRef {
            card_hash: hash,
            model_hash: m.content_hash(),
            law: &m.law.0,
            law_version: m.law_version,
            parameter: p,
        };
        for (k, p) in &m.parameters {
            it.r.record_model_parameter(format!("nasa9[{i}].{k}"), p.value, p.dims, mref(k));
        }
        if let Some((lo, hi)) = m.validity.bound("T") {
            bands.push((lo, hi));
            for (nm, v) in [("validity.T.lower", lo), ("validity.T.upper", hi)] {
                it.r.record_model_parameter(
                    format!("nasa9[{i}].{nm}"),
                    v,
                    TEMPERATURE_DIMS,
                    mref(nm),
                );
            }
        }
    }
    bands.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = bands.len();
    let Some((mut lo, mut hi)) = bands.iter().enumerate().find_map(|(i, &(a, b))| {
        (t0 >= a && (t0 < b || (i + 1 == n && t0 <= b)))
            .then_some((a, if i + 1 < n { b.next_down() } else { b }))
    }) else {
        return Err(it.invalid("initial temperature falls outside every NASA-9 region"));
    };
    for (_, c) in card.claims().claims_ordered() {
        if let Some((a, b)) = c
            .validity
            .bound("temperature")
            .or_else(|| c.validity.bound("T"))
        {
            (lo, hi) = (lo.max(a), hi.min(b));
        }
    }
    if !(lo.is_finite() && hi.is_finite() && lo > 0.0 && lo < hi && t0 >= lo && t0 <= hi) {
        return Err(it.invalid("gas card has no valid overlapping temperature band containing T0"));
    }
    Ok((lo, hi))
}

fn lower_piston_face(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let mut it = Scoped::new("interface", &ctx.interface.name);
    let DraftInterfaceKind::PistonFace { state } = &ctx.interface.kind else {
        return Err(it.invalid("expected PistonFace interface"));
    };
    let (mut area_opt, mut sign_opt) = (None, None);
    for (k, q) in state {
        match k.name.as_str() {
            "area" if area_opt.is_none() => area_opt = Some(*q),
            "sign" if sign_opt.is_none() => sign_opt = Some(*q),
            "area" | "sign" => return Err(it.invalid(format!("duplicate `{}`", k.name))),
            other => return Err(it.invalid(format!("unknown piston-face parameter `{other}`"))),
        }
    }
    let inp = ctx.draft.input();
    let (cat, fa, fb) = (
        &inp.catalog,
        ctx.interface.feature_a,
        ctx.interface.feature_b,
    );
    let port = |e, n| graph.resolve_entity_port(cat, e, n);
    let (mech, press, cav) = match (
        port(fa, "mechanics"),
        port(fb, "pressure"),
        port(fb, "mechanics"),
        port(fa, "pressure"),
    ) {
        (Some(m), Some(p), None, None) => (m, p, fb),
        (None, None, Some(m), Some(p)) => (m, p, fa),
        _ => return Err(it.invalid("piston-face requires `mechanics` and `pressure` ports")),
    };
    let area = if let Some(q) = area_opt {
        if !(q.dims == AREA_DIMS && q.value.is_finite() && q.value > 0.0) {
            return Err(it.invalid("piston-face `area` must be positive [m^2]"));
        }
        it.declare("state.area", q);
        q.value
    } else {
        let part = owning_part(inp, cav)
            .map(|i| &inp.parts[i])
            .ok_or_else(|| it.invalid("missing cavity part"))?;
        let DraftGeometrySource::Primitive(PrimitiveShape::Cylinder { radius, .. }) =
            &part.geometry
        else {
            return Err(it.invalid("piston-face without `area` requires cylinder cavity"));
        };
        if !(radius.value.is_finite() && radius.value > 0.0) {
            return Err(it.invalid("cavity cylinder radius must be positive finite"));
        }
        let a = PI * radius.value * radius.value;
        it.derive(part, "area", a, AREA_DIMS, "pi*r^2");
        a
    };
    let sign = match sign_opt {
        Some(q) if q.dims == Dims::NONE && q.value.abs().total_cmp(&1.0).is_eq() => {
            it.declare("state.sign", q);
            q.value
        }
        Some(_) => return Err(it.invalid("piston-face `sign` must be +1 or -1")),
        None => {
            it.declare("default_sign", QtyAny::dimensionless(1.0));
            1.0
        }
    };
    graph.attach_common_flow(mech, press, sign * area);
    Ok(it.interface_output(None))
}

#[cfg(test)]
mod tests {
    use std::{collections::BTreeSet, path::Path};

    use fs_gas::component::ClosedGasParcel;
    use fs_matdb::NormalizedMaterialCardPack;
    use fs_scenario::assembly::draft::*;
    use fs_scenario::assembly::{MaterialRef, ReferenceState};
    use fs_scenario::{EntityCatalog, EntityDeclaration, FrameId, Scenario};

    use super::*;
    use crate::lowering::owners::register_owner_defaults;
    use crate::lowering::{CardStore, LENGTH_DIMS};

    #[test]
    fn cavity_lowers_insulated_and_thermal_and_refuses_bad_geometry_or_band() {
        let p = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk");
        let pack = NormalizedMaterialCardPack::from_bytes(&std::fs::read(p).unwrap()).unwrap();
        let card = pack.card().clone();
        let lower_one = |geom, thermal, t0| {
            let mut cards = CardStore::new();
            let n2 = cards.insert_material(card.clone());
            let mut cat = EntityCatalog::new();
            let root = cat.declare(EntityDeclaration::assembly("cav")).unwrap();
            let eg = cat.declare(EntityDeclaration::part(root, "gas")).unwrap();
            let physics = if thermal {
                BTreeSet::from([PhysicsDomain::Fluid, PhysicsDomain::Thermal])
            } else {
                BTreeSet::from([PhysicsDomain::Fluid])
            };
            let mut env = DraftEnvironment::earth_lab_air();
            env.base.ambient_temperature = QtyAny::new(t0, TEMPERATURE_DIMS);
            env.base.ambient_pressure = QtyAny::new(90_000.0, P_DIMS);
            let sc = Scenario::new("cav", 1, env.base.clone());
            let mat = MaterialRef {
                card: n2,
                query_point: vec![],
                frame: FrameId(0),
                layer: None,
                unknown: vec![],
            };
            let part = DraftPartDecl {
                name: "gas".into(),
                entity: eg,
                kind: DraftPartKind::Cavity,
                geometry: DraftGeometrySource::Primitive(geom),
                materials: vec![mat],
                physics,
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::Lumped],
                    error_budget: FidelityErrorBudget {
                        relative_tol: 1e-6,
                        absolute_tol: QtyAny::dimensionless(1e-12),
                        max_work_units: Some(128),
                    },
                },
                uncertainty: vec![],
                initial: None,
            };
            let s = Dims([0, 0, 1, 0, 0, 0]);
            let time_plan = TimePlan {
                horizon: QtyAny::new(1e-3, s),
                step_policy: StepPolicy::Fixed {
                    dt: QtyAny::new(1e-4, s),
                },
                default_sample_rate: QtyAny::new(1e3, Dims([0, 0, -1, 0, 0, 0])),
                scheme: TimeScheme::Explicit,
            };
            let mk_src = |n: &str, p: &str, d| SourceDecl {
                name: n.into(),
                target: SourceTarget::Port {
                    entity: eg,
                    port: p.into(),
                },
                history: TimeHistory::Constant(QtyAny::new(0.0, d)),
            };
            let mut sources = vec![mk_src("flow", "pressure", Dims([3, 0, -1, 0, 0, 0]))];
            if thermal {
                sources.push(mk_src("heat", "thermal", Dims([2, 1, -3, 0, 0, 0])));
            }
            let inp = DraftAssemblyInput {
                name: "cav".into(),
                root,
                catalog: cat,
                reference_state: Some(ReferenceState::Assembled),
                environment: env,
                time_plan,
                parts: vec![part],
                joints: vec![],
                gear_meshes: vec![],
                interfaces: vec![],
                sources,
                controllers: vec![],
                probes: vec![],
            };
            let mut reg = LoweringRegistry::with_kernel_defaults();
            register_owner_defaults(&mut reg);
            reg.lower(&DraftAssembly::admit(inp, &sc).unwrap(), &cards)
        };
        let r = QtyAny::new(0.15, LENGTH_DIMS);
        let cyl = PrimitiveShape::Cylinder {
            radius: r,
            height: QtyAny::new(0.20, LENGTH_DIMS),
        };
        let iso = lower_one(cyl.clone(), false, 300.0).unwrap();
        iso.receipt.verify_completeness().unwrap();
        assert_eq!(iso.receipt.components[0].property_receipts.len(), 3);
        let gas = iso.node_map["gas"];
        let vol = iso
            .system
            .state::<ClosedGasParcel>(gas)
            .unwrap()
            .volume_m3();
        assert!((vol - PI * 0.15 * 0.15 * 0.20).abs() < 1e-15);
        assert!(iso.system.port(gas, "p").is_some() && iso.system.port(gas, "heat").is_none());
        let th = lower_one(cyl.clone(), true, 300.0).unwrap();
        th.receipt.verify_completeness().unwrap();
        assert!(th.system.port(th.node_map["gas"], "heat").is_some());
        let bad = |g, t| matches!(lower_one(g, false, t), Err(LoweringError::Invalid { .. }));
        assert!(bad(PrimitiveShape::Sphere { radius: r }, 300.0) && bad(cyl, 1500.0));
    }
}
