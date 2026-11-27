//! Normal Hertzian contact lowered onto `fs-dcontact::component::ContactInterface`.

use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_dcontact::{Obstacle, hunt_crossley_loss_from_restitution};
use fs_qty::Dims;
use fs_scenario::assembly::draft::{
    DraftGeometrySource, DraftInterfaceKind, DraftPartDecl, PrimitiveShape,
};

use crate::lowering::{
    FORCE_DIMS, HERTZ_CONTACT_FAMILY, InterfaceKindTag, InterfaceLoweringContext,
    InterfaceLoweringKey, LENGTH_DIMS, LoweredInterfaceOutput, LoweringError, LoweringGraph,
    LoweringRegistry, PRESSURE_DIMS, ProvenanceRecorder, VELOCITY_DIMS, attach_mech_endpoint,
    geometry_generator_label, owning_part_decl, resolve_geometry_fingerprint,
    route_dissipated_heat,
};

/// Register the card-less `HertzContact` interface constructor.
pub(crate) fn register(registry: &mut LoweringRegistry) {
    let key = InterfaceLoweringKey::new(InterfaceKindTag::HertzContact, HERTZ_CONTACT_FAMILY);
    registry.register_interface(key, lower_hertz_contact);
}

fn inv(scope: &str, detail: impl Into<String>) -> LoweringError {
    LoweringError::Invalid {
        item: scope.to_owned(),
        detail: detail.into(),
    }
}

fn inv_radius(geom: &DraftGeometrySource, scope: &str) -> Result<f64, LoweringError> {
    match geom {
        DraftGeometrySource::Primitive(PrimitiveShape::Sphere { radius })
            if radius.value.is_finite() && radius.value > 0.0 =>
        {
            Ok(1.0 / radius.value)
        }
        DraftGeometrySource::Primitive(PrimitiveShape::Box { .. }) => Ok(0.0),
        _ => Err(inv(scope, "HertzContact requires Sphere or Box primitive")),
    }
}

fn elastic_pair(
    ctx: &InterfaceLoweringContext<'_>,
    parts: (&DraftPartDecl, &DraftPartDecl),
    rec: &mut ProvenanceRecorder,
    scope: &str,
) -> Result<f64, LoweringError> {
    let env = ctx.query_env();
    let mut one = |part: &DraftPartDecl, side: &str| {
        let mr = part
            .materials
            .first()
            .ok_or_else(|| inv(scope, format!("{side} has no material")))?;
        let mut q = |p: &str| {
            rec.query_material_property(scope, ctx.cards, mr, p, &format!("{side}.{p}"), env)
                .map(|(s, _)| s)
        };
        let (e, nu) = (q("young_modulus")?, q("poisson_ratio")?);
        if e.dims != PRESSURE_DIMS || !(e.value.is_finite() && e.value > 0.0) {
            return Err(inv(scope, format!("{side}.young_modulus invalid [Pa]")));
        }
        if nu.dims != Dims::NONE || !(nu.value.is_finite() && nu.value > -1.0 && nu.value < 0.5) {
            return Err(inv(scope, format!("{side}.poisson_ratio out of (-1,0.5)")));
        }
        Ok((e.value, nu.value))
    };
    let ((ea, nua), (eb, nub)) = (one(parts.0, "side_a")?, one(parts.1, "side_b")?);
    Ok(1.0 / ((1.0 - nua * nua) / ea + (1.0 - nub * nub) / eb))
}

fn build_hertz(
    ctx: &InterfaceLoweringContext<'_>,
    (e_star, r_star): (f64, f64),
    rec: &mut ProvenanceRecorder,
    scope: &str,
) -> Result<(ContactInterface, f64, f64), LoweringError> {
    let DraftInterfaceKind::HertzContact { state } = &ctx.interface.kind else {
        return Err(inv(scope, "expected HertzContact interface"));
    };
    let get = |n| {
        state
            .iter()
            .find(|(k, _)| k.name == n)
            .map(|(_, q)| q.value)
    };
    let max_f = get("maximum_force").ok_or_else(|| inv(scope, "missing maximum_force"))?;
    let max_p =
        get("maximum_penetration").ok_or_else(|| inv(scope, "missing maximum_penetration"))?;
    let (gap, gap_src) = match get("gap") {
        Some(v) => (v, "state.gap"),
        None => (0.0, "not declared: 0"),
    };
    let (init_x, init_src) = match get("initial_closure") {
        Some(v) => (v, "state.initial_closure"),
        None => (0.0, "not declared: 0"),
    };
    let k_h = (4.0 / 3.0) * e_star * r_star.sqrt();
    let ef = "1/((1-nu_a^2)/E_a+(1-nu_b^2)/E_b)";
    for (k, v, d, f) in [
        ("effective_modulus", e_star, PRESSURE_DIMS, ef),
        ("hertz_exponent", 1.5, Dims::NONE, "Hertz spherical law"),
        ("maximum_force", max_f, FORCE_DIMS, "state.maximum_force"),
        (
            "maximum_penetration",
            max_p,
            LENGTH_DIMS,
            "state.maximum_penetration",
        ),
        ("gap", gap, LENGTH_DIMS, gap_src),
        ("initial_closure", init_x, LENGTH_DIMS, init_src),
    ] {
        rec.record_declared_field(k, v, d, scope, f);
    }
    let chi = match (get("restitution"), get("approach_speed")) {
        (Some(e), Some(v0)) => {
            rec.record_declared_field("restitution", e, Dims::NONE, scope, "state.restitution");
            rec.record_declared_field(
                "approach_speed",
                v0,
                VELOCITY_DIMS,
                scope,
                "state.approach_speed",
            );
            hunt_crossley_loss_from_restitution(e, v0).map_err(|err| inv(scope, err.to_string()))?
        }
        (None, None) => 0.0,
        _ => return Err(inv(scope, "restitution and approach_speed pair required")),
    };
    let nm = &ctx.interface.name;
    let prov = format!("hertz-contact:{nm}:E*={e_star:.6e},R*={r_star:.6e}");
    let mut law = Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], k_h, 1.5, prov)
        .map_err(|err| inv(scope, err.to_string()))?;
    if chi > 0.0 {
        law = law
            .with_internal_loss(chi)
            .map_err(|e| inv(scope, e.to_string()))?;
    }
    let solve = ContactSolve {
        max_iterations: 64,
        maximum_force_n: max_f,
        maximum_penetration_m: max_p,
        force_absolute_tolerance_n: 1.0e-12,
        force_relative_tolerance: 1.0e-12,
    };
    let iface = ContactInterface::new(law, solve)
        .and_then(ContactInterface::publishing_reaction)
        .map_err(|err| inv(scope, err.to_string()))?;
    Ok((iface, init_x, chi))
}

fn lower_hertz_contact(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let scope = format!("interface:{}", ctx.interface.name);
    let inp = ctx.draft.input();
    let (fa, fb) = (ctx.interface.feature_a, ctx.interface.feature_b);
    let part = |f, s| owning_part_decl(inp, f).ok_or_else(|| inv(&scope, s));
    let (pa, pb) = (
        part(fa, "feature_a has no part")?,
        part(fb, "feature_b has no part")?,
    );
    let mut rec = ProvenanceRecorder::default();
    let e_star = elastic_pair(ctx, (pa, pb), &mut rec, &scope)?;
    let (ir_a, ir_b) = (
        inv_radius(&pa.geometry, &scope)?,
        inv_radius(&pb.geometry, &scope)?,
    );
    let ir_star = ir_a + ir_b;
    if !(ir_star.is_finite() && ir_star > 0.0) {
        return Err(inv(&scope, "HertzContact refuses flat-on-flat geometry"));
    }
    let r_star = 1.0 / ir_star;
    let curved = if ir_a > 0.0 { pa } else { pb };
    rec.record_geometry_derived(
        "effective_radius",
        r_star,
        LENGTH_DIMS,
        resolve_geometry_fingerprint(&curved.geometry),
        geometry_generator_label(&curved.geometry),
        "1/(1/Ra+1/Rb)",
    );
    let (contact, init_x, chi) = build_hertz(ctx, (e_star, r_star), &mut rec, &scope)?;
    let name = format!("interface.{}", ctx.interface.name);
    let st = ContactState::at(init_x);
    let node = graph.add_on_line(&name, Some(ctx.interface.entity), contact, st)?;
    for (feat, pname) in [(fa, "a"), (fb, "b")] {
        let p = graph.port(node, pname)?;
        attach_mech_endpoint(graph, &inp.catalog, feat, inp.root, p, &scope)?;
    }
    let rn = graph.port(node, "reaction_n")?;
    graph.expose_port(ctx.interface.entity, "reaction_n", rn);
    if chi > 0.0 {
        let (h, t) = (graph.port(node, "heat")?, ctx.ambient_temperature_k());
        route_dissipated_heat(
            graph,
            &inp.catalog,
            inp.root,
            h,
            fa,
            Some(fb),
            t,
            &scope,
            &mut rec,
        )?;
    }
    let (property_receipts, parameter_sources) = rec.finish();
    Ok(LoweredInterfaceOutput {
        created_nodes: vec![name],
        property_receipts,
        parameter_sources,
    })
}
