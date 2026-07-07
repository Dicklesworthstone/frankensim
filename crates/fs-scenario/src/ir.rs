//! Canonical scenario IR: a deterministic s-expression encoding with
//! LOSSLESS round-trip (floats print in shortest-round-trip form; dims
//! travel as explicit SI exponent vectors, so no unit-string parsing is
//! involved). `write_ir` output is byte-stable — the ledger stores it as
//! a content-addressed artifact; `parse_ir` inverts it exactly.

use crate::ScenarioError;
use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
use crate::ensemble::{SpectrumModel, StochasticEnsemble};
use crate::frame::{Frame, FrameId, FrameMotion, FrameTree};
use crate::scenario::{
    Combination, ContactLaw, ContactModel, Environment, LoadCase, Scenario, Violation,
};
use crate::signal::{ChebProfile, Interp, TimeSignal};
use fs_cheb::Cheb1;
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use std::fmt::Write as _;

// ---------------------------------------------------------------- writer

fn w_qty(out: &mut String, q: &QtyAny) {
    let _ = write!(out, "(qty {}", q.value);
    for d in q.dims.0 {
        let _ = write!(out, " {d}");
    }
    out.push(')');
}

fn w_dims(out: &mut String, d: Dims) {
    out.push_str("(dims");
    for e in d.0 {
        let _ = write!(out, " {e}");
    }
    out.push(')');
}

fn w_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        if c == '"' || c == '\\' {
            out.push('\\');
        }
        out.push(c);
    }
    out.push('"');
}

fn w_vec3(out: &mut String, v: Vec3) {
    let _ = write!(out, "(vec {} {} {})", v.x, v.y, v.z);
}

fn w_floats(out: &mut String, head: &str, vs: &[f64]) {
    let _ = write!(out, "({head}");
    for v in vs {
        let _ = write!(out, " {v}");
    }
    out.push(')');
}

fn w_signal(out: &mut String, s: &TimeSignal) {
    match s {
        TimeSignal::Constant(q) => {
            out.push_str("(constant ");
            w_qty(out, q);
            out.push(')');
        }
        TimeSignal::Ramp {
            t_start,
            t_end,
            from,
            to,
        } => {
            let _ = write!(out, "(ramp {t_start} {t_end} ");
            w_qty(out, from);
            out.push(' ');
            w_qty(out, to);
            out.push(')');
        }
        TimeSignal::Table {
            times,
            values,
            dims,
            interp,
        } => {
            let tag = match interp {
                Interp::Linear => "linear",
                Interp::Hold => "hold",
            };
            let _ = write!(out, "(table {tag} ");
            w_dims(out, *dims);
            out.push(' ');
            w_floats(out, "times", times);
            out.push(' ');
            w_floats(out, "values", values);
            out.push(')');
        }
        TimeSignal::Chebfun(p) => {
            w_profile(out, "chebfun", p);
        }
    }
}

fn w_profile(out: &mut String, head: &str, p: &ChebProfile) {
    let (a, b) = p.cheb.domain();
    let _ = write!(out, "({head} ");
    w_dims(out, p.dims);
    let _ = write!(out, " {a} {b} ");
    w_floats(out, "coeffs", p.cheb.coeffs());
    out.push(')');
}

fn physics_tag(p: Physics) -> &'static str {
    match p {
        Physics::IncompressibleFlow => "incompressible-flow",
        Physics::Thermal => "thermal",
        Physics::Elasticity => "elasticity",
    }
}

fn kind_tag(k: BcKind) -> &'static str {
    match k {
        BcKind::Dirichlet => "dirichlet",
        BcKind::Neumann => "neumann",
        BcKind::Robin => "robin",
        BcKind::MassFlowInlet => "mass-flow-inlet",
        BcKind::PressureOutlet => "pressure-outlet",
        BcKind::WallNoSlip => "wall-no-slip",
        BcKind::WallSlip => "wall-slip",
        BcKind::Traction => "traction",
    }
}

fn w_bc(out: &mut String, bc: &BoundaryCondition) {
    out.push_str("(bc ");
    w_str(out, &bc.region);
    let _ = write!(
        out,
        " {} {} {} ",
        physics_tag(bc.physics),
        kind_tag(bc.kind),
        bc.frame
    );
    match &bc.value {
        None => out.push_str("none"),
        Some(BcValue::Uniform(q)) => {
            out.push_str("(uniform ");
            w_qty(out, q);
            out.push(')');
        }
        Some(BcValue::Signal(s)) => {
            out.push_str("(signal ");
            w_signal(out, s);
            out.push(')');
        }
        Some(BcValue::Profile(p)) => w_profile(out, "profile", p),
    }
    out.push(' ');
    match bc.compatibility {
        None => out.push_str("none"),
        Some(Compat::Incompressible) => out.push_str("incompressible"),
    }
    out.push(')');
}

fn w_frame(out: &mut String, f: &Frame) {
    let _ = write!(out, "(frame {} ", f.id.0);
    w_str(out, &f.name);
    let _ = write!(out, " {} ", f.parent.0);
    match &f.motion {
        FrameMotion::Fixed {
            orientation,
            translation,
        } => {
            let _ = write!(
                out,
                "(fixed (quat {} {} {} {}) ",
                orientation.w, orientation.x, orientation.y, orientation.z
            );
            w_vec3(out, *translation);
            out.push(')');
        }
        FrameMotion::Rotating { axis, center, rate } => {
            let _ = write!(out, "(rotating (vec {} {} {}) ", axis[0], axis[1], axis[2]);
            w_vec3(out, *center);
            out.push(' ');
            w_qty(out, rate);
            out.push(')');
        }
        FrameMotion::Tilt {
            axis,
            center,
            angle,
        } => {
            let _ = write!(out, "(tilt (vec {} {} {}) ", axis[0], axis[1], axis[2]);
            w_vec3(out, *center);
            out.push(' ');
            w_signal(out, angle);
            out.push(')');
        }
    }
    out.push(')');
}

fn w_ensemble(out: &mut String, e: &StochasticEnsemble) {
    out.push_str("(ensemble ");
    w_str(out, &e.name);
    let _ = write!(out, " {} {} ", e.seed, e.members);
    w_qty(out, &e.duration);
    out.push(' ');
    w_qty(out, &e.dt);
    out.push(' ');
    match &e.model {
        SpectrumModel::Dryden {
            sigma,
            length_scale,
            mean_speed,
        } => {
            out.push_str("(dryden ");
            w_qty(out, sigma);
            out.push(' ');
            w_qty(out, length_scale);
            out.push(' ');
            w_qty(out, mean_speed);
            out.push(')');
        }
        SpectrumModel::KanaiTajimi {
            s0,
            omega_g,
            zeta_g,
        } => {
            let _ = write!(out, "(kanai-tajimi {s0} ");
            w_qty(out, omega_g);
            let _ = write!(out, " {zeta_g})");
        }
        SpectrumModel::CarreauBand {
            eta_zero,
            eta_inf,
            lambda,
            n,
        } => {
            out.push_str("(carreau");
            for q in [
                &eta_zero[0],
                &eta_zero[1],
                &eta_inf[0],
                &eta_inf[1],
                &lambda[0],
                &lambda[1],
            ] {
                out.push(' ');
                w_qty(out, q);
            }
            let _ = write!(out, " {} {})", n[0], n[1]);
        }
    }
    out.push(')');
}

/// Serialize a scenario to canonical, byte-stable IR text.
#[must_use]
pub fn write_ir(s: &Scenario) -> String {
    let mut out = String::with_capacity(1024);
    out.push_str("(scenario ");
    w_str(&mut out, &s.name);
    let _ = write!(out, " {} (environment ", s.seed);
    for g in &s.environment.gravity {
        w_qty(&mut out, g);
        out.push(' ');
    }
    w_qty(&mut out, &s.environment.ambient_temperature);
    out.push(' ');
    w_qty(&mut out, &s.environment.ambient_pressure);
    out.push_str(") (frames");
    for f in &s.frames.frames {
        out.push(' ');
        w_frame(&mut out, f);
    }
    out.push_str(") (bcs");
    for bc in &s.base_bcs {
        out.push(' ');
        w_bc(&mut out, bc);
    }
    out.push_str(") (cases");
    for case in &s.cases {
        out.push_str(" (case ");
        w_str(&mut out, &case.name);
        for bc in &case.bcs {
            out.push(' ');
            w_bc(&mut out, bc);
        }
        out.push(')');
    }
    out.push_str(") (combos");
    for combo in &s.combinations {
        out.push_str(" (combo ");
        w_str(&mut out, &combo.name);
        for (case, factor) in &combo.terms {
            out.push_str(" (term ");
            w_str(&mut out, case);
            let _ = write!(out, " {factor})");
        }
        out.push(')');
    }
    out.push_str(") (ensembles");
    for e in &s.ensembles {
        out.push(' ');
        w_ensemble(&mut out, e);
    }
    out.push_str(") (contacts");
    for c in &s.contacts {
        out.push_str(" (contact ");
        w_str(&mut out, &c.region_a);
        out.push(' ');
        w_str(&mut out, &c.region_b);
        out.push(' ');
        match c.model {
            ContactModel::Frictionless => out.push_str("frictionless"),
            ContactModel::Tied => out.push_str("tied"),
            ContactModel::Coulomb { mu_s, mu_k } => {
                let _ = write!(out, "(coulomb {mu_s} {mu_k})");
            }
        }
        out.push(')');
    }
    out.push_str("))");
    out
}

// ---------------------------------------------------------------- parser

#[derive(Debug, Clone, PartialEq)]
enum Sx {
    Atom(String),
    Str(String),
    List(Vec<Sx>),
}

fn err(at: usize, what: &str) -> ScenarioError {
    ScenarioError::Parse {
        at,
        what: what.to_string(),
    }
}

fn parse_sx(text: &str) -> Result<Sx, ScenarioError> {
    let bytes = text.as_bytes();
    let mut pos = 0usize;
    let root = parse_one(bytes, &mut pos)?;
    skip_ws(bytes, &mut pos);
    if pos != bytes.len() {
        return Err(err(pos, "trailing bytes after the scenario form"));
    }
    Ok(root)
}

fn skip_ws(bytes: &[u8], pos: &mut usize) {
    while *pos < bytes.len() && bytes[*pos].is_ascii_whitespace() {
        *pos += 1;
    }
}

fn parse_one(bytes: &[u8], pos: &mut usize) -> Result<Sx, ScenarioError> {
    skip_ws(bytes, pos);
    match bytes.get(*pos) {
        None => Err(err(*pos, "unexpected end of input")),
        Some(b'(') => {
            *pos += 1;
            let mut items = Vec::new();
            loop {
                skip_ws(bytes, pos);
                match bytes.get(*pos) {
                    None => return Err(err(*pos, "unclosed list")),
                    Some(b')') => {
                        *pos += 1;
                        return Ok(Sx::List(items));
                    }
                    _ => items.push(parse_one(bytes, pos)?),
                }
            }
        }
        Some(b')') => Err(err(*pos, "unexpected ')'")),
        Some(b'"') => {
            *pos += 1;
            let mut s = String::new();
            loop {
                match bytes.get(*pos) {
                    None => return Err(err(*pos, "unterminated string")),
                    Some(b'\\') => {
                        let c = *bytes.get(*pos + 1).ok_or_else(|| err(*pos, "bad escape"))?;
                        s.push(c as char);
                        *pos += 2;
                    }
                    Some(b'"') => {
                        *pos += 1;
                        return Ok(Sx::Str(s));
                    }
                    Some(&c) => {
                        s.push(c as char);
                        *pos += 1;
                    }
                }
            }
        }
        Some(_) => {
            let start = *pos;
            while *pos < bytes.len()
                && !bytes[*pos].is_ascii_whitespace()
                && bytes[*pos] != b'('
                && bytes[*pos] != b')'
            {
                *pos += 1;
            }
            Ok(Sx::Atom(
                String::from_utf8_lossy(&bytes[start..*pos]).into_owned(),
            ))
        }
    }
}

// -------------------------------------------------------------- decoding

fn as_list<'a>(sx: &'a Sx, head: &str) -> Result<&'a [Sx], ScenarioError> {
    if let Sx::List(items) = sx
        && let Some(Sx::Atom(a)) = items.first()
        && a == head
    {
        return Ok(&items[1..]);
    }
    Err(err(0, &format!("expected ({head} ...) form")))
}

fn as_f64(sx: &Sx) -> Result<f64, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<f64>()
    {
        return Ok(v);
    }
    Err(err(0, "expected a number"))
}

fn as_str(sx: &Sx) -> Result<String, ScenarioError> {
    if let Sx::Str(s) = sx {
        return Ok(s.clone());
    }
    Err(err(0, "expected a string"))
}

fn as_atom(sx: &Sx) -> Result<&str, ScenarioError> {
    if let Sx::Atom(a) = sx {
        return Ok(a);
    }
    Err(err(0, "expected an atom"))
}

fn as_u64(sx: &Sx) -> Result<u64, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<u64>()
    {
        return Ok(v);
    }
    Err(err(0, "expected an unsigned integer"))
}

fn as_u32(sx: &Sx) -> Result<u32, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<u32>()
    {
        return Ok(v);
    }
    Err(err(0, "expected a u32"))
}

fn as_i8(sx: &Sx) -> Result<i8, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<i8>()
    {
        return Ok(v);
    }
    Err(err(0, "expected a small integer exponent"))
}

fn as_dims_at(items: &[Sx]) -> Result<Dims, ScenarioError> {
    if items.len() != 5 {
        return Err(err(0, "expected five dimension exponents"));
    }
    let mut d = [0i8; 5];
    for (slot, sx) in d.iter_mut().zip(items) {
        *slot = as_i8(sx)?;
    }
    Ok(Dims(d))
}

fn as_qty(sx: &Sx) -> Result<QtyAny, ScenarioError> {
    let items = as_list(sx, "qty")?;
    if items.len() != 6 {
        return Err(err(0, "qty needs value + five exponents"));
    }
    Ok(QtyAny::new(as_f64(&items[0])?, as_dims_at(&items[1..])?))
}

fn as_dims(sx: &Sx) -> Result<Dims, ScenarioError> {
    as_dims_at(as_list(sx, "dims")?)
}

fn as_floats(sx: &Sx, head: &str) -> Result<Vec<f64>, ScenarioError> {
    as_list(sx, head)?.iter().map(as_f64).collect()
}

fn as_vec3(sx: &Sx) -> Result<Vec3, ScenarioError> {
    let items = as_list(sx, "vec")?;
    if items.len() != 3 {
        return Err(err(0, "vec needs three components"));
    }
    Ok(Vec3::new(
        as_f64(&items[0])?,
        as_f64(&items[1])?,
        as_f64(&items[2])?,
    ))
}

fn as_profile(items: &[Sx]) -> Result<ChebProfile, ScenarioError> {
    if items.len() != 4 {
        return Err(err(0, "profile needs dims, domain, coeffs"));
    }
    let dims = as_dims(&items[0])?;
    let a = as_f64(&items[1])?;
    let b = as_f64(&items[2])?;
    let coeffs = as_floats(&items[3], "coeffs")?;
    Ok(ChebProfile {
        cheb: Cheb1::from_coeffs(a, b, coeffs),
        dims,
    })
}

fn as_signal(sx: &Sx) -> Result<TimeSignal, ScenarioError> {
    let Sx::List(items) = sx else {
        return Err(err(0, "expected a signal form"));
    };
    let head = as_atom(items.first().ok_or_else(|| err(0, "empty signal"))?)?;
    let rest = &items[1..];
    match head {
        "constant" => Ok(TimeSignal::Constant(as_qty(&rest[0])?)),
        "ramp" => {
            if rest.len() != 4 {
                return Err(err(0, "ramp needs t0 t1 from to"));
            }
            Ok(TimeSignal::Ramp {
                t_start: as_f64(&rest[0])?,
                t_end: as_f64(&rest[1])?,
                from: as_qty(&rest[2])?,
                to: as_qty(&rest[3])?,
            })
        }
        "table" => {
            if rest.len() != 4 {
                return Err(err(0, "table needs interp dims times values"));
            }
            let interp = match as_atom(&rest[0])? {
                "linear" => Interp::Linear,
                "hold" => Interp::Hold,
                other => return Err(err(0, &format!("unknown interp {other:?}"))),
            };
            Ok(TimeSignal::Table {
                interp,
                dims: as_dims(&rest[1])?,
                times: as_floats(&rest[2], "times")?,
                values: as_floats(&rest[3], "values")?,
            })
        }
        "chebfun" => Ok(TimeSignal::Chebfun(as_profile(rest)?)),
        other => Err(err(0, &format!("unknown signal {other:?}"))),
    }
}

fn as_physics(a: &str) -> Result<Physics, ScenarioError> {
    match a {
        "incompressible-flow" => Ok(Physics::IncompressibleFlow),
        "thermal" => Ok(Physics::Thermal),
        "elasticity" => Ok(Physics::Elasticity),
        other => Err(err(0, &format!("unknown physics {other:?}"))),
    }
}

fn as_kind(a: &str) -> Result<BcKind, ScenarioError> {
    match a {
        "dirichlet" => Ok(BcKind::Dirichlet),
        "neumann" => Ok(BcKind::Neumann),
        "robin" => Ok(BcKind::Robin),
        "mass-flow-inlet" => Ok(BcKind::MassFlowInlet),
        "pressure-outlet" => Ok(BcKind::PressureOutlet),
        "wall-no-slip" => Ok(BcKind::WallNoSlip),
        "wall-slip" => Ok(BcKind::WallSlip),
        "traction" => Ok(BcKind::Traction),
        other => Err(err(0, &format!("unknown bc kind {other:?}"))),
    }
}

fn as_bc(sx: &Sx) -> Result<BoundaryCondition, ScenarioError> {
    let items = as_list(sx, "bc")?;
    if items.len() != 6 {
        return Err(err(0, "bc needs region physics kind frame value compat"));
    }
    let value = match &items[4] {
        Sx::Atom(a) if a == "none" => None,
        Sx::List(inner) => {
            let head = as_atom(inner.first().ok_or_else(|| err(0, "empty bc value"))?)?;
            match head {
                "uniform" => Some(BcValue::Uniform(as_qty(&inner[1])?)),
                "signal" => Some(BcValue::Signal(as_signal(&inner[1])?)),
                "profile" => Some(BcValue::Profile(as_profile(&inner[1..])?)),
                other => return Err(err(0, &format!("unknown bc value {other:?}"))),
            }
        }
        _ => return Err(err(0, "bad bc value")),
    };
    let compatibility = match as_atom(&items[5])? {
        "none" => None,
        "incompressible" => Some(Compat::Incompressible),
        other => return Err(err(0, &format!("unknown compat {other:?}"))),
    };
    Ok(BoundaryCondition {
        region: as_str(&items[0])?,
        physics: as_physics(as_atom(&items[1])?)?,
        kind: as_kind(as_atom(&items[2])?)?,
        frame: as_u32(&items[3])?,
        value,
        compatibility,
    })
}

fn as_frame(sx: &Sx) -> Result<Frame, ScenarioError> {
    let items = as_list(sx, "frame")?;
    if items.len() != 4 {
        return Err(err(0, "frame needs id name parent motion"));
    }
    let Sx::List(motion_items) = &items[3] else {
        return Err(err(0, "bad frame motion"));
    };
    let head = as_atom(motion_items.first().ok_or_else(|| err(0, "empty motion"))?)?;
    let rest = &motion_items[1..];
    let motion = match head {
        "fixed" => {
            let q_items = as_list(&rest[0], "quat")?;
            if q_items.len() != 4 {
                return Err(err(0, "quat needs four components"));
            }
            FrameMotion::Fixed {
                orientation: Quat {
                    w: as_f64(&q_items[0])?,
                    x: as_f64(&q_items[1])?,
                    y: as_f64(&q_items[2])?,
                    z: as_f64(&q_items[3])?,
                },
                translation: as_vec3(&rest[1])?,
            }
        }
        "rotating" => {
            let axis = as_vec3(&rest[0])?;
            FrameMotion::Rotating {
                axis: [axis.x, axis.y, axis.z],
                center: as_vec3(&rest[1])?,
                rate: as_qty(&rest[2])?,
            }
        }
        "tilt" => {
            let axis = as_vec3(&rest[0])?;
            FrameMotion::Tilt {
                axis: [axis.x, axis.y, axis.z],
                center: as_vec3(&rest[1])?,
                angle: as_signal(&rest[2])?,
            }
        }
        other => return Err(err(0, &format!("unknown motion {other:?}"))),
    };
    Ok(Frame {
        id: FrameId(as_u32(&items[0])?),
        name: as_str(&items[1])?,
        parent: FrameId(as_u32(&items[2])?),
        motion,
    })
}

fn as_model(sx: &Sx) -> Result<SpectrumModel, ScenarioError> {
    let Sx::List(items) = sx else {
        return Err(err(0, "expected a spectrum model form"));
    };
    let head = as_atom(items.first().ok_or_else(|| err(0, "empty model"))?)?;
    let rest = &items[1..];
    match head {
        "dryden" => Ok(SpectrumModel::Dryden {
            sigma: as_qty(&rest[0])?,
            length_scale: as_qty(&rest[1])?,
            mean_speed: as_qty(&rest[2])?,
        }),
        "kanai-tajimi" => Ok(SpectrumModel::KanaiTajimi {
            s0: as_f64(&rest[0])?,
            omega_g: as_qty(&rest[1])?,
            zeta_g: as_f64(&rest[2])?,
        }),
        "carreau" => {
            if rest.len() != 8 {
                return Err(err(0, "carreau needs six qty bounds + two n bounds"));
            }
            Ok(SpectrumModel::CarreauBand {
                eta_zero: [as_qty(&rest[0])?, as_qty(&rest[1])?],
                eta_inf: [as_qty(&rest[2])?, as_qty(&rest[3])?],
                lambda: [as_qty(&rest[4])?, as_qty(&rest[5])?],
                n: [as_f64(&rest[6])?, as_f64(&rest[7])?],
            })
        }
        other => Err(err(0, &format!("unknown model {other:?}"))),
    }
}

fn as_ensemble(sx: &Sx) -> Result<StochasticEnsemble, ScenarioError> {
    let items = as_list(sx, "ensemble")?;
    if items.len() != 6 {
        return Err(err(0, "ensemble needs name seed members duration dt model"));
    }
    Ok(StochasticEnsemble {
        name: as_str(&items[0])?,
        seed: as_u64(&items[1])?,
        members: as_u32(&items[2])?,
        duration: as_qty(&items[3])?,
        dt: as_qty(&items[4])?,
        model: as_model(&items[5])?,
    })
}

fn as_contact(sx: &Sx) -> Result<ContactLaw, ScenarioError> {
    let contact_items = as_list(sx, "contact")?;
    if contact_items.len() != 3 {
        return Err(err(0, "contact needs two regions + model"));
    }
    let model = match &contact_items[2] {
        Sx::Atom(a) if a == "frictionless" => ContactModel::Frictionless,
        Sx::Atom(a) if a == "tied" => ContactModel::Tied,
        other @ Sx::List(_) => {
            let coulomb = as_list(other, "coulomb")?;
            if coulomb.len() != 2 {
                return Err(err(0, "coulomb needs mu_s + mu_k"));
            }
            ContactModel::Coulomb {
                mu_s: as_f64(&coulomb[0])?,
                mu_k: as_f64(&coulomb[1])?,
            }
        }
        _ => return Err(err(0, "bad contact model")),
    };
    Ok(ContactLaw {
        region_a: as_str(&contact_items[0])?,
        region_b: as_str(&contact_items[1])?,
        model,
    })
}

/// Parse canonical IR text back into a [`Scenario`].
///
/// # Errors
/// [`ScenarioError::Parse`] with a diagnosis for malformed input.
pub fn parse_ir(text: &str) -> Result<Scenario, ScenarioError> {
    let root = parse_sx(text)?;
    let items = as_list(&root, "scenario")?;
    if items.len() != 9 {
        return Err(err(0, "scenario needs name seed + seven sections"));
    }
    let env_items = as_list(&items[2], "environment")?;
    if env_items.len() != 5 {
        return Err(err(
            0,
            "environment needs gravity x3 + temperature + pressure",
        ));
    }
    let environment = Environment {
        gravity: [
            as_qty(&env_items[0])?,
            as_qty(&env_items[1])?,
            as_qty(&env_items[2])?,
        ],
        ambient_temperature: as_qty(&env_items[3])?,
        ambient_pressure: as_qty(&env_items[4])?,
    };
    let mut frames = FrameTree::new();
    for f in as_list(&items[3], "frames")? {
        frames.add(as_frame(f)?);
    }
    let base_bcs = as_list(&items[4], "bcs")?
        .iter()
        .map(as_bc)
        .collect::<Result<Vec<_>, _>>()?;
    let mut cases = Vec::new();
    for c in as_list(&items[5], "cases")? {
        let case_items = as_list(c, "case")?;
        let name = as_str(
            case_items
                .first()
                .ok_or_else(|| err(0, "case needs a name"))?,
        )?;
        let bcs = case_items[1..]
            .iter()
            .map(as_bc)
            .collect::<Result<Vec<_>, _>>()?;
        cases.push(LoadCase { name, bcs });
    }
    let mut combinations = Vec::new();
    for c in as_list(&items[6], "combos")? {
        let combo_items = as_list(c, "combo")?;
        let name = as_str(
            combo_items
                .first()
                .ok_or_else(|| err(0, "combo needs a name"))?,
        )?;
        let mut terms = Vec::new();
        for t in &combo_items[1..] {
            let term_items = as_list(t, "term")?;
            if term_items.len() != 2 {
                return Err(err(0, "term needs case + factor"));
            }
            terms.push((as_str(&term_items[0])?, as_f64(&term_items[1])?));
        }
        combinations.push(Combination { name, terms });
    }
    let ensembles = as_list(&items[7], "ensembles")?
        .iter()
        .map(as_ensemble)
        .collect::<Result<Vec<_>, _>>()?;
    let contacts = as_list(&items[8], "contacts")?
        .iter()
        .map(as_contact)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Scenario {
        name: as_str(&items[0])?,
        seed: as_u64(&items[1])?,
        frames,
        base_bcs,
        cases,
        combinations,
        ensembles,
        contacts,
        environment,
    })
}

/// Round-trip helper for lints/tests: violations if reparse ≠ original.
pub fn check_round_trip(s: &Scenario, out: &mut Vec<Violation>) {
    let text = write_ir(s);
    match parse_ir(&text) {
        Ok(back) if &back == s => {}
        Ok(_) => out.push(Violation {
            code: "ir-round-trip-drift",
            what: format!("scenario {:?} reparses to a different value", s.name),
            fix: "report this as an fs-scenario IR bug (canonical form must be lossless)"
                .to_string(),
        }),
        Err(e) => out.push(Violation {
            code: "ir-round-trip-parse",
            what: format!("scenario {:?} canonical IR failed to reparse: {e}", s.name),
            fix: "report this as an fs-scenario IR bug".to_string(),
        }),
    }
}
