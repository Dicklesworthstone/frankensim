//! Strict, bounded authoring for the existing CutFEM constructive-domain owner.
use super::*;
use fs_cutfem::csg3::{CsgBuilder3, CsgDomain3, CsgError3, CsgNode3, CsgOp3, MAX_CSG3_NODES};

const MAX_DEPTH: usize = 24;

fn construction(error: CsgError3) -> Failure {
    invalid(format!("constructive shape: {error}"))
}

fn visit(
    node: &Node,
    builder: &mut CsgBuilder3,
    depth: usize,
    nodes: &mut usize,
) -> Result<CsgNode3> {
    if depth > MAX_DEPTH || *nodes == MAX_CSG3_NODES {
        return Err(invalid("constructive shape exceeds 24 nesting levels or 128 nodes"));
    }
    *nodes += 1;
    let items = super::super::super::list(node, "constructive shape")?;
    let name = match items.first().map(|n| &n.kind) {
        Some(NodeKind::Symbol(name)) => name.as_str(),
        _ => return Err(invalid("constructive shape requires a named primitive or Boolean")),
    };
    match name {
        "half-space" => {
            let args = fields(node, name, &["normal", "offset-m"])?;
            builder.half_space(vector(args[0])?, scalar(args[1])?).map_err(construction)
        }
        "sphere" => {
            let args = fields(node, name, &["center-m", "radius-m"])?;
            builder.sphere(vector(args[0])?, scalar(args[1])?).map_err(construction)
        }
        "ellipsoid" => {
            let args = fields(node, name, &["center-m", "semi-axes-m"])?;
            builder.ellipsoid(vector(args[0])?, vector(args[1])?).map_err(construction)
        }
        "cylinder" => {
            let args = fields(node, name, &["axis", "center-m", "radius-m"])?;
            let axis = match &args[0].kind {
                NodeKind::Symbol(s) | NodeKind::Str(s) => match s.as_str() {
                    "x" => HeightAxis::X,
                    "y" => HeightAxis::Y,
                    "z" => HeightAxis::Z,
                    _ => return Err(invalid("cylinder axis must be x, y or z")),
                },
                _ => return Err(invalid("cylinder axis must be x, y or z")),
            };
            builder.cylinder(axis, vector(args[1])?, scalar(args[2])?).map_err(construction)
        }
        "box" => {
            let args = fields(node, name, &["center-m", "half-extents-m"])?;
            builder.box_region(vector(args[0])?, vector(args[1])?).map_err(construction)
        }
        "union" | "intersection" | "difference" => {
            let args = fields(node, name, &["blend-m", "left", "right"])?;
            let blend = scalar(args[0])?;
            let a = visit(args[1], builder, depth + 1, nodes)?;
            let b = visit(args[2], builder, depth + 1, nodes)?;
            let op = match name {
                "union" => CsgOp3::Union,
                "intersection" => CsgOp3::Intersection,
                _ => CsgOp3::Difference,
            };
            builder.combine(op, a, b, blend).map_err(construction)
        }
        _ => Err(invalid(format!("unsupported constructive primitive {name:?}"))),
    }
}

pub(super) fn parse(node: &Node) -> Result<CsgDomain3> {
    let mut builder = CsgBuilder3::new();
    let root = visit(node, &mut builder, 1, &mut 0)?;
    builder.finish(root).map_err(construction)
}

#[cfg(test)]
#[path = "csg_tests.rs"]
mod tests;
