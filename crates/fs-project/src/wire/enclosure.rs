//! Canonical enclosure declarations share the project's existing typed AST.
use super::*;
use crate::spec::{ConductionEnclosure, EnclosureRadiatingSurface, EnclosureViewFactorEvidence};

pub(super) fn lower(enclosure: &ConductionEnclosure) -> Result<Node, ProjectError> {
    let mut surfaces = vec![sym("surfaces")];
    for surface in &enclosure.surfaces {
        let mut row = vec![
            sym("surface"),
            kw("name"),
            text(&surface.name),
            kw("target"),
            text(&surface.target),
            kw("card"),
            text(&surface.card),
        ];
        if let Some(claim) = &surface.claim {
            row.extend([kw("claim"), text(claim)]);
        }
        row.extend([kw("query-temperature"), qty(surface.query_temperature)?]);
        surfaces.push(list(row));
    }
    let mut matrix = vec![sym("rows")];
    for factors in &enclosure.view_factors {
        let mut row = vec![sym("row")];
        row.extend(factors.iter().map(|&value| float(value)));
        matrix.push(list(row));
    }
    let evidence = match &enclosure.evidence {
        EnclosureViewFactorEvidence::Analytic { geometry } => {
            list(vec![sym("analytic"), kw("geometry"), text(geometry)])
        }
        EnclosureViewFactorEvidence::ExternalQmc {
            seed: value,
            samples,
            generator,
        } => {
            let samples = i64::try_from(*samples).map_err(|_| ProjectError {
                code: "fsim-enclosure-sample-count",
                detail: "external view-factor sample count exceeds the exact wire integer".into(),
                hint: "retain a positive sample count fitting i64".into(),
            })?;
            list(vec![
                sym("external-qmc"),
                kw("seed"),
                seed(*value),
                kw("samples"),
                int(samples),
                kw("generator"),
                text(generator),
            ])
        }
    };
    let mut declaration = vec![
        sym("enclosure-radiation"),
        kw("surfaces"),
        list(surfaces),
        kw("view-factors"),
        list(matrix),
        kw("evidence"),
        evidence,
        kw("closure-tolerance"),
        float(enclosure.closure_tolerance),
        kw("reciprocity-relative-tolerance"),
        float(enclosure.reciprocity_relative_tolerance),
        kw("max-iterations"),
        int(i64::from(enclosure.max_iterations)),
        kw("temperature-tolerance"),
        qty(enclosure.temperature_tolerance)?,
        kw("heat-tolerance"),
        qty(enclosure.heat_tolerance)?,
        kw("relaxation"),
        float(enclosure.relaxation),
    ];
    if let Some(model) = enclosure.spatial_model {
        declaration.extend([kw("spatial-model"), sym(model.slug())]);
    }
    if let Some(temperature) = enclosure.ambient_temperature {
        declaration.push(kw("ambient-temperature"));
        declaration.push(qty(temperature)?);
    }
    Ok(list(declaration))
}

fn malformed(out: &mut Vec<Violation>, what: impl Into<String>) {
    out.push(Violation {
        code: "project-enclosure-radiation-shape",
        what: what.into(),
        fix: "declare ordered surfaces, explicit matrix rows, source evidence and all convergence controls".into(),
    });
}

pub(super) fn read(node: &Node, out: &mut Vec<Violation>) -> Option<ConductionEnclosure> {
    let Some(("enclosure-radiation", body)) = section_name(node) else {
        malformed(
            out,
            "conduction.enclosure-radiation must be an enclosure-radiation declaration",
        );
        return None;
    };
    let pairs = read_pairs(
        body,
        "enclosure-radiation",
        &[
            "spatial-model",
            "ambient-temperature",
            "surfaces",
            "view-factors",
            "evidence",
            "closure-tolerance",
            "reciprocity-relative-tolerance",
            "max-iterations",
            "temperature-tolerance",
            "heat-tolerance",
            "relaxation",
        ],
        out,
    );
    let Some(("surfaces", nodes)) = field(&pairs, "surfaces").and_then(section_name) else {
        malformed(out, "enclosure-radiation.surfaces must be a surfaces list");
        return None;
    };
    let minimum = if field(&pairs, "ambient-temperature").is_some() { 1 } else { 2 };
    if !(minimum..=64).contains(&nodes.len()) {
        malformed(out, "enclosure-radiation requires at most 64 surfaces, at least two when closed or one with an explicit black environment");
        return None;
    }
    let mut surfaces = Vec::with_capacity(nodes.len());
    for node in nodes {
        let Some(("surface", body)) = section_name(node) else {
            malformed(out, "enclosure surfaces must be surface declarations");
            return None;
        };
        let fields = read_pairs(
            body,
            "enclosure-radiation.surface",
            &["name", "target", "card", "claim", "query-temperature"],
            out,
        );
        surfaces.push(EnclosureRadiatingSurface {
            name: expect_str(field(&fields, "name"), "enclosure.surface.name", out),
            target: expect_str(field(&fields, "target"), "enclosure.surface.target", out),
            card: expect_str(field(&fields, "card"), "enclosure.surface.card", out),
            claim: field(&fields, "claim")
                .map(|value| expect_str(Some(value), "enclosure.surface.claim", out)),
            query_temperature: expect_qty(
                field(&fields, "query-temperature"),
                "enclosure.surface.query-temperature",
                out,
            ),
        });
    }
    let Some(("rows", rows)) = field(&pairs, "view-factors").and_then(section_name) else {
        malformed(out, "enclosure view-factors must be a rows list");
        return None;
    };
    if rows.len() != surfaces.len() {
        malformed(
            out,
            "enclosure view-factor row count differs from its surface count",
        );
        return None;
    }
    let mut view_factors = Vec::with_capacity(rows.len());
    for row in rows {
        let Some(("row", values)) = section_name(row) else {
            malformed(out, "enclosure matrix rows must use the row head");
            return None;
        };
        if values.len() != surfaces.len() {
            malformed(
                out,
                "enclosure view-factor column count differs from its surface count",
            );
            return None;
        }
        view_factors.push(
            values
                .iter()
                .map(|node| expect_float(Some(node), "enclosure.view-factor", out))
                .collect(),
        );
    }
    let evidence = match field(&pairs, "evidence").and_then(section_name) {
        Some(("analytic", body)) => {
            let fields = read_pairs(body, "enclosure.evidence.analytic", &["geometry"], out);
            EnclosureViewFactorEvidence::Analytic {
                geometry: expect_str(
                    field(&fields, "geometry"),
                    "enclosure.evidence.geometry",
                    out,
                ),
            }
        }
        Some(("external-qmc", body)) => {
            let fields = read_pairs(
                body,
                "enclosure.evidence.external-qmc",
                &["seed", "samples", "generator"],
                out,
            );
            let seed = match field(&fields, "seed") {
                Some(Node {
                    kind: NodeKind::Seed(value),
                    ..
                }) => *value,
                _ => {
                    malformed(
                        out,
                        "external enclosure generation requires an explicit hexadecimal seed",
                    );
                    return None;
                }
            };
            let samples = match field(&fields, "samples") {
                Some(Node {
                    kind: NodeKind::Int(value),
                    ..
                }) if *value > 0 => *value as u64,
                _ => {
                    malformed(
                        out,
                        "external enclosure generation requires a positive exact sample count",
                    );
                    return None;
                }
            };
            EnclosureViewFactorEvidence::ExternalQmc {
                seed,
                samples,
                generator: expect_str(
                    field(&fields, "generator"),
                    "enclosure.evidence.generator",
                    out,
                ),
            }
        }
        _ => {
            malformed(
                out,
                "enclosure evidence must be explicit analytic or external-qmc evidence",
            );
            return None;
        }
    };
    let max_iterations = match field(&pairs, "max-iterations") {
        Some(Node {
            kind: NodeKind::Int(value),
            ..
        }) => match u32::try_from(*value) {
            Ok(value) if value > 0 => value,
            _ => {
                malformed(out, "enclosure max-iterations must be positive and fit u32");
                return None;
            }
        },
        _ => {
            malformed(out, "enclosure max-iterations must be an explicit integer");
            return None;
        }
    };
    Some(ConductionEnclosure {
        ambient_temperature: field(&pairs, "ambient-temperature")
            .map(|node| expect_qty(Some(node), "enclosure.ambient-temperature", out)),
        spatial_model: read_radiation_spatial_model(
            field(&pairs, "spatial-model"),
            "enclosure.spatial-model",
            out,
        ),
        surfaces,
        view_factors,
        evidence,
        closure_tolerance: expect_float(
            field(&pairs, "closure-tolerance"),
            "enclosure.closure-tolerance",
            out,
        ),
        reciprocity_relative_tolerance: expect_float(
            field(&pairs, "reciprocity-relative-tolerance"),
            "enclosure.reciprocity-relative-tolerance",
            out,
        ),
        max_iterations,
        temperature_tolerance: expect_qty(
            field(&pairs, "temperature-tolerance"),
            "enclosure.temperature-tolerance",
            out,
        ),
        heat_tolerance: expect_qty(
            field(&pairs, "heat-tolerance"),
            "enclosure.heat-tolerance",
            out,
        ),
        relaxation: expect_float(field(&pairs, "relaxation"), "enclosure.relaxation", out),
    })
}
