//! A `frankensim.shell-assembly.v1` table (e.g. data/instruments/geometry/guitar-esteso-1964-zimmer2024/shell-assembly-h12.tsv)
//! written as a v4 assembled body profile, choice for choice as `guitar_body_modes.rs` assembles it:
//! - every region is `orthotropic12` at grain 0 on its own axis, shear correction 5/6, at its part's thickness;
//! - consistent mass, drilling 1e-3, director-min-dot 0.9;
//! - no support (a free body);
//! - every triangle of a part flagged exterior (cavity) is an exterior (cavity) facet on its positive face;
//! - every tie is a full-rotation rigid edge tie, on the master's negative skin (side -1) or midsurface (0);
//! - the tie budget is 4 (6) per slave node plus 64 segments (slave DOFs), with 0.2 mm gap and pi angle tolerance.
//! What the table does not carry is explicit in `TableOptions`. Positions are translated and scaled into the v4 x/y
//! fractions of the performance body's length and width (the table's own bounding box); z stays in metres. Directors and
//! region axes are normalized to unit length (the table prints them rounded). Nothing else is changed or fitted.
use super::*;

/// The choices a shell-assembly table leaves to its user.
pub(crate) struct TableOptions<'a> {
    /// Squared angular-frequency window [rad^2/s^2] of the assembled body modes.
    pub window: (f64, f64),
    /// The v4 projection-work allowance.
    pub projection_work: usize,
    /// Each table pack name and the profile pack it reads: `@body` (the performance body pack) or a relative path.
    pub packs: &'a [(&'a str, &'a str)],
    /// The performance body record's thickness [m]: each part's thickness is written as a ratio to it.
    pub reference_thickness_m: f64,
    /// The component that carries the string bridges, and each string's bridge point (x, y) [m] in table coordinates,
    /// in string order. Each becomes the nearest node of that component, with a +z force direction.
    pub bridge_component: &'a str,
    pub bridges: &'a [(f64, f64)],
}

/// A converted table: the v4 profile, the body length and width [m] its x/y fractions refer to, how far each bridge point
/// was from the node it became [m], and the table's cavity row (volume [m^3], neck radius [m], neck length [m]) for the
/// performance's cavity record.
pub(crate) struct ConvertedTable {
    pub profile: String,
    pub length_m: f64,
    pub width_m: f64,
    pub bridge_snaps_m: Vec<f64>,
    pub cavity: Option<[f64; 3]>,
}

struct Part<'t> {
    name: &'t str,
    pack: &'t str,
    thickness_m: f64,
    axes: Vec<[f64; 3]>,
    nodes: Vec<[f64; 6]>,
    triangles: Vec<[usize; 4]>,
    exterior: bool,
    cavity: bool,
}

struct Tie<'t> {
    master: &'t str,
    slave: &'t str,
    side: i32,
    master_nodes: Vec<usize>,
    slave_nodes: Vec<usize>,
}

fn unit(v: [f64; 3]) -> Result<[f64; 3]> {
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if !(norm.is_finite() && norm > 0.0) {
        return Err("shell-assembly table direction is not a finite nonzero vector".into());
    }
    Ok([v[0] / norm, v[1] / norm, v[2] / norm])
}

fn floats<const N: usize>(fields: &[&str]) -> Result<[f64; N]> {
    if fields.len() != N {
        return Err(format!("shell-assembly table row needs {N} numbers").into());
    }
    let mut out = [0.0; N];
    for (slot, field) in out.iter_mut().zip(fields) {
        *slot = field.parse()?;
    }
    Ok(out)
}

/// Read the table and write it as a v4 assembled body profile (see the module documentation).
pub(crate) fn convert(table: &str, options: &TableOptions<'_>) -> Result<ConvertedTable> {
    let mut lines = table.lines();
    if lines.next() != Some("frankensim.shell-assembly.v1") {
        return Err("not a frankensim.shell-assembly.v1 table".into());
    }
    let (mut parts, mut ties, mut cavity) = (Vec::<Part<'_>>::new(), Vec::<Tie<'_>>::new(), None);
    let mut next = |what: &str| lines.next().ok_or_else(|| format!("shell-assembly table ends inside {what}"));
    while let Ok(line) = next("a record") {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "meta" | "drive" => {}
            "cavity" if cavity.is_none() => cavity = Some(floats::<3>(&f[1..])?),
            "component" if f.len() == 9 => {
                let (regions, nodes, triangles): (usize, usize, usize) = (f[4].parse()?, f[5].parse()?, f[6].parse()?);
                let mut part = Part {
                    name: f[1], pack: f[2], thickness_m: f[3].parse()?,
                    axes: Vec::with_capacity(regions), nodes: Vec::with_capacity(nodes),
                    triangles: Vec::with_capacity(triangles), exterior: f[7] == "1", cavity: f[8] == "1",
                };
                for _ in 0..regions {
                    let row: Vec<&str> = next("a region")?.split('\t').collect();
                    if row[0] != "region" { return Err("shell-assembly table expected a region row".into()); }
                    part.axes.push(unit(floats::<3>(&row[1..])?)?);
                }
                for _ in 0..nodes {
                    let row: Vec<&str> = next("a node")?.split('\t').collect();
                    if row[0] != "node" { return Err("shell-assembly table expected a node row".into()); }
                    let v = floats::<6>(&row[1..])?;
                    let d = unit([v[3], v[4], v[5]])?;
                    part.nodes.push([v[0], v[1], v[2], d[0], d[1], d[2]]);
                }
                for _ in 0..triangles {
                    let row: Vec<&str> = next("a triangle")?.split('\t').collect();
                    if row[0] != "tri" || row.len() != 5 { return Err("shell-assembly table expected a tri row".into()); }
                    let t = [row[1].parse()?, row[2].parse()?, row[3].parse()?, row[4].parse()?];
                    if t[..3].iter().any(|&n| n >= nodes) || t[3] >= regions {
                        return Err("shell-assembly table triangle names a missing node or region".into());
                    }
                    part.triangles.push(t);
                }
                parts.push(part);
            }
            "tie" if f.len() == 6 => {
                let mut list = |record: &str| -> Result<Vec<usize>> {
                    let row: Vec<&str> = next(record)?.split('\t').collect();
                    if row.len() != 2 || row[0] != record { return Err("shell-assembly table tie lacks its node list".into()); }
                    row[1].split(' ').map(|n| Ok(n.parse()?)).collect()
                };
                let (master_nodes, slave_nodes) = (list("master")?, list("slave")?);
                if master_nodes.len() != f[4].parse::<usize>()? || slave_nodes.len() != f[5].parse::<usize>()? {
                    return Err("shell-assembly table tie node counts disagree".into());
                }
                ties.push(Tie { master: f[1], slave: f[2], side: f[3].parse()?, master_nodes, slave_nodes });
            }
            _ => return Err(format!("unknown or malformed shell-assembly table record {}", f[0]).into()),
        }
    }
    let all = || parts.iter().flat_map(|p| &p.nodes);
    let bounds = |axis: usize| all().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), n| (lo.min(n[axis]), hi.max(n[axis])));
    let ((x0, x1), (y0, y1)) = (bounds(0), bounds(1));
    let (length_m, width_m) = (x1 - x0, y1 - y0);
    if !(length_m > 0.0 && width_m > 0.0) || !(options.reference_thickness_m > 0.0) {
        return Err("shell-assembly table needs a planform extent and a positive reference thickness".into());
    }
    let pack = |name: &str| options.packs.iter().find(|(table, _)| *table == name).map(|(_, profile)| *profile)
        .ok_or_else(|| format!("shell-assembly table pack {name} has no profile pack"));
    let index = |name: &str| parts.iter().position(|p| p.name == name)
        .ok_or_else(|| format!("shell-assembly table names an unknown part {name}"));

    use std::fmt::Write as _;
    let mut out = String::new();
    writeln!(out, "frankensim-string-body-v4\nwindow {} {}\nprojection-work {}\ncomponents {}",
        options.window.0, options.window.1, options.projection_work, parts.len())?;
    for part in &parts {
        writeln!(out, "component {}\ndirector-min-dot 0.9\nmass consistent\ndrilling 0.001\nregions {}",
            part.name, part.axes.len())?;
        let (profile_pack, ratio) = (pack(part.pack)?, part.thickness_m / options.reference_thickness_m);
        for (k, a) in part.axes.iter().enumerate() {
            writeln!(out, "region r{k} {profile_pack} orthotropic12 {ratio} 0 {} {} {} {}", 5.0 / 6.0, a[0], a[1], a[2])?;
        }
        writeln!(out, "nodes {}", part.nodes.len())?;
        for n in &part.nodes {
            writeln!(out, "node {} {} {} {} {} {}", (n[0] - x0) / length_m, (n[1] - y0) / width_m, n[2], n[3], n[4], n[5])?;
        }
        writeln!(out, "triangles {}", part.triangles.len())?;
        for t in &part.triangles {
            writeln!(out, "triangle {} {} {} r{}", t[0], t[1], t[2], t[3])?;
        }
        let facets = |on: bool| if on { part.triangles.len() } else { 0 };
        writeln!(out, "supports 0\npressure-surfaces {} {}", facets(part.exterior), facets(part.cavity))?;
        for (role, on) in [("exterior-facet", part.exterior), ("cavity-facet", part.cavity)] {
            for t in 0..facets(on) {
                writeln!(out, "{role} {t} positive")?;
            }
        }
        writeln!(out, "end-component")?;
    }
    let slaves: usize = ties.iter().map(|t| t.slave_nodes.len()).sum();
    writeln!(out, "tie-budget {} {} 0.0002 {}\nedge-ties {}", 4 * slaves + 64, 6 * slaves + 64, core::f64::consts::PI,
        ties.len())?;
    for tie in &ties {
        let offset = match tie.side {
            0 => "midsurface",
            -1 => "negative-skin",
            1 => "positive-skin",
            _ => return Err("shell-assembly table tie side must be -1, 0 or 1".into()),
        };
        index(tie.master)?;
        index(tie.slave)?;
        writeln!(out, "edge-tie {} {} {offset} full {} {}", tie.master, tie.slave, tie.master_nodes.len(),
            tie.slave_nodes.len())?;
        for n in &tie.master_nodes { writeln!(out, "master-node {n}")?; }
        for n in &tie.slave_nodes { writeln!(out, "slave-node {n}")?; }
    }
    let carrier = &parts[index(options.bridge_component)?];
    writeln!(out, "face-ties 0\nbridges {}", options.bridges.len())?;
    let mut bridge_snaps_m = Vec::with_capacity(options.bridges.len());
    for &(x, y) in options.bridges {
        let distance = |n: &[f64; 6]| ((n[0] - x).powi(2) + (n[1] - y).powi(2)).sqrt();
        let (node, snap) = carrier.nodes.iter().enumerate().map(|(i, n)| (i, distance(n)))
            .min_by(|a, b| a.1.total_cmp(&b.1)).ok_or("bridge component has no node")?;
        writeln!(out, "bridge {} {node} 0 0 1", options.bridge_component)?;
        bridge_snaps_m.push(snap);
    }
    Ok(ConvertedTable { profile: out, length_m, width_m, bridge_snaps_m, cavity })
}
