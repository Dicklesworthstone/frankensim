//! Authored, source-bound compliant shell attachments; all mechanics stay in
//! the shared material-interface and MITC owners.
use super::*;
use fs_couple::shear_shell::ResolvedShearShell;
use fs_couple::shear_shell::connections::{
    MaterialShellConnection, ShellInterfacePoint, ShellInterfaceSide,
    stiffness_per_area_dims, viscosity_per_area_dims,
};
use fs_matdb::NormalizedInterfacePack;
use fs_material::state_point::resolve_interface_state_point;
use fs_plate::shell::mitc::connections::ShellConnectionBudget;

const MAX_CONNECTIONS: usize = 64;
const MAX_POINTS: usize = 4096;

#[derive(Clone, Debug)]
enum Side {
    Region { name: String, points: Vec<ShellInterfacePoint> },
    Fixed,
}
#[derive(Clone, Debug)]
struct Declaration {
    name: String,
    pack: String,
    temperature: f64,
    omega: f64,
    stiffness: String,
    viscosity: String,
    area: f64,
    direction: [f64; 3],
    sides: [Side; 2],
}
#[derive(Clone, Default)]
pub(super) struct Inputs {
    declarations: Vec<Declaration>,
    budget: Option<ShellConnectionBudget>,
    packs: BTreeMap<String, NormalizedInterfacePack>,
    pack_bytes: usize,
}
impl core::fmt::Debug for Inputs {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        // Do not dump the complete normalized source claims in the render log.
        let hashes: Vec<_> = self.packs.iter()
            .map(|(path, pack)| (path, pack.content_hash())).collect();
        f.debug_struct("ShellInterfaces").field("declarations", &self.declarations)
            .field("budget", &self.budget).field("source_hashes", &hashes).finish()
    }
}
impl Inputs {
    pub(super) fn parse(rows: &mut Records<'_>, nodes: &[Node],
        triangles: &[[usize; 3]], regions: &[RegionInput]) -> Result<Self>
    {
        if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
            != Some("interfaces") {
            return Ok(Self::default());
        }
        let v = rows.row("interfaces", 3)?;
        let number = count(v[0], 1, MAX_CONNECTIONS)?;
        let max_terms = count(v[1], 1, MAX_POINTS)?;
        let max_work = count(v[2], 1, MAX_PROJECTION_WORK)?;
        let budget = ShellConnectionBudget { max_connections: number, max_terms, max_work };
        let mut declarations = Vec::<Declaration>::with_capacity(number);
        let mut points = 0usize;
        for _ in 0..number {
            let v = rows.row("interface", 10)?;
            token(v[0], 128)?;
            if declarations.iter().any(|d| d.name == v[0]) {
                return Err("shell interface names must be unique".into());
            }
            path(v[1])?;
            token(v[4], 128)?;
            token(v[5], 128)?;
            if v[4] == v[5] {
                return Err("shell stiffness and viscosity need distinct source properties".into());
            }
            let temperature = scalar(v[2])?;
            let omega = scalar(v[3])?;
            let area = scalar(v[6])?;
            if temperature <= 0.0 || omega < 0.0 || area <= 0.0 {
                return Err("shell interface needs positive kelvin/area and nonnegative angular frequency".into());
            }
            let direction = unit_direction(&v[7..10])?;
            let left = read_side(rows, nodes, triangles, regions, &mut points, max_terms)?;
            let right = read_side(rows, nodes, triangles, regions, &mut points, max_terms)?;
            if matches!((&left, &right), (Side::Fixed, Side::Fixed)) {
                return Err("shell interface needs at least one physical shell side".into());
            }
            declarations.push(Declaration { name: v[0].into(), pack: v[1].into(),
                temperature, omega, stiffness: v[4].into(), viscosity: v[5].into(),
                area, direction, sides: [left, right] });
        }
        Ok(Self { declarations, budget: Some(budget), packs: BTreeMap::new(), pack_bytes: 0 })
    }

    /// Conservative allowance known before any material artifact is loaded.
    pub(super) fn connection_count(&self) -> usize { self.declarations.len() }

    pub(super) fn preflight(&self, shell: &Spec, spec: &input::Spec) -> Result<()> {
        spec.network_dimensions_with_body_connections(true, self.connection_count())?;
        for declaration in &self.declarations {
            let center = |side: &Side| -> Result<Option<[f64; 3]>> {
                let Side::Region { points, .. } = side else { return Ok(None); };
                let mut result = [0.0; 3];
                for point in points {
                    let position = shell.position(point.node, spec);
                    for i in 0..3 {
                        result[i] += point.weight * (position[i] + point.offset_m[i]);
                    }
                }
                if result.iter().any(|x| !x.is_finite()) {
                    return Err("shell interface application point is unrepresentable".into());
                }
                Ok(Some(result))
            };
            // The shared owner checks the exact coincidence/virtual-work law
            // after binding; preflight checks representability before pack I/O.
            center(&declaration.sides[0])?;
            center(&declaration.sides[1])?;
        }
        Ok(())
    }

    pub(super) fn load(&mut self, gate: &CancelGate,
        loader: &mut impl FnMut(&str) -> Result<NormalizedInterfacePack>) -> Result<()>
    {
        let mut packs = BTreeMap::new();
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        for declaration in &self.declarations {
            poll(gate)?;
            if !packs.contains_key(&declaration.pack) {
                let pack = loader(&declaration.pack)?;
                poll(gate)?;
                let size = pack.to_bytes().len();
                if size > input::MAX_PACK_BYTES || size > remaining {
                    return Err("shell interface packs exceed the original material byte budgets".into());
                }
                remaining -= size;
                packs.insert(declaration.pack.clone(), pack);
            }
        }
        // All-or-nothing cold source loading; a refused replacement keeps the
        // previous complete source set rather than publishing a partial one.
        self.packs = packs;
        self.pack_bytes = input::MAX_TOTAL_PACK_BYTES - remaining;
        Ok(())
    }

    pub(super) fn has_source(&self, name: &str) -> bool {
        self.packs.contains_key(name)
    }

    pub(super) fn loaded_pack_bytes(&self) -> Result<usize> {
        if self.declarations.iter().any(|d| !self.packs.contains_key(&d.pack)) {
            return Err("shell interface packs must be loaded before bulk material construction".into());
        }
        Ok(self.pack_bytes)
    }

    pub(super) fn bind(&self, specimen: ResolvedShearShell, spec: &input::Spec,
        gate: &CancelGate) -> Result<ResolvedShearShell>
    {
        let Some(budget) = self.budget else { return Ok(specimen); };
        let mut inputs = Vec::with_capacity(self.declarations.len());
        for declaration in &self.declarations {
            poll(gate)?;
            let pack = self.packs.get(&declaration.pack)
                .ok_or("shell interface pack was not loaded through the source loader")?;
            let held = spec.material_states.get(&declaration.pack).map_or(&[][..], Vec::as_slice);
            let point = bridge_interfaces::source_query(pack, declaration.temperature,
                declaration.omega, held)?;
            let requirements = [
                Requirement::try_new(&declaration.stiffness, stiffness_per_area_dims(), Admit::Finite)?,
                Requirement::try_new(&declaration.viscosity, viscosity_per_area_dims(), Admit::Finite)?,
            ];
            let interface = resolve_interface_state_point(pack.card(), &point,
                &requirements, Pick::SingleClaimOnly)?;
            let side = |index: usize| match &declaration.sides[index] {
                Side::Region { name, points } => ShellInterfaceSide::Region {
                    name: name.clone(), points: points.clone(),
                },
                // Explicit zero-motion boundary. Its identity is the exact
                // ordered surface role authored in this interface card.
                Side::Fixed => ShellInterfaceSide::Fixed {
                    material: interface.surface_materials()[index].clone(),
                },
            };
            let left = side(0);
            let right = side(1);
            inputs.push(MaterialShellConnection {
                name: declaration.name.clone(), interface,
                stiffness_property: declaration.stiffness.clone(),
                viscosity_property: declaration.viscosity.clone(),
                area_m2: declaration.area, direction: declaration.direction, left, right,
            });
        }
        Ok(specimen.with_interfaces(inputs, budget, gate)?)
    }
}

fn read_side(rows: &mut Records<'_>, nodes: &[Node], triangles: &[[usize; 3]],
    regions: &[RegionInput], total: &mut usize, maximum: usize) -> Result<Side>
{
    match rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().nth(1)) {
        Some("fixed") => { rows.row("side", 1)?; Ok(Side::Fixed) }
        Some("region") => {
            let v = rows.row("side", 3)?;
            let name = v[1];
            let region = regions.iter().find(|r| r.name == name)
                .ok_or("shell interface refers to an unknown region")?;
            let number = count(v[2], 1, nodes.len())?;
            *total = total.checked_add(number).ok_or("shell interface point count overflow")?;
            if *total > maximum { return Err("shell interface exceeds the declared footprint budget".into()); }
            let mut points = Vec::<ShellInterfacePoint>::with_capacity(number);
            let mut sum = 0.0;
            for _ in 0..number {
                let v = rows.row("point", 5)?;
                let node = index(v[0], nodes.len())?;
                if points.iter().any(|p| p.node == node)
                    || !region.triangles.iter().any(|&i| triangles[i].contains(&node)) {
                    return Err("shell interface footprint needs distinct nodes belonging to its source region".into());
                }
                let weight = scalar(v[1])?;
                if weight <= 0.0 { return Err("shell interface footprint weights must be positive".into()); }
                sum += weight;
                points.push(ShellInterfacePoint { node, weight,
                    offset_m: [scalar(v[2])?, scalar(v[3])?, scalar(v[4])?] });
            }
            if !sum.is_finite() || (sum - 1.0).abs() > 1e-12 {
                return Err("shell interface footprint weights must sum to one".into());
            }
            Ok(Side::Region { name: name.into(), points })
        }
        _ => Err("shell interface side must be fixed or region NAME POINTS".into()),
    }
}
fn token(value: &str, maximum: usize) -> Result<()> {
    if value.is_empty() || value.len() > maximum || value.chars().any(char::is_control) {
        return Err("shell interface names must be bounded nonempty tokens".into());
    }
    Ok(())
}
fn path(value: &str) -> Result<()> {
    token(value, 512)?;
    if value.contains('\\') || value.contains(':')
        || value.split('/').any(|part| part.is_empty() || part == "." || part == "..") {
        return Err("shell interface path must stay relative to the performance directory".into());
    }
    Ok(())
}
