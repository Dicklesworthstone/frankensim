//! Source-bound Kelvin–Voigt mortars on the original assembled shell geometry.
//! Source resolution stays cold; every elastic and viscous operator remains in
//! the existing material, mortar and shared modal owners.
use super::*;
use fs_couple::plate_viscosity::{
    GLUE_NORMAL_STIFFNESS_PROPERTY, GLUE_NORMAL_VISCOSITY_PROPERTY,
    GLUE_SHEAR_STIFFNESS_1_PROPERTY, GLUE_SHEAR_STIFFNESS_2_PROPERTY,
    GLUE_SHEAR_VISCOSITY_1_PROPERTY, GLUE_SHEAR_VISCOSITY_2_PROPERTY, ResolvedGlueLine,
};
use fs_couple::shear_shell::assembly::MaterialShellComponent;
use fs_couple::shear_shell::connections::{stiffness_per_area_dims, viscosity_per_area_dims};
use fs_matdb::{NormalizedInterfacePack, SelectionPolicy};
use fs_material::state_point::{
    LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA, resolve_interface_state_point,
};

#[derive(Clone, Copy, Debug)]
enum Joint {
    Edge(usize),
    Face(usize),
}
#[derive(Clone, Debug)]
struct Declaration {
    joint: Joint,
    pack: String,
    regions: [String; 2],
    width_m: f64,
}
#[derive(Clone, Debug)]
struct Bound {
    source: ResolvedGlueLine,
    omega_band: (f64, f64),
    query: QueryPoint,
}
#[derive(Clone, Default)]
pub(super) struct Inputs {
    declarations: Vec<Declaration>,
    packs: BTreeMap<String, NormalizedInterfacePack>,
    bindings: Vec<Bound>,
    state: Option<(u64, u64)>,
}
impl core::fmt::Debug for Inputs {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let hashes: Vec<_> = self
            .packs
            .iter()
            .map(|(path, pack)| (path, pack.content_hash()))
            .collect();
        f.debug_struct("AssemblyGlue")
            .field("declarations", &self.declarations)
            .field("source_hashes", &hashes)
            .field("state", &self.state)
            .finish()
    }
}

impl Inputs {
    pub(super) fn read_edge(
        &mut self,
        rows: &mut Records<'_>,
        components: &[Component],
        edges: &[EdgeTie],
        index: usize,
    ) -> Result<()> {
        let Some(d) = read(rows, Joint::Edge(index))? else {
            return Ok(());
        };
        let tie = &edges[index];
        for (side, name) in [&tie.a, &tie.b].into_iter().zip(&d.regions) {
            let shell = &components[side.component].shell;
            let region = region(shell, name)?;
            if side.nodes.windows(2).any(|pair| {
                !region
                    .triangles
                    .iter()
                    .any(|&t| pair.iter().all(|node| shell.triangles[t].contains(node)))
            }) {
                return Err(
                    "every glue edge segment must belong to its named material region".into(),
                );
            }
        }
        self.declarations.push(d);
        Ok(())
    }

    pub(super) fn read_face(
        &mut self,
        rows: &mut Records<'_>,
        components: &[Component],
        faces: &[FaceTie],
        index: usize,
    ) -> Result<()> {
        let Some(d) = read(rows, Joint::Face(index))? else {
            return Ok(());
        };
        let tie = &faces[index];
        for ((component, triangles), name) in [
            (tie.master_component, &tie.master_tris),
            (tie.slave_component, &tie.slave_tris),
        ]
        .into_iter()
        .zip(&d.regions)
        {
            let region = region(&components[component].shell, name)?;
            if triangles.iter().any(|t| !region.triangles.contains(t)) {
                return Err(
                    "every glue face triangle must belong to its named material region".into(),
                );
            }
        }
        self.declarations.push(d);
        Ok(())
    }

    pub(super) fn load(
        &mut self,
        spec: &input::Spec,
        gate: &CancelGate,
        loader: &mut impl FnMut(&str) -> Result<NormalizedInterfacePack>,
    ) -> Result<()> {
        let mut packs = BTreeMap::<String, NormalizedInterfacePack>::new();
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        let mut bindings = Vec::with_capacity(self.declarations.len());
        for d in &self.declarations {
            poll(gate)?;
            if !packs.contains_key(&d.pack) {
                let pack = loader(&d.pack)?;
                poll(gate)?;
                let bytes = pack.to_bytes().len();
                if bytes > input::MAX_PACK_BYTES || bytes > remaining {
                    return Err("glue packs exceed the shared interface source byte budget".into());
                }
                remaining -= bytes;
                packs.insert(d.pack.clone(), pack);
            }
            let held = spec.material_states.get(&d.pack).map_or(&[][..], Vec::as_slice);
            let point = bridge_interfaces::source_query(&packs[&d.pack], spec.temperature,
                spec.omega, held)?;
            bindings.push(resolve(&packs[&d.pack], d, &point)?);
        }
        poll(gate)?;
        self.packs = packs;
        self.bindings = bindings;
        self.state = Some((spec.temperature.to_bits(), spec.omega.to_bits()));
        Ok(())
    }

    pub(super) fn has_source(&self, name: &str) -> bool {
        self.packs.contains_key(name)
    }

    pub(super) fn loaded_pack_bytes(&self) -> Result<usize> {
        if self.bindings.len() != self.declarations.len()
            || self
                .declarations
                .iter()
                .any(|d| !self.packs.contains_key(&d.pack))
        {
            return Err("load every glue interface before constructing material strings".into());
        }
        self.packs.values().try_fold(0usize, |n, p| {
            n.checked_add(p.to_bytes().len())
                .ok_or_else(|| "glue interface source byte count overflow".into())
        })
    }

    pub(super) fn bind(
        &self,
        spec: &input::Spec,
        components: &[MaterialShellComponent],
        edges: &mut [EdgeTie],
        faces: &mut [FaceTie],
        gate: &CancelGate,
    ) -> Result<()> {
        self.loaded_pack_bytes()?;
        if !self.declarations.is_empty()
            && self.state != Some((spec.temperature.to_bits(), spec.omega.to_bits()))
        {
            return Err(
                "glue interfaces must be reloaded at the current declared solid state".into(),
            );
        }
        for (d, bound) in self.declarations.iter().zip(&self.bindings) {
            poll(gate)?;
            let held = spec.material_states.get(&d.pack).map_or(&[][..], Vec::as_slice);
            let point = bridge_interfaces::source_query(&self.packs[&d.pack], spec.temperature,
                spec.omega, held)?;
            if point != bound.query {
                return Err("glue interfaces must be reloaded at the current declared environment".into());
            }
            let (a, b) = match d.joint {
                Joint::Edge(i) => (edges[i].a.component, edges[i].b.component),
                Joint::Face(i) => (faces[i].master_component, faces[i].slave_component),
            };
            let material = |component: usize, name: &str| {
                components[component]
                    .specimen
                    .regions()
                    .iter()
                    .find(|r| r.input().region.region.name == name)
                    .map(|r| r.input().region.material.material())
                    .ok_or("glue source region disappeared during material construction")
            };
            // Resolve again against the ACTUAL two bulk material identities.
            // A renamed pack or material swap cannot retain an unrelated joint.
            let source = ResolvedGlueLine::from_card(
                self.packs[&d.pack].card(),
                (material(a, &d.regions[0])?, material(b, &d.regions[1])?),
                &point,
                SelectionPolicy::SingleClaimOnly,
                d.width_m,
                true,
            )?;
            if source != bound.source {
                return Err(
                    "glue source selection changed between loading and material binding".into(),
                );
            }
            match d.joint {
                Joint::Edge(i) => edges[i].kind = TieKind::Compliant(source.glue_line),
                Joint::Face(i) => faces[i].kind = TieKind::Compliant(source.glue_line),
            }
        }
        Ok(())
    }

    pub(super) fn require_modes(&self, intervals: &[(f64, f64)], gate: &CancelGate) -> Result<()> {
        for bound in &self.bindings {
            for &(lo, hi) in intervals {
                poll(gate)?;
                let lo = fs_math::det::sqrt(lo);
                let hi = fs_math::det::sqrt(hi);
                if !lo.is_finite()
                    || !hi.is_finite()
                    || lo > hi
                    || lo < bound.omega_band.0
                    || hi > bound.omega_band.1
                {
                    return Err(
                        "assembled mode interval lies outside the original glue frequency support"
                            .into(),
                    );
                }
            }
        }
        Ok(())
    }
}

fn read(rows: &mut Records<'_>, joint: Joint) -> Result<Option<Declaration>> {
    if rows
        .lines
        .clone()
        .next()
        .and_then(|s| s.split_ascii_whitespace().next())
        != Some("glue")
    {
        return Ok(None);
    }
    let v = rows.row("glue", 4)?;
    if v[0].len() > 512
        || v[0].chars().any(char::is_control)
        || v[0].contains('\\')
        || v[0].contains(':')
        || v[0]
            .split('/')
            .any(|s| s.is_empty() || s == "." || s == "..")
    {
        return Err("glue source path must stay relative to the performance directory".into());
    }
    let width_m = scalar(v[3])?;
    if width_m <= 0.0 {
        return Err("glue needs an explicit positive physical line width in metres".into());
    }
    Ok(Some(Declaration {
        joint,
        pack: v[0].into(),
        regions: [v[1].into(), v[2].into()],
        width_m,
    }))
}

fn region<'a>(shell: &'a super::super::Spec, name: &str) -> Result<&'a RegionInput> {
    shell
        .regions
        .iter()
        .find(|r| r.name == name)
        .ok_or_else(|| "glue refers to an unknown material region on its ordered side".into())
}

fn resolve(pack: &NormalizedInterfacePack, d: &Declaration, point: &QueryPoint) -> Result<Bound> {
    let card = pack.card();
    if !card
        .claims_for(LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA)
        .is_empty()
    {
        return Err("assembled glue requires Kelvin-Voigt coefficients; a constant-loss anchor needs a different constitutive realization".into());
    }
    let mut requirements = Vec::with_capacity(6);
    for name in [
        GLUE_NORMAL_STIFFNESS_PROPERTY,
        GLUE_SHEAR_STIFFNESS_1_PROPERTY,
        GLUE_SHEAR_STIFFNESS_2_PROPERTY,
    ] {
        requirements.push(Requirement::try_new(
            name,
            stiffness_per_area_dims(),
            Admit::StrictlyPositive,
        )?);
    }
    for name in [
        GLUE_NORMAL_VISCOSITY_PROPERTY,
        GLUE_SHEAR_VISCOSITY_1_PROPERTY,
        GLUE_SHEAR_VISCOSITY_2_PROPERTY,
    ] {
        requirements.push(Requirement::try_new(
            name,
            viscosity_per_area_dims(),
            Admit::NonNegative,
        )?);
    }
    let state = resolve_interface_state_point(card, point, &requirements, Pick::SingleClaimOnly)?;
    let mut band = (0.0_f64, f64::MAX);
    for p in state.properties() {
        if !p.is_constant_at_fixed_environment() {
            return Err("glue coefficients must be constant at the resolved external conditions".into());
        }
        let (lo, hi) = p
            .angular_frequency_band_rad_s()?
            .ok_or("every glue coefficient needs explicit frequency applicability")?;
        band = (band.0.max(lo), band.1.min(hi));
    }
    if band.0 > band.1 {
        return Err("glue stiffness and viscosity frequency support must intersect".into());
    }
    let source = ResolvedGlueLine::from_card(
        card,
        (&state.surface_materials()[0], &state.surface_materials()[1]),
        point,
        SelectionPolicy::SingleClaimOnly,
        d.width_m,
        true,
    )?;
    Ok(Bound {
        source,
        omega_band: band,
        query: point.clone(),
    })
}
