//! Authored tetrahedral interior geometry, lowered through the original FE and
//! reciprocal material-shell owners. No frequency, overlap or gas coefficient
//! is inferred from an instrument name or accepted as a replacement for geometry.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::ShellCavityOptions;
use fs_couple::render::plate::impact::cavity::MAX_CAVITY_PRESSURE_MODES;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium, MAX_ORDER};
use fs_rep_mesh::TetComplex;
use std::collections::BTreeSet;

#[derive(Clone, Debug)]
pub(super) struct Spec {
    positions: Vec<[f64; 3]>,
    complex: TetComplex,
    order: usize,
    options: ShellCavityOptions,
    pressure_allowance: usize,
}

/// The optional volume follows a complete v3/v4 body. Existing body readers
/// still consume their entire prefix, and this reader consumes the entire tail.
pub(super) fn split(text: &str) -> Result<(&str, Option<Spec>)> {
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        if line.split_ascii_whitespace().next() == Some("interior-fe-v1") {
            return Ok((&text[..offset], Some(Spec::parse(&text[offset..])?)));
        }
        offset += line.len();
    }
    Ok((text, None))
}

fn count(value: &str, minimum: usize, maximum: usize) -> Result<usize> {
    let value: usize = value.parse()?;
    if !(minimum..=maximum).contains(&value) {
        return Err("interior FE count exceeds its explicit source allowance".into());
    }
    Ok(value)
}

impl Spec {
    fn parse(text: &str) -> Result<Self> {
        let mut r = Records { lines: text.lines(), line: 0 };
        let v = r.row("interior-fe-v1", 7)?;
        let order = count(v[0], 1, MAX_ORDER)?;
        let maximum_hz = scalar(v[1])?;
        let maximum_pressure_modes = count(v[2], 1, MAX_CAVITY_PRESSURE_MODES)?;
        let maximum_dofs = count(v[3], 4, 8192)?;
        let maximum_boundary_faces = count(v[4], 4, 8192)?;
        let maximum_projection_terms = count(v[5], 1, 1_000_000_000)?;
        let ritz_tol = scalar(v[6])?;
        if maximum_hz <= 0.0 || ritz_tol <= 0.0 || ritz_tol >= 1.0 {
            return Err("interior FE needs a positive frequency ceiling and relative eigen tolerance in (0,1)".into());
        }
        let options = ShellCavityOptions {
            maximum_hz,
            eigensolver: fs_modal::SliceOptions { ritz_tol, ..Default::default() },
            maximum_pressure_modes, maximum_dofs, maximum_boundary_faces,
            maximum_projection_terms, acoustic_budget: modal_budget(),
        };
        let n = count(r.row("air-nodes", 1)?[0], 4, maximum_dofs)?;
        let mut positions = Vec::with_capacity(n);
        for _ in 0..n {
            let v = r.row("air-node", 3)?;
            let p = [scalar(v[0])?, scalar(v[1])?, scalar(v[2])?];
            if !(0.0..=1.0).contains(&p[0]) || !(0.0..=1.0).contains(&p[1]) {
                return Err("interior x/y fractions must use the original body dimension convention".into());
            }
            positions.push(p);
        }
        let n = count(r.row("air-tets", 1)?[0], 1, 32_768)?;
        let mut tets = Vec::with_capacity(n);
        let mut cells = BTreeSet::new();
        let mut used = vec![false; positions.len()];
        for _ in 0..n {
            let v = r.row("air-tet", 4)?;
            let mut tet = [0u32; 4];
            for i in 0..4 {
                let node = count(v[i], 0, positions.len() - 1)?;
                tet[i] = node as u32;
                used[node] = true;
            }
            let mut key = tet;
            key.sort_unstable();
            if key.windows(2).any(|v| v[0] == v[1]) || !cells.insert(key) {
                return Err("interior tetrahedra need four distinct vertices and unique cells".into());
            }
            tets.push(tet);
        }
        if used.contains(&false) || r.lines.next().is_some() {
            return Err("interior FE must use every declared vertex and consume every record".into());
        }
        let complex = TetComplex::from_tets(positions.len(), tets);
        let dofs = fs_feec::SimplexSpace::new(&complex, order).ndof;
        if dofs > maximum_dofs {
            return Err("complete pressure FE space exceeds the authored degree-of-freedom allowance".into());
        }
        Ok(Self { positions, complex, order, options,
            pressure_allowance: maximum_pressure_modes.min(dofs) })
    }

    pub(super) fn pressure_allowance(&self) -> usize { self.pressure_allowance }

    fn physical_positions(&self, spec: &input::Spec) -> Vec<[f64; 3]> {
        self.positions.iter().map(|p|
            [p[0] * spec.body.length, p[1] * spec.body.width, p[2]]).collect()
    }

    pub(super) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        if spec.cavity.is_some() {
            return Err("interior FE and a separate lumped cavity would duplicate gas storage".into());
        }
        if self.options.maximum_hz >= 0.45 * f64::from(spec.rate) {
            return Err("interior pressure window exceeds the original acoustic Nyquist guard".into());
        }
        let positions = self.physical_positions(spec);
        if positions.iter().flatten().any(|v| !v.is_finite()) {
            return Err("scaled interior coordinates are unrepresentable".into());
        }
        for tet in &self.complex.tets {
            let origin = positions[tet[0] as usize];
            let [a,b,c] = [tet[1],tet[2],tet[3]].map(|i| {
                let p = positions[i as usize];
                [p[0]-origin[0],p[1]-origin[1],p[2]-origin[2]]
            });
            let volume = a[0]*(b[1]*c[2]-b[2]*c[1])
                - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]);
            if !volume.is_finite() || volume == 0.0 {
                return Err("interior FE contains a degenerate physical tetrahedron".into());
            }
        }
        Ok(())
    }

    pub(super) fn bind(&self, bank: ReciprocalStringBank, spec: &input::Spec,
        gate: &CancelGate) -> Result<ReciprocalStringBank>
    {
        poll(gate)?;
        self.preflight(spec)?;
        let positions = self.physical_positions(spec);
        let fe = HelmholtzFe::new(&self.complex, &positions, self.order, Medium {
            density_kg_m3: spec.gas.density,
            sound_speed_m_s: spec.gas.sound_speed,
        })?;
        poll(gate)?;
        // This owner verifies every wall facet and its outward pressure role,
        // retains the complete eigenwindow, and installs the reciprocal load.
        // It refuses missing walls or apertures; no closure surface is invented.
        Ok(bank.with_material_shell_cavity(&fe, self.options, gate)?)
    }
}
