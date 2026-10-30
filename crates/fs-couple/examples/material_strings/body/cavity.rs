//! Authored tetrahedral interior geometry, lowered through the original FE and
//! reciprocal material-shell owners. No frequency, overlap or gas coefficient
//! is inferred from an instrument name or accepted as a replacement for geometry.
//! Optional compact openings retain the original FE face IDs. Their cavity
//! coordinates bind before nonlinear attachments; the mandatory later shell
//! BEM stage supplies their joint exterior reaction before any render escapes.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::{
    ShellCavityOpening, ShellCavityOptions, ShellCavityThermalOptions,
    MAX_SHELL_CAVITY_OPENINGS,
};
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
    openings: Vec<ShellCavityOpening>,
    thermal: Option<ShellCavityThermalOptions>,
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
        if used.contains(&false) {
            return Err("interior FE must use every declared vertex".into());
        }
        let mut openings = Vec::new();
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
            == Some("air-openings")
        {
            let n = count(r.row("air-openings", 1)?[0], 1, MAX_SHELL_CAVITY_OPENINGS)?;
            let mut selected = BTreeSet::new();
            for _ in 0..n {
                let v = r.row("air-opening", 3)?;
                let internal_length_m = scalar(v[0])?;
                let internal_resistance_pa_s_m3 = scalar(v[1])?;
                let faces = count(v[2], 1, maximum_boundary_faces)?;
                if internal_length_m <= 0.0 || internal_resistance_pa_s_m3 < 0.0 {
                    return Err("air opening requires positive internal length and nonnegative internal resistance".into());
                }
                if faces > maximum_boundary_faces - selected.len() {
                    return Err("air openings exceed the original boundary-face allowance".into());
                }
                let mut boundary_faces = Vec::with_capacity(faces);
                for _ in 0..faces {
                    let v = r.row("air-opening-face", 3)?;
                    let face = [count(v[0], 0, positions.len()-1)? as u32,
                        count(v[1], 0, positions.len()-1)? as u32,
                        count(v[2], 0, positions.len()-1)? as u32];
                    let mut key = face;
                    key.sort_unstable();
                    if key.windows(2).any(|v| v[0] == v[1]) || !selected.insert(key) {
                        return Err("air openings need distinct original vertices and disjoint boundary faces".into());
                    }
                    // Canonical keys detect repetition only. Preserve authored
                    // FE indices/order for the original geometry/trace owner.
                    boundary_faces.push(face);
                }
                openings.push(ShellCavityOpening { boundary_faces, internal_length_m,
                    internal_resistance_pa_s_m3 });
            }
        }
        let thermal = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("air-thermal")
        {
            let v = r.row("air-thermal", 6)?;
            if v[0] != "isothermal-walls" {
                return Err("interior thermal loss requires an explicit isothermal solid-wall boundary".into());
            }
            let options = ShellCavityThermalOptions {
                maximum_decay_rate_per_s: scalar(v[1])?,
                maximum_thermal_modes: count(v[2], 1, 256)?,
                maximum_thermal_dofs: count(v[3], 1, 8192)?,
                maximum_projection_terms: count(v[4], 1, 1_000_000_000)?,
                maximum_unresolved_fraction: scalar(v[5])?,
            };
            options.validate()?;
            Some(options)
        } else { None };
        if r.lines.next().is_some() {
            return Err("interior FE must consume every record".into());
        }
        let complex = TetComplex::from_tets(positions.len(), tets);
        let dofs = fs_feec::SimplexSpace::new(&complex, order).ndof;
        if dofs > maximum_dofs {
            return Err("complete pressure FE space exceeds the authored degree-of-freedom allowance".into());
        }
        if thermal.is_some_and(|options| dofs > options.maximum_thermal_dofs) {
            return Err("original cavity FE exceeds its thermal degree-of-freedom allowance".into());
        }
        Ok(Self { positions, complex, order, options,
            pressure_allowance: maximum_pressure_modes.min(dofs), openings, thermal })
    }

    pub(super) fn pressure_allowance(&self) -> usize { self.pressure_allowance }
    pub(super) fn opening_count(&self) -> usize { self.openings.len() }
    pub(super) fn has_thermal_loss(&self) -> bool { self.thermal.is_some() }
    pub(super) fn thermal_memory_states(&self) -> usize {
        self.thermal.map_or(0, |options| options.maximum_thermal_modes)
    }

    fn physical_positions(&self, spec: &input::Spec) -> Vec<[f64; 3]> {
        self.positions.iter().map(|p|
            [p[0] * spec.body.length, p[1] * spec.body.width, p[2]]).collect()
    }

    pub(super) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        if spec.cavity.is_some() {
            return Err("interior FE and a separate lumped cavity would duplicate gas storage".into());
        }
        if spec.cavity_molecular.is_some() && self.thermal.is_none() {
            return Err("cavity_molecular requires the original air-thermal boundary and heat eigenwindow".into());
        }
        if let Some(thermal) = self.thermal {
            thermal.validate()?;
            if spec.prepared_linear.is_some() {
                return Err("interior thermal diffusion requires the original exact owner with gas memory".into());
            }
        }
        if self.options.maximum_hz >= 0.45 * f64::from(spec.rate) {
            return Err("interior pressure window exceeds the original acoustic Nyquist guard".into());
        }
        if !self.openings.is_empty() {
            let Some(input::BodyRadiation::ShellBem(radiation)) = spec.body_radiation else {
                return Err("interior air openings require the original shell BEM radiation owner".into());
            };
            if radiation.maximum_hz < self.options.maximum_hz {
                return Err("shell BEM band must cover the complete selected interior pressure window".into());
            }
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
        let gas = *spec.cavity_gas();
        let fe = HelmholtzFe::new(&self.complex, &positions, self.order, Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        })?;
        poll(gate)?;
        // The original owner verifies every wall and opening and constructs
        // the reciprocal ports before subsequent nonlinear attachments. The
        // authored source cannot return until its mandatory BEM stage binds
        // the complete body/opening exterior load.
        if let Some(molecular) = &spec.cavity_molecular {
            let thermal = self.thermal.ok_or("cavity_molecular lost its original air-thermal source")?;
            Ok(bank.with_material_shell_cavity_thermal_molecular(
                &fe, gas, self.options, self.openings.clone(), thermal,
                molecular.relaxation.clone(), molecular.options, gate,
            )?)
        } else {
            Ok(bank.with_material_shell_cavity_in_gas(
                &fe, gas, self.options, self.openings.clone(), self.thermal, gate,
            )?)
        }
    }
}


