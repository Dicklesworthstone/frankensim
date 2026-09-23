//! Bounded SI inputs for the existing snare bank, not a second instrument model.
use super::{SnareSet, HelicalWire, WireSpan, LineContact, ModalAcousticState, Error};
use fs_couple::render::plate::impact::string::StringStretching;
use std::io::Read;

const MAX_BYTES: u64 = 64 * 1024;

#[derive(Clone, Copy, Debug)]
pub struct Specification {
    pub engaged: SnareSet,
    pub disengaged_clearance_m: f64,
}
impl Specification {
    /// Every physical row is required, including explicit off/nonlinear choice.
    /// Neither a missing value nor an invalid supplied law gets a stock default.
    pub fn parse(text: &str) -> Result<Self, Error> {
        if text.len() as u64 > MAX_BYTES { return Err("snare specification exceeds 64 KiB".into()); }
        let mut rows = text.lines().enumerate().filter_map(|(i, line)| {
            let row = line.split('#').next().unwrap_or("").trim();
            (!row.is_empty()).then_some((i + 1, row))
        });
        if rows.next().map(|(_, r)| r) != Some("frankensim-snare-spec-v1") {
            return Err("snare specification needs frankensim-snare-spec-v1 header".into());
        }
        let (mut bank, mut coil, mut mechanics, mut contact, mut stretching) = (None, None, None, None, None);
        let mut seen = [false; 5];
        for (line, row) in rows {
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            let bad = || format!("snare specification line {line}: unknown, duplicate, malformed or nonfinite record");
            let kind = match fields[0] {
                "bank" => 0, "coil" => 1, "mechanics" => 2, "contact" => 3, "stretching" => 4,
                _ => return Err(bad().into()),
            };
            if seen[kind] { return Err(bad().into()); }
            seen[kind] = true;
            let number = |index: usize| -> Result<f64, Error> {
                let v = fields[index].parse::<f64>().map_err(|_| bad())?;
                if !v.is_finite() { return Err(bad().into()); }
                Ok(v)
            };
            match (kind, fields.len()) {
                (0, 6) => bank = Some((fields[1].parse::<usize>().map_err(|_| bad())?,
                    fields[2].parse::<usize>().map_err(|_| bad())?, fields[3].parse::<usize>().map_err(|_| bad())?,
                    number(4)?, number(5)?)),
                (1, 5) => coil = Some(HelicalWire { wire_radius_m: number(1)?, coil_radius_m: number(2)?,
                    pitch_m: number(3)?, density_kg_m3: number(4)? }),
                (2, 4) => mechanics = Some((number(1)?, number(2)?, number(3)?)),
                (3, 6) => contact = Some((number(1)?, number(2)?, number(3)?, number(4)?, number(5)?)),
                (4, 2) if fields[1] == "off" => stretching = Some(None),
                (4, 3) => stretching = Some(Some(StringStretching {
                    axial_rigidity_n: number(1)?, maximum_slope: number(2)?,
                })),
                _ => return Err(bad().into()),
            }
        }
        let (strands, modes_per_strand, contact_cells, length_m, width_m) = bank.ok_or("missing snare bank record")?;
        let (tension, bending, damping) = mechanics.ok_or("missing snare mechanics record")?;
        let (clearance, disengaged_clearance_m, stiffness, exponent, loss) = contact.ok_or("missing snare contact record")?;
        let specification = Self { engaged: SnareSet { strands, modes_per_strand, contact_cells, length_m, width_m,
            coil: coil.ok_or("missing snare coil record")?, tension_per_strand_n: tension,
            bending_per_strand_n_m2: bending, damping_per_s: damping, clearance_m: clearance,
            contact_stiffness_per_length: stiffness, contact_exponent: exponent, contact_internal_loss_s_m: loss,
            carrier: None,
            stretching: stretching.ok_or("missing snare stretching record; declare off or axial rigidity and slope limit")?,
        }, disengaged_clearance_m };
        specification.validate()?;
        Ok(specification)
    }
    pub fn load(path: &str) -> Result<Self, Error> {
        let mut text = String::new();
        std::fs::File::open(path)?.take(MAX_BYTES + 1).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    fn validate(self) -> Result<(), Error> {
        let s = self.engaged;
        s.mode_count()?; // Bound counts before any count-driven allocation.
        if !self.disengaged_clearance_m.is_finite() || self.disengaged_clearance_m <= s.clearance_m {
            return Err("disengaged snare clearance must be finite and greater than engaged clearance".into());
        }
        // Delegate material and contact admission to the original physical owners.
        // A single zero-motion strand is enough; no head mesh or eigensolve here.
        let span = WireSpan { endpoints_m: [[-0.5*s.length_m,0.0],[0.5*s.length_m,0.0]],
            linear_density_kg_m: s.coil.linear_density_kg_m()?, tension_n: s.tension_per_strand_n,
            bending_stiffness_n_m2: s.bending_per_strand_n_m2, damping_per_s: vec![s.damping_per_s;s.modes_per_strand] };
        let initial = vec![ModalAcousticState::default();s.modes_per_strand];
        if let Some(law) = s.stretching { span.stretching_body(initial, law)?; }
        else { span.body(initial)?; }
        LineContact::uniform(s.length_m,s.contact_cells,s.clearance_m,s.contact_stiffness_per_length,
            s.contact_exponent,s.contact_internal_loss_s_m,"caller-supplied effective snare law; not a calibration".into())?;
        Ok(())
    }
    pub fn select(self, disengaged: bool) -> SnareSet {
        if disengaged { SnareSet { clearance_m: self.disengaged_clearance_m, ..self.engaged } }
        else { self.engaged }
    }
}

pub fn option(args: &mut Vec<String>) -> Result<Option<String>, Error> {
    let mut path = None;
    let mut i = 0;
    while i < args.len() {
        if args[i] != "--snare-spec" { i += 1; continue; }
        if path.is_some() { return Err("--snare-spec may be supplied only once".into()); }
        let value = args.get(i+1).ok_or("--snare-spec needs an input file")?;
        if value.starts_with("--") { return Err("--snare-spec needs a file, not another option".into()); }
        path = Some(value.clone()); args.drain(i..i+2);
    }
    Ok(path)
}

/// Resolve once before numerical-image admission. None on non-snare commands;
/// never load a supplied snare file for a different instrument or ignore it.
pub fn select(path: Option<&str>, command: &str) -> Result<Option<SnareSet>, Error> {
    let disengaged = match command {
        "snare" | "snare-wav" | "snare-mic" => false,
        "snare-off" | "snare-off-wav" | "snare-off-mic" => true,
        _ => return if path.is_some() { Err("--snare-spec applies only to snare[-off][-wav|-mic]".into()) } else { Ok(None) },
    };
    match path {
        Some(path) => Ok(Some(Specification::load(path)?.select(disengaged))),
        None => Ok(Some(SnareSet::reference(disengaged))),
    }
}

#[cfg(test)]
mod tests;
