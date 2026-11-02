//! One physical source record shared by independent and coupled performances.
use super::Reader;

pub(super) struct PhysicalTone {
    pub voice: usize,
    pub port: usize,
    pub omega: f64,
    pub cosine: f64,
    pub sine: f64,
}

impl PhysicalTone {
    pub fn read(
        reader: &mut Reader<'_>,
        columns: &[&[Vec<f64>]],
        sample_rate_hz: u32,
        projection_terms: &mut usize,
    ) -> Result<Self, String> {
        let mut row = reader.row("tone").map_err(|e| e.to_string())?;
        let tone = Self {
            voice: row.parse().map_err(|e| e.to_string())?,
            port: row.parse().map_err(|e| e.to_string())?,
            omega: row.scalar().map_err(|e| e.to_string())?,
            cosine: row.scalar().map_err(|e| e.to_string())?,
            sine: row.scalar().map_err(|e| e.to_string())?,
        };
        row.finish().map_err(|e| e.to_string())?;
        let shapes = columns
            .get(tone.voice)
            .and_then(|ports| ports.get(tone.port))
            .ok_or("continuous tone requires an existing voice and physical port")?;
        if tone.omega <= 0.0 || tone.omega >= core::f64::consts::PI * f64::from(sample_rate_hz) {
            return Err("continuous force frequency must be positive and below the original clock's Nyquist limit".into());
        }
        let terms = shapes
            .len()
            .checked_mul(2)
            .ok_or("tone projection size overflow")?;
        *projection_terms = projection_terms
            .checked_sub(terms)
            .ok_or("continuous physical force projection exceeds compile_limits")?;
        Ok(tone)
    }
}
