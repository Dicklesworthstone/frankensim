// Edit primitive gas data in the native source. Constitutive evaluation and
// physical domain admission remain in fs_material::gas_species::input.
const standalone = new Set(['gas', 'cavity_gas', 'ambient', 'exterior', 'dg_exterior', 'dg_chimney_gas']);
const segments = new Set(['cylinder', 'cone', 'hole', 'duct_section']);
const numeric = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i;
const finite = value => typeof value === 'string' && numeric.test(value) && Number.isFinite(Number(value));
const sourceRows = source => source.replace(/^\ufeff/, '').split(/\r\n|\r|\n/)
  .map(line => line.trim()).filter(line => line && !line.startsWith('#'));
const resolvedSchema = 'frankensim-resolved-dg-performance-v1';

/** Locate only supported primitive records; aliases, packs and unknown forms stay untouched. */
export function nativeGasRegions(source) {
  const regions = [], rows = [];
  const records = sourceRows(source);
  const resolved = records[0] === resolvedSchema;
  // The finite-atmosphere owners read one original thermochemical record from
  // an attached file, without a performance keyword or a second gas model.
  const gasRecord = records.length === 1 && /^thermochemical(?:\s|$)/.test(records[0]);
  let line = 0;
  for (const match of source.matchAll(/[^\r\n]+|\r\n|\r|\n/g)) {
    if (/^[\r\n]/.test(match[0])) { line++; continue; }
    const bom = match.index === 0 && match[0].startsWith('\ufeff') ? 1 : 0;
    const tokens = [...match[0].slice(bom).matchAll(/[^ \t\v\f]+/g)];
    for (const token of tokens) token.index += bom;
    const words = tokens.map(token => token[0]);
    const key = words[0];
    if (key && !key.startsWith('#')) rows.push({start:match.index, line:line + 1, tokens, words});
    const bare = gasRecord && key === 'thermochemical';
    let first = bare ? 0 : standalone.has(key) ? 1 : -1;
    const inherited = segments.has(key);
    if (inherited) {
      const marker = words.indexOf('gas');
      if (marker > 0) first = marker + 1;
    }
    if (first < 0 || first >= words.length) continue;
    const fields = words.slice(first);
    const resolvedPrimary = resolved && key === 'gas';
    const supportedKinds = bare ? ['thermochemical']
      : resolvedPrimary ? ['dry-air', 'thermochemical'] : ['moist-air', 'mixture', 'thermochemical'];
    let kind = fields[0], count, consumed, values, pressureIndex;
    if (kind === 'dry-air') {
      if (!resolvedPrimary || !finite(fields[1]) || !finite(fields[2])) continue;
      consumed = 3;
      values = {kind, temperature:fields[1], pressure:fields[2], components:[]};
      pressureIndex = first + 2;
    } else if (kind === 'mixture' || kind === 'thermochemical') {
      count = Number(fields[4]); consumed = 5 + 2 * count;
      if (!/^\d+$/.test(fields[4] ?? '') || count < 1 || count > (kind === 'mixture' ? 6 : 128)
          || consumed > fields.length || !['mole', 'mass'].includes(fields[3])
          || !finite(fields[1]) || !finite(fields[2])) continue;
      const components = Array.from({length:count}, (_, i) => ({species:fields[5 + 2 * i], fraction:fields[6 + 2 * i]}));
      if (components.some(component => !finite(component.fraction))) continue;
      values = {kind, temperature:fields[1], pressure:fields[2], basis:fields[3], components};
      pressureIndex = first + 2;
    } else {
      consumed = inherited ? 2 : 3;
      if (fields.length < consumed || !fields.slice(0, consumed).every(finite)) continue;
      kind = 'moist-air';
      values = {kind, temperature:fields[0], pressure:inherited ? '' : fields[1],
        humidity:fields[inherited ? 1 : 2], components:[]};
      pressureIndex = inherited ? null : first + 1;
    }
    if (!supportedKinds.includes(kind)) continue;
    // duct_section may continue with its original loss declaration. Other
    // recognised records consume the rest of the line in their native readers.
    if (key !== 'duct_section' && consumed !== fields.length) continue;
    const start = match.index + tokens[first].index;
    const last = tokens[first + consumed - 1];
    const end = match.index + last.index + last[0].length;
    const pressureToken = pressureIndex === null ? null : tokens[pressureIndex];
    regions.push({source, start, end, line:line + 1, key, inherited, supportedKinds,
      label:`${bare ? 'Attached thermochemical gas' : key}${inherited ? ' gas' : ''} · line ${line + 1}`, values,
      pressureSpan:pressureToken ? [match.index + pressureToken.index, match.index + pressureToken.index + pressureToken[0].length] : null});
  }
  return [...regions, ...fieldStationRegions(source, rows, regions)].sort((a, b) => a.start - b.start);
}

// A station owns only T and fractions. Its column identities, x coordinate and
// common static pressure belong to the complete native bore declaration.
function fieldStationRegions(source, rows, regions) {
  const result = [];
  const count = (value, maximum, minimum = 1) => /^\d+$/.test(value ?? '')
    && Number(value) >= minimum && Number(value) <= maximum;
  const positive = value => finite(value) && Number(value) > 0;
  for (let i = 0; i < rows.length; i++) {
    const header = rows[i], words = header.words;
    if (words[0] !== 'dg_bore' || words[12] !== 'gas_field' || words[14] !== 'mole'
        || !count(words[1], 512) || !count(words[3], 512) || !count(words[4], 8)
        || !count(words[5], 64) || !count(words[10], 2048) || !count(words[11], 64)
        || !count(words[13], 4096, 2) || !count(words[15], 128)
        || !['open', 'closed', 'matched', 'baffled'].includes(words[2])
        || !words.slice(6, 10).every(positive)) continue;
    const segments = Number(words[1]), stations = Number(words[13]), size = Number(words[15]);
    if (words.length !== 16 + size || stations * size > 16384
        || Number(words[3]) < segments || Number(words[11]) < Number(words[4]) + 1
        || Number(words[6]) >= Number(words[7]) || Number(words[8]) >= 1 || Number(words[9]) >= 1) continue;
    const columns = words.slice(16);
    if (new Set(columns).size !== columns.length || columns.some(id => id.includes('#'))) continue;
    const geometry = rows.slice(i + 1, i + 1 + segments);
    if (geometry.length !== segments || geometry.some(({words:fields}) =>
      fields[0] !== 'dg_segment'
      || !((fields[1] === 'cylinder' && fields.length === 4) || (fields[1] === 'cone' && fields.length === 5))
      || !fields.slice(2).every(positive))) continue;
    const length = geometry.reduce((sum, row) => sum + Number(row.words.at(-1)), 0);
    const first = i + 1 + segments, samples = rows.slice(first, first + stations);
    if (!Number.isFinite(length) || samples.length !== stations
        || rows[first + stations]?.words[0] === 'gas_station') continue;
    let previous = -1, valid = true;
    for (const sample of samples) {
      const fields = sample.words, position = Number(fields[1]);
      if (fields[0] !== 'gas_station' || fields.length !== 3 + size
          || !finite(fields[1]) || position <= previous || !positive(fields[2])
          || fields.slice(3).some(value => !finite(value) || Number(value) < 0 || Number(value) > 1)
          || Math.abs(fields.slice(3).reduce((sum, value) => sum + Number(value), 0) - 1) > 1e-12) {
        valid = false; break;
      }
      previous = position;
    }
    if (!valid || Number(samples[0].words[1]) !== 0 || previous !== length) continue;
    const ambient = regions.filter(region => region.start < header.start
      && ['gas', 'ambient'].includes(region.key) && region.pressureSpan).at(-1);
    if (!ambient || !positive(ambient.values.pressure)) continue;
    samples.forEach((sample, index) => {
      const {tokens, words:fields} = sample;
      const last = tokens.at(-1);
      result.push({source, start:sample.start + tokens[2].index,
        end:sample.start + last.index + last[0].length, line:sample.line,
        key:'gas_station', inherited:true, supportedKinds:['thermochemical'], pressureSpan:null,
        label:'DG field line ' + header.line + ' · x = ' + fields[1] + ' m · station ' + (index + 1) + '/' + stations,
        field:{columns:[...columns], headerLine:header.line, position:fields[1], station:index + 1, stations},
        values:{kind:'thermochemical', temperature:fields[2], pressure:ambient.values.pressure, basis:'mole',
          components:columns.map((species, column) => ({species, fraction:fields[3 + column]}))}});
    });
  }
  return result;
}

function gasText(values, inherited) {
  const positive = (value, name) => {
    if (!finite(value) || Number(value) <= 0) throw new Error(`${name} must be a positive finite number.`);
    return value;
  };
  const temperature = positive(values.temperature, 'Temperature in kelvin');
  if (values.kind === 'dry-air') {
    return ['dry-air', temperature, positive(values.pressure, 'Absolute pressure in pascals')].join(' ');
  }
  if (values.kind === 'moist-air') {
    if (!finite(values.humidity) || Number(values.humidity) < 0 || Number(values.humidity) > 1) {
      throw new Error('Relative humidity must be a fraction between 0 and 1.');
    }
    return [temperature, ...(inherited ? [] : [positive(values.pressure, 'Absolute pressure in pascals')]), values.humidity].join(' ');
  }
  if (!['mixture', 'thermochemical'].includes(values.kind) || !['mole', 'mass'].includes(values.basis)) {
    throw new Error('Choose an explicit native gas model and fraction basis.');
  }
  const pressure = positive(values.pressure, 'Absolute pressure in pascals');
  const limit = values.kind === 'mixture' ? 6 : 128;
  if (!values.components?.length || values.components.length > limit) throw new Error(`Choose between 1 and ${limit} gas components.`);
  const names = new Set();
  let sum = 0;
  for (const component of values.components) {
    if (!component.species || /\s|#/.test(component.species) || names.has(component.species)) {
      throw new Error('Each species needs a distinct, nonempty native identity without spaces.');
    }
    names.add(component.species);
    if (!finite(component.fraction) || Number(component.fraction) < 0 || Number(component.fraction) > 1) {
      throw new Error('Each gas fraction must be a finite value between 0 and 1.');
    }
    sum += Number(component.fraction);
  }
  if (Math.abs(sum - 1) > 1e-12) throw new Error(`Gas fractions sum to ${sum}; enter fractions summing to 1. They are not normalized automatically.`);
  return [values.kind, temperature, pressure, values.basis, values.components.length,
    ...values.components.flatMap(component => [component.species, component.fraction])].join(' ');
}

/** Replace the selected source record, optionally changing every editable explicit pressure. */
export function editNativeGasRegion(source, region, values, commonPressure = false) {
  if (source !== region.source) throw new Error('The performance changed. Select the gas region again before applying.');
  if (!region.supportedKinds.includes(values.kind)) {
    throw new Error(`This native gas owner supports ${region.supportedKinds.join(' or ')}.`
      + (region.supportedKinds.includes('dry-air') ? ' A resolved flow cannot evolve a mixture or moist-air snapshot.' : ''));
  }
  if (region.supportedKinds.includes('dry-air')) {
    const hasEnergySolve = sourceRows(source).some(line => /^energy_solve(?:\s|$)/.test(line));
    if (values.kind === 'thermochemical' && !hasEnergySolve) {
      throw new Error('For resolved thermochemical gas, author an energy_solve TOLERANCE_J_PER_KG MAX_ITERATIONS record after fluid_limits in the performance editor. Its numerical limits are not inferred.');
    }
    if (values.kind === 'dry-air' && hasEnergySolve) {
      throw new Error('For resolved dry-air gas, remove the thermochemical energy_solve record in the performance editor before applying. It is not removed automatically.');
    }
  }
  if (region.field) {
    if (commonPressure || values.pressure !== region.values.pressure) {
      throw new Error('A gas-field station inherits static pressure. Select the ambient gas region to change it.');
    }
    if (values.basis !== 'mole' || values.components?.length !== region.field.columns.length
        || values.components.some((component, index) => component.species !== region.field.columns[index])) {
      throw new Error('Keep the gas-field mole basis and declared species columns in their original order. Change the complete header and source registry in the performance editor.');
    }
  }
  const text = gasText(values, region.inherited);
  const next = region.field ? [values.temperature, ...values.components.map(component => component.fraction)].join(' ') : text;
  const edits = [];
  if (JSON.stringify(values) !== JSON.stringify(region.values)) edits.push([region.start, region.end, next]);
  if (commonPressure) {
    if (!finite(values.pressure) || Number(values.pressure) <= 0) throw new Error('Enter the common absolute pressure in pascals.');
    for (const candidate of nativeGasRegions(source)) {
      if (!candidate.pressureSpan || candidate.start === region.start) continue;
      const [start, end] = candidate.pressureSpan;
      if (source.slice(start, end) !== values.pressure) edits.push([start, end, values.pressure]);
    }
  }
  for (const [start, end, text] of edits.sort((a, b) => b[0] - a[0])) source = source.slice(0, start) + text + source.slice(end);
  return source;
}
