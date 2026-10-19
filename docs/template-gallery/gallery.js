const $ = (selector, root = document) => root.querySelector(selector);
const state = {
  manifest: null,
  selected: null,
  rows: [],
  cache: new Map(),
  request: 0,
  changedOnly: false,
};

function el(tag, text, className) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (className) node.className = className;
  return node;
}

function svg(value, label) {
  const doc = new DOMParser().parseFromString(value, 'image/svg+xml');
  const root = doc.documentElement;
  if (root.localName !== 'svg' || doc.querySelector('parsererror'))
    throw new Error('Invalid drawing');
  root.querySelectorAll('script,foreignObject').forEach((n) => n.remove());
  for (const node of [root, ...root.querySelectorAll('*')]) {
    for (const attr of [...node.attributes]) {
      if (/^on/i.test(attr.name) || (/href$/.test(attr.name) && !attr.value.startsWith('#')))
        node.removeAttribute(attr.name);
    }
  }
  root.setAttribute('role', 'img');
  root.setAttribute('aria-label', label);
  return document.importNode(root, true);
}

async function getTemplate(summary) {
  if (!state.cache.has(summary.id)) {
    const request = (async () => {
      const response = await fetch(summary.data);
      if (!response.ok) throw new Error(`Could not load ${summary.id} (${response.status})`);
      const bytes = new Uint8Array(await response.arrayBuffer());
      // Some hosts decompress .gz responses via Content-Encoding themselves.
      if (bytes[0] !== 0x1f || bytes[1] !== 0x8b)
        return JSON.parse(new TextDecoder().decode(bytes));
      if (!('DecompressionStream' in window))
        throw new Error(
          'This gallery needs a browser with gzip DecompressionStream support (current Chrome, Firefox, Edge or Safari).',
        );
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      return JSON.parse(await new Response(stream).text());
    })();
    state.cache.set(summary.id, request);
    request.catch(() => state.cache.delete(summary.id));
  }
  return state.cache.get(summary.id);
}

function hasFlags(c) {
  return ['enabled', 'disabled'].some(
    (mode) =>
      c[mode] &&
      (c[mode].closeContacts ||
        c[mode].bondCrossings ||
        c[mode].sidechainAtomsInsideRings ||
        c[mode].stereoRoundTrip === false),
  );
}

function filteredRows() {
  const search = $('#search').value.trim().toLowerCase();
  const filter = $('#filter').value;
  const sidechain = $('#substituent').value;
  const rows = state.manifest.templates.filter((t) => {
    const numberMatch = /^\d+$/.test(search) && t.number === Number(search);
    const textMatch =
      !search ||
      t.id.toLowerCase().includes(search) ||
      t.referenceSmiles.toLowerCase().includes(search);
    return (
      (/^\d+$/.test(search) ? numberMatch : textMatch) &&
      (sidechain === 'all' || t.substituents.includes(sidechain)) &&
      (filter !== 'changed' || t.maxRmsd > 0.05) &&
      (filter !== 'flags' || t.flags > 0) &&
      (filter !== 'errors' || t.errors || t.limitation)
    );
  });
  const mode = $('#sort').value;
  rows.sort((a, b) => {
    if (mode === 'rmsd') return b.maxRmsd - a.maxRmsd || a.number - b.number;
    if (mode === 'flags') return b.flags - a.flags || a.number - b.number;
    if (mode === 'size') return b.atomCount - a.atomCount || a.number - b.number;
    return a.number - b.number;
  });
  return rows;
}

function updateList() {
  state.rows = filteredRows();
  const list = $('#template-list');
  list.replaceChildren();
  $('#result-count').textContent =
    `${state.rows.length} of ${state.manifest.templateCount} templates`;
  for (const row of state.rows) {
    const item = el('button', undefined, 'template-item');
    item.dataset.id = row.id;
    item.setAttribute('aria-current', String(row.id === state.selected));
    const left = el('div');
    left.append(
      el('div', row.id, 'item-id'),
      el('div', `${row.atomCount} atoms · ${row.caseCount} molecules`, 'item-size'),
    );
    const right = el('div');
    right.append(el('div', row.maxRmsd.toFixed(2), 'item-rmsd'));
    if (row.errors) right.append(el('div', `${row.errors} failure(s)`, 'item-flags'));
    else if (row.flags) right.append(el('div', `${row.flags} flagged panels`, 'item-flags'));
    item.append(left, right);
    item.addEventListener('click', () => select(row.id));
    list.append(item);
  }
  if (!state.rows.some((t) => t.id === state.selected)) {
    if (state.rows.length) select(state.rows[0].id);
    else {
      state.request++;
      $('#detail').replaceChildren(el('div', 'No templates match these filters.', 'empty'));
      $('#detail').setAttribute('aria-busy', 'false');
    }
  } else {
    select(state.selected, false);
  }
}

function metric(text, flagged) {
  return el('span', text, `metric${flagged ? ' flagged' : ''}`);
}

function panel(caseData, mode) {
  const value = caseData[mode];
  const panel = el('div', undefined, `depiction ${mode}`);
  const label = el('div', undefined, 'mode-label');
  label.append(
    el('span', undefined, 'mode-dot'),
    el('span', `Templates = ${mode === 'enabled' ? 'true' : 'false'}`),
  );
  const image = el('div', undefined, 'molecule-svg');
  image.append(svg(value.svg, `${caseData.label}, ring templates ${mode}`));
  const metrics = el('div', undefined, 'metrics');
  metrics.append(
    metric(`${value.closeContacts} close contacts`, value.closeContacts),
    metric(`${value.bondCrossings} crossings`, value.bondCrossings),
    metric(
      `${value.sidechainAtomsInsideRings} sidechain atoms inside rings`,
      value.sidechainAtomsInsideRings,
    ),
  );
  if (value.stereoRoundTrip !== null)
    metrics.append(
      metric(
        value.stereoRoundTrip ? 'Stereo round-trip OK' : 'Stereo round-trip changed',
        !value.stereoRoundTrip,
      ),
    );
  panel.append(label, image, metrics);
  return panel;
}

function download(text, filename, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = el('a');
  a.href = url;
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function caseCard(c, id, number) {
  const card = el('article', undefined, 'case-card');
  const header = el('div', undefined, 'case-header');
  header.append(el('h3', `${String(number).padStart(2, '0')} / ${c.label}`));
  if (c.error) header.append(el('span', 'Generation failure', 'pill error'));
  else {
    header.append(
      el('span', c.rmsd > 0.05 ? 'Changed' : 'Similar', `pill${c.rmsd > 0.05 ? ' changed' : ''}`),
    );
    if (hasFlags(c)) header.append(el('span', 'Inspect geometry', 'pill warn'));
    if (c.referenceEligible === false)
      header.append(el('span', 'Reference template ineligible', 'pill warn'));
    header.append(el('span', `Aligned RMSD ${c.rmsd.toFixed(3)}`, 'case-rmsd'));
  }
  card.append(header);
  if (c.error) card.append(el('div', c.error, 'error-content'));
  else {
    const pair = el('div', undefined, 'pair');
    pair.append(panel(c, 'enabled'), panel(c, 'disabled'));
    card.append(pair);
  }
  const details = el('details', undefined, 'molecule-details');
  details.append(el('summary', 'Molecule & comparison details'));
  details.append(el('code', c.smiles));
  if (c.sites?.length)
    details.append(
      el('p', `Attachment atom(s): ${c.sites.join(', ')} (1-based reference atom order).`),
    );
  if (c.scaffoldChange) details.append(el('p', c.scaffoldChange));
  if (typeof c.referenceEligible === 'boolean')
    details.append(
      el(
        'p',
        c.referenceEligible
          ? 'This molecule still matches the source template’s SMARTS. The library may choose another matching layout.'
          : 'The sidechain prevents a match to the source template’s SMARTS. Templates=true may use a different template or fall back to the template-free algorithm.',
      ),
    );
  if (c.scaffoldRmsd !== undefined)
    details.append(
      el(
        'p',
        `Scaffold RMSD: ${c.scaffoldRmsd.toFixed(3)}. Coordinates use a nominal bond length of 1.5.`,
      ),
    );
  const actions = el('div', undefined, 'actions');
  const copy = el('button', 'Copy SMILES');
  copy.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(c.smiles);
      copy.textContent = 'Copied';
    } catch {
      copy.textContent = 'Select the SMILES above to copy';
    }
  });
  const save = el('button', 'Download .smi');
  save.addEventListener('click', () =>
    download(`${c.smiles} ${id}-${number}\n`, `${id}-${number}.smi`, 'text/plain'),
  );
  actions.append(copy, save);
  details.append(actions);
  card.append(details);
  return card;
}

async function select(id, changeUrl = true) {
  state.selected = id;
  const token = ++state.request;
  const detail = $('#detail');
  detail.setAttribute('aria-busy', 'true');
  for (const b of document.querySelectorAll('.template-item'))
    b.setAttribute('aria-current', String(b.dataset.id === id));
  if (changeUrl) history.replaceState(null, '', `#${id}`);
  const summary = state.manifest.templates.find((t) => t.id === id);
  if (!summary) return;
  try {
    const row = await getTemplate(summary);
    if (token !== state.request) return;
    const fragment = $('#template-heading').content.cloneNode(true);
    $('.template-id', fragment).textContent = `RING TEMPLATE / ${id}`;
    $('h2', fragment).textContent = `${id} · ${row.atomCount}-atom scaffold`;
    $('.template-description', fragment).textContent =
      `${row.cases.length} molecules · ${row.referenceKind}. The enabled depiction uses the build’s full template library; a particular source layout is not forced.`;
    $('.source-svg', fragment).append(svg(row.sourceSvg, `${id} source template coordinates`));
    const reference = $('.reference-content', fragment);
    reference.append(el('code', row.referenceSmiles));
    const sourceLink = el(
      'a',
      row.referenceLine
        ? 'View public reference molecule ↗'
        : 'View Sketcher’s public template patch ↗',
    );
    sourceLink.href = row.referenceLine
      ? `${state.manifest.referenceRepository}/blob/${state.manifest.referenceCommit}/templates.smi#L${row.referenceLine}`
      : `https://github.com/rachelnwalker/sketcher/blob/${state.manifest.sketcherCommit}/external/rdkit_molecular_templates_pr_29.patch`;
    sourceLink.target = '_blank';
    sourceLink.rel = 'noreferrer';
    reference.append(sourceLink);
    const getData = el('button', 'Download template data');
    getData.addEventListener('click', () =>
      download(JSON.stringify(row, null, 2), `${id}.json`, 'application/json'),
    );
    reference.append(getData);
    const checkbox = $('#changed-only', fragment);
    checkbox.checked = state.changedOnly;
    checkbox.addEventListener('change', () => {
      state.changedOnly = checkbox.checked;
      select(id, false);
    });
    const position = state.rows.findIndex((t) => t.id === id);
    $('.previous', fragment).disabled = position <= 0;
    $('.next', fragment).disabled = position >= state.rows.length - 1;
    $('.previous', fragment).addEventListener('click', () => select(state.rows[position - 1].id));
    $('.next', fragment).addEventListener('click', () => select(state.rows[position + 1].id));
    const sidechain = $('#substituent').value;
    const cases = row.cases.filter(
      (c) =>
        (!state.changedOnly || c.rmsd > 0.05 || c.error) &&
        (sidechain === 'all' || c.substituent === sidechain),
    );
    $('.case-count', fragment).textContent = `${cases.length} of ${row.cases.length} molecules`;
    detail.replaceChildren(fragment);
    if (row.limitation) detail.append(el('div', row.limitation, 'warning-banner'));
    if (!cases.length)
      detail.append(
        el('div', 'No molecules match the current sidechain or difference filter.', 'empty'),
      );
    cases.forEach((c) => detail.append(caseCard(c, id, row.cases.indexOf(c) + 1)));
  } catch (error) {
    if (token === state.request) detail.replaceChildren(el('div', error.message, 'empty'));
  } finally {
    if (token === state.request) detail.setAttribute('aria-busy', 'false');
  }
}

function provenance() {
  const m = state.manifest;
  const box = $('#provenance-content');
  const dl = el('dl');
  const entries = {
    Build: m.runtimeLabel,
    'RDKit runtime': m.rdkitVersion,
    'Sketcher commit': m.sketcherCommit,
    'Template source commit': m.referenceCommit,
    'Runtime templates': String(m.runtimeTemplateCount),
    'Library SHA-256': m.library.sha256,
    Generated: new Date(m.generatedAt).toLocaleString(),
  };
  for (const [key, value] of Object.entries(entries))
    dl.append(el('dt', key), el('dd', value ?? 'Unavailable'));
  box.replaceChildren(dl);
  box.append(
    el(
      'p',
      'Both depictions start from the same sanitized molecule and atom order, with existing conformers removed. RDKit is forced as the coordinate backend. Template use is the only generation setting that changes.',
    ),
  );
  box.append(
    el(
      'p',
      'The disabled layout is rotated and translated onto the enabled scaffold without reflection or scaling. Both drawings share a canvas scale. A pair is marked changed when aligned RMSD exceeds 0.05.',
    ),
  );
  box.append(
    el(
      'p',
      'Flags count nonbonded close contacts below half the median bond length, proper bond intersections, and sidechain atoms inside scaffold ring polygons. They are heuristics for inspection, not chemical validity or drawing-quality verdicts.',
    ),
  );
  box.append(
    el(
      'p',
      'The manifest records the runtime library hash and the Sketcher source recipe separately. The installed native runtime is used directly; the source recipe is not substituted for a measured binary. All molecules come from public template references and generated sidechains.',
    ),
  );
}

async function start() {
  try {
    const response = await fetch('manifest.json');
    if (!response.ok)
      throw new Error(
        'No generated collection found. Run the gallery generator, then serve this directory over HTTP.',
      );
    state.manifest = await response.json();
    const m = state.manifest;
    $('#runtime').textContent = `${m.runtimeLabel} · RDKit ${m.rdkitVersion}`;
    $('#template-count').textContent = m.templateCount.toLocaleString();
    $('#molecule-count').textContent = m.totals.molecules.toLocaleString();
    $('#changed-count').textContent = m.totals.changedPairs.toLocaleString();
    const names = [...new Set(m.templates.flatMap((t) => t.substituents))].sort();
    for (const name of names) {
      const option = el('option', name);
      option.value = name;
      $('#substituent').append(option);
    }
    state.selected = m.templates.some((t) => t.id === location.hash.slice(1))
      ? location.hash.slice(1)
      : m.templates[0]?.id;
    if (!m.complete)
      $('.intro').append(
        el(
          'div',
          'Development sample: this is not the complete template collection.',
          'warning-banner',
        ),
      );
    provenance();
    updateList();
  } catch (error) {
    $('#detail').replaceChildren(el('div', error.message, 'empty'));
    $('#detail').setAttribute('aria-busy', 'false');
  }
}

for (const id of ['search', 'filter', 'sort', 'substituent'])
  $(`#${id}`).addEventListener(
    id === 'search' ? 'input' : 'change',
    () => state.manifest && updateList(),
  );
$('#build-details').addEventListener('click', () => $('#provenance').showModal());
$('#close-details').addEventListener('click', () => $('#provenance').close());
window.addEventListener('hashchange', () => {
  if (state.manifest?.templates.some((t) => t.id === location.hash.slice(1)))
    select(location.hash.slice(1), false);
});
document.addEventListener('keydown', (event) => {
  if (
    !state.manifest ||
    /INPUT|SELECT|TEXTAREA|BUTTON/.test(document.activeElement.tagName) ||
    $('#provenance').open
  )
    return;
  const position = state.rows.findIndex((t) => t.id === state.selected);
  const next =
    event.key === 'ArrowRight' ? position + 1 : event.key === 'ArrowLeft' ? position - 1 : -1;
  if (next >= 0 && next < state.rows.length) {
    event.preventDefault();
    select(state.rows[next].id);
  }
});
start();
