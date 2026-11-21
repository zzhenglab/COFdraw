/* Qualified drawing-plane diameter: the same geometry in the viewport and saved figures. */
"use strict";
(() => {
  const S = window.COFdrawStudio, $ = id => document.getElementById(id);
  const NS = 'http://www.w3.org/2000/svg';
  if (!S || !$('canvas-sheet')) return;
  const svg = document.createElementNS(NS, 'svg');
  svg.id = 'pore-measurement-overlay';
  svg.setAttribute('role', 'img');
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  svg.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;z-index:2;pointer-events:none;';
  svg.toggleAttribute('hidden', true);
  $('canvas-sheet').append(svg);
  let active = null, revision = 0;

  function node(tag, attributes = {}, text = null) {
    const value = document.createElementNS(NS, tag);
    for (const [key, attribute] of Object.entries(attributes)) value.setAttribute(key, String(attribute));
    if (text !== null) value.textContent = text;
    return value;
  }
  function clear() {
    revision++; active = null; svg.replaceChildren(); svg.toggleAttribute('hidden', true);
    $('pore-annotation-actions')?.remove();
    $('pore-annotation-note')?.remove();
  }
  function invalidate() {
    clear();
    $('pore-measurement').hidden = true;
  }
  function current(record) {
    return active === record && record.render === S.state.result && record.recipe === S.state.recipe &&
      !S.state.dirty && !S.state.inputsDirty && !S.state.jsonPending;
  }
  function note(message) {
    const p = document.createElement('p'); p.id = 'pore-annotation-note'; p.textContent = message;
    $('pore-measurement').append(p);
  }
  function transform(result, layout) {
    if (!layout || layout.display_projection) throw new Error('Use a flat chemical style to show the measured diameter on the figure.');
    const {annotation} = result, coordinates = layout.atom_coordinates, pixels = layout.atom_pixels;
    const scale = layout.pixels_per_unit;
    if (annotation?.coordinate_system !== 'drawing_xy_cartesian_y_up' || !coordinates?.length ||
        coordinates.length !== pixels?.length || !(scale > 0)) throw new Error('Render this drawing again to place its measurement.');
    for (const anchor of annotation.coordinate_anchors || []) {
      const xy = coordinates[anchor.atom_index];
      if (!xy || Math.hypot(xy[0] - anchor.xy[0], xy[1] - anchor.xy[1]) > 1e-6)
        throw new Error('The measured coordinates differ from this figure. Render and measure again.');
    }
    const offset = [pixels[0][0] - scale * coordinates[0][0], pixels[0][1] + scale * coordinates[0][1]];
    const project = xy => [scale * xy[0] + offset[0], -scale * xy[1] + offset[1]];
    // Reject a nonuniform/unsupported display transform rather than misplace a ruler.
    if (coordinates.some((xy, i) => Math.hypot(...project(xy).map((v, k) => v - pixels[i][k])) > .01))
      throw new Error('This visual transform cannot display a calibrated ruler. Use Standard style.');
    const points = annotation.endpoints?.map(project);
    if (points?.length !== 2 || !points.flat().every(Number.isFinite)) throw new Error('The measurement has no usable diameter endpoints.');
    return points;
  }
  function makeGroup(result, points, canvasSize) {
    const [a, b] = points, center = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
    const width = Math.max(1.5, canvasSize[0] / 500), font = Math.max(20, canvasSize[0] * .034);
    const length = Math.hypot(b[0] - a[0], b[1] - a[1]), head = Math.min(font * .6, length * .15);
    const color = '#205f70', group = node('g', {'data-cofdraw-pore-measurement': '1'});
    group.append(node('title', {}, `Pore opening · ${result.diameter_nm.toFixed(2)} nm`));
    group.append(node('line', {x1: a[0], y1: a[1], x2: b[0], y2: b[1], stroke: color,
      'stroke-width': width, 'data-measured-diameter': result.diameter_nm}));
    // Triangle tips coincide exactly with the reported diameter endpoints.
    group.append(node('path', {d: `M ${a[0]},${a[1]} L ${a[0]+head},${a[1]-head*.43} L ${a[0]+head},${a[1]+head*.43} Z M ${b[0]},${b[1]} L ${b[0]-head},${b[1]-head*.43} L ${b[0]-head},${b[1]+head*.43} Z`, fill: color}));
    const textStyle = {'text-anchor': 'middle', 'font-family': 'Arial, sans-serif', fill: color,
      stroke: '#ffffff', 'stroke-width': font * .23, 'stroke-linejoin': 'round', 'paint-order': 'stroke'};
    group.append(node('text', {...textStyle, x: center[0], y: center[1]-font*.62, 'font-size': font, 'font-weight': '600'}, `${result.diameter_nm.toFixed(2)} nm`));
    return group;
  }
  async function annotatedSVG(record) {
    if (!current(record)) throw new Error('The drawing changed. Measure it again before saving an annotated figure.');
    const response = await fetch(`.${record.render.downloads.svg}`);
    if (!response.ok) throw new Error('The drawing export is unavailable. Render it again.');
    const source = await response.text();
    if (!current(record)) throw new Error('The drawing changed while preparing the annotated figure.');
    const documentSVG = new DOMParser().parseFromString(source, 'image/svg+xml');
    if (documentSVG.querySelector('parsererror')) throw new Error('The drawing SVG could not be read.');
    const root = documentSVG.documentElement;
    if (root.localName !== 'svg') throw new Error('Expected an SVG drawing.');
    const metadata = node('metadata', {'data-cofdraw-measurement': '1'}, JSON.stringify({
      diameter_nm: record.result.diameter_nm, diameter_angstrom: record.result.diameter_angstrom,
      annotation: record.result.annotation, scope: record.result.scope, method: record.result.method,
      sampling_lower_bound: record.result.sampling_lower_bound,
      grid_diagonal_angstrom: record.result.grid_diagonal_angstrom,
      calibration_bond_angstrom: record.result.calibration_bond_angstrom}));
    root.append(documentSVG.importNode(metadata, true), documentSVG.importNode(record.group, true));
    return new XMLSerializer().serializeToString(root);
  }
  async function save(record, format) {
    const controls = $('pore-annotation-actions');
    if (!controls || controls.dataset.busy === 'true') return;
    controls.dataset.busy = 'true';
    controls.querySelectorAll('button').forEach(button => button.disabled = true);
    try {
      const source = await annotatedSVG(record);
      if (format === 'svg') await S.saveText('cofdraw-pore-annotated.svg', source, 'image/svg+xml');
      else {
        const url = URL.createObjectURL(new Blob([source], {type: 'image/svg+xml'}));
        try {
          const image = new Image();
          await new Promise((resolve, reject) => {image.onload = resolve; image.onerror = () => reject(new Error('Could not rasterize the annotated SVG.')); image.src = url;});
          if (!current(record)) throw new Error('The drawing changed. Measure it again before saving.');
          const canvas = document.createElement('canvas');
          [canvas.width, canvas.height] = record.canvasSize;
          canvas.getContext('2d').drawImage(image, 0, 0, canvas.width, canvas.height);
          const data = canvas.toDataURL('image/png');
          if (window.pywebview?.api?.save_image) await window.pywebview.api.save_image('cofdraw-pore-annotated.png', data);
          else {const link = document.createElement('a'); link.href = data; link.download = 'cofdraw-pore-annotated.png'; link.click();}
        } finally {URL.revokeObjectURL(url);}
      }
    } catch (error) {S.notice(error.message, true);}
    finally {controls.dataset.busy = 'false'; controls.querySelectorAll('button').forEach(button => button.disabled = !current(record));}
  }
  function show({result, audit, recipe}) {
    clear();
    if (S.state.dirty || S.state.inputsDirty || S.state.jsonPending || recipe !== S.state.recipe || audit !== S.state.result?.audit) return;
    try {
      const layout = audit.drawing_style?.layout, points = transform(result, layout);
      const canvasSize = layout.canvas_size, group = makeGroup(result, points, canvasSize);
      svg.setAttribute('viewBox', `0 0 ${canvasSize.join(' ')}`);
      svg.setAttribute('aria-label', `Pore opening, ${result.diameter_nm.toFixed(2)} nanometers`);
      svg.append(group); svg.toggleAttribute('hidden', false);
      active = {result, group, canvasSize, recipe, render: S.state.result, revision};
      const record = active, controls = document.createElement('div');
      controls.id = 'pore-annotation-actions'; controls.style.cssText = 'display:flex;gap:8px;flex-wrap:wrap;margin-top:10px';
      for (const format of ['svg', 'png']) {
        const button = document.createElement('button'); button.id = `save-pore-${format}`;
        button.className = 'button outline'; button.textContent = `Save annotated ${format.toUpperCase()}`;
        button.addEventListener('click', () => save(record, format)); controls.append(button);
      }
      const hide = document.createElement('button'); hide.className = 'button quiet'; hide.textContent = 'Remove arrow';
      hide.addEventListener('click', clear); controls.append(hide); $('pore-measurement').append(controls);
    } catch (error) {note(error.message);}
  }
  window.addEventListener('cofdraw:pore-measured', event => show(event.detail));
  window.addEventListener('cofdraw:rendered', invalidate);
  window.addEventListener('cofdraw:extentchanged', invalidate);
  window.addEventListener('cofdraw:layoutchanged', invalidate);
  window.addEventListener('cofdraw:state', () => {if (S.state.dirty || S.state.inputsDirty || S.state.jsonPending) invalidate();});
  window.COFdrawMeasure = {clear, get active() {return active;}};
})();
