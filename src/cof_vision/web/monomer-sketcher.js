/* Small local 2D graph editor. Chemical validation belongs to the RDKit API. */
"use strict";
class MonomerSketcher {
  constructor(svg, onChange = () => {}, onError = () => {}) {
    this.svg = svg; this.onChange = onChange; this.onError = onError; this.graph = {atoms: [], bonds: []};
    this.tool = "bond"; this.element = "C"; this.order = 1; this.selected = null;
    this.history = []; this.future = []; this.drag = null; this.scale = 35; this.center = [0, 0];
    svg.addEventListener("pointerdown", event => { try { this.down(event); } catch (error) { this.onError(error.message); } });
    svg.addEventListener("pointermove", event => this.move(event));
    svg.addEventListener("pointerup", event => { try { this.up(event); } catch (error) { this.drag = null; this.draw(); this.onError(error.message); } });
    svg.addEventListener("pointercancel", () => { this.drag = null; this.draw(); });
    svg.addEventListener("keydown", event => {
      if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === "z") { event.preventDefault(); event.shiftKey ? this.redo() : this.undo(); }
      if (event.key === "Delete" || event.key === "Backspace") { event.preventDefault(); try { this.removeSelected(); } catch (error) { this.onError(error.message); } }
      if (event.key === "Escape") { this.drag = null; this.selected = null; this.draw(); }
    });
    this.draw();
  }
  clone(value) { return JSON.parse(JSON.stringify(value)); }
  remember() { this.history.push(this.clone(this.graph)); if (this.history.length > 100) this.history.shift(); this.future = []; }
  changed() { this.draw(); this.onChange(this.exportGraph()); }
  exportGraph() { return this.clone(this.graph); }
  load(graph) {
    if (!graph || !Array.isArray(graph.atoms) || !Array.isArray(graph.bonds)) throw new Error("The editor needs atoms and bonds.");
    if (graph.atoms.length > 300) throw new Error("The monomer editor supports up to 300 atoms.");
    const atoms = graph.atoms.map(a => ({...a, element: a.element || a.symbol || "C", x: Number(a.x), y: Number(a.y), charge: Number(a.charge || 0)}));
    if (atoms.some(a => !Number.isFinite(a.x) || !Number.isFinite(a.y))) throw new Error("Atom coordinates must be finite.");
    const bonds = graph.bonds.map(b => ({...b, a: Number(b.a), b: Number(b.b), order: Number(b.order)}));
    if (bonds.some(b => !Number.isInteger(b.a) || !Number.isInteger(b.b) || b.a === b.b || !atoms[b.a] || !atoms[b.b] || ![1, 1.5, 2, 3].includes(b.order))) throw new Error("Invalid atom reference or bond order.");
    this.remember(); this.graph = {...this.clone(graph), atoms, bonds}; this.selected = null; this.fit(); this.changed();
  }
  clear() { this.remember(); this.graph = {atoms: [], bonds: []}; this.selected = null; this.fit(); this.changed(); }
  undo() { if (!this.history.length) return; this.future.push(this.clone(this.graph)); this.graph = this.history.pop(); this.selected = null; this.fit(); this.changed(); }
  redo() { if (!this.future.length) return; this.history.push(this.clone(this.graph)); this.graph = this.future.pop(); this.selected = null; this.fit(); this.changed(); }
  setTool(tool, value) { this.tool = tool; if (tool === "atom") this.element = value || "C"; if (tool === "bond") this.order = Number(value || 1); this.drag = null; this.draw(); }
  fit() {
    if (!this.graph.atoms.length) { this.scale = 35; this.center = [0, 0]; return; }
    const xs = this.graph.atoms.map(a => a.x), ys = this.graph.atoms.map(a => a.y);
    const width = Math.max(...xs) - Math.min(...xs), height = Math.max(...ys) - Math.min(...ys);
    this.center = [(Math.min(...xs) + Math.max(...xs)) / 2, (Math.min(...ys) + Math.max(...ys)) / 2];
    this.scale = Math.min(38, 540 / Math.max(3, width), 270 / Math.max(3, height));
  }
  screen(atom) { return [320 + (atom.x - this.center[0]) * this.scale, 180 - (atom.y - this.center[1]) * this.scale]; }
  position(event) {
    const p = new DOMPoint(event.clientX, event.clientY).matrixTransform(this.svg.getScreenCTM().inverse());
    return {x: (p.x - 320) / this.scale + this.center[0], y: (180 - p.y) / this.scale + this.center[1]};
  }
  atomAt(p) { let found = -1, distance = 13 / this.scale; this.graph.atoms.forEach((a, i) => { const d = Math.hypot(a.x - p.x, a.y - p.y); if (d < distance) { found = i; distance = d; } }); return found; }
  bondAt(p) {
    let found = -1, distance = 8 / this.scale;
    this.graph.bonds.forEach((b, i) => { const a = this.graph.atoms[b.a], c = this.graph.atoms[b.b]; const dx = c.x - a.x, dy = c.y - a.y, den = dx * dx + dy * dy; if (!den) return; const t = Math.max(0, Math.min(1, ((p.x - a.x) * dx + (p.y - a.y) * dy) / den)); const d = Math.hypot(p.x - a.x - t * dx, p.y - a.y - t * dy); if (d < distance) { found = i; distance = d; } }); return found;
  }
  hasStereo() { return this.graph.atoms.some(a => a.chiral_tag && !["CHI_UNSPECIFIED", "UNSPECIFIED", "0"].includes(String(a.chiral_tag))) || this.graph.bonds.some(b => (b.stereo && !["STEREONONE", "NONE", "0"].includes(String(b.stereo))) || (b.direction && !["NONE", "0"].includes(String(b.direction)))); }
  connectivityEditable() { if (this.hasStereo()) throw new Error("This monomer has assigned stereochemistry. Move atoms here, or edit its stereospecific SMILES / exported MOL in a full chemical editor to change connectivity."); }
  addAtom(p, element = "C") { if (this.graph.atoms.length >= 300) throw new Error("The editor limit is 300 atoms."); this.graph.atoms.push({element, x: p.x, y: p.y, charge: 0}); return this.graph.atoms.length - 1; }
  addBond(a, b, order = this.order) { if (a === b) return; const old = this.graph.bonds.find(x => (x.a === a && x.b === b) || (x.a === b && x.b === a)); if (old) old.order = order; else this.graph.bonds.push({a, b, order}); }
  deleteAtom(index) { this.graph.atoms.splice(index, 1); this.graph.bonds = this.graph.bonds.filter(b => b.a !== index && b.b !== index).map(b => ({...b, a: b.a > index ? b.a - 1 : b.a, b: b.b > index ? b.b - 1 : b.b})); }
  removeSelected() { if (!this.selected) return; this.connectivityEditable(); this.remember(); if (this.selected.kind === "atom") this.deleteAtom(this.selected.index); else this.graph.bonds.splice(this.selected.index, 1); this.selected = null; this.changed(); }
  charge(value) { if (this.selected?.kind !== "atom") return false; this.remember(); this.graph.atoms[this.selected.index].charge = Math.max(-3, Math.min(3, Number(value))); this.changed(); return true; }
  ring(p, existing, count, aromatic) {
    if (this.graph.atoms.length + count > 300) return;
    this.remember(); const radius = 1.5 / (2 * Math.sin(Math.PI / count));
    const anchor = existing >= 0 ? this.graph.atoms[existing] : p;
    const cx = existing >= 0 ? anchor.x + radius : p.x, cy = anchor.y;
    const ids = [];
    for (let i = 0; i < count; i++) { const angle = Math.PI - i * 2 * Math.PI / count; ids.push(i === 0 && existing >= 0 ? existing : this.addAtom({x: cx + radius * Math.cos(angle), y: cy + radius * Math.sin(angle)})); }
    for (let i = 0; i < count; i++) this.addBond(ids[i], ids[(i + 1) % count], aromatic && i % 2 === 0 ? 2 : 1);
    this.selected = {kind: "atom", index: ids[0]}; this.changed();
  }
  down(event) {
    if (event.button !== 0) return; event.preventDefault(); this.svg.focus();
    const p = this.position(event), atom = this.atomAt(p), bond = atom < 0 ? this.bondAt(p) : -1;
    if (this.tool !== "select") this.connectivityEditable();
    if (this.tool.startsWith("ring")) { this.ring(p, atom, this.tool === "ring5" ? 5 : 6, this.tool === "ring-aromatic"); return; }
    if (this.tool === "erase") { if (atom >= 0 || bond >= 0) { this.remember(); if (atom >= 0) this.deleteAtom(atom); else this.graph.bonds.splice(bond, 1); this.selected = null; this.changed(); } return; }
    if (this.tool === "atom") { this.remember(); if (atom >= 0) this.graph.atoms[atom].element = this.element; else this.addAtom(p, this.element); this.selected = {kind: "atom", index: atom >= 0 ? atom : this.graph.atoms.length - 1}; this.changed(); return; }
    if (this.tool === "select") { this.selected = atom >= 0 ? {kind: "atom", index: atom} : bond >= 0 ? {kind: "bond", index: bond} : null; if (atom >= 0) { this.remember(); this.drag = {kind: "move", atom}; this.svg.setPointerCapture(event.pointerId); } this.draw(); return; }
    if (bond >= 0) { this.remember(); this.graph.bonds[bond].order = this.order; this.selected = {kind: "bond", index: bond}; this.changed(); return; }
    this.drag = {kind: "bond", atom, start: atom >= 0 ? {...this.graph.atoms[atom]} : p, end: p}; this.svg.setPointerCapture(event.pointerId); this.draw();
  }
  move(event) {
    if (!this.drag) return; const p = this.position(event);
    if (this.drag.kind === "move") { const atom = this.graph.atoms[this.drag.atom]; atom.x = p.x; atom.y = p.y; this.draw(); return; }
    const hit = this.atomAt(p);
    if (hit >= 0 && hit !== this.drag.atom) this.drag.end = {...this.graph.atoms[hit]};
    else { const a = Math.round(Math.atan2(p.y - this.drag.start.y, p.x - this.drag.start.x) / (Math.PI / 6)) * Math.PI / 6; const length = event.altKey ? Math.max(.5, Math.hypot(p.x - this.drag.start.x, p.y - this.drag.start.y)) : 1.5; this.drag.end = {x: this.drag.start.x + length * Math.cos(a), y: this.drag.start.y + length * Math.sin(a)}; }
    this.draw();
  }
  up(event) {
    if (!this.drag) return; const drag = this.drag; this.drag = null;
    if (drag.kind === "move") { this.changed(); return; }
    let end = drag.end;
    if (Math.hypot(end.x - drag.start.x, end.y - drag.start.y) < .3) end = {x: drag.start.x + 1.299, y: drag.start.y + .75};
    this.remember(); const a = drag.atom >= 0 ? drag.atom : this.addAtom(drag.start); let b = this.atomAt(end); if (b < 0 || b === a) b = this.addAtom(end); this.addBond(a, b); this.selected = {kind: "atom", index: b}; this.changed();
    if (this.svg.hasPointerCapture(event.pointerId)) this.svg.releasePointerCapture(event.pointerId);
  }
  node(name, attrs = {}, text) { const e = document.createElementNS("http://www.w3.org/2000/svg", name); Object.entries(attrs).forEach(([k, v]) => e.setAttribute(k, v)); if (text !== undefined) e.textContent = text; return e; }
  draw() {
    const s = this.svg; s.replaceChildren();
    if (!this.graph.atoms.length && !this.drag) { s.append(this.node("text", {x:320,y:165,"text-anchor":"middle",class:"sketch-placeholder"}, "Choose a ring or draw a bond")); s.append(this.node("text", {x:320,y:188,"text-anchor":"middle",class:"sketch-placeholder small"}, "Click an atom to change its element. Drag to connect.")); }
    this.graph.bonds.forEach((b, index) => {
      const a = this.screen(this.graph.atoms[b.a]), c = this.screen(this.graph.atoms[b.b]); const d = Math.hypot(c[0]-a[0],c[1]-a[1]) || 1, nx = -(c[1]-a[1])/d, ny = (c[0]-a[0])/d;
      if (this.selected?.kind === "bond" && this.selected.index === index) s.append(this.node("line",{x1:a[0],y1:a[1],x2:c[0],y2:c[1],stroke:"#b8e2d5","stroke-width":12,"stroke-linecap":"round"}));
      const offsets = b.order === 3 ? [-4,0,4] : b.order === 2 || b.order === 1.5 ? [-2.5,2.5] : [0];
      offsets.forEach((o,j) => s.append(this.node("line",{x1:a[0]+nx*o,y1:a[1]+ny*o,x2:c[0]+nx*o,y2:c[1]+ny*o,stroke:"#243c3a","stroke-width":1.8,"stroke-linecap":"round",...(b.order===1.5&&j===1?{"stroke-dasharray":"4 3"}:{})})));
    });
    this.graph.atoms.forEach((atom,index) => {
      const [x,y] = this.screen(atom), selected = this.selected?.kind === "atom" && this.selected.index === index;
      if (selected) s.append(this.node("circle",{cx:x,cy:y,r:12,fill:"#d9f1e8",stroke:"#2b8770"}));
      const degree = this.graph.bonds.filter(b=>b.a===index||b.b===index).length;
      if (atom.element !== "C" || atom.charge || atom.isotope || !degree) {
        const label = `${atom.isotope || ""}${atom.element}${atom.charge ? (Math.abs(atom.charge) > 1 ? Math.abs(atom.charge) : "") + (atom.charge > 0 ? "+" : "−") : ""}`;
        s.append(this.node("rect",{x:x-8*label.length/2,y:y-10,width:8*label.length,height:20,fill:"white"}));
        s.append(this.node("text",{x,y:y+5,"text-anchor":"middle",fill:({N:"#3564a6",O:"#b84444",S:"#9a760d",F:"#39834a",Cl:"#39834a",Br:"#934533"})[atom.element]||"#243c3a",class:"sketch-atom"},label));
      } else if (selected) s.append(this.node("circle",{cx:x,cy:y,r:2.5,fill:"#206b56"}));
    });
    if (this.drag?.kind === "bond") { const a=this.screen(this.drag.start),b=this.screen(this.drag.end); s.append(this.node("line",{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:"#31a483","stroke-width":2,"stroke-dasharray":"5 4"})); }
  }
}
window.MonomerSketcher = MonomerSketcher;
