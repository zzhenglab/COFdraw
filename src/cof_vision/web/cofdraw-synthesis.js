/* Finite-precursor accessibility scores; never a score of COF synthesis success. */
"use strict";
(() => {
  const S=window.COFdrawStudio, $=id=>document.getElementById(id);
  const label=value=>value?.available&&Number.isFinite(value.score)?`SA Score ${value.score.toFixed(2)} / 10`:'SA Score unavailable';
  const scope='1 easier–10 harder. A molecular heuristic; not a synthesis route or COF reaction/yield prediction.';
  const card=document.createElement('section');card.id='precursor-sa-card';card.className='measurement-card';card.hidden=true;
  card.style.cssText='font-size:11px;padding:9px 12px;margin-bottom:12px;line-height:1.4';
  $('canvas-viewport').after(card);
  function update(){
    for(const block of ['a','b']){
      let el=$(`sa-${block}`);if(!el){el=document.createElement('p');el.id=`sa-${block}`;el.className='small-note';$(`inspector-${block}`).append(el);}
      const m=S.state.monomers[block],valid=m&&m.input===$(`smiles-${block}`).value.trim();
      el.textContent=valid?label(m.synthetic_accessibility):'SA Score updates after the monomer structure check.';el.title=scope;
    }
    const audit=S.state.result?.audit?.precursor_synthetic_accessibility;
    card.replaceChildren();card.hidden=!audit||S.state.dirty||S.state.inputsDirty||S.state.jsonPending;if(card.hidden)return;
    const title=document.createElement('strong');title.textContent='Declared precursor SA scores';card.append(title);
    const list=document.createElement('p');list.style.margin='3px 0';list.textContent=(audit.precursors||[]).map(p=>`${p.label||p.component}: ${label(p.synthetic_accessibility)}`).join(' · ');card.append(list);
    if(audit.complete&&Number.isFinite(audit.mean)&&Number.isFinite(audit.max)){
      const values=document.createElement('p');values.style.margin='3px 0';values.textContent=`Mean precursor score ${audit.mean.toFixed(2)} · Hardest precursor ${audit.max.toFixed(2)}`;card.append(values);
    }
    const note=document.createElement('small');note.textContent=scope;card.append(note);
  }
  window.addEventListener('cofdraw:state',update);window.addEventListener('cofdraw:rendered',update);update();
})();
