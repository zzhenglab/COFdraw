"""Check the displayed pores, precursor edits, project restore and saved drawing in a live app."""
import argparse
import json
import time
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as pw:
        chrome = Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe")
        browser = pw.chromium.launch(executable_path=str(chrome) if chrome.exists() else None)
        page = browser.new_page(viewport={"width": 1450, "height": 970}, accept_downloads=True)
        page.set_default_timeout(90000)
        errors = []
        page.on("pageerror", lambda exc: errors.append(str(exc)))
        def wait(expression, arg=None):
            deadline=time.monotonic()+85
            while not page.evaluate(expression, arg):
                if time.monotonic()>deadline:
                    raise AssertionError((expression, page.locator("#notice").inner_text(), page.locator("#extent-feedback").inner_text(), errors))
                page.wait_for_timeout(100)

        page.goto(args.base_url)
        wait("window.COFdrawStudio?.state.online && !document.querySelector('#assembly-count').disabled")

        def ready(example=None, face=None, seq=None):
            wait("""([example,face,seq])=>{
              const s=COFdrawStudio.state,a=s.recipe?.assembly,ex=document.querySelector('#assembly-example');
              return !s.busy&&!s.inputsDirty&&!s.dirty&&s.result&&a&&(!example||ex.value===example)&&
                (!face||a.pore_view===face)&&(!seq||a.linker_sequence===seq)&&
                document.querySelector('#framework-image').naturalWidth>0&&
                document.querySelector('#extent-options button:not(:disabled)');
            }""", arg=[example, face, seq])

        def extent(value):
            page.locator(f'#extent-options button[value="{value}"]').click()
            wait("""n=>COFdrawExtent.value===n&&!document.querySelector('#extent-bar').hasAttribute('aria-busy')
              && document.querySelector('#framework-image').complete && document.querySelector('#framework-image').naturalWidth>0""", arg=value)

        def three_pores():
            extent("three-pores")
            wait("""() => COFdrawExtent.visible?.summary?.pore_count===3 &&
                document.querySelector('#framework-image').getAttribute('src')==='.'+COFdrawExtent.visible.downloads.svg""")

        def save(ext, path):
            with page.expect_download() as event:
                page.locator("#exports").get_by_role("button", name=ext, exact=True).click()
            event.value.save_as(args.output / path)

        page.select_option("#assembly-count", "1")
        ready("cof-1")
        self_examples = page.locator("#assembly-example option").count() - 1
        page.select_option("#assembly-example", "bfbaeb-cof")
        ready("bfbaeb-cof", face="hexagon")
        page.select_option("#self-face", "triangle")
        ready("bfbaeb-cof", face="triangle")
        assert page.locator("#descriptor-block option").count() == 1
        three_pores()
        extent(2)
        page.screenshot(path=args.output / "self-kagome.png")
        save("SVG", "self-kagome.svg")
        extent(1)
        page.select_option("#assembly-example", "aba-chain")
        ready("aba-chain")
        assert page.locator("#measure-pore").is_disabled()
        assert not page.locator('#extent-options button[value="three-pores"]').count()
        assert "segment" in page.locator("#extent-options").inner_text().lower()
        page.select_option("#assembly-count", "3")
        ready("sioc-cof-1", seq="BCBCBC")
        mixed_examples = page.locator("#assembly-example option").count() - 1
        page.select_option("#assembly-arrangement", "BBB")
        ready("sioc-cof-1", seq="BBB")
        page.screenshot(path=args.output / "sioc-short-triangle.png")
        save("SVG", "sioc-short-triangle.svg")
        page.select_option("#assembly-arrangement", "CCC")
        ready("sioc-cof-1", seq="CCC")
        page.screenshot(path=args.output / "sioc-long-triangle.png")
        save("SVG", "sioc-long-triangle.svg")
        page.select_option("#assembly-arrangement", "BCBCBC")
        ready("sioc-cof-1", seq="BCBCBC")
        save("SVG", "sioc-hexagon.svg")
        extent(2)
        assert page.locator('#pore-type-legend i').count()==3
        page.screenshot(path=args.output / "sioc-neighboring-pores.png")
        save("SVG", "sioc-neighborhood.svg")
        save("MOL", "sioc-neighborhood.mol")
        page.locator("#color-pore-types").uncheck()
        wait("!document.querySelector('#extent-bar').hasAttribute('aria-busy') && COFdrawExtent.value===2")
        assert page.locator('#pore-type-legend i').count()==0
        page.locator("#color-pore-types").check()
        wait("document.querySelectorAll('#pore-type-legend i').length===3")
        three_pores()
        assert page.locator('#pore-type-legend i').count()==3
        with page.expect_download() as event:
            page.locator("#save-project").click()
        project = args.output / "sioc-project.json"
        event.value.save_as(project)
        assert json.loads(project.read_bytes())["drawing_extent"] == "three-pores"
        page.select_option("#assembly-example", "mixed-tapb")
        ready("mixed-tapb")
        page.locator('#project-file').set_input_files(project)
        ready("sioc-cof-1", seq="BCBCBC")
        wait("COFdrawExtent.value==='three-pores' && COFdrawExtent.visible?.summary?.pore_count===3")
        extent(1)
        page.get_by_role("button", name="Edit precursor C").click()
        page.fill("#sketch-smiles", "O=Cc1ccc(-c2ccc(C=O)c(F)c2)cc1")
        page.locator("#sketch-use").click()
        wait("!document.querySelector('#sketch-dialog').open && COFdrawStudio.state.recipe?.assembly?.custom_precursors?.length===1")
        ready()
        assert page.locator("#assembly-example").input_value() == ""
        extent(2)
        save("SVG", "edited-linker-neighborhood.svg")
        page.select_option("#assembly-count","2")
        page.locator("#cof-preset").select_option("tfb-tab-pore-v1")
        wait("""() => document.querySelector('#cof-preset').value==='tfb-tab-pore-v1' && !COFdrawStudio.state.busy &&
            !COFdrawStudio.state.inputsDirty && !COFdrawStudio.state.recipe?.assembly && !!COFdrawExtent.visible""")
        three_pores()
        save("PNG","three-connected-hexagons.png")
        save("SVG","three-connected-hexagons.svg")
        save("MOL","three-connected-hexagons.mol")
        page.select_option("#cof-preset","source:pyurea-source-layout-v1")
        wait("""() => COFdrawStudio.state.sourceId==='pyurea-source-layout-v1' && COFdrawExtent.value===1 &&
            !!COFdrawExtent.visible && !COFdrawStudio.state.busy""")
        assert not page.locator('#extent-options button[value="three-pores"]').count()
        page.select_option("#cof-preset","source:wiley-r198-tppa-so3h-pore-0002")
        wait("""() => COFdrawStudio.state.sourceId==='wiley-r198-tppa-so3h-pore-0002' &&
            !!COFdrawExtent.visible && !COFdrawStudio.state.busy""")
        three_pores()
        assert not errors, errors
        print(json.dumps({"self_examples":self_examples,"mixed_examples":mixed_examples,"script_errors":errors,
                          "passed":["self faces","chain","three mixed faces","neighboring pores","visible exports","restore","custom linker",
                                    "three connected hexagons","source switching"]}))
        browser.close()


if __name__ == "__main__":
    main()
